#include "InfernalBridge.h"

#include <cerrno>
#include <condition_variable>
#include <cstdint>
#include <cstdlib>
#include <csetjmp>
#include <cstring>
#include <deque>
#include <fstream>
#include <mutex>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/mman.h>
#include <sys/wait.h>
#include <fcntl.h>
#include <unistd.h>
#include <vector>

extern "C" {
int infernal_cmbuild_main(int argc, char **argv);
int infernal_cmcalibrate_main(int argc, char **argv);
void infernal_cmcalibrate_set_capture_to_memory(int enabled);
const char *infernal_cmcalibrate_get_captured_cm_text(size_t *opt_len);
}

namespace {

static bool fileExists(const std::string &p) {
    struct stat st;
    return stat(p.c_str(), &st) == 0;
}

static std::string makeShmName() {
    std::ostringstream os;
    os << "/mmseqs-infernal-" << static_cast<long long>(getpid()) << "-" << static_cast<long long>(std::rand());
    return os.str();
}

static bool createInMemoryFd(const std::string *content, int &fdOut, std::string &pathOut, std::string &error) {
    fdOut = -1;
    pathOut.clear();

    int fd = -1;
    std::string shmName;
    for (int attempt = 0; attempt < 16; ++attempt) {
        shmName = makeShmName();
        fd = shm_open(shmName.c_str(), O_CREAT | O_EXCL | O_RDWR, 0600);
        if (fd >= 0) {
            break;
        }
    }
    if (fd >= 0) {
        // Unlink immediately; object lives as long as fd is open.
        shm_unlink(shmName.c_str());
    } else {
        // Fallback for environments where shm_open is unavailable/restricted.
        char tmpTemplate[] = "/tmp/mmseqs-infernal-fd-XXXXXX";
        fd = mkstemp(tmpTemplate);
        if (fd < 0) {
            error = "failed to create anonymous in-memory fd";
            return false;
        }
        unlink(tmpTemplate);
    }

    const size_t n = (content != NULL) ? content->size() : 0;
    if (n > 0) {
        if (ftruncate(fd, static_cast<off_t>(n)) != 0) {
            error = "ftruncate failed for in-memory file";
            close(fd);
            return false;
        }
        size_t written = 0;
        while (written < n) {
            const ssize_t w = write(fd, content->data() + written, n - written);
            if (w <= 0) {
                error = "write failed for in-memory file";
                close(fd);
                return false;
            }
            written += static_cast<size_t>(w);
        }
        if (lseek(fd, 0, SEEK_SET) < 0) {
            error = "lseek failed for in-memory file";
            close(fd);
            return false;
        }
    }

    std::ostringstream fdPath;
    fdPath << "/dev/fd/" << fd;
    fdOut = fd;
    pathOut = fdPath.str();
    return true;
}

static thread_local bool gBridgeExitActive = false;
static thread_local int gBridgeExitCode = 1;
static thread_local std::jmp_buf gBridgeExitJmp;

extern "C" void infernal_bridge_exit(int code) {
    if (gBridgeExitActive) {
        gBridgeExitCode = code;
        std::longjmp(gBridgeExitJmp, 1);
    }
    std::exit(code);
}

static int runInfernalMainInProcess(int (*fn)(int, char **), const std::vector<std::string> &args) {
    std::vector<char *> argv;
    argv.reserve(args.size() + 1);
    for (size_t i = 0; i < args.size(); ++i) {
        argv.push_back(const_cast<char *>(args[i].c_str()));
    }
    argv.push_back(NULL);
    optind = 1;
    gBridgeExitCode = 1;
    gBridgeExitActive = true;
    if (setjmp(gBridgeExitJmp) == 0) {
        const int rc = fn(static_cast<int>(args.size()), argv.data());
        gBridgeExitActive = false;
        return rc;
    }
    gBridgeExitActive = false;
    return gBridgeExitCode;
}

} // namespace

namespace {

static bool readWholeFile(const std::string &path, std::string &out, std::string &error) {
    int fd = open(path.c_str(), O_RDONLY);
    if (fd < 0) {
        error = "failed to open CM file for readback";
        return false;
    }
    out.clear();
    char buf[65536];
    for (;;) {
        ssize_t r = read(fd, buf, sizeof(buf));
        if (r < 0) {
            if (errno == EINTR) continue;
            error = "read failed on CM file";
            close(fd);
            return false;
        }
        if (r == 0) break;
        out.append(buf, static_cast<size_t>(r));
    }
    close(fd);
    return true;
}

static bool runCmBuildInProcess(const std::string &stockholmText, std::string &cmText, std::string &error, bool calibrate) {
    char tmpTemplate[] = "/tmp/mmseqs-infernal-XXXXXX";
    char *tmpDir = mkdtemp(tmpTemplate);
    if (tmpDir == NULL) {
        error = "mkdtemp failed";
        return false;
    }
    const std::string dir(tmpDir);
    int stoFd = -1;
    std::string stoPath;
    const std::string cmPath = dir + "/model.cm";

    if (!createInMemoryFd(&stockholmText, stoFd, stoPath, error)) {
        rmdir(dir.c_str());
        return false;
    }

    std::vector<std::string> cmbArgs;
    cmbArgs.push_back("infernal-cmbuild");
    cmbArgs.push_back("-F");
    cmbArgs.push_back(cmPath);
    cmbArgs.push_back(stoPath);
    bool ok = (runInfernalMainInProcess(&infernal_cmbuild_main, cmbArgs) == 0);
    if (!ok) {
        error = "embedded infernal cmbuild failed";
        close(stoFd);
        unlink(cmPath.c_str());
        rmdir(dir.c_str());
        return false;
    }

    if (calibrate) {
        infernal_cmcalibrate_set_capture_to_memory(1);
        std::vector<std::string> calArgs;
        calArgs.push_back("infernal-cmcalibrate");
        calArgs.push_back("-L");
        calArgs.push_back("0.01");
        calArgs.push_back("--cpu");
        calArgs.push_back("1");
        calArgs.push_back(cmPath);
        if (runInfernalMainInProcess(&infernal_cmcalibrate_main, calArgs) != 0) {
            infernal_cmcalibrate_set_capture_to_memory(0);
            error = "embedded infernal cmcalibrate failed";
            close(stoFd);
            unlink(cmPath.c_str());
            rmdir(dir.c_str());
            return false;
        }
        size_t cmLen = 0;
        const char *calibratedCm = infernal_cmcalibrate_get_captured_cm_text(&cmLen);
        if (calibratedCm == NULL || cmLen == 0) {
            infernal_cmcalibrate_set_capture_to_memory(0);
            error = "embedded infernal did not return calibrated CM text";
            close(stoFd);
            unlink(cmPath.c_str());
            rmdir(dir.c_str());
            return false;
        }
        cmText.assign(calibratedCm, cmLen);
        infernal_cmcalibrate_set_capture_to_memory(0);
    } else {
        ok = readWholeFile(cmPath, cmText, error);
    }

    close(stoFd);
    unlink(cmPath.c_str());
    rmdir(dir.c_str());
    return ok;
}

static bool writeAllFd(int fd, const void *buf, size_t n) {
    const char *p = static_cast<const char *>(buf);
    while (n > 0) {
        ssize_t w = write(fd, p, n);
        if (w < 0) {
            if (errno == EINTR) continue;
            return false;
        }
        if (w == 0) return false;
        p += w;
        n -= static_cast<size_t>(w);
    }
    return true;
}

static bool readAllFd(int fd, void *buf, size_t n) {
    char *p = static_cast<char *>(buf);
    while (n > 0) {
        ssize_t r = read(fd, p, n);
        if (r < 0) {
            if (errno == EINTR) continue;
            return false;
        }
        if (r == 0) return false;
        p += r;
        n -= static_cast<size_t>(r);
    }
    return true;
}

} // namespace

namespace {

struct Worker {
    pid_t pid;
    int toWorker;   // parent writes, worker reads
    int fromWorker; // worker writes, parent reads
    std::mutex io;  // serialize send/recv on this worker
};

static void workerMainLoop(int fdIn, int fdOut) {
    for (;;) {
        uint8_t flag = 0;
        uint64_t len = 0;
        if (!readAllFd(fdIn, &flag, 1)) break;
        if (!readAllFd(fdIn, &len, sizeof(len))) break;
        std::string sto;
        if (len > 0) {
            if (len > (1ULL << 30)) break;
            sto.resize(static_cast<size_t>(len));
            if (!readAllFd(fdIn, &sto[0], len)) break;
        }
        std::string cmText;
        std::string err;
        bool ok = runCmBuildInProcess(sto, cmText, err, flag != 0);
        const std::string &payload = ok ? cmText : err;
        uint8_t outFlag = ok ? 1 : 0;
        uint64_t outLen = static_cast<uint64_t>(payload.size());
        if (!writeAllFd(fdOut, &outFlag, 1)) break;
        if (!writeAllFd(fdOut, &outLen, sizeof(outLen))) break;
        if (outLen > 0 && !writeAllFd(fdOut, payload.data(), payload.size())) break;
    }
    close(fdIn);
    close(fdOut);
}

} // namespace

namespace InfernalBridge {

struct WorkerPool {
    std::vector<Worker *> workers;
    std::mutex mu;
    std::condition_variable cv;
    std::deque<int> free;
};

bool isConfigured() {
    return true;
}

WorkerPool *startWorkerPool(int nWorkers) {
    if (nWorkers <= 0) return NULL;
    WorkerPool *pool = new WorkerPool();
    pool->workers.reserve(nWorkers);
    for (int i = 0; i < nWorkers; ++i) {
        int p2w[2];
        int w2p[2];
        if (pipe(p2w) != 0) {
            stopWorkerPool(pool);
            return NULL;
        }
        if (pipe(w2p) != 0) {
            close(p2w[0]); close(p2w[1]);
            stopWorkerPool(pool);
            return NULL;
        }
        pid_t pid = fork();
        if (pid < 0) {
            close(p2w[0]); close(p2w[1]);
            close(w2p[0]); close(w2p[1]);
            stopWorkerPool(pool);
            return NULL;
        }
        if (pid == 0) {
            // Child: close parent-side fds of our pipes + inherited fds of
            // any older siblings, then enter serial work loop.
            close(p2w[1]);
            close(w2p[0]);
            for (size_t j = 0; j < pool->workers.size(); ++j) {
                close(pool->workers[j]->toWorker);
                close(pool->workers[j]->fromWorker);
            }
            workerMainLoop(p2w[0], w2p[1]);
            _exit(0);
        }
        close(p2w[0]);
        close(w2p[1]);
        Worker *w = new Worker();
        w->pid = pid;
        w->toWorker = p2w[1];
        w->fromWorker = w2p[0];
        pool->workers.push_back(w);
        pool->free.push_back(static_cast<int>(pool->workers.size()) - 1);
    }
    return pool;
}

void stopWorkerPool(WorkerPool *pool) {
    if (pool == NULL) return;
    for (size_t i = 0; i < pool->workers.size(); ++i) {
        close(pool->workers[i]->toWorker); // EOF triggers worker exit
    }
    for (size_t i = 0; i < pool->workers.size(); ++i) {
        int status = 0;
        while (waitpid(pool->workers[i]->pid, &status, 0) < 0 && errno == EINTR) {}
        close(pool->workers[i]->fromWorker);
        delete pool->workers[i];
    }
    delete pool;
}

bool buildCmFromStockholmText(WorkerPool *pool,
                              const std::string &stockholmText,
                              std::string &cmText,
                              std::string &error,
                              bool calibrate) {
    if (pool == NULL || pool->workers.empty()) {
        error = "no infernal worker pool available";
        return false;
    }

    int idx = -1;
    {
        std::unique_lock<std::mutex> lk(pool->mu);
        pool->cv.wait(lk, [&]{ return !pool->free.empty(); });
        idx = pool->free.front();
        pool->free.pop_front();
    }
    Worker *w = pool->workers[idx];

    uint8_t flag = calibrate ? 1 : 0;
    uint64_t len = static_cast<uint64_t>(stockholmText.size());
    uint8_t rFlag = 0;
    uint64_t rLen = 0;
    std::string payload;
    bool ok = false;
    {
        std::lock_guard<std::mutex> iolk(w->io);
        ok = writeAllFd(w->toWorker, &flag, 1)
          && writeAllFd(w->toWorker, &len, sizeof(len))
          && (len == 0 || writeAllFd(w->toWorker, stockholmText.data(), stockholmText.size()));
        if (ok) {
            ok = readAllFd(w->fromWorker, &rFlag, 1)
              && readAllFd(w->fromWorker, &rLen, sizeof(rLen));
        }
        if (ok && rLen > 0) {
            if (rLen > (1ULL << 30)) {
                ok = false;
            } else {
                payload.resize(static_cast<size_t>(rLen));
                ok = readAllFd(w->fromWorker, &payload[0], payload.size());
            }
        }
    }

    {
        std::unique_lock<std::mutex> lk(pool->mu);
        pool->free.push_back(idx);
    }
    pool->cv.notify_one();

    if (!ok) {
        error = "infernal worker IPC failure";
        return false;
    }
    if (rFlag == 1) {
        cmText = std::move(payload);
        return true;
    }
    error = std::move(payload);
    return false;
}

} // namespace InfernalBridge
