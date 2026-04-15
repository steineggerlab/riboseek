#include "InfernalBridge.h"

#include <cstdlib>
#include <csetjmp>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/mman.h>
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

namespace InfernalBridge {

bool isConfigured() {
    return true;
}

bool buildCmFromStockholmText(const std::string &stockholmText, std::string &cmText, std::string &error) {
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
    if (runInfernalMainInProcess(&infernal_cmbuild_main, cmbArgs) != 0) {
        error = "embedded infernal cmbuild failed";
        close(stoFd);
        unlink(cmPath.c_str());
        rmdir(dir.c_str());
        return false;
    }

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

    close(stoFd);
    unlink(cmPath.c_str());
    rmdir(dir.c_str());
    return true;
}

} // namespace InfernalBridge
