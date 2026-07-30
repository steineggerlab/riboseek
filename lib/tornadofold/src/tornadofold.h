// tornadofold — wavefront (anti-diagonal) SIMD folder with exact Turner 2004 energy.
// Matrices stored diagonal-by-diagonal: D[s][i] holds cell (i, i+s). The O(N^3)
// multibranch bifurcation FMbif(i,j)=min_k FM(i,k)+FM1(k+1,j) becomes, per split
// span, a unit-stride vector add + min over i (NEON, 4x int32).
#pragma once
#include "energy.h"
#include <vector>
#include <string>
#include <algorithm>
#if (defined(__aarch64__) || defined(__ARM_NEON)) && !defined(DISABLE_NEON)
#include <arm_neon.h>
#define HAVE_NEON 1
#endif
#if defined(__wasm_simd128__) && !defined(DISABLE_NEON)
#include <wasm_simd128.h>
#define HAVE_WASM_SIMD 1
#endif
#if defined(__AVX512F__) && !defined(DISABLE_NEON)
#include <immintrin.h>
#define HAVE_AVX512 1
#endif
#if defined(__AVX2__) && !defined(DISABLE_NEON)
#include <immintrin.h>
#define HAVE_AVX2 1
#endif
#if defined(__SSE2__) && !defined(DISABLE_NEON)
#include <emmintrin.h>
#define HAVE_SSE2 1
#endif

namespace tornadofold {
using en::INF;

// Lower bound on any single internal-loop closure energy: the most stabilizing
// case is the best nearest-neighbour stack (-3.40 kcal/mol = -340 in centikcal
// units) across the Turner 2004 tables (stack/int11/int21/int22/general all
// >= -340). Used for an exact prune of the O(MAXLOOP^2) internal-loop scan.
static const int SINGLE_LOOP_LB = -340;

#ifdef HAVE_SSE2
// SSE2 has no signed 32-bit min (that is SSE4.1); emulate with compare + select.
static inline __m128i sse2_min_epi32(__m128i a, __m128i b) {
    __m128i lt = _mm_cmpgt_epi32(b, a); // 0xFFFFFFFF where a < b
    return _mm_or_si128(_mm_and_si128(lt, a), _mm_andnot_si128(lt, b));
}
#endif

struct TornadoFold {
    int n = 0;
    std::vector<int> b;
    std::string seq;
    en::EM em;
    // Flat span-major storage: cell (i, i+s) lives at off[s] + i. Replaces
    // vector<vector> to drop per-row heap overhead and the double-indirection
    // in the internal-loop V gather.
    std::vector<int32_t> Vf, FMf, FM1f;
    // FMbif is not stored as a full O(N^2) matrix: each span's bifurcation row is
    // consumed by FM (same span) and by V two spans later, so a 3-row rolling
    // buffer (cycled by pointer) is enough. Traceback recomputes it from FM/FM1.
    std::vector<int32_t> bifBuf;     // three n-length rows, cycled by pointer
    int32_t *bifCur, *bifP1, *bifP2; // FMbif rows for span s, s-1, s-2
    std::vector<int> off;            // off[s] = flat start of span s; size n+1
    std::vector<int32_t> Erow;
    static const int PAD = 8;

    int pt(int i, int j) const {
        return em.pt(i, j);
    }
    int eHairpin(int i, int j) const {
        return em.eHairpin(i, j);
    }
    int eIntLoop(int i, int j, int p, int q) const {
        return em.eIntLoop(i, j, p, q);
    }
    int mlStem(int p, int q) const {
        return em.mlStem(p, q);
    }
    int extStem(int k, int j) const {
        return em.extStem(k, j);
    }
    int mlClose(int i, int j) const {
        return em.mlClose(i, j);
    }

    void bifKernel(int s) {
        int m = n - s;
        int32_t* out = bifCur;
        for (int i = 0; i < m; ++i) {
            out[i] = INF;
        }
        for (int a = 0; a < s; ++a) {
            int bsp = s - 1 - a;
            const int32_t* pa = FMf.data() + off[a];
            const int32_t* pb = FM1f.data() + off[bsp] + (a + 1);
            int i = 0;
#ifdef HAVE_NEON
            for (; i + 4 <= m; i += 4) {
                int32x4_t va = vld1q_s32(pa + i), vb = vld1q_s32(pb + i);
                int32x4_t sum = vaddq_s32(va, vb);
                sum = vminq_s32(sum, vdupq_n_s32(2 * INF));
                int32x4_t cur = vld1q_s32(out + i);
                vst1q_s32(out + i, vminq_s32(cur, sum));
            }
#elif defined(HAVE_WASM_SIMD)
            for (; i + 4 <= m; i += 4) {
                v128_t va = wasm_v128_load(pa + i), vb = wasm_v128_load(pb + i);
                v128_t sum = wasm_i32x4_add(va, vb);
                sum = wasm_i32x4_min(sum, wasm_i32x4_splat(2 * INF));
                v128_t cur = wasm_v128_load(out + i);
                wasm_v128_store(out + i, wasm_i32x4_min(cur, sum));
            }
#elif defined(HAVE_AVX512)
            for (; i + 16 <= m; i += 16) {
                __m512i va = _mm512_loadu_si512((const void*)(pa + i));
                __m512i vb = _mm512_loadu_si512((const void*)(pb + i));
                __m512i sum = _mm512_add_epi32(va, vb);
                sum = _mm512_min_epi32(sum, _mm512_set1_epi32(2 * INF));
                __m512i cur = _mm512_loadu_si512((const void*)(out + i));
                _mm512_storeu_si512((void*)(out + i), _mm512_min_epi32(cur, sum));
            }
#elif defined(HAVE_AVX2)
            for (; i + 8 <= m; i += 8) {
                __m256i va = _mm256_loadu_si256((const __m256i*)(pa + i));
                __m256i vb = _mm256_loadu_si256((const __m256i*)(pb + i));
                __m256i sum = _mm256_add_epi32(va, vb);
                sum = _mm256_min_epi32(sum, _mm256_set1_epi32(2 * INF));
                __m256i cur = _mm256_loadu_si256((const __m256i*)(out + i));
                _mm256_storeu_si256((__m256i*)(out + i), _mm256_min_epi32(cur, sum));
            }
#elif defined(HAVE_SSE2)
            for (; i + 4 <= m; i += 4) {
                __m128i va = _mm_loadu_si128((const __m128i*)(pa + i));
                __m128i vb = _mm_loadu_si128((const __m128i*)(pb + i));
                __m128i sum = _mm_add_epi32(va, vb);
                sum = sse2_min_epi32(sum, _mm_set1_epi32(2 * INF));
                __m128i cur = _mm_loadu_si128((const __m128i*)(out + i));
                _mm_storeu_si128((__m128i*)(out + i), sse2_min_epi32(cur, sum));
            }
#endif
            // scalar tail (SIMD remainder), or the whole span when SIMD is off
            for (; i < m; ++i) {
                int va = pa[i], vb = pb[i];
                if (va >= INF || vb >= INF) {
                    continue;
                }
                int sm = va + vb;
                if (sm < out[i]) {
                    out[i] = sm;
                }
            }
        }
    }

    int fold(const std::string& str) {
        seq = str;
        n = (int)str.size();
        b.resize(n);
        for (int i = 0; i < n; ++i) {
            char c = str[i];
            b[i] = (c == 'A')   ? 0
                   : (c == 'C') ? 1
                   : (c == 'G') ? 2
                                : ((c == 'U' || c == 'T') ? 3 : 4);
        }
        em.set(seq, b.data(), n);
        off.resize(n + 1);
        off[0] = 0;
        for (int s = 0; s < n; ++s) {
            off[s + 1] = off[s] + (n - s);
        }
        const int tot = off[n];
        Vf.assign(tot + PAD, INF);
        FMf.assign(tot + PAD, INF);
        FM1f.assign(tot + PAD, INF);
        bifBuf.assign(3 * (size_t)n + PAD, INF);
        bifCur = bifBuf.data();
        bifP1 = bifBuf.data() + n;
        bifP2 = bifBuf.data() + 2 * (size_t)n;
        const int MLB = t2004::ML_BASE;
        for (int s = 1; s < n; ++s) {
            int m = n - s;
            for (int i = 0; i < m; ++i) {
                int j = i + s, p = pt(i, j), best = INF;
                if (p && s >= 4) {
                    best = eHairpin(i, j);
                    en::EM::IntPre pre = em.intPre(i, j); // (i,j)-invariant, hoisted
                    int maxp = std::min(j - 1, i + 31);
                    for (int a = i + 1; a <= maxp; ++a) {
                        int n1 = a - i - 1;
                        int minq = std::max(a + 1, j - 1 - (30 - n1));
                        for (int c = j - 1; c >= minq; --c) {
                            // minq bounds n1+n2 <= 30, so no explicit MAXLOOP check needed
                            int v2 = Vf[off[c - a] + a];
                            if (v2 >= INF) {
                                continue;
                            }
                            // Exact lower-bound prune: eIntLoop >= SINGLE_LOOP_LB,
                            // so if v2 + SINGLE_LOOP_LB already >= best this enclosed
                            // pair cannot improve V(i,j); skip its energy evaluation.
                            if (v2 + SINGLE_LOOP_LB >= best) {
                                continue;
                            }
                            int e = em.eIntLoop(pre, i, j, a, c);
                            if (e < INF && e + v2 < best) {
                                best = e + v2;
                            }
                        }
                    }
                    if (s - 2 >= 1) {
                        int bif = bifP2[i + 1];
                        if (bif < INF) {
                            int cl = mlClose(i, j);
                            if (bif + cl < best) {
                                best = bif + cl;
                            }
                        }
                    }
                }
                Vf[off[s] + i] = best;
            }
            for (int i = 0; i < m; ++i) {
                int j = i + s, f1 = INF;
                if (pt(i, j) && Vf[off[s] + i] < INF) {
                    f1 = Vf[off[s] + i] + mlStem(i, j);
                }
                int ext = FM1f[off[s - 1] + i];
                if (ext < INF && ext + MLB < f1) {
                    f1 = ext + MLB;
                }
                FM1f[off[s] + i] = f1;
            }
            bifKernel(s);
            for (int i = 0; i < m; ++i) {
                int fm = FM1f[off[s] + i];
                int u = FMf[off[s - 1] + (i + 1)];
                if (u < INF && u + MLB < fm) {
                    fm = u + MLB;
                }
                int bv = bifCur[i];
                if (bv < fm) {
                    fm = bv;
                }
                FMf[off[s] + i] = fm;
            }
            // cycle rolling FMbif rows: cur -> P1 -> P2 -> recycled -> cur
            int32_t* t = bifCur;
            bifCur = bifP2;
            bifP2 = bifP1;
            bifP1 = t;
        }
        Erow.assign(n, 0);
        for (int j = 0; j < n; ++j) {
            int e = (j > 0) ? Erow[j - 1] : 0;
            for (int k = 0; k <= j; ++k) {
                int v = Vf[off[j - k] + k];
                if (v >= INF) {
                    continue;
                }
                int pre = (k > 0) ? Erow[k - 1] : 0;
                int cand = pre + v + extStem(k, j);
                if (cand < e) {
                    e = cand;
                }
            }
            Erow[j] = e;
        }
        return n ? Erow[n - 1] : 0;
    }

    int getV(int i, int j) const {
        return Vf[off[j - i] + i];
    }
    int getFM(int i, int j) const {
        return FMf[off[j - i] + i];
    }
    int getFM1(int i, int j) const {
        return FM1f[off[j - i] + i];
    }
    int getFMbif(int i, int j) const {
        int v = INF;
        for (int k = i; k < j; ++k) {
            int a = getFM(i, k), c = getFM1(k + 1, j);
            if (a < INF && c < INF && a + c < v) {
                v = a + c;
            }
        }
        return v;
    }

    std::string traceback(int) {
        std::string dot(n, '.');
        int j = n - 1;
        while (j >= 0) {
            int e = Erow[j];
            if (j > 0 && Erow[j - 1] == e) {
                --j;
                continue;
            }
            bool found = false;
            for (int k = 0; k <= j; ++k) {
                int v = getV(k, j);
                if (v >= INF) {
                    continue;
                }
                int pre = (k > 0) ? Erow[k - 1] : 0;
                if (pre + v + extStem(k, j) == e) {
                    dot[k] = '(';
                    dot[j] = ')';
                    traceV(k, j, dot);
                    j = k - 1;
                    found = true;
                    break;
                }
            }
            if (!found) {
                --j;
            }
        }
        return dot;
    }
    void traceV(int i, int j, std::string& dot) {
        int v = getV(i, j);
        if (eHairpin(i, j) == v) {
            return;
        }
        int maxp = std::min(j - 1, i + 31);
        for (int a = i + 1; a <= maxp; ++a) {
            int n1 = a - i - 1;
            int minq = std::max(a + 1, j - 1 - (30 - n1));
            for (int c = j - 1; c >= minq; --c) {
                // minq bounds n1+n2 <= 30 (see fold())
                int v2 = getV(a, c);
                if (v2 >= INF) {
                    continue;
                }
                int e = eIntLoop(i, j, a, c);
                if (e < INF && e + v2 == v) {
                    dot[a] = '(';
                    dot[c] = ')';
                    traceV(a, c, dot);
                    return;
                }
            }
        }
        if ((j - 1) - (i + 1) >= 1) {
            int cl = mlClose(i, j);
            if (getFMbif(i + 1, j - 1) + cl == v) {
                traceFMbif(i + 1, j - 1, dot);
                return;
            }
        }
    }
    void traceFMbif(int i, int j, std::string& dot) {
        int val = getFMbif(i, j);
        for (int k = i; k < j; ++k) {
            int a = getFM(i, k), c = getFM1(k + 1, j);
            if (a >= INF || c >= INF) {
                continue;
            }
            if (a + c == val) {
                traceFM(i, k, dot);
                traceFM1(k + 1, j, dot);
                return;
            }
        }
    }
    void traceFM1(int i, int j, std::string& dot) {
        int f1 = getFM1(i, j);
        if (pt(i, j) && getV(i, j) < INF && getV(i, j) + mlStem(i, j) == f1) {
            dot[i] = '(';
            dot[j] = ')';
            traceV(i, j, dot);
            return;
        }
        if (j - 1 >= i && getFM1(i, j - 1) < INF && getFM1(i, j - 1) + t2004::ML_BASE == f1) {
            traceFM1(i, j - 1, dot);
        }
    }
    void traceFM(int i, int j, std::string& dot) {
        int fm = getFM(i, j);
        if (getFM1(i, j) == fm) {
            traceFM1(i, j, dot);
            return;
        }
        if (i + 1 <= j && getFM(i + 1, j) < INF && getFM(i + 1, j) + t2004::ML_BASE == fm) {
            traceFM(i + 1, j, dot);
            return;
        }
        traceFMbif(i, j, dot);
    }
};

} // namespace tornadofold
