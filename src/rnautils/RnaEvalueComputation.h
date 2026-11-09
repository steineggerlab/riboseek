#ifndef RIBOSEEK_RNA_EVALUE_COMPUTATION_H
#define RIBOSEEK_RNA_EVALUE_COMPUTATION_H

// riboseek's E-value model for the dinucleotide aligners (rnaalign, the realigner and the
// fragment merge), kept on the riboseek side of the tree (2026-10-19). It extends MMseqs2's
// EvalueComputation, which stays as upstream: MMseqs modules keep using the base class and
// its ALP statistics unchanged; only riboseek code constructs this one.
//   - calibrated constants for dinucleotide scoring (sequence and profile queries)
//   - the per-query sequence model (EvalueModelSeq4.h, legacy EvalueModel.h net)
//   - the free-lambda profile model (EvalueModelProf3.h)
//   - the target GC histogram / dinucleotide background those models read (TargetGcStat.h)

#include "EvalueComputation.h"
#include "EvalueModel.h"
#include "EvalueModelProf3.h"
// The runtime implements one profile model (S1n, 2026-10-20): 6 inputs, 2 outputs (lambda, mu). A header
// of any other shape (a 7-input model with lambda_KA, or a mu-only model with a lambda table, which has no
// NOUT) would read past the input array or mis-index the outputs, so refuse it at compile time.
static_assert(EvalueModelProf3::NIN == 6 && EvalueModelProf3::NOUT == 2,
              "EvalueModelProf3.h must be a 6-input, 2-output (lambda, mu) model: export it with export_prof3.py "
              "from a s4c --free-lambda --features cpg model");
#include "EvalueModelSeq4.h"
#include <vector>
#include <string>
#include <cstdio>
#include <utility>
#include <cstdlib>
#include <cmath>
#include <algorithm>

class RnaEvalueComputation : public EvalueComputation {
public:
    RnaEvalueComputation(size_t dbResCount, BaseMatrix *subMat, int gapOpen, int gapExtend,
                         bool bothStrands, bool useRnaCorrection)
        : EvalueComputation(dbResCount, subMat, gapOpen, gapExtend, bothStrands, useRnaCorrection),
          bothStrands(bothStrands), dinucProfile(false) {}

    // Which model scores this search: the profile model (a genuine profile query, built by
    // result2profile) or the sequence model. Set once, before any E-value is computed.
    void setProfileQuery(bool profileQuery) { dinucProfile = profileQuery; }

    // E-value for the RNA/dinucleotide aligners: the per-query model (below), the only E-value
    // riboseek reports since 2026-10-20 (the constant-parameter and ALP forms were removed).
    // Per-query cache of the model's output.
    //
    // lambda and mu depend only on the query's own features and the target GC
    // bin -- NOT on the score -- so for one query they are identical across
    // every alignment it produces. Recomputing the MLP per alignment costs
    // 2181 ns against 25 ns for reusing these: 88x, measured -O2 on the shipped
    // weights with four bins.
    //
    // The cache is OWNED BY THE CALLER, one per query profile, because a single
    // EvalueComputation is shared by every thread in rnaalign and anything
    // per-query stored on the shared object would be a data race. s_profile
    // already holds exactly the three inputs, so it owns one of these.
    struct QueryEvalueCache {
        enum { MAXBIN = 8 };
        double lam[MAXBIN];
        double mu[MAXBIN];
        double w[MAXBIN];               // bin weight; for the profile model, weight * K
        int nbin;
        double kGc, kLen, kNeff, kHet;  // the query this was filled for
        unsigned long serial;           // profile model: ProfileEvalueFeatures::serial filled for
        QueryEvalueCache() : nbin(-1), kGc(-1.0), kLen(-1.0), kNeff(-1.0), kHet(-1.0), serial(0) {}
        inline bool matches(double gc, double len, double neff, double het, int nb) const {
            return nbin == nb && kGc == gc && kLen == len && kNeff == neff && kHet == het;
        }
    };

    // Features of one profile query for the profile model
    // (EvalueModelProf3.h, exported by rawscores/repro/export_prof3.py). Filled
    // once per query in RnaSmithWaterman::ssw_init. Everything is computed from
    // the FORWARD-equivalent profile: on the reverse-complement copy that
    // splitstrand stores at an odd key, each column is permuted back through the
    // dinucleotide revcomp, so both copies get the forward strand's features --
    // what the model was trained on (forward copies, E doubled for strands).
    //   S        raw profile scores, L x 16, as stored (round(8 log2(p/bg)))
    //   gc       emission GC: mean over columns of P(first base G or C)
    //   cpg/tpa  CpG and TpA observed/expected of the mean dinucleotide emission
    //   neffFrac mean Neff byte / 255
    // serial changes with every query so a per-profile cache is never reused
    // across queries with coincidentally equal features.
    // het is the SEQUENCE model's 4th input (EvalueModelSeq4.h): within-query GC sd over
    // HETW-nt windows of the forward strand, filled for plain sequence queries; -1 when unset.
    struct ProfileEvalueFeatures {
        bool valid;
        int L;
        double gc, cpg, tpa, neffFrac;
        double het;
        std::vector<float> S;
        unsigned long serial;
        ProfileEvalueFeatures() : valid(false), L(0), gc(0.5), cpg(1.0), tpa(1.0), neffFrac(0.0), het(-1.0), serial(0) {}
    };

    // s5e_seq_boot.hetsd: population sd of GC over consecutive non-overlapping windows of
    // max(10, min(HETW, n/4)) nt starting at position 0 (a short tail is dropped); 0 for
    // n < 40 or fewer than 3 windows. isGc[i] says whether nucleotide i of the FORWARD
    // strand is G or C.
    static double windowGcSd(const std::vector<char> &isGc, int hetw) {
        const int n = (int)isGc.size();
        int w = std::max(10, std::min(hetw, n / 4));
        if (n < 40) return 0.0;
        std::vector<double> g;
        for (int i = 0; i + w <= n; i += w) {
            int c = 0;
            for (int j = i; j < i + w; j++) c += isGc[j];
            g.push_back((double)c / w);
        }
        if (g.size() < 3) return 0.0;
        double m = 0.0, v = 0.0;
        for (double x : g) m += x;
        m /= g.size();
        for (double x : g) v += (x - m) * (x - m);
        return std::sqrt(v / g.size());
    }

    // Dinucleotide state helpers, profile column order A C D E F G H I K L M N P Q R S
    // = GC AA CA TA CC AC TG GT CT AG GA TC GG AT CG TT.
    static const int *stateRevcomp() {      // state of the reverse complement dinucleotide
        static const int rc[16] = {0, 15, 6, 3, 12, 7, 2, 5, 9, 8, 11, 10, 4, 13, 14, 1};
        return rc;
    }
    static const int *stateFirstNuc() {     // 0=A 1=C 2=G 3=T
        static const int f[16] = {2, 0, 1, 3, 1, 0, 3, 2, 1, 0, 2, 3, 2, 0, 1, 3};
        return f;
    }
    static const int *stateSecondNuc() {
        static const int s2[16] = {1, 0, 0, 0, 1, 1, 2, 3, 3, 2, 0, 1, 2, 3, 2, 3};
        return s2;
    }
    // state of dinucleotide (n1, n2), nucleotides coded 0=A 1=C 2=G 3=T
    static int dinucState(int n1, int n2) {
        static const int t[16] = {1, 5, 9, 13,  2, 4, 14, 8,  10, 0, 12, 7,  3, 11, 6, 15};
        return t[n1 * 4 + n2];
    }

    // Fill f from L forward-equivalent score columns (16 each). Mirrors
    // profgc.decode ('fixed'), profgc.dinuc_oe and the emission GC in s5d_boot.py.
    static void fillProfileFeatures(ProfileEvalueFeatures &f, int L, double neffFrac) {
        f.L = L; f.neffFrac = neffFrac;
        const int *n1 = stateFirstNuc(), *n2 = stateSecondNuc();
        double m[16] = {0}, gc = 0.0;
        for (int j = 0; j < L; j++) {
            double p[16], sum = 0.0;
            for (int a = 0; a < 16; a++) {
                p[a] = EvalueModelProf3::DECODE_BG[a] * std::exp2(f.S[j * 16 + a] / 8.0);
                sum += p[a];
            }
            double g = 0.0;
            for (int a = 0; a < 16; a++) {
                p[a] /= sum; m[a] += p[a];
                if (n1[a] == 1 || n1[a] == 2) g += p[a];
            }
            gc += g;
        }
        if (L > 0) { gc /= L; for (int a = 0; a < 16; a++) m[a] /= L; }
        double f1[4] = {0}, f2[4] = {0};
        for (int a = 0; a < 16; a++) { f1[n1[a]] += m[a]; f2[n2[a]] += m[a]; }
        f.gc = gc;
        f.cpg = m[dinucState(1, 2)] / std::max(f1[1] * f2[2], 1e-300);
        f.tpa = m[dinucState(3, 0)] / std::max(f1[3] * f2[0], 1e-300);
        f.valid = (L > 0);
    }

    // S1n profile model (export_prof3.py on a s4c --free-lambda --features cpg model): no lambda table, no lambda_KA; the n -> 32 -> 32 -> 2 net
    // predicts (lambda, mu) per profile and target bin. Mean over seeds of each output, as s5d_boot.py does.
    static void prof3LamMu(const double *x, double &lam, double &mu) {
        using namespace EvalueModelProf3;
        double sl = 0.0, sm = 0.0;
        for (int s = 0; s < NSEED; s++) {
            const double *xm = XM + s * NIN, *xs = XS + s * NIN;
            const double *w0 = W0 + s * NIN * HID, *b0 = B0 + s * HID;
            const double *w1 = W1 + s * HID * HID, *b1 = B1 + s * HID;
            const double *w2 = W2 + s * HID * 2;
            double h0[HID], h1[HID];
            for (int j = 0; j < HID; j++) {
                double a = b0[j];
                for (int i = 0; i < NIN; i++) a += ((x[i] - xm[i]) / xs[i]) * w0[i * HID + j];
                h0[j] = a > 0 ? a : 0;
            }
            for (int j = 0; j < HID; j++) {
                double a = b1[j];
                for (int i = 0; i < HID; i++) a += h0[i] * w1[i * HID + j];
                h1[j] = a > 0 ? a : 0;
            }
            double o0 = B2[s * 2], o1 = B2[s * 2 + 1];
            for (int i = 0; i < HID; i++) { o0 += h1[i] * w2[i * 2]; o1 += h1[i] * w2[i * 2 + 1]; }
            sl += o0 * YS[s * 2] + YM[s * 2];
            sm += o1 * YS[s * 2 + 1] + YM[s * 2 + 1];
        }
        lam = sl / NSEED; mu = sm / NSEED;
    }
    // diag/kfit.py kpoly . KCOEF
    static double prof3K(double gc, double logLen, double tgc) {
        const double x0 = gc - 0.5, x1 = logLen - 2.5, x2 = tgc - 0.5;
        const double t[10] = {1.0, x0, x1, x2, x0 * x0, x1 * x1, x2 * x2, x0 * x1, x0 * x2, x1 * x2};
        double e = 0.0;
        for (int i = 0; i < 10; i++) e += EvalueModelProf3::KCOEF[i] * t[i];
        return std::exp(e);
    }

    // Per-query Gumbel parameters from the fitted MLP. Inputs are the query's own
    // composition and length plus one target-GC bin; E is summed over the bins.
    // See rawscores/repro/README.md. lambda and mu here describe ONE strand, so the
    // both-strands factor is applied by the caller exactly as for the constant form.
    static void mlpPredict(const double *x, int nin, const double *W0, const double *b0,
                           const double *W1, const double *b1, const double *W2,
                           const double *b2, const double *xm, const double *xs,
                           const double *ym, const double *ys, double &lam, double &mu) {
        double h0[32], h1[32], o[2];
        for (int j = 0; j < 32; j++) {
            double a = b0[j];
            for (int i = 0; i < nin; i++) a += ((x[i]-xm[i])/xs[i]) * W0[i*32+j];
            h0[j] = a > 0 ? a : 0;
        }
        for (int j = 0; j < 32; j++) {
            double a = b1[j];
            for (int i = 0; i < 32; i++) a += h0[i] * W1[i*32+j];
            h1[j] = a > 0 ? a : 0;
        }
        for (int j = 0; j < 2; j++) {
            double a = b2[j];
            for (int i = 0; i < 32; i++) a += h1[i] * W2[i*2+j];
            o[j] = a * ys[j] + ym[j];
        }
        lam = o[0]; mu = o[1];
    }

    // The database summary the model needs: (mean GC, weight) per bin. The weight
    // is the bin's RESIDUES divided by EvalueModel::R0_TARGET_LEN, the mean target
    // length of the nulls the model was fitted on -- not the bin's sequence count.
    // mu absorbed that length scale during fitting, so a count would assume every
    // database has targets of about that length; on a database of 3000-6000 nt
    // targets that under-reports E by 4.5-7.9x, measured. Four bins is enough.
    void setTargetGcHistogram(const std::vector<std::pair<double, double> > &h) { gcHist = h; }
    bool hasTargetGcHistogram() const { return gcHist.empty() == false; }

    // E-value from the profile model. Same database summary and
    // the same both-strands factor as nnEvalue; K and the residue weight are folded
    // into one per-bin weight in the cache.
    inline double prof3Evalue(double score, const ProfileEvalueFeatures &f, QueryEvalueCache *cache) const {
        const int nb = (int)gcHist.size();
        const double strand = bothStrands ? 2.0 : 1.0;
        QueryEvalueCache local;
        QueryEvalueCache *c = (cache != NULL && nb <= QueryEvalueCache::MAXBIN) ? cache : &local;
        if (!(c->nbin == nb && c->serial == f.serial && f.serial != 0)) {
            const double logLen = std::log10(std::max((double)f.L, 2.0));
            for (int b = 0; b < nb && b < QueryEvalueCache::MAXBIN; b++) {
                const double tgc = gcHist[b].first;
                const double x[6] = {f.gc, logLen, tgc, f.neffFrac, f.cpg, f.tpa};
                double lam, mu; prof3LamMu(x, lam, mu);
                c->lam[b] = std::max(lam, 1e-3); c->mu[b] = mu;
                c->w[b] = gcHist[b].second * prof3K(f.gc, logLen, tgc);
            }
            c->nbin = nb; c->serial = f.serial;
        }
        double total = 0.0;
        for (int b = 0; b < nb && b < QueryEvalueCache::MAXBIN; b++) {
            total += c->w[b] * -std::expm1(-std::exp(-c->lam[b] * (score - c->mu[b])));
        }
        return total * strand;
    }

    // Query features are ARGUMENTS, not object state: one EvalueComputation is
    // shared by every thread in rnaalign (declared outside the omp region and
    // handed to each per-thread RnaMatcher), so anything per-query stored on it
    // would be a data race. gcHist is written once before the region and only
    // read inside it, which is safe.
    // mean over seeds of the 4 -> 32 -> 32 -> 2 sequence network (EvalueModelSeq4.h)
    static void seq4Predict(const double *x, double &lam, double &mu) {
        using namespace EvalueModelSeq4;
        double sl = 0.0, sm = 0.0;
        for (int s = 0; s < NSEED; s++) {
            const double *xm = XM + s * NIN, *xs = XS + s * NIN;
            const double *w0 = W0 + s * NIN * HID, *b0 = B0 + s * HID;
            const double *w1 = W1 + s * HID * HID, *b1 = B1 + s * HID;
            const double *w2 = W2 + s * HID * 2, *b2 = B2 + s * 2;
            double h0[HID], h1[HID];
            for (int j = 0; j < HID; j++) {
                double a = b0[j];
                for (int i = 0; i < NIN; i++) a += ((x[i] - xm[i]) / xs[i]) * w0[i * HID + j];
                h0[j] = a > 0 ? a : 0;
            }
            for (int j = 0; j < HID; j++) {
                double a = b1[j];
                for (int i = 0; i < HID; i++) a += h0[i] * w1[i * HID + j];
                h1[j] = a > 0 ? a : 0;
            }
            double o[2];
            for (int j = 0; j < 2; j++) {
                double a = b2[j];
                for (int i = 0; i < HID; i++) a += h1[i] * w2[i * 2 + j];
                o[j] = a * YS[s * 2 + j] + YM[s * 2 + j];
            }
            sl += o[0]; sm += o[1];
        }
        lam = sl / NSEED; mu = sm / NSEED;
    }
    static double seq4K(double gc, double logLen, double tgc) {
        const double x0 = gc - 0.5, x1 = logLen - 2.5, x2 = tgc - 0.5;
        const double t[10] = {1.0, x0, x1, x2, x0 * x0, x1 * x1, x2 * x2, x0 * x1, x0 * x2, x1 * x2};
        double e = 0.0;
        for (int i = 0; i < 10; i++) e += EvalueModelSeq4::KCOEF[i] * t[i];
        return std::exp(e);
    }

    // het < 0 means the within-query GC sd is unknown; sequences then use the legacy
    // 3-input net, which needs no 4th input.
    inline double nnEvalue(double score, double seqLength, double queryGc,
                           double neffFrac, QueryEvalueCache *cache = NULL,
                           double het = -1.0) const {
        const int nb = (int)gcHist.size();
        const double strand = bothStrands ? 2.0 : 1.0;

        // Fast path: this query's (lambda, mu) are already in hand.
        if (cache != NULL && cache->matches(queryGc, seqLength, neffFrac, het, nb)) {
            double total = 0.0;
            for (int b = 0; b < nb; b++) {
                total += cache->w[b] * -std::expm1(-std::exp(-cache->lam[b] * (score - cache->mu[b])));
            }
            return total * strand;
        }

        const bool prof = dinucProfile;
        const int nin = prof ? EvalueModel::PROF_NIN : EvalueModel::SEQ_NIN;
        const double logLen = std::log10(std::max(seqLength, 2.0));
        const bool fill = (cache != NULL && nb <= QueryEvalueCache::MAXBIN);
        const bool seq4 = !prof && het >= 0.0;
        double total = 0.0;
        for (int b = 0; b < nb; b++) {
            double x[4] = { queryGc, logLen, gcHist[b].first, neffFrac };  // NOLINT
            double lam, mu;
            if (seq4) {
                // bin weight is residues / R0_TARGET_LEN; the model wants residues / its own R0
                x[3] = het;
                seq4Predict(x, lam, mu);
                if (lam < 1e-3) lam = 1e-3;
                const double w = gcHist[b].second * (EvalueModel::R0_TARGET_LEN / EvalueModelSeq4::R0)
                               * seq4K(queryGc, logLen, gcHist[b].first);
                if (fill) { cache->lam[b] = lam; cache->mu[b] = mu; cache->w[b] = w; }
                total += w * -std::expm1(-std::exp(-lam * (score - mu)));
                continue;
            }
            if (prof) {
                mlpPredict(x, nin, EvalueModel::PROF_W0, EvalueModel::PROF_b0,
                           EvalueModel::PROF_W1, EvalueModel::PROF_b1,
                           EvalueModel::PROF_W2, EvalueModel::PROF_b2,
                           EvalueModel::PROF_xm, EvalueModel::PROF_xs,
                           EvalueModel::PROF_ym, EvalueModel::PROF_ys, lam, mu);
            } else {
                mlpPredict(x, nin, EvalueModel::SEQ_W0, EvalueModel::SEQ_b0,
                           EvalueModel::SEQ_W1, EvalueModel::SEQ_b1,
                           EvalueModel::SEQ_W2, EvalueModel::SEQ_b2,
                           EvalueModel::SEQ_xm, EvalueModel::SEQ_xs,
                           EvalueModel::SEQ_ym, EvalueModel::SEQ_ys, lam, mu);
            }
            if (lam < 1e-3) lam = 1e-3;
            if (fill) { cache->lam[b] = lam; cache->mu[b] = mu; cache->w[b] = gcHist[b].second; }
            // P(score >= S) for a Gumbel, written as -expm1 for small-argument accuracy
            total += gcHist[b].second * -std::expm1(-std::exp(-lam * (score - mu)));
        }
        if (fill) {
            cache->nbin = nb; cache->kGc = queryGc; cache->kLen = seqLength; cache->kNeff = neffFrac; cache->kHet = het;
        }
        return total * strand;
    }

    // Every caller passes the query's features (RnaSmithWaterman fills them in ssw_init), and
    // the target GC histogram must have been set (rnaalign and the realigner stop if it cannot be
    // built), so there is no fallback form.
    inline double computeRnaEvalue(double score, double seqLength, double queryGc, double neffFrac,
                                   QueryEvalueCache *cache, const ProfileEvalueFeatures *pf) {
        if (dinucProfile && pf != NULL && pf->valid) {
            return prof3Evalue(score, *pf, cache);
        }
        return nnEvalue(score, seqLength, queryGc, neffFrac, cache, pf != NULL ? pf->het : -1.0);
    }

private:
    std::vector<std::pair<double, double> > gcHist;  // read-only inside omp regions
    // Copy of the base class's (private) strand flag; same value as the base was constructed with.
    const bool bothStrands;
    bool dinucProfile;
};

#endif // RIBOSEEK_RNA_EVALUE_COMPUTATION_H
