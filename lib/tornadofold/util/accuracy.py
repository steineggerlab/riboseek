#!/usr/bin/env python3
# Accuracy benchmark on ArchiveII: compare tornadofold vs RNAfold -d2 against
# reference structures. Reports micro- and macro-averaged sensitivity/PPV/F1.
#
# Point it at a directory of .ct or .bpseq reference files (the ArchiveII
# tarballs ship .bpseq):
#   ARCHIVEII=path/to/archiveII TORNADOFOLD=./tornadofold RNAFOLD=RNAfold \
#       python3 util/accuracy.py [maxlen]
import os, sys, glob, subprocess, re

DATADIR = os.environ.get("ARCHIVEII",
                         os.path.join(os.path.dirname(__file__), "..", "archiveII"))
TORNADOFOLD = os.environ.get("TORNADOFOLD", os.path.join(os.path.dirname(__file__), "..", "tornadofold"))
# RNAfold comparison is opt-in: only run it when RNAFOLD points at a binary.
RNAFOLD = os.environ.get("RNAFOLD")

def parse_bpseq(path):
    seq, pairs = [], {}
    for ln in open(path):
        t = ln.split()
        if len(t) < 3 or not t[0].isdigit(): continue
        i = int(t[0]); base = t[1].upper().replace('T','U'); j = int(t[2])
        seq.append(base if base in 'ACGU' else 'N')
        if j > i: pairs[i-1] = j-1
    return "".join(seq), pairs

def parse_ct(path):
    seq, pairs = [], {}
    with open(path) as fh:
        lines = fh.read().splitlines()
    # header line: "<len>  title"; body lines: idx base prev next pair natidx
    body = []
    started = False
    for ln in lines:
        t = ln.split()
        if not t: continue
        if not started:
            # first numeric-leading line is header
            if re.match(r'^\d+', t[0]):
                started = True
            continue
        if len(t) >= 5 and t[0].isdigit():
            body.append(t)
    for t in body:
        i = int(t[0]); base = t[1].upper().replace('T','U')
        j = int(t[4])
        seq.append(base if base in 'ACGU' else 'N')
        if j > i:
            pairs[i-1] = j-1
    return "".join(seq), pairs

def db_to_pairs(db):
    st, pr = [], {}
    for i,c in enumerate(db):
        if c=='(' : st.append(i)
        elif c==')':
            if st:
                k=st.pop(); pr[k]=i
    return pr

def load_dataset(maxlen=None):
    recs=[]
    files = sorted(glob.glob(os.path.join(DATADIR,"*.ct"))
                   + glob.glob(os.path.join(DATADIR,"*.bpseq")))
    for f in files:
        seq,pairs = parse_bpseq(f) if f.endswith(".bpseq") else parse_ct(f)
        if len(seq)==0: continue
        if maxlen and len(seq)>maxlen: continue
        recs.append((os.path.basename(f), seq, pairs))
    return recs

def run_tool_fasta(recs, cmd, kind):
    fa = "\n".join(f">{i}\n{seq}" for i,(_,seq,_) in enumerate(recs))+"\n"
    p = subprocess.run(cmd, input=fa, capture_output=True, text=True)
    out = p.stdout
    # extract dot-bracket lines in order
    structs=[]
    for ln in out.splitlines():
        s=ln.strip()
        # a structure line: first token is only structure chars
        tok=s.split()[0] if s else ""
        if tok and set(tok) <= set("().[]{}<>") and len(tok)>0:
            structs.append(tok)
    return structs

def score(recs, structs):
    # returns (micro dict, macro dict)
    TP=FP=FN=0
    mf1=ms=mp=0.0; nm=0
    per=[]
    for (name,seq,ref),db in zip(recs,structs):
        pred = db_to_pairs(db)
        refset = set((k,v) for k,v in ref.items())
        predset = set((k,v) for k,v in pred.items())
        tp=len(refset & predset); fp=len(predset-refset); fn=len(refset-predset)
        TP+=tp; FP+=fp; FN+=fn
        sens = tp/(tp+fn) if (tp+fn) else (1.0 if not predset else 0.0)
        ppv  = tp/(tp+fp) if (tp+fp) else (1.0 if not refset else 0.0)
        f1   = 2*sens*ppv/(sens+ppv) if (sens+ppv) else 0.0
        ms+=sens; mp+=ppv; mf1+=f1; nm+=1
        per.append((name,len(seq),sens,ppv,f1))
    micro_s = TP/(TP+FN) if (TP+FN) else 0
    micro_p = TP/(TP+FP) if (TP+FP) else 0
    micro_f = 2*micro_s*micro_p/(micro_s+micro_p) if (micro_s+micro_p) else 0
    macro = (ms/nm, mp/nm, mf1/nm)
    return (micro_s,micro_p,micro_f), macro, per

if __name__=="__main__":
    maxlen = int(sys.argv[1]) if len(sys.argv)>1 else None
    recs = load_dataset(maxlen)
    print(f"loaded {len(recs)} structures"+(f" (<= {maxlen} nt)" if maxlen else ""), flush=True)
    import time
    tools = [("tornadofold",[TORNADOFOLD])]
    if RNAFOLD:
        tools.append(("RNAfold -d2",[RNAFOLD,"--noPS","-d2"]))
    for label,cmd in tools:
        t0=time.time()
        structs = run_tool_fasta(recs, cmd, label)
        dt=time.time()-t0
        if len(structs)!=len(recs):
            print(f"  WARN {label}: got {len(structs)} structs for {len(recs)} recs")
        micro,macro,per = score(recs, structs[:len(recs)])
        print(f"\n{label}  ({dt:.1f}s)")
        print(f"  micro  sens={micro[0]:.3f}  ppv={micro[1]:.3f}  F1={micro[2]:.3f}")
        print(f"  macro  sens={macro[0]:.3f}  ppv={macro[1]:.3f}  F1={macro[2]:.3f}")
        # per-family macro F1
        fam={}
        for name,L,s,p,f in per:
            fk=name.split('_')[0]
            fam.setdefault(fk,[]).append(f)
        fl=" ".join(f"{k}={sum(v)/len(v):.2f}" for k,v in sorted(fam.items()))
        print(f"  family macroF1: {fl}", flush=True)
