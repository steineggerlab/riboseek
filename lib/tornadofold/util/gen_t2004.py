#!/usr/bin/env python3
"""Generate src/t2004.h directly from the NNDB Turner 2004 .dg bundle.

Reads the primary turner_2004 .dg source files (never a pre-baked header),
builds the nearest-neighbour tables, re-indexes them into this folder's compact
layout, and emits src/t2004.h. Also diffs every folding-reachable value against
the in-tree tables so any divergence is surfaced, not hidden.

  python3 util/gen_t2004.py --turner-source util/turner_2004.zip --out src/t2004.h
  python3 util/gen_t2004.py --turner-source util/turner_2004.zip --check-only

Layout: pair index = ViennaRNA pair type - 1 (CG..NN = 0..6); mismatch/dangle
base index bidx(N)=0, A..U=1..4; int22 canonical-only [6][6][4][4][4][4].
lxc defaults to ViennaRNA's precise 107.856 (the .dg file rounds it to 107.9);
pass --lxc 107.9 for pure .dg provenance (folds identically at real loop sizes).
"""
import argparse
import re
import shutil
import tempfile
import zipfile
from decimal import Decimal, ROUND_HALF_UP
from pathlib import Path

INF = 1_000_000
VIE_INF = 10_000_000
NT = "ACGU"
NT_INDEX = {b: i for i, b in enumerate(NT)}
CANON = ["CG", "GC", "GU", "UG", "AU", "UA"]  # ViennaRNA pairs 1..6
VALUE_RE = re.compile(r"^[+-]?(?:\d+\.?\d*|\.\d+|\.)$")


# --- .dg parsing (from the NNDB primary sources) ----------------------------

def scale(tok):
    if tok == ".":
        return VIE_INF
    return int((Decimal(tok) * 100).quantize(Decimal(1), rounding=ROUND_HALF_UP))


def reverse_pair(p):
    return p[::-1]


def read_grid_blocks(path, ncols):
    blocks, labels, rows = [], [], []
    for raw in path.read_text().splitlines():
        stripped = raw.strip()
        if not stripped or stripped.startswith("#"):
            if rows:
                blocks.append((labels, rows)); labels, rows = [], []
            continue
        values = [t for t in raw.split() if VALUE_RE.match(t)]
        if len(values) == ncols:
            rows.append(values)
        else:
            if rows:
                blocks.append((labels, rows)); labels, rows = [], []
            labels.append(stripped)
    if rows:
        blocks.append((labels, rows))
    return blocks


def scaled_grid(rows):
    return [[scale(t) for t in row] for row in rows]


def read_number_groups(path):
    groups, current = [], []
    for raw in path.read_text().splitlines():
        stripped = raw.strip()
        if not stripped or stripped.startswith("#"):
            if current:
                groups.append(current); current = []
            continue
        current.extend(stripped.split())
    if current:
        groups.append(current)
    return groups


def two_base_label(label):
    toks = label.split()
    return len(toks) == 2 and all(b in NT for t in toks for b in t)


def parse_pair_grids(path, key_from_labels):
    table = {}
    for labels, rows in read_grid_blocks(path, 4):
        table[key_from_labels(labels)] = scaled_grid(rows)
    return table


def stack_key(labels):
    top = next(l for l in labels if "X" in l)
    bottom = next(l for l in labels if "Y" in l)
    return top[top.index("X") - 1] + bottom[bottom.index("Y") - 1]


def parse_stack(path):
    return parse_pair_grids(path, stack_key)


def parse_mismatch(path):
    return parse_pair_grids(path, stack_key)


def parse_dangles(path):
    blocks = read_grid_blocks(path, 4)
    dangle3, dangle5 = {}, {}
    for index, (labels, rows) in enumerate(blocks):
        top_base = next(l for l in labels if "X" in l).replace("X", "")[0]
        plain_base = next(l for l in labels if "X" not in l and l in NT)
        values = [scale(t) for t in rows[0]]
        pair = plain_base + top_base
        (dangle3 if index < 16 else dangle5)[pair if index < 16 else reverse_pair(pair)] = values
    return dangle3, dangle5


def parse_specials(path):
    entries = []
    for raw in path.read_text().splitlines():
        s = raw.strip()
        if not s or s.startswith("#"):
            continue
        seq, energy = s.split()
        entries.append((seq, scale(energy)))
    entries.sort()
    return entries


def parse_loops(path):
    internal, bulge, hairpin = [VIE_INF], [VIE_INF], [VIE_INF]
    for raw in path.read_text().splitlines():
        s = raw.strip()
        if not s or s.startswith("#"):
            continue
        _, i_tok, b_tok, h_tok = s.split()
        internal.append(scale(i_tok)); bulge.append(scale(b_tok)); hairpin.append(scale(h_tok))
    return internal, bulge, hairpin


def parse_int11(path):
    table = {}
    for labels, rows in read_grid_blocks(path, 4):
        top, bottom = (l.split() for l in labels if two_base_label(l))
        table[(top[0] + bottom[0], top[1] + bottom[1])] = scaled_grid(rows)
    return table


def parse_int21(path):
    table = {}
    for labels, rows in read_grid_blocks(path, 4):
        pair_lines = [l.split() for l in labels if two_base_label(l)]
        top, bottom = pair_lines[0], pair_lines[1]
        panel_label = next(l for l in labels if "Y" in l)
        panel = panel_label[panel_label.index("Y") + 1]
        table[(top[0] + bottom[0], top[1] + bottom[1], panel)] = scaled_grid(rows)
    return table


def parse_int22(path):
    table = {}
    for labels, rows in read_grid_blocks(path, 16):
        top = next(l for l in labels if "X1" in l).split()
        bottom = next(l for l in labels if "X2" in l).split()
        table[(top[0] + bottom[0], top[-1] + bottom[-1])] = scaled_grid(rows)
    return table


# --- build ViennaRNA-layout [8]/[5] nested arrays ---------------------------

def pair_matrix(value_of):
    m = []
    for i in range(8):
        row = []
        for j in range(8):
            if i == 0 or j == 0:
                row.append(VIE_INF)
            elif i == 7 or j == 7:
                row.append(0)
            else:
                row.append(value_of(i, j))
        m.append(row)
    return m


def mismatch_block(pair_index, values_of):
    block = [[0] * 5 for _ in range(5)]
    for i in range(1, 5):
        for j in range(1, 5):
            block[i][j] = values_of(pair_index, i - 1, j - 1)
    return block


def build_mismatch(values_of):
    mm = [[[VIE_INF] * 5 for _ in range(5)]]                 # pair 0 sentinel
    mm += [mismatch_block(p, values_of) for p in range(1, 7)]
    mm.append([[0] * 5 for _ in range(5)])                   # pair 7 = NN
    return mm


def build_dangle(dangle):
    dg = [[VIE_INF] * 5]                                     # N sentinel
    dg += [[0] + dangle[label] for label in CANON]
    dg.append([0] * 5)                                       # NN
    return dg


def build_int11(table):
    def value(p1, p2, a, b):
        if not (1 <= p1 <= 6 and 1 <= p2 <= 6 and 1 <= a <= 4 and 1 <= b <= 4):
            return VIE_INF
        return table[(CANON[p1 - 1], reverse_pair(CANON[p2 - 1]))][a - 1][b - 1]
    return [[[[value(p1, p2, a, b) for b in range(5)] for a in range(5)]
             for p2 in range(8)] for p1 in range(8)]


def build_int21(table):
    def value(p1, p2, a3, a4, a5):
        if not (1 <= p1 <= 6 and 1 <= p2 <= 6 and 1 <= a3 <= 4 and 1 <= a4 <= 4 and 1 <= a5 <= 4):
            return VIE_INF
        grid = table[(CANON[p1 - 1], reverse_pair(CANON[p2 - 1]), NT[a4 - 1])]
        return grid[a3 - 1][a5 - 1]
    return [[[[[value(p1, p2, a3, a4, a5) for a5 in range(5)] for a4 in range(5)]
              for a3 in range(5)] for p2 in range(8)] for p1 in range(8)]


def build_int22(table):
    def value(p1, p2, a3, a4, a5, a6):
        if not (1 <= p1 <= 6 and 1 <= p2 <= 6 and min(a3, a4, a5, a6) >= 1 and max(a3, a4, a5, a6) <= 4):
            return VIE_INF
        grid = table[(CANON[p1 - 1], reverse_pair(CANON[p2 - 1]))]
        return grid[4 * (a3 - 1) + (a6 - 1)][4 * (a4 - 1) + (a5 - 1)]
    return [[[[[[value(p1, p2, a3, a4, a5, a6) for a6 in range(5)] for a5 in range(5)]
               for a4 in range(5)] for a3 in range(5)] for p2 in range(8)] for p1 in range(8)]


def build_vienna(src):
    groups = read_number_groups(src / "rna.miscloop.dg")
    stack_d = parse_stack(src / "rna.stack.dg")
    internal, bulge, hairpin = parse_loops(src / "rna.loop.dg")
    internal[2] = internal[3] = 100  # ViennaRNA convention (small internal loops)
    dangle3, dangle5 = parse_dangles(src / "rna.dangle.dg")

    def parsed_value(tab):
        return lambda p, i, j: tab[CANON[p - 1]][i][j]

    def reversed_value(tab):
        return lambda p, i, j: tab[reverse_pair(CANON[p - 1])][j][i]

    return {
        "stack": pair_matrix(lambda i, j: stack_d[CANON[i - 1]]
                             [NT_INDEX[CANON[j - 1][1]]][NT_INDEX[CANON[j - 1][0]]]),
        "mm_hairpin": build_mismatch(parsed_value(parse_mismatch(src / "rna.tstackh.dg"))),
        "mm_internal": build_mismatch(parsed_value(parse_mismatch(src / "rna.tstacki.dg"))),
        "mm_internal_1n": build_mismatch(parsed_value(parse_mismatch(src / "rna.tstacki1n.dg"))),
        "mm_internal_23": build_mismatch(parsed_value(parse_mismatch(src / "rna.tstacki23.dg"))),
        "mm_multi": build_mismatch(reversed_value(parse_mismatch(src / "rna.tstackm.dg"))),
        "mm_exterior": build_mismatch(reversed_value(parse_mismatch(src / "rna.tstack.dg"))),
        "dangle5": build_dangle(dangle5),
        "dangle3": build_dangle(dangle3),
        "hairpin": hairpin, "bulge": bulge, "internal": internal,
        "int11": build_int11(parse_int11(src / "rna.int11.dg")),
        "int21": build_int21(parse_int21(src / "rna.int21.dg")),
        "int22": build_int22(parse_int22(src / "rna.int22.dg")),
        "tri": parse_specials(src / "rna.triloop.dg"),
        "tetra": parse_specials(src / "rna.tloop.dg"),
        "hexa": parse_specials(src / "rna.hexaloop.dg"),
        "ML_closing": scale(groups[3][0]),
        "ML_BASE": scale(groups[3][1]),
        "ML_intern": -90,  # ViennaRNA convention
        "NINIO_max": scale(groups[1][0]),
        "NINIO_m": scale(groups[2][0]),
        "TerminalAU": scale(groups[7][0]),
        "lxc": float(Decimal(groups[0][0]) * 100),
    }


# --- re-index ViennaRNA [8]/[5] arrays into t2004.h [7]/[6] layout -----------

def reindex(v):
    t = {}
    t["stack"] = [[v["stack"][a + 1][b + 1] for b in range(7)] for a in range(7)]
    for name in ("mm_hairpin", "mm_internal", "mm_internal_1n", "mm_internal_23",
                 "mm_multi", "mm_exterior"):
        t[name] = [[[v[name][a + 1][x][y] for y in range(5)] for x in range(5)] for a in range(7)]
    for name in ("dangle5", "dangle3"):
        t[name] = [[v[name][a + 1][x] for x in range(5)] for a in range(7)]
    t["int11"] = [[[[v["int11"][a + 1][b + 1][x][y] for y in range(5)] for x in range(5)]
                   for b in range(7)] for a in range(7)]
    t["int21"] = [[[[[v["int21"][a + 1][b + 1][x][y][z] for z in range(5)] for y in range(5)]
                    for x in range(5)] for b in range(7)] for a in range(7)]
    t["int22"] = [[[[[[v["int22"][a + 1][b + 1][w + 1][x + 1][y + 1][z + 1]
                       for z in range(4)] for y in range(4)] for x in range(4)]
                    for w in range(4)] for b in range(6)] for a in range(6)]
    for k in ("hairpin", "bulge", "internal", "tri", "tetra", "hexa", "ML_closing",
              "ML_intern", "ML_BASE", "NINIO_max", "NINIO_m", "TerminalAU", "lxc"):
        t[k] = v[k]
    return t


# --- diff against the in-tree tables (folding-reachable entries only) --------

def flat(x):
    if isinstance(x, list):
        for e in x:
            yield from flat(e)
    else:
        yield x


def reshape(flat_list, shape):
    if len(shape) == 1:
        return list(flat_list)
    step = len(flat_list) // shape[0]
    return [reshape(flat_list[i * step:(i + 1) * step], shape[1:]) for i in range(shape[0])]


CUR_SHAPES = {
    "stack": (7, 7), "mm_hairpin": (7, 5, 5), "mm_internal": (7, 5, 5),
    "mm_internal_1n": (7, 5, 5), "mm_internal_23": (7, 5, 5), "mm_multi": (7, 5, 5),
    "mm_exterior": (7, 5, 5), "dangle5": (7, 5), "dangle3": (7, 5),
    "int11": (7, 7, 5, 5), "int21": (7, 7, 5, 5, 5), "int22": (6, 6, 4, 4, 4, 4),
    "hairpin": (31,), "bulge": (31,), "internal": (31,),
}


def parse_c_array(text, name):
    m = re.search(re.escape(name) + r"\s*(?:\[[^\]]*\])+\s*=\s*", text)
    if not m:
        raise KeyError(name)
    i = text.index("{", m.end()); depth = 0; j = i
    while j < len(text):
        depth += (text[j] == "{") - (text[j] == "}")
        if depth == 0:
            break
        j += 1
    py = re.sub(r"/\*.*?\*/", "", text[i:j + 1]).replace("{", "[").replace("}", "]")
    return eval(re.sub(r",\s*\]", "]", py))


def scalar(text, name, cast=int):
    return cast(re.search(name + r"\s*=\s*([-0-9.]+)", text).group(1))


def load_current(path):
    text = path.read_text()
    d = {n: reshape(list(flat(parse_c_array(text, n))), CUR_SHAPES[n]) for n in CUR_SHAPES}
    for k, cast in (("ML_closing", int), ("ML_intern", int), ("ML_BASE", int),
                    ("NINIO_max", int), ("NINIO_m", int), ("TerminalAU", int), ("lxc", float)):
        d[k] = scalar(text, k, cast)
    return d


def canon(name, arr):
    """Folding-reachable entries: canonical pairs (0..5), real bases (1..4)."""
    P, B = range(6), range(1, 5)
    if name == "stack":
        return [arr[a][b] for a in P for b in P]
    if name.startswith("mm_"):
        return [arr[a][x][y] for a in P for x in B for y in B]
    if name.startswith("dangle"):
        return [arr[a][x] for a in P for x in B]
    if name == "int11":
        return [arr[a][b][x][y] for a in P for b in P for x in B for y in B]
    if name == "int21":
        return [arr[a][b][x][y][z] for a in P for b in P for x in B for y in B for z in B]
    if name in ("hairpin", "bulge", "internal"):
        return arr[3:31]
    return list(flat(arr))  # int22 already canonical


def diff_canonical(new, cur):
    reps = []
    for name in CUR_SHAPES:
        a, b = canon(name, new[name]), canon(name, cur[name])
        d = [(x, y) for x, y in zip(a, b) if x != y]
        if len(a) != len(b):
            reps.append(f"  {name}: SHAPE new={len(a)} cur={len(b)}")
        elif d:
            reps.append(f"  {name}: {len(d)}/{len(a)} diffs, e.g. new/cur {d[:4]}")
    for k in ("ML_closing", "ML_intern", "ML_BASE", "NINIO_max", "NINIO_m", "TerminalAU", "lxc"):
        if new[k] != cur[k]:
            reps.append(f"  {k}: new={new[k]} cur={cur[k]}")
    return reps


# --- emit t2004.h -----------------------------------------------------------

def cflat(arr, int16=True):
    # Nested brace-initializer matching the array shape (so the header needs no
    # -Wno-missing-braces). Dead sentinel entries (VIE_INF) are clamped to fit
    # the target type: 0 for the int16 tables, our INF for the (int) loop arrays.
    if isinstance(arr, list):
        return "{" + ",".join(cflat(e, int16) for e in arr) + "}"
    return str(0 if (int16 and arr >= VIE_INF) else (INF if arr >= VIE_INF else arr))


def emit_t2004(new, lxc):
    L = [
        "// Turner 2004 energy tables for src/energy.h.",
        "// Generated by util/gen_t2004.py directly from the NNDB turner_2004 .dg bundle.",
        "// Layout: pair index = ViennaRNA pair type - 1; mismatch/dangle base index",
        "// bidx(N)=0, A..U=1..4; int22 canonical-only. Sentinel (NN pair, N base)",
        "// entries are folding-unreachable and set to 0. lxc is kept at ViennaRNA's",
        "// precise 107.856 (the .dg file rounds it to 107.9) for RNAfold-exact",
        "// large-loop extrapolation.",
        "#pragma once",
        "#include <cstdint>",
        "namespace t2004 {",
        "static const int INF = 1000000;",
        f"static const int16_t stack[7][7] = {cflat(new['stack'])};",
    ]
    for nm in ("hairpin", "bulge", "internal"):
        L.append(f"static const int {nm}[31] = {cflat(new[nm], int16=False)};")
    for nm in ("mm_hairpin", "mm_internal", "mm_internal_1n", "mm_internal_23",
               "mm_multi", "mm_exterior"):
        L.append(f"static const int16_t {nm}[7][5][5] = {cflat(new[nm])};")
    for nm in ("dangle5", "dangle3"):
        L.append(f"static const int16_t {nm}[7][5] = {cflat(new[nm])};")
    L.append(f"static const int16_t int11[7][7][5][5] = {cflat(new['int11'])};")
    L.append(f"static const int16_t int21[7][7][5][5][5] = {cflat(new['int21'])};")
    L.append(f"static const int16_t int22[6][6][4][4][4][4] = {cflat(new['int22'])};")
    L.append(f"static const int ML_BASE={new['ML_BASE']}, ML_closing={new['ML_closing']}, ML_intern={new['ML_intern']};")
    L.append(f"static const int NINIO_m={new['NINIO_m']}, NINIO_max={new['NINIO_max']};")
    L.append(f"static const int TerminalAU={new['TerminalAU']};")
    L.append(f"static const double lxc={lxc};")
    for struct, name, data in (("SPtetra", "tetra", new["tetra"]),
                               ("SPtri", "tri", new["tri"]),
                               ("SPhexa", "hexa", new["hexa"])):
        L.append("struct %s{const char* seq;int bonus;};" % struct)
        L.append("static const %s %s[]={%s};" % (
            struct, name, ",".join('{"%s",%d}' % (s, v) for s, v in data)))
        L.append(f"static const int {name}_n={len(data)};")
    L.append("} // namespace t2004")
    return "\n".join(L) + "\n"


def _find_dg_dir(root):
    for c in (root, root / "turner_2004"):
        if (c / "rna.stack.dg").is_file():
            return c
    matches = list(root.rglob("rna.stack.dg"))
    if matches:
        return matches[0].parent
    raise SystemExit(f"No .dg files under {root}")


def resolve_source(raw):
    path = Path(raw)
    if path.is_file() and path.suffix == ".zip":
        tmp = Path(tempfile.mkdtemp(prefix="turner2004_"))
        with zipfile.ZipFile(path) as z:
            z.extractall(tmp)
        return _find_dg_dir(tmp), tmp
    return _find_dg_dir(path), None


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--turner-source", required=True, help="turner_2004 .dg bundle (zip or dir)")
    ap.add_argument("--current", default="src/t2004.h", help="in-tree tables to diff against")
    ap.add_argument("--out", help="write t2004.h here (omit to only report the diff)")
    ap.add_argument("--lxc", type=float, default=107.856,
                    help="lxc override (default 107.856 = ViennaRNA/RNAfold; .dg rounds to 107.9)")
    ap.add_argument("--check-only", action="store_true")
    args = ap.parse_args()

    src, tmp = resolve_source(args.turner_source)
    try:
        new = reindex(build_vienna(src))
        cur_path = Path(args.current)
        if cur_path.is_file():
            reps = diff_canonical(new, load_current(cur_path))
            print("Canonical values IDENTICAL to current t2004.h." if not reps
                  else "DIFFERENCES from current t2004.h (folding-reachable):\n" + "\n".join(reps))
        if args.out and not args.check_only:
            Path(args.out).write_text(emit_t2004(new, args.lxc))
            print(f"wrote {args.out} (lxc={args.lxc})")
    finally:
        if tmp:
            shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    main()
