#!/usr/bin/env python3
"""tex_numbers.py DATA     (DATA = <workdir>/results, written by the workflow)
Every number of the RECOMB 2027 paper's Results and supplement: the rows of its tables as LaTeX (tab1, ablation, wgindex, classes,
bases) and one `key: value` line for each number of the text, grouped by section. Copy the values into recomb27.tex and
recomb27-supp.tex; stdlib only."""
import sys, csv, glob, statistics as st, collections

D = sys.argv[1]
rd = lambda p: list(csv.DictReader(open(p), delimiter="\t"))
rows = lambda pat: [r for p in sorted(glob.glob(f"{D}/{pat}")) for r in rd(p)]
num = lambda x: float(x) if x not in ("", "NA", "nan", None) else float("nan")
ok = lambda x: x == x
GiB = lambda kb: num(kb) / 2**20
rng = lambda v, f=".2f": f"{min(v):{f}}-{max(v):{f}}" if v else "-"
DEF = "mg60.s16"
bench = rows("bench.tsv"); sizes = rows("sizes.tsv")
B = {(r["rule"], r["chrom"], r["key"]): r for r in bench if ok(num(r["wall_s"]))}
SZ = collections.defaultdict(int)
for r in sizes: SZ[(r["rule"], r["chrom"], r["key"])] += int(r["bytes"])
CH = sorted({p.split("/")[-1].split(".")[0] for p in glob.glob(f"{D}/exp2/chr*.tsv")}, key=lambda c: int(c[3:]))

# ---- Table 1: building the indexes of the panel (noHG002, ten chromosomes) -------------------------------------------------------
pex = [B[("index", c, f"noHG002.{DEF}")] for c in CH]
vgr = [B[(r, c, "panel")] for c in CH for r in ("vg_gbz", "vg_dist", "vg_ri", "vg_hapl")]
print(f"tab1 \\pindex{{}} & {sum(SZ[('index', c, f'noHG002.{DEF}')] for c in CH) / 1e9:.2f} & {max(GiB(r['maxrss_kb']) for r in pex):.0f} & "
      f"{sum(num(r['wall_s']) for r in pex) / 60:,.0f} \\\\")
print(f"tab1 \\texttt{{vg haplotypes}} & {sum(SZ[(r, c, 'panel')] for c in CH for r in ('vg_dist', 'vg_ri', 'vg_hapl')) / 1e9:.2f} & "
      f"{max(GiB(r['maxrss_kb']) for r in vgr):.0f} & {sum(num(r['wall_s']) for r in vgr) / 60:,.0f} \\\\   (disk without the GBZ)")
n9 = [c for c in CH if c != "chr22"]
print(f"index nine chromosomes other than 22: pebble index {sum(num(B[('index', c, f'noHG002.{DEF}')]['wall_s']) for c in n9) / 60:.0f} min, "
      f"vg haplotypes {sum(num(B[(r, c, 'panel')]['wall_s']) for c in n9 for r in ('vg_gbz', 'vg_dist', 'vg_ri', 'vg_hapl')) / 60:.0f} min")

# ---- the index: colors of the records (Figure 5), the whole genome (Table 3) ------------------------------------------------------
tot = one = le10 = ge460 = 0
for p in sorted(glob.glob(f"{D}/exp1/chr*.rare.tsv")):
    L = open(p).read().split("\n"); n = int(L[0].split("n=")[1]); tot += n
    for l in L[2:]:
        f = l.split("\t")
        if len(f) != 2 or not f[0].isdigit(): break       # the histogram ends where "private per colour:" starts
        k, v = int(f[0]), int(f[1]); one += v * (k == 1); le10 += v * (k <= 10); ge460 += v * (k >= 460)
print(f"colors: {tot:,} records on the full graphs; private {one / tot:.1%}, <= 10 colors {le10 / tot:.1%}, >= 460 colors {ge460 / tot:.1%}")
wi = rd(f"{D}/wg/index.tsv")
for r in wi:
    print(f"wgindex {r['chrom'][3:]} & {num(r['gfa_bytes']) / 1e9:.1f} & {num(r['kmers']) / 1e6:.1f} & {num(r['kmers']) / num(r['windows']):.2f} & "
          f"{100 * num(r['windows_dark']) / num(r['windows']):.3f} & {num(r['bytes']) / 1e9:.2f} & {num(r['wall_s']) / 60:.0f} & {GiB(r['maxrss_kb']):.0f} \\\\")
s = lambda k: sum(num(r[k]) for r in wi)
print(f"wgindex total & {s('gfa_bytes') / 1e9:.1f} & {s('kmers') / 1e6:.1f} & {s('kmers') / s('windows'):.2f} & {100 * s('windows_dark') / s('windows'):.3f} & "
      f"{s('bytes') / 1e9:.2f} & {s('wall_s') / 60:.0f} & {max(GiB(r['maxrss_kb']) for r in wi):.0f} \\\\   (index/GFA {s('bytes') / s('gfa_bytes'):.1%})")

# ---- Table 2: the index parameters sweep (full graph, rank calls and the parents' cover at 1x and 10x) ------------------------------
idx = {(r["chrom"], r["cfg"]): r for r in rows("exp1/chr*.index.tsv")}
for (c, cfg), r in idx.items():     # a default index linked from another run has its time in bench.tsv, not in its own log
    if not ok(num(r["wall_s"])) and ("index", c, f"full.{cfg}") in B:
        r["wall_s"], r["maxrss_kb"] = B[("index", c, f"full.{cfg}")]["wall_s"], B[("index", c, f"full.{cfg}")]["maxrss_kb"]
e1 = rows("exp1/chr*.calls.tsv")
SWEEP = [("default ($\\Delta=60$, $s=16$)", DEF)] + [(f"$\\Delta={g}$", f"mg{g}.s16") for g in (0, 30, 150, 500, 1000, 2000)] + \
        [("$s=0$ (all junction $k$-mers)", "mg60.s0")] + [(f"$s={x}$", f"mg60.s{x}") for x in (12, 20, 24)] + \
        [("at most $10$ colors", f"{DEF}.max10"), ("private ($1$ color)", f"{DEF}.max1")]
for lab, cfg in SWEEP:
    ix = [idx[(c, cfg)] for c in CH if (c, cfg) in idx]
    rk = [r for r in e1 if r["cfg"] == cfg and r["kind"] == "rank"]
    own = [r for r in rk if r["sample"] in ("HG002", "HG005")]; par = [r for r in rk if r["sample"] in ("HG003", "HG004")]
    marg = min(num(r["margin"]) for r in own if r["cov"] == "1")
    cv = [r for r in e1 if r["cfg"] == cfg and r["kind"] == "cover" and r["cov"] == "10" and r["sample"] in ("HG003", "HG004")]
    tr = st.mean(num(r["h1" if r["sample"] == "HG003" else "h2"]) for r in cv); un = st.mean(num(r["h2" if r["sample"] == "HG003" else "h1"]) for r in cv)
    w = sum(num(r["wall_s"]) for r in ix) / 60
    print(f"ablation {lab} & {sum(num(r['bytes']) for r in ix) / 1e6:,.0f} & {'$<1$' if w < 1 else f'{w:.0f}'} & {max(GiB(r['maxrss_kb']) for r in ix):.0f} & "
          f"{sum(r['top1'] == 'True' for r in own)}/{len(own)} & {sum(r['top'] == 'HG002' for r in par)}/{len(par)} & {marg:.2f} & {100 * tr:.0f} / {100 * un:.1f} \\\\")
for f in ("0", "0.0001", "0.005"):
    cv = [r for r in e1 if r["kind"] == ("cover" if f == "0.0001" else f"coverF{f}") and r["cfg"] == DEF and r["cov"] == "10"]
    par = [r for r in cv if r["sample"] in ("HG003", "HG004")]
    print(f"cover threshold tau={f}: colors per window {st.mean(num(r['picks']) for r in cv):.1f}, transmitted {st.mean(num(r['h1' if r['sample'] == 'HG003' else 'h2']) for r in par):.1%}")

# ---- calling the donor (Figure 2): ten chromosomes, HiFi and short reads --------------------------------------------------------
e2 = [r for r in rows("exp2/chr*.tsv") if r["kind"] == "rank"]
for reads in ("hifi", "sr"):
    x = [r for r in e2 if r["reads"] == reads]
    right = [r for r in x if r["panel"] == "full" or r["sample"] == "HG005"]; absent = [r for r in x if r not in right]
    acc = lambda v: sum(num(r["D"]) >= 0.25 for r in v)
    print(f"calls {reads}: {len(x)} calls; expected sample first {sum(r['top1'] == 'True' for r in right)}/{len(right)}, accepted {acc(right)}/{len(right)}, "
          f"rejected {len(absent) - acc(absent)}/{len(absent)}")
    for lab, g in (("own reads", [r for r in right if r["sample"] in ("HG002", "HG005")]), ("parents", [r for r in right if r["sample"] in ("HG003", "HG004")]), ("absent", absent)):
        print(f"   {lab}: F {rng([num(r['donor_score'] if r['top1'] == 'True' else r['runner_score']) for r in g])} "
              f"D {rng([num(r['D']) for r in g])} pi {rng([num(r['pi']) for r in g])}")
    for cov in ("0.5", "30"):
        t = [num(r["wall_s"]) for (ru, c, k), r in B.items() if ru == "rank_call" and k.startswith(f"{reads}.") and f".{DEF}." in k and k.endswith(f".{cov}")]
        print(f"   median call at {cov}x: {st.median(t):.0f} s (n={len(t)}), peak {max(GiB(r['maxrss_kb']) for (ru, c, k), r in B.items() if ru == 'rank_call' and k.startswith(f'{reads}.') and f'.{DEF}.' in k):.1f} GiB" if t else f"   {cov}x: -")

# ---- the whole genome (Figure 3) ----------------------------------------------------------------------------------------------
for reads in ("hifi", "sr"):
    w = rd(f"{D}/wg/{reads}.tsv")
    for panel in sorted({r["panel"] for r in w}):
        for ch in ("wg", "wg_h5"):
            g = [r for r in w if r["panel"] == panel and r["chrom"] == ch]
            if not g: continue
            by = collections.defaultdict(list)
            for r in g: by[r["sample"]].append(r)
            print(f"wg {reads} {panel} {ch}: " + "; ".join(f"{s} call {','.join(sorted({r['top'] for r in v}))} D {rng([num(r['D']) for r in v])} "
                                                         f"pi {rng([num(r['pi_top']) for r in v])}" for s, v in sorted(by.items())))
        alone = [r for r in w if r["panel"] == panel and r["chrom"].startswith("chr")]
        print(f"wg {reads} {panel} one chromosome alone: accepted {sum(num(r['D']) >= 0.25 for r in alone)} of {len(alone)}")

# ---- personalized graphs (Figure 4): the cover against each other graph, read by read ------------------------------------------
pair = rows("exp3/chr*.paired.tsv")
NS = sorted({int(r["n"]) for r in pair})
for lab, pre in (("vg haplotypes", "vg.s{n}."), ("random", "rand.n{n}."), ("diploid", "vg.d32.")):
    for n in NS:
        v = [r for r in pair if r["n"] == str(n) and r["against"].startswith(pre.format(n=n))]
        if not v: continue
        d = [num(r["delta_cap100"]) for r in v]; worse = sorted((r["chrom"], r["cov"]) for r in v if num(r["delta_cap100"]) > 0)
        print(f"exp3 cover{n} vs {lab}: better {sum(x < 0 for x in d)} of {len(d)}, delta {min(d):+.2f}..{max(d):+.2f}; worse at {worse}")
fp = [num(r["delta_cap100"]) for r in pair if r["against"] == "fullpanel"]
print(f"exp3 cover - noHG002 graph: {min(fp):+.2f}..{max(fp):+.2f} edits per read")
arms = rows("exp3/chr*[0-9].tsv")
for m in ("gpeb_cover", "vg_sample", "vg_diploid", "full_panel"):
    print(f"exp3 MAPQ 60 {m}: {st.mean(num(r['mapq60']) for r in arms if r['method'] == m):.1%} of the mapped reads")
CLS, LAB = ["perfect", "small", "clipped", "sv", "unmapped"], {"gpeb_cover": "cover", "vg_sample": "\\texttt{vg haplotypes}",
           "vg_diploid": "\\texttt{vg haplotypes} diploid", "random": "random", "full_panel": "noHG002 graph"}
for m, lab in LAB.items():          # Table: read classes, mean over the chromosomes and coverages (each weighs the same)
    for n in sorted({int(r["n"] or 0) for r in arms if r["method"] == m}):
        x = [r for r in arms if r["method"] == m and int(r["n"] or 0) == n]
        v = [100 * st.mean(num(r[c]) for r in x) for c in CLS + ["mapq60"]] + [st.mean(num(r["edits_cap100"]) for r in x)]
        print(f"classes {lab} & {'--' if m == 'full_panel' else 2 if m == 'vg_diploid' else n} & " + " & ".join(f"{a:.1f}" for a in v) + " \\\\")
ex, px = rows("exp3/chr*.extra.tsv"), rows("exp3/chr*.paired_extra.tsv")
bp = lambda gid, k: sum(int(r[k]) for r in ex if r["gid"] == gid or r["gid"].startswith(gid + "."))
print(f"bases: {bp('fullpanel', 'bp') / 1e9:.1f} Gbp on each graph")
for n, other, lab in [(n, f"vg.s{n}", "\\texttt{vg haplotypes}") for n in NS] + [(32, "vg.d32", "\\texttt{vg haplotypes} diploid")]:
    pr = [r for r in px if r["n"] == str(n) and r["against"].startswith(other + ".")]
    both = sum(int(r["both_bp"]) for r in pr)
    for i, (name, gid, ck, kk) in enumerate([("cover", f"gpeb.cover{n}", "clip_bp", "core_bp"), (lab, other, "other_clip_bp", "other_core_bp")]):
        a, un = bp(gid, "bp"), bp(gid, "bp") - bp(gid, "bp_mapped"); clip, core = sum(int(r[ck]) for r in pr), sum(int(r[kk]) for r in pr)
        print(f"bases {n if i == 0 else ''} & {name} & {un / 1e6:.1f} & {100 * un / a:.2f} & {both / 1e9 if i == 0 else '':{'.2f' if i == 0 else ''}} & "
              f"{clip / 1e6:.1f} & {100 * clip / both:.2f} & {core / 1e6:.1f} & {100 * core / both:.2f} \\\\")

# ---- personalized graphs: sizes and cost --------------------------------------------------------------------------------------
gbz = {(r["chrom"], r["key"]): int(r["bytes"]) for r in sizes if r["rule"] in ("gfa2gbz", "vg_hap")}
for n in NS:
    rat = [b / gbz[(c, f"s{n}.ci1.{k.split('.')[-1]}")] for (c, k), b in gbz.items() if k.startswith(f"gpeb.cover{n}.c") and (c, f"s{n}.ci1.{k.split('.')[-1]}") in gbz]
    mi = lambda pre: st.median(sum(SZ[(r, c, f"{pre}.c{v}")] for c in n9 for r in ("g_dist", "g_minzip")) for v in ("0.5", "1", "2", "5", "10", "30")) / 1e9
    print(f"graphs N={n}: cover GBZ / vg GBZ median {st.median(rat) if rat else float('nan'):.2f}; distance+minimizer+zipcode, nine chromosomes: "
          f"cover {mi(f'gpeb.cover{n}'):.1f} GB, vg {mi(f'vg.s{n}.ci1'):.1f} GB")
def steps(m, n, c):
    g = f"vg.s{n}.ci1.c{c}" if m == "vg" else f"gpeb.cover{n}.c{c}"
    return {"panel index": [(r, "panel") for r in ("vg_gbz", "vg_dist", "vg_ri", "vg_hapl")] if m == "vg" else [("index", f"noHG002.{DEF}")],
            "graph": [("kmc", f"c{c}.ci1"), ("vg_hap", f"s{n}.ci1.c{c}")] if m == "vg" else [("gpeb_cover", f"cover{n}.c{c}"), ("gfa2gbz", g)],
            "distance": [("g_dist", g)], "minimizer+zipcode": [("g_minzip", g)], "mapping": [("giraffe", f"{g}.on{c}")]}
for grp, chs in (("nine chromosomes", n9), ("chr22", ["chr22"])):
    for m in ("vg", "cover"):
        for n, c in ((4, "0.5"), (32, "30")):
            row = {s: (sum(num(B[(r, ch, k)]["wall_s"]) for ch in chs for r, k in L) / 60, max(GiB(B[(r, ch, k)]["maxrss_kb"]) for ch in chs for r, k in L))
                   for s, L in steps(m, n, c).items()}
            print(f"cost {grp} {m} N={n} {c}x: " + ", ".join(f"{s} {t:.0f} min ({p:.0f} GiB)" for s, (t, p) in row.items()) + f"; total {sum(t for t, _ in row.values()):.0f} min")
cpu = lambda k, ch: num(B[("giraffe", ch, k)]["user_s"]) + num(B[("giraffe", ch, k)]["sys_s"])
print("cost mapping CPU cover / vg haplotypes, N=32 at 30x:", {ch: round(cpu("gpeb.cover32.c30.on30", ch) / cpu("vg.s32.ci1.c30.on30", ch), 1) for ch in CH})
