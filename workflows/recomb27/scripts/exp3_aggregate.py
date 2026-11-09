#!/usr/bin/env python3
"""exp3_aggregate.py OUT_PREFIX CHROM MAP_DIR      (MAP_DIR = exp3/<chrom>/map, holding <cov>x/<gid>.stats and .peread.gz)
OUT_PREFIX.tsv: one row per (graph, mapped rung) with the edits.py summary (read classes, MAPQ bins, edits).
OUT_PREFIX.paired.tsv: gpeb cover N against each other graph mapped with the SAME reads (vg sample N, vg diploid 32 for N = 2 and 32,
random N, full panel): mean difference in edits per read (cap 100, unmapped = 100; and uncapped over reads mapped by both), 95% bootstrap
over reads. Negative = gpeb has fewer edits."""
import sys, os, re, gzip, glob
import numpy as np
out, chrom, root = sys.argv[1], sys.argv[2], sys.argv[3]
rng = np.random.default_rng(0)

def arm(gid):
    m = re.match(r"vg\.(d32|s(\d+))\.ci(\d)\.c(.+)$", gid)
    if m: return ("vg_diploid" if m[1] == "d32" else "vg_sample", int(m[2] or 32), m[3], "", m[4])
    m = re.match(r"gpeb\.cover(\d+)\.c(.+)$", gid)
    if m: return ("gpeb_cover", int(m[1]), "", "", m[2])
    m = re.match(r"rand\.n(\d+)\.s(\d+)$", gid)
    if m: return ("random", int(m[1]), "", m[2], "")
    return ("full_panel", 0, "", "", "")

rows, per = [], {}
for p in sorted(glob.glob(f"{root}/*x/*.stats")):
    cov, gid = os.path.basename(os.path.dirname(p))[:-1], os.path.basename(p)[:-len(".stats")]
    s = dict(kv.split("=") for kv in open(p).read().split())
    method, n, ci, seed, gcov = arm(gid)
    rows.append(dict(chrom=chrom, cov=cov, gid=gid, method=method, n=n, ci=ci, seed=seed, graph_cov=gcov, **s))
    per[(cov, gid)] = np.loadtxt(gzip.open(p[:-len(".stats")] + ".peread.gz", "rt"), dtype=np.int64, ndmin=2)
if not rows: sys.exit(f"no .stats under {root}")
keys = list(rows[0])
with open(out + ".tsv", "w") as f:
    f.write("\t".join(keys) + "\n")
    for r in sorted(rows, key=lambda r: (float(r["cov"]), r["method"], r["n"], r["gid"])):
        f.write("\t".join(str(r.get(k, "")) for k in keys) + "\n")

def boot(d, B=1000):
    if len(d) == 0: return (float("nan"),) * 3
    bs = np.sort(np.concatenate([d[rng.integers(0, len(d), (50, len(d)))].mean(1) for _ in range(B // 50)]))   # batches: 30x chr1 is ~0.5M reads
    return d.mean(), bs[int(.025 * B)], bs[int(.975 * B) - 1]

with open(out + ".paired.tsv", "w") as f:
    f.write("chrom\tcov\tn\tagainst\tgpeb_cap100\tother_cap100\tdelta_cap100\tlo\thi\tboth_mapped\tdelta_nocap\tlo_nocap\thi_nocap\n")
    for (cov, gid), g in sorted(per.items(), key=lambda kv: (float(kv[0][0]), kv[0][1])):
        m = re.match(r"gpeb\.cover(\d+)\.c", gid)
        if not m: continue
        n = int(m[1])
        others = [o for o in per if o[0] == cov and (re.match(rf"vg\.s{n}\.ci\d\.c", o[1]) or (n in (2, 32) and o[1].startswith("vg.d32.")) or
                                                     o[1].startswith(f"rand.n{n}.") or o[1] == "fullpanel")]
        for o in sorted(others):
            v = per[o]
            if len(v) != len(g): continue
            a, b = np.minimum(g[:, 0], 100), np.minimum(v[:, 0], 100)
            d, lo, hi = boot(a - b)
            both = (g[:, 3] == 1) & (v[:, 3] == 1)
            dn, lon, hin = boot(g[both, 0] - v[both, 0])
            f.write(f"{chrom}\t{cov}\t{n}\t{o[1]}\t{a.mean():.4f}\t{b.mean():.4f}\t{d:.4f}\t{lo:.4f}\t{hi:.4f}\t{both.sum()}\t{dn:.4f}\t{lon:.4f}\t{hin:.4f}\n")
print(f"{len(rows)} graph x rung rows", file=sys.stderr)
