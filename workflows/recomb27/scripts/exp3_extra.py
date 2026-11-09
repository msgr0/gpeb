"""exp3_extra.py CHROM MAPDIR READSDIR OUT    (rule exp3_extra; same arms and pairs as exp3_aggregate.py)
Per graph x rung -> OUT.extra.tsv: edits capped at 100 over mapped reads, share of mapped reads with more than 100 edits,
the length of the unmapped reads (rung FASTQ, read in the order of the .peread.gz lines), and bases: all, mapped, and over
the mapped reads clipped (clip) and not matched inside the alignment (core = mismatched + inserted + deleted).
Paired cover N vs the same others as exp3_aggregate.py -> OUT.paired_extra.tsv, on the reads both graphs map,
with the 95% bootstrap interval: capped edits, and the share over 100 edits; and their bases, clip and core on each graph."""
import sys, os, re, glob, gzip
import numpy as np

CAP = 100
chrom, mapdir, readsdir, out = sys.argv[1:5]
rng = np.random.default_rng(0)

def lengths(fq):
    with open(fq, "rb") as f:
        return np.array([len(l) - 1 for i, l in enumerate(f) if i % 4 == 1], dtype=np.int64)

def boot(d, B=1000):
    if len(d) == 0: return (float("nan"),) * 3
    bs = np.sort(np.concatenate([d[rng.integers(0, len(d), (50, len(d)))].mean(1) for _ in range(B // 50)]))
    return d.mean(), bs[int(.025 * B)], bs[int(.975 * B) - 1]

per, rows = {}, []
for p in sorted(glob.glob(f"{mapdir}/*x/*.peread.gz")):
    cov, gid = os.path.basename(os.path.dirname(p))[:-1], os.path.basename(p)[:-len(".peread.gz")]
    per[(cov, gid)] = np.loadtxt(gzip.open(p, "rt"), dtype=np.int64, ndmin=2)
L = {cov: lengths(f"{readsdir}/HG002.G{cov}x.fq") for cov in {c for c, _ in per}}

for (cov, gid), v in sorted(per.items(), key=lambda kv: (float(kv[0][0]), kv[0][1])):
    assert len(v) == len(L[cov]), (cov, gid, len(v), len(L[cov]))
    m = v[:, 3] == 1; e = v[m, 0]; ul = L[cov][~m]
    rows.append(f"{chrom}\t{cov}\t{gid}\t{len(v)}\t{m.sum()}\t{np.minimum(e, CAP).mean():.4f}\t{(e > CAP).mean():.6f}\t"
                f"{(~m).sum()}\t{np.median(ul) if len(ul) else float('nan'):.0f}\t{ul.mean() if len(ul) else float('nan'):.0f}\t{np.median(L[cov]):.0f}\t"
                f"{L[cov].sum()}\t{L[cov][m].sum()}\t{v[m, 2].sum()}\t{v[m, 1].sum()}")
with open(out + ".extra.tsv", "w") as f:
    f.write("chrom\tcov\tgid\treads\tmapped_n\tedits_cap100_mapped\tover100_mapped\tunmapped_n\tunmapped_len_median\tunmapped_len_mean\treads_len_median\t"
            "bp\tbp_mapped\tclip_bp\tcore_bp\n")
    f.write("\n".join(rows) + "\n")

with open(out + ".paired_extra.tsv", "w") as f:
    f.write("chrom\tcov\tn\tagainst\tboth_mapped\tdelta_capmapped\tlo_capmapped\thi_capmapped\tdelta_over100\tlo_over100\thi_over100\t"
            "both_bp\tclip_bp\tother_clip_bp\tcore_bp\tother_core_bp\n")
    for (cov, gid), g in sorted(per.items(), key=lambda kv: (float(kv[0][0]), kv[0][1])):
        mm = re.match(r"gpeb\.cover(\d+)\.c", gid)
        if not mm: continue
        n = int(mm[1])
        others = [o for o in per if o[0] == cov and (re.match(rf"vg\.s{n}\.ci\d\.c", o[1]) or (n == 32 and o[1].startswith("vg.d32.")) or
                                                     o[1].startswith(f"rand.n{n}.") or o[1] == "fullpanel")]
        for o in sorted(others):
            v = per[o]
            both = (g[:, 3] == 1) & (v[:, 3] == 1)
            a, b = g[both, 0], v[both, 0]
            d1 = boot(np.minimum(a, CAP) - np.minimum(b, CAP))
            d2 = boot((a > CAP).astype(np.float64) - (b > CAP))
            f.write(f"{chrom}\t{cov}\t{n}\t{o[1]}\t{both.sum()}\t" + "\t".join(f"{x:.4f}" for x in d1) + "\t" + "\t".join(f"{x:.6f}" for x in d2) +
                    f"\t{L[cov][both].sum()}\t{g[both, 2].sum()}\t{v[both, 2].sum()}\t{g[both, 1].sum()}\t{v[both, 1].sum()}\n")
print(f"{chrom}: {len(rows)} graph x rung rows", file=sys.stderr)
