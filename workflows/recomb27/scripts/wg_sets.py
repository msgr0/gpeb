#!/usr/bin/env python3
"""wg_sets.py: the whole-genome calls (rules/wholegenome.smk), HiFi and short reads. gpeb's per-read table is too large to keep (a
short read ties ~350 colours, 7 kB a line), so each per-chromosome call is reduced as it streams out of gpeb:
  wg_sets.py reduce DIR < gpeb --classify stdout   DIR/reads.tsv.gz: read, hits, votes, set id, for the reads with a hit;
                                                   DIR/sets.tsv.gz: set id, the samples of the read's most-voted colours V(r), once each
  wg_sets.py conf OUT DIR...                       one conf.py row per reduced call (chrom = chrN, every read)
  wg_sets.py logs OUT DIR...                       the same rows from each call's gpeb log (its set-call[conf] line), no Python EM
  wg_sets.py merge MDIR [--min-hits T,...] DIR...  each read on the chromosome where it has the most hits (then votes, then the
                                                   first DIR), the genome-wide call: MDIR/conf.tsv, one conf.py row per threshold T
                                                   (default 1) on the reads with >= T hits there, chrom = wg for T = 1 and wg_h<T> else (a
                                                   short read with one or two hits can land on a chromosome by chance and hold one random
                                                   sample); MDIR/sets.tsv.gz (reads per sample set) and MDIR/chroms.tsv (reads per
                                                   chromosome), both for T = 1
The set of a read is the samples of V(r), taken from all_votes as conf.py does (top_colors and n_tied stop at ~4 kB)."""
import sys, os, gzip, re
from collections import Counter
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import conf


def reduce(d, inp=sys.stdin):
    ids, n, nhit = {}, 0, 0
    with gzip.open(f"{d}/reads.tsv.gz", "wt", compresslevel=1) as o:
        next(inp)
        for l in inp:
            n += 1
            x = l.split("\t", 10)
            if x[5] == "0": continue
            nhit += 1
            av = x[10].rstrip("\n")
            assert len(av) < 8100, (d, x[0], "all_votes may be cut short")
            # V(r) = the colours with the most votes, from all_votes: top_colors (and n_tied) stop at ~4 kB, 405 of 463 tied colours
            v = ",".join(sorted({c.split("#")[0] for c, _, k in (e.rpartition(":") for e in av.split(",")) if k == x[3]}))
            assert v, (d, x[0])
            o.write(f"{x[0]}\t{x[5]}\t{x[3]}\t{ids.setdefault(v, len(ids))}\n")
    with gzip.open(f"{d}/sets.tsv.gz", "wt", compresslevel=1) as o:
        for v, i in ids.items(): o.write(f"{i}\t{v}\n")
    print(f"reduce: {n} reads, {nhit} with hits, {len(ids)} distinct sample sets", file=sys.stderr)


def samples(d):   # set id -> its samples as a sorted tuple of interned names: ~1.3 kB for 150 samples, a frozenset of strings ~8 kB
    with gzip.open(f"{d}/sets.tsv.gz", "rt") as f:   # (conf.row only iterates and tests membership, so the call is the same)
        return {int(i): tuple(map(sys.intern, v.rstrip("\n").split(","))) for i, v in (l.split("\t") for l in f)}


def counts(d):
    s, c = samples(d), Counter()
    with gzip.open(f"{d}/reads.tsv.gz", "rt") as f:
        for l in f: c[int(l.rsplit("\t", 1)[1])] += 1
    out = Counter()
    for i, m in c.items(): out[s[i]] += m
    return out


def conf_one(d):
    return conf.row(d, counts(d))


CONF = re.compile(r"D (\S+)\s+D_lo (\S+)\s+hardest (\S+) \(a (\d+), b (\d+)\)\s+pi (\S+)\s+largest (\S+) (\S+), (\S+) (\S+)\s+\(EM (\d+) iterations, (\d+) read sets\)")


def log_row(d):   # the same row from gpeb's own log (its set-call[conf] line equals conf.py's; scores at 3 decimals)
    chrom, panel, cfg, leaf = d.rstrip("/").split("/")[-4:]
    sample, cov = re.match(r"rank\.(HG00\d)\.G([\d.]+)x(?:\.r0)?$", leaf).groups()
    t = open(f"{d}/log").read()
    rp = re.search(r"with-hits=(\d+)", t)[1]
    (top, st), (run, sr) = re.findall(r"set-call\[mode\]\s+(\S+)\s+(\S+)\s+\(", t)[:2]
    D, Dlo, h, a, b, pi, p1, v1, p2, v2, it, ns = CONF.search(t).groups()
    return "\t".join([chrom, panel, cfg, sample, cov, rp, ns, top, st, run, sr, D, Dlo, h, a, b, pi, p1, v1, p2, v2, it]) + "\n"


def merge(md, ds, ths=(1,)):
    best = {}
    for k, d in enumerate(ds):
        with gzip.open(f"{d}/reads.tsv.gz", "rt") as f:
            for l in f:
                r, h, v, i = l.split("\t")
                h, v = int(h), int(v); b = best.get(r)
                if b is None or (h, v) > b[:2]: best[r] = (h, v, k, int(i))   # one flat tuple a read: ~300 B with its key
    pick = {t: Counter((k, i) for h, _, k, i in best.values() if h >= t) for t in ths}; del best
    c, per = {t: Counter() for t in ths}, Counter()
    for k, d in enumerate(ds):     # one chromosome's sample sets in memory at a time
        s = samples(d)
        for t in ths:
            for (kk, i), m in pick[t].items():
                if kk == k:
                    c[t][s[i]] += m
                    if t == ths[0]: per[k] += m
    os.makedirs(md, exist_ok=True)
    with open(f"{md}/chroms.tsv", "w") as o:
        o.write("chrom\treads\n")
        for k, d in enumerate(ds): o.write(f"{d.rstrip('/').split('/')[-4]}\t{per[k]}\n")
    with gzip.open(f"{md}/sets.tsv.gz", "wt") as o:
        for ss, m in c[ths[0]].most_common(): o.write(f"{m}\t{','.join(sorted(ss))}\n")
    parts = md.rstrip("/").split("/")
    with open(f"{md}/conf.tsv", "w") as o:
        for t in ths:
            parts[-4] = "wg" if t == 1 else f"wg_h{t}"
            o.write(conf.row("/".join(parts), c[t]))


def selftest():
    import io, tempfile
    t = tempfile.mkdtemp()
    H = "#read\ttop_colors\tn_tied\tvotes\trunner_up\thits\tdisc_hits\tinter_n\tinter_colors\tn_masks\tall_votes\n"
    A = ["r1/1\tX#1,X#2\t2\t9\t2\t10\t10\t2\tX#1,X#2\t3\tX#1:9,X#2:9,Y#1:2", "r1/2\t-\t0\t0\t0\t0\t0\t0\t-\t0\t",
         "r2/1\tY#1\t1\t5\t1\t6\t6\t1\tY#1\t2\tY#1:5", "r3/1\tX#1,Y#2\t3\t4\t1\t5\t5\t1\tX#1\t2\tX#1:4,Y#2:4,Z#1:4,Y#1:1"]   # r3: top_colors cut
    B = ["r1/1\tZ#1\t1\t3\t1\t4\t4\t1\tZ#1\t2\tZ#1:3", "r1/2\tZ#1\t1\t3\t1\t4\t4\t1\tZ#1\t2\tZ#1:3",
         "r2/1\tZ#2\t1\t7\t1\t6\t6\t1\tZ#2\t2\tZ#2:7", "r3/1\t-\t0\t0\t0\t0\t0\t0\t-\t0\t"]
    ds = []
    for ch, ls in (("chrA", A), ("chrB", B)):
        ds.append(f"{t}/calls/{ch}/full/cfg/rank.HG002.G1x.r0"); os.makedirs(ds[-1])
        reduce(ds[-1], io.StringIO(H + "\n".join(ls) + "\n"))
    assert counts(ds[0]) == Counter({("X",): 1, ("Y",): 1, ("X", "Y", "Z"): 1})   # r3: its set from all_votes
    merge(f"{t}/calls/wg/full/cfg/rank.HG002.G1x.r0", ds, (1, 6))
    w, h = open(f"{t}/calls/wg/full/cfg/rank.HG002.G1x.r0/conf.tsv").read().splitlines()
    assert w.startswith("wg\tfull\tcfg\tHG002\t1\t4\t") and h.startswith("wg_h6\tfull\tcfg\tHG002\t1\t2\t"), (w, h)   # r1/1 (10 hits), r2 (6)
    m = f"{t}/calls/wg/full/cfg/rank.HG002.G1x.r0"
    assert open(f"{m}/chroms.tsv").read() == "chrom\treads\nchrA\t2\nchrB\t2\n"   # r1/1, r3 on A; r1/2, r2 (tied hits, more votes) on B
    got = {l.split("\t")[1].strip(): int(l.split("\t")[0]) for l in gzip.open(f"{m}/sets.tsv.gz", "rt")}
    assert got == {"X": 1, "X,Y,Z": 1, "Z": 2}, got
    assert open(f"{m}/conf.tsv").read().startswith("wg\tfull\tcfg\tHG002\t1\t4\t3\t")
    with open(f"{ds[0]}/log", "w") as f:   # gpeb's own lines (chr21 HG002 5x short reads, 2026-11-05)
        f.write("classify: reads=1557286 with-hits=933017 (59.9%) resolved-to-one=14210 (1.5% of hits) no-hit=624269\n"
                "set-call[mode]      HG002            0.959  (895203/933017 reads, 2 haplotypes)   <- call\n  private x\n"
                "set-call[mode]      HG01784          0.914  (853201/933017 reads, 2 haplotypes)\n"
                "set-call[conf]      D 0.7343  D_lo 0.7290  hardest HG01252 (a 54694, b 8379)  pi 0.5990  largest HG002 0.5990, GRCh38 0.0682  (EM 1292 iterations, 88427 read sets)\n")
    assert log_row(ds[0]) == "chrA\tfull\tcfg\tHG002\t1\t933017\t88427\tHG002\t0.959\tHG01784\t0.914\t0.7343\t0.7290\tHG01252\t54694\t8379\t0.5990\tHG002\t0.5990\tGRCh38\t0.0682\t1292\n", log_row(ds[0])
    print("selftest ok")


if __name__ == "__main__":
    a = sys.argv[1:]
    if a[0] == "--selftest": selftest()
    elif a[0] == "reduce": reduce(a[1])
    elif a[0] == "merge":
        if a[2] == "--min-hits": merge(a[1], a[4:], tuple(int(x) for x in a[3].split(",")))
        else: merge(a[1], a[2:])
    else:
        f = {"conf": conf_one, "logs": log_row}[a[0]]
        with Pool(int(os.environ.get("NPROC", 8))) as p, open(a[1], "w") as o:
            o.write("chrom\tpanel\tcfg\tsample\tcov\tr_plus\tn_sets\ttop\ts_top\trunner\ts_runner\tD\tD_lo\thardest\ta\tb\tpi_top\tpi1\tpi1_v\tpi2\tpi2_v\tem_iter\n")
            for r in p.imap_unordered(f, a[2:]): o.write(r); o.flush()
