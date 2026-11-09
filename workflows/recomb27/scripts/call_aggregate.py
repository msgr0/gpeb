#!/usr/bin/env python3
"""call_aggregate.py OUT.tsv EXPECT DIR [DIR ...]      one row per call
DIR = calls/<reads>/<chrom>/<panel>/<cfg>/<kind>.<sample>.G<cov>x, holding gpeb's log (and cover.sets.tsv for a cover call);
EXPECT = HG002:HG002,HG003:HG002,... (the donor whose sample a read set should call: the parents' reads call their child).
rank (gpeb --set-call mode --top 0): the top sample, the expected donor's share F of the reads with hits (donor_score) and its rank,
  the best sample other than the donor (runner) and the margin donor_score - runner_score; D, D_lo, hardest and pi from the
  set-call[conf] line: the margin of the call over its hardest competitor and the share of the reads the call explains (gpeb call_conf).
cover / coverF<f> (cover.sets.tsv): the share of the informative 200 kb windows whose picks hold DONOR#1 / DONOR#2 (95% block
  bootstrap, blocks of 5 windows), first-pick shares, the best other color, mean picks per window.
call_aggregate.py --selftest"""
import sys, os, re, random, statistics as st
random.seed(0)   # the block bootstrap of the cover rows

def cover(d, donor):
    W = {}
    for l in open(os.path.join(d, "cover.sets.tsv")):
        if l.startswith("window"): continue
        w, s, e, hits, unc, rk, col, gain, ua = l.rstrip("\n").split("\t")
        W.setdefault(int(w), [])
        if int(rk) > 0: W[int(w)].append((int(rk), col))
    inf = [w for w in sorted(W) if W[w]]
    if not inf: return None
    first = {w: min(W[w])[1] for w in inf}
    got = lambda c: [1 if any(p[1] == c for p in W[w]) else 0 for w in inf]
    cnt = {}
    for w in inf: cnt[first[w]] = cnt.get(first[w], 0) + 1
    oth = sorted(((v, c) for c, v in cnt.items() if not c.startswith(donor + "#")), reverse=True)
    return dict(h1=got(donor + "#1"), h2=got(donor + "#2"), oth=got(oth[0][1]) if oth else [0] * len(inf), oth_name=oth[0][1] if oth else "",
                f1=sum(first[w] == donor + "#1" for w in inf) / len(inf), f2=sum(first[w] == donor + "#2" for w in inf) / len(inf),
                picks=st.mean(len(W[w]) for w in inf), windows=len(inf))

def boot(s, B=2000, k=5):
    bl = [s[i:i + k] for i in range(0, len(s), k)]
    if not bl: return (float("nan"),) * 3
    bs = sorted(sum(x for b in (random.choice(bl) for _ in bl) for x in b) / max(1, sum(len(b) for b in bl)) for _ in range(B))
    return st.mean(s), bs[int(.025 * B)], bs[int(.975 * B)]

def rank(d, donor):
    txt = open(os.path.join(d, "log")).read()
    samp = {s: float(sc) for s, sc in re.findall(r"set-call\[(?:sample|mode)\]\s+(\S+)\s+([\d.]+)\s+\(", txt)}
    if not samp: return None
    ranked = sorted(samp.items(), key=lambda x: -x[1]); names = [s for s, _ in ranked]
    oth = [(s, v) for s, v in ranked if s != donor]
    ds = samp.get(donor, 0.0)
    c = re.search(r"set-call\[conf\]\s+D (\S+)\s+D_lo (\S+)\s+hardest (\S+) .*?\)\s+pi (\S+)", txt)
    rp = re.search(r"with-hits=(\d+)", txt)
    return dict(r_plus=int(rp[1]) if rp else "", top=names[0], top1=names[0] == donor, donor_score=ds,
                donor_rank=names.index(donor) + 1 if donor in names else "",
                runner=oth[0][0] if oth else "", runner_score=oth[0][1] if oth else float("nan"), margin=ds - oth[0][1] if oth else float("nan"),
                **(dict(D=float(c[1]), D_lo=float(c[2]), hardest=c[3], pi=float(c[4])) if c else {}))


HDR = ("chrom panel cfg reads sample donor cov kind windows picks h1 h1_lo h1_hi h2 h2_lo h2_hi h1_first h2_first best_other "
       "best_other_picked r_plus top top1 donor_score donor_rank runner runner_score margin D D_lo hardest pi").split()


def aggregate(out, expect, dirs):
    rows = []
    for d in dirs:
        d = d.rstrip("/"); parts = d.split("/")
        reads, chrom, panel, cfg = parts[-5:-1]
        kind, sample, cov = re.match(r"(.+)\.(HG\d+)\.G(.+)x$", parts[-1]).groups()     # kind may hold a dot: coverF0.0001
        donor = expect[sample]; r = dict(chrom=chrom, panel=panel, cfg=cfg, reads=reads, sample=sample, donor=donor, cov=cov, kind=kind)
        if kind.startswith("cover"):
            c = cover(d, donor)
            if not c: continue
            a, b = boot(c["h1"]), boot(c["h2"])
            r.update(windows=c["windows"], picks=c["picks"], h1=a[0], h1_lo=a[1], h1_hi=a[2], h2=b[0], h2_lo=b[1], h2_hi=b[2],
                     h1_first=c["f1"], h2_first=c["f2"], best_other=c["oth_name"], best_other_picked=st.mean(c["oth"]))
        else:
            k = rank(d, donor)
            if not k: continue
            r.update(k)
        rows.append(r)
    with open(out, "w") as o:
        o.write("\t".join(HDR) + "\n")
        for r in sorted(rows, key=lambda r: (r["chrom"], r["panel"], r["cfg"], r["reads"], r["sample"], float(r["cov"]), r["kind"])):
            o.write("\t".join(f"{r[h]:.4g}" if isinstance(r.get(h), float) else str(r.get(h, "")) for h in HDR) + "\n")
    print(f"{len(rows)} rows", file=sys.stderr)


def selftest():   # gpeb's own lines (chr21 HG002 5x short reads, 2026-11-05)
    import tempfile
    t = tempfile.mkdtemp(); d = f"{t}/calls/sr/chr21/full/mg60.s16/rank.HG002.G5x"; os.makedirs(d)
    with open(f"{d}/log", "w") as f:
        f.write("classify: reads=1557286 with-hits=933017 (59.9%) resolved-to-one=14210 (1.5% of hits) no-hit=624269\n"
                "set-call[mode]      HG002            0.959  (895203/933017 reads, 2 haplotypes)   <- call\n  private x\n"
                "set-call[mode]      HG01784          0.914  (853201/933017 reads, 2 haplotypes)\n"
                "set-call[conf]      D 0.7343  D_lo 0.7290  hardest HG01252 (a 54694, b 8379)  pi 0.5990  largest HG002 0.5990, GRCh38 0.0682  (EM 1292 iterations, 88427 read sets)\n")
    aggregate(f"{t}/out.tsv", {"HG002": "HG002"}, [d])
    r = dict(zip(*[l.rstrip("\n").split("\t") for l in open(f"{t}/out.tsv")]))
    assert (r["reads"], r["top"], r["donor_score"], r["runner"], r["D"], r["hardest"], r["pi"], r["r_plus"]) == \
        ("sr", "HG002", "0.959", "HG01784", "0.7343", "HG01252", "0.599", "933017"), r
    print("selftest ok")


if __name__ == "__main__":
    a = sys.argv[1:]
    if a[0] == "--selftest": selftest()
    else: aggregate(a[0], dict(x.split(":") for x in a[1].split(",")), a[2:])
