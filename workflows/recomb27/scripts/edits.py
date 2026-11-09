#!/usr/bin/env python3
"""edits.py READS_N GAF OUT_PEREAD FASTQ     (recomb27: v2's edits.py plus the largest indel, MAPQ and the read class)
Edits of a read = mismatches + inserted bases + deleted bases (GAF cs tag, case-insensitive) + query bases outside qs..qe.
An insertion at either end is a soft clip (giraffe writes clips as insertions): clip = leading + trailing insertion + outside bases.
maxindel = the longest insertion or deletion INSIDE the alignment (end clips excluded). Per read the record with fewest edits is kept.
Class: unmapped | perfect (0 edits) | sv (an inner indel >= 50 bp) | clipped (clip >= 50 bp, no sv) | small (any other edit).
OUT_PEREAD: one line per FASTQ read, in order: 'edits core clip mapped maxindel mapq' (unmapped = 100 100 0 0 0 0).
stdout: one key=value summary line (shares over all reads; MAPQ bins over mapped reads)."""
import sys, re
SV, CLIP, CAP = 50, 50, 100
TOK = re.compile(r":\d+|\*[A-Za-z][A-Za-z]|\+[A-Za-z]+|-[A-Za-z]+")

def record(t):
    L, qs, qe = int(t[1]), int(t[2]), int(t[3]); cs = next((x[5:] for x in t[12:] if x.startswith("cs:Z:")), None)
    if cs is None: return None
    m = ins = dele = lead = trail = big = 0; toks = TOK.findall(cs)
    for i, x in enumerate(toks):
        c, n = x[0], len(x) - 1
        if c == "*": m += 1
        elif c == "+":
            ins += n
            if i == 0: lead = n
            if i == len(toks) - 1: trail = n
            if 0 < i < len(toks) - 1: big = max(big, n)
        elif c == "-": dele += n; big = max(big, n)
    out = qs + (L - qe); clip = lead + trail + out
    tot = m + ins + dele + out
    return tot, tot - clip, clip, big, int(t[11])

def cls(r):
    if r is None: return "unmapped"
    tot, _, clip, big, _ = r
    return "perfect" if tot == 0 else "sv" if big >= SV else "clipped" if clip >= CLIP else "small"

def main(n_reads, path, out_per, fq):
    best = {}
    for l in open(path):
        if l[0] == "@": continue
        t = l.rstrip("\n").split("\t")
        if len(t) < 12 or t[5] == "*": continue
        r = record(t)
        if r is not None and (t[0] not in best or r[0] < best[t[0]][0]): best[t[0]] = r
    with open(fq) as g, open(out_per, "w") as o:
        for i, l in enumerate(g):
            if i % 4 == 0:
                r = best.get(l[1:].split()[0])
                o.write(f"{r[0]} {r[1]} {r[2]} 1 {r[3]} {r[4]}\n" if r else "100 100 0 0 0 0\n")
    v = list(best.values()); n = max(n_reads, 1); nm = max(len(v), 1)
    s = dict(reads=n_reads, mapped=len(v) / n)
    for c in ("perfect", "small", "clipped", "sv"): s[c] = sum(cls(r) == c for r in v) / n
    s["unmapped"] = 1 - len(v) / n
    q = [r[4] for r in v]
    s.update(mapq0=sum(x == 0 for x in q) / nm, mapq1_29=sum(0 < x < 30 for x in q) / nm, mapq30_59=sum(29 < x < 60 for x in q) / nm,
             mapq60=sum(x >= 60 for x in q) / nm, mapq_mean=sum(q) / nm,
             edits_cap100=(sum(min(r[0], CAP) for r in v) + CAP * (n_reads - len(v))) / n, edits_mapped=sum(r[0] for r in v) / nm)
    return s

def selftest():
    mk = lambda L, qs, qe, cs, mq=60: ["r", str(L), str(qs), str(qe), "+", ">1>2", "9", "0", "9", "0", "0", str(mq), f"cs:Z:{cs}"]
    assert record(mk(100, 0, 100, ":100")) == (0, 0, 0, 0, 60) and cls(record(mk(100, 0, 100, ":100"))) == "perfect"
    r = record(mk(100, 0, 100, ":9*CG:9*CG:32*TC:12+AATTCCATTC:27+G:46*AC:1")); assert r[0] == 15 and r[2] == 0 and r[3] == 10 and cls(r) == "small", r
    r = record(mk(100, 0, 100, "+GGATC:90-ACG:2")); assert r[0] == 8 and r[2] == 5 and r[3] == 3, r          # leading insertion is clip, not an indel
    r = record(mk(200, 0, 200, ":10-" + "A" * 60 + ":100")); assert r[3] == 60 and cls(r) == "sv", r           # 60 bp deletion inside
    r = record(mk(200, 0, 200, ":100+" + "C" * 70)); assert r[3] == 0 and r[2] == 70 and cls(r) == "clipped", r  # trailing 70 bp insertion = clip
    r = record(mk(100, 5, 90, ":85")); assert r[0] == 15 and r[2] == 15 and cls(r) == "small", r
    assert cls(None) == "unmapped"; print("edits.py selftest ok")

if __name__ == "__main__":
    if sys.argv[1] == "--selftest": selftest(); sys.exit(0)
    s = main(int(sys.argv[1]), sys.argv[2], sys.argv[3], sys.argv[4])
    print(" ".join(f"{k}={v:.6g}" if isinstance(v, float) else f"{k}={v}" for k, v in s.items()))
