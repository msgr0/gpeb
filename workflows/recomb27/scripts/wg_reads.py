#!/usr/bin/env python3
"""Whole-genome rungs from raw reads, with the hash of subset.py: u = crc32(name)/2^32, rung c = reads with u < c/depth.
  wg_reads.py sup UMAX OUT.fa.gz OUT.bases < FASTA       one uBAM (samtools fasta, one line per sequence): keep u < UMAX, count every base
  wg_reads.py rungs UMAX GENOME_LEN_FILE OUTPFX C1,C2,... SUP.fa.gz:BASES ...   depth = all bases / genome length;
                                                         OUTPFX.G<c>x.fa.gz per rung; fails if a rung needs u >= UMAX
Rungs are nested and the same read always lands in the same rungs, whichever uBAM it came from."""
import sys, zlib, gzip


def u(h):   # mates NAME/1 and NAME/2 (short reads) hash as NAME, so a pair always lands in the same rungs
    n = h[1:].split()[0]
    return zlib.crc32((n[:-2] if n[-2:] in ("/1", "/2") else n).encode()) / 2**32


def sup(umax, out, nout, inp=sys.stdin):
    bases = n = 0
    with gzip.open(out, "wt", compresslevel=1) as o:
        for h in inp:
            s = inp.readline()
            assert h[0] == ">" and s, "not one-line FASTA"
            bases += len(s) - 1; n += 1
            if u(h) < umax:
                o.write(h); o.write(s)
    with open(nout, "w") as f:
        print(f"{bases}\t{n}", file=f)


def rungs(glen, pfx, covs, sups):
    bases = sum(int(open(b).read().split()[0]) for _, b in sups)
    depth = bases / glen
    t = {c: float(c) / depth for c in covs}
    fh = {c: gzip.open(f"{pfx}.G{c}x.fa.gz", "wt", compresslevel=1) for c in covs}
    top = n = 0
    for fa, _ in sups:
        with gzip.open(fa, "rt") as f:
            for h in f:
                s = f.readline(); x = u(h); top = max(top, x); n += 1
                for c in covs:
                    if x < t[c]:
                        fh[c].write(h); fh[c].write(s)
    for o in fh.values():
        o.close()
    print(f"{bases} bases / {glen} = {depth:.2f}x; {n} reads in the first pass, u up to {top:.4f}; "
          + " ".join(f"{c}x:u<{t[c]:.4f}" for c in covs), file=sys.stderr)
    return depth, t, top


def selftest():
    import io, os, tempfile
    d = tempfile.mkdtemp()
    reads = [f">r{i}/ccs\n{'ACGT' * 25}\n" for i in range(2000)]          # 2000 reads x 100 bp = 200 kbp
    sup(0.5, f"{d}/a.fa.gz", f"{d}/a.n", io.StringIO("".join(reads[:1000])))
    sup(0.5, f"{d}/b.fa.gz", f"{d}/b.n", io.StringIO("".join(reads[1000:])))
    depth, t, top = rungs(100000, f"{d}/x", ["0.1", "0.5"], [(f"{d}/a.fa.gz", f"{d}/a.n"), (f"{d}/b.fa.gz", f"{d}/b.n")])
    assert depth == 2.0 and top < 0.5
    got = {c: [l for l in gzip.open(f"{d}/x.G{c}x.fa.gz", "rt") if l[0] == ">"] for c in ("0.1", "0.5")}
    want = {c: [r.split("\n")[0] + "\n" for r in reads if u(r) < float(c) / 2] for c in ("0.1", "0.5")}
    assert got == want and set(got["0.1"]) < set(got["0.5"]), "rungs differ from the hash, or are not nested"
    assert u(">p7/1 1:N\n") == u(">p7/2 2:N\n") == u(">p7\n") != u(">p8/1\n") and u(">m/ccs\n") != u(">m\n")   # mates together, HiFi untouched
    print("selftest ok")


if __name__ == "__main__":
    a = sys.argv[1:]
    if a[0] == "--selftest": selftest()
    elif a[0] == "sup": sup(float(a[1]), a[2], a[3])
    else:
        umax = float(a[1]); glen = int(open(a[2]).read()); pfx, covs = a[3], a[4].split(",")
        depth, t, top = rungs(glen, pfx, covs, [tuple(x.split(":")) for x in a[5:]])
        assert max(t.values()) <= umax, f"depth {depth:.1f}x too low: rung {max(covs, key=float)}x needs u < {max(t.values()):.4f} > {umax}"
