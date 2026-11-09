#!/usr/bin/env python3
"""subset.py IN CHROM_LEN OUTDIR SAMPLE COV [COV ...]     the nested coverage ladder of one chromosome's reads
IN = FASTQ (HiFi) or one-line FASTA (short reads, samtools fasta -n). Depth = all bases / CHROM_LEN; rung c keeps the reads with
u = crc32(name) / 2^32 < c / depth (one draw: rungs are nested, and a rung above the depth is every read). Mates share a name, so a
pair lands in the same rungs. Writes OUTDIR/SAMPLE.G<c>x.fq for FASTQ, OUTDIR/SAMPLE.G<c>x.fa.gz for FASTA.
subset.py --selftest"""
import sys, zlib, os, gzip


def ladder(inp, clen, out, sample, covs):
    fa = open(inp).read(1) == ">"; L = 2 if fa else 4      # lines per record
    with open(inp) as f:
        bases = sum(len(l) - 1 for i, l in enumerate(f) if i % L == 1)
    depth = bases / clen
    os.makedirs(out, exist_ok=True)
    fh = {c: gzip.open(f"{out}/{sample}.G{c}x.fa.gz", "wt", compresslevel=1) if fa else open(f"{out}/{sample}.G{c}x.fq", "w") for c in covs}
    with open(inp) as f:
        while h := f.readline():
            rec = h + "".join(f.readline() for _ in range(L - 1))
            u = zlib.crc32(h[1:].split()[0].encode()) / 2**32
            for c in covs:
                if u < float(c) / depth: fh[c].write(rec)
    for x in fh.values(): x.close()
    print(f"{sample}: {bases} bases / {clen} = {depth:.2f}x", file=sys.stderr)
    return depth


def selftest():
    import tempfile
    d = tempfile.mkdtemp()
    with open(f"{d}/in.fq", "w") as f:
        for i in range(1000): f.write(f"@r{i} x\n{'ACGT' * 25}\n+\n{'I' * 100}\n")     # 100 kbp
    assert ladder(f"{d}/in.fq", 50000, d, "S", ["0.5", "1", "9"]) == 2.0
    got = {c: [l.split()[0][1:] for i, l in enumerate(open(f"{d}/S.G{c}x.fq")) if i % 4 == 0] for c in ("0.5", "1", "9")}
    want = lambda c: [f"r{i}" for i in range(1000) if zlib.crc32(f"r{i}".encode()) / 2**32 < float(c) / 2]
    assert got["0.5"] == want("0.5") and set(got["0.5"]) < set(got["1"]) and len(got["9"]) == 1000, "rungs differ from the hash"
    print("selftest ok")


if __name__ == "__main__":
    a = sys.argv[1:]
    if a[0] == "--selftest": selftest()
    else: ladder(a[0], int(a[1]), a[2], a[3], a[4:])
