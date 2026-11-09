#!/usr/bin/env python3
"""bench.py WORKDIR BENCH.tsv SIZES.tsv
WORKDIR/bench/<rule>__<chrom>__<key>.time (GNU time -f '%e\\t%M\\t%U\\t%S'; a failed or killed command adds a 'Command ...' line first)
and WORKDIR/sizes/<rule>__<chrom>__<key>.tsv ('path\\tbytes' per file) -> one table each, replaced atomically (targets may finish together)."""
import sys, os, glob
w, bo, so = sys.argv[1:4]
with open(bo + ".tmp", "w") as o:
    o.write("rule\tchrom\tkey\twall_s\tmaxrss_kb\tuser_s\tsys_s\tstatus\n")
    for p in sorted(glob.glob(f"{w}/bench/*.time")):
        ls = [l.strip() for l in open(p) if l.strip()]
        if not ls: continue
        st = ls[0] if ls[0].startswith("Command") else "ok"
        v = ls[-1].split("\t")
        o.write("\t".join(os.path.basename(p)[:-5].split("__", 2) + (v if len(v) == 4 else ["NA"] * 4) + [st]) + "\n")
with open(so + ".tmp", "w") as o:
    o.write("rule\tchrom\tkey\tfile\tbytes\n")
    for p in sorted(glob.glob(f"{w}/sizes/*.tsv")):
        k = os.path.basename(p)[:-4].split("__", 2)
        for l in open(p):
            f, b = l.rstrip("\n").rsplit("\t", 1)
            o.write("\t".join(k + [os.path.basename(f), b]) + "\n")
os.replace(bo + ".tmp", bo); os.replace(so + ".tmp", so)
