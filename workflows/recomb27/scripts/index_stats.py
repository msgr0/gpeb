#!/usr/bin/env python3
"""index_stats.py WORKDIR CHROM PANEL OUT.tsv [CFG GFA]
One row per WORKDIR/idx/<chrom>.<panel>.<cfg>.pex (only CFG if given, with the bytes of its GFA in gfa_bytes: Table 3): what its build printed (logs/index.<chrom>.<panel>.<cfg>.log), its bytes, and the
GNU time line (bench/index__ or max_colours__). cfg = mg<gap>.s<syncmer>[.max<N>]."""
import sys, os, re, glob
w, chrom, panel, out = sys.argv[1:5]
only, gfa = (sys.argv[5], sys.argv[6]) if len(sys.argv) > 6 else (None, None)
PAT = {   # the stderr lines of gpeb --index / --max-colours (gpeb.c: serial ~8385, syncmer ~1351, index ~2254, phases ~2260, graph, max-colours ~2693)
    "serial": (r"serial: gfa_read ([\d.]+)s\s+bg_build ([\d.]+)s\s+cut ([\d.]+)s", ["gfa_read_s", "bg_build_s", "cut_s"]),
    "syncmer": (r"syncmer: s=\d+\s+(\d+)/(\d+) window k-mers kept \(([\d.]+)\), (\d+)/(\d+) windows left with none",
                ["sync_kept", "sync_cand", "sync_frac", "windows_dark", "windows_all"]),
    "index": (r"index: k=\d+ span=\d+ s=\d+ colors=(\d+) windows=(\d+) index=(\d+) kmers \((\d+) ambiguous dropped\)",
              ["colors", "windows", "kmers", "ambiguous_dropped"]),
    "phases": (r"phases: insert ([\d.]+)s\s+segment-scan ([\d.]+)s\s+walk-stream ([\d.]+)s\s+finalize ([\d.]+)s",
               ["insert_s", "segscan_s", "walk_s", "finalize_s"]),
    "graph": (r"graph: segments=(\d+) links=(\d+) walked_edges=(\d+)", ["segments", "links", "walked_edges"]),
    "maxcol": (r"max-colours: (\d+) -> (\d+) records", ["src_records", "kmers"]),
}
COLS = ["chrom", "panel", "cfg", "min_gap", "syncmer", "max_colours", "gfa_bytes", "bytes"] + [c for _, cs in PAT.values() for c in cs if c != "kmers"] + \
       ["kmers", "wall_s", "maxrss_kb", "status"]
rows = []
for p in sorted(glob.glob(f"{w}/idx/{chrom}.{panel}.*.pex")):
    cfg = os.path.basename(p)[len(f"{chrom}.{panel}."):-4]
    m = re.match(r"mg(\d+)\.s(\d+)(?:\.max(\d+))?$", cfg)
    if not m or (only and cfg != only): continue
    r = dict(chrom=chrom, panel=panel, cfg=cfg, min_gap=m[1], syncmer=m[2], max_colours=m[3] or "", bytes=os.path.getsize(p),
             gfa_bytes=os.path.getsize(gfa) if gfa else "")
    lg = f"{w}/logs/index.{chrom}.{panel}.{cfg}.log"
    txt = open(lg, errors="replace").read() if os.path.exists(lg) else ""
    for pat, cs in PAT.values():
        x = re.search(pat, txt)
        if x: r.update(zip(cs, x.groups()))
    t = f"{w}/bench/{'max_colours' if m[3] else 'index'}__{chrom}__{panel}.{'max' + m[3] if m[3] else cfg}.time"
    if os.path.exists(t):
        ls = [l.strip() for l in open(t) if l.strip()]
        v = ls[-1].split("\t")
        r.update(wall_s=v[0], maxrss_kb=v[1], status=ls[0] if ls[0].startswith("Command") else "ok")
    rows.append(r)
with open(out, "w") as o:
    o.write("\t".join(COLS) + "\n")
    for r in rows: o.write("\t".join(str(r.get(c, "NA")) for c in COLS) + "\n")
print(f"{len(rows)} indexes", file=sys.stderr)
