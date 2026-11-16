#!/usr/bin/env python3
"""Topology checks behind the counterexamples to the k-mer lemma of the pebble index, on a GFA 1.1 graph
(S, L and W lines). Nothing of the k-mer index is computed: the checks look only at the graph, its labels and
its walks. Pebbles are placed as the linear-time test multilin does: a walk pebbles the edge by which it enters
an oriented vertex of in-degree >= 2.

  A  two consecutive pebbles of a walk whose junctions are at most k-2 bases apart, so their windows share a k-mer.
  B  a pebbled edge whose window (k-1 bases on each side of the junction) differs between two walks;
     computed on the walks of the first --sample haplotypes.
  C  non-determinism: an oriented vertex with two out-neighbours whose labels start with the same base, so a
     k-mer can be read from the same locus along two branches (a record need not identify its pebbled edge).
  D  a walk that enters a merge by the same edge twice (a cycle), so the edge is pebbled twice; also whether
     the two windows are identical.

usage: python3 utils/topo_check.py graph.gfa [-k 31] [-j 8] [--sample 12]     (needs numpy)
test:  python3 utils/test_topo_check.py
"""
import argparse, re, sys
from multiprocessing import Pool
import numpy as np

COMP = bytes.maketrans(b"ACGTNacgtn", b"TGCANtgcan")
BASE = {65: 0, 67: 1, 71: 2, 84: 3, 97: 0, 99: 1, 103: 2, 116: 3}   # A C G T; anything else -> 4
G = {}                                                                # shared with forked workers


def load(path):
    ids, lens, seqs, la, oa, lb, ob, woff = [], [], [], [], [], [], [], []
    off = 0
    with open(path, "rb") as f:
        for line in f:
            t = line[:2]
            if t == b"S\t":
                p = line.rstrip(b"\n").split(b"\t", 3)
                ids.append(int(p[1])); seqs.append(p[2]); lens.append(len(p[2]))
            elif t == b"L\t":
                p = line.split(b"\t", 6)
                la.append(int(p[1])); oa.append(p[2] == b"-"); lb.append(int(p[3])); ob.append(p[4] == b"-")
            elif t == b"W\t":
                woff.append(off)
            off += len(line)
    ids = np.array(ids, dtype=np.int64)
    lo, hi = int(ids.min()), int(ids.max())
    idmap = np.full(hi - lo + 1, -1, dtype=np.int64); idmap[ids - lo] = np.arange(len(ids))
    n = len(ids)
    lens = np.array(lens, dtype=np.int64)
    starts = np.zeros(n + 1, dtype=np.int64); np.cumsum(lens, out=starts[1:])
    buf = b"".join(seqs)
    first = np.array([BASE.get(s[0], 4) for s in seqs], dtype=np.int8)
    last_c = np.array([BASE.get(s[-1], 4) for s in seqs], dtype=np.int8)
    last_c = np.where(last_c < 4, 3 - last_c, 4).astype(np.int8)     # complement of the last base (A<->T, C<->G)
    fb = np.empty(2 * n, dtype=np.int8); fb[0::2] = first; fb[1::2] = last_c
    a = 2 * idmap[np.array(la) - lo] + np.array(oa, dtype=np.int64)
    b = 2 * idmap[np.array(lb) - lo] + np.array(ob, dtype=np.int64)
    src = np.concatenate([a, b ^ 1]); dst = np.concatenate([b, a ^ 1])  # each link and its reverse
    e = np.unique(src * (2 * n) + dst)
    src, dst = e // (2 * n), e % (2 * n)
    indeg = np.bincount(dst, minlength=2 * n)
    G.update(path=path, ids=ids, lo=lo, idmap=idmap, lens=lens, starts=starts, buf=buf, fb=fb, indeg=indeg, n=n)
    return src, dst, woff


def determinism(src, dst):
    fb, indeg, n = G["fb"], G["indeg"], G["n"]
    key = src * 8 + fb[dst]
    u, cnt = np.unique(key, return_counts=True)
    bad = u[cnt >= 2]
    base = bad % 8
    acgt = bad[base < 4]
    conf = np.zeros(2 * n, dtype=bool); conf[acgt // 8] = True
    # same-first-base siblings at least one of which is a merge (so the edge into it is a pebble)
    sel = np.isin(key, acgt) & (indeg[dst] >= 2)
    merge_src = np.unique(src[sel])
    return dict(oriented_edges=len(src), merges=int((indeg >= 2).sum()),
                nondet_sources=int(len(np.unique(acgt // 8))), nondet_sources_N=int(len(np.unique(bad[base >= 4] // 8))),
                nondet_with_merge_sibling=int(len(merge_src))), conf, acgt


def oname(o):
    return "%d%s" % (G["ids"][o >> 1], "-" if o & 1 else "+")


def olabel(o):
    i = o >> 1; s = G["buf"][G["starts"][i]:G["starts"][i + 1]]
    return s.translate(COMP)[::-1] if o & 1 else s


def window(o, i, k):
    """the k-1 bases before and the k-1 bases from the start of step i of the walk o"""
    left, j, got = [], i - 1, 0
    while j >= 0 and got < k - 1:
        s = olabel(o[j]); left.append(s); got += len(s); j -= 1
    right, j, got = [], i, 0
    while j < len(o) and got < k - 1:
        s = olabel(o[j]); right.append(s); got += len(s); j += 1
    return b"".join(reversed(left))[-(k - 1):] + b"|" + b"".join(right)[:k - 1]


def walk(task):
    off, k, do_b = task
    with open(G["path"], "rb") as f:
        f.seek(off); p = f.readline().rstrip(b"\n").split(b"\t")
    name = b"#".join(p[1:4]).decode()
    toks = re.findall(rb"[<>]\d+", p[6])
    ids = np.fromiter((int(t[1:]) for t in toks), dtype=np.int64, count=len(toks))
    rev = np.fromiter((t[0] == 60 for t in toks), dtype=np.int64, count=len(toks))
    o = 2 * G["idmap"][ids - G["lo"]] + rev
    n2 = 2 * G["n"]
    pos = np.concatenate([[0], np.cumsum(G["lens"][o >> 1])[:-1]])
    jstep = np.nonzero(G["indeg"][o[1:]] >= 2)[0] + 1                   # steps entered by a pebbled edge
    J = pos[jstep]
    d = np.diff(J)
    close = np.nonzero(d <= k - 2)[0]
    r = dict(name=name, steps=len(o), pebbles=len(jstep), pairs=len(d), close=len(close),
             close_ex=(name, oname(o[jstep[close[0]] - 1]) + ">" + oname(o[jstep[close[0]]]),
                       oname(o[jstep[close[0] + 1] - 1]) + ">" + oname(o[jstep[close[0] + 1]]), int(d[close[0]])) if len(close) else None)
    # D: same edge into a merge traversed twice by this walk
    e = o[:-1] * n2 + o[1:]
    ue, inv, cnt = np.unique(e, return_inverse=True, return_counts=True)
    rep_nodes = int((np.unique(o, return_counts=True)[1] >= 2).sum())
    rep = np.nonzero((cnt >= 2) & (G["indeg"][ue % n2] >= 2))[0]
    same, ex = 0, None
    for q in rep:
        steps = np.nonzero(inv == q)[0][:2] + 1
        w1, w2 = window(o, steps[0], k), window(o, steps[1], k)
        if w1 == w2:
            same += 1
            if ex is None:
                ex = (name, oname(int(ue[q] // n2)) + ">" + oname(int(ue[q] % n2)), w1.decode())
    r.update(rep_nodes=rep_nodes, rep_edges=len(rep), rep_same_window=same, rep_ex=ex)
    # C: which siblings of non-deterministic vertices this walk takes
    cs = np.nonzero(G["conf"][o[:-1]])[0]
    r["taken"] = np.unique(o[cs] * n2 + o[cs + 1])
    # B: windows of every pebble (sampled walks only)
    if do_b:
        h = np.fromiter((hash(window(o, i, k)) for i in jstep), dtype=np.int64, count=len(jstep))
        ed = o[jstep - 1] * n2 + o[jstep]
        r["bw"] = np.unique(np.stack([ed, h], axis=1), axis=0)
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("gfa"); ap.add_argument("-k", type=int, default=31)
    ap.add_argument("-j", type=int, default=8); ap.add_argument("--sample", type=int, default=12)
    a = ap.parse_args()
    src, dst, woff = load(a.gfa)
    det, conf, acgt = determinism(src, dst)
    G["conf"] = conf
    print("graph", dict(vertices=G["n"], walks=len(woff), **det), flush=True)
    # sampled haplotypes for B: the first --sample distinct sample#hap in file order
    names = []
    with open(a.gfa, "rb") as f:
        for off in woff:
            f.seek(off); names.append(b"#".join(f.read(200).split(b"\t")[1:3]))
    keep = list(dict.fromkeys(names))[:a.sample]
    tasks = [(off, a.k, nm in keep) for off, nm in zip(woff, names)]
    tot = dict(steps=0, pebbles=0, pairs=0, close=0, rep_nodes_walks=0, rep_edges=0, rep_edges_walks=0, rep_same_window=0)
    close_ex = rep_ex = None; taken = []; bw = []
    with Pool(a.j) as pool:
        for r in pool.imap_unordered(walk, tasks, chunksize=4):
            for x in ("steps", "pebbles", "pairs", "close", "rep_edges", "rep_same_window"):
                tot[x] += r[x]
            tot["rep_nodes_walks"] += r["rep_nodes"] > 0
            tot["rep_edges_walks"] += r["rep_edges"] > 0
            close_ex = close_ex or r["close_ex"]; rep_ex = rep_ex or r["rep_ex"]
            taken.append(r["taken"])
            if "bw" in r: bw.append(r["bw"])
    n2 = 2 * G["n"]
    taken = np.unique(np.concatenate(taken)) if taken else np.zeros(0, np.int64)
    ts, td = taken // n2, taken % n2
    k2 = np.unique(ts * 8 + G["fb"][td], return_counts=True)
    both = k2[0][(k2[1] >= 2) & (k2[0] % 8 < 4)]                        # two same-base siblings both walked
    both_merge = np.unique(ts[np.isin(ts * 8 + G["fb"][td], both) & (G["indeg"][td] >= 2)])
    print("A", dict(pebbles=tot["pebbles"], consecutive_pairs=tot["pairs"], within_k_minus_2=tot["close"],
                    share=round(tot["close"] / max(1, tot["pairs"]), 4), example=close_ex), flush=True)
    if bw:
        bw = np.unique(np.concatenate(bw), axis=0)
        e, c = np.unique(bw[:, 0], return_counts=True)
        print("B", dict(sampled_haplotypes=len(keep), pebbled_edges=len(e), edges_with_2plus_windows=int((c >= 2).sum()),
                        share=round(float((c >= 2).mean()), 4)), flush=True)
    print("C", dict(nondet_sources_walked=int(len(np.unique(both // 8))), of_which_with_merge_sibling=int(len(both_merge))), flush=True)
    print("D", dict(walks_revisiting_a_vertex=tot["rep_nodes_walks"], walks_repeating_a_merge_edge=tot["rep_edges_walks"],
                    repeated_merge_edges=tot["rep_edges"], with_identical_windows=tot["rep_same_window"], example=rep_ex), flush=True)


if __name__ == "__main__":
    main()
