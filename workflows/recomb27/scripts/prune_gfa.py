#!/usr/bin/env python3
"""Panel without some samples: prune_gfa.py IN.gfa OUT.gfa SAMPLE [SAMPLE ...]
or with only some haplotypes:  prune_gfa.py IN.gfa OUT.gfa --keep NAME [NAME ...]   (NAME = SAMPLE or SAMPLE#HAP)
Drops the walks (W lines, or P lines named SAMPLE#...) of the listed samples, then keeps only the segments and links that a
remaining walk still uses, so the result is the graph the remaining haplotypes alone would make (not a graph with orphaned
HG005 nodes and edges).  Two passes over the file; the node and edge sets stay small (chr21: a few million)."""
import sys, re
inp, out, drop = sys.argv[1], sys.argv[2], set(sys.argv[3:])
keep = drop - {"--keep"} if "--keep" in drop else None
TOK = re.compile(r"([<>])(\d+)")
nodes, edges = set(), set()
dropped = kept = 0

def sample(t):
    return t[1] if t[0] == "W" else t[1].split("#")[0]

def gone(t):    # is this walk left out?
    if keep is None:
        return sample(t) in drop
    hap = f"{t[1]}#{t[2]}" if t[0] == "W" else "#".join(t[1].split("#")[:2])
    return sample(t) not in keep and hap not in keep

def edge(e1, e2):   # canonical form of the arc (e1 -> e2) and its reverse complement; e = 2*id + (forward)
    return min((e1 << 32) | e2, ((e2 ^ 1) << 32) | (e1 ^ 1))

with open(inp) as f:
    for l in f:
        if l[0] not in "WP":
            continue
        t = l.split("\t", 7)
        if gone(t):
            dropped += 1
            continue
        kept += 1
        seq = t[6] if t[0] == "W" else t[2]
        if t[0] == "P":
            seq = re.sub(r"(\d+)([+-])", lambda m: (">" if m[2] == "+" else "<") + m[1], seq.replace(",", ""))
        prev = None
        for o, i in TOK.findall(seq):
            e = 2 * int(i) + (o == ">")
            nodes.add(int(i))
            if prev is not None:
                edges.add(edge(prev, e))
            prev = e
print(f"{inp}: {kept} walks kept, {dropped} dropped; {len(nodes)} nodes, {len(edges)} edges used", file=sys.stderr)
ns = nl = 0
with open(inp) as f, open(out, "w") as o:
    for l in f:
        c = l[0]
        if c == "S":
            if int(l.split("\t", 2)[1]) in nodes:
                o.write(l); ns += 1
        elif c == "L":
            t = l.split("\t", 6)
            if edge(2 * int(t[1]) + (t[2] == "+"), 2 * int(t[3]) + (t[4] == "+")) in edges:
                o.write(l); nl += 1
        elif c in "WP":
            if not gone(l.split("\t", 3)):
                o.write(l)
        else:
            o.write(l)
print(f"{out}: {ns} segments, {nl} links", file=sys.stderr)
