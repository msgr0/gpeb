"""conf.py: the call's confidence over the sample sets of its reads, as gpeb's set-call[conf] line (gpeb.c call_conf), for the
whole-genome calls that wg_sets.py merges in Python. conf.py --selftest
Each read's set S(r) = the samples with a colour among its most-voted colours V(r); the call t = the sample holding the most reads. For every other sample j, a_j = reads t holds and j does not,
b_j = reads j holds and t does not (s_t - s_j = (a_j - b_j) / |R+|):
  D    = min_j (a_j - b_j) / (a_j + b_j), the call's margin over its hardest competitor among the reads that tell the two apart:
         0 = some sample is indistinguishable from the call, 1 = no sample holds a read the call does not, 1/2 = 3:1 for the call;
  D_lo = the same with each a_j / (a_j + b_j) replaced by its 95% Wilson lower bound (few discriminating reads -> lower);
  pi   = the share of the reads from each sample, by EM on the mixture P(r) = sum_{c in S(r)} pi_c; pi_t = the share the call explains.
row(DIR, counts): DIR = .../<chrom>/<panel>/<cfg>/rank.<sample>.G<cov>x names the call."""
import sys, re, math
from collections import Counter
import numpy as np


def wilson_lo(k, n, z=1.96):
    p = k / n
    return (p + z * z / (2 * n) - z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n))) / (1 + z * z / n)


def em(c, tol=1e-7, cap=5000):   # dense sets x samples matrix, fine up to ~1e5 distinct sets; sparse if a call outgrows it
    names = sorted({s for ss in c for s in ss}); ix = {s: i for i, s in enumerate(names)}
    M = np.zeros((len(c), len(names))); w = np.empty(len(c))   # float64: 0/1 entries, so the same sums as float32, without a cast per iteration
    for k, (ss, m) in enumerate(c.items()):
        M[k, [ix[s] for s in ss]] = 1; w[k] = m
    pi = np.full(len(names), 1 / len(names))
    for it in range(1, cap + 1):
        new = pi * (M.T @ (w / (M @ pi))) / w.sum()
        if np.abs(new - pi).max() < tol: return dict(zip(names, new)), it
        pi = new
    return dict(zip(names, pi)), cap


def conf(c):
    held, both = Counter(), Counter()
    for ss, m in c.items():
        for s in ss: held[s] += m
    t = max(sorted(held), key=held.get)
    for ss, m in c.items():
        if t in ss:
            for s in ss: both[s] += m
    D, Dlo, jmin, ab = 2.0, 2.0, "", (0, 0)
    for j in held:
        if j == t: continue
        a, b = held[t] - both[j], held[j] - both[j]
        d = (a - b) / (a + b) if a + b else 0.0
        lo = 2 * wilson_lo(a, a + b) - 1 if a + b else 0.0
        if d < D: D, jmin, ab = d, j, (a, b)
        Dlo = min(Dlo, lo)
    if D == 2.0: D = Dlo = 1.0   # no other sample holds any read
    pi, it = em(c)
    return t, held, D, Dlo, jmin, ab, pi, it




def row(d, c):   # the output line of the call in DIR d from its set counts c (wg_sets.py builds c from interned sets)
    chrom, panel, cfg, leaf = d.rstrip("/").split("/")[-4:]
    sample, cov = re.match(r"rank\.(HG00\d)\.G([\d.]+)x(?:\.r0)?$", leaf).groups()
    n = sum(c.values())
    t, held, D, Dlo, jmin, (a, b), pi, it = conf(c)
    r2 = sorted((s for s in held if s != t), key=lambda s: -held[s])[0]
    p = sorted(pi.items(), key=lambda kv: -kv[1])
    return (f"{chrom}\t{panel}\t{cfg}\t{sample}\t{cov}\t{n}\t{len(c)}\t{t}\t{held[t] / n:.4f}\t{r2}\t{held[r2] / n:.4f}\t{D:.4f}\t{Dlo:.4f}\t"
            f"{jmin}\t{a}\t{b}\t{pi[t]:.4f}\t{p[0][0]}\t{p[0][1]:.4f}\t{p[1][0]}\t{p[1][1]:.4f}\t{it}\n")


def selftest():
    c = Counter({frozenset("AB"): 10, frozenset("A"): 6, frozenset("B"): 2, frozenset("C"): 1})
    t, held, D, Dlo, jmin, ab, pi, it = conf(c)
    assert (t, jmin, ab, D) == ("A", "B", (6, 2), 0.5), (t, jmin, ab, D)    # B: a = 16 - 10, b = 12 - 10; C: 15/17
    assert Dlo < D and abs(sum(pi.values()) - 1) < 1e-6 and pi["A"] > pi["B"] > pi["C"] > 0, (Dlo, pi)
    c = Counter({frozenset("AB"): 5})                                        # B holds exactly A's reads: indistinguishable
    assert conf(c)[2] == 0.0
    print("selftest ok")


if __name__ == "__main__":
    sys.exit(selftest())
