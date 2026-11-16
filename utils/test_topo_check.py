# Builds one small graph per pattern of topo_check.py (k=3) and checks that it finds each one.
import ast, os, subprocess, sys, tempfile
here = os.path.dirname(os.path.abspath(__file__))
tmp = tempfile.mkdtemp()
def gfa(name, lab, E, W):
    ids = {v: i + 1 for i, v in enumerate(lab)}
    L = ["H\tVN:Z:1.1"] + ["S\t%d\t%s" % (ids[v], s) for v, s in lab.items()]
    L += ["L\t%d\t+\t%d\t+\t0M" % (ids[a], ids[b]) for a, b in E]
    L += ["W\t%s\t0\tchr\t0\t0\t%s" % (c, "".join(">%d" % ids[v] for v in w)) for c, w in W.items()]
    p = os.path.join(tmp, name + ".gfa"); open(p, "w").write("\n".join(L) + "\n"); return p
def run(p):
    out = subprocess.run([sys.executable, os.path.join(here, "topo_check.py"), p, "-k", "3", "-j", "2"],
                         check=True, capture_output=True, text=True).stdout
    return {l.split(" ", 1)[0]: ast.literal_eval(l.split(" ", 1)[1]) for l in out.splitlines()}
A = run(gfa("exA", dict(s="G", a="A", b="C", u="T", z="G", v="T", t="A"),
            [("s","a"),("s","b"),("a","u"),("b","u"),("u","v"),("u","z"),("z","v"),("v","t")],
            dict(c1="sauvt", c2="sbuvt", c3="sauzvt")))
assert A["A"]["within_k_minus_2"] == 2, A["A"]
B = run(gfa("exB", dict(s="A", u="GG", x="CC", v="T", y="A", z="C", t="GG"),
            [("s","u"),("s","x"),("u","v"),("x","v"),("v","y"),("v","z"),("y","t"),("z","t")],
            dict(c1="suvyt", c2="suvzt", c3="sxvyt")))
assert B["B"]["edges_with_2plus_windows"] == 1, B["B"]
C = run(gfa("exC", dict(s="G", u="A", y="C", v="TA", w="TC", t="GG"),
            [("s","u"),("s","y"),("u","v"),("y","v"),("u","w"),("v","t"),("w","t")],
            dict(c1="suvt", c2="suwt", c3="syvt")))
assert C["graph"]["nondet_sources"] == 1 and C["graph"]["nondet_with_merge_sibling"] == 1, C["graph"]
assert C["C"] == dict(nondet_sources_walked=1, of_which_with_merge_sibling=1), C["C"]
D = run(gfa("exD", {"1": "T", "2": "T", "3": "G", "4": "C", "5": "TG", "6": "CA", "7": "A"},
            [("1","2"),("2","3"),("2","4"),("3","5"),("4","5"),("5","6"),("6","5"),("6","7")],
            dict(c1="123567", c2="124567", c3="1235656567")))
assert D["D"]["repeated_merge_edges"] == 1 and D["D"]["with_identical_windows"] == 1, D["D"]
assert D["graph"]["nondet_sources"] == 0
print("topo_check finds all four patterns"); print(D["D"]["example"], C["C"], A["A"]["example"])
