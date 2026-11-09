# Included last by the Snakefile. Rules declare mem_gb (profiles/local caps their sum with --resources mem_gb); the Slurm executor
# reads mem_mb, so each rule's mem_gb is copied into mem_mb (harmless locally). It also picks the Alliance partition: by-core,
# large-memory nodes above 16 GB per core (Narval refused 150 GB on 8 cores in cpubase), tier b2 (<= 12 h) for the rules
# profiles/slurm gives more than 3 h, else b1; left alone, the executor picked cpubase_bynode_b1, which refuses 1-core jobs.
import inspect

def _raw(x):        # Snakemake 8 keeps the int or function; 9 wraps it in a Resource (the function in _evaluator, unset = value None)
    if type(x).__name__ != "Resource":
        return x
    ev = getattr(x, "_evaluator", None)
    return ev if ev is not None else x.value

def _apply(g, f):   # f(gb) with the same arguments as the mem_gb function (wildcards, and input when it asks for it)
    if not callable(g):
        return f(g or 0)
    if "input" in inspect.signature(g).parameters:
        return lambda wildcards, input: f(g(wildcards, input))
    return lambda wildcards: f(g(wildcards))

_LONG = {"g_dist", "giraffe", "g_minzip", "vg_hapl", "vg_dist", "build_index", "rank_call", "wg_merge"}

for _r in workflow.rules:
    _g = _raw(_r.resources.get("mem_gb"))
    if _g is not None:
        _r.resources["mem_mb"] = _apply(_g, lambda gb: 1024 * gb)
    _c = _raw(_r.resources.get("_cores"))
    _c = _c if isinstance(_c, int) else 1
    _tier = "b2" if _r.name in _LONG else "b1"
    _r.resources["slurm_partition"] = _apply(_g, lambda gb, c=_c, t=_tier: f"cpu{'large' if gb / max(1, c) > 16 else 'base'}_bycore_{t}")
