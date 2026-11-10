# Included last by the Snakefile. Rules declare mem_gb (profiles/local caps their sum with --resources mem_gb); the Slurm executor
# reads mem_mb, so each rule's mem_gb is copied into mem_mb (harmless locally).
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

for _r in workflow.rules:
    _g = _raw(_r.resources.get("mem_gb"))
    if _g is not None:
        _r.resources["mem_mb"] = _apply(_g, lambda gb: 1024 * gb)
