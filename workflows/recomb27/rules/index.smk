# Index study (paper: Results "The pebble index", Table 1's pebble-index side, supplement "Colors of the records" and "Index
# parameters sweep"). On the full graph: the default index and a one-at-a-time sweep around it (min-gap and syncmer rebuilds, two
# indexes cut from the default one that keep the k-mers held by at most 1 or 10 colors), each scored by the rank calls of HG002 and
# HG005 and by the cover of the parents' reads at 1x and 10x (cover-min-frac swept on the default index).
#   results/exp1/<chrom>.index.tsv   one row per index: the build's own counts, bytes, time and peak memory
#   results/exp1/<chrom>.calls.tsv   rank and cover rows (call_aggregate.py) for every index, sample and coverage
#   results/exp1/<chrom>.rare.tsv    the default index's records by number of colors (pex_rare)
X1 = config["exp1"]
CFGS = sorted({f"mg{g}.s{IX['syncmer']}" for g in X1["min_gaps"]} | {f"mg{IX['min_gap']}.s{s}" for s in X1["syncmers"]} | {D})
DERIVED = [f"{D}.max{n}" for n in X1["max_colours"]]

def exp1_units(chrom):
    u = [f"{W}/calls/hifi/{chrom}/full/{c}/rank.{s}.G{v}x/log" for c in CFGS + DERIVED for s in SM for v in X1["call_covs"]]
    u += [f"{W}/calls/hifi/{chrom}/full/{c}/cover.{s}.G{v}x/log" for c in CFGS + DERIVED for s in SM for v in X1["call_covs"]]
    u += [f"{W}/calls/hifi/{chrom}/full/{D}/coverF{f}.{s}.G{v}x/log" for f in X1["min_fracs"] for s in SM for v in X1["call_covs"]]
    return u

rule index_study:
    input: expand(f"{W}/results/exp1/{{chrom}}.{{t}}.tsv", chrom=CH, t=["calls", "index", "rare"])
    output: touch(f"{W}/results/exp1.done")
    shell: "{BENCH}"

rule exp1_calls:
    input: lambda wc: exp1_units(wc.chrom)
    output: f"{W}/results/exp1/{{chrom}}.calls.tsv"
    params: dirs = lambda wc, input: call_dirs(input)
    shell: "python3 -I {S}/call_aggregate.py {output} {EXPECT} {params.dirs}"

rule exp1_index:
    input: expand(f"{W}/idx/{{{{chrom}}}}.full.{{c}}.pex", c=CFGS + DERIVED)
    output: f"{W}/results/exp1/{{chrom}}.index.tsv"
    shell: "python3 -I {S}/index_stats.py {W} {wildcards.chrom} full {output}"

rule max_colours:   # the default index's records held by at most n colors (n = 1: the private ones), no rebuild
    input: gpeb = rules.gpeb.output.gpeb, idx = f"{W}/idx/{{chrom}}.{{panel}}.{D}.pex"
    output: f"{W}/idx/{{chrom}}.{{panel}}.{D}.max{{n}}.pex"
    resources: mem_gb = lambda wc: config["mem_gb_build"][wc.chrom] // 2
    params: b = bench("max_colours", "{panel}.max{n}"), sz = sizes("max_colours", "{panel}.max{n}")
    log: f"{W}/logs/index.{{chrom}}.{{panel}}.{D}.max{{n}}.log"
    shell: "{TIME} {params.b} {input.gpeb} --use-index {input.idx} --max-colours {wildcards.n} {output} 2> {log} && {STAT} {output} > {params.sz}"

rule rare:          # --max above any color count: the whole histogram (stderr); the record list itself is not kept
    input: rare = rules.gpeb.output.rare, idx = f"{W}/idx/{{chrom}}.full.{D}.pex"
    output: f"{W}/results/exp1/{{chrom}}.rare.tsv"
    resources: mem_gb = lambda wc: config["mem_gb_build"][wc.chrom] // 2
    shell: "{input.rare} --max 1000 {input.idx} > /dev/null 2> {output}"

def cover_frac(wc):  # cover = the configured --cover-min-frac; coverF<f> = the sweep
    return config["cover_min_frac"] if wc.kind == "cover" else wc.kind[len("coverF"):]

rule cover_call:    # per 200 kb window, the colors whose union carries the reads' hits; scored on the parents' transmitted haplotype
    input: gpeb = rules.gpeb.output.gpeb, idx = f"{W}/idx/{{chrom}}.{{panel}}.{{cfg}}.pex",
           reads = lambda wc: read_set(wc.reads, wc.chrom, wc.sample, wc.cov)
    output: f"{CALLS}/{{kind}}.{{sample}}.G{{cov}}x/log"
    wildcard_constraints: kind = r"cover(F[0-9.]+)?"
    threads: config["threads_query"]
    resources: mem_gb = call_mem
    params: d = lambda wc, output: os.path.dirname(output[0]), f = cover_frac,
            b = bench("cover_call", "{reads}.{panel}.{cfg}.{kind}.{sample}.{cov}")
    shell:
        "{TIME} {params.b} {input.gpeb} --threads {threads} --use-index {input.idx} --classify {IX[k]} {input.reads} --set-call cover "
        "--cover-max {config[call_cover_max]} --cover-window {config[cover_window]} --cover-min-frac {params.f} --mosaic-keep {REFS} "
        "--mosaic-out {params.d}/cover > /dev/null 2> {output}"
