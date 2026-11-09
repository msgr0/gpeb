# Personalized graphs (paper: Results "Personalized graphs with gpeb cover", Figure 4, Table 1's vg haplotypes side; supplement
# "Personalized graphs", read classes and bases). From each rung of HG002's HiFi reads, on noHG002: the cover with at most N colors
# per window, the graph of vg haplotypes with N haplotypes and its diploid sample; controls N random haplotypes and the whole
# noHG002 graph. The same reads are mapped on each graph with vg giraffe.
#   results/exp3/<chrom>.tsv           one row per (graph, rung): read classes, MAPQ bins, edits (edits.py)
#   results/exp3/<chrom>.paired.tsv    the cover with N colors against each other graph on the same reads, read by read
#   results/exp3/<chrom>.extra.tsv, .paired_extra.tsv   bases: unmapped, clipped and unmatched (exp3_extra.py)
# Graph ids: vg.<d32|sN>.ci<ci>.c<cov> and gpeb.cover<N>.c<cov> (built from rung cov, mapped with it), rand.n<N>.s<seed> and
# fullpanel (no reads; mapped with every rung).
E3, P3 = f"{W}/exp3", "noHG002"
G = E3 + "/{chrom}/graphs"
VP = W + "/vgpanel/{chrom}/panel"
VGCFG = {"d32": "--preset diploid", **{f"s{n}": f"--preset default --num-haplotypes {n}" for n in (2, 4, 8, 16, 32, 64)}}  # --preset FIRST

wildcard_constraints:
    gid = r"(vg|gpeb|rand)\.[A-Za-z0-9.]+|fullpanel"

def hifi(chrom, cov):
    return read_set("hifi", chrom, "HG002", cov)

def gids(chrom):
    """(graph id, rung mapped on it) for every arm of one chromosome"""
    u = []
    for c in COVS:
        u += [(f"vg.{v}.ci{ci}.c{c}", c) for v in config["vg_cfgs"] for ci in config["kmc_ci"]]
        u += [(f"gpeb.cover{n}.c{c}", c) for n in config["cover_caps"]]
        u += [(f"rand.n{n}.s{s}", c) for n in config["random_ns"] for s in config["random_seeds"]] + [("fullpanel", c)]
    return u

rule personalize:
    input: expand(f"{W}/results/exp3/{{chrom}}.{{t}}.tsv", chrom=CH, t=["paired", "extra"]),
           expand(f"{W}/idx/{{chrom}}.{P3}.{D}.pex", chrom=CH)      # Table 1: the pebble index of noHG002
    output: touch(f"{W}/results/exp3.done")
    shell: "{BENCH}"

rule exp3_chrom:
    input: maps = lambda wc: [f"{E3}/{wc.chrom}/map/{c}x/{g}.stats" for g, c in gids(wc.chrom)]
    output: f"{W}/results/exp3/{{chrom}}.tsv", f"{W}/results/exp3/{{chrom}}.paired.tsv"
    shell: "python3 -I {S}/exp3_aggregate.py {W}/results/exp3/{wildcards.chrom} {wildcards.chrom} {E3}/{wildcards.chrom}/map"

rule exp3_extra:
    input: rules.exp3_chrom.output, reads = lambda wc: [hifi(wc.chrom, c) for c in COVS]
    output: f"{W}/results/exp3/{{chrom}}.extra.tsv", f"{W}/results/exp3/{{chrom}}.paired_extra.tsv"
    resources: mem_gb = 16
    shell: "python3 -I {S}/exp3_extra.py {wildcards.chrom} {E3}/{wildcards.chrom}/map {W}/reads/hifi/{wildcards.chrom} {W}/results/exp3/{wildcards.chrom}"

# ---- the panel as vg needs it: GBZ, distance index, r-index, haplotype index (Table 1; one tool per rule: its own time and RSS) ----
rule vg_gbz:
    input: VGIN, gfa = f"{W}/panel/{{chrom}}.{P3}.gfa"
    output: f"{VP}.gbz"
    threads: config["threads_vg"]
    resources: mem_gb = lambda wc: config["mem_gb_vgpanel"][wc.chrom] // 2
    params: b = bench("vg_gbz", "panel"), sz = sizes("vg_gbz", "panel")
    log: f"{W}/logs/vg_gbz.{{chrom}}.log"
    shell: "{TIME} {params.b} {VG} gbwt -G {input.gfa} --gbz-format -g {output} --set-reference GRCh38 --set-reference CHM13 2> {log} && {STAT} {output} > {params.sz}"

rule vg_dist:
    input: VGIN, gbz = f"{VP}.gbz"
    output: f"{VP}.dist"
    threads: config["threads_vg"]
    resources: mem_gb = lambda wc: config["mem_gb_vgpanel"][wc.chrom]
    params: b = bench("vg_dist", "panel"), sz = sizes("vg_dist", "panel"), tmp = f"{W}/tmp/vg_dist.{{chrom}}"
    log: f"{W}/logs/vg_dist.{{chrom}}.log"
    shell: "mkdir -p {params.tmp} && {TIME} {params.b} {VG} index -j {output} -t {threads} -b {params.tmp} {input.gbz} 2> {log} && rm -rf {params.tmp} && {STAT} {output} > {params.sz}"

rule vg_ri:
    input: VGIN, gbz = f"{VP}.gbz"
    output: f"{VP}.ri"
    threads: config["threads_vg"]
    resources: mem_gb = lambda wc: config["mem_gb_vgpanel"][wc.chrom] // 2
    params: b = bench("vg_ri", "panel"), sz = sizes("vg_ri", "panel")
    log: f"{W}/logs/vg_ri.{{chrom}}.log"
    shell: "{TIME} {params.b} {VG} gbwt -Z {input.gbz} -r {output} 2> {log} && {STAT} {output} > {params.sz}"

rule vg_hapl:
    input: VGIN, gbz = f"{VP}.gbz", dist = f"{VP}.dist", ri = f"{VP}.ri"
    output: f"{VP}.hapl"
    threads: config["threads_vg"]
    resources: mem_gb = lambda wc: config["mem_gb_vgpanel"][wc.chrom] // 2
    params: b = bench("vg_hapl", "panel"), sz = sizes("vg_hapl", "panel")
    log: f"{W}/logs/vg_hapl.{{chrom}}.log"
    shell: "{TIME} {params.b} {VG} haplotypes -d {input.dist} -r {input.ri} -H {output} -t {threads} {input.gbz} 2> {log} && {STAT} {output} > {params.sz}"

# ---- vg arms: KMC counts of rung cov, then haplotype sampling ----------------------------------------------------------------
rule kmc:
    input: lambda wc: hifi(wc.chrom, wc.cov)
    output: f"{E3}/{{chrom}}/kmc/c{{cov}}.ci{{ci}}.kff"
    threads: 8
    resources: mem_gb = 16
    params: pre = lambda wc, output: output[0][:-4], tmp = f"{W}/tmp/kmc.{{chrom}}.{{cov}}.{{ci}}", b = bench("kmc", "c{cov}.ci{ci}")
    shell:
        "mkdir -p {params.tmp} && {TIME} {params.b} {KMC} -k29 -m16 -ci{wildcards.ci} -okff -t{threads} -fq {input} {params.pre} {params.tmp} > /dev/null "
        "&& rm -rf {params.tmp}"

rule vg_hap:        # --coverage = the rung, rounded and at least 1: vg's own estimate fails at low depth (its logs say "unreliable")
    input: VGIN, kff = f"{E3}/{{chrom}}/kmc/c{{cov}}.ci{{ci}}.kff", gbz = f"{VP}.gbz", dist = f"{VP}.dist", ri = f"{VP}.ri", hapl = f"{VP}.hapl"
    output: temp(f"{G}/vg.{{vcfg}}.ci{{ci}}.c{{cov}}.gbz")
    threads: config["threads_vg"]
    resources: mem_gb = lambda wc: config["mem_gb_vgpanel"][wc.chrom] // 4
    params: opts = lambda wc: VGCFG[wc.vcfg], k = lambda wc: max(1, round(float(wc.cov))),
            b = bench("vg_hap", "{vcfg}.ci{ci}.c{cov}"), sz = sizes("vg_hap", "{vcfg}.ci{ci}.c{cov}")
    log: f"{W}/logs/vg_hap.{{chrom}}.{{vcfg}}.ci{{ci}}.c{{cov}}.log"
    shell:
        "{LOCAL}D=$(L {input.dist}); {TIME} {params.b} {VG} haplotypes -t {threads} -v 2 -i {input.hapl} -d $D -r {input.ri} -k {input.kff} -g {output} "
        "{params.opts} --include-reference --coverage {params.k} {input.gbz} 2> {log} && {STAT} {output} > {params.sz}"

# ---- gpeb arms: the cover of rung cov on the noHG002 index, spliced into a GFA -------------------------------------------------
rule gpeb_cover:
    input: gpeb = rules.gpeb.output.gpeb, idx = f"{W}/idx/{{chrom}}.{P3}.{D}.pex", gfa = f"{W}/panel/{{chrom}}.{P3}.gfa",
           reads = lambda wc: hifi(wc.chrom, wc.cov)
    output: gfa = temp(f"{G}/gpeb.cover{{cap}}.c{{cov}}.gfa"), sets = f"{E3}/{{chrom}}/cover/cover{{cap}}.c{{cov}}.sets.tsv"
    threads: config["threads_query"]
    resources: mem_gb = lambda wc: config["mem_gb_build"][wc.chrom] // 2     # measured 9.7 GB on chr9
    params: pfx = lambda wc, output: output.sets[:-len(".sets.tsv")], b = bench("gpeb_cover", "cover{cap}.c{cov}")
    log: f"{W}/logs/gpeb_cover.{{chrom}}.cover{{cap}}.c{{cov}}.log"
    shell:
        "{TIME} {params.b} {input.gpeb} --threads {threads} --use-index {input.idx} {input.gfa} --classify {IX[k]} {input.reads} --set-call cover "
        "--cover-max {wildcards.cap} --cover-window {config[cover_window]} --cover-min-frac {config[cover_min_frac]} --mosaic-keep {REFS} "
        "--mosaic-out {params.pfx} --mosaic-gfa {output.gfa} > /dev/null 2> {log}"

# ---- random control: N whole haplotypes (references kept), drawn with a seed made of the chromosome, N and the seed -------------
rule random_gfa:
    input: f"{W}/panel/{{chrom}}.{P3}.gfa"
    output: gfa = temp(f"{G}/rand.n{{n}}.s{{seed}}.gfa"), picks = f"{E3}/{{chrom}}/random/n{{n}}.s{{seed}}.txt"
    params: refs = REFS.replace(",", " ")
    log: f"{W}/logs/random_gfa.{{chrom}}.n{{n}}.s{{seed}}.log"
    shell:
        """
        awk -F'\\t' '$1=="W" && $2!="CHM13" && $2!="GRCh38" {{print $2"#"$3}}' {input} | sort -u \
          | python3 -I -c 'import sys, random; h = sorted(l.strip() for l in sys.stdin); print("\\n".join(random.Random("{wildcards.chrom}:{wildcards.n}:{wildcards.seed}").sample(h, {wildcards.n})))' > {output.picks}
        python3 -I {S}/prune_gfa.py {input} {output.gfa} --keep {params.refs} $(cat {output.picks}) 2> {log}
        """

rule gfa2gbz:       # the cover and random graphs come as GFA; the conversion is part of their cost
    input: VGIN, gfa = f"{G}/{{gid}}.gfa"
    output: temp(f"{G}/{{gid}}.gbz")
    wildcard_constraints: gid = r"(gpeb|rand)\.[A-Za-z0-9.]+"
    threads: 4
    resources: mem_gb = 30
    params: b = bench("gfa2gbz", "{gid}"), sz = sizes("gfa2gbz", "{gid}")
    log: f"{W}/logs/gfa2gbz.{{chrom}}.{{gid}}.log"
    shell: "{TIME} {params.b} {VG} gbwt -G {input.gfa} --gbz-format -g {output} --set-reference GRCh38 --set-reference CHM13 2> {log} && {STAT} {output} > {params.sz}"

# ---- giraffe indexes of every graph (the whole noHG002 graph reuses the panel's GBZ and distance index), and the mapping -----------
def gbz_of(wc):
    return VP.format(chrom=wc.chrom) + ".gbz" if wc.gid == "fullpanel" else G.format(chrom=wc.chrom) + f"/{wc.gid}.gbz"

def dist_of(wc):
    return VP.format(chrom=wc.chrom) + ".dist" if wc.gid == "fullpanel" else G.format(chrom=wc.chrom) + f"/{wc.gid}.dist"

def arm_mem(table, chrom, gid, default):   # arm = the graph id without its rung (gpeb.cover32.c0.5 -> gpeb.cover32); chrom.arm > arm > chrom
    arm, m = re.sub(r"\.c[0-9.]+$", "", gid), config[table]
    return m.get(f"{chrom}.{arm}", m.get(arm, m.get(chrom, default)))

rule g_dist:
    input: VGIN, gbz = f"{G}/{{gid}}.gbz"
    output: temp(f"{G}/{{gid}}.dist")
    wildcard_constraints: gid = r"(vg|gpeb|rand)\.[A-Za-z0-9.]+"
    threads: 8
    resources: mem_gb = lambda wc: arm_mem("mem_gb_gdist", wc.chrom, wc.gid, 40)
    params: b = bench("g_dist", "{gid}"), sz = sizes("g_dist", "{gid}"), tmp = f"{W}/tmp/g_dist.{{chrom}}.{{gid}}"
    log: f"{W}/logs/g_dist.{{chrom}}.{{gid}}.log"
    shell: "mkdir -p {params.tmp} && {TIME} {params.b} {VG} index -j {output} -t {threads} -b {params.tmp} {input.gbz} 2> {log} && rm -rf {params.tmp} && {STAT} {output} > {params.sz}"

rule g_minzip:      # HiFi: k31 w50 weighted (giraffe's long-read defaults); one command writes the minimizer and zipcode indexes
    input: VGIN, gbz = gbz_of, dist = dist_of
    output: m = temp(f"{G}/{{gid}}.min"), z = temp(f"{G}/{{gid}}.zip")
    threads: 8
    resources: mem_gb = lambda wc: arm_mem("mem_gb_minzip", wc.chrom, wc.gid, 30)
    params: b = bench("g_minzip", "{gid}"), sz = sizes("g_minzip", "{gid}")
    log: f"{W}/logs/g_minzip.{{chrom}}.{{gid}}.log"
    shell: "{LOCAL}D=$(L {input.dist}); {TIME} {params.b} {VG} minimizer -d $D -z {output.z} -k 31 -w 50 -W -t {threads} -o {output.m} {input.gbz} 2> {log} && {STAT} {output} > {params.sz}"

rule giraffe:       # the GAF only lives inside the rule: edits.py keeps one line per read (edits core clip mapped maxindel mapq)
    input: VGIN, gbz = gbz_of, dist = dist_of, m = f"{G}/{{gid}}.min", z = f"{G}/{{gid}}.zip", reads = lambda wc: hifi(wc.chrom, wc.cov)
    output: per = f"{E3}/{{chrom}}/map/{{cov}}x/{{gid}}.peread.gz", stats = f"{E3}/{{chrom}}/map/{{cov}}x/{{gid}}.stats"
    threads: config["threads_map"]
    resources: mem_gb = lambda wc: arm_mem("mem_gb_giraffe", wc.chrom, wc.gid, 32)
    params: b = bench("giraffe", "{gid}.on{cov}"), gaf = f"{W}/tmp/giraffe.{{chrom}}.{{gid}}.{{cov}}.gaf"
    log: f"{W}/logs/giraffe.{{chrom}}.{{gid}}.{{cov}}.log"
    shell:
        """
        {LOCAL}D=$(L {input.dist})
        {TIME} {params.b} {VG} giraffe -b hifi -Z {input.gbz} -d $D -m {input.m} -z {input.z} -f {input.reads} -t {threads} -o gaf > {params.gaf} 2> {log}
        python3 -I {S}/edits.py $(( $(wc -l < {input.reads}) / 4 )) {params.gaf} {output.per}.tmp {input.reads} > {output.stats}
        gzip -c {output.per}.tmp > {output.per} && rm -f {output.per}.tmp {params.gaf}
        """
