# Shared rules: gpeb from this repository, the panels and their default index, the per-chromosome read sets (HiFi and short
# reads sliced from GIAB's aligned BAMs, each cut into a nested coverage ladder) and the rank call of a read set on an index.
# Every gpeb / vg / kmc command runs under GNU time: bench/<rule>__<chrom>__<key>.time = "wall_s  maxrss_kb  user_s  sys_s";
# every index or graph gets its size in sizes/<rule>__<chrom>__<key>.tsv ("path  bytes"); rule bench gathers both.
import math, os, re

W = os.path.abspath(config["workdir"])                         # relative paths: to the directory snakemake runs in
GFA = os.path.abspath(config["gfa_dir"])
CH, SM, COVS = config["chroms"], config["samples"], config["covs"]
PANELS, IX, REFS = config["panels"], config["index"], config["refs"]
D = f"mg{IX['min_gap']}.s{IX['syncmer']}"                      # the default index configuration
EXPECT = ",".join(f"{s}:{d}" for s, d in config["expect"].items())
REPO = os.path.normpath(os.path.join(workflow.basedir, "../.."))
S = os.path.join(workflow.basedir, "scripts")
VG = config["vg"] or f"{W}/bin/vg"                             # empty: the vg release named by vg_version, downloaded
VGIN = [] if config["vg"] else [f"{W}/bin/vg"]                  # ... and then an input of every rule that runs it
KMC = config["kmc"] or "kmc"
TIME = config["time_bin"] + " -f '%e\\t%M\\t%U\\t%S' -o"
# On a Slurm node vg mmaps its distance index, which crawls on a shared filesystem: copy it to the node's disk first.
# No-op where SLURM_TMPDIR is unset. Use: {LOCAL}D=$(L {input.dist}); vg ... -d $D   (cp -p: giraffe refuses a .dist newer than its .min)
LOCAL = 'L() { if [ -n "${SLURM_TMPDIR:-}" ]; then cp -p "$1" "$SLURM_TMPDIR/" && echo "$SLURM_TMPDIR/$(basename "$1")"; else echo "$1"; fi; }; '
bench = lambda name, key: f"{W}/bench/{name}__{{chrom}}__{key}.time"
sizes = lambda name, key: f"{W}/sizes/{name}__{{chrom}}__{key}.tsv"
STAT = "stat -c '%n\t%s'"

wildcard_constraints:
    chrom = "chr[0-9XY]+", sample = "HG[0-9]+", cov = r"[0-9.]+", panel = "|".join(PANELS), cfg = r"mg\d+\.s\d+(\.max\d+)?",
    reads = "hifi|sr", ci = r"\d", cap = r"\d+", vcfg = r"[ds]\d+", n = r"\d+", seed = r"\d+", mate = "[12]"

onstart:
    shell(f"mkdir -p {W}/bench {W}/sizes {W}/logs {W}/tmp")

# ---- tools --------------------------------------------------------------------------------------------------------------------
rule gpeb:          # gpeb and pex_rare from this repository at config gpeb_commit (git archive + the gfatools submodule)
    output: gpeb = f"{W}/bin/gpeb", rare = f"{W}/bin/pex_rare", commit = f"{W}/bin/COMMIT"
    params: src = f"{W}/src", rev = config["gpeb_commit"]
    threads: 4
    shell:
        """
        rm -rf {params.src} && mkdir -p {params.src}/libs/gfatools
        git -C {REPO} archive {params.rev} | tar -x -C {params.src}
        git -C {REPO}/libs/gfatools archive $(git -C {REPO} rev-parse {params.rev}:libs/gfatools) | tar -x -C {params.src}/libs/gfatools
        make -C {params.src} -j {threads} gpeb NATIVE=
        cp {params.src}/gpeb {output.gpeb} && cc -O2 -o {output.rare} {params.src}/utils/pex_rare.c
        {{ git -C {REPO} rev-parse {params.rev}; git -C {REPO} rev-parse {params.rev}:libs/gfatools; sha256sum {output.gpeb}; }} > {output.commit}
        """

rule vg:            # only when config vg is empty: the static release binary
    output: f"{W}/bin/vg"
    params: url = f"https://github.com/vgteam/vg/releases/download/{config['vg_version']}/vg"
    retries: 3
    shell: "curl -fsSL -o {output} {params.url} && chmod +x {output} && {output} version | head -1"

# ---- panels: the released GFA, or the graph the remaining haplotypes alone make ---------------------------------------------
rule panel_full:
    input: lambda wc: f"{GFA}/{wc.chrom}.gfa"
    output: f"{W}/panel/{{chrom}}.{{panel}}.gfa"
    wildcard_constraints: panel = "|".join(p for p, d in PANELS.items() if not d) or "^$"
    shell: "ln -sf {input} {output}"

rule panel_pruned:  # single-threaded, <= 2.4 GiB; chr21 ~15 min, chr2 ~75 min. protected(): a run that would redo it stops instead
    input: lambda wc: f"{GFA}/{wc.chrom}.gfa"
    output: protected(f"{W}/panel/{{chrom}}.{{panel}}.gfa")
    wildcard_constraints: panel = "|".join(p for p, d in PANELS.items() if d) or "^$"
    params: drop = lambda wc: " ".join(PANELS[wc.panel]), b = bench("panel_pruned", "{panel}")
    log: f"{W}/logs/panel.{{chrom}}.{{panel}}.log"
    shell: "{TIME} {params.b} python3 -I {S}/prune_gfa.py {input} {output}.part {params.drop} 2> {log} && mv {output}.part {output}"

# ---- the pebble index: the default one, and the min-gap / syncmer sweep of index.smk --------------------------------------------
def cfg_flags(cfg):
    mg, s = re.match(r"mg(\d+)\.s(\d+)$", cfg).groups()
    return f"--min-gap {mg} --syncmer {s}" + ("" if cfg == D else " --no-gfai")   # only the default build writes the GFA's .gfai

def index_mem(wc):  # the measured peaks scale with the GFA; min-gap 0 and syncmer 0 keep about twice the records
    return config["mem_gb_build"][wc.chrom] * (2 if re.match(r"mg0\.|.*\.s0$", wc.cfg) else 1)

rule build_index:
    input: gpeb = rules.gpeb.output.gpeb, gfa = f"{W}/panel/{{chrom}}.{{panel}}.gfa"
    output: f"{W}/idx/{{chrom}}.{{panel}}.{{cfg}}.pex"
    wildcard_constraints: cfg = r"mg\d+\.s\d+"
    threads: config["threads_build"]
    resources: mem_gb = index_mem
    params: flags = lambda wc: cfg_flags(wc.cfg), b = bench("index", "{panel}.{cfg}"), sz = sizes("index", "{panel}.{cfg}")
    log: f"{W}/logs/index.{{chrom}}.{{panel}}.{{cfg}}.log"
    shell:
        "{TIME} {params.b} {input.gpeb} --index {IX[k]} {output}.part {params.flags} --threads {threads} --mosaic-ref {IX[mosaic_ref]} "
        "{input.gfa} 2> {log} && mv {output}.part {output} && {STAT} {output} > {params.sz}"

# ---- reads: one chromosome sliced from GIAB's aligned BAM (no full download), then the nested ladder ---------------------------
# HiFi: CHM13v2.0 BAMs, kept as FASTQ (vg giraffe maps them in personalize.smk); short reads: GRCh38 BAMs, gzipped FASTA, mates
# named alike (samtools fasta -n) so that a pair lands in one rung, each mate classified as a read. Depth is measured on the BAM's
# reference. Rung c = the reads with crc32(name) / 2^32 < c / depth: one draw, nested.
def bam_url(reads, s):
    return f"{config['giab_base']}/{config[reads + '_bam'][s]}"

def read_set(reads, chrom, sample, cov):
    return f"{W}/reads/{reads}/{chrom}/{sample}.G{cov}x." + ("fq" if reads == "hifi" else "fa.gz")

rule chrlen:
    output: f"{W}/reads/{{reads}}/{{chrom}}.len"
    params: u = lambda wc: bam_url(wc.reads, SM[0])
    retries: 4
    resources: net = 1          # NCBI refuses many opens at once: at most `net` (profile) streams
    shell: "samtools view -H '{params.u}##idx##{params.u}.bai' | awk -v c=SN:{wildcards.chrom} '$1==\"@SQ\" && $2==c {{sub(\"LN:\",\"\",$3); print $3}}' > {output}"

rule slice:
    output: temp(f"{W}/reads/{{reads}}/{{chrom}}/{{sample}}.all")
    params: u = lambda wc: bam_url(wc.reads, wc.sample), conv = lambda wc: "fastq" if wc.reads == "hifi" else "fasta -n"
    threads: 2
    retries: 4                  # the remote seek fails now and then
    resources: net = 1
    log: f"{W}/logs/slice.{{reads}}.{{chrom}}.{{sample}}.log"
    shell: "samtools view -b -F 0x900 -@ {threads} '{params.u}##idx##{params.u}.bai' {wildcards.chrom} 2> {log} | samtools {params.conv} -@ {threads} - > {output} 2>> {log}"

rule ladder_hifi:
    input: reads = f"{W}/reads/hifi/{{chrom}}/{{sample}}.all", n = f"{W}/reads/hifi/{{chrom}}.len"
    output: [read_set("hifi", "{chrom}", "{sample}", c) for c in COVS]
    params: out = f"{W}/reads/hifi/{{chrom}}", covs = " ".join(COVS)
    log: f"{W}/logs/ladder.hifi.{{chrom}}.{{sample}}.log"
    shell: "python3 -I {S}/subset.py {input.reads} $(cat {input.n}) {params.out} {wildcards.sample} {params.covs} 2> {log}"

use rule ladder_hifi as ladder_sr with:
    input: reads = f"{W}/reads/sr/{{chrom}}/{{sample}}.all", n = f"{W}/reads/sr/{{chrom}}.len"
    output: [read_set("sr", "{chrom}", "{sample}", c) for c in COVS]
    params: out = f"{W}/reads/sr/{{chrom}}", covs = " ".join(COVS)
    log: f"{W}/logs/ladder.sr.{{chrom}}.{{sample}}.log"

# ---- the rank call: the sample of the panel the reads come from, with its margin D and share pi (gpeb --set-call mode) ---------
# run dir calls/<reads>/<chrom>/<panel>/<cfg>/rank.<sample>.G<cov>x; its log holds the set-call[mode] scores and the set-call[conf] line
CALLS = W + "/calls/{reads}/{chrom}/{panel}/{cfg}"

def call_mem(wc, input):    # measured <= 1.3 GB on chr9's default index (2.6x its .pex); the sweep's larger indexes need more
    return max(2, math.ceil(5 * os.path.getsize(input.idx) / 2**30) + 2)

rule rank_call:     # one thread: the call reads on one thread whatever --threads says
    input: gpeb = rules.gpeb.output.gpeb, idx = f"{W}/idx/{{chrom}}.{{panel}}.{{cfg}}.pex",
           reads = lambda wc: read_set(wc.reads, wc.chrom, wc.sample, wc.cov)
    output: f"{CALLS}/rank.{{sample}}.G{{cov}}x/log"
    resources: mem_gb = call_mem
    params: b = bench("rank_call", "{reads}.{panel}.{cfg}.{sample}.{cov}")
    shell:
        "{TIME} {params.b} {input.gpeb} --threads 1 --use-index {input.idx} --top 0 --classify {IX[k]} {input.reads} --set-call mode "
        "> /dev/null 2> {output}"

def call_dirs(units):
    return " ".join(os.path.dirname(f) for f in units)

# ---- bench.tsv and sizes.tsv: every timed command and every measured file of the run, rewritten when a target completes ---------
BENCH = f"python3 -I {S}/bench.py {W} {W}/results/bench.tsv {W}/results/sizes.tsv"
