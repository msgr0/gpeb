# Whole genome from unaligned reads (paper: Results "Whole genome", Figure 3; supplement "Calling the donor on the whole genome from
# unaligned reads", Table 3). The read sets of calls.smk are the reads an aligner put on each chromosome; here no aligner picks them:
# every read of GIAB's raw files (HiFi: the PacBio uBAMs; short reads: the Element FASTQs, each mate a read named NAME/1, NAME/2) is
# cut into whole-genome rungs by the same name hash, classified against the default index of each of chromosomes 1-22 and X, kept on
# the chromosome where it has the most hits, and called once over the pooled reads, on the full graphs and on noHG002.
#   results/wg/<reads>.tsv          conf.py rows: chrom = chrN, each chromosome's index alone on all the reads (gpeb's own call);
#                                   chrom = wg, the genome-wide call; wg_h<T>, the same on the reads with >= T hits on their chromosome
#   results/wg/<reads>.chroms.tsv   reads kept per chromosome, per call
#   results/wg/index.tsv            the 23 default indexes of the full graphs (Table 3)
WG, WC = config["wg_chroms"], config["wg_covs"]
R = W + "/wg"
WCALL = R + "/calls/{reads}/{chrom}/{panel}/" + D + "/rank.{sample}.G{cov}x"
WMERGE = R + "/calls/{reads}/wg/{panel}/" + D + "/rank.{sample}.G{cov}x"

def raw_files(reads, sample):       # stem -> URL (or path) of one raw file
    b = config["wg_base"] or config["giab_base"]
    if reads == "hifi":
        return {os.path.basename(p)[:-len(".bam")]: f"{b}/{p}" for p in config["wg_ubam"][sample]}
    return {f"R{m}": f"{b}/{config['wg_fastq'][sample]}_R{m}.fastq.gz" for m in (1, 2)}

def stream(wc):                     # the raw file as one-line FASTA on stdout
    u = raw_files(wc.reads, wc.sample)[wc.stem]
    if wc.reads == "hifi":
        return f"samtools fasta -@ 2 '{u}'"
    return f"curl -fsSL '{u}' | gzip -dc | awk -v m={wc.stem[1]} 'NR % 4 == 1 {{print \">\" substr($1, 2) \"/\" m; next}} NR % 4 == 2'"

def mem_by_cov(gb):                 # gb = per-rung list for WC
    return lambda wc: gb[WC.index(wc.cov)]

rule whole_genome:
    input: expand(f"{W}/results/wg/{{reads}}.{{t}}", reads=["hifi", "sr"], t=["tsv", "chroms.tsv"]), f"{W}/results/wg/index.tsv"
    output: touch(f"{W}/results/wg.done")
    shell: "{BENCH}"

rule wg_genome_len:     # CHM13v2.0: the summed @SQ lengths of a CHM13-aligned HiFi BAM (depth is bases / this, as in calls.smk)
    output: f"{R}/genome.len"
    params: u = bam_url("hifi", SM[0])
    retries: 4
    resources: net = 1
    shell: "samtools view -H '{params.u}' | awk '$1==\"@SQ\" {{for (i = 2; i <= NF; i++) if ($i ~ /^LN:/) s += substr($i, 4)}} END {{print s}}' > {output}"

rule wg_sup:            # one raw file streamed from NCBI, no copy on disk: its reads with u < wg_umax, and the bases of all its reads
    output: fa = f"{R}/reads/{{reads}}/{{sample}}/{{stem}}.sup.fa.gz", n = f"{R}/reads/{{reads}}/{{sample}}/{{stem}}.bases"
    params: cmd = stream
    threads: 3
    retries: 3
    resources: net = 1, mem_gb = 6
    log: f"{W}/logs/wg_sup.{{reads}}.{{sample}}.{{stem}}.log"
    shell: "{params.cmd} 2> {log} | python3 -I {S}/wg_reads.py sup {config[wg_umax]} {output.fa} {output.n} 2>> {log}"

rule wg_rungs:          # depth = all the sample's bases / genome length; nested whole-genome rungs cut from the first pass
    input: fa = lambda wc: [f"{R}/reads/{wc.reads}/{wc.sample}/{s}.sup.fa.gz" for s in raw_files(wc.reads, wc.sample)],
           n = lambda wc: [f"{R}/reads/{wc.reads}/{wc.sample}/{s}.bases" for s in raw_files(wc.reads, wc.sample)], glen = rules.wg_genome_len.output
    output: [f"{R}/reads/{{reads}}/{{sample}}.G{c}x.fa.gz" for c in WC]
    params: pfx = f"{R}/reads/{{reads}}/{{sample}}", covs = ",".join(WC), sups = lambda wc, input: " ".join(f"{a}:{b}" for a, b in zip(input.fa, input.n))
    resources: mem_gb = 6
    log: f"{W}/logs/wg_rungs.{{reads}}.{{sample}}.log"
    shell: "python3 -I {S}/wg_reads.py rungs {config[wg_umax]} {input.glen} {params.pfx} {params.covs} {params.sups} 2> {log}"

rule wg_call:           # one chromosome's index on every read of a rung; the per-read table is reduced as it streams (wg_sets.py)
    input: gpeb = rules.gpeb.output.gpeb, idx = f"{W}/idx/{{chrom}}.{{panel}}.{D}.pex", fa = f"{R}/reads/{{reads}}/{{sample}}.G{{cov}}x.fa.gz"
    output: r = f"{WCALL}/reads.tsv.gz", s = f"{WCALL}/sets.tsv.gz", log = f"{WCALL}/log"
    resources: mem_gb = mem_by_cov([4, 5, 6, 8])    # measured: <= 3.2 GB at 0.5x (short reads)
    params: d = WCALL, b = bench("wg_call", "{reads}.{panel}.{sample}.{cov}")
    shell:
        "{TIME} {params.b} {input.gpeb} --threads 1 --use-index {input.idx} --top 0 --classify {IX[k]} {input.fa} --set-call mode "
        "2> {output.log} | python3 -I {S}/wg_sets.py reduce {params.d} 2> {params.d}/reduce.log"

rule wg_merge:          # each read on its best chromosome; the genome-wide call at every wg_min_hits threshold
    input: lambda wc: [f"{R}/calls/{wc.reads}/{c}/{wc.panel}/{D}/rank.{wc.sample}.G{wc.cov}x/reads.tsv.gz" for c in WG]
    output: c = f"{WMERGE}/conf.tsv", n = f"{WMERGE}/chroms.tsv"
    resources: mem_gb = mem_by_cov([12, 18, 30, 60])  # measured on short reads: 7.9 GB at 0.5x; sets grow ~1.5x per doubling
    params: d = WMERGE, dirs = lambda wc, input: " ".join(os.path.dirname(f) for f in input), t = ",".join(map(str, config["wg_min_hits"])),
            b = f"{W}/bench/wg_merge__wg__{{reads}}.{{panel}}.{{sample}}.{{cov}}.time"
    shell: "{TIME} {params.b} python3 -I {S}/wg_sets.py merge {params.d} --min-hits {params.t} {params.dirs}"

rule wg_table:
    input: calls = lambda wc: expand(f"{WCALL}/log", reads=wc.reads, chrom=WG, panel=PANELS, sample=SM, cov=WC),
           merged = lambda wc: expand(f"{WMERGE}/conf.tsv", reads=wc.reads, panel=PANELS, sample=SM, cov=WC)
    output: t = f"{W}/results/wg/{{reads}}.tsv", c = f"{W}/results/wg/{{reads}}.chroms.tsv"
    params: dirs = lambda wc, input: " ".join(os.path.dirname(f) for f in input.calls)
    shell:
        """
        NPROC=1 python3 -I {S}/wg_sets.py logs {output.t} {params.dirs}
        cat {input.merged} >> {output.t}
        printf 'panel\\tsample\\tcov\\tchrom\\treads\\n' > {output.c}
        for f in {input.merged}; do d=$(dirname $f); p=$(echo $d | rev | cut -d/ -f3 | rev); s=$(basename $d | cut -d. -f2)
            v=$(basename $d | sed 's/.*\\.G\\(.*\\)x$/\\1/'); tail -n +2 $d/chroms.tsv | sed "s/^/$p\\t$s\\t$v\\t/" >> {output.c}; done
        """

rule wg_index:          # Table 3: the default index of the full graph of each chromosome, with the size of its GFA
    input: idx = expand(f"{W}/idx/{{chrom}}.full.{D}.pex", chrom=WG), gfa = expand(f"{W}/panel/{{chrom}}.full.gfa", chrom=WG)
    output: f"{W}/results/wg/index.tsv"
    params: tmp = f"{W}/tmp/wg_index"
    shell:
        """
        mkdir -p {params.tmp}
        for c in {WG}; do python3 -I {S}/index_stats.py {W} $c full {params.tmp}/$c.tsv {D} $(readlink -f {W}/panel/$c.full.gfa); done
        head -1 {params.tmp}/{WG[0]}.tsv > {output}; for c in {WG}; do tail -n +2 {params.tmp}/$c.tsv >> {output}; done
        """
