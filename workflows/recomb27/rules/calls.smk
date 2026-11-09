# Calling the donor (paper: Results "Calling the donor with gpeb call", Figure 2; supplement "Calling the donor", per chromosome).
# Every HiFi and short-read set of every chromosome, on the full graph (which holds HG002, HG005 and HG002 as the parents' child) and
# on noHG002 (which holds HG005 only), called on the default index.
#   results/exp2/<chrom>.tsv   one row per call (call_aggregate.py): the top sample, the expected donor's share F, the runner-up,
#                              the margin D over the hardest competitor and the share pi of the reads that the call explains

rule calls:
    input: expand(f"{W}/results/exp2/{{chrom}}.tsv", chrom=CH)
    output: touch(f"{W}/results/exp2.done")
    shell: "{BENCH}"

rule exp2_chrom:
    input: lambda wc: [f"{W}/calls/{r}/{wc.chrom}/{p}/{D}/rank.{s}.G{v}x/log" for r in ("hifi", "sr") for p in PANELS for s in SM for v in COVS]
    output: f"{W}/results/exp2/{{chrom}}.tsv"
    params: dirs = lambda wc, input: call_dirs(input)
    shell: "python3 -I {S}/call_aggregate.py {output} {EXPECT} {params.dirs}"
