# RECOMB 2027 workflow

This workflow reproduces every result of the paper *Colored Pebbling for Topology-Aware Indexing and Personalization of Pangenome
Graphs*. It starts from the HPRC release 2.1 graphs and the GIAB reads and ends with the tables in `<workdir>/results/`.

The paper uses two gpeb commands. **gpeb call** names the sample behind a set of reads (`gpeb --set-call mode`); **gpeb cover**
builds a personalized graph from them (`gpeb --set-call cover`).

| target | paper | results |
|---|---|---|
| `index_study` | Results "The pebble index"; supplement "Colors of the records" (Fig. 5), "Index parameters sweep" (Table 2) | `exp1/<chrom>.{index,calls,rare}.tsv` |
| `calls` | Results "Calling the donor" (Fig. 2); supplement "Calling the donor" (Fig. 6) | `exp2/<chrom>.tsv` |
| `whole_genome` | Results "Whole genome" (Fig. 3); supplement "Calling the donor on the whole genome" (Table 3) | `wg/{hifi,sr}.tsv`, `wg/{hifi,sr}.chroms.tsv`, `wg/index.tsv` |
| `personalize` | Results "Personalized graphs" (Fig. 4, Table 1); supplement "Personalized graphs" | `exp3/<chrom>.{tsv,paired.tsv,extra.tsv,paired_extra.tsv}` |
| every target | Table 1 and the costs in the text | `bench.tsv`, `sizes.tsv` |

## Before you start

Create the environment with `conda env create -f environment.yaml` and activate it (`gpeb-recomb27`). You also need git, because
the workflow builds gpeb from this repository, and network access to NCBI, because the reads are streamed from GIAB's FTP.

Put the HPRC graph, split per chromosome, in a folder as `<chrom>.gfa`. Then tell the workflow where things are, in
`config/config.yaml` or in a small overlay file of your own (`config/pasteur.yaml` is our machine yaml): `gfa_dir` is the folder of
the graphs and `workdir` is where everything is written. vg is downloaded unless `vg` points to a binary.

## Running

Every run has the same shape, from this directory:

```bash
snakemake <target> --profile profiles/local --configfile config/config.yaml config/<overlay>.yaml
```

Put the target first: `--configfile` takes several files and would swallow it. Add `-n` to see what would run. With no target, the
workflow runs all four. Set `cores` and `mem_gb` in `profiles/local/config.yaml` to your machine.

**`index_study`** builds the pebble index of each chromosome's full graph, then rebuilds it changing one parameter at a time. Every
index is judged by gpeb call on HG002's and HG005's reads, and by gpeb cover on the parents' reads.

**`calls`** runs gpeb call on every read set: HiFi and short reads of each sample, one chromosome at a time, at coverages from 0.5x
to 30x. It calls each set twice, on the full graph and on the graph without HG002, where the right answer is no longer there.

**`whole_genome`** does the same without an aligner. It takes every read of GIAB's raw files, finds the chromosome each read fits
best, and runs one gpeb call over the whole genome. It needs the graphs of all 23 chromosomes.

**`personalize`** builds personalized graphs for HG002 from the graph without HG002: gpeb cover, vg haplotypes, random haplotypes,
and the whole graph as a baseline. It maps HG002's reads on each with vg giraffe and compares them read by read. One job, chr22's
cover with 32 colors, needs 512 GB of memory; drop chr22 from `chroms` if your machine has less.

When the targets are done, `python3 scripts/tex_numbers.py <workdir>/results` prints the numbers of the paper's text and tables.

## A quick check

`config/smoke.yaml` shrinks the run to chromosome 21, two samples and two coverages:

```bash
snakemake calls index_study personalize --profile profiles/local --configfile config/config.yaml config/smoke.yaml
```

## On a Slurm cluster

Use `--profile profiles/slurm` instead (it needs `snakemake-executor-plugin-slurm`), after writing your account, and a partition
if your cluster needs one, in `profiles/slurm/config.yaml`. Each job then runs as its own Slurm job.
