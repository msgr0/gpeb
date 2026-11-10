# RECOMB 2027 workflow

Every result of the paper *Colored Pebbling for Topology-Aware Indexing and Personalization of Pangenome Graphs*, from the HPRC
release 2.1 graphs and GIAB reads to the tables in `<workdir>/results/`.

| target | paper | results |
|---|---|---|
| `index_study` | Results "The pebble index"; supplement "Colors of the records" (Fig. 5), "Index parameters sweep" (Table 2) | `exp1/<chrom>.{index,calls,rare}.tsv` |
| `calls` | Results "Calling the donor" (Fig. 2); supplement "Calling the donor" (Fig. 6) | `exp2/<chrom>.tsv` |
| `whole_genome` | Results "Whole genome" (Fig. 3); supplement "Calling the donor on the whole genome" (Table 3) | `wg/{hifi,sr}.tsv`, `wg/{hifi,sr}.chroms.tsv`, `wg/index.tsv` |
| `personalize` | Results "Personalized graphs" (Fig. 4, Table 1); supplement "Personalized graphs" (read classes, bases) | `exp3/<chrom>.{tsv,paired.tsv,extra.tsv,paired_extra.tsv}` |
| every target | Table 1 and the costs in the text | `bench.tsv` (time and peak memory of every command), `sizes.tsv` (every index and graph) |

## Setup

1. The environment: `conda env create -f environment.yaml && conda activate gpeb-recomb27` (Snakemake, samtools, kmc, Python with
   numpy, a C compiler, make, zlib, curl, GNU time). git must be on PATH: gpeb is built from this repository at `gpeb_commit`.
2. The graphs: the HPRC release 2.1 Minigraph-Cactus GRCh38 graph split per chromosome, one `<chrom>.gfa` per chromosome in
   `gfa_dir`.
3. The paths, in `config/config.yaml` or in an overlay of your own passed after it (`config/pasteur.yaml` is the paper's):
   `workdir` (every intermediate file and `results/`), `gfa_dir`, `vg` (empty: vg 1.78.0 is downloaded), `kmc` (empty: on PATH).
4. Network access to NCBI's GIAB FTP: the reads are streamed from it, never downloaded whole.

## Running

From this directory:

```bash
snakemake <target> --profile profiles/local --configfile config/config.yaml [config/<overlay>.yaml]
```

Add `-n` to list the jobs without running them. The target comes right after `snakemake`: `--configfile` takes several files and
would swallow it. Without a target, `all` runs the four below. `profiles/local` runs on one machine: set `cores` and the memory cap
(`resources: mem_gb`) to yours; each job's memory comes from the `mem_gb_*` keys of the config.

### `index_study`

On the full graph of each chromosome in `chroms`: the default index (`index`: k 31, min-gap 60, syncmer 16, mosaic reference
GRCh38) and a sweep around it, one parameter at a time (`exp1`: min-gap and syncmer rebuilds, indexes cut to the k-mers of at most
`max_colours` colors). Every index is scored by the rank calls of HG002 and HG005 and by the cover of the parents' reads at
`exp1.call_covs`; `--cover-min-frac` is swept on the default index (`exp1.min_fracs`).

```bash
snakemake index_study --profile profiles/local --configfile config/config.yaml
```

### `calls`

The donor call (`gpeb --top 0 --classify 31 <reads> --set-call mode`) of every read set on the default index of the full graph
and of the graph without HG002 (`noHG002`). The read sets: each sample of `samples`, HiFi (`hifi_bam`) and short reads (`sr_bam`),
one chromosome sliced from GIAB's aligned BAMs with samtools and cut into the coverages of `covs`. `expect` names the sample each
set should call.

```bash
snakemake calls --profile profiles/local --configfile config/config.yaml
```

### `whole_genome`

Every read of GIAB's raw files (`wg_ubam`: PacBio HiFi uBAMs; `wg_fastq`: Element FASTQs), cut into the whole-genome coverages of
`wg_covs`, classified against the default index of each chromosome of `wg_chroms`, kept on the chromosome where it has the most
hits and called once over the genome, on the full graphs and on noHG002; again on the reads with at least `wg_min_hits` hits.
Needs a GFA for each of the 23 chromosomes of `wg_chroms`. The reads stream from `giab_base`; `wg_base` points to a local copy
instead.

```bash
snakemake whole_genome --profile profiles/local --configfile config/config.yaml
```

### `personalize`

On noHG002, from HG002's HiFi reads at each coverage of `covs`: the personalized graph of `gpeb --set-call cover` with at most N
colors per window (`cover_caps`), the graphs of `vg haplotypes` from KMC k-mers (`vg_cfgs`: `d32` the diploid sample, `sN` N
haplotypes; `kmc_ci`), N random haplotypes (`random_ns`, `random_seeds`) and the whole noHG002 graph. The reads are mapped on
every graph with `vg giraffe -b hifi` and compared read by read.

```bash
snakemake personalize --profile profiles/local --configfile config/config.yaml
```

The distance index of chromosome 22's cover with 32 colors asks for 512 GB (`mem_gb_gdist`): lower `cover_caps`, drop chr22
from `chroms`, or run that job on a machine with the memory.

### The numbers of the paper

```bash
python3 scripts/tex_numbers.py <workdir>/results
```

prints, from the tables of the four targets, the numbers of the text and the tables.

## Quick check

`config/smoke.yaml`: chr21, HG002 and HG003, 0.5x and 1x, one sweep setting, N = 4.

```bash
snakemake calls index_study personalize --profile profiles/local --configfile config/config.yaml config/smoke.yaml
```

Self-tests of the scripts: `python3 scripts/<name>.py --selftest` for `subset.py`, `call_aggregate.py`, `wg_sets.py`, and
`python3 scripts/conf.py`.

## On Slurm

With `snakemake-executor-plugin-slurm` installed, every job becomes its own Slurm job:

```bash
snakemake <target> --profile profiles/slurm --configfile config/config.yaml config/<overlay>.yaml
```

Set `slurm_account` (and `slurm_partition`, if your cluster needs one) and the run times in `profiles/slurm/config.yaml` to your
cluster's. Run the controller in a login shell (`bash -lc`): the plugin needs `sacctmgr` on PATH.
