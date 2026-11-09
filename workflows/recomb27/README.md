# RECOMB 2027 workflow

Every result of the paper *Colored Pebbling for Topology-Aware Indexing and Personalization of Pangenome Graphs*, from the HPRC
release 2.1 graphs and GIAB reads to the tables in `<workdir>/results/`; `scripts/tex_numbers.py` prints from them the numbers of
the text and tables.

| target | paper | results |
|---|---|---|
| `index_study` | Results "The pebble index"; supplement "Colors of the records" (Fig. 5), "Index parameters sweep" (Table 2) | `exp1/<chrom>.{index,calls,rare}.tsv` |
| `calls` | Results "Calling the donor" (Fig. 2); supplement "Calling the donor" (Fig. 6) | `exp2/<chrom>.tsv` |
| `whole_genome` | Results "Whole genome" (Fig. 3); supplement "Calling the donor on the whole genome" (Table 3) | `wg/{hifi,sr}.tsv`, `wg/{hifi,sr}.chroms.tsv`, `wg/index.tsv` |
| `personalize` | Results "Personalized graphs" (Fig. 4, Table 1); supplement "Personalized graphs" (read classes, bases) | `exp3/<chrom>.{tsv,paired.tsv,extra.tsv,paired_extra.tsv}` |
| every target | Table 1 and the costs in the text | `bench.tsv` (time and peak memory of every command), `sizes.tsv` (every index and graph) |

`all` (the default) runs the four.

## What it does

- **gpeb** and `pex_rare` are built from this repository at `gpeb_commit` (default `HEAD`; `bin/COMMIT` records it and the binary's sha256).
- **Graphs**: per chromosome, the released GFA (`full`, 464 colors) and the graph without HG002's two haplotypes (`noHG002`,
  `scripts/prune_gfa.py`: drops the walks, then the segments and links no remaining walk uses).
- **Index**: `gpeb --index 31 --min-gap 60 --syncmer 16 --mosaic-ref GRCh38` (the default), plus the sweep of `exp1` in the config.
- **Per-chromosome reads**: GIAB HiFi (CHM13-aligned) and Element short reads (GRCh38-aligned), one chromosome sliced from the remote
  BAM with samtools, cut into a nested ladder 0.5-30x by a hash of the read name (`scripts/subset.py`).
- **Calls**: `gpeb --top 0 --classify 31 <reads> --set-call mode`; its log holds the call, every sample's share F, and the margin D and
  share pi of the call (`set-call[conf]`); `scripts/call_aggregate.py` makes one row per call.
- **Whole genome**: every read of the raw files (PacBio uBAMs, Element FASTQs, streamed from NCBI) cut into whole-genome rungs by the
  same hash (`scripts/wg_reads.py`), called against each of the 23 chromosomes' indexes, kept on its best chromosome and called once
  over all (`scripts/wg_sets.py`, `scripts/conf.py`), on every read and on the reads with >= 5 hits (`wg_min_hits`).
- **Personalized graphs**: on noHG002, from HG002's HiFi reads: `gpeb --set-call cover --cover-max N` spliced into a GFA, against
  `vg haplotypes` (KMC k-mers, N haplotypes or the diploid sample), N random haplotypes and the whole graph; every graph mapped with
  `vg giraffe -b hifi`, edits per read by `scripts/edits.py`, compared read by read by `scripts/exp3_aggregate.py` and `exp3_extra.py`.

## Requirements

- `environment.yaml` (`conda env create -f environment.yaml`): Snakemake, samtools, kmc, Python with numpy, a C compiler,
  make, zlib, curl and GNU time (`time_bin`); plus git (the build archives this repository) and network access to NCBI's GIAB FTP.
- vg 1.78.0: `vg` in the config points to the binary, or leave it empty to download the release.
- The HPRC release 2.1 Minigraph-Cactus GRCh38 graph split per chromosome, one `<chrom>.gfa` per chromosome in `gfa_dir`
  (chromosomes 1-22 and X).
- With Slurm: `snakemake-executor-plugin-slurm`.

## Running

From this directory, with `gfa_dir` and `workdir` set in `config/config.yaml` (or in an overlay of your own):

```bash
snakemake --profile profiles/local -n                     # what would run
snakemake --profile profiles/local                        # everything, on one machine (cores and memory in profiles/local)
snakemake calls --profile profiles/local                  # one target
```

The paper ran on pasteur (64 cores, 251 GB) and on the Alliance clusters, with the paths of `config/pasteur.yaml`,
`config/fir.yaml` and `config/narval.yaml`:

```bash
snakemake --profile profiles/local --configfile config/config.yaml config/pasteur.yaml
snakemake whole_genome --profile profiles/slurm --configfile config/config.yaml config/fir.yaml   # on Fir (run in `bash -lc`)
```

Targets come right after `snakemake`: `--configfile` and `--resources` take several values and would swallow them. Rules declare
`mem_gb` from measured peaks; `profiles/local` caps their sum, `rules/slurm.smk` turns them into Slurm requests.
Then:

```bash
python3 scripts/tex_numbers.py <workdir>/results
```

**A quick check** (`config/smoke.yaml`: chr21, HG002 and HG003, 0.5x and 1x, one sweep setting, N = 4), about an hour:

```bash
snakemake calls index_study personalize --profile profiles/local --configfile config/config.yaml config/smoke.yaml
```

Self-tests of the scripts: `python3 scripts/<name>.py --selftest` for `subset.py`, `call_aggregate.py`, `wg_sets.py`, and
`python3 scripts/conf.py`.

## Cost

Measured on pasteur (two 16-core EPYC 7301, 251 GB) and the Alliance clusters:

- The 23 indexes of the full graphs: 7 h in all, at most 31 min and 59 GiB each; pruning HG002 out of a graph: up to 75 min, 2.4 GiB.
- `vg haplotypes`' panel indexes of the ten chromosomes: 21 h, up to 170 GiB (chromosome 2's distance index).
- A HiFi call: from 12 s at 0.5x to 4.5 min at 30x, under 2.5 GiB; whole-genome calls: about 40 CPU hours per read type and graph.
- The distance index of chromosome 22's cover with 32 colors takes 3-4.5 h and up to 450 GiB at 5-30x: run it on a large-memory node
  (`mem_gb_gdist`); every other job fits in 230 GB.
- Disk: the whole-genome first pass keeps a fifth of every read (about 200 GB for the two read types); per-chromosome calls keep only
  their logs.
