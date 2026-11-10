![gpeb logo](logo.jpeg)

<!-- # gpeb -->

Graph pebbling for pangenome graphs: path inference and graph personalization.

From a pangenome GFA, `gpeb` finds the *pebbles* (the junctions that tell the
haplotype walks apart) and indexes their k-mers. With that index, **gpeb call** says
which sample a read set comes from, and **gpeb cover** builds a small personalized
graph for a sample from its reads.

## Install

**Prebuilt binary.** Each release carries `gpeb-linux-x86_64`, `gpeb-linux-aarch64`
(both static) and `gpeb-macos-arm64`. Download the one for your machine, then:

```bash
mv gpeb-linux-x86_64 gpeb && chmod +x gpeb
./gpeb --help
```

On macOS, a binary downloaded with a browser is quarantined; clear it once with
`xattr -d com.apple.quarantine gpeb`.

**From source.** Needs a C compiler, `make` and zlib.

```bash
git clone --recursive https://github.com/msgr0/gpeb gpeb && cd gpeb
make            # builds libs/gfatools/libgfa1.a, then ./gpeb
make test       # optional: the built-in self-test
```

If you cloned without `--recursive`, run `git submodule update --init libs/gfatools`
first. To build a binary that runs on other machines, use `make NATIVE=` (no
`-march=native`).

## gpeb call

The input is a GFA 1.1 with `W` lines (walks), such as an HPRC minigraph-cactus
chromosome graph.

**1. Pebble list.** One FASTA record per pebble: the junction edge, the walks it
belongs to and the `2 * span` bp window across it.

```bash
gpeb --span 60 graph.gfa > graph.peb
```

**2. Index, once per graph.** Writes a k-mer to haplotype-set index.

```bash
gpeb --index 31 --syncmer 16 --threads 16 graph-index.pex graph.gfa
```

**3. Classify reads** (FASTA or FASTQ, plain or gzipped). Per-read results go to
stdout as a TSV, the call for the whole read set to stderr. Passing the GFA as well
makes `gpeb` refuse an index built from a different graph.

```bash
gpeb --use-index graph-index.pex --classify 31 reads.fq.gz graph.gfa > reads.tsv
```

To call the sample behind the reads, take each read's most-voted haplotypes and rank
the samples:

```bash
gpeb --use-index graph-index.pex --top 0 --classify 31 reads.fq.gz --set-call mode > /dev/null
```

Its `set-call[conf]` line qualifies the call: `D`, the call's margin over the hardest
competing sample on the reads that tell the two apart (accept the call when D >= 1/4),
and `pi`, each sample's share of the reads.

### gpeb cover (personalized graph)

For a sample that is not in the panel, `--set-call cover` picks, per 200 kb window, the
few haplotypes that carry most of the reads' index hits, and writes the graph made of
them plus the references. The index must be built with positions on a reference colour
(`--mosaic-ref`):

```bash
gpeb --index 31 --syncmer 16 --threads 16 --mosaic-ref GRCh38 graph-index.pex graph.gfa

gpeb --threads 16 --use-index graph-index.pex graph.gfa --classify 31 reads.fq.gz \
     --set-call cover --mosaic-keep CHM13,GRCh38 \
     --mosaic-out sample --mosaic-gfa sample.gfa
```

This writes `sample.path.tsv` (the painting), `sample.sets.tsv` (the picks of every
window) and `sample.gfa`, a graph a short-read or HiFi aligner such as `vg giraffe`
can map to. Paired short reads: add `--mate R2.fq.gz`.

## Other commands

Not covered above; each is described by `./gpeb --help`.

- Cut modes: `--minimal`, `--acyclic`, `--early-cut`, `--full-graph`, `-c/--compress`,
  `--from-peb`
- Pebble filters: `--min-gap`, `--min-complexity`, `--colors`, `--color-sep`
- Index variants: `--keep-ambiguous`, `--keep-many`, `--fast-mask`, `--tol`,
  `--index-all`, `--subsample`, `--max-colours`, `--walk-dedup`
- Whole-set calls: `--set-call all|fraction|intersect|vote|sample|mode|mosaic|cover`, `--top`
- Cover tuning: `--cover-window`, `--cover-max`, `--cover-min-frac`, `--cover-min-gain`,
  `--cover-stray`, `--cover-rule`, `--cover-weight`, `--mosaic-check`
- Mosaic (two recombining haplotype tracks): `--set-call mosaic`, `--mosaic-block`,
  `--mosaic-reads`, `--mosaic-temper`, `--mosaic-exclude`
- Splicing a saved painting without reads: `--mosaic-splice`, `--make-gfai`, `--no-gfai`
- Splitting a large panel into one index per group: `--split-colors`, `--split-max-colors`,
  `--split-walk`, `--split-shared`, `--split-exclude`, `--split-dry-run`
- Diagnostics: `--kmer-stats`, `--pebble-ambiguity`, `--panel-jaccard`

## Reproducing the paper

[`workflows/recomb27`](workflows/recomb27) is the Snakemake workflow behind every result of the
RECOMB 2027 paper, from the HPRC release 2.1 graphs and GIAB reads; see its `README.md`.

## License

MIT, see `LICENSE`.
