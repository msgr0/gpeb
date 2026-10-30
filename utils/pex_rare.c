// pex_rare -- the rare k-mers of a .pex index, ranked by how few colours touch them.
//
// A .pex record is already k-mer -> colour mask (the key is the canonical 2-bit
// k-mer, A=0 C=1 G=2 T=3, first base in the high bits), so no inversion of a
// colour -> k-mers table is needed: the popcount of the mask IS the number of
// colours the k-mer touches. This keeps the records with 1..MAX colours and prints
// them sorted by that count (1 = private to one colour), ties grouped by colour set:
//
//   cc -O2 -o pex_rare pex_rare.c
//   ./pex_rare [--max 10] <file.pex> > rare.tsv      # PEX2..PEX5
//
// stdout: ncolours <TAB> kmer <TAB> comma-separated colour names
// stderr: records kept per ncolours, and the private (ncolours = 1) count per colour.
//
// Rarity is relative to what the index holds: junction candidates only (single-locus
// unless built with --keep-many, thinned by --syncmer), and only the colours of THIS
// index, so with --split-colors a k-mer private in one group may be shared in another.
// All kept records are held in memory (8 + 8*ceil(nc/64) bytes each).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static uint32_t NW;

static int cmp(const void *a, const void *b) {
  const uint64_t *x = a, *y = b;                // record = key, mask[NW]
  uint32_t px = 0, py = 0;
  for (uint32_t w = 0; w < NW; w++) { px += __builtin_popcountll(x[1 + w]); py += __builtin_popcountll(y[1 + w]); }
  if (px != py) return px < py ? -1 : 1;
  for (uint32_t w = 0; w < NW; w++) if (x[1 + w] != y[1 + w]) return x[1 + w] < y[1 + w] ? -1 : 1;
  return x[0] < y[0] ? -1 : x[0] > y[0];
}

#define RD(p, n) do { if (fread(p, 1, n, f) != (size_t)(n)) { fprintf(stderr, "pex_rare: truncated header\n"); return 1; } } while (0)

int main(int argc, char **argv) {
  uint32_t max = 10;
  int a = 1;
  if (a + 1 < argc && !strcmp(argv[a], "--max")) { max = (uint32_t)atoi(argv[a + 1]); a += 2; }
  if (a >= argc) { fprintf(stderr, "usage: pex_rare [--max N] file.pex > rare.tsv\n"); return 1; }
  FILE *f = fopen(argv[a], "rb"); if (!f) { perror("open"); return 1; }
  char m[4]; RD(m, 4);
  int v5 = !memcmp(m, "PEX5", 4), v4 = v5 || !memcmp(m, "PEX4", 4), v3 = v4 || !memcmp(m, "PEX3", 4);
  if (!v3 && memcmp(m, "PEX2", 4)) { fprintf(stderr, "pex_rare: bad magic\n"); return 1; }
  uint64_t u64; uint32_t sl, k, u32, sk = 0, nc;
  RD(&u64, 8); RD(&u64, 8); RD(&sl, 4); fseek(f, sl, SEEK_CUR);
  RD(&k, 4); RD(&u32, 4);                               // k, span
  if (v3) RD(&sk, 4);
  if (v4) RD(&u32, 4);                                  // flags
  if (v5) { RD(&u32, 4); RD(&u64, 8); RD(&sl, 4); fseek(f, sl, SEEK_CUR); }
  RD(&nc, 4);
  char **name = calloc(nc, sizeof *name);
  for (uint32_t c = 0; c < nc; c++) { RD(&sl, 4); name[c] = calloc(sl + 1, 1); RD(name[c], sl); }
  uint64_t n; RD(&n, 8);
  NW = (nc + 63) / 64; if (!NW) NW = 1;
  fprintf(stderr, "k=%u s=%u nc=%u n=%llu\n", k, sk, nc, (unsigned long long)n);

  size_t rw = 1 + NW, nk = 0, cap = 1 << 20;
  uint64_t *rec = malloc(cap * rw * 8), *buf = malloc(rw * 8);
  for (uint64_t i = 0; i < n; i++) {
    if (fread(buf, 8, rw, f) != rw) { fprintf(stderr, "pex_rare: short read at %llu\n", (unsigned long long)i); return 1; }
    if (v5) fseek(f, 2, SEEK_CUR);                      // uint16 block id
    uint32_t pc = 0; for (uint32_t w = 0; w < NW; w++) pc += __builtin_popcountll(buf[1 + w]);
    if (!pc || pc > max) continue;
    if (nk == cap) { cap *= 2; rec = realloc(rec, cap * rw * 8); if (!rec) { perror("realloc"); return 1; } }
    memcpy(rec + nk++ * rw, buf, rw * 8);
  }
  fclose(f);
  qsort(rec, nk, rw * 8, cmp);

  uint64_t *hist = calloc(max + 1, 8), *priv = calloc(nc, 8);
  char *s = malloc(k + 1); s[k] = 0;
  for (size_t i = 0; i < nk; i++) {
    const uint64_t *r = rec + i * rw;
    uint32_t pc = 0; for (uint32_t w = 0; w < NW; w++) pc += __builtin_popcountll(r[1 + w]);
    hist[pc]++;
    for (uint32_t j = 0; j < k; j++) s[j] = "ACGT"[r[0] >> 2 * (k - 1 - j) & 3];
    printf("%u\t%s\t", pc, s);
    for (uint32_t c = 0, first = 1; c < nc; c++)
      if (r[1 + c / 64] >> (c % 64) & 1) { printf(first ? "%s" : ",%s", name[c]); first = 0; if (pc == 1) priv[c]++; }
    putchar('\n');
  }
  fprintf(stderr, "ncolours\tkmers\n");
  for (uint32_t c = 1; c <= max; c++) fprintf(stderr, "%u\t%llu\n", c, (unsigned long long)hist[c]);
  fprintf(stderr, "private per colour:\n");
  for (uint32_t c = 0; c < nc; c++) fprintf(stderr, "%s\t%llu\n", name[c], (unsigned long long)priv[c]);
  return 0;
}
