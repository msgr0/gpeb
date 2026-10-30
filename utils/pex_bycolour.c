// pex_bycolour -- a .pex index turned around: colour -> the rare k-mers it touches.
//
// A .pex record is k-mer -> colour mask. This writes the inverse as a FASTA-like
// file, one record per colour (index order): the header is the colour name and the
// single sequence line is the comma-separated list of the k-mers whose mask has that
// colour, rarest first (fewest colours touching the k-mer), ties in key order.
//
//   >HG002#1
//   ACGT...,GGCA...,...
//
// Only k-mers touched by at most MAX colours are listed (default 10): the full
// inverse holds mean-popcount * records entries (hundreds of times the index on a
// 462-colour panel). `--max 1` lists each colour's private k-mers. A colour with none
// still gets its header and an empty line. k-mers are canonical, A=0 C=1 G=2 T=3.
//
//   cc -O2 -o pex_bycolour pex_bycolour.c
//   ./pex_bycolour [--max 10] <file.pex> > by_colour.fa      # PEX2..PEX5
//
// All kept records are held in memory (8 + 8*ceil(nc/64) bytes each). See pex_rare.c
// for the same records as k-mer -> colours rows.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

static uint32_t NW;

static uint32_t popc(const uint64_t *r) {            // record = key, mask[NW]
  uint32_t p = 0; for (uint32_t w = 0; w < NW; w++) p += __builtin_popcountll(r[1 + w]);
  return p;
}
static int cmp(const void *a, const void *b) {
  const uint64_t *x = a, *y = b;
  uint32_t px = popc(x), py = popc(y);
  if (px != py) return px < py ? -1 : 1;
  return x[0] < y[0] ? -1 : x[0] > y[0];
}

#define RD(p, n) do { if (fread(p, 1, n, f) != (size_t)(n)) { fprintf(stderr, "pex_bycolour: truncated header\n"); return 1; } } while (0)

int main(int argc, char **argv) {
  uint32_t max = 10;
  int a = 1;
  if (a + 1 < argc && !strcmp(argv[a], "--max")) { max = (uint32_t)atoi(argv[a + 1]); a += 2; }
  if (a >= argc) { fprintf(stderr, "usage: pex_bycolour [--max N] file.pex > by_colour.fa\n"); return 1; }
  FILE *f = fopen(argv[a], "rb"); if (!f) { perror("open"); return 1; }
  char m[4]; RD(m, 4);
  int v5 = !memcmp(m, "PEX5", 4), v4 = v5 || !memcmp(m, "PEX4", 4), v3 = v4 || !memcmp(m, "PEX3", 4);
  if (!v3 && memcmp(m, "PEX2", 4)) { fprintf(stderr, "pex_bycolour: bad magic\n"); return 1; }
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
    if (fread(buf, 8, rw, f) != rw) { fprintf(stderr, "pex_bycolour: short read at %llu\n", (unsigned long long)i); return 1; }
    if (v5) fseek(f, 2, SEEK_CUR);                      // uint16 block id
    uint32_t pc = popc(buf);
    if (!pc || pc > max) continue;
    if (nk == cap) { cap *= 2; rec = realloc(rec, cap * rw * 8); if (!rec) { perror("realloc"); return 1; } }
    memcpy(rec + nk++ * rw, buf, rw * 8);
  }
  fclose(f);
  qsort(rec, nk, rw * 8, cmp);

  setvbuf(stdout, NULL, _IOFBF, 1 << 20);
  char *s = malloc(k + 1); s[k] = 0;
  uint64_t entries = 0;
  for (uint32_t c = 0; c < nc; c++) {
    printf(">%s\n", name[c]);
    for (size_t i = 0, first = 1; i < nk; i++) {
      const uint64_t *r = rec + i * rw;
      if (!(r[1 + c / 64] >> (c % 64) & 1)) continue;
      for (uint32_t j = 0; j < k; j++) s[j] = "ACGT"[r[0] >> 2 * (k - 1 - j) & 3];
      printf(first ? "%s" : ",%s", s); first = 0; entries++;
    }
    putchar('\n');
  }
  fprintf(stderr, "%zu k-mers with 1..%u colours, %llu colour entries over %u colours\n",
          nk, max, (unsigned long long)entries, nc);
  return 0;
}
