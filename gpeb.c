//
//  gpeb.c -- minimal cycle-aware guardian-list pebbling for bidirected GFA.
//
//  Build (tool):  cc -O2 -Ilibs gpeb.c libs/gfatools/libgfa1.a -lz -lpthread -lm -o gpeb
//  Run:           ./gpeb in.gfa > out.fa
//  Self-check:    cc -DGPEB_TEST -O2 -Ilibs gpeb.c libs/gfatools/libgfa1.a -lz -lpthread -lm -o t && ./t
//
// Must precede every include: -std=c11 (strict ISO) hides POSIX clock_gettime/
// CLOCK_MONOTONIC from the glibc headers, and the build fails.
#define _POSIX_C_SOURCE 200809L
#ifdef __APPLE__
#define _DARWIN_C_SOURCE         // with _POSIX_C_SOURCE alone macOS hides the nanosecond stat fields
#endif

#include <gfatools/gfa.h>
#include <gfatools/gfa-priv.h>   // gfa_arc_a, gfa_n_vtx, gfa_comp_table
#include <errno.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>                // --classify reads a plain or gzipped FASTA/FASTQ
#include <pthread.h>             // --threads: parallel walk/segment streaming
#include <time.h>                // phase timing (--index build breakdown)
#include <sys/stat.h>            // the .gfai: the size and modification time of the GFA it describes
#include <math.h>                // --set-call mosaic: exp/log of the emission table

#ifdef __APPLE__
#define ST_MTIM(sb) ((sb).st_mtimespec)
#define ST_ATIM(sb) ((sb).st_atimespec)
#else
#define ST_MTIM(sb) ((sb).st_mtim)
#define ST_ATIM(sb) ((sb).st_atim)
#endif

#define GSPAN   30               // default --span window (bp)
#define UNREACH UINT32_MAX

// bg_t  = bidirected graph (orientation-preserving "st-graph")
// Node ids are ORIENTED vertices v = (segment_id << 1) | strand, strand 0='+',
// 1='-', with n == gfa_n_vtx(gfa) == 2*n_seg. sname[]/seq[] stay indexed by
// SEGMENT (v>>1); a '-' vertex is reverse-complemented at emit time.
typedef struct {
  uint32_t   n;      // n_nodes: oriented vertex count (2*n_seg)
  uint32_t **adj;    // adj[v]  = oriented out-neighbour ids
  uint32_t  *adjn;   // adjn[v] = length of adj[v]
  uint32_t **pv;     // pv[k]   = oriented vertices of walk k (paths)
  uint32_t  *pn;     // pn[k]   = length of pv[k] (path_n)
  uint32_t   np;     // n_paths: number of walks
  char     **sname;  // sname[seg] = segment name (borrowed from gfa)
  char     **seq;    // seq[seg]   = segment sequence (borrowed from gfa)
  uint32_t  *slen;   // slen[seg]  = its length, precomputed (nsl() runs once per k-mer hit)
  char     **pname;  // pname[k]   = color name, interned into cpool (borrowed)
  char     **cpool;  // cpool[c]   = the distinct color names, OWNED
  uint32_t   ncpool; // number of distinct colors
  uint32_t  *psrc;   // psrc[k]    = index of walk k in the GFA's walk/path table; differs from k under --colors
} bg_t;

// gpb_t = guardian pebble: one ranked placement on a directed edge head->tail.
// Per-placement (not deduped): a cyclic walk crossing an edge twice yields two.
// Color is the walk index `pid` (one color per walk); rank is per-walk.
typedef struct {
  uint32_t head, tail; // oriented edge endpoints
  uint32_t rank;       // 1-based rank among this walk's pebbles, in walk order
  uint32_t pid;        // walk that placed it (color + span lookup)
  uint32_t pos;        // index of head in pv[pid] (exact edge slot, cycle-safe)
} gpb_t;

// gset_t = guardian set (dynamic array of gpb_t)
typedef struct { gpb_t *p; uint32_t n, cap; } gset_t;

// bg_build(walks, full_graph, colors):
//   walks: 1 = W-lines, 0 = P-lines.
//   full_graph: 0 builds reachability from the union of walk edges only (G_P), so a
//     merge no sample crosses plants no guardian; 1 uses all L-lines (more guardians).
//   colors (--colors FILE): whitelist, one color per line, restricting the build to
//     those colors (like pre-splitting the GFA); NULL keeps all. A line naming a color
//     selects it exactly; a line with no '#' selects every haplotype of that sample.
//     Blank lines and lines starting with '#' are ignored.
//   Names are computed over ALL walks first, so a color keeps the name a full-panel
//     build gives it (the "#hap" suffix follows the whole graph, not the surviving
//     subset); otherwise existing indexes would stop matching by name.
//
// gpeb parameters:
//   min_gap (bp): after emitting a pebble, skip further ones until the walk advances
//     >= min_gap bp. Only thins the output; pebbling decisions are unchanged. 0 = off.
//   min_complexity (0..1): drop pebbles whose span-bp window has a distinct-4-mer
//     fraction below this (tandem repeats, homopolymers). Applied before ranking. 0 = off.
//   early (--early-cut): slide each merge guardian back along the forced corridor
//     (in/out-degree 1) to the divergence edge instead of the re-merge edge. Composes
//     with every merge mode. The slide refuses to cross a vertex on a cycle: starting
//     the segment earlier could let the walk re-enter it and fire a revisit-cut the late
//     anchor never sees, changing the guardian count.
//     INVARIANT: sliding never changes which merges are cut, so the count is the same
//     and guardians only move backward (min_gap and min_complexity are position-sensitive
//     and may drop a different subset).
//   fast (default): cut at any in-degree >= 2 vertex; no dominators, O(V+E). Sound on
//     any graph but not minimal: a merge can have indeg >= 2 in the graph while the
//     anchor cannot reach one of its in-edges.
//   minimal (--minimal): the last-edge lemma test: is an in-edge source other than the
//     walk's own reachable from the anchor while AVOIDING the merge (the exclusion keeps
//     it correct on a cycle)? Three tiers: indeg < 2 in O(1); an O(1) screen off the
//     super-source dominator tree; then a backward search (back_reach). Exact and
//     minimal on any graph.
//   acyclic (--acyclic, DAG only): decide the merge by exact path counting clamped at 2,
//     one forward sweep per walk over a topological numbering (see swp_t). Same guardians
//     as --minimal on a DAG; refuses a cyclic graph. Overrides fast.

// ---- oriented-vertex sequence helpers --------------------------------------

static inline const char *nseq(const bg_t *g, uint32_t v) {
  return (v < g->n) ? g->seq[v >> 1] : NULL;
}
static inline size_t nsl(const bg_t *g, uint32_t v) {
  return (v < g->n) ? g->slen[v >> 1] : 0;
}

// fputs_sanitized: write s, replacing '|' with repl (path names may contain '|')
static void wsan(const char *s, char repl, FILE *out) {
  if (!s) return;
  for (; *s; s++) fputc(*s == '|' ? repl : *s, out);
}

// fputs_node_name: "<seg_name><+|->"
static void wname(const bg_t *g, uint32_t v, FILE *out) {
  const char *nm = (v < g->n) ? g->sname[v >> 1] : "?";
  fputs(nm, out);
  fputc((v & 1) ? '-' : '+', out);
}

// write_oriented_suffix: `take` bp from the END of v's oriented sequence
static void wsuf(const bg_t *g, uint32_t v, uint32_t take, char *dst) {
  const char *s = nseq(g, v);
  if (!s || take == 0) return;
  size_t L = strlen(s);
  if (!(v & 1)) memcpy(dst, s + (L - take), take);
  else for (uint32_t k = 0; k < take; k++)
    dst[k] = (char)gfa_comp_table[(uint8_t)s[take - 1 - k]];
}

// write_oriented_prefix: `take` bp from the START of v's oriented sequence
static void wpre(const bg_t *g, uint32_t v, uint32_t take, char *dst) {
  const char *s = nseq(g, v);
  if (!s || take == 0) return;
  size_t L = strlen(s);
  if (!(v & 1)) memcpy(dst, s, take);
  else for (uint32_t k = 0; k < take; k++)
    dst[k] = (char)gfa_comp_table[(uint8_t)s[L - 1 - k]];
}

// build_win: fill `buf` (capacity >= 2*span) with the `span`-bp window flanking
// the edge at walk slot i -> i+1 (left = suffix ending at head, right = prefix
// from tail), crossing node boundaries. Returns the byte length. Shared by the
// FASTA dumps, the complexity filter and the k-mer passes so all see the same
// sequence. lf_out (may be NULL) gets the junction offset: it sits at buf[lf].
static size_t build_win(const bg_t *g, uint32_t pid, uint32_t i, uint32_t span,
                        char *buf, size_t *lf_out) {
  const uint32_t *path = g->pv[pid];
  uint32_t n = g->pn[pid];

  size_t lf = 0;                                       // left: end of head back
  for (int32_t j = (int32_t)i; j >= 0 && lf < span; j--) {
    size_t L = nsl(g, path[j]);
    if (L == 0) continue;
    size_t take = (L < span - lf) ? L : (span - lf);
    wsuf(g, path[j], (uint32_t)take, buf + (span - lf - take));
    lf += take;
  }
  memmove(buf, buf + (span - lf), lf);

  size_t rf = 0;                                       // right: start of tail fwd
  for (uint32_t j = i + 1; j < n && rf < span; j++) {
    size_t L = nsl(g, path[j]);
    if (L == 0) continue;
    size_t take = (L < span - rf) ? L : (span - rf);
    wpre(g, path[j], (uint32_t)take, buf + lf + rf);
    rf += take;
  }
  if (lf_out) *lf_out = lf;
  return lf + rf;
}

// kmer_complexity: fraction of DISTINCT 4-mers over s (distinct/total). Low for
// tandem repeats / homopolymers (telomere, satellite, microsat), ~1 for ordinary
// sequence. Non-ACGT breaks the k-mer run. Returns 1.0 when too short to judge.
static double kmer_complexity(const char *s, size_t len) {
  if (len < 4) return 1.0;
  uint8_t seen[32];                                    // 256-bit set over 4-mers
  memset(seen, 0, sizeof seen);
  uint32_t distinct = 0, total = 0, kmer = 0;
  int run = 0;
  for (size_t i = 0; i < len; i++) {
    int b;
    switch (s[i]) {
      case 'A': case 'a': b = 0; break;
      case 'C': case 'c': b = 1; break;
      case 'G': case 'g': b = 2; break;
      case 'T': case 't': b = 3; break;
      default:  b = -1;
    }
    if (b < 0) { run = 0; continue; }
    kmer = ((kmer << 2) | (uint32_t)b) & 0xFF;
    if (++run >= 4) {
      total++;
      uint8_t m = (uint8_t)(1u << (kmer & 7));
      if (!(seen[kmer >> 3] & m)) { seen[kmer >> 3] |= m; distinct++; }
    }
  }
  return total ? (double)distinct / (double)total : 1.0;
}

// ---- Cooper-Harvey-Kennedy dominators ----------------------------------------
// Iterative, O(V+E) -- no recursion (a whole-chromosome backbone is hundreds of
// thousands of levels deep).

// dom_intersect: two-finger climb up the partial idom tree via postorder nums.
static uint32_t isect(uint32_t a, uint32_t b, const uint32_t *idom,
                      const uint32_t *post) {
  while (a != b) {
    while (post[a] < post[b]) a = idom[a];
    while (post[b] < post[a]) b = idom[b];
  }
  return a;
}

// ---- guardian pebbling ------------------------------------------------------

static void gadd(gset_t *gs, uint32_t head, uint32_t tail, uint32_t rank,
                 uint32_t pid, uint32_t pos) {
  if (gs->n == gs->cap) {
    gs->cap = gs->cap ? gs->cap * 2 : 16;
    gs->p = realloc(gs->p, gs->cap * sizeof(gpb_t));
  }
  gpb_t *q = &gs->p[gs->n++];
  q->head = head; q->tail = tail; q->rank = rank; q->pid = pid; q->pos = pos;
}

// peb_read: rebuild a guardian set from a .peb instead of re-running the cut. The
// index build needs only (pid, pos) of each pebble, so one cut can serve many index
// builds. head/tail are recomputed from pv[pid][pos..pos+1] rather than parsed, so they
// cannot disagree with the walk.
static gset_t *peb_read(const char *path, const bg_t *g) {
  gzFile fp = gzopen(path, "r");
  if (!fp) { fprintf(stderr, "peb: cannot open %s\n", path); return NULL; }
  gset_t *gs = calloc(1, sizeof(gset_t));
  size_t cap = 1 << 16; char *line = malloc(cap);
  size_t nline = 0, bad = 0;
  while (gzgets(fp, line, (int)cap)) {
    if (line[0] != '>') continue;
    nline++;
    // >>head>tail|hm|color|rank|pid|pos -- take the last two fields.
    char *f[8]; int nf = 0;
    for (char *s = line; *s && nf < 8; ) {
      f[nf++] = s;
      char *e = strchr(s, '|');
      if (!e) break;
      *e = 0; s = e + 1;
    }
    if (nf < 6) { bad++; continue; }
    uint32_t rank = (uint32_t)strtoul(f[3], NULL, 10);
    uint32_t pid  = (uint32_t)strtoul(f[4], NULL, 10);
    uint32_t pos  = (uint32_t)strtoul(f[5], NULL, 10);
    // Bounds are checked: a .peb from another graph would index out of range.
    // Same guardian invariant as gdump: pos+1 < pn.
    if (pid >= g->np || pos + 1 >= g->pn[pid]) { bad++; continue; }
    gadd(gs, g->pv[pid][pos], g->pv[pid][pos + 1], rank, pid, pos);
  }
  free(line); gzclose(fp);
  if (bad) {
    fprintf(stderr, "peb: %s: %zu of %zu records unusable -- refusing\n",
            path, bad, nline);
    fprintf(stderr, "     a .peb needs the |pid|pos fields and must come from"
                    " THIS graph; re-cut it with a current gpeb\n");
    free(gs->p); free(gs); return NULL;
  }
  fprintf(stderr, "peb: %s -> %u guardians (cut not re-run)\n", path, gs->n);
  return gs;
}

// back_reach: does `a` reach `w` WITHOUT passing through `blocked`? Backward search
// from w over the predecessor CSR; the caller passes the merge as `blocked`, since a route
// through the merge is not a simple path to it. Returns 1 = reaches, 0 = does not.
// Two O(1) prunes off the super-source dominator tree, valid only if `a` is S-reachable
// (`domok`):
//   gid[u] == UNREACH             => a cannot reach u: drop
//   u dominates a, a on no cycle  => a cannot reach u (it would close a cycle): drop
// There is deliberately no "a dominates u => a reaches u" prune: that route may pass
// through `blocked`.
// ponytail: unbounded, so one search can cost O(V+E).
static int back_reach(uint32_t a, uint32_t w, uint32_t blocked,
                      const uint32_t *phead, const uint32_t *pred,
                      const uint32_t *gid, const uint32_t *tin, const uint32_t *tout,
                      const uint8_t *incyc, uint32_t *lvis, uint32_t gen,
                      uint32_t *stk, uint64_t *visited) {
  uint32_t sp = 0, seen = 0;
  // Both prunes assume `a` is S-reachable. If every vertex has an in-edge there is no
  // source, gid[] is UNREACH everywhere and the prunes would discard the whole search (an
  // under-cut), so skip them.
  int domok = (gid[a] != UNREACH);
  if (w == blocked) return 0;
  stk[sp++] = w; lvis[w] = gen;
  while (sp) {
    uint32_t u = stk[--sp];
    if (u == a) { *visited += seen; return 1; }
    for (uint32_t pi = phead[u]; pi < phead[u + 1]; pi++) {
      uint32_t x = pred[pi];
      if (x == a) { *visited += seen; return 1; }   // must precede the prunes: a dominates itself
      if (x == blocked) continue;     // a route through the merge is not a simple path to it
      if (lvis[x] == gen) continue;
      if (domok && gid[x] == UNREACH) continue;
      if (domok && !incyc[a] && tin[x] <= tin[a] && tout[a] <= tout[x]) continue;
      seen++; lvis[x] = gen; stk[sp++] = x;
    }
  }
  *visited += seen;
  return 0;
}

// scc_mark: incyc[v] = 1 iff v lies on a cycle (nontrivial SCC or self-loop).
// Iterative Tarjan.
static void scc_mark(const bg_t *g, uint8_t *incyc) {
  uint32_t n = g->n, UN = UINT32_MAX;
  uint32_t *idx = malloc(n * sizeof(uint32_t)), *low = malloc(n * sizeof(uint32_t));
  uint32_t *cst = malloc(n * sizeof(uint32_t)); // Tarjan component stack
  uint32_t *dfs = malloc(n * sizeof(uint32_t)); // explicit DFS stack: vertices
  uint32_t *it  = malloc(n * sizeof(uint32_t)); // ...and per-frame next-succ index
  uint32_t *buf = malloc(n * sizeof(uint32_t)); // scratch to size a popped SCC
  uint8_t  *onc = calloc(n, 1);
  for (uint32_t i = 0; i < n; i++) { idx[i] = UN; incyc[i] = 0; }
  uint32_t cnt = 0, csp = 0;
  for (uint32_t s = 0; s < n; s++) {
    if (idx[s] != UN) continue;
    uint32_t dsp = 0;
    dfs[0] = s; it[0] = 0; idx[s] = low[s] = cnt++; cst[csp++] = s; onc[s] = 1;
    for (;;) {
      uint32_t v = dfs[dsp];
      if (it[dsp] < g->adjn[v]) {
        uint32_t w = g->adj[v][it[dsp]++];
        if (w >= n) continue;
        if (w == v) { incyc[v] = 1; continue; }         // self-loop
        if (idx[w] == UN) {
          idx[w] = low[w] = cnt++; cst[csp++] = w; onc[w] = 1;
          dfs[++dsp] = w; it[dsp] = 0;
        } else if (onc[w] && idx[w] < low[v]) low[v] = idx[w];
      } else {
        if (low[v] == idx[v]) {
          uint32_t sz = 0, w;
          do { w = cst[--csp]; onc[w] = 0; buf[sz++] = w; } while (w != v);
          if (sz > 1) for (uint32_t i = 0; i < sz; i++) incyc[buf[i]] = 1;
        }
        if (dsp == 0) break;
        uint32_t par = dfs[dsp - 1];
        if (low[v] < low[par]) low[par] = low[v];
        dsp--;
      }
    }
  }
  free(idx); free(low); free(cst); free(dfs); free(it); free(buf); free(onc);
}

// pred_csr: predecessor lists of g, self-loops dropped. The preds of v are
// pred[phead[v] .. phead[v+1]), so its in-degree is phead[v+1] - phead[v].
static void pred_csr(const bg_t *g, uint32_t **phead_out, uint32_t **pred_out) {
  uint32_t n = g->n, *phead = calloc(n + 1, sizeof(uint32_t));
  for (uint32_t v = 0; v < n; v++)
    for (uint32_t a = 0; a < g->adjn[v]; a++) { uint32_t w = g->adj[v][a];
      if (w < n && w != v) phead[w + 1]++; }
  for (uint32_t v = 0; v < n; v++) phead[v + 1] += phead[v];
  uint32_t *pred = malloc((phead[n] ? phead[n] : 1) * sizeof(uint32_t));
  uint32_t *fill = malloc((n ? n : 1) * sizeof(uint32_t));
  memcpy(fill, phead, n * sizeof(uint32_t));
  for (uint32_t v = 0; v < n; v++)
    for (uint32_t a = 0; a < g->adjn[v]; a++) { uint32_t w = g->adj[v][a];
      if (w < n && w != v) pred[fill[w]++] = v; }
  free(fill);
  *phead_out = phead; *pred_out = pred;
}

// gdom_build: ONE super-source dominator tree over G_P, built once. Virtual source
// S=n links to every in-degree-0 vertex. Fills idom[0..n-1] (a source's idom is S;
// unreachable => UNREACH) and idom-tree enter/exit times, so "u dominates v" ==
// tin[u]<=tin[v] && tout[v]<=tout[u]. `ghead` is pred_csr's phead of g (marks the
// sources). Iterative CHK fixpoint, O(V+E).
static void gdom_build(const bg_t *g, const uint32_t *ghead,
                       uint32_t *idom, uint32_t *tin, uint32_t *tout) {
  uint32_t n = g->n, S = n, N = n + 1;
  uint32_t nsrc = 0;
  for (uint32_t v = 0; v < n; v++) if (ghead[v + 1] == ghead[v]) nsrc++;
  uint32_t *srcs = malloc((nsrc ? nsrc : 1) * sizeof(uint32_t));
  { uint32_t j = 0; for (uint32_t v = 0; v < n; v++) if (ghead[v + 1] == ghead[v]) srcs[j++] = v; }
  #define GD_DEG(v) ((v) == S ? nsrc : g->adjn[(v)])
  #define GD_ADJ(v, a) ((v) == S ? srcs[(a)] : g->adj[(v)][(a)])

  uint32_t *order = malloc(N * sizeof(uint32_t)), *post = malloc(N * sizeof(uint32_t));
  uint32_t *stk = malloc(N * sizeof(uint32_t)), *it = malloc(N * sizeof(uint32_t));
  uint8_t  *vis = calloc(N, 1);
  uint32_t R = 0, sp = 0; vis[S] = 1; stk[0] = S; it[0] = 0; sp = 1;
  while (sp > 0) { uint32_t v = stk[sp - 1];
    if (it[sp - 1] < GD_DEG(v)) { uint32_t w = GD_ADJ(v, it[sp - 1]); it[sp - 1]++;
      if (w < N && !vis[w]) { vis[w] = 1; stk[sp] = w; it[sp] = 0; sp++; } }
    else { post[v] = R; order[R++] = v; sp--; } }

  uint32_t *phead = calloc(N + 1, sizeof(uint32_t)), *fill = malloc(N * sizeof(uint32_t));
  for (uint32_t i = 0; i < R; i++) { uint32_t v = order[i];
    for (uint32_t a = 0; a < GD_DEG(v); a++) { uint32_t w = GD_ADJ(v, a);
      if (w < N && vis[w] && w != v) phead[w + 1]++; } }
  for (uint32_t i = 0; i < N; i++) phead[i + 1] += phead[i];
  uint32_t tot = phead[N]; uint32_t *pred = malloc((tot ? tot : 1) * sizeof(uint32_t));
  memcpy(fill, phead, N * sizeof(uint32_t));
  for (uint32_t i = 0; i < R; i++) { uint32_t v = order[i];
    for (uint32_t a = 0; a < GD_DEG(v); a++) { uint32_t w = GD_ADJ(v, a);
      if (w < N && vis[w] && w != v) pred[fill[w]++] = v; } }

  uint32_t *id = malloc(N * sizeof(uint32_t));
  for (uint32_t v = 0; v < N; v++) id[v] = UNREACH;  // S-unreachable stay UNREACH
  id[S] = S;
  int changed = 1;
  while (changed) { changed = 0;
    for (int32_t i = (int32_t)R - 1; i >= 0; i--) { uint32_t b = order[i];
      if (b == S) continue;
      uint32_t ni = UNREACH;
      for (uint32_t pi = phead[b]; pi < phead[b + 1]; pi++) { uint32_t p = pred[pi];
        if (id[p] == UNREACH) continue;
	ni = (ni == UNREACH) ? p : isect(p, ni, id, post); }
      if (ni != UNREACH && id[b] != ni) { id[b] = ni; changed = 1; } } }

  uint32_t *chead = calloc(N + 1, sizeof(uint32_t)), *child = malloc(N * sizeof(uint32_t));
  for (uint32_t i = 0; i < R; i++) { uint32_t x = order[i]; if (x != S) chead[id[x] + 1]++; }
  for (uint32_t i = 0; i < N; i++) chead[i + 1] += chead[i];
  memcpy(fill, chead, N * sizeof(uint32_t));
  for (uint32_t i = 0; i < R; i++) { uint32_t x = order[i]; if (x != S) child[fill[id[x]]++] = x; }
  uint32_t *ti = calloc(N, sizeof(uint32_t)), *to = calloc(N, sizeof(uint32_t));
  uint32_t timer = 0; sp = 0; stk[0] = S; it[0] = chead[S]; sp = 1; ti[S] = timer++;
  while (sp > 0) { uint32_t v = stk[sp - 1];
    if (it[sp - 1] < chead[v + 1]) { uint32_t c = child[it[sp - 1]++];
      ti[c] = timer++; stk[sp] = c; it[sp] = chead[c]; sp++; }
    else { to[v] = timer++; sp--; } }

  for (uint32_t v = 0; v < n; v++) { idom[v] = id[v]; tin[v] = ti[v]; tout[v] = to[v]; }

  free(srcs); free(order); free(post); free(stk); free(it); free(vis);
  free(phead); free(fill); free(pred); free(id); free(chead); free(child); free(ti); free(to);
  #undef GD_DEG
  #undef GD_ADJ
}

static uint32_t uf_find(uint32_t *uf, uint32_t x) {   // union-find, path halving
  while (uf[x] != x) { uf[x] = uf[uf[x]]; x = uf[x]; }
  return x;
}

// ---- forward single-sweep path counting (--acyclic) -------------------------
//
// The merge test on a DAG is a path count clamped at 2. Vertices are numbered
// topologically and swept ONCE per walk, carrying cnt[u] = min(2, #paths anchor->u).
// Re-anchoring is free: a vertex below the new anchor is dropped by one comparison
// (topo[w] >= A), and one above it has not been swept yet, so nothing is recomputed.
// A walk costs one sweep of its topological span however many guardians it collects.
// This mirrors the paper's backward sweep forwards, the direction gpeb's scan runs in
// (hence the pebble sits on the edge entering the merge, not the one leaving the
// divergence).
typedef struct {
  uint32_t *topo, *order;
  uint32_t *phead, *pred;
  uint8_t  *cnt;
  uint32_t A, pos;
} swp_t;

static void swp_free(swp_t *s) {
  if (!s) return;
  free(s->topo); free(s->order); free(s->phead); free(s->pred); free(s->cnt); free(s);
}

// Kahn's algorithm one weakly-connected component at a time, so each component is a
// contiguous block of `order` and a walk sweeps only its own component's span. Returns
// NULL if the graph has a cycle (no topological order exists).
static swp_t *swp_build(const bg_t *g) {
  uint32_t n = g->n, *phead, *pred;
  pred_csr(g, &phead, &pred);
  uint32_t *ind = malloc((n ? n : 1) * sizeof(uint32_t));
  uint32_t *uf  = malloc((n ? n : 1) * sizeof(uint32_t));
  for (uint32_t v = 0; v < n; v++) uf[v] = v;
  for (uint32_t v = 0; v < n; v++) {
    ind[v] = phead[v + 1] - phead[v];
    for (uint32_t pi = phead[v]; pi < phead[v + 1]; pi++) {
      uint32_t rv = uf_find(uf, pred[pi]), rw = uf_find(uf, v); if (rv != rw) uf[rv] = rw; }
  }

  uint32_t *cstart = calloc(n + 1, sizeof(uint32_t));      // vertices bucketed by WCC
  uint32_t *byc = malloc(n * sizeof(uint32_t));
  for (uint32_t v = 0; v < n; v++) cstart[uf_find(uf, v) + 1]++;
  for (uint32_t v = 0; v < n; v++) cstart[v + 1] += cstart[v];
  { uint32_t *fill = malloc(n * sizeof(uint32_t));
    memcpy(fill, cstart, n * sizeof(uint32_t));
    for (uint32_t v = 0; v < n; v++) byc[fill[uf_find(uf, v)]++] = v;
    free(fill); }

  uint32_t *order = malloc(n * sizeof(uint32_t)), R = 0, qh = 0;
  for (uint32_t r = 0; r < n; r++) {
    if (cstart[r + 1] == cstart[r]) continue;
    for (uint32_t i = cstart[r]; i < cstart[r + 1]; i++)
      if (!ind[byc[i]]) order[R++] = byc[i];
    while (qh < R) {
      uint32_t v = order[qh++];
      for (uint32_t a = 0; a < g->adjn[v]; a++) { uint32_t w = g->adj[v][a];
        if (w < n && w != v && --ind[w] == 0) order[R++] = w; }
    }
  }
  free(cstart); free(byc); free(uf); free(ind);
  if (R < n) { free(order); free(phead); free(pred); return NULL; } // a cycle held vertices back

  swp_t *s = calloc(1, sizeof(swp_t));
  s->order = order;
  s->topo = malloc(n * sizeof(uint32_t));
  for (uint32_t i = 0; i < n; i++) s->topo[order[i]] = i;
  s->cnt = malloc(n);
  s->phead = phead; s->pred = pred;
  return s;
}

// Re-anchor. Nothing is cleared: cnt[] is only read at positions in [A, pos), and
// swp_count writes every one of them after this reset.
static inline void swp_reset(swp_t *s, uint32_t v) {
  s->A = s->topo[v]; s->cnt[v] = 1; s->pos = s->A + 1;
}

// min(2, #paths anchor->x). Targets must arrive in non-decreasing topological
// order, which a walk does: its edges point forward in any topological numbering.
static inline uint32_t swp_count(swp_t *s, uint32_t x) {
  uint32_t tx = s->topo[x];
  while (s->pos <= tx) {
    uint32_t u = s->order[s->pos], c = 0;
    for (uint32_t pi = s->phead[u]; pi < s->phead[u + 1]; pi++) {
      uint32_t w = s->pred[pi];
      if (s->topo[w] < s->A) continue;                     // below the anchor: no route
      c += s->cnt[w];
      if (c >= 2) { c = 2; break; }
    }
    s->cnt[u] = (uint8_t)c;
    s->pos++;
  }
  return s->cnt[x];
}

// walk_colors: color[w] for every walk. Walks sharing a pname POINTER are one
// color -- bg_build interns names -- numbered in first-seen order. *cname (may be
// NULL) gets the distinct names in that order, owned by the caller.
static uint32_t *walk_colors(const bg_t *g, const char ***cname, uint32_t *nc_out) {
  uint32_t np = g->np, nc = 0;
  uint32_t *color = malloc((np ? np : 1) * sizeof(uint32_t));
  const char **cptr = malloc((np ? np : 1) * sizeof(char *));
  for (uint32_t w = 0; w < np; w++) {
    const char *s = g->pname ? g->pname[w] : NULL;
    uint32_t c = 0; while (c < nc && cptr[c] != s) c++;
    if (c == nc) cptr[nc++] = s;
    color[w] = c;
  }
  if (cname) *cname = cptr; else free(cptr);
  *nc_out = nc;
  return color;
}

gset_t *gpeb(const bg_t *g, int fast, uint32_t min_gap, uint32_t span,
             double min_complexity, int early, int acyclic) {
  gset_t *gs = calloc(1, sizeof(gset_t));
  if (!g || g->np == 0 || !g->adj) return gs;
  uint32_t n = g->n;
  if (acyclic) fast = 0;

  // --acyclic: one forward sweep per walk over a topological numbering; DAG-only,
  // enforced.
  swp_t *sw = NULL;
  if (acyclic && !(sw = swp_build(g))) {
    fprintf(stderr, "gpeb: --acyclic needs an acyclic graph, and this one has a"
                    " cycle. Use --minimal instead.\n");
    free(gs); return NULL;
  }

  // --minimal: by the last-edge lemma the merge test is reachability from the anchor to
  // an in-edge source other than the walk's own, AVOIDING the merge; a local search answers
  // it instead of a dominator row per anchor.
  int useloc = !fast && !acyclic;
  uint8_t  *seen = calloc(n, 1);
  uint32_t *sl   = malloc(n * sizeof(uint32_t));
  char *wbuf = (min_complexity > 0.0) ? malloc(2 * (span ? span : 1)) : NULL;

  // Predecessor CSR: in-degree feeds the fast rule, the indeg<2 exit of --minimal and
  // --early-cut's corridor test; pred feeds the reduction. --acyclic keeps its own.
  uint32_t *phead = NULL, *pred = NULL;
  if (fast || useloc || early) pred_csr(g, &phead, &pred);
  #define INDEG(v) (phead[(v) + 1] - phead[(v)])
  // incyc marks cycle vertices: --early-cut will not slide across one, and the local
  // search's prune is valid only for an anchor on no cycle.
  uint8_t  *incyc = (early || useloc) ? malloc(n) : NULL;
  if (incyc) scc_mark(g, incyc);
  // One super-source dominator tree, built once, for the O(1) in-edge screen and the
  // local search's prunes.
  uint32_t *gid = NULL, *gtin = NULL, *gtout = NULL;
  if (useloc) {
    gid = malloc(n * sizeof(uint32_t));
    gtin = malloc(n * sizeof(uint32_t));
    gtout = malloc(n * sizeof(uint32_t));
    gdom_build(g, phead, gid, gtin, gtout);
  }
  // --minimal census (reported on stderr): which branch decided each merge test.
  uint64_t nb_nomerge = 0, nb_screen = 0, nb_local = 0, nb_visited = 0;
  uint32_t *lvis = NULL, *lstk = NULL, lgen = 0;
  if (useloc) {
    lvis = calloc(n, sizeof(uint32_t)); lstk = malloc((size_t)n * sizeof(uint32_t));
  }

  // Color aggregation: walks that share a sample name are one color, with a
  // single continuous rank sequence. Walks are processed color-major (counting
  // sort) so each color's pebbles are contiguous and rank-ascending.
  uint32_t np = g->np, nc = 0;
  uint32_t *color = walk_colors(g, NULL, &nc);
  uint32_t *crank = calloc(nc ? nc : 1, sizeof(uint32_t)); // per-color 1-based rank
  uint32_t *ord   = malloc(np * sizeof(uint32_t));
  uint32_t *coff  = calloc(nc + 1, sizeof(uint32_t));
  for (uint32_t k = 0; k < np; k++) coff[color[k] + 1]++;
  for (uint32_t c = 0; c < nc; c++) coff[c + 1] += coff[c];
  for (uint32_t k = 0; k < np; k++) ord[coff[color[k]]++] = k; // stable: walks grouped by color

  for (uint32_t oi = 0; oi < np; oi++) {
    uint32_t k = ord[oi];
    const uint32_t *p = g->pv[k];
    uint32_t pn = g->pn[k];
    if (pn < 1) continue;
    uint32_t c = color[k];

    uint32_t anchor = p[0];
    uint32_t apos   = 0;                      // ...and its position along the walk
    uint32_t ns = 0;
    if (anchor < n) { seen[anchor] = 1; sl[ns++] = anchor; }
    if (sw && anchor < n) swp_reset(sw, anchor);

    // bp offset up to and including the current head; tracked only when coalescing.
    uint32_t wbp = min_gap ? (uint32_t)nsl(g, p[0]) : 0;
    uint32_t last_bp = 0; int have = 0;       // last emitted pebble's bp (this walk)

    for (uint32_t s = 1; s < pn; s++) {
      uint32_t nxt = p[s];
      int rep = (nxt < n) && seen[nxt];       // revisit => not simple
      int mul;                                // >=2 routes from anchor?
      if (nxt >= n) mul = 1;
      else if (fast) mul = INDEG(nxt) >= 2;   // conservative: any merge in the graph
      else if (useloc) {
        if (INDEG(nxt) < 2) { mul = 0; nb_nomerge++; }
        else {
          // Screen the in-edges in O(1) off the dominator tree; only what it leaves open goes
          // to the backward search. By the last-edge lemma the question is pure reachability: is
          // an in-edge source other than the walk's own reachable from the anchor?
          //   w == anchor                             => yes, CUT
          //   w dominates anchor, anchor on no cycle  => anchor cannot reach w: drop that edge
          // The !incyc guard keeps the second rule exact: w dominating the anchor gives a path
          // w->anchor, so anchor->w would close a cycle through the anchor; without the guard it
          // under-cuts. If every other in-edge is dropped, no second route exists: NO CUT.
          int yes = 0, unknown = 0;
          for (uint32_t pi = phead[nxt]; pi < phead[nxt + 1]; pi++) {
            uint32_t w = pred[pi];
            if (w == p[s - 1]) continue;
            if (w == anchor) { yes = 1; break; }   // direct edge + the walk's own arrival: two simple paths; must precede the prune
            // No other O(1) "yes" rule: the anchor dominating w says it reaches w, not that it
            // does so AVOIDING nxt.
            if (gid[anchor] != UNREACH && !incyc[anchor] &&
                gtin[w] <= gtin[anchor] && gtout[anchor] <= gtout[w]) continue;
              // gid[anchor] guard: see back_reach (a source-less graph leaves the whole tree UNREACH).
            unknown = 1;
          }
          if (yes || !unknown) { mul = yes; nb_screen++; }
          else {
            mul = 0; nb_local++;
            for (uint32_t pi = phead[nxt]; pi < phead[nxt + 1] && !mul; pi++)
              if (pred[pi] != p[s - 1])
                mul = back_reach(anchor, pred[pi], nxt, phead, pred, gid, gtin, gtout,
                                 incyc, lvis, ++lgen, lstk, &nb_visited);
          }
        }
      }
      else                                    // --acyclic: running sweep, no recompute
        mul = (anchor < n) && swp_count(sw, nxt) >= 2;
      if (rep || mul) {
        // Guardian goes on edge (p[gp], p[gp+1]); the next segment anchors at p[gp+1].
        // Default: cut as LATE as possible, on the edge entering the merge (gp = s-1).
        uint32_t gp = s - 1;

        // --early-cut: slide the guardian back along the forced corridor (in-degree 1 and
        // out-degree 1) preceding the merge, to the divergence edge. Safe in every merge mode:
        //  (1) merge decisions are unaffected: every corridor vertex has outdeg 1, so from the
        //      slid-back anchor the only route forward lands on the merge and simple paths
        //      biject with those from the merge. The absorbed span is fast-forwarded below, not
        //      re-tested.
        //  (2) revisit decisions are unaffected only if no corridor vertex lies on a cycle
        //      (the !incyc guard): sliding starts the segment earlier, so `seen` covers more of
        //      the walk and a cycle vertex could fire a revisit-cut the late anchor never sees.
        // INVARIANT (asserted in GPEB_TEST on random DAGs and cyclic graphs): sliding never
        // changes WHICH merges are cut, only moves guardians backward. With min_gap or
        // min_complexity set the kept subset may differ (both are position-sensitive).
        if (early && mul && !rep) {
          int32_t i = (int32_t)s - 1;
          while (i > (int32_t)apos && INDEG(p[i]) == 1 && g->adjn[p[i]] == 1
                 && !incyc[p[i]]) i--;
          if ((uint32_t)i + 1 < s) gp = (uint32_t)i;   // slide back to the corridor entrance
        }

        uint32_t ghead = p[gp], gtail = p[gp + 1];
        int emit = (min_gap == 0 || !have || wbp - last_bp >= min_gap);
        if (emit && wbuf) {
          size_t wl = build_win(g, k, gp, span, wbuf, NULL);
          if (kmer_complexity(wbuf, wl) < min_complexity) emit = 0;
        }
        if (emit) {
          gadd(gs, ghead, gtail, ++crank[c], k, gp); // ranked per color, not per walk
          last_bp = wbp; have = 1;
        }
        for (uint32_t i = 0; i < ns; i++) seen[sl[i]] = 0;
        ns = 0;
        anchor = gtail; apos = gp + 1;        // re-anchor at the guardian's tail
        // The sweep re-anchors at the MERGE, not at gtail: it is single-pass and cannot go
        // back, and from either end of a forced corridor the simple paths onward are in
        // bijection.
        if (sw && nxt < n) swp_reset(sw, nxt);
        if (anchor < n) { seen[anchor] = 1; sl[ns++] = anchor; }
        // Forced from the new anchor (corridor + merge): absorb without re-testing.
        for (uint32_t j = apos + 1; j <= s; j++)
          if (p[j] < n && !seen[p[j]]) { seen[p[j]] = 1; sl[ns++] = p[j]; }
      } else if (nxt < n && !seen[nxt]) {
        seen[nxt] = 1; sl[ns++] = nxt;
      }
      if (min_gap && nxt < n) wbp += (uint32_t)nsl(g, nxt);
    }
    for (uint32_t i = 0; i < ns; i++) seen[sl[i]] = 0;
  }

  if (useloc && n >= (1u << 16)) {     // skip toy graphs: GPEB_TEST calls gpeb() thousands of times
    uint64_t dec = nb_nomerge + nb_screen + nb_local;
    fprintf(stderr, "cut: %llu merge tests -- %llu indeg<2, %llu in-edge screen (%.2f%%),"
            " %llu local search (%.2f%%, %llu vertices)\n",
            (unsigned long long)dec, (unsigned long long)nb_nomerge,
            (unsigned long long)nb_screen, dec ? 100.0 * nb_screen / dec : 0.0,
            (unsigned long long)nb_local, dec ? 100.0 * nb_local / dec : 0.0,
            (unsigned long long)nb_visited);
  }
  #undef INDEG
  free(lvis); free(lstk);
  free(gid); free(gtin); free(gtout);
  free(phead); free(pred);
  free(seen); free(sl); free(wbuf); free(incyc);
  swp_free(sw);
  free(color); free(crank); free(ord); free(coff);
  return gs;
}

void gpeb_free(gset_t *gs) { if (gs) { free(gs->p); free(gs); } }

// gdump: one FASTA record per pebble, walk order,
//   >>head>tail|hm|<walk_name>|<rank>  +  span-bp window across the junction.
int gdump(const bg_t *g, const gset_t *gs, FILE *out, uint32_t span) {
  if (!g || !gs || !out) return -1;
  char *buf = malloc(2 * (span ? span : 1));
  for (uint32_t i = 0; i < gs->n; i++) {
    const gpb_t *p = &gs->p[i];
    const char *cn = (g->pname && p->pid < g->np && g->pname[p->pid])
                         ? g->pname[p->pid] : "?";
    fputs(">>", out);
    wname(g, p->head, out);
    fputc('>', out);
    wname(g, p->tail, out);
    fputs("|hm|", out);
    wsan(cn, '_', out);
    // pid/pos close the record: --from-peb needs them to rebuild the guardian set.
    // (color, rank) cannot name the walk when one color spans several walks.
    fprintf(out, "|%u|%u|%u\n", p->rank, p->pid, p->pos);
    // exact edge slot => right window even if a cycle re-traverses the edge.
    // guardian invariant: pos = s-1 with s in [1,pn-1], so pos+1 < pn always.
    if (p->pid < g->np && p->pos + 1 < g->pn[p->pid])
      fwrite(buf, 1, build_win(g, p->pid, p->pos, span, buf, NULL), out);
    fputc('\n', out);
  }
  free(buf);
  return 0;
}

// qsort comparator: group pebbles by edge (head, then tail).
static int cmp_edge(const void *a, const void *b) {
  const gpb_t *x = a, *y = b;
  if (x->head != y->head) return x->head < y->head ? -1 : 1;
  if (x->tail != y->tail) return x->tail < y->tail ? -1 : 1;
  return 0;
}

// gdump_c: one record per (edge, span-window), listing every color that shares it.
// Sorts gs in place by edge, then sub-groups each edge-run by window BYTES (a sub-span
// head/tail spills into path-specific neighbours, so same edge != same window). Header:
//   >>head>tail|np|N|hm|name|rank(|hm|name|rank)*
// ponytail: O(R^2) window compares per edge-run (R = colors on it); sort by window if an
// edge is ever shared by thousands.
int gdump_c(const bg_t *g, gset_t *gs, FILE *out, uint32_t span) {
  if (!g || !gs || !out) return -1;
  qsort(gs->p, gs->n, sizeof(gpb_t), cmp_edge);

  uint32_t wcap = 2 * (span ? span : 1);
  char    *win  = NULL;   // the current edge-run's R windows, wcap bytes each
  size_t  *wl   = NULL;   // their lengths
  uint8_t *done = NULL;   // sub-grouping marks
  uint32_t cap  = 0;

  for (uint32_t lo = 0; lo < gs->n; ) {
    uint32_t hi = lo + 1;
    while (hi < gs->n && gs->p[hi].head == gs->p[lo].head
                      && gs->p[hi].tail == gs->p[lo].tail) hi++;
    uint32_t R = hi - lo;
    if (R > cap) {
      cap  = R;
      win  = realloc(win,  (size_t)cap * wcap);
      wl   = realloc(wl,   cap * sizeof *wl);
      done = realloc(done, cap);
    }
    for (uint32_t k = 0; k < R; k++) {
      const gpb_t *p = &gs->p[lo + k];
      done[k] = 0;
      wl[k] = (p->pid < g->np && p->pos + 1 < g->pn[p->pid])
                ? build_win(g, p->pid, p->pos, span, win + (size_t)k * wcap, NULL) : 0;
    }
    for (uint32_t k = 0; k < R; k++) {
      if (done[k]) continue;
      char *wk = win + (size_t)k * wcap;
      uint32_t np = 0;
      for (uint32_t j = k; j < R; j++)
        if (!done[j] && wl[j] == wl[k] && memcmp(win + (size_t)j * wcap, wk, wl[k]) == 0)
          np++;
      const gpb_t *pk = &gs->p[lo + k];
      fputs(">>", out); wname(g, pk->head, out);
      fputc('>', out);  wname(g, pk->tail, out);
      fprintf(out, "|np|%u", np);
      for (uint32_t j = k; j < R; j++) {
        if (done[j] || wl[j] != wl[k] || memcmp(win + (size_t)j * wcap, wk, wl[k]) != 0)
          continue;
        done[j] = 1;
        const gpb_t *pj = &gs->p[lo + j];
        const char *cn = (g->pname && pj->pid < g->np && g->pname[pj->pid])
                             ? g->pname[pj->pid] : "?";
        fputs("|hm|", out); wsan(cn, '_', out); fprintf(out, "|%u", pj->rank);
      }
      fputc('\n', out);
      fwrite(wk, 1, wl[k], out);
      fputc('\n', out);
    }
    lo = hi;
  }
  free(win); free(wl); free(done);
  return 0;
}

// ---- k-mer locus uniqueness of guardian windows -----------------------------
// How often does a guardian window carry a k-mer that pins a SINGLE place in the
// graph? A k-mer occurring elsewhere admits a false-positive match.
//
// "Elsewhere" is counted in LOCI, not occurrences: a guardian crossed by 24 haplotypes
// has each k-mer 24 times, always at the same graph position, which is one locus. An
// occurrence is keyed by (oriented vertex, offset) of the k-mer's first base.
//
// Only the k-mers of the guardian windows are inserted; the walks are then streamed
// once and each k-mer looked up. Memory O(#guardian k-mers), time O(total walk bp).


// Monotonic wall-clock seconds, for the per-phase timings of an index build.
static double now_s(void) {
  struct timespec t;
  clock_gettime(CLOCK_MONOTONIC, &t);
  return (double)t.tv_sec + (double)t.tv_nsec * 1e-9;
}

// ---- color sets -------------------------------------------------------------
// One bit per color over nw = ceil(nc/64) words; every operation is a word loop. For
// nc <= 64 the layout equals a single uint64_t mask, so .pex files need no version bump.
static inline uint32_t cm_nw(uint32_t nc) { return nc ? (nc + 63u) / 64u : 1u; }
static inline void cm_set(uint64_t *m, uint32_t c) { m[c >> 6] |= 1ULL << (c & 63); }
static inline int  cm_get(const uint64_t *m, uint32_t c) { return (m[c >> 6] >> (c & 63)) & 1; }
static inline void cm_and(uint64_t *a, const uint64_t *b, uint32_t nw)
  { for (uint32_t z = 0; z < nw; z++) a[z] &= b[z]; }
static inline void cm_cpy(uint64_t *a, const uint64_t *b, uint32_t nw)
  { memcpy(a, b, (size_t)nw * sizeof(uint64_t)); }
static inline int  cm_eq(const uint64_t *a, const uint64_t *b, uint32_t nw)
  { return !memcmp(a, b, (size_t)nw * sizeof(uint64_t)); }
static inline int  cm_any(const uint64_t *m, uint32_t nw)
  { for (uint32_t z = 0; z < nw; z++) if (m[z]) return 1; return 0; }
static inline uint32_t cm_pop(const uint64_t *m, uint32_t nw)
  { uint32_t t = 0; for (uint32_t z = 0; z < nw; z++) t += (uint32_t)__builtin_popcountll(m[z]); return t; }

#define KH_EMPTY 0
#define KH_WANT  1     // in a guardian window, not yet located in the walks
#define KH_ONE   2     // seen at exactly one locus so far
#define KH_MANY  3     // seen at two or more distinct loci
#define KH_IDX   4     // kept in the classifiable index, carries a color set
#define KH_DROP  5     // was KH_IDX, removed by --subsample; distinct from KH_MANY (ambiguous)

typedef struct { uint64_t *key, *loc; uint8_t *st; size_t mask, n; } kh_t;

// `loc` is needed only while locating k-mers in the walks (the first locus seen, so a
// second one can flip the entry to KH_MANY). Tables that never locate (window-dedup sets,
// indexes loaded from a .pex) leave it NULL, saving 8 bytes per slot.
static kh_t *kh_new(size_t want, int with_loc) {
  size_t cap = 16; while (cap < want * 2) cap <<= 1;   // <=50% load
  kh_t *h = calloc(1, sizeof(kh_t));
  h->key = malloc(cap * sizeof(uint64_t));
  h->loc = with_loc ? malloc(cap * sizeof(uint64_t)) : NULL;
  h->st  = calloc(cap, 1);
  h->mask = cap - 1;
  return h;
}
static void kh_free(kh_t *h) { if (h) { free(h->key); free(h->loc); free(h->st); free(h); } }

static inline uint64_t kh_mix(uint64_t x) {          // splitmix64 finalizer
  x ^= x >> 30; x *= 0xbf58476d1ce4e5b9ULL;
  x ^= x >> 27; x *= 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

// FNV-1a, chainable: pass FNV_SEED to start, the previous result to continue.
#define FNV_SEED 1469598103934665603ULL
static uint64_t fnv1a(const void *p, size_t n, uint64_t h) {
  const unsigned char *b = p;
  for (size_t i = 0; i < n; i++) { h ^= b[i]; h *= 1099511628211ULL; }
  return h;
}

// Grow and rehash. Linear probing loops forever on a full table, so the load factor is
// enforced (the number of window k-mers is not known in advance).
static void kh_grow(kh_t *h) {
  size_t ocap = h->mask + 1, ncap = ocap * 2;
  uint64_t *ok = h->key, *ol = h->loc; uint8_t *os = h->st;
  h->key = malloc(ncap * sizeof(uint64_t));
  h->loc = ol ? malloc(ncap * sizeof(uint64_t)) : NULL;
  h->st  = calloc(ncap, 1);
  h->mask = ncap - 1;
  for (size_t j = 0; j < ocap; j++) {
    if (os[j] == KH_EMPTY) continue;
    size_t i = kh_mix(ok[j]) & h->mask;
    while (h->st[i] != KH_EMPTY) i = (i + 1) & h->mask;
    h->key[i] = ok[j]; h->st[i] = os[j];
    if (ol) h->loc[i] = ol[j];
  }
  free(ok); free(ol); free(os);
}

// Slot of `key`, inserting it as KH_WANT when `insert`; SIZE_MAX if absent. The EMPTY
// check is an atomic load so lookups are well-defined while other threads update other
// slots' state; relaxed suffices because a slot is never re-emptied.
static size_t kh_slot(kh_t *h, uint64_t key, int insert) {
  for (;;) {
    size_t i = kh_mix(key) & h->mask;
    for (;;) {
      if (__atomic_load_n(&h->st[i], __ATOMIC_RELAXED) == KH_EMPTY) {
        if (!insert) return (size_t)-1;
        if ((h->n + 1) * 2 > h->mask + 1) break;     // >50% load: grow, then retry
        h->key[i] = key; h->st[i] = KH_WANT; h->n++;
        if (h->loc) h->loc[i] = 0;
        return i;
      }
      if (h->key[i] == key) return i;
      i = (i + 1) & h->mask;
    }
    kh_grow(h);
  }
}

static inline int kh_base(char c) {
  switch (c) { case 'A': case 'a': return 0; case 'C': case 'c': return 1;
               case 'G': case 'g': return 2; case 'T': case 't': return 3; }
  return -1;
}
// canonical k-mer (min of forward and reverse complement) as a 2-bit packing.
static inline uint64_t kh_canon(uint64_t f, uint64_t r) { return f < r ? f : r; }

// ---- closed syncmers, streamed one base at a time ---------------------------
// A guardian junction yields k-1 straddling k-mers, highly redundant: a read holding
// any ONE has crossed the junction. Keeping a fixed fraction (closed syncmers, density
// about 2/(k-s+1)) shrinks the index. The selector must be a function of the k-mer
// ALONE so the read side can reproduce it: keep the k-mer when the minimum of its k-s+1
// s-mers falls at either end. Two properties matter:
//
//  strand symmetry. Keys are CANONICAL, so a strand-dependent predicate would let walks
//    crossing a segment in opposite orientations select different k-mers under one key
//    and record only some colors. Hashing each s-mer canonically fixes that: reverse-
//    complementing a k-mer maps s-mer position p to k-s-p, so the end pair {0, k-s} maps
//    to itself and an RC-invariant s-mer hash makes the predicate RC-invariant. One
//    forward scan then suffices.
//
//  no junction goes dark. In any k-s+1 consecutive k-mers, let the minimum s-mer over
//    their combined range sit at m. If m <= k-s the k-mer starting at m has it as its
//    first s-mer, else the k-mer ending at m has it as its last; either way one is
//    selected. A junction offers k-1 >= k-s+1 of them, so every guardian keeps at least
//    one (kmer_insert audits the windows that come out empty anyway).
//
// The end test compares VALUES, not the argmin: "leftmost wins" would break symmetry when
// the minimum value occurs at position 0 and again inside. The deque is trimmed at the
// front only after the new entry goes in, so it transiently holds ns+1 entries; with
// ns = k-s+1 <= 30 (k <= 31, s >= 2) that fits in 32.
#define SYNC_RING 32

typedef struct {
  uint32_t s, ns;                 // s-mer length; window k-s+1 (ns == 0: filter off)
  uint32_t sh;                    // 2*(s-1)
  uint64_t smask;
  uint64_t sf, sr;                // rolling s-mer, forward and reverse complement
  int      run;                   // valid bases since the last non-ACGT
  int64_t  j;                     // index of the last s-mer pushed; -1 = none yet
  uint64_t hbuf[SYNC_RING];       // its hash, addressed by j & (SYNC_RING-1)
  uint64_t dqh[SYNC_RING];        // monotonic deque: the hashes...
  int64_t  dqj[SYNC_RING];        // ...and the s-mer index each came from
  uint32_t dqa, dqb;              // the deque is [dqa, dqb), both mod SYNC_RING
} sync_t;

static inline void sync_reset(sync_t *y) {
  y->sf = y->sr = 0; y->run = 0; y->j = -1; y->dqa = y->dqb = 0;
}

// s == 0 leaves the filter OFF and sync_push always selects, so callers keep one code
// path.
static inline void sync_init(sync_t *y, uint32_t k, uint32_t s) {
  y->s = s; y->ns = s ? k - s + 1 : 0;
  y->sh = s ? 2 * (s - 1) : 0;
  y->smask = s ? ((1ULL << (2 * s)) - 1) : 0;
  sync_reset(y);
}

// Feed one 2-bit base; returns whether the k-mer ENDING here is selected. The caller
// tracks whether k bases have been seen; that agrees with the window filling (ns s-mers
// = ns + s - 1 = k bases).
static inline int sync_push(sync_t *y, int b) {
  if (!y->ns) return 1;                               // filter off
  y->sf = ((y->sf << 2) | (uint64_t)b) & y->smask;
  y->sr = (y->sr >> 2) | ((uint64_t)(3 - b) << y->sh);
  // Saturate at s: `run` only answers "is the first s-mer complete".
  if (y->run < (int)y->s && ++y->run < (int)y->s) return 0;
  uint64_t h = kh_mix(y->sf < y->sr ? y->sf : y->sr);
  int64_t  j = ++y->j;
  y->hbuf[j & (SYNC_RING - 1)] = h;
  // Sliding-window minimum. Popping on >= keeps the newest of equal hashes, which
  // is safe only because the end test below compares VALUES and never positions.
  while (y->dqb != y->dqa && y->dqh[(y->dqb - 1) & (SYNC_RING - 1)] >= h) y->dqb--;
  y->dqh[y->dqb & (SYNC_RING - 1)] = h;
  y->dqj[y->dqb & (SYNC_RING - 1)] = j;
  y->dqb++;
  int64_t first = j - (int64_t)y->ns + 1;
  if (first < 0) return 0;                            // fewer than k bases so far
  while (y->dqj[y->dqa & (SYNC_RING - 1)] < first) y->dqa++;
  uint64_t mn = y->dqh[y->dqa & (SYNC_RING - 1)];
  return mn == h || mn == y->hbuf[first & (SYNC_RING - 1)];
}

// kroll_t: the rolling canonical k-mer every scan in this file runs -- 2-bit
// forward and reverse-complement packing, reset on non-ACGT, plus the syncmer
// predicate (s == 0: every k-mer selects).
typedef struct { uint64_t f, r, kmask; uint32_t k, sh; int run; sync_t y; } kroll_t;

static inline void kr_reset(kroll_t *q) { q->f = q->r = 0; q->run = 0; sync_reset(&q->y); }
static inline void kr_init(kroll_t *q, uint32_t k, uint32_t s) {
  q->k = k; q->kmask = (1ULL << (2 * k)) - 1; q->sh = 2 * (k - 1);
  sync_init(&q->y, k, s); kr_reset(q);
}
// Feed one base. -1: no complete k-mer ends here (non-ACGT, or fewer than k bases
// since the last one); 0: one does but the syncmer filter drops it; 1: selected.
// kr_key is then its canonical key.
static inline int kr_push(kroll_t *q, char c) {
  int b = kh_base(c);
  if (b < 0) { kr_reset(q); return -1; }
  q->f = ((q->f << 2) | (uint64_t)b) & q->kmask;
  q->r = (q->r >> 2) | ((uint64_t)(3 - b) << q->sh);
  int sel = sync_push(&q->y, b);
  if (q->run < (int)q->k) q->run++;                  // saturates: only "run >= k" is asked
  return q->run < (int)q->k ? -1 : sel;
}
static inline uint64_t kr_key(const kroll_t *q) { return kh_canon(q->f, q->r); }

// kstat_run: for the guardian set, how many window k-mers pin a single locus.
// Windows are DEDUPED first: a window shared by many colors is one index entry, so
// the rates describe what a k-mer index would actually hold. `absent` must come out
// 0 -- a window's k-mer occurs, by construction, in the walk the window came from.
typedef struct { size_t nwin, tot_km, uniq_km, win_uniq, absent; } kstat_t;

// kmer_insert: insert every guardian-window k-mer as KH_WANT; kmer_locate_color then
// streams the walks once to mark each KH_ONE (one locus) or KH_MANY (ambiguous). Shared
// by --kmer-stats and --classify. *nwin_out gets the distinct-window count. NULL if k is
// out of range.
// `all_kmers`: index EVERY k-mer of the walks instead of the guardian windows (the
// ablation against a complete colored index).
// `junc_only`: keep only k-mers STRADDLING the junction. One wholly in a flank belongs to
// the segment, not the edge, and pidx_build would discard it; dropping it here avoids
// inserting, rehashing and looking it up once per crossing walk, and makes the index
// cost independent of span (k-1 straddlers per window at span >= k-1). --kmer-stats
// passes 0: its question concerns the whole window.
static kh_t *kmer_insert(const bg_t *g, const gset_t *gs, uint32_t span, uint32_t k,
                         size_t *nwin_out, int all_kmers, int junc_only, uint32_t sk) {
  if (k < 2 || k > 31) { fprintf(stderr, "kmer: k must be in [2,31]\n"); return NULL; }
  if (sk && (sk < 2 || sk >= k)) {
    fprintf(stderr, "kmer: --syncmer s must be in [2,%u) for k=%u\n", k, k); return NULL;
  }
  char *buf = malloc(2 * (span ? span : 1));
  kroll_t q; kr_init(&q, k, sk);
  size_t ncand = 0, nkept = 0, ndark = 0;

  // Dedupe windows by content hash and collect their k-mers.
  kh_t *wseen = kh_new(gs->n ? gs->n : 1, 0);         // distinct window hashes
  kh_t *km    = kh_new(gs->n ? (size_t)gs->n * 4 : 16, 1);  // locates: needs loc
  size_t nwin = 0;
  if (all_kmers) {
    for (uint32_t w = 0; w < g->np; w++) {
      kr_reset(&q);
      for (uint32_t si = 0; si < g->pn[w]; si++) {
        uint32_t v = g->pv[w][si];
        const char *sq = nseq(g, v); size_t L = nsl(g, v);
        for (size_t o = 0; o < L; o++) {
          char c = (v & 1) ? (char)gfa_comp_table[(uint8_t)sq[L - 1 - o]] : sq[o];
          if (kr_push(&q, c) == 1) kh_slot(km, kr_key(&q), 1);
        }
      }
    }
  }
  for (uint32_t i = 0; i < gs->n && !all_kmers; i++) {
    const gpb_t *p = &gs->p[i];
    if (p->pid >= g->np || p->pos + 1 >= g->pn[p->pid]) continue;
    size_t lf = 0;                                    // junction sits at buf[lf]
    size_t wl = build_win(g, p->pid, p->pos, span, buf, &lf);
    // The dedup key is the WHOLE window even under junc_only, so the window count does not
    // depend on it.
    uint64_t wh = fnv1a(buf, wl, FNV_SEED);
    size_t ws = kh_slot(wseen, wh, 1);
    if (wseen->st[ws] != KH_WANT) continue;
    wseen->st[ws] = KH_ONE; nwin++;
    size_t cand = 0, kept = 0;
    kr_reset(&q);
    for (size_t j = 0; j < wl; j++) {
      int sel = kr_push(&q, buf[j]);
      if (sel < 0) continue;                          // a k-mer ends here, so j+1 >= k
      size_t start = j + 1 - k;
      if (junc_only && !(start + k > lf && start < lf)) continue;
      cand++;
      if (!sel) continue;
      kept++;
      kh_slot(km, kr_key(&q), 1);
    }
    ncand += cand; nkept += kept;
    if (cand && !kept) ndark++;
  }

  // A window offering all k-1 straddling k-mers keeps at least one (k-1 >= k-s+1), so
  // `dark` counts the windows that do not offer all of them: truncated at the end of a
  // walk, or with a broken ACGT run. Printed so the guarantee is audited on every build.
  if (sk)
    fprintf(stderr, "syncmer: s=%u  %zu/%zu window k-mers kept (%.4f), "
                    "%zu/%zu windows left with none (%.4f%%)\n",
            sk, nkept, ncand, ncand ? (double)nkept / (double)ncand : 0.0,
            ndark, nwin, nwin ? 100.0 * (double)ndark / (double)nwin : 0.0);

  free(buf); kh_free(wseen);
  if (nwin_out) *nwin_out = nwin;
  return km;
}

// Pass 2: stream every walk, keyed by graph position, and settle two things for each
// k-mer of the table:
//   locus -- KH_ONE if every occurrence is at the same graph position, else KH_MANY.
//   color -- which haplotypes contain it, read straight off the walks. Deriving it from
//            the guardian edge is unsound: two bubble branches agreeing over the k-mer
//            share one locus, so a read on the branch without a guardian would be
//            excluded by that edge's color set and an intersection would delete its own
//            color. From the walks, a read that is a substring of haplotype c and holds
//            k-mer K implies c is in K's mask.
// Colors are written unconditionally and the verdict is applied once at the end
// (pidx_build). `cmask` NULL collects no colors (--kmer-stats); otherwise color[w] is
// walk w's color index and nw the mask width.

// ---- parallel locate/color ----------------------------------------------------
// The passes below (walks in kmer_locate_color, segments in kmer_scan_segments) are
// read-mostly against `km`: every lookup is insert=0 and most positions miss the
// candidate table. Only a hit updates km->st/loc and cmask, so striped spinlocks
// serialise just that rare write. A global array is fine: one index is built at a time.
#define SLOCK_BITS 16
#define SLOCK_N (1u << SLOCK_BITS)
static unsigned char g_slock[SLOCK_N];

static inline void slot_lock(size_t sl) {
  uint32_t idx = (uint32_t)(kh_mix((uint64_t)sl) & (SLOCK_N - 1));
  for (unsigned char zero = 0;
       !__atomic_compare_exchange_n(&g_slock[idx], &zero, (unsigned char)1, 0,
                                    __ATOMIC_ACQUIRE, __ATOMIC_RELAXED); zero = 0) {}
}
static inline void slot_unlock(size_t sl) {
  uint32_t idx = (uint32_t)(kh_mix((uint64_t)sl) & (SLOCK_N - 1));
  __atomic_store_n(&g_slock[idx], (unsigned char)0, __ATOMIC_RELEASE);
}

// kh_mark: the KH_WANT/ONE/MANY transition, lock-free on the hot path. Once a k-mer's
// first locus is recorded, a later hit can only confirm it (nothing to write) or flip ONE
// to MANY, which is absorbing. So:
//   MANY: relaxed load and return.
//   ONE : `loc` is written once, before the release-store of ONE, so after an acquire-load
//         sees ONE it is race-free to read. A mismatched locus stores MANY with an atomic
//         store: idempotent and absorbing.
//   WANT: the only real race (two threads, two loci, one `loc` to claim); taken under the
//         striped lock with the state re-checked. Once per distinct k-mer.
// The verdict is order-independent (MANY iff the occurrences span >= 2 loci) and color
// masks are ORs, so the table and the .pex are byte-identical at every thread count.
// `st` uses the atomic builtins everywhere because kh_slot's lock-free probe reads the
// same byte concurrently.
// `mt` == 0 (single-threaded caller) skips the atomics: a locked RMW costs even without
// contention, and the segment pass does one per mask word. The test is a loop-invariant,
// well-predicted branch.
static inline void kh_mark(kh_t *km, size_t sl, uint64_t loc, int mt) {
  if (!mt) {
    uint8_t s0 = km->st[sl];
    if (s0 == KH_MANY) return;
    if (s0 == KH_WANT) { km->loc[sl] = loc; km->st[sl] = KH_ONE; }
    else if (km->loc[sl] != loc) km->st[sl] = KH_MANY;
    return;
  }
  uint8_t st = __atomic_load_n(&km->st[sl], __ATOMIC_ACQUIRE);
  if (st == KH_MANY) return;
  if (st == KH_ONE) {
    if (km->loc[sl] != loc)
      __atomic_store_n(&km->st[sl], (uint8_t)KH_MANY, __ATOMIC_RELAXED);
    return;
  }
  slot_lock(sl);
  st = __atomic_load_n(&km->st[sl], __ATOMIC_RELAXED);
  if (st == KH_WANT) {
    km->loc[sl] = loc;
    __atomic_store_n(&km->st[sl], (uint8_t)KH_ONE, __ATOMIC_RELEASE);
  } else if (st == KH_ONE && km->loc[sl] != loc) {
    __atomic_store_n(&km->st[sl], (uint8_t)KH_MANY, __ATOMIC_RELAXED);
  }
  slot_unlock(sl);
}

// Color bits are ORs: commutative and idempotent, so no lock. The common case is
// "already set": a relaxed load and no write, avoiding cache-line contention.
static inline void cm_or_atomic(uint64_t *w, uint64_t bits, int mt) {
  if (!mt) { *w |= bits; return; }                    // serial: see kh_mark
  if ((__atomic_load_n(w, __ATOMIC_RELAXED) & bits) != bits)
    __atomic_fetch_or(w, bits, __ATOMIC_RELAXED);
}

// par_for: fn(ctx, i) for every i < n. With nthreads >= 2 the work is claimed off
// a shared counter, the calling thread working too; below that it is the plain
// serial loop, zero threading overhead.
typedef void (*par_fn)(const void *ctx, uint32_t i);
typedef struct { par_fn fn; const void *ctx; uint32_t n; unsigned int next; } par_t;

static void *par_worker(void *arg) {
  par_t *p = arg;
  for (uint32_t i; (i = __atomic_fetch_add(&p->next, 1, __ATOMIC_RELAXED)) < p->n; )
    p->fn(p->ctx, i);
  return NULL;
}

static void par_for(uint32_t n, int nthreads, par_fn fn, const void *ctx) {
  if (nthreads < 2) { for (uint32_t i = 0; i < n; i++) fn(ctx, i); return; }
  par_t p = { fn, ctx, n, 0 };
  pthread_t *th = malloc((size_t)(nthreads - 1) * sizeof(pthread_t));
  int nth = 0;
  // Work is claimed from a shared counter, so a thread that cannot be created only
  // slows the loop down: the calling thread and the ones that started finish it.
  if (th) for (; nth < nthreads - 1; nth++)
    if (pthread_create(&th[nth], NULL, par_worker, &p)) break;
  par_worker(&p);
  for (int t = 0; t < nth; t++) pthread_join(th[t], NULL);
  free(th);
}

typedef struct {
  const bg_t *g; kh_t *km; uint32_t k, nw;
  uint64_t *cmask; const uint32_t *color; int skip_interior;
  uint32_t sk;                                        // --syncmer s (0 = every k-mer)
  int mt;                                             // >1 thread: use atomics
} walk_job_t;

// The physical place of one base read by a walk: its segment and FORWARD offset,
// so a base reads the same whichever strand the walk takes through it.
static inline uint64_t end_place(const bg_t *g, uint32_t v, uint32_t o) {
  uint32_t fo = (v & 1) ? (uint32_t)nsl(g, v) - 1 - o : o;
  return ((uint64_t)(v >> 1) << 32) | fo;
}

// kmer_locus: locus key of one k-mer occurrence, from the walk's vertex and offset at
// its first (v0,o0) and last (v1,o1) base. The key is the first base of the FORWARD
// reading: the first base when the walk enters on a '+' vertex, the last when on a '-'
// one, so a walk reading a place backward gets the same key (otherwise one place spanning
// several segments would give two loci and a spurious KH_MANY). A k-mer whose ends lie on
// vertices of opposite strand (an inversion breakpoint) has no forward reading; it takes
// the smaller end, also direction-free.
static inline uint64_t kmer_locus(const bg_t *g, uint32_t v0, uint32_t o0,
                                  uint32_t v1, uint32_t o1) {
  uint64_t e0 = end_place(g, v0, o0), e1 = end_place(g, v1, o1);
  if ((v0 & 1) == (v1 & 1)) return (v0 & 1) ? e1 : e0;
  return e0 < e1 ? e0 : e1;
}

static void walk_process_one(const void *ctx, uint32_t w) {
  const walk_job_t *j = ctx;
  const bg_t *g = j->g; kh_t *km = j->km;
  uint32_t k = j->k, nw = j->nw;
  uint32_t cw = 0; uint64_t bit = 0;
  if (j->cmask) { cw = j->color[w] >> 6; bit = 1ULL << (j->color[w] & 63); }
  uint32_t qv[32], qo[32]; int qh = 0;                // ring of the last k positions
  kroll_t q; kr_init(&q, k, j->sk);
  for (uint32_t s = 0; s < g->pn[w]; s++) {
    uint32_t v = g->pv[w][s];
    const char *sq = nseq(g, v);
    size_t L = nsl(g, v);
    for (size_t o = 0; o < L; o++) {
      char c = (v & 1) ? (char)gfa_comp_table[(uint8_t)sq[L - 1 - o]] : sq[o];
      qv[qh & 31] = v; qo[qh & 31] = (uint32_t)o; qh++;
      // The syncmer filter goes first, ahead of the interior skip and the probe: it rejects
      // ~87.5% of positions, and each rejection saves a hash probe into a table far larger
      // than cache.
      if (kr_push(&q, c) != 1) continue;              // no k-mer, or not a syncmer
      if (j->skip_interior && o >= k - 1) continue;   // deferred to the segment pass
      size_t sl = kh_slot(km, kr_key(&q), 0);
      if (sl == (size_t)-1) continue;                 // not a candidate k-mer
      // Direction-free locus: see kmer_locus.
      int z0 = (qh - (int)k) & 31, z1 = (qh - 1) & 31;
      uint64_t loc = kmer_locus(g, qv[z0], qo[z0], qv[z1], qo[z1]);
      if (j->cmask) cm_or_atomic(&j->cmask[(size_t)sl * nw + cw], bit, j->mt);
      kh_mark(km, sl, loc, j->mt);
    }
  }
}

// `skip_interior`: a k-mer wholly inside one segment has the same value and locus
// whichever walk visits it, so pidx_build settles those once per SEGMENT
// (kmer_scan_segments) and skips them here. `o >= k - 1` means this vertex alone supplied
// all k bases: a plain position count, orientation-independent and immune to a short cycle
// revisiting the vertex (unlike comparing start and end vertex ids). kstat_run passes 0:
// it has no segment pass.
// The table must NOT grow here: lookups use insert=0, so `cmask` (sized from km->mask
// by the caller) stays valid.
//
// Walks, not colors, are the unit of work: one color can span several walks and walk
// lengths vary widely, so claiming walk-by-walk off a shared counter balances load.
static void kmer_locate_color(const bg_t *g, kh_t *km, uint32_t k,
                              uint64_t *cmask, const uint32_t *color, uint32_t nw,
                              int skip_interior, int nthreads, uint32_t sk) {
  walk_job_t j = { g, km, k, nw, cmask, color, skip_interior, sk, nthreads >= 2 };
  par_for(g->np, nthreads, walk_process_one, &j);
}

// segcolor[seg] = union of colors of every walk traversing seg in either orientation,
// i.e. the color set of any interior k-mer of that segment. One pass over the walks'
// vertex lists rather than their bases.
static uint64_t *build_segcolor(const bg_t *g, const uint32_t *color, uint32_t nw) {
  uint32_t nseg = g->n / 2;
  uint64_t *segcolor = calloc((size_t)nseg * nw, sizeof(uint64_t));
  for (uint32_t w = 0; w < g->np; w++) {
    uint32_t cw = color[w] >> 6;
    uint64_t bit = 1ULL << (color[w] & 63);
    for (uint32_t s = 0; s < g->pn[w]; s++)
      segcolor[(size_t)(g->pv[w][s] >> 1) * nw + cw] |= bit;
  }
  return segcolor;
}

// The interior half of the locate+color pass, once per segment. A k-mer wholly inside a
// segment has the same value and locus (seg, offset) whichever walk and direction visits
// it, so scanning the stored sequence forward once covers every traversal. A hit still
// runs the full KH_WANT/ONE/MANY transition: an interior sequence equal to another
// candidate's value is a genuine ambiguity.
typedef struct {
  const bg_t *g; kh_t *km; uint32_t k, nw;
  uint64_t *cmask; const uint64_t *segcolor;
  uint32_t sk;
  int mt;
} seg_job_t;

static void seg_process_one(const void *ctx, uint32_t seg) {
  const seg_job_t *j = ctx;
  uint32_t k = j->k, nw = j->nw;
  const uint64_t *sc = j->segcolor + (size_t)seg * nw;
  if (!cm_any(sc, nw)) return;                        // no walk ever traverses it
  const char *sq = j->g->seq[seg];
  uint32_t L = j->g->slen[seg];
  if (!sq || L < k) return;                           // no interior k-mer fits
  kroll_t q; kr_init(&q, k, j->sk);
  for (uint32_t o = 0; o < L; o++) {
    if (kr_push(&q, sq[o]) != 1) continue;
    size_t sl = kh_slot(j->km, kr_key(&q), 0);
    if (sl == (size_t)-1) continue;
    // Same locus as walk_process_one: read forward inside one segment, the first base.
    uint64_t loc = ((uint64_t)seg << 32) | (uint32_t)(o + 1 - k);
    if (j->cmask) { uint64_t *dst = j->cmask + sl * nw;
                     for (uint32_t wi = 0; wi < nw; wi++)
                       if (sc[wi]) cm_or_atomic(&dst[wi], sc[wi], j->mt); }
    kh_mark(j->km, sl, loc, j->mt);
  }
}

// Segments, not colors, are the unit of work here: a segment carries no single
// color (many can share it) and its own length is what actually varies the work.
static void kmer_scan_segments(const bg_t *g, kh_t *km, uint32_t k, uint64_t *cmask,
                               const uint64_t *segcolor, uint32_t nw, int nthreads,
                               uint32_t sk) {
  seg_job_t j = { g, km, k, nw, cmask, segcolor, sk, nthreads >= 2 };
  par_for(g->n / 2, nthreads, seg_process_one, &j);
}

static kstat_t kstat_run(const bg_t *g, const gset_t *gs, uint32_t span, uint32_t k,
                         uint32_t sk) {
  kstat_t R = {0};
  size_t nwin = 0;
  // Both passes, no colors, no segment shortcut, one thread: it asks about WINDOWS, not
  // the index (a diagnostic, not a build).
  kh_t *km = kmer_insert(g, gs, span, k, &nwin, 0, 0, sk);
  if (!km) return R;
  kmer_locate_color(g, km, k, NULL, NULL, 0, 0, 1, sk);
  char *buf = malloc(2 * (span ? span : 1));
  kroll_t q; kr_init(&q, k, sk);

  // Re-walk the deduped windows and score them against the located k-mers.
  kh_t *wseen2 = kh_new(gs->n ? gs->n : 1, 0);
  size_t win_with_uniq = 0, tot_km = 0, uniq_km = 0, absent = 0;
  for (uint32_t i = 0; i < gs->n; i++) {
    const gpb_t *p = &gs->p[i];
    if (p->pid >= g->np || p->pos + 1 >= g->pn[p->pid]) continue;
    size_t wl = build_win(g, p->pid, p->pos, span, buf, NULL);
    uint64_t wh = fnv1a(buf, wl, FNV_SEED);
    size_t ws = kh_slot(wseen2, wh, 1);
    if (wseen2->st[ws] != KH_WANT) continue;
    wseen2->st[ws] = KH_ONE;
    int has_uniq = 0;
    kr_reset(&q);
    for (size_t j = 0; j < wl; j++) {
      if (kr_push(&q, buf[j]) != 1) continue;
      tot_km++;
      size_t sl = kh_slot(km, kr_key(&q), 0);
      uint8_t st = (sl == (size_t)-1) ? KH_WANT : km->st[sl];
      if (st == KH_ONE) { uniq_km++; has_uniq = 1; }
      else if (st == KH_WANT) absent++;
    }
    if (has_uniq) win_with_uniq++;
  }

  R.nwin = nwin; R.tot_km = tot_km; R.uniq_km = uniq_km;
  R.win_uniq = win_with_uniq; R.absent = absent;
  free(buf); kh_free(wseen2); kh_free(km);
  return R;
}

static void kstat(const bg_t *g, const gset_t *gs, uint32_t span, uint32_t k,
                  uint32_t sk, FILE *out) {
  kstat_t R = kstat_run(g, gs, span, k, sk);
  fprintf(out, "k=%u span=%u guardians=%u windows=%zu kmers=%zu\n",
          k, span, gs->n, R.nwin, R.tot_km);
  fprintf(out, "  single-locus kmers : %zu (%.1f%%)\n",
          R.uniq_km, R.tot_km ? 100.0 * (double)R.uniq_km / (double)R.tot_km : 0.0);
  fprintf(out, "  windows with >=1   : %zu (%.1f%%)  <- indexable without ambiguity\n",
          R.win_uniq, R.nwin ? 100.0 * (double)R.win_uniq / (double)R.nwin : 0.0);
  if (R.absent)
    fprintf(out, "  unlocated kmers    : %zu (should be 0; a window k-mer must occur "
                 "in the walk it came from)\n", R.absent);
}

// ---- read -> color classification against guardian k-mers -------------------
// The index holds single-locus guardian k-mers, each with the bitmask of colors whose
// guardian windows contain it (a window shared by 24 haplotypes cannot separate those 24
// but excludes every other color). A read votes +1 for each color in the mask of each of
// its k-mers found in the index; the winner is the arg-max. `disc` counts hits whose mask
// is a proper subset of all colors, i.e. those that discriminate.

// set-level strategies (--set-call); SET_ALL reports all three for comparison.
#define SET_ALL       0
#define SET_FRACTION  1
#define SET_INTERSECT 2
#define SET_VOTE      3
#define SET_SAMPLE    4
#define SET_MOSAIC    5      // paint the set with two recombining tracks (needs positions)
#define SET_COVER     6      // per window of the chromosome, the few colours whose union covers the hits

// pidx_t = the classifiable index: junction k-mers -> mask of colors containing
// them, plus the color names. Built from a graph (pidx_build) or loaded from a
// .pex file (pex_read); `run_classify` only ever sees this, so both paths share
// one code path. `own` marks names this struct must free (the .pex case).
typedef struct {
  kh_t     *km;
  uint64_t *cmask;                // nw words per slot: cmask + slot*nw
  const char **cname;
  uint32_t  nc, nw, k, span;
  uint32_t  sk;                   // --syncmer s; 0 = every junction k-mer (PEX2)
  size_t    nwin, nkmer;
  int       own;
  uint64_t  src_hash, src_size;   // identity of the GFA it was built from
  uint32_t  flags;                // PEX_F_*: how the masks were built (PEX4)
  char     *src_path;
  // Positions (--mosaic-ref, PEX5): the block of the reference walk each record sits in.
  uint16_t *blk;                  // per SLOT of km; BLK_NONE = no position. NULL = none stored
  uint32_t  blk_bp;               // block size in bp
  uint64_t  blk_len;              // end coordinate of the reference walk(s)
  char     *blk_ref;              // the colour whose walk(s) give the coordinates
} pidx_t;
#define PEX_F_FASTMASK 1u        // masks from pebble windows (--fast-mask), not the walks
#define BLK_NONE 0xFFFFu         // a record with no position (off the reference, no neighbour)
#define BLK_MAX  0xFFFEu         // the largest block id; a chromosome needs length/block < this

static void pidx_free(pidx_t *ix) {
  if (!ix) return;
  kh_free(ix->km); free(ix->cmask);
  if (ix->own && ix->cname) {
    for (uint32_t c = 0; c < ix->nc; c++) free((void *)ix->cname[c]);
  }
  free((void *)ix->cname); free(ix->src_path); free(ix->blk); free(ix->blk_ref); free(ix);
}

static int peb_valid(const bg_t *g, const gpb_t *p)
  { return p->pid < g->np && p->pos + 1 < g->pn[p->pid]; }

// peb_masks: the NON-SATURATED mask of every candidate key: the union of colors of every
// distinct guardian window that proposes it (a window's colors being every haplotype with
// a pebble on that exact content), ORed into `nsat` (nw words per slot of `km`). No walk
// is streamed. `norig`, if given, gets the exact number of distinct origin windows per
// key; nprop/nkeys/nmulti (each may be NULL) the proposal counts. Used by
// --pebble-ambiguity and --fast-mask.
static void peb_masks(const bg_t *g, const gset_t *gs, uint32_t span, uint32_t k,
                      uint32_t sk, int junc_only, kh_t *km, const uint32_t *color,
                      uint32_t nw, uint64_t *nsat, uint32_t *norig,
                      size_t *nprop, size_t *nkeys, size_t *nmulti) {
  char *buf = malloc(2 * (span ? span : 1));
  kroll_t q; kr_init(&q, k, sk);
  // Pass A: hash every pebble's window and grow the window table to its final size before
  // any buffer indexes it (kh_grow does not know about them). The hash is kept per pebble
  // for the next passes.
  uint64_t *ph = malloc((gs->n ? gs->n : 1) * sizeof(uint64_t));
  kh_t *wt = kh_new(gs->n ? gs->n : 1, 0);
  for (uint32_t i = 0; i < gs->n; i++) {
    const gpb_t *p = &gs->p[i];
    if (!peb_valid(g, p)) continue;
    size_t lf = 0;
    size_t wl = build_win(g, p->pid, p->pos, span, buf, &lf);
    uint64_t wh = fnv1a(buf, wl, FNV_SEED);
    ph[i] = wh;
    kh_slot(wt, wh, 1);
  }

  // Pass B: every pebble ORs its color into its window's mask. Must complete before any
  // key reads a mask, since pebbles sharing a content arrive in any order.
  uint64_t *wcol = calloc((wt->mask + 1) * (size_t)nw, sizeof(uint64_t));
  for (uint32_t i = 0; i < gs->n; i++) {
    const gpb_t *p = &gs->p[i];
    if (!peb_valid(g, p)) continue;
    size_t ws = kh_slot(wt, ph[i], 0);
    uint32_t c = color[p->pid];
    wcol[ws * nw + (c >> 6)] |= 1ULL << (c & 63);
  }

  // Pass C: candidates are extracted once per distinct window content, as in
  // kmer_insert, each ORing that window's final mask into its non-saturated mask. `lastw`
  // (1 + the last window that proposed the key) counts each (key, window) pair once even
  // when a key repeats inside a window; it is needed only for the origin counts, since
  // masks are ORs.
  uint32_t *lastw = norig ? calloc(km->mask + 1, sizeof(uint32_t)) : NULL;
  for (uint32_t i = 0; i < gs->n; i++) {
    const gpb_t *p = &gs->p[i];
    if (!peb_valid(g, p)) continue;
    size_t ws = kh_slot(wt, ph[i], 0);
    if (wt->st[ws] != KH_WANT) continue;                // candidates already pulled from this content
    wt->st[ws] = KH_ONE;
    size_t lf = 0;
    size_t wl = build_win(g, p->pid, p->pos, span, buf, &lf);
    const uint64_t *wc = wcol + ws * nw;
    kr_reset(&q);
    for (size_t j = 0; j < wl; j++) {
      int sel = kr_push(&q, buf[j]);
      if (sel < 0) continue;
      size_t start = j + 1 - k;
      if (junc_only && !(start + k > lf && start < lf)) continue;
      if (!sel) continue;
      size_t sl = kh_slot(km, kr_key(&q), 0);
      if (sl == (size_t)-1) continue;                   // not a candidate: km came from other windows
      for (uint32_t z = 0; z < nw; z++) nsat[sl * nw + z] |= wc[z];
      if (!norig) continue;
      if (nprop) (*nprop)++;
      if (lastw[sl] == (uint32_t)ws + 1) continue;
      lastw[sl] = (uint32_t)ws + 1;
      if (norig[sl] == 0) { if (nkeys) (*nkeys)++; }
      else if (norig[sl] == 1) { if (nmulti) (*nmulti)++; }  // crossing from 1 to 2 origins
      norig[sl]++;
    }
  }
  free(buf); free(ph); kh_free(wt); free(wcol); free(lastw);
}

// KH_MANY keys, by WHERE their occurrences sit; each is placed by kmer_locus (the same
// key kh_mark uses). Class 0 is a consistency check and must be empty.
//   0 strand    -- every occurrence is one physical place, yet KH_MANY: must be 0
//   1 segment   -- >=2 places, all inside one segment: a tandem repeat, or
//                  branches agreeing over the k-mer but ending apart in it
//   2 sibling   -- >=2 segments, every pair adjacent or sharing a neighbour, e.g.
//                  the branches of one bubble ending in the same bases
//   3 distant   -- up to PAMB_LCAP places, some pair not close in the graph
//   4 highcopy  -- more than PAMB_LCAP distinct places
#define PAMB_NC 5
#define PAMB_LCAP 8
static const char *pamb_cname[PAMB_NC] = { "strand", "segment", "sibling", "distant", "highcopy" };

// The physical-place pass over the walks, for KH_MANY keys only. `mid` maps a
// key slot to its KH_MANY id (UINT32_MAX otherwise); each id keeps up to
// PAMB_LCAP distinct places, and `ovf` records that there were more.
typedef struct {
  const bg_t *g; const kh_t *km; uint32_t k, sk;
  const uint32_t *mid; uint64_t *place; uint8_t *nplace, *ovf;
} ploc_t;

static void pamb_place_one(const void *ctx, uint32_t w) {
  const ploc_t *j = ctx;
  const bg_t *g = j->g; uint32_t k = j->k;
  uint32_t qv[32], qo[32]; int qh = 0;                // ring of the last k positions
  kroll_t q; kr_init(&q, k, j->sk);
  for (uint32_t s = 0; s < g->pn[w]; s++) {
    uint32_t v = g->pv[w][s];
    const char *sq = nseq(g, v);
    size_t L = nsl(g, v);
    for (size_t o = 0; o < L; o++) {
      char c = (v & 1) ? (char)gfa_comp_table[(uint8_t)sq[L - 1 - o]] : sq[o];
      qv[qh & 31] = v; qo[qh & 31] = (uint32_t)o; qh++;
      if (kr_push(&q, c) != 1) continue;
      size_t sl = kh_slot((kh_t *)j->km, kr_key(&q), 0);
      if (sl == (size_t)-1) continue;
      uint32_t id = j->mid[sl];
      if (id == UINT32_MAX) continue;
      // Exactly as walk_process_one keys the locus.
      int z0 = (qh - (int)k) & 31, z1 = (qh - 1) & 31;
      uint64_t pl = kmer_locus(g, qv[z0], qo[z0], qv[z1], qo[z1]);
      slot_lock(id);
      uint64_t *P = j->place + (size_t)id * PAMB_LCAP;
      uint8_t np = j->nplace[id], seen = 0;
      for (uint8_t x = 0; x < np; x++) if (P[x] == pl) { seen = 1; break; }
      if (!seen) {
        if (np < PAMB_LCAP) { P[np] = pl; j->nplace[id] = np + 1; }
        else j->ovf[id] = 1;
      }
      slot_unlock(id);
    }
  }
}

// Two segments are close if adjacent, or if they share an out-neighbour (both
// orientations of each). g->adj holds only the arcs the walks actually take, so
// adjacency is directed: seen from a but maybe not from b. seg_close tests it
// one way; seg_near asks both ways, so a key's class never depends on the order
// its places were recorded in (threads make that order vary).
static int seg_close(const bg_t *g, uint32_t a, uint32_t b) {
  for (uint32_t va = 2 * a; va <= 2 * a + 1; va++)
    for (uint32_t x = 0; x < g->adjn[va]; x++) {
      uint32_t na = g->adj[va][x] >> 1;
      if (na == b) return 1;
      for (uint32_t vb = 2 * b; vb <= 2 * b + 1; vb++)
        for (uint32_t y = 0; y < g->adjn[vb]; y++)
          if ((g->adj[vb][y] >> 1) == na) return 1;
    }
  return 0;
}
static int seg_near(const bg_t *g, uint32_t a, uint32_t b)
  { return seg_close(g, a, b) || seg_close(g, b, a); }

static int pamb_class(const bg_t *g, const uint64_t *P, uint8_t np, uint8_t ovf) {
  if (ovf) return 4;
  if (np <= 1) return 0;
  int one_seg = 1;
  for (uint8_t x = 1; x < np; x++) if ((P[x] >> 32) != (P[0] >> 32)) { one_seg = 0; break; }
  if (one_seg) return 1;
  for (uint8_t x = 0; x < np; x++)
    for (uint8_t y = x + 1; y < np; y++) {
      uint32_t a = (uint32_t)(P[x] >> 32), b = (uint32_t)(P[y] >> 32);
      if (a != b && !seg_near(g, a, b)) return 3;
    }
  return 2;
}


// many_classes: the locus class (an index into pamb_cname) of every KH_MANY slot
// of `km`, from one more stream over the walks restricted to those keys; 0xff
// for every other slot. Shared by --pebble-ambiguity, which reports the classes,
// and by an index build with --keep-many, which keeps some of them.
static uint8_t *many_classes(const bg_t *g, const kh_t *km, uint32_t k, uint32_t sk,
                             int nthreads) {
  size_t nks = km->mask + 1;
  uint32_t *mid = malloc(nks * sizeof(uint32_t));
  uint32_t nmany = 0;
  for (size_t s = 0; s < nks; s++) mid[s] = km->st[s] == KH_MANY ? nmany++ : UINT32_MAX;
  ploc_t pj = { g, km, k, sk, mid,
                calloc((size_t)(nmany ? nmany : 1) * PAMB_LCAP, sizeof(uint64_t)),
                calloc(nmany ? nmany : 1, 1), calloc(nmany ? nmany : 1, 1) };
  if (nmany) par_for(g->np, nthreads, pamb_place_one, &pj);
  uint8_t *cls = malloc(nks ? nks : 1);
  memset(cls, 0xff, nks ? nks : 1);
  for (size_t s = 0; s < nks; s++)
    if (mid[s] != UINT32_MAX)
      cls[s] = (uint8_t)pamb_class(g, pj.place + (size_t)mid[s] * PAMB_LCAP,
                                   pj.nplace[mid[s]], pj.ovf[mid[s]]);
  free(mid); free(pj.place); free(pj.nplace); free(pj.ovf);
  return cls;
}

// --keep-many: one bit per class; --keep-ambiguous sets them all.
#define KEEP_ALL ((1u << PAMB_NC) - 1)
static int parse_keep_many(const char *v, uint32_t *out) {
  char buf[256]; snprintf(buf, sizeof buf, "%s", v);
  for (char *t = strtok(buf, ","); t; t = strtok(NULL, ",")) {
    if (!strcmp(t, "all")) { *out = KEEP_ALL; continue; }
    int c = 0; while (c < PAMB_NC && strcmp(t, pamb_cname[c])) c++;
    if (c == PAMB_NC) {
      fprintf(stderr, "--keep-many: unknown class '%s' (strand, segment, sibling,"
                      " distant, highcopy, all)\n", t);
      return 1;
    }
    *out |= 1u << c;
  }
  return 0;
}

// `fast_mask` (--fast-mask): every candidate is indexed with its non-saturated mask
// (peb_masks) and nothing else: no segment or walk stream, hence no KH_ONE/KH_MANY
// verdict. Not sound: a read can lose its own color when one of its k-mers' masks misses
// it, so queries default to a tolerant intersection (--tol, see run_classify).
static pidx_t *pidx_build(const bg_t *g, const gset_t *gs, uint32_t span, uint32_t k,
                          uint32_t keep_many, int fast_mask, int all_kmers, int nthreads,
                          uint32_t sk) {
  const char **cptr; uint32_t nc;
  uint32_t *color = walk_colors(g, &cptr, &nc);       // same rule as gpeb()
  uint32_t nw = cm_nw(nc);

  size_t nwin = 0;
  double t_a = now_s();
  kh_t *km = kmer_insert(g, gs, span, k, &nwin, all_kmers, !all_kmers, sk);
  double t_ins = now_s() - t_a;
  if (!km) { free(color); free(cptr); return NULL; }

  // Colors are collected during the locate pass; the KH_ONE/KH_MANY verdict is applied at
  // the end. cmask is sized after the insert pass, when the table stops growing; neither
  // locate pass inserts.
  uint64_t *cmask = calloc((size_t)(km->mask + 1) * nw, sizeof(uint64_t));

  if (fast_mask) {
    t_a = now_s();
    peb_masks(g, gs, span, k, sk, !all_kmers, km, color, nw, cmask, NULL, NULL, NULL, NULL);
    double t_peb = now_s() - t_a;
    size_t kept = 0;
    for (size_t s = 0; s <= km->mask; s++)
      if (km->st[s] != KH_EMPTY && cm_any(cmask + s * nw, nw)) { km->st[s] = KH_IDX; kept++; }
    fprintf(stderr, "index: k=%u span=%u s=%u colors=%u windows=%zu index=%zu kmers"
                    " (--fast-mask: pebble-window masks, no locus check)\n",
            k, span, sk, nc, nwin, kept);
    fprintf(stderr, "phases: insert %.2fs  peb-masks %.2fs  (nw=%u word%s/mask)\n",
            t_ins, t_peb, nw, nw == 1 ? "" : "s");
    pidx_t *ix = calloc(1, sizeof(pidx_t));
    ix->km = km; ix->cmask = cmask; ix->cname = cptr; ix->nc = nc; ix->nw = nw;
    ix->k = k; ix->span = span; ix->sk = sk; ix->nwin = nwin; ix->nkmer = kept; ix->own = 0;
    ix->flags = PEX_F_FASTMASK;
    free(color);
    return ix;
  }

  // Interior k-mers first, once per segment; the walk stream then fills in only k-mers
  // that straddle two segments. Every candidate is still visited once either way.
  uint64_t *segcolor = build_segcolor(g, color, nw);
  t_a = now_s();
  kmer_scan_segments(g, km, k, cmask, segcolor, nw, nthreads, sk);
  double t_seg = now_s() - t_a; t_a = now_s();
  free(segcolor);
  kmer_locate_color(g, km, k, cmask, color, nw, 1, nthreads, sk);
  double t_walk = now_s() - t_a; t_a = now_s();

  // KH_IDX marks the slots the index holds. The verdict must test the state: every k-mer
  // the walks contain now carries colors, ambiguous ones included, so cm_any(mask) alone
  // would keep them all. The mask test is a second condition (a k-mer no color contains
  // indexes nothing). Masks on dropped slots are never read: pex_write, run_classify and
  // panel_nn gate on KH_IDX.
  // `keep_many` (--keep-many / --keep-ambiguous) also keeps the KH_MANY keys of the given
  // locus classes. Sound because a KH_MANY key's mask is read off the walks and holds every
  // haplotype carrying it at ANY locus. Keeping only some classes needs each key's class
  // (many_classes); keeping all does not.
  uint8_t *cls = (keep_many && keep_many != KEEP_ALL) ? many_classes(g, km, k, sk, nthreads) : NULL;
  size_t kept = 0, dropped = 0, kept_many[PAMB_NC] = {0}, kept_many_all = 0;
  for (size_t s = 0; s <= km->mask; s++)
    if (km->st[s] != KH_EMPTY) {
      int many = km->st[s] == KH_MANY;
      int keep = km->st[s] == KH_ONE ||
                 (many && (keep_many == KEEP_ALL || (cls && ((keep_many >> cls[s]) & 1))));
      if (keep && cm_any(cmask + s * nw, nw)) {
        if (many) { kept_many_all++; if (cls) kept_many[cls[s]]++; }
        km->st[s] = KH_IDX; kept++;
      } else dropped++;
    }
  free(cls);
  if (keep_many) {
    fprintf(stderr, "keep-many: %zu KH_MANY keys indexed", kept_many_all);
    if (keep_many == KEEP_ALL) fprintf(stderr, " (all classes)");
    else for (int c = 0; c < PAMB_NC; c++)
      if ((keep_many >> c) & 1) fprintf(stderr, " %s=%zu", pamb_cname[c], kept_many[c]);
    fputc('\n', stderr);
  }
  double t_fin = now_s() - t_a;
  fprintf(stderr, "index: k=%u span=%u s=%u colors=%u windows=%zu index=%zu kmers"
                  " (%zu ambiguous dropped)\n", k, span, sk, nc, nwin, kept, dropped);
  // Per-phase wall clock (the segment and walk passes scale differently).
  fprintf(stderr, "phases: insert %.2fs  segment-scan %.2fs  walk-stream %.2fs"
                  "  finalize %.2fs  (nw=%u word%s/mask)\n",
          t_ins, t_seg, t_walk, t_fin, nw, nw == 1 ? "" : "s");

  pidx_t *ix = calloc(1, sizeof(pidx_t));
  ix->km = km; ix->cmask = cmask; ix->cname = cptr; ix->nc = nc; ix->nw = nw;
  ix->k = k; ix->span = span; ix->sk = sk; ix->nwin = nwin; ix->nkmer = kept; ix->own = 0;
  free(color);
  return ix;
}

// ---- pebble ambiguity: windows proposing the same candidate k-mer ----------
// Two colour masks per key: saturated (`ix->cmask`, read off the walks: every haplotype
// whose walk contains the key; what gpeb indexes) and non-saturated (union of the colours
// of every distinct window proposing it). The latter must be a subset of the former;
// `nbad` counts violations and must be 0.
// `ix->cmask` also holds the saturated mask at KH_MANY slots, but only when the index is
// built WITHOUT --keep-ambiguous (otherwise KH_MANY is promoted to KH_IDX).
// Windows are grouped by CONTENT hash as in kmer_insert, so distinct loci sharing window
// text collapse into one "origin".
// Bins of distinct origin windows per key: 1, 2, 3-5, 6-9, >=10.
#define PAMB_NB 5
static const char *pamb_bname[PAMB_NB] = { "1", "2", "3-5", "6-9", ">=10" };
static int pamb_bin(uint32_t no) { return no == 1 ? 0 : no == 2 ? 1 : no <= 5 ? 2 : no <= 9 ? 3 : 4; }

typedef struct {
  size_t nprop, nkeys, nmulti;          // proposals, distinct keys, keys from >=2 windows
  size_t nbad;                          // non-saturated NOT a subset of saturated: must be 0
  // [state][bin]: state 0 = single-locus (KH_IDX, indexed), 1 = KH_MANY (dropped)
  size_t n[2][PAMB_NB], diff[2][PAMB_NB];
  uint64_t miss_sum[2][PAMB_NB];        // colors in saturated, absent from non-saturated
  uint32_t miss_max[2][PAMB_NB];
  uint32_t omax[2];
  size_t cn[PAMB_NC][PAMB_NB], cdiff[PAMB_NC][PAMB_NB];   // KH_MANY by locus class
} pamb_t;

static pamb_t pebble_ambiguity_run(const bg_t *g, const gset_t *gs, const pidx_t *ix,
                                    int junc_only, int nthreads) {
  pamb_t R; memset(&R, 0, sizeof R);
  uint32_t k = ix->k, span = ix->span, sk = ix->sk, nw = ix->nw;
  uint32_t nc2;
  uint32_t *color = walk_colors(g, NULL, &nc2);

  size_t nks = ix->km->mask + 1;                        // key-table slots
  uint64_t *nsat = calloc(nks * (size_t)nw, sizeof(uint64_t));
  uint32_t *norig = calloc(nks, sizeof(uint32_t));
  peb_masks(g, gs, span, k, sk, junc_only, ix->km, color, nw, nsat, norig,
            &R.nprop, &R.nkeys, &R.nmulti);

  // Pass D: the locus class of every KH_MANY key.
  uint8_t *cls = many_classes(g, ix->km, k, sk, nthreads);

  // Report, keyed on the real build's own verdict per slot.
  for (size_t s = 0; s < nks; s++) {
    uint32_t no = norig[s];
    if (!no) continue;                                  // never a candidate here
    uint8_t st = ix->km->st[s];
    int si = st == KH_MANY ? 1 : st == KH_IDX ? 0 : -1;
    if (si < 0) continue;                               // no walk carries it: nothing to compare
    int b = pamb_bin(no);
    const uint64_t *sat = ix->cmask + s * nw, *ns = nsat + s * nw;
    R.n[si][b]++;
    if (no > R.omax[si]) R.omax[si] = no;
    uint32_t miss = 0; int bad = 0;
    for (uint32_t z = 0; z < nw; z++) {
      miss += (uint32_t)__builtin_popcountll(sat[z] & ~ns[z]);
      if (ns[z] & ~sat[z]) bad = 1;
    }
    if (bad) R.nbad++;
    if (miss) {
      R.diff[si][b]++;
      R.miss_sum[si][b] += miss;
      if (miss > R.miss_max[si][b]) R.miss_max[si][b] = miss;
    }
    if (si == 1) {
      int c = cls[s];
      R.cn[c][b]++;
      if (miss) R.cdiff[c][b]++;
    }
  }

  free(color);
  free(nsat); free(norig); free(cls);
  return R;
}

static void pebble_ambiguity(const bg_t *g, const gset_t *gs, const pidx_t *ix,
                              int junc_only, int nthreads, FILE *out) {
  pamb_t R = pebble_ambiguity_run(g, gs, ix, junc_only, nthreads);
  fprintf(out, "pebble-ambiguity: k=%u span=%u s=%u\n", ix->k, ix->span, ix->sk);
  fprintf(out, "  candidate k-mers: %zu distinct keys touched, %zu (key,window) proposals\n",
          R.nkeys, R.nprop);
  fprintf(out, "  proposed by exactly 1 window: %zu (%.4f%%)\n",
          R.nkeys - R.nmulti, R.nkeys ? 100.0 * (double)(R.nkeys - R.nmulti) / (double)R.nkeys : 0.0);
  fprintf(out, "  proposed by >=2 distinct windows: %zu (%.4f%%)\n",
          R.nmulti, R.nkeys ? 100.0 * (double)R.nmulti / (double)R.nkeys : 0.0);
  fprintf(out, "  saturated = mask read off the walks (what gpeb indexes); non-saturated ="
               " union of the proposing windows' colors\n");
  fprintf(out, "  non-saturated not a subset of saturated: %zu keys (must be 0)\n", R.nbad);
  static const char *sname[2] = { "KH_ONE ", "KH_MANY" };
  static const char *oname[2] = { "1 window   ", ">=2 windows" };
  // 1 vs >=2 split: bin 0 against the sum of the others.
  for (int si = 0; si < 2; si++)
    for (int mi = 0; mi < 2; mi++) {
      size_t n = 0, d = 0; uint64_t ms = 0; uint32_t mx = 0;
      for (int b = mi ? 1 : 0; b < (mi ? PAMB_NB : 1); b++) {
        n += R.n[si][b]; d += R.diff[si][b]; ms += R.miss_sum[si][b];
        if (R.miss_max[si][b] > mx) mx = R.miss_max[si][b];
      }
      fprintf(out, "  %s %s: n=%zu differ=%zu (%.4f%%) missing colors per differing key"
                   " mean %.3f max %u\n", sname[si], oname[mi], n, d,
              n ? 100.0 * (double)d / (double)n : 0.0,
              d ? (double)ms / (double)d : 0.0, mx);
    }
  for (int si = 0; si < 2; si++) {
    fprintf(out, "  origin windows per key, %s (max %u):\n", sname[si], R.omax[si]);
    for (int b = 0; b < PAMB_NB; b++) {
      size_t n = R.n[si][b], d = R.diff[si][b];
      fprintf(out, "    %s windows=%-4s n=%zu differ=%zu (%.4f%%) missing colors per"
                   " differing key mean %.3f max %u\n", sname[si], pamb_bname[b], n, d,
              n ? 100.0 * (double)d / (double)n : 0.0,
              d ? (double)R.miss_sum[si][b] / (double)d : 0.0, R.miss_max[si][b]);
    }
  }
  fprintf(out, "  KH_MANY by where the occurrences sit (place = the locus key,"
               " kmer_locus; up to %d places tracked):\n", PAMB_LCAP);
  for (int c = 0; c < PAMB_NC; c++) {
    size_t n = 0, d = 0;
    for (int b = 0; b < PAMB_NB; b++) { n += R.cn[c][b]; d += R.cdiff[c][b]; }
    fprintf(out, "    KH_MANY class=%-8s n=%zu differ=%zu (%.4f%%)  by windows", pamb_cname[c],
            n, d, n ? 100.0 * (double)d / (double)n : 0.0);
    for (int b = 0; b < PAMB_NB; b++) fprintf(out, " %s=%zu/%zu", pamb_bname[b], R.cn[c][b], R.cdiff[c][b]);
    fputc('\n', out);
  }
  fprintf(out, "  (a single-locus key can still differ: a branch with no guardian of its own\n"
               "   contributes no window, but its walk reaches the locus all the same)\n");
}

// ---- .pex : the built index, serialized -------------------------------------
// Caching the index avoids the walk streaming that dominates a classify run. The file
// carries the source GFA's identity so a stale index is refused:
//   "PEX3" | src_hash | src_size | (len,src_path) | k | span | s | nc
//          | (len,bytes)*nc | n | (key,mask)*n
// `s` is the syncmer parameter. PEX2 lacks it and reads back as s == 0 (every junction
// k-mer indexed); only the writer moved to PEX3.
#define PEX_MAGIC  "PEX3"
#define PEX_MAGIC2 "PEX2"
#define PEX_MAGIC4 "PEX4"
#define PEX_MAGIC5 "PEX5"
// PEX4 = PEX3 plus a `flags` word after `s` (PEX_F_*). It is written ONLY when a
// flag is set, i.e. for a --fast-mask index, so every normal index stays PEX3
// byte for byte and older readers keep reading it; a reader that cannot tell a
// fast-mask index apart would query it with the wrong default tolerance.
// PEX5 = PEX4 plus POSITIONS (--mosaic-ref): after `flags`, blk_bp (u32), blk_len
// (u64) and the reference colour (len, bytes); then every record carries a uint16
// block id after its mask. A new magic rather than a flag bit, because an older
// reader ignores flag bits it does not know and would misparse the records
// instead of refusing the file.
// FNV-1a over the raw file bytes (the compressed bytes for a .gz): what is identified is
// the file, not its expansion.
static uint64_t file_hash(const char *path, uint64_t *size_out) {
  FILE *f = fopen(path, "rb");
  if (!f) return 0;
  uint64_t h = FNV_SEED, n = 0;
  unsigned char buf[1 << 16];
  size_t got;
  while ((got = fread(buf, 1, sizeof buf, f)) > 0) { n += got; h = fnv1a(buf, got, h); }
  fclose(f);
  if (size_out) *size_out = n;
  return h;
}

static int pex_write(const char *path, const pidx_t *ix, const char *src) {
  // Written next to the target and renamed over it, so a full disk or a crash never
  // leaves a truncated index under the final name.
  char tmp[4096];
  if (snprintf(tmp, sizeof tmp, "%s.tmp", path) >= (int)sizeof tmp) {
    fprintf(stderr, "pex: path too long: %s\n", path); return -1;
  }
  FILE *f = fopen(tmp, "wb");
  if (!f) { fprintf(stderr, "pex: cannot write %s\n", tmp); return -1; }
  // `src` NULL: the index was DERIVED from another one (subsampling), so the source
  // identity is inherited and the GFA need not be present.
  uint64_t hsz = 0, hh = 0;
  const char *sp = src;
  if (src) hh = file_hash(src, &hsz);
  else { hh = ix->src_hash; hsz = ix->src_size; sp = ix->src_path; }
  uint32_t sl = (uint32_t)(sp ? strlen(sp) : 0);
  int v5 = ix->blk != NULL;
  fwrite(v5 ? PEX_MAGIC5 : ix->flags ? PEX_MAGIC4 : PEX_MAGIC, 1, 4, f);
  fwrite(&hh, sizeof hh, 1, f);
  fwrite(&hsz, sizeof hsz, 1, f);
  fwrite(&sl, sizeof sl, 1, f);
  if (sl) fwrite(sp, 1, sl, f);
  fwrite(&ix->k, sizeof ix->k, 1, f);
  fwrite(&ix->span, sizeof ix->span, 1, f);
  fwrite(&ix->sk, sizeof ix->sk, 1, f);
  if (ix->flags || v5) fwrite(&ix->flags, sizeof ix->flags, 1, f);
  if (v5) {
    uint32_t rl = (uint32_t)(ix->blk_ref ? strlen(ix->blk_ref) : 0);
    fwrite(&ix->blk_bp, sizeof ix->blk_bp, 1, f);
    fwrite(&ix->blk_len, sizeof ix->blk_len, 1, f);
    fwrite(&rl, sizeof rl, 1, f);
    if (rl) fwrite(ix->blk_ref, 1, rl, f);
  }
  fwrite(&ix->nc, sizeof ix->nc, 1, f);
  for (uint32_t c = 0; c < ix->nc; c++) {
    const char *s = ix->cname[c] ? ix->cname[c] : "?";
    uint32_t L = (uint32_t)strlen(s);
    fwrite(&L, sizeof L, 1, f); fwrite(s, 1, L, f);
  }
  uint64_t n = 0;
  for (size_t s = 0; s <= ix->km->mask; s++)
    if (ix->km->st[s] == KH_IDX) n++;
  fwrite(&n, sizeof n, 1, f);
  for (size_t s = 0; s <= ix->km->mask; s++)
    if (ix->km->st[s] == KH_IDX) {
      fwrite(&ix->km->key[s], sizeof(uint64_t), 1, f);
      fwrite(ix->cmask + s * ix->nw, sizeof(uint64_t), ix->nw, f);
      if (v5) fwrite(&ix->blk[s], sizeof(uint16_t), 1, f);
    }
  // ferror is sticky: it reports any fwrite above that failed (e.g. ENOSPC); fflush and
  // fclose catch the data still buffered.
  int bad = ferror(f) | fflush(f);
  bad |= fclose(f);
  if (!bad && rename(tmp, path)) bad = 1;
  if (bad) {
    fprintf(stderr, "pex: failed writing %s (disk full or not writable?)\n", path);
    remove(tmp);
    return -1;
  }
  fprintf(stderr, "pex: wrote %s (%llu kmers, %u colors, k=%u span=%u s=%u,"
                  " %zu B/record, src=%s%s%s)\n",
          path, (unsigned long long)n, ix->nc, ix->k, ix->span, ix->sk,
          sizeof(uint64_t) * (1 + (size_t)ix->nw) + (v5 ? sizeof(uint16_t) : 0),
          sp ? sp : "?", v5 ? ", positions on " : "", v5 ? ix->blk_ref : "");
  return 0;
}

static pidx_t *pex_read(const char *path) {
  FILE *f = fopen(path, "rb");
  if (!f) { fprintf(stderr, "pex: cannot read %s\n", path); return NULL; }
  char magic[4];
  int v3 = 0, v4 = 0, v5 = 0;
  if (fread(magic, 1, 4, f) != 4 ||
      (!(v5 = !memcmp(magic, PEX_MAGIC5, 4)) && !(v4 = !memcmp(magic, PEX_MAGIC4, 4)) &&
       !(v3 = !memcmp(magic, PEX_MAGIC, 4)) && memcmp(magic, PEX_MAGIC2, 4))) {
      fprintf(stderr, "pex: %s is not a %s, %s, %s or %s file"
                      " (older indexes must be rebuilt)\n",
            path, PEX_MAGIC5, PEX_MAGIC4, PEX_MAGIC, PEX_MAGIC2); fclose(f); return NULL;
  }
  if (v5) v4 = 1;                                 // PEX5 is PEX4 plus positions
  if (v4) v3 = 1;                                 // PEX4 is PEX3 plus flags
  pidx_t *ix = calloc(1, sizeof(pidx_t)); ix->own = 1;
  uint32_t sl = 0;
  if (fread(&ix->src_hash, sizeof ix->src_hash, 1, f) != 1 ||
      fread(&ix->src_size, sizeof ix->src_size, 1, f) != 1 ||
      fread(&sl, sizeof sl, 1, f) != 1) goto bad;
  if (sl) {
    ix->src_path = malloc((size_t)sl + 1);
    if (fread(ix->src_path, 1, sl, f) != sl) goto bad;
    ix->src_path[sl] = 0;
  }
  if (fread(&ix->k, sizeof ix->k, 1, f) != 1 ||
      fread(&ix->span, sizeof ix->span, 1, f) != 1) goto bad;
  ix->sk = 0;                                     // PEX2: no filter, every k-mer
  if (v3 && fread(&ix->sk, sizeof ix->sk, 1, f) != 1) goto bad;
  if (v4 && fread(&ix->flags, sizeof ix->flags, 1, f) != 1) goto bad;
  if (v5) {
    uint32_t rl = 0;
    if (fread(&ix->blk_bp, sizeof ix->blk_bp, 1, f) != 1 ||
        fread(&ix->blk_len, sizeof ix->blk_len, 1, f) != 1 ||
        fread(&rl, sizeof rl, 1, f) != 1 || !ix->blk_bp) goto bad;
    ix->blk_ref = malloc((size_t)rl + 1);
    if (fread(ix->blk_ref, 1, rl, f) != rl) goto bad;
    ix->blk_ref[rl] = 0;
  }
  if (fread(&ix->nc, sizeof ix->nc, 1, f) != 1) goto bad;
  ix->cname = calloc(ix->nc ? ix->nc : 1, sizeof(char *));
  for (uint32_t c = 0; c < ix->nc; c++) {
    uint32_t L;
    if (fread(&L, sizeof L, 1, f) != 1) goto bad;
    char *s = malloc((size_t)L + 1);
    if (fread(s, 1, L, f) != L) { free(s); goto bad; }
    s[L] = 0; ix->cname[c] = s;
  }
  ix->nw = cm_nw(ix->nc);
  if (ix->sk && (ix->sk < 2 || ix->sk >= ix->k)) {
    fprintf(stderr, "pex: %s declares s=%u for k=%u, which is out of range\n",
            path, ix->sk, ix->k);
    goto bad;
  }
  uint64_t n;
  if (fread(&n, sizeof n, 1, f) != 1) goto bad;
  ix->km = kh_new(n ? (size_t)n : 1, 0);   // already located; loc would be dead
  ix->cmask = calloc((size_t)(ix->km->mask + 1) * ix->nw, sizeof(uint64_t));
  if (v5) {
    ix->blk = malloc((ix->km->mask + 1) * sizeof(uint16_t));
    for (size_t z = 0; z <= ix->km->mask; z++) ix->blk[z] = BLK_NONE;
  }
  for (uint64_t i = 0; i < n; i++) {
    uint64_t key;
    if (fread(&key, sizeof key, 1, f) != 1) goto bad;
    size_t sl = kh_slot(ix->km, key, 1);          // the table never grows past here
    if (fread(ix->cmask + sl * ix->nw, sizeof(uint64_t), ix->nw, f) != ix->nw) goto bad;
    if (v5 && fread(&ix->blk[sl], sizeof(uint16_t), 1, f) != 1) goto bad;
    ix->km->st[sl] = KH_IDX;
  }
  fclose(f);
  ix->nkmer = (size_t)n;
  fprintf(stderr, "pex: loaded %s (%llu kmers, %u colors, k=%u span=%u s=%u%s",
          path, (unsigned long long)n, ix->nc, ix->k, ix->span, ix->sk,
          ix->flags & PEX_F_FASTMASK ? ", --fast-mask" : "");
  if (v5) fprintf(stderr, ", positions on %s in %u bp blocks", ix->blk_ref, ix->blk_bp);
  fprintf(stderr, ")\n");
  return ix;
bad:
  fprintf(stderr, "pex: %s is truncated or corrupt\n", path);
  fclose(f); pidx_free(ix); return NULL;
}

// ---- --subsample : thin a built index, without rebuilding it ----------------
// Keeps a key iff kh_mix(key ^ seed) < threshold, so for one seed the keys kept at rate
// r/2 are a subset of those kept at r (nested ladder). The seed matters: kh_mix also
// places keys in kh_slot, so an unseeded test would correlate retention with table
// position; it also provides replicates.
// Records are marked KH_DROP rather than deleted: pex_write emits only KH_IDX slots, and
// the table keeps its shape.
// it afterwards.
static size_t pex_subsample(pidx_t *ix, double rate, uint64_t seed) {
  // rate >= 1 must keep everything EXACTLY: a UINT64_MAX threshold would still drop the
  // key hashing to the maximum.
  int keep_all = (rate >= 1.0);
  uint64_t thr = keep_all ? 0 : (uint64_t)(rate * 18446744073709551616.0);
  size_t in = 0, kept = 0;
  for (size_t s = 0; s <= ix->km->mask; s++) {
    if (ix->km->st[s] != KH_IDX) continue;
    in++;
    if (keep_all || kh_mix(ix->km->key[s] ^ seed) < thr) { kept++; continue; }
    ix->km->st[s] = KH_DROP;
  }
  ix->nkmer = kept;
  fprintf(stderr, "subsample: %zu -> %zu records (rate requested %.6g, achieved"
                  " %.6g, seed %llu)\n",
          in, kept, rate, in ? (double)kept / (double)in : 0.0,
          (unsigned long long)seed);
  return kept;
}

// ---- panel geometry ---------------------------------------------------------
// Jaccard similarity of each colour's nearest neighbour, measured on the index itself.
// The score a colour MISSING from the panel can reach rises with it: near-duplicate
// members leave no gap between a true call and an impostor. One pass over the index, no
// reads. Pairs are taken only over colours set in a mask, not over all |C|^2.
static void panel_nn(const pidx_t *ix, uint32_t *nn, double *jac) {
  uint32_t nc = ix->nc, nw = ix->nw;
  uint64_t *cnt = calloc(nc ? nc : 1, sizeof(uint64_t));
  uint64_t *both = calloc((size_t)nc * nc + 1, sizeof(uint64_t));
  uint32_t *bits = malloc((nc ? nc : 1) * sizeof(uint32_t));
  for (size_t s = 0; s <= ix->km->mask; s++) {
    if (ix->km->st[s] != KH_IDX) continue;
    const uint64_t *m = ix->cmask + s * nw;
    uint32_t nb = 0;
    for (uint32_t w = 0; w < nw; w++)
      for (uint64_t b = m[w]; b; b &= b - 1)
        bits[nb++] = (w << 6) + (uint32_t)__builtin_ctzll(b);
    for (uint32_t i = 0; i < nb; i++) {
      cnt[bits[i]]++;
      for (uint32_t j = i + 1; j < nb; j++) both[(size_t)bits[i] * nc + bits[j]]++;
    }
  }
  for (uint32_t a = 0; a < nc; a++) {
    nn[a] = a; jac[a] = 0.0;
    for (uint32_t b = 0; b < nc; b++) {
      if (b == a) continue;
      uint64_t in = a < b ? both[(size_t)a * nc + b] : both[(size_t)b * nc + a];
      uint64_t un = cnt[a] + cnt[b] - in;
      double j = un ? (double)in / (double)un : 0.0;
      if (j > jac[a]) { jac[a] = j; nn[a] = b; }
    }
  }
  free(cnt); free(both); free(bits);
}

// `topn`: ranks printed per set-call block; 0 = all, down to the first zero, so absence
// from the list means exactly zero.
// `tol` < 0 = not given: 1 (strict intersection) for a normal index, 0.9 for a
// --fast-mask one. Below 1 a read's set is every colour present in at least
// ceil(tol * hits) of its hit k-mers' masks (read off `votes`); at 1 it is the plain
// intersection and the code path is the untouched one.
// One whole line, EOL kept, into *line (grown as needed); 0 at the end of the file.
// gzgets alone stops at the buffer size, which would split a long-read line in two.
static int gz_line(gzFile fp, char **line, size_t *cap) {
  size_t len = 0;
  for (;;) {
    if (len + 2 > *cap) { *cap *= 2; *line = realloc(*line, *cap); if (!*line) { perror("realloc"); exit(1); } }
    if (!gzgets(fp, *line + len, (int)(*cap - len))) { if (!len) return 0; break; }
    len += strlen(*line + len);
    if (len && (*line)[len - 1] == '\n') break;
  }
  return 1;
}

static void run_classify(const pidx_t *ix, const char *reads, int setcall, int topn,
                         double tol, FILE *out) {
  if (tol < 0) tol = (ix->flags & PEX_F_FASTMASK) ? 0.9 : 1.0;
  if (tol < 1.0)
    fprintf(stderr, "classify: tolerant intersection --tol %.3g%s\n", tol,
            (ix->flags & PEX_F_FASTMASK) ? " (--fast-mask index)" : "");
  kh_t *km = ix->km;
  uint64_t *cmask = ix->cmask;
  const char **cptr = ix->cname;
  uint32_t nc = ix->nc, nw = ix->nw;
  // The query side applies the SAME syncmer predicate as the build: selection is a
  // property of the k-mer up to reverse complement, so a k-mer the index skipped can
  // never be in the table and need not be probed.
  kroll_t q; kr_init(&q, ix->k, ix->sk);

  // Colour -> sample. A colour is "SAMPLE#HAP" when the sample has several haplotypes in
  // the graph and plain "SAMPLE" otherwise (bg_build's usplit rule). A diploid read set
  // splits credit between its donor's haplotypes; summing afterwards would double-count
  // reads holding both, so the per-sample union must be taken per read, in the loop below.
  uint32_t ns = 0;
  uint32_t *csamp = malloc((nc ? nc : 1) * sizeof(uint32_t));
  const char **sname = malloc((nc ? nc : 1) * sizeof(char *));
  uint32_t *slen_ = malloc((nc ? nc : 1) * sizeof(uint32_t));
  for (uint32_t c = 0; c < nc; c++) {
    const char *nm = cptr[c] ? cptr[c] : "?";
    const char *h = strchr(nm, '#');
    uint32_t L = h ? (uint32_t)(h - nm) : (uint32_t)strlen(nm);
    uint32_t sx = 0;
    while (sx < ns && !(slen_[sx] == L && !strncmp(sname[sx], nm, L))) sx++;
    if (sx == ns) { sname[ns] = nm; slen_[ns] = L; ns++; }
    csamp[c] = sx;
  }
  uint32_t *scontains = calloc(ns ? ns : 1, sizeof(uint32_t));
  // Per read, HOW MANY of a sample's colours the intersection holds. Exactly one: the
  // read pins that haplotype (private); more: it identifies the donor but not the copy
  // (shared). Counted directly, since marginals determine it only for two haplotypes.
  uint32_t *scnt      = calloc(ns ? ns : 1, sizeof(uint32_t));
  uint32_t *slast     = calloc(ns ? ns : 1, sizeof(uint32_t));
  uint32_t *sshar     = calloc(ns ? ns : 1, sizeof(uint32_t));
  uint32_t *cpriv     = calloc(nc ? nc : 1, sizeof(uint32_t));

  gzFile fp = gzopen(reads, "r");
  if (!fp) { fprintf(stderr, "classify: cannot open %s\n", reads); return; }
  fprintf(out, "#read\ttop_colors\tn_tied\tvotes\trunner_up\thits\tdisc_hits"
                 "\tinter_n\tinter_colors\tn_masks\tall_votes\n");

  uint64_t *all = calloc(nw, sizeof(uint64_t));       // every color, for `disc`
  for (uint32_t c = 0; c < nc; c++) cm_set(all, c);
  size_t cap = 1 << 16; char *line = malloc(cap);
  char name[256] = "?", tie[4096], allv[8192], isec[4096]; int in_seq = 0, fastq = 0;
  int fq_hdr = 1;                                     // in FASTQ, next line is a header
  uint32_t *votes = calloc(nc, sizeof(uint32_t));
  uint32_t *contains = calloc(nc, sizeof(uint32_t));  // reads whose intersection has c
  uint64_t *vtotal = calloc(nc, sizeof(uint64_t));    // raw k-mer hits over the set
  uint64_t *ginter = calloc(nw, sizeof(uint64_t));    // AND of per-read intersections
  // Which index k-mers the read set hit (one byte per slot), to report how much of each
  // colour's signature was covered: separates "not enough data" from "not this colour".
  uint8_t *seen = calloc(km->mask + 1, 1);
  int ginit = 0;
  size_t nread = 0, nclass = 0, nohit = 0, nresolved = 0, nsig = 0;
  size_t hits = 0, disc = 0;
  // Voting asks which colour is named most often; the path-signature argument says the
  // right operator is INTERSECTION: the guardian edges of a colour intersect, over their
  // colour sets, to exactly that colour (given no path is a subpath of another). Both
  // are accumulated. `nmask` counts the DISTINCT colour sets a read saw (k-mers of one
  // window share a set).
  uint64_t *inter = calloc(nw, sizeof(uint64_t)); int have_inter = 0;
  uint64_t *masks = malloc(64 * (size_t)nw * sizeof(uint64_t)); int nmask = 0;

  // Flush the current read. Print EVERY colour at the maximum vote: picking one among
  // ties would be arbitrary, and callers need "resolved to one" vs "narrowed to a set".
  #define CLS_FLUSH() do { if (in_seq) { nread++; \
      uint32_t v1 = 0, v2 = 0, ntie = 0; \
      for (uint32_t c = 0; c < nc; c++) \
        if (votes[c] > v1) { v2 = v1; v1 = votes[c]; } \
        else if (votes[c] > v2) v2 = votes[c]; \
      if (v1) { nclass++; \
        if (tol < 1.0) { \
          double tx = tol * (double)hits; uint32_t thr = (uint32_t)tx; \
          if ((double)thr < tx - 1e-9) thr++;            /* ceil, without libm */ \
          memset(inter, 0, (size_t)nw * sizeof(uint64_t)); \
          for (uint32_t c = 0; c < nc; c++) if (votes[c] >= thr) cm_set(inter, c); } \
        size_t o = 0; \
        for (uint32_t c = 0; c < nc && o + 64 < sizeof tie; c++) \
          if (votes[c] == v1) { ntie++; \
            o += (size_t)snprintf(tie + o, sizeof tie - o, "%s%s", \
                                  o ? "," : "", cptr[c] ? cptr[c] : "?"); } \
        if (ntie == 1) nresolved++; \
        if (!ginit) { cm_cpy(ginter, inter, nw); ginit = 1; } \
        else cm_and(ginter, inter, nw); \
        for (uint32_t c = 0; c < nc; c++) vtotal[c] += votes[c]; \
        size_t q = 0; \
        for (uint32_t c = 0; c < nc && q + 64 < sizeof allv; c++) \
          if (votes[c]) q += (size_t)snprintf(allv + q, sizeof allv - q, "%s%s:%u", \
                                  q ? "," : "", cptr[c] ? cptr[c] : "?", votes[c]); \
        uint32_t ni = cm_pop(inter, nw); \
        size_t z = 0; \
        for (uint32_t c = 0; c < nc; c++) if (cm_get(inter, c)) { \
          contains[c]++; scnt[csamp[c]]++; slast[csamp[c]] = c; \
          if (z + 64 < sizeof isec) z += (size_t)snprintf(isec + z, sizeof isec - z, \
                                  "%s%s", z ? "," : "", cptr[c] ? cptr[c] : "?"); } \
        for (uint32_t sx = 0; sx < ns; sx++) \
          if (scnt[sx]) { scontains[sx]++; \
            if (scnt[sx] == 1) cpriv[slast[sx]]++; else sshar[sx]++; \
            scnt[sx] = 0; } \
        if (ni == 1) nsig++; \
        fprintf(out, "%s\t%s\t%u\t%u\t%u\t%zu\t%zu\t%u\t%s\t%d\t%s\n", \
                name, tie, ntie, v1, v2, hits, disc, \
                ni, ni ? isec : "-", nmask, q ? allv : "-"); } \
      else { nohit++; fprintf(out, "%s\t-\t0\t0\t0\t0\t0\t0\t-\t0\t-\n", name); } \
      memset(votes, 0, nc * sizeof(uint32_t)); hits = disc = 0; \
      memset(inter, 0, (size_t)nw * sizeof(uint64_t)); have_inter = 0; nmask = 0; } \
    kr_reset(&q); } while (0)

  while (gz_line(fp, &line, &cap)) {
    size_t L = strlen(line);
    while (L && (line[L-1] == '\n' || line[L-1] == '\r')) line[--L] = 0;
    if (L == 0) continue;
    // A quality line may itself start with '@', so in FASTQ a header is recognised by
    // POSITION in the record, not by its first byte.
    if (line[0] == '>' || (line[0] == '@' && (!fastq || fq_hdr))) {
      if (line[0] == '@') fastq = 1;
      fq_hdr = 0;
      CLS_FLUSH();
      size_t j = 0; while (j + 1 < sizeof name && line[j+1] && line[j+1] != ' '
                           && line[j+1] != '\t') { name[j] = line[j+1]; j++; }
      name[j] = 0; in_seq = 1;
      continue;
    }
    if (fastq && line[0] == '+') {                    // skip the quality line
      gz_line(fp, &line, &cap);           // discard the quality line, however long
      fq_hdr = 1;                                     // ...the next line opens a record
      continue;
    }
    if (!in_seq) continue;
    for (size_t j = 0; j < L; j++) {
      if (kr_push(&q, line[j]) != 1) continue;       // no k-mer, or the index cannot hold it
      size_t sl = kh_slot(km, kr_key(&q), 0);
      if (sl == (size_t)-1 || km->st[sl] != KH_IDX) continue;   // absent, or dropped
      seen[sl] = 1;
      const uint64_t *m = cmask + sl * nw;
      hits++;
      if (!cm_eq(m, all, nw)) disc++;
      // Set bits only. This runs per HIT, so scanning all nc colors would make the
      // panel size -- not the read length -- the cost of classifying a read.
      for (uint32_t z = 0; z < nw; z++)
        for (uint64_t b = m[z]; b; b &= b - 1)
          votes[(z << 6) + (uint32_t)__builtin_ctzll(b)]++;
      if (!have_inter) { cm_cpy(inter, m, nw); have_inter = 1; }
      else cm_and(inter, m, nw);
      int seen = 0;
      for (int z = 0; z < nmask; z++)
        if (cm_eq(masks + (size_t)z * nw, m, nw)) { seen = 1; break; }
      if (!seen && nmask < 64) { cm_cpy(masks + (size_t)nmask * nw, m, nw); nmask++; }
    }
  }
  CLS_FLUSH();
  #undef CLS_FLUSH
  fprintf(stderr, "classify: reads=%zu with-hits=%zu (%.1f%%) resolved-to-one=%zu"
                  " (%.1f%% of hits) no-hit=%zu\n",
          nread, nclass, nread ? 100.0 * (double)nclass / (double)nread : 0.0,
          nresolved, nclass ? 100.0 * (double)nresolved / (double)nclass : 0.0, nohit);
  fprintf(stderr, "classify: signature intersection singleton=%zu (%.1f%% of hits)\n",
          nsig, nclass ? 100.0 * (double)nsig / (double)nclass : 0.0);

  // Set-level call: one colour for the WHOLE read set (SET_ALL runs every strategy).
  //  intersect -- AND of the per-read intersections. Exact on clean reads, but ONE
  //    corrupted read empties it (a sequencing error turning one allele's junction k-mer
  //    into the other's), so it gets worse with more data.
  //  fraction  -- rank by the fraction of reads whose intersection contains the colour.
  //    Sound (the true colour ranks first); whether it sits far above the rest is a
  //    property of the PANEL (near-duplicate members), reported below as nearest-neighbour
  //    J. A corrupted read costs 1/n instead of everything. The margin separates
  //    in-panel (~0.9) from absent samples (~0) better than the score. The default.
  //  vote      -- rank by raw junction-k-mer hits per colour, no intersection.
  if (nclass) {
    uint8_t *shown = calloc(nc, 1);
    if (setcall == SET_ALL || setcall == SET_INTERSECT) {
      uint32_t ni = cm_pop(ginter, nw);
      fprintf(stderr, "set-call[intersect] %u color(s)%s", ni, ni ? ": " : " -- EMPTY "
              "(a corrupted read killed it; use fraction)");
      for (uint32_t c = 0; c < nc; c++)
        if (cm_get(ginter, c)) fprintf(stderr, "%s", cptr[c] ? cptr[c] : "?"),
                               fprintf(stderr, " ");
      fprintf(stderr, "\n");
    }
    // Per color: how many of the k-mers carrying it did the set cover?
    uint64_t *sig = calloc(nc ? nc : 1, sizeof(uint64_t));
    uint64_t *cov = calloc(nc ? nc : 1, sizeof(uint64_t));
    for (size_t s2 = 0; s2 <= km->mask; s2++) {
      if (km->st[s2] != KH_IDX) continue;
      const uint64_t *m = cmask + s2 * nw;
      for (uint32_t w = 0; w < nw; w++)
        for (uint64_t b = m[w]; b; b &= b - 1) {
          uint32_t c = (w << 6) + (uint32_t)__builtin_ctzll(b);
          sig[c]++; if (seen[s2]) cov[c]++;
        }
    }
    if (setcall == SET_ALL || setcall == SET_FRACTION) {
      memset(shown, 0, nc);
      for (uint32_t rank = 0; rank < (topn ? (uint32_t)topn : nc); rank++) {
        uint32_t best = nc, bv = 0;
        for (uint32_t c = 0; c < nc; c++)
          if (!shown[c] && contains[c] >= bv) { bv = contains[c]; best = c; }
        if (best == nc || bv == 0) break;
        shown[best] = 1;
        fprintf(stderr, "set-call[fraction]  %-16s %.3f  (%u/%zu reads)"
                        "  signature %llu/%llu (%.1f%%)%s\n",
                cptr[best] ? cptr[best] : "?", (double)bv / (double)nclass,
                bv, nclass, (unsigned long long)cov[best], (unsigned long long)sig[best],
                sig[best] ? 100.0 * (double)cov[best] / (double)sig[best] : 0.0,
                rank ? "" : "   <- call");
      }
    }
    // Per SAMPLE: share of reads whose intersection contains ANY haplotype of the donor.
    // [fraction] ranks assemblies, this ranks individuals; neither replaces the other.
    if (setcall == SET_ALL || setcall == SET_SAMPLE) {
      uint8_t *sshown = calloc(ns ? ns : 1, 1);
      for (uint32_t rank = 0; rank < (topn ? (uint32_t)topn : ns); rank++) {
        uint32_t best = ns, bv = 0;
        for (uint32_t sx = 0; sx < ns; sx++)
          if (!sshown[sx] && scontains[sx] >= bv) { bv = scontains[sx]; best = sx; }
        if (best == ns || bv == 0) break;
        sshown[best] = 1;
        uint32_t nh = 0;
        for (uint32_t c = 0; c < nc; c++) if (csamp[c] == best) nh++;
        fprintf(stderr, "set-call[sample]    %-16.*s %.3f  (%u/%zu reads, %u haplotype%s)%s\n",
                (int)slen_[best], sname[best], (double)bv / (double)nclass,
                bv, nclass, nh, nh == 1 ? "" : "s", rank ? "" : "   <- call");
        if (nh > 1) {                       // the private/shared split, per copy
          fprintf(stderr, "                      private");
          for (uint32_t c = 0; c < nc; c++)
            if (csamp[c] == best)
              fprintf(stderr, " %s %.3f", cptr[c] ? cptr[c] : "?",
                      (double)cpriv[c] / (double)nclass);
          fprintf(stderr, "   shared %.3f\n", (double)sshar[best] / (double)nclass);
        }
      }
      free(sshown);
    }
    if (setcall == SET_ALL || setcall == SET_VOTE) {
      memset(shown, 0, nc);
      uint64_t tot = 0;
      for (uint32_t c = 0; c < nc; c++) tot += vtotal[c];
      for (uint32_t rank = 0; rank < (topn ? (uint32_t)topn : nc); rank++) {
        uint32_t best = nc; uint64_t bv = 0;
        for (uint32_t c = 0; c < nc; c++)
          if (!shown[c] && vtotal[c] >= bv) { bv = vtotal[c]; best = c; }
        if (best == nc || bv == 0) break;
        shown[best] = 1;
        fprintf(stderr, "set-call[vote]      %-16s %llu hits (%.1f%% of all)%s\n",
                cptr[best] ? cptr[best] : "?", (unsigned long long)bv,
                tot ? 100.0 * (double)bv / (double)tot : 0.0, rank ? "" : "   <- call");
      }
    }
    // Qualify the call with the panel's geometry: no score can separate a colour from a
    // very close neighbour, so report J rather than let a caller invent a threshold.
    if (nc > 1) {
      uint32_t *nn = malloc(nc * sizeof(uint32_t));
      double *jc = malloc(nc * sizeof(double));
      panel_nn(ix, nn, jc);
      uint32_t worst = 0;
      for (uint32_t c = 1; c < nc; c++) if (jc[c] > jc[worst]) worst = c;
      uint32_t best = nc, bv = 0;                      // the color just called
      for (uint32_t c = 0; c < nc; c++) if (contains[c] >= bv) { bv = contains[c]; best = c; }
      if (best < nc)
        fprintf(stderr, "panel: %s's nearest neighbour is %s at J=%.3f\n",
                cptr[best] ? cptr[best] : "?", cptr[nn[best]] ? cptr[nn[best]] : "?", jc[best]);
      fprintf(stderr, "panel: closest pair anywhere is %s ~ %s at J=%.3f -- the score a color"
                      " NOT in the panel can reach rises with this, so read the ranking\n",
              cptr[worst] ? cptr[worst] : "?", cptr[nn[worst]] ? cptr[nn[worst]] : "?", jc[worst]);
      free(nn); free(jc);
    }
    free(shown); free(sig); free(cov);
  }
  free(csamp); free((void *)sname); free(slen_); free(scontains);
  free(scnt); free(slast); free(sshar); free(cpriv);
  free(line); free(votes); free(contains); free(vtotal);
  free(all); free(ginter); free(inter); free(masks); free(seen); gzclose(fp);
}

// ---- positions: a block id for every single-locus record (--mosaic-ref) -------
//
// A record's locus is ix->km->loc[slot] = (segment << 32) | forward offset of the first
// base of the k-mer's forward reading (kmer_locus). Projecting it onto a walk of ONE
// reference colour gives a coordinate, and coordinate / B a block. The projection is per
// SEGMENT, so one pass over the walks serves every record:
//   - a segment the reference walk(s) visit gets the coordinate of its first base
//     (PK_REF; PK_REFMUL if visited at several places -- the first visit is kept);
//   - any other segment (a bubble allele) inherits the coordinate of the junction it
//     hangs from: the END of the last reference segment before it along the first walk
//     that visits it, or the START of that segment when the walk runs against the
//     reference. This places it within the bubble's width;
//   - a segment no walk reaches from a reference segment has none (PK_NONE).
// Offsets inside a segment are added: a v1.1 segment reaches 3.2 Mbp.
#define PK_NONE   0
#define PK_REF    1
#define PK_INH    2
#define PK_REFMUL 3

typedef struct {
  uint32_t  nseg;
  uint32_t *coord;                     // UINT32_MAX = none; PK_REF: first base on the
                                       // reference; PK_INH: the junction coordinate
  uint8_t  *kind, *rev;                // rev: the reference's first visit runs on '-'
  uint64_t  len;                       // end coordinate of the reference walk(s)
  uint32_t  nrefw;                     // reference walks (W-lines of the colour)
  uint64_t  nref_steps, ninh_diff;     // ninh_diff: a later walk placed a segment more
                                       // than a block away from where the first did
} proj_t;

static void proj_free(proj_t *P) {
  if (!P) return;
  free(P->coord); free(P->kind); free(P->rev); free(P);
}

static proj_t *proj_build(const gfa_t *gfa, const bg_t *g, const char *refname, uint32_t B) {
  uint32_t nseg = gfa->n_seg;
  uint8_t *isref = calloc(g->np ? g->np : 1, 1);
  uint32_t nrefw = 0;
  for (uint32_t w = 0; w < g->np; w++)
    if (g->pname[w] && !strcmp(g->pname[w], refname)) { isref[w] = 1; nrefw++; }
  if (!nrefw) {
    fprintf(stderr, "--mosaic-ref: no walk of colour \"%s\" (the first colour here is \"%s\")\n",
            refname, g->ncpool ? g->cpool[0] : "?");
    free(isref); return NULL;
  }
  // The W-lines of the reference colour share ONE coordinate system only if they are
  // fragments of one contig; two contigs would overlap and mix their blocks.
  const char *ctg = NULL;
  for (uint32_t w = 0; w < g->np; w++) {
    if (!isref[w]) continue;
    const char *c = gfa->sseq[gfa->walk[g->psrc[w]].snid].name;
    if (ctg && strcmp(ctg, c)) {
      fprintf(stderr, "--mosaic-ref: colour \"%s\" has walks on several contigs (%s, %s); the"
                      " coordinates would overlap. Build from a one-chromosome graph\n",
              refname, ctg, c);
      free(isref); return NULL;
    }
    ctg = c;
  }
  proj_t *P = calloc(1, sizeof *P);
  P->nseg = nseg; P->nrefw = nrefw;
  P->coord = malloc((size_t)nseg * sizeof(uint32_t));
  memset(P->coord, 0xff, (size_t)nseg * sizeof(uint32_t));
  P->kind = calloc(nseg ? nseg : 1, 1); P->rev = calloc(nseg ? nseg : 1, 1);
  // Pass 1: the reference walk(s) give every segment they visit its coordinate. A W-line
  // carries its own start, so a reference in several fragments lands in ONE system.
  for (uint32_t w = 0; w < g->np; w++) {
    if (!isref[w]) continue;
    int64_t st = gfa->walk[g->psrc[w]].st;
    uint64_t pos = st > 0 ? (uint64_t)st : 0;
    for (uint32_t s = 0; s < g->pn[w]; s++) {
      uint32_t v = g->pv[w][s], seg = v >> 1;
      P->nref_steps++;
      if (P->kind[seg] == PK_NONE) { P->coord[seg] = (uint32_t)pos; P->kind[seg] = PK_REF; P->rev[seg] = v & 1; }
      else P->kind[seg] = PK_REFMUL;
      pos += g->slen[seg];
    }
    if (pos > P->len) P->len = pos;
  }
  // Pass 2: every other walk, in file order. The first value a segment is given is
  // kept; later walks only vote on whether it was stable.
  for (uint32_t w = 0; w < g->np; w++) {
    if (isref[w]) continue;
    int have = 0, same = 1; uint32_t aseg = 0;
    for (uint32_t s = 0; s < g->pn[w]; s++) {
      uint32_t v = g->pv[w][s], seg = v >> 1; uint8_t kd = P->kind[seg];
      if (kd == PK_REF || kd == PK_REFMUL) {
        have = 1; aseg = seg; same = ((v & 1) == P->rev[seg]);
        continue;
      }
      if (!have) continue;                       // before the walk's first reference step
      uint32_t val = same ? P->coord[aseg] + g->slen[aseg] : P->coord[aseg];
      if (kd == PK_NONE) { P->coord[seg] = val; P->kind[seg] = PK_INH; }
      else {
        uint32_t o = P->coord[seg], d = o > val ? o - val : val - o;
        if (d > B) P->ninh_diff++;
      }
    }
  }
  free(isref);
  return P;
}

// The coordinate of a locus (segment, forward offset), or UINT32_MAX.
static inline uint32_t proj_coord(const proj_t *P, const bg_t *g, uint64_t loc) {
  uint32_t seg = (uint32_t)(loc >> 32), fo = (uint32_t)loc;
  uint8_t kd = P->kind[seg];
  if (kd == PK_NONE) return UINT32_MAX;
  uint32_t c = P->coord[seg];
  if (kd == PK_REF || kd == PK_REFMUL) c += P->rev[seg] ? (g->slen[seg] - 1 - fo) : fo;
  return c;
}

// Give every record of a freshly built index its block. Only a SOUND single-locus
// index has a position to give: a k-mer seen at several loci has no one place, and a
// --fast-mask or --index-all index never located its k-mers.
static int pidx_set_blocks(pidx_t *ix, const gfa_t *gfa, const bg_t *g,
                           const char *refname, uint32_t B) {
  if (!ix->km->loc || (ix->flags & PEX_F_FASTMASK)) {
    fprintf(stderr, "--mosaic-ref: this index never located its k-mers\n"); return -1;
  }
  if (!B) { fprintf(stderr, "--mosaic-block must be >= 1\n"); return -1; }
  double t0 = now_s();
  proj_t *P = proj_build(gfa, g, refname, B);
  if (!P) return -1;
  size_t ns = ix->km->mask + 1, nrec = 0, nnone = 0, nref = 0, ninh = 0;
  uint16_t *blk = malloc(ns * sizeof(uint16_t));
  for (size_t s = 0; s < ns; s++) blk[s] = BLK_NONE;
  for (size_t s = 0; s < ns; s++) {
    if (ix->km->st[s] != KH_IDX) continue;
    nrec++;
    uint64_t loc = ix->km->loc[s];
    uint32_t c = proj_coord(P, g, loc);
    if (c == UINT32_MAX) { nnone++; continue; }
    uint64_t b = c / B;
    if (b > BLK_MAX) {
      fprintf(stderr, "--mosaic-ref: block %llu does not fit a uint16 id; use a larger"
                      " --mosaic-block than %u\n", (unsigned long long)b, B);
      free(blk); proj_free(P); return -1;
    }
    blk[s] = (uint16_t)b;
    if (P->kind[(uint32_t)(loc >> 32)] == PK_INH) ninh++; else nref++;
  }
  uint32_t nsr = 0, nsi = 0, nsn = 0;
  for (uint32_t z = 0; z < P->nseg; z++) {
    if (P->kind[z] == PK_NONE) nsn++; else if (P->kind[z] == PK_INH) nsi++; else nsr++;
  }
  fprintf(stderr, "positions: %u walk(s) of %s span %llu bp -> %llu blocks of %u bp; "
                  "segments: %u on the reference, %u inherit, %u unplaced\n",
          P->nrefw, refname, (unsigned long long)P->len,
          (unsigned long long)((P->len + B - 1) / B), B, nsr, nsi, nsn);
  fprintf(stderr, "positions: %zu records: %zu on reference segments, %zu inherited, "
                  "%zu without a position (%.1fs)\n", nrec, nref, ninh, nnone, now_s() - t0);
  if (P->ninh_diff)
    fprintf(stderr, "positions: %llu walk steps placed an inherited segment more than a block"
                    " from where the first walk did\n", (unsigned long long)P->ninh_diff);
  free(ix->blk); free(ix->blk_ref);
  ix->blk = blk; ix->blk_bp = B; ix->blk_len = P->len; ix->blk_ref = strdup(refname);
  proj_free(P);
  return 0;
}

// ---- mosaic: paint a read set with two haplotype tracks (--set-call mosaic) ----
//
//
// Every read says, for each colour c, how many of its index hits c does NOT have,
// m(c) = hits - votes(c). A diploid pair of colours explains it better or worse, and a
// Li-Stephens chain over blocks of the reference lets the best pair change along the
// chromosome.
//
// Model. An observation is a read, or a read PAIR with votes and hits of the two mates
// summed (one haplotype, one place). Copying colour c, it has probability
//     p(c) = eps + (1 - eps) exp(-lambda m(c))
// (eps caps the cost of a read no colour explains -- a chimera, a misplaced read;
// lambda < 1 because the hits of one variant are correlated); a diploid pair (a, b)
// emits (p(a) + p(b)) / 2. Per block step each copy stays with probability 1 - rho or
// jumps to a uniformly drawn colour. The transition factorises, so the Viterbi path
// over the C x C ORDERED pairs costs O(C^2) per block:
//   V'(a,b) = E(a,b) + max( V(a,b) + 2 log(1-rho),
//                           max_a' V(a',b) + log(1-rho) + log(rho/C),   copy 1 jumps
//                           max_b' V(a,b') + log(1-rho) + log(rho/C),   copy 2 jumps
//                           max V + 2 log(rho/C) )                      both jump
// An observation lives in the MEDIAN block of its index hits (robust to stray hits far
// away). Observations that no colour tells from another (every m equal) carry no
// information and are not stored; a block with none still propagates the chain.
#define NOCOL 0xFFFFu            // a track holds no colour here (cover painting)
#define MR_AUTO  0
#define MR_HIFI  1
#define MR_SHORT 2

typedef struct {
  double lambda, rho, eps;       // < 0: the preset of the read type
  int reads;                     // MR_*
  const char *mate;              // --mate: the second file of read pairs
  const char *out;               // --mosaic-out PREFIX; NULL = the painting goes to stdout
  const char *gfa_out;           // --mosaic-gfa: splice the painting into this graph
  int threads;                   // --threads: the cover's scan of the reads (the other modes are serial)
  const char *gfa_in;            // the GFA --mosaic-gfa splices into, when the caller has not parsed it (read after the painting, for its colours only)
  uint64_t gfa_hash;             // ...its FNV-1a, when the staleness check computed it (gfa_hashed): the .gfai is checked against it too
  int gfa_hashed;
  const char *exclude;           // --mosaic-exclude: colours the painting may not use
  const char *keep;              // --mosaic-keep: colours kept whole next to the tracks
  const char *panel;             // label for the painting's header
  // --set-call cover
  uint32_t cwin;                 // --cover-window: bp per window (0 = the whole chromosome)
  int      cauto;                // --cover-window auto: windows of a fixed amount of evidence, not of bp
  uint64_t chits;                // --cover-hits: the hits an automatic window holds, at least
  uint32_t cmax;                 // --cover-max: most colours picked per window
  double   cfrac;                // --cover-min-frac: stop below this share of the window's hits
  double   cgain;                // --cover-min-gain: ...or below this many hits
  int      cstray;               // --cover-stray: blocks from its read's median a hit may lie at
  int      crule;                // CR_GREEDY | CR_SUPPORT
  int      cweight;              // CW_COUNT | CW_PRESENCE
  const char *check, *check_mate;  // --mosaic-check: an independent read set to score the graph on
  int      temper;               // --mosaic-temper: TEMPER_AUTO | TEMPER_OFF | TEMPER_GLOBAL | TEMPER_BLOCK
  double   temper_n;             // --mosaic-temper-n: observations per block (<= 0: 58 per 10 kb)
} mosaic_opt_t;
#define CR_GREEDY   0            // each pick covers what the earlier ones left: the union objective
#define CR_SUPPORT  1            // CONTROL: the N colours with the most hits each, ignoring overlap
#define CW_COUNT    0            // a record weighs the times the reads hit it (capped)
#define CW_PRESENCE 1            // a record weighs 1: alleles, not depth
// --mosaic-temper. A block's evidence grows with its reads, the cost of a switch does
// not, so at high coverage a few discordant pairs buy a switch. Tempering multiplies a
// block's log-emission by w <= 1, bringing the evidence back to N observations per block
// (N = 58 per 10 kb, the value the preset was fitted at).
#define TEMPER_AUTO   0          // the default: global for the short preset, off for HiFi
#define TEMPER_OFF    1          // no tempering (the painting before tempering, byte for byte)
#define TEMPER_GLOBAL 2          // one w = min(1, N / median n_obs per block) for the whole painting
#define TEMPER_BLOCK  3          // w_t = min(1, N / n_obs of block t)
#define TEMPER_N_10KB 58.0       // N for 10 kb blocks; scaled with the block size

// Presets chosen on held-out reads; the score is flat in rho above 0.1.
static void mosaic_preset(int reads, double *lam, double *rho, double *eps) {
  if (reads == MR_SHORT) { *lam = 4.0; *rho = 0.1; *eps = 0.01; }
  else                   { *lam = 0.5; *rho = 0.1; *eps = 0.001; }
}

// FASTA or FASTQ, plain or gzipped, one record at a time. A record is delimited by its
// POSITION (FASTQ: four lines; FASTA: up to the next '>'), never by a leading byte --
// a quality line may start with '@'. Lines of any length.
typedef struct {
  gzFile fp;
  char *line; size_t lcap;
  char *seq;  size_t scap, slen;
  char name[256];
  int have;                      // `line` holds the next FASTA header, not yet consumed
} rdr_t;

static int rdr_open(rdr_t *r, const char *path) {
  memset(r, 0, sizeof *r);
  r->fp = gzopen(path, "r");
  if (!r->fp) { fprintf(stderr, "mosaic: cannot open %s\n", path); return -1; }
  gzbuffer(r->fp, 1 << 20);
  return 0;
}
static void rdr_close(rdr_t *r) { if (r->fp) gzclose(r->fp); free(r->line); free(r->seq); }

// One line into r->line, EOL stripped; its length, or -1 at the end of the file.
static long rdr_line(rdr_t *r) {
  size_t len = 0;
  for (;;) {
    if (len + 2 > r->lcap) { r->lcap = r->lcap ? 2 * r->lcap : (1u << 16); r->line = realloc(r->line, r->lcap); }
    if (!gzgets(r->fp, r->line + len, (int)(r->lcap - len))) { if (!len) return -1; break; }
    len += strlen(r->line + len);
    if (len && r->line[len - 1] == '\n') break;
  }
  while (len && (r->line[len - 1] == '\n' || r->line[len - 1] == '\r')) r->line[--len] = 0;
  return (long)len;
}

static void rdr_cat(rdr_t *r, const char *s, size_t L) {
  if (r->slen + L + 1 > r->scap) { r->scap = (r->slen + L + 1) * 2; r->seq = realloc(r->seq, r->scap); }
  memcpy(r->seq + r->slen, s, L); r->slen += L; r->seq[r->slen] = 0;
}

// 1 = a record is in r->name / r->seq, 0 = end of file, -1 = malformed.
static int rdr_next(rdr_t *r) {
  long L;
  if (r->have) L = (long)strlen(r->line);
  else do { L = rdr_line(r); } while (L == 0);          // blank lines between records
  if (L < 0) return 0;
  r->have = 0;
  char h = r->line[0];
  if (h != '>' && h != '@') {
    fprintf(stderr, "mosaic: not a FASTA/FASTQ header: %.40s\n", r->line); return -1;
  }
  size_t j = 0;
  while (j + 1 < sizeof r->name && r->line[j + 1] && r->line[j + 1] != ' ' && r->line[j + 1] != '\t')
    { r->name[j] = r->line[j + 1]; j++; }
  r->name[j] = 0;
  r->slen = 0; if (r->seq) r->seq[0] = 0;
  if (h == '@') {
    L = rdr_line(r);
    if (L < 0) { fprintf(stderr, "mosaic: truncated FASTQ\n"); return -1; }
    rdr_cat(r, r->line, (size_t)L);
    if (rdr_line(r) < 0 || rdr_line(r) < 0) {            // the '+' line and the quality
      fprintf(stderr, "mosaic: truncated FASTQ\n"); return -1;
    }
  } else {
    for (;;) {
      L = rdr_line(r);
      if (L < 0) break;
      if (L && r->line[0] == '>') { r->have = 1; break; }
      rdr_cat(r, r->line, (size_t)L);
    }
  }
  return 1;
}

// The two files of a pair name the same fragment, up to a trailing /1 or /2.
static int same_pair(const char *a, const char *b) {
  size_t la = strlen(a), lb = strlen(b);
  if (la > 2 && a[la - 2] == '/' && (a[la - 1] == '1' || a[la - 1] == '2')) la -= 2;
  if (lb > 2 && b[lb - 2] == '/' && (b[lb - 1] == '1' || b[lb - 1] == '2')) lb -= 2;
  return la == lb && !memcmp(a, b, la);
}

static int cmp_u16(const void *x, const void *y) {
  uint16_t a = *(const uint16_t *)x, b = *(const uint16_t *)y;
  return (a > b) - (a < b);
}

// The evidence of a read set, one row per INFORMATIVE observation.
typedef struct {
  uint32_t nc, nB;
  uint32_t *nobs_blk, *ninf_blk;       // per block: placed observations, informative ones
  size_t    n, cap;                    // informative observations kept
  uint16_t *blk;                       // [n]    the observation's block
  uint16_t *m;                         // [n*nc] hits the colour does not have, capped at 65535
  uint64_t  n_obs, n_nohit, n_unplaced, n_placed, n_inf, bases;
  int       want_m;                    // with hit counts too, also keep the per-observation vectors
} mobs_t;

static void mobs_free(mobs_t *O) {
  free(O->nobs_blk); free(O->ninf_blk); free(O->blk); free(O->m);
}

// Scan the reads exactly as run_classify does -- same rolling k-mer, same syncmer
// predicate, same lookup -- but keep, per observation, the hit count, the votes of
// every colour, and the blocks of the hits.
//
// With `hitcnt` (the cover mode) the per-observation vectors are not kept at all: what
// is kept is, per index record, how many times the reads hit it. A hit more than
// `stray` blocks from its read's median block is dropped there (a repeat copy the index
// calls single-locus can put a hit far away).
static int mosaic_collect(const pidx_t *ix, const char *reads, const char *mate, mobs_t *O,
                          uint16_t *hitcnt, int stray) {
  uint32_t nc = ix->nc, nw = ix->nw;
  kroll_t q; kr_init(&q, ix->k, ix->sk);
  rdr_t r1, r2; int paired = mate != NULL;
  if (rdr_open(&r1, reads)) return -1;
  if (paired && rdr_open(&r2, mate)) { rdr_close(&r1); return -1; }
  uint32_t *votes = calloc(nc ? nc : 1, sizeof(uint32_t));
  size_t bcap = 1u << 12, nb = 0;
  uint16_t *bl = malloc(bcap * sizeof(uint16_t)), *bs = malloc(bcap * sizeof(uint16_t));
  uint32_t *hs = hitcnt ? malloc(bcap * sizeof(uint32_t)) : NULL;   // the slot of each placed hit
  int rc = 0;
  for (;;) {
    int a = rdr_next(&r1);
    if (a < 0) { rc = -1; break; }
    if (a == 0) {
      if (paired && rdr_next(&r2) == 1) {
        fprintf(stderr, "mosaic: --mate has more records than the reads\n"); rc = -1;
      }
      break;
    }
    if (paired) {
      int b = rdr_next(&r2);
      if (b != 1) { fprintf(stderr, "mosaic: --mate ended early or is malformed\n"); rc = -1; break; }
      if (!same_pair(r1.name, r2.name)) {
        fprintf(stderr, "mosaic: record %llu: the pair is %s / %s, not one fragment\n",
                (unsigned long long)O->n_obs + 1, r1.name, r2.name);
        rc = -1; break;
      }
    }
    O->n_obs++;
    uint32_t hits = 0; nb = 0;
    memset(votes, 0, (size_t)nc * sizeof(uint32_t));
    for (int mt = 0; mt < (paired ? 2 : 1); mt++) {
      const rdr_t *r = mt ? &r2 : &r1;
      kr_reset(&q);
      O->bases += r->slen;
      for (size_t j = 0; j < r->slen; j++) {
        if (kr_push(&q, r->seq[j]) != 1) continue;
        size_t sl = kh_slot(ix->km, kr_key(&q), 0);
        if (sl == (size_t)-1 || ix->km->st[sl] != KH_IDX) continue;
        hits++;
        if (!hitcnt || O->want_m) {
          const uint64_t *mk = ix->cmask + sl * nw;
          for (uint32_t z = 0; z < nw; z++)
            for (uint64_t b = mk[z]; b; b &= b - 1) votes[(z << 6) + (uint32_t)__builtin_ctzll(b)]++;
        }
        if (ix->blk[sl] != BLK_NONE) {
          if (nb == bcap) {
            bcap *= 2; bl = realloc(bl, bcap * sizeof(uint16_t)); bs = realloc(bs, bcap * sizeof(uint16_t));
            if (hs) hs = realloc(hs, bcap * sizeof(uint32_t));
          }
          if (hs) hs[nb] = (uint32_t)sl;
          bl[nb++] = ix->blk[sl];
        }
      }
    }
    if (!hits) { O->n_nohit++; continue; }
    if (!nb)   { O->n_unplaced++; continue; }
    memcpy(bs, bl, nb * sizeof(uint16_t));
    qsort(bs, nb, sizeof(uint16_t), cmp_u16);
    uint32_t b = bs[(nb - 1) / 2];                       // the lower median
    if (b >= O->nB) b = O->nB - 1;                       // a coordinate at the very end
    O->n_placed++; O->nobs_blk[b]++;
    if (hitcnt) {
      for (size_t i = 0; i < nb; i++) {
        int d = (int)bl[i] - (int)b;
        if (d < -stray || d > stray) continue;
        if (hitcnt[hs[i]] < 65535) hitcnt[hs[i]]++;
        O->ninf_blk[bl[i] < O->nB ? bl[i] : O->nB - 1]++;   // here: the hits each block received
      }
      if (!O->want_m) continue;
    }
    uint32_t mn = UINT32_MAX, mx = 0;
    for (uint32_t c = 0; c < nc; c++) {
      uint32_t m = hits - votes[c];
      if (m < mn) mn = m;
      if (m > mx) mx = m;
    }
    if (mn == mx) continue;                              // nothing here tells a colour from another
    if (O->n == O->cap) {
      O->cap = O->cap ? 2 * O->cap : (1u << 14);
      O->blk = realloc(O->blk, O->cap * sizeof(uint16_t));
      O->m   = realloc(O->m, O->cap * (size_t)nc * sizeof(uint16_t));
    }
    O->blk[O->n] = (uint16_t)b;
    uint16_t *row = O->m + O->n * (size_t)nc;
    for (uint32_t c = 0; c < nc; c++) { uint32_t m = hits - votes[c]; row[c] = (uint16_t)(m > 65535 ? 65535 : m); }
    O->n++; O->n_inf++; if (!hitcnt) O->ninf_blk[b]++;
  }
  free(votes); free(bl); free(bs); free(hs);
  rdr_close(&r1); if (paired) rdr_close(&r2);
  return rc;
}

// The cover's scan of the reads on several threads. In cover mode the scan only COUNTS
// (hits per index record, placed reads, hits per block), so it commutes: reads are batched
// and cut into chunks for par_for; each chunk keeps its own block counters (merged after
// the batch) and per-record counts use a saturating atomic increment, giving
// min(total, 65535) in any order. Same result as mosaic_collect; not used by the mosaic mode.
typedef struct { uint64_t n_nohit, n_unplaced, n_placed, bases; uint32_t *nobs_blk, *ninf_blk; } cres_t;
typedef struct {
  const pidx_t *ix; uint16_t *hitcnt; int stray, per; uint32_t nB;     // per: reads per fragment (1 or 2)
  const char *sq; const size_t *off; const uint32_t *len;               // the batch's reads, one after another
  const uint32_t *chunk0; cres_t *res;                                  // chunk c = fragments [chunk0[c], chunk0[c+1])
} ccx_t;

static inline void sat_inc16(uint16_t *p) {
  uint16_t v = __atomic_load_n(p, __ATOMIC_RELAXED);
  while (v < 65535 && !__atomic_compare_exchange_n(p, &v, (uint16_t)(v + 1), 1, __ATOMIC_RELAXED, __ATOMIC_RELAXED)) ;
}

static void ccollect_chunk(const void *vx, uint32_t c) {
  const ccx_t *X = vx; const pidx_t *ix = X->ix; cres_t *R = &X->res[c];
  kroll_t q; kr_init(&q, ix->k, ix->sk);
  size_t bcap = 1u << 12, nb;
  uint16_t *bl = malloc(bcap * sizeof(uint16_t)), *bs = malloc(bcap * sizeof(uint16_t));
  uint32_t *hs = malloc(bcap * sizeof(uint32_t));
  for (uint32_t f = X->chunk0[c]; f < X->chunk0[c + 1]; f++) {
    uint32_t hits = 0; nb = 0;
    for (int mt = 0; mt < X->per; mt++) {
      size_t ri = (size_t)f * (size_t)X->per + (size_t)mt;
      const char *seq = X->sq + X->off[ri]; uint32_t L = X->len[ri];
      kr_reset(&q); R->bases += L;
      for (size_t j = 0; j < L; j++) {
        if (kr_push(&q, seq[j]) != 1) continue;
        size_t sl = kh_slot(ix->km, kr_key(&q), 0);
        if (sl == (size_t)-1 || ix->km->st[sl] != KH_IDX) continue;
        hits++;
        if (ix->blk[sl] != BLK_NONE) {
          if (nb == bcap) {
            bcap *= 2; bl = realloc(bl, bcap * sizeof(uint16_t)); bs = realloc(bs, bcap * sizeof(uint16_t));
            hs = realloc(hs, bcap * sizeof(uint32_t));
          }
          hs[nb] = (uint32_t)sl; bl[nb++] = ix->blk[sl];
        }
      }
    }
    if (!hits) { R->n_nohit++; continue; }
    if (!nb)   { R->n_unplaced++; continue; }
    memcpy(bs, bl, nb * sizeof(uint16_t));
    qsort(bs, nb, sizeof(uint16_t), cmp_u16);
    uint32_t b = bs[(nb - 1) / 2];                       // the lower median
    if (b >= X->nB) b = X->nB - 1;
    R->n_placed++; R->nobs_blk[b]++;
    for (size_t i = 0; i < nb; i++) {
      int d = (int)bl[i] - (int)b;
      if (d < -X->stray || d > X->stray) continue;
      sat_inc16(&X->hitcnt[hs[i]]);
      R->ninf_blk[bl[i] < X->nB ? bl[i] : X->nB - 1]++;
    }
  }
  free(bl); free(bs); free(hs);
}

static int mosaic_collect_par(const pidx_t *ix, const char *reads, const char *mate, mobs_t *O,
                              uint16_t *hitcnt, int stray, int nthreads) {
  rdr_t r1, r2; int paired = mate != NULL, per = paired ? 2 : 1;
  if (rdr_open(&r1, reads)) return -1;
  if (paired && rdr_open(&r2, mate)) { rdr_close(&r1); return -1; }
  const size_t BATCH = (size_t)64 << 20, CHUNK = (size_t)1 << 20;        // bases per batch, per chunk
  size_t sqcap = BATCH + ((size_t)1 << 20), sqn = 0, ncap = 1u << 16, nr = 0;
  char *sq = malloc(sqcap); size_t *off = malloc(ncap * sizeof(size_t)); uint32_t *len = malloc(ncap * sizeof(uint32_t));
  uint32_t nB = O->nB, ncc = 0;
  uint32_t *chunk0 = NULL; cres_t *res = NULL; uint32_t rescap = 0;
  int rc = 0, done = 0;
  while (!done && !rc) {
    sqn = 0; nr = 0;
    while (sqn < BATCH) {
      int a = rdr_next(&r1);
      if (a < 0) { rc = -1; break; }
      if (a == 0) {
        if (paired && rdr_next(&r2) == 1) { fprintf(stderr, "mosaic: --mate has more records than the reads\n"); rc = -1; }
        done = 1; break;
      }
      if (paired) {
        int b = rdr_next(&r2);
        if (b != 1) { fprintf(stderr, "mosaic: --mate ended early or is malformed\n"); rc = -1; break; }
        if (!same_pair(r1.name, r2.name)) {
          fprintf(stderr, "mosaic: record %llu: the pair is %s / %s, not one fragment\n",
                  (unsigned long long)O->n_obs + nr / 2 + 1, r1.name, r2.name);
          rc = -1; break;
        }
      }
      for (int mt = 0; mt < per; mt++) {
        const rdr_t *r = mt ? &r2 : &r1;
        if (nr == ncap) { ncap *= 2; off = realloc(off, ncap * sizeof(size_t)); len = realloc(len, ncap * sizeof(uint32_t)); }
        if (sqn + r->slen > sqcap) { sqcap = (sqn + r->slen) * 2; sq = realloc(sq, sqcap); }
        memcpy(sq + sqn, r->seq, r->slen); off[nr] = sqn; len[nr] = (uint32_t)r->slen; nr++; sqn += r->slen;
      }
    }
    if (rc || !nr) break;
    uint32_t nfrag = (uint32_t)(nr / (size_t)per);
    // chunk boundaries: about CHUNK bases each, whole fragments
    ncc = 0; size_t acc = 0;
    uint32_t cc_cap = nfrag + 2;
    chunk0 = realloc(chunk0, (cc_cap + 1) * sizeof(uint32_t));
    chunk0[0] = 0;
    for (uint32_t f = 0; f < nfrag; f++) {
      for (int mt = 0; mt < per; mt++) acc += len[(size_t)f * per + mt];
      if (acc >= CHUNK || f + 1 == nfrag) { chunk0[++ncc] = f + 1; acc = 0; }
    }
    if (ncc > rescap) {
      res = realloc(res, ncc * sizeof(cres_t));
      for (uint32_t c = rescap; c < ncc; c++) {
        res[c].nobs_blk = calloc(nB, sizeof(uint32_t)); res[c].ninf_blk = calloc(nB, sizeof(uint32_t));
      }
      rescap = ncc;
    }
    for (uint32_t c = 0; c < ncc; c++) {
      uint32_t *a = res[c].nobs_blk, *b = res[c].ninf_blk;
      memset(&res[c], 0, sizeof(cres_t)); res[c].nobs_blk = a; res[c].ninf_blk = b;
      memset(a, 0, nB * sizeof(uint32_t)); memset(b, 0, nB * sizeof(uint32_t));
    }
    ccx_t X = { ix, hitcnt, stray, per, nB, sq, off, len, chunk0, res };
    par_for(ncc, nthreads, ccollect_chunk, &X);
    O->n_obs += nfrag;
    for (uint32_t c = 0; c < ncc; c++) {
      O->n_nohit += res[c].n_nohit; O->n_unplaced += res[c].n_unplaced; O->n_placed += res[c].n_placed;
      O->bases += res[c].bases;
      for (uint32_t b = 0; b < nB; b++) { O->nobs_blk[b] += res[c].nobs_blk[b]; O->ninf_blk[b] += res[c].ninf_blk[b]; }
    }
  }
  for (uint32_t c = 0; c < rescap; c++) { free(res[c].nobs_blk); free(res[c].ninf_blk); }
  free(res); free(chunk0); free(sq); free(off); free(len);
  rdr_close(&r1); if (paired) rdr_close(&r2);
  return rc;
}

// The Viterbi painting. `cs` lists the Cs colours the two tracks may use (the others
// are excluded); returns path[nB][2] as indices into the FULL colour list.
// temper / temper_n / temper_w: TEMPER_OFF, GLOBAL (every block's E times temper_w) or
// BLOCK (E times min(1, temper_n / n_obs of the block)), resolved by run_mosaic.
static uint16_t *mosaic_paint(const mobs_t *O, const uint32_t *cs, uint32_t Cs,
                              double lam, double rho, double eps,
                              int temper, double temper_n, double temper_w, double *score_out) {
  uint32_t nB = O->nB, nc = O->nc;
  // p(m) is a table, and so is log((p(m1) + p(m2)) / 2): an observation costs one
  // lookup per ordered pair. The table stops where p has reached eps to 1e-12.
  uint32_t cap = 8;
  while (cap < 4096 && (1.0 - eps) * exp(-lam * (double)cap) > 1e-12 * eps) cap++;
  uint32_t T = cap + 1;
  double *tab = malloc(T * sizeof(double)), *lt = malloc((size_t)T * T * sizeof(double));
  for (uint32_t i = 0; i < T; i++) tab[i] = eps + (1.0 - eps) * exp(-lam * (double)i);
  for (uint32_t i = 0; i < T; i++)
    for (uint32_t j = 0; j < T; j++) lt[(size_t)i * T + j] = log(0.5 * (tab[i] + tab[j]));

  // observations by block (counting sort)
  size_t *bs = calloc((size_t)nB + 1, sizeof(size_t));
  for (size_t r = 0; r < O->n; r++) bs[O->blk[r] + 1]++;
  for (uint32_t b = 0; b < nB; b++) bs[b + 1] += bs[b];
  uint32_t *ord = malloc((O->n ? O->n : 1) * sizeof(uint32_t));
  { size_t *fill = malloc(((size_t)nB + 1) * sizeof(size_t));
    memcpy(fill, bs, ((size_t)nB + 1) * sizeof(size_t));
    for (size_t r = 0; r < O->n; r++) ord[fill[O->blk[r]]++] = (uint32_t)r;
    free(fill); }

  size_t cells = (size_t)Cs * Cs;
  uint8_t  *code = malloc((size_t)nB * cells);          // which transition won, per state
  uint16_t *argrow = malloc((size_t)nB * Cs * sizeof(uint16_t));
  uint16_t *argcol = malloc((size_t)nB * Cs * sizeof(uint16_t));
  uint16_t *glob = malloc((size_t)nB * 2 * sizeof(uint16_t));
  double *V = malloc(cells * sizeof(double)), *E = malloc(cells * sizeof(double));
  double *rm = malloc(Cs * sizeof(double)), *cm = malloc(Cs * sizeof(double));
  uint16_t *ar = malloc(Cs * sizeof(uint16_t)), *ac = malloc(Cs * sizeof(uint16_t));
  uint32_t *mm = malloc(Cs * sizeof(uint32_t));
  double ls = log1p(-rho), lj = log(rho / (double)Cs);
  double s2 = 2.0 * ls, sj = ls + lj, g2 = 2.0 * lj;

  for (uint32_t t = 0; t < nB; t++) {
    int has = bs[t + 1] > bs[t];
    if (has) {                                          // E(a,b) = sum over the block's observations
      memset(E, 0, cells * sizeof(double));
      for (size_t i = bs[t]; i < bs[t + 1]; i++) {
        const uint16_t *mr = O->m + (size_t)ord[i] * nc;
        for (uint32_t a = 0; a < Cs; a++) { uint32_t v = mr[cs[a]]; mm[a] = v < cap ? v : cap; }
        for (uint32_t a = 0; a < Cs; a++) {
          const double *lrow = lt + (size_t)mm[a] * T;
          double *Ea = E + (size_t)a * Cs;
          for (uint32_t b = a; b < Cs; b++) Ea[b] += lrow[mm[b]];
        }
      }
      for (uint32_t a = 1; a < Cs; a++)                 // mirror the upper triangle
        for (uint32_t b = 0; b < a; b++) E[(size_t)a * Cs + b] = E[(size_t)b * Cs + a];
      // Tempering scales E, not the transitions: scaling log(1-rho) and log(rho/Cs) by
      // 1/w has the same argmax in exact arithmetic, but in floating point it picks
      // another pair in 0.2-0.4% of the blocks of real paintings.
      double wt = 1.0;
      if (temper == TEMPER_GLOBAL) wt = temper_w;
      else if (temper == TEMPER_BLOCK && (double)O->nobs_blk[t] > temper_n)
        wt = temper_n / (double)O->nobs_blk[t];
      if (wt != 1.0) for (size_t z = 0; z < cells; z++) E[z] *= wt;
    }
    if (t == 0) {
      double v0 = -2.0 * log((double)Cs);
      for (size_t z = 0; z < cells; z++) V[z] = v0 + (has ? E[z] : 0.0);
      continue;
    }
    // row/column maxima of the previous V, first index on ties
    for (uint32_t b = 0; b < Cs; b++) { rm[b] = V[b]; ar[b] = 0; }
    for (uint32_t a = 1; a < Cs; a++)
      for (uint32_t b = 0; b < Cs; b++) { double x = V[(size_t)a * Cs + b]; if (x > rm[b]) { rm[b] = x; ar[b] = (uint16_t)a; } }
    for (uint32_t a = 0; a < Cs; a++) {
      const double *Va = V + (size_t)a * Cs;
      double mxv = Va[0]; uint16_t am = 0;
      for (uint32_t b = 1; b < Cs; b++) if (Va[b] > mxv) { mxv = Va[b]; am = (uint16_t)b; }
      cm[a] = mxv; ac[a] = am;
    }
    uint32_t gb = 0;
    for (uint32_t b = 1; b < Cs; b++) if (rm[b] > rm[gb]) gb = b;
    double g = rm[gb]; uint32_t ga = ar[gb];
    double c3 = g + g2;
    uint8_t *cd = code + (size_t)t * cells;
    for (uint32_t a = 0; a < Cs; a++) {
      double c2 = cm[a] + sj;
      for (uint32_t b = 0; b < Cs; b++) {
        size_t z = (size_t)a * Cs + b;
        double best = V[z] + s2; uint8_t k = 0, c;
        double cand = rm[b] + sj;
        if (cand > best) { best = cand; k = 1; }
        if (c2 > best) { best = c2; k = 2; }
        if (c3 > best) { best = c3; k = 3; }
        c = k; cd[z] = c;
        V[z] = has ? best + E[z] : best;
      }
    }
    memcpy(argrow + (size_t)t * Cs, ar, Cs * sizeof(uint16_t));
    memcpy(argcol + (size_t)t * Cs, ac, Cs * sizeof(uint16_t));
    glob[2 * (size_t)t] = (uint16_t)ga; glob[2 * (size_t)t + 1] = (uint16_t)gb;
  }
  size_t fl = 0;
  for (size_t z = 1; z < cells; z++) if (V[z] > V[fl]) fl = z;
  if (score_out) *score_out = V[fl];
  uint32_t a = (uint32_t)(fl / Cs), b = (uint32_t)(fl % Cs);
  uint16_t *path = malloc((size_t)nB * 2 * sizeof(uint16_t));
  path[2 * (size_t)(nB - 1)] = (uint16_t)cs[a]; path[2 * (size_t)(nB - 1) + 1] = (uint16_t)cs[b];
  for (uint32_t t = nB - 1; t > 0; t--) {
    uint8_t k = code[(size_t)t * cells + (size_t)a * Cs + b];
    if (k == 1)      a = argrow[(size_t)t * Cs + b];
    else if (k == 2) b = argcol[(size_t)t * Cs + a];
    else if (k == 3) { a = glob[2 * (size_t)t]; b = glob[2 * (size_t)t + 1]; }
    path[2 * (size_t)(t - 1)] = (uint16_t)cs[a]; path[2 * (size_t)(t - 1) + 1] = (uint16_t)cs[b];
  }
  free(tab); free(lt); free(bs); free(ord); free(code); free(argrow); free(argcol); free(glob);
  free(V); free(E); free(rm); free(cm); free(ar); free(ac); free(mm);
  return path;
}

// ---- mosaic: splice the painting into a graph (--mosaic-gfa) -------------------
//
// A painting says "from block 120 to 140, copy HG00733#1, then HG01106#2". A graph needs
// that as a WALK, so the tracks are cut and joined where two haplotypes are guaranteed
// to meet: at a node of the reference walk that both of them visit, in the same
// orientation. Joined there, the walk is an edge-valid path -- x's steps up to the node,
// then y's steps after it -- and nothing is invented.
//
// Every W-line of a colour is first ORIENTED along the reference (a W-line that runs
// mostly against it is reversed) and each of its steps
// gets a key
//     2 * (reference position of the last anchor at or before it) + (1 if not itself an
//     anchor)
// (2*position-1 for steps before the first anchor), so a bubble stays with the anchor in
// front of it. A switch x -> y at coordinate X takes the first shared anchor at or after
// X, looking up to 200 kb ahead and then back to the previous cut; a "strong" node, one
// where x's walk keeps going in the W-line it entered and y's keeps going after it,
// beats a "weak" one. A tract of x is every step of x whose key lies between its two
// cuts. If the two share no anchor the cut is a plain coordinate cut and the fragments
// are NOT joined. A W-line with no anchor at all has no coordinate: it is kept only when
// the track is one colour, and dropped (and counted) otherwise.
//
// Output: the reference colours named by --mosaic-keep (default: the coordinate
// reference) and the two tracks (sample `mosaic`, haplotype 1 and 2, one W-line `frag<k>`
// per maximal joined run), over the subgraph those walks induce -- the segments they
// visit and every edge between two of them -- which is what a graph aligner is given.
static int64_t splice_fwd = 200000;    // --mosaic-fwd: how far ahead of a switch to look for a shared node

typedef struct {                       // the coordinate reference
  uint32_t  nref;                      // distinct reference segments (first visit)
  int32_t  *kof;                       // segment -> k, or -1
  int64_t  *rp;                        // k -> coordinate of the first visit
  uint8_t  *rs;                        // k -> strand of the first visit (1 = '-')
  uint32_t *pos_to_k;                  // k in coordinate order
  int64_t  *spos;                      // rp[pos_to_k[i]], the sorted coordinates
  int64_t   length;                    // end coordinate
  uint64_t  repeats;                   // repeat visits dropped
} sref_t;

typedef struct { int64_t p; uint32_t k; } pk_t;
static int pk_cmp(const void *x, const void *y) {
  const pk_t *a = x, *b = y;
  if (a->p != b->p) return a->p < b->p ? -1 : 1;
  return (a->k > b->k) - (a->k < b->k);
}

static void sref_free(sref_t *R) {
  if (!R) return;
  free(R->kof); free(R->rp); free(R->rs); free(R->pos_to_k); free(R->spos); free(R);
}

static sref_t *sref_build(const gfa_t *gfa, const bg_t *g, const char *refname) {
  uint64_t nsteps = 0;
  for (uint32_t w = 0; w < g->np; w++)
    if (g->pname[w] && !strcmp(g->pname[w], refname)) nsteps += g->pn[w];
  if (!nsteps) { fprintf(stderr, "mosaic: no walk of colour \"%s\"\n", refname); return NULL; }
  sref_t *R = calloc(1, sizeof *R);
  R->kof = malloc((size_t)gfa->n_seg * sizeof(int32_t));
  for (uint32_t z = 0; z < gfa->n_seg; z++) R->kof[z] = -1;
  R->rp = malloc(nsteps * sizeof(int64_t)); R->rs = malloc(nsteps);
  for (uint32_t w = 0; w < g->np; w++) {
    if (!g->pname[w] || strcmp(g->pname[w], refname)) continue;
    int64_t st = gfa->walk[g->psrc[w]].st, pos = st > 0 ? st : 0;
    for (uint32_t s = 0; s < g->pn[w]; s++) {
      uint32_t v = g->pv[w][s], seg = v >> 1;
      if (R->kof[seg] < 0) { R->kof[seg] = (int32_t)R->nref; R->rp[R->nref] = pos; R->rs[R->nref] = v & 1; R->nref++; }
      else R->repeats++;
      pos += g->slen[seg];
    }
    if (pos > R->length) R->length = pos;
  }
  pk_t *a = malloc((size_t)R->nref * sizeof(pk_t));
  for (uint32_t k = 0; k < R->nref; k++) { a[k].p = R->rp[k]; a[k].k = k; }
  qsort(a, R->nref, sizeof(pk_t), pk_cmp);
  R->pos_to_k = malloc((size_t)R->nref * sizeof(uint32_t));
  R->spos = malloc((size_t)R->nref * sizeof(int64_t));
  for (uint32_t i = 0; i < R->nref; i++) { R->pos_to_k[i] = a[i].k; R->spos[i] = a[i].p; }
  free(a);
  return R;
}

// ---- the .gfai: a GFA's W-lines in chunks, so a splice reads only what it uses ----------
//
// A pangenome GFA is almost all W-lines, and a painting uses most colours but only a few
// haplotypes' worth of steps, so parsing everything made the splice the memory peak.
// --index --mosaic-ref writes <gfa>.gfai next to the GFA: every W-line in chunks of `ch`
// steps taken ALONG THE REFERENCE, with what swl_build computes for it -- whether it runs
// against the reference, its first anchor, the key (pf2) each chunk starts from and the
// range of keys in it. The splice reads the S- and L-lines and the colours it keeps
// whole, and every other W-line chunk by chunk, only where a cut or a tract asks for its
// keys; the output graph is identical. A .gfa.gz, a GFA without a .gfai, or one whose size
// or mtime no longer match is read whole.
#define GFAI_MAGIC   "GFAI"
#define GFAI_VERSION 1u
static uint32_t gfai_chunk = 1024;     // --gfai-chunk: steps per chunk when a .gfai is written
static int gfai_off = 0;               // --no-gfai: read the whole GFA even next to a .gfai

typedef struct numtab_s numtab_t;      // the step parser's name -> id table (defined with the GFA reader)

typedef struct {                       // a chunk: steps [c*ch, c*ch + ch) of a W-line, along the reference
  uint64_t off;                        // its bytes in the file, in file order: [off, off + len)
  uint32_t len;
  int64_t  base;                       // swl_build's `base` where the chunk starts
  int64_t  kmin, kmax;                 // the smallest and largest key in it
} gfai_ck_t;

typedef struct {                       // a W-line
  uint64_t loff, llen;                 // the line: the offset of its 'W', its length without the newline
  uint32_t n, first;                   // steps; the first anchor along the reference
  uint8_t  flip, mono, has;            // as swl_build finds them: reversed, keys never decrease, has an anchor
  int64_t  rp_first;                   // the reference coordinate of the first anchor
  uint64_t bp;                         // its sequence
  uint64_t ck0; uint32_t nck;          // its chunks in the table
  uint32_t hlen; char *hdr;            // the line from 'W' to the tab before the steps, inclusive
} gfai_w_t;

typedef struct {
  uint64_t  size, hash;                // the GFA it describes: size, FNV-1a of its bytes (as a .pex records it)
  int64_t   mt_s, mt_ns;               // ...and its modification time
  char     *ref;                       // the coordinate reference the keys are of
  uint32_t  ch, n_seg, nref;
  uint32_t  nw; gfai_w_t *w;           // every W-line, in file order
  uint64_t  nck; gfai_ck_t *ck;
  // while splicing
  FILE     *f;                         // the GFA
  numtab_t *nt; const gfa_t *gfa; const sref_t *R;
  char     *buf; size_t bcap;
  uint64_t  nload, bload;              // chunks read, bytes read
  int       err;                       // a chunk did not parse as the .gfai said: the GFA changed
} gfai_t;

static int gfai_steps(gfai_t *X, const gfai_ck_t *K, uint32_t *out, uint32_t cnt);   // with the GFA reader

typedef struct scol_s scol_t;

typedef struct {                       // one W-line of a colour, oriented along the reference
  uint32_t  w;                         // its walk in bg_t
  uint32_t  n;                         // steps
  uint32_t *v;                         // oriented vertices (a copy, reversed if flipped); NULL when lazy
  int64_t  *pf2;                       // per-step key; NULL if the W-line has no anchor, or is lazy
  uint8_t   mono;                      // pf2 never decreases: the steps of a window are one interval, found by binary search
  uint8_t   has;                       // the W-line has an anchor (it has keys)
  uint8_t   lazy;                      // its steps are read from the GFA chunk by chunk (.gfai)
  uint32_t  nanc;
  uint32_t *anc_idx;                   // steps that are anchors...
  int32_t  *anc_k;                     // ...and their reference index
  uint64_t  bp;                        // its sequence
  // lazy: chunk c holds steps [c*ch, c*ch + ch) once read
  gfai_t   *X; const gfai_w_t *xw;
  uint32_t **lv; int64_t **lk;         // per chunk: vertices and keys, NULL until read
  scol_t   *C; uint32_t wi;            // the colour and the index of this W-line in it
} swl_t;

typedef struct { int32_t k, w; uint32_t i; int8_t s; } occ_e_t;

struct scol_s {                        // a colour the painting uses
  uint32_t c;                          // its colour id in the index
  uint32_t nwl; swl_t *wl;
  int32_t  *ow; uint32_t *oi; int8_t *os;   // per reference k: W-line, step, strand (lazy)
  // A colour with a lazy W-line keeps the same three in a hash of the reference segments its
  // chunks have shown so far: they are complete for every key range that was asked for.
  uint8_t  lazy;
  occ_e_t *oh; size_t ohm, ohn;
};

static void swl_build(swl_t *L, const bg_t *g, const sref_t *R, uint32_t w) {
  uint32_t n = g->pn[w];
  memset(L, 0, sizeof *L);
  L->w = w; L->n = n;
  L->v = malloc((n ? n : 1) * sizeof(uint32_t));
  if (n) memcpy(L->v, g->pv[w], n * sizeof(uint32_t));
  uint64_t nanc = 0, nsame = 0;
  for (uint32_t i = 0; i < n; i++) {
    int32_t k = R->kof[L->v[i] >> 1];
    if (k >= 0) { nanc++; if ((L->v[i] & 1) == R->rs[k]) nsame++; }
    L->bp += g->slen[L->v[i] >> 1];
  }
  if (nanc && 2 * nsame < nanc) {                  // mostly against the reference: reverse it
    for (uint32_t i = 0, j = n; i < j; ) { j--; uint32_t t = L->v[i]; L->v[i] = L->v[j] ^ 1u; L->v[j] = t ^ 1u; i++; }
  }
  L->nanc = 0;
  L->anc_idx = malloc((nanc ? nanc : 1) * sizeof(uint32_t));
  L->anc_k = malloc((nanc ? nanc : 1) * sizeof(int32_t));
  for (uint32_t i = 0; i < n; i++) {
    int32_t k = R->kof[L->v[i] >> 1];
    if (k >= 0) { L->anc_idx[L->nanc] = i; L->anc_k[L->nanc] = k; L->nanc++; }
  }
  if (!L->nanc) return;
  L->has = 1;
  L->pf2 = malloc((size_t)n * sizeof(int64_t));
  uint32_t first = L->anc_idx[0], ai = 0;
  int64_t base = 0;
  for (uint32_t i = 0; i < n; i++) {
    int anc = ai < L->nanc && L->anc_idx[ai] == i;
    if (anc) { base = R->rp[L->anc_k[ai]]; ai++; }
    L->pf2[i] = i < first ? 2 * R->rp[L->anc_k[0]] - 1 : 2 * base + (anc ? 0 : 1);
  }
  L->mono = 1;
  for (uint32_t i = 1; i < n; i++) if (L->pf2[i] < L->pf2[i - 1]) { L->mono = 0; break; }
}

// A W-line of the .gfai: nothing is read until a step of it is asked for.
static void swl_lazy(swl_t *L, gfai_t *X, const gfai_w_t *W, uint32_t w, scol_t *C, uint32_t wi) {
  memset(L, 0, sizeof *L);
  L->w = w; L->n = W->n; L->mono = W->mono; L->has = W->has; L->bp = W->bp; L->lazy = 1;
  L->X = X; L->xw = W; L->C = C; L->wi = wi;
  L->lv = calloc(W->nck ? W->nck : 1, sizeof(uint32_t *));
  L->lk = calloc(W->nck ? W->nck : 1, sizeof(int64_t *));
}

// The reference segment k is visited by W-line w of C at step i, strand s. What scol_occ keeps
// is the visit of the LARGEST W-line (ties: the later one) and its FIRST visit there, which
// does not depend on the order the visits arrive in.
static void occ_put(scol_t *C, int32_t k, uint32_t w, uint32_t i, int8_t s) {
  if (2 * (C->ohn + 1) > C->ohm) {
    size_t m = C->ohm ? 2 * C->ohm : 1024;
    occ_e_t *h = malloc(m * sizeof(occ_e_t));
    for (size_t z = 0; z < m; z++) h[z].k = -1;
    for (size_t z = 0; z < C->ohm; z++) if (C->oh[z].k >= 0) {
      size_t q = ((uint32_t)C->oh[z].k * 2654435761u) & (m - 1);
      while (h[q].k >= 0) q = (q + 1) & (m - 1);
      h[q] = C->oh[z];
    }
    free(C->oh); C->oh = h; C->ohm = m;
  }
  size_t q = ((uint32_t)k * 2654435761u) & (C->ohm - 1);
  while (C->oh[q].k >= 0 && C->oh[q].k != k) q = (q + 1) & (C->ohm - 1);
  occ_e_t *e = &C->oh[q];
  if (e->k < 0) { e->k = k; e->w = (int32_t)w; e->i = i; e->s = s; C->ohn++; return; }
  if ((uint32_t)e->w == w) { if (i < e->i) { e->i = i; e->s = s; } return; }
  uint32_t nn = C->wl[w].n, ne = C->wl[e->w].n;
  if (nn > ne || (nn == ne && w > (uint32_t)e->w)) { e->w = (int32_t)w; e->i = i; e->s = s; }
}

// The W-line of C holding the reference segment k, or -1; *i, *s: its step and strand there.
static int occ_get(const scol_t *C, uint32_t k, uint32_t *i, int8_t *s) {
  if (!C->lazy) {
    int w = C->ow[k];
    if (w >= 0) { *i = C->oi[k]; *s = C->os[k]; }
    return w;
  }
  if (!C->ohm) return -1;
  size_t q = (k * 2654435761u) & (C->ohm - 1);
  while (C->oh[q].k >= 0) {
    if ((uint32_t)C->oh[q].k == k) { *i = C->oh[q].i; *s = C->oh[q].s; return C->oh[q].w; }
    q = (q + 1) & (C->ohm - 1);
  }
  return -1;
}

// Read chunk c of a lazy W-line: its vertices along the reference and their keys, computed as
// swl_build does from the `base` the chunk starts with; its anchors go into the colour's hash.
static void swl_load(swl_t *L, uint32_t c) {
  gfai_t *X = L->X; const gfai_w_t *W = L->xw; const sref_t *R = X->R;
  uint32_t s0 = c * X->ch, cnt = L->n - s0 < X->ch ? L->n - s0 : X->ch;
  uint32_t *v = calloc(cnt, sizeof(uint32_t)); int64_t *p = calloc(cnt, sizeof(int64_t));
  L->lv[c] = v; L->lk[c] = p;
  const gfai_ck_t *K = &X->ck[W->ck0 + c];
  if (gfai_steps(X, K, v, cnt)) { X->err = 1; memset(v, 0, cnt * sizeof(uint32_t)); return; }
  if (W->flip)
    for (uint32_t i = 0, j = cnt; i < j; ) { j--; uint32_t t = v[i]; v[i] = v[j] ^ 1u; v[j] = t ^ 1u; i++; }
  int64_t base = K->base;
  for (uint32_t j = 0; j < cnt; j++) {
    uint32_t i = s0 + j;
    int32_t k = R->kof[v[j] >> 1];
    if (k >= 0) { base = R->rp[k]; occ_put(L->C, k, L->wi, i, (int8_t)(v[j] & 1)); }
    if (W->has) p[j] = i < W->first ? 2 * W->rp_first - 1 : 2 * base + (k >= 0 ? 0 : 1);
  }
}

static inline uint32_t swl_v(swl_t *L, uint32_t i) {
  if (!L->lazy) return L->v[i];
  uint32_t c = i / L->X->ch;
  if (!L->lv[c]) swl_load(L, c);
  return L->lv[c][i - c * L->X->ch];
}
static inline int64_t swl_k(swl_t *L, uint32_t i) {
  if (!L->lazy) return L->pf2[i];
  uint32_t c = i / L->X->ch;
  if (!L->lk[c]) swl_load(L, c);
  return L->lk[c][i - c * L->X->ch];
}

// Read every chunk of a lazy W-line that holds a key in [lo, hi].
static void swl_need(swl_t *L, int64_t lo, int64_t hi) {
  if (!L->lazy || !L->has) return;
  const gfai_ck_t *K = &L->X->ck[L->xw->ck0];
  uint32_t a = 0, nk = L->xw->nck;
  if (L->mono) {                                   // the chunks' key ranges follow each other: start at the first that reaches lo
    uint32_t hi2 = nk;
    while (a < hi2) { uint32_t m = a + (hi2 - a) / 2; if (K[m].kmax < lo) a = m + 1; else hi2 = m; }
  }
  for (uint32_t c = a; c < nk; c++) {
    if (L->mono && K[c].kmin > hi) break;
    if (K[c].kmax >= lo && K[c].kmin <= hi && !L->lv[c]) swl_load(L, c);
  }
}

// The first step of a monotone W-line whose key is > x (n if none).
static uint32_t swl_upper(swl_t *L, int64_t x) {
  if (!L->lazy) {
    uint32_t lo = 0, hi = L->n;
    while (lo < hi) { uint32_t mid = lo + (hi - lo) / 2; if (L->pf2[mid] > x) hi = mid; else lo = mid + 1; }
    return lo;
  }
  const gfai_ck_t *K = &L->X->ck[L->xw->ck0];
  uint32_t lo = 0, hi = L->xw->nck;
  while (lo < hi) { uint32_t mid = lo + (hi - lo) / 2; if (K[mid].kmax > x) hi = mid; else lo = mid + 1; }
  if (lo == L->xw->nck) return L->n;
  uint32_t c = lo, s0 = c * L->X->ch, cnt = L->n - s0 < L->X->ch ? L->n - s0 : L->X->ch;
  if (!L->lk[c]) swl_load(L, c);
  uint32_t a = 0, b = cnt;
  while (a < b) { uint32_t mid = a + (b - a) / 2; if (L->lk[c][mid] > x) b = mid; else a = mid + 1; }
  return s0 + a;
}

static void scol_need(scol_t *C, int64_t lo, int64_t hi) {
  if (C->lazy) for (uint32_t i = 0; i < C->nwl; i++) swl_need(&C->wl[i], lo, hi);
}

static void swl_free(swl_t *L) {
  free(L->v); free(L->pf2); free(L->anc_idx); free(L->anc_k);
  if (L->lazy) {
    for (uint32_t c = 0; c < L->xw->nck; c++) { free(L->lv[c]); free(L->lk[c]); }
    free(L->lv); free(L->lk);
  }
}

static void scol_free(scol_t *C) {
  for (uint32_t i = 0; i < C->nwl; i++) swl_free(&C->wl[i]);
  free(C->wl); free(C->ow); free(C->oi); free(C->os); free(C->oh);
}

// Where each reference segment sits in this colour: the LARGEST W-line that holds it
// (ties: the later one), and the FIRST step of that W-line that visits it.
static void scol_occ(const sref_t *R, scol_t *C) {
  if (C->ow) return;
  C->ow = malloc((size_t)R->nref * sizeof(int32_t));
  C->oi = malloc((size_t)R->nref * sizeof(uint32_t));
  C->os = malloc(R->nref ? R->nref : 1);
  for (uint32_t k = 0; k < R->nref; k++) C->ow[k] = -1;
  uint32_t *ord = malloc((C->nwl ? C->nwl : 1) * sizeof(uint32_t));
  for (uint32_t i = 0; i < C->nwl; i++) ord[i] = i;
  for (uint32_t a = 1; a < C->nwl; a++)             // insertion sort by length, stable
    for (uint32_t b = a; b > 0 && C->wl[ord[b - 1]].n > C->wl[ord[b]].n; b--)
      { uint32_t t = ord[b]; ord[b] = ord[b - 1]; ord[b - 1] = t; }
  for (uint32_t oi = 0; oi < C->nwl; oi++) {
    uint32_t wi = ord[oi]; const swl_t *L = &C->wl[wi];
    for (uint32_t j = L->nanc; j-- > 0; ) {         // last to first: the first visit wins
      int32_t k = L->anc_k[j];
      C->ow[k] = (int32_t)wi; C->oi[k] = L->anc_idx[j]; C->os[k] = (int8_t)(L->v[L->anc_idx[j]] & 1);
    }
  }
  free(ord);
}

typedef struct {                       // where a switch is made
  int      node;                       // 1: at a shared reference node; 0: a plain coordinate cut
  int64_t  b;                          // key bound: steps with key <= b belong to the earlier tract
  int64_t  X, pos;                     // the switch coordinate; the node's coordinate
  int      wlx, wly;                   // the W-line of x / y holding the node
  uint32_t ix, iy;                     // its step in each
  int      strong;
} cut_t;

static size_t lb64(const int64_t *a, size_t n, int64_t x) {      // first i with a[i] >= x
  size_t lo = 0, hi = n;
  while (lo < hi) { size_t m = (lo + hi) / 2; if (a[m] < x) lo = m + 1; else hi = m; }
  return lo;
}
static size_t ub64(const int64_t *a, size_t n, int64_t x) {      // first i with a[i] > x
  size_t lo = 0, hi = n;
  while (lo < hi) { size_t m = (lo + hi) / 2; if (a[m] <= x) lo = m + 1; else hi = m; }
  return lo;
}

// The cut node of the switch x -> y at X, or 0. Forward first (up to splice_fwd, and short
// of the next switch), then backward down to the previous cut; within a window a strong
// candidate anywhere beats a weak one.
static int find_cut(const sref_t *R, scol_t *cx, scol_t *cy, int64_t X, int64_t Xnext,
                    int64_t pos_prev, int main_wl_x, int64_t i_in, cut_t *out) {
  if (!cx->lazy) scol_occ(R, cx);
  if (!cy->lazy) scol_occ(R, cy);
  size_t n = R->nref;
  size_t lo_fwd = lb64(R->spos, n, X);
  int64_t hi_x = Xnext < 0 ? X + splice_fwd : (X + splice_fwd < Xnext - 1 ? X + splice_fwd : Xnext - 1);
  size_t hi_fwd = ub64(R->spos, n, hi_x);
  size_t lo_bwd = ub64(R->spos, n, pos_prev), hi_bwd = lo_fwd;
  for (int dir = 1; dir >= -1; dir -= 2) {
    size_t a = dir > 0 ? lo_fwd : lo_bwd, b = dir > 0 ? hi_fwd : hi_bwd;
    // Every visit of a reference segment carries the key 2 * its coordinate: reading the chunks
    // with keys in the window makes the occupancy of every segment of the window complete.
    if (b > a) { scol_need(cx, 2 * R->spos[a], 2 * R->spos[b - 1]); scol_need(cy, 2 * R->spos[a], 2 * R->spos[b - 1]); }
    long fs = -1, fb = -1, ls = -1, lbs = -1;
    for (size_t i = a; i < b; i++) {
      uint32_t k = R->pos_to_k[i], oix = 0, oiy = 0; int8_t osx = 0, osy = 0;
      int wx = occ_get(cx, k, &oix, &osx), wy = occ_get(cy, k, &oiy, &osy);
      if (wx < 0 || wy < 0 || osx != osy) continue;
      int strong = (main_wl_x < 0 || (wx == main_wl_x && (int64_t)oix > i_in)) &&
                   oiy + 1 < cy->wl[wy].n;
      if (fb < 0) fb = (long)i;
      lbs = (long)i;
      if (strong) { if (fs < 0) fs = (long)i; ls = (long)i; }
    }
    long pick = dir > 0 ? (fs >= 0 ? fs : fb) : (ls >= 0 ? ls : lbs);
    int strong = dir > 0 ? fs >= 0 : ls >= 0;
    if (pick < 0) continue;
    uint32_t k = R->pos_to_k[pick]; int8_t sx, sy;
    out->node = 1; out->pos = R->rp[k]; out->b = 2 * R->rp[k]; out->X = X; out->strong = strong;
    out->wlx = occ_get(cx, k, &out->ix, &sx); out->wly = occ_get(cy, k, &out->iy, &sy);
    return 1;
  }
  return 0;
}

typedef struct { int32_t slot; int64_t X0, X1; } tract_t;     // slot < 0: the track holds no colour here
typedef struct { uint32_t slot; int wl; uint32_t s, e; int sc, ec; int64_t key; int tract; } piece_t;
typedef struct { uint32_t *v; uint32_t n; uint64_t bp; } frag_t;

typedef struct { const piece_t *p; size_t i; } pord_t;
static int pord_cmp(const void *x, const void *y) {
  const pord_t *a = x, *b = y;
  if (a->p->key != b->p->key) return a->p->key < b->p->key ? -1 : 1;
  if (a->p->tract != b->p->tract) return a->p->tract < b->p->tract ? -1 : 1;
  return (a->i > b->i) - (a->i < b->i);
}

typedef struct {
  uint32_t tracts, ncut_node, ncut_coord, ngap, nstrong, nfrag;
  uint64_t bp, dropped_bp, dropped_wl;
  int64_t  far_cut;                    // the largest |cut - switch| distance
} strk_t;

static int splice_track(const sref_t *R, scol_t *cols, const tract_t *tr, uint32_t m,
                        const bg_t *g, frag_t **frags, strk_t *st) {
  // The steps of a lazy W-line are read through swl_v / swl_k, chunk by chunk, and only the
  // chunks whose keys can fall in a tract are asked for: what is never asked for is never read.
  memset(st, 0, sizeof *st); st->tracts = m;
  cut_t *cuts = calloc(m > 1 ? m - 1 : 1, sizeof(cut_t));
  int64_t pos_prev = -1, i_in = -1; int main_wl = -1;
  for (uint32_t j = 0; j + 1 < m; j++) {
    int64_t X = tr[j].X1, Xn = j + 2 < m ? tr[j + 1].X1 : -1;
    cut_t c; memset(&c, 0, sizeof c);
    if (tr[j].slot < 0 || tr[j + 1].slot < 0) {                // next to a stretch with no colour: no one to join
      cuts[j].node = 0; cuts[j].b = 2 * X - 1; cuts[j].X = X;
      pos_prev = X - 1; main_wl = -1; i_in = -1;
      st->ngap++;
    } else if (find_cut(R, &cols[tr[j].slot], &cols[tr[j + 1].slot], X, Xn, pos_prev, main_wl, i_in, &c)) {
      cuts[j] = c; pos_prev = c.pos; main_wl = c.wly; i_in = c.iy;
      st->ncut_node++; st->nstrong += c.strong;
      int64_t d = c.pos > X ? c.pos - X : X - c.pos;
      if (d > st->far_cut) st->far_cut = d;
    } else {
      cuts[j].node = 0; cuts[j].b = 2 * X - 1; cuts[j].X = X;
      pos_prev = X - 1; main_wl = -1; i_in = -1;
      st->ncut_coord++;
    }
  }
  size_t np = 0, pcap = 64; piece_t *pc = malloc(pcap * sizeof(piece_t));
  uint32_t maxn = 0;
  for (uint32_t j = 0; j < m; j++)
    for (uint32_t wi = 0; tr[j].slot >= 0 && wi < cols[tr[j].slot].nwl; wi++)
      if (cols[tr[j].slot].wl[wi].n > maxn) maxn = cols[tr[j].slot].wl[wi].n;
  uint8_t *inr = calloc(maxn ? maxn : 1, 1);               // all zero between uses: only [a, b) is ever written
  for (uint32_t j = 0; j < m; j++) {
    if (tr[j].slot < 0) continue;
    scol_t *cx = &cols[tr[j].slot];
    int64_t Lb = j == 0 ? INT64_MIN : cuts[j - 1].b, Hb = j == m - 1 ? INT64_MAX : cuts[j].b;
    const cut_t *ci = j > 0 ? &cuts[j - 1] : NULL, *co = j + 1 < m ? &cuts[j] : NULL;
    for (uint32_t wi = 0; wi < cx->nwl; wi++) {
      swl_t *L = &cx->wl[wi];
      if (!L->has) continue;
      // The steps in range are those with Lb < pf2 <= Hb. When pf2 never decreases they are one interval
      // [t0, t1), found by binary search, instead of a pass over the whole W-line for every tract.
      uint32_t t0 = 0, t1 = L->n;
      if (L->mono && !L->lazy) {
        uint32_t lo = 0, hi = L->n;                         // t0: first step with pf2 > Lb
        while (lo < hi) { uint32_t mid = lo + (hi - lo) / 2; if (L->pf2[mid] > Lb) hi = mid; else lo = mid + 1; }
        t0 = lo; lo = t0; hi = L->n;                        // t1: first step with pf2 > Hb
        while (lo < hi) { uint32_t mid = lo + (hi - lo) / 2; if (L->pf2[mid] > Hb) hi = mid; else lo = mid + 1; }
        t1 = lo;
      } else if (L->mono) {                                 // the same two searches, over the chunks first
        t0 = swl_upper(L, Lb); t1 = swl_upper(L, Hb);
        if (t1 < t0) t1 = t0;
      }
      uint32_t a = t0, b = t1;                              // the part of inr that can be nonzero
      if (ci && ci->node && ci->wly == (int)wi) { if (ci->iy < a) a = ci->iy; if (ci->iy + 2 > b) b = ci->iy + 2 < L->n ? ci->iy + 2 : L->n; }
      if (co && co->node && co->wlx == (int)wi) { if (co->ix < a) a = co->ix; if (co->ix + 1 > b) b = co->ix + 1; }
      if (b < a) b = a;
      if (!L->lazy) for (uint32_t i = t0; i < t1; i++) inr[i] = L->pf2[i] > Lb && L->pf2[i] <= Hb;
      else if (L->mono) for (uint32_t i = t0; i < t1; i++) { int64_t q = swl_k(L, i); inr[i] = q > Lb && q <= Hb; }
      else {                                                // keys go back and forth: only the chunks whose keys reach the tract
        const gfai_ck_t *K = &L->X->ck[L->xw->ck0];
        for (uint32_t c = 0; c < L->xw->nck; c++) {
          if (K[c].kmax <= Lb || K[c].kmin > Hb) continue;
          uint32_t s0 = c * L->X->ch, s1 = s0 + L->X->ch < L->n ? s0 + L->X->ch : L->n;
          for (uint32_t i = s0; i < s1; i++) { int64_t q = swl_k(L, i); inr[i] = q > Lb && q <= Hb; }
        }
      }
      if (ci && ci->node && ci->wly == (int)wi) { inr[ci->iy] = 0; if (ci->iy + 1 < L->n) inr[ci->iy + 1] = 1; }
      if (co && co->node && co->wlx == (int)wi) inr[co->ix] = 1;
      for (uint32_t i = a; i < b; ) {
        if (!inr[i]) { i++; continue; }
        uint32_t s = i, e = i;
        while (e + 1 < L->n && inr[e + 1]) e++;
        i = e + 1;
        int sc = -1, ec = -1;
        if (ci && ci->node && ci->wly == (int)wi && s <= ci->iy + 1 && ci->iy + 1 <= e) { s = ci->iy + 1; sc = (int)j - 1; }
        if (co && co->node && co->wlx == (int)wi && s <= co->ix && co->ix <= e) { e = co->ix; ec = (int)j; }
        if (e < s) continue;
        if (np == pcap) { pcap *= 2; pc = realloc(pc, pcap * sizeof(piece_t)); }
        pc[np++] = (piece_t){ tr[j].slot, (int)wi, s, e, sc, ec, swl_k(L, s), (int)j };
      }
      memset(inr + a, 0, b - a);
    }
  }
  for (uint32_t j = 0; j < m; j++) {
    if (tr[j].slot < 0) continue;
    int seen = 0;
    for (uint32_t z = 0; z < j; z++) if (tr[z].slot == tr[j].slot) seen = 1;
    if (seen) continue;
    const scol_t *cx = &cols[tr[j].slot];
    for (uint32_t wi = 0; wi < cx->nwl; wi++) {
      const swl_t *L = &cx->wl[wi];
      if (L->has) continue;
      if (m == 1) {
        if (np == pcap) { pcap *= 2; pc = realloc(pc, pcap * sizeof(piece_t)); }
        pc[np++] = (piece_t){ tr[j].slot, (int)wi, 0, L->n - 1, -1, -1, INT64_MIN / 2 + (int64_t)wi, (int)j };
      } else {
        st->dropped_wl++;
        st->dropped_bp += L->bp;
      }
    }
  }
  free(inr);
  int *by_start = malloc((m > 1 ? m - 1 : 1) * sizeof(int)), *by_end = malloc((m > 1 ? m - 1 : 1) * sizeof(int));
  for (uint32_t j = 0; j + 1 < m; j++) by_start[j] = by_end[j] = -1;
  for (size_t i = 0; i < np; i++) {
    if (pc[i].sc >= 0) by_start[pc[i].sc] = (int)i;
    if (pc[i].ec >= 0) by_end[pc[i].ec] = (int)i;
  }
  int *nxt = malloc((np ? np : 1) * sizeof(int)); uint8_t *prv = calloc(np ? np : 1, 1);
  for (size_t i = 0; i < np; i++) nxt[i] = -1;
  for (uint32_t j = 0; j + 1 < m; j++)
    if (by_end[j] >= 0 && by_start[j] >= 0) { nxt[by_end[j]] = by_start[j]; prv[by_start[j]] = 1; }
  pord_t *po = malloc((np ? np : 1) * sizeof(pord_t)); size_t nh = 0;
  for (size_t i = 0; i < np; i++) if (!prv[i]) { po[nh].p = &pc[i]; po[nh].i = i; nh++; }
  qsort(po, nh, sizeof(pord_t), pord_cmp);
  frag_t *fr = calloc(nh ? nh : 1, sizeof(frag_t));
  int bad = 0;
  for (size_t h = 0; h < nh; h++) {
    size_t tot = 0;
    for (int i = (int)po[h].i; i >= 0; i = nxt[i]) tot += pc[i].e - pc[i].s + 1;
    fr[h].v = malloc(tot * sizeof(uint32_t));
    int prev = -1;
    for (int i = (int)po[h].i; i >= 0; prev = i, i = nxt[i]) {
      swl_t *L = &cols[pc[i].slot].wl[pc[i].wl];
      if (prev >= 0) {                                       // the junction is one and the same node
        swl_t *A = &cols[pc[prev].slot].wl[pc[prev].wl];
        if (swl_v(A, pc[prev].e) != swl_v(L, pc[i].s - 1)) bad++;
      }
      for (uint32_t s = pc[i].s; s <= pc[i].e; s++) {
        uint32_t x = swl_v(L, s);
        fr[h].v[fr[h].n++] = x; fr[h].bp += g->slen[x >> 1];
      }
    }
    st->bp += fr[h].bp;
  }
  st->nfrag = (uint32_t)nh;
  *frags = fr;
  free(po); free(nxt); free(prv); free(by_start); free(by_end); free(pc); free(cuts);
  if (bad) { fprintf(stderr, "mosaic: splice: %d joined fragments do not meet at one node\n", bad); return -1; }
  return 0;
}

// Does the graph have the edge u -> v (oriented vertices)?
static int has_edge(const gfa_t *gfa, uint32_t u, uint32_t v) {
  const gfa_arc_t *a = gfa_arc_a(gfa, u);
  uint32_t n = gfa_arc_n(gfa, u);
  for (uint32_t i = 0; i < n; i++) if (a[i].w == v) return 1;
  return 0;
}

// Splice the painting and write the graph. `path` is [nB][2] colour ids of the index.
// X: the .gfai of the GFA, whose W-lines that `g` holds without steps are read from it; NULL = none.
static int mosaic_splice(const gfa_t *gfa, const bg_t *g, const pidx_t *ix, const uint16_t *path,
                         uint32_t nB, uint32_t T, const char *outfile, const char *keep_csv, gfai_t *X) {
  double t0 = now_s();
  sref_t *R = sref_build(gfa, g, ix->blk_ref);
  if (!R) return 1;
  fprintf(stderr, "splice: reference %s: %u segments, length %lld bp, %llu repeat visits dropped\n",
          ix->blk_ref, R->nref, (long long)R->length, (unsigned long long)R->repeats);
  if (X && (R->nref != X->nref || gfa->n_seg != X->n_seg)) {
    fprintf(stderr, "splice: the .gfai was written for %u segments and %u reference segments, this graph has %u and %u\n",
            X->n_seg, X->nref, gfa->n_seg, R->nref);
    sref_free(R); return 1;
  }
  if (X) X->R = R;
  uint32_t nc = ix->nc, nslot = 0;
  int32_t *slot_of = malloc((nc ? nc : 1) * sizeof(int32_t));
  for (uint32_t c = 0; c < nc; c++) slot_of[c] = -1;
  for (uint32_t t = 0; t < nB; t++) for (uint32_t k = 0; k < T; k++) {
    uint16_t c = path[T * (size_t)t + k];
    if (c != NOCOL && slot_of[c] < 0) slot_of[c] = (int32_t)nslot++;
  }
  scol_t *cols = calloc(nslot, sizeof(scol_t));
  for (uint32_t c = 0; c < nc; c++) if (slot_of[c] >= 0) {
    scol_t *C = &cols[slot_of[c]]; C->c = c;
    uint32_t nw = 0;
    for (uint32_t w = 0; w < g->np; w++) if (g->pname[w] && !strcmp(g->pname[w], ix->cname[c])) nw++;
    if (!nw) { fprintf(stderr, "splice: colour %s of the painting has no walk in this graph\n", ix->cname[c]); return 1; }
    C->wl = calloc(nw, sizeof(swl_t));
    for (uint32_t w = 0; w < g->np; w++) {
      if (!g->pname[w] || strcmp(g->pname[w], ix->cname[c])) continue;
      const gfai_w_t *xw = X && !g->pn[w] && g->psrc[w] < X->nw ? &X->w[g->psrc[w]] : NULL;
      if (xw && xw->n) { swl_lazy(&C->wl[C->nwl], X, xw, w, C, C->nwl); C->lazy = 1; }
      else swl_build(&C->wl[C->nwl], g, R, w);
      C->nwl++;
    }
    if (C->lazy)
      for (uint32_t wi = 0; wi < C->nwl; wi++) {
        const swl_t *L = &C->wl[wi];
        if (L->lazy) continue;
        for (uint32_t j = 0; j < L->nanc; j++)
          occ_put(C, L->anc_k[j], wi, L->anc_idx[j], (int8_t)(L->v[L->anc_idx[j]] & 1));
      }
  }
  fprintf(stderr, "splice: %u colours loaded (%.1fs)\n", nslot, now_s() - t0);

  frag_t **frags = calloc(T, sizeof(frag_t *)); strk_t *stat = calloc(T, sizeof(strk_t)); int rc = 0;
  uint64_t all_bp = 0;
  for (uint32_t k = 0; k < T && !rc; k++) {
    tract_t *tr = malloc((size_t)nB * sizeof(tract_t)); uint32_t m = 0;
    for (uint32_t t = 0; t < nB; t++) {
      uint16_t pc = path[T * (size_t)t + k];
      int32_t slot = pc == NOCOL ? -1 : slot_of[pc];
      int64_t s = (int64_t)t * ix->blk_bp, e = s + ix->blk_bp;
      if (e > (int64_t)ix->blk_len) e = (int64_t)ix->blk_len;
      if (m && tr[m - 1].slot == slot) tr[m - 1].X1 = e;
      else { tr[m].slot = slot; tr[m].X0 = s; tr[m].X1 = e; m++; }
    }
    rc = splice_track(R, cols, tr, m, g, &frags[k], &stat[k]);
    free(tr);
    if (!rc) {
      all_bp += stat[k].bp;
      fprintf(stderr, "splice: track %u: %u tracts, %u cuts at a shared node (%u strong, farthest %lld bp from "
                      "its switch), %u at a coordinate, %u next to a gap; %u fragments, %llu bp = %.3f x the reference\n",
              k + 1, stat[k].tracts, stat[k].ncut_node, stat[k].nstrong, (long long)stat[k].far_cut,
              stat[k].ncut_coord, stat[k].ngap, stat[k].nfrag, (unsigned long long)stat[k].bp,
              R->length ? (double)stat[k].bp / (double)R->length : 0.0);
    }
    if (!rc && stat[k].dropped_wl)
      fprintf(stderr, "splice: track %u: %llu W-lines without a reference node dropped (%llu bp)\n",
              k + 1, (unsigned long long)stat[k].dropped_wl, (unsigned long long)stat[k].dropped_bp);
  }
  if (!rc && T > 2)
    fprintf(stderr, "splice: the %u tracks hold %.2f haplotypes' worth of sequence in all\n", T,
            R->length ? (double)all_bp / (double)R->length : 0.0);
  char **kn = NULL; uint32_t nk = 0;
  if (!rc) {
    const char *csv = keep_csv ? keep_csv : ix->blk_ref;
    for (const char *p = csv; *p; ) {
      const char *e = strchr(p, ','); size_t n = e ? (size_t)(e - p) : strlen(p);
      if (n) { kn = realloc(kn, (nk + 1) * sizeof(char *)); kn[nk] = strndup(p, n); nk++; }
      if (!e) break;
      p = e + 1;
    }
    for (uint32_t i = 0; i < nk; i++) {
      int hit = 0;
      for (uint32_t w = 0; w < g->np; w++) if (g->pname[w] && !strcmp(g->pname[w], kn[i])) hit = 1;
      if (!hit) { fprintf(stderr, "--mosaic-keep: no walk of colour \"%s\"\n", kn[i]); rc = 1; }
    }
  }
  if (!rc && X) {
    uint64_t lw = 0, all = 0;
    for (uint32_t s = 0; s < nslot; s++)
      for (uint32_t wi = 0; wi < cols[s].nwl; wi++) if (cols[s].wl[wi].lazy) { lw++; all += cols[s].wl[wi].xw->nck; }
    if (X->err) {
      fprintf(stderr, "splice: a chunk of the GFA does not read as its .gfai says: the GFA changed after the .gfai"
                      " was written; delete the .gfai, or pass --no-gfai\n");
      rc = 1;
    } else
      fprintf(stderr, "splice: .gfai: %llu of the %llu chunks of %llu W-lines read, %.1f MB of steps\n",
              (unsigned long long)X->nload, (unsigned long long)all, (unsigned long long)lw, X->bload / 1e6);
  }
  // every consecutive pair of every output walk must be an edge of the graph
  if (!rc) {
    uint64_t bad = 0, steps = 0;
    for (uint32_t k = 0; k < T; k++)
      for (uint32_t f = 0; f < stat[k].nfrag; f++)
        for (uint32_t i = 0; i + 1 < frags[k][f].n; i++) {
          steps++;
          if (!has_edge(gfa, frags[k][f].v[i], frags[k][f].v[i + 1])) bad++;
        }
    if (bad) { fprintf(stderr, "splice: %llu of %llu consecutive steps are not edges of the graph\n",
                       (unsigned long long)bad, (unsigned long long)steps); rc = 1; }
    else fprintf(stderr, "splice: all %llu steps of the tracks are edges of the graph\n", (unsigned long long)steps);
  }
  if (!rc) {
    FILE *f = fopen(outfile, "w");
    if (!f) { fprintf(stderr, "splice: cannot write %s\n", outfile); rc = 1; }
    else {
      uint8_t *seen = calloc(((size_t)gfa->n_seg >> 3) + 1, 1);
      #define SEEN(z) (seen[(z) >> 3] & (1u << ((z) & 7)))
      #define MARK(z) (seen[(z) >> 3] |= (uint8_t)(1u << ((z) & 7)))
      for (uint32_t k = 0; k < T; k++)
        for (uint32_t fi = 0; fi < stat[k].nfrag; fi++)
          for (uint32_t i = 0; i < frags[k][fi].n; i++) MARK(frags[k][fi].v[i] >> 1);
      for (uint32_t w = 0; w < g->np; w++) {
        int keepw = 0;
        for (uint32_t i = 0; i < nk; i++) if (g->pname[w] && !strcmp(g->pname[w], kn[i])) keepw = 1;
        if (!keepw) continue;
        for (uint32_t s = 0; s < g->pn[w]; s++) MARK(g->pv[w][s] >> 1);
      }
      fprintf(f, "H\tVN:Z:1.1\n");
      uint64_t nseg = 0, nlink = 0;
      for (uint32_t z = 0; z < gfa->n_seg; z++)
        if (SEEN(z)) {
          fputs("S\t", f); fputs(g->sname[z], f); putc('\t', f); fputs(g->seq[z] ? g->seq[z] : "*", f); putc('\n', f);
          nseg++;
        }
      for (uint64_t a = 0; a < gfa->n_arc; a++) {
        const gfa_arc_t *r = &gfa->arc[a];
        if (r->del || r->comp) continue;
        uint32_t u = (uint32_t)(r->v_lv >> 33), v = r->w >> 1;
        if (!SEEN(u) || !SEEN(v)) continue;
        char ob[16]; snprintf(ob, sizeof ob, "%dM\n", r->ov < r->ow ? r->ov : r->ow);
        fputs("L\t", f); fputs(g->sname[u], f); putc('\t', f); putc("+-"[r->v_lv >> 32 & 1], f); putc('\t', f);
        fputs(g->sname[v], f); putc('\t', f); putc("+-"[r->w & 1], f); putc('\t', f); fputs(ob, f);
        nlink++;
      }
      for (uint32_t w = 0; w < g->np; w++) {
        int keepw = 0;
        for (uint32_t i = 0; i < nk; i++) if (g->pname[w] && !strcmp(g->pname[w], kn[i])) keepw = 1;
        if (!keepw) continue;
        const gfa_walk_t *wk = &gfa->walk[g->psrc[w]];
        fprintf(f, "W\t%s\t%d\t%s\t", wk->sample, wk->hap, gfa->sseq[wk->snid].name);
        if (wk->st >= 0 && wk->en >= 0) fprintf(f, "%lld\t%lld\t", (long long)wk->st, (long long)wk->en);
        else fprintf(f, "*\t*\t");
        for (uint32_t s = 0; s < g->pn[w]; s++)
          fprintf(f, "%c%s", "><"[g->pv[w][s] & 1], g->sname[g->pv[w][s] >> 1]);
        fputc('\n', f);
      }
      for (uint32_t k = 0; k < T; k++)
        for (uint32_t fi = 0; fi < stat[k].nfrag; fi++) {
          fprintf(f, "W\tmosaic\t%u\tfrag%u\t0\t%llu\t", k + 1, fi, (unsigned long long)frags[k][fi].bp);
          for (uint32_t i = 0; i < frags[k][fi].n; i++)
            fprintf(f, "%c%s", "><"[frags[k][fi].v[i] & 1], g->sname[frags[k][fi].v[i] >> 1]);
          fputc('\n', f);
        }
      fclose(f); free(seen);
      #undef SEEN
      #undef MARK
      fprintf(stderr, "splice: wrote %s: %llu segments, %llu links, %u kept colour(s) + %u tracks (%.1fs)\n",
              outfile, (unsigned long long)nseg, (unsigned long long)nlink, nk, T, now_s() - t0);
    }
  }
  for (uint32_t k = 0; k < T; k++) if (frags[k]) {
    for (uint32_t f = 0; f < stat[k].nfrag; f++) free(frags[k][f].v);
    free(frags[k]);
  }
  free(frags); free(stat);
  for (uint32_t i = 0; i < nk; i++) free(kn[i]);
  free(kn);
  for (uint32_t s = 0; s < nslot; s++) scol_free(&cols[s]);
  free(cols); free(slot_of); sref_free(R);
  return rc;
}

// A colour is excluded by its full name (HG002#1) or by its sample (HG002), as --colors does.
static int mosaic_listed(const char *csv, const char *name) {
  if (!csv) return 0;
  const char *h = strchr(name, '#');
  size_t L = strlen(name), Ls = h ? (size_t)(h - name) : L;
  for (const char *p = csv; *p; ) {
    const char *e = strchr(p, ','); size_t n = e ? (size_t)(e - p) : strlen(p);
    if ((n == L && !strncmp(p, name, n)) || (n == Ls && !strncmp(p, name, n))) return 1;
    if (!e) break;
    p = e + 1;
  }
  return 0;
}

static int mcmp_u32(const void *x, const void *y) {
  uint32_t a = *(const uint32_t *)x, b = *(const uint32_t *)y;
  return (a > b) - (a < b);
}

// --mosaic-temper: how much evidence the painting has against the N it is tempered to.
// chat = the median n_obs over the blocks with data (mean of the two middle ones at an
// even count) / N; w = min(1, 1/chat), 1 without data. `lt1` = share of the blocks with
// data that block tempering scales down (n_obs > N).
typedef struct { double N, median, chat, w, lt1; } temper_t;
static temper_t mosaic_temper_stats(const mobs_t *O, double N) {
  temper_t s = { N, 0.0, 0.0, 1.0, 0.0 };
  uint32_t *v = malloc((O->nB ? O->nB : 1) * sizeof(uint32_t)), n = 0, over = 0;
  for (uint32_t t = 0; t < O->nB; t++)
    if (O->nobs_blk[t]) { v[n++] = O->nobs_blk[t]; over += (double)O->nobs_blk[t] > N; }
  qsort(v, n, sizeof(uint32_t), mcmp_u32);
  if (n) s.median = n & 1 ? (double)v[n / 2] : 0.5 * ((double)v[n / 2 - 1] + (double)v[n / 2]);
  if (N > 0) s.chat = s.median / N;
  if (s.chat > 1.0) s.w = 1.0 / s.chat;
  if (n) s.lt1 = (double)over / (double)n;
  free(v);
  return s;
}

// What the painting looks like: tracts (maximal runs of one colour in one track),
// switches, the colours it uses, and how much of it sits on the reference colour.
static void mosaic_report(const pidx_t *ix, const mobs_t *O, const uint16_t *path) {
  uint32_t nB = O->nB, nc = ix->nc;
  uint32_t *tl = malloc((size_t)nB * 2 * sizeof(uint32_t)), ntl = 0;
  uint64_t sw[2] = {0, 0}, onref = 0, data_blocks = 0, same = 0;
  uint8_t *used = calloc(nc ? nc : 1, 1);
  for (int k = 0; k < 2; k++) {
    uint32_t run = 1;
    for (uint32_t t = 0; t < nB; t++) {
      used[path[2 * (size_t)t + k]] = 1;
      if (t && path[2 * (size_t)t + k] != path[2 * (size_t)(t - 1) + k]) { tl[ntl++] = run; run = 1; sw[k]++; }
      else if (t) run++;
    }
    tl[ntl++] = run;
  }
  qsort(tl, ntl, sizeof(uint32_t), mcmp_u32);
  uint64_t tot = 0, acc = 0, n50 = 0;
  for (uint32_t i = 0; i < ntl; i++) tot += tl[i];
  for (uint32_t i = ntl; i-- > 0 && !n50; ) { acc += tl[i]; if (2 * acc >= tot) n50 = tl[i]; }
  uint32_t nused = 0, nsamp = 0;
  for (uint32_t c = 0; c < nc; c++) if (used[c]) {
    nused++;
    const char *a = ix->cname[c] ? ix->cname[c] : "?", *h = strchr(a, '#');
    size_t La = h ? (size_t)(h - a) : strlen(a);
    int seen = 0;
    for (uint32_t d = 0; d < c && !seen; d++) if (used[d]) {
      const char *b = ix->cname[d] ? ix->cname[d] : "?", *hb = strchr(b, '#');
      size_t Lb = hb ? (size_t)(hb - b) : strlen(b);
      seen = La == Lb && !strncmp(a, b, La);
    }
    if (!seen) nsamp++;
  }
  for (uint32_t t = 0; t < nB; t++) {
    if (O->nobs_blk[t]) data_blocks++;
    if (path[2 * (size_t)t] == path[2 * (size_t)t + 1]) same++;
    for (int k = 0; k < 2; k++)
      if (ix->blk_ref && ix->cname[path[2 * (size_t)t + k]] &&
          !strcmp(ix->cname[path[2 * (size_t)t + k]], ix->blk_ref)) onref++;
  }
  fprintf(stderr, "mosaic: %u blocks of %u bp (%llu with data); %llu + %llu switches, "
                  "%u tracts, median %u / N50 %llu blocks\n",
          nB, ix->blk_bp, (unsigned long long)data_blocks, (unsigned long long)sw[0],
          (unsigned long long)sw[1], ntl, ntl ? tl[ntl / 2] : 0, (unsigned long long)n50);
  fprintf(stderr, "mosaic: %u colours (%u samples) used; both tracks the same colour in %.1f%% of "
                  "the blocks; %.1f%% of the track-blocks are on %s\n",
          nused, nsamp, 100.0 * (double)same / (double)nB,
          100.0 * (double)onref / (2.0 * (double)nB), ix->blk_ref ? ix->blk_ref : "?");
  // The most common pair. A donor that is in the panel paints as one pair throughout; a donor
  // that is not paints as a different pair every few blocks, and no pair covers much.
  uint32_t *pc = calloc((size_t)nc * nc, sizeof(uint32_t)), bestpair = 0; uint64_t ndat = 0;
  for (uint32_t t = 0; t < nB; t++) {
    if (!O->nobs_blk[t]) continue;
    uint32_t a = path[2 * (size_t)t], b = path[2 * (size_t)t + 1];
    if (a > b) { uint32_t x = a; a = b; b = x; }
    ndat++;
    if (++pc[(size_t)a * nc + b] > pc[bestpair]) bestpair = a * nc + b;
  }
  if (ndat)
    fprintf(stderr, "mosaic: most common pair %s + %s in %.1f%% of the %llu blocks with data\n",
            ix->cname[bestpair / nc], ix->cname[bestpair % nc],
            100.0 * (double)pc[bestpair] / (double)ndat, (unsigned long long)ndat);
  free(pc);
  free(tl); free(used);
}

// The painting, in the format the splice and the mapping driver of the study read: every
// block of the chromosome, the colour of each of the T tracks ("-" where a track holds none),
// and how many observations (or, for the cover, hits) the block held.
static void path_write(FILE *f, const pidx_t *ix, const mobs_t *O, const uint16_t *path, uint32_t T,
                       const char *hdr, const char *lastcol) {
  fprintf(f, "# %s\n", hdr);
  fprintf(f, "block\tstart\tend");
  for (uint32_t k = 1; k <= T; k++) fprintf(f, "\ttrack%u", k);
  fprintf(f, "\tn_obs\t%s\n", lastcol);
  for (uint32_t b = 0; b < O->nB; b++) {
    uint64_t s = (uint64_t)b * ix->blk_bp, e = s + ix->blk_bp;
    if (e > ix->blk_len) e = ix->blk_len;
    fprintf(f, "%u\t%llu\t%llu", b, (unsigned long long)s, (unsigned long long)e);
    for (uint32_t k = 0; k < T; k++) {
      uint16_t c = path[(size_t)b * T + k];
      fprintf(f, "\t%s", c == NOCOL ? "-" : ix->cname[c]);
    }
    fprintf(f, "\t%u\t%u\n", O->nobs_blk[b], O->ninf_blk[b]);
  }
}

static void mosaic_write_path(FILE *f, const pidx_t *ix, const mobs_t *O, const uint16_t *path,
                              double lam, double rho, double eps, int reads, const char *panel,
                              int temper, const temper_t *ts) {
  char hdr[1024];
  int k = snprintf(hdr, sizeof hdr, "panel=%s reads=%s obs=%llu B=%u coord=%s lambda=%g rho=%g eps=%g",
                   panel ? panel : "-", reads == MR_SHORT ? "short" : "hifi",
                   (unsigned long long)O->n_placed, ix->blk_bp, ix->blk_ref, lam, rho, eps);
  if (temper == TEMPER_GLOBAL && k > 0 && (size_t)k < sizeof hdr)
    snprintf(hdr + k, sizeof hdr - k, " temper=global N=%g chat=%.4f w=%.6f", ts->N, ts->chat, ts->w);
  else if (temper == TEMPER_BLOCK && k > 0 && (size_t)k < sizeof hdr)
    snprintf(hdr + k, sizeof hdr - k, " temper=block N=%g chat=%.4f tempered=%.4f", ts->N, ts->chat, ts->lt1);
  path_write(f, ix, O, path, 2, hdr, "n_inf");
}

// ---- cover: per window, the few colours whose UNION covers the reads' hits ------
//
// (--set-call cover.) A graph aligner needs every allele a read carries to be a path in the
// graph, and the graph is the union of the chosen haplotypes, so what counts for a window is how
// many of the alleles the reads show are present SOMEWHERE in the chosen set, not which colour
// carries them. An index hit is such an allele (a k-mer straddling a variant junction that the
// read has); every hit has a block, so evidence is additive and each window is decided alone,
// with no transition model:
//
//     covered(S) = sum over index records r hit n_r times with mask_r & S != 0 of n_r
//
// The always-present colours (the reference) cover what they carry first; the greedy then adds
// the colour covering the most uncovered hits. This is weighted maximum coverage (monotone,
// submodular), so gains never increase and the stop rule can read them: a colour is added only
// while it brings at least `min_gain` hits and `min_frac` of the window's hits. Nothing forces
// every read or record to be covered.
//
// The score is a screen, not a ranking: it does not transfer from whole-haplotype graphs to
// graphs of a few colours per window, whose alleles must also lie on one haplotype path.
#define NOCOL 0xFFFFu

typedef struct { uint32_t cmax; double min_gain, min_frac; int rule; } cover_par_t;
typedef struct { uint64_t hits, unc0, unc_end; } cover_stat_t;

// One window. `slot[i]`, `wt[i]` (i < n): the index records hit in it and their weights;
// `base`: the colours already in the graph; `allowed`: the colours that may be picked;
// `prev`: the colours picked in the previous window, which win exact ties (a pair of
// identical haplotypes then does not flip from window to window). Returns the picks in
// order, with the hits each one newly covered.
static uint32_t cover_window(const uint64_t *cmask, uint32_t nw, uint32_t nc, const uint32_t *slot,
                             const uint32_t *wt, uint32_t n, const uint64_t *base,
                             const uint8_t *allowed, const uint8_t *prev, const cover_par_t *P,
                             uint32_t *pick, uint64_t *pgain, cover_stat_t *st) {
  uint64_t *gain = calloc(nc ? nc : 1, sizeof(uint64_t));
  uint8_t *cov = malloc(n ? n : 1), *chosen = calloc(nc ? nc : 1, 1);
  uint64_t H = 0, U = 0;
  for (uint32_t i = 0; i < n; i++) {
    const uint64_t *m = cmask + (size_t)slot[i] * nw;
    H += wt[i];
    int c0 = 0;
    for (uint32_t z = 0; z < nw && !c0; z++) c0 = (m[z] & base[z]) != 0;
    cov[i] = (uint8_t)c0;
    if (c0) continue;
    U += wt[i];
    for (uint32_t z = 0; z < nw; z++)
      for (uint64_t b = m[z]; b; b &= b - 1) gain[(z << 6) + (uint32_t)__builtin_ctzll(b)] += wt[i];
  }
  st->hits = H; st->unc0 = U;
  uint32_t k = 0;
  while (k < P->cmax) {
    uint32_t best = nc; uint64_t bk = 0;                       // key = 2 * gain + the stay bonus
    for (uint32_t c = 0; c < nc; c++) {
      if (chosen[c] || (allowed && !allowed[c])) continue;
      uint64_t key = 2 * gain[c] + (prev && prev[c] ? 1u : 0u);
      if (best == nc || key > bk) { best = c; bk = key; }
    }
    if (best == nc || gain[best] == 0) break;
    uint64_t g = gain[best];
    if (P->rule == CR_GREEDY && ((double)g < P->min_gain || (double)g < P->min_frac * (double)H)) break;
    pick[k] = best; pgain[k] = g; chosen[best] = 1; k++;
    if (P->rule != CR_GREEDY) continue;
    for (uint32_t i = 0; i < n; i++) {
      if (cov[i]) continue;
      const uint64_t *m = cmask + (size_t)slot[i] * nw;
      if (!((m[best >> 6] >> (best & 63)) & 1)) continue;
      cov[i] = 1; U -= wt[i];
      for (uint32_t z = 0; z < nw; z++)
        for (uint64_t b = m[z]; b; b &= b - 1) gain[(z << 6) + (uint32_t)__builtin_ctzll(b)] -= wt[i];
    }
  }
  if (P->rule != CR_GREEDY) {
    U = 0;
    for (uint32_t i = 0; i < n; i++) {
      const uint64_t *m = cmask + (size_t)slot[i] * nw;
      int c = 0;
      for (uint32_t z = 0; z < nw && !c; z++) c = (m[z] & base[z]) != 0;
      for (uint32_t j = 0; j < k && !c; j++) c = (m[pick[j] >> 6] >> (pick[j] & 63)) & 1;
      if (!c) U += wt[i];
    }
  }
  st->unc_end = U;
  free(gain); free(cov); free(chosen);
  return k;
}

typedef struct {
  uint32_t nB, nwin, T, cmax;      // blocks, windows, tracks, picks per window
  uint32_t *wb;                    // [nwin+1]: the first block of each window (windows are runs of whole blocks)
  uint16_t *path;                  // [nB][T]: the colour of each track, NOCOL where it is empty
  uint32_t *nsel;                  // [nwin]
  uint16_t *sel;                   // [nwin][cmax]  the picks, in order
  uint64_t *gain;                  // [nwin][cmax]  the hits each newly covered
  cover_stat_t *st;                // [nwin]
  uint32_t wcap;                   // the cap of a record's weight
} cover_t;

static void cover_free(cover_t *C) {
  if (!C) return;
  free(C->wb); free(C->path); free(C->nsel); free(C->sel); free(C->gain); free(C->st); free(C);
}

// The windows: runs of whole blocks. Fixed: `W` blocks each. Automatic: a window is closed once it
// holds `target` hits (or `maxb` blocks), so the evidence a choice is made from does not depend
// on the coverage. A last window with less than half the target is merged into the one before.
static uint32_t *cover_windows(const uint32_t *hits_blk, uint32_t nB, uint32_t W, int automatic,
                               uint64_t target, uint32_t maxb, uint32_t *nwin_out) {
  uint32_t *wb = malloc(((size_t)nB + 2) * sizeof(uint32_t)); uint32_t nw = 0;
  wb[0] = 0;
  if (!automatic) {
    for (uint32_t b = W; b < nB; b += W) wb[++nw] = b;
    wb[++nw] = nB;
  } else {
    uint64_t cum = 0; uint32_t start = 0;
    for (uint32_t b = 0; b < nB; b++) {
      cum += hits_blk[b];
      if (cum >= target || b + 1 - start >= maxb) {
        if (b + 1 < nB) wb[++nw] = b + 1;
        cum = 0; start = b + 1;
      }
    }
    if (cum > 0 && cum * 2 < target && nw > 0) nw--;
    wb[++nw] = nB;
  }
  *nwin_out = nw;
  return wb;
}

// Every window of the chromosome, then the tracks: a colour picked in two neighbouring
// windows keeps its track, so a track is cut only where its colour really changes.
static cover_t *mosaic_cover(const pidx_t *ix, const uint16_t *hitcnt, const uint32_t *hits_blk, uint32_t nB,
                             const mosaic_opt_t *opt, const uint64_t *base, const uint8_t *allowed) {
  uint32_t nc = ix->nc, nw = ix->nw;
  uint32_t W = opt->cwin ? (opt->cwin + ix->blk_bp - 1) / ix->blk_bp : nB;
  if (W < 1) W = 1;
  if (W > nB) W = nB;
  uint32_t nwin = 0;
  uint32_t *wb = cover_windows(hits_blk, nB, W, opt->cauto, opt->chits, 100, &nwin);
  uint32_t *winof = malloc((size_t)nB * sizeof(uint32_t));
  for (uint32_t w = 0; w < nwin; w++) for (uint32_t b = wb[w]; b < wb[w + 1]; b++) winof[b] = w;
  size_t nslots = ix->km->mask + 1, nrec = 0;
  double sum = 0;
  for (size_t s = 0; s < nslots; s++)
    if (ix->km->st[s] == KH_IDX && hitcnt[s] && ix->blk[s] != BLK_NONE) { nrec++; sum += hitcnt[s]; }
  // A record hit hundreds of times (a collapsed repeat) must not decide a window alone.
  uint32_t wcap = nrec ? (uint32_t)(3.0 * sum / (double)nrec + 0.999) : 4;
  if (wcap < 4) wcap = 4;
  size_t *ws = calloc((size_t)nwin + 1, sizeof(size_t));
  for (size_t s = 0; s < nslots; s++)
    if (ix->km->st[s] == KH_IDX && hitcnt[s] && ix->blk[s] != BLK_NONE) {
      uint32_t b = ix->blk[s]; if (b >= nB) b = nB - 1;
      ws[winof[b] + 1]++;
    }
  for (uint32_t w = 0; w < nwin; w++) ws[w + 1] += ws[w];
  uint32_t *wslot = malloc((nrec ? nrec : 1) * sizeof(uint32_t)), *wwt = malloc((nrec ? nrec : 1) * sizeof(uint32_t));
  { size_t *fill = malloc(((size_t)nwin + 1) * sizeof(size_t));
    memcpy(fill, ws, ((size_t)nwin + 1) * sizeof(size_t));
    for (size_t s = 0; s < nslots; s++)
      if (ix->km->st[s] == KH_IDX && hitcnt[s] && ix->blk[s] != BLK_NONE) {
        uint32_t b = ix->blk[s]; if (b >= nB) b = nB - 1;
        size_t o = fill[winof[b]]++;
        wslot[o] = (uint32_t)s;
        wwt[o] = opt->cweight == CW_PRESENCE ? 1u : (hitcnt[s] < wcap ? hitcnt[s] : wcap);
      }
    free(fill); }

  cover_t *C = calloc(1, sizeof *C);
  C->nB = nB; C->nwin = nwin; C->wb = wb; C->cmax = opt->cmax; C->wcap = wcap;
  C->nsel = calloc(nwin, sizeof(uint32_t));
  C->sel = calloc((size_t)nwin * opt->cmax, sizeof(uint16_t));
  C->gain = calloc((size_t)nwin * opt->cmax, sizeof(uint64_t));
  C->st = calloc(nwin, sizeof(cover_stat_t));
  cover_par_t P = { opt->cmax, opt->cgain, opt->cfrac, opt->crule };
  uint8_t *prev = calloc(nc ? nc : 1, 1);
  uint32_t *pk = malloc(opt->cmax * sizeof(uint32_t)); uint64_t *pg = malloc(opt->cmax * sizeof(uint64_t));
  for (uint32_t w = 0; w < nwin; w++) {
    uint32_t n = (uint32_t)(ws[w + 1] - ws[w]);
    uint32_t k = cover_window(ix->cmask, nw, nc, wslot + ws[w], wwt + ws[w], n, base, allowed, prev, &P,
                              pk, pg, &C->st[w]);
    C->nsel[w] = k;
    memset(prev, 0, nc ? nc : 1);
    for (uint32_t j = 0; j < k; j++) {
      C->sel[(size_t)w * opt->cmax + j] = (uint16_t)pk[j];
      C->gain[(size_t)w * opt->cmax + j] = pg[j];
      prev[pk[j]] = 1;
    }
    if (k > C->T) C->T = k;
  }
  free(pk); free(pg); free(prev); free(ws); free(wslot); free(wwt);
  uint32_t T = C->T ? C->T : 1;
  C->T = T;
  C->path = malloc((size_t)nB * T * sizeof(uint16_t));
  for (size_t z = 0; z < (size_t)nB * T; z++) C->path[z] = NOCOL;
  uint16_t *slotcol = malloc(T * sizeof(uint16_t));
  for (uint32_t t = 0; t < T; t++) slotcol[t] = NOCOL;
  for (uint32_t w = 0; w < nwin; w++) {
    const uint16_t *S = C->sel + (size_t)w * opt->cmax; uint32_t k = C->nsel[w];
    uint8_t placed[64] = {0};                                  // k <= cmax <= 64
    for (uint32_t t = 0; t < T; t++) {
      int keep = 0;
      for (uint32_t j = 0; j < k && slotcol[t] != NOCOL; j++) if (S[j] == slotcol[t]) { keep = 1; placed[j] = 1; }
      if (!keep) slotcol[t] = NOCOL;
    }
    for (uint32_t j = 0; j < k; j++) {
      if (placed[j]) continue;
      uint32_t t = 0; while (t < T && slotcol[t] != NOCOL) t++;
      slotcol[t] = S[j];
    }
    for (uint32_t b = wb[w]; b < wb[w + 1]; b++)
      for (uint32_t t = 0; t < T; t++) C->path[(size_t)b * T + t] = slotcol[t];
  }
  free(slotcol); free(winof);
  return C;
}

// What the sets look like, and the gain curve: the share of the hits that no colour of
// the graph carries after 0, 1, 2, ... picks, which is the quantity the objective lowers.
static void cover_report(const pidx_t *ix, const mobs_t *O, const cover_t *C, const mosaic_opt_t *opt) {
  uint32_t nc = ix->nc;
  uint64_t H = 0, U0 = 0, Uend = 0, nwhit = 0, kept = 0, picks = 0;
  uint64_t *Uk = calloc((size_t)opt->cmax + 1, sizeof(uint64_t));
  uint32_t *hist = calloc((size_t)opt->cmax + 1, sizeof(uint32_t));
  uint8_t *used = calloc(nc ? nc : 1, 1);
  for (uint32_t w = 0; w < C->nwin; w++) {
    const cover_stat_t *s = &C->st[w];
    if (!s->hits) continue;
    nwhit++; H += s->hits; U0 += s->unc0; Uend += s->unc_end;
    uint64_t u = s->unc0;
    for (uint32_t j = 0; j <= opt->cmax; j++) {
      Uk[j] += u;
      if (j < C->nsel[w]) u -= C->gain[(size_t)w * opt->cmax + j];
    }
    hist[C->nsel[w]]++; picks += C->nsel[w];
    for (uint32_t j = 0; j < C->nsel[w]; j++) {
      uint16_t c = C->sel[(size_t)w * opt->cmax + j]; used[c] = 1;
      if (w) for (uint32_t i = 0; i < C->nsel[w - 1]; i++) if (C->sel[(size_t)(w - 1) * opt->cmax + i] == c) { kept++; break; }
    }
  }
  uint32_t nused = 0, nsamp = 0;
  for (uint32_t c = 0; c < nc; c++) if (used[c]) {
    nused++;
    const char *a = ix->cname[c], *h = strchr(a, '#'); size_t La = h ? (size_t)(h - a) : strlen(a);
    int seen = 0;
    for (uint32_t d = 0; d < c && !seen; d++) if (used[d]) {
      const char *b = ix->cname[d], *hb = strchr(b, '#'); size_t Lb = hb ? (size_t)(hb - b) : strlen(b);
      seen = La == Lb && !strncmp(a, b, La);
    }
    if (!seen) nsamp++;
  }
  uint32_t wmin = UINT32_MAX, wmax = 0;
  for (uint32_t w = 0; w < C->nwin; w++) {
    uint32_t l = C->wb[w + 1] - C->wb[w];
    if (l < wmin) wmin = l;
    if (l > wmax) wmax = l;
  }
  fprintf(stderr, "cover: %u windows (%.0f kb on average, %u to %u kb), %llu with hits; %u tracks; "
                  "%.2f colours per window (%llu picks)\n",
          C->nwin, (double)C->nB * ix->blk_bp / 1e3 / (double)C->nwin, wmin * ix->blk_bp / 1000, wmax * ix->blk_bp / 1000,
          (unsigned long long)nwhit, C->T, nwhit ? (double)picks / (double)nwhit : 0.0, (unsigned long long)picks);
  fprintf(stderr, "cover: colours per window:");
  for (uint32_t j = 0; j <= opt->cmax; j++) if (hist[j]) fprintf(stderr, " %u:%u", j, hist[j]);
  fprintf(stderr, "\ncover: hits carried by no colour of the graph, after 0,1,2,... picks:");
  for (uint32_t j = 0; j <= opt->cmax && j < 9; j++)
    fprintf(stderr, " %.2f%%", H ? 100.0 * (double)Uk[j] / (double)H : 0.0);
  fprintf(stderr, "\ncover: %.2f%% of the hits are carried by the reference alone; %.3f%% stay uncovered at the end\n",
          H ? 100.0 * (double)(H - U0) / (double)H : 0.0, H ? 100.0 * (double)Uend / (double)H : 0.0);
  fprintf(stderr, "cover: %u colours (%u samples) used; %.1f%% of the picks were also picked in the previous window\n",
          nused, nsamp, picks ? 100.0 * (double)kept / (double)picks : 0.0);
  (void)O;
  free(Uk); free(hist); free(used);
}

static void cover_write_sets(FILE *f, const pidx_t *ix, const cover_t *C) {
  fprintf(f, "window\tstart\tend\thits\tuncovered_by_base\trank\tcolour\tgain\tuncovered_after\n");
  for (uint32_t w = 0; w < C->nwin; w++) {
    uint64_t s = (uint64_t)C->wb[w] * ix->blk_bp, e = (uint64_t)C->wb[w + 1] * ix->blk_bp;
    if (e > ix->blk_len) e = ix->blk_len;
    uint64_t u = C->st[w].unc0;
    if (!C->nsel[w]) fprintf(f, "%u\t%llu\t%llu\t%llu\t%llu\t0\t-\t0\t%llu\n", w, (unsigned long long)s,
                             (unsigned long long)e, (unsigned long long)C->st[w].hits,
                             (unsigned long long)C->st[w].unc0, (unsigned long long)u);
    for (uint32_t j = 0; j < C->nsel[w]; j++) {
      uint64_t g = C->gain[(size_t)w * C->cmax + j];
      u -= g;
      fprintf(f, "%u\t%llu\t%llu\t%llu\t%llu\t%u\t%s\t%llu\t%llu\n", w, (unsigned long long)s,
              (unsigned long long)e, (unsigned long long)C->st[w].hits, (unsigned long long)C->st[w].unc0, j + 1,
              ix->cname[C->sel[(size_t)w * C->cmax + j]], (unsigned long long)g, (unsigned long long)u);
    }
  }
}

// The colours that are in every output graph (--mosaic-keep, default the coordinate
// reference), as a mask: they carry what they carry before anything is chosen.
static uint64_t *mosaic_base(const pidx_t *ix, const mosaic_opt_t *opt, uint32_t *nbase_out) {
  uint64_t *base = calloc(ix->nw, sizeof(uint64_t)); uint32_t nb = 0;
  const char *keep = opt->keep ? opt->keep : ix->blk_ref;
  for (uint32_t c = 0; c < ix->nc; c++)
    if (mosaic_listed(keep, ix->cname[c] ? ix->cname[c] : "?")) { cm_set(base, c); nb++; }
  if (nbase_out) *nbase_out = nb;
  return base;
}

// How much of an INDEPENDENT read set does the graph carry? The graph of a painting is its
// reference colours plus, in each block, the colours of its tracks. A hit is carried if one of
// those colours has the k-mer; the share of uncarried hits predicts the edits a graph aligner
// makes. Pass reads NOT used to choose the colours, or small windows look better than they are.
static void mosaic_check1(const pidx_t *ix, const char *reads, const char *mate,
                          const uint16_t *path, uint32_t T, uint32_t nB, const uint64_t *base, int stray,
                          const char *tag) {
  mobs_t O; memset(&O, 0, sizeof O);
  O.nc = ix->nc; O.nB = nB; O.want_m = 1;
  O.nobs_blk = calloc(nB, sizeof(uint32_t)); O.ninf_blk = calloc(nB, sizeof(uint32_t));
  uint16_t *hc = calloc(ix->km->mask + 1, sizeof(uint16_t));
  if (mosaic_collect(ix, reads, mate, &O, hc, stray)) { mobs_free(&O); free(hc); return; }
  // The best single colour of the graph, per read: how many of the read's hits does the colour
  // of the graph that fits it best NOT carry? A read whose alleles lie on one haplotype path of
  // the graph is aligned without a gap; one whose alleles are split between two colours is
  // not, however many of them the union carries. (Reads that no colour tells apart add 0.)
  double mm_graph = 0, mm_base = 0;
  for (size_t i = 0; i < O.n; i++) {
    const uint16_t *row = O.m + i * (size_t)ix->nc; uint32_t b = O.blk[i];
    uint32_t bg = UINT32_MAX, bb = UINT32_MAX;
    for (uint32_t c = 0; c < ix->nc; c++)
      if ((base[c >> 6] >> (c & 63)) & 1) { if (row[c] < bb) bb = row[c]; }
    bg = bb;
    for (uint32_t k = 0; k < T; k++) {
      uint16_t c = path[(size_t)b * T + k];
      if (c != NOCOL && row[c] < bg) bg = row[c];
    }
    mm_base += bb == UINT32_MAX ? 0 : bb; mm_graph += bg == UINT32_MAX ? 0 : bg;
  }
  uint64_t tot = 0, miss_base = 0, miss = 0;
  uint32_t nw = ix->nw;
  for (size_t sl = 0; sl <= ix->km->mask; sl++) {
    if (ix->km->st[sl] != KH_IDX || !hc[sl] || ix->blk[sl] == BLK_NONE) continue;
    uint32_t b = ix->blk[sl]; if (b >= nB) b = nB - 1;
    const uint64_t *m = ix->cmask + sl * nw;
    tot += hc[sl];
    int c0 = 0;
    for (uint32_t z = 0; z < nw && !c0; z++) c0 = (m[z] & base[z]) != 0;
    if (c0) continue;
    miss_base += hc[sl];
    int c1 = 0;
    for (uint32_t k = 0; k < T && !c1; k++) {
      uint16_t c = path[(size_t)b * T + k];
      c1 = c != NOCOL && ((m[c >> 6] >> (c & 63)) & 1);
    }
    if (!c1) miss += hc[sl];
  }
  double rd = O.n_placed ? (double)O.n_placed : 1.0;
  fprintf(stderr, "check[%s]: %s: %llu reads, %llu placed hits; carried by no colour of the graph: "
                  "reference alone %.3f%% (%.2f per read), this graph %.3f%% (%.2f per read); "
                  "best single colour misses %.2f hits per read (reference alone %.2f)\n",
          tag, reads, (unsigned long long)O.n_obs, (unsigned long long)tot,
          tot ? 100.0 * (double)miss_base / (double)tot : 0.0, (double)miss_base / rd,
          tot ? 100.0 * (double)miss / (double)tot : 0.0, (double)miss / rd, mm_graph / rd, mm_base / rd);
  mobs_free(&O); free(hc);
}

static void mosaic_check(const pidx_t *ix, const mosaic_opt_t *opt, const uint16_t *path, uint32_t T,
                         uint32_t nB, const uint64_t *base, int stray, const char *tag) {
  char *rl = strdup(opt->check), *ml = opt->check_mate ? strdup(opt->check_mate) : NULL;
  char *rs = NULL, *ms = NULL;
  char *r = strtok_r(rl, ",", &rs), *m = ml ? strtok_r(ml, ",", &ms) : NULL;
  for (; r; r = strtok_r(NULL, ",", &rs), m = ml ? strtok_r(NULL, ",", &ms) : NULL) {
    if (ml && !m) { fprintf(stderr, "check: --mosaic-check-mate has fewer files than --mosaic-check\n"); break; }
    mosaic_check1(ix, r, m, path, T, nB, base, stray, tag);
  }
  free(rl); free(ml);
}

static int splice_out(const mosaic_opt_t *opt, const gfa_t *gfa, const bg_t *g, const pidx_t *ix,
                      const uint16_t *path, uint32_t nB, uint32_t T);

static int run_mosaic(const pidx_t *ix, const char *reads, const mosaic_opt_t *opt,
                      const gfa_t *gfa, const bg_t *g) {
  if (!ix->blk) {
    fprintf(stderr, "--set-call mosaic needs positions in the index: build it with"
                    " --mosaic-ref COLOUR\n");
    return 1;
  }
  uint32_t nc = ix->nc;
  mobs_t O; memset(&O, 0, sizeof O);
  O.nc = nc; O.nB = (uint32_t)((ix->blk_len + ix->blk_bp - 1) / ix->blk_bp);
  if (!O.nB) { fprintf(stderr, "mosaic: the index covers no reference\n"); return 1; }
  O.nobs_blk = calloc(O.nB, sizeof(uint32_t)); O.ninf_blk = calloc(O.nB, sizeof(uint32_t));
  double t0 = now_s();
  if (mosaic_collect(ix, reads, opt->mate, &O, NULL, 0)) { mobs_free(&O); return 1; }
  fprintf(stderr, "mosaic: %llu %s: %llu without a hit, %llu with hits but no position, %llu placed, "
                  "%llu informative (%.1fs)\n",
          (unsigned long long)O.n_obs, opt->mate ? "pairs" : "reads",
          (unsigned long long)O.n_nohit, (unsigned long long)O.n_unplaced,
          (unsigned long long)O.n_placed, (unsigned long long)O.n_inf, now_s() - t0);
  if (!O.n_inf) { fprintf(stderr, "mosaic: no informative observation, nothing to paint\n"); mobs_free(&O); return 1; }

  int rt = opt->reads;
  if (rt == MR_AUTO) {
    double mean = (double)O.bases / ((double)O.n_obs * (opt->mate ? 2.0 : 1.0));
    rt = mean >= 1000.0 ? MR_HIFI : MR_SHORT;
    fprintf(stderr, "mosaic: mean read length %.0f bp -> %s preset\n", mean, rt == MR_HIFI ? "hifi" : "short");
  }
  double lam, rho, eps;
  mosaic_preset(rt, &lam, &rho, &eps);
  if (opt->lambda >= 0) lam = opt->lambda;
  if (opt->rho >= 0)    rho = opt->rho;
  if (opt->eps >= 0)    eps = opt->eps;
  if (!(lam > 0 && rho > 0 && rho < 1 && eps > 0 && eps < 1)) {
    fprintf(stderr, "mosaic: need lambda > 0, 0 < rho < 1, 0 < eps < 1\n"); mobs_free(&O); return 1;
  }
  uint32_t *cs = malloc((nc ? nc : 1) * sizeof(uint32_t)), Cs = 0;
  for (uint32_t c = 0; c < nc; c++)
    if (!mosaic_listed(opt->exclude, ix->cname[c] ? ix->cname[c] : "?")) cs[Cs++] = c;
  if (!Cs) { fprintf(stderr, "mosaic: --mosaic-exclude leaves no colour\n"); free(cs); mobs_free(&O); return 1; }
  double need = (double)O.nB * (double)Cs * (double)Cs;
  if (need > 3.2e10) {
    fprintf(stderr, "mosaic: %u blocks x %u^2 colour pairs needs %.0f GB of back-pointers; "
                    "narrow the panel (--mosaic-exclude) or use a larger --mosaic-block\n",
            O.nB, Cs, need / 1e9);
    free(cs); mobs_free(&O); return 1;
  }
  fprintf(stderr, "mosaic: lambda %g rho %g eps %g, %u colours%s, %u blocks (%.2f GB of back-pointers)\n",
          lam, rho, eps, Cs, Cs < nc ? " allowed" : "", O.nB, need / 1e9);
  // --mosaic-temper: the short preset is tempered globally unless told otherwise; N scales
  // with the index's block size (n_obs grows with it).
  int temper = opt->temper == TEMPER_AUTO ? (rt == MR_SHORT ? TEMPER_GLOBAL : TEMPER_OFF) : opt->temper;
  double temper_n = opt->temper_n > 0 ? opt->temper_n : TEMPER_N_10KB * (double)ix->blk_bp / 10000.0;
  temper_t ts = mosaic_temper_stats(&O, temper_n);
  if (temper == TEMPER_GLOBAL)
    fprintf(stderr, "mosaic: temper global, N %g, median %.1f observations per block -> w %.6f\n",
            ts.N, ts.median, ts.w);
  else if (temper == TEMPER_BLOCK)
    fprintf(stderr, "mosaic: temper block, N %g, median %.1f observations per block, %.1f%% of the "
                    "blocks with data scaled down\n", ts.N, ts.median, 100.0 * ts.lt1);
  t0 = now_s();
  double score = 0;
  uint16_t *path = mosaic_paint(&O, cs, Cs, lam, rho, eps, temper, temper_n, ts.w, &score);
  fprintf(stderr, "mosaic: Viterbi log-score %.1f (%.1fs)\n", score, now_s() - t0);
  mosaic_report(ix, &O, path);
  if (opt->check) {
    uint64_t *base = mosaic_base(ix, opt, NULL);
    mosaic_check(ix, opt, path, 2, O.nB, base, rt == MR_HIFI ? 5 : 2, "mosaic");
    free(base);
  }

  FILE *out = stdout; char pp[4096];
  if (opt->out) {
    snprintf(pp, sizeof pp, "%s.path.tsv", opt->out);
    out = fopen(pp, "w");
    if (!out) { fprintf(stderr, "mosaic: cannot write %s\n", pp); free(cs); free(path); mobs_free(&O); return 1; }
  }
  mosaic_write_path(out, ix, &O, path, lam, rho, eps, rt, opt->panel, temper, &ts);
  if (opt->out) { fclose(out); fprintf(stderr, "mosaic: painting written to %s\n", pp); }
  int rcs = 0;
  if (opt->gfa_out) rcs = splice_out(opt, gfa, g, ix, path, O.nB, 2);
  free(cs); free(path); mobs_free(&O);
  return rcs;
}


// The mean length of the first reads decides the preset when --mosaic-reads is not given.
static int guess_reads(const char *path, int paired) {
  if (paired) return MR_SHORT;
  rdr_t r;
  if (rdr_open(&r, path)) return MR_SHORT;
  uint64_t bases = 0, n = 0;
  while (n < 200 && rdr_next(&r) == 1) { bases += r.slen; n++; }
  rdr_close(&r);
  return n && bases / n >= 1000 ? MR_HIFI : MR_SHORT;
}

static int run_cover(const pidx_t *ix, const char *reads, const mosaic_opt_t *opt,
                     const gfa_t *gfa, const bg_t *g) {
  if (!ix->blk) {
    fprintf(stderr, "--set-call cover needs positions in the index: build it with"
                    " --mosaic-ref COLOUR\n");
    return 1;
  }
  uint32_t nc = ix->nc, nw = ix->nw;
  if (opt->cmax < 1 || opt->cmax > 64) { fprintf(stderr, "--cover-max: 1 to 64\n"); return 1; }
  if (opt->cfrac < 0 || opt->cfrac >= 1 || opt->cgain < 0) {
    fprintf(stderr, "--cover-min-frac must be in [0, 1) and --cover-min-gain >= 0\n"); return 1;
  }
  int rt = opt->reads == MR_AUTO ? guess_reads(reads, opt->mate != NULL) : opt->reads;
  int stray = opt->cstray >= 0 ? opt->cstray : (rt == MR_HIFI ? 5 : 2);
  mobs_t O; memset(&O, 0, sizeof O);
  O.nc = nc; O.nB = (uint32_t)((ix->blk_len + ix->blk_bp - 1) / ix->blk_bp);
  if (!O.nB) { fprintf(stderr, "cover: the index covers no reference\n"); return 1; }
  O.nobs_blk = calloc(O.nB, sizeof(uint32_t)); O.ninf_blk = calloc(O.nB, sizeof(uint32_t));
  uint16_t *hitcnt = calloc(ix->km->mask + 1, sizeof(uint16_t));
  double t0 = now_s();
  if (opt->threads > 1 ? mosaic_collect_par(ix, reads, opt->mate, &O, hitcnt, stray, opt->threads)
                       : mosaic_collect(ix, reads, opt->mate, &O, hitcnt, stray)) { mobs_free(&O); free(hitcnt); return 1; }
  fprintf(stderr, "cover: %llu %s (%s reads): %llu without a hit, %llu with hits but no position, %llu placed"
                  "; hits kept within %d blocks of their read's median (%.1fs)\n",
          (unsigned long long)O.n_obs, opt->mate ? "pairs" : "reads", rt == MR_HIFI ? "hifi" : "short",
          (unsigned long long)O.n_nohit, (unsigned long long)O.n_unplaced, (unsigned long long)O.n_placed,
          stray, now_s() - t0);
  if (!O.n_placed) { fprintf(stderr, "cover: no placed read, nothing to cover\n"); mobs_free(&O); free(hitcnt); return 1; }

  uint32_t nbase = 0;
  uint64_t *base = mosaic_base(ix, opt, &nbase);
  const char *keep = opt->keep ? opt->keep : ix->blk_ref;
  uint8_t *allowed = malloc(nc ? nc : 1);
  for (uint32_t c = 0; c < nc; c++) allowed[c] = !mosaic_listed(opt->exclude, ix->cname[c] ? ix->cname[c] : "?");
  (void)nw;
  if (!nbase) { fprintf(stderr, "cover: --mosaic-keep names no colour of the index\n"); mobs_free(&O); free(hitcnt); free(base); free(allowed); return 1; }
  t0 = now_s();
  cover_t *C = mosaic_cover(ix, hitcnt, O.ninf_blk, O.nB, opt, base, allowed);
  fprintf(stderr, "cover: %u colour(s) already in the graph (%s); at most %u per window, stop below %g hits or "
                  "%.3g of the window's hits; weight %s (capped at %u), rule %s (%.1fs)\n",
          nbase, keep, opt->cmax, opt->cgain, opt->cfrac, opt->cweight == CW_PRESENCE ? "presence" : "count",
          C->wcap, opt->crule == CR_SUPPORT ? "support (control)" : "greedy union", now_s() - t0);
  cover_report(ix, &O, C, opt);
  if (opt->check) {
    mosaic_check(ix, opt, C->path, C->T, O.nB, base, stray, "cover");
  }

  char hdr[1024];
  char win[64];
  if (opt->cauto) snprintf(win, sizeof win, "auto:%llu", (unsigned long long)opt->chits);
  else snprintf(win, sizeof win, "%u", opt->cwin ? opt->cwin : (uint32_t)(ix->blk_len));
  snprintf(hdr, sizeof hdr, "panel=%s reads=%s obs=%llu B=%u coord=%s solver=cover window=%s max=%u minfrac=%g "
                            "mingain=%g stray=%d rule=%s weight=%s",
           opt->panel ? opt->panel : "-", rt == MR_SHORT ? "short" : "hifi", (unsigned long long)O.n_placed,
           ix->blk_bp, ix->blk_ref, win, opt->cmax, opt->cfrac, opt->cgain, stray,
           opt->crule == CR_SUPPORT ? "support" : "greedy", opt->cweight == CW_PRESENCE ? "presence" : "count");
  FILE *out = stdout; char pp[4096]; int rcs = 0;
  if (opt->out) {
    snprintf(pp, sizeof pp, "%s.path.tsv", opt->out);
    out = fopen(pp, "w");
    if (!out) { fprintf(stderr, "cover: cannot write %s\n", pp); rcs = 1; }
  }
  if (!rcs) {
    path_write(out, ix, &O, C->path, C->T, hdr, "n_hits");
    if (opt->out) {
      fclose(out); fprintf(stderr, "cover: painting written to %s\n", pp);
      snprintf(pp, sizeof pp, "%s.sets.tsv", opt->out);
      FILE *sf = fopen(pp, "w");
      if (sf) { cover_write_sets(sf, ix, C); fclose(sf); fprintf(stderr, "cover: the picks of every window written to %s\n", pp); }
    }
    if (opt->gfa_out) rcs = splice_out(opt, gfa, g, ix, C->path, O.nB, C->T);
  }
  cover_free(C); mobs_free(&O); free(hitcnt); free(base); free(allowed);
  return rcs;
}


// --mosaic-score: score a painting that is already on disk (any number of tracks) on
// --mosaic-check reads, with the index only: no model is run.
static int run_score_file(const pidx_t *ix, const char *painting, const mosaic_opt_t *opt) {
  if (!ix->blk) { fprintf(stderr, "--mosaic-score needs an index with positions\n"); return 1; }
  FILE *f = fopen(painting, "r");
  if (!f) { fprintf(stderr, "--mosaic-score: cannot read %s\n", painting); return 1; }
  char *line = NULL; size_t cap = 0; ssize_t L;
  uint32_t T = 0, nB = 0, rcap = 1024; uint16_t *path = NULL; int rc = 0;
  while ((L = getline(&line, &cap, f)) > 0) {
    while (L && (line[L - 1] == '\n' || line[L - 1] == '\r')) line[--L] = 0;
    if (line[0] == '#' || !L) continue;
    if (!strncmp(line, "block\t", 6)) {
      T = 0;
      for (char *t = line; t; ) { if (!strncmp(t, "track", 5)) T++; char *e = strchr(t, '\t'); t = e ? e + 1 : NULL; }
      continue;
    }
    if (!T) T = 2;
    if (T > 64) { rc = 1; break; }
    if (!path) path = malloc((size_t)rcap * T * sizeof(uint16_t));
    char *fld[3 + 64 + 4]; uint32_t nf = 0;
    for (char *t = line; nf < 3 + T + 2; ) { fld[nf++] = t; char *e = strchr(t, '\t'); if (!e) break; *e = 0; t = e + 1; }
    if (nf < 3 + T) { rc = 1; break; }
    if (nB == rcap) { rcap *= 2; path = realloc(path, (size_t)rcap * T * sizeof(uint16_t)); }
    for (uint32_t k = 0; k < T; k++) {
      if (!strcmp(fld[3 + k], "-")) { path[(size_t)nB * T + k] = NOCOL; continue; }
      uint32_t c = 0;
      while (c < ix->nc && strcmp(ix->cname[c], fld[3 + k])) c++;
      if (c == ix->nc) {                                  // SAMPLE#HAP for a lone haplotype, as in --mosaic-splice
        const char *h = strchr(fld[3 + k], '#');
        if (h) { c = 0; while (c < ix->nc && (strncmp(ix->cname[c], fld[3 + k], (size_t)(h - fld[3 + k])) ||
                                              ix->cname[c][h - fld[3 + k]])) c++; }
      }
      if (c == ix->nc) { fprintf(stderr, "--mosaic-score: colour %s is not in the index\n", fld[3 + k]); rc = 1; break; }
      path[(size_t)nB * T + k] = (uint16_t)c;
    }
    if (rc) break;
    nB++;
  }
  free(line); fclose(f);
  uint32_t nBix = (uint32_t)((ix->blk_len + ix->blk_bp - 1) / ix->blk_bp);
  if (!rc && (nB != nBix || !T)) {
    fprintf(stderr, "--mosaic-score: %s has %u blocks, the index %u\n", painting, nB, nBix); rc = 1;
  }
  if (!rc) {
    uint64_t *base = mosaic_base(ix, opt, NULL);
    mosaic_check(ix, opt, path, T, nB, base, opt->cstray >= 0 ? opt->cstray : 5, "score");
    free(base);
  }
  free(path);
  return rc;
}


static int cmp_u32(const void *a, const void *b) {
  uint32_t x = *(const uint32_t *)a, y = *(const uint32_t *)b;
  return (x > y) - (x < y);
}

// path-union adjacency (G_P): edges are exactly the consecutive pairs of the
// walks, deduped per vertex. L-lines that no walk crosses are excluded, so the
// dominator pass never plants a guardian at a merge no sample reaches. Needs
// g->n / g->np / g->pv / g->pn already set.
static void build_padj(bg_t *g) {
  g->adjn = calloc(g->n, sizeof(uint32_t));
  g->adj  = calloc(g->n, sizeof(uint32_t *));
  uint32_t *cap = calloc(g->n, sizeof(uint32_t));
  for (uint32_t k = 0; k < g->np; k++)
    for (uint32_t s = 1; s < g->pn[k]; s++) cap[g->pv[k][s - 1]]++;
  for (uint32_t v = 0; v < g->n; v++)
    if (cap[v]) g->adj[v] = malloc(cap[v] * sizeof(uint32_t));
  for (uint32_t k = 0; k < g->np; k++)
    for (uint32_t s = 1; s < g->pn[k]; s++) {
      uint32_t u = g->pv[k][s - 1];
      g->adj[u][g->adjn[u]++] = g->pv[k][s];
    }
  for (uint32_t v = 0; v < g->n; v++) {
    if (g->adjn[v] > 1) {
      qsort(g->adj[v], g->adjn[v], sizeof(uint32_t), cmp_u32);
      uint32_t w = 1;
      for (uint32_t r = 1; r < g->adjn[v]; r++)
        if (g->adj[v][r] != g->adj[v][w - 1]) g->adj[v][w++] = g->adj[v][r];
      g->adjn[v] = w;
    }
    if (g->adjn[v] < cap[v])
      g->adj[v] = g->adjn[v]
                    ? realloc(g->adj[v], g->adjn[v] * sizeof(uint32_t))
                    : (free(g->adj[v]), NULL);
  }
  free(cap);
}


// wl_t: the --colors whitelist, read once and matched against the color names
// bg_build has already computed over the WHOLE graph.
typedef struct { char **e; uint32_t n; } wl_t;

void wl_free(wl_t *w) {
  if (!w) return;
  for (uint32_t i = 0; i < w->n; i++) free(w->e[i]);
  free(w->e); free(w);
}

wl_t *wl_read(const char *path) {
  FILE *f = fopen(path, "r");
  if (!f) { fprintf(stderr, "--colors: cannot read %s\n", path); return NULL; }
  wl_t *w = calloc(1, sizeof(wl_t));
  uint32_t cap = 64; w->e = malloc(cap * sizeof(char *));
  char line[1024];
  while (fgets(line, sizeof line, f)) {
    char *a = line, *b;
    while (*a == ' ' || *a == '\t') a++;
    b = a + strlen(a);
    while (b > a && (b[-1] == '\n' || b[-1] == '\r' || b[-1] == ' ' || b[-1] == '\t')) b--;
    *b = 0;
    if (!*a || *a == '#') continue;         // blank, or a comment: no color starts with '#'
    if (w->n == cap) { cap *= 2; w->e = realloc(w->e, cap * sizeof(char *)); }
    w->e[w->n++] = strdup(a);
  }
  fclose(f);
  return w;
}

// A whitelist entry matches a color name either exactly, or -- when the entry
// carries no '#' -- as the sample of a "sample#hap" name. So HG01109 selects
// HG01109#1 and HG01109#2 as well as a bare HG01109, while HG01109#1 selects
// only that one. Prefix matching is deliberately NOT done: it would let a typo
// silently widen the selection.
static int name_hit(const char *e, const char *name) {
  if (!strcmp(e, name)) return 1;
  size_t L = strlen(e);
  return !strchr(e, '#') && !strncmp(e, name, L) && name[L] == '#';
}

static int wl_hit(const wl_t *w, const char *name) {
  for (uint32_t i = 0; i < w->n; i++) if (name_hit(w->e[i], name)) return 1;
  return 0;
}


// ===========================================================================
// --split-colors: partition the panel and build one index per group, reading the GFA once.
//
// Two caps, because either alone is insufficient. The colour cap keeps a group at nw=1 (one
// 64-bit mask word per k-mer); the walk cap is what actually bounds RSS and is the binding one
// on big chromosomes.
//
// Whole samples stay together: splitting a diploid sample across groups would put its two
// haplotypes in different indexes, where neither can act as the other's near neighbour.
//
// The group lists name SAMPLES, and the whitelist matcher expands a '#'-free entry to every
// haplotype of that sample. --split-shared names colours (or samples) replicated into every group.
// ===========================================================================
typedef struct { char **name; uint64_t *bp; uint32_t *nhap; uint32_t n; } smp_t;

static void smp_free(smp_t *s) {
  if (!s) return;
  for (uint32_t i = 0; i < s->n; i++) free(s->name[i]);
  free(s->name); free(s->bp); free(s->nhap); free(s);
}

// Samples of the panel with their walk span and haplotype count, straight off
// the W lines. No graph is built: the partition only needs these three columns.
static smp_t *scan_samples(gfa_walk_t *src, uint32_t ns) {
  smp_t *s = calloc(1, sizeof(smp_t));
  s->name = malloc((ns ? ns : 1) * sizeof(char *));
  s->bp   = calloc(ns ? ns : 1, sizeof(uint64_t));
  s->nhap = calloc(ns ? ns : 1, sizeof(uint32_t));
  int32_t *seen_hap = malloc((ns ? ns : 1) * 8 * sizeof(int32_t));   // <=8 haps/sample
  for (uint32_t k = 0; k < ns; k++) {
    const char *nm = src[k].sample ? src[k].sample : "?";
    uint32_t i = 0; while (i < s->n && strcmp(s->name[i], nm)) i++;
    if (i == s->n) { s->name[i] = strdup(nm); s->n++; }
    s->bp[i] += (uint64_t)(src[k].en - src[k].st);
    int dup = 0;
    for (uint32_t h = 0; h < s->nhap[i] && h < 8; h++)
      if (seen_hap[i * 8 + h] == src[k].hap) { dup = 1; break; }
    if (!dup && s->nhap[i] < 8) seen_hap[i * 8 + s->nhap[i]++] = src[k].hap;
  }
  free(seen_hap);
  return s;
}

static wl_t *shared_list(const char *csv) {
  wl_t *w = calloc(1, sizeof(wl_t));
  if (!csv || !*csv) { w->e = malloc(sizeof(char *)); return w; }
  uint32_t cap = 8; w->e = malloc(cap * sizeof(char *));
  const char *p = csv;
  while (*p) {
    const char *q = strchr(p, ',');
    size_t L = q ? (size_t)(q - p) : strlen(p);
    if (L) {
      if (w->n == cap) { cap *= 2; w->e = realloc(w->e, cap * sizeof(char *)); }
      w->e[w->n] = malloc(L + 1); memcpy(w->e[w->n], p, L); w->e[w->n][L] = 0; w->n++;
    }
    if (!q) break;
    p = q + 1;
  }
  return w;
}

static int is_shared(const wl_t *sh, const char *sample) {
  for (uint32_t i = 0; i < sh->n; i++) {
    const char *e = sh->e[i];
    if (!strcmp(e, sample)) return 1;
    if (!strncmp(e, sample, strlen(sample)) && e[strlen(sample)] == '#') return 1;
  }
  return 0;
}

static int bg_quiet = 0;                               // the splice builds with its own whitelist: no --colors report
static int bg_nopadj = 0;                              // the splice reads walks only: no path-union adjacency (pebbling's)
bg_t *bg_build(gfa_t *gfa, int walks, int full_graph, int colsep,
               const wl_t *colors) {
  bg_t *g = calloc(1, sizeof(bg_t));
  g->n = gfa_n_vtx(gfa);
  gfa_walk_t *src = walks ? gfa->walk : gfa->path;
  uint32_t ns = walks ? gfa->n_walk : gfa->n_path;

  // borrowed seg views, indexed by SEGMENT id (v>>1); not owned.
  g->sname = malloc(gfa->n_seg * sizeof(char *));
  g->seq   = malloc(gfa->n_seg * sizeof(char *));
  g->slen  = malloc(gfa->n_seg * sizeof(uint32_t));
  for (uint32_t i = 0; i < gfa->n_seg; i++) {
    g->sname[i] = gfa->seg[i].name;
    g->seq[i]   = gfa->seg[i].seq;
    g->slen[i]  = gfa->seg[i].seq ? (uint32_t)strlen(gfa->seg[i].seq) : 0;
  }
  // A color is a (sample, haplotype) PAIR. gfatools interns only the sample, so naming a walk
  // after it alone would merge the two haplotypes of a diploid sample. The haplotype is appended
  // only where the sample is ambiguous by itself, so panels that already encode it in the sample
  // name keep their names. Names are INTERNED: colors are grouped by pointer identity downstream.
  g->pname = ns ? malloc(ns * sizeof(char *)) : NULL;
  g->cpool = ns ? malloc(ns * sizeof(char *)) : NULL;
  g->ncpool = 0;
  if (ns) {
    const char **usm = malloc(ns * sizeof(char *));    // distinct sample pointers
    int32_t *uhap = malloc(ns * sizeof(int32_t));      // the first hap seen for it
    uint8_t *usplit = calloc(ns, 1);                   // ...and whether others differ
    uint32_t nu = 0;
    for (uint32_t k = 0; k < ns; k++) {
      uint32_t u = 0; while (u < nu && usm[u] != src[k].sample) u++;
      if (u == nu) { usm[nu] = src[k].sample; uhap[nu] = src[k].hap; nu++; }
      else if (uhap[u] != src[k].hap) usplit[u] = 1;
    }
    // With `colsep`, the color is the name truncated at its first occurrence (a P-line graph names
    // paths per CONTIG). Two interned names can then land on one color, so the pool is keyed by
    // the resulting STRING rather than by the source pointer.
    char *nmbuf = NULL;
    for (uint32_t k = 0; k < ns; k++) {
      const char *sm = src[k].sample ? src[k].sample : "?";
      uint32_t u = 0; while (u < nu && usm[u] != src[k].sample) u++;
      int32_t key = (u < nu && usplit[u]) ? src[k].hap : INT32_MIN;
      size_t cut = strlen(sm);
      if (colsep) { const char *e = strchr(sm, colsep); if (e) cut = (size_t)(e - sm); }
      size_t L = cut + 24;
      nmbuf = realloc(nmbuf, L);
      memcpy(nmbuf, sm, cut); nmbuf[cut] = 0;
      if (key != INT32_MIN) snprintf(nmbuf + cut, L - cut, "#%d", (int)src[k].hap);
      uint32_t c = 0;
      while (c < g->ncpool && strcmp(g->cpool[c], nmbuf)) c++;
      if (c == g->ncpool) {
        g->cpool[c] = strdup(nmbuf); g->ncpool++;
      }
      g->pname[k] = g->cpool[c];
    }
    free(nmbuf);
    free(usm); free(uhap); free(usplit);
  }

  // --colors: keep only the whitelisted colors. Applied AFTER interning, so names match a
  // full-panel build; the pool is then compacted, which drops nc and with it nw = ceil(nc/64).
  uint8_t *keepc = NULL;                               // per color of the full pool
  uint8_t *keepw = NULL;                               // per walk; NULL = keep all
  uint32_t nkeep = ns;
  uint32_t nc_all = g->ncpool;
  if (colors && g->ncpool) {
    const wl_t *w = colors;
    keepc = calloc(g->ncpool, 1);
    for (uint32_t c = 0; c < g->ncpool; c++) keepc[c] = (uint8_t)wl_hit(w, g->cpool[c]);
    uint32_t nc2 = 0;
    for (uint32_t c = 0; c < g->ncpool; c++) if (keepc[c]) nc2++;
    if (!nc2) {
      fprintf(stderr, "--colors: none of the %u entries matched a color of this "
                      "graph (first color here is \"%s\")\n", w->n, g->cpool[0]);
      free(keepc); free(g->pname);
      for (uint32_t c = 0; c < g->ncpool; c++) free(g->cpool[c]);
      free(g->cpool); free(g); return NULL;
    }
    // Report entries that matched nothing: almost always a naming mistake.
    uint32_t miss = 0;
    for (uint32_t i = 0; i < w->n; i++) {
      int hit = 0;
      for (uint32_t c = 0; c < g->ncpool && !hit; c++)
        if (keepc[c] && name_hit(w->e[i], g->cpool[c])) hit = 1;
      if (!hit) { if (miss < 5 && !bg_quiet) fprintf(stderr, "--colors: no match for \"%s\"\n", w->e[i]); miss++; }
    }
    if (miss > 5 && !bg_quiet) fprintf(stderr, "--colors: ...and %u more unmatched entries\n", miss - 5);
    // Compact the pool, freeing the dropped names, and repoint every walk.
    char **np = malloc(nc2 * sizeof(char *));
    uint32_t *map = malloc(g->ncpool * sizeof(uint32_t));
    uint32_t j = 0;
    for (uint32_t c = 0; c < g->ncpool; c++) {
      if (keepc[c]) { map[c] = j; np[j++] = g->cpool[c]; }
      else { map[c] = UINT32_MAX; free(g->cpool[c]); }
    }
    // pname points into cpool, so translate through the pointer, not an index.
    keepw = malloc(ns ? ns : 1);
    nkeep = 0;
    for (uint32_t k = 0; k < ns; k++) {
      uint32_t c = 0; while (c < g->ncpool && g->cpool[c] != g->pname[k]) c++;
      keepw[k] = (c < g->ncpool && map[c] != UINT32_MAX);
      if (keepw[k]) nkeep++;
    }
    free(g->cpool); g->cpool = np; g->ncpool = nc2;
    free(map);
    if (!bg_quiet)
      fprintf(stderr, "--colors: %u of %u colors kept, %u of %u walks\n",
              nc2, nc_all, nkeep, ns);
  }
  free(keepc);

  // walks: oriented vertex ids verbatim. Built first: G_P adjacency needs them. Under --colors
  // only the kept walks are materialised, so every later pass sees exactly the sub-panel.
  g->np = nkeep;
  g->pv = nkeep ? malloc(nkeep * sizeof(uint32_t *)) : NULL;
  g->pn = nkeep ? malloc(nkeep * sizeof(uint32_t))   : NULL;
  char **kept_name = nkeep ? malloc(nkeep * sizeof(char *)) : NULL;
  g->psrc = nkeep ? malloc(nkeep * sizeof(uint32_t)) : NULL;
  uint32_t kk = 0;
  for (uint32_t k = 0; k < ns; k++) {
    if (keepw && !keepw[k]) continue;
    gfa_walk_t *p = &src[k];
    g->psrc[kk] = k;
    kept_name[kk] = g->pname[k];
    g->pn[kk] = p->n_v;
    if (p->n_v > 0) {
      g->pv[kk] = malloc(p->n_v * sizeof(uint32_t));
      for (int32_t s = 0; s < p->n_v; s++) g->pv[kk][s] = p->v[s];
    } else g->pv[kk] = NULL;
    kk++;
  }
  if (keepw) { free(g->pname); g->pname = kept_name; free(keepw); }
  else free(kept_name);

  if (full_graph) {
    g->adj  = malloc(g->n * sizeof(uint32_t *));
    g->adjn = malloc(g->n * sizeof(uint32_t));
    for (uint32_t v = 0; v < g->n; v++) {
      gfa_arc_t *a = gfa_arc_a(gfa, v);
      uint32_t na = gfa_arc_n(gfa, v);
      g->adjn[v] = na;
      if (na) {
        g->adj[v] = malloc(na * sizeof(uint32_t));
        for (uint32_t k = 0; k < na; k++) g->adj[v][k] = a[k].w;
      } else g->adj[v] = NULL;
    }
  } else if (!bg_nopadj) {
    build_padj(g);
  }
  return g;
}

void bg_free(bg_t *g) {
  if (!g) return;
  if (g->adj) { for (uint32_t v = 0; v < g->n; v++) free(g->adj[v]); free(g->adj); }
  free(g->adjn);
  if (g->pv) { for (uint32_t k = 0; k < g->np; k++) free(g->pv[k]); free(g->pv); }
  free(g->pn); free(g->psrc);
  free(g->sname); free(g->seq); free(g->slen); free(g->pname);  // borrowed contents
  if (g->cpool) { for (uint32_t c = 0; c < g->ncpool; c++) free(g->cpool[c]); }
  free(g->cpool);                                      // color names are OWNED
  free(g);
}

// ---- the graph of a splice: the walks the painting needs, and no others ------------
//
// A splice reads few of a panel's colours, but the generic reader parses every W-line (most of
// the bytes of the GFA) with one name -> id hash lookup per step. This reader gives EVERY W-line
// its record (sample, haplotype, contig, start, end), so bg_build names the colours as for the
// whole graph, and parses the steps only of the samples in `want`; the other W-lines keep an
// empty walk. S and L lines go through gfatools' own parsers. An unwanted W-line also loses its
// optional tags, which nothing in gpeb reads.
int gfa_parse_S(gfa_t *g, char *s);   // gfatools/gfa-io.c: weak definitions, no prototype in its headers
int gfa_parse_L(gfa_t *g, char *s);
int gfa_parse_W(gfa_t *g, char *s);

typedef struct { char *s; size_t l, m; } lbuf_t;
static void lbuf_add(lbuf_t *b, const char *p, size_t n) {
  if (b->l + n + 1 > b->m) { b->m = (b->l + n + 1) * 3 / 2 + 4096; b->s = realloc(b->s, b->m); }
  memcpy(b->s + b->l, p, n); b->l += n; b->s[b->l] = 0;
}

// a W-line is wanted when some name is its sample or a colour of it (SAMPLE#HAP)
static int sample_wanted(char *const *want, uint32_t nw, const char *sm, size_t L) {
  if (!want) return 1;                                   // NULL: every sample
  for (uint32_t i = 0; i < nw; i++) {
    size_t le = strlen(want[i]);
    if (le == L && !memcmp(want[i], sm, L)) return 1;
    if (le > L && !memcmp(want[i], sm, L) && want[i][L] == '#') return 1;
  }
  return 0;
}

int gfa_aux_parse(char *s, uint8_t **data, int *max);

// Segment names that are plain integers in strictly increasing order (every vg and minigraph-cactus
// GFA) are found by a search in a sorted array that starts from the previous step's position: a
// walk goes through nearly consecutive ids, so most steps cost one comparison. Anything else goes
// to gfa_name2id.
struct numtab_s { uint32_t n_seg; uint64_t *nm; int ok; uint32_t hint; };

static void numtab_build(numtab_t *t, const gfa_t *g) {
  free(t->nm); t->nm = NULL; t->ok = 0; t->hint = 0; t->n_seg = g->n_seg;
  if (!g->n_seg) return;
  t->nm = malloc((size_t)g->n_seg * sizeof(uint64_t));
  for (uint32_t z = 0; z < g->n_seg; z++) {
    const char *s = g->seg[z].name; uint64_t x = 0; size_t L = 0;
    if (!s || !*s) goto fail;
    for (; s[L]; L++) { if (s[L] < '0' || s[L] > '9' || L >= 18) goto fail; x = x * 10 + (uint64_t)(s[L] - '0'); }
    if (s[0] == '0' && L > 1) goto fail;                 // "007" and "7" are different names
    if (z && x <= t->nm[z - 1]) goto fail;               // a binary search needs them strictly increasing
    t->nm[z] = x;
  }
  t->ok = 1; return;
fail:
  free(t->nm); t->nm = NULL;
}

static int32_t numtab_find(numtab_t *t, uint64_t x) {
  const uint64_t *a = t->nm; uint32_t n = t->n_seg, h = t->hint, lo, hi, step = 1;
  if (h >= n) h = n - 1;
  if (a[h] == x) return (int32_t)h;
  if (a[h] < x) {                                        // gallop up from the last hit: a[lo - 1] < x
    lo = h + 1; hi = lo + step < n ? lo + step : n;
    while (hi < n && a[hi - 1] < x) { lo = hi; step <<= 1; hi = lo + step < n ? lo + step : n; }
  } else {                                               // gallop down: a[hi] > x
    hi = h; lo = hi > step ? hi - step : 0;
    while (lo > 0 && a[lo] > x) { hi = lo; step <<= 1; lo = hi > step ? hi - step : 0; }
  }
  while (lo < hi) { uint32_t mid = lo + (hi - lo) / 2; if (a[mid] < x) lo = mid + 1; else hi = mid; }
  if (lo < n && a[lo] == x) { t->hint = lo; return (int32_t)lo; }
  return -1;
}

// gfatools' gfa_parse_W for a wanted W-line, with the steps resolved by numtab_find.
static int gfa_parse_W_fast(gfa_t *g, char *s, numtab_t *nt) {
  char *p = s + 2, *fld[5];
  for (int i = 0; i < 5; i++) {
    fld[i] = p; char *e = strchr(p, '\t');
    if (!e) return -1;
    *e = 0; p = e + 1;
  }
  char *steps = p, *se = strchr(steps, '\t');
  size_t sl = se ? (size_t)(se - steps) : strlen(steps);
  char *rest = se ? se + 1 : NULL, *end = steps + sl;
  gfa_walk_t t; GFA_BZERO(&t, 1);
  t.hap = atoi(fld[1]); t.st = atol(fld[3]); t.en = atol(fld[4]);
  size_t cap = sl / 8 + 16, n = 0;
  uint32_t *v = malloc(cap * sizeof(uint32_t));
  for (char *q = steps; q < end; ) {
    char o = *q, *name = q + 1, *e = name;
    while (e < end && *e != '>' && *e != '<') e++;
    size_t L = (size_t)(e - name); int32_t seg = -1, done = 0;
    if (nt->ok && L > 0 && L < 19 && !(name[0] == '0' && L > 1)) {
      uint64_t x = 0; size_t k = 0;
      for (; k < L; k++) { unsigned d = (unsigned char)name[k] - '0'; if (d > 9) break; x = x * 10 + d; }
      if (k == L) { seg = numtab_find(nt, x); done = 1; }
    }
    if (!done) { char sv = *e; *e = 0; seg = gfa_name2id(g, name); *e = sv; }
    if (seg >= 0) {
      if (n == cap) { cap += cap >> 1; v = realloc(v, cap * sizeof(uint32_t)); }
      v[n++] = (uint32_t)seg << 1 | (o == '<');
    }
    q = e;
  }
  t.n_v = (int32_t)n; t.v = v;
  int l_aux, m_aux = 0; uint8_t *aux = 0;
  l_aux = gfa_aux_parse(rest, &aux, &m_aux);
  t.sample = gfa_sample_add(g, fld[0]);
  t.snid = gfa_sseq_add(g, fld[2]);
  if (l_aux > 0) t.aux.m_aux = m_aux, t.aux.l_aux = l_aux, t.aux.aux = aux;
  else if (aux) free(aux);
  GFA_GROW(gfa_walk_t, g->walk, g->n_walk, g->m_walk);
  g->walk[g->n_walk++] = t;
  return 0;
}

// One line of a plain GFA, as gfa_read_walks and the .gfai reader take it: S and L through
// gfatools' parsers, a W-line through gfa_parse_W_fast, or as a header with no steps (`wskip`).
// 1: a FASTA-style line, so not a plain GFA.
static int gfa_take_line(gfa_t *g, lbuf_t *ln, int wskip, numtab_t *nt) {
  if (ln->l >= 3 && ln->s[1] == '\t') {
    char t = ln->s[0];
    if (t == 'S') gfa_parse_S(g, ln->s);
    else if (t == 'L') gfa_parse_L(g, ln->s);
    else if (t == 'W') {
      if (wskip) gfa_parse_W(g, ln->s);
      else {
        if (nt->n_seg != g->n_seg) numtab_build(nt, g);  // first walk, or a segment came after one: (re)build
        gfa_parse_W_fast(g, ln->s, nt);
      }
    }
  } else if (ln->l && ln->s[0] == '>') return 1;
  return 0;
}

// NULL when the file is not a plain GFA (a FASTA-style segment line): the caller falls back to gfa_read.
static gfa_t *gfa_read_walks(const char *path, char *const *want, uint32_t nw, uint64_t *nskip) {
  gzFile fp = gzopen(path, "r");
  if (!fp) return NULL;
  gzbuffer(fp, 1 << 20);                                 // the default 8 KB buffer means very many small read() calls
  enum { CH = 1 << 22 };
  char *buf = malloc(CH); size_t n = 0, i = 0;
  gfa_t *g = gfa_init(); lbuf_t ln = {0, 0, 0};
  numtab_t nt = {0, NULL, 0, 0};
  int eof = 0, bad = 0; *nskip = 0;
  for (;;) {
    ln.l = 0; int skip = 0, eol = 0, wskip = 0;
    if (ln.s) ln.s[0] = 0;
    while (!eol) {
      if (i >= n) {
        if (eof) break;
        int k = gzread(fp, buf, CH);
        if (k <= 0) { eof = 1; break; }
        n = (size_t)k; i = 0;
      }
      char *p = buf + i; size_t av = n - i;
      if (skip) {
        char *e = memchr(p, '\n', av);
        if (e) { i += (size_t)(e - p) + 1; eol = 1; } else i = n;
        continue;
      }
      if (!ln.l && av >= 2 && p[0] == 'W' && p[1] == '\t') {
        size_t lim = av < 4096 ? av : 4096, h = 0, t2 = 0; int nt = 0;
        for (size_t j = 0; j < lim && p[j] != '\n'; j++)
          if (p[j] == '\t') { nt++; if (nt == 2) t2 = j; if (nt == 6) { h = j; break; } }
        if (h && !sample_wanted(want, nw, p + 2, t2 - 2)) {   // header complete and not wanted: keep the header only
          lbuf_add(&ln, p, h + 1); i += h + 1; skip = 1; wskip = 1; (*nskip)++;
          continue;
        }
      }
      char *e = memchr(p, '\n', av); size_t take = e ? (size_t)(e - p) : av;
      lbuf_add(&ln, p, take); i += take + (e ? 1 : 0);
      if (e) eol = 1;
    }
    if (gfa_take_line(g, &ln, wskip, &nt)) { bad = 1; break; }
    if (!eol) break;
  }
  free(ln.s); free(buf); free(nt.nm); gzclose(fp);
  if (bad) { gfa_destroy(g); return NULL; }
  gfa_finalize(g);
  return g;
}


static void gfai_free(gfai_t *X) {
  if (!X) return;
  for (uint32_t i = 0; i < X->nw; i++) free(X->w[i].hdr);
  free(X->w); free(X->ck); free(X->ref); free(X->buf);
  if (X->nt) { free(X->nt->nm); free(X->nt); }
  if (X->f) fclose(X->f);
  free(X);
}

static int file_is_gz(const char *path) {
  FILE *f = fopen(path, "rb");
  if (!f) return 0;
  int a = fgetc(f), b = fgetc(f);
  fclose(f);
  return a == 0x1f && b == 0x8b;
}

// The segment of a step named [name, e), found as gfa_parse_W_fast finds it; -1 if none.
// *e is written and restored: the caller's buffer has a byte to spare after the steps.
static int32_t step_id(const gfa_t *g, numtab_t *nt, char *name, char *e) {
  size_t L = (size_t)(e - name);
  if (nt->ok && L > 0 && L < 19 && !(name[0] == '0' && L > 1)) {
    uint64_t x = 0; size_t k = 0;
    for (; k < L; k++) { unsigned d = (unsigned char)name[k] - '0'; if (d > 9) break; x = x * 10 + d; }
    if (k == L) return numtab_find(nt, x);
  }
  char sv = *e; *e = 0;
  int32_t seg = gfa_name2id(g, name);
  *e = sv;
  return seg;
}

// The `cnt` steps of chunk K, in file order. -1 if the bytes are not exactly that many steps.
static int gfai_steps(gfai_t *X, const gfai_ck_t *K, uint32_t *out, uint32_t cnt) {
  if ((size_t)K->len + 1 > X->bcap) { X->bcap = (size_t)K->len + 1; X->buf = realloc(X->buf, X->bcap); }
  if (fseeko(X->f, (off_t)K->off, SEEK_SET) || fread(X->buf, 1, K->len, X->f) != K->len) return -1;
  X->nload++; X->bload += K->len;
  char *q = X->buf, *end = q + K->len;
  *end = 0;
  uint32_t n = 0;
  while (q < end) {
    char o = *q;
    if (o != '>' && o != '<') return -1;
    char *name = q + 1, *e = name;
    while (e < end && *e != '>' && *e != '<') e++;
    int32_t seg = step_id(X->gfa, X->nt, name, e);
    if (seg < 0 || n == cnt) return -1;
    out[n++] = (uint32_t)seg << 1 | (o == '<');
    q = e;
  }
  return n == cnt ? 0 : -1;
}

static int fw_u64(FILE *f, uint64_t x) { return fwrite(&x, 8, 1, f) == 1 ? 0 : -1; }
static int fw_u32(FILE *f, uint32_t x) { return fwrite(&x, 4, 1, f) == 1 ? 0 : -1; }
static int fr_u64(FILE *f, uint64_t *x) { return fread(x, 8, 1, f) == 1 ? 0 : -1; }
static int fr_u32(FILE *f, uint32_t *x) { return fread(x, 4, 1, f) == 1 ? 0 : -1; }

// Write <gfa_path>.gfai: the GFA `gfa` was parsed from, with the keys of the walks of `refname`
// in `g`. Every step of every W-line is read again from the file and checked against `gfa`, so
// the chunks parse at splice time into exactly the steps the whole-file reader would give.
// 0: written. 1: not written, and said why -- the caller goes on without it.
static int gfai_write(const char *gfa_path, const gfa_t *gfa, const bg_t *g, const char *refname) {
  double t0 = now_s();
  if (!gfa->n_walk) { fprintf(stderr, "gfai: %s has no W-lines: no .gfai\n", gfa_path); return 1; }
  if (file_is_gz(gfa_path)) {
    fprintf(stderr, "gfai: %s is compressed, and a compressed file cannot be read at an offset: no .gfai"
                    " (a splice will read it whole)\n", gfa_path);
    return 1;
  }
  struct stat sb;
  if (stat(gfa_path, &sb)) { fprintf(stderr, "gfai: cannot stat %s: no .gfai\n", gfa_path); return 1; }
  sref_t *R = sref_build(gfa, g, refname);
  if (!R) return 1;
  uint32_t ch = gfai_chunk, nw = gfa->n_walk;
  gfai_w_t *W = calloc(nw, sizeof(gfai_w_t));
  uint64_t nck = 0, nmono = 0, nflip = 0;
  for (uint32_t i = 0; i < nw; i++) { W[i].n = (uint32_t)gfa->walk[i].n_v; W[i].nck = (W[i].n + ch - 1) / ch; W[i].ck0 = nck; nck += W[i].nck; }
  gfai_ck_t *CK = calloc(nck ? nck : 1, sizeof(gfai_ck_t));
  // what swl_build computes, chunk by chunk
  for (uint32_t i = 0; i < nw; i++) {
    const gfa_walk_t *wk = &gfa->walk[i]; gfai_w_t *Wi = &W[i]; uint32_t n = Wi->n;
    uint64_t nanc = 0, nsame = 0;
    for (uint32_t s = 0; s < n; s++) {
      uint32_t v = wk->v[s]; int32_t k = R->kof[v >> 1];
      if (k >= 0) { nanc++; if ((v & 1) == R->rs[k]) nsame++; }
      Wi->bp += g->slen[v >> 1];
    }
    Wi->flip = nanc && 2 * nsame < nanc; Wi->has = nanc > 0; Wi->mono = 1;
    nflip += Wi->flip;
    for (uint32_t j = 0; j < n && Wi->has; j++) {
      uint32_t v = Wi->flip ? wk->v[n - 1 - j] ^ 1u : wk->v[j]; int32_t k = R->kof[v >> 1];
      if (k >= 0) { Wi->first = j; Wi->rp_first = R->rp[k]; break; }
    }
    int64_t base = 0, prev = INT64_MIN;
    for (uint32_t j = 0; j < n; j++) {
      gfai_ck_t *K = &CK[Wi->ck0 + j / ch];
      if (j % ch == 0) { K->base = base; K->kmin = INT64_MAX; K->kmax = INT64_MIN; }
      uint32_t v = Wi->flip ? wk->v[n - 1 - j] ^ 1u : wk->v[j]; int32_t k = R->kof[v >> 1];
      if (k >= 0) base = R->rp[k];
      if (!Wi->has) continue;
      int64_t key = j < Wi->first ? 2 * Wi->rp_first - 1 : 2 * base + (k >= 0 ? 0 : 1);
      if (key < K->kmin) K->kmin = key;
      if (key > K->kmax) K->kmax = key;
      if (key < prev) Wi->mono = 0;
      prev = key;
    }
    nmono += Wi->has && Wi->mono;
  }
  // the byte offsets, from the file itself: every step is checked against the parsed walk
  numtab_t nt = {0, NULL, 0, 0}; numtab_build(&nt, gfa);
  FILE *f = fopen(gfa_path, "rb");
  enum { BS = 1 << 23 };
  char *blk = malloc(BS), hdr[4096], nm[256];
  uint64_t boff = 0, loff = 0, *bo = NULL;
  uint32_t wi = 0, ntab = 0, hl = 0, nl = 0, x = 0;
  int st = 0, instep = 0, bad = !f; char so = 0;
  const char *why = !f ? "cannot read it" : NULL;
  gfai_w_t *Wc = NULL;
  #define GFAI_BAD(m) do { bad = 1; why = (m); } while (0)
  // finish the step being read: it is step x - 1 of W-line wi
  #define GFAI_STEP() do { if (instep) { int32_t sg = step_id(gfa, &nt, nm, nm + nl);                     \
      if (sg < 0 || x - 1 >= Wc->n || ((uint32_t)sg << 1 | (so == '<')) != gfa->walk[wi].v[x - 1])         \
        { GFAI_BAD("a step reads differently from the parsed walk"); }                                     \
      instep = 0; } } while (0)
  // record the offset of file step x when it bounds a chunk; GFAI_END: the end of the steps
  #define GFAI_BOUND(o) do { uint32_t n_ = Wc->n;                                                        \
      if (!Wc->flip) { if (x % ch == 0) bo[x / ch] = (o); }                                             \
      else { if ((n_ - x) % ch == 0) bo[(n_ - x) / ch] = (o); if (x == 0) bo[Wc->nck] = (o); } } while (0)
  #define GFAI_END(o) do { if (Wc->flip) bo[0] = (o); else bo[Wc->nck] = (o); } while (0)
  #define GFAI_LINE(end) do { if (x != Wc->n) GFAI_BAD("a W-line has a different number of steps");    \
      else { Wc->llen = (end) - loff;                                                                   \
        for (uint32_t c_ = 0; c_ < Wc->nck; c_++) {                                                     \
          uint64_t s_ = Wc->flip ? bo[c_ + 1] : bo[c_], e_ = Wc->flip ? bo[c_] : bo[c_ + 1];            \
          if (s_ == UINT64_MAX || e_ == UINT64_MAX || e_ <= s_ || e_ - s_ > UINT32_MAX)                 \
            { GFAI_BAD("a chunk has no bounds"); break; }                                               \
          CK[Wc->ck0 + c_].off = s_; CK[Wc->ck0 + c_].len = (uint32_t)(e_ - s_); } }                    \
      free(bo); bo = NULL; wi++; st = 0; } while (0)
  size_t got;
  while (!bad && (got = fread(blk, 1, BS, f)) > 0) {
    for (size_t p = 0; p < got && !bad; ) {
      char c = blk[p];
      switch (st) {
      case 0:                                              // the start of a line
        if (c == 'W') { st = 5; loff = boff + p; p++; }
        else if (c == '\n') p++;
        else st = 1;
        break;
      case 1: {                                            // any other line: to its end
        char *e = memchr(blk + p, '\n', got - p);
        if (e) { p = (size_t)(e - blk) + 1; st = 0; } else p = got;
        break; }
      case 5:                                              // "W" at the start of a line
        if (c != '\t') { st = 1; break; }
        if (wi >= nw) { GFAI_BAD("more W-lines than the parsed graph"); break; }
        Wc = &W[wi]; Wc->loff = loff; hl = 0; hdr[hl++] = 'W'; hdr[hl++] = '\t'; ntab = 1; st = 2; p++;
        break;
      case 2:                                              // the header, to the 6th tab
        if (c == '\n' || hl >= sizeof hdr) { GFAI_BAD("a W-line header is malformed or too long"); break; }
        hdr[hl++] = c; p++;
        if (c == '\t' && ++ntab == 6) {
          size_t sl = strcspn(hdr + 2, "\t");
          const char *sm = gfa->walk[wi].sample ? gfa->walk[wi].sample : "";
          if (strlen(sm) != sl || memcmp(sm, hdr + 2, sl)) { GFAI_BAD("a W-line's sample differs from the parsed graph"); break; }
          Wc->hlen = hl; Wc->hdr = malloc(hl + 1); memcpy(Wc->hdr, hdr, hl); Wc->hdr[hl] = 0;
          bo = malloc((Wc->nck + 1) * sizeof(uint64_t));
          for (uint32_t z = 0; z <= Wc->nck; z++) bo[z] = UINT64_MAX;
          x = 0; instep = 0; st = 3;
        }
        break;
      case 3:                                              // the steps
        if (c == '>' || c == '<') {
          GFAI_STEP(); if (bad) break;
          GFAI_BOUND(boff + p); so = c; nl = 0; instep = 1; x++; p++;
        } else if (c == '\t' || c == '\n') {
          GFAI_STEP(); if (bad) break;
          GFAI_END(boff + p);
          if (c == '\n') { GFAI_LINE(boff + p); p++; } else { st = 4; p++; }
        } else {
          if (!instep || nl + 1 >= sizeof nm) { GFAI_BAD("a step name is malformed or too long"); break; }
          nm[nl++] = c; p++;
        }
        break;
      case 4: {                                            // the tags after the steps: to the end of the line
        char *e = memchr(blk + p, '\n', got - p);
        if (e) { GFAI_LINE(boff + (uint64_t)(e - blk)); p = (size_t)(e - blk) + 1; } else p = got;
        break; }
      }
    }
    boff += got;
  }
  if (!bad && st == 3) { GFAI_STEP(); if (!bad) { GFAI_END(boff); GFAI_LINE(boff); } }   // the file ends in a W-line
  else if (!bad && st == 4) GFAI_LINE(boff);
  else if (!bad && st == 2) GFAI_BAD("the file ends in a W-line header");
  if (!bad && wi != nw) GFAI_BAD("fewer W-lines than the parsed graph");
  #undef GFAI_BAD
  #undef GFAI_STEP
  #undef GFAI_BOUND
  #undef GFAI_END
  #undef GFAI_LINE
  free(bo); free(blk); free(nt.nm);
  if (f) fclose(f);
  uint64_t fsz = 0, fh = bad ? 0 : file_hash(gfa_path, &fsz);
  if (!bad && fsz != (uint64_t)sb.st_size) { bad = 1; why = "the GFA changed while it was read"; }
  size_t pl = strlen(gfa_path);
  char *out = malloc(pl + 6), *tmp = malloc(pl + 10);
  snprintf(out, pl + 6, "%s.gfai", gfa_path); snprintf(tmp, pl + 10, "%s.gfai.tmp", gfa_path);
  FILE *fo = NULL;
  if (!bad && !(fo = fopen(tmp, "wb"))) { bad = 1; why = "cannot write it next to the GFA"; }
  if (!bad) {
    int e = 0;
    uint32_t rl = (uint32_t)strlen(refname);
    e |= fwrite(GFAI_MAGIC, 1, 4, fo) != 4;
    e |= fw_u32(fo, GFAI_VERSION); e |= fw_u64(fo, fsz); e |= fw_u64(fo, fh);
    e |= fw_u64(fo, (uint64_t)ST_MTIM(sb).tv_sec); e |= fw_u64(fo, (uint64_t)ST_MTIM(sb).tv_nsec);
    e |= fw_u32(fo, ch); e |= fw_u32(fo, gfa->n_seg); e |= fw_u32(fo, R->nref);
    e |= fw_u32(fo, rl); e |= fwrite(refname, 1, rl, fo) != rl;
    e |= fw_u32(fo, nw); e |= fw_u64(fo, nck);
    for (uint32_t i = 0; i < nw && !e; i++) {
      const gfai_w_t *Wi = &W[i]; uint8_t fl[4] = { Wi->flip, Wi->mono, Wi->has, 0 };
      e |= fw_u64(fo, Wi->loff); e |= fw_u64(fo, Wi->llen); e |= fw_u32(fo, Wi->n); e |= fw_u32(fo, Wi->first);
      e |= fwrite(fl, 1, 4, fo) != 4; e |= fw_u64(fo, (uint64_t)Wi->rp_first); e |= fw_u64(fo, Wi->bp);
      e |= fw_u64(fo, Wi->ck0); e |= fw_u32(fo, Wi->nck); e |= fw_u32(fo, Wi->hlen);
      e |= fwrite(Wi->hdr, 1, Wi->hlen, fo) != Wi->hlen;
    }
    for (uint64_t c = 0; c < nck && !e; c++) {
      const gfai_ck_t *K = &CK[c];
      e |= fw_u64(fo, K->off); e |= fw_u32(fo, K->len);
      e |= fw_u64(fo, (uint64_t)K->base); e |= fw_u64(fo, (uint64_t)K->kmin); e |= fw_u64(fo, (uint64_t)K->kmax);
    }
    e |= fclose(fo) != 0;
    if (e || rename(tmp, out)) { remove(tmp); bad = 1; why = "cannot write it next to the GFA"; }
  }
  if (bad) fprintf(stderr, "gfai: no %s: %s\n", out, why ? why : "?");
  else {
    struct stat ob; uint64_t osz = stat(out, &ob) ? 0 : (uint64_t)ob.st_size;
    fprintf(stderr, "gfai: wrote %s: %u W-lines (%llu against %s, %llu with keys that never decrease), %llu chunks of %u"
                    " steps, %.1f MB (%.1fs)\n", out, nw, (unsigned long long)nflip, refname,
            (unsigned long long)nmono, (unsigned long long)nck, ch, osz / 1e6, now_s() - t0);
  }
  for (uint32_t i = 0; i < nw; i++) free(W[i].hdr);
  free(W); free(CK); free(out); free(tmp); sref_free(R);
  return bad;
}

// The .gfai of `gfa_path`, if it can be used for a splice on `refname`; NULL otherwise, with
// a line saying why (the splice then reads the whole GFA, as without one). `hash`: the GFA's
// FNV-1a when the caller has computed it anyway (have_hash), checked too.
static gfai_t *gfai_open(const char *gfa_path, const char *refname, uint64_t hash, int have_hash) {
  if (gfai_off) { fprintf(stderr, "splice: --no-gfai: reading the whole GFA\n"); return NULL; }
  if (file_is_gz(gfa_path)) { fprintf(stderr, "splice: %s is compressed: reading the whole GFA\n", gfa_path); return NULL; }
  size_t pl = strlen(gfa_path);
  char *path = malloc(pl + 6); snprintf(path, pl + 6, "%s.gfai", gfa_path);
  FILE *f = fopen(path, "rb");
  if (!f) { fprintf(stderr, "splice: no %s: reading the whole GFA\n", path); free(path); return NULL; }
  gfai_t *X = calloc(1, sizeof *X);
  const char *why = NULL; char mg[4]; uint32_t ver = 0, rl = 0; uint64_t ms = 0, mn = 0;
  struct stat sb;
  if (fread(mg, 1, 4, f) != 4 || memcmp(mg, GFAI_MAGIC, 4) || fr_u32(f, &ver)) why = "not a .gfai";
  else if (ver != GFAI_VERSION) why = "another version of the format";
  else if (fr_u64(f, &X->size) || fr_u64(f, &X->hash) || fr_u64(f, &ms) || fr_u64(f, &mn) ||
           fr_u32(f, &X->ch) || fr_u32(f, &X->n_seg) || fr_u32(f, &X->nref) || fr_u32(f, &rl) || rl > 4096)
    why = "truncated";
  if (!why) {
    X->mt_s = (int64_t)ms; X->mt_ns = (int64_t)mn;
    X->ref = malloc(rl + 1);
    if (fread(X->ref, 1, rl, f) != rl) why = "truncated";
    X->ref[rl] = 0;
  }
  if (!why && stat(gfa_path, &sb)) why = "cannot stat the GFA";
  if (!why && ((uint64_t)sb.st_size != X->size || (int64_t)ST_MTIM(sb).tv_sec != X->mt_s ||
               (int64_t)ST_MTIM(sb).tv_nsec != X->mt_ns))
    why = "written for another version of the GFA (its size or modification time changed)";
  if (!why && have_hash && hash != X->hash) why = "written for another version of the GFA (its contents changed)";
  if (!why && strcmp(X->ref, refname)) why = "its keys are of another reference";
  if (!why && (fr_u32(f, &X->nw) || fr_u64(f, &X->nck) || !X->ch)) why = "truncated";
  if (!why) {
    X->w = calloc(X->nw ? X->nw : 1, sizeof(gfai_w_t));
    for (uint32_t i = 0; i < X->nw && !why; i++) {
      gfai_w_t *W = &X->w[i]; uint8_t fl[4]; uint64_t rp = 0;
      if (fr_u64(f, &W->loff) || fr_u64(f, &W->llen) || fr_u32(f, &W->n) || fr_u32(f, &W->first) ||
          fread(fl, 1, 4, f) != 4 || fr_u64(f, &rp) || fr_u64(f, &W->bp) || fr_u64(f, &W->ck0) ||
          fr_u32(f, &W->nck) || fr_u32(f, &W->hlen) || W->hlen > 4096 || W->ck0 + W->nck > X->nck)
        { why = "truncated"; break; }
      W->flip = fl[0]; W->mono = fl[1]; W->has = fl[2]; W->rp_first = (int64_t)rp;
      W->hdr = malloc(W->hlen + 1);
      if (fread(W->hdr, 1, W->hlen, f) != W->hlen) why = "truncated";
      W->hdr[W->hlen] = 0;
    }
  }
  if (!why) {
    X->ck = malloc((X->nck ? X->nck : 1) * sizeof(gfai_ck_t));
    for (uint64_t c = 0; c < X->nck && !why; c++) {
      gfai_ck_t *K = &X->ck[c]; uint64_t b = 0, lo = 0, hi = 0;
      if (fr_u64(f, &K->off) || fr_u32(f, &K->len) || fr_u64(f, &b) || fr_u64(f, &lo) || fr_u64(f, &hi)) why = "truncated";
      K->base = (int64_t)b; K->kmin = (int64_t)lo; K->kmax = (int64_t)hi;
    }
  }
  fclose(f);
  if (!why && !(X->f = fopen(gfa_path, "rb"))) why = "cannot read the GFA";
  if (why) { fprintf(stderr, "splice: %s: %s; reading the whole GFA\n", path, why); free(path); gfai_free(X); return NULL; }
  fprintf(stderr, "splice: using %s (%u W-lines in %llu chunks of %u steps)\n", path, X->nw,
          (unsigned long long)X->nck, X->ch);
  free(path);
  return X;
}

// Bytes [a, b) of the GFA: whole lines, each taken as gfa_read_walks takes it. -1 on a read error.
static int gfai_lines(FILE *f, uint64_t a, uint64_t b, gfa_t *g, numtab_t *nt, lbuf_t *ln,
                      char *blk, size_t bsz, int *bad) {
  if (a >= b) return 0;
  if (fseeko(f, (off_t)a, SEEK_SET)) return -1;
  ln->l = 0; if (ln->s) ln->s[0] = 0;
  uint64_t left = b - a;
  while (left && !*bad) {
    size_t want = left < bsz ? (size_t)left : bsz, got = fread(blk, 1, want, f);
    if (!got) return -1;
    left -= got;
    for (size_t i = 0; i < got && !*bad; ) {
      char *e = memchr(blk + i, '\n', got - i);
      size_t take = e ? (size_t)(e - (blk + i)) : got - i;
      lbuf_add(ln, blk + i, take); i += take;
      if (e) { i++; if (gfa_take_line(g, ln, 0, nt)) *bad = 1; ln->l = 0; ln->s[0] = 0; }
    }
  }
  if (ln->l && !*bad) { if (gfa_take_line(g, ln, 0, nt)) *bad = 1; ln->l = 0; }
  return 0;
}

// The GFA through its .gfai: every line but the W-lines as in the file; the W-lines of the
// samples in `full` whole, every other one as its header (no steps), so the walks the splice
// reads lazily stay empty in the graph. *nskip: the W-lines of no sample in `all` (what
// gfa_read_walks(all) would skip); *nfull: the W-lines parsed whole.
static gfa_t *gfa_read_gfai(gfai_t *X, char *const *full, uint32_t nf, char *const *all, uint32_t na,
                            uint64_t *nskip, uint64_t *nfull) {
  gfa_t *g = gfa_init(); numtab_t nt = {0, NULL, 0, 0}; lbuf_t ln = {0, 0, 0};
  enum { BS = 1 << 22 };
  char *blk = malloc(BS);
  uint64_t pos = 0; int bad = 0, io = 0;
  *nskip = 0; *nfull = 0;
  for (uint32_t i = 0; i <= X->nw && !bad && !io; i++) {
    uint64_t upto = i < X->nw ? X->w[i].loff : X->size;
    if (gfai_lines(X->f, pos, upto, g, &nt, &ln, blk, BS, &bad)) { io = 1; break; }
    if (i == X->nw || bad) break;
    const gfai_w_t *W = &X->w[i];
    const char *sm = W->hdr + 2; size_t sl = strcspn(sm, "\t");
    if (!sample_wanted(all, na, sm, sl)) (*nskip)++;
    ln.l = 0;
    if (sample_wanted(full, nf, sm, sl)) {
      if (ln.m < W->llen + 1) { ln.m = W->llen + 1; ln.s = realloc(ln.s, ln.m); }
      if (fseeko(X->f, (off_t)W->loff, SEEK_SET) || fread(ln.s, 1, W->llen, X->f) != W->llen) { io = 1; break; }
      ln.l = W->llen; ln.s[ln.l] = 0;
      if (ln.l < 2 || ln.s[0] != 'W' || ln.s[1] != '\t') { io = 1; break; }
      (*nfull)++;
      if (gfa_take_line(g, &ln, 0, &nt)) bad = 1;
    } else {
      lbuf_add(&ln, W->hdr, W->hlen);
      if (gfa_take_line(g, &ln, 1, &nt)) bad = 1;
    }
    pos = W->loff + W->llen + 1;
  }
  free(ln.s); free(blk); free(nt.nm);
  if (bad || io || g->n_walk != X->nw) { gfa_destroy(g); return NULL; }
  gfa_finalize(g);
  X->nt = calloc(1, sizeof(numtab_t)); numtab_build(X->nt, g); X->gfa = g;
  return g;
}

// The names and the samples they belong to (SAMPLE for SAMPLE#HAP), once each.
static char **want_expand(char *const *names, uint32_t nn, uint32_t *nw_out) {
  char **want = malloc(2 * (nn ? nn : 1) * sizeof(char *)); uint32_t nw = 0;
  for (uint32_t i = 0; i < nn; i++) {
    const char *h = strchr(names[i], '#');
    for (int pass = 0; pass < 2; pass++) {
      if (pass && !h) break;
      size_t L = pass ? (size_t)(h - names[i]) : strlen(names[i]); uint32_t j = 0;
      while (j < nw && !(strlen(want[j]) == L && !memcmp(want[j], names[i], L))) j++;
      if (j == nw) want[nw++] = strndup(names[i], L);
    }
  }
  *nw_out = nw;
  return want;
}

// `who` names the option in the messages. names: the colours the splice will read (the
// painting's, the kept ones, the reference); the samples they belong to are added here.
// With a .gfai (*X), only `full` (the kept colours and the reference) is parsed whole and the
// other W-lines are read later, chunk by chunk; if the GFA cannot be read through it, *X is
// dropped (and freed) and the GFA is read as without one.
static int splice_graph_open(const char *who, const char *gfa_path, char *const *names, uint32_t nn,
                             char *const *full, uint32_t nfull, gfai_t **X, gfa_t **gfa, bg_t **g) {
  double t0 = now_s();
  uint32_t nw = 0, nwf = 0;
  char **want = want_expand(names, nn, &nw), **wantf = *X ? want_expand(full, nfull, &nwf) : NULL;
  uint64_t nskip = 0, nparsed = 0;
  gfa_t *gf = NULL;
  if (*X) {
    gf = gfa_read_gfai(*X, wantf, nwf, want, nw, &nskip, &nparsed);
    if (!gf) {
      fprintf(stderr, "splice: the GFA does not read as its .gfai says; reading the whole GFA\n");
      gfai_free(*X); *X = NULL;
    }
  }
  if (!gf) {
    gf = gfa_read_walks(gfa_path, want, nw, &nskip);
    if (!gf) gf = gfa_read(gfa_path);                    // not a plain GFA, or unreadable: the generic reader decides
    if (gf) nparsed = gf->n_walk - nskip;
  }
  int rc = 0; bg_t *bg = NULL;
  if (!gf || !gf->n_walk) { fprintf(stderr, "%s: %s has no W-lines\n", who, gfa_path); rc = 1; }
  else {
    double tr = now_s() - t0;
    wl_t wl = { want, nw };
    bg_quiet = 1; bg_nopadj = 1;
    bg = bg_build(gf, 1, 0, 0, nskip ? &wl : NULL);
    bg_quiet = 0; bg_nopadj = 0;
    if (!bg) rc = 1;
    else if (*X)
      fprintf(stderr, "splice: graph read in %.1fs (%llu of %u W-lines parsed whole, %llu more read from the .gfai as"
                      " needed, for %u colours), built in %.1fs\n", tr, (unsigned long long)nparsed, gf->n_walk,
              (unsigned long long)(gf->n_walk - nskip - nparsed), nn, now_s() - t0 - tr);
    else fprintf(stderr, "splice: graph read in %.1fs (%u of %u W-lines parsed for %u colours), built in %.1fs\n", tr,
                 (uint32_t)(gf->n_walk - nskip), gf->n_walk, nn, now_s() - t0 - tr);
  }
  for (uint32_t i = 0; i < nw; i++) free(want[i]);
  for (uint32_t i = 0; i < nwf; i++) free(wantf[i]);
  free(want); free(wantf);
  if (rc) { if (gf) gfa_destroy(gf); if (*X) (*X)->gfa = NULL; return rc; }
  *gfa = gf; *g = bg;
  return 0;
}

static void splice_names_add(char ***v, uint32_t *n, const char *s, size_t L) {
  for (uint32_t i = 0; i < *n; i++) if (strlen((*v)[i]) == L && !memcmp((*v)[i], s, L)) return;
  *v = realloc(*v, (*n + 1) * sizeof(char *)); (*v)[(*n)++] = strndup(s, L);
}
static void splice_names_csv(char ***v, uint32_t *n, const char *csv) {
  for (const char *p = csv; *p; ) {
    const char *e = strchr(p, ','); size_t L = e ? (size_t)(e - p) : strlen(p);
    if (L) splice_names_add(v, n, p, L);
    if (!e) break;
    p = e + 1;
  }
}

// the graph for a painting that is still in memory: parse it for the painting's colours only
static int splice_out(const mosaic_opt_t *opt, const gfa_t *gfa, const bg_t *g, const pidx_t *ix,
                      const uint16_t *path, uint32_t nB, uint32_t T) {
  if (gfa) return mosaic_splice(gfa, g, ix, path, nB, T, opt->gfa_out, opt->keep, NULL);
  if (!opt->gfa_in) { fprintf(stderr, "--mosaic-gfa needs the GFA the index was built from, as input\n"); return 1; }
  char **names = NULL, **full = NULL; uint32_t nn = 0, nf = 0;
  uint8_t *used = calloc(ix->nc ? ix->nc : 1, 1);
  for (size_t t = 0; t < (size_t)nB * T; t++) if (path[t] != NOCOL) used[path[t]] = 1;
  for (uint32_t c = 0; c < ix->nc; c++) if (used[c] && ix->cname[c]) splice_names_add(&names, &nn, ix->cname[c], strlen(ix->cname[c]));
  free(used);
  splice_names_csv(&names, &nn, opt->keep ? opt->keep : ix->blk_ref);
  splice_names_add(&names, &nn, ix->blk_ref, strlen(ix->blk_ref));
  splice_names_csv(&full, &nf, opt->keep ? opt->keep : ix->blk_ref);
  splice_names_add(&full, &nf, ix->blk_ref, strlen(ix->blk_ref));
  gfa_t *gf = NULL; bg_t *bg = NULL;
  gfai_t *X = gfai_open(opt->gfa_in, ix->blk_ref, opt->gfa_hash, opt->gfa_hashed);
  int rc = splice_graph_open("--mosaic-gfa", opt->gfa_in, names, nn, full, nf, &X, &gf, &bg);
  for (uint32_t i = 0; i < nn; i++) free(names[i]);
  for (uint32_t i = 0; i < nf; i++) free(full[i]);
  free(names); free(full);
  if (rc) { gfai_free(X); return rc; }
  rc = mosaic_splice(gf, bg, ix, path, nB, T, opt->gfa_out, opt->keep, X);
  bg_free(bg); gfa_destroy(gf); gfai_free(X);
  return rc;
}

// the colour names a path.tsv uses, as run_splice_file will read it (the fields of the track columns)
// (ref: the coordinate reference of its header, as run_splice_file reads it -- the last coord=)
static int painting_names(const char *painting, char ***names, uint32_t *nn, char *ref, size_t rcap) {
  FILE *f = fopen(painting, "r");
  if (!f) { fprintf(stderr, "--mosaic-splice: cannot read %s\n", painting); return 1; }
  char *line = NULL; size_t cap = 0; ssize_t L; uint32_t T = 0;
  while ((L = getline(&line, &cap, f)) > 0) {
    while (L && (line[L - 1] == '\n' || line[L - 1] == '\r')) line[--L] = 0;
    if (line[0] == '#') {
      for (char *t = strtok(line + 1, " \t"); t; t = strtok(NULL, " \t"))
        if (!strncmp(t, "coord=", 6)) { splice_names_add(names, nn, t + 6, strlen(t + 6)); snprintf(ref, rcap, "%s", t + 6); }
      continue;
    }
    if (!strncmp(line, "block\t", 6)) {
      T = 0;
      for (char *t = line; t; ) { if (!strncmp(t, "track", 5)) T++; char *e = strchr(t, '\t'); t = e ? e + 1 : NULL; }
      continue;
    }
    if (!L) continue;
    if (!T) T = 2;
    char *t = line; uint32_t fi = 0;
    while (t && fi < 3 + T) {
      char *e = strchr(t, '\t'); size_t FL = e ? (size_t)(e - t) : strlen(t);
      if (fi >= 3 && FL && !(FL == 1 && *t == '-')) splice_names_add(names, nn, t, FL);
      t = e ? e + 1 : NULL; fi++;
    }
  }
  free(line); fclose(f);
  return 0;
}

// --mosaic-splice: splice a painting that is already on disk (path.tsv, as written by
// --set-call mosaic) into a graph. No reads, no index: the reference colour and the
// block size come from the painting's header.
static int run_splice_file(const char *painting, const char *gfa_path, const char *out,
                           const char *keep) {
  char **pn = NULL, **full = NULL; uint32_t pnn = 0, nf = 0; char pref[256] = "";
  if (painting_names(painting, &pn, &pnn, pref, sizeof pref)) return 1;
  splice_names_csv(&pn, &pnn, keep ? keep : "");
  splice_names_csv(&full, &nf, keep ? keep : pref);
  if (pref[0]) splice_names_add(&full, &nf, pref, strlen(pref));
  gfa_t *gfa = NULL; bg_t *g = NULL;
  gfai_t *X = pref[0] ? gfai_open(gfa_path, pref, 0, 0) : NULL;
  int rco = splice_graph_open("--mosaic-splice", gfa_path, pn, pnn, full, nf, &X, &gfa, &g);
  for (uint32_t i = 0; i < pnn; i++) free(pn[i]);
  for (uint32_t i = 0; i < nf; i++) free(full[i]);
  free(pn); free(full);
  if (rco) { gfai_free(X); return 1; }
  FILE *f = fopen(painting, "r");
  if (!f) { fprintf(stderr, "--mosaic-splice: cannot read %s\n", painting); return 1; }
  char *line = NULL; size_t cap = 0; ssize_t L;
  char ref[256] = ""; uint32_t B = 0, T = 0;
  uint32_t nB = 0, rcap = 1024; uint16_t *path = NULL;
  uint64_t prev_end = 0, first = 1, blk_len = 0; int rc = 0;
  while ((L = getline(&line, &cap, f)) > 0) {
    while (L && (line[L - 1] == '\n' || line[L - 1] == '\r')) line[--L] = 0;
    if (line[0] == '#') {
      for (char *t = strtok(line + 1, " \t"); t; t = strtok(NULL, " \t")) {
        if (!strncmp(t, "coord=", 6)) snprintf(ref, sizeof ref, "%s", t + 6);
        else if (!strncmp(t, "B=", 2)) B = (uint32_t)strtoul(t + 2, NULL, 10);
      }
      continue;
    }
    if (!strncmp(line, "block\t", 6)) {                  // the header names the tracks: track1 .. trackT
      T = 0;
      for (char *t = line; t; ) {
        if (!strncmp(t, "track", 5)) T++;
        char *e = strchr(t, '\t'); t = e ? e + 1 : NULL;
      }
      continue;
    }
    if (!L) continue;
    if (!T) T = 2;                                         // no header: the two-track format
    if (!path) path = malloc((size_t)rcap * T * sizeof(uint16_t));
    char *fld[3 + 64 + 4]; uint32_t nf = 0;
    if (T > 64) { fprintf(stderr, "--mosaic-splice: more than 64 tracks\n"); rc = 1; break; }
    for (char *t = line; nf < 3 + T + 2; ) { fld[nf++] = t; char *e = strchr(t, '\t'); if (!e) break; *e = 0; t = e + 1; }
    if (nf < 3 + T) { fprintf(stderr, "--mosaic-splice: malformed row in %s\n", painting); rc = 1; break; }
    uint64_t s0 = strtoull(fld[1], NULL, 10), e0 = strtoull(fld[2], NULL, 10);
    if (!first && s0 != prev_end) {
      fprintf(stderr, "--mosaic-splice: blocks are not contiguous at block %s\n", fld[0]); rc = 1; break;
    }
    first = 0; prev_end = e0; blk_len = e0;
    if (nB == rcap) { rcap *= 2; path = realloc(path, (size_t)rcap * T * sizeof(uint16_t)); }
    for (uint32_t k = 0; k < T; k++) {
      if (!strcmp(fld[3 + k], "-")) { path[(size_t)nB * T + k] = NOCOL; continue; }
      uint32_t c = 0;
      while (c < g->ncpool && strcmp(g->cpool[c], fld[3 + k])) c++;
      if (c == g->ncpool) {
        // gpeb names a colour SAMPLE#HAP only where the sample has several haplotypes; a
        // painting made with another convention may say HG002#1 for a lone haplotype.
        const char *h = strchr(fld[3 + k], '#');
        if (h) { c = 0; while (c < g->ncpool && (strncmp(g->cpool[c], fld[3 + k], (size_t)(h - fld[3 + k])) ||
                                                  g->cpool[c][h - fld[3 + k]])) c++; }
      }
      if (c == g->ncpool) { fprintf(stderr, "--mosaic-splice: colour %s is not in the graph\n", fld[3 + k]); rc = 1; break; }
      path[(size_t)nB * T + k] = (uint16_t)c;
    }
    if (rc) break;
    nB++;
  }
  free(line); fclose(f);
  if (!rc && (!nB || !ref[0] || !B)) {
    fprintf(stderr, "--mosaic-splice: %s needs a header with coord= and B=, and at least one block\n", painting); rc = 1;
  }
  if (!rc) {
    pidx_t fake; memset(&fake, 0, sizeof fake);
    fake.nc = g->ncpool; fake.cname = (const char **)g->cpool;
    fake.blk_ref = ref; fake.blk_bp = B; fake.blk_len = blk_len;
    rc = mosaic_splice(gfa, g, &fake, path, nB, T, out, keep, X);
  }
  free(path); bg_free(g); gfa_destroy(gfa); gfai_free(X);
  return rc;
}


#ifdef GPEB_TEST
#include <assert.h>
#include <fcntl.h>               // AT_FDCWD: the .gfai test ages its GFA
#include <math.h>
// brute oracle: >=2 distinct simple paths src->x (DFS, capped at 2).
static int csp(const bg_t *g, uint32_t c, uint32_t d, uint8_t *on, int cap) {
  if (c == d) return 1;
  int t = 0; on[c] = 1;
  for (uint32_t i = 0; i < g->adjn[c] && t < cap; i++) {
    uint32_t w = g->adj[c][i];
    if (w < g->n && !on[w]) t += csp(g, w, d, on, cap);
  }
  on[c] = 0; return t;
}
// bg_t carries a precomputed length per segment; derive it so a fixture that forgets
// `slen` cannot leave a NULL array behind for the k-mer path.
static void tslen(bg_t *g, uint32_t nseg) {
  g->slen = malloc(nseg * sizeof(uint32_t));
  for (uint32_t i = 0; i < nseg; i++)
    g->slen[i] = g->seq[i] ? (uint32_t)strlen(g->seq[i]) : 0;
}

// ---- closed-syncmer predicate ----------------------------------------------
// A wrong verdict from the streaming selector is silent: it drops a junction or records a
// k-mer under a canonical key with only half its colors. Hence the brute-force oracle below.

// The oracle: recompute every s-mer of one k-mer from scratch, no rolling state.
static int sync_naive(const uint8_t *b, uint32_t k, uint32_t s) {
  uint32_t ns = k - s + 1;
  uint64_t smask = (1ULL << (2 * s)) - 1, best = ~0ULL, h0 = 0, hl = 0;
  for (uint32_t p = 0; p < ns; p++) {
    uint64_t fs = 0, rs = 0;
    for (uint32_t j = 0; j < s; j++) {
      fs = ((fs << 2) | b[p + j]) & smask;
      rs = (rs >> 2) | ((uint64_t)(3 - b[p + j]) << (2 * (s - 1)));
    }
    uint64_t h = kh_mix(fs < rs ? fs : rs);
    if (h < best) best = h;
    if (p == 0) h0 = h;
    if (p == ns - 1) hl = h;
  }
  return best == h0 || best == hl;
}

static void test_syncmer(void) {
  uint64_t rng = 0x9e3779b97f4a7c15ULL;
  #define RND() (rng = kh_mix(rng + 0x9e3779b97f4a7c15ULL))
  const uint32_t N = 200000;
  uint8_t *seq = malloc(N);
  for (uint32_t i = 0; i < N; i++) seq[i] = (uint8_t)(RND() & 3);

  const uint32_t ks[] = {31, 31, 31, 21, 15}, ss[] = {16, 12, 24, 11, 7};
  for (uint32_t t = 0; t < 5; t++) {
    uint32_t k = ks[t], s = ss[t], ns = k - s + 1;
    sync_t y; sync_init(&y, k, s);
    size_t sel = 0, tot = 0; int64_t last = -1; int64_t maxgap = 0;
    for (uint32_t i = 0; i < N; i++) {
      int got = sync_push(&y, seq[i]);
      if (i + 1 < k) { assert(!got); continue; }      // no k-mer yet: never selects
      int want = sync_naive(seq + i + 1 - k, k, s);
      assert(got == want);                            // (1) streaming == oracle
      tot++; sel += (size_t)got;
      if (got) { if (last >= 0 && (int64_t)i - last > maxgap) maxgap = (int64_t)i - last;
                 last = (int64_t)i; }
    }
    // (2) no junction goes dark: consecutive selections are never more than
    // k-s+1 apart, so any window of that many k-mers holds at least one.
    assert(maxgap > 0 && maxgap <= (int64_t)ns);
    // (3) density is 2/(k-s+1). 200k positions puts the sampling error well
    // inside 10%, so a predicate that drifted would be caught, not tolerated.
    double d = (double)sel / (double)tot, want = 2.0 / (double)ns;
    assert(fabs(d - want) < 0.10 * want);

    // (4) strand symmetry -- the property that makes a canonical key sound.
    for (uint32_t trial = 0; trial < 20000; trial++) {
      uint8_t fw[32], rc[32];
      for (uint32_t j = 0; j < k; j++) fw[j] = (uint8_t)(RND() & 3);
      for (uint32_t j = 0; j < k; j++) rc[j] = (uint8_t)(3 - fw[k - 1 - j]);
      assert(sync_naive(fw, k, s) == sync_naive(rc, k, s));
    }
  }

  // s == 0 is the filter switched off: every position selects, unconditionally.
  sync_t off; sync_init(&off, 31, 0);
  for (uint32_t i = 0; i < 1000; i++) assert(sync_push(&off, (int)(RND() & 3)));

  // A reset mid-stream must leave no state behind: the same suffix scanned after
  // a reset has to give the same verdicts as scanning it from the start. This is
  // what a non-ACGT base does, and what every caller does between walks.
  {
    uint32_t k = 31, s = 16;
    sync_t a, b2; sync_init(&a, k, s); sync_init(&b2, k, s);
    for (uint32_t i = 0; i < 5000; i++) (void)sync_push(&a, seq[i]);
    sync_reset(&a);
    for (uint32_t i = 0; i < 5000; i++) {
      int ga = sync_push(&a, seq[i]), gb = sync_push(&b2, seq[i]);
      assert(ga == gb);
    }
  }
  free(seq);
  #undef RND
}

// ref_trace: the guardian trace a scan produces when MULTI is decided by BRUTE
// FORCE -- csp counts simple paths and is valid on cyclic graphs. Independent of
// dominators, of reachability labels and of the last-edge reduction, so it is the
// only oracle that can arbitrate between them. Mirrors gpeb()'s scan exactly:
// cut at the edge entering the merge, re-anchor there, clear `seen`.
static uint32_t ref_trace(const bg_t *g, const uint32_t *p, uint32_t pn,
                          uint32_t *hd, uint32_t *tl, uint8_t *on, uint8_t *seen) {
  uint32_t ng = 0, anchor = p[0];
  memset(seen, 0, g->n); seen[anchor] = 1;
  for (uint32_t s = 1; s < pn; s++) {
    uint32_t v = p[s];
    int r = seen[v];
    int m = csp(g, anchor, v, on, 2) >= 2;
    if (r || m) {
      hd[ng] = p[s - 1]; tl[ng] = v; ng++;
      memset(seen, 0, g->n); anchor = v; seen[v] = 1;
    } else seen[v] = 1;
  }
  return ng;
}

// same_trace: guardian sets agree as (head, tail) sequences, in order.
static int same_trace(const gset_t *gs, const uint32_t *hd, const uint32_t *tl,
                      uint32_t ng) {
  if (gs->n != ng) return 0;
  for (uint32_t i = 0; i < ng; i++)
    if (gs->p[i].head != hd[i] || gs->p[i].tail != tl[i]) return 0;
  return 1;
}

// The .gfai's lazy W-lines against swl_build on the whole graph: every vertex, every key,
// the occupancy, and the chunk searches, on W-lines that run forward, against the reference,
// through an inversion (keys that go back), with an insertion, a repeat, and none of its nodes.
static void test_gfai(void) {
  const char *gp = "gpeb_gfai_selftest.gfa";
  FILE *f = fopen(gp, "w");
  assert(f);
  fprintf(f, "H\tVN:Z:1.1\n");
  for (int i = 0; i < 10; i++) fprintf(f, "S\t%d\t%.*s\n", i + 1, 7 + i, "ACGTACGTACGTACGTACGTAC");
  fprintf(f, "S\t11\tGGGG\nS\t12\tTT\nS\t13\tCCC\n");
  fprintf(f, "W\tREF\t0\tchr\t0\t115\t>1>2>3>4>5>6>7>8>9>10\n"
             "W\tA\t1\tchr\t0\t119\t>1>2>11>3>4>5>6>7>8>9>10\tXX:Z:tag\n"
             "W\tB\t1\tchr\t0\t115\t<10<9<8<7<6<5<4<3<2<1\n"
             "W\tC\t1\tchr\t0\t115\t>1>2>3<6<5<4>7>8>9>10\n"
             "W\tD\t1\tchr\t0\t40\t>1>2>3>4\n"
             "W\tD\t1\tchr\t60\t115\t>7>8>9>10\n"
             "W\tD\t1\tchr2\t0\t5\t>12>13\n"
             "W\tE\t1\tchr\t0\t70\t>1>2>3>2>3>4>5\n"
             "W\tF\t2\tchr\t0\t0\t");                      // no newline at the end, and no steps
  fclose(f);
  gfa_t *gfa = gfa_read(gp);
  assert(gfa && gfa->n_walk == 9);
  bg_quiet = 1; bg_nopadj = 1;
  bg_t *bg = bg_build(gfa, 1, 0, 0, NULL);
  uint32_t ch0 = gfai_chunk;
  for (uint32_t ch = 1; ch <= 4; ch++) {
    gfai_chunk = ch;
    assert(gfai_write(gp, gfa, bg, "REF") == 0);
    gfai_t *X = gfai_open(gp, "REF", 0, 0);
    assert(X && X->nw == 9 && X->ch == ch);
    char *full[1] = { "REF" }; uint64_t nskip = 0, nfull = 0;
    gfa_t *g2 = gfa_read_gfai(X, full, 1, NULL, 0, &nskip, &nfull);
    assert(g2 && nfull == 1 && g2->n_walk == gfa->n_walk && g2->n_seg == gfa->n_seg);
    bg_t *b2 = bg_build(g2, 1, 0, 0, NULL);
    sref_t *R = sref_build(gfa, bg, "REF"), *R2 = sref_build(g2, b2, "REF");
    assert(R && R2 && R->nref == R2->nref);
    X->R = R2;
    for (uint32_t c = 0; c < b2->ncpool; c++) {
      if (!strcmp(b2->cpool[c], "REF")) continue;
      scol_t M, Z; memset(&M, 0, sizeof M); memset(&Z, 0, sizeof Z);
      uint32_t nw = 0;
      for (uint32_t w = 0; w < b2->np; w++) if (b2->pname[w] == b2->cpool[c]) nw++;
      M.wl = calloc(nw, sizeof(swl_t)); Z.wl = calloc(nw, sizeof(swl_t));
      for (uint32_t w = 0; w < b2->np; w++) {
        if (b2->pname[w] != b2->cpool[c]) continue;
        assert(b2->pn[w] == 0 && bg->pn[w] == gfa->walk[w].n_v);
        swl_build(&M.wl[M.nwl++], bg, R, w);
        const gfai_w_t *xw = &X->w[b2->psrc[w]];
        if (xw->n) { swl_lazy(&Z.wl[Z.nwl], X, xw, w, &Z, Z.nwl); Z.lazy = 1; } else swl_build(&Z.wl[Z.nwl], b2, R2, w);
        Z.nwl++;
      }
      for (uint32_t wi = 0; wi < nw; wi++) {
        swl_t *a = &M.wl[wi], *z = &Z.wl[wi];
        assert(a->n == z->n && a->has == z->has && a->bp == z->bp && (!a->has || a->mono == z->mono));
        if (z->lazy && z->has && z->mono)                      // the chunked search == the binary search
          for (int64_t x = -3; x < 2 * 120; x++) assert(swl_upper(z, x) == swl_upper(a, x));
        if (z->lazy && z->has) {                               // swl_need reads every chunk a key range touches
          for (int64_t lo = 0; lo < 2 * 120; lo += 17) {
            swl_need(z, lo, lo + 30);
            for (uint32_t i = 0; i < a->n; i++)
              if (a->pf2[i] >= lo && a->pf2[i] <= lo + 30) assert(z->lv[i / ch] != NULL);
          }
        }
        for (uint32_t i = 0; i < a->n; i++) {
          assert(swl_v(z, i) == a->v[i]);
          if (a->has) assert(swl_k(z, i) == a->pf2[i]);
        }
      }
      if (Z.lazy) {                                            // every chunk read: the hash == scol_occ
        scol_occ(R, &M);
        for (uint32_t k = 0; k < R->nref; k++) {
          uint32_t ia = 0, iz = 0; int8_t sa = 0, sz = 0;
          int wa = occ_get(&M, k, &ia, &sa), wz = occ_get(&Z, k, &iz, &sz);
          assert(wa == wz && (wa < 0 || (ia == iz && sa == sz)));
        }
      }
      scol_free(&M); scol_free(&Z);
    }
    assert(!X->err && X->nload > 0);
    sref_free(R); sref_free(R2); bg_free(b2); gfa_destroy(g2); gfai_free(X);
  }
  gfai_chunk = ch0;
  // refused: keys of another reference, and a GFA that changed after its .gfai
  assert(!gfai_open(gp, "A", 0, 0));
  { struct stat sb; assert(!stat(gp, &sb));
    struct timespec tv[2] = { ST_ATIM(sb), ST_MTIM(sb) }; tv[1].tv_sec += 10;
    assert(!utimensat(AT_FDCWD, gp, tv, 0));
    assert(!gfai_open(gp, "REF", 0, 0)); }
  bg_quiet = 0; bg_nopadj = 0;
  char gi[64]; snprintf(gi, sizeof gi, "%s.gfai", gp);
  remove(gi); remove(gp);
  bg_free(bg); gfa_destroy(gfa);
  printf("[ok] .gfai: lazy W-lines (forward, reversed, inverted, insertion, repeat, off the reference, empty) read"
         " as swl_build reads them at chunks of 1-4 steps; occupancy and chunk searches exact; stale ones refused\n");
}

int main(void) {
  test_syncmer();
  test_gfai();
  // Aggregation: two walks with the SAME sample -> one color, continuous ranks.
  uint32_t d0[] = {1,2}, d1[] = {3}, d2[] = {3};
  uint32_t *adj2[4] = {d0,d1,d2,NULL}; uint32_t adjn2[4] = {2,1,1,0};
  uint32_t w0[] = {0,1,3}, w1[] = {0,2,3};   // each hits merge vertex 3 -> 1 pebble
  uint32_t *pv[2] = {w0,w1}; uint32_t pn2[2] = {3,3};
  char samp[] = "HG";                        // same pointer for both walks
  char *pnames[2] = {samp, samp};
  bg_t ag = {0}; ag.n = 4; ag.adj = adj2; ag.adjn = adjn2;
  ag.np = 2; ag.pv = pv; ag.pn = pn2; ag.pname = pnames;
  gset_t *g2 = gpeb(&ag, 1, 0, 0, 0.0, 0, 0); // fast mode, no coalescing/filter
  assert(g2->n == 2);
  assert(g2->p[0].pid != g2->p[1].pid);       // two different walks
  assert(g2->p[0].rank == 1 && g2->p[1].rank == 2); // ranks continuous across the color
  gpeb_free(g2);

  // Lever A/B: a cycle {1,2} on the walk 0->1->2->3. Vertex 1 has indeg 2 (from 0
  // and the back-edge 2->1), so --dag over-cuts there; but from anchor 0 there is
  // exactly one SIMPLE path to each vertex, so exact/scc/cache place NO pebble.
  uint32_t c0[] = {1}, c1[] = {2}, c2[] = {1,3};
  uint32_t *cadj[4] = {c0,c1,c2,NULL}; uint32_t cadjn[4] = {1,1,2,0};
  uint8_t cyc[4]; bg_t cg = {0}; cg.n = 4; cg.adj = cadj; cg.adjn = cadjn;
  scc_mark(&cg, cyc);
  assert(cyc[1] && cyc[2] && !cyc[0] && !cyc[3]);       // {1,2} is the only cycle
  uint32_t cw[] = {0,1,2,3}; uint32_t *cpv[1] = {cw}; uint32_t cpn[1] = {4};
  char csamp[] = "S"; char *cnm[1] = {csamp};
  cg.np = 1; cg.pv = cpv; cg.pn = cpn; cg.pname = cnm;
  gset_t *ge  = gpeb(&cg, 0, 0, 0, 0.0, 0, 0); assert(ge->n  == 0); gpeb_free(ge);  // exact
  gset_t *gd  = gpeb(&cg, 1, 0, 0, 0.0, 0, 0); assert(gd->n  == 1); gpeb_free(gd);  // --dag over-cuts

  // --acyclic: exact DAG trace by the forward single-sweep path count. Node 3 has
  // indeg 2 (from 1 and 2), but 2 is a SEPARATE source unreachable from the walk's
  // anchor 0, so the count from 0 is 1: exact/--acyclic place NO pebble, --fast over-cuts.
  uint32_t k0[]={1}, k1[]={3}, k2[]={3};
  uint32_t *kadj[4]={k0,k1,k2,NULL}; uint32_t kadjn[4]={1,1,1,0};
  bg_t kg={0}; kg.n=4; kg.adj=kadj; kg.adjn=kadjn;
  { swp_t *ks = swp_build(&kg); assert(ks); swp_reset(ks, 0);
    assert(swp_count(ks,1)==1 && swp_count(ks,3)==1);   // 3 reached once from 0
    swp_free(ks);
    uint32_t d0[]={1,2}, d1[]={3}, d2[]={3};            // diamond: 3 has 2 paths
    uint32_t *dadj[4]={d0,d1,d2,NULL}; uint32_t dadjn[4]={2,1,1,0};
    bg_t dg={0}; dg.n=4; dg.adj=dadj; dg.adjn=dadjn;
    swp_t *ds = swp_build(&dg); assert(ds); swp_reset(ds, 0);
    assert(swp_count(ds,3)==2); swp_free(ds); }         // clamped multi-path count

  // A cyclic graph has no topological order, so --acyclic refuses it.
  { uint32_t c0[]={1}, c1[]={0};
    uint32_t *cy2[2]={c0,c1}; uint32_t cyn2[2]={1,1};
    bg_t cy={0}; cy.n=2; cy.adj=cy2; cy.adjn=cyn2;
    assert(swp_build(&cy) == NULL);
    uint32_t cw[]={0,1}; uint32_t *cpv[1]={cw}; uint32_t cpn[1]={2};
    char csm[]="C"; char *cnm[1]={csm};
    cy.np=1; cy.pv=cpv; cy.pn=cpn; cy.pname=cnm;
    fprintf(stderr, "  (expected refusal) ");
    assert(gpeb(&cy, 0,0,0,0.0,0,1) == NULL); }
  printf("[ok] --acyclic refuses cyclic input instead of answering\n");
  uint32_t kw[]={0,1,3}; uint32_t *kpv[1]={kw}; uint32_t kpn[1]={3};
  char ksmp[]="K"; char *knm[1]={ksmp};
  kg.np=1; kg.pv=kpv; kg.pn=kpn; kg.pname=knm;
  gset_t *gct = gpeb(&kg, 0,0,0,0.0,0,1); assert(gct->n==0); gpeb_free(gct); // --acyclic exact
  gset_t *gcx = gpeb(&kg, 0,0,0,0.0,0,0); assert(gcx->n==0); gpeb_free(gcx); // --minimal agrees
  gset_t *gcd = gpeb(&kg, 1,0,0,0.0,0,0); assert(gcd->n==1); gpeb_free(gcd); // --fast over-cuts

  // gdom_build: one super-source dominator tree. On kg the sources are {0,2};
  // vertex 3 is entered from 1 and from the separate source 2, so idom(3)=S.
  // Validates gdom_build's idom + tin/tout.
  { uint32_t gi[4], gt[4], go[4];
    uint32_t *ph, *pr; pred_csr(&kg, &ph, &pr);
    gdom_build(&kg, ph, gi, gt, go);
    free(ph); free(pr);
    assert(gi[1]==0 && gi[0]==4 && gi[2]==4 && gi[3]==4); // idom: 1<-0; {0,2,3}<-S(=n=4)
    assert(gt[0] <= gt[1] && go[1] <= go[0]);             // 0 dominates 1
    assert(!(gt[0] <= gt[3] && go[3] <= go[0])); }        // 0 does NOT dominate 3

  // --acyclic on random DAGs: the sweep equals brute force from EVERY anchor, and
  // the whole trace (re-anchoring included) equals --minimal's.
  for (int trial = 0; trial < 2000; trial++) {
    uint32_t nn = 2 + rand() % 7;
    uint32_t **A = malloc(nn * sizeof(uint32_t*));
    uint32_t *An = calloc(nn, sizeof(uint32_t));
    for (uint32_t u = 0; u < nn; u++) {
      A[u] = malloc(nn * sizeof(uint32_t));
      for (uint32_t v = u + 1; v < nn; v++)          // u<v only => acyclic (DAG)
        if ((rand() & 3) == 0) A[u][An[u]++] = v;
    }
    bg_t r = {0}; r.n = nn; r.adj = A; r.adjn = An;
    swp_t *sw = swp_build(&r); assert(sw);           // random DAG: never refused
    uint8_t *on = calloc(nn, 1);
    for (uint32_t src = 0; src < nn; src++) {
      swp_reset(sw, src);
      for (uint32_t i = sw->topo[src]; i < nn; i++) { // targets in topological order
        uint32_t x = sw->order[i];
        int bc = csp(&r, src, x, on, 2); if (bc > 2) bc = 2; // csp clamps loosely
        assert(swp_count(sw, x) == (uint32_t)bc);           // vs min(2,#paths)
      }
    }
    free(on); swp_free(sw);
    { uint32_t maxlen = 2 * nn, pl = 0, cur = (uint32_t)(rand() % nn);
      uint32_t *pw = malloc(maxlen * sizeof(uint32_t)); pw[pl++] = cur;
      while (An[cur] > 0 && pl < maxlen) { cur = A[cur][rand() % An[cur]]; pw[pl++] = cur; }
      if (pl >= 2) {                                  // same trace as exact, cuts and all
        uint32_t *cpv[1]={pw}; uint32_t cpn[1]={pl}; char csm[]="C"; char *cnm[1]={csm};
        r.np=1; r.pv=cpv; r.pn=cpn; r.pname=cnm;
        gset_t *gc = gpeb(&r, 0,0,0,0.0,0,1), *gx = gpeb(&r, 0,0,0,0.0,0,0);
        assert(gc && gx && gc->n == gx->n);
        for (uint32_t a = 0; a < gc->n; a++)
          assert(gc->p[a].head == gx->p[a].head && gc->p[a].tail == gx->p[a].tail);
        gpeb_free(gc); gpeb_free(gx);
        r.np = 0; r.pv = NULL; r.pn = NULL; r.pname = NULL;
      }
      free(pw); }
    for (uint32_t u = 0; u < nn; u++) free(A[u]);
    free(A); free(An);
  }
  printf("[ok] --acyclic sweep == brute path count, and == --minimal trace, on random DAGs\n");

  // EXACT vs BRUTE FORCE on random CYCLIC walks: the reference is a simple-path count, so it
  // arbitrates the last-edge reduction rather than assuming it.
  { long checked = 0; srand(23);
    for (int trial = 0; trial < 4000; trial++) {
      uint32_t nn = 3 + rand() % 8;
      uint32_t **A = malloc(nn * sizeof(uint32_t*));
      uint32_t *An = calloc(nn, sizeof(uint32_t));
      for (uint32_t u = 0; u < nn; u++) { A[u] = malloc(nn * sizeof(uint32_t));
        for (uint32_t v = 0; v < nn; v++)
          if (v != u && (rand() & 3) == 0) A[u][An[u]++] = v; }   // cycles allowed
      bg_t r = {0}; r.n = nn; r.adj = A; r.adjn = An;
      uint32_t maxlen = 3 * nn;
      uint32_t *pw = malloc(maxlen * sizeof(uint32_t)), pl = 0, cur = (uint32_t)(rand() % nn);
      pw[pl++] = cur;
      while (An[cur] > 0 && pl < maxlen) { cur = A[cur][rand() % An[cur]]; pw[pl++] = cur; }
      if (pl >= 2) {
        uint32_t *pv1[1] = {pw}; uint32_t pn1[1] = {pl};
        char sm[] = "X"; char *nm1[1] = {sm};
        r.np = 1; r.pv = pv1; r.pn = pn1; r.pname = nm1;
        uint8_t *on = calloc(nn, 1), *seen = malloc(nn);
        uint32_t *hd = malloc(maxlen * sizeof(uint32_t)), *tl = malloc(maxlen * sizeof(uint32_t));
        uint32_t ng = ref_trace(&r, pw, pl, hd, tl, on, seen);
        gset_t *ee = gpeb(&r, 0,0,0,0.0,0,0);   // --minimal
        assert(same_trace(ee, hd, tl, ng));
        gpeb_free(ee);
        free(on); free(seen); free(hd); free(tl);
        checked++;
      }
      free(pw); for (uint32_t u = 0; u < nn; u++) free(A[u]); free(A); free(An);
    }
    printf("[ok] --minimal == brute-force simple-path trace on %ld random CYCLIC walks\n", checked);
  }

  // Default vs --minimal guardian SETS (paper, Proposition prop:containment). The
  // default cuts wherever --minimal does on a merge, but each mode empties S at its
  // own cuts, so on a cyclic graph a revisit cut can fire in one and not the other.
  // Fixed witness: 0->1->2->3->0 plus 4->2, walked 0,1,2,3,0. The default cuts at 2
  // (in-degree 2) and re-anchors there, so the return to 0 is no revisit; --minimal
  // sees that 4 is unreachable from 0, keeps going, and cuts on the revisit of 0.
  { uint32_t x0[]={1}, x1[]={2}, x2[]={3}, x3[]={0}, x4[]={2};
    uint32_t *xadj[5]={x0,x1,x2,x3,x4}; uint32_t xadjn[5]={1,1,1,1,1};
    uint32_t xw[]={0,1,2,3,0}; uint32_t *xpv[1]={xw}; uint32_t xpn[1]={5};
    char xsm[]="X"; char *xnm[1]={xsm};
    bg_t xg={0}; xg.n=5; xg.adj=xadj; xg.adjn=xadjn;
    xg.np=1; xg.pv=xpv; xg.pn=xpn; xg.pname=xnm;
    gset_t *lin = gpeb(&xg, 1,0,0,0.0,0,0), *dom = gpeb(&xg, 0,0,0,0.0,0,0);
    assert(lin->n == 1 && lin->p[0].pos == 1);   // (1,2): the merge, unreachable 2nd in-edge
    assert(dom->n == 1 && dom->p[0].pos == 3);   // (3,0): the revisit of the anchor
    gpeb_free(lin); gpeb_free(dom); }

  // ...and what DOES hold on any graph. Per walk: |minimal| <= |default|, and every
  // --minimal guardian the default lacks is a revisit onto its own anchor a, with
  // indeg(a) = 1, paired with its own default-only guardian earlier in that segment
  // such that every vertex after the partner's head, up to a, has in-degree 1.
  { long checked = 0, only = 0; srand(29);
    for (int trial = 0; trial < 4000; trial++) {
      uint32_t nn = 3 + rand() % 8;
      uint32_t **A = malloc(nn * sizeof(uint32_t*));
      uint32_t *An = calloc(nn, sizeof(uint32_t)), *ideg = calloc(nn, sizeof(uint32_t));
      for (uint32_t u = 0; u < nn; u++) { A[u] = malloc(nn * sizeof(uint32_t));
        for (uint32_t v = 0; v < nn; v++)
          if (v != u && (rand() & 3) == 0) { A[u][An[u]++] = v; ideg[v]++; } }
      bg_t r = {0}; r.n = nn; r.adj = A; r.adjn = An;
      uint32_t maxlen = 3 * nn;
      uint32_t *pw = malloc(maxlen * sizeof(uint32_t)), pl = 0, cur = (uint32_t)(rand() % nn);
      pw[pl++] = cur;
      while (An[cur] > 0 && pl < maxlen) { cur = A[cur][rand() % An[cur]]; pw[pl++] = cur; }
      if (pl >= 2) {
        uint32_t *pv1[1] = {pw}; uint32_t pn1[1] = {pl};
        char sm[] = "N"; char *nm1[1] = {sm};
        r.np = 1; r.pv = pv1; r.pn = pn1; r.pname = nm1;
        gset_t *lin = gpeb(&r, 1,0,0,0.0,0,0), *dom = gpeb(&r, 0,0,0,0.0,0,0);
        assert(dom->n <= lin->n);
        uint8_t *inl = calloc(pl, 1), *ind = calloc(pl, 1), *used = calloc(pl, 1);
        for (uint32_t i = 0; i < lin->n; i++) inl[lin->p[i].pos] = 1;
        for (uint32_t i = 0; i < dom->n; i++) ind[dom->p[i].pos] = 1;
        uint32_t alpha = 0;                      // --minimal's anchor position
        for (uint32_t g = 0; g + 1 < pl; g++) {
          if (!ind[g]) continue;
          uint32_t s = g + 1;
          if (!inl[g]) {
            only++;
            assert(pw[s] == pw[alpha] && ideg[pw[s]] == 1);
            int found = 0;
            for (int32_t x = (int32_t)s - 2; x >= (int32_t)alpha && !found; x--) {
              if (!inl[x] || ind[x] || used[x]) continue;
              int chain = 1;
              for (uint32_t j = (uint32_t)x + 2; j <= s; j++) if (ideg[pw[j]] != 1) { chain = 0; break; }
              if (chain) { used[x] = 1; found = 1; }
            }
            assert(found);
          }
          alpha = s;
        }
        free(inl); free(ind); free(used); gpeb_free(lin); gpeb_free(dom);
        checked++;
      }
      free(pw); for (uint32_t u = 0; u < nn; u++) free(A[u]); free(A); free(An); free(ideg);
    }
    assert(only > 0);                            // the random graphs do exercise the pairing
    printf("[ok] default vs --minimal on %ld random cyclic walks: |minimal| <= |default|; all %ld"
           " minimal-only guardians are anchor revisits paired with a default-only cut\n",
           checked, only);
  }

  // --early-cut, in EVERY mode. Fixed graph first: vertex 0 diverges (0->1, 0->5),
  // the walked branch is the forced corridor 1->2->3 (each in/out-degree 1), and 4
  // is the merge (in from 3 and 5). Late cuts on the edge entering the merge (3->4);
  // early slides the whole corridor back to the divergence edge (0->1). Exactly one
  // cut either way -- that equality is the invariant, checked at scale below.
  { uint32_t e0[]={1,5}, e1[]={2}, e2[]={3}, e3[]={4}, e4[]={6}, e5[]={4};
    uint32_t *eadj[7]={e0,e1,e2,e3,e4,e5,NULL}; uint32_t eadjn[7]={2,1,1,1,1,1,0};
    uint32_t ew[]={0,1,2,3,4};                 // walk down the corridor branch
    uint32_t *epv[1]={ew}; uint32_t epn[1]={5};
    char esm[]="E"; char *enm[1]={esm};
    bg_t eg={0}; eg.n=7; eg.adj=eadj; eg.adjn=eadjn;
    eg.np=1; eg.pv=epv; eg.pn=epn; eg.pname=enm;
    // (fast, acyclic) x (late, early); --minimal is the all-zero row.
    int modes[][2] = { {1,0}, {0,0}, {0,1} };
    for (size_t m = 0; m < sizeof modes / sizeof modes[0]; m++) {
      int fl = modes[m][0], ac = modes[m][1];
      gset_t *late  = gpeb(&eg, fl,0,0,0.0, 0,ac);
      gset_t *earl  = gpeb(&eg, fl,0,0,0.0, 1,ac);
      assert(late->n == 1 && earl->n == 1);
      assert(late->p[0].head == 3 && late->p[0].tail == 4);   // at the merge
      assert(earl->p[0].head == 0 && earl->p[0].tail == 1);   // at the divergence
      gpeb_free(late); gpeb_free(earl);
    }
  }

  // Randomized: the --early-cut INVARIANT, on random DAGs and random CYCLIC graphs,
  // in every mode. Sliding may only move a guardian BACKWARD along a forced
  // corridor -- it must never change which merges are cut. So with min_gap=0 and
  // min_complexity=0 (both position-sensitive, hence excluded): the guardian count
  // is identical, and pairwise (pebbles stay in walk order) each early guardian
  // sits at or before its late counterpart.
  { long checked = 0; srand(23);
    for (int trial = 0; trial < 4000; trial++) {
      uint32_t nn = 3 + rand() % 9;
      int cyclic = trial & 1;                  // half DAGs, half cyclic
      uint32_t **A = malloc(nn * sizeof(uint32_t*));
      uint32_t *An = calloc(nn, sizeof(uint32_t));
      for (uint32_t u = 0; u < nn; u++) { A[u] = malloc(nn * sizeof(uint32_t));
        for (uint32_t v = cyclic ? 0 : u + 1; v < nn; v++)
          if (v != u && (rand() & 3) == 0) A[u][An[u]++] = v; }
      bg_t r = {0}; r.n = nn; r.adj = A; r.adjn = An;
      uint32_t maxlen = 3 * nn;
      uint32_t *pw = malloc(maxlen * sizeof(uint32_t)), pl = 0, cur = (uint32_t)(rand() % nn);
      pw[pl++] = cur;
      while (An[cur] > 0 && pl < maxlen) { cur = A[cur][rand() % An[cur]]; pw[pl++] = cur; }
      if (pl >= 2) {
        uint32_t *pv1[1] = {pw}; uint32_t pn1[1] = {pl};
        char sm[] = "V"; char *nm1[1] = {sm};
        r.np = 1; r.pv = pv1; r.pn = pn1; r.pname = nm1;
        // (fast, acyclic): --acyclic is a DAG-only rule, so it is compared
        // against itself (late vs early), never across modes.
        int modes[][2] = { {1,0}, {0,0}, {0,1} };
        for (size_t m = 0; m < sizeof modes / sizeof modes[0]; m++) {
          int fl=modes[m][0], ac=modes[m][1];
          if (ac && cyclic) continue;          // --acyclic undefined on cycles
          gset_t *late = gpeb(&r, fl,0,0,0.0, 0,ac);
          gset_t *earl = gpeb(&r, fl,0,0,0.0, 1,ac);
          assert(late->n == earl->n);          // same merges cut, always
          for (uint32_t a = 0; a < late->n; a++)
            assert(earl->p[a].pos <= late->p[a].pos);   // only ever slides backward
          gpeb_free(late); gpeb_free(earl);
        }
        checked++;
      }
      free(pw); for (uint32_t u = 0; u < nn; u++) free(A[u]); free(A); free(An);
    }
    printf("[ok] --early-cut invariant (count equal, pos monotone) on %ld random DAG+cyclic graphs, all modes\n", checked);
  }

  // k-mer LOCUS uniqueness (--kmer-stats). Two alleles s2/s3 between a shared
  // upstream s1 and a shared downstream s4; walk Z re-plays X's exact junction
  // sequence (s2 | s4[:30]) at a different place, as s5|s6. So every k-mer of X's
  // window occurs at two loci, while Y's window keeps the k-mers that touch s3.
  // With span=30, k=21: Y contributes 10 k-mers wholly inside s3 plus 20 crossing
  // the s3|s4 junction = 30 single-locus k-mers, and X contributes none -- so
  // exactly one of the two windows is unambiguously indexable. This also pins the
  // locus semantics: the shared s4 prefix must count as ambiguous (s4 and s6), and
  // `absent` must be 0.
  { char q1[61], q2[31], q3[31], q4[61], q5[31], q6[31];
    unsigned rs = 12345;                              // tiny LCG: reproducible seq
    #define RB() ("ACGT"[(rs = rs * 1103515245u + 12345u) >> 24 & 3])
    for (int i = 0; i < 60; i++) q1[i] = RB(); q1[60] = 0;
    for (int i = 0; i < 30; i++) q2[i] = RB(); q2[30] = 0;
    for (int i = 0; i < 30; i++) q3[i] = RB(); q3[30] = 0;
    for (int i = 0; i < 60; i++) q4[i] = RB(); q4[60] = 0;
    #undef RB
    // The two alleles must differ at the base adjacent to the junction, else the
    // k-mer starting there is the same in both and the count is off by one -- a
    // property of the test sequence, not of the code. Forcing it keeps the
    // expected value exact instead of sequence-dependent.
    q3[29] = (q2[29] == 'A') ? 'C' : 'A';
    memcpy(q5, q2, 31);                               // s5 == s2
    memcpy(q6, q4, 30); q6[30] = 0;                   // s6 == s4[:30]
    char *sq[6]  = {q1, q2, q3, q4, q5, q6};
    char *snm[6] = {(char*)"1",(char*)"2",(char*)"3",(char*)"4",(char*)"5",(char*)"6"};
    uint32_t v0[]={2}, v2[]={6}, v4[]={6}, v8[]={10}; // forward vertices: seg*2
    uint32_t *kadj2[12]={v0,NULL,v2,NULL,v4,NULL,NULL,NULL,v8,NULL,NULL,NULL};
    uint32_t kadjn2[12]={1,0,1,0,1,0,0,0,1,0,0,0};
    uint32_t wX[]={0,2,6}, wY[]={0,4,6}, wZ[]={8,10};
    uint32_t *kpv2[3]={wX,wY,wZ}; uint32_t kpn2[3]={3,3,2};
    char nX[]="X", nY[]="Y", nZ[]="Z"; char *knm2[3]={nX,nY,nZ};
    bg_t kb={0}; kb.n=12; kb.adj=kadj2; kb.adjn=kadjn2;
    kb.np=3; kb.pv=kpv2; kb.pn=kpn2; kb.pname=knm2; kb.seq=sq; kb.sname=snm;
    tslen(&kb, 6);
    // vertex 6 (=s4) has in-degree 2 -> one guardian on each allele's in-edge.
    gset_t *kg2 = gpeb(&kb, 1, 0, 30, 0.0, 0, 0);
    assert(kg2->n == 2);
    kstat_t R = kstat_run(&kb, kg2, 30, 21, 0);
    assert(R.nwin == 2);          // the two windows are distinct
    assert(R.absent == 0);        // every window k-mer located in its own walk
    assert(R.tot_km == 80);       // 2 windows x (2*30 - 21 + 1) k-mers
    assert(R.uniq_km == 30);      // only Y's s3-touching k-mers pin one locus
    assert(R.win_uniq == 1);      // ...so only Y's window is indexable
    gpeb_free(kg2);
  }

  // Color sets WIDER THAN ONE WORD. A chain of NB bubbles with NCOL walks, walk w
  // taking allele bit i of w at bubble i: every junction k-mer's color set is then
  // a known NCOL-bit pattern, and with NCOL > 64 some of them carry a color a
  // single uint64_t mask could not even represent. Checks the word/bit arithmetic
  // in the build, the .pex record width, and the set-bit iteration the per-hit
  // vote loop relies on.
  { enum { NB = 7, NCOL = 100, SEGLEN = 40, KK = 21, NSEG = 3 * NB + 1 };
    char sq[NSEG][SEGLEN + 1], snb[NSEG][8], *sqp[NSEG], *snm[NSEG];
    unsigned rs = 987654321u;
    #define RB() ("ACGT"[(rs = rs * 1103515245u + 12345u) >> 24 & 3])
    for (uint32_t s2 = 0; s2 < NSEG; s2++) {
      for (int i = 0; i < SEGLEN; i++) sq[s2][i] = RB();
      sq[s2][SEGLEN] = 0; sqp[s2] = sq[s2];
      snprintf(snb[s2], sizeof snb[s2], "%u", s2); snm[s2] = snb[s2];
    }
    #undef RB
    // The alleles must differ at the base ADJACENT to the junction, else the k-mers
    // covering only that tail are shared, land at two loci and drop out -- the same
    // off-by-one the locus test above pins, here it would silently shrink the index.
    for (uint32_t i = 0; i < NB; i++) {
      char *A = sq[3*i+1], *B = sq[3*i+2];
      if (A[SEGLEN-1] == B[SEGLEN-1]) B[SEGLEN-1] = (A[SEGLEN-1] == 'A') ? 'C' : 'A';
    }
    uint32_t nv = 2 * NSEG;
    uint32_t **adj = calloc(nv, sizeof(uint32_t *)), *adjn = calloc(nv, sizeof(uint32_t));
    uint32_t (*ab)[2] = calloc(nv, sizeof(*ab));
    for (uint32_t i = 0; i < NB; i++) {              // X_i -> {A_i,B_i} -> X_{i+1}
      uint32_t X = 6*i, A = 6*i+2, B = 6*i+4, X2 = 6*i+6;
      ab[X][0] = A;  ab[X][1] = B;  adj[X] = ab[X];  adjn[X] = 2;
      ab[A][0] = X2; adj[A] = ab[A]; adjn[A] = 1;
      ab[B][0] = X2; adj[B] = ab[B]; adjn[B] = 1;
    }
    uint32_t wv[NCOL][2*NB+1], *pv2[NCOL], pn3[NCOL];
    char cnm[NCOL][8], *pnm[NCOL];
    for (uint32_t w = 0; w < NCOL; w++) {
      uint32_t j = 0;
      for (uint32_t i = 0; i < NB; i++) {
        wv[w][j++] = 6*i;
        wv[w][j++] = ((w >> i) & 1) ? 6*i+4 : 6*i+2;
      }
      wv[w][j++] = 6*NB;
      pv2[w] = wv[w]; pn3[w] = j;
      snprintf(cnm[w], sizeof cnm[w], "c%u", w); pnm[w] = cnm[w];
    }
    bg_t wb = {0}; wb.n = nv; wb.adj = adj; wb.adjn = adjn;
    wb.np = NCOL; wb.pv = pv2; wb.pn = pn3; wb.pname = pnm; wb.seq = sqp; wb.sname = snm;
    tslen(&wb, nv / 2);
    gset_t *wg = gpeb(&wb, 1, 0, SEGLEN, 0.0, 0, 0);
    assert(wg->n == NCOL * NB);                      // one merge per bubble per walk
    uint64_t exp_pre[2*NB][2];                       // the only color sets that exist
    memset(exp_pre, 0, sizeof exp_pre);
    for (uint32_t i = 0; i < NB; i++)
      for (uint32_t w = 0; w < NCOL; w++) cm_set(exp_pre[2*i + ((w >> i) & 1)], w);
    pidx_t *ix = pidx_build(&wb, wg, SEGLEN, KK, 0, 0, 0, 1, 0);
    assert(ix && ix->nw == 2 && ix->nc == NCOL);     // 100 colors -> two words
    assert(ix->nkmer == (size_t)2 * NB * (KK - 1));  // k-1 junction k-mers per allele

    // --threads must not change the answer. This graph's bubbles share segments
    // across up to NCOL walks, so a 100-color, 8-thread build hits the exact
    // contention --threads exists to handle: many threads racing the same slot.
    for (int nt = 2; nt <= 8; nt *= 2) {
      pidx_t *tx = pidx_build(&wb, wg, SEGLEN, KK, 0, 0, 0, nt, 0);
      assert(tx && tx->nw == ix->nw && tx->nc == ix->nc && tx->nkmer == ix->nkmer);
      size_t matched = 0;
      for (size_t s2 = 0; s2 <= ix->km->mask; s2++) {
        if (ix->km->st[s2] != KH_IDX) continue;
        size_t t = kh_slot(tx->km, ix->km->key[s2], 0);
        assert(t != (size_t)-1 && tx->km->st[t] == KH_IDX);
        assert(cm_eq(tx->cmask + t * tx->nw, ix->cmask + s2 * ix->nw, ix->nw));
        matched++;
      }
      assert(matched == ix->nkmer);
      pidx_free(tx);
    }
    printf("[ok] --threads 2/4/8 == --threads 1, same (key,mask) set over %d colors\n", NCOL);

    // --syncmer end to end. Two things have to hold and neither is implied by the
    // unit tests above: the filtered index must be a strict SUBSET of the
    // unfiltered one carrying identical masks (so the filter drops records and
    // changes nothing else), and every allele must survive it -- an allele whose
    // k-mers were all filtered away is a guardian gone dark, which is the one
    // failure that would make the index quietly less sensitive rather than smaller.
    {
      uint32_t sk = KK / 2;
      pidx_t *sx = pidx_build(&wb, wg, SEGLEN, KK, 0, 0, 0, 1, sk);
      assert(sx && sx->sk == sk && sx->nc == ix->nc && sx->nw == ix->nw);
      assert(sx->nkmer > 0 && sx->nkmer < ix->nkmer);
      uint8_t *live = calloc(2 * NB, 1);
      for (size_t s2 = 0; s2 <= sx->km->mask; s2++) {
        if (sx->km->st[s2] != KH_IDX) continue;
        size_t t = kh_slot(ix->km, sx->km->key[s2], 0);
        assert(t != (size_t)-1 && ix->km->st[t] == KH_IDX);          // subset
        assert(cm_eq(sx->cmask + s2 * sx->nw, ix->cmask + t * ix->nw, ix->nw));
        uint8_t bs[32];                                              // and selected
        for (uint32_t z = 0; z < KK; z++)
          bs[z] = (uint8_t)((sx->km->key[s2] >> (2 * (KK - 1 - z))) & 3);
        assert(sync_naive(bs, KK, sk));
        for (int e = 0; e < 2 * NB; e++)
          if (cm_eq(sx->cmask + s2 * sx->nw, exp_pre[e], sx->nw)) live[e] = 1;
      }
      for (int e = 0; e < 2 * NB; e++) assert(live[e]);               // no allele lost
      printf("[ok] --syncmer s=%u: %zu of %zu records kept (%.3f), all %d alleles"
             " still represented\n", sk, sx->nkmer, ix->nkmer,
             (double)sx->nkmer / (double)ix->nkmer, 2 * NB);
      free(live); pidx_free(sx);
    }

    // --subsample. Three properties, none of them implied by "it drops records":
    // rate 1 is an identity, a lower rate NESTS inside a higher one at a fixed
    // seed, and the masks of whatever survives are untouched. The ladder is only
    // an ablation if all three hold -- otherwise a kink in the curve could be the
    // sampler rather than the index.
    {
      pidx_t *one = pidx_build(&wb, wg, SEGLEN, KK, 0, 0, 0, 1, 0);
      pidx_t *h   = pidx_build(&wb, wg, SEGLEN, KK, 0, 0, 0, 1, 0);
      pidx_t *q   = pidx_build(&wb, wg, SEGLEN, KK, 0, 0, 0, 1, 0);
      size_t n0 = ix->nkmer;
      assert(pex_subsample(one, 1.0, 7) == n0);         // rate >= 1 is an identity
      size_t nh = pex_subsample(h, 0.5,  7);
      size_t nq = pex_subsample(q, 0.25, 7);
      assert(nq < nh && nh < n0);
      // 280 records, so the sampling s.d. at r=0.5 is ~8.4 and these bounds are
      // over 5 sigma out. They catch a broken threshold, not an unlucky draw.
      assert(nh > n0 * 35 / 100 && nh < n0 * 65 / 100);
      assert(nq > n0 * 13 / 100 && nq < n0 * 38 / 100);
      for (size_t s2 = 0; s2 <= q->km->mask; s2++) {
        if (q->km->st[s2] != KH_IDX) continue;
        size_t a = kh_slot(h->km, q->km->key[s2], 0);
        assert(a != (size_t)-1 && h->km->st[a] == KH_IDX);      // 0.25 nests in 0.5
        size_t b2 = kh_slot(ix->km, q->km->key[s2], 0);
        assert(b2 != (size_t)-1 && ix->km->st[b2] == KH_IDX);   // ...and in the whole
        assert(cm_eq(q->cmask + s2 * q->nw, ix->cmask + b2 * ix->nw, ix->nw));
      }
      printf("[ok] --subsample: %zu -> %zu -> %zu records, nested, masks intact\n",
             n0, nh, nq);
      pidx_free(one); pidx_free(h); pidx_free(q);
    }

    int wide = 0;
    for (size_t s2 = 0; s2 <= ix->km->mask; s2++) {
      if (ix->km->st[s2] != KH_IDX) continue;
      const uint64_t *m = ix->cmask + s2 * ix->nw;
      int known = 0;
      for (int e = 0; e < 2*NB; e++) if (cm_eq(m, exp_pre[e], ix->nw)) { known = 1; break; }
      assert(known);
      if (m[1]) wide = 1;                            // a color >= 64 is actually set
    }
    assert(wide);

    // Set-bit iteration (the per-hit vote loop) must agree with a full color scan.
    for (int e = 0; e < 2*NB; e++) {
      uint32_t a = 0, b = 0;
      for (uint32_t c = 0; c < NCOL; c++) if (cm_get(exp_pre[e], c)) a += c + 1;
      for (uint32_t z = 0; z < 2; z++)
        for (uint64_t bt = exp_pre[e][z]; bt; bt &= bt - 1)
          b += (z << 6) + (uint32_t)__builtin_ctzll(bt) + 1;
      assert(a == b);
    }

    const char *tp = "gpeb_pex_selftest.tmp";        // .pex carries nw words/record
    assert(pex_write(tp, ix, NULL) == 0);
    pidx_t *rd = pex_read(tp);
    assert(rd && rd->nw == ix->nw && rd->nc == ix->nc && rd->nkmer == ix->nkmer);
    assert(ix->km->loc != NULL);    // the built table located k-mers, so it has loc
    assert(rd->km->loc == NULL);    // a loaded one never will -- 8 bytes/slot saved
    for (size_t s2 = 0; s2 <= ix->km->mask; s2++) {
      if (ix->km->st[s2] != KH_IDX) continue;
      size_t t = kh_slot(rd->km, ix->km->key[s2], 0);
      assert(t != (size_t)-1 && rd->km->st[t] == KH_IDX);
      assert(cm_eq(rd->cmask + t * rd->nw, ix->cmask + s2 * ix->nw, ix->nw));
    }
    // Panel geometry is exactly predictable here: a color's k-mer set is one allele
    // per bubble, KK-1 junction k-mers each, so two colors agreeing on b of the NB
    // bubbles share 20b of 140. The closest possible neighbour differs in one bit,
    // giving 120/160 = 0.75, and every color below NCOL has such a neighbour.
    { uint32_t *nn = malloc(NCOL * sizeof(uint32_t));
      double *jc = malloc(NCOL * sizeof(double));
      panel_nn(ix, nn, jc);
      for (uint32_t c = 0; c < NCOL; c++) {
        assert(fabs(jc[c] - 0.75) < 1e-9);              // one differing bubble
        assert(__builtin_popcount(c ^ nn[c]) == 1);     // ...and that is the neighbour
      }
      free(nn); free(jc);
      printf("[ok] panel geometry: nearest neighbour exact over %d colors\n", NCOL); }

    pidx_free(rd); remove(tp);

    // FASTQ records are delimited by POSITION, not by a leading '@': a quality line may
    // start with '@' too.
    const char *fq = "gpeb_fq_selftest.tmp", *ft = "gpeb_fqout_selftest.tmp";
    FILE *qf = fopen(fq, "w");
    assert(qf);
    for (int q = 0; q < 3; q++)                      // quality is all '@'
      fprintf(qf, "@r%d\nACGTACGTACGTACGTACGTACGT\n+\n@@@@@@@@@@@@@@@@@@@@@@@@\n", q);
    fclose(qf);
    FILE *of = fopen(ft, "w");
    assert(of);
    run_classify(ix, fq, SET_ALL, 5, -1, of);
    fclose(of);
    int nl = 0;
    FILE *rf = fopen(ft, "r");
    assert(rf);
    for (int ch; (ch = fgetc(rf)) != EOF; ) if (ch == '\n') nl++;
    fclose(rf);
    assert(nl == 4);                                 // header + one line per read
    remove(fq); remove(ft);

    pidx_free(ix); gpeb_free(wg);
    free(adj); free(adjn); free(ab);
    printf("[ok] %d-color index: masks exact over 2 words, .pex round-trip clean\n", NCOL);
  }

  // Complexity: tandem/homopolymer windows score low, diverse sequence high.
  const char *rep = "AACCCTAACCCTAACCCTAACCCTAACCCT"; // telomere-like
  const char *div = "ACGTTGCAAGTCCTGACTTAGGCATCGATG"; // diverse
  assert(kmer_complexity(rep, strlen(rep)) < 0.5);
  assert(kmer_complexity(div, strlen(div)) > 0.9);

  // G_P adjacency: diamond walks {0,1,3},{0,2,3} -> 0:{1,2} 1:{3} 2:{3} 3:{}.
  // vertex 0 is shared by both walks with distinct successors (dedup keeps both);
  // vertex 3 is a sink (never a source). Edge 4+ (no walk) must not appear.
  uint32_t pw0[] = {0,1,3}, pw1[] = {0,2,3};
  uint32_t *ppv[2] = {pw0, pw1}; uint32_t ppn[2] = {3,3};
  bg_t pg = {0}; pg.n = 5; pg.np = 2; pg.pv = ppv; pg.pn = ppn;
  build_padj(&pg);
  assert(pg.adjn[0]==2 && pg.adj[0][0]==1 && pg.adj[0][1]==2);
  assert(pg.adjn[1]==1 && pg.adj[1][0]==3);
  assert(pg.adjn[2]==1 && pg.adj[2][0]==3);
  assert(pg.adjn[3]==0 && pg.adjn[4]==0);
  for (uint32_t v = 0; v < pg.n; v++) free(pg.adj[v]);
  free(pg.adj); free(pg.adjn);

  // cover_window against brute force. Random records (a mask over 9 colours, a weight) and
  // a reference colour that is always in the graph. The greedy must add up exactly -- the
  // gains it reports are the marginal gains, they never increase, and base coverage plus
  // gains is the coverage -- keep at least 1 - 1/e of the best 3 colours (the bound a
  // monotone submodular objective guarantees; in practice it is almost always optimal),
  // and obey the stop rule and the tie rule.
  { srand(99); long ntest = 0, optimal = 0; const uint32_t NC = 9;
    for (int trial = 0; trial < 600; trial++) {
      uint32_t n = 20 + (uint32_t)(rand() % 60);
      uint64_t *cm = malloc(n * sizeof(uint64_t)); uint32_t *sl = malloc(n * sizeof(uint32_t)), *wt = malloc(n * sizeof(uint32_t));
      uint64_t total = 0;
      for (uint32_t i = 0; i < n; i++) {
        uint64_t m = 0; while (!m) m = (uint64_t)(rand() & ((1 << NC) - 1));
        cm[i] = m; sl[i] = i; wt[i] = 1 + (uint32_t)(rand() % 5); total += wt[i];
      }
      uint64_t base = 1;                                      // colour 0 is the reference
      uint8_t allowed[9]; for (uint32_t c = 0; c < NC; c++) allowed[c] = 1;
      cover_par_t P = { 3, 0.0, 0.0, CR_GREEDY }; cover_stat_t st; uint32_t pick[3]; uint64_t pg[3];
      uint32_t k = cover_window(cm, 1, NC, sl, wt, n, &base, allowed, NULL, &P, pick, pg, &st);
      uint64_t cov0 = 0;
      for (uint32_t i = 0; i < n; i++) if (cm[i] & base) cov0 += wt[i];
      assert(st.hits == total && st.unc0 == total - cov0);
      uint64_t cg = cov0, S = base;
      for (uint32_t j = 0; j < k; j++) {
        if (j) assert(pg[j] <= pg[j - 1]);                    // gains never increase
        uint64_t marg = 0;                                    // the marginal gain, recomputed
        for (uint32_t i = 0; i < n; i++) if (!(cm[i] & S) && (cm[i] >> pick[j] & 1)) marg += wt[i];
        assert(marg == pg[j]);
        S |= 1ULL << pick[j]; cg += marg;
        for (uint32_t c = 0; c < NC; c++) if (c != pick[j]) { uint64_t m2 = 0;   // nobody beat the pick
          for (uint32_t i = 0; i < n; i++) if (!(cm[i] & (S & ~(1ULL << pick[j])) ) && (cm[i] >> c & 1)) m2 += wt[i];
          if (j == 0) assert(m2 <= pg[0]); }
      }
      assert(st.unc_end == total - cg);
      uint64_t best = 0;                                      // brute force over every 3 non-reference colours
      for (uint32_t a = 1; a < NC; a++) for (uint32_t b = a; b < NC; b++) for (uint32_t c = b; c < NC; c++) {
        uint64_t T2 = base | 1ULL << a | 1ULL << b | 1ULL << c, v = 0;
        for (uint32_t i = 0; i < n; i++) if (cm[i] & T2) v += wt[i];
        if (v > best) best = v;
      }
      assert(cg <= best);
      assert((double)(cg - cov0) >= 0.632 * (double)(best - cov0) - 1e-9);   // 1 - 1/e of the gain
      optimal += cg == best; ntest++;
      // stop rule: a floor of 1 hit more than the first gain allows no pick at all
      cover_par_t P2 = { 3, (double)pg[0] + 1.0, 0.0, CR_GREEDY };
      if (k) assert(cover_window(cm, 1, NC, sl, wt, n, &base, allowed, NULL, &P2, pick, pg, &st) == 0);
      free(cm); free(sl); free(wt);
    }
    // ties go to the colour picked in the previous window: records hit only by colours 1 and 2 alike
    { uint64_t cm2[2] = { 1ULL << 1 | 1ULL << 2, 1ULL << 1 | 1ULL << 2 }; uint32_t sl2[2] = {0, 1}, wt2[2] = {3, 3};
      uint64_t base2 = 1; uint8_t prev[9] = {0}, allowed[9]; for (int c = 0; c < 9; c++) allowed[c] = 1;
      cover_par_t P = { 1, 0.0, 0.0, CR_GREEDY }; cover_stat_t st; uint32_t pick[1]; uint64_t pg[1];
      assert(cover_window(cm2, 1, 9, sl2, wt2, 2, &base2, allowed, NULL, &P, pick, pg, &st) == 1 && pick[0] == 1);
      prev[2] = 1;
      assert(cover_window(cm2, 1, 9, sl2, wt2, 2, &base2, allowed, prev, &P, pick, pg, &st) == 1 && pick[0] == 2);
      // the control picks the colours with the most hits each, overlap and all
      uint64_t cm3[3] = { 1ULL << 1 | 1ULL << 2, 1ULL << 1 | 1ULL << 2, 1ULL << 3 }; uint32_t sl3[3] = {0, 1, 2}, wt3[3] = {5, 5, 4};
      cover_par_t Pc = { 2, 0.0, 0.0, CR_SUPPORT }, Pg = { 2, 0.0, 0.0, CR_GREEDY }; uint32_t pk2[2]; uint64_t pg2[2];
      assert(cover_window(cm3, 1, 9, sl3, wt3, 3, &base2, allowed, NULL, &Pc, pk2, pg2, &st) == 2 && pk2[0] == 1 && pk2[1] == 2);
      assert(st.unc_end == 4);                                // colour 3's records stay uncovered: the control overlaps
      assert(cover_window(cm3, 1, 9, sl3, wt3, 3, &base2, allowed, NULL, &Pg, pk2, pg2, &st) == 2 && pk2[0] == 1 && pk2[1] == 3);
      assert(st.unc_end == 0); }
    printf("[ok] cover_window: %ld random windows, gains exact and non-increasing, >= 1-1/e of the best 3 colours"
           " (optimal in %ld), stop and tie rules\n", ntest, optimal); }

  // --set-call mosaic, end to end, on a world whose answer is known. A reference and
  // five haplotypes over 130 SNP sites, one every 60 bp. The donor is diploid: copy 1
  // is H1 up to site 65 and H3 after it (ONE recombination), copy 2 is H5 throughout.
  // Error-free 100 bp reads from both copies are placed by the index alone, and the
  // painting must be {H1,H5} on the left and {H3,H5} on the right; the spliced tracks
  // must be one edge-valid walk each, as long as the reference.
  { enum { NS = 130, SG = 59, NH = 6, RL = 100, BLK = 500 };
    const int LEN = (NS + 1) * SG + NS;
    static char sseq[NS + 1][SG + 1]; char refal[NS], altal[NS]; uint8_t al[NH][NS];
    const char *nm[NH] = { "REF", "H1", "H2", "H3", "H4", "H5" };
    unsigned rs = 2024;
    #define RB3() ("ACGT"[(rs = rs * 1103515245u + 12345u) >> 24 & 3])
    for (int i = 0; i <= NS; i++) { for (int j = 0; j < SG; j++) sseq[i][j] = RB3(); sseq[i][SG] = 0; }
    for (int i = 0; i < NS; i++) {
      refal[i] = RB3(); altal[i] = "ACGT"[(strchr("ACGT", refal[i]) - "ACGT" + 1) & 3];
    }
    for (int h = 0; h < NH; h++)
      for (int i = 0; i < NS; i++) { rs = rs * 1103515245u + 12345u; al[h][i] = h ? (uint8_t)((rs >> 24) & 1) : 0; }
    #undef RB3
    const char *gp = "gpeb_mosaic_selftest.gfa";
    FILE *gf = fopen(gp, "w");
    assert(gf);
    fprintf(gf, "H\tVN:Z:1.1\n");
    for (int i = 0; i <= NS; i++) fprintf(gf, "S\ts%d\t%s\n", i, sseq[i]);
    for (int i = 0; i < NS; i++) {
      fprintf(gf, "S\ta%d\t%c\nS\tb%d\t%c\n", i, refal[i], i, altal[i]);
      fprintf(gf, "L\ts%d\t+\ta%d\t+\t0M\nL\ts%d\t+\tb%d\t+\t0M\n"
                  "L\ta%d\t+\ts%d\t+\t0M\nL\tb%d\t+\ts%d\t+\t0M\n", i, i, i, i, i, i + 1, i, i + 1);
    }
    for (int h = 0; h < NH; h++) {
      fprintf(gf, "W\t%s\t0\tchr\t0\t%d\t>s0", nm[h], LEN);
      for (int i = 0; i < NS; i++) fprintf(gf, ">%c%d>s%d", al[h][i] ? 'b' : 'a', i, i + 1);
      fprintf(gf, "\n");
    }
    fclose(gf);
    gfa_t *gfa = gfa_read(gp);
    assert(gfa && gfa->n_walk == NH);
    bg_t *bg = bg_build(gfa, 1, 0, 0, NULL);
    gset_t *gs = gpeb(bg, 1, GSPAN, GSPAN, 0.0, 0, 0);
    pidx_t *ix = pidx_build(bg, gs, GSPAN, 21, 0, 0, 0, 1, 0);
    assert(ix && ix->nkmer > 0);
    assert(pidx_set_blocks(ix, gfa, bg, "REF", BLK) == 0);
    assert(ix->blk_len == (uint64_t)LEN && !strcmp(ix->blk_ref, "REF"));
    // Positions: every window of the reference that is an index record must carry a block
    // within the k-mer's own span, [p, p+20] -- the check that pins the projection.
    { char *rseq = malloc((size_t)LEN + 1); int o = 0;
      for (int i = 0; i <= NS; i++) { memcpy(rseq + o, sseq[i], SG); o += SG; if (i < NS) rseq[o++] = refal[i]; }
      kroll_t q; kr_init(&q, 21, 0);
      int seen = 0;
      for (int p = 0; p < LEN; p++) {
        if (kr_push(&q, rseq[p]) != 1) continue;
        size_t sl = kh_slot(ix->km, kr_key(&q), 0);
        if (sl == (size_t)-1 || ix->km->st[sl] != KH_IDX || ix->blk[sl] == BLK_NONE) continue;
        int lo = (p - 20) / BLK, hi = p / BLK;
        assert((int)ix->blk[sl] >= lo - 1 && (int)ix->blk[sl] <= hi + 1);
        seen++;
      }
      assert(seen > NS);
      free(rseq); }
    // the block ids survive a .pex round trip, as a PEX5 file
    { const char *tp = "gpeb_pex5_selftest.tmp";
      assert(pex_write(tp, ix, NULL) == 0);
      pidx_t *rd = pex_read(tp);
      assert(rd && rd->blk && rd->blk_bp == BLK && rd->blk_len == (uint64_t)LEN && !strcmp(rd->blk_ref, "REF"));
      size_t nchk = 0;
      for (size_t z = 0; z <= ix->km->mask; z++) {
        if (ix->km->st[z] != KH_IDX) continue;
        size_t t = kh_slot(rd->km, ix->km->key[z], 0);
        assert(t != (size_t)-1 && rd->km->st[t] == KH_IDX && rd->blk[t] == ix->blk[z]);
        nchk++;
      }
      assert(nchk == ix->nkmer);
      pidx_free(rd); remove(tp); }
    // the donor's two copies, and error-free reads from both
    char *cp[2];
    for (int c = 0; c < 2; c++) {
      cp[c] = malloc((size_t)LEN + 1); int o = 0;
      for (int i = 0; i <= NS; i++) {
        memcpy(cp[c] + o, sseq[i], SG); o += SG;
        if (i < NS) { int h = c ? 5 : (i < NS / 2 ? 1 : 3); cp[c][o++] = al[h][i] ? altal[i] : refal[i]; }
      }
      cp[c][o] = 0; assert(o == LEN);
    }
    const char *rp = "gpeb_mosaic_selftest.fq";
    FILE *rf = fopen(rp, "w");
    assert(rf);
    for (int c = 0; c < 2; c++)
      for (int st = 0; st + RL <= LEN; st += 7)
        fprintf(rf, "@c%d_%d\n%.*s\n+\n%.*s\n", c, st, RL, cp[c] + st, RL,
                "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII");
    fclose(rf);
    mosaic_opt_t mo; memset(&mo, 0, sizeof mo);
    mo.lambda = mo.rho = mo.eps = -1; mo.reads = MR_SHORT; mo.out = "gpeb_mosaic_selftest";
    mo.gfa_out = "gpeb_mosaic_selftest.out.gfa"; mo.keep = "REF"; mo.panel = "selftest";
    assert(run_mosaic(ix, rp, &mo, gfa, bg) == 0);
    FILE *pf = fopen("gpeb_mosaic_selftest.path.tsv", "r");
    assert(pf);
    char ln[512]; int nrow = 0, left_ok = 0, nleft = 0, right_ok = 0, nright = 0;
    while (fgets(ln, sizeof ln, pf)) {
      if (ln[0] == '#' || !strncmp(ln, "block", 5)) continue;
      char *f[8]; int nf = 0;
      for (char *t = strtok(ln, "\t\n"); t && nf < 8; t = strtok(NULL, "\t\n")) f[nf++] = t;
      assert(nf >= 6);
      int b = atoi(f[0]); nrow++;
      int has5 = !strcmp(f[3], "H5") || !strcmp(f[4], "H5");
      int has1 = !strcmp(f[3], "H1") || !strcmp(f[4], "H1");
      int has3 = !strcmp(f[3], "H3") || !strcmp(f[4], "H3");
      if (b <= 5)  { nleft++;  left_ok  += has5 && has1; }     // the switch is in block 7
      if (b >= 9)  { nright++; right_ok += has5 && has3; }
    }
    fclose(pf);
    assert(nrow == (LEN + BLK - 1) / BLK);
    assert(left_ok == nleft && right_ok == nright);
    // the splice: one edge-valid walk per track, as long as the reference
    FILE *of = fopen("gpeb_mosaic_selftest.out.gfa", "r");
    assert(of);
    char *wl = malloc(1 << 20); int nw1 = 0, nw2 = 0;
    while (fgets(wl, 1 << 20, of)) {
      if (strncmp(wl, "W\tmosaic\t", 9)) continue;
      char *t = strtok(wl, "\t"); (void)t;
      char *sm = strtok(NULL, "\t"), *hp = strtok(NULL, "\t"), *sq = strtok(NULL, "\t"), *st0 = strtok(NULL, "\t"),
           *en = strtok(NULL, "\t"); (void)sm; (void)sq; (void)st0;
      assert(atoi(en) == LEN);
      if (atoi(hp) == 1) nw1++; else nw2++;
    }
    fclose(of); free(wl);
    assert(nw1 == 1 && nw2 == 1);
    // splicing the SAVED painting gives the same graph
    assert(run_splice_file("gpeb_mosaic_selftest.path.tsv", gp, "gpeb_mosaic_selftest.out2.gfa", "REF") == 0);
    { FILE *a = fopen("gpeb_mosaic_selftest.out.gfa", "rb"), *b = fopen("gpeb_mosaic_selftest.out2.gfa", "rb");
      assert(a && b); int ca, cb, same = 1;
      do { ca = fgetc(a); cb = fgetc(b); if (ca != cb) same = 0; } while (ca != EOF && cb != EOF);
      assert(same); fclose(a); fclose(b); }
    printf("[ok] --set-call mosaic: painted %d blocks ({H1,H5} left of the recombination, {H3,H5} right of it),"
           " spliced into one edge-valid walk per track\n", nrow);

    // --set-call cover on the same reads: windows of 8 blocks (4 kb), so the recombination is
    // at the edge of the first one. The donor's alleles are covered by {H1, H5} on the left
    // and {H3, H5} on the right, and the stop rule must not ask for a third colour.
    { mosaic_opt_t co; memset(&co, 0, sizeof co);
      co.lambda = co.rho = co.eps = -1; co.reads = MR_SHORT; co.out = "gpeb_cover_selftest";
      co.gfa_out = "gpeb_cover_selftest.out.gfa"; co.keep = "REF"; co.panel = "selftest";
      co.cwin = 8 * BLK; co.cmax = 6; co.cfrac = 0.02; co.cgain = 3; co.cstray = -1;
      co.check = rp;                                          // in-sample here: the check must run
      assert(run_cover(ix, rp, &co, gfa, bg) == 0);
      FILE *cf = fopen("gpeb_cover_selftest.path.tsv", "r");
      assert(cf);
      char *hl = malloc(1 << 16); int ntr = 0, w0 = 0, w1 = 0, rows = 0;
      while (fgets(hl, 1 << 16, cf)) {
        if (hl[0] == '#') continue;
        char *fl[32]; int nf2 = 0;
        for (char *t = strtok(hl, "\t\n"); t && nf2 < 32; t = strtok(NULL, "\t\n")) fl[nf2++] = t;
        if (!strcmp(fl[0], "block")) { for (int z = 0; z < nf2; z++) ntr += !strncmp(fl[z], "track", 5); continue; }
        int b = atoi(fl[0]); rows++;
        int h1 = 0, h3 = 0, h5 = 0, n = 0;
        for (int z = 0; z < ntr; z++) { n += strcmp(fl[3 + z], "-") != 0;
          h1 += !strcmp(fl[3 + z], "H1"); h3 += !strcmp(fl[3 + z], "H3"); h5 += !strcmp(fl[3 + z], "H5"); }
        assert(n <= 3);
        if (b == 1) w0 = h1 && h5;                            // block 1 is in the first window
        if (b == 14) w1 = h3 && h5;                           // block 14 is in the second
      }
      fclose(cf); free(hl);
      assert(rows == (LEN + BLK - 1) / BLK && ntr >= 2 && w0 && w1);
      // the saved painting (with "-" for an empty track) splices to the same graph
      assert(run_splice_file("gpeb_cover_selftest.path.tsv", gp, "gpeb_cover_selftest.out2.gfa", "REF") == 0);
      FILE *a = fopen("gpeb_cover_selftest.out.gfa", "rb"), *b2 = fopen("gpeb_cover_selftest.out2.gfa", "rb");
      assert(a && b2); int ca, cb, same = 1;
      do { ca = fgetc(a); cb = fgetc(b2); if (ca != cb) same = 0; } while (ca != EOF && cb != EOF);
      assert(same); fclose(a); fclose(b2);
      printf("[ok] --set-call cover: %d tracks, {H1,H5} picked left of the recombination and {H3,H5} right of it,"
             " no third colour asked for; the saved painting splices to the same graph\n", ntr);
      // the parallel scan of the reads and the graph read AFTER the painting, for its colours only,
      // give the very same painting and graph as the serial scan and the whole-graph parse
      { mosaic_opt_t c2 = co; c2.out = "gpeb_cover_selftest2"; c2.gfa_out = "gpeb_cover_selftest2.out.gfa";
        c2.check = NULL; c2.threads = 4; c2.gfa_in = gp;
        assert(run_cover(ix, rp, &c2, NULL, NULL) == 0);
        const char *f1[2] = { "gpeb_cover_selftest.path.tsv", "gpeb_cover_selftest.out.gfa" };
        const char *f2[2] = { "gpeb_cover_selftest2.path.tsv", "gpeb_cover_selftest2.out.gfa" };
        for (int z = 0; z < 2; z++) {
          FILE *fa = fopen(f1[z], "rb"), *fb = fopen(f2[z], "rb");
          assert(fa && fb); int xa, xb, same2 = 1;
          do { xa = fgetc(fa); xb = fgetc(fb); if (xa != xb) same2 = 0; } while (xa != EOF && xb != EOF);
          assert(same2); fclose(fa); fclose(fb);
          if (z) remove(f2[z]); else { remove(f2[z]); remove("gpeb_cover_selftest2.sets.tsv"); }
        }
        // the reader keeps every W-line's header but parses the steps of the wanted samples only
        char *want[3] = { "H1", "H5", "REF" }; uint64_t nskip = 0;
        gfa_t *lg = gfa_read_walks(gp, want, 3, &nskip);
        assert(lg && nskip == (uint64_t)lg->n_walk - 3 && lg->n_walk == gfa->n_walk);
        uint32_t nfull = 0, nempty = 0;
        for (uint32_t w = 0; w < lg->n_walk; w++) { if (lg->walk[w].n_v > 0) nfull++; else nempty++; }
        assert(nfull == 3 && nempty == nskip);
        for (uint32_t w = 0; w < lg->n_walk; w++) if (lg->walk[w].n_v > 0)
          for (uint32_t w2 = 0; w2 < gfa->n_walk; w2++)
            if (!strcmp(lg->walk[w].sample, gfa->walk[w2].sample) && lg->walk[w].hap == gfa->walk[w2].hap) {
              assert(lg->walk[w].n_v == gfa->walk[w2].n_v);
              assert(!memcmp(lg->walk[w].v, gfa->walk[w2].v, (size_t)lg->walk[w].n_v * sizeof(uint32_t)));
            }
        gfa_destroy(lg);
        printf("[ok] cover: parallel scan and graph read after the painting give the same painting and graph;"
               " the reader skipped the steps of %llu of %u W-lines and kept the others exact\n",
               (unsigned long long)nskip, gfa->n_walk);
        // with a .gfai next to the GFA, the splice after the painting reads the painted colours chunk
        // by chunk, and writes the very same graph
        uint32_t ch0 = gfai_chunk; gfai_chunk = 2;
        assert(gfai_write(gp, gfa, bg, "REF") == 0);
        gfai_chunk = ch0;
        mosaic_opt_t c3 = co; c3.out = "gpeb_cover_selftest3"; c3.gfa_out = "gpeb_cover_selftest3.out.gfa";
        c3.check = NULL; c3.gfa_in = gp;
        assert(run_cover(ix, rp, &c3, NULL, NULL) == 0);
        { FILE *fa = fopen("gpeb_cover_selftest.out.gfa", "rb"), *fb = fopen("gpeb_cover_selftest3.out.gfa", "rb");
          assert(fa && fb); int xa, xb, same3 = 1;
          do { xa = fgetc(fa); xb = fgetc(fb); if (xa != xb) same3 = 0; } while (xa != EOF && xb != EOF);
          assert(same3); fclose(fa); fclose(fb); }
        { char gi[256]; snprintf(gi, sizeof gi, "%s.gfai", gp); remove(gi); }
        remove("gpeb_cover_selftest3.out.gfa"); remove("gpeb_cover_selftest3.path.tsv"); remove("gpeb_cover_selftest3.sets.tsv");
        printf("[ok] cover: with a .gfai the splice after the painting writes the same graph\n"); }
      remove("gpeb_cover_selftest.path.tsv"); remove("gpeb_cover_selftest.sets.tsv");
      remove("gpeb_cover_selftest.out.gfa"); remove("gpeb_cover_selftest.out2.gfa"); }

    remove(gp); remove(rp); remove("gpeb_mosaic_selftest.path.tsv");
    remove("gpeb_mosaic_selftest.out.gfa"); remove("gpeb_mosaic_selftest.out2.gfa");
    free(cp[0]); free(cp[1]);
    pidx_free(ix); gpeb_free(gs); bg_free(bg); gfa_destroy(gfa);
  }

  // --mosaic-temper: the median, chat and w by hand. Even and odd counts (blocks without
  // data ignored), w capped at 1, no data at all, the share block tempering scales down.
  { mobs_t T; memset(&T, 0, sizeof T);
    uint32_t even[7] = { 10, 0, 20, 30, 40, 50, 60 };          // median (30 + 40) / 2
    T.nB = 7; T.nobs_blk = even;
    temper_t s = mosaic_temper_stats(&T, 58.0);
    assert(fabs(s.median - 35.0) < 1e-12 && fabs(s.chat - 35.0 / 58.0) < 1e-12 && s.w == 1.0);
    assert(fabs(s.lt1 - 1.0 / 6.0) < 1e-12);                   // only 60 > 58
    uint32_t odd[7] = { 10, 20, 30, 40, 50, 60, 2000 };        // median 40: chat 4, w 1/4
    T.nobs_blk = odd;
    s = mosaic_temper_stats(&T, 10.0);
    assert(s.median == 40.0 && s.chat == 4.0 && s.w == 0.25 && fabs(s.lt1 - 6.0 / 7.0) < 1e-12);
    uint32_t none[3] = { 0, 0, 0 };
    T.nB = 3; T.nobs_blk = none;
    s = mosaic_temper_stats(&T, 58.0);
    assert(s.median == 0.0 && s.w == 1.0 && s.lt1 == 0.0);
    printf("[ok] --mosaic-temper: median, chat, w and the tempered share by hand\n");
  }

  printf("ok\n");
  return 0;
}
#else
// One line per flag; the rationale lives in the docs.
static void usage(const char *prog, FILE *out) {
  fprintf(out,
"usage: %s [options] <in.gfa> > out.fa\n"
"\n"
"Ranked guardian-pebbling of a bidirected GFA. One FASTA record per pebble:\n"
"  >>head>tail|hm|<sample>|<rank>|<pid>|<pos>, then the 2*span-bp junction window.\n"
"\n"
"cut mode, at most one (default: cut at any in-degree>=2 vertex; sound, not\n"
"always minimal):\n"
"  --minimal            exact, minimal guardian set on any graph\n"
"  --acyclic            exact on a DAG by one sweep per walk; refuses a cyclic graph\n"
"\n"
"cut options:\n"
"  --early-cut          cut at the divergence edge instead of the merge (any mode)\n"
"  --full-graph         reachability over all L-lines, not only walked edges\n"
"  --color-sep C        a color is the path name up to the first C ('.' or '#')\n"
"  --colors FILE        keep only these colors (one per line; a name without '#'\n"
"                       selects every haplotype of that sample)\n"
"  --min-gap N          skip pebbles closer than N bp along a walk (default %u)\n"
"  --min-complexity F   drop windows with distinct-4-mer fraction < F\n"
"  --span N             flank bp on each side of the junction (default %u)\n"
"  -c, --compress       one record per (edge, window), listing every color\n"
"  --from-peb F.peb     take the guardians from a .peb (plain/.gz), skip the cut\n"
"\n"
"index and classify:\n"
"  --index K OUT.pex    build the junction k-mer index, write it, stop\n"
"  --classify K READS   color FASTA/FASTQ reads (plain/.gz); per-read TSV on stdout\n"
"  --use-index F.pex    classify (or --subsample) against a saved index; pass the\n"
"                       GFA too to refuse a stale one\n"
"  --syncmer S          index only closed syncmers at s-mer length S (0 = off)\n"
"  --keep-ambiguous     also index k-mers found at several loci (all classes)\n"
"  --keep-many LIST     also index the KH_MANY k-mers of these locus classes, e.g.\n"
"                       segment,sibling,highcopy (strand, segment, sibling, distant,\n"
"                       highcopy; all = --keep-ambiguous); one more walk pass\n"
"  --fast-mask          color each k-mer from its pebble windows only, skipping the\n"
"                       walk stream and the locus check; unsound, since masks\n"
"                       can miss colors, so its queries default to --tol 0.9\n"
"  --tol F              keep a color in a read's set if >= F of its k-mer masks hold\n"
"                       it (default 1 = strict intersection; 0.9 for --fast-mask)\n"
"  --index-all          ABLATION: index every walk k-mer, not only junction ones\n"
"  --threads N          threads for the index build's streaming passes (default 1)\n"
"  --set-call MODE      whole-set call: all (default), fraction, intersect, vote,\n"
"                       sample, mosaic, cover\n"
"  --top N              ranks printed per set-call block (default 5, 0 = all)\n"
"  --subsample R OUT    with --use-index: keep a seeded, nested fraction R, stop\n"
"  --subsample-seed N   another --subsample draw\n"
"\n"
"mosaic: two haplotype tracks that recombine along the chromosome, fitted to a read set\n"
"  --mosaic-ref COLOUR  with --index: store a block id per record, projected on COLOUR's\n"
"                       walk(s) (PEX5, +2 bytes/record); needs a sound single-locus index.\n"
"                       Also writes GFA.gfai next to a plain GFA: where its W-lines' steps\n"
"                       sit, so a splice reads only the ones it uses (a .gfa.gz, or a GFA\n"
"                       without a valid .gfai, is read whole)\n"
"  --mosaic-block N     block size in bp (default 10000)\n"
"  --set-call mosaic    with --use-index (or --classify): paint the reads; the painting\n"
"                       (path.tsv) goes to stdout or --mosaic-out\n"
"  --mate FILE          second file of read pairs: one observation per pair\n"
"  --mosaic-reads T     hifi | short: the preset of lambda/rho/eps (default: by read length)\n"
"  --mosaic-lambda X, --mosaic-rho X, --mosaic-eps X   override the preset\n"
"  --mosaic-exclude L   colours or samples the tracks may not use (e.g. the references)\n"
"  --mosaic-out PFX     write PFX.path.tsv instead of stdout\n"
"  --mosaic-gfa OUT     splice the painting into a graph (needs the GFA as input)\n"
"  --mosaic-keep L      colours kept whole in OUT next to the two tracks (default: the\n"
"                       coordinate reference)\n"
"  --mosaic-fwd N       how far ahead of a switch to look for a shared node (default 200000)\n"
"  --mosaic-temper M    global | block | off: scale each block's evidence down to\n"
"                       --mosaic-temper-n observations, so more coverage does not buy\n"
"                       spurious switches (default: global for short reads, off for HiFi)\n"
"  --mosaic-temper-n N  observations per block the preset was fitted at\n"
"                       (default 58 per 10 kb of block)\n"
"  --skip-stale-check   with --use-index and the GFA: do not hash the GFA to check that the index\n"
"                       was built from it (2 s on a chr21 panel); the graph is then trusted\n"
"  --threads N          also the cover's scan of the reads (--set-call cover), same result\n"
"  --mosaic-splice P    splice a saved painting P (path.tsv) into the GFA, no reads;\n"
"                       with --mosaic-gfa OUT [--mosaic-keep L]\n"
"  --make-gfai COLOUR   write GFA.gfai alone, for an index built before it existed\n"
"                       (COLOUR = the index's --mosaic-ref)\n"
"  --no-gfai            splice reading the whole GFA even when it has a .gfai\n"
"  --gfai-chunk N       steps per chunk of a .gfai being written (default 1024)\n"
"\n"
"cover: per window of the chromosome, the few colours whose UNION carries the most of the\n"
"reads' index hits; N tracks instead of two, no model of how colours change. Same index,\n"
"--mate, --mosaic-ref/-out/-gfa/-keep/-exclude/-fwd/-reads as above (--mosaic-keep also\n"
"names the colours already in the graph, which cover what they carry before any pick)\n"
"  --set-call cover     with --use-index (or --classify); path.tsv has one column per track\n"
"                       and PFX.sets.tsv the picks and gains of every window\n"
"  --mosaic-check R     score the graph on an INDEPENDENT read set R (--mosaic-check-mate R2):\n"
"                       the share of its hits that no colour of the graph carries, and how many\n"
"                       hits the best single colour of the graph misses per read\n"
"  --mosaic-score P     with --use-index and --mosaic-check: score a saved painting P, no model run\n"
"  --cover-window BP    window size (default 200000; 0 = the whole chromosome; auto = EXPERIMENTAL,\n"
"                       windows of whole blocks that each hold --cover-hits hits)\n"
"  --cover-hits N       the hits an automatic window holds, at least (default 4000)\n"
"  --cover-max N        most colours per window (default 32)\n"
"  --cover-min-frac F   add a colour only while it covers >= F of the window's hits (default 0.0005)\n"
"  --cover-min-gain G   ...and at least G hits (default 2)\n"
"  --cover-stray R      drop hits farther than R blocks from their read's median (5 HiFi, 2 short)\n"
"  --cover-rule R       greedy (default) | support: CONTROL, the N colours with most hits each\n"
"  --cover-weight W     count (default, capped) | presence: a record weighs 1\n"
"\n"
"split the panel into one index per group, parsing the GFA once:\n"
"  --split-colors PFX   write PFX.g0.pex, PFX.g1.pex, ... and PFX.groups.tsv\n"
"  -k K                 k-mer length for --split-colors\n"
"  --split-max-colors N colors per group (default 64)\n"
"  --split-walk BP      walk bases per group (default 4e9, 0 = no cap)\n"
"  --split-shared LIST  comma-separated colors/samples copied into every group\n"
"  --split-exclude LIST comma-separated samples left out of the panel (leave-one-out)\n"
"  --split-dry-run      print the partition, build nothing\n"
"\n"
"diagnostics (report instead of FASTA):\n"
"  --kmer-stats K       how many windows carry a single-locus k-mer\n"
"  --pebble-ambiguity K candidate k-mers proposed by several windows, and whether\n"
"                       the window-derived mask matches the walk-derived one\n"
"  -h, --help           this help\n",
    prog, (unsigned)GSPAN, (unsigned)GSPAN);
}

#define MAX_SPAN    100000       // --span cap: a window is 2*span bp, held per pebble
#define MAX_THREADS 64           // --threads cap

// optval: the value after an option, refusing to swallow the NEXT option (so a forgotten
// value is not mistaken for a file name). A bare "-" is still a legal value.
static const char *optval(int argc, char **argv, int *i, const char *flag) {
  if (*i + 1 >= argc) {
    fprintf(stderr, "%s: missing value\n", flag);
    exit(1);
  }
  const char *v = argv[++*i];
  if (v[0] == '-' && v[1]) {
    fprintf(stderr, "%s: expected a value but got the option %s\n", flag, v);
    exit(1);
  }
  return v;
}

// optuint: an option's value as an integer in [lo, hi]; digits only, no wrap-around.
static unsigned long optuint(int argc, char **argv, int *i, const char *flag,
                             unsigned long lo, unsigned long hi) {
  const char *v = optval(argc, argv, i, flag);
  char *end; errno = 0;
  unsigned long x = strtoul(v, &end, 10);
  if (end == v || *end || errno || x < lo || x > hi) {
    fprintf(stderr, "%s: expected an integer in [%lu, %lu], got '%s'\n", flag, lo, hi, v);
    exit(1);
  }
  return x;
}


// Build one index per group. Returns 0 on success. The GFA is parsed by the caller and
// reused for every group.
static int split_build(gfa_t *gfa, const char *prefix, int walks, int full_graph,
                       int colsep, uint32_t maxcol, uint64_t maxwalk,
                       const char *shared_csv, const char *exclude_csv, const char *in,
                       int fast, uint32_t min_gap, uint32_t span,
                       double min_complexity, int early, int acyclic,
                       uint32_t ik, uint32_t keep_many, int fast_mask,
                       int all_kmers, int nthreads, uint32_t sk, int dry) {
  gfa_walk_t *src = walks ? gfa->walk : gfa->path;
  uint32_t ns = walks ? gfa->n_walk : gfa->n_path;
  if (!ns) { fprintf(stderr, "--split-colors: the GFA has no W or P lines\n"); return 1; }

  smp_t *sm = scan_samples(src, ns);
  wl_t *sh = shared_list(shared_csv);
  uint32_t nshared_col = 0;
  uint64_t shared_bp = 0;
  for (uint32_t i = 0; i < sm->n; i++)
    if (is_shared(sh, sm->name[i])) { nshared_col += sm->nhap[i]; shared_bp += sm->bp[i]; }
  if (sh->n && !nshared_col)
    fprintf(stderr, "--split-shared: no sample matched; groups will not share any colour\n");
  // --split-exclude: samples left out of the panel altogether, for leave-one-out.
  // A sample both shared and excluded is refused rather than silently resolved.
  wl_t *ex = shared_list(exclude_csv);
  uint32_t nex = 0;
  for (uint32_t i = 0; i < sm->n; i++)
    if (is_shared(ex, sm->name[i])) {
      if (is_shared(sh, sm->name[i])) {
        fprintf(stderr, "--split-exclude: %s is also in --split-shared\n", sm->name[i]);
        smp_free(sm); wl_free(sh); wl_free(ex); return 1;
      }
      nex++;
    }
  if (ex->n) fprintf(stderr, "--split-exclude: %u of %u samples left out\n", nex, sm->n);
  if (ex->n && !nex) fprintf(stderr, "--split-exclude: no sample matched\n");

  // Greedy left-to-right over the samples, in name order so the partition is
  // reproducible: close a group when adding the next sample would break either
  // cap. Sorting by size would balance better but would make the grouping
  // depend on the panel's composition, and every index already built would
  // stop corresponding to a group.
  uint32_t *order = malloc(sm->n * sizeof(uint32_t));
  uint32_t nord = 0;
  for (uint32_t i = 0; i < sm->n; i++)
    if (!is_shared(sh, sm->name[i]) && !is_shared(ex, sm->name[i])) order[nord++] = i;
  for (uint32_t a = 0; a + 1 < nord; a++)                    // insertion sort by name
    for (uint32_t b = a + 1; b < nord; b++)
      if (strcmp(sm->name[order[b]], sm->name[order[a]]) < 0) {
        uint32_t t = order[a]; order[a] = order[b]; order[b] = t;
      }
  if (maxcol <= nshared_col) {
    fprintf(stderr, "--split-max-colors %u leaves no room for the %u shared colours\n",
            maxcol, nshared_col);
    free(order); smp_free(sm); wl_free(sh); wl_free(ex); return 1;
  }

  char path[4096];
  snprintf(path, sizeof path, "%s.groups.tsv", prefix);
  FILE *mf = fopen(path, "w");
  if (mf) fprintf(mf, "group\tpex\tsamples\tcolors\twalk_bp\n");

  uint32_t gi = 0, i = 0;
  int rc = 0;
  while (i < nord && !rc) {
    // Accumulate one group.
    uint32_t ncol = nshared_col, nsam = 0;
    uint64_t bp = shared_bp;
    uint32_t start = i;
    while (i < nord) {
      uint32_t sidx = order[i];
      uint32_t c2 = ncol + sm->nhap[sidx];
      uint64_t b2 = bp + sm->bp[sidx];
      if (nsam && (c2 > maxcol || (maxwalk && b2 > maxwalk))) break;
      ncol = c2; bp = b2; nsam++; i++;
    }
    // Whitelist: the shared colours plus this group's samples, by name.
    wl_t *w = calloc(1, sizeof(wl_t));
    w->e = malloc((sh->n + nsam) * sizeof(char *));
    for (uint32_t z = 0; z < sh->n; z++) w->e[w->n++] = strdup(sh->e[z]);
    for (uint32_t z = start; z < i; z++) w->e[w->n++] = strdup(sm->name[order[z]]);

    snprintf(path, sizeof path, "%s.g%u.pex", prefix, gi);
    fprintf(stderr, "=== group g%u: %u samples, %u colours, %.2f Gbp -> %s\n",
            gi, nsam, ncol, (double)bp / 1e9, path);
    if (dry) {                      // --split-dry-run: report the partition only
      if (mf) fprintf(mf, "g%u\t%s\t%u\t%u\t%llu\n", gi, path, nsam, ncol,
                      (unsigned long long)bp);
      wl_free(w); gi++; continue;
    }
    bg_t *g = bg_build(gfa, walks, full_graph, colsep, w);
    if (!g) { wl_free(w); rc = 1; break; }
    gset_t *gs = gpeb(g, fast, min_gap, span, min_complexity, early, acyclic);
    if (!gs) { bg_free(g); wl_free(w); rc = 1; break; }
    pidx_t *ix = pidx_build(g, gs, span, ik, keep_many, fast_mask, all_kmers, nthreads, sk);
    if (!ix) rc = 1;
    else {
      if (pex_write(path, ix, in)) rc = 1;
      pidx_free(ix);
    }
    if (mf && !rc) fprintf(mf, "g%u\t%s\t%u\t%u\t%llu\n", gi, path, nsam, ncol,
                           (unsigned long long)bp);
    gpeb_free(gs); bg_free(g); wl_free(w);
    gi++;
  }
  if (mf) fclose(mf);
  if (!rc) fprintf(stderr, "--split-colors: %u groups %s, manifest %s.groups.tsv\n",
                   gi, dry ? "planned (dry run, nothing built)" : "written", prefix);
  free(order); smp_free(sm); wl_free(sh); wl_free(ex);
  return rc;
}

int main(int argc, char **argv) {
  const char *in = NULL; int compress = 0, early = 0, full_graph = 0;
  int minimal = 0;                 // default (neither): the conservative superset
  int colsep = 0;                            // --color-sep: cut path names here
  const char *colors = NULL;                                      // --colors
  const char *split_pfx = NULL, *split_shared = NULL;             // --split-colors
  const char *split_exclude = NULL;                               // --split-exclude
  uint32_t split_maxcol = 64;                                     // --split-max-colors
  int split_dry = 0;                                              // --split-dry-run
  uint64_t split_walk = 4000000000ULL;                            // --split-walk
  int acyclic = 0; uint32_t kstats = 0, pamb = 0;
  uint32_t ck = 0; const char *creads = NULL; int fast_mask = 0; uint32_t keep_many = 0;  // --classify
  int all_kmers = 0;                                              // --index-all
  uint32_t syncs = 0;                                             // --syncmer S
  const char *subout = NULL; double subrate = 0.0;                // --subsample
  uint64_t subseed = 0x5bf03635d9e2a1c7ULL;                       // --subsample-seed
  int setcall = SET_ALL;                                          // --set-call
  const char *mref = NULL; uint32_t mblock = 10000;               // --mosaic-ref/-block
  const char *msplice = NULL, *mscore = NULL;                     // --mosaic-splice, --mosaic-score
  const char *mkgfai = NULL;                                      // --make-gfai COLOUR
  mosaic_opt_t mo; memset(&mo, 0, sizeof mo);                     // --set-call mosaic
  mo.lambda = mo.rho = mo.eps = -1; mo.reads = MR_AUTO;
  mo.cwin = 200000; mo.cauto = 0; mo.chits = 4000; mo.cmax = 32; mo.cfrac = 0.0005; mo.cgain = 2; mo.cstray = -1;   // --set-call cover
  int cover_opt = 0;                                              // a --cover-* option was given
  int topn = 5;                                                   // --top
  double tol = -1;                          // --tol; <0 = per-index default
  const char *mkindex = NULL, *usepex = NULL; uint32_t ik = 0;    // --index/--use-index
  const char *frompeb = NULL;                                     // --from-peb
  int nthreads = 1;                                               // --threads
  int skip_stale = 0;                                             // --skip-stale-check
  uint32_t min_gap = GSPAN;                  // default: no overlapping span windows
  uint32_t span = GSPAN;                     // flank bp each side (--span)
  double min_complexity = 0.0;              // default: no low-complexity filter
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help")) {
      usage(argv[0], stdout); return 0;
    }
    else if (!strcmp(argv[i], "--minimal")) minimal = 1;
    else if (!strcmp(argv[i], "--full-graph")) full_graph = 1;
    else if (!strcmp(argv[i], "--color-sep"))
      colsep = optval(argc, argv, &i, "--color-sep")[0];  // truncate names here
    else if (!strcmp(argv[i], "--acyclic")) acyclic = 1;
    else if (!strcmp(argv[i], "-c") || !strcmp(argv[i], "--compress"))
      compress = 1;                           // one record per (edge, window)
    else if (!strcmp(argv[i], "--early-cut"))
      early = 1;                              // cut at the divergence edge, not the merge
    else if (!strcmp(argv[i], "--colors")) colors = optval(argc, argv, &i, "--colors");
    else if (!strcmp(argv[i], "--split-colors"))
      split_pfx = optval(argc, argv, &i, "--split-colors");
    else if (!strcmp(argv[i], "--split-max-colors"))
      split_maxcol = (uint32_t)strtoul(optval(argc, argv, &i, "--split-max-colors"), NULL, 10);
    else if (!strcmp(argv[i], "--split-walk"))
      split_walk = strtoull(optval(argc, argv, &i, "--split-walk"), NULL, 10);
    else if (!strcmp(argv[i], "--split-dry-run")) split_dry = 1;
    else if (!strcmp(argv[i], "--split-shared"))
      split_shared = optval(argc, argv, &i, "--split-shared");
    else if (!strcmp(argv[i], "--split-exclude"))
      split_exclude = optval(argc, argv, &i, "--split-exclude");
    else if (!strcmp(argv[i], "--min-gap"))
      min_gap = (uint32_t)strtoul(optval(argc, argv, &i, "--min-gap"), NULL, 10);
    else if (!strcmp(argv[i], "--min-complexity"))
      min_complexity = strtod(optval(argc, argv, &i, "--min-complexity"), NULL);
    else if (!strcmp(argv[i], "--kmer-stats"))
      kstats = (uint32_t)strtoul(optval(argc, argv, &i, "--kmer-stats"), NULL, 10);
    else if (!strcmp(argv[i], "--pebble-ambiguity"))
      pamb = (uint32_t)strtoul(optval(argc, argv, &i, "--pebble-ambiguity"), NULL, 10);
    else if (!strcmp(argv[i], "--span"))
      span = (uint32_t)optuint(argc, argv, &i, "--span", 1, MAX_SPAN);
    else if (!strcmp(argv[i], "--classify")) {
      ck = (uint32_t)strtoul(optval(argc, argv, &i, "--classify"), NULL, 10);
      creads = optval(argc, argv, &i, "--classify");
    }
    else if (!strcmp(argv[i], "--keep-ambiguous")) keep_many = KEEP_ALL;
    else if (!strcmp(argv[i], "--keep-many")) {
      if (parse_keep_many(optval(argc, argv, &i, "--keep-many"), &keep_many)) return 1;
    }
    else if (!strcmp(argv[i], "--fast-mask")) fast_mask = 1;
    else if (!strcmp(argv[i], "--index-all")) all_kmers = 1;
    else if (!strcmp(argv[i], "--syncmer"))
      syncs = (uint32_t)strtoul(optval(argc, argv, &i, "--syncmer"), NULL, 10);
    else if (!strcmp(argv[i], "--subsample")) {
      subrate = strtod(optval(argc, argv, &i, "--subsample"), NULL);
      subout = optval(argc, argv, &i, "--subsample");
    }
    else if (!strcmp(argv[i], "--tol")) {
      tol = strtod(optval(argc, argv, &i, "--tol"), NULL);
      if (!(tol > 0 && tol <= 1)) { fprintf(stderr, "--tol must be in (0, 1]\n"); return 1; }
    }
    else if (!strcmp(argv[i], "--top"))
      topn = (int)strtol(optval(argc, argv, &i, "--top"), NULL, 10);
    else if (!strcmp(argv[i], "--subsample-seed"))
      subseed = strtoull(optval(argc, argv, &i, "--subsample-seed"), NULL, 10);
    else if (!strcmp(argv[i], "--index")) {
      ik = (uint32_t)strtoul(optval(argc, argv, &i, "--index"), NULL, 10);
      mkindex = optval(argc, argv, &i, "--index");
    }
    else if (!strcmp(argv[i], "-k"))          // k alone: --split-colors has no OUT
      ik = (uint32_t)strtoul(optval(argc, argv, &i, "-k"), NULL, 10);
    else if (!strcmp(argv[i], "--use-index"))
      usepex = optval(argc, argv, &i, "--use-index");
    else if (!strcmp(argv[i], "--from-peb"))
      frompeb = optval(argc, argv, &i, "--from-peb");
    else if (!strcmp(argv[i], "--skip-stale-check")) skip_stale = 1;
    else if (!strcmp(argv[i], "--threads"))
      nthreads = (int)optuint(argc, argv, &i, "--threads", 1, MAX_THREADS);
    else if (!strcmp(argv[i], "--set-call")) {
      const char *m = optval(argc, argv, &i, "--set-call");
      if      (!strcmp(m, "all"))       setcall = SET_ALL;
      else if (!strcmp(m, "fraction"))  setcall = SET_FRACTION;
      else if (!strcmp(m, "intersect")) setcall = SET_INTERSECT;
      else if (!strcmp(m, "vote"))      setcall = SET_VOTE;
      else if (!strcmp(m, "sample"))    setcall = SET_SAMPLE;
      else if (!strcmp(m, "mosaic"))    setcall = SET_MOSAIC;
      else if (!strcmp(m, "cover"))     setcall = SET_COVER;
      else { fprintf(stderr, "--set-call: all|fraction|intersect|vote|sample|mosaic|cover\n"); return 1; }
    }
    else if (!strcmp(argv[i], "--mosaic-ref")) mref = optval(argc, argv, &i, "--mosaic-ref");
    else if (!strcmp(argv[i], "--mosaic-block"))
      mblock = (uint32_t)strtoul(optval(argc, argv, &i, "--mosaic-block"), NULL, 10);
    else if (!strcmp(argv[i], "--mate")) mo.mate = optval(argc, argv, &i, "--mate");
    else if (!strcmp(argv[i], "--mosaic-reads")) {
      const char *m = optval(argc, argv, &i, "--mosaic-reads");
      if      (!strcmp(m, "hifi"))  mo.reads = MR_HIFI;
      else if (!strcmp(m, "short")) mo.reads = MR_SHORT;
      else { fprintf(stderr, "--mosaic-reads: hifi|short\n"); return 1; }
    }
    else if (!strcmp(argv[i], "--mosaic-lambda"))
      mo.lambda = strtod(optval(argc, argv, &i, "--mosaic-lambda"), NULL);
    else if (!strcmp(argv[i], "--mosaic-rho"))
      mo.rho = strtod(optval(argc, argv, &i, "--mosaic-rho"), NULL);
    else if (!strcmp(argv[i], "--mosaic-eps"))
      mo.eps = strtod(optval(argc, argv, &i, "--mosaic-eps"), NULL);
    else if (!strcmp(argv[i], "--mosaic-exclude")) mo.exclude = optval(argc, argv, &i, "--mosaic-exclude");
    else if (!strcmp(argv[i], "--mosaic-out")) mo.out = optval(argc, argv, &i, "--mosaic-out");
    else if (!strcmp(argv[i], "--mosaic-gfa")) mo.gfa_out = optval(argc, argv, &i, "--mosaic-gfa");
    else if (!strcmp(argv[i], "--mosaic-splice")) msplice = optval(argc, argv, &i, "--mosaic-splice");
    else if (!strcmp(argv[i], "--no-gfai")) gfai_off = 1;
    else if (!strcmp(argv[i], "--make-gfai")) mkgfai = optval(argc, argv, &i, "--make-gfai");
    else if (!strcmp(argv[i], "--gfai-chunk")) {
      long v = strtol(optval(argc, argv, &i, "--gfai-chunk"), NULL, 10);
      if (v < 1 || v > 1 << 20) { fprintf(stderr, "--gfai-chunk: 1 to 1048576 steps\n"); return 1; }
      gfai_chunk = (uint32_t)v; }
    else if (!strcmp(argv[i], "--mosaic-score")) mscore = optval(argc, argv, &i, "--mosaic-score");
    else if (!strcmp(argv[i], "--mosaic-check")) mo.check = optval(argc, argv, &i, "--mosaic-check");
    else if (!strcmp(argv[i], "--mosaic-check-mate")) mo.check_mate = optval(argc, argv, &i, "--mosaic-check-mate");
    else if (!strcmp(argv[i], "--cover-window")) {
      const char *w = optval(argc, argv, &i, "--cover-window"); cover_opt = 1;
      if (!strcmp(w, "auto")) mo.cauto = 1;
      else { mo.cwin = (uint32_t)strtoul(w, NULL, 10); mo.cauto = 0; } }
    else if (!strcmp(argv[i], "--cover-hits")) {
      mo.chits = strtoull(optval(argc, argv, &i, "--cover-hits"), NULL, 10); cover_opt = 1; }
    else if (!strcmp(argv[i], "--cover-max")) {
      mo.cmax = (uint32_t)strtoul(optval(argc, argv, &i, "--cover-max"), NULL, 10); cover_opt = 1; }
    else if (!strcmp(argv[i], "--cover-min-frac")) {
      mo.cfrac = strtod(optval(argc, argv, &i, "--cover-min-frac"), NULL); cover_opt = 1; }
    else if (!strcmp(argv[i], "--cover-min-gain")) {
      mo.cgain = strtod(optval(argc, argv, &i, "--cover-min-gain"), NULL); cover_opt = 1; }
    else if (!strcmp(argv[i], "--cover-stray")) {
      mo.cstray = (int)strtol(optval(argc, argv, &i, "--cover-stray"), NULL, 10); cover_opt = 1; }
    else if (!strcmp(argv[i], "--cover-rule")) {
      const char *m = optval(argc, argv, &i, "--cover-rule"); cover_opt = 1;
      if      (!strcmp(m, "greedy"))  mo.crule = CR_GREEDY;
      else if (!strcmp(m, "support")) mo.crule = CR_SUPPORT;
      else { fprintf(stderr, "--cover-rule: greedy|support\n"); return 1; }
    }
    else if (!strcmp(argv[i], "--cover-weight")) {
      const char *m = optval(argc, argv, &i, "--cover-weight"); cover_opt = 1;
      if      (!strcmp(m, "count"))    mo.cweight = CW_COUNT;
      else if (!strcmp(m, "presence")) mo.cweight = CW_PRESENCE;
      else { fprintf(stderr, "--cover-weight: count|presence\n"); return 1; }
    }
    else if (!strcmp(argv[i], "--mosaic-fwd"))
      splice_fwd = (int64_t)strtoll(optval(argc, argv, &i, "--mosaic-fwd"), NULL, 10);
    else if (!strcmp(argv[i], "--mosaic-keep")) mo.keep = optval(argc, argv, &i, "--mosaic-keep");
    else if (!strcmp(argv[i], "--mosaic-temper")) {
      const char *m = optval(argc, argv, &i, "--mosaic-temper");
      if      (!strcmp(m, "off"))    mo.temper = TEMPER_OFF;
      else if (!strcmp(m, "global")) mo.temper = TEMPER_GLOBAL;
      else if (!strcmp(m, "block"))  mo.temper = TEMPER_BLOCK;
      else { fprintf(stderr, "--mosaic-temper: global, block or off\n"); return 1; }
    }
    else if (!strcmp(argv[i], "--mosaic-temper-n"))
      mo.temper_n = strtod(optval(argc, argv, &i, "--mosaic-temper-n"), NULL);
    // An unrecognised option is an error, not the input filename.
    else if (argv[i][0] == '-' && argv[i][1]) {
      fprintf(stderr, "unknown option: %s (try --help)\n", argv[i]);
      return 1;
    }
    else in = argv[i];
  }
  // Cut modes are one choice, not flags that stack: the default, --minimal or --acyclic.
  if (minimal && acyclic) {
    fprintf(stderr, "--minimal and --acyclic are alternative cut modes; pass one\n");
    return 1;
  }
  int fast = !minimal && !acyclic;
  if (mref && (keep_many || fast_mask || all_kmers || split_pfx)) {
    fprintf(stderr, "--mosaic-ref needs a sound single-locus, single-index build: not with"
                    " --keep-many, --keep-ambiguous, --fast-mask, --index-all or --split-colors\n");
    return 1;
  }
  if (mref && usepex) {
    fprintf(stderr, "--mosaic-ref is a build option: pass it with --index, to a GFA\n"); return 1;
  }
  if (msplice) {                              // splice a saved painting, nothing else
    if (!in || !mo.gfa_out) {
      fprintf(stderr, "--mosaic-splice PATH.tsv needs --mosaic-gfa OUT and the GFA as input\n"); return 1;
    }
    return run_splice_file(msplice, in, mo.gfa_out, mo.keep);
  }
  if (mkgfai) {                               // the .gfai of a GFA whose index is already built
    if (!in || usepex || mkindex) { fprintf(stderr, "--make-gfai COLOUR takes the GFA alone as input\n"); return 1; }
    uint64_t nsk = 0;
    gfa_t *gf = gfa_read_walks(in, NULL, 0, &nsk);       // every W-line, with gpeb's fast step parser
    if (!gf) gf = gfa_read(in);
    if (!gf) { fprintf(stderr, "failed to read GFA: %s\n", in); return 1; }
    bg_quiet = 1; bg_nopadj = 1;
    bg_t *bg = bg_build(gf, 1, 0, 0, NULL);
    bg_quiet = 0; bg_nopadj = 0;
    int rcg = !bg || gfai_write(in, gf, bg, mkgfai);
    bg_free(bg); gfa_destroy(gf);
    return rcg;
  }
  if ((setcall == SET_MOSAIC || setcall == SET_COVER) && !creads) {
    fprintf(stderr, "--set-call mosaic and cover need --classify K READS\n"); return 1;
  }
  if ((mo.mate || mo.gfa_out || mo.out || mo.exclude || mo.keep || mo.reads || mo.lambda >= 0 ||
       mo.rho >= 0 || mo.eps >= 0 || mo.check || mo.check_mate) && setcall != SET_MOSAIC && setcall != SET_COVER && !mscore) {
    fprintf(stderr, "--mate and the --mosaic-* options belong to --set-call mosaic and cover\n"); return 1;
  }
  if ((mo.lambda >= 0 || mo.rho >= 0 || mo.eps >= 0) && setcall == SET_COVER) {
    fprintf(stderr, "--mosaic-lambda, -rho and -eps belong to --set-call mosaic; cover has --cover-*\n"); return 1;
  }
  if (cover_opt && setcall != SET_COVER) {
    fprintf(stderr, "the --cover-* options belong to --set-call cover\n"); return 1;
  }
  if (mo.gfa_out && !(in && (usepex || mref))) {
    fprintf(stderr, "--mosaic-gfa needs the GFA the index was built from, as input\n"); return 1;
  }
  // Classifying from a prebuilt .pex skips the graph entirely -- that is the point
  // of the cache, so no GFA is required (or read) on this path.
  if (usepex) {
    if (!creads && !subout && !mscore) {
      fprintf(stderr, "--use-index needs --classify, --subsample or --mosaic-score\n"); return 1;
    }
    if (subout && !(subrate > 0.0)) {
      fprintf(stderr, "--subsample: rate must be > 0 (got %g)\n", subrate); return 1;
    }
    pidx_t *ix = pex_read(usepex);
    if (!ix) return 1;
    if (ck && ck != ix->k)
      fprintf(stderr, "note: --classify k=%u ignored, the index was built with k=%u\n",
              ck, ix->k);
    if (syncs && syncs != ix->sk)
      fprintf(stderr, "note: --syncmer %u ignored, the index was built with s=%u"
                      " (the filter is baked into which k-mers it holds)\n",
              syncs, ix->sk);
    // Staleness. Checkable only against the graph, so passing it is what enables
    // the check -- and omitting it is the deliberate way to skip one. A mismatch
    // is fatal: a stale index answers with the colors of a graph that no longer
    // exists, which is worse than not answering.
    if (in && skip_stale) {
      fprintf(stderr, "pex: staleness check skipped (--skip-stale-check): %s is trusted to be the graph the index was built from\n", in);
    } else if (in) {
      uint64_t sz = 0, h = file_hash(in, &sz);
      if (h != ix->src_hash || sz != ix->src_size) {
        fprintf(stderr, "pex: %s is STALE -- built from %s\n", usepex,
                ix->src_path ? ix->src_path : "?");
        if (sz != ix->src_size)
          fprintf(stderr, "     that was %llu bytes; %s is %llu bytes now\n",
                  (unsigned long long)ix->src_size, in, (unsigned long long)sz);
        else
          fprintf(stderr, "     same size (%llu bytes) but the contents differ\n",
                  (unsigned long long)sz);
        fprintf(stderr, "     rebuild:  gpeb --index %u <out.pex> %s\n", ix->k, in);
        pidx_free(ix); return 1;
      }
      fprintf(stderr, "pex: source matches %s\n", in);
      mo.gfa_hash = h; mo.gfa_hashed = 1;
    } else {
      fprintf(stderr, "pex: staleness not checked (pass the GFA to verify it)\n");
    }
    // After the staleness check on purpose: thinning a stale index would just
    // produce a second stale index, with the mismatch one step further from view.
    if (subout) {
      pex_subsample(ix, subrate, subseed);
      int rcs = pex_write(subout, ix, NULL) ? 1 : 0;   // NULL: inherit the source
      pidx_free(ix);
      return rcs;
    }
    if (mscore) {
      if (!mo.check) { fprintf(stderr, "--mosaic-score needs --mosaic-check READS\n"); pidx_free(ix); return 1; }
      int rcs = run_score_file(ix, mscore, &mo);
      pidx_free(ix);
      return rcs;
    }
    if (setcall == SET_MOSAIC || setcall == SET_COVER) {
      const char *bn = strrchr(usepex, '/');
      mo.panel = bn ? bn + 1 : usepex;
      mo.threads = nthreads;
      int rcm = 0;
      mo.gfa_in = in;                         // the splice parses the graph AFTER the painting, for the painting's colours only
      rcm = setcall == SET_COVER ? run_cover(ix, creads, &mo, NULL, NULL) : run_mosaic(ix, creads, &mo, NULL, NULL);
      pidx_free(ix);
      return rcm;
    }
    run_classify(ix, creads, setcall, topn, tol, stdout);
    pidx_free(ix);
    return 0;
  }
  if (!in) { usage(argv[0], stderr); return 1; }
  if (fast_mask && all_kmers) {
    fprintf(stderr, "--fast-mask colors k-mers from their windows; --index-all has none\n");
    return 1;
  }
  double t_p = now_s();
  gfa_t *gfa = gfa_read(in);
  double t_parse = now_s() - t_p; t_p = now_s();
  if (!gfa) { fprintf(stderr, "failed to read GFA: %s\n", in); return 1; }
  int walks = gfa->n_walk > 0;              // auto: W-lines if present, else P-lines
  if (split_pfx) {                          // one index per group, GFA parsed once
    if (split_maxcol < 1) { fprintf(stderr, "--split-max-colors must be >= 1\n");
                            gfa_destroy(gfa); return 1; }
    if (!ik && !split_dry) {
      fprintf(stderr, "--split-colors needs the k-mer length: -k 31\n");
      gfa_destroy(gfa); return 1;
    }
    int rc2 = split_build(gfa, split_pfx, walks, full_graph, colsep, split_maxcol,
                          split_walk, split_shared, split_exclude, in, fast, min_gap, span,
                          min_complexity, early, acyclic,
                          ik, keep_many, fast_mask, all_kmers, nthreads, syncs, split_dry);
    gfa_destroy(gfa);
    return rc2;
  }
  wl_t *wl = colors ? wl_read(colors) : NULL;
  if (colors && !wl) { gfa_destroy(gfa); return 1; }
  bg_t *g = bg_build(gfa, walks, full_graph, colsep, wl);
  wl_free(wl);
  if (!g) { gfa_destroy(gfa); return 1; }
  double t_bg = now_s() - t_p; t_p = now_s();
  gset_t *gs = frompeb ? peb_read(frompeb, g)
                       : gpeb(g, fast, min_gap, span, min_complexity, early,
                              acyclic);
  if (!gs) { bg_free(g); gfa_destroy(gfa); return 1; }
  // These three are SERIAL -- --threads does not touch them -- so an index build's
  // remaining headroom is here, not in the streaming passes.
  if (mkindex || creads)
    fprintf(stderr, "serial: gfa_read %.2fs  bg_build %.2fs  cut %.2fs\n",
            t_parse, t_bg, now_s() - t_p);
  int rc = 0;
  if (mkindex) {                            // --index: build and cache, then stop
    pidx_t *ix = pidx_build(g, gs, span, ik, keep_many, fast_mask, all_kmers, nthreads, syncs);
    if (!ix) rc = 1;
    else {
      if (mref && pidx_set_blocks(ix, gfa, g, mref, mblock)) rc = 1;
      else rc = pex_write(mkindex, ix, in) ? 1 : 0;
      // the .gfai a splice on this index's positions reads the GFA through (a warning if it cannot be written)
      if (!rc && mref && !gfai_off && walks && !colsep) gfai_write(in, gfa, g, mref);
      pidx_free(ix);
    }
  }
  else if (creads) {                        // --classify straight from the GFA
    pidx_t *ix = pidx_build(g, gs, span, ck, keep_many, fast_mask, all_kmers, nthreads, syncs);
    if (!ix) rc = 1;
    else {
      if (mref && pidx_set_blocks(ix, gfa, g, mref, mblock)) rc = 1;
      else if (setcall == SET_MOSAIC || setcall == SET_COVER) {
        mo.panel = in;
        rc = setcall == SET_COVER ? run_cover(ix, creads, &mo, gfa, g) : run_mosaic(ix, creads, &mo, gfa, g);
      }
      else run_classify(ix, creads, setcall, topn, tol, stdout);
      pidx_free(ix);
    }
  }
  else if (kstats)   kstat(g, gs, span, kstats, syncs, stdout);  // analysis, not a dump
  else if (pamb) {                          // --pebble-ambiguity: analysis, not a dump
    // Never --keep-ambiguous here: it would promote KH_MANY slots to KH_IDX and
    // erase the very distinction this report is about.
    pidx_t *ix = pidx_build(g, gs, span, pamb, 0, 0, all_kmers, nthreads, syncs);
    if (!ix) rc = 1; else { pebble_ambiguity(g, gs, ix, !all_kmers, nthreads, stderr); pidx_free(ix); }
  }
  else if (compress) gdump_c(g, gs, stdout, span);
  else               gdump(g, gs, stdout, span);
  gpeb_free(gs); bg_free(g); gfa_destroy(gfa);
  return rc;
}
#endif
