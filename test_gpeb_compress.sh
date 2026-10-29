#!/usr/bin/env sh
# test_gpeb_compress.sh — proves gpeb -c MUST group by (edge, sequence), not by
# edge alone.
#
# Setup: head node 3 is 3bp (< span), so the span window's left flank spills into
# the PREDECESSOR, which differs by path (node 1 = A.. for walk X, node 2 = C..
# for walk Y). Both walks cross the same edge 3+>4+ (4 is a merge, so it gets
# pebbled), but the two windows differ. An edge-only key would collapse them into
# one record and silently drop one window; edge+seq must keep them separate.
set -eu
cd "$(dirname "$0")"
make gpeb >/dev/null
[ -x ./gpeb ] || { echo "FAIL: gpeb did not build"; exit 1; }

rep() { i=0; while [ "$i" -lt "$1" ]; do printf %s "$2"; i=$((i + 1)); done; }
gfa=$(mktemp); trap 'rm -f "$gfa"' EXIT
{
  printf 'H\tVN:Z:1.1\n'
  printf 'S\t1\t%s\n' "$(rep 60 A)"      # X's predecessor
  printf 'S\t2\t%s\n' "$(rep 60 C)"      # Y's predecessor
  printf 'S\t3\tGGG\n'                    # head: 3bp, < span
  printf 'S\t4\t%s\n' "$(rep 60 T)"      # tail
  printf 'S\t5\t%s\n' "$(rep 60 A)"      # 2nd source into 4 -> 4 is a merge
  for e in '1 3' '2 3' '3 4' '5 4'; do
    set -- $e; printf 'L\t%s\t+\t%s\t+\t0M\n' "$1" "$2"
  done
  printf 'W\tX\t0\tc\t0\t123\t>1>3>4\n'
  printf 'W\tY\t0\tc\t0\t123\t>2>3>4\n'
} > "$gfa"

# --full-graph: node 5 is the 2nd entrance to node 4 but is in no walk, so the
# merge at 4 only exists in the full L-line graph, not in G_P (the default).
out=$(./gpeb --full-graph --min-gap 0 -c "$gfa" 2>/dev/null)

# gpeb (edge+seq): edge 3+>4+ must be TWO records ...
n=$(printf '%s\n' "$out" | grep -c '^>>3+>4+')
# ... with TWO distinct window sequences.
seqs=$(printf '%s\n' "$out" | grep -A1 '^>>3+>4+' | grep -vE '^(>>|--)' | sort -u | wc -l)
# edge-ONLY key (the wrong grouping) would see just ONE >>3+>4+ edge -> collapse
# both into one record, losing a window.
edge=$(printf '%s\n' "$out" | grep '^>>3+>4+' | sed 's/|.*//' | sort -u | wc -l)

fail=0
[ "$n" -eq 2 ]    || { echo "FAIL: expected 2 records for 3+>4+, got $n"; fail=1; }
[ "$seqs" -eq 2 ] || { echo "FAIL: expected 2 distinct windows, got $seqs"; fail=1; }
[ "$edge" -eq 1 ] || { echo "FAIL: expected 1 distinct edge string, got $edge"; fail=1; }

if [ "$fail" -eq 0 ]; then
  echo "PASS: edge 3+>4+ kept as 2 records / 2 windows; edge-only key (1 edge) would wrongly merge."
fi
exit "$fail"
