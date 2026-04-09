#!/bin/bash
# Measure peak memory usage of CONGA
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
CONGA_DIR="$(dirname "$SCRIPT_DIR")"

# Run CONGA and measure peak RSS using /usr/bin/time on macOS
/usr/bin/time -l "$CONGA_DIR/conga" \
    --input "$SCRIPT_DIR/test.bam" \
    --ref "$SCRIPT_DIR/ref.fa" \
    --sonic "$SCRIPT_DIR/test.sonic" \
    --dels "$SCRIPT_DIR/known_dels.bed" \
    --dups "$SCRIPT_DIR/known_dups.bed" \
    --mappability "$SCRIPT_DIR/mappability.bed" \
    --exclude "$SCRIPT_DIR/low_map_exclude.bed" \
    --out "$SCRIPT_DIR/output/mem_test" \
    --first-chr 0 --last-chr 4 2>&1 | grep -E "(maximum resident|peak memory)"
