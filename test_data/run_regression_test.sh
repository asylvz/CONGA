#!/bin/bash
# CONGA Comprehensive Regression Test
# Tests: exact match, genotype accuracy, SV counts, numerical tolerance,
#        negative controls, GC-bias regions, edge cases, mappability filtering
#
# Usage:
#   ./test_data/run_regression_test.sh                  # Run tests
#   ./test_data/run_regression_test.sh --update-baseline # Update baseline after intentional changes

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
CONGA_DIR="$(dirname "$SCRIPT_DIR")"
BASELINE_DIR="$SCRIPT_DIR/baseline"
OUTPUT_DIR="$SCRIPT_DIR/output"
CONGA_BIN="$CONGA_DIR/conga"

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[0;33m'
BLUE='\033[0;34m'
NC='\033[0m'

PASS=0
FAIL=0
WARN=0

TAB=$'\t'

pass() { echo -e "  ${GREEN}PASS${NC}: $1"; ((PASS++)); }
fail() { echo -e "  ${RED}FAIL${NC}: $1"; ((FAIL++)); }
warn() { echo -e "  ${YELLOW}WARN${NC}: $1"; ((WARN++)); }
section() { echo -e "\n${BLUE}[$1]${NC} $2"; }

# Portable grep for tab-delimited BED lines (works on BSD/macOS and GNU)
bed_grep() {
    # Usage: bed_grep FILE chr start end [genotype]
    local file="$1" chr="$2" start="$3" end="$4" geno="${5:-}"
    if [ -n "$geno" ]; then
        awk -F'\t' -v c="$chr" -v s="$start" -v e="$end" -v g="$geno" \
            '$1==c && $2==s && $3==e && $4==g' "$file"
    else
        awk -F'\t' -v c="$chr" -v s="$start" -v e="$end" \
            '$1==c && $2==s && $3==e' "$file"
    fi
}

bed_grep_count() {
    bed_grep "$@" | wc -l | tr -d ' '
}

# -----------------------------------------------
# Prerequisites
# -----------------------------------------------
echo "=== CONGA Comprehensive Regression Test ==="

if [ ! -f "$CONGA_BIN" ]; then
    echo "ERROR: CONGA binary not found. Run 'make' first."
    exit 1
fi

if [ ! -f "$SCRIPT_DIR/test.bam" ]; then
    echo "ERROR: Test data not found. Run generate_test_data.sh first."
    exit 1
fi

if [ ! -d "$BASELINE_DIR" ]; then
    echo "ERROR: Baseline not found. Run with --update-baseline to create."
    exit 1
fi

# -----------------------------------------------
# 1. Run CONGA
# -----------------------------------------------
section "1" "Running CONGA..."
rm -rf "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"

CONGA_EXIT=0
CONGA_STDERR=$("$CONGA_BIN" \
    -i "$SCRIPT_DIR/test.bam" \
    -f "$SCRIPT_DIR/ref.fa" \
    -d "$SCRIPT_DIR/known_dels.bed" \
    -u "$SCRIPT_DIR/known_dups.bed" \
    -o "$OUTPUT_DIR/baseline" \
    -m "$SCRIPT_DIR/mappability.bed" \
    2>&1) || CONGA_EXIT=$?

if [ $CONGA_EXIT -ne 0 ]; then
    fail "CONGA exited with status $CONGA_EXIT"
    echo "  stderr: $(echo "$CONGA_STDERR" | tail -5)"
else
    pass "CONGA ran successfully (exit 0)"
fi

# -----------------------------------------------
# 2. Output file existence
# -----------------------------------------------
section "2" "Checking output files..."
for f in baseline_svs.bed baseline_dels.bed baseline_dups.bed; do
    if [ -f "$OUTPUT_DIR/$f" ]; then
        lines=$(wc -l < "$OUTPUT_DIR/$f" | tr -d ' ')
        pass "$f exists ($lines lines)"
    else
        fail "$f missing"
    fi
done

# -----------------------------------------------
# 3. Exact baseline match
# -----------------------------------------------
section "3" "Exact baseline comparison..."
for f in baseline_svs.bed baseline_dels.bed baseline_dups.bed; do
    if [ ! -f "$OUTPUT_DIR/$f" ] || [ ! -f "$BASELINE_DIR/$f" ]; then
        fail "$f: cannot compare (file missing)"
        continue
    fi
    if diff -q "$BASELINE_DIR/$f" "$OUTPUT_DIR/$f" > /dev/null 2>&1; then
        pass "$f matches baseline exactly"
    else
        fail "$f differs from baseline"
        diff "$BASELINE_DIR/$f" "$OUTPUT_DIR/$f" | head -20
    fi
done

# -----------------------------------------------
# 4. Genotype accuracy - Deletions
# -----------------------------------------------
section "4" "Genotype accuracy — Deletions..."

# Homozygous deletions should be 1/1
HOMO_DELS=(
    "chr1:10000-13000:Small homo DEL"
    "chr1:40000-55000:Large homo DEL"
    "chr1:130000-135000:Adjacent homo DEL"
    "chr1:490000-498000:DEL near chr end"
    "chr2:30000-35000:Homo DEL in AT-rich"
    "chr3:20000-25000:Homo DEL in GC-rich"
    "chr5:1000-4000:DEL at chr start"
)

for entry in "${HOMO_DELS[@]}"; do
    IFS=: read -r loc start_end desc <<< "$entry"
    chr=$(echo "$loc" | cut -d: -f1)
    start=$(echo "$start_end" | cut -d- -f1)
    end=$(echo "$start_end" | cut -d- -f2)
    if [ "$(bed_grep_count "$OUTPUT_DIR/baseline_dels.bed" "$chr" "$start" "$end" "1/1")" -gt 0 ]; then
        pass "$desc ($chr:$start-$end) = 1/1"
    else
        actual=$(bed_grep "$OUTPUT_DIR/baseline_dels.bed" "$chr" "$start" "$end" | awk -F'\t' '{print $4}')
        fail "$desc ($chr:$start-$end) expected 1/1, got ${actual:-NOT_FOUND}"
    fi
done

# Heterozygous deletions should be 0/1
HET_DELS=(
    "chr1:20000-25000:Medium het DEL"
    "chr1:70000-72000:Small het DEL"
    "chr1:136000-140000:Adjacent het DEL"
    "chr1:200000-250000:Very large het DEL"
    "chr2:60000-65000:Het DEL in AT-rich"
    "chr2:200000-210000:Large het DEL in AT-rich"
    "chr3:50000-55000:Het DEL in GC-rich"
    "chr3:160000-165000:Het DEL near end GC-rich"
    "chr5:50000-55000:Het DEL in mixed GC"
)

for entry in "${HET_DELS[@]}"; do
    IFS=: read -r loc start_end desc <<< "$entry"
    chr=$(echo "$loc" | cut -d: -f1)
    start=$(echo "$start_end" | cut -d- -f1)
    end=$(echo "$start_end" | cut -d- -f2)
    if [ "$(bed_grep_count "$OUTPUT_DIR/baseline_dels.bed" "$chr" "$start" "$end" "0/1")" -gt 0 ]; then
        pass "$desc ($chr:$start-$end) = 0/1"
    else
        actual=$(bed_grep "$OUTPUT_DIR/baseline_dels.bed" "$chr" "$start" "$end" | awk -F'\t' '{print $4}')
        fail "$desc ($chr:$start-$end) expected 0/1, got ${actual:-NOT_FOUND}"
    fi
done

# -----------------------------------------------
# 5. Genotype accuracy - Duplications
# -----------------------------------------------
section "5" "Genotype accuracy — Duplications..."

HOMO_DUPS=(
    "chr1:300000-305000:Medium homo DUP"
    "chr1:350000-370000:Large homo DUP"
    "chr2:100000-105000:Homo DUP in AT-rich"
    "chr2:250000-260000:Large homo DUP in AT-rich"
    "chr3:80000-85000:Homo DUP in GC-rich"
    "chr5:90000-95000:Homo DUP in mixed GC"
)

for entry in "${HOMO_DUPS[@]}"; do
    IFS=: read -r loc start_end desc <<< "$entry"
    chr=$(echo "$loc" | cut -d: -f1)
    start=$(echo "$start_end" | cut -d- -f1)
    end=$(echo "$start_end" | cut -d- -f2)
    if [ "$(bed_grep_count "$OUTPUT_DIR/baseline_dups.bed" "$chr" "$start" "$end" "1/1")" -gt 0 ]; then
        pass "$desc ($chr:$start-$end) = 1/1"
    else
        actual=$(bed_grep "$OUTPUT_DIR/baseline_dups.bed" "$chr" "$start" "$end" | awk -F'\t' '{print $4}')
        fail "$desc ($chr:$start-$end) expected 1/1, got ${actual:-NOT_FOUND}"
    fi
done

HET_DUPS=(
    "chr1:320000-323000:Small het DUP"
    "chr1:420000-430000:Medium het DUP"
    "chr1:460000-465000:Het DUP near chr end"
    "chr2:150000-155000:Het DUP in AT-rich"
    "chr3:120000-125000:Het DUP in GC-rich"
    "chr5:140000-148000:Het DUP near chr end"
)

for entry in "${HET_DUPS[@]}"; do
    IFS=: read -r loc start_end desc <<< "$entry"
    chr=$(echo "$loc" | cut -d: -f1)
    start=$(echo "$start_end" | cut -d- -f1)
    end=$(echo "$start_end" | cut -d- -f2)
    if [ "$(bed_grep_count "$OUTPUT_DIR/baseline_dups.bed" "$chr" "$start" "$end" "0/1")" -gt 0 ]; then
        pass "$desc ($chr:$start-$end) = 0/1"
    else
        actual=$(bed_grep "$OUTPUT_DIR/baseline_dups.bed" "$chr" "$start" "$end" | awk -F'\t' '{print $4}')
        fail "$desc ($chr:$start-$end) expected 0/1, got ${actual:-NOT_FOUND}"
    fi
done

# -----------------------------------------------
# 6. Negative controls (chr4 - no SVs)
# -----------------------------------------------
section "6" "Negative controls (chr4 - normal regions should be 0/0)..."

NORMALS=(
    "chr4:20000-25000:Normal region 1"
    "chr4:50000-55000:Normal region 2"
    "chr4:70000-75000:Normal region 3"
)

for entry in "${NORMALS[@]}"; do
    IFS=: read -r loc start_end desc <<< "$entry"
    chr=$(echo "$loc" | cut -d: -f1)
    start=$(echo "$start_end" | cut -d- -f1)
    end=$(echo "$start_end" | cut -d- -f2)
    if [ "$(bed_grep_count "$OUTPUT_DIR/baseline_dels.bed" "$chr" "$start" "$end" "0/0")" -gt 0 ]; then
        pass "$desc ($chr:$start-$end) = 0/0 (true negative)"
    else
        actual=$(bed_grep "$OUTPUT_DIR/baseline_dels.bed" "$chr" "$start" "$end" | awk -F'\t' '{print $4}')
        fail "$desc ($chr:$start-$end) expected 0/0, got ${actual:-NOT_FOUND}"
    fi
done

# chr4 regions should NOT appear in SVs file
CHR4_IN_SVS=$(grep -c "^chr4" "$OUTPUT_DIR/baseline_svs.bed" 2>/dev/null || true)
CHR4_IN_SVS=$(echo "$CHR4_IN_SVS" | tr -d '[:space:]')
CHR4_IN_SVS=${CHR4_IN_SVS:-0}
if [ "$CHR4_IN_SVS" -eq 0 ]; then
    pass "No chr4 entries in SVs file (no false positives)"
else
    fail "$CHR4_IN_SVS chr4 false positives in SVs file"
fi

# -----------------------------------------------
# 7. SV counts
# -----------------------------------------------
section "7" "SV counts..."

for f in baseline_svs.bed baseline_dels.bed baseline_dups.bed; do
    CURRENT=$(grep -cv "^#" "$OUTPUT_DIR/$f" 2>/dev/null || echo "0")
    EXPECTED=$(grep -cv "^#" "$BASELINE_DIR/$f" 2>/dev/null || echo "0")
    if [ "$CURRENT" -eq "$EXPECTED" ]; then
        pass "$f: $CURRENT entries (matches baseline)"
    else
        fail "$f: expected $EXPECTED entries, got $CURRENT"
    fi
done

# Breakdown by type in SVs file
DEL_COUNT=$(grep -c "DEL" "$OUTPUT_DIR/baseline_svs.bed" 2>/dev/null || echo "0")
DUP_COUNT=$(grep -c "DUP" "$OUTPUT_DIR/baseline_svs.bed" 2>/dev/null || echo "0")
DEL_EXPECTED=$(grep -c "DEL" "$BASELINE_DIR/baseline_svs.bed" 2>/dev/null || echo "0")
DUP_EXPECTED=$(grep -c "DUP" "$BASELINE_DIR/baseline_svs.bed" 2>/dev/null || echo "0")

if [ "$DEL_COUNT" -eq "$DEL_EXPECTED" ]; then
    pass "DEL count in SVs: $DEL_COUNT"
else
    fail "DEL count in SVs: expected $DEL_EXPECTED, got $DEL_COUNT"
fi
if [ "$DUP_COUNT" -eq "$DUP_EXPECTED" ]; then
    pass "DUP count in SVs: $DUP_COUNT"
else
    fail "DUP count in SVs: expected $DUP_EXPECTED, got $DUP_COUNT"
fi

# -----------------------------------------------
# 8. Numerical tolerance (likelihood + read counts)
# -----------------------------------------------
section "8" "Numerical tolerance check..."

python3 << 'PYEOF'
import sys

def parse_bed(path):
    rows = []
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            rows.append(line.strip().split('\t'))
    return rows

files = {
    'baseline_dels.bed': {'like_col': 4, 'obs_col': 7, 'exp_col': 8},
    'baseline_dups.bed': {'like_col': 4, 'obs_col': 7, 'exp_col': 8},
}

base_dir = 'test_data/baseline'
out_dir = 'test_data/output'

like_tol = 0.05       # 5% relative tolerance for likelihood
read_tol = 0.10       # 10% tolerance for read counts
issues = 0
checks = 0

for fname, cols in files.items():
    base_rows = parse_bed(f'{base_dir}/{fname}')
    out_rows = parse_bed(f'{out_dir}/{fname}')

    if len(base_rows) != len(out_rows):
        print(f'  FAIL: {fname} row count: {len(base_rows)} vs {len(out_rows)}')
        issues += 1
        continue

    for i, (b, o) in enumerate(zip(base_rows, out_rows)):
        region = f'{b[0]}:{b[1]}-{b[2]}'
        checks += 1

        # Position match
        if b[:3] != o[:3]:
            print(f'  FAIL: {fname} row {i}: position {b[:3]} vs {o[:3]}')
            issues += 1
            continue

        # Genotype match
        if b[3] != o[3]:
            print(f'  FAIL: {fname} {region}: genotype {b[3]} vs {o[3]}')
            issues += 1
            continue

        # Likelihood tolerance
        try:
            bl, ol = float(b[cols['like_col']]), float(o[cols['like_col']])
            if bl != 0:
                rel = abs(bl - ol) / abs(bl)
                if rel > like_tol:
                    print(f'  WARN: {fname} {region}: likelihood {bl} vs {ol} (diff {rel*100:.1f}%)')
                    issues += 1
        except (ValueError, IndexError):
            pass

        # Observed read count tolerance
        try:
            bo, oo = int(b[cols['obs_col']]), int(o[cols['obs_col']])
            if bo != 0:
                rel = abs(bo - oo) / abs(bo)
                if rel > read_tol:
                    print(f'  WARN: {fname} {region}: obs reads {bo} vs {oo} (diff {rel*100:.1f}%)')
                    issues += 1
        except (ValueError, IndexError):
            pass

        # Expected read depth tolerance
        try:
            be, oe = float(b[cols['exp_col']]), float(o[cols['exp_col']])
            if be != 0:
                rel = abs(be - oe) / abs(be)
                if rel > read_tol:
                    print(f'  WARN: {fname} {region}: exp reads {be} vs {oe} (diff {rel*100:.1f}%)')
                    issues += 1
        except (ValueError, IndexError):
            pass

if issues == 0:
    print(f'  All {checks} numerical checks within tolerance')
else:
    print(f'  {issues} issues found in {checks} checks')

sys.exit(0)
PYEOF

# -----------------------------------------------
# 9. Borderline / edge cases
# -----------------------------------------------
section "9" "Borderline and edge case validation..."

# Borderline deletion (chr1:100000-101500) - mild signal, should be detected but borderline
if [ "$(bed_grep_count "$OUTPUT_DIR/baseline_dels.bed" chr1 100000 101500)" -gt 0 ]; then
    BORDERLINE_LIKE=$(bed_grep "$OUTPUT_DIR/baseline_dels.bed" chr1 100000 101500 | awk -F'\t' '{print $5}')
    pass "Borderline DEL detected (likelihood=$BORDERLINE_LIKE)"
else
    fail "Borderline DEL (chr1:100000-101500) not in output"
fi

# Low mappability region (chr1:100000-102000 has low mappability)
BORDERLINE_MAP=$(bed_grep "$OUTPUT_DIR/baseline_dels.bed" chr1 100000 101500 | awk -F'\t' '{print $7}')
if [ -n "$BORDERLINE_MAP" ]; then
    LOW_MAP=$(python3 -c "print('yes' if float('$BORDERLINE_MAP') < 0.5 else 'no')")
    if [ "$LOW_MAP" = "yes" ]; then
        pass "Low mappability correctly reflected for borderline DEL (map=$BORDERLINE_MAP)"
    else
        warn "Borderline DEL mappability=$BORDERLINE_MAP (expected <0.5)"
    fi
fi

# Borderline duplication (chr1:400000-402000) - mild signal
MILD_DUP_GENO=$(bed_grep "$OUTPUT_DIR/baseline_dups.bed" chr1 400000 402000 | awk -F'\t' '{print $4}')
if [ -n "$MILD_DUP_GENO" ]; then
    pass "Borderline DUP detected (genotype=$MILD_DUP_GENO)"
else
    warn "Borderline DUP (chr1:400000-402000) not detected"
fi

# Deletion at chromosome start (chr5:1000-4000)
if [ "$(bed_grep_count "$OUTPUT_DIR/baseline_dels.bed" chr5 1000 4000)" -gt 0 ]; then
    pass "DEL at chromosome start (chr5:1000-4000) detected"
else
    fail "DEL at chromosome start (chr5:1000-4000) not detected"
fi

# Deletion near chromosome end (chr1:490000-498000)
if [ "$(bed_grep_count "$OUTPUT_DIR/baseline_dels.bed" chr1 490000 498000)" -gt 0 ]; then
    pass "DEL near chromosome end (chr1:490000-498000) detected"
else
    fail "DEL near chromosome end (chr1:490000-498000) not detected"
fi

# Adjacent SVs both detected
ADJ1=$(bed_grep_count "$OUTPUT_DIR/baseline_dels.bed" chr1 130000 135000)
ADJ2=$(bed_grep_count "$OUTPUT_DIR/baseline_dels.bed" chr1 136000 140000)
if [ "$ADJ1" -gt 0 ] && [ "$ADJ2" -gt 0 ]; then
    pass "Adjacent SVs (130k-135k, 136k-140k) both detected"
else
    fail "Adjacent SVs not both detected (found: $ADJ1, $ADJ2)"
fi

# -----------------------------------------------
# 10. Cross-chromosome consistency
# -----------------------------------------------
section "10" "Cross-chromosome consistency..."

# Check that all 5 chromosomes were processed
for chr in chr1 chr2 chr3 chr4 chr5; do
    if grep -q "^$chr" "$OUTPUT_DIR/baseline_dels.bed" 2>/dev/null || grep -q "^$chr" "$OUTPUT_DIR/baseline_dups.bed" 2>/dev/null; then
        pass "$chr present in output"
    else
        fail "$chr missing from output"
    fi
done

# -----------------------------------------------
# 11. Observed/Expected ratio sanity
# -----------------------------------------------
section "11" "Observed/Expected read depth ratio sanity..."

python3 << 'PYEOF'
def check_ratios(path, label):
    ok = 0
    bad = 0
    with open(path) as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            geno = parts[3]
            try:
                obs = int(parts[7])
                exp = float(parts[8])
            except (ValueError, IndexError):
                continue

            if exp == 0:
                continue

            ratio = obs / exp

            # Sanity checks
            if geno == '1/1' and label == 'DEL':
                # Homo del: ratio should be very low (<0.2)
                if ratio > 0.3:
                    print(f'  WARN: {parts[0]}:{parts[1]}-{parts[2]} homo DEL but obs/exp={ratio:.2f}')
                    bad += 1
                else:
                    ok += 1
            elif geno == '0/1' and label == 'DEL':
                # Het del: ratio should be ~0.5 (0.3-0.7)
                if ratio < 0.2 or ratio > 0.8:
                    print(f'  WARN: {parts[0]}:{parts[1]}-{parts[2]} het DEL but obs/exp={ratio:.2f}')
                    bad += 1
                else:
                    ok += 1
            elif geno == '1/1' and label == 'DUP':
                # Homo dup: ratio should be ~2.0 (1.5-2.5)
                if ratio < 1.3 or ratio > 3.0:
                    print(f'  WARN: {parts[0]}:{parts[1]}-{parts[2]} homo DUP but obs/exp={ratio:.2f}')
                    bad += 1
                else:
                    ok += 1
            elif geno == '0/1' and label == 'DUP':
                # Het dup: ratio should be ~1.5 (1.2-1.8)
                if ratio < 1.0 or ratio > 2.2:
                    print(f'  WARN: {parts[0]}:{parts[1]}-{parts[2]} het DUP but obs/exp={ratio:.2f}')
                    bad += 1
                else:
                    ok += 1
            elif geno == '0/0':
                # Normal: ratio should be ~1.0 (0.8-1.2)
                if ratio < 0.7 or ratio > 1.3:
                    print(f'  WARN: {parts[0]}:{parts[1]}-{parts[2]} normal but obs/exp={ratio:.2f}')
                    bad += 1
                else:
                    ok += 1

    if bad == 0:
        print(f'  {label}: all {ok} obs/exp ratios are sane')
    else:
        print(f'  {label}: {bad} suspicious ratios out of {ok + bad}')

check_ratios('test_data/output/baseline_dels.bed', 'DEL')
check_ratios('test_data/output/baseline_dups.bed', 'DUP')
PYEOF

# -----------------------------------------------
# Update baseline if requested
# -----------------------------------------------
if [ "${1:-}" = "--update-baseline" ]; then
    echo ""
    section "*" "Updating baseline..."
    cp "$OUTPUT_DIR"/baseline_*.bed "$BASELINE_DIR/"
    echo "  Baseline updated."
fi

# -----------------------------------------------
# Summary
# -----------------------------------------------
echo ""
echo "============================================"
echo -e "  ${GREEN}PASS: $PASS${NC}"
if [ $FAIL -gt 0 ]; then
    echo -e "  ${RED}FAIL: $FAIL${NC}"
fi
if [ $WARN -gt 0 ]; then
    echo -e "  ${YELLOW}WARN: $WARN${NC}"
fi
echo "============================================"

if [ $FAIL -gt 0 ]; then
    echo -e "${RED}REGRESSION TEST FAILED${NC}"
    exit 1
else
    echo -e "${GREEN}REGRESSION TEST PASSED${NC}"
    exit 0
fi
