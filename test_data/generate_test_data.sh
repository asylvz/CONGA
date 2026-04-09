#!/bin/bash
# Generate comprehensive synthetic test data for CONGA regression testing
# Covers: various SV sizes, coverage levels, GC content, edge cases, split-reads

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
CONGA_DIR="$(dirname "$SCRIPT_DIR")"
DATA_DIR="$SCRIPT_DIR"

SEED=42

echo "=== CONGA Comprehensive Test Data Generator ==="
echo "Output directory: $DATA_DIR"

# -----------------------------------------------
# 1. Generate synthetic reference genome (FASTA)
#    - chr1: 500kb, mixed GC content
#    - chr2: 300kb, AT-rich regions
#    - chr3: 200kb, GC-rich regions
#    - chr4: 100kb, normal (no SVs - negative control)
#    - chr5: 150kb, edge cases
# -----------------------------------------------
echo "[1/7] Generating reference genome (5 chromosomes, 1.25 Mb total)..."

python3 << 'PYEOF'
import random
random.seed(42)

chroms = {
    'chr1': 500000,
    'chr2': 300000,
    'chr3': 200000,
    'chr4': 100000,
    'chr5': 150000,
}

def gen_seq(length, gc_bias=0.5):
    """Generate sequence with controllable GC content."""
    bases_gc = 'GC'
    bases_at = 'AT'
    seq = []
    for _ in range(length):
        if random.random() < gc_bias:
            seq.append(random.choice(bases_gc))
        else:
            seq.append(random.choice(bases_at))
    return ''.join(seq)

with open('test_data/ref.fa', 'w') as f:
    for chrom, length in chroms.items():
        if chrom == 'chr2':
            # AT-rich: 70% AT
            seq = gen_seq(length, gc_bias=0.3)
        elif chrom == 'chr3':
            # GC-rich: 70% GC
            seq = gen_seq(length, gc_bias=0.7)
        elif chrom == 'chr5':
            # Mixed: alternating GC-rich and AT-rich blocks (10kb each)
            seq = ''
            for i in range(0, length, 10000):
                block_len = min(10000, length - i)
                if (i // 10000) % 2 == 0:
                    seq += gen_seq(block_len, gc_bias=0.7)
                else:
                    seq += gen_seq(block_len, gc_bias=0.3)
        else:
            seq = gen_seq(length, gc_bias=0.5)

        f.write(f'>{chrom}\n')
        for i in range(0, length, 80):
            f.write(seq[i:i+80] + '\n')

total = sum(chroms.values())
print(f'Generated reference: {len(chroms)} chromosomes, {total:,} bp total')
PYEOF

samtools faidx "$DATA_DIR/ref.fa"
echo "  Reference indexed."

# -----------------------------------------------
# 2. Generate synthetic reads with realistic patterns
# -----------------------------------------------
echo "[2/7] Generating synthetic reads..."

python3 << 'PYEOF'
import random
random.seed(42)

# Load reference
ref_seqs = {}
chrom = None
with open('test_data/ref.fa') as f:
    for line in f:
        line = line.strip()
        if line.startswith('>'):
            chrom = line[1:].split()[0]
            ref_seqs[chrom] = []
        else:
            ref_seqs[chrom].append(line)
for c in ref_seqs:
    ref_seqs[c] = ''.join(ref_seqs[c])

# =============================================
# SV definitions — comprehensive test coverage
# =============================================
# Format: (chrom, start, end, sv_type, description)
#
# sv_type controls read generation:
#   'homo_del'  -> 0 copies (skip all reads)
#   'het_del'   -> 1 copy   (skip ~50% reads)
#   'mild_del'  -> slight reduction (~30% fewer reads, borderline call)
#   'homo_dup'  -> 4 copies (3x extra reads)
#   'het_dup'   -> 3 copies (1.5x reads)
#   'mild_dup'  -> slight increase (~30% more reads, borderline call)
#   'normal'    -> no change (false positive test)

sv_regions = [
    # ---- chr1: Various deletion sizes ----
    ('chr1',  10000,  13000, 'homo_del',  'Small homozygous deletion (3kb)'),
    ('chr1',  20000,  25000, 'het_del',   'Medium heterozygous deletion (5kb)'),
    ('chr1',  40000,  55000, 'homo_del',  'Large homozygous deletion (15kb)'),
    ('chr1',  70000,  72000, 'het_del',   'Small heterozygous deletion (2kb)'),
    ('chr1', 100000, 101500, 'mild_del',  'Borderline deletion (1.5kb, mild signal)'),
    # Adjacent SVs (tests boundary handling)
    ('chr1', 130000, 135000, 'homo_del',  'Adjacent DEL part 1'),
    ('chr1', 136000, 140000, 'het_del',   'Adjacent DEL part 2 (1kb gap)'),
    # Large SV
    ('chr1', 200000, 250000, 'het_del',   'Very large heterozygous deletion (50kb)'),
    # Near chromosome end
    ('chr1', 490000, 498000, 'homo_del',  'Deletion near chromosome end'),

    # ---- chr1: Various duplication sizes ----
    ('chr1', 300000, 305000, 'homo_dup',  'Medium homozygous duplication (5kb)'),
    ('chr1', 320000, 323000, 'het_dup',   'Small heterozygous duplication (3kb)'),
    ('chr1', 350000, 370000, 'homo_dup',  'Large homozygous duplication (20kb)'),
    ('chr1', 400000, 402000, 'mild_dup',  'Borderline duplication (2kb, mild signal)'),
    ('chr1', 420000, 430000, 'het_dup',   'Medium heterozygous duplication (10kb)'),
    ('chr1', 460000, 465000, 'het_dup',   'Duplication near chromosome end region'),

    # ---- chr2: AT-rich region SVs ----
    ('chr2',  30000,  35000, 'homo_del',  'Homo DEL in AT-rich region'),
    ('chr2',  60000,  65000, 'het_del',   'Het DEL in AT-rich region'),
    ('chr2', 100000, 105000, 'homo_dup',  'Homo DUP in AT-rich region'),
    ('chr2', 150000, 155000, 'het_dup',   'Het DUP in AT-rich region'),
    ('chr2', 200000, 210000, 'het_del',   'Large het DEL in AT-rich (10kb)'),
    ('chr2', 250000, 260000, 'homo_dup',  'Large homo DUP in AT-rich (10kb)'),

    # ---- chr3: GC-rich region SVs ----
    ('chr3',  20000,  25000, 'homo_del',  'Homo DEL in GC-rich region'),
    ('chr3',  50000,  55000, 'het_del',   'Het DEL in GC-rich region'),
    ('chr3',  80000,  85000, 'homo_dup',  'Homo DUP in GC-rich region'),
    ('chr3', 120000, 125000, 'het_dup',   'Het DUP in GC-rich region'),
    ('chr3', 160000, 165000, 'het_del',   'Het DEL near end GC-rich (5kb)'),

    # ---- chr4: No SVs (negative control) ----
    # All known_SVs BED entries for chr4 will be normal regions
    # CONGA should call 0/0 for these
    ('chr4',  20000,  25000, 'normal',    'Normal region 1 (should be 0/0)'),
    ('chr4',  50000,  55000, 'normal',    'Normal region 2 (should be 0/0)'),
    ('chr4',  70000,  75000, 'normal',    'Normal region 3 (should be 0/0)'),

    # ---- chr5: Edge cases (mixed GC blocks) ----
    ('chr5',   1000,   4000, 'homo_del',  'DEL at chromosome start'),
    ('chr5',  50000,  55000, 'het_del',   'DEL spanning GC/AT boundary'),
    ('chr5',  90000,  95000, 'homo_dup',  'DUP in mixed GC region'),
    ('chr5', 140000, 148000, 'het_dup',   'DUP near chromosome end'),
]

# Coverage target per base (controls read density)
# ~10x average coverage for a bioinformatics low-coverage scenario
coverage_target = {
    'chr1': 10,
    'chr2': 8,   # slightly lower coverage
    'chr3': 12,  # slightly higher coverage
    'chr4': 10,
    'chr5': 6,   # low coverage (ancient DNA simulation)
}

def mutate_read(seq, error_rate=0.005):
    bases = 'ACGT'
    result = list(seq)
    for i in range(len(result)):
        if random.random() < error_rate:
            result[i] = random.choice([b for b in bases if b != result[i]])
    return ''.join(result)

def get_sv_at_pos(chrom, pos, read_len):
    """Check if position falls within an SV region."""
    for sv in sv_regions:
        sv_chr, sv_start, sv_end, sv_type, _ = sv
        if sv_chr == chrom and sv_start <= pos <= sv_end - read_len:
            return sv_type
    return 'normal'

# Different read lengths to test min_read_length filtering
read_lengths = [75, 100, 100, 100, 150]  # weighted toward 100

with open('test_data/reads.sam', 'w') as f:
    # SAM header
    for chrom, seq in ref_seqs.items():
        f.write(f'@SQ\tSN:{chrom}\tLN:{len(seq)}\n')
    f.write('@RG\tID:test\tSM:TestSample\tLB:lib1\tPL:ILLUMINA\n')

    read_id = 0
    for chrom, seq in ref_seqs.items():
        chr_len = len(seq)
        cov = coverage_target[chrom]
        read_len_base = 100
        num_reads = int(chr_len * cov / read_len_base)

        for _ in range(num_reads):
            read_len = random.choice(read_lengths)
            pos = random.randint(0, chr_len - read_len - 1)

            sv_type = get_sv_at_pos(chrom, pos, read_len)

            # Modify read generation based on SV type
            if sv_type == 'homo_del':
                continue  # skip all reads
            elif sv_type == 'het_del':
                if random.random() < 0.5:
                    continue  # skip ~50%
            elif sv_type == 'mild_del':
                if random.random() < 0.3:
                    continue  # skip ~30%

            read_seq = seq[pos:pos + read_len]
            if len(read_seq) < read_len:
                continue

            read_seq = mutate_read(read_seq)
            qual = ''.join(chr(random.randint(25, 40) + 33) for _ in range(read_len))

            # Vary mapping quality realistically
            mapq = random.choices([0, 10, 20, 30, 40, 50, 60],
                                   weights=[1, 2, 5, 15, 25, 30, 22])[0]

            f.write(f'read_{read_id}\t0\t{chrom}\t{pos+1}\t{mapq}\t{read_len}M\t*\t0\t0\t{read_seq}\t{qual}\tRG:Z:test\n')
            read_id += 1

            # Add extra reads for duplications
            if sv_type == 'homo_dup':
                # ~2x total: add 1 extra read
                for _ in range(1):
                    read_seq2 = mutate_read(seq[pos:pos + read_len])
                    qual2 = ''.join(chr(random.randint(25, 40) + 33) for _ in range(read_len))
                    mapq2 = random.choices([30, 40, 50, 60], weights=[15, 25, 30, 30])[0]
                    f.write(f'read_{read_id}\t0\t{chrom}\t{pos+1}\t{mapq2}\t{read_len}M\t*\t0\t0\t{read_seq2}\t{qual2}\tRG:Z:test\n')
                    read_id += 1
            elif sv_type == 'het_dup':
                if random.random() < 0.5:
                    read_seq2 = mutate_read(seq[pos:pos + read_len])
                    qual2 = ''.join(chr(random.randint(25, 40) + 33) for _ in range(read_len))
                    mapq2 = random.choices([30, 40, 50, 60], weights=[15, 25, 30, 30])[0]
                    f.write(f'read_{read_id}\t0\t{chrom}\t{pos+1}\t{mapq2}\t{read_len}M\t*\t0\t0\t{read_seq2}\t{qual2}\tRG:Z:test\n')
                    read_id += 1
            elif sv_type == 'mild_dup':
                if random.random() < 0.3:
                    read_seq2 = mutate_read(seq[pos:pos + read_len])
                    qual2 = ''.join(chr(random.randint(25, 40) + 33) for _ in range(read_len))
                    mapq2 = random.choices([30, 40, 50, 60], weights=[15, 25, 30, 30])[0]
                    f.write(f'read_{read_id}\t0\t{chrom}\t{pos+1}\t{mapq2}\t{read_len}M\t*\t0\t0\t{read_seq2}\t{qual2}\tRG:Z:test\n')
                    read_id += 1

print(f'Generated {read_id:,} reads across 5 chromosomes')
PYEOF

samtools view -bS "$DATA_DIR/reads.sam" | samtools sort -o "$DATA_DIR/test.bam"
samtools index "$DATA_DIR/test.bam"
rm "$DATA_DIR/reads.sam"
echo "  BAM file created and indexed."

# -----------------------------------------------
# 3. Generate BED files for known SVs
# -----------------------------------------------
echo "[3/7] Generating known SV BED files..."

cat > "$DATA_DIR/known_dels.bed" << 'EOF'
chr1	10000	13000
chr1	20000	25000
chr1	40000	55000
chr1	70000	72000
chr1	100000	101500
chr1	130000	135000
chr1	136000	140000
chr1	200000	250000
chr1	490000	498000
chr2	30000	35000
chr2	60000	65000
chr2	200000	210000
chr3	20000	25000
chr3	50000	55000
chr3	160000	165000
chr4	20000	25000
chr4	50000	55000
chr4	70000	75000
chr5	1000	4000
chr5	50000	55000
EOF

cat > "$DATA_DIR/known_dups.bed" << 'EOF'
chr1	300000	305000
chr1	320000	323000
chr1	350000	370000
chr1	400000	402000
chr1	420000	430000
chr1	460000	465000
chr2	100000	105000
chr2	150000	155000
chr2	250000	260000
chr3	80000	85000
chr3	120000	125000
chr5	90000	95000
chr5	140000	148000
EOF

echo "  $(wc -l < "$DATA_DIR/known_dels.bed" | tr -d ' ') DELs and $(wc -l < "$DATA_DIR/known_dups.bed" | tr -d ' ') DUPs defined."

# -----------------------------------------------
# 4. Generate SONIC file
# -----------------------------------------------
echo "[4/7] Generating SONIC file..."

# Gaps file
cat > "$DATA_DIR/gaps.bed" << 'EOF'
chr1	0	500
chr2	0	500
chr3	0	500
chr4	0	500
chr5	0	500
EOF

# RepeatMasker .out format
cat > "$DATA_DIR/repeats.out" << 'EOF'
   SW  perc perc perc  query      position in query           matching       repeat              position in  repeat
score  div. del. ins.  sequence    begin     end    (left)    repeat         class/family         begin  end (left)   ID

  200   5.0  0.0  0.0  chr1         1000     1100 (498900) +  AluSx          SINE/Alu               1    100  (200)    1
  180   6.0  0.0  0.0  chr2         2000     2100 (297900) +  L1             LINE/L1                1    100  (200)    2
  150   4.0  0.0  0.0  chr3         3000     3100 (196900) +  AluY           SINE/Alu               1    100  (200)    3
EOF

# Segmental duplications
cat > "$DATA_DIR/segdups.bed" << 'EOF'
chr1	300000	305000
chr1	350000	370000
chr3	80000	85000
EOF

"$CONGA_DIR/sonic/sonic" --ref "$DATA_DIR/ref.fa" \
        --gaps "$DATA_DIR/gaps.bed" \
        --reps "$DATA_DIR/repeats.out" \
        --dups "$DATA_DIR/segdups.bed" \
        --make-sonic "$DATA_DIR/test.sonic" 2>&1 || {
    echo "WARNING: sonic build failed."
    exit 1
}

echo "  SONIC file created."

# -----------------------------------------------
# 5. Generate mappability file
# -----------------------------------------------
echo "[5/7] Generating mappability file..."

python3 << 'PYEOF'
import random
random.seed(43)

chroms = {'chr1': 500000, 'chr2': 300000, 'chr3': 200000, 'chr4': 100000, 'chr5': 150000}
window = 500

# Define some low-mappability regions for testing
low_map_regions = {
    'chr1': [(100000, 102000)],    # overlaps with borderline deletion
    'chr3': [(80000, 82000)],      # overlaps with GC-rich DUP
    'chr5': [(90000, 92000)],      # overlaps with edge case DUP
}

with open('test_data/mappability.bed', 'w') as f:
    for chrom, length in chroms.items():
        for start in range(0, length, window):
            end = min(start + window, length)
            # Check if in low-mappability region
            is_low = False
            for lr in low_map_regions.get(chrom, []):
                if start < lr[1] and end > lr[0]:
                    is_low = True
                    break
            if is_low:
                score = round(random.uniform(0.1, 0.4), 4)
            else:
                score = round(random.uniform(0.7, 1.0), 4)
            f.write(f'{chrom}\t{start}\t{end}\t{score}\n')

print('Mappability file created (with low-map regions for edge cases)')
PYEOF

echo "  Mappability file created."

# -----------------------------------------------
# 6. Generate low-mappability exclusion regions
# -----------------------------------------------
echo "[6/7] Generating low-mappability exclusion file..."

cat > "$DATA_DIR/low_map_exclude.bed" << 'EOF'
chr1	100000	102000
chr3	80000	82000
chr5	90000	92000
EOF

echo "  Low-mappability exclusion file created."

# -----------------------------------------------
# 7. Verify all files
# -----------------------------------------------
echo "[7/7] Verifying test data..."

EXPECTED_FILES="ref.fa ref.fa.fai test.bam test.bam.bai known_dels.bed known_dups.bed test.sonic mappability.bed low_map_exclude.bed"
ALL_OK=true

for f in $EXPECTED_FILES; do
    if [ -f "$DATA_DIR/$f" ]; then
        SIZE=$(ls -lh "$DATA_DIR/$f" | awk '{print $5}')
        echo "  OK: $f ($SIZE)"
    else
        echo "  MISSING: $f"
        ALL_OK=false
    fi
done

# Stats
echo ""
samtools flagstat "$DATA_DIR/test.bam" 2>/dev/null | head -3
echo ""

if $ALL_OK; then
    echo "=== All test data generated successfully ==="
else
    echo "=== Some files are missing! ==="
    exit 1
fi
