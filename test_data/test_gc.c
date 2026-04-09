#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

/*
 * Standalone test for GC content computation.
 * Replicates the exact SONIC algorithm:
 *   - char_count starts at 1 (SONIC sets char_count=1 after FASTA header)
 *   - First window covers 99 bases (bases 0-98), subsequent windows 100 bases each
 *   - GC% = (G+C count) * 100 / 100  (N bases count in denominator but not numerator)
 *   - Truncation to char (not rounding) during storage
 *   - Partial last window is dropped (not written to profile)
 *   - Callers use round() on retrieval via get_gc_content
 */

#define GC_WINDOW 100

/* This is the function we're testing - matches genome.c which replicates SONIC exactly */
static char* compute_gc_profile(const char *seq, int seq_len, int *profile_len)
{
	int num_windows = (seq_len / GC_WINDOW) + 1;
	char *profile = (char*) calloc(num_windows, sizeof(char));
	int gc = 0;
	int char_count = 1; /* Match SONIC: char_count starts at 1 after FASTA header */
	int window_id = 0;

	int i;
	for(i = 0; i < seq_len; i++)
	{
		char c = seq[i];
		if(c == 'G' || c == 'g' || c == 'C' || c == 'c')
			gc++;
		char_count++;

		if(char_count % GC_WINDOW == 0)
		{
			profile[window_id] = (char)(100.0 * gc / GC_WINDOW);
			window_id++;
			gc = 0;
		}
	}

	/* Match SONIC: partial last window is dropped (not written) */

	*profile_len = window_id;
	return profile;
}

/* Replicate sonic_get_gc_content: returns float GC% for a range */
static float get_gc_content(const char *profile, int pos_start, int pos_end)
{
	int window_count = 0;
	float gc_content = 0.0;
	int start_gc = pos_start;

	while(start_gc < pos_end)
	{
		window_count++;
		gc_content += (float) profile[start_gc / GC_WINDOW];
		start_gc += GC_WINDOW;
	}

	if(gc_content != 0.0)
		return gc_content / window_count;
	else
		return 0.0;
}

static int tests_passed = 0;
static int tests_failed = 0;

#define CHECK(cond, msg) do { \
	if(cond) { tests_passed++; printf("  PASS: %s\n", msg); } \
	else { tests_failed++; printf("  FAIL: %s\n", msg); } \
} while(0)

#define CHECK_INT(val, expected, msg) do { \
	int _v = (val), _e = (expected); \
	if(_v == _e) { tests_passed++; printf("  PASS: %s (got %d)\n", msg, _v); } \
	else { tests_failed++; printf("  FAIL: %s (expected %d, got %d)\n", msg, _e, _v); } \
} while(0)

#define CHECK_FLOAT(val, expected, tol, msg) do { \
	float _v = (val), _e = (expected); \
	if(fabs(_v - _e) <= tol) { tests_passed++; printf("  PASS: %s (got %.2f)\n", msg, _v); } \
	else { tests_failed++; printf("  FAIL: %s (expected %.2f, got %.2f)\n", msg, _e, _v); } \
} while(0)

/* Helper: build a sequence with given GC fraction */
static char* make_seq(int length, float gc_fraction)
{
	char *seq = (char*) malloc(length + 1);
	int gc_count = (int)(length * gc_fraction);
	int i;
	for(i = 0; i < length; i++)
	{
		if(i < gc_count / 2)
			seq[i] = 'G';
		else if(i < gc_count)
			seq[i] = 'C';
		else
			seq[i] = 'A';
	}
	seq[length] = '\0';
	return seq;
}

void test_basic_gc()
{
	printf("\n[1] Basic GC content tests\n");
	/*
	 * SONIC behavior: char_count=1, first window fires at char_count=100.
	 * So first window covers bases 0-98 (99 bases), denominator is always GC_WINDOW=100.
	 */

	/* 100bp all GC: window 0 covers 0-98 (99 G/C), base 99 is partial → dropped */
	char *seq = make_seq(100, 1.0);
	int plen;
	char *profile = compute_gc_profile(seq, 100, &plen);
	CHECK_INT(plen, 1, "100bp all-GC: 1 window");
	CHECK_INT((int)profile[0], 99, "100bp all-GC: first window has 99 bases → GC=99");
	free(seq); free(profile);

	/* 100bp all AT */
	seq = make_seq(100, 0.0);
	profile = compute_gc_profile(seq, 100, &plen);
	CHECK_INT((int)profile[0], 0, "100bp all-AT: GC=0%");
	free(seq); free(profile);

	/* 100bp 50% GC (G*25 + C*25 + A*50): all GC in bases 0-49, window covers 0-98 */
	seq = make_seq(100, 0.5);
	profile = compute_gc_profile(seq, 100, &plen);
	CHECK_INT((int)profile[0], 50, "100bp 50% GC: GC=50% (all GC within first 99 bases)");
	free(seq); free(profile);

	/* 100bp 40% GC */
	seq = make_seq(100, 0.4);
	profile = compute_gc_profile(seq, 100, &plen);
	CHECK_INT((int)profile[0], 40, "100bp 40% GC: GC=40%");
	free(seq); free(profile);

	/* 100bp 73% GC */
	seq = make_seq(100, 0.73);
	profile = compute_gc_profile(seq, 100, &plen);
	CHECK_INT((int)profile[0], 73, "100bp 73% GC: GC=73%");
	free(seq); free(profile);
}

void test_multiple_windows()
{
	printf("\n[2] Multiple window tests\n");

	/* 200bp: first 100bp = 30% GC, second 100bp = 70% GC
	 * SONIC windows: [0-98] (99 bases), [99-198] (100 bases), base 199 dropped
	 */
	char seq[201];
	int i;

	for(i = 0; i < 30; i++) seq[i] = 'G';
	for(i = 30; i < 100; i++) seq[i] = 'A';
	for(i = 100; i < 170; i++) seq[i] = 'C';
	for(i = 170; i < 200; i++) seq[i] = 'T';
	seq[200] = '\0';

	int plen;
	char *profile = compute_gc_profile(seq, 200, &plen);
	CHECK_INT(plen, 2, "200bp: 2 windows");
	/* Window 0: bases 0-98 = 30G + 69A = 30 GC */
	CHECK_INT((int)profile[0], 30, "Window 0: 30% GC");
	/* Window 1: bases 99-198 = A + C*70 + T*29 = 70 GC */
	CHECK_INT((int)profile[1], 70, "Window 1: 70% GC");

	float gc;
	gc = get_gc_content(profile, 0, 100);
	CHECK_FLOAT(gc, 30.0, 0.01, "get_gc_content(0, 100) = 30%");

	gc = get_gc_content(profile, 100, 200);
	CHECK_FLOAT(gc, 70.0, 0.01, "get_gc_content(100, 200) = 70%");

	gc = get_gc_content(profile, 50, 150);
	CHECK_FLOAT(gc, 30.0, 0.01, "get_gc_content(50, 150) = tile 0 only = 30%");

	free(profile);
}

void test_n_bases()
{
	printf("\n[3] N base handling tests\n");

	/* 100bp: 50 N's + 25 G's + 25 T's
	 * Window 0 covers 0-98: N*50 + G*25 + T*24 = 25 GC */
	char seq[101];
	int i;
	for(i = 0; i < 50; i++) seq[i] = 'N';
	for(i = 50; i < 75; i++) seq[i] = 'G';
	for(i = 75; i < 100; i++) seq[i] = 'T';
	seq[100] = '\0';

	int plen;
	char *profile = compute_gc_profile(seq, 100, &plen);
	CHECK_INT((int)profile[0], 25, "50N + 25G + 25T: GC=25% (denominator=100, not 50)");
	free(profile);

	/* All N's */
	for(i = 0; i < 100; i++) seq[i] = 'N';
	seq[100] = '\0';
	profile = compute_gc_profile(seq, 100, &plen);
	CHECK_INT((int)profile[0], 0, "100 N's: GC=0%");
	free(profile);

	/* 90 N's + 10 C's: Window 0 covers 0-98: N*90 + C*9 = 9 GC
	 * (base 99 = C but it falls outside the first 99-base window) */
	for(i = 0; i < 90; i++) seq[i] = 'N';
	for(i = 90; i < 100; i++) seq[i] = 'C';
	seq[100] = '\0';
	profile = compute_gc_profile(seq, 100, &plen);
	CHECK_INT((int)profile[0], 9, "90N + 10C: GC=9% (base 99 outside first window)");
	free(profile);
}

void test_truncation()
{
	printf("\n[4] Truncation vs rounding tests\n");

	/* 33 GC bases, all within 0-32 → within first 99-base window */
	char seq[101];
	int i;

	for(i = 0; i < 33; i++) seq[i] = 'G';
	for(i = 33; i < 100; i++) seq[i] = 'A';
	seq[100] = '\0';

	int plen;
	char *profile = compute_gc_profile(seq, 100, &plen);
	int stored = (int)profile[0];
	CHECK_INT(stored, 33, "33/100 stored as 33 (truncation)");

	float gc = get_gc_content(profile, 0, 100);
	int gc_val = (int) round(gc);
	CHECK_INT(gc_val, 33, "round(33.0) = 33");
	free(profile);
}

void test_sonic_window_boundaries()
{
	printf("\n[5] SONIC window boundary tests (char_count=1 offset)\n");

	/* Key SONIC behavior: first window is 99 bases, rest are 100 bases.
	 * Partial last window is dropped.
	 *
	 * For N bases: number of windows = floor((N+1)/100)
	 *   N=99:  floor(100/100)=1 (exactly 1 window of 99 bases)
	 *   N=100: floor(101/100)=1 (99 bases in window, 1 base dropped)
	 *   N=199: floor(200/100)=2
	 *   N=200: floor(201/100)=2
	 *   N=299: floor(300/100)=3
	 *   N=300: floor(301/100)=3
	 *   N=50:  floor(51/100)=0  (too short for any window)
	 */

	int plen;
	char seq[301];
	int i;

	/* 99bp: exactly fills first window */
	for(i = 0; i < 99; i++) seq[i] = 'G';
	seq[99] = '\0';
	char *profile = compute_gc_profile(seq, 99, &plen);
	CHECK_INT(plen, 1, "99bp: exactly 1 window (99 bases)");
	CHECK_INT((int)profile[0], 99, "99bp all-GC: 99 GC / 100 denominator = 99");
	free(profile);

	/* 199bp: 2 windows, first=99 bases, second=100 bases */
	for(i = 0; i < 199; i++) seq[i] = 'G';
	seq[199] = '\0';
	profile = compute_gc_profile(seq, 199, &plen);
	CHECK_INT(plen, 2, "199bp: 2 windows");
	CHECK_INT((int)profile[0], 99, "199bp window 0: 99 GC");
	CHECK_INT((int)profile[1], 100, "199bp window 1: 100 GC");
	free(profile);

	/* 200bp: still 2 windows (base 199 is partial → dropped) */
	for(i = 0; i < 200; i++) seq[i] = 'G';
	seq[200] = '\0';
	profile = compute_gc_profile(seq, 200, &plen);
	CHECK_INT(plen, 2, "200bp: 2 windows (base 199 dropped)");
	free(profile);

	/* 50bp: too short, 0 windows */
	for(i = 0; i < 50; i++) seq[i] = 'G';
	seq[50] = '\0';
	profile = compute_gc_profile(seq, 50, &plen);
	CHECK_INT(plen, 0, "50bp: 0 windows (too short for first window)");
	free(profile);

	/* 300bp: 3 windows */
	for(i = 0; i < 300; i++) seq[i] = 'G';
	seq[300] = '\0';
	profile = compute_gc_profile(seq, 300, &plen);
	CHECK_INT(plen, 3, "300bp: 3 windows");
	CHECK_INT((int)profile[0], 99, "300bp window 0: 99 (99 bases)");
	CHECK_INT((int)profile[1], 100, "300bp window 1: 100 (100 bases)");
	CHECK_INT((int)profile[2], 100, "300bp window 2: 100 (100 bases)");
	free(profile);

	/* Verify with non-uniform: 300bp with distinct tiles
	 * tile 0: 20 G + 80 A (bases 0-99, but window 0 covers 0-98)
	 * tile 1: 60 C + 40 T (bases 100-199, window 1 covers 99-198)
	 * tile 2: 90 G + 10 A (bases 200-299, window 2 covers 199-298)
	 */
	for(i = 0; i < 20; i++) seq[i] = 'G';
	for(i = 20; i < 100; i++) seq[i] = 'A';
	for(i = 100; i < 160; i++) seq[i] = 'C';
	for(i = 160; i < 200; i++) seq[i] = 'T';
	for(i = 200; i < 290; i++) seq[i] = 'G';
	for(i = 290; i < 300; i++) seq[i] = 'A';
	seq[300] = '\0';

	profile = compute_gc_profile(seq, 300, &plen);
	CHECK_INT(plen, 3, "300bp non-uniform: 3 windows");
	/* Window 0: bases 0-98 → 20G + 79A = 20 GC → 20 */
	CHECK_INT((int)profile[0], 20, "Non-uniform window 0: 20% (20G in 0-98)");
	/* Window 1: bases 99-198 → A + C*60 + T*39 = 60 GC → 60 */
	CHECK_INT((int)profile[1], 60, "Non-uniform window 1: 60% (A + 60C + 39T)");
	/* Window 2: bases 199-298 → T + G*90 + A*9 = 90 GC → 90 */
	CHECK_INT((int)profile[2], 90, "Non-uniform window 2: 90% (T + 90G + 9A)");
	free(profile);
}

void test_tile_indexing()
{
	printf("\n[6] Tile indexing tests (get_gc_content lookup)\n");

	/* Use the non-uniform 300bp profile from above */
	char seq[301];
	int i;
	for(i = 0; i < 20; i++) seq[i] = 'G';
	for(i = 20; i < 100; i++) seq[i] = 'A';
	for(i = 100; i < 160; i++) seq[i] = 'C';
	for(i = 160; i < 200; i++) seq[i] = 'T';
	for(i = 200; i < 290; i++) seq[i] = 'G';
	for(i = 290; i < 300; i++) seq[i] = 'A';
	seq[300] = '\0';

	int plen;
	char *profile = compute_gc_profile(seq, 300, &plen);

	/* get_gc_content steps by 100 from start_gc, reads profile[start_gc/100] */
	float gc;
	gc = get_gc_content(profile, 0, 100);
	CHECK_FLOAT(gc, 20.0, 0.01, "pos 0-99 -> tile 0 = 20%");

	gc = get_gc_content(profile, 50, 150);
	CHECK_FLOAT(gc, 20.0, 0.01, "pos 50-149 -> tile0 only (step: 50->150 stops) = 20%");

	gc = get_gc_content(profile, 99, 199);
	CHECK_FLOAT(gc, 20.0, 0.01, "pos 99-198 -> tile0 only (step: 99->199 stops) = 20%");

	gc = get_gc_content(profile, 150, 250);
	CHECK_FLOAT(gc, 60.0, 0.01, "pos 150-249 -> tile1 only (step: 150->250 stops) = 60%");

	gc = get_gc_content(profile, 0, 200);
	CHECK_FLOAT(gc, 40.0, 0.01, "pos 0-199 -> avg(tile0+tile1) = (20+60)/2 = 40%");

	gc = get_gc_content(profile, 100, 300);
	CHECK_FLOAT(gc, 75.0, 0.01, "pos 100-299 -> avg(tile1+tile2) = (60+90)/2 = 75%");

	/* Single-tile lookup (how CONGA actually calls it) */
	int gc_val = (int) round(get_gc_content(profile, 0, 100));
	CHECK_INT(gc_val, 20, "Single-tile lookup pos=0: gc_val=20");

	gc_val = (int) round(get_gc_content(profile, 155, 255));
	CHECK_INT(gc_val, 60, "Single-tile lookup pos=155: tile1 only = 60");

	free(profile);
}

void test_case_insensitivity()
{
	printf("\n[7] Case insensitivity tests\n");

	/* 100bp: g*25+c*25+a*25+t*25. Window 0 covers 0-98: g*25+c*25+a*24 = 50 GC */
	char seq[101];
	int i;
	for(i = 0; i < 25; i++) seq[i] = 'g';
	for(i = 25; i < 50; i++) seq[i] = 'c';
	for(i = 50; i < 75; i++) seq[i] = 'a';
	for(i = 75; i < 100; i++) seq[i] = 't';
	seq[100] = '\0';

	int plen;
	char *profile = compute_gc_profile(seq, 100, &plen);
	CHECK_INT((int)profile[0], 50, "Lowercase gcat: GC=50%");
	free(profile);

	/* Mixed case */
	for(i = 0; i < 30; i++) seq[i] = (i % 2 == 0) ? 'G' : 'g';
	for(i = 30; i < 100; i++) seq[i] = (i % 2 == 0) ? 'A' : 'a';
	seq[100] = '\0';
	profile = compute_gc_profile(seq, 100, &plen);
	CHECK_INT((int)profile[0], 30, "Mixed case GgAa: GC=30%");
	free(profile);
}

void test_real_sequence()
{
	printf("\n[8] Real-like sequence tests\n");

	/* ATCGATCG pattern (50% GC). 400bp = 4 windows.
	 * Window 0: bases 0-98 (99 bases). 12 full repeats (96 bases, 48 GC) + "ATC" (1 GC) = 49 GC
	 * Window 1-3: each 100 bases, starts mid-pattern, still ~50 GC each */
	char seq[401];
	int i;
	const char *pattern = "ATCGATCG";
	int plen_pat = 8;
	for(i = 0; i < 400; i++)
		seq[i] = pattern[i % plen_pat];
	seq[400] = '\0';

	int plen;
	char *profile = compute_gc_profile(seq, 400, &plen);
	CHECK_INT(plen, 4, "400bp: 4 windows");
	CHECK_INT((int)profile[0], 49, "ATCGATCG: tile 0 = 49% (99 bases: 12*4+1=49 GC)");
	CHECK_INT((int)profile[1], 50, "ATCGATCG: tile 1 = 50%");
	CHECK_INT((int)profile[2], 50, "ATCGATCG: tile 2 = 50%");
	CHECK_INT((int)profile[3], 50, "ATCGATCG: tile 3 = 50%");

	float gc = get_gc_content(profile, 250, 350);
	int gc_val = (int) round(gc);
	CHECK_INT(gc_val, 50, "get_gc_content for ATCGATCG region (tile 2) = 50");
	free(profile);
}

void test_chromosome_scale()
{
	printf("\n[9] Chromosome-scale tests\n");

	/* Simulate a small chromosome (10000bp) with uniform 40% GC */
	int len = 10000;
	char *seq = (char*) malloc(len + 1);
	int i;
	for(i = 0; i < 4000; i++) seq[i] = 'G';
	for(i = 4000; i < len; i++) seq[i] = 'A';
	seq[len] = '\0';

	int plen;
	char *profile = compute_gc_profile(seq, len, &plen);
	/* floor((10000+1)/100) = 100 windows */
	CHECK_INT(plen, 100, "10000bp: 100 windows");

	/* Window 0: bases 0-98, all G → 99 GC → 99 */
	CHECK_INT((int)profile[0], 99, "10000bp window 0: 99 (99 G bases)");

	/* Windows 1-39: bases 99-3998, all G → 100 GC → 100 */
	CHECK_INT((int)profile[1], 100, "10000bp window 1: 100 (all G)");
	CHECK_INT((int)profile[39], 100, "10000bp window 39: 100 (all G)");

	/* Window 40: bases 3999-4098. seq[3999]=G, seq[4000-4098]=A → 1 GC → 1 */
	CHECK_INT((int)profile[40], 1, "10000bp window 40: 1 (boundary G/A)");

	/* Windows 41-99: all A → 0 */
	CHECK_INT((int)profile[41], 0, "10000bp window 41: 0 (all A)");
	CHECK_INT((int)profile[99], 0, "10000bp window 99: 0 (all A)");

	/* Verify get_gc_content for SV-like region */
	float gc = get_gc_content(profile, 500, 600);
	CHECK_INT((int)round(gc), 100, "get_gc_content in GC-rich region = 100");

	gc = get_gc_content(profile, 5000, 5100);
	CHECK_INT((int)round(gc), 0, "get_gc_content in AT-rich region = 0");

	free(seq);
	free(profile);
}

void test_consistency_with_sonic()
{
	printf("\n[10] SONIC compatibility verification\n");

	/* SONIC behavior details:
	 * 1. char_count starts at 1 → first window = 99 bases
	 * 2. stored = (char)(100.0 * gc / GC_WINDOW) → truncation
	 * 3. Partial last window dropped
	 * 4. get_gc_content: identical stepping logic
	 *
	 * Test with a full 100-base window (standard case) */
	int plen;
	char seq[101];
	int i;

	/* 47 GC out of 99 bases (first window) */
	for(i = 0; i < 47; i++) seq[i] = 'G';
	for(i = 47; i < 100; i++) seq[i] = 'A';
	seq[100] = '\0';
	char *profile = compute_gc_profile(seq, 100, &plen);
	CHECK_INT((int)profile[0], 47, "Full window 47/100: exact = 47");
	free(profile);

	/* Verify the 1-base offset: 200bp where base 98 and 99 differ.
	 * This is where SONIC's offset matters most. */
	char seq2[201];
	for(i = 0; i < 98; i++) seq2[i] = 'A';
	seq2[98] = 'G';   /* This G is IN window 0 (base 98, last of 99 bases) */
	seq2[99] = 'C';   /* This C is IN window 1 (base 99, first of 100 bases) */
	for(i = 100; i < 200; i++) seq2[i] = 'A';
	seq2[200] = '\0';

	profile = compute_gc_profile(seq2, 200, &plen);
	CHECK_INT(plen, 2, "200bp offset test: 2 windows");
	CHECK_INT((int)profile[0], 1, "Offset test window 0: G at base 98 counted → 1 GC");
	CHECK_INT((int)profile[1], 1, "Offset test window 1: C at base 99 counted → 1 GC");
	free(profile);

	/* Verify SONIC's lookup function is byte-identical */
	char test_profile[] = {40, 60, 80};
	float gc = get_gc_content(test_profile, 0, 100);
	CHECK_FLOAT(gc, 40.0, 0.01, "Lookup profile[0]=40 → 40.0");
	gc = get_gc_content(test_profile, 0, 200);
	CHECK_FLOAT(gc, 50.0, 0.01, "Lookup avg(40,60) = 50.0");
	gc = get_gc_content(test_profile, 0, 300);
	CHECK_FLOAT(gc, 60.0, 0.01, "Lookup avg(40,60,80) = 60.0");
}

int main()
{
	printf("=== GC Content Computation Tests (SONIC-compatible) ===\n");

	test_basic_gc();
	test_multiple_windows();
	test_n_bases();
	test_truncation();
	test_sonic_window_boundaries();
	test_tile_indexing();
	test_case_insensitivity();
	test_real_sequence();
	test_chromosome_scale();
	test_consistency_with_sonic();

	printf("\n============================================\n");
	printf("  PASS: %d, FAIL: %d\n", tests_passed, tests_failed);
	printf("============================================\n");

	if(tests_failed > 0)
	{
		printf("GC TESTS FAILED\n");
		return 1;
	}
	printf("ALL GC TESTS PASSED\n");
	return 0;
}
