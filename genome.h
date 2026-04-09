#ifndef GENOME_H_
#define GENOME_H_

#include <htslib/sam.h>
#include <htslib/faidx.h>
#include <stdbool.h>

#define GC_WINDOW 100

typedef struct _repeat_info
{
	char *repeat_type;
	char *repeat_class;
} repeat_info;

typedef struct _repeat_interval
{
	int start;
	int end;
	repeat_info *repeat_item;
} repeat_interval;

typedef struct _genome_info
{
	int number_of_chromosomes;
	char **chromosome_names;
	int *chromosome_lengths;

	char **chromosome_gc_profile;

	int *number_of_repeats_in_chromosome;
	repeat_interval **reps;

	int last_chromosome_index; /* cache for genome_find_chromosome_index */
} genome_info;

/* Load chromosome names and lengths from BAM header */
genome_info* genome_load_from_bam(bam_hdr_t *bam_header);

/* Compute GC profile for a single chromosome from reference sequence */
void genome_compute_gc_profile(genome_info *genome, const char *ref_seq, int seq_len, int chr_index);

/* Get average GC content for a genomic range (replicates sonic_get_gc_content) */
float genome_get_gc_content(genome_info *genome, const char *chr_name, int pos_start, int pos_end);

/* Check if a position overlaps a satellite region (replicates sonic_is_satellite) */
float genome_is_satellite(genome_info *genome, const char *chr_name, int pos_start, int pos_end);

/* Load repeat regions from a BED file (chr start end repeat_type repeat_class) */
void genome_load_repeats(genome_info *genome, const char *reps_file);

/* Find chromosome index by name */
int genome_find_chromosome_index(genome_info *genome, const char *chr_name);

/* Free genome_info */
void genome_free(genome_info *genome);

#endif /* GENOME_H_ */
