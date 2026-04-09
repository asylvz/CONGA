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

	int last_chromosome_index;
} genome_info;

genome_info* genome_load_from_bam(bam_hdr_t *bam_header);
void genome_compute_gc_profile(genome_info *genome, const char *ref_seq, int seq_len, int chr_index);
float genome_get_gc_content(genome_info *genome, const char *chr_name, int pos_start, int pos_end);
float genome_is_satellite(genome_info *genome, const char *chr_name, int pos_start, int pos_end);
void genome_load_repeats(genome_info *genome, const char *reps_file);
int genome_find_chromosome_index(genome_info *genome, const char *chr_name);
void genome_free(genome_info *genome);

#endif /* GENOME_H_ */
