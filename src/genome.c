#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "genome.h"

extern long long memUsage;

static void* genome_getMem(size_t size)
{
	void *ret = malloc(size);
	if(ret == NULL)
	{
		fprintf(stderr, "Cannot allocate memory (%zu bytes)\n", size);
		exit(EXIT_FAILURE);
	}
	memUsage += size;
	return ret;
}

genome_info* genome_load_from_bam(bam_hdr_t *bam_header)
{
	int i;

	genome_info *genome = (genome_info*) genome_getMem(sizeof(genome_info));
	genome->number_of_chromosomes = bam_header->n_targets;
	genome->chromosome_names = (char**) genome_getMem(sizeof(char*) * genome->number_of_chromosomes);
	genome->chromosome_lengths = (int*) genome_getMem(sizeof(int) * genome->number_of_chromosomes);
	genome->chromosome_gc_profile = (char**) genome_getMem(sizeof(char*) * genome->number_of_chromosomes);
	genome->number_of_repeats_in_chromosome = (int*) genome_getMem(sizeof(int) * genome->number_of_chromosomes);
	genome->reps = (repeat_interval**) genome_getMem(sizeof(repeat_interval*) * genome->number_of_chromosomes);

	for(i = 0; i < genome->number_of_chromosomes; i++)
	{
		genome->chromosome_names[i] = (char*) genome_getMem(strlen(bam_header->target_name[i]) + 1);
		strcpy(genome->chromosome_names[i], bam_header->target_name[i]);
		genome->chromosome_lengths[i] = bam_header->target_len[i];
		genome->chromosome_gc_profile[i] = NULL;
		genome->number_of_repeats_in_chromosome[i] = 0;
		genome->reps[i] = NULL;
	}

	genome->last_chromosome_index = -1;

	return genome;
}

/* GC profile: 100bp non-overlapping windows, partial last window dropped */
void genome_compute_gc_profile(genome_info *genome, const char *ref_seq, int seq_len, int chr_index)
{
	int num_windows = (seq_len / GC_WINDOW) + 1;
	char *profile = (char*) calloc(num_windows, sizeof(char));
	int gc = 0;
	int char_count = 1;
	int window_id = 0;
	int i;

	for(i = 0; i < seq_len; i++)
	{
		char c = ref_seq[i];
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

	if(genome->chromosome_gc_profile[chr_index] != NULL)
		free(genome->chromosome_gc_profile[chr_index]);

	genome->chromosome_gc_profile[chr_index] = profile;
}

/* Average GC content over a genomic range */
float genome_get_gc_content(genome_info *genome, const char *chr_name, int pos_start, int pos_end)
{
	int chr_index = genome_find_chromosome_index(genome, chr_name);
	if(chr_index == -1 || genome->chromosome_gc_profile[chr_index] == NULL)
		return 0.0;

	char *profile = genome->chromosome_gc_profile[chr_index];
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

/* Check whether two intervals overlap */
static int interval_intersects(int pos_start, int pos_end, int start, int end)
{
	if(pos_start >= start && pos_end < end)
		return 1;
	if(pos_start <= start && pos_end > end)
		return 1;
	if(pos_start <= start && pos_end >= start)
		return 1;
	if(pos_start <= end && pos_end > end)
		return 1;
	return 0;
}

/* Fraction of interval1 that overlaps with interval2 */
static float intersection_fraction(int start_1, int end_1, int start_2, int end_2)
{
	if(start_1 >= start_2 && end_1 <= end_2)
		return 1.0;
	if(start_1 >= start_2 && end_1 >= end_2)
		return (float)(end_2 - start_1) / (float)(end_1 - start_1);
	if(start_1 <= start_2 && end_1 <= end_2)
		return (float)(end_1 - start_2) / (float)(end_1 - start_1);
	if(start_1 <= start_2 && end_1 <= end_2)
		return (float)(end_2 - start_2) / (float)(end_1 - start_1);
	return 0.0;
}

/* Binary search for an overlapping repeat interval */
static repeat_interval* genome_intersect(genome_info *genome, const char *chr_name, int pos_start, int pos_end)
{
	int start, end, med;
	int interval_count;
	repeat_interval *interval_list;
	int chr_index;

	chr_index = genome_find_chromosome_index(genome, chr_name);
	if(chr_index == -1)
		return NULL;

	interval_list = genome->reps[chr_index];
	interval_count = genome->number_of_repeats_in_chromosome[chr_index];

	if(interval_list == NULL || interval_count == 0)
		return NULL;

	start = 0;
	end = interval_count - 1;
	med = (start + end) / 2;

	while(1)
	{
		if(start > end)
			return NULL;

		if(interval_intersects(pos_start, pos_end, interval_list[med].start, interval_list[med].end))
			return &interval_list[med];

		if(start == med || end == med)
		{
			if(interval_intersects(pos_start, pos_end, interval_list[start].start, interval_list[start].end))
				return &interval_list[start];
			else if(interval_intersects(pos_start, pos_end, interval_list[end].start, interval_list[end].end))
				return &interval_list[end];
			return NULL;
		}

		else if(pos_start < interval_list[med].start)
		{
			end = med;
			med = (start + end) / 2;
		}

		else
		{
			start = med;
			med = (start + end) / 2;
		}
	}

	return NULL;
}

/* Return satellite overlap fraction for a query interval, 0 if not satellite */
float genome_is_satellite(genome_info *genome, const char *chr_name, int pos_start, int pos_end)
{
	repeat_interval *interval;
	char *is_satellite;

	interval = genome_intersect(genome, chr_name, pos_start, pos_end);

	if(interval == NULL)
		return 0;

	is_satellite = strstr(interval->repeat_item->repeat_class, "Satel");
	if(is_satellite != NULL)
		return intersection_fraction(pos_start, pos_end, interval->start, interval->end);

	is_satellite = strstr(interval->repeat_item->repeat_type, "Satel");
	if(is_satellite != NULL)
		return intersection_fraction(pos_start, pos_end, interval->start, interval->end);

	return 0;
}

static int compare_repeat(const void *a, const void *b)
{
	const repeat_interval *ra = (const repeat_interval*) a;
	const repeat_interval *rb = (const repeat_interval*) b;

	if(ra->start != rb->start)
		return ra->start - rb->start;
	return ra->end - rb->end;
}

/* Load repeat regions from a tab-delimited file: chr start end type class */
void genome_load_repeats(genome_info *genome, const char *reps_file)
{
	FILE *fp;
	char line[1024];
	char *chr_name, *start_str, *end_str, *rtype, *rclass;
	int i, chr_index;

	fp = fopen(reps_file, "r");
	if(fp == NULL)
	{
		fprintf(stderr, "[CONGA ERROR] Cannot open repeats file: %s\n", reps_file);
		return;
	}

	/* First pass: count entries per chromosome */
	while(fgets(line, 1024, fp) != NULL)
	{
		if(line[0] == '#' || line[0] == '\n')
			continue;

		chr_name = strtok(line, "\t\r\n ");
		if(chr_name == NULL) continue;

		chr_index = genome_find_chromosome_index(genome, chr_name);
		if(chr_index == -1) continue;

		genome->number_of_repeats_in_chromosome[chr_index]++;
	}

	/* Allocate arrays */
	for(i = 0; i < genome->number_of_chromosomes; i++)
	{
		if(genome->number_of_repeats_in_chromosome[i] > 0)
		{
			genome->reps[i] = (repeat_interval*) genome_getMem(
				sizeof(repeat_interval) * genome->number_of_repeats_in_chromosome[i]);
			genome->number_of_repeats_in_chromosome[i] = 0;
		}
	}

	/* Second pass: load entries */
	rewind(fp);
	while(fgets(line, 1024, fp) != NULL)
	{
		if(line[0] == '#' || line[0] == '\n')
			continue;

		chr_name = strtok(line, "\t\r\n ");
		if(chr_name == NULL) continue;

		start_str = strtok(NULL, "\t\r\n ");
		if(start_str == NULL) continue;

		end_str = strtok(NULL, "\t\r\n ");
		if(end_str == NULL) continue;

		rtype = strtok(NULL, "\t\r\n ");
		if(rtype == NULL) continue;

		rclass = strtok(NULL, "\t\r\n ");
		if(rclass == NULL) continue;

		chr_index = genome_find_chromosome_index(genome, chr_name);
		if(chr_index == -1) continue;

		int idx = genome->number_of_repeats_in_chromosome[chr_index];
		genome->reps[chr_index][idx].start = atoi(start_str);
		genome->reps[chr_index][idx].end = atoi(end_str);

		genome->reps[chr_index][idx].repeat_item = (repeat_info*) genome_getMem(sizeof(repeat_info));
		genome->reps[chr_index][idx].repeat_item->repeat_type = (char*) genome_getMem(strlen(rtype) + 1);
		strcpy(genome->reps[chr_index][idx].repeat_item->repeat_type, rtype);
		genome->reps[chr_index][idx].repeat_item->repeat_class = (char*) genome_getMem(strlen(rclass) + 1);
		strcpy(genome->reps[chr_index][idx].repeat_item->repeat_class, rclass);

		genome->number_of_repeats_in_chromosome[chr_index]++;
	}

	fclose(fp);

	/* Sort per-chromosome arrays by position for binary search */
	for(i = 0; i < genome->number_of_chromosomes; i++)
	{
		if(genome->number_of_repeats_in_chromosome[i] > 0)
		{
			qsort(genome->reps[i], genome->number_of_repeats_in_chromosome[i],
				sizeof(repeat_interval), compare_repeat);
		}
	}

	int total = 0;
	for(i = 0; i < genome->number_of_chromosomes; i++)
		total += genome->number_of_repeats_in_chromosome[i];
	fprintf(stderr, "[CONGA INFO] Loaded %d repeat regions from %s\n", total, reps_file);
}

/* Find chromosome index by name, with cache for sequential access */
int genome_find_chromosome_index(genome_info *genome, const char *chr_name)
{
	int i;

	if(genome->last_chromosome_index != -1 &&
	   strcmp(genome->chromosome_names[genome->last_chromosome_index], chr_name) == 0)
		return genome->last_chromosome_index;

	for(i = 0; i < genome->number_of_chromosomes; i++)
	{
		if(strcmp(genome->chromosome_names[i], chr_name) == 0)
		{
			genome->last_chromosome_index = i;
			return i;
		}
	}
	return -1;
}

void genome_free(genome_info *genome)
{
	int i, j;
	if(genome == NULL) return;

	for(i = 0; i < genome->number_of_chromosomes; i++)
	{
		if(genome->chromosome_names[i] != NULL)
			free(genome->chromosome_names[i]);
		if(genome->chromosome_gc_profile[i] != NULL)
			free(genome->chromosome_gc_profile[i]);
		if(genome->reps[i] != NULL)
		{
			for(j = 0; j < genome->number_of_repeats_in_chromosome[i]; j++)
			{
				if(genome->reps[i][j].repeat_item != NULL)
				{
					free(genome->reps[i][j].repeat_item->repeat_type);
					free(genome->reps[i][j].repeat_item->repeat_class);
					free(genome->reps[i][j].repeat_item);
				}
			}
			free(genome->reps[i]);
		}
	}

	free(genome->chromosome_names);
	free(genome->chromosome_lengths);
	free(genome->chromosome_gc_profile);
	free(genome->number_of_repeats_in_chromosome);
	free(genome->reps);
	free(genome);
}
