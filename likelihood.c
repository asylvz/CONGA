#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "likelihood.h"
#include "free.h"
#include "kmer.h"

#define hashtable_size 10000
int sv_count;
int del_count;
int dup_count;

int check_rp_intersection(svs arr[], int locus, int start_array, int end_array)
{
	if(end_array >= start_array)
	{
		int mid = start_array + (end_array - start_array) / 2;
		//fprintf(stderr,"Start:%d - End:%d - Mid:%d - loc-start:%d - loc-end:%d, Split: %d\n",start_array, end_array, mid, arr[mid].start, arr[mid].end, locus);

		if(locus > arr[mid].start - WRONGMAP_WINDOW && locus < arr[mid].end + WRONGMAP_WINDOW)
		{
			return mid;
		}
		else if(locus > arr[mid].start)
		{
			return check_rp_intersection(arr, locus, mid + 1, end_array);
		}
		else if(locus < arr[mid].start)
		{
			return check_rp_intersection(arr, locus, start_array, mid - 1);
		}
	}
	return -1;
}


void count_ReadPairs()
{
	SplitRow *splitRowPtr;
	SplitsInfo *aPtr = all_split_reads;

	int split_del_cnt = 0, split_dup_cnt = 0;
	int del_present = 0, del_border_present = 0, dup_present = 0, all_rp = 0, del_left, del_right, dup_left, dup_right, del_left_out, del_right_out;
	int del_start, del_end, dup_start, dup_end, i;

	int del_interval[50000];
	int dup_interval[50000];

	if(aPtr != NULL)
	{
		splitRowPtr = aPtr->head;
		while (splitRowPtr != NULL)
		{
			if(splitRowPtr->svType == DUPLICATION)
			{
				for(i = 0; i < dup_count; i++)
				{
					dup_start = all_svs_dup[i].start;
					dup_end = all_svs_dup[i].end;

					if(splitRowPtr->locMapLeftEnd >= (dup_start - WRONGMAP_WINDOW_DEL) && splitRowPtr->locMapLeftEnd <= (dup_end + WRONGMAP_WINDOW_DEL)
							&& splitRowPtr->locMapRightStart <= (dup_end + WRONGMAP_WINDOW_DEL) && splitRowPtr->locMapRightStart >= (dup_start - WRONGMAP_WINDOW_DEL))
					{
						dup_present++;
						all_svs_dup[i].rp++;
					}
				}
			}
			else if(splitRowPtr->svType == DELETION)
			{
				for(i = 0; i < del_count; i++)
				{
					del_start = all_svs_del[i].start;
					del_end = all_svs_del[i].end;

					if(splitRowPtr->locMapLeftEnd <= (del_start + WRONGMAP_WINDOW) && splitRowPtr->locMapLeftEnd >= (del_start - WRONGMAP_WINDOW_DEL)
							&& splitRowPtr->locMapRightStart >= (del_end - WRONGMAP_WINDOW) && splitRowPtr->locMapRightStart <= (del_end + WRONGMAP_WINDOW_DEL))
					{
						del_border_present++;
						all_svs_del[i].border_rp++;
					}
				}
			}
			else
				fprintf(stderr,"ERROR %c", splitRowPtr->svType);
			splitRowPtr = splitRowPtr->next;
		}
	}
	fprintf(stderr,"%d DELS and %d DUPS overlap with a known loci using %d bps wrong-map window\n", del_border_present, dup_present, WRONGMAP_WINDOW_DEL);
}

double lpoisson(int observed, double lambda)
{
	double penalty;

	penalty = 0.01;
	if  (lambda == 0.0)
		lambda = penalty;
	//fprintf(stderr,"%d - %f - %f - %f - %f - %f\n", observed, log(lambda), observed * log(lambda), lambda, lgamma(observed + 1), observed * log(lambda) - lambda - lgamma(observed + 1));
	return observed * log(lambda) - lambda - lgamma(observed + 1);
}


void calculate_likelihood_CNV(bam_info *in_bam, parameters *params, svs arr[], int count, char* chr_name, char type)
{
	int gc_val, i;
	double lhomo, lhete, lnone, score;
	float expected_rd_filtered = 0, expected_rd_unfiltered = 0, expected_kmer = 0.0, lambda;
	int observed_rd_filtered = 0, totalReadCount_kmer = 0, observed_rd_unfiltered = 0;
	double mappability_score = 0;


	for( i = arr[count].start; i < arr[count].end; i++)
	{
		gc_val = ( int)round ( sonic_get_gc_content(params->this_sonic, arr[count].chr_name, i, i + WINDOWSLIDE));

		expected_rd_unfiltered += in_bam->expected_rd_unfiltered[gc_val];
		observed_rd_unfiltered += in_bam->rd_unfiltered[i];

		if(params->mappability_file != NULL)
			mappability_score += in_bam->mappability[i];
	}
	//fprintf(stderr,"%d - %f ----- %d - %f\n", observed_rd_filtered, expected_rd_filtered, observed_rd_unfiltered, expected_rd_unfiltered);

	if(params->mappability_file != NULL)
		arr[count].mappability = mappability_score / (double) (arr[count].end - arr[count].start);

	if(!params->no_kmer)
	{
		for(i = arr[count].start; i < arr[count].end; i += KMERWINDOWSLIDE)
		{
			gc_val = (int)round ( sonic_get_gc_content(params->this_sonic, chr_name, i, i + KMERWINDOWSIZE));
			expected_kmer += in_bam->expected_kmer[gc_val];
			totalReadCount_kmer += in_bam->kmer[i];
		}
		arr[count].k_mer = totalReadCount_kmer;
		arr[count].expected_kmer = expected_kmer;

		if(type == DELETION)
		{
			lhomo = lpoisson(totalReadCount_kmer, 0.0);
			lhete = lpoisson(totalReadCount_kmer, 0.5 * expected_kmer);
			lnone = lpoisson(totalReadCount_kmer, expected_kmer);

			all_svs_del[count].likelihood_kmer = max(lhomo, lhete) / lnone;
		}
		else
		{
			lhomo = lpoisson(totalReadCount_kmer, 2 * expected_kmer);
			lhete = lpoisson(totalReadCount_kmer, 1.5 * expected_kmer);
			lnone = lpoisson(totalReadCount_kmer, expected_kmer);

			all_svs_dup[count].likelihood_kmer = max(lhomo, lhete) / lnone;
		}
		//fprintf(stderr,"%d - %lf - %lf\n", totalReadCount_kmer, expected_kmer, arr[count].likelihood_kmer);
	}

	if(type == DELETION)
	{
		lhomo = lpoisson(observed_rd_unfiltered, 0.0);
		lhete = lpoisson(observed_rd_unfiltered, 0.5 * expected_rd_unfiltered);
		lnone = lpoisson(observed_rd_unfiltered, expected_rd_unfiltered);

		score = max(lhomo, lhete) / lnone;
		fprintf(stderr,"%lf - %lf - %lf - %lf\n", lhomo, lhete, lnone, score);

		arr[count].likelihood_unfiltered = score;
		arr[count].observed_rd_all = observed_rd_unfiltered;
		arr[count].expected_rd_all = expected_rd_unfiltered;


		//fprintf(stderr,"%s\t%d\t%d\t%.2lf\t%.2f\t%f\t%f\t%f\t%ld\t%.2f\n", all_svs_del[count].chr_name, all_svs_del[count].start, all_svs_del[count].end, all_svs_del[count].del_likelihood, all_svs_del[count].copy_number,lhomo, lhete,lnone, arr[count].depth, expectedReadCount);

	}
	else if(type == DUPLICATION)
	{
		lhomo = lpoisson(observed_rd_unfiltered, 2 * expected_rd_unfiltered);
		lhete = lpoisson(observed_rd_unfiltered, 1.5 * expected_rd_unfiltered);
		lnone = lpoisson(observed_rd_unfiltered, expected_rd_unfiltered);

		score = max(lhomo, lhete) / lnone;

		arr[count].likelihood_unfiltered = score;
		arr[count].observed_rd_all = observed_rd_unfiltered;
		arr[count].expected_rd_all = expected_rd_unfiltered;
	}
}

void calculate_expected_CN( bam_info *in_bam, parameters *params, svs arr[], int count, char type)
{
	int pos, gc_val, i;
	float totalReadCount, expectedReadCount;

	totalReadCount = 0;
	expectedReadCount = 0;


	for( i = arr[count].start; i < arr[count].end; i++)
	{
		gc_val = ( int)round ( sonic_get_gc_content(params->this_sonic, arr[count].chr_name, i, i + WINDOWSLIDE));
		expectedReadCount += in_bam->expected_rd_unfiltered[gc_val];

		totalReadCount += ( float)in_bam->rd_unfiltered[i];
	}
	//if(type == DUPLICATION)
		//fprintf(stderr,"%lf - %lf\n", totalReadCount, expectedReadCount);

	arr[count].copy_number = ( float)( 2 * totalReadCount) / ( float)( expectedReadCount);
}


void output_SVs( parameters *params, FILE* fpSVs, FILE* fp_del, FILE* fp_dup)
{
	int count;
	int sv_cnt_dup = 0;
	int sv_cnt_del = 0;

	//qsort( all_svs, sv_count, sizeof( int), compare_start_pos);
	for( count = 0; count < del_count; count++)
	{
		fprintf(fp_del,"%s\t%d\t%d\t%.1f\t%.2f\t%d\t%.2lf\n", all_svs_del[count].chr_name, all_svs_del[count].start, all_svs_del[count].end, all_svs_del[count].copy_number, all_svs_del[count].likelihood_unfiltered, all_svs_del[count].border_rp, all_svs_del[count].mappability);

		if(params->mappability_file != NULL)
		{
			if(all_svs_del[count].likelihood_unfiltered < 0.5 && all_svs_del[count].mappability > 0.5)
			{
				fprintf(fpSVs,"%s\t%d\t%d\tDEL\t%.1f\t%.2f\t%d\t%.2lf\n", all_svs_del[count].chr_name, all_svs_del[count].start, all_svs_del[count].end, all_svs_del[count].copy_number, all_svs_del[count].likelihood_unfiltered, all_svs_del[count].border_rp, all_svs_del[count].mappability);
				sv_cnt_del++;
			}
		}
		else
		{
			if(all_svs_del[count].likelihood_unfiltered < 0.5)
			{
				fprintf(fpSVs,"%s\t%d\t%d\tDEL\t%.1f\t%.2f\t%d\n", all_svs_del[count].chr_name, all_svs_del[count].start, all_svs_del[count].end, all_svs_del[count].copy_number, all_svs_del[count].likelihood_unfiltered, all_svs_del[count].border_rp);
				sv_cnt_del++;
			}
		}
	}

	for( count = 0; count < dup_count; count++)
	{
		fprintf(fp_dup,"%s\t%d\t%d\t%.1f\t%.2lf\t%d\t%.2lf\n", all_svs_dup[count].chr_name, all_svs_dup[count].start, all_svs_dup[count].end, all_svs_dup[count].copy_number, all_svs_dup[count].likelihood_unfiltered, all_svs_dup[count].rp, all_svs_dup[count].mappability);

		if(!params->no_sr)
		{
			if(all_svs_dup[count].likelihood_unfiltered < 100 && all_svs_dup[count].rp > 10)
			{
				fprintf(fpSVs,"%s\t%d\t%d\tDUP\t%.1f\t%.2lf\t%d\t%.2lf\n", all_svs_dup[count].chr_name, all_svs_dup[count].start, all_svs_dup[count].end, all_svs_dup[count].copy_number, all_svs_dup[count].likelihood_unfiltered, all_svs_dup[count].rp, all_svs_dup[count].mappability);
				sv_cnt_dup++;
			}
		}
		else
		{
			if(all_svs_dup[count].likelihood_unfiltered < 1)
			{
				fprintf(fpSVs,"%s\t%d\t%d\tDUP\t%.2f\t%.2lf\t%d\t%.2lf\n", all_svs_dup[count].chr_name, all_svs_dup[count].start, all_svs_dup[count].end, all_svs_dup[count].copy_number, all_svs_dup[count].likelihood_unfiltered, all_svs_dup[count].rp, all_svs_dup[count].mappability);
				sv_cnt_dup++;
			}
		}
	}

	fprintf(stderr,"\nFound %d DELs - %d DUPs\n\n",sv_cnt_del, sv_cnt_dup);

	total_dels += sv_cnt_del;
	total_dups += sv_cnt_dup;
}

void find_depths( bam_info *in_bam, parameters *params, char* chr_name, int chr_index)
{
	int count;

	for( count = 0; count < del_count; count++)
	{
		calculate_expected_CN( in_bam, params, all_svs_del, count, DELETION);
		calculate_likelihood_CNV( in_bam, params, all_svs_del, count, chr_name, DELETION);
	}
	for( count = 0; count < dup_count; count++)
	{
		calculate_expected_CN( in_bam, params, all_svs_dup, count, DUPLICATION);
		calculate_likelihood_CNV( in_bam, params, all_svs_dup, count, chr_name, DUPLICATION);
	}
}


void find_SVs( bam_info *in_bam, parameters *params, FILE* fp_del, FILE* fp_dup, FILE* fp_SVs, char* chr_name, int chr_index)
{
	sv_count = 0;
	del_count = 0;
	dup_count = 0;
	int kmer_hash_size = 0;
	long total_kmers = 0;

	fprintf(stderr,"\nLoading known SVs ");
	load_known_SVs( &all_svs_del, &all_svs_dup, params, chr_name, &del_count, &dup_count);
	fprintf( stderr, "(%d DELS, %d DUPS in chromosome %s - larger than the threshold %d)\n", del_count, dup_count, chr_name, params->min_sv_size);

	//Sort the known SVs
	qsort( all_svs_del, del_count, sizeof(svs), compare_start_pos);
	qsort( all_svs_dup, dup_count, sizeof(svs), compare_start_pos);

	sv_count = del_count + dup_count;

	if(sv_count == 0)
	{
		free(in_bam->rd_unfiltered);
		return;
	}

	if(params->low_map_regions != NULL)
	{
		//fprintf(stderr,"Loading Low mappability regions\n");
		check_low_mappability(params, all_svs_del, all_svs_dup, chr_name, del_count, dup_count);
	}

	if(!params->no_sr)
	{
		count_ReadPairs();
		free_splits(in_bam);
	}

	/*if(!params->no_kmer)
	{
		// Read the fastq file
		fprintf(stderr,"\nReading K-MERS\n");
		kmer_hash_size = read_kmer_jellyfish(params);

		//fprintf(stderr,"\nHash size %d\n", kmer_hash_size);
		if(kmer_hash_size < MINKMERHASHSIZE)
		{
			free(in_bam->rd_unfiltered);
			return;
		}

		init_kmer_per_chr(in_bam, params, chr_index);

		fprintf(stderr,"-->calculating k-mer counts");
		total_kmers = calc_kmer_counts(in_bam, params, chr_index);
		fprintf(stderr," (%li kmers)\n", total_kmers);

		free_hash_table_kmer(params);

		fprintf(stderr,"-->calculating expected counts\n");
		calc_expected_kmer(in_bam, params, chr_index);

	}*/

	//Check mappability
	fprintf(stderr,"Finding mappability for each region\n");
	if(params->mappability_file != NULL)
	{
		init_mappability_per_chr( in_bam, params, chr_index);
		load_mappability_regions( in_bam, params, chr_name);
	}

	fprintf(stderr,"\nCalculating Likelihoods\n");
	find_depths(in_bam, params, chr_name, chr_index);

	free(in_bam->rd_unfiltered);

	if(params->mappability_file != NULL)
		free(in_bam->mappability);

	/*if(!params->no_kmer)
	{
		free(in_bam->kmer);
		in_bam->kmer = NULL;
	}*/
	//fprintf(stderr,"Outputting\n");
	output_SVs(params, fp_SVs, fp_del, fp_dup);

	free_SVs(all_svs_del, del_count);
	free_SVs(all_svs_dup, dup_count);
}
