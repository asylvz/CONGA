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
	int del_start, del_end, i;

	int del_interval[50000];

	if(aPtr != NULL)
	{
		splitRowPtr = aPtr->head;
		while (splitRowPtr != NULL)
		{

			del_left = check_rp_intersection(all_svs_del, splitRowPtr->locMapLeftStart, 0, del_count);
			del_right = check_rp_intersection(all_svs_del, splitRowPtr->locMapRightStart, 0, del_count);


			dup_left = check_rp_intersection(all_svs_dup, splitRowPtr->locMapLeftStart, 0, dup_count);
			dup_right = check_rp_intersection(all_svs_dup, splitRowPtr->locMapRightStart, 0, dup_count);

			if(splitRowPtr->svType == DELETION)
			{
				split_del_cnt++;

				if((del_left != -1) && (del_right == del_left))
				{
					del_start = all_svs_del[del_left].start;
					del_end = all_svs_del[del_left].end;

					del_interval[del_present] = del_left;

					all_svs_del[del_left].rp++;
					del_present++;
				}
			}
			else if(splitRowPtr->svType == DUPLICATION)
			{
				split_dup_cnt++;

				if((dup_left != -1) && (dup_left == dup_right))
				{
					all_svs_dup[dup_left].rp++;
					dup_present++;
				}
			}
			else
				fprintf(stderr,"\nERRORRRR\n");

			splitRowPtr = splitRowPtr->next;
		}
	}

	splitRowPtr = aPtr->head;
	while (splitRowPtr != NULL)
	{
		for(i = 0; i < del_present; i++)
		{
			del_left = del_interval[i];
			del_start = all_svs_del[del_left].start;
			del_end = all_svs_del[del_left].end;

			if(splitRowPtr->locMapLeftEnd <= del_start && splitRowPtr->locMapLeftEnd >= (del_start - WRONGMAP_WINDOW_DEL)
					&& splitRowPtr->locMapRightStart >= del_end && splitRowPtr->locMapRightStart <= (del_end + WRONGMAP_WINDOW_DEL))
			{
				del_border_present++;
				all_svs_del[del_left].border_rp++;
			}
		}
		splitRowPtr = splitRowPtr->next;
	}
	fprintf(stderr,"%d DELS (%d around breakpoints) and %d DUPS overlap with a known loci using %d bps wrong-map window\n", del_present, del_border_present, dup_present, WRONGMAP_WINDOW);
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


void calculate_likelihood_CNV(bam_info *in_bam, parameters *params, svs arr[], int count, float *mean_kmer, char* chr_name, char type)
{
	int gc_val, i;
	double lhomo, lhete, lnone, score, expected_kmer;
	float expectedReadCount = 0, lambda;
	int totalReadCount = 0, totalReadCount_kmer;


	if(arr[count].low_mappability == true)
	{
			arr[count].likelihood = 9999;
			return;
	}

	for( i = arr[count].start; i < arr[count].end; i++)
	{
		gc_val = ( int)round ( sonic_get_gc_content(params->this_sonic, arr[count].chr_name, i, i + WINDOWSLIDE));
		expectedReadCount += in_bam->mean_rd_per_gc[gc_val];
		totalReadCount += in_bam->read_depth_per_chr[i];
	}
	//arr[count].depth = (long) totalReadCount;

	if(!params->no_kmer)
	{
		totalReadCount_kmer = 0;
		expected_kmer = 0;
		for(i = arr[count].start; i < arr[count].end; i += KMERWINDOWSLIDE)
		{
			gc_val = (int)round ( sonic_get_gc_content(params->this_sonic, chr_name, i, i + KMERWINDOWSLIDE));
			expected_kmer += mean_kmer[gc_val];
			totalReadCount_kmer += calculate_kmers(params, chr_name, i, i + KMERWINDOWSLIDE);
		}
		arr[count].k_mer = totalReadCount_kmer;
		arr[count].expected_kmer = expected_kmer;

		lhomo = lpoisson(totalReadCount_kmer, 0.0);
		lhete = lpoisson(totalReadCount_kmer, 0.5 * expected_kmer);
		lnone = lpoisson(totalReadCount_kmer, expected_kmer);

		if(type == DELETION)
			all_svs_del[count].likelihood_kmer = max(lhomo, lhete) / lnone;
		else
			all_svs_dup[count].likelihood_kmer = max(lhomo, lhete) / lnone;

		//fprintf(stderr,"%d - %lf - %lf\n", totalReadCount_kmer, expected_kmer, arr[count].likelihood_kmer);
	}

	if(type == DELETION)
	{

		lhomo = lpoisson(totalReadCount, 0.0);
		lhete = lpoisson(totalReadCount, 0.5 * expectedReadCount);
		lnone = lpoisson(totalReadCount, expectedReadCount);

		score = max(lhomo, lhete) / lnone;
		arr[count].likelihood = score;

		//fprintf(stderr,"%s\t%d\t%d\t%.2lf\t%.2f\t%f\t%f\t%f\t%ld\t%.2f\n", all_svs_del[count].chr_name, all_svs_del[count].start, all_svs_del[count].end, all_svs_del[count].del_likelihood, all_svs_del[count].copy_number,lhomo, lhete,lnone, arr[count].depth, expectedReadCount);

	}
	else if(type == DUPLICATION)
	{
		lhomo = lpoisson(totalReadCount, 2 * expectedReadCount);
		lhete = lpoisson(totalReadCount, 1.5 * expectedReadCount);
		lnone = lpoisson(totalReadCount, expectedReadCount);

		score = max(lhomo, lhete) / lnone;
		arr[count].likelihood = score;
	}
}

void calculate_expected_CN( bam_info *in_bam, parameters *params, svs arr[], int count)
{
	int pos, gc_val, i;
	float totalReadCount, expectedReadCount;

	totalReadCount = 0;
	expectedReadCount = 0;

	for( i = arr[count].start; i < arr[count].end; i++)
	{
		gc_val = ( int)round ( sonic_get_gc_content(params->this_sonic, arr[count].chr_name, i, i + WINDOWSLIDE));
		expectedReadCount += in_bam->mean_rd_per_gc[gc_val];

		totalReadCount += ( float)in_bam->read_depth_per_chr[i];

	}
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
		//if(all_svs_del[count].del_likelihood > params->rd_threshold && (all_svs_del[count].copy_number < 0.3 || all_svs_del[count].rp >= params->rp_support))
		//if(all_svs_del[count].del_likelihood > params->rd_threshold && all_svs_del[count].rp < 10 && all_svs_del[count].copy_number <= 0.4)
		fprintf(fp_del,"%s\t%d\t%d\t%.2lf\t%.2f\t%d\t%d\t%d\t%f\t%f\n", all_svs_del[count].chr_name, all_svs_del[count].start, all_svs_del[count].end, all_svs_del[count].likelihood, all_svs_del[count].copy_number, all_svs_del[count].rp, all_svs_del[count].border_rp, all_svs_del[count].k_mer, all_svs_del[count].likelihood_kmer, all_svs_del[count].expected_kmer);
		if(all_svs_del[count].likelihood < 0.5 && all_svs_del[count].likelihood_kmer < 100)
		{
			fprintf(fpSVs,"%s\t%d\t%d\tDEL\t%.2lf\t%.1f\t%d\t%d\t%f\n", all_svs_del[count].chr_name, all_svs_del[count].start, all_svs_del[count].end, all_svs_del[count].likelihood, all_svs_del[count].copy_number, all_svs_del[count].border_rp, all_svs_del[count].k_mer, all_svs_del[count].likelihood_kmer);
			sv_cnt_del++;
		}
	}

	for( count = 0; count < dup_count; count++)
	{
		//if(all_svs_dup[count].dup_likelihood > (params->rd_threshold / 10) && all_svs_dup[count].rp >= params->rp_support)
		fprintf(fp_dup,"%s\t%d\t%d\t%.2lf\t%.2f\t%d\t%d\t%f\n", all_svs_dup[count].chr_name, all_svs_dup[count].start, all_svs_dup[count].end, all_svs_dup[count].likelihood, all_svs_dup[count].copy_number, all_svs_dup[count].rp, all_svs_dup[count].k_mer, all_svs_dup[count].likelihood_kmer);

		if(all_svs_dup[count].likelihood < 0.5)
		{
			fprintf(fpSVs,"%s\t%d\t%d\tDUP\t%.2lf\t%.1f\t%d\t%d\t%f\n", all_svs_dup[count].chr_name, all_svs_dup[count].start, all_svs_dup[count].end, all_svs_dup[count].likelihood, all_svs_dup[count].copy_number, all_svs_dup[count].rp, all_svs_dup[count].k_mer, all_svs_dup[count].likelihood_kmer);
			sv_cnt_dup++;
		}
	}

	fprintf(stderr,"\nFound %d DELs - %d DUPs\n\n",sv_cnt_del, sv_cnt_dup);

	total_dels += sv_cnt_del;
	total_dups += sv_cnt_dup;
}

void find_depths( bam_info *in_bam, parameters *params, char* chr_name, int chr_index)
{
	int count;
	float mean_kmer[101];

	if(!params->no_kmer)
		calc_mean_kmer(params, chr_name, chr_index, &mean_kmer[0]);

	for( count = 0; count < del_count; count++)
	{
		calculate_expected_CN( in_bam, params, all_svs_del, count);
		calculate_likelihood_CNV( in_bam, params, all_svs_del, count,  &mean_kmer[0], chr_name, DELETION);
	}
	for( count = 0; count < dup_count; count++)
	{
		calculate_expected_CN( in_bam, params, all_svs_dup, count);
		calculate_likelihood_CNV( in_bam, params, all_svs_dup, count, &mean_kmer[0], chr_name, DUPLICATION);
	}
}


void count_kmers_find_likelihoods(parameters *params, char* chr_name, int chr_index)
{
	int count, var_size, i;
	int return_value = -1;
	int max_size = params->this_sonic->chromosome_lengths[chr_index] - 1001;
	int k_mer_count = 0, rand_loci;
	long sum_kmers = 0;
	double expected_kmer_per_base, lhomo, lhete, lnone, observed, expected_kmer;
	float mean_kmer[101];
	int gc_val;
	int totalReadCount = 0;

	//Calculate expected k_mer
	calc_mean_kmer(params, chr_name, chr_index, &mean_kmer[0]);

	for( count = 0; count < del_count; count ++)
	{
		totalReadCount = 0;
		expected_kmer = 0;

		for(i = all_svs_del[count].start; i < all_svs_del[count].end; i+= KMERWINDOWSIZE)
		{
			gc_val = ( int)round ( sonic_get_gc_content(params->this_sonic, chr_name, i, i + KMERWINDOWSLIDE));
			expected_kmer += mean_kmer[gc_val];
			totalReadCount += calculate_kmers(params, chr_name, i, i + KMERWINDOWSLIDE);
		}


		//k_mer_count = calculate_kmers(params, chr_name, all_svs_del[count].start, all_svs_del[count].end);
		all_svs_del[count].k_mer = totalReadCount;
		var_size = all_svs_del[count].end - all_svs_del[count].start;
		observed = (double) k_mer_count * ((double) KMERWINDOWSLIDE / var_size);

		lhomo = lpoisson(totalReadCount, 0.0);
		lhete = lpoisson(totalReadCount, 0.5 * expected_kmer);
		lnone = lpoisson(totalReadCount, expected_kmer);

		all_svs_del[count].likelihood_kmer = max(lhomo, lhete) / lnone;

		//fprintf(stderr,"%d - %lf - %lf\n", totalReadCount, expected_kmer, all_svs_del[count].del_likelihood_kmer);
	}

	for( count = 0; count < dup_count; count++)
	{
		k_mer_count = calculate_kmers(params, chr_name, all_svs_dup[count].start, all_svs_dup[count].end);
		all_svs_dup[count].k_mer = k_mer_count;
		var_size = all_svs_dup[count].end - all_svs_dup[count].start;
		observed = (double) k_mer_count * ((double) KMERWINDOWSLIDE / var_size);

		lhomo = lpoisson(observed, 2 * expected_kmer_per_base);
		lhete = lpoisson(observed, 1.5 * expected_kmer_per_base);
		lnone = lpoisson(observed, expected_kmer_per_base);

		all_svs_dup[count].likelihood_kmer = max(lhomo, lhete) / lnone;
	}
}


void find_SVs( bam_info *in_bam, parameters *params, FILE* fp_del, FILE* fp_dup, FILE* fp_SVs, char* chr_name, int chr_index)
{
	sv_count = 0;
	del_count = 0;
	dup_count = 0;

	fprintf(stderr,"\nLoading SVs\n");
	load_known_SVs( &all_svs_del, &all_svs_dup, params, chr_name, &del_count, &dup_count);
	fprintf( stderr, "%d DELS, %d DUPS in the input file for chromosome %s (larger than the threshold of %d)\n", del_count, dup_count, chr_name, params->min_sv_size);

	//Sort the known SVs
	qsort( all_svs_del, del_count, sizeof(svs), compare_start_pos);
	qsort( all_svs_dup, dup_count, sizeof(svs), compare_start_pos);

	sv_count = del_count + dup_count;

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

	if(!params->no_kmer)
	{
		// Read the fastq file

		fprintf(stderr,"Reading k-mer counts\n");
		read_kmer_jellyfish(params);
	}

	fprintf(stderr,"Finding depths\n");
	find_depths(in_bam, params, chr_name, chr_index);
	free( in_bam->read_depth_per_chr);
	if(!params->no_kmer)
		free_hash_table_kmer(params);

	//fprintf(stderr,"Counting k-mers (hash size=%ld)\n", params->hash_size_kmer);
	//count_kmers_find_likelihoods(params, chr_name, chr_index);
	//free_hash_table_kmer(params);

	//Filter Low Mappability Regions

	output_SVs(params, fp_SVs, fp_del, fp_dup);

	free_SVs(all_svs_del, del_count);
	free_SVs(all_svs_dup, dup_count);
}
