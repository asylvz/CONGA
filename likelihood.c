#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "likelihood.h"
#include "free.h"

#define hashtable_size 10000
int sv_count;
int del_count;
int dup_count;
svs* all_svs_del;
svs* all_svs_dup;

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

	int del_interval[10000];

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

			if(splitRowPtr->locMapLeftEnd < del_start && splitRowPtr->locMapLeftEnd > del_start - WRONGMAP_WINDOW_DEL
					&& splitRowPtr->locMapRightStart > del_end && splitRowPtr->locMapRightStart < del_end + WRONGMAP_WINDOW_DEL)
			{
				del_border_present++;
				all_svs_del[del_left].border_rp++;
			}
		}
		splitRowPtr = splitRowPtr->next;
	}

	//fprintf(stderr,"%d DELS and %d DUPS found by split-reads\n", split_del_cnt, split_dup_cnt);
	fprintf(stderr,"%d DELS (%d around breakpoints) and %d DUPS overlap with a known loci using %d bps wrong-map window\n", del_present, del_border_present, dup_present, WRONGMAP_WINDOW);
}


void calculate_likelihood_CNV(bam_info *in_bam, parameters *params, svs arr[], int count)
{
	int gc_val, i;
	float expectedReadCount = 0, lambda;
	int totalReadCount = 0;

	for( i = arr[count].start; i < arr[count].end; i++)
	{
		gc_val = ( int)round ( sonic_get_gc_content(params->this_sonic, arr[count].chr_name, i, i + WINDOWSLIDE));
		expectedReadCount += in_bam->mean_rd_per_gc[gc_val];
		totalReadCount += in_bam->read_depth_per_chr[i];
	}
	arr[count].depth = (long) totalReadCount;
	for( i = 0; i < 10; i++)
	{
		if( i > 0)
		{
			lambda = (((float)i / ( float)2 ) * expectedReadCount);
			arr[count].cnv_probability[i] =  exp( totalReadCount * log( lambda) - lambda - lgamma( totalReadCount + 1));
		}
		else if( i == 0)
		{
			if(totalReadCount < (0.2) * expectedReadCount)
				arr[count].cnv_probability[i] = 1;
			else
			{
				lambda = ((0.01) * expectedReadCount);
				arr[count].cnv_probability[i] = exp( totalReadCount * log( lambda) - lambda - lgamma( totalReadCount + 1));
			}
		}
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


void output_SVs( parameters *params, FILE* fp_del, FILE* fp_dup, FILE* fpSVs, char* chr_name)
{
	int count;
	int sv_cnt_dup = 0;
	int sv_cnt_del = 0;
	double maxRDlikelihoodDel, maxRDlikelihoodDup;
	double sumLikelihoodDup, sumLikelihoodDel, sumLikelihoodNorm;
	double epsilon = 1e-5; /* was 1e-200 before */

	for( count = 0; count < del_count; count++)
	{
		sumLikelihoodNorm = all_svs_del[count].cnv_probability[2];
		sumLikelihoodDel = all_svs_del[count].cnv_probability[0] + all_svs_del[count].cnv_probability[1];
		sumLikelihoodDup = all_svs_del[count].cnv_probability[3] + all_svs_del[count].cnv_probability[4] +
				all_svs_del[count].cnv_probability[5] + all_svs_del[count].cnv_probability[6];

		all_svs_del[count].del_likelihood = sumLikelihoodDel / (double)( sumLikelihoodDup + sumLikelihoodNorm + epsilon);
		fprintf(fp_del,"%s\t%d\t%d\t%.2lf\t%.2f\t%d\t%d\n", all_svs_del[count].chr_name, all_svs_del[count].start, all_svs_del[count].end, all_svs_del[count].del_likelihood, all_svs_del[count].copy_number, all_svs_del[count].rp, all_svs_del[count].border_rp);
	}

	for( count = 0; count < dup_count; count++)
	{
		sumLikelihoodNorm = all_svs_dup[count].cnv_probability[2];
		sumLikelihoodDel = all_svs_dup[count].cnv_probability[0] + all_svs_dup[count].cnv_probability[1];
		sumLikelihoodDup = all_svs_dup[count].cnv_probability[3] + all_svs_dup[count].cnv_probability[4] +
				all_svs_dup[count].cnv_probability[5] + all_svs_dup[count].cnv_probability[6];

		all_svs_dup[count].dup_likelihood = sumLikelihoodDup / (double)( sumLikelihoodDel + sumLikelihoodNorm + epsilon);
		fprintf(fp_dup,"%s\t%d\t%d\t%.2lf\t%.2f\t%d\n", all_svs_dup[count].chr_name, all_svs_dup[count].start, all_svs_dup[count].end, all_svs_dup[count].dup_likelihood, all_svs_dup[count].copy_number, all_svs_dup[count].rp);
	}

	//qsort( all_svs, sv_count, sizeof( int), compare_start_pos);
	for( count = 0; count < del_count; count++)
	{
		//if(all_svs_del[count].del_likelihood > params->rd_threshold && (all_svs_del[count].copy_number < 0.3 || all_svs_del[count].rp >= params->rp_support))
		if(all_svs_del[count].del_likelihood > params->rd_threshold && all_svs_del[count].rp < 10 && all_svs_del[count].copy_number <= 0.4)
		{
			fprintf(fpSVs,"%s\t%d\t%d\tDEL\t%.2lf\t%.1f\t%d\n", all_svs_del[count].chr_name, all_svs_del[count].start, all_svs_del[count].end, all_svs_del[count].del_likelihood, all_svs_del[count].copy_number, all_svs_del[count].rp);
			sv_cnt_del++;
		}
	}

	for( count = 0; count < dup_count; count++)
	{
		if(all_svs_dup[count].dup_likelihood > (params->rd_threshold / 10) && all_svs_dup[count].rp >= params->rp_support)
		{
			fprintf(fpSVs,"%s\t%d\t%d\tDUP\t%.2lf\t%.1f\t%d\n", all_svs_dup[count].chr_name, all_svs_dup[count].start, all_svs_dup[count].end, all_svs_dup[count].dup_likelihood, all_svs_dup[count].copy_number, all_svs_dup[count].rp);
			sv_cnt_dup++;
		}
	}

	fprintf(stderr,"\nFound %d DELs - %d DUPs\n\n",sv_cnt_del, sv_cnt_dup);

	total_dels += sv_cnt_del;
	total_dups += sv_cnt_dup;
}

void find_depths( bam_info *in_bam, parameters *params)
{
	int count;

	for( count = 0; count < del_count; count++)
	{
		calculate_expected_CN( in_bam, params, all_svs_del, count);
		calculate_likelihood_CNV( in_bam, params, all_svs_del, count);
	}
	for( count = 0; count < dup_count; count++)
	{
		calculate_expected_CN( in_bam, params, all_svs_dup, count);
		calculate_likelihood_CNV( in_bam, params, all_svs_dup, count);
	}
}


void find_SVs( bam_info *in_bam, parameters *params, FILE* fp_del, FILE* fp_dup, FILE* fp_SVs, char* chr_name)
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

	count_ReadPairs();

	find_depths(in_bam, params);
	output_SVs(params, fp_del, fp_dup, fp_SVs, chr_name);

	free_SVs(all_svs_del,del_count);
	free_SVs(all_svs_dup,dup_count);
}
