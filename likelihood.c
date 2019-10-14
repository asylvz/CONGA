#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "likelihood.h"

#define hashtable_size 10000
int sv_count;
int del_count;
int dup_count;
svs* all_svs;

void calculateLikelihoodCNV(bam_info *in_bam, parameters *params, int count)
{
	int gc_val, i;
	float expectedReadCount = 0, lambda;
	int totalReadCount = 0;


	for( i = all_svs[count].start; i < all_svs[count].end; i++)
	{
		gc_val = ( int)round ( sonic_get_gc_content(params->this_sonic, all_svs[count].chr_name, i, i + WINDOWSLIDE));
		expectedReadCount += in_bam->mean_rd_per_gc[gc_val];
		totalReadCount += in_bam->read_depth_per_chr[i];

	}
	all_svs[count].depth = ( long)totalReadCount;
	for( i = 0; i < 10; i++)
	{
		if( i > 0)
		{
			lambda = ( ( ( float)i / ( float)2 ) * expectedReadCount);
			all_svs[count].cnv_probability[i] =  exp( totalReadCount * log( lambda) - lambda - lgamma( totalReadCount + 1));
			//fprintf(stderr, "%lf\n", all_svs[count].cnv_probability[i]);
		}
		else if( i == 0)
		{
			if( totalReadCount < ( 0.2) * expectedReadCount)
				all_svs[count].cnv_probability[i] = 1;
			else
			{
				lambda = ( ( 0.01) * expectedReadCount);
				all_svs[count].cnv_probability[i] = exp( totalReadCount * log( lambda) - lambda - lgamma( totalReadCount + 1));
				//fprintf(stderr, "%lf\n", all_svs[count].cnv_probability[i]);
				//fprintf(stderr, "i = %d read cnt= %d expected=%f - CNV=%lf\n", i, totalReadCount, expectedReadCount, cnvProb[count][i]);
			}
		}
	}
}

void calculateExpectedCN( bam_info *in_bam, parameters *params, int count)
{
	int pos, gc_val, i;
	float totalReadCount, expectedReadCount;

	totalReadCount = 0;
	expectedReadCount = 0;

	for( i = all_svs[count].start; i < all_svs[count].end; i++)
	{
		gc_val = ( int)round ( sonic_get_gc_content(params->this_sonic, all_svs[count].chr_name, i, i + WINDOWSLIDE));
		expectedReadCount += in_bam->mean_rd_per_gc[gc_val];

		totalReadCount += ( float)in_bam->read_depth_per_chr[i];

	}
	all_svs[count].copy_number = ( float)( 2 * totalReadCount) / ( float)( expectedReadCount);
	//fprintf(stderr, "total=%f - expected=%f - CN =%f\n", totalReadCount, expectedReadCount, vars.copy_number);

}


void output_SVs( parameters *params, FILE* fp_del, FILE* fp_dup, FILE* fpSVs, char* chr_name)
{
	int count;
	int sv_cnt_dup = 0;
	int sv_cnt_del = 0;
	double maxRDlikelihoodDel, maxRDlikelihoodDup;
	double sumLikelihoodDup, sumLikelihoodDel, sumLikelihoodNorm;
	double epsilon = 1e-5; /* was 1e-200 before */

	for( count = 0; count < sv_count; count++)
	{
		sumLikelihoodNorm = all_svs[count].cnv_probability[2];
		sumLikelihoodDel = all_svs[count].cnv_probability[0] + all_svs[count].cnv_probability[1];
		sumLikelihoodDup = all_svs[count].cnv_probability[3] + all_svs[count].cnv_probability[4] +
				all_svs[count].cnv_probability[5] + all_svs[count].cnv_probability[6];

		if(all_svs[count].SV_type == DUPLICATION)
		{
			all_svs[count].dup_likelihood = sumLikelihoodDup / (double)( sumLikelihoodDel + sumLikelihoodNorm + epsilon);
			fprintf(fp_dup,"%s\t%d\t%d\t%.2lf\n", all_svs[count].chr_name, all_svs[count].start, all_svs[count].end, all_svs[count].dup_likelihood);
		}
		else
		{
			all_svs[count].del_likelihood = sumLikelihoodDel / (double)( sumLikelihoodDup + sumLikelihoodNorm + epsilon);
			fprintf(fp_del,"%s\t%d\t%d\t%.2lf\n", all_svs[count].chr_name, all_svs[count].start, all_svs[count].end, all_svs[count].del_likelihood);
		}
	}

	//Sort the SVs
	//qsort( all_svs, sv_count, sizeof( int), compare_start_pos);

	for( count = 0; count < sv_count; count++)
	{
		//fprintf(stderr,"%s\t%d\t%d\t%c\t%.2lf\t%.2lf\t%.2f\n", all_svs[count].chr_name, all_svs[count].start, all_svs[count].end, all_svs[count].SV_type, all_svs[count].del_likelihood, all_svs[count].dup_likelihood,all_svs[count].copy_number);
		if(all_svs[count].SV_type == DUPLICATION && all_svs[count].dup_likelihood > params->rd_threshold * 5 && all_svs[count].copy_number > 3)
		{
			fprintf(fpSVs,"%s\t%d\t%d\tDUP\t%.2lf\t%.1f\n", all_svs[count].chr_name, all_svs[count].start, all_svs[count].end, all_svs[count].dup_likelihood, all_svs[count].copy_number);
			sv_cnt_dup++;
		}
		else if(all_svs[count].SV_type == DELETION && all_svs[count].del_likelihood > params->rd_threshold && all_svs[count].copy_number < 0.3)
		{
			fprintf(fpSVs,"%s\t%d\t%d\tDEL\t%.2lf\t%.1f\n", all_svs[count].chr_name, all_svs[count].start, all_svs[count].end, all_svs[count].del_likelihood, all_svs[count].copy_number);
			sv_cnt_del++;
		}
	}
	fprintf(stderr,"  (%d DELs - %d DUPs)\n",sv_cnt_del, sv_cnt_dup);

	total_dels += sv_cnt_del;
	total_dups += sv_cnt_dup;

	//Free the SVs
	for( count = 0; count < sv_count; count++)
	{
		if(all_svs[count].chr_name != NULL)
			free(all_svs[count].chr_name);
	}
	free(all_svs);
}

void find_depths( bam_info *in_bam, parameters *params)
{
	int count;

	for( count = 0; count < sv_count; count++)
	{
		calculateExpectedCN( in_bam, params, count);
		calculateLikelihoodCNV( in_bam, params, count);
	}
}


void find_SVs( bam_info *in_bam, parameters *params, FILE* fp_del, FILE* fp_dup, FILE* fp_SVs, char* chr_name)
{
	sv_count = 0;
	del_count = 0;
	dup_count = 0;

	load_known_SVs( &all_svs, params, chr_name, &del_count, &dup_count);
	//fprintf( stderr, "\n\n%d DELS, %d DUPS in the input for chromosome %s\n\n", del_count, dup_count, chr_name);

	sv_count = del_count + dup_count;

	find_depths(in_bam, params);
	output_SVs(params, fp_del, fp_dup, fp_SVs, chr_name);

}
