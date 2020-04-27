#include "read_distribution.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <htslib/faidx.h>
#include <stdbool.h>

#define STR_SIZE 255


void init_rd_per_chr( bam_info* in_bam, parameters* param, int chr_index)
{
	in_bam->read_count = 0;
	in_bam->read_depth_per_chr = ( short*) getMem( sizeof( short) * ( param->this_sonic->chromosome_lengths[chr_index]));

	memset (in_bam->read_depth_per_chr, 0, (param->this_sonic->chromosome_lengths[chr_index] * sizeof(short)));
}


void calc_mu_per_chr( bam_info* in_bam, int chromosome_length)
{
	int i;
	long rd_cnt = 0, window_total = 0;
	double cov;

	for( i = 0; i < chromosome_length; i++)
	{
		rd_cnt += (long) in_bam->read_depth_per_chr[i];
		window_total++;
	}
	/* Calculate mu values */
	in_bam->mean = ( double)(rd_cnt) / ( double)window_total;

	fprintf( logFile, "Read Count:%li  Window count:%li mean=%f\n", rd_cnt, window_total, in_bam->mean);

	//cov = ( double) rd_cnt * (in_bam->libraries[lib_index]->read_length) / ( double) ( window_total * WINDOWSIZE );
	//fprintf( logFile,"Coverage is %.0lfx\n", round(cov));

}


void calc_mean_per_chr( parameters *params, bam_info* in_bam, int chr_index)
{
	int lib_index, i, gc_val = -1, window_per_gc[101], end;
	long rd_per_gc[101];

	calc_mu_per_chr( in_bam, params->this_sonic->chromosome_lengths[chr_index]);

	/* Calculate mu_GC values */
	for( i = 0; i < 101; i++)
	{
		rd_per_gc[i] = 0;
		window_per_gc[i] = 0;
	}

	for( i = 0; i < params->this_sonic->chromosome_lengths[chr_index]; i++)
	{
		if((i + WINDOWSLIDE) < params->this_sonic->chromosome_lengths[chr_index])
			end = i + WINDOWSLIDE;
		else
			end = params->this_sonic->chromosome_lengths[chr_index];

		gc_val = (int) round ( sonic_get_gc_content( params->this_sonic, params->this_sonic->chromosome_names[chr_index], i, end));
		rd_per_gc[gc_val] += ( long) in_bam->read_depth_per_chr[i];
		window_per_gc[gc_val]++;
	}

	in_bam->mean_rd_per_gc[0] = 0.0;
	for( i = 1; i < 101; i++)
	{
		in_bam->mean_rd_per_gc[i] = ( float)rd_per_gc[i] / ( window_per_gc[i]);
		if( isnanf( in_bam->mean_rd_per_gc[i]) || isinff( ( in_bam->mean_rd_per_gc[i])) == -1
				|| isinff( ( in_bam->mean_rd_per_gc[i])) == 1 )
			in_bam->mean_rd_per_gc[i] = 0;
		//fprintf(stderr,"GC = %d - RD=%ld\tWINDOW=%d\tMEAN=%f\n", i, rd_per_gc[i],window_per_gc[i], in_bam->mean_rd_per_gc[i]);
	}
}
