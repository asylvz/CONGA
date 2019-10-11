#include "common.h"
#include <math.h>
#include "free.h"
#include "likelihood.h"


void get_sample_name(bam_info* in_bam, char* header_text)
{
	char *tmp_header = NULL;
	set_str( &( tmp_header), header_text);
	char* p = strtok( tmp_header, "\t\n");
	char sample_name_buffer[1024];

	while( p != NULL)
	{
		/* If the current token has "SM" as the first two characters,
			we have found our Sample Name */
		if( p[0] == 'S' && p[1] == 'M')
		{
			/* Get the Sample Name */
			strncpy( sample_name_buffer, p + 3, strlen( p) - 3);

			/* Add the NULL terminator */
			sample_name_buffer[strlen( p) - 3] = '\0';

			/* Exit loop */
			break;
		}
		p = strtok( NULL, "\t\n");
	}

	set_str( &( in_bam->sample_name), sample_name_buffer);
	free( tmp_header);
}

void count_reads_bam( bam_info* in_bam, parameters* params)
{
	bam1_core_t bam_alignment_core;
	bam1_t* bam_alignment = bam_init1();

	while( bam_itr_next( in_bam->bam_file, in_bam->iter, bam_alignment) > 0)
	{
		bam_alignment_core = bam_alignment->core;

		/* Increase the read depth and read count for RD filtering */
		in_bam->read_depth_per_chr[bam_alignment_core.pos]++;
		in_bam->read_count++;
	}
	bam_destroy1( bam_alignment);
}


void read_bam( bam_info* in_bam, parameters *params)
{
	int i, bam_index, chr_index, chr_index_bam, return_value, not_in_bam = 0;
	char svfile_del[MAX_SEQ], svfile_dup[MAX_SEQ], svfile[MAX_SEQ];
	FILE *fpDel = NULL, *fpDup = NULL, *fpSVs = NULL;
	char chr[50];

	sprintf( svfile, "%s%s_svs.bed", params->outdir, params->outprefix);
	fprintf( stderr, "\nOutput SV file: %s\n", svfile);
	fpSVs = safe_fopen( svfile,"w");
	fprintf(fpSVs,"#CHR\tSTART_SV\tEND_SV\tSV_TYPE\tLIKELIHOOD\tCOPY_NUMBER\n");

	sprintf( svfile_del, "%s%s_dels.bed", params->outdir, params->outprefix);
	fprintf( stderr, "Output Del file: %s\n", svfile_del);
	fpDel = safe_fopen( svfile_del,"w");

	sprintf( svfile_dup, "%s%s_dups.bed", params->outdir, params->outprefix);
	fprintf( stderr, "Output DUP file: %s\n", svfile_dup);
	fpDup = safe_fopen( svfile_dup,"w");


	for( chr_index = 0; chr_index < params->this_sonic->number_of_chromosomes; chr_index++)
	{
		if (chr_index < params->first_chrom)
			chr_index = params->first_chrom;

		if (chr_index > params->last_chrom)
		{
			chr_index = params->this_sonic->number_of_chromosomes;
			continue;
		}

		/* HTS implementation */
		in_bam->bam_file = safe_hts_open( params->bam_file, "r");

		/* Read in BAM header information */
		in_bam->bam_header = bam_hdr_read( ( in_bam->bam_file->fp).bgzf);

		/* Load the bam index file */
		in_bam->bam_file_index = sam_index_load( in_bam->bam_file, params->bam_file);
		if( in_bam->bam_file_index == NULL)
		{
			fprintf( stderr, "Error: Sam Index cannot be loaded (sam_index_load)\n");
			exit( 1);
		}

		chr_index_bam = find_chr_index_bam( params->this_sonic->chromosome_names[chr_index], in_bam->bam_header);
		not_in_bam = 0;
		if( chr_index_bam == -1)
		{
			fprintf( stderr, "\nCannot find chromosome name %s in BAM/CRAM %s", params->this_sonic->chromosome_names[chr_index], in_bam->sample_name);
			not_in_bam = 1;
			continue;
		}

		in_bam->iter = sam_itr_queryi( in_bam->bam_file_index, chr_index_bam, 0, params->this_sonic->chromosome_lengths[chr_index]);
		if( in_bam->iter == NULL)
		{
			fprintf( stderr, "Error: Iterator cannot be loaded (bam_itr_queryi)\n");
			exit( 1);
		}

		fprintf( stderr, "\n");
		fprintf( stderr, "Reading BAM [%s] - Chromosome: %s", in_bam->sample_name, in_bam->bam_header->target_name[chr_index_bam]);

		strcpy(chr, in_bam->bam_header->target_name[chr_index_bam]);

		/* Initialize the read depth and read count */
		init_rd_per_chr( in_bam, params, chr_index);

		/* Read bam file for this chromosome */
		count_reads_bam( in_bam, params);

		/* Mean value (mu) calculation */
		calc_mean_per_chr( params, in_bam, chr_index);

		/* Close the BAM file */
		return_value = hts_close( in_bam->bam_file);
		if( return_value != 0)
		{
			fprintf( stderr, "Error closing BAM file\n");
			exit( 1);
		}
		/* Free the bam related files */
		sam_itr_destroy( in_bam->iter);
		bam_hdr_destroy( in_bam->bam_header);
		hts_idx_destroy(in_bam->bam_file_index);

		if( not_in_bam == 1)
			continue;

		find_SVs( in_bam, params, fpDel, fpDup, fpSVs, chr);
	}
	fprintf( stderr, "\n");
	fclose( fpDel);
	fclose( fpDup);
	fclose( fpSVs);
}
