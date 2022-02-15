#include "common.h"
#include <math.h>
#include "free.h"
#include "likelihood.h"
#include "split_read.h"


struct SplitsInfo *all_split_reads = NULL;

SplitRow *createSplitRow(int locMapLeftStart, int locMapLeftEnd, char orientationLeft,
		int locMapRightStart, int locMapRightEnd, char orientationRight, char svType)
{
	SplitRow *newRow = (SplitRow*) getMem(sizeof(SplitRow));

	newRow->next = NULL;

	newRow->locMapLeftEnd = locMapLeftEnd;
	newRow->locMapLeftStart = locMapLeftStart;
	newRow->locMapRightStart = locMapRightStart;
	newRow->locMapRightEnd = locMapRightEnd;
	newRow->orientationLeft = orientationLeft;
	newRow->orientationRight = orientationRight;

	newRow->svType = svType;

	return newRow;
}

SplitRow *determine_SvType( splitRead *ptrSplitRead, posMapSplitRead *ptrPosMapSplit)
{
	int pos1_1, pos1_2, pos2_1, pos2_2, i = 0;

	/* Length of read A(left) and read B(right) */
	int lengthRead, lengthSplit;
	bool is_split_beginning = false;

	if(strstr(ptrSplitRead->readName, "_read2") != NULL)
		is_split_beginning = true;

	/* If soft clip is at the end */
	lengthSplit = ptrSplitRead->read_length / 2;
	lengthRead = ptrSplitRead->read_length - lengthSplit;

	if( ptrSplitRead->pos < ptrPosMapSplit->posMap)
	{
		pos1_1 = ptrSplitRead->pos;
		pos1_2 = ptrSplitRead->pos + lengthRead;
		pos2_1 = ptrPosMapSplit->posMap;
		pos2_2 = ptrPosMapSplit->posMap + lengthSplit;
	}
	else if( ptrPosMapSplit->posMap < ptrSplitRead->pos)
	{
		pos1_1 = ptrPosMapSplit->posMap;
		pos1_2 = ptrPosMapSplit->posMap + lengthSplit;
		pos2_1 = ptrSplitRead->pos;
		pos2_2 = ptrSplitRead->pos + lengthRead;
	}

	if( pos1_2 >= pos2_1)
		return NULL;

	if( ptrSplitRead->orient == FORWARD && ptrPosMapSplit->orient == FORWARD)
	{
		if( (( ptrSplitRead->pos < ptrPosMapSplit->posMap) && is_split_beginning == false) ||
				(( ptrSplitRead->pos > ptrPosMapSplit->posMap) && is_split_beginning == true))
		{
			SplitRow * newRow = createSplitRow (pos1_1, pos1_2, FORWARD, pos2_1, pos2_2, FORWARD, DELETION);
			return newRow;
		}
		else if( (( ptrPosMapSplit->posMap < ptrSplitRead->pos) && is_split_beginning == false) ||
				(( ptrPosMapSplit->posMap > ptrSplitRead->pos) && is_split_beginning == true))
		{
			SplitRow * newRow = createSplitRow (pos1_1, pos1_2, FORWARD, pos2_1, pos2_2, FORWARD, DUPLICATION);
			return newRow;
		}
	}
	return NULL;
}

int read_SplitReads(splitRead *ptrSoftClip, parameters *params, int chr_index)
{
	float is_satellite = 0.0;
	SplitRow *newRow;
	posMapSplitRead *ptrPosMapSoftClip = NULL;

	all_split_reads = ( SplitsInfo *) getMem( sizeof( struct SplitsInfo));
	all_split_reads->size = 0;
	all_split_reads->head = NULL;
	all_split_reads->tail = NULL;

	while( ptrSoftClip != NULL)
	{
		ptrPosMapSoftClip = ptrSoftClip->ptrSplitMap;
		while( ptrPosMapSoftClip != NULL)
		{
			is_satellite = sonic_is_satellite( params->this_sonic, ptrSoftClip->chromosome_name, ptrSoftClip->pos, ptrSoftClip->pos + 1 ) +
					sonic_is_satellite( params->this_sonic, ptrSoftClip->chromosome_name, ptrPosMapSoftClip->posMap, ptrPosMapSoftClip->posMap + 1);

			newRow = NULL;
			if ( is_satellite == 0 && ptrSoftClip->qual > params->mq_threshold && strcmp(ptrSoftClip->chromosome_name, params->this_sonic->chromosome_names[chr_index]) == 0
					&& ptrPosMapSoftClip->mapq > params->mq_threshold && ptrSoftClip->pos > 0 && ptrPosMapSoftClip->posMap > 0
					&& ptrSoftClip->pos < params->this_sonic->chromosome_lengths[chr_index] && ptrPosMapSoftClip->posMap < params->this_sonic->chromosome_lengths[chr_index])
			{
				newRow = determine_SvType(ptrSoftClip, ptrPosMapSoftClip);
				if( newRow != NULL)
				{
					/* For Deletion */
					if( newRow->svType == DELETION)
					{
						newRow->orientationLeft = FORWARD;
						newRow->orientationRight = REVERSE;
						newRow->locMapLeftStart -= SOFTCLIP_WRONGMAP_WINDOW;
						newRow->locMapLeftEnd -= SOFTCLIP_WRONGMAP_WINDOW;
						newRow->locMapRightStart += SOFTCLIP_WRONGMAP_WINDOW;
						newRow->locMapRightEnd += SOFTCLIP_WRONGMAP_WINDOW;
					}
					/* For Duplication */
					else if(newRow->svType == DUPLICATION)
					{
						newRow->orientationLeft = REVERSE;
						newRow->orientationRight = FORWARD;
						newRow->locMapLeftStart -= SOFTCLIP_WRONGMAP_WINDOW;
						newRow->locMapLeftEnd -= SOFTCLIP_WRONGMAP_WINDOW;
						newRow->locMapRightStart += SOFTCLIP_WRONGMAP_WINDOW;
						newRow->locMapRightEnd += SOFTCLIP_WRONGMAP_WINDOW;
					}
					else
						newRow = NULL;
				}

				if( newRow == NULL)
					;//fprintf( stderr, "ERROR loading divet from bam (soft clip)\n");
				else
				{
					if( all_split_reads->head == NULL || all_split_reads->tail == NULL)
					{
						all_split_reads->head = newRow;
						all_split_reads->tail = newRow;
					}
					else
					{
						all_split_reads->tail->next = newRow;
						all_split_reads->tail = newRow;
					}
					all_split_reads->size++;
				}
			}
			ptrPosMapSoftClip = ptrPosMapSoftClip->next;
		}
		ptrSoftClip = ptrSoftClip->next;
	}
	fprintf(stderr,"\nCONGA paired %d single-end reads\n", all_split_reads->size);
	return all_split_reads->size;
}

int write_sequences(parameters *params, bam1_t* bam_alignment, FILE* fp, int base_count)
{
	int i, k = 0;
	int hash[4];

	bam1_core_t bam_alignment_core = bam_alignment->core;
	char str[bam_alignment_core.l_qseq + 1];

	if( bam_alignment_core.pos == 0)
		return -1;


	for(i = 0; i < bam_alignment_core.l_qseq; i++)
	{
		if( bam_seqi( bam_get_seq( bam_alignment), i) == 1)
			str[k] = 'A';
		else if( bam_seqi( bam_get_seq( bam_alignment), i) == 2)
			str[k] = 'C';
		else if(bam_seqi( bam_get_seq( bam_alignment), i) == 4)
			str[k] = 'G';
		else if( bam_seqi( bam_get_seq( bam_alignment), i) == 8)
			str[k] = 'T';
		//else //In case of N...
		//return -1;
		k++;
	}
	base_count += k;

	str[k] = '\0';

	if(str != NULL)
		fprintf(fp,"%s",str);

	return base_count;
}

void count_reads_bam( bam_info* in_bam, parameters* params, int chr_index, int* base_count_bam )
{
	bam1_core_t bam_alignment_core;
	char seq_file[MAX_SEQ];
	int return_type;
	int cnt_reads = 0;

	bam1_t* bam_alignment = bam_init1();

	while( sam_itr_next( in_bam->bam_file, in_bam->iter, bam_alignment) > 0)
	{
		bam_alignment_core = bam_alignment->core;

		if(sonic_is_satellite( params->this_sonic, params->this_sonic->chromosome_names[chr_index], bam_alignment_core.pos, bam_alignment_core.pos + 20) == 0
				&& bam_alignment_core.qual > params->mq_threshold && is_proper( bam_alignment_core.flag))
		{
			if( !params->no_sr && params->dup_file && bam_alignment_core.l_qseq > params->min_read_length)
			{
				return_type = find_split_reads( in_bam, params, bam_alignment, chr_index);
			}
		}

		in_bam->read_depth[bam_alignment_core.pos]++;
		in_bam->total_read_count_unfiltered++;
		cnt_reads++;

	}
	fprintf(stderr," (%d reads, %ld split-reads)\n", cnt_reads, split_read_count);
	bam_destroy1( bam_alignment);
}


void read_bam( bam_info* in_bam, parameters *params)
{
	int i, bam_index, chr_index, chr_index_bam, return_value, not_in_bam = 0;
	char svfile_del[MAX_SEQ], svfile_dup[MAX_SEQ], svfile[MAX_SEQ];
	FILE *fpDel = NULL, *fpDup = NULL, *fpSVs = NULL;
	int base_count_bam = 0;
	long bp_cnt;

	sprintf( svfile, "%s%s_svs.bed", params->outdir, params->outprefix);
	fprintf( stderr, "\nOutput SV file: %s\n", svfile);
	fpSVs = safe_fopen( svfile,"w");
	fprintf(fpSVs,"#CHR\tSTART_SV\tEND_SV\tSV_TYPE\tCOPY_NUMBER\tLIKELIHOOD\tREAD_PAIR\tMAPPABILITY\n");

	if(params->del_file)
	{
		sprintf( svfile_del, "%s%s_dels.bed", params->outdir, params->outprefix);
		fprintf( stderr, "Output Del file: %s\n", svfile_del);
		fpDel = safe_fopen( svfile_del,"w");
		fprintf(fpDel,"#CHR\tSTART_SV\tEND_SV\tCOPY_NUMBER\tLIKELIHOOD\tREAD_PAIR\tMAPPABILITY\tOBSERVED_READS\tEXPECTED_READS\tNOSV_LIKELIHOOD\tHETEROZYGOUS_LIKELIHOOD\tHOMOZYGOUS_LIKELIHOOD\n");
	}
	if(params->dup_file)
	{
		sprintf( svfile_dup, "%s%s_dups.bed", params->outdir, params->outprefix);
		fprintf( stderr, "Output DUP file: %s\n", svfile_dup);
		fpDup = safe_fopen( svfile_dup,"w");
		fprintf(fpDup,"#CHR\tSTART_SV\tEND_SV\tCOPY_NUMBER\tLIKELIHOOD\tREAD_PAIR\tMAPPABILITYOBSERVED_READS\tEXPECTED_READS\tNOSV_LIKELIHOOD\tHETEROZYGOUS_LIKELIHOOD\tHOMOZYGOUS_LIKELIHOOD\n");
	}

	/* HTS implementation */
	in_bam->bam_file = safe_hts_open( params->bam_file, "r");

	/* Read in BAM header information */
	in_bam->bam_header = sam_hdr_read( in_bam->bam_file);

	/* Load the bam index file */
	in_bam->bam_file_index = sam_index_load( in_bam->bam_file, params->bam_file);
	if( in_bam->bam_file_index == NULL)
	{
		fprintf( stderr, "Error: Sam Index cannot be loaded (sam_index_load)\n");
		exit( 1);
	}

	/* Extract the Sample Name from the header text */
	get_sample_name( in_bam, in_bam->bam_header->text);

	for( chr_index = 0; chr_index < params->this_sonic->number_of_chromosomes; chr_index++)
	{
		if (chr_index < params->first_chrom)
			chr_index = params->first_chrom;

		if (chr_index > params->last_chrom)
		{
			chr_index = params->this_sonic->number_of_chromosomes;
			continue;
		}

		if( strstr( params->this_sonic->chromosome_names[chr_index], "X") != NULL || strstr( params->this_sonic->chromosome_names[chr_index], "Y") != NULL)
			continue;


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

		/* Initialize the read depth and read count */
		init_rd_per_chr( in_bam, params, chr_index);

		if(!params->no_sr && params->dup_file)
		{
			fprintf( stderr, "\nReading the Reference Genome");
			bp_cnt = readReferenceSeq(params, chr_index);
			build_hash_table(params->ref_seq, bp_cnt, params->hash_size, HASH_COUNT);
			create_hash_table(params, bp_cnt);
		}

		/* Read bam file for this chromosome */
		fprintf(stderr,"\n-->counting reads");
		count_reads_bam( in_bam, params, chr_index, &base_count_bam);

		if(!params->no_sr && params->dup_file)
			free_hash_table(params);

		/* Mean value (mu) calculation */
		calc_mean_per_chr( params, in_bam, chr_index);

		/* Free the bam related files */
		sam_itr_destroy( in_bam->iter);

		if( not_in_bam == 1)
			continue;

		//Load Split-Reads
		if(!params->no_sr && params->dup_file)
		{
			fprintf(stderr, "Loading Split-Reads\n");
			read_SplitReads(in_bam->listSplitRead, params, chr_index);
		}

		//fprintf( stderr, "\nLikelihood Estimation\n");
		find_SVs( in_bam, params, fpDel, fpDup, fpSVs, in_bam->bam_header->target_name[chr_index_bam], chr_index);
	}

	/* Close the BAM file */
	return_value = hts_close( in_bam->bam_file);
	if( return_value != 0)
	{
		fprintf( stderr, "Error closing BAM file\n");
		exit( 1);
	}

	bam_hdr_destroy(in_bam->bam_header);
	hts_idx_destroy(in_bam->bam_file_index);

	fprintf( stderr, "\n");

	if(fpDel)
		fclose( fpDel);
	if(fpDup)
		fclose( fpDup);
	fclose( fpSVs);
}
