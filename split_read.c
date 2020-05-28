/*
 * split_read.c
 *
 *  Created on: Nov 18, 2019
 *      Author: tardis
 */

#include "split_read.h"
#include <htslib/faidx.h>
#include <time.h>

#define HASH_COUNT 0
#define HASH_BUILD 1

#define SR_LOOKAHEAD 100000


int **hash_table_array;
int *hash_table_count;
int *hash_table_iter;
float total_hash_build_time = 0;

long split_read_count = 0;

void free_hash_table(parameters *params)
{
	int i;
	int hash_size = params->hash_size;

	for(i = 0; i < hash_size; i++){
		if(hash_table_count[i] != 0)
			free (hash_table_array[i]);
	}
	free(hash_table_array);

	free(params->ref_seq);
	params->ref_seq = NULL;
}

unsigned int hash_function_next( unsigned int prev_hash, unsigned int mask, const char next_char)
{
	/* this strictly assumes HASHKMERLEN < 16 */
	return (((prev_hash & mask) << 2) |  ((next_char & 0x6) >> 1));
}

unsigned int hash_function_ref( char *str)
{
	/* this strictly assumes HASHKMERLEN < 16 and is_kmer_valid is already called and returned TRUE */

	int i = 0;
	unsigned int val = 0; unsigned int numericVal = 0;
	while(i < HASHKMERLEN)
	{
		numericVal = (str[i++] & 0x6) >> 1;
		val = (val << 2) | numericVal;
	}
	return val;
}


void init_hash_count(parameters *params)
{
	params->hash_size = pow (4, HASHKMERLEN);
	hash_table_count = (int *) getMem( params->hash_size * sizeof (int));
	hash_table_iter = (int *) getMem( params->hash_size * sizeof (int));
	params->ref_seq = NULL;
}

int is_kmer_valid (char *str){

	int i, l;
	l = strlen(str);

	if (l < HASHKMERLEN)
		return 0;

	for (i = 0; i < HASHKMERLEN; i++)
		if (str[i] != 'A' && str[i] != 'C' && str[i] != 'G' && str[i] != 'T')
			return 0;

	return 1;
}

posMapSplitRead *almostPerfect_match_seq_ref( parameters *params, int chr_index, char *str, int pos)
{
	int i, index, posMapSize, posMap[10000], hammingDisMap[10000];
	char orient[10000];// orient of the mapping
	int dist, reverseMatch;
	lociInRef *ptr;
	char seed[HASHKMERLEN + 1];
	char *strRev;
	posMapSplitRead *tmpSplitMap = NULL, *returnPtr = NULL;

	int *hash_ptr;
	int num_hits;
	int cnt_hits;
	int str_length;
	int dist_max;

	returnPtr = NULL;

	str_length = strlen(str);
	dist_max = 0.05 * str_length;

	if( str_length < HASHKMERLEN)
		return NULL;

	strncpy( seed, str, HASHKMERLEN);
	seed[HASHKMERLEN] = '\0';

	if ( is_kmer_valid (seed) )
	{
		index = hash_function_ref( seed);
		num_hits = hash_table_count[index];
		hash_ptr = hash_table_array[index];
	}
	else
	{
		num_hits = 0;
		hash_ptr = NULL;
	}

	posMapSize = 0;

	for (cnt_hits = 0; cnt_hits < num_hits; cnt_hits++)
	{
		if ( abs (hash_ptr[cnt_hits] - pos) < SR_LOOKAHEAD)
		{
			dist = hammingDistance( &( params->ref_seq[hash_ptr[cnt_hits]]), str, str_length);
			if( dist <= dist_max)
			{
				posMap[posMapSize] = hash_ptr[cnt_hits];
				orient[posMapSize] = FORWARD;
				hammingDisMap[posMapSize] = dist;
				posMapSize = posMapSize + 1;
			}
		}
	}

	if( posMapSize < MAX_MAPPING)
	{
		strRev = ( char *)getMem( (str_length + 1) * sizeof( char));
		for( i = 0; i < str_length; i++)
		{
			if( str[i] == 'A')
				strRev[str_length - i - 1] = 'T';
			else if( str[i] == 'T')
				strRev[str_length - i - 1] = 'A';
			else if( str[i] == 'G')
				strRev[str_length - i - 1] = 'C';
			else if( str[i] == 'C')
				strRev[str_length - i - 1] = 'G';
			else if( str[i] == 'N')
				strRev[str_length - i - 1] = 'N';
		}

		strRev[str_length] = '\0';
		strncpy( seed, strRev, HASHKMERLEN);
		seed[HASHKMERLEN] = '\0';

		if (is_kmer_valid( seed))
		{
			index = hash_function_ref( seed);
			num_hits = hash_table_count[index];
			hash_ptr = hash_table_array[index];
		}
		else {
			hash_ptr = NULL;
			num_hits = 0;
		}

		reverseMatch = 0;
		for (cnt_hits = 0; cnt_hits < num_hits; cnt_hits++)
		{
			if ( abs (hash_ptr[cnt_hits] - pos) < SR_LOOKAHEAD)
			{
				dist = hammingDistance( &( params->ref_seq[hash_ptr[cnt_hits]]), strRev, str_length);
				if( dist <= dist_max)
				{
					posMap[posMapSize] = hash_ptr[cnt_hits];
					orient[posMapSize] = REVERSE;
					hammingDisMap[posMapSize] = dist;
					posMapSize = posMapSize + 1;
					reverseMatch = 1;
				}
			}
			if( posMapSize > MAX_MAPPING)
				break;
		}
		freeMem( strRev, str_length+1);
	}

	//posMapSize < 12 &&
	if( posMapSize > 0 && posMapSize < MAX_MAPPING)
	{
		for( i = 0; i < posMapSize; i++)
		{
			tmpSplitMap = ( posMapSplitRead *) getMem( sizeof( posMapSplitRead));
			tmpSplitMap->posMap = posMap[i];
			tmpSplitMap->orient = orient[i];

			//if( posMapSize < 11)
			tmpSplitMap->mapq = 60 / posMapSize;
			//else
			//posMapSize = 0;

			tmpSplitMap->next = returnPtr;
			returnPtr = tmpSplitMap;
		}
	}

	return returnPtr;
}

int find_split_reads( bam_info* in_bam, parameters* params, bam1_t* bam_alignment, int chr_index)
{
	uint8_t *tmp;
	int return_type, i;
	float avgPhredQual = 0;
	char str[512], str2[512];
	uint8_t *a_qual, * a_qual2;

	bam1_core_t bam_alignment_core = bam_alignment->core;

	if( bam_alignment_core.pos == 0)
		return -1;
	else
	{
		/* Split the read */
		splitRead *newEl = ( splitRead *) getMem( sizeof( splitRead));

		newEl->readName = NULL;
		set_str( &(newEl->readName), bam_alignment->data);

		/* Get the name of the chromosome */
		newEl->chromosome_name = NULL;
		set_str( &(newEl->chromosome_name), params->this_sonic->chromosome_names[chr_index]);

		newEl->pos = bam_alignment_core.pos;
		newEl->qual = bam_alignment_core.qual;
		newEl->orient = FORWARD;
		newEl->split_start = (bam_alignment_core.l_qseq / 2);
		newEl->read_length = bam_alignment_core.l_qseq;

		a_qual = bam_get_qual( bam_alignment);

		for( i = newEl->split_start; i < bam_alignment_core.l_qseq; i++)
			avgPhredQual = avgPhredQual + a_qual[i];

		avgPhredQual = ( float)avgPhredQual / ( float)(bam_alignment_core.l_qseq - newEl->split_start);

		newEl->avgPhredQualSplitRead = ( int)floorf( avgPhredQual);

		if(newEl->avgPhredQualSplitRead < params->mq_threshold)
		{
			free(newEl->chromosome_name);
			free(newEl->readName);
			free(newEl);
			return -1;
		}


		//str = (char *) getMem( ( bam_alignment_core.l_qseq / 2 + 1) * sizeof( char));
		newEl->ptrSplitMap = NULL;

		int k = 0;
		for( i = newEl->read_length / 2; i < newEl->read_length; i++)
		{
			if( bam_seqi( bam_get_seq( bam_alignment), i) == 1)
				str[k] = 'A';
			else if( bam_seqi( bam_get_seq( bam_alignment), i) == 2)
				str[k] = 'C';
			else if(bam_seqi( bam_get_seq( bam_alignment), i) == 4)
				str[k] = 'G';
			else if( bam_seqi( bam_get_seq( bam_alignment), i) == 8)
				str[k] = 'T';
			else if( bam_seqi( bam_get_seq( bam_alignment), i) == 15)
				str[k] = 'N';
			k++;
		}
		str[k] = '\0';

		if(str != NULL)
		{
			newEl->ptrSplitMap = almostPerfect_match_seq_ref( params, chr_index, str, newEl->pos);
			//free(str);
		}

		newEl->next = in_bam->listSplitRead;
		in_bam->listSplitRead = newEl;

		split_read_count++;

		//Split the other part also

		char read_name_2[1000];

		splitRead *newEl2 = ( splitRead *) getMem( sizeof( splitRead));

		sprintf( read_name_2, "%s_read2", bam_alignment->data);
		newEl2->readName = NULL;
		set_str( &(newEl2->readName), read_name_2);

		/* Get the name of the chromosome */
		newEl2->chromosome_name = NULL;
		set_str( &(newEl2->chromosome_name), params->this_sonic->chromosome_names[chr_index]);

		newEl2->pos = bam_alignment_core.pos + (bam_alignment_core.l_qseq / 2);
		newEl2->qual = bam_alignment_core.qual;
		newEl2->orient = FORWARD;
		newEl2->split_start = 0;
		newEl2->read_length = bam_alignment_core.l_qseq;

		a_qual2 = bam_get_qual( bam_alignment);

		for( i = 0; i < bam_alignment_core.l_qseq / 2; i++)
			avgPhredQual = avgPhredQual + a_qual2[i];

		avgPhredQual = ( float)avgPhredQual / ( float)(bam_alignment_core.l_qseq / 2);

		newEl2->avgPhredQualSplitRead = ( int)floorf( avgPhredQual);

		if(newEl2->avgPhredQualSplitRead < params->mq_threshold)
		{
			free(newEl2->chromosome_name);
			free(newEl2->readName);
			free(newEl2);
			return -1;
		}

		newEl2->ptrSplitMap = NULL;

		k = 0;
		for( i = 0; i < newEl2->read_length / 2; i++)
		{
			if( bam_seqi( bam_get_seq( bam_alignment), i) == 1)
				str2[k] = 'A';
			else if( bam_seqi( bam_get_seq( bam_alignment), i) == 2)
				str2[k] = 'C';
			else if(bam_seqi( bam_get_seq( bam_alignment), i) == 4)
				str2[k] = 'G';
			else if( bam_seqi( bam_get_seq( bam_alignment), i) == 8)
				str2[k] = 'T';
			else if( bam_seqi( bam_get_seq( bam_alignment), i) == 15)
				str2[k] = 'N';
			k++;
		}
		str2[k] = '\0';

		if(str2 != NULL)
		{
			newEl2->ptrSplitMap = almostPerfect_match_seq_ref( params, chr_index, str2, newEl2->pos);
			//free(str);
		}

		newEl2->next = in_bam->listSplitRead;
		in_bam->listSplitRead = newEl2;

		split_read_count++;

		return RETURN_SUCCESS;
	}
}


void build_hash_table(const char *ref, int len, int hash_size, int mode)
{
	int i = 0, j = 0;
	char seed[HASHKMERLEN + 1];
	int hash_val;
	int mask;
	int len_limit;

	mask = (hash_size - 1) >> 2;

	if (mode == HASH_COUNT)
		memset (hash_table_count, 0, hash_size * sizeof (int));
	else
		memset (hash_table_iter, 0, hash_size * sizeof (int));


	len_limit = len - HASHKMERLEN + 1;

	/* get first kmer */
	while (!is_dna_letter(ref[i])) i++;

	strncpy (seed, ref+i, HASHKMERLEN);
	seed[HASHKMERLEN] = 0;

	while (!is_kmer_valid (seed)){
		i++;
		strncpy (seed, ref+i, HASHKMERLEN);
		seed[HASHKMERLEN] = 0;
	}

	hash_val = hash_function_ref (seed);
	if (mode == HASH_COUNT)
	{
		( hash_table_count[hash_val])++;
	}
	else if (hash_table_count[hash_val] != 0)
	{
		hash_table_array[hash_val][hash_table_iter[hash_val]] = i;
		( hash_table_iter[hash_val])++;
	}

	j = i + HASHKMERLEN;

	while (i < len_limit)
	{
		i++; //shift
		if (is_dna_letter(ref[j]))
		{
			hash_val = hash_function_next( hash_val, mask, ref[j++]);
			if (mode == HASH_COUNT)
				( hash_table_count[hash_val])++;
			else if (hash_table_count[hash_val] != 0)
			{
				hash_table_array[hash_val][hash_table_iter[hash_val]] = i;
				( hash_table_iter[hash_val])++;
			}
		}
		else
		{
			/* recover from non-ACGT */
			while (!is_dna_letter(ref[j]) && i < len_limit) { i++; j++; }
			if (i >= len_limit)
				break;

			strncpy (seed, ref+i, HASHKMERLEN);
			seed[HASHKMERLEN] = 0;

			while (!is_kmer_valid (seed)){
				i++;
				strncpy (seed, ref+i, HASHKMERLEN);
				seed[HASHKMERLEN] = 0;
			}

			j = i + HASHKMERLEN;
			hash_val = hash_function_ref (seed);

			if (mode == HASH_COUNT)
				( hash_table_count[hash_val])++;
			else if (hash_table_count[hash_val] != 0)
			{
				hash_table_array[hash_val][hash_table_iter[hash_val]] = i;
				( hash_table_iter[hash_val])++;
			}
		}
	}
}

void init_hash_table(parameters *params)
{
	int i;
	int hash_size = params->hash_size;
	hash_table_array = (int **) getMem (hash_size * sizeof (int *));

	for (i = 0; i < hash_size; i++)
	{
		if (hash_table_count[i] != 0 && hash_table_count[i] < MAX_SR_HIT)
			hash_table_array[i] = (int *) getMem (hash_table_count[i] * sizeof (int));
		else
		{
			hash_table_count[i] = 0;
			hash_table_array[i] = NULL;
		}
	}
}

void create_hash_table( parameters *params, int len)
{
	init_hash_table( params);
	build_hash_table( params->ref_seq, len, params->hash_size, HASH_BUILD);
}



/*void map_split_reads( bam_info* in_bam, parameters* params, int chr_index)
{
	int i, tmp;
	char *str;
	splitRead *ptrSplitRead;

	ptrSplitRead = in_bam->listSplitRead;
	while( ptrSplitRead != NULL)
	{
		if( ptrSplitRead->avgPhredQualSplitRead > params->mq_threshold && ptrSplitRead->split_sequence != NULL)
		{
			//&& ptrSplitRead->numSoftClipInConcordance > 1
			tmp = ptrSplitRead->read_length - ptrSplitRead->split_start;
			if (tmp > 0)
			{
				str = ( char *) getMem( sizeof( char) * (tmp + 1));
				strncpy( str, &( ptrSplitRead->split_sequence[ptrSplitRead->split_start]), tmp);
				str[tmp] = '\0';

				ptrSplitRead->ptrSplitMap = almostPerfect_match_seq_ref( params, chr_index, str, ptrSplitRead->pos);
				//fprintf(stderr,"DONEE - %d\n",z);
				if( str != NULL)
					free( str);
			}
		}
		ptrSplitRead = ptrSplitRead->next;
	}
}*/


/* Read the reference genome */
void readReferenceSeq( parameters *params, int chr_index)
{
	int i, min, max, loc_length;
	char *ref_seq;
	long bp_cnt = 0;
	faidx_t* ref_fai;

	if (params->ref_seq != NULL)
	{
		fprintf (stderr, "Reference genome is already loaded.\n");
		return;
	}

	min = 0, max = 999;
	ref_fai = fai_load( params->ref_genome);

	params->ref_seq = ( char *) getMem( (params->this_sonic->chromosome_lengths[chr_index] + 1) * sizeof( char));

	while ( max < params->this_sonic->chromosome_lengths[chr_index])
	{
		ref_seq = faidx_fetch_seq( ref_fai, params->this_sonic->chromosome_names[chr_index], min, max, &loc_length);

		for( i = 0; i < loc_length; i++)
		{
			/* can we do this faster with memcpy? */
			if( bp_cnt < params->this_sonic->chromosome_lengths[chr_index])
				params->ref_seq[bp_cnt] = toupper( ref_seq[i]);
			bp_cnt++;
		}
		if( bp_cnt >= params->this_sonic->chromosome_lengths[chr_index])
			break;

		min += loc_length;
		max += loc_length;
		free( ref_seq);
	}
	fai_destroy( ref_fai);

	params->ref_seq[bp_cnt] = '\0';

	build_hash_table(params->ref_seq, bp_cnt, params->hash_size, HASH_COUNT);

	create_hash_table(params, bp_cnt);
}
