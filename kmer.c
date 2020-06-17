/*
 * kmer.c
 *
 *  Created on: Apr 8, 2020
 *      Author: tardis
 */
#include <stdio.h>
#include "kmer.h"
#include "common.h"
#include <htslib/faidx.h>
#include <math.h>
#include "mhash.h"

svs* all_svs_del;
svs* all_svs_dup;

HashInfo **hash_table_kmer;

int isPrime(int n)
{
	int i;
	for(i = 2; i <= sqrt(n); i++)
		if( n % i == 0 ) return 0;
	return 1;
}
int firstPrime(int n)
{
	while( !isPrime(n) )
		n++;
	return n;
}

void free_hash_table_kmer(parameters *params)
{
	int i;

	HashInfo *tmp = NULL, *tmp_next = NULL;
	for(i = 0; i < params->hash_size_kmer; i++)
	{
		tmp = hash_table_kmer[i];
		while(tmp != NULL)
		{
			tmp_next = tmp->next;
			free(tmp);
			tmp = tmp_next;
		}
	}
	free(hash_table_kmer);
	hash_table_kmer = NULL;
}

int is_kmer_valid_likelihood (char *str)
{
	int i, l = 0;
	l = (int) strlen(str);

	if (l != KMER)
		return 0;

	for (i = 0; i < KMER; i++)
	{
		if (str[i] != 'A' && str[i] != 'C' && str[i] != 'G' && str[i] != 'T')
			return 0;
	}

	return 1;
}

int legitimateSeed(char* seed)
{
	int i;
	for (i = 0; i < KMER; i++)
	{
		if(!is_dna_letter(seed[i]))
			return 0;
	}
	return 1;
}

/* Read the reference genome */
int read_kmer_jellyfish( parameters *params)
{
	int i, j, mer_count = 0;
	char *kmer_seq, *return_value, *tmp;
	unsigned int hash_val1, hash_val2;
	FILE* fp_reads = NULL;
	char seqFile[MAX_SEQ];
	int hash[4];
	char line[64];
	int hash_size, kmer_freq;
	HashInfo *newEl = NULL;

	//tmp_file = fopen("dnm_seqs", "w");
	sprintf( seqFile, "mer_counts.fa");
	fp_reads = fopen(seqFile, "r");
	if (fp_reads == NULL){
		printf("Could not open file %s", seqFile);
		exit(1);
	}
	while(!feof(fp_reads))
	{
		return_value = fgets(line, 64, fp_reads);
		if (line == NULL)
			continue;

		/* If the read line is empty */
		int isEmpty = 1;
		int len = strlen (line);
		for (i = 0; i < len; i++)
			if (!isspace (line[i]))
			{
				isEmpty = 0;
				break;
			}
		if (isEmpty)
			continue;

		tmp = strtok (line, ROW_DELIMITERS);

		if(strstr(tmp,">") != NULL)
			mer_count++;
		line[0] = '\0';
	}
	rewind(fp_reads);
	//fprintf(stderr,"There are %d lines", mer_count);

	hash_size = firstPrime(mer_count * 1.2);
	params->hash_size_kmer = hash_size;

	hash_table_kmer = (HashInfo **) getMem ((hash_size + 1) * sizeof (HashInfo*));

	for(i = 0; i < hash_size; i++)
		hash_table_kmer[i] = NULL;

	i = 0;

	// Loops through the frequency file
	while(!feof(fp_reads))
	{
		// This is executed once we get the frequency after this if
		// First we get the frequency, then the sequence
		if(i != 0 && (i % 2) == 0)
		{
			MurmurHash3_x86_128(kmer_seq, strlen(kmer_seq), 42, hash);
			hash_val1 = (unsigned) hash[0] % hash_size;

			MurmurHash3_x86_128(kmer_seq, strlen(kmer_seq), 11, hash);
			hash_val2 = (unsigned) hash[0] % hash_size;

			//if(hash_val1 == hash_val2)
			//fprintf(stderr,"PROBLEMM IN HASH - Same values (kmer.c)\n");

			newEl = (HashInfo *) getMem(sizeof(HashInfo));
			newEl->freq = kmer_freq;
			newEl->hash2 = hash_val2;

			//fprintf(tmp_file,"%s\t%ld\t%u\t%u\t%d\n", kmer_seq, strlen(kmer_seq), hash_val1, hash_val2, newEl->freq);
			newEl->next = hash_table_kmer[hash_val1];
			hash_table_kmer[hash_val1] = newEl;

			free(kmer_seq);
			kmer_seq = NULL;
		}

		return_value = fgets(line, 64, fp_reads);

		if (line == NULL)
			continue;

		/* If the read line is empty */
		int isEmpty = 1;
		int len = strlen (line);
		for (j = 0; j < len; j++)
			if (!isspace (line[j]))
			{
				isEmpty = 0;
				break;
			}
		if (isEmpty)
			continue;

		if((i % 2) == 0)
		{
			tmp = strtok (line, ROW_DELIMITERS_MERS);
			kmer_freq = (short) atoi(tmp);
			//fprintf(stderr,"%d\t%s\n",i, kmer_freq);
		}
		else
		{
			tmp = strtok (line, ROW_DELIMITERS);
			kmer_seq = NULL;
			set_str(&kmer_seq, tmp);
		}
		line[0] = '\0';
		i++;
	}
	fclose(fp_reads);

	return hash_size;
}


char* read_ref_seq( parameters *params, char* chr_name, int start, int end)
{
	int i, min, max, loc_length, chr_index;
	char *ref_seq = NULL;
	faidx_t* ref_fai;

	chr_index = sonic_refind_chromosome_index( params->this_sonic, chr_name);

	ref_fai = fai_load( params->ref_genome);
	ref_seq = faidx_fetch_seq( ref_fai, params->this_sonic->chromosome_names[chr_index], start, end, &loc_length);

	fai_destroy(ref_fai);

	return ref_seq;
}


int kmer_count_interval(long hash_size_kmer, int variation_length, char* seq)
{
	int i, j;
	char *seq_tmp = NULL, *seq_tmp_rev = NULL;
	int hash[4];
	unsigned int hash_val1, hash_val2;
	HashInfo *tmp = NULL;
	int kmer_count_forward = 0, kmer_count_reverse = 0;
	int is_kmer_valid = 0;


	for(i = 0; i < variation_length - KMER; i += KMERSLIDE)
	{
		// For forward strand
		seq_tmp = substring(seq, i, KMER);

		is_kmer_valid = is_kmer_valid_likelihood(seq_tmp);
		if(!is_kmer_valid)
		{
			if(seq_tmp != NULL)
				free(seq_tmp);
			continue;
		}
		//fprintf(stderr,"HEREE\n");
		//fprintf(stderr,"%s\n",seq_tmp);

		MurmurHash3_x86_128(seq_tmp, strlen(seq_tmp), 42, hash);
		hash_val1 = (unsigned) hash[0] % hash_size_kmer;

		MurmurHash3_x86_128(seq_tmp, strlen(seq_tmp), 11, hash);
		hash_val2 = (unsigned) hash[0] % hash_size_kmer;
		//fprintf(stderr,"%u - %u\n", hash_val1, hash_val2);

		tmp = hash_table_kmer[hash_val1];
		while(tmp != NULL)
		{
			if(tmp->hash2 == hash_val2)
			{
				kmer_count_forward += tmp->freq;
				break;
			}
			tmp = tmp->next;
		}

		//For reverse strand
		seq_tmp_rev = reverseComplement(seq_tmp);

		is_kmer_valid = is_kmer_valid_likelihood(seq_tmp_rev);
		if(!is_kmer_valid)
		{
			if(seq_tmp_rev != NULL)
			{
				free(seq_tmp_rev);
				seq_tmp_rev = NULL;
			}
			if(seq_tmp != NULL)
			{
				free(seq_tmp);
				seq_tmp = NULL;
			}
			continue;
		}
		//fprintf(stderr,"REV = %s\n",seq_tmp_rev);
		MurmurHash3_x86_128(seq_tmp_rev, strlen(seq_tmp_rev), 42, hash);
		hash_val1 = (unsigned) hash[0] % hash_size_kmer;

		MurmurHash3_x86_128(seq_tmp_rev, strlen(seq_tmp_rev), 11, hash);
		hash_val2 = (unsigned) hash[0] % hash_size_kmer;

		free(seq_tmp);
		seq_tmp = NULL;

		free(seq_tmp_rev);
		seq_tmp_rev = NULL;

		tmp = hash_table_kmer[hash_val1];
		while(tmp != NULL)
		{
			if(tmp->hash2 == hash_val2)
			{
				kmer_count_reverse += tmp->freq;
				break;
			}
			tmp = tmp->next;
		}
	}
	return (kmer_count_forward + kmer_count_reverse);
}

void init_kmer_per_chr( bam_info* in_bam, parameters* param, int chr_index)
{
	in_bam->kmer = ( short*) getMem( sizeof( short) * ( param->this_sonic->chromosome_lengths[chr_index]));
	memset (in_bam->kmer, 0, (param->this_sonic->chromosome_lengths[chr_index] * sizeof(short)));
}

long calc_kmer_counts(bam_info *in_bam, parameters *params, int chr_index)
{
	int i, j, k, end;
	short tmp = 0;
	long total_kmers = 0;
	char* seq = NULL;
	//char *sub_seq = NULL;
	char sub_seq[KMERWINDOWSIZE + 10];

	seq = read_ref(params, chr_index);

	for( i = 0; i < params->this_sonic->chromosome_lengths[chr_index]; i += KMERWINDOWSLIDE)
	{
		end = i;

		if((i + KMERWINDOWSIZE) < params->this_sonic->chromosome_lengths[chr_index])
			end = i + KMERWINDOWSIZE;
		else
			end = params->this_sonic->chromosome_lengths[chr_index];

		//sub_seq = ( char*) getMem( sizeof( char) * ((end - i) + 1));

		k = 0;
		for(j = i; j < end; j++)
		{
			sub_seq[k] = seq[j];
			k++;
		}
		sub_seq[k] = '\0';

		tmp = 0;
		tmp = (short) kmer_count_interval(params->hash_size_kmer, end - i + 1, sub_seq);
		total_kmers += (long) tmp;

		for(j = 0; j < (end - i); j++)
			in_bam->kmer[i + j] = tmp;

		//sub_seq[0] = '\0';
		memset(&sub_seq[0], 0, sizeof(sub_seq));
	}
	//fprintf(stderr,"here3");

	if(seq != NULL)
		free(seq);
	seq = NULL;
	//fprintf(stderr,"here4");

	return total_kmers;
}

void calc_expected_kmer(bam_info *in_bam, parameters *params, int chr_index)
{
	int i, gc_val = -1, window_per_gc[101], end;
	long rd_per_gc[101];

	/* Calculate mu_GC values */
	for( i = 0; i < 101; i++)
	{
		rd_per_gc[i] = 0;
		window_per_gc[i] = 0;
		in_bam->expected_kmer[i] = 0.0;
	}

	for( i = 0; i < params->this_sonic->chromosome_lengths[chr_index]; i += KMERWINDOWSLIDE)
	{
		if((i + KMERWINDOWSIZE) < params->this_sonic->chromosome_lengths[chr_index])
			end = i + KMERWINDOWSIZE;
		else
			end = params->this_sonic->chromosome_lengths[chr_index];

		gc_val = (int) round (sonic_get_gc_content(params->this_sonic, params->this_sonic->chromosome_names[chr_index], i, end));
		rd_per_gc[gc_val] += (long) in_bam->kmer[i];
		window_per_gc[gc_val]++;
	}

	for( i = 1; i < 101; i++)
	{
		in_bam->expected_kmer[i] = ( float)rd_per_gc[i] / ( window_per_gc[i]);
		if( isnanf( in_bam->expected_kmer[i]) || isinff( ( in_bam->expected_kmer[i])) == -1
				|| isinff( ( in_bam->expected_kmer[i])) == 1 )
			in_bam->expected_kmer[i] = 0;
		//fprintf(stderr,"GC = %d - RD=%ld\tWINDOW=%d\tMEAN=%f\n", i, rd_per_gc[i],window_per_gc[i], in_bam->expected_kmer[i]);
	}
}


void read_fastq(parameters* params)
{
	FILE* fp_reads = NULL;
	char seqFile[MAX_SEQ];

	fp_reads = fopen(params->fastq, "r");
	if (fp_reads == NULL){
		printf("Could not open file %s",params->fastq);
		exit(1);
	}
}

