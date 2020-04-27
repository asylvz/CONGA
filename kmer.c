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
	int hash_size = params->hash_size_kmer;
	HashInfo *tmp, *tmp_next;
	for(i = 0; i < hash_size; i++)
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
}

int is_kmer_valid_likelihood (char *str)
{

	int i, l;
	l = (int) strlen(str);

	if (l != KMER)
		return 0;

	for (i = 0; i < KMER; i++)
		if (str[i] != 'A' && str[i] != 'C' && str[i] != 'G' && str[i] != 'T')
			return 0;

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
void read_kmer_jellyfish( parameters *params)
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

	hash_table_kmer = (HashInfo **) getMem (hash_size * sizeof (HashInfo*));

	for(i = 0; i < hash_size; i++)
		hash_table_kmer[i] = NULL;

	i = 0;
	while(!feof(fp_reads))
	{
		if(i != 0 && (i % 2) == 0)
		{
			MurmurHash3_x86_128(kmer_seq, strlen(kmer_seq),42, hash);
			hash_val1 = (unsigned) hash[0] % hash_size;

			MurmurHash3_x86_128(kmer_seq, strlen(kmer_seq),11, hash);
			hash_val2 = (unsigned) hash[0] % hash_size;

			//if(hash_val1 == hash_val2)
				//fprintf(stderr,"PROBLEMM IN HASH - Same values (kmer.c)\n");

			newEl = (HashInfo *) getMem(sizeof(HashInfo));
			newEl->freq = kmer_freq;
			newEl->hash2 = hash_val2;

			//fprintf(stderr,"%s\t%ld\t%u\t%u\t%d\n", kmer_seq, strlen(kmer_seq), hash_val1, hash_val2, newEl->freq);
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

	/*
	for(i = 0; i < hash_size; i++)
	{
		HashInfo *tmp = hash_table_kmer[i];
		while(tmp != NULL)
		{
			fprintf(stderr,"%u\t %d\n",tmp->hash2, tmp->freq);
			tmp = tmp->next;
		}
	}*/
}



/* Read the reference genome
void read_kmer_seqs( parameters *params, int base_count, int mode)
{
	int i, j;
	unsigned int hash_val;
	FILE* fp_reads = NULL;
	char seqFile[MAX_SEQ];
	int mask, len_limit;
	char seed[KMER + 1];
	int hash[4];
	int *hash_table_kmer_iter;

	if (mode == HASH_COUNT)
	{
		hash_table_kmer_count = (int *) getMem( params->hash_size_kmer * sizeof (int));

		sprintf( seqFile, "%s%s_seqs.fa", params->outdir, params->outprefix);
		params->kmer_seq = (char *) getMem( (base_count + 1) * sizeof( char));

		memset (hash_table_kmer_count, 0, params->hash_size_kmer * sizeof (int));

		fp_reads = fopen(seqFile, "r");
		if (fp_reads == NULL){
			printf("Could not open file %s",seqFile);
			exit(1);
		}
		len_limit = fread(params->kmer_seq, sizeof(char), base_count, fp_reads);
		fclose(fp_reads);
		params->kmer_seq[base_count] = '\0';
	}

	else
	{
		hash_table_kmer_iter = (int *) getMem( params->hash_size_kmer * sizeof (int));
		hash_table_kmer_array = (HashInfo **) getMem (params->hash_size_kmer * sizeof (HashInfo *));

		for (i = 0; i < params->hash_size_kmer; i++)
		{
			//if( hash_table_kmer_count[i] > 100)
			//fprintf(stderr,"%d\n", hash_table_kmer_count[i]);
			if (hash_table_kmer_count[i] != 0 && hash_table_kmer_count[i] < MAX_KHASH_HIT)
			{
				//fprintf(stderr,"%d\n", hash_table_kmer_count[i]);
				hash_table_kmer_array[i] = (HashInfo *) getMem (hash_table_kmer_count[i] * sizeof (HashInfo));

				for(j = 0; j < hash_table_kmer_count[i]; j++)
					hash_table_kmer_array[i][j].freq = 0;
			}
			else
			{
				hash_table_kmer_count[i] = 0;
				hash_table_kmer_array[i] = NULL;

			}
		}
		memset (hash_table_kmer_iter, 0, params->hash_size_kmer * sizeof (int));
	}
	//fprintf(stderr,"Here\n");
	i = 0;
	while (i < base_count)
	{
		strncpy (seed, params->kmer_seq + i, KMER);
		seed[KMER] = '\0';

		MurmurHash3_x86_128(seed, strlen(seed),42, hash);
		hash_val = (unsigned) hash[0] % params->hash_size_kmer;
		//fprintf(stderr,"%s - %u - lim:%ld\n", seed, hash_val, params->hash_size_kmer);

		if(hash_val < 0)
		{
			printf("ERRRRROOR HASH VALUE\n");
			exit(1);
		}
		if (mode == HASH_COUNT)
		{
			(hash_table_kmer_count[hash_val])++;
		}
		else if (hash_table_kmer_count[hash_val] != 0)
		{

			hash_table_kmer_array[hash_val][hash_table_kmer_iter[hash_val]].pos = i;
			hash_table_kmer_array[hash_val][hash_table_kmer_iter[hash_val]].freq++;
			(hash_table_kmer_iter[hash_val])++;
		}

		i+=10; //shift
	}
	if (mode != HASH_COUNT)
		free(hash_table_kmer_iter);
}*/


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

int query_jellyfish(char* str, int count, char type)
{
	int k = 0, z, freq = 0, freq_sum = 0;
	FILE *fd;
	fd = popen(str, "r");
	char* ptr_tok = NULL;
	char buffer[256];
	size_t chread;

	if (!fd) return -1;

	while ((chread = fread(buffer, 1, sizeof(buffer), fd)) != 0)
	{
		k++;
		if(k == JELLYWINDOWSIZE)
			break;
		//fprintf(stderr,"Line: %s\n\n", comout);
		ptr_tok = strtok(buffer, " ");
		z = 0;
		freq = 0;
		while (ptr_tok != NULL)
		{
			if(z == 1)
			{
				//fprintf(stderr,"String: %s\t", ptr_tok);
				freq = atoi(ptr_tok);
				break;
			}
			ptr_tok = strtok(NULL, " ");
			z++;
		}
		freq_sum += freq;
	}
	pclose(fd);
	if(type == DELETION)
		all_svs_del[count].k_mer += freq_sum;
	else
		all_svs_dup[count].k_mer += freq_sum;

	return 1;
}
int calculate_kmers_jellyfish(parameters *params, char* chr_name, int count, int total_variant, char type)
{
	int variation_length, i, j = 0, freq;
	char* seq = NULL, *seq_tmp = NULL;
	char cmd_jelly[(JELLYWINDOWSIZE * (KMER + 2))];

	if(type == DELETION)
	{
		variation_length = all_svs_del[count].end - all_svs_del[count].start + 1;
		seq = read_ref_seq(params, chr_name, all_svs_del[count].start, all_svs_del[count].end);
	}
	else
	{
		variation_length = all_svs_dup[count].end - all_svs_dup[count].start + 1;
		seq = read_ref_seq(params, chr_name, all_svs_dup[count].start, all_svs_dup[count].end);
	}

	fprintf(stderr,"%d of %d (size = %d)\n", count, total_variant, variation_length);

	for(i = 0; i < variation_length - KMER; i += 10)
	{
		if(j > JELLYWINDOWSIZE || i == 0)
		{
			if(j > JELLYWINDOWSIZE)
			{
				//fprintf(stderr,"%s\n\n", cmd_jelly);
				query_jellyfish(cmd_jelly, count, type);
			}
			cmd_jelly[0] = '\0';
			sprintf(cmd_jelly, "jellyfish-2.3.0/bin/jellyfish query mer_counts.jf");
			j = 0;
		}
		else
		{
			seq_tmp = substring(seq, i, KMER);
			//fprintf(stderr,"%d - %s\n", j, seq_tmp);
			if(is_kmer_valid_likelihood(seq_tmp) != 0)
			{
				strncat(cmd_jelly, seq_tmp, KMER + 1);
				j++;
			}
			else
			{
				//fprintf(stderr,"ERROR SIZE=%d\t%s",(int) strlen(seq_tmp), seq_tmp);
				continue;
			}
		}
	}
	if(seq != NULL)
	{
		free(seq);
		seq = NULL;
	}
	return 1;
}


int calculate_kmers(parameters *params, char* chr_name, int start, int end)
{
	int variation_length, i, j;
	char* seq = NULL, *seq_tmp = NULL, *seq_tmp_rev = NULL;
	int hash[4];
	unsigned int hash_val1, hash_val2;
	HashInfo *tmp;
	int k_mer_count = 0;


	variation_length = end - start + 1;
	seq = read_ref_seq(params, chr_name, start, end);


	//fprintf(stderr,"%d of %d (size = %d)\n", count, total_variant, variation_length);

	for(i = 0; i < variation_length - KMER; i+=10)
	{
		for(j = 0; j < 2; j++)
		{
			if(j == 0)
			{
				seq_tmp = substring(seq, i, KMER);
				if(seq_tmp == NULL)
					continue;

				if(!is_kmer_valid_likelihood(seq_tmp))
					continue;

				MurmurHash3_x86_128(seq_tmp, strlen(seq_tmp), 42, hash);
				hash_val1 = (unsigned) hash[0] % params->hash_size_kmer;

				MurmurHash3_x86_128(seq_tmp, strlen(seq_tmp), 11, hash);
				hash_val2 = (unsigned) hash[0] % params->hash_size_kmer;
			}
			else
			{
				seq_tmp_rev = reverseComplement(seq_tmp);

				if(seq_tmp_rev == NULL)
					continue;

				if(!is_kmer_valid_likelihood(seq_tmp_rev))
					continue;

				MurmurHash3_x86_128(seq_tmp_rev, strlen(seq_tmp_rev), 42, hash);
				hash_val1 = (unsigned) hash[0] % params->hash_size_kmer;

				MurmurHash3_x86_128(seq_tmp_rev, strlen(seq_tmp_rev), 11, hash);
				hash_val2 = (unsigned) hash[0] % params->hash_size_kmer;

				free(seq_tmp);
				seq_tmp = NULL;

				free(seq_tmp_rev);
				seq_tmp_rev = NULL;
			}

			tmp = hash_table_kmer[hash_val1];
			while(tmp != NULL)
			{
				//fprintf(stderr,"%d\t%u\n", tmp->freq,tmp->hash2);
				if(tmp->hash2 == hash_val2)
				{
					k_mer_count += tmp->freq;
					break;
				}
				tmp = tmp->next;
			}
		}
	}
	if(seq != NULL)
	{
		free(seq);
		seq = NULL;
	}
	return k_mer_count;
}

void calc_mean_kmer( parameters *params, char* chr_name, int chr_index, float *mean_kmer_per_gc)
{
	int i, gc_val = -1, window_per_gc[101], end;
	long rd_per_gc[101];

	/* Calculate mu_GC values */
	for( i = 0; i < 101; i++)
	{
		rd_per_gc[i] = 0;
		window_per_gc[i] = 0;
	}

	for( i = 0; i < params->this_sonic->chromosome_lengths[chr_index]; i += 100)
	{
		if((i + KMERWINDOWSLIDE) < params->this_sonic->chromosome_lengths[chr_index])
			end = i + KMERWINDOWSLIDE;
		else
			end = params->this_sonic->chromosome_lengths[chr_index];

		gc_val = (int) round (sonic_get_gc_content(params->this_sonic, params->this_sonic->chromosome_names[chr_index], i, end));
		rd_per_gc[gc_val] += (long) calculate_kmers(params, chr_name, i, end);
		window_per_gc[gc_val]++;
	}

	mean_kmer_per_gc[0] = 0.0;
	for( i = 1; i < 101; i++)
	{
		mean_kmer_per_gc[i] = ( float)rd_per_gc[i] / ( window_per_gc[i]);
		if( isnanf( mean_kmer_per_gc[i]) || isinff( ( mean_kmer_per_gc[i])) == -1
				|| isinff( ( mean_kmer_per_gc[i])) == 1 )
			mean_kmer_per_gc[i] = 0;
		fprintf(stderr,"GC = %d - RD=%ld\tWINDOW=%d\tMEAN=%f\n", i, rd_per_gc[i],window_per_gc[i], mean_kmer_per_gc[i]);
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

