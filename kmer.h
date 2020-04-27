/*
 * kmer.h
 *
 *  Created on: Apr 8, 2020
 *      Author: tardis
 */

#ifndef KMER_H_
#define KMER_H_

#include "svs.h"

extern svs* all_svs_del;
extern svs* all_svs_dup;

#define JELLYWINDOWSIZE 1200
#define KMERWINDOWSIZE 1000
#define KMERWINDOWSLIDE 1000
#define MAX_SEQ 1000
#define HASH_COUNT 0
#define HASH_BUILD 1
#define MAX_KHASH_HIT 5


typedef struct HashInfo
{
	unsigned int hash2;
	short freq;

	struct HashInfo *next;
}HashInfo;

void read_kmer_jellyfish( parameters *params);
void free_hash_table_kmer(parameters *params);
char* read_ref_seq( parameters *params, char* chr_name, int start, int end);
int query_jellyfish(char* str, int count, char type);
int is_kmer_valid_likelihood (char *str);
//void init_hash_count_kmer(parameters *params);
int calculate_kmers_jellyfish(parameters *params, char* chr_name, int count, int total_variant, char type);
void read_kmer_seqs( parameters *params, int base_count, int mode);
void init_hash_count_kmer(parameters *params);
int calculate_kmers(parameters *params, char* chr_name, int start, int end);
void calc_mean_kmer( parameters *params, char* chr_name, int chr_index, float *mean_kmer_per_gc);

#endif /* KMER_H_ */
