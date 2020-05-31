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
#define KMERWINDOWSIZE 100
#define KMERWINDOWSLIDE 10
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

int read_kmer_jellyfish( parameters *params);
void free_hash_table_kmer(parameters *params);
char* read_ref_seq( parameters *params, char* chr_name, int start, int end);
int query_jellyfish(char* str, int count, char type);
int is_kmer_valid_likelihood (char *str);
int calculate_kmers_jellyfish(parameters *params, char* chr_name, int count, int total_variant, char type);
void read_kmer_seqs( parameters *params, int base_count, int mode);
void init_hash_count_kmer(parameters *params);
int kmer_count_interval(parameters *params, int start, int end);
void init_kmer_per_chr( bam_info* in_bam, parameters* param, int chr_index);
void calc_expected_kmer(bam_info *in_bam, parameters *params, int chr_index);
void calc_kmer_counts(bam_info *in_bam, parameters *params, int chr_index);

#endif /* KMER_H_ */
