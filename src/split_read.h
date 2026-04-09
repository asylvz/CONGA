/*
 * split_read.h
 *
 *  Created on: Nov 18, 2019
 *      Author: tardis
 */

#ifndef SPLIT_READ_H_
#define SPLIT_READ_H_

#include "common.h"

#define MAX_SR_HIT 50000
#define MAX_MAPPING 100
#define HASH_COUNT 0
#define HASH_BUILD 1

#define HASHKMERLEN 10

extern int *hash_table_count;
extern int *hash_table_iter;

typedef struct lociInRef
{
	int pos; /* pos in the referenceSeqInterest array (keeps the reference genome for each chromosome) */
	struct lociInRef *next;
}lociInRef;

typedef struct posMapSplitRead{
	int posMap;
	char orient;
	int mapq; // the mapq calculate based on number of mappings for this softclip + the distance

	struct posMapSplitRead *next;
}posMapSplitRead;

typedef struct splitRead
{
	char *readName;
	char *chromosome_name;
	int pos;
	char orient;
	int qual;
	int opCount;
	int split_start;
	int read_length;
	//char *split_sequence;
	int avgPhredQualSplitRead; // the average phred quality of the basepairs which have be cut (the soft clip part).

	struct splitRead *next;
	posMapSplitRead *ptrSplitMap; // a linked list of all the positions that the soft clipped part of the read maps in the reference genome (chromosome_name:windowStart-windowEnd)
}splitRead;

int find_split_reads( bam_info* in_bam, parameters* params, bam1_t* bam_alignment, int chr_index);
void map_split_reads( bam_info* in_bam, parameters* params, int chr_index);
void init_hash_count(parameters *params);
void free_hash_table(parameters *params);
void build_hash_table(const char *ref, int len, int hash_size, int mode);
void create_hash_table( parameters *params, int len);

#endif /* SPLIT_READ_H_ */
