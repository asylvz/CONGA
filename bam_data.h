#ifndef __PROCESSBAM
#define __PROCESSBAM

/* htslib headers */
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <stdbool.h>
#include "common.h"

/* Maximum sequence/quality length */
#define MAX_SEQ 1000

/* Gender of the bam file */
enum gender{ MALE, FEMALE};

typedef struct _bam_info
{
	int read_count; /* total number of reads in this library */
	short* read_depth_per_chr; /* read depth */
	float mean;
	float mean_rd_per_gc[101]; /* GC percentages, i.e., GC[13]=323 means 343 windows have GC of 13% */

	htsFile* bam_file; /* file pointer to the BAM file */
	hts_idx_t* bam_file_index;
	hts_itr_t *iter;
	bam_hdr_t* bam_header;

	enum gender sample_gender; /* gender of the sample */
	char* sample_name; /* name of the sample, parsed from SM in the BAM header */
} bam_info;


/* Function Prototypes */
void load_bam( bam_info* in_bam, char* path);
void print_bam( bam_info* in_bam);
void read_bam( bam_info* in_bam, parameters *params);


/* BAM Utility functions */
void get_sample_name( bam_info* in_bam, char* header_text);


#endif
