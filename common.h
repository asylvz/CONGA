#ifndef __COMMON
#define __COMMON

#include <stdio.h>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <zlib.h>
#include <stdbool.h>
#include "sonic/sonic.h"


#define DELETION 'D'
#define DUPLICATION 'E'

#define KMER 33
#define KMERSLIDE 5
#define MINKMERHASHSIZE 20

#define LEFT 'L'
#define RIGHT 'R'
#define NONE 'N'
#define FORWARD 'F'
#define REVERSE 'R'

/* Exit Codes */
#define EXIT_SUCCESS 0
#define EXIT_COMMON 1
#define EXIT_MAXBAMS 2
#define EXIT_PARAM_ERROR 3
#define EXIT_EXTERNAL_PROG_ERROR 4
#define EXIT_FILE_OPEN_ERROR 5
#define EXIT_READGROUP 6
#define EXIT_SONIC 7

/* Return Codes */
#define RETURN_SUCCESS 1
#define RETURN_ERROR 0

#define MAX_BAMS 256
#define MAXLISTBRKPOINTINTR 10000000;

/* Maximum filename length */
#define MAX_LENGTH 1024

// Track memory usage
extern long long memUsage;
extern FILE *logFile; //Defined in tardis.c

/* Gender of the bam file */
enum gender{ MALE, FEMALE};

extern struct SplitsInfo *all_split_reads;

typedef struct SplitsInfo
{
	int size;
	struct SplitRow *head;
	struct SplitRow *tail;

}SplitsInfo;

typedef struct _params
{
	char* ref_genome; /* path to reference genome - fasta */
	char* outdir;
	char *low_map_regions; /*Regions with low mappability is excluded from consideration - a BED input is required */
	char *dup_file;
	char *del_file;
	char* bam_file; /* the actual list that holds all bam file paths after tokenization */
	char* outprefix; /* prefix for the output files */
	int min_sv_size; /* minimum size of known SV loci. Smaller sized SVs than this value are discarded */
	int min_read_length;
	int first_chrom; /*the first chromosome as indexed in the ref to be computer for. 0 by default*/
	int last_chrom; /*the last chromosome as indexed in the ref to be computer for. ref->chrom_count by default*/
	float c_score; /* Likelihood score threshold */
	int mq_threshold; /* Minimum mapping quality */
	int rp_support; /* Minimum number of read-pairs needed to support duplications */
	char *ref_seq; /* reference sequence per chromosome */
	int hash_size; /* size of the hash table for split read mapping */
	char *sonic_file; /* SONIC file name */
	char *mappability_file; /* a file that contains mappability information for various intervals */
	int load_sonic; /*load SONIC file*/
	int no_sr; /* Don't use split-read */
	char *sonic_info; /* SONIC reference information string for building */
	sonic *this_sonic; /* SONIC */
} parameters;

typedef struct _bam_info
{
	int total_read_count_unfiltered; /* total number of reads that are > some mapping quality */
	short* read_depth; /* read depth low qual */
	float* mappability; /* mappability value for each base */
	float mean;
	float expected_read_depth[101];

	htsFile* bam_file; /* file pointer to the BAM file */
	hts_idx_t* bam_file_index;
	hts_itr_t *iter;
	bam_hdr_t* bam_header;

	enum gender sample_gender; /* gender of the sample */
	char* sample_name; /* name of the sample, parsed from SM in the BAM header */
	struct splitRead *listSplitRead;
} bam_info;

typedef struct SplitRow
{
	int locMapLeftEnd;
	int locMapLeftStart;
	char orientationLeft;
	int locMapRightStart;
	int locMapRightEnd;
	int startPosition;
	int endPosition;
	char orientationRight;
	char svType;

	struct SplitRow *next;
} SplitRow;


/* Parameter related TARDIS functions */
void init_params( parameters**);
void print_params( parameters*);

/* FILE opening and error printing functions. For opening regular and BAM/SAM
 files safely */
void print_error( char*);
FILE* safe_fopen( char* path, char* mode);
gzFile safe_fopen_gz( char* path, char* mode);
htsFile* safe_hts_open( char* path, char* mode);

/* String functions */
void set_str( char **target, char *source); /* Even safer than strncpy */
void reverse_string( char* str);

/* Misc. Utility */
int compare_start_pos (const void *a, const void *b);
int compare_size_int( const void* p, const void* q);
int find_chr_index_bam( char* chromosome_name, bam_hdr_t* bam_header);
int max( int x, int y);
int min( int x, int y);
int hammingDistance( char *str1, char *str2, int len);
void get_working_directory(parameters *params);
int is_dna_letter( char base);
int is_proper( int flag);
void get_sample_name(bam_info* in_bam, char* header_text);
char *substring(char *string, int position, int length);
char* reverseComplement( char* str);
long readReferenceSeq( parameters *params, int chr_index);
char* get_refseq( parameters *params, char* chr_name, int start, int end);
char* read_ref( parameters *params, int chr_index);

// Memory allocation/tracking functions
void* getMem( size_t size);
void freeMem( void* ptr, size_t size);
double getMemUsage();


#endif
