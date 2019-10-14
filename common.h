#ifndef __COMMON
#define __COMMON

#include <stdio.h>
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <zlib.h>
#include <stdbool.h>
#include "sonic/sonic.h"

//#define MAIN_DELETION_CLUSTER
#define DELETION 'D'
#define DUPLICATION 'E'

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

typedef struct _params
{
	char* ref_genome; /* path to reference genome - fasta */
	char* outdir;
	char *dup_file;
	char *del_file;
	char* bam_file; /* the actual list that holds all bam file paths after tokenization */
	char* outprefix; /* prefix for the output files */
	int min_sv_size;
	int load_sonic; /*load SONIC file*/
	char *sonic_info; /* SONIC reference information string for building */
	int first_chrom; /*the first chromosome as indexed in the ref to be computer for. 0 by default*/
	int last_chrom; /*the last chromosome as indexed in the ref to be computer for. ref->chrom_count by default*/
	int rd_threshold; /* Threshold is used in RD filtering, calWeight() in vh_setcover.c */
	int mq_threshold; /* Minimum mapping quality */
	char *sonic_file; /* SONIC file name */
	sonic *this_sonic; /* SONIC */
} parameters;

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
int vh_cmprReadNameStr (const void *a, const void *b);
void get_working_directory(parameters *params);

// Memory allocation/tracking functions
void* getMem( size_t size);
void freeMem( void* ptr, size_t size);
double getMemUsage();


#endif
