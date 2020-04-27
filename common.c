#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>

#include <htslib/sam.h>
#include <htslib/hts.h>
#include "common.h"
#include "svs.h"

// Track memory usage
long long memUsage = 0;

void init_params( parameters** params)
{
	int i;

	/* initialize parameters */
	*params = ( parameters*) getMem( sizeof( parameters));
	( *params)->ref_genome = NULL;
	( *params)->dup_file = NULL;
	( *params)->del_file = NULL;
	( *params)->low_map_regions = NULL;
	( *params)->sonic_file = NULL;
	( *params)->this_sonic = NULL;
	( *params)->outprefix = NULL;
	( *params)->outdir = NULL;
	( *params)->bam_file = NULL;
	( *params)->bam_file = ( char*) getMem( sizeof( char));
	( *params)->min_sv_size = 0;
	( *params)->rp_support = 0;
	( *params)->first_chrom = 0;
	( *params)->last_chrom = -1;
	( *params)->load_sonic = 0;
	( *params)->sonic_info = NULL;
	( *params)->ref_seq = NULL;
	//( *params)->kmer_seq = NULL;
	( *params)->hash_size = 0;
	( *params)->min_read_length = 0;
}


void get_working_directory(parameters *params)
{
	char *directory;
	char *prefix;
	int i;

	directory = strrchr(params->outprefix, '/');
	prefix = NULL;
	if (directory != NULL)
		set_str( &prefix, directory + 1);
	if ( prefix != NULL)
		fprintf (stderr, "prefix: %s\n", prefix);

	if (directory == NULL)
	{
		set_str(&(params->outdir), "");
		return;
	}

	i = 0;
	while (params->outprefix+i != directory)
		i++;

	params->outdir = (char *) getMem((i+2) * sizeof(char));
	memcpy(params->outdir, params->outprefix, i*sizeof(char));
	params->outdir[i]='/';
	params->outdir[i+1] = 0;

	if ( prefix != NULL)
		fprintf (stderr, "prefix2: %s\n", prefix);

	if (prefix!=NULL)
	{
		free( params->outprefix);
		params->outprefix = NULL;
		set_str( &(params->outprefix), prefix);
		free (prefix);
	}
}


void print_params( parameters* params)
{
	int i;
	printf("\n");

	printf( "%-30s%s\n","BAM input:",params->bam_file);
	fprintf( logFile,"%-30s%s\n","BAM input:",params->bam_file);
	printf( "%-30s%s\n","Reference genome:", params->ref_genome);
	printf( "%-30s%s\n","SONIC file:", params->sonic_file);

	fprintf( logFile, "%-30s%s\n","Reference genome:", params->ref_genome);
	fprintf( logFile, "%-30s%s\n","SONIC file:", params->sonic_file);
	fprintf( logFile, "%-30s%d\n","First chrom:", params->first_chrom);
	fprintf( logFile, "%-30s%d\n","Last chrom:", params->last_chrom);
}

void print_error( char* msg)
{
	/* print error message and exit */
	fprintf( stderr, "\n%s\n", msg);
	fprintf( stderr, "Invoke parameter -h for help.\n");
	exit( EXIT_COMMON);
}


FILE* safe_fopen( char* path, char* mode)
{
	/* Safe file open. Try to open a file; exit if file does not exist */
	FILE* file;
	char err[500];

	file = fopen( path, mode);  
	if( !file)
	{
		sprintf( err, "[CONGA INPUT ERROR] Unable to open file %s in %s mode.", path, mode[0]=='w' ? "write" : "read");
		print_error( err);

	}
	return file;
}

gzFile safe_fopen_gz( char* path, char* mode)
{
	/* Safe file open. Try to open a file; exit if file does not exist */
	gzFile file;
	char err[500];

	file = gzopen( path, mode);  
	if( !file)
	{
		sprintf( err, "[CONGA INPUT ERROR] Unable to open file %s in %s mode.", path, mode[0]=='w' ? "write" : "read");
		print_error( err);		
	}
	return file;
}

htsFile* safe_hts_open( char* path, char* mode)
{
	htsFile* bam_file;
	char err[500];

	bam_file = hts_open( path, mode);
	if( !bam_file)
	{
		sprintf( err, "[CONGA INPUT ERROR] Unable to open file %s in %s mode.", path, mode[0]=='w' ? "write" : "read");
		print_error( err);
	}

	return bam_file;
}



/* Even safer than strncpy as it dynamically allocates space for the string if
 there hasn't been already */
void set_str( char** target, char* source)
{
	if( *target != NULL)
	{
		free( ( *target));
	}

	if (source != NULL)
	{
		( *target) = ( char*) getMem( sizeof( char) * ( strlen( source) + 1));
		strncpy( ( *target), source, ( strlen( source) + 1));
	}
	else
	{
		(*target) = NULL;
	}
}


/* Reverse a given string */
void reverse_string( char* str)
{
	int i;
	char swap;
	int len = strlen( str);

	for( i = 0; i < len / 2; i++)
	{
		swap = str[i];
		str[i] = str[len - i - 1];
		str[len - i - 1] = swap;
	}
}


int compare_start_pos (const void *a, const void *b)
{
	int num1 = ((svs*)a)->start;
	int num2 = ((svs*)b)->start;

	int num1_end = ((svs*)a)->end;
	int num2_end = ((svs*)b)->end;

	if(num1 > num2)
		return 1;
	else if(num1 == num2)
		return num1_end - num2_end;
	else
		return -1;

	//return num1 - num2;
}

int compare_size_int( const void* p, const void* q)
{
	int i = *( const int*) p;
	int j = *( const int*) q;

	if( i < j)
	{
		return -1;
	}
	else if( i == j)
	{
		return 0;
	}
	else
	{
		return 1;
	}
}

void* getMem( size_t size)
{
	void* ret;

	ret = malloc( size);
	if( ret == NULL)
	{
		fprintf( stderr, "Cannot allocate memory. Currently addressed memory = %0.2f MB, requested memory = %0.2f MB.\nCheck the available main memory.\n", getMemUsage(), ( float) ( size / 1048576.0));
		exit( 0);
	}

	memUsage = memUsage + size;
	return ret;
}

void freeMem( void* ptr, size_t size)
{
	memUsage = memUsage - size;
	free( ptr);
}

double getMemUsage()
{
	return memUsage / 1048576.0;
}

int max( int x, int y)
{
	if( x < y)
		return y;
	else
		return x;
}

int min( int x, int y)
{
	if( x < y)
		return x;
	else
		return y;
}

int hammingDistance( char *str1, char *str2, int len)
{
	int dist = 0, i;
	for( i = 0; i < len; i++)
	{
		if( str1[i] != str2[i])
			dist++;
	}
	return dist;
}

int find_chr_index_bam( char* chromosome_name, bam_hdr_t* bam_header)
{
	int i, len;
	char *tmp;

	for(i = 0; i < bam_header->n_targets; i++)
	{
		if( strcmp( chromosome_name, bam_header->target_name[i]) == 0)
			return i;
	}
	return -1;
}

int vh_cmprReadNameStr (const void *a, const void *b)
{
	return strcmp (*(char **) a, *(char **) b);
}

/* check if ACGT */
int is_dna_letter( char base)
{
	if (base == 'A')
		return 1;
	if (base == 'C')
		return 1;
	if (base == 'G')
		return 1;
	if (base == 'T')
		return 1;
	return 0;
}

int is_proper( int flag)
{
	if ( (flag & BAM_FSECONDARY) == 0 && (flag & BAM_FSUPPLEMENTARY) == 0 && (flag & BAM_FDUP) == 0 && (flag & BAM_FQCFAIL) == 0)
		return 1;

	return 0;
}

void get_sample_name(bam_info* in_bam, char* header_text)
{
	char *tmp_header = NULL;
	set_str( &( tmp_header), header_text);
	char* p = strtok( tmp_header, "\t\n");
	char sample_name_buffer[1024];

	while( p != NULL)
	{
		/* If the current token has "SM" as the first two characters,
			we have found our Sample Name */
		if( p[0] == 'S' && p[1] == 'M')
		{
			/* Get the Sample Name */
			strncpy( sample_name_buffer, p + 3, strlen( p) - 3);

			/* Add the NULL terminator */
			sample_name_buffer[strlen( p) - 3] = '\0';

			/* Exit loop */
			break;
		}
		p = strtok( NULL, "\t\n");
	}

	set_str( &( in_bam->sample_name), sample_name_buffer);
	free( tmp_header);
}

char *substring(char *string, int position, int length)
{
	char *pointer;
	int c;

	pointer = (char*) getMem(sizeof(char) * length + 1);

	if (pointer == NULL)
	{
		printf("Unable to allocate memory.\n");
		exit(1);
	}
	//Adds one space character at the beginning for jellyfish -don't forget-

	//*(pointer) = ' ';
	for (c = 0 ; c < length ; c++)
	{
		*(pointer + c) = *(string + position - 1);
		string++;
	}

	*(pointer + c) = '\0';

	return pointer;
}

/* Return the complement of a base */
char complement_char( char base)
{
	switch( base)
	{
	case 'A':
		return 'T';
		break;
	case 'C':
		return 'G';
		break;
	case 'G':
		return 'C';
		break;
	case 'T':
		return 'A';
		break;
	default:
		return 'N';
		break;
	}
	return 'X';
}


char* reverseComplement( char* str)
{
	int i;
	char* str2 = NULL;
	char tmp;

	set_str( &str2, str);
	reverse_string(str2);

	for(i = 0; i < strlen(str2); i++)
	{
		tmp = complement_char(str2[i]);
		str2[i] = tmp;
	}
	return str2;
}
