#include "svdepth.h"

#include <stdio.h>
#include <time.h>

#include "bam_data.h"
#include "common.h"
#include "cmdline.h"
#include "sonic/sonic.h"
#include "free.h"
#include "read_distribution.h"

FILE *logFile = NULL;
int total_dels = 0;
int total_dups = 0;

int main( int argc, char** argv)
{
	bam_info* in_bam;
	parameters* params;
	int return_value;
	char* username;
	int i;
	time_t rawtime;
	struct tm * timeinfo;

	time ( &rawtime);
	timeinfo = localtime( &rawtime);

	/* Keeping simple logs in svdepth.log file */
	logFile = safe_fopen ("svdepth.log", "w");
	fprintf( logFile, "#CreationDate=%d.%d.%d\n\n", timeinfo->tm_year+1900, timeinfo->tm_mon+1, timeinfo->tm_mday);

	/* Set program parameters */
	init_params( &params);

	/* Parse command line arguments */	
	return_value = parse_command_line( argc, argv, params);

	print_params( params);

	/* Load SONIC */
	params->this_sonic = sonic_load(params->sonic_file);

	if (params->last_chrom < params->first_chrom){
		params->last_chrom = params->this_sonic->number_of_chromosomes - 1;
	}

	/* Read BAM files */
	in_bam = ( bam_info*) getMem( sizeof( bam_info));
	in_bam->sample_name = NULL;

	htsFile* bam_file = safe_hts_open( params->bam_file, "r");

	/* Read in BAM header information */
	bam_hdr_t* bam_header = bam_hdr_read( ( bam_file->fp).bgzf);

	/* Extract the Sample Name from the header text */
	get_sample_name( in_bam, bam_header->text);

	return_value = hts_close( bam_file);

	//Find the depth of SVs
	read_bam(in_bam, params);

	//Free the data structures
	free_DS(in_bam, params);

	username = ( char*) getMem( MAX_SEQ * sizeof( char));
	getlogin_r( username, (MAX_SEQ - 1));
	fprintf( stderr, "\nThank you %s. I found %d DELs and %d DUPs. Hope to see you again...\n", username, total_dels, total_dups);

	fclose( logFile);

	free( username);
	return EXIT_SUCCESS;
}
