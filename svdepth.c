#include "svdepth.h"
#include <stdio.h>
#include <time.h>

#include "bam_data.h"
#include "cmdline.h"
#include "sonic/sonic.h"
#include "free.h"
#include "read_distribution.h"
#include "split_read.h"
#include "mhash.h"
#include "kmer.h"

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

	/* Keeping simple logs in conga.log file */
	logFile = safe_fopen ("conga.log", "w");
	fprintf( logFile, "#CreationDate=%d.%d.%d\n\n", timeinfo->tm_year+1900, timeinfo->tm_mon+1, timeinfo->tm_mday);
/*
	char x[100] = ">;7>?>2?>?>?@@?@A;0?(B@D@@C@??>9+@B()+B;D>2A9@<2'>@(9='(7'=D;,0<A,=%081-*5*";

	int sum = 0;
	for(i = 0; i < strlen(x); i++)
	{
		int b = x[i] - 33;
		fprintf(stderr,"%c - %d - %lf\n", x[i], b, pow(10,(double) (-b)/10));
		sum += x[i] - 33;
	}

	double a = pow(10,(double) (-((double) sum/i))/10);
	fprintf(stderr,"%d\t%d\t%f\n",sum, sum/i, a);*/


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
	in_bam->listSplitRead = NULL;

	if(!params->no_sr)
		init_hash_count( params);

	//Find the depth of SVs
	read_bam(in_bam, params);

	//Free the data structures
	free_DS(in_bam, params);

	username = ( char*) getMem( MAX_SEQ * sizeof( char));
	getlogin_r( username, (MAX_SEQ - 1));
	fprintf( stderr, "\nThank you %s. I found %d DELs and %d DUPs. Hope to see you again...\n", username, total_dels, total_dups);

	fclose(logFile);

	free(username);
	return EXIT_SUCCESS;
}
