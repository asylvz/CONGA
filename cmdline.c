#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cmdline.h"

#include "svdepth.h"

#define FIRST_CHROM 10001
#define LAST_CHROM 10002


int parse_command_line( int argc, char** argv, parameters* params)
{
	int index;
	int o;
	static int load_sonic = 0;
	static int do_remap = 0;
	char *threshold = NULL, *mapping_qual = NULL, *rp_support = NULL;

	static struct option long_options[] = 
	{
			{"input"  , required_argument,   0, 'i'},
			{"ref"    , required_argument,   0, 'f'},
			{"dups"    , required_argument,   0, 'u'},
			{"dels"    , required_argument,   0, 'd'},
			{"help"   , no_argument,         0, 'h'},
			{"version", no_argument,         0, 'v'},
			{"sv-size"    , required_argument,	 0, 'l'},
			{"min-read-length"    , required_argument,	 0, 'b'},
			{"out"    , required_argument,	 0, 'o'},
			{"sonic"    , required_argument,	 0, 's'},
			{"sonic-info"    , required_argument,	 0, 'n'},
			{"first-chrom", required_argument, 0, FIRST_CHROM},
			{"last-chrom", required_argument, 0, LAST_CHROM},
			{"rd", required_argument, 0, 'a'},
			{"rp", required_argument, 0, 'j'},
			{"mq", required_argument, 0, 'e'},
			{0        , 0,                   0,  0 }
	};

	if( argc == 1)
	{
		print_help();
		return 0;
	}

	while( ( o = getopt_long( argc, argv, "hvb:i:f:g:d:r:o:m:c:s:a:e:n:j:k:u", long_options, &index)) != -1)
	{
		switch( o)
		{

		case 'i':
			set_str( &( params->bam_file), optarg);
			break;

		case 'f':
			set_str( &( params->ref_genome), optarg);
			break;

		case 'd':
			set_str( &( params->del_file), optarg);
			break;

		case 'u':
			set_str( &( params->dup_file), optarg);
			break;

		case 's':
			set_str( &( params->sonic_file), optarg);
			load_sonic = 1;
			break;

		case 'n':
			set_str( &( params->sonic_info), optarg);
			break;

		case 'o':
			set_str( &( params->outprefix), optarg);
			break;

		case 'l':
			params->min_sv_size = atoi( optarg);
			break;

		case 'b':
			params->min_read_length = atoi( optarg);
			break;

		case 'h':
			print_help();
			return 0;
			break;

		case FIRST_CHROM:
			params->first_chrom = atoi(optarg);
			break;

		case LAST_CHROM:
			params->last_chrom = atoi(optarg);
			break;

		case 'a':
			set_str( &( threshold), optarg);
			break;

		case 'e':
			set_str( &( mapping_qual), optarg);
			break;

		case 'j':
			set_str( &( rp_support), optarg);
			break;

		case 'v':
			fprintf( stderr, "\n...SvDepth....\n");
			fprintf( stderr, "Version %s\n\tLast update: %s, build date: %s\n\n", SVDEPTH_VERSION, SVDEPTH_UPDATE, BUILD_DATE);
			return 0;
			break;
		}
	}

	/* check if outprefix is given */
	if( params->outprefix == NULL)
	{
		fprintf( stderr, "[SvDepth CMDLINE ERROR] Please enter the output file name prefix using the --out option.\n");
		return EXIT_PARAM_ERROR;
	}

	/* check if --ref   is invoked */
	if( params->ref_genome == NULL)
	{
		fprintf( stderr, "[SvDepth CMDLINE ERROR] Please enter reference genome file (FASTA) using the --ref option.\n");
		return EXIT_PARAM_ERROR;
	}


	/* check if --sonic  is invoked */
	if( params->sonic_file == NULL && load_sonic)
	{
		fprintf( stderr, "[SvDepth CMDLINE ERROR] Please enter the SONIC file (BED) using the --sonic option.\n");
		return EXIT_PARAM_ERROR;
	}

	if( params->min_sv_size <= 0)
	{
		params->min_sv_size = 1000;
		fprintf( stderr, "Minimum size of an SV is set to %d\n", params->min_sv_size);
	}

	if( threshold == NULL)
		params->rd_threshold = 1000;
	else
	{
		params->rd_threshold = atoi(threshold);
		free( threshold);
	}

	if( rp_support == NULL)
		params->rp_support = 2;
	else
	{
		params->rp_support = atoi(rp_support);
		free( threshold);
	}

	if( mapping_qual == NULL)
		params->mq_threshold = 5;
	else
	{
		params->mq_threshold = atoi(mapping_qual);
		free( mapping_qual);
	}

	if (load_sonic)
		params->load_sonic = load_sonic;

	if ( params->sonic_info == NULL)
		set_str( &(params->sonic_info), params->ref_genome);

	get_working_directory(params);
	fprintf(stderr, "[SvDepth INFO] Working directory: %s\n", params->outdir);

	return RETURN_SUCCESS;

}

void print_help( void)
{  
	fprintf( stdout, "\n... SvDepth ...\n");
	fprintf( stdout, "Version %s\n\tLast update: %s, build date: %s\n\n", SVDEPTH_VERSION, SVDEPTH_UPDATE, BUILD_DATE);
	fprintf( stdout, "\tParameters:\n\n");
	fprintf( stdout, "\t--input [BAM file ]        : Input file in sorted and indexed BAM format.\n");
	fprintf( stdout, "\t--out   [output prefix]    : Prefix for the output file names.\n");
	fprintf( stdout, "\t--ref   [reference genome] : Reference genome in FASTA format.\n");
	fprintf( stdout, "\t--sonic [sonic file]       : SONIC file that contains assembly annotations.\n");
	fprintf( stdout, "\t--dels [BED file]          : Known deletion SVs in BED format\n");
	fprintf( stdout, "\t--dups [BED file]          : Known duplication SVs in BED format\n");

	fprintf( stdout, "\n\n\tInformation:\n");
	fprintf( stdout, "\t--version                  : Print version and exit.\n");
	fprintf( stdout, "\t--help                     : Print this help screen and exit.\n\n");
}
