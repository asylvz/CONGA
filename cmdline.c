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
	static int load_sonic = 0, no_sr = 0, no_kmer = 0;
	static int do_remap = 0;
	char *min_rd_support = NULL, *min_mapping_qual = NULL, *min_rp_support = NULL;

	static struct option long_options[] = 
	{
			{"rd", required_argument, 0, 'a'},
			{"min-read-length"    , required_argument,	 0, 'b'},
			{"dels"    , required_argument,   0, 'd'},
			{"mq", required_argument, 0, 'e'},
			{"ref"    , required_argument,   0, 'f'},
			{"fastq"    , required_argument,   0, 'g'},
			{"help"   , no_argument,         0, 'h'},
			{"input"  , required_argument,   0, 'i'},
			{"rp", required_argument, 0, 'j'},
			{"min-sv-size"    , required_argument,	 0, 'l'},
			{"mappability"    , required_argument,	 0, 'm'},
			{"sonic-info"    , required_argument,	 0, 'n'},
			{"out"    , required_argument,	 0, 'o'},
			{"sonic"    , required_argument,	 0, 's'},
			{"dups"    , required_argument,   0, 'u'},
			{"version", no_argument,         0, 'v'},
			{"exclude", required_argument, 0,'x'},
			{"no-sr", no_argument, &no_sr, 1},
			{"no-kmer", no_argument, &no_kmer, 1},
			{"first-chrom", required_argument, 0, FIRST_CHROM},
			{"last-chrom", required_argument, 0, LAST_CHROM},
			{0        , 0,                   0,  0 }
	};

	if( argc == 1)
	{
		print_help();
		return 0;
	}

	while( ( o = getopt_long( argc, argv, "hvb:i:f:g:d:r:o:m:c:s:a:e:n:j:k:u:x", long_options, &index)) != -1)
	{
		switch( o)
		{
		case 'a':
			set_str( &( min_rd_support), optarg);
			break;

		case 'b':
			params->min_read_length = atoi( optarg);
			break;

		case 'd':
			set_str( &( params->del_file), optarg);
			break;

		case 'e':
			set_str( &( min_mapping_qual), optarg);
			break;

		case 'f':
			set_str( &( params->ref_genome), optarg);
			break;

		case 'g':
			set_str( &( params->fastq), optarg);
			break;

		case 'h':
			print_help();
			return 0;
			break;

		case 'i':
			set_str( &( params->bam_file), optarg);
			break;

		case 'j':
			set_str( &( min_rp_support), optarg);
			break;

		case 'l':
			params->min_sv_size = atoi( optarg);
			break;

		case 'm':
			set_str( &( params->mappability_file), optarg);
			break;

		case 'n':
			set_str( &( params->sonic_info), optarg);
			break;

		case 'o':
			set_str( &( params->outprefix), optarg);
			break;

		case 's':
			set_str( &( params->sonic_file), optarg);
			load_sonic = 1;
			break;

		case 'u':
			set_str( &( params->dup_file), optarg);
			break;

		case 'v':
			fprintf( stderr, "\n...CONGA....\n");
			fprintf( stderr, "Version %s\n\tLast update: %s, build date: %s\n\n", CONGA_VERSION, CONGA_UPDATE, BUILD_DATE);
			return 0;
			break;

		case 'x':
			set_str( &( params->low_map_regions), optarg);
			break;

		case FIRST_CHROM:
			params->first_chrom = atoi(optarg);
			break;

		case LAST_CHROM:
			params->last_chrom = atoi(optarg);
			break;
		}
	}

	params->no_sr = no_sr;
	params->no_kmer = no_kmer;

	/* check if outprefix is given */
	if( params->outprefix == NULL)
	{
		fprintf( stderr, "[CONGA CMDLINE ERROR] Please enter the output file name prefix using the --out option.\n");
		return EXIT_PARAM_ERROR;
	}

	/* check if --ref   is invoked */
	if( params->ref_genome == NULL)
	{
		fprintf( stderr, "[CONGA CMDLINE ERROR] Please enter reference genome file (FASTA) using the --ref option.\n");
		return EXIT_PARAM_ERROR;
	}


	/* check if --sonic  is invoked */
	if( params->sonic_file == NULL && load_sonic)
	{
		fprintf( stderr, "[CONGA CMDLINE ERROR] Please enter the SONIC file (BED) using the --sonic option.\n");
		return EXIT_PARAM_ERROR;
	}

	if( params->min_sv_size <= 0)
	{
		params->min_sv_size = 1000;
		fprintf( stderr, "Minimum size of an SV is set to %d\n", params->min_sv_size);
	}

	if( params->min_read_length <= 0)
	{
		params->min_read_length = 60;
		fprintf( stderr, "Minimum size of a read is set to %d\n", params->min_read_length);
	}

	if( min_rd_support == NULL)
		params->rd_threshold = 1000;
	else
	{
		params->rd_threshold = atoi(min_rd_support);
		free( min_rd_support);
	}

	if( min_rp_support == NULL)
		params->rp_support = 20;
	else
	{
		params->rp_support = atoi(min_rp_support);
		free( min_rd_support);
	}

	if( min_mapping_qual == NULL)
		params->mq_threshold = 5;
	else
	{
		params->mq_threshold = atoi(min_mapping_qual);
		free( min_mapping_qual);
	}

	if (load_sonic)
		params->load_sonic = load_sonic;

	if ( params->sonic_info == NULL)
		set_str( &(params->sonic_info), params->ref_genome);

	get_working_directory(params);
	fprintf(stderr, "[CONGA INFO] Working directory: %s\n", params->outdir);

	return RETURN_SUCCESS;

}

void print_help( void)
{  
	fprintf( stdout, "\n... CONGA (COpy Number Genotyping in Ancient genomes) ...\n");
	fprintf( stdout, "Version %s\n\tLast update: %s, build date: %s\n\n", CONGA_VERSION, CONGA_UPDATE, BUILD_DATE);
	fprintf( stdout, "\tParameters:\n\n");
	fprintf( stdout, "\t--input [BAM file ]        	: Input file in sorted and indexed BAM format.\n");
	fprintf( stdout, "\t--out   [output prefix]    	: Prefix for the output file names.\n");
	fprintf( stdout, "\t--ref   [reference genome] 	: Reference genome in FASTA format.\n");
	fprintf( stdout, "\t--sonic [sonic file]       	: SONIC file that contains assembly annotations.\n");
	fprintf( stdout, "\t--dels [BED file]          	: Known deletion SVs in BED format\n");
	fprintf( stdout, "\t--dups [BED file]          	: Known duplication SVs in BED format\n");
	fprintf( stdout, "\t--mapability [BED file]     : Mappability file in BED format\n");
	fprintf( stdout, "\t--min-read-length [int]    	: Minimum length of a read to be processed for RP (default: 60 bps)\n");
	fprintf( stdout, "\t--min-sv-size [int]    		: Minimum length of a read to be processed (default: 1000 bps)\n");
	fprintf( stdout, "\t--no-sr                     : Split read mapping is disabled\n");
	fprintf( stdout, "\t--no-kmer                   : k-mer likelihood calculation is disabled (faster)\n");

	fprintf( stdout, "\n\n\tInformation:\n");
	fprintf( stdout, "\t--version                  : Print version and exit.\n");
	fprintf( stdout, "\t--help                     : Print this help screen and exit.\n\n");
}
