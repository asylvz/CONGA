#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cmdline.h"
#include "svdepth.h"

#define FIRST_CHROM 10001
#define LAST_CHROM 10002


int parse_cmd_line( int argc, char** argv, parameters* params)
{
	int index;
	int o;
	static int load_sonic = 0, no_sr = 1;
	static int do_remap = 0;
	char *min_rd_support = NULL, *min_mapping_qual = NULL, *min_rp_support = NULL;

	static struct option long_options[] = 
	{
			{"c-score", required_argument, 0, 'a'},
			{"min-read-length"    , required_argument,	 0, 'b'},
			{"dels"    , required_argument,   0, 'd'},
			{"min-mapq", required_argument, 0, 'e'},
			{"ref"    , required_argument,   0, 'f'},
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
			{"first-chr", required_argument, 0, FIRST_CHROM},
			{"last-chr", required_argument, 0, LAST_CHROM},
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

		case 'h':
			print_help();
			return EXIT_SUCCESS;
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
			//fprintf( stderr, "\n...CONGA....\n");
			fprintf( stderr, "\n\tCONGA Version %s\n\tLast update: %s, build date: %s\n", CONGA_VERSION, CONGA_UPDATE, BUILD_DATE);
			fprintf(stderr,"\tFor more information, check https://github.com/asylvz/CONGA\n\n");
			return EXIT_SUCCESS;
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
		params->c_score = 0.5;
	else
	{
		params->c_score = atof(min_rd_support);
		free( min_rd_support);
	}

	if( min_rp_support == NULL)
	{
		params->rp_support = 10;
		params->no_sr = 1;
	}
	else
	{
		params->rp_support = atoi(min_rp_support);
		params->no_sr = 0;
		free( min_rd_support);
	}

	if( min_mapping_qual == NULL)
		params->mq_threshold = -1;
	else
	{
		params->mq_threshold = atoi(min_mapping_qual);
		free( min_mapping_qual);
	}

	if (load_sonic)
		params->load_sonic = load_sonic;

	if ( params->sonic_info == NULL)
		set_str( &(params->sonic_info), params->ref_genome);

	//params->no_sr = no_sr;
	get_working_directory(params);
	fprintf(stderr, "[CONGA INFO] Working directory: %s\n", params->outdir);

	return RETURN_SUCCESS;

}

void print_help( void)
{  
	fprintf( stdout, "\n\t... CONGA (COpy Number variation Genotyping in Ancient genomes) ...\n\n");
	fprintf( stdout, "\tVersion %s\n\tLast update: %s, build date: %s\n\n", CONGA_VERSION, CONGA_UPDATE, BUILD_DATE);
	fprintf( stdout, "\tParameters:\n");
	fprintf( stdout, "\t--input [BAM file]        : Input file in sorted and indexed BAM format (required).\n");
	fprintf( stdout, "\t--out [output prefix]     : Prefix for the output file names (required).\n");
	fprintf( stdout, "\t--ref [reference genome]  : Reference genome in FASTA format (required).\n");
	fprintf( stdout, "\t--sonic [sonic file]      : SONIC file that contains assembly annotations (required).\n");
	fprintf( stdout, "\t--dels [BED file]         : Known deletion SVs in BED format\n");
	fprintf( stdout, "\t--dups [BED file]         : Known duplication SVs in BED format\n");
	fprintf( stdout, "\t--first-chr [chr index]   : The index of the first chromosome for genotyping in your BAM\n");
	fprintf( stdout, "\t--last-chr [chr index]    : The index of the last chromosome for genotyping in your BAM\n");
	fprintf( stdout, "\t--mapability [BED file]   : Mappability file in BED format\n");
	fprintf( stdout, "\t--min-read-length [INT]   : Minimum length of a read to be processed for RP (default: 60 bps)\n");
	fprintf( stdout, "\t--min-sv-size [INT]       : Minimum length of a CNV (default: 1000 bps)\n");
	fprintf( stdout, "\t--min-mapq [INT]          : Minimum mapping quality filter for reads (default: no-filter)\n");
	fprintf( stdout, "\t--c-score [FLOAT]         : Minimum c-score to filter variants (More conservative with lower values, default: 0.5).\n");
	fprintf( stdout, "\t--rp [INT]                : Enable split-read and set minimum read-pair support for a duplication (Suggested for >5x only).");

	fprintf( stdout, "\n\n\tInformation:\n");
	fprintf( stdout, "\t--version                 : Print version and exit.\n");
	fprintf( stdout, "\t--help                    : Print this help screen and exit.\n\n");
	fprintf(stderr,"\n\t* For more information, please consult https://github.com/asylvz/CONGA\n\n");
}
