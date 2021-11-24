CONGA (COpy Number variation Genotyping in Ancient genomes) 
======

CONGA is a genotyping algorithm for Copy Number Variations (large deletions and duplications) in ancient genomes. It is tailored for calling homozygous and heterozygous CNV genotypes at low depths of coverage using read-depth and read-pair information from a BAM file with Illumina short reads.

CONGA is developed and tested using Linux operating system... 

Please feel free to send me an e-mail (asoylev@gmail.com), or better yet open an issue for your questions.


Requirements
============

 * zlib   	(http://www.zlib.net)
 * htslib	(included as submodule; http://htslib.org/)
 * sonic  	(included as submodule; https://github.com/calkan/sonic)

Fetching CONGA
===============

	git clone https://github.com/asylvz/CONGA --recursive

Compilation
===========

Type:

	make libs
	make

Running CONGA
===========================

	conga -i myinput.bam --ref human_g1k_v37.fasta --sonic human_g1k_v37.sonic  \
		--dels known_dels.bed --dups known_dups.bed --out myoutput


SONIC file (annotations container)
==================================

SONIC files for some human and mouse genome reference versions are available at external repo: https://github.com/BilkentCompGen/sonic-prebuilt

 * human_g1k_v37.sonic: SONIC file for Human Reference Genome GRCh37 (1000 Genomes Project version)
	* Also download the reference genome at: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz. 
 * ucsc_hg19.sonic: SONIC file for the human reference genome, UCSC version build hg19.
	* Also download the reference genome at: http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz. Deflate the tar archive and concatenate all chromosomes into a single FASTA file.
 * ucsc_hg38.sonic: SONIC file for the human reference genome build 38.
	* Also download the reference genome at: http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz. Deflate the tar archive and concatenate all chromosomes into a single FASTA file.

Building the SONIC file
=======================

If you are working with a different reference genome please refer to the SONIC development repository: https://github.com/calkan/sonic/

The README.md file includes documentation on how to obtain the necessary files for different genomes from the UCSC Genome Browser.

Example Genotype file
=====================
	1	668630		850204
	1	963826		974172
	1	1171539		1179729
	1	1249799		1265722
	1	2374226		2379823
	...

* The columns are "Chromosome Name" (TAB) "Start Position of a CNV" (TAB) "End Postion of a CNV"
* This file should be seperate for duplications and deletions if both are to be genotyped.


All parameters
==============

	--input 		[BAM file]         : Input files in sorted and indexed BAM format. (required)
	--out   		[output prefix]    : Prefix for the output file names. (required)
	--ref   		[reference genome] : Reference genome in FASTA format. (required)
	--sonic 		[sonic file]       : SONIC file that contains assembly annotations. (required)
	--dels          	[bed file]         : Known deletion SVs in bed format\n");
	--dups          	[bed file]         : Known duplication SVs in bed format\n");
	--mappability   	[bed file]         : Mappability file in BED format
	--first-chr     	[chromosome index] : The index of the first chromosome for genotyping in your BAM.
	--last-chr      	[chromosome index] : The index of the last chromosome for genotyping in your BAM.
	--min-read-length	[integer]	   : Minimum length of a read to be processed for RP (default: 60 bps)
	--min-sv-size		[integer]	   : Minimum length of a CNV (default: 1000 bps)
	--no-sr					   : Split read mapping is disabled
	
	Information:
	--version                  : Print version and exit.
	--help                     : Print this help screen and exit.


Citation
========
Currently under review...
