SvDepth
======

Requirements
============

 * zlib   (http://www.zlib.net)
 * htslib (included as submodule; http://htslib.org/)
 * sonic  (included as submodule; https://github.com/calkan/sonic)

Fetching SvDepth
===============

	

Compilation
===========

Type:

	make libs
	make

Running SvDepth
===========================

	svdepth -i myinput.bam --ref human_g1k_v37.fasta --sonic human_g1k_v37.sonic  \
		--dels known_dels.bed --dups known_dups.bed --out myoutput


SONIC file (annotations container)
==================================

SONIC files are available under aux/

 * human_g1k_v37.sonic: SONIC file for Human Reference Genome GRCh37 (1000 Genomes Project version)
 	* Also download the reference genome at: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz. 
 * ucsc_hg19.sonic: SONIC file for the human reference genome, UCSC version build hg19.
	* Also download the reference genome at: http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz. Deflate the tar archive and concatenate all chromosomes into a single FASTA file.
 * ucsc_hg38.sonic: SONIC file for the human reference genome build 38.
	* Also download the reference genome at: http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chromFa.tar.gz. Deflate the tar archive and concatenate all chromosomes into a single FASTA file.

Make sure that the same reference was used to align the reads beforehand (BAM file) and to create the SONIC file. The SONIC files and the reference FASTA files linked above are compatible.


All parameters
==============

	--input [BAM file ]        : Input files in sorted and indexed BAM format.
	--out   [output prefix]    : Prefix for the output file names.
	--ref   [reference genome] : Reference genome in FASTA format.
	--sonic [sonic file]       : SONIC file that contains assembly annotations.
	--dels  [bed file]         : Known deletion SVs in bed format\n");
	--dups  [bed file]         : Known duplication SVs in bed format\n");

	
	Information:
	--version                  : Print version and exit.
	--help                     : Print this help screen and exit.
