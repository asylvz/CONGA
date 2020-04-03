CONGA
======


Requirements
============

 * zlib   (http://www.zlib.net)
 * htslib (included as submodule; http://htslib.org/)
 * sonic  (included as submodule; https://github.com/calkan/sonic)

Fetching CONGA
===============

	git clone https://github.com/asylvz/CONGA --recursive

Compilation
===========

Type:

	make libs
	make

Running SvDepth
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
 * mm9.sonic: SONIC file for the mouse reference genome version mm9.
	* Also download the reference genome at: http://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/chromFa.tar.gz. Deflate the tar archive and concatenate all chromosomes into a single FASTA file.
 * mm10.sonic: SONIC file for the mouse reference genome version mm10.
	* Also download the reference genome at: http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz. Deflate the tar archive and concatenate all chromosomes into a single FASTA file.

Building the SONIC file
=======================

If you are working with a different reference genome please refer to the SONIC development repository: https://github.com/calkan/sonic/

The README.md file includes documentation on how to obtain the necessary files for different genomes from the UCSC Genome Browser.


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
