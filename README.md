# CONGA 
(**CO**py **N**umber variation **G**enotyping in **A**ncient genomes and low-coverage sequencing data) 

CONGA is a genotyping algorithm for Copy Number Variations (large deletions and duplications) in ancient genomes. It is tailored for calling homozygous and heterozygous CNV genotypes at low depths of coverage using read-depth and read-pair information from a BAM file with Illumina short single-end reads.

Please feel free to send me an e-mail (asoylev@gmail.com), or better yet open an issue for your questions.


## Requirements

 * htslib	(included as submodule; http://htslib.org/)
 	* libbz2 and liblzma are required by htslib (e.g., "sudo apt-get install libbz2-dev" and "sudo apt-get install liblzma-dev")
 * sonic  	(included as submodule; https://github.com/calkan/sonic)

CONGA is developed and tested using Linux operating system... 

## Downloading and running

	git clone https://github.com/asylvz/CONGA --recursive
	cd CONGA && make libs
	make

	./conga -i myinput.bam --ref human_g1k_v37.fasta --sonic human_g1k_v37.sonic  \
		--dels known_dels.bed --dups known_dups.bed --out myoutput


## SONIC file (required)

You need to input a SONIC file as input to CONGA (--sonic). This file contains some annotation based on the reference genome that you use. You can use one of the already created ones from: https://github.com/BilkentCompGen/sonic-prebuilt

 * human_g1k_v37.sonic: SONIC file for Human Reference Genome GRCh37 (1000 Genomes Project version)
 * ucsc_hg19.sonic: SONIC file for the human reference genome, UCSC version build hg19.
 * ucsc_hg38.sonic: SONIC file for the human reference genome build 38.

If you are working with a different reference genome, you need to create the SONIC file yourself. This is a straightforward process; please refer to the SONIC development repository: https://github.com/calkan/sonic/


## Sample Genotype file (required)

	1	668630		850204
	1	963826		974172
	1	1171539		1179729
	1	1249799		1265722
	1	2374226		2379823
	...

* The columns are "Chromosome Name" (TAB) "Start Position of a CNV" (TAB) "End Postion of a CNV"
* This file should be seperate for duplications and deletions if both are to be genotyped.


## Sample Mappability File (optional)

	1	63913643	63913648	0.2
	1	63913648	63913649	0.25
	1	63913649	63913653	0.5
	1	63913653	63913659	0.333333
	...

Using a mappability file (--mappability) increases the accuracy of CONGA's predictions. We used the 100-mer mappability file from http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeMapability/ and converted the bigWig file into a BED file using "bigWigToBedGraph".

* The columns are "Chromosome Name" (TAB) "Start Position of a CNV" (TAB) "End Postion of a CNV" (TAB) "Mappability value"
    * Note that the mappability value should be between [0,1], where lower values indicate lower mappability intervals, i.e., repeat-rich regions, etc. 


## All parameters

	--input 		[BAM file]         : Input files in sorted and indexed BAM format. (required)
	--out   		[output prefix]    : Prefix for the output file names. (required)
	--ref   		[reference genome] : Reference genome in FASTA format. (required)
	--sonic 		[sonic file]       : SONIC file that contains assembly annotations. (required)
	--dels          	[bed file]         : Known deletion SVs in bed format
	--dups          	[bed file]         : Known duplication SVs in bed format
	--mappability   	[bed file]         : Mappability file in BED format
	--first-chr     	[chromosome index] : The index of the first chromosome for genotyping in your BAM.
	--last-chr      	[chromosome index] : The index of the last chromosome for genotyping in your BAM.
	--min-read-length	[integer]	   : Minimum length of a read to be processed for RP (default: 60 bps)
	--min-sv-size		[integer]	   : Minimum length of a CNV (default: 1000 bps)
	--min-mapq		[integer]	   : Minimum mapping quality threshold for reads (default: -1)
	--c-score               [float]            : Minimum c-score to filter variants (More conservative with lower values, default: 0.5).
	--rp                    [integer]          : Enable split-read and set minimum read-pair support for a duplication (Suggested for >5x only).
	
	Information:
	--version                  		   : Print version and exit.
	--help 		                           : Print this help screen and exit.


## Citation

*Arda Soylev, Sevim Seda Cokoglu, Dilek Koptekin, Can Alkan, and Mehmet Somel. "CONGA: Copy number variation genotyping in ancient genomes and low-coverage sequencing data." bioRxiv (2021).*
