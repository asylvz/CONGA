# CONGA 
(**CO**py **N**umber variation **G**enotyping in **A**ncient genomes and low-coverage sequencing data) 

CONGA is a genotyping algorithm for Copy Number Variations (large deletions and duplications) in ancient genomes. It is tailored for calling homozygous and heterozygous CNV genotypes at low depths of coverage using read-depth and read-pair information from a BAM file with Illumina short single-end reads.

Please feel free to send me an e-mail (asoylev@gmail.com), or better yet open an issue for your questions.


## Installation

### Bioconda (recommended)

	conda install -c bioconda conga

### From source

#### Requirements

 * htslib	(included as submodule; http://htslib.org/)
 	* libbz2, liblzma, libcurl are required by htslib

*Installing development libraries (requires sudo access):* ***"sudo apt-get install zlib1g-dev libbz2-dev liblzma-dev libcurl4-openssl-dev"***

#### Build

	git clone https://github.com/asylvz/CONGA --recursive
	cd CONGA && make libs && make

	./conga --input myinput.bam --ref human_g1k_v37.fasta \
		--dels known_dels.bed --dups known_dups.bed --out myoutput

* Optionally, you can provide a repeats file to filter out reads in satellite regions: `--reps repeats.bed` (only used with `--rp`)

## Docker Usage

	docker pull asylvz/conga

	docker run --user=$UID -v /home/projects/conga:/input -v /home/projects/conga:/output asylvz/conga --input /input/myinput.bam --ref /input/human_g1k_v37.fasta --dels /input/known_dels.bed --dups /input/known_dups.bed --out /output/mydockertest


## Output

CONGA produces three output files based on the `--out` prefix:

* `<prefix>_svs.bed` — Filtered SV calls (deletions and duplications combined), filtered by c-score and mappability.
* `<prefix>_dels.bed` — All deletion genotyping results with detailed metrics (observed/expected read depth).
* `<prefix>_dups.bed` — All duplication genotyping results with detailed metrics (observed/expected read depth).


## Sample Genotype file (required)

### You can use the "svcalls.sh" script under /scripts to generate CNV calls from the 1K Phase 3 SV call set

	1	668630		850204
	1	963826		974172
	1	1171539		1179729
	1	1249799		1265722
	1	2374226		2379823
	...

* The columns are "Chromosome Name" (TAB) "Start Position of a CNV" (TAB) "End Position of a CNV"
* Provide separate files for deletions (`--dels`) and duplications (`--dups`).


## Mappability file (optional)

	1	63913643	63913648	0.2
	1	63913648	63913649	0.25
	1	63913649	63913653	0.5
	1	63913653	63913659	0.333333
	...

Using a mappability file (--mappability) increases the accuracy of CONGA's predictions. We used the 100-mer mappability file from http://hgdownload.cse.ucsc.edu/goldenpath/hg19/encodeDCC/wgEncodeMapability/ and converted the bigWig file into a BED file using "bigWigToBedGraph".

* The columns are "Chromosome Name" (TAB) "Start" (TAB) "End" (TAB) "Mappability value"
    * Note that the mappability value should be between [0,1], where lower values indicate lower mappability intervals, i.e., repeat-rich regions, etc. 


## Repeats file (optional, only used with --rp)

You can optionally provide a repeats file to filter out reads in satellite regions using `--reps`. This is only effective when split-read mode is enabled with `--rp`. The file should be tab-delimited with 5 columns:

	chr1	10000	20000	(CATTC)n	Satellite
	chr1	50000	60000	L1ME1		LINE/L1

* The columns are "Chromosome Name" (TAB) "Start" (TAB) "End" (TAB) "Repeat Type" (TAB) "Repeat Class"
* This can be generated from RepeatMasker output. All repeats are loaded; satellite filtering (looking for "Satel" in type or class) is done at runtime.


## Parameters

**Required:**

	--input            Input file in sorted and indexed BAM format.
	--out              Prefix for the output file names.
	--ref              Reference genome in FASTA format.
	--dels             Known deletion SVs in BED format.
	--dups             Known duplication SVs in BED format.

**Optional:**

	--mappability      Mappability file in BED format.
	--reps             Repeat regions file to filter out satellite regions (only used with --rp).
	--rp               Enable split-read and set read-pair support threshold for duplications (suggested for >5x).
	--first-chr        Index of the first chromosome for genotyping (default: all).
	--last-chr         Index of the last chromosome for genotyping (default: all).
	--min-read-length  Minimum read length for read-pair analysis (default: 60).
	--min-sv-size      Minimum CNV length (default: 1000).
	--min-mapq         Minimum mapping quality threshold (default: no filter).
	--c-score          Minimum c-score for filtering; lower is more conservative (default: 0.5).

**Info:**

	--version          Print version and exit.
	--help             Print help screen and exit.


## Citation

*Arda Söylev, Sevim Seda Çokoglu, Dilek Koptekin, Can Alkan, and Mehmet Somel. "CONGA: Copy number variation genotyping in ancient genomes and low-coverage sequencing data." PLOS Computational Biology 18, no. 12 (2022): e1010788.* https://doi.org/10.1371/journal.pcbi.1010788
