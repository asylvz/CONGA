
This script can be used to generate CNV call set using 1000 Genomes Phase 3 SVs to be used to genotype CNVs using CONGA.

- If you want to generate deletions of size >1000 bps, use:
./svcalls 1000 or ./svcalls 1000 del 

- For duplications of size 1000, use:
./svcalls 1000 dup

* 1000 can be replaced with any value.
* You may need to use "chmod 755 svcalls.sh" first
* Note that it removes mobile element deletions (e.g., ALU, L1)

Please cite "Arda Söylev, Sevim Seda Çokoglu, Dilek Koptekin, Can Alkan, and Mehmet Somel. "CONGA: Copy number variation genotyping in ancient genomes and low-coverage sequencing data." PLOS Computational Biology 18, no. 12 (2022): e1010788. https://doi.org/10.1371/journal.pcbi.1010788"


You can open an issue in Github or e-mail me: asoylev@gmail.com for any issues
