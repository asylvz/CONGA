#ifndef READ_DISTRIBUTION_H_
#define READ_DISTRIBUTION_H_

#include "bam_data.h"
#include "common.h"

#define WINDOWSIZE 100
#define WINDOWSLIDE 100
#define STDEVLIMIT 5.0

/* Number of chromosomoes used in read depth analysis */
#define RDCHRCOUNT 24
#define CHRX 22
#define CHRY 23

void init_rd_per_chr( bam_info* in_bam, parameters* param, int chr_index);
void calc_mean_per_chr( parameters *params, bam_info* in_bam, int chr_index);
#endif /* READ_DISTRIBUTION_H_ */
