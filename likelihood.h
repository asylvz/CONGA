
#ifndef __WEIGHTED_SETCOVER_MULTICOLOR__
#define __WEIGHTED_SETCOVER_MULTICOLOR__

#include <stdlib.h>
#include <stdio.h>
#include "bam_data.h"
#include "svs.h"
#include "read_distribution.h"

#define inf 10000000 // a constant to represent infinity
#define TRUE 1
#define FALSE 0
#define true 1
#define false 0

#define WRONGMAP_WINDOW 50
#define WRONGMAP_WINDOW_DEL 1000

extern int total_dels;
extern int total_dups;

void find_SVs( bam_info *in_bam, parameters *params, FILE* fp_del, FILE* fp_dup, FILE* fp_SVs, char* chr_name);
#endif
