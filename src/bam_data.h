#ifndef __PROCESSBAM
#define __PROCESSBAM

/* htslib headers */
#include <htslib/sam.h>
#include <htslib/hts.h>
#include <stdbool.h>
#include "common.h"

extern long split_read_count;

/* Maximum sequence/quality length */
#define MAX_SEQ 1000
#define SOFTCLIP_WRONGMAP_WINDOW 50


/* Function Prototypes */
void read_bam( bam_info* in_bam, parameters *params);

#endif
