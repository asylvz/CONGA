

#ifndef FREE_H_
#define FREE_H_

#include "bam_data.h"
#include "svs.h"

void free_DS(bam_info* in_bam, parameters *params);
void free_SVs(svs* sv_all, int sv_count);
void free_splits(bam_info* in_bam);


#endif /* FREE_H_ */
