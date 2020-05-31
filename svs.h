
#ifndef SVS_H_
#define SVS_H_

#include "common.h"

#define ROW_DELIMITERS " \t\r\n"
#define ROW_DELIMITERS_MERS " \t\r\n>"

typedef struct _svs
{
	int id;
	char* chr_name;
	int start;
	int end;
	char SV_type;
	//int observed_rd_filtered;
	int observed_rd_all;
	//float expected_rd_filtered;
	float expected_rd_all;
	int rp;
	int border_rp;
	int k_mer;
	double expected_kmer;
	float copy_number;
	//double likelihood_filtered;
	double likelihood_kmer;
	double likelihood_unfiltered;
	double mappability;
	bool low_mappability;
}svs;

void load_mappability_regions(bam_info* in_bam, parameters *params, char* chr);
void check_low_mappability(parameters *params, svs* vars_del, svs* vars_dup, char* chr, int del_count, int dup_count);
int load_known_SVs(svs** vars_del, svs** vars_dup, parameters *params, char* chr, int* del_count, int* dup_count);
#endif /* SVS_H_ */
