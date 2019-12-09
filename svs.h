
#ifndef SVS_H_
#define SVS_H_

#include "common.h"

#define ROW_DELIMITERS " \t\r\n"

typedef struct _svs
{
	int id;
	char* chr_name;
	int start;
	int end;
	char SV_type;
	long depth;
	int rp;
	float copy_number;
	double del_likelihood;
	double dup_likelihood;
	double cnv_probability[10];
}svs;


int load_known_SVs(svs** vars_del, svs** vars_dup, parameters *params, char* chr, int* del_count, int* dup_count);
#endif /* SVS_H_ */
