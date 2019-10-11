
#include <stdio.h>
#include "free.h"

void free_SVs(svs* sv_all, int sv_count)
{
	int count;
	for( count = 0; count < sv_count; count++)
	{
		//fprintf(stderr,"\n%d\n", count);
		if(sv_all[count].chr_name != NULL)
			free(sv_all[count].chr_name);
	}
	free(sv_all);
}

void free_DS(bam_info* in_bam, parameters *params)
{
	/* Free the read depth array*/
	free( in_bam->read_depth_per_chr);

	/* Free bams and related libraries */
	free( in_bam->sample_name);
	free( in_bam);

	/* Free params struct */
	free( params->bam_file);
	free( params->del_file);
	free( params->dup_file);
	free( params->ref_genome);
	free( params->outprefix);
	free( params->outdir);
	free( params->sonic_file);
	free( params);
}
