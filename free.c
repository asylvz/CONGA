
#include <stdio.h>
#include "free.h"
#include "split_read.h"


void free_splits(bam_info* in_bam)
{
	int i;
	splitRead *sfcPtr, *sfcPtrNext;

	SplitsInfo *aPtr;
	SplitRow *tmp, *tmp_next;
	posMapSplitRead *ptrSplitMapPtr, *ptrSplitMapNext;

	sfcPtr = in_bam->listSplitRead;
	while( sfcPtr != NULL)
	{
		sfcPtrNext = sfcPtr->next;
		if( sfcPtr->readName != NULL)
			free( sfcPtr->readName);
		if( sfcPtr->chromosome_name != NULL)
			free( sfcPtr->chromosome_name);


		ptrSplitMapPtr = sfcPtr->ptrSplitMap;
		while(ptrSplitMapPtr != NULL)
		{
			ptrSplitMapNext = ptrSplitMapPtr->next;

			if(ptrSplitMapPtr != NULL)
				free( ptrSplitMapPtr);

			ptrSplitMapPtr = ptrSplitMapNext;
		}

		free( sfcPtr);
		sfcPtr = sfcPtrNext;
	}
	in_bam->listSplitRead = NULL;


	aPtr = all_split_reads;
	if( aPtr != NULL)
	{
		tmp = aPtr->head;
		while( tmp != NULL)
		{
			tmp_next = tmp->next;
			free( tmp);
			tmp = tmp_next;
		}

		aPtr->head = NULL;
		aPtr->tail = NULL;

		free( aPtr);
	}
	all_split_reads = NULL;
}

void free_SVs(svs* sv_all, int sv_count)
{
	int count;
	for(count = 0; count < sv_count; count++)
	{
		if(sv_all[count].chr_name != NULL)
			free(sv_all[count].chr_name);
	}
	free(sv_all);
}

void free_DS(bam_info* in_bam, parameters *params)
{
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
	free( params->sonic_info);
	free( params->mappability_file);
	free( params);

	if(hash_table_count != NULL)
		free(hash_table_count);

	if(hash_table_iter != NULL)
		free(hash_table_iter);
}
