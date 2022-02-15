
#include <stdio.h>
#include <stdlib.h>
#include "svs.h"


int load_known_SVs(svs** vars_del, svs** vars_dup, parameters *params, char* chr, int* del_count, int* dup_count)
{
	int i;
	char* return_value;
	char line[512];
	char *sv_start, *sv_end, *chr_name;
	int start_sv, end_sv;
	FILE *sv_file_del = NULL, *sv_file_dup = NULL;

	if(params->del_file != NULL)
		sv_file_del = safe_fopen(params->del_file, "r");

	if(params->dup_file != NULL)
		sv_file_dup = safe_fopen(params->dup_file, "r");

	//Count the lines first

	if(params->del_file != NULL)
	{
		while(!feof(sv_file_del))
		{
			return_value = fgets(line, 512, sv_file_del);

			if (line == NULL)
				continue;

			// If the read line is empty
			int isEmpty = 1;
			int len = strlen (line);
			for (i = 0; i < len; i++)
			{
				if (!isspace (line[i]))
				{
					isEmpty = 0;
					break;
				}
			}
			if (isEmpty)
				continue;

			chr_name = strtok (line, ROW_DELIMITERS);

			sv_start = strtok (NULL, ROW_DELIMITERS);
			start_sv = atoi(sv_start);

			sv_end = strtok (NULL, ROW_DELIMITERS);
			end_sv = atoi(sv_end);

			if(strcmp(chr_name, chr) != 0 || (end_sv - start_sv) < params->min_sv_size)
			{
				line[0] = '\0';
				continue;
			}

			(*del_count)++;

			line[0] = '\0';

		}
		rewind(sv_file_del);
	}

	if(params->dup_file != NULL)
	{
		while(!feof(sv_file_dup))
		{
			return_value = fgets(line, 512, sv_file_dup);

			if (line == NULL)
				continue;

			// If the read line is empty
			int isEmpty = 1;
			int len = strlen (line);
			for (i = 0; i < len; i++)
			{
				if (!isspace (line[i]))
				{
					isEmpty = 0;
					break;
				}
			}
			if (isEmpty)
				continue;

			chr_name = strtok (line, ROW_DELIMITERS);

			sv_start = strtok (NULL, ROW_DELIMITERS);
			start_sv = atoi(sv_start);

			sv_end = strtok (NULL, ROW_DELIMITERS);
			end_sv = atoi(sv_end);

			if(strcmp(chr_name, chr) != 0 || (end_sv - start_sv) < params->min_sv_size)
			{
				line[0] = '\0';
				continue;
			}

			(*dup_count)++;

			line[0] = '\0';
		}
		rewind(sv_file_dup);
	}

	if(params->del_file != NULL)
		(*vars_del) = ( svs*) getMem( sizeof( svs) * ((*del_count) + 1));

	if(params->dup_file != NULL)
		(*vars_dup) = ( svs*) getMem( sizeof( svs) * ((*dup_count) + 1));

	//fprintf(stderr,"%d dels and %d dups available\n",(*del_count), (*dup_count));
	int cnt = 0;

	if(params->del_file != NULL)
	{
		while(!feof(sv_file_del))
		{
			return_value = fgets (line, 512, sv_file_del);

			if (line == NULL)
				continue;
			// If the read line is empty
			int isEmpty = 1;
			int len = strlen (line);
			for (i = 0; i < len; i++)
			{
				if (!isspace (line[i]))
				{
					isEmpty = 0;
					break;
				}
			}
			if (isEmpty)
				continue;

			chr_name = strtok (line, ROW_DELIMITERS);

			sv_start = strtok (NULL, ROW_DELIMITERS);
			start_sv = atoi(sv_start);

			sv_end = strtok (NULL, ROW_DELIMITERS);
			end_sv = atoi(sv_end);

			if(strcmp(chr_name, chr) != 0 || (end_sv - start_sv) < params->min_sv_size)
			{
				line[0] = '\0';
				continue;
			}


			(*vars_del)[cnt].chr_name = NULL;
			(*vars_del)[cnt].chr_name = ( char*) getMem( sizeof( char) * ( strlen( chr_name) + 1));
			strncpy( (*vars_del)[cnt].chr_name, chr_name, ( strlen( chr_name) + 1));

			(*vars_del)[cnt].id = cnt;
			(*vars_del)[cnt].start = start_sv;
			(*vars_del)[cnt].end = end_sv;
			(*vars_del)[cnt].SV_type = DELETION;
			(*vars_del)[cnt].rp = 0;
			(*vars_del)[cnt].observed_rd_sv = 0;
			(*vars_del)[cnt].expected_rd_sv = 0.0;
			(*vars_del)[cnt].border_rp = 0;
			(*vars_del)[cnt].low_mappability = false;

			cnt++;

			line[0] = '\0';
		}
		fclose(sv_file_del);
	}

	cnt = 0;
	if(params->dup_file != NULL)
	{
		while(!feof(sv_file_dup))
		{
			return_value = fgets (line, 512, sv_file_dup);

			if (line == NULL)
				continue;
			// If the read line is empty
			int isEmpty = 1;
			int len = strlen (line);
			for (i = 0; i < len; i++)
			{
				if (!isspace (line[i]))
				{
					isEmpty = 0;
					break;
				}
			}
			if (isEmpty)
				continue;

			chr_name = strtok (line, ROW_DELIMITERS);

			sv_start = strtok (NULL, ROW_DELIMITERS);
			start_sv = atoi(sv_start);

			sv_end = strtok (NULL, ROW_DELIMITERS);
			end_sv = atoi(sv_end);

			if(strcmp(chr_name, chr) != 0 || (end_sv - start_sv) < params->min_sv_size)
			{
				line[0] = '\0';
				continue;
			}

			//fprintf(stderr, "%s - %s - %d - %d - %d\n", chr, chr_name, start_sv, end_sv, end_sv-start_sv);

			(*vars_dup)[cnt].chr_name = NULL;
			(*vars_dup)[cnt].chr_name = ( char*) getMem( sizeof( char) * ( strlen( chr_name) + 1));
			strncpy( (*vars_dup)[cnt].chr_name, chr_name, ( strlen( chr_name) + 1));

			(*vars_dup)[cnt].id = cnt;
			(*vars_dup)[cnt].start = start_sv;
			(*vars_dup)[cnt].end = end_sv;
			(*vars_dup)[cnt].SV_type = DUPLICATION;
			(*vars_dup)[cnt].rp = 0;
			(*vars_dup)[cnt].observed_rd_sv = 0;
			(*vars_dup)[cnt].expected_rd_sv = 0.0;
			(*vars_dup)[cnt].border_rp = 0;
			(*vars_dup)[cnt].low_mappability = false;

			cnt++;

			line[0] = '\0';
		}
		fclose(sv_file_dup);
	}
	return 1;
}


void check_low_mappability(parameters *params, svs* vars_del, svs* vars_dup, char* chr, int del_count, int dup_count)
{
	int i;
	char* return_value;
	char line[512];
	char *sv_start, *sv_end, *chr_name;
	int start_sv, end_sv;
	FILE *low_map_file;
	low_map_file = safe_fopen(params->low_map_regions, "r");
	int cnt_del = 0; int cnt_dup = 0;

	//Count the lines first
	while(!feof(low_map_file))
	{
		return_value = fgets(line, 512, low_map_file);

		if (line == NULL)
			continue;

		// If the read line is empty
		int isEmpty = 1;
		int len = strlen (line);
		for (i = 0; i < len; i++)
		{
			if (!isspace (line[i]))
			{
				isEmpty = 0;
				break;
			}
		}
		if (isEmpty)
			continue;

		chr_name = strtok (line, ROW_DELIMITERS);

		sv_start = strtok (NULL, ROW_DELIMITERS);
		start_sv = atoi(sv_start);

		sv_end = strtok (NULL, ROW_DELIMITERS);
		end_sv = atoi(sv_end);

		if(strcmp(chr_name, chr) != 0)
		{
			line[0] = '\0';
			continue;
		}
		if(params->del_file)
		{
			for(i = 0; i < del_count; i++)
			{
				if((start_sv >= vars_del[i].start && start_sv <= vars_del[i].end) || (end_sv >= vars_del[i].start && end_sv <= vars_del[i].end))
				{
					vars_del[i].low_mappability = true;
					cnt_del++;
				}
			}
		}
		if(!params->dup_file)
		{
			for(i = 0; i < dup_count; i++)
			{
				if((start_sv >= vars_dup[i].start && start_sv <= vars_dup[i].end) || (end_sv >= vars_dup[i].start && end_sv <= vars_dup[i].end))
				{
					vars_dup[i].low_mappability = true;
					cnt_dup++;
				}
			}
		}
		line[0] = '\0';
	}
	fclose(low_map_file);
	fprintf(stderr,"Filtered Low Mappability Regions (%d dels and %d dups)\n", cnt_del, cnt_dup);
}

void load_mappability_regions(bam_info* in_bam, parameters *params, char* chr)
{
	int i;
	char* return_value;
	char line[512];
	char *sv_start, *sv_end, *chr_name, *mappability_tmp;
	int start_sv, end_sv;
	float mappability;
	FILE *mappability_file;
	mappability_file = safe_fopen(params->mappability_file, "r");

	//Count the lines first
	while(!feof(mappability_file))
	{
		return_value = fgets(line, 512, mappability_file);

		if (line == NULL)
			continue;

		// If the read line is empty
		int isEmpty = 1;
		int len = strlen (line);
		for (i = 0; i < len; i++)
		{
			if (!isspace (line[i]))
			{
				isEmpty = 0;
				break;
			}
		}
		if (isEmpty)
			continue;

		chr_name = strtok (line, ROW_DELIMITERS);

		if(strcmp(chr_name, chr) != 0)
		{
			line[0] = '\0';
			continue;
		}

		sv_start = strtok (NULL, ROW_DELIMITERS);
		start_sv = atoi(sv_start);

		sv_end = strtok (NULL, ROW_DELIMITERS);
		end_sv = atoi(sv_end);

		mappability_tmp = strtok (NULL, ROW_DELIMITERS);
		mappability = atof(mappability_tmp);

		int i;
		for(i = start_sv; i <= end_sv; i++)
		{
			in_bam->mappability[i] = mappability;
		}

		line[0] = '\0';
	}
	fclose(mappability_file);
	//fprintf(stderr,"Added mappability information for each position\n");
}

