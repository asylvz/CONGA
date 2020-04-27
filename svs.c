
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
	FILE *sv_file_del, *sv_file_dup;
	sv_file_del = safe_fopen(params->del_file, "r");
	sv_file_dup = safe_fopen(params->dup_file, "r");

	//Count the lines first
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
	//fprintf(stderr,"%d\n", (*del_count));

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


	(*vars_del) = ( svs*) getMem( sizeof( svs) * ((*del_count) + 1));
	(*vars_dup) = ( svs*) getMem( sizeof( svs) * ((*dup_count) + 1));

	//fprintf(stderr,"%d dels and %d dups available\n",(*del_count), (*dup_count));
	int cnt = 0;
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
		(*vars_del)[cnt].k_mer = 0;
		(*vars_del)[cnt].expected_kmer = 0.0;
		(*vars_del)[cnt].border_rp = 0;
		(*vars_del)[cnt].low_mappability = false;

		cnt++;

		line[0] = '\0';
	}
	fclose(sv_file_del);

	cnt = 0;
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
		(*vars_dup)[cnt].k_mer = 0;
		(*vars_dup)[cnt].border_rp = 0;
		(*vars_dup)[cnt].low_mappability = false;

		cnt++;

		line[0] = '\0';
	}
	fclose(sv_file_dup);

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
		for(i = 0; i < del_count; i++)
		{
			if((start_sv >= vars_del[i].start && start_sv <= vars_del[i].end) || (start_sv >= vars_del[i].start && start_sv <= vars_del[i].end))
			{
				vars_del[i].low_mappability = true;
				cnt_del++;
			}
		}
		for(i = 0; i < dup_count; i++)
		{
			if((start_sv >= vars_dup[i].start && start_sv <= vars_dup[i].end) || (start_sv >= vars_dup[i].start && start_sv <= vars_dup[i].end))
			{
				vars_dup[i].low_mappability = true;
				cnt_dup++;
			}
		}
		line[0] = '\0';
	}
	fprintf(stderr,"Filtered Low Mappability Regions (%d dels and %d dups)\n", cnt_del, cnt_dup);
}
