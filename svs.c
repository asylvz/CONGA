
#include <stdio.h>
#include <stdlib.h>
#include "svs.h"
#include "common.h"


int load_known_SVs(svs** vars, char* input_del, char* input_dup, char* chr, int* del_count, int* dup_count)
{
	int i;
	char* return_value;
	char line[255]; //shouldn't this be longer?
	char *sv_start, *sv_end, *chr_name;
	int start_sv, end_sv;
	FILE *sv_file_del, *sv_file_dup;
	sv_file_del = safe_fopen(input_del, "r");
	sv_file_dup = safe_fopen(input_dup, "r");

	//Count the lines first
	while(!feof(sv_file_del))
	{
		return_value = fgets(line, 255, sv_file_del);

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

		if(strcmp(chr_name, chr) == 0)
			(*del_count)++;

		line[0] = '\0';
	}
	rewind(sv_file_del);

	while(!feof(sv_file_dup))
	{
		return_value = fgets(line, 255, sv_file_dup);

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

		if(strcmp(chr_name, chr) == 0)
			(*dup_count)++;

		line[0] = '\0';
	}
	rewind(sv_file_dup);

	(*vars) = ( svs*) getMem( sizeof( svs) * ((*del_count) + (*dup_count) + 1));

	int cnt = 0;
	while(!feof(sv_file_del))
	{
		return_value = fgets (line, 255, sv_file_del);

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


		(*vars)[cnt].chr_name = NULL;
		(*vars)[cnt].chr_name = ( char*) getMem( sizeof( char) * ( strlen( chr_name) + 1));
		strncpy( (*vars)[cnt].chr_name, chr_name, ( strlen( chr_name) + 1));

		(*vars)[cnt].id = cnt;
		(*vars)[cnt].start = start_sv;
		(*vars)[cnt].end = end_sv;
		(*vars)[cnt].SV_type = DELETION;

		cnt++;

		line[0] = '\0';
	}
	fclose(sv_file_del);

	while(!feof(sv_file_dup))
	{
		return_value = fgets (line, 255, sv_file_dup);

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


		(*vars)[cnt].chr_name = NULL;
		(*vars)[cnt].chr_name = ( char*) getMem( sizeof( char) * ( strlen( chr_name) + 1));
		strncpy( (*vars)[cnt].chr_name, chr_name, ( strlen( chr_name) + 1));

		(*vars)[cnt].id = cnt;
		(*vars)[cnt].start = start_sv;
		(*vars)[cnt].end = end_sv;
		(*vars)[cnt].SV_type = DUPLICATION;

		cnt++;

		line[0] = '\0';
	}
	fclose(sv_file_dup);

	return 1;
}
