#ifndef __COMMANDLINE
#define __COMMANDLINE

#include "common.h"

int parse_cmd_line( int, char**, parameters*);
void parse_bam_list( parameters** params);
void print_help( void);

#endif
