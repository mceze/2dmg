//
//  2dmg_io.h
//  2dmg
//
//  Created by Marco Ceze on 11/11/14.
//  Copyright (c) 2014 Marco Ceze. All rights reserved.
//

#ifndef _dmg__dmg_io_h
#define _dmg__dmg_io_h

#include "2dmg_struct.h"
#include "2dmg_def.h"


/******************************************************************/
/* function:  mg_scan_n_num */
/* scan "n" doubles from a string  */
int mg_scan_n_num( const char *line, int *n, int *vi, double *vr);

/******************************************************************/
/* function:  mg_read_bgri_file */
/* reads a boundary discretization file  */
int mg_read_bgri_file(const char *InFile, mg_Mesh *Mesh);

/******************************************************************/
/* function:  mg_parse_input_line */
int mg_parse_input_line(char line[], char **pkey,
                        char **pvalue);

/******************************************************************/
/* function:  mg_read_input_file */
int mg_read_input_file(char const ParamFile[]);

/******************************************************************/
/* function:  mg_read_input_file */
int mg_get_input_char(char const ParamName[], char **pvalue);

#endif
