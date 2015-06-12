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

/******************************************************************/
/* function: mg_mesh_2_matlab */
/* converts mesh to matlab format */
int mg_mesh_2_matlab(mg_Mesh *Mesh, mg_Front *Front, char *FileName);

/******************************************************************/
/* function: mg_write_mesh */
/* writes mesh to a file with connectivities and boundary information */
int mg_write_mesh(mg_Mesh *Mesh, char *FileName);

/******************************************************************/
/* function: mg_read_mesh */
/* reads mesh from file with connectivities and boundary information */
int mg_read_mesh(mg_Mesh **pMesh, char *FileName);

/******************************************************************/
/* function: mg_read_mesh */
/* reads mesh from file with connectivities and boundary information */
int mg_read_geo(mg_Geometry **pGeo, char *FileName);

#endif
