//
//  2dmg_utils.h
//  2dmg
//
//  Created by Marco Ceze on 11/11/14.
//  Copyright (c) 2014 Marco Ceze. All rights reserved.
//

#ifndef ___dmg___dmg_utils__
#define ___dmg___dmg_utils__

#include "2dmg_def.h"
#include "2dmg_struct.h"

/******************************************************************/
/* function:  mg_alloc*/
/* wrapper for malloc with error handling*/
int mg_alloc( void **pchunk, int n, int size);

/******************************************************************/
/* function:  mg_free*/
/* initializes a mg_List*/
void mg_free( void *chunk);

/******************************************************************/
/* function:  mg_alloc2*/
/* wrapper for 2d malloc with error handling*/
int mg_alloc2( void ***pchunk, int n1, int n2, int size);

/******************************************************************/
/* function:  mg_free2*/
/* wrapper for 2d free*/
void mg_free2(void **chunk);

/******************************************************************/
/* function:  mg_realloc*/
/* wrapper for 2d realloc*/
int mg_realloc( void **pchunk, const int n, const int size);

/******************************************************************/
/* function:  mg_init_list*/
/* initializes a mg_List*/
void mg_init_list(mg_List *list);

/******************************************************************/
/* function: mg_create_mesh */
/* creates and initilizes a mesh structure */
int mg_create_mesh(mg_Mesh **pMesh);

/******************************************************************/
/* function: mg_destroy_mesh */
/* deallocates the memory occupied by the mesh structure */
void mg_destroy_mesh(mg_Mesh *Mesh);

/******************************************************************/
/* function: mg_init_face */
/* initializes a face structure */
void mg_init_face(mg_FaceData *face);

/******************************************************************/
/* function: mg_free_front_face */
/* frees a front face structure */
void mg_free_front_face(struct mg_FrontFace *FFace);

/******************************************************************/
/* function: mg_add_2_ord_set */
/* adds "entry" to a ordered "set" and keeps it ordered */
int mg_add_2_ord_set(const int entry, int *set_size, int **set,
                     int **orig_rank, bool AllowRepeat);

/******************************************************************/
/* function: mg_rm_frm_ord_set */
/* removes entry from set */
int mg_rm_frm_ord_set(const int entry, int *set_size, int **set,
                      int n2Rm, int *nRmd);

/******************************************************************/
/* function: mg_build_connectivity */
/* builds Node2Elem and Node2Face connectivities */
int mg_build_connectivity(mg_Mesh *Mesh);

/******************************************************************/
/* function: mg_init_ord_data_list */
/* initializes an ordered data list*/
void mg_init_ord_data_list(mg_OrderedDataList *List, int DataSize);

/******************************************************************/
/* function: mg_free_ord_data_list */
/* frees the memory ocuppied by an ordered data list*/
void mg_free_ord_data_list(mg_OrderedDataList *List);

/******************************************************************/
/* function: mg_rm_frm_ord_data_list */
/* removes "data" corresponding to "entry" from "List" */
int mg_rm_frm_ord_data_list(int entry, mg_OrderedDataList *List);

/******************************************************************/
/* function: mg_add_2_ord_data_list */
/* adds "data" corresponding to "entry" to "list". It outputs the
 data rank position in the list. It will substitute data if
 AllowRepeat == false */
int mg_add_2_ord_data_list(int entry, void **Data,
                           mg_OrderedDataList *List, int *prank,
                           bool AllowRepeat);

/******************************************************************/
/* function: mg_edges_intersect */
/* checks if edges defined by *X0 and *X1 intersect, if so,
 xint receives the intersection point*/
bool mg_edges_intersect(double X0[4], double X1[4], double *xint);

/******************************************************************/
/* function: mg_coord_inside_elem */
/* checks if a coordinate is inside an element */
bool mg_coord_inside_elem(mg_Mesh *Mesh, int elem, double coord[2]);

/******************************************************************/
/* function: mg_find_elem_frm_coord */
/* finds element containing a point given by its coordinates  */
int mg_find_elem_frm_coord(mg_Mesh *Mesh, int elem_start,
                           double coord[2], int *pelem);

/******************************************************************/
/* function: mg_check_exist */
/* checks if "num" exists in "vector" and returns "index" */
void mg_check_exist(int num, int size, int *vector, int *index);

/******************************************************************/
/* function: mg_mesh_2_matlab */
/* converts mesh to matlab format */
int mg_mesh_2_matlab(mg_Mesh *Mesh, mg_Front *Front, char *FileName);

/******************************************************************/
/* function: mg_write_mesh */
/* writes mesh to a file with connectivities and boundary information */
int mg_write_mesh(mg_Mesh *Mesh, char *FileName);


#endif /* defined(___dmg___dmg_utils__) */
