//
//  2dmg_math.h
//  2dmg
//
//  Created by Marco Ceze on 11/11/14.
//  https://github.com/mceze/2dmg
//

#ifndef ___dmg___dmg_math__
#define ___dmg___dmg_math__

#include <stdio.h>
#include "2dmg_struct.h"
#include "2dmg_def.h"
#include "2dmg_metric_struct.h"


/******************************************************************/
/* function:  mg_limited_pair */
/* Pairs (x,y)=>z uniquely using ylim key.
 The maximum value for y is ylim. x >= 0; 0 =< y <= ylim.  */
int mg_limited_pair(const int x, const int y, int *z, const int ylim);

/******************************************************************/
/* function:  mg_limited_pair_inv */
/* Decodes z=>(x,y) uniquely using ylim key.
 The maximum value for y is ylim. x >= 0; 0 =< y <= ylim.  */
int mg_limited_pair_inv(int *x, int *y, const int z, const int ylim);

/******************************************************************/
/* function:  mg_binary_search */
/* searches for "target" in set between positions "begin" and "end".
 set[rank] == target if "target" is returned err_OK, "rank" contains the
 rightful position to insert "target" when returned err_NOT_FOUND.
 */
int mg_binary_search(const int target, const int *set, int begin,
                     int end, int *rank);

/******************************************************************/
/* function: mg_circumcircle */
/* builds a circulcircle from 3 coordinates listed in "coord" */
int mg_circumcircle(double coord[3][2], double center[2], double *radius);

/******************************************************************/
/* function: mg_get_metric */
/* gets metric value at a (x,y) */
int mg_get_metric(mg_Metric *Metric, double *x, double *y, int np, double *M);

/******************************************************************/
/* function: mg_metric_dist */
/* computes metric distance between 2 points */
int mg_metric_dist(mg_Metric *Metric, int order, double *coord,
                   double *dist);

/******************************************************************/
/* function: mg_metric_length */
/* computes metric length of a segment */
int mg_metric_length(mg_Metric *Metric, mg_Segment *Segment, int order,
                     double *length);

/******************************************************************/
/* function: mg_mxm */
/* product of 2 full matrices*/
void
mg_mxm(int m, int n, int l, double *a, double *b, double *c);

/******************************************************************/
/* function: mg_eig2 */
/* eigenvectors and eigenvalues for 2by2 matrices*/
int
mg_eig2(const double M[4], double V[4], double lambda[2]);

/******************************************************************/
/* function: mg_circumellipse */
int
mg_circumellipse(double *coord, mg_Ellipse *Ellipse);

/******************************************************************/
/* function: mg_inside_ellipse */
bool
mg_inside_ellipse(double *coord, mg_Ellipse *Ellipse);

/******************************************************************/
/* function: mg_ellipse_frm_face_p */
/* creates a steiner ellipse using a front face and a point */
int mg_ellipse_frm_face_p(mg_Mesh *Mesh, mg_FaceData *face,
                           double *Popt, mg_Ellipse *Ellipse);

#endif /* defined(___dmg___dmg_math__) */
