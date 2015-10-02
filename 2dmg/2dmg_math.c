//
//  2dmg_math.c
//  2dmg
//
//  Created by Marco Ceze on 11/11/14.
//  https://github.com/mceze/2dmg
//

#include "2dmg_math.h"
#include "2dmg_def.h"
#include "2dmg_metric_analytic.h"
#include "2dmg_struct.h"
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp.h>

/******************************************************************/
/* function:  mg_limited_pair */
/* Pairs (x,y)=>z uniquely using ylim key. 
 The maximum value for y is ylim. x >= 0; 0 =< y <= ylim.  */
int mg_limited_pair(const int x, const int y, int *z, const int ylim)
{
  int ys;
  ys = ylim+1;
  if (x < 0 || y < 0 || y > ylim) return error(err_INPUT_ERROR);
  (*z) = floor(x/ys)*ys*ys+y*ys+x%ys;
  return err_OK;
}

/******************************************************************/
/* function:  mg_limited_pair_inv */
/* Decodes z=>(x,y) uniquely using ylim key. 
 The maximum value for y is ylim. x >= 0; 0 =< y <= ylim.  */
int mg_limited_pair_inv(int *x, int *y, const int z, const int ylim)
{
  int ys;
  ys = ylim+1;
  (*x) = floor(z/(ys*ys))*ys+(z%(ys*ys))%ys;
  (*y) = floor((z%(ys*ys))/ys);
  return err_OK;
}

/******************************************************************/
/* function:  mg_binary_search */
/* searches for "target" in set between positions "begin" and "end".
 set[rank] == target if "target" is returned err_OK, "rank" contains the
 rightful position to insert "target" when returned err_NOT_FOUND.
 */
int mg_binary_search(const int target, const int *set, int begin,
                     int end, int *rank)
{
  int center;
  int size;
  
  if (end < begin)
    return err_INPUT_ERROR;
  
  size = end-begin+1;
  if (size == 1){
    if (target == set[begin]){
      if (rank != NULL)
        (*rank) = begin;
      return err_OK;
    }
    else{
      if (rank != NULL){
        if (target < set[begin])
          (*rank) = begin-1;
        else
          (*rank) = end+1;
      }
      return err_NOT_FOUND;
    }
  }
  
  if (target == set[begin]){
    if (rank != NULL)
      (*rank) = begin;
    return err_OK;
  }
  
  if (target == set[end]){
    if (rank != NULL)
      (*rank) = end;
    return err_OK;
  }
  
  //if within range
  if (target > set[end]){
    if (rank != NULL)
      (*rank) = end+1;
    return err_NOT_FOUND;
  }
  else if (target < set[begin]){
    if (rank != NULL)
      (*rank) = begin-1;
    return err_NOT_FOUND;
  }
  else{
    while (size > 1){
      center = begin + size/2-1;
      if (target > set[center]){
        //in top portion of set
        begin = center+1;
      }
      else {
        end = center;
      }
      size = end-begin+1;
    }
    if (target == set[begin]){
      if (rank != NULL)
        (*rank) = begin;
      return err_OK;
    }
    else{
      if (target > set[begin]){
        if (rank != NULL)
          (*rank) = begin+1;
      }
      else{
        if (rank != NULL)
          (*rank) = begin;
      }
      return err_NOT_FOUND;
    }
  }
  
  return err_OK;
}

/******************************************************************/
/* function: mg_circumcircle */
/* builds a circulcircle from 3 coordinates listed in "coord" */
int mg_circumcircle(double coord[3][2], double center[2], double *radius)
{
  double det, A[2], B[2], C[2], A2[2], B2[2], C2[2], D;
  
  A[0] = coord[0][0];
  B[0] = coord[1][0];
  C[0] = coord[2][0];
  A[1] = coord[0][1];
  B[1] = coord[1][1];
  C[1] = coord[2][1];
  
  //check if nodes are colinear
  det = A[0]*B[1]+B[0]*C[1]+A[1]*C[0];
  det-= (C[0]*B[1]+A[1]*B[0]+A[0]*C[1]);
  if (fabs(det) < MEPS) return error(err_INPUT_ERROR);
  
  A2[0] = A[0]*A[0];
  A2[1] = A[1]*A[1];
  B2[0] = B[0]*B[0];
  B2[1] = B[1]*B[1];
  C2[0] = C[0]*C[0];
  C2[1] = C[1]*C[1];
  
  D = 2.*(A[0]*(B[1]-C[1])+B[0]*(C[1]-A[1])+C[0]*(A[1]-B[1]));
  center[0] = ((A2[0]+A2[1])*(B[1]-C[1])+(B2[0]+B2[1])*(C[1]-A[1])+(C2[0]+C2[1])*(A[1]-B[1]))/D;
  center[1] = ((A2[0]+A2[1])*(C[0]-B[0])+(B2[0]+B2[1])*(A[0]-C[0])+(C2[0]+C2[1])*(B[0]-A[0]))/D;
  //radius computed from definition
  (*radius) = sqrt((A[0]-center[0])*(A[0]-center[0])+(A[1]-center[1])*(A[1]-center[1]));
  
  
  return err_OK;
}

/******************************************************************/
/* function: mg_get_metric */
/* gets metric value at a (x,y) */
int mg_get_metric(mg_Metric *Metric, double *x, double *y, int np, double *M)
{
  switch (Metric->type) {
    case mge_Metric_Uniform:
      mg_metric_uniform(x, y, np, M);
      break;
    case mge_Metric_Analitic1:
      mg_metric_linx(x, y, np, M);
      break;
    case mge_Metric_Analitic2:
      mg_metric_expx(x, y, np, M);
      break;
    case mge_Metric_Analitic3:
      mg_metric_sqx(x, y, np, M);
      break;
    default:
      error(err_NOT_SUPPORTED);
      break;
  }
  return err_OK;
}

/******************************************************************/
/* function: mg_metric_dist */
/* computes metric distance between 2 points */
int mg_metric_dist(mg_Metric *Metric, int order, double *coord,
                   double *dist)
{
  int ierr;
  gsl_integration_glfixed_table *gltable;
  double dx = coord[1]-coord[0], dy = coord[3]-coord[2];
  double ab[2], M[3], x, y, ds2, xgl, wgl;
  ab[0] = dx;
  ab[1] = dy;
  int ip;
  
  if ((gltable = gsl_integration_glfixed_table_alloc(order)) == NULL)
    return error(err_GSL_ERROR);
  
  (*dist) = 0.0;
  for (ip = 0; ip < gltable->n; ip++) {
    gsl_integration_glfixed_point(0.0, 1.0, ip, &xgl, &wgl, gltable);
    x = coord[0]+xgl*dx;
    y = coord[2]+xgl*dy;
    call(mg_get_metric(Metric, &x, &y, 1, M));
    ds2   = metriclen(ab,M);
    (*dist) += wgl*sqrt(ds2);
  }
  
  gsl_integration_glfixed_table_free(gltable);
  
  return err_OK;
}

/******************************************************************/
/* function: mg_metric_length */
/* computes metric length of a segment */
int mg_metric_length(mg_Metric *Metric, mg_Segment *Segment, int order,
                     double *length)
{
  int ierr, iq;
  gsl_integration_glfixed_table *gltable;
  double tq, wq, xq, yq, dxq, dyq, Mq[3], ab[2], dl2;
  
  (*length) = 0.0;
  
  //integration table
  if ((gltable = gsl_integration_glfixed_table_alloc(order)) == NULL)
    return error(err_GSL_ERROR);
  for (iq = 0; iq < gltable->n; iq++) {
    //get parametric coordinate for quadrature point
    gsl_integration_glfixed_point(0.0, 1.0, iq, &tq, &wq, gltable);
    //evaluate global coordinate and tangent
    xq = gsl_interp_eval(Segment->interp[0],Segment->s,
                         Segment->Coord+0*Segment->nPoint, tq,
                         Segment->accel[0]);
    dxq = gsl_interp_eval_deriv(Segment->interp[0],Segment->s,
                                Segment->Coord+0*Segment->nPoint, tq,
                                Segment->accel[0]);
    yq = gsl_interp_eval(Segment->interp[1],Segment->s,
                         Segment->Coord+1*Segment->nPoint, tq,
                         Segment->accel[1]);
    dyq = gsl_interp_eval_deriv(Segment->interp[1],Segment->s,
                                Segment->Coord+1*Segment->nPoint, tq,
                                Segment->accel[1]);
    //matrix value at quadrature point
    call(mg_get_metric(Metric, &xq, &yq, 1, Mq));
    ab[0] = dxq;
    ab[1] = dyq;
    //dl2 = ab^T*M*ab;
    dl2 = metriclen(ab, Mq);
    (*length) += wq*sqrt(dl2);
  }
  
  gsl_integration_glfixed_table_free(gltable);
  
  return err_OK;
}

/******************************************************************/
/* function: mg_lu */
/* decomposes full matrix as A = L*U*/
static int
mg_lu(int n, double *a, double *l, double *u)
{
  
  int i, j, m;
  
  for(i = 0; i < n; i++) for(j = 0; j < n; j++) l[i*n+j] = u[i*n+j] = 0.0;
  
  for(i = 0; i < n; i++) l[i*n+i] = 1.0;
  
  for(m = 0; m < n; m++){
    for(i = m; i < n; i++){
      u[m*n+i] = a[m*n+i];
      for(j = 0; j < m; j++){
        u[m*n+i] -= l[m*n+j]*u[j*n+i];
      }
    }
    for(i = m+1; i < n; i++){
      l[i*n+m] = a[i*n+m];
      for (j = 0; j < m; j++){
        l[i*n+m] -= l[i*n+j]*u[j*n+m];
      }
      if (fabs(u[m*n+m]) < MEPS)
        return error(err_SINGULAR);
      l[i*n+m] = l[i*n+m]/u[m*n+m];
    }
  }
  
  return err_OK;
}

/******************************************************************/
/* function: mg_luinv */
/* inverts L and U factors*/
static void
mg_luinv(int n, double *l, double *u, double *linv, double *uinv)
{
  int i, j, k;
  double sum;
  //Computing inverse of U
  for(j = n-1; j >= 0; j--){
    for(i = n-1; i >= 0; i--){
      if(i == j) uinv[i*n+j] = 1.0/u[i*n+j];
      else {
        sum = 0.0;
        for(k = n-1; k > i; k--) sum += u[i*n+k]*uinv[k*n+j];
        uinv[i*n+j] = -sum/u[i*n+i];
      }
    }
  }
  
  //Computing inverse of L
  for(j = 0; j < n; j++){
    for(i = 0; i < n; i++){
      if(i == j) linv[i*n+j] = 1.0/l[i*n+j];
      else {
        sum = 0.0;
        for(k = 0; k < i; k++) sum += l[i*n+k]*linv[k*n+j];
        linv[i*n+j] = -sum/l[i*n+i];
      }
    }
  }
}
/******************************************************************/
/* function: mg_mxm */
/* product of 2 full matrices*/
void
mg_mxm(int m, int n, int l, double *a, double *b, double *c)
{
  int i, j, k;
  
  for (i = 0; i < m; i++) {
    for (k = 0; k < l; k++) {
      c[i*l+k] = 0.0;
      for (j = 0; j < n; j++) {
        c[i*l+k] += a[i*n+j]*b[j*l+k];
      }
    }
  }
}

/******************************************************************/
/* function: mg_circumellipse */
int
mg_circumellipse(double *coord, mg_Ellipse *Ellipse)
{
  //builds Steiner circumellipse
  double x[3], y[3], M[4], L[4], U[4], Linv[4], Uinv[4];
  int i;
  double T, D, L1, L2, t, V1norm, V2norm;
  double Cp[4], Mt[4];
  double sxx, sxy, syy, t1, t2, Kt[4];
  
  for (i = 0; i < 3; i++) {
    x[i] = coord[i*2];
    y[i] = coord[i*2+1];
  }
  //centroid
  Ellipse->Ot[0] = (x[0]+x[1]+x[2])/3.0;
  Ellipse->Ot[1] = (y[0]+y[1]+y[2])/3.0;
  
  //transformation matrix
  M[0] = x[1]-x[0];
  M[1] = x[2]-x[0];
  M[2] = y[1]-y[0];
  M[3] = y[2]-y[0];
  
  //inertia coeeficients
  sxy = (x[1]-x[0])*(y[1]-y[0]) + (x[2]-x[0])*(y[2]-y[0]) + (x[2]-x[1])*(y[2]-y[1]);
  sxx = (x[1]-x[0])*(x[1]-x[0]) + (x[2]-x[0])*(x[2]-x[0]) + (x[2]-x[1])*(x[2]-x[1]);
  syy = (y[1]-y[0])*(y[1]-y[0]) + (y[2]-y[0])*(y[2]-y[0]) + (y[2]-y[1])*(y[2]-y[1]);
  
  //principal directions
  T = sxx+syy;
  D = sxx*syy-sxy*sxy;
  //eigenvalues
  t = sqrt(0.25*T*T-D);
  L1 = 0.5*T + t;
  L2 = 0.5*T - t;
  
  if (fabs(sxy) < MEPS){
    Ellipse->V[0] = 1.0;
    Ellipse->V[1] = 0.0;
    Ellipse->V[2] = 0.0;
    Ellipse->V[3] = 1.0;
  }
  else {
    V1norm = sqrt((L1-syy)*(L1-syy)+sxy*sxy);
    V2norm = sqrt((L2-syy)*(L2-syy)+sxy*sxy);
    Ellipse->V[0] = (L1-syy)/V1norm;
    Ellipse->V[1] = (L2-syy)/V2norm;
    Ellipse->V[2] = sxy/V1norm;
    Ellipse->V[3] = sxy/V2norm;
  }
  //inverse rotation matrix
  mg_lu(2, Ellipse->V, L, U);
  mg_luinv(2, L, U, Linv, Uinv);
  mg_mxm(2, 2, 2, Uinv, Linv, Ellipse->Vinv);
  //get semi-axis lengths
  mg_mxm(2, 2, 2, Ellipse->Vinv, M, Mt);
  for (i = 0; i < 4; i++) Mt[i] /= SQRT3;
  t1 = atan(-Mt[0]*SQRT3/(2.0*Mt[1]-Mt[0]));
  t2 = atan(-Mt[2]*SQRT3/(2.0*Mt[3]-Mt[2]));
  Kt[0] = cos(t1)-sin(t1)/SQRT3;
  Kt[1] = cos(t2)-sin(t2)/SQRT3;
  Kt[2] = 2.0*sin(t1)/SQRT3;
  Kt[3] = 2.0*sin(t2)/SQRT3;
  mg_mxm(2, 2, 2, M, Kt, Cp);
  for (i = 0; i < 4; i++) Cp[i] /= SQRT3;
  mg_mxm(2, 2, 2, Ellipse->Vinv, Cp, Kt);
  Ellipse->rho[0] = max(fabs(Kt[1]), fabs(Kt[3]));
  Ellipse->rho[1] = max(fabs(Kt[0]), fabs(Kt[2]));
  
  return err_OK;
}

/******************************************************************/
/* function: mg_inside_ellipse */
bool
mg_inside_ellipse(double *coord, mg_Ellipse *Ellipse)
{
  bool inside;
  double coord_p[2], ct[2], r;
  
  //shift coordinates by ellipse's center
  coord_p[0] = coord[0]-Ellipse->Ot[0];
  coord_p[1] = coord[1]-Ellipse->Ot[1];
  
  //rotate reference
  mg_mxm(2, 2, 1, Ellipse->Vinv, coord_p, ct);
  
  //scale directions
  ct[0] /= Ellipse->rho[0];
  ct[1] /= Ellipse->rho[1];
  
  r = sqrt(ct[0]*ct[0]+ct[1]*ct[1]);
  
  inside = (r <= 1.0);
  
  return inside;
}

/******************************************************************/
/* function: mg_ellipse_frm_face_p */
/* creates a steiner ellipse using a front face and a point */
int mg_ellipse_frm_face_p(mg_Mesh *Mesh, mg_FaceData *face,
                           double *Popt, mg_Ellipse *Ellipse)
{
  int ierr, dim = Mesh->Dim, i, d;
  double coord[6];
  int *node = face->node;
  
  for (i = 0; i < 2; i++)
    for (d = 0; d < dim; d++)
      coord[i*dim+d] = Mesh->Coord[node[i]*dim+d];
  coord[4] = Popt[0];
  coord[5] = Popt[1];
  
  call(mg_circumellipse(coord, Ellipse));
  
  return err_OK;
}

