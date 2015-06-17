//
//  2dmg_math.c
//  2dmg
//
//  Created by Marco Ceze on 11/11/14.
//  Copyright (c) 2014 Marco Ceze. All rights reserved.
//

#include "2dmg_math.h"
#include "2dmg_def.h"
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
/* function: mg_metric_dist */
/* computes metric distance between 2 points */
int mg_metric_dist(mg_Metric *Metric, int order, double *coord,
                   double *dist)
{
  gsl_integration_glfixed_table *gltable;
  double dx = coord[1]-coord[0], dy = coord[3]-coord[2];
  double ab[2], M[4], x, y, ds2, xgl, wgl;
  ab[0] = dx;
  ab[1] = dy;
  int ip;
  
  if ((gltable = gsl_integration_glfixed_table_alloc(order)) == NULL)
    return error(err_GSL_ERROR);
  
  switch (Metric->type) {
    case mge_Metric_Analitic1:
      (*dist) = 0.0;
      for (ip = 0; ip < gltable->n; ip++) {
        gsl_integration_glfixed_point(0.0, 1.0, ip, &xgl, &wgl, gltable);
        x = coord[0]+xgl*dx;
        y = coord[2]+xgl*dy;
        M[0] = 1.0+exp(-pow((x-0.5)/(0.1), 2.0)); M[1] = 0.0;
        M[2] = 0.0; M[3] = 1.0;
        ds2   = ab[0]*(M[0]*ab[0]+M[1]*ab[1])+ab[1]*(M[2]*ab[0]+M[3]*ab[1]);
        (*dist) += wgl*sqrt(ds2);
      }
      break;
    default:
      error(err_NOT_SUPPORTED);
      break;
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
  
  
  return err_OK;
}

