//
//  2dmg_geo.c
//  2dmg
//
//  Created by Marco Ceze on 6/10/15.
//  Copyright (c) 2015 Marco Ceze. All rights reserved.
//

#include "2dmg_geo.h"

/******************************************************************/
/* function:  mg_create_geo */
/* creates a mg_Geometry structure with "nBoundary" boundaries  */
int mg_create_geo(mg_Geometry **pGeo, int nBoundary, int nPoint,
                  int Dim)
{
  int ierr, i;
  mg_Segment *Boundary;
  
  call(mg_alloc((void**)pGeo, 1, sizeof(mg_Geometry)));
  (*pGeo)->nPoint = nPoint;
  (*pGeo)->Dim = Dim;
  call(mg_alloc((void**)&((*pGeo)->Coord), nPoint*Dim,
                sizeof(double)));
  (*pGeo)->nBoundary = nBoundary;
  call(mg_alloc((void**)&((*pGeo)->Boundary), nBoundary,
                sizeof(mg_Segment)));
  //init boundaries
  for (i = 0; i < nBoundary; i++) {
    call(mg_alloc((void**)&Boundary, 1,sizeof(mg_Segment)));
    (*pGeo)->Boundary[i] = Boundary;
    (*pGeo)->Boundary[i]->nPoint = 0;
    (*pGeo)->Boundary[i]->Coord  = NULL;
    (*pGeo)->Boundary[i]->s  = NULL;
    (*pGeo)->Boundary[i]->interp = NULL;
    (*pGeo)->Boundary[i]->accel  = NULL;
    (*pGeo)->Boundary[i]->interp_type = -1;//unitialized
    call(mg_alloc((void**)&(Boundary->Name),MAXSTRLEN,sizeof(char)));
  }
  
  return err_OK;
}

/******************************************************************/
/* function:  mg_destroy_geo */
/* destroys a mg_Geometry structure with "nBoundary" boundaries  */
void mg_destroy_geo(mg_Geometry *Geo)
{
  int i, d;
  mg_free((void*)Geo->Coord);
  for (i = 0; i < Geo->nBoundary; i++) {
    mg_free((void*)Geo->Boundary[i]->Coord);
    mg_free((void*)Geo->Boundary[i]->s);
    mg_free((void*)Geo->Boundary[i]->Name);
    if (Geo->Boundary[i]->interp != NULL){
      for (d = 0; d < Geo->Dim;d++){
        gsl_interp_free(Geo->Boundary[i]->interp[d]);
        gsl_interp_accel_free(Geo->Boundary[i]->accel[d]);
      }
      mg_free((void*)Geo->Boundary[i]->interp);
      mg_free((void*)Geo->Boundary[i]->accel);
    }
    mg_free((void*)Geo->Boundary[i]);
  }
  mg_free((void*)Geo->Boundary);
  mg_free((void*)Geo);
}

/******************************************************************/
/* function:  mg_init_segment */
/* initializes the interpolant of a mg_Segment  */
int mg_init_segment(mg_Geometry *Geo, int iseg)
{
  int ierr, i;
  mg_Segment *Seg = Geo->Boundary[iseg];
  
  switch (Seg->interp_type) {
    case mge_Linear:
      for (i = 0; i < Geo->Dim; i++) {
        Seg->interp[i] = gsl_interp_alloc(gsl_interp_linear,Seg->nPoint);
        gsl_interp_init(Seg->interp[i], Seg->s,Seg->Coord
                        +i*Seg->nPoint, Seg->nPoint);
        Seg->accel[i] = gsl_interp_accel_alloc();
      }
      break;
    case mge_CSpline:
      for (i = 0; i < Geo->Dim; i++) {
        Seg->interp[i] = gsl_interp_alloc(gsl_interp_cspline,Seg->nPoint);
        gsl_interp_init(Seg->interp[i], Seg->s,Seg->Coord
                        +(i*Seg->nPoint), Seg->nPoint);
        Seg->accel[i] = gsl_interp_accel_alloc();
      }
      break;
    default:
      return error(err_NOT_SUPPORTED);
      break;
  }
  
  
  return err_OK;
}
