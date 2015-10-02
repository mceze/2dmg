//
//  main.c
//  test_interp
//
//  Created by Marco Ceze on 6/10/15.
//  
//

#include <gsl/gsl_interp.h>
#include <stdio.h>
#include "2dmg_def.h"
#include "2dmg_utils.h"
#include "2dmg_struct.h"
#include "2dmg_geo.h"
#include "2dmg_io.h"

int main(int argc, const char * argv[]) {
  int ierr, d,b, i, np = 15;
  double xi,yi,dxi,dyi,t=.25, coord[]={0.0,1.0,0.0,1.0}, dist;
  double *tref, scale;
  mg_Geometry *Geo;
  mg_Metric *Metric;
  
  call(mg_read_geo(&Geo, "box.geo"));
  b=0;
  d=0;
  xi = gsl_interp_eval(Geo->Boundary[b]->interp[d],
                       Geo->Boundary[b]->s,
                       Geo->Boundary[b]->Coord+d*Geo->Boundary[b]->nPoint, t,
                       Geo->Boundary[b]->accel[d]);
  dxi = gsl_interp_eval_deriv(Geo->Boundary[b]->interp[d],
                              Geo->Boundary[b]->s,
                              Geo->Boundary[b]->Coord+d*Geo->Boundary[b]->nPoint, t,
                              Geo->Boundary[b]->accel[d]);
  
  d=1;
  yi = gsl_interp_eval(Geo->Boundary[b]->interp[d],
                       Geo->Boundary[b]->s,
                       Geo->Boundary[b]->Coord+d*Geo->Boundary[b]->nPoint, t,
                       Geo->Boundary[b]->accel[d]);
  dyi = gsl_interp_eval_deriv(Geo->Boundary[b]->interp[d],
                              Geo->Boundary[b]->s,
                              Geo->Boundary[b]->Coord+d*Geo->Boundary[b]->nPoint, t,
                              Geo->Boundary[b]->accel[d]);
  
  Metric = malloc(sizeof(mg_Metric));
  Metric->type = mge_Metric_Analitic2;
  Metric->order = 8;
//  Metric.type = mge_Metric_Uniform;
//  Metric.order = 1;
  call(mg_create_mesh(&Metric->BGMesh));
  Metric->BGMesh->Dim = 2;
//  call(mg_metric_dist(&Metric, 1, coord, &dist));
//  printf("dist: %1.12e\n",dist);
//  call(mg_metric_dist(&Metric, 2, coord, &dist));
//  printf("dist: %1.12e\n",dist);
//  call(mg_metric_dist(&Metric, 4, coord, &dist));
//  printf("dist: %1.12e\n",dist);
//  call(mg_metric_dist(&Metric, 8, coord, &dist));
//  printf("dist: %1.12e\n",dist);
//  call(mg_metric_dist(&Metric, 16, coord, &dist));
//  printf("dist: %1.12e\n",dist);
//  call(mg_metric_dist(&Metric, 32, coord, &dist));
//  printf("dist: %1.12e\n",dist);
//  call(mg_metric_dist(&Metric, 64, coord, &dist));
//  printf("dist: %1.12e\n",dist);

  
//  call(mg_metric_length(&Metric, Geo->Boundary[b], 16, &dist));
//  printf("dist: %1.12e\n",dist);
//  call(mg_metric_length(&Metric, Geo->Boundary[b], 32, &dist));
//  printf("dist: %1.12e\n",dist);
  call(mg_metric_length(Metric, Geo->Boundary[b], 64, &dist));
  
  call(mg_mesh_segment(Geo->Boundary[b], Metric, np, &scale, &tref));
  for (i = 0; i < np; i++) {
    t=tref[i];
    d=0;
    xi = gsl_interp_eval(Geo->Boundary[b]->interp[d],
                         Geo->Boundary[b]->s,
                         Geo->Boundary[b]->Coord+d*Geo->Boundary[b]->nPoint, t,
                         Geo->Boundary[b]->accel[d]);
    d=1;
    yi = gsl_interp_eval(Geo->Boundary[b]->interp[d],
                         Geo->Boundary[b]->s,
                         Geo->Boundary[b]->Coord+d*Geo->Boundary[b]->nPoint, t,
                         Geo->Boundary[b]->accel[d]);
    printf("%1.8e %1.8e %1.8e\n",tref[i],xi,yi);
  }

  mg_free((void*)tref);
  mg_destroy_mesh(Metric->BGMesh);
  mg_free((void*)Metric);
  mg_destroy_geo(Geo);
  
  return 0;
}
