//
//  2dmg_metric_struct.h
//  2dmg
//
//  Created by Marco Ceze on 6/26/15.
//  Copyright (c) 2015 Marco Ceze. All rights reserved.
//

#ifndef _dmg__dmg_metric_struct_h
#define _dmg__dmg_metric_struct_h

#include <gsl/gsl_interp.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas.h>
#include "2dmg_struct.h"

/******************************************************************/
/* enumerators for metric types */
enum mge_Metric {
  mge_Metric_Uniform,
  mge_Metric_Analitic1,
  mge_Metric_Analitic2,
  mge_Metric_Analitic3,
  mge_Metric_Last
};
static char *mge_MetricName[mge_Metric_Last] = {
  "MetricUniform",
  "MetricAnalytic1",
  "MetricAnalytic2",
  "MetricAnalytic3"
};

/******************************************************************/
/* mesh structure */
typedef struct
{
  enum mge_Metric type;
  mg_Mesh *BGMesh;
  int order; //interpolation order (Lagrange basis)
  double *M;
}
mg_Metric;

/******************************************************************/
/* structure: mg_gsl_multimin_params */
/* parameters for optimizer*/
typedef struct
{
  mg_Segment *Segment;
  mg_Metric *Metric;
  double *scale;
  
}
mg_gsl_multimin_params;

#endif
