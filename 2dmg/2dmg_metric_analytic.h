//
//  2dmg_metric_analytic.h
//  2dmg
//
//  Created by Marco Ceze on 6/17/15.
//  https://github.com/mceze/2dmg
//

#ifndef _____dmg_metric_analytic__
#define _____dmg_metric_analytic__

#include <stdio.h>
#include <stdlib.h>
#include "2dmg_def.h"
#include "2dmg_struct.h"
#include "2dmg_utils.h"
#include "2dmg_math.h"
#include "2dmg_plot.h"

/******************************************************************/
/* function:  mg_metric_uniform */
/* receives array of x,y and returns metric values at each point  */
void mg_metric_uniform(double *x, double *y, int np, double *M);

/******************************************************************/
/* function:  mg_metric_expx */
/* receives array of x,y and returns metric values at each point  */
void mg_metric_expx(double *x, double *y, int np, double *M);

/******************************************************************/
/* function:  mg_metric_sqx */
/* receives array of x,y and returns metric values at each point  */
void mg_metric_sqx(double *x, double *y, int np, double *M);

/******************************************************************/
/* function:  mg_metric_linx */
/* receives array of x,y and returns metric values at each point  */
void mg_metric_linx(double *x, double *y, int np, double *M);

#endif /* defined(_____dmg_metric_analytic__) */
