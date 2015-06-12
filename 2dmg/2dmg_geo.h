//
//  2dmg_geo.h
//  2dmg
//
//  Created by Marco Ceze on 6/10/15.
//  Copyright (c) 2015 Marco Ceze. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include "2dmg_def.h"
#include "2dmg_struct.h"
#include "2dmg_utils.h"
#include "2dmg_math.h"
#include "2dmg_plot.h"

#ifndef _dmg__dmg_geo_h
#define _dmg__dmg_geo_h

/******************************************************************/
/* function:  mg_create_geo */
/* creates a mg_Geometry structure with "nBoundary" boundaries  */
int mg_create_geo(mg_Geometry **pGeo, int nBoundary, int nPoint,
                  int Dim);

/******************************************************************/
/* function:  mg_destroy_geo */
/* destroys a mg_Geometry structure with "nBoundary" boundaries  */
void mg_destroy_geo(mg_Geometry *Geo);

/******************************************************************/
/* function:  mg_init_segment */
/* initializes the interpolant of a mg_Segment  */
int mg_init_segment(mg_Geometry *Geo, int iseg);

#endif
