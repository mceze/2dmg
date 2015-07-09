//
//  2dmg_metric.h
//  2dmg
//
//  Created by Marco Ceze on 3/1/15.
//  Copyright (c) 2015 Marco Ceze. All rights reserved.
//

#ifndef _dmg__dmg_metric_h
#define _dmg__dmg_metric_h

#include "2dmg_struct.h"

/******************************************************************/
/* metric field */
enum mg_MetricType {
  mg_MetricAnalitycTest,
  mg_MetricFEData,
  mg_MetricLast
};
static char *mg_MetricTypeName[mg_MetricLast] = {
  "AnalyticTest",
  "FEData"
};

typedef struct
{
  enum mg_MetricType type;
  bool isotropic;
  double *M;
  struct mg_Mesh *Mesh;
}
mg_MetricField;

#endif
