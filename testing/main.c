//
//  main.c
//  testing
//
//  Created by Marco Ceze on 7/24/15.
//  Copyright (c) 2015 Marco Ceze. All rights reserved.
//

#include "2dmg_def.h"
#include "2dmg_utils.h"
#include "2dmg_struct.h"
#include "2dmg_geo.h"
#include "2dmg_io.h"

int main(int argc, const char * argv[]) {
  int ierr;
  bool inside;
  double coord[6], rho[2], V[4], Vinv[4], C[2], p[2];
  mg_Ellipse Ellipse;
  
  
  coord[0] = -4.6;
  coord[1] = 0.5;
  coord[2] = -1.0;
  coord[3] = 0.0;
  coord[4] = 7.0;
  coord[5] = 5.0;
  
  call(mg_circumellipse(coord, &Ellipse));
  
  p[0] = 1.0;
  p[1] = 0.74;
  inside = mg_inside_ellipse(p, &Ellipse);
  
  
  return 0;
}
