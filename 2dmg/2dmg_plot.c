//
//  2dmg_plot.c
//  2dmg
//
//  Created by Marco Ceze on 11/16/14.
//  Copyright (c) 2014 Marco Ceze. All rights reserved.
//

#include "2dmg_plot.h"
#include "2dmg_def.h"

/******************************************************************/
/* function:  mg_plot_mesh */
int mg_plot_mesh(mg_Mesh *Mesh)
{
  int f, node0, node1, n, dim = Mesh->Dim;
  double x[3], y[3], xlim[2], ylim[2];
  
  plsdev("xwin");
  // Initialize plplot
  plinit();
  plfont( 2 );
  pladv( 0 );
  //set viewport as full window -10%
  plvpor( 0.05, 0.95, 0.05, 0.95 );
  //get xmin, xmax, ymin, ymax
  xlim[0] = xlim[1] = Mesh->Coord[0];
  ylim[0] = ylim[1] = Mesh->Coord[1];
  for (n = 1; n < Mesh->nNode; n++){
    if (Mesh->Coord[n*dim+0] < xlim[0])
      xlim[0] = Mesh->Coord[n*dim+0];
    if (Mesh->Coord[n*dim+0] > xlim[1])
      xlim[1] = Mesh->Coord[n*dim+0];
    if (Mesh->Coord[n*dim+1] < ylim[0])
      ylim[0] = Mesh->Coord[n*dim+1];
    if (Mesh->Coord[n*dim+1] > ylim[1])
      ylim[1] = Mesh->Coord[n*dim+1];
  }
//  xlim[0] = -0.1; xlim[1] = 1.1;
//  ylim[0] = -0.55; ylim[1] = 0.55;
  //set plot limits
  plwind(xlim[0], xlim[1], ylim[0], ylim[1] );
  plcol0( 1 );
  //plot faces
  for (f = 0; f < Mesh->nFace; f++) {
    node0 = Mesh->Face[f]->node[0];
    node1 = Mesh->Face[f]->node[1];
    x[0] = Mesh->Coord[node0*dim];
    x[1] = Mesh->Coord[node1*dim];
    y[0] = Mesh->Coord[node0*dim+1];
    y[1] = Mesh->Coord[node1*dim+1];
    plline( 2, x, y );
  }
  
  plend();
  
  return err_OK;
}
