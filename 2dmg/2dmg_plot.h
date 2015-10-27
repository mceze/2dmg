//
//  2dmg_plot.h
//  2dmg
//
//  Created by Marco Ceze on 11/16/14.
//  https://github.com/mceze/2dmg
//

#ifndef ___dmg___dmg_plot__
#define ___dmg___dmg_plot__

#include "2dmg_def.h"
#include "plConfig.h"
#include "plplot.h"
#include "plevent.h"

/******************************************************************/
/* mesh plot structure */
typedef struct
{
  double range[4]; //viewport range
  bool tree_on;
  bool ellipse_on;
  bool elem_on;
}
mg_MeshPlot;

/******************************************************************/
/* function:  mg_init_plot_mesh */
void mg_init_plot_mesh(mg_Mesh *Mesh, mg_MeshPlot *PlotInfo);

/******************************************************************/
/* function:  mg_plot_ellipse */
int mg_plot_ellipse(mg_Ellipse *ellipse);

/******************************************************************/
/* function:  mg_plot_mesh */
int mg_plot_mesh(mg_Mesh *Mesh, mg_MeshPlot *PlotInfo);

/******************************************************************/
/* function:  mg_close_plot_mesh */
void mg_close_plot_mesh(void);

/******************************************************************/
/* function:  mg_show_mesh */
int mg_show_mesh(mg_Mesh *Mesh);

#endif /* defined(___dmg___dmg_plot__) */
