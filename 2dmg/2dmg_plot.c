//
//  2dmg_plot.c
//  2dmg
//
//  Created by Marco Ceze on 11/16/14.
//  https://github.com/mceze/2dmg
//

#include "2dmg_plot.h"
#include "2dmg_utils.h"
#include "2dmg_math.h"
#include "2dmg_def.h"


/******************************************************************/
/* function:  mg_init_plot_mesh */
void mg_init_plot_mesh(mg_Mesh *Mesh, mg_MeshPlot *PlotInfo)
{
  int n, dim = Mesh->Dim;
  double xlim[2], ylim[2];
  
  plsetopt("geometry","800x800");
  plsdev("xwin");
  // Initialize plplot
  plinit();
  plfont( 1 );
  pladv( 0 );
  plvsta();
  //set viewport as full window -10%
  plvpor( 0.05, 0.95, 0.05, 0.95 );
  plvasp(1.0);
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
  plscolbg(110,110,110);
  
  PlotInfo->range[0] = xlim[0];
  PlotInfo->range[1] = xlim[1];
  PlotInfo->range[2] = ylim[0];
  PlotInfo->range[3] = ylim[1];
  PlotInfo->tree_on    = false;
  PlotInfo->ellipse_on = false;
  PlotInfo->elem_on    = false;
  
  //set plot limits
  plwind(xlim[0], xlim[1], ylim[0], ylim[1] );
  plschr(0.0,0.5);
}

/******************************************************************/
/* function:  mg_plot_ellipse */
int mg_plot_ellipse(mg_Ellipse *ellipse)
{
  int ierr, i;
  int n = 50;
  double xplt[50], yplt[50];
  double cosval[50] = {1.0000000000e+00, 9.9179001382e-01,\
    9.6729486304e-01, 9.2691675735e-01, 8.7131870412e-01, \
    8.0141362187e-01, 7.1834935010e-01, 6.2348980186e-01, \
    5.1839256831e-01, 4.0478334312e-01, 2.8452758663e-01, \
    1.5959989503e-01, 3.2051577572e-02, -9.6023025908e-02,\
    -2.2252093396e-01, -3.4536505442e-01, -4.6253829024e-01,\
    -5.7211666012e-01, -6.7230089026e-01, -7.6144595837e-01,\
    -8.3808810489e-01, -9.0096886790e-01, -9.4905574701e-01,\
    -9.8155915699e-01, -9.9794539275e-01, -9.9794539275e-01,\
    -9.8155915699e-01, -9.4905574701e-01, -9.0096886790e-01,\
    -8.3808810489e-01, -7.6144595837e-01, -6.7230089026e-01,\
    -5.7211666012e-01, -4.6253829024e-01, -3.4536505442e-01,\
    -2.2252093396e-01, -9.6023025908e-02, 3.2051577572e-02, \
    1.5959989503e-01, 2.8452758663e-01, 4.0478334312e-01,\
    5.1839256831e-01, 6.2348980186e-01, 7.1834935010e-01,\
    8.0141362187e-01, 8.7131870412e-01, 9.2691675735e-01,\
    9.6729486304e-01, 9.9179001382e-01, 1.0000000000e+00};
  double sinval[50] = {0.0000000000e+00, 1.2787716168e-01,\
    2.5365458391e-01, 3.7526700488e-01, 4.9071755200e-01,\
    5.9811053049e-01, 6.9568255060e-01, 7.8183148247e-01,\
    8.5514276301e-01, 9.1441262302e-01, 9.5866785304e-01,\
    9.8718178341e-01, 9.9948621620e-01, 9.9537911295e-01,\
    9.7492791218e-01, 9.3846842205e-01, 8.8659930637e-01,\
    8.2017225460e-01, 7.4027799708e-01, 6.4822839531e-01,\
    5.4553490121e-01, 4.3388373912e-01, 3.1510821802e-01,\
    1.9115862870e-01, 6.4070219981e-02, -6.4070219981e-02,\
    -1.9115862870e-01, -3.1510821802e-01, -4.3388373912e-01,\
    -5.4553490121e-01, -6.4822839531e-01, -7.4027799708e-01,\
    -8.2017225460e-01, -8.8659930637e-01, -9.3846842205e-01,\
    -9.7492791218e-01, -9.9537911295e-01, -9.9948621620e-01,\
    -9.8718178341e-01, -9.5866785304e-01, -9.1441262302e-01,\
    -8.5514276301e-01, -7.8183148247e-01, -6.9568255060e-01,\
    -5.9811053049e-01, -4.9071755200e-01, -3.7526700488e-01,\
    -2.5365458391e-01, -1.2787716168e-01, -2.4492935983e-16};
  
  double R[4], S[4], T[4], x[2], xp[2];
  //rotation
  R[0] = ellipse->V[0];
  R[1] = -ellipse->V[2];
  R[2] = ellipse->V[2];
  R[3] = ellipse->V[0];
  //scaling
  S[0] = ellipse->rho[0];
  S[1] = S[2] = 0.0;
  S[3] = ellipse->rho[1];
  mg_mxm(2, 2, 2, R, S, T);
  
  for (i = 0; i < n; i++) {
    x[0] = cosval[i];
    x[1] = sinval[i];
    mg_mxm(2, 2, 1, T, x, xp);
    xplt[i] = xp[0]+ellipse->Ot[0];
    yplt[i] = xp[1]+ellipse->Ot[1];
  }
  plcol0( 3 );
  plline( n, xplt, yplt );
  
  return err_OK;
}

/******************************************************************/
/* function:  mg_plot_mesh */
int mg_plot_mesh(mg_Mesh *Mesh, mg_MeshPlot *PlotInfo)
{
  int ierr, f, node0, node1, dim = Mesh->Dim, n, e, in, node;
  double x[3], y[3], xlim[2],ylim[2], coord[6];
  char text[10];
  mg_Ellipse Ellipse;
  
  plflush();
  plclear();
  plcol0( 1 );
  
  xlim[0] = PlotInfo->range[0];
  xlim[1] = PlotInfo->range[1];
  ylim[0] = PlotInfo->range[2];
  ylim[1] = PlotInfo->range[3];
  
  //set plot limits
  plwind(xlim[0], xlim[1], ylim[0], ylim[1] );
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
  
  if (PlotInfo->tree_on){
    if (Mesh->QuadTree != NULL)
      call(mg_plot_branch(Mesh->QuadTree));
  }
  
  if (PlotInfo->ellipse_on) {
    for (e = 0; e < Mesh->nElem; e++){
      for (in = 0; in < Mesh->Elem[e].nNode; in++) {
        node = Mesh->Elem[e].node[in];
        coord[in*Mesh->Dim] = Mesh->Coord[node*Mesh->Dim];
        coord[in*Mesh->Dim+1] = Mesh->Coord[node*Mesh->Dim+1];
      }
      call(mg_circumellipse(coord, &Ellipse));
      call(mg_plot_ellipse(&Ellipse));
    }
  }
  
  if (PlotInfo->elem_on) {
    plcol0( 4 );
    for (e = 0; e < Mesh->nElem; e++){
      coord[0] = 0.0;
      coord[1] = 0.0;
      for (in = 0; in < Mesh->Elem[e].nNode; in++) {
        node = Mesh->Elem[e].node[in];
        coord[0] += Mesh->Coord[node*Mesh->Dim]/Mesh->Elem[e].nNode;
        coord[1] += Mesh->Coord[node*Mesh->Dim+1]/Mesh->Elem[e].nNode;
      }
      sprintf(text, "%d",e);
      plptex(coord[0],coord[1],0.0,0.0,0.5,text);
    }
  }
  
  return err_OK;
}

/******************************************************************/
/* function:  mg_close_plot_mesh */
void mg_close_plot_mesh(void)
{
    plend();
}

/******************************************************************/
/* function:  mg_show_mesh */
int mg_show_mesh(mg_Mesh *Mesh)
{
  int ierr, n;
  bool open = true, refresh = false;
  double range[4];
  int cursorval;
  mg_MeshPlot PlotInfo;
  
  static PLGraphicsIn gin;
  
  mg_init_plot_mesh(Mesh, &PlotInfo);
  
  //plot mesh
  call(mg_plot_mesh(Mesh, &PlotInfo));
  
  while (open) {
//    cursorval = plGetCursor( &gin );
//    printf("cursorval = %d type: %d wx: %1.2e wy: %1.2e dx: %1.2e dy: %1.2e key: %s button: %d\n",
//           cursorval, gin.type,gin.wX,gin.wY,gin.dX,gin.dY,gin.string, gin.button);
    //zoom
    if (strcmp(gin.string,"z")==0){
      while (1) {
        printf("Left click top left corner of zoom window\n");
        cursorval = plGetCursor( &gin );
        if (gin.button == 1){
          PlotInfo.range[0] = gin.wX;
          PlotInfo.range[3] = gin.wY;
          while (1) {
            printf("Left click bottom right corner of zoom window\n");
            cursorval = plGetCursor( &gin );
            if (gin.button == 1){
              PlotInfo.range[1] = gin.wX;
              PlotInfo.range[2] = gin.wY;
              break;
            }
          }
          break;
        }
      }
      refresh = true;
    }
    //plot qtree
    if (strcmp(gin.string,"t")==0){
      PlotInfo.tree_on = !PlotInfo.tree_on;
      refresh = true;
    }
    
    //plot ellipses
    if (strcmp(gin.string,"e")==0){
      PlotInfo.ellipse_on = !PlotInfo.ellipse_on;
      refresh = true;
    }
    
    //show elem numbers
    if (strcmp(gin.string,"s")==0){
      PlotInfo.elem_on = !PlotInfo.elem_on;
      refresh = true;
    }
    
    //reset range
    if (strcmp(gin.string,"r")==0){
      for (n = 1; n < Mesh->nNode; n++){
        if (Mesh->Coord[n*Mesh->Dim+0] < PlotInfo.range[0])
          PlotInfo.range[0] = Mesh->Coord[n*Mesh->Dim+0];
        if (Mesh->Coord[n*Mesh->Dim+0] > PlotInfo.range[1])
          PlotInfo.range[1] = Mesh->Coord[n*Mesh->Dim+0];
        if (Mesh->Coord[n*Mesh->Dim+1] < PlotInfo.range[2])
          PlotInfo.range[2] = Mesh->Coord[n*Mesh->Dim+1];
        if (Mesh->Coord[n*Mesh->Dim+1] > PlotInfo.range[3])
          PlotInfo.range[3] = Mesh->Coord[n*Mesh->Dim+1];
      }
      refresh = true;
    }
    //quit
    if (strcmp(gin.string,"q")==0){
      break;
    }
    
    if (refresh){
      call(mg_plot_mesh(Mesh, &PlotInfo));
      refresh = false;
    }
  }
  
  mg_close_plot_mesh();
  
  return err_OK;
}

