//
//  2dmg_qtree.c
//  2dmg
//
//  Created by Marco Ceze on 12/19/14.
//  
//

#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <memory.h>
#include "2dmg_qtree.h"
#include "2dmg_error.h"
#include "plConfig.h"
#include "plplot.h"

/******************************************************************/
/* function:  mg_init_branch */
int mg_init_branch(mg_qtree *branch)
{
  int i;
  
  branch->c[0] = branch->c[1] = 0.0;
  branch->ds[0] = branch->ds[1] = 0.0;
  branch->capacity = MAXCAPACITY;
  branch->n_entry = 0;
  //leave allocated space
  if ((branch->data = malloc(branch->capacity*sizeof(mg_data_entry))) == NULL)
    return error(err_MEMORY_ERROR);
  //loop over children and set them as NULL
  for (i = 0; i < 4; i++)
    branch->child[i] = NULL;
  
  return err_OK;
}

/******************************************************************/
/* function:  mg_destroy_branch */
void mg_destroy_branch(mg_qtree *branch)
{
  int i;
  //free up children first
  if (branch->child[0] != NULL){
    for (i = 0; i < 4; i++) {
      mg_destroy_branch(branch->child[i]);
    }
  }
  //free data storage
  free(branch->data);
  branch->data = NULL;
}

/******************************************************************/
/* function:  mg_branch_qtree */
int mg_branch_qtree(mg_qtree *qtree)
{
  int i, ierr;
  double wx[4] = {-1.,1.,-1.,1.}, wy[4] = {-1.,-1.,1.,1.};
  mg_qtree *branch;
  
  for (i = 0; i < 4; i++) {
    branch = malloc(sizeof(mg_qtree));
    call(mg_init_branch(branch));
    qtree->child[i] = branch;
    //assign half dimensions and centroids
    qtree->child[i]->ds[0] = 0.5*qtree->ds[0];
    qtree->child[i]->ds[1] = 0.5*qtree->ds[1];
    qtree->child[i]->c[0] = qtree->c[0]+wx[i]*qtree->child[i]->ds[0];
    qtree->child[i]->c[1] = qtree->c[1]+wy[i]*qtree->child[i]->ds[1];
  }
  
  return err_OK;
}

/******************************************************************/
/* function:  mg_add_qtree_entry */
int mg_add_qtree_entry(double coord[2], void **data, mg_qtree *qtree)
{
  int ierr, i, quad;
  
  //first check if within bounds
  for (i = 0; i < 2; i++)
    if (coord[i] < qtree->c[i]-qtree->ds[i] ||
        coord[i] > qtree->c[i]+qtree->ds[i])
      return error(err_OUT_OF_BOUNDS);
  
  
  if (qtree->n_entry+1 <= qtree->capacity &&
      qtree->child[0] == NULL) {
    //add entry to bin
    qtree->n_entry++;
    qtree->data[qtree->n_entry-1].coord[0] = coord[0];
    qtree->data[qtree->n_entry-1].coord[1] = coord[1];
    //store only pointer
    qtree->data[qtree->n_entry-1].data = data;
  }
  else {
    if (qtree->child[0] == NULL){
      //branch out children
      call(mg_branch_qtree(qtree));
    }
    //add new entry one of the children
    quad = quadrant(coord, qtree->c);
    call(mg_add_qtree_entry(coord, data, qtree->child[quad]));
    //loop over old entries and redistribute amongst the children
    for (i = 0; i < qtree->n_entry; i++){
      quad = quadrant(qtree->data[i].coord, qtree->c);
      call(mg_add_qtree_entry(qtree->data[i].coord, &qtree->data[i].data,
                              qtree->child[quad]));
    }
    qtree->n_entry = 0;
  }
  
  
  return err_OK;
}

/******************************************************************/
/* function:  mg_find_branch */
int mg_find_branch(mg_qtree *trunk, mg_qtree **pbranch, double coord[2])
{
  int ierr, i, quad;
  mg_qtree *new_trunk;

  //first check if within bounds
  for (i = 0; i < 2; i++)
    if (coord[i] < trunk->c[i]-trunk->ds[i] ||
        coord[i] > trunk->c[i]+trunk->ds[i])
      return error(err_OUT_OF_BOUNDS);
  //check if trunk has branches
  if (trunk->child[0] == NULL)
    (*pbranch) = trunk; //got to minimal subdivision
  else {
    quad = quadrant(coord, trunk->c);
    new_trunk = trunk->child[quad];
    call(mg_find_branch(new_trunk, pbranch, coord));
  }
  
  return err_OK;
}

/******************************************************************/
/* function:  mg_find_entry */
int mg_find_entry(mg_qtree *trunk, double coord[2], void *target_data,
                  double tol, mg_data_entry **data)
{
  int ierr, i;
  bool found = false;
  mg_qtree *branch;
  
  //find branch containing coord
  call(mg_find_branch(trunk, &branch, coord));

  //check if match will be exact or by distance
  if (target_data != NULL ) {
    //looking for exact match by pointer
    for (i = 0; i < branch->n_entry; i++) {
      if (branch->data[i].data == target_data) {
        (*data) = branch->data+i;
        found = true;
        break;
      }
    }
  }
  else {
    for (i = 0; i < branch->n_entry; i++) {
      if (distance(coord, branch->data[i].coord) <= tol) {
        (*data) = branch->data+i;
        found = true;
        break;
      }
    }
  }
  if (!found) return err_NOT_FOUND;
  
  return err_OK;
}

/******************************************************************/
/* function:  mg_del_entry */
int mg_del_entry(mg_qtree *branch, void *data)
{
  int i;
  
  for (i = 0; i < branch->n_entry; i++) {
    if (branch->data[i].data == data) {
      branch->data[i].data = NULL;
      break;
    }
  }
  //ensure continuity of data array
  memmove(branch->data+i,branch->data+i+1,i+1-branch->n_entry);
  branch->n_entry--;
  
  return err_OK;
}

/******************************************************************/
/* function:  mg_plot_branch */
int mg_plot_branch(mg_qtree *branch)
{
  int i, ierr;
  double x[5], y[5], xp[MAXCAPACITY], yp[MAXCAPACITY];
  
  x[0] = branch->c[0]-branch->ds[0];
  x[1] = branch->c[0]+branch->ds[0];
  x[2] = branch->c[0]+branch->ds[0];
  x[3] = branch->c[0]-branch->ds[0];
  x[4] = branch->c[0]-branch->ds[0];
  
  y[0] = branch->c[1]-branch->ds[1];
  y[1] = branch->c[1]-branch->ds[1];
  y[2] = branch->c[1]+branch->ds[1];
  y[3] = branch->c[1]+branch->ds[1];
  y[4] = branch->c[1]-branch->ds[1];
  
  plcol0( 2 );
  plline( 5, x, y );
  
  //plot points
  for (i = 0; i < branch->n_entry; i++) {
    xp[i] = branch->data[i].coord[0];
    yp[i] = branch->data[i].coord[1];
  }
  if (branch->n_entry > 0){
    plcol0( 1 );
    plpoin(branch->n_entry, xp, yp, 22);
  }
  if (branch->child[0] != NULL)
    for (i = 0; i < 4; i++)
      call(mg_plot_branch(branch->child[i]));
  
  
  return err_OK;
}

/******************************************************************/
/* function:  mg_plot_qtree */
int mg_plot_qtree(mg_qtree *qtree, int strm)
{
  int ierr;
  double xlim[2], ylim[2];
  
  plsdev("aqt");
  // Initialize plplot
  plsstrm (	strm);
  plinit();
  plfont( 2 );
  pladv( 0 );
//  plbop (	);
  //set viewport as full window -10%
  plvpor( 0.05, 0.95, 0.05, 0.95 );
  //set plot limits
  xlim[0] = qtree->c[0]-qtree->ds[0];
  xlim[1] = qtree->c[0]+qtree->ds[0];
  ylim[0] = qtree->c[1]-qtree->ds[1];
  ylim[1] = qtree->c[1]+qtree->ds[1];
  
  plwind(xlim[0], xlim[1], ylim[0], ylim[1] );
//  plenv0(xlim[0], xlim[1], ylim[0], ylim[1],2,-2);
  
  
  call(mg_plot_branch(qtree));
  
  plend();
  //pleop (	);
  //plspause( 1 );

  return err_OK;
}


