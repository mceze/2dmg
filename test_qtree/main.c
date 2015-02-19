//
//  main.c
//  test_qtree
//
//  Created by Marco Ceze on 1/15/15.
//  Copyright (c) 2015 Marco Ceze. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include "2dmg_qtree.h"
#include "2dmg_error.h"


int main(int argc, const char * argv[]) {
  int ierr;
  double x[2];
  mg_qtree *qtree, *branch;
  mg_data_entry *D;
  
  qtree = malloc(sizeof(mg_qtree));
  call(mg_init_branch(qtree));
  qtree->c[0] = 0.0;
  qtree->c[1] = 0.0;
  qtree->ds[0] = 1.0;
  qtree->ds[1] = 1.0;
  
  x[0] = 0.2; x[1] = 0.3;
  call(mg_add_qtree_entry(x, NULL, qtree));
  x[0] = 0.2; x[1] = 0.35;
  call(mg_add_qtree_entry(x, NULL, qtree));
  call(mg_plot_qtree(qtree,0));
  
  x[0] = 0.4; x[1] = 0.2;
  call(mg_add_qtree_entry(x, NULL, qtree));
  call(mg_plot_qtree(qtree,0));
  
  x[0] = 0.2; x[1] = -0.3;
  call(mg_add_qtree_entry(x, NULL, qtree));
//  call(mg_plot_qtree(qtree));
  
  x[0] = 0.25; x[1] = 0.6;
  call(mg_add_qtree_entry(x, NULL, qtree));
//  call(mg_plot_qtree(qtree));
  x[0] = -0.2; x[1] = -0.3;
  call(mg_add_qtree_entry(x, NULL, qtree));
//  call(mg_plot_qtree(qtree));
  x[0] = 0.8; x[1] = 0.3;
  call(mg_add_qtree_entry(x, (void*)qtree, qtree));
  //call(mg_plot_qtree(qtree));
  x[0] = 0.21; x[1] = 0.3;
  call(mg_find_branch(qtree, &branch, x));
  x[0] = 0.8; x[1] = 0.3;
  call(mg_find_entry(qtree, x, qtree, 1e-12, &D));
  branch = (mg_qtree*)&(D->data[0]);
  
  mg_destroy_branch(qtree);
  free(qtree);
  
  return err_OK;
}
