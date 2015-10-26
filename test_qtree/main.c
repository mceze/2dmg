//
//  main.c
//  test_qtree
//
//  Created by Marco Ceze on 1/15/15.
//  
//

#include <stdio.h>
#include <stdlib.h>
#include "2dmg_qtree.h"
#include "2dmg_error.h"
#include "plConfig.h"
#include "plplot.h"


int main(int argc, const char * argv[]) {
  int ierr, np = 100, i;
  double x[2],r;
  mg_qtree *qtree, *branch;
  mg_data_entry *D;
  
  qtree = malloc(sizeof(mg_qtree));
  call(mg_init_branch(qtree));
  qtree->c[0] = 0.0;
  qtree->c[1] = 0.0;
  qtree->ds[0] = 6.0;
  qtree->ds[1] = 6.0;
  
  call(mg_plot_qtree(qtree,0));
  
  for (i = 0; i < np; i++) {
    x[0] = x[1] = 10.0;
    while (fabs(x[0])>=qtree->ds[0] || fabs(x[1])>=qtree->ds[1]){
      r = 10.0;
      while (r > 1.0){
        x[0] = ((rand()%10000)/10000.0-0.5)*2.0;
        x[1] = ((rand()%10000)/10000.0-0.5)*2.0;
        r = x[0]*x[0] + x[1]*x[1];
      }
      x[0] = x[0]*sqrt(-2.0*log(r)/r);
      x[1] = x[1]*sqrt(-2.0*log(r)/r);
    }
    call(mg_add_qtree_entry(x, NULL, qtree));
    call(mg_plot_branch(qtree));
  }
  //plot tree
  //call(mg_plot_qtree(qtree,0));
  plend();
  
  mg_destroy_branch(qtree);
  free(qtree);
  
  return err_OK;
}
