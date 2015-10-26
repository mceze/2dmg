//
//  main.c
//  plot_mesh
//
//  Created by Marco Ceze on 5/26/15.
//  
//

#include <gsl/gsl_interp.h>
#include <stdio.h>
#include "2dmg_def.h"
#include "2dmg_utils.h"
#include "2dmg_io.h"
#include "2dmg_plot.h"

int main(int argc, const char * argv[])
{
  int ierr;
  mg_Mesh *Mesh;
  
  call(mg_read_mesh(&Mesh, "rae2822.gri"));
  
  call(mg_show_mesh(Mesh));
  
  mg_destroy_mesh(Mesh);
  
  return err_OK;
}
