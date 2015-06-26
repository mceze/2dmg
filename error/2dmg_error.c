//
//  2dmg_error.c
//  2dmg
//
//  Created by Marco Ceze on 1/13/15.
//  
//

#include "2dmg_error.h"

/******************************************************************/
/* function:  mg_error */
/* reports the file, line and call  where
 the error happened */
int error_report( char *file, int line, char *call, int ierr)
{
  if (ierr == err_OK) return err_OK;
  printf("**************Error: %d**************",ierr);
  printf("\nFile: %s\nLine: %d\nCall: %s\n",
         file, line, call);
  fflush(stdout);
  return ierr;
}