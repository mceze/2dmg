//
//  2dmg_error.h
//  2dmg
//
//  Created by Marco Ceze on 1/13/15.
//  
//

#ifndef ___dmg__error__
#define ___dmg__error__

#include <stdio.h>

/******************************************************************/
/* Error codes */
#define err_OK               0
#define err_READWRITE_ERROR -1
#define err_INCONSISTENCY   -2
#define err_INPUT_ERROR     -3
#define err_INCOMPATIBLE    -4
#define err_MEMORY_ERROR    -5
#define err_NOT_FOUND       -6
#define err_OUT_OF_BOUNDS   -7
#define err_LOGIC_ERROR     -8
#define err_NOT_SUPPORTED   -9
#define err_MESH_ERROR      -10
#define err_HSEARCH_ERROR   -11
#define err_GSL_ERROR       -12
#define err_SINGULAR        -13

/******************************************************************/
/* Error Macro: used to report error occurrences */
#define error(X) (error_report( __FILE__, __LINE__, #X, (X)))
#define call(X) do{ierr = error(X); if(ierr != err_OK) return ierr;} while(0)

/******************************************************************/
/* function:  error_report*/
/* reports the file, line and call  where the error happened*/
int error_report( char *file, int line, char *call, int ierr);

#endif /* defined(___dmg__error__) */
