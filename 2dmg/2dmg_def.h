//
//  2dmg_def.h
//  2dmg
//
//  Created by Marco Ceze on 11/11/14.
//  Copyright (c) 2014 Marco Ceze. All rights reserved.
//

#ifndef _dmg__dmg_def_h
#define _dmg__dmg_def_h

/******************************************************************/
/* standard headers used */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stddef.h>
#include <math.h>
#include <stdbool.h>
#include <errno.h>
#include <ctype.h>
#include "2dmg_io.h"
#include "2dmg_error.h"

#define _GNU_SOURCE
#include <search.h>

/******************************************************************/
/* constant macros */
#define MEPS              1e-16 // machine precision
#define MAXLINELEN        400 //maximum line length
#define MAXSTRLEN         100 //maximum string length
#define MAXLONGLINELEN    1000  // maximum long line length, in characters
#define MAXFACENNODE      25
#define HOLLOWNEIGHTAG   -2147483648 //minimum integer value
#define LEFTNEIGHINDEX    0 //left tag
#define RIGHTNEIGHINDEX   1 //right tag
#define HALFSQRT3         0.866025403784439
#define SQRT3             1.73205080756888
#define NPARAMLIST        100 //hash table size

/******************************************************************/
/* Useful macros */
#define max(a,b)  (((a) > (b)) ? (a) : (b))
#define min(a,b)  (((a) < (b)) ? (a) : (b))
#define sign(a) (((a) < 0.0) ? -1 : 1)
#define swap(a,b, t)  {t = a; a = b; b = t;}
#define cycle3(a,b,c, t)  {t = a; a = b; b = c; c = t;}
#define cycle4(a,b,c,d, t)  {t = a; a = b; b = c; c = d; d = t;}
#define nextincycle(i,ni)  i+1-(int)(i/(ni-1))*(i+1)
#define previncycle(i,ni)  i-1+(1-(int)((ni-1+i)/ni))*ni;

#endif
