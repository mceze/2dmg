//
//  2dmg_math.h
//  2dmg
//
//  Created by Marco Ceze on 11/11/14.
//  Copyright (c) 2014 Marco Ceze. All rights reserved.
//

#ifndef ___dmg___dmg_math__
#define ___dmg___dmg_math__

#include <stdio.h>


/******************************************************************/
/* function:  mg_limited_pair */
/* Pairs (x,y)=>z uniquely using ylim key.
 The maximum value for y is ylim. x >= 0; 0 =< y <= ylim.  */
int mg_limited_pair(const int x, const int y, int *z, const int ylim);

/******************************************************************/
/* function:  mg_limited_pair_inv */
/* Decodes z=>(x,y) uniquely using ylim key.
 The maximum value for y is ylim. x >= 0; 0 =< y <= ylim.  */
int mg_limited_pair_inv (int *x, int *y, const int z, const int ylim);

/******************************************************************/
/* function:  mg_binary_search */
/* searches for "target" in set between positions "begin" and "end".
 set[rank] == target if "target" is returned err_OK, "rank" contains the
 rightful position to insert "target" when returned err_NOT_FOUND.
 */
int mg_binary_search(const int target, const int *set, int begin,
                     int end, int *rank);

/******************************************************************/
/* function: mg_circumcircle */
/* builds a circulcircle from 3 coordinates listed in "coord" */
int mg_circumcircle(double coord[3][2], double center[2], double *radius);

#endif /* defined(___dmg___dmg_math__) */
