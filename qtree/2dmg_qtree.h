//
//  2dmg_qtree.h
//  2dmg
//
//  Created by Marco Ceze on 12/19/14.
//  
//

#ifndef ___dmg__qtree__
#define ___dmg__qtree__

#include <stdio.h>
#include <math.h>

//macros for children numbering
#define SOUTHWEST 0
#define SOUTHEAST 1
#define NORTHWEST 2
#define NORTHEAST 3
#define quadrant(cd,c)  1*(cd[0]>c[0])+2*(cd[1]>c[1])
#define distance(a,b) sqrt((a[0]-b[0])*(a[0]-b[0])+(a[1]-b[1])*(a[1]-b[1]))
//hardcode max capacity
#define MAXCAPACITY 2



/******************************************************************/
/* dataentry */
typedef struct
{
  double coord[2]; //coordinates of data entry
  void *data;    //pointer to data chunk
}
mg_data_entry;

/******************************************************************/
/* Tree structure */
struct mg_qtree
{
  int capacity; // max number of entry in each container
  int n_entry;  // number of entries (<= capacity)
  double c[2], ds[2]; //center and half dimension
  mg_data_entry *data; //array of "capacity" size with data and its coords.
  struct mg_qtree *child[4]; //array of children ordered lexicographically
};
typedef struct mg_qtree mg_qtree;

/******************************************************************/
/* function:  mg_init_branch */
int mg_init_branch(mg_qtree *branch);

/******************************************************************/
/* function:  mg_destroy_branch */
void mg_destroy_branch(mg_qtree *branch);

/******************************************************************/
/* function:  mg_branch_qtree */
int mg_branch_qtree(mg_qtree *qtree);

/******************************************************************/
/* function:  mg_add_qtree_entry */
int mg_add_qtree_entry(double coord[2], void **data, mg_qtree *qtree);

/******************************************************************/
/* function:  mg_find_branch */
int mg_find_branch(mg_qtree *trunk, mg_qtree **pbranch, double coord[2]);

/******************************************************************/
/* function:  mg_find_entry */
int mg_find_entry(mg_qtree *trunk, double coord[2], void *target_data,
                  double tol, mg_data_entry **data);

/******************************************************************/
/* function:  mg_plot_branch */
int mg_plot_branch(mg_qtree *branch);

/******************************************************************/
/* function:  mg_plot_qtree */
int mg_plot_qtree(mg_qtree *qtree, int strm);



#endif /* defined(___dmg__qtree__) */
