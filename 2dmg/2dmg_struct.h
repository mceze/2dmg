//
//  2dmg_struct.h
//  2dmg
//
//  Created by Marco Ceze on 11/11/14.
//  https://github.com/mceze/2dmg
//

#ifndef _dmg__dmg_struct_h
#define _dmg__dmg_struct_h

#include <gsl/gsl_interp.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas.h>

/******************************************************************/
/* List structure */
typedef struct
{
  int nItem;
  int *Item;
}
mg_List;

/******************************************************************/
/* mesh component stack */
typedef struct
{
  mg_List *Elem, *Face, *Node;
}
mg_MeshComponentStack;

/******************************************************************/
/* facedata structure*/
typedef struct
{
  int nNode; //number of nodes constituting this face
  int *node; //oriented list of nodes
  int elem[2]; //[0] -> Left; [1] -> Right
  double *normal; //store unit normal vector
  double *centroid; //store centroid
  double area;    //store area/length (3D/2D)
  double Marea;   //metric area/length (3D/2D)
}
mg_FaceData;

/******************************************************************/
/* elemdata structure */
typedef struct
{
  int nNode; //number of nodes and faces (same number)
  int *node; //oriented list of nodes
  int *face; //oriented list of faces
  int *nbor; //oriented list of neighbors
}
mg_ElemData;

/******************************************************************/
/* Ordered data list structure */
typedef struct
{
  int nEntry;      //number of entries
  int *Entry;      //ordered array of entries
  void **Data;     /* array of pointers to data corresponding entry.
                    e.g.: (Data+i)->ID == Entry[i] */
  int DataSize;
}
mg_OrderedDataList;

/******************************************************************/
/* front face structure */
struct mg_FrontFace
{
  int ID; //global number for this face
  int iloop; //loop index within front
  mg_FaceData *face; //pointer to face in mesh
  struct mg_FrontFace *next, *prev; // links to keep loop contigous
};
typedef struct mg_FrontFace mg_FrontFace;

/******************************************************************/
/* loop structure: defines a closed poligon in the mesh */
typedef struct
{
  mg_OrderedDataList *FacesInLoop;
  mg_FrontFace *head, *tail;
}
mg_Loop;

/******************************************************************/
/* front structure: single structure containing possibly more than
 one loop (front) */
typedef struct
{
  int nloop;
  mg_Loop **loop;
}
mg_Front;

/******************************************************************/
/* enumerators for boundary types */
enum mge_GeoInterp {
  mge_Linear,
  mge_Polynomial,
  mge_CSpline,
  mge_CSplinePeriodic,
  mge_Akima,
  mge_AkimaPeriodic,
  mge_GeoInterpLast
};
static char *mge_GeoInterpName[mge_GeoInterpLast] = {
  "Linear",
  "Polynomial",
  "CSpline",
  "CSplinePeriodic",
  "Akima",
  "AkimaPeriodic"
};

/******************************************************************/
/* segment structure */
typedef struct
{
  int nPoint, *Point;
  double *Coord;//first x, then y, then z
  double *s; // (0 <= s >= 1) parametric coordinate
  double Length, Lm; //physical and metric length
  //Geometry interpolation
  gsl_interp **interp;
  gsl_interp_accel **accel;
  //metric interpolation
  //interpolation object for metric values (stored as SPD)
  double *Mij;
  gsl_interp **Mij_interp;
  gsl_interp_accel **Mij_accel;
  enum mge_GeoInterp interp_type;
  char *Name;
}
mg_Segment;

/******************************************************************/
/* geometry structure */
typedef struct
{
  int nBoundary, nPoint, Dim;
  double *Coord;
  mg_Segment **Boundary;
}
mg_Geometry;

/******************************************************************/
/* list structure */
struct mg_Item {
  int Id;
  struct mg_Item *next;
};

/******************************************************************/
/* mesh structure */
typedef struct
{
  int nNode, nElem, nFace, nBfg, Dim;
  int *nBface;
  char **BNames;
  double *Coord;
  mg_ElemData *Elem;
  mg_FaceData **Face;//storing only pointers to mg_FaceData structures
  mg_List *Node2Elem, *Node2Face;
  mg_MeshComponentStack *Stack;
  
}
mg_Mesh;


#endif
