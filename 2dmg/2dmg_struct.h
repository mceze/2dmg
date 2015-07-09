//
//  2dmg_struct.h
//  2dmg
//
//  Created by Marco Ceze on 11/11/14.
//  Copyright (c) 2014 Marco Ceze. All rights reserved.
//

#ifndef _dmg__dmg_struct_h
#define _dmg__dmg_struct_h

#include "2dmg_metric.h"

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
 one loop (fronts) */
typedef struct
{
  int nloop;
  mg_Loop **loop;
}
mg_Front;

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
  mg_MetricField *Metric;
}
mg_Mesh;


#endif
