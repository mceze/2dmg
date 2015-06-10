//
//  2dmg_utils.c
//  2dmg
//
//  Created by Marco Ceze on 11/11/14.
//  Copyright (c) 2014 Marco Ceze. All rights reserved.
//

#include "2dmg_def.h"
#include "2dmg_utils.h"
#include "2dmg_struct.h"
#include "2dmg_math.h"

/******************************************************************/
/* function:  mg_alloc*/
/* wrapper for malloc with error handling*/
int mg_alloc( void **pchunk, int n, int size)
{
  int totalsize;
  
  (*pchunk) = NULL;
  totalsize = (int) n*size;

  if (totalsize == 0) return err_OK;
  if (totalsize < 0){
    printf("Warning, requesting allocation of negative memory size.\n");
    return err_OK;
  }
  if (((*pchunk) = (void *)malloc(totalsize)) == NULL)
    return error(err_MEMORY_ERROR);
  
  return err_OK;
}

/******************************************************************/
/* function:  mg_free*/
/* wrapper for free*/
void mg_free( void *chunk)
{
  if (chunk == NULL) return;
  free( (void *)chunk);
}

/******************************************************************/
/* function:  mg_alloc2*/
/* wrapper for 2d malloc with error handling*/
int mg_alloc2( void ***pchunk, int n1, int n2, int size)
{
  char *temp;
  int i, totalsize;
  
  (*pchunk) = NULL;
  totalsize = (int) n1*n2*size;
  
  if (totalsize == 0) return err_OK;
  
  if (totalsize < 0){
    printf("Warning, requesting allocation of negative memory size.\n");
    return err_OK;
  }

  if ((temp = (char *)malloc( totalsize)) == NULL)
    return error(err_MEMORY_ERROR);
  if (((*pchunk) = (void **)malloc( n1*sizeof(char *) )) == NULL)
    return error(err_MEMORY_ERROR);
  for(i = 0; i<n1; i++)
    (*pchunk)[i] = temp + (int) i*n2*size;
  
  return err_OK;
}

/******************************************************************/
/* function:  mg_free2*/
/* wrapper for 2d free*/
void mg_free2(void **chunk)
{
  if (chunk == NULL) return;
  free( (void * ) chunk[0]);
  free( (void **) chunk   );
}

/******************************************************************/
/* function:  mg_realloc*/
/* wrapper for 2d realloc*/
int mg_realloc( void **pchunk, const int n, const int size)
{
  int merr;
  int totalsize;
  
  if ((*pchunk) == NULL)
    return error(mg_alloc(pchunk, n, size));
  
  totalsize = (int) n*size;
  
  if (totalsize <= 0) {
    mg_free((*pchunk));
    (*pchunk) = NULL;
    return err_OK;
  }
  
  if (((*pchunk) = (char *)realloc( (void *)(*pchunk), totalsize)) == NULL) {
    merr = errno;
    printf("In mg_realloc: totalsize = %d, errno = %d\n", totalsize, merr);
    if (merr == ENOMEM) printf("ENOMEM\n");
    if (merr == EAGAIN) printf("EAGAIN\n");
    return error(err_MEMORY_ERROR);
  }
  
  return err_OK;
}

/******************************************************************/
/* function:  mg_init_list*/
/* initializes a mg_List*/
void mg_init_list(mg_List *list)
{
  list->nItem = 0;
  list->Item  = NULL;
}

/******************************************************************/
/* function: mg_create_mesh */
/* creates and initilizes a mesh structure */
int mg_create_mesh(mg_Mesh **pMesh)
{
  int ierr;
  call(mg_alloc((void**)pMesh, 1, sizeof(mg_Mesh)));
  (*pMesh)->Dim    = 0;
  (*pMesh)->nBfg   = 0;
  (*pMesh)->nElem  = 0;
  (*pMesh)->nFace  = 0;
  (*pMesh)->nNode  = 0;
  (*pMesh)->BNames = NULL;
  (*pMesh)->Elem = NULL;
  (*pMesh)->Face = NULL;
  (*pMesh)->Node2Elem = NULL;
  (*pMesh)->Coord = NULL;
  (*pMesh)->Node2Face = NULL;
  call(mg_alloc((void**)&(*pMesh)->Stack, 1, sizeof(mg_MeshComponentStack)));
  call(mg_alloc((void**)&(*pMesh)->Stack->Elem, 1, sizeof(mg_List)));
  mg_init_list(&(*pMesh)->Stack->Elem[0]);
  call(mg_alloc((void**)&(*pMesh)->Stack->Face, 1, sizeof(mg_List)));
  mg_init_list(&(*pMesh)->Stack->Face[0]);
  call(mg_alloc((void**)&(*pMesh)->Stack->Node, 1, sizeof(mg_List)));
  mg_init_list(&(*pMesh)->Stack->Node[0]);
  
  return err_OK;
}

/******************************************************************/
/* function: mg_destroy_mesh */
/* deallocates the memory occupied by the mesh structure */
void mg_destroy_mesh(mg_Mesh *Mesh)
{
  int i;
  
  //nBFaces
  mg_free((void*)Mesh->nBface);
  mg_free2((void**)Mesh->BNames);
  //Coord
  mg_free((void*)Mesh->Coord);
  //Elem
  for (i = 0; i < Mesh->nElem; i++) {
    mg_free((void*)Mesh->Elem[i].node);
    mg_free((void*)Mesh->Elem[i].nbor);
    mg_free((void*)Mesh->Elem[i].face);
  }
  mg_free((void*)Mesh->Elem);
  //faces
  for (i = 0; i < Mesh->nFace; i++) {
    mg_free((void*)Mesh->Face[i]->node);
    mg_free((void*)Mesh->Face[i]->normal);
    mg_free((void*)Mesh->Face[i]->centroid);
    Mesh->Face[i]->node = NULL;
    Mesh->Face[i]->normal = Mesh->Face[i]->centroid = NULL;
    mg_free((void*)Mesh->Face[i]);
  }
  mg_free((void*)Mesh->Face);
  //destroy connectivities
  for (i = 0; i < Mesh->nNode; i++) {
    if (Mesh->Node2Elem != NULL)
      mg_free((void*)Mesh->Node2Elem[i].Item);
    if (Mesh->Node2Face != NULL)
      mg_free((void*)Mesh->Node2Face[i].Item);
  }
  mg_free((void*)Mesh->Node2Elem);
  mg_free((void*)Mesh->Node2Face);
  
  //destroy stack
  mg_free((void*)Mesh->Stack->Elem->Item);
  mg_free((void*)Mesh->Stack->Elem);
  mg_free((void*)Mesh->Stack->Face->Item);
  mg_free((void*)Mesh->Stack->Face);
  mg_free((void*)Mesh->Stack->Node->Item);
  mg_free((void*)Mesh->Stack->Node);
  mg_free((void*)Mesh->Stack);
  
  
  mg_free((void*)Mesh);
}

/******************************************************************/
/* function: mg_init_face */
/* initializes a face structure */
void mg_init_face(mg_FaceData *face)
{
  face->nNode = 0;
  face->node = NULL;
  face->normal = NULL;
  face->centroid = NULL;
  face->elem[LEFTNEIGHINDEX] = -1;
  face->elem[RIGHTNEIGHINDEX] = -1;
  face->area = -1.0;
}

/******************************************************************/
/* function: mg_free_front_face */
/* frees a front face structure */
void mg_free_front_face(mg_FrontFace *FFace)
{
  if (FFace!=NULL){
    FFace->ID = -1;
    FFace->iloop = -1;
    FFace->next = FFace->prev = NULL;
    FFace->face = NULL;
    mg_free((void*)FFace);
    FFace = NULL;
  }
}

/******************************************************************/
/* function: mg_add_2_ord_set */
/* adds "entry" to a ordered "set" and keeps it ordered */
int mg_add_2_ord_set(const int entry, int *set_size, int **set,
                     int **orig_rank, bool AllowRepeat)
{
  int ierr, rank, movesize, src = 0, dest = 0;
  
  rank = 0;
  if ((*set_size) == 0)
    ierr = err_NOT_FOUND;
  else
    ierr = mg_binary_search(entry, (*set), 0, (*set_size)-1, &rank);
  
  if (ierr == err_NOT_FOUND || (ierr == err_OK && AllowRepeat)){
    ierr = error(mg_realloc((void**)&(*set), (*set_size)+1, sizeof(int)));
    if (ierr != err_OK) return ierr;
    if (orig_rank != NULL){
      ierr = error(mg_realloc((void**)&(*orig_rank), (*set_size)+1,
                                 sizeof(int)));
      if (ierr != err_OK) return ierr;
    }
    //make sure to keep the crescent order
    if (rank == (*set_size)){
      movesize = 0;
    }
    else if (rank == -1) {
      //move all one position forward
      rank = 0;
      src = 0;
      dest = 1;
      movesize = (*set_size);
    }
    else {
      src = rank;
      dest = rank+1;
      movesize = (*set_size)-rank;
    }
    if (movesize > 0) {
      if (memmove((*set)+dest,(*set)+src,movesize*sizeof(int)) == NULL)
        return error(err_MEMORY_ERROR);
      if (orig_rank != NULL)
        if (memmove((*orig_rank)+dest,(*orig_rank)+src,
                    movesize*sizeof(int)) == NULL)
          return error(err_MEMORY_ERROR);
    }
    if (orig_rank != NULL)
      orig_rank[0][rank] = (*set_size);
    set[0][rank] = entry;
    (*set_size)++;
  }
  else if (ierr != err_OK) return ierr;
  
  return err_OK;
}

/******************************************************************/
/* function: mg_rm_frm_ord_set */
/* removes entry from set */
int mg_rm_frm_ord_set(const int entry, int *set_size, int **set,
                      int n2Rm, int *nRmd)
{
  int ierr, rank = 0, dest, src, movesize;
  
  (*nRmd) = 0;
  if ((*set_size) == 0){
    (*set) = NULL;
    ierr = err_NOT_FOUND;
  }
  else
    ierr = mg_binary_search(entry, (*set), 0, (*set_size)-1, &rank);
  if (ierr == err_NOT_FOUND)
    return err_OK;
  else if (ierr == err_OK){
    //back-track from rank until set[rank] != entry;
    if ((*set)[rank] != entry) return error(err_LOGIC_ERROR);
    while (rank >=0 && (*set)[rank] == entry)
      rank--;
    rank++;
    dest = rank;
    while ((*set)[rank] == entry && n2Rm >= 0 && rank < (*set_size)) {
      rank++;
      n2Rm--;
    }
    (*nRmd) = rank-dest;
    src = rank;
    movesize = (*set_size)-src;
    if (memmove((*set)+dest,(*set)+src, movesize*sizeof(int)) == NULL)
      return error(err_MEMORY_ERROR);
    //now trim array
    (*set_size) -= (*nRmd);
    call(mg_realloc((void**)&(*set), (*set_size), sizeof(int)));
  }
  else return error(ierr);
  
  return err_OK;
}

/******************************************************************/
/* function: mg_build_connectivity */
/* builds Node2Elem and Node2Face connectivities */
int mg_build_connectivity(mg_Mesh *Mesh)
{
  int ierr, elem, node, in, face, istack;
  
  //Node2Elem
  if (Mesh->Node2Elem == NULL){
    call(mg_alloc((void**)&Mesh->Node2Elem, Mesh->nNode, sizeof(mg_List)));
    for (node = 0; node < Mesh->nNode; node++)
      mg_init_list(&Mesh->Node2Elem[node]);
  }
  istack = 0;
  for (elem = 0; elem < Mesh->nElem; elem++) {
    //skip if on stack
    if (Mesh->Stack->Elem->nItem >0)
      if (elem == Mesh->Stack->Elem->Item[istack]){
        istack++;
        continue;
      }
    for (in = 0; in < Mesh->Elem[elem].nNode; in++) {
      node = Mesh->Elem[elem].node[in];
      call(mg_add_2_ord_set(elem, &Mesh->Node2Elem[node].nItem,
                            &Mesh->Node2Elem[node].Item, NULL, false));
    }
  }
  //Node2Face
  if (Mesh->Node2Face == NULL){
    call(mg_alloc((void**)&Mesh->Node2Face, Mesh->nNode, sizeof(mg_List)));
    for (node = 0; node < Mesh->nNode; node++)
      mg_init_list(&Mesh->Node2Face[node]);
  }
  istack = 0;
  for (face = 0; face < Mesh->nFace; face++) {
    //skip if on stack
    if (Mesh->Stack->Face->nItem >0)
      if (face == Mesh->Stack->Face->Item[istack]){
        istack++;
        continue;
      }
    for (in = 0; in < Mesh->Face[face]->nNode; in++) {
      node = Mesh->Face[face]->node[in];
      call(mg_add_2_ord_set(face, &Mesh->Node2Face[node].nItem,
                            &Mesh->Node2Face[node].Item, NULL, false));
    }
  }
  
  return err_OK;
}

/******************************************************************/
/* function: mg_init_ord_data_list */
/* initializes an ordered data list*/
void mg_init_ord_data_list(mg_OrderedDataList *List, int DataSize)
{
  List->nEntry = 0;
  List->Entry  = NULL;
  List->Data  = NULL;
  List->DataSize = DataSize;
}

/******************************************************************/
/* function: mg_free_ord_data_list */
/* frees the memory ocuppied by an ordered data list*/
void mg_free_ord_data_list(mg_OrderedDataList *List)
{
  int i;
  for (i = 0; i < List->nEntry; i++)
    List->Data[i] = NULL;
  List->nEntry = 0;
  mg_free((void*)List->Entry);
  mg_free((void*)List->Data);
  List->DataSize = 0;
}

/******************************************************************/
/* function: mg_rm_frm_ord_data_list */
/* removes "data" corresponding to "entry" from "List" */
int mg_rm_frm_ord_data_list(int entry, mg_OrderedDataList *List)
{
  int ierr, rank, movesize, DataSize, src, dest;
  
  DataSize = List->DataSize;
  call(mg_binary_search(entry, List->Entry, 0, List->nEntry-1,
                        &rank));
  movesize = List->nEntry-1-rank;
  src = rank+1;
  dest = rank;
  if (movesize > 0) {
    if (memmove(List->Entry+dest,List->Entry+src,
                movesize*sizeof(int)) == NULL)
      return error(err_MEMORY_ERROR);
    if (memmove(List->Data+dest,List->Data+src,
                movesize*DataSize) == NULL)
      return error(err_MEMORY_ERROR);
  }
  List->nEntry--;
  call(mg_realloc((void**)&List->Entry, List->nEntry, sizeof(int)));
  call(mg_realloc((void**)&List->Data, List->nEntry,DataSize));
  
  return err_OK;
}

/******************************************************************/
/* function: mg_add_2_ord_data_list */
/* adds "data" corresponding to "entry" to "list". It outputs the
 data rank position in the list. It will substitute data if 
 AllowRepeat == false */
int mg_add_2_ord_data_list(int entry, void **Data,
                           mg_OrderedDataList *List, int *prank,
                           bool AllowRepeat)
{
  int ierr, rank, movesize, src = 0, dest = 0, DataSize;
  DataSize = List->DataSize;
  rank = 0;
  if (List->nEntry == 0)
    ierr = err_NOT_FOUND;
  else
    ierr = mg_binary_search(entry, List->Entry, 0, List->nEntry-1,
                            &rank);
  
  if (ierr == err_NOT_FOUND || (ierr == err_OK && AllowRepeat)){
    call(mg_realloc((void**)&List->Entry, List->nEntry+1,sizeof(int)));
    call(mg_realloc((void**)&List->Data, List->nEntry+1,DataSize));
    //make sure to keep the crescent order
    if (rank == List->nEntry){
      movesize = 0;
    }
    else if (rank == -1) {
      //move all one position forward
      rank = 0;
      src = 0;
      dest = 1;
      movesize = List->nEntry;
    }
    else {
      src = rank;
      dest = rank+1;
      movesize = List->nEntry-rank;
    }
    if (movesize > 0) {
      if (memmove(List->Entry+dest,List->Entry+src,
                  movesize*sizeof(int)) == NULL)
        return error(err_MEMORY_ERROR);
      if (memmove(List->Data+dest,List->Data+src,
                  movesize*DataSize) == NULL)
        return error(err_MEMORY_ERROR);
    }
    
    List->Data[rank] = Data[0];
    List->Entry[rank] = entry;
    List->nEntry++;
  }
  else if (ierr != err_OK) return ierr;
  
  if (prank != NULL)
    (*prank) = rank;
  
  return err_OK;
}

/******************************************************************/
/* function: mg_coord_inside_elem */
/* checks if a coordinate is inside an element */
bool mg_coord_inside_elem(mg_Mesh *Mesh, int elem, double coord[2])
{
  bool inside;
  int nleft = 0, sgn, i, d, dim = Mesh->Dim;
  double proj;
  mg_FaceData *face;
  //loop around element and if coord is to the left of all faces,
  //element contains coord.
  for (i = 0; i < Mesh->Elem[elem].nNode; i++) {
    face = Mesh->Face[Mesh->Elem[elem].face[i]];
    if (elem == face->elem[LEFTNEIGHINDEX]) sgn = 1;
    else sgn = -1;
    for (proj = 0.0, d = 0; d < dim; d++)
      proj+=sgn*face->normal[d]*(coord[d]-face->centroid[d]);
    if (proj > 1.0e-6) nleft++;
  }
  inside = (nleft == 3);
  
  return inside;
}

/******************************************************************/
/* function: mg_edges_intersect */
/* checks if edges defined by *X0 and *X1 intersect, if so,
 xint receives the intersection point*/
bool mg_edges_intersect(double X0[4], double X1[4], double *xint)
{
  double det, m11, m12, m21, m22, b0, b1, x00, x01, x10, x11;
  double y00, y01, y10, y11, qsi0, qsi1;
  
  x00 = X0[0];
  x01 = X0[2];
  x10 = X1[0];
  x11 = X1[2];
  
  y00 = X0[1];
  y01 = X0[3];
  y10 = X1[1];
  y11 = X1[3];
  
  m11 = x01-x00;
  m12 = -x11+x10;
  m21 = y01-y00;
  m22 = -y11+y10;
  
  det = m11*m22-m21*m12;
  //check if edges are approximately parallel
  if (fabs(det) <= 1e-7) return false;
  
  b0 = -x00+x10;
  b1 = -y00+y10;
  
  qsi0 = (m22*b0-m12*b1)/det;
  qsi1 = (-m21*b0+m11*b1)/det;
  
  if ((qsi0 >= 0.0 && qsi0 <= 1.0) &&
      (qsi1 >= 0.0 && qsi1 <= 1.0)){
    xint[0] = x00*(1.0-qsi0)+x01*qsi0;
    xint[1] = y00*(1.0-qsi0)+y01*qsi0;
    return true;
  }
  else
    return false;
}

/******************************************************************/
/* function: mg_find_elem_frm_coord */
/* finds element containing a point given by its coordinates  */
int mg_find_elem_frm_coord(mg_Mesh *Mesh, int elem_start,
                           double coord[2], int *pelem)
{
  int ierr, iface, i, dim, node0, node1;
  double X0[4], X1[4], Xint[2];
  bool intersect = false;
  
  dim = Mesh->Dim;
  
  X0[0] = coord[0]; X0[1] = coord[1];
  X0[2] = X0[3] = 0.0;
  //get elem centroid
  for (i = 0; i < Mesh->Elem[elem_start].nNode; i++) {
    node0 = Mesh->Elem[elem_start].node[i];
    X0[2] += Mesh->Coord[node0*dim+0];
    X0[3] += Mesh->Coord[node0*dim+1];
  }
  X0[2] /= Mesh->Elem[elem_start].nNode;
  X0[3] /= Mesh->Elem[elem_start].nNode;
//  printf("%d %1.6e %1.6e %1.6e %1.6e\n",elem_start,X0[0],X0[1],X0[2],X0[3]);
  
  for (iface = 0; iface < Mesh->Elem[elem_start].nNode; iface++) {
    node0 = Mesh->Face[Mesh->Elem[elem_start].face[iface]]->node[0];
    node1 = Mesh->Face[Mesh->Elem[elem_start].face[iface]]->node[1];
    X1[0] = Mesh->Coord[node0*dim+0];
    X1[1] = Mesh->Coord[node0*dim+1];
    X1[2] = Mesh->Coord[node1*dim+0];
    X1[3] = Mesh->Coord[node1*dim+1];
    if ((intersect = mg_edges_intersect(X0, X1, Xint)) == true)
      break;
  }
  //if edges intersect, change elem_start to its neighbor accross iface
  //and continue recurrence
  if (intersect) {
    if (Mesh->Elem[elem_start].nbor[iface] < 0) {
      (*pelem) = elem_start;
      return err_NOT_FOUND;
    }
    elem_start = Mesh->Elem[elem_start].nbor[iface];
    call(mg_find_elem_frm_coord(Mesh, elem_start, coord, pelem));
  }
  else {
    //found the element but let's verify
    if (!mg_coord_inside_elem(Mesh, elem_start, coord))
      return error(err_LOGIC_ERROR);
    (*pelem) = elem_start;
  }
  return err_OK;
}

/******************************************************************/
/* function: mg_check_exist */
/* checks if "num" exists in "vector" and returns "index" */
void mg_check_exist(int num, int size, int *vector, int *index)
{
  int i;
  (*index) = -1;
  
  for (i = 0; i < size; i++){
    if (num == vector[i]){
      (*index) = i;
      break;
    }
  }
}

/******************************************************************/
/* function: mg_value_2_enum */
/* matches a enumerator with its name */
int
mg_value_2_enum(const char value[], char *EName[], int ELast,
                int *ival)
{
  int k;
  for (k=0; k<ELast; k++){
    if (strcmp(EName[k], value) == 0){
      *ival = k;
      return err_OK;
    }
  }
  
  return err_NOT_FOUND;
}

