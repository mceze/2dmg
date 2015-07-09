//
//  2dmg.c
//  2dmg
//
//  Created by Marco Ceze on 11/11/14.
//  https://github.com/mceze/2dmg
//

#include <stdlib.h>
#include "2dmg_def.h"
#include "2dmg_struct.h"
#include "2dmg_utils.h"
#include "2dmg_math.h"
#include "2dmg_geo.h"
#include "2dmg_plot.h"

/******************************************************************/
/* function: mg_calc_face_info */
/* calculates face properties */
int mg_calc_face_info(mg_Mesh *Mesh)
{
  int ierr, f, i, *node, in, istack;
  double delta[3];
  mg_FaceData *face;
  
  istack = 0;
  for (f = 0; f < Mesh->nFace+Mesh->Stack->Face->nItem; f++) {
    face = Mesh->Face[f];
    //if face normal is null, assume all face data is stale or uninitialized
    //do not calculate if on stack
    if (Mesh->Stack->Face->nItem >0)
      if (f == Mesh->Stack->Face->Item[istack]){
        istack++;
        continue;
      }
    if (face->normal == NULL) {
      node = face->node;
      call(mg_alloc((void**)&face->normal, Mesh->Dim, sizeof(double)));
      call(mg_alloc((void**)&face->centroid, Mesh->Dim, sizeof(double)));
      switch (Mesh->Dim) {
        case 2:
          face->area = 0.0;
          for (i=0;i<Mesh->Dim;i++){
            delta[i] = Mesh->Coord[node[1]*Mesh->Dim+i]-Mesh->Coord[node[0]*Mesh->Dim+i];
            face->area += delta[i]*delta[i];
            face->centroid[i] = 0.0;
            for (in = 0; in < face->nNode; in++)
              face->centroid[i] += Mesh->Coord[node[in]*Mesh->Dim+i]/face->nNode;
          }
          face->area = sqrt(face->area);
          face->normal[0] = -delta[1]/face->area;
          face->normal[1] = delta[0]/face->area;
          break;
        case 3:
          return error(err_NOT_SUPPORTED);
          break;
        default:
          break;
      }
    }
  }
  
  return err_OK;
}

/******************************************************************/
/* function: mg_free_loop */
/* frees the memory stored in loop */
int mg_free_loop(mg_Loop *Loop)
{
  mg_FrontFace *FFace, *FFaceNext;
  FFace = Loop->head;
  while (FFace != Loop->tail) {
    FFaceNext = FFace->next;
    mg_free_front_face(FFace);
    FFace = FFaceNext;
  }
  mg_free_front_face(FFace);//destroy the tail also
  Loop->head = NULL;
  Loop->tail = NULL;
  mg_free_ord_data_list(Loop->FacesInLoop);
  
  return err_OK;
}

/******************************************************************/
/* function: mg_build_loop */
/* builds a closed loop of front faces oriented counter-clockwise*/
int mg_build_loop(mg_Mesh *Mesh, int faceID, mg_Loop *Loop, int iloop)
{
  int ierr, f, nborface = -1, ID, node1;
  bool open, foundnext;
  mg_FrontFace *FFace, *NextFFace;
  mg_OrderedDataList *FacesInLoop;
  
  //initialize list of faces in loop
  //(this facilitates finding if a face is in a loop)
  call(mg_alloc((void**)&FacesInLoop, 1, sizeof(mg_OrderedDataList)));
  Loop->FacesInLoop = FacesInLoop;
  mg_init_ord_data_list(Loop->FacesInLoop, sizeof(mg_FrontFace));
  //create loop's head
  call(mg_alloc((void**)&FFace, 1, sizeof(mg_FrontFace)));
  FFace->ID = faceID;
  FFace->iloop = iloop;
  //store pointer to mesh face
  FFace->face = Mesh->Face[faceID];
  FFace->prev = FFace->next = NULL;
  //start at head and follow path until it closes
  Loop->head = FFace;
  Loop->tail = NULL;
  open = true;
  while (open) {
    ID = FFace->ID;
    //add face to list
    call(mg_add_2_ord_data_list(ID, (void**)&FFace,
                                Loop->FacesInLoop, NULL, false));
    foundnext = false;
    node1 = FFace->face->node[1];
    for (f = 0; f < Mesh->Node2Face[node1].nItem; f++) {
      nborface = Mesh->Node2Face[node1].Item[f];
      if (nborface == ID) continue;
      //find neighboring face for which the unmeshed domain is to its left
      if (Mesh->Face[nborface]->elem[LEFTNEIGHINDEX] == HOLLOWNEIGHTAG){
        foundnext = true;
        break;
      }
    }
    //should have found the next front face
    if (!foundnext)
      return error(err_LOGIC_ERROR);
    call(mg_alloc((void**)&NextFFace, 1, sizeof(mg_FrontFace)));
    if (nborface == -1)
      return error(err_LOGIC_ERROR);
    NextFFace->ID = nborface;
    NextFFace->iloop = iloop;
    //store pointer to mesh face
    NextFFace->face = Mesh->Face[nborface];
    NextFFace->prev = FFace;
    //check if loop is closed
    if (NextFFace->ID == Loop->head->ID){
      open = false;
      Loop->tail = NextFFace->prev;
    }
    else{
      NextFFace->next = NULL;
      //link with previous face
      FFace->next = NextFFace;
      //point to next face
      FFace = NextFFace;
    }
  }
  Loop->head->prev = Loop->tail;
  Loop->tail->next = Loop->head;
  NextFFace = FFace = NULL;
  
  return err_OK;
}

/******************************************************************/
/* function: mg_create_front */
/* creates a front - a collection of loops of front faces */
int mg_create_front (mg_Mesh *Mesh, mg_Front *Front)
{
  int ierr, faceID, iloop, istack;
  bool foundseed;
  mg_Loop *Loop;
  
  //initialize front
  Front->nloop = 0;
  Front->loop = NULL;
  //loop over faces, pick a seed face and generate loop.
  //stop when can't find any more seeds
  istack = 0;
  for (faceID = 0; faceID < Mesh->nFace; faceID++) {
    //do not consider if on stack
    if (Mesh->Stack->Face->nItem >0)
      if (faceID == Mesh->Stack->Face->Item[istack]){
        istack++;
        continue;
      }
    if (Mesh->Face[faceID]->elem[LEFTNEIGHINDEX] == HOLLOWNEIGHTAG){
      //check if it is the first loop
      if (Front->nloop == 0){
        call(mg_alloc((void**)&Loop, 1, sizeof(mg_Loop)));
        call(mg_build_loop(Mesh, faceID, Loop, 0));
      }
      else {
        foundseed = true;
        //loop over loops check if faceID is part of any other loop
        for (iloop = 0; iloop < Front->nloop; iloop++) {
          if (Front->loop[iloop]->FacesInLoop->nEntry == 0) continue;
          ierr = mg_binary_search(faceID, Front->loop[iloop]->FacesInLoop->Entry,0,
                                  Front->loop[iloop]->FacesInLoop->nEntry-1, NULL);
          if (ierr == err_OK) foundseed = false;
          else if (ierr != err_NOT_FOUND) return error(ierr);
        }
        if (foundseed){
          call(mg_alloc((void**)&Loop, 1, sizeof(mg_Loop)));
          call(mg_build_loop(Mesh, faceID, Loop, Front->nloop));
        }
        else continue;
      }
      Front->nloop++;
      call(mg_realloc((void**)&Front->loop,Front->nloop, sizeof(mg_Loop)));
      Front->loop[Front->nloop-1] = Loop;
    }
  }
  
  return err_OK;
}

/******************************************************************/
//   FUNCTION Definition: xf_VerifyFront2D
/* NOTE: this function is for verifying the assumptions of the front
 that ensure proper propagation of the front into the domain. It is
 supposed to be used only for debugging */
int xf_VerifyFront2D(mg_Front *Front)
{
  int ierr, loopID, *nfacevisit;
  int iface, rank;
  mg_Loop *Loop;
  mg_FrontFace *FFace, *FFtest;
  
  for (loopID = 0; loopID < Front->nloop; loopID++) {
    Loop = Front->loop[loopID];
    if (Loop->FacesInLoop->nEntry == 0) {
      //printf("Loop %d is empty.\n",loopID);
      continue;
    }
    /*go around loop and verify if it closes and if all (and only those)
     faces are present in the list*/
    call(mg_alloc((void**)&nfacevisit,Loop->FacesInLoop->nEntry,
                  sizeof(int)));
    for (iface = 0; iface < Loop->FacesInLoop->nEntry; iface++)
      nfacevisit[iface] = 0;
    FFace = Loop->head;
    iface = 0;
    call(mg_binary_search(FFace->ID,Loop->FacesInLoop->Entry,0,
                          Loop->FacesInLoop->nEntry-1, &rank));
    FFtest = (mg_FrontFace*)Loop->FacesInLoop->Data[rank];
    if (FFtest != FFace){
      printf("Inconsistent FrontFace pointers in face list for face %d.\n",
             FFace->ID);
      return error(err_LOGIC_ERROR);
    }
    FFace = FFace->next;
    nfacevisit[rank]++;
    iface++;
    while (FFace != Loop->head && iface < Loop->FacesInLoop->nEntry) {
      call(mg_binary_search(FFace->ID,Loop->FacesInLoop->Entry,0,
                            Loop->FacesInLoop->nEntry-1, &rank));
      FFtest = (mg_FrontFace*)Loop->FacesInLoop->Data[rank];
      if (FFtest != FFace){
        printf("Inconsistent FrontFace pointers in face list for face %d.\n",
               FFace->ID);
        return error(err_LOGIC_ERROR);
      }
      nfacevisit[rank]++;
      if (nfacevisit[rank] != 1) {
        printf("Face %d occurs more than once in loop.\n", FFace->ID);
        return error(err_LOGIC_ERROR);
      }
      if (FFace->face->elem[LEFTNEIGHINDEX] != HOLLOWNEIGHTAG) {
        printf("Face %d left side is not the hollow region.\n", FFace->ID);
        return error(err_LOGIC_ERROR);
      }
      FFace = FFace->next;
      iface++;
    }
    mg_free((void*)nfacevisit);
    if (iface != Loop->FacesInLoop->nEntry){
      printf("Number of faces in Loop %d is inconsistent with its list.\n",
             loopID);
      return error(err_LOGIC_ERROR);
    }
  }
  
  return err_OK;
}

/******************************************************************/
/* function: mg_find_seed_face */
/* selects a face to advance from front */
int mg_find_seed_face(mg_Front *Front, mg_FrontFace **pSeedFace)
{
  int iloop;
  bool FoundSeed = false;
  mg_FrontFace *FFace;
  
  if (Front->nloop < 1) return error(err_INPUT_ERROR);
  
  //for now, hard code minimum length criterion
  for (iloop = 0; iloop < Front->nloop; iloop++) {
    FFace = Front->loop[iloop]->head;
    if (FFace == NULL) continue;//empty loop
    if (!FoundSeed){
      (*pSeedFace) = FFace;
      FoundSeed = true;
    }
    while (FFace != Front->loop[iloop]->tail) {
      if (FFace->face->area < (*pSeedFace)->face->area) {
        (*pSeedFace) = FFace;
      }
      FFace = FFace->next;
    }
  }
  
  return err_OK;
}

/******************************************************************/
/* function: mg_find_node_in_front */
/* finsd node with "nodeID0" in the front */
int mg_find_node_in_front(int nodeID0, mg_Mesh *Mesh, mg_Front *Front,
                          mg_OrderedDataList **pNode2FFace)
{
  int ierr, iface, faceID, iloop, rank;
  mg_Loop *Loop;
  mg_OrderedDataList *Node2FFace;
  mg_FrontFace *FFace;
  
  call(mg_alloc((void**)&Node2FFace, 1, sizeof(mg_OrderedDataList)));
  mg_init_ord_data_list(Node2FFace, sizeof(mg_FrontFace));
  
  for (iface = 0; iface < Mesh->Node2Face[nodeID0].nItem; iface++) {
    faceID = Mesh->Node2Face[nodeID0].Item[iface];
    for (iloop = 0; iloop < Front->nloop; iloop++) {
      Loop = Front->loop[iloop];
      if (Loop->FacesInLoop->nEntry == 0) continue;
      ierr = mg_binary_search(faceID, Loop->FacesInLoop->Entry, 0,
                              Loop->FacesInLoop->nEntry-1, &rank);
      if (ierr == err_OK){
        FFace = (mg_FrontFace*)(Loop->FacesInLoop->Data[rank]);
        call(mg_add_2_ord_data_list(faceID, (void**)&FFace, Node2FFace,
                                    NULL, true));
      }
      else if (ierr != err_NOT_FOUND) return error(ierr);
    }
  }
  if (Node2FFace->nEntry == 0)
    mg_free_ord_data_list(Node2FFace);
  (*pNode2FFace) = Node2FFace;
  
  return err_OK;
}

/******************************************************************/
/* function: mg_comp_mesh_size*/
/* computes the mesh size for face FFace*/
int mg_comp_mesh_size(mg_Mesh *Mesh, mg_FrontFace *FFace, double **prho,
                      bool *isotropic)
{
  int ierr, nface, inode, iface, nodeID, faceID;
  
  //for now, make it isotropic and solely based on edge length
  (*isotropic) = true;
  call(mg_alloc((void**)prho, 1, sizeof(double)));
  //compute rho as the average length of all faces connected to FFace
  //start at -area/-1 so it does not double count FFace
  (*prho)[0] = -FFace->face->area;
  nface = -1;
  for (inode = 0; inode < FFace->face->nNode; inode++) {
    nodeID = FFace->face->node[inode];
    for (iface = 0; iface < Mesh->Node2Face[nodeID].nItem; iface++) {
      faceID = Mesh->Node2Face[nodeID].Item[iface];
      (*prho)[0] += Mesh->Face[faceID]->area;
      nface++;
    }
  }
  (*prho)[0] /= nface;
  
  return err_OK;
}

/******************************************************************/
/* function: mg_nodes_frnt_dist */
/* selects nodes in front within a distance "rho" */
int mg_nodes_frnt_dist(mg_Mesh *Mesh, mg_Front *Front,
                       mg_FrontFace *SelfFace, double *rho, double c,
                       bool IsoFlag, mg_List *NodeList)
{
  int ierr, iloop, nodeID0, inode0, inode1, nodeID1, d, dim;
  double distance, delta;
  mg_FrontFace *FFace;
  
  if (IsoFlag == false) return error(err_NOT_SUPPORTED);
  
  NodeList->nItem = 0;
  NodeList->Item = NULL;
  dim = Mesh->Dim;
  //loop over front and compute distances
  for (inode0 = 0; inode0 < SelfFace->face->nNode; inode0++){
    nodeID0 = SelfFace->face->node[inode0];
    for (iloop = 0; iloop < Front->nloop; iloop++) {
      if (Front->loop+iloop == NULL) continue;
      FFace = Front->loop[iloop]->head;
      while (FFace != Front->loop[iloop]->tail) {
        if (FFace != SelfFace) {//skip self face
          for (inode1 = 0; inode1 < FFace->face->nNode; inode1++) {
            nodeID1 = FFace->face->node[inode1];
            distance = 0.0;
            for (d = 0; d < dim; d++) {
              delta = Mesh->Coord[nodeID1*dim+d]-Mesh->Coord[nodeID0*dim+d];
              distance += delta*delta;
            }
            //check if within distance
            if (distance <= c*(*rho)*c*(*rho))
              call(mg_add_2_ord_set(nodeID1, &NodeList->nItem,
                                    &NodeList->Item, NULL, false));
          }
        }
        FFace = FFace->next;
      }
    }
  }
  
  return err_OK;
}

/******************************************************************/
/* function: mg_build_circle_frm_face */
/* builds a circumcircle from a face and a node and outputs
 its radius*/
int mg_build_circle_frm_face(mg_Mesh *Mesh, mg_FaceData *face,
                             double *NCoord, double center[2],
                             double *radius)
{
  int ierr, dim = Mesh->Dim;
  double coord[3][2];
  
  coord[0][0] = Mesh->Coord[face->node[0]*dim];
  coord[0][1] = Mesh->Coord[face->node[0]*dim+1];
  coord[1][0] = Mesh->Coord[face->node[1]*dim];
  coord[1][1] = Mesh->Coord[face->node[1]*dim+1];
  coord[2][0] = NCoord[0];
  coord[2][1] = NCoord[1];
  call(mg_circumcircle(coord, center, radius));
  
  return err_OK;
}

/******************************************************************/
/* function: mg_tri_frm_face_node */
/* builds a triangle from a face and a node and adds it to the
 Mesh*/
int mg_tri_frm_face_node(mg_Mesh *Mesh, mg_FrontFace *FFace, int nodeID2,
                         double *coord)
{
  int ierr, elemID0, in, nodeID, faceID0, faceID1, faceID2, f, faceID, i, t, idx;
  int nbor;
  bool face0new, face1new, newnode, NodeIsFromStack, Face0IsFromStack;
  bool Face1IsFromStack;
  mg_FaceData *face;
  
  /* Steps:
   1- add element to mesh
   2- add node to mesh if node is new
   3- update Node2Elem
   4- add faces to mesh
   5- update Node2Face
   
   |       2 (this is nodeID2 received by this function)
   |       *
   |    1 / \ 0
   |     / e \
   |   0*-----*1
   |       2 (front face)
   
   6- Front gets update outside from FFace */
  
  /*************************************************************************/
  //element ID: add element to mesh, get id from element stack if not empty,
  //            otherwise, add element to mesh.
  /*************************************************************************/
  Mesh->nElem++;
  if (Mesh->Stack->Elem->nItem == 0) {
    //empty stack, need to allocate space
    elemID0 = Mesh->nElem-1;
    call(mg_realloc((void**)&Mesh->Elem, Mesh->nElem, sizeof(mg_ElemData)));
    Mesh->Elem[elemID0].nNode = FFace->face->nNode+1;
    call(mg_alloc((void**)&Mesh->Elem[elemID0].node, Mesh->Elem[elemID0].nNode,
                  sizeof(int)));
    //face info (note: nnode == nface == number of neighbors)
    call(mg_alloc((void**)&Mesh->Elem[elemID0].face, Mesh->Elem[elemID0].nNode,
                  sizeof(int)));
    call(mg_alloc((void**)&Mesh->Elem[elemID0].nbor, Mesh->Elem[elemID0].nNode,
                  sizeof(int)));
  }
  else {
    //get next element number in stack and update stack
    elemID0 = Mesh->Stack->Elem->Item[0];
    call(mg_rm_frm_ord_set(elemID0, &(Mesh->Stack->Elem->nItem),
                           &(Mesh->Stack->Elem->Item), 1, &t));
    if (t != 1) return error(err_MESH_ERROR);
  }
  //elem2node info
  for (in = 0; in < FFace->face->nNode; in++){
    nodeID = FFace->face->node[in];
    Mesh->Elem[elemID0].node[in] = nodeID;
  }
  Mesh->Elem[elemID0].node[in] = nodeID2;
  
  Mesh->Elem[elemID0].face[2] = faceID2 = FFace->ID;
  Mesh->Face[faceID2]->elem[LEFTNEIGHINDEX] = elemID0;
  
  /*************************************************************************/
  //Faces and node: check if node came from stack, if so, node is
  //                considered new and we need to use new coordinates.
  //                if node did not come from stack but it is new,
  //                allocate space for it.
  /*************************************************************************/
  face0new = face1new = false;
  call(mg_rm_frm_ord_set(nodeID2, &(Mesh->Stack->Node->nItem),
                         &(Mesh->Stack->Node->Item), 1, &t));
  NodeIsFromStack = (t==1);
  newnode = (NodeIsFromStack || nodeID2 >= Mesh->nNode);
  //Now, take care of the node and faces
  if (!newnode) {
    //existing node
    //should not pass in coord if node is not new
    if (coord != NULL) return error(err_INPUT_ERROR);
    //check if either face0 or face1 are connected to nodeID2
    faceID0 = faceID1 = -1;
    for (f = 0; f < Mesh->Node2Face[nodeID2].nItem; f++) {
      faceID = Mesh->Node2Face[nodeID2].Item[f];
      for (i = 0; i < 2; i++){
        if (nodeID2 == Mesh->Face[faceID]->node[i]){
          nodeID = Mesh->Face[faceID]->node[!i];//get other node
          if (nodeID == FFace->face->node[0])
            faceID1 = faceID;
          if (nodeID == FFace->face->node[1])
            faceID0 = faceID;
        }
      }
    }
    Face0IsFromStack = false;
    if (faceID0 == -1) {
      //faceID0 is not connected to nodeID2, this is a new face,
      //so check stack to see if we need to allocate space for it
      if (Mesh->Stack->Face->nItem == 0)
        faceID0 = Mesh->nFace;
      else {//stack is not empty
        faceID0 = Mesh->Stack->Face->Item[0];
        call(mg_rm_frm_ord_set(faceID0, &(Mesh->Stack->Face->nItem),
                               &(Mesh->Stack->Face->Item), 1, &t));
        if (t != 1) return error(err_MESH_ERROR);
        Face0IsFromStack = true;
      }
      Mesh->nFace++;
      face0new = true;
    }
    Face1IsFromStack = false;
    if (faceID1 == -1) {
      //faceID1 is not connected to nodeID2, this is a new face,
      //so check stack to see if we need to allocate space for it
      if (Mesh->Stack->Face->nItem == 0)
        faceID1 = Mesh->nFace;
      else {//stack is not empty
        faceID1 = Mesh->Stack->Face->Item[0];
        call(mg_rm_frm_ord_set(faceID1, &(Mesh->Stack->Face->nItem),
                               &(Mesh->Stack->Face->Item), 1, &t));
        if (t != 1) return error(err_MESH_ERROR);
        Face1IsFromStack = true;
      }
      Mesh->nFace++;
      face1new = true;
    }
  }
  else {
    //new node
    //add node coordinates
    if (coord == NULL) return error(err_INPUT_ERROR);
    //new node can only come from stack or be sequentially added
    if (!NodeIsFromStack && nodeID2 != Mesh->nNode)
      return error(err_INPUT_ERROR);
    
    Mesh->nNode++;
    if (NodeIsFromStack == false) {
      //create space for coordinatess
      //if node was not on stack, create space for it
      call(mg_realloc((void**)&Mesh->Coord, Mesh->nNode*Mesh->Dim, sizeof(double)));
      //create space for Node2Elem and Node2Face connectivities
      call(mg_realloc((void**)&Mesh->Node2Elem, Mesh->nNode, sizeof(mg_List)));
      call(mg_realloc((void**)&Mesh->Node2Face, Mesh->nNode, sizeof(mg_List)));
      //initialize connectivities
      mg_init_list(&Mesh->Node2Elem[nodeID2]);
      mg_init_list(&Mesh->Node2Face[nodeID2]);
    }
    for (i = 0; i < Mesh->Dim; i++)
      Mesh->Coord[nodeID2*Mesh->Dim+i] = coord[i];
    call(mg_add_2_ord_set(elemID0, &Mesh->Node2Elem[nodeID2].nItem,
                          &Mesh->Node2Elem[nodeID2].Item, NULL, false));
    //new faces, check stack for both
    Face0IsFromStack = false;
    if (Mesh->Stack->Face->nItem == 0)
      faceID0 = Mesh->nFace;
    else {//stack is not empty
      faceID0 = Mesh->Stack->Face->Item[0];
      call(mg_rm_frm_ord_set(faceID0, &(Mesh->Stack->Face->nItem),
                             &(Mesh->Stack->Face->Item), 1, &t));
      if (t != 1) return error(err_MESH_ERROR);
      Face0IsFromStack = true;
    }
    Mesh->nFace++;
    face0new = true;
    
    Face1IsFromStack = false;
    if (Mesh->Stack->Face->nItem == 0)
      faceID1 = Mesh->nFace;
    else {//stack is not empty
      faceID1 = Mesh->Stack->Face->Item[0];
      call(mg_rm_frm_ord_set(faceID1, &(Mesh->Stack->Face->nItem),
                             &(Mesh->Stack->Face->Item), 1, &t));
      if (t != 1) return error(err_MESH_ERROR);
      Face1IsFromStack = true;
    }
    Mesh->nFace++;
    face1new = true;
    
    call(mg_add_2_ord_set(faceID0, &Mesh->Node2Face[nodeID2].nItem,
                          &Mesh->Node2Face[nodeID2].Item, NULL, false));
    call(mg_add_2_ord_set(faceID1, &Mesh->Node2Face[nodeID2].nItem,
                          &Mesh->Node2Face[nodeID2].Item, NULL, false));
  }
  //face IDs (faceID2 is assigned above)
  Mesh->Elem[elemID0].face[1] = faceID1;
  Mesh->Elem[elemID0].face[0] = faceID0;
  //now, take care of new faces
  if (!Face0IsFromStack || !Face1IsFromStack)//if either face is not from stack, allocate space
    call(mg_realloc((void**)&Mesh->Face, Mesh->nFace, sizeof(mg_FaceData)));
  if (face0new){
    if (!Face0IsFromStack){
      //create face structure and store pointer
      call(mg_alloc((void**)&face, 1, sizeof(mg_FaceData)));
      Mesh->Face[faceID0] = face;
      call(mg_alloc((void**)&Mesh->Face[faceID0]->node, 2, sizeof(int)));
    }
    Mesh->Face[faceID0]->nNode = 2;
    Mesh->Face[faceID0]->node[0] = nodeID2;
    Mesh->Face[faceID0]->node[1] = FFace->face->node[1];
    Mesh->Face[faceID0]->elem[LEFTNEIGHINDEX] = HOLLOWNEIGHTAG;
    Mesh->Face[faceID0]->elem[RIGHTNEIGHINDEX] = elemID0;
    Mesh->Elem[elemID0].nbor[0] = Mesh->Face[faceID0]->elem[LEFTNEIGHINDEX];
    //face info will be updated later
    //note: if face came from the stack, these will be NULL already
    Mesh->Face[faceID0]->normal = NULL;
    Mesh->Face[faceID0]->centroid = NULL;
    Mesh->Face[faceID0]->area = -1.0;
    for (i = 0; i < 2; i++) {
      nodeID = Mesh->Face[faceID0]->node[i];
      call(mg_add_2_ord_set(faceID0, &Mesh->Node2Face[nodeID].nItem,
                            &Mesh->Node2Face[nodeID].Item, NULL, false));
    }
  }
  else {
    //update elemental connectivity
    Mesh->Face[faceID0]->elem[LEFTNEIGHINDEX] = elemID0;
    Mesh->Elem[elemID0].nbor[0] = Mesh->Face[faceID0]->elem[RIGHTNEIGHINDEX];
    nbor = Mesh->Face[faceID0]->elem[RIGHTNEIGHINDEX];
    if (nbor >= 0) { //not a boundary
      mg_check_exist(faceID0, Mesh->Elem[nbor].nNode, Mesh->Elem[nbor].face, &idx);
      Mesh->Elem[nbor].nbor[idx] = elemID0;
    }
  }
  if (face1new){
    if (!Face1IsFromStack){
      //create face structure and store pointer
      call(mg_alloc((void**)&face, 1, sizeof(mg_FaceData)));
      Mesh->Face[faceID1] = face;
      call(mg_alloc((void**)&Mesh->Face[faceID1]->node, 2, sizeof(int)));
    }
    Mesh->Face[faceID1]->nNode = 2;
    Mesh->Face[faceID1]->node[1] = nodeID2;
    Mesh->Face[faceID1]->node[0] = FFace->face->node[0];
    Mesh->Face[faceID1]->elem[LEFTNEIGHINDEX] = HOLLOWNEIGHTAG;
    Mesh->Face[faceID1]->elem[RIGHTNEIGHINDEX] = elemID0;
    Mesh->Elem[elemID0].nbor[1] = Mesh->Face[faceID1]->elem[LEFTNEIGHINDEX];
    //face info will be updated later
    Mesh->Face[faceID1]->normal = NULL;
    Mesh->Face[faceID1]->centroid = NULL;
    Mesh->Face[faceID1]->area = -1.0;
    for (i = 0; i < 2; i++) {
      nodeID = Mesh->Face[faceID1]->node[i];
      call(mg_add_2_ord_set(faceID1, &Mesh->Node2Face[nodeID].nItem,
                            &Mesh->Node2Face[nodeID].Item, NULL, false));
    }
  }
  else {
    //update elemental connectivity
    Mesh->Face[faceID1]->elem[LEFTNEIGHINDEX] = elemID0;
    Mesh->Elem[elemID0].nbor[1] = Mesh->Face[faceID1]->elem[RIGHTNEIGHINDEX];
    nbor = Mesh->Face[faceID1]->elem[RIGHTNEIGHINDEX];
    if (nbor >= 0) { //not a boundary
      mg_check_exist(faceID1, Mesh->Elem[nbor].nNode, Mesh->Elem[nbor].face, &i);
      if (i < 0) return error(err_LOGIC_ERROR);
      Mesh->Elem[nbor].nbor[i] = elemID0;
    }
  }
  //calculate normals,areas and centroids for new faces
  //only calculates the uninitialized values (new faces)
  call(mg_calc_face_info(Mesh));
  
  Mesh->Face[faceID2]->elem[LEFTNEIGHINDEX] = elemID0;
  Mesh->Elem[elemID0].nbor[2] = Mesh->Face[faceID2]->elem[RIGHTNEIGHINDEX];
  if ((nbor = Mesh->Face[faceID2]->elem[RIGHTNEIGHINDEX]) >= 0) { //not a boundary
    mg_check_exist(faceID2, Mesh->Elem[nbor].nNode, Mesh->Elem[nbor].face, &i);
    if (i < 0)
      return error(err_LOGIC_ERROR);
    Mesh->Elem[nbor].nbor[i] = elemID0;
  }
  
  return err_OK;
}

/******************************************************************/
/* function: mg_update_front */
/* updates front */
int mg_update_front(mg_Mesh *Mesh, mg_Front *Front, mg_FrontFace *FF2)
{
  int ierr, elemID0, faceID0, faceID1, nodeID2, f, faceID = 0, faceID2;
  int iloop, oldloop, newloop, leftloop, rightloop;
  bool face0new, face1new, node2new, sameloop;
  mg_FrontFace *FF1=NULL, *FF0=NULL, *FFace=NULL, *FF0Next=NULL;
  mg_FrontFace *FF2Prev=NULL, *FF2Next=NULL, *FF1Next=NULL;
  mg_FrontFace *FF0Prev=NULL, *FF1Prev=NULL;
  mg_Loop *Loop, *LoopNew;
  mg_OrderedDataList *Node2FFace;
  
  /* FFace should not be a valid front face at this point, i.e.,
   left side of face should not point to HOLLOWNEIGHTAG
   */
  elemID0 = FF2->face->elem[LEFTNEIGHINDEX];
  faceID2 = FF2->ID;
  
  if (elemID0 == HOLLOWNEIGHTAG) return error(err_INPUT_ERROR);
  //FFace should be face 2 of elem
  if (Mesh->Elem[elemID0].face[2] != FF2->ID)
    return error(err_LOGIC_ERROR);
  faceID0 = Mesh->Elem[elemID0].face[0];
  faceID1 = Mesh->Elem[elemID0].face[1];
  face0new = face1new = false;
  if (Mesh->Face[faceID0]->elem[LEFTNEIGHINDEX] == HOLLOWNEIGHTAG)
    face0new = true;
  if (Mesh->Face[faceID1]->elem[LEFTNEIGHINDEX] == HOLLOWNEIGHTAG)
    face1new = true;
  //check if node2 has just been created
  nodeID2 = Mesh->Elem[elemID0].node[2];
  node2new = false;
  if (Mesh->Node2Face[nodeID2].nItem == 2)
    node2new = true;
  
  if (node2new) {
    //both faces must be new also, let's check
    if (!face0new || !face1new) return error(err_LOGIC_ERROR);
    //create a FrontFace structure for faceID0 and faceID1
    call(mg_alloc((void**)&FF0, 1, sizeof(mg_FrontFace)));
    call(mg_alloc((void**)&FF1, 1, sizeof(mg_FrontFace)));
    FF2Prev = FF2->prev;
    FF2Next = FF2->next;
    //setup face0
    FF0->ID = faceID0;
    FF0->iloop = FF2->iloop;
    FF0->face = Mesh->Face[faceID0];
    FF0->next = FF2->next;
    FF0->prev = FF1;
    //setup face1
    FF1->ID = faceID1;
    FF1->iloop = FF2->iloop;
    FF1->face = Mesh->Face[faceID1];
    FF1->next = FF0;
    FF1->prev = FF2->prev;
    //stitch loop
    FF2Prev->next = FF1;
    FF2Next->prev = FF0;
    Loop = Front->loop[FF2->iloop];
    call(mg_rm_frm_ord_data_list(faceID2, Loop->FacesInLoop));
    call(mg_add_2_ord_data_list(FF0->ID, (void**)&FF0, Loop->FacesInLoop,
                                NULL, false));
    call(mg_add_2_ord_data_list(FF1->ID, (void**)&FF1, Loop->FacesInLoop,
                                NULL, false));
    Loop->head = FF0;
    Loop->tail = FF1;
    mg_free_front_face(FF2);
  }
  else {
    //one new face
    if ((face0new && !face1new) || (face1new && !face0new)){
      Loop = Front->loop[FF2->iloop];
      if (face0new){// face 0 is new
        //create a FrontFace structure for faceID0
        call(mg_alloc((void**)&FF0, 1, sizeof(mg_FrontFace)));
        FF1 = FF2->prev;
        FF2Next = FF2->next;
        FF1Prev = FF1->prev;
        //setup face0
        FF0->ID = faceID0;
        FF0->iloop = FF2->iloop;
        FF0->face = Mesh->Face[faceID0];
        FF0->next = FF2->next;
        FF0->prev = FF1->prev;
        //stitch loop
        FF2Next->prev = FF0;
        FF1Prev->next = FF0;
        //remove faces from loop
        call(mg_rm_frm_ord_data_list(faceID1, Loop->FacesInLoop));
        call(mg_rm_frm_ord_data_list(faceID2, Loop->FacesInLoop));
        call(mg_add_2_ord_data_list(FF0->ID, (void**)&FF0, Loop->FacesInLoop,
                                    NULL, false));
        Loop->head = FF0;
        Loop->tail = FF0->prev;
        mg_free_front_face(FF1);
        mg_free_front_face(FF2);
      }
      else {//face 1 is new
        //create a FrontFace structure for faceID1
        call(mg_alloc((void**)&FF1, 1, sizeof(mg_FrontFace)));
        FF0 = FF2->next;
        FF0Next = FF0->next;
        FF2Prev = FF2->prev;
        //setup face1
        FF1->ID = faceID1;
        FF1->iloop = FF2->iloop;
        FF1->face = Mesh->Face[faceID1];
        FF1->next = FF0->next;
        FF1->prev = FF2->prev;
        //stitch loop
        FF0Next->prev = FF1;
        FF2Prev->next = FF1;
        //remove faces from loop
        call(mg_rm_frm_ord_data_list(faceID0, Loop->FacesInLoop));
        call(mg_rm_frm_ord_data_list(faceID2, Loop->FacesInLoop));
        call(mg_add_2_ord_data_list(FF1->ID, (void**)&FF1, Loop->FacesInLoop,
                                    NULL, false));
        Loop->head = FF1;
        Loop->tail = FF1->prev;
        mg_free_front_face(FF0);
        mg_free_front_face(FF2);
      }
    }
    else if (face0new && face1new) {//no new node
      //check if merging 2 loops or spliting one loop
      //check which loops contain nodeID2
      FF2Prev = FF2->prev;
      FF2Next = FF2->next;
      call(mg_find_node_in_front(nodeID2, Mesh, Front, &Node2FFace));
      sameloop = false;
      for (f = 0; f < Node2FFace->nEntry; f++){
        FFace = (mg_FrontFace*)(Node2FFace->Data[f]);
        faceID = FFace->ID;
        if (FFace->iloop == FF2->iloop) {
          sameloop = true;
          break;
        }
      }
      call(mg_alloc((void**)&FF0, 1, sizeof(mg_FrontFace)));
      call(mg_alloc((void**)&FF1, 1, sizeof(mg_FrontFace)));
      if (sameloop) { //we are splitting the loop in 2
        /* Diagram of splitting a loop into 2:
         f1Next     f0Prev
         *----------*----------*
         |         / \         |
         |    f1ID/ ^ \f0ID    |
         |       /  |  \       |
         *------*-------*------*
         f2Prev   f2ID    f2Next
         */
        //First step: Remove faceID2 from loop
        leftloop = FF2->iloop;
        rightloop = Front->nloop;
        if (Mesh->Face[faceID]->node[0] == nodeID2){
          FF1Next = FFace;
          FF0Prev = FF1Next->prev;
        }
        if (Mesh->Face[faceID]->node[1] == nodeID2){
          FF0Prev = FFace;
          FF1Next = FF0Prev->next;
        }
        //Work on left loop
        //setup face1
        FF1->ID = faceID1;
        FF1->iloop = leftloop;
        FF1->face = Mesh->Face[faceID1];
        FF1->next = FF1Next;
        FF1->prev = FF2Prev;
        //stitch left loop
        FF2Prev->next = FF1;
        FF1Next->prev = FF1;
        Loop = Front->loop[leftloop];
        Loop->head = FF1;
        Loop->tail = FF2Prev;
        call(mg_rm_frm_ord_data_list(faceID2, Loop->FacesInLoop));
        call(mg_add_2_ord_data_list(FF1->ID, (void**)&FF1, Loop->FacesInLoop,
                                    NULL, false));
        //Work on right loop
        //setup face1
        FF0->ID = faceID0;
        FF0->iloop = rightloop;
        FF0->face = Mesh->Face[faceID0];
        FF0->next = FF2Next;
        FF0->prev = FF0Prev;
        //stitch right loop
        FF0Prev->next = FF0;
        FF2Next->prev = FF0;
        call(mg_alloc((void**)&LoopNew, 1, sizeof(mg_Loop)));
        call(mg_alloc((void**)&LoopNew->FacesInLoop, 1,
                      sizeof(mg_OrderedDataList)));
        mg_init_ord_data_list(LoopNew->FacesInLoop, sizeof(mg_FrontFace));
        LoopNew->head = FF0;
        LoopNew->tail = FF0Prev;
        call(mg_add_2_ord_data_list(FF0->ID, (void**)&FF0,LoopNew->FacesInLoop,
                                    NULL, false));
        //loop around loopnew and change appropriate values
        FFace = FF0->next;
        while (FFace != LoopNew->head) {
          FFace->iloop = rightloop;
          //remove from left loop
          call(mg_rm_frm_ord_data_list(FFace->ID, Loop->FacesInLoop));
          //add to new loop
          call(mg_add_2_ord_data_list(FFace->ID, (void**)&FFace,LoopNew->FacesInLoop,
                                      NULL, false));
          FFace = FFace->next;
        }
        //add new loop to front
        Front->nloop++;
        call(mg_realloc((void**)&Front->loop, Front->nloop, sizeof(mg_Loop)));
        Front->loop[Front->nloop-1] = LoopNew;
        mg_free_front_face(FF2);
      }
      else {//merging 2 loops
        if (Node2FFace->nEntry != 2) return error(err_LOGIC_ERROR);
        FFace = (mg_FrontFace*)(Node2FFace->Data[0]);
        faceID = FFace->ID;
        if (Mesh->Face[faceID]->node[0] == nodeID2){
          FF1Next = FFace;
          FF0Prev = FF1Next->prev;
        }
        else if (Mesh->Face[faceID]->node[1] == nodeID2){
          FF0Prev = FFace;
          FF1Next = FF0Prev->next;
        }
        else return error(err_LOGIC_ERROR);
        oldloop = FF2->iloop;
        newloop = FF1Next->iloop;
        //setup FF0
        FF0->ID = faceID0;
        FF0->iloop = newloop;
        FF0->face = Mesh->Face[faceID0];
        FF0->prev = FF0Prev;
        FF0->next = FF2Next;
        //setup FF1
        FF1->ID = faceID1;
        FF1->iloop = newloop;
        FF1->face = Mesh->Face[faceID1];
        FF1->prev = FF2Prev;
        FF1->next = FF1Next;
        //stitch loops
        FF2Prev->next = FF1;
        FF1Next->prev = FF1;
        FF2Next->prev = FF0;
        FF0Prev->next = FF0;
        //clean up oldloop
        Loop = Front->loop[oldloop];
        LoopNew = Front->loop[newloop];
        LoopNew->head = FF1;
        LoopNew->tail = FF1->prev;
        Loop->head = Loop->tail = NULL;
        call(mg_rm_frm_ord_data_list(faceID2, Loop->FacesInLoop));
        mg_free_front_face(FF2);
        //add faces to new loop
        call(mg_add_2_ord_data_list(FF0->ID, (void**)&FF0, LoopNew->FacesInLoop,
                                    NULL, false));
        call(mg_add_2_ord_data_list(FF1->ID, (void**)&FF1, LoopNew->FacesInLoop,
                                    NULL, false));
        FFace = FF2Next;
        while (FFace != FF1) {
          call(mg_rm_frm_ord_data_list(FFace->ID, Loop->FacesInLoop));
          FFace->iloop = newloop;
          call(mg_add_2_ord_data_list(FFace->ID, (void**)&FFace,
                                      LoopNew->FacesInLoop,NULL, false));
          FFace = FFace->next;
        }
        mg_free_ord_data_list(Loop->FacesInLoop);
        mg_free_ord_data_list(Node2FFace);
        mg_free(Node2FFace);
      }
    }
    else {
      //no new node and no new faces means the loop disappears
      iloop = FF2->iloop;
      Loop = Front->loop[iloop];
      if (Loop->FacesInLoop->nEntry != 3) return error(err_LOGIC_ERROR);
      FF2Next = FF2->next;
      if (FF2Next->ID != faceID0) return error(err_LOGIC_ERROR);
      FF2Prev = FF2->prev;
      if (FF2Prev->ID != faceID1) return error(err_LOGIC_ERROR);
      //now let's remove the loop from the front
      mg_free_loop(Loop);
    }
  }
  
  return err_OK;
}

/******************************************************************/
/* function: mg_bld_tri_frm_close_pts */
/* build as many triangles possible from list of nearby points */
int mg_bld_tri_frm_close_pts(mg_Mesh *Mesh, mg_Front *Front,
                             mg_FrontFace *SelfFace, double *rhomax,
                             bool IsoFlag, mg_List *CloseNodes,
                             bool *success, double *radii)
{
  int ierr, in, nodeID0, d, dim, nodeID1, iloop;
  double proj, radius, center[2], dist2, xint[2], X0[4], X1[4];
  bool intersect=false, first;
  mg_FaceData *face = SelfFace->face;
  mg_List NodesOnLeft;
  mg_FrontFace *FFace;
  mg_Loop *Loop;
  
  if (!IsoFlag) return error(err_NOT_SUPPORTED);
  
  dim = Mesh->Dim;
  NodesOnLeft.nItem = 0;
  NodesOnLeft.Item = NULL;
  (*success) = false;
  //get nodes on left side of selfface
  for (in = 0; in < CloseNodes->nItem; in++) {
    nodeID0 = CloseNodes->Item[in];
    proj = 0.0;
    for (d = 0; d < dim; d++){
      proj+=face->normal[d]*(Mesh->Coord[nodeID0*dim+d]-face->centroid[d]);
    }
    if (proj > 1e-5) {
      call(mg_add_2_ord_set(nodeID0, &NodesOnLeft.nItem, &NodesOnLeft.Item,
                            NULL, false));
    }
  }
  //loop over nodes on left side and build the smallest (in radius) possible triangle
  if (NodesOnLeft.nItem > 0) {
    nodeID0 = NodesOnLeft.Item[0];
    call(mg_build_circle_frm_face(Mesh, face, Mesh->Coord+nodeID0*dim, center,
                                  &radius));
    for (in = 1; in < NodesOnLeft.nItem; in++) {
      nodeID1 = NodesOnLeft.Item[in];
      dist2 = (Mesh->Coord[nodeID1*dim]-center[0])*(Mesh->Coord[nodeID1*dim]-center[0]);
      dist2+= (Mesh->Coord[nodeID1*dim+1]-center[1])*(Mesh->Coord[nodeID1*dim+1]-center[1]);
      if (dist2 < radius*radius) {
        nodeID0 = nodeID1;
        call(mg_build_circle_frm_face(Mesh, face, Mesh->Coord+nodeID0*dim, center,
                                      &radius));
      }
    }
    //check is radius is acceptable
    if (radius <= (*rhomax)) {
      for (iloop = 0; iloop < Front->nloop; iloop++) {
        Loop = Front->loop[iloop];
        if (Loop->FacesInLoop->nEntry == 0) continue;
        FFace = Loop->head;
        first = true;
        intersect = false;
        while ((first || FFace != Loop->head) && !intersect) {
          first = false;
          if (FFace != SelfFace && nodeID0 != FFace->face->node[0] && nodeID0 != FFace->face->node[1]){
            //compute intersection between front face and triangle height vector
            X0[0] = Mesh->Coord[FFace->face->node[0]*dim];
            X0[1] = Mesh->Coord[FFace->face->node[0]*dim+1];
            X0[2] = Mesh->Coord[FFace->face->node[1]*dim];
            X0[3] = Mesh->Coord[FFace->face->node[1]*dim+1];
            //trinagle height vector
            X1[0] = SelfFace->face->centroid[0];
            X1[1] = SelfFace->face->centroid[1];
            X1[2] = Mesh->Coord[nodeID0*dim];
            X1[3] = Mesh->Coord[nodeID0*dim+1];
            intersect = mg_edges_intersect(X0,X1, xint);
          }
          FFace = FFace->next;
        }
        if (intersect)
          break;
      }
      if (!intersect){
        //build triangle and update front
        call(mg_tri_frm_face_node(Mesh, SelfFace, nodeID0, NULL));
        call(mg_update_front(Mesh, Front, SelfFace));
        (*success) = true;
      }
    }
  }
  
  //substitute CloseNodes by nodes on the left (not the best way to do this)
  CloseNodes->nItem = NodesOnLeft.nItem;
  mg_free((void*)CloseNodes->Item);
  CloseNodes->Item = NodesOnLeft.Item;
  
  (*radii) = radius;//assuming isotropic here
  
  return err_OK;
}

/******************************************************************/
/* function: mg_find_face_in_frt */
/* finds a face with "faceID" in "Front" */
int mg_find_face_in_frt(mg_Front *Front, int faceID, mg_FrontFace **pFFace)
{
  int ierr, iloop, rank;
  mg_FrontFace *FFace = NULL;
  for (iloop = 0; iloop < Front->nloop; iloop++) {
    if (Front->loop[iloop]->FacesInLoop->nEntry == 0) continue;
    ierr = mg_binary_search(faceID, Front->loop[iloop]->FacesInLoop->Entry, 0,
                            Front->loop[iloop]->FacesInLoop->nEntry-1, &rank);
    if (ierr == err_OK) {//found face
      FFace = (mg_FrontFace*)Front->loop[iloop]->FacesInLoop->Data[rank];
      break;
    }
    else if (ierr != err_NOT_FOUND) return error(ierr);
  }
  if (FFace == NULL) return error(err_NOT_FOUND);
  
  (*pFFace) = FFace;
  
  return err_OK;
}

/******************************************************************/
/* function: mg_rm_broken_elems */
/* removes elements from mesh */
int mg_rm_broken_elems(mg_Mesh *Mesh, mg_Front *Front,
                       mg_List *BrokenElems,
                       mg_List *CandidateNodes)
{
  int ierr, e, elem, f, faceID, nborID, fIDX2rm[3], nf2rm, fIDX2kp[3], nf2kp;
  int idx, t, n, nodeID, node2rm;
  mg_FaceData *face;
  mg_FrontFace *FFace = NULL, *FFacePrev, *FFaceNext, *FFaceNew, *FFaceStale;
  
  //note: the first nElemInFront of this list are ordered
  for (e = 0; e < BrokenElems->nItem; e++) {
    elem = BrokenElems->Item[e];
    //loop over faces and fix connectivities
    nf2rm = nf2kp = 0;
    for (f = 0; f < Mesh->Elem[elem].nNode; f++) {
      faceID = Mesh->Elem[elem].face[f];
      //set pointer to face
      face = Mesh->Face[faceID];
      if (face->elem[LEFTNEIGHINDEX] == HOLLOWNEIGHTAG){
        //this is a front face, mark index for removal
        fIDX2rm[nf2rm] = f;
        nf2rm++;
      }
      else {
        fIDX2kp[nf2kp] = f;
        nf2kp++;
      }
    }
    if (nf2rm == 3 || nf2kp == 3) return error(err_LOGIC_ERROR);
    if (nf2rm == 0 || nf2kp == 0) return error(err_LOGIC_ERROR);
    /***************************************************************************/
    //two distinct cases
    if (nf2rm == 1) {//one face to remove
      //stitch front
      faceID = Mesh->Elem[elem].face[fIDX2rm[0]];
      call(mg_find_face_in_frt(Front, faceID, &FFace));
      call(mg_rm_frm_ord_data_list(faceID, Front->loop[FFace->iloop]->FacesInLoop));
      FFaceNext = FFace->next;
      FFacePrev = FFace->prev;
      //make sure we keep correct orientation of the loop
      //note: the nodal order for the front face is such that the
      //normal point to the hollow
      f = Mesh->Elem[elem].face[fIDX2kp[1]];
      if (Mesh->Face[faceID]->node[0] == Mesh->Face[f]->node[0] ||
          Mesh->Face[faceID]->node[0] == Mesh->Face[f]->node[1]) {
        swap(fIDX2kp[0], fIDX2kp[1], t);
      }
      FFace->ID = Mesh->Elem[elem].face[fIDX2kp[0]];
      FFace->face = Mesh->Face[FFace->ID];
      call(mg_alloc((void**)&FFaceNew, 1, sizeof(mg_FrontFace)));
      FFace->next = FFaceNew;
      FFaceNew->prev = FFace;
      FFaceNew->next = FFaceNext;
      FFaceNext->prev = FFaceNew;
      FFaceNew->iloop = FFace->iloop;
      //reset tail and head
      Front->loop[FFace->iloop]->head = FFaceNew;
      Front->loop[FFace->iloop]->tail = FFace;
      //new front face is the second face to keep
      FFaceNew->ID = Mesh->Elem[elem].face[fIDX2kp[1]];
      FFaceNew->face = Mesh->Face[FFaceNew->ID];
      call(mg_add_2_ord_data_list(FFace->ID, (void**)&FFace,
                                  Front->loop[FFace->iloop]->FacesInLoop,
                                  NULL, false));
      call(mg_add_2_ord_data_list(FFaceNew->ID, (void**)&FFaceNew,
                                  Front->loop[FFaceNew->iloop]->FacesInLoop,
                                  NULL, false));
      
      node2rm = -1;
      //add new front node (opposite to fIDX2rm[0]) to list of candidate nodes
      nodeID = Mesh->Elem[elem].node[fIDX2rm[0]];
      call(mg_add_2_ord_set(nodeID, &CandidateNodes->nItem,&CandidateNodes->Item,
                            NULL, false));
    }
    else {//2faces and one node to remove
      //stitch front
      faceID = Mesh->Elem[elem].face[fIDX2rm[0]];
      call(mg_find_face_in_frt(Front, faceID, &FFace));
      call(mg_rm_frm_ord_data_list(faceID, Front->loop[FFace->iloop]->FacesInLoop));
      FFaceNext = FFace->next;
      FFaceStale = FFace->prev;
      if (FFaceStale->ID != Mesh->Elem[elem].face[fIDX2rm[1]])
        return error(err_LOGIC_ERROR);
      call(mg_rm_frm_ord_data_list(FFaceStale->ID,
                                   Front->loop[FFaceStale->iloop]->FacesInLoop));
      FFacePrev = FFaceStale->prev;
      FFacePrev->next = FFace;
      FFace->prev = FFacePrev;
      mg_free_front_face(FFaceStale);
      //substitute old front face by only face to keep
      FFace->ID = Mesh->Elem[elem].face[fIDX2kp[0]];
      FFace->face = Mesh->Face[FFace->ID];
      call(mg_add_2_ord_data_list(FFace->ID, (void**)&FFace,
                                  Front->loop[FFace->iloop]->FacesInLoop,
                                  NULL, false));
      //NOTE: node to remove is the one opposite to the face to keep
      node2rm = Mesh->Elem[elem].node[fIDX2kp[0]];
      //remove node from candidate nodes if it was a candidate
      call(mg_rm_frm_ord_set(node2rm,&CandidateNodes->nItem,&CandidateNodes->Item, 1, &t));
      //add nodes from face to keep to list of candidate nodes
      faceID = Mesh->Elem[elem].face[fIDX2kp[0]];
      for (n = 0; n < Mesh->Face[faceID]->nNode; n++) {
        nodeID = Mesh->Face[faceID]->node[n];
        call(mg_add_2_ord_set(nodeID, &CandidateNodes->nItem,
                              &CandidateNodes->Item, NULL, false));
      }
    }
    /***************************************************************************/
    //update face connections and possibly swap nodes and reverse normal
    for (f = 0; f < nf2kp; f++) {
      faceID = Mesh->Elem[elem].face[fIDX2kp[f]];
      face = Mesh->Face[faceID];
      if (face->elem[LEFTNEIGHINDEX] == elem) {
        face->elem[LEFTNEIGHINDEX] = HOLLOWNEIGHTAG;
      }
      else if (face->elem[RIGHTNEIGHINDEX] == elem){
        //need to flip normal and and swap nodes
        swap(face->node[0], face->node[1], t);
        face->normal[0] *= -1.0;
        face->normal[1] *= -1.0;
      }//has to be either left or right
      else return error(err_LOGIC_ERROR);
      //update neighbor lists in neighbors
      nborID = Mesh->Elem[elem].nbor[fIDX2kp[f]];
      if (nborID >= 0) { //not a boundary
        mg_check_exist(faceID, 3, Mesh->Elem[nborID].face, &idx);
        if (idx < 0) return error(err_MESH_ERROR);
        Mesh->Elem[nborID].nbor[idx] = HOLLOWNEIGHTAG;
      }
    }
    
    //update node2face
    for (f = 0; f < nf2rm; f++) {
      faceID = Mesh->Elem[elem].face[fIDX2rm[f]];
      for (n = 0; n < Mesh->Face[faceID]->nNode; n++) {
        nodeID = Mesh->Face[faceID]->node[n];
        call(mg_rm_frm_ord_set(faceID, &Mesh->Node2Face[nodeID].nItem,
                               &Mesh->Node2Face[nodeID].Item, 1, &t));
        if (t != 1) error(err_MESH_ERROR);
        mg_free((void*)Mesh->Face[faceID]->centroid);
        mg_free((void*)Mesh->Face[faceID]->normal);
        Mesh->Face[faceID]->centroid = NULL;
        Mesh->Face[faceID]->normal = NULL;
      }
    }
    //update stack of removed mesh components
    //add elem to stack
    call(mg_add_2_ord_set(elem, &Mesh->Stack->Elem->nItem,&Mesh->Stack->Elem->Item,
                          NULL, false));
    //reduce number of valid elements
    Mesh->nElem--;
    //add node to remove (if any) to stack
    if (node2rm >= 0){
      call(mg_add_2_ord_set(node2rm, &Mesh->Stack->Node->nItem,&Mesh->Stack->Node->Item,
                            NULL, false));
      Mesh->nNode--;
    }
    //add faces to remove to stack
    for (f = 0; f < nf2rm; f++) {
      call(mg_add_2_ord_set(Mesh->Elem[elem].face[fIDX2rm[f]],
                            &Mesh->Stack->Face->nItem,&Mesh->Stack->Face->Item,
                            NULL, false));
      Mesh->nFace--;
    }
  }//BrokenElems.nItem
  
  BrokenElems->nItem = 0;
  mg_free((void*)BrokenElems->Item);
  
  return err_OK;
}

/******************************************************************/
/* function: mg_neigh_srch_brkn_tri */
/* searches for broken triangles via a neighbor search */
int mg_neigh_srch_brkn_tri(mg_Mesh *Mesh, mg_List *BrokenTri,
                           double *newcoord, double *rhomax)
{
  int ierr, e, elem, n, nbor, oppnode, dim, d;
  int *BrokenTriList;
  double *coord, center[2], radius, dist2;
  mg_FaceData *face;
  
  /* Local node arrangement diagram:
   |     2
   |     *
   |  1 / \ 0
   |   / e \
   | 0*-----*1
   |     2
   */
  
  call(mg_alloc((void**)&BrokenTriList, BrokenTri->nItem, sizeof(int)));
  memcpy(BrokenTriList, BrokenTri->Item, BrokenTri->nItem*sizeof(int));
  dim = Mesh->Dim;
  for (e = 0; e < BrokenTri->nItem; e++) {
    //note: this list will keep expanding while new broken triangles are found
    elem = BrokenTriList[e];
    //check elem's neighbors
    for (n = 0; n < Mesh->Elem[elem].nNode; n++) {
      //note: number of neighbors is the same as number of nodes
      nbor = Mesh->Elem[elem].nbor[n];
      if (nbor >= 0){//neighbor is not a boundary
        //check if nbor has been listed before
        ierr = mg_binary_search(nbor, BrokenTri->Item, 0, BrokenTri->nItem-1,
                                NULL);
        if (ierr == err_NOT_FOUND){
          //get first face in neighbor and its opposing node to compute the circumcircle
          face = Mesh->Face[Mesh->Elem[nbor].face[0]];
          oppnode = Mesh->Elem[nbor].node[0]; //look at diagram above
          coord = Mesh->Coord+oppnode*dim;
          call(mg_build_circle_frm_face(Mesh, face, coord, center, &radius));
          dist2 = 0.0;
          for (d = 0; d < dim; d++)
            dist2 += (newcoord[d]-center[d])*(newcoord[d]-center[d]);
          if (dist2 <= radius*radius) {
            //newcoord is inside nbor's circumcircle
            //add to ordered list to facilitate search
            call(mg_add_2_ord_set(nbor, &BrokenTri->nItem, &BrokenTri->Item,
                                  NULL, false));
            //non-ordered list
            call(mg_realloc((void**)&BrokenTriList, BrokenTri->nItem, sizeof(int)));
            //add new item at end of list
            BrokenTriList[BrokenTri->nItem-1] = nbor;
            //keep track of maximum allowable circumradius
            if (radius > (*rhomax)) (*rhomax) = radius;
          }
        }
        else if (ierr != err_OK) return error(ierr);
      }
    }
  }
  
  mg_free((void*)BrokenTri->Item);
  //note: the first elements are adjacent to front
  BrokenTri->Item = BrokenTriList;
  
  return err_OK;
}

/******************************************************************/
/* function: mg_cand_nds_2_cand_fcs */
/* converts a list of candidate nodes to candidate faces for
 forming triangles */
int mg_cand_nds_2_cand_fcs(mg_Mesh *Mesh, mg_Front *Front,
                           mg_List *CandidateNodes,
                           mg_OrderedDataList **pCandidateFaces)
{
  int ierr, inode, nodeID, iface, othernodeID;
  mg_OrderedDataList *Node2FFace = NULL, *CandidateFaces;
  mg_FrontFace *FFace;
  
  call(mg_alloc((void**)&CandidateFaces, 1, sizeof(mg_OrderedDataList)));
  mg_init_ord_data_list(CandidateFaces, sizeof(mg_FrontFace));
  
  for (inode = 0; inode < CandidateNodes->nItem; inode++) {
    nodeID = CandidateNodes->Item[inode];
    call(mg_find_node_in_front(nodeID, Mesh, Front, &Node2FFace));
    for (iface = 0; iface < Node2FFace->nEntry; iface++) {
      FFace = (mg_FrontFace*)Node2FFace->Data[iface];
      if (FFace->face->node[0] == nodeID)
        othernodeID = FFace->face->node[1];
      else
        othernodeID = FFace->face->node[0];
      //if othernodeID is a candidate, then faceID is a candidate
      ierr = mg_binary_search(othernodeID, CandidateNodes->Item, 0,
                              CandidateNodes->nItem-1, NULL);
      if (ierr == err_NOT_FOUND){
        call(mg_add_2_ord_data_list(FFace->ID, (void**)&FFace,
                                    CandidateFaces, NULL, false));
      }
      else if (ierr != err_OK) return error(ierr);
    }
    if (Node2FFace->nEntry > 0)//non empty
      mg_free_ord_data_list(Node2FFace);
  }
  
  (*pCandidateFaces) = CandidateFaces;
  
  return err_OK;
}

/******************************************************************/
/* function: mg_add_new_node */
/* adds a node to the mesh that forms a triangle with "ActiveFace"*/
int
mg_add_new_node(mg_Mesh *Mesh, mg_Front *Front, mg_FrontFace *ActiveFace,
                mg_List *CandidateNodes, bool isoflag, double rhomax,
                double c, bool *success)
{
  int ierr, newnodeID, d, iloop, elem, idx, oppnode, dim, nBrokenTriInFront;
  int icface, nsuccess, t, n;
  double UpperRBound, LowerRBound, newcoord[3], *coord, center[2], radius, dist2;
  double newrhomax, xint[2], X0[4], X1[4];
  bool NodeFromStack, first;
  mg_List *BrokenTri, *BrokenFFace;
  mg_Loop *Loop;
  mg_FrontFace *FFace;
  mg_FaceData *gface;
  mg_OrderedDataList *CandidateFaces;
  
  if (!isoflag) return error(err_NOT_SUPPORTED);
  
  UpperRBound = c*rhomax;//circumradius for closest candidate node
  LowerRBound = ActiveFace->face->area/2;//inscribed radius for edge
  if (LowerRBound > UpperRBound) return error(err_LOGIC_ERROR);
  //new node coordinates
  dim = Mesh->Dim;
  //call(mg_alloc((void**)&newcoord, dim, sizeof(double)));
  for (d = 0; d < dim; d++)
    newcoord[d] = ActiveFace->face->centroid[d]+(HALFSQRT3*UpperRBound+SQRT3*LowerRBound)/(2.0)*ActiveFace->face->normal[d];
  //loop over front and check for broken triangles
  call(mg_alloc((void**)&BrokenTri,1,sizeof(mg_List)));
  mg_init_list(&BrokenTri[0]);
  call(mg_alloc((void**)&BrokenFFace,1,sizeof(mg_List)));
  mg_init_list(&BrokenFFace[0]);
  for (iloop = 0; iloop < Front->nloop; iloop++) {
    Loop = Front->loop[iloop];
    if (Loop->FacesInLoop->nEntry == 0) continue;
    FFace = Loop->head;
    first = true;
    while (first || FFace != Loop->head) {
      first = false;
      if (FFace != ActiveFace){
        elem = FFace->face->elem[RIGHTNEIGHINDEX];
        if (elem >= 0) {//not a boundary
          mg_check_exist(FFace->ID, Mesh->Elem[elem].nNode,
                         Mesh->Elem[elem].face, &idx);
          if (idx < 0) return error(err_MESH_ERROR);
          //point to opposite node to FFace
          oppnode = Mesh->Elem[elem].node[idx];
          coord = Mesh->Coord+oppnode*dim;
          call(mg_build_circle_frm_face(Mesh, FFace->face, coord, center, &radius));
          dist2 = 0.0;
          for (d = 0; d < dim; d++)dist2 += (newcoord[d]-center[d])*(newcoord[d]-center[d]);
          if (dist2 <= radius*radius) {
            //elem is not delaunay anymore due to new point
            call(mg_add_2_ord_set(elem,&BrokenTri->nItem,&BrokenTri->Item,NULL, false));
          }
        }
        else {
          //compute intersection between front face and triangle height vector
          X0[0] = Mesh->Coord[FFace->face->node[0]*dim];
          X0[1] = Mesh->Coord[FFace->face->node[0]*dim+1];
          X0[2] = Mesh->Coord[FFace->face->node[1]*dim];
          X0[3] = Mesh->Coord[FFace->face->node[1]*dim+1];
          //trinagle height vector
          X1[0] = ActiveFace->face->centroid[0];
          X1[1] = ActiveFace->face->centroid[1];
          X1[2] = newcoord[0];
          X1[3] = newcoord[1];
          if (mg_edges_intersect(X0,X1, xint)) {
            call(mg_add_2_ord_set(FFace->ID,&BrokenFFace->nItem,&BrokenFFace->Item,NULL, false));
            if (FFace->face->area > rhomax)
              rhomax = FFace->face->area/2.0;
            for (n = 0; n < FFace->face->nNode; n++)
              call(mg_add_2_ord_set(FFace->face->node[n],&CandidateNodes->nItem,&CandidateNodes->Item,
                                    NULL, false));
          }
        }
      }
      FFace = FFace->next;
    }
  }
  
  if (BrokenFFace->nItem > 0){
    call(mg_bld_tri_frm_close_pts(Mesh, Front, ActiveFace, &rhomax, isoflag,
                                  CandidateNodes, success, &radius));
    if ((*success)){
      BrokenFFace->nItem = 0;
      mg_free((void*)BrokenFFace->Item);
      return err_OK;
    }
  }
  
  NodeFromStack = false;
  if (Mesh->Stack->Node->nItem > 0){
    newnodeID = Mesh->Stack->Node->Item[0];
    NodeFromStack = true;
  }
  else
    newnodeID = Mesh->nNode;
  
  nBrokenTriInFront = BrokenTri->nItem;
  if (nBrokenTriInFront == 0 && BrokenFFace->nItem == 0) {//no bronken triangles, accept point and update front
    //build triangle and update front
    call(mg_tri_frm_face_node(Mesh, ActiveFace, newnodeID, newcoord));
    call(mg_update_front(Mesh, Front, ActiveFace));
    (*success) = true;
  }
  else {
    newrhomax = rhomax;
    //initiate neighbor search
    call(mg_neigh_srch_brkn_tri(Mesh, BrokenTri, newcoord, &newrhomax));
    if (newrhomax > rhomax) {
      //update list of close nodes
      call(mg_nodes_frnt_dist(Mesh, Front, ActiveFace, &newrhomax,2.0,
                              isoflag,CandidateNodes));
    }
    /*remove intersected triangles, update list of candidate nodes, and include
     removed mesh components in the mesh stack*/
    call(mg_rm_broken_elems(Mesh, Front, BrokenTri, CandidateNodes));
    //convert candidate nodes into cadidate front faces
    call(mg_cand_nds_2_cand_fcs(Mesh, Front, CandidateNodes, &CandidateFaces));
    //loop over candidate faces and build triangles
    
    nsuccess = 0;
    for (icface = 0; icface < CandidateFaces->nEntry; icface++) {
      FFace = (mg_FrontFace*)CandidateFaces->Data[icface];
      gface = FFace->face;
      call(mg_build_circle_frm_face(Mesh, gface, newcoord, center, &radius));
      if (radius < newrhomax){
        //form circle
        call(mg_tri_frm_face_node(Mesh, FFace, newnodeID, (nsuccess == 0)?
                                  newcoord:NULL));
        if (nsuccess == 0 && NodeFromStack == true){
          //first success, it means newnode is not new anymore
          call(mg_rm_frm_ord_set(newnodeID, &Mesh->Stack->Node->nItem,
                                 &Mesh->Stack->Node->Item, 1, &t));
          if (t != 1) return error(err_MESH_ERROR);
        }
        nsuccess++;
        call(mg_update_front(Mesh, Front, FFace));
      }
    }
    if (nsuccess > 0) (*success) = true;
  }
  
  return err_OK;
}


/******************************************************************/
/* function: mg_advance_front */
/* advances the mesh front */
int mg_advance_front(mg_Mesh *Mesh, mg_Front *Front)
{
  int ierr;
  double *rhomax, radius;
  bool isoflag, success;
  mg_FrontFace *SeedFace;
  mg_List *CloseNodes;
  
  call(mg_find_seed_face(Front, &SeedFace));
  call(mg_comp_mesh_size(Mesh, SeedFace, &rhomax, &isoflag));
  call(mg_alloc((void**)&CloseNodes, 1, sizeof(mg_List)));
  call(mg_nodes_frnt_dist(Mesh, Front, SeedFace, rhomax, 2.0,
                          isoflag,CloseNodes));
  call(mg_bld_tri_frm_close_pts(Mesh, Front, SeedFace, rhomax, isoflag,
                                CloseNodes, &success, &radius));
  if (!success) {
    call(mg_add_new_node(Mesh, Front, SeedFace, CloseNodes, isoflag,
                         rhomax[0], 2.0,&success));
    if (!success) return error(err_MESH_ERROR);
  }
  
  CloseNodes[0].nItem = 0;
  mg_free((void*)CloseNodes[0].Item);
  mg_free((void*)CloseNodes);
  mg_free((void*)rhomax);
  
  return err_OK;
}

/******************************************************************/
/* function: mg_front_empty */
/* checks if front is empty */
bool mg_front_empty(mg_Front *Front)
{
  int iloop;
  mg_Loop *Loop;
  
  for (iloop = 0; iloop < Front->nloop; iloop++){
    Loop = Front->loop[iloop];
    if (Loop != NULL)
      if (Loop->head != NULL)
        return false;
  }
  
  return true;
}

/******************************************************************/
/* function: mg_prealloc_msh_comp */
/* preallocates mesh components to avoid too much reallocation */
int mg_prealloc_msh_comp(mg_Mesh *Mesh, int nElem, int nFace, int nNode)
{
  int ierr, i, k;
  mg_FaceData *GFace;
  
  if (Mesh->Stack->Elem->nItem != 0 || Mesh->Stack->Elem->Item != NULL)
    return error(err_INPUT_ERROR);
  
  if (nElem > Mesh->nElem){
    call(mg_realloc((void**)&Mesh->Elem, nElem, sizeof(mg_ElemData)));
    call(mg_realloc((void**)&Mesh->Stack->Elem->Item, nElem-Mesh->nElem,
                    sizeof(int)));
    for (i = Mesh->nElem; i < nElem; i++) {
      Mesh->Elem[i].nNode = 3;
      call(mg_alloc((void**)&Mesh->Elem[i].node, Mesh->Elem[i].nNode,
                    sizeof(int)));
      call(mg_alloc((void**)&Mesh->Elem[i].face, Mesh->Elem[i].nNode,
                    sizeof(int)));
      call(mg_alloc((void**)&Mesh->Elem[i].nbor, Mesh->Elem[i].nNode,
                    sizeof(int)));
      k = Mesh->Stack->Elem->nItem;
      Mesh->Stack->Elem->Item[k] = i;
      Mesh->Stack->Elem->nItem++;
    }
  }
  
  if (Mesh->Stack->Face->nItem != 0 || Mesh->Stack->Face->Item != NULL)
    return error(err_INPUT_ERROR);
  
  if (nFace > Mesh->nFace){
    call(mg_realloc((void**)&Mesh->Face, nFace, sizeof(mg_FaceData)));
    call(mg_realloc((void**)&Mesh->Stack->Face->Item, nFace-Mesh->nFace,
                    sizeof(int)));
    for (i = Mesh->nFace; i < nFace; i++) {
      call(mg_alloc((void**)&GFace, 1, sizeof(mg_FaceData)));
      Mesh->Face[i] = GFace;
      Mesh->Face[i]->area = -1.0;
      Mesh->Face[i]->centroid = NULL;
      Mesh->Face[i]->normal = NULL;
      Mesh->Face[i]->nNode = 2;
      call(mg_alloc((void**)&Mesh->Face[i]->node, Mesh->Face[i]->nNode,
                    sizeof(int)));
      k = Mesh->Stack->Face->nItem;
      Mesh->Stack->Face->Item[k] = i;
      Mesh->Stack->Face->nItem++;
    }
  }
  
  if (Mesh->Stack->Node->nItem != 0 || Mesh->Stack->Node->Item != NULL)
    return error(err_INPUT_ERROR);
  
  if (nNode > Mesh->nNode){
    call(mg_realloc((void**)&Mesh->Coord, nNode*Mesh->Dim,
                    sizeof(double)));
    call(mg_realloc((void**)&Mesh->Stack->Node->Item, nNode-Mesh->nNode,
                    sizeof(int)));
    call(mg_realloc((void**)&Mesh->Node2Elem, nNode, sizeof(mg_List)));
    call(mg_realloc((void**)&Mesh->Node2Face, nNode, sizeof(mg_List)));
    for (i = Mesh->nNode; i < nNode; i++) {
      mg_init_list(&Mesh->Node2Elem[i]);
      k = Mesh->Stack->Node->nItem;
      Mesh->Stack->Node->Item[k] = i;
      Mesh->Stack->Node->nItem++;
    }
  }
  
  return err_OK;
}

/******************************************************************/
/* Main program */
int main(int argc, char *argv[])
{
  int ierr, len, i;
  char ParFile[MAXSTRLEN], *InFile, *OutFile,*pext;
  char MeshName[MAXSTRLEN];
  mg_Mesh *Mesh;
  mg_Front Front;
  mg_Geometry *Geo;
  
  /* Check number of arguments */
  if( argc != 2 ){
    printf("Usage:\n");
    printf("2dmg <parfile>\n");
    printf("\n");
    printf("Where <parfile> is the name of the parameter file.\n");
    printf("\n");
    return error(err_INPUT_ERROR);
  }
  
  // Parameter file name
  strcpy(ParFile, argv[1]);
  //parse input file an put in hash table
  call(mg_read_input_file(ParFile));
  
  
  //Read boundary discretization or geometry file
  call(mg_get_input_char("BoundaryMesh", &InFile));
  len = (int)strlen(InFile);
  if (strcmp(InFile, "None") != 0) {
    if (len < 5) {
      printf("Invalid BoundaryMesh filename");
      return err_INPUT_ERROR;
    }
    pext = InFile+len-5;
    if (strncmp(pext,".bgri",5) == 0) {
      call(mg_create_mesh(&Mesh));
      call(mg_read_bgri_file(InFile, Mesh));
    }
    else
      return error(err_NOT_SUPPORTED);
  }
  else {//read geometry file
    call(mg_get_input_char("GeometryFile", &InFile));
    call(mg_read_geo(&Geo, InFile));
    //mesh boundary and create initial front
    mg_Metric *Metric;
    int nNodeInSeg[5]={10,12,9,7,55};
    Metric = malloc(sizeof(mg_Metric));
    Metric->type = mge_Metric_Uniform;
    Metric->order = 8;
    //  Metric.type = mge_Metric_Uniform;
    //  Metric.order = 1;
    call(mg_create_mesh(&Metric->BGMesh));
    Metric->BGMesh->Dim = 2;
    call(mg_create_mesh(&Mesh));
    call(mg_create_bmesh_from_geo(Geo, Metric, nNodeInSeg, Mesh, &Front));
    mg_destroy_mesh(Metric->BGMesh);
    mg_free((void*)Metric);
  }
  
  //fill in face information
  call(mg_calc_face_info(Mesh));
  //fill in connectivities
  call(mg_build_connectivity(Mesh));
  //call(mg_prealloc_msh_comp(Mesh, 3, 15, 20));
  //create front
  call(mg_create_front(Mesh, &Front));
  call(mg_mesh_2_matlab(Mesh, &Front, "mesh_initial.m"));
  i = 0;
  while (!mg_front_empty(&Front)){
    printf("it = %d nElem = %d\n",i,Mesh->nElem);
//    sprintf(MeshName, "mesh_at_%d.m",i);
//    call(mg_mesh_2_matlab(Mesh, &Front, MeshName));
    //advance front
    ierr=error(mg_advance_front(Mesh, &Front));
    if (ierr != err_OK) {
      call(mg_mesh_2_matlab(Mesh, &Front,"mesh_error.m"));
      exit(0);
    }
    //DEBUGGING
    //call(xf_VerifyFront2D(&Front));
    i++;
  }
  //call(mg_plot_mesh(Mesh));
  call(mg_get_input_char("OutputMesh", &OutFile));
  call(mg_write_mesh(Mesh, OutFile));
  call(mg_mesh_2_matlab(Mesh, &Front, "mesh_final.m"));
  printf("Number of triangles: %d\nDone.\n",Mesh->nElem);
  
  mg_destroy_mesh(Mesh);
  //destroy hash table
  hdestroy();
  
  return err_OK;
}