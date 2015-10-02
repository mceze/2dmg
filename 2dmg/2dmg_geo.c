//
//  2dmg_geo.c
//  2dmg
//
//  Created by Marco Ceze on 6/10/15.
//  https://github.com/mceze/2dmg
//

#include "2dmg_geo.h"

/******************************************************************/
/* function:  mg_create_geo */
/* creates a mg_Geometry structure with "nBoundary" boundaries  */
int mg_create_geo(mg_Geometry **pGeo, int nBoundary, int nPoint,
                  int Dim)
{
  int ierr, i;
  mg_Segment *Boundary;
  
  call(mg_alloc((void**)pGeo, 1, sizeof(mg_Geometry)));
  (*pGeo)->nPoint = nPoint;
  (*pGeo)->Dim = Dim;
  call(mg_alloc((void**)&((*pGeo)->Coord), nPoint*Dim,
                sizeof(double)));
  (*pGeo)->nBoundary = nBoundary;
  call(mg_alloc((void**)&((*pGeo)->Boundary), nBoundary,
                sizeof(mg_Segment)));
  //init boundaries
  for (i = 0; i < nBoundary; i++) {
    call(mg_alloc((void**)&Boundary, 1,sizeof(mg_Segment)));
    (*pGeo)->Boundary[i] = Boundary;
    (*pGeo)->Boundary[i]->nPoint = 0;
    (*pGeo)->Boundary[i]->Point  = NULL;
    (*pGeo)->Boundary[i]->Coord  = NULL;
    (*pGeo)->Boundary[i]->s  = NULL;
    (*pGeo)->Boundary[i]->interp = NULL;
    (*pGeo)->Boundary[i]->accel  = NULL;
    (*pGeo)->Boundary[i]->Mij = NULL;
    (*pGeo)->Boundary[i]->Mij_interp = NULL;
    (*pGeo)->Boundary[i]->Mij_accel  = NULL;
    (*pGeo)->Boundary[i]->interp_type = -1;//unitialized
    call(mg_alloc((void**)&((*pGeo)->Boundary[i]->Name),MAXSTRLEN,sizeof(char)));
    (*pGeo)->Boundary[i]->Length = -1.0;//uninitialized
    (*pGeo)->Boundary[i]->Lm = -1.0;//uninitialized
  }
  
  return err_OK;
}

/******************************************************************/
/* function:  mg_destroy_geo */
/* destroys a mg_Geometry structure with "nBoundary" boundaries  */
void mg_destroy_geo(mg_Geometry *Geo)
{
  int i, d;
  mg_free((void*)Geo->Coord);
  for (i = 0; i < Geo->nBoundary; i++) {
    mg_free((void*)Geo->Boundary[i]->Coord);
    mg_free((void*)Geo->Boundary[i]->s);
    mg_free((void*)Geo->Boundary[i]->Name);
    if (Geo->Boundary[i]->interp != NULL){
      for (d = 0; d < Geo->Dim;d++){
        gsl_interp_free(Geo->Boundary[i]->interp[d]);
        gsl_interp_accel_free(Geo->Boundary[i]->accel[d]);
      }
      mg_free((void*)Geo->Boundary[i]->Point);
      mg_free((void*)Geo->Boundary[i]->interp);
      mg_free((void*)Geo->Boundary[i]->accel);
    }
    mg_free((void*)Geo->Boundary[i]);
  }
  mg_free((void*)Geo->Boundary);
  mg_free((void*)Geo);
}

/******************************************************************/
/* function:  mg_init_segment */
/* initializes the interpolant of a mg_Segment  */
int mg_init_segment(mg_Geometry *Geo, int iseg)
{
  int i;
  mg_Segment *Seg = Geo->Boundary[iseg];
  
  switch (Seg->interp_type) {
    case mge_Linear:
      for (i = 0; i < Geo->Dim; i++) {
        Seg->interp[i] = gsl_interp_alloc(gsl_interp_linear,Seg->nPoint);
        gsl_interp_init(Seg->interp[i], Seg->s,Seg->Coord
                        +i*Seg->nPoint, Seg->nPoint);
        Seg->accel[i] = gsl_interp_accel_alloc();
      }
      break;
    case mge_CSpline:
      for (i = 0; i < Geo->Dim; i++) {
        Seg->interp[i] = gsl_interp_alloc(gsl_interp_cspline,Seg->nPoint);
        gsl_interp_init(Seg->interp[i], Seg->s,Seg->Coord
                        +(i*Seg->nPoint), Seg->nPoint);
        Seg->accel[i] = gsl_interp_accel_alloc();
      }
      break;
    case mge_CSplinePeriodic:
      for (i = 0; i < Geo->Dim; i++) {
        Seg->interp[i] = gsl_interp_alloc(gsl_interp_cspline_periodic,Seg->nPoint);
        gsl_interp_init(Seg->interp[i], Seg->s,Seg->Coord
                        +(i*Seg->nPoint), Seg->nPoint);
        Seg->accel[i] = gsl_interp_accel_alloc();
      }
      break;
    case mge_Polynomial:
      for (i = 0; i < Geo->Dim; i++) {
        Seg->interp[i] = gsl_interp_alloc(gsl_interp_polynomial,Seg->nPoint);
        gsl_interp_init(Seg->interp[i], Seg->s,Seg->Coord
                        +(i*Seg->nPoint), Seg->nPoint);
        Seg->accel[i] = gsl_interp_accel_alloc();
      }
      break;
    default:
      return error(err_NOT_SUPPORTED);
      break;
  }
  
  
  return err_OK;
}

/******************************************************************/
/* function:  mg_create_bmesh_from_geo */
/* initializes the interpolant of a mg_Segment  */
int mg_create_bmesh_from_geo(mg_Geometry *Geo, mg_Metric *Metric,
                             int *nNodeInSeg, mg_Mesh *Mesh,
                             mg_Front *Front)
{
  int ierr, iseg, inode, in, d, iface, elem, nrmd;
  int nose_node, tail_node, i, n, ref_node_id;
  mg_List *SegList;
  mg_Segment *Seg;
  mg_FaceData *Face;
  struct mg_Item *seg_root, *seg_curr;
  double scale, *t;
  
  call(mg_alloc((void**)&SegList, 1, sizeof(mg_List)));
  //initialize list with contiguous set of segments
  SegList->nItem = Geo->nBoundary;
  call(mg_alloc((void**)&SegList->Item,SegList->nItem,sizeof(int)));
  for (iseg = 0; iseg < Geo->nBoundary; iseg++)
    SegList->Item[iseg] = iseg;
  
  //create space for nodes in mesh
  //assume empty mesh
  Mesh->nNode = 0;
  Mesh->Dim = Geo->Dim;
  
  //initialize face counter
  Mesh->nFace = 0;
  Mesh->nBfg  = Geo->nBoundary;
  //array of number of boundary faces
  call(mg_alloc((void**)&Mesh->nBface, Mesh->nBfg, sizeof(int)));
  //array of boudary group names
  call(mg_alloc2((void***)&Mesh->BNames, Mesh->nBfg, MAXSTRLEN, sizeof(char)));
  
  inode = 0;
  iface = 0;
  //while list of segments is not empty
  while (SegList->nItem > 0) {
    //always get first in list
    iseg = SegList->Item[0];
    Seg = Geo->Boundary[iseg];
    memcpy(Mesh->BNames[iseg], Seg->Name,strlen(Seg->Name)*sizeof(char));
    //check if segment is a closed loop
    if (Seg->Point[0] == Seg->Point[Seg->nPoint-1]) {
      //segment closes on itself
      call(mg_mesh_segment(Seg,Metric,nNodeInSeg[iseg],&scale,&t));
      nNodeInSeg[iseg]--;//subtract last repeated node
      //update node list
      Mesh->nNode += nNodeInSeg[iseg];
      call(mg_realloc((void**)&Mesh->Coord, Mesh->nNode*Mesh->Dim,
                      sizeof(double)));
      //allocate space for faces
      call(mg_realloc((void**)&Mesh->Face, Mesh->nFace+nNodeInSeg[iseg],
                      sizeof(mg_FaceData)));
      //get global coordinates and create loop. 
      for (in = 0; in < nNodeInSeg[iseg]; in++) {
        for (d = 0; d < Mesh->Dim; d++) {
          Mesh->Coord[(inode+in)*Mesh->Dim+d] = gsl_interp_eval(Seg->interp[d],Seg->s,
                                                                Seg->Coord+d*Seg->nPoint,t[in],
                                                                Seg->accel[d]);
        }
        call(mg_alloc((void**)&Face, 1, sizeof(mg_FaceData)));
        mg_init_face(Face);
        Face->nNode = 2;
        call(mg_alloc((void**)&Face->node, Face->nNode, sizeof(int)));
        Face->node[0] = inode+in;
        if (in < nNodeInSeg[iseg]-1)
          Face->node[1] = inode+in+1;
        else
          Face->node[1] = inode;//close last face with first node
        Face->elem[LEFTNEIGHINDEX] = HOLLOWNEIGHTAG;
        //encode face number and boundary group into one integer
        call(mg_limited_pair(in, iseg, &elem, Mesh->nBfg));
        Face->elem[RIGHTNEIGHINDEX] = -(elem+1);
        Mesh->Face[iface+in] = Face;
      }
      //these are equal because the this segment closes on itself
      inode += nNodeInSeg[iseg];
      iface += nNodeInSeg[iseg];
      Mesh->nBface[iseg] = nNodeInSeg[iseg];
      Mesh->nFace += Mesh->nBface[iseg];
      mg_free((void*)t);
      //remove segment from list
      call(mg_rm_frm_ord_set(iseg, &SegList->nItem, &SegList->Item, 1, &nrmd));
      if (nrmd != 1) return error(err_LOGIC_ERROR);
    }
    else {
      //start stitching segments
      seg_root = (struct mg_Item*)malloc(sizeof(struct mg_Item));
      seg_curr  = (struct mg_Item*)malloc(sizeof(struct mg_Item));
      seg_root->Id = iseg;
      seg_root->next = seg_curr;
      seg_curr->Id   = -1;
      seg_curr->next = NULL;
      //initial geometry node id
      nose_node = Geo->Boundary[seg_root->Id]->Point[0];
      n = Geo->Boundary[seg_root->Id]->nPoint;
      tail_node = Geo->Boundary[seg_root->Id]->Point[n-1];
      //remove root segment from list and pick next
      call(mg_rm_frm_ord_set(iseg, &SegList->nItem, &SegList->Item, 1, &nrmd));
      if (nrmd != 1) return error(err_LOGIC_ERROR);
      //loop around geometry until get back to root
      while (tail_node != nose_node) {
        for (i = 0; i < SegList->nItem; i++) {
          iseg = SegList->Item[i];
          if (Geo->Boundary[iseg]->Point[0] == tail_node) {
            //set new tail_node
            n = Geo->Boundary[iseg]->nPoint;
            tail_node = Geo->Boundary[iseg]->Point[n-1];
            seg_curr->Id = iseg;
            seg_curr->next  = (struct mg_Item*)malloc(sizeof(struct mg_Item));
            seg_curr = seg_curr->next;
            seg_curr->Id = -1;
            seg_curr->next = NULL;
            //remove newly found segment from list
            call(mg_rm_frm_ord_set(iseg, &SegList->nItem, &SegList->Item, 1, &nrmd));
            if (nrmd != 1) return error(err_LOGIC_ERROR);
            break;//restart loop over list
          }
        }
      }
      ref_node_id = inode;
      //now mesh each segment and discard appropriate nodes.
      seg_curr = seg_root;
      while (seg_curr->next != NULL){
        iseg = seg_curr->Id;
        Seg = Geo->Boundary[iseg];
        call(mg_mesh_segment(Seg,Metric,nNodeInSeg[iseg],&scale,&t));
        //allocate space for new nodes
        nNodeInSeg[iseg]--;
        Mesh->nNode += nNodeInSeg[iseg];
        call(mg_realloc((void**)&Mesh->Coord, Mesh->nNode*Mesh->Dim,
                        sizeof(double)));
        //allocate space for faces
        call(mg_realloc((void**)&Mesh->Face, Mesh->nFace+nNodeInSeg[iseg],
                        sizeof(mg_FaceData)));
        nose_node = Mesh->nNode;//next node to be created
        //get global coordinates and create faces.
        //discard last node as it is the last node of previous segment
//        printf("%s:\n",Seg->Name);
        for (in = 0; in < nNodeInSeg[iseg]; in++) {
          for (d = 0; d < Mesh->Dim; d++) {
            Mesh->Coord[(inode+in)*Mesh->Dim+d] = gsl_interp_eval(Seg->interp[d],Seg->s,
                                                                  Seg->Coord+d*Seg->nPoint,t[in],
                                                                  Seg->accel[d]);
//            printf("%1.8f ",Mesh->Coord[(inode+in)*Mesh->Dim+d]);
          }
//          printf("\n");
          call(mg_alloc((void**)&Face, 1, sizeof(mg_FaceData)));
          mg_init_face(Face);
          Face->nNode = 2;
          call(mg_alloc((void**)&Face->node, Face->nNode, sizeof(int)));
          Face->node[0] = inode+in;
          Face->node[1] = inode+in+1;
          Face->elem[LEFTNEIGHINDEX] = HOLLOWNEIGHTAG;
          //encode face number and boundary group into one integer
          call(mg_limited_pair(in, iseg, &elem, Mesh->nBfg));
          Face->elem[RIGHTNEIGHINDEX] = -(elem+1);
          Mesh->Face[iface+in] = Face;
        }
        //these are equal because the this segment closes on itself
        inode += nNodeInSeg[iseg];
        iface += nNodeInSeg[iseg];
        Mesh->nBface[iseg] = nNodeInSeg[iseg];
        Mesh->nFace += Mesh->nBface[iseg];
        mg_free((void*)t);
        seg_curr = seg_curr->next;
      }
      //fix last node id
      Face->node[1] = ref_node_id;
      
      mg_free_linked_list(seg_root);
    }
  }
  
  
  mg_free((void*)SegList->Item);
  
  return err_OK;
}
