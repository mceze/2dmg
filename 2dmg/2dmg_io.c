//
//  2dmg_io.c
//  2dmg
//
//  Created by Marco Ceze on 11/11/14.
//  https://github.com/mceze/2dmg
//

#define _GNU_SOURCE
#include <search.h>
#include "2dmg_def.h"
#include "2dmg_struct.h"
#include "2dmg_utils.h"
#include "2dmg_math.h"
#include "2dmg_geo.h"
#include "2dmg_io.h"

/******************************************************************/
/* function:  mg_scan_n_num */
/* scan "n" numbers from a string  */
int mg_scan_n_num( const char *line, int *n, int *vi, double *vr)
{
  int ierr, pos, i, k;
  int ndelim, iv=0;
  const char delim[] = " ";
  char value[MAXSTRLEN];
  bool indelim, indelimcur;
  
  /* return 0 if null */
  if (line == NULL){
    (*n) = 0;
    return err_OK;
  }
  
  /* Error if vi, vr both null or both not null */
  if ( ((vi == NULL) && (vr == NULL)) ||  ((vi != NULL) && (vr != NULL)))
    return error(err_INPUT_ERROR);
  
  /* Get number of delim characters */
  ndelim = (int)strlen(delim);
  
  /* Loop until hit start of token # ntok */
  k = 0;
  pos = 0;
  indelim = true;
  while ((line[pos] != '\0') && (line[pos] != '\n')){
    indelimcur = false;
    for (i=0; i<ndelim; i++)
      if (line[pos] == delim[i]){
        indelimcur = true;
        break;
      }
    
    if ((!indelimcur) && (indelim)){
      // hit start of token
      indelim = false;
      iv = 0;
    }
    if ((indelimcur) && (!indelim)){
      // hit end of token
      value[iv] = '\0';
      if (vi != NULL)
        ierr = sscanf(value, "%d", vi+k);
      else
        ierr = sscanf(value, "%lf", vr+k);
      if (ierr != 1) return error(err_READWRITE_ERROR);
      k++;
    }
    indelim = indelimcur;
    if (!indelim){
      value[iv] = line[pos];
      iv++;
    }
    pos++;
  }
  
  if (!indelim){
    value[iv] = '\0';
    if (vi != NULL)
      ierr = sscanf(value, "%d", vi+k);
    else
      ierr = sscanf(value, "%lf", vr+k);
    if (ierr != 1) return error(err_READWRITE_ERROR);
    k++;
  }
  (*n) = k;
  
  return err_OK;
}

/******************************************************************/
/* function:  mg_read_bgri_file */
/* reads a boundary discretization file  */
int mg_read_bgri_file(const char *InFile, mg_Mesh *Mesh)
{
  int ierr, nnode, nface, dim, i, n, nbfg, k, nk, nj, e;
  int node[3], j;
  char line[MAXLONGLINELEN];
  FILE *bgri = NULL;
  mg_FaceData *Face;
  
  printf("Reading %s file:\n",InFile);
  //open boundary mesh file
  bgri = fopen(InFile, "r");
  if (!bgri) return error(err_READWRITE_ERROR);
  rewind(bgri);
  
  if (fscanf(bgri, "%d %d %d\n",&nnode,&nface,&dim) != 3)
    return error(err_READWRITE_ERROR);
  if (dim != 2) return error(err_NOT_SUPPORTED);
  
  Mesh->nNode = nnode;
  Mesh->Dim = dim;
  //allocate and read coordinates
  call(mg_alloc((void**)&Mesh->Coord, Mesh->nNode*Mesh->Dim, sizeof(double)));
  for (i = 0; i < Mesh->nNode; i++) {
    if (fgets(line, MAXLONGLINELEN, bgri) == NULL)
      return error(err_READWRITE_ERROR);
    call(mg_scan_n_num(line, &k, NULL, Mesh->Coord+i*dim));
    if (k != dim) return error(err_READWRITE_ERROR);
  }
  //get number of boundary groups
  if (fgets(line, MAXLONGLINELEN, bgri) == NULL)
    return error(err_READWRITE_ERROR);
  if (sscanf(line, "%d",&nbfg) != 1)
    return error(err_READWRITE_ERROR);
  Mesh->nBfg = nbfg;
  call(mg_alloc2((void***)&Mesh->BNames, nbfg, MAXSTRLEN, sizeof(char)));
  //allocate and read boundaries
  call(mg_alloc((void**)&Mesh->Face, nface, sizeof(mg_FaceData)));
  call(mg_alloc((void**)&Mesh->nBface, nbfg, sizeof(int)));
  Mesh->nFace = nface;
  n = 0;
  for (i = 0; i < nbfg; i++) {
    if (fgets(line, MAXLONGLINELEN, bgri) == NULL)
      return error(err_READWRITE_ERROR);
    if (sscanf(line, "%d %d %s",&nk,&nj,Mesh->BNames[i]) != 3)
      return error(err_READWRITE_ERROR);
    Mesh->nBface[i] = nk;
    if (nj != 2) return error(err_NOT_SUPPORTED);//for now
    printf("Read BGroup: %s\n",Mesh->BNames[i]);
    for (k = 0; k < nk; k++) {
      if (fgets(line, MAXLONGLINELEN, bgri) == NULL)
        return error(err_READWRITE_ERROR);
      call(mg_alloc((void**)&Face, 1, sizeof(mg_FaceData)));
      mg_init_face(Face);
      //nodes
      call(mg_alloc((void**)&Face->node, Mesh->Dim, sizeof(int)));
      call(mg_scan_n_num(line, &j, node, NULL));
      if (j != nj)return error(err_READWRITE_ERROR);
      Face->nNode = nj;
      //convert to 0-based numbering
      for (j=0;j<nj;j++) Face->node[j] = node[j]-1;
      //left element will be the interior
      Face->elem[LEFTNEIGHINDEX] = HOLLOWNEIGHTAG;
      //encode face number and boundary group into one integer
      call(mg_limited_pair(k, i, &e, Mesh->nBfg));
      //right side is the boundary
      Face->elem[RIGHTNEIGHINDEX] = -(e+1);
      Mesh->Face[n] = Face;
      n++;
    }
  }
  if (n != nface){
    printf("File indicates nface = %d but I only read %d.\n",nface,n);
    return error(err_READWRITE_ERROR);
  }
  
  fclose(bgri);
  
  return err_OK;
}

/******************************************************************/
/* function:  mg_parse_input_line */
int mg_parse_input_line(char line[MAXLINELEN], char **pkey,
                        char **pvalue)
{
  int ierr, i, len;
  char *ptr;
  bool delim_found = false;
  len = (int)strlen(line);
  
  for (i = 0; i < len; i++) {
    if (line[i] == '='){
      delim_found = true;
      break;
    }
  }
  if (i < 1 || !delim_found)
    return error(err_READWRITE_ERROR);
  call(mg_alloc((void**)pkey, MAXSTRLEN, sizeof(char)));
  call(mg_alloc((void**)pvalue, MAXSTRLEN, sizeof(char)));
  line[i] = '\0';
  ptr = (char*)line+i+1;
  sscanf(line, "%s",(*pkey));
  sscanf(ptr, "%s",(*pvalue));
  
  return err_OK;
}


/******************************************************************/
/* function:  mg_read_input_file */
int mg_read_input_file(char const ParamFile[])
{
  int ierr, ikey;
  char line[MAXLINELEN], *value;
  ENTRY item, *entry;
  FILE *fid;
  
  if ((fid = fopen(ParamFile, "r")) == NULL)
    return error(err_READWRITE_ERROR);
  rewind(fid);
  
  ierr = hcreate(NPARAMLIST);
  if (ierr == 0) return error(err_MEMORY_ERROR);
  ikey = 0;
  while (!feof(fid)) {
    fgets(line, MAXSTRLEN, fid);
    //check if line is a comment
    if (strncmp(line, "#",1) == 0 ||
        strncmp(line, "%",1) == 0 ||
        strncmp(line, "//",2) == 0 ||
        line[0] == '\n')
      continue;
    //parse key and value
    call(mg_parse_input_line(line, &item.key, &value));
    //add entry to hash
    item.data = value;
    if ((entry = hsearch(item, ENTER)) == NULL)
      return error(err_HSEARCH_ERROR);
    
    ikey++;
  }
  
  return err_OK;
}

/******************************************************************/
/* function:  mg_read_input_file */
int mg_get_input_char(char const ParamName[], char **pvalue)
{
  ENTRY *e, target;
  
  target.key = malloc(MAXSTRLEN*sizeof(char));
  sprintf(target.key, "%s",ParamName);
  if ((e = hsearch(target, FIND)) == NULL)
    return error(err_HSEARCH_ERROR);
  (*pvalue) = (char*)e->data;
  free(target.key);
  
  return err_OK;
}

/******************************************************************/
/* function: mg_mesh_2_matlab */
/* converts mesh to matlab format */
int mg_mesh_2_matlab(mg_Mesh *Mesh, mg_Front *Front, char *FileName)
{
  int nodeID, d, faceID, elemID, loopID;
  mg_FrontFace *FFace;
  FILE *fid;
  
  if ((fid=fopen(FileName ,"w"))==NULL)
    return error(err_READWRITE_ERROR);
  
  //print node array
  fprintf(fid,"coord=[");
  for (nodeID = 0; nodeID < Mesh->nNode; nodeID++){
    for (d = 0; d < Mesh->Dim; d++)
      fprintf(fid, "%1.8e ",Mesh->Coord[nodeID*Mesh->Dim+d]);
    fprintf(fid, "\n");
  }
  fprintf(fid, "];\n");
  
  //print face array
  fprintf(fid,"face=[");
  for (faceID = 0; faceID < Mesh->nFace; faceID++) {
    for (d = 0; d < Mesh->Face[faceID]->nNode; d++)
      fprintf(fid, "%d ",Mesh->Face[faceID]->node[d]+1);
    fprintf(fid, "\n");
  }
  fprintf(fid, "];\n");
  
  //print element array
  fprintf(fid,"elem=[");
  for (elemID = 0; elemID < Mesh->nElem; elemID++) {
    for (d = 0; d < Mesh->Elem[elemID].nNode; d++)
      fprintf(fid, "%d ",Mesh->Elem[elemID].node[d]+1);
    fprintf(fid, "\n");
  }
  fprintf(fid, "];\n");
  
  //print front
  fprintf(fid,"front=[");
  for (loopID = 0; loopID < Front->nloop; loopID++) {
    if (Front->loop[loopID]->FacesInLoop->nEntry == 0)
      continue;
    FFace = Front->loop[loopID]->head;
    if (FFace == NULL) continue;
    for (d = 0; d < FFace->face->nNode; d++)
      fprintf(fid, "%d ",FFace->face->node[d]+1);
    fprintf(fid, "\n");
    FFace = FFace->next;
    while (FFace != Front->loop[loopID]->head) {
      for (d = 0; d < FFace->face->nNode; d++)
        fprintf(fid, "%d ",FFace->face->node[d]+1);
      fprintf(fid, "\n");
      FFace = FFace->next;
    }
  }
  fprintf(fid, "];\n");
  
  fclose(fid);
  
  
  return err_OK;
}

/******************************************************************/
/* function: mg_write_mesh */
/* writes mesh to a file with connectivities and boundary information */
int mg_write_mesh(mg_Mesh *Mesh, char *FileName)
{
  int i, d;
  FILE *fid;
  
  if ((fid = fopen(FileName, "w")) == NULL)
    return error(err_READWRITE_ERROR);
  //write header
  fprintf(fid, "%% Dim nNode nFace nElem nBfg\n");
  fprintf(fid, "%d %d %d %d %d\n", Mesh->Dim, Mesh->nNode, Mesh->nFace, Mesh->nElem, Mesh->nBfg);
  //write nodal coordinates
  fprintf(fid, "%% Node coordinates\n");
  for (i = 0; i < Mesh->nNode; i++) {
    for (d = 0; d < Mesh->Dim; d++)
      fprintf(fid, "%1.12e ",Mesh->Coord[i*Mesh->Dim+d]);
    fprintf(fid, "\n");
  }
  //write boundary grou information
  fprintf(fid, "%% BGroup nBface\n");
  for (i = 0; i < Mesh->nBfg; i++) {
    fprintf(fid, "%s %d\n",Mesh->BNames[i],Mesh->nBface[i]);
  }
  //write element-to-node connectivity
  fprintf(fid, "%% Element to node connectivity\n");
  for (i = 0; i < Mesh->nElem; i++){
    for (d = 0; d < Mesh->Elem[i].nNode; d++)
      fprintf(fid, "%d ",Mesh->Elem[i].node[d]);
    fprintf(fid, "\n");
  }
  //Face connectivity
  fprintf(fid, "%% n0 n1 eL eR\n");
  for (i = 0; i < Mesh->nFace; i++) {
    fprintf(fid, "%d %d %d %d\n",Mesh->Face[i]->node[0],
            Mesh->Face[i]->node[1],Mesh->Face[i]->elem[LEFTNEIGHINDEX],
            Mesh->Face[i]->elem[RIGHTNEIGHINDEX]);
  }
  
  fclose(fid);
  
  
  return err_OK;
}

/******************************************************************/
/* function: mg_read_mesh */
/* reads mesh from file with connectivities and boundary information */
int mg_read_mesh(mg_Mesh **pMesh, char *FileName)
{
  int ierr, n, vi[5], i, j, elemL, elemR;
  char line[MAXLINELEN];
  mg_Mesh *Mesh;
  mg_FaceData *Face;
  FILE *fid = fopen(FileName, "r");
  
  call(mg_create_mesh(&Mesh));
  
  //get mesh info
  fgets(line, MAXLINELEN, fid);
  while (line[0] == '%') {//skip comments
    fgets(line, MAXLINELEN, fid);
  }
  call(mg_scan_n_num(line, &n, vi, NULL));
  if (n != 5) return error(err_READWRITE_ERROR);
  Mesh->Dim   = vi[0];
  Mesh->nNode = vi[1];
  Mesh->nFace = vi[2];
  Mesh->nElem = vi[3];
  Mesh->nBfg  = vi[4];
  call(mg_alloc((void **)&Mesh->Coord, Mesh->nNode*Mesh->Dim,
                sizeof(double)));
  call(mg_alloc((void**)&Mesh->Elem, Mesh->nElem, sizeof(mg_ElemData)));
  call(mg_alloc2((void ***)&Mesh->BNames, Mesh->nBfg, MAXSTRLEN,
                 sizeof(char)));
  call(mg_alloc((void **)&Mesh->nBface, Mesh->nBfg,
                sizeof(int)));
  
  //get coordinates
  fgets(line, MAXLINELEN, fid);
  while (line[0] == '%') {//skip comments
    fgets(line, MAXLINELEN, fid);
  }
  for (i = 0 ; i < Mesh->nNode; i++) {
    call(mg_scan_n_num(line, &n, NULL, Mesh->Coord+i*Mesh->Dim));
    fgets(line, MAXLINELEN, fid);
    if (n != Mesh->Dim) return error(err_READWRITE_ERROR);
  }
  //get boundary groups
  while (line[0] == '%') {//skip comments
    fgets(line, MAXLINELEN, fid);
  }
  for (i = 0; i < Mesh->nBfg; i++) {
    sscanf(line,"%s %d",Mesh->BNames[i],&Mesh->nBface[i]);
    fgets(line, MAXLINELEN, fid);
  }
  //element-to-node connectivity
  while (line[0] == '%') {//skip comments
    fgets(line, MAXLINELEN, fid);
  }
  for (i = 0; i < Mesh->nElem; i++) {
    call(mg_scan_n_num(line, &n, vi, NULL));
    //for now, complain if elem is not a triangle
    if (n != 3) return error(err_NOT_SUPPORTED);
    Mesh->Elem[i].nNode = n;
    call(mg_alloc((void**)&(Mesh->Elem[i].face), Mesh->Elem[i].nNode,
                  sizeof(int)));
    call(mg_alloc((void**)&(Mesh->Elem[i].nbor), Mesh->Elem[i].nNode,
                  sizeof(int)));
    call(mg_alloc((void**)&(Mesh->Elem[i].node), Mesh->Elem[i].nNode,
                  sizeof(int)));
    for (j = 0; j < Mesh->Elem[i].nNode; j++)
      Mesh->Elem[i].node[j] = vi[j];
    fgets(line, MAXLINELEN, fid);
  }
  //face information
  call(mg_alloc((void**)&Mesh->Face, Mesh->nFace, sizeof(mg_FaceData)));
  while (line[0] == '%') {//skip comments
    fgets(line, MAXLINELEN, fid);
  }
  for (i = 0; i < Mesh->nFace; i++) {
    call(mg_scan_n_num(line, &n, vi, NULL));
    //may have to change this check because of multinode faces in the future
    if (n != 4) return error(err_READWRITE_ERROR);
    call(mg_alloc((void**)&Face, 1, sizeof(mg_FaceData)));
    mg_init_face(Face);
    Mesh->Face[i] = Face;
    Face->nNode = n-2;//number of entries read minus element numbers
    call(mg_alloc((void**)&Face->node, Face->nNode, sizeof(int)));
    for (j = 0; j < Face->nNode; j++)
      Face->node[j] = vi[j];
    elemL = vi[n-2];
    elemR = vi[n-1];
    Face->elem[LEFTNEIGHINDEX] = elemL;
    Face->elem[RIGHTNEIGHINDEX] = elemR;
    Mesh->Face[i] = Face;
    //figure out which face this is
    //left element
    if (elemL >= 0){
      for (j = 0; j < Mesh->Elem[elemL].nNode; j++) {
        if (Mesh->Elem[elemL].node[j] == Face->node[0]) {
          n = nextincycle(j, Mesh->Elem[elemL].nNode);
          if (Face->node[Face->nNode-1] != Mesh->Elem[elemL].node[n])
            return error(err_LOGIC_ERROR);
          n = previncycle(j, Mesh->Elem[elemL].nNode);
          Mesh->Elem[elemL].face[n] = i;
          Mesh->Elem[elemL].nbor[n] = elemR;
        }
      }
    }
    //right element
    if (elemR >= 0){
      for (j = 0; j < Mesh->Elem[elemR].nNode; j++) {
        if (Mesh->Elem[elemR].node[j] == Face->node[0]) {
          n = previncycle(j, Mesh->Elem[elemR].nNode);
          if (Face->node[Face->nNode-1] != Mesh->Elem[elemR].node[n])
            return error(err_LOGIC_ERROR);
          n = nextincycle(j, Mesh->Elem[elemR].nNode);
          Mesh->Elem[elemR].face[n] = i;
          Mesh->Elem[elemR].nbor[n] = elemL;
        }
      }
    }
    fgets(line, MAXLINELEN, fid);
  }
  
  
  (*pMesh) = Mesh;
  fclose(fid);
  
  return err_OK;
}

/******************************************************************/
/* function: mg_read_mesh */
/* reads mesh from file with connectivities and boundary information */
int mg_read_geo(mg_Geometry **pGeo, char *FileName)
{
  int ierr, n, vi[5], i, nB, id, nP, d;
  double ds;
  char line[MAXLINELEN], type[MAXSTRLEN];
  FILE *fid = fopen(FileName, "r");
  
  //get geo info
  fgets(line, MAXLINELEN, fid);
  while (line[0] == '%') {//skip comments
    fgets(line, MAXLINELEN, fid);
  }
  call(mg_scan_n_num(line, &n, vi, NULL));
  if (n != 3) return error(err_READWRITE_ERROR);
  //create geo structure
  call(mg_create_geo(pGeo, vi[2], vi[1], vi[0]));
  //read coordinates
  fgets(line, MAXLINELEN, fid);
  while (line[0] == '%') {//skip comments
    fgets(line, MAXLINELEN, fid);
  }
  for (i = 0; i < (*pGeo)->nPoint; i++) {
    call(mg_scan_n_num(line, &n, NULL, (*pGeo)->Coord+i*(*pGeo)->Dim));
    fgets(line, MAXLINELEN, fid);
    while (line[0] == '%') {//skip comments
      fgets(line, MAXLINELEN, fid);
    }
  }
  //read segments
  for (nB = 0; nB < (*pGeo)->nBoundary; nB++) {
    while (line[0] == '%') {//skip comments
      fgets(line, MAXLINELEN, fid);
    }
    sscanf(line, "%s %s %d",(*pGeo)->Boundary[nB]->Name,type,&nP);
    (*pGeo)->Boundary[nB]->nPoint = nP;
    ds = 1.0/(nP-1);
    call(mg_value_2_enum(type, mge_GeoInterpName,(int)mge_GeoInterpLast,
                         (int*)&((*pGeo)->Boundary[nB]->interp_type)));
    //allocate and read points
    call(mg_alloc((void**)&((*pGeo)->Boundary[nB]->Coord),
                  (*pGeo)->Boundary[nB]->nPoint*(*pGeo)->Dim,
                  sizeof(double)));
    call(mg_alloc((void**)&((*pGeo)->Boundary[nB]->s),
                  (*pGeo)->Boundary[nB]->nPoint,sizeof(double)));
    call(mg_alloc((void**)&((*pGeo)->Boundary[nB]->Point),
                  (*pGeo)->Boundary[nB]->nPoint,sizeof(int)));
    call(mg_alloc((void**)&((*pGeo)->Boundary[nB]->interp),(*pGeo)->Dim,
                  sizeof(gsl_interp)));
    call(mg_alloc((void**)&((*pGeo)->Boundary[nB]->accel),(*pGeo)->Dim,
                  sizeof(gsl_interp_accel)));
    for (i = 0; i < nP; i++) {
      fscanf(fid, "%d\n",&id);
      id--;
      (*pGeo)->Boundary[nB]->Point[i] = id;
      for (d = 0; d < (*pGeo)->Dim; d++){
        (*pGeo)->Boundary[nB]->Coord[d*nP+i] = (*pGeo)->Coord[id*(*pGeo)->Dim+d];
      }
      (*pGeo)->Boundary[nB]->s[i] = i*ds;
    }
    //initialize interpolants
    call(mg_init_segment((*pGeo), nB));
    
    fgets(line, MAXLINELEN, fid);
  }
  
  fclose(fid);
  
  return err_OK;
}

