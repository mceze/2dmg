//
//  2dmg_io.c
//  2dmg
//
//  Created by Marco Ceze on 11/11/14.
//  Copyright (c) 2014 Marco Ceze. All rights reserved.
//

#define _GNU_SOURCE
#include <search.h>
#include "2dmg_def.h"
#include "2dmg_struct.h"
#include "2dmg_utils.h"
#include "2dmg_math.h"
#include "2dmg_io.h"

/******************************************************************/
/* function:  mg_scan_n_num */
/* scan "n" numbers from a string  */
int mg_scan_n_num( const char *line, int *n, int *vi, double *vr)
{
  int ierr, pos, i;
  int ndelim, iv=0;
  const char delim[] = " \t";
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
  (*n) = 0;
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
        ierr = sscanf(value, "%d", vi+(*n));
      else
        ierr = sscanf(value, "%lf", vr+(*n));
      if (ierr != 1) return error(err_READWRITE_ERROR);
      (*n)++;
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
      ierr = sscanf(value, "%d", vi+(*n));
    else
      ierr = sscanf(value, "%lf", vr+(*n));
    if (ierr != 1) return error(err_READWRITE_ERROR);
    (*n)++;
  }
  
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

