#include <assert.h>
#include <stdint.h>

#include <gkyl_alloc.h>
#include <gkyl_rect_grid.h>
#include <gkyl_rect_grid_priv.h>

void
gkyl_rect_grid_init(struct gkyl_rect_grid *grid, int ndim,
  const double *lower, const double *upper, const int *cells)
{
//  // MF 2023/07/07: commenting this out because it causes seg faults in g2.
//  *grid = (struct gkyl_rect_grid) { };
  
  grid->ndim = ndim;  
  grid->cellVolume = 1.0;
  for (int i=0; i<ndim; ++i) {
    grid->lower[i] = lower[i];
    grid->upper[i] = upper[i];
    grid->cells[i] = cells[i];
    grid->dx[i] = (upper[i]-lower[i])/cells[i];
    grid->cellVolume *= grid->dx[i];
  }
}

int* gkyl_find_cell(struct gkyl_rect_grid *grid, const double *point, bool pickLower, const int **knownIdx){
  int *cellIdx;
  int nDim = grid->ndim;
  int searchNum = 0;
  int *searchDim, *idxOut;
  int *dimTrans[GKYL_MAX_DIM];
  bool allLessEq = true;
  int *iStart, *iEnd, *iMid, *iNew, *cells;
  double lowerInDir[nDim], upperInDir[nDim];
  
  printf("Start gkyl_find_cell\n");
  cellIdx = (int*)malloc(nDim*sizeof(int));
  searchDim = (int*)malloc(nDim*sizeof(int));
  //lowerInDir = (int*)malloc(nDim*sizeof(int));
  //upperInDir = (int*)malloc(nDim*sizeof(int));
  for(int d=0; d<nDim; d++){
    dimTrans[d] = (int*)malloc(1*sizeof(int));
    if(knownIdx[d]==NULL){
      searchDim[searchNum] = d;
      *dimTrans[d] = searchNum;
      searchNum = searchNum + 1;
    }else{
      dimTrans[d] = NULL;      
      printf("In else, d=%i, dimTrans[d]=%i, test:%i\n", d, dimTrans[d], NULL);
      cellIdx[d] = *knownIdx[d];
    }
    if(dimTrans[d]==NULL){
      printf("searchDim[d]: %i, dimTrans[d] is NULL, searchNum: %i\n",searchDim[searchNum-1]+1, searchNum);
    }else{
      printf("searchDim[d]: %i, dimTrans[d]: %i, searchNum: %i\n",searchDim[searchNum-1]+1, *dimTrans[d]+1, searchNum);
    }
  }

  iStart = (int*)malloc(searchNum*sizeof(int));
  iEnd = (int*)malloc(searchNum*sizeof(int));
  iMid = (int*)malloc(searchNum*sizeof(int));
  iNew = (int*)malloc(searchNum*sizeof(int));
  cells = grid -> cells;
  for(int d=0; d<searchNum; d++){
    iStart[d] = 1;
    iEnd[d] = cells[searchDim[d]];
    iMid[d] = 0;
    iNew[d] = 0;
  }
  /* Below we use a binary search. That is, if the i-th coordinate of the point in
   * question is below(above) the ith-coordinate of the current mid-point, we search
   * the lower(upper) half along that direction in the next iteration.
   */
  while (allLessEq){
    for(int d=0; d<searchNum; d++){
      iMid[d] = (iStart[d]+iEnd[d])/2;//Integer division intentional
      printf("dim=%i\n",d);
      printf(" Mid[d]=%i,",iMid[d]);
      printf(" Start[d]=%i,",iStart[d]);
      printf(" End[d]=%i,\n",iEnd[d]);
    }

    if(isInCell(grid, point, iMid, dimTrans, knownIdx)){
      //Check if neighboring cells also contain this point.
      //TO DO...

      for(int d=0; d<searchNum; d++){
	cellIdx[searchDim[d]] = iMid[d];
        //cellIdx[d] = knownIdx[d]==NULL ? idxOut[dimTrans[d]] : knownIdx[d];
	//printf("d=%i,mid=%i,idx=%i\n",d,iMid[d],cellIdx[d]);
      }
      break;
    }else{
      InDir(grid, iMid, dimTrans, knownIdx, lowerInDir, upperInDir);
      for(int d=0; d<searchNum; d++){
	if(point[searchDim[d]] < lowerInDir[searchDim[d]]){
	  iEnd[d] = iMid[d]-1;
	}else if(point[searchDim[d]] > upperInDir[searchDim[d]]){
	  iStart[d] = iMid[d]+1;
	}
      }
    }

    for(int d=0; d<searchNum; d++){
      if(iStart[d]>iEnd[d]){
	allLessEq = false;
	break;
      }
    }
  }

  printf("End gkyl_find_cell\n");
  return cellIdx;
}

bool isInCell(const struct gkyl_rect_grid *grid, const double *pIn, int *iIn, int *dimTrans[1], const int **knownIdx){
  double eps = 1.0e-14;
  int checkIdx;
  bool inCell = true;
  int nDim;
  nDim = grid -> ndim;
  double lowerInDir[nDim], upperInDir[nDim];
  //lowerInDir = (int*)malloc(nDim*sizeof(int));
  //upperInDir = (int*)malloc(nDim*sizeof(int));
  for(int d=0;d<nDim; d++){
    lowerInDir[d]=0;
    upperInDir[d]=0;
  }
  InDir(grid, iIn, dimTrans, knownIdx, lowerInDir, upperInDir);
  for(int d=0;d<nDim; d++){
    printf("lowerInDir[%i] = %f \n",d,lowerInDir[d]);
  }
  for(int d=0; d<nDim; d++){
    printf("d=%i, point[%i] = %f\n",d,d,pIn[d]);
    printf("lowerInDir[%i] = %f ",d,lowerInDir[d]);
    printf("upperInDir[%i] = %f\n",d,upperInDir[d]);
    if(lowerInDir[d]-eps>pIn[d] || upperInDir[d]+eps<pIn[d]){
      inCell = false;
      break;
    }
  }
  return inCell;
}

void InDir(const struct gkyl_rect_grid *grid, int *iIn, int *dimTrans[1], const int **knownIdx, double lowerInDir[], double upperInDir[]){
  int nDim, checkIdx;
  const double *dx, *lower;
  nDim = grid -> ndim;
  double testarr[nDim];
  setbuf(stdout, NULL);
  dx = grid -> dx;
  lower = grid -> lower;
  for(int d=0; d<nDim; d++){
    checkIdx = knownIdx[d]==NULL ? iIn[*dimTrans[d]] : *knownIdx[d];
    printf("lowerInDir before calc: %f\n",lowerInDir[d]);
    lowerInDir[d] = lower[d]+(checkIdx-1)*dx[d];
    testarr[d] = lower[d]+(checkIdx-1)*dx[d];
    upperInDir[d] = lower[d]+(checkIdx)*dx[d];
    printf("d=%i, checkIDx= %i\n",d,checkIdx);
    printf("lower[%i]=%f, dx[%i]=%f\n",d,lower[d],d,dx[d]);
    printf("lowerInDir[d]=%f, calc:%f upperInDir[d]=%f\n", lowerInDir[d],lower[d]+(checkIdx-1)*dx[d],upperInDir[d]);
    printf("all lowerindir %f,%f,%f,%f\n",lowerInDir[0],lowerInDir[1],lowerInDir[2],lowerInDir[3]);
    printf("all testarr %f,%f,%f,%f\n",testarr[0],testarr[1],testarr[2],testarr[3]);
  }
  for(int d=0;d<nDim; d++){
    printf("lowerInDir[%i] = %f ",d,lowerInDir[d]);
    printf("testarr[%i] = %f \n",d,testarr[d]);
  }
}

void
gkyl_rect_grid_write(const struct gkyl_rect_grid *grid, FILE *fp)
{
  // dimension and shape are written as 64 bit integers
  uint64_t ndim = grid->ndim;
  uint64_t cells[GKYL_MAX_DIM];
  for (int d=0; d<grid->ndim; ++d)
    cells[d] = grid->cells[d];
  
  fwrite(&ndim, sizeof(uint64_t), 1, fp);
  fwrite(cells, sizeof(uint64_t), grid->ndim, fp);
  fwrite(grid->lower, sizeof(double), grid->ndim, fp);
  fwrite(grid->upper, sizeof(double), grid->ndim, fp);
}

bool
gkyl_rect_grid_read(struct gkyl_rect_grid *grid, FILE *fp)
{
  uint64_t ndim = grid->ndim;
  uint64_t cells64[GKYL_MAX_DIM];

  if (1 != fread(&ndim, sizeof(uint64_t), 1, fp))
    return false;
  if (ndim != fread(cells64, sizeof(uint64_t), ndim, fp))
    return false;

  double lower[GKYL_MAX_DIM], upper[GKYL_MAX_DIM];
  if (ndim != fread(lower, sizeof(double), ndim, fp))
    return false;
  if (ndim != fread(upper, sizeof(double), ndim, fp))
    return false;

  // copy into regular int array
  int cells[GKYL_MAX_DIM];
  for (int d=0; d<ndim; ++d) cells[d] = cells64[d];

  gkyl_rect_grid_init(grid, ndim, lower, upper, cells);

  return true;
}
