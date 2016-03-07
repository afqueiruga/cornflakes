#ifndef _SPATIAL_HASH_H
#define _SPATIAL_HASH_H

#define SPATIAL_HASH_MAXDIM 3

typedef double real_t;

typedef struct spatialhash_t {
  int dim;
  int Npart;

  real_t start[SPATIAL_HASH_MAXDIM];
  real_t h[SPATIAL_HASH_MAXDIM];
  int Ncell[SPATIAL_HASH_MAXDIM];
  
  int * cells;
  int * list;
} spatialhash_t;

void SpatialHash_init(spatialhash_t * sh, int Npart, int dim,
		      real_t *start, real_t * end,real_t *h);
void SpatialHash_destroy(spatialhash_t * sh);

void SpatialHash_Push(spatialhash_t * sh, int A, real_t * x);
void Build_New_Hash(spatialhash_t * sh, int npart,int dim, real_t * x, real_t cutoff);
void SpatialHash_Scanall(spatialhash_t * sh, real_t * x,void act(int,int) );
void SpatialHash_ScanPt(spatialhash_t * sh, real_t * x, void act(int,int) );

void SpatialHash_Scan_Area(spatialhash_t * sh, int * hs, int A, void act(int,int), int half );
void SpatialHash_Scan_Cell(spatialhash_t * sh, void act(int,int),
			   int a,
			   int * hs, int * off, int filter);





#endif
