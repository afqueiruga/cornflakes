#ifndef __TARGET_H
#define __TARGET_H

#include "cfmat.h"
#include "cfdata.h"

/* A polymorphic class for targets with different backends */
typedef struct _TARGET_VTABLE_t _TARGET_VTABLE_t;
typedef struct target_t {
  int rank;
  cfdata_t *R;
  cfmat_t *K;
} target_t;

/* The helper constructor */
void Target_New(target_t * self, int onum,
		kernel_t * ke, hypergraph_t * hg, int ndof,
		char * backend);
void Target_New_From_Ptr(target_t * self, int rank,  void * payload);
/* The member methods */
real_t * Target_Place(target_t * self, int n, int * dofs, real_t * vals);
void Target_Destroy(target_t * self);
void Target_Wipe(target_t * self);
void Target_Finalize(target_t * self);


#if 0
/*
 * The constructor for the built in type, meant to be used with 
 * Scipy in python-cornflakes
 */
void Target_Default_New(target_t * self, int onum,
			kernel_t * ke, hypergraph_t * hg, int ndof);
void Target_Default_From_Array(target_t * self, int rank, int ndof,
			       real_t * V, int * II, int * JJ);
typedef struct Target_Default_data_t {
  int nalloc;
  int own;
  real_t * V;
  int * II;
  int * JJ;

  real_t * Viter;
  int * IIiter;
  int * JJiter;
} Target_Default_data_t;
#define Target_Default_Data(x) ((Target_Default_data_t*)((x)->data))
#endif
#endif
