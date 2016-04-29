#ifndef __TARGET_H
#define __TARGET_H

#include "kernel.h"
#include "hypergraph.h"

/* A polymorphic class for targets with different backends */
typedef struct _TARGET_VTABLE_t _TARGET_VTABLE_t;

typedef struct target_t {
  int rank;
  int N;
  void * data;
  const _TARGET_VTABLE_t * vtable;
} target_t;

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

/* The helper constructor */
void Target_New(target_t * self, int onum,
		kernel_t * ke, hypergraph_t * hg, int ndof,
		char * backend);

/* The member methods */
real_t * Target_Place(target_t * self, int n, int * dofs, real_t * vals);
void Target_Destroy(target_t * self);
void Target_Wipe(target_t * self);

#endif
