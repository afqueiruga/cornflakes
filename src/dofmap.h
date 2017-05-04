#ifndef __H_DOFMAP
#define __H_DOFMAP

#include "hypergraph.h"

/* Vtable forward declaration*/
typedef struct _DOFMAP_VTABLE_t _DOFMAP_VTABLE_t;

/* Structure. Using a union for now... Not extensible*/
typedef struct Dofmap {
  //enum dofmap_type type;
  union {
    struct strided_t{
      int stride;
      int offset;
    }  strided ;
    struct tabled_t{
      int stride;
      int * table;
      int offset;
    } tabled;
  } U;
  const _DOFMAP_VTABLE_t * vtable;
} dofmap_t;

/* Constructors */
void Dofmap_Strided(dofmap_t * dm, int stride, int offset);
void Dofmap_Tabled(dofmap_t * dm, int Nentry, int stride, int * table, int offset);
/* Methods */
void Dofmap_Get_List(dofmap_t * dm, int nvert, hypervertex_t * Vs, int * dofs, int * ndofs);
/* Methods in the vtable */
void Dofmap_Get(dofmap_t * dm, hypervertex_t V, int * dofs, int * ndofs);
int  Dofmap_Max_Len(dofmap_t * dm);
void Dofmap_Destroy(dofmap_t * dm);


#endif
