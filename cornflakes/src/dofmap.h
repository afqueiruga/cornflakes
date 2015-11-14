#ifndef __H_DOFMAP
#define __H_DOFMAP

#include "hypergraph.h"

enum dofmap_type  {
  DOFMAP_STRIDED,
  DOFMAP_TABLE
} ;


/* Vtable forward declaration*/
typedef struct _DOFMAP_VTABLE_t _DOFMAP_VTABLE_t;

/* Structure. Using a union for now... Not extensible*/
typedef struct dofmap_t {
  //enum dofmap_type type;
  union {
    struct strided_t{
      int stride;
      int offset;
    }  strided ;
    struct tabled_t{
      int stride;
      int * table;
    } tabled;
  } U;
  const _DOFMAP_VTABLE_t * vtable;
} dofmap_t;

/* Constructors */
void Dofmap_Strided(dofmap_t * dm, int stride, int offset);
void Dofmap_Tabled(dofmap_t * dm, int stride, int * table);
/* Methods */
void Dofmap_Get(dofmap_t * dm, hypervertex_t V, int * dofs, int * ndofs);
int  Dofmap_Max_Len(dofmap_t * dm);


#endif
