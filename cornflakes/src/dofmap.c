#include "dofmap.h"

#include <stdlib.h>

struct _DOFMAP_VTABLE_t {
  int (*Max_Len)(dofmap_t * self);
  void (*Get)(dofmap_t * self, hypervertex_t V, int * dofs, int *ndofs);
  void (*Destroy)(dofmap_t * self);
};

/*
 * Class interface
 */
void Dofmap_Get(dofmap_t * dm, hypervertex_t V, int * dofs, int * ndofs) {
  dm->vtable->Get(dm,V,dofs,ndofs);
}
int Dofmap_Max_Len(dofmap_t * dm) {
  return dm->vtable->Max_Len(dm);
}
void Dofmap_Get_List(dofmap_t * dm, int nvert, hypervertex_t * Vs, int * dofs, int * ndofs) {
  int i;
  int * dofiter = dofs;
  *ndofs = 0;
  int nd;
  for(i=0; i<nvert; i++) {
    Dofmap_Get(dm, Vs[i], dofiter, &nd);
    dofiter += nd;
    *ndofs += nd;
  }
}
void Dofmap_Destroy(dofmap_t * dm) {
  dm->vtable->Destroy(dm);
}


/*
 * Strided implementation
 */
int Dofmap_Strided_Max_Len(dofmap_t * dm) {
  return dm->U.strided.stride;
}
void Dofmap_Strided_Get(dofmap_t * dm, hypervertex_t V, int * dofs, int * ndofs) {
  int i;
  *ndofs = dm->U.strided.stride;
  for(i=0;i<dm->U.strided.stride; i++) {
    dofs[i] =  dm->U.strided.offset + dm->U.strided.stride * ( (int)V ) + i;
  }
}
void Dofmap_Strided_Destroy(dofmap_t * dm) {

}
const _DOFMAP_VTABLE_t Dofmap_Strided_vtable = {
  .Max_Len = Dofmap_Strided_Max_Len,
  .Get = Dofmap_Strided_Get,
  .Destroy = Dofmap_Strided_Destroy
};
void Dofmap_Strided(dofmap_t * dm, int stride, int offset) {
  dm->vtable = &Dofmap_Strided_vtable;
  dm->U.strided.stride = stride;
  dm->U.strided.offset = offset;
}


/*
 * Tabled implementation.
 */
int Dofmap_Tabled_Max_Len(dofmap_t * dm) {
  return dm->U.tabled.stride;
}
void Dofmap_Tabled_Get(dofmap_t * dm, hypervertex_t V, int * dofs, int *ndofs) {
  *ndofs = dm->U.tabled.stride;
  int i;
  for(i=0;i<dm->U.tabled.stride; i++) {
    dofs[i] = dm->U.tabled.table[dm->U.tabled.stride*( (int)V ) + i];
  }
}
void Dofmap_Tabled_Destroy(dofmap_t * dm) {
  free(dm->U.tabled.table);
}
const _DOFMAP_VTABLE_t Dofmap_Tabled_vtable = {
  .Max_Len = Dofmap_Tabled_Max_Len,
  .Get = Dofmap_Tabled_Get,
  .Destroy = Dofmap_Tabled_Destroy
};
void Dofmap_Tabled(dofmap_t * dm, int Nentry, int stride, int * table) {
  dm->vtable = &Dofmap_Tabled_vtable;
  dm->U.tabled.stride = stride;
  dm->U.tabled.table = malloc( sizeof(int) * Nentry * stride );
  int i;
  for(i=0;i<Nentry*stride;i++) dm->U.tabled.table[i] = table[i];
}
