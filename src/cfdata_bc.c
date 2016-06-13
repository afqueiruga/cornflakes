#include "cfdata_bc.h"

#include <stdlib.h>
#include <stdio.h>

#define data(self) (CFData_BC_Data(self))
real_t * CFData_BC_Place(cfdata_t * self,
			    int n, int * dofs, real_t * ker_out) {
  int i;
  index_t dofs_sub[n];
  real_t vals_sub[n];
  int n_sub = 0;
  for(i=0;i<n;i++) {
    index_t mapped = IndexMap_Get(data(self)->map, dofs[i]);
    if(mapped >= 0) {
      dofs_sub[n_sub] = mapped;
      vals_sub[n_sub] = ker_out[i];
      n_sub++;
    }
  }
  data(self)->R->vtable->Place(data(self)->R, n_sub,dofs_sub,vals_sub);
  return ker_out + n;
}
void CFData_BC_Scatter(cfdata_t * self, real_t * src) {
  /* Also doesn't make sense */
  printf("CFData_BC_Scatter is unimplemented : non-concrete data\n");
}
void CFData_BC_Copy(cfdata_t *self, cfdata_t * src) {
  /* Nope, don't do this one either! */
  printf("CFData_BC_Copy is unimplemented : non-concrete data\n");
}
void CFData_BC_Destroy(cfdata_t * self) {
  /* I don't own anything other than my internet struct of ptrs */
  free(data(self));
}
void CFData_BC_Wipe(cfdata_t * self) {
  CFData_Wipe(data(self)->R);
}
void CFData_BC_Finalize(cfdata_t * self) {
  CFData_Finalize(data(self)->R);
}
void CFData_BC_Get_Values(cfdata_t * self, int ndof,int *dofs, real_t * vals) {
  /* Actually, This routine doesn't make sense! */
  printf("CFData_BC_Get_Values is unimplemented : non-concrete data\n");
  /* CFData_Get_Values(data(self)->R,ndof_sub,dofs_sub,vals); */
}
void CFData_BC_Get_Ptr(cfdata_t * self, real_t **ptr) {
  /* This one also doesn't make sense */
  printf("CFData_BC_Get_Ptr is unimplemented : non-concrete data\n");
}
void CFData_BC_Release_Ptr(cfdata_t * self, real_t **ptr) {
  /* This one also doesn't make sense */
  printf("CFData_BC_Release_Ptr is unimplemented : non-concrete data\n");
}
void CFData_BC_Print(cfdata_t * self) {
  CFData_Print(data(self)->R);
}

/* vtable */
const _CFDATA_VTABLE_t cfdata_bc_vtable = {
  .Get_Values = CFData_BC_Get_Values,
  .Scatter = CFData_BC_Scatter,
  .Copy = CFData_BC_Copy,
  .Place = &CFData_BC_Place,
  .Destroy = &CFData_BC_Destroy,
  .Wipe = &CFData_BC_Wipe,
  .Finalize = &CFData_BC_Finalize,
  .Get_Ptr = &CFData_BC_Get_Ptr,
  .Release_Ptr = &CFData_BC_Release_Ptr,
  .Print = CFData_BC_Print
};

void CFData_BC_New(cfdata_t * self,
		   cfdata_t * R, indexmap_t * map) {
  self->vtable = &cfdata_bc_vtable;
  self->own = 0;
  self->N = map->end - map->start;

  self->data = malloc(sizeof(cfdata_bc_data_t));
  data(self)->map = map;
  data(self)->R = R;
}
