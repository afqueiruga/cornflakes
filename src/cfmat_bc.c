#include "cfmat_bc.h"

#include <stdlib.h>

#define data(x) CFMat_BC_Data(x)

real_t * CFMat_BC_Place(cfmat_t * self, int n, int * dofs, real_t * vals) {
  int i,j;
  // I need to extract the square block to call CFMat_Place on the sub mat
  // and then add the rectangular contribution to R
  int mapped[n];
  int n_sub = 0;
  for(i=0;i<n;i++) {
    mapped[i] = IndexMap_Get(data(self)->map, i);
    if(mapped[i]>=0) n_sub++;
  }
  if(n_sub==n) {
    // No BCs in this block. Nothing happens to R
    data(self)->K->vtable->Place(data(self)->K, n, dofs, vals);    
  } else {
    // Need to trim
    // Ri += Kij uj
    /* int nR = n-n_sub; */
    real_t R_vals[n_sub];
    index_t dofs_sub[n_sub];
    for(i=0;i<n_sub;i++) {
      R_vals[i] = 0;
      /* for(j=0;j< */
    }
    data(self)->R->vtable->Place(data(self)->R, n_sub,dofs_sub,R_vals);

    // K += K_sub
    real_t vals_sub[n_sub*n_sub];
    
    data(self)->K->vtable->Place(data(self)->K, n_sub, dofs_sub, vals_sub);
  }
  return vals + n*n;
}
void CFMat_BC_Wipe(cfmat_t * self) {
  CFMat_Wipe(data(self)->K);
  CFData_Wipe(data(self)->R);
}
void CFMat_BC_Finalize(cfmat_t * self) {
  CFMat_Finalize(data(self)->K);
  CFData_Finalize(data(self)->R);
}
void CFMat_BC_Destroy(cfmat_t * self) {
  // I don't own K and R, but I do own the map
  // Nope, I don't even own the map
  free(self->data);
  //IndexMap_Destroy(&self->map);
}

const _CFMAT_VTABLE_t CFMat_BC_vtable = {
  .Place = &CFMat_BC_Place,
  .Destroy = &CFMat_BC_Destroy,
  .Wipe = &CFMat_BC_Wipe,
  .Finalize = &CFMat_BC_Finalize
};

void CFMat_BC_New(cfmat_t * self,
		  cfmat_t * K, cfdata_t * R, cfdata_t * u,
		  indexmap_t* imap ) {
  self->vtable = &CFMat_BC_vtable;
  self->data = malloc(sizeof(cfmat_bc_data_t));
  data(self)->map = imap;
  data(self)->K = K;
  data(self)->R = R;
  data(self)->u = u;
}
