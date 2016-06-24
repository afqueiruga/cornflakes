#include "cfmat_bc.h"

#include <stdlib.h>

#define data(x) CFMat_BC_Data(x)

real_t * CFMat_BC_Place(cfmat_t * self, int n, int * dofs, real_t * vals) {
  int i,j;
  // I need to extract the square block to call CFMat_Place on the sub mat
  // and then add the rectangular contribution to R
  index_t mapped[n];
  index_t dofs_sub[n];
  index_t i_sub[n], i_bc[n];
  int n_sub = 0, n_bc=0;
  for(i=0;i<n;i++) {
    mapped[i] = IndexMap_Get(data(self)->map, dofs[i]);
    /* printf("%d ",mapped[i]); */
    if(mapped[i]>=0) {
      dofs_sub[n_sub] = mapped[i];
      i_sub[n_sub] = i;
      n_sub++;
    } else {
      i_bc[n_bc] = i;
      n_bc++;
    }
  }
  /* printf("\n"); */
  /* for(i=0;i<n_sub;i++) printf("%d ",dofs_sub[i]); printf("\n"); */
  /* for(i=0;i<n_bc;i++) printf("%d ",i_bc[i]); printf("\n"); */
  if(n_sub==n) {
    // No BCs in this block. Nothing happens to R
    // But, they still map to new indices
    data(self)->K->vtable->Place(data(self)->K, n, dofs_sub, vals);    
  } else {
    // Need to trim
    // Fill up trimmed dofs
    /* index_t dofs_sub[n_sub]; */
    
    // Ri += Kij uj
    /* Calculate Ri += Kij uj */
    real_t R_vals[n_sub];
    real_t ubar[n_bc];
    int bcdof[n_bc];
    for(i=0;i<n_bc;i++) bcdof[i] = dofs[i_bc[i]];
    CFData_Get_Values(data(self)->u, n_bc, bcdof ,ubar);
    for(i=0;i<n_sub;i++) {
      R_vals[i] = 0;
      for(j=0;j<n_bc;j++) {
	R_vals[i] += vals[i*n + i_bc[j]] * ubar[j] ;
      }
    }
    /* for(i=0;i<n_sub;i++) printf("%lf ",R_vals[i]); printf("\n"); */
    data(self)->R->vtable->Place(data(self)->R, n_sub,dofs_sub,R_vals);


    /* Copy up K to shrink it */
    // K += K_sub
    real_t vals_sub[n_sub*n_sub];
    for(i=0;i<n_sub;i++) {
      for(j=0;j<n_sub;j++) {
	vals_sub[i*n_sub+j] = vals[i_sub[i]*n+i_sub[j]];
      }
    }
    data(self)->K->vtable->Place(data(self)->K, n_sub, dofs_sub, vals_sub);
  }
  return vals + n*n;
}
void CFMat_BC_Set_Value(cfmat_t * self, int i, int j, real_t v) {
  // TODO THIS IS WHERE YOU ARE WORKING
  int mi,mj;
  mi = IndexMap_Get(data(self)->map, i);
  if(mi>=0) {
    mj = IndexMap_Get(data(self)->map, j);
    if(mj>=0) {
      // Add to Kij
      CFMat_Set_Value(data(self)->K, mi,mj,v);
    } else {

    }
  } else {
    // Add to Ri += Kij uj
    mj = IndexMap_Get(data(self)->map, j);
    if(mj>=0) {
      real_t ubar;
      CFData_Get_Values(data(self)->u,1,&i,&ubar);
      ubar *= v;
      CFData_Place(data(self)->R, 1,&mj,&ubar);
    }
  }
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
  .Set_Value = CFMat_BC_Set_Value,
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
  self->N = imap->end - imap->start;
  self->own = 0;
}
