#include "target.h"

struct _TARGET_VTABLE_t {
  real_t * (*Place)(target_t * self, int n, int * dofs, real_t * vals);
  void (*Destroy)(target_t * self);
  void (*Wipe)(target_t * self);
};

/*
 * Class interface
 */
real_t * Target_Place(target_t * self, int n, int * dofs, real_t * vals) {
  return self->vtable->Place(self,n,dofs,vals);
}
void Target_Destroy(target_t * self) {
  self->vtable->Destroy(self);
}
void Target_Wipe(target_t * self) {
  self->vtable->Wipe(self);
}


/*
 * The default implementation
 */
real_t * Target_Default_Place(target_t * self,
			      int n, int * dofs, real_t * ker_out) {
  int i,j;
  switch(self->rank) {
  case 2:
    // Fill this block
    for(i=0;i<n;i++) {
      for(j=0;j<n;j++) {
	self->data->IIiter[n*i + j ] = dofs[i];
	self->data->JJiter[n*i + j ] = dofs[j];
	self->data->Viter [n*i + j ] = ker_out[ n*i + j];
      }
    }
    // Advance our iterators.
    self->data->IIiter += n*n;
    self->data->JJiter += n*n;
    self->data->Viter += n*n;
    return ker_out + n*n;
    break;
  case 1:
    for(i=0;i<n;i++) {
      self->data->V[dofs[i]] += ker_out[i];
    }
    return ker_out + n;
    break;
  default:
    for(i=0;i<n;i++) {
      self->data->V[0] += ker_out[i];
    }
    return ker_out + 1;
    break;
  }
}
void Target_Default_Destroy(target_t * self) {
  free(self->data->V);
  if(self->rank >= 2) {
    free(self->data->II);
    free(self->data->JJ);
  }
  free(self->data);
}
void Target_Default_Wipe(target_t * self) {
  int i;
  switch(self->rank) {
  case 0:
    self->data->V[0] = 0.0;
    break;
  case 1:
    for(i=0;i<self->N;i++) {
      self->data->V[i]=0.0;
    }
    break;
  default:
    self->data->Viter = self->data->V;
    self->data->IIiter = self->data->II;
    self->data->JJiter = self->data->JJ;
  }
}

const _TARGET_VTABLE_t Table_Default_vtable = {
  .Place = &Target_Default_Place;
  .Destroy = &Target_Default_Destroy;
  .Wipe = &Target_Default_Wipe;
};
void Target_Default_New(target_t * self, int ioutp,
			kernel_t * ke, hypergraph_t * hg, int ndof) {
  int matsize;
  
  self->vtable = Table_Default_vtable;
  self->data = malloc(sizeof(struct Target_Default_data_t));

  self->rank = ke->outp[ioutp].rank;
  self->N = ndof;
  switch(self->rank) {
  case 0:
    self->data->V = malloc( sizeof(real_t) );
    break;
  case 1:
    self->data->rank = 1;
    self->data->V = malloc( ndof*sizeof(real_t) );
    break;
  case 2:
    matsize = 0;
    for( i=0; i<hg->n_types; i++ ) {
      oplen = kernel_outp_len(ke,ke->outp+j,hg->he[i].l_edge);
      matsize += hg->he[i].n_edge * oplen;
    }
    self->data->V = malloc( matsize * sizeof(real_t) );
    self->data->II = malloc( matsize * sizeof(int) );
    self->data->JJ = malloc( matsize * sizeof(int) );
      
    break;
  default:
    printf("Cannot handle higher rank\n!");
  }
}
