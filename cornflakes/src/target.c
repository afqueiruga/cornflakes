#include "target.h"

#include "stdlib.h"
#include "stdio.h"


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
 * Helper method for making a target with a specific backend
 */
void Target_New(target_t * self, int onum,
		kernel_t * ke, hypergraph_t * hg, int ndof,
		char * backend) {
  switch(backend[0]) {
  case 'd':
  default:
    Target_Default_New(self,onum,ke,hg,ndof);
  }
}
	

/*
 * The default implementation
 */
#define data(x) Target_Default_Data(x)
real_t * Target_Default_Place(target_t * self,
			      int n, int * dofs, real_t * ker_out) {
  int i,j;
  switch(self->rank) {
  case 2:
    // Fill this block
    for(i=0;i<n;i++) {
      for(j=0;j<n;j++) {
	data(self)->IIiter[n*i + j ] = dofs[i];
	data(self)->JJiter[n*i + j ] = dofs[j];
	data(self)->Viter [n*i + j ] = ker_out[ n*i + j];
      }
    }
    // Advance our iterators.
    data(self)->IIiter += n*n;
    data(self)->JJiter += n*n;
    data(self)->Viter += n*n;
    return ker_out + n*n;
    break;
  case 1:
    for(i=0;i<n;i++) {
      data(self)->V[dofs[i]] += ker_out[i];
    }
    return ker_out + n;
    break;
  default:
    for(i=0;i<n;i++) {
      data(self)->V[0] += ker_out[i];
    }
    return ker_out + 1;
    break;
  }
}
void Target_Default_Destroy(target_t * self) {
  if(data(self)->own) {
    free(data(self)->V);
    if(self->rank >= 2) {
      free(data(self)->II);
      free(data(self)->JJ);
    }
  }
  free(data(self));
}
void Target_Default_Wipe(target_t * self) {
  int i;
  switch(self->rank) {
  case 0:
    data(self)->V[0] = 0.0;
    break;
  case 1:
    for(i=0;i<self->N;i++) {
      data(self)->V[i]=0.0;
    }
    break;
  default:
    data(self)->Viter = data(self)->V;
    data(self)->IIiter = data(self)->II;
    data(self)->JJiter = data(self)->JJ;
  }
}

const _TARGET_VTABLE_t Table_Default_vtable = {
  .Place = &Target_Default_Place,
  .Destroy = &Target_Default_Destroy,
  .Wipe = &Target_Default_Wipe
};
void Target_Default_New(target_t * self, int onum,
			kernel_t * ke, hypergraph_t * hg, int ndof) {
  int matsize, oplen;
  int i;
  self->vtable = &Table_Default_vtable;
  self->data = malloc(sizeof(struct Target_Default_data_t));

  self->rank = ke->outp[onum].rank;
  self->N = ndof;
  data(self)->own = 1;
  switch(self->rank) {
  case 0:
    data(self)->V = malloc( sizeof(real_t) );
    break;
  case 1:
    data(self)->V = malloc( ndof*sizeof(real_t) );
    break;
  case 2:
    matsize = 0;
    for( i=0; i<hg->n_types; i++ ) {
      oplen = kernel_outp_len(ke,ke->outp+onum,hg->he[i].l_edge);
      matsize += hg->he[i].n_edge * oplen;
    }
    data(self)->V = malloc( matsize * sizeof(real_t) );
    data(self)->II = malloc( matsize * sizeof(int) );
    data(self)->JJ = malloc( matsize * sizeof(int) );
      
    break;
  default:
    printf("Cannot handle higher rank\n!");
  }
  // For good measure:
  Target_Default_Wipe(self);
}
void Target_Default_From_Array(target_t * self, int rank, int ndof,
			       real_t * V, int * II, int * JJ) {
  self->vtable = &Table_Default_vtable;
  self->data = malloc(sizeof(struct Target_Default_data_t));

  self->rank = rank;
  self->N = ndof;

  data(self)->own = 0;
  
  data(self)->V = V;
  if(self->rank > 1) {
    data(self)->II = II;
    data(self)->JJ = JJ;
  }
  Target_Default_Wipe(self);
}
