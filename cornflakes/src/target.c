#include "target.h"

struct _TARGET_VTABLE_t {
  void (*Place)(target_t * self, int n, int * dofs, real_t * vals);
  void (*Destroy)(target_t * self);
  void (*Wipe)(target_t * self);
};

/*
 * Class interface
 */
void Target_Place(target_t * self, int n, int * dofs, real_t * vals) {
  self->vtable->Place(self,n,dofs,vals);
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

void Target_Default_Place(target_t * self, int n, int * dofs, real_t * vals) {
  
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
  
}

const _TARGET_VTABLE_t Table_Default_vtable = {
  .Place = &Target_Default_Place;
  .Destroy = &Target_Default_Destroy;
  .Wipe = &Target_Default_Wipe;
};
void Target_Default_New(target_t * self, int ioutp,
			kernel_t * ke, hypergraph_t * hg, int ndof) {
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
