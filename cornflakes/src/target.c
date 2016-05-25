#include "target.h"

#include "stdlib.h"
#include "stdio.h"


/*
 * Class interface
 */
real_t * Target_Place(target_t * self, int n, int * dofs, real_t * vals) {
  switch(self->rank) {
  case 2:
    return CFMat_Place(self->K,n,dofs,vals);
    break;
  default: //0,1
    return CFData_Place(self->R,n,dofs,vals);
  }
}
void Target_Destroy(target_t * self) {
    switch(self->rank) {
  case 2:
    return CFMat_Destroy(self->K);
    break;
  default: //0,1
    return CFData_Destroy(self->R);
  }
}
void Target_Wipe(target_t * self) {
  switch(self->rank) {
  case 2:
    CFMat_Wipe(self->K);
    break;
  default: //0,1
    CFData_Wipe(self->R);
  }
}
void Target_Finalize(target_t * self) {
  switch(self->rank) {
  case 2:
    CFMat_Finalize(self->K);
    break;
  default: //0,1
    CFData_Finalize(self->R);
  }
}
#if 0

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

const _TARGET_VTABLE_t Table_Default_vtable = {
  .Place = &Target_Default_Place,
  .Destroy = &Target_Default_Destroy,
  .Wipe = &Target_Default_Wipe,
  .Finalize = &Target_Default_Finalize
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
#endif
