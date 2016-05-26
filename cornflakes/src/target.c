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


void Target_New(target_t * self, int onum,
		kernel_t * ke, hypergraph_t * hg, int ndof,
		char * backend) {
  /* unimplemented */
}
void Target_New_From_Ptr(target_t * self, int rank,  void * payload) {
  self->rank = rank;
  switch(self->rank) {
  case 2:
    self->K = payload;
    break;
  case 1:
    self->R = payload;
  }
}
