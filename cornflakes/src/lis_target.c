#ifdef USE_LIS
#include "lis_target.h"

void Target_LIS_From_Obj(target_t * self, int rank, void * obj) {
  
}
#define data(x) (Target_LIS_Data(x))
void Target_LIS_Destroy(target_t * self) {
  if(data(self)->own) {
    switch(self->rank) {
    case 0:
    case 1:
      lis_vector_destroy((data(self)->R));
      break;
    case 2:
      lis_matrix_unset((data(self)->K));
      lis_matrix_destroy((data(self)->K));
    break;
    default:
      //
      break;
    }
  }
}

real_t * Target_LIS_Place(target_t * self,
			  int n, int *dofs, real_t * ker_out)
{
  switch(self->rank) {
  case 2:
    for(int i=0;i<n;i++) {
      for(int j=0;j<n;j++) {
	lis_matrix_set_value(LIS_ADD_VALUE, dofs[i],dofs[j], ker_out[n*i+j], data(self)->K);
      }
    }
    return ker_out + n*n;
  case 1:
    lis_vector_set_values(LIS_ADD_VALUE, n,dofs, ker_out, data(self)->R);
    return ker_out + n;
  default:
    lis_vector_set_value(LIS_ADD_VALUE, 0, ker_out[0], data(self)->R);

    return ker_out + 1;
  }
}
#endif
