#include "cfmat.h"

real_t * CFMat_Place(cfmat_t * self, int n, int * dofs, real_t * vals) {
  return self->vtable->Place(self,n,dofs,vals);
}
void CFMat_Set_Value(cfmat_t * self, int i, int j, real_t v) {
  self->vtable->Set_Value(self,i,j,v);
}
void CFMat_Wipe(cfmat_t * self) {
  self->vtable->Wipe(self);
}
void CFMat_Finalize(cfmat_t * self) {
  self->vtable->Finalize(self);
}
void CFMat_Destroy(cfmat_t * self) {
  self->vtable->Destroy(self);
}
