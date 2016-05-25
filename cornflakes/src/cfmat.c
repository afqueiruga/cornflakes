#include "cfmat.h"

void CFMat_Place(cfmat_t * self, int n, int * dofs, real_t * vals) {
  self->vtable->Place(self,n,dofs,vals);
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
