#include "cfmat.h"

real_t * CFMat_Place(cfmat_t * self, int ln, int * ldofs,int rn, int * rdofs,  real_t * vals) {
  return self->vtable->Place(self,ln,ldofs,rn,rdofs,vals);
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
