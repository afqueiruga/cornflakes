#include "cfdata.h"

void CFData_Get_Values(cfdata_t * self, int ndof,int * dofs, real_t * vals) {
  self->vtable->Get_Values(self,ndof,dofs,vals);
}
real_t * CFData_Place(cfdata_t * self, int n, int * dofs, real_t * vals) {
  return self->vtable->Place(self,n,dofs,vals);
}
void CFData_Scatter(cfdata_t * self, real_t * src) {
  self->vtable->Scatter(self,src);
}
void CFData_Wipe(cfdata_t * self) {
  self->vtable->Wipe(self);
}
void CFData_Finalize(cfdata_t * self) {
  self->vtable->Finalize(self);
}
void CFData_Destroy(cfdata_t * self) {
  self->vtable->Destroy(self);
}
void CFData_Get_Ptr(cfdata_t * self, real_t **dat) {
  self->vtable->Get_Ptr(self,dat);
}
void CFData_Release_Ptr(cfdata_t * self, real_t **dat) {
  self->vtable->Get_Ptr(self,dat);
}
