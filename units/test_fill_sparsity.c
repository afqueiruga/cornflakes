#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "cornflakes.h"

/* Just copy-and-pasted in a simple kernel */
void linear_spring_eval_wr(int l_edge, const real_t * restrict in, real_t * restrict out) {
  // Not needed for this test.
}
void kmap_linear_spring0(int * edge, int l_edge, int * verts, int * n, int * dim) {
  int i;
  *dim = 2;
  *n = (int)((int)(2));
  if(edge) {
    for( i=0 ; i<*n ; i++ ) {
      verts[i] = edge[(int)(0) + 1*i];
    }
  }
}
void kmap_linear_spring1(int * edge, int l_edge, int * verts, int * n, int * dim) {
  int i;
  *dim = 2;
  *n = (int)(1);
  if(edge) verts[0]=0;
}
kernel_t kernel_linear_spring = {
  .nmap = 2,
  .maps = {kmap_linear_spring0,kmap_linear_spring1},
  .ninp = 3,
  .inp = {{ .field_number=0, .map_num=0, .name="x" },{ .field_number=1, .map_num=0, .name="v" },{ .field_number=2, .map_num=1, .name="params" }},
  
  .noutp = 2,
  .outp = {{ .rank = 1, .nmap = 1, .map_nums = { 0 } },{ .rank = 2, .nmap = 1, .map_nums = { 0 } }},
  .eval=linear_spring_eval_wr,
  .name="linear_spring"
};





int main(int argc, char **argv) {
  /* The problem set up. I don't need any actually data for this test. */
  hypergraph_t Bonds;
  dofmap_t NodeVec;
  dofmap_t Global;
  int ndof = 2*5;
  
  /* Make a graph */
  Hypergraph_Alloc(&Bonds, 1);
#define ADDPAIR(i,j) {int P[2] = {i,j}; Hypergraph_Push_Edge(&Bonds, 2, P);}
  ADDPAIR(0,1);
  ADDPAIR(1,2);
  ADDPAIR(0,3);
  ADDPAIR(2,3);
#undef ADDPAIR

  /* Setup dmaps */
  Dofmap_Strided(&NodeVec, 2,0);
  Dofmap_Strided(&Global, 2,0);

  /* Set up some matrices and targets */
  cfmat_t K;
  cfdata_t R;
  target_t att[2];
  CFData_Default_New(&R,ndof);
  CFMat_Default_New(&K,1,&kernel_linear_spring,&Bonds,ndof);
  Target_New_From_Ptr(att+0, 1,&R);
  Target_New_From_Ptr(att+1, 2,&K);
  
  /* Call fill_sparsity finally */
  dofmap_t *dms[2] = {&NodeVec,&Global};
  fill_sparsity(&kernel_linear_spring,&Bonds,dms, att);

  
  Dofmap_Destroy(&NodeVec);
  Dofmap_Destroy(&Global);
  CFMat_Destroy(&K);
  CFData_Destroy(&R);
  Hypergraph_Destroy(&Bonds);
  return 0;
}
