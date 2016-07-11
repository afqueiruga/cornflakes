#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

#include "cornflakes.h"


/* Just copy-and-pasted in a simple kernel */
void linear_spring_eval_wr(int l_edge, const real_t * restrict in, real_t * restrict out) {
  out[0] = in[0];
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
void kmap_linear_spring2(int * edge, int l_edge, int * verts, int * n, int * dim) {
  int i;
  *dim = 1;
  *n = (int)(1);
  if(edge) verts[0]=0;
}
kernel_t kernel_linear_spring = {
  .nmap = 3,
  .maps = {kmap_linear_spring0,kmap_linear_spring1,kmap_linear_spring2},
  .ninp = 3,
  .inp = {{ .field_number=0, .map_num=0, .name="x" },{ .field_number=1, .map_num=0, .name="v" },{ .field_number=2, .map_num=1, .name="params" }},
  
  .noutp = 1,
  .outp = {{ .rank = 1, .nmap = 1, .map_nums = { 2 } }},
  .eval=linear_spring_eval_wr,
  .name="linear_spring"
};



int main(int argc, char **argv) {
  /* The problem set up. I need data objects for this test, but they don't
     need to be filled with anything meaningful. */
  hypergraph_t Bonds, htrue;
  dofmap_t NodeVec, Global;
  cfdata_t x,v,params;
  int ndof = 2*5;
  
  /* Make a graph */
  Hypergraph_Alloc(&Bonds, 1);
#define ADDPAIR(i,j) {int P[2] = {i,j}; Hypergraph_Push_Edge(&Bonds, 2, P);}
  ADDPAIR(0,1);
  ADDPAIR(1,2);
  ADDPAIR(0,3);
  ADDPAIR(2,3);
  ADDPAIR(0,4);
#undef ADDPAIR

  /* Setup dmaps */
  Dofmap_Strided(&NodeVec, 2,0);
  Dofmap_Strided(&Global, 2,0);
  
  /* Initialize CFData */
  CFData_Default_New(&x,ndof);
  CFData_Default_New(&v,ndof);
  CFData_Default_New(&params,2);
  CFData_Wipe(&x);
  {int i=0;real_t v=1.0;CFData_Place(&x, 1,&i, &v);}
  
  /* Call filter */
  dofmap_t *dms[2] = {&NodeVec,&Global};
  cfdata_t *data[3] = {&x,&v,&params};
  filter(&kernel_linear_spring,&Bonds,dms,data, &htrue);

  /* Make the checks */
  assert(htrue.n_types == 1);
  assert(htrue.he[0].n_edge == 3);

  Dofmap_Destroy(&NodeVec);
  Dofmap_Destroy(&Global);
  CFData_Destroy(&x);
  CFData_Destroy(&v);
  CFData_Destroy(&params);

  Hypergraph_Destroy(&Bonds);
  Hypergraph_Destroy(&htrue);
  return 0;
}
