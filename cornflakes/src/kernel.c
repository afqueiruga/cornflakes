#include "kernel.h"

#include "math.h"

int kernel_dof_len(dof_t * d, int l_edge) {
  if(d->v_start<0) {
    /* Global (no vertices) */
    return d->dim;
  } else {
    if(d->v_end<0) {
      /* Variable-length (vertices from start to the end of the edge )*/
      return d->dim*(l_edge-d->v_start);
    } else {
      /* Fixed-length (the standard case) */
      return d->dim*(d->v_end-d->v_start);
    }
  }   
}


int kernel_inp_len(kernel_t * ke, int l_edge) {
  int i, ndof=0;
  for(i=0;i<ke->ninp;i++) {
    ndof += kernel_dof_len(ke->inp + i, l_edge);
  }
  return ndof;
}
int kernel_outp_len(outp_t * ou, int l_edge) {
  int i, len;
  if(ou->rank==0) return 1;
  for(i=0;i<ou->ndof;i++) {
    len += kernel_dof_len(ou->dofs + i, l_edge);
  }
  if(ou->rank==1)
    return len;
  else
    return len*len;
  }
int kernel_outps_len(kernel_t * ke,int l_edge) {
  int i,j, len=0;
  for(i=0;i<ke->noutp;i++) {
    len += kernel_outp_len(ke->outps+i, l_edge);
  }
  return len;
}

