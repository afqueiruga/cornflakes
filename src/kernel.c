#include "kernel.h"

#include <math.h>
#include <stdlib.h>



int kernel_map_len(k_map_t  km, int l_edge) {
  int dim;
  int nsel;
  km(NULL,l_edge,NULL, &nsel, &dim);
  return nsel*dim;
}

int kernel_inps_len(kernel_t * ke, int l_edge) {
  int i, ndof=0;
  for(i=0;i<ke->ninp;i++) {
    ndof += kernel_map_len(ke->maps[ ke->inp[i].map_num ], l_edge);
  }
  return ndof;
}

int kernel_outp_ndof(kernel_t * ke, outp_t * ou, int l_edge) {
  int i, len=0;
  if(ou->rank==0) return 1;
  for(i=0;i<ou->nmap;i++) {
    len += kernel_map_len(ke->maps[ ou->map_nums[i] ], l_edge);
  }
  return len;
}
int kernel_outp_len(kernel_t * ke, outp_t * ou, int l_edge) {
  int i, len=0;
  if(ou->rank==0) return 1;
  for(i=0;i<ou->nmap;i++) {
    len += kernel_map_len(ke->maps[ ou->map_nums[i] ], l_edge);
  }
  if(ou->rank==1)
    return len;
  else
    return len*len;
}

int kernel_outps_len(kernel_t * ke, int l_edge) {
  int i,j, len=0;
  for(i=0;i<ke->noutp;i++) {
    len += kernel_outp_len(ke, ke->outp+i, l_edge);
  }
  return len;
}
