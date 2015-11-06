#ifndef __KERNEL_H
#define __KERNEL_H

typedef double real_t;

#define KERNEL_OUT_DOF_MAX 5
#define KERNEL_IN_DOF_MAX 5

enum dof_loc {
  LOC_GLOBAL,
  LOC_NODE,
  LOC_EDGE
};

typedef struct dof_t {
  enum dof_loc loc;
  int len;
} dof_t;

typedef struct kernel_t {
  int ninp;
  dof_t inp[KERNEL_IN_DOF_MAX];
  int noutp;
  dof_t outp[KERNEL_IN_DOF_MAX];
  void (*assem)(real_t* in, real_t*out);
} kernel_t ;


int kernel_inp_len(kernel_t * ke, int l_edge);
int kernel_outp_len(kernel_t * ke, int l_edge);



extern kernel_t particle_kernel_strct ;

#endif
