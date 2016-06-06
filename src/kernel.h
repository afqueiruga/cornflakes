#ifndef __KERNEL_H
#define __KERNEL_H

typedef double real_t;

#define KERNEL_MAP_MAX 10
#define KERNEL_INP_MAX 10
#define KERNEL_OUTP_MAX 10
#define KERNEL_OUT_MAP_MAX 10
// User only gets 31 characters to name
#define KERNEL_NAME_MAX 32


/* 
 * NOMENCLATURE: Prefixes for indexing and dimensions:
 * g_ : Global ordering amongst all processors
 * l_ : Local ordering on this processor
 * k_ : Kernel-level ordering for computation (small numbers)
 * g  : geometry (no _ ) Should be 1,2, or 3  (or even 4 O_o)
 * t  : topology (no _ ) Meh....
 */

/*
 * The k_map typedef. It's a function that takes in the edge, and tells the assembler
 * Which vertices its field lives on, plus their dimension
 * One of them should look like this:

void kmapVec (..) {
*n = 2;
  verts[0] = edge[0];
  verts[1] = edge[1];
  *dim = 2;
}
-or-
void kmapVec (..) {
  *n = (l_edge-1) / 2;
  for( int i=  0 ; i < *n ; i++ ) {
    verts[i] = edge[ start +  i*stride ];
  }
  *dim = 3
}
 */
typedef void (*k_map_t)(int * edge, int l_edge, int * verts, int * n, int * dim);

typedef struct inp_t {
  int field_number;
  int map_num;
  char name[KERNEL_NAME_MAX];
} inp_t;
/* 
 * Data type for output naming.
 * Kernels define multiple outputs to save computation time
 * With tangents, conjugate vectors, etc.
 */
typedef struct outp_t {
  int rank;
  int nmap;
  int map_nums[KERNEL_OUT_MAP_MAX];
  char name[KERNEL_NAME_MAX];
} outp_t;

typedef struct kernel_t {
  int nmap;
  k_map_t maps[KERNEL_MAP_MAX];
  
  int ninp;
  inp_t inp[KERNEL_INP_MAX];
  int noutp;
  outp_t outp[KERNEL_OUTP_MAX];
  
  char name[KERNEL_NAME_MAX];

  void (*eval)(int len, const real_t * restrict in, real_t * restrict out);
} kernel_t;

int kernel_outps_len(kernel_t * ke, int l_edge);
int kernel_outp_len(kernel_t * ke, outp_t * ou, int l_edge);
int kernel_inps_len(kernel_t * ke, int l_edge);
int kernel_outp_ndof(kernel_t * ke, outp_t * ou, int l_edge);

#endif
