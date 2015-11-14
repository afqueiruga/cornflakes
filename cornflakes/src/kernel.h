#ifndef __KERNEL_H
#define __KERNEL_H

typedef double real_t;

#define KERNEL_OUTP_MAX 5
#define KERNEL_IN_DOF_MAX 5
#define KERNEL_NAME_MAX 10

/* DEPRECATED
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
*/


/* 
 * NOMENCLATURE: Prefixes for indexing and dimensions:
 * g_ : Global ordering amongst all processors
 * l_ : Local ordering on this processor
 * k_ : Kernel-level ordering for computation (small numbers)
 * g  : geometry (no _ ) Should be 1,2,3 
 * t  : topology (no _ ) Meh....
 */

/*
 * Data type for DOF naming.
 */
typedef struct dof_t {
  int field_number; //Which location in **data is it in?
  int dim; // How many dofs are associated with a vertex?
  int v_start; // What's the position of the first vertex in the edge?
               // -1 Means it is a global DOF (Not associated with the vertices).
               // hardcoded optimization for the super-common thing of having global
               // parameters passed in, like TIME.
  int v_end; // What's the position of the last vertex in the edge?
             // -1 means the field takes in a variable number of vertices! To the end.
  char name[KERNEL_NAME_MAX]; // The name of the field for info grabbing.
} dof_t;

/* 
 * Data type for output naming.
 * Kernels define multiple outputs to save computation time
 * With tangents, conjugate vectors, etc.
 */
typedef struct outp_t {
  int rank; // Scalar, vector, or matrix
  int k_len; // How long is the kernel calculation // VARIABLE LENGTH?????
  // HOW TO MAP IT OUT? // SHOULD JUST BE IN THE MAPPING FUNCTION
} outp_t;

/*
 * Kernel data type
 */
typedef struct kernel_t {
  int ninp;
  dof_t inp[KERNEL_IN_DOF_MAX];

  int noutp;
  outp_t outp[KERNEL_OUTP_MAX]; // The kernels can have multiple outputs
  
  char name[KERNEL_NAME_MAX]; // Identifier for infoing

  void (*eval)(int * len, real_t * in, real_t * out); // The payload. Eval the kernel.
} kernel_t;


int kernel_inp_len(kernel_t * ke, int l_edge);
int kernel_outp_len(kernel_t * ke, int l_edge);

#endif
