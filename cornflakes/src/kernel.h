#ifndef __KERNEL_H
#define __KERNEL_H

typedef double real_t;

#define KERNEL_MAP_MAX 5
#define KERNEL_INP_MAX 5
#define KERNEL_OUTP_MAX 5
#define KERNEL_OUT_MAP_MAX 5
// User only gets 31 characters to name
#define KERNEL_NAME_MAX 32

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
 * g  : geometry (no _ ) Should be 1,2, or 3  (or even 4 O_o)
 * t  : topology (no _ ) Meh....
 */

/*
 * Data type for DOF naming.
 */
/*
typedef struct dof_t {
  int field_number; //Which location in **data and **dmap is it in?
  int dim; // How many dofs are associated with a vertex?
  int v_start; // What's the position of the first vertex in the edge?
               // -1 Means it is a global DOF (Not associated with the vertices).
               // hardcoded optimization for the super-common thing of having global
               // parameters passed in, like TIME.
  int v_end; // What's the position of the last vertex in the edge?
             // -1 means the field takes in a variable number of vertices! To the end.
  char name[KERNEL_NAME_MAX]; // The name of the field for info grabbing.
} dof_t;
*/
/* 
 * Data type for output naming.
 * Kernels define multiple outputs to save computation time
 * With tangents, conjugate vectors, etc.
 */
/*
typedef struct outp_t {
  int rank; // Scalar, vector, or matrix
  //int k_len; // How long is the kernel calculation // VARIABLE LENGTH?????
  // Same dof_t struct is used, so the hyperedge vertices have to store sufficient info
  int ndof;
  // for the output mapping as well.
  dof_t dofs[KERNEL_OUT_DOF_MAX];
} outp_t;
*/
/*
 * Kernel data type
 */
/*
typedef struct kernel_t {
  // Each input is directly a dof_t
  int ninp;
  dof_t inp[KERNEL_IN_DOF_MAX];

  // The kernels can have multiple outputs
  // The outputs need their own types
  int noutp;
  outp_t outps[KERNEL_OUTP_MAX]; 
  
  char name[KERNEL_NAME_MAX]; // Identifier for infoing

  void (*eval)(int len, real_t * in, real_t * out); // The payload. Eval the kernel.
} kernel_t;
*/


/*
 * Redoing it again...
 */

typedef struct k_map_t {
  int dim;
  int v_start;
  int v_end;
} k_map_t;

typedef struct inp_t {
  int field_number;
  int map_num;
  char name[KERNEL_NAME_MAX];
} inp_t;

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

int kernel_map_len(k_map_t * km, int l_edge);
int kernel_inps_len(kernel_t * ke, int l_edge);
int kernel_outp_ndof(kernel_t * ke, outp_t * ou, int l_edge);
int kernel_outp_len(kernel_t * ke, outp_t * ou, int l_edge);
int kernel_outps_len(kernel_t * ke, int l_edge);

#endif
