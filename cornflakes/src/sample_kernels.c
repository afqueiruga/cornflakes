#include "sample_kernels.h"



void peri2d_eval(int ninp, real_t * in, real_t * out) {

}

kernel_t peri2d_kernel = {
  .ninp = 3,
  .inp = {
    { .field_number = 0, .dim = 2, .v_start=0, .v_end=1, .name = "x" },
    { .field_number = 1, .dim = 2, .v_start=0, .v_end=1, .name = "v" },
    { .field_number = 2, .dim = 2, .v_start=0, .v_end=1, .name = "X" },
    { .field_number = 3, .dim = 1, .v_start=2, .v_end=2, .name = "alpha" },
    { .field_number = 4, .dim = 2, .v_start=0, .v_end=1, .name = "param" }
  },
  .noutp = 3,
  .outps = {
    { .rank = 1, .ndof=1, .dofs = {
	{ .field_number = 0, .dim = 2, .v_start=0, .v_end=1, .name = "fx" } }
    },
    { .rank = 2, .ndof=1, .dofs = {
	{ .field_number = 0, .dim = 2, .v_start=0, .v_end=1, .name = "kx" } }
    },
    { .rank = 2, .ndof=1, .dofs = {
	{ .field_number = 0, .dim = 2, .v_start=0, .v_end=1, .name = "kv" } }
    },
    { .rank = -1, .ndof = -1 },
    { .rank = -1, .ndof = -1 },
  },
  .eval = peri2d_eval,
  .name = {"peri2d"}
};



