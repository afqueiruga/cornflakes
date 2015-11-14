#include "kernel.h"

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
  .outp = {
    { .rank = 1, .len = 4 },
    { .rank = 2, .len = 4 },
    { .rank = 2, .len = 4 },
    { .rank = -1, .len = -1 },
    { .rank = -1, .len = -1 },
  },
  .name = {"peri2d"}
};
