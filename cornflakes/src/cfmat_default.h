#ifndef __CFMAT_DEFAULT_H
#define __CFMAT_DEFAULT_H

typedef struct CFMat_Default_data_t {
  int nalloc;
  real_t * V;
  int * II;
  int * JJ;

  real_t * Viter;
  int * IIiter;
  int * JJiter;
} CFMat_Default_data_t;

#endif
