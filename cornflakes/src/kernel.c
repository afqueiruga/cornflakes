#include "kernel.h"

#include "math.h"


int kernel_inp_len(kernel_t * ke,int l_edge) {
  int i;/*
  int len_loc_in = 0;
  for(i=0; i<ke->ninp; i++) {
    if(ke->inp[i].loc == LOC_NODE) {
      len_loc_in += ke->inp[i].len*l_edge;
    } else {
      len_loc_in += ke->inp[i].len;
    }
  }
  return len_loc_in;*/ return 0;
}
int kernel_outp_len(kernel_t * ke,int l_edge) {
  int i;/*
  int len_loc_out = 0;
  for(i=0; i<ke->noutp; i++) {
    if(ke->outp[i].loc == LOC_NODE) {
      len_loc_out += ke->outp[i].len*l_edge;
    } else {
      len_loc_out += ke->outp[i].len;
    }
  }
  return len_loc_out;*/ return 0;
}
