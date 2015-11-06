#include "kernel.h"

#include "math.h"


int kernel_inp_len(kernel_t * ke,int l_edge) {
  int i;
  int len_loc_in = 0;
  for(i=0; i<ke->ninp; i++) {
    if(ke->inp[i].loc == LOC_NODE) {
      len_loc_in += ke->inp[i].len*l_edge;
    } else {
      len_loc_in += ke->inp[i].len;
    }
  }
  return len_loc_in;
}
int kernel_outp_len(kernel_t * ke,int l_edge) {
  int i;
  int len_loc_out = 0;
  for(i=0; i<ke->noutp; i++) {
    if(ke->outp[i].loc == LOC_NODE) {
      len_loc_out += ke->outp[i].len*l_edge;
    } else {
      len_loc_out += ke->outp[i].len;
    }
  }
  return len_loc_out;
}





void eval_peri(/*Inputs:*/
real_t * x,
real_t * v,
real_t * X,
real_t * params
,
/*Outputs:*/
double * F,
double * K)
{
/* Evalaute F:*/
F[0] = params[0]*(x[0] - x[2])*(-sqrt(pow(-X[0] + X[2], 2) + pow(-X[1] + X[3], 2)) + sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)))/sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2));
F[1] = params[0]*(x[1] - x[3])*(-sqrt(pow(-X[0] + X[2], 2) + pow(-X[1] + X[3], 2)) + sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)))/sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2));
F[2] = params[0]*(-x[0] + x[2])*(-sqrt(pow(-X[0] + X[2], 2) + pow(-X[1] + X[3], 2)) + sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)))/sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2));
F[3] = params[0]*(-x[1] + x[3])*(-sqrt(pow(-X[0] + X[2], 2) + pow(-X[1] + X[3], 2)) + sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)))/sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2));

/* Evalaute K:*/
K[0] = params[0]*pow(x[0] - x[2], 2)*sqrt(pow(X[0] - X[2], 2) + pow(X[1] - X[3], 2))/pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 3.0L/2.0L) - params[0]*sqrt(pow(X[0] - X[2], 2) + pow(X[1] - X[3], 2))/sqrt(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2)) + params[0];
K[1] = params[0]*(x[0]*x[1] - x[0]*x[3] - x[1]*x[2] + x[2]*x[3])*sqrt(pow(X[0], 2) - 2*X[0]*X[2] + pow(X[1], 2) - 2*X[1]*X[3] + pow(X[2], 2) + pow(X[3], 2))/pow(pow(x[0], 2) - 2*x[0]*x[2] + pow(x[1], 2) - 2*x[1]*x[3] + pow(x[2], 2) + pow(x[3], 2), 3.0L/2.0L);
K[2] = -params[0]*pow(x[0] - x[2], 2)*sqrt(pow(X[0] - X[2], 2) + pow(X[1] - X[3], 2))/pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 3.0L/2.0L) + params[0]*sqrt(pow(X[0] - X[2], 2) + pow(X[1] - X[3], 2))/sqrt(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2)) - params[0];
K[3] = -params[0]*(x[0]*x[1] - x[0]*x[3] - x[1]*x[2] + x[2]*x[3])*sqrt(pow(X[0], 2) - 2*X[0]*X[2] + pow(X[1], 2) - 2*X[1]*X[3] + pow(X[2], 2) + pow(X[3], 2))/pow(pow(x[0], 2) - 2*x[0]*x[2] + pow(x[1], 2) - 2*x[1]*x[3] + pow(x[2], 2) + pow(x[3], 2), 3.0L/2.0L);
K[4] = params[0]*(x[0]*x[1] - x[0]*x[3] - x[1]*x[2] + x[2]*x[3])*sqrt(pow(X[0], 2) - 2*X[0]*X[2] + pow(X[1], 2) - 2*X[1]*X[3] + pow(X[2], 2) + pow(X[3], 2))/pow(pow(x[0], 2) - 2*x[0]*x[2] + pow(x[1], 2) - 2*x[1]*x[3] + pow(x[2], 2) + pow(x[3], 2), 3.0L/2.0L);
K[5] = params[0]*pow(x[1] - x[3], 2)*sqrt(pow(X[0] - X[2], 2) + pow(X[1] - X[3], 2))/pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 3.0L/2.0L) - params[0]*sqrt(pow(X[0] - X[2], 2) + pow(X[1] - X[3], 2))/sqrt(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2)) + params[0];
K[6] = -params[0]*(x[0]*x[1] - x[0]*x[3] - x[1]*x[2] + x[2]*x[3])*sqrt(pow(X[0], 2) - 2*X[0]*X[2] + pow(X[1], 2) - 2*X[1]*X[3] + pow(X[2], 2) + pow(X[3], 2))/pow(pow(x[0], 2) - 2*x[0]*x[2] + pow(x[1], 2) - 2*x[1]*x[3] + pow(x[2], 2) + pow(x[3], 2), 3.0L/2.0L);
K[7] = -params[0]*pow(x[1] - x[3], 2)*sqrt(pow(X[0] - X[2], 2) + pow(X[1] - X[3], 2))/pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 3.0L/2.0L) + params[0]*sqrt(pow(X[0] - X[2], 2) + pow(X[1] - X[3], 2))/sqrt(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2)) - params[0];
K[8] = -params[0]*pow(x[0] - x[2], 2)*sqrt(pow(X[0] - X[2], 2) + pow(X[1] - X[3], 2))/pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 3.0L/2.0L) + params[0]*sqrt(pow(X[0] - X[2], 2) + pow(X[1] - X[3], 2))/sqrt(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2)) - params[0];
K[9] = -params[0]*(x[0]*x[1] - x[0]*x[3] - x[1]*x[2] + x[2]*x[3])*sqrt(pow(X[0], 2) - 2*X[0]*X[2] + pow(X[1], 2) - 2*X[1]*X[3] + pow(X[2], 2) + pow(X[3], 2))/pow(pow(x[0], 2) - 2*x[0]*x[2] + pow(x[1], 2) - 2*x[1]*x[3] + pow(x[2], 2) + pow(x[3], 2), 3.0L/2.0L);
K[10] = params[0]*pow(x[0] - x[2], 2)*sqrt(pow(X[0] - X[2], 2) + pow(X[1] - X[3], 2))/pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 3.0L/2.0L) - params[0]*sqrt(pow(X[0] - X[2], 2) + pow(X[1] - X[3], 2))/sqrt(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2)) + params[0];
K[11] = params[0]*(x[0]*x[1] - x[0]*x[3] - x[1]*x[2] + x[2]*x[3])*sqrt(pow(X[0], 2) - 2*X[0]*X[2] + pow(X[1], 2) - 2*X[1]*X[3] + pow(X[2], 2) + pow(X[3], 2))/pow(pow(x[0], 2) - 2*x[0]*x[2] + pow(x[1], 2) - 2*x[1]*x[3] + pow(x[2], 2) + pow(x[3], 2), 3.0L/2.0L);
K[12] = -params[0]*(x[0]*x[1] - x[0]*x[3] - x[1]*x[2] + x[2]*x[3])*sqrt(pow(X[0], 2) - 2*X[0]*X[2] + pow(X[1], 2) - 2*X[1]*X[3] + pow(X[2], 2) + pow(X[3], 2))/pow(pow(x[0], 2) - 2*x[0]*x[2] + pow(x[1], 2) - 2*x[1]*x[3] + pow(x[2], 2) + pow(x[3], 2), 3.0L/2.0L);
K[13] = -params[0]*pow(x[1] - x[3], 2)*sqrt(pow(X[0] - X[2], 2) + pow(X[1] - X[3], 2))/pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 3.0L/2.0L) + params[0]*sqrt(pow(X[0] - X[2], 2) + pow(X[1] - X[3], 2))/sqrt(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2)) - params[0];
K[14] = params[0]*(x[0]*x[1] - x[0]*x[3] - x[1]*x[2] + x[2]*x[3])*sqrt(pow(X[0], 2) - 2*X[0]*X[2] + pow(X[1], 2) - 2*X[1]*X[3] + pow(X[2], 2) + pow(X[3], 2))/pow(pow(x[0], 2) - 2*x[0]*x[2] + pow(x[1], 2) - 2*x[1]*x[3] + pow(x[2], 2) + pow(x[3], 2), 3.0L/2.0L);
K[15] = params[0]*pow(x[1] - x[3], 2)*sqrt(pow(X[0] - X[2], 2) + pow(X[1] - X[3], 2))/pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 3.0L/2.0L) - params[0]*sqrt(pow(X[0] - X[2], 2) + pow(X[1] - X[3], 2))/sqrt(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2)) + params[0];
}





void particle_kernel_calc(real_t* in,real_t* out) {
  eval_peri(
	    in,
	    in+4,
	    in+8,
	    in+12,

	    out,
	    out+4
	    );
	    
}

kernel_t particle_kernel_strct = {
  .ninp = 4,
  .inp = { {LOC_NODE,2},
	   {LOC_NODE,2},
	   {LOC_NODE,2},
	   {LOC_GLOBAL,1},
	   {LOC_GLOBAL,0} },
  .noutp = 1,
  .outp = { {LOC_NODE,2},
	    {LOC_GLOBAL,0},
	    {LOC_GLOBAL,0},
	    {LOC_GLOBAL,0},
	    {LOC_GLOBAL,0} },
  .assem = particle_kernel_calc
};

