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
double * KX,
double * KV)
{
/* Evalaute F:*/
F[0] = -params[0]*(x[0] - x[2])*(-sqrt(pow(-X[0] + X[2], 2) + pow(-X[1] + X[3], 2)) + sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)))/sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) - params[1]*(-x[0] + x[2])*((-v[0] + v[2])*(-x[0] + x[2]) + (-v[1] + v[3])*(-x[1] + x[3]))/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2));
F[1] = -params[0]*(x[1] - x[3])*(-sqrt(pow(-X[0] + X[2], 2) + pow(-X[1] + X[3], 2)) + sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)))/sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) - params[1]*(-x[1] + x[3])*((-v[0] + v[2])*(-x[0] + x[2]) + (-v[1] + v[3])*(-x[1] + x[3]))/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2));
F[2] = -params[0]*(-x[0] + x[2])*(-sqrt(pow(-X[0] + X[2], 2) + pow(-X[1] + X[3], 2)) + sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)))/sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) + params[1]*(-x[0] + x[2])*((-v[0] + v[2])*(-x[0] + x[2]) + (-v[1] + v[3])*(-x[1] + x[3]))/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2));
F[3] = -params[0]*(-x[1] + x[3])*(-sqrt(pow(-X[0] + X[2], 2) + pow(-X[1] + X[3], 2)) + sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)))/sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) + params[1]*(-x[1] + x[3])*((-v[0] + v[2])*(-x[0] + x[2]) + (-v[1] + v[3])*(-x[1] + x[3]))/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2));

/* Evalaute KX:*/
KX[0] = (-params[0]*pow(x[0] - x[2], 2)*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 7.0L/2.0L)*(sqrt(pow(X[0] - X[2], 2) + pow(X[1] - X[3], 2)) - sqrt(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2))) + params[0]*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 9.0L/2.0L)*(sqrt(pow(X[0] - X[2], 2) + pow(X[1] - X[3], 2)) - sqrt(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2))) - 2*params[1]*pow(x[0] - x[2], 2)*((v[0] - v[2])*(x[0] - x[2]) + (v[1] - v[3])*(x[1] - x[3]))*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 3) + pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 4)*(-params[0]*pow(x[0] - x[2], 2) + params[1]*(v[0] - v[2])*(x[0] - x[2]) + params[1]*((v[0] - v[2])*(x[0] - x[2]) + (v[1] - v[3])*(x[1] - x[3]))))/pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 5);
KX[1] = (x[0] - x[2])*(-params[0]*(x[1] - x[3])*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 3)*(sqrt(pow(X[0] - X[2], 2) + pow(X[1] - X[3], 2)) - sqrt(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2))) - 2*params[1]*(x[1] - x[3])*((v[0] - v[2])*(x[0] - x[2]) + (v[1] - v[3])*(x[1] - x[3]))*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 5.0L/2.0L) + (-params[0]*(x[1] - x[3]) + params[1]*(v[1] - v[3]))*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 7.0L/2.0L))/pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 9.0L/2.0L);
KX[2] = (params[0]*pow(x[0] - x[2], 2)*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 7.0L/2.0L)*(sqrt(pow(X[0] - X[2], 2) + pow(X[1] - X[3], 2)) - sqrt(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2))) - params[0]*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 9.0L/2.0L)*(sqrt(pow(X[0] - X[2], 2) + pow(X[1] - X[3], 2)) - sqrt(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2))) + 2*params[1]*pow(x[0] - x[2], 2)*((v[0] - v[2])*(x[0] - x[2]) + (v[1] - v[3])*(x[1] - x[3]))*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 3) + pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 4)*(params[0]*pow(x[0] - x[2], 2) - params[1]*(v[0] - v[2])*(x[0] - x[2]) - params[1]*((v[0] - v[2])*(x[0] - x[2]) + (v[1] - v[3])*(x[1] - x[3]))))/pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 5);
KX[3] = (x[0] - x[2])*(params[0]*(x[1] - x[3])*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 3)*(sqrt(pow(X[0] - X[2], 2) + pow(X[1] - X[3], 2)) - sqrt(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2))) + 2*params[1]*(x[1] - x[3])*((v[0] - v[2])*(x[0] - x[2]) + (v[1] - v[3])*(x[1] - x[3]))*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 5.0L/2.0L) + (params[0]*(x[1] - x[3]) - params[1]*(v[1] - v[3]))*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 7.0L/2.0L))/pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 9.0L/2.0L);
KX[4] = (x[1] - x[3])*(-params[0]*(x[0] - x[2])*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 3)*(sqrt(pow(X[0] - X[2], 2) + pow(X[1] - X[3], 2)) - sqrt(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2))) - 2*params[1]*(x[0] - x[2])*((v[0] - v[2])*(x[0] - x[2]) + (v[1] - v[3])*(x[1] - x[3]))*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 5.0L/2.0L) + (-params[0]*(x[0] - x[2]) + params[1]*(v[0] - v[2]))*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 7.0L/2.0L))/pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 9.0L/2.0L);
KX[5] = (-params[0]*pow(x[1] - x[3], 2)*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 7.0L/2.0L)*(sqrt(pow(X[0] - X[2], 2) + pow(X[1] - X[3], 2)) - sqrt(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2))) + params[0]*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 9.0L/2.0L)*(sqrt(pow(X[0] - X[2], 2) + pow(X[1] - X[3], 2)) - sqrt(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2))) - 2*params[1]*pow(x[1] - x[3], 2)*((v[0] - v[2])*(x[0] - x[2]) + (v[1] - v[3])*(x[1] - x[3]))*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 3) + pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 4)*(-params[0]*pow(x[1] - x[3], 2) + params[1]*(v[1] - v[3])*(x[1] - x[3]) + params[1]*((v[0] - v[2])*(x[0] - x[2]) + (v[1] - v[3])*(x[1] - x[3]))))/pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 5);
KX[6] = (x[1] - x[3])*(params[0]*(x[0] - x[2])*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 3)*(sqrt(pow(X[0] - X[2], 2) + pow(X[1] - X[3], 2)) - sqrt(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2))) + 2*params[1]*(x[0] - x[2])*((v[0] - v[2])*(x[0] - x[2]) + (v[1] - v[3])*(x[1] - x[3]))*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 5.0L/2.0L) + (params[0]*(x[0] - x[2]) - params[1]*(v[0] - v[2]))*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 7.0L/2.0L))/pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 9.0L/2.0L);
KX[7] = (params[0]*pow(x[1] - x[3], 2)*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 7.0L/2.0L)*(sqrt(pow(X[0] - X[2], 2) + pow(X[1] - X[3], 2)) - sqrt(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2))) - params[0]*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 9.0L/2.0L)*(sqrt(pow(X[0] - X[2], 2) + pow(X[1] - X[3], 2)) - sqrt(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2))) + 2*params[1]*pow(x[1] - x[3], 2)*((v[0] - v[2])*(x[0] - x[2]) + (v[1] - v[3])*(x[1] - x[3]))*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 3) + pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 4)*(params[0]*pow(x[1] - x[3], 2) - params[1]*(v[1] - v[3])*(x[1] - x[3]) - params[1]*((v[0] - v[2])*(x[0] - x[2]) + (v[1] - v[3])*(x[1] - x[3]))))/pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 5);
KX[8] = (params[0]*pow(x[0] - x[2], 2)*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 7.0L/2.0L)*(sqrt(pow(X[0] - X[2], 2) + pow(X[1] - X[3], 2)) - sqrt(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2))) - params[0]*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 9.0L/2.0L)*(sqrt(pow(X[0] - X[2], 2) + pow(X[1] - X[3], 2)) - sqrt(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2))) + 2*params[1]*pow(x[0] - x[2], 2)*((v[0] - v[2])*(x[0] - x[2]) + (v[1] - v[3])*(x[1] - x[3]))*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 3) + pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 4)*(params[0]*pow(x[0] - x[2], 2) - params[1]*(v[0] - v[2])*(x[0] - x[2]) - params[1]*((v[0] - v[2])*(x[0] - x[2]) + (v[1] - v[3])*(x[1] - x[3]))))/pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 5);
KX[9] = (x[0] - x[2])*(params[0]*(x[1] - x[3])*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 3)*(sqrt(pow(X[0] - X[2], 2) + pow(X[1] - X[3], 2)) - sqrt(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2))) + 2*params[1]*(x[1] - x[3])*((v[0] - v[2])*(x[0] - x[2]) + (v[1] - v[3])*(x[1] - x[3]))*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 5.0L/2.0L) + (params[0]*(x[1] - x[3]) - params[1]*(v[1] - v[3]))*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 7.0L/2.0L))/pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 9.0L/2.0L);
KX[10] = (-params[0]*pow(x[0] - x[2], 2)*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 7.0L/2.0L)*(sqrt(pow(X[0] - X[2], 2) + pow(X[1] - X[3], 2)) - sqrt(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2))) + params[0]*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 9.0L/2.0L)*(sqrt(pow(X[0] - X[2], 2) + pow(X[1] - X[3], 2)) - sqrt(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2))) - 2*params[1]*pow(x[0] - x[2], 2)*((v[0] - v[2])*(x[0] - x[2]) + (v[1] - v[3])*(x[1] - x[3]))*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 3) + pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 4)*(-params[0]*pow(x[0] - x[2], 2) + params[1]*(v[0] - v[2])*(x[0] - x[2]) + params[1]*((v[0] - v[2])*(x[0] - x[2]) + (v[1] - v[3])*(x[1] - x[3]))))/pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 5);
KX[11] = (x[0] - x[2])*(-params[0]*(x[1] - x[3])*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 3)*(sqrt(pow(X[0] - X[2], 2) + pow(X[1] - X[3], 2)) - sqrt(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2))) - 2*params[1]*(x[1] - x[3])*((v[0] - v[2])*(x[0] - x[2]) + (v[1] - v[3])*(x[1] - x[3]))*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 5.0L/2.0L) + (-params[0]*(x[1] - x[3]) + params[1]*(v[1] - v[3]))*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 7.0L/2.0L))/pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 9.0L/2.0L);
KX[12] = (x[1] - x[3])*(params[0]*(x[0] - x[2])*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 3)*(sqrt(pow(X[0] - X[2], 2) + pow(X[1] - X[3], 2)) - sqrt(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2))) + 2*params[1]*(x[0] - x[2])*((v[0] - v[2])*(x[0] - x[2]) + (v[1] - v[3])*(x[1] - x[3]))*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 5.0L/2.0L) + (params[0]*(x[0] - x[2]) - params[1]*(v[0] - v[2]))*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 7.0L/2.0L))/pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 9.0L/2.0L);
KX[13] = (params[0]*pow(x[1] - x[3], 2)*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 7.0L/2.0L)*(sqrt(pow(X[0] - X[2], 2) + pow(X[1] - X[3], 2)) - sqrt(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2))) - params[0]*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 9.0L/2.0L)*(sqrt(pow(X[0] - X[2], 2) + pow(X[1] - X[3], 2)) - sqrt(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2))) + 2*params[1]*pow(x[1] - x[3], 2)*((v[0] - v[2])*(x[0] - x[2]) + (v[1] - v[3])*(x[1] - x[3]))*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 3) + pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 4)*(params[0]*pow(x[1] - x[3], 2) - params[1]*(v[1] - v[3])*(x[1] - x[3]) - params[1]*((v[0] - v[2])*(x[0] - x[2]) + (v[1] - v[3])*(x[1] - x[3]))))/pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 5);
KX[14] = (x[1] - x[3])*(-params[0]*(x[0] - x[2])*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 3)*(sqrt(pow(X[0] - X[2], 2) + pow(X[1] - X[3], 2)) - sqrt(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2))) - 2*params[1]*(x[0] - x[2])*((v[0] - v[2])*(x[0] - x[2]) + (v[1] - v[3])*(x[1] - x[3]))*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 5.0L/2.0L) + (-params[0]*(x[0] - x[2]) + params[1]*(v[0] - v[2]))*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 7.0L/2.0L))/pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 9.0L/2.0L);
KX[15] = (-params[0]*pow(x[1] - x[3], 2)*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 7.0L/2.0L)*(sqrt(pow(X[0] - X[2], 2) + pow(X[1] - X[3], 2)) - sqrt(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2))) + params[0]*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 9.0L/2.0L)*(sqrt(pow(X[0] - X[2], 2) + pow(X[1] - X[3], 2)) - sqrt(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2))) - 2*params[1]*pow(x[1] - x[3], 2)*((v[0] - v[2])*(x[0] - x[2]) + (v[1] - v[3])*(x[1] - x[3]))*pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 3) + pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 4)*(-params[0]*pow(x[1] - x[3], 2) + params[1]*(v[1] - v[3])*(x[1] - x[3]) + params[1]*((v[0] - v[2])*(x[0] - x[2]) + (v[1] - v[3])*(x[1] - x[3]))))/pow(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2), 5);

/* Evalaute KV:*/
KV[0] = params[1]*pow(x[0] - x[2], 2)/(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2));
KV[1] = params[1]*(x[0] - x[2])*(x[1] - x[3])/(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2));
KV[2] = -params[1]*pow(x[0] - x[2], 2)/(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2));
KV[3] = -params[1]*(x[0] - x[2])*(x[1] - x[3])/(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2));
KV[4] = params[1]*(x[0] - x[2])*(x[1] - x[3])/(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2));
KV[5] = params[1]*pow(x[1] - x[3], 2)/(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2));
KV[6] = -params[1]*(x[0] - x[2])*(x[1] - x[3])/(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2));
KV[7] = -params[1]*pow(x[1] - x[3], 2)/(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2));
KV[8] = -params[1]*pow(x[0] - x[2], 2)/(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2));
KV[9] = -params[1]*(x[0] - x[2])*(x[1] - x[3])/(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2));
KV[10] = params[1]*pow(x[0] - x[2], 2)/(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2));
KV[11] = params[1]*(x[0] - x[2])*(x[1] - x[3])/(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2));
KV[12] = -params[1]*(x[0] - x[2])*(x[1] - x[3])/(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2));
KV[13] = -params[1]*pow(x[1] - x[3], 2)/(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2));
KV[14] = params[1]*(x[0] - x[2])*(x[1] - x[3])/(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2));
KV[15] = params[1]*pow(x[1] - x[3], 2)/(pow(x[0] - x[2], 2) + pow(x[1] - x[3], 2));
}

void particle_kernel_calc(real_t* in,real_t* out) {
  eval_peri(
	    in,
	    in+4,
	    in+8,
	    in+12,

	    out,
	    out+4,
	    out+4+16
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


void eval_peri3d(/*Inputs:*/
real_t * x,
real_t * v,
real_t * X,
real_t * params
,
/*Outputs:*/
double * F)
{
/* Evalaute F:*/
F[0] = -1.0*params[0]*(x[0] - x[3])*(-sqrt(pow(-X[0] + X[3], 2) + pow(-X[1] + X[4], 2) + pow(-X[2] + X[5], 2)) + sqrt(pow(-x[0] + x[3], 2) + pow(-x[1] + x[4], 2) + pow(-x[2] + x[5], 2)))/sqrt(pow(-x[0] + x[3], 2) + pow(-x[1] + x[4], 2) + pow(-x[2] + x[5], 2)) - params[1]*(-x[0] + x[3])*((-v[0] + v[3])*(-x[0] + x[3]) + (-v[1] + v[4])*(-x[1] + x[4]) + (-v[2] + v[5])*(-x[2] + x[5]))/(pow(-x[0] + x[3], 2) + pow(-x[1] + x[4], 2) + pow(-x[2] + x[5], 2));
F[1] = -1.0*params[0]*(x[1] - x[4])*(-sqrt(pow(-X[0] + X[3], 2) + pow(-X[1] + X[4], 2) + pow(-X[2] + X[5], 2)) + sqrt(pow(-x[0] + x[3], 2) + pow(-x[1] + x[4], 2) + pow(-x[2] + x[5], 2)))/sqrt(pow(-x[0] + x[3], 2) + pow(-x[1] + x[4], 2) + pow(-x[2] + x[5], 2)) - params[1]*(-x[1] + x[4])*((-v[0] + v[3])*(-x[0] + x[3]) + (-v[1] + v[4])*(-x[1] + x[4]) + (-v[2] + v[5])*(-x[2] + x[5]))/(pow(-x[0] + x[3], 2) + pow(-x[1] + x[4], 2) + pow(-x[2] + x[5], 2));
F[2] = -1.0*params[0]*(x[2] - x[5])*(-sqrt(pow(-X[0] + X[3], 2) + pow(-X[1] + X[4], 2) + pow(-X[2] + X[5], 2)) + sqrt(pow(-x[0] + x[3], 2) + pow(-x[1] + x[4], 2) + pow(-x[2] + x[5], 2)))/sqrt(pow(-x[0] + x[3], 2) + pow(-x[1] + x[4], 2) + pow(-x[2] + x[5], 2)) - params[1]*(-x[2] + x[5])*((-v[0] + v[3])*(-x[0] + x[3]) + (-v[1] + v[4])*(-x[1] + x[4]) + (-v[2] + v[5])*(-x[2] + x[5]))/(pow(-x[0] + x[3], 2) + pow(-x[1] + x[4], 2) + pow(-x[2] + x[5], 2));
F[3] = -1.0*params[0]*(-x[0] + x[3])*(-sqrt(pow(-X[0] + X[3], 2) + pow(-X[1] + X[4], 2) + pow(-X[2] + X[5], 2)) + sqrt(pow(-x[0] + x[3], 2) + pow(-x[1] + x[4], 2) + pow(-x[2] + x[5], 2)))/sqrt(pow(-x[0] + x[3], 2) + pow(-x[1] + x[4], 2) + pow(-x[2] + x[5], 2)) + params[1]*(-x[0] + x[3])*((-v[0] + v[3])*(-x[0] + x[3]) + (-v[1] + v[4])*(-x[1] + x[4]) + (-v[2] + v[5])*(-x[2] + x[5]))/(pow(-x[0] + x[3], 2) + pow(-x[1] + x[4], 2) + pow(-x[2] + x[5], 2));
F[4] = -1.0*params[0]*(-x[1] + x[4])*(-sqrt(pow(-X[0] + X[3], 2) + pow(-X[1] + X[4], 2) + pow(-X[2] + X[5], 2)) + sqrt(pow(-x[0] + x[3], 2) + pow(-x[1] + x[4], 2) + pow(-x[2] + x[5], 2)))/sqrt(pow(-x[0] + x[3], 2) + pow(-x[1] + x[4], 2) + pow(-x[2] + x[5], 2)) + params[1]*(-x[1] + x[4])*((-v[0] + v[3])*(-x[0] + x[3]) + (-v[1] + v[4])*(-x[1] + x[4]) + (-v[2] + v[5])*(-x[2] + x[5]))/(pow(-x[0] + x[3], 2) + pow(-x[1] + x[4], 2) + pow(-x[2] + x[5], 2));
F[5] = -1.0*params[0]*(-x[2] + x[5])*(-sqrt(pow(-X[0] + X[3], 2) + pow(-X[1] + X[4], 2) + pow(-X[2] + X[5], 2)) + sqrt(pow(-x[0] + x[3], 2) + pow(-x[1] + x[4], 2) + pow(-x[2] + x[5], 2)))/sqrt(pow(-x[0] + x[3], 2) + pow(-x[1] + x[4], 2) + pow(-x[2] + x[5], 2)) + params[1]*(-x[2] + x[5])*((-v[0] + v[3])*(-x[0] + x[3]) + (-v[1] + v[4])*(-x[1] + x[4]) + (-v[2] + v[5])*(-x[2] + x[5]))/(pow(-x[0] + x[3], 2) + pow(-x[1] + x[4], 2) + pow(-x[2] + x[5], 2));
}



void particle3d_kernel_calc(real_t* in,real_t* out) {
  eval_peri3d(
	    in,
	    in+6,
	    in+12,
	    in+18,

	    out
	    
	    );
	    
}

kernel_t particle3d_kernel_strct = {
  .ninp = 4,
  .inp = { {LOC_NODE,3},
	   {LOC_NODE,3},
	   {LOC_NODE,3},
	   {LOC_GLOBAL,2},
	   {LOC_GLOBAL,0} },
  .noutp = 1,
  .outp = { {LOC_NODE,3},
	    {LOC_GLOBAL,0},
	    {LOC_GLOBAL,0},
	    {LOC_GLOBAL,0},
	    {LOC_GLOBAL,0} },
  .assem = particle3d_kernel_calc
};
