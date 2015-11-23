#include "sample_peri.h"
#include "math.h"

void peri_eval(/*Inputs:*/
real_t * x,
real_t * v,
real_t * X,
real_t * alpha,
real_t * params
,
/*Outputs:*/
double * F,
double * KX,
double * KV)
{

/*Declarations:*/


/*Intermediates:*/


/*Final calcs:*/
/* Evalaute F:*/
F[0] = -1.0*alpha[0]*params[0]*(x[0] - x[2])*(-sqrt(pow(-X[0] + X[2], 2) + pow(-X[1] + X[3], 2)) + sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)))/sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) - alpha[0]*params[1]*(-x[0] + x[2])*((-v[0] + v[2])*(-x[0] + x[2]) + (-v[1] + v[3])*(-x[1] + x[3]))/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2));
F[1] = -1.0*alpha[0]*params[0]*(x[1] - x[3])*(-sqrt(pow(-X[0] + X[2], 2) + pow(-X[1] + X[3], 2)) + sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)))/sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) - alpha[0]*params[1]*(-x[1] + x[3])*((-v[0] + v[2])*(-x[0] + x[2]) + (-v[1] + v[3])*(-x[1] + x[3]))/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2));
F[2] = -1.0*alpha[0]*params[0]*(-x[0] + x[2])*(-sqrt(pow(-X[0] + X[2], 2) + pow(-X[1] + X[3], 2)) + sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)))/sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) + alpha[0]*params[1]*(-x[0] + x[2])*((-v[0] + v[2])*(-x[0] + x[2]) + (-v[1] + v[3])*(-x[1] + x[3]))/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2));
F[3] = -1.0*alpha[0]*params[0]*(-x[1] + x[3])*(-sqrt(pow(-X[0] + X[2], 2) + pow(-X[1] + X[3], 2)) + sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)))/sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) + alpha[0]*params[1]*(-x[1] + x[3])*((-v[0] + v[2])*(-x[0] + x[2]) + (-v[1] + v[3])*(-x[1] + x[3]))/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2));

/* Evalaute KX:*/
KX[0] = -1.0*alpha[0]*params[0]*(-x[0] + x[2])*(x[0] - x[2])*(-sqrt(pow(-X[0] + X[2], 2) + pow(-X[1] + X[3], 2)) + sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)))/pow(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2), 3.0L/2.0L) - 1.0*alpha[0]*params[0]*pow(x[0] - x[2], 2)/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) - 1.0*alpha[0]*params[0]*(-sqrt(pow(-X[0] + X[2], 2) + pow(-X[1] + X[3], 2)) + sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)))/sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) - alpha[0]*params[1]*(v[0] - v[2])*(-x[0] + x[2])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) - alpha[0]*params[1]*(-2*x[0] + 2*x[2])*(-x[0] + x[2])*((-v[0] + v[2])*(-x[0] + x[2]) + (-v[1] + v[3])*(-x[1] + x[3]))/pow(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2), 2) + alpha[0]*params[1]*((-v[0] + v[2])*(-x[0] + x[2]) + (-v[1] + v[3])*(-x[1] + x[3]))/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2));
KX[1] = -1.0*alpha[0]*params[0]*(x[0] - x[2])*(-x[1] + x[3])*(-sqrt(pow(-X[0] + X[2], 2) + pow(-X[1] + X[3], 2)) + sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)))/pow(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2), 3.0L/2.0L) - 1.0*alpha[0]*params[0]*(x[0] - x[2])*(x[1] - x[3])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) - alpha[0]*params[1]*(v[1] - v[3])*(-x[0] + x[2])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) - alpha[0]*params[1]*(-x[0] + x[2])*(-2*x[1] + 2*x[3])*((-v[0] + v[2])*(-x[0] + x[2]) + (-v[1] + v[3])*(-x[1] + x[3]))/pow(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2), 2);
KX[2] = -1.0*alpha[0]*params[0]*(-x[0] + x[2])*(x[0] - x[2])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) - 1.0*alpha[0]*params[0]*pow(x[0] - x[2], 2)*(-sqrt(pow(-X[0] + X[2], 2) + pow(-X[1] + X[3], 2)) + sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)))/pow(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2), 3.0L/2.0L) + 1.0*alpha[0]*params[0]*(-sqrt(pow(-X[0] + X[2], 2) + pow(-X[1] + X[3], 2)) + sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)))/sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) - alpha[0]*params[1]*(-v[0] + v[2])*(-x[0] + x[2])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) - alpha[0]*params[1]*(-x[0] + x[2])*(2*x[0] - 2*x[2])*((-v[0] + v[2])*(-x[0] + x[2]) + (-v[1] + v[3])*(-x[1] + x[3]))/pow(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2), 2) - alpha[0]*params[1]*((-v[0] + v[2])*(-x[0] + x[2]) + (-v[1] + v[3])*(-x[1] + x[3]))/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2));
KX[3] = -1.0*alpha[0]*params[0]*(x[0] - x[2])*(-x[1] + x[3])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) - 1.0*alpha[0]*params[0]*(x[0] - x[2])*(x[1] - x[3])*(-sqrt(pow(-X[0] + X[2], 2) + pow(-X[1] + X[3], 2)) + sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)))/pow(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2), 3.0L/2.0L) - alpha[0]*params[1]*(-v[1] + v[3])*(-x[0] + x[2])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) - alpha[0]*params[1]*(-x[0] + x[2])*(2*x[1] - 2*x[3])*((-v[0] + v[2])*(-x[0] + x[2]) + (-v[1] + v[3])*(-x[1] + x[3]))/pow(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2), 2);
KX[4] = -1.0*alpha[0]*params[0]*(-x[0] + x[2])*(x[1] - x[3])*(-sqrt(pow(-X[0] + X[2], 2) + pow(-X[1] + X[3], 2)) + sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)))/pow(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2), 3.0L/2.0L) - 1.0*alpha[0]*params[0]*(x[0] - x[2])*(x[1] - x[3])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) - alpha[0]*params[1]*(v[0] - v[2])*(-x[1] + x[3])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) - alpha[0]*params[1]*(-2*x[0] + 2*x[2])*(-x[1] + x[3])*((-v[0] + v[2])*(-x[0] + x[2]) + (-v[1] + v[3])*(-x[1] + x[3]))/pow(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2), 2);
KX[5] = -1.0*alpha[0]*params[0]*(-x[1] + x[3])*(x[1] - x[3])*(-sqrt(pow(-X[0] + X[2], 2) + pow(-X[1] + X[3], 2)) + sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)))/pow(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2), 3.0L/2.0L) - 1.0*alpha[0]*params[0]*pow(x[1] - x[3], 2)/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) - 1.0*alpha[0]*params[0]*(-sqrt(pow(-X[0] + X[2], 2) + pow(-X[1] + X[3], 2)) + sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)))/sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) - alpha[0]*params[1]*(v[1] - v[3])*(-x[1] + x[3])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) - alpha[0]*params[1]*(-2*x[1] + 2*x[3])*(-x[1] + x[3])*((-v[0] + v[2])*(-x[0] + x[2]) + (-v[1] + v[3])*(-x[1] + x[3]))/pow(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2), 2) + alpha[0]*params[1]*((-v[0] + v[2])*(-x[0] + x[2]) + (-v[1] + v[3])*(-x[1] + x[3]))/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2));
KX[6] = -1.0*alpha[0]*params[0]*(-x[0] + x[2])*(x[1] - x[3])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) - 1.0*alpha[0]*params[0]*(x[0] - x[2])*(x[1] - x[3])*(-sqrt(pow(-X[0] + X[2], 2) + pow(-X[1] + X[3], 2)) + sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)))/pow(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2), 3.0L/2.0L) - alpha[0]*params[1]*(-v[0] + v[2])*(-x[1] + x[3])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) - alpha[0]*params[1]*(2*x[0] - 2*x[2])*(-x[1] + x[3])*((-v[0] + v[2])*(-x[0] + x[2]) + (-v[1] + v[3])*(-x[1] + x[3]))/pow(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2), 2);
KX[7] = -1.0*alpha[0]*params[0]*(-x[1] + x[3])*(x[1] - x[3])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) - 1.0*alpha[0]*params[0]*pow(x[1] - x[3], 2)*(-sqrt(pow(-X[0] + X[2], 2) + pow(-X[1] + X[3], 2)) + sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)))/pow(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2), 3.0L/2.0L) + 1.0*alpha[0]*params[0]*(-sqrt(pow(-X[0] + X[2], 2) + pow(-X[1] + X[3], 2)) + sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)))/sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) - alpha[0]*params[1]*(-v[1] + v[3])*(-x[1] + x[3])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) - alpha[0]*params[1]*(-x[1] + x[3])*(2*x[1] - 2*x[3])*((-v[0] + v[2])*(-x[0] + x[2]) + (-v[1] + v[3])*(-x[1] + x[3]))/pow(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2), 2) - alpha[0]*params[1]*((-v[0] + v[2])*(-x[0] + x[2]) + (-v[1] + v[3])*(-x[1] + x[3]))/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2));
KX[8] = -1.0*alpha[0]*params[0]*pow(-x[0] + x[2], 2)*(-sqrt(pow(-X[0] + X[2], 2) + pow(-X[1] + X[3], 2)) + sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)))/pow(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2), 3.0L/2.0L) - 1.0*alpha[0]*params[0]*(-x[0] + x[2])*(x[0] - x[2])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) + 1.0*alpha[0]*params[0]*(-sqrt(pow(-X[0] + X[2], 2) + pow(-X[1] + X[3], 2)) + sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)))/sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) + alpha[0]*params[1]*(v[0] - v[2])*(-x[0] + x[2])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) + alpha[0]*params[1]*(-2*x[0] + 2*x[2])*(-x[0] + x[2])*((-v[0] + v[2])*(-x[0] + x[2]) + (-v[1] + v[3])*(-x[1] + x[3]))/pow(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2), 2) - alpha[0]*params[1]*((-v[0] + v[2])*(-x[0] + x[2]) + (-v[1] + v[3])*(-x[1] + x[3]))/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2));
KX[9] = -1.0*alpha[0]*params[0]*(-x[0] + x[2])*(-x[1] + x[3])*(-sqrt(pow(-X[0] + X[2], 2) + pow(-X[1] + X[3], 2)) + sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)))/pow(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2), 3.0L/2.0L) - 1.0*alpha[0]*params[0]*(-x[0] + x[2])*(x[1] - x[3])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) + alpha[0]*params[1]*(v[1] - v[3])*(-x[0] + x[2])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) + alpha[0]*params[1]*(-x[0] + x[2])*(-2*x[1] + 2*x[3])*((-v[0] + v[2])*(-x[0] + x[2]) + (-v[1] + v[3])*(-x[1] + x[3]))/pow(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2), 2);
KX[10] = -1.0*alpha[0]*params[0]*pow(-x[0] + x[2], 2)/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) - 1.0*alpha[0]*params[0]*(-x[0] + x[2])*(x[0] - x[2])*(-sqrt(pow(-X[0] + X[2], 2) + pow(-X[1] + X[3], 2)) + sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)))/pow(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2), 3.0L/2.0L) - 1.0*alpha[0]*params[0]*(-sqrt(pow(-X[0] + X[2], 2) + pow(-X[1] + X[3], 2)) + sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)))/sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) + alpha[0]*params[1]*(-v[0] + v[2])*(-x[0] + x[2])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) + alpha[0]*params[1]*(-x[0] + x[2])*(2*x[0] - 2*x[2])*((-v[0] + v[2])*(-x[0] + x[2]) + (-v[1] + v[3])*(-x[1] + x[3]))/pow(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2), 2) + alpha[0]*params[1]*((-v[0] + v[2])*(-x[0] + x[2]) + (-v[1] + v[3])*(-x[1] + x[3]))/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2));
KX[11] = -1.0*alpha[0]*params[0]*(-x[0] + x[2])*(-x[1] + x[3])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) - 1.0*alpha[0]*params[0]*(-x[0] + x[2])*(x[1] - x[3])*(-sqrt(pow(-X[0] + X[2], 2) + pow(-X[1] + X[3], 2)) + sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)))/pow(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2), 3.0L/2.0L) + alpha[0]*params[1]*(-v[1] + v[3])*(-x[0] + x[2])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) + alpha[0]*params[1]*(-x[0] + x[2])*(2*x[1] - 2*x[3])*((-v[0] + v[2])*(-x[0] + x[2]) + (-v[1] + v[3])*(-x[1] + x[3]))/pow(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2), 2);
KX[12] = -1.0*alpha[0]*params[0]*(-x[0] + x[2])*(-x[1] + x[3])*(-sqrt(pow(-X[0] + X[2], 2) + pow(-X[1] + X[3], 2)) + sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)))/pow(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2), 3.0L/2.0L) - 1.0*alpha[0]*params[0]*(x[0] - x[2])*(-x[1] + x[3])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) + alpha[0]*params[1]*(v[0] - v[2])*(-x[1] + x[3])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) + alpha[0]*params[1]*(-2*x[0] + 2*x[2])*(-x[1] + x[3])*((-v[0] + v[2])*(-x[0] + x[2]) + (-v[1] + v[3])*(-x[1] + x[3]))/pow(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2), 2);
KX[13] = -1.0*alpha[0]*params[0]*pow(-x[1] + x[3], 2)*(-sqrt(pow(-X[0] + X[2], 2) + pow(-X[1] + X[3], 2)) + sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)))/pow(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2), 3.0L/2.0L) - 1.0*alpha[0]*params[0]*(-x[1] + x[3])*(x[1] - x[3])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) + 1.0*alpha[0]*params[0]*(-sqrt(pow(-X[0] + X[2], 2) + pow(-X[1] + X[3], 2)) + sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)))/sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) + alpha[0]*params[1]*(v[1] - v[3])*(-x[1] + x[3])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) + alpha[0]*params[1]*(-2*x[1] + 2*x[3])*(-x[1] + x[3])*((-v[0] + v[2])*(-x[0] + x[2]) + (-v[1] + v[3])*(-x[1] + x[3]))/pow(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2), 2) - alpha[0]*params[1]*((-v[0] + v[2])*(-x[0] + x[2]) + (-v[1] + v[3])*(-x[1] + x[3]))/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2));
KX[14] = -1.0*alpha[0]*params[0]*(-x[0] + x[2])*(-x[1] + x[3])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) - 1.0*alpha[0]*params[0]*(x[0] - x[2])*(-x[1] + x[3])*(-sqrt(pow(-X[0] + X[2], 2) + pow(-X[1] + X[3], 2)) + sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)))/pow(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2), 3.0L/2.0L) + alpha[0]*params[1]*(-v[0] + v[2])*(-x[1] + x[3])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) + alpha[0]*params[1]*(2*x[0] - 2*x[2])*(-x[1] + x[3])*((-v[0] + v[2])*(-x[0] + x[2]) + (-v[1] + v[3])*(-x[1] + x[3]))/pow(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2), 2);
KX[15] = -1.0*alpha[0]*params[0]*pow(-x[1] + x[3], 2)/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) - 1.0*alpha[0]*params[0]*(-x[1] + x[3])*(x[1] - x[3])*(-sqrt(pow(-X[0] + X[2], 2) + pow(-X[1] + X[3], 2)) + sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)))/pow(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2), 3.0L/2.0L) - 1.0*alpha[0]*params[0]*(-sqrt(pow(-X[0] + X[2], 2) + pow(-X[1] + X[3], 2)) + sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)))/sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) + alpha[0]*params[1]*(-v[1] + v[3])*(-x[1] + x[3])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)) + alpha[0]*params[1]*(-x[1] + x[3])*(2*x[1] - 2*x[3])*((-v[0] + v[2])*(-x[0] + x[2]) + (-v[1] + v[3])*(-x[1] + x[3]))/pow(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2), 2) + alpha[0]*params[1]*((-v[0] + v[2])*(-x[0] + x[2]) + (-v[1] + v[3])*(-x[1] + x[3]))/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2));

/* Evalaute KV:*/
KV[0] = -alpha[0]*params[1]*(-x[0] + x[2])*(x[0] - x[2])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2));
KV[1] = -alpha[0]*params[1]*(-x[0] + x[2])*(x[1] - x[3])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2));
KV[2] = -alpha[0]*params[1]*pow(-x[0] + x[2], 2)/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2));
KV[3] = -alpha[0]*params[1]*(-x[0] + x[2])*(-x[1] + x[3])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2));
KV[4] = -alpha[0]*params[1]*(x[0] - x[2])*(-x[1] + x[3])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2));
KV[5] = -alpha[0]*params[1]*(-x[1] + x[3])*(x[1] - x[3])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2));
KV[6] = -alpha[0]*params[1]*(-x[0] + x[2])*(-x[1] + x[3])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2));
KV[7] = -alpha[0]*params[1]*pow(-x[1] + x[3], 2)/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2));
KV[8] = alpha[0]*params[1]*(-x[0] + x[2])*(x[0] - x[2])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2));
KV[9] = alpha[0]*params[1]*(-x[0] + x[2])*(x[1] - x[3])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2));
KV[10] = alpha[0]*params[1]*pow(-x[0] + x[2], 2)/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2));
KV[11] = alpha[0]*params[1]*(-x[0] + x[2])*(-x[1] + x[3])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2));
KV[12] = alpha[0]*params[1]*(x[0] - x[2])*(-x[1] + x[3])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2));
KV[13] = alpha[0]*params[1]*(-x[1] + x[3])*(x[1] - x[3])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2));
KV[14] = alpha[0]*params[1]*(-x[0] + x[2])*(-x[1] + x[3])/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2));
KV[15] = alpha[0]*params[1]*pow(-x[1] + x[3], 2)/(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2));

}



void peri_eval_wr(int ninp, real_t * in, real_t * out) {
   peri_eval(in+0, in+4, in+8, in+12, in+13, 
out+0, out+4, out+8);
}


kernel_t kern_peri = {
  .nmap = 3,
  .maps = {
    { .dim = 2, .v_start = 0, .v_end = 2 },
    { .dim = 1, .v_start = 2, .v_end = 3 },
    { .dim = 2, .v_start= -1, .v_end = -1 }
  },

  .ninp = 5,
  .inp = {
    { .field_number = 0, .map_num = 0, .name = "x" },
    { .field_number = 1, .map_num = 0, .name = "v" },
    { .field_number = 2, .map_num = 0, .name = "X" },
    { .field_number = 3, .map_num = 1, .name = "alpha" },
    { .field_number = 4, .map_num = 2, .name = "params" }
  },

  .noutp = 3,
  .outp = {
    {.rank = 1, .nmap = 1, .map_nums = { 0 }, .name = "F" },
    {.rank = 2, .nmap = 1, .map_nums = { 0 }, .name = "KX" },
    {.rank = 2, .nmap = 1, .map_nums = { 0 }, .name = "KV" }
  },

  .eval = peri_eval_wr,
  .name = "peri"
};



#if 0
kernel_t kern_peri = {
.ninp = 5,
.inp = {
{ .field_number = 0, .dim = 2, .v_start = 0, .v_end = 2, .name = "x" },
{ .field_number = 1, .dim = 2, .v_start = 0, .v_end = 2, .name = "v" },
{ .field_number = 2, .dim = 2, .v_start = 0, .v_end = 2, .name = "X" },
{ .field_number = 3, .dim = 1, .v_start = 2, .v_end = 3, .name = "alpha" },
{ .field_number = 4, .dim = 1, .v_start = 2, .v_end = 3, .name = "params" }
},

.noutp = 3,
.outps = {
{ .rank = 1, .ndof = 1, .dofs = { { .field_number = 0, .dim = 2, .v_start = 0, .v_end = 2, .name = "F" } } },
{ .rank = 2, .ndof = 1, .dofs = { { .field_number = 0, .dim = 2, .v_start = 0, .v_end = 2, .name = "KX" } } },
{ .rank = 2, .ndof = 1, .dofs = { { .field_number = 0, .dim = 2, .v_start = 0, .v_end = 2, .name = "KV" } } },

},

.eval=peri_eval_wr,
.name="peri"
};
#endif
