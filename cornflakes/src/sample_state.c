#include "sample_state.h"
#include "math.h"
void state_eval(int l_edge,
                /*Inputs:*/
                const real_t * restrict y,
                const real_t * restrict v,
                const real_t * restrict x,
                const real_t * restrict param,
                /*Outputs:*/
                real_t * restrict F)
{
    real_t theta[1];
    {
        int i;
        for(i=0; i<1; i++) theta[i]=0.0;
    }
    int i;
    for(i=1; i<l_edge; i++) {
        /* Evaluation of theta */
        theta[0] += -1.0*param[0]*((x[0] - x[2*i])*(y[0] - y[2*i]) + (x[1] - x[2*i + 1])*(y[1] - y[2*i + 1]))*(sqrt(pow(x[0] - x[2*i], 2) + pow(x[1] - x[2*i + 1], 2)) - sqrt(pow(y[0] - y[2*i], 2) + pow(y[1] - y[2*i + 1], 2)))/((pow(x[0] - x[2*i], 2) + pow(x[1] - x[2*i + 1], 2))*sqrt(pow(y[0] - y[2*i], 2) + pow(y[1] - y[2*i + 1], 2)));
    }

    for(i=1; i<l_edge; i++) {
        /* Evaluation of F */
        F[0] += 2.0*pow(param[0], -9.0)*(y[0] - y[2*i])*(3.0*M_E*M_PI*pow(param[0], 6)*sqrt(pow(x[0] - x[2*i], 2) + pow(x[1] - x[2*i + 1], 2))*sqrt(pow(y[0] - y[2*i], 2) + pow(y[1] - y[2*i + 1], 2))*(sqrt(pow(x[0] - x[2*i], 2) + pow(x[1] - x[2*i + 1], 2)) - sqrt(pow(y[0] - y[2*i], 2) + pow(y[1] - y[2*i + 1], 2)))*(12.0*pow(param[2], 2) + 6.0*param[2] - 6.0) + 4.0*pow(param[0], 4.0)*param[1]*theta[0]*(param[2] + 1)*(7.0*param[2] - 2.0)*((x[0] - x[2*i])*(y[0] - y[2*i]) + (x[1] - x[2*i + 1])*(y[1] - y[2*i + 1])))/(pow(M_PI, 2)*(param[2] + 1)*(pow(x[0] - x[2*i], 2) + pow(x[1] - x[2*i + 1], 2))*(pow(y[0] - y[2*i], 2) + pow(y[1] - y[2*i + 1], 2))*(12.0*pow(param[2], 2) + 6.0*param[2] - 6.0));
        F[1] += 2.0*pow(param[0], -9.0)*(y[1] - y[2*i + 1])*(3.0*M_E*M_PI*pow(param[0], 6)*sqrt(pow(x[0] - x[2*i], 2) + pow(x[1] - x[2*i + 1], 2))*sqrt(pow(y[0] - y[2*i], 2) + pow(y[1] - y[2*i + 1], 2))*(sqrt(pow(x[0] - x[2*i], 2) + pow(x[1] - x[2*i + 1], 2)) - sqrt(pow(y[0] - y[2*i], 2) + pow(y[1] - y[2*i + 1], 2)))*(12.0*pow(param[2], 2) + 6.0*param[2] - 6.0) + 4.0*pow(param[0], 4.0)*param[1]*theta[0]*(param[2] + 1)*(7.0*param[2] - 2.0)*((x[0] - x[2*i])*(y[0] - y[2*i]) + (x[1] - x[2*i + 1])*(y[1] - y[2*i + 1])))/(pow(M_PI, 2)*(param[2] + 1)*(pow(x[0] - x[2*i], 2) + pow(x[1] - x[2*i + 1], 2))*(pow(y[0] - y[2*i], 2) + pow(y[1] - y[2*i + 1], 2))*(12.0*pow(param[2], 2) + 6.0*param[2] - 6.0));
        /* Evaluation of F */
        F[2*i] -= 2.0*pow(param[0], -9.0)*(y[0] - y[2*i])*(3.0*M_E*M_PI*pow(param[0], 6)*sqrt(pow(x[0] - x[2*i], 2) + pow(x[1] - x[2*i + 1], 2))*sqrt(pow(y[0] - y[2*i], 2) + pow(y[1] - y[2*i + 1], 2))*(sqrt(pow(x[0] - x[2*i], 2) + pow(x[1] - x[2*i + 1], 2)) - sqrt(pow(y[0] - y[2*i], 2) + pow(y[1] - y[2*i + 1], 2)))*(12.0*pow(param[2], 2) + 6.0*param[2] - 6.0) + 4.0*pow(param[0], 4.0)*param[1]*theta[0]*(param[2] + 1)*(7.0*param[2] - 2.0)*((x[0] - x[2*i])*(y[0] - y[2*i]) + (x[1] - x[2*i + 1])*(y[1] - y[2*i + 1])))/(pow(M_PI, 2)*(param[2] + 1)*(pow(x[0] - x[2*i], 2) + pow(x[1] - x[2*i + 1], 2))*(pow(y[0] - y[2*i], 2) + pow(y[1] - y[2*i + 1], 2))*(12.0*pow(param[2], 2) + 6.0*param[2] - 6.0));
        F[2*i + 1] -= 2.0*pow(param[0], -9.0)*(y[1] - y[2*i + 1])*(3.0*M_E*M_PI*pow(param[0], 6)*sqrt(pow(x[0] - x[2*i], 2) + pow(x[1] - x[2*i + 1], 2))*sqrt(pow(y[0] - y[2*i], 2) + pow(y[1] - y[2*i + 1], 2))*(sqrt(pow(x[0] - x[2*i], 2) + pow(x[1] - x[2*i + 1], 2)) - sqrt(pow(y[0] - y[2*i], 2) + pow(y[1] - y[2*i + 1], 2)))*(12.0*pow(param[2], 2) + 6.0*param[2] - 6.0) + 4.0*pow(param[0], 4.0)*param[1]*theta[0]*(param[2] + 1)*(7.0*param[2] - 2.0)*((x[0] - x[2*i])*(y[0] - y[2*i]) + (x[1] - x[2*i + 1])*(y[1] - y[2*i + 1])))/(pow(M_PI, 2)*(param[2] + 1)*(pow(x[0] - x[2*i], 2) + pow(x[1] - x[2*i + 1], 2))*(pow(y[0] - y[2*i], 2) + pow(y[1] - y[2*i + 1], 2))*(12.0*pow(param[2], 2) + 6.0*param[2] - 6.0));
    }
}

void state_eval_wr(int l_edge, const real_t * restrict in, real_t * restrict out) {
    state_eval(l_edge,in, in+2*(l_edge-0), in+2*(l_edge-0)+2*(l_edge-0), in+2*(l_edge-0)+2*(l_edge-0)+2*(l_edge-0),
               out);
}

kernel_t kern_state = {
    .nmap = 2,
    .maps = {
        { .dim=4, .v_start=-1, .v_end=-1 },
        { .dim=2, .v_start=0, .v_end=-1 }
    },

    .ninp = 4,
    .inp = {
        { .field_number=0, .map_num=1, .name="y" },
        { .field_number=1, .map_num=1, .name="v" },
        { .field_number=2, .map_num=1, .name="x" },
        { .field_number=3, .map_num=0, .name="param" }
    },

    .noutp = 1,
    .outp = {
        { .rank = 1, .nmap = 1, .map_nums = { 1 } }
    },
    .eval=state_eval_wr,
    .name="state"
};

