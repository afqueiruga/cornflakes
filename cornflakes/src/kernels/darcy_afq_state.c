#include "darcy_afq_state.h"

#include "math.h"

void darcy_afq_eval(int l_edge,
                    /*Inputs:*/
                    const real_t * restrict p,
                    const real_t * restrict x,
                    const real_t * restrict param,
                    /*Outputs:*/
                    real_t * restrict FQ,
                    real_t * restrict KQ)
{
    /* Evaluation of FQ */
    FQ[0] = 2*param[1]*param[6]*(1.0 - sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2))/param[0])*(-p[0] + p[1])/(param[0]*param[3]*sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)));
    FQ[1] = -2*param[1]*param[2]*(1.0 - sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2))/param[0])*(-p[0] + p[1])/(param[0]*param[7]*sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)));

    /* Evaluation of KQ */
    KQ[0] = -2*param[1]*param[6]*(1.0 - sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2))/param[0])/(param[0]*param[3]*sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)));
    KQ[1] = 2*param[1]*param[6]*(1.0 - sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2))/param[0])/(param[0]*param[3]*sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)));
    KQ[2] = 2*param[1]*param[2]*(1.0 - sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2))/param[0])/(param[0]*param[7]*sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)));
    KQ[3] = -2*param[1]*param[2]*(1.0 - sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2))/param[0])/(param[0]*param[7]*sqrt(pow(-x[0] + x[2], 2) + pow(-x[1] + x[3], 2)));
}

void darcy_afq_eval_wr(int l_edge, const real_t * restrict in, real_t * restrict out) {
    darcy_afq_eval(l_edge,in, in+2, in+2+4,
                   out, out+2);
}

kernel_t kernel_darcy_afq = {
    .nmap = 3,
    .maps = {
        { .dim=1, .v_start=0, .v_end=2 },
        { .dim=2, .v_start=0, .v_end=2 },
        { .dim=4, .v_start=0, .v_end=2 }
    },

    .ninp = 3,
    .inp = {
        { .field_number=0, .map_num=0, .name="p" },
        { .field_number=1, .map_num=1, .name="x" },
        { .field_number=2, .map_num=2, .name="param" }
    },

    .noutp = 2,
    .outp = {
        { .rank = 1, .nmap = 1, .map_nums = { 0 } },
        { .rank = 2, .nmap = 1, .map_nums = { 0 } }
    },
    .eval=darcy_afq_eval_wr,
    .name="darcy_afq"
};

