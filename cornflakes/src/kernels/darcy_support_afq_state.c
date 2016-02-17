#include "darcy_support_afq_state.h"

#include "math.h"

void darcy_support_afq_eval(int l_edge,
                            /*Inputs:*/
                            const real_t * restrict p,
                            const real_t * restrict x,
                            const real_t * restrict param,
                            /*Outputs:*/
                            real_t * restrict SUP)
{

    if( pow(1.0 - sqrt(pow(-0.288675134594813*sqrt(param[6]) - x[0] + x[2], 2) + pow(-0.288675134594813*sqrt(param[6]) - x[1] + x[3], 2))/param[0], 2) > 0.0 ) {
        /* Evaluation of SUP */
        SUP[0] += 0.25*param[6]*pow(1.0 - sqrt(pow(-0.288675134594813*sqrt(param[6]) - x[0] + x[2], 2) + pow(-0.288675134594813*sqrt(param[6]) - x[1] + x[3], 2))/param[0], 2);
        SUP[1] += 0.25*param[6]*pow(1.0 - sqrt(pow(-0.288675134594813*sqrt(param[6]) - x[0] + x[2], 2) + pow(-0.288675134594813*sqrt(param[6]) - x[1] + x[3], 2))/param[0], 2);
    } else {

    }



    if( pow(1.0 - sqrt(pow(-0.288675134594813*sqrt(param[6]) - x[0] + x[2], 2) + pow(0.288675134594813*sqrt(param[6]) - x[1] + x[3], 2))/param[0], 2) > 0.0 ) {
        /* Evaluation of SUP */
        SUP[0] += 0.25*param[6]*pow(1.0 - sqrt(pow(-0.288675134594813*sqrt(param[6]) - x[0] + x[2], 2) + pow(0.288675134594813*sqrt(param[6]) - x[1] + x[3], 2))/param[0], 2);
        SUP[1] += 0.25*param[6]*pow(1.0 - sqrt(pow(-0.288675134594813*sqrt(param[6]) - x[0] + x[2], 2) + pow(0.288675134594813*sqrt(param[6]) - x[1] + x[3], 2))/param[0], 2);
    } else {

    }



    if( pow(1.0 - sqrt(pow(-0.288675134594813*sqrt(param[6]) - x[1] + x[3], 2) + pow(0.288675134594813*sqrt(param[6]) - x[0] + x[2], 2))/param[0], 2) > 0.0 ) {
        /* Evaluation of SUP */
        SUP[0] += 0.25*param[6]*pow(1.0 - sqrt(pow(-0.288675134594813*sqrt(param[6]) - x[1] + x[3], 2) + pow(0.288675134594813*sqrt(param[6]) - x[0] + x[2], 2))/param[0], 2);
        SUP[1] += 0.25*param[6]*pow(1.0 - sqrt(pow(-0.288675134594813*sqrt(param[6]) - x[1] + x[3], 2) + pow(0.288675134594813*sqrt(param[6]) - x[0] + x[2], 2))/param[0], 2);
    } else {

    }



    if( pow(1.0 - sqrt(pow(0.288675134594813*sqrt(param[6]) - x[0] + x[2], 2) + pow(0.288675134594813*sqrt(param[6]) - x[1] + x[3], 2))/param[0], 2) > 0.0 ) {
        /* Evaluation of SUP */
        SUP[0] += 0.25*param[6]*pow(1.0 - sqrt(pow(0.288675134594813*sqrt(param[6]) - x[0] + x[2], 2) + pow(0.288675134594813*sqrt(param[6]) - x[1] + x[3], 2))/param[0], 2);
        SUP[1] += 0.25*param[6]*pow(1.0 - sqrt(pow(0.288675134594813*sqrt(param[6]) - x[0] + x[2], 2) + pow(0.288675134594813*sqrt(param[6]) - x[1] + x[3], 2))/param[0], 2);
    } else {

    }

}

void darcy_support_afq_eval_wr(int l_edge, const real_t * restrict in, real_t * restrict out) {
    darcy_support_afq_eval(l_edge,in, in+2, in+2+4,
                           out);
}

kernel_t kernel_darcy_support_afq = {
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

    .noutp = 1,
    .outp = {
        { .rank = 1, .nmap = 1, .map_nums = { 0 } }
    },
    .eval=darcy_support_afq_eval_wr,
    .name="darcy_support_afq"
};

