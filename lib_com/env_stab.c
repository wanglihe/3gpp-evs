/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"
#include <math.h>

/*--------------------------------------------------------------------------*
 * Local constants
 *--------------------------------------------------------------------------*/

#define ENV_STAB_SMO_HO 10                      /* number of hangover frames when switching from music to speech state */

/*--------------------------------------------------------------------------*/
/*  Function  env_stability()                                               */
/*  ~~~~~~~~~~~~~~~~~~~~~~~~~                                               */
/*                                                                          */
/*  Envelope stability measure                                              */
/*--------------------------------------------------------------------------*/

float env_stability(
    const short *ynrm,         /*i  : Norm vector for current frame          */
    const short nb_sfm,        /*i  : Number of sub-bands                    */
    short *mem_norm,     /*i/o: Norm vector memory from past frame     */
    float *mem_env_delta /*i/o: Envelope stability memory for smoothing*/
)
{
    float env_delta;
    float env_stab;
    float tmp,tmp_stab;
    short i;

    /* Calculate envelope stability parameter */
    env_delta = 0;
    for (i = 0; i < nb_sfm; i++)
    {
        tmp = (float)(mem_norm[i] - ynrm[i]);
        env_delta += tmp * tmp;
        mem_norm[i] = ynrm[i];
    }
    env_delta = (float) sqrt(env_delta / nb_sfm);

    env_delta = ENV_SMOOTH_FAC * env_delta + (1-ENV_SMOOTH_FAC) **mem_env_delta;
    *mem_env_delta = env_delta;

    tmp_stab = (float)fabs(env_delta - M_STAB_TBL);

    /* Table lookup for smooth transitions */
    i = (short)usquant(tmp_stab, &tmp, 0, D_STAB_TBL, L_STAB_TBL);
    if(env_delta < M_STAB_TBL)
    {
        env_stab = 1.0f - stab_trans[i];
    }
    else
    {
        env_stab = stab_trans[i];
    }

    return env_stab;
}

/*--------------------------------------------------------------------------*
 * env_stab_smo_fx()
 *
 *
 *--------------------------------------------------------------------------*/

float env_stab_smo(
    float env_stab,            /*i  : env_stab value                         */
    float *env_stab_state_p,   /*i/o: env_stab state probabilities           */
    short *ho_cnt              /*i/o: hangover counter for speech state      */
)
{
    short state, prev_state;
    float maxval, pp[NUM_ENV_STAB_PLC_STATES], pa[NUM_ENV_STAB_PLC_STATES];
    /* get previous state */
    prev_state = maximum(env_stab_state_p,NUM_ENV_STAB_PLC_STATES,&maxval);

    /* assume two states: speech(0), music(1) */
    /* set a posteriori likelihoods for the two states according to env_stab */
    env_stab = (env_stab - stab_trans[L_STAB_TBL-1])/(1-2*stab_trans[L_STAB_TBL-1]);
    pp[0] = 1.0f-env_stab;
    pp[1] = env_stab;

    /* calculate a priori likelihoods */
    pa[0] = dotp(env_stab_tp[0],env_stab_state_p,NUM_ENV_STAB_PLC_STATES);
    pa[1] = dotp(env_stab_tp[1],env_stab_state_p,NUM_ENV_STAB_PLC_STATES);

    /* multiply elementwise with a posteriori likelihoods */
    v_mult(pa,pp,env_stab_state_p,NUM_ENV_STAB_PLC_STATES);

    /* renormalize state probabilities */
    v_multc(env_stab_state_p,1.0f/sum_f(env_stab_state_p,NUM_ENV_STAB_PLC_STATES),env_stab_state_p,NUM_ENV_STAB_PLC_STATES);

    /* find maximum index as return value */
    state = maximum(env_stab_state_p,NUM_ENV_STAB_PLC_STATES,&maxval);

    /* apply some hangover for speech */
    if (state==0 && prev_state==1)
    {
        *ho_cnt=ENV_STAB_SMO_HO;
    }
    if (*ho_cnt>0)
    {
        pp[0]=1;
        pp[1]=0;
        (*ho_cnt)--;
    }

    return state;
}
