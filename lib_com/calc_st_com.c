/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"



/*----------------------------------------------------------------------------
 * calc_rc0_h()
 *
 * computes 1st parcor from composed filter impulse response
 *---------------------------------------------------------------------------*/

static void calc_rc0_h(
    const float *h,       /* i  : impulse response of composed filter */
    float *rc0      /* o  : 1st parcor */
)
{
    float acf0, acf1;
    float temp, temp2;
    const float *ptrs;
    int i;

    /* computation of the autocorrelation function acf */
    temp = (float) 0.;
    for (i = 0; i < LONG_H_ST; i++)
    {
        temp += h[i] * h[i];
    }
    acf0 = temp;

    temp = (float) 0.;
    ptrs = h;
    for (i = 0; i < LONG_H_ST - 1; i++)
    {
        temp2 = *ptrs++;
        temp += temp2 * (*ptrs);
    }
    acf1 = temp;

    /* Initialisation of the calculation */
    if (acf0 == (float) 0.)
    {
        *rc0 = (float) 0.;
        return;
    }

    /* Compute 1st parcor */
    if (acf0 < (float) fabs (acf1))
    {
        *rc0 = (float) 0.0;
        return;
    }
    *rc0 = -acf1 / acf0;

    return;
}


/*----------------------------------------------------------------------------
 * calc_st_filt()
 *
 * computes impulse response of A(gamma2) / A(gamma1)
 * controls gain : computation of energy impulse response as
 *                 SUMn  (abs (h[n])) and computes parcor0
 *---------------------------------------------------------------------------- */

void calc_st_filt(
    const float *apond2,      /* i  : coefficients of numerator               */
    const float *apond1,      /* i  : coefficients of denominator             */
    float *parcor0,     /* o  : 1st parcor calcul. on composed filter   */
    float *sig_ltp_ptr, /* i/o: input of 1/A(gamma1) : scaled by 1/g0   */
    float *mem_zero,    /* i/o: All zero memory                         */
    const short L_subfr,      /* i  : the length of subframe                  */
    const short extl          /* i  : extension layer info                    */

)
{
    float h[LONG_H_ST];
    float g0, temp;
    int i;

    /* compute i.r. of composed filter apond2 / apond1 */
    if( extl == SWB_TBE )
    {
        syn_filt( apond1, LPC_SHB_ORDER, apond2, h, LONG_H_ST, mem_zero, 0 );
    }
    else
    {
        syn_filt( apond1, M, apond2, h, LONG_H_ST, mem_zero, 0 );
    }

    /* compute 1st parcor */
    calc_rc0_h( h, parcor0 );

    /* compute g0 */
    g0 = (float) 0.;
    for (i = 0; i < LONG_H_ST; i++)
    {
        g0 += (float) fabs (h[i]);
    }

    /* Scale signal input of 1/A(gamma1) */
    if (g0 > (float) 1.)
    {
        temp = (float) 1. / g0;

        for (i = 0; i < L_subfr; i++)
        {
            sig_ltp_ptr[i] = sig_ltp_ptr[i] * temp;
        }
    }

    return;
}

/*----------------------------------------------------------------------------
 * filt_mu()
 *
 * tilt filtering with : (1 + mu z-1) * (1/1-|mu|)
 *      computes y[n] = (1/1-|mu|) (x[n]+mu*x[n-1])
 *---------------------------------------------------------------------------*/

void filt_mu(
    const float *sig_in,      /* i  : signal (beginning at sample -1) */
    float *sig_out,     /* o  : output signal                   */
    const float parcor0,      /* i  : parcor0 (mu = parcor0 * gamma3) */
    const short L_subfr,      /* i  : the length of subframe          */
    const short extl          /* i  : extension layer info            */
)
{
    short n;
    float mu, ga, temp;
    const float *ptrs;

    if( extl == SWB_TBE )
    {
        if(parcor0 > 0.0f)
        {
            mu = parcor0 * GAMMA3_PLUS_WB;
        }
        else
        {
            mu = parcor0 * GAMMA3_MINUS_WB;
        }
    }
    else
    {
        if (parcor0 > 0.0f)
        {
            mu = parcor0 * GAMMA3_PLUS;
        }
        else
        {
            mu = parcor0 * GAMMA3_MINUS;
        }
    }

    ga = (float) 1. / ((float) 1. - (float) fabs (mu));

    ptrs = sig_in;                /* points on sig_in(-1) */

    for (n = 0; n < L_subfr; n++)
    {
        temp = mu * (*ptrs++);
        temp += (*ptrs);
        sig_out[n] = ga * temp;
    }

    return;
}




/*----------------------------------------------------------------------------
 * scale_st()
 *
 * control of the subframe gain
 * gain[n] = AGC_FAC_FX * gain[n-1] + (1 - AGC_FAC_FX) g_in/g_out
 *---------------------------------------------------------------------------*/

void scale_st(
    const float *sig_in,      /* i  : postfilter input signal         */
    float *sig_out,     /* i/o: postfilter output signal        */
    float *gain_prec,   /* i/o: last value of gain for subframe */
    const short L_subfr,      /* i  : the length of subframe          */
    const short extl          /* i  : extension layer info            */
)
{
    int i;
    float gain_in, gain_out;
    float g0, gain;
    float agc_fac1_para = 0.0f;
    float agc_fac_para = 0.0f;

    if( extl == SWB_TBE )
    {
        agc_fac1_para = AGC_FAC1_WB;
        agc_fac_para = AGC_FAC_WB;
    }
    else
    {
        agc_fac1_para = AGC_FAC1;
        agc_fac_para = AGC_FAC;
    }

    /* compute input gain */
    gain_in = (float) 0.;
    for (i = 0; i < L_subfr; i++)
    {
        gain_in += (float) fabs (sig_in[i]);
    }

    if ( gain_in == 0.0f )
    {
        g0 = 0.0f;
    }
    else
    {
        /* Compute output gain */
        gain_out = 0.0f;
        for (i = 0; i < L_subfr; i++)
        {
            gain_out += (float) fabs (sig_out[i]);
        }

        if (gain_out == 0.0f)
        {
            *gain_prec = 0.0f;
            return;
        }

        g0 = gain_in / gain_out;
        g0 *= agc_fac1_para;
    }

    /* compute gain(n) = AGC_FAC gain(n-1) + (1-AGC_FAC)gain_in/gain_out */
    /* sig_out(n) = gain(n) sig_out(n) */
    gain = *gain_prec;
    for (i = 0; i < L_subfr; i++)
    {
        gain *= agc_fac_para;
        gain += g0;
        sig_out[i] *= gain;
    }

    *gain_prec = gain;

    return;
}

void blend_subfr2( float *sigIn1, float *sigIn2, float *sigOut)
{

    float fac1 = 1.f - (1.f / L_SUBFR);
    float fac2 = 0.f + (1.f / L_SUBFR);
    float step = 1.f / (L_SUBFR/2);
    int i;

    for(i=0; i<L_SUBFR/2; i++)
    {
        sigOut[i] = fac1 * sigIn1[i] + fac2 * sigIn2[i];
        fac1 -= step;
        fac2 += step;
    }

    return;
}
