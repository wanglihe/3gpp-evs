/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"

/*---------------------------------------------------------------------*
 * Local functions
 *---------------------------------------------------------------------*/

static void hp400_12k8(float signal[], const short lg,  float mem[]);
static void filt_6k_8k(float signal[], const short lg, float mem[]);

/*---------------------------------------------------------------------*
 * hf_cod_init()
 *
 *
 *---------------------------------------------------------------------*/

void hf_cod_init(
    float *mem_hp400_enc,         /* o: memory of hp 400 Hz filter   */
    float *mem_hf1_enc,           /* o: HF band-pass filter memory   */
    float *mem_syn_hf_enc,        /* o: HF synthesis memory          */
    float *mem_hf2_enc,           /* o: HF band-pass filter memory   */
    float *gain_alpha             /* o: smoothing gain for transitions between active and inactive frames */
)
{
    set_f( mem_hp400_enc, 0, 4 );
    set_f( mem_hf1_enc, 0, L_FIR-1 );
    set_f( mem_syn_hf_enc, 0, M );
    set_f( mem_hf2_enc, 0, L_FIR-1 );
    *gain_alpha = 1.0;

    return;
}


/*---------------------------------------------------------------------*
 * hf_cod()
 *
 *
 *---------------------------------------------------------------------*/

void hf_cod(
    const long  core_brate,                     /* i  : core bitrate                 */
    const float *speech16k,                     /* i  : original speech at 16 kHz    */
    const float Aq[],                           /* i  : quantized Aq                 */
    const float exc[],                          /* i  : excitation at 12.8 kHz       */
    float synth[],                        /* i  : 12.8kHz synthesis signal     */
    short *seed2_enc,                     /* i/o: random seed for HF noise gen */
    float *mem_hp400_enc,                 /* i/o: memory of hp 400 Hz filter   */
    float *mem_syn_hf_enc,                /* i/o: HF synthesis memory          */
    float *mem_hf1_enc,                   /* i/o: HF band-pass filter memory   */
    float *mem_hf2_enc,                   /* i/o: HF band-pass filter memory   */
    const short *dtxHangoverCount,
    float *gain_alpha,                    /* i/o: smoothing gain for transitions between active and inactive frames */
    short *hf_gain                        /* o  :  HF gain to be transmitted to decoder */
)
{
    short i;
    float ener_hf, ener_exc, ener_input, fac, HF_syn[L_SUBFR16k], tmp, ener, scale;
    float Ap[M16k+1];
    float HF_SP[L_SUBFR16k];
    float HF_est_gain;
    float HF_calc_gain;
    float HF_corr_gain;
    short HF_gain_ind;
    float dist_min, dist;
    float HF[L_SUBFR16k];                        /* o  : HF excitation                */

    /* Original speech signal as reference for high band gain quantisation */
    for (i = 0; i < L_SUBFR16k; i++)
    {
        HF_SP[i] = speech16k[i];
    }

    /*-----------------------------------------------------------------*
     * generate white noise vector
     *-----------------------------------------------------------------*/

    for (i=0; i<L_SUBFR16k; i++)
    {
        HF[i] = (float)own_random(seed2_enc);
    }

    /*-----------------------------------------------------------------*
     * calculate energy scaling factor so that white noise would have the
     * same energy as exc12k8
     *-----------------------------------------------------------------*/

    ener_exc = 0.01f;
    for (i=0; i<L_SUBFR; i++)
    {
        ener_exc += exc[i]*exc[i];
    }

    ener_hf = 0.01f;
    for (i=0; i<L_SUBFR16k; i++)
    {
        ener_hf += HF[i]*HF[i];
    }

    scale = (float)(sqrt(ener_exc/ener_hf));

    for (i=0; i<L_SUBFR16k; i++)
    {
        HF[i] *= scale;
    }

    /*-----------------------------------------------------------------*
     * calculate energy scaling factor to respect tilt of synth12k8
     * (tilt: 1=voiced, -1=unvoiced)
     *-----------------------------------------------------------------*/

    hp400_12k8( synth, L_SUBFR, mem_hp400_enc );

    ener = 0.001f;
    tmp = 0.001f;
    for (i=1; i<L_SUBFR; i++)
    {
        ener += synth[i]*synth[i];  /* ener = r[0] */
        tmp += synth[i]*synth[i-1]; /* tmp  = r[1] */
    }
    fac = tmp/ener;

    HF_est_gain = (float)(1.0f - fac);
    if( core_brate == SID_1k75 || core_brate == FRAME_NO_DATA )
    {
        /* emphasize HF noise in CNG */
        HF_est_gain *= 1.25f; /* full alignment with G.722 and AMR-WB */
    }

    if (HF_est_gain < 0.1f)
    {
        HF_est_gain = 0.1f;
    }

    if (HF_est_gain > 1.0f) /* this condition is not in G.722.2, but in AMR-WB!*/
    {
        HF_est_gain = 1.0f;
    }

    /*-----------------------------------------------------------------*
     * synthesis of noise: 4.8kHz..5.6kHz --> 6kHz..7kHz
     *-----------------------------------------------------------------*/

    weight_a( Aq, Ap, 0.6f, M );
    syn_filt( Ap, M, HF, HF_syn, L_SUBFR16k, mem_syn_hf_enc, 1 );

    /*-----------------------------------------------------------------*
     * high pass filtering (0.9375ms of delay = 15 samples@16k)
     *-----------------------------------------------------------------*/

    filt_6k_8k( HF_syn, L_SUBFR16k, mem_hf1_enc );

    /* filtering of the original signal */
    filt_6k_8k( HF_SP, L_SUBFR16k, mem_hf2_enc );

    /* check the gain difference */
    ener_input = 0.01f;
    ener_hf = 0.01f;
    for ( i=0; i<L_SUBFR16k; i++ )
    {
        ener_input += HF_SP[i]*HF_SP[i];
        ener_hf += HF_syn[i]*HF_syn[i];
    }

    HF_calc_gain = (float)sqrt(ener_input/ener_hf);

    /* set energy of HF synthesis to energy of original HF:
       cross-fade between HF levels in active and inactive frame in hangover period */

    *gain_alpha *= (float)(10-(*dtxHangoverCount))/7.0f;
    if ((10-(*dtxHangoverCount)) > 6)
    {
        *gain_alpha = 1.0f;
    }

    HF_corr_gain = (*gain_alpha)*HF_calc_gain + (1.0f-(*gain_alpha))*HF_est_gain;
    HF_corr_gain /= 2.0f;  /* to stay in alignlent with AMR-WB legacy decoder where decoded gain is multiplied by 2 */

    /* Quantize the correction gain */
    dist_min = 100000.0f;
    HF_gain_ind = 0;
    for ( i = 0; i < 16; i++ )
    {
        dist = (HF_corr_gain-HP_gain[i])*(HF_corr_gain-HP_gain[i]);
        if (dist_min > dist)
        {
            dist_min = dist;
            HF_gain_ind = i;
        }
    }

    *hf_gain = HF_gain_ind;

    return;
}

/*-----------------------------------------------------------------------*
 * Function hp400_12k8()                                                 *
 *                                                                       *
 * 2nd order Cheb2 high pass filter with cut off frequency at 400 Hz.    *
 * Optimized for fixed-point to get the following frequency response  :  *
 *                                                                       *
 *  frequency  :   0Hz   100Hz  200Hz  300Hz  400Hz  630Hz  1.5kHz  3kHz *
 *  dB loss  :   -infdB  -30dB  -20dB  -10dB  -3dB   +6dB    +1dB    0dB *
 *                                                                       *
 * Algorithm  :                                                          *
 *                                                                       *
 *  y[i] = b[0]*x[i] + b[1]*x[i-1] + b[2]*x[i-2]                         *
 *                   + a[1]*y[i-1] + a[2]*y[i-2];                        *
 *                                                                       *
 *  short b[3] = {3660, -7320,  3660};       in Q12                      *
 *  short a[3] = {4096,  7320, -3540};       in Q12                      *
 *                                                                       *
 *  float -->   b[3] = {0.893554687, -1.787109375,  0.893554687};        *
 *              a[3] = {1.000000000,  1.787109375, -0.864257812};        *
 *-----------------------------------------------------------------------*/

static void hp400_12k8(
    float signal[],  /* i/o: signal            */
    const short lg,        /* i  : lenght of signal  */
    float mem[]      /* i/o: filter memory [4] */
)
{
    short i;
    float x0, x1, x2;
    float yy0, yy1, y2;

    yy1 = mem[0];
    y2 = mem[1];
    x0 = mem[2];
    x1 = mem[3];
    for(i=0; i<lg; i++)
    {
        x2 = x1;
        x1 = x0;
        x0 = signal[i];
        yy0 = yy1*a_hp400[1] + y2*a_hp400[2] + x0*b_hp400[0] + x1*b_hp400[1] + x2*b_hp400[2];

        signal[i] = yy0;
        y2 = yy1;
        yy1 = yy0;
    }

    mem[0] = yy1;
    mem[1] = y2;
    mem[2] = x0;
    mem[3] = x1;

    return;
}

/*-------------------------------------------------------------------*
 * filt_6k_7k:
 *
 * 15th order band pass 6kHz to 7kHz FIR filter.
 *
 * frequency  :4kHz   5kHz  5.5kHz  6kHz  6.5kHz 7kHz  7.5kHz  8kHz
 * dB loss  : -60dB  -45dB  -13dB   -3dB   0dB   -3dB  -13dB  -45dB
 * (gain=4.0)
 *-------------------------------------------------------------------*/

static void filt_6k_8k(
    float signal[],  /* i/o: signal        */
    const short lg,        /* i  : signal length */
    float mem[]      /* i/o: filter memory */
)
{
    short i, j;
    float s, x[L_FRAME48k/NB_SUBFR+(L_FIR-1)];

    for( i=0; i<(L_FIR-1); i++ )
    {
        x[i] = mem[i];
    }

    for(i=0; i<lg; i++)
    {
        x[i+(L_FIR-1)] = signal[i];
    }

    for(i=0; i<lg; i++)
    {
        s = 0.0;
        for( j=0; j<L_FIR; j++ )
        {
            s += x[i+j] * fir_6k_8k[j];
        }

        signal[i] = (float)(s*1.0f);
    }

    for(i=0; i<(L_FIR-1); i++)
    {
        mem[i] = x[i+lg];
    }

    return;
}


