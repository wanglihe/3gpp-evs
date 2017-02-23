/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <stdlib.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"

/*-----------------------------------------------------------------*
 * Local functions
 *-----------------------------------------------------------------*/

static void create_random_vector( float output[], const short length, short seed[] );
static void flip_spectrum( const float input[], float output[], const short length );
static void Hilbert_transform( float tmp_R[], float tmp_I[], float *tmpi_R, float *tmpi_I, const short length, const short HB_stage_id );

/*-------------------------------------------------------------------*
 * swb_tbe_reset()
 *
 * Reset the SWB TBE encoder
 *-------------------------------------------------------------------*/

void swb_tbe_reset(
    float mem_csfilt[],
    float mem_genSHBexc_filt_down_shb[],
    float state_lpc_syn[],
    float syn_overlap[],
    float state_syn_shbexc[],
    float *tbe_demph,
    float *tbe_premph,
    float mem_stp_swb[],
    float *gain_prec_swb
)
{
    set_f( mem_csfilt, 0, 2);
    set_f( mem_genSHBexc_filt_down_shb, 0.0f, (2*ALLPASSSECTIONS_STEEP+1) );
    set_f( state_lpc_syn, 0.0f, LPC_SHB_ORDER );
    set_f( syn_overlap, 0.0f, L_SHB_LAHEAD );
    set_f( state_syn_shbexc, 0.0f, L_SHB_LAHEAD );
    *tbe_demph = 0.0f;
    *tbe_premph = 0.0f;
    set_f(mem_stp_swb, 0, LPC_SHB_ORDER);
    *gain_prec_swb = 1.0f;

    return;
}

/*-------------------------------------------------------------------*
 * swb_tbe_reset_synth()
 *
 * Reset the extra parameters needed for synthesis of the SWB TBE output
 *-------------------------------------------------------------------*/

void swb_tbe_reset_synth(
    float genSHBsynth_Hilbert_Mem[],
    float genSHBsynth_state_lsyn_filt_shb_local[]
)
{
    set_f( genSHBsynth_Hilbert_Mem, 0.0f, HILBERT_MEM_SIZE );
    set_f( genSHBsynth_state_lsyn_filt_shb_local, 0.0f, 2 * ALLPASSSECTIONS_STEEP );

    return;
}

/*-------------------------------------------------------------------*
 * tbe_celp_exc_offset()
 *
 * Compute tbe bwe celp excitation offset
 *-------------------------------------------------------------------*/
short tbe_celp_exc_offset(
    const short T0,               /* i  : Integer pitch */
    const short T0_frac           /* i  : Fractional part of the pitch */
)
{
    short offset;
    offset = T0 * HIBND_ACB_L_FAC
             + (int) ((float) T0_frac * 0.25f * HIBND_ACB_L_FAC + 2 * HIBND_ACB_L_FAC + 0.5f)
             - 2 * HIBND_ACB_L_FAC;

    return offset;
}
/*-------------------------------------------------------------------*
 * flip_and_downmix_generic()
 *
 * flips the spectrum and downmixes the signals, lpf if needed
 *-------------------------------------------------------------------*/

void flip_and_downmix_generic(
    float input[],            /* i  : input spectrum          */
    float output[],           /* o  : output  spectrum        */
    const short length,             /* i  : length of spectra       */
    float mem1_ext[],         /* i/o: Hilbert filter memory   */
    float mem2_ext[],         /* i/o: memory                  */
    float mem3_ext[],         /* i/o: memory                  */
    short *phase_state        /* i/o: Phase state in case frequency isn't multiple of 50 Hz */
)
{
    short i, j;
    float tmp[L_FRAME32k + HILBERT_ORDER1];
    float tmpi_R[L_FRAME32k];
    float tmpi_I[L_FRAME32k];
    float tmpi2_R[L_FRAME32k + HILBERT_ORDER2];
    float tmpi2_I[L_FRAME32k + HILBERT_ORDER2];
    float tmp_R[L_FRAME32k + HILBERT_ORDER2];
    float tmp_I[L_FRAME32k + HILBERT_ORDER2];
    short k, period;
    float recip_period;
    float local_negsin_table[L_FRAME16k];
    float local_cos_table[L_FRAME16k];

    period = 17;    /* == (short) (32000.0f / 1850.0f + 0.5f); */

    recip_period = 256.0f / (float) period;
    for( i=0; i<period; i++ )
    {
        k = (short) (i * recip_period + 0.5f);
        if( k <= 64 )
        {
            local_negsin_table[i] = -sincos_t[k];
            local_cos_table[i] = sincos_t[64 - k];
        }
        else if( k <= 128 )
        {
            local_negsin_table[i] = -sincos_t[128 - k];
            local_cos_table[i] = -sincos_t[k - 64];
        }
        else if( k <= 192 )
        {
            local_negsin_table[i] = sincos_t[k - 128];
            local_cos_table[i] = -sincos_t[192 - k];
        }
        else
        {
            local_negsin_table[i] = sincos_t[256 - k];
            local_cos_table[i] = sincos_t[k - 192];
        }
    }

    for( i = 0; i < length; i = i+2 )
    {
        input[i] = -input[i];
    }

    mvr2r( input, tmp + HILBERT_ORDER1, length );

    mvr2r( mem1_ext, tmp, HILBERT_ORDER1 );

    /* Hilber transform stage - 0 */
    Hilbert_transform( tmp,                /* i: Real component of HB */
                       tmp,                /* i: Imag component of HB */
                       tmpi_R,             /* o: Real component of HB */
                       tmpi_I,             /* o: Imag. component of HB */
                       length,             /* i: length of the spectra */
                       0);                 /* i: HB transform stage */

    mvr2r( mem2_ext, tmpi2_R, HILBERT_ORDER2 );
    mvr2r( mem3_ext, tmpi2_I, HILBERT_ORDER2 );

    /* Hilber transform stage - 1 */
    Hilbert_transform( tmpi_R,             /* i: Real component of HB */
                       tmpi_I,             /* i: Imag component of HB */
                       tmpi2_R,            /* o: Real component of HB */
                       tmpi2_I,            /* o: Imag. component of HB */
                       length,             /* i: length of the spectra */
                       1);                 /* i: HB transform stage */

    mvr2r( tmp + length, mem1_ext, HILBERT_ORDER1 );
    mvr2r( mem2_ext+HILBERT_ORDER2, tmp_R, HILBERT_ORDER2 );
    mvr2r( mem3_ext+HILBERT_ORDER2, tmp_I, HILBERT_ORDER2 );

    /* Hilber transform stage - 2 */
    Hilbert_transform( tmpi2_R,            /* i: Real component of HB */
                       tmpi2_I,            /* i: Imag component of HB */
                       tmpi_R,             /* o: Real component of HB */
                       tmpi_I,             /* o: Imag. component of HB */
                       length,             /* i: length of the spectra */
                       2);                 /* i: HB transform stage */

    mvr2r( tmpi2_R + length, mem2_ext, HILBERT_ORDER2 );
    mvr2r( tmpi2_I + length, mem3_ext, HILBERT_ORDER2 );

    /* Hilber transform stage - 3 */
    Hilbert_transform( tmpi_R,             /* i: Real component of HB */
                       tmpi_I,             /* i: Imag component of HB */
                       tmp_R,              /* o: Real component of HB */
                       tmp_I,              /* o: Imag. component of HB */
                       length,             /* i: length of the spectra */
                       3);                 /* i: HB transform stage */

    mvr2r( tmp_R + length, mem2_ext+HILBERT_ORDER2, HILBERT_ORDER2 );
    mvr2r( tmp_I + length, mem3_ext+HILBERT_ORDER2, HILBERT_ORDER2 );

    if( *phase_state >= period )
    {
        *phase_state = 0;
    }

    for( i=0, j=*phase_state; i < length; )
    {
        for( ; (j < period) && (i < length); j++, i++ )
        {
            output[i] = tmp_R[i + HILBERT_ORDER2] * local_cos_table[j] + tmp_I[i + HILBERT_ORDER2] * local_negsin_table[j];
        }

        if( j >= period )
        {
            j = 0;
        }
    }

    *phase_state = j;

    return;
}

/*----------------------------------------------
 * Hilbert_transform()
 *
 * Hilbert transform
 *------------------------------------------------*/

static void Hilbert_transform(
    float tmp_R[],             /* i: Real component of HB */
    float tmp_I[],             /* i: Real component of HB */
    float tmpi_R[],            /* o: Real component of HB */
    float tmpi_I[],            /* o: Imag. component of HB */
    const short length,              /* i: input length */
    const short HB_stage_id          /* i: HB transform stage */
)
{
    short i, hb_filter_stage, offset;

    hb_filter_stage = 2*HB_stage_id;
    offset = (HB_stage_id == 0) ? 1 : 0;

    if (HB_stage_id == 0 || HB_stage_id == 2)
    {
        for( i=0; i<length; i++ )
        {
            tmpi_R[i] =   tmp_R[i + 4] * Hilbert_coeffs[hb_filter_stage][0 + offset]
                          + tmp_R[i + 2] * Hilbert_coeffs[hb_filter_stage][2 + offset]
                          + tmp_R[i]     * Hilbert_coeffs[hb_filter_stage][4 + offset];

            tmpi_I[i] =   tmp_I[i + 4 + offset] * Hilbert_coeffs[hb_filter_stage+1][0]
                          + tmp_I[i + 2 + offset] * Hilbert_coeffs[hb_filter_stage+1][2]
                          + tmp_I[i + offset]     * Hilbert_coeffs[hb_filter_stage+1][4];
        }
    }
    else if (HB_stage_id == 1 || HB_stage_id == 3)
    {
        for( i=0; i<length; i++ )
        {
            tmpi_R[i+4] =    tmp_R[i]
                             - tmpi_R[i + 2] * Hilbert_coeffs[hb_filter_stage][2]
                             - tmpi_R[i]     * Hilbert_coeffs[hb_filter_stage][4];

            tmpi_I[i+4] =    tmp_I[i]
                             - tmpi_I[i + 2] * Hilbert_coeffs[hb_filter_stage+1][2]
                             - tmpi_I[i]     * Hilbert_coeffs[hb_filter_stage+1][4];
        }
    }

    return;
}

/*-------------------------------------------------------------------*
 * flip_spectrum()
 *
 *
 *-------------------------------------------------------------------*/

void flip_spectrum(
    const float input[],        /* i  : input spectrum          */
    float output[],       /* o  : output  spectrum        */
    const short length          /* i  : vector length           */
)
{
    short i;

    for( i = 0; i < length; i = i+2 )
    {
        output[i] = -input[i];
        output[i+1] = input[i+1];
    }

    return;
}

/*-------------------------------------------------------------------*
 * flip_spectrum_and_decimby4()
 *
 *
 *-------------------------------------------------------------------*/

void flip_spectrum_and_decimby4(
    const float input[],        /* i  : input spectrum          */
    float output[],       /* o  : output  spectrum        */
    const short length,         /* i  : vector length           */
    float mem1[],         /* i/o  : memory                */
    float mem2[],         /* i/o  : memory                */
    const short ramp_flag       /*i: flag to trigger slow ramp-up of output following change of core (HQ to ACELP or 12k8 to 16k ACELP) */
)
{
    short i;
    float factor, tmp[L_FRAME16k/2];
    float input_change[L_FRAME16k];

    if( ramp_flag )
    {
        factor = 4.0f / length;
        for( i = 0; i < length/4; i = i+2 )
        {
            input_change[i] = -input[i] * (i * factor);
            input_change[i+1] = input[i+1] * ((i + 1.0f) * factor);
        }
    }
    else
    {
        i = 0;
    }

    for( ; i < length; i = i+2 )
    {
        input_change[i] = -input[i];
        input_change[i+1] = input[i+1];
    }

    Decimate_allpass_steep( input_change, mem1, L_FRAME16k, tmp);
    Decimate_allpass_steep( tmp, mem2, L_FRAME16k/2, output);

    return;
}

/*-------------------------------------------------------------------*
 * GenShapedWBExcitation()
 *
 * Synthesize spectrally shaped highband excitation signal for the wideband
 *-------------------------------------------------------------------*/

void GenShapedWBExcitation(
    float *excSHB,                     /* o   : synthesized shaped shb exctiation     */
    const float *lpc_shb,                    /* i   : lpc coefficients                      */
    float *exc4kWhtnd,                 /* o   : whitened synthesized shb excitation   */
    float *mem_csfilt,                 /* i/o : memory                                */
    float *mem_genSHBexc_filt_down1,   /* i/o : memory                                */
    float *mem_genSHBexc_filt_down2,   /* i/o : memory                                */
    float *mem_genSHBexc_filt_down3,   /* i/o : memory                                */
    float *state_lpc_syn,              /* i/o : memory                                */
    const short coder_type,                  /* i   : coding type                           */
    const float *bwe_exc_extended,           /* i   : bandwidth extended exciatation        */
    short bwe_seed[],                  /* i/o : random number generator seed          */
    const float voice_factors[],             /* i   : voicing factor                        */
    const short uv_flag,                     /* i   : unvoiced flag                         */
    const short igf_flag
)
{
    short i, j, k;
    float wht_fil_mem [ LPC_WHTN_ORDER_WB ];
    float lpc_whtn[ LPC_WHTN_ORDER_WB + 1];
    float R[LPC_WHTN_ORDER_WB + 2];
    float excTmp[ L_FRAME16k];
    float excTmp2[ L_FRAME16k / 4];
    float exc4k[ L_FRAME16k / 4];
    float pow1, pow2, scale;
    float excNoisyEnv[ L_FRAME16k / 4];
    float csfilt_num2[1] = {0.05f};
    float csfilt_den2[2] = {1.0f, -0.96f};
    float temp1, temp2;
    float ervec[LPC_WHTN_ORDER_WB+2];
    float tmp_vfac;
    float avg_voice_fac = 0.25f*sum_f(voice_factors, NB_SUBFR);

    if( igf_flag && (coder_type == VOICED || avg_voice_fac > 0.35f) )
    {
        csfilt_num2[0] = 0.2f;
        csfilt_den2[1] = -0.8f;
    }
    else if( igf_flag && (coder_type == UNVOICED || avg_voice_fac < 0.2f) )
    {
        csfilt_num2[0] = 0.01f;
        csfilt_den2[1] = -0.99f;
    }
    set_f( wht_fil_mem, 0, LPC_WHTN_ORDER_WB );

    Decimate_allpass_steep( bwe_exc_extended, mem_genSHBexc_filt_down1, L_FRAME32k, excTmp );

    flip_spectrum_and_decimby4( excTmp, exc4k, L_FRAME16k, mem_genSHBexc_filt_down2, mem_genSHBexc_filt_down3, 0 );

    if( uv_flag )
    {
        /* unvoiced signal */
        create_random_vector( exc4kWhtnd, L_FRAME16k/4, bwe_seed );
    }
    else
    {
        autocorr(exc4k, R, LPC_WHTN_ORDER_WB+1, L_FRAME16k/4, win_flatten_4k, 0, 1, 1);

        /* Ensure R[0] isn't zero when entering Levinson Durbin */
        R[0] = max(R[0], 1.0e-8f);
        for( i = 0; i <= LPC_WHTN_ORDER_WB; i++ )
        {
            R[i] = R[i] * wac[i];
        }
        lev_dur( lpc_whtn, R, LPC_WHTN_ORDER_WB, ervec );

        fir( exc4k, lpc_whtn, exc4kWhtnd, wht_fil_mem, L_FRAME16k/4, LPC_WHTN_ORDER_WB, 0);

        /* Ensure pow1 is greater than zero when computing normalization */
        for( i=0, pow1=0.00001f; i<L_FRAME16k/4; i++ )
        {
            excTmp2[i] = (float)(fabs(exc4kWhtnd[i]));
            pow1 += exc4kWhtnd[i] * exc4kWhtnd[i];
        }

        for( i=0; i<L_FRAME16k/4; i++ )
        {
            excNoisyEnv[i] = *mem_csfilt + csfilt_num2[0] * excTmp2[i];
            *mem_csfilt = -csfilt_den2[1] * excNoisyEnv[i];
        }

        create_random_vector(exc4k, L_FRAME16k/4, bwe_seed);

        /* Ensure pow2 is greater than zero when computing normalization */
        for( i=0, pow2=0.00001f; i<L_FRAME16k/4; i++ )
        {
            exc4k[i] *= excNoisyEnv[i];
            pow2 += exc4k[i] * exc4k[i];
        }

        if( coder_type == UNVOICED || ( igf_flag && avg_voice_fac < 0.2f) )
        {
            scale = sqrt(pow1/pow2);
            if ((pow2)==0) scale = 0;

            for( i=0; i<L_FRAME16k/4; i++ )
            {
                exc4kWhtnd[i] = exc4k[i] * scale;
            }
        }
        else
        {
            for( i=0, k=0; i<4; i++ )
            {

                if( igf_flag && coder_type == VOICED )
                {
                    tmp_vfac = 2*voice_factors[i];
                    tmp_vfac = min(1, tmp_vfac);
                }
                else
                {
                    tmp_vfac = voice_factors[i];

                }

                temp1 = root_a( tmp_vfac );
                temp2 = root_a_over_b( pow1 * (1.0f - tmp_vfac), pow2 );


                for( j=0; j<L_FRAME16k/16; j++,k++ )
                {
                    exc4kWhtnd[k] = temp1 * exc4kWhtnd[k] + temp2 * exc4k[k];
                }
            }
        }
    }

    syn_filt( lpc_shb, LPC_SHB_ORDER_WB, exc4kWhtnd, excSHB, L_FRAME16k/4, state_lpc_syn, 1 );

    return;
}

/*-------------------------------------------------------------------*
 * GenWBSynth()
 *
 * Generate 16 KHz sampled highband component from synthesized highband
 *-------------------------------------------------------------------*/

void GenWBSynth(
    const float *input_synspeech,           /* i  : input synthesized speech    */
    float *shb_syn_speech_16k,        /* o  : output highband compnent    */
    float *state_lsyn_filt_shb1,      /* i/o: memory                      */
    float *state_lsyn_filt_shb2       /* i/o: memory                      */
)
{
    float speech_buf_16k1[L_FRAME16k], speech_buf_16k2[L_FRAME16k];

    Interpolate_allpass_steep( input_synspeech, state_lsyn_filt_shb1, L_FRAME16k / 4, speech_buf_16k1);

    Interpolate_allpass_steep( speech_buf_16k1, state_lsyn_filt_shb2, L_FRAME16k / 2, speech_buf_16k2);

    flip_spectrum( speech_buf_16k2, shb_syn_speech_16k, L_FRAME16k );

    return;
}

/*-------------------------------------------------------------------*
 * PostShortTerm()
 *
 * Short term processing
 *-------------------------------------------------------------------*/

void PostShortTerm(
    float *sig_in,             /* i  : input signal (pointer to current subframe */
    float *lpccoeff,           /* i  : LPC coefficients for current subframe */
    float *sig_out,            /* o  : postfiltered output */
    float *mem_stp,            /* i/o: postfilter memory*/
    float *ptr_mem_stp,        /* i/o: pointer to postfilter memory*/
    float *ptr_gain_prec,      /* i/o: for gain adjustment*/
    float *mem_zero,           /* i/o: null memory to compute h_st*/
    const float formant_fac          /* i  : Strength of post-filter [0,1] */
)
{
    float apond1[LPC_SHB_ORDER+1];        /* denominator coeff.*/
    float apond2[LONG_H_ST];              /* numerator coeff.  */
    float sig_ltp[L_SUBFR16k+1];          /* residual signal   */
    float parcor0;
    float g1,g2;

    set_f( apond1, 0, LPC_SHB_ORDER+1 );
    set_f( apond2, 0, LONG_H_ST );
    set_f( sig_ltp, 0, L_SUBFR16k+1 );

    /* Obtain post-filter weights  */
    g1 = GAMMA0 + GAMMA_SHARP * formant_fac;
    g2 = GAMMA0 - GAMMA_SHARP * formant_fac;

    /* Compute weighted LPC coefficients */
    weight_a( lpccoeff, apond1, g1, LPC_SHB_ORDER );
    weight_a( lpccoeff, apond2, g2, LPC_SHB_ORDER );

    /* Compute A(gamma2) residual */
    residu( apond2, LPC_SHB_ORDER, sig_in, sig_ltp+1, L_SUBFR16k );

    /* Save last output of 1/A(gamma1)  */
    sig_ltp[0] = *ptr_mem_stp;

    /* Control short term pst filter gain and compute parcor0   */
    calc_st_filt( apond2, apond1, &parcor0, sig_ltp+1, mem_zero, L_SUBFR16k, SWB_TBE );

    /* 1/A(gamma1) filtering, mem_stp is updated */
    syn_filt( apond1, LPC_SHB_ORDER,sig_ltp+1, sig_ltp+1, L_SUBFR16k, mem_stp, 1 );

    /* (1 + mu z-1) tilt filtering */
    filt_mu( sig_ltp, sig_out, parcor0, L_SUBFR16k, SWB_TBE );

    /* gain control */
    scale_st( sig_in, sig_out, ptr_gain_prec, L_SUBFR16k, SWB_TBE );

    return;
}

/*-------------------------------------------------------------------*
 * swb_formant_fac()
 *
 * Find strength of adaptive formant postfilter using tilt of the high
 * band. The 2nd lpc coefficient is used as a tilt approximation.
 *-------------------------------------------------------------------*/

float swb_formant_fac(       /* o  : Formant filter strength [0,1]    */
    const float lpc_shb2,    /* i  : 2nd HB LPC coefficient           */
    float *tilt_mem          /* i/o: Tilt smoothing memory            */
)
{
    float formant_fac;
    float tmp;

    /* Smoothen tilt value */
    tmp = 0.5f * (float)fabs(lpc_shb2) + 0.5f * *tilt_mem;
    *tilt_mem = tmp;

    /* Map to PF strength */
    formant_fac = (tmp - SWB_TILT_LOW)*SWB_TILT_DELTA;
    if (formant_fac > 1.0f)
    {
        formant_fac = 1.0f;
    }
    else if (formant_fac < 0.0f)
    {
        formant_fac = 0.0f;
    }

    formant_fac = 1.0f - 0.5f*formant_fac;

    return formant_fac;
}

/*-------------------------------------------------------------------*
 * GenShapedSHBExcitation()
 *
 * Synthesize spectrally shaped highband excitation signal
 *-------------------------------------------------------------------*/

void GenShapedSHBExcitation(
    float *excSHB,                     /* o  : synthesized shaped shb excitation */
    const float *lpc_shb,                    /* i  : lpc coefficients */
    float *White_exc16k_FB,            /* o  : white excitation for the Fullband extension */
    float *mem_csfilt,                 /* i/o: memory */
    float *mem_genSHBexc_filt_down_shb,/* i/o: memory */
    float *state_lpc_syn,              /* i/o: memory */
    const short coder_type,                  /* i  : coding type */
    const float *bwe_exc_extended,           /* i  : bandwidth extended excitation */
    short bwe_seed[],                  /* i/o: random number generator seed */
    float voice_factors[],             /* i  : voicing factor*/
    const short extl,                        /* i  : extension layer */
    float *tbe_demph,                  /* i/o: de-emphasis memory                      */
    float *tbe_premph,                 /* i/o: pre-emphasis memory                     */
    float *lpc_shb_sf,                 /* i:   LP coefficients                         */
    float *shb_ener_sf,
    float *shb_res_gshape,
    float *shb_res,
    short *vf_ind,
    const float formant_fac,                  /* i   : Formant sharpening factor [0..1] */
    float fb_state_lpc_syn[],           /* i/o: memory */
    float *fb_tbe_demph,                /* i/o: fb de-emphasis memory                   */
    const long bitrate,                       /* i  : bitrate                                 */
    const short prev_bfi                      /* i  : previous frame was concealed            */
)
{
    short i, j, k;
    float wht_fil_mem[LPC_WHTN_ORDER];
    float lpc_whtn[LPC_WHTN_ORDER + 1];
    float R[LPC_WHTN_ORDER + 2];
    float exc32k[L_FRAME32k], exc16k[L_FRAME16k];
    float pow1, pow2, scale, temp1, temp2;
    float excTmp2[L_FRAME16k];
    short nbSubFr;
    float excNoisyEnv[ L_FRAME16k];
    float csfilt_num2[1] = {0.2f};
    float csfilt_den2[2] = {1.0f, -0.8f};
    float varEnvShape;
    float ervec[LPC_WHTN_ORDER+2];
    float exc16kWhtnd[L_FRAME16k];
    float temp = 0.0f;
    float *White_exc16k;
    float voiceFacEst[NB_SUBFR16k];
    float syn_shb_ener_sf[4], tempSHB[80];
    float zero_mem[LPC_SHB_ORDER];
    float vf_tmp;
    float White_exc16k_FB_temp[L_FRAME16k];
    float fb_deemph_fac = 0.48f;

    set_f( zero_mem, 0, LPC_SHB_ORDER);
    set_f( wht_fil_mem, 0, LPC_WHTN_ORDER );

    for(i = 0; i < L_FRAME32k; i++)
    {
        exc32k[i] = ((i%2)==0)?(-bwe_exc_extended[i]):(bwe_exc_extended[i]);
    }

    /* Decimate by 2 */
    Decimate_allpass_steep( exc32k, mem_genSHBexc_filt_down_shb, 2*L_FRAME16k, exc16k );

    autocorr( exc16k, R, LPC_WHTN_ORDER+1, L_FRAME16k, win_flatten, 0, 1, 1 );

    /* Ensure R[0] isn't zero when entering Levinson-Durbin */
    R[0] = max(R[0], 1.0e-8f);
    for( i = 0; i <= LPC_WHTN_ORDER; i++ )
    {
        R[i] = R[i] * wac[i];
    }

    /* Ensure R[0] isn't zero when entering Levinson-Durbin */
    R[0] += 1.0e-8f;
    lev_dur( lpc_whtn, R, LPC_WHTN_ORDER, ervec );

    fir( exc16k, lpc_whtn, exc16kWhtnd, wht_fil_mem, L_FRAME16k, LPC_WHTN_ORDER, 0 );

    if( bitrate >= ACELP_24k40 )
    {
        for(i = 0; i < L_FRAME16k; i++)
        {
            exc16kWhtnd[i] *= shb_res_gshape[(short)(i/80)];
        }
    }

    for( k=0, pow1=0.00001f; k<L_FRAME16k; k++ )
    {
        excTmp2[k ] = (float)(fabs(exc16kWhtnd[k]));
        pow1 += exc16kWhtnd[k] * exc16kWhtnd[k];
    }
    if( bitrate <= ACELP_13k20 && bitrate >= ACELP_7k20)
    {
        varEnvShape = mean(voice_factors, 4);
    }
    else
    {
        varEnvShape = mean(voice_factors, 5);
    }

    if ( extl == FB_TBE)
    {
        fb_deemph_fac = max((0.68f - (float)pow(varEnvShape, 3)), 0.48f);
    }

    varEnvShape = 1.09875f - 0.49875f * varEnvShape;
    varEnvShape = min( max(varEnvShape, 0.6f), 0.999f);
    csfilt_num2[0] = 1.0f - varEnvShape;
    csfilt_den2[1] = - varEnvShape;

    if (*mem_csfilt == 0 && ( bitrate == ACELP_9k60 || bitrate == ACELP_16k40 || bitrate == ACELP_24k40 ) )
    {
        /* pre-init smoothing avoid energy drop outs */
        float tmp_scale = 0;
        for (i=0; i<L_SUBFR16k/4; i++)
        {
            tmp_scale += excTmp2[i];
        }

        /* don't apply for FB in case the FB start-frame was potentially lost - White_exc16k is very sensitive to enery mismatch between enc - dec */
        /* rather stick to the more conservative approach, to avoid potential clippings */
        if( !(prev_bfi && extl == FB_TBE) )
        {
            /* use weak smoothing for 1st frame after switching to make filter recover more quickly */
            varEnvShape = 0.8f;
            csfilt_num2[0] = 1.0f - varEnvShape;
            csfilt_den2[1] = - varEnvShape;
        }
        *mem_csfilt = varEnvShape*(tmp_scale/(L_SUBFR16k/4));
    }

    /* Track the low band envelope */
    for( k = 0; k < L_FRAME16k; k++ )
    {
        excNoisyEnv[k] = *mem_csfilt + csfilt_num2[0] * excTmp2[k];
        *mem_csfilt = -csfilt_den2[1] * excNoisyEnv[k];
    }

    White_exc16k = exc16k;

    create_random_vector( White_exc16k, 256, bwe_seed );
    create_random_vector( White_exc16k + 256, L_FRAME16k - 256, bwe_seed );

    for( k=0, pow2=0.00001f; k<L_FRAME16k; k++ )
    {
        White_exc16k[k] *= excNoisyEnv[k];
        pow2 += White_exc16k[k] * White_exc16k[k];
    }

    if( bitrate >= ACELP_24k40 )
    {
        if( *vf_ind == 20 ) /* encoder side */
        {
            Estimate_mix_factors( shb_res, exc16kWhtnd, White_exc16k, pow1, pow2, voiceFacEst, vf_ind );
            temp = (voiceFacEst[0] > 0.7f)? 1.0f : 0.8f;
        }
        else /* decoder side */
        {
            temp = ((*vf_ind * 0.125f) > 0.7f)? 1.0f : 0.8f;
        }
        for(i = 0; i < NB_SUBFR16k; i++)
        {
            voice_factors[i] *=  temp;
        }
    }

    mvr2r( White_exc16k, White_exc16k_FB, L_FRAME16k );
    deemph( White_exc16k, PREEMPH_FAC, L_FRAME16k, tbe_demph );

    if( coder_type == UNVOICED )
    {
        scale = sqrt( pow1/pow2 );
        if ((pow2)==0) scale = 0;
        for( k=0; k<L_FRAME16k; k++ )
        {
            exc16kWhtnd[k] = White_exc16k[k] * scale;
        }

        preemph( exc16kWhtnd, PREEMPH_FAC, L_FRAME16k, tbe_premph );
    }
    else
    {
        nbSubFr = ( bitrate < ACELP_24k40 )? NB_SUBFR : NB_SUBFR16k;
        for( i = 0, k = 0; i < nbSubFr; i++ )
        {
            if(coder_type == VOICED && (bitrate < ACELP_24k40) )
            {
                temp = sqrt( voice_factors[i] );
                temp1 = sqrt(temp);
                temp2 = sqrt( (pow1 * (1.0f - temp))/pow2 );
                if ((pow2)==0) temp2 = 0;
            }
            else
            {
                /* Adjust noise mixing for formant sharpening filter */
                vf_tmp = SWB_NOISE_MIX_FAC * formant_fac;
                vf_tmp = voice_factors[i] * (1.0f - vf_tmp);

                temp1 = sqrt(vf_tmp);
                temp2 = sqrt((pow1 * (1.0f - vf_tmp))/pow2);
                if ((pow2)==0) temp2 = 0;
            }

            for( j=0; j<L_FRAME16k/nbSubFr; j++, k++ )
            {
                exc16kWhtnd[k] = temp1 * exc16kWhtnd[k] + temp2 * White_exc16k[k];
            }

            temp = sqrt( 1.0f - voice_factors[i] );
            temp = PREEMPH_FAC*temp/(temp1+temp);

            preemph( &exc16kWhtnd[i*L_FRAME16k/nbSubFr], temp, L_FRAME16k/nbSubFr, tbe_premph );
        }
    }

    if ( bitrate < ACELP_24k40 )
    {
        syn_filt( lpc_shb, LPC_SHB_ORDER, exc16kWhtnd, excSHB, L_FRAME16k, state_lpc_syn, 1 );
    }
    else
    {
        set_f( zero_mem, 0, LPC_SHB_ORDER);
        syn_filt( lpc_shb_sf,                      LPC_SHB_ORDER, exc16kWhtnd,    tempSHB, 80, zero_mem, 1 );
        syn_shb_ener_sf[0] = 0.125f * sum2_f(tempSHB, 80);
        syn_filt( lpc_shb_sf+(LPC_SHB_ORDER+1),   LPC_SHB_ORDER, exc16kWhtnd+ 80, tempSHB, 80, zero_mem, 1 );
        syn_shb_ener_sf[1] = 0.125f * sum2_f(tempSHB, 80);
        syn_filt( lpc_shb_sf+2*(LPC_SHB_ORDER+1), LPC_SHB_ORDER, exc16kWhtnd+160, tempSHB, 80, zero_mem, 1 );
        syn_shb_ener_sf[2] = 0.125f * sum2_f(tempSHB, 80);
        syn_filt( lpc_shb_sf+3*(LPC_SHB_ORDER+1), LPC_SHB_ORDER, exc16kWhtnd+240, tempSHB, 80, zero_mem, 1 );
        syn_shb_ener_sf[3] = 0.125f * sum2_f(tempSHB, 80);
        if(bitrate <= ACELP_32k)
        {
            tempSHB[0] = (float)(shb_ener_sf[0])/(syn_shb_ener_sf[0]+syn_shb_ener_sf[1]+syn_shb_ener_sf[2]+syn_shb_ener_sf[3]) ;
            for(i = 0; i < L_FRAME16k; i++)
            {
                exc16kWhtnd[i] = exc16kWhtnd[i] * sqrt(tempSHB[0]);
            }
        }




        syn_filt( lpc_shb_sf,                     LPC_SHB_ORDER, exc16kWhtnd,     excSHB,     80, state_lpc_syn, 1 );
        syn_filt( lpc_shb_sf+(LPC_SHB_ORDER+1),   LPC_SHB_ORDER, exc16kWhtnd+ 80, excSHB+ 80, 80, state_lpc_syn, 1 );
        syn_filt( lpc_shb_sf+2*(LPC_SHB_ORDER+1), LPC_SHB_ORDER, exc16kWhtnd+160, excSHB+160, 80, state_lpc_syn, 1 );
        syn_filt( lpc_shb_sf+3*(LPC_SHB_ORDER+1), LPC_SHB_ORDER, exc16kWhtnd+240, excSHB+240, 80, state_lpc_syn, 1 );
    }

    if ( extl == FB_TBE)
    {

        syn_filt( lpc_shb, LPC_SHB_ORDER, White_exc16k_FB, White_exc16k_FB_temp, L_FRAME16k, fb_state_lpc_syn, 1 );

        for( i=0; i<L_FRAME16k; i++ )
        {
            White_exc16k_FB_temp[i] *= cos_fb_exc[i%32];
        }

        flip_spectrum( White_exc16k_FB_temp, White_exc16k_FB, L_FRAME16k );

        deemph( White_exc16k_FB, fb_deemph_fac, L_FRAME16k, fb_tbe_demph );
    }
    else
    {
        for( i=0; i<L_FRAME16k; i++ )
        {
            White_exc16k_FB[i] = 0.0f;
        }
    }

    return;
}

/*-------------------------------------------------------------------*
 * GenSHBSynth()
 *
 * Generate 32 KHz sampled highband component from synthesized highband
 *-------------------------------------------------------------------*/

void GenSHBSynth(
    const float *input_synspeech,               /* i  : input synthesized speech    */
    float *shb_syn_speech_32k,            /* o  : output highband component   */
    float Hilbert_Mem[],                  /* i/o: memory                      */
    float state_lsyn_filt_shb_local[],    /* i/o: memory                      */
    const short L_frame,                        /* i   : ACELP frame length         */
    short *syn_dm_phase
)
{
    float speech_buf_32k[L_FRAME32k];
    short i;

    Interpolate_allpass_steep( input_synspeech, state_lsyn_filt_shb_local, L_FRAME16k, speech_buf_32k );

    if( L_frame == L_FRAME )
    {
        flip_and_downmix_generic( speech_buf_32k, shb_syn_speech_32k, L_FRAME32k, Hilbert_Mem,
                                  Hilbert_Mem + HILBERT_ORDER1, Hilbert_Mem + (HILBERT_ORDER1+2*HILBERT_ORDER2), syn_dm_phase );
    }
    else
    {
        for( i = 0; i < L_FRAME32k; i++ )
        {
            shb_syn_speech_32k[i] = ((i%2)==0)?(-speech_buf_32k[i]):(speech_buf_32k[i]);
        }
    }

    return;
}

/*-------------------------------------------------------------------*
 * ScaleShapedSHB()
 *
 *
 *-------------------------------------------------------------------*/

void ScaleShapedSHB(
    const short length,             /* i  : SHB overlap length      */
    float *synSHB,            /* i/o: synthesized shb signal  */
    float *overlap,           /* i/o: buffer for overlap-add  */
    const float *subgain,           /* i  : subframe gain           */
    const float frame_gain,         /* i  : frame gain              */
    const float *win,               /* i  : window                  */
    const float *subwin             /* i  : subframes window        */
)
{
    const short *skip;
    short i, j, k, l_shb_lahead, l_frame;
    short join_length, num_join;
    float mod_syn[L_FRAME16k+L_SHB_LAHEAD], sum_gain;

    /* initilaization */
    l_frame = L_FRAME16k;
    l_shb_lahead = L_SHB_LAHEAD;
    skip = skip_bands_SWB_TBE;

    if( length == SHB_OVERLAP_LEN/2 )
    {
        skip = skip_bands_WB_TBE;
        l_frame = L_FRAME16k/4;
        l_shb_lahead = L_SHB_LAHEAD/4;
    }

    /* apply gain for each subframe, and store noise output signal using overlap-add */
    set_f( mod_syn, 0, l_frame+l_shb_lahead );

    if( length == SHB_OVERLAP_LEN/2 )
    {
        sum_gain = 0;
        for( k=0; k<length/2; k++ )
        {
            sum_gain = subwin[2*k+2]*subgain[0];
            mod_syn[skip[0]+k] = synSHB[skip[0]+k] * sum_gain;
            mod_syn[skip[0]+k+length/2] = synSHB[skip[0]+k+length/2] * subgain[0];
        }
        for( i=1; i<NUM_SHB_SUBFR/2; i++ )
        {
            for( k=0; k<length; k++ )
            {
                sum_gain = subwin[k+1]*subgain[i] + subwin[length-k-1]*subgain[i-1];
                mod_syn[skip[i]+k] = synSHB[skip[i]+k] * sum_gain;
            }
        }
        for( k=0; k<length/2; k++ )
        {
            sum_gain = subwin[length-2*k-2]*subgain[i-1];
            mod_syn[skip[i]+k] = synSHB[skip[i]+k] * sum_gain;
        }
    }
    else
    {
        num_join = NUM_SHB_SUBFR/NUM_SHB_SUBGAINS;
        join_length = num_join*length;
        for (k = 0, j = 0; k < length; k++)
        {
            mod_syn[j] = synSHB[j]*subwin[k+1]*subgain[0];
            j++;
        }
        for (i = 0; i < NUM_SHB_SUBGAINS-1; i++)
        {
            for (k = 0; k < join_length - length; k++)
            {
                mod_syn[j] = synSHB[j]*subgain[i*num_join];
                j++;
            }

            for (k = 0; k < length; k++)
            {
                mod_syn[j] = synSHB[j]*(subwin[length-k-1]*subgain[i*num_join] + subwin[k+1]*subgain[(i+1)*num_join]);
                j++;
            }
        }
        for (k = 0; k < join_length - length; k++)
        {
            mod_syn[j] = synSHB[j]*subgain[(NUM_SHB_SUBGAINS-1)*num_join];
            j++;
        }
        for (k = 0; k < length; k++)
        {
            mod_syn[j] = synSHB[j]*subwin[length-k-1]*subgain[(NUM_SHB_SUBGAINS-1)*num_join];
            j++;
        }
    }

    for( i=0; i<l_shb_lahead; i++ )
    {
        synSHB[i] = mod_syn[i] * win[i] * frame_gain;
        synSHB[i] += overlap[i];
        synSHB[i+l_shb_lahead] = mod_syn[i] * frame_gain;
    }

    for( ; i<l_frame; i++ )
    {
        synSHB[i] = mod_syn[i] * frame_gain;
    }

    for( ; i<l_frame+l_shb_lahead; i++ )
    {
        overlap[i-l_frame] = mod_syn[i] * win[l_frame+l_shb_lahead-1-i] * frame_gain;
    }

    return;
}

/*-------------------------------------------------------------------*
 * non_linearity()
 *
 * Apply a non linearity to the SHB excitation
 * -------------------------------------------------------------------*/

void non_linearity(
    const float input[],                /* i  : input signal    */
    float output[],               /* o  : output signal   */
    float old_bwe_exc_extended[], /* i/o: memory bugffer  */
    const short length,                 /* i  : input length    */
    float *prev_scale             /* i/o: memory          */
    ,short  coder_type,           /* i  : Coder Type          */
    float   *voice_factors,       /* i  : Voice Factors       */
    const short   L_frame			    /* i  : ACELP frame length  */
)
{
    short i,j;
    float max = 0.0;
    float scale, temp;
    float scale_step;
    float *p_out;

    short en_abs = 0;
    float v_fac = 0, ths;
    short nframes;


    if (L_frame == L_FRAME16k)
    {
        nframes = 5;
        ths = 0.87f;
    }
    else
    {
        nframes = 4;
        ths = 0.94f;
    }


    for (i=0; i<nframes; i++)
    {
        v_fac += voice_factors[i];

    }
    v_fac /= nframes;

    if (coder_type == VOICED && v_fac > ths )
    {
        en_abs = 1;
    }





    p_out = output + NL_BUFF_OFFSET;   /* NL_BUFF_OFFSET = 12 */
    /* update buffer memory */
    mvr2r( old_bwe_exc_extended, output, NL_BUFF_OFFSET );

    for (i=j=0; i<length/2; i++)
    {
        if ((temp = (float) fabs(input[i])) > max)
        {
            max = temp;
            j = i;
        }
    }

    if (max > 1.0f)
    {
        scale = 0.67f / max;
    }
    else
    {
        scale = 0.67f;
    }


    if ( *prev_scale <= 0.0 || *prev_scale > 1024.0f * scale )
    {
        scale_step = 1.0;
        *prev_scale = scale;
    }
    else
    {
        scale_step = 1.0f;
        if(j != 0)
        {
            scale_step = (float) exp(1.0f / (float) j * (float) log(scale / *prev_scale));
        }
    }

    for (i=0; i<length/2; i++)
    {
        if (input[i] >= 0.0)
        {
            *p_out++ = (input[i] * input[i]) **prev_scale;
        }
        else
        {
            if (en_abs)
            {
                *p_out++ =  1.0f * (input[i] * input[i]) **prev_scale;
            }
            else
            {
                *p_out++ =  -1.0f * (input[i] * input[i]) **prev_scale;
            }

        }

        if (i < j)
        {
            *prev_scale *= scale_step;
        }
    }

    max = 0.0f;
    for (i=j=length/2; i<length; i++)
    {
        if ((temp = (float) fabs(input[i])) > max)
        {
            max = temp;
            j = i;
        }
    }

    if (max > 1.0f)
    {
        scale = 0.67f / max;
    }
    else
    {
        scale = 0.67f;
    }


    if ( *prev_scale <= 0.0 || *prev_scale > 1024.0f * scale )
    {
        scale_step = 1.0;
        *prev_scale = scale;
    }
    else
    {
        scale_step = 1.0f;
        if(j != length/2)
        {
            scale_step = (float) exp(1.0f / (float) (j - length/2) * (float) log(scale / *prev_scale));
        }
    }

    for (i=length/2; i<length; i++)
    {
        if (input[i] >= 0.0)
        {
            *p_out++ = (input[i] * input[i]) **prev_scale;
        }
        else
        {
            if (en_abs)
            {
                *p_out++ =  1.0f * (input[i] * input[i]) **prev_scale;
            }
            else
            {
                *p_out++ =  -1.0f * (input[i] * input[i]) **prev_scale;
            }
        }

        if (i < j)
        {
            *prev_scale *= scale_step;
        }
    }

    /* update buffer memory */
    mvr2r( output + L_FRAME32k, old_bwe_exc_extended, NL_BUFF_OFFSET );

    return;
}


/*-------------------------------------------------------------------*
 * create_random_vector()
 *
 * creates random number vector
 * -------------------------------------------------------------------*/

void create_random_vector(
    float output[],         /* o  : output random vector      */
    const short length,           /* i  : length of random vector   */
    short seed[]            /* i/o: start seed                */
)
{
    short i, j, k;
    float scale1, scale2;

    j = (short) (own_random(&seed[0]) * 0.0078f);
    j = abs(j) & 0xff;
    k = (short) (own_random(&seed[1]) * 0.0078f);
    k = abs(k) & 0xff;

    while( k==j )
    {
        k = (short)(own_random(&seed[1]) * 0.0078f);
        k = abs(k) & 0xff;
    }

    if( own_random(&seed[0]) < 0 )
    {
        scale1 = -563.154f;     /* -200.00f * 0.35f/0.1243f; */
    }
    else
    {
        scale1 = 563.154f;      /* 200.00f * 0.35f/0.1243f; */
    }

    if( own_random(&seed[1]) < 0 )
    {
        scale2 = -225.261f;     /* -80.00f * 0.35f/0.1243f; */
    }
    else
    {
        scale2 = 225.261f;      /* 80.00f * 0.35f/0.1243f; */
    }

    for( i=0; i<length; i++, j++, k++ )
    {
        j &= 0xff;
        k &= 0xff;
        output[i] = scale1 * gaus_dico_swb[j] + scale2 * gaus_dico_swb[k];
    }

    return;
}


/*-------------------------------------------------------------------*
 * interp_code_5over2()
 *
 * Used to interpolate the excitation from the core sample rate
 * of 12.8 kHz to 32 kHz.
 * Simple linear interpolator - No need for precision here.
 *-------------------------------------------------------------------*/

void interp_code_5over2(
    const float inp_code[],         /* i  : input vector                */
    float interp_code[],      /* o  : output vector               */
    const short inp_length          /* i  : length of input vector      */
)
{
    short i, kk, kkp1;
    float factor_i[5] = {0.2f, 0.6f, 1.0f, 0.6f, 0.2f};
    float factor_j[5] = {0.8f, 0.4f, 0.0f, 0.4f, 0.8f};

    interp_code[0] = inp_code[0];
    interp_code[1] = inp_code[0] * factor_i[3] + inp_code[1] * factor_j[3];
    interp_code[2] = inp_code[0] * factor_i[4] + inp_code[1] * factor_j[4];

    for (i=3, kk=1, kkp1=2; i < (inp_length - 2) * HIBND_ACB_L_FAC; i+=5, kk++, kkp1++)
    {
        interp_code[i] = inp_code[kk] * factor_j[0] + inp_code[kkp1] * factor_i[0];
        interp_code[i+1] = inp_code[kk] * factor_j[1] + inp_code[kkp1] * factor_i[1];
        interp_code[i+2] = inp_code[kkp1] * factor_i[2];
        kk++;
        kkp1++;
        interp_code[i+3] = inp_code[kk] * factor_i[3] + inp_code[kkp1] * factor_j[3];
        interp_code[i+4] = inp_code[kk] * factor_i[4] + inp_code[kkp1] * factor_j[4];
    }

    interp_code[i] = inp_code[kk] * factor_j[0];
    interp_code[i+1] = inp_code[kk] * factor_j[1];

    return;
}

/*-------------------------------------------------------------------*
 * interp_code_4over2()
 *
 * Used to interpolate the excitation from the core sample rate
 * of 16 kHz to 32 kHz.
 * Simple linear interpolator - No need for precision here.
 *-------------------------------------------------------------------*/

void interp_code_4over2(
    const float inp_code[],         /* i  : input vector                */
    float interp_code[],      /* o  : output vector               */
    const short inp_length          /* i  : length of input vector      */
)
{
    short   i,j;
    for (i=j=0; i<inp_length-1; i++, j+=2)
    {
        interp_code[j] = inp_code[i];
        interp_code[j+1] = inp_code[i] * 0.5f + inp_code[i+1] * 0.5f;
    }

    interp_code[j] = inp_code[i];
    interp_code[j+1] = inp_code[i] * 0.5f;

    return;
}

/*-------------------------------------------------------------------*
 * fb_tbe_reset_synth()
 *
 * Reset the extra parameters needed for synthesis of the FB TBE output
 *-------------------------------------------------------------------*/

void fb_tbe_reset_synth(
    float fbbwe_hpf_mem[][4],
    float *prev_fbbwe_ratio
)
{
    set_f( fbbwe_hpf_mem[0], 0, 4 );
    set_f( fbbwe_hpf_mem[1], 0, 4 );
    set_f( fbbwe_hpf_mem[2], 0, 4 );
    set_f( fbbwe_hpf_mem[3], 0, 4 );
    *prev_fbbwe_ratio = 1.0f;

    return;
}

/*-------------------------------------------------------------------*
 * wb_tbe_extras_reset()
 *
 * Reset the extra parameters only required for WB TBE encoding
 *-------------------------------------------------------------------*/

void wb_tbe_extras_reset(
    float mem_genSHBexc_filt_down_wb2[],
    float mem_genSHBexc_filt_down_wb3[]
)
{
    set_f( mem_genSHBexc_filt_down_wb2, 0.0f, (2*ALLPASSSECTIONS_STEEP+1) );
    set_f( mem_genSHBexc_filt_down_wb3, 0.0f, (2*ALLPASSSECTIONS_STEEP+1) );

    return;
}

/*-------------------------------------------------------------------*
 * wb_tbe_extras_reset_synth()
 *
 * Reset the extra parameters only required for WB TBE synthesis
 *-------------------------------------------------------------------*/

void wb_tbe_extras_reset_synth(
    float state_lsyn_filt_shb[],
    float state_lsyn_filt_dwn_shb[],
    float mem_resamp_HB[]
)
{
    set_f( state_lsyn_filt_shb, 0.0f, 2 * ALLPASSSECTIONS_STEEP );
    set_f( state_lsyn_filt_dwn_shb, 0.0f, 2 * ALLPASSSECTIONS_STEEP );
    set_f( mem_resamp_HB, 0.0f, INTERP_3_1_MEM_LEN );

    return;
}

/*-------------------------------------------------------------------*
 * elliptic_bpf_48k_generic()
 *
 * 18th-order elliptic bandpass filter at 14.0 to 20 kHz sampled at 48 kHz
 * Implemented as 3 fourth order sections cascaded.
 *-------------------------------------------------------------------*/

void elliptic_bpf_48k_generic(
    const float input[],                /* i  : input signal                            */
    float output[],               /* o  : output signal                           */
    float memory[][4],            /* i/o: 4 arrays of 4 for memory                */
    const float full_band_bpf[][5]      /* i  : filter coefficients b0,b1,b2,a0,a1,a2   */
)
{
    short i;
    float tmp[L_FRAME48k], tmp2[L_FRAME48k];

    tmp[0] = memory[0][0] * full_band_bpf[0][4] + memory[0][1] * full_band_bpf[0][3] + memory[0][2] * full_band_bpf[0][2] + memory[0][3] * full_band_bpf[0][1] + input[0] * full_band_bpf[0][0]
             - full_band_bpf[3][1] * memory[1][3] - full_band_bpf[3][2] * memory[1][2] - full_band_bpf[3][3] * memory[1][1] - full_band_bpf[3][4] * memory[1][0];
    tmp[1] = memory[0][1] * full_band_bpf[0][4] + memory[0][2] * full_band_bpf[0][3] + memory[0][3] * full_band_bpf[0][2] + input[0] * full_band_bpf[0][1] + input[1] * full_band_bpf[0][0]
             - full_band_bpf[3][1] * tmp[0] - full_band_bpf[3][2] * memory[1][3] - full_band_bpf[3][3] * memory[1][2] - full_band_bpf[3][4] * memory[1][1];
    tmp[2] = memory[0][2] * full_band_bpf[0][4] + memory[0][3] * full_band_bpf[0][3] + input[0] * full_band_bpf[0][2] + input[1] * full_band_bpf[0][1] + input[2] * full_band_bpf[0][0]
             - full_band_bpf[3][1] * tmp[1] - full_band_bpf[3][2] * tmp[0] - full_band_bpf[3][3] * memory[1][3] - full_band_bpf[3][4] * memory[1][2];
    tmp[3] = memory[0][3] * full_band_bpf[0][4] + input[0] * full_band_bpf[0][3] + input[1] * full_band_bpf[0][2] + input[2] * full_band_bpf[0][1] + input[3] * full_band_bpf[0][0]
             - full_band_bpf[3][1] * tmp[2] - full_band_bpf[3][2] * tmp[1] - full_band_bpf[3][3] * tmp[0] - full_band_bpf[3][4] * memory[1][3];

    for( i=4; i<L_FRAME48k; i++ )
    {
        tmp[i] = input[i-4] * full_band_bpf[0][4] + input[i-3] * full_band_bpf[0][3] + input[i-2] * full_band_bpf[0][2] + input[i-1] * full_band_bpf[0][1] + input[i] * full_band_bpf[0][0]
                 - full_band_bpf[3][1] * tmp[i-1] - full_band_bpf[3][2] * tmp[i-2] - full_band_bpf[3][3] * tmp[i-3] - full_band_bpf[3][4] * tmp[i-4];
    }

    memory[0][0] = input[L_FRAME48k-4];
    memory[0][1] = input[L_FRAME48k-3];
    memory[0][2] = input[L_FRAME48k-2];
    memory[0][3] = input[L_FRAME48k-1];

    tmp2[0] = memory[1][0] * full_band_bpf[1][4] + memory[1][1] * full_band_bpf[1][3] + memory[1][2] * full_band_bpf[1][2] + memory[1][3] * full_band_bpf[1][1] + tmp[0] * full_band_bpf[1][0]
              - full_band_bpf[4][1] * memory[2][3] - full_band_bpf[4][2] * memory[2][2] - full_band_bpf[4][3] * memory[2][1] - full_band_bpf[4][4] * memory[2][0];
    tmp2[1] = memory[1][1] * full_band_bpf[1][4] + memory[1][2] * full_band_bpf[1][3] + memory[1][3] * full_band_bpf[1][2] + tmp[0] * full_band_bpf[1][1] + tmp[1] * full_band_bpf[1][0]
              - full_band_bpf[4][1] * tmp2[0] - full_band_bpf[4][2] * memory[2][3] - full_band_bpf[4][3] * memory[2][2] - full_band_bpf[4][4] * memory[2][1];
    tmp2[2] = memory[1][2] * full_band_bpf[1][4] + memory[1][3] * full_band_bpf[1][3] + tmp[0] * full_band_bpf[1][2] + tmp[1] * full_band_bpf[1][1] + tmp[2] * full_band_bpf[1][0]
              - full_band_bpf[4][1] * tmp2[1] - full_band_bpf[4][2] * tmp2[0] - full_band_bpf[4][3] * memory[2][3] - full_band_bpf[4][4] * memory[2][2];
    tmp2[3] = memory[1][3] * full_band_bpf[1][4] + tmp[0] * full_band_bpf[1][3] + tmp[1] * full_band_bpf[1][2] + tmp[2] * full_band_bpf[1][1] + tmp[3] * full_band_bpf[1][0]
              - full_band_bpf[4][1] * tmp2[2] - full_band_bpf[4][2] * tmp2[1] - full_band_bpf[4][3] * tmp2[0] - full_band_bpf[4][4] * memory[2][3];

    for( i=4; i<L_FRAME48k; i++ )
    {
        tmp2[i] = tmp[i-4] * full_band_bpf[1][4] + tmp[i-3] * full_band_bpf[1][3] + tmp[i-2] * full_band_bpf[1][2] + tmp[i-1] * full_band_bpf[1][1] + tmp[i] * full_band_bpf[1][0]
                  - full_band_bpf[4][1] * tmp2[i-1] - full_band_bpf[4][2] * tmp2[i-2] - full_band_bpf[4][3] * tmp2[i-3] - full_band_bpf[4][4] * tmp2[i-4];
    }

    memory[1][0] = tmp[L_FRAME48k-4];
    memory[1][1] = tmp[L_FRAME48k-3];
    memory[1][2] = tmp[L_FRAME48k-2];
    memory[1][3] = tmp[L_FRAME48k-1];

    output[0] = memory[2][0] * full_band_bpf[2][4] + memory[2][1] * full_band_bpf[2][3] + memory[2][2] * full_band_bpf[2][2] + memory[2][3] * full_band_bpf[2][1] + tmp2[0] * full_band_bpf[2][0]
                - full_band_bpf[5][1] * memory[3][3] - full_band_bpf[5][2] * memory[3][2] - full_band_bpf[5][3] * memory[3][1] - full_band_bpf[5][4] * memory[3][0];
    output[1] = memory[2][1] * full_band_bpf[2][4] + memory[2][2] * full_band_bpf[2][3] + memory[2][3] * full_band_bpf[2][2] + tmp2[0] * full_band_bpf[2][1] + tmp2[1] * full_band_bpf[2][0]
                - full_band_bpf[5][1] * output[0] - full_band_bpf[5][2] * memory[3][3] - full_band_bpf[5][3] * memory[3][2] - full_band_bpf[5][4] * memory[3][1];
    output[2] = memory[2][2] * full_band_bpf[2][4] + memory[2][3] * full_band_bpf[2][3] + tmp2[0] * full_band_bpf[2][2] + tmp2[1] * full_band_bpf[2][1] + tmp2[2] * full_band_bpf[2][0]
                - full_band_bpf[5][1] * output[1] - full_band_bpf[5][2] * output[0] - full_band_bpf[5][3] * memory[3][3] - full_band_bpf[5][4] * memory[3][2];
    output[3] = memory[2][3] * full_band_bpf[2][4] + tmp2[0] * full_band_bpf[2][3] + tmp2[1] * full_band_bpf[2][2] + tmp2[2] * full_band_bpf[2][1] + tmp2[3] * full_band_bpf[2][0]
                - full_band_bpf[5][1] * output[2] - full_band_bpf[5][2] * output[1] - full_band_bpf[5][3] * output[0] - full_band_bpf[5][4] * memory[3][3];

    for( i=4; i<L_FRAME48k; i++ )
    {
        output[i] = tmp2[i-4] * full_band_bpf[2][4] + tmp2[i-3] * full_band_bpf[2][3] + tmp2[i-2] * full_band_bpf[2][2] + tmp2[i-1] * full_band_bpf[2][1] + tmp2[i] * full_band_bpf[2][0]
                    - full_band_bpf[5][1] * output[i-1] - full_band_bpf[5][2] * output[i-2] - full_band_bpf[5][3] * output[i-3] - full_band_bpf[5][4] * output[i-4];
    }

    memory[2][0] = tmp2[L_FRAME48k-4];
    memory[2][1] = tmp2[L_FRAME48k-3];
    memory[2][2] = tmp2[L_FRAME48k-2];
    memory[2][3] = tmp2[L_FRAME48k-1];

    memory[3][0] = output[L_FRAME48k-4];
    memory[3][1] = output[L_FRAME48k-3];
    memory[3][2] = output[L_FRAME48k-2];
    memory[3][3] = output[L_FRAME48k-1];

    return;
}


/*-------------------------------------------------------------------*
 * synthesise_fb_high_band()
 *
 * Creates the highband output for full band  - 14.0 to 20 kHz
 * Using the energy shaped white excitation signal from the SWB BWE.
 * The excitation signal input is sampled at 16kHz and so is upsampled
 * to 48 kHz first.
 * Uses a complementary split filter to code the two regions from
 * 14kHz to 16kHz and 16 kHz to 20 kHz.
 * One of 16 tilt filters is also applied afterwards to further
 * refine the spectral shape of the fullband signal.
 * The tilt is specified in dB per kHz. N.B. Only negative values are
 * accomodated.
 *-------------------------------------------------------------------*/

void synthesise_fb_high_band(
    const float excitation_in[],    /* i  : full band excitation                                */
    float output[],           /* o  : high band speech - 14.0 to 20 kHz                   */
    const float fb_exc_energy,		/* i  : full band excitation energy                         */
    const float ratio,				/* i  : energy ratio		                                */
    const short L_frame,			/* i  : ACELP frame length                                  */
    const short bfi,                /* i  : fec flag			                                */
    float *prev_fbbwe_ratio,  /* o  : previous frame energy for FEC                       */
    float bpf_memory[][4]     /* i/o: memory for elliptic bpf 48k                         */
)
{
    short i, j;
    float excitation_in_interp3[L_FRAME48k];
    float tmp[L_FRAME48k];
    float temp1, ratio2;

    /* Interpolate the white energy shaped gaussian excitation from 16 kHz to 48 kHz with zeros */
    /* white excitation from DC to 8 kHz resampled to produce DC to 24 kHz excitation.          */
    for( i=0, j=0; i<L_FRAME48k; i+=3, j++ )
    {
        excitation_in_interp3[i] = 3.0f * excitation_in[j];
        excitation_in_interp3[i+1] = 0.0f;
        excitation_in_interp3[i+2] = 0.0f;
    }

    if( L_frame == L_FRAME16k )
    {
        /* for 16kHz ACELP core */
        elliptic_bpf_48k_generic(excitation_in_interp3, tmp, bpf_memory, full_band_bpf_3);
    }
    else
    {
        /* for 12.8kHz ACELP core */
        elliptic_bpf_48k_generic( excitation_in_interp3, tmp, bpf_memory, full_band_bpf_1 );
    }
    temp1 = sum2_f( tmp, L_FRAME48k ) + 0.001f;
    ratio2 = ratio * sqrt( fb_exc_energy / temp1 );

    if( !bfi )
    {
        *prev_fbbwe_ratio = ratio;
    }
    else
    {
        *prev_fbbwe_ratio = ratio*0.5f;
    }
    for( i=0; i<L_FRAME48k; i++ )
    {
        output[i] = tmp[i]*ratio2;
    }

    return;
}

/*-------------------------------------------------------------------*
 * Estimate_mix_factors()                                            *
 *                                                                   *
 * Estimate mix factors for SHB excitation generation                *
 *-------------------------------------------------------------------*/

void Estimate_mix_factors(
    const float *shb_res,                 /* i  : SHB LP residual */
    const float *exc16kWhtnd,             /* i  : SHB transformed low band excitation */
    const float *White_exc16k,            /* i  : Modulated envelope shaped white noise  */
    const float pow1,                     /* i  : SHB exc. power for normalization */
    const float pow2,                     /* i  : White noise excitation for normalization */
    float *vf_modified,             /* o  : Estimated voice factors */
    short *vf_ind                   /* o  : voice factors VQ index */
)
{
    float shb_res_local[L_FRAME16k], WN_exc_local[L_FRAME16k];
    float pow3, temp_p1_p2, temp_p1_p3;
    float temp_numer1[L_FRAME16k], temp_numer2[L_FRAME16k];
    short i, length;

    mvr2r(shb_res, shb_res_local, L_FRAME16k);
    mvr2r(White_exc16k, WN_exc_local, L_FRAME16k);

    pow3 = dotp(shb_res_local, shb_res_local, L_FRAME16k);

    pow3 += 0.00001f;
    temp_p1_p2 = (float)sqrt(pow1/pow2);
    temp_p1_p3 = (float)sqrt(pow1/pow3);


    for(i = 0; i < L_FRAME16k; i++)
    {
        WN_exc_local[i] *= temp_p1_p2;
        shb_res_local[i] *= temp_p1_p3;
    }
    for(i = 0; i < L_FRAME16k; i++)
    {
        temp_numer1[i] = shb_res_local[i] - WN_exc_local[i];
        temp_numer2[i] = exc16kWhtnd[i] - WN_exc_local[i];
    }

    length = L_FRAME16k;
    for(i = 0; i < 1; i++)
    {
        temp_p1_p2 = dotp(temp_numer1+i*length, temp_numer2+i*length, length);
        temp_p1_p3 = dotp(temp_numer2+i*length, temp_numer2+i*length, length);
        vf_modified[i] = min( max( (temp_p1_p2 / temp_p1_p3), 0.1f), 0.99f);
    }

    *vf_ind = usquant(vf_modified[0], &temp_p1_p2, 0.125, 0.125, 1<<NUM_BITS_SHB_VF);
    set_f(vf_modified, temp_p1_p2, NB_SUBFR16k);

    return;
}

/*-------------------------------------------------------------------*
 * prep_tbe_exc()                                                    *
 *                                                                   *
 * Prepare TBE excitation                                            *
 *-------------------------------------------------------------------*/

void prep_tbe_exc(
    const short L_frame,            /* i  : length of the frame         */
    const short i_subfr,            /* i  : subframe index              */
    const float gain_pit,           /* i  : Pitch gain                  */
    const float gain_code,          /* i  : algebraic codebook gain     */
    const float code[],             /* i  : algebraic excitation        */
    const float voice_fac,          /* i  : voicing factor              */
    float *voice_factors,     /* o  : TBE voicing factor          */
    float bwe_exc[],          /* i/o: excitation for TBE          */
    const float gain_preQ,          /* i  : prequantizer excitation gain*/
    const float code_preQ[],        /* i  : prequantizer excitation     */
    const short T0,                 /* i  : integer pitch variables     */
    const short coder_type,         /* i  : coding type                 */
    const long  core_brate          /* i  : core bitrate                */
)
{
    short i;
    float tmp_code[L_SUBFR * HIBND_ACB_L_FAC];
    float tmp_code_preInt[L_SUBFR];
    float tmp = 1.0f;

    *voice_factors = VF_0th_PARAM + VF_1st_PARAM * voice_fac + VF_2nd_PARAM * voice_fac * voice_fac;

    if( (coder_type == VOICED || T0 > 115.5f) && core_brate > ACELP_8k00 )
    {
        tmp = 1.0f;
        *voice_factors *= tmp;
    }

    *voice_factors = min( max(0.000001f, *voice_factors), 0.999999f);

    if( L_frame == L_FRAME )
    {
        interp_code_5over2( code, tmp_code, L_SUBFR );

        for( i = 0; i < L_SUBFR * HIBND_ACB_L_FAC;  i++ )
        {
            bwe_exc[i + i_subfr * HIBND_ACB_L_FAC] = gain_pit * bwe_exc[i + i_subfr * HIBND_ACB_L_FAC] +
                    gain_code * tmp_code[i];
        }
    }
    else
    {
        for( i = 0; i < L_SUBFR;  i++ )
        {
            tmp_code_preInt[i] = gain_code * code[i] + 2 * gain_preQ * code_preQ[i];
        }

        interp_code_4over2( tmp_code_preInt, tmp_code, L_SUBFR );

        for( i = 0; i < L_SUBFR * 2;  i++ )
        {
            bwe_exc[i + i_subfr * 2] = gain_pit * bwe_exc[i + i_subfr * 2] + tmp_code[i];
        }
    }

    return;
}


/*-------------------------------------------------------------------*
 * get_tbe_bits()                                                    *
 *                                                                   *
 * Determine TBE bit consumption per frame from bit rate             *
 *-------------------------------------------------------------------*/

short get_tbe_bits(
    short bitrate,
    short bandwidth,
    short rf_mode
)
{
    short i, bits = 0;

    if( rf_mode )
    {
        /* TBE bits for core, primary frame */
        if( bandwidth == WB && bitrate == ACELP_13k20 )
        {
            /* Gain frame: 4, Gain shapes: 0, and LSFs: 2 */
            bits = NUM_BITS_SHB_FrameGain_LBR_WB + NUM_BITS_LBR_WB_LSF;
        }
        else if( bandwidth == SWB && bitrate == ACELP_13k20 )
        {
            /* Gain frame: 5, Gain shapes: 5, and lowrate LSFs: 8 */
            bits = NUM_BITS_SHB_FRAMEGAIN + NUM_BITS_SHB_SUBGAINS + 8;
        }
    }
    else
    {
        if( bandwidth == WB && bitrate == ACELP_9k60 )
        {
            bits = NUM_BITS_LBR_WB_LSF + NUM_BITS_SHB_FrameGain_LBR_WB;
        }
        else if( bandwidth == SWB || bandwidth == FB )
        {
            if( bitrate == ACELP_9k60 )
            {
                bits = NUM_BITS_SHB_FRAMEGAIN + NUM_BITS_SHB_SUBGAINS + 8;
            }
            else if( bitrate >= ACELP_13k20 && bitrate <= ACELP_32k )
            {
                bits = NUM_BITS_SHB_SUBGAINS + NUM_BITS_SHB_FRAMEGAIN + NUM_LSF_GRID_BITS + MIRROR_POINT_BITS;

                for( i=0; i<NUM_Q_LSF; i++ )
                {
                    bits += lsf_q_num_bits[i];
                }
            }

            if( bitrate >= ACELP_24k40 )
            {
                bits += NUM_BITS_SHB_ENER_SF + NUM_BITS_SHB_VF + NUM_BITS_SHB_RES_GS*NB_SUBFR16k;
            }

            if( bandwidth == SWB && (bitrate == ACELP_16k40 || bitrate == ACELP_24k40) )
            {
                bits += BITS_TEC + BITS_TFA;
            }

            if( bandwidth == FB )
            {
                /* fullband slope */
                bits += 4;
            }
        }
    }

    return bits;
}
