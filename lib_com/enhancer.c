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

static void agc2( const float *sig_in, float *sig_out,const short l_trm );

/*---------------------------------------------------------------------*
 * enhancer()
 *
 * Enhancement of the excitation signal before synthesis
 *---------------------------------------------------------------------*/

void enhancer(
    const short codec_mode,       /* i  : flag indicating Codec Mode              */
    const long  core_brate,       /* i  : core bitrate                            */
    const short cbk_index,        /* i  :                                         */
    const short Opt_AMR_WB,       /* i  : flag indicating AMR-WB IO mode          */
    const short coder_type,       /* i  : coding type                             */
    const short L_frame,          /* i  : frame size                              */
    const float voice_fac,        /* i  : subframe voicing estimation             */
    const float stab_fac,         /* i  : LP filter stablility measure            */
    const float norm_gain_code,   /* i  : normalized innovative cb. gain          */
    const float gain_inov,        /* i  : gain of the unscaled innovation         */
    float *gc_threshold,    /* i/o: code threshold                          */
    float *code,            /* i/o: innovation                              */
    float *pt_exc2,         /* i/o: adapt. excitation/total exc.            */
    const float gain_pit,         /* i  : Quantized pitch gain                    */
    float *dispMem          /* i/o: Phase dispersion algorithm memory       */
)
{
    float tmp, gain_code, new_norm_gain_code, fac;
    short i;
    float pit_sharp;
    float excp[L_SUBFR];

    pit_sharp = gain_pit;

    /*-----------------------------------------------------------------*
     * Phase dispersion
     *
     * Enhance noise at low bit rates
     *-----------------------------------------------------------------*/

    i = 2;            /* no dispersion   */
    if( Opt_AMR_WB )
    {
        if ( core_brate <= ACELP_6k60 )
        {
            i = 0;        /* high dispersion  */
        }
        else if ( core_brate <= ACELP_8k85 )
        {
            i = 1;        /* low dispersion  */
        }
    }
    else if( codec_mode == MODE1 && coder_type != UNVOICED )
    {
        if ( core_brate <= ACELP_7k20 )
        {
            i = 0;        /* high dispersion  */
        }
        else if ( ( coder_type == GENERIC || coder_type == TRANSITION || coder_type == AUDIO || coder_type == INACTIVE ) && core_brate <= ACELP_9k60 )
        {
            i = 1;        /* low dispersion  */
        }
    }
    else if( codec_mode == MODE2 )
    {
        if( ((coder_type!=VOICED) && cbk_index<=2) || ((coder_type==UNVOICED) && L_frame==L_FRAME && cbk_index<=10) || ((coder_type==UNVOICED) && L_frame==L_FRAME16k && cbk_index<=14))
        {
            i = 0;        /* high dispersion  */
        }
        else if( (coder_type!=VOICED) && (cbk_index<=7) )
        {
            i = 1;        /* low dispersion  */
        }
    }

    phase_dispersion( norm_gain_code, gain_pit, code, i, dispMem );

    /*------------------------------------------------------------
     * Noise enhancer
     *
     * Enhance excitation on noise (modify code gain). If signal is noisy and LPC filter is stable,
     * move code gain 1.5 dB towards its threshold. This decreases by 3 dB noise energy variation.
     *-----------------------------------------------------------*/

    if ( norm_gain_code < *gc_threshold )
    {
        new_norm_gain_code = (float)(norm_gain_code * 1.19f);
        if ( new_norm_gain_code > *gc_threshold )
        {
            new_norm_gain_code = *gc_threshold;
        }
    }
    else
    {
        new_norm_gain_code = (float)(norm_gain_code / 1.19f);
        if ( new_norm_gain_code < *gc_threshold )
        {
            new_norm_gain_code = *gc_threshold;
        }
    }
    *gc_threshold = new_norm_gain_code;

    /* calculate new code gain */
    fac = stab_fac * (0.5f * (1.0f - voice_fac));  /* 1 = unvoiced, 0 = voiced */
    gain_code = fac * new_norm_gain_code + (1.0f - fac) * norm_gain_code;
    gain_code *= gain_inov;

    for (i=0; i<L_SUBFR; i++)
    {
        code[i] *= gain_code;
    }

    /*------------------------------------------------------------*
     * Pitch enhancer
     *
     * Enhance excitation on voiced (HP filtering of code). On voiced signal, filter code by a smooth HP
     * filter to decrease the energy of code at low frequency
     *------------------------------------------------------------*/

    if( !Opt_AMR_WB && codec_mode == MODE1 && coder_type == UNVOICED )
    {
        mvr2r( code, pt_exc2, L_SUBFR );
    }
    else
    {
        if ( Opt_AMR_WB && ( core_brate == ACELP_8k85 || core_brate == ACELP_6k60 ) )
        {
            pit_sharp = gain_pit;
            if( pit_sharp > 1.0 )
            {
                pit_sharp = 1.0;
            }

            if ( pit_sharp > 0.5 )
            {
                for (i = 0; i < L_SUBFR; i++)
                {
                    excp[i] = pt_exc2[i] * pit_sharp * 0.25f;
                }
            }
        }

        /*-----------------------------------------------------------------
         * Do a simple noncasual "sharpening": effectively an FIR
         * filter with coefs [-tmp 1.0 -tmp] where tmp = 0 ... 0.25
         * This is applied to code and added to exc2
         *-----------------------------------------------------------------*/
        if( L_frame == L_FRAME16k )
        {
            tmp = (float)(0.150f*(1.0f+voice_fac));     /* 0.30=voiced, 0=unvoiced */
        }
        else
        {
            tmp = (float)(0.125f * (1.0f + voice_fac));     /* 0.25=voiced, 0=unvoiced */
        }
        pt_exc2[0] += code[0] - (tmp * code[1]);
        for ( i=1; i<L_SUBFR-1; i++ )
        {
            pt_exc2[i] += code[i] - (tmp*code[i-1]) - (tmp*code[i+1]);
        }
        pt_exc2[L_SUBFR-1] += code[L_SUBFR-1] - (tmp*code[L_SUBFR-2]);

        if ( Opt_AMR_WB && ( core_brate == ACELP_8k85 || core_brate == ACELP_6k60 ) )
        {
            if ( pit_sharp > 0.5f )
            {
                for (i = 0; i < L_SUBFR; i++)
                {
                    excp[i] += pt_exc2[i];
                }

                agc2( pt_exc2, excp, L_SUBFR );
                mvr2r( excp, pt_exc2, L_SUBFR );
            }
        }
    }

    return;
}

/*-----------------------------------------------------------------------*
 * agc2()
 *
 * Adaptive gain control
 *-----------------------------------------------------------------------*/

static void agc2(
    const float *sig_in,   /* i  : postfilter input signal  */
    float *sig_out,  /* i/o: postfilter output signal */
    const short l_trm      /* i  : subframe size            */
)
{
    short i;
    float gain_in, gain_out;
    float g0, gain;


    gain_out = 0.0f;
    for(i=0; i<l_trm; i++)
    {
        gain_out += sig_out[i]*sig_out[i];
    }

    if (gain_out == 0.0f)
    {
        return;
    }

    gain_in = 0.0f;
    for(i=0; i<l_trm; i++)
    {
        gain_in += sig_in[i]*sig_in[i];
    }
    if (gain_in == 0.0f)
    {
        g0 = 0.0f;
    }
    else
    {
        g0 = (float) sqrt(gain_in / gain_out);
    }

    gain = g0;
    for (i=0; i<l_trm; i++)
    {
        sig_out[i] *= gain;
    }

    return;
}
