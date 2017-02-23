/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "prot.h"


/*---------------------------------------------------------------------*
 * Local constants
 *---------------------------------------------------------------------*/

#define TILT_COMP_LIM       0.75f

/*---------------------------------------------------------*
 * Local functions
 *---------------------------------------------------------*/

static float calc_tilt(const float *x, const short len);

/*--------------------------------------------------------------------*
 * stat_noise_uv_mod()
 *
 * Modifies excitation signal in stationary noise segments
 *--------------------------------------------------------------------*/

void stat_noise_uv_mod(
    const short coder_type,       /* i  : coder type                           */
    float noisiness,        /* i  : noisiness parameter                  */
    const float *lsp_old,         /* i  : old LSP vector at 4th sfr            */
    const float *lsp_new,         /* i  : LSP vector at 4th sfr                */
    const float *lsp_mid,         /* i  : LSP vector at 2nd sfr                */
    float *Aq,              /* o  : A(z) quantized for the 4 subframes   */
    float *exc2,            /* o  : excitation buffer                    */
    const short bfi,              /* i  : bad frame indicator                  */
    float *ge_sm,           /* i/o: smoothed excitation gain             */
    short *uv_count,        /* i/o: unvoiced counter                     */
    short *act_count,       /* i/o: activation counter                   */
    float lspold_s[],       /* i/o: old LSP                              */
    short *noimix_seed,     /* i/o: mixture seed                         */
    float *st_min_alpha,    /* i/o: minimum alpha                        */
    float *exc_pe,          /* i/o: memory of the preemphasis filter     */
    const long  bitrate,          /* i  : core bitrate                         */
    const short bwidth            /* i  : Sampling rate                        */
)
{
    short i, k;
    short i_subfr;
    float exctilt;
    float vare, ge, randval;
    float alpha, min_alpha;
    float lspnew_s[M], oldlsp_mix[M], midlsp_mix[M], newlsp_mix[M];
    float beta;
    float noimix_fac;

    alpha = 1.0f;
    min_alpha = 0.5f;

    /*---------------------------------------------------------*
     * Update minimum mixing factor alpha
     *---------------------------------------------------------*/

    /* Activate modifications for WB/SWB <= 9.6kbps for NB only at 9.6kbps */
    if ( coder_type == INACTIVE && ( bitrate == ACELP_9k60 || (bitrate < ACELP_9k60 && bwidth > NB) ) )
    {
        if ( !bfi )
        {
            min_alpha = max( noisiness/31.0f * 0.5f + 0.5f, *st_min_alpha - 0.05f );
            *st_min_alpha = min_alpha;
        }
        else
        {
            min_alpha = *st_min_alpha;
        }
    }

    /*---------------------------------------------------------*
     * Mix excitation signal with random noise
     *---------------------------------------------------------*/

    /* Activate modifications for WB/SWB <= 9.6kbps for NB only at 9.6kbps */
    if ( coder_type == INACTIVE && (bitrate == ACELP_9k60 || (bitrate < ACELP_9k60 && bwidth > NB ) ) )
    {
        /* preemphasize the excitation signal with its tilt */
        if( min_alpha < TILT_COMP_LIM )
        {
            for (i_subfr=0; i_subfr<L_FRAME; i_subfr+=L_SUBFR)
            {
                exctilt = calc_tilt( &exc2[i_subfr], L_SUBFR );
                exctilt = (TILT_COMP_LIM - min_alpha)/(TILT_COMP_LIM - 0.5f) * exctilt;

                preemph( &exc2[i_subfr], exctilt, L_SUBFR, exc_pe );
            }
        }

        /* set the mixing factor alpha */
        (*uv_count)++;
        if (*uv_count <= START_NG)
        {
            *act_count = 3;
            alpha = 1;
            mvr2r(lsp_new, lspold_s, M);

        }
        else
        {
            *act_count = 0;
            if (*uv_count > FULL_NG)
            {
                *uv_count = FULL_NG;
            }

            alpha = 1 + ((float)*uv_count - START_NG)/((float)FULL_NG - START_NG) * (min_alpha - 1.0f);
        }

        /* calculate lowpass-filtered excitation gain */
        vare = 0.01f;
        for (i=0; i<L_FRAME; i++)
        {
            vare += exc2[i] * exc2[i];
        }

        ge = (float)sqrt(vare/(float)L_FRAME);
        if(*uv_count == 1)
        {
            *ge_sm = ge;
        }
        else
        {
            *ge_sm = ISP_SMOOTHING_QUANT_A1 **ge_sm + (1-ISP_SMOOTHING_QUANT_A1) * ge;
        }

        /* generate mixture of excitation and noise */
        beta = 2*(alpha -0.5f);
        noimix_fac = (beta + *ge_sm / ge * (1 - beta))/(float)sqrt(alpha*alpha + (1-alpha)*(1-alpha));

        for (i=0; i<L_FRAME; i++)
        {
            randval = ge * (float)sqrt(12.0f) * ((float)own_random(noimix_seed)/65536.0f);
            exc2[i] = noimix_fac * ( exc2[i] * alpha + randval * (1-alpha) );
        }

        /* generate low-pass filtered version of LSP coefficients */
        for ( k=0; k<M; k++ )
        {
            lspnew_s[k] = (float) ISP_SMOOTHING_QUANT_A1 * lspold_s[k] + (float) (1-ISP_SMOOTHING_QUANT_A1) * lsp_new[k];
        }

        for ( i=0; i<M; i++ )
        {
            oldlsp_mix[i] = beta * lsp_old[i] + (1-beta) * lspold_s[i];
            midlsp_mix[i] = beta * lsp_mid[i] + (1-beta) * 0.5f * (lspold_s[i] + lspnew_s[i]);
            newlsp_mix[i] = beta * lsp_new[i] + (1-beta) * lspnew_s[i];
        }

        /* redo the interpolation of LSP coefficients and recalculte A(z) */
        int_lsp4( L_FRAME, oldlsp_mix, midlsp_mix, newlsp_mix, Aq, M, 0 );

        mvr2r( lspnew_s, lspold_s, M );
    }
    else
    {
        /* active signal - reset counters */
        (*act_count)++;
        if (*act_count > 3)
        {
            *act_count = 3;
            *uv_count = 0;
        }
    }

    return;
}

/*---------------------------------------------------------------------------*
 * calc_tilt()
 *
 * Calculate spectral tilt by means of 1st-order LP analysis
 *---------------------------------------------------------------------------*/

static float calc_tilt(
    const float *x,            /* i  : Signal input   */
    const short len            /* i  : lenght         */
)
{
    int i;
    float r0, r1;

    r0 = 0;
    r1 = 0;

    for ( i=0; i<len-1; i++ )
    {
        r0 += x[i]*x[i];
        r1 += x[i]*x[i+1];
    }

    return r1/r0;
}
