/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "prot.h"

/*-------------------------------------------------------------------*
 * Local constants
 *-------------------------------------------------------------------*/

#define GAIN_PIT_MAX      1.2f
#define ACELP_GAINS_CONST 0.8f    /* adaptive codebook gain constraint */


/*-------------------------------------------------------------------*
 * Local functions
 *-------------------------------------------------------------------*/

static float adpt_enr( const short codec_mode, const float *exc, const float *h1, float *y1, const short L_subfr, float *gain, float *g_corr,
                       const short clip_gain, const float *xn, float *xn2, short use_prev_sf_pit_gain );

/*-------------------------------------------------------------------*
 * lp_filt_exc_enc()
 *
 * Low-pass filtering of the adaptive excitation
 * Innovation target construction
 * Gain quantization limitation
 *-------------------------------------------------------------------*/

short lp_filt_exc_enc(
    const short codec_mode,     /* i  : codec mode                                  */
    const long  core_brate,     /* i  : core bitrate                                */
    const short Opt_AMR_WB,     /* i  : flag indicating AMR-WB IO mode              */
    const short coder_type,     /* i  : coding type                                 */
    const short i_subfr,        /* i  : subframe index                              */
    float *exc,           /* i/o: pointer to excitation signal frame          */
    const float *h1,            /* i  : weighted filter input response              */
    const float *xn,            /* i  : target vector                               */
    float *y1,            /* o  : zero-memory filtered adaptive excitation    */
    float *xn2,           /* o  : target vector for innovation search         */
    const short L_subfr,        /* i  : length of vectors for gain quantization     */
    const short L_frame,        /* i  : frame size                                  */
    float *g_corr,        /* o  : ACELP correlation values                    */
    const short clip_gain,      /* i  : adaptive gain clipping flag                 */
    float *gain_pit,      /* o  : adaptive excitation gain                    */
    short *lp_flag        /* i/o  : mode selection                            */
)
{
    float ener, ener_tmp, gain1, gain2, g_corr2[2], exc_tmp[L_FRAME16k], xn2_tmp[L_FRAME16k];
    float y1_tmp[L_FRAME16k];
    short select, i;
    short use_prev_sf_pit_gain;

    gain1 = 0.0f;
    gain2 = 0.0f;
    ener = 0.0f;
    ener_tmp = 0.0f;
    use_prev_sf_pit_gain = 0;

    /*-----------------------------------------------------------------*
     * Select LP filtering flag
     *-----------------------------------------------------------------*/

    if( codec_mode == MODE1 )
    {
        if ( ( Opt_AMR_WB || coder_type == GENERIC || coder_type == TRANSITION ) && core_brate < ACELP_11k60 )
        {
            *lp_flag = LOW_PASS;
        }
        else if ( core_brate >= ACELP_11k60 && coder_type != AUDIO )
        {
            *lp_flag = NORMAL_OPERATION;
        }
        else
        {
            *lp_flag = FULL_BAND;
        }
    }

    /*----------------------------------------------------------------*
     * Find energy of the fixed cb. target with respect to adaptive
     *   exc. filtering
     * Find flag about splitting the gain quantizer in half
     * - find the target energy if adaptive exc. is not filtered
     * - filter the adaptive excitation and find the target energy
     *   if the exc. is filtered.
     *----------------------------------------------------------------*/

    if( codec_mode == MODE2 && coder_type == 100 )
    {
        use_prev_sf_pit_gain = 1;
    }

    if( *lp_flag == FULL_BAND || *lp_flag == NORMAL_OPERATION )
    {
        if( use_prev_sf_pit_gain == 1 )
        {
            ener = adpt_enr( codec_mode, &exc[i_subfr], h1, y1, L_subfr, gain_pit, g_corr, clip_gain, xn, xn2, use_prev_sf_pit_gain );
        }
        else
        {
            ener = adpt_enr( codec_mode, &exc[i_subfr], h1, y1, L_subfr, &gain1, g_corr, clip_gain, xn, xn2, use_prev_sf_pit_gain );
        }
    }

    /*----------------------------------------------------------------*
     * Find energy of the fixed cb. target with respect to adaptive
     *   exc. filtering
     * filter the adaptive excitation and find the target energy
     *   if the exc. is filtered.
     *----------------------------------------------------------------*/

    if( *lp_flag == LOW_PASS || *lp_flag == NORMAL_OPERATION )
    {
        if( codec_mode == MODE2 && L_frame == L_FRAME16k )
        {
            for (i=0; i<L_subfr; i++)
            {
                exc_tmp[i] = (float)(0.21f * exc[i-1+i_subfr] + 0.58f * exc[i+i_subfr] + 0.21f * exc[i+1+i_subfr]);
            }
        }
        else
        {
            for ( i=0; i<L_subfr; i++ )
            {
                exc_tmp[i] = (float)(0.18f * exc[i-1+i_subfr] + 0.64f * exc[i+i_subfr] + 0.18f * exc[i+1+i_subfr]);
            }
        }

        if( use_prev_sf_pit_gain == 1 )
        {
            ener_tmp = adpt_enr( codec_mode, exc_tmp, h1, y1_tmp, L_subfr, gain_pit, g_corr2, clip_gain, xn, xn2_tmp, use_prev_sf_pit_gain );
        }
        else
        {
            ener_tmp = adpt_enr( codec_mode, exc_tmp, h1, y1_tmp, L_subfr, &gain2, g_corr2, clip_gain, xn, xn2_tmp, use_prev_sf_pit_gain );
        }
    }

    /*-----------------------------------------------------------------*
     * use the best prediction (minimize quadratic error)
     *-----------------------------------------------------------------*/

    if( ( (ener_tmp < ener) && (*lp_flag == NORMAL_OPERATION) ) || (*lp_flag == LOW_PASS) )
    {
        /* use the LP filter for pitch excitation prediction */
        select = LOW_PASS;
        mvr2r(exc_tmp, &exc[i_subfr], L_subfr);
        mvr2r(y1_tmp, y1, L_subfr);
        mvr2r( xn2_tmp, xn2, L_subfr );

        if( use_prev_sf_pit_gain == 0 )
        {
            *gain_pit = gain2;
            g_corr[0] = g_corr2[0];
            g_corr[1] = g_corr2[1];
        }
    }
    else
    {
        /* no LP filter used for pitch excitation prediction */
        select = FULL_BAND;
        if( use_prev_sf_pit_gain == 0 )
        {
            *gain_pit = gain1;
        }
    }

    return(select);
}

/*-------------------------------------------------------------------*
 * adpt_enr()
 *
 * Find the adaptive excitation energy
 * This serves to decide about the filtering of the adaptive excitation
 *-------------------------------------------------------------------*/

static float adpt_enr(
    const short codec_mode,  /* i  : codec mode                             */
    const float *exc,        /* i  : Excitation vector                      */
    const float *h1,         /* i  : impuls response                        */
    float *y1,         /* o  : zero-memory filtered adpt. excitation  */
    const short L_subfr,     /* i  : vector length                          */
    float *gain,       /* o  : subframe adaptive gain                 */
    float *g_corr,     /* o  : correlations for adptive gain          */
    const short clip_gain,   /* i  : adaptive gain clipping flag            */
    const float *xn,         /* i  : adaptive codebook target               */
    float *xn2,        /* o  : algebraic codebook target              */
    short use_prev_sf_pit_gain /* i : flag to use prev sf pitch gain or not */
)
{
    float ener;

    conv( exc, h1, y1, L_subfr );
    if( use_prev_sf_pit_gain == 0 )
    {
        *gain = corr_xy1( xn, y1, g_corr, L_subfr, codec_mode==MODE2 );

        /* clip gain, if necessary to avoid problems at decoder */
        if ( clip_gain == 1 && *gain > 0.95f )
        {
            *gain = 0.95f;
        }

        if( clip_gain == 2 && *gain > 0.65f )
        {
            *gain = 0.65f;
        }
    }

    /* find energy of new target xn2[] */
    updt_tar( xn, xn2, y1, *gain, L_subfr );
    ener = dotp( xn2, xn2, L_subfr );

    return ener;
}

/*-------------------------------------------------------------------*
 * corr_xy1()
 *
 * Find the correlations between the target xn[] and the filtered adaptive
 * codebook excitation y1[]. ( <y1,y1>  and -2<xn,y1> )
 *-------------------------------------------------------------------*/

float corr_xy1(            /* o  : pitch gain  (0..GAIN_PIT_MAX)         */
    const float xn[],      /* i  : target signal                         */
    const float y1[],      /* i  : filtered adaptive codebook excitation */
    float g_corr[],  /* o  : correlations <y1,y1>  and -2<xn,y1>   */
    const short L_subfr,   /* i  : vector length                          */
    const short norm_flag  /* i : flag for constraining pitch contribution */
)
{
    float temp1, temp2, gain,gain_p_snr;

    /*-----------------------------------------------------------------*
     * Find the ACELP correlations and the pitch gain
     * (for current subframe)
     *-----------------------------------------------------------------*/

    /* Compute scalar product <xn[],y1[]> */
    temp1 = dotp(xn, y1, L_subfr);

    /* Compute scalar product <y1[],y1[]> */
    temp2 = dotp(y1, y1, L_subfr) + 0.01f;

    g_corr[0] = temp2;
    g_corr[1] = -2.0f*temp1 + 0.01f;

    /* find pitch gain and bound it by [0,GAIN_PIT_MAX] */
    if( norm_flag )
    {
        gain = (temp1 + 0.01f)/temp2;
    }
    else
    {
        gain = temp1/temp2;
    }

    if( gain < 0.0f )
    {
        gain = 0.0f;
    }
    if( gain > GAIN_PIT_MAX )
    {
        gain = GAIN_PIT_MAX;
    }

    /*Limit the energy of pitch contribution*/
    if( norm_flag )
    {
        /* Compute scalar product <xn[],xn[]> */
        temp1 = dotp(xn, xn, L_subfr);
        gain_p_snr = ACELP_GAINS_CONST*sqrt(temp1/temp2);

        if( gain>gain_p_snr )
        {
            gain = gain_p_snr;
        }
    }

    return gain;
}
