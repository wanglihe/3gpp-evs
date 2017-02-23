/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "prot.h"


/*-------------------------------------------------------------------*
 * cb_shape()
 *
 * pre-emphasis, pitch sharpening and formant sharpening of the algebraic codebook
 *-------------------------------------------------------------------*/

void cb_shape(
    const short preemphFlag,          /* i  : flag for pre-emphasis                           */
    const short pitchFlag,            /* i  : flag for pitch sharpening                       */
    const short scramblingFlag,       /* i  : flag for phase scrambling                       */
    const short formantFlag,          /* i  : flag for formant sharpening                     */
    const short formantTiltFlag,      /* i  : flag for formant tilt                           */
    const float g1,                   /* i  : formant sharpening numerator weighting          */
    const float g2,                   /* i  : formant sharpening denominator weighting        */
    const float *p_Aq,                /* i  : LP filter coefficients                          */
    float *code,                /* i/o: signal to shape                                 */
    const float tilt_code,            /* i  : tilt of code                                    */
    const float pt_pitch              /* i  : pointer to current subframe fractional pitch    */
)
{
    float buff[M+L_SUBFR], A_num[M+1], A_den[M+1];
    float tmp, tilt;
    short i, round_T0;

    /* pre-emphasize the algebraic codebook */
    if ( preemphFlag )
    {
        tmp = 0.0f;
        preemph( code, tilt_code, L_SUBFR, &tmp );
    }

    /* pitch sharpening */
    if ( pitchFlag )
    {
        round_T0 = (short)(pt_pitch + 0.4f);
        for (i = round_T0; i < L_SUBFR; i++)
        {
            code[i] += code[i - round_T0] * PIT_SHARP;
        }
    }

    /* phase scrambling filter */
    if ( scramblingFlag )
    {
        buff[0] = code[0];
        for (i = 1; i < L_SUBFR; i++)
        {
            buff[i] = code[i];
            code[i] = 0.7f * buff[i] + buff[i-1] - 0.7f * code[i-1];
        }
    }

    /* formant sharpening (only on active speech)  */
    if ( formantFlag || formantTiltFlag )
    {
        weight_a( p_Aq, A_num, g1, M );
        weight_a( p_Aq, A_den, g2, M );

        set_f( buff, 0, M+L_SUBFR );

        /* formant tilt */
        if( formantTiltFlag )
        {
            mvr2r( A_num, buff+M, M+1 );
            syn_filt( A_den, M, buff+M, buff+M, L_SUBFR, buff, 0 );
            tilt = get_gain( buff+M+1, buff+M, L_SUBFR-1, NULL );
            tmp = 0.0;
            preemph( code, 0.5f*tilt_code-0.25f*tilt, L_SUBFR, &tmp );
        }
        else
        {
            mvr2r( code, buff+M, L_SUBFR );
            residu( A_num, M, buff+M, code, L_SUBFR );
            syn_filt( A_den, M, code, code, L_SUBFR, buff, 0 );
        }
    }

    return;
}
