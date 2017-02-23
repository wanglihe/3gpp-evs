/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"


/*----------------------------------------------------------------------------------*
* fer_energy()
*
* Estimation of pitch-synchronous (voiced sounds) or half-frame energy
*----------------------------------------------------------------------------------*/

void fer_energy(
    const int L_frame,    /* i  : frame length                           */
    const short clas,     /* i  : frame classification                   */
    const float *synth,   /* i  : synthesized speech at Fs = 12k8 Hz     */
    const float pitch,    /* i  : pitch period                           */
    float *enr,     /* o  : pitch-synchronous or half_frame energy */
    const short offset    /* i  : speech pointer offset (0 or L_frame)   */
)
{
    short        len;
    const float  *pt_synth;

    if( clas == VOICED_CLAS || clas == ONSET || clas == SIN_ONSET ) /* Voiced or Onset current frame */
    {
        len = (short)( pitch + 0.5f );    /* pitch value */

        pt_synth = synth;
        if( offset != 0 )
        {
            pt_synth = synth + L_frame - len;
        }

        emaximum( pt_synth, len, enr );  /* pitch synchronous E */
    }
    else
    {
        pt_synth = synth;
        if( offset != 0 )
        {
            pt_synth = synth + L_frame/2;
        }

        *enr = dotp( pt_synth, pt_synth, L_frame/2 );
        *enr /= (float)(L_frame/2);

    }
    return;
}


/*------------------------------------------------------------------------*
 * frame_energy()
 *
 * Compute pitch-synchronous energy at the frame end
 *------------------------------------------------------------------------*/

float frame_energy(
    const short L_frame,           /* i: length of the frame                            */
    const float *pitch,            /* i: pitch values for each subframe                 */
    const float *speech,           /* i: pointer to speech signal for E computation     */
    const float lp_speech,         /* i: long-term active speech energy average         */
    float       *frame_ener        /* o: pitch-synchronous energy at frame end          */
)
{
    float enern;
    const float *pt1;
    short len;
    float dotProd;

    len = (short)( 0.5f * (pitch[2] + pitch[3]) + 0.5f );
    if( len < L_SUBFR )
    {
        len *= 2;
    }

    pt1 = speech + L_frame - len;

    dotProd = dotp( pt1, pt1, len );
    if ( 0 == dotProd )
    {
        *frame_ener = MIN_LOG_VAL_60dB;
    }
    else
    {
        *frame_ener = 10.0f * (float) log10( dotProd / (float) len );
    }
    enern = *frame_ener - lp_speech;

    return enern;
}
