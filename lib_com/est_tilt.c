/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "prot.h"


/*-------------------------------------------------------------------*
 * est_tilt()
 *
 * Estimate spectral tilt based on the relative E of adaptive
 * and innovative excitations
 *-------------------------------------------------------------------*/

float est_tilt(             /* o  : tilt of the code             */
    const float *adpt_exc,  /* i  : adaptive excitation vector   */
    const float gain_pit,   /* i  : adaptive gain                */
    const float *fixe_exc,  /* i  : algebraic exctitation vector */
    const float gain_code,  /* i  : algebraic code gain          */
    float *voice_fac, /* o  : voicing factor               */
    const short L_subfr,    /* i : subframe size                 */
    const short flag_tilt   /* i : flag for special tilt         */
)
{
    float ener, tmp, tilt_code;

    ener = dotp( adpt_exc, adpt_exc, L_subfr );
    ener *= gain_pit * gain_pit;                    /* energy of pitch excitation */

    tmp = dotp( fixe_exc, fixe_exc, L_subfr );
    tmp *= gain_code * gain_code;                   /* energy of innovative code excitation */

    /* find voice factor (1=voiced, -1=unvoiced) */
    *voice_fac = (float)((ener - tmp) / (ener + tmp + 0.01f));

    /* find tilt of code for next subframe */
    if (flag_tilt==0)
    {
        /*Between 0 (=unvoiced) and 0.5 (=voiced)*/
        tilt_code = (float)(0.25f*(1.0f + *voice_fac));
    }
    else if (flag_tilt==1)
    {
        /*Between 0.25 (=unvoiced) and 0.5 (=voiced)*/
        tilt_code = (float)(0.25f + (*voice_fac+1.0f)*0.125f);
    }
    else
    {
        /*Between 0.28 (=unvoiced) and 0.56 (=voiced)*/
        tilt_code = (float)(0.28f + (*voice_fac+1.0f)*0.14f);
    }

    return tilt_code;
}

