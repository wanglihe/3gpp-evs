/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"


/*-------------------------------------------------------------------*
 * Local constants
 *-------------------------------------------------------------------*/

#define ATT_LENGHT            64
#define ATT_SEG_LEN           (L_FRAME/ATT_LENGHT)
#define INV_ATT_SEG_LEN       (1.0f/ATT_SEG_LEN)
#define INV_L_FRAME           (1.0f/L_FRAME)


/*-------------------------------------------------------------------*
 * pre_echo_att()
 *
 * Attenuation of the pre-echo when encoder specifies an attack
 *-------------------------------------------------------------------*/

void pre_echo_att(
    float *Last_frame_ener,     /* i/o: Energy of the last frame         */
    float *exc,                 /* i/o: Excitation of the current frame  */
    const short attack_flag,          /* i  : flag signalling attack encoded by AC mode (GSC) */
    const short last_coder_type       /* i  : Last coder type                  */
)
{
    float etmp;
    float etmp1;
    float finc[ATT_LENGHT], ratio;
    short attack_pos, i;

    if ( attack_flag == 1 && last_coder_type == AUDIO)
    {
        /*-------------------------------------------------------------------------*
         * Find where the onset (attack) occurs by computing the energy per section
         * The inverse weighting aims to favor the first maxima in case of
         * gradual onset
         *-------------------------------------------------------------------------*/

        for(i = 0; i < ATT_LENGHT; i++)
        {
            finc[i] = sum2_f( exc + i*ATT_SEG_LEN, ATT_SEG_LEN )*((float)(ATT_LENGHT-i)/(ATT_LENGHT));
        }

        etmp = -1;
        attack_pos = maximum(finc, ATT_LENGHT, &etmp);

        /* Scaled the maximum energy and allowed 6 dB increase*/
        etmp *= INV_ATT_SEG_LEN;
        etmp1 = etmp;
        *Last_frame_ener *= 4.0f;

        /* If the maximum normalized energy > last frame energy + 6dB */
        if( etmp > *Last_frame_ener && attack_pos > 0 )
        {
            /* Find the average energy before the attack */
            etmp = sum_f( finc, attack_pos ) +  0.01f;
            etmp /= (attack_pos*ATT_SEG_LEN);

            /* Find the correction factor and apply it before the attack */
            ratio = (float)sqrt(*Last_frame_ener/etmp);


            /* Pre-echo atttenuation should never increase the energy */
            ratio = min(ratio, 1.0f);
            for(i = 0; i < attack_pos*ATT_SEG_LEN; i++)
            {
                exc[i] *= ratio;
            }
        }
        *Last_frame_ener = etmp1;
    }
    else
    {
        /*-------------------------------------------------------*
         * In normal cases, just compute the energy of the frame
         *-------------------------------------------------------*/

        etmp = sum2_f( exc, L_FRAME ) + 0.01f;

        etmp *= INV_L_FRAME;
        *Last_frame_ener = etmp;
    }

    return;
}
