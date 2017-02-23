/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "prot.h"

/*-------------------------------------------------------------------*
 * frame_spec_dif_cor_rate()
 *
 *
 *-------------------------------------------------------------------*/

void frame_spec_dif_cor_rate(
    float spec_amp[],               /*(i) spectral amplitude*/
    float pre_spec_low_dif[],       /*(io) low spectrum different*/
    float f_tonality_rate[]         /*(o) tonality rate*/
)
{
    int i;
    float spec_low_dif[59];
    float tmp,spec_low_dif_cor_rate,spec_dif_cor_rate;
    float m,dx,dy;


    for(i=0; i<PRE_SPEC_DIF_NUM; i++)
    {
        tmp = spec_amp[i+6] - spec_amp[i+5];
        if(tmp < 0)
        {
            spec_low_dif[i] = 0;
        }
        else
        {
            spec_low_dif[i] = tmp;
        }
    }
    m  = 0;
    dx = 0;
    dy = 0;
    for(i=0; i<PRE_SPEC_DIF_NUM; i++)
    {
        m  += spec_low_dif[i]*pre_spec_low_dif[i];
        dx += spec_low_dif[i]*spec_low_dif[i];
        dy += pre_spec_low_dif[i]*pre_spec_low_dif[i];
    }

    /* 1073741.824 = 0.001 * 32768 * 32768 */
    spec_low_dif_cor_rate = (float)(m/sqrt((dx*dy+1073741.824f)));

    spec_dif_cor_rate  = spec_low_dif_cor_rate;
    f_tonality_rate[0] = spec_dif_cor_rate;
    f_tonality_rate[1] = f_tonality_rate[1]*0.96f + spec_dif_cor_rate*0.04f;
    f_tonality_rate[2] = f_tonality_rate[2]*0.90f + spec_dif_cor_rate*0.1f;

    for(i=0; i<PRE_SPEC_DIF_NUM; i++)
    {
        pre_spec_low_dif[i] = spec_low_dif[i];
    }

    return;
}
