/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "prot.h"


/*-------------------------------------------------------------------*
 * ltd_stable()
 *
 *
 *-------------------------------------------------------------------*/

void ltd_stable(
    float frames_power[],       /* i/o: energy of several frames    */
    float ltd_stable_rate[],    /* o  : time-domain stable rate     */
    float frame_energy,         /* i  : current frame energy        */
    int frameloop               /* i  : number of  frames           */
)
{
    int   i;
    float tmp;
    float mid_frame_amp[28];
    float seg_amp;
    float dif,apow;

    frames_power[0] = (float)(sqrt(frame_energy)+0.001);

    if(frameloop<3)
    {
        for(i=1; i<POWER_NUM; i++)
        {
            frames_power[i] = 	frames_power[0];
        }
    }

    for(i=0; i<20; i++)
    {
        mid_frame_amp[i] = frames_power[2*i] + frames_power[2*i+1];
    }
    seg_amp = 0;
    for(i=0; i<20; i++)
    {
        seg_amp += mid_frame_amp[i];
    }
    seg_amp = seg_amp/20;
    dif  = 0;
    apow = 0;
    for(i=0; i<20; i++)
    {
        tmp = mid_frame_amp[i] - seg_amp;
        dif += tmp*tmp;
        apow += mid_frame_amp[i]*mid_frame_amp[i];
    }
    tmp = dif/(apow+0.0001f);
    ltd_stable_rate[0] = dif/(apow+FLT_MIN);


    seg_amp = 0;
    for(i=0; i<14; i++)
    {
        seg_amp += mid_frame_amp[i];
    }
    seg_amp = seg_amp/14;
    dif  = 0;
    apow = 0;
    for(i=0; i<14; i++)
    {
        tmp = mid_frame_amp[i] - seg_amp;
        dif += tmp*tmp;
        apow += mid_frame_amp[i]*mid_frame_amp[i];
    }
    ltd_stable_rate[1] = dif/(apow+0.0001f);

    seg_amp = 0;
    for(i=0; i<8; i++)
    {
        seg_amp += mid_frame_amp[i];
    }
    seg_amp = seg_amp/8;
    dif  = 0;
    apow = 0;
    for(i=0; i<8; i++)
    {
        tmp = mid_frame_amp[i] - seg_amp;
        dif += tmp*tmp;
        apow += mid_frame_amp[i]*mid_frame_amp[i];
    }
    ltd_stable_rate[2] = (float)(dif/(apow+0.0001f));
    ltd_stable_rate[3] = 0.90f*ltd_stable_rate[3] + 0.1f*ltd_stable_rate[2];

    for(i=POWER_NUM-1; i>0; i--)
    {
        frames_power[i] = frames_power[i-1];
    }

    return;
}
