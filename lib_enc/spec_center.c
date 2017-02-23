/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "prot.h"
#include "rom_enc.h"
#include "cnst.h"


/*-------------------------------------------------------------------*
 * spec_center()
 *
 *
 *-------------------------------------------------------------------*/

void spec_center(
    float spec_power[],               /*(i) energy of sub-band divided uniformly*/
    float sp_center[],                /*(o) spectral centroid*/
    int bandwidth                     /*(i) band width*/
)
{
    int i;
    float t_sp_center,frame_power;

    short sp_center_band_num = SP_CENTER_BAND_NUM_TAB[bandwidth];
    float sp_center_mem;

    t_sp_center = 0;
    frame_power = 0;

    for(i=1; i<sp_center_band_num; i++)
    {
        t_sp_center += spec_power[i]*(i);
        frame_power += 	spec_power[i];
    }
    sp_center[3] = (float)( (t_sp_center+DELTA1[bandwidth])/(frame_power+DELTA2[bandwidth]));

    sp_center_mem = 0;
    frame_power = 0;

    for(i=0; i<10; i++)
    {
        sp_center_mem += spec_power[i]*(i+1);
        frame_power += 	spec_power[i];
    }

    /* 107374184.f = 32768.f * 32768.f * 0.1 */
    t_sp_center = (float)( (sp_center_mem+107374184.f)/(frame_power+107374184.f));

    sp_center[0]= 0.7f*sp_center[0]+0.3f*t_sp_center;
    sp_center[2]= t_sp_center;

    if(bandwidth == CLDFBVAD_NB_ID)
    {
        t_sp_center = (float)( (sp_center_mem)/(frame_power+FLT_MIN));
        sp_center[0]= 0.7f*sp_center[0]+0.3f*t_sp_center;
        sp_center[2]= t_sp_center;
    }

    return;
}
