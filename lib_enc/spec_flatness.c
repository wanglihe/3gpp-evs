/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "prot.h"

#include <assert.h>


/*-------------------------------------------------------------------*
 * spec_flatness()
 *
 *
 *-------------------------------------------------------------------*/

void spec_flatness(
    float spec_amp[],               /*(i) spectral amplitude*/
    float smooth_spec_amp[],        /*(i) smoothed spectral amplitude*/
    float sSFM[]                    /*(o) spectral flatness rate*/
)
{
    int i;
    double prods,sums;

    float SFM;
    for(i=MIN_AMP_ID; i<=MAX_AMP_ID; i++)
    {
        smooth_spec_amp[i-MIN_AMP_ID] = smooth_spec_amp[i-MIN_AMP_ID]*0.7f + spec_amp[i]*0.3f;
    }
    /* sSFM1 */
    prods = 1;
    sums = 0;
    assert(MIN_AMP_ID<=5);
    for(i=(5-MIN_AMP_ID); i<(20-MIN_AMP_ID); i++)
    {
        prods = smooth_spec_amp[i]*prods;
        sums  = sums + smooth_spec_amp[i];
    }

    if(prods>0)
    {
        prods = pow(prods,1.0/15);
    }
    else
    {
        prods=0;
    }
    sums = sums/15;

    SFM = (float)((prods+3276.8f)/(sums+3276.8f));
    sSFM[0] = 0.85f*sSFM[0] + 0.15f*SFM;

    /* sSFM2 */
    prods = 1;
    sums = 0;
    for(i=(20-MIN_AMP_ID); i<(40-MIN_AMP_ID); i++)
    {
        prods = smooth_spec_amp[i]*prods;
        sums  = sums + smooth_spec_amp[i];
    }

    if(prods>0)
    {
        prods = pow(prods,1.0/20);
    }
    else
    {
        prods =0;
    }
    sums = sums/20;
    SFM= (float)((prods+3276.8f)/(sums+3276.8f));
    sSFM[1] = 0.85f*sSFM[1] + 0.15f*SFM;
    /* sSFM3 */
    prods = 1;
    sums  = 0;
    for(i=(40-MIN_AMP_ID); i<=(MAX_AMP_ID-MIN_AMP_ID); i++)
    {
        prods = smooth_spec_amp[i]*prods;
        sums  = sums + smooth_spec_amp[i];
    }

    if(prods>0)
    {
        prods = pow(prods,0.04);
    }
    else
    {
        prods =0;
    }
    sums = sums/25;
    SFM = (float)((prods+3276.8f)/(sums+3276.8f));
    sSFM[2] = 0.85f*sSFM[2] + 0.15f*SFM;


    return;
}
