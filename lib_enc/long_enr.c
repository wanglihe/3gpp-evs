/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "prot.h"


/*-------------------------------------------------------------------*
 * long_enr()
 *
 * Compute relative energy, long-term average total noise energy and total active speech energy
 *-------------------------------------------------------------------*/

void long_enr(
    Encoder_State *st,                  /* i/o: encoder state structure                  */
    const float Etot,                 /* i  : total channel energy                     */
    const short localVAD_HE_SAD,      /* i  : HE-SAD flag without hangover             */
    short high_lpn_flag
)
{
    float tmp;
    /*-----------------------------------------------------------------*
     * Compute long term estimate of total noise energy
     * and total active speech energy
     *-----------------------------------------------------------------*/

    if( st->ini_frame < 4 )
    {
        st->lp_noise = st->totalNoise;
        tmp = st->lp_noise + 10.0f;

        if( st->lp_speech < tmp )
        {
            st->lp_speech = tmp;
        }
    }
    else
    {
        if ( st->ini_frame < 150 )
        {
            st->lp_noise = 0.95f * st->lp_noise + 0.05f * st->totalNoise;
        }
        else
        {
            st->lp_noise = 0.98f * st->lp_noise + 0.02f * st->totalNoise;
        }

        if ( localVAD_HE_SAD && !high_lpn_flag )
        {
            if ( ( st->lp_speech - Etot ) < 10.0f )
            {
                st->lp_speech = 0.98f * st->lp_speech + 0.02f * Etot;
            }
            else
            {
                st->lp_speech = st->lp_speech - 0.05f;
            }
        }
    }

    /*-----------------------------------------------------------------*
     * Initialize parameters for energy tracking and signal dynamics
     *-----------------------------------------------------------------*/
    return;
}
