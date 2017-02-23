/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <stdlib.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"

/*-------------------------------------------------------------------*
 * vad_param_updt()
 *
 * Update parameters used by the VAD or DTX
 *--------------------------------------------------------------------*/

void vad_param_updt(
    Encoder_State *st,              /* i/o: encoder state structure                               */
    const short pitch[3],         /* i  : open loop pitch lag for each half-frame               */
    const float voicing[3],       /* i  : maximum normalized correlation for each half-frame    */
    const float corr_shift,       /* i  : correlation shift                                     */
    const float A[]               /* i  : A(z) unquantized for the 4 subframes                  */
)
{
    float refl[M+1];
    short tmp_active_flag;


    if( !st->Opt_AMR_WB )
    {
        /* fix explanation: after function dtx_fx, the "vad_flag"
           parameter can not be used for SID scheduling purposes any
           longer as dtx_fx can schedules active frames even if the
           initial analyzed vad_flag is 0 ) in the worst case without
           the fix an active frame could be classified as SID frame,
           quite/very unlikley though
        */
        tmp_active_flag = 0;

        if( (st->core_brate != SID_2k40) && (st->core_brate != 0)  ) /* Note,  core_brate_fx can be -1 */
        {
            tmp_active_flag = 1; ; /* requires active coding according  to dtx_fx logic  */
        }

        if( st->Opt_DTX_ON != 0  &&  tmp_active_flag == 0  &&  st->ini_frame > 3)

        {
            /* update the counter of consecutive inactive frames in DTX */
            st->consec_inactive++;
            if( st->consec_inactive > 5 )
            {
                st->consec_inactive = 5;
            }

            if( st->consec_inactive == 5 )
            {
                /* compute spectral tilt parameter */
                a2rc( &A[1], refl, M );

                if( st->spectral_tilt_reset == 1 )
                {
                    st->spectral_tilt_reset = 0;
                    st->running_avg = refl[0];
                    st->ra_deltasum = 0;
                }

                st->ra_deltasum += (0.80f * st->running_avg + 0.20f * refl[0]) - st->running_avg;
                st->running_avg = 0.80f * st->running_avg + 0.20f * refl[0];

                if( fabs(st->ra_deltasum) > 0.2f )
                {
                    st->spectral_tilt_reset = 1;
                    st->running_avg = 0;
                    st->ra_deltasum = 0;
                    st->trigger_SID = 1;
                }
            }
        }
        else
        {
            st->trigger_SID = 0;
            st->consec_inactive = 0;
        }

        if( st->trigger_SID == 1 )
        {
            if( st->cng_cnt >= 8 )
            {
                /* Declare SID frame due to spectral tilt changes */
                st->cnt_SID = 1;
                st->core_brate = SID_2k40;
                st->trigger_SID = 0;
            }
            else if ( st->core_brate == SID_2k40 )
            {
                /* SID fame has already been declared before */
                st->trigger_SID = 0;
            }
        }
    }


    if( (voicing[0] + voicing[1] + voicing[2]) / 3 + corr_shift > 0.65 &&
            (short)(abs(pitch[0] - st->pitO) + abs(pitch[1] - pitch[0]) + abs(pitch[2] - pitch[1])) / 3 < 14 )
    {
        (st->voiced_burst)++;
    }
    else
    {
        st->voiced_burst = 0;
    }

    /* Update previous voicing value for next frame use */
    st->voicing_old = (voicing[0] + voicing[1] + voicing[2]) / 3 + corr_shift;

    return;
}

