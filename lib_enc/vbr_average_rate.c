/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"


/*------------------------------------------------------------------*
 * Local constants
 *------------------------------------------------------------------*/

#define RATEWIN         600     /* length of the rate control window. This is 600 active speech frames. This equals roughly 12s of active speech */

/*------------------------------------------------------------------*
 * update_average_rate()
 *
 * SC-VBR update average data rate
 *------------------------------------------------------------------*/

void update_average_rate(
    Encoder_State *st           /* i/o: encoder state structure */
)
{
    float avratetarg;           /* target rate for next RATEWIN active frames */
    float target;               /* target set by VBR_ADR_MAX_TARGET*RATEWIN*10 */

    if ( st->numactive == RATEWIN )  /* goes into rate control only the numactive ==RATEWIN. So rate control is triggered after each RATEWIN avtive frames */
    {
        /* after 1000 blocks of RATEWIN frames, we change the way we control the average rate by using
           st->global_avr_rate=0.99*st->global_avr_rate+0.01*st->sum_of_rates. This will avoid
           veriables growing indefinitely while providing a good long term average rate */

        if ( st->global_frame_cnt < 1000 )
        {
            st->global_frame_cnt++;
            st->global_avr_rate = (st->global_avr_rate * (st->global_frame_cnt-1) + st->sum_of_rates) / st->global_frame_cnt;
        }
        else
        {
            st->global_avr_rate = 0.01f * st->sum_of_rates + 0.99f * st->global_avr_rate;
        }

        if ( st->sum_of_rates == 0 )
        {
            st->sum_of_rates = (float) (RATEWIN * VBR_ADR_MAX_TARGET * 10);
        }

        target = VBR_ADR_MAX_TARGET * 10 * RATEWIN;

        if ( target < st->global_avr_rate )  /* Action is taken to reduce the averge rate. Only initiated if the global rate > target rate */
        {
            /* Check the vad snr values to table the noisey/not noisey decision */

            if ( st->SNR_THLD < 67 )   /* Currently in QFF mode. The bumpup thresholds are slightly relaxed for noisy speech. */
            {
                /*  Increase the threshold so the the bumpup procedure is done using the noisy thresholds.
                    Use 3.5 steps to quickly ramp up the rate control to reduce the settling time */
                st->SNR_THLD += 3.5f;
            }
            else if ( st->mode_QQF == 0 && st->sum_of_rates > target ) /* Now SNR_THLD is in the max allowed. Sill the global average is higher and
                                                                          last RATEWIN frames have a higher agerage than the target rate. Now slightly
                                                                          more aggresive rate control is used by changing the mode to QQF. Still the
                                                                          same strict bumpups (more bumpups,higher rate) are used. */
            {
                /* Kick in QQF mode */
                st->mode_QQF = 1;
            }
            else if ( st->sum_of_rates > target )  /* Actions (1) and (2) are not sufficient to control the rate. Still the last RATEWIN active
                                                      frames have a higher average rate than the target rate. More aggresive rate control is
                                                      needed. At this point the rate_control flag is set. This will enable the more relaxed
                                                      bump up thresholds (less bump ups->reduced rate)*/
            {
                /* Relaxed bump ups are used */
                st->rate_control = 1;

                /* This will be triggered only if the gloabl average rate is considerablly higher than the target rate.
                   Keep a higher threshold to avoid short term rate increases over the target rate. */
                if ( st->global_avr_rate > (target+420.0f) ) /* Last resort rate control. This is a safer rate control mechanism by increasing NELPS */
                {
                    st->Last_Resort = 1;    /* compute based on a larger window as the last resort */
                }
                else
                {
                    st->Last_Resort = 0;
                }
            }
            else if ( st->sum_of_rates < target ) /* If the average rate of last RATEWIN frames is controlled by above actions, disable the most
                                                     aggresive rate control mechanisms. Still keep QQF mode as the global rate is not under
                                                     the target rate*/
            {
                st->Last_Resort = 0;
                st->mode_QQF = 1;
                st->rate_control = 0;
            }
        }
        else
        {
            /* floding back to lesser and leser aggresive rate control mechanisms gradually if global rate is under control */
            st->Last_Resort = 0;

            if ( st->rate_control == 1 )
            {
                st->rate_control = 0;
            }
            else if ( st->mode_QQF == 1) /* now rate control is not active and still the global rate is below the target. so go to QFF mode */
            {
                st->mode_QQF = 0;
            }
            else
            {
                if ( st->SNR_THLD >= 60 )
                {
                    st->SNR_THLD -= 1.5f;
                }
                else
                {
                    st->SNR_THLD = 60.0f;
                }
            }
        }

        if ( st->global_avr_rate < target-120 )  /* In QFF mode and global rate is less than target rate-0.2kbps. We can send some Q frames
                                                    to F frames to improve the quality */
        {
            /* kick in bouncing back from Q to F */
            st->Q_to_F = 1;

            /* average rate for next 600ms = global_rate * 2 - rate of the past RATEWIN active frames */
            avratetarg = (float)((RATEWIN * 10) * 2 * VBR_ADR_MAX_TARGET - st->global_avr_rate);

            /* compute the percentage of frames that needed to be sent to F. st->pattern_m is computed as % val * 1000. eg. if % is 10%, then
               st->pattern_m=100 . Later this value is used in voiced.enc to bump up 10% of PPP frames to F frames. */
            st->pattern_m = (short)(1000 * (avratetarg - 6.15f * RATEWIN * 10)/(10 * RATEWIN * 0.1f) );

            if ( st->pattern_m < 0 )
            {
                st->pattern_m = 0;      /* no bump up will ever happen */
            }

            if ( st->pattern_m > 1000 )
            {
                st->pattern_m = 1000;   /* 10% of bump ups */
            }

            st->patterncount = 0;
        }
        else
        {
            st->Q_to_F = 0;
        }

        st->sum_of_rates = 0;
        st->numactive = 0;
    }

    st->numactive++;

    /* sum the total number of bits (in kbytes) * 10 here */
    st->sum_of_rates += (st->core_brate / 1000.0f) * 10;

    return;
}


