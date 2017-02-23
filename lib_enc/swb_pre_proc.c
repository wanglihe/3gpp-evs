/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <stdlib.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"
#include "rom_enc.h"


/*-------------------------------------------------------------------*
 * Local constants
 *
 *-------------------------------------------------------------------*/

#define CLDFB_NO_CHANNELS_HB 20


/*-------------------------------------------------------------------*
 * wb_pre_proc()
 *
 * - Resampling of input signal when input signal sampling rate
 *   is above 16kHz
 * - Common WB TBE and WB BWE pre-processing
 *-------------------------------------------------------------------*/

void wb_pre_proc(
    Encoder_State *st,                  /* i/o: encoder state structure             */
    const float *new_inp_resamp16k,   /* i  : original input signal               */
    float *hb_speech            /* o  : HB target signal (6-8kHz) at 16kHz  */
)
{
    short Sample_Delay_WB_BWE, ramp_flag;
    float old_input[NS2SA(16000, DELAY_FD_BWE_ENC_NS + DELAY_FIR_RESAMPL_NS) + L_FRAME16k];
    float *highband_new_speech, highband_old_speech [(L_LOOK_12k8 + L_SUBFR + L_FRAME) * 5/16];
    short fSwitchFromIO = 0;

    if ( (st->last_total_brate == ACELP_6k60) ||
            (st->last_total_brate == ACELP_8k85) ||
            (st->last_total_brate == ACELP_12k65) ||
            (st->last_total_brate == ACELP_14k25) ||
            (st->last_total_brate == ACELP_15k85) ||
            (st->last_total_brate >= ACELP_18k25 && st->last_total_brate <= ACELP_23k85) )
    {
        fSwitchFromIO = 1;
    }

    set_f( old_input, 0, NS2SA(16000, DELAY_FD_BWE_ENC_12k8_NS + DELAY_FIR_RESAMPL_NS) + L_FRAME16k );

    if ( st->extl == WB_BWE || st->extl == WB_TBE || st->igf )
    {
        ramp_flag = 0;
        if( (st->last_extl != WB_TBE && st->last_extl != WB_BWE && !st->igf) || (st->igf && fSwitchFromIO) )
        {
            ramp_flag = 1;
        }

        if ( !st->ppp_mode)
        {
            flip_spectrum_and_decimby4( new_inp_resamp16k, hb_speech, L_FRAME16k, st->decim_state1, st->decim_state2, ramp_flag );

            if( st->extl != WB_TBE )
            {
                /* Update the previous wideband speech buffer in case of a WB_BWE frame - this code is in wb_tbe_enc */
                Sample_Delay_WB_BWE = (L_LOOK_12k8 + L_SUBFR) * 5/16;

                highband_new_speech = highband_old_speech + Sample_Delay_WB_BWE;
                mvr2r( hb_speech, highband_new_speech, L_FRAME16k / 4 );
                mvr2r( highband_old_speech + L_FRAME16k / 4, st->old_speech_wb, Sample_Delay_WB_BWE );
            }
        }
    }
    else
    {
        set_f( st->decim_state1, 0.0f, (2*ALLPASSSECTIONS_STEEP+1) );
        set_f( st->decim_state2, 0.0f, (2*ALLPASSSECTIONS_STEEP+1) );
        set_f( st->old_speech_wb, 0.0f, (L_LOOK_12k8 + L_SUBFR) * 5/16 );
    }

    /* st->old_input_wb and st->old_wtda_wb must be updated each frame, or there are often some clicks during WB TBE <-> WB BWE switching */
    if ( (st->extl != WB_BWE || (st->extl == WB_BWE && st->total_brate <= ACELP_8k00)) && !st->ppp_mode )
    {
        Sample_Delay_WB_BWE = NS2SA( 16000, DELAY_FD_BWE_ENC_12k8_NS );

        mvr2r( new_inp_resamp16k, &old_input[Sample_Delay_WB_BWE], L_FRAME16k );
        mvr2r( st->old_input_wb, old_input, Sample_Delay_WB_BWE );
        mvr2r( new_inp_resamp16k + L_FRAME16k - Sample_Delay_WB_BWE, st->old_input_wb, Sample_Delay_WB_BWE );
        if ((st->extl != SWB_BWE) && (st->extl != FB_BWE))
        {
            mvr2r( old_input, st->old_wtda_swb, L_FRAME16k );
        }
    }

    return;
}


/*-------------------------------------------------------------------*
 * swb_pre_proc()
 *
 * - Calculate the 6 to 14 kHz (or 7.5 - 15.5 kHz) SHB target signal
 *   for SWB TBE or SWB BWE coding
 * - Common SWB TBE and SWB BWE pre-processing
 *-------------------------------------------------------------------*/

void swb_pre_proc(
    Encoder_State *st,               /* i/o: encoder state structure                */
    const float *input,            /* i  : original input signal                  */
    float *new_swb_speech,   /* o  : original input signal at 32kHz         */
    float *shb_speech,       /* o  : SHB target signal (6-14kHz) at 16kHz   */
    float realBuffer[CLDFB_NO_COL_MAX][CLDFB_NO_CHANNELS_MAX], /* i : real buffer */
    float imagBuffer[CLDFB_NO_COL_MAX][CLDFB_NO_CHANNELS_MAX]  /* i : imag buffer */
)
{
    short Sample_Delay_SWB_BWE, inner_frame, delay;
    long inner_Fs;
    float old_input[NS2SA(48000, DELAY_FD_BWE_ENC_NS + DELAY_FIR_RESAMPL_NS) + L_FRAME48k];
    float spchTmp[640];
    short i, j;
    short startB, endB;
    float *realBufferFlipped[CLDFB_NO_COL_MAX];
    float *imagBufferFlipped[CLDFB_NO_COL_MAX];
    float realBufferTmp[CLDFB_NO_COL_MAX][20];
    float imagBufferTmp[CLDFB_NO_COL_MAX][20];
    short ts, nB, uB;
    float sign;

    for( j=0; j < CLDFB_NO_COL_MAX; j++ )
    {
        set_f( realBufferTmp[j], 0, 20 );
        set_f( imagBufferTmp[j], 0, 20 );
        realBufferFlipped[j] = realBufferTmp[j];
        imagBufferFlipped[j] = imagBufferTmp[j];
    }

    set_f( old_input, 0.0f, NS2SA(48000, DELAY_FD_BWE_ENC_12k8_NS + DELAY_FIR_RESAMPL_NS) + L_FRAME48k );

    if( st->input_Fs == 32000 )
    {
        mvr2r( input, new_swb_speech, L_FRAME32k );
        if( st->last_extl != SWB_BWE && st->last_extl != FB_BWE && st->extl != SWB_BWE_HIGHRATE)
        {
            Sample_Delay_SWB_BWE = NS2SA( 32000, DELAY_FD_BWE_ENC_12k8_NS + DELAY_FIR_RESAMPL_NS );
            mvr2r( st->old_fdbwe_speech, &old_input[Sample_Delay_SWB_BWE], L_FRAME32k );
            set_f( old_input, 0, Sample_Delay_SWB_BWE );
            mvr2r( st->old_fdbwe_speech + L_FRAME32k - Sample_Delay_SWB_BWE, st->old_input, Sample_Delay_SWB_BWE );
            mvr2r( old_input, st->old_wtda_swb, L_FRAME32k );
        }

        if( st->extl != SWB_BWE && st->extl != FB_BWE )
        {
            mvr2r( input, st->old_fdbwe_speech, L_FRAME32k );
        }
    }
    else  /* 48 kHz */
    {
        if( st->codec_mode == MODE1 )
        {
            if( st->extl != SWB_BWE && st->extl != FB_BWE && st->core == ACELP_CORE)
            {
                /* move the resampling out of the TDBWE path as new_swb_speech is not needed for TDBWE. */
                mvr2r( input, st->old_fdbwe_speech, L_FRAME48k );
            }
            else
            {
                if( st->last_extl != SWB_BWE && st->last_extl != FB_BWE )
                {
                    /* resample 48 kHz to 32kHz */
                    if( st->last_bwidth == FB )
                    {
                        inner_frame = L_FRAME48k;
                        inner_Fs = 48000;
                        mvr2r( st->old_fdbwe_speech, new_swb_speech, L_FRAME48k );
                    }
                    else
                    {
                        inner_frame = L_FRAME32k;
                        inner_Fs = 32000;
                        decimate_2_over_3_allpass( st->old_fdbwe_speech, L_FRAME48k, new_swb_speech, st->dec_2_over_3_mem, allpass_poles_3_ov_2,
                                                   decimate_3_ov_2_lowpass_num, decimate_3_ov_2_lowpass_den, st->dec_2_over_3_mem_lp );
                    }

                    Sample_Delay_SWB_BWE = NS2SA( inner_Fs, DELAY_FD_BWE_ENC_12k8_NS + DELAY_FIR_RESAMPL_NS );
                    mvr2r( new_swb_speech, &old_input[Sample_Delay_SWB_BWE], inner_frame );
                    set_f( old_input, 0, Sample_Delay_SWB_BWE );
                    mvr2r( new_swb_speech + inner_frame - Sample_Delay_SWB_BWE, st->old_input, Sample_Delay_SWB_BWE );
                    mvr2r( old_input, st->old_wtda_swb, inner_frame );
                }
                /* resample 48 kHz to 32kHz */
                if( st->bwidth == FB )
                {
                    mvr2r( input, new_swb_speech, L_FRAME48k );
                }
                else
                {
                    decimate_2_over_3_allpass( input, L_FRAME48k, new_swb_speech, st->dec_2_over_3_mem, allpass_poles_3_ov_2,
                                               decimate_3_ov_2_lowpass_num, decimate_3_ov_2_lowpass_den, st->dec_2_over_3_mem_lp );
                }
            }
        }
        else
        {
            /* resample 48 kHz to 32kHz */
            if( st->bwidth == FB )
            {
                mvr2r( input, new_swb_speech, L_FRAME48k );
            }
            else
            {
                decimate_2_over_3_allpass( input, L_FRAME48k, new_swb_speech, st->dec_2_over_3_mem, allpass_poles_3_ov_2,
                                           decimate_3_ov_2_lowpass_num, decimate_3_ov_2_lowpass_den, st->dec_2_over_3_mem_lp );
            }
        }
    }

    if( ( st->core == ACELP_CORE && st->extl != SWB_BWE_HIGHRATE && st->extl != FB_BWE_HIGHRATE ) ||
            ( ( st->total_brate == ACELP_9k60 || st->rf_mode ) && st->bwidth == SWB ) )
    {
        if( st->L_frame == L_FRAME )
        {
            startB= 34;
            endB= 14;
            for( ts = 0; ts < CLDFB_NO_COL_MAX; ts++ )
            {
                for( nB = startB, uB=0; nB > endB; nB--,uB++ )
                {
                    sign = (ts%2) ? 1.0f : -1.0f;
                    realBufferFlipped[ts][uB] = -sign*realBuffer[ts][nB];
                    imagBufferFlipped[ts][uB] =  sign*imagBuffer[ts][nB];
                }
            }
        }
        else
        {
            startB = 39;
            endB = 19;
            for( ts = 0; ts < CLDFB_NO_COL_MAX; ts++ )
            {
                for( nB = startB, uB=0; nB > endB; nB--,uB++ )
                {
                    realBufferFlipped[ts][uB] = -realBuffer[ts][nB];
                    imagBufferFlipped[ts][uB] =  imagBuffer[ts][nB];
                }
            }
        }

        {
            float CldfbHB = 0;
            for (nB = 0; nB < 10; nB++)
            {
                for (ts = 0; ts < CLDFB_NO_COL_MAX; ts++)
                {
                    CldfbHB +=  (realBufferFlipped[ts][nB] * realBufferFlipped[ts][nB] + imagBufferFlipped[ts][nB] * imagBufferFlipped[ts][nB]);
                }
            }
            if( CldfbHB <= 0 )
            {
                CldfbHB = 1.0f;
            }
            st->cldfbHBLT = 0.9f * st->cldfbHBLT + 0.1f * ( 0.221462f /*=1/log10(32768)*/ * (log10(CldfbHB) - 1.0f) );
        }
        cldfbSynthesis( realBufferFlipped, imagBufferFlipped, shb_speech, -1, st->cldfbSynTd );

        if( st->extl != WB_TBE && st->extl != SWB_TBE && st->extl != FB_TBE )
        {
            /* Update the previous superwideband speech buffer in case of a SWB_BWE frame - this code is in swb_tbe_enc */
            delay = L_LOOK_16k + L_SUBFR16k;
            mvr2r( shb_speech + L_FRAME16k - delay, st->old_speech_shb, delay );
        }
    }
    else
    {
        if( st->bwidth == FB || st->core == ACELP_CORE)
        {
            set_f( st->old_speech_shb, 0, L_LOOK_16k + L_SUBFR16k );
            set_f( shb_speech, 0, L_FRAME16k );   /* shb_speech for FB/SWB BWE_HIGHRATE is not used at 64kbps */
        }
        else
        {
            /* flip the spectrm */
            mvr2r( new_swb_speech, spchTmp, L_FRAME32k );

            for( i = 0; i < L_FRAME32k; i = i+2 )
            {
                spchTmp[i] = -spchTmp[i];
            }

            Decimate_allpass_steep( spchTmp, st->state_ana_filt_shb, L_FRAME32k, shb_speech );
            mvr2r( shb_speech + L_FRAME16k - (L_LOOK_16k + L_SUBFR16k), st->old_speech_shb, L_LOOK_16k + L_SUBFR16k );
        }

        /* Reset CLDFB synthesis buffer */
        set_f( st->cldfbSynTd->cldfb_state, 0.0f, st->cldfbSynTd->p_filter_length + st->cldfbSynTd->no_channels*st->cldfbSynTd->no_col );
    }

    /* Memory reset to compensate for 0.9375 ms offset when transitioning from IO to SWB */
    if( st->last_extl == -1 )
    {
        delay = NS2SA(st->input_Fs, DELAY_FIR_RESAMPL_NS);
        for( i = 0; i < delay; i++ )
        {
            shb_speech[i] = (float)i * (0.03f * shb_speech[2*delay-1-i]);
        }
    }

    return;
}
