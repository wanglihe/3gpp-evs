/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <assert.h>
#include "options.h"
#include "cnst.h"
#include "rom_enc.h"
#include "rom_com.h"
#include "prot.h"

/*-------------------------------------------------------------------*
 * Local constants
 *-------------------------------------------------------------------*/

#define ALPHA_BWD           0.75f
#define BWD_LT_THRESH       0.6f

#define BWD_COUNT_MAX       100
#define BWD_COUNT_WIDER_BW  10

#define CLDFB_ENER_OFFSET 1.6f

/*-------------------------------------------------------------------*
 * bw_detect()
 *
 * bandwidth detector
 *-------------------------------------------------------------------*/

void bw_detect(
    Encoder_State *st,              /* i/o: Encoder State       */
    const float signal_in[],        /* i  : input signal        */
    const short localVAD,
    float *enerBuffer
)
{
    short i, j, k, bw_max, bin_width, n_bins;
    float spect[BWD_TOTAL_WIDTH], in_win[BWD_TOTAL_WIDTH];
    float spect_bin[BWD_N_BINS_MAX];
    float cldfb_bin[9];
    const float *pt, *pt1;
    float max_NB, max_WB, max_SWB, max_FB, mean_NB, mean_WB, mean_SWB, mean_FB;
    short cldfb_bin_width = 4;

    if( st->input_Fs > 8000 )
    {
        if( enerBuffer != NULL )
        {
            float ScalFac = 1.0f;

            ScalFac =  1/ ( st->cldfbAnaEnc->scale * st->cldfbAnaEnc->scale * 8.f);
            set_f( cldfb_bin, 0.001f, 9 );

            /* NB: 1.2 - 2.8 kHz, 4 cldfb-bands*/
            cldfb_bin[0] += sum_f( &(enerBuffer[3]), cldfb_bin_width );

            /* WB: 4.4 - 7.2 kHz, 8 cldfb-bands, mid band(14) counted twice */
            if( st->input_Fs >= 16000 )
            {
                cldfb_bin[1] += sum_f( &(enerBuffer[11]), cldfb_bin_width );
                cldfb_bin[2] += sum_f( &(enerBuffer[14]), cldfb_bin_width );
            }

            /* SWB: 9.2 - 15.6 kHz, 16 cldfb-bands */
            if( st->input_Fs >= 32000 )
            {
                cldfb_bin[3] += sum_f( &(enerBuffer[23]), cldfb_bin_width );
                cldfb_bin[4] += sum_f( &(enerBuffer[27]), cldfb_bin_width );
                cldfb_bin[5] += sum_f( &(enerBuffer[31]), cldfb_bin_width );
                cldfb_bin[6] += sum_f( &(enerBuffer[35]), cldfb_bin_width );
            }

            /* FB:  16.8 - 20.0 kHz, 8 cldfb-bands */
            if( st->input_Fs >= 48000 )
            {
                cldfb_bin[7] += sum_f( &(enerBuffer[42]), cldfb_bin_width );
                cldfb_bin[8] += sum_f( &(enerBuffer[46]), cldfb_bin_width );
            }

            for (i=0; i<9; i++)
            {
                cldfb_bin[i] = (float)log10( cldfb_bin[i]* ScalFac );  /* see formula used in perform_noise_estimation_enc() for CNG */
            }

        }
        else
        {

            /* set width of a speactral bin (corresponds to 1.5kHz) */
            if( st->input_Fs == 16000 )
            {
                bw_max = WB;
                bin_width = 60;
                n_bins = 5;                    /* spectrum to 7.5 kHz */
            }
            else if( st->input_Fs == 32000 )
            {
                bw_max = SWB;
                bin_width = 30;
                n_bins = 10;                    /* spectrum to 15 kHz */
            }
            else  /* st->input_Fs == 48000 */
            {
                bw_max = FB;
                bin_width = 20;
                n_bins = BWD_N_BINS_MAX;        /* spectrum to 19.5 kHz */
            }

            /*---------------------------------------------------------------------*
             * windowing of the input signal
             *---------------------------------------------------------------------*/

            pt = signal_in;
            pt1 = hann_window_320;

            /* 1st half of the window */
            for( i=0; i<BWD_TOTAL_WIDTH/2; i++ )
            {
                in_win[i] = *pt++ **pt1++;
            }
            pt1--;

            /* 2nd half of the window */
            for( ; i<BWD_TOTAL_WIDTH; i++ )
            {
                in_win[i] = *pt++ **pt1--;
            }

            /*---------------------------------------------------------------------*
             * tranform into frequency domain
             *---------------------------------------------------------------------*/

            edct( in_win, spect, BWD_TOTAL_WIDTH );

            /*---------------------------------------------------------------------*
             * compute energy per spectral bins
             *---------------------------------------------------------------------*/

            set_f( spect_bin, 0.001f, n_bins );

            for( k=0; k<=bw_max; k++ )
            {
                for( i=bwd_start_bin[k]; i<=bwd_end_bin[k]; i++ )
                {
                    for( j=0; j<bin_width; j++ )
                    {
                        spect_bin[i] += spect[i*bin_width + j]*spect[i*bin_width + j];
                    }
                    spect_bin[i] = (float)log10(spect_bin[i]);
                }
            }
        }

        if( enerBuffer != NULL )
        {
            /* cldfb detections */
            mean_NB = mean( cldfb_bin, 1 );                       /* NB: 1.2 - 2.8 kHz, 4 cldfb-bands (1 bin)   */
            maximum       ( cldfb_bin, 1 , &max_NB );
            mean_WB = mean( cldfb_bin + 1, 2 );                   /* WB: 4.4 - 7.2 kHz, 8 cldfb-bands (2 bins)  */
            maximum       ( cldfb_bin + 1, 2 , &max_WB );

            mean_NB  += CLDFB_ENER_OFFSET;
            max_NB   += CLDFB_ENER_OFFSET;
            mean_WB  += CLDFB_ENER_OFFSET;
            max_WB   += CLDFB_ENER_OFFSET;

            if( st->input_Fs == 16000 )
            {
                /* for 16kHz sampled inputs, do not check SWB & FB */
                mean_SWB = 0.0f;
                max_SWB = 0.0f;
                mean_FB = 0.0f;
                max_FB = 0.0f;
            }
            else if( st->input_Fs == 32000 )
            {
                /* for 32kHz sampled inputs, do not check FB */
                mean_FB = 0.0f;
                max_FB = 0.0f;
                mean_SWB = mean( cldfb_bin + 3, 4 );               /*  SWB: 9.2 - 15.6 kHz, 16 cldfb-bands (4 bins) */
                maximum        ( cldfb_bin + 3, 4 , &max_SWB );
                mean_SWB += CLDFB_ENER_OFFSET;
                max_SWB  += CLDFB_ENER_OFFSET;
            }
            else
            {
                mean_SWB = mean( cldfb_bin + 3, 4 ) ;              /*  SWB: 9.2 - 15.6 kHz, 16 cldfb-bands (4 bins) */
                maximum        ( cldfb_bin + 3, 4 , &max_SWB );
                mean_FB = mean ( cldfb_bin + 7, 2 );               /*  FB: 16.8 - 20.0 kHz, 8 cldfb-bands (2 bins) */
                maximum        ( cldfb_bin + 7, 2 , &max_FB );

                mean_SWB += CLDFB_ENER_OFFSET;
                max_SWB  += CLDFB_ENER_OFFSET;
                mean_FB  += CLDFB_ENER_OFFSET;
                max_FB   += CLDFB_ENER_OFFSET;
            }

        }
        else
        {
            mean_NB = mean( spect_bin + bwd_start_bin[0], bwd_end_bin[0]-bwd_start_bin[0]+1 );            /* NB:  1.5-3.0kHz (1 bin)   */
            maximum       ( spect_bin + bwd_start_bin[0], bwd_end_bin[0]-bwd_start_bin[0]+1, &max_NB );
            mean_WB = mean( spect_bin + bwd_start_bin[1], bwd_end_bin[1]-bwd_start_bin[1]+1 );            /* WB:  4.5-7.5kHz (2 bins)  */
            maximum(        spect_bin + bwd_start_bin[1], bwd_end_bin[1]-bwd_start_bin[1]+1, &max_WB );

            if( st->input_Fs == 16000 )
            {
                /* for 16kHz sampled inputs, do not check SWB & FB */
                mean_SWB = 0.0f;
                max_SWB = 0.0f;
                mean_FB = 0.0f;
                max_FB = 0.0f;
            }
            else if( st->input_Fs == 32000 )
            {
                mean_SWB = mean( spect_bin + bwd_start_bin[2], bwd_end_bin[2]-bwd_start_bin[2]+1 );       /* SWB: 9.0-15.0kHz (4 bins) */
                maximum( spect_bin + bwd_start_bin[2], bwd_end_bin[2]-bwd_start_bin[2]+1, &max_SWB );

                /* for 32kHz sampled inputs, do not check FB */
                mean_FB = 0.0f;
                max_FB = 0.0f;
            }
            else
            {
                mean_SWB = mean( spect_bin + bwd_start_bin[2], bwd_end_bin[2]-bwd_start_bin[2]+1 );       /* SWB: 9.0-15.0kHz (4 bins) */
                maximum( spect_bin + bwd_start_bin[2], bwd_end_bin[2]-bwd_start_bin[2]+1, &max_SWB );
                mean_FB = mean( spect_bin + bwd_start_bin[3], bwd_end_bin[3]-bwd_start_bin[3]+1 );        /* FB: 16.5-19.5kHz (2 bins) */
                maximum( spect_bin + bwd_start_bin[3], bwd_end_bin[3]-bwd_start_bin[3]+1, &max_FB );
            }
        }

        /*---------------------------------------------------------------------*
         * update LT counters and energies
         *---------------------------------------------------------------------*/

        if( localVAD || st->lp_noise > 30 )
        {
            st->lt_mean_NB = ALPHA_BWD * st->lt_mean_NB + (1-ALPHA_BWD) * mean_NB;
            st->lt_mean_WB = ALPHA_BWD * st->lt_mean_WB + (1-ALPHA_BWD) * mean_WB;
            st->lt_mean_SWB = ALPHA_BWD * st->lt_mean_SWB + (1-ALPHA_BWD) * mean_SWB;

            if( enerBuffer != NULL )
            {
                if( 0.9f * max_WB > BWD_LT_THRESH * st->lt_mean_NB )
                {
                    if( 2.5f * max_WB > max_NB )
                    {
                        st->count_WB++;
                    }
                }
                else
                {
                    if( 3.5f * mean_WB < mean_NB )
                    {
                        st->count_WB--;
                    }
                }

                if( 0.83f * max_SWB > BWD_LT_THRESH * st->lt_mean_WB && max_WB > BWD_LT_THRESH * st->lt_mean_NB )
                {
                    if( 2 * max_SWB > max_WB )
                    {
                        st->count_SWB++;
                    }
                }
                else
                {
                    if( 3 * mean_SWB < mean_WB )
                    {
                        st->count_SWB--;
                    }
                }

                if( max_FB > BWD_LT_THRESH * st->lt_mean_SWB && 0.83f * max_SWB > BWD_LT_THRESH * st->lt_mean_WB && max_WB > BWD_LT_THRESH * st->lt_mean_NB )
                {
                    if( 3 * max_FB > max_SWB )
                    {
                        st->count_FB++;
                    }
                }
                else
                {
                    if( 4.1f * mean_FB < mean_SWB )
                    {
                        st->count_FB--;
                    }
                }

            }
            else
            {
                if( max_WB > BWD_LT_THRESH * st->lt_mean_NB )
                {
                    if( 2 * max_WB > max_NB )
                    {
                        st->count_WB++;
                    }
                }
                else
                {
                    if( 2.6f * mean_WB < mean_NB )
                    {
                        st->count_WB--;
                    }
                }

                if( max_SWB > BWD_LT_THRESH * st->lt_mean_WB && max_WB > BWD_LT_THRESH * st->lt_mean_NB )
                {
                    if( 2 * max_SWB > max_WB )
                    {
                        st->count_SWB++;
                    }
                }
                else
                {
                    if( 3 * mean_SWB < mean_WB )
                    {
                        st->count_SWB--;
                    }
                }

                if( max_FB > BWD_LT_THRESH * st->lt_mean_SWB && max_SWB > BWD_LT_THRESH * st->lt_mean_WB && max_WB > BWD_LT_THRESH * st->lt_mean_NB )
                {
                    if( 2 * max_FB > max_SWB )
                    {
                        st->count_FB++;
                    }
                }
                else
                {
                    if( 3 * mean_FB < mean_SWB )
                    {
                        st->count_FB--;
                    }
                }
            }

            st->count_WB  = min(st->count_WB,BWD_COUNT_MAX);
            st->count_SWB = min(st->count_SWB,BWD_COUNT_MAX);
            st->count_FB  = min(st->count_FB,BWD_COUNT_MAX);
            st->count_WB  = max(st->count_WB,0);
            st->count_SWB = max(st->count_SWB,0);
            st->count_FB  = max(st->count_FB,0);

            /*---------------------------------------------------------------------*
             * check against thresholds
             * detect a band-width change
             *---------------------------------------------------------------------*/

            /* switching to a higher BW */
            if( st->last_input_bwidth == NB )
            {
                if( st->count_WB > BWD_COUNT_WIDER_BW )
                {
                    st->input_bwidth = WB;
                    st->count_WB = BWD_COUNT_MAX;

                    if( st->count_SWB > BWD_COUNT_WIDER_BW )
                    {
                        st->input_bwidth = SWB;
                        st->count_SWB = BWD_COUNT_MAX;

                        if( st->count_FB > BWD_COUNT_WIDER_BW )
                        {
                            st->input_bwidth = FB;
                            st->count_FB = BWD_COUNT_MAX;
                        }
                    }
                }
            }

            if( st->last_input_bwidth == WB && st->input_Fs > 16000 )
            {
                if( st->count_SWB > BWD_COUNT_WIDER_BW )
                {
                    st->input_bwidth = SWB;
                    st->count_SWB = BWD_COUNT_MAX;

                    if( st->count_FB > BWD_COUNT_WIDER_BW )
                    {
                        st->input_bwidth = FB;
                        st->count_FB = BWD_COUNT_MAX;
                    }
                }
            }

            if( st->last_input_bwidth == SWB && st->input_Fs > 32000 )
            {
                if( st->count_FB > BWD_COUNT_WIDER_BW )
                {
                    st->input_bwidth = FB;
                    st->count_FB = BWD_COUNT_MAX;
                }
            }

            /* switching to a lower BW */
            if( st->last_input_bwidth == FB )
            {
                if( st->count_FB < 10 )
                {
                    st->input_bwidth = SWB;
                    st->count_FB = 0;
                }
                if( st->count_SWB < 10 )
                {
                    st->input_bwidth = WB;
                    st->count_SWB = 0;
                    st->count_FB = 0;
                }
                if( st->count_WB < 10 )
                {
                    st->input_bwidth = NB;
                    st->count_WB = 0;
                    st->count_SWB = 0;
                    st->count_FB = 0;
                }
            }

            if( st->last_input_bwidth == SWB )
            {
                if( st->count_SWB < 10 )
                {
                    st->input_bwidth = WB;
                    st->count_SWB = 0;
                    st->count_FB = 0;
                }
                if( st->count_WB < 10 )
                {
                    st->input_bwidth = NB;
                    st->count_WB = 0;
                    st->count_SWB = 0;
                    st->count_FB = 0;
                }
            }

            if( st->last_input_bwidth == WB )
            {
                if( st->count_WB < 10 )
                {
                    st->input_bwidth = NB;
                    st->count_WB = 0;
                    st->count_SWB = 0;
                    st->count_FB = 0;
                }
            }
        }
    }

    /* verify that maximum encoded bandwidth (specified on the command line) is not exceeded */
    if( st->input_bwidth > st->max_bwidth )
    {
        st->input_bwidth = st->max_bwidth;
    }


    /* Set and limit the encoded bandwidth */
    if ( st->codec_mode == MODE1 )
    {
        long total_brate;

        st->bwidth = st->input_bwidth;
        total_brate = st->total_brate;
        if ( total_brate <= ACELP_9k60 && st->bwidth != NB && st->bwidth != WB )
        {
            st->bwidth = WB;
        }
        else if ( total_brate >= ACELP_13k20 && total_brate <= ACELP_16k40 && st->bwidth > SWB )
        {
            st->bwidth = SWB;
        }
        else if ( total_brate >= ACELP_32k && st->bwidth < WB )
        {
            st->bwidth = WB;
        }
    }
    else
    {
        short n, bits_frame_nominal, tmpBandwidthMin;

        bits_frame_nominal = st->total_brate / 50;
        for( n=0; n<FRAME_SIZE_NB; n++ )
        {
            if( FrameSizeConfig[n].frame_bits == bits_frame_nominal )
            {
                break;
            }
        }
        if( n == FRAME_SIZE_NB )
        {
            assert(!"Bitrate not supported: not part of EVS");
        }

        tmpBandwidthMin = FrameSizeConfig[n].bandwidth_min;

        if( st->rf_mode )
        {
            tmpBandwidthMin = WB;
        }

        st->bwidth = max(min(st->input_bwidth, FrameSizeConfig[n].bandwidth_max), tmpBandwidthMin);
    }

    return;
}
