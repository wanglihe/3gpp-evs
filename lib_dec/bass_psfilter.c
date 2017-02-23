/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <stdlib.h>
#include "prot.h"
#include "cnst.h"
#include "rom_com.h"
#include "options.h"
#include <assert.h>

/*---------------------------------------------------------------------*
 * Local constants
 *---------------------------------------------------------------------*/

#define NBPSF_L_EXTRA 120

#define BPF_STOP_STOPBAND_16 16


/*---------------------------------------------------------------------*
 * Local functions
 *---------------------------------------------------------------------*/

static short Pit_track( float syn[], short T );

/*---------------------------------------------------------------------*
 * bass_psfilter_init()
 *
 * Initialisation of postfiltering variables
 *---------------------------------------------------------------------*/

void bass_psfilter_init(
    float old_syn[],        /* o  : Old synthesis buffer 1        */
    float *mem_deemph_err,  /* o  : Error deemphasis memory       */
    float *lp_ener          /* o  : long_term error signal energy */
)
{
    /* post-filter memories */
    *mem_deemph_err = 0.0f;
    *lp_ener = 0.0f;
    set_f( old_syn, 0, NBPSF_PIT_MAX );

    return;
}

/*---------------------------------------------------------------------*
 * bass_psfilter()
 *
 * Perform low-frequency postfiltering
 *---------------------------------------------------------------------*/

void bass_psfilter(
    const short Opt_AMR_WB,      /* i  : AMR-WB IO flag                     */
    const float synth_in[],      /* i  : input synthesis (at 16kHz)         */
    const short L_frame,         /* i  : length of the last frame           */
    const float pitch_buf[],     /* i  : pitch for every subfr [0,1,2,3]    */
    float old_syn[],       /* i/o: NBPSF_PIT_MAX                      */
    float *mem_deemph_err, /* i/o: Error deemphasis memory            */
    float *lp_ener,        /* i/o: long_term error signal energy      */
    const short bpf_off,         /* i  : do not use BPF when set to 1       */
    float v_stab,          /* i  : stability factor                   */
    float *v_stab_smooth,  /* i/o: smoothed stability factor          */
    float *mem_mean_pit,   /* i/o: average pitch memory               */
    short *Track_on_hist,  /* i/o: History of half frame usage        */
    short *vibrato_hist,   /* i/o: History of frames declared as vibrato*/
    float *psf_att,        /* i/o: Post filter attenuation factor     */
    const short coder_type,      /* i  : coder_type                         */
    float bpf_noise_buf[]  /* o  : BPF error signal (at int_fs)       */
)
{
    short i, i_subfr, T;
    float gain, alpha, corr, ener;
    float syn_buf[NBPSF_PIT_MAX+L_FRAME16k+NBPSF_PIT_MAX], *syn;
    float syn2_buf[L_FRAME16k], *syn2;
    float err[L_HALFR16k];
    short T_sf[NB_SUBFR16k];
    short nb_subfr, T_update;
    float delta_v_stab;
    short dist_pit_diff, idx_pit_min, idx_pit_max, vibrato, Track_on;
    float loc_pit_max, loc_pit_min, diff_pit;
    float TrackOnR, vibratR, alp_tmp;

    Track_on = 0;
    T_update = 01;
    vibrato = 0;
    nb_subfr = L_frame/L_SUBFR;

    /*-------------------------------------------------------
     *  Initialize pointers to various synthesis buffers
     *
     *   |--------------------syn_buf--------------------------------|
     *   |-----old_syn-----|-------synth_in--------|----extrapol---- |
     *   |--NBPSF_PIT_MAX--| sf1 | sf2 | sf3 | sf4 |--NBPSF_PIT_MAX--|
     *                     |------syn2_buf---------|
     *                     |------L_frame----------|
     *                     |----bpf_noise_buf------|
     *
     *-------------------------------------------------------*/

    mvr2r( old_syn, syn_buf, NBPSF_PIT_MAX );
    mvr2r( synth_in, syn_buf + NBPSF_PIT_MAX, L_frame );

    if( !(pitch_buf == NULL || bpf_off) )
    {
        for(i = L_TRACK_HIST-1; i > 0; i--)
        {
            mem_mean_pit[i] = mem_mean_pit[i-1];
        }

        mem_mean_pit[i] = mean( pitch_buf, nb_subfr );
    }

    idx_pit_min = minimum(mem_mean_pit, L_TRACK_HIST, &loc_pit_min);
    idx_pit_max = maximum(mem_mean_pit, L_TRACK_HIST, &loc_pit_max);

    dist_pit_diff = (short)abs(idx_pit_max - idx_pit_min);
    diff_pit = loc_pit_max-loc_pit_min;

    if( L_frame == L_FRAME16k )
    {
        diff_pit *= 0.8f;
    }

    if( coder_type != INACTIVE && diff_pit >= 2.0f  && diff_pit < 10 && dist_pit_diff>= 3)
    {
        vibrato = 1;
    }

    TrackOnR = (float)sum_s(Track_on_hist,L_TRACK_HIST)/L_TRACK_HIST;
    vibratR = (float)sum_s(vibrato_hist,L_TRACK_HIST);
    alp_tmp = 1.0f-TrackOnR;

    if( vibrato )
    {
        vibratR = vibratR * vibratR * -0.009f + 1.0f;
        alp_tmp *= vibratR;
    }

    alp_tmp = max( 0.1f, alp_tmp );

    if( alp_tmp > *psf_att )
    {
        *psf_att = 0.05f * alp_tmp + 0.95f * *psf_att;
    }
    else
    {
        *psf_att = 0.4f * alp_tmp + 0.6f * *psf_att;
    }

    if( pitch_buf == NULL || bpf_off )
    {
        /* do not use BPF for HQ core */
        T_update = 80;
        set_s( T_sf, 0, 5 );
        syn = &syn_buf[NBPSF_PIT_MAX+L_frame];
        for (i=0; i<NBPSF_PIT_MAX; i++)
        {
            syn[i] = syn[i-T_update];
        }
    }
    else
    {
        /* extrapolation of synth_in */
        for( i=0; i<nb_subfr; i++ )
        {
            /* copy subframe pitch values [1,2,3,4] of the current frame */
            T_sf[i] = (short)(pitch_buf[i] + 0.5f);

            if(L_frame == L_FRAME16k)
            {
                /* Safety check, should be very rare */
                if( T_sf[i] > NBPSF_PIT_MAX )
                {
                    T_sf[i] = NBPSF_PIT_MAX;
                }
            }
            else
            {
                if( T_sf[i] > PIT_MAX )
                {
                    T_sf[i] = PIT_MAX;
                }
            }
        }

        T = T_sf[nb_subfr-1];
        syn = &syn_buf[NBPSF_PIT_MAX+L_frame];
        for( i=0; i<NBPSF_PIT_MAX; i++ )
        {
            syn[i] = syn[i-T];
        }
    }

    for( i_subfr=0; i_subfr<L_frame; i_subfr+=L_SUBFR )
    {
        T = T_sf[i_subfr/L_SUBFR];

        /* Pitch limitation already done in the previous loop */
        syn = &syn_buf[NBPSF_PIT_MAX+i_subfr];
        syn2 = &syn2_buf[i_subfr];

        if( T != 0 )
        {
            if( T>=PIT_MIN && Opt_AMR_WB )
            {
                T = Pit_track( syn, T );

                if( T != T_sf[i_subfr/L_SUBFR] )
                {
                    Track_on = 1;
                }
            }

            /* symetric pitch prediction : phase is opposite between harmonic */
            for( i=0; i<L_SUBFR; i++ )
            {
                syn2[i] = 0.5f*(syn[i-T]+syn[i+T]);
            }

            /* gain of prediction */
            corr = 0.0001f;
            for (i=0; i<L_SUBFR; i++)
            {
                corr += syn[i] * syn2[i];
            }

            ener = 0.0001f;
            for (i=0; i<L_SUBFR; i++)
            {
                ener += syn2[i] * syn2[i];
            }
            gain = corr / ener;

            /* error of prediction for noise estimator */
            for (i=0; i<L_SUBFR; i++)
            {
                err[i] = syn[i]-(gain*syn2[i]);
            }

            /* alpha = post-filtering factor (0=OFF, 1=100%) */
            ener = (float)(ener + pow( 10.0, 0.1 * *lp_ener) );       /* limit alpha (post-filtering) when noise is high */
            alpha = corr / ener;
            if (alpha > 0.5f)
            {
                alpha = 0.5f;
            }

            alpha *= (*psf_att);
            if (alpha > 0.3f && Track_on)
            {
                alpha = 0.3f;
            }
            else if (alpha > 0.4f && vibrato )
            {
                alpha = 0.4f;
            }

            if( alpha < 0.0f )
            {
                alpha = 0.0f;
            }

            *v_stab_smooth = 0.8f * v_stab + 0.2f * (*v_stab_smooth);
            delta_v_stab = (float) fabs(*v_stab_smooth - v_stab);
            v_stab = (float) pow(v_stab, 0.5f);
            alpha = (1.0f + 0.15f*v_stab - 2.0f*delta_v_stab) * alpha;

            for (i=0; i<L_SUBFR; i++)
            {
                syn2[i] = alpha*(syn[i]-syn2[i]);
            }
        }

        else
        {
            /* symetric pitch prediction : phase is opposite between harmonic */
            for( i=0; i<L_SUBFR; i++ )
            {
                syn2[i] = 0.5f*(syn[i-T_update]+syn[i+T_update]);
            }

            if( coder_type == AUDIO )  /* GSC mode without temporal component */
            {
                Track_on = 1;
            }
            else
            {
                Track_on = 0;
            }

            /* gain of prediction */
            corr = 0.0001f;
            for (i=0; i<L_SUBFR; i++)
            {
                corr += syn[i] * syn2[i];
            }
            ener = 0.0001f;
            for (i=0; i<L_SUBFR; i++)
            {
                ener += syn2[i] * syn2[i];
            }
            gain = corr / ener;

            /* error of prediction for noise estimator */
            for (i=0; i<L_SUBFR; i++)
            {
                err[i] = syn[i]-(gain*syn2[i]);
            }

            /* alpha = post-filtering factor (0=OFF, 1=100%) */
            alpha = 0.0f;

            *v_stab_smooth = 0.8f * v_stab + 0.2f * (*v_stab_smooth);

            set_f( syn2, 0.0f, L_SUBFR );
        }

        /* low-frequency noise estimator (lp_ener is average in dB) */
        deemph( err, 0.9f, L_SUBFR, mem_deemph_err );   /* +20dB at 50Hz */

        ener = 0.0001f;
        for( i=0; i<L_SUBFR; i++ )
        {
            ener += err[i]*err[i];
        }
        ener = (float)(10.0*log10(ener));
        *lp_ener = (float)(0.99f * *lp_ener + 0.01f * ener);

        /* just write out the error signal */
        mvr2r( syn2, bpf_noise_buf + i_subfr, L_SUBFR );
    }

    /*-------------------------------------------------------*
     * update memory for next frame
     *-------------------------------------------------------*/

    for( i = L_TRACK_HIST-1; i > 0; i-- )
    {
        Track_on_hist[i] = Track_on_hist[i-1];
        vibrato_hist[i] = vibrato_hist[i-1];
    }

    Track_on_hist[i] = Track_on;
    vibrato_hist[i] = vibrato;

    mvr2r( syn_buf+L_frame, old_syn, NBPSF_PIT_MAX );

    return;
}


/*---------------------------------------------------------------------*
 * Pit_track()
 *
 * Perform pitch tracking and test pitch/2 to avoid continuous pitch doubling
 *---------------------------------------------------------------------*/

static short Pit_track(/* o :  Pitch                                   */
    float syn[],       /* i :  synthesis [-NBPSF_PIT_MAX..L_HALFR16k]  */
    short T            /* i :  pitch period (>= PIT_MIN)               */
)
{
    short i, T2;
    float tmp, corr, ener, cn;
    float *v1, *v2;

    /*----------------------------------------------------------------*
     * Test pitch/2 to avoid continuous pitch doubling
     * (short pitch is limited to PIT_MIN (34 = 376Hz) by the encoder
     *----------------------------------------------------------------*/

    T2 = T>>1;

    v1 = &syn[-NBPSF_L_EXTRA];
    v2 = &syn[-T2-NBPSF_L_EXTRA];

    ener = 0.01f;
    for (i=0; i<L_HALFR16k+NBPSF_L_EXTRA; i++)
    {
        ener += v1[i]*v1[i];
    }
    corr = 0.01f;
    for (i=0; i<L_HALFR16k+NBPSF_L_EXTRA; i++)
    {
        corr += v1[i]*v2[i];
    }
    tmp = 0.01f;
    for (i=0; i<L_HALFR16k+NBPSF_L_EXTRA; i++)
    {
        tmp += v2[i]*v2[i];
    }

    /* cn = normalized correlation of pitch/2 */
    cn = corr / (float)sqrt(ener*tmp);
    if (cn > 0.95f)
    {
        T = T2;
    }

    return(T);
}


/*---------------------------------------------------------------------*
 * addBassPostFilter()
 *
 * Add BPF component in cldfb domain
 *---------------------------------------------------------------------*/

void addBassPostFilter(
    const float *harm_timeIn,
    int samplesToProcess,
    float **rAnalysis,
    float **iAnalysis,
    HANDLE_CLDFB_FILTER_BANK cldfb
)
{
    float *tmp_R[CLDFB_NO_COL_MAX];
    float *tmp_I[CLDFB_NO_COL_MAX];
    float cldfbBufferReal[CLDFB_NO_COL_MAX][20];
    float cldfbBufferImag[CLDFB_NO_COL_MAX][20];
    short i, b;
    int maxBand;
    const float *weights;
    int nCol = cldfb->no_col;
    int nColToProcess = nCol;
    int nChan = cldfb->no_channels;

    if( samplesToProcess > -1 )
    {
        nColToProcess = ((samplesToProcess + cldfb->no_channels - 1) / cldfb->no_channels);
    }

    assert(nCol == 16);

    weights = bpf_weights_16;

    if( nChan > BPF_STOP_STOPBAND_16 )
    {
        maxBand = BPF_STOP_STOPBAND_16;
    }
    else
    {
        maxBand = nChan;
    }

    for (i=0; i < nColToProcess; i++)
    {
        tmp_R[i] = cldfbBufferReal[i];
        tmp_I[i] = cldfbBufferImag[i];
    }

    /* do the analysis of filtered signal */
    cldfbAnalysis( harm_timeIn, tmp_R, tmp_I, samplesToProcess, cldfb );

    /* now do the subtraction */
    for( i = 0; i < nColToProcess; i++ )
    {
        /* loop over low frequency bands */
        for( b = 0; b < maxBand; b++ )
        {
            rAnalysis[i][b] -= weights[b] * tmp_R[i][b];
            iAnalysis[i][b] -= weights[b] * tmp_I[i][b];
        }
    }

    return;
}
