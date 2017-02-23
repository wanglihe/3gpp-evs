/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"

/*-------------------------------------------------------------------*
 * Local constants
 *-------------------------------------------------------------------*/

#define MAX_SNR1       45.0f
#define INV_MAX_SNR    (1.f / (MAX_SNR1-1.0f))   /* max. SNR considered for noise subtraction in voiced segments */
#define MAX_SNR_SNR1   (MAX_SNR1 * INV_MAX_SNR)  /* 45 * (1 / (MAX_SNR1-1)) */


#define BIN_1KHZ       (short)(1000.f/BIN_16kdct)
#define BIN_2KHZ       (short)(2000.f/BIN_16kdct)
#define BIN_4KHZ       (short)(4000.f/BIN_16kdct)

#define MAX_GN_R       0.2f
#define ALPH           1.00f
#define BET            (1.925f-ALPH)
#define MAXX           5

/*-------------------------------------------------------------------*
 * Local functions
 *-------------------------------------------------------------------*/

static void spectrum_mod_dct( float data[], const float lf_E[], float lf_EO[],
                              const float noiseE[], const float minGain, float lp_gbin[],
                              const short music_flag, short min_band, const float MAX_GN, const short max_band );

static void analy_sp_dct( const float dct_in[], float dct_buf[], float *fr_bands, float *lf_E, float *etot );

static void find_enr_dct( const float data[], float band[], float *ptE, float *Etot, const short min_band,
                          const short max_band, float *Bin_E, const float bin_freq );

/*------------------------------------------------------------------------*
 * LD_music_post_filter()
 *
 * Music post-filter
 *------------------------------------------------------------------------*/

void LD_music_post_filter
(
    const float dtc_in[],         /* i   : input synthesis                            */
    float dtc_out[],        /* o   : output synthesis                           */
    const  long core_brate,       /* i   : core bitrate                               */
    short *last_music_flag, /* i/o : Previous music detection ouptut            */
    float *thresh,          /* i/o : Detection thresold                         */
    short *nb_thr_1,        /* i/o : Number of consecutives frames of level 1   */
    short *nb_thr_3,        /* i/o : Number of consecutives frames of level 3   */
    float *lt_diff_etot,    /* i/o : Long term total energy variation           */
    float *mem_etot,        /* i/o : Total energy memory                        */
    const float min_ns_gain,      /* i   : minimum gain for inter-harm noise red.     */
    float bckr[],           /* i/o : per band bckgnd. noise energy estimate     */
    float lf_EO[],          /* i/o : old per bin E for previous half frame      */
    float lp_gbin[],        /* i/o : smoothed suppression gain, per FFT bin     */
    float *filt_lfE,        /* i   : post filter weighting coefficient          */
    short *last_nonfull_music,/* i : Number of frames sinces la "speech like" frame */
    const short coder_type        /* i   : Coder type : -1 in case of IO              */
    ,const short Last_coder_type   /* i   : last Coder type                            */
)
{
    float DCT_buf[DCT_L_POST];
    float fr_bands[MBANDS_GN_LD];
    float lf_E[VOIC_BINS_HR];
    float etot, ftmp;
    short i, j, k;
    short min_band = 0;
    float local_min_gain = min_ns_gain;
    short music_flag2 = 0;
    float max_val;
    float max_ovf_2k, max_ovf_4k, max_ovf_6k;
    float min_g_2k, min_g_4k, min_g_6k;
    float m_ave, old_ftmp;
    float old_ftmp_1;
    float tmp_lfE[DCT_L_POST];
    float MAX_GN = MAX_GN_R;
    short MAX_band = MBANDS_GN_LD;

    /*------------------------------------------------------------------------*
     * Frequency analysis
     *------------------------------------------------------------------------*/

    analy_sp_dct( dtc_in, DCT_buf, fr_bands, lf_E, &etot );

    /*------------------------------------------------------------------------*
     * Find signal classification
     *------------------------------------------------------------------------*/

    music_flag2 = stab_est( etot, lt_diff_etot, mem_etot,
                            nb_thr_3, nb_thr_1, thresh, last_music_flag, 1 );

    if ( core_brate < ACELP_6k60 || Last_coder_type != AUDIO )
    {
        /* do not perform music improvement on SID frames */
        music_flag2 = 0;
    }

    if( music_flag2 < 4 )
    {
        *last_nonfull_music = 0;
    }
    else
    {
        (*last_nonfull_music)++;
    }

    *last_nonfull_music = min( 51, *last_nonfull_music );

    /*------------------------------------------------------------------------*
     * Remapping of bands
     * Section to "remap" the minimal band and the minimum gain for our needs
     *------------------------------------------------------------------------*/

    if( music_flag2 > 3)
    {
        min_band = 2;
        local_min_gain = 0.25119f;
    }
    else if( music_flag2 == 3)
    {
        min_band = 3;
        local_min_gain = 0.25119f;
    }
    else if( music_flag2 == 2)
    {
        min_band = 4;
        local_min_gain = 0.35481f;
    }
    else if( music_flag2 == 1)
    {
        min_band = 4;
        local_min_gain = 0.50119f;
    }

    min_band += 4;

    MAX_GN = 0.1f;
    if( core_brate > ACELP_9k60 )
    {
        /* overshoot not allowed, since GSC already matches the energy */
        MAX_GN = 0.0f;
    }

    if( coder_type == AUDIO )
    {
        /* with GSC we know for sure that we are in music */
        min_band = min( min_band, 3 );
    }

    /*------------------------------------------------------------------------*
     * Approximation of the inter-harmonic noise level
     * - sort the bin energy
     * - compupte the average energy per band excluding the maximum energy bin
     *------------------------------------------------------------------------*/

    j = 0;
    for (i = 0; i < MBANDS_GN_LD; i++)
    {
        m_ave = 0.0f;
        max_val = 0.0f;

        for( k=j; k < mfreq_bindiv_LD[i]+j; k++ )
        {
            m_ave += lf_E[k];
            max_val = max(max_val, lf_E[k]);
        }

        m_ave -= max_val;
        m_ave /= mfreq_bindiv_LD[i] - 1;

        bckr[i] = m_ave*sc_qnoise[i];
        j += mfreq_bindiv_LD[i];
    }

    i = maximum(lf_E, DCT_L_POST, &m_ave);

    /*------------------------------------------------------------------------*
     * - Normalisation of the energy vector between [0.72, 5], with the form of pow(x,4)
     * - Simple LP filtering along the frequency domain
     * - LT averaging with the past and in function of the stability factor
     *------------------------------------------------------------------------*/

    m_ave = ALPH/lf_E[i];

    ftmp = lf_E[0]*m_ave+BET;
    ftmp *= ftmp;
    ftmp *= ftmp;
    ftmp *= ftmp;
    ftmp = min(ftmp, MAXX);
    old_ftmp = ftmp;

    ftmp = lf_E[1]*m_ave+BET;
    ftmp *= ftmp;
    ftmp *= ftmp;
    ftmp *= ftmp;
    ftmp = min(ftmp, MAXX);
    old_ftmp_1 = ftmp;
    tmp_lfE[0] = 0.5f*old_ftmp + 0.5f*ftmp;

    for(i = 1; i < DCT_L_POST-1; i++)
    {
        tmp_lfE[i] = 0.333f*old_ftmp + 0.333f*old_ftmp_1;
        old_ftmp = old_ftmp_1;
        ftmp = lf_E[i+1]*m_ave+BET;
        ftmp *= ftmp;
        ftmp *= ftmp;
        ftmp *= ftmp;
        old_ftmp_1 = min(ftmp, MAXX);
        tmp_lfE[i] += 0.333f * old_ftmp_1;
    }

    ftmp = lf_E[i] * m_ave + BET;
    ftmp *= ftmp;
    ftmp *= ftmp;
    ftmp *= ftmp;
    ftmp = min(ftmp, MAXX);
    tmp_lfE[i] = 0.5f * old_ftmp + 0.5f * ftmp;

    for(i = 0; i < BIN_4KHZ; i++)
    {
        filt_lfE[i] = tmp_lfE[i] * 0.05f + 0.95f * filt_lfE[i];
    }

    for(; i < DCT_L_POST; i++)
    {
        filt_lfE[i] = tmp_lfE[i] * 0.15f + 0.85f * filt_lfE[i];
    }

    /*------------------------------------------------------------------------*
     * Reduce inter-harmonic noise with SNR based method
     * Second stage of spectral shaping modification based on the pow(x,4) energy spectrum
     *------------------------------------------------------------------------*/

    if( coder_type == AUDIO )
    {
        MAX_band = 16;
    }

    spectrum_mod_dct( DCT_buf, lf_E, lf_EO, bckr, local_min_gain, lp_gbin, music_flag2, min_band, MAX_GN, MAX_band );

    i = 0;
    if( music_flag2 >= 1 )
    {
        for(i = 0; i < BIN_1KHZ; i++)
        {
            ftmp = min(1.0f,  filt_lfE[i]);
            DCT_buf[i] *= ftmp;
        }
    }

    if( *last_nonfull_music > 40 )
    {
        max_ovf_2k = 1.25f;
        max_ovf_4k = 1.5f;
        max_ovf_6k = 1.5f;

        min_g_2k = 0.0f;
        min_g_4k = 0.0f;
        min_g_6k = 0.0f;

        if( coder_type == AUDIO )
        {
            max_ovf_2k = 1.0f;
            max_ovf_4k = 1.1f;
            max_ovf_6k = 1.25f;

            min_g_2k = 0.75f;
            min_g_4k = 0.5f;
            min_g_6k = 0.5f;

            if( core_brate  > ACELP_9k60 )
            {
                max_ovf_4k = 1.0f;
                max_ovf_6k = 1.15f;

                min_g_2k = 0.925f;
                min_g_4k = 0.825f;
                min_g_6k = 0.75f;
            }
        }
        else if( core_brate  >= ACELP_12k65 )
        {
            max_ovf_2k = 1.0f;
            max_ovf_4k = 1.25f;

            if( core_brate > ACELP_15k85 )
            {
                max_ovf_4k = 1.0f;
                max_ovf_6k = 1.25f;

                min_g_2k = 0.75f;
                min_g_4k = 0.5f;
                min_g_6k = 0.5f;
            }
        }

        for(; i < BIN_2KHZ; i++)
        {
            ftmp = min(max_ovf_2k,  filt_lfE[i]);
            ftmp = max(min_g_2k, ftmp);
            DCT_buf[i] *= ftmp;
        }

        for(; i < BIN_4KHZ; i++)
        {
            ftmp = min(max_ovf_4k,  filt_lfE[i]);
            ftmp = max(min_g_4k, ftmp);
            DCT_buf[i] *= ftmp;
        }

        if( coder_type != AUDIO || core_brate  > ACELP_8k85 )
        {
            /* Do not modify HF when coded with GSC at LR, because the spectrum is just noise */
            for(; i < DCT_L_POST; i++)
            {
                ftmp = min(max_ovf_6k,  filt_lfE[i]);
                ftmp = max(min_g_6k, ftmp);
                DCT_buf[i] *= ftmp;
            }
        }
    }
    else if( *last_nonfull_music > 25 )
    {
        /* When unsure on content type only slight clean-up allowed, no overshoot allowed */
        for(; i < DCT_L_POST; i++)
        {
            ftmp = min(1.f,  filt_lfE[i]);
            DCT_buf[i] *= ftmp;
        }
    }

    /* reconstruction of the enhanced synthesis */
    mvr2r( DCT_buf, dtc_out, DCT_L_POST );
}

/*---------------------------------------------------------------------------*
 * spectrum_mod_dct()
 *
 * spectrum enhancement according to the output of signal_type_clas()
 *---------------------------------------------------------------------------*/

static void spectrum_mod_dct(
    float data[],       /* i/o: DCT spectrum                                */
    const float lf_E[],       /* i:   per bin E for first 46 bins (without DC)    */
    float lf_EO[],      /* i/o: old per bin E for previous half frame       */
    const float noiseE[],     /* i:   per band background noise energy estimate   */
    const float minGain,      /* i:   minimum suppression gain                    */
    float lp_gbin[],    /* i/o: Smoothed suppression gain, per FFT bin      */
    const short music_flag,   /* i:   music ? 1:0                                 */
    short min_band,
    const float MAX_GN,
    const short MAX_band
)
{
    float maxNoise;
    float binE[VOIC_BINS_HR];
    float gain, minE;
    float freq;
    float slope;
    float inv_noise[MBANDS_GN_LD];
    float *pt_gbin, alpha, tmpN;
    short i, cnt;
    float *pt;
    float tmp, shift;
    const float *pt2;

    gain = 0.0f;

    /*-----------------------------------------------------------------*
     * Compute the inverse of noise
     *-----------------------------------------------------------------*/

    for (i=0; i<=MBANDS_GN_LD-1; i++)
    {
        inv_noise[i] = 1.0f / noiseE[i];
    }

    /*----------------------------------------------------------------------*
     * Perform noise reduction for 1 frames
     *----------------------------------------------------------------------*/

    for (i=0; i < VOIC_BINS_HR; i++)
    {
        binE[i] = 0.3f * lf_EO[i] + 0.7f * lf_E[i];
    }
    mvr2r(lf_E, lf_EO, VOIC_BINS_HR); /* update */

    /*----------------------------------------------------------------------*
     * Find the maximum noise in a critical band
     *----------------------------------------------------------------------*/

    maxNoise = 0.0f;
    for (i=0; i<=MBANDS_GN_LD-1; i++)
    {
        set_max(&maxNoise, noiseE[i]);
    }

    /* pointer initialization */
    pt = &data[0];

    /*-----------------------------------------------------------------*
     * Initialization for active speech frames or VAD hangover frames,
     * (exclude Clean speech)
     *-----------------------------------------------------------------*/

    if ( music_flag != 0 ) /* prevent subtraction on clean speech */
    {
        if( maxNoise <= 10.0f)
        {
            minE = .5625f;
        }
        else
        {
            minE = minGain*minGain;
        }

        pt2 = binE;
        freq = 0;

        pt_gbin = lp_gbin;
        tmp = INV_MAX_SNR;
        slope =  tmp - tmp * minE;
        shift =  MAX_SNR_SNR1 * minE - tmp;
        for (i=0; i < min_band; i++)
        {
            for (; freq <= mfreq_loc_LD[i]; freq += BIN_16kdct)
            {
                pt2++;
                pt++;
                *pt_gbin++ = 1.0f;
            }
        }

        /*-----------------------------------------------------------------*
         * Per Frequency MODE1_BIN processing
         * For highly voiced and highly pitched speech, use per bin
         * subtraction in low frequencies (maximum up to 3700 Hz,
         * first 17 critical bands)
         *-----------------------------------------------------------------*/

        for (; i < MAX_band/*MBANDS_GN_LD*/; i++)
        {

            tmp = INV_MAX_SNR_tab[i];
            slope =  tmp - tmp * minE;
            shift =  MAX_SNR_SNR1_tab[i] * minE - tmp;

            tmpN = slope * inv_noise[i];
            tmp = 0.0f;
            cnt = 0;
            while (freq <= mfreq_loc_LD[i])
            {
                gain = 1.0f;

                if (noiseE[i] >= 0.5f)            /* Do not alter if noise E very low */
                {
                    gain = tmpN **pt2 + shift; /* limits: [x,y] = {[1, minE], [MAX_SNR1, 1]}, */
                }

                pt2++;

                if (gain < minE)
                {
                    gain = minE;
                }

                if (gain > 1.0f+MAX_GN)
                {
                    gain = 1.0f+MAX_GN;
                }

                /* the gain smoothing control: stronger lp filtering for lower gains */
                alpha = 1.0f - (float)sqrt(gain);

                *pt_gbin = gain + alpha **pt_gbin;
                *pt++ *= *pt_gbin;
                cnt++;
                freq += BIN_16kdct;
                pt_gbin++;
            }
        }
    }
    else
    {
        freq = BIN_16kdct;
        pt_gbin = lp_gbin;
        for (i=0; i < MBANDS_GN_LD; i++)
        {
            for (; freq <= mfreq_loc_LD[i]; freq += BIN_16kdct)
            {
                *pt_gbin = 0.9f* *pt_gbin + 0.1f;
                pt_gbin++;
            }
        }
    }

    return;
}

/*----------------------------------------------------------------------------------*
 * analy_sp_dct()
 *
 * Spectral analysis of the current synthesized frame
 *----------------------------------------------------------------------------------*/

static void analy_sp_dct(
    const float dct_in[], /* i: input DCT spectrum                  */
    float dct_buf[],      /* i: output DCT spectrum                 */
    float *fr_bands,      /* o: energy in critical frequency bands */
    float *lf_E,          /* o: per bin E for first...             */
    float *etot           /* o: total input energy                 */
)
{
    float Bin_E[DCT_L_POST];

    *etot = 0.0f;
    mvr2r( dct_in, dct_buf, DCT_L_POST );

    /*-----------------------------------------------------------------*
     * Windowing
     * Compute spectrum
     * find energy per critical frequency band and total energy in dB
     *-----------------------------------------------------------------*/

    find_enr_dct( dct_buf, fr_bands, lf_E, etot, 0, MBANDS_GN_LD, Bin_E, BIN_16kdct );

    /* find average log total energy over both half-frames */
    *etot = 10.0f * (float)log10(*etot) - 3.0103f;

    return;
}

/*------------------------------------------------------------------------*
 * find_enr_dct)
 *
 * find input signal energy for each critical band and first 74 LF bins
 * The energy is normalized by the number of frequency bins in a channel
 *------------------------------------------------------------------------*/

void find_enr_dct(
    const float data[],    /* i  : fft result, for the format see fft_rel.c */
    float band[],    /* o  : per band energy                          */
    float *ptE,      /* o  : per bin energy  for low frequencies      */
    float *Etot,     /* o  : total energy                             */
    const short min_band,  /* i  : minimum critical band                    */
    const short max_band,  /* i  : maximum critical band                    */
    float *Bin_E,    /* o  : Per bin energy                           */
    const float bin_freq   /* i  : Number of frequency bins                 */
)
{
    short i, cnt;
    float freq, tmp;
    const float *ptR;

    ptR = &data[0];         /* pointer to first real coefficient */
    freq = 0;
    for( i=0; i < max_band; i++ )
    {
        band[i] = 0.0f;
        cnt = 0;
        while( freq <= mfreq_loc_LD[i] )
        {
            /* energy  */
            *ptE = *ptR **ptR;

            /* normalization - corresponds to FFT normalization by 2/L_FFT */
            *ptE *= 1.0f / (DCT_L_POST);

            *Bin_E = *ptE;
            Bin_E++;
            band[i] += *ptE++;
            ptR++;

            freq += bin_freq;
            cnt++;
        }

        /* normalization per frequency bin */
        band[i] /= cnt;

        if ( band[i] < E_MIN )
        {
            band[i] = E_MIN;
        }
    }

    /*-----------------------------------------------------------------*
     * Find the total energy over the input bandwidth
     *-----------------------------------------------------------------*/

    tmp = *Etot;
    for( i = min_band; i <= NB_LIMIT_BAND; i++ )
    {
        /* total channel energy */
        tmp += band[i];
    }

    *Etot = tmp;

    return;
}

/*------------------------------------------------------------------------*
 * Prep_music_postP()
 *
 * Performs the steps needed to do the music post processing
 *------------------------------------------------------------------------*/

void Prep_music_postP(
    float exc_buffer_in[],  /* i/o: excitation buffer */
    float dct_buffer_out[], /* o  : DCT output buffer */
    float filt_lfE[],       /* i/o: long term spectrum energy */
    const short last_core,        /* i  : last core  */
    const float *pitch_buf,       /* i  : current frame pitch information */
    float *LDm_enh_lp_gbin  /* o  : smoothed suppression gain, per bin FFT */
)
{
    short i;
    float *pt1;
    const float *pt2;
    short s_pit, fr_pit;

    s_pit = (short)(pitch_buf[3]);
    fr_pit = (short)((pitch_buf[3] - s_pit) * 4.0f);

    /*------------------------------------------------------------*
     * Extrapolation of the last future part and windowing
     *------------------------------------------------------------*/

    if( last_core == HQ_CORE )
    {
        set_f( filt_lfE, 1.0f, DCT_L_POST );
        set_f( LDm_enh_lp_gbin, 1.0f, VOIC_BINS_HR );
        pt1 = exc_buffer_in + OFFSET2 - 1;
        pt2 = pt1 + (short)(pitch_buf[0] + 0.5f);
        for( i = 0; i < OFFSET2; i++ )
        {
            *pt1 = *pt2;
            pt1--;
            pt2--;
        }
    }

    pt1 = exc_buffer_in + DCT_L_POST - OFFSET2;
    pred_lt4( pt1, pt1, s_pit, fr_pit, OFFSET2, inter4_2, L_INTERPOL2, PIT_UP_SAMP );

    pt2 = post_dct_wind;
    for( i = 0; i < OFFSET2; i++ )
    {
        *pt1 *= *pt2;
        pt1++;
        pt2++;
    }

    pt1 = exc_buffer_in;
    pt2--;
    for( i = 0; i < OFFSET2; i++ )
    {
        *pt1 *= *pt2;
        pt1++;
        pt2--;
    }

    edct( exc_buffer_in, dct_buffer_out, DCT_L_POST );

    return;
}

/*------------------------------------------------------------------------*
 * Post_music_postP()
 *
 * Going back from frequency to time domain from the enhanced spectrum
 * to retreive the aligned excitation and redo the synthesis
 *------------------------------------------------------------------------*/

void Post_music_postP(
    float dct_buffer_in[],  /* i/o: excitation buffer */
    float exc_buffer_out[], /* o  : DCT output buffer */
    float *exc2,            /* i/o: Current excitation to be overwriten */
    const float *mem_tmp,         /* i  : previous frame synthesis memory     */
    float *st_mem_syn2,     /* i/o: current frame synthesis memory      */
    const float *Aq,              /* i  : LPC filter coefficients             */
    float *syn              /* i/o: 12k8 synthesis                      */
)
{
    edct( dct_buffer_in, exc_buffer_out, DCT_L_POST );

    mvr2r( exc_buffer_out + OFFSET2, exc2, L_FRAME );
    mvr2r( mem_tmp, st_mem_syn2, M );

    syn_12k8( L_FRAME, Aq, exc2, syn, st_mem_syn2, 1 );

    return;
}

