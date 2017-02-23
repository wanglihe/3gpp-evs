/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <math.h>
#include "cnst.h"
#include "prot.h"

/*-----------------------------------------------------------------*
 * Local constants
 *-----------------------------------------------------------------*/

#define ALPHA           0.1f

#define COR_MIN8        0.65f
#define COR_MAX8        0.7f
#define TH_PC8          14
#define TH_EPS8         10.4f
#define TH_STA8         5.0e5f
#define K_8             0.0091f
#define C_8             0.3185f
#define ALPHA_MAX8      0.999f
#define THR_NCHAR8      1.0f

#define COR_MIN16       0.52f
#define COR_MAX16       0.85f
#define TH_PC16         4
#define TH_EPS16        1.6f
/* 10.0e5f causes problems with music - as the noise estimate starts to track the music */
#define TH_STA16        3.5e5f
#define K_16            0.0245f
#define C_16            -0.235f
#define ALPHA_MAX16     0.99f
#define THR_NCHAR16     1.0f

#define TH_PC           12
#define TH_PRED         0.8f
#define THR_SPDIV       5

#define HC_CNT          20           /* limit for harm corr count      */
#define HC_CNT_SLOW     80           /* limit for harm corr count slow */
#define BCKR_SLOW_UPDATE_SCALE 0.1f  /* step size for slow bckr update */

#define HE_LT_THR1      10.0f
#define HE_LT_THR2      30.0f
#define HE_LT_CNT       30

/*-----------------------------------------------------------------*
 * noise_est_init()
 *
 * Initialization of Noise estimator
 *-----------------------------------------------------------------*/

void noise_est_init(
    float *totalNoise,            /* o  : noise estimate over all critical bands     */
    short *first_noise_updt,      /* o  : noise update initialization flag           */
    float bckr[],                 /* o  : per band background noise energy estimate  */
    float enrO[],                 /* o  : per band old input energy                  */
    float ave_enr[],              /* o  : per band long-term average energies        */
    short *pitO,                  /* o  : open-loop pitch values from preceed. frame */
    short *aEn,                   /* o  : noise adaptation hangover counter          */
    short *st_harm_cor_cnt,       /* i/o: 1st harm correlation timer                */
    short *bg_cnt,                /* i/o: pause burst length counter                */
    float *lt_tn_track,
    float *lt_tn_dist,
    float *lt_Ellp_dist,
    float *lt_haco_ev,
    short *low_tn_track_cnt
    ,float *Etot_st_est,
    float *Etot_sq_st_est
)
{
    short i;

    for( i=0; i<NB_BANDS; i++ )
    {
        enrO[i] = E_MIN;
        bckr[i] = E_MIN;
        ave_enr[i] = E_MIN;
    }

    pitO[0] = 0;
    *totalNoise = 0.0f;
    *first_noise_updt = 0;

    *aEn = 6;
    *st_harm_cor_cnt = 0;
    *bg_cnt = 0;
    *lt_tn_track=0.20f;
    *lt_tn_dist=0.0f;
    *lt_Ellp_dist=0.0f;
    *lt_haco_ev=0.4f;
    *low_tn_track_cnt=0;

    *Etot_st_est=20.0f;
    *Etot_sq_st_est=400.0f;

    return;
}

/*-----------------------------------------------------------------*
 * noise_est_pre()
 *
 * Track energy and signal dynamics
 *-----------------------------------------------------------------*/

void noise_est_pre(
    const float Etot,           /* i  : Energy of current frame  */
    const short ini_frame,      /* i  : Frame number (init)      */
    float *Etot_l,        /* i/o: Track energy from below  */
    float *Etot_h,        /* i/o: Track energy from above  */
    float *Etot_l_lp,     /* i/o: Smoothed low energy      */
    float *Etot_last,     /* i/o: Energy of last frame     */
    float *Etot_v_h2,      /* i/o: Energy variations        */
    float *sign_dyn_lp,   /* i/o: Smoother signal dynamics */
    short harm_cor_cnt
    ,float *Etot_lp
)
{
    if ( ini_frame <= 1 )
    {
        *Etot_lp   = Etot;
        *Etot_h = Etot;
        *Etot_l = Etot;
        *Etot_l_lp = Etot;
        *Etot_last = Etot;
        *Etot_v_h2 = 0.0f;
        *sign_dyn_lp = 0.0f;
    }
    else
    {
        *Etot_lp = 0.20f*Etot + 0.80f* *Etot_lp;
        *Etot_h -= 0.04f;
        if ( Etot > *Etot_h )
        {
            *Etot_h = Etot;
        }
        *Etot_l += 0.08f;

        /* Could even be higher but it also delays first entry to DTX */
        if ( harm_cor_cnt > 50 )
        {
            if ( ini_frame < min(150,MAX_FRAME_COUNTER-1) &&
                    (Etot_h - Etot_lp) < 3.0f )
            {
                *Etot_l += min(2,(*Etot_last-*Etot_l)*0.1f);
            }

            /* Avoids large steps in short active segments */
            if ( *Etot_last - *Etot_l > HE_LT_THR2 && harm_cor_cnt > 250 )
            {
                *Etot_l += (*Etot_last-*Etot_l)*0.02f;
            }
            else if ( (*Etot_last - *Etot_l) > HE_LT_THR1 )
            {
                *Etot_l += 0.08f;
            }
        }
        if ( Etot < *Etot_l )
        {
            *Etot_l = Etot;
        }

        if ( ini_frame < 100 && *Etot_l < *Etot_l_lp )
        {
            *Etot_l_lp = 0.1f * *Etot_l + (1.0f - 0.1f) **Etot_l_lp;
        }
        else if ( ( harm_cor_cnt > HE_LT_CNT && (*Etot_last - *Etot_l) > HE_LT_THR2 ) || ( harm_cor_cnt > HE_LT_CNT && (ini_frame < 150 ) ) ||
                  ( (*Etot_l_lp - *Etot_l) > HE_LT_THR2 ) )
        {
            *Etot_l_lp = 0.03f * *Etot_l + (1.0f - 0.03f) **Etot_l_lp;
        }
        else
        {
            *Etot_l_lp = 0.02f * *Etot_l + (1.0f - 0.02f) **Etot_l_lp;
        }
        *sign_dyn_lp = 0.1f * (*Etot_h - *Etot_l) + (1.0f - 0.1f) **sign_dyn_lp;
    }

    return;
}

/*-----------------------------------------------------------------*
 * noise_est_down()
 *
 * Down-ward noise udatation routine
 * Total Noise computation, relative frame Energy computation
 * Noise energy update - here, the energy is updated only if it is
 * decreasing to improve noise suppression. Otherwise, the noise
 * update is done on noise-only frames and this decision is made in
 * nois_est() later in this file.
 *-----------------------------------------------------------------*/

void noise_est_down(
    const float fr_bands[],        /* i  : per band input energy (contains 2 vectors) */
    float bckr[],            /* i/o: per band background noise energy estimate  */
    float tmpN[],            /* o  : temporary noise update                     */
    float enr[],             /* o  : averaged energy over both subframes        */
    const short min_band,          /* i  : minimum critical band                      */
    const short max_band,          /* i  : maximum critical band                      */
    float *totalNoise,       /* o  : noise estimate over all critical bands     */
    const float Etot,              /* i  : Energy of current frame                    */
    float *Etot_last,        /* i/o: Energy of last frame                       */
    float *Etot_v_h2         /* i/o: Energy variaions of noise frames           */
)
{
    const float *pt1, *pt2;
    short i;
    float Etot_v;

    /*-----------------------------------------------------------------*
     * Estimate total noise energy
     *-----------------------------------------------------------------*/

    *totalNoise = 0.0f;
    for( i = min_band; i <= max_band; i++ )
    {
        *totalNoise += bckr[i];
    }
    *totalNoise = 10.0f * (float)log10( *totalNoise );

    /*-----------------------------------------------------------------*
     * Average energy per frame for each frequency band
     *-----------------------------------------------------------------*/

    pt1 = fr_bands;
    pt2 = fr_bands + NB_BANDS;

    for( i=0 ; i < NB_BANDS; i++ )
    {
        enr[i] = 0.5f * ( *pt1++ + *pt2++ );
    }

    /*-----------------------------------------------------------------*
     * Background noise energy update
     *-----------------------------------------------------------------*/

    for( i=0; i< NB_BANDS; i++ )
    {
        tmpN[i] = (1-ALPHA) * bckr[i] + ALPHA * enr[i];
        if( tmpN[i] < bckr[i] )
        {
            bckr[i] = tmpN[i];          /* Defend to increase noise estimate: keep as it is or decrease  */
        }
    }

    /*------------------------------------------------------------------*
     * Energy variation update
     *------------------------------------------------------------------*/

    Etot_v = (float) fabs(*Etot_last - Etot);
    *Etot_v_h2 = (1.0f-0.02f) **Etot_v_h2 + 0.02f * min(3.0f,Etot_v);
    if (*Etot_v_h2 < 0.1f)
    {
        *Etot_v_h2 = 0.1f;
    }

    return;
}

/*-----------------------------------------------------------------*
 * noise_est()
 *
 * Noise energy estimation (noise energy is updated in case of noise-only frame)
 *-----------------------------------------------------------------*/

void noise_est(
    Encoder_State *st,                    /* i/o: encoder state structure                                */
    const float tmpN[],                 /* i  : temporary noise update                                 */
    const short *pitch,                 /* i  : open-loop pitch values for each half-frame             */
    const float *voicing,               /* i  : normalized correlation for all half-frames             */
    const float *epsP,                  /* i  : LP prediction error energies                           */
    const float Etot,                   /* i  : total channel E                                        */
    const float relE,                   /* i  : relative frame energy                                  */
    const float corr_shift,             /* i  : normalized correlation correction                      */
    const float enr[],                  /* i  : averaged energy over both subframes                    */
    float fr_bands[],             /* i  : spectrum per critical bands of the current frame       */
    float *cor_map_sum,           /* o  : sum of correlation map from mult-harm analysis         */
    float *sp_div,                /* o  : soectral diversity feature                             */
    float *non_staX,              /* o  : non-stationarity for sp/mus classifier                 */
    short *loc_harm,              /* o  : multi-harmonicity flag for UV classifier               */
    const float *lf_E,                  /* i  : per bin energy  for low frequencies                    */
    short *st_harm_cor_cnt,       /* i/o : 1st harm correlation timer                            */
    const float Etot_l_lp,              /* i   : Smoothed low energy                                   */
    float *sp_floor               /* o  : noise floor estimate                                   */
)
{
    short i, tmp_pc, pc, spec_div, noise_char;
    float alpha, th_eps, th_sta, non_sta, cor_min, cor_max;
    float non_sta2, alpha2, sum_num, sum_den, *pt1, *pt2, ftemp, ftemp2, nchar_thr;
    float updt_step,log_enr;

    short aE_bgd,sd1_bgd,bg_bgd2;
    short tn_ini;

    float epsP_0_2,epsP_0_2_ad,epsP_0_2_ad_lp_max;
    float epsP_2_16,epsP_2_16_dlp,epsP_2_16_dlp_max;
    short PAU,BG_1,NEW_POS_BG;

    float haco_ev_max;
    float Etot_l_lp_thr;
    float comb_ahc_epsP,comb_hcm_epsP;

    short enr_bgd,cns_bgd,lp_bgd,ns_mask;
    short lt_haco_mask, bg_haco_mask;
    short SD_1,bg_bgd3,PD_1,PD_2,PD_3,PD_4,PD_5;

    float tmp_enr,tmp_ave,tmp_ave2;
    float non_staB;

    float lim_Etot;

    /*-----------------------------------------------------------------*
     * Initialization
     *-----------------------------------------------------------------*/

    st->ener_RAT = 10.0f * (float)log10( mean(lf_E, 8) );
    st->ener_RAT /= (Etot + 0.01f);

    if( st->ener_RAT < 0.0f )
    {
        st->ener_RAT = 0.0f;
    }

    if( st->ener_RAT > 1.0 )
    {
        st->ener_RAT = 1.0f;
    }

    /*-----------------------------------------------------------------*
     * Set the threshold for eps & non_sta based on input sampling rate
     * The reason is that in case of 8kHz sampling input, there is nothing
     * between 4kHz-6.4kHz. In noisy conditions, this makes a fast
     * transition even in noise-only parts, hence producing a "higher
     * order" spectral envelope => the epsP ratio is much less effective.
     *-----------------------------------------------------------------*/

    if( st->input_bwidth != NB )
    {
        /* WB input */
        th_eps = TH_EPS16;
        th_sta = TH_STA16;
        cor_min = COR_MIN16;
        cor_max = COR_MAX16;
    }
    else
    {
        /* NB input */
        th_eps = TH_EPS8;
        th_sta = TH_STA8;
        cor_min = COR_MIN8;
        cor_max = COR_MAX8;
    }

    /*-----------------------------------------------------------------*
     * Estimation of pitch stationarity
     *-----------------------------------------------------------------*/

    pc = (short)( abs(pitch[0] - st->pitO) + abs(pitch[1] - pitch[0]) );

    if ( ( (voicing[0] + voicing[1] + voicing[2]) / 3.0f + corr_shift ) < cor_min )
    {
        /* low correlation -> probably inactive signal */
        tmp_pc = TH_PC;
    }
    else
    {
        tmp_pc = pc;
    }

    st->pitO = pitch[1];

    /*-----------------------------------------------------------------*
     * Multi-harmonic analysis
     *-----------------------------------------------------------------*/

    *loc_harm = multi_harm( st->Bin_E, st->old_S, st->cor_map, &st->multi_harm_limit, st->total_brate,
                            st->bwidth, &st->cor_strong_limit, &st->mean_avr_dyn, &st->last_sw_dyn, cor_map_sum, sp_floor );

    /*-----------------------------------------------------------------*
     * Detection of frames with non-stationary spectral content
     *-----------------------------------------------------------------*/

    /* weighted sum of spectral changes per critical bands */
    sum_num = 0;
    sum_den = 0;

    pt1 = fr_bands + 10;
    pt2 = st->fr_bands2 + 10;
    for (i = 10; i <= st->max_band; i++)
    {
        if (*pt1 > *pt2)
        {
            sum_num += *pt1 **pt1 / *pt2;
            sum_den += *pt1;
        }
        else
        {
            sum_num += *pt2 **pt2 / *pt1;
            sum_den += *pt2;
        }

        pt1++;
        pt2++;
    }

    /* calculation of spectral diversity */
    if (sum_num > THR_SPDIV * sum_den)
    {
        spec_div = 1;
    }
    else
    {
        spec_div = 0;
    }

    *sp_div = sum_num / (sum_den + 1e-5f);

    /*-----------------------------------------------------------------*
     * Detection of frames with high energy content in high frequencies
     *-----------------------------------------------------------------*/

    /* calculation of energy in first 10 critical bands */
    ftemp = sum_f( &fr_bands[st->min_band], 10 - st->min_band );

    /* calculation of energy in the rest of bands */
    ftemp2 = sum_f( &fr_bands[10], st->max_band - 10 + 1 );

    if ( ftemp < 1e2 || ftemp2 < 1e2 )
    {
        ftemp2 = 0;
    }
    else
    {
        ftemp2 /= ftemp;
    }

    if ( ftemp2 > 10 )
    {
        ftemp2 = 10;
    }

    /* update LT value of the final parameter */
    st->noise_char = M_ALPHA * st->noise_char + (1-M_ALPHA) * ftemp2;

    if( st->input_bwidth == NB )
    {
        nchar_thr = THR_NCHAR8;
    }
    else
    {
        nchar_thr = THR_NCHAR16;
    }

    if( st->noise_char > nchar_thr )
    {
        noise_char = 1;
    }
    else
    {
        noise_char = 0;
    }

    /* save the 2 last spectra per crit. bands for the future */
    mvr2r( st->fr_bands1, st->fr_bands2, NB_BANDS );
    mvr2r( fr_bands+NB_BANDS, st->fr_bands1, NB_BANDS );

    /*-----------------------------------------------------------------*
     * Non-stationarity estimation for each band (handicap high E frames in average computing)
     *-----------------------------------------------------------------*/

    /* set averaging factor */
    ftemp = relE;
    if( ftemp < 0.0f )
    {
        ftemp = 0.0f;
    }

    alpha =  0.064f * ftemp + 0.75f;

    if( alpha > 0.999f )
    {
        alpha = 0.999f;
    }

    /* during significant attacks, replace LT energy by the */
    /* current energy - this will cause non_sta2 failures to occur in */
    /* different frames than non_sta failures */
    alpha2 = alpha;
    if ( spec_div > 0)
    {
        alpha2 = 0.0f;
    }

    /* calculate non-stationarity */
    non_sta = 1.0f;
    non_sta2 = 1.0f;
    *non_staX = 0.0f;
    non_staB =0.0f;
    for( i = st->min_band; i <= st->max_band; i++ )
    {
        /* + 1.0f added to reduce sencitivity to non stationarity in low energies  */
        tmp_enr = enr[i] + 1.0f;
        if( non_sta <= th_sta )                               /* Just to limit the saturation */
        {
            tmp_ave =  st->ave_enr[i] + 1.0f;
            if( tmp_enr > tmp_ave )
            {
                non_sta = non_sta * (tmp_enr / tmp_ave );    /* non_stationarity measure  */
            }
            else
            {
                non_sta = non_sta * (tmp_ave / tmp_enr );    /* non_stationarity measure  */
            }
        }
        st->ave_enr[i] = alpha * st->ave_enr[i] + (1-alpha) * enr[i]; /* update long-term average */

        /* calculation of another non-stationarity measure (following attacks) */
        if( non_sta2 <= th_sta )
        {
            tmp_ave2 =  st->ave_enr2[i] + 1.0f;
            if( tmp_enr > tmp_ave2 )
            {
                non_sta2 = non_sta2 * ( tmp_enr / tmp_ave2 );
            }
            else
            {
                non_sta2 = non_sta2 * (tmp_ave2 / tmp_enr );
            }
        }

        st->ave_enr2[i] = alpha2 * st->ave_enr2[i] + (1-alpha2) * enr[i];

        /* calculate non-stationarity feature for speech/music classifier */
        if( i >= START_BAND_SPMUS && i < NB_BANDS_SPMUS+START_BAND_SPMUS )
        {
            log_enr = (float)log(enr[i]);
            if( log_enr > st->past_log_enr[i-START_BAND_SPMUS] )
            {
                *non_staX += log_enr - st->past_log_enr[i-START_BAND_SPMUS];
            }
            else
            {
                *non_staX += st->past_log_enr[i-START_BAND_SPMUS] - log_enr;
            }

            st->past_log_enr[i-START_BAND_SPMUS] = log_enr;
        }
        /* calculate non-stationarity feature relative background */
        if (st->ini_frame < 100)
        {
            /* During init don't include updates */
            if ( i >= 2 && i <= 16 )
            {
                non_staB += (float)fabs(log(enr[i] + 1.0f) -
                                        log(E_MIN + 1.0f));
            }
        }
        else
        {
            /* After init compare with background estimate */
            if ( i >= 2 && i <= 16 )
            {
                non_staB += (float)fabs(log(enr[i] + 1.0f) -
                                        log(st->bckr[i] + 1.0f));
            }
        }
        if (non_staB >= 128)
        {
            non_staB = 32767.0/256.0f;
        }
    }
    if ( Etot < -5.0f )
    {
        non_sta  = 1.0f;
        non_sta2 = 1.0f;
    }

    lim_Etot = max(20.0f,Etot);

    if ( st->ini_frame < 150 )
    {
        /* Allow use of quicker filter during init - if needed */
        st->Etot_st_est = 0.25f * lim_Etot + (1.0f-0.25F) * st->Etot_st_est;
        st->Etot_sq_st_est = 0.25f * lim_Etot * lim_Etot + (1.0f-0.25f) * st->Etot_sq_st_est;
    }
    else
    {
        st->Etot_st_est = 0.25f * lim_Etot + (1.0f-0.25f) * st->Etot_st_est;
        st->Etot_sq_st_est = 0.25f * lim_Etot * lim_Etot + (1.0f-0.25f) * st->Etot_sq_st_est;
    }
    /*-----------------------------------------------------------------*
     * Count frames since last correlation or harmonic event
     *-----------------------------------------------------------------*/
    if ( Etot > 0 &&
            (*loc_harm > 0 || 0.5f * (voicing[0] + voicing[1]) > 0.85f)  )

    {
        st->harm_cor_cnt = 0;
    }
    else
    {
        st->harm_cor_cnt += 1;
    }

    if ( st->harm_cor_cnt > 1 && ( ( Etot < 15.0f )  ||
                                   ( st->ini_frame>10 &&
                                     ( Etot-st->Etot_lp ) > 7.0f )
                                 )
       )
    {
        st->harm_cor_cnt = 1;
    }
    if ( st->harm_cor_cnt > 1 &&
            Etot > 30.0f &&
            (st->Etot_sq_st_est - st->Etot_st_est*st->Etot_st_est) > 8.0f)
    {
        st->harm_cor_cnt = max(1,(short) round_f( (float) st->harm_cor_cnt / 4.0f )) ;
    }

    /*-----------------------------------------------------------------*
     * Energy-based pause-length counter
     *-----------------------------------------------------------------*/

    if ( st->bg_cnt >= 0 && (Etot - st->Etot_l_lp) > 5 )
    {
        /* probably speech burst */
        st->bg_cnt = -1;
    }
    else
    {
        if ( st->bg_cnt == -1 && (Etot - st->Etot_l_lp) < 5 )
        {
            /* probably start of speech pause */
            st->bg_cnt = 0;
        }
    }

    if (st->bg_cnt >= 0)
    {
        st->bg_cnt += 1;
    }

    /*-----------------------------------------------------------------*
     * Linear predition efficiency 0 to 2 order
     *-----------------------------------------------------------------*/

    epsP_0_2 = max(0 , min(8, epsP[0] / epsP[2]));
    st->epsP_0_2_lp = 0.15f * epsP_0_2 + (1.0f-0.15f) * st->epsP_0_2_lp;
    epsP_0_2_ad = (float) fabs(epsP_0_2 - st->epsP_0_2_lp );
    if (epsP_0_2_ad < st->epsP_0_2_ad_lp)
    {
        st->epsP_0_2_ad_lp = 0.1f * epsP_0_2_ad + (1.0f - 0.1f) * st->epsP_0_2_ad_lp;
    }
    else
    {
        st->epsP_0_2_ad_lp = 0.2f * epsP_0_2_ad + (1.0f - 0.2f) * st->epsP_0_2_ad_lp;
    }
    epsP_0_2_ad_lp_max = max(epsP_0_2_ad,st->epsP_0_2_ad_lp);

    /*-----------------------------------------------------------------*
     * Linear predition efficiency 2 to 16 order
     *-----------------------------------------------------------------*/

    epsP_2_16 = max(0 , min(8, epsP[2] / epsP[16]));
    if (epsP_2_16 > st->epsP_2_16_lp)
    {
        st->epsP_2_16_lp = 0.2f * epsP_2_16 + (1.0f-0.2f) * st->epsP_2_16_lp;
    }
    else
    {
        st->epsP_2_16_lp = 0.03f * epsP_2_16 + (1.0f-0.03f) * st->epsP_2_16_lp;
    }
    st->epsP_2_16_lp2 = 0.02f * epsP_2_16 + (1.0f-0.02f) * st->epsP_2_16_lp2;

    epsP_2_16_dlp = st->epsP_2_16_lp-st->epsP_2_16_lp2;

    if (epsP_2_16_dlp < st->epsP_2_16_dlp_lp2 )
    {
        st->epsP_2_16_dlp_lp2 = 0.02f * epsP_2_16_dlp + (1.0f-0.02f) * st->epsP_2_16_dlp_lp2;
    }
    else
    {
        st->epsP_2_16_dlp_lp2 = 0.05f * epsP_2_16_dlp + (1.0f-0.05f) * st->epsP_2_16_dlp_lp2;
    }

    epsP_2_16_dlp_max  = max(epsP_2_16_dlp,st->epsP_2_16_dlp_lp2);

    /*-----------------------------------------------------------------*
     * long term extensions of frame features
     *-----------------------------------------------------------------*/

    st->lt_tn_track = 0.03f* (Etot - st->totalNoise < 10) + 0.97f*st->lt_tn_track;
    st->lt_tn_dist = 0.03f* (Etot - st->totalNoise) + 0.97f*st->lt_tn_dist;
    st->lt_Ellp_dist = 0.03f* (Etot - st->Etot_l_lp) + 0.97f*st->lt_Ellp_dist;

    if (st->harm_cor_cnt == 0)
    {
        st->lt_haco_ev = 0.03f + 0.97f*st->lt_haco_ev;
    }
    else
    {
        st->lt_haco_ev = 0.99f*st->lt_haco_ev;
    }

    if (st->lt_tn_track<0.05f)
    {
        st->low_tn_track_cnt++;
    }
    else
    {
        st->low_tn_track_cnt=0;
    }


    /* update of the long-term non-stationarity measure (between 0 and 1) */
    if ( (non_sta > th_sta) || (*loc_harm > 0) )
    {
        st->act_pred = M_GAMMA * st->act_pred + (1-M_GAMMA) * 1;
    }
    else
    {
        st->act_pred = M_GAMMA * st->act_pred + (1-M_GAMMA) * 0;
    }

    /*-----------------------------------------------------------------*
     * Increment/decrement counter for enabling background noise update
     *-----------------------------------------------------------------*/
    if( ( (*st_harm_cor_cnt < 3*HC_CNT_SLOW ) && ( ( non_sta > th_sta ) ||
            ( tmp_pc < TH_PC ) ||
            ( noise_char > 0) )
        ) ||
            ( (st->ini_frame > 150) && (Etot - Etot_l_lp) > 10 ) ||
            ( 0.5f * (voicing[0]+voicing[1]) > cor_max ) ||
            ( epsP[2] / epsP[16] > th_eps ) ||
            ( *loc_harm > 0) ||
            ((st->act_pred > 0.8f) && (non_sta2 > th_sta))
      )

    {
        /* active signal present - increment counter */
        st->aEn = st->aEn + 2;
    }
    else
    {
        /* background noise present - decrement counter */
        st->aEn = st->aEn - 1;
    }

    if( st->aEn > 6 )
    {
        st->aEn = 6;
    }
    else if ( st->aEn < 0 )
    {
        st->aEn = 0;
    }

    /*--------------------------------------------------------------*
     * Background noise update
     * (if this is not noise-only frame, bckr has been already updated downwards in nois_est_down())
     *--------------------------------------------------------------*/
    /* Additional detectors */

    comb_ahc_epsP = max(max(st->act_pred,st->lt_haco_ev),epsP_2_16_dlp);
    comb_hcm_epsP = max(max(st->lt_haco_ev,epsP_2_16_dlp_max),epsP_0_2_ad_lp_max);

    haco_ev_max = max(*st_harm_cor_cnt==0,st->lt_haco_ev);
    Etot_l_lp_thr = st->Etot_l_lp + (1.5f + 1.5f * (st->Etot_lp<50.0f))*st->Etot_v_h2;


    enr_bgd = Etot < Etot_l_lp_thr;
    cns_bgd = (epsP_0_2 > 7.95f) && (non_sta< 1e3f);
    lp_bgd  = epsP_2_16_dlp_max < 0.10f;
    ns_mask = non_sta < 1e5f;
    lt_haco_mask = st->lt_haco_ev < 0.5f;
    bg_haco_mask = haco_ev_max < 0.4f;

    SD_1 = ( (epsP_0_2_ad > 0.5f) && (epsP_0_2 > 7.95f) );

    bg_bgd3 = enr_bgd || ( ( cns_bgd || lp_bgd ) && ns_mask && lt_haco_mask && SD_1==0 );

    PD_1 = (epsP_2_16_dlp_max < 0.10f ) ;
    PD_2 = (epsP_0_2_ad_lp_max < 0.10f ) ;
    PD_3 = (comb_ahc_epsP < 0.85f );
    PD_4 = comb_ahc_epsP < 0.15f;
    PD_5 =  comb_hcm_epsP < 0.30f;

    BG_1 = ( (SD_1==0) || (Etot < Etot_l_lp_thr) )  && bg_haco_mask && (st->act_pred < 0.85f) && (st->Etot_lp < 50.0f);

    PAU = (st->aEn==0) || ( (Etot < 55.0f) && (SD_1==0) && ( ( PD_3 && (PD_1 || PD_2 ) ) ||  ( PD_4 || PD_5 ) ) );

    NEW_POS_BG = (PAU | BG_1) & bg_bgd3;

    /* Original silence detector works in most cases */
    aE_bgd = st->aEn == 0;

    /* When the signal dynamics is high and the energy is close to the background estimate */
    sd1_bgd = (st->sign_dyn_lp > 15) && (Etot - st->Etot_l_lp ) < 2*st->Etot_v_h2 && st->harm_cor_cnt > 20;

    /* init conditions steadily dropping act_pred and/or lt_haco_ev */
    tn_ini = st->ini_frame < 150 && st->harm_cor_cnt > 5 &&
             (Etot-st->Etot_lp) < 7 &&
             ( (st->act_pred < 0.59f && st->lt_haco_ev <0.23f ) ||
               st->act_pred < 0.38f ||
               st->lt_haco_ev < 0.15f ||
               non_staB < 50.0f ||
               aE_bgd
               || ( Etot < 42.0f
                    && st->harm_cor_cnt > 10
                    && st->lt_haco_ev <0.35f
                    && st->act_pred <0.8f )
             );
    /* Energy close to the background estimate serves as a mask for other background detectors */
    bg_bgd2 = Etot < Etot_l_lp_thr || tn_ini ;

    updt_step=0.0f;
    if (( bg_bgd2 && ( aE_bgd || sd1_bgd || st->lt_tn_track >0.90f || NEW_POS_BG ) ) ||
            tn_ini )
    {
        if( ( ( st->act_pred < 0.85f )  &&
                aE_bgd &&
                ( st->lt_Ellp_dist < 10  || sd1_bgd ) && st->lt_tn_dist<40 &&
                ( ( Etot - st->totalNoise ) < 10.0f ) ) ||
                ( st->first_noise_updt == 0 && st->harm_cor_cnt > 80 && aE_bgd && st->lt_aEn_zero > 0.5f ) ||
                ( tn_ini && ( aE_bgd || non_staB < 10.0 || st->harm_cor_cnt > 80 ) )
          )
        {
            updt_step=1.0f;
            st->first_noise_updt = 1;
            for( i=0; i< NB_BANDS; i++ )
            {
                st->bckr[i] = tmpN[i];
            }
        }
        else if ( ( ( st->act_pred < 0.80f ) && ( aE_bgd  || PAU )  &&  st->lt_haco_ev < 0.10f )  ||
                  ( ( st->act_pred < 0.70f ) && ( aE_bgd || non_staB < 17.0f ) && PAU &&  st->lt_haco_ev < 0.15f ) ||
                  ( st->harm_cor_cnt > 80 && st->totalNoise > 5.0f && Etot < max(1.0f,Etot_l_lp + 1.5f* st->Etot_v_h2) ) ||
                  ( st->harm_cor_cnt > 50 && st->first_noise_updt > 30 && aE_bgd  && st->lt_aEn_zero>0.5f ) ||
                  tn_ini
                )
        {
            updt_step=0.1f;
            if ( !aE_bgd &&
                    st->harm_cor_cnt < 50 &&
                    ( st->act_pred > 0.6f ||
                      ( !tn_ini && Etot_l_lp - st->totalNoise < 10.0f && non_staB > 8.0f ) ) )
            {
                updt_step=0.01f;
            }
            if (updt_step > 0.0f )
            {
                st->first_noise_updt = 1;
                for( i=0; i< NB_BANDS; i++ )
                {
                    st->bckr[i] = st->bckr[i] + updt_step * (tmpN[i]-st->bckr[i]);
                }
            }
        }
        else if (aE_bgd || st->harm_cor_cnt > 100 )
        {
            (  st->first_noise_updt) += 1;
        }
    }
    else
    {
        /* If in music lower bckr to drop further */
        if ( st->low_tn_track_cnt > 300 && st->lt_haco_ev >0.9f && st->totalNoise > 0.0f)
        {
            updt_step=-0.02f;
            for( i=0; i< NB_BANDS; i++ )
            {
                if (st->bckr[i] > 2*E_MIN)
                {
                    st->bckr[i] = 0.98f*st->bckr[i];
                }
            }
        }
    }
    st->lt_aEn_zero = 0.2f * (st->aEn==0) + (1-0.2f)*st->lt_aEn_zero;
    return ;
}
