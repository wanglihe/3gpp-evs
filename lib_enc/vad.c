/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_enc.h"


/*-----------------------------------------------------------------*
 * Local constants
 *-----------------------------------------------------------------*/

#define HANGOVER_LONG             10         /* Hangover for CNG */
#define HANGOVER_LONG_HE          20         /* Hangover of CNG */
#define HANGOVER_LONG_MUSIC       20         /* Hangover of CNG */
#define HANGOVER_LONG_NB          8          /* Hangover for CNG */
#define ACTIVE_FRAMES             3          /* Number of consecutive active SPEECH frames necessary to trigger HO */

#define TH16_2                    35.0f      /* long-term SNR that separates the curves for clean speech and noisy speech */
#define TH8_1                     20.0f      /* long-term SNR that separates the curves for clean speech and noisy speech     */
#define TH16_2_NFLAG              35.0f
#define TH8_1_NFLAG               35.0f


#define SNR_OUTLIER_WGHT_1        1.00f
#define SNR_OUTLIER_WGHT_2        1.01f
#define SNR_OUTLIER_WGHT_3        1.02f
#define OUTLIER_THR_1             10.0f
#define OUTLIER_THR_2             6.0f
#define MAX_SNR_OUTLIER_IND       17
#define MAX_SNR_OUTLIER_1         10.0f
#define MAX_SNR_OUTLIER_2         25.0f
#define MAX_SNR_OUTLIER_3         50.0f

/*---------------------------------------------------------------------*
 * wb_vad_init()
 *
 * VAD initializations
 *---------------------------------------------------------------------*/

void wb_vad_init(
    short *nb_active_frames,      /* o  : nb of consecutive active speech frames                    */
    short *hangover_cnt,
    float *lp_speech,             /* o  : long-term active speech level                             */
    short *nb_active_frames_he,   /* o  : nb of consecutive active speech frames                    */
    short *hangover_cnt_he,
    float *bcg_flux,              /* o  : background noise fluctuation                              */
    short *soft_hangover,         /* o  : soft hangover counter                                     */
    short *voiced_burst,          /* o  : consecutive voiced speech counter                         */
    short *bcg_flux_init,         /* o  : initialization period for noise fluctuation estimation    */
    short *nb_active_frames_he1,  /* o  : nb of consecutive active speech frames 1                  */
    short *hangover_cnt_he1,
    long  *vad_flag_reg_H,
    long  *vad_flag_reg_L,
    long  *vad_prim_reg,
    short *vad_flag_cnt_50,
    short *vad_prim_cnt_16,
    short *hangover_cnt_dtx,
    short *flag_noisy_speech_snr
    ,short *hangover_cnt_music     /* o  : counter of VAD DTX Music hangover frames                  */
)
{
    *hangover_cnt = 0;                    /* Hangover counter initialized to 0                      */
    *nb_active_frames = ACTIVE_FRAMES;    /* The counter of SPEECH frames necessary to trigger HO   */
    /* is set to max (-> start with hangover)                 */
    *lp_speech = 45.0f;                   /* Initialize the long-term active speech level in dB     */
    *flag_noisy_speech_snr=0;

    *vad_flag_reg_H = 0L;
    *vad_flag_reg_L = 0L;
    *vad_prim_reg = 0L;
    *vad_flag_cnt_50 = 0;
    *vad_prim_cnt_16 = 0;

    /* By default one should not start with a hangover */
    *hangover_cnt_dtx = HANGOVER_LONG;             /* hangover for DTX                              */
    *hangover_cnt_music = HANGOVER_LONG_MUSIC;     /* hangover for DTX                              */

    *hangover_cnt_he = 0;                 /* Hangover counter initialized to 0                      */
    *nb_active_frames_he = ACTIVE_FRAMES; /* The counter of SPEECH frames necessary to trigger HO   */
    *bcg_flux = 70;
    *soft_hangover = 0;
    *voiced_burst = 0;
    *bcg_flux_init = 50;
    *nb_active_frames_he1 = ACTIVE_FRAMES;
    *hangover_cnt_he1 = 0;

    return;
}

/*-----------------------------------------------------------------*
 * sing_thr_snr_acc()
 *
 * accumulate snr_sum with significance thresholds
 *-----------------------------------------------------------------*/

static void sign_thr_snr_acc(
    float *snr_sum,
    float snr,
    float sign_thr,
    float min_snr
)
{
    if( snr >= sign_thr )
    {
        *snr_sum = *snr_sum + snr;
    }
    else
    {
        *snr_sum = *snr_sum + min_snr;
    }

    return;
}

/*-----------------------------------------------------------------*
 * dtx_hangover_addition()
 *
 * accumulate snr_sum with significance thresholds
 *-----------------------------------------------------------------*/

short dtx_hangover_addition(
    Encoder_State *st,                /* i/o: encoder state structure */
    const short localVAD,
    const short vad_flag,
    const float lp_snr,
    const short cldfb_subtraction,
    short *vad_hover_flag
)
{
    short hangover_short_dtx, flag_dtx;
    short ho_limit_clean;

    flag_dtx = 0;

    /* Determine initial hangover length */
    hangover_short_dtx = 2; /* was 1 */
    if ( ( lp_snr < 16.0f && st->input_bwidth != NB ) ||
            st->prim_act_he > 0.95f )
    {
        hangover_short_dtx = 3; /* was 2 */
    }

    /* Adjust hangover according to activity history */
    if (st->vad_prim_cnt_16 > 12 ) /* 12 requires roughly > 80% primary activity */
    {
        hangover_short_dtx = hangover_short_dtx + 2;
    }

    if (st->vad_flag_cnt_50 >40  ) /* 40 requires roughtly > 80% flag activity */
    {
        hangover_short_dtx = hangover_short_dtx + 5;
    }

    /* Keep hangover_short lower than maximum hangover count */
    if (hangover_short_dtx >  HANGOVER_LONG-1)
    {
        hangover_short_dtx = HANGOVER_LONG-1;
    }

    /* Only allow short HO if not sufficient active frames */
    ho_limit_clean = 3;
    if (st->core == AMR_WB_CORE )
    {
        ho_limit_clean = 2;
    }

    if ( st->input_bwidth != NB && st->core != AMR_WB_CORE && lp_snr > 25.0f )
    {
        ho_limit_clean = 2;
    }

    if ( ho_limit_clean != 0 )
    {
        if ( (hangover_short_dtx > ho_limit_clean) && ( ( st->vad_prim_cnt_16 < 7 ) || ( lp_snr > 16 && st->prim_act_he < 0.85 ) ) )
        {
            hangover_short_dtx = ho_limit_clean;
        }
    }


    /* hangover adjustment from combined FFT + CLDFBVAD */
    if (st->core != AMR_WB_CORE)
    {
        hangover_short_dtx = hangover_short_dtx - cldfb_subtraction;
        if( hangover_short_dtx < 0 )
        {
            hangover_short_dtx = 0;
        }
    }
    if ( vad_flag == 1 ) /* Speech present */
    {
        flag_dtx = 1;

        /* Add hangover after sufficient # of active frames or sufficient activity during last second */
        if (st->nb_active_frames >= ACTIVE_FRAMES || st->vad_flag_cnt_50 >45 ) /* 45 requires roughtly > 90% flag activity */
        {
            st->hangover_cnt_dtx = 0;
        }

        /* inside HO period */
        if( st->hangover_cnt_dtx < HANGOVER_LONG && st->hangover_cnt_dtx != 0 )
        {
            st->hangover_cnt_dtx++;
        }

        st->hangover_terminate_flag = 0;

        /* Music hangover when music detected */
        if( st->prim_act_he > 0.98f && st->Etot_lp > 40 && st->vad_prim_cnt_16 > 14 && st->vad_flag_cnt_50 > 48 )
        {
            st->hangover_cnt_music = 0;
        }

        /* inside music HO period */
        if (st->hangover_cnt_music < HANGOVER_LONG_MUSIC && st->hangover_cnt_music != 0 )
        {
            st->hangover_cnt_music++;
        }
    }
    else
    {
        /* Reset the counter of speech frames necessary to start hangover algorithm */
        if(st->hangover_cnt_dtx < HANGOVER_LONG )    /* inside HO period */
        {
            st->hangover_cnt_dtx++;
        }
        if(st->hangover_cnt_music < HANGOVER_LONG_MUSIC )    /* inside music HO period */
        {
            st->hangover_cnt_music++;
        }
        /* fast terminate DTX hangover if st->hangover_terminate_flag is set */
        if ( st->hangover_terminate_flag == 1 )
        {
            st->hangover_cnt = HANGOVER_LONG;
            st->hangover_cnt_dtx = HANGOVER_LONG;
            st->hangover_terminate_flag = 0;
            /* only shorten music hangover when low energy frames */
            if ( st->Etot_lp < 20.0f )
            {
                st->hangover_cnt_music = HANGOVER_LONG_MUSIC;
            }
        }

        if( st->hangover_cnt_dtx <= hangover_short_dtx )  /* "hard" hangover  */
        {
            flag_dtx = 1;
        }

        if( st->hangover_cnt_music <= 15 )  /* "hard" hangover  */
        {
            flag_dtx = 1;
        }

    }


    if ( flag_dtx != 0 && localVAD == 0 )
    {
        *vad_hover_flag = 1;
    }

    return flag_dtx ;
}


/*-----------------------------------------------------------------*
 * wb_vad()
 *
 * Voice Activity Detector
 *-----------------------------------------------------------------*/

short wb_vad(
    Encoder_State *st,                    /* i/o: encoder state structure                    */
    const float fr_bands[],             /* i  : per band input energy (contains 2 vectors) */
    short *localVAD,
    short *noisy_speech_HO,      /* o  : SC-VBR noisy speech HO flag                */
    short *clean_speech_HO,       /* o  : SC-VBR clean speech HO flag                */
    short *NB_speech_HO,          /* o  : SC-VBR NB speech HO flag                   */
    float *snr_sum_he,            /* o  : Output snr_sum as weighted spectral measure */
    short *localVAD_HE_SAD,       /* o  : HE_SAD decision without hangovers          */
    short *flag_noisy_speech_snr  /* o:  */
)
{
    short i, j, flag, hangover_short;
    float snr[NB_BANDS], snr_sum, thr1, thr2, lp_snr, nk, nc, th_clean;
    const float *pt1, *pt2, *pt3;
    float min_snr, sign_thr;
    float fr_enr;
    float ftmp, ftmp1;
    float mssnr = 0;
    float snr_sumt;
    float vad_thr;
    short hangover_hd;
    short snr_idx;
    float delta1, delta2, delta3;
    short flag_he1;
    float mssnr_hov;
    short stmp;
    float msnr;
    float snr_outlier;
    short snr_outlier_index;
    float accum_ener_L, accum_ener_H;
    float delta4;
    float snr18=1.0f, snr19=1.0f;
    short nb_sig_snr;
    float nv;
    float snr_sum_HE_SAD;
    float sign_thr_HE_SAD,min_snr_HE_SAD;
    float nv_ofs;
    float thr1_ol;
    float snr_sum_ol;
    snr_outlier = 0;
    snr_outlier_index = 0;
    accum_ener_L = 0;
    accum_ener_H = 0;

    if( st->input_bwidth == NB )
    {
        st->min_band = 1;
        st->max_band = 16;
    }
    else
    {
        st->min_band = 0;
        st->max_band = 19;
    }

    /*---------------------------------------------------------------------*
     * set SNR thresholds depending on the input rate
     *---------------------------------------------------------------------*/

    if(st->max_band == 19 )  /* WB input */
    {
        nk = 0.1f;
        nc = 16.1f;
        nv = 2.05f;
        nv_ofs = 1.65f;
        th_clean = TH16_2;
        if ( st->input_bwidth == WB )
        {
            sign_thr = 1.3f;
            min_snr  = 0.8f;
        }
        else
        {
            sign_thr = 1.75f;
            min_snr  = 0.25f;
        }
        sign_thr_HE_SAD = 2.5f;
        min_snr_HE_SAD  = 0.2f;
    }
    else                    /* NB input */
    {
        nk = 0.10f;
        nc = 16.0f;
        nv = 4.00f;         /* Was 4.5f but trunkated to 4.00 used when converted to short */
        nv_ofs = 1.15f;
        th_clean = TH8_1;
        sign_thr = 1.75f;
        min_snr  = 0.25f;

        sign_thr_HE_SAD = 2.65f;
        min_snr_HE_SAD  = 0.05f;
    }

    hangover_short = 0;


    if( st->Opt_SC_VBR )
    {
        *noisy_speech_HO = 0;
        *clean_speech_HO = 0;
        *NB_speech_HO = 0;
    }

    /*---------------------------------------------------------------------*
     * compute SNR for each band & total
     *---------------------------------------------------------------------*/

    pt1 = fr_bands;
    pt2 = fr_bands + NB_BANDS;
    snr_sum = 0.0f;
    *snr_sum_he = 0.0f;
    snr_sumt = 0;
    mssnr_hov = 0;
    snr_sum_HE_SAD = 0.0f;
    lp_snr = st->lp_speech - st->lp_noise;

    if ( lp_snr > 24.0f )
    {
        snr_idx = 0;
    }
    else if ( lp_snr > 18 )
    {
        snr_idx = 1;
    }
    else
    {
        snr_idx = 2;
    }

    if ( snr_idx == 0 )
    {
        stmp = 6;
        delta1 = 0.0f;
        delta2 = 0.0f;
        delta3 = 0.0f;
        delta4 = 0.0f;
        vad_thr = 2.4f*lp_snr - 42.2f;
        vad_thr = min(vad_thr, 80);
    }
    else if ( snr_idx == 1)
    {
        stmp = 6;
        delta1 = 0.1f;
        delta2 = 0.2f;
        delta3 = 0.2f;
        delta4 = 0.2f;
        vad_thr = 2.4f*lp_snr - 40.2f;
        vad_thr = min(vad_thr, 80);
    }
    else
    {
        stmp = 9;
        delta1 = 0.2f;
        delta2 = 0.4f;
        delta3 = 0.3f;
        delta4 = 0.4f;
        vad_thr = 2.5f*lp_snr - 10.0f;
        vad_thr = max(vad_thr, 1);
    }
    pt3 = st->bckr;
    nb_sig_snr = 20;

    for( i=st->min_band; i<=st->max_band; i++ )
    {
        ftmp  = *pt1++;
        ftmp1 = *pt2++;
        fr_enr = ( 0.2f * st->enrO[i] + 0.4f * ftmp + 0.4f * ftmp1 );

        if (ftmp > ftmp1)
        {
            snr[i] = ( 0.2f * st->enrO[i] + 0.4f * ftmp + 0.4f * ftmp1 ) / *pt3++;
        }
        else
        {
            snr[i] = ( 0.2f * st->enrO[i] + 0.3f * ftmp + 0.5f * ftmp1 ) / *pt3++;
        }

        if ( snr[i] < 2.0f )
        {
            nb_sig_snr--;
        }

        if ( snr[i] < 1 )
        {
            snr[i] = 1;
        }

        snr[i] = (float)log10(snr[i]);
        snr_sumt += snr[i];
        if (i < 2)
        {
            ftmp = snr[i] + delta1;
        }
        else if (i < 7)
        {
            ftmp = snr[i] + delta2;
        }
        else if (i <18)
        {
            ftmp = snr[i] + delta3;
        }
        else
        {
            ftmp = snr[i] + delta4;
        }
        ftmp1 = ftmp;
        if ( i < 7 )
        {
            ftmp1 = ftmp + 0.4f;
        }
        ftmp = min(ftmp, 2.0f);
        ftmp1 = min(ftmp1, 2.0f);
        msnr = 1;
        for (j=0; j<stmp; j++)
        {
            msnr *= ftmp;
        }
        mssnr += msnr;
        if ( i == 18 )
        {
            snr18 = msnr;
        }
        else if ( i == 19 )
        {
            snr19 = msnr;
        }
        msnr = 1;
        for (j=0; j<stmp; j++)
        {
            msnr *= ftmp1;
        }
        mssnr_hov += msnr;
        snr[i] = fr_enr / st->bckr[i];

        sign_thr_snr_acc( &snr_sum_HE_SAD, snr[i], sign_thr_HE_SAD, min_snr_HE_SAD );
        sign_thr_snr_acc( &snr_sum,snr[i], sign_thr, min_snr );

        /* To make snr[] compatible with older versions where snr[i] >= 1
           also this could be removed if this no longer is a requriement */
        if( snr[i] < 1.0f )
        {
            snr[i] = 1.0f;
        }
        /* accumulate background noise energy in bands [0-2]  and in bands [3-19]*/
        if( i < 3 )
        {
            accum_ener_L = accum_ener_L + st->bckr[i];
        }
        else
        {
            accum_ener_H = accum_ener_H + st->bckr[i];
        }

        /* identify the outlier band */
        if( snr[i] > snr_outlier )
        {
            snr_outlier = snr[i];
            snr_outlier_index = i;
        }
    }

    if ( (st->max_band == 19) && (snr[18] > 5.0f) && (snr[19] > 5.0f) )
    {
        ftmp = (mssnr + 3*(snr18 + snr19)) * 0.77f;
        if ( ftmp > mssnr )
        {
            mssnr = ftmp;
        }
    }
    else if ( snr_idx != 0 && nb_sig_snr > 13 )
    {
        if ( 2.5f*lp_snr - 15.5f > 0 )
        {
            mssnr += 2.5f*lp_snr - 15.5f;
        }
    }


    /* Separate SNR_SUM modification to */
    snr_sum_ol = snr_sum;
    if(st->max_band == 19 && snr_outlier < MAX_SNR_OUTLIER_3 && snr_outlier_index > 3 && snr_outlier_index < MAX_SNR_OUTLIER_IND)  /* Update the total SNR only for WB signals */
    {
        if( (accum_ener_L > OUTLIER_THR_1 * accum_ener_H ) || (snr_outlier < MAX_SNR_OUTLIER_1) )
        {
            snr_sum_ol = SNR_OUTLIER_WGHT_1 * (snr_sum_ol - snr_outlier);
        }
        else if( (accum_ener_L > OUTLIER_THR_2 * accum_ener_H ) || (snr_outlier < MAX_SNR_OUTLIER_2) )
        {
            snr_sum_ol = SNR_OUTLIER_WGHT_2 * (snr_sum_ol - snr_outlier);
        }
        else
        {
            snr_sum_ol = SNR_OUTLIER_WGHT_3 * (snr_sum_ol - snr_outlier);
        }
    }

    st->snr_sum_vad = 0.5f * st->snr_sum_vad + 0.5f * snr_sum_ol;

    snr_sum_ol = 10.0f * (float)log10( snr_sum_ol );
    snr_sum = snr_sum_ol; /* for NB no outlier modification */

    snr_sum_HE_SAD = 10.0f * (float)log10( snr_sum_HE_SAD );
    *snr_sum_he=snr_sum_HE_SAD;

    /*---------------------------------------------------------------------*
     * compute thr1 for SAD decision
     *---------------------------------------------------------------------*/

    lp_snr = st->lp_speech - st->lp_noise;            /* long-term SNR */

    if (lp_snr < st->sign_dyn_lp)
    {
        lp_snr +=1;

        if (lp_snr > st->sign_dyn_lp)
        {
            lp_snr = st->sign_dyn_lp;
        }
    }

    thr1 = nk * lp_snr + nc + nv * ( st->Etot_v_h2 - nv_ofs);              /* Linear function for noisy speech */

    if (lp_snr > 20.0f )
    {
        thr1 = thr1 + 0.3f * (lp_snr - 20.0f);
        if ( st->max_band==16 && lp_snr > 40 && thr1 > 24.1f && st->lp_speech < 45.0f )
        {
            thr1 = 24.1f;
        }
    }

    /*---------------------------------------------------------------------*
     * WB input
     * SNR threshold computing
     * Hangover control & final VAD decision
     *---------------------------------------------------------------------*/

    if( st->input_bwidth != NB )
    {
        /* Outlier Detection first calculates thr1_ol and snr_sum_ol instead of
           thr1 and snr_sum */

        thr1_ol = thr1;
        if( lp_snr < th_clean )
        {
            hangover_short = 4;
            if( ( snr_outlier_index <= 4 && (st->last_coder_type > UNVOICED) && !st->Opt_SC_VBR ) ||
                    ( snr_outlier_index <= 4 && (st->last_7k2_coder_type > UNVOICED) && st->Opt_SC_VBR ) )
            {
                thr1_ol = thr1 - 1.0f ;
                snr_sum_ol = 10.0f * (float)log10( st->snr_sum_vad );
            }
            else if ( ( (st->last_coder_type <= UNVOICED) && (snr_outlier < MAX_SNR_OUTLIER_2)  && !st->Opt_SC_VBR ) ||
                      ( (st->last_7k2_coder_type <= UNVOICED) && (snr_outlier < MAX_SNR_OUTLIER_2)  && st->Opt_SC_VBR ) )

            {
                thr1_ol = thr1 + (float)(1.0f - 0.04f * snr_outlier);
            }
            else
            {
                thr1_ol = thr1 + max(0, (float)(0.6f - 0.01f * snr_outlier));
            }
        }
        else
        {
            if( st->Opt_SC_VBR )
            {
                hangover_short = 3;
            }
            else
            {
                hangover_short = 3;
            }
        }

        /* The use of outlier detection had been removed by accident at some point */
        snr_sum = snr_sum_ol;
        thr1 = thr1_ol;

        /* DTX HANGOVER ADDITION MOVED TO pre_proc() */

        flag_he1 = 0;
        *localVAD = 0;
        if ( mssnr > vad_thr )
        {
            *localVAD = 1;                /* he1 primary decision */
            flag_he1 = 1;
            st->nb_active_frames_he1++;   /* Counter of consecutive active speech frames */
            if ( st->nb_active_frames_he1 >= ACTIVE_FRAMES )
            {
                st->nb_active_frames_he1 = ACTIVE_FRAMES;
                st->hangover_cnt_he1 = 0;   /* Reset the counter of hangover frames after at least "active_frames" speech frames */
            }
            /* inside HO period */
            if ( st->hangover_cnt_he1 < HANGOVER_LONG_HE && st->hangover_cnt_he1 != 0 )
            {
                st->hangover_cnt_he1++;
            }

            if ( st->soft_hangover > 0 )
            {
                st->soft_hangover--;
            }
        }
        else
        {
            /* Reset the counter of speech frames necessary to start hangover algorithm */
            st->nb_active_frames_he1 = 0;
        }


        if ( st->voiced_burst > 3 )
        {
            if ( st->bcg_flux < 40 )
            {
                st->soft_hangover = hangover_sf_tbl[snr_idx+3];
            }
            else
            {
                st->soft_hangover = hangover_sf_tbl[snr_idx];
            }
        }


        hangover_hd = hangover_hd_tbl[snr_idx];

        if ( st->bcg_flux < 40 )
        {
            hangover_hd = (hangover_hd>>1) + 1;
        }


        if ( flag_he1 == 0 && st->soft_hangover > 0 )
        {
            if ( mssnr_hov > vad_thr )
            {
                flag_he1 = 1;
                st->soft_hangover--;
            }
            else
            {
                st->soft_hangover = 0;
            }

            if ( st->soft_hangover < 0 )
            {
                st->soft_hangover = 0;
            }
        }

        if ( flag_he1 == 0 && st->hangover_cnt_he1 < hangover_hd && st->soft_hangover == 0 )
        {
            flag_he1 = 1;
            st->hangover_cnt_he1++;
        }

        /* Calculate background stationarity */
        if ( flag_he1 == 0 && st->first_noise_updt > 0 )
        {
            if ( snr_sumt > st->bcg_flux )
            {
                if ( st->bcg_flux_init-- > 0 )
                {
                    if ( snr_sumt > st->bcg_flux+50 )
                    {
                        st->bcg_flux = 0.9f * st->bcg_flux + (1-0.9f)*(st->bcg_flux+50);
                    }
                    else
                    {
                        st->bcg_flux = 0.9f * st->bcg_flux + (1-0.9f)*snr_sumt;
                    }
                }
                else
                {
                    if ( snr_sumt > st->bcg_flux+10 )
                    {
                        st->bcg_flux = 0.99f * st->bcg_flux + (1-0.99f)*(st->bcg_flux+10);
                    }
                    else
                    {
                        st->bcg_flux = 0.99f * st->bcg_flux + (1-0.99f)*snr_sumt;
                    }
                }
            }
            else
            {
                if ( st->bcg_flux_init-- > 0 )
                {
                    if ( snr_sumt < st->bcg_flux-30 )
                    {
                        st->bcg_flux = 0.95f * st->bcg_flux + (1-0.95f)*(st->bcg_flux-30);
                    }
                    else
                    {
                        st->bcg_flux = 0.95f * st->bcg_flux + (1-0.95f)*snr_sumt;
                    }
                }
                else
                {
                    if ( snr_sumt < st->bcg_flux-10 )
                    {
                        st->bcg_flux = 0.9992f * st->bcg_flux + (1-0.9992f)*(st->bcg_flux-10);
                    }
                    else
                    {
                        st->bcg_flux = 0.9992f * st->bcg_flux + (1-0.9992f)*snr_sumt;
                    }
                }
            }

            if ( st->bcg_flux_init < 0 )
            {
                st->bcg_flux_init = 0;
            }
        }

        flag = 0;
        *localVAD = 0;

        if ( snr_sum > thr1 && flag_he1 == 1 ) /* Speech present */
        {
            flag = 1;
            *localVAD = 1;
            st->nb_active_frames++;   /* Counter of consecutive active speech frames */
            if (st->nb_active_frames >= ACTIVE_FRAMES )
            {
                st->nb_active_frames = ACTIVE_FRAMES;
                st->hangover_cnt = 0;   /* Reset the counter of hangover frames after at least "active_frames" speech frames */
            }

            /* inside HO period */
            if(st->hangover_cnt < HANGOVER_LONG && st->hangover_cnt != 0 )
            {
                st->hangover_cnt++;
            }
        }
        else
        {
            /* Reset the counter of speech frames necessary to start hangover algorithm */
            st->nb_active_frames = 0;
            if(st->hangover_cnt < HANGOVER_LONG )    /* inside HO period */
            {
                st->hangover_cnt++;
            }

            if(st->hangover_cnt <= hangover_short )  /* "hard" hangover  */
            {
                /* send the extra 3 HO frames to NELP */
                if ( (lp_snr < th_clean) && (st->Opt_SC_VBR) && (st->hangover_cnt >= 2) )
                {
                    *noisy_speech_HO = 1;
                }

                if ( (lp_snr >= th_clean) && (st->Opt_SC_VBR) && (st->hangover_cnt >= 2) )
                {
                    *clean_speech_HO = 1;
                }

                flag = 1;
            }

        }

        /* localVAD and vad_flag for HE-SAD - in parallel with normal localVAD and vad_flag */
        *localVAD_HE_SAD = 0;
        if( snr_sum_HE_SAD > thr1 && (flag_he1 == 1) ) /* Speech present */
        {
            *localVAD_HE_SAD = 1;
        }
    }

    /*---------------------------------------------------------------------*
     * NB input
     * SNR threshold computing
     * Hangover control & final VAD decision
     *---------------------------------------------------------------------*/

    else                                          /* NB input */
    {
        /* Add localVAD_HE_SAD also for NB operation for use with speech music classifier */
        *localVAD_HE_SAD = 0;
        if (snr_sum_HE_SAD > thr1 )
        {
            *localVAD_HE_SAD = 1;
        }

        *localVAD        = 0;       /* init needed in NB,           otherwise it can be undefined          */
        if ( snr_sum > thr1 )                     /* Speech present */
        {
            st->nb_active_frames++;                /* Counter of consecutive active speech frames */
            if (st->nb_active_frames >= ACTIVE_FRAMES )
            {
                st->nb_active_frames = ACTIVE_FRAMES;
                st->hangover_cnt = 0;                /* Reset the counter of hangover frames after at least "active_frames" speech frames */
            }

            *localVAD = 1;
        }
        else
        {
            st->nb_active_frames = 0;    /* Reset the counter of speech frames necessary to start hangover algorithm */
        }

        if(st->hangover_cnt < HANGOVER_LONG_NB )
        {
            st->hangover_cnt++;
            if( lp_snr < 19.0f )                  /* very low SNR */
            {
                thr1 -= 5.2f;
            }
            else if( lp_snr < 35.0f )             /* low SNR */
            {
                thr1 -= 2.0f;
            }
        }

        if (st->Opt_DTX_ON)
        {
            if (lp_snr < th_clean)
            {
                thr2 = thr1 - 1.10f;

            }
            else
            {
                thr2 = thr1 - 1.5f;
            }
        }
        else
        {
            if (lp_snr < th_clean)
            {
                thr2 = thr1 - 1.3f;
            }
            else
            {
                thr2 = thr1 - 1.5f;
            }
        }

        flag = 0;
        if ( snr_sum > thr1 )                      /* Speech present */
        {
            flag = 1;
        }

        if (  (snr_sum < thr1) && (snr_sum > thr2) )   /* Speech present */
        {
            flag = 1;
            *localVAD = 0;
            *NB_speech_HO = 1;
        }

        /* Need to handle the case when switching from WB -> NB */
    }



    if( st->input_bwidth != NB )
    {
        *flag_noisy_speech_snr = (lp_snr < TH16_2_NFLAG );          /*original threshold: 35dB*/
    }
    else
    {
        *flag_noisy_speech_snr = (lp_snr < TH8_1_NFLAG );    /*original threshold: 20dB, not yet tested!*/
    }

    /* SC-VBR */
    st->vadsnr = snr_sum;
    st->vadnoise = thr1;

    /* Updates */
    st->prim_act_quick = 0.2f * (*localVAD)  + (1.0f -0.2f)* st->prim_act_quick;
    st->prim_act_slow = 0.01f * (*localVAD)  + (1.0f-0.01f)* st->prim_act_slow;
    if (st->prim_act_quick <= st->prim_act_slow)
    {
        st->prim_act = 0.1f * st->prim_act_quick + (1.0f-0.1f)* st->prim_act;
    }
    else
    {
        st->prim_act = 0.1f * st->prim_act_slow + (1.0f-0.1f)* st->prim_act;
    }

    st->prim_act_quick_he = 0.2f * *localVAD_HE_SAD + (1.0f -0.2f)* st->prim_act_quick_he;
    st->prim_act_slow_he = 0.01f * *localVAD_HE_SAD + (1.0f-0.01f)* st->prim_act_slow_he;

    if (st->prim_act_quick_he <= st->prim_act_slow_he)
    {
        st->prim_act_he = 0.1f * st->prim_act_quick_he + (1.0f-0.1f)* st->prim_act_he;
    }
    else
    {
        st->prim_act_he = 0.1f * st->prim_act_slow_he + (1.0f-0.1f)* st->prim_act_he;
    }


    if ((st->vad_flag_reg_H & (long) 0x40000L) != 0) /* 0x4000L = 0x01L << 18 */
    {
        st->vad_flag_cnt_50 = st->vad_flag_cnt_50-1;
    }

    st->vad_flag_reg_H = (st->vad_flag_reg_H & (long) 0x3fffffffL ) << 1;

    if ( ( st->vad_flag_reg_L & (long) 0x40000000L) != 0)
    {
        st->vad_flag_reg_H = st->vad_flag_reg_H |  0x01L;
    }

    st->vad_flag_reg_L = (st->vad_flag_reg_L & (long) 0x3fffffffL ) << 1;

    if ( flag ) /* should not include the extra DTX hangover */
    {
        st->vad_flag_reg_L = st->vad_flag_reg_L | 0x01L;
        st->vad_flag_cnt_50 = st->vad_flag_cnt_50+1;
    }

    if ((st->vad_prim_reg & (long) 0x8000L ) != 0) /* 0x8000L = 1L << 15 */
    {
        st->vad_prim_cnt_16 = st->vad_prim_cnt_16-1;
    }

    st->vad_prim_reg = (st->vad_prim_reg & (long) 0x3fffffffL ) << 1;

    if(*localVAD)
    {
        st->vad_prim_reg = st->vad_prim_reg | 0x01L;
        st->vad_prim_cnt_16 = st->vad_prim_cnt_16+1;
    }



    return flag;
}
