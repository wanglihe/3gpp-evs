/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_enc.h"
#include "rom_com.h"    /* Common static table prototypes         */


/*---------------------------------------------------------------------*
 * Local constants
 *---------------------------------------------------------------------*/

#define ATT_NSEG              32
#define ATT_SEG_LEN           (L_FRAME/ATT_NSEG)
#define ATT_3LSUB_POS         (3 * ATT_NSEG / NB_SUBFR)
#define ATT_3LSUB_POS_16k     (short)((4.0f * ATT_NSEG / (float)NB_SUBFR16k) + 0.5f)

/*---------------------------------------------------------------------*
 * Local functions
 *---------------------------------------------------------------------*/

static short sp_mus_classif_1st( Encoder_State *st, const short localVAD, const short pitch[3], const float voicing[3], const float lsp_new[M],
                                 const float cor_map_sum, const float epsP[M+1], const float PS[], float non_sta, float relE,
                                 float *voi_fv, float *cor_map_sum_fv, float *LPCErr, short *high_lpn_flag );

static void music_mixed_classif_improv( Encoder_State *st, const float *new_inp, short *sp_aud_decision1,
                                        const short vad_flag, const float *voicing, const float *epsP, const float etot,
                                        const float old_cor, const float cor_map_sum );

static void spec_analysis( float *Bin_E, float *p2v_map );

static void flux( float *Bin_E, float *p2v_map, float *old_Bin_E, float *buf_flux, short attack_hangover, float dec_mov );

static void tonal_dist( float *p2v_map, float *buf_pkh, float *buf_Ntonal, float *buf_Ntonal2, float *buf_Ntonal_lf );

static short mode_decision( Encoder_State *st, short len, float *dec_mov, float *buf_flux, float *buf_epsP_tilt, float *buf_pkh, float *buf_cor_map_sum,
                            float *buf_Ntonal, float *buf_Ntonal2, float *buf_Ntonal_lf, float *buf_dlp );


static void sp_mus_classif_2nd( Encoder_State *st, const short sp_aud_decision1, short *sp_aud_decision2, const short pitch[],
                                const float Etot, short *coder_type, short *attack_flag,
                                const float *inp, const short localVAD, const short vad_flag );

static void var_cor_calc( const float old_corr, float *mold_corr, float var_cor_t[], short *high_stable_cor );

static short attack_det( const float *inp, const short last_clas, const short localVAD, const short coder_type, const int total_brate );

static void tonal_context_improv( Encoder_State *st, const float PS[], short *sp_aud_decision1, short *sp_aud_decision2,
                                  const short vad_flag, const short pitch[3], const float voicing[3], const float voi_fv,
                                  const float cor_map_sum_fv, const float LPCErr );
static void order_spectrum( float *vec, short len );

static void detect_sparseness( Encoder_State *st, const short localVAD_HE_SAD, short *sp_aud_decision1,
                               short *sp_aud_decision2, const float voi_fv );


/*---------------------------------------------------------------------*
 * speech_music_classif()
 *
 * Speech/music classification
 *
 * The following technologies are used based on the outcome of the sp/mus classifier
 * sp_aud_decision1  sp_aud_decision2
 *       0                 0             use ACELP (+TD BWE)
 *       1                 0             use ACELP (+FD BWE) or HQ/LR-MDCT depending on bitrate
 *       1                 1             use GSC (+FD BWE) or HQ/LR-MDCT depending on bitrate
 *
 *       0                 1             exceptionally use GSC (+FD BWE) instead of LR-MDCT at 13.2 kbps (WB/SWB) for sparse spectra
 *---------------------------------------------------------------------*/

void speech_music_classif(         /* o  : 1st stage decision (1-music, 0-speech or noise) */
    Encoder_State *st,           /* i/o: state structure                                 */
    short *sp_aud_decision0, /* o  : 1st stage speech/music decision                 */
    short *sp_aud_decision1, /* o  : 1st stage speech/music decision                 */
    short *sp_aud_decision2, /* o  : 2nd stage speech/music decision                 */
    const float *new_inp,          /* i  : new input signal                                */
    const float *inp,              /* i  : input signal to locate attach position          */
    const short vad_flag,
    const short localVAD,
    const short localVAD_HE_SAD,   /* i  : HE-SAD flag without hangover                    */
    const short pitch[3],          /* i  : open-loop pitch estimate in three subframes     */
    const float voicing[3],        /* i  : voicing estimate in three subframes             */
    const float lsp_new[M],        /* i  : LSPs in current frame                           */
    const float cor_map_sum,       /* i  : correlation map sum (from multi-harmonic anal.) */
    const float epsP[M+1],         /* i  : LP prediciton error                             */
    const float PS[],              /* i  : energy spectrum                                 */
    const float Etot,              /* i  : total frame energy                              */
    const float old_cor,           /* i  : max correlation from previous frame             */
    short *coder_type,       /* i/o: coding type                                     */
    short *attack_flag,      /* o  : flag to indicate if attack is to be treated by TC or GSC */
    const float non_sta,           /* i  : unbound non-stationarity for sp/mus classifier */
    const float relE,              /* i  : relative frame energy */
    short *high_lpn_flag,
    const short flag_spitch        /* i  : flag to indicate very short stable pitch */
)
{
    float voi_fv, cor_map_sum_fv, LPCErr;

    /* 1st stage speech/music classification based on the GMM model */
    *sp_aud_decision1 = sp_mus_classif_1st( st, localVAD_HE_SAD, pitch, voicing, lsp_new, cor_map_sum,
                                            epsP, PS, non_sta, relE, &voi_fv, &cor_map_sum_fv, &LPCErr, high_lpn_flag );

    if( st->codec_mode == MODE1 || st->sr_core == 12800)
    {

        /* Improvement of the 1st stage decision for mixed/music content */
        if ( !st->Opt_SC_VBR && (st->total_brate != ACELP_24k40) )
        {
            music_mixed_classif_improv( st, new_inp, sp_aud_decision1, vad_flag, voicing, epsP, Etot, old_cor, cor_map_sum );
        }

        *sp_aud_decision0 = *sp_aud_decision1;

        /* 2nd stage speech/music classification (rewrite music to speech in onsets) */
        *sp_aud_decision2 = *sp_aud_decision1;

        if( st->bwidth > NB )
        {
            sp_mus_classif_2nd( st, *sp_aud_decision1, sp_aud_decision2, pitch, Etot, coder_type, attack_flag,
                                inp, localVAD, vad_flag );

            if ( flag_spitch && st->bwidth == WB && st->total_brate < ACELP_13k20 )
            {
                /* avoid switch to AUDIO/MUSIC class for very short stable high pitch
                   and/or stable pitch with high correlation at low bitrates*/
                *sp_aud_decision2 = 0;
            }
        }

        /* Context-based improvement of 1st and 2nd stage decision on stable tonal signals */
        if( !st->Opt_SC_VBR && (st->total_brate != ACELP_24k40) )
        {
            tonal_context_improv( st, PS, sp_aud_decision1, sp_aud_decision2, vad_flag, pitch, voicing,
                                  voi_fv, cor_map_sum_fv, LPCErr );
        }

        /* Avoid using LR-MDCT on sparse spectra, use GSC instead at 13.2 kbps (WB/SWB) */
        if ( !st->Opt_SC_VBR && st->total_brate == 13200 && vad_flag == 1 &&
                ( st->bwidth == WB || st->bwidth == SWB ) )
        {
            detect_sparseness( st, localVAD_HE_SAD, sp_aud_decision1, sp_aud_decision2, voi_fv );
        }

        /* override speech/music classification to ACELP when background noise level reaches certain level */
        /* this is a patch against mis-classifications during active noisy speech segments */
        if ( st->lp_noise > 12.0f )
        {
            *sp_aud_decision1 = 0;
            *sp_aud_decision2 = 0;
        }

        /* set GSC noisy speech flag on unvoiced SWB segments */
        st->GSC_noisy_speech = 0;
        if ( vad_flag == 1 && st->total_brate >= ACELP_13k20 && st->total_brate < ACELP_24k40 &&
                st->lp_noise > 12.0f  && *sp_aud_decision1 == 0 && st->bwidth >= SWB &&
                st->coder_type_raw == UNVOICED )
        {
            st->GSC_noisy_speech = 1;
        }

        /* Select AUDIO frames */
        if ( st->codec_mode == MODE1 && (*sp_aud_decision2 || st->GSC_noisy_speech) )
        {
            *coder_type = AUDIO;
            st->noise_lev = NOISE_LEVEL_SP0;
        }

    }

    return;
}


/*---------------------------------------------------------------------*
 * sp_mus_classif_1st()
 *
 * 1st stage speech/music classification (based on the GMM model)
 *---------------------------------------------------------------------*/

static short sp_mus_classif_1st(   /* o  : decision flag (1-music, 0-speech or noise)      */
    Encoder_State *st,           /* i/o: state structure                                 */
    const short localVAD,
    const short pitch[3],          /* i  : open-loop pitch estimate in three subframes     */
    const float voicing[3],        /* i  : voicing estimate in three subframes             */
    const float lsp_new[M],        /* i  : LSPs in current frame                           */
    const float cor_map_sum,       /* i  : correlation map sum (from multi-harmonic anal.) */
    const float epsP[M+1],         /* i  : LP prediciton error                             */
    const float PS[],              /* i  : energy spectrum                                 */
    float non_sta,           /* i  : unbound non-stationarity                        */
    float relE,              /* i  : relative frame energy                           */
    float *voi_fv,           /* o  : scaled voicing feature                          */
    float *cor_map_sum_fv,   /* o  : scaled correlation map feature                  */
    float *LPCErr,           /* o  : scaled LP prediction error feature              */
    short *high_lpn_flag

)
{
    short i, k, p, dec, vad;
    float dlp, ftmp, lepsP1, sum_PS, ps_diff, ps_sta, wrelE, wdrop, wght, mx;
    float FV[N_FEATURES], *pFV = FV, PS_norm[128], dPS[128], lsp[M];
    float pys, pym, xm[N_FEATURES], py, lps = 0, lpm = 0;
    const float *pSF;
    float pyn, lpn = 0;

    /*------------------------------------------------------------------*
     * Initialization
     *------------------------------------------------------------------*/

    vad = localVAD;

    /*------------------------------------------------------------------*
     * Preparation of the feature vector
     *------------------------------------------------------------------*/

    /* [0] OL pitch */
    if ( st->tc_cnt == 1 || st->tc_cnt == 2 )
    {
        *pFV++ = (float)pitch[2];
    }
    else
    {
        *pFV++ = (float)(pitch[0] + pitch[1] + pitch[2]) / 3.0f;
    }

    /* [1] voicing */
    if ( st->tc_cnt == 1 || st->tc_cnt == 2 )
    {
        *pFV++ = voicing[2];
    }
    else
    {
        *pFV++ = (float)(voicing[0] + voicing[1] + voicing[2]) / 3.0f;
    }

    /* [2,3,4,5,6] LSFs */
    mvr2r( lsp_new, lsp, M );

    ftmp = (float)acos(lsp[1]);
    *pFV++ = ftmp + st->last_lsp[1];
    st->last_lsp[1] = ftmp;

    ftmp = (float)acos(lsp[2]);
    *pFV++ = ftmp + st->last_lsp[2];
    st->last_lsp[2] = ftmp;

    ftmp = (float)acos(lsp[3]);
    *pFV++ = ftmp + st->last_lsp[3];
    st->last_lsp[3] = ftmp;

    ftmp = (float)acos(lsp[4]);
    *pFV++ = ftmp + st->last_lsp[4];
    st->last_lsp[4] = ftmp;

    ftmp = (float)acos(lsp[5]);
    *pFV++ = ftmp + st->last_lsp[5];
    st->last_lsp[5] = ftmp;

    /* [7] cor_map_sum */
    *pFV++ = cor_map_sum + st->last_cor_map_sum;
    st->last_cor_map_sum = cor_map_sum;

    /* [8] non_sta */
    *pFV++ = non_sta + st->last_non_sta;
    st->last_non_sta = non_sta;

    /* [9] epsP */
    if ( st->bwidth == NB )
    {
        /* do not take into account (statistics are too different) */
        *pFV++ = -1.647f;
    }
    else
    {
        lepsP1 = (float)log(epsP[1] + 1e-5f);
        ftmp = (float)log(epsP[13]) - lepsP1;
        *pFV++ = ftmp + st->past_epsP2;
        st->past_epsP2 = ftmp;
    }

    /* calculation of differential normalized power spectrum */
    sum_PS = 1e-5f;
    for ( i = LOWEST_FBIN; i < HIGHEST_FBIN; i++ )
    {
        sum_PS += PS[i];
    }

    for ( i = LOWEST_FBIN; i < HIGHEST_FBIN; i++ )
    {
        PS_norm[i] = PS[i] / sum_PS;
        dPS[i] = (float)fabs(PS_norm[i] - st->past_PS[i-LOWEST_FBIN]);
    }

    /* [10] ps_diff (spectral difference) */
    ps_diff = 0;
    for ( i = LOWEST_FBIN; i < HIGHEST_FBIN; i++ )
    {
        ps_diff += dPS[i];
    }

    ps_diff = (float)log(ps_diff + 1e-5f);
    *pFV++ = ps_diff + st->past_ps_diff;
    st->past_ps_diff = ps_diff;

    /* [11] ps_sta (spectral stationarity) */
    ps_sta = 0;
    for ( i = LOWEST_FBIN; i < HIGHEST_FBIN; i++ )
    {
        mx = PS_norm[i] > st->past_PS[i-LOWEST_FBIN] ? PS_norm[i] : st->past_PS[i-LOWEST_FBIN];
        ps_sta += mx / (dPS[i] + 1e-5f);
    }

    *pFV++ = (float)log(ps_sta + 1e-5f);
    mvr2r( &PS_norm[LOWEST_FBIN], st->past_PS, HIGHEST_FBIN-LOWEST_FBIN );

    /*------------------------------------------------------------------*
     * Scaling of the feature vector
     *------------------------------------------------------------------*/

    pFV = FV;
    if ( st->bwidth == NB )
    {
        pSF = SF_8k;
    }
    else
    {
        pSF = SF;
    }

    for ( i=0; i<N_FEATURES; i++, pFV++, pSF += 2 )
    {
        *pFV = pSF[0] **pFV + pSF[1];
    }

    /* store some scaled parameters for later correction of the 1st stage speech/music classification */
    *voi_fv = FV[1];
    *cor_map_sum_fv = FV[7];
    *LPCErr = FV[9];

    /*------------------------------------------------------------------*
     * Calculation of posterior probability
     * Log-probability
     *------------------------------------------------------------------*/

    pys = pym = pyn = 1e-5f;

    /* run loop for all mixtures (for each mixture, calculate the probability of speech or noise and the probability of music) */
    for ( k = 0; k < N_MIXTURES; k++ )
    {
        /* active frames - calculate the probability of speech */
        for ( p = 0; p < N_FEATURES; p++ )
        {
            xm[p] = FV[p] - m_speech[k*N_FEATURES+p];
        }

        py = lvm_speech[k] + dot_product_mat(xm, &invV_speech[k*N_FEATURES*N_FEATURES], N_FEATURES );
        pys += (float)exp(py);
        /* inactive frames - calculate the probability of noise */
        for ( p = 0; p < N_FEATURES; p++ )
        {
            xm[p] = FV[p] - m_noise[k*N_FEATURES+p];
        }

        py = lvm_noise[k] + dot_product_mat(xm, &invV_noise[k*N_FEATURES*N_FEATURES], N_FEATURES );
        pyn += (float)exp(py);

        /* either active or inactive frames - calculate the probability of music */
        for ( p = 0; p < N_FEATURES; p++ )
        {
            xm[p] = FV[p] - m_music[k*N_FEATURES+p];
        }

        py = lvm_music[k] + dot_product_mat(xm, &invV_music[k*N_FEATURES*N_FEATURES], N_FEATURES );
        pym += (float)exp(py);
    }

    /* calculate log-probability */
    lps = (float)log(pys) - 0.5f * N_FEATURES * (float)log(2*EVS_PI);
    lpm = (float)log(pym) - 0.5f * N_FEATURES * (float)log(2*EVS_PI);
    lpn = (float)log(pyn) - 0.5f * N_FEATURES * (float)log(2*EVS_PI);

    *high_lpn_flag = 0;
    if ( lpn > lps && lpn > lpm )
    {
        *high_lpn_flag = 1;
    }
    if ( !vad )
    {
        /* artificially increase log-probability of noise */
        lps = lpn * 1.2f;
    }

    st->lpm = lpm;
    st->lps = lps;

    /* determine HQ Generic speech class */
    if ( lps > lpm + 0.5f )
    {
        st->hq_generic_speech_class = 1;
    }
    else
    {
        st->hq_generic_speech_class = 0;
    }

    /*------------------------------------------------------------------*
     * State machine (sp_mus_state < 0 .. inactive, > 0 .. entry, = 0 .. active )
     *------------------------------------------------------------------*/

    if ( vad )
    {
        if ( relE < -20 || (lps <= -5 && lpm <= -5) )
        {
            if ( st->sp_mus_state > 0 )
            {
                if ( st->sp_mus_state < HANG_LEN )
                {
                    /* energy is too low but we are in entry period -> reset the inactive counter to allow new entry later */
                    st->inact_cnt = 0;
                }

                /* energy is too low -> we are going to instable state */
                st->sp_mus_state = 0;
            }
            else if ( st->sp_mus_state > -HANG_LEN )
            {
                /* energy is still too low -> we are still in instable state */
                st->sp_mus_state--;
            }
        }
        else if ( st->sp_mus_state <= 0 )
        {
            if ( st->inact_cnt == 0 )
            {

                st->sp_mus_state = 1;
            }
            else
            {

                st->sp_mus_state = HANG_LEN;
            }

            st->inact_cnt = 12;
        }
        else if ( st->sp_mus_state > 0 && st->sp_mus_state < HANG_LEN )
        {
            /* we are inside an entry period -> increment the counter of entry frames */
            st->sp_mus_state++;
        }

        if ( st->sp_mus_state < 0 && st->inact_cnt > 0 )
        {
            st->inact_cnt--;
        }
    }
    else
    {
        if ( st->sp_mus_state > 0 && st->sp_mus_state < HANG_LEN )
        {
            st->inact_cnt = 0;
        }
        else if ( st->inact_cnt > 0 )
        {
            st->inact_cnt--;
        }

        if ( st->sp_mus_state > 0 && st->sp_mus_state < HANG_LEN )
        {

            st->sp_mus_state = -HANG_LEN;
        }
        else if ( st->sp_mus_state > 0 )
        {

            st->sp_mus_state = -1;
        }
        else if ( st->sp_mus_state > -HANG_LEN )
        {
            /* we are in inactive state */
            st->sp_mus_state--;
        }

    }

    /*------------------------------------------------------------------*
     * Decision without hangover
     * Weighted decision
     *------------------------------------------------------------------*/

    /* decision without hangover (0 - speech/noise, 1 - music) */
    dec = lpm > lps;
    dlp = lpm - lps;

    if ( !vad )
    {
        dec = 0;
        dlp = 0;
    }

    /* calculate weight based on relE (close to 0.01 in low-E regions, close to 1 in high-E regions) */
    wrelE = 1.0f + relE/15;

    if ( wrelE > 1.0f )
    {
        wrelE = 1.0f;
    }
    else if ( wrelE < 0.01f )
    {
        wrelE = 0.01f;
    }

    /* calculate weight based on drops of dlp (close to 1 during sudden drops of dlp, close to 0 otherwise) */
    if ( dlp < 0 && dlp < st->past_dlp[0] )
    {
        if ( st->past_dlp[0] > 0 )
        {
            st->wdrop = -dlp;
        }
        else
        {
            st->wdrop += st->past_dlp[0] - dlp;
        }
    }
    else
    {
        st->wdrop = 0;
    }

    wdrop = st->wdrop/20;

    if ( wdrop > 1.0f )
    {
        wdrop = 1.0f;
    }
    else if ( wdrop < 0.1f )
    {
        wdrop = 0.1f;
    }

    /* combine weights into one */
    wght = wrelE * wdrop;
    if ( wght < 0.01f )
    {
        wght = 0.01f;
    }

    /* calculate weighted decision */
    st->wdlp_0_95_sp = wght * dlp + (1 - wght) * st->wdlp_0_95_sp;

    if ( st->sp_mus_state == -HANG_LEN )
    {
        st->wdlp_0_95_sp = 0;
    }

    /*------------------------------------------------------------------*
     * Final speech/music decision
     *------------------------------------------------------------------*/

    if ( !vad && st->sp_mus_state == -HANG_LEN )
    {
        /* inactive state */
        dec = 0;
    }
    else if ( st->sp_mus_state <= 0 )
    {
        /* transition from active to inactive state or instable state */
        dec = st->past_dec[0];
    }
    else if ( st->sp_mus_state > 0 && st->sp_mus_state < HANG_LEN )
    {
        /* entry state -> final decision is calculated based on weighted average of past non-binary decisions */
        ftmp = w[st->sp_mus_state-1][0] * dlp;
        ftmp += dotp( &w[st->sp_mus_state-1][1], st->past_dlp, HANG_LEN-1 );
        dec = ftmp > 2.0f;
    }
    else
    {
        /* stable active state */
        if ( st->wdlp_0_95_sp > 0 && st->past_dec[0] == 0 && st->past_dec[1] == 0 && st->past_dec[2] == 0 )
        {
            /* switching from speech to music */
            dec = 1;
        }
        else if ( st->past_dec[0] == 1 && st->wdlp_0_95_sp < 0 )
        {
            /* switching from music to speech */
            dec = 0;
        }
        else
        {
            dec = st->past_dec[0];
        }
    }

    /*------------------------------------------------------------------*
     * Updates
     *------------------------------------------------------------------*/

    /* update buffer of past non-binary decisions */
    mvr2r( &st->past_dlp[0], &st->past_dlp[1], HANG_LEN-2 );
    st->past_dlp[0] = dlp;

    /* update buffer of past binary decisions */
    mvs2s( &st->past_dec[0], &st->past_dec[1], HANG_LEN-2 );
    st->past_dec[0] = dec;

    return dec;
}


/*---------------------------------------------------------------------*
 * sp_mus_classif_2nd()
 *
 * 2nd stage speech/music classifier (convert music to speech for onsets)
 *---------------------------------------------------------------------*/

static void sp_mus_classif_2nd(
    Encoder_State *st,                  /* i/o: encoder state structure                */
    const short sp_aud_decision1,     /* i  : 1st stage decision flag                */
    short *sp_aud_decision2,    /* o  : 2nd stage decision flag                */
    const short pitch[3],             /* i  : open-loop pitch estimate in three subframes     */
    const float Etot,                 /* i  : total frame energy                     */
    short *coder_type,          /* i/o: coder type                             */
    short *attack_flag,         /* i/o: attack flag (GSC or TC)                */
    const float *inp,                 /* i  : input signal                           */
    const short localVAD,
    const short vad_flag
)
{
    short attack;

    /* initialization */
    *attack_flag = 0;

    /* signal stability estimation */
    stab_est( Etot, st->gsc_lt_diff_etot, &st->gsc_mem_etot,
              &st->gsc_nb_thr_3, &st->gsc_nb_thr_1, st->gsc_thres, &st->gsc_last_music_flag, vad_flag );

    /* calculate variance of correlation */
    var_cor_calc( st->old_corr, &st->mold_corr, st->var_cor_t, &st->high_stable_cor );

    /* attack detection */
    attack = attack_det( inp, st->clas, localVAD, *coder_type, st->total_brate );

    /* change decision from music hto speech in certain special cases */
    if( sp_aud_decision1 == 1 )
    {
        if( st->ener_RAT < 0.18f && st->lt_dec_thres > 15.0f )
        {
            /* strong music decision but almost no content below 1kHz */
            *sp_aud_decision2 = 0;
        }
        else if( st->high_stable_cor && pitch[0] >= 130 )
        {
            /* prevent GSC in highly correlated signal with low energy variation */
            /* this is basically a patch against bassoon-type of music */
            *sp_aud_decision2 = 0;

            if ( st->codec_mode == MODE1 && *coder_type == TRANSITION )
            {
                *coder_type = GENERIC;
            }
        }
        else if( st->gsc_lt_diff_etot[MAX_LT-1] > 4.5f &&
                 (st->gsc_lt_diff_etot[MAX_LT-1] - st->gsc_lt_diff_etot[MAX_LT-2] > 10.0f) )
        {
            if ( st->tc_cnt == 1 )
            {
                /* do TC coding instead of GC/VC if onset has been already declared before */
                *sp_aud_decision2 = 0;
                if ( st->codec_mode == MODE1 )
                {
                    *coder_type = TRANSITION;
                }
            }
            else
            {
                if( attack >= ATT_3LSUB_POS )
                {
                    /* do TC coding if attack is located in the last subframe */
                    *sp_aud_decision2 = 0;
                    *attack_flag = 1;
                    if ( st->codec_mode == MODE1 )
                    {
                        *coder_type = TRANSITION;
                    }
                }
                else if( attack >= ATT_SEG_LEN/2 )
                {
                    /* do GSC coding if attack is located after the first quarter of the first subframe */
                    /* (pre-echo will be treated at the decoder side) */
                    *attack_flag = 1;
                }
            }
        }
    }
    else if ( localVAD == 1 && *coder_type == GENERIC &&
              ( (attack >= ATT_3LSUB_POS && st->total_brate < ACELP_24k40) ||
                (attack >= ATT_3LSUB_POS_16k && st->total_brate >= ACELP_24k40 && st->total_brate < ACELP_48k) )
            )
    {
        /* do TC coding an attack is located in the last subframe */
        *attack_flag = 1;
        if ( st->codec_mode == MODE1 )
        {
            *coder_type = TRANSITION;
        }
    }

    return;
}

/*---------------------------------------------------------------------*
 * var_cor_calc()
 *
 * Calculate variance of correlation
 *---------------------------------------------------------------------*/

static void var_cor_calc(
    const float old_corr,
    float *mold_corr,
    float var_cor_t[],
    short *high_stable_cor
)
{
    short i;
    float var_cor;

    /* update buffer of old correlation values */
    for( i = VAR_COR_LEN-1; i > 0; i-- )
    {
        var_cor_t[i] = var_cor_t[i-1];
    }
    var_cor_t[i] = old_corr;

    /* calculate variance of correlation */
    var_cor = var( var_cor_t, VAR_COR_LEN );

    /* set flag in case of highly-correlated stable signal */
    if( *mold_corr > 0.8f && var_cor < 5e-4f )
    {
        *high_stable_cor = 1;
    }
    else
    {
        *high_stable_cor = 0;
    }

    /* update average correlation */
    *mold_corr = 0.1f * old_corr + 0.9f * *mold_corr;

    return;
}

/*---------------------------------------------------------------------*
 * attack_det()
 *
 * Attack detection
 *---------------------------------------------------------------------*/

static short attack_det(
    const float *inp,                 /* i  : input signal                           */
    const short last_clas,            /* i  : last signal clas                       */
    const short localVAD,
    const short coder_type,           /* i  : coder type                             */
    const int   total_brate           /* i  : total bit-rate                         */
)
{
    short i, attack;
    float etmp, etmp2, finc[ATT_NSEG];
    short att_3lsub_pos;

    att_3lsub_pos = ATT_3LSUB_POS;
    if( total_brate >= ACELP_24k40 )
    {
        att_3lsub_pos = ATT_3LSUB_POS_16k;
    }

    /* compute energy per section */
    for( i=0; i<ATT_NSEG; i++ )
    {
        finc[i] = sum2_f( inp + i*ATT_SEG_LEN, ATT_SEG_LEN );
    }

    attack = maximum( finc, ATT_NSEG, &etmp );

    if( localVAD == 1 && coder_type == GENERIC )
    {
        /* compute mean energy in the first three subframes */
        etmp = mean( finc, att_3lsub_pos );

        /* compute mean energy after the attack */
        etmp2 = mean( finc + attack, ATT_NSEG - attack );

        /* and compare them */
        if( etmp * 8 > etmp2 )
        {
            /* stop, if the attack is not sufficiently strong */
            attack = 0;
        }

        if( last_clas == VOICED_CLAS && etmp * 20 > etmp2 )
        {
            /* stop, if the signal was voiced and the attack is not sufficiently strong */
            attack = 0;
        }

        /* compare wrt. other sections (reduces miss-classification) */
        if( attack > 0 )
        {
            etmp2 = finc[attack];

            for( i=2; i<att_3lsub_pos-2; i++ )
            {
                if( finc[i] * 2.0f > etmp2 )
                {
                    /* stop, if the attack is not sufficiently strong */
                    attack = 0;
                    break;
                }
            }
        }
    }

    /* compare wrt. other sections (reduces miss-classification) */
    else if( attack > 0 )
    {
        etmp2 = finc[attack];

        for( i=2; i<att_3lsub_pos-2; i++ )
        {
            if( i != attack && finc[i] * 1.3f > etmp2 )
            {
                /* stop, if the attack is not sufficiently strong */
                attack = 0;
                break;
            }
        }
    }

    return attack;
}


/*------------------------------------------------------------------------*
 * music_mixed_classif_improv()
 *
 * Improve 1st stage speech/music decision for mixed&music signals
 *------------------------------------------------------------------------*/

static void music_mixed_classif_improv(
    Encoder_State *st,                        /* i/o: Encoder state structure                         */
    const float *new_inp,                   /* i  : new input signal                                */
    short *sp_aud_decision1,          /* i/o: improved 1st stage speech/music decision        */
    const short vad_flag,
    const float *voicing,                   /* i  : voicing estimate                                */
    const float *epsP,                      /* i  : LP prediction error                             */
    const float etot,                       /* i  : total frame energy                              */
    const float old_cor,                    /* i  : normalized correlation                          */
    const float cor_map_sum                 /* i  : correlation map sum                             */
)
{
    short i, dec, len, percus_flag;
    float p2v_map[128], ftmp, ftmp1, lt_diff, log_max_spl, epsP_tilt, max_spl;

    /* find sample with maximum absolute amplitude */
    max_spl = 0;
    for ( i=0; i<L_FRAME; i++ )
    {
        if ( fabs(new_inp[i]) > max_spl )
        {
            max_spl = (float)fabs(new_inp[i]);
        }
    }

    /* music is considered only appearing in high SNR condition and active signal */
    if ( vad_flag == 0 || st->lp_speech - st->lp_noise < 25 )
    {
        st->dec_mov = 0.5f;
        st->dec_mov1 = 0.5f;

        if ( vad_flag == 0 )
        {
            st->onset_cnt = 0;
        }

        return;
    }

    st->onset_cnt++;

    if ( st->onset_cnt > 9 )
    {
        st->onset_cnt = 9;
    }

    if ( st->onset_cnt == 1 )
    {
        set_f( st->buf_flux, -100, BUF_LEN );
    }

    /* spectral analysis */
    spec_analysis( st->Bin_E, p2v_map );

    /* percussive music detection */
    log_max_spl = 20 * (float)log(max_spl + 0.0001f);
    lt_diff = log_max_spl - st->mov_log_max_spl;

    for ( i=0; i<3; i++ )
    {
        st->buf_etot[i] = st->buf_etot[i+1];
    }
    st->buf_etot[i] = etot;

    percus_flag = 0;
    if ( st->buf_etot[1] - st->buf_etot[0] > 6 && st->buf_etot[2] < st->buf_etot[1] && st->buf_etot[1] - st->lp_speech > 3 )
    {
        if ( st->buf_etot[1] - st->buf_etot[3] > 3 && st->buf_etot[3] < st->buf_etot[2] && 0.5f * (0.5f * (voicing[0] + voicing[1]) + old_cor) < 0.75f )
        {
            if ( st->dec_mov > 0.8f )
            {
                percus_flag = 1;
            }
            else if ( old_cor < 0.75f && voicing[0] < 0.75f && voicing[1] < 0.75f && st->old_lt_diff[0] > 10 )
            {
                percus_flag = 1;
            }
        }
    }

    /* sound attack detection */
    if ( st->buf_etot[3] - st->buf_etot[2] > 6 && st->dec_mov > 0.9f && etot - st->lp_speech > 5 && st->old_lt_diff[0] > 5 )
    {
        st->attack_hangover = 3;
    }

    if ( voicing[0] > 0.9f && voicing[1] > 0.9f )
    {
        if ( log_max_spl > st->mov_log_max_spl )
        {
            st->mov_log_max_spl = 0.75f * st->mov_log_max_spl + (1 - 0.75f) * log_max_spl;
        }
        else
        {
            st->mov_log_max_spl = 0.995f * st->mov_log_max_spl + (1 - 0.995f) * log_max_spl;
        }
    }

    st->old_lt_diff[0] = st->old_lt_diff[1];
    st->old_lt_diff[1] = lt_diff;

    /* calculate and buffer spectral energy fluctuation */
    flux( st->Bin_E, p2v_map, st->old_Bin_E, st->buf_flux, st->attack_hangover, st->dec_mov );

    st->attack_hangover--;
    if ( st->attack_hangover < 0 )
    {
        st->attack_hangover = 0;
    }

    /* identify flux buffer status */
    len = 0;
    for ( i = BUF_LEN-1; i >= 0 && st->buf_flux[i] >= 0; i-- )
    {
        len++;
    }

    /* reset flux buffer if percussive music is detected */
    if ( percus_flag == 1 )
    {
        set_f( &st->buf_flux[BUF_LEN-len], 5, len );
    }

    /* calculate and buffer the tilt of residual LP analysis energies */
    ftmp = 0.00001f;
    ftmp1 = 0;
    for ( i=1; i<16; i++ )
    {
        ftmp += epsP[i] * epsP[i];
        ftmp1 += epsP[i] * epsP[i+1];
    }

    epsP_tilt = ftmp1/ftmp;

    for ( i=0; i<BUF_LEN-1; i++ )
    {
        st->buf_epsP_tilt[i] = st->buf_epsP_tilt[i+1];
    }
    st->buf_epsP_tilt[i] = epsP_tilt;

    /* calculate and buffer highband spectral peakness */
    tonal_dist( p2v_map, st->buf_pkh, st->buf_Ntonal, st->buf_Ntonal2, st->buf_Ntonal_lf );

    /* buffer sum of correlation map */
    for ( i=0; i<BUF_LEN-1; i++ )
    {
        st->buf_cor_map_sum[i] = st->buf_cor_map_sum[i+1];
    }
    st->buf_cor_map_sum[i] = cor_map_sum;

    /* buffer voicing metric */
    for ( i=0; i<9; i++ )
    {
        st->buf_dlp[i] = st->buf_dlp[i+1];
    }
    st->buf_dlp[i] = st->lps - st->lpm;

    /* classification */
    dec = mode_decision( st, len, &st->dec_mov, st->buf_flux, st->buf_epsP_tilt, st->buf_pkh, st->buf_cor_map_sum,
                         st->buf_Ntonal, st->buf_Ntonal2, st->buf_Ntonal_lf, st->buf_dlp );

    /* update long term moving average of the classification decisions */
    if ( len > 30 )
    {
        st->dec_mov = 0.97f * st->dec_mov + (1 - 0.97f) * dec;
        st->dec_mov1 = 0.97f * st->dec_mov1 + (1 - 0.97f) * dec;
    }

    /* update long-term unvoiced counter */
    if ( (st->coder_type_raw == UNVOICED || st->coder_type_raw == INACTIVE) && etot > 1.5f && st->buf_Ntonal2[59] < 2 )
    {
        st->UV_cnt1 -= 8;
    }
    else
    {
        st->UV_cnt1++;
    }

    if ( st->UV_cnt1 > 300 )
    {
        st->UV_cnt1 = 300;
    }
    else if ( st->UV_cnt1 < 0 )
    {
        st->UV_cnt1 = 0;
    }

    st->LT_UV_cnt1 = 0.9f * st->LT_UV_cnt1 + 0.1f * st->UV_cnt1;

    /* revert classification decision due to long-term unvoiced counter */
    if ( dec == 1 && st->dec_mov1 < 0.2f && st->LT_UV_cnt1 < 200 )
    {
        dec = 0;
    }

    /* overwrite 1st stage speech/music decision to music */
    if ( dec == 1 )
    {
        *sp_aud_decision1 = 1;
    }

    return;
}


/*---------------------------------------------------------------------*
 * spec_analysis()
 *
 * Spectral analysis for mixed/music classification improvement
 *---------------------------------------------------------------------*/

static void spec_analysis(
    float *Bin_E,                     /* i  : log energy spectrum of the current frame        */
    float *p2v_map                    /* o  : spectral peakiness map                          */
)
{
    short i, k, m;
    float peak[L_FFT/4+1];
    float valley[L_FFT/4+1];
    short peak_idx[L_FFT/4+1];
    short valey_idx[L_FFT/4+1];
    float p2v[L_FFT/4+1];

    /* find spectral peaks */
    k = 0;
    for ( i=1; i<L_FFT/2-2; i++ )
    {
        if ( Bin_E[i] > Bin_E[i-1] && Bin_E[i] > Bin_E[i+1] )
        {
            peak[k] = Bin_E[i];
            peak_idx[k] = i;
            k++;
        }
    }
    assert(k+1<L_FFT/4+1);
    peak_idx[k] = -1;
    peak_idx[k+1] = -1;

    if ( k == 0 )
    {
        for ( i=0; i<L_FFT/2-1; i++ )
        {
            p2v_map[i] = 0;
        }

        return;
    }

    /* find spectral valleys */
    m = 0;
    if ( Bin_E[0] < Bin_E[1] )
    {
        valley[0] = Bin_E[0];
        valey_idx[0] = 0;
        m++;
    }

    k = L_FFT/2-2;
    for ( i=L_FFT/2-3; i >= 0 && Bin_E[i+1] > Bin_E[i]; i-- )
    {
        k = i;
    }

    for ( i=1; i<k; i++ )
    {
        if ( Bin_E[i] < Bin_E[i-1] && Bin_E[i] < Bin_E[i+1] )
        {
            valley[m] = Bin_E[i];
            valey_idx[m] = i;
            m++;
        }
    }

    valley[m] = Bin_E[k];
    valey_idx[m] = k;

    /* find spectral peak to valley distances */
    k = 0;
    for (i=0; i<m; i++)
    {
        if ( peak_idx[k] > valey_idx[i] && peak_idx[k] < valey_idx[i+1] )
        {
            p2v[k] = 2*peak[k] - valley[i] - valley[i+1];
            k++;
        }
    }

    for ( i=0; i<L_FFT/2-1; i++ )
    {
        p2v_map[i] = 0;
    }

    for ( i=0; i<k; i++ )
    {
        p2v_map[peak_idx[i]] = p2v[i];
    }

    return;
}

/*---------------------------------------------------------------------*
 * flux()
 *
 * Calculation of spectral flux
 *---------------------------------------------------------------------*/

static void flux(
    float *Bin_E,                     /* i  : log energy spectrum of the current frame        */
    float *p2v_map,                   /* i  : spectral peakiness map                          */
    float *old_Bin_E,                 /* i/o: log energy spectrum of the frame 60ms ago       */
    float *buf_flux,                  /* i/o: buffer storing spectral energy fluctuation      */
    short attack_hangover,            /* i/o: hangover preventing flux buffering              */
    float dec_mov                     /* i/o: moving average of classifier decision           */
)
{
    short i;
    float *pt1,*pt2,*pt3,*pt4,*pt5,*pt6;
    float flux;
    short cnt;

    /* calculate flux */
    flux = 0 ;
    cnt = 0;
    for ( i=0; i<N_OLD_BIN_E; i++ )
    {
        if ( p2v_map[i] != 0 )
        {
            flux += (float)fabs( Bin_E[i] - old_Bin_E[i] );
            cnt++;
        }
    }

    if ( cnt == 0 )
    {
        flux = 5;
    }
    else
    {
        flux = flux/(float)cnt;
    }

    if ( flux > 20 && dec_mov > 0.8f )
    {
        flux = 20;
    }

    /* update old Bin_E buffer */
    pt1 = old_Bin_E;
    pt2 = old_Bin_E + N_OLD_BIN_E;
    pt3 = Bin_E;
    pt4 = old_Bin_E + N_OLD_BIN_E;
    pt5 = old_Bin_E + 2*N_OLD_BIN_E;
    pt6 = old_Bin_E + 2*N_OLD_BIN_E;

    for ( i=0; i<N_OLD_BIN_E; i++ )
    {
        *pt1++ = *pt2++;
        *pt4++ = *pt5++;
        *pt6++ = *pt3++;
    }
    /* update flux buffer */
    if ( attack_hangover <= 0 )
    {
        for ( i=0; i<BUF_LEN-1; i++ )
        {
            buf_flux[i] = buf_flux[i+1];
        }

        buf_flux[i] = flux;
    }

    return;
}

/*---------------------------------------------------------------------*
 * tonal_dist()
 *
 * Calculation of spectral distance
 *---------------------------------------------------------------------*/

static void tonal_dist(
    float *p2v_map,                   /* i  : spectral peakiness map                          */
    float *buf_pkh,                   /* i/o: buffer storing highband spectral peakiness      */
    float *buf_Ntonal,                /* i/o: buffer storing No.of 1st spectral tone          */
    float *buf_Ntonal2,               /* i/o: buffer storing No.of 2nd spectral tone          */
    float *buf_Ntonal_lf              /* i/o: buffer storing low band spectral tone ratio     */
)
{
    short i;
    float pk;
    short Ntonal;
    short Ntonal2;
    short Ntonal_lf;

    /* find number of tonals, number of tonals at low-band,
    spectral peakiness at high-band */
    pk = 0;
    Ntonal = 0;
    Ntonal2 = 0;
    Ntonal_lf = 0;
    for ( i=0; i<64; i++ )
    {
        if ( p2v_map[i] > 55 )
        {
            Ntonal++;
        }

        if ( p2v_map[i] > 80 )
        {
            Ntonal2++;
            Ntonal_lf++;
        }
    }

    for ( i=64; i<127; i++ )
    {
        if ( p2v_map[i] != 0 )
        {
            pk += p2v_map[i];
        }

        if ( p2v_map[i] > 55 )
        {
            Ntonal++;
        }

        if ( p2v_map[i] > 80 )
        {
            Ntonal2++;
        }
    }

    /* update buffers */
    for ( i=0; i<BUF_LEN-1; i++ )
    {
        buf_pkh[i] = buf_pkh[i+1];
        buf_Ntonal[i] = buf_Ntonal[i+1];
        buf_Ntonal2[i] = buf_Ntonal2[i+1];
        buf_Ntonal_lf[i] = buf_Ntonal_lf[i+1];
    }

    buf_pkh[i] = pk;
    buf_Ntonal[i] = (float)Ntonal;
    buf_Ntonal2[i] = (float)Ntonal2;
    buf_Ntonal_lf[i] = (float)Ntonal_lf;
}

/*---------------------------------------------------------------------*
 * mode_decision()
 *
 * Decision about internal mode of the mixed/music classifier improvement
 *---------------------------------------------------------------------*/

static short mode_decision(
    Encoder_State *st,
    short len,                         /* i  : buffering status                                */
    float *dec_mov,                    /* i/o: moving average of classifier decision           */
    float *buf_flux,                   /* i  : buffer storing spectral energy fluctuation      */
    float *buf_epsP_tilt,              /* i  : buffer storing LP prediciton error tilt         */
    float *buf_pkh,                    /* i  : buffer storing highband spectral peakiness      */
    float *buf_cor_map_sum,            /* i  : buffer storing correlation map sum              */
    float *buf_Ntonal,                 /* i  : buffer storing No.of 1st spectral tone          */
    float *buf_Ntonal2,                /* i  : buffer storing No.of 2nd spectral tone          */
    float *buf_Ntonal_lf,              /* i  : buffer storing low band spectral tone ratio     */
    float *buf_dlp                     /* i  : buffer storing voicing estimate                 */
)
{
    short mode;
    short i;
    short voiced_cnt;
    float M_pkh;
    float M_cor_map_sum;
    float M_Ntonal;
    float M_flux;
    float V_epsP_tilt;
    float lf_Ntonal_ratio;

    mode = *dec_mov > 0.5f;

    if ( len <= 5 )
    {
        return ( mode );
    }
    else if ( len < 10 )
    {
        M_pkh = mean( buf_pkh+BUF_LEN-len, len );
        M_cor_map_sum = mean( buf_cor_map_sum+BUF_LEN-len, len );
        M_Ntonal = mean( buf_Ntonal+BUF_LEN-len, len );
        V_epsP_tilt = var( buf_epsP_tilt+BUF_LEN-len, len );

        voiced_cnt = 0;
        for ( i=9; i>3; i-- )
        {
            if ( buf_dlp[i] > 0.0f )
            {
                voiced_cnt++;
            }
        }

        if ( (M_pkh > 1100 || V_epsP_tilt < 0.00008f || M_cor_map_sum > 100) && voiced_cnt < 4 )
        {
            mode = 1;
        }
        else if ( M_Ntonal > 27 && voiced_cnt < 4 )
        {
            mode = 1;
        }
    }
    else
    {
        voiced_cnt = 0;
        for ( i=0; i<10; i++ )
        {
            if ( buf_dlp[i] > 0.0f )
            {
                voiced_cnt++;
            }
        }

        M_flux = mean( &buf_flux[BUF_LEN-10], 10 );
        M_pkh = mean( buf_pkh+BUF_LEN-10, 10 );
        M_cor_map_sum = mean( buf_cor_map_sum+BUF_LEN-10, 10 );
        V_epsP_tilt = var( buf_epsP_tilt+BUF_LEN-10, 10 );

        if ( (M_flux < 8.5f || (V_epsP_tilt < 0.001f && M_flux < 12.0f) || M_pkh > 1050 || M_cor_map_sum > 100)
                && voiced_cnt < 3 && mean( &buf_flux[55], 5 ) < 15 )
        {
            mode = 1;
            *dec_mov = 1;
            return ( mode );
        }

        if ( M_flux > 16.0f || (M_flux > 15 && voiced_cnt > 2) || mean( &buf_flux[55], 5 ) > 19.0f || (buf_flux[59] >= 20 && st->lps-st->lpm > 0) )
        {
            *dec_mov = 0;
            mode = 0;
            return ( mode );
        }

        for ( i=10; i<len; i++ )
        {
            M_flux = mean( &buf_flux[BUF_LEN-i], i );
            M_pkh = mean( buf_pkh+BUF_LEN-i, i );
            M_cor_map_sum = mean( buf_cor_map_sum+BUF_LEN-i, i );
            V_epsP_tilt = var( buf_epsP_tilt+BUF_LEN-i, i );

            if ( ((M_flux < 12+0.05f*(len-10) && mean( &buf_flux[BUF_LEN-10], 10 ) < 15) || V_epsP_tilt < 0.0001f+0.000018f*(len-10)
                    || M_pkh > 1050-5.0f*(len-10) || M_cor_map_sum > 95-0.3f*(len-10)) && voiced_cnt < 3 )
            {
                mode = 1;
                return( mode );
            }
        }

        if ( len == BUF_LEN )
        {
            M_Ntonal = mean( buf_Ntonal, BUF_LEN );
            lf_Ntonal_ratio = sum_f(buf_Ntonal_lf, BUF_LEN)/(sum_f(buf_Ntonal2, BUF_LEN) + 0.0001f);

            if ( M_Ntonal > 18 || lf_Ntonal_ratio < 0.2f )
            {
                mode = 1;
            }
            else if ( M_Ntonal < 1 )
            {
                mode = 0;
            }
        }
    }

    return( mode );
}

/*----------------------------------------------------------------------------------*
 * tonal_context_improv()
 *
 * Context-based improvement of 1st/2nd stage speech/music decision on stable tonal signals
 *----------------------------------------------------------------------------------*/

static void tonal_context_improv(
    Encoder_State *st,               /* i/o: encoder state structure                       */
    const float PS[],              /* i  : energy spectrum                               */
    short *sp_aud_decision1, /* i/o: 1st stage speech/music decision               */
    short *sp_aud_decision2, /* i/o: 2nd stage speech/music decision               */
    const short vad_flag,
    const short pitch[3],          /* i  : open-loop pitch estimate in three subframes   */
    const float voicing[3],        /* i  : voicing estimate in three subframes           */
    const float voi_fv,            /* i  : scaled voicing feature                        */
    const float cor_map_sum_fv,    /* i  : scaled correlation map feature                */
    const float LPCErr             /* i  : scaled LP prediction error feature            */
)
{
    short lt_pitch_diff;
    float sort_max, sort_avg, sort_val[80];
    float tonality, tonality1, tonality2, tonality3, t2, t3, tL, err, cor, dft;

    /* reset in case of codec mode swithing */
    if( st->last_codec_mode == MODE2 )
    {
        set_f( st->tonality2_buf, 0, HANG_LEN_INIT );
        set_f( st->tonality3_buf, 0, HANG_LEN_INIT );
        set_f( st->LPCErr_buf, 0, HANG_LEN_INIT );
        st->lt_music_hangover = 0;
        st->lt_music_state    = 0;
        st->lt_speech_state   = 0;
        st->lt_speech_hangover= 0;
    }

    /* estimate maximum tonality in bands [0-1 kHz], [1-2kHz] and [2-4kHz] */
    mvr2r( PS, sort_val, 80 );

    /* tonality in band 0-1 kHz */
    v_sort(sort_val, 0, 19);
    sort_max = sort_val[19];
    sort_avg = sum_f(&sort_val[0], 10);
    tonality1 = sort_max / sort_avg;

    /* tonality in band 1-2 kHz */
    v_sort(sort_val, 20, 39);
    sort_max = sort_val[39];
    sort_avg = sum_f(&sort_val[20], 10);
    tonality2 = sort_max / sort_avg;

    /* tonality in band 2-4 kHz */
    v_sort(sort_val, 40, 79);
    sort_max = sort_val[79];
    sort_avg = sum_f(&sort_val[40], 20);
    tonality3 = sort_max / sort_avg;

    tonality = max(max(tonality1, tonality2), tonality3);

    if( st->hangover_cnt == 10 && vad_flag == 1 )
    {
        /* long-term voicing parameter */
        st->lt_voicing = 0.1f * st->lt_voicing + 0.9f * *voicing;

        /* long-term correlation value */
        st->lt_corr = 0.1f * st->lt_corr + 0.9f * st->old_corr;

        /* long-term tonality measure */
        st->lt_tonality = 0.1f * st->lt_tonality + 0.9f * tonality;
    }
    else
    {
        /* long-term voicing parameter */
        st->lt_voicing = 0.7f * st->lt_voicing + 0.3f * *voicing;

        /* long-term correlation value */
        st->lt_corr = 0.7f * st->lt_corr + 0.3f * st->old_corr;

        /* long-term tonality measure */
        st->lt_tonality = 0.5f * st->lt_tonality + 0.5f * tonality;
    }

    /* pitch difference w.r.t to past 3 frames */
    lt_pitch_diff =  (short)abs(st->lt_corr_pitch[0] - pitch[0]);
    lt_pitch_diff += (short)abs(st->lt_corr_pitch[1] - pitch[0]);
    lt_pitch_diff += (short)abs(st->lt_corr_pitch[2] - pitch[0]);

    st->lt_corr_pitch[0] = st->lt_corr_pitch[1];
    st->lt_corr_pitch[1] = st->lt_corr_pitch[2];
    st->lt_corr_pitch[2] = pitch[0];

    st->lt_old_mode[0] = st->lt_old_mode[1];
    st->lt_old_mode[1] = st->lt_old_mode[2];

    if ( *sp_aud_decision1 == 1 &&
            ( min(min(tonality1, tonality2), tonality3) > 50.0f ) &&
            ( tonality1+tonality2 > 200.0f && tonality2 + tonality3 > 200.0f && tonality1 + tonality3 > 200.0f ) &&
            ( st->lt_tonality < 20000.0f ) &&
            ( ( st->lt_tonality > 1000 && max(st->lt_voicing, *voicing) > 0.99f ) ||
              ( st->lt_tonality > 1500 && st->lt_corr > 0.99f ) ||
              ( st->lt_tonality > 3000 && st->lowrate_pitchGain > 0.96f ) ||
              ( lt_pitch_diff == 0 && st->lowrate_pitchGain > 0.89f ) ) )
    {
        if( sum_s(st->lt_old_mode, 2) < 2 )
        {
            /* probably speech - change the decision to speech */
            *sp_aud_decision1 = 0;
            *sp_aud_decision2 = 0;

            if( st->lt_hangover == 0 )
            {
                st->lt_hangover = 6;
            }
        }
    }
    else
    {
        /* not speech, but still in the hangover period - change the decision to speech */
        if( st->lt_hangover > 0 )
        {
            *sp_aud_decision1 = 0;
            *sp_aud_decision2 = 0;
            st->lt_hangover--;
        }
    }

    /* calculate standard deviation of log-tonality */
    mvr2r( st->tonality2_buf + 1, st->tonality2_buf, HANG_LEN_INIT - 1 );
    st->tonality2_buf[HANG_LEN_INIT - 1] = 0.2f*(float)log10(tonality2);
    t2 = std_dev( st->tonality2_buf, HANG_LEN_INIT );

    mvr2r( st->tonality3_buf + 1, st->tonality3_buf, HANG_LEN_INIT - 1 );
    st->tonality3_buf[HANG_LEN_INIT - 1] = 0.2f*(float)log10(tonality3);
    t3 = std_dev( st->tonality3_buf, HANG_LEN_INIT );

    tL  = 0.2f*(float)log10(st->lt_tonality);

    /* calculate standard deviation of residual LP energy */
    mvr2r( st->LPCErr_buf + 1, st->LPCErr_buf, HANG_LEN_INIT - 1 );
    st->LPCErr_buf[HANG_LEN_INIT - 1] = LPCErr;
    err = std_dev( st->LPCErr_buf, HANG_LEN_INIT );

    cor = max( voi_fv - cor_map_sum_fv, 0.0f );
    dft = 0.2f * (float) fabs( log10(tonality2) - log10(tonality3) );

    /* state machine for strong music */
    if( *sp_aud_decision1 == 1 && st->lt_music_state == 0 && st->lt_music_hangover == 0 &&
            t2 < 0.54f && t2 > 0.26f && t3 > 0.22f && tL < 0.54f && tL > 0.26f && err > 0.5f )
    {
        st->lt_music_state = 1;
        st->lt_music_hangover = 6;
    }
    else if( st->lt_music_state == 1 && st->lt_music_hangover == 0 &&
             t2 < 0.34 && t3 < 0.26f && tL < 0.45f)
    {
        st->lt_music_state = 0;
        st->lt_music_hangover = 6;
    }

    if( st->lt_music_hangover > 0 )
    {
        st->lt_music_hangover--;
    }

    /* state machine for strong speech */
    if( *sp_aud_decision1 == 1 && st->lt_speech_state == 0 && st->lt_speech_hangover == 0 &&
            cor > 0.40f && dft < 0.1f && voi_fv > 2 * cor_map_sum_fv + 0.12f &&
            t2 < cor && t3 < cor && tL < cor && cor_map_sum_fv < cor && voi_fv > cor && voi_fv > 0.76f )
    {
        st->lt_speech_state = 1;
        st->lt_speech_hangover = 6;
    }
    else if( st->lt_speech_state == 1 && st->lt_speech_hangover == 0 && cor < 0.40f )
    {
        st->lt_speech_state = 0;
        st->lt_speech_hangover = 6;
    }

    if( st->lt_speech_hangover > 0 )
    {
        st->lt_speech_hangover--;
    }

    /* final decision */
    if ( *sp_aud_decision1 == 1 && st->lt_speech_state == 1 )
    {
        /* strong speech - probably error in speech/music classification */
        *sp_aud_decision1 = 0;
        *sp_aud_decision2 = 0;
    }
    else if ( *sp_aud_decision1 == 0 && st->lt_music_state == 1 )
    {
        /* strong music - probably error in speech/music classification */
        *sp_aud_decision1 = 1;
        *sp_aud_decision2 = 1;
    }

    /* update the buffer of past decisions */
    st->lt_old_mode[2] = *sp_aud_decision1;

    return;
}

/*---------------------------------------------------------------------*
 * detect_sparseness()
 *
 *
 *---------------------------------------------------------------------*/

static void detect_sparseness(
    Encoder_State *st,               /* i/o: encoder state structure                */
    const short localVAD_HE_SAD,   /* i  : HE-SAD flag without hangover           */
    short *sp_aud_decision1, /* i/o: 1st stage speech/music decision        */
    short *sp_aud_decision2, /* i/o: 2nd stage speech/music decision        */
    const float voi_fv             /* i  : scaled voicing feature                 */
)
{
    float sum;
    float ftmp;
    float ftmp1;
    float S1[128];
    short i,j;
    short hb_sp_high_flag = 0;
    short lb_sp_high_flag = 0;
    float sumh;
    float sparse;
    float tmp_buf[4];
    float Mlpe = 0.0f;
    float Mv = 0.0f;
    float Msp;

    mvr2r( st->Bin_E, S1, 128 );

    sum = 0;
    for ( i=0; i<80; i++ )
    {
        if ( S1[i] < 0 )
        {
            S1[i] = 0;
        }
        sum += S1[i];
    }

    sumh = 0;
    for ( i=80; i<128; i++ )
    {
        if ( S1[i] < 0 )
        {
            S1[i] = 0;
        }
        sumh += S1[i];
    }

    sum += sumh;

    /* order spectral from max to min */
    order_spectrum(S1, 128);

    /* calculate spectral sparseness in the range 0 - 6.4 kHz */
    j = 0;
    ftmp = 0.0f;
    ftmp1 = 0.75f * sum;
    for ( i=0; i<128; i++ )
    {
        ftmp += S1[i];
        if ( ftmp > ftmp1 )
        {
            j = i;
            break;
        }
    }

    for ( i=0; i<HANG_LEN_INIT-1; i++ )
    {
        st->sparse_buf[i] = st->sparse_buf[i+1];
    }

    sparse = (float)j;
    st->sparse_buf[i] = sparse;

    if ( st->bwidth == WB )
    {
        Msp = mean(st->sparse_buf, 8);

        /* find long-term smoothed sparseness */
        if ( st->last_vad_spa == 0 )
        {
            set_f( &st->sparse_buf[0], sparse, HANG_LEN_INIT-1 );
            st->LT_sparse = sparse;
        }
        else
        {
            set_f(tmp_buf, 0.0f, 4);

            for ( i=0; i<HANG_LEN_INIT; i++ )
            {
                for ( j=0; j<4; j++ )
                {
                    if ( st->sparse_buf[i] > tmp_buf[j] )
                    {
                        mvr2r(&tmp_buf[j], &tmp_buf[j+1], 3-j);
                        tmp_buf[j] = st->sparse_buf[i];
                        break;
                    }
                }
            }

            ftmp = 0.25f*(HANG_LEN_INIT*Msp - sum_f(tmp_buf, 4)) - st->LT_sparse;

            st->LT_sparse = st->LT_sparse + 0.25f * ftmp;
        }

        /* find high-band sparseness */
        mvr2r(st->Bin_E+80, S1, 48);
        order_spectrum(S1, 48);

        for ( i=0; i<HANG_LEN_INIT-1; i++ )
        {
            st->hf_spar_buf[i] = st->hf_spar_buf[i+1];
        }
        st->hf_spar_buf[i] = sum_f(S1, 5)/(sumh + 0.1f);
        if ( mean(st->hf_spar_buf, 8) > 0.2f )
        {
            hb_sp_high_flag = 1;
        }

        /* find low-band sparseness */
        mvr2r(st->Bin_E, S1, 60);
        order_spectrum(S1, 60);

        if ( sum_f(S1, 5)/sum_f(S1,60) > 0.18f )
        {
            lb_sp_high_flag = 1;
        }

        /* find smoothed linear prediction efficiency */
        for ( i=0; i<7; i++ )
        {
            st->lpe_buf[i] = st->lpe_buf[i+1];
        }

        st->lpe_buf[i] = st->past_epsP2;
        Mlpe = mean(st->lpe_buf, 8);

        /* find smoothed voicing */
        for ( i=0; i<HANG_LEN_INIT-1; i++ )
        {
            st->voicing_buf[i] = st->voicing_buf[i+1];
        }

        st->voicing_buf[i] = voi_fv;
        Mv = mean(st->voicing_buf, 8);
    }

    /* avoid using LR-MDCT on sparse spectra */
    if ( *sp_aud_decision1 == 1 )
    {
        if ( st->bwidth == WB )
        {
            ftmp = 90;
        }
        else
        {
            ftmp = 91;
        }
        if ( sparse > ftmp )
        {
            *sp_aud_decision1 = 0;
            *sp_aud_decision2 = 1;
            st->gsc_hangover = 1;
        }
        else if ( st->gsc_hangover == 1 )
        {
            if ( sparse > 85 )
            {
                *sp_aud_decision1 = 0;
                *sp_aud_decision2 = 1;
            }
            else if ( fabs(sparse - mean(&st->sparse_buf[HANG_LEN_INIT-1-st->gsc_cnt], st->gsc_cnt)) < 7.0f )
            {
                *sp_aud_decision1 = 0;
                *sp_aud_decision2 = 1;
            }
        }

        if ( st->bwidth == WB )
        {
            if ( st->LT_sparse > 60 && sparse > 50 && Mlpe < -1.3f && Mv > 0.85f &&
                    lb_sp_high_flag == 0 && ( (hb_sp_high_flag == 0 && sumh > 0.15f * sum) || sumh <= 0.15f * sum ) )
            {
                *sp_aud_decision1 = 0;
                *sp_aud_decision2 = 1;
                st->gsc_hangover = 1;
            }
            else if ( st->gsc_hangover == 1 && !( *sp_aud_decision1 == 0 && *sp_aud_decision2 == 1) )
            {
                if ( fabs(sparse - mean(&st->sparse_buf[HANG_LEN_INIT-1-st->gsc_cnt], st->gsc_cnt)) < 7.0f )
                {
                    *sp_aud_decision1 = 0;
                    *sp_aud_decision2 = 1;
                }
            }
        }
    }

    /* update the counter of consecutive GSC frames with sparse spectrum */
    if ( *sp_aud_decision1 == 0 && *sp_aud_decision2 == 1 )
    {
        (st->gsc_cnt)++;
        if ( st->gsc_cnt > 7 )
        {
            st->gsc_cnt = 7;
        }
    }
    else
    {
        st->gsc_cnt = 0;
        st->gsc_hangover = 0;
    }

    st->last_vad_spa = localVAD_HE_SAD;

    return;
}


/*---------------------------------------------------------------------*
 * order_spectrum()
 *
 *
 *---------------------------------------------------------------------*/

static void order_spectrum(
    float *vec,
    short len
)
{
    short i, j;
    float *pts,*pte;
    float max, min;
    short imax, imin;

    pts = &vec[0];
    pte = &vec[len-1];
    for ( i=0; i<len/2; i++ )
    {
        max = vec[i];
        min = vec[i];
        imax = i;
        imin = i;
        for ( j=i; j<len-i; j++ )
        {
            if ( vec[j] > max )
            {
                max = vec[j];
                imax = j;
            }
            else
            {
                if ( vec[j] < min )
                {
                    min = vec[j];
                    imin = j;
                }
            }
        }

        vec[imax] = *pts;
        vec[imin] = *pte;
        *pts++ = max;
        *pte-- = min;
    }

    return;
}

