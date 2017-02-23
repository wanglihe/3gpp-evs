/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "rom_enc.h"
#include "rom_com.h"
#include "prot.h"

/*-------------------------------------------------------------------*
 * pre_proc()
 *
 * Pre-processing (spectral analysis, LP analysis, VAD, OL pitch calculation, coder mode selection, ...)
 *--------------------------------------------------------------------*/

void pre_proc(
    Encoder_State *st,                      /* i/o: encoder state structure                  */
    const short input_frame,              /* i  : frame length                             */
    const float signal_in[],              /* i  : new samples                              */
    float old_inp_12k8[],           /* i/o: buffer of old input signal               */
    float old_inp_16k[],            /* i/o: buffer of old input signal @ 16kHz       */
    float **inp,                    /* o  : ptr. to inp. signal in the current frame */
    short *sp_aud_decision1,        /* o  : 1st stage speech/music classification    */
    short *sp_aud_decision2,        /* o  : 2nd stage speech/music classification    */
    float fr_bands[2*NB_BANDS],     /* o  : energy in frequency bands                */
    short *vad_flag,
    short *localVAD,
    float *Etot,                    /* o  : total energy                             */
    float *ener,                    /* o  : residual energy from Levinson-Durbin     */
    short pitch[3],                 /* o  : open-loop pitch values for quantiz.      */
    float voicing[3],               /* o  : OL maximum normalized correlation        */
    float A[NB_SUBFR16k*(M+1)],     /* o  : A(z) unquantized for the 4 subframes     */
    float Aw[NB_SUBFR16k*(M+1)],    /* o  : weighted A(z) unquantized for subframes  */
    float epsP[M+1],                /* o  : LP prediction errors                     */
    float lsp_new[M],               /* o  : LSPs at the end of the frame             */
    float lsp_mid[M],               /* o  : LSPs in the middle of the frame          */
    short *coder_type,              /* o  : coder type                               */
    short *sharpFlag,               /* o  : formant sharpening flag                  */
    short *vad_hover_flag,
    short *attack_flag,             /* o  : flag signalling attack encoded by AC mode (GSC)    */
    float *new_inp_resamp16k,       /* o  : new input signal @16kHz, non pre-emphasised, used by the WB TBE/BWE */
    short *Voicing_flag,            /* o  : voicing flag for HQ FEC                  */
    float realBuffer[CLDFB_NO_COL_MAX][CLDFB_NO_CHANNELS_MAX], /* i/o : real buffer */
    float imagBuffer[CLDFB_NO_COL_MAX][CLDFB_NO_CHANNELS_MAX], /* i/o : imag buffer */
    short *hq_core_type           /* o  : HQ core type                             */
)
{
    short delay;
    float *inp_12k8, *new_inp_12k8, *inp_16k, *new_inp_16k;   /* pointers to current frame and new data */
    float old_wsp[L_WSP], *wsp;                               /* weighted input signal buffer         */
    float pitch_fr[NB_SUBFR];                                 /* fractional pitch values */
    float voicing_fr[NB_SUBFR];                               /* fractional pitch gains               */
    float lf_E[2*VOIC_BINS];                                  /* per bin spectrum energy in lf        */
    float tmpN[NB_BANDS];                                     /* Temporary noise update               */
    float tmpE[NB_BANDS];                                     /* Temporary averaged energy of 2 sf.   */
    float ee[2];                                              /* Spectral tilt                        */
    float corr_shift;                                         /* correlation shift                    */
    float relE;                                               /* frame relative energy                */
    short loc_harm;                                           /* harmonicity flag                     */
    float cor_map_sum, sp_div, PS[128];                       /* speech/music clasif. parameters      */
    short L_look;                                             /* length of look-ahead                 */
    float snr_sum_he;                                         /* HE SAD parameters                    */
    short localVAD_HE_SAD;                                    /* HE SAD parameters                    */
    short vad_flag_dtx;                                       /* HE-SAD flag with additional DTX HO   */
    short vad_flag_cldfb;
    float old_cor;
    short uc_clas;
    float hp_E[2];                                            /* Energy in HF                         */
    short noisy_speech_HO, clean_speech_HO, NB_speech_HO;     /* SC-VBR HO flags                      */
    float non_staX;                                           /* unbound non-stationarity for sp/mus clas. */
    short sr_core_tmp;
    short L_frame_tmp;
    short flag_spitch;
    float lsf_new[M], stab_fac;
    float band_energies[2*NB_BANDS];                          /* energy in critical bands without minimum noise floor E_MIN */
    float enerBuffer[CLDFB_NO_CHANNELS_MAX];
    float currFlatness;
    short high_lpn_flag;
    short cldfb_addition = 0;
    short alw_pitch_lag_12k8[2];
    float alw_voicing[2];
    float fft_buff[2*L_FFT];
    float sp_floor;
    short sp_aud_decision0;
    short last_core_orig;



    /*------------------------------------------------------------------*
     * Initialization
     *------------------------------------------------------------------*/

    localVAD_HE_SAD = 0;
    NB_speech_HO = 0;
    clean_speech_HO = 0;
    noisy_speech_HO = 0;
    snr_sum_he = 0;
    currFlatness = 0;

    *vad_hover_flag = 0;
    uc_clas = VOICED_CLAS;
    *sp_aud_decision1 = 0;
    *sp_aud_decision2 = 0;
    *coder_type = GENERIC;
    st->noise_lev = NOISE_LEVEL_SP0;
    *attack_flag = 0;

    st->bump_up = 0;
    st->ppp_mode = 0;
    st->nelp_mode = 0;
    st->avoid_HQ_VBR_NB = 0;

    L_look = L_LOOK_12k8;                                       /* lookahead at 12.8kHz */

    new_inp_12k8 = old_inp_12k8 + L_INP_MEM;                    /* pointer to new samples of the input signal in 12.8kHz core */
    inp_12k8 = new_inp_12k8 - L_look;                           /* pointer to the current frame of input signal in 12.8kHz core */
    mvr2r( st->old_inp_12k8, old_inp_12k8, L_INP_MEM );

    mvr2r( st->old_wsp, old_wsp, L_WSP_MEM );
    wsp = old_wsp + L_WSP_MEM;                                  /* pointer to the current frame of weighted signal in 12.8kHz core */

    old_cor = st->old_corr;                                     /* save old_cor for speech/music classifier */

    st->rf_mode = st->Opt_RF_ON;

    last_core_orig = st->last_core;

    /*--------------------------------------------------------------*
     * Cldfb analysis
     *---------------------------------------------------------------*/

    st->prevEnergyHF = st->currEnergyHF;

    analysisCldfbEncoder( st, signal_in, input_frame, realBuffer, imagBuffer, enerBuffer );

    /*----------------------------------------------------------------*
     * Change the sampling frequency to 12.8 kHz
     *----------------------------------------------------------------*/

    modify_Fs( signal_in, input_frame, st->input_Fs, new_inp_12k8, 12800, st->mem_decim,(st->max_bwidth == NB) );

    /* save input resampled at 12.8kHz, non-preemhasised */
    mvr2r( new_inp_12k8, st->buf_speech_enc+L_FRAME32k, L_FRAME );

    /*------------------------------------------------------------------*
     * Perform fixed preemphasis (12.8 kHz signal) through 1 - g*z^-1
     *-----------------------------------------------------------------*/

    preemph( new_inp_12k8, PREEMPH_FAC, L_FRAME, &st->mem_preemph );

    /*-------------------------------------------------------------------------*
     * Spectral analysis
     *--------------------------------------------------------------------------*/

    analy_sp( inp_12k8, st->Bin_E, st->Bin_E_old, fr_bands, lf_E, Etot, st->min_band, st->max_band, band_energies, PS, fft_buff );

    /*----------------------------------------------------------------*
     * SAD (1-signal, 0-noise)
     *----------------------------------------------------------------*/

    noise_est_pre( *Etot, st->ini_frame, &st->Etot_l, &st->Etot_h, &st->Etot_l_lp, &st->Etot_last,
                   &st->Etot_v_h2, &st->sign_dyn_lp, st->harm_cor_cnt, &st->Etot_lp );

    *vad_flag = wb_vad( st, fr_bands, localVAD, &noisy_speech_HO, &clean_speech_HO, &NB_speech_HO, &snr_sum_he, &localVAD_HE_SAD, &(st->flag_noisy_speech_snr) );

    vad_flag_cldfb = vad_proc( realBuffer, imagBuffer, enerBuffer, st->cldfbAnaEnc->no_channels, &st->vad_st, &cldfb_addition,*vad_flag );


    /* Combine decisions from SADS */
    if ( *vad_flag == 1 && vad_flag_cldfb == 0 )
    {
        *localVAD = 0;
    }

    *vad_flag = vad_flag_cldfb;

    vad_flag_dtx = dtx_hangover_addition( st, *localVAD, *vad_flag, st->lp_speech-st->lp_noise, cldfb_addition, vad_hover_flag );

    /*----------------------------------------------------------------*
     * NB/WB/SWB/FB bandwidth detector
     *----------------------------------------------------------------*/

    bw_detect( st, signal_in, *localVAD, enerBuffer );

    /*----------------------------------------------------------------*
     * Noise energy down-ward update and total noise energy estimation
     * Long-term energies and relative frame energy updates
     * Correlation correction as a function of total noise level
     *----------------------------------------------------------------*/

    noise_est_down( fr_bands, st->bckr, tmpN, tmpE, st->min_band, st->max_band, &st->totalNoise, *Etot, &st->Etot_last, &st->Etot_v_h2 );

    relE = *Etot - st->lp_speech;

    corr_shift = correlation_shift( st->totalNoise );

    /*----------------------------------------------------------------*
     * FD-CNG Noise Estimator
     *----------------------------------------------------------------*/

    resetFdCngEnc ( st );

    perform_noise_estimation_enc( band_energies, enerBuffer, st->hFdCngEnc );

    /*-----------------------------------------------------------------*
     * Select SID or FRAME_NO_DATA frame if DTX enabled
     *-----------------------------------------------------------------*/

    dtx( st, vad_flag_dtx, inp_12k8 );

    /*----------------------------------------------------------------*
     * Adjust FD-CNG Noise Estimator
     *----------------------------------------------------------------*/

    if( (st->last_total_brate != st->total_brate) || (st->last_bwidth != st->bwidth) )
    {
        configureFdCngEnc( st->hFdCngEnc, st->bwidth, st->rf_mode&&st->total_brate==13200?9600:st->total_brate );
    }
    if ( st->hFdCngEnc != NULL && st->Opt_DTX_ON )
    {
        AdjustFirstSID( st->hFdCngEnc->hFdCngCom->npart, st->hFdCngEnc->msPeriodog, st->hFdCngEnc->energy_ho,
                        st->hFdCngEnc->msNoiseEst, st->hFdCngEnc->msNoiseEst_old, &(st->hFdCngEnc->hFdCngCom->active_frame_counter), st );
    }

    /*----------------------------------------------------------------*
     * Reconfigure Mode 2
     *----------------------------------------------------------------*/

    if ( st->codec_mode == MODE2 )
    {
        SetModeIndex( st, st->total_brate, st->bwidth );
    }

    calcLoEnvCheckCorrHiLo( st->cldfbAnaEnc->no_col, freqTable, st->tecEnc.loBuffer, st->tecEnc.loTempEnv,
                            st->tecEnc.loTempEnv_ns, st->tecEnc.hiTempEnv, &(st->tecEnc.corrFlag) );

    /*---------------------------------------------------------------*
     * Time Domain Transient Detector
     *---------------------------------------------------------------*/

    if( st->tcx10Enabled || st->tcx20Enabled )
    {
        RunTransientDetection( signal_in, input_frame, &st->transientDetection );

        currFlatness = GetTCXAvgTemporalFlatnessMeasure( &st->transientDetection, NSUBBLOCKS, 0 );
    }

    /*----------------------------------------------------------------*
     * LP analysis
     *----------------------------------------------------------------*/

    alw_pitch_lag_12k8[0] = st->old_pitch_la;
    alw_pitch_lag_12k8[1] = st->old_pitch_la;
    alw_voicing[0] = st->old_voicing_la;
    alw_voicing[1] = st->old_voicing_la;

    analy_lp( inp_12k8, L_FRAME, L_look, ener, A, epsP, lsp_new, lsp_mid, st->lsp_old1, alw_pitch_lag_12k8, alw_voicing, 12800 );

    lsp2lsf( lsp_new, lsf_new, M, 12800 );
    stab_fac = lsf_stab( lsf_new, st->lsf_old1, 0, L_FRAME );
    mvr2r( lsf_new, st->lsf_old1, M );

    /*----------------------------------------------------------------*
     * Compute weighted input (for OL pitch analysis)
     * OL pitch analysis
     * stable high pitch detection
     * 1/4 pitch precision improvement
     *----------------------------------------------------------------*/

    find_wsp( L_FRAME, L_SUBFR, NB_SUBFR, A, Aw, inp_12k8, TILT_FAC, wsp, &st->mem_wsp, GAMMA1, L_look );

    if( *vad_flag == 0 )
    {
        /* reset the OL pitch tracker memories during inactive frames */
        pitch_ol_init( &st->old_thres, &st->old_pitch, &st->delta_pit, &st->old_corr );
    }

    pitch_ol( pitch, voicing, &st->old_pitch, &st->old_corr, corr_shift, &st->old_thres, &st->delta_pit,
              st->old_wsp2, wsp, st->mem_decim2, relE, L_look, st->clas, st->input_bwidth, st->Opt_SC_VBR );

    /* Updates for adaptive lag window memory */
    st->old_pitch_la = pitch[2];
    st->old_voicing_la = voicing[2];

    /* Detection of very short stable pitch period (MODE1 bit-rates) */
    StableHighPitchDetect( &flag_spitch, pitch, voicing, st->Bin_E, wsp, *localVAD, &st->voicing_sm, &st->voicing0_sm,
                           &st->LF_EnergyRatio_sm, &st->predecision_flag, &st->diff_sm, &st->energy_sm );

    /* 1/4 pitch precision improvement */
    if( st->total_brate <= ACELP_24k40 )
    {
        pitch_ol2( PIT_MIN_EXTEND, pitch[0], &pitch_fr[0], &voicing_fr[0], 0, wsp, 7 );
        pitch_ol2( PIT_MIN_EXTEND, pitch[0], &pitch_fr[1], &voicing_fr[1], L_SUBFR, wsp, 7 );
        pitch_ol2( PIT_MIN_EXTEND, pitch[1], &pitch_fr[2], &voicing_fr[2], 2*L_SUBFR, wsp, 7 );
        pitch_ol2( PIT_MIN_EXTEND, pitch[1], &pitch_fr[3], &voicing_fr[3], 3*L_SUBFR, wsp, 7 );
    }
    else
    {
        pitch_fr[0] = pitch[0];
        pitch_fr[1] = pitch[0];
        pitch_fr[2] = pitch[1];
        pitch_fr[3] = pitch[1];

        voicing_fr[0] = voicing[0];
        voicing_fr[1] = voicing[0];
        voicing_fr[2] = voicing[1];
        voicing_fr[3] = voicing[1];
    }

    /*------------------------------------------------------------------*
     * Update estimated noise energy and voicing cut-off frequency
     *-----------------------------------------------------------------*/

    noise_est( st, tmpN, pitch, voicing, epsP, *Etot, relE, corr_shift, tmpE, fr_bands, &cor_map_sum,
               &sp_div, &non_staX, &loc_harm, lf_E, &st->harm_cor_cnt, st->Etot_l_lp, &sp_floor );

    /*------------------------------------------------------------------*
     * Update parameters used in the VAD and DTX
     *-----------------------------------------------------------------*/

    vad_param_updt( st, pitch, voicing, corr_shift, A );

    /*-----------------------------------------------------------------*
     * Find spectral tilt
     * UC and VC frame selection
     *-----------------------------------------------------------------*/

    find_tilt( fr_bands, st->bckr, ee, pitch, voicing, lf_E, corr_shift, st->input_bwidth,
               st->max_band, hp_E, st->codec_mode, &(st->bckr_tilt_lt), st->Opt_SC_VBR);

    *coder_type = find_uv( st, pitch_fr, voicing_fr, voicing, inp_12k8, *localVAD, ee,
                           corr_shift, relE, *Etot, hp_E, &flag_spitch, st->voicing_sm, last_core_orig );

    /*-----------------------------------------------------------------*
     * channel aware mode configuration                                *
     *-----------------------------------------------------------------*/

    if( !st->Opt_RF_ON )
    {
        st->rf_mode = 0;
        st->rf_target_bits_write = 0;
    }
    else if( st->rf_mode && st->core_brate != FRAME_NO_DATA && st->core_brate != SID_2k40 )
    {
        /* the RF config is for (n- fec_offset)th frame that will be packed along with the n-th frame bistream */
        st->rf_mode = 1;
        st->codec_mode = MODE2;

        st->rf_target_bits_write = st->rf_targetbits_buff[st->rf_fec_offset];
    }
    else
    {
        st->rf_mode = 0;
        st->codec_mode = MODE1;
        st->rf_indx_frametype[0] = RF_NO_DATA;
        st->rf_targetbits_buff[0] = 6;  /* rf_mode: 1, rf_frame_type: 3, and fec_offset: 2 */
    }

    /*-----------------------------------------------------------------*
     * Signal classification for FEC
     * TC frame selection
     *-----------------------------------------------------------------*/

    st->clas = signal_clas( st, coder_type, voicing, inp_12k8, *localVAD, pitch, ee, relE, L_look, &uc_clas );

    st->Local_VAD = *localVAD;

    /*----------------------------------------------------------------*
     * Speech/music classification
     * AC frame selection
     *----------------------------------------------------------------*/

    speech_music_classif( st, &sp_aud_decision0, sp_aud_decision1, sp_aud_decision2, new_inp_12k8, inp_12k8, *vad_flag,
                          *localVAD, localVAD_HE_SAD, pitch, voicing, lsp_new, cor_map_sum, epsP, PS,
                          *Etot, old_cor, coder_type, attack_flag, non_staX, relE, &high_lpn_flag, flag_spitch );

    long_enr( st, *Etot, localVAD_HE_SAD, high_lpn_flag );

    /*----------------------------------------------------------------*
     * Final VAD correction ( when HE-SAD is used instead of the normal VAD,
     * rewrite the VAD flag by VAD flag  with DTX hangover for further processing)
     *----------------------------------------------------------------*/

    if( st->Opt_DTX_ON )
    {
        *vad_flag = vad_flag_dtx;
    }

    /*----------------------------------------------------------------*
     * Selection of internal ACELP Fs (12.8 kHz or 16 kHz)
     *----------------------------------------------------------------*/

    if( st->codec_mode == MODE1 )
    {
        if( st->core_brate == FRAME_NO_DATA )
        {
            /* prevent "L_frame" changes in CNG segments */
            st->L_frame = st->last_L_frame;
        }
        else if( st->core_brate == SID_2k40 && st->bwidth >= WB && st->first_CNG && st->act_cnt2 < MIN_ACT_CNG_UPD )
        {
            /* prevent "L_frame" changes in SID frame after short segment of active frames */
            st->L_frame = st->last_CNG_L_frame;
        }
        else if( ( st->core_brate == SID_2k40 && st->total_brate >= ACELP_9k60 && ((st->bwidth == WB && !( st->total_brate == ACELP_13k20 && st->cng_type == FD_CNG )) || (st->cng_type == LP_CNG && st->bwidth > WB && st->total_brate >= ACELP_16k40)) ) ||
                 ( st->total_brate > ACELP_24k40 && st->total_brate < HQ_96k ) || ( st->total_brate == ACELP_24k40 && st->bwidth >= WB ) )
        {
            st->L_frame = L_FRAME16k;
        }
        else
        {
            st->L_frame = L_FRAME;
        }

        if( st->ini_frame == 0 )
        {
            /* avoid switching of internal ACELP Fs in the very first frame */
            st->last_L_frame = st->L_frame;
        }

        if( st->L_frame == L_FRAME )
        {
            st->gamma = GAMMA1;
            st->preemph_fac = PREEMPH_FAC;
        }
        else
        {
            st->gamma = GAMMA16k;
            st->preemph_fac = PREEMPH_FAC_16k;
        }

        st->sr_core = 50*st->L_frame;
        st->encoderLookahead_enc = NS2SA(st->sr_core, ACELP_LOOK_NS);
        st->encoderPastSamples_enc = (st->L_frame*9)>>4;
    }

    /*-----------------------------------------------------------------*
     * coder_type rewriting in case of switching
     * IC frames selection
     * enforce TC frames in case of switching
     *-----------------------------------------------------------------*/

    if( st->codec_mode == MODE1 )
    {
        /* enforce TRANSITION frames */
        if( st->last_L_frame != st->L_frame && st->core_brate != FRAME_NO_DATA && st->core_brate != SID_2k40 && (st->coder_type_raw != VOICED) )
        {
            /* enforce TC frame in case of ACELP@12k8 <-> ACELP@16k core switching */
            *coder_type = TRANSITION;
        }
        else if( st->last_core == HQ_CORE || st->last_core == TCX_10_CORE ||  st->last_core == TCX_20_CORE )
        {
            /* enforce TC frame in case of HQ -> ACELP core switching */
            *coder_type = TRANSITION;
        }
        else if( st->last_core_brate <= SID_2k40 && st->cng_type == FD_CNG )
        {
            /* enforce TC frame in case of FD_CNG -> ACELP switching (past excitation not available) */
            *coder_type = TRANSITION;
        }

        /* select INACTIVE frames */
        else if( st->total_brate <= ACELP_24k40 && *vad_flag == 0 )
        {
            /* inactive frames will be coded by GSC technology */
            /* except for the VBR mode. VBR mode uses NELP for that */
            if( !( st->Opt_SC_VBR && vad_flag_dtx) )
            {
                *coder_type = INACTIVE;
                st->noise_lev = NOISE_LEVEL_SP3;
            }
        }
        else if( st->total_brate > ACELP_24k40 &&
                 ( (*vad_flag == 0 && st->bwidth >= SWB && st->max_bwidth >= SWB) || (*localVAD == 0 && (st->bwidth <= WB || st->max_bwidth <= WB)) )
               )
        {
            /* inactive frames will be coded by AVQ technology */
            *coder_type = INACTIVE;
        }
    }
    else /* st->codec_mode == MODE2 */
    {
        if( !(*vad_flag) )
        {
            *coder_type = INACTIVE;
        }
        else if( *coder_type > GENERIC )
        {
            *coder_type = GENERIC;
        }
    }

    /*---------------------------------------------------------------*
     * SC-VBR - decision about PPP/NELP mode
     *---------------------------------------------------------------*/

    if( st->Opt_SC_VBR )
    {
        set_ppp_mode( st, coder_type, noisy_speech_HO, clean_speech_HO, NB_speech_HO,
                      *localVAD, localVAD_HE_SAD, vad_flag, pitch, *sp_aud_decision1 );
    }

    if ( !st->Opt_AMR_WB && !st->rf_mode )
    {
        if ( st->total_brate == ACELP_13k20 || st->total_brate == ACELP_32k )
        {
            st->mdct_sw_enable = MODE1;
        }
        else if ( ACELP_16k40 <= st->total_brate && st->total_brate <= ACELP_24k40 )
        {
            st->mdct_sw_enable = MODE2;
        }
    }

    if( st->codec_mode == MODE1 )
    {
        /*---------------------------------------------------------------------*
         * Decision matrix (selection of technologies)
         *---------------------------------------------------------------------*/

        decision_matrix_enc( st, *sp_aud_decision1, *sp_aud_decision2, *coder_type, *vad_flag, hq_core_type );

        /* HQ_CORE/TCX_20_CORE decision */
        if( st->core == HQ_CORE ) /* Decision matrix decided for MDCT coding */
        {
            if( (st->bwidth == SWB || st->bwidth == FB) && st->total_brate == ACELP_32k )
            {
                /* Select MDCT Core */
                st->core = mdct_classifier( fft_buff, st, *vad_flag, enerBuffer );
            }
            if( st->total_brate == ACELP_13k20 && st->bwidth != FB )
            {
                MDCT_selector( st, sp_floor, *Etot, cor_map_sum, voicing, enerBuffer, *vad_flag );
            }
        }
        else
        {
            MDCT_selector_reset( st );
        }

        /* Switch to MODE2 if TCX_20_CORE */
        if( st->core == TCX_20_CORE )
        {
            st->codec_mode = MODE2;

            if ( st->last_codec_mode == MODE1 )
            {
                int last_total_brate = st->last_total_brate;
                st->last_total_brate = -1;
                SetModeIndex( st, st->total_brate, st->bwidth );
                st->last_total_brate = last_total_brate;
            }
            else
            {
                SetModeIndex( st, st->total_brate, st->bwidth );

                st->sr_core = getCoreSamplerateMode2( st->total_brate, st->bwidth, st->rf_mode );
                st->L_frame = st->sr_core / 50;
                st->encoderLookahead_enc = NS2SA(st->sr_core, ACELP_LOOK_NS);
                st->encoderPastSamples_enc = (st->L_frame*9)>>4;

                if ( st->sr_core == 12800 )
                {
                    st->preemph_fac = PREEMPH_FAC;
                    st->gamma = GAMMA1;
                }
                else
                {
                    st->preemph_fac = PREEMPH_FAC_16k;
                    st->gamma = GAMMA16k;
                }

                st->igf = getIgfPresent( st->total_brate, st->bwidth, st->rf_mode );
            }

            *coder_type = st->coder_type_raw;

            if( *vad_flag == 0 )
            {
                *coder_type = INACTIVE;
            }
            else if( *coder_type > GENERIC )
            {
                *coder_type = GENERIC;
            }

            st->mdct_sw = MODE1;
        }
    }

    /*-----------------------------------------------------------------*
     * Update of ACELP harmonicity counter (used in ACELP transform codebook @32kbps)
     *-----------------------------------------------------------------*/

    if( st->total_brate == ACELP_32k && loc_harm == 1 && cor_map_sum > 50 &&
            st->clas == VOICED_CLAS && *coder_type == GENERIC )
    {
        st->last_harm_flag_acelp++;

        if( st->last_harm_flag_acelp > 10 )
        {
            st->last_harm_flag_acelp = 10;
        }
    }
    else
    {
        st->last_harm_flag_acelp = 0;
    }

    /*-----------------------------------------------------------------*
     * Update audio frames counter (used for UV decision)
     *-----------------------------------------------------------------*/

    if( *coder_type == AUDIO )
    {
        st->audio_frame_cnt += AUDIO_COUNTER_STEP;
    }
    else if (*coder_type != INACTIVE)
    {
        st->audio_frame_cnt--;
    }

    if( st->audio_frame_cnt > AUDIO_COUNTER_MAX )
    {
        st->audio_frame_cnt = AUDIO_COUNTER_MAX;
    }

    if( st->audio_frame_cnt < 0 )
    {
        st->audio_frame_cnt = 0;
    }

    /*-----------------------------------------------------------------*
     * Set formant sharpening flag
     *-----------------------------------------------------------------*/

    *sharpFlag = 0;
    if( *coder_type == TRANSITION )
    {
        if( ( st->total_brate > ACELP_48k && st->bwidth < SWB ) ||                                     /* Deactivate for core bitrates higher than 48.0 kb/s */
                ( st->total_brate >= ACELP_13k20 && st->total_brate <= ACELP_16k40 ) ||                    /* Deactivate for bitrates <13.2, 16.4> kb/s (this is basically due to lack of signaling configurations */
                ( st->total_brate > ACELP_16k40 && st->lp_noise > FORMANT_SHARPENING_NOISE_THRESHOLD )  )  /* Deactivate for bitrates >= 24.4 kb/s if the long-term noise level exceeds 34 dB */
        {
            *sharpFlag = 0;
        }
        else
        {
            *sharpFlag = 1;
        }
    }

    if( *coder_type == GENERIC || *coder_type == VOICED )
    {
        if( *vad_hover_flag ||
                ( st->total_brate > ACELP_48k && st->bwidth < SWB ) ||               /* Deactivate for core bitrates higher than 48.0 kb/s */
                ( st->total_brate >= ACELP_13k20 && st->lp_noise > FORMANT_SHARPENING_NOISE_THRESHOLD
                  && st->total_brate > CNA_MAX_BRATE)  ) /* Deactivate for bitrates >= 13.2 kb/s if the long-term noise level exceeds 34 dB */
        {
            *sharpFlag = 0;
        }
        else
        {
            *sharpFlag = 1;
        }
    }

    /* channel-aware mode - due to lack of signalling bit, sharpFlag is 1 always in RF mode */
    if( st->rf_mode && ( *coder_type == VOICED || *coder_type == GENERIC ) )
    {
        *sharpFlag = 1;
    }

    /*-----------------------------------------------------------------*
     * Set voicing flag for HQ FEC
     *-----------------------------------------------------------------*/

    if ( *sp_aud_decision1 == 0 && ( *coder_type == VOICED || *coder_type == GENERIC ) )
    {
        *Voicing_flag = 1;
    }
    else
    {
        *Voicing_flag = 0;
    }

    /*---------------------------------------------------------------*
     * Preprocessing at other sampling frequency rate (16/25.6/32kHz)
     *----------------------------------------------------------------*/

    sr_core_tmp = (st->codec_mode == MODE1)? 16000 : max(16000,st->sr_core);      /* indicates the ACELP sampling rate */
    L_frame_tmp = (st->codec_mode == MODE1)? L_FRAME16k : max(L_FRAME16k,st->L_frame);

    L_look = NS2SA(sr_core_tmp, ACELP_LOOK_NS);             /* lookahead at other sampling rate (16kHz, 25.6kHz, 32kHz) */

    new_inp_16k = old_inp_16k + L_INP_MEM;                  /* pointer to new samples of the input signal in 16kHz core */
    inp_16k = new_inp_16k - L_look;                         /* pointer to the current frame of input signal in 16kHz core */
    mvr2r( st->old_inp_16k, old_inp_16k, L_INP_MEM );

    /*---------------------------------------------------------------*
     * Change the sampling frequency to 16/25.6/32 kHz
     *----------------------------------------------------------------*/

    if( st->input_Fs == sr_core_tmp )
    {
        /* no resampling needed, only delay adjustment to account for the FIR resampling delay */
        delay = NS2SA(st->input_Fs, DELAY_FIR_RESAMPL_NS);
        mvr2r( st->mem_decim16k + delay, new_inp_16k, delay );
        mvr2r( signal_in, new_inp_16k + delay, input_frame - delay );
        mvr2r( signal_in + input_frame - 2*delay, st->mem_decim16k, 2*delay );
    }
    else if( st->input_Fs == 32000 || st->input_Fs == 48000 )
    {
        modify_Fs( signal_in, input_frame, st->input_Fs, new_inp_16k, sr_core_tmp, st->mem_decim16k, 0 );
    }
    else    /* keep memories up-to-date in case of bit-rate switching */
    {
        /* no resampling needed, only delay adjustment to account for the FIR resampling delay */
        delay = NS2SA(st->input_Fs, DELAY_FIR_RESAMPL_NS);
        mvr2r( st->mem_decim16k + delay, new_inp_16k, delay );
        mvr2r( signal_in, new_inp_16k + delay, input_frame - delay );
        mvr2r( signal_in + input_frame - 2*delay, st->mem_decim16k, 2*delay );
    }

    if( sr_core_tmp == 16000 )
    {
        /* save input resampled at 16kHz, non-preemhasised */
        mvr2r( new_inp_16k, new_inp_resamp16k, L_FRAME16k );
    }
    else if( sr_core_tmp > 16000 )
    {
        /* reset the buffer, the signal is needed for WB BWEs */
        set_f( new_inp_resamp16k, 0.0f, L_FRAME16k );
    }

    /*------------------------------------------------------------------*
     * Perform fixed preemphasis (16kHz signal) through 1 - g*z^-1
     *-----------------------------------------------------------------*/

    if( (st->tcxonly == 0) || ( st->codec_mode == MODE1 && st->input_Fs > 8000 ) )
    {
        st->mem_preemph_enc = new_inp_16k[L_frame_tmp-1];
    }

    if( st->input_Fs > 8000 && sr_core_tmp == 16000 )
    {
        preemph( new_inp_16k, PREEMPH_FAC_16k, L_FRAME16k, &(st->mem_preemph16k) );
    }
    else if( st->input_Fs > 8000 )    /* keep memory up-to-date in case of bit-rate switching */
    {
        st->mem_preemph16k = new_inp_16k[L_frame_tmp-1];
    }

    /*-----------------------------------------------------------------*
     * Redo LP analysis at 16kHz if ACELP@16k core was selected
     * update buffers
     *-----------------------------------------------------------------*/

    if( ( ((st->tcxonly == 0) || !(st->core_brate != FRAME_NO_DATA || st->core_brate != SID_2k40)) && st->L_frame == L_FRAME16k && st->codec_mode == MODE2 ) ||
            ( st->L_frame == L_FRAME16k && st->codec_mode == MODE1 ) )
    {
        /* update signal buffers */
        mvr2r( new_inp_resamp16k, st->buf_speech_enc+L_FRAME16k, L_FRAME16k );
        mvr2r( new_inp_16k, st->buf_speech_enc_pe+L_FRAME16k, L_FRAME16k );

        /*--------------------------------------------------------------*
         * LPC analysis
         *---------------------------------------------------------------*/

        if( st->last_L_frame == L_FRAME && st->codec_mode == MODE1 )
        {
            /* this is just an approximation, but it is sufficient */
            mvr2r( st->lsp_old1, st->lspold_enc, M );
        }

        analy_lp( inp_16k, L_FRAME16k, L_look, ener, A, epsP, lsp_new, lsp_mid, st->lspold_enc, pitch, voicing, 16000 );

        /*--------------------------------------------------------------*
        * Compute Weighted Input
        *---------------------------------------------------------------*/

        if( st->codec_mode == MODE2 )
        {
            find_wsp( L_FRAME16k, L_SUBFR, st->nb_subfr, A, Aw, st->speech_enc_pe, PREEMPH_FAC_16k, st->wspeech_enc, &st->mem_wsp_enc, st->gamma, L_LOOK_16k );
        }
        else
        {
            weight_a_subfr( NB_SUBFR16k, A, Aw, GAMMA16k, M );
        }

    }
    else
    {
        /* update signal buffers */
        mvr2r( new_inp_12k8, st->buf_speech_enc_pe+st->L_frame, L_FRAME );
        mvr2r( st->buf_speech_enc+L_FRAME32k, st->buf_speech_enc+st->L_frame, L_FRAME );

        if( st->tcxonly == 0 )
        {
            mvr2r( wsp, st->wspeech_enc, L_FRAME + L_LOOK_12k8 );
        }
    }

    /*-----------------------------------------------------------------*
     * ACELP/TCX20 Switching Decision
     *-----------------------------------------------------------------*/

    if( st->codec_mode == MODE2 )
    {
        if( st->core_brate != FRAME_NO_DATA && st->core_brate != SID_2k40 && st->tcxonly == 0 )
        {
            core_acelp_tcx20_switching( st, *vad_flag, sp_aud_decision0, non_staX, pitch, pitch_fr,
                                        voicing_fr, currFlatness, lsp_mid, stab_fac );
        }

        if( st->mdct_sw_enable == MODE2 && !st->rf_mode )
        {
            if (st->core == TCX_20_CORE) /* Switching only possible from TCX_20 frames, not from TCX_10 frames */
            {
                /* Select MDCT Core */
                if( (st->bwidth == SWB || st->bwidth == FB) && st->total_brate == ACELP_24k40 )
                {
                    st->core = mdct_classifier( fft_buff, st, *vad_flag, enerBuffer );
                }

                if( st->total_brate == ACELP_16k40 && st->bwidth != FB )
                {
                    MDCT_selector( st, sp_floor, *Etot, cor_map_sum, voicing, enerBuffer, *vad_flag );
                }
            }
            else
            {
                MDCT_selector_reset( st );
            }

            /* Do the switching that was decided in the MDCT selector */
            if( st->core == HQ_CORE )
            {
                st->codec_mode = MODE1;
                st->mdct_sw = MODE2;
            }
            else if( st->last_codec_mode == MODE1 && st->last_core == HQ_CORE )
            {
                int L_frame_old = st->last_L_frame;
                st->last_L_frame = st->L_frame;
                SetModeIndex( st, st->total_brate, st->bwidth );
                st->last_L_frame = L_frame_old;
            }
        }

        /*--------------------------------------------------------------*
         * TCX mode decision
         *---------------------------------------------------------------*/

        SetTCXModeInfo( st, &st->transientDetection, &st->tcx_cfg.tcx_curr_overlap_mode );
    }

    /*-----------------------------------------------------------------*
     * Updates
     *-----------------------------------------------------------------*/

    /* update old weighted speech buffer - for OL pitch analysis */
    mvr2r( &old_wsp[L_FRAME], st->old_wsp, L_WSP_MEM );

    /* update old input signal buffer */
    mvr2r( &old_inp_12k8[L_FRAME], st->old_inp_12k8, L_INP_MEM );

    /* update old input signal @16kHz buffer */
    if( st->input_Fs > 8000 && sr_core_tmp == 16000 )
    {
        mvr2r( &old_inp_16k[L_frame_tmp], st->old_inp_16k, L_INP_MEM );
    }
    else if( st->input_Fs > 8000 )
    {
        lerp( st->old_inp_12k8+L_INP_MEM-L_INP_MEM*4/5, st->old_inp_16k, L_INP_MEM, L_INP_MEM*4/5 );
    }

    if( (sr_core_tmp == 16000) && st->tcxonly && st->codec_mode == MODE2 )
    {
        /* copy input resampled at 16kHz, non-preemhasised */
        mvr2r( new_inp_resamp16k, new_inp_16k, L_FRAME16k );
    }

    /* update of old per-band energy spectrum */
    mvr2r( fr_bands + NB_BANDS, st->enrO, NB_BANDS );

    /* set the pointer of the current frame for the ACELP core */
    if ( st->L_frame == L_FRAME )
    {
        *inp = inp_12k8;
    }
    else
    {
        *inp = inp_16k;
    }


    return;
}
