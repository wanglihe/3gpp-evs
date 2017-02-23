/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_enc.h"
#include "rom_com.h"
#include "prot.h"

/*-------------------------------------------------------------------*
 * amr_wb_enc()
 *
 * AMR-WB encoder
 *-------------------------------------------------------------------*/

void amr_wb_enc(
    Encoder_State *st,                      /* i/o: encoder state structure */
    const short *input_sp,          /* i  : input signal            */
    const short n_samples           /* i  : number of input samples */
)
{
    short i, delay, harm_flag;
    float old_inp[L_INP_12k8], *new_inp, *inp;                /* buffer of old input signal           */
    float old_inp_16k[L_INP_12k8+L_SUBFR], *inp_16k, *new_inp_16k;/* buffer of old input signal @16kHz*/
    float old_exc[L_EXC], *exc;                               /* excitation signal buffer             */
    float old_wsp[L_WSP], *wsp;                               /* weighted input signal buffer         */
    short input_frame;                                        /* frame length at input sampling freq. */
    float fr_bands[2*NB_BANDS];                               /* energy in frequency bands            */
    float lf_E[2*VOIC_BINS];                                  /* per bin spectrum energy in lf        */
    float tmpN[NB_BANDS];                                     /* temporary noise update               */
    float tmpE[NB_BANDS], PS[L_FFT/2];                        /* temporary averaged energy of 2 sf.   */
    float corr_shift;                                         /* correlation shift                    */
    float relE;                                               /* frame relative energy                */
    float non_staX, cor_map_sum, sp_div;
    short vad_flag;
    short localVAD;
    float Etot;                                               /* total energy                         */
    float ener;                                               /* residual energy from Levinson-Durbin */
    short pitch[3];                                           /* open-loop pitch values               */
    float voicing[3];                                         /* open-loop pitch gains                */
    float A[NB_SUBFR*(M+1)];                                  /* A(z) unquantized for the 4 subframes */
    float Aw[NB_SUBFR*(M+1)];                                 /* weigted A(z) unquant. for 4 subframes*/
    float epsP[M+1];                                          /* LP prediction errors                 */
    float isp_new[M];                                         /* ISPs at the end of the frame         */
    float isf_new[M];                                         /* ISFs at the end of the frame         */
    float isp_tmp[M];
    float Aq[NB_SUBFR*(M+1)];                                 /* A(z) quantized for the 4 subframes   */
    float syn[L_FRAME];                                       /* synthesis vector                     */
    float res[L_FRAME];                                       /* residual signal for FER protection   */
    float exc2[L_FRAME];                                      /* enhanced excitation                  */
    float pitch_buf[NB_SUBFR];                                /* floating pitch for each subframe     */
    float dummy_buf[L_FRAME32k];                              /* dummy buffer - no usage              */
    float snr_sum_he;
    short allow_cn_step;
    short localVAD_HE_SAD;
    short tmps;
    short vad_flag_dtx;
    short vad_hover_flag;
    short coder_type;
    short hf_gain[NB_SUBFR];
    short high_lpn_flag;
    float lp_bckr, hp_bckr;
    float q_env[NUM_ENV_CNG];
    short sid_bw = 0;
    float exc3[L_FRAME];
    float fft_buff[2*L_FFT];
    float sp_floor;
    float tmp;


    /*------------------------------------------------------------------*
     * Initialization
     *------------------------------------------------------------------*/

    st->L_frame = L_FRAME;
    st->gamma = GAMMA1;
    st->core = AMR_WB_CORE;
    st->core_brate = st->total_brate;
    st->input_bwidth = st->last_input_bwidth;
    st->bwidth = st->last_bwidth;
    st->extl = -1;
    coder_type = GENERIC;
    input_frame = (short)(st->input_Fs / 50);                   /* frame length of the input signal */
    st->encoderPastSamples_enc = (L_FRAME*9)/16;
    st->encoderLookahead_enc = L_LOOK_12k8;

    st->bpf_off = 0;
    if( st->last_core == HQ_CORE || st->last_codec_mode == MODE2 )
    {
        st->bpf_off = 1;
    }

    st->igf = 0;

    /* Updates in case of EVS primary mode -> AMR-WB IO mode switching */
    if( st->last_core != AMR_WB_CORE )
    {
        updt_IO_switch_enc( st, input_frame );
    }

    /* Updates in case of HQ -> AMR-WB IO switching */
    core_switching_pre_enc( st, &(st->LPDmem), NULL, NULL );

    set_s( hf_gain, 0, NB_SUBFR );

    set_f( old_inp, 0.0f, L_INP_12k8 );
    exc = old_exc + L_EXC_MEM;                                  /* pointer to excitation signal in the current frame */
    mvr2r( st->LPDmem.old_exc, old_exc, L_EXC_MEM );


    new_inp = old_inp + L_INP_MEM;                              /* pointer to new samples of the input signal */
    inp = new_inp - L_LOOK_12k8;                                /* pointer to current frame of input signal */
    wsp = old_wsp + L_WSP_MEM;                                  /* pointer to current frame of weighted signal */

    mvr2r( st->old_inp_12k8, old_inp, L_INP_MEM );
    mvr2r( st->old_wsp, old_wsp, L_WSP_MEM );

    new_inp_16k = old_inp_16k + L_INP_MEM;                      /* pointer to new samples of the input signal in 16kHz core */
    inp_16k = new_inp_16k - L_LOOK_16k;                         /* pointer to the current frame of input signal in 16kHz core */
    mvr2r( st->old_inp_16k, old_inp_16k, L_INP_MEM );

    /* in case of switching, reset AMR-WB BWE memories */
    if( st->total_brate == ACELP_23k85 && st->last_core_brate != ACELP_23k85 )
    {
        hf_cod_init( st->mem_hp400_enc, st->mem_hf_enc, st->mem_syn_hf_enc, st->mem_hf2_enc, &st->gain_alpha );
    }

    /*----------------------------------------------------------------*
     * set input samples buffer
     *----------------------------------------------------------------*/

    /* get delay to synchronize ACELP and MDCT frame */
    delay = NS2SA(st->input_Fs, DELAY_FIR_RESAMPL_NS);

    mvr2r( st->input - delay, st->old_input_signal, input_frame+delay );

    /*----------------------------------------------------------------*
     * Buffering of input signal
     * (convert 'short' input data to 'float')
     * HP filtering
     *----------------------------------------------------------------*/

    mvs2r( input_sp, st->input, n_samples );

    for( i = n_samples; i < input_frame; i++ )
    {
        st->input[i] = 0;
    }

    hp20( st->input, input_frame, st->mem_hp20_in, st->input_Fs );

    /*-----------------------------------------------------------------*
     * switching from ACELP@16k core to AMR-WB IO mode
     *-----------------------------------------------------------------*/

    st->rate_switching_reset = 0;

    if( st->last_core != AMR_WB_CORE && st->last_L_frame == L_FRAME16k && st->last_core != HQ_CORE)
    {
        /* in case of switching, do not apply BPF */
        st->bpf_off = 1;
        /* convert old quantized LSP vector */
        st->rate_switching_reset = lsp_convert_poly( st->lsp_old, L_FRAME, 1 );

        /* convert old quantized LSF vector */
        lsp2lsf( st->lsp_old, st->lsf_old, M, INT_FS_12k8 );

        /* Reset LPC mem */
        mvr2r( GEWB_Ave, st->mem_AR, M );
        set_zero( st->mem_MA, M );

        /* update synthesis filter memories */
        synth_mem_updt2( L_FRAME, st->last_L_frame, st->LPDmem.old_exc, st->LPDmem.mem_syn_r, st->mem_syn1, st->LPDmem.mem_syn, ENC );
        mvr2r( st->LPDmem.old_exc, old_exc, L_EXC_MEM );
        mvr2r( st->mem_syn1, st->LPDmem.mem_syn2, M );
        mvr2r( st->LPDmem.mem_syn2, st->LPDmem.mem_syn3, M );

        /* lsp -> isp */
        mvr2r( stable_ISP, isp_tmp, M );
        lsp2isp( st->lsp_old, st->lsp_old, isp_tmp, M );

    }

    /* update buffer of old subframe pitch values */
    if( st->last_L_frame != L_FRAME )
    {
        if( st->last_L_frame == L_FRAME32k )
        {
            tmp = (float)12800/(float)32000;
        }
        else if( st->last_L_frame == 512 )
        {
            tmp = (float)12800/(float)25600;
        }
        else /* st->last_L_frame == L_FRAME16k */
        {
            tmp = (float)12800/(float)16000;
        }

        for( i=NB_SUBFR16k-NB_SUBFR; i<NB_SUBFR16k; i++ )
        {
            st->old_pitch_buf[i-1] = tmp * st->old_pitch_buf[i];
        }

        for( i=2*NB_SUBFR16k-NB_SUBFR; i<2*NB_SUBFR16k; i++ )
        {
            st->old_pitch_buf[i-2] = tmp * st->old_pitch_buf[i];
        }
    }
    if( st->last_bwidth == NB && st->ini_frame != 0 )
    {
        st->rate_switching_reset = 1;
    }

    /*----------------------------------------------------------------*
     * Change the sampling frequency to 12.8 kHz
     *----------------------------------------------------------------*/

    modify_Fs( st->input, input_frame, st->input_Fs, new_inp, 12800, st->mem_decim, 0 );

    /* update signal buffer */
    mvr2r( new_inp, st->buf_speech_enc+L_FRAME, L_FRAME );

    /*------------------------------------------------------------------*
     * Perform fixed preemphasis through 1 - g*z^-1
     *-----------------------------------------------------------------*/

    preemph( new_inp, PREEMPH_FAC, L_FRAME, &st->mem_preemph);

    /*----------------------------------------------------------------*
     * Compute spectrum, find energy per critical frequency band
     * Track energy and signal dynamics
     * Detect NB spectrum in a 16kHz-sampled input
     *----------------------------------------------------------------*/

    analy_sp( inp, st->Bin_E, st->Bin_E_old, fr_bands, lf_E, &Etot, st->min_band, st->max_band, dummy_buf, PS, fft_buff );

    noise_est_pre( Etot, st->ini_frame, &st->Etot_l, &st->Etot_h, &st->Etot_l_lp, &st->Etot_last, &st->Etot_v_h2, &st->sign_dyn_lp, st->harm_cor_cnt, &st->Etot_lp );

    /*----------------------------------------------------------------*
     * VAD
     *----------------------------------------------------------------*/

    vad_flag = wb_vad( st, fr_bands, &localVAD, &tmps, &tmps, &tmps, &snr_sum_he, &localVAD_HE_SAD, &st->flag_noisy_speech_snr );

    if( vad_flag == 0 )
    {
        coder_type = INACTIVE;
    }

    /* apply DTX hangover for CNG analysis */
    vad_flag_dtx = dtx_hangover_addition( st, localVAD, vad_flag, st->lp_speech-st->lp_noise, 0, &vad_hover_flag );


    /*-----------------------------------------------------------------*
     * Select SID or FRAME_NO_DATA frame if DTX enabled
     *-----------------------------------------------------------------*/

    if ( st->last_core != AMR_WB_CORE )
    {
        st->fd_cng_reset_flag = 1;
    }
    else if ( st->fd_cng_reset_flag > 0 && st->fd_cng_reset_flag < 10 )
    {
        st->fd_cng_reset_flag++;
    }
    else
    {
        st->fd_cng_reset_flag = 0;
    }

    dtx( st, vad_flag_dtx, inp );

    /*----------------------------------------------------------------*
     * Noise energy down-ward update and total noise energy estimation
     * Long-term energies and relative frame energy updates
     * Correlation correction as a function of total noise level
     *----------------------------------------------------------------*/

    noise_est_down( fr_bands, st->bckr, tmpN, tmpE, st->min_band, st->max_band, &st->totalNoise, Etot, &st->Etot_last, &st->Etot_v_h2 );

    high_lpn_flag = 0;
    long_enr( st, Etot, localVAD_HE_SAD, high_lpn_flag );
    relE = Etot - st->lp_speech;

    if( st->bwidth != NB )
    {
        lp_bckr = mean( st->bckr, 10 );
    }
    else
    {
        lp_bckr = mean( st->bckr+1, 9 );
    }
    hp_bckr = 0.5f * (st->bckr[st->max_band-1] + st->bckr[st->max_band]);
    st->bckr_tilt_lt = 0.9f * st->bckr_tilt_lt + 0.1f * lp_bckr / hp_bckr;

    corr_shift = correlation_shift(st->totalNoise);

    /*----------------------------------------------------------------*
     * WB, SWB and FB bandwidth detector
     *----------------------------------------------------------------*/

    bw_detect( st, st->input, localVAD, NULL );

    /* in AMR_WB IO, limit the maximum band-width to WB */
    if( st->bwidth > WB )
    {
        st->bwidth = WB;
    }

    /*----------------------------------------------------------------*
     * Perform LP analysis
     * Compute weighted inp
     * Perform open-loop pitch analysis
     * Perform 1/4 pitch precision improvement
     *----------------------------------------------------------------*/

    if ( vad_flag == 0 )
    {
        /* reset the OL pitch tracker memories during inactive frames */
        pitch_ol_init( &st->old_thres, &st->old_pitch, &st->delta_pit, &st->old_corr ) ;
    }

    /* LP analysis */
    analy_lp_AMR_WB( inp, &ener, A, epsP, isp_new, st->lsp_old1, isf_new, st->old_pitch_la, st->old_voicing_la );

    /* compute weighted input */
    find_wsp( L_FRAME, L_SUBFR, NB_SUBFR, A, Aw, inp, TILT_FAC, wsp, &st->mem_wsp, GAMMA1, L_LOOK_12k8 );

    /* open-loop pitch analysis */
    pitch_ol( pitch,voicing, &st->old_pitch, &st->old_corr, corr_shift, &st->old_thres,
              &st->delta_pit, st->old_wsp2, wsp, st->mem_decim2, relE, L_LOOK_12k8, 0, st->bwidth, 0 );

    st->old_pitch_la = pitch[2];
    st->old_voicing_la = voicing[2];


    vad_param_updt( st, pitch, voicing, corr_shift, A );

    /*------------------------------------------------------------------*
     * Update estimated noise energy and voicing cut-off frequency
     *-----------------------------------------------------------------*/

    noise_est( st, tmpN, pitch, voicing, epsP, Etot, relE, corr_shift, tmpE, fr_bands, &cor_map_sum,
               &sp_div, &non_staX, &harm_flag, lf_E, &st->harm_cor_cnt, st->Etot_l_lp, &sp_floor );

    /*----------------------------------------------------------------*
     * Change the sampling frequency to 16 kHz,
     *   input@16kHz needed for AMR-WB IO BWE @23.85kbps
     *----------------------------------------------------------------*/

    if ( st->input_Fs == 16000 )
    {
        /* no resampling needed, only delay adjustement to account for the FIR resampling delay */
        tmps = NS2SA(16000, DELAY_FIR_RESAMPL_NS);
        mvr2r( st->mem_decim16k + tmps, new_inp_16k, tmps );
        mvr2r( st->input, new_inp_16k + tmps, input_frame - tmps );
        mvr2r( st->input + input_frame - 2*tmps, st->mem_decim16k, 2*tmps );
    }
    else if( st->input_Fs == 32000 || st->input_Fs == 48000 )
    {
        modify_Fs( st->input, input_frame, st->input_Fs, new_inp_16k, 16000, st->mem_decim16k, 0 );
    }

    /*----------------------------------------------------------------*
     * Encoding of SID frames
     *----------------------------------------------------------------*/

    if ( st->core_brate == SID_1k75 || st->core_brate == FRAME_NO_DATA )
    {
        /* encode CNG parameters */
        CNG_enc( st, L_FRAME, Aq, inp, ener, isp_new, isf_new , &allow_cn_step, st->burst_ho_cnt, q_env, &sid_bw, st->exc_mem2 );

        /* comfort noise generation */
        CNG_exc( st->core_brate, L_FRAME, &st->Enew, &st->cng_seed, exc, exc2, &st->lp_ener,
                 st->last_core_brate, &st->first_CNG, &st->cng_ener_seed, dummy_buf, allow_cn_step, &st->last_allow_cn_step, st->num_ho,
                 q_env, st->lp_env, st->old_env, st->exc_mem, st->exc_mem1, &sid_bw, &st->cng_ener_seed1, exc3, st->Opt_AMR_WB );

        if ( st->first_CNG == 0 )
        {
            st->first_CNG = 1;
        }
        /* synthesis */
        syn_12k8( L_FRAME, Aq, exc2, dummy_buf, st->LPDmem.mem_syn3, 1 );  /* dummy_buf = temporary buffer to handle syn1[] */

        /* reset the encoder */
        CNG_reset_enc( st, &(st->LPDmem), pitch_buf, dummy_buf+L_FRAME, 0 );

        /* update st->mem_syn1 for ACELP core switching */
        mvr2r( st->LPDmem.mem_syn3, st->mem_syn1, M );

        /* update ACELP core synthesis filter memory */
        mvr2r( st->LPDmem.mem_syn3, st->LPDmem.mem_syn, M );

        /* update old synthesis buffer - needed for ACELP internal sampling rate switching */
        mvr2r( dummy_buf + L_FRAME - L_SYN_MEM, st->LPDmem.mem_syn_r, L_SYN_MEM );

        /* Update MODE2 core switching memory */
        deemph( dummy_buf, PREEMPH_FAC, L_FRAME, &(st->LPDmem.syn[M]) );
        mvr2r( dummy_buf+L_FRAME-M-1, st->LPDmem.syn, M+1 );
    }

    /*----------------------------------------------------------------*
     * Encoding of all other frames
     *----------------------------------------------------------------*/

    else
    {
        /*-----------------------------------------------------------------*
         * After inactive period, use the most up-to-date ISPs
         *-----------------------------------------------------------------*/

        if( st->last_core_brate == FRAME_NO_DATA || st->last_core_brate == SID_1k75 )
        {
            mvr2r( st->lspCNG, st->lsp_old, M );
            isp2isf( st->lspCNG, st->lsf_old, M, INT_FS_12k8 );
            set_f( old_exc, 0, L_EXC_MEM );
        }

        /*-----------------------------------------------------------------*
         * ISF Quantization and interpolation
         *-----------------------------------------------------------------*/

        isf_enc_amr_wb( st, isf_new, isp_new, Aq, &st->stab_fac );

        /*---------------------------------------------------------------*
         * Calculation of LP residual (filtering through A[z] filter)
         *---------------------------------------------------------------*/

        calc_residu( inp, res, Aq, L_FRAME );
        st->burst_ho_cnt = 0;

        /*------------------------------------------------------------*
         * Encode excitation
         *------------------------------------------------------------*/

        encod_amr_wb( st, &(st->LPDmem), inp, Aw, Aq, pitch, voicing, res, syn, exc, exc2, pitch_buf, hf_gain, inp_16k );


        /* update st->mem_syn1 for ACELP core switching */
        mvr2r( st->LPDmem.mem_syn, st->mem_syn1, M );

        /* update old synthesis buffer - needed for ACELP internal sampling rate switching */
        mvr2r( syn + L_FRAME - L_SYN_MEM, st->LPDmem.mem_syn_r, L_SYN_MEM );

        /* Update MODE2 core switching memory */
        mvr2r( syn, dummy_buf, L_FRAME );
        deemph( dummy_buf, PREEMPH_FAC, L_FRAME, &(st->LPDmem.syn[M]) );
        mvr2r( dummy_buf+L_FRAME-M-1, st->LPDmem.syn, M+1 );

        /*--------------------------------------------------------------------------------------*
         * Write VAD information into the bitstream in AMR-WB IO mode
         *--------------------------------------------------------------------------------------*/

        push_indice( st, IND_VAD_FLAG, vad_flag, 1 );

    }


    /*-----------------------------------------------------------------*
     * Updates
     *-----------------------------------------------------------------*/

    /* update old weighted speech buffer - for OL pitch analysis */
    mvr2r( &old_wsp[L_FRAME], st->old_wsp, L_WSP_MEM );

    /* update old input signal buffer */
    mvr2r( &old_inp[L_FRAME], st->old_inp_12k8, L_INP_MEM );

    /* update old input signal @16kHz buffer */
    if( st->input_Fs > 8000 )
    {
        mvr2r( &old_inp_16k[L_FRAME16k], st->old_inp_16k, L_INP_MEM );
    }

    /* update of old per-band energy spectrum */
    mvr2r( fr_bands + NB_BANDS, st->enrO, NB_BANDS );

    /* update the last bandwidth */
    st->last_input_bwidth = st->input_bwidth;
    st->last_bwidth = st->bwidth;

    /* update signal buffers */
    mvr2r( new_inp, st->buf_speech_enc_pe+L_FRAME, L_FRAME );
    mvr2r( wsp, st->buf_wspeech_enc+L_FRAME+L_SUBFR, L_FRAME + L_LOOK_12k8 );

    updt_enc( st, L_FRAME, coder_type, old_exc, pitch_buf, 0, Aq, isf_new, isp_new, dummy_buf );

    core_encode_update( st );

    /* update main codec parameters */
    st->last_extl = -1;
    st->last_core = st->core;
    st->last_L_frame = L_FRAME;
    st->last_core_brate = st->core_brate;
    st->last_total_brate = st->total_brate;
    st->Etot_last = Etot;
    st->last_coder_type_raw = st->coder_type_raw;
    st->last_codec_mode = st->codec_mode;

    /* Increase the counter of initialization frames */
    if( st->ini_frame < MAX_FRAME_COUNTER )
    {
        (st->ini_frame)++;
    }

    if( st->core_brate > SID_1k75 )
    {
        st->last_active_brate = st->total_brate;
    }

    if ( st->core_brate > SID_1k75 && st->first_CNG )
    {
        if( st->act_cnt >= BUF_DEC_RATE )
        {
            st->act_cnt = 0;
        }

        st->act_cnt++;

        if( st->act_cnt == BUF_DEC_RATE && st->ho_hist_size > 0 )
        {
            st->ho_hist_size--;
        }

        if( ++(st->act_cnt2) >= MIN_ACT_CNG_UPD )
        {
            st->act_cnt2 = MIN_ACT_CNG_UPD;
        }
    }




    return;
}
