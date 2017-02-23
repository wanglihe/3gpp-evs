/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <assert.h>
#include <string.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"


/*-----------------------------------------------------------------------*
 * Local functions
 *-----------------------------------------------------------------------*/

static void init_tcx( Encoder_State *st, int L_frame_old );
static void init_sig_buffers( Encoder_State *st, const int L_frame_old, const short L_subfr );
static void init_core_sig_ana( Encoder_State *st );
static void init_acelp( Encoder_State *st, int L_frame_old );
static void init_modes( Encoder_State *st );


/*-----------------------------------------------------------------------*
 * init_coder_ace_plus()
 *
 * Initialization of state variables
 *-----------------------------------------------------------------------*/

void init_coder_ace_plus( Encoder_State *st )
{
    int L_frame_old; /*keep old frame size for switching */
    short L_subfr;

    /* Bitrate */
    st->tcxonly = getTcxonly(st->total_brate);

    /* Core Sampling Rate */
    st->sr_core = getCoreSamplerateMode2(st->total_brate, st->bwidth, st->rf_mode );
    st->fscale  = sr2fscale(st->sr_core);

    /* Narrowband? */
    st->narrowBand = (st->bwidth == NB)?1:0;

    /* Core Framing */
    L_frame_old = st->last_L_frame;
    st->L_frame = st->sr_core / 50;
    st->L_frame_past = -1;

    st->L_frameTCX = st->input_Fs / 50;

    if( st->L_frame == L_FRAME16k && st->total_brate <= ACELP_32k )
    {
        st->nb_subfr = NB_SUBFR16k;
    }
    else
    {
        st->nb_subfr = NB_SUBFR;
    }
    L_subfr = st->L_frame/st->nb_subfr;

    /* Core Lookahead */
    st->encoderLookahead_enc = NS2SA(st->sr_core, ACELP_LOOK_NS);
    st->encoderLookahead_FB = NS2SA(st->input_Fs, ACELP_LOOK_NS);

    if( st->ini_frame == 0 )
    {
        st->acelpFramesCount = 0;
        st->prevTempFlatness = 1.0f;
    }

    /* Initialize TBE */
    st->prev_coder_type = GENERIC;
    set_f( st->prev_lsf_diff, 0.5f, LPC_SHB_ORDER-2 );
    st->prev_tilt_para = 0.0f;
    set_zero( st->cur_sub_Aq, M+1 );

    st->currEnergyHF = 0;

    /* Initialize LPC analysis/quantization */
    if( st->sr_core <= 16000 && st->tcxonly == 0 )
    {
        st->lpcQuantization = 1;
    }
    else
    {
        st->lpcQuantization = 0;
    }

    st->next_force_safety_net = 0;
    if( (st->last_L_frame != st->L_frame) || (st->last_core == AMR_WB_CORE) || (st->last_core == HQ_CORE) )
    {
        set_f( st->mem_MA, 0.0f, M );
        mvr2r( GEWB_Ave, st->mem_AR, M );
    }

    /* Initialize IGF */
    memset( &st->hIGFEnc, 0, sizeof(st->hIGFEnc) );
    st->hIGFEnc.infoStopFrequency = -1;
    if( st->igf )
    {
        IGFEncSetMode( &st->hIGFEnc, st->total_brate, st->bwidth, st->rf_mode );
    }

    /* Initialize TCX */
    init_tcx( st, L_frame_old );

    /* Initialize Core Signal Analysis Module */
    init_core_sig_ana( st );

    /* Initialize Signal Buffers */
    init_sig_buffers( st, L_frame_old, L_subfr );
    /* Initialize ACELP */
    init_acelp( st, L_frame_old );
    if( st->ini_frame == 0 )
    {
        st->tec_tfa = 0;
    }
    if( st->tec_tfa == 0 )
    {
        resetTecEnc( &(st->tecEnc), 0);
    }
    else
    {
        resetTecEnc( &(st->tecEnc), 1);
    }

    if( st->bwidth == SWB && (st->total_brate == ACELP_16k40 || st->total_brate == ACELP_24k40) )
    {
        st->tec_tfa = 1;
    }
    else
    {
        st->tec_tfa = 0;
    }

    st->tec_flag = 0;
    st->tfa_flag = 0;

    /* Initialize DTX */
    if( st->ini_frame == 0 )
    {

        vad_init(&st->vad_st);
    }

    if( st->total_brate == ACELP_9k60 || st->total_brate == ACELP_16k40 || st->total_brate == ACELP_24k40 )
    {
        st->glr = 1;
    }
    else
    {
        st->glr = 0;
    }

    st->glr_reset = 0;

    /* Initialize ACELP/TCX Modes */
    init_modes( st );

    /* Init I/O */


    /* Adaptive BPF */
    set_zero( st->mem_bpf, 2*L_FILT16k );
    set_zero( st->mem_error_bpf, 2*L_FILT16k );

    if( st->total_brate >= HQ_48k )
    {
        st->enablePlcWaveadjust = 1;
    }
    else
    {
        st->enablePlcWaveadjust = 0;
    }

    open_PLC_ENC_EVS( &st->plcExt, st->sr_core );

    st->glr_idx[0] = 0;
    st->glr_idx[1] = 0;
    st->mean_gc[0] = 0.0f;
    st->mean_gc[1] = 0.0f;
    st->prev_lsf4_mean = 0.0f;
    st->last_stab_fac = 0.0f;

    return;
}


/*-----------------------------------------------------------------------*
 * init_tcx()
 *
 * Initialization of TCX
 *-----------------------------------------------------------------------*/

static void init_tcx(
    Encoder_State *st,
    int L_frame_old
)
{
    short mdctWindowLength;
    short mdctWindowLengthFB;

    /* Share the memories for 2xTCX10/4xTCX5 and for TCX20 */
    st->spectrum[0] = st->spectrum_long;
    st->spectrum[1] = st->spectrum_long + N_TCX10_MAX;

    st->tcx_cfg.tcx5Size = NS2SA(st->sr_core, FRAME_SIZE_NS/4); /* Always 5 ms */
    st->tcx_cfg.tcx5SizeFB = NS2SA(st->input_Fs, FRAME_SIZE_NS/4); /* Always 5 ms */

    st->tcx_cfg.tcx_mdct_window_length_old = st->tcx_cfg.tcx_mdct_window_length;
    mdctWindowLength = getMdctWindowLength(st->fscale);
    mdctWindowLengthFB = mdctWindowLength * st->input_Fs / st->sr_core;

    /* Initialize the TCX MDCT window */
    /*Symmetric window = sinus LD window*/
    st->tcx_cfg.tcx_mdct_window_delay  = mdctWindowLength;
    st->tcx_cfg.tcx_mdct_window_delayFB  = mdctWindowLengthFB;
    st->tcx_cfg.tcx_mdct_window_length = mdctWindowLength;
    st->tcx_cfg.tcx_mdct_window_lengthFB  = mdctWindowLengthFB;

    mdct_window_sine( st->tcx_cfg.tcx_mdct_window, st->tcx_cfg.tcx_mdct_window_length );
    mdct_window_sine( st->tcx_cfg.tcx_mdct_windowFB, st->tcx_cfg.tcx_mdct_window_lengthFB );

    mdct_window_sine( st->tcx_cfg.tcx_mdct_window_half, st->tcx_cfg.tcx_mdct_window_length/2 );
    mdct_window_sine( st->tcx_cfg.tcx_mdct_window_halfFB, st->tcx_cfg.tcx_mdct_window_lengthFB/2 );

    /*ALDO windows for MODE2*/
    mdct_window_aldo(st->tcx_cfg.tcx_aldo_window_1, st->tcx_cfg.tcx_aldo_window_2, st->L_frame);
    mdct_window_aldo(st->tcx_cfg.tcx_aldo_window_1_FB, st->tcx_cfg.tcx_aldo_window_2_FB, NS2SA(st->input_Fs, FRAME_SIZE_NS));
    st->tcx_cfg.tcx_aldo_window_1_trunc = st->tcx_cfg.tcx_aldo_window_1 + NS2SA(st->sr_core, N_ZERO_MDCT_NS);
    st->tcx_cfg.tcx_aldo_window_1_FB_trunc = st->tcx_cfg.tcx_aldo_window_1_FB + NS2SA(st->input_Fs, N_ZERO_MDCT_NS);

    /*1.25ms transition window for ACELP->TCX*/
    st->tcx_cfg.tcx_mdct_window_trans_length = NS2SA(st->sr_core, ACELP_TCX_TRANS_NS);
    mdct_window_sine( st->tcx_cfg.tcx_mdct_window_trans, st->tcx_cfg.tcx_mdct_window_trans_length );
    st->tcx_cfg.tcx_mdct_window_trans_lengthFB = NS2SA(st->input_Fs, ACELP_TCX_TRANS_NS);
    mdct_window_sine( st->tcx_cfg.tcx_mdct_window_transFB, st->tcx_cfg.tcx_mdct_window_trans_lengthFB );

    /*Mid-OLA*/
    /*compute minimum length for "half" window: lookahead - 5ms. It must be also multiple of 2*/
    st->tcx_cfg.tcx_mdct_window_half_length=2*((st->encoderLookahead_enc-(int)(0.005f*st->sr_core+0.5f))>>1);
    st->tcx_cfg.tcx_mdct_window_half_lengthFB=2*((st->encoderLookahead_FB-(int)(0.005f*st->input_Fs+0.5f))>>1);
    assert( (st->tcx_cfg.tcx_mdct_window_half_length>16) && "Half window can not be large enough!");

    mdct_window_sine( st->tcx_cfg.tcx_mdct_window_half, st->tcx_cfg.tcx_mdct_window_half_length );
    mdct_window_sine( st->tcx_cfg.tcx_mdct_window_halfFB, st->tcx_cfg.tcx_mdct_window_half_lengthFB );

    /* minimum overlap 1.25 ms */
    st->tcx_cfg.tcx_mdct_window_min_length = st->sr_core / 800;
    st->tcx_cfg.tcx_mdct_window_min_lengthFB = st->input_Fs / 800;
    mdct_window_sine(st->tcx_cfg.tcx_mdct_window_minimum, st->tcx_cfg.tcx_mdct_window_min_length);
    mdct_window_sine(st->tcx_cfg.tcx_mdct_window_minimumFB, st->tcx_cfg.tcx_mdct_window_min_lengthFB);

    /* TCX Offset */
    st->tcx_cfg.tcx_offset = (st->tcx_cfg.tcx_mdct_window_delay>>1);
    st->tcx_cfg.tcx_offsetFB = (st->tcx_cfg.tcx_mdct_window_delayFB>>1);
    /*<0 rectangular transition with optimized window size = L_frame+L_frame/4*/
    st->tcx_cfg.lfacNext = st->tcx_cfg.tcx_offset - st->L_frame/4;
    st->tcx_cfg.lfacNextFB = st->tcx_cfg.tcx_offsetFB - st->L_frameTCX/4;

    if( st->ini_frame == 0 )
    {
        st->tcx_cfg.tcx_curr_overlap_mode = st->tcx_cfg.tcx_last_overlap_mode = ALDO_WINDOW;
    }
    /* Init TCX target bits correction factor */
    st->LPDmem.tcx_target_bits_fac = 1.0f;

    st->measuredBwRatio  = 1.0f;
    st->noiseTiltFactor  = 0.5625f;
    st->noiseLevelMemory = 0;

    /*SQ deadzone & memory quantization*/
    st->tcx_cfg.sq_rounding = 0.375f; /*deadzone of 1.25->rounding=1-1.25/2 (No deadzone=0.5)*/
    set_i( st->memQuantZeros, 0, L_FRAME_PLUS );

    /* TCX rate loop */
    st->tcx_cfg.tcxRateLoopOpt = (st->tcxonly) ? 2 : 0;

    /* TCX bandwidth */
    st->tcx_cfg.bandwidth = getTcxBandwidth(st->bwidth);

    /* set number of coded lines */
    st->tcx_cfg.tcx_coded_lines = getNumTcxCodedLines(st->bwidth);

    /* TNS in TCX */
    st->tcx_cfg.pCurrentTnsConfig = NULL;
    st->tcx_cfg.fIsTNSAllowed = getTnsAllowed( st->total_brate, st->igf );

    if( st->tcx_cfg.fIsTNSAllowed )
    {
        InitTnsConfigs( bwMode2fs[st->bwidth], st->tcx_cfg.tcx_coded_lines, st->tcx_cfg.tnsConfig, (&st->hIGFEnc)->infoStopFrequency, st->total_brate);
    }

    /* TCX-LTP */
    st->tcxltp = getTcxLtp(st->sr_core);

    if( st->ini_frame == 0 )
    {
        st->tcxltp_pitch_int_past = st->L_frame;
        st->tcxltp_pitch_fr_past = 0;
        st->tcxltp_gain_past = 0.f;
        st->tcxltp_norm_corr_past = 0.f;
    }
    else if( st->L_frame!=L_frame_old && !((st->total_brate==ACELP_16k40||st->total_brate==ACELP_24k40)&&(st->total_brate==st->last_total_brate)&&(st->last_bwidth==st->bwidth)) )
    {
        int pitres, pitres_old;
        float pit, pit_old;

        if ( L_frame_old%160==0 )
        {
            pitres_old = 6;
        }
        else
        {
            pitres_old = 4;
        }
        pit_old = (float)st->tcxltp_pitch_int_past + (float)st->tcxltp_pitch_fr_past/(float)pitres_old;
        if (st->L_frame%160==0 )
        {
            pitres = 6;
        }
        else
        {
            pitres = 4;
        }
        pit = pit_old * (float)st->L_frame/(float)L_frame_old;
        st->tcxltp_pitch_int_past = (int)pit;
        st->tcxltp_pitch_fr_past = (int)( (pit-(float)st->tcxltp_pitch_int_past)*(float)pitres );
    }

    /* Context HM*/
    st->tcx_cfg.ctx_hm = getCtxHm( st->total_brate, st->rf_mode );

    /* Residual Coding*/
    st->tcx_cfg.resq = getResq(st->total_brate);
    st->tcx_cfg.tcxRateLoopOpt = (st->tcx_cfg.resq && !st->tcxonly) ? 1 : st->tcx_cfg.tcxRateLoopOpt;

    st->tcx_lpc_shaped_ari = getTcxLpcShapedAri( st->total_brate, st->bwidth, st->rf_mode );

    return;

}

/*-----------------------------------------------------------------------*
 * init_sig_buffers()
 *
 * Initialization of signal buffers
 *-----------------------------------------------------------------------*/

static void init_sig_buffers(
    Encoder_State *st,
    const int L_frame_old,
    const short L_subfr
)
{
    /* Encoder Past Samples at encoder-sampling-rate */
    st->encoderPastSamples_enc = (st->L_frame*9)/16;

    /* Initialize Signal Buffers and Pointers at encoder-sampling-rate */
    if ( st->ini_frame == 0 )
    {
        set_zero(st->buf_speech_enc, L_PAST_MAX_32k+L_FRAME32k+L_NEXT_MAX_32k);
        set_zero(st->buf_speech_enc_pe, L_PAST_MAX_32k+L_FRAME32k+L_NEXT_MAX_32k);
        set_zero(st->buf_speech_ltp, L_PAST_MAX_32k+L_FRAME32k+L_NEXT_MAX_32k);
        set_zero(st->buf_wspeech_enc, L_FRAME16k+L_SUBFR+L_FRAME16k+L_NEXT_MAX_16k);
    }
    else if ( st->L_frame!=L_frame_old && !((st->total_brate==16400||st->total_brate==24400)&&(st->total_brate==st->last_total_brate)&&(st->last_bwidth==st->bwidth)) )
    {
        lerp( st->buf_speech_enc, st->buf_speech_enc, st->L_frame, L_frame_old );

        if( (st->last_core != TCX_20_CORE) && (st->last_core != TCX_10_CORE) )
        {
            mvr2r( st->buf_speech_enc, st->buf_speech_ltp, st->L_frame );
        }

        mvr2r( st->old_wsp, st->buf_wspeech_enc+st->L_frame + L_SUBFR-L_WSP_MEM,L_WSP_MEM);

        /*Resamp buffers needed only for ACELP*/
        if( st->L_frame == L_FRAME && !st->tcxonly )
        {
            mvr2r( st->old_inp_12k8, st->buf_speech_enc_pe+st->L_frame-L_INP_MEM,L_INP_MEM);
        }
        else if( st->L_frame == L_FRAME16k && !st->tcxonly )
        {
            lerp( st->buf_wspeech_enc+st->L_frame + L_SUBFR-L_WSP_MEM, st->buf_wspeech_enc+st->L_frame + L_SUBFR-310, 310, L_WSP_MEM );
            mvr2r( st->old_inp_16k, st->buf_speech_enc_pe+st->L_frame-L_INP_MEM,L_INP_MEM);
        }

        st->mem_preemph_enc = st->buf_speech_enc[st->L_frame-1];
        st->mem_wsp_enc = st->buf_wspeech_enc[st->L_frame+L_SUBFR-1];

    }
    /*coming from TCXonly modes*/
    else if( !st->tcxonly && st->last_total_brate>ACELP_32k )
    {
        mvr2r( st->old_wsp, st->buf_wspeech_enc+st->L_frame + L_SUBFR-L_WSP_MEM,L_WSP_MEM);

        /*Resamp buffers needed only for ACELP*/
        if( st->L_frame == L_FRAME16k)
        {
            lerp( st->buf_wspeech_enc+st->L_frame + L_SUBFR-L_WSP_MEM, st->buf_wspeech_enc+st->L_frame + L_SUBFR-310, 310, L_WSP_MEM );
        }
        st->LPDmem.mem_w0 = 0;
        st->mem_wsp_enc = st->buf_wspeech_enc[st->L_frame+L_SUBFR-1];
    }

    st->new_speech_enc      = st->buf_speech_enc      + st->encoderPastSamples_enc    + st->encoderLookahead_enc;
    st->new_speech_enc_pe   = st->buf_speech_enc_pe   + st->encoderPastSamples_enc    + st->encoderLookahead_enc;
    st->new_speech_ltp      = st->buf_speech_ltp      + st->encoderPastSamples_enc    + st->encoderLookahead_enc;
    st->new_speech_TCX      = st->input_buff + L_FRAME48k + NS2SA(48000, DELAY_FIR_RESAMPL_NS) - NS2SA(st->input_Fs, DELAY_FIR_RESAMPL_NS);

    st->speech_enc          = st->buf_speech_enc      + st->encoderPastSamples_enc;
    st->speech_enc_pe       = st->buf_speech_enc_pe   + st->encoderPastSamples_enc;
    st->speech_ltp          = st->buf_speech_ltp      + st->encoderPastSamples_enc;
    st->speech_TCX          = st->new_speech_TCX      - st->encoderLookahead_FB;

    st->wspeech_enc         = st->buf_wspeech_enc     + st->L_frame + L_subfr;

    if( st->ini_frame == 0 || st->L_frame != L_frame_old || st->last_codec_mode == MODE1 )
    {
        set_zero( st->buf_synth, OLD_SYNTH_SIZE_ENC+L_FRAME32k );
    }
    st->synth               = st->buf_synth           + st->L_frame + L_subfr;

    return;

}


/*-----------------------------------------------------------------------*
 * init_core_sig_ana()
 *
 *
 *-----------------------------------------------------------------------*/

static void init_core_sig_ana(
    Encoder_State *st
)
{
    /* Pre-emphasis factor and memory */
    if (st->fscale < (16000*FSCALE_DENOM)/12800)
    {
        st->preemph_fac = PREEMPH_FAC; /*WB*/
    }
    else if (st->fscale < (24000*FSCALE_DENOM)/12800)
    {
        st->preemph_fac = PREEMPH_FAC_16k; /*WB*/
    }
    else
    {
        st->preemph_fac = PREEMPH_FAC_SWB; /*SWB*/
    }

    st->tcx_cfg.preemph_fac = st->preemph_fac;

    if ( st->sr_core==16000 )
    {
        st->gamma = GAMMA16k;
    }
    else
    {
        st->gamma = GAMMA1;
    }


    if( st->narrowBand )
    {
        st->min_band = 1;
        st->max_band = 16;
    }
    else
    {
        st->min_band = 0;
        st->max_band = 19;
    }

    return;
}


/*-----------------------------------------------------------------------*
 * init_acelp()
 *
 *
 *-----------------------------------------------------------------------*/

static void init_acelp(
    Encoder_State *st,
    int L_frame_old
)
{
    short mem_syn_r_size_old;
    short mem_syn_r_size_new;

    /* Init pitch lag */
    st->pit_res_max = initPitchLagParameters(st->sr_core, &st->pit_min, &st->pit_fr1, &st->pit_fr1b, &st->pit_fr2, &st->pit_max);

    /* Init LPDmem */
    if (st->ini_frame == 0)
    {
        set_zero( st->LPDmem.syn, 1+M );

        set_zero( st->LPDmem.Txnq, L_FRAME32k/2+64 );
        st->LPDmem.acelp_zir = st->LPDmem.Txnq + (st->L_frame/2);
        set_zero( st->LPDmem.mem_syn_r, L_SYN_MEM );
    }
    else /*Rate switching*/
    {
        if( st->last_core == ACELP_CORE )
        {
            lerp( st->LPDmem.Txnq,st->LPDmem.Txnq, st->L_frame/2, L_frame_old/2 );
        }
        else
        {
            lerp( st->LPDmem.Txnq, st->LPDmem.Txnq, st->tcx_cfg.tcx_mdct_window_length, st->tcx_cfg.tcx_mdct_window_length_old );
        }
        st->LPDmem.acelp_zir = st->LPDmem.Txnq + (st->L_frame/2);

        /* Rate switching */
        if( st->last_codec_mode == MODE1 )
        {
            mvr2r( st->mem_syn1, st->LPDmem.mem_syn2, M );
            set_zero( st->LPDmem.Txnq, L_FRAME32k/2+64 );
            set_zero( st->LPDmem.syn, M );
        }

        if( st->last_core == AMR_WB_CORE )
        {
            st->next_force_safety_net = 1;
            st->last_core = ACELP_CORE;
        }

        if( st->last_codec_mode == MODE1 && st->last_core == HQ_CORE )
        {
            /*Reset of ACELP memories*/
            st->next_force_safety_net = 1;
            st->rate_switching_reset = 1;
            st->LPDmem.tilt_code = TILT_CODE;
            set_zero( st->LPDmem.old_exc, L_EXC_MEM );
            set_zero( st->LPDmem.syn, 1+M );
            st->LPDmem.mem_w0 = 0.0f;
            set_zero( st->LPDmem.mem_syn, M );
            set_zero( st->LPDmem.mem_syn2, M );

            /* unquantized LPC*/
            if ( !((st->total_brate == ACELP_16k40 || st->total_brate == ACELP_24k40) && st->total_brate == st->last_total_brate && st->last_bwidth == st->bwidth) )
            {
                mvr2r( st->lsp_old1, st->lspold_enc, M ); /*lsp old @12.8kHz*/
                if( st->L_frame == L_FRAME16k )
                {
                    lsp_convert_poly( st->lspold_enc, st->L_frame, 0 );
                }
            }
            mvr2r( st->lspold_enc, st->lsp_old, M ); /*used unquantized values for mid-LSF Q*/
            lsp2lsf( st->lsp_old, st->lsf_old, M, st->sr_core );

            st->last_core = TCX_20_CORE;

            st->tcx_cfg.last_aldo=1;  /*It was previously ALDO*/
            st->tcx_cfg.tcx_curr_overlap_mode = ALDO_WINDOW;

            /*ALDO overlap windowed past: also used in MODE1 but for other MDCT-FB*/
            set_f( st->old_out, 0, st->L_frame );

        }
        else
        {
            if( st->L_frame != L_frame_old && st->L_frame <= L_FRAME16k && L_frame_old <= L_FRAME16k ) /* Rate switching between 12.8 and 16 kHz*/
            {
                float tmp, A[M+1], Ap[M+1],tmp_buf[M+1];

                /* convert quantized LSP vector */
                st->rate_switching_reset = lsp_convert_poly( st->lsp_old, st->L_frame, 0 );
                lsp2lsf( st->lsp_old, st->lsf_old, M, st->sr_core );

                if( st->L_frame == L_FRAME16k )
                {
                    mvr2r( st->lsp_old, st->lspold_enc, M );
                }
                else
                {
                    mvr2r( st->lsp_old1, st->lspold_enc, M );
                }

                synth_mem_updt2( st->L_frame, st->last_L_frame, st->LPDmem.old_exc, st->LPDmem.mem_syn_r, st->LPDmem.mem_syn2, st->LPDmem.mem_syn, ENC );

                /* Update wsyn */
                lsp2a_stab( st->lsp_old, A, M );
                weight_a( A, Ap, GAMMA1, M );
                tmp = 0.f;
                tmp_buf[0] = 0.f;
                mvr2r( st->LPDmem.mem_syn2, tmp_buf+1, M );
                deemph( tmp_buf+1, st->preemph_fac, M, &tmp );
                residu( Ap, M, tmp_buf+M, &tmp, 1 );
                st->LPDmem.mem_w0 = st->wspeech_enc[-1] - tmp;
            }
            else if( st->L_frame != L_frame_old )   /* Rate switching involving TCX only modes */
            {
                /*Partial reset of ACELP memories*/
                st->next_force_safety_net = 1;
                st->rate_switching_reset = 1;

                /*reset partly some memories*/
                st->LPDmem.tilt_code = TILT_CODE;
                set_zero( st->LPDmem.old_exc, L_EXC_MEM );

                /*Resamp others memories*/
                /*Size of LPC syn memory*/
                mem_syn_r_size_old = (int)(1.25*L_frame_old/20.f);
                mem_syn_r_size_new = (int)(1.25*st->L_frame/20.f);
                lerp( st->LPDmem.mem_syn_r+L_SYN_MEM-mem_syn_r_size_old, st->LPDmem.mem_syn_r+L_SYN_MEM-mem_syn_r_size_new, mem_syn_r_size_new, mem_syn_r_size_old );
                mvr2r( st->LPDmem.mem_syn_r+L_SYN_MEM-M, st->LPDmem.mem_syn, M );
                mvr2r( st->LPDmem.mem_syn, st->LPDmem.mem_syn2, M );

                /*Untouched memories : st->LPDmem.syn & st->LPDmem.mem_w0*/
                st->LPDmem.mem_w0 = 0;

                /* unquantized LPC*/
                mvr2r( st->lsp_old1, st->lspold_enc, M ); /*lsp old @12.8kHz*/
                if( st->L_frame == L_FRAME16k )
                {
                    lsp_convert_poly( st->lspold_enc, st->L_frame, 0 );
                }
                mvr2r( st->lspold_enc, st->lsp_old, M ); /*used unquantized values for mid-LSF Q*/
                lsp2lsf( st->lsp_old, st->lsf_old, M, st->sr_core );
            }
            /* necessary in BASOP only, due to different representations of st->lsf_old */
            /* else if ( !st->tcxonly && (st->L_frame == L_FRAME16k) && (st->last_total_brate > ACELP_32k) ) */
            /* { */
            /* 	lsp2lsf( st->lsp_old, st->lsf_old, M, st->sr_core ); */
            /* } */
        }
    }

    if( st->last_bwidth == NB && st->bwidth != NB && st->ini_frame != 0 )
    {
        st->rate_switching_reset = 1;
    }

    /* Post-processing */
    set_zero( st->dispMem, 8 );
    st->LPDmem.gc_threshold = 0.0f;

    /* Pulse Search configuration */
    st->acelp_autocorr = 1;

    /*Use for 12.8 kHz sampling rate and low bitrates, the conventional pulse search->better SNR*/
    if( ((st->total_brate <= ACELP_9k60 || st->rf_mode == 1) && (st->sr_core == 12800)) )
    {
        st->acelp_autocorr = 0;
    }

    /*BPF parameters for adjusting gain in function of background noise*/
    if( st->codec_mode == MODE2 )
    {
        st->pst_lp_ener = 0.0f;
        if( st->last_codec_mode == MODE1 )
        {
            st->pst_mem_deemp_err = 0.0f;
        }
    }


    return;
}

/*-----------------------------------------------------------------------*
 * init_modes()
 *
 *
 *-----------------------------------------------------------------------*/

static void init_modes(
    Encoder_State *st
)
{
    short n;

    /* Restrict ACE/TCX20/TCX10 mode */
    st->restrictedMode = getRestrictedMode(st->total_brate, st->Opt_AMR_WB);
    st->acelpEnabled = (st->restrictedMode & 1) == 1;
    st->tcx20Enabled = (st->restrictedMode & 2) == 2;
    st->tcx10Enabled = (st->restrictedMode & 4) == 4;

    /* TCX mode (TCX20 TCX10_10 or NO_TCX) */
    st->tcxMode = NO_TCX;

    /*st->bits_frame_nominal = (int)( (float)st->L_frame * (float)FSCALE_DENOM * (float)st->total_brate / ( (float)st->fscale * 12800.0f ) );*/
    st->bits_frame_nominal = (int)( (float)st->L_frame/(float)st->fscale * (float)FSCALE_DENOM/128.0f * (float)st->total_brate/100.0f + 0.49f );

    if( st->Opt_AMR_WB )
    {
        st->bits_frame      = st->bits_frame_nominal;
        st->bits_frame_core = st->bits_frame_nominal;
        st->frame_size_index = 0;
    }
    else
    {
        for (n=0; n<FRAME_SIZE_NB; n++)
        {
            if (FrameSizeConfig[n].frame_bits==st->bits_frame_nominal)
            {
                st->frame_size_index = n;
                st->bits_frame       = FrameSizeConfig[n].frame_bits;
                st->bits_frame_core  = FrameSizeConfig[n].frame_net_bits;
                break;
            }
        }
        if (n==FRAME_SIZE_NB)
        {
            assert(!"Bitrate not supported: not part of EVS");
        }
    }

    /* Reconfigure core */
    core_coder_reconfig( st );

    return;
}
