/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <assert.h>
#include "options.h"
#include "stat_com.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"

/*-----------------------------------------------------------------------*
 * open_decoder_LPD()
 *
 * Initialization of state variables
 *-----------------------------------------------------------------------*/

void open_decoder_LPD(
    Decoder_State *st,
    const int bit_rate,
    const int bandwidth
)
{
    short i;
    short mdctWindowLength;
    short mem_syn_r_size_old;
    short mem_syn_r_size_new;
    short fscaleFB;
    short mdctWindowLengthFB;
    short encoderLookahead, encoderLookaheadFB;

    st->total_brate = bit_rate;
    st->fscale_old  = st->fscale;
    st->sr_core     = getCoreSamplerateMode2(st->total_brate, bandwidth, st->rf_flag);
    st->fscale      = sr2fscale(st->sr_core);
    fscaleFB        = sr2fscale(st->output_Fs);

    /* initializing variables for frame lengths etc. right in the beginning */
    st->L_frame  = st->sr_core / 50;
    st->L_frameTCX = st->output_Fs / 50;
    if (st->ini_frame == 0)
    {
        st->last_L_frame = st->L_frame_past = st->L_frame;
        st->L_frameTCX_past  = st->L_frameTCX;
    }

    st->tcxonly = getTcxonly(bit_rate);

    /* the TD TCX PLC in MODE1 still runs with 80ms subframes */
    if ( (st->L_frame == L_FRAME16k && st->total_brate <= ACELP_32k) || (st->tcxonly && (st->sr_core == 32000 || st->sr_core == 16000)))
    {
        st->nb_subfr = NB_SUBFR16k;
    }
    else
    {
        st->nb_subfr = NB_SUBFR;
    }
    st->bits_frame  = (int) ( ((float)st->L_frame/(float)st->fscale) * (float)FSCALE_DENOM/128.0f * (float)bit_rate /100.0f +0.49f);

    st->TcxBandwidth   = getTcxBandwidth(bandwidth);
    st->narrowBand     = (bandwidth == NB)?1:0;
    encoderLookahead   = (L_LOOK_12k8*st->fscale)/FSCALE_DENOM;
    encoderLookaheadFB = (L_LOOK_12k8*fscaleFB)/FSCALE_DENOM;

    st->pit_res_max = initPitchLagParameters(st->sr_core, &st->pit_min, &st->pit_fr1, &st->pit_fr1b, &st->pit_fr2, &st->pit_max);
    if ( st->ini_frame == 0)
    {
        st->pit_res_max_past = st->pit_res_max;
    }
    st->pit_max_TCX =  st->pit_max * st->output_Fs / st->sr_core;
    st->pit_min_TCX =  st->pit_min * st->output_Fs / st->sr_core;

    /*Preemphasis param*/
    if (st->fscale < (16000*FSCALE_DENOM)/12800)
    {
        st->preemph_fac = PREEMPH_FAC; /*NB*/
    }
    else if (st->fscale < (24000*FSCALE_DENOM)/12800)
    {
        st->preemph_fac = PREEMPH_FAC_16k; /*WB*/
    }
    else
    {
        st->preemph_fac = PREEMPH_FAC_SWB; /*SWB*/
    }

    if( st->sr_core==16000 )
    {
        st->gamma = GAMMA16k;
    }
    else
    {
        st->gamma = GAMMA1;
    }

    /* LPC quantization */
    if ( st->sr_core <= 16000 && st->tcxonly==0 )
    {
        st->lpcQuantization = 1;
    }
    else
    {
        st->lpcQuantization = 0;
    }

    if ( st->tcxonly==0 )
    {
        st->numlpc = 1;
    }
    else
    {
        st->numlpc = 2;
    }

    /* Initialize TBE */
    st->prev_coder_type = GENERIC;
    set_f( st->prev_lsf_diff, 0.5f, LPC_SHB_ORDER-2 );
    st->prev_tilt_para = 0.0f;
    set_zero( st->cur_sub_Aq, M+1 );

    /*TCX config*/
    st->tcx_cfg.preemph_fac=st->preemph_fac;
    st->tcx_cfg.tcx_mdct_window_length_old=st->tcx_cfg.tcx_mdct_window_length;
    mdctWindowLength = getMdctWindowLength(st->fscale);
    mdctWindowLengthFB = getMdctWindowLength(fscaleFB);

    /* Initialize the TCX MDCT window */
    st->tcx_cfg.tcx_mdct_window_delay = mdctWindowLength;
    st->tcx_cfg.tcx_mdct_window_delayFB = mdctWindowLengthFB;
    st->tcx_cfg.tcx_mdct_window_length = mdctWindowLength;
    st->tcx_cfg.tcx_mdct_window_lengthFB = mdctWindowLengthFB;

    mdct_window_sine( st->tcx_cfg.tcx_mdct_window, st->tcx_cfg.tcx_mdct_window_length );
    mdct_window_sine( st->tcx_cfg.tcx_mdct_windowFB, st->tcx_cfg.tcx_mdct_window_lengthFB );
    set_zero(st->tcx_cfg.tcx_mdct_window + st->tcx_cfg.tcx_mdct_window_length, L_MDCT_OVLP_MAX - st->tcx_cfg.tcx_mdct_window_length);
    set_zero(st->tcx_cfg.tcx_mdct_windowFB + st->tcx_cfg.tcx_mdct_window_lengthFB, L_MDCT_OVLP_MAX - st->tcx_cfg.tcx_mdct_window_lengthFB);

    mdct_window_sine( st->tcx_cfg.tcx_mdct_window_half, st->tcx_cfg.tcx_mdct_window_length/2 );
    mdct_window_sine( st->tcx_cfg.tcx_mdct_window_halfFB, st->tcx_cfg.tcx_mdct_window_lengthFB/2 );

    /* minimum overlap 1.25 ms */
    st->tcx_cfg.tcx_mdct_window_min_length = st->sr_core / 800;
    st->tcx_cfg.tcx_mdct_window_min_lengthFB = st->output_Fs / 800;

    mdct_window_sine( st->tcx_cfg.tcx_mdct_window_minimum, st->tcx_cfg.tcx_mdct_window_min_length );
    mdct_window_sine( st->tcx_cfg.tcx_mdct_window_minimumFB, st->tcx_cfg.tcx_mdct_window_min_lengthFB );

    /*ALDO windows for MODE2*/
    mdct_window_aldo(st->tcx_cfg.tcx_aldo_window_1, st->tcx_cfg.tcx_aldo_window_2, st->L_frame);
    mdct_window_aldo(st->tcx_cfg.tcx_aldo_window_1_FB, st->tcx_cfg.tcx_aldo_window_2_FB, NS2SA(st->output_Fs, FRAME_SIZE_NS));
    st->tcx_cfg.tcx_aldo_window_1_trunc = st->tcx_cfg.tcx_aldo_window_1 + NS2SA(st->sr_core, N_ZERO_MDCT_NS);
    st->tcx_cfg.tcx_aldo_window_1_FB_trunc = st->tcx_cfg.tcx_aldo_window_1_FB + NS2SA(st->output_Fs, N_ZERO_MDCT_NS);

    /*1.25ms transition window for ACELP->TCX*/
    mdct_window_sine( st->tcx_cfg.tcx_mdct_window_trans, (st->L_frame>>1)-st->tcx_cfg.tcx_mdct_window_length);
    mdct_window_sine( st->tcx_cfg.tcx_mdct_window_transFB, (st->L_frameTCX>>1)-st->tcx_cfg.tcx_mdct_window_lengthFB);

    /*Mid-OLA*/
    /*compute minimum length for "half" window: lookahead - 5ms. It must be also multiple of 2*/
    st->tcx_cfg.tcx_mdct_window_half_length = 2*((encoderLookahead-(int)(0.005f*st->sr_core+0.5f))>>1);
    st->tcx_cfg.tcx_mdct_window_half_lengthFB = 2*((encoderLookaheadFB-(int)(0.005f*st->output_Fs+0.5f))>>1);
    assert( (st->tcx_cfg.tcx_mdct_window_half_length>16) && "Half window can not be large enough!");

    mdct_window_sine( st->tcx_cfg.tcx_mdct_window_half, st->tcx_cfg.tcx_mdct_window_half_length );
    mdct_window_sine( st->tcx_cfg.tcx_mdct_window_halfFB, st->tcx_cfg.tcx_mdct_window_half_lengthFB );
    if (st->ini_frame == 0)
    {
        st->tcx_cfg.tcx_last_overlap_mode = st->tcx_cfg.tcx_curr_overlap_mode = ALDO_WINDOW;
    }

    /* TCX Offset */
    st->tcx_cfg.tcx_offset = (st->tcx_cfg.tcx_mdct_window_delay>>1);
    st->tcx_cfg.tcx_offsetFB = (st->tcx_cfg.tcx_mdct_window_delayFB>>1);
    /*<0 rectangular transition with optimized window size = L_frame+L_frame/4*/
    st->tcx_cfg.lfacNext = st->tcx_cfg.tcx_offset - st->L_frame/4;
    st->tcx_cfg.lfacNextFB = st->tcx_cfg.tcx_offsetFB - st->L_frameTCX/4;

    /* set number of coded lines */
    st->tcx_cfg.tcx_coded_lines = getNumTcxCodedLines(bandwidth);

    /* TNS in TCX */
    st->tcx_cfg.pCurrentTnsConfig = NULL;
    st->tcx_cfg.fIsTNSAllowed = getTnsAllowed( st->total_brate, st->igf );

    if (st->tcx_cfg.fIsTNSAllowed)
    {
        InitTnsConfigs( bwMode2fs[bandwidth], st->tcx_cfg.tcx_coded_lines, st->tcx_cfg.tnsConfig, st->hIGFDec.infoIGFStopFreq, st->total_brate );
    }

    resetTecDec( &(st->tecDec) );

    /* Initialize decoder delay */
    /*Constraint for adaptive BPF, otherwise parameter estimation and post-processing not time aligned*/
    if( st->tcxonly == 0 )
    {
        i = st->tcx_cfg.lfacNext>0?st->tcx_cfg.lfacNext:0;
        assert( 0==i );
    }

    st->flag_cna = 0;
    if( st->ini_frame == 0 )
    {
        st->last_flag_cna = 0;
    }

    /* Static vectors to zero */
    if (st->ini_frame == 0)
    {
        st->last_is_cng = 0;

        st->rate_switching_reset = 0;
        set_zero(st->old_syn_Overl, L_FRAME32k/2);

        set_zero(st->syn_Overl_TDAC, L_FRAME32k/2);
        set_zero(st->syn_OverlFB, L_FRAME_MAX/2);
        set_zero(st->syn_Overl_TDACFB, L_FRAME_MAX/2);

        set_zero(st->syn_Overl, L_FRAME32k/2);

        set_zero(st->old_synth, OLD_SYNTH_INTERNAL_DEC);

        set_zero( st->synth_history, L_PROT48k + L_FRAME_MAX);

        set_zero( st->syn, M+1 );

        set_zero( st->mem_syn_r, L_SYN_MEM );

        mem_syn_r_size_old = 0;         /* just to avoid MSVC warnings */
        mem_syn_r_size_new = 0;         /* just to avoid MSVC warnings */

        st->con_tcx = 0;
    }
    else
    {
        /* Reset old_synth in case of core sampling rate switching and Mode 1/2 switching*/
        if( (st->L_frame != st->last_L_frame) || (st->last_codec_mode == MODE1) )
        {
            set_zero(st->old_synth, OLD_SYNTH_INTERNAL_DEC);
        }

        /*Compute size of old and new memories*/
        mem_syn_r_size_old=(int)(1.25*st->last_L_frame/20.f);
        mem_syn_r_size_new=(int)(1.25*st->L_frame/20.f);

        /*Reset LPC mem*/
        if( (st->L_frame != st->last_L_frame) || (st->last_core == AMR_WB_CORE) || (st->last_core == HQ_CORE) )
        {
            set_zero( st->mem_MA, M );
            if( st->sr_core == 16000 )
            {
                mvr2r( GEWB2_Ave, st->mem_AR, M );
            }
            else
            {
                mvr2r( GEWB_Ave, st->mem_AR, M );
            }
        }

        /*Mode 1/2 switching*/
        if( st->last_codec_mode == MODE1 )
        {
            mvr2r( st->lsp_old, st->lspold_uw, M );
            mvr2r( st->lsf_old, st->lsfold_uw, M );
            set_zero( st->syn, M );
        }
        if( st->last_core == AMR_WB_CORE )
        {
            st->last_core = ACELP_CORE;
            st->last_core_bfi = ACELP_CORE;
        }

        if( st->last_codec_mode == MODE1 && st->last_core == ACELP_CORE )
        {
            /* Switching from Mode 1 ACELP */
            st->last_core_bfi = ACELP_CORE;

            /*PLC*/
            if( st->prev_bfi != 0 )
            {
                float *w;
                short W1,nz,delay_comp;

                W1 = st->tcx_cfg.tcx_mdct_window_lengthFB;
                w = st->tcx_cfg.tcx_mdct_windowFB;

                nz = NS2SA(st->output_Fs, N_ZERO_MDCT_NS);
                delay_comp = NS2SA(st->output_Fs, DELAY_CLDFB_NS); /*CLDFB delay*/

                mvr2r(st->fer_samples+delay_comp, st->syn_OverlFB, st->L_frameTCX/2);
                lerp (st->fer_samples+delay_comp, st->syn_Overl  , st->L_frame   /2, st->L_frameTCX/2); /*ACELP(bfi)->TCX(rect)*/
                /*old_out needed for MODE1 routine and syn_Overl_TDAC for MODE2 routine*/
                set_f(st->old_out, 0, nz);
                mvr2r(st->fer_samples+delay_comp, st->old_out+nz,W1);

                for (i = 0; i < W1; i++)
                {
                    st->old_out[i+nz] *= w[W1-1-i]*w[W1-1-i];
                }
                set_f(&st->old_out[W1+nz], 0, nz);

                lerp (st->old_out     , st->old_outLB       , st->L_frame     , st->L_frameTCX);

                mvr2r(st->old_out  +nz, st->syn_Overl_TDACFB, st->L_frameTCX/2);
                nz = NS2SA(st->sr_core, N_ZERO_MDCT_NS);
                mvr2r(st->old_outLB+nz, st->syn_Overl_TDAC  , st->L_frame   /2);
            }
        }

        if ( st->last_codec_mode == MODE2 &&
                st->L_frame != st->last_L_frame &&
                ( ( st->m_frame_type == SID_FRAME && st->last_core > ACELP_CORE ) ||
                  (st->last_core > ACELP_CORE && st->core > ACELP_CORE) || st->prev_bfi ) )
        {
            lerp(st->old_outLB, st->old_outLB, st->L_frame, st->last_L_frame);
        }

        /* Rate switching */
        if( st->last_codec_mode == MODE1 && st->last_core == HQ_CORE )
        {
            /* Switching from MDCT */

            /*Reset of ACELP memories*/
            st->rate_switching_reset=1;
            st->tilt_code = TILT_CODE;
            set_zero(st->old_exc, L_EXC_MEM_DEC);
            set_zero(st->syn, 1+M);
            set_zero(st->mem_syn2, M);

            /*OLA -> zero */
            set_zero(st->old_syn_Overl, L_FRAME32k/2);
            set_zero(st->syn_Overl_TDAC, L_FRAME32k/2);
            set_zero(st->syn_Overl, L_FRAME32k/2);
            set_zero(st->syn_Overl_TDACFB, L_FRAME_MAX/2);
            mvr2r(st->old_out+NS2SA(st->output_Fs, N_ZERO_MDCT_NS), st->syn_OverlFB, st->tcx_cfg.tcx_mdct_window_lengthFB);
            st->tcx_cfg.last_aldo=1;  /*It was previously ALDO*/
            st->tcx_cfg.tcx_curr_overlap_mode = ALDO_WINDOW;
            /*OLA for Mode 2 TCX always reset in Mode switching cases*/
            set_f( st->old_outLB, 0, st->L_frame );

            st->last_core_bfi = TCX_20_CORE;
            st->pfstat.on=0;
            /* reset CLDFB memories */
            cldfb_reset_memory( st->cldfbAna );
            cldfb_reset_memory( st->cldfbBPF );
            cldfb_reset_memory( st->cldfbSyn );
        }
        else if( (st->L_frame!=st->last_L_frame) && (st->L_frame<=L_FRAME16k) && (st->last_L_frame <=L_FRAME16k)) /* Rate switching between 12.8 and 16 kHz*/
        {
            /*Interpolation of ACELP memories*/

            /* convert quantized LSP vector */
            st->rate_switching_reset= lsp_convert_poly( st->lsp_old, st->L_frame, 0 );
            lsp2a_stab( st->lsp_old, st->old_Aq_12_8, M );

            lsp2lsf( st->lsp_old, st->lsf_old, M, st->sr_core );
            mvr2r( st->lsp_old, st->lspold_uw, M );
            mvr2r( st->lsf_old, st->lsfold_uw, M );

            if( !st->last_con_tcx )
            {
                synth_mem_updt2( st->L_frame, st->last_L_frame, st->old_exc, st->mem_syn_r, st->mem_syn2, NULL, DEC );
            }

            /*mem of deemphasis stayed unchanged.*/
        }
        else if( st->L_frame!=st->last_L_frame ) /* Rate switching involving TCX only modes */
        {
            /*Partial reset of ACELP memories*/
            st->rate_switching_reset = 1;

            /*reset partly some memories*/
            st->tilt_code = TILT_CODE;
            if( !st->last_con_tcx )
            {
                set_zero( st->old_exc, L_EXC_MEM_DEC );
            }
            set_zero( st->old_Aq_12_8, M+1 );

            /*Resamp others memories*/
            /*Size of LPC syn memory*/
            lerp( st->mem_syn_r+L_SYN_MEM-mem_syn_r_size_old, st->mem_syn_r+L_SYN_MEM-mem_syn_r_size_new, mem_syn_r_size_new, mem_syn_r_size_old );
            mvr2r( st->mem_syn_r+L_SYN_MEM-M, st->mem_syn2, M );
        }
        /* update of lsf_old only needed in BASOP */
        /* else if( !st->tcxonly && (st->L_frame == L_FRAME16k) && (st->last_total_brate > ACELP_32k) ) */
        /* { */
        /*     lsp2lsf( st->lsp_old, st->lsf_old, M, st->sr_core ); */
        /* } */
    }

    if( st->last_bwidth == NB && st->bwidth != NB && st->ini_frame != 0 )
    {
        st->rate_switching_reset=1;
    }

    st->old_synth_len = 2*st->L_frame;
    st->old_synth_lenFB = 2*st->L_frameTCX;

    /* bass pf reset */
    st->bpf_gain_param = 0;
    set_f( st->pst_old_syn, 0, NBPSF_PIT_MAX );

    /* Formant postfilter */
    if (st->ini_frame == 0)
    {
        /* do nothing */
    }
    else if( st->last_codec_mode == MODE2 )
    {
        if( !st->tcxonly )
        {
            if( st->pfstat.on )
            {
                lerp( st->pfstat.mem_stp+L_SYN_MEM-mem_syn_r_size_old, st->pfstat.mem_stp+L_SYN_MEM-mem_syn_r_size_new, mem_syn_r_size_new, mem_syn_r_size_old );
                lerp( st->pfstat.mem_pf_in+L_SYN_MEM-mem_syn_r_size_old, st->pfstat.mem_pf_in+L_SYN_MEM-mem_syn_r_size_new, mem_syn_r_size_new, mem_syn_r_size_old );
            }
            else
            {
                set_zero( st->pfstat.mem_stp, L_SYN_MEM );
                set_zero( st->pfstat.mem_pf_in, L_SYN_MEM );
                st->pfstat.reset = 1.f;
                st->pfstat.gain_prec = 1.f;
            }
        }
        else if( st->pfstat.on )
        {
            lerp( st->pfstat.mem_stp+L_SYN_MEM-mem_syn_r_size_old, st->pfstat.mem_stp+L_SYN_MEM-mem_syn_r_size_new, mem_syn_r_size_new, mem_syn_r_size_old );
            lerp( st->pfstat.mem_pf_in+L_SYN_MEM-mem_syn_r_size_old, st->pfstat.mem_pf_in+L_SYN_MEM-mem_syn_r_size_new, mem_syn_r_size_new, mem_syn_r_size_old );
        }
    }
    else
    {
        /*codec mode switching*/

        /*reset post-filter except for Narrowband*/
        if ( ((short)(st->output_Fs / 50)) != L_FRAME8k )
        {
            st->pfstat.reset = 1;
            if( st->pfstat.on != 0 )
            {
                st->pfstat.reset = 0;
                lerp( st->pfstat.mem_stp+L_SYN_MEM-mem_syn_r_size_old, st->pfstat.mem_stp+L_SYN_MEM-mem_syn_r_size_new, mem_syn_r_size_new, mem_syn_r_size_old );
                lerp( st->pfstat.mem_pf_in+L_SYN_MEM-mem_syn_r_size_old, st->pfstat.mem_pf_in+L_SYN_MEM-mem_syn_r_size_new, mem_syn_r_size_new, mem_syn_r_size_old );
            }
        }
        else
        {
            /*feed last value old_synth as it is used for pre-emphasis mem*/
            st->old_synth[st->old_synth_len-1]=st->syn[M];
            st->pst_old_syn[NBPSF_PIT_MAX-1]=st->syn[M];
        }
    }

    /* lsf and lsp initialization */
    if (st->ini_frame == 0 )
    {
        mvr2r( st->lsp_old, st->lspold_uw, M );
        mvr2r( st->lsf_old, st->lsfold_uw, M );

        set_zero( st->lsf_cng, M );
    }

    st->seed_tcx_plc = 21845;
    st->past_gpit = 0.0f;
    st->past_gcode = 0.0f;
    st->gc_threshold = 0.0f;

    lsf2lsp(st->lsf_cng, st->lspold_cng, M, INT_FS_12k8);
    lsp2a_stab(st->lspold_cng, st->Aq_cng, M);
    st->plcBackgroundNoiseUpdated = 0;
    mvr2r(st->lsf_old, st->lsf_q_cng, M);
    mvr2r(st->lsf_old, st->old_lsf_q_cng, M);
    mvr2r(st->lsp_old, st->lsp_q_cng, M);
    mvr2r(st->lsp_old, st->old_lsp_q_cng, M);
    set_zero(st->mem_syn_unv_back, M);
    st->last_gain_syn_deemph = 1.f;

    if( st->last_codec_mode == MODE1 || st->ini_frame == 0 )
    {
        /* this assumes that MODE1 fades out in the frequency domain -
           otherwise some data from MODE1 would be needed here */
        st->last_concealed_gain_syn_deemph = 1.f;
        st->conceal_eof_gain = 1.0f;
    }
    /* Post processing */
    set_zero(st->mem_Aq, NB_SUBFR16k*(M+1));

    st->lp_ener_bfi = 60.0f;
    if( st->ini_frame == 0 )
    {
        st->prev_bfi = 0;
        st->last_core_bfi = -1;
        st->tcxConceal_recalc_exc = 0;
    }
    st->prev_old_bfi = 0;

    st->nbLostCmpt = 0;
    st->noise_filling_index = 0;

    mvr2r(st->lsf_old, st->lsf_adaptive_mean, M);
    mvr2r(st->lsf_old, st->lsfoldbfi0, M);
    mvr2r(st->lsf_old, st->lsfoldbfi1, M);

    st->clas_dec  = UNVOICED_CLAS;

    if (!st->last_con_tcx)
    {
        st->old_enr_LP = 0.0f;              /* LP filter E of last good voiced frame or local LP filter E in TD TCX PLC */
    }

    if (st->prev_bfi)
    {
        /* calculate energy at the end of the previous frame */
        if( st->core == ACELP_CORE && st->last_core == HQ_CORE )
        {
            fer_energy( st->L_frameTCX, UNVOICED_CLAS, st->previoussynth, -1, &st->enr_old, 1 );
        }
    }
    else
    {
        st->last_good = UNVOICED_CLAS;      /* last good received frame for concealment */
        st->enr_old = 0.0f;                 /* energy at the end of the previous frame */
    }
    st->lp_gainc = 0.0f;
    st->lp_gainp = 0.0f;

    st->prev_widow_left_rect = 0;

    st->CngLevelBackgroundTrace_bfi  = PLC_MIN_CNG_LEV;
    st->NoiseLevelIndex_bfi          = PLC_MIN_STAT_BUFF_SIZE-1;
    st->CurrLevelIndex_bfi           = 0;
    st->LastFrameLevel_bfi           = PLC_MIN_CNG_LEV;
    set_f( st->NoiseLevelMemory_bfi, PLC_MIN_CNG_LEV, PLC_MIN_STAT_BUFF_SIZE);

    st->cummulative_damping_tcx = 1.0f;
    st->cummulative_damping = 1.0f;

    for( i=0; i<2*NB_SUBFR16k+2; i++ )
    {
        st->old_pitch_buf[i] = (float)st->pit_min;
    }

    for( i=0; i<2*NB_SUBFR16k+2; i++ )
    {
        st->mem_pitch_gain[i] = 1.f;
    }

    st->old_fpitch   = (float)st->pit_min;
    st->old_fpitchFB = (float)st->pit_min_TCX;

    st->rate_switching_init = 1;

    st->reset_mem_AR = 0;

    /* For phase dispersion */
    set_zero(st->dispMem,8);

    st->voice_fac = -1; /* purely unvoiced  */

    /* TCX-LTP */
    st->tcxltp = getTcxLtp(st->sr_core);
    if( st->ini_frame == 0 || st->last_codec_mode == MODE1 )
    {
        st->tcxltp_pitch_int = st->pit_max;
        st->tcxltp_pitch_fr = 0;
        st->tcxltp_last_gain_unmodified = 0.f;
        if (st->ini_frame == 0)
        {
            set_f( st->tcxltp_mem_in, 0.0f, TCXLTP_MAX_DELAY );
            set_f( st->tcxltp_mem_out, 0.0f, L_FRAME48k );
            st->tcxltp_pitch_int_post_prev = 0;
            st->tcxltp_pitch_fr_post_prev = 0;
            st->tcxltp_gain_post_prev = 0.f;
            st->tcxltp_filt_idx_prev = -1;
        }
    }

    st->pst_mem_deemp_err = 0.0f;
    st->tcx_cfg.ctx_hm = getCtxHm(st->total_brate, st->rf_flag);
    st->last_ctx_hm_enabled = 0;

    st->tcx_cfg.resq = getResq(st->total_brate);
    st->tcx_cfg.sq_rounding = 0.375f; /*deadzone of 1.25->rounding=1-1.25/2 (No deadzone=0.5)*/

    st->tcx_lpc_shaped_ari = getTcxLpcShapedAri( st->total_brate, bandwidth, st->rf_flag );

    st->envWeighted = 0;

    if( st->tcxonly )
    {
        st->p_bpf_noise_buf = NULL;
    }
    else
    {
        st->p_bpf_noise_buf = st->bpf_noise_buf;
    }

    if( bandwidth == SWB && (st->total_brate == ACELP_16k40 || st->total_brate == ACELP_24k40) )
    {
        st->tec_tfa = 1;
    }
    else
    {
        st->tec_tfa = 0;
    }

    st->tec_flag = 0;
    st->tfa_flag = 0;

    /* needed in decoder to read the bitstream */
    st->enableGplc = 0;

    st->flagGuidedAcelp = 0;
    st->T0_4th = L_SUBFR;
    st->guidedT0 = st->T0_4th;


    if( st->total_brate >= HQ_48k )
    {
        st->enablePlcWaveadjust = 1;
        if( st->ini_frame == 0 || st->last_total_brate < HQ_48k || st->last_codec_mode == MODE1 || st->force_lpd_reset )
        {
            concealment_init( st->L_frameTCX, &st->plcInfo );
        }
    }
    else
    {
        st->enablePlcWaveadjust = 0;
    }

    /* PLC: [TCX: Tonal Concealment] */
    st->tonalMDCTconceal.nScaleFactors = 0;
    st->tonalMDCTconceal.nSamples = 0;
    st->tonalMDCTconceal.lastPcmOut = 0x0;
    st->tonalMDCTconceal.lastBlockData.tonalConcealmentActive = 0;
    st->tonalMDCTconceal.lastBlockData.nSamples = 0;

    TonalMDCTConceal_Init( &st->tonalMDCTconceal,
                           st->L_frameTCX,
                           st->L_frame, FDNS_NPTS,
                           &st->tcx_cfg
                         );

    st->last_tns_active = 0;
    st->second_last_tns_active = 0;
    st->second_last_core = -1;
    st->tcxltp_second_last_pitch = st->old_fpitch;
    st->tcxltp_third_last_pitch  = st->old_fpitch;


    if( st->total_brate == ACELP_9k60 || st->total_brate == ACELP_16k40 || st->total_brate == ACELP_24k40 )
    {
        st->dec_glr = 1;
    }
    else
    {
        st->dec_glr = 0;
    }

    st->dec_glr_idx = 0;

    st->enableTcxLpc = 1;
    st->VAD = 0;
    st->old_gaintcx_bfi=0.0f;
    st->tcx_cfg.na_scale = 1.f;
    st->tcx_hm_LtpPitchLag = -1;
    st->tcxltp_gain = 0.0f;
    return;
}

