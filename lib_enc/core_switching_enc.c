/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "rom_enc.h"
#include "rom_com.h"
#include "prot.h"
#include "math.h"
#include "assert.h"

/*---------------------------------------------------------------------*
 * core_switching_pre_enc()
 *
 * Preprocessing (preparing) for ACELP/HQ core switching
 *---------------------------------------------------------------------*/

void core_switching_pre_enc(
    Encoder_State *st,            /* i/o: encoder state structure           */
    LPD_state *mem,           /* i/o: encoder state structure           */
    const float *old_inp_12k8,    /* i  : old input signal @12.8kHz         */
    const float *old_inp_16k      /* i  : old input signal @16kHz           */
)
{
    short Sample_Delay_HP, Sample_Delay_LP;


    /* Codec mode switching */
    if( st->last_codec_mode == MODE2 )
    {
        st->mem_deemph = st->LPDmem.syn[M];

        mvr2r( mem->mem_syn2, st->mem_syn1, M );

        st->igf = 0;

        if( st->last_core != ACELP_CORE )
        {
            /* reset BWE memories */
            set_f( st->old_bwe_exc, 0, PIT16k_MAX*2 );
            st->bwe_non_lin_prev_scale = 0.0f;
        }
        set_f( st->old_syn_12k8_16k, 0, NS2SA(16000, DELAY_FD_BWE_ENC_NS) );

        if( st->last_core == TCX_20_CORE || st->last_core == TCX_10_CORE )
        {
            st->last_core = HQ_CORE;


            set_f( st->last_ni_gain, 0, BANDS_MAX );
            set_f( st->last_env, 0, BANDS_MAX );
            st->last_max_pos_pulse = 0;

            st->mode_count = 0;
            st->mode_count1 = 0;

            set_s( st->prev_SWB_peak_pos, 0, SPT_SHORTEN_SBNUM );
            st->prev_frm_hfe2 = 0;
            st->prev_stab_hfe2 = 0;

            /*ALDO overlap windowed past: also used in MODE2 but for other MDCT-LB*/
            set_f( st->old_out, 0, L_FRAME32k );

        }
        if( st->L_frame == L_FRAME16k && st->last_L_frame==L_FRAME )
        {
            mvr2r( st->lsp_old, st->lsp_old16k, M );
            st->rate_switching_reset_16kHz = lsp_convert_poly( st->lsp_old16k, L_FRAME16k, 0 );
        }

        st->use_acelp_preq = 0;
    }

    if( st->last_core == -1 && st->core == HQ_CORE )
    {
        /* very first frame is HQ_CORE */
        st->last_core = HQ_CORE;
    }

    if( st->core == HQ_CORE && (st->last_core == ACELP_CORE || st->last_core == AMR_WB_CORE) ) /* HQ init */
    {

        set_f( st->last_ni_gain, 0, BANDS_MAX );
        set_f( st->last_env, 0, BANDS_MAX );
        st->last_max_pos_pulse = 0;

        st->mode_count = 0;
        st->mode_count1 = 0;

        set_s( st->prev_SWB_peak_pos, 0, SPT_SHORTEN_SBNUM );
        st->prev_frm_hfe2 = 0;
        st->prev_stab_hfe2 = 0;

        set_f( st->old_out, 0, L_FRAME32k );
    }

    /* Here we only handle cases where last_ppp and last_nelp not updated when coming from CodecB or other cores
       within ACELP_CORE if switching from another bitarate to vbr, last_ppp and last_nelp is always updated in the previous frame */
    if( sub(st->core, ACELP_CORE) == 0 && ( st->last_core != ACELP_CORE  || st->last_codec_mode == MODE2))
    {
        st->last_last_ppp_mode = 0;
        st->last_ppp_mode = 0;
        st->last_nelp_mode =0;
    }

    /* Handle state reset of stat_noise_uv_mod memory */
    if( st->core == ACELP_CORE && (st->last_core != ACELP_CORE || st->last_codec_mode == MODE2 || st->last_total_brate <= PPP_NELP_2k80 ) )
    {
        st->act_count = 3;
        st->uv_count = 0;
    }

    if( (st->core == ACELP_CORE || st->core == AMR_WB_CORE) && st->last_core == HQ_CORE )
    {
        /* Reset the ACELP core in case of HQ->ACELP core switching */

        if(st->L_frame == L_FRAME16k )
        {
            mvr2r( TRWB2_Ave, st->lsf_old, M ); /* init of LSP */
            lsf2lsp( st->lsf_old, st->lsp_old, M, INT_FS_16k );
        }
        else
        {
            mvr2r( TRWB_Ave, st->lsf_old, M ); /* init of LSP */
            lsf2lsp( st->lsf_old, st->lsp_old, M, INT_FS_12k8 );
        }

        st->mem_deemph = 0;
        st->LPDmem.syn[M] = 0;
        set_f( mem->mem_syn2, 0.0f, M );
        set_f( mem->mem_syn, 0.0f, M );
        set_f( st->mem_syn1, 0.0f, M );
        st->Nb_ACELP_frames = 0;

        /* Reset ACELP parameters */
        set_zero( st->mem_MA, M );
        mvr2r( GEWB_Ave, st->mem_AR, M );
        mem->mem_w0 = 0.0f;
        mem->tilt_code = 0.0f;
        init_gp_clip(st->clip_var);
        mem->gc_threshold = 0.0f;
        set_f( st->dispMem, 0, 8 );

        st->last_coder_type = GENERIC;

        mvr2r( st->old_pitch_buf + st->L_frame/L_SUBFR, st->old_pitch_buf, st->L_frame/L_SUBFR );
        set_f( st->old_pitch_buf + st->L_frame/L_SUBFR, L_SUBFR, st->L_frame/L_SUBFR );
        /* Reset old ACELP buffers */
        set_f( mem->old_exc, 0, L_EXC_MEM );
        set_f( st->old_bwe_exc, 0, PIT16k_MAX*2 );


        /* reset BWE memories */
        st->bwe_non_lin_prev_scale = 0.0;
        set_f( st->old_syn_12k8_16k, 0, NS2SA(16000, DELAY_FD_BWE_ENC_NS) );
    }

    if( st->input_Fs >= 16000 && st->last_extl != WB_BWE && st->extl == WB_BWE )
    {
        if( st->last_extl != SWB_BWE && st->last_extl != FB_BWE )
        {
            st->prev_mode = NORMAL;
            st->modeCount = 0;
        }

        st->prev_L_swb_norm1 = 8;
    }

    if( ( st->input_Fs >= 32000 && st->last_extl != SWB_BWE && st->extl == SWB_BWE ) ||
            ( st->input_Fs >= 48000 && st->last_extl != FB_BWE && st->extl == FB_BWE ) )
    {
        /* we are switching to SWB BWE - reset SWB BWE buffers */
        if( st->L_frame == L_FRAME )
        {
            Sample_Delay_HP = NS2SA( 16000, ACELP_LOOK_NS + DELAY_FD_BWE_ENC_12k8_NS + DELAY_FIR_RESAMPL_NS - DELAY_CLDFB_NS );
            Sample_Delay_LP = NS2SA( 12800, ACELP_LOOK_NS + DELAY_FD_BWE_ENC_12k8_NS );

            mvr2r( old_inp_12k8 + L_INP_MEM + L_FRAME - Sample_Delay_LP, st->old_input_lp, Sample_Delay_LP );
        }
        else
        {
            Sample_Delay_HP = NS2SA( 16000, ACELP_LOOK_NS + DELAY_FD_BWE_ENC_16k_NS + DELAY_FIR_RESAMPL_NS - DELAY_CLDFB_NS );
            Sample_Delay_LP = NS2SA( 16000, ACELP_LOOK_NS + DELAY_FD_BWE_ENC_16k_NS );

            mvr2r( old_inp_16k + L_INP_MEM + L_FRAME - Sample_Delay_LP, st->old_input_lp, Sample_Delay_LP );
        }
        mvr2r( st->old_speech_shb + L_LOOK_16k + L_SUBFR16k - Sample_Delay_HP, st->new_input_hp, Sample_Delay_HP );

        if (st->last_extl != WB_BWE)
        {
            st->prev_mode = NORMAL;
            st->modeCount = 0;
        }
        st->EnergyLF = 0.0f;
        st->prev_L_swb_norm1 = 8;
    }



    return;
}


/*---------------------------------------------------------------------*
 * core_switching_post_enc()
 *
 * Postprocessing for ACELP/HQ core switching
 *---------------------------------------------------------------------*/

void core_switching_post_enc(
    Encoder_State *st,            /* i/o: encoder state structure             */
    const float *old_inp_12k8,  /* i  : old input signal @12.8kHz           */
    const float *old_inp_16k,   /* i  : old input signal @16kHz             */
    const short pitch[3],       /* i  : open-loop pitch values for quantiz. */
    const float voicing[3],     /* i  : Open-loop pitch gains               */
    const float A[]             /* i  : unquant. LP filter coefs.           */
)
{
    short T_op[3];



    mvs2s( pitch, T_op, 3 );


    if( st->core == HQ_CORE )
    {
        st->use_acelp_preq = 0;

        if( (st->last_core == ACELP_CORE || st->last_core == AMR_WB_CORE) )  /* core switching ==> CELP subframe encoding */
        {
            acelp_core_switch_enc( st, &(st->LPDmem), old_inp_12k8 + L_INP_MEM - NS2SA(INT_FS_12k8, ACELP_LOOK_NS),
                                   old_inp_16k + L_INP_MEM - NS2SA(INT_FS_16k, ACELP_LOOK_NS), T_op, voicing, A );
        }

        st->bwe_non_lin_prev_scale = 0.0;
        st->mem_deemph_old_syn = 0.0f;
    }
    else
    {

        /* reset SWB TBE buffers */
        if( st->extl == WB_TBE && st->last_extl != WB_TBE )
        {
            wb_tbe_extras_reset( st->mem_genSHBexc_filt_down_wb2, st->mem_genSHBexc_filt_down_wb3 );

            if( st->last_extl != WB_BWE )
            {
                set_f( st->decim_state1, 0.0f, (2*ALLPASSSECTIONS_STEEP+1) );
                set_f( st->decim_state2, 0.0f, (2*ALLPASSSECTIONS_STEEP+1) );
            }

            set_f( st->state_syn_shbexc, 0, L_SHB_LAHEAD/4 );
            set_f( st->syn_overlap, 0, L_SHB_LAHEAD );
            set_f( st->mem_csfilt, 0, 2 );
        }

        if( (st->extl == SWB_TBE || st->extl == FB_TBE) &&
                ( st->last_core == HQ_CORE || st->L_frame != st->last_L_frame || (st->last_extl != SWB_TBE && st->last_extl != FB_TBE) )
          )
        {
            set_f( st->state_ana_filt_shb, 0.0f, (2*ALLPASSSECTIONS_STEEP+1) );
            set_f( st->old_speech_shb, 0.0f, L_LOOK_16k + L_SUBFR16k );

            swb_tbe_reset( st->mem_csfilt, st->mem_genSHBexc_filt_down_shb, st->state_lpc_syn, st->syn_overlap,
                           st->state_syn_shbexc, &(st->tbe_demph), &(st->tbe_premph), st->mem_stp_swb, &(st->gain_prec_swb) );

            set_f( st->dec_2_over_3_mem,0.0f, 12 );
            set_f( st->dec_2_over_3_mem_lp,0.0f, 6 );

        }
        else if( (st->extl == SWB_TBE || st->extl == FB_TBE) &&
                 ( (st->last_total_brate != st->total_brate) || (st->last_bwidth != st->bwidth) ||
                   (st->last_codec_mode != MODE1) || (st->rf_mode_last != st->rf_mode) ) )
        {
            set_f( st->state_lpc_syn, 0.0f, LPC_SHB_ORDER );
            set_f( st->state_syn_shbexc, 0.0f, L_SHB_LAHEAD );
            set_f( st->mem_stp_swb, 0.0f, LPC_SHB_ORDER );
            set_f( st->mem_zero_swb, 0, LPC_SHB_ORDER );
            st->gain_prec_swb = 1.0f;
        }

        /* Interp_3_2 CNG buffers reset */

        if( st->extl == FB_TBE && ( st->last_extl != FB_TBE || st->L_frame != st->last_L_frame ) )
        {
            set_f(st->fb_state_lpc_syn, 0, LPC_SHB_ORDER);
            st->fb_tbe_demph = 0;
            fb_tbe_reset_enc( st->elliptic_bpf_2_48k_mem, &st->prev_fb_energy );
        }
    }



    return;
}


/*---------------------------------------------------------------------*
 * core_switching_hq_prepare_enc()
 *
 * Preprocessing in the first HQ frame after ACELP frame
 * - modify bit allocation for HQ core removing CELP subframe budget
 * - update st->old_wtda to modify windows at the encoder
 *---------------------------------------------------------------------*/

void core_switching_hq_prepare_enc(
    Encoder_State *st,                /* i/o: encoder state structure */
    short *num_bits,          /* i/o: bit budget update       */
    const short input_frame,        /* i  : frame length            */
    float *wtda_audio,
    const float *audio
)
{
    short delta,  Loverlapp, i;
    short n;
    long  cbrate;

    /* set multiplication factor according to the sampling rate */
    delta = 1;
    if( input_frame == L_FRAME16k )
    {
        delta = 2;
    }
    else if( input_frame == L_FRAME32k )
    {
        delta = 4;
    }
    else if( input_frame == L_FRAME48k )
    {
        delta = 6;
    }

    /* set switching frame bit-rate */
    if( st->last_L_frame == L_FRAME)
    {
        if( st->core_brate > ACELP_24k40 )
        {
            cbrate = ACELP_24k40;
        }
        else
        {
            cbrate = st->core_brate;
        }
        /* subtract ACELP switching frame bits */
        if( st->core_brate >= ACELP_11k60 )
        {
            /* subtract one bit for LP filtering flag */
            (*num_bits)--;
        }

        *num_bits -= ACB_bits_tbl[BIT_ALLOC_IDX(cbrate, GENERIC, 0, 0)];     /* pitch bits */
        *num_bits -= gain_bits_tbl[BIT_ALLOC_IDX(cbrate, TRANSITION, 0, 0)]; /* gain bits */
        *num_bits -= FCB_bits_tbl[BIT_ALLOC_IDX(cbrate, GENERIC, 0, 0)];     /* FCB bits  */
    }
    else  /* L_frame == L_FRAME16k */
    {
        if( st->core_brate <= ACELP_8k00 )
        {
            cbrate = ACELP_8k00;
        }
        else if( st->core_brate <= ACELP_14k80 )
        {
            cbrate = ACELP_14k80;
        }
        else
        {
            cbrate = min(st->core_brate,ACELP_22k60);
        }

        /* subtract ACELP switching frame bits */
        if( st->core_brate >= ACELP_11k60 )
        {
            /* subtract one bit for LP filtering flag */
            (*num_bits)--;
        }


        *num_bits -= ACB_bits_16kHz_tbl[BIT_ALLOC_IDX_16KHZ(cbrate, GENERIC, 0, 0)];     /* pitch bits */
        *num_bits -= gain_bits_16kHz_tbl[BIT_ALLOC_IDX_16KHZ(cbrate, GENERIC, 0, 0)];    /* gain bits */
        *num_bits -= FCB_bits_16kHz_tbl[BIT_ALLOC_IDX_16KHZ(cbrate, GENERIC, 0, 0)];     /* FCB bits  */
    }

    /* subtract BWE bits */
    if ( !((inner_frame_tbl[st->bwidth]==L_FRAME16k && st->last_L_frame==L_FRAME16k)|| inner_frame_tbl[st->bwidth] == L_FRAME8k) )
    {
        *num_bits -= (NOOFGAINBITS1 + AUDIODELAYBITS);
    }


    n = ((float)input_frame * N_ZERO_MDCT_NS/FRAME_SIZE_NS);

    /* Transition window at the encoder */
    Loverlapp = delta*SWITCH_OVERLAP_8k;
    for( i = 0; i <n; i++ )
    {
        wtda_audio[i+input_frame/2] = - audio[n-i-1];
    }

    for( i = n; i <input_frame/2-Loverlapp; i++ )
    {
        wtda_audio[i+input_frame/2] = - audio[n-i-1];

    }


    for( i = input_frame/2-Loverlapp; i <input_frame/2; i++ )
    {
        wtda_audio[i+input_frame/2] =  - audio[n-i-1] *(float)cos((i+1-input_frame/2+Loverlapp)*EVS_PI/(2*(Loverlapp+1)));  /* win=cos() */
    }

    /* reset state of old_out if switching */
    set_f( st->old_out, 0.0f, L_FRAME32k );

    return;
}

