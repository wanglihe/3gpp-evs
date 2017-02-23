/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"

/*-------------------------------------------------------------------*
* Local functions
*-------------------------------------------------------------------*/

static void configure_core_coder( Encoder_State *st, short *coder_type, const short localVAD );

static void writeFrameHeader( Encoder_State *st );

static void initFrameHeader( Encoder_State *st );

/*-------------------------------------------------------------------*
 * evs_enc()
 *
 * Principal encoder routine
 *-------------------------------------------------------------------*/

void evs_enc(
    Encoder_State *st,                                 /* i/o: encoder state structure  */
    const short *data,                               /* i  : input signal             */
    const short n_samples                            /* i  : number of input samples  */
)
{
    short i, input_frame, delay;
    float old_inp_12k8[L_INP_12k8], *inp;            /* buffer of input signal @ 12k8            */
    float old_inp_16k[L_INP];                        /* buffer of input signal @ 16kHz           */
    short sp_aud_decision1;                          /* 1st stage speech/music classification    */
    short sp_aud_decision2;                          /* 2nd stage speech/music classification    */
    float fr_bands[2*NB_BANDS];                      /* energy in frequency bands                */
    short vad_flag;
    short localVAD;
    float Etot;                                      /* total energy; correlation shift          */
    float ener;                                      /* residual energy from Levinson-Durbin     */
    short pitch[3];                                  /* open-loop pitch values for quantization  */
    float voicing[3];                                /* OL maximum normalized correlation        */
    float A[NB_SUBFR16k*(M+1)];                      /* A(z) unquantized for subframes           */
    float Aw[NB_SUBFR16k*(M+1)];                     /* weighted A(z) unquantized for subframes  */
    float epsP[M+1];                                 /* LP prediction errors                     */
    float lsp_new[M];                                /* LSPs at the end of the frame             */
    float lsp_mid[M];                                /* ISPs in the middle of the frame          */
    short coder_type;                                /* coder type                               */
    short sharpFlag;                                 /* formant sharpening flag                  */
    short vad_hover_flag;
    short hq_core_type;                              /* HQ core type (HQ, or LR-MDCT)            */
    short attack_flag;                               /* flag signalling attack encoded by the AC mode (GSC) */
    float new_inp_resamp16k[L_FRAME16k];             /* new input signal @16kHz, non pre-emphasised, used by the WB TBE/BWE */
    float old_syn_12k8_16k[L_FRAME16k];              /* ACELP core synthesis at 12.8kHz or 16kHz to be used by the SWB BWE */
    float shb_speech[L_FRAME16k];
    float hb_speech[L_FRAME16k/4];
    float new_swb_speech[L_FRAME48k];
    float bwe_exc_extended[L_FRAME32k+NL_BUFF_OFFSET];
    float voice_factors[NB_SUBFR16k];
    float fb_exc[L_FRAME16k];
    short Voicing_flag;
    float pitch_buf[NB_SUBFR16k];
    short unbits;
    short padBits;
    float realBuffer[CLDFB_NO_COL_MAX][CLDFB_NO_CHANNELS_MAX]; /* real buffer */
    float imagBuffer[CLDFB_NO_COL_MAX][CLDFB_NO_CHANNELS_MAX]; /* imag buffer */


    /*------------------------------------------------------------------*
     * Initializiation
     *-----------------------------------------------------------------*/

    input_frame = (short)(st->input_Fs / 50);
    st->core = -1;
    st->extl = -1;
    st->core_brate = -1;
    st->input_bwidth = st->last_input_bwidth;
    st->bwidth = st->last_bwidth;
    hq_core_type = -1;
    unbits = 0;

    st->bits_frame_core = 0;              /* For getting bit consumption in core coder */
    st->lp_cng_mode2 = 0;
    st->mdct_sw_enable = 0;
    st->mdct_sw = 0;
    st->rate_switching_reset = 0;

    /*----------------------------------------------------------------*
     * set input samples buffer
     *----------------------------------------------------------------*/

    /* get delay to synchronize ACELP and MDCT frame */
    delay = NS2SA(st->input_Fs, DELAY_FIR_RESAMPL_NS);

    mvr2r( st->input - delay, st->old_input_signal, input_frame+delay );

    /*----------------------------------------------------------------*
     * convert 'short' input data to 'float'
     *----------------------------------------------------------------*/

    for( i=0; i<n_samples; i++ )
    {
        st->input[i] = (float)data[i];
    }

    if( n_samples < input_frame )
    {
        set_f( st->input + n_samples, 0.0f, input_frame - n_samples );
    }

    /*----------------------------------------------------------------*
     * HP filtering
     *----------------------------------------------------------------*/

    hp20( st->input, input_frame, st->mem_hp20_in, st->input_Fs );

    /*----------------------------------------------------------------*
     * Updates in case of AMR-WB IO mode -> EVS primary mode switching
     *----------------------------------------------------------------*/

    if( st->last_core == AMR_WB_CORE )
    {
        updt_IO_switch_enc( st, input_frame );
        set_f( st->old_speech_shb, 0, L_LOOK_16k + L_SUBFR16k );
        cldfb_reset_memory( st->cldfbAnaEnc );
        cldfb_reset_memory( st->cldfbSynTd );
    }

    /*---------------------------------------------------------------------*
     * Pre-processing
     *---------------------------------------------------------------------*/

    pre_proc( st, input_frame, st->input, old_inp_12k8, old_inp_16k, &inp, &sp_aud_decision1,
              &sp_aud_decision2, fr_bands, &vad_flag, &localVAD, &Etot, &ener, pitch, voicing,
              A, Aw, epsP, lsp_new, lsp_mid, &coder_type, &sharpFlag, &vad_hover_flag,
              &attack_flag, new_inp_resamp16k, &Voicing_flag, realBuffer, imagBuffer, &hq_core_type );

    st->sharpFlag = sharpFlag;

    if( st->mdct_sw == MODE2 )
    {

        st->bits_frame_nominal = st->total_brate / 50;
        initFrameHeader( st );

        writeFrameHeader( st );

        if( (st->total_brate > ACELP_24k40 && st->total_brate < HQ_96k) || (st->total_brate == ACELP_24k40 && st->bwidth >= WB) )
        {
            st->L_frame = L_FRAME16k;
            st->gamma = GAMMA16k;
            st->preemph_fac = PREEMPH_FAC_16k;

            weight_a_subfr( NB_SUBFR16k, A, Aw, GAMMA16k, M );

            if( st->last_L_frame == L_FRAME && st->ini_frame != 0 )
            {
                /* this is just an approximation, but it is sufficient */
                mvr2r( st->lsp_old1, st->lspold_enc, M );
            }
        }
        else
        {
            st->L_frame = L_FRAME;
            st->gamma = GAMMA1;
            st->preemph_fac = PREEMPH_FAC;
        }

        st->sr_core = 50*st->L_frame;
        st->core_brate = st->total_brate;

        st->igf = 0;
        hq_core_type = NORMAL_HQ_CORE;

        if( (st->bwidth == SWB || st->bwidth == WB) && st->total_brate <= LRMDCT_CROSSOVER_POINT )
        {
            /* note that FB (bit-rate >= 24400bps) is always coded with NORMAL_HQ_CORE */
            hq_core_type = LOW_RATE_HQ_CORE;
        }
        else if( st->bwidth == NB )
        {
            hq_core_type = LOW_RATE_HQ_CORE;
        }
    }

    /*---------------------------------------------------------------------*
     * Encoding
     *---------------------------------------------------------------------*/

    if( st->codec_mode == MODE1 )
    {
        /* write signalling info into the bitstream */
        signalling_enc( st, coder_type, sharpFlag );

        /*---------------------------------------------------------------------*
         * Preprocessing (preparing) for ACELP/HQ core switching
         *---------------------------------------------------------------------*/

        core_switching_pre_enc( st, &(st->LPDmem), old_inp_12k8, old_inp_16k );

        /*---------------------------------------------------------------------*
         * ACELP core encoding
         *---------------------------------------------------------------------*/

        if( st->core == ACELP_CORE )
        {
            acelp_core_enc( st, &st->LPDmem, inp, vad_flag, ener, pitch, voicing,
                            A, Aw, epsP, lsp_new, lsp_mid, coder_type, sharpFlag, vad_hover_flag, attack_flag,
                            bwe_exc_extended, voice_factors, old_syn_12k8_16k, pitch_buf, &unbits );
        }

        /*---------------------------------------------------------------------*
         * HQ core encoding
         *---------------------------------------------------------------------*/

        if( st->core == HQ_CORE )
        {
            hq_core_enc( st, st->input - delay, input_frame, hq_core_type, Voicing_flag );
        }

        /*---------------------------------------------------------------------*
        * Postprocessing for ACELP/HQ core switching
        *---------------------------------------------------------------------*/

        core_switching_post_enc( st, old_inp_12k8, old_inp_16k, pitch, voicing, A );
    }

    else   /* MODE2 */
    {

        /*----------------------------------------------------------------*
         * Configuration of core coder/SID
         * Write Frame Header
         *----------------------------------------------------------------*/

        configure_core_coder( st, &coder_type, localVAD );

        if( st->mdct_sw != MODE1 )
        {
            writeFrameHeader( st );
        }

        /*----------------------------------------------------------------*
         * Core-Coder
         *----------------------------------------------------------------*/

        /* Call main encoding function */
        enc_acelp_tcx_main( old_inp_16k + L_INP_MEM, st, coder_type, pitch, voicing, Aw, lsp_new,
                            lsp_mid, st->hFdCngEnc, bwe_exc_extended, voice_factors, pitch_buf,
                            vad_hover_flag );

        /*---------------------------------------------------------------------*
         * Postprocessing for Mode 1/2 switching
         *---------------------------------------------------------------------*/
        /* TBE for Mode 2 interface */
        if( st->igf && st->core_brate > SID_2k40 )
        {
            if( st->core == ACELP_CORE )
            {
                switch (st->bwidth)
                {
                case WB:
                    st->extl = WB_TBE;
                    st->extl_brate = WB_TBE_0k35;
                    break;

                case SWB:
                    st->extl = SWB_TBE;
                    st->extl_brate = SWB_TBE_1k6;
                    break;

                case FB:
                    st->extl = FB_TBE;
                    st->extl_brate = FB_TBE_1k8;
                    break;
                }
            }
            else
            {
                coder_type = -1;
                st->extl = IGF_BWE;
                st->extl_brate = 0;
            }

            st->core_brate = st->total_brate - st->extl_brate;

            if( st->tec_tfa == 1 )
            {
                st->core_brate -= BITS_TEC;
                st->core_brate -= BITS_TFA;
            }
        }

        /*----------------------------------------------------------------*
         * Complete Bitstream Writing
         *----------------------------------------------------------------*/

        /* Pad the bitstream with zeros and byte-alignment*/
        if( st->igf && st->core == ACELP_CORE && st->core_brate > SID_2k40 )
        {
            padBits = ((st->bits_frame+7)/8)*8 - (st->nb_bits_tot + (st->rf_target_bits_write - ((st->rf_mode==1)?1:0) ) + get_tbe_bits(st->total_brate, st->bwidth, st->rf_mode ));
        }
        else
        {
            padBits = ((st->bits_frame+7)/8)*8 - (st->nb_bits_tot + (st->rf_target_bits_write - ((st->rf_mode==1)?1:0) ));
        }

        for( i = 0; i<padBits; i++ )
        {
            push_next_indice(st, 0, 1);
        }

    }

    /*---------------------------------------------------------------------*
     * WB TBE encoding
     * WB BWE encoding
     *---------------------------------------------------------------------*/


    if ( st->input_Fs >= 16000 && st->bwidth < SWB )
    {
        /* Common pre-processing for WB TBE and WB BWE */
        wb_pre_proc( st, new_inp_resamp16k, hb_speech );
    }

    if ( st->extl == WB_TBE )
    {
        /* WB TBE encoder */
        wb_tbe_enc( st, coder_type, hb_speech, bwe_exc_extended, voice_factors, pitch_buf, voicing );

        if( st->codec_mode == MODE2 )
        {
            tbe_write_bitstream( st );
        }
    }
    else if ( st->extl == WB_BWE )
    {
        /* WB BWE encoder */
        wb_bwe_enc( st, new_inp_resamp16k, coder_type );
    }

    /*---------------------------------------------------------------------*
     * SWB(FB) TBE encoding
     * SWB(FB) BWE encoding
     *---------------------------------------------------------------------*/

    if( !st->Opt_SC_VBR && st->input_Fs >= 32000 )
    {
        /* Common pre-processing for SWB(FB) TBE and SWB(FB) BWE */
        swb_pre_proc( st, st->input, new_swb_speech, shb_speech, realBuffer, imagBuffer );
    }
    else if( st->input_Fs >= 32000 )
    {
        set_f( st->old_speech_shb, 0.0f, L_LOOK_16k + L_SUBFR16k );
        set_f( shb_speech, 0.0f, L_FRAME16k );
    }

    /* SWB TBE encoder */
    if ( st->extl == SWB_TBE || st->extl == FB_TBE || (st->igf && st->core == ACELP_CORE && st->extl != WB_TBE) )
    {
        if( st->core_brate != FRAME_NO_DATA && st->core_brate != SID_2k40 )
        {
            swb_tbe_enc( st, coder_type, shb_speech, bwe_exc_extended, voice_factors, fb_exc, voicing, pitch_buf );

            if ( st->extl == FB_TBE )
            {
                /* FB TBE encoder */
                fb_tbe_enc( st, st->input, fb_exc );
            }

            if( st->codec_mode == MODE2 )
            {
                if( st->tec_tfa == 1 )
                {
                    tecEnc_TBE( &(st->tecEnc.corrFlag), voicing, coder_type);

                    if( coder_type == INACTIVE )
                    {
                        st->tec_flag = 0;
                        st->tecEnc.corrFlag = 0;
                    }
                    st->tfa_flag = tfaEnc_TBE( st->tfa_enr, st->last_core, voicing, pitch_buf );
                    set_TEC_TFA_code( st->tecEnc.corrFlag, &st->tec_flag, &st->tfa_flag );
                }
                else
                {
                    st->tec_flag = 0;
                    st->tecEnc.corrFlag = 0;
                    st->tfa_flag = 0;
                }

                tbe_write_bitstream( st );
            }
        }
    }
    else if( st->extl == SWB_BWE || st->extl == FB_BWE )
    {
        /* SWB(FB) BWE encoder */
        swb_bwe_enc( st, old_inp_12k8, old_inp_16k, old_syn_12k8_16k, new_swb_speech, shb_speech, coder_type );
    }
    else if( st->extl == SWB_BWE_HIGHRATE || st->extl == FB_BWE_HIGHRATE )
    {
        swb_bwe_enc_hr( st, st->input - delay, input_frame, coder_type, unbits );
    }


    /*---------------------------------------------------------------------*
     * SWB DTX/CNG encoding
     *---------------------------------------------------------------------*/

    if ( st->Opt_DTX_ON && input_frame >= L_FRAME32k )
    {
        /* SHB DTX/CNG encoder */
        swb_CNG_enc( st, shb_speech, old_syn_12k8_16k );
    }


    /*---------------------------------------------------------------------*
     * Channel-aware mode - write signaling information into the bit-stream
     *---------------------------------------------------------------------*/

    signalling_enc_rf( st );


    /*---------------------------------------------------------------------*
     * Updates - MODE1
     *---------------------------------------------------------------------*/


    st->last_sr_core = st->sr_core;
    st->last_codec_mode = st->codec_mode;
    st->last_L_frame = st->L_frame;
    st->last_core = st->core;

    st->last_total_brate = st->total_brate;
    st->last_core_brate = st->core_brate;
    st->last_extl = st->extl;
    st->last_input_bwidth = st->input_bwidth;
    st->last_bwidth = st->bwidth;
    st->Etot_last = Etot;
    st->last_coder_type_raw = st->coder_type_raw;

    if( st->core_brate > SID_2k40 )
    {
        st->last_active_brate = st->total_brate;
    }

    if ( st->core == HQ_CORE )
    {
        /* in the HQ core, coder_type is not used so it could have been set to anything */
        st->prev_coder_type = GENERIC;
    }
    else
    {
        st->prev_coder_type = coder_type;
    }

    if( st->core_brate > SID_2k40 && st->first_CNG == 1 )
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

    if ( st->core_brate <= SID_2k40 && st->first_CNG == 0 && st->cng_type == LP_CNG )
    {
        st->first_CNG = 1;
    }

    /*-----------------------------------------------------------------*
     * Increase the counter of initialization frames
     * Limit the max number of init. frames
     *-----------------------------------------------------------------*/

    if( st->ini_frame < MAX_FRAME_COUNTER )
    {
        (st->ini_frame)++;
    }

    /* synchronisation of CNG seeds */
    if( st->core_brate != FRAME_NO_DATA && st->core_brate != SID_2k40 )
    {
        own_random( &(st->cng_seed) );
        own_random( &(st->cng_ener_seed) );
    }


    /*---------------------------------------------------------------------*
     * Updates - MODE2
     *---------------------------------------------------------------------*/

    if( st->mdct_sw == MODE2 )
    {
        st->codec_mode = MODE2;
        st->sr_core = getCoreSamplerateMode2( st->total_brate, st->bwidth, st->rf_mode );
        st->L_frame = st->sr_core / 50;
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

    /* update FER clas */
    st->last_clas = st->clas;

    /* Update Core */
    core_encode_update( st );

    if( st->mdct_sw == MODE1 )
    {
        st->codec_mode = MODE1;
    }

    if( st->lp_cng_mode2 )
    {
        st->codec_mode = MODE2;
    }

    /* RF mode updates */
    if( st->rf_mode )
    {
        if (st->rf_frame_type == RF_NELP )
        {
            st->last_nelp_mode = 1;
        }
        else
        {
            st->last_nelp_mode = 0;
        }
    }
    st->rf_mode_last = st->rf_mode;
    if( st->Opt_RF_ON )
    {
        st->L_frame = L_FRAME;
        st->rf_mode = 1;
    }




    return;
}


/*-------------------------------------------------------------------*
 * initFrameHeader()
 *
 * Init Mode 2 frame header
 *-------------------------------------------------------------------*/

static void initFrameHeader(
    Encoder_State *st             /* i/o: encoder state structure  */
)
{

    int n;

    if( st->core_brate == SID_2k40 )
    {
        /*Get size of frame*/
        st->bits_frame       = FRAME_2_4;
        st->bits_frame_core += FRAME_2_4-4;     /*1 bit for SID on/off + 2 bits for bandwith in case of SID + 1 bit CNG type */
        st->frame_size_index = 2;
    }
    else if( st->core_brate == FRAME_NO_DATA )
    {
        st->bits_frame       = FRAME_0;
        st->bits_frame_core += st->bits_frame;
        st->frame_size_index = 0;
    }
    else
    {
        for( n=0; n<FRAME_SIZE_NB; n++ )
        {
            if( FrameSizeConfig[n].frame_bits == st->bits_frame_nominal )
            {
                st->frame_size_index = n;
                st->bits_frame = FrameSizeConfig[n].frame_bits;
                st->bits_frame_core = FrameSizeConfig[n].frame_net_bits;
                break;
            }
        }
    }

    return;
}

/*-------------------------------------------------------------------*
 * writeFrameHeader()
 *
 * Write Mode 2 frame header
 *-------------------------------------------------------------------*/

static void writeFrameHeader(
    Encoder_State *st             /* i/o: encoder state structure  */
)
{

    if( st->core_brate != FRAME_NO_DATA )
    {
        /* SID flag at 2.4kbps */
        if( st->core_brate == SID_2k40 )
        {
            if ( st->cng_type == FD_CNG )
            {
                /* write SID/CNG type flag */
                push_next_indice( st, 1, 1 );

                /* write bandwidth mode */
                push_next_indice( st, st->bwidth, 2 );

                /* write L_frame */
                if( st->L_frame == L_FRAME )
                {
                    push_next_indice( st, 0, 1 );
                }
                else
                {
                    push_next_indice( st, 1, 1 );
                }
            }
        }
        else /* active frames */
        {
            if( st->rf_mode == 0 )
            {
                push_next_indice( st, st->bwidth-FrameSizeConfig[st->frame_size_index].bandwidth_min, FrameSizeConfig[st->frame_size_index].bandwidth_bits);
            }
        }

        /* Write reserved bit */
        if( FrameSizeConfig[st->frame_size_index].reserved_bits && st->rf_mode == 0)
        {
            push_next_indice( st, 0, FrameSizeConfig[st->frame_size_index].reserved_bits );
        }
    }

    return;
}


/*------------------------------------------------------------------------*
 * Configuration of core coder/SID
 *------------------------------------------------------------------------*/

static void configure_core_coder(
    Encoder_State *st,            /* i/o: encoder state structure       */
    short *coder_type,    /* i  : coder type                    */
    const short localVAD
)
{
    initFrameHeader( st );

    if( st->core_brate != SID_2k40 && st->core_brate != FRAME_NO_DATA )
    {
        if( st->tcxonly )
        {
            *coder_type = GENERIC;
        }

        st->tcx_cfg.coder_type = *coder_type;


        if( !st->tcxonly && !localVAD && st->tcx_cfg.coder_type == GENERIC )
        {
            st->tcx_cfg.coder_type = UNVOICED;
        }
    }

    st->igf = getIgfPresent( st->total_brate, st->bwidth, st->rf_mode );

    if( st->core_brate != SID_2k40 && st->core_brate != FRAME_NO_DATA )
    {
        st->core_brate = st->total_brate;
    }

    return;
}
