/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <assert.h>
#include "stat_dec.h"
#include "prot.h"
#include "rom_com.h"
#include "options.h"


/*-------------------------------------------------------------------*
* decode_frame_type()
*
*
*--------------------------------------------------------------------*/

static void decode_frame_type
(
    Decoder_State *st
)
{
    int frame_size_index, n, total_brate;

    frame_size_index = 0;
    total_brate = st->total_brate;

    /* Get Frame Type (NULL,SID,ACTIVE) and Frame Mode (2kbps, 4kbps,...) */

    if (st->mdct_sw == MODE1)
    {
        st->m_frame_type = ACTIVE_FRAME;

        for (n=0; n<FRAME_SIZE_NB; n++)
        {
            if (FrameSizeConfig[n].frame_bits == st->total_brate/50)
            {
                frame_size_index = n;
                break;
            }
        }

    }
    else
    {

        /* ZERO Frame */
        if( st->total_brate == FRAME_NO_DATA )
        {
            st->bwidth = st->last_bwidth;
            st->m_frame_type = ZERO_FRAME;
        }

        /* SID frame */
        else if( st->total_brate == SID_2k40 )
        {
            unsigned short frame_len_indicator;
            st->cng_type = get_next_indice(st, 1);
            if( st->cng_type != FD_CNG )
            {
                st->BER_detect = 1;
                st->cng_type = FD_CNG;
            }
            st->m_frame_type = SID_FRAME;
            frame_size_index = 1;
            st->bwidth = get_next_indice(st, 2);

            frame_len_indicator = get_next_indice(st, 1);
            if( st->bwidth == NB )
            {
                if( frame_len_indicator )
                {
                    st->BER_detect = 1;
                }
                frame_len_indicator = 0;
            }
            if( frame_len_indicator == 0 )
            {
                st->L_frame = L_FRAME;
                st->total_brate = 9600;
            }
            else
            {
                st->L_frame = L_FRAME16k;
                if ( st->last_total_brate==16400 || st->last_total_brate==24400 )
                {
                    st->total_brate = st->last_total_brate;
                }
                else
                {
                    st->total_brate = 16400;
                }
            }

            for (n=0; n<FRAME_SIZE_NB; n++)
            {
                if (FrameSizeConfig[n].frame_bits == st->total_brate/50)
                {
                    frame_size_index = n;
                    break;
                }
            }
        }
        /* EVS MODES */
        else
        {
            /* Get Frame mode */
            st->m_frame_type = ACTIVE_FRAME;

            for( n=0; n<FRAME_SIZE_NB; n++ )
            {
                if( FrameSizeConfig[n].frame_bits == st->total_brate/50 )
                {
                    frame_size_index =  n;
                    break;
                }
            }


            if (st->rf_flag == 0)
            {
                /* Get bandwidth mode */
                st->bwidth = get_next_indice(st, FrameSizeConfig[frame_size_index].bandwidth_bits);
                st->bwidth += FrameSizeConfig[frame_size_index].bandwidth_min;
            }
            else
            {
                st->bwidth += FrameSizeConfig[frame_size_index].bandwidth_min;
            }

            if (st->bwidth > FB)
            {
                st->bwidth = FB;
                st->BER_detect = 1;
            }

            if (st->bwidth > SWB && st->total_brate < ACELP_16k40)
            {
                st->bwidth = SWB;
                st->BER_detect = 1;
            }

            /* Get reserved bit */
            if (FrameSizeConfig[frame_size_index].reserved_bits && st->rf_flag == 0)
            {
                int dummy = get_next_indice(st, 1);
                if (dummy != 0)
                {
                    st->BER_detect = 1;
                }
                assert( FrameSizeConfig[frame_size_index].reserved_bits == 1);
            }
        }
    }

    st->rate_switching_init = 0;

    if( st->last_codec_mode != MODE2 || !st->BER_detect )
    {
        /* Mode  or Rate Change */
        if( (st->m_frame_type == ACTIVE_FRAME || st->m_frame_type == SID_FRAME) &&  ( (st->total_brate != st->last_total_brate) || (st->bwidth!=st->last_bwidth ) || (st->last_codec_mode == MODE1) || (st->rf_flag !=st->rf_flag_last) || st->force_lpd_reset) )
        {

            st->rate_switching_init = 1;
            /* Reconf Core */
            mode_switch_decoder_LPD( st, st->bwidth, st->total_brate, frame_size_index );

            /* Reconf. CLDFB: check if the CLDFB works on the right sample rate */
            if ( (st->cldfbAna->no_channels * st->cldfbAna->no_col) != st->L_frame )
            {
                resampleCldfb (st->cldfbAna, (st->L_frame * 50));
                if( st->L_frame <= L_FRAME16k )
                {
                    resampleCldfb (st->cldfbBPF, (st->L_frame * 50));
                }
            }

            if (st->bwidth == NB)
            {
                short nBand_nb = (8000*st->cldfbSyn->no_channels / st->output_Fs);
                st->cldfbSyn->bandsToZero =  st->cldfbSyn->no_channels - nBand_nb;
            }
            else
            {
                st->cldfbSyn->bandsToZero = 0;
            }

            /*Reconf Frequency-domain based CNG*/
            configureFdCngDec( st->hFdCngDec, st->bwidth, st->rf_flag==1&&st->total_brate==13200?9600:st->total_brate, st->L_frame );
            if ( st->last_L_frame!=st->L_frame && st->L_frame<=320 && st->last_L_frame<=320 )
            {
                lerp( st->hFdCngDec->hFdCngCom->olapBufferSynth2, st->hFdCngDec->hFdCngCom->olapBufferSynth2, st->L_frame*2, st->last_L_frame*2 );
                if ( st->m_frame_type==SID_FRAME && st->hFdCngDec->hFdCngCom->frame_type_previous!= ACTIVE_FRAME )
                {
                    lerp( st->hFdCngDec->hFdCngCom->olapBufferSynth, st->hFdCngDec->hFdCngCom->olapBufferSynth, st->L_frame*2, st->last_L_frame*2 );
                    if( st->L_frame==L_FRAME)
                    {
                        for (n=0; n < st->L_frame*2; n++)
                        {
                            st->hFdCngDec->hFdCngCom->olapBufferSynth[n] = st->hFdCngDec->hFdCngCom->olapBufferSynth[n]*1.25f;
                        }
                    }
                    else
                    {
                        for (n=0; n < st->L_frame*2; n++)
                        {
                            st->hFdCngDec->hFdCngCom->olapBufferSynth[n] = st->hFdCngDec->hFdCngCom->olapBufferSynth[n]/1.25f;
                        }
                    }
                }
            }
            if ( st->bwidth!=st->last_bwidth )
            {
                st->hFdCngDec->hFdCngCom->msFrCnt_init_counter = 0;
                st->hFdCngDec->hFdCngCom->init_old = FLT_MAX;
            }

            if( st->tcxonly )
            {
                st->p_bpf_noise_buf = NULL;
            }
            else
            {
                st->p_bpf_noise_buf = st->bpf_noise_buf;
            }
        }
    }

    st->total_brate = total_brate;

    return;

}


/*-------------------------------------------------------------------*
* dec_acelp_tcx_frame()
*
* Main decoding function
*--------------------------------------------------------------------*/

void dec_acelp_tcx_frame(
    Decoder_State *st,                /* i/o: encoder state structure             */
    short *coder_type,        /* o  : coder type                          */
    short *concealWholeFrame, /* i/o: concealment flag                    */
    float *output,            /* o  : synthesis                           */
    float *bpf_noise_buf,     /* i/o: BPF noise buffer                    */
    float *pcmbufFB,
    float bwe_exc_extended[], /* i/o: bandwidth extended excitation       */
    float *voice_factors,     /* o  : voicing factors                     */
    float pitch_buf[]         /* o  : floating pitch for each subframe    */
)
{
    short               i;
    int                 start_bit_pos;
    short               tmp;
    short bitsRead;
    int param[DEC_NPRM_DIV*NB_DIV];

    float old_bwe_exc[(PIT16k_MAX + (L_FRAME16k + 1) + L_SUBFR16k) * 2]; /* excitation buffer */
    float *ptr_bwe_exc;              /* pointer to BWE excitation signal in the current frame */


    start_bit_pos = st->next_bit_pos;
    if( st->rf_flag == 1 )
    {
        start_bit_pos -= 2;
    }

    /* -------------------------------------------------------------- */
    /* IDENTIFY FRAME TYPE                                            */
    /* -------------------------------------------------------------- */

    st->m_old_frame_type = st->m_frame_type;

    if( *concealWholeFrame == 0 )
    {

        unsigned char m_frame_type = st->m_frame_type;
        short bwidth = st->bwidth;
        short cng_type = st->cng_type;
        short L_frame = st->L_frame;
        long total_brate = st->last_total_brate;

        decode_frame_type( st );

        st->force_lpd_reset = 0;

        if( (st->last_codec_mode) == MODE2 && (st->BER_detect || (st->prev_bfi && st->m_frame_type == ZERO_FRAME && st->m_old_frame_type == ACTIVE_FRAME)))
        {
            /* Copy back parameters from previous frame, because there is a high risk they are corrupt
             * Do concealment with configuration used in previous frame                                */
            st->m_frame_type = m_frame_type;
            st->bwidth = bwidth;
            st->cng_type = cng_type;
            st->L_frame = L_frame;
            if( st->ini_frame != 0 )
            {
                st->total_brate = total_brate;
            }

            *concealWholeFrame = 1;
            st->m_decodeMode = DEC_CONCEALMENT_EXT;
            st->BER_detect = 0;

            if( (st->bwidth != st->last_bwidth ) || (st->rf_flag != st->rf_flag_last) || (st->total_brate != st->last_total_brate) )
            {
                st->force_lpd_reset = 1;
            }

            st->core_brate = st->last_core_brate;
            st->bfi = 1;
            if(st->ini_frame == 0 )
            {
                st->tcx_cfg.tcx_coded_lines = getNumTcxCodedLines(SWB);
            }
        }
        else
        {
            st->core_brate = st->total_brate;
            bpf_noise_buf = st->p_bpf_noise_buf;
        }
    }

    if( *concealWholeFrame != 0 )
    {
        /* add two zero bytes for arithmetic coder flush */
        for( i=0; i<8*2; i++ )
        {
            st->bit_stream[i] = 0;
        }
    }


    if( !(st->m_frame_type == SID_FRAME || st->m_frame_type == ZERO_FRAME) )
    {

        /* -------------------------------------------------------------- */
        /* DECODE CORE                                                    */
        /* -------------------------------------------------------------- */

        if( *concealWholeFrame )
        {

            tmp = 0; /*to avoid empty counting */
        }

        tmp = st->total_brate/50 - (st->next_bit_pos - start_bit_pos);

        bitsRead = 0;

        /* update old BWE excitation buffer */
        set_f( old_bwe_exc + PIT16k_MAX * 2, 0.f, ((L_FRAME16k + 1) + L_SUBFR16k) * 2 );
        ptr_bwe_exc = old_bwe_exc + PIT16k_MAX * 2;
        mvr2r( st->old_bwe_exc, old_bwe_exc, PIT16k_MAX * 2 );

        /* Decode the LPD data */
        if( st->m_decodeMode == DEC_NO_FRAM_LOSS )
        {
            decoder_LPD( output, pcmbufFB, &tmp, st, bpf_noise_buf, 0,
                         &bitsRead, coder_type, param, pitch_buf, voice_factors, ptr_bwe_exc );
            if( !st->rate_switching_init && (st->last_codec_mode) == MODE2 && !(st->use_partial_copy && st->rf_frame_type >= RF_TCXFD && st->rf_frame_type <= RF_TCXTD2) && st->bfi )
            {
                *concealWholeFrame = 1;
                st->m_decodeMode = DEC_CONCEALMENT_EXT;
                st->BER_detect = 0;
            }

        }
        else if( st->m_decodeMode == DEC_CONCEALMENT_EXT )
        {
            decoder_LPD( output, pcmbufFB, NULL, st, bpf_noise_buf, 1, /* bfi - st->bfi can be 0 here - MODE2 stays in PLC when DTX appears after a loss */
                         &bitsRead, coder_type, NULL, pitch_buf, voice_factors, ptr_bwe_exc );
        }


        if( ( !st->bfi && (st->prev_bfi || st->prev_use_partial_copy )) || ((st->last_vbr_hw_BWE_disable_dec==1) && (st->vbr_hw_BWE_disable_dec==0)) )
        {
            st->bwe_non_lin_prev_scale = 0.0f;
            set_f( st->old_bwe_exc_extended, 0.0f, NL_BUFF_OFFSET );
        }

        if( st->core == ACELP_CORE && st->igf && st->con_tcx == 0 )
        {
            non_linearity( ptr_bwe_exc, bwe_exc_extended, st->old_bwe_exc_extended, L_FRAME32k, &st->bwe_non_lin_prev_scale, *coder_type, voice_factors, st->L_frame );

            /* update the old BWE exe memory */
            mvr2r( &old_bwe_exc[L_FRAME32k], st->old_bwe_exc, PIT16k_MAX * 2 );
        }
        else
        {
            set_f( st->old_bwe_exc_extended, 0, NL_BUFF_OFFSET );
            set_f( st->old_bwe_exc, 0, PIT16k_MAX * 2 );   /* reset old non_linear exc during igf frames */
            st->bwe_non_lin_prev_scale = 0.0f;
        }

        /* for ACELP mode, skip core data to read TD-BWE side info */
        if( (!st->bfi) && st->core == ACELP_CORE && st->total_brate > 0)
        {
            /* target bs-position "-2", because of earlier "start_bit_pos -= 2;", which are included in "st->rf_target_bits"*/
            /* from "-2" to "-3" as flag-bit not considered in rf_target_bits */
            if (st->rf_flag)
            {
                get_next_indice_tmp(st, start_bit_pos + st->total_brate/50 - st->rf_target_bits - 3 - get_tbe_bits(st->total_brate, st->bwidth, st->rf_flag) - st->next_bit_pos);
            }
            else
            {
                get_next_indice_tmp(st, start_bit_pos + st->total_brate/50 - st->rf_target_bits - get_tbe_bits(st->total_brate, st->bwidth, st->rf_flag) - st->next_bit_pos);
            }
            tbe_read_bitstream(st);
        }

        if( *concealWholeFrame )
        {
            /*"LPD dec - All BFI"*/

            tmp = 0; /*to avoid empty counting */
        }

        /* updates */
        st->last_voice_factor = voice_factors[st->nb_subfr-1];
        st->last_coder_type = *coder_type;
    }
    else
    {


        if( st->m_frame_type == SID_FRAME )
        {
            /* Decode the FD-CNG bitstream */
            FdCng_decodeSID( st );
        }

        /* updates */
        st->last_voice_factor = 0;
        st->last_coder_type = INACTIVE;


    }

    return;
}
