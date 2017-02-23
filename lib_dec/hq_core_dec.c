/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"

/*--------------------------------------------------------------------------
 *  hq_core_dec()
 *
 *  HQ core decoder
 *--------------------------------------------------------------------------*/

void hq_core_dec(
    Decoder_State *st,                     /* i/o: decoder state structure            */
    float synth[],                 /* o  : output synthesis                   */
    const short output_frame,            /* i  : output frame length                */
    const short hq_core_type,            /* i  : HQ core type                       */
    const short core_switching_flag      /* i  : ACELP->HQ switching frame flag     */
)
{
    short num_bits, is_transient, hqswb_clas, inner_frame;
    short i, j, flag_uv, num_Sb, nb_sfm;
    short ynrm[NB_SFM], num_bands_p[MAX_SB_NB];
    float ener_match, mean_en_high;
    float t_audio_q[L_FRAME48k], wtda_audio[2*L_FRAME48k];
    short delay_comp;
    float normq[NB_SFM], SWB_fenv[SWB_FENV+DIM_FB];
    const short *sfmsize, *sfm_start, *sfm_end;
    float gapsynth[L_FRAME48k];


    /*--------------------------------------------------------------------------
     * Initializations
     *--------------------------------------------------------------------------*/

    set_f( t_audio_q, 0.0f, L_FRAME48k );
    set_f( gapsynth, 0.0f, L_FRAME48k );
    set_s( num_bands_p, 0, MAX_SB_NB );
    set_s( ynrm, 39, NB_SFM );              /* Initialize to the smallest value */
    mean_en_high = 0.0f;
    sfm_start = sfm_end = NULL;
    num_Sb = nb_sfm = 0;

    if (st->tcx_cfg.tcx_curr_overlap_mode == FULL_OVERLAP)
    {
        st->tcx_cfg.tcx_last_overlap_mode = ALDO_WINDOW;
    }
    else
    {
        st->tcx_cfg.tcx_last_overlap_mode = st->tcx_cfg.tcx_curr_overlap_mode;
    }
    st->tcx_cfg.tcx_curr_overlap_mode = ALDO_WINDOW;

    /*--------------------------------------------------------------------------
     * Find the number of bits for transform-domain coding
     *--------------------------------------------------------------------------*/

    /* set the total bit-budget */
    num_bits = (short)(st->total_brate / 50);

    if( !st->bfi )
    {
        if( core_switching_flag )
        {
            /* Preprocessing in the first HQ frame after ACELP frame */
            core_switching_hq_prepare_dec( st, &num_bits, output_frame );

            /* During ACELP->HQ core switching, limit the HQ core bitrate to 48kbps */
            if( num_bits > HQ_48k / 50 )
            {
                num_bits = (short)(HQ_48k / 50);
            }
        }

        /* subtract signalling bits */
        num_bits -= st->next_bit_pos;

        /* set FEC parameters */
        flag_uv = 1 - st->HqVoicing;

        /* subtract the number of bits for pitch & gain at higher bitrates */
        if ( !(core_switching_flag) && st->core_brate > MINIMUM_RATE_TO_ENCODE_VOICING_FLAG )
        {
            st->HqVoicing = (short) get_next_indice( st, 1 );
            num_bits -= 1;
        }
        else
        {
            st->HqVoicing = 0;
            if ( st->core_brate > MINIMUM_RATE_TO_ENCODE_VOICING_FLAG )
            {
                st->HqVoicing = 1;
            }
        }
    }
    else
    {
        flag_uv = 0;
    }

    /* set inner frame (== coded bandwidth) length */
    inner_frame = inner_frame_tbl[st->bwidth];

    if( st->bfi == 0)
    {
        if( output_frame >= L_FRAME16k )
        {
            st->ph_ecu_HqVoicing = st->HqVoicing;
        }
        else
        {
            st->ph_ecu_HqVoicing = 0;
        }
    }

    if( output_frame == L_FRAME8k )
    {
        hq_configure_bfi( &nb_sfm, &num_Sb, num_bands_p, &sfmsize, &sfm_start, &sfm_end );
    }

    /*--------------------------------------------------------------------------
     * transform-domain decoding
     *--------------------------------------------------------------------------*/

    if( st->bfi )
    {
        is_transient = st->old_is_transient[0];
        if( output_frame >= L_FRAME16k )
        {
            hq_ecu( st->prev_good_synth, t_audio_q, &st->time_offs, st->X_sav, &st->num_p, st->plocs, st->plocsi, st->env_stab,
                    &st->last_fec, st->ph_ecu_HqVoicing, &st->ph_ecu_active, gapsynth, st->prev_bfi, st->old_is_transient,
                    st->mag_chg_1st, st->Xavg, &st->beta_mute, output_frame, st );
        }
        else
        {
            HQ_FEC_processing( st, t_audio_q, is_transient, st->ynrm_values, st->r_p_values, num_Sb, nb_sfm, num_bands_p,
                               output_frame, sfm_start, sfm_end );
        }

        st->old_is_transient[2] = st->old_is_transient[1];
        st->old_is_transient[1] = st->old_is_transient[0];

        if( output_frame >= L_FRAME16k )
        {
            /* keep st->previoussynth updated as in FEC_HQ_pitch_analysis but no LP analysis */
            delay_comp = NS2SA(st->output_Fs, DELAY_CLDFB_NS);
            mvr2r( st->previoussynth + delay_comp, st->previoussynth, output_frame - delay_comp );
            mvr2r( st->delay_buf_out, st->previoussynth + output_frame - delay_comp, delay_comp );

            flag_uv = 1;                                  /*  disable costly pitch out synthesis in bfi frame  */
            st->HqVoicing = 1-flag_uv;                    /*  safety setting for HQ->ACELP switch logic           */
            set_f( st->fer_samples, 0.0f, L_FRAME48k );   /*  safety, create a known signal state for HQ->ACELP switch logic */
        }
    }
    else
    {
        if( hq_core_type == LOW_RATE_HQ_CORE )
        {
            if(st->prev_bfi == 1 )
            {
                set_f( st->last_ni_gain, 0, BANDS_MAX );
                set_f( st->last_env, 0, BANDS_MAX );
                st->last_max_pos_pulse = 0;
            }

            /* HQ low rate decoder */
            hq_lr_dec( st, t_audio_q, inner_frame, num_bits, &is_transient );

            hqswb_clas = is_transient;
        }
        else
        {
            /* HQ high rate decoder */
            hq_hr_dec( st, t_audio_q, inner_frame, num_bits, ynrm, &is_transient, &hqswb_clas, SWB_fenv );
        }

        /* scaling (coefficients are in nominal level) */
        if( output_frame != NORM_MDCT_FACTOR )
        {
            ener_match = (float)sqrt((float)output_frame/(float)NORM_MDCT_FACTOR);

            for( i=0; i<inner_frame; i++ )
            {
                t_audio_q[i] *= ener_match;
            }
        }

        HQ_FEC_Mem_update( st, t_audio_q, normq, ynrm, num_bands_p, is_transient, hqswb_clas,
                           core_switching_flag, nb_sfm, num_Sb, &mean_en_high,
                           hq_core_type, output_frame );
    }

    /*--------------------------------------------------------------------------
     * Attenuate HFs in case of band-width switching (from higher BW to lower BW)
     *--------------------------------------------------------------------------*/

    /* attenuate HFs in case of band-width switching */
    if( st->bws_cnt1 > 0 )
    {
        ener_match = (float)st->bws_cnt1 / (float)N_NS2W_FRAMES;

        if( is_transient )
        {
            for( i = 0; i < NUM_TIME_SWITCHING_BLOCKS; i++ )
            {
                for( j=inner_frame_tbl[st->bwidth-1]/NUM_TIME_SWITCHING_BLOCKS; j<inner_frame/NUM_TIME_SWITCHING_BLOCKS; j++ )
                {
                    t_audio_q[i*inner_frame/NUM_TIME_SWITCHING_BLOCKS + j] *= ener_match;
                }
            }
        }
        else
        {
            for( i=inner_frame_tbl[st->bwidth-1]; i<inner_frame; i++ )
            {
                t_audio_q[i] *= ener_match;
            }
        }
    }

    /* WB/SWB bandwidth switching */
    if( is_transient )
    {
        mvr2r( t_audio_q + 240, st->t_audio_q, 80 );
    }
    else
    {
        mvr2r( t_audio_q, st->t_audio_q, L_FRAME );
    }

    /*--------------------------------------------------------------------------
     * Inverse transform
     * Overlap-add
     * Pre-echo reduction
     *--------------------------------------------------------------------------*/

    if( output_frame == L_FRAME8k || st->bfi == 0)
    {
        if((output_frame != inner_frame) && (st->bfi == 1))
        {
            inverse_transform( t_audio_q, wtda_audio, is_transient, output_frame, output_frame );
        }
        else
        {
            inverse_transform( t_audio_q, wtda_audio, is_transient, output_frame, inner_frame );
        }
    }
    if( output_frame == L_FRAME8k )
    {
        if( st->bfi == 0 && st->prev_bfi == 0 )
        {
            mvr2r( st->old_out + (short)(N_ZERO_MDCT_NS*output_frame/FRAME_SIZE_NS), st->prev_oldauOut, output_frame - (short)(N_ZERO_MDCT_NS*output_frame/FRAME_SIZE_NS) );
        }
        else if( st->prev_bfi == 1 )
        {
            set_f( st->prev_oldauOut, 0.0f, output_frame );
        }
        if( (st->prev_bfi == 1 || st->bfi == 1) && !st->old_is_transient[2] && st->last_core == HQ_CORE && st->last_codec_mode == MODE1)
        {
            time_domain_FEC_HQ( st, wtda_audio, synth, mean_en_high, output_frame );
        }
        else
        {
            window_ola( wtda_audio, synth, st->old_out, output_frame, st->tcx_cfg.tcx_last_overlap_mode, st->tcx_cfg.tcx_curr_overlap_mode, st->prev_bfi, st->oldHqVoicing, st->oldgapsynth );
            st->phase_mat_next = 0;
        }

        if( (!st->bfi && !st->prev_bfi) || (!(output_frame >= L_FRAME16k)) )
        {
            preecho_sb( st->core_brate, wtda_audio, synth, output_frame, &st->memfilt_lb, &st->mean_prev_hb, &st->smoothmem,
                        &st->mean_prev, &st->mean_prev_nc, &st->wmold_hb, &st->prevflag, &st->pastpre, st->bwidth );
        }

    }
    else
    {
        if( st->bfi && output_frame >= L_FRAME16k)
        {
            window_ola( t_audio_q, synth, st->old_out, output_frame, ALDO_WINDOW, ALDO_WINDOW, st->prev_bfi && !st->ph_ecu_active, st->oldHqVoicing, st->oldgapsynth );
        }
        else
        {
            window_ola( wtda_audio, synth, st->old_out, output_frame, st->tcx_cfg.tcx_last_overlap_mode, st->tcx_cfg.tcx_curr_overlap_mode, st->prev_bfi && !st->ph_ecu_active, st->oldHqVoicing, st->oldgapsynth );
        }

        if( (!st->bfi && !st->prev_bfi) || (!(output_frame >= L_FRAME16k)) )
        {
            preecho_sb( st->core_brate, wtda_audio, synth, output_frame, &st->memfilt_lb, &st->mean_prev_hb, &st->smoothmem,
                        &st->mean_prev, &st->mean_prev_nc, &st->wmold_hb, &st->prevflag, &st->pastpre, st->bwidth );
        }

    }

    if (!st->bfi
            && st->prev_bfi
            && st->last_total_brate >= HQ_48k
            && st->last_codec_mode == MODE2
            && (st->last_core_bfi == TCX_20_CORE || st->last_core_bfi == TCX_10_CORE)
            && st->plcInfo.concealment_method == TCX_NONTONAL
            && st->plcInfo.nbLostCmpt < 4 )
    {
        waveform_adj2(
            st->tonalMDCTconceal.secondLastPcmOut,
            synth,
            st->plcInfo.data_noise,
            &st->plcInfo.outx_new_n1,
            &st->plcInfo.nsapp_gain,
            &st->plcInfo.nsapp_gain_n,
            &st->plcInfo.recovery_gain,
            st->plcInfo.step_concealgain,
            st->plcInfo.Pitch,
            st->plcInfo.FrameSize,
            0,
            st->plcInfo.nbLostCmpt + 1,
            st->bfi);
    }

    if ( output_frame >= L_FRAME16k )
    {
        if( st->ph_ecu_HqVoicing )
        {
            st->oldHqVoicing = 1;
            mvr2r( gapsynth, st->oldgapsynth, L_FRAME48k );
        }
        else
        {
            st->oldHqVoicing = 0;
        }
    }
    else
    {
        st->oldHqVoicing = 0;
    }

    if( st->nbLostCmpt == FRAMECTTOSTART_MDCT )
    {
        st->HqVoicing = 0;
    }

    if( output_frame == L_FRAME8k )
    {
        mvr2r( wtda_audio, st->oldIMDCTout, L_FRAME8k/2 );
        mvr2r( st->old_auOut_2fr + output_frame, st->old_auOut_2fr, output_frame );
        mvr2r( synth, st->old_auOut_2fr + output_frame, output_frame );
    }


    /* update buffer of old subframe pitch values */
    if( st->last_core == HQ_CORE && st->L_frame != st->last_L_frame )
    {
        set_f( &st->old_pitch_buf[st->L_frame/L_SUBFR], (float)L_SUBFR, st->L_frame/L_SUBFR );
    }
    mvr2r( &st->old_pitch_buf[st->L_frame/L_SUBFR], st->old_pitch_buf, st->L_frame/L_SUBFR );
    set_f( &st->old_pitch_buf[st->L_frame/L_SUBFR], (float)L_SUBFR, st->L_frame/L_SUBFR );
    mvr2r( &st->mem_pitch_gain[2], &st->mem_pitch_gain[st->L_frame/L_SUBFR+2], st->L_frame/L_SUBFR );
    set_zero( &st->mem_pitch_gain[2], st->L_frame/L_SUBFR );


    return;
}

