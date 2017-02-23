/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"

/*-------------------------------------------------------------------*
 * decod_audio()
 *
 * Decode audio (AC) frames
 *-------------------------------------------------------------------*/

void decod_audio(
    Decoder_State *st,                  /* i/o: decoder static memory                     */
    float dct_epit[],           /* o  : GSC excitation in DCT domain              */
    const float *Aq,                  /* i  : LP filter coefficient                     */
    const short coder_type,           /* i  : coding type                               */
    float *tmp_noise,           /* o  : long term temporary noise energy          */
    float *pitch_buf,           /* o  : floating pitch values for each subframe   */
    float *voice_factors,       /* o  : voicing factors                           */
    float *exc,                 /* i/o: adapt. excitation exc                     */
    float *exc2,                /* i/o: adapt. excitation/total exc               */
    float *bwe_exc,             /* o  : excitation for SWB TBE                    */
    float *lsf_new,             /* i  : ISFs at the end of the frame              */
    float *gain_buf
)
{
    short tmp_nb_bits_tot, pit_band_idx;
    float code[L_SUBFR];
    short Diff_len, nb_subfr, i;
    short nb_frame_flg;
    float Es_pred = 0.0f;
    short Len, max_len;
    short attack_flag;
    float low_pit;
    short last_bin;
    short nbits;
    float exc_wo_nf[L_FRAME];

    short nb_bits; /* number of bits */
    int   indice;   /* parameter indices to read */

    /*---------------------------------------------------------------*
     * Initialization
     *---------------------------------------------------------------*/

    Diff_len  = 0;

    /* decode GSC attack flag (used to reduce possible pre-echo) */
    attack_flag = (short) get_next_indice( st, 1 );

    /* decode GSC SWB speech flag */
    if( coder_type != INACTIVE && st->total_brate >= ACELP_13k20 )
    {
        st->GSC_noisy_speech = (short) get_next_indice( st, 1 );
    }

    /* safety check in case of bit errors */
    if( st->GSC_noisy_speech && st->bwidth != SWB )
    {
        st->BER_detect = 1;
        st->GSC_noisy_speech = 0;
    }

    /*---------------------------------------------------------------*
     * Decode energy dynamics
     *---------------------------------------------------------------*/

    if( st->GSC_noisy_speech )
    {
        nb_subfr = NB_SUBFR;
        st->cor_strong_limit = 0;
        st->noise_lev = NOISE_LEVEL_SP3;
    }
    else
    {
        if( st->core_brate <= ACELP_8k00 )
        {
            st->noise_lev = (short)get_next_indice( st, 2 ) + NOISE_LEVEL_SP2;
        }
        else
        {
            st->noise_lev = (short)get_next_indice( st, 3 ) +  NOISE_LEVEL_SP0;
        }

        /*---------------------------------------------------------------*
         * Decode number of subframes
         *---------------------------------------------------------------*/

        st->cor_strong_limit = 1;
        nb_subfr = SWNB_SUBFR;
        if( st->core_brate >= ACELP_9k60 )
        {
            nbits = 1;

            nb_frame_flg = (short)get_next_indice( st, nbits );

            if( (nb_frame_flg & 0x1) == 0)
            {
                nb_subfr = 2*SWNB_SUBFR;
                st->cor_strong_limit = 0;
            }
        }
    }

    /*---------------------------------------------------------------*
     * Decode the last band where the adaptive (pitch) contribution is significant
     *---------------------------------------------------------------*/

    if( st->core_brate < CFREQ_BITRATE )
    {
        if( st->core_brate < ACELP_9k60 && coder_type == INACTIVE )
        {
            nbits = 1;
        }
        else
        {
            nbits = 3;
        }
    }
    else
    {
        nbits = 4;
    }

    if( st->core_brate < ACELP_9k60 && coder_type != INACTIVE )
    {
        pit_band_idx = 1;
    }
    else
    {
        pit_band_idx = (short)get_next_indice( st, nbits );
    }

    if( pit_band_idx != 0 )
    {
        if( st->core_brate < ACELP_9k60 )
        {
            pit_band_idx = 7+BAND1k2;  /* At low rate, if pitch model is chosen, then for to be use on extented and constant frequency range */
        }
        else
        {
            pit_band_idx += BAND1k2;
        }

        /* detect bit errors in the bitstream */
        if( pit_band_idx > 13 ) /* The maximum decodable index is 10 + BAND1k2 (3) = 13 */
        {
            pit_band_idx = 13;
            st->BER_detect = 1;
        }

        Diff_len = (short)(mfreq_loc[pit_band_idx]/BIN_SIZE);
    }

    st->Last_GSC_pit_band_idx = pit_band_idx;

    /*--------------------------------------------------------------------------------------*
     * Decode adaptive (pitch) excitation contribution
     * Reset unvaluable part of the adaptive (pitch) excitation contribution
     *--------------------------------------------------------------------------------------*/

    if( pit_band_idx > BAND1k2 )
    {
        /*---------------------------------------------------------------*
         * Decode adaptive (pitch) excitation contribution
         *---------------------------------------------------------------*/

        if( st->GSC_noisy_speech && nb_subfr == NB_SUBFR )
        {
            nb_bits = Es_pred_bits_tbl[BIT_ALLOC_IDX(st->core_brate, GENERIC, -1, -1)];
            indice = (short)get_next_indice( st, nb_bits );
            Es_pred_dec(&Es_pred, indice, nb_bits, 0);
        }

        dec_pit_exc( st, L_FRAME, Aq, Es_pred, pitch_buf, code, exc, nb_subfr, gain_buf );

        if( st->core_brate < ACELP_9k60 )
        {
            minimum( pitch_buf, L_FRAME>>6, &low_pit);

            if( low_pit < 64 )
            {
                pit_band_idx = 9+BAND1k2;
                if(st->bwidth == NB)
                {
                    pit_band_idx = 7+BAND1k2;
                }
            }
            else if( low_pit < 128 )
            {
                pit_band_idx = 5+BAND1k2;
            }
            else
            {
                pit_band_idx = 3+BAND1k2;
            }

            Diff_len = (short)(mfreq_loc[pit_band_idx]/BIN_SIZE);
            st->Last_GSC_pit_band_idx = pit_band_idx;
        }

        /*---------------------------------------------------------------*
         * DCT transform
         *---------------------------------------------------------------*/

        edct( exc, dct_epit, L_FRAME );

        /*---------------------------------------------------------------*
         * Reset unvaluable part of the adaptive (pitch) excitation contribution
         *---------------------------------------------------------------*/

        max_len = L_FRAME - Diff_len;
        if(st->bwidth == NB)
        {
            max_len = 160-Diff_len;
        }

        Len = 80;
        if( max_len < 80 )
        {
            Len = max_len;
        }

        if(st->core_brate == ACELP_8k00 && st->bwidth != NB )
        {
            for (i=0; i < max_len; i++)
            {
                dct_epit[i+Diff_len] = 0.0f;
            }
        }
        else
        {
            for (i = 0; i < Len; i++)
            {
                dct_epit[i+Diff_len] *= sm_table[i];
            }
            for (; i < max_len; i++)
            {
                dct_epit[i+Diff_len] = 0.0f;
            }
        }

        st->bfi_pitch = (short)(mean(pitch_buf, nb_subfr)+0.5f);
        st->bfi_pitch_frame = L_FRAME;

        Diff_len++;
        st->bpf_off = 0;
    }
    else
    {
        /* No adaptive (pitch) excitation contribution */
        st->bpf_off = 1;
        set_f( dct_epit, 0.0f, L_FRAME );
        set_f( pitch_buf, (float)L_SUBFR, NB_SUBFR );

        set_f( gain_buf, 0.f, NB_SUBFR16k);

        st->bfi_pitch = L_SUBFR;
        st->bfi_pitch_frame = L_FRAME;
        st->lp_gainp = 0.0f;
        st->lp_gainc = 0.0f;
        st->tilt_code = 0;
        pit_band_idx = 0;
        Diff_len = 0;
    }

    /*--------------------------------------------------------------------------------------*
     * GSC decoder
     *--------------------------------------------------------------------------------------*/

    /* find the current total number of bits used */
    tmp_nb_bits_tot = st->next_bit_pos;

    if( st->extl_brate > 0 )
    {
        /* subtract 1 bit for TBE/BWE BWE flag (bit counted in extl_brate) */
        tmp_nb_bits_tot--;
    }

    if( coder_type == INACTIVE && st->core_brate <= ACELP_9k60 )
    {
        tmp_nb_bits_tot += 5;  /* for noisiness */
    }

    gsc_dec( st, dct_epit, pit_band_idx, Diff_len, tmp_nb_bits_tot, nb_subfr, coder_type, &last_bin, lsf_new, exc_wo_nf, tmp_noise );

    /*--------------------------------------------------------------------------------------*
     * iDCT transform
     *--------------------------------------------------------------------------------------*/

    edct( dct_epit, exc, L_FRAME );
    edct( exc_wo_nf, exc_wo_nf, L_FRAME );

    /*----------------------------------------------------------------------*
     * Remove potential pre-echo in case an onset has been detected
     *----------------------------------------------------------------------*/

    pre_echo_att( &st->Last_frame_ener, exc, attack_flag, st->last_coder_type );

    /*--------------------------------------------------------------------------------------*
     * Update BWE excitation
     *--------------------------------------------------------------------------------------*/

    set_f( voice_factors, 0.0f, NB_SUBFR16k );
    interp_code_5over2( exc, bwe_exc, L_FRAME );

    /*--------------------------------------------------------------------------------------*
     * Updates
     *--------------------------------------------------------------------------------------*/

    mvr2r( exc, exc2, L_FRAME );
    mvr2r( exc_wo_nf, exc, L_FRAME );

    /*--------------------------------------------------------------------------------------*
     * Channel aware mode parameters
     *--------------------------------------------------------------------------------------*/

    set_f( st->tilt_code_dec, 0, NB_SUBFR16k );

    return;
}



/*-------------------------------------------------------------------*
 * gsc_dec()
 *
 * Generic audio signal decoder
 *-------------------------------------------------------------------*/

void gsc_dec(
    Decoder_State *st,              /* i/o: State structure                                       */
    float exc_dct_in[],     /* i/o: dct of pitch-only excitation / total excitation       */
    const short pit_band_idx,     /* i  : bin position of the cut-off frequency                 */
    const short Diff_len,         /* i  : Lenght of the difference signal (before pure spectral)*/
    const short bits_used,        /* i  : Number of bit used before frequency Q                 */
    const short nb_subfr,         /* i  : Number of subframe considered                         */
    const short coder_type,       /* i  : coding type                                           */
    short *last_bin,        /* i  : last bin of bit allocation                            */
    float *lsf_new,         /* i  : ISFs at the end of the frame                          */
    float *exc_wo_nf,       /* o  : excitation (in f domain) without noisefill            */
    float *tmp_noise        /* o  : long-term noise energy                                */
)
{
    short i, j, bit, nb_subbands, pvq_len;
    short bitallocation_band[MBANDS_GN];
    short bitallocation_exc[2];
    float Ener_per_bd_iQ[MBANDS_GN];
    short max_ener_band[MBANDS_GN];
    float exc_diffQ[L_FRAME];
    short bits_per_bands[MBANDS_GN]; /*Q3*/
    float concat_out[L_FRAME];
    short npulses[NB_SFM];
    short maxpulse[NB_SFM];
    float mean_gain;
    short Mbands_gn = 16;
    short seed_init;

    /*--------------------------------------------------------------------------------------*
     * Initialization
     *--------------------------------------------------------------------------------------*/

    bit = bits_used;
    set_f( exc_diffQ, 0.0f, L_FRAME );

    /*--------------------------------------------------------------------------------------*
     * Gain decoding
     *--------------------------------------------------------------------------------------*/

    if( st->bfi || st->BER_detect )
    {
        /* copy old gain */
        mvr2r( st->old_y_gain, Ener_per_bd_iQ, Mbands_gn );
        mean_gain = st->lp_gainc/10.0f;
        for( i=0; i<Mbands_gn; i++ )
        {
            Ener_per_bd_iQ[i] += mean_gain;
        }

        st->lp_gainc *= 0.98f;
    }
    else
    {
        mean_gain = gsc_gaindec( st, Ener_per_bd_iQ, st->core_brate, st->old_y_gain, coder_type, st->bwidth );

        st->lp_gainc  = 10.0f * mean_gain;
    }

    *tmp_noise = st->lp_gainc;

    *last_bin = 0;
    if( st->core_brate == ACELP_8k00 && st->bwidth != NB )
    {
        bitallocation_exc[0] = 0;
        bitallocation_exc[1] = 0;
    }

    set_s( bitallocation_band, 0, MBANDS_GN );

    if( st->bfi || st->BER_detect )
    {
        /*--------------------------------------------------------------------------------------*
         * Copy old spectrum
         * reduce spectral dynamic
         * save spectrum
         *--------------------------------------------------------------------------------------*/

        if(st->last_good == INACTIVE_CLAS || st->Last_GSC_noisy_speech_flag == 1)
        {
            for( i=0; i<L_FRAME; i++ )
            {
                st->Last_GSC_spectrum[i] = 0.8f * own_random(&st->seed_tcx)/32768.0f+0.2f*st->Last_GSC_spectrum[i];
            }
            mvr2r( st->Last_GSC_spectrum, exc_diffQ, L_FRAME );
        }

        mvr2r( st->Last_GSC_spectrum, exc_diffQ, L_FRAME );

        for( i=0; i<L_FRAME; i++ )
        {
            st->Last_GSC_spectrum[i] *= 0.75f;
        }
    }
    else
    {
        /*--------------------------------------------------------------------------------------*
         * PVQ decoder
         *--------------------------------------------------------------------------------------*/

        bands_and_bit_alloc( st->cor_strong_limit, st->noise_lev, st->core_brate, Diff_len, bit, &bit, Ener_per_bd_iQ,
                             max_ener_band, bits_per_bands, &nb_subbands, NULL, NULL, &pvq_len, coder_type,
                             st->bwidth, st->GSC_noisy_speech );

        set_s( npulses, 0, NB_SFM );

        pvq_core_dec( st, gsc_sfm_start, gsc_sfm_end, gsc_sfm_size, concat_out, bit, nb_subbands, bits_per_bands, NULL, npulses, maxpulse, ACELP_CORE );

        seed_init = 0;
        for( j = 0; j < nb_subbands; j++ )
        {
            mvr2r( concat_out+j*16, exc_diffQ + max_ener_band[j]*16, 16);

            if( max_ener_band[j] > *last_bin )
            {
                *last_bin = max_ener_band[j];
            }

            bitallocation_band[max_ener_band[j]] = 1;

            seed_init += npulses[j];
        }
        if( st->last_coder_type != AUDIO /* First audio frame */
                && st->last_coder_type != UNVOICED  ) /* last_coder_type == INACTIVE is overwritten in update_dec to UNVOICED */
        {
            for(j = 0; j < nb_subbands*16; j++)
            {
                if( concat_out[j] > 0 )
                {
                    seed_init = (short)((int)seed_init<<3);
                }
                if( concat_out[j] < 0 )
                {
                    seed_init += 3;
                }
            }

            st->seed_tcx = seed_init;
        }

        if( st->core_brate == ACELP_8k00 && st->bwidth != NB )
        {
            if( exc_diffQ[L_FRAME8k - 2] != 0 )
            {
                bitallocation_exc[0] = 1;
            }

            if( exc_diffQ[L_FRAME8k - 1] != 0 )
            {
                bitallocation_exc[1] = 1;
            }
        }

        mvr2r( exc_diffQ, st->Last_GSC_spectrum, L_FRAME );

        /*--------------------------------------------------------------------------------------*
         * Skip adaptive (pitch) contribution frequency band (no noise added over the time contribution)
         * Find x pulses between 1.6-3.2kHz to code in the spectrum of the residual signal
         * Gain is based on the inter-correlation gain between the pulses found and residual signal
         *--------------------------------------------------------------------------------------*/

        freq_dnw_scaling( st->cor_strong_limit, coder_type, st->noise_lev, st->core_brate, exc_diffQ );
    }

    /*--------------------------------------------------------------------------------------*
     * Estimate noise level
     *--------------------------------------------------------------------------------------*/

    highband_exc_dct_in( st->core_brate, mfreq_bindiv_loc, *last_bin, Diff_len, st->noise_lev, pit_band_idx, exc_diffQ,
                         &st->seed_tcx, Ener_per_bd_iQ, nb_subfr, exc_dct_in, st->last_coder_type, bitallocation_band, lsf_new,
                         st->last_exc_dct_in, &st->last_ener, st->last_bitallocation_band, bitallocation_exc, st->bfi, coder_type,
                         st->bwidth, exc_wo_nf, st->GSC_noisy_speech
                         ,st->lt_ener_per_band
                       );

    exc_dct_in[0] = 0;

    return;
}
