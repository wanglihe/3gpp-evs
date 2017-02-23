/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"

/*-------------------------------------------------------------------*
 * Local functions
 *-------------------------------------------------------------------*/

static float edyn( const float *vec, const short lvec );

static void gsc_enc( Encoder_State *st, float res_dct_in[], float exc_dct_in[], const short Diff_len,
                     const short bits_used, const short nb_subfr, const short coder_type,
                     float *lsf_new, float* exc_wo_nf, float *tmp_noise );

/*-------------------------------------------------------------------*
 * encod_audio()
 *
 * Encode audio (AC) frames
 *-------------------------------------------------------------------*/

void encod_audio(
    Encoder_State *st,                   /* i/o: State structure                                  */
    LPD_state *mem,                  /* i/o: encoder memories                                 */
    const float speech[],              /* i  : input speech                                     */
    const float Aw[],                  /* i  : weighted A(z) unquantized for subframes          */
    const float Aq[],                  /* i  : 12k8 Lp coefficient                              */
    const short T_op[],                /* i  : open loop pitch                                  */
    const float voicing[],             /* i  : voicing                                          */
    const float *res,                  /* i  : residual signal                                  */
    float *synth,                /* i/o: core synthesis                                   */
    float *exc,                  /* i/o: current non-enhanced excitation                  */
    float *pitch_buf,            /* i/o: floating pitch values for each subframe          */
    float *voice_factors,        /* o  : voicing factors                                  */
    float *bwe_exc,              /* o  : excitation for SWB TBE                           */
    const short attack_flag,           /* i  : Flag that point to an attack coded with AC mode (GSC)    */
    const short coder_type,            /* i  : coding type                                      */
    float *lsf_new,              /* i  : current frame ISF vector                         */
    float *tmp_noise             /* o  : long-term noise energy                           */
)
{
    const float *p_Aq;
    short i, i_subfr, nb_subfr, last_pit_bin;
    short T0_tmp, T0_frac_tmp, nb_subfr_flag;
    short tmp_nb_bits_tot;
    float Es_pred;
    float dct_res[L_FRAME], dct_epit[L_FRAME];
    float m_mean;
    float exc_wo_nf[L_FRAME];
    short nb_bits;      /*number of bits*/
    int   indice;       /*parameter indices to write*/

    m_mean = 0.0f;
    tmp_nb_bits_tot = 0;

    T0_tmp = 64;
    T0_frac_tmp = 0;
    mvr2r( mem->mem_syn, st->mem_syn_tmp, M );
    st->mem_w0_tmp = mem->mem_w0;
    Es_pred = 0;

    /*---------------------------------------------------------------*
     * Encode GSC attack flag (used to reduce possible pre-echo)
     * Encode GSC SWB speech flag
     *---------------------------------------------------------------*/

    push_indice( st, IND_GSC_ATTACK, attack_flag, 1 );

    if( coder_type != INACTIVE && st->total_brate >= ACELP_13k20 )
    {
        push_indice( st,IND_GSC_SWB_SPEECH, st->GSC_noisy_speech, 1);
    }

    /*---------------------------------------------------------------*
     * Find and encode the number of subframes
     *---------------------------------------------------------------*/

    if ( st->core_brate >= ACELP_9k60 && st->core_brate <= ACELP_13k20 )
    {
        for( i = 0; i < 5; i++ )
        {
            if( fabs(st->gsc_lt_diff_etot[MAX_LT-i-1]) > 6.0f && st->cor_strong_limit == 1 )
            {
                st->cor_strong_limit = 0;
            }
        }
    }

    if( st->GSC_noisy_speech )
    {
        nb_subfr = NB_SUBFR;
        st->cor_strong_limit = 0;
        nb_subfr_flag = 1;
    }
    else
    {
        if( (st->cor_strong_limit == 0 || coder_type == INACTIVE) && st->core_brate >= ACELP_9k60 )
        {
            nb_subfr = 2;
            nb_subfr_flag = 0;
            st->cor_strong_limit = 0;
        }
        else
        {
            nb_subfr = SWNB_SUBFR;
            nb_subfr_flag = 1;
        }

        if( st->core_brate >= ACELP_9k60 )
        {
            /* nb_subfr_flag can only have the value 0 or 1 */
            push_indice( st, IND_HF_NOISE, nb_subfr_flag, 1);
        }
    }

    /*---------------------------------------------------------------*
     * Compute adaptive (pitch) excitation contribution
     *---------------------------------------------------------------*/

    if( st->GSC_noisy_speech && nb_subfr == NB_SUBFR )
    {
        nb_bits = Es_pred_bits_tbl[BIT_ALLOC_IDX(st->core_brate, GENERIC, -1, -1)];
        Es_pred_enc( &Es_pred, &indice, L_FRAME, L_SUBFR, res, voicing, nb_bits, 0 );
        push_indice( st, IND_ES_PRED, indice, nb_bits );
    }

    enc_pit_exc( st, mem, speech, Aw, Aq, Es_pred, T_op, voicing, res, synth, exc, &T0_tmp, &T0_frac_tmp, pitch_buf, nb_subfr, &st->lt_gpitch );

    /*---------------------------------------------------------------*
     * DCT transform
     *---------------------------------------------------------------*/

    edct( exc, dct_epit, L_FRAME );
    edct( res, dct_res, L_FRAME );

    /*---------------------------------------------------------------*
     * Calculate energy dynamics
     *---------------------------------------------------------------*/

    for( i = 7; i < 15; i++ )
    {
        m_mean += edyn( dct_res+i*16, 16 );
    }
    m_mean *= 0.125f;

    if( m_mean > st->mid_dyn )
    {
        st->mid_dyn = 0.2f * st->mid_dyn + 0.8f * m_mean;
    }
    else
    {
        st->mid_dyn = 0.6f * st->mid_dyn + 0.4f * m_mean;
    }

    if( coder_type != INACTIVE )
    {
        st->noise_lev = (NOISE_LEVEL_SP3+1) - usquant(st->mid_dyn, &m_mean, MIN_DYNAMIC, DYNAMIC_RANGE/GSC_NF_STEPS, GSC_NF_STEPS);
        if( st->noise_lev > NOISE_LEVEL_SP3)
        {
            st->noise_lev = NOISE_LEVEL_SP3;
        }
    }

    st->past_dyn_dec = st->noise_lev;

    if( st->core_brate <= ACELP_8k00 )
    {
        if( st->noise_lev <= NOISE_LEVEL_SP2 )
        {
            st->noise_lev = NOISE_LEVEL_SP2;
        }

        push_indice( st, IND_NOISE_LEVEL, st->noise_lev - NOISE_LEVEL_SP2, 2 );
    }
    else if( st->GSC_noisy_speech )
    {
        st->noise_lev = NOISE_LEVEL_SP3;
    }
    else
    {
        push_indice( st, IND_NOISE_LEVEL, st->noise_lev - NOISE_LEVEL_SP0, 3 );
    }

    /*---------------------------------------------------------------*
     * Find and encode the last band where the adaptive (pitch) contribution is significant
     *---------------------------------------------------------------*/

    last_pit_bin = Pit_exc_contribution_len( st, dct_res, dct_epit, pitch_buf, &st->pit_exc_hangover, coder_type );

    if( last_pit_bin == 0 )
    {
        mem->tilt_code = 0.0f;
    }
    else
    {
        last_pit_bin++;
    }

    /*--------------------------------------------------------------------------------------*
     * GSC encoder
     *--------------------------------------------------------------------------------------*/

    /* Find the current total number of bits used */
    tmp_nb_bits_tot = st->nb_bits_tot;

    if( st->extl_brate > 0 )
    {
        /* subtract 1 bit for TBE/BWE BWE flag (bit counted in extl_brate) */
        tmp_nb_bits_tot--;
    }

    if( coder_type == INACTIVE && st->core_brate <= ACELP_9k60 )
    {
        /* add 5 bits for noisiness */
        tmp_nb_bits_tot += 5;
    }

    gsc_enc( st, dct_res, dct_epit, last_pit_bin, tmp_nb_bits_tot, nb_subfr, coder_type, lsf_new, exc_wo_nf, tmp_noise );

    /*--------------------------------------------------------------------------------------*
     * iDCT transform
     *--------------------------------------------------------------------------------------*/

    edct( dct_epit, exc, L_FRAME );
    edct( exc_wo_nf, exc_wo_nf, L_FRAME );


    /*--------------------------------------------------------------------------------------*
     * Remove potential pre-echo in case an onset has been detected
     *--------------------------------------------------------------------------------------*/

    pre_echo_att( &st->Last_frame_ener, exc, attack_flag, st->last_coder_type );

    /*--------------------------------------------------------------------------------------*
     * Update BWE excitation
     *--------------------------------------------------------------------------------------*/

    set_f( voice_factors, 0.0f, NB_SUBFR);
    interp_code_5over2( exc, bwe_exc, L_FRAME );

    /*--------------------------------------------------------------------------------------*
     * Synthesis
     *--------------------------------------------------------------------------------------*/

    p_Aq = Aq;
    for ( i_subfr=0; i_subfr<L_FRAME; i_subfr+=L_SUBFR )
    {
        syn_filt( p_Aq, M, &exc_wo_nf[i_subfr], &synth[i_subfr], L_SUBFR, mem->mem_syn, 1 );
        p_Aq += (M+1);
    }

    /*--------------------------------------------------------------------------------------*
     * Updates
     *--------------------------------------------------------------------------------------*/

    mem->mem_w0 = st->mem_w0_tmp;
    mvr2r( exc_wo_nf, exc, L_FRAME );

    return;
}

/*-------------------------------------------------------------------*
 * gsc_enc()
 *
 * Generic audio signal encoder
 *-------------------------------------------------------------------*/

static void gsc_enc(
    Encoder_State *st,              /* i/o: State structure                                       */
    float res_dct_in[],     /* i  : dct of residual signal                                */
    float exc_dct_in[],     /* i/o: dct of pitch-only excitation / total excitation       */
    const short Diff_len,         /* i  : Lenght of the difference signal (before pure spectral)*/
    const short bits_used,        /* i  : Number of bit used before frequency Q                 */
    const short nb_subfr,         /* i  : Number of subframe considered                         */
    const short coder_type,       /* i  : coding type                                           */
    float *lsf_new,         /* i  : ISFs at the end of the frame                          */
    float *exc_wo_nf,       /* o  : excitation (in f domain) without noisefill            */
    float *tmp_noise        /* o  : long-term noise energy                                */
)
{
    short i;
    float exc_diffQ[L_FRAME];
    float exc_diff[L_FRAME];
    short bit;
    short nb_subbands;
    short pvq_len;
    short bits_per_bands[MBANDS_GN]; /*Q3*/
    short tmp_band;
    float concat_in[L_FRAME];
    float concat_out[L_FRAME];
    short max_ener_band[MBANDS_GN], j;
    float Ener_per_bd_iQ[MBANDS_GN];
    short last_bin;
    short bitallocation_band[MBANDS_GN];
    short bitallocation_exc[2];
    short npulses[NB_SFM];
    short maxpulse[NB_SFM];
    float mean_gain;
    short seed_init;

    /*--------------------------------------------------------------------------------------*
     * Initialization
     *--------------------------------------------------------------------------------------*/

    bit = bits_used;
    set_f( exc_diffQ, 0.0f, L_FRAME );

    /*--------------------------------------------------------------------------------------*
     * Calculate the difference between the residual spectrum and the spectrum of adaptive excitation
     * (non valuable temporal content present in exc_dct_in is already zeroed)
     *--------------------------------------------------------------------------------------*/

    v_sub( res_dct_in, exc_dct_in, exc_diff, L_FRAME );
    exc_diff[0] = 0;

    /*--------------------------------------------------------------------------------------*
     * Multiply the difference spectrum with the normalized spectral shape of the residual signal
     * This improves the stability of the differnece spectrum since the spectral shape of the
     * residual signal is less suseptible to rapid changes than the difference spectrum
     *--------------------------------------------------------------------------------------*/

    if( Diff_len == 0 )
    {
        tmp_band = 0;
    }
    else
    {
        tmp_band = st->mem_last_pit_band;
    }

    Ener_per_band_comp( exc_diff, Ener_per_bd_iQ, MBANDS_GN, 1 );

    /*--------------------------------------------------------------------------------------*
     * Gain quantizaion
     *--------------------------------------------------------------------------------------*/

    mean_gain = gsc_gainQ( st, Ener_per_bd_iQ, Ener_per_bd_iQ, st->core_brate, coder_type, st->bwidth );
    *tmp_noise = 10.0f * mean_gain;

    /*--------------------------------------------------------------------------------------*
     * PVQ encoder
     *--------------------------------------------------------------------------------------*/

    bands_and_bit_alloc( st->cor_strong_limit, st->noise_lev, st->core_brate, Diff_len, bit, &bit, Ener_per_bd_iQ,
                         max_ener_band, bits_per_bands, &nb_subbands, exc_diff, concat_in, &pvq_len,
                         coder_type, st->bwidth, st->GSC_noisy_speech );

    set_s( npulses, 0, NB_SFM );
    bit -= pvq_core_enc(st, concat_in, concat_out, bit, nb_subbands, gsc_sfm_start, gsc_sfm_end, gsc_sfm_size, bits_per_bands, NULL, npulses, maxpulse, ACELP_CORE );

    /* write unused bits */
    while( bit > 0 )
    {
        i = min( bit, 16 );
        push_indice( st, IND_UNUSED, 0, i );
        bit -= i;
    }

    /* Reorder Q bands */
    last_bin = 0;
    set_s( bitallocation_band, 0, MBANDS_GN );
    seed_init = 0;

    for( j = 0; j < nb_subbands; j++ )
    {
        mvr2r( concat_out+j*16, exc_diffQ + max_ener_band[j]*16, 16 );

        if( max_ener_band[j] > last_bin )
        {
            last_bin = max_ener_band[j];
        }

        bitallocation_band[max_ener_band[j]] = 1;

        seed_init += npulses[j];
    }
    if( st->last_coder_type != AUDIO /* First audio frame */
            && st->last_coder_type != UNVOICED  ) /* last_coder_type == INACTIVE is overwritten in update_enc to UNVOICED */
    {
        for( j = 0; j < nb_subbands*16; j++ )
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
        bitallocation_exc[0] = 0;
        bitallocation_exc[1] = 0;

        if( exc_diffQ[L_FRAME8k - 2] != 0 )
        {
            bitallocation_exc[0] = 1;
        }

        if( exc_diffQ[L_FRAME8k - 1] != 0 )
        {
            bitallocation_exc[1] = 1;
        }
    }

    /*--------------------------------------------------------------------------------------*
     * Skip adaptive (pitch) contribution frequency band (no noise added over the adaptive (pitch) contribution)
     * Find x pulses between 1.6-3.2kHz to code in the spectrum of the residual signal
     * Gain is based on the inter-correlation gain between the pulses found and residual signal
     *--------------------------------------------------------------------------------------*/

    freq_dnw_scaling( st->cor_strong_limit, coder_type, st->noise_lev, st->core_brate, exc_diffQ );

    /*--------------------------------------------------------------------------------------*
     * Estimate noise level
     *--------------------------------------------------------------------------------------*/

    highband_exc_dct_in( st->core_brate, mfreq_bindiv_loc, last_bin, Diff_len, st->noise_lev, tmp_band, exc_diffQ,
                         &st->seed_tcx, Ener_per_bd_iQ, nb_subfr, exc_dct_in, st->last_coder_type, bitallocation_band, lsf_new,
                         st->last_exc_dct_in, &st->last_ener, st->last_bitallocation_band, bitallocation_exc, 0, coder_type,
                         st->bwidth, exc_wo_nf, st->GSC_noisy_speech, NULL );

    exc_dct_in[0] = 0;

    return;
}

/*---------------------------------------------------------------------*
 * edyn()
 *
 * Calculate energy dynamics in a vector (ratio of energy maximum to energy mean)
 *---------------------------------------------------------------------*/

static float edyn(         /* o  : ratio of max to mean    */
    const float *vec,      /* i  : input vector            */
    const short lvec       /* i  : length of input vector  */
)
{
    short j;
    float temp, ener_max, ener_mean, dyn;

    ener_mean = 1.0f;
    ener_max = 1.0f;

    for( j=0; j<lvec; j++ )
    {
        temp = vec[j] * vec[j];

        if( temp > ener_max )
        {
            ener_max = temp;
        }
        ener_mean += temp;
    }
    ener_mean /= lvec;
    dyn = 10.0f * (ener_max / ener_mean);

    return dyn;
}
