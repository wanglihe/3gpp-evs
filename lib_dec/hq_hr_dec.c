/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"

/*--------------------------------------------------------------------------*
 * hq_pred_hb_bws()
 *
 * HQ core HB band-width switching handling
 *--------------------------------------------------------------------------*/

static void hq_pred_hb_bws(
    Decoder_State *st,            /* i/o: decoder state structure             */
    const short *ynrm,          /* i  : norm quantization index vector      */
    const short length,         /* i  : frame length                        */
    const short hqswb_clas,     /* i  : HQ SWB class                        */
    const float *SWB_fenv       /* i  : SWB frequency envelopes             */
)
{
    short i;

    /* SWB switching to WB */
    if ( length >= L_FRAME32k )      /* wb switch to swb */
    {
        /* calculate the switching parameters */
        if( (hqswb_clas != HQ_GEN_SWB && st->core_brate <= HQ_32k) || st->core_brate > HQ_32k )
        {
            st->prev_ener_shb = 0.0f;
            for( i=25; i<SFM_N_HARM; i++ )
            {
                st->prev_ener_shb += dicn[ynrm[i]];
            }
            st->prev_ener_shb /= 6;
        }
        else
        {
            st->prev_ener_shb = 0.0f;
            for( i=0; i<SWB_FENV-3; i++ )
            {
                st->prev_ener_shb += SWB_fenv[i];
            }
            st->prev_ener_shb /= (SWB_FENV-3);
        }
    }

    if( st->last_inner_frame >= L_FRAME32k )
    {
        set_f( st->prev_SWB_fenv, st->prev_ener_shb, SWB_FENV );
    }

    return;
}


/*--------------------------------------------------------------------------*
 * hq_hr_dec()
 *
 * HQ high rate decoding routine
 *--------------------------------------------------------------------------*/

void hq_hr_dec(
    Decoder_State *st,                    /* i/o: decoder state structure                 */
    float *t_audio_q,             /* o  : transform-domain coefficients           */
    const short length,                 /* i  : frame length                            */
    short num_bits,               /* i  : number of available bits                */
    short *ynrm,                  /* o  : norm quantization index vector          */
    short *is_transient,          /* o  : transient flag                          */
    short *hqswb_clas,            /* o  : HQ SWB class                            */
    float *SWB_fenv               /* o  : SWB frequency envelopes                 */
)
{
    short nb_sfm;
    short sum, hcode_l;
    const short *sfmsize, *sfm_start, *sfm_end;
    short num_sfm, numnrmibits;
    short nf_idx;
    short normqlg2[NB_SFM], R[NB_SFM];
    short pulses[NB_SFM], maxpulse[NB_SFM];
    float env_stab;
    short Rsubband[NB_SFM]; /*Q3*/
    short start_norm, Npeaks = 0;
    float noise_level[HVQ_BWE_NOISE_BANDS];
    short peak_idx[HVQ_MAX_PEAKS_32k];
    short hq_generic_offset;
    short num_env_bands;
    short hq_generic_exc_clas = 0;
    short core_sfm;
    short har_freq_est1, har_freq_est2;
    short flag_dis;
    const short *subband_search_offset;
    short wBands[2];
    short b_delta_env;
    short n_band;

    /*------------------------------------------------------------------*
     * Initializations
     *------------------------------------------------------------------*/
    set_s( pulses, 0, NB_SFM );
    set_s( maxpulse, 0, NB_SFM );
    flag_dis = 1;
    har_freq_est1 = 0;
    har_freq_est2 = 0;
    set_s( peak_idx, 0, HVQ_MAX_PEAKS_32k);

    /*------------------------------------------------------------------*
     * Decode classification
     *------------------------------------------------------------------*/

    num_bits -= hq_classifier_dec( st, st->core_brate, length, is_transient, hqswb_clas );


    /*------------------------------------------------------------------*
     * set quantization parameters
     *------------------------------------------------------------------*/

    hq_configure( length, *hqswb_clas, st->core_brate, &num_sfm, &nb_sfm, &start_norm,
                  &num_env_bands, &numnrmibits, &hq_generic_offset, &sfmsize, &sfm_start, &sfm_end );

    /*------------------------------------------------------------------*
     * Unpacking bit-stream
     *------------------------------------------------------------------*/

    nf_idx = 0;
    if( !*is_transient && *hqswb_clas != HQ_HVQ && !(length == L_FRAME16k && st->core_brate == HQ_32k))
    {
        nf_idx = (short)get_next_indice( st, 2 );
    }


    /*------------------------------------------------------------------*
     * Decode envelope
     *------------------------------------------------------------------*/

    hcode_l = decode_envelope_indices( st, start_norm, num_env_bands, numnrmibits, ynrm, NORMAL_HQ_CORE, *is_transient );
    num_bits -= hcode_l + NORM0_BITS + FLAGS_BITS;

    dequantize_norms( st, start_norm, num_env_bands, *is_transient, ynrm, normqlg2 );

    if ( *hqswb_clas == HQ_GEN_SWB || *hqswb_clas == HQ_GEN_FB )
    {
        hq_generic_exc_clas = swb_bwe_gain_deq( st, HQ_CORE, NULL, SWB_fenv, st->core_brate == HQ_32k, *hqswb_clas );
        if (hq_generic_exc_clas == HQ_GENERIC_SP_EXC)
        {
            num_bits++;        /* conditional 1 bit saving for representing FD3 BWE excitation class */
        }
        map_hq_generic_fenv_norm(*hqswb_clas, SWB_fenv, ynrm, normqlg2, num_env_bands, nb_sfm, hq_generic_offset);
    }

    env_stab = 0;
    if( *hqswb_clas == HQ_HVQ )
    {
        st->mem_env_delta = 0;
    }
    else if( length == L_FRAME32k )
    {
        env_stab = env_stability( ynrm, SFM_N_ENV_STAB, st->mem_norm, &st->mem_env_delta );
    }
    else
    {
        st->mem_norm[0] = 31;
        st->mem_env_delta = 0;
    }

    if( *hqswb_clas == HQ_HVQ )
    {
        st->env_stab = 1.0f;        /* stable by definition */
    }
    else
    {
        if( length == L_FRAME32k )
        {
            st->env_stab = env_stab;    /* calculated stability */
        }
        else
        {
            st->env_stab = env_stability( ynrm, SFM_N_ENV_STAB_WB, st->mem_norm_hqfec, &st->mem_env_delta_hqfec );
        }
    }
    st->env_stab_plc=env_stab_smo(min(st->env_stab,1.0f-stab_trans[L_STAB_TBL-1]),st->env_stab_state_p,&st->envstabplc_hocnt);

    /*------------------------------------------------------------------*
     * Bit allocation
     *------------------------------------------------------------------*/

    hq_bit_allocation( st->core_brate, length, *hqswb_clas, &num_bits, normqlg2, nb_sfm, sfmsize, noise_level,
                       R, Rsubband, &sum, &core_sfm, num_env_bands );

    if( st->bws_cnt1 > 0 && *hqswb_clas == HQ_GEN_SWB && st->core_brate == HQ_24k40 )
    {
        if(st->L_frame == L_FRAME16k )
        {
            for (n_band = 0; n_band < 4; n_band++)
            {
                SWB_fenv[n_band] *= (float)st->bws_cnt1 / (float)N_NS2W_FRAMES;
            }
        }

        for (n_band = 4; n_band < SWB_FENV; n_band++)
        {
            SWB_fenv[n_band] *= (float)st->bws_cnt1 / (float)N_NS2W_FRAMES;
        }
    }

    if ( *hqswb_clas == HQ_GEN_SWB || *hqswb_clas == HQ_GEN_FB )
    {
        b_delta_env = get_nor_delta_hf(st, ynrm, Rsubband, num_env_bands, nb_sfm, core_sfm );
        sum -= b_delta_env;
    }

    /*------------------------------------------------------------------*
     * Decode spectral fine structure using HVQ/PVQ
     *------------------------------------------------------------------*/

    if( *hqswb_clas == HQ_HVQ )
    {
        hvq_dec( st, num_bits, st->core_brate, ynrm, R, noise_level, peak_idx, &Npeaks, t_audio_q, st->core );
    }
    else
    {
        pvq_core_dec(st, sfm_start, sfm_end, sfmsize, t_audio_q, sum, nb_sfm, Rsubband, R, pulses, maxpulse, HQ_CORE );
    }

    if( *hqswb_clas == HQ_HVQ || *hqswb_clas == HQ_HARMONIC )
    {
        subband_search_offset = subband_search_offsets_13p2kbps_Har;
        wBands[0] = SWB_SB_BW_LEN0_16KBPS_HAR;
        wBands[1] = SWB_SB_BW_LEN1_16KBPS_HAR;

        har_est( t_audio_q, 300, &har_freq_est1, &har_freq_est2, &flag_dis, &st->prev_frm_hfe2, subband_search_offset, wBands, &st->prev_stab_hfe2 );

        st->prev_frm_hfe2 = har_freq_est2;
    }

    if( *hqswb_clas != HQ_HARMONIC || *hqswb_clas != HQ_HVQ || flag_dis == 0 )
    {
        st->prev_frm_hfe2 = 0; /*reset*/
        st->prev_stab_hfe2 = 0;/*reset*/
    }

    /*------------------------------------------------------------------*
     * Spectral filling
     *------------------------------------------------------------------*/
    fill_spectrum( t_audio_q, R, *is_transient, ynrm, SWB_fenv, hq_generic_offset, nf_idx, length, env_stab,
                   &st->no_att_hangover, &st->energy_lt, &st->hq_generic_seed, hq_generic_exc_clas,
                   core_sfm, *hqswb_clas, noise_level, st->core_brate, st->prev_noise_level, &(st->prev_R), st->prev_coeff_out,
                   peak_idx, Npeaks, pulses, st->old_is_transient[0], st->prev_normq, st->prev_env, st->prev_bfi,
                   sfmsize, sfm_start, sfm_end, &st->prev_L_swb_norm, st->prev_hqswb_clas, num_sfm, num_env_bands );

    enforce_zero_for_min_envelope( *hqswb_clas, ynrm, t_audio_q, nb_sfm, sfm_start, sfm_end );

    if( *is_transient )
    {
        de_interleave_spectrum( t_audio_q, length );
    }

    /*------------------------------------------------------------------*
     * WB/SWB bandwidth switching
     *------------------------------------------------------------------*/

    hq_pred_hb_bws( st, ynrm, length, *hqswb_clas, SWB_fenv );

    /* update */
    st->prev_hqswb_clas = *hqswb_clas;

    return;
}
