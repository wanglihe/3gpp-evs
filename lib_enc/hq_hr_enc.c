/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_enc.h"
#include "rom_com.h"

/*--------------------------------------------------------------------------*
 * hq_hr_enc()
 *
 * HQ high rate encoding routine
 *--------------------------------------------------------------------------*/

void hq_hr_enc(
    Encoder_State *st,                /* i/o: encoder state structure         */
    float *t_audio,           /* i/o: transform-domain coefficients   */
    const short length,             /* i  : length of spectrum              */
    short *num_bits,          /* i  : number of available bits        */
    const short is_transient        /* i  : transient flag                  */
)
{
    short nb_sfm;
    short sum, hcode_l;
    short difidx[NB_SFM];
    short normqlg2[NB_SFM], ynrm[NB_SFM];
    short nf_idx;
    short LCmode;
    short shape_bits, num_sfm, numnrmibits;
    short hqswb_clas;
    short num_env_bands;
    short Npeaks, start_norm;
    short difidx_org[NB_SFM];
    short R[NB_SFM];
    short peaks[HVQ_MAX_PEAKS];
    const short *sfmsize, *sfm_start, *sfm_end;
    short npulses[NB_SFM], maxpulse[NB_SFM];
    short Rsubband[NB_SFM]; /* Q3 */
    float t_audio_q[L_FRAME48k];
    float nf_gains[HVQ_NF_GROUPS], pe_gains[HVQ_NF_GROUPS];
    float noise_level[HVQ_BWE_NOISE_BANDS];
    short hq_generic_offset;
    float hq_generic_fenv[HQ_FB_FENV];
    short hq_generic_exc_clas = 0;
    short core_sfm;
    short har_freq_est1 = 0, har_freq_est2 = 0;
    short flag_dis = 1;
    const short *subband_search_offset;
    short wBands[2];
    short b_delta_env;

    /*------------------------------------------------------------------*
     * Initializations
     *------------------------------------------------------------------*/

    Npeaks = 0;
    set_s( npulses, 0, NB_SFM );
    set_s( maxpulse, 0, NB_SFM );
    set_s( difidx_org, 0, NB_SFM );
    set_f( t_audio_q, 0.0f, L_FRAME48k );
    set_f( nf_gains, 0.0f, HVQ_NF_GROUPS );
    set_f( pe_gains, 0.0f, HVQ_NF_GROUPS );

    /*------------------------------------------------------------------*
     * Classification
     *------------------------------------------------------------------*/

    *num_bits -= hq_classifier_enc( st, length, t_audio, is_transient, &Npeaks, peaks, pe_gains, nf_gains, &hqswb_clas );

    /*------------------------------------------------------------------*
     * Set quantization parameters
     *------------------------------------------------------------------*/

    hq_configure( length, hqswb_clas, st->core_brate, &num_sfm, &nb_sfm, &start_norm, &num_env_bands, &numnrmibits, &hq_generic_offset,
                  &sfmsize, &sfm_start, &sfm_end );

    /*------------------------------------------------------------------*
     * Transient frame handling
     *------------------------------------------------------------------*/

    /* Interleave MLT coefficients of 4 sub-vectors in case of transient frame */
    if( is_transient )
    {
        interleave_spectrum( t_audio, length );
    }

    /*------------------------------------------------------------------*
     * Scalar quantization of norms
     * Encode norm indices
     *------------------------------------------------------------------*/

    /* calculate and quantize norms */
    calc_norm( t_audio, ynrm, normqlg2, start_norm, num_env_bands, sfmsize, sfm_start );

    /* create differential code of quantized norm indices */
    diff_envelope_coding( is_transient, num_env_bands, start_norm, ynrm, normqlg2, difidx );

    /* Find norm coding mode and calculate number of bits */
    hcode_l = encode_envelope_indices( st, num_env_bands, numnrmibits, difidx, &LCmode, 0, NORMAL_HQ_CORE, is_transient );
    *num_bits -= hcode_l + NORM0_BITS + FLAGS_BITS;

    /* Encode norm indices */
    encode_envelope_indices( st, num_env_bands, numnrmibits, difidx, &LCmode, 1, NORMAL_HQ_CORE, is_transient );

    /*------------------------------------------------------------------*
     * HQ GENERIC BWE encoding
     *------------------------------------------------------------------*/

    if ( hqswb_clas == HQ_GEN_SWB || hqswb_clas == HQ_GEN_FB )
    {
        hq_generic_hf_encoding( t_audio, hq_generic_fenv, hq_generic_offset, st, &hq_generic_exc_clas );
        if (hq_generic_exc_clas == HQ_GENERIC_SP_EXC)
        {
            (*num_bits)++;        /* conditional 1 bit saving for representing FD3 BWE excitation class */
        }
        map_hq_generic_fenv_norm(hqswb_clas, hq_generic_fenv, ynrm, normqlg2, num_env_bands, nb_sfm, hq_generic_offset);
    }

    /*------------------------------------------------------------------*
     * Bit allocation
     *------------------------------------------------------------------*/

    hq_bit_allocation( st->core_brate, length, hqswb_clas, num_bits, normqlg2,
                       nb_sfm, sfmsize, noise_level, R, Rsubband, &sum,
                       &core_sfm, num_env_bands);

    /*------------------------------------------------------------------*
     * Normalize coefficients with quantized norms
     *------------------------------------------------------------------*/
    if( hqswb_clas != HQ_HVQ )
    {
        if (hqswb_clas == HQ_GEN_SWB || hqswb_clas == HQ_GEN_FB)
        {
            b_delta_env = calc_nor_delta_hf( st, t_audio, ynrm, Rsubband, num_env_bands, nb_sfm, sfmsize, sfm_start, core_sfm );
            sum -= b_delta_env;
        }
        normalizecoefs( t_audio, ynrm, nb_sfm, sfm_start, sfm_end );
    }

    /*------------------------------------------------------------------*
     * Quantize/code spectral fine structure using PVQ or HVQ
     *------------------------------------------------------------------*/

    if( hqswb_clas == HQ_HVQ )
    {
        sum = hvq_enc( st, st->core_brate, *num_bits, Npeaks, ynrm, R, peaks, nf_gains,
                       noise_level, pe_gains, t_audio, t_audio_q );
        *num_bits -= sum;
    }
    else
    {
        shape_bits = pvq_core_enc( st, t_audio, t_audio_q, sum, nb_sfm, sfm_start, sfm_end, sfmsize, Rsubband, R, npulses, maxpulse, HQ_CORE );
        *num_bits += (sum - shape_bits);
    }

    if( hqswb_clas == HQ_HVQ || hqswb_clas == HQ_HARMONIC )
    {
        subband_search_offset = subband_search_offsets_13p2kbps_Har;
        wBands[0] = SWB_SB_BW_LEN0_16KBPS_HAR;
        wBands[1] = SWB_SB_BW_LEN1_16KBPS_HAR;

        har_est( t_audio_q, 300, &har_freq_est1, &har_freq_est2, &flag_dis, &st->prev_frm_hfe2, subband_search_offset, wBands, &st->prev_stab_hfe2 );

        st->prev_frm_hfe2 = har_freq_est2;
    }

    if( hqswb_clas != HQ_HARMONIC || hqswb_clas != HQ_HVQ || flag_dis == 0)
    {
        st->prev_frm_hfe2 = 0; /*reset*/
        st->prev_stab_hfe2 = 0; /*reset*/
    }

    nf_idx = 0;
    if ( !is_transient && hqswb_clas != HQ_HVQ && !(length == L_FRAME16k && st->core_brate == HQ_32k) )
    {
        if (hqswb_clas == HQ_GEN_SWB || hqswb_clas == HQ_GEN_FB)
        {
            nf_idx = noise_adjust( t_audio, R, sfm_start, sfm_end, max(core_sfm,num_env_bands-1));
            push_indice( st, IND_NF_IDX, nf_idx, 2 );
        }
        else
        {
            nf_idx = noise_adjust( t_audio, R, sfm_start, sfm_end, core_sfm );
            push_indice( st, IND_NF_IDX, nf_idx, 2 );
        }
    }


    /* updates */
    st->prev_hqswb_clas = hqswb_clas;

    return;
}
