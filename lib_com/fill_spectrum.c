/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <stdlib.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"

/*--------------------------------------------------------------------------*
 * fill_spectrum()
 *
 * Apply spectral filling by
 * - filling zero-bit bands below BWE region
 * - applying BWE above transition frequency
 *--------------------------------------------------------------------------*/

void fill_spectrum(
    float *coeff,               /* i/o: normalized MLT spectrum / nf spectrum                 */
    short *R,                   /* i  : number of pulses per band                             */
    const short is_transient,         /* i  : transient flag                                        */
    short norm[],               /* i  : quantization indices for norms                        */
    const float *hq_generic_fenv,     /* i  : HQ GENERIC envelope                                   */
    const short hq_generic_offset,    /* i  : HQ GENERIC offset                                     */
    const short nf_idx,               /* i  : noise fill index                                      */
    const short length,               /* i  : Length of spectrum (32 or 48 kHz)                     */
    const float env_stab,             /* i  : Envelope stability measure [0..1]                     */
    short *no_att_hangover,     /* i/o: Frame counter for attenuation hangover                */
    float *energy_lt,           /* i/o: Long-term energy measure for transient detection      */
    short *bwe_seed,            /* i/o: random seed for generating BWE input                  */
    const short hq_generic_exc_clas,  /* i  : HQ generic hf excitation class                        */
    const short core_sfm,             /* i  : index of the end band for core                        */
    short HQ_mode,              /* i  : HQ mode                                               */
    float noise_level[],        /* i  : noise levels for harmonic modes                       */
    long  core_brate,           /* i  : target bit-rate                                       */
    float prev_noise_level[],   /* i/o: noise factor in previous frame                        */
    short *prev_R,              /* i/o: bit allocation info. in previous frame                */
    float *prev_coeff_out,      /* i/o: decoded spectrum in previous frame                    */
    const short *peak_idx,            /* i  : HVQ peak indices                                      */
    const short Npeaks,               /* i  : Number of HVQ peaks                                   */
    const short *npulses,             /* i  : Number of assigned pulses per band                    */
    short prev_is_transient,    /* i  : previous transient flag                               */
    float *prev_normq,          /* i  : previous norms                                        */
    float *prev_env,            /* i  : previous noise envelopes                              */
    short prev_bfi,             /* i  : previous bad frame indicator                          */
    const short *sfmsize,             /* i  : Length of bands                                       */
    const short *sfm_start,           /* i  : Start of bands                                        */
    const short *sfm_end,             /* i  : End of bands                                          */
    short *prev_L_swb_norm,     /* i/o: last normalize length for harmonic mode               */
    short prev_hq_mode,         /* i  : previous HQ mode                                      */
    const short num_sfm,              /* i  : Number of bands                                       */
    const short num_env_bands         /* i  : Number of envelope bands                              */
)
{
    float CodeBook[FREQ_LENGTH];
    short cb_size = 0;
    short last_sfm;
    float CodeBook_mod[FREQ_LENGTH];
    float norm_adj[NB_SFM];
    short high_sfm = 23;
    short flag_32K_env_hangover;
    short bin_th = 0;
    short peak_pos[L_HARMONIC_EXC];
    short bwe_peaks[L_FRAME48k];
    float normq_v[NB_SFM];
    float coeff_fine[L_FRAME48k];
    float coeff_out[L_FRAME48k];

    set_s( peak_pos, 0, L_HARMONIC_EXC );
    set_s( bwe_peaks, 0, L_FRAME48k );
    set_f( norm_adj, 1.0f, num_sfm );
    set_f( coeff_out, 0.0f, L_FRAME48k );

    if ( HQ_mode == HQ_TRANSIENT )
    {
        last_sfm = num_sfm-1;
    }
    else if (HQ_mode == HQ_GEN_SWB || HQ_mode == HQ_GEN_FB)
    {
        last_sfm = max(core_sfm,num_env_bands-1);
    }
    else
    {
        last_sfm = core_sfm;
    }

    if ( HQ_mode == HQ_HARMONIC )
    {
        high_sfm = (core_brate == HQ_24k40) ? HVQ_THRES_SFM_24k-1 : HVQ_THRES_SFM_32k-1;
        if( last_sfm < high_sfm )
        {
            last_sfm = high_sfm;
        }
    }
    else if ( HQ_mode == HQ_HVQ )
    {
        bin_th = sfm_end[last_sfm];
    }

    /* Transient analysis for envelope stability measure */
    if( length == L_FRAME32k )
    {
        env_stab_transient_detect( is_transient, length, norm, no_att_hangover, energy_lt, HQ_mode, bin_th, coeff );
    }

    if (  length == L_FRAME16k || ((length == L_FRAME32k && HQ_mode != HQ_HARMONIC && HQ_mode != HQ_HVQ) && *no_att_hangover == 0) )
    {
        /* Norm adjustment function */
        env_adj( npulses, length, last_sfm, norm_adj, env_stab, sfmsize );
    }

    flag_32K_env_hangover = ( length == L_FRAME32k && ( (env_stab < 0.5f && *no_att_hangover == 0) || HQ_mode == HQ_HVQ ) );

    /*----------------------------------------------------------------*
     * Build noise-fill codebook
     *----------------------------------------------------------------*/

    if ( HQ_mode != HQ_HVQ )
    {
        cb_size = build_nf_codebook( flag_32K_env_hangover, coeff, sfm_start, sfmsize, sfm_end, last_sfm, R, CodeBook, CodeBook_mod );
    }

    /*----------------------------------------------------------------*
     * Prepare fine structure for Harmonic and HVQ
     *----------------------------------------------------------------*/

    if ( HQ_mode == HQ_HARMONIC )
    {
        harm_bwe_fine( R, last_sfm, high_sfm, num_sfm, norm, sfm_start, sfm_end, prev_L_swb_norm, coeff, coeff_out, coeff_fine );
    }
    else if ( HQ_mode == HQ_HVQ )
    {
        hvq_bwe_fine( last_sfm, num_sfm, sfm_end, peak_idx, Npeaks, peak_pos, prev_L_swb_norm, coeff, bwe_peaks, coeff_fine );
    }

    /*----------------------------------------------------------------*
     * Apply noise-fill
     *----------------------------------------------------------------*/

    if( HQ_mode != HQ_HVQ && cb_size > 0 )
    {
        apply_noisefill_HQ( R, length, flag_32K_env_hangover, core_brate, last_sfm, CodeBook, CodeBook_mod, cb_size, sfm_start, sfm_end, sfmsize, coeff );
    }

    /*----------------------------------------------------------------*
     * Normal mode BWE
     *----------------------------------------------------------------*/

    if ( HQ_mode == HQ_NORMAL )
    {
        hq_fold_bwe( last_sfm, sfm_end, num_sfm, coeff );
    }

    /*----------------------------------------------------------------*
     * Apply noise-fill adjustment
     *----------------------------------------------------------------*/

    if( (length >= L_FRAME32k || core_brate > HQ_32k || core_brate < HQ_24k40) && HQ_mode != HQ_HVQ )
    {
        apply_nf_gain( nf_idx, last_sfm, R, sfm_start, sfm_end, coeff );
    }

    /*----------------------------------------------------------------*
     * Prepare fine strucutre for HQ GENERIC
     *----------------------------------------------------------------*/
    if ( HQ_mode == HQ_GEN_SWB || HQ_mode == HQ_GEN_FB )
    {
        hq_generic_fine( coeff, last_sfm, sfm_start, sfm_end, bwe_seed, coeff_fine );
    }

    /*----------------------------------------------------------------*
     * Apply envelope
     *----------------------------------------------------------------*/

    if ( HQ_mode != HQ_HARMONIC && HQ_mode != HQ_HVQ )
    {
        apply_envelope( coeff, norm, norm_adj, num_sfm, last_sfm, HQ_mode, length, sfm_start, sfm_end, normq_v, coeff_out, coeff_fine );
    }

    /*----------------------------------------------------------------*
     * Harmonic BWE, HVQ BWE and HQ SWB BWE
     *----------------------------------------------------------------*/

    if( HQ_mode == HQ_HARMONIC )
    {
        harm_bwe( coeff_fine, coeff, num_sfm, sfm_start, sfm_end, last_sfm, high_sfm, R, prev_hq_mode, norm, noise_level, prev_noise_level, bwe_seed, coeff_out );
    }
    else if ( HQ_mode == HQ_HVQ )
    {
        hvq_bwe( coeff, coeff_fine, sfm_start, sfm_end, sfmsize, last_sfm,
                 prev_hq_mode, bwe_peaks, bin_th, num_sfm, core_brate, R, norm,
                 noise_level, prev_noise_level, bwe_seed, coeff_out );
    }
    else if ( HQ_mode == HQ_GEN_SWB || HQ_mode == HQ_GEN_FB )
    {
        hq_generic_bwe( HQ_mode, coeff_fine, hq_generic_fenv, coeff_out, hq_generic_offset, prev_L_swb_norm,
                        hq_generic_exc_clas, sfm_end, num_sfm, num_env_bands, R );
    }

    /*----------------------------------------------------------------*
     * HQ WB BWE refinements
     *----------------------------------------------------------------*/

    if( length == L_FRAME16k && core_brate == HQ_32k )
    {
        hq_wb_nf_bwe( coeff, is_transient, prev_bfi, normq_v, num_sfm, sfm_start, sfm_end, sfmsize, last_sfm, R,
                      prev_is_transient, prev_normq, prev_env, bwe_seed, prev_coeff_out, prev_R, coeff_out );
    }

    /*----------------------------------------------------------------*
     * Update memories
     *----------------------------------------------------------------*/

    if ( HQ_mode != HQ_HARMONIC && HQ_mode != HQ_HVQ )
    {
        prev_noise_level[0] = 0.1f;
        prev_noise_level[1] = 0.1f;
    }
    if ( !(length == L_FRAME16k && core_brate == HQ_32k ) )
    {
        set_f( prev_env, 0, SFM_N_WB );
        set_f( prev_normq, 0, SFM_N_WB );
    }

    if ( length == L_FRAME32k && core_brate <= HQ_32k )
    {
        *prev_R = R[SFM_N_WB-1];
        mvr2r( coeff_out + L_FRAME16k - L_HQ_WB_BWE, prev_coeff_out, L_HQ_WB_BWE );
    }

    mvr2r( coeff_out, coeff, L_FRAME48k );

    return;
}
