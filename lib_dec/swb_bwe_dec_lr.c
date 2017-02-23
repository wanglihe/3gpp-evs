/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"
#include "stat_com.h"



/*-------------------------------------------------------------------*
 * DecodeSWBGenericParameters()
 *
 * Decoding of generic subband coding parameters
 *-------------------------------------------------------------------*/

static void DecodeSWBGenericParameters(
    Decoder_State *st,                    /* i/o: decoder state structure   */
    short *lagIndices,            /* o  : lowband index for each subband              */
    const short nBands_search,          /* i  : number of subbnads for SSearch              */
    const short BANDS,                  /* i  : total number of subbands per frame          */
    const short *p2a_flags,             /* i  : HF tonal flag                               */
    const short hq_swb_clas             /* i  : mode of operation HQ_NORMAL or HQ_HARMONIC  */
)
{
    short sb;

    /* lag index for each subband (except last two) */
    for (sb = 0; sb < nBands_search; sb++)
    {
        if( hq_swb_clas == HQ_HARMONIC )
        {
            lagIndices[sb] = (short)get_next_indice( st, bits_lagIndices_mode0_Har[sb]);
        }
        else
        {
            if( p2a_flags[BANDS-NB_SWB_SUBBANDS+sb] == 0 )
            {
                lagIndices[sb] = (short)get_next_indice( st, bits_lagIndices[sb]);
            }
            else
            {
                lagIndices[sb] = 0;
            }
        }
    }

    return;
}


/*-------------------------------------------------------------------*
 * DecodeSWBSubbands()
 *
 * Main routine for generic SWB coding
 *
 * High-frequency subbands are replicated based on the lowband signal using a lowband index denoting
 * the selected lowband subband as well as linear and logarithmic domain gains
 *-------------------------------------------------------------------*/

static void DecodeSWBSubbands(
    Decoder_State *st,            /* i/o: decoder state structure   */
    float  *spectra,              /* i/o: MDCT domain spectrum             */
    const short fLenLow,                /* i  : lowband length                   */
    const short fLenHigh,               /* i  : highband length                  */
    const short nBands,                 /* i  : number of subbands               */
    const short *sbWidth,               /* i  : subband lengths                  */
    short *lagIndices,
    float *lagGains,              /* i  : first gain for each subband      */
    short BANDS,                  /* i  : number subbands per frame        */
    short *band_start,            /* i  : band start of each SB            */
    short *band_end,              /* i  : band end of each SB              */
    float *band_energy,           /* i  : band energy of each SB           */
    short *p2a_flags,             /* i  : HF tonal indicator               */
    const short hqswb_clas,             /* i  : class information                */
    const short har_bands,              /* i  : number of LF harmonic bands      */
    const short *subband_search_offset,
    short *prev_frm_hfe2,
    short *prev_stab_hfe2,
    short band_width[],           /* i   :  subband band widths                */
    const short *subband_offsets,       /* i   :  subband offsets for sparse filling */
    const float spectra_ni[],           /* i   :  core coder with sparseness filled  */
    short *ni_seed                /* i/o: random seed for search buffer NI     */
)
{
    short i;
    short k;
    float sspectra[L_FRAME32k];
    float sspectra_ni[L_FRAME32k],sspectra_diff[L_FRAME32k], th_g[NB_SWB_SUBBANDS];
    float ss_min=1.0f,g ,be_tonal[SWB_HAR_RAN1], xSynth_har[L_FRAME32k];
    GainItem pk_sf[(NB_SWB_SUBBANDS)*8];
    short lagIndices_real[NB_SWB_SUBBANDS];
    short pul_res[NB_SWB_SUBBANDS],cnt,imin;
    short har_freq_est1 = 0;
    short har_freq_est2 = 0;
    short flag_dis = 1;
    short pos_max_hfe2=0;

    set_s(pul_res,0,NB_SWB_SUBBANDS);
    set_f( xSynth_har, 0.0f, fLenHigh );


    if( hqswb_clas == HQ_HARMONIC )
    {
        /* Harmonic Structure analysis */
        pos_max_hfe2 = har_est( spectra, fLenLow, &har_freq_est1, &har_freq_est2, &flag_dis, prev_frm_hfe2, subband_search_offset, sbWidth, prev_stab_hfe2 );
        /* Spectrum normalization for the corecoder */
        noise_extr_corcod( spectra, spectra_ni, sspectra, sspectra_diff, sspectra_ni, fLenLow, st->prev_hqswb_clas, &st->prev_ni_ratio );
        /* Harmonic Structure analysis */

        if(flag_dis == 0)
        {
            if(har_freq_est2 != SWB_HAR_RAN1 || har_freq_est2 != *prev_frm_hfe2)
            {
                har_freq_est2 +=lagIndices[0];
            }
        }
        /*Generate HF noise*/
        genhf_noise( sspectra_diff, xSynth_har, sspectra, BANDS, har_bands, har_freq_est2, pos_max_hfe2,pul_res,pk_sf,fLenLow,
                     fLenHigh, sbWidth, lagIndices, subband_offsets, subband_search_offset );

        imin =(short) get_next_indice( st, 2);
        g=(float) pow (10.0f,gain_table[imin]);
        /* tonal energy estimation*/
        ton_ene_est( xSynth_har, be_tonal, band_energy, band_start, band_end, band_width, fLenLow, fLenHigh, BANDS, har_bands, g, pk_sf, pul_res );

        /*HF Spectrum Generation*/
        Gettonl_scalfact( xSynth_har, spectra_ni, fLenLow, fLenHigh, har_bands, BANDS, band_energy, band_start, band_end, p2a_flags, be_tonal, pk_sf, pul_res);
        if(flag_dis == 0)
        {
            *prev_frm_hfe2 = 0;
        }
        else
        {
            *prev_frm_hfe2 = har_freq_est2;
        }

        for( k = har_bands; k<BANDS; k++ )
        {
            for( i = band_start[k]; i <= band_end[k]; i++ )
            {
                spectra[i] = xSynth_har[i-fLenLow];
            }
        }
    }
    else if(hqswb_clas == HQ_NORMAL)
    {
        /* Spectrum normalization for the corecoder */
        ss_min = spectrumsmooth_noiseton(spectra,spectra_ni,sspectra,sspectra_diff,sspectra_ni,fLenLow, ni_seed);

        convert_lagIndices_pls2smp( (short *)lagIndices, nBands, lagIndices_real, sspectra, sbWidth, fLenLow );

        for (k = 0; k < nBands; k++)
        {
            if( p2a_flags[BANDS-NB_SWB_SUBBANDS+k] == 1 )
            {
                lagIndices_real[k] = 0;
            }
        }
        /*get levels for missing bands*/
        GetlagGains( sspectra_ni, &band_energy[BANDS-NB_SWB_SUBBANDS], nBands, sbWidth, lagIndices_real, fLenLow, lagGains );
        for(k=0; k<NB_SWB_SUBBANDS; k++)
        {
            lagGains[k] *= 0.9f;
        }
        cnt = 0;
        for(k=0; k<NB_SWB_SUBBANDS; k++)
        {
            th_g[k] = 0.0f;
            if(p2a_flags[BANDS-NB_SWB_SUBBANDS+k] == 0)
            {
                th_g[k] = lagGains[k]*ss_min;
                cnt++;
            }
        }
        /* Construct spectrum */
        GetSynthesizedSpecThinOut( sspectra_ni, xSynth_har, nBands, sbWidth, lagIndices_real, lagGains, fLenLow );

        /*Level adjustment for the missing bands*/
        noiseinj_hf( xSynth_har, th_g, band_energy, st->prev_En_sb, p2a_flags, BANDS, band_start, band_end, fLenLow );
        /* xSynth is the reconstructed high-band */
        for( k = BANDS-NB_SWB_SUBBANDS; k<BANDS; k++ )
        {
            if( p2a_flags[k] == 0 )
            {
                for( i = band_start[k]; i <= band_end[k]; i++ )
                {
                    spectra[i] = xSynth_har[i-fLenLow];
                }
            }
            else
            {
                for( i = band_start[k]; i <= band_end[k]; i++ )
                {
                    spectra[i] = spectra_ni[i];
                }
            }
        }
    }

    return;
}


/*-------------------------------------------------------------------*
 * swb_bwe_dec_lr()
 *
 * Main decoding routine of SWB BWE for the LR MDCT core
 *-------------------------------------------------------------------*/

void swb_bwe_dec_lr(
    Decoder_State *st,                 /* i/o: decoder state structure   */
    const float m_core[],            /* i  : lowband synthesis                       */
    float m[],                 /* o  : highband synthesis with lowband zeroed  */
    const long  total_brate,         /* i  : total bitrate for selecting subband pattern */
    short BANDS,               /* i  : Number subbands/Frame */
    short *band_start,         /* i  : Band Start of each SB */
    short *band_end,           /* i  : Band end of each SB */
    float *band_energy,        /* i  : BAnd energy of each SB */
    short *p2a_flags,          /* i  : HF tonal Indicator */
    const short hqswb_clas,          /* i  : class information */
    short lowlength,           /* i  : Lowband Length */
    short highlength,          /* i  : Highband Length */
    const short har_bands,           /* i  : Number of LF harmonic bands */
    short *prev_frm_hfe2,
    short *prev_stab_hfe2,
    short band_width[],        /* i  : subband bandwidth                           */
    const float y2_ni[],             /* i/o: Sparse filled corecoder                     */
    short *ni_seed             /* i/o: random seed */
)
{
    short k;
    short nBands;
    short nBands_search;
    short wBands[NB_SWB_SUBBANDS];
    short lagIndices[NB_SWB_SUBBANDS];
    float lagGains[NB_SWB_SUBBANDS];
    short swb_lowband, swb_highband;
    const short *subband_search_offset;

    const short *subband_offsets;

    subband_search_offset = subband_search_offsets_13p2kbps_Har;
    subband_offsets = subband_offsets_sub5_13p2kbps_Har;
    hf_parinitiz(total_brate,hqswb_clas,lowlength,highlength,wBands,&subband_search_offset,&subband_offsets,&nBands,&nBands_search,&swb_lowband,&swb_highband);
    /* Decoding of the SWB parameters */
    DecodeSWBGenericParameters( st, lagIndices, nBands_search, BANDS, p2a_flags, hqswb_clas );


    /* Copy WB synthesis for SWB decoding */
    mvr2r( m_core, m, swb_lowband + swb_highband );

    /* Generic subband processing */
    DecodeSWBSubbands(
        st, m, swb_lowband, swb_highband, nBands, wBands, lagIndices,
        lagGains, BANDS, band_start, band_end, band_energy, p2a_flags, hqswb_clas, har_bands, subband_search_offset,
        prev_frm_hfe2, prev_stab_hfe2, band_width,subband_offsets, y2_ni, ni_seed);

    /* Smoothen highest frequencies */
    m[swb_lowband+swb_highband-1] *= 0.0625f;
    m[swb_lowband+swb_highband-2] *= 0.125f;
    m[swb_lowband+swb_highband-3] *= 0.25f;
    m[swb_lowband+swb_highband-4] *= 0.5f;

    /* Set frequencies below 6.4 kHz to zero */
    if( hqswb_clas == HQ_NORMAL )
    {
        for (k=0; k < swb_lowband; k++)
        {
            m[k] = 0.0f;
        }
    }

    return;
}

