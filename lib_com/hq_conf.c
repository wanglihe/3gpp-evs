/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"


/*--------------------------------------------------------------------------*
 * hq_configure()
 *
 * Configuration routine for HQ mode
 *--------------------------------------------------------------------------*/

void hq_configure(
    const short length,              /* i  : Frame length                      */
    const short hqswb_clas,          /* i  : HQ SWB class                      */
    const long core_brate,           /* i  : Codec bitrate                     */
    short *num_sfm,            /* o  : Total number of subbands          */
    short *nb_sfm,             /* o  : Total number of coded bands       */
    short *start_norm,         /* o  : First norm to be SDE encoded      */
    short *num_env_bands,      /* o  : Number coded envelope bands       */
    short *numnrmibits,        /* o  : Number of bits in fall-back norm encoding   */
    short *hq_generic_offset,  /* o  : Freq offset for HQ GENERIC        */
    short const **sfmsize,     /* o  : Subband bandwidths                */
    short const **sfm_start,   /* o  : Subband start coefficients        */
    short const **sfm_end      /* o  : Subband end coefficients          */
)
{
    *start_norm = 0;

    if( length == L_FRAME48k )
    {
        if ( hqswb_clas == HQ_GEN_FB )
        {
            *num_sfm = NB_SFM;
            *sfmsize = band_len;
            *sfm_start = band_start;
            *sfm_end = band_end;

            if( core_brate == HQ_32k )
            {
                *hq_generic_offset = HQ_GENERIC_FOFFSET_32K;
            }
            else if ( core_brate == HQ_16k40 || core_brate == HQ_24k40 )
            {
                *hq_generic_offset = HQ_GENERIC_FOFFSET_24K4;
            }

            /* setting start frequency of FD BWE */
            if( core_brate == HQ_32k )
            {
                *num_env_bands = SFM_N_STA_10k;
            }
            else if( core_brate == HQ_16k40 || core_brate == HQ_24k40 )
            {
                *num_env_bands = SFM_N_STA_8k;
            }

            *nb_sfm = *num_sfm;
        }
        else
        {
            if(hqswb_clas == HQ_HARMONIC)
            {
                *num_sfm = SFM_N_HARM_FB;
                *nb_sfm = SFM_N_HARM_FB;
                *num_env_bands = SFM_N_HARM_FB;

                *sfmsize = band_len_harm;
                *sfm_start = band_start_harm;
                *sfm_end = band_end_harm;
            }
            else if(hqswb_clas == HQ_HVQ)
            {
                if( core_brate == HQ_24k40 )
                {
                    *num_sfm = SFM_N_HARM_FB;
                    *nb_sfm = HVQ_THRES_SFM_24k;
                    *num_env_bands = *num_sfm - *nb_sfm;

                    *sfmsize = band_len_harm;
                    *sfm_start = band_start_harm;
                    *sfm_end = band_end_harm;
                    *start_norm = HVQ_THRES_SFM_24k;
                }
                else
                {
                    *num_sfm = SFM_N_HARM_FB;
                    *nb_sfm = HVQ_THRES_SFM_32k;
                    *num_env_bands = *num_sfm - *nb_sfm;

                    *sfmsize = band_len_harm;
                    *sfm_start = band_start_harm;
                    *start_norm = HVQ_THRES_SFM_32k;
                    *sfm_end = band_end_harm;
                }
            }
            else
            {
                *num_sfm = NB_SFM;
                *nb_sfm = *num_sfm;
                *num_env_bands = NB_SFM;

                *sfmsize = band_len;
                *sfm_start = band_start;
                *sfm_end = band_end;
            }
        }
    }
    else if( length == L_FRAME32k )
    {
        if( hqswb_clas == HQ_HARMONIC )
        {
            *num_sfm = SFM_N_HARM;
            *nb_sfm = SFM_N_HARM;
            *num_env_bands = SFM_N_HARM;

            *sfmsize = band_len_harm;
            *sfm_start = band_start_harm;
            *sfm_end = band_end_harm;
        }
        else if ( hqswb_clas == HQ_HVQ )
        {
            if( core_brate == HQ_24k40 )
            {
                *num_sfm = SFM_N_HARM;
                *nb_sfm = HVQ_THRES_SFM_24k;
                *num_env_bands = *num_sfm - *nb_sfm;

                *sfmsize = band_len_harm;
                *sfm_start = band_start_harm;
                *sfm_end = band_end_harm;
                *start_norm = HVQ_THRES_SFM_24k;
            }
            else
            {
                *num_sfm = SFM_N_HARM;
                *nb_sfm = HVQ_THRES_SFM_32k;
                *num_env_bands = *num_sfm - *nb_sfm;

                *sfmsize = band_len_harm;
                *sfm_start = band_start_harm;
                *start_norm = HVQ_THRES_SFM_32k;
                *sfm_end = band_end_harm;
            }
        }
        else if ( hqswb_clas == HQ_GEN_SWB )
        {
            *num_sfm = SFM_N_SWB;
            *sfmsize = band_len;
            *sfm_start = band_start;
            *sfm_end = band_end;

            if( core_brate == HQ_32k )
            {
                *hq_generic_offset = HQ_GENERIC_FOFFSET_32K;
            }
            else if ( core_brate == HQ_24k40 )
            {
                *hq_generic_offset = HQ_GENERIC_FOFFSET_24K4;
            }

            /* setting start frequency of FD BWE */
            if( core_brate == HQ_32k )
            {
                *num_env_bands = SFM_N_STA_10k;
            }
            else if( core_brate == HQ_24k40 )
            {
                *num_env_bands = SFM_N_STA_8k;
            }

            *nb_sfm = *num_sfm;
        }
        else
        {
            /* HQ_NORMAL and HQ_TRANSIENT */
            *num_sfm = SFM_N_SWB;
            *nb_sfm = *num_sfm;
            *num_env_bands = SFM_N_SWB;

            *sfmsize = band_len;
            *sfm_start = band_start;
            *sfm_end = band_end;
        }
    }
    else
    {
        *num_sfm = SFM_N_WB;
        *nb_sfm = *num_sfm;
        *num_env_bands = SFM_N_WB;

        *sfmsize = band_len_wb;
        *sfm_start = band_start_wb;
        *sfm_end = band_end_wb;
    }

    *numnrmibits = (*num_env_bands - 1) * NORMI_BITS;

    return;
}
