/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"
#include "basop_util.h"
#include "basop_proto_func.h"

/*--------------------------------------------------------------------------*
 * mdct_spectrum_denorm()
 *
 *
 *--------------------------------------------------------------------------*/

void mdct_spectrum_denorm(
    const int   inp_vector[],
    float y2[],
    const short band_start[],
    const short band_end[],
    const short band_width[],
    const float band_energy[],
    const int   npulses[],
    const short bands,
    const float ld_slope,
    const float pd_thresh
)
{
    short i, k;
    float Eyy, gamma, pd, gain_tweak;

    for (k = 0; k < bands; k++)
    {
        Eyy = 0;
        for (i = band_start[k]; i <= band_end[k]; i++)
        {
            Eyy += (float) inp_vector[i] * inp_vector[i];
        }

        if (Eyy > 0.0f)
        {
            /* Set gamma to be pulse gain which results in perfect quantized subband energy */
            gamma = (float) sqrt (pow (2.0f, band_energy[k]) / Eyy);

            /* Adjust gamma based on pulse density (0 bit MSE gain estimator) */
            pd = (float) npulses[k] / band_width[k];
            if (pd < pd_thresh)
            {
                gain_tweak = (float) pow (2.0f, (ld_slope * log2_f (pd / pd_thresh)));
                gamma *= gain_tweak;
            }

            for (i = band_start[k]; i <= band_end[k]; i++)
            {
                y2[i] = gamma * inp_vector[i];
            }
        }
    }

    return;
}
/*--------------------------------------------------------------------------*
 * hq2_core_configure()
 *
 *
 *--------------------------------------------------------------------------*/

void hq2_core_configure(
    const short frame_length,               /* i  : frame length                            */
    const short num_bits,                   /* i  : number of bits                          */
    const short is_transient,               /* i  : transient flag                          */
    short *bands,
    short *length,
    short band_width[],
    short band_start[],
    short band_end[],
    Word32 *L_qint,
    Word16 *eref_fx,
    Word16 *bit_alloc_weight_fx,
    short *gqlevs,
    short *Ngq,
    short *p2a_bands,
    float *p2a_th,
    float *pd_thresh,
    float *ld_slope,
    float *ni_coef,
    float *ni_pd_th,
    long  bwe_br
)
{
    const Xcore_Config *xcore_config;

    short i, k;
    short bands_sh;

    if( frame_length == L_FRAME8k )
    {
        if( is_transient )
        {
            if (num_bits <= ACELP_7k20 / 50)
            {
                xcore_config = &xcore_config_8kHz_007200bps_short;
            }
            else if (num_bits <= ACELP_8k00 / 50)
            {
                xcore_config = &xcore_config_8kHz_008000bps_short;
            }
            else if (num_bits <= ACELP_13k20 / 50)
            {
                xcore_config = &xcore_config_8kHz_013200bps_short;
            }
            else
            {
                xcore_config = &xcore_config_8kHz_016400bps_short;
            }
        }
        else
        {
            if (num_bits <= ACELP_7k20 / 50)
            {
                xcore_config = &xcore_config_8kHz_007200bps_long;
            }
            else if (num_bits <= ACELP_8k00 / 50)
            {
                xcore_config = &xcore_config_8kHz_008000bps_long;
            }
            else if (num_bits <= ACELP_13k20 / 50)
            {
                xcore_config = &xcore_config_8kHz_013200bps_long;
            }
            else
            {
                xcore_config = &xcore_config_8kHz_016400bps_long;
            }
        }
    }
    else if( frame_length == L_FRAME16k )
    {
        if( is_transient )
        {
            if (num_bits <= ACELP_13k20 / 50)
            {
                xcore_config = &xcore_config_16kHz_013200bps_short;
            }
            else
            {
                xcore_config = &xcore_config_16kHz_016400bps_short;
            }
        }
        else
        {
            if (num_bits <= ACELP_13k20 / 50)
            {
                xcore_config = &xcore_config_16kHz_013200bps_long;
            }
            else
            {
                xcore_config = &xcore_config_16kHz_016400bps_long;
            }
        }
    }
    else    /* (bwidth == SWB) */
    {
        if( is_transient )
        {
            if (bwe_br == ACELP_13k20)
            {
                xcore_config = &xcore_config_32kHz_013200bps_short;
            }
            else
            {
                xcore_config = &xcore_config_32kHz_016400bps_short;
            }
        }
        else
        {
            if (bwe_br == ACELP_13k20)
            {
                xcore_config = &xcore_config_32kHz_013200bps_long;
            }
            else
            {
                xcore_config = &xcore_config_32kHz_016400bps_long;
            }
        }
    }

    *bands = xcore_config->bands;
    *length = xcore_config->bw;
    *L_qint = xcore_config->L_qint;
    *eref_fx = xcore_config->eref_fx;
    *bit_alloc_weight_fx = xcore_config->bit_alloc_weight_fx;
    *gqlevs = xcore_config->gqlevs;
    *Ngq = xcore_config->Ngq;

    *p2a_bands = xcore_config->p2a_bands;
    *p2a_th = xcore_config->p2a_th;

    *pd_thresh = xcore_config->pd_thresh;
    *ld_slope = xcore_config->ld_slope;
    *ni_coef = xcore_config->ni_coef;
    *ni_pd_th = xcore_config->ni_pd_th;

    mvs2s (xcore_config->band_width, band_width, *bands);

    /* Expand band_width[] table for short windows */
    if( is_transient )
    {
        bands_sh = *bands;
        *bands = NUM_TIME_SWITCHING_BLOCKS * bands_sh;
        *length *= NUM_TIME_SWITCHING_BLOCKS;

        for( i = 1; i <= 3; i++ )
        {
            for( k = 0; k < bands_sh; k++ )
            {
                band_width[i * bands_sh + k] = band_width[k];
            }
        }
    }

    /* Formulate band_start and band_end tables from band_width table */
    band_start[0] = 0;
    band_end[0] = band_width[0] - 1;
    for( k = 1; k < *bands; k++ )
    {
        band_start[k] = band_start[k - 1] + band_width[k - 1];
        band_end[k] = band_start[k] + band_width[k] - 1;
    }


    return;
}

/*--------------------------------------------------------------------------*
 * reverse_transient_frame_energies()
 *
 *
 *--------------------------------------------------------------------------*/

void reverse_transient_frame_energies(
    float band_energy[],          /* o  : band energies       */
    const short bands                   /* i  : number of bands     */
)
{
    short k, k1, k2;
    float be;

    k1 = bands/4;
    k2 = bands/2-1;
    for( k = 0; k < bands/8; k++ )
    {
        be = band_energy[k1];
        band_energy[k1] = band_energy[k2];
        band_energy[k2] = be;
        k1++, k2--;
    }

    k1 = 3*bands/4;
    k2 = bands-1;
    for( k = 0; k < bands/8; k++ )
    {
        be = band_energy[k1];
        band_energy[k1] = band_energy[k2];
        band_energy[k2] = be;
        k1++, k2--;
    }

    return;
}

void bit_allocation_second_fx(
    Word32 *Rk,
    Word32 *Rk_sort,
    Word16  BANDS,
    const Word16 *band_width,
    Word16 *k_sort,
    Word16 *k_num,
    const Word16 *p2a_flags,
    const Word16  p2a_bands,
    const Word16 *last_bitalloc,
    const Word16  input_frame
)
{
    Word16 k, k2 = 0;
    Word16 ever_bits[BANDS_MAX], ever_sort[BANDS_MAX];/*Q12 */
    Word16 class_flag = 0;
    Word16 rk_temp = 32767, ever_temp = 32767;/*Q12 */
    Word16 exp;
    Word16 tmp;
    Word32 L_tmp;

    for (k = 0; k < BANDS; k++)
    {
        if((( sub(k_sort[k],sub(BANDS,p2a_bands)) >= 0 )&&( sub(p2a_flags[k_sort[k]],1) == 0 )) ||
                (( sub(k_sort[k],sub(BANDS,2)) >= 0 )&&( sub(last_bitalloc[sub(k_sort[k], sub(BANDS,2))], 1) == 0 )))
        {
            exp = norm_s(band_width[k_sort[k]]);
            tmp = shl(band_width[k_sort[k]],exp);/*Q(exp) */
            tmp = div_s(16384,tmp);/*Q(15+14-exp = 29-exp) */
            L_tmp = Mult_32_16(Rk_sort[k],tmp);/* Q(16+29-exp-15 = 30-exp) */
            tmp = sub(18,exp);
            ever_bits[k] = extract_l(L_shr(L_tmp,tmp));/*Q12 */
            if( sub(ever_bits[k],rk_temp) < 0 )
            {
                rk_temp = ever_bits[k];
                k2 = k;
            }
            class_flag = 1;
        }
    }
    if( class_flag ==0 || sub(input_frame,L_FRAME8k) == 0)
    {
        for(k = 0; k < BANDS; k++)
        {
            if( sub(k_sort[k],sub(BANDS,p2a_bands)) < 0 && Rk_sort[k] > 0 )
            {
                exp = norm_s(band_width[k_sort[k]]);
                tmp = shl(band_width[k_sort[k]],exp);/*Q(exp) */
                tmp = div_s(16384,tmp);/*Q(15+14-exp = 29-exp) */
                L_tmp = Mult_32_16(Rk_sort[k],tmp);/* Q(16+29-exp-15 = 30-exp) */
                tmp = sub(18,exp);
                ever_sort[k] = extract_l(L_shr(L_tmp,tmp));/*Q12 */
                IF(sub(ever_sort[k],ever_temp) < 0)
                {
                    ever_temp = ever_sort[k];
                    k2 = k;
                }
            }
        }
    }

    k_num[0] = k2;
    if(sub(k_sort[k2],sub(BANDS,1)) == 0)
    {
        for (k = 0; k < BANDS; k++)
        {
            if(sub(k_sort[k],sub(k_sort[k2],1)) == 0)
            {
                k_num[1] = k;
            }
        }
    }
    else if(k_sort[k2] == 0)
    {
        for (k = 0; k < BANDS; k++)
        {
            if(sub(k_sort[k],add(k_sort[k2],1)) == 0)
            {
                k_num[1] = k;
            }
        }
    }
    else
    {
        if ( L_sub( Rk[sub(k_sort[k2],1)],Rk[add(k_sort[k2],1)] ) < 0 )
        {
            for (k = 0; k < BANDS; k++)
            {
                if(sub(k_sort[k],sub(k_sort[k2],1)) == 0)
                {
                    k_num[1] = k;
                }
            }
        }
        else
        {
            for (k = 0; k < BANDS; k++)
            {
                if(sub(k_sort[k],add(k_sort[k2],1)) == 0)
                {
                    k_num[1] = k;
                }
            }
        }
    }
}

/*--------------------------------------------------------------------------*
 * spt_shorten_domain_pre()
 *
 * Compute shorten subband if previous frame has spectral peak.
 *--------------------------------------------------------------------------*/

void spt_shorten_domain_pre(
    const short band_start[],           /* i:   Starting position of sub band             */
    const short band_end[],             /* i:   End position of sub band                  */
    const short prev_SWB_peak_pos[],    /* i:   Spectral peak                             */
    const short BANDS,                  /* i:   total number of bands                     */
    const long  bwe_br,                 /* i:   bitrate information                       */
    short       new_band_start[],       /* o:   Starting position of new shorten sub band */
    short       new_band_end[],         /* o:   End position of new shorten sub band      */
    short       new_band_width[]        /* o:   new sub band bandwidth                    */
)
{
    int j;
    int k;
    int kpos;

    short new_band_width_half;
    const short *p_bw_SPT_tbl; /* pointer of bw_SPT_tbl */

    p_bw_SPT_tbl = bw_SPT_tbl[0];
    if( bwe_br == HQ_16k40 )
    {
        p_bw_SPT_tbl = bw_SPT_tbl[1];
    }

    kpos = 0;
    j = 0;
    for(k=BANDS-SPT_SHORTEN_SBNUM; k<BANDS; k++)
    {
        if ( prev_SWB_peak_pos[kpos] != 0)
        {
            new_band_width[j] = p_bw_SPT_tbl[j];

            /*shorten the bandwidth for pulse resolution*/
            new_band_width_half = new_band_width[j]/2;
            new_band_start[j] = prev_SWB_peak_pos[kpos] - new_band_width_half;
            new_band_end[j]   = prev_SWB_peak_pos[kpos] + new_band_width_half;

            if( new_band_start[j] < band_start[k] )
            {
                new_band_start[j] = band_start[k];
                new_band_end[j] = new_band_start[j] + (new_band_width[j] - 1);
            }
            else if ( new_band_end[j] > band_end[k] )
            {
                new_band_end[j]   = band_end[k];
                new_band_start[j] = new_band_end[j] - (new_band_width[j] - 1);
            }
        }
        else
        {
            new_band_width[j] = p_bw_SPT_tbl[j];

            /*shorten the bandwidth for pulse resolution*/
            new_band_width_half = new_band_width[j]/2;
            new_band_start[j] = (band_start[k]+band_end[k])/2 - new_band_width_half;
            new_band_end[j]   = (band_start[k]+band_end[k])/2 + new_band_width_half;
        }

        kpos++;
        j++;
    }

    return;
}

/*--------------------------------------------------------------------------*
 * spt_shorten_domain_band_save()
 *
 * Store the original subband information
 *--------------------------------------------------------------------------*/

void spt_shorten_domain_band_save(
    const short bands,                  /* i:   total subband                */
    const short band_start[],           /* i:   starting position of subband */
    const short band_end[],             /* i:   end position of subband      */
    const short band_width[],           /* i:   band width of subband        */
    short org_band_start[],       /* o:   starting position of subband */
    short org_band_end[],         /* o:   end position of subband      */
    short org_band_width[]        /* o:   band width of subband        */
)
{
    int k;
    int kpos;

    kpos = 0;
    for(k=bands-SPT_SHORTEN_SBNUM; k<bands; k++)
    {
        org_band_start[kpos] = band_start[k];
        org_band_end[kpos]   = band_end[k];
        org_band_width[kpos] = band_width[k];
        kpos++;
    }

    return;
}

/*--------------------------------------------------------------------------*
 * spt_shorten_domain_band_restore()
 *
 * Restrore the subband information
 *--------------------------------------------------------------------------*/

void spt_shorten_domain_band_restore(
    const short bands,                  /* i:   total subband                */
    short band_start[],           /* i/o: starting position of subband */
    short band_end[],             /* i/o: end position of subband      */
    short band_width[],           /* i/o: band width of subband        */
    const short org_band_start[],       /* o:   starting position of subband */
    const short org_band_end[],         /* o:   end position of subband      */
    const short org_band_width[]        /* o:   band width of subband        */
)
{
    int k;
    int kpos;

    kpos = 0;
    for(k=bands-SPT_SHORTEN_SBNUM; k<bands; k++)
    {
        band_start[k] = org_band_start[kpos];
        band_end[k]   = org_band_end[kpos];
        band_width[k] = org_band_width[kpos];
        kpos++;
    }

    return;
}

/*--------------------------------------------------------------------------*
 * spt_swb_peakpos_tmp_save
 *
 * Save Peak position for every higher subband
 *--------------------------------------------------------------------------*/

void spt_swb_peakpos_tmp_save(
    const float y2[],                   /* i:   coded spectral information   */
    const short bands,                  /* i:   total number of bands        */
    const short band_start[],           /* i:   starting position of subband */
    const short band_end[],             /* i:   end position of subband      */
    short prev_SWB_peak_pos_tmp[] /* o:   spectral peaks               */
)
{

    short   i, j, k;
    float peak_max;

    j = 0;
    for(k=bands-SPT_SHORTEN_SBNUM; k<bands; k++)
    {
        peak_max = 0;
        prev_SWB_peak_pos_tmp[j] = 0;
        for(i=band_start[k]; i<=band_end[k]; i++)
        {
            if( peak_max < fabs(y2[i]) )
            {
                peak_max =(float) fabs(y2[i]);
                prev_SWB_peak_pos_tmp[j] = i;
            }
        }
        j++;
    }
    return;
}
