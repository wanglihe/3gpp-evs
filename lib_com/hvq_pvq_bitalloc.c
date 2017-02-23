/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"

/*--------------------------------------------------------------------------*/
/*  Function  hvq_pvq_bitalloc                                              */
/*  ~~~~~~~~~~~~~~~~~~~~~~~~                                                */
/*                                                                          */
/*  Calculate the number of PVQ bands to code and allocate bits based on    */
/*  the number of available bits.                                           */
/*--------------------------------------------------------------------------*/

short hvq_pvq_bitalloc(
    short num_bits,               /* i/o: Number of available bits (including gain bits) */
    const short brate,                  /* i  : bitrate                     */
    const short bwidth,                 /* i  : Encoded bandwidth           */
    const short *ynrm,                  /* i  : Envelope coefficients       */
    const int   manE_peak,              /* i  : Peak energy mantissa        */
    const short expE_peak,              /* i  : Peak energy exponent        */
    short *Rk,                    /* o  : bit allocation for concatenated vector */
    short *R,                     /* i/o: Global bit allocation       */
    short *sel_bands,             /* o  : Selected bands for encoding */
    short *n_sel_bands            /* o  : No. of selected bands for encoding */
)
{
    short num_bands, k, band_max_bits;
    short k_max;
    short k_start;
    float E_max;
    int   E_max5;
    short expo;
    short align;
    int acc;
    float tmp;
    int   env_mean;
    short reciprocal;
    short num_sfm;
    short bit_cost;

    if (bwidth == FB)
    {
        num_sfm = SFM_N_HARM_FB;
    }
    else
    {
        num_sfm = SFM_N_HARM;
    }

    if (brate == HQ_24k40)
    {
        band_max_bits = HVQ_BAND_MAX_BITS_24k;
        k_start = HVQ_THRES_SFM_24k;
        if (bwidth == FB)
        {
            reciprocal = 2731;
        }
        else
        {
            reciprocal = 3277;
        }
    }
    else
    {
        band_max_bits = HVQ_BAND_MAX_BITS_32k;
        k_start = HVQ_THRES_SFM_32k;
        if (bwidth == FB)
        {
            reciprocal = 3641;
        }
        else
        {
            reciprocal = 4681;
        }
    }

    num_bands = num_bits / band_max_bits;
    num_bits = num_bits - num_bands*band_max_bits;

    if (num_bits >= HVQ_NEW_BAND_BIT_THR)
    {
        num_bands++;
    }
    else
    {
        num_bits += band_max_bits;
    }

    /* safety check in case of bit errors */
    if (num_bands < 1)
    {
        return 0;
    }

    *n_sel_bands = 0;
    env_mean = 0;
    E_max = 0;
    k_max = k_start;
    for ( k = k_start; k < num_sfm; k++ )
    {
        tmp = dicn[ynrm[k]];
        env_mean += ynrm[k];
        if( tmp > E_max )
        {
            E_max = tmp;
            k_max = k;
        }
    }
    env_mean = 2*((int)env_mean * (int)reciprocal);

    if(band_len_harm[k_max] == 96)
    {
        bit_cost = 61;
    }
    else
    {
        QuantaPerDsDirac(band_len_harm[k_max], 1, hBitsN, &bit_cost);
    }


    expo = max(0, ynrm[k_max]-1) >> 1;
    E_max5 = E_max5_tbl[ynrm[k_max]];
    align = expo - expE_peak;
    align = align + (19 - 14) - (31 - 2*12);
    if (align < 0)
    {
        acc = E_max5 - (manE_peak >> -align);
    }
    else
    {
        acc = (E_max5 >> align) - manE_peak;
    }
    if ( acc > 0
            && ((env_mean - (ynrm[k_max]<<16)) > 0x30000L)
            && num_bands > 1
            && (num_bits - HVQ_PVQ_GAIN_BITS)<<3 >= bit_cost )
    {
        sel_bands[*n_sel_bands] = k_max;
        (*n_sel_bands)++;
        R[k_max] = 1; /* Mark that the band has been encoded for fill_spectrum */
    }

    /* Allocate bits */
    for (k = 0; k < num_bands-1; k++)
    {
        Rk[k] = band_max_bits-HVQ_PVQ_GAIN_BITS;
    }
    Rk[num_bands-1] = num_bits-HVQ_PVQ_GAIN_BITS;

    return num_bands;
}

