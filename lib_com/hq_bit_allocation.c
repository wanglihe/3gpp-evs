/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"


/*--------------------------------------------------------------------------*
 * hq_bit_allocation()
 *
 * Assign bits for HQ fine structure coding with PVQ
 *--------------------------------------------------------------------------*/

void hq_bit_allocation(
    const long core_brate,      /* i  : Core bit-rate                      */
    const short length,         /* i  : Frame length                       */
    const short hqswb_clas,     /* i  : HQ class                           */
    short *num_bits,      /* i/o: Remaining bit budget               */
    const short *normqlg2,      /* i  : Quantized norms                    */
    const short nb_sfm,         /* i  : Number sub bands to be encoded     */
    const short *sfmsize,       /* i  : Sub band bandwidths                */
    float *noise_level,   /* o  : HVQ noise level                    */
    short *R,             /* o  : Bit allocation per sub band        */
    short *Rsubband,      /* o  : Fractional bit allocation (Q3)     */
    short *sum,           /* o  : Sum of allocated shape bits        */
    short *core_sfm,      /* o  : Last coded band in core            */
    const short num_env_bands   /* i  : Number sub bands to be encoded for HQ_SWB_BWE  */
)
{
    short sfm_limit = nb_sfm;
    short i;
    short idx[NB_SFM];
    short wnorm[NB_SFM];
    short avrg_wnorm;
    short E_low;
    short E_hb_mean;
    short E_max;
    short i_max;


    set_s( R, 0, NB_SFM );
    for( i = 0; i < nb_sfm; i++ )
    {
        idx[i] = i;
    }
    if( hqswb_clas != HQ_TRANSIENT && hqswb_clas != HQ_HVQ && !(length == L_FRAME16k && core_brate == HQ_32k))
    {
        /* 'nf_idx' 2-bits index written later */
        *num_bits -= 2;
    }

    if ( hqswb_clas == HQ_GEN_SWB || hqswb_clas == HQ_GEN_FB )
    {
        if (core_brate == HQ_32k)
        {
            *num_bits -= HQ_GENERIC_SWB_NBITS2;
        }
        else
        {
            *num_bits -= HQ_GENERIC_SWB_NBITS;
        }
        if (length == L_FRAME48k)
        {
            *num_bits -= HQ_GENERIC_FB_NBITS;
        }
    }

    if( (length == L_FRAME48k)  &&  (hqswb_clas != HQ_HARMONIC)  &&  (hqswb_clas != HQ_HVQ) )
    {
        map_quant_weight( normqlg2, wnorm, hqswb_clas == HQ_TRANSIENT );
    }
    else
    {
        mvs2s( normqlg2, wnorm, nb_sfm );
    }


    if( hqswb_clas == HQ_HARMONIC )
    {
        /* classification and limit bandwidth for bit allocation */
        sfm_limit -= 2;
        limit_band_noise_level_calc( wnorm, &sfm_limit, core_brate, noise_level );

        /* Detect important band in high frequency region */
        E_low = sum_s(wnorm, SFM_G1);
        i_max = 0;
        E_max = MIN16B;
        E_hb_mean = 0;
        for( i = SFM_G1; i < nb_sfm; i++)
        {
            E_hb_mean += wnorm[i];
            if( wnorm[i] > E_max)
            {
                E_max = wnorm[i];
                i_max = i;
            }
        }
        E_hb_mean = E_hb_mean >> 4; /* Truncated division by SFM_G1 */
        if ( E_max * 15 < E_low || E_max * 0.67 < E_hb_mean || i_max < sfm_limit )
        {
            i_max = 0;
        }

        set_s( wnorm + sfm_limit, -20, nb_sfm - sfm_limit );
        if ( i_max > 0)
        {
            wnorm[i_max] = E_max;
        }

    }


    if( hqswb_clas == HQ_HVQ )
    {
        *sum = 0;
    }
    else if ( hqswb_clas == HQ_GEN_SWB || (hqswb_clas == HQ_TRANSIENT && length == L_FRAME32k && core_brate <= HQ_32k) )
    {
        *sum = BitAllocF( wnorm, core_brate, *num_bits, nb_sfm, R, Rsubband, hqswb_clas, num_env_bands );
    }
    else if( length == L_FRAME16k && core_brate == HQ_32k )
    {
        if( hqswb_clas != HQ_TRANSIENT )
        {
            avrg_wnorm = wnorm[10];
            for( i=11; i<18; i++ )
            {
                avrg_wnorm += wnorm[i];
            }

            avrg_wnorm /= 8;
            for( i=0; i<4; i++ )
            {
                if( wnorm[i] < avrg_wnorm )
                {
                    wnorm[i] = avrg_wnorm;
                }
            }

            /* Estimate number of bits per band */
            *sum = BitAllocWB( wnorm, *num_bits, nb_sfm, R, Rsubband );
        }
        else
        {
            reordvct(wnorm, nb_sfm, idx);
            bitalloc( wnorm, idx, *num_bits, nb_sfm, QBIT_MAX2, R, sfmsize, hqswb_clas );
            bitallocsum( R, nb_sfm, sum, Rsubband, *num_bits, length, sfmsize );
        }
    }
    else
    {
        reordvct( wnorm, nb_sfm, idx );

        /* enlarge the wnorm value so that more bits can be allocated to (sfm_limit/2 ~ sfm_limit) range */
        if( hqswb_clas == HQ_HARMONIC )
        {
            for( i=sfm_limit/2; i<sfm_limit; i++ )
            {
                wnorm[i] = wnorm[sfm_limit/2 - 1];
            }
        }

        bitalloc( wnorm, idx, *num_bits, nb_sfm, QBIT_MAX2, R, sfmsize, hqswb_clas );
        bitallocsum( R, nb_sfm, sum, Rsubband, *num_bits, length, sfmsize );
    }

    /* Find last coded core band */
    *core_sfm = nb_sfm - 1;
    if( hqswb_clas == HQ_NORMAL || hqswb_clas == HQ_HARMONIC )
    {
        *core_sfm = find_last_band( R, nb_sfm );
    }
    else if (hqswb_clas == HQ_GEN_SWB || hqswb_clas == HQ_GEN_FB )
    {
        *core_sfm = find_last_band( R, nb_sfm );
        if (*core_sfm < num_env_bands)
        {
            *core_sfm = num_env_bands-1;
        }
    }

    *num_bits -= *sum;

    return;
}
