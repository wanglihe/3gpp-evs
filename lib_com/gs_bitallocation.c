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

static float Find_bit_frac( const short nb_band, const short remaining_bits );

/*-------------------------------------------------------------------*
 * bands_and_bit_alloc()
 *
 * AC mode (GSC) bands and bits allocation
 *-------------------------------------------------------------------*/

void bands_and_bit_alloc(
    const short cor_strong_limit, /* i  : HF correlation                                        */
    const short noise_lev,        /* i  : dwn scaling factor                                    */
    const long  core_brate,       /* i  : core bit rate                                         */
    const short Diff_len,         /* i  : Lenght of the difference signal (before pure spectral)*/
    const short bits_used,        /* i  : Number of bit used before frequency Q                 */
    short *bit,             /* i/o: Number of bit allowed for frequency quantization      */
    float *Ener_per_bd_iQ,  /* i/o: Quantized energy vector                               */
    short *max_ener_band,   /* o  : Sorted order                                          */
    short *bits_per_bands_s,/* i/o: Number of bit allowed per allowed subband (Q3)        */
    short *nb_subbands,     /* o  : Number of subband allowed                             */
    const float *exc_diff,        /* i  : Difference signal to quantize (encoder side only)     */
    float *concat_in,       /* o  : Concatened PVQ's input vector (encoder side only)     */
    short *pvq_len,         /* o  : Number of bin covered with the PVQ                    */
    const short coder_type,       /* i  : coding type                                           */
    const short bwidth,           /* i  : input signal bandwidth                                */
    const short GSC_noisy_speech  /* i  : GSC noisy speech flag                                 */
)
{
    short bandoffset, i, j, nb_bands_max, bit_new_bands, bit_tmp, st_band, nb_bands;
    float bit_fracf, etmp;
    float sum_bit;
    float ener_vec[MBANDS_GN];
    short nb_tot_bands;
    short bit_index, bit_index_mem, imax;
    short pos, band;
    float SWB_bit_budget;
    float bits_per_bands[MBANDS_GN];

    /* initializations */
    nb_tot_bands = 16;
    set_f( bits_per_bands, 0.0f, MBANDS_GN );

    /* To adapt current energy band to PVQ freq band for sorting*/
    ener_vec[0] = Ener_per_bd_iQ[0] + Ener_per_bd_iQ[1];
    mvr2r( Ener_per_bd_iQ + 1, ener_vec, 15 );
    ener_vec[15] = ener_vec[14];

    for(i = 0; i < MBANDS_GN; i++)
    {
        ener_vec[i] = (float)((short)(ener_vec[i]*4096.f+0.5f));
    }

    /*------------------------------------------------------------------------
     * Determination of the number of bits available to the frequency domain
     * Allocation of a maximum number of band to be encoded
     *-----------------------------------------------------------------------*/

    nb_bands_max = nb_tot_bands;
    bit_new_bands = 5;

    bit_index = BRATE2IDX(core_brate)*17;
    bit_index_mem = bit_index;

    if( (coder_type == AUDIO || coder_type == INACTIVE) && bwidth == NB )
    {
        if(core_brate >= ACELP_9k60)
        {
            *bit = (short)(core_brate*(1.0f/50) + 0.5f) - bits_used - 25;
        }
        else
        {
            *bit = (short)(core_brate*(1.0f/50) + 0.5f) - bits_used - 21;
        }

        nb_tot_bands = 10;
    }
    else
    {
        *bit = (short)(core_brate*(1.0f/50) + 0.5f) - bits_used - GSC_freq_bits[bit_index];
    }

    if( GSC_noisy_speech )
    {
        SWB_bit_budget = *bit;
        nb_bands = 5;
        st_band = nb_bands;

        set_f( bits_per_bands, 0, MBANDS_GN );

        bit_fracf = Find_bit_frac(nb_bands, SWB_bit_budget);  /* Supplementary bits distributed only on first bands */

        nb_tot_bands = nb_bands_max - 6;

        if( nb_tot_bands > 16 )
        {
            nb_tot_bands = 16;
        }

        for(j = 0; j < 2; j++)
        {
            i = j;
            max_ener_band[j] = i;
            ener_vec[i] = 0;
        }

        for(; j < nb_bands; j++)
        {
            i = maximum(ener_vec, nb_tot_bands, &etmp);
            max_ener_band[j] = i;
            ener_vec[i] = 0;
        }

        set_f( bits_per_bands, bit_fracf, nb_bands );
    }
    else
    {
        bit_index++;
        bit_tmp = *bit - GSC_freq_bits[bit_index];
        bit_index++;
        nb_bands_max += GSC_freq_bits[bit_index];
        bit_index++;

        *pvq_len = 112;
        st_band = 7;

        if( core_brate <= ACELP_9k60 )
        {
            *pvq_len = 80;
            st_band = 5;
            if( Diff_len == 0 )
            {
                nb_bands_max += 2;
                bit_tmp -= 13;
            }
        }
        else if( Diff_len == 0 )
        {
            nb_bands_max += 2;
            bit_tmp -= 17;
        }

        nb_bands = *pvq_len/16;

        /*------------------------------------------------------------------------
         * Adjustement of the maximum number of bands in function of the
         * dynamics of the spectrum (more or less speech like)
         *-----------------------------------------------------------------------*/

        if( coder_type == INACTIVE || noise_lev >= NOISE_LEVEL_SP3 )
        {
            /* Probably classification error -> concentrate bits on LF */
            if( core_brate >= ACELP_8k00 )
            {
                nb_bands_max = nb_bands+1;
            }
            else
            {
                nb_bands_max = nb_bands;
            }
        }
        else if( noise_lev >= NOISE_LEVEL_SP2 ||
                 (core_brate <= ACELP_13k20 && core_brate >= ACELP_9k60 && cor_strong_limit == 0) )  /* Very low dynamic, tend to speech, do not try to code HF at all */
        {
            nb_bands_max -= 2;
        }
        else if( noise_lev >= NOISE_LEVEL_SP1) /* Very low dynamic, tend to speech, code less HF */
        {
            nb_bands_max -= 1;
        }

        if( bwidth == NB && nb_bands_max > 10 )
        {
            nb_bands_max = 10;
        }

        /*------------------------------------------------------------------------
         * Find extra number of band to code according to bit rate availables
         *-----------------------------------------------------------------------*/

        while ( bit_tmp >= bit_new_bands && nb_bands <= nb_bands_max - 1 )
        {
            bit_tmp -= bit_new_bands;
            nb_bands++;
        }

        /*------------------------------------------------------------------------
         * Fractional bits to distribute on the first x bands
         *-----------------------------------------------------------------------*/

        bit_fracf = Find_bit_frac(st_band, bit_tmp);  /* Supplementary bits distributed only on first bands */

        /*------------------------------------------------------------------------
         * Complete the bit allocation per frequency band
         *-----------------------------------------------------------------------*/

        imax = 5;
        if( core_brate > ACELP_9k60 )
        {
            imax = 7;
        }

        for(i = 0; i < imax; i++)
        {
            bits_per_bands[i] = GSC_freq_bits[bit_index] + bit_fracf;
            bit_index++;
        }

        if( Diff_len == 0 )
        {
            bit_index = bit_index_mem+10;
            for( i = 0; i < 7; i++ )
            {
                bits_per_bands[i] += GSC_freq_bits[bit_index];
                bit_index++;
            }
        }

        /*--------------------------------------------------------------------------
         * Complete the bit allocation per frequency band for 16kHz high brate mode
         *--------------------------------------------------------------------------*/

        for( j = st_band; j < nb_bands; j++ )
        {
            bits_per_bands[j] = bit_new_bands;
        }

        /*--------------------------------------------------------------------------
         * Compute a maximum band (band offset) for the search on maximal energy
         * This is function of the spectral dynamic and the bitrate
         *--------------------------------------------------------------------------*/

        bandoffset = nb_tot_bands - (nb_bands + 2);

        if( noise_lev <= NOISE_LEVEL_SP1a )
        {
            bandoffset--;
        }
        else if ( (core_brate <= ACELP_13k20 && (coder_type == INACTIVE || noise_lev >= NOISE_LEVEL_SP3)) ||
                  (core_brate <= ACELP_13k20 && core_brate >= ACELP_9k60 && cor_strong_limit == 0) )
        {
            bandoffset++;
        }

        if( bandoffset < 0 )
        {
            bandoffset = 0;
        }

        /*--------------------------------------------------------------------------
         * Initiazed sorted vector
         * For the first x bands to be included in th final sorted vector
         * Sort the remaining bands in decrease energy order
         *--------------------------------------------------------------------------*/

        for(j = 0; j < nb_tot_bands; j++)
        {
            max_ener_band[j] = -10;
        }

        for(j = 0; j < st_band; j++)
        {
            max_ener_band[j] = j;
            ener_vec[j] = -10;
        }

        pos = st_band;
        for(; j < nb_bands; j++)
        {
            i = maximum(ener_vec, nb_tot_bands-bandoffset, &etmp);
            if(i > pos)
            {
                pos = i;
            }
            max_ener_band[j] = i;
            ener_vec[i] = -10;
        }

        /* re-allocate bits to the frames such that the highest band with allocated bits is higher than the threshold */
        if( nb_tot_bands - bandoffset > nb_bands && ( pos > 7 && core_brate == ACELP_8k00 ) && bwidth == WB )
        {
            band = nb_tot_bands-bandoffset-nb_bands;
            for( j=0; j<band; j++ )
            {
                i = maximum( ener_vec, nb_tot_bands-bandoffset, &etmp );
                max_ener_band[nb_bands+j] = i;
                ener_vec[i] = -10;
                bits_per_bands[nb_bands+j] = 5;
            }
            nb_bands += band;

            bit_tmp = 5*band;
            if( band <= 2 )
            {
                for(j = st_band-1; j < nb_bands; j++)
                {
                    bits_per_bands[j] += 1;
                }
                bit_tmp += nb_bands - st_band + 1;
            }

            i = 0;
            j = 0;
            while( bit_tmp > 0 )
            {
                bits_per_bands[j] -= 1;

                if ( j == st_band - 1 - i )
                {
                    j = 0;
                }
                else
                {
                    ++j;
                }

                if( j == 0 && i < st_band - 1)
                {
                    i++;
                }

                bit_tmp -= 1;
            }
        }

        /*--------------------------------------------------------------------------
         * Bit sum verification for GSC inactive at very high rate
         * The maximum number of bits per band of length 16 is 112
         * Redistribute the overage bits if needed
         *--------------------------------------------------------------------------*/

        sum_bit = 0;
        j = 0;
        for( i = 0; i < nb_bands; i++ )
        {
            if( bits_per_bands[i] > 112 )
            {
                sum_bit += bits_per_bands[i] - 112;
                bits_per_bands[i] = 112;
                j = i+1;
            }

            /* safety check for overage bit reallocation */
            else if( bits_per_bands[i] + sum_bit/3 > 112 )
            {
                j = i+1;
            }
        }

        if( sum_bit != 0 )
        {
            sum_bit /= (nb_bands - j);
            for( i = j; i < nb_bands; i++ )
            {
                bits_per_bands[i] += sum_bit;
            }
        }
    }

    /*--------------------------------------------------------------------------
     * second step of bit sum verification, normally sum_bit == *bit
     *--------------------------------------------------------------------------*/

    sum_bit = 0.00f;
    for( i = 0; i < nb_bands; i++ )
    {
        bits_per_bands[i] = (float)floor(bits_per_bands[i]);
        sum_bit += bits_per_bands[i];
    }

    if( *bit > sum_bit )
    {
        i = nb_bands-1;
        while(*bit > sum_bit)
        {
            bits_per_bands[i]++;
            sum_bit++;
            i--;
            if(i==0)
            {
                i = nb_bands-1;
            }
        }
    }

    /*--------------------------------------------------------------------------
     * Recompute the real number/length of frequency bands to encode
     *--------------------------------------------------------------------------*/

    *nb_subbands = nb_bands;
    *pvq_len = *nb_subbands*16;

    /*--------------------------------------------------------------------------
     * Concatenate bands (encoder only)
     *--------------------------------------------------------------------------*/

    if( exc_diff != NULL )
    {
        for( j = 0; j < nb_bands; j++ )
        {
            mvr2r( exc_diff + max_ener_band[j]*16, concat_in+j*16, 16 );
        }
    }
    for( j = 0; j < nb_bands; j++ )
    {
        bits_per_bands_s[j] = ((short)bits_per_bands[j]) << 3;
    }

    return;
}


/*-------------------------------------------------------------------*
 * Find_bit_frac()
 *
 * Computes the fraction of the remaining bit budget to allocate to the bands
 *-------------------------------------------------------------------*/

static float Find_bit_frac(
    const short nb_band,
    const short remaining_bits
)
{
    float var_out;
    short inv_bandQ15;
    int L_num;

    inv_bandQ15 = 6553;
    if (nb_band == 7)
    {
        inv_bandQ15 = 4681;
    }
    L_num = inv_bandQ15*remaining_bits;
    L_num *= 8;
    var_out = L_num/262144.0f;

    return (var_out);
}
