/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_enc.h"
#include "rom_com.h"
#include "prot.h"


static void find_enr( const float data[], float band[], float *ptE, float *Etot, const short min_band,
                      const short max_band, float *Bin_E, const short bin_freq, float *band_ener );

/*-------------------------------------------------------------------*
 * analy_sp()
 *
 * Spectral analysis
 *-------------------------------------------------------------------*/

void analy_sp(
    float *speech,     /* i  : speech buffer                                    */
    float *Bin_E,      /* o  : per bin log energy spectrum                      */
    float *Bin_E_old,  /* o  : per bin log energy spectrum for mid-frame        */
    float *fr_bands,   /* o  : per band energy spectrum (2 analyses)            */
    float lf_E[],      /* o  : per bin E for first VOIC_BINS bins (without DC)  */
    float *Etot,       /* o  : total input energy                               */
    const short min_band,    /* i  : minimum critical band                            */
    const short max_band,    /* i  : maximum critical band                            */
    float *band_ener,  /* o: energy in critical frequency bands without minimum noise floor E_MIN */
    float *PS          /* o  : Per bin energy spectrum                          */
    ,float *fft_buff    /* o  : FFT coefficients                                 */
)
{
    short i_subfr, i;
    float *pt_bands, *pt_fft, *pt;
    const float *pt1;
    /*-----------------------------------------------------------------*
     * Perform two spectral analyses
     * Find energy per critical frequency band and total energy in dB
     *-----------------------------------------------------------------*/

    pt_bands = fr_bands;
    pt_fft = fft_buff;
    *Etot = 0.0f;

    for( i_subfr = 0; i_subfr <= 1; i_subfr++ )
    {
        /* set pointer to the beginning of the signal for spectral analysis */
        if (i_subfr == 0)
        {
            /* set the pointer for first analysis window */
            pt = speech + 3*(L_SUBFR/2) - L_FFT/2;
        }
        else
        {
            /* set the pointer for second analysis window */
            pt = speech + 7*(L_SUBFR/2) - L_FFT/2;
        }

        /* 1st half of the window */
        pt1 = sqrt_han_window;
        for( i=0; i<L_FFT/2; i++ )
        {
            pt_fft[i] = *pt++ **pt1++;
        }
        /* 2nd half of the window */
        for( i=L_FFT/2; i<L_FFT; i++ )
        {
            pt_fft[i] = *pt++ **pt1--;
        }

        /* compute the spectrum */
        fft_rel( pt_fft, L_FFT, LOG2_L_FFT );

        /* find energy per critical band */
        find_enr( pt_fft, pt_bands, lf_E + i_subfr * VOIC_BINS, Etot, min_band,
                  max_band, Bin_E + i_subfr * (L_FFT/2), (short)BIN, band_ener + i_subfr*NB_BANDS );

        pt_bands += NB_BANDS;
        pt_fft += L_FFT;
    }

    /* Average total log energy over both half-frames */
    *Etot = 10.0f * (float)log10( 0.5f * *Etot );

    /* Per-bin log-energy spectrum */
    Bin_E[L_FFT/2-1] = Bin_E[L_FFT/2-2];
    Bin_E[L_FFT-1] = Bin_E[L_FFT-2];

    for ( i=0; i<L_FFT/2; i++ )
    {
        Bin_E_old[i] = Bin_E[i];
        PS[i] = (Bin_E[i] + 1e-5f + Bin_E[i+L_FFT/2] + 1e-5f)/2.0f;
        Bin_E[i] = (float)(10.0f * log(PS[i]));
    }

    return;
}


/*------------------------------------------------------------------------*
 * find_enr()
 *
 * find input signal energy for each critical band and first 74 LF bins
 * The energy is normalized by the number of frequency bins in a channel
 *------------------------------------------------------------------------*/

static void find_enr(
    const float data[],    /* i  : fft result, for the format see fft_rel.c */
    float band[],    /* o  : per band energy                          */
    float *ptE,      /* o  : per bin energy  for low frequencies      */
    float *Etot,     /* i/o: total energy                             */
    const short min_band,  /* i  : minimum critical band                    */
    const short max_band,  /* i  : maximum critical band                    */
    float *Bin_E,    /* o  : Per bin energy                           */
    const short bin_freq,  /* i  : Number of frequency bins                 */
    float *band_ener /* o  : per band energy without E_MIN            */
)
{
    short i, cnt;
    float freq, tmp;
    const float *ptR, *ptI;
    short voic_band;
    float norm_val;

    norm_val = 4.0f / (L_FFT*L_FFT);

    voic_band = VOIC_BAND_8k;
    if ( bin_freq == 50 )
    {
        voic_band = VOIC_BAND;
    }

    ptR = &data[1];         /* first real  */
    ptI = &data[L_FFT-1];   /* first imaginary */

    /* for low frequency bins, save per bin energy for the use in find_tilt() */
    freq = bin_freq;
    for( i=0; i < voic_band; i++ )                      /* up to maximum allowed voiced critical band */
    {
        band[i] = 0.0f;
        cnt = 0;
        while( freq <= crit_bands[i] )
        {
            *ptE = *ptR **ptR + *ptI **ptI;           /* energy  */
            *ptE *= norm_val;                           /* normalization - corresponds to FFT normalization by 2/L_FFT */
            *Bin_E++ = *ptE;
            band[i] += *ptE++;
            ptR++;
            ptI--;

            freq += bin_freq;
            cnt++;
        }

        band[i] *= inv_tbl[cnt];                        /* normalization per frequency bin */

        band_ener[i] = band[i];                         /* per band energy without E_MIN   */

        if ( band[i] < E_MIN )
        {
            band[i] = E_MIN;
        }
    }

    /* continue computing the energy per critical band for higher frequencies */
    if ( bin_freq == 50 )
    {
        for( i = voic_band; i < NB_BANDS; i++ )
        {
            band[i] = 0.0f;
            cnt = 0;
            while( freq <= crit_bands[i] )
            {
                *Bin_E = *ptR **ptR + *ptI **ptI;
                *Bin_E *= norm_val;
                band[i] += *Bin_E;
                Bin_E++;
                ptR++;
                ptI--;

                freq += bin_freq;
                cnt++;
            }

            band[i] *= inv_tbl[cnt];                        /* normalization per frequency bin */

            band_ener[i] = band[i];                         /* per band energy without E_MIN   */

            if( band[i] < E_MIN )
            {
                band[i] = E_MIN;
            }
        }
    }

    /* find the total log energy */
    tmp = *Etot;
    for( i = min_band; i <= max_band; i++ )
    {
        tmp += band[i];
    }

    *Etot = tmp;

    return;
}
