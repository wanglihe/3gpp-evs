/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <math.h>
#include "options.h"
#include "prot.h"
#include "rom_com.h"

/*--------------------------------------------------------------------------*
 * Local functions
 *--------------------------------------------------------------------------*/

static void overlap_hq_bwe( const float *hq_swb_overlap_buf, float *coeff_out1, const short n_swb_overlap_offset, const short n_swb_overlap,
                            const short *R, const short num_env_bands, const short num_sfm, const short *sfm_end );

/*--------------------------------------------------------------------------*
 * hq_swb_harmonic_calc_norm_envelop()
 *
 * Calculate normalization envelop
 *--------------------------------------------------------------------------*/

void hq_swb_harmonic_calc_norm_envelop(
    float *SWB_signal,           /* i  : input signal                */
    float *envelope,             /* o  : output envelope             */
    int   L_swb_norm,            /* i  : length of normaliztion      */
    int   SWB_flength            /* i  : length of input signal      */
)
{
    int lookback;
    int env_index;
    int n_freq;
    int n_lag_now;
    int n_lag;
    int i;

    lookback = L_swb_norm/2;
    env_index = 0;
    for (n_freq = 0; n_freq < lookback; n_freq++)
    {
        n_lag_now = lookback+n_freq;

        /* Apply MA filter */
        envelope[env_index] = EPSILON;
        for (n_lag=0; n_lag<n_lag_now; n_lag++)
        {
            envelope[env_index] += (float)fabs(SWB_signal[n_lag]);
        }
        env_index++;
    }

    n_lag_now = L_swb_norm;

    for (n_freq = 0; n_freq < SWB_flength-L_swb_norm; n_freq++)
    {
        /* Apply MA filter */
        envelope[env_index] = EPSILON;
        for (n_lag=0; n_lag<n_lag_now; n_lag++)
        {
            envelope[env_index] += (float)fabs(SWB_signal[n_freq+n_lag]);
        }
        env_index++;
    }

    for (n_freq = SWB_flength-L_swb_norm, i = 0; n_freq<SWB_flength-lookback; n_freq++, i++)
    {
        n_lag_now = L_swb_norm-i;

        /* Apply MA filter */
        envelope[env_index] = EPSILON;
        for (n_lag=0; n_lag<n_lag_now; n_lag++)
        {
            envelope[env_index] += (float)fabs(SWB_signal[n_freq+n_lag]);
        }
        env_index++;
    }

    return;
}

/*--------------------------------------------------------------------------*
 * noise_level_calc()
 *
 * Calculate noise level and limited band
 *--------------------------------------------------------------------------*/

void limit_band_noise_level_calc(
    short *wnorm,             /* i  : reordered norm of sub-vectors        */
    short *limit,             /* o  : highest band of bit allocation       */
    long core_brate,         /* o  : bit rate                             */
    float *noise_level        /* o  : noise level                          */
)
{
    float ener_limit, ener_sum, fact;
    short i;
    short nb_sfm;

    nb_sfm = *limit;
    ener_limit = 1e-5f;
    *noise_level = 0.0f;
    for(i=0; i< 10; i++)
    {
        ener_limit += wnorm[i];
        *noise_level += (float)abs(wnorm[i+1] - wnorm[i]);
    }
    ener_sum = ener_limit;

    for(i=10; i<(nb_sfm-1); i++)
    {
        ener_sum += wnorm[i];
        *noise_level += (float)abs(wnorm[i+1] - wnorm[i]);
    }
    ener_sum += wnorm[nb_sfm-1];

    if( core_brate == HQ_24k40 )
    {
        fact = 0.885f;
    }
    else
    {
        fact = 0.942f;
    }

    i = 9;
    while ( ener_limit < ener_sum * fact && i+1 < nb_sfm )
    {
        ener_limit += wnorm[++i];
    }
    *limit = i;

    /* calculate noise level for spectrum filling */
    *noise_level /= ener_sum;
    if(*noise_level < 0)
    {
        *noise_level = 0.25f;
    }

    *noise_level = 0.25f - *noise_level;

    if( *noise_level < 0.0f )
    {
        *noise_level = 0.0f;
    }

    return;
}

/*--------------------------------------------------------------------------*
 * build_nf_codebook()
 *
 * Build noise-fill codebook for HQ mode
 *--------------------------------------------------------------------------*/

short build_nf_codebook(         /* o  : Number of coefficients in nf codebook */
    const short flag_32K_env_ho, /* i  : Envelope attenuation hangover flag */
    const float *coeff,          /* i  : Coded spectral coefficients   */
    const short *sfm_start,      /* i  : Subband start indices         */
    const short *sfmsize,        /* i  : Subband widths                */
    const short *sfm_end,        /* i  : Subband end indices           */
    const short last_sfm,        /* i  : Last coded  band              */
    const short *R,              /* i  : Per-band bit allocation       */
    float *CodeBook,       /* o  : Noise-fill codebook           */
    float *CodeBook_mod    /* o  : Densified noise-fill codebook */
)
{
    short sfm_base;
    short sfm;
    short E_cb_vec;
    short i,j;
    short cb_size;

    /* Build codebook */
    cb_size = 0;
    for ( sfm = 0; sfm <= last_sfm; sfm++ )
    {
        if (R[sfm] != 0)
        {
            if (flag_32K_env_ho)
            {
                /* Build compressed (+/- 1) noise-fill codebook */
                sfm_base = sfm_start[sfm];
                for (i = 0; i < sfmsize[sfm]/8; i++)
                {
                    E_cb_vec = 0;
                    for (j = sfm_base+i*8; j < sfm_base+(i+1)*8; j++)
                    {
                        if (coeff[j] > 0.0f)
                        {
                            CodeBook_mod[cb_size] = 1.0f;
                            E_cb_vec++;
                        }
                        else if (coeff[j] < 0.0f)
                        {
                            CodeBook_mod[cb_size] = -1.0f;
                            E_cb_vec++;
                        }
                        else
                        {
                            CodeBook_mod[cb_size] = 0.0f;
                        }
                        cb_size++;
                    }

                    if (E_cb_vec < 2)
                    {
                        cb_size -= 8;
                    }
                }
            }
            else
            {
                for (j = sfm_start[sfm]; j < sfm_end[sfm]; j++)
                {
                    CodeBook[cb_size] = coeff[j];
                    cb_size++;
                }
            }
        }
    }

    if (flag_32K_env_ho)
    {
        for (j = 0; j < cb_size; j++)
        {
            if (CodeBook_mod[j] != 0.0f)
            {
                /* Densify codebook */
                CodeBook[j] = sign(CodeBook_mod[j])*(CodeBook_mod[j]*CodeBook_mod[j] + CodeBook_mod[cb_size-j-1]*CodeBook_mod[cb_size-j-1]);
            }
            else
            {
                CodeBook[j] = CodeBook_mod[cb_size-j-1];
            }
        }
    }

    return cb_size;
}

/*--------------------------------------------------------------------------*
 * find_last_band()
 *
 * Find the last band which has bits allocated
 *--------------------------------------------------------------------------*/

short find_last_band(           /* o  : index of last band              */
    const short *bitalloc,      /* i  : bit allocation                  */
    const short nb_sfm          /* i  : number of possibly coded bands  */
)
{
    short sfm, core_sfm;

    core_sfm = nb_sfm-1;

    for (sfm = nb_sfm-1; sfm >= 0; sfm--)
    {
        if ( bitalloc[sfm] != 0 )
        {
            core_sfm = sfm;
            break;
        }
    }

    return core_sfm;
}

/*--------------------------------------------------------------------------*
 * apply_noisefill_HQ()
 *
 * Inject noise in non-coded bands
 *--------------------------------------------------------------------------*/

void apply_noisefill_HQ(
    const short *R,             /* i  : bit allocation                  */
    const short length,         /* i  : input frame length              */
    const short flag_32K_env_ho,/* i  : envelope stability hangover flag*/
    const long  core_brate,     /* i  : core bit rate                   */
    const short last_sfm,       /* i  : last coded subband              */
    const float *CodeBook,      /* i  : Noise-fill codebook             */
    const float *CodeBook_mod,  /* i  : Densified noise-fill codebook   */
    const short cb_size,        /* i  : Codebook length                 */
    const short *sfm_start,     /* i  : Subband start coefficient       */
    const short *sfm_end,       /* i  : Subband end coefficient         */
    const short *sfmsize,       /* i  : Subband band width              */
    float *coeff          /* i/o: coded/noisefilled spectrum      */
)
{
    short sfm;
    short cb_pos;
    float E_cb_vec;
    float E_corr;
    float cb_buff[64];
    short i,j;

    if( length >= L_FRAME32k || core_brate > HQ_32k || core_brate < HQ_24k40 )
    {
        /* Read from codebook */
        cb_pos = 0;

        for (sfm = 0; sfm <= last_sfm; sfm++)
        {
            if (R[sfm] == 0)
            {
                if (flag_32K_env_ho)
                {
                    E_cb_vec = 0.0f;
                    if (sfm < 20)
                    {
                        for (i = 0; i < sfmsize[sfm]; i++)
                        {
                            cb_buff[i] = CodeBook_mod[cb_pos];
                            E_cb_vec += cb_buff[i]*cb_buff[i];
                            cb_pos++;
                            if(cb_pos >= cb_size)
                            {
                                cb_pos = 0;
                            }
                        }
                    }
                    else
                    {
                        for (i = 0; i < sfmsize[sfm]; i++)
                        {
                            cb_buff[i] = CodeBook[cb_pos];
                            E_cb_vec += cb_buff[i]*cb_buff[i];
                            cb_pos++;
                            if(cb_pos >= cb_size)
                            {
                                cb_pos = 0;
                            }
                        }
                    }

                    E_corr = E_cb_vec / ((float) sfmsize[sfm]);
                    E_corr = 1.0f / (float)sqrt(E_corr);

                    for (j = sfm_start[sfm]; j < sfm_end[sfm]; j++)
                    {
                        coeff[j] = cb_buff[j - sfm_start[sfm]] * E_corr;
                    }
                }
                else
                {
                    for (j = sfm_start[sfm]; j < sfm_end[sfm]; j++)
                    {
                        coeff[j] = CodeBook[cb_pos];
                        cb_pos++;
                        cb_pos = (cb_pos>=cb_size) ? 0 : cb_pos;
                    }
                }
            }
        }
    }

    return;
}


/*--------------------------------------------------------------------------*
 * harm_bwe_fine()
 *
 * Prepare harmonic BWE fine structure
 *--------------------------------------------------------------------------*/

void harm_bwe_fine(
    const short *R,                 /* i  : bit allocation                              */
    const short last_sfm,           /* i  : last coded subband                          */
    const short high_sfm,           /* i  : higher transition band to BWE               */
    const short num_sfm,            /* i  : total number of bands                       */
    const short *norm,              /* i  : quantization indices for norms              */
    const short *sfm_start,         /* i  : Subband start coefficient                   */
    const short *sfm_end,           /* i  : Subband end coefficient                     */
    short *prev_L_swb_norm,   /* i/o: last normalize length                       */
    float *coeff,             /* i/o: coded/noisefilled normalized spectrum       */
    float *coeff_out,         /* o  : coded/noisefilled spectrum                  */
    float *coeff_fine         /* o  : BWE fine structure                          */
)
{
    short sfm;
    short i;
    float normq;
    float SWB_signal[L_HARMONIC_EXC];
    float envelope[L_HARMONIC_EXC];
    float *src, *dst, *end;

    short norm_width = 64;

    /* shape the spectrum */
    for (sfm = 0; sfm <= last_sfm; sfm++)
    {
        if( R[sfm] != 0 )
        {
            normq = dicn[norm[sfm]];
            for (i = sfm_start[sfm]; i < sfm_end[sfm]; i++)
            {
                coeff_out[i] = coeff[i]*normq;
            }
        }
        else
        {
            for (i = sfm_start[sfm]; i < sfm_end[sfm]; i++)
            {
                coeff_out[i] = 0.0f;
            }
        }
    }

    /* excitation replication */
    mvr2r( coeff_out, SWB_signal, L_HARMONIC_EXC );
    calc_normal_length( HQ_CORE, coeff_out, HQ_HARMONIC, -1, &norm_width, prev_L_swb_norm );
    hq_swb_harmonic_calc_norm_envelop( SWB_signal, envelope, norm_width, L_HARMONIC_EXC );

    /* Normalize with envelope */
    for (i = 0; i < L_HARMONIC_EXC; i++)
    {
        SWB_signal[i] = SWB_signal[i] / envelope[i];
    }

    dst = coeff_fine + sfm_end[last_sfm];
    end = coeff_fine + sfm_end[num_sfm-1];

    if ( (sfm_end[last_sfm] - sfm_end[high_sfm]) <= L_HARMONIC_EXC - START_EXC )
    {
        src = SWB_signal + START_EXC + (sfm_end[last_sfm] - sfm_end[high_sfm]);
    }
    else
    {
        src = SWB_signal + L_HARMONIC_EXC - 1;
    }

    while (dst < end)
    {
        while (dst < end && src < &SWB_signal[L_HARMONIC_EXC])
        {
            *dst++ = *src++;
        }
        src --;

        while (dst < end && src >= &SWB_signal[START_EXC])
        {
            *dst++ = *src--;
        }
        src++;
    }

    return;
}

/*--------------------------------------------------------------------------*
 * hvq_bwe_fine()
 *
 * Prepare HVQ BWE fine structure
 *--------------------------------------------------------------------------*/

void hvq_bwe_fine(
    const short last_sfm,           /* i  : last coded subband                          */
    const short num_sfm,            /* i  : total number of bands                       */
    const short *sfm_end,           /* i  : Subband end coefficient                     */
    const short *peak_idx,          /* i  : Peak index                                  */
    const short Npeaks,             /* i  : Number of peaks                             */
    short *peak_pos,          /* i/o: Peak positions                              */
    short *prev_L_swb_norm,   /* i/o: last normalize length                       */
    float *coeff,             /* i/o: coded/noisefilled normalized spectrum       */
    short *bwe_peaks,         /* o  : Positions of peaks in BWE                   */
    float *coeff_fine         /* o  : HVQ BWE fine structure                      */
)
{
    short i, j;
    float SWB_signal[L_HARMONIC_EXC];
    float envelope[L_HARMONIC_EXC];
    float *src, *dst, *end;
    short *peak_dst, *peak_src;
    short norm_width = 64;

    /* excitation replication */
    mvr2r( coeff, SWB_signal, L_HARMONIC_EXC );
    calc_normal_length( HQ_CORE, coeff, HQ_HVQ, -1, &norm_width, prev_L_swb_norm );

    hq_swb_harmonic_calc_norm_envelop( SWB_signal, envelope, norm_width, L_HARMONIC_EXC );

    /* Normalize with envelope */
    for (i = 0; i < L_HARMONIC_EXC; i++)
    {
        SWB_signal[i] = SWB_signal[i] / envelope[i];
    }

    dst = coeff_fine;
    end = coeff_fine + sfm_end[num_sfm-1] - sfm_end[last_sfm];

    src = SWB_signal + START_EXC;
    peak_src = peak_pos + START_EXC;

    for(i = 0; i < Npeaks; i++)
    {
        if ( peak_idx[i] < L_HARMONIC_EXC )
        {
            peak_pos[peak_idx[i]] = 1;
        }
    }

    i = L_HARMONIC_EXC-1;
    while( i-- > 0 )
    {
        if( peak_pos[i] == 1 )
        {
            break;
        }
    }

    if( i < 180 )
    {
        i = 180;
    }

    for( j = L_HARMONIC_EXC-1; j > i+1; j-- )
    {
        SWB_signal[j] = 0.0f;
    }

    peak_dst = bwe_peaks + sfm_end[last_sfm];
    while ( dst < end )
    {
        while ( dst < end && src < &SWB_signal[L_HARMONIC_EXC] )
        {
            *dst++ = *src++;
            *peak_dst++ = *peak_src++;
        }
        peak_src--;
        src --;

        while (dst < end && src >= &SWB_signal[START_EXC])
        {
            *dst++ = *src--;
            *peak_dst++ = *peak_src--;
        }
        peak_src++;
        src++;
    }

    return;
}

/*--------------------------------------------------------------------------*
 * hq_fold_bwe()
 *
 * HQ mode folding BWE
 *--------------------------------------------------------------------------*/

void hq_fold_bwe(
    const short last_sfm,           /* i  : last coded subband                      */
    const short *sfm_end,           /* i  : Subband end coefficient                 */
    const short num_sfm,            /* i  : Number of subbands                      */
    float *coeff              /* i/o: coded/noisefilled normalized spectrum   */
)
{
    short low_coeff;
    short first_coeff;
    float *src, *dst, *end;

    /* Find replication range for BWE */
    low_coeff = sfm_end[last_sfm] >> 1;
    src = coeff + sfm_end[last_sfm] - 1;

    first_coeff = sfm_end[last_sfm];
    dst = coeff + sfm_end[last_sfm];
    end = coeff + sfm_end[num_sfm-1];

    /* Generate BWE with spectral folding */
    while (dst < end)
    {
        while (dst < end && src >= &coeff[low_coeff])
        {
            *dst++ = *src--;
        }

        src++;

        while (dst < end && src < &coeff[first_coeff])
        {
            *dst++ = *src++;
        }
    }

    return;
}

/*--------------------------------------------------------------------------*
 * apply_nf_gain()
 *
 * Apply noise fill gain
 *--------------------------------------------------------------------------*/

void apply_nf_gain(
    const short nf_idx,             /* i  : noise fill gain index                   */
    const short last_sfm,           /* i  : last coded subband                      */
    const short *R,                 /* i  : bit allocation                          */
    const short *sfm_start,         /* i  : Subband start coefficient               */
    const short *sfm_end,           /* i  : Subband end coefficient                 */
    float *coeff              /* i/o: coded/noisefilled normalized spectrum   */
)
{
    short sfm;
    short j;
    float nf_scale;

    nf_scale = 1.0f / (1 << nf_idx);
    for (sfm = 0; sfm <= last_sfm; sfm++)
    {
        if (R[sfm] == 0)
        {
            for (j = sfm_start[sfm]; j < sfm_end[sfm]; j++)
            {
                /* Scale NoiseFill */
                coeff[j] = coeff[j] * nf_scale;
            }
        }
    }

    return;
}

/*--------------------------------------------------------------------------*
 * hq_generic_fine()
 *
 * Prepare HQ SWB BWE fine structure
 *--------------------------------------------------------------------------*/

void hq_generic_fine(
    float *coeff,             /* i  : coded/noisefilled normalized spectrum   */
    const short last_sfm,           /* i  : Last coded band                         */
    const short *sfm_start,         /* i  : Subband start coefficient               */
    const short *sfm_end,           /* i  : Subband end coefficient                 */
    short *bwe_seed,          /* i/o: random seed for generating BWE input    */
    float *coeff_out1         /* o  : HQ SWB BWE input                        */
)
{
    short sfm;
    short i;
    float multi;

    multi = 1.0f;
    for (sfm = 0; sfm <= last_sfm; sfm++)
    {
        for (i = sfm_start[sfm]; i < sfm_end[sfm]; i++)
        {
            if (coeff[i]==0.f)
            {
                coeff_out1[i] = (float)multi*( own_random(bwe_seed)>0 ? 1.0f: -1.0f)*0.5f;
            }
            else
            {
                coeff_out1[i] = coeff[i];
            }
        }
    }

    return;
}

/*--------------------------------------------------------------------------*
 * harm_bwe()
 *
 * HQ Harmonic BWE
 *--------------------------------------------------------------------------*/

void harm_bwe(
    const float *coeff_fine,        /* i  : fine structure for BWE                  */
    const float *coeff,             /* i  : coded/noisefilled normalized spectrum   */
    const short num_sfm,            /* i  : Number of subbands                      */
    const short *sfm_start,         /* i  : Subband start coefficient               */
    const short *sfm_end,           /* i  : Subband end coefficient                 */
    const short last_sfm,           /* i  : last coded subband                      */
    const short high_sfm,           /* i  : higher transition band to BWE           */
    const short *R,                 /* i  : bit allocation                          */
    const short prev_hq_mode,       /* i  : previous hq mode                        */
    short *norm,              /* i/o: quantization indices for norms          */
    float *noise_level,       /* i/o: noise levels for harmonic modes         */
    float *prev_noise_level,  /* i/o: noise factor in previous frame          */
    short *bwe_seed,          /* i/o: random seed for generating BWE input    */
    float *coeff_out          /* o  : coded/noisefilled spectrum              */
)
{
    short i, j;
    short sfm;
    float normq;
    float norm_adj;
    float E_L;

    float alfa = 0.5f;
    float alpha, beta;
    short idx;
    float fac;
    float *src, *dst;

    for (sfm = 0; sfm <= last_sfm; sfm++)
    {
        if (R[sfm] == 0)
        {
            normq = dicn[norm[sfm]];
            for (i = sfm_start[sfm]; i < sfm_end[sfm]; i++)
            {
                coeff_out[i] = coeff[i]*normq;
            }
        }
    }
    noise_level[1] = noise_level[0];

    /* shaping the BWE spectrum further by envelopes and noise factors */
    noise_level[0] = 0.9f*prev_noise_level[0] + 0.1f*noise_level[0];
    noise_level[1] = 0.9f*prev_noise_level[1] + 0.1f*noise_level[1];

    if( prev_hq_mode == HQ_NORMAL || prev_hq_mode == HQ_GEN_SWB)
    {
        if( noise_level[0] < 0.25f )
        {
            noise_level[0] *= 4.0f;
        }

        if( noise_level[1] < 0.25f )
        {
            noise_level[1] *= 4.0f;
        }
    }

    E_L = EPSILON;
    for( i=last_sfm+1; i < num_sfm; i++ )
    {
        E_L = EPSILON;
        for (j = sfm_start[i]; j < sfm_end[i]; j++)
        {
            E_L += coeff_fine[j] * coeff_fine[j];
        }
        E_L = (float)sqrt((sfm_end[i] - sfm_start[i])/E_L);

        normq = dicn[norm[i]];
        norm_adj = normq*E_L;
        alfa = (i > 27) ? noise_level[1] : noise_level[0];
        beta = (float)sqrt(1.0f - alfa);
        alpha = (float)sqrt(alfa) * 0.5f;

        for(sfm = sfm_start[i]; sfm < sfm_end[i]; sfm++)
        {
            coeff_out[sfm] = (beta*coeff_fine[sfm]*norm_adj + alpha*own_random(bwe_seed)/32768.0f*normq);
        }
    }

    prev_noise_level[0] = noise_level[0];
    prev_noise_level[1] = noise_level[1];
    idx = 16;
    src = &coeff_out[sfm_end[high_sfm] + L_HARMONIC_EXC - START_EXC];
    dst = src-1;

    for( i = 0; i < idx; i++ )
    {
        fac = i / (2.0f * idx);
        *src++ *= 0.5f + fac;
        *dst-- *= 0.5f + fac;
    }
    if( num_sfm ==33 )
    {
        set_f(&coeff_out[800], 0, 160);
    }
    return;
}

/*--------------------------------------------------------------------------*
 * HVQ_bwe()
 *
 * HQ HVQ BWE
 *--------------------------------------------------------------------------*/

void hvq_bwe(
    const float *coeff,             /* i  : coded/noisefilled spectrum              */
    const float *coeff_fine,        /* i  : BWE fine structure                      */
    const short *sfm_start,         /* i  : Subband start coefficient               */
    const short *sfm_end,           /* i  : Subband end coefficient                 */
    const short *sfm_len,           /* i  : Subband length                          */
    const short last_sfm,           /* i  : last coded subband                      */
    const short prev_hq_mode,       /* i  : previous hq mode                        */
    const short *bwe_peaks,         /* i  : HVQ bwe peaks                           */
    const short bin_th,             /* i  : HVQ transition bin                      */
    const short num_sfm,            /* i  : Number of bands                         */
    const long  core_brate,         /* i  : Core bit-rate                           */
    const short *R,                 /* i  : Bit allocation                          */
    short *norm,              /* i/o: quantization indices for norms          */
    float *noise_level,       /* i/o: noise levels for harmonic modes         */
    float *prev_noise_level,  /* i/o: noise factor in previous frame          */
    short *bwe_seed,          /* i/o: random seed for generating BWE input    */
    float *coeff_out          /* o  : coded/noisefilled spectrum              */
)
{
    short i, j;
    short N;
    float normq;
    float E_L;

    short bwe_noise_th = 0;
    short peak_band, low, high, sel_norm;
    short norm_ind;
    float tmp_norm = 0;
    short idx;
    float fac;
    float *src, *dst;
    short istart, iend;
    short offset = sfm_end[last_sfm];

    mvr2r( coeff, coeff_out, L_FRAME48k );

    bwe_noise_th = bin_th+(sfm_end[num_sfm-1] - bin_th)/HVQ_BWE_NOISE_BANDS;
    logqnorm(&coeff_out[sfm_start[last_sfm]], &norm[last_sfm], 40, sfm_len[last_sfm], thren);

    /* shaping the BWE spectrum further by envelopes and noise factors */
    noise_level[0] = 0.9f*prev_noise_level[0] + 0.1f*noise_level[0];
    noise_level[1] = 0.9f*prev_noise_level[1] + 0.1f*noise_level[1];

    if( prev_hq_mode == HQ_NORMAL || prev_hq_mode == HQ_GEN_SWB)
    {
        if( noise_level[0] < 0.25f )
        {
            noise_level[0] *= 4.0f;
        }

        if( noise_level[1] < 0.25f )
        {
            noise_level[1] *= 4.0f;
        }
    }

    norm_ind = last_sfm+1;
    if (core_brate == HQ_24k40)
    {
        peak_band = 0;
        E_L = EPSILON;
        for (i = sfm_start[norm_ind]; i < sfm_end[norm_ind+1]; i++)
        {
            if (bwe_peaks[i])
            {
                peak_band = 1;
            }
            E_L += coeff_fine[i-offset] * coeff_fine[i-offset];
        }
        E_L = (float)sqrt((sfm_end[norm_ind+1] - sfm_start[norm_ind])/E_L);

        normq = 0.1f*dicn[norm[norm_ind-1]] + 0.8f*dicn[norm[norm_ind]] + 0.1f*dicn[norm[norm_ind+1]];
        tmp_norm = 0.1f*dicn[norm[norm_ind]] + 0.8f*dicn[norm[norm_ind+1]] + 0.1f*dicn[norm[norm_ind+2]];

        istart = sfm_start[norm_ind];
        iend = istart + sfm_len[norm_ind]/2;
        for (i = istart; i < iend; i++)
        {
            coeff_out[i] = ((1.0f - noise_level[i/bwe_noise_th])*coeff_fine[i-offset]*E_L + noise_level[i/bwe_noise_th]*own_random(bwe_seed)/32768.0f)*normq;
        }

        j = 0;
        N = sfm_len[norm_ind]/2+sfm_len[norm_ind+1]/2-1;
        istart = iend;
        iend = sfm_start[norm_ind+1] + sfm_len[norm_ind+1]/2;
        for (i = istart; i < iend; i++)
        {
            coeff_out[i] = ((float)(N-j)/N*normq + (float)j/N*tmp_norm)*((1.0f - noise_level[i/bwe_noise_th])*coeff_fine[i-offset]*E_L + noise_level[i/bwe_noise_th]*own_random(bwe_seed)/32768.0f);
            j++;
        }

        istart = iend;
        iend = sfm_end[norm_ind+1];
        for (i = istart; i < iend; i++)
        {
            coeff_out[i] = ((1.0f - noise_level[i/bwe_noise_th])*coeff_fine[i-offset]*E_L + noise_level[i/bwe_noise_th]*own_random(bwe_seed)/32768.0f)*tmp_norm;
        }

        norm_ind += 2;
    }

    for ( ; norm_ind < num_sfm; norm_ind++)
    {
        if ( R[norm_ind] == 0 )
        {
            peak_band = 0;
            E_L = EPSILON;

            for (i = sfm_start[norm_ind]; i < sfm_end[norm_ind]; i++)
            {
                if (bwe_peaks[i])
                {
                    peak_band = 1;
                    break;
                }
            }

            istart = sfm_start[norm_ind];
            iend = sfm_end[norm_ind];

            if ( peak_band == 1 && norm_ind > last_sfm+1 && norm_ind < num_sfm-1 )
            {
                istart -= sfm_len[norm_ind-1]/2;
                iend += sfm_len[norm_ind+1]/2;
            }

            for (i = istart; i < iend; i++)
            {
                E_L += coeff_fine[i-offset] * coeff_fine[i-offset];
            }
            E_L = (float)sqrt((iend - istart)/E_L);

            if ( peak_band )
            {
                if ( norm_ind+1 > num_sfm-1 )
                {
                    normq = 0.15f*dicn[norm[norm_ind-1]] + 0.85f*dicn[norm[norm_ind]];
                }
                else
                {
                    normq = 0.1f*dicn[norm[norm_ind-1]] + 0.8f*dicn[norm[norm_ind]] + 0.1f*dicn[norm[norm_ind+1]];
                }
            }
            else
            {
                low = norm_ind;
                high = min(norm_ind+1, num_sfm-1);
                sel_norm = norm[norm_ind-1];
                for (j = low; j <= high; j++)
                {
                    if (norm[j] > sel_norm)
                    {
                        sel_norm = norm[j];
                    }
                }
                normq = dicn[sel_norm];
            }

            for (i = sfm_start[norm_ind]; i < sfm_end[norm_ind]; i++)
            {
                coeff_out[i] = ((1.0f - noise_level[i/bwe_noise_th])*coeff_fine[i-offset]*E_L + noise_level[i/bwe_noise_th]*own_random(bwe_seed)/32768.0f)*normq;
            }
        }
        else /* R[norm_ind] > 0 */
        {
            for (i = sfm_start[norm_ind]; i < sfm_end[norm_ind]; i++)
            {
                coeff_out[i] = coeff[i]; /* Scaling already applied */
            }
        }

    }

    prev_noise_level[0] = noise_level[0];
    prev_noise_level[1] = noise_level[1];
    idx = 16;
    src = &coeff_out[sfm_end[last_sfm] + L_HARMONIC_EXC - START_EXC];
    dst = src-1;

    for( i = 0; i < idx; i++ )
    {
        fac = i/(2.0f * idx);
        *src++ *= 0.5f + fac;
        *dst-- *= 0.5f + fac;
    }

    return;
}
/*-------------------------------------------------------------------*
* hvq_concat_bands()
*
* Compute the band limits for concatenated bands for PVQ target signal in HVQ
*--------------------------------------------------------------------------*/
void hvq_concat_bands
(
    const short pvq_bands,          /* i  : Number of bands in concatenated PVQ target  */
    const short *sel_bnds,          /* i  : Array of selected high bands                */
    const short n_sel_bnds,         /* i  : Number of selected high bands               */
    short *hvq_band_start,    /* i  : Band start indices                          */
    short *hvq_band_width,    /* i  : Band widths                                 */
    short *hvq_band_end       /* i  : Band end indices                            */
)
{
    short k;
    short s;

    s = 0;

    for (k = 0; k < pvq_bands; k++)
    {

        if( k >= pvq_bands - n_sel_bnds)
        {
            hvq_band_start[k] = hvq_band_end[k-1];
            hvq_band_width[k] = band_len_harm[sel_bnds[s]];
            hvq_band_end[k]   = hvq_band_end[k-1] + band_len_harm[sel_bnds[s]];
            s++;
        }
        else
        {
            hvq_band_start[k] = k * HVQ_PVQ_COEFS;
            hvq_band_width[k] = HVQ_PVQ_COEFS;
            hvq_band_end[k]   = (k + 1) * HVQ_PVQ_COEFS;
        }
    }

    return;
}


/*--------------------------------------------------------------------------*
 * map_hq_generic_fenv_norm()
 *
 * mapping high frequency envelope to high band norm
 *--------------------------------------------------------------------------*/
void map_hq_generic_fenv_norm(
    const short hqswb_clas,
    const float *hq_generic_fenv,
    short *ynrm,
    short *normqlg2,
    const short num_env_bands,
    const short nb_sfm,
    const short hq_generic_offset)
{
    float env_fl[17];
    short i;

    set_f( env_fl, 0, 17 );
    if (hq_generic_offset == 144)
    {
        env_fl[0] = hq_generic_fenv[1];
        env_fl[1] = hq_generic_fenv[2]*0.6640625f+hq_generic_fenv[3]*0.3359375f;
        env_fl[2] = hq_generic_fenv[3]*0.6640625f+hq_generic_fenv[4]*0.3359375f;
        env_fl[3] = hq_generic_fenv[4]*0.3359375f+hq_generic_fenv[5]*0.6640625f;
        env_fl[4] = hq_generic_fenv[5]*0.3359375f+hq_generic_fenv[6]*0.6640625f;
        env_fl[5] = hq_generic_fenv[7];
        env_fl[6] = hq_generic_fenv[8]*0.75f+hq_generic_fenv[9]*0.25f;
        env_fl[7] = hq_generic_fenv[9]*0.75f+hq_generic_fenv[10]*0.25f;
        env_fl[8] = hq_generic_fenv[10]*0.25f+hq_generic_fenv[11]*0.75f;
    }
    else
    {
        env_fl[0] = hq_generic_fenv[0]*0.3359375f+hq_generic_fenv[1]*0.6640625f;
        env_fl[1] = hq_generic_fenv[1]*0.3359375f+hq_generic_fenv[2]*0.6640625f;
        env_fl[2] = hq_generic_fenv[3];
        env_fl[3] = hq_generic_fenv[4]*0.6640625f+hq_generic_fenv[5]*0.3359375f;
        env_fl[4] = hq_generic_fenv[5]*0.6640625f+hq_generic_fenv[6]*0.3359375f;
        env_fl[5] = hq_generic_fenv[6]*0.3359375f+hq_generic_fenv[7]*0.6640625f;
        env_fl[6] = hq_generic_fenv[7]*0.3359375f+hq_generic_fenv[8]*0.6640625f;
        env_fl[7] = hq_generic_fenv[8]*0.3359375f+hq_generic_fenv[9]*0.6640625f;
        env_fl[8] = hq_generic_fenv[9]*0.3359375f+hq_generic_fenv[10]*0.6640625f;
        env_fl[9] = hq_generic_fenv[10]*0.25f+hq_generic_fenv[11]*0.75f;
        env_fl[10] = hq_generic_fenv[12];
        env_fl[11] = hq_generic_fenv[13];
    }

    if (hqswb_clas == HQ_GEN_FB)
    {
        if (hq_generic_offset == 144)
        {
            env_fl[9] = hq_generic_fenv[12];
            env_fl[10] = hq_generic_fenv[12]*0.25f+hq_generic_fenv[13]*0.75f;
            env_fl[11] = hq_generic_fenv[13]*0.5f+hq_generic_fenv[14]*0.5f;
            env_fl[12] = hq_generic_fenv[14];
            env_fl[13] = hq_generic_fenv[14];
        }
        else
        {
            env_fl[12] = hq_generic_fenv[14];
            env_fl[13] = hq_generic_fenv[14]*0.25f+hq_generic_fenv[15]*0.75f;
            env_fl[14] = hq_generic_fenv[15]*0.5f+hq_generic_fenv[16]*0.5f;
            env_fl[15] = hq_generic_fenv[16];
            env_fl[16] = hq_generic_fenv[16];
        }
    }

    logqnorm_2( env_fl,40, num_env_bands, nb_sfm, ynrm+num_env_bands, normqlg2+num_env_bands, thren);

    for(i=num_env_bands; i<nb_sfm; ++i)
    {
        normqlg2[i] = dicnlg2[min(ynrm[i]+10,39)];
    }

    return;
}

static void update_rsubband(const short nb_sfm,
                            short *Rsubband,
                            short b_add_bits_denv
                           )
{
    short i;

    /* updating bit allocation */
    while (b_add_bits_denv>0)
    {
        i=nb_sfm-1;
        while(b_add_bits_denv>0 && i>=0)
        {
            if (Rsubband[i]>24)
            {
                Rsubband[i] -= 8;
                b_add_bits_denv--;
            }
            i--;
        }
    }

    return;
}

short get_nor_delta_hf(
    Decoder_State *st,
    short *ynrm,
    short *Rsubband,
    const short num_env_bands,
    const short nb_sfm,
    const short core_sfm
)
{
    short i;
    short delta,bitsforDelta,add_bits_denv;

    add_bits_denv = 0;
    if (core_sfm >= num_env_bands)
    {
        bitsforDelta = (short)get_next_indice(st,2);
        bitsforDelta += 2;
        add_bits_denv += 2;

        for(i=num_env_bands; i<nb_sfm; ++i)
        {
            if (Rsubband[i]!=0)
            {
                delta = (short)get_next_indice(st,bitsforDelta);
                ynrm[i] += delta - (1<<(bitsforDelta-1));

                /* safety check in case of bit errors */
                if ( ynrm[i] < 0 || ynrm[i] > 39 )
                {
                    ynrm[i] = 39;
                    st->BER_detect = 1;
                }
                add_bits_denv += bitsforDelta;
            }
        }

        update_rsubband(nb_sfm, Rsubband,add_bits_denv);
    }
    return add_bits_denv;
}

short calc_nor_delta_hf(
    Encoder_State *st,
    const float *t_audio,
    short *ynrm,
    short *Rsubband,
    const short num_env_bands,
    const short nb_sfm,
    const short *sfmsize,
    const short *sfm_start,
    const short core_sfm
)
{
    short i;
    short ynrm_t[44],normqlg2_t[44];
    short delta,max_delta,min_delta,bitsforDelta,add_bits_denv;


    short temp_num=0;

    max_delta=-100;
    calc_norm( t_audio, ynrm_t, normqlg2_t, 0, nb_sfm, sfmsize, sfm_start );
    add_bits_denv = 0;
    for(i=num_env_bands; i<nb_sfm; ++i)
    {
        if (Rsubband[i]!=0)
        {
            delta = ynrm_t[i] - ynrm[i];
            if (delta > 0)
            {
                delta += 1;
            }
            else
            {
                delta = -delta;
            }
            if (delta>max_delta)
            {
                max_delta = delta;
            }
        }
    }
    if (core_sfm >= num_env_bands)
    {
        if (max_delta < 16)
        {
            bitsforDelta = 2;
            while(max_delta>=2)
            {
                bitsforDelta++;
                max_delta >>= 1;
            }
        }
        else
        {
            bitsforDelta = 5;
        }
        max_delta = (1<<(bitsforDelta-1))-1;
        min_delta = (max_delta+1)*(-1);

        /* updating norm & storing delta norm */
        add_bits_denv = 2;
        push_indice( st, IND_DELTA_ENV_HQ, bitsforDelta-2 , 2 );
        for(i=num_env_bands; i<nb_sfm; ++i)
        {
            if (Rsubband[i]!=0)
            {
                delta = ynrm_t[i] - ynrm[i];
                if (delta > max_delta)
                {
                    delta = max_delta;
                }
                else if (delta < min_delta)
                {
                    delta = min_delta;
                }
                push_indice( st, IND_DELTA_ENV_HQ, delta - min_delta , bitsforDelta );
                ynrm[i] += delta;
                add_bits_denv += bitsforDelta;


                temp_num++;
            }
        }

        /* updating bit allocation */
        update_rsubband(nb_sfm, Rsubband,add_bits_denv);

    }

    return add_bits_denv;
}

/*-------------------------------------------------------------------*
* hq_generic_bwe()
*
* HQ GENERIC BWE
*--------------------------------------------------------------------------*/
void hq_generic_bwe(
    const short HQ_mode,                  /* i  : HQ mode                                     */
    float *coeff_out1,              /* i/o: BWE input & temporary buffer                */
    const float *hq_generic_fenv,         /* i  : SWB frequency envelopes                     */
    float *coeff_out,               /* o  : SWB signal in MDCT domain                   */
    const short hq_generic_offset,        /* i  : frequency offset for representing hq generic*/
    short *prev_L_swb_norm,         /* i/o: last normalize length                       */
    const short hq_generic_exc_clas,      /* i  : hq generic hf excitation class              */
    const short *sfm_end,                 /* i  : End of bands                                */
    const short num_sfm,                  /* i  : Number of bands                             */
    const short num_env_bands,            /* i  : Number of coded envelope bands              */
    const short *R                        /* i  : Bit allocation                              */
)
{
    short n_swb_overlap_offset, n_swb_overlap;
    float hq_swb_overlap_buf[640];

    n_swb_overlap_offset = swb_bwe_subband[0] + hq_generic_offset;

    n_swb_overlap = sfm_end[num_env_bands-1] - n_swb_overlap_offset;
    mvr2r( &coeff_out[n_swb_overlap_offset], hq_swb_overlap_buf, n_swb_overlap + sfm_end[num_sfm-1] - sfm_end[num_env_bands-1] );

    hq_generic_hf_decoding( HQ_mode, coeff_out1, hq_generic_fenv, coeff_out, hq_generic_offset, prev_L_swb_norm, hq_generic_exc_clas, R );

    overlap_hq_bwe( hq_swb_overlap_buf, coeff_out, n_swb_overlap_offset, n_swb_overlap, R, num_env_bands, num_sfm, sfm_end );

    return;
}


/*--------------------------------------------------------------------------*
 * hq_wb_nf_bwe()
 *
 * HQ WB noisefill and BWE
 *--------------------------------------------------------------------------*/

void hq_wb_nf_bwe(
    const float *coeff,             /* i  : coded/noisefilled normalized spectrum   */
    const short is_transient,       /* i  : is transient flag                       */
    const short prev_bfi,           /* i  : previous bad frame indicator            */
    const float *normq_v,           /* i  : norms                                   */
    const short num_sfm,            /* i  : Number of subbands                      */
    const short *sfm_start,         /* i  : Subband start coefficient               */
    const short *sfm_end,           /* i  : Subband end coefficient                 */
    const short *sfmsize,           /* i  : Subband band width                      */
    const short last_sfm,           /* i  : last coded subband                      */
    const short *R,                 /* i  : bit allocation                          */
    const short prev_is_transient,  /* i  : previous transient flag                 */
    float *prev_normq,        /* i/o: previous norms                          */
    float *prev_env,          /* i/o: previous noise envelopes                */
    short *bwe_seed,          /* i/o: random seed for generating BWE input    */
    float *prev_coeff_out,    /* i/o: decoded spectrum in previous frame      */
    short *prev_R,            /* i/o: bit allocation info. in previous frame  */
    float *coeff_out          /* o  : coded/noisefilled spectrum              */
)
{
    short i;
    short sfm;
    short total_bit;
    short num;

    float bitalloc_var;
    float sharp;
    float mean;
    float peak;
    float fabs_coeff_out;
    float harm_para;
    float alfa = 0.5;
    float env;
    float step;
    float min_coef;
    float avrg_norm;
    float prev_avrg_norm;

    if( is_transient == 0 )
    {
        if( prev_bfi == 1 )
        {
            mvr2r(normq_v, prev_normq, SFM_N_WB);
        }

        /* the variance of bit allocation */
        total_bit = 0;
        bitalloc_var = 0.0f;
        for (sfm = 8; sfm <= last_sfm; sfm++)
        {
            bitalloc_var += (float)abs(R[sfm] - R[sfm-1]);
            total_bit += R[sfm];
        }
        bitalloc_var = (last_sfm > 8 && total_bit > 0) ? (bitalloc_var / total_bit) : 0;

        /* calculate the peak-average ratio of saturable subbands */
        num = 0;
        sharp = EPSILON;
        for(sfm = last_sfm; sfm >= 8; sfm--)
        {
            if(R[sfm] >= rat[sfm]*sfmsize[sfm])
            {
                peak = 0.0f;
                mean = EPSILON;
                for (i = sfm_start[sfm]; i < sfm_end[sfm]; i++)
                {
                    fabs_coeff_out = (float)fabs(coeff_out[i]);
                    mean += fabs_coeff_out;
                    if(fabs_coeff_out > peak)
                    {
                        peak = fabs_coeff_out;
                    }
                }
                sharp += sfmsize[sfm]*peak/mean;
                num ++;
            }
        }

        sharp = (num != 0) ? 2.0f*num/sharp : 1.0f;
        harm_para = sharp;
        if( last_sfm == 0 )
        {
            step = 0;
        }
        else
        {
            step = 5.0f * sharp / last_sfm;
        }
        alfa = 2.5f;

        /* fill noise for the insaturable subbands */
        for( sfm = 0; sfm < num_sfm; sfm++ )
        {
            env = 0.0f;
            if(R[sfm] != 0 && R[sfm] < 1.5f*sfmsize[sfm])
            {
                /* calculate the energy of the undecoded coefficients */
                peak = 0.0f;
                min_coef = FLT_MAX;
                env = normq_v[sfm]*normq_v[sfm]*sfmsize[sfm];
                for (i = sfm_start[sfm]; i < sfm_end[sfm]; i++)
                {
                    fabs_coeff_out = (float)fabs(coeff_out[i]);
                    if(fabs_coeff_out < min_coef && coeff_out[i] != 0)
                    {
                        min_coef = fabs_coeff_out;
                    }
                    if(fabs_coeff_out > peak)
                    {
                        peak = fabs_coeff_out;
                    }
                    env -= coeff_out[i]*coeff_out[i];
                }

                if (env > 0 )
                {
                    if(sfm == 0)
                    {
                        avrg_norm = normq_v[0] + normq_v[1] + normq_v[2];
                        prev_avrg_norm = prev_normq[0] + prev_normq[1] + prev_normq[2];
                    }
                    else if (sfm == 25)
                    {
                        avrg_norm = normq_v[23] + normq_v[24] + normq_v[25];
                        prev_avrg_norm = prev_normq[23] + prev_normq[24] + prev_normq[25];
                    }
                    else
                    {
                        avrg_norm = normq_v[sfm-1] + normq_v[sfm] + normq_v[sfm+1];
                        prev_avrg_norm = prev_normq[sfm-1] + prev_normq[sfm] + prev_normq[sfm+1];
                    }

                    if(bitalloc_var > 0.3f || 4.0f*normq_v[sfm] < peak)
                    {
                        /* calculate the noise magnitude of harmonic signal */
                        env = (float)(avrg_norm*harm_para *sqrt(env/sfmsize[sfm])/peak);
                    }
                    else
                    {
                        /* calculate the noise magnitude of normal signal */
                        env = sharp *(float)sqrt(env/sfmsize[sfm]);
                        if(alfa*normq_v[sfm] < peak)
                        {
                            env *= env/peak;
                        }
                        sharp += step;
                    }
                    if(env > 0.5f*min_coef)
                    {
                        env = 0.5f*min_coef;
                    }

                    if(prev_bfi == 1)
                    {
                        prev_env[sfm] = env;
                    }
                    /* smooth the noise magnitudes between inter-frame */
                    if(prev_avrg_norm > 0.5f*avrg_norm && prev_avrg_norm < 2.0f*avrg_norm && prev_is_transient == 0)
                    {
                        env = 0.5f*env + 0.5f*prev_env[sfm];
                    }

                    for (i = sfm_start[sfm]; i < sfm_end[sfm]; i++)
                    {
                        if (coeff[i] == 0)
                        {
                            coeff_out[i] = (float)( own_random(bwe_seed))/32768.0f;
                            coeff_out[i] *= env;
                        }
                    }
                }
                else
                {
                    env = 0.0f;
                }
            }
            else if(R[sfm] == 0)
            {
                /* fill random noise for 0 bit subbands */
                for (i = sfm_start[sfm]; i < sfm_end[sfm]; i++)
                {
                    if(coeff[i] == 0)
                    {
                        coeff_out[i] = (float)( own_random(bwe_seed))/32768.0f;
                        coeff_out[i] *= normq_v[sfm];
                    }
                }

                env = normq_v[sfm];
            }
            if(sfm == SFM_N_WB-1 && prev_is_transient == 0 && prev_normq[sfm] > 0.5f*normq_v[sfm] && prev_normq[sfm] < 2.0f*normq_v[sfm] && bitalloc_var <= 0.3f )
            {
                float *p_prev_coeff_out = prev_coeff_out;
                for (i = sfm_start[sfm]+12; i < sfm_end[sfm]; i++)
                {
                    if( fabs(coeff_out[i]) > 4.0f * fabs(*p_prev_coeff_out) ||
                            fabs(coeff_out[i]) < 0.25f * fabs(*p_prev_coeff_out) ||
                            (R[sfm] * (*prev_R) == 0 && R[sfm] + (*prev_R) != 0) )
                    {
                        coeff_out[i] = (coeff_out[i] > 0) ? (float)(0.5f*(fabs(coeff_out[i]) + fabs(*p_prev_coeff_out))) : (float)(-0.5f*(fabs(coeff_out[i]) + fabs(*p_prev_coeff_out)));
                    }
                    p_prev_coeff_out++;
                }
            }

            prev_env[sfm] = env;
        }
    }
    else
    {
        /* fill random noise for 0 bit subbands of transient frame */
        for(sfm = 0; sfm < num_sfm; sfm++)
        {
            if( R[sfm] == 0 )
            {
                for (i = sfm_start[sfm]; i < sfm_end[sfm]; i++)
                {
                    coeff_out[i] = (float)( own_random(bwe_seed))/32768.0f;
                    coeff_out[i] *= normq_v[sfm];
                }
            }
        }

        set_f( prev_env, 0, SFM_N_WB );
    }

    mvr2r(normq_v, prev_normq, SFM_N_WB);
    mvr2r( coeff_out + L_FRAME16k - L_HQ_WB_BWE, prev_coeff_out, L_HQ_WB_BWE );
    *prev_R = R[SFM_N_WB-1];

    return;
}


/*--------------------------------------------------------------------------*
 * enforce_zero_for_min_envelope()
 *
 * Detect minimum level of envelope and set corresponding bands to zero
 *--------------------------------------------------------------------------*/

void enforce_zero_for_min_envelope(
    const short hqswb_clas,     /* i  : HQ coding mode                     */
    const short *ynrm,          /* i  : Envelope indices                   */
    float *coefsq,        /* i/o: Quantized spectrum/zeroed spectrum */
    short nb_sfm,         /* i  : Number of coded sub bands          */
    const short *sfm_start,     /* i  : Sub band start indices             */
    const short *sfm_end        /* i  : Sub band end indices               */
)
{
    short i, j;

    /* prevent non-zero output for all-zero input */
    if( hqswb_clas != HQ_HVQ )
    {
        if( ynrm[0] == 31 )
        {
            for( j = sfm_start[0]; j < sfm_end[0]; j++ )
            {
                coefsq[j] = 0.0f;
            }
        }

        for( i = 1; i < nb_sfm; i++ )
        {
            if( ynrm[i] == 39 )
            {
                for( j = sfm_start[i]; j < sfm_end[i]; j++ )
                {
                    coefsq[j] = 0.0f;
                }
            }
        }
    }

    return;
}


/*--------------------------------------------------------------------------*
 * overlap_hq_bwe()
 *
 * Overlapping at the boundary between HQ core and BWE
 *--------------------------------------------------------------------------*/

static void overlap_hq_bwe(
    const float *hq_swb_overlap_buf,    /* i  : spectrum from HQ core                   */
    float *coeff_out,             /* i/o: spectrum from BWE, overlapped output    */
    const short n_swb_overlap_offset,   /* i  : starting offset of overlapping          */
    const short n_swb_overlap,          /* i  : length of overlapping                   */
    const short *R,                     /* i  : Bit allocation                          */
    const short num_env_bands,          /* i  : Number of coded envelope bands          */
    const short num_sfm,                /* i  : Number of bands                         */
    const short *sfm_end                /* i  : Band end indices                        */
)
{
    short i;
    float step;
    float weighting;
    short n_band;

    if (R[num_env_bands-1] != 0)
    {
        mvr2r( hq_swb_overlap_buf, &coeff_out[n_swb_overlap_offset], n_swb_overlap );
    }
    else
    {
        /*weighting = 0.8f;*/
        step = 1.0f / (float)n_swb_overlap;
        weighting = 1.0f;
        for (i = 0; i < n_swb_overlap; i++)
        {
            coeff_out[n_swb_overlap_offset+i] = hq_swb_overlap_buf[i] * weighting + coeff_out[n_swb_overlap_offset+i] * (1.0f - weighting);
            weighting -= step;
        }
    }


    for (n_band = num_env_bands; n_band < num_sfm; n_band++)
    {
        if (R[n_band] !=0 )
        {
            for(i=sfm_end[n_band-1]; i<sfm_end[n_band]; ++i)
            {
                coeff_out[i] = hq_swb_overlap_buf[i-n_swb_overlap_offset];
            }
        }
    }

    return;
}


/*--------------------------------------------------------------------------*
 * apply_envelope()
 *
 * Apply spectral envelope with envelope adjustments
 *--------------------------------------------------------------------------*/

void apply_envelope(
    const float *coeff,             /* i/o: Coded/noisefilled normalized spectrum   */
    const short *norm,              /* i  : Envelope                                */
    const float *norm_adj,          /* i  : Envelope adjustment                     */
    const short num_sfm,            /* i  : Total number of bands                   */
    const short last_sfm,           /* i  : Last coded band                         */
    const short HQ_mode,            /* i  : HQ mode                                 */
    const short length,             /* i  : Frame length                            */
    const short *sfm_start,         /* i  : Sub band start indices                  */
    const short *sfm_end,           /* i  : Sub band end indices                    */
    float *normq_v,           /* o  : Envelope with adjustment                */
    float *coeff_out,         /* o  : coded/noisefilled spectrum              */
    float *coeff_out1         /* o  : noisefilled spectrum for HQ SWE BWE     */
)
{
    short i;
    short sfm;
    float normq;
    short len;

    len = num_sfm;
    if( HQ_mode == HQ_GEN_SWB || HQ_mode == HQ_GEN_FB )
    {
        len = last_sfm + 1;
    }

    if( length == L_FRAME16k )
    {
        for (sfm = 0; sfm < num_sfm; sfm++)
        {
            normq_v[sfm] = dicn[norm[sfm]];
            normq = normq_v[sfm] * norm_adj[sfm];

            for (i = sfm_start[sfm]; i < sfm_end[sfm]; i++)
            {
                coeff_out[i] = coeff[i]*normq;
            }
        }
    }
    else
    {
        for (sfm = 0; sfm < len; sfm++)
        {
            normq_v[sfm] = dicn[norm[sfm]];
            normq_v[sfm] *= norm_adj[sfm];

            normq = normq_v[sfm];
            for (i = sfm_start[sfm]; i < sfm_end[sfm]; i++)
            {
                coeff_out[i] = coeff[i]*normq;
            }
        }

        if ( HQ_mode == HQ_GEN_SWB || HQ_mode == HQ_GEN_FB )
        {
            for (sfm = 0; sfm <= last_sfm; sfm++)
            {
                normq = normq_v[sfm];
                for (i = sfm_start[sfm]; i < sfm_end[sfm]; i++)
                {
                    coeff_out1[i] = coeff_out1[i]*normq;
                }
            }
        }
    }

    return;
}

/*-----------------------------------------------------------------------------
 * floating_point_add:
 *
 * Add two floating point numbers in integer representation: x <- x + y
 * Ported from BASOP code to ensure interoperability
 *----------------------------------------------------------------------------*/
void floating_point_add(
    int   *mx,  /* io: mantissa of the addend Q31 */
    short *ex,  /* io: exponent of the addend Q0  */
    const int    my,  /* i:  mantissa of the adder Q31  */
    const short  ey   /* i:  exponent of the adder Q0   */
)
{
    Word32 accX, accY;
    Word16 align, expo;

    accX = *mx >> 1;
    accY =  my >> 1;
    align = *ex - ey;
    if (align < 0)
    {
        if( align > -32) /* If align < -32, (accY >> (-align) = 0 */
        {
            accX = accX + (accY >> (-align));
        }
    }
    else
    {
        if( align < 32) /* If align > 32, (accX >> align) = 0 */
        {
            accX = accY + (accX >> align);
        }
        else
        {
            accX = accY;
        }
        *ex = ey;
    }

    expo = norm_l(accX); /* aligned to BASOP */
    *mx = accX << expo;
    *ex = *ex + expo - 1;

    return;
}

