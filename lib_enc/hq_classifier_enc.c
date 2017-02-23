/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"

/*-----------------------------------------------------------------*
 * Local constants
 *-----------------------------------------------------------------*/

#define SHARP_DIST_THRES            22.2f

/*-----------------------------------------------------------------*
 * Local functions
 *-----------------------------------------------------------------*/

static void hvq_classifier( const float *input, short *prev_Npeaks, short *prev_peaks,
                            short *hqswb_clas, short *Npeaks, short *peaks, const long core_brate,
                            const short last_core, float *nf_gains, short *hvq_hangover, float *pe_gains );

/*--------------------------------------------------------------------------*
 * hq_classifier_enc()
 *
 * HQ mode selector (decision_matrix)
 *--------------------------------------------------------------------------*/

short hq_classifier_enc(            /* o  : Consumed bits                   */
    Encoder_State *st,              /* i/o: encoder state structure         */
    const short length,             /* i  : Frame length                    */
    const float *coefs,             /* i  : Spectral coefficients           */
    const short is_transient,       /* i  : Transient flag                  */
    short *Npeaks,            /* o  : Number of identified peaks      */
    short *peaks,             /* o  : Peak indices                    */
    float *pe_gains,          /* o  : Peak gains                      */
    float *nf_gains,          /* o  : Noise-fill gains                */
    short *hqswb_clas         /* o  : HQ class                        */
)
{
    short bits;

    *hqswb_clas = HQ_NORMAL;
    if( is_transient )
    {
        *hqswb_clas = HQ_TRANSIENT;
    }
    else if( length == L_FRAME32k && st->core_brate <= HQ_32k && st->bwidth == st->last_bwidth )
    {
        /* Detect HQ_HARMONIC mode */
        *hqswb_clas = peak_avrg_ratio( st->total_brate, coefs, NUMC_N+96, &st->mode_count, &st->mode_count1 );

        /* Detect HQ_HVQ mode */
        hvq_classifier( coefs, &st->prev_Npeaks, st->prev_peaks, hqswb_clas, Npeaks, peaks, st->core_brate, st->last_core,
                        nf_gains, &st->hvq_hangover, pe_gains );
    }
    else if ( length == L_FRAME48k && st->core_brate <= HQ_32k && st->bwidth == st->last_bwidth )
    {
        /* Detect HQ_HARMONIC mode */
        *hqswb_clas = peak_avrg_ratio( st->total_brate, coefs, NUMC_N+96, &st->mode_count, &st->mode_count1 );
        /* Detect HQ_HVQ mode */
        hvq_classifier( coefs, &st->prev_Npeaks, st->prev_peaks, hqswb_clas, Npeaks, peaks, st->core_brate, st->last_core,
                        nf_gains, &st->hvq_hangover, pe_gains );
    }
    if( length == L_FRAME48k && st->core_brate <= HQ_32k && *hqswb_clas == HQ_NORMAL)
    {
        *hqswb_clas = HQ_GEN_FB;
    }
    /* set the required number of bits for signalling */
    if( length >= L_FRAME32k && st->core_brate <= HQ_32k )
    {
        bits = 2;
    }
    else
    {
        bits = 1;
    }

    /* write signalling info to the bitstream */
    if ( length == L_FRAME48k && st->core_brate <= HQ_32k )
    {
        if ( *hqswb_clas >= HQ_GEN_SWB)
        {
            push_indice( st, IND_HQ_SWB_CLAS, *hqswb_clas - 5, bits );
        }
        else
        {
            push_indice( st, IND_HQ_SWB_CLAS, *hqswb_clas, bits );
        }
    }
    else
    {
        push_indice( st, IND_HQ_SWB_CLAS, *hqswb_clas, bits );
    }

    if ( *hqswb_clas == HQ_NORMAL && length == L_FRAME32k && st->core_brate <= HQ_32k )
    {
        *hqswb_clas = HQ_GEN_SWB;
    }

    return bits;
}

/*--------------------------------------------------------------------------*
 * peak_avrg_ratio()
 *
 * Classify the input signal and decide if it has a harmonic structure
 *--------------------------------------------------------------------------*/
short peak_avrg_ratio(
    const long total_brate,
    const float *input_hi,          /* i  : input signal               */
    const short N,                  /* i  : number of coefficients     */
    short *mode_count,        /* i/o: HQ_HARMONIC mode count     */
    short *mode_count1        /* i/o: HQ_NORMAL mode count       */
)
{
    float mean, peak, sharp;
    short i, j, q, k, k1, hqswb_clas;
    float input_abs[L_FRAME32k];

    for ( i = 96; i < N; i++)
    {
        input_abs[i] = (float) fabs(input_hi[i]);
    }

    hqswb_clas = HQ_NORMAL;

    k = 0;
    k1 = 0;
    q = 96;
    for( i = 3; i < 17; i++ )
    {
        peak = 0.0f;
        mean = EPSILON;
        for(j = 0; j < 32; j++, q++)
        {
            mean += input_abs[q];

            if ( input_abs[q] > peak )
            {
                peak = input_abs[q];
            }
        }

        sharp = 32 * peak / mean;

        if( i < 8 )
        {
            if( sharp > 4.5 )
            {
                k += 1;
            }
        }
        else
        {
            if( sharp > 3.6 && peak > 10 )
            {
                k1 += 1;
            }
        }
    }

    if( k + k1 >= 10 && k1 > 5 )
    {
        if( *mode_count < 8 )
        {
            (*mode_count)++;
        }

        if( *mode_count1 > 0 )
        {
            (*mode_count1)--;
        }
    }
    else
    {
        if( *mode_count > 0 )
        {
            (*mode_count)--;
        }

        if( *mode_count1 < 8 )
        {
            (*mode_count1)++;
        }
    }
    if( (k + k1 >= 5 && k1 > 2 && total_brate == HQ_24k40) || (((k + k1 >= 10 && k1 > 5) || *mode_count >= 5) && *mode_count1 < 5) )
    {
        hqswb_clas = HQ_HARMONIC;
    }

    return hqswb_clas;
}


/*--------------------------------------------------------------------------*
 * hvq_classifier()
 *
 * Classification of spectral content for HQ_HVQ mode
 *--------------------------------------------------------------------------*/

static void hvq_classifier(
    const float *input,             /* i  : input signal                */
    short *prev_Npeaks,       /* i/o: Peak number memory          */
    short *prev_peaks,        /* i/o: Peak indices memory         */
    short *hqswb_clas,        /* i/o: HQ class                    */
    short *Npeaks,            /* o  : Number of peaks             */
    short *peaks,             /* o  : Peak indices                */
    const long  core_brate,         /* i  : Core bit-rate               */
    const short last_core,          /* i  : Last core used              */
    float *nf_gains,          /* o  : Noisefloor gains            */
    short *hvq_hangover,      /* i/o: Mode-switch hangover        */
    float *pe_gains           /* o  : peak gains                  */
)
{
    const float *p_adj;
    float sharp_dist;
    float nf, pe, d, peak, thr_tmp, m;
    float input_abs[L_FRAME32k], thr[L_FRAME16k];
    float pe_mean[HVQ_NSUB_32k], nf_mean[HVQ_NSUB_32k];
    float sharp[HVQ_NSUB_32k];

    short num_sharp_bands, i, j, k, q, peak_th, nsub, pindx, N, offset;
    short num_peak_cands, high, low;
    short peak_cand_idx[HVQ_THRES_BIN_32k], avail_peaks[HVQ_NSUB_32k];

    if( *hqswb_clas == HQ_HARMONIC && last_core != ACELP_CORE && last_core != AMR_WB_CORE )
    {
        set_f(thr, 0.0f, L_FRAME16k);
        if( core_brate == HQ_24k40 )
        {
            nsub = HVQ_NSUB_24k;
        }
        else
        {
            nsub = HVQ_NSUB_32k;
        }

        N = nsub * HVQ_BW;

        for( i = 0; i < N; i++ )
        {
            input_abs[i] = (float) fabs(input[i]);
        }

        *Npeaks = 0;
        nf = 800;
        pe = 800;
        num_sharp_bands = 0;
        k = 0;
        q = 0;
        sharp_dist = 0;

        /* Find peak threshold */
        for( i = 0; i < nsub; i++ )
        {
            peak = 0.0f;
            nf_mean[i] = EPSILON;
            pe_mean[i] = EPSILON;
            for( j = 0; j < HVQ_BW; j++, q++ )
            {
                d = input_abs[q];

                if( d > nf )
                {
                    nf = HVQ_NF_WEIGHT1 * nf + (1 - HVQ_NF_WEIGHT1) * d;
                }
                else
                {
                    nf = HVQ_NF_WEIGHT2 * nf + (1 - HVQ_NF_WEIGHT2) * d;
                }

                if( d > pe )
                {
                    pe = HVQ_PE_WEIGHT1 * pe + (1 - HVQ_PE_WEIGHT1) * d;
                }
                else
                {
                    pe = HVQ_PE_WEIGHT2 * pe + (1 - HVQ_PE_WEIGHT2) * d;
                }

                nf_mean[i] += nf;
                pe_mean[i] += pe;

                if( d > peak )
                {
                    peak = d;
                }
            }

            nf_mean[i] /= HVQ_BW;
            pe_mean[i] /= HVQ_BW;


            thr_tmp = (float)pow( pe_mean[i]/nf_mean[i], HVQ_THR_POW ) * nf_mean[i];
            set_f( &thr[k], thr_tmp, HVQ_BW );
            k += HVQ_BW;

            sharp[i] = peak / nf_mean[i];
            sharp_dist += sharp[i] - HVQ_SHARP_THRES;

            if( sharp[i] > HVQ_SHARP_THRES )
            {
                num_sharp_bands++;
            }
        }

        /* Estimate noise floor gains */
        offset = nsub % 2;
        for( i = 0; i < 2*(nsub/2); i++ )
        {
            nf_gains[(2*i+1)/nsub] += nf_mean[i+offset];
            pe_gains[(2*i+1)/nsub] += pe_mean[i+offset];
        }

        for( i = 0; i < HVQ_NF_GROUPS; i++ )
        {
            nf_gains[i] /= nsub/HVQ_NF_GROUPS;
            pe_gains[i] /= nsub/HVQ_NF_GROUPS;
        }

        /* Allocate available peaks */
        for( i = 0; i < nsub; i++ )
        {
            avail_peaks[i] = HVQ_PA_PEAKS_SHARP1;
            if( nf_mean[i] < nf_gains[(2*i+1)/nsub] * HVQ_PA_FAC )
            {
                if( sharp[i] < HVQ_PA_SHARP_THRES3 )
                {
                    avail_peaks[i] = HVQ_PA_PEAKS_SHARP3;
                }
                else if( sharp[i] < HVQ_PA_SHARP_THRES2 )
                {
                    avail_peaks[i] = HVQ_PA_PEAKS_SHARP2;
                }
            }
        }

        /* Adjust threshold around previous peaks */
        for( i = 0; i < *prev_Npeaks; i++ )
        {
            j = prev_peaks[i] - 2;
            k = prev_peaks[i] + 2;
            p_adj = hvq_thr_adj;

            for( q = j; q < k; q++ )
            {
                thr[q] *= *p_adj++;
            }
        }

        num_peak_cands = 0;

        /* Remove everything below threshold for peak search */
        input_abs[0] = 0;
        input_abs[1] = 0;
        input_abs[N-2] = 0;
        input_abs[N-1] = 0;
        for( i = 0; i < N-2; i++ )
        {
            if( input_abs[i] < thr[i] )
            {
                input_abs[i] = 0;
            }
            else
            {
                input_abs[num_peak_cands] = input_abs[i];
                peak_cand_idx[num_peak_cands] = i;
                num_peak_cands++;
            }
        }

        if( core_brate == HQ_24k40 )
        {
            peak_th = HVQ_MAX_PEAKS_24k_CLAS;
        }
        else
        {
            peak_th = HVQ_MAX_PEAKS_32k;
        }

        /* Find peaks */
        pindx = maximum( input_abs, num_peak_cands, &m );
        i = 0;

        while( m > 0 && i < peak_th+1)
        {
            if( avail_peaks[peak_cand_idx[pindx]/HVQ_BW] > 0 )
            {
                peaks[i++] = peak_cand_idx[pindx];
                avail_peaks[peak_cand_idx[pindx]/HVQ_BW]--;
            }

            j = pindx - 2;
            k = pindx + 2;

            if( j < 0 )
            {
                j = 0;
            }

            if( k > num_peak_cands-1 )
            {
                k = num_peak_cands-1;
            }

            low = peak_cand_idx[pindx] - 2;
            high = peak_cand_idx[pindx] + 2;

            if( low < 0 )
            {
                low = 0;
            }

            if( high > N-1 )
            {
                high = N-1;
            }

            for( q = j; q <= pindx; q++ )
            {
                if( peak_cand_idx[q] >= low )
                {
                    peak_cand_idx[q] = 0;
                    input_abs[q] = 0;
                }
            }

            for( q = pindx + 1; q <= k; q++ )
            {
                if( peak_cand_idx[q] <= high )
                {
                    peak_cand_idx[q] = 0;
                    input_abs[q] = 0;
                }
            }

            pindx = maximum(input_abs, num_peak_cands, &m);
        }

        *Npeaks = i;

        /* decision about HQ_HVQ mode */
        if (*Npeaks > HVQ_MIN_PEAKS)
        {
            if( num_sharp_bands > nsub - 3 && *Npeaks <= peak_th )
            {
                sharp_dist /= nsub;
                if( sharp_dist <= SHARP_DIST_THRES && *hvq_hangover < 0 )
                {
                    (*hvq_hangover)++;
                }
                else
                {
                    *hqswb_clas = HQ_HVQ;
                    *hvq_hangover = 2;
                }

                /* update memory */
                *prev_Npeaks = *Npeaks;
                mvs2s( peaks, prev_peaks, *Npeaks );
            }
            else
            {
                if( *hvq_hangover > 0 )
                {
                    *hqswb_clas = HQ_HVQ;
                    (*hvq_hangover)--;
                }
                else
                {
                    *hvq_hangover = -1;
                }
            }
        }
        else
        {
            *hvq_hangover = -1;
        }


        if( core_brate == HQ_24k40 )
        {
            *Npeaks = min( HVQ_MAX_PEAKS_24k, *Npeaks );
        }
        else
        {
            *Npeaks = min( HVQ_MAX_PEAKS_32k, *Npeaks );
        }
    }
    else
    {
        *prev_Npeaks = 0;
        *hvq_hangover = 0;
    }

    return;
}


