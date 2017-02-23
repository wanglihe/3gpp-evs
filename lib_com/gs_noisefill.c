/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"


/*-------------------------------------------------------------------*
 * gs_noisf()
 *
 * Noise fill-in function
 *-------------------------------------------------------------------*/

static void gs_noisf(
    const short Start_BIN,    /* i  : First bin for noise fill */
    const short NB_Qbins,     /* i  : Number of bin per band   */
    const float Noise_fac,    /* i  : Noise level              */
    const float *y_norm,      /* i  : Quantized pulses         */
    float *exc_diffQ,   /* o  : Quantized pulses with noise added */
    short *seed_tcx,    /* i  : Random generator seed    */
    const short coder_type    /* i  : coder type               */
)
{
    float ftmp;
    short i, k;
    short NB_zer = NB_Qbins/2;

    if( coder_type == INACTIVE )
    {
        NB_zer = 2;
    }

    /*----------------------------------------------*
     * noise fill-in on unquantized subvector       *
     * injected only from 1066Hz to 6400Hz.         *
     *----------------------------------------------*/

    for( k=Start_BIN; k<NB_Qbins + Start_BIN; k+=NB_zer )
    {
        ftmp = 0.0;
        for(i=k; i<k+NB_zer; i++)
        {
            exc_diffQ[i] = (float)y_norm[i];
            ftmp += exc_diffQ[i]*exc_diffQ[i];
        }

        if (ftmp <.5)
        {
            for(i=k; i<k+NB_zer; i++)
            {
                exc_diffQ[i] += Noise_fac*((float)own_random(seed_tcx)/32768.0f);
            }
        }
        else
        {
            /* This is added only to keep the seed in sync between different compilers */
            for(i=k; i<k+NB_zer; i++)
            {
                own_random(seed_tcx);
            }
        }

    }

    return;
}

/*-------------------------------------------------------------------*
 * EstimateNoiseLevel_inner()
 *
 * Estimate noise level from the power spectrum
 *-------------------------------------------------------------------*/

static void EstimateNoiseLevel_inner(
    float *noisepb,     /* o  : Noise per band */
    const long  bitrate,      /* i  : Bitrate of the codec */
    const short i_band,       /* i  : First band to compute the noise  */
    const short Mbands_gn     /* i  : number of bands                  */
)
{
    short i;
    float noise_offset;

    noise_offset = 0.25f;
    if( bitrate > ACELP_24k40 )
    {
        noise_offset = .2f;
    }
    else if( bitrate >= ACELP_22k60 )
    {
        noise_offset = .3f;
    }
    else if( bitrate >= ACELP_9k60 )
    {
        noise_offset = 0.35f;
    }
    else
    {
        noise_offset = .4f;
    }

    set_f( noisepb + i_band, noise_offset, Mbands_gn - i_band );

    for( i = i_band; i < 5; i++ )
    {
        if( noisepb[i] > 0.2f )
        {
            noisepb[i] = 0.2f;
        }
    }

    return;
}


/*-------------------------------------------------------------------*
 * EstimateNoiseLevel()
 *
 *
 *-------------------------------------------------------------------*/

static void EstimateNoiseLevel(
    float *noisepb,             /* o  : Noise per band                          */
    const long  bitrate,              /* i  : Bitrate of the codec                    */
    const short Diff_len,             /* i  : number of bin before cut-off frequency  */
    const short Mbands_gn,            /* i  : number of bands                         */
    const short coder_type,           /* i  : coder type                              */
    const short noise_lev,            /* i  : pulses dynamic                          */
    const short pit_band_idx,         /* i  : bin position of the cut-off frequency   */
    short last_bin,             /* i  : the last bin of bit allocation          */
    short bwidth
)
{
    short i_band;

    i_band = 0;

    if( Diff_len < L_FRAME )
    {
        EstimateNoiseLevel_inner(noisepb, bitrate, i_band, Mbands_gn);

        if( coder_type != INACTIVE )
        {
            if( (bitrate == ACELP_8k00 && last_bin > 8) && bwidth != NB)
            {
                while(Mbands_gn > i_band )
                {
                    noisepb[i_band] *= 2.0f;
                    i_band++;
                }
            }
            else
            {
                while(pit_band_idx > i_band )
                {
                    noisepb[i_band] /= 2.0f;
                    i_band++;
                }
            }
        }
    }

    if ( (coder_type == INACTIVE || noise_lev >= NOISE_LEVEL_SP3) )
    {
        for( i_band = 9; i_band < Mbands_gn; i_band++ )
        {
            noisepb[i_band] *= 1.15f;
        }
    }

    return;
}


/*-------------------------------------------------------------------*
 * Apply_NoiseFill()
 *
 *
 *-------------------------------------------------------------------*/

static void Apply_NoiseFill(
    float *exc_diffQ,             /* i/o: Noise per band                          */
    short *seed_tcx,              /* i  : Seed for noise                          */
    const float *noisepb,               /* i  : Noise per band                          */
    const short Diff_len,               /* i  : number of bin before cut-off frequency  */
    const short Mbands_gn,              /* i  : number of bands                         */
    const short coder_type,             /* i  : coder type                              */
    const short *freq_nsbin_per_band    /* i  : bin per bands tables                    */
)
{
    short StartBin, NB_Qbins, i_band;
    StartBin = 0;
    NB_Qbins  = 0;

    for( i_band = 0; i_band < Mbands_gn; i_band++ )
    {
        StartBin += NB_Qbins;
        NB_Qbins = freq_nsbin_per_band[i_band];

        if( Diff_len < L_FRAME )
        {
            gs_noisf( StartBin, NB_Qbins, noisepb[i_band], exc_diffQ, exc_diffQ, seed_tcx, coder_type );
        }
    }

    return;
}


/*-------------------------------------------------------------------*
 * freq_dnw_scaling()
 *
 *
 *-------------------------------------------------------------------*/

void freq_dnw_scaling(
    const short cor_strong_limit, /* i  : HF correlation                */
    const short coder_type,       /* i  : coder type                    */
    const short noise_lev,        /* i  : Noise level                   */
    const long  core_brate,       /* i  : Core bitrate                  */
    float fy_norm[]         /* i/o: Frequency quantized parameter */
)
{
    float sc_dyn;
    short start_sc, i;

    sc_dyn = 1.0f;
    start_sc = L_FRAME;

    if( core_brate <= ACELP_8k00 && coder_type == INACTIVE )
    {
        sc_dyn *= .15f;
        start_sc = 64;
    }
    else if( coder_type == INACTIVE )
    {
        sc_dyn *= .25f;
        start_sc = 80;
    }
    else
    {
        sc_dyn = (float)(NOISE_LEVEL_SP3 - noise_lev)/10.0f + 0.4f;
        start_sc = 112 + (NOISE_LEVEL_SP3 - noise_lev) * 16;

        if( noise_lev == NOISE_LEVEL_SP0 )
        {
            start_sc = L_FRAME;
        }
    }

    for(i = start_sc; i < L_FRAME; i++)
    {
        fy_norm[i] *= sc_dyn;
    }

    if( (core_brate < ACELP_13k20 && cor_strong_limit  == 0) || core_brate < ACELP_9k60)
    {
        for(i = 160; i < L_FRAME; i++)
        {
            if( fy_norm[i] > 1.0f )
            {
                fy_norm[i] = 1.0f;
            }

            if( fy_norm[i] < -1.0f )
            {
                fy_norm[i] = -1.0f;
            }
        }
    }
    else if( core_brate < ACELP_22k60 )
    {
        for(i = 160; i < L_FRAME; i++)
        {
            if( fy_norm[i] > 1.5f )
            {
                fy_norm[i] = 1.5f;
            }

            if( fy_norm[i] < -1.5f )
            {
                fy_norm[i] = -1.5f;
            }
        }
    }

    return;
}


/*-------------------------------------------------------------------*
 * Decreas_freqPeak()
 *
 *
 *-------------------------------------------------------------------*/

static void Decreas_freqPeak(
    float *lsf_new,                        /* i  : ISFs at the end of the frame                          */
    float *exc_diffQ,                      /* i/o: frequency coefficients of per band                    */
    float rat                              /* i  : threshold of ratio between consecutive lsf_new_diff   */
)
{
    short i, j, k;
    short last_bin = 0;
    short pos = 0;
    float *src;
    float avrg, max;
    float lsf_new_diff[M];

    for(j=1; j<(M-1); j++)
    {
        lsf_new_diff[j] = lsf_new[j] - lsf_new[j-1];
    }

    avrg = 0.0f;
    /* This is to prevent a possible div by 0 in the '*(src) = (*src > 0) ?...'
       loop. The value of 'max' is not important because it will be mutiplied
       by 'avrg' and the result will be close to 0. The 'fabs(*src)/max'
       div by 0 error will be avoided. */
    max = 0.001f;
    for(i=160; i<L_FRAME; i++)
    {
        if(fabs(exc_diffQ[i]) > max)
        {
            max = (float)fabs(exc_diffQ[i]);
            pos = i;
        }
        avrg += (float)fabs(exc_diffQ[i]);
    }
    avrg /= 96;

    for(i=0; i<(M-1); i++)
    {
        if(lsf_new[i] > 4000)
        {
            last_bin = i;
            break;
        }
    }

    for(i=last_bin; i<14; i++)
    {
        if(lsf_new_diff[i] < rat*lsf_new_diff[i-1])
        {
            src = &exc_diffQ[(i-1)*16];
            for(j=0; j<2; j++)
            {
                for(k=0; k<16; k++)
                {
                    if(fabs(*src) > 2.0f*avrg)
                    {
                        *(src) = (*src > 0) ? (float)(avrg*(2.0f-fabs(*src)/max)) : (float)(-avrg*(2.0f-fabs(*src)/max));
                    }
                    src++;
                }
            }
        }
    }

    if(fabs(exc_diffQ[pos]) == max && max > 4.0f*avrg)
    {
        for(i=pos-1; i<pos+2; i++)
        {
            exc_diffQ[pos] *=0.5f;
        }
    }

    return;
}


/*-------------------------------------------------------------------*
 * envelop_modify()
 *
 *
 *-------------------------------------------------------------------*/

static void envelop_modify(
    float *exc_diffQ,                /* i/o: frequency coefficients of per band      */
    short *seed_tcx,                 /* i  : Seed for noise                          */
    short last_bin,                  /* i  : last bin of bit allocation              */
    float *Ener_per_bd_iQ            /* i  : Quantized energy of targeted vector     */
)
{
    short i, j, end_band;
    float Ener, Ener1, *src;
    float weight = 1.0f;

    end_band = L_FRAME;
    Ener = 0.1f;
    for(i=last_bin*16; i<end_band; i++)
    {
        Ener += exc_diffQ[i] * exc_diffQ[i];
    }
    Ener = (float)sqrt((end_band - last_bin*16)/Ener);

    weight = 0.5f;

    src = &exc_diffQ[16*last_bin];
    for(i=last_bin; i<last_bin+4; i++)
    {
        Ener1 = (float)(0.4f*pow(10, Ener_per_bd_iQ[i+1]));
        for(j=0; j<16; j++)
        {
            *src = Ener1*(weight*(*src)*Ener + (1.0f-weight)*own_random(seed_tcx)/32768.0f);
            src++;
        }
    }

    Ener1 = (float)(0.4f*pow(10, Ener_per_bd_iQ[15]));

    src = &exc_diffQ[224];
    for(j=0; j<32; j++)
    {
        *src = Ener1*(weight*(*src)*Ener + (1.0f-weight)*own_random(seed_tcx)/32768.0f);
        src++;
    }

    return;
}


/*-------------------------------------------------------------------*
 * highband_exc_dct_in()
 *
 *
 *-------------------------------------------------------------------*/

void highband_exc_dct_in(
    const long core_brate,                 /* i  : core bitrate                            */
    const short *mfreq_bindiv_loc,         /* i  : bin per bands tables                    */
    short last_bin,                  /* i  : last bin of bit allocation              */
    short Diff_len,                  /* i  : number of bin before cut-off frequency  */
    short noise_lev,                 /* i  : pulses dynamic                          */
    short pit_band_idx,              /* i  : bin position of the cut-off frequency   */
    float *exc_diffQ,                /* i  : frequency coefficients of per band      */
    short *seed_tcx,                 /* i  : Seed for noise                          */
    float *Ener_per_bd_iQ,           /* i  : Quantized energy of targeted vector     */
    short nb_subfr,                  /* i  : Number of subframe considered           */
    float *exc_dct_in,               /* o  : dct of residual signal                  */
    short last_coder_type,           /* i  : coding type of last frame               */
    short *bitallocation_band,       /* i  : bit allocation flag of each band        */
    float *lsf_new,                  /* i  : LSFs at the end of the frame            */
    float *last_exc_dct_in,          /* i  : dct of residual signal of last frame    */
    float *last_ener,                /* i  : frequency energy  of last frame         */
    short *last_bitallocation_band,  /* i  : bit allocation flag of each band  of last frame   */
    short *bitallocation_exc,        /* i  : flag of decoded coefficients            */
    short bfi,                       /* i  : bad frame indicator                     */
    const short coder_type,                /* i  : coder type                              */
    short bwidth,
    float *exc_wo_nf,                /* o  : temporal excitation (in f domain) without noisefill   */
    const short GSC_noisy_speech
    ,float *lt_ener_per_band           /* i/o: Average per band energy */
)
{
    short i, j;
    short MAX_Bin = 0;
    short last_bin_tmp;
    float noisepb[MBANDS_GN];
    float Ener_per_bd_yQ[MBANDS_GN];
    float *src, *dst, *end;
    float ener = 0.0f;
    short length_bin, bwe_flag = 0;

    for( j=10; j<MBANDS_GN; j++ )
    {
        ener += (float)pow(10, Ener_per_bd_iQ[j]);
    }

    if( core_brate == ACELP_8k00 && bwidth != NB )
    {
        if(last_coder_type != AUDIO)
        {
            *last_ener = ener;
        }

        if((last_bin > 8 || Diff_len != 0) && last_coder_type == AUDIO)
        {
            MAX_Bin = 10;
            bwe_flag = 1;
        }
        else
        {
            MAX_Bin = 15;
        }

        last_bin_tmp = last_bin;
        if(last_bin < MAX_Bin)
        {
            last_bin = MAX_Bin;
        }
        last_bin += 1;
    }
    else
    {
        last_bin = MBANDS_GN;
        last_bin_tmp = last_bin;
    }

    if( bfi )
    {
        set_f( noisepb, 0.4f, MBANDS_GN );
    }
    else
    {
        EstimateNoiseLevel( noisepb, core_brate, Diff_len, last_bin, coder_type, noise_lev, pit_band_idx, last_bin_tmp, bwidth );
    }

    if( exc_wo_nf != NULL )
    {
        mvr2r( exc_diffQ, exc_wo_nf, L_FRAME );
    }

    if( GSC_noisy_speech && !bfi )
    {
        set_f( noisepb, 0.1f, MBANDS_GN );
    }

    Apply_NoiseFill( exc_diffQ, seed_tcx, noisepb, Diff_len, last_bin, coder_type, mfreq_bindiv_loc );

    /*--------------------------------------------------------------------------------------*
     * Quantize average gain
     * Substract Q averaged gain
     * VQ of remaining gain per band
     *--------------------------------------------------------------------------------------*/

    if( core_brate == ACELP_8k00 && bwidth != NB )
    {
        Ener_per_band_comp( exc_diffQ, Ener_per_bd_yQ, last_bin+1, 0 );
    }
    else
    {
        Ener_per_band_comp( exc_diffQ, Ener_per_bd_yQ, MBANDS_GN, 1 );

        if( nb_subfr < 4 )
        {
            for( i = L_FRAME-16; i < L_FRAME; i++ )
            {
                exc_diffQ[i] *= (0.067f * i - 15.0f);
            }
        }
    }
    /*--------------------------------------------------------------------------------------*
     * Check potential energy excitation overshoot
     *--------------------------------------------------------------------------------------*/
    if( bfi )
    {
        if (GSC_noisy_speech == 0 &&  coder_type > UNVOICED) /* Here coder_type == last_coder_type because of the bfi */
        {
            for( i=0; i<last_bin; i++ )
            {
                Ener_per_bd_iQ[i]= min( Ener_per_bd_iQ[i], (lt_ener_per_band[i]-0.0376f)- Ener_per_bd_yQ[i]);
                lt_ener_per_band[i] -= 0.0188f;
            }
            for( ; i<MBANDS_GN; i++ )
            {
                Ener_per_bd_iQ[i]= min( Ener_per_bd_iQ[i], (lt_ener_per_band[i]-0.0376f));
                lt_ener_per_band[i] -= 0.0188f;
            }
        }
        else
        {
            for( i=0; i<last_bin; i++ )
            {
                Ener_per_bd_iQ[i]= min( Ener_per_bd_iQ[i], (lt_ener_per_band[i]+0.3f)- Ener_per_bd_yQ[i]);
                lt_ener_per_band[i] -= 0.0188f;
            }
            for(; i<MBANDS_GN; i++ )
            {
                Ener_per_bd_iQ[i]= min( Ener_per_bd_iQ[i], (lt_ener_per_band[i]+0.3f));
                lt_ener_per_band[i] -= 0.0188f;
            }
        }
    }
    /*--------------------------------------------------------------------------------------*
     * Apply decoded gain onto the difference signal
     *--------------------------------------------------------------------------------------*/

    if( GSC_noisy_speech )
    {
        for( i= 0; i < L_FRAME; i++ )
        {
            exc_diffQ[i] *= 0.9f;
        }
    }
    Comp_and_apply_gain( exc_diffQ, Ener_per_bd_iQ, Ener_per_bd_yQ, last_bin, 0 );

    if( exc_wo_nf != NULL )
    {
        Comp_and_apply_gain( exc_wo_nf, Ener_per_bd_iQ, Ener_per_bd_yQ, last_bin, 1 );
        v_add( exc_dct_in, exc_wo_nf, exc_wo_nf, L_FRAME );
    }

    /*--------------------------------------------------------------------------------------*
     * add the correction layer to the LF bins,
     * and add the quantized pulses or the noise for the higher part of the spectrum
     * (non valuable temporal content already zeroed)
     * DC is Zeroed
     *--------------------------------------------------------------------------------------*/

    v_add( exc_dct_in, exc_diffQ, exc_dct_in, L_FRAME );

    if( core_brate == ACELP_8k00 && bwidth != NB )
    {
        if( bwe_flag == 1 )
        {
            last_bin -= 1;
            src = &exc_diffQ[L_FRAME-1];
            dst = &exc_dct_in[MAX_Bin*16-1];
            end = &exc_diffQ[last_bin*16-1];

            while (src > end)
            {
                *src-- = *dst--;
            }

            if( (bitallocation_exc[0] != 0 || bitallocation_exc[1] != 0) && core_brate == ACELP_8k00 )
            {
                exc_diffQ[160] = 0.0f;
            }

            envelop_modify(exc_diffQ, seed_tcx, last_bin, Ener_per_bd_iQ);

            mvr2r( &exc_diffQ[last_bin*16], &exc_dct_in[last_bin*16], L_FRAME-last_bin*16 );
        }

        if( nb_subfr < 4 )
        {
            for( i = L_FRAME-16; i < L_FRAME; i++ )
            {
                exc_dct_in[i] *= (0.067f*i-15.f);
            }
        }

        if( ener < 2*(*last_ener) && ener > 0.5f*(*last_ener) )
        {
            length_bin = 6;
            if(last_coder_type != AUDIO)
            {
                set_s( last_bitallocation_band, 0, 6 );
                mvr2r( &exc_dct_in[(4+length_bin)*16], &last_exc_dct_in[(4+length_bin)*16], length_bin*16 );
            }

            for(i=4; i<(4+length_bin); i++)
            {
                if( !(bitallocation_band[i] == 0 && last_bitallocation_band[i-4] == 0))
                {
                    src = &exc_dct_in[(i+length_bin)*16];
                    dst = &last_exc_dct_in[(i+length_bin)*16];
                    for(j=0; j<16; j++)
                    {
                        if(fabs(*src) > 3.0f*fabs(*dst))
                        {
                            *src = (*src > 0) ? (float)(0.5f*(*src + fabs(*dst))) : (float)(0.5f*(*src - fabs(*dst)));
                        }
                        else if(fabs(*dst) > 3.0f*fabs(*src))
                        {
                            *src = (*src > 0) ? (float)(0.7f*(*src) + 0.3f*fabs(*dst)) : (float)(0.7f*(*src) - 0.3f*fabs(*dst));
                        }
                        src++;
                        dst++;
                    }
                }
            }
        }

        if( bwe_flag == 1 )
        {
            Decreas_freqPeak( lsf_new, exc_dct_in, 0.3f );
        }
        else
        {
            Decreas_freqPeak( lsf_new, exc_dct_in, 0.5f );
        }
    }

    mvr2r( &exc_dct_in[64], &last_exc_dct_in[64], L_FRAME-64 );
    mvs2s( &bitallocation_band[4], last_bitallocation_band, 6 );
    *last_ener = ener;

    return;
}
