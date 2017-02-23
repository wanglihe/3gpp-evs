/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <math.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"
#include "stat_com.h"

/*--------------------------------------------------------------------------*
 * GetSubbandCorrIndex2_har()
 *
 * Finds the index of best correlation between highband (*inBuf) and lowband (*predBuf) for current subband of length (fLen)
 *--------------------------------------------------------------------------*/

static short GetSubbandCorrIndex2_har( /* o  : best correlation index                   */
    const float *inBuf,                /* i  : target buffer (i.e., highband signal)    */
    const float *predBuf,              /* i  : prediction buffer (i.e., lowband)        */
    const short fLen,                  /* i  : window length                            */
    const short maxLag,                /* i  : search length                            */
    const GainItem *G_item,            /* i  : nZero nonzero components (value+index)   */
    const short nZero,                 /* i  : number of nonzero components used in     */
    short *prev_frame_bstindx    /* i  : previous frame best Indices              */
)
{
    short bestIdx, i, j;
    float corr, energy, corr_sq;
    float lagCorr_sq, lagEnergy, eOld;
    short   absPos;
    short N1, N2;


    absPos = 0;
    bestIdx = 0;
    lagCorr_sq = 0.0f;
    lagEnergy = 1e30f;

    for( i = 0, energy = 0.0f; i < fLen - 1; i++, predBuf++ )
    {
        energy += *predBuf **predBuf;
    }

    predBuf -= fLen-1;
    eOld = 0.0f;
    N1 = max(0, *prev_frame_bstindx-maxLag/2);
    if( *prev_frame_bstindx < 0 )
    {
        N2 = maxLag-1;
    }
    else
    {
        N2 = min(maxLag-1, *prev_frame_bstindx+maxLag/2);
    }

    predBuf +=N1;

    /* find the best lag */
    for( i = N1; i <= N2; i++ )
    {
        corr = 0.0f;

        /* get the energy, remove the old and update with the new energy index */
        for(j = 0, energy = 0.0f; j < fLen; j++, predBuf++)
        {
            energy += *predBuf **predBuf;
        }

        predBuf -= fLen;

        /* energy to be removed in the next lag */
        eOld = *predBuf;
        eOld *= eOld;

        /* get cross-correlation */
        if( energy )
        {
            for(j = 0, corr = 0.0f; j < nZero; j++)
            {
                corr += inBuf[G_item[j].gainIndex]* predBuf[G_item[j].gainIndex];
            }

            corr_sq = corr*corr;

            if( (double)lagCorr_sq*(double)energy < (double)corr_sq*(double)lagEnergy )
            {
                bestIdx = i;
                lagCorr_sq = corr_sq;
                lagEnergy = energy;
            }
        }

        predBuf++;
        absPos++;
    }

    if( lagCorr_sq == 0.0f && *prev_frame_bstindx < 0 )
    {
        bestIdx = 0;
    }
    else
    {
        if( lagCorr_sq == 0.0f )
        {
            bestIdx = *prev_frame_bstindx;
        }
    }

    *prev_frame_bstindx = bestIdx;

    return bestIdx;
}



/*--------------------------------------------------------------------------*
 * getswbindices_har()
 *
 * Finds the pulse index of best correlation between highband (*yos) and lowband (*y2) for two groups of length sbLen
 *--------------------------------------------------------------------------*/

static void getswbindices_har(
    float *yos,                         /* i  : original input spectrum */
    float  *y2,                         /* i  : decoded spectrum */
    const short nBands_search,                /* i  : number of bands */
    short *lagIndices,                  /* o  : pulse index */
    short *prev_frame_bstindx,          /* i/o: prev frame index */
    const short swb_lowband,                  /* i  : length of the LF spectrum */
    const short *subband_offsets,
    const short *sbWidth,
    const short *subband_search_offset

)
{
    GainItem Nbiggest[(NB_SWB_SUBBANDS_HAR_SEARCH_SB)*N_NBIGGEST_PULSEARCH];
    const float *refBuf,*predBuf;
    short search_offset[NB_SWB_SUBBANDS_HAR_SEARCH_SB];
    short nlags[NB_SWB_SUBBANDS_HAR_SEARCH_SB];
    short j, k, sb;
    short n_nbiggestsearch[NB_SWB_SUBBANDS_HAR];
    float low_freqsgnl[L_FRAME32k];

    /* Initializations */
    predBuf = y2;

    /* Get the number of HF groups for performing Similarity search */
    for( sb = 0; sb < nBands_search; sb++ )
    {
        nlags[sb] = (short)pow(2, bits_lagIndices_mode0_Har[sb]);
    }

    for( sb=0; sb<nBands_search; sb++ )
    {
        j = sb * N_NBIGGEST_PULSEARCH ;
        refBuf = yos+swb_lowband + subband_offsets[sb];

        /* Find NBiggest samples in HF Groups */
        FindNBiggest2_simple( refBuf, Nbiggest+j, sbWidth[sb], &n_nbiggestsearch[sb],N_NBIGGEST_PULSEARCH );
        search_offset[sb] = subband_search_offset[sb];
    }

    /* Similarity Search for the HF spectrum */
    for( sb=0; sb<nBands_search; sb++ )
    {
        k = 0;
        if( sb == 0 )
        {
            /* copy SSmoothed LF Spectrum */
            mvr2r( predBuf + search_offset[sb] - nlags[sb]/2, low_freqsgnl, sbWidth[sb] + nlags[sb] );
        }
        else
        {
            /* copy SSmoothed LF Spectrum */
            for( j = search_offset[sb] + nlags[sb]/2; j > search_offset[sb] - sbWidth[sb] - nlags[sb]/2; j-- )
            {
                low_freqsgnl[k] = predBuf[j];
                k++;
            }
        }
        /* correlation b/w HF spectrum Group1 of length sbLen and decoded LF spectrum */
        lagIndices[sb] = GetSubbandCorrIndex2_har( yos + swb_lowband + subband_offsets[sb], low_freqsgnl,
                         sbWidth[sb], nlags[sb], Nbiggest+(sb*N_NBIGGEST_PULSEARCH),
                         n_nbiggestsearch[sb], &prev_frame_bstindx[sb] );
    }

    return;
}

static short GetSubbandCorrIndex2_pulsestep(
    const float *inBuf,
    const float *predBuf,
    const float *predBufMa,
    const short fLen,
    const short maxLag,
    const GainItem *G_item,
    const short nZero,
    short ssearch_buflim,
    const float *predBuf_ni
)
{
    short bestIdx, i, j;
    float corr, energy, corr_sq;
    float lagCorr_sq, lagEnergy;
    short absPos;
    short ib_flag = 0;


    absPos = 0;
    bestIdx = -1;
    lagCorr_sq = 0.0f;
    lagEnergy = 1e30f;

    /* Get the initial energy for zero lag */
    while( *predBuf == 0.0f && absPos < ssearch_buflim )
    {
        predBuf++;
        absPos++;
        predBuf_ni++;
    }

    if( absPos == ssearch_buflim )
    {
        predBuf--;
        absPos--;
        predBuf_ni--;
    }
    for( i = 0, energy = 0.0f; i < fLen; i++, predBuf_ni++ )
    {
        energy += *predBuf_ni **predBuf_ni;
    }

    predBuf_ni -= fLen;
    lagEnergy = energy;

    /* Find the best lag */
    for( i = 0; i < maxLag; i++ )
    {
        corr = 0.0f;

        /* Get the energy, remove the old and update with the new energy index */
        for(j = 0, energy = 0.0f; j < fLen; j++, predBuf_ni++)
        {
            energy += *predBuf_ni **predBuf_ni;
        }
        predBuf_ni -= fLen;

        /* Get cross-correlation */
        if( energy )
        {
            for(j = 0, corr = 0.0f; j < nZero; j++)
            {
                corr += inBuf[G_item[j].gainIndex]* predBufMa[G_item[j].gainIndex];
            }

            corr_sq = corr*corr;

            if( (lagCorr_sq == 0.0f && corr_sq == 0.0f) || (double)lagCorr_sq*(double)energy < (double)corr_sq*(double)lagEnergy )
            {
                bestIdx = i;
                lagCorr_sq = corr_sq;
                lagEnergy = energy;
            }
        }

        predBuf++;
        absPos++;
        predBuf_ni++;
        while( *predBuf == 0.0f && absPos < ssearch_buflim )
        {
            predBuf++;
            absPos++;
            predBuf_ni++;
        }

        if( absPos >= ssearch_buflim )
        {
            if( bestIdx == -1 )
            {
                ib_flag = 1;
            }

            break;
        }
    }
    if( ib_flag )
    {
        bestIdx = 0;
    }
    return bestIdx;
}


/*--------------------------------------------------------------------------*
 * GetSWBIndices()
 *
 * Finds the subband lags for subband coding. Each lag corresponds to
 * the section that is copied and scaled from the low-frequency band
 * to form the high-frequency subband.
 *
 * The subbands are searched in two steps to reduce the complexity of
 * the search relative to full search. In first step, a most representative
 * subband is selected based on the similarity of the subbands.Full search
 * is performed for the selected subband.Based on the lagIndice value for
 * this subband, the region of highest interest is selected. A partial
 * search is performed for rest of the subbands around this region of the
 * low-frequency range.
 *--------------------------------------------------------------------------*/

static void GetSWBIndices(
    const float *predBuf,               /* i  : low-frequency band                  */
    const float *targetBuf,             /* i  : high-frequency band                 */
    const short nBands_search,          /* i  : number of search subbands           */
    const short *sbWidth,               /* i  : subband lengths                     */
    short *lagIndices,            /* o  : selected lags for subband coding    */
    const short predBufLen,             /* i  : low-frequency band length           */
    GainItem *Nbiggest,              /* o  : most representative region          */
    const short *subband_offsets,       /* o  : N biggest components                */
    float *predBuf_ni             /* i  : low-frequency band filled noise     */
)
{
    const float *refBuf;
    short search_offset[NB_SWB_SUBBANDS];
    short nlags[NB_SWB_SUBBANDS];
    short sbLen;
    short j, sb;
    short n_nbiggestsearch[NB_SWB_SUBBANDS];
    short ssearch_buflim;
    float sspectra_ma[L_FRAME32k];

    /* Initializations */
    for (sb = 0; sb < nBands_search; sb++)
    {
        j = sb * N_NBIGGEST_PULSEARCH;
        sbLen = sbWidth[sb];
        refBuf = targetBuf + subband_offsets[sb];
        FindNBiggest2_simple(refBuf, Nbiggest+j, sbLen, &n_nbiggestsearch[sb],N_NBIGGEST_PULSEARCH);
    }

    /* Selection of most representative subband (full search) */
    for ( sb = 0; sb < nBands_search; sb++ )
    {
        nlags[sb] = (short)pow(2, bits_lagIndices[sb]);
    }

    for ( sb = 0; sb < nBands_search; sb++ )
    {
        search_offset[sb] = subband_search_offsets[sb];
    }
    sspectra_ma[0] = (predBuf_ni[0]+predBuf_ni[1])/2.0f;

    for( sb=1; sb < predBufLen-1; sb++ )
    {
        sspectra_ma[sb] = (predBuf_ni[sb-1] + predBuf_ni[sb] + predBuf_ni[sb+1])/3.0f;
    }

    sspectra_ma[sb] = (predBuf_ni[sb-1]+predBuf_ni[sb])/2.0f;
    /* Partial search for rest of subbands except the last which is fixed */
    for( sb = 0; sb < nBands_search; sb++ )
    {
        sbLen = sbWidth[sb];
        ssearch_buflim = predBufLen - (sbLen + search_offset[sb]);
        lagIndices[sb] = GetSubbandCorrIndex2_pulsestep( targetBuf + subband_offsets[sb], predBuf + search_offset[sb],
                         sspectra_ma + search_offset[sb], sbLen, nlags[sb], Nbiggest+(sb*N_NBIGGEST_PULSEARCH),n_nbiggestsearch[sb],
                         ssearch_buflim, predBuf_ni + search_offset[sb]);
    }

    return;
}



static void gethar_noisegn(
    Encoder_State *st,
    float spectra[],
    float noise_flr[],
    float xSynth_har[],
    const short sbWidth[],
    const short lagIndices[],
    const short bands,
    const short har_bands,
    const short fLenLow,
    const short fLenHigh,
    const short subband_offsets[],
    const short subband_search_offset[],
    short band_start[],
    short band_end[],
    short band_width[],
    float band_energy[],
    float be_tonal[],
    float *sspectra,
    const short har_freq_est2,
    const short pos_max_hfe2,
    short *pul_res,
    GainItem pk_sf[]
)
{
    GainItem get_pk[N_NBIGGEST_SEARCH_LRG_B];
    short n_nbiggestsearch,imin,gqlevs;
    short i;
    float hfspec[L_FRAME32k];
    float g,g1,g2,dmin,d;

    /*Generate HF noise*/
    genhf_noise( noise_flr, xSynth_har, sspectra, bands, har_bands, har_freq_est2, pos_max_hfe2, pul_res, pk_sf, fLenLow,
                 fLenHigh, sbWidth, lagIndices, subband_offsets, subband_search_offset );

    mvr2r(spectra+fLenLow,hfspec,fLenHigh);
    FindNBiggest2_simple(hfspec, get_pk, fLenHigh, &n_nbiggestsearch,N_NBIGGEST_SEARCH_LRG_B);
    for(i=0; i<n_nbiggestsearch; i++)
    {
        hfspec[get_pk[i].gainIndex]=0.0f;
    }

    g1= sum2_f(hfspec, fLenHigh);
    g2= sum2_f(xSynth_har,fLenHigh);
    imin = 0;
    if( g1 != 0.0 && g2 != 0.0 )
    {
        g = (float) log10(sqrt(g1/g2));
        gqlevs = 4;
        dmin = FLT_MAX;
        imin = 0;

        for (i = 0; i < gqlevs; i++)
        {
            d = (float) fabs (g - gain_table[i]);
            if (d < dmin)
            {
                dmin = d;
                imin = i;
            }
        }
    }
    push_indice( st,IND_NOISEG, imin, 2);
    g=(float) pow (10.0f,gain_table[imin]);
    /*Tonal energy estimation*/
    ton_ene_est(xSynth_har,be_tonal,band_energy,band_start,band_end,band_width,fLenLow,fLenHigh,bands,har_bands,g,pk_sf,pul_res);

    return;
}
/*--------------------------------------------------------------------------*
 * EncodeSWBSubbands()
 *
 * Main routine for generic SWB coding. High-frequency subband
 * replicated based on the lowband signal. A search is perform
 * find lowband indices denoting the selected lowband subband.
 *--------------------------------------------------------------------------*/

static void EncodeSWBSubbands(
    Encoder_State *st,                        /* i/o: encoder state structure      */
    float *spectra,                   /* i/o: MDCT domain spectrum                        */
    const short fLenLow,                    /* i  : lowband length                              */
    const short fLenHigh,                   /* i  : highband length                             */
    const short nBands,                     /* i  : number of subbands                          */
    const short nBands_search,              /* i  : number of subbands to be searched for BWE   */
    const short *sbWidth,                   /* i  : subband lengths                             */
    const short *subband_offsets,           /* i  : Subband offset for BWE                      */
    short *lagIndices,                /* o  : lowband index for each subband              */
    float *lagGains,                  /* o  : first gain for each subband                 */
    short BANDS,                      /* i  : noise estimate from WB part                 */
    short *band_start,                /* i  : Number subbands/Frame                       */
    short *band_end,                  /* i  : Band Start of each SB                       */
    float *band_energy,               /* i  : Band end of each SB                         */
    const short *p2a_flags,                 /* i  : BAnd energy of each SB                      */
    const short hqswb_clas,                 /* i  : lowband synthesis                           */
    short *prev_frm_index,            /* i  : clas information                            */
    const short har_bands,                  /* i/o: Index of the previous Frame                 */
    const short *subband_search_offset,     /* i  : Number of harmonic LF bands                 */
    short *prev_frm_hfe2,
    short *prev_stab_hfe2,
    short band_width[],
    const float spectra_ni[],               /* i  : coded MDCT domain spectrum + noise          */
    short *ni_seed                    /* i/o: random seed for search buffer NI            */
)
{
    GainItem Nbiggest[NB_SWB_SUBBANDS*N_NBIGGEST_PULSEARCH];
    int i, k;
    float sspectra[L_FRAME32k];
    float sspectra_ni[L_FRAME32k],sspectra_diff[L_FRAME32k],be_tonal[SWB_HAR_RAN1], xSynth_har[L_FRAME32k];
    float ss_min=1.0f,th_g[NB_SWB_SUBBANDS];
    GainItem pk_sf[(NB_SWB_SUBBANDS)*8];
    short pul_res[NB_SWB_SUBBANDS],cnt;
    short har_freq_est1=0,har_freq_est2=0;
    short  flag_dis = 1;
    short pos_max_hfe2=0;

    set_f(sspectra, 0.0f, fLenLow);
    set_f(sspectra_ni, 0.0f, fLenLow);
    set_f(xSynth_har, 0.0f, L_FRAME32k);
    set_s(pul_res,0,NB_SWB_SUBBANDS);


    if( hqswb_clas == HQ_HARMONIC )
    {
        /* Harmonic Structure analysis */
        pos_max_hfe2 = har_est( spectra, fLenLow, &har_freq_est1, &har_freq_est2, &flag_dis, prev_frm_hfe2, subband_search_offset, sbWidth, prev_stab_hfe2 );
        /* Spectrum normalization for the corecoder */
        noise_extr_corcod( spectra, spectra_ni, sspectra, sspectra_diff, sspectra_ni, fLenLow, st->prev_hqswb_clas, &st->prev_ni_ratio );

        /* Find best indices for each group */
        getswbindices_har(spectra, sspectra_ni,nBands_search, lagIndices, prev_frm_index,fLenLow,subband_offsets, sbWidth , subband_search_offset);
        /* Write the indices into the bitstream */
        for (k = 0; k < nBands_search; k++)
        {
            push_indice( st,IND_LAGINDICES, lagIndices[k], bits_lagIndices_mode0_Har[k]);
        }

        if(flag_dis == 0)
        {
            if(har_freq_est2 != SWB_HAR_RAN1 || har_freq_est2 != *prev_frm_hfe2)
            {
                har_freq_est2 += lagIndices[0];
            }
        }

        /*noise generation*/
        gethar_noisegn( st, spectra, sspectra_diff, xSynth_har, sbWidth, lagIndices, BANDS, har_bands, fLenLow, fLenHigh, subband_offsets, subband_search_offset,
                        band_start, band_end, band_width, band_energy, be_tonal, sspectra, har_freq_est2, pos_max_hfe2, pul_res, pk_sf );

        /*HF Spectrum Generation*/
        Gettonl_scalfact( xSynth_har, spectra_ni, fLenLow, fLenHigh, har_bands, BANDS, band_energy, band_start, band_end, p2a_flags, be_tonal, pk_sf, pul_res );

        if(flag_dis == 0)
        {
            *prev_frm_hfe2 = 0;
        }
        else
        {
            *prev_frm_hfe2 = har_freq_est2;
        }

        for( k = har_bands; k < BANDS; k++ )
        {
            for( i = band_start[k]; i <= band_end[k]; i++ )
            {
                spectra[i] = xSynth_har[i-fLenLow];
            }
        }
    }
    else
    {
        /* Spectrum normalization for the corecoder*/
        ss_min = spectrumsmooth_noiseton(spectra,spectra_ni,sspectra,sspectra_diff,sspectra_ni,fLenLow, ni_seed);

        /* Get lag indices */
        GetSWBIndices( sspectra, spectra + fLenLow, nBands, sbWidth, lagIndices, fLenLow,
                       Nbiggest, subband_offsets, sspectra );
        /* Bitstream operations */
        for (k = 0; k < nBands; k++)
        {
            if( p2a_flags[BANDS-NB_SWB_SUBBANDS+k] == 0 )
            {
                push_indice( st,IND_LAGINDICES,lagIndices[k], bits_lagIndices[k]);
            }
            else
            {
                lagIndices[k] = 0;
            }
        }
        convert_lagIndices_pls2smp( lagIndices, nBands, lagIndices, sspectra, sbWidth, fLenLow );
        /*get levels for missing bands*/
        GetlagGains( sspectra_ni, &band_energy[BANDS-NB_SWB_SUBBANDS], nBands, sbWidth, lagIndices, fLenLow, lagGains );
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
        GetSynthesizedSpecThinOut(sspectra_ni, xSynth_har, nBands, sbWidth, lagIndices, lagGains, fLenLow );
        /*Level adjustment for the missing bands*/
        noiseinj_hf( xSynth_har, th_g, band_energy, st->prev_En_sb, p2a_flags, BANDS, band_start, band_end, fLenLow );

        for( k = BANDS - NB_SWB_SUBBANDS; k < BANDS; k++ )
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

/*--------------------------------------------------------------------------*
 * swb_bwe_enc_lr()
 *
 * Main encoding routine of  SWB BWE for the LR MDCT core
 *--------------------------------------------------------------------------*/

void swb_bwe_enc_lr(
    Encoder_State *st,               /* i/o: encoder state structure      */
    const float m_core[],          /* i  : lowband synthesis                            */
    const float m_orig[],          /* i/o: scaled orig signal (MDCT)                    */
    float m[],               /* o  : highband synthesis with lowband zeroed       */
    const long  total_brate,       /* i  : total bitrate for selecting subband pattern  */
    short BANDS,             /* i  : Total number of Subbands in a frame          */
    short *band_start,       /* i  : band start of each SB                        */
    short *band_end,         /* i  : band end of each SB                          */
    float *band_energy,      /* i  : band_energy of each SB                       */
    short *p2a_flags,        /* i  : HF tonal indicator                           */
    const short hqswb_clas,        /* i  : HQ_NORMAL or HQ_HARMONIC mode                */
    short lowlength,         /* i  : lowband length                               */
    short highlength,        /* i  : highband length                              */
    short *prev_frm_index,   /* i/o: previous frame lag index for harmonic mode   */
    const short har_bands,         /* i  : Number of LF harmonic bands                  */
    short *prev_frm_hfe2,
    short *prev_stab_hfe2,
    short band_width[],      /* i  : subband bandwidths                           */
    const float y2_ni[],           /* i/o: Sparse filled core coder                     */
    short *ni_seed           /* i/o: random seed for search buffer NI             */
)
{
    short k;
    short nBands;
    short nBands_search;
    short wBands[NB_SWB_SUBBANDS];
    short lagIndices[NB_SWB_SUBBANDS];
    float lagGains[NB_SWB_SUBBANDS];
    const short *subband_offsets;
    short swb_lowband, swb_highband;
    const short *subband_search_offset;
    subband_search_offset = subband_search_offsets_13p2kbps_Har;
    subband_offsets = subband_offsets_sub5_13p2kbps_Har;
    hf_parinitiz(total_brate,hqswb_clas,lowlength,highlength,wBands,&subband_search_offset,&subband_offsets,&nBands,&nBands_search,&swb_lowband,&swb_highband);
    /* Prepare m[], low part from WB core, high part from 32k input */
    mvr2r( m_core, m, swb_lowband );
    mvr2r( m_orig + swb_lowband, m + swb_lowband, swb_highband );

    /* SWB BWE encoding */
    EncodeSWBSubbands( st, m, swb_lowband, swb_highband, nBands, nBands_search, wBands, subband_offsets,
                       lagIndices, lagGains, BANDS, band_start, band_end, band_energy, p2a_flags,
                       hqswb_clas, prev_frm_index, har_bands, subband_search_offset, prev_frm_hfe2, prev_stab_hfe2,band_width, y2_ni, ni_seed );
    m[swb_lowband+swb_highband-1] *= 0.0625f;
    m[swb_lowband+swb_highband-2] *= 0.125f;
    m[swb_lowband+swb_highband-3] *= 0.25f;
    m[swb_lowband+swb_highband-4] *= 0.5f;

    /* set frequencies below 6.4 kHz to zero */
    if(hqswb_clas == HQ_NORMAL)
    {
        for (k=0; k < swb_lowband; k++)
        {
            m[k] = 0.0f;
        }
    }

    return;
}

