/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <assert.h>
#include <math.h>
#include "options.h"
#include "prot.h"
#include "rom_com.h"



/*-------------------------------------------------------------------
 * Local functions
 *-------------------------------------------------------------------*/

static void mhvals( int d, float * m );


/*-------------------------------------------------------------------
 * createFdCngCom()
 *
 * Create an instance of type FD_CNG_COM
 *-------------------------------------------------------------------*/

void createFdCngCom(
    HANDLE_FD_CNG_COM * hFdCngCom
)
{
    HANDLE_FD_CNG_COM hs;

    /* Allocate memory */
    hs = (HANDLE_FD_CNG_COM) calloc(1, sizeof (FD_CNG_COM));

    *hFdCngCom = hs;

    return;
}


/*-------------------------------------------------------------------
 * initFdCngCom()
 *
 *
 *-------------------------------------------------------------------*/

void initFdCngCom(
    HANDLE_FD_CNG_COM hs,   /* i/o: Contains the variables related to the FD-based CNG process */
    float scale
)
{

    /* Calculate FFT scaling factor */
    hs->scalingFactor = 1 / (scale*scale*8.f);

    /* Initialize the overlap-add */
    set_f( hs->timeDomainBuffer, 0.0f, L_FRAME16k );
    set_f( hs->olapBufferAna, 0.0f, FFTLEN );
    set_f( hs->olapBufferSynth, 0.0f, FFTLEN );
    set_f( hs->olapBufferSynth2, 0.0f, FFTLEN );

    /* Initialize the comfort noise generation */
    set_f( hs->fftBuffer, 0.0f, FFTLEN );
    set_f( hs->cngNoiseLevel, 0.0f, FFTCLDFBLEN );

    /* Initialize quantizer */
    set_f( hs->sidNoiseEst, 0.0f, NPART );
    set_f( hs->A_cng, 0.0f, M+1 );
    hs->A_cng[0] = 1.f;

    /* Set some counters and flags */
    hs->inactive_frame_counter  = 0; /* Either SID or zero frames */
    hs->active_frame_counter    = 0;
    hs->frame_type_previous     = ACTIVE_FRAME;
    hs->flag_noisy_speech       = 0;
    hs->likelihood_noisy_speech = 0.f;

    /* Initialize noise estimation algorithm */
    set_f( hs->periodog, 0.0f, PERIODOGLEN );
    mhvals(MSNUMSUBFR*MSSUBFRLEN, &(hs->msM_win));
    mhvals(MSSUBFRLEN, &(hs->msM_subwin));
    set_f( hs->msPeriodogSum, 0.0f, 2 );
    set_f( hs->msPsdSum, 0.0f, 2 );
    set_f( hs->msSlope, 0.0f, 2 );
    set_f( hs->msQeqInvAv, 0.0f, 2 );
    hs->msFrCnt_init_counter = 0;
    hs->msFrCnt_init_thresh  = 1;
    hs->init_old             = 0;
    hs->offsetflag           = 0;
    hs->msFrCnt              = MSSUBFRLEN;
    hs->msMinBufferPtr       = 0;
    set_f( hs->msAlphaCor, 0.3f, 2 );

    return;
}


/*-------------------------------------------------------------------
 * deleteFdCngCom()
 *
 * Delete an instance of type FD_CNG_COM
 *-------------------------------------------------------------------*/

void deleteFdCngCom(
    HANDLE_FD_CNG_COM * hFdCngCom /* i/o: Contains the variables related to the FD-based CNG process */
)
{
    HANDLE_FD_CNG_COM hsCom = *hFdCngCom;
    if (hsCom != NULL)
    {
        free(hsCom);
        *hFdCngCom = NULL;
    }

    return;
}


/*-------------------------------------------------------------------
 * initPartitions()
 *
 * Initialize the spectral partitioning
 *-------------------------------------------------------------------*/

void initPartitions(
    const int * part_in,
    int npart_in,
    int startBand,
    int stopBand,
    int * part_out,
    int * npart_out,
    int * midband,
    float * psize,
    float * psize_inv,
    int stopBandFR
)
{
    int i, j, len_out;

    if (part_in != NULL)
    {
        if (stopBandFR > startBand)
        {
            len_out = stopBandFR - startBand;                                               /*part_out*/
            for(i = 0 ; i < len_out ; i++)
            {
                part_out[i] = i;
            }
        }
        else
        {
            len_out = 0;
        }                                                                                 /*npart_in,part_out*/
        for(j=0 ; j<npart_in && part_in[j]<stopBand ; j++)
        {
            if (part_in[j] >= stopBandFR && part_in[j] >= startBand)
            {
                part_out[len_out++] = part_in[j] - startBand;
            }
        }
    }
    else
    {
        len_out = stopBand - startBand;                                                   /*part_out*/
        for (i = 0 ; i < len_out ; i++)
        {
            part_out[i] = i;
        }
    }
    *npart_out = len_out;
    getmidbands(part_out, len_out, midband, psize, psize_inv);

    return;
}


/*-------------------------------------------------------------------
 * compress_range()
 *
 * Apply some dynamic range compression based on the log
 *-------------------------------------------------------------------*/

void compress_range(
    float* in,
    float* out,
    int len
)
{
    float* ptrIn = in;
    float* ptrOut = out;
    int i;

    /* out = log2( 1 + in ) */
    for(i = 0 ; i < len ; i++)
    {
        *ptrOut = (float)log10(*ptrIn + 1.f);
        ptrIn++;
        ptrOut++;
    }
    v_multc(out, 1.f/(float)log10(2.f), out, len);

    /* Quantize to simulate a fixed-point representation 6Q9 */
    v_multc(out, CNG_LOG_SCALING, out, len);
    for(ptrOut = out ; ptrOut < out + len ; ptrOut++)
    {
        *ptrOut = (float)((int)(*ptrOut+0.5f));
        if ( *ptrOut == 0.f )
        {
            *ptrOut = 1.f;
        }
    }
    v_multc(out, 1./CNG_LOG_SCALING, out, len);

    return;
}


/*-------------------------------------------------------------------
 * expand_range()
 *
 * Apply some dynamic range expansion to undo the compression
 *-------------------------------------------------------------------*/

void expand_range(
    float* in,
    float* out,
    int len
)
{
    float* ptrIn = in;
    float* ptrOut = out;
    int i;

    /* out = (2^(in) - 1) */
    for(i = 0 ; i < len ; i++)
    {
        *ptrOut = (float)pow(2.f,*ptrIn) - 1.f;
        if ( *ptrOut < 0.0003385080526823181f )
        {
            *ptrOut = 0.0003385080526823181f;
        }
        ptrIn++;
        ptrOut++;
    }

    return;
}


/*-------------------------------------------------------------------
 * minimum_statistics()
 *
 * Noise estimation using Minimum Statistics (MS)
 *-------------------------------------------------------------------*/

void minimum_statistics(
    int len,                     /* i  : Vector length */
    int lenFFT,                   /* i  : Length of the FFT part of the vectors */
    float * psize,
    float * msPeriodog,            /* i  : Periodograms */
    float * msNoiseFloor,
    float * msNoiseEst,            /* o  : Noise estimates */
    float * msAlpha,
    float * msPsd,
    float * msPsdFirstMoment,
    float * msPsdSecondMoment,
    float * msMinBuf,
    float * msBminWin,
    float * msBminSubWin,
    float * msCurrentMin,
    float * msCurrentMinOut,
    float * msCurrentMinSubWindow,
    int   * msLocalMinFlag,
    int   * msNewMinFlag,
    float * msPeriodogBuf,
    int   * msPeriodogBufPtr,
    HANDLE_FD_CNG_COM st           /* i/o: FD_CNG structure containing all buffers and variables */
)
{
    float msM_win          = st->msM_win;
    float msM_subwin       = st->msM_subwin;
    float * msPsdSum       = st->msPsdSum;
    float * msPeriodogSum  = st->msPeriodogSum;
    float slope;
    float * ptr;
    float msAlphaCorAlpha  = MSALPHACORALPHA;
    float msAlphaCorAlpha2 = 1.f-MSALPHACORALPHA;

    short i,j,k;
    float scalar, scalar2, scalar3;
    float snr, BminCorr, QeqInv, QeqInvAv;
    float beta;
    float msAlphaHatMin2;
    int   len2 = MSNUMSUBFR*len;
    int   current_len;
    short start, stop, cnt;
    int   totsize;

    /* No minimum statistics at initialization */
    if (st->msFrCnt_init_counter < st->msFrCnt_init_thresh)
    {
        mvr2r(msPeriodog, msPsd, len);
        mvr2r(msPeriodog, msNoiseFloor, len);
        mvr2r(msPeriodog, msNoiseEst, len);
        mvr2r(msPeriodog, msPsdFirstMoment, len);
        set_f( msPsdSecondMoment, 0.0f, len);
        msPeriodogSum[0] = dotp(msPeriodog, psize, lenFFT);
        msPsdSum[0] = msPeriodogSum[0];
        if (lenFFT<len)
        {
            msPeriodogSum[1] = dotp(msPeriodog+lenFFT, psize+lenFFT, len-lenFFT);
            msPsdSum[1] = msPeriodogSum[1];
        }

        /* Increment frame counter at initialization */
        /* Some frames are sometimes zero at initialization => ignore them */
        if (msPeriodog[0]<st->init_old)
        {
            set_f(msCurrentMinOut, FLT_MAX, len);
            set_f(msCurrentMin, FLT_MAX, len);
            set_f(msMinBuf, FLT_MAX, len2);
            set_f(msCurrentMinSubWindow, FLT_MAX, len);
            st->msFrCnt_init_counter++;
        }
        st->init_old = msPeriodog[0];
    }

    else
    {

        /* Consider the FFT and CLDFB bands separately
           - first iteration for FFT bins,
           - second one for CLDFB bands in SWB mode */
        start = 0;
        stop  = lenFFT;
        totsize = st->stopFFTbin-st->startBand;
        cnt   = 0;                                                                         /*msAlphaCor*/
        while (stop > start)
        {
            current_len = stop-start;

            /* Compute smoothed correction factor for PSD smoothing */
            msPeriodogSum[cnt] = dotp(msPeriodog+start, psize+start, current_len);
            scalar      = msPeriodogSum[cnt]*msPeriodogSum[cnt] + DELTA;
            scalar2     = msPsdSum[cnt] - msPeriodogSum[cnt];
            scalar      = max(scalar / (scalar + scalar2*scalar2) , MSALPHACORMAX);
            st->msAlphaCor[cnt] =   msAlphaCorAlpha*st->msAlphaCor[cnt] + msAlphaCorAlpha2*scalar;

            /* Compute SNR */
            snr = dotp(msNoiseFloor+start, psize+start, current_len);
            snr = (msPsdSum[cnt] + DELTA) / (snr + DELTA);
            snr = (float)pow(snr, MSSNREXP);
            msAlphaHatMin2 = min( MSALPHAHATMIN, snr );
            scalar = MSALPHAMAX * st->msAlphaCor[cnt];                                      /*msAlpha,msPsd,msPeriodog,msNoiseFloor*/
            for(j=start ; j<stop ; j++)
            {
                /* Compute optimal smoothing parameter for PSD estimation */
                scalar2    = msNoiseFloor[j] + DELTA;
                scalar2   *= scalar2;
                scalar3    = msPsd[j] - msNoiseFloor[j];
                msAlpha[j] = max( (scalar*scalar2) / (scalar2 + scalar3*scalar3) , msAlphaHatMin2);

                /* Compute the PSD (smoothed periodogram) in each band */
                msPsd[j] = msAlpha[j] * msPsd[j] + (1.f-msAlpha[j]) * msPeriodog[j];
            }
            msPsdSum[cnt] = dotp(msPsd+start, psize+start, current_len);
            QeqInvAv = 0;
            scalar  = ((float)(MSNUMSUBFR*MSSUBFRLEN)-1.f)*(1.f-msM_win);
            scalar2 = ((float)MSSUBFRLEN-1.f)*(1.f-msM_subwin);                             /*msAlpha,msPsd,msPsdFirstMoment,msPsdSecondMoment,msNoiseFloor,msBminSubWin,msBminWin,psize*/
            for(j=start ; j<stop ; j++)
            {
                /* Compute variance of PSD */
                beta    = min( msAlpha[j]*msAlpha[j], MSBETAMAX );
                scalar3 = msPsd[j]-msPsdFirstMoment[j];
                msPsdFirstMoment[j]  = beta*msPsdFirstMoment[j] + (1.f - beta)*msPsd[j];
                msPsdSecondMoment[j] = beta*msPsdSecondMoment[j] + (1.f - beta)*scalar3*scalar3;
                /* Compute inverse of amount of degrees of freedom */
                QeqInv = min( (msPsdSecondMoment[j] + DELTA) / (2.f * msNoiseFloor[j]*msNoiseFloor[j] + DELTA), MSQEQINVMAX );
                QeqInvAv += QeqInv * psize[j];

                /* Compute bias correction Bmin */
                msBminWin[j]    = 1.f + scalar *QeqInv / (0.5f - msM_win   *QeqInv);
                msBminSubWin[j] = 1.f + scalar2*QeqInv / (0.5f - msM_subwin*QeqInv);
            }
            QeqInvAv /= totsize;
            st->msQeqInvAv[cnt] = QeqInvAv;

            /* New minimum? */
            BminCorr  = 1.f + MSAV * (float)sqrt(QeqInvAv);                                 /*msPsd,msBminWin,msNewMinFlag,msCurrentMin,msCurrentMinSubWindow*/
            for(j=start ; j<stop; j++)
            {
                scalar  = BminCorr * msPsd[j];
                scalar2 = scalar * msBminWin[j];
                if (scalar2 < msCurrentMin[j])
                {
                    msNewMinFlag[j]          = 1;
                    msCurrentMin[j]          = scalar2;
                    msCurrentMinSubWindow[j] = scalar * msBminSubWin[j];
                }
                else
                {
                    msNewMinFlag[j]          = 0;
                }
            }

            /* This is used later to identify local minima */
            if (st->msFrCnt >= MSSUBFRLEN)
            {
                i = 0;
                while (i < 3)
                {
                    if (st->msQeqInvAv[cnt] < msQeqInvAv_thresh[i])
                    {
                        break;
                    }
                    else
                    {
                        i++;
                    }
                }
                st->msSlope[cnt] = msNoiseSlopeMax[i];
            }

            /* Consider the FFT and CLDFB bands separately */
            start = stop;
            stop  = len;
            totsize = st->stopBand - st->stopFFTbin;
            cnt++;
        } /*while (stop > start)*/

        /* Update minimum between sub windows */
        if( st->msFrCnt > 1 && st->msFrCnt < MSSUBFRLEN )
        {
            /*msNewMinFlag,msCurrentMinSubWindow,msCurrentMinOut*/
            for(j=0 ; j<len; j++)
            {
                if (msNewMinFlag[j] > 0)
                {
                    msLocalMinFlag[j] = 1;
                }
                if (msCurrentMinSubWindow[j] < msCurrentMinOut[j])
                {
                    msCurrentMinOut[j] = msCurrentMinSubWindow[j];
                }
            }
            /* Get the current noise floor */
            mvr2r(msCurrentMinOut, msNoiseFloor, len);
        }

        /* sub window complete */
        else
        {
            if (st->msFrCnt >= MSSUBFRLEN)
            {
                /* Collect buffers */
                mvr2r(msCurrentMinSubWindow, msMinBuf+len*st->msMinBufferPtr, len);

                /* Compute minimum among all buffers */
                mvr2r(msMinBuf, msCurrentMinOut, len);
                ptr = msMinBuf + len;
                for (i=1 ; i<MSNUMSUBFR ; i++)
                {
                    /*msCurrentMinOut*/
                    for(j=0 ; j<len; j++)
                    {
                        if (*ptr < msCurrentMinOut[j])
                        {
                            msCurrentMinOut[j] = *ptr;
                        }
                        ptr++;
                    }
                }

                /* Take over local minima */
                slope = st->msSlope[0];                                                       /*msLocalMinFlag,msNewMinFlag,msCurrentMinSubWindow,msCurrentMinOut*/
                for(j=0 ; j<len ; j++)
                {
                    if (j==lenFFT)
                    {
                        slope = st->msSlope[1];
                    }
                    if ( msLocalMinFlag[j] && !msNewMinFlag[j] &&
                            msCurrentMinSubWindow[j]<slope*msCurrentMinOut[j] &&
                            msCurrentMinSubWindow[j]>msCurrentMinOut[j] )
                    {
                        msCurrentMinOut[j] = msCurrentMinSubWindow[j];
                        i = j;
                        for(k=0 ; k<MSNUMSUBFR; k++)
                        {
                            msMinBuf[i] = msCurrentMinOut[j];
                            i += len;
                        }
                    }
                }

                /* Reset */
                set_i  ( msLocalMinFlag, 0, len);
                set_f( msCurrentMin, FLT_MAX, len);

                /* Get the current noise floor */
                mvr2r(msCurrentMinOut, msNoiseFloor, len);
            }
        }

        /* Detect sudden offsets based on the FFT bins (core bandwidth) */
        if (msPsdSum[0]>50.f*msPeriodogSum[0])
        {
            if (st->offsetflag>0)
            {
                mvr2r(msPeriodog, msPsd, len);
                mvr2r(msPeriodog, msCurrentMinOut, len);
                set_f( st->msAlphaCor, 1.0f, cnt);
                set_f( msAlpha, 0.0f, len);
                mvr2r(msPeriodog, msPsdFirstMoment, len);
                set_f( msPsdSecondMoment, 0.0f, len);
                msPsdSum[0] = dotp(msPeriodog, psize, lenFFT);
                if (lenFFT<len)
                {
                    msPsdSum[1] = dotp(msPeriodog+lenFFT, psize+lenFFT, len-lenFFT);
                }
            }
            st->offsetflag = 1;
        }
        else
        {
            st->offsetflag = 0;
        }


        /* Increment frame counter */
        if (st->msFrCnt == MSSUBFRLEN)
        {
            st->msFrCnt = 1;
            st->msMinBufferPtr++;
            if (st->msMinBufferPtr==MSNUMSUBFR)
            {
                st->msMinBufferPtr = 0;
            }
        }
        else
        {
            (st->msFrCnt)++;
        }

        {
            /* Smooth noise estimate during CNG phases */                                   /*msNoiseEst,msNoiseFloor*/
            for(j=0 ; j<len ; j++)
            {
                msNoiseEst[j] = 0.95f*msNoiseEst[j] + 0.05f*msNoiseFloor[j];
            }
        }
    }

    {
        /* Collect buffers */
        mvr2r(msPeriodog, msPeriodogBuf+len*(*msPeriodogBufPtr), len);
        (*msPeriodogBufPtr)++;
        if ((*msPeriodogBufPtr)==MSBUFLEN)
        {
            (*msPeriodogBufPtr) = 0;
        }

        /* Upper limit the noise floors with the averaged input energy */                 /*msNoiseEst*/
        for(j=0 ; j<len ; j++)
        {
            scalar = msPeriodogBuf[j];
            for(i=j+len; i<MSBUFLEN*len; i+=len)
            {
                scalar += msPeriodogBuf[i];
            }                                                                              /*division by a constant = multiplication by its (constant) inverse */
            scalar *= (1.f/MSBUFLEN);
            if (msNoiseEst[j]>scalar)
            {
                msNoiseEst[j] = scalar;
            }
            assert(msNoiseEst[j] >= 0);
        }
    }

    return;
}


/*-------------------------------------------------------------------
 * apply_scale()
 *
 * Apply bitrate-dependent scale
 *-------------------------------------------------------------------*/

void apply_scale(
    float *scale,
    int bandwidth,
    int bitrate
)
{
    int i;
    int scaleTableSize = sizeof (scaleTable) / sizeof (scaleTable[0]);

    for (i = 0 ; i < scaleTableSize ; i++)
    {
        if ( (bandwidth == scaleTable[i].bwmode) &&
                (bitrate >= scaleTable[i].bitrateFrom) &&
                (bitrate < scaleTable[i].bitrateTo) )
        {
            break;
        }
    }

    *scale += scaleTable[i].scale;

    return;
}


/*-------------------------------------------------------------------
 * bandcombinepow()
 *
 * Compute the power for each partition
 *-------------------------------------------------------------------*/

void bandcombinepow(
    float* bandpow,     /* i  : Power for each band                       */
    int    nband,       /* i  : Number of bands                           */
    int*   part,        /* i  : Partition upper boundaries (band indices starting from 0) */
    int    npart,       /* i  : Number of partitions                      */
    float* psize_inv,   /* i  : Inverse partition sizes                   */
    float* partpow      /* o  : Power for each partition                  */
)
{

    int   i, p;
    float temp;

    if (nband == npart)
    {
        mvr2r(bandpow, partpow, nband);
    }
    else
    {
        /* Compute the power in each partition */
        i = 0;                                                                            /*part,partpow,psize_inv*/
        for (p = 0; p < npart; p++)
        {
            /* Arithmetic averaging of power for all bins in partition */
            temp = 0;
            for ( ; i <= part[p]; i++)
            {
                temp += bandpow[i];
            }
            partpow[p] = temp*psize_inv[p];
        }
    }

    return;
}


/*-------------------------------------------------------------------
 * scalebands()
 *
 * Scale partitions (with smoothing)
 *-------------------------------------------------------------------*/

void scalebands(
    float* partpow,           /* i  : Power for each partition */
    int*   part,              /* i  : Partition upper boundaries (band indices starting from 0) */
    int    npart,             /* i  : Number of partitions */
    int*   midband,           /* i  : Central band of each partition */
    int    nFFTpart,          /* i  : Number of FFT partitions */
    int    nband,             /* i  : Number of bands */
    float* bandpow,           /* o  : Power for each band */
    short  flag_fft_en
)
{
    int   i, j=0, nint, startBand, startPart, stopPart;
    float val, delta = 0.f;


    /* Interpolate the bin/band-wise levels from the partition levels */
    if (nband == npart)
    {
        mvr2r(partpow, bandpow, npart);
    }
    else
    {
        startBand = 0;
        startPart = 0;
        stopPart  = nFFTpart;
        while (startBand < nband)
        {
            if (flag_fft_en || startPart>=nFFTpart)
            {

                /* first half partition */
                j = startPart;
                val = partpow[j];
                for (i = startBand; i <= midband[j]; i++)
                {
                    bandpow[i] = val;
                }
                j++;

                delta = 1;
                /* inner partitions */
                for ( ; j < stopPart ; j++)
                {
                    nint = midband[j] - midband[j-1];
                    /* log-linear interpolation */                                                                        /* Only one new LOG needs to be computed per loop iteration */
                    delta = (float)exp( (log(partpow[j]+DELTA) - log(partpow[j-1]+DELTA)) * normReciprocal[nint] );
                    val = partpow[j-1];
                    for ( ; i < midband[j]; i++)
                    {
                        val *= delta;
                        bandpow[i] = val;
                    }
                    bandpow[i++] = partpow[j];
                }
                if ( delta>1.f )
                {
                    delta = 1.f;
                }

                /* last half partition */
                val = partpow[stopPart-1];
                for ( ; i <= part[stopPart-1]; i++)
                {
                    val *= delta;
                    bandpow[i] = val;
                }
            }
            startBand = part[stopPart-1]+1;
            startPart = stopPart;
            stopPart  = npart;
        }
    }

    return;
}


/*-------------------------------------------------------------------
 * getmidbands()
 *
 * Get central band for each partition
 *-------------------------------------------------------------------*/

void getmidbands(
    int*   part,              /* i  : Partition upper boundaries (band indices starting from 0) */
    int    npart,             /* i  : Number of partitions */
    int*   midband,           /* o  : Central band of each partition */
    float* psize,             /* o  : Partition sizes */
    float* psize_inv          /* o  : Inverse of partition sizes */
)
{
    int   j;


    /* first half partition */
    midband[0] = part[0];
    psize[0] = (float)part[0] + 1.f;
    psize_inv[0] = normReciprocal[part[0] + 1];
    /* inner partitions */                                                              /*part,midband,psize_inv*/
    for (j = 1; j < npart; j++)
    {
        midband[j] = (part[j-1]+1 + part[j]) >> 1;
        psize[j]   = (float)(part[j] - part[j-1]);
        psize_inv[j] = normReciprocal[part[j] - part[j-1]];
    }

    return;
}


/*-------------------------------------------------------------------
 * AnalysisSTFT()
 *
 * STFT analysis filterbank
 *-------------------------------------------------------------------*/

void AnalysisSTFT(
    const float *  timeDomainInput,
    float *  fftBuffer,          /* o  : FFT bins */
    HANDLE_FD_CNG_COM st        /* i/o: FD_CNG structure containing all buffers and variables */
)
{
    float *olapBuffer = st->olapBufferAna;
    const float *olapWin = st->olapWinAna;

    /* Shift and cascade for overlap-add */
    mvr2r(olapBuffer+st->frameSize, olapBuffer, st->fftlen-st->frameSize);
    mvr2r(timeDomainInput, olapBuffer+st->fftlen-st->frameSize, st->frameSize);

    /* Window the signal */
    v_mult(olapBuffer, olapWin, fftBuffer, st->fftlen);


    /* Perform FFT */
    RFFTN(fftBuffer, st->fftSineTab, st->fftlen, -1);

    return;
}


/*-------------------------------------------------------------------
 * SynthesisSTFT()
 *
 * STFT synthesis filterbank
 *-------------------------------------------------------------------*/

void SynthesisSTFT(
    float * fftBuffer,                    /* i  : FFT bins */
    float * timeDomainOutput,
    float * olapBuffer,
    const float * olapWin,
    int tcx_transition,
    HANDLE_FD_CNG_COM st                  /* i/o: FD_CNG structure containing all buffers and variables */
)
{
    short i;
    float buf[M+1+320], tmp;

    /* Perform IFFT */
    RFFTN(fftBuffer, st->fftSineTab, st->fftlen, 1);

    /* Perform overlap-add */
    mvr2r(olapBuffer+st->frameSize, olapBuffer, st->frameSize);
    set_f( olapBuffer+st->frameSize, 0.0f, st->frameSize);                       /*olapBuffer, fftBuffer, olapWin*/
    if ( tcx_transition )
    {
        for (i=0 ; i<5*st->frameSize/4 ; i++)
        {
            olapBuffer[i] = fftBuffer[i];
        }
    }
    else
    {
        for (i=st->frameSize/4 ; i<3*st->frameSize/4 ; i++)
        {
            olapBuffer[i] += fftBuffer[i] * olapWin[i-st->frameSize/4];
        }
        for ( ; i<5*st->frameSize/4 ; i++)
        {
            olapBuffer[i] = fftBuffer[i];
        }
    }
    for ( ; i<7*st->frameSize/4 ; i++)
    {
        olapBuffer[i] = fftBuffer[i] * olapWin[i-3*st->frameSize/4];
    }

    for ( ; i<st->fftlen ; i++)
    {
        olapBuffer[i] = 0;
    }

    /* Get time-domain signal */
    v_multc( olapBuffer+st->frameSize/4, (float)(st->fftlen/2), timeDomainOutput, st->frameSize);

    /* Get excitation */
    v_multc( olapBuffer+st->frameSize/4-(M+1), (float)(st->fftlen/2), buf, M+1+st->frameSize );
    tmp = buf[0];
    preemph( buf+1, PREEMPH_FAC, M+st->frameSize, &tmp );
    residu( st->A_cng, M, buf+1+M, st->exc_cng, st->frameSize );

    return;
}


/*-------------------------------------------------------------------
 * mhvals()
 *
 * Compute some values used in the bias correction of the minimum statistics algorithm
 *-------------------------------------------------------------------*/

static void mhvals(
    int d,
    float * m
)
{
    int i, j;
    float qi, qj, q;
    int len = sizeof(d_array)/sizeof(int);

    i = 0;
    for (i=0 ; i<len ; i++)
    {
        if (d<=d_array[i])
        {
            break;
        }
    }
    if (i==len)
    {
        i = len-1;
        j = i;
    }
    else
    {
        j = i-1;
    }
    if (d==d_array[i])
    {
        *m = m_array[i];
    }
    else
    {
        qj = (float)sqrt((float)d_array[i-1]);    /*interpolate using sqrt(d)*/
        qi = (float)sqrt((float)d_array[i]);
        q = (float)sqrt((float)d);
        *m = m_array[i] + (qi*qj/q-qj)*(m_array[j]-m_array[i])/(qi-qj);
    }

    return;
}

/*-------------------------------------------------------------------
 * rand_gauss()
 *
 * Random generator with Gaussian distribution with mean 0 and std 1
 *-------------------------------------------------------------------*/

void rand_gauss(float *x, short *seed)
{

    float temp;
    temp = (float)own_random(seed);
    temp += (float)own_random(seed);
    temp += (float)own_random(seed);
    temp *= OUTMAX_INV;

    *x = temp;

    return;
}


/*-------------------------------------------------------------------
 * lpc_from_spectrum()
 *
 *
 *-------------------------------------------------------------------*/

void lpc_from_spectrum(
    float* powspec,
    int start,
    int stop,
    int fftlen,
    const float *fftSineTab,
    float *A,
    float preemph_fac
)
{
    int i;
    float r[32], nf;
    float fftBuffer[FFTLEN], *ptr, *pti;

    /* Power Spectrum */
    ptr = fftBuffer;
    pti = fftBuffer+1;
    nf = 1e-3f;
    for ( i = 0; i < start; i++ )
    {
        *ptr = nf;
        *pti = 0.f;
        ptr += 2;
        pti += 2;
    }
    for ( ; i < stop; i++ )
    {
        *ptr = max( nf, powspec[i-start] );
        *pti = 0.f;
        ptr += 2;
        pti += 2;
    }
    for ( ; i < fftlen/2; i++ )
    {
        *ptr = nf;
        *pti = 0.f;
        ptr += 2;
        pti += 2;
    }
    fftBuffer[1] = nf;

    /* Pre-emphasis */
    ptr = fftBuffer;
    for ( i = 0; i < fftlen/2; i++ )
    {
        *ptr *= (1.f+preemph_fac*preemph_fac-2.0f*preemph_fac*(float)cos(-2.0f*EVS_PI*(float)i/(float)fftlen));

        ptr += 2;
    }
    fftBuffer[1] *= (1.f+preemph_fac*preemph_fac+2.0f*preemph_fac);

    /* Autocorrelation */
    RFFTN(fftBuffer, fftSineTab, fftlen, 1);
    for ( i = 0; i <= M; i++ )
    {
        r[i] = fftBuffer[i] * (fftlen/2) * (fftlen/2);
    }
    if ( r[0] < 100.f )
    {
        r[0] = 100.f;
    }

    r[0] *= 1.0005f;

    /* LPC */
    lev_dur( A, r, M, NULL );

    return;

}


/*-------------------------------------------------------------------
 * FdCng_exc()
 *
 * Generate FD-CNG as LP excitation
 *-------------------------------------------------------------------*/

void FdCng_exc(
    HANDLE_FD_CNG_COM hs,
    short *CNG_mode,
    short L_frame,
    float *lsp_old,
    short first_CNG,
    float *lspCNG,
    float *Aq,                    /* o:   LPC coeffs */
    float *lsp_new,               /* o:   lsp  */
    float *lsf_new,               /* o:   lsf  */
    float *exc,                   /* o:   LP excitation   */
    float *exc2,                  /* o:   LP excitation   */
    float *bwe_exc                /* o:   LP excitation for BWE */
)
{
    short i;
    *CNG_mode = -1;

    /*Get excitation */
    for(i=0; i<L_frame/L_SUBFR; i++)
    {
        mvr2r( hs->A_cng, Aq+i*(M+1), M+1 );
    }
    a2lsp_stab( Aq, lsp_new, lsp_old );
    if( first_CNG == 0 )
    {
        mvr2r( lsp_old, lspCNG, M );
    }
    for( i=0; i<M; i++ )
    {
        /* AR low-pass filter  */
        lspCNG[i] = CNG_ISF_FACT * lspCNG[i] + (1-CNG_ISF_FACT) * lsp_new[i];
    }
    if (L_frame == L_FRAME16k)
    {
        lsp2lsf( lsp_new, lsf_new, M, INT_FS_16k );
    }
    else
    {
        lsp2lsf( lsp_new, lsf_new, M, INT_FS_12k8 );
    }
    mvr2r( hs->exc_cng, exc, L_frame );
    mvr2r( hs->exc_cng, exc2, L_frame );
    if( L_frame == L_FRAME )
    {
        interp_code_5over2( exc2, bwe_exc, L_frame );
    }
    else
    {
        interp_code_4over2( exc2, bwe_exc, L_frame );
    }

    return;
}
