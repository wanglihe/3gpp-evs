/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <assert.h>
#include <math.h>
#include "typedef.h"
#include "options.h"
#include "prot.h"
#include "rom_com.h"

/*-------------------------------------------------------------------
 * createFdCngDec()
 *
 * Create an instance of type FD_CNG
 *-------------------------------------------------------------------*/

void createFdCngDec(
    HANDLE_FD_CNG_DEC* hFdCngDec
)
{
    HANDLE_FD_CNG_DEC hs;

    /* Allocate memory */
    hs = (HANDLE_FD_CNG_DEC) calloc(1, sizeof (FD_CNG_DEC));


    createFdCngCom(&(hs->hFdCngCom));

    *hFdCngDec = hs;

    return;
}


/*-------------------------------------------------------------------
 * initFdCngDec()
 *
 * Initialize an instance of type FD_CNG
 *-------------------------------------------------------------------*/

void initFdCngDec(
    HANDLE_FD_CNG_DEC hs,   /* i/o: Contains the variables related to the FD-based CNG process */
    float scale
)
{
    /* Initialize common */
    initFdCngCom( hs->hFdCngCom, scale );

    /* Set some counters and flags */
    hs->flag_dtx_mode           = 0;
    hs->lp_noise                = -20.f;
    hs->lp_speech               = 25.f;

    /* Initialize noise estimation algorithm */
    set_f( hs->bandNoiseShape, 0.0f, FFTLEN2 );
    set_f( hs->partNoiseShape, 0.0f, NPART );
    set_f( hs->msPeriodog, 0.0f, NPART_SHAPING );
    set_f( hs->msAlpha, 0.0f, NPART_SHAPING );
    set_f( hs->msBminWin, 0.0f, NPART_SHAPING );
    set_f( hs->msBminSubWin, 0.0f, NPART_SHAPING );
    set_f( hs->msPsd, 0.0f, NPART_SHAPING );
    set_f( hs->msNoiseFloor, 0.0f, NPART_SHAPING );
    set_f( hs->msNoiseEst, 0.0f, NPART_SHAPING );
    set_f( hs->msMinBuf, FLT_MAX, MSNUMSUBFR*NPART_SHAPING );
    set_f( hs->msCurrentMin, FLT_MAX, NPART_SHAPING );
    set_f( hs->msCurrentMinOut, FLT_MAX, NPART_SHAPING );
    set_f( hs->msCurrentMinSubWindow, FLT_MAX, NPART_SHAPING );
    set_i( hs->msLocalMinFlag, 0, NPART_SHAPING );
    set_i( hs->msNewMinFlag, 0, NPART_SHAPING );
    set_f( hs->msPsdFirstMoment, 0.0f, NPART_SHAPING );
    set_f( hs->msPsdSecondMoment, 0.0f, NPART_SHAPING );
    hs->msPeriodogBufPtr = 0;
    set_f( hs->msPeriodogBuf, 0.0f, MSBUFLEN*NPART_SHAPING );
    set_f( hs->msLogPeriodog, 0.0f, NPART_SHAPING );
    set_f( hs->msLogNoiseEst, 0.0f, NPART_SHAPING );

    return;
}


/*-------------------------------------------------------------------
 * configureFdCngDec()
 *
 * Configure an instance of type FD_CNG
 *-------------------------------------------------------------------*/

void configureFdCngDec(
    HANDLE_FD_CNG_DEC hsDec,   /* i/o: Contains the variables related to the FD-based CNG process */
    short bandwidth,
    int bitrate,
    short L_frame
)
{
    int j, stopBandFR;
    HANDLE_FD_CNG_COM hsCom = hsDec->hFdCngCom;

    hsCom->CngBandwidth = bandwidth;
    if ( hsCom->CngBandwidth == FB )
    {
        hsCom->CngBandwidth = SWB;
    }
    if ( bitrate != FRAME_NO_DATA && bitrate != SID_2k40 )
    {
        hsCom->CngBitrate = bitrate;
    }
    hsCom->numSlots = 16;

    /* NB configuration */
    if ( bandwidth == NB )
    {
        hsCom->FdCngSetup = FdCngSetup_nb;
        hsCom->numCoreBands = 16;
        hsCom->regularStopBand = 16;
    }

    /* WB configuration */
    else if ( bandwidth == WB )
    {
        /* FFT 6.4kHz, no CLDFB */
        if ( hsCom->CngBitrate <= ACELP_8k00 && L_frame==L_FRAME )
        {
            hsCom->FdCngSetup = FdCngSetup_wb1;
            hsCom->numCoreBands = 16;
            hsCom->regularStopBand = 16;
        }
        /* FFT 6.4kHz, CLDFB 8.0kHz */
        else if ( hsCom->CngBitrate <= ACELP_13k20 || L_frame==L_FRAME )
        {
            hsCom->FdCngSetup = FdCngSetup_wb2;
            hsCom->numCoreBands = 16;
            hsCom->regularStopBand = 20;
            if ( L_frame==L_FRAME16k )
            {
                hsCom->FdCngSetup = FdCngSetup_wb2;
                hsCom->numCoreBands = 20;
                hsCom->regularStopBand = 20;
                hsCom->FdCngSetup.fftlen = 640;
                hsCom->FdCngSetup.stopFFTbin = 256;
            }
        }
        /* FFT 8.0kHz, no CLDFB */
        else
        {
            hsCom->FdCngSetup = FdCngSetup_wb3;
            hsCom->numCoreBands = 20;
            hsCom->regularStopBand = 20;
        }
    }

    /* SWB/FB configuration */
    else
    {
        /* FFT 6.4kHz, CLDFB 14kHz */
        if ( L_frame==L_FRAME )
        {
            hsCom->FdCngSetup = FdCngSetup_swb1;
            hsCom->numCoreBands = 16;
            hsCom->regularStopBand = 35;
        }
        /* FFT 8.0kHz, CLDFB 16kHz */
        else
        {
            hsCom->FdCngSetup = FdCngSetup_swb2;
            hsCom->numCoreBands = 20;
            hsCom->regularStopBand = 40;
        }
    }
    hsCom->fftlen = hsCom->FdCngSetup.fftlen;
    hsCom->stopFFTbin = hsCom->FdCngSetup.stopFFTbin;

    /* Configure the SID quantizer and the Comfort Noise Generator */

    hsCom->startBand = 2;
    hsCom->stopBand = hsCom->FdCngSetup.sidPartitions[hsCom->FdCngSetup.numPartitions-1] + 1;
    initPartitions(hsCom->FdCngSetup.sidPartitions, hsCom->FdCngSetup.numPartitions, hsCom->startBand, hsCom->stopBand, hsCom->part, &hsCom->npart, hsCom->midband, hsCom->psize, hsCom->psize_inv, 0);
    if ( hsCom->stopFFTbin == 160 )
    {
        hsCom->nFFTpart = 17;
    }
    else if ( hsCom->stopFFTbin == 256 )
    {
        hsCom->nFFTpart = 20;
    }
    else
    {
        hsCom->nFFTpart = 21;
    }
    hsCom->nCLDFBpart = hsCom->npart - hsCom->nFFTpart;
    for(j=0; j<hsCom->nCLDFBpart; j++)
    {
        hsCom->CLDFBpart[j] = hsCom->part[j+hsCom->nFFTpart] - (hsCom->stopFFTbin-hsCom->startBand);
        hsCom->CLDFBpsize_inv[j] = hsCom->psize_inv[j+hsCom->nFFTpart];
    }

    stopBandFR = (int)floor( 1000.f/*Hz*/ / 25.f/*Hz/Bin*/ );
    if ( stopBandFR > hsCom->stopFFTbin )
    {
        stopBandFR = hsCom->stopFFTbin;
    }
    initPartitions(hsCom->FdCngSetup.shapingPartitions, hsCom->FdCngSetup.numShapingPartitions,
                   hsCom->startBand, hsCom->stopFFTbin, hsDec->part_shaping, &hsDec->npart_shaping, hsDec->midband_shaping,
                   hsDec->psize_shaping, hsDec->psize_inv_shaping, stopBandFR );

    hsDec->nFFTpart_shaping = hsDec->npart_shaping;

    switch (hsCom->fftlen)
    {
    case 512:
        hsCom->fftSineTab = NULL;
        hsCom->olapWinAna = olapWinAna512;
        hsCom->olapWinSyn = olapWinSyn256;
        break;
    case 640:
        hsCom->fftSineTab = fftSineTab640;
        hsCom->olapWinAna = olapWinAna640;
        hsCom->olapWinSyn = olapWinSyn320;
        break;
    default:
        assert(!"Unsupported FFT length for FD-based CNG");
        break;
    }
    hsCom->frameSize = hsCom->fftlen >> 1;

    return;
}


/*-------------------------------------------------------------------
 * deleteFdCngDec()
 *
 * Delete the instance of type FD_CNG
 *-------------------------------------------------------------------*/

void deleteFdCngDec(
    HANDLE_FD_CNG_DEC * hFdCngDec
)
{

    HANDLE_FD_CNG_DEC hsDec = *hFdCngDec;
    if (hsDec != NULL)
    {
        deleteFdCngCom(&(hsDec->hFdCngCom));
        free(hsDec);
        *hFdCngDec = NULL;
    }

    return;
}


/*-------------------------------------------------------------------
 * ApplyFdCng()
 *
 * Apply the CLDFB-based CNG at the decoder
 *-------------------------------------------------------------------*/

void ApplyFdCng(
    float * timeDomainInput,
    float ** realBuffer,          /* i/o: Real part of the buffer */
    float ** imagBuffer,          /* i/o: Imaginary part of the buffer */
    HANDLE_FD_CNG_DEC st,         /* i/o: FD_CNG structure containing all buffers and variables */
    unsigned char m_frame_type,   /* i  : Type of frame at the decoder side */
    Decoder_State *stdec,
    const int concealWholeFrame,
    short is_music
)
{
    float* cngNoiseLevel = st->hFdCngCom->cngNoiseLevel;
    float* sidNoiseEst = st->hFdCngCom->sidNoiseEst;
    int j,k;
    float factor;
    float lsp_cng[M];

    if( st->hFdCngCom->frame_type_previous == ACTIVE_FRAME )
    {
        st->hFdCngCom->inactive_frame_counter = 0;
    }

    switch (m_frame_type)
    {

    case ACTIVE_FRAME:
        /**************************
        * ACTIVE_FRAME at DECODER *
        **************************/

        st->hFdCngCom->inactive_frame_counter = 0;
        st->hFdCngCom->sid_frame_counter = 0;
        /* set noise estimation inactive during concealment, as no update with noise generated by concealment should be performed. */
        /* set noise estimation inactive when we have bit errors, as no update with noise generated by corrupt frame (biterror) should be performed. */
        if ( concealWholeFrame==0 &&
                *timeDomainInput<FLT_MAX &&
                *timeDomainInput>(-FLT_MAX) &&
                *(timeDomainInput+st->hFdCngCom->frameSize-1)<FLT_MAX &&
                *(timeDomainInput+st->hFdCngCom->frameSize-1)>(-FLT_MAX)
                &&
                !(!st->flag_dtx_mode&&stdec->VAD)
                && !(stdec->cng_type==LP_CNG && st->flag_dtx_mode)
                && ( is_music == 0 )
                && (!stdec->BER_detect)
           )
        {
            /* Perform noise estimation at the decoder */
            perform_noise_estimation_dec(timeDomainInput, st );

            /* Update the shaping parameters */
            scalebands(st->msNoiseEst, st->part_shaping, st->nFFTpart_shaping, st->midband_shaping, st->nFFTpart_shaping, st->hFdCngCom->stopFFTbin-st->hFdCngCom->startBand, st->bandNoiseShape, 1);

            /* Update CNG levels */
            if ( st->flag_dtx_mode && stdec->cng_type == FD_CNG )
            {
                bandcombinepow(st->bandNoiseShape, st->hFdCngCom->stopFFTbin-st->hFdCngCom->startBand, st->hFdCngCom->part, st->hFdCngCom->nFFTpart, st->hFdCngCom->psize_inv, st->partNoiseShape); /* This needs to be done only once per inactive phase */

                j = 0;
                for(k=0 ; k<st->hFdCngCom->nFFTpart ; k++)
                {
                    factor = (st->hFdCngCom->sidNoiseEst[k]+DELTA)/(st->partNoiseShape[k]+DELTA);
                    for(; j<=st->hFdCngCom->part[k] ; j++)
                    {
                        cngNoiseLevel[j] = st->bandNoiseShape[j] * factor;
                    }
                }
            }
            else
            {
                /* This sets the new CNG levels until a SID update overwrites it */
                mvr2r(st->bandNoiseShape, cngNoiseLevel, st->hFdCngCom->stopFFTbin-st->hFdCngCom->startBand); /* This sets the new CNG levels until a SID update overwrites it */
            }

            stdec->cngTDLevel = (float)sqrt( (sum_f(cngNoiseLevel, st->hFdCngCom->stopFFTbin - st->hFdCngCom->startBand) / 2 * st->hFdCngCom->fftlen) / stdec->L_frame);
        }
        if ((concealWholeFrame==1) && (stdec->nbLostCmpt==1)
                && sum_f(cngNoiseLevel+st->hFdCngCom->startBand, st->hFdCngCom->stopFFTbin - st->hFdCngCom->startBand)>0.01f)
        {
            /* update lsf cng estimate for concealment. Do that during concealment, in order to avoid addition clean channel complexity*/
            lpc_from_spectrum( cngNoiseLevel, st->hFdCngCom->startBand, st->hFdCngCom->stopFFTbin, st->hFdCngCom->fftlen, st->hFdCngCom->fftSineTab, st->hFdCngCom->A_cng, 0 );

            a2lsp_stab( st->hFdCngCom->A_cng, lsp_cng, stdec->lspold_cng );
            mvr2r( lsp_cng, stdec->lspold_cng, M );

            lsp2lsf( lsp_cng, stdec->lsf_cng, M, stdec->sr_core );
            stdec->plcBackgroundNoiseUpdated = 1;
        }
        break;

    case SID_FRAME:
        st->flag_dtx_mode = 1;

    case ZERO_FRAME:

        if( stdec!=NULL && stdec->cng_type == LP_CNG )
        {
            /* Perform noise estimation on inactive phase at the decoder */
            perform_noise_estimation_dec(timeDomainInput, st );

            /* Update the shaping parameters */
            scalebands(st->msNoiseEst, st->part_shaping, st->nFFTpart_shaping, st->midband_shaping, st->nFFTpart_shaping, st->hFdCngCom->stopFFTbin-st->hFdCngCom->startBand, st->bandNoiseShape, 1);

            /* This sets the new CNG levels until a SID update overwrites it */
            mvr2r(st->bandNoiseShape, cngNoiseLevel, st->hFdCngCom->stopFFTbin-st->hFdCngCom->startBand); /* This sets the new CNG levels until a SID update overwrites it */

            stdec->cngTDLevel = (float)sqrt( (sum_f(cngNoiseLevel, st->hFdCngCom->stopFFTbin - st->hFdCngCom->startBand) / 2 * st->hFdCngCom->fftlen) / stdec->L_frame);
            break;
        }

        st->hFdCngCom->inactive_frame_counter++;

        /*************************************
        * SID_FRAME or ZERO_FRAME at DECODER *
        *************************************/

        /* Detect first non-active frame */
        if (st->hFdCngCom->inactive_frame_counter == 1)
        {
            /* Compute the fine spectral structure of the comfort noise shape using the decoder-side noise estimates */
            bandcombinepow(st->bandNoiseShape, st->hFdCngCom->stopFFTbin-st->hFdCngCom->startBand, st->hFdCngCom->part, st->hFdCngCom->nFFTpart, st->hFdCngCom->psize_inv, st->partNoiseShape);
        }
        if (m_frame_type == SID_FRAME)
        {
            if (st->hFdCngCom->msFrCnt_init_counter < st->hFdCngCom->msFrCnt_init_thresh)
            {
                /* At initialization, interpolate the bin/band-wise levels from the partition levels */
                scalebands(sidNoiseEst, st->hFdCngCom->part, st->hFdCngCom->npart, st->hFdCngCom->midband, st->hFdCngCom->nFFTpart, st->hFdCngCom->stopBand-st->hFdCngCom->startBand, cngNoiseLevel, 1);
            }
            else
            {
                /* Interpolate the CLDFB band levels from the SID (partition) levels */
                if (st->hFdCngCom->regularStopBand>st->hFdCngCom->numCoreBands)
                {
                    scalebands(sidNoiseEst, st->hFdCngCom->part, st->hFdCngCom->npart, st->hFdCngCom->midband, st->hFdCngCom->nFFTpart, st->hFdCngCom->stopBand-st->hFdCngCom->startBand, cngNoiseLevel, 0);
                }
                /* Shape the SID noise levels in each FFT bin */
                j = 0;
                for(k=0 ; k<st->hFdCngCom->nFFTpart ; k++)
                {
                    factor = (sidNoiseEst[k]+DELTA)/(st->partNoiseShape[k]+DELTA);
                    for(; j<=st->hFdCngCom->part[k] ; j++)
                    {
                        cngNoiseLevel[j] = st->bandNoiseShape[j] * factor;
                    }
                }
            }
        }

        if( stdec->codec_mode == MODE2 )
        {
            /* Generate comfort noise during SID or zero frames */
            generate_comfort_noise_dec(realBuffer, imagBuffer, stdec);
        }

        break;

    default:
        break;
    }

    return;
}


/*-------------------------------------------------------------------
 * perform_noise_estimation_dec()
 *
 * Perform noise estimation at the decoder
 *-------------------------------------------------------------------*/

void perform_noise_estimation_dec(
    const float * timeDomainInput,
    HANDLE_FD_CNG_DEC st  /* i/o: FD_CNG structure containing all buffers and variables */
)
{
    float * ptr_r;
    float * ptr_i;
    int   startBand    = st->hFdCngCom->startBand;
    int   stopFFTbin   = st->hFdCngCom->stopFFTbin;
    float * fftBuffer  = st->hFdCngCom->fftBuffer;
    float * periodog   = st->hFdCngCom->periodog;
    float * ptr_per    = periodog;
    float * msPeriodog = st->msPeriodog;
    float * msNoiseEst = st->msNoiseEst;

    float * msLogPeriodog = st->msLogPeriodog;
    float * msLogNoiseEst = st->msLogNoiseEst;

    int   * part       = st->part_shaping;
    int     npart      = st->npart_shaping;
    int     nFFTpart   = st->nFFTpart_shaping;
    float * psize      = st->psize_shaping;
    float * psize_inv  = st->psize_inv_shaping;

    /* Perform STFT analysis */
    AnalysisSTFT(timeDomainInput, fftBuffer, st->hFdCngCom);

    /* Compute the squared magnitude in each FFT bin */
    if( startBand == 0 )
    {
        (*ptr_per) = fftBuffer[0]*fftBuffer[0]; /* DC component */
        ptr_per++;
        ptr_r = fftBuffer + 2;
    }
    else
    {
        ptr_r = fftBuffer + 2*startBand;
    }

    ptr_i = ptr_r+1;

    for( ; ptr_per < periodog+stopFFTbin-startBand ; ptr_per++)
    {
        (*ptr_per) = (*ptr_r)*(*ptr_r) + (*ptr_i)*(*ptr_i);
        ptr_r += 2;
        ptr_i += 2;
    }
    /* Nyquist frequency is discarded */

    /* Rescale to get energy/sample: it should be 2*(1/N)*(2/N), parseval relation with 1/N,*2 for nrg computed till Nyquist only, 2/N as windowed samples correspond to half a frame*/
    v_multc( periodog, 4.f/(float)(st->hFdCngCom->fftlen*st->hFdCngCom->fftlen), periodog, stopFFTbin-startBand);

    /* Adjust to the desired frequency resolution by averaging over spectral partitions for SID transmission */
    bandcombinepow( periodog, stopFFTbin-startBand, part, npart, psize_inv, msPeriodog );

    /* Compress MS inputs */
    compress_range( msPeriodog, msLogPeriodog, npart );

    /* Call the minimum statistics routine for noise estimation */
    minimum_statistics( npart, nFFTpart, psize, msLogPeriodog, st->msNoiseFloor, msLogNoiseEst, st->msAlpha, st->msPsd, st->msPsdFirstMoment, st->msPsdSecondMoment,
                        st->msMinBuf, st->msBminWin, st->msBminSubWin, st->msCurrentMin, st->msCurrentMinOut, st->msCurrentMinSubWindow,
                        st->msLocalMinFlag, st->msNewMinFlag, st->msPeriodogBuf, &(st->msPeriodogBufPtr), st->hFdCngCom );

    /* Expand MS outputs */
    expand_range(msLogNoiseEst, msNoiseEst, npart);

    return;
}


/*-------------------------------------------------------------------
 * FdCng_decodeSID()
 *
 * Decode the FD-CNG bitstream
 *-------------------------------------------------------------------*/

void FdCng_decodeSID(
    Decoder_State *corest             /* i/o: decoder state structure */
)
{
    int N;
    float* sidNoiseEst;
    float gain;
    float preemph_fac;

    int i, index;
    float v[32];
    int indices[32];
    HANDLE_FD_CNG_COM st;

    st = (corest->hFdCngDec)->hFdCngCom;

    sidNoiseEst = st->sidNoiseEst;
    preemph_fac = corest->preemph_fac;

    N = st->npart;
    gain = 0.0f;
    st->sid_frame_counter++;

    /* Read bitstream */
    for ( i=0; i<stages_37bits; i++ )
    {
        indices[i] = get_next_indice(corest, bits_37bits[i]);
    }

    index = get_next_indice( corest, 7 );

    /* MSVQ decoder */
    msvq_dec( cdk_37bits, NULL, NULL, stages_37bits, N, maxN_37bits, indices, v, NULL );

    /* Decode gain */
    gain = ((float)index-60.f)/1.5f;

    /* Apply gain and undo log */
    for ( i=0; i<N; i++ )
    {
        sidNoiseEst[i] = (float)pow( 10.f, (v[i]+gain) / 10.f );
    }

    /* NB last band energy compensation */

    if (st->CngBandwidth == NB)
    {
        sidNoiseEst[N-1] *= NB_LAST_BAND_SCALE;
    }

    if ( st->CngBandwidth == SWB && st->CngBitrate <= ACELP_13k20 )
    {
        sidNoiseEst[N-1] *= SWB_13k2_LAST_BAND_SCALE;
    }

    scalebands(sidNoiseEst, st->part, st->npart, st->midband, st->nFFTpart, st->stopBand-st->startBand, st->cngNoiseLevel, 1);

    lpc_from_spectrum(st->cngNoiseLevel, st->startBand, st->stopFFTbin, st->fftlen, st->fftSineTab, st->A_cng, preemph_fac );

    return;
}


/*-------------------------------------------------------------------
 * noisy_speech_detection()
 *
 *
 *-------------------------------------------------------------------*/

void noisy_speech_detection(
    const short vad,
    const float * ftimeInPtr, /* i  : input time-domain frame                  */
    const int     frameSize,  /* i  : frame size                               */
    const float * msNoiseEst, /* i  : noise estimate over all critical bands   */
    const float * psize,      /* i  : partition sizes                          */
    const int     nFFTpart,   /* i  : Number of partitions taken into account  */
    float *lp_noise,          /* i/o: long term total Noise energy average     */
    float *lp_speech,         /* i/o: long term active speech energy average   */
    short *flag_noisy_speech
)
{
    float tmp;

    if( vad == 0 )
    {
        tmp = dotp(msNoiseEst, psize, nFFTpart);
        *lp_noise = 0.995f * *lp_noise + 0.005f * 10.f*(float)log10( tmp + DELTA);
    }
    else
    {
        tmp = dotp(ftimeInPtr, ftimeInPtr, frameSize) * 2.f/frameSize;
        *lp_speech = 0.995f * *lp_speech + 0.005f * 10.f*(float)log10( tmp + DELTA );
    }

    tmp = *lp_speech - 45.f;
    if (*lp_noise<tmp)
    {
        *lp_noise = tmp;
    }

    *flag_noisy_speech = (*lp_speech-*lp_noise) < 28.f;

    return;
}


/*-------------------------------------------------------------------
 * generate_comfort_noise_dec()
 *
 * Generate the comfort noise based on the target noise level
 *-------------------------------------------------------------------*/

void generate_comfort_noise_dec(
    float ** bufferReal,        /* o  : Real part of input bands      */
    float ** bufferImag,        /* o  : Imaginary part of input bands */
    Decoder_State *stdec
)
{
    short i;
    short j;
    float * ptr_r;
    float * ptr_i;
    HANDLE_FD_CNG_DEC std = stdec->hFdCngDec;
    HANDLE_FD_CNG_COM st = std->hFdCngCom;
    float * cngNoiseLevel = st->cngNoiseLevel;
    float * ptr_level = cngNoiseLevel;
    short * seed = &(st->seed);
    float scale = 1.f;
    float scaleCldfb = CLDFB_SCALING / st->scalingFactor;
    float * fftBuffer = st->fftBuffer;
    float * timeDomainOutput = st->timeDomainBuffer;
    float preemph_fac = stdec->preemph_fac;
    int tcx_transition = 0;
    float enr, att;
    short bwidth = stdec->bwidth;
    short CNG_mode = stdec->CNG_mode;

    /*
      Generate Gaussian random noise in real and imaginary parts of the FFT bins
      Amplitudes are adjusted to the estimated noise level cngNoiseLevel in each bin
    */
    if (st->startBand==0)
    {
        rand_gauss(&fftBuffer[0], seed);
        fftBuffer[0] *= (float)sqrt(scale **ptr_level); /* DC component in FFT */
        ptr_level++;
        ptr_r = fftBuffer + 2;
    }
    else
    {
        fftBuffer[0] = 0.f;
        set_f( fftBuffer+2, 0.0f, 2*(st->startBand-1) );
        ptr_r = fftBuffer + 2*st->startBand;
    }
    ptr_i = ptr_r + 1;
    for( ; ptr_level < cngNoiseLevel+st->stopFFTbin-st->startBand ; ptr_level++)
    {
        /* Real part in FFT bins */
        rand_gauss(ptr_r, seed);
        (*ptr_r) *= (float)sqrt((scale **ptr_level)*0.5f);
        ptr_r += 2;
        /* Imaginary part in FFT bins */
        rand_gauss(ptr_i, seed);
        (*ptr_i) *= (float)sqrt((scale **ptr_level)*0.5f);
        ptr_i += 2;
    }

    /* Remaining FFT bins are set to zero */
    set_f( fftBuffer+2*st->stopFFTbin, 0.0f, st->fftlen-2*st->stopFFTbin);

    /* Nyquist frequency is discarded */
    fftBuffer[1] = 0.f;

    /* If previous frame is active, reset the overlap-add buffer */
    if( st->frame_type_previous==ACTIVE_FRAME )
    {
        set_f( st->olapBufferSynth, 0.0f, st->fftlen);
        if( (stdec->core > 0 && stdec->codec_mode == MODE2) || stdec->codec_mode == MODE1 )
        {
            tcx_transition = 1;
        }
    }

    /* Perform STFT synthesis */
    SynthesisSTFT(fftBuffer, timeDomainOutput, st->olapBufferSynth, st->olapWinSyn, tcx_transition, st );

    /* update CNG excitation energy for LP_CNG */

    /* calculate the residual signal energy */
    enr = dotp( st->exc_cng, st->exc_cng, st->frameSize ) / st->frameSize;

    /* convert log2 of residual signal energy */
    enr = (float)log10( enr + 0.1f ) / (float)log10( 2.0f );

    /* decrease the energy in case of WB input */
    if( bwidth != NB )
    {
        if( bwidth == WB )
        {
            if( CNG_mode >= 0 )
            {
                /* Bitrate adapted attenuation */
                att = ENR_ATT[CNG_mode];
            }
            else
            {
                /* Use least attenuation for higher bitrates */
                att = ENR_ATT[4];
            }
        }
        else
        {
            att = 1.5f;
        }

        enr -= att;
    }

    stdec->lp_ener = 0.8f * stdec->lp_ener + 0.2f * pow( 2.0f, enr );

    /*
      Generate Gaussian random noise in real and imaginary parts of the CLDFB bands
      Amplitudes are adjusted to the estimated noise level cngNoiseLevel in each band
    */
    if (bufferReal!=NULL && st->numCoreBands < st->regularStopBand)
    {
        for(j=st->numCoreBands ; j<st->regularStopBand ; j++)
        {
            for(i=0 ; i<st->numSlots ; i++)
            {
                /* Real part in CLDFB band */
                rand_gauss(&bufferReal[i][j], seed);
                bufferReal[i][j] *= (float)sqrt((scaleCldfb **ptr_level)*0.5f);
                /* Imaginary part in CLDFB band */
                rand_gauss(&bufferImag[i][j], seed);
                bufferImag[i][j] *= (float)sqrt((scaleCldfb **ptr_level)*0.5f);
            }
            ptr_level++;
        }
    }

    /* Overlap-add when previous frame is active */
    if( st->frame_type_previous == ACTIVE_FRAME && stdec->codec_mode == MODE2 )
    {
        float noise[2048], old_exc_ener = 0.f, gain = 0.f, tmp;
        int N = st->frameSize;
        short seed = st->seed;
        float *old_exc, old_Aq[M+1], *old_syn_pe, old_syn;

        if( stdec->core > ACELP_CORE )
        {
            tcx_windowing_synthesis_current_frame(  timeDomainOutput, stdec->tcx_cfg.tcx_mdct_window, /*Keep sine windows for limiting Time modulation*/
                                                    stdec->tcx_cfg.tcx_mdct_window_half, stdec->tcx_cfg.tcx_mdct_window_minimum, stdec->tcx_cfg.tcx_mdct_window_length,
                                                    stdec->tcx_cfg.tcx_mdct_window_half_length, stdec->tcx_cfg.tcx_mdct_window_min_length, 0,
                                                    stdec->tcx_cfg.tcx_last_overlap_mode==ALDO_WINDOW?FULL_OVERLAP:stdec->tcx_cfg.tcx_last_overlap_mode,
                                                    NULL, NULL, NULL, NULL, NULL, N/2, stdec->tcx_cfg.tcx_offset<0?-stdec->tcx_cfg.tcx_offset:0,
                                                    1, 0, 0 );

            if( stdec->tcx_cfg.last_aldo )
            {
                for (i=0; i<(st->frameSize-NS2SA(stdec->sr_core, N_ZERO_MDCT_NS)); i++)
                {
                    timeDomainOutput[i] += stdec->old_outLB[i+NS2SA(stdec->sr_core, N_ZERO_MDCT_NS)];
                }
            }
            else
            {
                tcx_windowing_synthesis_past_frame( stdec->syn_Overl,
                                                    stdec->tcx_cfg.tcx_mdct_window,
                                                    stdec->tcx_cfg.tcx_mdct_window_half,
                                                    stdec->tcx_cfg.tcx_mdct_window_minimum, stdec->tcx_cfg.tcx_mdct_window_length, stdec->tcx_cfg.tcx_mdct_window_half_length,
                                                    stdec->tcx_cfg.tcx_mdct_window_min_length, stdec->tcx_cfg.tcx_last_overlap_mode );

                for (i=0; i<stdec->tcx_cfg.tcx_mdct_window_length; i++)
                {
                    timeDomainOutput[i] += stdec->syn_Overl[i];
                }
            }
        }
        else
        {

            mvr2r(stdec->old_Aq_12_8, old_Aq,M+1);
            old_exc = stdec->old_exc+L_EXC_MEM_DEC-(N/2);
            old_syn_pe = stdec->mem_syn2;
            old_syn = stdec->syn[M];

            for ( i=0; i<N/2; i++ )
            {
                old_exc_ener += old_exc[i]*old_exc[i];
            }

            old_exc_ener = (float)sqrt( old_exc_ener / (float)(N/2) );

            for ( i=0; i<N; i++ )
            {
                rand_gauss(&(noise[i]), &(seed));
                gain += noise[i]*noise[i];
            }

            gain = old_exc_ener / (float)sqrt( gain / (float)N );

            for ( i=0; i<N; i++ )
            {
                noise[i] *= gain;
            }

            syn_filt( old_Aq, M,noise, noise, N, old_syn_pe, 0);

            tmp = old_syn;

            deemph( noise, preemph_fac, N, &tmp );

            for ( i = 0; i < N/2 ; i++ )
            {
                timeDomainOutput[i] += noise[i] * st->olapWinSyn[N/2+i];
            }
        }
    }

    return;
}


/*-------------------------------------------------------------------
 * generate_comfort_noise_dec_hf()
 *
 * Generate the comfort noise based on the target noise level for the CLDFB part
 *-------------------------------------------------------------------*/

void generate_comfort_noise_dec_hf(
    float ** bufferReal,        /* o  : Real part of input bands      */
    float ** bufferImag,        /* o  : Imaginary part of input bands */
    Decoder_State *stdec
)
{
    short i,j;
    float * ptr_level;
    HANDLE_FD_CNG_COM st = stdec->hFdCngDec->hFdCngCom;
    short * seed = &(st->seed);
    float scale = CLDFB_SCALING / st->scalingFactor;

    ptr_level=st->cngNoiseLevel+st->stopFFTbin-st->startBand;
    /*
      Generate Gaussian random noise in real and imaginary parts of the CLDFB bands
      Amplitudes are adjusted to the estimated noise level cngNoiseLevel in each band
    */
    if (st->numCoreBands < st->regularStopBand)
    {
        for(j=st->numCoreBands ; j<st->regularStopBand ; j++)
        {
            for(i=0 ; i<st->numSlots ; i++)
            {
                /* Real part in CLDFB band */
                rand_gauss(&bufferReal[i][j], seed);
                bufferReal[i][j] *= (float)sqrt((scale **ptr_level)*0.5f);
                /* Imaginary part in CLDFB band */
                rand_gauss(&bufferImag[i][j], seed);
                bufferImag[i][j] *= (float)sqrt((scale **ptr_level)*0.5f);
            }
            ptr_level++;
        }
    }

    return;
}


/*-------------------------------------------------------------------
 * generate_masking_noise()
 *
 * Generate additional comfort noise (kind of noise filling)
 *-------------------------------------------------------------------*/

void generate_masking_noise(
    float * timeDomainBuffer,     /* i/o: time-domain signal */
    HANDLE_FD_CNG_COM st,         /* i/o: FD_CNG structure containing all buffers and variables */
    short length,
    short core
)
{
    float * cngNoiseLevel = st->cngNoiseLevel;
    float * fftBuffer     = st->fftBuffer;
    short i;
    float maskingNoise[L_FRAME16k];
    float * ptr_r;
    float * ptr_i;
    float * ptr_level = cngNoiseLevel;
    int startBand = st->startBand;
    short * seed = &(st->seed);
    float scale = 1.f;
    int scaleTableSize;

    if( core != AMR_WB_CORE )
    {
        /* Compute additional CN level */
        scaleTableSize = sizeof (scaleTable_cn_only) / sizeof (scaleTable_cn_only[0]);

        for (i = 0 ; i < scaleTableSize ; i++)
        {
            if ( (st->CngBandwidth == scaleTable_cn_only[i].bwmode) &&
                    (st->CngBitrate >= scaleTable_cn_only[i].bitrateFrom) &&
                    (st->CngBitrate < scaleTable_cn_only[i].bitrateTo) )
            {
                break;
            }
        }

        scale *= (float)pow( 10.f,-scaleTable_cn_only[i].scale/10.f ) - 1.f;
    }
    else
    {
        /* Compute additional CN level */
        scaleTableSize = sizeof (scaleTable_cn_only_amrwbio) / sizeof (scaleTable_cn_only_amrwbio[0]);

        for (i = 0 ; i < scaleTableSize ; i++)
        {
            if (st->CngBitrate >= scaleTable_cn_only_amrwbio[i][0])
            {
                break;
            }
        }

        if(i<scaleTableSize)
        {
            scale *= (float)pow( 10.f,-scaleTable_cn_only_amrwbio[i][1]/10.f ) - 1.f;
        }
        else
        {
            scale=0.f;
        }
    }

    /* Exclude clean speech */
    scale *= st->likelihood_noisy_speech;

    /*
      Generate Gaussian random noise in real and imaginary parts of the FFT bins
      Amplitudes are adjusted to the estimated noise level cngNoiseLevel in each bin
    */
    if (startBand==0)
    {
        rand_gauss(&fftBuffer[0], seed);
        ptr_r = fftBuffer + 2;
        fftBuffer[0] *= (float)sqrt(scale **ptr_level ); /* DC component in FFT */
        ptr_level++;
    }
    else
    {
        fftBuffer[0] = 0.f;
        set_f( fftBuffer+2, 0.0f, 2*(startBand-1) );
        ptr_r = fftBuffer + 2*startBand;
    }
    ptr_i = ptr_r + 1;

    for( ; ptr_level < cngNoiseLevel+st->stopFFTbin-startBand ; ptr_level++)
    {
        /* Real part in FFT bins */
        rand_gauss(ptr_r, seed);
        (*ptr_r) *= (float)sqrt((scale **ptr_level )*0.5f);
        ptr_r += 2;
        /* Imaginary part in FFT bins */
        rand_gauss(ptr_i, seed);
        (*ptr_i) *= (float)sqrt((scale **ptr_level )*0.5f);
        ptr_i += 2;
    }

    /* Remaining FFT bins are set to zero */
    set_f( fftBuffer+2*st->stopFFTbin, 0.0f, st->fftlen-2*st->stopFFTbin);
    /* Nyquist frequency is discarded */
    fftBuffer[1] = 0.f;

    /* Perform STFT synthesis */
    SynthesisSTFT( fftBuffer, maskingNoise, st->olapBufferSynth2, st->olapWinSyn, 0, st );

    /* Add some comfort noise on top of decoded signal */
    v_add(maskingNoise, timeDomainBuffer, timeDomainBuffer, min(st->frameSize,length));

    return;
}


/*-------------------------------------------------------------------
 * generate_masking_noise_update_seed()
 *
 * Update seed for scenarios where generate_masking_noise() is
 * not called based on signal statistics
 *-------------------------------------------------------------------*/

void generate_masking_noise_update_seed(
    HANDLE_FD_CNG_COM st          /* i/o: FD_CNG structure containing all buffers and variables */
)
{
    float * cngNoiseLevel = st->cngNoiseLevel;
    float * ptr_level = cngNoiseLevel;
    int startBand = st->startBand;
    short * seed = &(st->seed);
    float tmp = 0;

    /*
      Generate Gaussian random noise in real and imaginary parts of the FFT bins
      Amplitudes are adjusted to the estimated noise level cngNoiseLevel in each bin
    */
    if (startBand==0)
    {
        rand_gauss(&tmp, seed);
        ptr_level++;
    }

    for( ; ptr_level < cngNoiseLevel+st->stopFFTbin-startBand ; ptr_level++)
    {
        /* Real part in FFT bins */
        rand_gauss(&tmp, seed);
        rand_gauss(&tmp, seed);
    }

    return;
}


/*-------------------------------------------------------------------
 * generate_masking_noise_mdct()
 *
 * Generate additional comfort noise (kind of noise filling)
 *-------------------------------------------------------------------*/

void generate_masking_noise_mdct(
    float * mdctBuffer,         /* i/o: time-domain signal */
    HANDLE_FD_CNG_COM st        /* i/o: FD_CNG structure containing all buffers and variables */
)
{

    float * cngNoiseLevel = st->cngNoiseLevel;
    short i;
    float maskingNoise[2*L_FRAME16k];
    float * ptr_r;

    float * ptr_level = cngNoiseLevel;
    int startBand = st->startBand;
    short * seed = &(st->seed);
    float scale = 1.f;
    int scaleTableSize;

    /* Compute additional CN level */
    scaleTableSize = sizeof (scaleTable_cn_only) / sizeof (scaleTable_cn_only[0]);

    for (i = 0 ; i < scaleTableSize ; i++)
    {
        if ( (st->CngBandwidth == scaleTable_cn_only[i].bwmode) &&
                (st->CngBitrate >= scaleTable_cn_only[i].bitrateFrom) &&
                (st->CngBitrate < scaleTable_cn_only[i].bitrateTo) )
        {
            break;
        }
    }

    scale *= (float)pow( 10.f,-scaleTable_cn_only[i].scale/10.f ) - 1.f;

    /* Exclude clean speech */
    scale *= st->likelihood_noisy_speech;

    /*
      Generate Gaussian random noise in real and imaginary parts of the FFT bins
      Amplitudes are adjusted to the estimated noise level cngNoiseLevel in each bin
    */
    if (startBand==0)
    {
        rand_gauss(&maskingNoise[0], seed);
        maskingNoise[0] *= (float)sqrt(scale **ptr_level * 0.5f ); /* DC component in FFT */
        ptr_level++;
        ptr_r = maskingNoise + 1;
    }
    else
    {
        maskingNoise[0] = 0.f;
        set_f( maskingNoise+1, 0.0f, (startBand-1) );
        ptr_r = maskingNoise + startBand;
    }

    for( ; ptr_level < cngNoiseLevel+st->stopFFTbin-startBand ; ptr_level++)
    {
        /* MDCT bins */
        rand_gauss(ptr_r, seed);
        (*ptr_r) *= (float)sqrt(scale **ptr_level * 0.5f );
        ptr_r += 1;

    }

    /*re-normalization of energy level: M/sqrt(2)*/
    v_multc(maskingNoise, (float)sqrt(NORM_MDCT_FACTOR), maskingNoise,st->stopFFTbin);

    /* Add some comfort noise on top of decoded signal */
    v_add(maskingNoise, mdctBuffer, mdctBuffer, st->stopFFTbin);

    return;
}
