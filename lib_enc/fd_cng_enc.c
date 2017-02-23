/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <assert.h>
#include <math.h>
#include "rom_enc.h"
#include "rom_com.h"
#include "prot.h"
#include "stat_enc.h"
#include "options.h"

/*-------------------------------------------------------------------*
* createFdCngEnc()
*
*
*-------------------------------------------------------------------*/

void createFdCngEnc(HANDLE_FD_CNG_ENC* hFdCngEnc)
{
    HANDLE_FD_CNG_ENC hs;

    /* Allocate memory */
    hs = (HANDLE_FD_CNG_ENC) calloc(1, sizeof (FD_CNG_ENC));

    createFdCngCom(&(hs->hFdCngCom));

    *hFdCngEnc = hs;

    return;
}

/*-------------------------------------------------------------------*
* initFdCngEnc()
*
* Initialize FD_CNG
*-------------------------------------------------------------------*/

void initFdCngEnc(
    HANDLE_FD_CNG_ENC hsEnc,  /* i/o: Contains the variables related to the FD-based CNG process */
    int input_Fs,             /* i: input signal sampling frequency in Hz */
    float scale               /* i: scaling factor */
)
{
    int j;
    HANDLE_FD_CNG_COM hsCom = hsEnc->hFdCngCom;

    /* Initialize common */

    initFdCngCom( hsCom, scale );

    /* Configure the Noise Estimator */

    hsCom->numSlots          = 16;
    hsCom->numCoreBands      = 16;
    hsCom->regularStopBand   = input_Fs/800;
    if ( hsCom->regularStopBand > 40 )
    {
        hsCom->regularStopBand = 40;
    }

    hsCom->startBand = 2;
    if ( hsCom->regularStopBand == 10 )
    {
        hsCom->stopFFTbin = 160;
        hsCom->stopBand = 160;
        hsCom->nFFTpart = 17;
    }
    else
    {
        hsCom->stopFFTbin = 256;
        hsCom->stopBand = hsCom->regularStopBand - hsCom->numCoreBands + hsCom->stopFFTbin;
        hsCom->nFFTpart = 20;
    }

    initPartitions( sidparts_encoder_noise_est, sizeof(sidparts_encoder_noise_est)/sizeof(int), hsCom->startBand, hsCom->stopBand, hsCom->part, &hsCom->npart, hsCom->midband, hsCom->psize, hsCom->psize_inv, 0);

    hsCom->nCLDFBpart = hsCom->npart - hsCom->nFFTpart;
    for(j=0; j<hsCom->nCLDFBpart; j++)
    {
        hsCom->CLDFBpart[j] = hsCom->part[j+hsCom->nFFTpart] - (256-hsCom->startBand);
        hsCom->CLDFBpsize_inv[j] = hsCom->psize_inv[j+hsCom->nFFTpart];
    }

    /* Initialize the Noise Estimator */

    set_f( hsEnc->msPeriodog, 0.0f, NPART );
    set_f( hsEnc->msAlpha, 0.0f, NPART );
    set_f( hsEnc->msBminWin, 0.0f, NPART );
    set_f( hsEnc->msBminSubWin, 0.0f, NPART );
    set_f( hsEnc->msPsd, 0.0f, NPART );
    set_f( hsEnc->msNoiseFloor, 0.0f, NPART );
    set_f( hsEnc->msNoiseEst, 0.0f, NPART );
    set_f( hsEnc->energy_ho, 0.0f, NPART );
    set_f( hsEnc->msNoiseEst_old, 0.0f, NPART );
    set_f( hsEnc->msMinBuf, FLT_MAX, MSNUMSUBFR*NPART );
    set_f( hsEnc->msCurrentMin, FLT_MAX, NPART );
    set_f( hsEnc->msCurrentMinOut, FLT_MAX, NPART );
    set_f( hsEnc->msCurrentMinSubWindow, FLT_MAX, NPART );
    set_i( hsEnc->msLocalMinFlag, 0, NPART );
    set_i( hsEnc->msNewMinFlag, 0, NPART );
    set_f( hsEnc->msPsdFirstMoment, 0.0f, NPART );
    set_f( hsEnc->msPsdSecondMoment, 0.0f, NPART );
    hsEnc->msPeriodogBufPtr = 0;
    set_f( hsEnc->msPeriodogBuf, 0.0f, MSBUFLEN*NPART );
    set_f( hsEnc->msLogPeriodog, 0.0f, NPART );
    set_f( hsEnc->msLogNoiseEst, 0.0f, NPART );

    return;
}

/*-------------------------------------------------------------------*
* configureFdCngEnc()
*
* Configure FD_CNG
*-------------------------------------------------------------------*/

void configureFdCngEnc(
    HANDLE_FD_CNG_ENC hsEnc,   /* i/o: Contains the variables related to the FD-based CNG process */
    short bandwidth,
    int bitrate
)
{
    HANDLE_FD_CNG_COM hsCom = hsEnc->hFdCngCom;
    float   psizeDec[NPART];
    float   psize_invDec[NPART];

    hsCom->CngBandwidth = bandwidth;
    if ( hsCom->CngBandwidth == FB )
    {
        hsCom->CngBandwidth = SWB;
    }
    hsCom->CngBitrate = bitrate;

    /* NB configuration */
    if ( bandwidth == NB )
    {
        hsCom->FdCngSetup = FdCngSetup_nb;
    }

    /* WB configuration */
    else if ( bandwidth == WB )
    {
        /* FFT 6.4kHz, no CLDFB */
        if ( bitrate <= ACELP_8k00 )
        {
            hsCom->FdCngSetup = FdCngSetup_wb1;
        }
        /* FFT 6.4kHz, CLDFB 8.0kHz */
        else if ( bitrate <= ACELP_13k20 )
        {
            hsCom->FdCngSetup = FdCngSetup_wb2;
        }
        /* FFT 8.0kHz, no CLDFB */
        else
        {
            hsCom->FdCngSetup = FdCngSetup_wb3;
        }
    }

    /* SWB/FB configuration */
    else
    {
        /* FFT 6.4kHz, CLDFB 14kHz */
        if ( bitrate <= ACELP_13k20 )
        {
            hsCom->FdCngSetup = FdCngSetup_swb1;
        }
        /* FFT 8.0kHz, CLDFB 16kHz */
        else
        {
            hsCom->FdCngSetup = FdCngSetup_swb2;
        }
    }
    hsCom->fftlen = hsCom->FdCngSetup.fftlen;
    hsEnc->stopFFTbinDec = hsCom->FdCngSetup.stopFFTbin;

    /* Configure the SID quantizer and the Confort Noise Generator */

    hsEnc->startBandDec = hsCom->startBand;
    hsEnc->stopBandDec = hsCom->FdCngSetup.sidPartitions[hsCom->FdCngSetup.numPartitions-1] + 1;
    initPartitions( hsCom->FdCngSetup.sidPartitions, hsCom->FdCngSetup.numPartitions, hsEnc->startBandDec, hsEnc->stopBandDec,
                    hsEnc->partDec, &hsEnc->npartDec, hsEnc->midbandDec, psizeDec, psize_invDec, 0 );

    if ( hsEnc->stopFFTbinDec == 160 )
    {
        hsEnc->nFFTpartDec = 17;
    }
    else if ( hsEnc->stopFFTbinDec == 256 )
    {
        hsEnc->nFFTpartDec = 20;
    }
    else
    {
        hsEnc->nFFTpartDec = 21;
    }

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


/*-------------------------------------------------------------------*
* deleteFdCngEnc()
*
* Delete the instance of type FD_CNG
*-------------------------------------------------------------------*/

void deleteFdCngEnc(
    HANDLE_FD_CNG_ENC *hFdCngEnc
)
{

    HANDLE_FD_CNG_ENC hsEnc = *hFdCngEnc;
    if (hsEnc != NULL)
    {
        deleteFdCngCom(&(hsEnc->hFdCngCom));
        free(hsEnc);
        *hFdCngEnc = NULL;
    }

    return;
}

/*-------------------------------------------------------------------*
* resetFdCngEnc()
*
* Reset the instance of type FD_CNG
*-------------------------------------------------------------------*/

void resetFdCngEnc(
    Encoder_State * st
)
{
    int n;
    float totalNoiseIncrease;

    /* Detect fast increase of totalNoise */
    totalNoiseIncrease = st->totalNoise - st->last_totalNoise;
    st->last_totalNoise = st->totalNoise;
    if ( totalNoiseIncrease > 0 )
    {
        if ( st->totalNoise_increase_len == TOTALNOISE_HIST_SIZE )
        {
            for ( n = 0; n < TOTALNOISE_HIST_SIZE - 1; n++ )
            {
                st->totalNoise_increase_hist[n] = st->totalNoise_increase_hist[n+1];
            }
            st->totalNoise_increase_hist[TOTALNOISE_HIST_SIZE-1] = totalNoiseIncrease;
        }
        else
        {
            st->totalNoise_increase_hist[st->totalNoise_increase_len] = totalNoiseIncrease;
            st->totalNoise_increase_len++;
        }
    }
    else
    {
        st->totalNoise_increase_len = 0;
    }
    totalNoiseIncrease = 0.f;
    for ( n = 0; n < st->totalNoise_increase_len; n++ )
    {
        totalNoiseIncrease += st->totalNoise_increase_hist[n];
    }
    if (
        (totalNoiseIncrease > 5 && st->totalNoise_increase_len == TOTALNOISE_HIST_SIZE && st->ini_frame > 150 )
        ||
        ( st->input_bwidth > st->last_input_bwidth )
        ||
        ( st->last_core == AMR_WB_CORE )
    )
    {
        st->fd_cng_reset_flag = 1;
        st->hFdCngEnc->hFdCngCom->msFrCnt_init_counter = 0;
        st->hFdCngEnc->hFdCngCom->init_old = FLT_MAX;
    }
    else if ( st->fd_cng_reset_flag > 0 && st->fd_cng_reset_flag < 10 )
    {
        st->fd_cng_reset_flag++;
    }
    else
    {
        st->fd_cng_reset_flag = 0;
    }

    return;
}


/*-------------------------------------------------------------------*
* perform_noise_estimation_enc()
*
* Perform noise estimation
*-------------------------------------------------------------------*/

void perform_noise_estimation_enc(
    float *band_energies,         /* i: energy in critical bands without minimum noise floor E_MIN */
    float *enerBuffer,
    HANDLE_FD_CNG_ENC st          /* i/o: CNG structure containing all buffers and variables */
)
{
    short i, j;
    int   numCoreBands = st->hFdCngCom->numCoreBands;
    int   regularStopBand = st->hFdCngCom->regularStopBand;
    int   numSlots     = st->hFdCngCom->numSlots;
    float numSlots_inv = 1.f/(float)numSlots; /*enough if done only once*/
    float * periodog   = st->hFdCngCom->periodog;
    float * ptr_per    = periodog;
    int     npart      = st->hFdCngCom->npart;
    int     nFFTpart   = st->hFdCngCom->nFFTpart;
    float * psize      = st->hFdCngCom->psize;
    float * msPeriodog = st->msPeriodog;
    float * msNoiseEst = st->msNoiseEst;

    float * msLogPeriodog = st->msLogPeriodog;
    float * msLogNoiseEst = st->msLogNoiseEst;


    /* preemphasis compensation and grouping of per bin energies into msPeriodog */
    for(i=0; i<nFFTpart; i++)
    {
        msPeriodog[i] = 0.5f*(band_energies[i] + band_energies[i+NB_BANDS]);
        msPeriodog[i] *= preemphCompensation[i];
    }

    /* Adjust to the desired time resolution by averaging the periodograms over the time slots */
    for(j=numCoreBands ; j<regularStopBand ; j++)
    {
        (*ptr_per) = (enerBuffer[j] * numSlots_inv * st->hFdCngCom->scalingFactor);
        ptr_per++;
    }

    /* Adjust filterbank to the desired frequency resolution by averaging over spectral partitions for SID transmission */
    if (numCoreBands < regularStopBand)
    {
        bandcombinepow(periodog, regularStopBand-numCoreBands, st->hFdCngCom->CLDFBpart, st->hFdCngCom->nCLDFBpart, st->hFdCngCom->CLDFBpsize_inv, &msPeriodog[nFFTpart]);
    }

    /* Compress MS inputs */
    compress_range(msPeriodog, msLogPeriodog, npart);

    /* Call the minimum statistics routine for noise estimation */
    minimum_statistics( npart, nFFTpart, psize, msLogPeriodog, st->msNoiseFloor, msLogNoiseEst,
                        st->msAlpha, st->msPsd, st->msPsdFirstMoment, st->msPsdSecondMoment,
                        st->msMinBuf, st->msBminWin, st->msBminSubWin,
                        st->msCurrentMin, st->msCurrentMinOut, st->msCurrentMinSubWindow,
                        st->msLocalMinFlag, st->msNewMinFlag, st->msPeriodogBuf, &(st->msPeriodogBufPtr), st->hFdCngCom );

    /* Expand MS outputs */
    expand_range( msLogNoiseEst, msNoiseEst, npart );


    return;
}


/*-------------------------------------------------------------------*
* AdjustFirstSID()
*
* Adjust the noise estimator at the beginning of each CNG phase (encoder-side)
*-------------------------------------------------------------------*/

void AdjustFirstSID(
    int npart,
    float * msPeriodog,
    float * energy_ho,
    float * msNoiseEst,
    float * msNoiseEst_old,
    short * active_frame_counter,
    Encoder_State *stcod
)
{

    float lambda;
    int i;

    if ( stcod->cnt_SID == 1 && stcod->last_core_brate > SID_2k40 )
    {
        /* Detect the hangover period and the first SID frame at the beginning of each CNG phase */

        /* Average input energy over hangover period */
        mvr2r(msPeriodog, energy_ho, npart); /*First hangover frame*/
        /* Set first SID to current input level but add some smoothing */
        lambda = (float)pow( 0.96f, (float)(*active_frame_counter+1) );
        v_multc( msNoiseEst_old, lambda, msNoiseEst_old, npart );
        v_multc( energy_ho, 1-lambda, energy_ho, npart );

        v_add( msNoiseEst_old, energy_ho, energy_ho, npart );
        for (i=0; i<npart; i++)
        {
            if (msNoiseEst[i]>energy_ho[i])
            {
                msNoiseEst[i] = energy_ho[i];
            }
        }
        *active_frame_counter = 0;
    }

    if ( stcod->core_brate != SID_2k40 && stcod->core_brate != FRAME_NO_DATA )
    {
        (*active_frame_counter)++; /* Count the number of active frames in a row */
    }
    else
    {
        mvr2r( msNoiseEst, msNoiseEst_old, npart); /* Store the noise estimate obtained in the CNG phases */
    }

    return;
}

/*-------------------------------------------------------------------*
* FdCng_encodeSID()
*
* Generate a bitstream out of the partition levels
*-------------------------------------------------------------------*/

void FdCng_encodeSID(
    HANDLE_FD_CNG_ENC stenc,        /* i/o: CNG structure containing all buffers and variables */
    Encoder_State *corest,
    float  preemph_fac
)
{
    int N;
    HANDLE_FD_CNG_COM st = stenc->hFdCngCom;
    float* E = stenc->msNoiseEst;
    float gain;
    int i, index;
    float v[32], e;
    int indices[32];
    float w[32];
    (void)preemph_fac;


    /* Init */
    N = stenc->npartDec;

    /* Convert to LOG */
    e = 0.f;
    for ( i=0; i<N; i++ )
    {
        v[i] = 10.f*(float)log10( E[i] + 1e-4f );
        e += v[i];
    }

    /* Normalize MSVQ input */
    gain = 0.f;
    for ( i=N_GAIN_MIN; i<N_GAIN_MAX; i++ )
    {
        gain += v[i];
    }

    gain /= (float)(N_GAIN_MAX-N_GAIN_MIN);

    for ( i=0; i<N; i++ )
    {
        v[i] -= gain;
    }

    /* MSVQ encoder */
    set_f( w, 1.0f, N );

    msvq_enc( cdk_37bits, NULL, NULL, v, levels_37bits, maxC_37bits, stages_37bits, w, N, maxN_37bits, indices );

    /* MSVQ decoder */
    msvq_dec( cdk_37bits, NULL, NULL, stages_37bits, N, maxN_37bits, indices, v, NULL );

    /* Compute gain */
    gain = 0.f;
    for ( i=0; i<N; i++ )
    {
        gain += v[i];
    }

    gain = ( e - gain ) / (float)N;

    /* Apply bitrate-dependant scale */
    apply_scale( &gain, st->CngBandwidth, st->CngBitrate );

    /* Quantize gain */
    index = (short) floor( gain*1.5f + 60.f + 0.5f );

    if (index < 0)
    {
        index = 0;
    }

    if (index > 127)
    {
        index = 127;
    }

    gain = ((float)index-60.f)/1.5f;

    /* Apply gain and undo log */
    for ( i=0; i<N; i++ )
    {

        st->sidNoiseEst[i] = (float)pow( 10.f, (v[i]+gain)/10.f );
    }

    /* NB last band energy compensation */
    if (st->CngBandwidth == NB)
    {
        st->sidNoiseEst[N-1] *= NB_LAST_BAND_SCALE;
    }

    if ( st->CngBandwidth == SWB && st->CngBitrate <= ACELP_13k20 )
    {
        st->sidNoiseEst[N-1] *= SWB_13k2_LAST_BAND_SCALE;
    }

    /* Write bitstream */
    if ( corest->codec_mode == MODE2 )
    {
        for ( i=0; i<stages_37bits; i++ )
        {
            push_next_indice(corest, indices[i], bits_37bits[i]);
        }

        push_next_indice(corest, index, 7);
    }
    else
    {
        push_indice( corest, IND_SID_TYPE, 1, 1 );
        push_indice( corest, IND_ACELP_16KHZ, corest->bwidth, 2 );
        push_indice( corest, IND_ACELP_16KHZ, corest->L_frame == L_FRAME16k ? 1 : 0, 1 );

        for ( i=0; i<stages_37bits; i++ )
        {
            push_indice( corest, IND_LSF, indices[i], bits_37bits[i] );
        }

        push_indice( corest, IND_ENERGY, index, 7 );
    }

    /* Interpolate the bin/band-wise levels from the partition levels */
    scalebands(st->sidNoiseEst, stenc->partDec, stenc->npartDec, stenc->midbandDec, stenc->nFFTpartDec, stenc->stopBandDec-stenc->startBandDec, st->cngNoiseLevel, 1);

    lpc_from_spectrum(st->cngNoiseLevel, stenc->startBandDec, stenc->stopFFTbinDec, st->fftlen, st->fftSineTab, st->A_cng, preemph_fac );

    return;
}


void generate_comfort_noise_enc( Encoder_State *stcod )
{
    short i;
    float * ptr_r;
    float * ptr_i;
    HANDLE_FD_CNG_ENC stenc = stcod->hFdCngEnc;
    HANDLE_FD_CNG_COM st = stenc->hFdCngCom;
    float * cngNoiseLevel = st->cngNoiseLevel;
    float * ptr_level = cngNoiseLevel;
    short * seed = &(st->seed);
    float scale = 1.f;
    float * fftBuffer = st->fftBuffer;
    float * timeDomainOutput = st->timeDomainBuffer;
    float preemph_fac = stcod->preemph_fac;
    int tcx_transition = 0;
    float enr, att;
    short bwidth = stcod->bwidth;
    short CNG_mode = stcod->CNG_mode;

    /*
      Generate Gaussian random noise in real and imaginary parts of the FFT bins
      Amplitudes are adjusted to the estimated noise level cngNoiseLevel in each bin
    */
    if (stenc->startBandDec==0)
    {
        rand_gauss(&fftBuffer[0], seed);
        fftBuffer[0] *= (float)sqrt(scale **ptr_level); /* DC component in FFT */
        ptr_level++;
        ptr_r = fftBuffer + 2;
    }
    else
    {
        fftBuffer[0] = 0.f;
        set_f( fftBuffer+2, 0.0f, 2*(stenc->startBandDec-1) );
        ptr_r = fftBuffer + 2*stenc->startBandDec;
    }
    ptr_i = ptr_r + 1;
    for( ; ptr_level < cngNoiseLevel+stenc->stopFFTbinDec-stenc->startBandDec ; ptr_level++)
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
    set_f( fftBuffer+2*stenc->stopFFTbinDec, 0.0f, st->fftlen-2*stenc->stopFFTbinDec);
    /* Nyquist frequency is discarded */
    fftBuffer[1] = 0.f;

    /* If previous frame is active, reset the overlap-add buffer */
    if( stcod->last_core_brate > SID_2k40 )
    {
        set_f( st->olapBufferSynth, 0.0f, st->fftlen);
        if( (stcod->last_core>0 && stcod->codec_mode == MODE2)||stcod->codec_mode == MODE1 )
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

    stcod->lp_ener = 0.8f * stcod->lp_ener + 0.2f * pow( 2.0f, enr );


    /* Overlap-add when previous frame is active */
    if( stcod->last_core_brate > SID_2k40 && stcod->codec_mode == MODE2 )
    {
        float noise[2048], old_exc_ener = 0.f, gain = 0.f, tmp;
        int N = st->frameSize;
        short seed = st->seed;
        float *old_exc, old_Aq[M+1], *old_syn_pe, old_syn;

        if( stcod->last_core > 0 )
        {
            tcx_windowing_synthesis_current_frame(  timeDomainOutput,
                                                    stcod->tcx_cfg.tcx_mdct_window, /*Keep sine windows for limiting Time modulation*/
                                                    stcod->tcx_cfg.tcx_mdct_window_half,
                                                    stcod->tcx_cfg.tcx_mdct_window_minimum,
                                                    stcod->tcx_cfg.tcx_mdct_window_length,
                                                    stcod->tcx_cfg.tcx_mdct_window_half_length,
                                                    stcod->tcx_cfg.tcx_mdct_window_min_length, 0,
                                                    stcod->tcx_cfg.tcx_last_overlap_mode==ALDO_WINDOW?FULL_OVERLAP:stcod->tcx_cfg.tcx_last_overlap_mode,
                                                    NULL, NULL, NULL, NULL, NULL, N/2, stcod->tcx_cfg.tcx_offset<0?-stcod->tcx_cfg.tcx_offset:0,
                                                    1, 0, 0 );

            if(stcod->tcx_cfg.last_aldo)
            {
                for (i=0; i<st->frameSize; i++)
                {
                    timeDomainOutput[i] += stcod->old_out[i+NS2SA(stcod->sr_core, N_ZERO_MDCT_NS)];
                }
            }
            else
            {
                tcx_windowing_synthesis_past_frame( stcod->LPDmem.Txnq,
                                                    stcod->tcx_cfg.tcx_aldo_window_1_trunc,
                                                    stcod->tcx_cfg.tcx_mdct_window_half,
                                                    stcod->tcx_cfg.tcx_mdct_window_minimum,
                                                    stcod->tcx_cfg.tcx_mdct_window_length,
                                                    stcod->tcx_cfg.tcx_mdct_window_half_length,
                                                    stcod->tcx_cfg.tcx_mdct_window_min_length,
                                                    stcod->tcx_cfg.tcx_last_overlap_mode );

                for (i=0; i<stcod->tcx_cfg.tcx_mdct_window_length; i++)
                {
                    timeDomainOutput[i] += stcod->LPDmem.Txnq[i];
                }
            }
        }
        else
        {
            lsp2a_stab( stcod->lsp_old, old_Aq, M );
            old_exc = stcod->LPDmem.old_exc+L_EXC_MEM-(N/2);
            old_syn_pe = stcod->LPDmem.mem_syn2;
            old_syn = stcod->LPDmem.syn[M];
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
