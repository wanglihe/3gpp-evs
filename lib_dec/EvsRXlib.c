/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "options.h"
#include "prot.h"
#include "EvsRXlib.h"
#include "jbm_jb4sb.h"
#include "jbm_pcmdsp_apa.h"
#include "jbm_pcmdsp_fifo.h"
#include "cnst.h"


struct EVS_RX
{
    unsigned int             nSamplesFrame;
    Decoder_State           *st;
    JB4_HANDLE               hJBM;
    unsigned int             lastDecodedWasActive;
    PCMDSP_APA_HANDLE        hTimeScaler;
    PCMDSP_FIFO_HANDLE       hFifoAfterTimeScaler;
#ifdef SUPPORT_JBM_TRACEFILE
    FILE                    *jbmTraceFile;
#endif
};

/* function to check if a frame contains a SID */
static int isSidFrame( unsigned int size );


/* Opens the EVS Receiver instance. */
EVS_RX_ERROR EVS_RX_Open(EVS_RX_HANDLE* phEvsRX,
                         Decoder_State *st,
                         Word16 jbmSafetyMargin)
{
    EVS_RX_HANDLE hEvsRX;
    uint16_t wss, css;

    *phEvsRX = NULL;

    /* Create EVS Receiver handle */
    *phEvsRX = (EVS_RX_HANDLE) calloc(1, sizeof(struct EVS_RX) );
    if ( *phEvsRX == NULL )
    {
        return EVS_RX_MEMORY_ERROR;
    }
    hEvsRX = *phEvsRX;

    hEvsRX->st = st;
    /* do not use codec for time stretching (PLC) before initialization with first received frame */
    st->codec_mode = 0;

    /* open JBM */
    hEvsRX->hJBM = 0;
    if( JB4_Create(&(hEvsRX->hJBM)) != 0)
    {
        EVS_RX_Close(phEvsRX);
        return EVS_RX_INIT_ERROR;
    }

    /* init JBM */
    if(JB4_Init(hEvsRX->hJBM, jbmSafetyMargin) != 0)
    {
        EVS_RX_Close(phEvsRX);
        return EVS_RX_INIT_ERROR;
    }


    hEvsRX->lastDecodedWasActive = 0;
    hEvsRX->nSamplesFrame = st->output_Fs / 50;

    if(st->output_Fs == 8000)
    {
        wss = 1;
        css = 1;
    }
    else if(st->output_Fs == 16000)
    {
        wss = 2;
        css = 1;
    }
    else if(st->output_Fs == 32000)
    {
        wss = 4;
        css = 2;
    }
    else if(st->output_Fs == 48000)
    {
        wss = 6;
        css = 3;
    }
    else
    {
        assert(0 || "unknown sample rate!");
        wss = css = 1; /* just to avoid compiler warning */
    }

    /* initialize time scaler and FIFO after time scaler */
    if( apa_init( &hEvsRX->hTimeScaler ) != 0 ||
            apa_set_rate( hEvsRX->hTimeScaler, st->output_Fs, 1 ) != 0 ||
            apa_set_complexity_options( hEvsRX->hTimeScaler, wss, css) != 0 ||
            apa_set_quality( hEvsRX->hTimeScaler, 1, 4, 4 ) != 0 ||
            pcmdsp_fifo_create( &hEvsRX->hFifoAfterTimeScaler ) != 0 ||
            pcmdsp_fifo_init( hEvsRX->hFifoAfterTimeScaler, st->output_Fs * 4 / 50 /* 4 frames */, 1, 2 /* Word16 */ ) != 0 )
    {
        EVS_RX_Close(phEvsRX);
        return EVS_RX_TIMESCALER_ERROR;
    }

    return EVS_RX_NO_ERROR;
}

#ifdef SUPPORT_JBM_TRACEFILE
/* Sets the name of the JBM trace file which will be created. */
EVS_RX_ERROR
EVS_RX_SetJbmTraceFileName(EVS_RX_HANDLE hEvsRX,
                           const char *jbmTraceFileName)
{
    /* JBM trace file writing is only done for EVS testing and is not instrumented. */
    if( hEvsRX->jbmTraceFile )
        fclose( hEvsRX->jbmTraceFile );
    if( jbmTraceFileName != NULL )
    {
        hEvsRX->jbmTraceFile = fopen( jbmTraceFileName, "w" );
        if( !hEvsRX->jbmTraceFile )
        {
            return EVS_RX_WRONG_PARAMS;
        }
        fprintf( hEvsRX->jbmTraceFile, "#rtpSeqNo;rtpTs;rcvTime;playTime;active\n" );
    }
    return EVS_RX_NO_ERROR;
}
#endif

/* Feeds one frame into the receiver. */
EVS_RX_ERROR
EVS_RX_FeedFrame(EVS_RX_HANDLE hEvsRX,
                 unsigned char *au,
                 unsigned int auSize,
                 unsigned short rtpSequenceNumber,
                 unsigned long rtpTimeStamp,
                 unsigned int rcvTime_ms)
{
    JB4_DATAUNIT_HANDLE dataUnit;
    int16_t partialCopyFrameType, partialCopyOffset;
    int result;

    assert( auSize != 0 );
    assert( (auSize + 7) / 8 <= MAX_AU_SIZE );

    /* check if frame contains a partial copy and get its offset */
    evs_dec_previewFrame(au, auSize, &partialCopyFrameType, &partialCopyOffset);

    /* create data unit for primary copy in the frame */
    dataUnit = JB4_AllocDataUnit(hEvsRX->hJBM);
    memcpy(dataUnit->data, au, (auSize + 7) / 8);
    dataUnit->dataSize = auSize;
    dataUnit->duration = 20;
    dataUnit->sequenceNumber = rtpSequenceNumber;
    dataUnit->silenceIndicator = isSidFrame( dataUnit->dataSize );
    dataUnit->timeScale = 1000;
    dataUnit->rcvTime = rcvTime_ms;
    dataUnit->timeStamp = rtpTimeStamp;
    dataUnit->partial_frame = 0;
    dataUnit->partialCopyOffset = partialCopyOffset;

    /* add the frame to the JBM */
    result = JB4_PushDataUnit(hEvsRX->hJBM, dataUnit, rcvTime_ms);
    if(result != 0)
    {
        return EVS_RX_JBM_ERROR;
    }

    if(partialCopyFrameType != RF_NO_DATA && partialCopyOffset != 0)
    {
        /* create data unit for partial copy in the frame */
        dataUnit = JB4_AllocDataUnit(hEvsRX->hJBM);
        memcpy(dataUnit->data, au, (auSize + 7) / 8);
        dataUnit->dataSize = auSize;
        dataUnit->duration = 20;
        dataUnit->sequenceNumber = rtpSequenceNumber;
        dataUnit->silenceIndicator = 0; /* there are no partial copies for SID frames */
        dataUnit->timeScale = 1000;
        dataUnit->rcvTime = rcvTime_ms;
        dataUnit->timeStamp = rtpTimeStamp - partialCopyOffset * dataUnit->duration;
        dataUnit->partial_frame = 1;
        dataUnit->partialCopyOffset = partialCopyOffset;

        /* add the frame to the JBM */
        result = JB4_PushDataUnit(hEvsRX->hJBM, dataUnit, rcvTime_ms);
        if(result != 0)
        {
            return EVS_RX_JBM_ERROR;
        }
    }
    return EVS_RX_NO_ERROR;
}


/* Retrieves one frame of output PCM data. */
EVS_RX_ERROR
EVS_RX_GetSamples(EVS_RX_HANDLE hEvsRX,
                  unsigned int*  nOutSamples,
                  Word16        *pcmBuf,
                  unsigned int   pcmBufSize,
                  unsigned int  systemTimestamp_ms
                 )
{
    Decoder_State *st;
    unsigned int soundCardFrameSize, extBufferedSamples;
    uint32_t extBufferedTime_ms, scale, maxScaling;
    uint16_t nTimeScalerOutSamples;
    int timeScalingDone, result;
    JB4_DATAUNIT_HANDLE dataUnit;
    float output[3 * L_FRAME48k];       /* 'float' buffer for output synthesis */

    assert(hEvsRX->nSamplesFrame <= pcmBufSize);
    assert(hEvsRX->nSamplesFrame <= APA_BUF);

    st = hEvsRX->st;
    soundCardFrameSize = hEvsRX->nSamplesFrame;
    timeScalingDone = 0;


    /* make sure that the FIFO after decoder/scaler contains at least one sound card frame (i.e. 20ms) */
    while( pcmdsp_fifo_nReadableSamples( hEvsRX->hFifoAfterTimeScaler ) < soundCardFrameSize )
    {
        extBufferedSamples = pcmdsp_fifo_nReadableSamples( hEvsRX->hFifoAfterTimeScaler );
        extBufferedTime_ms = extBufferedSamples * 1000 / st->output_Fs;
        dataUnit = NULL;
        /* pop one access unit from the jitter buffer */
        result = JB4_PopDataUnit(hEvsRX->hJBM, systemTimestamp_ms, extBufferedTime_ms, &dataUnit, &scale, &maxScaling);
        if(result != 0)
        {
            return EVS_RX_JBM_ERROR;
        }
        maxScaling = maxScaling * st->output_Fs / 1000;
        /* avoid time scaling multiple times in one sound card slot */
        if( scale != 100U )
        {
            if( timeScalingDone )
                scale = 100;
            else
                timeScalingDone = 1;
        }

        /* copy bitstream into decoder state */
        if(dataUnit)
        {
            if( st->codec_mode != 0 )
            {
                read_indices_from_djb( st, dataUnit->data, dataUnit->dataSize, (dataUnit->partial_frame==TRUE)? 1:0, dataUnit->nextCoderType );
                if(dataUnit->partial_frame != 0)
                {
                    st->codec_mode = MODE2;
                    st->use_partial_copy = 1;
                }
            }
            else /* initialize decoder with first received frame */
            {
                /* initialize, since this is needed within read_indices_from_djb, to correctly set st->last_codec_mode */
                st->ini_frame = 0;
                st->prev_use_partial_copy = 0;
                /* initialize st->last_codec_mode, since this is needed for init_decoder() */
                read_indices_from_djb( st, dataUnit->data, dataUnit->dataSize, 0, 0 );

                assert(st->codec_mode != 0);
                init_decoder( st );
                /* parse frame again because init_decoder() overwrites st->total_brate */
                read_indices_from_djb( st, dataUnit->data, dataUnit->dataSize, 0, 0 );

            }
        }
        else if( st->codec_mode != 0 )
        {
            read_indices_from_djb( st, NULL, 0, 0, 0 );
        }

        /* run the main decoding routine */
        if( st->codec_mode == MODE1 )
        {
            if( st->Opt_AMR_WB )
            {
                amr_wb_dec( st, output );
            }
            else
            {
                evs_dec( st, output, FRAMEMODE_NORMAL );
            }
        }
        else if( st->codec_mode == MODE2 )
        {
            if(st->bfi == 0)
            {
                evs_dec(st, output, FRAMEMODE_NORMAL);   /* FRAMEMODE_NORMAL */
            }
            else if ( st->bfi == 2 )
            {
                evs_dec(st, output, FRAMEMODE_FUTURE);   /* FRAMEMODE_FUTURE */
            }
            else /* conceal */
            {
                evs_dec(st, output, FRAMEMODE_MISSING);
            }
        }
        /* convert 'float' output data to 'short' */
        if( st->codec_mode == MODE1 || st->codec_mode == MODE2 )
        {
            syn_output( output, hEvsRX->nSamplesFrame, pcmBuf );

            /* increase the counter of initialization frames */
            if( st->ini_frame < MAX_FRAME_COUNTER )
            {
                st->ini_frame++;
            }
        }
        else /* codec mode to use not known yet */
        {
            set_s( pcmBuf, 0, hEvsRX->nSamplesFrame );
        }

        if(dataUnit)
        {
            if(dataUnit->partial_frame != 0)
            {
                hEvsRX->lastDecodedWasActive = 1;
            }
            else
            {
                hEvsRX->lastDecodedWasActive = !dataUnit->silenceIndicator;
            }
            /* data unit memory is no longer used */
            JB4_FreeDataUnit(hEvsRX->hJBM, dataUnit);
        }

        /* limit scale to range supported by time scaler */
        if( scale < APA_MIN_SCALE )
            scale = APA_MIN_SCALE;
        else if( scale > APA_MAX_SCALE )
            scale = APA_MAX_SCALE;
        /* apply time scaling on decoded/concealed samples */
        if( apa_set_scale( hEvsRX->hTimeScaler, scale ) != 0 )
        {
            return EVS_RX_TIMESCALER_ERROR;
        }
        result = apa_exec( hEvsRX->hTimeScaler, pcmBuf, hEvsRX->nSamplesFrame,
                           maxScaling, pcmBuf, &nTimeScalerOutSamples );
        if( result != 0 )
        {
            return EVS_RX_TIMESCALER_ERROR;
        }
        assert(nTimeScalerOutSamples <= pcmBufSize);
        assert(nTimeScalerOutSamples <= APA_BUF);
        (void)pcmBufSize;
        /* append scaled samples to FIFO */
        if( pcmdsp_fifo_write( hEvsRX->hFifoAfterTimeScaler,
                               (uint8_t*)pcmBuf, nTimeScalerOutSamples ) != 0 )
        {
            return EVS_RX_TIMESCALER_ERROR;
        }
#ifdef SUPPORT_JBM_TRACEFILE
        /* write JBM trace file entry */
        /* JBM trace file writing is only done for EVS testing and is not instrumented. */
        if( hEvsRX->jbmTraceFile )
        {
            /* the first sample of the decoded/concealed frame will be played after the samples in the ring buffer */
            double playTime = systemTimestamp_ms + extBufferedSamples * 1000.0 / st->output_Fs;
            /* rtpSeqNo;rtpTs;rcvTime;playTime;active\n */
            if( dataUnit )
            {
                if(dataUnit->partial_frame == 1)
                {
                    fprintf( hEvsRX->jbmTraceFile, "%d;%d;%d;%f;%d;%d\n",
                             -1, -1, -1, playTime, hEvsRX->lastDecodedWasActive, dataUnit->partialCopyOffset );
                }
                else
                {
                    fprintf( hEvsRX->jbmTraceFile, "%u;%u;%u;%f;%d\n",
                             dataUnit->sequenceNumber, dataUnit->timeStamp, dataUnit->rcvTime,
                             playTime, hEvsRX->lastDecodedWasActive );
                }

            }
            else
            {
                fprintf( hEvsRX->jbmTraceFile, "%d;%d;%d;%f;%d\n",
                         -1, -1, -1,
                         playTime, hEvsRX->lastDecodedWasActive );
            }
        }
#endif
    }

    /* fetch one frame for the sound card from FIFO */
    *nOutSamples = soundCardFrameSize;
    if( pcmdsp_fifo_read( hEvsRX->hFifoAfterTimeScaler, *nOutSamples, (uint8_t*)pcmBuf ) != 0 )
    {
        return EVS_RX_TIMESCALER_ERROR;
    }
    return EVS_RX_NO_ERROR;
}

EVS_RX_ERROR
EVS_RX_Get_FEC_offset( EVS_RX_HANDLE hEvsRX, short *offset
                       , short *FEC_hi
                     )
{

    *offset = JB4_getFECoffset(hEvsRX->hJBM);
    *FEC_hi = JB4_FECoffset(hEvsRX->hJBM);
    return EVS_RX_NO_ERROR;

}




/* Returns 1 if the jitter buffer is empty, otherwise 0. */
unsigned int
EVS_RX_IsEmpty(EVS_RX_HANDLE hEvsRX )
{
    unsigned int isEmpty;

    isEmpty = 0;
    if(JB4_bufferedDataUnits(hEvsRX->hJBM) == 0U)
        isEmpty = 1;

    return isEmpty;
}

/* Closes the receiver instance. */
EVS_RX_ERROR
EVS_RX_Close(EVS_RX_HANDLE* phRX )
{
    /* Free all memory */
    if( phRX == NULL || *phRX == NULL )
    {
        return EVS_RX_NO_ERROR;
    }

    destroy_decoder( (*phRX)->st );

    if( (*phRX)->hJBM )
        JB4_Destroy( &(*phRX)->hJBM );


    if( (*phRX)->hTimeScaler )
        apa_exit( &(*phRX)->hTimeScaler );

    if( (*phRX)->hFifoAfterTimeScaler )
        pcmdsp_fifo_destroy( &(*phRX)->hFifoAfterTimeScaler );

#ifdef SUPPORT_JBM_TRACEFILE
    if( (*phRX)->jbmTraceFile )
        fclose( (*phRX)->jbmTraceFile );
#endif

    free( *phRX );
    *phRX = NULL;
    phRX = NULL;

    return EVS_RX_NO_ERROR;
}

/* function to check if a frame contains a SID */
static int isSidFrame( unsigned int size )
{
    int ret = 0;
    if(size == SID_1k75 / 50)
    {
        ret = 1; /* AMR-WB SID */
    }
    else if(size == SID_2k40 / 50)
    {
        ret = 1; /* EVS SID */
    }
    return ret;
}
