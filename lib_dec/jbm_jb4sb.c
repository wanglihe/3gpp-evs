/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/


/*! \file jbm_jb4sb.c EVS Jitter Buffer Management Interface */

/* system headers */
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "options.h"
/* instrumentation */
/* local headers */
#include "jbm_jb4_circularbuffer.h"
#include "jbm_jb4_inputbuffer.h"
#include "jbm_jb4_jmf.h"
#include "jbm_jb4sb.h"
#include "prot.h"


#define JB4_MIN(a,b)  ((a)>(b) ? (b) : (a))
#define JB4_MAX(a,b)  ((a)>(b) ? (a) : (b))

#define MAXOFFSET 10

/*! Calculates the difference between two RTP timestamps - the diff is positive, if B 'later', negative otherwise */
static int32_t JB4_rtpTimeStampDiff( uint32_t tsA, uint32_t tsB );
/* function to calculate different options for the target playout delay */
static void JB4_targetPlayoutDelay( const JB4_HANDLE h, uint32_t *targetMin,
                                    uint32_t *targetMax, uint32_t *targetDtx, uint32_t *targetStartUp );
/*! function to do playout adaptation before playing the next data unit */
/*! In case of time shrinking, data units will be dropped before the next data unit to play is returned and
 *  in case of time stretching a empty data unit is returned and the frame should be concealed.
 *  @param[in] now current system time
 *  @param[out] dataUnit the next data unit to play
 *  @param[out] scale the scale in percent used as target for time scaling of the returned data unit
 *  @param[out] maxScaling the maximum allowed external time scaling */
static int JB4_adaptPlayout( JB4_HANDLE h, uint32_t sysTime, uint32_t extBufferedTime,
                             JB4_DATAUNIT_HANDLE *pDataUnit, uint32_t *scale, uint32_t *maxScaling );
/*! function to do playout adaptation before playing the first data unit */
/*! @param[in] now current system time
 *  @param[out] prebuffer true, if the data unit should be prebuffered */
static void JB4_adaptFirstPlayout( JB4_HANDLE h, uint32_t sysTime, bool_t *prebuffer );
/*! function for playout adaptation while active (no DTX) */
static void JB4_adaptActivePlayout( JB4_HANDLE h, uint32_t sysTime,
                                    uint32_t extBufferedTime, uint32_t *scale, uint32_t *maxScaling );
/*! function for playout adaptation while DTX */
static void JB4_adaptDtxPlayout( JB4_HANDLE h, uint32_t sysTime, bool_t *stretchTime );
/*! function to look into the buffer and check if it makes sense to drop a data unit */
/*! @param[out] dropEarly true, if a data unit could be dropped early
 *  @param[out] buffered the buffered time span in timeScale units
 *  @return true, if a data unit could be dropped */
static int JB4_inspectBufferForDropping( const JB4_HANDLE h, bool_t *dropEarly, uint32_t *buffered );
/* function to look into the buffer and check if it makes sense to drop a data unit during DTX */
static int JB4_checkDtxDropping( const JB4_HANDLE h );
/*! function to estimate the short term jitter  */
static void JB4_estimateShortTermJitter( JB4_HANDLE h, uint32_t rcvTime, uint32_t rtpTimeStamp );
/*! function to pop a data unit from the buffer */
static void JB4_popFromBuffer( JB4_HANDLE h, uint32_t sysTime, JB4_DATAUNIT_HANDLE *pDataUnit );
/*! function to drop a data unit from the buffer - updates nShrinked */
static void JB4_dropFromBuffer( JB4_HANDLE h, uint32_t sysTime );
/*! function to calculate the playout delay based on the current jitter */
/*! @param[in] playTime the system time when the data unit will be played
 *  @param[in] timeStamp the time stamp of the data unit to played
 *  @param[out] delay the calculated playout delay */
static int JB4_playoutDelay( const JB4_HANDLE h, uint32_t playTime,
                             uint32_t rtpTimeStamp, uint32_t *delay );
/*! function to update lastPlayoutDelay and lastTargetTime after popFromBuffer() */
static void JB4_updateLastTimingMembers( JB4_HANDLE h, uint32_t playTime, uint32_t rtpTimeStamp );
/*! function to compare the RTP time stamps of two data units: newElement==arrayElement ? 0 : (newElement>arrayElement ? +1 : -1) */
static int JB4_inputBufferCompareFunction( const JB4_INPUTBUFFER_ELEMENT newElement,
        const JB4_INPUTBUFFER_ELEMENT arrayElement, bool_t *replaceWithNewElementIfEqual );


/*! Jitter Buffer Management Interface */
struct JB4
{
    /*! @name statistics for user */
    /*@{ */
    /*! the number of late lost data units */
    uint32_t                     nLateLost;
    /*! the number of data units that were available (not NULL) at playout time */
    uint32_t                     nAvailablePopped;
    /*! the number of data units that were not available (NULL) at playout time */
    uint32_t                     nUnavailablePopped;
    /*! the number of unavailable pops since the last available one - used as temp value for nLost and nStretched */
    uint32_t                     nLostOrStretched;
    /*! the number of data units that were lost at playout time */
    uint32_t                     nLost;
    /*! the number of empty data units inserted for playout adaptation */
    uint32_t                     nStretched;
    /*! the number of data units dropped for playout adaptation */
    /*! This function counts all time shrinking events, no matter if a dropped data unit was actually available. */
    uint32_t                     nShrinked;
    /*! the number of data units that were returned to create comfort noice (including NULL) */
    uint32_t                     nComfortNoice;
    /*! the number of jitter induced concealment operations (as defined in 3GPP TS 26.114) */
    uint32_t                     jitterInducedConcealments;
    /*! the target playout delay of the last returned data unit */
    uint32_t                     targetPlayoutDelay;
    /*! the target playout time of the last returned data unit */
    uint32_t                     lastTargetTime;
    /*@} */
    /*! @name internal configuration values - do not change!!! */
    /*@{ */
    /*! internal time scale for all calculations */
    int                          timeScale;
    /*! internal frame duration in timeScale units */
    uint32_t                     frameDuration;
    /*@} */
    /*! @name jitter buffer configuration values */
    /*@{ */
    /*! the allowed delay reserve in addition to network jitter to reduce late-loss [milliseconds] */
    Word32                       safetyMargin;
    /*@} */
    /*! @name data for short term jitter estimation */
    /*@{ */
    /*! short term jitter measure FIFO */
    JB4_JMF_HANDLE               stJmf;
    /*! FIFO of short term jitter values */
    JB4_CIRCULARBUFFER_HANDLE    stJitterFifo;
    /*! FIFO of RTP time stamps for the values stored in stJitterFifo */
    JB4_CIRCULARBUFFER_HANDLE    stTimeStampFifo;
    /*! short term jitter */
    uint32_t                     stJitter;
    /*@} */
    /*! @name jitter buffer data */
    /*@{ */
    /*! true, if a data unit was already popped from the buffer */
    bool_t                       firstDataUnitPopped;
    /*! system time of the previous JB4_PopDataUnit() call */
    uint32_t                     prevPopSysTime;
    /*! RTP timestamp of the last played/dropped data unit that was actually available */
    uint32_t                     lastReturnedTs;
    /*! true, if the last popped data unit contained no active signal, i.e. silence -> hint for DTX */
    bool_t                       lastPoppedWasSilence;
    /*! the playout time minus the minimum offset of the last played data unit in microseconds */
    int32_t                      lastPlayoutOffset;
    /*! RTP time stamp of the next data unit that is expected to be fetched from the buffer */
    uint32_t                     nextExpectedTs;
    Word16                       rfOffset2Active;
    Word16                       rfOffset3Active;
    Word16                       rfOffset5Active;
    Word16                       rfOffset7Active;
    Word32                       rfDelay;
    /*! long term jitter measure FIFO */
    JB4_JMF_HANDLE               ltJmf;

    uint32_t  FecOffWinLen;
    uint32_t  FecOffWin[10];
    uint32_t  optimum_offset;

    float  netLossRate;
    Word32  nPartialCopiesUsed;
    Word32 last_nLost;
    Word32 last_ntot;

    uint32_t  totWin;
    bool_t                     pre_partial_frame;
    /*@} */

    /*! @name members to store the data units */
    /*@{ */
    /*! the data unit buffer */
    JB4_INPUTBUFFER_HANDLE       inputBuffer;
    struct JB4_DATAUNIT          memorySlots[MAX_JBM_SLOTS];
    JB4_DATAUNIT_HANDLE          freeMemorySlots[MAX_JBM_SLOTS];
    unsigned int                 nFreeMemorySlots;
    /*@} */
}; /* JB4 */


int JB4_Create( JB4_HANDLE *ph )
{
    JB4_HANDLE h = calloc( 1, sizeof( struct JB4 ) );
    short iter;

    /* statistics for user */
    h->nLateLost                 = 0;
    h->nAvailablePopped          = 0;
    h->nUnavailablePopped        = 0;
    h->nLostOrStretched          = 0;
    h->nLost                     = 0;
    h->nStretched                = 0;
    h->nShrinked                 = 0;
    h->nComfortNoice             = 0;
    h->jitterInducedConcealments = 0;
    h->targetPlayoutDelay        = 0;
    h->lastTargetTime            = 0;
    /* internal configuration values - do not change!!! */
    h->timeScale                 = 0;
    h->frameDuration             = 0;
    /* jitter buffer configuration values: done in JB4_Init() */
    /* short term jitter evaluation */
    JB4_JMF_Create( &h->stJmf );
    JB4_CIRCULARBUFFER_Create( &h->stJitterFifo );
    JB4_CIRCULARBUFFER_Create( &h->stTimeStampFifo );
    h->stJitter                  = 0;
    /* jitter buffer data */
    h->firstDataUnitPopped       = false;
    h->prevPopSysTime            = 0;
    h->lastReturnedTs            = 0;
    h->lastPoppedWasSilence      = false;
    h->lastPlayoutOffset         = 0;
    h->nextExpectedTs            = 0;
    h->rfOffset2Active           = 0;
    h->rfOffset3Active           = 0;
    h->rfOffset5Active           = 0;
    h->rfOffset7Active           = 0;
    h->rfDelay                   = 0;
    JB4_JMF_Create( &h->ltJmf );
    h->pre_partial_frame         = 0;

    h->FecOffWinLen               = 0;
    for (iter = 0; iter < 10; iter++ )
    {
        h->FecOffWin[iter] = 0;
    }
    h->optimum_offset             = 3;
    h->totWin                     = 0;
    h->netLossRate                = 0.0f;
    move32();
    h->nPartialCopiesUsed         = 0;
    move32();
    h->last_nLost                 = 0;
    move32();
    h->last_ntot                  = 0;
    move32();

    /* members to store the data units */
    JB4_INPUTBUFFER_Create( &h->inputBuffer );
    /* allocate memory for data units */
    for(iter = 0; iter < MAX_JBM_SLOTS; ++iter)
    {
        h->memorySlots[iter].data = malloc(MAX_AU_SIZE);
        h->freeMemorySlots[iter] = &h->memorySlots[iter];
    }
    h->nFreeMemorySlots = MAX_JBM_SLOTS;
    *ph = h;
    return 0;
}

void JB4_Destroy( JB4_HANDLE *ph )
{
    JB4_HANDLE h;
    unsigned int i;

    if( !ph )
    {
        return;
    }
    h = *ph;
    if( !h )
    {
        return;
    }

    JB4_JMF_Destroy( &h->stJmf );
    JB4_CIRCULARBUFFER_Destroy( &h->stJitterFifo );
    JB4_CIRCULARBUFFER_Destroy( &h->stTimeStampFifo );
    JB4_JMF_Destroy( &h->ltJmf );
    JB4_INPUTBUFFER_Destroy( &h->inputBuffer );

    for(i = 0; i < MAX_JBM_SLOTS; ++i)
    {
        free(h->memorySlots[i].data);
    }

    free( h );
    *ph = NULL;

}

int JB4_Init( JB4_HANDLE h, Word16 safetyMargin )
{
    unsigned int ltJmfSize, stFifoSize, stJmfSize, stJmfAllowedLateLoss;
    unsigned int inputBufferCapacity;

    /* internal timescale is 1000, frame duration is 20ms */
    h->timeScale            = 1000; /* ms */
    h->frameDuration        = 20;   /* ms */

    /* jitter buffer configuration values */
    h->safetyMargin         = safetyMargin;

    /* long term jitter measure FIFO: 500 frames and 10s */
    ltJmfSize = 10000;
    JB4_JMF_Init( h->ltJmf, h->timeScale, ltJmfSize / 20, ltJmfSize, 1000 );
    /* short term jitter evaluation */
    stFifoSize           = 200;
    stJmfSize            = 50;
    stJmfAllowedLateLoss = 60; /* 6%, e.g. ignore three packets out of 50 */
    JB4_CIRCULARBUFFER_Init( h->stJitterFifo, stFifoSize );
    JB4_CIRCULARBUFFER_Init( h->stTimeStampFifo, stFifoSize );
    JB4_JMF_Init( h->stJmf, h->timeScale,
                  stJmfSize, h->timeScale /* 1s */, 1000 - stJmfAllowedLateLoss );

    inputBufferCapacity = MAX_JBM_SLOTS - 2;
    JB4_INPUTBUFFER_Init( h->inputBuffer, inputBufferCapacity, JB4_inputBufferCompareFunction );
    return 0;
}

/* Returns a memory slot to store a new data unit */
JB4_DATAUNIT_HANDLE JB4_AllocDataUnit( JB4_HANDLE h )
{
    JB4_DATAUNIT_HANDLE dataUnit;
    while(h->nFreeMemorySlots == 0)
    {
        assert(JB4_INPUTBUFFER_IsEmpty(h->inputBuffer) == 0);
        JB4_dropFromBuffer(h, 0);
    }

    --h->nFreeMemorySlots;
    dataUnit = h->freeMemorySlots[h->nFreeMemorySlots];
    h->freeMemorySlots[h->nFreeMemorySlots] = NULL;
    assert(dataUnit != NULL);
    return dataUnit;
}

/* Notifies the JBM that a data unit is no longer used and the memory can be reused */
void JB4_FreeDataUnit( JB4_HANDLE h, JB4_DATAUNIT_HANDLE dataUnit )
{
    assert(dataUnit != NULL);
    assert(h->nFreeMemorySlots < MAX_JBM_SLOTS);
    h->freeMemorySlots[h->nFreeMemorySlots] = dataUnit;
    h->nFreeMemorySlots++;
}

int JB4_PushDataUnit( JB4_HANDLE h, JB4_DATAUNIT_HANDLE dataUnit, uint32_t rcvTime )
{
    JB4_DATAUNIT_HANDLE droppedDataUnit = NULL;

    assert( dataUnit->duration == h->frameDuration );
    assert( dataUnit->timeScale == (unsigned int)h->timeScale );

    /* ignore frames from too far in future (3 seconds) */
    if( h->firstDataUnitPopped && JB4_rtpTimeStampDiff( h->lastReturnedTs, dataUnit->timeStamp ) >=
            (int32_t) (50 * 3 * dataUnit->duration) )
    {
        JB4_FreeDataUnit(h, dataUnit);
        return 0;
    }

    /* reserve space for one element to add: drop oldest if buffer is full */
    while(JB4_INPUTBUFFER_IsFull(h->inputBuffer))
    {
        JB4_dropFromBuffer(h, rcvTime);
    }
    assert(JB4_INPUTBUFFER_IsFull(h->inputBuffer) == 0);

    /* do statistics on partial copy offset using active primary copies to
     * avoid unexpected resets because RF_NO_DATA partial copies are dropped before JBM */
    if(dataUnit->silenceIndicator == 0 && dataUnit->partial_frame == 0)
    {
        if(dataUnit->partialCopyOffset == 0)
        {
            if(h->rfOffset2Active > 0)
                --h->rfOffset2Active;
            if(h->rfOffset3Active > 0)
                --h->rfOffset3Active;
            if(h->rfOffset5Active > 0)
                --h->rfOffset5Active;
            if(h->rfOffset7Active > 0)
                --h->rfOffset7Active;
        }
        else if(dataUnit->partialCopyOffset == 2)
        {
            h->rfOffset2Active = 100;
            h->rfOffset3Active = 0;
            h->rfOffset5Active = 0;
            h->rfOffset7Active = 0;
        }
        else if(dataUnit->partialCopyOffset == 3)
        {
            h->rfOffset2Active = 0;
            h->rfOffset3Active = 100;
            h->rfOffset5Active = 0;
            h->rfOffset7Active = 0;
        }
        else if(dataUnit->partialCopyOffset == 5)
        {
            h->rfOffset2Active = 0;
            h->rfOffset3Active = 0;
            h->rfOffset5Active = 100;
            h->rfOffset7Active = 0;
        }
        else if(dataUnit->partialCopyOffset == 7)
        {
            h->rfOffset2Active = 0;
            h->rfOffset3Active = 0;
            h->rfOffset5Active = 0;
            h->rfOffset7Active = 100;
        }
    }

    if(dataUnit->partial_frame != 0)
    {
        /* check for "real" late loss: a frame with higher/same timestamp was already returned to be fed into decoder */
        if( h->firstDataUnitPopped && JB4_rtpTimeStampDiff( h->lastReturnedTs, dataUnit->timeStamp ) <= 0 )
        {
            JB4_FreeDataUnit(h, dataUnit);
            return 0;
        }

        /* drop partial copy if the missing frame was already concealed */
        if( h->firstDataUnitPopped )
        {
            if( dataUnit->partialCopyOffset <= 3 && JB4_rtpTimeStampDiff( h->nextExpectedTs, dataUnit->timeStamp ) < 0)
            {
                JB4_FreeDataUnit(h, dataUnit);
                return 0;
            }
            else if( dataUnit->partialCopyOffset == 5 && JB4_rtpTimeStampDiff( h->nextExpectedTs, dataUnit->timeStamp ) < -40)
            {
                JB4_FreeDataUnit(h, dataUnit);
                return 0;
            }
            else if( dataUnit->partialCopyOffset == 7 && JB4_rtpTimeStampDiff( h->nextExpectedTs, dataUnit->timeStamp ) < -80)
            {
                JB4_FreeDataUnit(h, dataUnit);
                return 0;
            }
        }

        /* try to store partial copy - will be dropped if primary copy already available */
        if(JB4_INPUTBUFFER_Enque( h->inputBuffer, dataUnit, (void**)&droppedDataUnit ) == 0)
        {
            /* partial copy is useful, consider it in long-term jitter estimation */
            if(dataUnit->partialCopyOffset <= 3)
            {
                JB4_JMF_PushPacket( h->ltJmf, rcvTime, dataUnit->timeStamp );
            }
        }
        else
        {
            JB4_FreeDataUnit(h, dataUnit);
        }
        if(droppedDataUnit != NULL)
        {
            JB4_FreeDataUnit(h, droppedDataUnit);
        }
    }
    else
    {
        /* calculate jitter */
        JB4_JMF_PushPacket( h->ltJmf, rcvTime, dataUnit->timeStamp );
        JB4_estimateShortTermJitter( h, rcvTime, dataUnit->timeStamp );
        /* check for "real" late loss: a frame with higher/same timestamp was already returned to be fed into decoder */
        if( h->firstDataUnitPopped && JB4_rtpTimeStampDiff( h->lastReturnedTs, dataUnit->timeStamp ) <= 0 )
        {
            if( !dataUnit->silenceIndicator )
            {
                ++h->nLateLost;
                /* deletion of a speech frame because it arrived at the JBM too late */
                ++h->jitterInducedConcealments;
            }
            JB4_FreeDataUnit(h, dataUnit);
            return 0;
        }
        /* store data unit */
        if(JB4_INPUTBUFFER_Enque( h->inputBuffer, dataUnit, (void**)&droppedDataUnit) != 0)
        {
            JB4_FreeDataUnit(h, dataUnit);
        }
        if(droppedDataUnit != NULL)
        {
            JB4_FreeDataUnit(h, droppedDataUnit);
        }
    }
    return 0;
}



int JB4_getFECoffset(JB4_HANDLE h)
{
    return (int)h->optimum_offset;
}
short JB4_FECoffset(JB4_HANDLE h)
{
    if ( h->netLossRate <  0.05 )
    {
        return (short)0;
    }
    else
    {
        return (short)1;
    }
}



int JB4_PopDataUnit( JB4_HANDLE h, uint32_t sysTime, uint32_t extBufferedTime,
                     JB4_DATAUNIT_HANDLE *pDataUnit, uint32_t *scale, uint32_t *maxScaling )
{
    int ret;

    assert( sysTime >= h->prevPopSysTime );
    if( sysTime > h->prevPopSysTime + 20 )
    {
        h->lastPlayoutOffset += 20;
    }
    h->prevPopSysTime = sysTime;


    ret = JB4_adaptPlayout( h, sysTime, extBufferedTime, pDataUnit, scale, maxScaling );


    return ret;
}

/* Calculates the difference between two RTP timestamps - the diff is positive, if B 'later', negative otherwise */
static int32_t JB4_rtpTimeStampDiff( uint32_t tsA, uint32_t tsB )
{
    int32_t ret;
    /* do not dare to inline this function, casting to int32_t is important here! */
    ret = (int32_t)(tsB - tsA);
    return ret;
}

/* function to get the number of data units contained in the buffer */
unsigned int JB4_bufferedDataUnits( const JB4_HANDLE h )
{
    return JB4_INPUTBUFFER_Size( h->inputBuffer );
}


/*****************************************************************************
 **************************** private functions ******************************
 *****************************************************************************/


/* function to calculate different options for the target playout delay */
static void JB4_targetPlayoutDelay( const JB4_HANDLE h, uint32_t *targetMin,
                                    uint32_t *targetMax, uint32_t *targetDtx, uint32_t *targetStartUp )
{
    uint32_t ltJitter, extraDelayReserve;
    /* adapt target delay to partial copy offset */
    extraDelayReserve = 0;
    h->rfDelay = 0;
    if(h->rfOffset7Active != 0)
    {
        h->rfDelay = 140;
    }
    else if(h->rfOffset5Active != 0)
    {
        h->rfDelay = 100;
    }
    else if(h->rfOffset2Active == 0 && h->rfOffset3Active == 0)
    {
        /* keep some delay reserve for RF-off */
        extraDelayReserve = 15;
    }

    /* get estimated long term jitter */
    if( JB4_JMF_Jitter( h->ltJmf, &ltJitter ) == 0 )
    {
        /* combine long term and short term jitter to calculate target delay values */
        *targetMax     = h->stJitter + h->safetyMargin + h->rfDelay;
        *targetMin     = JB4_MIN( ltJitter + 20 + h->rfDelay + extraDelayReserve, *targetMax );
        *targetDtx     = JB4_MIN( ltJitter + extraDelayReserve, h->stJitter );
        *targetStartUp = ( *targetMin + *targetMax + extraDelayReserve / 4) / 2;
    }
    else
    {
        /* combine long term and short term jitter to calculate target delay values */
        *targetMax     = h->safetyMargin;
        *targetMin     = JB4_MIN( 20, *targetMax );
        *targetDtx     = 0;
        *targetStartUp = ( *targetMin + *targetMax ) / 2;
    }
    if(*targetStartUp < 60)
    {
        *targetStartUp = 60;
    }
}

/* function to do playout adaptation before playing the next data unit */
static int JB4_adaptPlayout( JB4_HANDLE h, uint32_t sysTime, uint32_t extBufferedTime,
                             JB4_DATAUNIT_HANDLE *pDataUnit, uint32_t *scale, uint32_t *maxScaling )
{
    bool_t stretchTime;
    /* reset scale */
    if( scale == NULL || maxScaling == NULL )
    {
        return -1;
    }
    *scale = 100;
    *maxScaling = 0;
    stretchTime = false;

    /* switch type of current playout (first one, active, DTX) */
    if( !h->firstDataUnitPopped )
    {
        JB4_adaptFirstPlayout( h, sysTime, &stretchTime );
    }
    else if( h->lastPoppedWasSilence )
    {
        JB4_adaptDtxPlayout( h, sysTime, &stretchTime );
    }
    else
    {
        JB4_adaptActivePlayout( h, sysTime, extBufferedTime, scale, maxScaling );
    }

    /* time shrinking done if needed, now do time stretching or pop data unit to play */
    if( stretchTime )
    {
        /* return empty data unit */
        *pDataUnit = NULL;
        if( h->firstDataUnitPopped )
        {
            ++h->nUnavailablePopped;
            if( !h->lastPoppedWasSilence )
            {
                ++h->nStretched;
                /* jitter-induced insertion (e.g. buffer underflow) */
                ++h->jitterInducedConcealments;
            }
        }
        /* add one frame to last playout delay */
        h->lastPlayoutOffset += h->frameDuration;
    }
    else
    {
        /* return next data unit from buffer */
        JB4_popFromBuffer( h, sysTime, pDataUnit );
    }

    return 0;
}

/* function for playout adaptation while active (no DTX) */
static void JB4_adaptActivePlayout( JB4_HANDLE h, uint32_t sysTime,
                                    uint32_t extBufferedTime, uint32_t *scale, uint32_t *maxScaling )
{
    JB4_DATAUNIT_HANDLE nextDataUnit;
    bool_t convertToLateLoss, dropEarly;
    uint32_t targetMin, targetMax, targetDtx, targetStartUp, targetMaxStretch;
    uint32_t currPlayoutDelay, gap, buffered;
    uint32_t dropGapMax, dropRateMin, dropRateMax, rate;
    int32_t minOffTicks, tsDiffToNextDataUnit;

    JB4_targetPlayoutDelay( h, &targetMin, &targetMax, &targetDtx, &targetStartUp );
    if( JB4_JMF_MinOffset( h->ltJmf, &minOffTicks ) != 0 )
    {
        return;
    }
    h->targetPlayoutDelay = ( targetMin + targetMax ) / 2;

    convertToLateLoss = false;
    dropEarly         = false;
    dropGapMax        = 200;
    dropRateMin       = 5;
    dropRateMax       = 200; /* 20% */

    /* calculate current playout delay */
    currPlayoutDelay = h->lastPlayoutOffset - minOffTicks + extBufferedTime;
    if( !JB4_INPUTBUFFER_IsEmpty( h->inputBuffer ) )
    {
        nextDataUnit = (JB4_DATAUNIT_HANDLE)JB4_INPUTBUFFER_Front( h->inputBuffer );
        tsDiffToNextDataUnit = JB4_rtpTimeStampDiff( h->nextExpectedTs, nextDataUnit->timeStamp );
        if( tsDiffToNextDataUnit < 0 )
        {
            convertToLateLoss = true;
            /* time stretching is expected -> increase playout delay to allow dropping the late frame */
            currPlayoutDelay -= tsDiffToNextDataUnit;
            currPlayoutDelay += 1;
        }
    }

    /*  decided between shrinking/stretching */
    if( currPlayoutDelay > targetMax )   /* time shrinking */
    {
        gap = currPlayoutDelay - h->targetPlayoutDelay;
        /* check if gap is positive and dropping is allowed
         * and buffer contains enough time (ignoring one frame) */
        if( gap > 0 &&
                JB4_inspectBufferForDropping( h, &dropEarly, &buffered ) == 0 &&
                ( convertToLateLoss ||
                  ( buffered + h->frameDuration + extBufferedTime ) > targetMax ) )
        {
            if( convertToLateLoss )
            {
                JB4_dropFromBuffer( h, sysTime );
            }
            else if( dropEarly )
            {
                JB4_dropFromBuffer( h, sysTime );
                ++h->nLostOrStretched;
            }
            else
            {
                /* limit gap to [gapMin,gapMax] and calculate current drop rate from gap */
                rate = JB4_MIN( (uint32_t)(gap), dropGapMax ) *
                       ( dropRateMax - dropRateMin ) / dropGapMax + dropRateMin;
                *scale = ( 1000 - rate ) / 10;
                *maxScaling = currPlayoutDelay - targetMax;
            }
        }
    }
    else   /* time stretching */
    {
        uint32_t delayWithClearedExternalBuffer;
        /* Stretching only makes sense if we win one additional frame in the input buffer.
         * If too much additional delay would be required to do so, then do not scale.
         * Also make sure that the delay doesn't increase too much. */
        delayWithClearedExternalBuffer = currPlayoutDelay - extBufferedTime + h->frameDuration;
        targetMaxStretch = targetMax - h->frameDuration;
        if( delayWithClearedExternalBuffer + h->frameDuration <= targetMaxStretch &&
                currPlayoutDelay < targetMaxStretch && currPlayoutDelay < (uint32_t)(110 + h->rfDelay / 4))
        {
            *scale = 120;
            *maxScaling = targetMaxStretch - currPlayoutDelay;
        }
    }

}

/* function for playout adaptation while DTX */
static void JB4_adaptDtxPlayout( JB4_HANDLE h, uint32_t sysTime, bool_t *stretchTime )
{
    JB4_DATAUNIT_HANDLE firstDu;
    uint32_t firstTs;
    uint32_t targetMin, targetMax, targetDtx, targetStartUp;
    uint32_t currPlayoutDelay, headRoom;
    int32_t minOffTicks, tsDiffToNextDataUnit;

    JB4_targetPlayoutDelay( h, &targetMin, &targetMax, &targetDtx, &targetStartUp );
    if( JB4_JMF_MinOffset( h->ltJmf, &minOffTicks ) != 0 )
    {
        return;
    }

    /* calculate current playout delay */
    currPlayoutDelay = h->lastPlayoutOffset - minOffTicks;

    /* check for startup after DTX */
    if( !JB4_INPUTBUFFER_IsEmpty( h->inputBuffer ) )
    {
        firstDu = (JB4_DATAUNIT_HANDLE)JB4_INPUTBUFFER_Front( h->inputBuffer );
        firstTs = firstDu->timeStamp;

        tsDiffToNextDataUnit = JB4_rtpTimeStampDiff( h->nextExpectedTs, firstTs );
        /* check if the next available data unit should already be used (time stamp order) */
        if( tsDiffToNextDataUnit > 0 )
        {
            /* time stretching is expected -> increase playout delay */
            currPlayoutDelay += tsDiffToNextDataUnit;
        }
        if( !firstDu->silenceIndicator )
        {
            /* recalculate playout delay based on first buffered data unit */
            JB4_playoutDelay( h, sysTime, firstTs, &currPlayoutDelay );
            /* check if the next available data unit should already be used (time stamp order) */
            if( tsDiffToNextDataUnit > 0 )
            {
                /* time stretching is expected -> increase playout delay */
                currPlayoutDelay += tsDiffToNextDataUnit;
            }
            h->targetPlayoutDelay = targetStartUp;
            headRoom = 600 * h->frameDuration / 1000;
            /*  decided between shrinking/stretching */
            if( currPlayoutDelay > targetStartUp + headRoom )   /* time shrinking */
            {
                if( JB4_checkDtxDropping( h ) )
                {
                    JB4_dropFromBuffer( h, sysTime );
                }
            }
            else if( currPlayoutDelay + headRoom < targetStartUp )   /* time stretching */
            {
                *stretchTime = true;
            }
            return;
        }
    }

    /* adapt while DTX */
    h->targetPlayoutDelay = targetDtx;

    /*  decided between shrinking/stretching */
    if( currPlayoutDelay >= targetDtx + h->frameDuration )   /* time shrinking */
    {
        if( JB4_checkDtxDropping( h ) )
        {
            JB4_dropFromBuffer( h, sysTime );
        }
    }
    else if( currPlayoutDelay + 500 * h->frameDuration / 1000 < targetDtx )   /* time stretching */
    {
        *stretchTime = true;
    }

}

/* function to do playout adaptation before playing the first data unit */
static void JB4_adaptFirstPlayout( JB4_HANDLE h, uint32_t sysTime, bool_t *prebuffer )
{
    uint32_t currPlayoutDelay;
    JB4_DATAUNIT_HANDLE firstDu;
    uint32_t targetMin, targetMax, targetDtx, targetStartUp;
    if( JB4_INPUTBUFFER_IsEmpty( h->inputBuffer ) )
    {
        *prebuffer = true;
        return;
    }
    JB4_targetPlayoutDelay( h, &targetMin, &targetMax, &targetDtx, &targetStartUp );
    if(targetStartUp < h->frameDuration)
    {
        return;
    }
    /* calculate delay if first data unit would be played now */
    firstDu = (JB4_DATAUNIT_HANDLE)JB4_INPUTBUFFER_Front( h->inputBuffer );
    if( JB4_playoutDelay( h, sysTime, firstDu->timeStamp, &currPlayoutDelay ) != 0 )
    {
        *prebuffer = true;
        return;
    }
    if( currPlayoutDelay + h->frameDuration / 2 < targetStartUp )   /* time stretching */
    {
        *prebuffer = true;
    }
    else   /* no adaptation, start playout */
    {
        *prebuffer = false;
    }
}

/* function to look into the buffer and check if it makes sense to drop a data unit */
static int JB4_inspectBufferForDropping( const JB4_HANDLE h, bool_t *dropEarly, uint32_t *buffered )
{
    unsigned int inputBufferSize;
    int16_t seqNrDiff;
    int32_t bufferedTs;
    uint32_t firstTs;
    uint64_t beginTs, endTs;
    JB4_DATAUNIT_HANDLE firstDu, secondDu, lastDu;

    assert( !h->lastPoppedWasSilence );
    *dropEarly   = false;
    *buffered    = 0;
    inputBufferSize = JB4_INPUTBUFFER_Size( h->inputBuffer );
    if( inputBufferSize == 0U )
    {
        return -1;
    }

    firstDu = (JB4_DATAUNIT_HANDLE)JB4_INPUTBUFFER_Front( h->inputBuffer );
    firstTs = firstDu->timeStamp;
    /* check for loss: sequence number diff is exactly 0 in the valid case */
    if( h->firstDataUnitPopped )
    {
        seqNrDiff = JB4_rtpTimeStampDiff( h->nextExpectedTs, firstTs ) /
                    (int32_t)(h->frameDuration);
    }
    else
    {
        seqNrDiff = 0;
    }
    if( seqNrDiff <= 0 )
    {
        /* preview data unit to play after dropping */
        if( inputBufferSize <= 1U )
        {
            /* data unit to play missing, avoid drop followed by concealment */
            return -1;
        }
        secondDu = JB4_INPUTBUFFER_Element( h->inputBuffer, 1 );
        if( firstTs + h->frameDuration != secondDu->timeStamp )
        {
            /* data unit to play is not available, avoid drop followed by concealment */
            return -1;
        }
        /* calculate buffered time span */
        bufferedTs = 0;
    }
    else if( seqNrDiff == 2 )
    {
        /* data unit to play is not available, avoid dropping followed by concealment */
        return -1;
    }
    else   /* seqNoDiff == 1 || seqNoDiff > 2 */
    {
        /* first data unit is not available -> drop it early to avoid concealment
         * This is very aggressive: ignores the maximum drop rate (50% drop and 50% concealment for adjacent lost),
         * but on the other hand, dropping sounds better than concealment. */
        *dropEarly = true;
        /* data unit to drop (first one) is lost */
        bufferedTs = 0;
    }

    /* add time stamp difference of last and first actually buffered data unit */
    if( inputBufferSize == 1U )
    {
        bufferedTs += h->frameDuration;
    }
    else
    {
        lastDu = (JB4_DATAUNIT_HANDLE)JB4_INPUTBUFFER_Back( h->inputBuffer );
        beginTs = firstTs;
        endTs   = lastDu->timeStamp + h->frameDuration;
        /* check for RTP time stamp wrap around */
        if( endTs < beginTs )
        {
            endTs = endTs + 0xFFFFFFFF;
        }
        bufferedTs += (int32_t)(endTs - beginTs);
    }

    /* the result should not be negative */
    if( bufferedTs < 0 )
    {
        return -1;
    }
    *buffered = bufferedTs;

    return 0;
}

/* function to look into the buffer and check if it makes sense to drop a data unit */
static int JB4_checkDtxDropping( const JB4_HANDLE h )
{
    unsigned int inputBufferSize;
    int16_t seqNrDiff;
    JB4_DATAUNIT_HANDLE firstDu;
    int droppingAllowed;


    assert( h->firstDataUnitPopped );
    assert( h->lastPoppedWasSilence );
    droppingAllowed = 1;
    inputBufferSize = JB4_INPUTBUFFER_Size( h->inputBuffer );
    if( inputBufferSize > 0U )
    {
        firstDu = (JB4_DATAUNIT_HANDLE)JB4_INPUTBUFFER_Front( h->inputBuffer );
        /* check for loss: sequence number diff is exactly 0 in the valid case */
        seqNrDiff = JB4_rtpTimeStampDiff( h->nextExpectedTs, firstDu->timeStamp ) /
                    (int32_t)(h->frameDuration);
        if( seqNrDiff <= 0 )
        {
            /* no not drop first active frame */
            droppingAllowed = 0;
        }
    }
    /* else: buffer empty, allow dropping FRAME_NO_DATA */

    return droppingAllowed;
}

/* function to estimate the short term jitter */
static void JB4_estimateShortTermJitter( JB4_HANDLE h, uint32_t rcvTime, uint32_t rtpTimeStamp )
{
    uint32_t jitter, duration, maxDuration;
    int32_t minTime, maxTime;
    JB4_CIRCULARBUFFER_ELEMENT maxElement, dequedElement;


    jitter = 0;
    JB4_JMF_PushPacket( h->stJmf, rcvTime, rtpTimeStamp );
    /* save delta delay */
    if( JB4_JMF_Jitter( h->stJmf, &jitter ) == 0 )
    {
        /* compensate difference between both offsets */
        int32_t stOffset, ltOffset;
        JB4_JMF_MinOffset( h->stJmf, &stOffset );
        JB4_JMF_MinOffset( h->ltJmf, &ltOffset );
        jitter += stOffset - ltOffset;
        assert( (int)jitter >= 0 );
        if( JB4_CIRCULARBUFFER_IsFull( h->stJitterFifo ) )
        {
            JB4_CIRCULARBUFFER_Deque( h->stJitterFifo, &dequedElement );
            JB4_CIRCULARBUFFER_Deque( h->stTimeStampFifo, &dequedElement );
        }
        JB4_CIRCULARBUFFER_Enque( h->stJitterFifo, jitter );
        JB4_CIRCULARBUFFER_Enque( h->stTimeStampFifo, rtpTimeStamp );

        /* check for duration and discard first entries if too long */
        minTime = JB4_CIRCULARBUFFER_Front( h->stTimeStampFifo );
        maxTime = JB4_CIRCULARBUFFER_Back( h->stTimeStampFifo );
        if( maxTime > minTime )
        {
            duration = maxTime - minTime;
            maxDuration = 4 * h->timeScale;
            while( duration > maxDuration )
            {
                JB4_CIRCULARBUFFER_Deque( h->stJitterFifo, &dequedElement );
                JB4_CIRCULARBUFFER_Deque( h->stTimeStampFifo, &dequedElement );
                minTime = JB4_CIRCULARBUFFER_Front( h->stTimeStampFifo );
                if( maxTime <= minTime )
                {
                    break;
                }
                duration = maxTime - minTime;
            }
        }
    }

    /* update h->stJitter */
    if( !JB4_CIRCULARBUFFER_IsEmpty( h->stJitterFifo ) )
    {
        JB4_CIRCULARBUFFER_Max( h->stJitterFifo, &maxElement );
        /* round up to full frame duration */
        h->stJitter = (uint32_t)ceil( (double)( maxElement ) / h->frameDuration ) *
                      h->frameDuration;
    }

}

/* function to pop a data unit from the buffer */
static void JB4_popFromBuffer( JB4_HANDLE h, uint32_t sysTime, JB4_DATAUNIT_HANDLE *pDataUnit )
{
    JB4_DATAUNIT_HANDLE nextDataUnit;
    uint32_t nStretched;
    int32_t  tsDiff;

    JB4_DATAUNIT_HANDLE tempDataUnit;
    unsigned int readlen ;
    unsigned short i;
    int frameoffset;
    unsigned int maxval;

    Word32 lost, total_rec ;



    JB4_DATAUNIT_HANDLE partialCopyDu;
    unsigned int searchpos, endpos;


    /* check if a data unit is available */
    if( JB4_INPUTBUFFER_IsEmpty( h->inputBuffer ) )
    {
        /* no data unit available */
        *pDataUnit = NULL;
        h->nextExpectedTs += h->frameDuration;
        if( h->lastPoppedWasSilence )
        {
            ++h->nComfortNoice;
        }
        else
        {
            ++h->nUnavailablePopped;
            ++h->nLostOrStretched;
        }

        return;
    }

    /* preview next data unit in sequence order */
    nextDataUnit = (JB4_DATAUNIT_HANDLE)JB4_INPUTBUFFER_Front( h->inputBuffer );

    /* check if this is the first data unit */
    if( !h->firstDataUnitPopped )
    {
        h->firstDataUnitPopped = true;
        /* adjust sequence numbers to avoid handling first packet as loss */
        h->nextExpectedTs = nextDataUnit->timeStamp;
    }

    /* check if the next available data unit should already be used (time stamp order) */
    tsDiff = JB4_rtpTimeStampDiff( nextDataUnit->timeStamp, h->nextExpectedTs );

    h->totWin += 1;
    if ( ( h->totWin > 3000) || ( h->FecOffWinLen > 100 ) )
    {
        maxval = h->FecOffWin[1];
        h->optimum_offset = 1;
        for( i = 2; i < MAXOFFSET ; i++ )
        {
            if ( h->FecOffWin[i] > maxval )
            {
                maxval = h->FecOffWin[i] ;
                h->optimum_offset = i ;
            }
            h->FecOffWin[i] = 0;
        }
        h->FecOffWin[0] = 0;
        h->FecOffWin[1] = 0;
        h->totWin = 0;
        h->FecOffWinLen = 0;


        lost =h->nLost+ h->nPartialCopiesUsed - h->last_nLost ;
        total_rec =  h->nAvailablePopped + h->nUnavailablePopped - h->last_ntot ;

        if ( lost != 0 && total_rec != 0 )
        {
            h->netLossRate = (float)lost/(float)total_rec;
        }
        else
        {
            h->netLossRate = 0.0f;
        }
        h->last_nLost = L_add(h->nLost, h->nPartialCopiesUsed);
        h->last_ntot = L_add(h->nAvailablePopped , h->nUnavailablePopped);

    }

    if( tsDiff < 0 )
    {
        readlen = JB4_INPUTBUFFER_Size( h->inputBuffer );
        for ( i=0; i < readlen; i++)
        {
            tempDataUnit = (JB4_DATAUNIT_HANDLE)JB4_INPUTBUFFER_Element( h->inputBuffer, i );
            if ( ! tempDataUnit->partial_frame && !h->lastPoppedWasSilence )
            {
                frameoffset = JB4_rtpTimeStampDiff(  h->nextExpectedTs, tempDataUnit->timeStamp )/20 ;

                if ( frameoffset > 0  && frameoffset < MAXOFFSET )
                {
                    h->FecOffWin[frameoffset] += 1;
                }
            }
        }
        h->FecOffWinLen += 1;

        /* next expected data unit is missing
         * -> conceal network loss, do time stretching or create comfort noise */
        *pDataUnit = NULL;

        /* update statistics */
        h->nextExpectedTs += h->frameDuration;
        if( h->lastPoppedWasSilence )
        {
            ++h->nComfortNoice;
        }
        else
        {
            ++h->nUnavailablePopped;
            ++h->nLostOrStretched;
        }
        return;
    }

    /* fetch the next data unit from buffer */
    *pDataUnit = nextDataUnit;
    nextDataUnit->nextCoderType = INACTIVE;
    if ( h->pre_partial_frame || nextDataUnit->partial_frame )
    {
        if ( nextDataUnit->partial_frame )
        {
            h->pre_partial_frame = 1;
        }
        else if ( h->pre_partial_frame )
        {
            h->pre_partial_frame = 0;
        }

        endpos = JB4_INPUTBUFFER_Size(h->inputBuffer);
        for(searchpos = 0; searchpos < endpos; searchpos++)
        {
            partialCopyDu = (JB4_DATAUNIT_HANDLE)JB4_INPUTBUFFER_Element(h->inputBuffer, searchpos);
            if ( partialCopyDu->timeStamp == nextDataUnit->timeStamp + partialCopyDu->duration )
            {
                get_NextCoderType( partialCopyDu->data, &nextDataUnit->nextCoderType);
                break;
            }
        }
    }
    JB4_INPUTBUFFER_Deque( h->inputBuffer, (void**)pDataUnit );

    if ( nextDataUnit->partial_frame )
    {
        h->nPartialCopiesUsed += 1;

        readlen = JB4_INPUTBUFFER_Size( h->inputBuffer );
        for ( i=0; i < readlen; i++)
        {
            tempDataUnit = (JB4_DATAUNIT_HANDLE)JB4_INPUTBUFFER_Element( h->inputBuffer, i );
            if ( ! tempDataUnit->partial_frame && !h->lastPoppedWasSilence )
            {
                frameoffset = JB4_rtpTimeStampDiff(  h->nextExpectedTs, tempDataUnit->timeStamp )/20 ;

                if ( frameoffset > 0  && frameoffset < MAXOFFSET )
                {
                    h->FecOffWin[frameoffset] += 1;
                }
            }
        }
        h->FecOffWinLen += 1;
    }

    /* update statistics */
    if( h->nLostOrStretched != 0U )
    {
        assert( h->lastPoppedWasSilence == false );
        /* separate concealments since last available pop in lost and stretched */
        nStretched = tsDiff / h->frameDuration;
        assert( h->nLostOrStretched >= nStretched );
        h->nLost                     += h->nLostOrStretched - nStretched;
        /* jitter-induced insertion (e.g. buffer underflow) */
        h->jitterInducedConcealments += nStretched;
        h->nStretched += nStretched;
        h->nLostOrStretched = 0;
    }
    h->lastReturnedTs        = nextDataUnit->timeStamp;
    JB4_updateLastTimingMembers( h, sysTime, nextDataUnit->timeStamp );
    h->nextExpectedTs        = nextDataUnit->timeStamp + h->frameDuration;
    if( nextDataUnit->silenceIndicator )
    {
        h->lastPoppedWasSilence = true;
        ++h->nComfortNoice;
    }
    else
    {
        h->lastPoppedWasSilence = false;
        ++h->nAvailablePopped;
    }
}

/* function to drop a data unit from the buffer - updates nShrinked */
static void JB4_dropFromBuffer( JB4_HANDLE h, uint32_t sysTime )
{
    JB4_DATAUNIT_HANDLE nextDataUnit, dataUnit;
    int32_t  tsDiff;
    uint32_t nStretched;
    (void)sysTime;

    /* check if a data unit is available */
    if( JB4_INPUTBUFFER_IsEmpty( h->inputBuffer ) )
    {
        return;
    }
    /* preview next data unit in sequence order */
    nextDataUnit = (JB4_DATAUNIT_HANDLE)JB4_INPUTBUFFER_Front( h->inputBuffer );

    /* check if this is the first data unit */
    if( !h->firstDataUnitPopped )
    {
        h->firstDataUnitPopped = true;
        /* adjust sequence numbers to avoid handling first packet as loss */
        h->nextExpectedTs = nextDataUnit->timeStamp;
    }

    /* check if the next available data unit should already be used (time stamp order) */
    tsDiff = JB4_rtpTimeStampDiff( nextDataUnit->timeStamp, h->nextExpectedTs );
    if( tsDiff < 0 )
    {
        /* next expected data unit is missing, remember this data unit as popped,
         * but do not count it as lost, because it will not be concealed */
        h->nextExpectedTs += h->frameDuration;
        /* substract one frame from last playout delay */
        h->lastPlayoutOffset -= h->frameDuration;
        if( !h->lastPoppedWasSilence )
        {
            ++h->nShrinked;
            /* modification of the output timeline due to link loss */
            ++h->nUnavailablePopped;
            ++h->nLostOrStretched;
        }
        if( h->lastTargetTime != 0U )
        {
            h->lastTargetTime += h->frameDuration;
        }
        return;
    }

    /* fetch the next data unit from buffer */
    JB4_INPUTBUFFER_Deque( h->inputBuffer, (void *)&dataUnit );
    /* update statistics */
    if( h->nLostOrStretched != 0U )
    {
        assert( h->lastPoppedWasSilence == false );
        /* separate concealments since last available pop in lost and stretched */
        nStretched = tsDiff / h->frameDuration;
        assert( h->nLostOrStretched >= nStretched );

        /* convert stretching followed by shrinking to late-loss */
        if( nStretched > 0U )
        {
            --nStretched;
            ++h->nLateLost;
            h->nLost += h->nLostOrStretched - nStretched;
            /* jitter-induced insertion (e.g. buffer underflow) */
            h->jitterInducedConcealments += nStretched;
            if( !dataUnit->silenceIndicator )
            {
                /* JBM induced removal of a speech frame (intentional frame dropping) */
                ++h->jitterInducedConcealments;
            }
            h->nStretched += nStretched;
        }
        else
        {
            h->nLost += h->nLostOrStretched;
            ++h->nShrinked;
            if( !dataUnit->silenceIndicator )
            {
                /* JBM induced removal of a speech frame (intentional frame dropping) */
                ++h->jitterInducedConcealments;
            }
        }
        h->nLostOrStretched = 0;
    }
    else
    {
        if( !dataUnit->silenceIndicator )
        {
            ++h->nShrinked;
            /* JBM induced removal of a speech frame (intentional frame dropping) */
            ++h->jitterInducedConcealments;
        }
    }

    h->lastReturnedTs        = dataUnit->timeStamp;
    h->lastPoppedWasSilence  = dataUnit->silenceIndicator;
    h->nextExpectedTs        = dataUnit->timeStamp + h->frameDuration;

    /* substract one frame from last playout delay */
    h->lastPlayoutOffset -= h->frameDuration;
    if( h->lastTargetTime != 0U )
        h->lastTargetTime += h->frameDuration;
    JB4_FreeDataUnit(h, dataUnit);
}

/* function to calculate the playout delay based on the current jitter */
static int JB4_playoutDelay( const JB4_HANDLE h, uint32_t playTime, uint32_t rtpTimeStamp, uint32_t *delay )
{
    int32_t minOffTicks;

    if( JB4_JMF_MinOffset( h->ltJmf, &minOffTicks ) != 0 )
    {
        return -1;
    }

    *delay = playTime - minOffTicks - rtpTimeStamp;

    return 0;
}

/* function to update lastPlayoutDelay and lastTargetTime after popFromBuffer() */
static void JB4_updateLastTimingMembers( JB4_HANDLE h, uint32_t playTime,
        uint32_t rtpTimeStamp )
{
    int32_t minOffTicks;

    if( JB4_JMF_MinOffset( h->ltJmf, &minOffTicks ) != 0 )
    {
        return;
    }

    /* playoutDelay = playTime - minOffset - timeStamp */
    h->lastPlayoutOffset = playTime - rtpTimeStamp;
    /* targetTime = minOffset + timeStamp + targetDelay */
    h->lastTargetTime = (uint32_t)( minOffTicks + rtpTimeStamp + h->targetPlayoutDelay );

}

/* function to compare the RTP time stamps of two data units: newElement==arrayElement ? 0 : (newElement>arrayElement ? +1 : -1) */
static int JB4_inputBufferCompareFunction( const JB4_INPUTBUFFER_ELEMENT newElement,
        const JB4_INPUTBUFFER_ELEMENT arrayElement, bool_t *replaceWithNewElementIfEqual )
{
    JB4_DATAUNIT_HANDLE newDataUnit, arrayDataUnit;
    int32_t diff;
    int result;

    *replaceWithNewElementIfEqual = 0;
    newDataUnit   = (JB4_DATAUNIT_HANDLE)newElement;
    arrayDataUnit = (JB4_DATAUNIT_HANDLE)arrayElement;
    diff = JB4_rtpTimeStampDiff( arrayDataUnit->timeStamp, newDataUnit->timeStamp );
    if( diff > 0 )
    {
        result = 1;
    }
    else if( diff < 0 )
    {
        result = -1;
    }
    else   /* equal timestamps */
    {
        result = 0;
        if(newDataUnit->partial_frame == 0 && arrayDataUnit->partial_frame == 1)
        {
            /* replace partial copy with primary copy */
            *replaceWithNewElementIfEqual = 1;
        }
        else if(newDataUnit->partial_frame == arrayDataUnit->partial_frame &&
                newDataUnit->dataSize > arrayDataUnit->dataSize)
        {
            /* if both are primary or partial: take the one with higher size (e.g. higher bitrate) */
            *replaceWithNewElementIfEqual = 1;
        }
    }
    return result;
}
