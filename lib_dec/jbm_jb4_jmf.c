/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

/** \file jbm_jb4_jmf.cpp jitter measure fifo - a fifo used for windowed measure of network status */

/* system includes */
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include "options.h"
/* local includes */
#include "jbm_jb4_jmf.h"
#include "jbm_jb4_circularbuffer.h"
/* instrumentation */


/** jitter measure fifo - a fifo used for windowed measure of network status */
struct JB4_JMF
{
    /** scale of system time and RTP time stamps */
    int                       timeScale;
    /** the window size of the fifo as time in sysTimeScale */
    unsigned int              maxWindowDuration;
    /** considered fraction in 1/1000 units, e.g. 900 ignores 10% of the highest samples */
    unsigned int              consideredFraction;

    /** fifo containing the delay entries (ordered by receive time) */
    JB4_CIRCULARBUFFER_HANDLE fifo;
    /** fifo containing the offset entries (ordered by receive time) */
    JB4_CIRCULARBUFFER_HANDLE offsetFifo;
    /** fifo containing the RTP times of the values in offsetFifo (ordered by receive time) */
    JB4_CIRCULARBUFFER_HANDLE timeStampFifo;
    /** flag if the first packet was already pushed */
    int                       firstPacketPushed;
    /** last packets system time in microseconds */
    int32_t                   lastSysTime;
    /** RTP time stamp of the last pushed packet */
    int32_t                   lastRtpTimeStamp;
    /** last packets calculated delay value */
    int32_t                   lastDelay;
    /** number of elements to ignore for percentile calculation - value set within init */
    int16_t                   nElementsToIgnore;
};


/** helper function to add an entry at back of the buffer */
static void JB4_JMF_pushBack( JB4_JMF_HANDLE h, int32_t delay, int32_t offset, uint32_t time );
/** helper function to remove an entry from the front of the buffer */
static void JB4_JMF_popFront( JB4_JMF_HANDLE h );


int JB4_JMF_Create( JB4_JMF_HANDLE *ph )
{
    JB4_JMF_HANDLE h = malloc( sizeof( struct JB4_JMF ) );


    JB4_CIRCULARBUFFER_Create( &h->fifo );
    JB4_CIRCULARBUFFER_Create( &h->offsetFifo );
    JB4_CIRCULARBUFFER_Create( &h->timeStampFifo );

    h->timeScale             = 1000;
    h->consideredFraction    = 1000;
    h->firstPacketPushed     = 0;
    h->lastSysTime           = 0;
    h->lastRtpTimeStamp      = 0;
    h->lastDelay             = 0;
    h->nElementsToIgnore     = 0;

    *ph = h;

    return 0;
}

void JB4_JMF_Destroy( JB4_JMF_HANDLE *ph )
{
    JB4_JMF_HANDLE h;

    if( !ph )
    {
        return;
    }
    h = *ph;
    if( !h )
    {
        return;
    }

    JB4_CIRCULARBUFFER_Destroy( &h->fifo );
    JB4_CIRCULARBUFFER_Destroy( &h->offsetFifo );
    JB4_CIRCULARBUFFER_Destroy( &h->timeStampFifo );

    free( h );
    *ph = NULL;
}

/* function to set the window size of the fifo and the fraction which will be considered */
int JB4_JMF_Init( JB4_JMF_HANDLE h, int timeScale, unsigned int windowSize,
                  unsigned int windowDuration, unsigned int consideredFraction )
{

    /* check parameters */
    if( windowSize != 0U && consideredFraction * windowSize / 1000 < 2 )
    {
        return -1;
    }
    if( consideredFraction > 1000 )
    {
        return -1;
    }

    /* store values */
    h->timeScale          = timeScale;
    h->maxWindowDuration  = windowDuration;
    h->consideredFraction = consideredFraction;

    JB4_CIRCULARBUFFER_Init( h->fifo, windowSize );
    JB4_CIRCULARBUFFER_Init( h->offsetFifo, windowSize );
    JB4_CIRCULARBUFFER_Init( h->timeStampFifo, windowSize );

    h->nElementsToIgnore = windowSize * ( 1000 - consideredFraction ) / 1000;
    return 0;
}

/* function to calculate delay for the current packet */
int JB4_JMF_PushPacket( JB4_JMF_HANDLE h, uint32_t sysTime, uint32_t rtpTimeStamp )
{
    int32_t rtpTimeDiff, sysTimeDiff;
    int32_t offset, delay;


    /* check if this is the first entry */
    if( h->firstPacketPushed == 0 )
    {
        h->firstPacketPushed = 1;
        h->lastSysTime       = sysTime;
        h->lastRtpTimeStamp  = rtpTimeStamp;
        return 0;
    }

    rtpTimeDiff = rtpTimeStamp - h->lastRtpTimeStamp;
    sysTimeDiff = sysTime - h->lastSysTime;
    offset      = sysTime - rtpTimeStamp;

    /* get the delay (yes, signed!!!!) */
    delay = sysTimeDiff - rtpTimeDiff + h->lastDelay;

    /* remember old values */
    h->lastSysTime      = sysTime;
    h->lastRtpTimeStamp = rtpTimeStamp;
    /* reset delay if absolute value is greater than 60s
     * to avoid overflow caused by clockdrift */
    if( delay > 60 * h->timeScale || delay < -60 * h->timeScale )
    {
        h->lastDelay = 0;
    }
    else
    {
        h->lastDelay = delay;
    }

    JB4_JMF_pushBack( h, delay, offset, rtpTimeStamp );

    return 0;
}

/* function to get the current jitter */
int JB4_JMF_Jitter( const JB4_JMF_HANDLE h, uint32_t *jitter )
{
    JB4_CIRCULARBUFFER_ELEMENT min, percentile;

    /* sanity check (must not be empty) and return invalid result if there is only one entry */
    if( JB4_CIRCULARBUFFER_Size( h->fifo ) < 2U )
    {
        return -1;
    }

    JB4_CIRCULARBUFFER_MinAndPercentile( h->fifo, h->nElementsToIgnore, &min, &percentile );

    /* return the difference between the highest considered and the smallest value */
    *jitter = percentile - min;
    assert( percentile >= min );

    return 0;
}

/* function to get the minimum offset between received time and time stamp of all entries in the fifo */
int JB4_JMF_MinOffset( const JB4_JMF_HANDLE h, int32_t *offset )
{
    JB4_CIRCULARBUFFER_ELEMENT min;

    if( JB4_CIRCULARBUFFER_IsEmpty( h->offsetFifo ) )
    {
        return -1;
    }

    JB4_CIRCULARBUFFER_Min( h->offsetFifo, &min );

    *offset = min;

    return 0;
}


/*****************************************************************************
 **************************** private functions ******************************
 *****************************************************************************/

/* helper function to add entry at back of the buffer */
static void JB4_JMF_pushBack( JB4_JMF_HANDLE h, int32_t delay, int32_t offset, uint32_t time )
{
    int32_t minTime, maxTime;
    uint32_t duration;


    /* check for size and discard first entry if too big */
    if( JB4_CIRCULARBUFFER_IsFull( h->fifo ) )
    {
        JB4_JMF_popFront( h );
    }

    /* push back new entry */
    JB4_CIRCULARBUFFER_Enque( h->fifo, delay );
    JB4_CIRCULARBUFFER_Enque( h->offsetFifo, offset );
    JB4_CIRCULARBUFFER_Enque( h->timeStampFifo, time );

    /* check for duration and discard first entries if too long */
    minTime = JB4_CIRCULARBUFFER_Front( h->timeStampFifo );
    maxTime = JB4_CIRCULARBUFFER_Back( h->timeStampFifo );
    if( maxTime > minTime )
    {
        duration = maxTime - minTime;
        while( duration > h->maxWindowDuration )
        {
            JB4_JMF_popFront( h );
            minTime = JB4_CIRCULARBUFFER_Front( h->timeStampFifo );
            if( maxTime <= minTime )
            {
                break;
            }
            duration = maxTime - minTime;
        }
    }

}

/* helper function to remove an entry from the front of the buffer */
static void JB4_JMF_popFront( JB4_JMF_HANDLE h )
{
    JB4_CIRCULARBUFFER_ELEMENT tmpElement;


    /* try to remove one element - fails if empty */
    if( JB4_CIRCULARBUFFER_Deque( h->fifo, &tmpElement ) != 0 )
    {
        return;
    }
    /* also remove offset entry */
    JB4_CIRCULARBUFFER_Deque( h->offsetFifo, &tmpElement );
    JB4_CIRCULARBUFFER_Deque( h->timeStampFifo, &tmpElement );

}

