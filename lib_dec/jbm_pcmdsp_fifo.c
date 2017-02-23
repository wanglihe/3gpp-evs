/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

/*! @file jbm_pcmdsp_fifo.c Ringbuffer (FIFO) with fixed capacity for audio samples */

/* system headers */
#include <stdlib.h>
#include <string.h>
#include "options.h"
/* instrumentation */
/* local headers */
#include "jbm_pcmdsp_fifo.h"


/** Ringbuffer (FIFO) with fixed capacity for audio samples. */
struct PCMDSP_FIFO
{
    /** number of currently stored samples per channel */
    unsigned int size;
    /** maximum allowed number of samples per channel */
    unsigned int capacity;
    /** sample size in bytes per channel */
    unsigned int nBytesPerSampleSet;

    /** begin of the FIFO data (pointer to bytes) */
    uint8_t *dataBegin;
    /** end of the FIFO data (pointer to bytes) */
    uint8_t *dataEnd;
    /** position of next write operation (pointer to bytes) */
    uint8_t *dataWriteIterator;
    /** position of next read operation (pointer to bytes) */
    uint8_t *dataReadIterator;
};


/* Creates a FIFO. */
int pcmdsp_fifo_create( PCMDSP_FIFO_HANDLE *ph )
{
    PCMDSP_FIFO_HANDLE h = malloc( sizeof( struct PCMDSP_FIFO ) );

    h->size               = 0;
    h->capacity           = 0;
    h->nBytesPerSampleSet = 0;
    h->dataBegin          = NULL;
    h->dataEnd            = NULL;
    h->dataWriteIterator  = NULL;
    h->dataReadIterator   = NULL;

    *ph = h;

    return 0;
}

/* Destroys the FIFO. */
void pcmdsp_fifo_destroy( PCMDSP_FIFO_HANDLE *ph )
{
    PCMDSP_FIFO_HANDLE h;

    if( !ph )
        return;
    h = *ph;
    if( !h )
        return;

    if( h->dataBegin )
        free( h->dataBegin );
    free( h );
    *ph = NULL;

}

/* Initializes the FIFO with a fixed maximum allowed number audio samples. */
int pcmdsp_fifo_init( PCMDSP_FIFO_HANDLE h, unsigned int nSamples,
                      unsigned int nChannels, unsigned int nBytesPerSample )
{
    unsigned int nDataBytes;

    h->capacity           = nSamples;
    h->nBytesPerSampleSet = nChannels * nBytesPerSample;
    nDataBytes            = nSamples * h->nBytesPerSampleSet;
    h->dataBegin          = malloc(nDataBytes);
    h->dataEnd            = h->dataBegin + nDataBytes;
    h->dataWriteIterator  = h->dataBegin;
    h->dataReadIterator   = h->dataBegin;

    return 0;
}

/* Writes the given audio data to the FIFO. */
int pcmdsp_fifo_write( PCMDSP_FIFO_HANDLE h, const uint8_t *samples, unsigned int nSamplesPerChannel )
{
    unsigned int nBytesToWrite;

    /* check for empty input buffer */
    if( nSamplesPerChannel == 0U )
        return 0;
    /* check, if enough space left */
    if( nSamplesPerChannel > h->capacity - h->size )
        return -1;

    nBytesToWrite = nSamplesPerChannel * h->nBytesPerSampleSet;
    if( h->dataWriteIterator + nBytesToWrite > h->dataEnd )
    {
        /* wrap around: writing two parts */
        unsigned int bytesOfFirstPart, secondSize;
        bytesOfFirstPart = h->dataEnd - h->dataWriteIterator;
        secondSize       = nBytesToWrite - bytesOfFirstPart;
        memcpy( h->dataWriteIterator, samples, bytesOfFirstPart );
        memcpy( h->dataBegin, samples + bytesOfFirstPart, secondSize );
        h->dataWriteIterator = h->dataBegin + secondSize;
    }
    else
    {
        /* no wrap around: simple write */
        memcpy( h->dataWriteIterator, samples, nBytesToWrite );
        h->dataWriteIterator += nBytesToWrite;
    }
    h->size += nSamplesPerChannel;

    return 0;
}

/* Reads the given number of audio samples from the FIFO. */
int pcmdsp_fifo_read( PCMDSP_FIFO_HANDLE h, unsigned int nSamplesPerChannel, uint8_t *samples )
{
    unsigned int nBytesToRead;

    /* check for empty output buffer */
    if( nSamplesPerChannel == 0U )
        return 0;
    /* check, if enough bytes readable */
    if( nSamplesPerChannel > h->size )
        return -1;

    nBytesToRead = nSamplesPerChannel * h->nBytesPerSampleSet;
    if( h->dataReadIterator + nBytesToRead > h->dataEnd )
    {
        /* wrap around: reading two parts */
        unsigned int bytesOfFirstPart, nBytesOfSecondPart;
        bytesOfFirstPart   = h->dataEnd - h->dataReadIterator;
        nBytesOfSecondPart = nBytesToRead - bytesOfFirstPart;
        memcpy( samples, h->dataReadIterator, bytesOfFirstPart );
        memcpy( samples + bytesOfFirstPart, h->dataBegin, nBytesOfSecondPart );
        h->dataReadIterator = h->dataBegin + nBytesOfSecondPart;
    }
    else
    {
        /* no wrap around: simple read */
        memcpy( samples, h->dataReadIterator, nBytesToRead );
        h->dataReadIterator += nBytesToRead;
    }
    h->size -= nSamplesPerChannel;

    return 0;
}

/* Returns the number of samples per channel that can be read (number of currently stored samples). */
unsigned int pcmdsp_fifo_nReadableSamples( const PCMDSP_FIFO_HANDLE h )
{
    return h->size;
}

