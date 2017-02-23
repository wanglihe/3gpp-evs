/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

/*! @file jbm_pcmdsp_fifo.h Ringbuffer (FIFO) with fixed capacity for audio samples. */

#ifndef JBM_PCMDSP_FIFO_H
#define JBM_PCMDSP_FIFO_H JBM_PCMDSP_FIFO_H

/* local headers */
#include "jbm_types.h"

/** handle for FIFO with fixed capacity */
typedef struct PCMDSP_FIFO *PCMDSP_FIFO_HANDLE;

/** Creates a FIFO.
 * @param[out] ph pointer to created handle
 * @return 0 if succeeded */
int pcmdsp_fifo_create( PCMDSP_FIFO_HANDLE *ph );
/** Destroys the FIFO. */
void pcmdsp_fifo_destroy( PCMDSP_FIFO_HANDLE *ph );
/** Initializes the FIFO with a fixed maximum allowed number of audio samples.
 * @param[in] nSamples maximum allowed number of samples per channel (capacity)
 * @param[in] nChannels number of audio channels
 * @param[in] nBytesPerSample size in bytes per sample per channel
 * @return 0 if succeeded */
int pcmdsp_fifo_init( PCMDSP_FIFO_HANDLE h, unsigned int nSamples,
                      unsigned int nChannels, unsigned int nBytesPerSample );

/** Writes the given audio data to the FIFO.
 * @param[in] samples pointer to audio samples to append
 * @param[in] nSamplesPerChannel the number of samples per channel to append
 * @return 0 if succeeded */
int pcmdsp_fifo_write( PCMDSP_FIFO_HANDLE h, const uint8_t *samples, unsigned int nSamplesPerChannel );
/** Reads the given number of audio samples from the FIFO.
 * @param[in] nSamplesPerChannel the number of samples per channel to read from the FIFO
 * @param[in] samples pointer where the audio samples will be copied to
 * @return 0 if succeeded */
int pcmdsp_fifo_read( PCMDSP_FIFO_HANDLE h, unsigned int nSamplesPerChannel, uint8_t *samples );

/** Returns the number of samples per channel that can be read (number of currently stored samples). */
unsigned int pcmdsp_fifo_nReadableSamples( const PCMDSP_FIFO_HANDLE h );

#endif /* JBM_PCMDSP_FIFO_H */
