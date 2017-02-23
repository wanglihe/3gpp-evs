/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#ifndef EvsRXLIB_H
#define EvsRXLIB_H

/* local headers */
#include "stat_dec.h"

/*
 * ENUMS
 */

/* Receiver error enums */
typedef enum _EVS_RX_ERROR
{
    EVS_RX_NO_ERROR          = 0x0000,
    EVS_RX_MEMORY_ERROR      = 0x0001,
    EVS_RX_WRONG_PARAMS      = 0x0002,
    EVS_RX_INIT_ERROR        = 0x0003,
    EVS_RX_RECEIVER_ERROR    = 0x0004,
    EVS_RX_DECODER_ERROR     = 0x0005,
    EVS_RX_JBM_ERROR         = 0x0006,
    EVS_RX_TIMESCALER_ERROR  = 0x0007,
    EVS_RX_NOT_IMPLEMENTED   = 0x0010

} EVS_RX_ERROR;


/*
 * Structures
 */

typedef struct EVS_RX*     EVS_RX_HANDLE;

/*
 * Functions
 */

/*! Opens the EVS Receiver instance. */
EVS_RX_ERROR
EVS_RX_Open(EVS_RX_HANDLE* phEvsRX,
            Decoder_State *st,
            Word16 jbmSafetyMargin);

/*! Sets the name of the JBM trace file which will be created. */
EVS_RX_ERROR
EVS_RX_SetJbmTraceFileName(EVS_RX_HANDLE hEvsRX,
                           const char *jbmTraceFileName);

/*! Feeds one frame into the receiver. */
EVS_RX_ERROR
EVS_RX_FeedFrame(EVS_RX_HANDLE hEvsRX,
                 unsigned char *au,
                 unsigned int   auSize,
                 unsigned short rtpSequenceNumber,
                 unsigned long  rtpTimeStamp,
                 unsigned int   rcvTime_ms);

/*! Retrieves one frame of output PCM data. */
EVS_RX_ERROR
EVS_RX_GetSamples(EVS_RX_HANDLE hEvsRX,
                  unsigned int* nOutSamples,
                  Word16       *pcmBuf,
                  unsigned int  pcmBufSize,
                  unsigned int  systemTimestamp_ms
                 );

EVS_RX_ERROR
EVS_RX_Get_FEC_offset( EVS_RX_HANDLE hEvsRX, short* offset
                       , short *FEC_hi
                     );


/*! Returns 1 if the jitter buffer is empty, otherwise 0. */
/*  Intended for flushing at the end of the main loop but not during normal operation! */
unsigned int
EVS_RX_IsEmpty(EVS_RX_HANDLE hEvsRX );

/*! Closes the receiver instance. */
EVS_RX_ERROR
EVS_RX_Close(EVS_RX_HANDLE* phEvsRX );

#endif
