/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

/** \file jbm_jb4sb.h EVS Jitter Buffer Management Interface */

#ifndef JBM_JB4SB_H
#define JBM_JB4SB_H JBM_JB4SB_H

#include "jbm_types.h"

/** handle for jitter buffer */
typedef struct JB4* JB4_HANDLE;

/** jitter buffer data units (access unit together with RTP seqNo, timestamp, ...) */
struct JB4_DATAUNIT
{
    /** the RTP sequence number (16 bits) */
    uint16_t sequenceNumber;
    /** the RTP time stamp (32 bits) of this chunk in timeScale() units */
    uint32_t timeStamp;
    /** the duration of this chunk in timeScale() units */
    uint32_t duration;
    /** the RTP time scale, which is used for timeStamp() and duration() */
    uint32_t timeScale;
    /** the receive time of the RTP packet in milliseconds */
    uint32_t rcvTime;
    /** true, if the data unit contains only silence */
    bool_t   silenceIndicator;

    /** the binary encoded access unit */
    uint8_t *data;
    /** the size of the binary encoded access unit [bits] */
    uint32_t dataSize;

    /** identify if the data unit has a partial copy of a previous frame */
    bool_t partial_frame;
    /** offset of the partial copy contained in that frame or zero */
    int16_t partialCopyOffset;
    int16_t nextCoderType;
};
/** handle for jitter buffer data units */
typedef struct JB4_DATAUNIT* JB4_DATAUNIT_HANDLE;


int JB4_Create( JB4_HANDLE *ph );
void JB4_Destroy( JB4_HANDLE *ph );

int JB4_Init( JB4_HANDLE h, Word16 safetyMargin );

/** Returns a memory slot to store a new data unit */
JB4_DATAUNIT_HANDLE JB4_AllocDataUnit( JB4_HANDLE h );
/** Notifies the JBM that a data unit is no longer used and the memory can be reused */
void JB4_FreeDataUnit( JB4_HANDLE h, JB4_DATAUNIT_HANDLE dataUnit );

int JB4_PushDataUnit( JB4_HANDLE h, JB4_DATAUNIT_HANDLE dataUnit, uint32_t rcvTime );
int JB4_PopDataUnit( JB4_HANDLE h, uint32_t sysTime, uint32_t extBufferedTime,
                     JB4_DATAUNIT_HANDLE *pDataUnit, uint32_t *scale, uint32_t *maxScaling );

int JB4_getFECoffset(JB4_HANDLE h);

short JB4_FECoffset(JB4_HANDLE h);


/** function to get the number of data units contained in the buffer */
unsigned int JB4_bufferedDataUnits( const JB4_HANDLE h );

#endif /* JBM_JB4SB_H */
