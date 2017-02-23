/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

/** \file jbm_jb4_jmf.h jitter measure fifo - a fifo used for windowed measure of network status */

#ifndef JBM_JB4_JMF_H
#define JBM_JB4_JMF_H JBM_JB4_JMF_H

#include "jbm_types.h"

/** handle for jitter measure fifo - a fifo used for windowed measure of network status */
typedef struct JB4_JMF *JB4_JMF_HANDLE;

/**@name functions to manage the fifo */
/**@{ */
int JB4_JMF_Create( JB4_JMF_HANDLE *ph );
void JB4_JMF_Destroy( JB4_JMF_HANDLE *ph );
/** function to set the window size of the fifo and the fraction which will be considered */
/** @param timeScale scale of system time and RTP time stamps
 *  @param windowSize the window size of the fifo in number of packets
 *  @param windowDuration the window size of the fifo as time in sysTimeScale
 *  @param consideredFraction the considered fraction in 1/1000 units, e.g. 900 ignores 10% of the highest samples
 *  @return 0 on success */
int JB4_JMF_Init( JB4_JMF_HANDLE h, int timeScale, unsigned int windowSize,
                  unsigned int windowDuration, unsigned int consideredFraction );
/**@} */

/**@name functions to push packets and get the current jitter rate */
/**@{ */
/** function to calculate jitter for the current packet */
int JB4_JMF_PushPacket( JB4_JMF_HANDLE h, uint32_t sysTime, uint32_t rtpTimeStamp );
/** function to get the current jitter */
int JB4_JMF_Jitter( const JB4_JMF_HANDLE h, uint32_t *jitter );
/** function to get the minimum offset between received time and time stamp of all entries in the fifo */
/*! This value is the offset of the fastest transmitted packet of all packets currently
 *  contained in the fifo.
 *  @param[out] offset the minimum offset in microseconds */
int JB4_JMF_MinOffset( const JB4_JMF_HANDLE h, int32_t *offset );
/**@} */

#endif /* JBM_JB4_JMF_H */
