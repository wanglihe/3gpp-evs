/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#ifndef G192_H
#define G192_H G192_H

#include <stdio.h>

/*
 * ENUMS
 */

/* error enums */

typedef enum _G192_ERROR
{
    G192_NO_ERROR          = 0x0000,
    G192_MEMORY_ERROR      = 0x0001,
    G192_WRONG_PARAMS      = 0x0002,
    G192_INIT_ERROR        = 0x0003,
    G192_WRITE_ERROR       = 0x0004,
    G192_READ_ERROR        = 0x0005,
    G192_FILE_NOT_FOUND    = 0x0006,
    G192_NOT_IMPLEMENTED   = 0x0010,
    G192_NOT_INITIALIZED   = 0x0100,
    G192_UNKNOWN_ERROR     = 0x1000
} G192_ERROR;

/*
 * Structures
 */

/* main handle */
struct __G192;
typedef struct __G192 * G192_HANDLE;

/*
 * Functions
 */

G192_ERROR
G192_Reader_Open(G192_HANDLE* phG192, FILE * filename);

G192_ERROR
G192_ReadVoipFrame_compact(G192_HANDLE const hG192,
                           unsigned char * const serial,
                           short * const num_bits,
                           unsigned short * const rtpSequenceNumber,
                           unsigned int * const rtpTimeStamp,
                           unsigned int * const rcvTime_ms);

G192_ERROR
G192_ReadVoipFrame_short(G192_HANDLE const hG192,
                         short * const serial,
                         short * const num_bits,
                         unsigned short * const rtpSequenceNumber,
                         unsigned int * const rtpTimeStamp,
                         unsigned int * const rcvTime_ms);

G192_ERROR
G192_Reader_Close(G192_HANDLE* phG192);

#endif /* G192_H */
