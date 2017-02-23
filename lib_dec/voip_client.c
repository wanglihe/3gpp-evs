/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include "options.h"
#include "prot.h"
#include "cnst.h"
#include "rom_com.h"
#include "EvsRXlib.h"
#include "g192.h"


/*------------------------------------------------------------------------------------------*
 * Global variables
 *------------------------------------------------------------------------------------------*/
extern long frame;                 /* Counter of frames */

/*------------------------------------------------------------------------------------------*
 * Local constants
 *------------------------------------------------------------------------------------------*/

#define PCMBUFSIZE            1920

/*------------------------------------------------------------------------------------------*
 * decodeVoip()
 *
 * Main function for EVS decoder with VOIP mode
 *------------------------------------------------------------------------------------------*/

int decodeVoip(
    Decoder_State *st,
    FILE *f_stream,
    FILE *f_synth,
#ifdef SUPPORT_JBM_TRACEFILE
    const char *jbmTraceFileName,
#endif
    const char *jbmFECoffsetFileName, /* : Output file  for Optimum FEC offset        */
    const short quietMode
)
{
    /* input/output */
    G192_HANDLE  g192 = NULL;
    G192_ERROR   g192err;

    short optimum_offset,FEC_hi;
    FILE *f_offset = NULL;

    /* main loop */
    unsigned int nextPacketRcvTime_ms = 0;
    unsigned int systemTime_ms = 0;

    EVS_RX_HANDLE hRX = NULL;
    EVS_RX_ERROR rxerr = EVS_RX_NO_ERROR;
    Word16 jbmSafetyMargin = 60; /* allowed delay reserve in addition to network jitter to reduce late-loss [milliseconds] */
    short dec_delay, zero_pad;

    unsigned char   au[2560];
    short           auSize;
    unsigned short  rtpSequenceNumber;
    unsigned int    rtpTimeStamp;

    Word16       pcmBuf[3 * L_FRAME48k] = {0};
    unsigned int pcmBufSize  = 3 * L_FRAME48k;

    /* open input file */
    g192err = G192_Reader_Open(&g192, f_stream);
    if(g192err != G192_NO_ERROR)
    {
        fprintf(stderr,"error in G192_Reader_Open(): %d\n", g192err);
        return -1;
    }


    if(jbmFECoffsetFileName)
    {
        f_offset =  fopen( jbmFECoffsetFileName, "w+" );
        if(f_offset == NULL)
        {
            fprintf(stderr,"unable to open CA offset file: %s\n", jbmFECoffsetFileName);
            G192_Reader_Close(&g192);
            return -1;
        }
    }

    /* initialize receiver (wraps decoder) */
    rxerr = EVS_RX_Open(&hRX, st, jbmSafetyMargin);
    if(rxerr)
    {
        fprintf(stderr,"unable to open receiver\n");
        G192_Reader_Close(&g192);
        if( f_offset ) fclose (f_offset);
        return -1;
    }
#ifdef SUPPORT_JBM_TRACEFILE
    rxerr = EVS_RX_SetJbmTraceFileName(hRX, jbmTraceFileName);
    if(rxerr)
    {
        fprintf(stderr,"unable to set JBM trace file name: %s\n", jbmTraceFileName);
        G192_Reader_Close(&g192);
        if( f_offset ) fclose (f_offset);
        EVS_RX_Close(&hRX);
        return -1;
    }
#endif


    /* calculate the delay compensation to have the decoded signal aligned with the original input signal */
    /* the number of first output samples will be reduced by this amount */
    dec_delay = NS2SA(st->output_Fs, get_delay(DEC, st->output_Fs));
    zero_pad = dec_delay;


    /* read first packet */
    g192err = G192_ReadVoipFrame_compact(g192, au, &auSize,
                                         &rtpSequenceNumber, &rtpTimeStamp, &nextPacketRcvTime_ms);
    if(g192err != G192_NO_ERROR)
    {
        fprintf(stderr,"failed to read first RTP packet\n");
        G192_Reader_Close(&g192);
        if( f_offset ) fclose (f_offset);
        EVS_RX_Close(&hRX);
        return -1;
    }

    if( quietMode == 0 )
    {
        fprintf( stdout, "\n------ Running the decoder ------\n\n" );
        fprintf( stdout, "Frames processed:       " );
    }
    else
    {
        fprintf( stdout, "\n-- Start the decoder (quiet mode) --\n\n" );
    }

    /* main receiving/decoding loop */
    for( ; ; )
    {
        unsigned int nSamples = 0;
        /* read all packets with a receive time smaller than the system time */
        while( nextPacketRcvTime_ms <= systemTime_ms )
        {
            /* feed the previous read packet into the receiver now */
            rxerr = EVS_RX_FeedFrame(hRX, au, auSize, rtpSequenceNumber, rtpTimeStamp,
                                     nextPacketRcvTime_ms);
            if (rxerr != EVS_RX_NO_ERROR)
            {
                printf("\nerror in feeding access unit: %8x", rxerr);
                G192_Reader_Close(&g192);
                if( f_offset ) fclose (f_offset);
                EVS_RX_Close(&hRX);
                return -1;
            }
            /* read the next packet */
            g192err = G192_ReadVoipFrame_compact(g192, au, &auSize,
                                                 &rtpSequenceNumber, &rtpTimeStamp, &nextPacketRcvTime_ms);
            if(g192err == G192_READ_ERROR)
            {
                /* finished reading */
                nextPacketRcvTime_ms = (unsigned int) -1;
            }
            else if(g192err != G192_NO_ERROR)
            {
                fprintf(stderr,"failed to read RTP packet\n");
                G192_Reader_Close(&g192);
                if( f_offset ) fclose (f_offset);
                EVS_RX_Close(&hRX);
                return -1;
            }
        }

        /* we are finished when all packets have been received and jitter buffer is empty */
        if( nextPacketRcvTime_ms == (unsigned int)(-1) && EVS_RX_IsEmpty(hRX) )
            break;

        /* decode and get samples */
        rxerr = EVS_RX_GetSamples(hRX, &nSamples, pcmBuf, pcmBufSize, systemTime_ms
                                 );


        EVS_RX_Get_FEC_offset(hRX, &optimum_offset, &FEC_hi);

        if ( st->writeFECoffset == 1 && f_offset )
        {
            if ( FEC_hi == 1)
            {
                fprintf( f_offset, "HI " );
            }
            else
            {
                fprintf( f_offset, "LO " );
            }

            if ( optimum_offset == 1 || optimum_offset == 2 )
            {
                optimum_offset =2;
            }
            else if ( optimum_offset == 3 || optimum_offset == 4 )
            {
                optimum_offset = 3;
            }
            else if ( optimum_offset == 5 || optimum_offset == 6 )
            {
                optimum_offset = 5;
            }
            else if ( optimum_offset >= 7)
            {
                optimum_offset = 7;
            }

            fprintf( f_offset, "%d\n", optimum_offset );
        }

        if(rxerr != EVS_RX_NO_ERROR)
        {
            printf("\nerror in getting samples: %8x", rxerr);
            G192_Reader_Close(&g192);
            if( f_offset ) fclose(f_offset);
            EVS_RX_Close(&hRX);
            return -1;
        }

        if ( dec_delay == 0 )
        {
            fwrite( pcmBuf, sizeof(Word16), nSamples, f_synth );
        }
        else
        {
            if ( dec_delay <= (short)nSamples )
            {
                fwrite( pcmBuf + dec_delay, sizeof(Word16), nSamples - dec_delay, f_synth );
                dec_delay = 0;
            }
            else
            {
                dec_delay -= nSamples;
            }
        }

        frame++;
        if( quietMode == 0 )
        {
            fprintf( stdout, "%-8ld\b\b\b\b\b\b\b\b", frame );
        }
        systemTime_ms += 20;
    }

    /* add zeros at the end to have equal length of synthesized signals */
    set_s( pcmBuf, 0, zero_pad );
    fwrite( pcmBuf, sizeof(Word16), zero_pad, f_synth );

    if (quietMode == 0)
    {
        fprintf( stdout, "\n\n" );
        fprintf( stdout, "Decoding finished\n\n" );
    }
    else
    {
        fprintf( stdout, "Decoding of %ld frames finished\n\n", frame );
    }
    fprintf( stdout, "\n\n" );

    /* free memory etc. */
    G192_Reader_Close(&g192);
    EVS_RX_Close(&hRX);
    return 0;
}


