/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "mime.h"


/*---------------------------------------------------------------------*
 * Local functions
 *---------------------------------------------------------------------*/

static void usage_dec(void);
static char *to_upper( char *str );


/*---------------------------------------------------------------------*
 * io_ini_dec()
 *
 * Processing of command line parameters
 *---------------------------------------------------------------------*/

void io_ini_dec(
    const int argc,              /* i  : command line arguments number             */
    char *argv[],           /* i  : command line arguments                    */
    FILE **f_stream,        /* o  : input bitstream file                      */
    FILE **f_synth,         /* o  : output synthesis file                     */
    short *quietMode,          /* o  : limited printouts                         */
    short *noDelayCmp,         /* o  : turn off delay compensation               */
    Decoder_State *st,               /* o  : Decoder static variables structure        */
#ifdef SUPPORT_JBM_TRACEFILE
    char **jbmTraceFileName,  /* o  : VOIP tracefilename                        */
#endif
    char **jbmFECoffsetFileName /* : Output file  for Optimum FEC offset       */
)
{
    short i;
    Word16 evs_magic, amrwb_magic;

    print_disclaimer( stderr );

    st->writeFECoffset = 0;

    /*-----------------------------------------------------------------*
     * Initialization
     *-----------------------------------------------------------------*/

    i = 1;
    *f_synth = NULL;
    *f_stream = NULL;
    *quietMode = 0;
    *noDelayCmp = 0;

    st->codec_mode = 0; /* unknown before first frame */
    st->Opt_AMR_WB = 0;
    st->Opt_VOIP = 0;

    st->bitstreamformat = G192;
    st->amrwb_rfc4867_flag = -1;

    if ( argc <= 1 )
    {
        usage_dec();
    }

    /*-----------------------------------------------------------------*
     * Optional input arguments
     *-----------------------------------------------------------------*/

    while ( i < argc-3 )
    {
        /*-----------------------------------------------------------------*
         * VOIP mode
         *-----------------------------------------------------------------*/

        if ( strcmp( to_upper(argv[i]), "-VOIP") == 0)
        {
            st->Opt_VOIP = 1;
            i += 1;
        }

#ifdef SUPPORT_JBM_TRACEFILE
        /*-----------------------------------------------------------------*
         * VOIP tracefile
         *-----------------------------------------------------------------*/

        else if ( strcmp( to_upper(argv[i]), "-TRACEFILE" ) == 0 )
        {
            *jbmTraceFileName = argv[i+1];
            i = i + 2;
        }
#endif
        /*-----------------------------------------------------------------*
        * FEC offset file
        *-----------------------------------------------------------------*/

        else if ( strcmp( to_upper(argv[i]), "-FEC_CFG_FILE" ) == 0 )
        {
            st->writeFECoffset = 1;
            *jbmFECoffsetFileName = argv[i+1];
            i = i + 2;
        }

        /*-----------------------------------------------------------------*
         * Quiet mode
         *-----------------------------------------------------------------*/

        else if ( strcmp( to_upper(argv[i]), "-Q" ) == 0 )
        {
            *quietMode = 1;
            i++;
        }

        /*-----------------------------------------------------------------*
         * deactivate delay compensation
         *-----------------------------------------------------------------*/

        else if( strcmp( to_upper(argv[i]), "-NO_DELAY_CMP" ) == 0 )
        {
            *noDelayCmp = 1;
            i++;
        }

        /*-----------------------------------------------------------------*
        * MIME input file format
        *-----------------------------------------------------------------*/
        else if( strcmp( to_upper(argv[i]), "-MIME" ) == 0 )
        {
            st->bitstreamformat = MIME;
            st->amrwb_rfc4867_flag = 0;
            i++;
        }

        /*-----------------------------------------------------------------*
         * Option not recognized
         *-----------------------------------------------------------------*/

        else
        {
            {
                fprintf(stderr, "Error: Unknown option %s\n\n", argv[i]);
                usage_dec();
            }
        }

    } /* end of while  */


    /*-----------------------------------------------------------------*
     * Mandatory input arguments
     *-----------------------------------------------------------------*/

    /*-----------------------------------------------------------------*
     * Output sampling frequency
     *-----------------------------------------------------------------*/

    if( i < argc - 2 )
    {
        st->output_Fs = (int)atoi( argv[i] ) * 1000;

        if( st->output_Fs != 8000 && st->output_Fs != 16000 && st->output_Fs != 32000 && st->output_Fs != 48000 )
        {
            fprintf(stderr, "Error: %d kHz is not a supported sampling rate\n\n", atoi( argv[i] ) );
            usage_dec();
        }

        i++;
    }
    else
    {
        fprintf (stderr, "Error: Sampling rate is not specified\n\n");
        usage_dec();
    }

    /*-----------------------------------------------------------------*
     * Input bitstream file
     *-----------------------------------------------------------------*/

    if( i < argc - 1 )
    {
        if ( (*f_stream = fopen(argv[i], "rb")) == NULL)
        {
            fprintf(stderr,"Error: input bitstream file %s cannot be opened\n\n", argv[i]);
            usage_dec();
        }
        /* If MIME/storage format selected, scan for the magic number at the beginning of the bitstream file */
        if( st->bitstreamformat == MIME )
        {
            char buf[13];
            evs_magic   = 1 ;
            amrwb_magic = 1;

            if(NULL == fgets(buf, 13, *f_stream))
            {
                fprintf(stderr,"Error: input bitstream file %s cannot be read\n\n", argv[i]);
                usage_dec();
            }
            /* verify AMRWB magic number */
            if ( strncmp(buf, AMRWB_MAGIC_NUMBER, strlen(AMRWB_MAGIC_NUMBER)))
            {
                amrwb_magic = 0;
            }
            if ( strncmp(buf, EVS_MAGIC_NUMBER, strlen(EVS_MAGIC_NUMBER)))
            {
                evs_magic = 0;
            }
            if( evs_magic != 0 )
            {
                if ((fread(&buf,sizeof(char), 4, *f_stream) != 4 ) || !((buf[3] == 1) && (buf[2] == 0) && (buf[1] == 0) &&  (buf[0] == 0)) )
                {
                    fprintf(stderr, "Error: input bitstream file %s specifies unsupported number of evs audio channels\n\n",argv[i]);
                    usage_dec();
                }
            }

            if( evs_magic == 0 &&  amrwb_magic == 0 )
            {
                /* no valid MIME magic number  */
                fprintf(stderr, "Error: input bitstream file %s specifies unsupported MIME magic number (%13s) \n\n",argv[i],buf );
                usage_dec();
            }

            if( evs_magic )
            {
                fprintf( stderr, "Found MIME Magic number %s\n", EVS_MAGIC_NUMBER );
                st->amrwb_rfc4867_flag = 0;
            }
            else
            {
                fprintf( stderr, "Found MIME Magic number %s\n",AMRWB_MAGIC_NUMBER );
                st->amrwb_rfc4867_flag = 1;
                st->Opt_AMR_WB = 1;    /*  needed in case first initial RFC4867 frames/ToCs are lost */
            }
        }

        else if( st->Opt_VOIP == 0 )
        {
            /* G.192 format ....  preread the G.192 sync header */
            unsigned short  utmp;
            if ( fread( &utmp, sizeof(unsigned short), 1, *f_stream ) != 1 )
            {
                /* error during pre-reading */
                if( ferror( *f_stream ) )
                {
                    fprintf(stderr, "Error: input G.192 bitstream file %s , can not be read  \n\n",argv[i] );
                }
                else
                {
                    fprintf(stderr, "Error: input G.192 bitstream file %s , has zero size, can not be read  \n\n",argv[i] );
                }
                usage_dec();
            }
            if( utmp != SYNC_GOOD_FRAME && utmp != SYNC_BAD_FRAME )
            {
                /* check for a valid first G.192 synch  word in Sync Header  */
                fprintf(stderr, "Error: input bitstream file %s does not have a valid G.192 synch word value \n\n",argv[i]);
                usage_dec();
            }
            /* now rewind the G.192 bitstream file */
            fseek( *f_stream , 0L, SEEK_SET );
        }
        /*  JBM format */

        fprintf( stdout, "Input bitstream file:   %s\n", argv[i]);

        i++;
    }
    else
    {
        fprintf (stderr, "Error: no input bitstream file specified\n\n");
        usage_dec();
    }

    /*-----------------------------------------------------------------*
     * Output synthesis file
     *-----------------------------------------------------------------*/

    if( i < argc )
    {
        if ( (*f_synth = fopen(argv[i], "wb")) == NULL )
        {
            fprintf( stderr, "Error: ouput synthesis file %s cannot be opened\n\n", argv[i] );
            usage_dec();
        }

        fprintf( stdout, "Output synthesis file:  %s\n", argv[i] );
        i++;
    }
    else
    {
        fprintf( stderr, "Error: no output synthesis file specified\n\n" );
        usage_dec();
    }

    fprintf( stdout, "\n" );

    if( !st->Opt_VOIP )
    {
        /*-----------------------------------------------------------------*
         * Read information from bitstream
         *-----------------------------------------------------------------*/
        st->ini_frame = 0; /* initialize, since this is needed within read_indices, to correctly set st->last_codec_mode */
        if( st->bitstreamformat == G192 )
        {
            read_indices( st, *f_stream, 1 );
        }
        else
        {
            st->total_brate=0;  /* make sure total_brate is deterministic  even if there are no MIME ToCs */
            read_indices_mime( st, *f_stream, 1 );   /* rew_flag == 1 ,  checks only very first  frame    */
            if( st->amrwb_rfc4867_flag != 0 )
            {
                fseek(*f_stream,strlen(AMRWB_MAGIC_NUMBER), SEEK_SET);    /* restart after 9 bytes */
            }
            else
            {
                fseek(*f_stream,strlen(EVS_MAGIC_NUMBER)+4, SEEK_SET); /* restart after  16 bytes */
            }
        }

        /*-----------------------------------------------------------------*
         * Print info on screen
         *-----------------------------------------------------------------*/

        /*-----------------------------------------------------------------*
         * Print output sampling frequency
         *-----------------------------------------------------------------*/

        fprintf( stdout, "Output sampling rate:   %d Hz\n", st->output_Fs );

        /*-----------------------------------------------------------------*
         * Print bitrate
         *-----------------------------------------------------------------*/

        fprintf( stdout, "Bitrate:                %.2f kbps\n", (float)st->total_brate/1000 );

        if( st->total_brate <= 0 )
        {
            if( st->bitstreamformat == G192 )
            {
                fprintf( stdout, "Active Bitrate not identified in bitstream file \n" );
            }
            else /* MIME */
            {
                fprintf( stdout, "Active Bitrate not identified from first MIME frame \n" );
            }
        }

    }
    return;
}


/*---------------------------------------------------------------------*
 * to_upper()
 *
 * Capitalize all letters of a string.
 * (normally to_upper() function would be used but it does not work in Unix)
 *---------------------------------------------------------------------*/

static char *to_upper( char *str )
{
    short i;
    char *p = str;

    i = 0;
    while (str[i] != 0)
    {
        if (str[i] >= 'a' && str[i] <= 'z') str[i] -= 0x20;
        i++;
    }

    return p;
}

static void usage_dec( void )
{
    fprintf(stdout,"Usage : EVS_dec.exe [Options] Fs bitstream_file output_file\n\n");

    fprintf(stdout,"Mandatory parameters:\n");
    fprintf(stdout,"---------------------\n");
    fprintf(stdout,"Fs                  : Output sampling rate in kHz (8, 16, 32 or 48)\n");
    fprintf(stdout,"bitstream_file      : Input bitstream filename or RTP packet filename (in VOIP mode)\n");
    fprintf(stdout,"output_file         : Output speech filename \n\n");

    fprintf(stdout,"Options:\n");
    fprintf(stdout,"--------\n");

    fprintf(stdout, "-VOIP              : VOIP mode\n");
#ifdef SUPPORT_JBM_TRACEFILE
    fprintf(stdout, "-Tracefile TF      : Generate trace file named TF\n");
#endif
    fprintf(stdout, "-no_delay_cmp      : Turn off delay compensation\n");
    fprintf(stdout, "-fec_cfg_file      : Optimal channel aware configuration computed by the JBM   \n");
    fprintf(stdout, "                     as described in Section 6.3.1 of TS26.448. The output is \n");
    fprintf(stdout, "                     written into a .txt file. Each line contains the FER indicator \n");
    fprintf(stdout, "                     (HI|LO) and optimal FEC offset. \n");

    fprintf(stdout, "-mime              : Mime bitstream file format\n");
    fprintf(stdout, "                     The decoder may read both TS26.445 Annex.2.6 and RFC4867 Mime Storage\n");
    fprintf(stdout, "                     Format files, the magic word in the mime file is used to determine\n");
    fprintf(stdout, "                     which of the two supported formats is in use.\n");
    fprintf(stdout, "                     default bitstream file format is G.192\n");
    fprintf(stdout, "-q                 : Quiet mode, no frame counter\n");
    fprintf(stdout, "                     default is OFF\n");
    fprintf(stdout, "\n");

    exit(-1);
}

