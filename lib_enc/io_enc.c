/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"
#include "mime.h"


/*-------------------------------------------------------------------*
 * Local functions
 *-------------------------------------------------------------------*/

static void usage_enc(void);
static char *to_upper( char *str );

/*---------------------------------------------------------------------*
 * io_ini_enc()
 *
 * Processing of command line parameters
 *---------------------------------------------------------------------*/

void io_ini_enc(
    const int   argc,                /* i  : command line arguments number             */
    char  *argv[],             /* i  : command line arguments                    */
    FILE  **f_input,           /* o  : input signal file                         */
    FILE  **f_stream,          /* o  : output bitstream file                     */
    FILE  **f_rate,            /* o  : bitrate switching profile (0 if N/A)      */
    FILE  **f_bwidth,          /* o  : bandwidth switching profile (0 if N/A)    */
    FILE  **f_rf,              /* o  : channel aware configuration file          */
    short  *quietMode,         /* o  : limit printouts                           */
    short  *noDelayCmp,        /* o  : turn off delay compensation               */
    Encoder_State *st                  /* o  : state structure                           */
)
{
    short i = 1, j;
    int tmp;
    short max_bwidth_user = -1;
    char max_bwidth_string[4];
    char stmp[FILENAME_MAX];

    short rf_file_attempt = 1;
    char first_char;
    char rf_file_name[300];

    print_disclaimer( stderr );

    /*-----------------------------------------------------------------*
     * Initialization
     *-----------------------------------------------------------------*/

    *f_input = NULL;
    *f_stream = NULL;
    *quietMode = 0;
    *noDelayCmp = 0;

    st->input_Fs = 16000;
    st->total_brate = ACELP_12k65;

    st->Opt_AMR_WB = 0;
    st->Opt_DTX_ON = 0;
    st->Opt_RF_ON = 0;
    st->rf_fec_offset = 0;
    st->rf_fec_indicator = 1;
    st->max_bwidth = SWB;
    st->interval_SID = FIXED_SID_RATE;
    st->var_SID_rate_flag = 1;
    st->Opt_SC_VBR = 0;
    st->last_Opt_SC_VBR =0;
    st->bitstreamformat = G192;

    if ( argc < 5 )
    {
        usage_enc();
    }

    /*-----------------------------------------------------------------*
     * Optional input arguments
     *-----------------------------------------------------------------*/

    while ( i < argc-4 )
    {
        /*-----------------------------------------------------------------*
         * Bandwidth limitation
         *-----------------------------------------------------------------*/

        if ( strcmp( to_upper(argv[i]), "-MAX_BAND") == 0 )
        {
            strncpy( stmp, argv[i+1], sizeof(stmp) );

            if (strcmp( to_upper(stmp), "-NB") == 0 || strcmp( to_upper(stmp), "NB") == 0)
            {
                st->max_bwidth = max_bwidth_user = NB;
            }
            else if (strcmp( to_upper(stmp), "-WB") == 0 || strcmp( to_upper(stmp), "WB") == 0)
            {
                st->max_bwidth = max_bwidth_user = WB;
            }
            else if (strcmp( to_upper(stmp), "-SWB") == 0 || strcmp( to_upper(stmp), "SWB") == 0)
            {
                st->max_bwidth = max_bwidth_user = SWB;
            }
            else if (strcmp( to_upper(stmp), "-FB") == 0 || strcmp( to_upper(stmp), "FB") == 0)
            {
                st->max_bwidth = max_bwidth_user = FB;
            }
            else if ( (*f_bwidth = fopen(argv[i+1], "rb")) == NULL )
            {
                fprintf(stderr, "Error: incorect bandwidth specification or the bandwidth profile file could not be opened: %s\n\n", argv[i+1]);
                usage_enc();
            }

            /* read first bandwidth value from the profile file (just to check validity) */
            if ( *f_bwidth != NULL )
            {
                if ( fscanf( *f_bwidth, "%d %s", &tmp, stmp) == EOF )
                {
                    fprintf (stderr,"Error: cannot read the bandwidth profile file\n\n");
                    usage_enc();
                }

                rewind(*f_bwidth);

                if ( (st->max_bwidth = CONV_BWIDTH( stmp )) == -1 )
                {
                    fprintf (stderr,"Error: incorrect bandwidth specified (only NB, WB, SWB and FB are supported)\n\n");
                    usage_enc();
                }
                else
                {
                    fprintf(stdout, "Bandwidth switching file: %s\n", argv[i+1]);
                }
            }
            else
            {
                fprintf(stdout, "Maximum encoded bandwidth: %s\n", stmp);
            }

            i += 2;
        }

        /*-----------------------------------------------------------------*
         * Quiet mode
         *-----------------------------------------------------------------*/

        else if( strcmp( to_upper(argv[i]), "-Q" ) == 0 )
        {
            i++;
            *quietMode = 1;
        }

        /*-----------------------------------------------------------------*
         * DTX/CNG
         *-----------------------------------------------------------------*/

        else if( strcmp( to_upper(argv[i]), "-DTX" ) == 0 )
        {
            i++;
            if( i < argc-4 )
            {
                if( sscanf(argv[i], "%d", &tmp) <= 0 )
                {
                    tmp = FIXED_SID_RATE;
                }
                else
                {
                    i++;
                }
            }
            else
            {
                tmp = FIXED_SID_RATE;
            }

            st->Opt_DTX_ON = 1;

            if( tmp == 0 )
            {
                st->var_SID_rate_flag = 1;
                st->interval_SID = 0;
            }
            else if (tmp >= 3 && tmp <= 100)
            {
                st->var_SID_rate_flag = 0;
                st->interval_SID = (short)tmp;
            }
            else
            {
                fprintf(stderr, "Error: Incorrect SID update interval specified: %d (supported 3-100)\n\n", tmp);
                usage_enc();
            }
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
         * Activate channel-aware mode
         *-----------------------------------------------------------------*/

        else if ( strcmp( to_upper(argv[i]), "-RF" ) == 0 )
        {
            st->Opt_RF_ON = 1;
            i++;
            if( i < argc-4 )
            {
                rf_file_attempt = 1;
                sscanf(argv[i], "%s", rf_file_name);
                if (strcmp( to_upper(argv[i]), "LO") == 0)
                {
                    st->rf_fec_indicator = 0;
                    i++;
                    rf_file_attempt = 0;
                }
                else if (strcmp( to_upper(argv[i]), "HI") == 0)
                {
                    st->rf_fec_indicator = 1;
                    i++;
                    rf_file_attempt = 0;
                }

                if ( ( sscanf(argv[i], "%d", &tmp) == 1 ) && ( i < argc-4 ) )
                {
                    if ( tmp == 0 )
                    {
                        st->Opt_RF_ON = 0;
                        st->rf_fec_offset = 0;
                        i++;
                        rf_file_attempt = 0;
                    }
                    else if( tmp == 2 || tmp == 3 || tmp == 5 || tmp == 7 )
                    {
                        st->rf_fec_offset = tmp;
                        i += 1;
                        rf_file_attempt = 0;
                    }
                    else
                    {
                        fprintf (stderr,"Error: incorrect FEC offset specified; RF offset can be 2, 3, 5, or 7. \n\n");
                        usage_enc();
                    }
                }
                else
                {
                    st->rf_fec_offset = FEC_OFFSET;
                }
                if ( rf_file_attempt == 1 )
                {
                    first_char = rf_file_name[0];
                    if ( first_char != '-' )
                    {
                        if ( (*f_rf = fopen(rf_file_name, "r")) == NULL )
                        {
                            fprintf(stderr, "Error: Incorrect specification of channal aware configuration or the CA configuration file could not be opened: %s\n\n", argv[i]);
                            usage_enc();
                        }
                        i++;
                    }
                }
            }
            else
            {
                st->rf_fec_indicator = 1;
                st->rf_fec_offset = FEC_OFFSET;
            }
        }

        /*-----------------------------------------------------------------*
        * MIME output file format
        *-----------------------------------------------------------------*/
        else if( strcmp( to_upper(argv[i]), "-MIME" ) == 0 )
        {
            st->bitstreamformat = MIME;
            i++;
        }

        /*-----------------------------------------------------------------*
         * Option not recognized
         *-----------------------------------------------------------------*/

        else
        {
            fprintf(stderr, "Error: option not recognized, %s\n\n", argv[i]);
            usage_enc();
        }

    } /* end of while  */

    /*-----------------------------------------------------------------*
     * Mandatory input arguments
     *-----------------------------------------------------------------*/

    /*-----------------------------------------------------------------*
     * Bitrate
     *-----------------------------------------------------------------*/

    if ( i < argc-2 )
    {
        /* check if profile file has been entered instead of a fixed bitrate */
        if (sscanf(argv[i], "%d", &tmp) != 1)
        {
            if ( (*f_rate = fopen(argv[i], "rb")) == NULL )
            {
                fprintf(stderr, "Error: bitrate profile file %s could not be opened\n\n", argv[i]);
                usage_enc();
            }

            /* read first bitrate value from the profile file (just to check validity)*/
            if ( fread( &st->total_brate, 4, 1, *f_rate ) != 1 && feof(*f_rate) )
            {
                fprintf (stderr,"Error: cannot read the bitrate profile file\n\n");
                usage_enc();
            }

            rewind(*f_rate);

            fprintf(stdout, "Bitrate switching file: %s\n", argv[i]);
        }
        else
        {
            st->total_brate = tmp;

            /* channel-aware mode is supported only at 13.20 */
            if( st->Opt_RF_ON && st->total_brate != ACELP_13k20 )
            {
                st->Opt_RF_ON = 0;
            }
        }

        /* SC-VBR at 5.90 kbps */
        if ( st->total_brate == ACELP_5k90 )
        {
            st->Opt_SC_VBR = 1;
            st->last_Opt_SC_VBR = st->Opt_SC_VBR;
            st->total_brate = ACELP_7k20;
            if ( st->max_bwidth != NB )
            {
                st->max_bwidth = WB;
            }
        }

        /* check if the entered bitrate is supported */
        j = 0;
        while ( j < SIZE_BRATE_TBL && st->total_brate != brate_tbl[j] )
        {
            j++;
        }

        /* AMR-WB IO mode/EVS primary mode determination */
        if ( j >= SIZE_BRATE_TBL )
        {
            switch ( st->total_brate )
            {
            case ACELP_6k60 :
            case ACELP_8k85 :
            case ACELP_12k65 :
            case ACELP_14k25 :
            case ACELP_15k85 :
            case ACELP_18k25 :
            case ACELP_19k85 :
            case ACELP_23k05 :
            case ACELP_23k85 :
                break;
            default :
            {
                fprintf(stderr, "Error: Incorrect bitrate specification: %ld\n\n", st->total_brate );
                usage_enc();
            }
            break;
            }

            st->Opt_AMR_WB = 1;
        }
        else
        {
            st->Opt_AMR_WB = 0;
        }

        i++;
    }
    else
    {
        fprintf(stderr, "Error: no bitrate specified\n\n");
        usage_enc();
    }

    /*-----------------------------------------------------------------*
     * Input sampling frequency
     *-----------------------------------------------------------------*/

    if( i < argc-2 )
    {
        st->input_Fs = (int) atoi( argv[i] ) * 1000;
        if( st->Opt_RF_ON && st->input_Fs == 8000 )
        {
            /* Reset RF parameters if NB input_Fs;
               For FB input @13.2kbps, the codec is shifted down to SWB below. */
            st->Opt_RF_ON = 0;
            st->rf_fec_offset = 0;
        }
        if( (st->input_Fs != 8000) && ( st->input_Fs != 16000) && (st->input_Fs != 32000) && (st->input_Fs != 48000) )
        {
            fprintf(stderr, "Error: %d kHz is not a supported sampling rate\n\n", atoi( argv[i] ) );
            usage_enc();
        }

        i++;
    }
    else
    {
        fprintf(stderr, "Error: no input sampling frequency specified\n\n");
        usage_enc();
    }
    /*-----------------------------------------------------------------*
     * Input file
     *-----------------------------------------------------------------*/

    if( i < argc-1 )
    {
        if ( (*f_input = fopen(argv[i], "rb")) == NULL )
        {
            fprintf(stderr, "Error: input audio file %s could not be opened\n\n", argv[i]);
            usage_enc();
        }

        fprintf(stdout, "Input audio file:       %s\n", argv[i]);
        i++;
    }
    else
    {
        fprintf(stderr, "Error: no input file specified\n\n");
        usage_enc();
    }

    /*-----------------------------------------------------------------*
     * Output bitstream file
     *-----------------------------------------------------------------*/

    if( i < argc )
    {
        if ( (*f_stream = fopen(argv[i], "wb")) == NULL)
        {
            fprintf(stderr, "Error: output bitstream file %s could not be opened\n\n", argv[i]);
            usage_enc();
        }

        fprintf(stdout, "Output bitstream file:  %s\n", argv[i]);
        i++;
        /* If MIME/storage format selected, write the magic number at the beginning of the bitstream file */
        if( st->bitstreamformat == MIME )
        {
            char buf[4];
            fwrite(EVS_MAGIC_NUMBER, sizeof(char), strlen(EVS_MAGIC_NUMBER), *f_stream);
            buf[0] = buf[1] = buf[2] = 0;
            buf[3] = 1;   /* one channel */
            fwrite(&buf, sizeof(char), 4, *f_stream);
            fprintf(stdout, "Output bitstream file format: MIME");

        }
    }
    else
    {
        fprintf(stderr, "Error: no output bitstream file specified\n\n");
        usage_enc();
    }

    fprintf( stdout, "\n" );


    /* Prevent st->max_bwidth from being higher than Fs/2 */
    if ( st->input_Fs == 8000 && st->max_bwidth > NB )
    {
        st->max_bwidth = NB;
    }
    else if ( st->input_Fs == 16000 && st->max_bwidth > WB )
    {
        st->max_bwidth = WB;
    }
    else if ( st->input_Fs == 32000 && st->max_bwidth > SWB )
    {
        st->max_bwidth = SWB;
    }

    if ( st->Opt_AMR_WB )
    {
        st->codec_mode = MODE1;
    }
    else
    {
        switch ( st->total_brate )
        {
        case 5900:
            st->codec_mode = MODE1;
            break;
        case 7200:
            st->codec_mode = MODE1;
            break;
        case 8000:
            st->codec_mode = MODE1;
            break;
        case 9600:
            st->codec_mode = MODE2;
            break;
        case 13200:
            st->codec_mode = MODE1;
            break;
        case 16400:
            st->codec_mode = MODE2;
            break;
        case 24400:
            st->codec_mode = MODE2;
            break;
        case 32000:
            st->codec_mode = MODE1;
            break;
        case 48000:
            st->codec_mode = MODE2;
            break;
        case 64000:
            st->codec_mode = MODE1;
            break;
        case 96000:
            st->codec_mode = MODE2;
            break;
        case 128000:
            st->codec_mode = MODE2;
            break;
        }
    }

    if( st->total_brate == 13200 && st->Opt_RF_ON == 1)
    {
        st->codec_mode = MODE2;
    }

    st->last_codec_mode = st->codec_mode;

    /*-----------------------------------------------------------------*
     * Print info on screen
     *-----------------------------------------------------------------*/
    /*-----------------------------------------------------------------*
     * Print input signal sampling frequency
     *-----------------------------------------------------------------*/

    fprintf(stdout, "Input sampling rate:    %d Hz\n", st->input_Fs);

    /*-----------------------------------------------------------------*
     * Print bitrate
     *-----------------------------------------------------------------*/

    if ( st->Opt_SC_VBR )
    {
        fprintf(stdout, "Average bitrate:        %.2f kbps\n", (float)ACELP_5k90/1000);
    }
    else
    {
        fprintf(stdout, "Bitrate:                %.2f kbps\n", (float)st->total_brate/1000);
    }

    /*-----------------------------------------------------------------*
     * Print CNG update interval, if DTX is activated
     *-----------------------------------------------------------------*/

    if ( st->Opt_DTX_ON )
    {
        if( st->var_SID_rate_flag )
        {
            fprintf(stdout, "DTX:                    ON, variable CNG update interval\n");
        }
        else
        {
            fprintf(stdout, "DTX:                    ON, CNG update interval = %d frames\n", st->interval_SID);
        }
    }

    /*-----------------------------------------------------------------*
     * Print channel-aware mode info
     *-----------------------------------------------------------------*/

    if ( st->Opt_RF_ON )
    {
        if ( *f_rf == NULL )
        {
            fprintf(stdout, "Channel-aware mode:     ON, FEC indicator : %s  FEC offset: %d \n", (st->rf_fec_indicator==0)?"LO":"HI", st->rf_fec_offset);
        }
        else
        {
            fprintf(stdout, "Channel-aware mode:     ON, Channel-aware config file name:  %s  \n", rf_file_name );
        }
    }

    /*-----------------------------------------------------------------*
     * Print potential limitation of audio bandwidth
     *-----------------------------------------------------------------*/
    switch( st->max_bwidth )
    {
    case NB:
        strncpy(max_bwidth_string, "NB\0", sizeof(max_bwidth_string));
        break;
    case WB:
        strncpy(max_bwidth_string, "WB\0", sizeof(max_bwidth_string));
        break;
    case SWB:
        strncpy(max_bwidth_string, "SWB\0", sizeof(max_bwidth_string));
        break;
    case FB:
        strncpy(max_bwidth_string, "FB\0", sizeof(max_bwidth_string));
        break;
    default:
        memset(max_bwidth_string, 0, sizeof(max_bwidth_string));
        break;
    }

    if ( st->Opt_SC_VBR && !st->Opt_DTX_ON)
    {
        fprintf(stderr, "\nError: SC-VBR 5900 bit/s not supported without DTX\n\n");
        exit(-1);
    }

    if ( (max_bwidth_user != -1) && (st->max_bwidth != max_bwidth_user) )
    {
        fprintf(stdout, "\nBandwidth limited to %s.\n", max_bwidth_string);
    }
    if( (max_bwidth_user == -1) && (st->max_bwidth < FB) && (st->input_Fs == 48000) )
    {
        fprintf(stdout, "\nBandwidth limited to %s. To enable FB coding, please use -max_band FB.\n", max_bwidth_string);
    }
    if( (st->max_bwidth == FB) && (st->total_brate < ACELP_16k40) )
    {
        fprintf(stdout, "\nFB coding not supported below %.2f kbps. ", ACELP_16k40 / 1000.f);
        if( st->total_brate < ACELP_9k60 )
        {
            fprintf(stdout, "Switching to WB.\n");
        }
        else
        {
            fprintf(stdout, "Switching to SWB.\n");
        }
    }
    if( (st->max_bwidth == SWB) && (st->total_brate < ACELP_9k60) )
    {
        fprintf(stdout, "\nSWB coding not supported below %.2f kbps. Switching to WB.", ACELP_9k60 / 1000.f);
    }
    /* in case of 8kHz input sampling or "-max_band NB", require the total bitrate to be below 24.40 kbps */
    if ( ((st->max_bwidth == NB) || (st->input_Fs == 8000)) && (st->total_brate > ACELP_24k40) )
    {
        fprintf(stderr, "\nError: Unsupported mode NB %ld bit/s, NB mode supports rates 5900-24400 bit/s\n\n", st->total_brate);
        usage_enc();
    }


    fprintf(stdout, "\n");
    return;
}

/*---------------------------------------------------------------------*
 * read_next_rfparam()
 *
 * Read next channel aware configuration parameters
 *      p: RF FEC indicator HI or LO, and
 *      o: RF FEC offset in number of frame (2, 3, 5, or 7)
 *---------------------------------------------------------------------*/

void read_next_rfparam(short *rf_fec_offset,     /* o: rf offset                         */
                       short *rf_fec_indicator,  /* o: rf FEC indicator                  */
                       FILE *f_rf               /* i: file pointer to read parameters   */
                      )
{
    char rline[10], str[4];
    short tmp;

    /* initialize */
    *rf_fec_offset = 0;
    *rf_fec_indicator = 1;

    if( f_rf != NULL )
    {
        while ( fgets(rline, 10, f_rf) == NULL && feof(f_rf) )
        {
            rewind(f_rf);
        }
    }
    else
    {
        return;
    }

    if ( sscanf (rline,"%s %hd",str,&tmp) != 2)
    {
        fprintf(stderr, "Error in the RF configuration file. There is no proper configuration line.\n");
        exit(-1);
    }

    /* Read RF FEC indicator */
    if ( strcmp( to_upper(str), "HI" ) == 0 )
    {
        *rf_fec_indicator = 1;
    }
    else if ( strcmp( to_upper(str), "LO" ) == 0 )
    {
        *rf_fec_indicator = 0;
    }
    else
    {
        fprintf(stderr, " Incorrect FEC indicator string. Exiting the encoder.\n");
        exit(-1);
    }

    /* Read RF FEC offset */
    if( tmp == 0 || tmp == 2 || tmp == 3 || tmp == 5 || tmp == 7 )
    {
        *rf_fec_offset = tmp;
    }
    else
    {
        fprintf (stderr,"Error: incorrect FEC offset specified in the RF configration file; RF offset can be 2, 3, 5, or 7. \n");
        exit(-1);
    }
}


/*---------------------------------------------------------------------*
 * read_next_brate()
 *
 * Read next bitrate from profile file (only if invoked on cmd line)
 *---------------------------------------------------------------------*/

void read_next_brate(
    long  *total_brate,           /* i/o: total bitrate                             */
    const int last_total_brate,   /* i  : last total bitrate                        */
    FILE  *f_rate,                /* i  : bitrate switching profile (0 if N/A)      */
    int   input_Fs,               /* i  : input sampling frequency                  */
    short *Opt_AMR_WB,            /* i  : flag indicating AMR-WB IO mode            */
    short *Opt_SC_VBR,            /* i/o: SC-VBR flag                               */
    short *codec_mode             /* i/o: Mode 1 or 2                               */
)
{
    short j;

    /* read next bitrate value from the profile file */
    if( f_rate != NULL )
    {
        while ( fread( total_brate, 4, 1, f_rate ) != 1 && feof(f_rate) )
        {
            rewind(f_rate);
        }
    }
    else
    {
        *total_brate = last_total_brate;
    }

    /* SC-VBR at 5.90 kbps */
    if ( *total_brate == ACELP_5k90 )
    {
        *Opt_SC_VBR = 1;
        *total_brate = ACELP_7k20;
    }
    else
    {
        *Opt_SC_VBR = 0;
    }

    /* check if the entered bitrate is supported */
    j = 0;
    while ( j < SIZE_BRATE_TBL && *total_brate != brate_tbl[j] )
    {
        j++;
    }

    /* AMR-WB IO mode/EVS primary mode determination */
    if ( j >= SIZE_BRATE_TBL )
    {
        switch (*total_brate)
        {
        case ACELP_6k60 :
        case ACELP_8k85 :
        case ACELP_12k65 :
        case ACELP_14k25 :
        case ACELP_15k85 :
        case ACELP_18k25 :
        case ACELP_19k85 :
        case ACELP_23k05 :
        case ACELP_23k85 :
            break;
        default :
        {
            fprintf(stderr, "Error: Incorrect bitrate specification: %ld\n\n", *total_brate );
            usage_enc();
        }
        break;
        }

        *Opt_AMR_WB = 1;
    }
    else
    {
        *Opt_AMR_WB = 0;
    }

    /* in case of 8kHz signal, limit the total bitrate to 24.40 kbps */
    if ( input_Fs == 8000 && *total_brate > ACELP_24k40 )
    {
        *total_brate = ACELP_24k40;
    }

    if ( *Opt_AMR_WB )
    {
        *codec_mode = MODE1;
    }
    else
    {
        switch ( *total_brate )
        {
        case 2800:
            *codec_mode = MODE1;
            break;
        case 5900:
            *codec_mode = MODE1;
            break;
        case 7200:
            *codec_mode = MODE1;
            break;
        case 8000:
            *codec_mode = MODE1;
            break;
        case 9600:
            *codec_mode = MODE2;
            break;
        case 13200:
            *codec_mode = MODE1;
            break;
        case 16400:
            *codec_mode = MODE2;
            break;
        case 24400:
            *codec_mode = MODE2;
            break;
        case 32000:
            *codec_mode = MODE1;
            break;
        case 48000:
            *codec_mode = MODE2;
            break;
        case 64000:
            *codec_mode = MODE1;
            break;
        case 96000:
            *codec_mode = MODE2;
            break;
        case 128000:
            *codec_mode = MODE2;
            break;
        }
    }


    return;
}


/*---------------------------------------------------------------------*
 * read_next_bwidth()
 *
 * Read next bandwidth from profile file (only if invoked on cmd line)
 *---------------------------------------------------------------------*/

void read_next_bwidth(
    short *max_bwidth,            /* i/o: maximum encoded bandwidth                 */
    FILE  *f_bwidth,              /* i  : bandwidth switching profile (0 if N/A)    */
    long  *bwidth_profile_cnt,    /* i/o: counter of frames for bandwidth switching profile file */
    int    input_Fs               /* i  : input sampling rate                       */
)
{
    int res;
    char stmp[4];

    if ( *bwidth_profile_cnt == 0 )
    {
        /* read next bandwidth value and number of frames from the profile file */
        while ( (res = fscanf( f_bwidth, "%ld %3s", bwidth_profile_cnt, stmp)) != 2 && feof(f_bwidth))
        {
            rewind(f_bwidth);
        }

        (*bwidth_profile_cnt)--;

        if (strcmp( to_upper(stmp), "NB") == 0 )
        {
            *max_bwidth = NB;
        }
        else if (strcmp( to_upper(stmp), "WB") == 0 )
        {
            *max_bwidth = WB;
        }
        else if (strcmp( to_upper(stmp), "SWB") == 0 )
        {
            *max_bwidth = SWB;
        }
        else if (strcmp( to_upper(stmp), "FB") == 0 )
        {
            *max_bwidth = FB;
        }

        /* Prevent st->max_bwidth from being higher than Fs/2 */
        if ( input_Fs == 8000 && *max_bwidth > NB )
        {
            *max_bwidth = NB;
        }
        else if ( input_Fs == 16000 && *max_bwidth > WB )
        {
            *max_bwidth = WB;
        }
        else if ( input_Fs == 32000 && *max_bwidth > SWB )
        {
            *max_bwidth = SWB;
        }


    }
    else
    {
        /* current profile still active, only decrease the counter */
        (*bwidth_profile_cnt)--;
    }

    return;
}



/*---------------------------------------------------------------------*
 * to_upper()
 *
 * Capitalize all letters of a string.
 * (normally, _strupr() function would be used but it does not work in Unix)
 *---------------------------------------------------------------------*/

static char *to_upper( char *str )
{
    short i;

    i = 0;
    while (str[i] != 0)
    {
        if (str[i] >= 'a' && str[i] <= 'z') str[i] += 'A' - 'a';
        i++;
    }

    return str;
}

static void usage_enc( void )
{
    fprintf(stdout, "Usage: EVS_cod.exe [Options] R Fs input_file bitstream_file\n\n");

    fprintf(stdout, "Mandatory parameters:\n");
    fprintf(stdout, "---------------------\n");
    fprintf(stdout, "R                : Bitrate in bps, \n");
    fprintf(stdout, "                   for EVS native modes R = (5900*, 7200, 8000, 9600, 13200, 16400,\n");
    fprintf(stdout, "                                             24400, 32000, 48000, 64000, 96000, 128000) \n");
    fprintf(stdout, "                                             *VBR mode (average bitrate),\n");
    fprintf(stdout, "                   for AMR-WB IO modes R =  (6600, 8850, 12650, 14250, 15850, 18250,\n");
    fprintf(stdout, "                                             19850, 23050, 23850) \n");
    fprintf(stdout, "                   Alternatively, R can be a bitrate switching file which consits of R values\n");
    fprintf(stdout, "                   indicating the bitrate for each frame in bit/s. These values are stored in\n");
    fprintf(stdout, "                   binary format using 4 bytes per value\n");
    fprintf(stdout, "Fs               : Input sampling rate in kHz, Fs = (8, 16, 32 or 48) \n");

    fprintf(stdout, "input_file       : Input signal filename \n");
    fprintf(stdout, "bitstream_file   : Output bitstream filename \n\n");

    fprintf(stdout, "Options:\n");
    fprintf(stdout, "--------\n");
    fprintf(stdout, "-q               : Quiet mode, no frame counters\n");
    fprintf(stdout, "                 : default is deactivated\n");
    fprintf(stdout, "-dtx D           : Activate DTX mode, D = (0, 3-100) is the SID update rate\n");
    fprintf(stdout, "                   where 0 = adaptive, 3-100 = fixed in number of frames,\n");
    fprintf(stdout, "                   default is deactivated\n");
    fprintf(stdout, "-dtx             : Activate DTX mode with a SID update rate of 8 frames\n");
    fprintf(stdout, "-rf p o          : Activate channel-aware mode for WB and SWB signal at 13.2kbps, \n");
    fprintf(stdout, "                   where FEC indicator, p: LO or HI, and FEC offset, o: 2, 3, 5, or 7 in number of frames.\n");
    fprintf(stdout, "                   Alternatively p and o can be replaced by a rf configuration file with each line  \n");
    fprintf(stdout, "                   contains the values of p and o separated by a space, \n");
    fprintf(stdout, "                   default is deactivated \n");
    fprintf(stdout, "-max_band B      : Activate bandwidth limitation, B = (NB, WB, SWB or FB)\n");
    fprintf(stdout, "                   alternatively, B can be a text file where each line contains\n");
    fprintf(stdout, "                   \"nb_frames B\"\n");
    fprintf(stdout, "-no_delay_cmp    : Turn off delay compensation\n");

    fprintf(stdout, "-mime            : Mime output bitstream file format\n");
    fprintf(stdout, "                   The encoder produces TS26.445 Annex.2.6 Mime Storage Format, (not RFC4867 Mime Format).\n");
    fprintf(stdout, "                   default output bitstream file format is G.192\n");
    fprintf(stdout, "\n");

    exit(-1);
}
