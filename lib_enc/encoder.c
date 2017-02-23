/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <stdio.h>
#include <wchar.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include <sys/stat.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"
#include "g192.h"

/*------------------------------------------------------------------------------------------*
 * Global variables
 *------------------------------------------------------------------------------------------*/
long frame = 0;                 /* Counter of frames */



/*------------------------------------------------------------------------------------------*
 * Main encoder function
 *------------------------------------------------------------------------------------------*/

int main( int argc, char** argv )
{

    FILE *f_input;                                        /* MODE1 - input signal file */
    FILE *f_stream;                                       /* MODE1 - output bitstream file */
    FILE *f_rate = NULL;                                  /* MODE1 - bitrate switching profile file */
    FILE *f_bwidth = NULL;                                /* MODE1 - bandwidth switching profile file */
    FILE *f_rf = NULL;                                    /* Channel aware configuration file */
    long bwidth_profile_cnt = 0;                          /* MODE1 - counter of frames for bandwidth switching profile file */
    short tmps, input_frame, enc_delay, n_samples;
    short quietMode = 0;
    short noDelayCmp = 0;
    short data[L_FRAME48k];                               /* 'short' buffer for input signal */
    Indice ind_list[MAX_NUM_INDICES];                     /* list of indices */
    Encoder_State *st;                                    /* MODE1 - encoder state structure */
    UWord8 pFrame[(MAX_BITS_PER_FRAME + 7) >> 3];
    Word16 pFrame_size = 0;
    short Opt_RF_ON_loc, rf_fec_offset_loc;






    /*------------------------------------------------------------------------------------------*
     * Allocate memory for static variables
     * Processing of command-line parameters
     * Encoder initialization
     *------------------------------------------------------------------------------------------*/

    if ( (st = (Encoder_State *) malloc( sizeof(Encoder_State) ) ) == NULL )
    {
        fprintf(stderr, "Can not allocate memory for encoder state structure\n");
        exit(-1);
    }

    io_ini_enc( argc, argv, &f_input, &f_stream, &f_rate, &f_bwidth,
                &f_rf, &quietMode, &noDelayCmp, st );

    Opt_RF_ON_loc = st->Opt_RF_ON;
    rf_fec_offset_loc = st->rf_fec_offset;

    st->ind_list = ind_list;
    init_encoder( st );

    input_frame = (short)(st->input_Fs / 50);

    /*------------------------------------------------------------------------------------------*
     * Compensate for encoder delay (bitstream aligned with input signal)
     * Compensate for the rest of codec delay (local synthesis aligned with decoded signal and original signal)
     *------------------------------------------------------------------------------------------*/
    enc_delay = NS2SA( st->input_Fs, get_delay(ENC, st->input_Fs) );
    if ( noDelayCmp == 0 )
    {
        /* read samples and throw them away */
        if( (tmps = (short)fread(data, sizeof(short), enc_delay, f_input)) != enc_delay )
        {
            fprintf(stderr, "Can not read the data from input file\n");
            exit(-1);
        }
        /* the number of first output samples will be reduced by this amount */
    }


    /*------------------------------------------------------------------------------------------*
     * Loop for every frame of input data
     * - Read the input data
     * - Process switching files
     * - Run the encoder
     * - Write the parameters into output bitstream file
     *------------------------------------------------------------------------------------------*/

    if( quietMode == 0 )
    {
        fprintf( stdout, "\n------ Running the encoder ------\n\n" );
        fprintf( stdout, "Frames processed:       " );
    }
    else
    {
        fprintf( stdout, "\n-- Start the encoder (quiet mode) --\n\n" );
    }

    while( (n_samples = (short)fread(data, sizeof(short), input_frame, f_input)) > 0 )
    {
        if(f_rf != NULL)
        {
            read_next_rfparam( &st->rf_fec_offset, &st->rf_fec_indicator, f_rf);
            rf_fec_offset_loc = st->rf_fec_offset;
        }

        if (f_rate != NULL)
        {
            /* read next bitrate from profile file (only if invoked on the cmd line) */
            read_next_brate( &st->total_brate, st->last_total_brate,
                             f_rate, st->input_Fs, &st->Opt_AMR_WB, &st->Opt_SC_VBR, &st->codec_mode);
        }

        if (f_bwidth != NULL)
        {
            /* read next bandwidth from profile file (only if invoked on the cmd line) */
            read_next_bwidth( &st->max_bwidth, f_bwidth, &bwidth_profile_cnt, st->input_Fs );
        }

        if( ( st->Opt_RF_ON && ( st->total_brate != ACELP_13k20 ||  st->input_Fs == 8000 || st->max_bwidth == NB ) ) || st->rf_fec_offset == 0 )
        {
            if( st->total_brate == ACELP_13k20 )
            {
                st->codec_mode = MODE1;
                reset_rf_indices(st);
            }
            st->Opt_RF_ON = 0;
            st->rf_fec_offset = 0;
        }

        if( Opt_RF_ON_loc && rf_fec_offset_loc != 0 && L_sub( st->total_brate, ACELP_13k20 ) == 0 && L_sub( st->input_Fs, 8000 ) != 0 && st->max_bwidth != NB )
        {
            st->codec_mode = MODE2;
            if(st->Opt_RF_ON == 0)
            {
                reset_rf_indices(st);
            }
            st->Opt_RF_ON = 1;
            st->rf_fec_offset = rf_fec_offset_loc;
        }

        /* in case of 8kHz sampling rate or when in "max_band NB" mode, limit the total bitrate to 24.40 kbps */
        if ( ((st->input_Fs == 8000)|| (st->max_bwidth == NB)) && (st->total_brate > ACELP_24k40) )
        {
            st->total_brate = ACELP_24k40;
            st->codec_mode = MODE2;
        }


        /* run the main encoding routine */

        if ( st->Opt_AMR_WB )
        {
            amr_wb_enc( st, data, n_samples );
        }
        else
        {
            evs_enc( st, data, n_samples );
        }

        /* pack indices into serialized payload format */
        if( st->bitstreamformat == MIME )
        {
            indices_to_serial( st, pFrame, &pFrame_size );
        }

        /* write indices into bitstream file */
        write_indices( st, f_stream
                       , pFrame, pFrame_size
                     );

        fflush( stderr );

        frame++;
        if( quietMode == 0 )
        {
            fprintf( stdout, "%-8ld\b\b\b\b\b\b\b\b", frame );
        }
    }
    if( quietMode == 0 )
    {
        fprintf( stdout, "\n\n" );
        fprintf( stdout, "Encoding finished\n\n" );
    }
    else
    {
        fprintf( stdout, "Encoding of %ld frames finished\n\n", frame );
    }





    /*------------------------------------------------------------------------------------------*
     * Close files and free ressources
     *------------------------------------------------------------------------------------------*/

    fclose( f_input );
    fclose( f_stream );
    if ( f_rate ) fclose ( f_rate );
    if ( f_bwidth ) fclose ( f_bwidth );
    destroy_encoder( st );
    free( st );


    return 0;
}


