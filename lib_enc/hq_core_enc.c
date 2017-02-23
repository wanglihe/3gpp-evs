/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"

/*--------------------------------------------------------------------------
 * hq_core_enc()
 *
 * HQ core encoder
 *--------------------------------------------------------------------------*/

void hq_core_enc(
    Encoder_State *st,             /* i/o: encoder state structure */
    const float *audio,          /* i  : input audio signal      */
    const short input_frame_orig,/* i  : frame length            */
    const short hq_core_type,    /* i  : HQ core type            */
    const short Voicing_flag     /* i  : HQ core voicing flag    */
)
{
    short i, is_transient, num_bits, extra_unused;
    float wtda_audio[2 * L_FRAME48k];
    float t_audio[L_FRAME48k];
    short inner_frame, input_frame;
    float ener_match;


    set_f( t_audio, 0, L_FRAME48k );
    st->Nb_ACELP_frames = 0;

    /* set input_frame length */
    input_frame = input_frame_orig;

    st->tcx_cfg.tcx_last_overlap_mode = st->tcx_cfg.tcx_curr_overlap_mode;
    st->tcx_cfg.tcx_curr_overlap_mode = ALDO_WINDOW;

    /*--------------------------------------------------------------------------
     * Preprocessing in the first HQ frame after ACELP frame
     * Find the number of bits for PVQ coding
     * Write signalling information
     *--------------------------------------------------------------------------*/

    num_bits = (short)(st->total_brate / 50);
    extra_unused = 0;

    /*--------------------------------------------------------------------------
     * Detect signal transition
     *--------------------------------------------------------------------------*/

    is_transient = detect_transient( audio, st, input_frame, HQ_CORE );

    /*--------------------------------------------------------------------------
     * Windowing and time-domain aliasing
     * DCT transform
     *--------------------------------------------------------------------------*/

    wtda( audio, wtda_audio, NULL, st->tcx_cfg.tcx_last_overlap_mode, st->tcx_cfg.tcx_curr_overlap_mode, input_frame );

    if ( st->last_core == ACELP_CORE || st->last_core == AMR_WB_CORE )
    {
        /* Preprocessing in the first HQ frame after ACELP frame */
        core_switching_hq_prepare_enc( st, &num_bits, input_frame, wtda_audio, audio );

        /* During ACELP->HQ core switching, limit the HQ core bitrate to 48kbps */
        if( num_bits > HQ_48k / 50 )
        {
            extra_unused = num_bits - (short)(HQ_48k / 50);
            num_bits = (short)(HQ_48k / 50);
        }
    }
    /* subtract signalling bits */
    num_bits -= st->nb_bits_tot;
    direct_transform( wtda_audio, t_audio, is_transient, input_frame );

    /* scale coefficients to their nominal level (8kHz) */
    if( input_frame != NORM_MDCT_FACTOR )
    {
        ener_match = (float)sqrt((float)NORM_MDCT_FACTOR/(float)input_frame);

        for( i=0; i<input_frame; i++ )
        {
            t_audio[i] *= ener_match;
        }
    }

    /* limit encoded band-width according to the command-line OR BWD limitation */
    inner_frame = inner_frame_tbl[st->bwidth];

    if( input_frame > inner_frame )
    {
        if( is_transient )
        {
            for( i = 1; i < NUM_TIME_SWITCHING_BLOCKS; i++ )
            {
                mvr2r( t_audio + i*input_frame/NUM_TIME_SWITCHING_BLOCKS, t_audio + i*inner_frame/NUM_TIME_SWITCHING_BLOCKS, inner_frame/NUM_TIME_SWITCHING_BLOCKS );
            }
        }

        set_f( t_audio + inner_frame, 0.0f, input_frame - inner_frame );
    }

    /*--------------------------------------------------------------------------
     * Classify whether to put extra bits for FER mitigation
     *--------------------------------------------------------------------------*/

    if( st->last_core == HQ_CORE && st->core_brate >  MINIMUM_RATE_TO_ENCODE_VOICING_FLAG )
    {
        if ( Voicing_flag > 0 )
        {
            push_indice( st, IND_HQ_VOICING_FLAG, 1, 1 );
            num_bits -= 1;
        }
        else
        {
            push_indice( st, IND_HQ_VOICING_FLAG, 0, 1 );
            num_bits -= 1;
        }
    }

    /*--------------------------------------------------------------------------
     * Transform-domain encoding
     *--------------------------------------------------------------------------*/

    if( hq_core_type == LOW_RATE_HQ_CORE )
    {
        /* HQ low rate encoder */
        hq_lr_enc( st, t_audio, inner_frame, &num_bits, is_transient );
    }
    else
    {
        /* HQ high rate encoder */
        hq_hr_enc( st, t_audio, inner_frame, &num_bits, is_transient );
    }

    /* write all unused bits to the bitstream */
    num_bits += extra_unused;

    while( num_bits >= 16 )
    {
        push_indice( st, IND_UNUSED, 0, 16 );
        num_bits -= 16;
    }

    if( num_bits != 0 )
    {
        push_indice( st, IND_UNUSED, 0, num_bits );
    }




    return;
}

