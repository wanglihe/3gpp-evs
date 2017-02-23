/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"


/*----------------------------------------------------------------------*
 * inov_decode()
 *
 * Decode the algebraic innovation and do pitch sharpening
 *----------------------------------------------------------------------*/

void inov_decode(
    Decoder_State *st,                        /* i/o: decoder state structure   */
    const long  core_brate,                 /* i  : core bitrate                                */
    const short Opt_AMR_WB,                 /* i  : flag indicating AMR-WB IO mode              */
    const short L_frame,                    /* i  : length of the frame                         */
    const short coder_type,                 /* i  : coding type                                 */
    const short sharpFlag,                  /* i  : formant sharpening flag                     */
    const short i_subfr,                    /* i  : subframe index                              */
    const short tc_subfr,                   /* i  : TC subframe index                           */
    const float *p_Aq,                      /* i  : LP filter coefficients                      */
    const float tilt_code,                  /* i  : tilt of the excitation of previous subframe */
    const float pt_pitch,                   /* i  : pointer to current subframe fractional pitch*/
    float *code                      /* o  : algebraic excitation                        */
)
{
    short nBits;
    float g1, g2;

    if ( L_frame == L_FRAME )
    {
        g1 = FORMANT_SHARPENING_G1;
        g2 = FORMANT_SHARPENING_G2;
    }
    else
    {
        g1 = FORMANT_SHARPENING_G1_16k;
        g2 = FORMANT_SHARPENING_G2_16k;
    }

    if ( !Opt_AMR_WB )
    {
        if( L_frame == L_FRAME )
        {
            nBits = FCB_bits_tbl[BIT_ALLOC_IDX(core_brate, coder_type, i_subfr, TC_SUBFR2IDX(tc_subfr))];
        }
        else  /* L_frame == L_FRAME16k */
        {
            nBits = FCB_bits_16kHz_tbl[BIT_ALLOC_IDX_16KHZ(core_brate, coder_type, i_subfr, TC_SUBFR2IDX_16KHZ(tc_subfr))];
        }

        if( nBits == 7 )
        {
            dec_acelp_1t64( st, code );
        }
        else if( nBits == 12 )
        {
            dec_acelp_2t32( st, code );
        }
        else
        {
            dec_acelp_4t64( st, nBits, code, Opt_AMR_WB );
        }
    }
    else
    {
        if ( core_brate == ACELP_6k60 )
        {
            dec_acelp_2t32( st, code );
        }
        else if ( core_brate == ACELP_8k85 )
        {
            dec_acelp_4t64( st, 20, code, Opt_AMR_WB );
        }
        else if ( core_brate == ACELP_12k65)
        {
            dec_acelp_4t64( st, 36, code, Opt_AMR_WB );
        }
        else if ( core_brate == ACELP_14k25)
        {
            dec_acelp_4t64( st, 44, code, Opt_AMR_WB );
        }
        else if ( core_brate == ACELP_15k85)
        {
            dec_acelp_4t64( st, 52, code, Opt_AMR_WB );
        }
        else if ( core_brate == ACELP_18k25)
        {
            dec_acelp_4t64( st, 64, code, Opt_AMR_WB );
        }
        else if ( core_brate == ACELP_19k85)
        {
            dec_acelp_4t64( st, 72, code, Opt_AMR_WB );
        }
        else
        {
            dec_acelp_4t64( st, 88, code, Opt_AMR_WB );
        }
    }

    cb_shape( 1, 1, 0, sharpFlag, 0, g1, g2, p_Aq, code, tilt_code, pt_pitch );

    return;
}
