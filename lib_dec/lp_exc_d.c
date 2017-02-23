/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "prot.h"

/*---------------------------------------------------------------------*
 * lp_filt_exc_dec()
 *
 * Low-pass filtering of the adaptive exctitation
 *---------------------------------------------------------------------*/

void lp_filt_exc_dec(
    Decoder_State *st,             /* i/o: decoder state structure                 */
    const short codec_mode,     /* i : coder mode                               */
    const long  core_brate,     /* i  : core bitrate                            */
    const short Opt_AMR_WB,     /* i  : flag indicating AMR-WB IO mode          */
    const short coder_type,     /* i  : coding type                             */
    const short i_subfr,        /* i  : subframe index                          */
    const short L_subfr,        /* i  : subframe size                           */
    const short L_frame,        /* i  : frame size                              */
    short lp_flag,        /* i  : operation mode signalling               */
    float *exc            /* i/o: pointer to the excitation signal frame  */
)
{
    short i;
    float code[L_FRAME];

    /*-----------------------------------------------------------------*
     * Select LP filtering of the adaptive excitation
     *-----------------------------------------------------------------*/

    if( codec_mode == MODE1 )
    {
        if ( ( Opt_AMR_WB || coder_type == GENERIC || coder_type == TRANSITION ) && core_brate < ACELP_11k60 )
        {
            lp_flag = LOW_PASS;
        }
        else if ( core_brate >= ACELP_11k60 )
        {
            lp_flag = (short)get_next_indice( st, 1 );
        }
        else
        {
            lp_flag = FULL_BAND;
        }
    }

    /*--------------------------------------------------------------------*
     * Find pitch excitation with LP filter
     *--------------------------------------------------------------------*/

    if ( lp_flag == LOW_PASS )
    {
        /* pointer positionning to avoid doing it inside the loop */
        if( codec_mode == MODE2 && L_frame == L_FRAME16k )
        {
            for (i=0; i<L_subfr; i++)
            {
                code[i] = (float)(0.21f * exc[i-1+i_subfr] + 0.58f * exc[i+i_subfr] + 0.21f * exc[i+1+i_subfr]);
            }
        }
        else
        {
            for (i=0; i<L_subfr; i++)
            {
                code[i] = (float)(0.18f * exc[i-1+i_subfr] + 0.64f * exc[i+i_subfr] + 0.18f * exc[i+1+i_subfr]);
            }
        }

        mvr2r(code, &exc[i_subfr], L_subfr);
    }

    return;
}
