/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"


/*---------------------------------------------------------------------*
 * inov_encode()
 *
 * Encode the algebraic innovation
 *---------------------------------------------------------------------*/

void inov_encode(
    Encoder_State *st,              /* i/o: encoder state structure                         */
    const long  core_brate,       /* i  : core bitrate                                    */
    const short Opt_AMR_WB,       /* i  : flag indicating AMR-WB IO mode                  */
    const short L_frame,          /* i  : length of the frame                             */
    const short last_L_frame,     /* i  : length of the last frame                        */
    const short coder_type,       /* i  : coding type                                     */
    const short bwidth,           /* i  : input signal bandwidth                          */
    const short sharpFlag,        /* i  : formant sharpening flag                         */
    const short i_subfr,          /* i  : subframe index                                  */
    const short tc_subfr,         /* i  : TC subframe index                               */
    const float *p_Aq,            /* i  : LP filter coefficients                          */
    const float gain_pit,         /* i  : adaptive excitation gain                        */
    float *cn,              /* i/o: target vector in residual domain                */
    const float *exc,             /* i  : pointer to excitation signal frame              */
    float *h1,              /* i/o: weighted filter input response                  */
    const float tilt_code,        /* i  : tilt of the excitation of previous subframe     */
    const float pt_pitch,         /* i  : pointer to current subframe fractional pitch    */
    const float *xn2,             /* i  : target vector for innovation search             */
    float *code,            /* o  : algebraic excitation                            */
    float *y2,              /* o  : zero-memory filtered algebraic excitation       */
    short *unbits           /* o  : number of unused bits for  EVS_PI               */
)
{
    float dn[L_SUBFR];            /* Correlation between xn2 and h1 */
    short nBits, cmpl_flag;
    short k;
    float g1, g2;
    float cn2[L_SUBFR];
    float Rw[L_SUBFR];
    short i, acelpautoc;

    if( L_frame == L_FRAME )
    {
        g1 = FORMANT_SHARPENING_G1;
        g2 = FORMANT_SHARPENING_G2;
    }
    else
    {
        g1 = FORMANT_SHARPENING_G1_16k;
        g2 = FORMANT_SHARPENING_G2_16k;
    }

    /*----------------------------------------------------------------*
     * Update target vector for codebook search in residual domain
     * Preemphasize the impulse response and include fixed-gain pitch contribution into impulse resp. h1[] (pitch sharpenning)
     * Correlation between target xn2[] and impulse response h1[]
     *----------------------------------------------------------------*/

    if( core_brate > ACELP_13k20 && !Opt_AMR_WB )
    {
        acelpautoc = 1;

        cb_shape( 1, 1, 0, sharpFlag, 0, g1, g2, p_Aq, h1, tilt_code, pt_pitch );

        corr_xh( h1, Rw, h1, L_SUBFR );

        for( k=0; k<L_SUBFR; k++ )
        {
            cn2[k] = xn2[k];

            for( i=0; i<k; i++ )
            {
                cn2[k] -= cn2[i] * h1[k-i];
            }
        }

        E_ACELP_toeplitz_mul( Rw, cn2, dn );
        mvr2r( cn2, cn, L_SUBFR );
    }
    else
    {
        acelpautoc = 0;

        updt_tar( cn, cn, &exc[i_subfr], gain_pit, L_SUBFR );

        cb_shape( 1, 1, 0, sharpFlag, 0, g1, g2, p_Aq, h1, tilt_code, pt_pitch );

        corr_xh( xn2, dn, h1, L_SUBFR );
    }

    /*-----------------------------------------------------------------*
     * Set complexity reduction flag to limit the number of iterations
     * in algebraic innovation search
     *-----------------------------------------------------------------*/

    cmpl_flag = 0;

    if( L_frame == L_FRAME && coder_type == TRANSITION )
    {
        if( core_brate == ACELP_8k00 && i_subfr == 0 && tc_subfr < L_SUBFR )
        {
            cmpl_flag = 3;
        }

        if( core_brate == ACELP_11k60 && ( (i_subfr == 0 && tc_subfr < L_SUBFR) || tc_subfr == TC_0_0 || (i_subfr == 3*L_SUBFR && tc_subfr == TC_0_64)) )
        {
            cmpl_flag = 3;
        }

        if( (core_brate == ACELP_13k20 || core_brate == ACELP_12k15 ) && ( (i_subfr == 0 && tc_subfr < L_SUBFR) || tc_subfr <= TC_0_64 ) )
        {
            cmpl_flag = 3;
        }
    }

    if( L_frame == L_FRAME16k )
    {
        if( core_brate <= ACELP_32k )
        {
            cmpl_flag = 4;

            if( coder_type == TRANSITION && bwidth > WB )
            {
                if( i_subfr <= L_SUBFR )
                {
                    cmpl_flag -= 1;
                }
                else
                {
                    cmpl_flag -= 2;
                }
            }
        }
        else if( core_brate <= ACELP_48k )
        {
            cmpl_flag = 3;

            if( coder_type == TRANSITION )
            {
                if( i_subfr <= L_SUBFR )
                {
                    cmpl_flag -= 1;
                }
                else
                {
                    cmpl_flag -= 2;
                }
            }
        }
        else
        {
            cmpl_flag = 4;

            if( coder_type == TRANSITION )
            {
                if( i_subfr <= L_SUBFR )
                {
                    cmpl_flag -= 1;
                }
                else
                {
                    cmpl_flag -= 2;
                }
            }
        }

        if( coder_type == INACTIVE )
        {
            cmpl_flag = 4;
        }
    }

    /* reduce number of iterations in a frame where there is an internal sampling rate switch in order not to increase the WC complexity  */
    if( L_frame != last_L_frame && core_brate > ACELP_13k20 && (core_brate < ACELP_32k || bwidth == WB) )
    {
        if( cmpl_flag > 1 )
        {
            cmpl_flag--;
        }
    }

    /*-----------------------------------------------------------------*
     * Find and encode the algebraic innovation
     *-----------------------------------------------------------------*/

    set_f( y2, 0, L_SUBFR );

    if( !Opt_AMR_WB )
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
            acelp_1t64( st, dn, h1, code, y2 );
        }
        else if( nBits == 12 )
        {
            acelp_2t32( st, dn, h1, code, y2 );
        }
        else
        {
            *unbits += acelp_4t64( st, dn, cn, h1, Rw, acelpautoc, code, y2, nBits, cmpl_flag, Opt_AMR_WB );
        }
    }
    else
    {
        if (core_brate == ACELP_6k60 )
        {
            acelp_2t32( st, dn, h1, code, y2 );
        }
        else if( (core_brate == ACELP_8k85) )
        {
            acelp_4t64( st, dn, cn, h1, Rw, acelpautoc, code, y2, 20, cmpl_flag, Opt_AMR_WB );
        }
        else if( core_brate == ACELP_12k65 )
        {
            acelp_4t64( st, dn, cn, h1, Rw, acelpautoc, code, y2, 36, cmpl_flag, Opt_AMR_WB );
        }
        else if( core_brate == ACELP_14k25 )
        {
            acelp_4t64( st, dn, cn, h1, Rw, acelpautoc, code, y2, 44, cmpl_flag, Opt_AMR_WB );
        }
        else if( core_brate == ACELP_15k85 )
        {
            acelp_4t64( st, dn, cn, h1, Rw, acelpautoc, code, y2, 52, cmpl_flag, Opt_AMR_WB );
        }
        else if( core_brate == ACELP_18k25 )
        {
            acelp_4t64( st, dn, cn, h1, Rw, acelpautoc, code, y2, 64, cmpl_flag, Opt_AMR_WB );
        }
        else if( core_brate == ACELP_19k85 )
        {
            acelp_4t64( st, dn, cn, h1, Rw, acelpautoc, code, y2, 72, cmpl_flag, Opt_AMR_WB );
        }
        else if( core_brate == ACELP_23k05 )
        {
            acelp_4t64( st, dn, cn, h1, Rw, acelpautoc, code, y2, 88, cmpl_flag, Opt_AMR_WB );
        }
        else if( core_brate == ACELP_23k85 )
        {
            acelp_4t64( st, dn, cn, h1, Rw, acelpautoc, code, y2, 88, cmpl_flag, Opt_AMR_WB );
        }
    }

    /*----------------------------------------------------------------*
     * Pitch sharpening
     *----------------------------------------------------------------*/

    cb_shape( 1, 1, 0, sharpFlag, 0, g1, g2, p_Aq, code, tilt_code, pt_pitch );

    return;
}
