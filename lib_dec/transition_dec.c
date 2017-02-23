/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"

/*----------------------------------------------------------------------*
 * Local functions
 *----------------------------------------------------------------------*/

static void tc_dec( Decoder_State *st, const short L_frame, float exc[], short *T0, short *T0_frac, const short i_subfr,
                    const short tc_subfr, short *position, const long core_brate, float bwe_exc[] );

/*-------------------------------------------------------------------*
 * transition_dec()
 *
 * Principal function for TC decoding
 *-------------------------------------------------------------------*/

void transition_dec(
    Decoder_State *st,          /* i/o: decoder state structure                 */
    const long  core_brate,   /* i  : core bitrate                            */
    const short Opt_AMR_WB,   /* i  : flag indicating AMR-WB IO mode          */
    const short L_frame,      /* i  : length of the frame                     */
    const short i_subfr,      /* i  : subframe index                          */
    const short coder_type,   /* i  : coder type                              */
    const short tc_subfr,     /* i  : TC subframe index                       */
    short *Jopt_flag,   /* i  : joint optimization flag                 */
    float *exc,         /* o  : excitation signal                       */
    short *T0,          /* o  : close loop integer pitch                */
    short *T0_frac,     /* o  : close loop fractional part of the pitch */
    short *T0_min,      /* i/o: delta search min for sf 2 & 4           */
    short *T0_max,      /* i/o: delta search max for sf 2 & 4           */
    float **pt_pitch,   /* o  : floating pitch values                   */
    short *position,    /* i/o: first glottal impulse position in frame */
    float *bwe_exc      /* o  : excitation for SWB TBE                  */
)
{
    short i, pit_flag, pit_start, pit_limit, index, nBits;
    short limit_flag;
    short offset;

    /* Set limit_flag to 0 for restrained limits, and 1 for extended limits */
    limit_flag = 0;

    /*---------------------------------------------------------------------*
     * zero adaptive contribution (glottal shape codebook search not
     *                             in first subframe(s) )
     *---------------------------------------------------------------------*/

    if( tc_subfr > i_subfr+TC_0_192 )
    {
        set_f(&exc[i_subfr], 0, L_SUBFR);

        if( L_frame == L_FRAME )
        {
            set_f(&bwe_exc[i_subfr*HIBND_ACB_L_FAC], 0, (short) (L_SUBFR*HIBND_ACB_L_FAC));         /* set past excitation buffer to 0 */
        }
        else
        {
            set_f(&bwe_exc[i_subfr*2], 0, (short) (L_SUBFR*2));         /* set past excitation buffer to 0 */
        }

        *T0 = L_SUBFR;
        *T0_frac = 0;
        **pt_pitch = (float)L_SUBFR;
    }

    /*---------------------------------------------------------------------*
     * glottal shape codebook search
     *---------------------------------------------------------------------*/

    else if( (tc_subfr-i_subfr >= 0) && (tc_subfr-i_subfr <= TC_0_192) )
    {
        set_f( exc-L_EXC_MEM, 0, L_EXC_MEM );         /* set past excitation buffer to 0 */

        if( L_frame == L_FRAME )
        {
            set_f( bwe_exc-PIT_MAX*HIBND_ACB_L_FAC, 0, PIT_MAX*HIBND_ACB_L_FAC);         /* set past excitation buffer to 0 */
        }
        else
        {
            set_f( bwe_exc-PIT16k_MAX*2, 0, PIT16k_MAX*2);         /* set past excitation buffer to 0 */
        }

        /* glottal shape codebook contribution construction */
        tc_dec( st, L_frame, exc, T0, T0_frac, i_subfr, tc_subfr, position, core_brate, bwe_exc );

        **pt_pitch = (float)(*T0) + (float)(*T0_frac)/4.0f;   /* save subframe pitch values  */
        *Jopt_flag = 1;
    }

    /*---------------------------------------------------------------------*
     * Regular ACELP Decoding using GENERIC type decoder
     * (all subframes following subframe with glottal shape codebook seach)
     * - search the position of the 2nd glottal impulse in case that the first
     *   one is in the 1st subframe (different adaptive contribution
     *   construction and the pitch period coding is used)
     *---------------------------------------------------------------------*/

    else if (tc_subfr < i_subfr)
    {
        if( L_frame == L_FRAME )
        {
            *Jopt_flag = 1;

            if( (i_subfr - tc_subfr >= L_SUBFR) && (i_subfr - tc_subfr <= L_SUBFR + TC_0_192) )
            {
                pit_flag = 0;
            }
            else
            {
                pit_flag = L_SUBFR;
            }

            if( tc_subfr == TC_0_0 )
            {
                if( i_subfr == L_SUBFR )
                {
                    limit_T0( L_FRAME, 8, pit_flag, limit_flag, *T0, 0, T0_min, T0_max );
                }

                pit_flag = 1;
            }

            /*-----------------------------------------------------------------*
             * get number of bits for pitch decoding
             *-----------------------------------------------------------------*/

            nBits = ACB_bits_tbl[BIT_ALLOC_IDX(core_brate, TRANSITION, i_subfr, TC_SUBFR2IDX(tc_subfr))];

            /*------------------------------------------------------------*
             * first glottal impulse is in the 1st subframe
             *------------------------------------------------------------*/

            if( (i_subfr == L_SUBFR) && (tc_subfr >= TC_0_128)  )
            {
                /*--------------------------------------------------------*
                 * second glottal impulse is in the 3rd or 4th subframe
                 * - build exc[] in 2nd subframe
                 *--------------------------------------------------------*/

                *T0 = 2*L_SUBFR;
                *T0_frac = 0;
                *Jopt_flag = 0;

                /* set adaptive part of exciation for curent subframe to 0 */
                set_f( &exc[i_subfr], 0, (short)(L_SUBFR+1) );

                set_f( &bwe_exc[i_subfr*HIBND_ACB_L_FAC], 0, (short)(L_SUBFR*HIBND_ACB_L_FAC) );
            }
            else if( (i_subfr == L_SUBFR) && (tc_subfr == TC_0_64) )
            {
                /*--------------------------------------------------------*
                 * second glottal impulse is in the 2nd subframe,
                 * - build exc[] in 2nd subframe
                 *--------------------------------------------------------*/

                if( PIT_MIN > (*position) )
                {
                    pit_start = L_SUBFR - (*position);
                }
                else
                {
                    pit_start = PIT_MIN;
                }

                if( pit_start < PIT_MIN )
                {
                    pit_start = PIT_MIN;
                }

                pit_limit = 2*pit_start + (*position);

                /* 7 bit pitch DECODER */
                index = (short)get_next_indice( st, nBits );

                *T0 = (short) (floor( pit_start + index/2 ));
                *T0_frac = (index - (*T0 - pit_start)*2 ) * 2;
                limit_T0( L_FRAME, 8, pit_flag, limit_flag, *T0, 0, T0_min, T0_max );   /* find T0_min and T0_max */

                /* Find the adaptive codebook vector - ACELP long-term prediction   */
                pred_lt4( &exc[i_subfr], &exc[i_subfr], *T0, *T0_frac, L_SUBFR+1, inter4_2, L_INTERPOL2, PIT_UP_SAMP );

                offset = tbe_celp_exc_offset(*T0, *T0_frac);
                for (i=0; i<L_SUBFR * HIBND_ACB_L_FAC; i++)
                {
                    bwe_exc[i + i_subfr * HIBND_ACB_L_FAC] = bwe_exc[i + i_subfr * HIBND_ACB_L_FAC - offset ];
                }

            }
            else if( (i_subfr == 2*L_SUBFR) && (tc_subfr == TC_0_128) )
            {
                /*--------------------------------------------------------*
                 * second glottal impulse is in the 3rd subframe
                 * - build exc[] in 3rd subframe
                 *--------------------------------------------------------*/

                /* 7bit pitch DECODER */
                pit_start = 2*L_SUBFR - (*position);

                index = (short)get_next_indice( st, nBits );

                *T0 = (short) (floor( pit_start + (short)(index/2) ));
                *T0_frac = ( index - (*T0 - pit_start)*2 ) * 2;

                limit_T0( L_FRAME, 8, pit_flag, limit_flag, *T0, 0, T0_min, T0_max );   /* find T0_min and T0_max */

                /* Find the adaptive codebook vector. ACELP long-term prediction */
                pred_lt4( &exc[i_subfr], &exc[i_subfr], *T0, *T0_frac, L_SUBFR+1, inter4_2, L_INTERPOL2, PIT_UP_SAMP );

                offset = tbe_celp_exc_offset(*T0, *T0_frac);
                for (i=0; i<L_SUBFR * HIBND_ACB_L_FAC; i++)
                {
                    bwe_exc[i + i_subfr * HIBND_ACB_L_FAC] = bwe_exc[i + i_subfr * HIBND_ACB_L_FAC - offset ];
                }

            }
            else if( (i_subfr == 2*L_SUBFR) && (tc_subfr == TC_0_192) )
            {
                /*--------------------------------------------------------*
                 * second glottal impulse is in the 4th subframe
                 * - build exc[] in 3rd subframe
                 *--------------------------------------------------------*/

                *T0 = 4*L_SUBFR;
                *T0_frac = 0;
                *Jopt_flag = 0;

                /* set adaptive part of exciation for curent subframe to 0 */
                set_f( &exc[i_subfr], 0, (short)(L_SUBFR+1) );

                set_f( &bwe_exc[i_subfr*HIBND_ACB_L_FAC], 0, (short)(L_SUBFR*HIBND_ACB_L_FAC) );
            }
            else if( (i_subfr == 3*L_SUBFR) && (tc_subfr == TC_0_192) )
            {
                /*--------------------------------------------------------*
                 * second glottal impulse is in the 4th subframe
                 * - build exc[] in 4th subframe
                 *--------------------------------------------------------*/

                pit_start = 3*L_SUBFR - (*position);
                pit_limit = 2*L_FRAME - PIT_MAX - 2*(*position) - 2;

                index = (short)get_next_indice( st, nBits );

                if( index < (pit_limit-pit_start)*2)
                {
                    *T0 = (short) (floor( pit_start + (index/2) ));
                    *T0_frac = (index - ((*T0) - pit_start)*2 )*2;
                }
                else
                {
                    *T0 = index + pit_limit - (pit_limit-pit_start)*2;
                    *T0_frac = 0;
                }

                /* biterror detection mechanism */
                if( ((*T0<<2) + *T0_frac) > (PIT_MAX<<2)+2 )
                {
                    *T0 = L_SUBFR;
                    *T0_frac = 0;
                    st->BER_detect = 1;
                }

                /* Find the adaptive codebook vector. ACELP long-term prediction   */
                pred_lt4( &exc[i_subfr], &exc[i_subfr], *T0, *T0_frac, L_SUBFR+1, inter4_2, L_INTERPOL2, PIT_UP_SAMP );

                offset = tbe_celp_exc_offset(*T0, *T0_frac);
                for (i=0; i<L_SUBFR * HIBND_ACB_L_FAC; i++)
                {
                    bwe_exc[i + i_subfr * HIBND_ACB_L_FAC] = bwe_exc[i + i_subfr * HIBND_ACB_L_FAC - offset ];
                }
            }
            else if( (i_subfr == 3*L_SUBFR) && (tc_subfr == TC_0_128) )
            {
                /*--------------------------------------------------------*
                 * second glottal impulse in the 3rd subframe
                 * build exc[] in 4th subframe
                 *--------------------------------------------------------*/

                index = (short)get_next_indice( st, nBits );

                delta_pit_dec( 2, index, T0, T0_frac, *T0_min );

                /* Find the adaptive codebook vector. ACELP long-term prediction   */
                pred_lt4( &exc[i_subfr], &exc[i_subfr], *T0, *T0_frac, L_SUBFR+1, inter4_2, L_INTERPOL2, PIT_UP_SAMP);

                offset = tbe_celp_exc_offset(*T0, *T0_frac);
                for (i=0; i<L_SUBFR * HIBND_ACB_L_FAC; i++)
                {
                    bwe_exc[i + i_subfr * HIBND_ACB_L_FAC] = bwe_exc[i + i_subfr * HIBND_ACB_L_FAC - offset ];
                }
            }

            /*------------------------------------------------------------*
             * first glottal impulse is NOT in the 1st subframe,
             * or two impulses are in the 1st subframe
             *------------------------------------------------------------*/
            else
            {
                index = (short)get_next_indice( st, nBits );

                pit_Q_dec( 0, index, nBits, 8, pit_flag, limit_flag, T0, T0_frac, T0_min, T0_max, &st->BER_detect );

                /* Find the adaptive codebook vector */
                pred_lt4( &exc[i_subfr], &exc[i_subfr], *T0, *T0_frac, L_SUBFR+1, inter4_2, L_INTERPOL2, PIT_UP_SAMP );

                offset = tbe_celp_exc_offset(*T0, *T0_frac);
                for (i=0; i<L_SUBFR * HIBND_ACB_L_FAC; i++)
                {
                    bwe_exc[i + i_subfr * HIBND_ACB_L_FAC] = bwe_exc[i + i_subfr * HIBND_ACB_L_FAC - offset];
                }
            }

            /*-----------------------------------------------------------------*
             * LP filtering of the adaptive excitation (if non-zero)
             *-----------------------------------------------------------------*/

            if( *Jopt_flag )
            {
                lp_filt_exc_dec( st, MODE1,core_brate, Opt_AMR_WB, coder_type, i_subfr, L_SUBFR, L_frame, 0, exc );
            }

            /*---------------------------------------------------------------------*
             * fill the pitch buffer - needed for post-processing and FEC_clas_estim()
             *---------------------------------------------------------------------*/

            **pt_pitch = (float)(*T0) + (float)(*T0_frac)/4.0f;   /* save subframe pitch values  */
            if( (tc_subfr >= 2*L_SUBFR) && (i_subfr == 3*L_SUBFR) )
            {
                (*pt_pitch) -= 3;
                **pt_pitch = (float)(*T0) + (float)(*T0_frac)/4.0f;
                (*pt_pitch)++;
                **pt_pitch = (float)(*T0) + (float)(*T0_frac)/4.0f;
                (*pt_pitch)++;
                **pt_pitch = (float)(*T0) + (float)(*T0_frac)/4.0f;
                (*pt_pitch)++;
            }
            else if( (tc_subfr == L_SUBFR) && (i_subfr == 2*L_SUBFR) )
            {
                (*pt_pitch) -= 2;
                **pt_pitch = (float)(*T0) + (float)(*T0_frac)/4.0f;
                (*pt_pitch)++;
                **pt_pitch = (float)(*T0) + (float)(*T0_frac)/4.0f;
                (*pt_pitch)++;
            }
            else if( (tc_subfr == TC_0_64) && (i_subfr == L_SUBFR) )
            {
                (*pt_pitch) -= 1;
                **pt_pitch = (float)(*T0) + (float)(*T0_frac)/4.0f;
                (*pt_pitch)++;
            }
            else if( (tc_subfr == TC_0_128) && (i_subfr == 2*L_SUBFR) )
            {
                (*pt_pitch) -= 2;
                **pt_pitch = (float)(*T0) + (float)(*T0_frac)/4.0f;
                (*pt_pitch)++;
                **pt_pitch = (float)(*T0) + (float)(*T0_frac)/4.0f;
                (*pt_pitch)++;
            }
            else if( (tc_subfr == TC_0_192) && (i_subfr == 3*L_SUBFR) )
            {
                (*pt_pitch) -= 3;
                **pt_pitch = (float)(*T0) + (float)(*T0_frac)/4.0f;
                (*pt_pitch)++;
                **pt_pitch = (float)(*T0) + (float)(*T0_frac)/4.0f;
                (*pt_pitch)++;
                **pt_pitch = (float)(*T0) + (float)(*T0_frac)/4.0f;
                (*pt_pitch)++;
            }

        }
        else  /* L_frame == L_FRAME16k */
        {
            if( i_subfr >= 2*L_SUBFR )
            {
                limit_flag = 1;
            }

            if( i_subfr - tc_subfr == L_SUBFR )
            {
                limit_T0( L_FRAME16k, 8, 0, limit_flag, *T0, *T0_frac, T0_min, T0_max );   /* find T0_min and T0_max */
            }

            /*-----------------------------------------------------------------*
             * get number of bits and index for pitch decoding
             *-----------------------------------------------------------------*/

            nBits = ACB_bits_16kHz_tbl[BIT_ALLOC_IDX_16KHZ(core_brate, coder_type, i_subfr, TC_SUBFR2IDX_16KHZ(tc_subfr))];

            index = (short)get_next_indice( st, nBits );

            /*-----------------------------------------------------------------*
             * Find adaptive part of excitation, encode pitch period
             *-----------------------------------------------------------------*/

            if( nBits == 10 )
            {
                pit16k_Q_dec( index, nBits, limit_flag, T0, T0_frac, T0_min, T0_max, &st->BER_detect );
            }
            else if( nBits == 8 )     /* tc_subfr==0 && i_subfr==L_SUBFR */
            {
                /*-----------------------------------------------------------------------------*
                 * The pitch range is encoded absolutely with 8 bits and is divided as follows:
                 *   PIT16k_MIN  to PIT16k_FR2_TC0_2SUBFR-1 resolution 1/4 (frac = 0,1,2 or 3)
                 *   PIT16k_FR2_TC0_2SUBFR to 2*L_SUBFR     resolution 1/2 (frac = 0 or 2)
                 *-----------------------------------------------------------------------------*/

                if( index < (PIT16k_FR2_TC0_2SUBFR-PIT16k_MIN)*4 )
                {
                    *T0 = PIT16k_MIN + (index/4);
                    *T0_frac = index - (*T0 - PIT16k_MIN)*4;
                }
                else
                {
                    index -=  (PIT16k_FR2_TC0_2SUBFR-PIT16k_MIN)*4;
                    *T0 = PIT16k_FR2_TC0_2SUBFR + (index/2);
                    *T0_frac = index - (*T0 - PIT16k_FR2_TC0_2SUBFR)*2;
                    (*T0_frac) *= 2;
                }

                /* biterror detection mechanism */
                if( ((*T0<<2) + *T0_frac) > ((2*L_SUBFR)<<2) )
                {
                    *T0 = L_SUBFR;
                    *T0_frac = 0;
                    st->BER_detect = 1;
                }
            }
            else if( nBits == 6 )
            {
                delta_pit_dec( 4, index, T0, T0_frac, *T0_min );
            }
            if( nBits == 6 )
            {
                limit_T0( L_FRAME16k, 8, L_SUBFR, limit_flag, *T0, *T0_frac, T0_min, T0_max );   /* find T0_min and T0_max */
            }

            /*-----------------------------------------------------------------*
             * - find the adaptive codebook vector
             * - LP filtering of the adaptive excitation (if non-zero)
             *-----------------------------------------------------------------*/

            if( (i_subfr == L_SUBFR) && (*T0 == 2*L_SUBFR) )
            {
                /* no adaptive excitation in the second subframe */
                set_f( &exc[i_subfr], 0, L_SUBFR+1 );

                get_next_indice( st, 1 );      /* this bit is actually not needed */

                set_f( &bwe_exc[i_subfr * 2], 0, L_SUBFR * 2 );
            }
            else
            {
                /* Find the adaptive codebook vector - ACELP long-term prediction */
                pred_lt4( &exc[i_subfr], &exc[i_subfr], *T0, *T0_frac, L_SUBFR+1, inter4_2, L_INTERPOL2, PIT_UP_SAMP);

                for (i=0; i<L_SUBFR * 2; i++)
                {
                    bwe_exc[i + i_subfr * 2] = bwe_exc[i + i_subfr * 2 - *T0 * 2 - (int) ((float) *T0_frac * 0.5f + 4 + 0.5f) + 4];
                }

                lp_filt_exc_dec( st, MODE1,core_brate, Opt_AMR_WB, coder_type, i_subfr, L_SUBFR, L_frame, 0, exc );

                *Jopt_flag = 1;
            }

            **pt_pitch = (float)(*T0) + (float)(*T0_frac)/4.0f;   /* save subframe pitch values  */

            /*---------------------------------------------------------------------*
             * fill the pitch buffer - needed for post-processing and FEC_clas_estim()
             *---------------------------------------------------------------------*/

            if( (i_subfr - tc_subfr == L_SUBFR) || (tc_subfr==0 && i_subfr==2*L_SUBFR) )
            {
                index = i_subfr/L_SUBFR;
                (*pt_pitch) -= index;

                for( i=0; i<index; i++ )
                {
                    **pt_pitch = (float)(*T0) + (float)(*T0_frac)/4.0f;
                    (*pt_pitch)++;
                }
            }
        }
    }

    return;
}

/*---------------------------------------------------------------------*
 *  tc_dec()
 *
 *  Principal function for TC decoding.
 *  - constructs glottal codebook contribution
 *  - uses pitch sharpening
 *  - uses gain_trans
 *---------------------------------------------------------------------*/

static void tc_dec(
    Decoder_State *st,            /* i/o: decoder state structure   */
    const short L_frame,        /* i  : length of the frame                         */
    float exc[],          /* o  : glottal codebook contribution               */
    short *T0,            /* o  : close-loop pitch period                     */
    short *T0_frac,       /* o  : close-loop pitch period - fractional part   */
    const short i_subfr,        /* i  : subframe index                              */
    const short tc_subfr,       /* i  : TC subframe index                           */
    short *position,      /* o  : first glottal impulse position in frame     */
    const long  core_brate,     /* i  : core bitrate                                */
    float bwe_exc[]       /* o  : excitation for SWB TBE                   */
)
{
    short i, imp_shape, imp_pos, imp_sign, imp_gain, nBits;
    float gain_trans;
    short index;

    /*----------------------------------------------------------------*
     * find the number of bits
     *----------------------------------------------------------------*/

    if( L_frame == L_FRAME )
    {
        nBits = ACB_bits_tbl[BIT_ALLOC_IDX(core_brate, TRANSITION, i_subfr, TC_SUBFR2IDX(tc_subfr))];
    }
    else  /* L_frame == L_FRAME16k */
    {
        nBits = ACB_bits_16kHz_tbl[BIT_ALLOC_IDX_16KHZ(core_brate, TRANSITION, i_subfr, TC_SUBFR2IDX_16KHZ(tc_subfr))];
    }

    /*----------------------------------------------------------------*
     * decode parameter T0 (pitch period)
     *----------------------------------------------------------------*/

    if( L_frame == L_FRAME )
    {
        if ( ((i_subfr == 0) && ((tc_subfr == 0) || (tc_subfr == TC_0_64) || (tc_subfr == TC_0_128) || (tc_subfr == TC_0_192) )) || (tc_subfr == L_SUBFR) )
        {
            *T0 = L_SUBFR;
            *T0_frac = 0;
        }
        else if( (tc_subfr == 3*L_SUBFR) )
        {
            i = (short)get_next_indice( st, nBits );

            if( nBits == 9 )
            {
                abs_pit_dec( 4, i, 0, T0, T0_frac );
            }
            else
            {
                abs_pit_dec( 2, i, 0, T0, T0_frac );
            }
        }
        else
        {
            i = (short)get_next_indice( st, nBits );

            if( i == 0 )
            {
                *T0 = L_SUBFR;
                *T0_frac = 0;
            }
            else
            {
                if( tc_subfr == TC_0_0 )
                {
                    delta_pit_dec( 2, i, T0, T0_frac, PIT_MIN-1 );
                }
                else
                {
                    delta_pit_dec( 0, i, T0, T0_frac, PIT_MIN-1 );
                }
            }
        }
    }
    else  /* L_frame == L_FRAME16k */
    {
        i = (short)get_next_indice( st, nBits );

        if( nBits == 10 )
        {
            if( i < (PIT16k_FR2_EXTEND_10b-PIT16k_MIN_EXTEND)*4 )
            {
                *T0 = PIT16k_MIN_EXTEND + (i/4);
                *T0_frac = i - ((*T0 - PIT16k_MIN_EXTEND)*4);
            }
            else
            {
                index =  i - (PIT16k_FR2_EXTEND_10b-PIT16k_MIN_EXTEND)*4;
                *T0 = PIT16k_FR2_EXTEND_10b + (index/2);
                *T0_frac = index - (*T0 - PIT16k_FR2_EXTEND_10b)*2;
                (*T0_frac) *= 2;
            }
        }
        else if( nBits == 6 )
        {
            *T0 = PIT16k_MIN + (i/2);
            *T0_frac = i - (*T0 - PIT16k_MIN)*2;
            *T0_frac *= 2;
        }
    }

    /*----------------------------------------------------------------*
     * decode other TC parameters
     *----------------------------------------------------------------*/

    imp_shape = (short)get_next_indice( st, 3 );
    imp_pos = (short)get_next_indice( st, 6 );
    imp_sign = (short)get_next_indice( st, 1 );
    imp_gain = (short)get_next_indice( st, 3 );

    /*----------------------------------------------------------------*
     * restore gain_trans
     * build glottal codebook contribution
     *----------------------------------------------------------------*/

    gain_trans = tbl_gain_trans_tc[imp_gain];

    if( imp_sign == 0 )
    {
        gain_trans *= -1;
    }

    /* build glottal codebook contribution */
    set_f( &exc[i_subfr], 0, L_SUBFR );

    for( i = (imp_pos-L_IMPULSE2); i <= (imp_pos+L_IMPULSE2); i++ )
    {
        if( (i >= 0) && (i < L_SUBFR) )
        {
            exc[i+i_subfr] = glottal_cdbk[(imp_shape)*L_IMPULSE+i-imp_pos+L_IMPULSE2]*gain_trans;
        }
    }

    /*--------------------------------------------------------------*
     * adapt. search of the second impulse in the same subframe
     * (when appears)
     *--------------------------------------------------------------*/

    pred_lt4_tc( exc, *T0, *T0_frac, inter4_2, imp_pos, i_subfr );

    if( L_frame == L_FRAME )
    {
        interp_code_5over2(&exc[i_subfr], &bwe_exc[i_subfr * HIBND_ACB_L_FAC], L_SUBFR);
    }
    else
    {
        interp_code_4over2(&exc[i_subfr], &bwe_exc[i_subfr * 2], L_SUBFR);
    }

    *position = imp_pos + i_subfr;

    return;
}
/*-------------------------------------------------------------------*
 * tc_classif()
 *
 * TC subframe classification decoding
 *-------------------------------------------------------------------*/

short tc_classif(
    Decoder_State *st,      /* i/o: decoder state structure           */
    const short L_frame   /* i  : length of the frame               */
)
{
    short tc_subfr, indice;

    if( L_frame == L_FRAME )
    {
        if ( get_next_indice( st, 1 ) )
        {
            tc_subfr = TC_0_0;
        }
        else
        {
            if ( get_next_indice( st, 1 ) )
            {
                tc_subfr = 0;

                if ( get_next_indice( st, 1 ) )
                {
                    tc_subfr = TC_0_192;
                }
                else
                {
                    if ( get_next_indice( st, 1 ) )
                    {
                        tc_subfr = TC_0_64;
                    }
                    else
                    {
                        tc_subfr = TC_0_128;
                    }
                }
            }
            else
            {
                if ( get_next_indice( st, 1 ) )
                {
                    tc_subfr = L_SUBFR;
                }
                else
                {
                    if ( get_next_indice( st, 1 ) )
                    {
                        tc_subfr = 2*L_SUBFR;
                    }
                    else
                    {
                        tc_subfr = 3*L_SUBFR;
                    }
                }
            }
        }
    }
    else  /* L_frame == L_FRAME16k */
    {
        indice = (short) get_next_indice( st, 2 );

        if( indice < 3 )
        {
            tc_subfr = indice * L_SUBFR;
        }
        else
        {
            if( get_next_indice( st, 1 ) == 0 )
            {
                tc_subfr = 3*L_SUBFR;
            }
            else
            {
                tc_subfr = 4*L_SUBFR;
            }
        }
    }

    return( tc_subfr );
}
