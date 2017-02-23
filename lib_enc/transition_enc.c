/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"

/*-----------------------------------------------------------------*
 * Local functions
 *-----------------------------------------------------------------*/

static void gain_trans_enc( float *gain_trans, float exc[], short *quant_index, short *quant_sign );

static void tc_enc( Encoder_State *st, const long core_brate, const short L_frame, const short i_subfr, short *tc_subfr, short *position,
                    const float *h1, const float *xn, float *exc, float *y1, short *T0_min, short *T0_max, short *T0,
                    short *T0_frac, float *gain_pit, float g_corr[], float *bwe_exc );

/*-----------------------------------------------------------------*
 * transition_enc()
 *
 * Principal function for adaptive excitation construction in TC
 *-----------------------------------------------------------------*/

void transition_enc(
    Encoder_State *st,             /* i/o: encoder state structure      */
    const long  core_brate,     /* i  : core bitrate                                */
    const short L_frame,        /* i  : length of the frame                         */
    const short coder_type,     /* i  : coding type                                 */
    const short i_subfr,        /* i  : subframe index                              */
    short *tc_subfr,      /* i/o: TC subframe index                           */
    short *Jopt_flag,     /* i  : joint optimization flag                     */
    short *position,      /* i/o: maximum of residual signal index            */
    const float voicing[],      /* i  : normalized correlations (from OL pitch)     */
    const short T_op[],         /* i  : open loop pitch estimates in current frame  */
    short *T0,            /* i/o: close loop integer pitch                    */
    short *T0_frac,       /* i/o: close loop fractional part of the pitch     */
    short *T0_min,        /* i/o: lower limit for close-loop search           */
    short *T0_max,        /* i/o: higher limit for close-loop search          */
    float *exc,           /* i/o: pointer to excitation signal frame          */
    float *y1,            /* o  : zero-memory filtered adaptive excitation    */
    const float *res,           /* i  : pointer to the LP residual signal frame     */
    const float *h1,            /* i  : weighted filter input response              */
    const float *xn,            /* i  : target vector                               */
    float *xn2,           /* o  : target vector for innovation search         */
    float *gp_cl,         /* i/o: memory of gain of pitch clipping algorithm  */
    float *gain_pit,      /* o  : adaptive excitation gain                    */
    float *g_corr,        /* o  : ACELP correlation values                    */
    short *clip_gain,     /* i/o: adaptive gain clipping flag                 */
    float **pt_pitch,     /* o  : floating pitch values                       */
    float *bwe_exc        /* o  : excitation for SWB TBE                      */
)
{
    short i, pit_flag, pit_start, pit_limit, index, nBits;
    float temp;
    short limit_flag, mult_Top;
    short lp_select, lp_flag;
    int offset=0;
    /* set limit_flag to 0 for restrained limits, and 1 for extended limits */
    limit_flag = 0;

    pit_start = PIT_MIN;

    /*-----------------------------------------------------------------*
     * TC: subrame determination for glottal shape search
     * -------------------------------------------------------
     * tc_subfr == 0         - TC in 1st subframe
     * tc_subfr == TC_0_0    - TC in 1st subframe + information about T0
     * tc_subfr == L_SUBFR   - TC in 2nd subframe
     * tc_subfr == 2*L_SUBFR - TC in 3rd subframe
     * tc_subfr == 3*L_SUBFR - TC in 4th subframe
     *-----------------------------------------------------------------*/

    if( i_subfr == 0 )
    {
        if( *tc_subfr == 3*L_SUBFR )
        {
            if( L_frame == L_FRAME )
            {
                *position = emaximum( res + 3*L_SUBFR,min(T_op[0]+2,L_SUBFR), &temp ) + 3*L_SUBFR;
                *tc_subfr = 3*L_SUBFR;
            }
            else /* L_frame == L_FRAME16k */
            {
                *position = emaximum( res + 4*L_SUBFR,min(T_op[0]+2,L_SUBFR), &temp ) + 4*L_SUBFR;
                *tc_subfr = 4*L_SUBFR;
            }
        }
        else
        {
            *position = emaximum( res, (short)(T_op[0]+2), &temp );

            /* correction in case of possibly wrong T_op (double-pitch values) */
            if( (L_frame == L_FRAME && T_op[0] > 2*PIT_MIN) ||
                    (L_frame == L_FRAME16k && T_op[0] > 2*PIT16k_MIN)
              )
            {
                short position_tmp, len;
                float aver, temp2;

                len = (short)(T_op[0]/2+2);

                position_tmp = emaximum( res, len, &temp2 );
                aver = dotp( res, res, len) + 0.01f;
                aver = sqrt(aver/len);

                temp = sqrt(temp);
                temp2 = sqrt(temp2);

                if( temp2 > 0.8f * temp && aver < 0.25f * temp )
                {
                    *position = position_tmp;
                }
            }

            *tc_subfr = (short) floor((*position)/L_SUBFR)*L_SUBFR;
        }

        mult_Top = 1;
        if( limit_flag == 0 )
        {
            if( L_frame == L_FRAME && T_op[1] < PIT_MIN )
            {
                mult_Top = 2;
            }

            if( L_frame == L_FRAME16k && T_op[1] < PIT16k_MIN )
            {
                mult_Top = 2;
            }
        }

        limit_T0( L_frame, 8, 0, limit_flag, mult_Top*T_op[1], 0, T0_min, T0_max );
    }

    /*-----------------------------------------------------------------*
     * zero adaptive excitation signal construction
     *-----------------------------------------------------------------*/

    if ( *tc_subfr > i_subfr )
    {
        *gain_pit = 0.0f;
        *clip_gain = 0;
        g_corr[0] = 1.0f;
        g_corr[1] = 1.0f;

        set_f(&exc[i_subfr], 0, L_SUBFR);  /* set excitation for current subrame to 0 */

        if( L_frame == L_FRAME )
        {
            set_f(&bwe_exc[i_subfr*HIBND_ACB_L_FAC], 0, (short) (L_SUBFR*HIBND_ACB_L_FAC));         /* set past excitation buffer to 0 */
        }
        else
        {
            set_f(&bwe_exc[i_subfr*2], 0, (short) (L_SUBFR*2));         /* set past excitation buffer to 0 */
        }

        set_f(y1, 0, L_SUBFR);             /* set filtered adaptive excitation to 0 */
        mvr2r(xn, xn2, L_SUBFR);           /* target vector for codebook search */
        *T0 = L_SUBFR;
        *T0_frac = 0;

        **pt_pitch = (float)(*T0) + (float)(*T0_frac)/4.0f;   /* save subframe pitch values  */
    }

    /*-----------------------------------------------------------------*
     * glottal codebook contribution construction
     *-----------------------------------------------------------------*/

    else if ( *tc_subfr == i_subfr )
    {
        if( L_frame == L_FRAME )
        {
            set_f( bwe_exc-PIT_MAX*HIBND_ACB_L_FAC, 0, PIT_MAX*HIBND_ACB_L_FAC);         /* set past excitation buffer to 0 */
        }
        else
        {
            set_f( bwe_exc-PIT16k_MAX*2, 0, PIT16k_MAX*2);         /* set past excitation buffer to 0 */
        }

        tc_enc( st, core_brate, L_frame, i_subfr, tc_subfr, position, h1, xn, exc, y1, T0_min, T0_max, T0, T0_frac, gain_pit, g_corr, bwe_exc );

        *clip_gain = gp_clip( core_brate, voicing, i_subfr, coder_type, xn, gp_cl );
        updt_tar( xn, xn2, y1, *gain_pit, L_SUBFR );

        **pt_pitch = (float)(*T0) + (float)(*T0_frac)/4.0f;   /* save subframe pitch values  */
        *Jopt_flag = 1;
    }

    /*--------------------------------------------------------------*
     * other subframes -> GENERIC encoding type,
     * standard adaptive excitation contribution
     * - exemption only in case when first glottal impulse is
     * in the 1st subframe and the second one in 2nd subframe
     * and later
     *--------------------------------------------------------------*/

    else if ( *tc_subfr < i_subfr)
    {
        if( L_frame == L_FRAME )
        {
            *Jopt_flag = 1;

            /* pit_flag for T0 bits number coding determination */
            if( ((i_subfr - *tc_subfr) == L_SUBFR) || ((i_subfr - *tc_subfr) == L_SUBFR-TC_0_0) )
            {
                pit_flag = 0;
            }
            else
            {
                pit_flag = L_SUBFR;
            }
            if( *tc_subfr == TC_0_0 )
            {
                if( i_subfr == L_SUBFR )
                {
                    limit_T0( L_FRAME, 8, pit_flag, limit_flag, *T0, 0, T0_min, T0_max );
                }
                pit_flag = 1;
            }

            /*----------------------------------------------------------*
             * if tc_subfr==0, change tc_subfr corresponding to the
             * second glot. impulse position
             *----------------------------------------------------------*/

            if( (*tc_subfr == 0) && (i_subfr == L_SUBFR) )
            {
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

                /* Find the closed loop pitch period */
                *T0 = pitch_fr4( &exc[i_subfr], xn, h1, *T0_min, *T0_max, T0_frac, pit_flag, limit_flag, pit_start, pit_limit, L_FRAME, L_SUBFR );

                offset = tbe_celp_exc_offset(*T0, *T0_frac);
                for (i=0; i<L_SUBFR * HIBND_ACB_L_FAC; i++)
                {
                    bwe_exc[i + i_subfr * HIBND_ACB_L_FAC] = bwe_exc[i + i_subfr * HIBND_ACB_L_FAC - offset ];
                }

                if( (*T0)>(2*L_SUBFR-(*position)) )
                {
                    if( (*T0) + (*position) >= 3*L_SUBFR)
                    {
                        /* second glottal impulse is in the 4th subframe */
                        *tc_subfr = TC_0_192;
                    }
                    else
                    {
                        /* second glottal impulse is in the 3rd subframe */
                        *tc_subfr = TC_0_128;
                    }
                }
                else if( (*tc_subfr == 0) && (i_subfr == L_SUBFR) )
                {
                    /* second glottal impulse is in the 2nd subframe */
                    *tc_subfr = TC_0_64;
                }
            }

            /*-----------------------------------------------------------------*
             * get number of bits for pitch encoding
             *-----------------------------------------------------------------*/

            nBits = ACB_bits_tbl[BIT_ALLOC_IDX(core_brate, TRANSITION, i_subfr, TC_SUBFR2IDX(*tc_subfr))];

            /*-----------------------------------------------------------------*
             * Find adaptive part of excitation, encode pitch period
             *-----------------------------------------------------------------*/

            /* first glottal impulse is in the 1st subrame */
            if( (i_subfr == L_SUBFR) && (*tc_subfr >= TC_0_128) )
            {
                /*--------------------------------------------------------*
                 * second glottal impulse is in the 3rd or 4th subframe
                 * - build exc[] in 2nd subframe
                 *--------------------------------------------------------*/

                *T0 = 2*L_SUBFR;
                *T0_frac = 0;
                *Jopt_flag = 0;

                set_f( &exc[i_subfr], 0, (short)(L_SUBFR+1) );
                set_f( &bwe_exc[i_subfr*HIBND_ACB_L_FAC], 0, (short)(L_SUBFR*HIBND_ACB_L_FAC) );
            }
            else if( (i_subfr == L_SUBFR) && (*tc_subfr == TC_0_64) )
            {
                /*--------------------------------------------------------*
                 * second glottal impulse is in the 2nd subframe,
                 * - build exc[] in 2nd subframe
                 *--------------------------------------------------------*/

                if( (*T0)+(*position) < L_SUBFR )
                {
                    /* impulse must be in the 2nd subframe (not in 1st) */
                    *T0 = L_SUBFR - (*position);
                    *T0_frac = 0;
                }

                if( (*T0)+(*position) >= 2*L_SUBFR )
                {
                    /* impulse must be in the 2nd subframe (not in 3rd) */
                    *T0 = 2*L_SUBFR-1-(*position);
                    *T0_frac = 2;
                }

                limit_T0( L_FRAME, 8, pit_flag, limit_flag, *T0, 0, T0_min, T0_max ); /* find T0_min and T0_max for delta search */

                /* 7bit ENCODER */
                index = (*T0-pit_start)*2 + *T0_frac/2;
                push_indice( st, IND_PITCH, index, nBits );

                /* Find the adaptive codebook vector - ACELP long-term prediction */
                pred_lt4( &exc[i_subfr], &exc[i_subfr], *T0, *T0_frac, L_SUBFR+1, inter4_2, L_INTERPOL2, PIT_UP_SAMP);

                offset = tbe_celp_exc_offset(*T0, *T0_frac);
                for (i=0; i<L_SUBFR * HIBND_ACB_L_FAC; i++)
                {
                    bwe_exc[i + i_subfr * HIBND_ACB_L_FAC] = bwe_exc[i + i_subfr * HIBND_ACB_L_FAC - offset ];
                }
            }
            else if( (i_subfr == 2*L_SUBFR) && (*tc_subfr == TC_0_128) )
            {
                /*--------------------------------------------------------*
                 * second glottal impulse is in the 3rd subframe
                 * - build exc[] in 3rd subframe
                 *--------------------------------------------------------*/

                pit_start = 2*L_SUBFR - (*position);

                pit_flag = 0;

                *T0 = pitch_fr4( &exc[i_subfr], xn, h1, *T0_min, *T0_max, T0_frac, pit_flag, limit_flag, pit_start, 3*L_SUBFR, L_FRAME, L_SUBFR );

                if( (*T0)+(*position) < 2*L_SUBFR )
                {
                    /* impulse must be in the 3rd subframe (not in 2nd) */
                    *T0 = 2*L_SUBFR - (*position);
                    *T0_frac = 0;
                }
                if( (*T0)+(*position) >= 3*L_SUBFR )
                {
                    /* impulse must be in the 3rd subframe (not in 4th) */
                    *T0 = 3*L_SUBFR - 1 - (*position);
                    *T0_frac = 2;
                }

                limit_T0( L_FRAME, 8, pit_flag, limit_flag, *T0, 0, T0_min, T0_max ); /* find T0_min and T0_max for delta search */

                index = (*T0-pit_start)*2 + *T0_frac/2;
                push_indice( st, IND_PITCH, index, nBits );

                /* Find the adaptive codebook vector - ACELP long-term prediction */
                pred_lt4( &exc[i_subfr], &exc[i_subfr], *T0, *T0_frac, L_SUBFR+1, inter4_2, L_INTERPOL2, PIT_UP_SAMP );

                offset = tbe_celp_exc_offset(*T0, *T0_frac);
                for (i=0; i<L_SUBFR * HIBND_ACB_L_FAC; i++)
                {
                    bwe_exc[i + i_subfr * HIBND_ACB_L_FAC] = bwe_exc[i + i_subfr * HIBND_ACB_L_FAC - offset ];
                }
            }
            else if( (i_subfr == 2*L_SUBFR) && (*tc_subfr == TC_0_192) )
            {
                /*--------------------------------------------------------*
                 * second glottal impulse is in the 4th subframe
                 * - build exc[] in 3rd subframe
                 *--------------------------------------------------------*/

                *T0 = 4*L_SUBFR;
                *T0_frac = 0;
                *Jopt_flag = 0;

                set_f( &exc[i_subfr], 0, (short)(L_SUBFR+1) );
                set_f( &bwe_exc[i_subfr*HIBND_ACB_L_FAC], 0, (short)(L_SUBFR*HIBND_ACB_L_FAC) );
            }
            else if( (i_subfr == 3*L_SUBFR) && (*tc_subfr == TC_0_192) )
            {
                /*--------------------------------------------------------*
                 * second glottal impulse is in the 4th subframe
                 * - build exc[] in 4th subframe
                 *--------------------------------------------------------*/
                /* always T0_frac = 0 */
                pit_flag = 0;

                *T0 = pitch_fr4( &exc[i_subfr], xn, h1, *T0_min, *T0_max, T0_frac, pit_flag, limit_flag, PIT_MIN, PIT_FR1_8b, L_FRAME, L_SUBFR );

                if( (*T0)+(*position) < 3*L_SUBFR )
                {
                    /* impulse must be in the 4th subframe (not in 3rd) */
                    *T0 = 3*L_SUBFR - (*position);
                    *T0_frac = 0;
                }

                pit_start = 3*L_SUBFR - (*position);
                pit_limit = 2*L_FRAME - PIT_MAX - 2 - 2*(*position);

                if( *T0 < pit_limit )
                {
                    index = (*T0-pit_start)*2 + *T0_frac/2;
                }
                else
                {
                    index = *T0 - pit_limit + (pit_limit-pit_start)*2;
                    *T0_frac = 0;
                }

                push_indice( st, IND_PITCH, index, nBits );

                /* Find the adaptive codebook vector - ACELP long-term prediction */
                pred_lt4( &exc[i_subfr], &exc[i_subfr], *T0, *T0_frac, L_SUBFR+1, inter4_2, L_INTERPOL2, PIT_UP_SAMP);

                offset = tbe_celp_exc_offset(*T0, *T0_frac);
                for (i=0; i<L_SUBFR * HIBND_ACB_L_FAC; i++)
                {
                    bwe_exc[i + i_subfr * HIBND_ACB_L_FAC] = bwe_exc[i + i_subfr * HIBND_ACB_L_FAC - offset ];
                }
            }
            else if( (i_subfr == 3*L_SUBFR) && (*tc_subfr == TC_0_128) )
            {
                /*--------------------------------------------------------*
                 * second glottal impulse in the 3rd subframe
                 * build exc[] in 4th subframe
                 *--------------------------------------------------------*/

                pit_flag = L_SUBFR;
                *T0 = pitch_fr4( &exc[i_subfr], xn, h1, *T0_min, *T0_max, T0_frac, pit_flag, limit_flag, PIT_MIN, PIT_FR1_8b, L_FRAME, L_SUBFR );
                index = delta_pit_enc( 2, *T0, *T0_frac, *T0_min );
                push_indice( st, IND_PITCH, index, nBits );

                /* Find the adaptive codebook vector - ACELP long-term prediction */
                pred_lt4( &exc[i_subfr], &exc[i_subfr], *T0, *T0_frac, L_SUBFR+1, inter4_2, L_INTERPOL2, PIT_UP_SAMP );

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
                if( nBits == 8 || nBits == 5 )
                {
                    if( !((*tc_subfr == 0) && (i_subfr == L_SUBFR)) )
                    {
                        *T0 = pitch_fr4( &exc[i_subfr], xn, h1, *T0_min, *T0_max, T0_frac, pit_flag, limit_flag, PIT_MIN, PIT_FR1_8b, L_FRAME, L_SUBFR );
                    }
                }
                else
                {
                    *T0 = pitch_fr4( &exc[i_subfr], xn, h1, *T0_min, *T0_max, T0_frac, pit_flag, limit_flag, PIT_FR2_9b, PIT_FR1_9b, L_FRAME, L_SUBFR );
                }

                pit_Q_enc( st, 0, nBits, 8, pit_flag, limit_flag, *T0, *T0_frac, T0_min, T0_max );

                /* Find the adaptive codebook vector - ACELP long-term prediction */
                pred_lt4( &exc[i_subfr], &exc[i_subfr], *T0, *T0_frac, L_SUBFR+1, inter4_2, L_INTERPOL2, PIT_UP_SAMP );

                offset = tbe_celp_exc_offset(*T0, *T0_frac);
                for (i=0; i<L_SUBFR * HIBND_ACB_L_FAC; i++)
                {
                    bwe_exc[i + i_subfr * HIBND_ACB_L_FAC] = bwe_exc[i + i_subfr * HIBND_ACB_L_FAC - offset ];
                }
            }

            if( *Jopt_flag == 0 )
            {
                /* adaptive/TC excitation is zero */
                mvr2r( xn, xn2, L_SUBFR );
                g_corr[0] = 0.0f;
                g_corr[1] = 0.0f;
                *clip_gain = 0;
            }
            else
            {
                *clip_gain = gp_clip( core_brate, voicing, i_subfr, coder_type, xn, gp_cl );

                lp_select = lp_filt_exc_enc( MODE1, core_brate, 0, coder_type, i_subfr, exc, h1,
                                             xn, y1, xn2, L_SUBFR, L_frame, g_corr, *clip_gain, gain_pit, &lp_flag );

                if( lp_flag == NORMAL_OPERATION )
                {
                    push_indice( st, IND_LP_FILT_SELECT, lp_select, 1 );
                }
            }

            **pt_pitch = (float)(*T0) + (float)(*T0_frac)/4.0f;   /* save subframe pitch values  */

            /*---------------------------------------------------------------------*
             * fill the pitch buffer - needed for post-processing
             *---------------------------------------------------------------------*/

            if( (*tc_subfr >= 2*L_SUBFR) && (i_subfr == 3*L_SUBFR) )
            {
                (*pt_pitch) -= 3;
                **pt_pitch = (float)(*T0) + (float)(*T0_frac)/4.0f;
                (*pt_pitch)++;
                **pt_pitch = (float)(*T0) + (float)(*T0_frac)/4.0f;
                (*pt_pitch)++;
                **pt_pitch = (float)(*T0) + (float)(*T0_frac)/4.0f;
                (*pt_pitch)++;
            }
            else if( (*tc_subfr == L_SUBFR) && (i_subfr == 2*L_SUBFR) )
            {
                (*pt_pitch) -= 2;
                **pt_pitch = (float)(*T0) + (float)(*T0_frac)/4.0f;
                (*pt_pitch)++;
                **pt_pitch = (float)(*T0) + (float)(*T0_frac)/4.0f;
                (*pt_pitch)++;
            }
            else if( (*tc_subfr == TC_0_64) && (i_subfr == L_SUBFR) )
            {
                (*pt_pitch) -= 1;
                **pt_pitch = (float)(*T0) + (float)(*T0_frac)/4.0f;
                (*pt_pitch)++;
            }
            else if( (*tc_subfr == TC_0_128) && (i_subfr == 2*L_SUBFR) )
            {
                (*pt_pitch) -= 2;
                **pt_pitch = (float)(*T0) + (float)(*T0_frac)/4.0f;
                (*pt_pitch)++;
                **pt_pitch = (float)(*T0) + (float)(*T0_frac)/4.0f;
                (*pt_pitch)++;
            }
            else if( (*tc_subfr == TC_0_192) && (i_subfr == 3*L_SUBFR) )
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

            if( i_subfr <= 2*L_SUBFR )
            {
                if( i_subfr < 2*L_SUBFR )
                {
                    mult_Top = 1;
                    if( T_op[0] < PIT16k_MIN )
                    {
                        mult_Top = 2;
                    }
                    limit_T0( L_FRAME16k, 8, 0, limit_flag, mult_Top*T_op[0], 0, T0_min, T0_max );  /* TC0 second subfr. */
                }
                else
                {
                    limit_T0( L_FRAME16k, 8, 0, limit_flag, T_op[1], 0, T0_min, T0_max );  /* TC0 third subfr., or TC64 third subfr. */
                }
            }

            /*-----------------------------------------------------------------*
             * get number of bits for pitch encoding
             *-----------------------------------------------------------------*/

            nBits = ACB_bits_16kHz_tbl[BIT_ALLOC_IDX_16KHZ(core_brate, coder_type, i_subfr, TC_SUBFR2IDX_16KHZ(*tc_subfr))];

            /*-----------------------------------------------------------------*
             * Find adaptive part of excitation, encode pitch period
             *-----------------------------------------------------------------*/

            if( nBits == 10 )
            {
                *T0 = pitch_fr4( &exc[i_subfr], xn, h1, *T0_min, *T0_max, T0_frac, 0, limit_flag, PIT16k_FR2_EXTEND_10b, PIT16k_MAX, L_frame, L_SUBFR );
                pit16k_Q_enc( st, nBits, limit_flag, *T0, *T0_frac, T0_min, T0_max );
            }
            else if( nBits == 8 )     /* tc_subfr==0 && i_subfr==L_SUBFR */
            {
                /*-----------------------------------------------------------------------------*
                 * The pitch range is encoded absolutely with 8 bits and is divided as follows:
                 *   PIT16k_MIN  to PIT16k_FR2_TC0_2SUBFR-1 resolution 1/4 (frac = 0,1,2 or 3)
                 *   PIT16k_FR2_TC0_2SUBFR to 2*L_SUBFR     resolution 1/2 (frac = 0 or 2)
                 *-----------------------------------------------------------------------------*/

                *T0 = pitch_fr4( &exc[i_subfr], xn, h1, *T0_min, *T0_max, T0_frac, 0, limit_flag, PIT16k_FR2_TC0_2SUBFR, 2*L_SUBFR, L_frame, L_SUBFR );

                if( *T0_max > 2*L_SUBFR )
                {
                    *T0 = 2*L_SUBFR;
                    *T0_frac = 0;
                }

                if( *T0 < PIT16k_FR2_TC0_2SUBFR )
                {
                    index = (*T0)*4 + (*T0_frac) - (PIT16k_MIN*4);
                }
                else
                {
                    index = (*T0)*2 + ((*T0_frac)>>1) - (PIT16k_FR2_TC0_2SUBFR*2) + ((PIT16k_FR2_TC0_2SUBFR-PIT16k_MIN)*4);
                }

                push_indice( st, IND_PITCH, index, nBits );
            }
            else if( nBits == 6 )
            {
                /* delta search */
                *T0 = pitch_fr4( &exc[i_subfr], xn, h1, *T0_min, *T0_max, T0_frac, L_SUBFR, limit_flag, PIT16k_FR2_EXTEND_9b, PIT16k_FR1_EXTEND_9b, L_frame, L_SUBFR );

                index = delta_pit_enc( 4, *T0, *T0_frac, *T0_min );

                push_indice( st, IND_PITCH, index, nBits );
            }
            if( nBits == 6 )
            {
                limit_T0( L_FRAME16k, 8, L_SUBFR, limit_flag, *T0, *T0_frac, T0_min, T0_max );
            }

            /*-----------------------------------------------------------------*
             * - gain clipping test to avoid unstable synthesis
             * - LP filtering of the adaptive excitation
             * - codebook target computation
             *-----------------------------------------------------------------*/

            if( (i_subfr == L_SUBFR) && (*T0 == 2*L_SUBFR) )
            {
                *gain_pit = 0.0f;
                *clip_gain = 0;
                g_corr[0] = 0.01f;
                g_corr[1] = 0.01f;
                *Jopt_flag = 0;

                set_f( &exc[i_subfr], 0, L_SUBFR+1 );   /* set excitation for current subrame to 0 */

                push_indice( st, IND_LP_FILT_SELECT, 0, 1 );      /* this bit is actually not needed */

                mvr2r( xn, xn2, L_SUBFR );              /* target vector for codebook search */
                set_f( y1, 0, L_SUBFR );                /* set filtered adaptive excitation to 0 */
                set_f( &bwe_exc[i_subfr * 2], 0, L_SUBFR * 2 );
            }
            else
            {
                /* Find the adaptive codebook vector - ACELP long-term prediction */
                pred_lt4( &exc[i_subfr], &exc[i_subfr], *T0, *T0_frac, L_SUBFR+1, inter4_2, L_INTERPOL2, PIT_UP_SAMP );

                for (i=0; i<L_SUBFR * 2; i++)
                {
                    bwe_exc[i + i_subfr * 2] = bwe_exc[i + i_subfr * 2 - *T0 * 2 - (int) ((float) *T0_frac * 0.5f + 4 + 0.5f) + 4];
                }

                *clip_gain = gp_clip( core_brate, voicing, i_subfr, coder_type, xn, gp_cl );

                lp_select = lp_filt_exc_enc( MODE1, core_brate, 0, coder_type, i_subfr, exc, h1,
                                             xn, y1, xn2, L_SUBFR, L_frame, g_corr, *clip_gain, gain_pit, &lp_flag );

                if( lp_flag == NORMAL_OPERATION )
                {
                    push_indice( st, IND_LP_FILT_SELECT, lp_select, 1 );
                }

                *Jopt_flag = 1;
            }

            **pt_pitch = (float)(*T0) + (float)(*T0_frac)/4.0f;   /* save subframe pitch value  */

            /*---------------------------------------------------------------------*
             * fill the pitch buffer - needed for post-processing
             *---------------------------------------------------------------------*/

            if( (i_subfr - *tc_subfr == L_SUBFR) || (*tc_subfr==0 && i_subfr==2*L_SUBFR) )
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



/*-------------------------------------------------------------------------------------------*
 * tc_enc()
 *
 * Principal function for transition coding (TC) in encoder.
 * Glottal codebook contribution part:
 *
 *           |----|             |----|                                      xn
 *  imp_pos->||   |  imp_shape->| g1 |                                       |
 *           | |  |             | g2 |     ----   exc  |---|  y1   ----      |
 *           |  | |-------------|    |----|gain|-------| h |------|gain|----(-)---> xn2
 *           |   ||             | gn |     ----        |---|       ----
 *           |----|             |----|
 *          codebook          excitation gain_trans    h_orig     gain_pit
 *
 *-------------------------------------------------------------------------------------------*/

static void tc_enc(
    Encoder_State *st,            /* i/o: encoder state structure                 */
    const long  core_brate,     /* i  : core bitrate                            */
    const short L_frame,        /* i  : length of the frame                     */
    const short i_subfr,        /* i  : subrame index                           */
    short *tc_subfr,      /* i/o: TC subframe index                       */
    short *position,      /* i/o: index of the residual signal maximum    */
    const float *h1,            /* i  : weighted filter input response          */
    const float *xn,            /* i  : target signal                           */
    float *exc,           /* o  : glottal codebook contribution           */
    float *y1,            /* o  : filtered glottal codebook contribution  */
    short *T0_min,        /* o  : lower pitch limit                       */
    short *T0_max,        /* o  : higher pitch limit                      */
    short *T0,            /* o  : close loop integer pitch                */
    short *T0_frac,       /* o  : close loop fractional part of the pitch */
    float *gain_pit,      /* o  : pitch gain  (0..GAIN_PIT_MAX)           */
    float g_corr[],       /* o  : correlations <y1,y1>  and -2<xn,y1>     */
    float *bwe_exc        /* i/o: excitation for SWB TBE                  */
)
{
    short imp_shape, imp_pos, imp_sign, imp_gain, index, nBits;
    float gain_trans;

    imp_pos = *position-i_subfr;

    /*-----------------------------------------------------------------*
     * get number of bits for pitch encoding
     *-----------------------------------------------------------------*/

    if( L_frame == L_FRAME )
    {
        nBits = ACB_bits_tbl[BIT_ALLOC_IDX(core_brate, TRANSITION, i_subfr, TC_SUBFR2IDX(*tc_subfr))];
    }
    else  /* L_frame == L_FRAME16k */
    {
        nBits = ACB_bits_16kHz_tbl[BIT_ALLOC_IDX_16KHZ(core_brate, TRANSITION, i_subfr, TC_SUBFR2IDX_16KHZ(*tc_subfr))];
    }

    /*--------------------------------------------------------------*
     * Closed loop pitch search
     *--------------------------------------------------------------*/

    *T0_frac = 0;
    if( L_frame == L_FRAME )
    {
        if( (*T0_min <= L_SUBFR) || (*tc_subfr == 3*L_SUBFR) )
        {
            if( nBits == 9 )
            {
                *T0 = pitch_fr4( &exc[i_subfr], xn, h1, *T0_min, *T0_max, T0_frac, 0, 0, PIT_FR2_9b, PIT_FR1_9b, L_FRAME, L_SUBFR );
            }
            else if( nBits == 6 )
            {
                *T0 = pitch_fr4( &exc[i_subfr], xn, h1, *T0_min, *T0_max, T0_frac, 0, 0, PIT_MIN, L_SUBFR, L_FRAME, L_SUBFR );
            }
            else
            {
                *T0 = pitch_fr4( &exc[i_subfr], xn, h1, *T0_min, *T0_max, T0_frac, 0, 0, PIT_MAX, PIT_MIN, L_FRAME, L_SUBFR );
            }
        }
        else
        {
            *T0 = L_SUBFR;
        }

        if( (*tc_subfr == L_SUBFR) && (*T0 < L_SUBFR) )
        {
            *T0 = L_SUBFR;
        }
    }
    else  /* L_frame == L_FRAME16k */
    {
        if( nBits == 10 )
        {
            *T0 = pitch_fr4( &exc[i_subfr], xn, h1, *T0_min, *T0_max, T0_frac, 0, 1, PIT16k_FR2_EXTEND_10b, PIT16k_MAX, L_FRAME16k, L_SUBFR );
        }
        else if( nBits == 6 )
        {
            /* T0_frac with 1/2 sample resolution */
            *T0 = pitch_fr4( &exc[i_subfr], xn, h1, *T0_min, *T0_max, T0_frac, 0, 0, PIT16k_MIN, L_SUBFR, L_FRAME16k, L_SUBFR );

            if( *T0 > L_SUBFR )
            {
                *T0 = L_SUBFR;
                *T0_frac = 0;
            }
        }
    }

    /* set tc_subfr to TC_0_0 */
    if( i_subfr == 0 && L_frame == L_FRAME && ( *T0 < L_SUBFR || *tc_subfr == 3*L_SUBFR ) )
    {
        *tc_subfr = TC_0_0;
    }

    /*--------------------------------------------------------------*
     * Builds glottal codebook contribution
     *--------------------------------------------------------------*/

    set_impulse( xn, h1, &exc[i_subfr], y1, &imp_shape, &imp_pos, &gain_trans );

    /*--------------------------------------------------------------*
     * quantize gain_trans and scale glottal codebook contribution
     *--------------------------------------------------------------*/

    gain_trans_enc( &gain_trans, &exc[i_subfr], &imp_gain, &imp_sign );

    /* set past excitation buffer to zeros */
    set_f( exc-L_EXC_MEM, 0, L_EXC_MEM );

    /*--------------------------------------------------------------*
     * adapt. search of the second impulse in the same subframe
     * (when appears)
     *--------------------------------------------------------------*/

    pred_lt4_tc( exc, *T0, *T0_frac, inter4_2, imp_pos, i_subfr);

    if( L_frame == L_FRAME )
    {
        interp_code_5over2(&exc[i_subfr], &bwe_exc[i_subfr * HIBND_ACB_L_FAC], L_SUBFR);
    }
    else
    {
        interp_code_4over2(&exc[i_subfr], &bwe_exc[i_subfr * 2], L_SUBFR);
    }

    /*--------------------------------------------------------------*
     * compute glottal-shape codebook excitation
     *--------------------------------------------------------------*/

    /* create filtered glottal codebook contribution */
    conv( &exc[i_subfr], h1, y1, L_SUBFR );

    /* gain_pit computation */
    *gain_pit = corr_xy1( xn, y1, g_corr, L_SUBFR,0);

    /*--------------------------------------------------------------*
     * Encode parameters and write indices
     *--------------------------------------------------------------*/

    if( L_frame == L_FRAME )
    {
        if( ( (i_subfr != 0) || (*tc_subfr == TC_0_0) )
                && (*tc_subfr != L_SUBFR))
        {
            /* write pitch index */
            if( (*T0 >= L_SUBFR) && (*tc_subfr != 3*L_SUBFR) )
            {
                push_indice( st, IND_PITCH, 0, nBits );
            }
            else if( *tc_subfr == 3*L_SUBFR )
            {
                if( nBits == 9 )
                {
                    index = abs_pit_enc( 4, 0, *T0, *T0_frac );
                }
                else
                {
                    index = abs_pit_enc( 2, 0, *T0, *T0_frac );
                }

                push_indice( st, IND_PITCH, index, nBits );

                limit_T0( L_FRAME, 8, 0, 0, *T0, 0, T0_min, T0_max );
            }
            else
            {
                if( nBits == 6 )
                {
                    index = delta_pit_enc( 2, *T0, *T0_frac, PIT_MIN-1 );
                    push_indice( st, IND_PITCH, index, nBits );
                }
                else
                {
                    index = delta_pit_enc( 0, *T0, *T0_frac, PIT_MIN-1 );
                    push_indice( st, IND_PITCH, index, nBits );
                }
            }
        }
    }
    else  /* L_frame == L_FRAME16k */
    {
        if( nBits == 10 )
        {
            pit16k_Q_enc( st, nBits, 1, *T0, *T0_frac, T0_min, T0_max );
        }
        else if( nBits == 6 )
        {
            index = 2*(*T0 - PIT16k_MIN) + *T0_frac/2;
            push_indice( st, IND_PITCH, index, nBits );
        }
    }

    push_indice( st, IND_TC_IMP_SHAPE, imp_shape, 3 );
    push_indice( st, IND_TC_IMP_POS, imp_pos, 6 );
    push_indice( st, IND_TC_IMP_SIGN, imp_sign, 1 );
    push_indice( st, IND_TC_IMP_GAIN, imp_gain, 3 );

    *position = imp_pos + i_subfr;

    return;
}

/*-----------------------------------------------------------------*
 * gain_trans_enc()
 *
 * Quantize gain_trans of TC (gains of glottal impulses).
 * - Uses scalar quantization prototypes tbl_gain_trans_tc[N_GAIN_TC].
 * - Gains the glottal codebook contibution signal.
 *-----------------------------------------------------------------*/

static void gain_trans_enc(
    float       *gain_trans,    /*  i/o: gain for TC                   */
    float       exc[],          /*  i/o: glottal codebook contribution */
    short       *quant_index,   /*  o  : index of quantized gain_trans */
    short       *quant_sign     /*  o  : sign of quantized gain_trans  */
)
{
    short i;


    if( *gain_trans < 0 )
    {
        *gain_trans *= -1;
        *quant_sign = 0;
    }
    else
    {
        *quant_sign = 1;
    }

    *quant_index = N_GAIN_TC-1;

    for( i=0; i<N_GAIN_TC-1; i++ )
    {
        if( *gain_trans < tbl_gain_trans_tc[i] )
        {
            *quant_index = i;
            break;
        }
    }

    /* restore gain_trans */
    *gain_trans = tbl_gain_trans_tc[i];

    if( *quant_sign == 0 )
    {
        *gain_trans *= -1;
    }

    for( i=0; i < L_SUBFR; i++ )
    {
        exc[i] *= (*gain_trans);
    }

    return;
}
