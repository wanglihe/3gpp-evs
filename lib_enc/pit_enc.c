/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "cnst.h"
#include "prot.h"
#include "rom_enc.h"
#include "rom_com.h"


/*------------------------------------------------------------------*
 * pit_encode()
 *
 * Close-loop pitch lag search and pitch lag quantization
 * Adaptive excitation construction
 *------------------------------------------------------------------*/

float pit_encode(                 /* o  : Fractional pitch for each subframe          */
    Encoder_State *st,               /* i/o: encoder state structure                     */
    const long  core_brate,       /* i  : core bitrate                                */
    const short Opt_AMR_WB,       /* i  : flag indicating AMR-WB IO mode              */
    const short L_frame,          /* i  : length of the frame                         */
    const short coder_type,       /* i  : coding type                                 */
    short *limit_flag,      /* i/o: restrained(0) or extended(1) Q limits       */
    const short i_subfr,          /* i  : subframe index                              */
    float *exc,             /* i/o: pointer to excitation signal frame          */
    const short L_subfr,          /* i  : subframe length                             */
    const short *T_op,            /* i  : open loop pitch estimates in current frame  */
    short *T0_min,          /* i/o: lower limit for close-loop search           */
    short *T0_max,          /* i/o: higher limit for close-loop search          */
    short *T0,              /* i/o: close loop integer pitch                    */
    short *T0_frac,         /* i/o: close loop fractional part of the pitch     */
    const float *h1,              /* i  : weighted filter input response              */
    const float *xn               /* i  : target vector                               */
)
{
    float pitch;
    short pit_flag, delta, mult_Top, nBits;

    /*----------------------------------------------------------------*
     * Initializations
     *----------------------------------------------------------------*/

    /* Set pit_flag to 0 for every subframe with absolute pitch search */
    pit_flag = i_subfr;
    if( i_subfr == 2*L_SUBFR )
    {
        pit_flag = 0;
    }

    /*-----------------------------------------------------------------*
     * Limit range of pitch search
     * Fractional pitch search
     * Pitch quantization
     *-----------------------------------------------------------------*/

    mult_Top = 1;

    if( !Opt_AMR_WB )
    {
        /*----------------------------------------------------------------*
         * Set limit_flag to 0 for restrained limits, and 1 for extended limits
         *----------------------------------------------------------------*/

        if( i_subfr == 0 )
        {
            *limit_flag = 1;
            if( coder_type == VOICED )
            {
                *limit_flag = 2;     /* double-extended limits */
            }

            if( coder_type == GENERIC && core_brate == ACELP_7k20 )
            {
                *limit_flag = 0;
            }
        }
        else if( i_subfr == 2*L_SUBFR && coder_type == GENERIC && core_brate <= ACELP_13k20 )
        {
            if( *T0 > (PIT_FR1_EXTEND_8b + PIT_MIN)>>1 )
            {
                *limit_flag = 0;
            }
        }

        /* check the minimum pitch value */
        if( *limit_flag == 0 )
        {
            if( ( i_subfr == 0 && T_op[0] < PIT_MIN ) ||
                    ( i_subfr == 2*L_SUBFR && T_op[1] < PIT_MIN ) )
            {
                mult_Top = 2;
            }
        }

        /*-------------------------------------------------------*
         *  Retrieve the number of Q bits
         *-------------------------------------------------------*/

        nBits = 0;
        if( coder_type != AUDIO )
        {
            /* find the number of bits */
            if( L_frame == L_FRAME )
            {
                nBits = ACB_bits_tbl[BIT_ALLOC_IDX(core_brate, coder_type, i_subfr, 0)];
            }
            else  /* L_frame == L_FRAME16k */
            {
                nBits = ACB_bits_16kHz_tbl[BIT_ALLOC_IDX_16KHZ(core_brate, coder_type, i_subfr, 0)];
            }
        }

        if( coder_type == AUDIO )
        {
            /*-------------------------------------------------------*
             *  Pitch encoding in AUDIO mode
             *  (both ACELP@12k8 and ACELP@16k cores)
             *-------------------------------------------------------*/

            delta = 4;

            if( L_subfr == L_frame/2 && i_subfr != 0 )
            {
                pit_flag = L_SUBFR;
            }

            if( pit_flag == 0 )
            {
                nBits = 10;
            }
            else
            {
                nBits = 6;
            }

            /* pitch lag search limitation */
            if( i_subfr == 0 )
            {
                limit_T0( L_FRAME, delta, pit_flag, *limit_flag, mult_Top*T_op[0], 0, T0_min, T0_max );
            }
            else if( i_subfr == 2*L_SUBFR && pit_flag == 0 )
            {
                limit_T0( L_FRAME, delta, pit_flag, *limit_flag, mult_Top*T_op[1], 0, T0_min, T0_max );
            }

            /* search and encode the closed loop pitch period */
            *T0 = pitch_fr4( &exc[i_subfr], xn, h1, *T0_min, *T0_max, T0_frac, pit_flag, *limit_flag, PIT_MAX, PIT_MAX, L_FRAME, L_subfr );

            pit_Q_enc( st, 0, nBits, delta, pit_flag, *limit_flag, *T0, *T0_frac, T0_min, T0_max );
        }
        else if( coder_type == VOICED )
        {
            /*-------------------------------------------------------*
             *  Pitch encoding in VOICED mode (ACELP@12k8 core only)
             *-------------------------------------------------------*/

            delta = 4;

            if( i_subfr == 2*L_SUBFR )
            {
                pit_flag = i_subfr;
            }

            /* pitch lag search limitation */
            if( i_subfr == 0 )
            {
                limit_T0( L_FRAME, delta, pit_flag, *limit_flag, mult_Top*T_op[0], 0, T0_min, T0_max );
            }

            /* search and encode the closed loop pitch period */
            if( nBits == 9 || nBits == 5 )
            {
                *T0 = pitch_fr4( &exc[i_subfr], xn, h1, *T0_min, *T0_max, T0_frac, pit_flag, *limit_flag, PIT_FR2_DOUBLEEXTEND_9b, PIT_FR1_DOUBLEEXTEND_9b, L_FRAME, L_SUBFR );
            }
            else if( nBits == 10 )
            {
                *T0 = pitch_fr4( &exc[i_subfr], xn, h1, *T0_min, *T0_max, T0_frac, pit_flag, *limit_flag, PIT_MAX, PIT_MAX, L_FRAME, L_SUBFR );
            }

            pit_Q_enc( st, 0, nBits, delta, pit_flag, *limit_flag, *T0, *T0_frac, T0_min, T0_max );
        }
        else
        {
            /*-------------------------------------------------------*
             *  Pitch encoding in GENERIC mode
             *  (both ACELP@12k8 and ACELP@16k cores)
             *-------------------------------------------------------*/

            delta = 8;

            /* pitch lag search limitation */
            if( i_subfr == 0 )
            {
                limit_T0( L_frame, delta, pit_flag, *limit_flag, mult_Top*T_op[0], 0, T0_min, T0_max );
            }
            else if( i_subfr == 2*L_SUBFR )
            {
                limit_T0( L_frame, delta, pit_flag, *limit_flag, mult_Top*T_op[1], 0, T0_min, T0_max );
            }

            /* search and encode the closed loop pitch period */
            if( L_frame == L_FRAME )
            {
                if( nBits == 8 || nBits == 5 )
                {
                    if( *limit_flag == 0 )
                    {
                        *T0 = pitch_fr4( &exc[i_subfr], xn, h1, *T0_min, *T0_max, T0_frac, pit_flag, *limit_flag, PIT_MIN, PIT_FR1_8b, L_FRAME, L_SUBFR );
                    }
                    else
                    {
                        *T0 = pitch_fr4( &exc[i_subfr], xn, h1, *T0_min, *T0_max, T0_frac, pit_flag, *limit_flag, PIT_MIN_EXTEND, PIT_FR1_EXTEND_8b, L_FRAME, L_SUBFR );
                    }
                }
                else if( nBits == 9 || nBits == 6 )
                {
                    if( *limit_flag == 0 )
                    {
                        *T0 = pitch_fr4( &exc[i_subfr], xn, h1, *T0_min, *T0_max, T0_frac, pit_flag, *limit_flag, PIT_FR2_9b, PIT_FR1_9b, L_FRAME, L_SUBFR );
                    }
                    else
                    {
                        *T0 = pitch_fr4( &exc[i_subfr], xn, h1, *T0_min, *T0_max, T0_frac, pit_flag, *limit_flag, PIT_FR2_EXTEND_9b, PIT_FR1_EXTEND_9b, L_FRAME, L_SUBFR );
                    }
                }
                else if( nBits == 10 )
                {
                    *T0 = pitch_fr4( &exc[i_subfr], xn, h1, *T0_min, *T0_max, T0_frac, pit_flag, *limit_flag, PIT_MAX, PIT_MAX, L_FRAME, L_SUBFR );
                }

                pit_Q_enc( st, 0, nBits, delta, pit_flag, *limit_flag, *T0, *T0_frac, T0_min, T0_max );
            }
            else  /* L_frame == L_FRAME16k */
            {
                if( nBits == 9 || nBits == 6 )
                {
                    *T0 = pitch_fr4( &exc[i_subfr], xn, h1, *T0_min, *T0_max, T0_frac, pit_flag, *limit_flag, PIT16k_FR2_EXTEND_9b, PIT16k_FR1_EXTEND_9b, L_FRAME16k, L_SUBFR );
                }
                else if( nBits == 10 )
                {
                    *T0 = pitch_fr4( &exc[i_subfr], xn, h1, *T0_min, *T0_max, T0_frac, pit_flag, *limit_flag, PIT16k_FR2_EXTEND_10b, PIT16k_MAX, L_FRAME16k, L_SUBFR );
                }

                pit16k_Q_enc( st, nBits, *limit_flag, *T0, *T0_frac, T0_min, T0_max );
            }
        }
    }

    /*-------------------------------------------------------*
     *  Pitch encoding in AMR-WB IO mode
     *-------------------------------------------------------*/

    else
    {
        delta = 8;
        *limit_flag = 0;

        if( core_brate == ACELP_6k60 )
        {
            nBits = 5;

            /* pitch lag search limitation */
            if( i_subfr == 0 )
            {
                limit_T0( L_FRAME, delta, pit_flag, *limit_flag, mult_Top*T_op[0], 0, T0_min, T0_max );
                nBits = 8;
            }

            if( i_subfr == 2*L_SUBFR )
            {
                /* rewrite pit_flag - it must not be zero */
                pit_flag = i_subfr;
            }

            /* search and encode the closed loop pitch period */
            *T0 = pitch_fr4( &exc[i_subfr], xn, h1, *T0_min, *T0_max, T0_frac, pit_flag, *limit_flag, PIT_MIN, PIT_FR1_8b, L_FRAME, L_SUBFR );
        }
        else if( core_brate == ACELP_8k85 )
        {
            nBits = 5;

            /* pitch lag search limitation */
            if( i_subfr == 0 )
            {
                limit_T0( L_FRAME, delta, pit_flag, *limit_flag, mult_Top*T_op[0], 0, T0_min, T0_max );
                nBits = 8;
            }
            else if( i_subfr == 2*L_SUBFR )
            {
                limit_T0( L_FRAME, delta, pit_flag, *limit_flag, mult_Top*T_op[1], 0, T0_min, T0_max );
                nBits = 8;
            }

            /* search and encode the closed loop pitch period */
            *T0 = pitch_fr4( &exc[i_subfr], xn, h1, *T0_min, *T0_max, T0_frac, pit_flag, *limit_flag, PIT_MIN, PIT_FR1_8b, L_FRAME, L_SUBFR );
        }
        else
        {
            nBits = 6;

            /* pitch lag search limitation */
            if( i_subfr == 0 )
            {
                limit_T0( L_FRAME, delta, pit_flag, *limit_flag, mult_Top*T_op[0], 0, T0_min, T0_max );
                nBits = 9;
            }
            else if( i_subfr == 2*L_SUBFR )
            {
                limit_T0( L_FRAME, delta, pit_flag, *limit_flag, mult_Top*T_op[1], 0, T0_min, T0_max );
                nBits = 9;
            }
            else
            {
                limit_T0( L_FRAME, delta, pit_flag, 0, *T0, 0, T0_min, T0_max );         /* T0_frac==0 to keep IO with AMR-WB */
            }

            /* search and encode the closed loop pitch period */
            *T0 = pitch_fr4( &exc[i_subfr], xn, h1, *T0_min, *T0_max, T0_frac, pit_flag, *limit_flag, PIT_FR2_9b, PIT_FR1_9b, L_FRAME, L_SUBFR );
        }

        pit_Q_enc( st, 1, nBits, delta, pit_flag, *limit_flag, *T0, *T0_frac, T0_min, T0_max );
    }

    /*-------------------------------------------------------*
     * Compute floating pitch output
     *-------------------------------------------------------*/

    pitch = (float)(*T0) + (float)(*T0_frac)/4.0f;   /* save subframe pitch values  */

    return pitch;
}

/*-------------------------------------------------------------------*
 * pitch_fr4()
 *
 * Find the closed loop pitch period with 1/4 subsample resolution.
 *-------------------------------------------------------------------*/

short pitch_fr4(            /* o  : chosen integer pitch lag                   */
    const float exc[],      /* i  : excitation buffer                          */
    const float xn[],       /* i  : target signal                              */
    const float h[],        /* i  : weighted synthesis filter impulse response */
    const short t0_min,     /* i  : minimum value in the searched range.       */
    const short t0_max,     /* i  : maximum value in the searched range.       */
    short *pit_frac,  /* o  : chosen fraction (0, 1, 2 or 3)             */
    const short i_subfr,    /* i  : flag to first subframe                     */
    const short limit_flag, /* i  : flag for limits (0=restrained, 1=extended) */
    const short t0_fr2,     /* i  : minimum value for resolution 1/2           */
    const short t0_fr1,     /* i  : minimum value for resolution 1             */
    const short L_frame,    /* i  : length of the frame                        */
    const short L_subfr     /* i  : size of subframe                           */
)
{
    short i, fraction, step;
    short t0, t1, t_min, t_max, pit_min;
    float cor_max, max, temp;
    float *corr, corr_v[15+2*L_INTERPOL1+1];

    /* initialization */
    if( limit_flag == 0 )
    {
        if( L_frame == L_FRAME )
        {
            pit_min = PIT_MIN;
        }
        else /* L_frame == L_FRAME16k */
        {
            pit_min = PIT16k_MIN;
        }
    }
    else
    {
        if( L_frame == L_FRAME )
        {
            pit_min = PIT_MIN_EXTEND;
            if( limit_flag == 2 )
            {
                pit_min = PIT_MIN_DOUBLEEXTEND;
            }
        }
        else /* L_frame == L_FRAME16k */
        {
            pit_min = PIT16k_MIN_EXTEND;
        }
    }


    /*-----------------------------------------------------------------*
     * - Find interval to compute normalized correlation
     * - allocate memory to normalized correlation vector
     * - Compute normalized correlation between target and filtered
     *   excitation
     *-----------------------------------------------------------------*/

    t_min = t0_min - L_INTERPOL1;
    t_max = t0_max + L_INTERPOL1;
    corr  = &corr_v[-t_min];      /* corr[t_min..t_max] */

    norm_corr( exc, xn, h, t_min, t_max, corr, L_subfr );

    /*-----------------------------------------------------------------*
     * Find integer pitch
     *-----------------------------------------------------------------*/

    max = corr[t0_min];
    t0  = t0_min;

    for( i=t0_min+1; i<=t0_max; i++ )
    {
        if( corr[i] >= max)
        {
            max = corr[i];
            t0 = i;
        }
    }

    if( t0_fr1 == pit_min )
    {
        /* don't search fraction (for 7b/4b quant) */
        if((i_subfr == 0) && (t0 >= t0_fr2))
        {
            i = (t0>>1)*2;       /* 2 samples resolution */
            if( (i+2) > PIT_MAX )
            {
                i -= 2;
            }
            if( corr[i] > corr[i+2] )
            {
                t0 = i;
            }
            else
            {
                t0 = i+2;
            }
        }

        *pit_frac = 0;

        return(t0);
    }
    if( (i_subfr == 0) && (t0 >= t0_fr1) )
    {
        *pit_frac = 0;

        return(t0);
    }

    /*------------------------------------------------------------------*
     * Search fractionnal pitch with 1/4 subsample resolution.
     * search the fractions around t0 and choose the one which maximizes
     * the interpolated normalized correlation.
     *-----------------------------------------------------------------*/

    t1 = t0;
    step = 1;                /* 1/4 subsample resolution */
    fraction = 1;
    if( ((i_subfr == 0) && (t0 >= t0_fr2)) || (t0_fr2 == pit_min) )
    {
        step = 2;            /* 1/2 subsample resolution */
        fraction = 2;
    }

    if( t0 == t0_min )        /* Limit case */
    {
        fraction = 0;
        cor_max = interpolation( &corr[t0], E_ROM_inter4_1, fraction, PIT_UP_SAMP, 4 );
    }
    else /* Process negative fractions */
    {
        t0--;
        cor_max = interpolation( &corr[t0], E_ROM_inter4_1, fraction, PIT_UP_SAMP, 4 );
        for( i=(fraction+step); i<=3; i=i+step )
        {
            temp = interpolation( &corr[t0], E_ROM_inter4_1, i, PIT_UP_SAMP, 4 );
            if (temp > cor_max)
            {
                cor_max = temp;
                fraction = i;
            }
        }
    }

    for( i=0; i<=3; i=i+step )    /* Process positive fractions */
    {
        temp = interpolation( &corr[t1], E_ROM_inter4_1, i, PIT_UP_SAMP, 4 );
        if( temp > cor_max )
        {
            cor_max = temp;
            fraction = i;
            t0 = t1;
        }
    }

    *pit_frac = fraction;

    return (t0);
}

/*-------------------------------------------------------------------*
 * norm_corr()
 *
 * Find the normalized correlation between the target vector and the
 * filtered past excitation (correlation between target and filtered
 * excitation divided by the square root of energy of filtered
 * excitation)
 *---------------------------------------------------------------------*/

void norm_corr(
    const float exc[],        /* i  : excitation buffer                          */
    const float xn[],         /* i  : target signal                              */
    const float h[],          /* i  : weighted synthesis filter impulse response */
    const short t_min,        /* i  : minimum value of searched range            */
    const short t_max,        /* i  : maximum value of searched range            */
    float corr_norm[],  /* o  : normalized correlation                     */
    const short L_subfr       /* i  : subframe size                              */
)
{
    short t, j, k;
    float excf[L_FRAME];      /* filtered past excitation */    /* length up to L_FRAME in GSC */
    float alp, ps, norm;

    k = - t_min;

    /*-----------------------------------------------------------------*
     * compute the filtered excitation for the first delay t_min
     *-----------------------------------------------------------------*/

    conv( &exc[k], h, excf, L_subfr );

    /*-----------------------------------------------------------------*
     * loop for every possible period
     *-----------------------------------------------------------------*/

    for (t = t_min; t <= t_max; t++)
    {
        /* Compute correlation between xn[] and excf[] */

        ps = 0.0f;
        for (j = 0; j < L_subfr; ++j)
        {
            ps += xn[j]*excf[j];
        }

        /* Compute 1/sqrt(energie of excf[]) */

        alp = 0.01f;
        for (j = 0; j < L_subfr; ++j)
        {
            alp += excf[j]*excf[j];
        }
        norm = inv_sqrt(alp);

        /* Normalize correlation = correlation * (1/sqrt(energie)) */
        corr_norm[t] = ps*norm;

        /* update the filtered excitation excf[] for the next iteration */
        if (t != t_max)
        {
            k--;
            for (j = L_subfr-1; j > 0; j--)
            {
                excf[j] = excf[j-1] + exc[k]*h[j];
            }
            excf[0] = exc[k];
        }
    }
    return;
}

/*-------------------------------------------------------------------*
 * abs_pit_enc()
 *
 * Encode pitch lag absolutely with resolution for shortest pitches
 * depending on parameter 'fr_step':
 * fr_step = 2: pitch range encoded with 8 bits
 * fr_step = 4: pitch range encoded with 9 bits
 *-------------------------------------------------------------------*/

short abs_pit_enc(            /* o  : pitch index                         */
    const short fr_steps,     /* i  : fractional resolution step          */
    const short limit_flag,   /* i  : restrained(0) or extended(1) limits */
    const short T0,           /* i  : integer pitch lag                   */
    const short T0_frac       /* i  : pitch fraction                      */
)
{
    short pitch_index;

    if( limit_flag == 0 )
    {
        if( fr_steps == 2 )
        {
            /*-----------------------------------------------------------------*
             * The pitch range is encoded absolutely with 8 bits
             * and is divided as follows:
             *   PIT_MIN to PIT_FR1_8b-1  resolution 1/2 (frac = 0 or 2)
             *   PIT_FR1_8b to PIT_MAX    resolution 1   (frac = 0)
             *-----------------------------------------------------------------*/
            if( T0 < PIT_FR1_8b )
            {
                pitch_index = T0*2 + (T0_frac>>1) - (PIT_MIN*2);
            }
            else
            {
                pitch_index = T0 - PIT_FR1_8b + ((PIT_FR1_8b-PIT_MIN)*2);
            }

        }
        else if( fr_steps == 4 )
        {
            /*-------------------------------------------------------------------*
             * The pitch range is encoded absolutely with 9 bits
             * and is divided as follows:
             *   PIT_MIN    to PIT_FR2_9b-1  resolution 1/4 (frac = 0,1,2 or 3)
             *   PIT_FR2_9b to PIT_FR1_9b-1  resolution 1/2 (frac = 0 or 2)
             *   PIT_FR1_9b to PIT_MAX       resolution 1   (frac = 0)
             *-------------------------------------------------------------------*/
            if( T0 < PIT_FR2_9b )
            {
                pitch_index = T0*4 + T0_frac - (PIT_MIN*4);
            }
            else if( T0 < PIT_FR1_9b )
            {
                pitch_index = T0*2 + (T0_frac>>1) - (PIT_FR2_9b*2) + ((PIT_FR2_9b-PIT_MIN)*4);
            }
            else
            {
                pitch_index = T0 - PIT_FR1_9b + ((PIT_FR2_9b-PIT_MIN)*4) + ((PIT_FR1_9b-PIT_FR2_9b)*2);
            }

        }
        else  /* fr_step == 0 */
        {
            /* not used in the codec */
            pitch_index = 0;
        }
    }
    else if( limit_flag == 1 )   /* extended Q range */
    {
        if( fr_steps == 2 )
        {
            /*-----------------------------------------------------------------*
             * The pitch range is encoded absolutely with 8 bits
             * and is divided as follows:
             *   PIT_MIN_EXTEND to PIT_FR1_EXTEND_8b-1  resolution 1/2 (frac = 0 or 2)
             *   PIT_FR1_EXTEND_8b to PIT_MAX           resolution 1   (frac = 0)
             *-----------------------------------------------------------------*/

            if( T0 < PIT_FR1_EXTEND_8b )
            {
                pitch_index = T0*2 + (T0_frac>>1) - (PIT_MIN_EXTEND*2);
            }
            else
            {
                pitch_index = T0 - PIT_FR1_EXTEND_8b + ((PIT_FR1_EXTEND_8b-PIT_MIN_EXTEND)*2);
            }

        }
        else if( fr_steps == 4 )
        {
            /*-------------------------------------------------------------------*
             * The pitch range is encoded absolutely with 9 bits
             * and is divided as follows:
             *   PIT_MIN_EXTEND    to PIT_FR2__EXTEND9b-1  resolution 1/4 (frac = 0,1,2 or 3)
             *   PIT_FR2_EXTEND_9b to PIT_FR1__EXTEND9b-1  resolution 1/2 (frac = 0 or 2)
             *   PIT_FR1_EXTEND_9b to PIT_MAX              resolution 1   (frac = 0)
             *-------------------------------------------------------------------*/

            if( T0 < PIT_FR2_EXTEND_9b )
            {
                pitch_index = T0*4 + T0_frac - (PIT_MIN_EXTEND*4);
            }
            else if( T0 < PIT_FR1_EXTEND_9b )
            {
                pitch_index = T0*2 + (T0_frac>>1) - (PIT_FR2_EXTEND_9b*2) + ((PIT_FR2_EXTEND_9b-PIT_MIN_EXTEND)*4);
            }
            else
            {
                pitch_index = T0 - PIT_FR1_EXTEND_9b + ((PIT_FR2_EXTEND_9b-PIT_MIN_EXTEND)*4) + ((PIT_FR1_EXTEND_9b-PIT_FR2_EXTEND_9b)*2);
            }

        }
        else  /* fr_step == 0 */
        {
            /* not used in the codec */
            pitch_index = 0;
        }
    }
    else  /* double-extended Q range */
    {
        if( fr_steps == 2 )
        {
            /*-----------------------------------------------------------------*
             * The pitch range is encoded absolutely with 8 bits
             * and is divided as follows:
             *   PIT_MIN_DOUBLEEXTEND    to PIT_FR1_DOUBLEEXTEND_8b-1  resolution 1/2 (frac = 0 or 2)
             *   PIT_FR1_DOUBLEEXTEND_8b to PIT_MAX                    resolution 1   (frac = 0)
             *-----------------------------------------------------------------*/

            if( T0 < PIT_FR1_DOUBLEEXTEND_8b )
            {
                pitch_index = T0*2 + (T0_frac>>1) - (PIT_MIN_DOUBLEEXTEND*2);
            }
            else
            {
                pitch_index = T0 - PIT_FR1_DOUBLEEXTEND_8b + ((PIT_FR1_DOUBLEEXTEND_8b-PIT_MIN_DOUBLEEXTEND)*2);
            }
        }
        else if( fr_steps == 4 )
        {
            /*-------------------------------------------------------------------*
             * The pitch range is encoded absolutely with 9 bits
             * and is divided as follows:
             *   PIT_MIN_DOUBLEEXTEND    to PIT_FR2_DOUBLEEXTEND9b-1  resolution 1/4 (frac = 0,1,2 or 3)
             *   PIT_FR2_DOUBLEEXTEND_9b to PIT_FR1_DOOBLEEXTEND9b-1  resolution 1/2 (frac = 0 or 2)
             *   PIT_FR1_DOUBLEEXTEND_9b to PIT_MAX                   resolution 1   (frac = 0)
             *-------------------------------------------------------------------*/

            if( T0 < PIT_FR2_DOUBLEEXTEND_9b )
            {
                pitch_index = T0*4 + T0_frac - (PIT_MIN_DOUBLEEXTEND*4);
            }
            else if( T0 < PIT_FR1_DOUBLEEXTEND_9b )
            {
                pitch_index = T0*2 + (T0_frac>>1) - (PIT_FR2_DOUBLEEXTEND_9b*2) + ((PIT_FR2_DOUBLEEXTEND_9b-PIT_MIN_DOUBLEEXTEND)*4);
            }
            else
            {
                pitch_index = T0 - PIT_FR1_DOUBLEEXTEND_9b + ((PIT_FR2_DOUBLEEXTEND_9b-PIT_MIN_DOUBLEEXTEND)*4) + ((PIT_FR1_DOUBLEEXTEND_9b-PIT_FR2_DOUBLEEXTEND_9b)*2);
            }
        }
        else  /* fr_step == 0 */
        {
            /* not used in the codec */
            pitch_index = 0;
        }
    }

    return pitch_index;
}

/*-------------------------------------------------------------------*
 * delta_pit_enc()
 *
 * Encode pitch lag differentially from T0_min to T0_max
 * with resolution depending on parameter 'fr_step':
 * fr_step = 0: resolution 1   (frac = 0), or
 * fr_step = 2: resolution 1/2 (frac = 0 or 2), or
 * fr_step = 4: resolution 1/4 (frac = 0, 1, 2, or 3)
 *-------------------------------------------------------------------*/

short delta_pit_enc(          /* o  : pitch index                             */
    const short fr_steps,     /* i  : fractional resolution step              */
    const short T0,           /* i  : integer pitch lag                       */
    const short T0_frac,      /* i  : pitch fraction                          */
    const short T0_min        /* i  : delta search min                        */
)
{
    short pitch_index;

    if( fr_steps == 0 )
    {
        pitch_index = T0 - T0_min;
    }
    else if( fr_steps == 2 )
    {
        pitch_index = (T0 - T0_min) * 2 + (T0_frac>>1);
    }
    else /* fr_steps == 4 */
    {
        pitch_index = (T0 - T0_min) * 4 + T0_frac;
    }

    return pitch_index;
}

/*-------------------------------------------------------------------*
 * pit_Q_enc()
 *
 * Encode subframe pitch lag
 *-------------------------------------------------------------------*/

void pit_Q_enc(
    Encoder_State *st,          /* i/o: encoder state structure      */
    const short Opt_AMR_WB,   /* i  : flag indicating AMR-WB IO mode          */
    const short nBits,        /* i  : # of Q bits                             */
    const short delta,        /* i  : Half the CL searched interval           */
    const short pit_flag,     /* i  : absolute(0) or delta(1) pitch Q         */
    const short limit_flag,   /* i  : restrained(0) or extended(1) Q limits   */
    const short T0,           /* i  : integer pitch lag                       */
    const short T0_frac,      /* i  : pitch fraction                          */
    short *T0_min,      /* i/o: delta search min                        */
    short *T0_max       /* o  : delta search max                        */
)
{
    short pitch_index;

    if( nBits == 10 )         /* absolute encoding with 10 bits */
    {
        if( limit_flag == 0 )
        {
            pitch_index = T0*4 + T0_frac - (PIT_MIN*4);
        }
        else if( limit_flag == 1 )
        {
            pitch_index = T0*4 + T0_frac - (PIT_MIN_EXTEND*4);
        }
        else  /* limit_flag == 2 */
        {
            pitch_index = T0*4 + T0_frac - (PIT_MIN_DOUBLEEXTEND*4);
        }
    }
    else if( nBits == 9 )     /* absolute encoding with 9 bits */
    {
        pitch_index = abs_pit_enc( 4, limit_flag, T0, T0_frac );

        /* find T0_min and T0_max for delta search */
        if( Opt_AMR_WB )
        {
            limit_T0( L_FRAME, delta, pit_flag, 0, T0, 0, T0_min, T0_max );         /* T0_frac==0 to keep IO with AMR-WB */
        }
    }
    else if( nBits == 8 )     /* absolute encoding with 8 bits */
    {
        pitch_index = abs_pit_enc( 2, limit_flag, T0, T0_frac );

        /* find T0_min and T0_max for delta search */
        if( Opt_AMR_WB )
        {
            limit_T0( L_FRAME, delta, pit_flag, 0, T0, 0, T0_min, T0_max );         /* T0_frac==0 to keep IO with AMR-WB */
        }
    }
    else if( nBits == 6 )     /* relative encoding with 6 bits */
    {
        pitch_index = delta_pit_enc( 4, T0, T0_frac, *T0_min );
    }
    else if( nBits == 5 )     /* relative encoding with 5 bits */
    {
        if( delta == 8 )
        {
            pitch_index = delta_pit_enc( 2, T0, T0_frac, *T0_min );
        }
        else  /* delta == 4 */
        {
            pitch_index = delta_pit_enc( 4, T0, T0_frac, *T0_min );
        }
    }
    else  /* nBits == 4 ) */  /* relative encoding with 4 bits */
    {
        if( delta == 8 )
        {
            pitch_index = delta_pit_enc( 0, T0, T0_frac, *T0_min );
        }
        else  /* delta == 4 */
        {
            pitch_index = delta_pit_enc( 2, T0, T0_frac, *T0_min );
        }
    }

    if( !Opt_AMR_WB )
    {
        /* find T0_min and T0_max for delta search */
        limit_T0( L_FRAME, delta, L_SUBFR, limit_flag, T0, T0_frac, T0_min, T0_max );
    }

    push_indice( st, IND_PITCH, pitch_index, nBits );

    return;
}

/*-------------------------------------------------------------------*
 * pit16k_Q_enc()
 *
 * Encode subframe pitch lag @16kHz core
 *-------------------------------------------------------------------*/

void pit16k_Q_enc(
    Encoder_State *st,          /* i/o: encoder state structure      */
    const short nBits,        /* i  : # of Q bits                             */
    const short limit_flag,   /* i  : restrained(0) or extended(1) Q limits   */
    const short T0,           /* i  : integer pitch lag                       */
    const short T0_frac,      /* i  : pitch fraction                          */
    short *T0_min,      /* i/o: delta search min                        */
    short *T0_max       /* o  : delta search max                        */
)
{
    short pitch_index;

    if( nBits == 10 )         /* absolute encoding with 10 bits */
    {
        {
            if( T0 < PIT16k_FR2_EXTEND_10b )
            {
                pitch_index = T0*4 + T0_frac - (PIT16k_MIN_EXTEND*4);
            }
            else
            {
                pitch_index = T0*2 + (T0_frac>>1) - (PIT16k_FR2_EXTEND_10b*2) + ((PIT16k_FR2_EXTEND_10b-PIT16k_MIN_EXTEND)*4);
            }
        }

        push_indice( st, IND_PITCH, pitch_index, nBits );
    }
    else if( nBits == 9 )     /* absolute encoding with 9 bits */
    {
        {
            /*-------------------------------------------------------------------*
             * The pitch range is encoded absolutely with 9 bits
             * and is divided as follows:
             *   PIT16k_EXTEND_MIN    to PIT16k_FR2_EXTEND_9b-1  resolution 1/4 (frac = 0,1,2 or 3)
             *   PIT16k_FR2_EXTEND_9b to PIT16k_FR1_EXTEND_9b-1  resolution 1/2 (frac = 0 or 2)
             *   PIT16k_FR1_EXTEND_9b to PIT16k_MAX              resolution 1   (frac = 0)
             *-------------------------------------------------------------------*/

            if( T0 < PIT16k_FR2_EXTEND_9b )
            {
                pitch_index = T0*4 + T0_frac - (PIT16k_MIN_EXTEND*4);
            }
            else if( T0 < PIT16k_FR1_EXTEND_9b )
            {
                pitch_index = T0*2 + (T0_frac>>1) - (PIT16k_FR2_EXTEND_9b*2) + ((PIT16k_FR2_EXTEND_9b-PIT16k_MIN_EXTEND)*4);
            }
            else
            {
                pitch_index = T0 - PIT16k_FR1_EXTEND_9b + ((PIT16k_FR2_EXTEND_9b-PIT16k_MIN_EXTEND)*4) + ((PIT16k_FR1_EXTEND_9b-PIT16k_FR2_EXTEND_9b)*2);
            }
        }

        push_indice( st, IND_PITCH, pitch_index, 9 );
    }
    else  /* nBits == 6 */    /* relative encoding with 6 bits */
    {
        pitch_index = (T0 - *T0_min) * 4 + T0_frac;

        push_indice( st, IND_PITCH, pitch_index, nBits );
    }

    limit_T0( L_FRAME16k, 8, L_SUBFR, limit_flag, T0, T0_frac, T0_min, T0_max );

    return;
}


/*------------------------------------------------------------------*
 * limit_T0_voiced2:
 *
 *
 *------------------------------------------------------------------*/

static void limit_T0_voiced2(
    int res,
    const short *T_op,
    int *T0_min,
    int *T0_min_frac,
    int *T0_max,
    int *T0_max_frac,
    short pit_min,
    short pit_max,
    int i_subfr
)
{
    int t, temp1, temp2;

    /* Lower-bound */
    if (i_subfr == 0)
    {
        temp1 = (T_op[0]*res) - 32;
    }
    else
    {
        temp1 = (T_op[1]*res) - 32;
    }

    if (T_op[0]<T_op[1])
    {
        t = (T_op[0]*res) - 16;
    }
    else
    {
        t = (T_op[1]*res) - 16;
    }

    if (temp1<t)
    {
        temp1 = t;
    }

    temp2 = temp1 / res;
    *T0_min = temp2;
    *T0_min_frac = temp1 - temp2*res;

    if ( *T0_min < pit_min)
    {
        *T0_min = pit_min;
        *T0_min_frac = 0;
    }

    /* Higher-bound */
    temp1 = (*T0_min*res) + *T0_min_frac + 64 - 1;

    if (T_op[0]<T_op[1])
    {
        t = (T_op[1]*res) + 16 + res - 1;
    }
    else
    {
        t = (T_op[0]*res) + 16 + res - 1;
    }

    if (temp1>t)
    {
        temp1 = t;
    }

    temp2 = temp1 / res;
    *T0_max = temp2;
    *T0_max_frac = temp1 - temp2*res;

    if ( *T0_max > pit_max)
    {
        *T0_max = pit_max;
        *T0_max_frac = res - 1;
        temp1 = (*T0_max*res) - 64 + res;
        temp2 = temp1 / res;
        *T0_min = temp2;
        *T0_min_frac = temp1 - temp2*res;
    }

    return;
}


/*------------------------------------------------------------------*
 * Mode2_pit_encode:
 *
 * Close-loop pitch lag search and pitch lag quantization
 * Adaptive excitation construction
 *------------------------------------------------------------------*/

void Mode2_pit_encode(
    short coder_type,   /* i  : coding model                               */
    short i_subfr,      /* i  : subframe index                             */
    int **pt_indice,  /* i/o: quantization indices pointer               */
    float *exc,         /* i/o: pointer to excitation signal frame         */
    const short *T_op,        /* i  : open loop pitch estimates in current frame */
    int *T0_min,      /* i/o: lower limit for close-loop search          */
    int *T0_min_frac, /* i/o: lower limit for close-loop search          */
    int *T0_max,      /* i/o: higher limit for close-loop search         */
    int *T0_max_frac, /* i/o: higher limit for close-loop search         */
    int *T0,          /* i/o: close loop integer pitch                   */
    int *T0_frac,     /* i/o: close loop fractional part of the pitch    */
    int *T0_res,      /* i/o: close loop pitch resolution                */
    float *h1,          /* i  : weighted filter impulse response           */
    float *xn,          /* i  : target vector                              */
    int pit_min,
    int pit_fr1,
    int pit_fr1b,
    int pit_fr2,
    int pit_max,
    int pit_res_max)
{
    int pit_flag;

    /* Pitch flag */
    pit_flag = i_subfr;

    if (i_subfr == (2*L_SUBFR))
    {
        pit_flag = 0;
    }

    /*-----------------------------------------------------------------*
     *  - Limit range of pitch search
     *  - Fractional pitch search
     *  - Pitch quantization
     *-----------------------------------------------------------------*/

    if(coder_type == 0) /*Unvoiced Coding do nothing*/
    {
        *T0 = L_SUBFR;
        *T0_frac = 0;
        *T0_res = 1;
    }
    else if(coder_type == 1) /* 8/4/4/4 (EVS) */
    {
        if (i_subfr == 0)
        {
            limit_T0_voiced( 4, pit_res_max>>1, T_op[0], 0, 1, T0_min, T0_min_frac, T0_max, T0_max_frac, pit_min, pit_max );
        }
        else
        {
            limit_T0_voiced( 4, pit_res_max>>1, *T0, *T0_frac, *T0_res, T0_min, T0_min_frac, T0_max, T0_max_frac, pit_min, pit_max );
        }

        *T0 = E_GAIN_closed_loop_search(exc, xn, h1, *T0_min, *T0_min_frac, *T0_max, *T0_max_frac, pit_res_max>>1, T0_frac, T0_res, pit_res_max,
                                        i_subfr, pit_min, pit_min, pit_fr1b, L_SUBFR );

        if (i_subfr == 0)
        {
            Mode2_abs_pit_enc( *T0, *T0_frac, pt_indice, pit_min, pit_fr1b, pit_min, pit_res_max );
        }
        else
        {
            Mode2_delta_pit_enc( *T0, *T0_frac, (pit_res_max>>1), *T0_min, *T0_min_frac, pt_indice );
        }

    }
    else if(coder_type == 2) /* 8/5/8/5 (EVS) */
    {

        if (i_subfr == 0)
        {
            limit_T0_voiced( 5, pit_res_max>>1, T_op[0], 0, 1, T0_min, T0_min_frac, T0_max, T0_max_frac, pit_min, pit_max );
        }
        else if(i_subfr == 2*L_SUBFR)
        {
            limit_T0_voiced( 5, pit_res_max>>1, T_op[1], 0, 1, T0_min, T0_min_frac, T0_max, T0_max_frac, pit_min, pit_max );
        }
        else
        {
            limit_T0_voiced( 5, pit_res_max>>1, *T0, *T0_frac, *T0_res, T0_min, T0_min_frac, T0_max, T0_max_frac, pit_min, pit_max );
        }

        *T0 = E_GAIN_closed_loop_search(exc, xn, h1, *T0_min, *T0_min_frac, *T0_max, *T0_max_frac, pit_res_max>>1, T0_frac, T0_res, pit_res_max,
                                        pit_flag, pit_min, pit_min, pit_fr1b, L_SUBFR);

        if (pit_flag == 0)
        {
            Mode2_abs_pit_enc( *T0, *T0_frac, pt_indice, pit_min, pit_fr1b, pit_min, pit_res_max );
        }
        else
        {
            Mode2_delta_pit_enc( *T0, *T0_frac, (pit_res_max>>1), *T0_min, *T0_min_frac, pt_indice );
        }
    }
    else if(coder_type == 3) /* 9/6/6/6 (HRs- VC) */
    {
        int pit_res_max2 = pit_res_max;

        if ( pit_min==PIT_MIN_16k )
        {

            pit_res_max2 = pit_res_max >> 1;
        }

        if (i_subfr == 0)
        {

            limit_T0_voiced2( pit_res_max2, T_op, T0_min, T0_min_frac, T0_max, T0_max_frac, pit_min, pit_max, i_subfr );
        }
        else
        {
            limit_T0_voiced( 6, pit_res_max2, *T0, 0, 1, T0_min, T0_min_frac, T0_max, T0_max_frac, pit_min, pit_max);
        }

        *T0 = E_GAIN_closed_loop_search(exc, xn, h1, *T0_min, *T0_min_frac, *T0_max, *T0_max_frac, pit_res_max2, T0_frac, T0_res, pit_res_max,
                                        i_subfr, pit_min, pit_fr2, pit_fr1, L_SUBFR);

        if (i_subfr == 0) /* if 1st subframe */
        {
            Mode2_abs_pit_enc( *T0, *T0_frac, pt_indice, pit_min, pit_fr1, pit_fr2, pit_res_max );
        }
        else
        {
            Mode2_delta_pit_enc( *T0, *T0_frac, pit_res_max2, *T0_min, *T0_min_frac, pt_indice );
        }
    }
    else if(coder_type == 4) /* 9/6/9/6 (AMRWB) */
    {
        int pit_res_max2 = pit_res_max;

        if ( pit_min==PIT_MIN_16k )
        {

            pit_res_max2 = pit_res_max >> 1;
        }


        if ( (i_subfr == 0) || (i_subfr == 2*L_SUBFR) )
        {

            limit_T0_voiced2( pit_res_max2, T_op, T0_min, T0_min_frac, T0_max, T0_max_frac, pit_min, pit_max, i_subfr );
        }
        else
        {
            limit_T0_voiced( 6, pit_res_max2, *T0, 0, 1, T0_min, T0_min_frac, T0_max, T0_max_frac, pit_min, pit_max);
        }

        *T0 = E_GAIN_closed_loop_search(exc, xn, h1, *T0_min, *T0_min_frac, *T0_max, *T0_max_frac, pit_res_max2, T0_frac, T0_res, pit_res_max,
                                        pit_flag, pit_min, pit_fr2, pit_fr1, L_SUBFR);

        if (pit_flag == 0) /* if 1st/3rd/5th subframe */
        {
            Mode2_abs_pit_enc( *T0, *T0_frac, pt_indice, pit_min, pit_fr1, pit_fr2, pit_res_max );
        }
        else /* if subframe 2 or 4 */
        {
            Mode2_delta_pit_enc( *T0, *T0_frac, pit_res_max2, *T0_min, *T0_min_frac, pt_indice );
        }
    }
    else if(coder_type == 8) /* 8/5/5/5 (RF all pred mode) */
    {
        if (i_subfr == 0)
        {
            limit_T0_voiced( 5, pit_res_max>>1, T_op[0], 0, 1, T0_min, T0_min_frac, T0_max, T0_max_frac, pit_min, pit_max );
        }
        else
        {
            limit_T0_voiced( 5, pit_res_max>>1, *T0, *T0_frac, *T0_res, T0_min, T0_min_frac, T0_max, T0_max_frac, pit_min, pit_max );
        }
        *T0 = E_GAIN_closed_loop_search(exc, xn, h1, *T0_min, *T0_min_frac, *T0_max, *T0_max_frac, pit_res_max>>1, T0_frac, T0_res, pit_res_max,
                                        i_subfr, pit_min, pit_min, pit_fr1b, L_SUBFR );

        if (i_subfr == 0)
        {
            Mode2_abs_pit_enc( *T0, *T0_frac, pt_indice, pit_min, pit_fr1b, pit_min, pit_res_max );
        }
        else
        {
            Mode2_delta_pit_enc( *T0, *T0_frac, (pit_res_max>>1), *T0_min, *T0_min_frac, pt_indice );
        }
    }
    else if(coder_type == 9) /* 8/0/8/0 (RF mode Gen pred) */
    {
        if (i_subfr == 0)
        {
            limit_T0_voiced( 4, pit_res_max>>1, T_op[0], 0, 1, T0_min, T0_min_frac, T0_max, T0_max_frac, pit_min, pit_max );
        }
        else
        {
            limit_T0_voiced( 4, pit_res_max>>1, *T0, *T0_frac, *T0_res, T0_min, T0_min_frac, T0_max, T0_max_frac, pit_min, pit_max );
        }
        *T0 = E_GAIN_closed_loop_search(exc, xn, h1, *T0_min, *T0_min_frac, *T0_max, *T0_max_frac, pit_res_max>>1, T0_frac, T0_res, pit_res_max,
                                        i_subfr, pit_min, pit_min, pit_fr1b, L_SUBFR );

        if (i_subfr == 0)
        {
            Mode2_abs_pit_enc( *T0, *T0_frac, pt_indice, pit_min, pit_fr1b, pit_min, pit_res_max );
        }
        else
        {
            Mode2_delta_pit_enc( *T0, *T0_frac, (pit_res_max>>1), *T0_min, *T0_min_frac, pt_indice );
        }
    }
    else
    {
        assert(0);
    }

    return;
}


/*-------------------------------------------------------------------*
 * Mode2_abs_pit_enc:
 *
 * Encode pitch lag absolutely
 *-------------------------------------------------------------------*/

void Mode2_abs_pit_enc(
    short T0,          /* i  : integer pitch lag              */
    int T0_frac,     /* i  : pitch fraction                 */
    int **pt_indice, /* i/o: pointer to Vector of Q indexes */
    short pit_min,
    short pit_fr1,
    short pit_fr2,
    short pit_res_max
)
{
    short pit_res_max_half;

    pit_res_max_half= pit_res_max>>1;

    if (T0 < pit_fr2)
    {
        **pt_indice = T0*pit_res_max + T0_frac - (pit_min*pit_res_max);
    }
    else if (T0 < pit_fr1)
    {
        **pt_indice = T0*pit_res_max_half + T0_frac - (pit_fr2*pit_res_max_half) + ((pit_fr2-pit_min)*pit_res_max);

    }
    else
    {
        **pt_indice = T0 - pit_fr1 + ((pit_fr2-pit_min)*pit_res_max) + ((pit_fr1-pit_fr2)*pit_res_max_half);
    }

    (*pt_indice)++;

    return;
}


/*-------------------------------------------------------------------*
 * Mode2_delta_pit_enc:
 *
 * Encode pitch lag differentially
 *-------------------------------------------------------------------*/

void Mode2_delta_pit_enc(
    short T0,          /* i  : integer pitch lag              */
    int T0_frac,     /* i  : pitch fraction                 */
    int T0_res,      /* i  : pitch resolution               */
    short T0_min,      /* i/o: delta search min               */
    short T0_min_frac, /* i/o: delta search min               */
    int **pt_indice  /* i/o: pointer to Vector of Q indexes */
)
{

    **pt_indice = (T0 - T0_min) * T0_res + T0_frac - T0_min_frac;

    (*pt_indice)++;

    return;
}
