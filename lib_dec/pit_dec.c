/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"
#include "options.h"


/*----------------------------------------------------------*
 * pit_decode()
 *
 * Decode OL pitch lag
 *----------------------------------------------------------*/

float pit_decode(                 /* o  : floating pitch value                    */
    Decoder_State *st,              /* i/o: decoder state structure                 */
    const long  core_brate,       /* i  : core bitrate                            */
    const short Opt_AMR_WB,       /* i  : flag indicating AMR-WB IO mode          */
    const short L_frame,          /* i  : length of the frame                     */
    short i_subfr,          /* i  : subframe index                          */
    const short coder_type,       /* i  : coding type                             */
    short *limit_flag,      /* i/o: restrained(0) or extended(1) Q limits   */
    short *T0,              /* o  : close loop integer pitch                */
    short *T0_frac,         /* o  : close loop fractional part of the pitch */
    short *T0_min,          /* i/o: delta search min                        */
    short *T0_max,          /* i/o: delta search max                        */
    const short L_subfr           /* i  : subframe length                         */
)
{
    float pitch;
    short pitch_index, nBits, pit_flag;

    pitch_index = 0;

    /*----------------------------------------------------------------*
     * Set pit_flag = 0 for every subframe with absolute pitch search
     *----------------------------------------------------------------*/

    pit_flag = i_subfr;
    if( i_subfr == 2*L_SUBFR )
    {
        pit_flag = 0;
    }

    /*-------------------------------------------------------*
     * Retrieve the pitch index
     *-------------------------------------------------------*/

    if( !Opt_AMR_WB )
    {
        /*----------------------------------------------------------------*
         * pitch Q: Set limit_flag to 0 for restrained limits, and 1 for extended limits
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

            pitch_index = (short)get_next_indice( st, nBits );
        }

        /*-------------------------------------------------------*
         * Pitch decoding in AUDIO mode
         * (both ACELP@12k8 and ACELP@16k cores)
         *-------------------------------------------------------*/

        if( coder_type == AUDIO )
        {
            if( L_subfr == L_FRAME/2 && i_subfr != 0 )
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

            pitch_index = (short)get_next_indice( st, nBits );

            if( L_subfr == L_FRAME/2 && i_subfr != 0 && pitch_index >= 32 )  /* safety check in case of bit errors */
            {
                pitch_index = pitch_index>>1;
                st->BER_detect = 1;
            }

            pit_Q_dec( 0, pitch_index, nBits, 4, pit_flag, *limit_flag, T0, T0_frac, T0_min, T0_max, &st->BER_detect );
        }
        else if( coder_type == VOICED )
        {
            /*-------------------------------------------------------*
             * Pitch decoding in VOICED mode
             * (ACELP@12k8 core only)
             *-------------------------------------------------------*/

            if( i_subfr == 2*L_SUBFR )
            {
                pit_flag = i_subfr;
            }

            pit_Q_dec( 0, pitch_index, nBits, 4, pit_flag, *limit_flag, T0, T0_frac, T0_min, T0_max, &st->BER_detect );
        }
        else
        {
            /*-------------------------------------------------------*
             *  Pitch decoding in GENERIC mode
             * (both ACELP@12k8 and ACELP@16k cores)
             *-------------------------------------------------------*/

            if( L_frame == L_FRAME )
            {
                pit_Q_dec( 0, pitch_index, nBits, 8, pit_flag, *limit_flag, T0, T0_frac, T0_min, T0_max, &st->BER_detect );
            }
            else
            {
                pit16k_Q_dec( pitch_index, nBits, *limit_flag, T0, T0_frac, T0_min, T0_max, &st->BER_detect );
            }
        }
    }

    /*-------------------------------------------------------*
     *  Pitch decoding in AMR-WB IO mode
     *-------------------------------------------------------*/

    else
    {
        *limit_flag = 0;

        if( i_subfr == 0 || ( i_subfr == 2*L_SUBFR && core_brate == ACELP_8k85 ) )
        {
            nBits = 8;
        }
        else
        {
            nBits = 5;
        }

        if( core_brate > ACELP_8k85 )
        {
            nBits = 6;

            if( i_subfr == 0 || i_subfr == 2*L_SUBFR )
            {
                nBits = 9;
            }
        }

        pitch_index = (short)get_next_indice( st, nBits );

        pit_Q_dec( 1, pitch_index, nBits, 8, pit_flag, *limit_flag, T0, T0_frac, T0_min, T0_max, &st->BER_detect );
    }

    /*-------------------------------------------------------*
     * Compute floating pitch output
     *-------------------------------------------------------*/

    pitch = (float)(*T0) + (float)(*T0_frac)/4.0f;

    return pitch;
}


/*---------------------------------------------------------------------*
 * abs_pit_dec()
 *
 * Decode the absolute pitch
 *---------------------------------------------------------------------*/

void abs_pit_dec(
    const short fr_steps,    /* i:   fractional resolution steps (0, 2, 4)    */
    short pitch_index, /* i:   pitch index                              */
    const short limit_flag,   /* i  : restrained(0) or extended(1) limits     */
    short *T0,         /* o:   integer pitch lag                        */
    short *T0_frac     /* o:   pitch fraction                           */
)
{
    if( limit_flag == 0 )
    {
        if( fr_steps == 2 )
        {
            if( pitch_index < (PIT_FR1_8b-PIT_MIN)*2 )
            {
                *T0 = PIT_MIN + (pitch_index/2);
                *T0_frac = pitch_index - ((*T0 - PIT_MIN)*2);
                *T0_frac *= 2;
            }
            else
            {
                *T0 = pitch_index + PIT_FR1_8b - ((PIT_FR1_8b-PIT_MIN)*2);
                *T0_frac = 0;
            }
        }
        else if( fr_steps == 4 )
        {
            if( pitch_index < (PIT_FR2_9b-PIT_MIN)*4 )
            {
                *T0 = PIT_MIN + (pitch_index/4);
                *T0_frac = pitch_index - (*T0 - PIT_MIN)*4;
            }
            else if( pitch_index < ( (PIT_FR2_9b-PIT_MIN)*4 + (PIT_FR1_9b-PIT_FR2_9b)*2) )
            {
                pitch_index -=  (PIT_FR2_9b-PIT_MIN)*4;
                *T0 = PIT_FR2_9b + (pitch_index/2);
                *T0_frac = pitch_index - (*T0 - PIT_FR2_9b)*2;
                (*T0_frac) *= 2;
            }
            else
            {
                *T0 = pitch_index + PIT_FR1_9b - ((PIT_FR2_9b-PIT_MIN)*4) - ((PIT_FR1_9b-PIT_FR2_9b)*2);
                *T0_frac = 0;
            }
        }
        else  /* fr_steps == 0 */
        {
            /* not used in the codec */
        }
    }
    else if( limit_flag == 1 )   /* extended Q range */
    {
        if( fr_steps == 2 )
        {
            if( pitch_index < (PIT_FR1_EXTEND_8b-PIT_MIN_EXTEND)*2 )
            {
                *T0 = PIT_MIN_EXTEND + (pitch_index/2);
                *T0_frac = pitch_index - ((*T0 - PIT_MIN_EXTEND)*2);
                *T0_frac *= 2;
            }
            else
            {
                *T0 = pitch_index + PIT_FR1_EXTEND_8b - ((PIT_FR1_EXTEND_8b-PIT_MIN_EXTEND)*2);
                *T0_frac = 0;
            }
        }
        else if( fr_steps == 4 )
        {
            if( pitch_index < (PIT_FR2_EXTEND_9b-PIT_MIN_EXTEND)*4 )
            {
                *T0 = PIT_MIN_EXTEND + (pitch_index/4);
                *T0_frac = pitch_index - (*T0 - PIT_MIN_EXTEND)*4;
            }
            else if( pitch_index < ( (PIT_FR2_EXTEND_9b-PIT_MIN_EXTEND)*4 + (PIT_FR1_EXTEND_9b-PIT_FR2_EXTEND_9b)*2) )
            {
                pitch_index -=  (PIT_FR2_EXTEND_9b-PIT_MIN_EXTEND)*4;
                *T0 = PIT_FR2_EXTEND_9b + (pitch_index/2);
                *T0_frac = pitch_index - (*T0 - PIT_FR2_EXTEND_9b)*2;
                (*T0_frac) *= 2;
            }
            else
            {
                *T0 = pitch_index + PIT_FR1_EXTEND_9b - ((PIT_FR2_EXTEND_9b-PIT_MIN_EXTEND)*4) - ((PIT_FR1_EXTEND_9b-PIT_FR2_EXTEND_9b)*2);
                *T0_frac = 0;
            }
        }
        else  /* fr_steps == 0 */
        {
            /* not used in the codec */
        }
    }
    else  /* limit_flag == 2 */
    {
        if( fr_steps == 2 )
        {
            if( pitch_index < (PIT_FR1_DOUBLEEXTEND_8b-PIT_MIN_DOUBLEEXTEND)*2 )
            {
                *T0 = PIT_MIN_DOUBLEEXTEND + (pitch_index/2);
                *T0_frac = pitch_index - ((*T0 - PIT_MIN_DOUBLEEXTEND)*2);
                *T0_frac *= 2;
            }
            else
            {
                *T0 = pitch_index + PIT_FR1_DOUBLEEXTEND_8b - ((PIT_FR1_DOUBLEEXTEND_8b-PIT_MIN_DOUBLEEXTEND)*2);
                *T0_frac = 0;
            }
        }
        else if( fr_steps == 4 )
        {
            if( pitch_index < (PIT_FR2_DOUBLEEXTEND_9b-PIT_MIN_DOUBLEEXTEND)*4 )
            {
                *T0 = PIT_MIN_DOUBLEEXTEND + (pitch_index/4);
                *T0_frac = pitch_index - (*T0 - PIT_MIN_DOUBLEEXTEND)*4;
            }
            else if( pitch_index < ( (PIT_FR2_DOUBLEEXTEND_9b-PIT_MIN_DOUBLEEXTEND)*4 + (PIT_FR1_DOUBLEEXTEND_9b-PIT_FR2_DOUBLEEXTEND_9b)*2) )
            {
                pitch_index -=  (PIT_FR2_DOUBLEEXTEND_9b-PIT_MIN_DOUBLEEXTEND)*4;
                *T0 = PIT_FR2_DOUBLEEXTEND_9b + (pitch_index/2);
                *T0_frac = pitch_index - (*T0 - PIT_FR2_DOUBLEEXTEND_9b)*2;
                (*T0_frac) *= 2;
            }
            else
            {
                *T0 = pitch_index + PIT_FR1_DOUBLEEXTEND_9b - ((PIT_FR2_DOUBLEEXTEND_9b-PIT_MIN_DOUBLEEXTEND)*4) - ((PIT_FR1_DOUBLEEXTEND_9b-PIT_FR2_DOUBLEEXTEND_9b)*2);
                *T0_frac = 0;
            }
        }
        else  /* fr_steps == 0 */
        {
            /* not used in the codec */
        }
    }

    return;
}


/*---------------------------------------------------------------------*
 * delta_pit_dec()
 *
 * Decode delta pitch
 *---------------------------------------------------------------------*/

void delta_pit_dec(
    const short fr_steps,    /* i  : fractional resolution steps (0, 2, 4)    */
    const short pitch_index, /* i  : pitch index                              */
    short *T0,         /* o  : integer pitch lag                        */
    short *T0_frac,    /* o  : pitch fraction                           */
    const short T0_min       /* i  : delta search min                         */
)
{
    if( fr_steps == 0 )
    {
        *T0 = T0_min + pitch_index;
        *T0_frac = 0;
    }
    else if( fr_steps == 2 )
    {
        *T0 = T0_min + pitch_index/2;
        *T0_frac = pitch_index - (*T0 - T0_min)*2;
        *T0_frac *= 2;
    }
    else if( fr_steps == 4 )
    {
        *T0 = T0_min + pitch_index/4;
        *T0_frac = pitch_index - (*T0 - T0_min)*4;
    }

    return;
}


/*-------------------------------------------------*
 * pit_Q_dec()
 *
 * pitch decoding
 *-------------------------------------------------*/

void pit_Q_dec(
    const short Opt_AMR_WB,   /* i  : flag indicating AMR-WB IO mode          */
    const short pitch_index,  /* i  : pitch index                             */
    const short nBits,        /* i  : # of Q bits                             */
    const short delta,        /* i  : Half the CL searched interval           */
    const short pit_flag,     /* i  : absolute(0) or delta(1) pitch Q         */
    const short limit_flag,   /* i  : restrained(0) or extended(1) Q limits   */
    short *T0,          /* o  : integer pitch lag                       */
    short *T0_frac,     /* o  : pitch fraction                          */
    short *T0_min,      /* i/o: delta search min                        */
    short *T0_max       /* i/o: delta search max                        */
    ,short *BER_detect   /* o  : BER detect flag                         */
)
{
    if( nBits == 10 )         /* absolute decoding with 10 bits */
    {
        if( limit_flag == 0 )
        {
            *T0 = PIT_MIN + (pitch_index/4);
            *T0_frac = pitch_index - (*T0 - PIT_MIN)*4;
        }
        else if( limit_flag == 1 )
        {
            *T0 = PIT_MIN_EXTEND + (pitch_index/4);
            *T0_frac = pitch_index - (*T0 - PIT_MIN_EXTEND)*4;
        }
        else  /* limit_flag == 2 */
        {
            *T0 = PIT_MIN_DOUBLEEXTEND + (pitch_index/4);
            *T0_frac = pitch_index - (*T0 - PIT_MIN_DOUBLEEXTEND)*4;
        }
    }
    else if( nBits == 9 )     /* absolute decoding with 9 bits */
    {
        abs_pit_dec( 4, pitch_index, limit_flag, T0, T0_frac );

        /* find T0_min and T0_max for delta search */
        if( Opt_AMR_WB )
        {
            limit_T0( L_FRAME, delta, pit_flag, 0, *T0, 0, T0_min, T0_max );        /* T0_frac==0 to keep IO with AMR-WB */
        }
    }
    else if( nBits == 8 )     /* absolute decoding with 8 bits */
    {
        abs_pit_dec( 2, pitch_index, limit_flag, T0, T0_frac );

        /* find T0_min and T0_max for delta search */
        if( Opt_AMR_WB )
        {
            limit_T0( L_FRAME, delta, pit_flag, 0, *T0, 0, T0_min, T0_max );        /* T0_frac==0 to keep IO with AMR-WB */
        }
    }
    else if( nBits == 6 )     /* relative decoding with 6 bits */
    {
        delta_pit_dec( 4, pitch_index, T0, T0_frac, *T0_min );
    }
    else if( nBits == 5 )     /* relative decoding with 5 bits */
    {
        if( delta == 8 )
        {
            delta_pit_dec( 2, pitch_index, T0, T0_frac, *T0_min );
        }
        else  /* delta == 4 */
        {
            delta_pit_dec( 4, pitch_index, T0, T0_frac, *T0_min );
        }
    }
    else /* nBits == 4 */     /* relative decoding with 4 bits */
    {
        if( delta == 8 )
        {
            delta_pit_dec( 0, pitch_index, T0, T0_frac, *T0_min );
        }
        else  /* delta == 4 */
        {
            delta_pit_dec( 2, pitch_index, T0, T0_frac, *T0_min );
        }
    }

    /* biterror detection mechanism */
    if( ((*T0<<2) + *T0_frac) > (PIT_MAX<<2)+2 && pit_flag == 0 && !Opt_AMR_WB )
    {
        *T0 = L_SUBFR;
        *T0_frac = 0;
        *BER_detect = 1;
    }

    if( !Opt_AMR_WB )
    {
        /* find T0_min and T0_max for delta search */
        limit_T0( L_FRAME, delta, L_SUBFR, limit_flag, *T0, *T0_frac, T0_min, T0_max );
    }

    return;
}

/*-------------------------------------------------*
  * pit16k_Q_dec()
  *
  * pitch decoding @16kHz core
  *-------------------------------------------------*/

void pit16k_Q_dec(
    const short pitch_index,  /* i  : pitch index                             */
    const short nBits,        /* i  : # of Q bits                             */
    const short limit_flag,   /* i  : restrained(0) or extended(1) limits     */
    short *T0,          /* o  : integer pitch lag                       */
    short *T0_frac,     /* o  : pitch fraction                          */
    short *T0_min,      /* i/o: delta search min                        */
    short *T0_max       /* i/o: delta search max                        */
    ,short *BER_detect   /* o  : BER detect flag                         */
)
{
    short index;

    if( nBits == 10 )         /* absolute decoding with 10 bits */
    {
        {
            if( pitch_index < (PIT16k_FR2_EXTEND_10b-PIT16k_MIN_EXTEND)*4 )
            {
                *T0 = PIT16k_MIN_EXTEND + (pitch_index/4);
                *T0_frac = pitch_index - ((*T0 - PIT16k_MIN_EXTEND)*4);
            }
            else
            {
                index =  pitch_index - (PIT16k_FR2_EXTEND_10b-PIT16k_MIN_EXTEND)*4;
                *T0 = PIT16k_FR2_EXTEND_10b + (index/2);
                *T0_frac = index - (*T0 - PIT16k_FR2_EXTEND_10b)*2;
                (*T0_frac) *= 2;
            }
        }

    }
    else if( nBits == 9 )     /* absolute decoding with 9 bits */
    {
        {
            if (pitch_index < (PIT16k_FR2_EXTEND_9b-PIT16k_MIN_EXTEND)*4)
            {
                *T0 = PIT16k_MIN_EXTEND + (pitch_index/4);
                *T0_frac = pitch_index - (*T0 - PIT16k_MIN_EXTEND)*4;
            }
            else if (pitch_index < ( (PIT16k_FR2_EXTEND_9b-PIT16k_MIN_EXTEND)*4 + (PIT16k_FR1_EXTEND_9b-PIT16k_FR2_EXTEND_9b)*2) )
            {
                index =  pitch_index - (PIT16k_FR2_EXTEND_9b-PIT16k_MIN_EXTEND)*4;
                *T0 = PIT16k_FR2_EXTEND_9b + (index/2);
                *T0_frac = index - (*T0 - PIT16k_FR2_EXTEND_9b)*2;
                (*T0_frac) *= 2;
            }
            else
            {
                *T0 = pitch_index + PIT16k_FR1_EXTEND_9b - ((PIT16k_FR2_EXTEND_9b-PIT16k_MIN_EXTEND)*4) - ((PIT16k_FR1_EXTEND_9b-PIT16k_FR2_EXTEND_9b)*2);
                *T0_frac = 0;
            }
        }
    }
    else  /* nBits == 6 */    /* relative decoding with 6 bits */
    {
        delta_pit_dec( 4, pitch_index, T0, T0_frac, *T0_min );
    }

    /* biterror detection mechanism */
    if( ((*T0<<2) + *T0_frac) > (PIT16k_MAX<<2) && nBits >= 9 )
    {
        *T0 = L_SUBFR;
        *T0_frac = 0;
        *BER_detect = 1;
    }

    /* find T0_min and T0_max for delta search */
    limit_T0( L_FRAME16k, 8, L_SUBFR, limit_flag, *T0, *T0_frac, T0_min, T0_max );

    return;
}


/*----------------------------------------------------------*
 * Mode2_pit_decode
 *
 * Decode pitch lag
 *----------------------------------------------------------*/

float Mode2_pit_decode(       /* o:   floating pitch value                      */
    const short coder_type,   /* i:   coding model                              */
    short i_subfr,      /* i:   subframe index                            */
    int L_subfr,
    int **pt_indice,  /* i/o: quantization indices pointer              */
    int *T0,          /* i/o:   close loop integer pitch                */
    int *T0_frac,     /* o:   close loop fractional part of the pitch   */
    int *T0_res,      /* i/o: pitch resolution                          */
    int *T0_min,      /* i/o: lower limit for close-loop search         */
    int *T0_min_frac, /* i/o: lower limit for close-loop search         */
    int *T0_max,      /* i/o: higher limit for close-loop search        */
    int *T0_max_frac, /* i/o: higher limit for close-loop search        */
    int  pit_min,
    int  pit_fr1,
    int  pit_fr1b,
    int  pit_fr2,
    int  pit_max,
    int  pit_res_max
)
{
    float pitch;

    if( coder_type == 0 ) /*Unvoiced Coding do nothing*/
    {

        *T0 = L_subfr;
        *T0_frac = 0;
        *T0_res = 1;
    }
    else if(coder_type == 1) /* 8/4/4/4 (EVS) */
    {

        if (i_subfr == 0)
        {
            Mode2_abs_pit_dec( T0, T0_frac, T0_res, pt_indice, pit_min, pit_fr1b, pit_min, pit_res_max);
        }
        else
        {
            limit_T0_voiced( 4, pit_res_max>>1, *T0, *T0_frac, *T0_res, T0_min, T0_min_frac, T0_max, T0_max_frac, pit_min, pit_max );
            *T0_res = (pit_res_max>>1);
            Mode2_delta_pit_dec( T0, T0_frac, *T0_res, T0_min, T0_min_frac, pt_indice);
        }
    }
    else if(coder_type == 2) /* 8/5/8/5 (EVS) */
    {

        if ( ( i_subfr == 0 ) || ( i_subfr == 2*L_subfr ) )
        {
            Mode2_abs_pit_dec( T0, T0_frac, T0_res, pt_indice, pit_min, pit_fr1b, pit_min, pit_res_max);
        }
        else
        {
            limit_T0_voiced( 5, pit_res_max>>1, *T0, *T0_frac, *T0_res, T0_min, T0_min_frac, T0_max, T0_max_frac, pit_min, pit_max );
            *T0_res = (pit_res_max>>1);
            Mode2_delta_pit_dec( T0, T0_frac, *T0_res, T0_min, T0_min_frac, pt_indice);
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
            Mode2_abs_pit_dec( T0, T0_frac, T0_res, pt_indice, pit_min, pit_fr1, pit_fr2, pit_res_max);
        }
        else
        {
            limit_T0_voiced( 6, pit_res_max2, *T0, 0, 1, T0_min, T0_min_frac, T0_max, T0_max_frac, pit_min, pit_max );
            *T0_res = pit_res_max2;
            Mode2_delta_pit_dec( T0, T0_frac, *T0_res, T0_min, T0_min_frac, pt_indice);
        }
    }
    else if(coder_type == 4) /* 9/6/9/6 (AMRWB) */
    {
        int pit_res_max2 = pit_res_max;
        if ( pit_min==PIT_MIN_16k )
        {
            pit_res_max2 = pit_res_max >> 1;
        }

        if ( ( i_subfr == 0 ) || ( i_subfr == 2*L_subfr ) )
        {
            Mode2_abs_pit_dec( T0, T0_frac, T0_res, pt_indice, pit_min, pit_fr1, pit_fr2, pit_res_max);
        }
        else
        {
            limit_T0_voiced( 6, pit_res_max2, *T0, 0, 1, T0_min, T0_min_frac, T0_max, T0_max_frac, pit_min, pit_max );
            *T0_res = pit_res_max2;
            Mode2_delta_pit_dec( T0, T0_frac, *T0_res, T0_min, T0_min_frac, pt_indice);
        }
    }
    else if(coder_type == 8) /* 8/5/5/5 (RF all pred mode) */
    {

        if (i_subfr == 0)
        {
            Mode2_abs_pit_dec( T0, T0_frac, T0_res, pt_indice, pit_min, pit_fr1b, pit_min, pit_res_max);
        }
        else
        {
            limit_T0_voiced( 5, pit_res_max>>1, *T0, *T0_frac, *T0_res, T0_min, T0_min_frac, T0_max, T0_max_frac, pit_min, pit_max );
            *T0_res = (pit_res_max>>1);
            Mode2_delta_pit_dec( T0, T0_frac, *T0_res, T0_min, T0_min_frac, pt_indice);
        }
    }
    else if(coder_type == 9) /* 8/0/8/0 (RF gen pred mode) */
    {

        if (i_subfr == 0)
        {
            Mode2_abs_pit_dec( T0, T0_frac, T0_res, pt_indice, pit_min, pit_fr1b, pit_min, pit_res_max);
        }
        else
        {
            limit_T0_voiced( 4, pit_res_max>>1, *T0, *T0_frac, *T0_res, T0_min, T0_min_frac, T0_max, T0_max_frac, pit_min, pit_max );
            *T0_res = (pit_res_max>>1);
            Mode2_delta_pit_dec( T0, T0_frac, *T0_res, T0_min, T0_min_frac, pt_indice);
        }
    }
    else /*RCELP 8 bits pitch delay*/
    {
        assert(0);
    }

    /*-------------------------------------------------------*
     * Compute floating pitch output
     *-------------------------------------------------------*/

    pitch = (float)(*T0) + (float)(*T0_frac)/(float)(*T0_res);   /* save subframe pitch values  */

    return pitch;
}


/*---------------------------------------------------------------------*
 * Mode2_abs_pit_dec
 *
 * Decode the absolute pitch
 *---------------------------------------------------------------------*/
void Mode2_abs_pit_dec(
    int *T0,         /* o:   integer pitch lag              */
    int *T0_frac,    /* o:   pitch fraction                 */
    int *T0_res,     /* o:   pitch resolution               */
    int **pt_indice, /* i/o: pointer to Vector of Q indexes */
    int pit_min,
    int pit_fr1,
    int pit_fr2,
    int pit_res_max
)
{
    int index;
    int pit_res_max_half;

    index = **pt_indice;
    (*pt_indice)++;
    pit_res_max_half = pit_res_max>>1;

    if (index < (pit_fr2-pit_min)*pit_res_max)
    {

        *T0 = pit_min + (index/pit_res_max);

        *T0_frac = index - (*T0 - pit_min)*pit_res_max;
        *T0_res = pit_res_max;
    }
    else if (index < ( (pit_fr2-pit_min)*pit_res_max + (pit_fr1-pit_fr2)*pit_res_max_half) )
    {
        index -=  (pit_fr2-pit_min)*pit_res_max;

        *T0 = pit_fr2 + (index/pit_res_max_half);

        *T0_frac = index - (*T0 - pit_fr2)*pit_res_max_half;
        *T0_res = pit_res_max_half;
    }
    else
    {
        *T0 = index + pit_fr1 - ((pit_fr2-pit_min)*pit_res_max) - ((pit_fr1-pit_fr2)*pit_res_max_half);
        *T0_frac = 0;
        *T0_res = 1;
    }

    return;
}


/*---------------------------------------------------------------------*
 * Routine Mode2_delta_pit_dec()
 *
 * Decode delta pitch
 *---------------------------------------------------------------------*/
void Mode2_delta_pit_dec(
    int       *T0,          /* o:   integer pitch lag              */
    int       *T0_frac,     /* o:   pitch fraction                 */
    int       T0_res,       /* i:   pitch resolution               */
    int       *T0_min,      /* i/o: delta search min               */
    int       *T0_min_frac, /* i: delta search min                 */
    int       **pt_indice   /* i/o: pointer to Vector of Q indexes */
)
{
    int index;

    index = **pt_indice;
    (*pt_indice)++;

    *T0 = *T0_min + (index+*T0_min_frac)/T0_res;

    *T0_frac = index+*T0_min_frac - (*T0 - *T0_min)*T0_res;

    return;
}
