/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "prot.h"

/*-------------------------------------------------*
 * Local constants
 *-------------------------------------------------*/

#define LIMIT_PIT_REL_LOWER       2       /* delta interval to extend pitch coding in relative Q */
#define LIMIT_PIT_REL_UPPER       0

/*-------------------------------------------------*
 * limit_T0()
 *
 * Close-loop pitch lag search limitation
 *-------------------------------------------------*/

void limit_T0(
    const short L_frame,      /* i  : length of the frame                                  */
    const short delta,        /* i  : Half the close-loop searched interval                */
    const short pit_flag,     /* i  : selecting absolute(0) or delta(1) pitch quantization */
    const short limit_flag,   /* i  : flag for Q limits (0=restrained, 1=extended)         */
    const short T0,           /* i  : rough pitch estimate around which the search is done */
    const short T0_frac,      /* i  : pitch estimate fractional part                       */
    short *T0_min,      /* o  : lower pitch limit                                    */
    short *T0_max       /* o  : higher pitch limit                                   */
)
{
    short delta2, T1;
    short pit_min, pit_max;

    if( limit_flag == 0 )  /* restrained Q limits */
    {
        /* set limits */
        if( L_frame == L_FRAME )
        {
            pit_max = PIT_MAX;
            pit_min = PIT_MIN;
        }
        else /* L_frame == L_FRAME16k */
        {
            pit_max = PIT16k_MAX;
            pit_min = PIT16k_MIN;
        }

        delta2 = 2 * delta - 1;

        T1 = T0;
        if( T0_frac >= 2 )
        {
            T1++;
        }
        *T0_min = T1 - delta;

        if( *T0_min < pit_min )
        {
            *T0_min = pit_min;
        }
        *T0_max = *T0_min + delta2;

        if( *T0_max > pit_max )
        {
            *T0_max = pit_max;
            *T0_min = *T0_max - delta2;
        }
    }
    else  /* extended Q limits */
    {

        /* set limits */
        if( L_frame == L_FRAME )
        {
            pit_max = PIT_MAX;
            pit_min = PIT_MIN_EXTEND;

            if( limit_flag == 2 )
            {
                pit_min = PIT_MIN_DOUBLEEXTEND;
            }
        }
        else /* L_frame == L_FRAME16k */
        {
            pit_max = PIT16k_MAX;
            pit_min = PIT16k_MIN_EXTEND;
        }

        delta2 = 2 * delta - 1;

        T1 = T0;
        if( T0_frac >= 2 )
        {
            T1++;
        }
        *T0_min = T1 - delta;

        if( pit_flag == 0 )
        {
            /* subframes with absolute search: keep Q range */
            if( *T0_min < pit_min )
            {
                *T0_min = pit_min;
            }
            *T0_max = *T0_min + delta2;

            if( *T0_max > pit_max )
            {
                *T0_max = pit_max;
                *T0_min = *T0_max - delta2;
            }
        }
        else
        {
            /* subframes with relative search: extend Q range */
            if( *T0_min < pit_min - LIMIT_PIT_REL_LOWER )
            {
                *T0_min = pit_min - LIMIT_PIT_REL_LOWER;
            }

            if( *T0_min < L_INTERPOL )
            {
                *T0_min = L_INTERPOL ;
            }
            *T0_max = *T0_min + delta2;

            if( *T0_max > pit_max + LIMIT_PIT_REL_UPPER )
            {
                *T0_max = pit_max + LIMIT_PIT_REL_UPPER;
                *T0_min = *T0_max - delta2;
            }
        }
    }

    return;
}



/*-------------------------------------------------*
* Routine limit_T0_voiced()
*
* Close-loop pitch lag search limitation
*-------------------------------------------------*/

void limit_T0_voiced(
    int nbits,
    int res,
    int T0,            /* i  : rough pitch estimate around which the search is done */
    int T0_frac,       /* i  : pitch estimate fractional part                       */
    int T0_res,        /* i  : pitch resolution                                     */
    int *T0_min,       /* o  : lower pitch limit                                    */
    int *T0_min_frac,  /* o  : lower pitch limit                                    */
    int *T0_max,       /* o  : higher pitch limit                                   */
    int *T0_max_frac,  /* o  : higher pitch limit                                   */
    int pit_min,       /* i  : Minimum pitch lag                                    */
    int pit_max        /* i  : Maximum pitch lag                                    */
)
{
    short T1, temp1, temp2;


    /* Mid-point */
    T1 = T0;
    if( (T0_res > 1) && (T0_frac >= (T0_res>>1)) )
    {
        T1++;
    }

    /* Lower-bound */
    temp1 = (T1*res) - (1<<(nbits-1));
    temp2 = temp1 / res;
    *T0_min = temp2;
    *T0_min_frac = temp1 - temp2*res;
    if ( *T0_min < pit_min)
    {
        *T0_min = pit_min;
        *T0_min_frac = 0;
    }

    /* Higher-bound */
    temp1 = (*T0_min*res) + *T0_min_frac + (1<<nbits) - 1;
    temp2 = temp1 / res;
    *T0_max = temp2;
    *T0_max_frac = temp1 - temp2*res;
    if ( *T0_max > pit_max)
    {
        *T0_max = pit_max;
        *T0_max_frac = res - 1;
        temp1 = (*T0_max*res) + *T0_max_frac - (1<<nbits) + 1;
        temp2 = temp1 / res;
        *T0_min = temp2;
        *T0_min_frac = temp1 - temp2*res;
    }

}


