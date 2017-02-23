/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "prot.h"

/*------------------------------------------------------------------
 * weight_a()
 *
 * Weighting of LP filter coefficients, ap[i] = a[i] * (gamma^i)
 *------------------------------------------------------------------*/

void weight_a(
    const float *a,     /* i  : LP filter coefficients          */
    float *ap,    /* o  : weighted LP filter coefficients */
    const float gamma,  /* i  : weighting factor                */
    const short m       /* i  : order of LP filter              */
)
{
    float f;
    short i;

    ap[0] = a[0];
    f = gamma;

    for (i = 1; i <= m; i++)
    {
        ap[i] = f*a[i];
        f *= gamma;
    }

    return;
}



/*------------------------------------------------------------------
 * weight_a_subfr()
 *
 * Weighting of LP filter coefficients for multiple subframes,
 * ap[i] = a[i] * (gamma^i)
 *------------------------------------------------------------------*/

void weight_a_subfr(
    const short nb_subfr,   /* i  : number of subframes             */
    const float *A,         /* i  : LP filter coefficients          */
    float *Aw,        /* o  : weighted LP filter coefficients */
    const float gamma,      /* i  : weighting factor                */
    const short m           /* i  : order of LP filter              */
)
{
    short i, j;
    float tmp;

    for( j=0; j<nb_subfr; j++ )
    {
        Aw[j*(m+1)] = A[j*(m+1)];
    }

    /* Smoothing aka spreading aka masking envelope generation */
    tmp = gamma;
    for( i = 1; i < m+1; i++ )
    {
        for( j=0; j<nb_subfr; j++ )
        {
            Aw[i+j*(m+1)] = A[i+j*(m+1)] * tmp;
        }
        tmp *= gamma;
    }

    return;
}
