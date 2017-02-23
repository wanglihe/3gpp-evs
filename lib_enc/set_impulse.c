/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "cnst.h"
#include "rom_com.h"
#include "prot.h"

/*-----------------------------------------------------------------*
 * Local constant
 *-----------------------------------------------------------------*/
#define INPOL 4 /* +- range in samples for impulse position searching */

/*-----------------------------------------------------------------*
 * Local function prototype
 *-----------------------------------------------------------------*/
static void convolve_tc(const float g[], const float h[], float y[], const short L_1, const short L_2);
static void correlate_tc(const float *x, float *y, const float *h, const short start, const short L_1,
                         const short L_2);
static void convolve_tc2(const float g[], const float h[], float y[], const short pos_max);

/*---------------------------------------------------------------------------------------*
 * Function  set_impulse() for TC                                                        *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                                                        *
 * Builds glottal codebook contribution based on glotal impulses positions finding.      *
 *                                                                                       *
 * Returns a position of the glotal impulse center and                                   *
 * a number of the glotal impulse shape.                                                 *
 *                                                                                       *
 *               |----|              |----|                                   xn         *
 *     imp_pos-> ||   |  imp_shape-> | g1 |                                   |          *
 *               | |  |              | g2 |  exc    |---|    y1  ----         |          *
 *               |  | |--------------|    |---------| h |-------|gain|-------(-)---> xn2 *
 *               |   ||              | gn |         |---|        ----                    *
 *               |----|              |----|                                              *
 *              codebook           excitation       h_orig       gain                    *
 *                                                                                       *
 *                                                                                       *
 *                                 nominator      dd    <xn,y1>*<xn,y1>                  *
 * Searching criterion: maximize ------------- = ---- = -----------------                *
 *                                denominator     rr        <y1,y1>                      *
 *                                                                                       *
 * Notice: gain = gain_trans * gain_pit (computed in trans_enc() function)               *
 *                                                                                       *
 *---------------------------------------------------------------------------------------*/

void set_impulse(
    const float xn[],       /* i  : target signal                                   */
    const float h_orig[],   /* i  : impulse response of weighted synthesis filter   */
    float exc[],      /* o  : adaptive codebook excitation                    */
    float y1[],       /* o  : filtered adaptive codebook excitation           */
    short *imp_shape, /* o  : adaptive codebook index                         */
    short *imp_pos,   /* o  : position of the glotal impulse center index     */
    float *gain_trans /* o  : transition gain                                 */
)
{
    float rr[L_SUBFR];  /*  criterion: nominator coeficients   */
    float dd[L_SUBFR];  /*  criterion: denominator coeficients */
    float gh[L_SUBFR];  /*  convolution of 'g' and 'h' filters */
    float krit, krit_max;
    short i, j, m;
    short start1, start2, end1;


    krit_max = -1.0e+12f;

    /* main loop */
    /*  impulse  */
    for( m = 0; m < NUM_IMPULSE; m++ )
    {
        /* set searching ranges */
        if( *imp_pos<L_SUBFR-INPOL )
        {
            end1 = *imp_pos+INPOL;
        }
        else
        {
            end1 = L_SUBFR;
        }
        if( *imp_pos>INPOL )
        {
            start1 = *imp_pos-INPOL;
        }
        else
        {
            start1 = 0;
        }
        if( start1>L_IMPULSE2 )
        {
            start2 = start1;
        }
        else
        {
            start2 = L_IMPULSE2;
        }

        /*-----------------------------------------------------*
         *   nominator & DEnominator, gh=conv(g,h)
         *-----------------------------------------------------*/
        if( start1<L_IMPULSE2 )
        {
            rr[start1] = 0;
            dd[start1] = 0;
            convolve_tc(&glottal_cdbk[m*L_IMPULSE+L_IMPULSE2-start1], &h_orig[0], gh, (short) (L_IMPULSE-L_IMPULSE2+start1), L_SUBFR);

            /* nominator & Denominator row <0> */
            for( i=0; i < L_SUBFR; i++ )
            {
                rr[start1] += gh[i]*gh[i];
                dd[start1] += gh[i]*xn[i];
            }
            for( i=start1+1; i<L_IMPULSE2; i++ )
            {
                rr[i] = 0;
                dd[i] = 0;
                /* Denominator rows <1,L_IMPULSE2-1> */
                for( j=L_SUBFR-1; j > 0; j-- )
                {
                    gh[j]  = gh[j-1] + glottal_cdbk[m*L_IMPULSE+L_IMPULSE2-i]*h_orig[j];
                    rr[i] += gh[j]*gh[j];
                    dd[i] += gh[j]*xn[j];
                }
                gh[0]  = glottal_cdbk[m*L_IMPULSE+L_IMPULSE2-i]*h_orig[0];
                rr[i] += gh[0]*gh[0];
                dd[i] += gh[0]*xn[0];
                /* move rr and dd into rr[i] and dd[i] */
            }
            /* complete convolution(excitation,h_orig) */
            for( j=L_SUBFR-1; j > 0; j-- )
            {
                gh[j] = gh[j-1] + glottal_cdbk[m*L_IMPULSE]*h_orig[j];
            }
        }
        else
        {
            convolve_tc( &glottal_cdbk[m*L_IMPULSE], &h_orig[0], gh, L_IMPULSE, L_SUBFR );
        }
        if( end1>=start2 )
        {
            /* Denominator row <L_SUBFR-1> */
            rr[L_SUBFR-1] = 0;
            for (j=0; j <= L_IMPULSE2; j++)
            {
                rr[L_SUBFR-1] += gh[j]*gh[j];
            }
            /* move rr into rr[L_SUBFFR-1 */
            /* Denominator rows <L_IMPULSE2,L_SUBFR-2> */
            for( i=L_SUBFR-2; i >= start2; i-- )
            {
                rr[i] = rr[i+1] + gh[L_SUBFR+L_IMPULSE2-1-i]*gh[L_SUBFR+L_IMPULSE2-1-i];
            }

            /* nominator rows <L_IMPULSE2,L_SUBFR-1> */
            correlate_tc( xn, &dd[L_IMPULSE2], gh, (short) (start2-L_IMPULSE2), L_SUBFR, (short)(end1-L_IMPULSE2) );
        }

        /*------------------------------------------------------*
         *    maxim. criterion
         *------------------------------------------------------*/
        for( i=start1; i < end1; i++ )
        {
            krit = (float)(dd[i]*dd[i])/rr[i];
            if( krit > krit_max )
            {
                krit_max = krit;
                *imp_pos = i;
                *imp_shape = m;
            }
        }
    }

    /*--------------------------------------------------------*
     *    Build the excitation using found codeword
     *--------------------------------------------------------*/

    set_f(exc, 0, L_SUBFR);
    set_f(y1, 0, L_SUBFR);
    for( i=(*imp_pos-L_IMPULSE2); i<=(*imp_pos+L_IMPULSE2); i++ )
    {
        if( (i >= 0) && (i < L_SUBFR) )
        {
            exc[i] = glottal_cdbk[(*imp_shape)*L_IMPULSE+i-(*imp_pos)+L_IMPULSE2];
        }
    }

    /*------------------------------------------------------*
     *    Form filtered excitation, find gain_trans
     *------------------------------------------------------*/

    convolve_tc2( exc, h_orig, y1, *imp_pos );

    /* Find the ACELP correlations and the pitch gain (for current subframe) */
    *gain_trans = dotp( xn, y1, L_SUBFR )/(dotp( y1, y1, L_SUBFR ) + 0.01f);

    return;

}

/*-------------------------------------------------------------------*
 * convolve_tc:
 *
 *   convolution for different vectors' lengths
 *-------------------------------------------------------------------*/
static void convolve_tc(
    const float g[],    /* i  : input vector                              */
    const float h[],    /* i  : impulse response (or second input vector) */
    float y[],    /* o  : output vetor (result of convolution)      */
    const short L_1,    /* i  : vector h size                             */
    const short L_2     /* i  : vector g size                             */
)
{
    float temp;
    short   i, n;


    for( n = 0; n < L_2; n++ )
    {
        temp = g[0] * h[n];
        for( i = 1; i < ((n<L_1) ? (n+1) : L_1) ; i++ )
        {
            temp += g[i] * h[n-i];
        }
        y[n] = temp;
    }
    return;
}

/*-------------------------------------------------------------------*
 * convolve_tc2:
 *
 *   convolution for one vector with only L_IMPULSE nonzero coefficients
 *-------------------------------------------------------------------*/
static void convolve_tc2(
    const float g[],        /* i  : input vector                              */
    const float h[],        /* i  : impulse response (or second input vector) */
    float y[],        /* o  : output vetor (result of convolution)      */
    const short pos_max     /* o  : artificial impulse position               */
)

{
    float temp;
    short i, n;
    short i_start, i_end;
    short i_end2;


    i_start = pos_max-L_IMPULSE2;

    if( i_start<0 )
    {
        i_start = 0;
    }

    i_end = pos_max+L_IMPULSE;
    if (i_end > L_SUBFR)
    {
        i_end = L_SUBFR;
    }
    for( n = i_start; n < L_SUBFR; n++ )
    {
        temp = g[0] * h[n];
        i_end2 = ((n<=i_end) ? (n+1) : i_end);
        for( i = 1; i < i_end2; i++ )
        {
            temp += g[i] * h[n-i];
        }
        y[n] = temp;
    }
    return;
}

/*-------------------------------------------------------------------*
 * correlate_tc:
 *
 *   correlation for different vectors' lengths
 *-------------------------------------------------------------------*/
static void correlate_tc(
    const float *x,     /* i: target signal                                   */
    float       *y,     /* o: correlation between x[] and h[]                 */
    const float *h,     /* i: impulse response (of weighted synthesis filter) */
    const short start,  /* i: index of iterest                                */
    const short L_1,    /* i: vector size                                     */
    const short L_2     /* i: index of interest                               */
)
{
    short i, j;
    float s;


    for( i = start; i < L_2; i++ )
    {
        s = 0.0f;
        for( j = i; j < L_1; j++ )
        {
            s += x[j]*h[j-i];
        }
        y[i] = s;
    }

    return;
}
