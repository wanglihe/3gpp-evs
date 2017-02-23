/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"

/*--------------------------------------------------------------------*
 * residu()
 *
 * Compute the LP residual by filtering the input speech through A(z)
 *--------------------------------------------------------------------*/

void residu(
    const float *a,  /* i  : LP filter coefficients           */
    const short m,   /* i  : order of LP filter               */
    const float *x,  /* i  : input signal (usually speech)    */
    float *y,  /* o  : output signal (usually residual) */
    const short l    /* i  : size of filtering                */
)
{
    float s;
    short i, j;

    for (i = 0; i < l; i++)
    {
        s = x[i];
        for (j = 1; j <= m; j++)
        {
            s += a[j]*x[i-j];
        }
        y[i] = s;
    }

    return;
}

/*--------------------------------------------------------------------*
 * calc_residu()
 *
 * Compute the LP residual by filtering the input through A(z) in all subframes
 *--------------------------------------------------------------------*/

void calc_residu(
    const float *speech,       /* i  : weighted speech signal                    */
    float *res,          /* o  : residual signal                           */
    const float *p_Aq,         /* i  : quantized LP filter coefficients          */
    const short L_frame        /* i  : size of frame                             */
)
{
    short i_subfr;

    for( i_subfr = 0; i_subfr < L_frame; i_subfr += L_SUBFR )
    {
        /* calculate the residual signal */
        residu( p_Aq, M, &speech[i_subfr], &res[i_subfr], L_SUBFR );

        /* next subframe */
        p_Aq += (M+1);
    }
    return;
}
