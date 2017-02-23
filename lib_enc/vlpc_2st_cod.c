/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "cnst.h"
#include "prot.h"

/*------------------------------------------------------------------*
* vlpc_2st_cod()
*
*
*------------------------------------------------------------------*/

int vlpc_2st_cod(   /* output: number of allocated bits        */
    const float *lsf, /* input:  normalized vector to quantize   */
    float *lsfq,      /* i/o:    i:1st stage   o:1st+2nd stage   */
    int *indx,        /* output: index[] (4 bits per words)      */
    int mode,         /* input:  0=abs, >0=rel                   */
    float sr_core
)
{
    int    i, nbits;
    float  w[M], x[M], tmp;
    int nq, xq[M];
    float scale = sr_core/INT_FS_12k8;
    float scaleinv = 1.f / scale;

    /* 0 bit with true weighting: save 0.5 bit */

    lsf_weight_2st( lsf, w, 1, sr_core );

    for (i=0; i<M; i++)
    {
        x[i] = (lsf[i] - lsfq[i])*scaleinv;
    }

    for (i=0; i<M; i++)
    {
        x[i] /= w[i];
    }

    tmp = 0.0f;

    for (i=0; i<M; i++)
    {
        tmp += x[i]*x[i];
    }

    if (tmp < 8.0f)
    {
        indx[0] = 0;
        indx[1] = 0;

        return(10);       /* 2*(2+3) */
    }

    /* weighting from the 1st stage */
    lsf_weight_2st( lsfq, w, mode, sr_core );

    /* find lsf residual */

    for (i=0; i<M; i++)
    {
        x[i] = (lsf[i] - lsfq[i])*scaleinv;
    }

    /* scale the residual */

    for (i=0; i<M; i++)
    {
        x[i] /= w[i];
    }

    /* quantize */
    AVQ_cod_lpc( x, xq, indx, 2 );

    /* quantized lsf */

    for (i=0; i<M; i++)
    {
        lsfq[i] += scale*(w[i]*(float)xq[i]);
    }

    /* total number of bits using entropic code to index the quantizer number */
    nbits = 0;

    for (i=0; i<2; i++)
    {
        nq = indx[i];
        nbits += (2+(nq*4));             /* 2 bits to specify Q2,Q3,Q4,ext */

        if (nq > 6) nbits += nq-3;       /* unary code (Q7=1110, ...) */
        else if (nq > 4) nbits += nq-4;  /* Q5=0, Q6=10 */
        else if (nq == 0) nbits += 3;    /* Q0=110 */
    }

    /* reorder */
    v_sort( lsfq, 0, M-1 );
    reorder_lsf( lsfq, LSF_GAP, M, sr_core );

    return( nbits );
}
