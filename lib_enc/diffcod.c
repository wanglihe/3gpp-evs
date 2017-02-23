/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "prot.h"
#include "rom_com.h"


/*--------------------------------------------------------------------------*/
/*  Function  diffcod                                                       */
/*  ~~~~~~~~~~~~~~~~~                                                       */
/*                                                                          */
/*  Differential coding for indices of quantized norms                      */
/*--------------------------------------------------------------------------*/
/*  short    *normqlg2  (o)   quantized norms in log2                       */
/*  short    N          (i)   number of sub-vectors                         */
/*  short    *y         (i/o) indices of quantized norms                    */
/*  short    *difidx    (o)   differential code                             */
/*--------------------------------------------------------------------------*/

void diffcod(
    const short N,
    short *y,
    short *difidx
)
{
    short i, k, r;

    for (i=N-1; i>0; i--)
    {
        r = i - 1;
        k = y[i] - y[r];
        if (k<(-15))
        {
            y[r] = y[i] + 15;
        }
    }

    for (i=1; i<N; i++)
    {
        r = i - 1;
        k = y[i] - y[r];
        if (k>16)
        {
            k = 16;
            y[i] = y[r] + 16;
        }
        difidx[r] = k + 15;
    }

    return;
}




/*--------------------------------------------------------------------------
 * diffcod_lrmdct()
 *
 * Differential coding for indices of quantized norms
 *--------------------------------------------------------------------------*/

void diffcod_lrmdct(
    const short N,                  /* i  : number of sub-vectors       */
    const int be_ref,               /* o  : band energy reference */
    int *y,                   /* i/o: indices of quantized norms */
    int *difidx,              /* o  : differential code */
    const short is_transient        /* i  : transient flag  */
)
{
    short i, m, r;
    int k;
    int thr_l,thr_h;

    if(is_transient)
    {
        thr_l=-15;
        thr_h=16;
    }
    else
    {
        thr_l=-32;
        thr_h=31;
    }

    difidx[0] = y[0] - be_ref;
    if(difidx[0]>thr_h)
    {
        difidx[0] = thr_h;
        y[0] = be_ref + thr_h;
    }

    if(difidx[0]<thr_l)
    {
        difidx[0] = thr_l;
        y[0] = be_ref + thr_l;
    }

    m = N - 1;
    for (i=m; i>0; i--)
    {
        r = i - 1;
        k = y[i] - y[r];
        if (k<thr_l)
        {
            y[r] = y[i] - thr_l;
        }
    }

    for (i=1; i<N; i++)
    {
        r = i - 1;
        k = y[i] - y[r];
        if (k>thr_h)
        {
            k = thr_h;
            y[i] = y[r] + thr_h;
        }
        difidx[i] = k;
    }

    return;
}

