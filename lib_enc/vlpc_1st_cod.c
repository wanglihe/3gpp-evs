/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <assert.h>
#include "cnst.h"
#include "prot.h"


extern float const dico_lsf_abs_8b[];

/*------------------------------------------------------------------*
* lsf_weight()
*
* outputs only the weightings, doesn't do anything with the lsfq
*------------------------------------------------------------------*/

static void lsf_weight(
    const float *lsfq,
    float *w,
    float sr_core
)
{
    int i;
    float inv_di0, inv_di1;
    float scale = sr_core/INT_FS_12k8;
    float freq_max = sr_core/2.f;

    /* Verify, that M is pair, otherwise adapt exit of loop below */
    assert ((M & 1) == 0);

    /* weighting function */
    inv_di0 = scale / lsfq[0];
    for (i=1; i<(M-2); i+=2)
    {
        inv_di1 = scale / (lsfq[i] - lsfq[i-1]);
        w[i-1] = inv_di0 + inv_di1;

        inv_di0 = scale / (lsfq[i+1] - lsfq[i]);
        w[i] = inv_di1 + inv_di0;
    }
    inv_di1 = scale / (lsfq[i] - lsfq[i-1]);
    w[i-1] = inv_di0 + inv_di1;

    inv_di0 = scale / (freq_max - lsfq[i]);
    w[i] = inv_di1 + inv_di0;

    return;
}

/*------------------------------------------------------------------*
* vlpc_1st_cod()
*
*
*------------------------------------------------------------------*/

int vlpc_1st_cod(       /* output: codebook index                  */
    const float *lsf,     /* input:  vector to quantize              */
    float *lsfq,          /* i/o:    i:prediction   o:quantized lsf  */
    float sr_core
    ,float *wout           /* o: lsf weights */
)
{
    int    i, j, index;
    float  w[M], x[M];
    float dist_min, dist, temp;
    const float *p_dico;
    int hit=0;
    float scale = sr_core/INT_FS_12k8;
    float scaleinv = 1.f / scale;

    /* weighting */
    lsf_weight(lsf, w, sr_core );

    mvr2r( w, wout, M );

    /* remove lsf prediction/means */

    for( i=0; i<M; i++ )
    {
        x[i] = (lsf[i] - lsfq[i])*scaleinv;
    }

    dist_min = 1.0e30f;
    p_dico = dico_lsf_abs_8b;
    index = 0;

    for (i = 0; i < 256; i++)
    {
        dist = 0.0;

        for (j = 0; j < M; j++)
        {
            temp = x[j] - *p_dico++;
            dist += w[j] * temp * temp;
        }

        if (dist < dist_min)
        {
            dist_min = dist;
            index = i;
            hit++; /*just for testing*/
        }
    }

    /* quantized vector */
    p_dico = &dico_lsf_abs_8b[index * M];

    for (j = 0; j < M; j++)
    {
        lsfq[j] += scale **p_dico++; /* += cause it's differential */
    }

    assert(index < 256);
    assert(hit>0);

    return index;
}
