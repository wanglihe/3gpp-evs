/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <assert.h>
#include <math.h>
#include "cnst.h"
#include "prot.h"


/*------------------------------------------------------------------*
* lsf_weight_2st()
*
*
*------------------------------------------------------------------*/

void lsf_weight_2st(
    const float *lsfq,
    float *w,
    int mode,
    float sr_core
)
{
    int i;
    float d[M+1];
    float freq_max = sr_core / 2.f;
    float freq_div = freq_max / (float)M;

    /* compute lsf distance */
    d[0] = lsfq[0];
    for (i=1; i<M; i++)
    {
        d[i] = lsfq[i] - lsfq[i-1];
    }
    d[M] = freq_max - lsfq[M-1];

    /* weighting function */
    for( i=0; i<M; i++ )
    {
        assert(d[i]>0);

        if (mode == 0)
        {
            w[i] = (float) (60.0f / (freq_div/sqrt(d[i]*d[i+1])));    /* abs */
        }
        else if (mode == 1)
        {
            w[i] = (float) (65.0f / (freq_div/sqrt(d[i]*d[i+1])));    /* mid */
        }
        else
        {
            w[i] = (float) (63.0f / (freq_div/sqrt(d[i]*d[i+1])));    /* rel2 */
        }
    }

    return;
}
