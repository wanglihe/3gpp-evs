/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/


#include <math.h>
#include "options.h"
#include "prot.h"

/*-------------------------------------------------------------------
 * ham_cos_window()
 *
 *
 *-------------------------------------------------------------------*/

void ham_cos_window(
    float *fh,
    int n1,
    int n2
)
{
    float cc, cte;
    int i;
    cte = PI2/(float)(2*n1 - 1);
    cc = 0.0f;
    for (i = 0; i < n1; i++)
    {
        *fh++ = 0.54f - 0.46f * (float)cos(cc);
        cc += cte;
    }
    cte = PI2/(float)(4*n2 - 1);
    cc = 0.0f;
    for (i = 0; i < n2; i++)
    {
        *fh++ = (float)cos(cc);
        cc += cte;
    }
    return;
}
