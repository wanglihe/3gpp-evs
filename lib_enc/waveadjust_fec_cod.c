/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "prot.h"


/*-------------------------------------------------------------------*
 * SFM_Cal()
 *
 *
 *--------------------------------------------------------------------*/

float SFM_Cal(
    float const fcoef[],
    int n
)
{
    int i, k;
    double geoMean = 0, arithMean = 0, SFM;
    double tmp[4] = {1, 1, 1, 1};

    for (k = 0; k < 4; k++)
    {

        for (i = k*(n>>2); i < (k+1)*(n>>2); i++)
        {
            tmp[k] *= (fabs(fcoef[i]) + EPSILON);
            arithMean += fabs(fcoef[i]);
        }
    }
    geoMean = log(tmp[0])+log(tmp[1])+log(tmp[2])+log(tmp[3]);
    geoMean = geoMean/n;
    geoMean = exp(geoMean);
    arithMean /= n;
    SFM = (geoMean+EPSILON)/(arithMean+EPSILON);

    return (float)SFM;
}
