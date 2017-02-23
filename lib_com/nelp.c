/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include "math.h"
#include "rom_com.h"
#include "prot.h"

/*-------------------------------------------------------------------*
 * dequantize_uvg()
 *
 * Dequantize unvoiced gains
 *--------------------------------------------------------------------*/

float dequantize_uvg(
    int   iG1,                   /* i: gain 1 index      */
    int   *iG2,                  /* i: gain 2 index      */
    float *G,                    /* o: quantized gain    */
    short bwidth
)
{
    short i, k;
    const float (*UVG1CB)[2] = NULL;
    const float (*UVG2CB1)[5] = NULL;
    const float (*UVG2CB2)[5] = NULL;

    if ( bwidth == NB )
    {
        UVG1CB = UVG1CB_NB;
        UVG2CB1 = UVG2CB1_NB;
        UVG2CB2 = UVG2CB2_NB;
    }

    else if ( bwidth == WB || bwidth == SWB )
    {
        UVG1CB = UVG1CB_WB;
        UVG2CB1 = UVG2CB1_WB;
        UVG2CB2 = UVG2CB2_WB;
    }

    for( i=0; i<2; i++)
    {
        for (k=0; k<5; k++)
        {
            if(i==0)
            {
                G[i*5+k] = (float) pow(10.0, UVG1CB[iG1][i]) * UVG2CB1[iG2[i]][k];
            }
            else if (i==1)
            {
                G[i*5+k] = (float) pow(10.0, UVG1CB[iG1][i]) * UVG2CB2[iG2[i]][k];
            }
        }
    }
    return(0.0);
}

/*-------------------------------------------------------------------*
 * generate_nelp_excitation()
 *
 * Generate excitation for NELP coding.
 *--------------------------------------------------------------------*/

void generate_nelp_excitation(
    short *seed,            /* i/o: random number seed    */
    float *Gains,           /* i  : excitation gains      */
    float *output,          /* o  : excitation output     */
    float gain_fac          /* i  : gain factor           */
)
{
    short i,len,j;
    float tmp[31], tmp1[31], tmpf;
    short k1,k2, I[31], tmpi;

    for (i=0; i<10; i++)
    {
        if (i==9)
        {
            len=31;
        }
        else
        {
            len=25;
        }

        for (j=0; j<len; j++)
        {
            tmp[j]=((*seed)=521*(*seed)+259)/32768.0f;
            tmp1[j]=ABSVAL(tmp[j]);
            I[j]=j;
        }
        for (k1=0; k1<len-1; k1++)
        {
            for (k2=k1+1; k2<len; k2++)
            {
                if (tmp1[k2]>tmp1[k1])
                {
                    tmpi=I[k2];
                    tmpf=tmp1[k2];
                    tmp1[k2]=tmp1[k1];
                    I[k2]=I[k1];
                    tmp1[k1]=tmpf;
                    I[k1]=tmpi;
                }
            }
        }

        /*using a factor of 1.37 to compensate for the ~ 2.5 ( or 2.73) dB diff between this scheme and EVS-UV */
        for (j=0; j<(short) rint_new(len/4.0f); j++)
        {
            output[i*25+I[j]]=(float) (Gains[i]*sqrt(3.0f)*tmp[I[j]]*gain_fac);
        }
        for (; j<len; j++)
        {
            output[i*25+I[j]]=0;
        }
    }

}
