/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "prot.h"
#include "rom_enc.h"




/*-------------------------------------------------------------------*
 * fft16()
 *
 *
 *-------------------------------------------------------------------*/

static void fft16(float *r_samp, float *i_samp)
{
    int i,j,N,Nv2,nm1,k;
    float tmpr[16],tmpi[16];
    float r1,s1,r2,s2;


    for (i = 0; i < 16; i++)
    {
        tmpr[i] = r_samp[i]*M_inr[i] - i_samp[i]*M_ini[i];
        tmpi[i] = r_samp[i]*M_ini[i] + i_samp[i]*M_inr[i];
    }

    for (i = 0; i < 8; i++)
    {
        s1 = tmpr[i] - tmpr[8+i];
        r1 = tmpi[i] - tmpi[8+i];

        tmpr[i] = tmpr[i] + tmpr[8+i];
        tmpi[i] = tmpi[i] + tmpi[8+i];

        tmpr[i+8] = s1*M_r[i] - r1*M_i[i];
        tmpi[i+8] = s1*M_i[i] + r1*M_r[i];
    }

    for (i = 0; i < 4; i++)
    {
        s1 = tmpr[i] - tmpr[4+i];
        r1 = tmpi[i] - tmpi[4+i];

        tmpr[i] = tmpr[i] + tmpr[4+i];
        tmpi[i] = tmpi[i] + tmpi[4+i];

        tmpr[i+4] = s1*M_r[2*i] - r1*M_i[2*i];
        tmpi[i+4] = s1*M_i[2*i] + r1*M_r[2*i];

    }
    for (i = 0; i < 4; i++)
    {
        s1 = tmpr[i+8] - tmpr[12+i];
        r1 = tmpi[i+8] - tmpi[12+i];

        tmpr[i+8] = tmpr[i+8] + tmpr[12+i];
        tmpi[i+8] = tmpi[i+8] + tmpi[12+i];

        tmpr[i+12] = s1*M_r[2*i] - r1*M_i[2*i];
        tmpi[i+12] = s1*M_i[2*i] + r1*M_r[2*i];

    }


    for (i = 0; i< 16; i=i+4)
    {
        s1 = tmpr[i] - tmpr[2+i];
        r1 = tmpi[i] - tmpi[2+i];
        s2 = tmpr[i+1] - tmpr[3+i];
        r2 = tmpi[i+1] - tmpi[3+i];

        tmpr[i] = tmpr[i] + tmpr[2+i];
        tmpi[i] = tmpi[i] + tmpi[2+i];
        tmpr[i+1] = tmpr[i+1] + tmpr[3+i];
        tmpi[i+1] = tmpi[i+1] + tmpi[3+i];

        tmpr[i+2] = s1*M_r[0] - r1*M_i[0];
        tmpi[i+2] = s1*M_i[0] + r1*M_r[0];
        tmpr[i+3] = s2*M_r[4] - r2*M_i[4];
        tmpi[i+3] = s2*M_i[4] + r2*M_r[4];
    }

    for (i = 0; i < 16; i= i+2)
    {
        s1 = tmpr[i] - tmpr[1+i];
        r1 = tmpi[i] - tmpi[1+i];

        tmpr[i] = tmpr[i] + tmpr[1+i];
        tmpi[i] = tmpi[i] + tmpi[1+i];

        tmpr[i+1] = s1*M_r[0] - r1*M_i[0];
        tmpi[i+1] = s1*M_i[0] + r1*M_r[0];
    }


    N = 16;
    Nv2=N>>1;
    nm1=N-1;
    j=0;
    for(i=0; i<nm1; i++)
    {
        if(i<j)
        {
            r_samp[j] = tmpr[i]*M_Wr[j] - tmpi[i]*M_Wi[j];
            i_samp[j] = tmpr[i]*M_Wi[j] + tmpi[i]*M_Wr[j];
            r_samp[i] = tmpr[j]*M_Wr[i] - tmpi[j]*M_Wi[i];
            i_samp[i] = tmpr[j]*M_Wi[i] + tmpi[j]*M_Wr[i];
        }
        else if (i==j)
        {
            r_samp[i] = tmpr[i]*M_Wr[i] - tmpi[i]*M_Wi[i];
            i_samp[i] = tmpr[i]*M_Wi[i] + tmpi[i]*M_Wr[i];
        }

        k=Nv2;
        while(k<=j)
        {
            j-=k;
            k>>=1;
        }
        j+=k;
    }

    r_samp[15] = tmpr[15]*M_Wr[15] - tmpi[15]*M_Wi[15];
    i_samp[15] = tmpr[15]*M_Wi[15] + tmpi[15]*M_Wr[15];

    return;
}

/*-------------------------------------------------------------------*
 * subband_FFT()
 *
 *
 *-------------------------------------------------------------------*/

void subband_FFT(
    float Sr[16][60],            /*(i) real part of the cldfb*/
    float Si[16][60],            /*(i) imag part of the cldfb*/
    float *spec_amp              /*(o) spectral amplitude*/
)
{
    int i,j;
    float tmpr[16],tmpi[16];
    float ptmp[16],tmp1;
    for(i=0; i<10; i++)
    {
        for(j=0; j<16; j++)
        {
            tmpr[j] = Sr[j][i];
            tmpi[j] = Si[j][i];
        }

        fft16(tmpr,tmpi);

        for(j=0; j<16; j++)
        {
            ptmp[j] = tmpr[j]*tmpr[j] + tmpi[j]*tmpi[j];
        }
        if(i%2==0)
        {
            for(j=0; j<8; j++)
            {
                tmp1 = ptmp[j] + ptmp[15-j] ;
                spec_amp[i*8+j] = (float)sqrt(tmp1);
            }
        }
        else
        {
            for(j=0; j<8; j++)
            {
                tmp1 = ptmp[j] + ptmp[15-j] ;
                spec_amp[i*8+7-j] =(float)sqrt(tmp1);
            }
        }
    }

    return;
}




