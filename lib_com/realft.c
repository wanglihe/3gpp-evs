/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "cnst.h"
#include "options.h"
#include "prot.h"

/*-------------------------------------------------------------------*
 * four1()
 *
 *  From "numerical recipes in C".
 *  Replace data by its DFT, if isign is input as 1; or replace data
 *  by nn times its inverse-DFT, if isign is input as -1.
 *  data is a complex array of length nn, input as a real
 *  array data[1...2nn]. nn must be an integer power of 2
 *-------------------------------------------------------------------*/

static void four1(
    float *data,     /* i/o: data array   .......... */
    short nn,        /* i  : length of data array    */
    short isign      /* i  : sign +1 or -1           */
)
{
    short n,mmax,m,j,istep,i;
    float wtemp,wr,wpr,wpi,wi,theta;
    float tempr,tempi;

    n=nn << 1;
    j=1;

    /* this is the bit-reversal section of the routine */
    for (i=1; i<n; i+=2)
    {
        if (j > i)
        {
            /* exchange the two complex numbers */
            SWAP(data[j],data[i]);
            SWAP(data[j+1],data[i+1]);
        }
        m=n >> 1;
        while (m >= 2 && j > m)
        {
            j -= m;
            m >>= 1;
        }
        j += m;
    }
    mmax=2;
    /* here begins the Danielson-Lanczos section of the routine */
    /* Outer loop executed log2(nn) times */
    while (n > mmax)
    {
        istep=2*mmax;
        /* initialization for the trigonometric recurrence */
        theta=(float) (6.28318530717959/(isign*mmax));
        wtemp=(float) (sin(0.5f*theta));
        wpr = -2.0f*wtemp*wtemp;
        wpi=(float) sin(theta);
        wr=1.0f;
        wi=0.0f;
        /* here are the two nested loops */
        for (m=1; m<mmax; m+=2)
        {
            for (i=m; i<=n; i+=istep)
            {
                /* this is Danielson-Lanczos formula */
                j=i+mmax;
                tempr=wr*data[j]-wi*data[j+1];
                tempi=wr*data[j+1]+wi*data[j];
                data[j]=data[i]-tempr;
                data[j+1]=data[i+1]-tempi;
                data[i] += tempr;
                data[i+1] += tempi;
            }
            /* trigonometric recurrence */
            wr=(wtemp=wr)*wpr-wi*wpi+wr;
            wi=wi*wpr+wtemp*wpi+wi;
        }
        mmax=istep;
    }

    return;
}

/*-------------------------------------------------------------------------*
 * realft()
 *
 * from "numerical recipes in C".
 * Calculates the Fourier Transform of a set of 2*n real-valued data points.
 * Replaces this data (which is stored in the array data[1..2n]) by the
 * positive frequancy half of its complex Fourier Transform. The real-valued
 * first and last components of the complex transform are returned as elements
 * data[1] and data[2] respectively. n must be a power of 2. This routine
 * also calculates the inverse transform of a complex data array if it is the
 * tranform of real data. (Results in this case must be multiplied by 1/n.)
 *--------------------------------------------------------------------------*/

void realft(
    float *data,               /* i/o: data array   .......... */
    short n,                   /* i  : length of data array    */
    short isign                /* i  : sign +1 or -1           */
)
{
    short i,i1,i2,i3,i4,n2p3;
    float c1=0.5,c2,h1r,h1i,h2r,h2i;
    float wr,wi,wpr,wpi,wtemp,theta;

    theta=(float) (EVS_PI/(float) n);
    if (isign == 1)
    {
        /* the forward transorm here */
        c2 = -0.5;
        four1(data,n,1);
    }
    else
    {
        /* otherwise set up for the inverse transform */
        c2=0.5;
        theta = -theta;
    }
    wtemp=(float) sin(0.5f*theta);
    wpr = -2.0f*wtemp*wtemp;
    wpi=(float) sin(theta);
    wr=1.0f+wpr;
    wi=wpi;
    n2p3=2*n+3;
    /* case i=1 done separately below */
    for (i=2; i<=n/2; i++)
    {
        i4=1+(i3=n2p3-(i2=1+(i1=i+i-1)));
        /* the two separate transforms are separated out of data */
        h1r=c1*(data[i1]+data[i3]);
        h1i=c1*(data[i2]-data[i4]);
        h2r = -c2*(data[i2]+data[i4]);
        h2i=c2*(data[i1]-data[i3]);
        /* here they are recombined to form the true transform
        of the original real data */
        data[i1]=h1r+wr*h2r-wi*h2i;
        data[i2]=h1i+wr*h2i+wi*h2r;
        data[i3]=h1r-wr*h2r+wi*h2i;
        data[i4] = -h1i+wr*h2i+wi*h2r;
        /* the recurrence */
        wr=(wtemp=wr)*wpr-wi*wpi+wr;
        wi=wi*wpr+wtemp*wpi+wi;
    }
    if (isign == 1)
    {
        /* squeeze the first and the last data together to get them
        all within the original data */
        data[1] = (h1r=data[1])+data[2];
        data[2] = h1r-data[2];
    }
    else
    {
        /* this is the inverse transform for the case isign=-1 */
        data[1]=c1*((h1r=data[1])+data[2]);
        data[2]=c1*(h1r-data[2]);
        four1(data,n,-1);
    }

    return;
}
