/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "prot.h"
#include "rom_com.h"


/*-----------------------------------------------------------------*
 * Local constants
 *-----------------------------------------------------------------*/

#define N_MAX_FFT     1024
#define INV_SQR2  0.70710676908493f

/*---------------------------------------------------------------------*
 * ifft_rel()
 *
 * Calculate the inverse FFT of a real signal
 *
 * Based on the FORTRAN code from the article "Real-valued Fast Fourier Transform Algorithms"
 * by Sorensen, ... in IEEE Trans. on ASSP, Vol. ASSP-35, No. June 6th 1987.
 *
 * Input: the io[] signal containing the spectrum in the following order :
 *
 * Re[0], Re[1], ..  Re[n/2], Im[n/2-1], .. Im[1]
 *---------------------------------------------------------------------*/

void ifft_rel(
    float io[],  /* i/o: input/output vector   */
    const short n,     /* i  : vector length         */
    const short m      /* i  : log2 of vector length */
)
{
    short i, j, k;
    short step;
    short n2, n4, n8, i0;
    short is, id;
    float *x,*xi0, *xi1, *xi2, *xi3, *xi4, *xup1, *xdn6, *xup3, *xdn8;
    float xt;
    float r1;
    float t1, t2, t3, t4, t5;
    float cc1, cc3, ss1, ss3;
    const float *s, *s3, *c, *c3;
    const short *idx;
    float temp[512];
    float n_inv;

    n_inv = 1.0f/n;

    /*-----------------------------------------------------------------*
     * IFFT
     *-----------------------------------------------------------------*/

    x = &io[-1];
    n2 = 2*n;
    for (k=1; k<m; k++)
    {
        is = 0;
        id = n2;
        n2 = n2 >> 1;
        n4 = n2 >> 2;
        n8 = n4 >> 1;
        while (is < n-1)
        {
            xi1 = x + is + 1;
            xi2 = xi1 + n4;
            xi3 = xi2 + n4;
            xi4 = xi3 + n4;

            for (i=is; i<n; i+= id)
            {
                t1 = *xi1 - *xi3;
                *xi1 += *xi3;
                *xi2 = 2.0f**xi2;
                *xi3 = t1- 2.0f**xi4;
                *xi4 = t1 + 2.0f**xi4;
                if (n4 != 1)
                {
                    t1 = (*(xi2+n8) - *(xi1+n8))*INV_SQR2;
                    t2 = (*(xi4+n8) + *(xi3+n8))*INV_SQR2;

                    *(xi1+n8) += *(xi2+n8);
                    *(xi2+n8) = *(xi4+n8) - *(xi3+n8);
                    *(xi3+n8) = (float)(2.0f * (-t2-t1));
                    *(xi4+n8) = (float)(2.0f * (-t2+t1));
                }
                xi1 += id;
                xi2 += id;
                xi3 += id;
                xi4 += id;
            }
            is = 2*id - n2;
            id  = 4*id;
        }
        step  = N_MAX_FFT/n2;

        s = sincos_t_ext + step;
        c = s + N_MAX_FFT/4;
        s3 = sincos_t_ext + 3*step;
        c3 = s3 + N_MAX_FFT/4;
        for (j=2; j<=n8; j++)
        {
            cc1 = *c ;
            ss1 = *s;
            cc3 = *c3;
            ss3 = *s3;

            is  = 0;
            id  = 2 * n2;

            c += step;
            s += step;

            c3 += 3*step;
            s3 += 3*step;
            while (is < n-1)
            {
                xup1 = x + j + is;
                xup3 = xup1 + 2*n4;
                xdn6 = xup3 - 2*j +2;
                xdn8 = xdn6 + 2*n4;

                for (i=is; i<n; i+=id)
                {
                    t1 = *xup1 - *xdn6;
                    *xup1 = *xup1 + *xdn6;
                    xup1 += n4;
                    xdn6 -= n4;

                    t2 = *xdn6 - *xup1;
                    *xdn6 = *xup1 + *xdn6;

                    xdn6 += n4;
                    t3 = *xdn8 + *xup3;
                    *xdn6  = *xdn8 - *xup3;

                    xup3 += n4;
                    xdn8 -= n4;

                    t4 =  *xup3 + *xdn8;
                    *xup1= *xup3 - *xdn8;

                    t5 = t1 - t4;
                    t1 = t1 + t4;
                    t4 = t2 - t3;
                    t2 = t2 + t3;
                    *xup3 = t1*cc3 - t2*ss3;
                    xup3 -= n4;
                    *xup3 = t5*cc1 + t4*ss1;
                    *xdn8 = -t4*cc1 + t5*ss1;

                    xdn8 += n4;
                    *xdn8 = t2*cc3 + t1*ss3;

                    xup1 -= n4;
                    xup1 += id;
                    xup3 += id;
                    xdn6 += id;
                    xdn8 += id;
                }
                is = 2*id - n2;
                id = 4*id;
            }
        }
    }

    /*-----------------------------------------------------------------*
     * Length two butterflies
     *-----------------------------------------------------------------*/

    is = 1;
    id = 4;
    while (is < n)
    {
        xi0 = x + is ;
        xi1 = xi0 + 1;

        for (i0=is; i0<=n; i0+=id)
        {
            r1 = *xi0;
            *xi0= r1 + *xi1;
            *xi1 = r1 - *xi1;
            xi0 += id;
            xi1 += id;
        }
        is = 2*id - 1;
        id = 4*id;
    }

    /*-----------------------------------------------------------------*
     * Digit reverse counter
     *-----------------------------------------------------------------*/

    idx = fft256_read_indexes;
    xi0 = temp-1;
    if (n == 128)
    {
        for (i=0; i<n; i++)
        {
            j = *idx++;
            temp[i] = x[1+(j>>1)];
        }
    }
    else if (n == 256)
    {
        for (i=0; i<n; i++)
        {
            j = *idx++;
            temp[i] = x[1+j];
        }
    }
    else if (n == 512)
    {
        for (i=0; i<256; i++)
        {
            j = *idx++;
            temp[i] = x[1+2*j];
            temp[i+256] = x[2+2*j];
        }
    }
    else
    {
        xi0 = x;
        j = 1;
        for (i=1; i<n; i++)
        {
            if (i < j)
            {
                xt = x[j];
                x[j] = x[i];
                x[i] = xt;
            }
            k = n >> 1;
            while (k < j)
            {
                j = j - k;
                k = k >> 1;
            }
            j = j + k;
        }
    }

    /*-----------------------------------------------------------------*
     * Normalization
     *-----------------------------------------------------------------*/

    for (i=1; i<=n; i++)
    {
        x[i] = xi0[i] * n_inv;
    }

    return;
}
