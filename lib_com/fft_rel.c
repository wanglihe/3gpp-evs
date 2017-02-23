/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "prot.h"
#include "rom_com.h"

/*---------------------------------------------------------------------*
 * Local constants
 *---------------------------------------------------------------------*/

#define N_MAX_FFT        1024
#define N_MAX_DIV2  (N_MAX_FFT>>1)
#define N_MAX_DIV4  (N_MAX_DIV2>>1)

/*---------------------------------------------------------------------*
 *  fft_rel()
 *
 *  Computes the split-radix FFT in place for the real-valued
 *  signal x of length n.  The algorithm has been ported from
 *  the Fortran code of [1].
 *
 *  The function  needs sine and cosine tables t_sin and t_cos,
 *  and the constant N_MAX_FFT.  The table  entries  are defined as
 *  sin(2*pi*i) and cos(2*pi*i) for i = 0, 1, ..., N_MAX_FFT-1. The
 *  implementation  assumes  that any entry  will not be needed
 *  outside the tables. Therefore, N_MAX_FFT and n must be properly
 *  set.  The function has been tested  with the values n = 16,
 *  32, 64, 128, 256, and N_MAX_FFT = 1280.
 *
 *  References
 *  [1] H.V. Sorensen,  D.L. Jones, M.T. Heideman, C.S. Burrus,
 *      "Real-valued fast  Fourier transform  algorithm,"  IEEE
 *      Trans. on Signal Processing,  Vol.35, No.6, pp 849-863,
 *      1987.
 *
 *  OUTPUT
 *      x[0:n-1]  Transform coeffients in the order re[0], re[1],
 *                ..., re[n/2], im[n/2-1], ..., im[1].
 *---------------------------------------------------------------------*/

void fft_rel(
    float x[],  /* i/o: input/output vector    */
    const short n,    /* i  : vector length          */
    const short m     /* i  : log2 of vector length  */
)
{
    short i, j, k, n1, n2, n4;
    short step;
    float xt, t1, t2;
    float *x0, *x1, *x2;
    float *xi2, *xi3, *xi4, *xi1;
    const float *s, *c;
    const short *idx;

    /* !!!! NMAX = 256 is hardcoded here  (similar optimizations should be done for NMAX > 256) !!! */

    float *x2even, *x2odd;
    float temp[512];

    if ( n == 128 || n == 256 || n == 512 )
    {
        idx = fft256_read_indexes;

        /* Combined Digit reverse counter & Length two butterflies */
        if (n == 128)
        {
            x2 = temp;
            for (i = 0; i < 64; i++)
            {
                j = *idx++;
                k = *idx++;

                *x2++ = x[j>>1] + x[k>>1];
                *x2++ = x[j>>1] - x[k>>1];
            }
        }
        else if (n == 256)
        {
            x2 = temp;
            for (i = 0; i < 128; i++)
            {
                j = *idx++;
                k = *idx++;

                *x2++ = x[j] + x[k];
                *x2++ = x[j] - x[k];
            }
        }
        else if (n == 512)
        {
            x2even = temp;
            x2odd = temp + 256;

            for (i = 0; i < 128; i++)
            {
                j = 2 * *idx++;
                k = 2 * *idx++;

                *x2even++ = x[j] + x[k];
                *x2even++ = x[j] - x[k];
                *x2odd++ = x[++j] + x[++k];
                *x2odd++ = x[j] - x[k];
            }
        }

        /*-----------------------------------------------------------------*
         * 1st Stage Loop has been Unrolled because n4 is '1' and that
         * allows the elimination of the 'for_ (j = 1; j < n4; j++)' loop
         * and the associated pointers initialization.
         * Also, it allows to Put the Data from 'temp' back into 'x' due
         * to the previous Combined Digit Reverse and Length two butterflies
         *-----------------------------------------------------------------*/

        /*for_ (k = 2; k < 3; k++)*/
        {
            x0 = temp;
            x1 = x0 + 2;
            x2 = x;

            for (i = 0; i < n; i += 4)
            {
                *x2++ = *x0++ + *x1;    /* x[i] = xt + x[i+n2];    */
                *x2++ = *x0;
                *x2++ = *--x0 - *x1++;  /* x[i+n2] = xt - x[i+n2];      */
                *x2++ = -*x1;         /* x[i+n2+n4] = -x[i+n2+n4];     */

                x0 += 4;
                x1 += 3; /* x1 has already advanced */
            }
        }
    }
    else
    {
        /*-----------------------------------------------------------------*
         * Digit reverse counter
         *-----------------------------------------------------------------*/

        j = 0;
        x0 = &x[0];
        for (i = 0; i < n-1; i++)
        {
            if (i < j)
            {
                xt   = x[j];
                x[j] = *x0;
                *x0  = xt;
            }
            x0++;
            k = n/2;
            while (k <= j)
            {
                j -= k;
                k  = k>>1;
            }
            j += k;
        }

        /*-----------------------------------------------------------------*
         * Length two butterflies
         *-----------------------------------------------------------------*/

        x0 = &x[0];
        x1 = &x[1];
        for (i = 0; i < n/2; i++)
        {
            *x1 = *x0   - *x1;
            *x0 = *x0*2 - *x1;

            x0++;
            x0++;
            x1++;
            x1++;
        }

        /*-----------------------------------------------------------------*
         * 1st Stage Loop has been Unrolled because n4 is '1' and that
         * allows the elimination of the 'for_ (j = 1; j < n4; j++)' loop
         * and the associated pointers initialization.
         *-----------------------------------------------------------------*/

        /* for_ (k = 2; k < 3; k++) */
        {
            x0 = x;
            x1 = x0 + 2;

            for (i = 0; i < n; i += 4)
            {
                *x1 = *x0   - *x1;      /* x[i+n2] = xt - x[i+n2];      */
                *x0 = *x0*2 - *x1++;    /* x[i] = xt + x[i+n2];    */
                *x1 = -*x1;             /* x[i+n2+n4] = -x[i+n2+n4];     */

                x0 += 4;
                x1 += 3; /* x1 has already advanced */
            }
        }
    }

    /*-----------------------------------------------------------------*
     * Other butterflies
     *
     * The implementation described in [1] has been changed by using
     * table lookup for evaluating sine and cosine functions.  The
     * variable ind and its increment step are needed to access table
     * entries.  Note that this implementation assumes n4 to be so
     * small that ind will never exceed the table.  Thus the input
     * argument n and the constant N_MAX_FFT must be set properly.
     *-----------------------------------------------------------------*/

    n4 = 1;
    n2 = 2;
    n1 = 4;

    step = N_MAX_DIV4;

    for (k = 3; k <= m; k++)
    {
        step >>= 1;
        n4 <<= 1;
        n2 <<= 1;
        n1 <<= 1;

        x0 = x;
        x1 = x0 + n2;
        x2 = x1 + n4;

        for (i = 0; i < n; i += n1)
        {
            *x1 = *x0   - *x1;      /* x[i+n2] = xt - x[i+n2];      */
            *x0 = *x0*2 - *x1;      /* x[i] = xt + x[i+n2];    */
            *x2 = -*x2;             /* x[i+n2+n4] = -x[i+n2+n4];     */

            s = sincos_t_ext;
            c = s + N_MAX_FFT/4;          /* 1024/4 = 256, 256/4=64 */
            xi1 = x0;
            xi3 = xi1 + n2;
            xi2 = xi3;
            x0 += n1;
            xi4 = x0;

            for (j = 1; j < n4; j++)
            {
                xi3++;
                xi1++;
                xi4--;
                xi2--;
                c += step;
                s+= step;          /* autoincrement by ar0 */

                t1 = *xi3**c + *xi4**s;    /* t1 = *xi3**(pt_c+ind) + *xi4**(pt_s+ind);   */
                t2 = *xi3**s - *xi4**c;    /* t2 = *xi3**(pt_s+ind) - *xi4**(pt_c+ind);     */

                *xi4 =  *xi2 - t2;
                *xi2 =  *xi1 - t1;
                *xi1 =  *xi1*2 - *xi2;
                *xi3 = -2*t2 - *xi4;
            }

            x1 += n1;
            x2 += n1;
        }
    }

    return;
}
