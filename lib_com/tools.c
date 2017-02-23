/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>
#include <ctype.h>
#include <math.h>
#include <assert.h>
#include "options.h"
#include "prot.h"


/*------------------------------------------------------------------*
 * own_random()
 *
 * Random generator
 *------------------------------------------------------------------*/

short own_random(       /* o  : output random value */
    short *seed          /* i/o: random seed         */
)
{
    *seed = (short) (*seed * 31821L + 13849L);

    return(*seed);
}

/*---------------------------------------------------------------------
 * sign()
 *
 *---------------------------------------------------------------------*/

float sign(             /* o  : sign of x (+1/-1) */
    const float x       /* i  : input value of x  */
)
{
    if (x < 0.0f)
    {
        return -1.0f;
    }
    else
    {
        return 1.0f;
    }
}

/*---------------------------------------------------------------------
 * log2_f()
 *
 *---------------------------------------------------------------------*/

float log2_f(            /* o  : logarithm2 of x */
    const float x        /* i  : input value of x   */
)
{
    return (float)(log(x)/log(2.0f));
}

short norm_ul( unsigned int UL_var1)
{
    short  var_out;

    if (UL_var1 == 0)
    {
        var_out = 0;
    }
    else
    {
        for (var_out = 0; UL_var1 < (unsigned int) 0x80000000U; var_out++)
        {
            UL_var1 <<= 1;
        }
    }
    BASOP_CHECK();
    return (var_out);
}





/*---------------------------------------------------------------------
 * sum_s()
 * sum_f()
 *
 *---------------------------------------------------------------------*/

short sum_s(           /* o  : sum of all vector elements            */
    const short *vec,  /* i  : input vector                          */
    const short lvec   /* i  : length of input vector                */
)
{
    short i;
    short tmp;

    tmp = 0;
    for( i=0; i<lvec; i++ )
    {
        tmp += vec[i];
    }

    return tmp;
}


float sum_f(           /* o  : sum of all vector elements            */
    const float *vec,  /* i  : input vector                          */
    const short lvec   /* i  : length of input vector                */
)
{
    short i;
    float tmp;

    tmp = 0.0f;
    for( i=0; i<lvec; i++ )
    {
        tmp += vec[i];
    }

    return tmp;
}

/*----------------------------------------------------------------------
 * sum2_f()
 *
 *---------------------------------------------------------------------*/

float sum2_f(          /* o  : sum of all squared vector elements    */
    const float *vec,  /* i  : input vector                          */
    const short lvec   /* i  : length of input vector                */
)
{
    short i;
    float tmp;

    tmp = 0.0f;
    for( i=0; i<lvec; i++ )
    {
        tmp += vec[i] * vec[i];
    }

    return tmp;
}

/*-------------------------------------------------------------------*
 * set_c()
 * set_s()
 * set_f()
 * set_i()
 *
 * Set the vector elements to a value
 *-------------------------------------------------------------------*/

void set_c(
    char y[],  /* i/o: Vector to set                       */
    const char a,    /* i  : Value to set the vector to          */
    const short N     /* i  : Length of the vector                */
)
{
    short i;

    for (i=0 ; i<N ; i++)
    {
        y[i] = a;
    }

    return;
}


void set_s(
    short y[],  /* i/o: Vector to set                       */
    const short a,    /* i  : Value to set the vector to          */
    const short N     /* i  : Length of the vector                */
)
{
    short i;

    for (i=0 ; i<N ; i++)
    {
        y[i] = a;
    }

    return;
}

void set_i(
    int y[],  /* i/o: Vector to set                       */
    const int a,    /* i  : Value to set the vector to          */
    const short N     /* i  : Length of the vector                */
)
{
    short i;

    for (i=0 ; i<N ; i++)
    {
        y[i] = a;
    }

    return;
}

void set_f(
    float y[],  /* i/o: Vector to set                       */
    const float a,    /* i  : Value to set the vector to          */
    const short N     /* i  : Lenght of the vector                */
)
{
    short i ;

    for (i=0 ; i<N ; i++)
    {
        y[i] = a;
    }

    return;
}

/*---------------------------------------------------------------------*
 * set_zero()
 *
 * Set a vector vec[] of dimension lvec to zero
 *---------------------------------------------------------------------*/

void set_zero(
    float *vec, /* o  : input vector         */
    int lvec  /* i  : length of the vector */
)
{
    int i;

    for (i = 0; i < lvec; i++)
    {
        *vec++ = 0.0f;
    }

    return;
}


/*---------------------------------------------------------------------*
 * mvr2r()
 * mvs2s()
 * mvr2s()
 * mvs2r()
 * mvi2i()
 *
 * Transfer the contents of vector x[] to vector y[]
 *---------------------------------------------------------------------*/

void mvr2r(
    const float x[],  /* i  : input vector  */
    float y[],  /* o  : output vector */
    const short n     /* i  : vector size   */
)
{
    short i;

    if ( n <= 0 )
    {
        /* cannot transfer vectors with size 0 */
        return;
    }

    if (y < x)
    {
        for (i = 0; i < n; i++)
        {
            y[i] = x[i];
        }
    }
    else
    {
        for (i = n-1; i >= 0; i--)
        {
            y[i] = x[i];
        }
    }

    return;
}

void mvs2s(
    const short x[],  /* i  : input vector  */
    short y[],  /* o  : output vector */
    const short n     /* i  : vector size   */
)
{
    short i;

    if ( n <= 0 )
    {
        /* cannot transfer vectors with size 0 */
        return;
    }

    if (y < x)
    {
        for (i = 0; i < n; i++)
        {
            y[i] = x[i];
        }
    }
    else
    {
        for (i = n-1; i >= 0; i--)
        {
            y[i] = x[i];
        }
    }

    return;
}

unsigned int mvr2s(
    const float x[],  /* i  : input vector  */
    short y[],  /* o  : output vector */
    const short n     /* i  : vector size   */
)
{
    short i;
    float temp;
    unsigned int noClipping = 0;

    if ( n <= 0 )
    {
        /* cannot transfer vectors with size 0 */
        return 0;
    }

    if ((void*)y < (const void*)x)
    {
        for (i = 0; i < n; i++)
        {
            temp = x[i];
            temp = (float)floor(temp + 0.5f);

            if (temp >  32767.0f )
            {
                temp =  32767.0f;
                noClipping++;
            }
            else if (temp < -32768.0f )
            {
                temp = -32768.0f;
                noClipping++;
            }

            y[i] = (short)temp;
        }
    }
    else
    {
        for (i = n-1; i >= 0; i--)
        {
            temp = x[i];
            temp = (float)floor(temp + 0.5f);

            if (temp >  32767.0f )
            {
                temp =  32767.0f;
                noClipping++;
            }
            else if (temp < -32768.0f )
            {
                temp = -32768.0f;
                noClipping++;
            }

            y[i] = (short)temp;
        }
    }

    return noClipping;
}

void mvs2r(
    const short x[],  /* i  : input vector  */
    float y[],  /* o  : output vector */
    const short n     /* i  : vector size   */
)
{
    short i;

    if ( n <= 0 )
    {
        /* cannot transfer vectors with size 0 */
        return;
    }

    if ((void*)y < (const void*)x)
    {
        for (i = 0; i < n; i++)
        {
            y[i] = (float)x[i];
        }
    }
    else
    {
        for (i = n-1; i >= 0; i--)
        {
            y[i] = (float)x[i];
        }
    }

    return;
}


void mvi2i(
    const int x[],  /* i  : input vector  */
    int y[],  /* o  : output vector */
    const int n     /* i  : vector size   */
)
{
    int i;

    if ( n <= 0 )
    {
        /* no need to transfer vectors with size 0 */
        return;
    }

    if (y < x)
    {
        for (i = 0; i < n; i++)
        {
            y[i] = x[i];
        }
    }
    else
    {
        for (i = n-1; i >= 0; i--)
        {
            y[i] = x[i];
        }
    }

    return;
}

/*---------------------------------------------------------------------*
 * maximum()
 *
 * Find index and value of the maximum in a vector
 *---------------------------------------------------------------------*/

short maximum(         /* o  : index of the maximum value in the input vector */
    const float *vec,  /* i  : input vector                                   */
    const short lvec,  /* i  : length of input vector                         */
    float *max   /* o  : maximum value in the input vector              */
)
{
    short j, ind;
    float tmp;

    ind = 0;
    tmp = vec[0];

    for ( j=1; j<lvec; j++ )
    {
        if( vec[j] > tmp )
        {
            ind = j;
            tmp = vec[j];
        }
    }

    if ( max != 0 )
    {
        *max = tmp;
    }

    return ind;
}

/*---------------------------------------------------------------------*
 * minimum()
 *
 * Find index of a minimum in a vector
 *---------------------------------------------------------------------*/

short minimum(         /* o  : index of the minimum value in the input vector */
    const float *vec,  /* i  : input vector                                   */
    const short lvec,  /* i  : length of input vector                         */
    float *min   /* o  : minimum value in the input vector              */
)
{
    short j, ind;
    float tmp;

    ind = 0;
    tmp = vec[0];

    for ( j=1 ; j<lvec ; j++ )
    {
        if( vec[j] < tmp )
        {
            ind = j;
            tmp = vec[j];
        }
    }

    if ( min != 0 )
    {
        *min = tmp;
    }

    return ind;
}

/*---------------------------------------------------------------------*
 * emaximum()
 *
 * Find index of a maximum energy in a vector
 *---------------------------------------------------------------------*/

short emaximum(            /* o  : return index with max energy value in vector   */
    const float *vec,      /* i  : input vector                                   */
    const short lvec,      /* i  : length of input vector                         */
    float *ener_max  /* o  : maximum energy value                           */
)
{
    short j, ind;
    float temp;

    *ener_max = 0.0f;
    ind = 0;

    for( j=0; j<lvec; j++ )
    {
        temp = vec[j] * vec[j];

        if( temp > *ener_max )
        {
            ind = j;
            *ener_max = temp;
        }
    }

    return ind;
}


/*---------------------------------------------------------------------*
 * mean()
 *
 * Find the mean of the vector
 *---------------------------------------------------------------------*/

float mean(            /* o  : mean of vector                         */
    const float *vec,  /* i  : input vector                           */
    const short lvec   /* i  : length of input vector                 */
)
{
    float tmp;

    tmp = sum_f( vec, lvec ) / (float)lvec;

    return tmp;
}

/*---------------------------------------------------------------------*
 * dotp()
 *
 * Dot product of vector x[] and vector y[]
 *---------------------------------------------------------------------*/

float dotp(            /* o  : dot product of x[] and y[]    */
    const float  x[],  /* i  : vector x[]                    */
    const float  y[],  /* i  : vector y[]                    */
    const short  n     /* i  : vector length                 */
)
{
    short i;
    float suma;

    suma = x[0] * y[0];

    for(i=1; i<n; i++)
    {
        suma += x[i] * y[i];
    }

    return suma;
}



/*---------------------------------------------------------------------*
 * inv_sqrt()
 *
 * Find the inverse square root of the input value
 *---------------------------------------------------------------------*/

float inv_sqrt(    /* o  : inverse square root of input value */
    const float x  /* i  : input value                        */
)
{
    return (float)(1.0 / sqrt(x));
}

/*-------------------------------------------------------------------*
 * conv()
 *
 * Convolution between vectors x[] and h[] written to y[]
 * All vectors are of length L. Only the first L samples of the
 * convolution are considered.
 *-------------------------------------------------------------------*/

void conv(
    const float x[],   /* i  : input vector                              */
    const float h[],   /* i  : impulse response (or second input vector) */
    float y[],   /* o  : output vetor (result of convolution)      */
    const short L      /* i  : vector size                               */
)
{
    float temp;
    short i, n;
    for (n = 0; n < L; n++)
    {
        temp = x[0] * h[n];
        for (i = 1; i <= n; i++)
        {
            temp += x[i] * h[n-i];
        }
        y[n] = temp;
    }

    return;
}

/*-------------------------------------------------------------------*
 * fir()
 *
 * FIR filtering of vector x[] with filter having impulse response h[]
 * written to y[]
 * The input vector has length L and the FIR filter has an order of K, i.e.
 * K+1 coefficients. The memory of the input signal is provided in the vector mem[]
 * which has K values
 * The maximum length of the input signal is L_FRAME32k and the maximum order
 * of the FIR filter is L_FILT_MAX
 *-------------------------------------------------------------------*/

void fir(
    const float x[],   /* i  : input vector                              */
    const float h[],   /* i  : impulse response of the FIR filter        */
    float y[],   /* o  : output vector (result of filtering)       */
    float mem[], /* i/o: memory of the input signal (L samples)    */
    const short L,     /* i  : input vector size                         */
    const short K,     /* i  : order of the FIR filter (K+1 coefs.)      */
    const short upd    /* i  : 1 = update the memory, 0 = not            */
)
{
    float buf_in[L_FRAME48k+60], buf_out[L_FRAME48k], s;
    short i, j;

    /* prepare the input buffer (copy and update memory) */
    mvr2r( mem, buf_in, K );
    mvr2r( x, buf_in + K, L );

    if ( upd )
    {
        mvr2r( buf_in + L, mem, K );
    }

    /* do the filtering */
    for ( i = 0; i < L; i++ )
    {
        s = buf_in[K+i] * h[0];

        for (j = 1; j <= K; j++)
        {
            s += h[j] * buf_in[K+i-j];
        }

        buf_out[i] = s;
    }

    /* copy to the output buffer */
    mvr2r( buf_out, y, L );

    return;
}

/*-------------------------------------------------------------------*
 * v_add()
 *
 * Addition of two vectors sample by sample
 *-------------------------------------------------------------------*/

void v_add(
    const float x1[],   /* i  : Input vector 1                       */
    const float x2[],   /* i  : Input vector 2                       */
    float y[],    /* o  : Output vector that contains vector 1 + vector 2  */
    const short N       /* i  : Vector lenght                                    */
)
{
    short i ;

    for (i=0 ; i<N ; i++)
    {
        y[i] = x1[i] + x2[i] ;
    }

    return;
}

/*-------------------------------------------------------------------*
 * v_sub()
 *
 * Subtraction of two vectors sample by sample
 *-------------------------------------------------------------------*/

void v_sub(
    const float x1[],   /* i  : Input vector 1                                   */
    const float x2[],   /* i  : Input vector 2                                   */
    float y[],    /* o  : Output vector that contains vector 1 - vector 2  */
    const short N       /* i  : Vector lenght                                    */
)
{
    short i ;

    for (i=0 ; i<N ; i++)
    {
        y[i] = x1[i] - x2[i] ;
    }

    return;
}

/*-------------------------------------------------------------------*
 * v_mult()
 *
 * Multiplication of two vectors
 *-------------------------------------------------------------------*/

void v_mult(
    const float x1[],   /* i  : Input vector 1                                   */
    const float x2[],   /* i  : Input vector 2                                   */
    float y[],    /* o  : Output vector that contains vector 1 .* vector 2 */
    const short N       /* i  : Vector lenght                                    */
)
{
    short i ;

    for (i=0 ; i<N ; i++)
    {
        y[i] = x1[i] * x2[i] ;
    }

    return;
}

/*-------------------------------------------------------------------*
 * v_multc()
 *
 * Multiplication of vector by constant
 *-------------------------------------------------------------------*/

void v_multc(
    const float x[],    /* i  : Input vector                                     */
    const float c,      /* i  : Constant                                         */
    float y[],    /* o  : Output vector that contains c*x                  */
    const short N       /* i  : Vector length                                    */
)
{
    short i ;

    for (i=0 ; i<N ; i++)
    {
        y[i] = c * x[i];
    }

    return;
}


/*-------------------------------------------------------------------*
 * squant()
 *
 * Scalar quantizer according to MMSE criterion (nearest neighbour in Euclidean space)
 *
 * Searches a given codebook to find the nearest neighbour in Euclidean space.
 * Index of the winning codeword and the winning codeword itself are returned.
 *-------------------------------------------------------------------*/

int squant(             /* o: index of the winning codeword   */
    const float x,      /* i: scalar value to quantize        */
    float *xq,    /* o: quantized value                 */
    const float cb[],   /* i: codebook                        */
    const int   cbsize  /* i: codebook size                   */
)
{
    float dist, mindist, tmp;
    int   c, idx;
    idx = 0;
    mindist = 1e16f;

    for ( c = 0; c < cbsize; c++)
    {
        dist = 0.0f;
        tmp = x - cb[c];
        dist += tmp*tmp;
        if (dist < mindist)
        {
            mindist = dist;
            idx = c;
        }
    }

    *xq = cb[idx];

    return idx;
}


/*-------------------------------------------------------------------*
 * usquant()
 *
 * Uniform scalar quantizer according to MMSE criterion
 * (nearest neighbour in Euclidean space)
 *
 * Applies quantization based on scale and round operations.
 * Index of the winning codeword and the winning codeword itself are returned.
 *-------------------------------------------------------------------*/

short usquant(          /* o: index of the winning codeword   */
    const float x,      /* i: scalar value to quantize        */
    float *xq,    /* o: quantized value                 */
    const float qlow,   /* i: lowest codebook entry (index 0) */
    const float delta,  /* i: quantization step               */
    const short cbsize  /* i: codebook size                   */
)
{
    short idx;

    idx = (short) max(0.f, min(cbsize - 1, ( (x - qlow)/delta + 0.5f)));
    *xq = idx*delta + qlow;

    return idx;
}

/*-------------------------------------------------------------------*
 * usdequant()
 *
 * Uniform scalar de-quantizer routine
 *
 * Applies de-quantization based on scale and round operations.
 *-------------------------------------------------------------------*/

float usdequant(
    const int idx,      /* i: quantizer index                 */
    const float qlow,   /* i: lowest codebook entry (index 0) */
    const float delta   /* i: quantization step               */
)
{
    float g;

    g = idx * delta + qlow;

    return( g );
}


/*-------------------------------------------------------------------*
 * vquant()
 *
 * Vector quantizer according to MMSE criterion (nearest neighbour in Euclidean space)
 *
 * Searches a given codebook to find the nearest neighbour in Euclidean space.
 * Index of the winning codevector and the winning vector itself are returned.
 *-------------------------------------------------------------------*/

int vquant(                   /* o: index of the winning codevector */
    float x[],          /* i: vector to quantize              */
    const float x_mean[],     /* i: vector mean to subtract (0 if none) */
    float xq[],         /* o: quantized vector                */
    const float cb[],         /* i: codebook                        */
    const int   dim,          /* i: dimension of codebook vectors   */
    const int   cbsize        /* i: codebook size                   */
)
{
    float dist, mindist, tmp;
    int   c, d, idx, j;

    idx = 0;
    mindist = 1e16f;

    if (x_mean != 0)
    {
        for( d = 0; d < dim; d++ )
        {
            x[d] -= x_mean[d];
        }
    }

    j = 0;
    for ( c = 0; c < cbsize; c++ )
    {
        dist = 0.0f;
        for( d = 0; d < dim; d++ )
        {
            tmp = x[d] - cb[j++];
            dist += tmp*tmp;
        }

        if ( dist < mindist )
        {
            mindist = dist;
            idx = c;
        }
    }

    if ( xq == 0 )
    {
        return idx;
    }

    j =  idx*dim;
    for ( d = 0; d < dim; d++)
    {
        xq[d] = cb[j++];
    }

    if (x_mean != 0)
    {
        for( d = 0; d < dim; d++)
        {
            xq[d] += x_mean[d];
        }
    }

    return idx;
}

/*-------------------------------------------------------------------*
 * w_vquant()
 *
 * Vector quantizer according to MMSE criterion (nearest neighbour in Euclidean space)
 *
 * Searches a given codebook to find the nearest neighbour in Euclidean space.
 * Weights are put on the error for each vector element.
 * Index of the winning codevector and the winning vector itself are returned.
 *-------------------------------------------------------------------*/

int w_vquant(                 /* o: index of the winning codevector */
    float x[],          /* i: vector to quantize              */
    const float x_mean[],     /* i: vector mean to subtract (0 if none) */
    const short weights[],    /* i: error weights                   */
    float xq[],         /* o: quantized vector                */
    const float cb[],         /* i: codebook                        */
    const int   dim,          /* i: dimension of codebook vectors   */
    const int   cbsize,       /* i: codebook size                   */
    const short rev_vect      /* i: reverse codebook vectors        */
)
{
    float dist, mindist, tmp;
    int   c, d, idx, j, k;

    idx = 0;
    mindist = 1e16f;

    if (x_mean != 0)
    {
        for( d = 0; d < dim; d++)
        {
            x[d] -= x_mean[d];
        }
    }

    j = 0;
    if (rev_vect)
    {
        k = dim-1;
        for ( c = 0; c < cbsize; c++)
        {
            dist = 0.0f;

            for( d = k; d >= 0; d--)
            {
                tmp = x[d] - cb[j++];
                dist += weights[d]*(tmp*tmp);
            }

            if (dist < mindist)
            {
                mindist = dist;
                idx = c;
            }
        }

        if (xq == 0)
        {
            return idx;
        }

        j =  idx*dim;
        for ( d = k; d >= 0; d--)
        {
            xq[d] = cb[j++];
        }
    }
    else
    {
        for ( c = 0; c < cbsize; c++)
        {
            dist = 0.0f;

            for( d = 0; d < dim; d++)
            {
                tmp = x[d] - cb[j++];
                dist += weights[d]*(tmp*tmp);
            }

            if (dist < mindist)
            {
                mindist = dist;
                idx = c;
            }
        }

        if (xq == 0)
        {
            return idx;
        }

        j =  idx*dim;
        for ( d = 0; d < dim; d++)
        {
            xq[d] = cb[j++];
        }
    }

    if (x_mean != 0)
    {
        for( d = 0; d < dim; d++)
        {
            xq[d] += x_mean[d];
        }
    }

    return idx;
}


/*----------------------------------------------------------------------------------*
 * v_sort()
 *
 * Sorting of vectors. This is very fast with almost ordered vectors.
 *----------------------------------------------------------------------------------*/

void v_sort(
    float *r,    /* i/o: Vector to be sorted in place */
    const short lo,    /* i  : Low limit of sorting range   */
    const short up     /* I  : High limit of sorting range  */
)
{
    short i, j;
    float tempr;

    for ( i=up-1; i>=lo; i-- )
    {
        tempr = r[i];
        for ( j=i+1; j<=up && (tempr>r[j]); j++ )
        {
            r[j-1] = r[j];
        }

        r[j-1] = tempr;
    }

    return;
}

/*---------------------------------------------------------------------*
 * var()
 *
 * Calculate the variance of a vector
 *---------------------------------------------------------------------*/

float var(                /* o: variance of vector                    */
    const float *x,        /* i: input vector                          */
    const int len          /* i: length of inputvector                 */
)
{
    float m;
    float v;
    short i;

    m = mean(x, (const short)len);

    v = 0.0f;
    for (i = 0; i < len; i++)
    {
        v += (x[i] - m)*(x[i] - m);
    }
    v /= len;

    return v;
}


/*---------------------------------------------------------------------*
 * std_dev()
 *
 * Calculate the standard deviation of a vector
 *---------------------------------------------------------------------*/

float std_dev(             /* o: standard deviation                    */
    const float *x,        /* i: input vector                          */
    const int len          /* i: length of the input vector            */
)
{
    short i;
    float std;

    std = 1e-16f;
    for( i = 0; i < len; i++)
    {
        std += x[i] * x[i];
    }

    std = (float)sqrt( std / len );

    return std;
}



/*---------------------------------------------------------------------*
 * dot_product_mat()
 *
 * Calculates dot product of type x'*A*x, where x is column vector of size m,
 * and A is square matrix of size m*m
 *---------------------------------------------------------------------*/

float dot_product_mat(    /* o  : the dot product x'*A*x        */
    const float  *x,      /* i  : vector x                      */
    const float  *A,      /* i  : matrix A                      */
    const short  m        /* i  : vector & matrix size          */
)
{
    short i,j;
    float suma, tmp_sum;
    const float *pt_x, *pt_A;

    pt_A = A;
    suma = 0;

    for(i=0; i<m; i++)
    {
        tmp_sum = 0;
        pt_x = x;
        for (j=0; j<m; j++)
        {
            tmp_sum += *pt_x++ **pt_A++;
        }

        suma += x[i] * tmp_sum;
    }

    return suma;
}

/*--------------------------------------------------------------------------------*
 * polezero_filter()
 *
 * Y(Z)=X(Z)(b[0]+b[1]z^(-1)+..+b[L]z^(-L))/(a[0]+a[1]z^(-1)+..+a[M]z^(-M))
 *          mem[n]=x[n]+cp[0]mem[n-1]+..+cp[M-1]mem[n-M], where cp[i]=-a[i+1]/a[0]
 *          y[n]=cz[0]mem[n]+cz[1]mem[n-1]+..+cz[L]mem[n-L], where cz[i]=b[i]/a[0]
 *          mem={mem[n-K] mem[n-K+1] . . . . mem[n-2] mem[n-1]}, where K=max(L,M)
 *
 * a[0] must be equal to 1.0f!
 *---------------------------------------------------------------------------------*/

void polezero_filter (
    const float *in,            /* i  : input vector                              */
    float *out,           /* o  : output vector                             */
    const short N,              /* i  : input vector size                         */
    const float *b,             /* i  : numerator coefficients                    */
    const float *a,             /* i  : denominator coefficients                  */
    const short order,          /* i  : filter order                              */
    float *mem            /* i/o: filter memory                             */
)
{
    short i, j, k;


    for (i=0; i<order; i++)
    {
        out[i] = in[i] * b[0];
        for (j = 0; j < i; j++)
        {
            out[i] += in[i-1-j] * b[j+1] - out[i-1-j] * a[j+1];
        }

        for (k = order-1; j < order; j++, k--)
        {
            out[i] += mem[k] * b[j+1] - mem[k+order] * a[j+1];
        }
    }

    for (; i < N; i++)
    {
        out[i] = in[i] * b[0];
        for (j = 0; j < order; j++)
        {
            out[i] += in[i-1-j] * b[j+1] - out[i-1-j] * a[j+1];
        }
    }

    for (i=0; i<order; i++)
    {
        mem[i] = in[N - order + i];
        mem[i+order] = out[N - order + i];
    }

    return;
}


static
float fleft_shift(float input, int shift)
{
    return (input * (float) pow(2.0, (double) shift));
}

static
float fright_shift(float input, int shift)
{
    return (input * (float) pow(0.5, (double) shift));
}



/*--------------------------------------------------------------------------------*
 * root_a()
 *
 * Implements a quadratic approximation to sqrt(a)
 * Firstly, a is normalized to lie between 0.25 & 1.0
 * by shifting the input left or right by an even number of
 * shifts. Even shifts represent powers of 4 which, after
 * the sqrt, can easily be converted to powers of 2 and are
 * easily dealt with.
 * At the heart of the algorithm is a quadratic
 * approximation of the curve sqrt(a) for 0.25 <= a <= 1.0.
 * Sqrt(a) approx = 0.27 + 1.0127 * a - 0.2864 * a^2
 *
 *---------------------------------------------------------------------------------*/

float root_a(float a)
{
    short int shift_a;
    float mod_a;
    float approx;

    if (a <= 0.0f)
    {
        return 0.0;
    }


    /* This next piece of code implements a "norm" function */
    /* and returns the shift needed to scale "a" to have a  */
    /* 1 in the (MSB-1) position. This is equivalent to     */
    /* giving a value between 0.5 & 1.0.                    */
    mod_a = a;

    shift_a = 0;
    while (mod_a > 1.0)
    {
        mod_a /= 2.0;
        shift_a--;
    }

    while (mod_a < 0.5)
    {
        mod_a *= 2.0;
        shift_a++;
    }


    shift_a &= 0xfffe;
    mod_a = fleft_shift(a, (int) shift_a);

    approx = 0.27f + 1.0127f * mod_a - 0.2864f * mod_a * mod_a;

    approx = fright_shift(approx, (int) (shift_a >> 1));

    return (approx);
}

/*--------------------------------------------------------------------------------*
 * root_a_over_b()
 *
 * Implements an approximation to sqrt(a/b)
 * Firstly a & b are normalized to lie between 0.25 & 1.0
 * by shifting the inputs left or right by an even number
 * of shifts.
 * Even shifts represent powers of 4 which, after the sqrt,
 * become powers of 2 and are easily dealt with.
 * At the heart of the algorithm is an approximation of the
 * curve sqrt(a/b) for 0.25 <= a <= 1.0 & 0.25 <= b <= 1.0.
 * Given the value of b, the 2nd order coefficients p0, p1
 * & p2 can be determined so that...
 * Sqrt(a/b) approx = p0 + p1 * a + p2 * a^2
 * where p0 approx =  0.7176 - 0.8815 * b + 0.4429 * b^2
 *       p1 approx =  2.6908 - 3.3056 * b + 1.6608 * b^2
 *       p2 approx = -0.7609 + 0.9346 * b - 0.4695 * b^2
 *
 *---------------------------------------------------------------------------------*/

float root_a_over_b(float a, float b)
{
    short int shift_a, shift_b, shift;
    float mod_a, mod_b;
    float p2 = -0.7609f;
    float p1 = 2.6908f;
    float p0 = 0.7176f;
    float b_sqr;
    float approx;

    if ((a <= 0.0f) || (b <= 0.0f))
    {
        return 0.0;
    }

    a += 0x00000001;
    b += 0x00000001;

    /* This next piece of code implements a "norm" function */
    /* and returns the shift needed to scale "a" to have a  */
    /* 1 in the (MSB-1) position. This is equivalent to     */
    /* giving a value between 0.5 & 1.0.                    */
    mod_a = a;

    shift_a = 0;
    while (mod_a > 1.0)
    {
        mod_a /= 2.0;
        shift_a--;
    }

    while (mod_a < 0.5)
    {
        mod_a *= 2.0;
        shift_a++;
    }

    shift_a &= 0xfffe;
    mod_a = fleft_shift(a, (int) shift_a);

    /* This next piece of code implements a "norm" function */
    /* and returns the shift needed to scale "b" to have a  */
    /* 1 in the (MSB-1) position. This is equivalent to     */
    /* giving a value between 0.5 & 1.0.                    */
    mod_b = b;

    shift_b = 0;
    while (mod_b > 1.0)
    {
        mod_b /= 2.0;
        shift_b--;
    }

    while (mod_b < 0.5)
    {
        mod_b *= 2.0;
        shift_b++;
    }

    shift_b &= 0xfffe;
    mod_b = fleft_shift(b, (int) shift_b);

    shift = (shift_b - shift_a) >> 1;

    b_sqr = mod_b * mod_b;

    p2 += 0.9346f * mod_b + -0.4695f * b_sqr;
    p1 += -3.3056f * mod_b + 1.6608f * b_sqr;
    p0 += -0.8815f * mod_b + 0.4429f * b_sqr;

    approx = p0 + p1 * mod_a + p2 * mod_a * mod_a;

    approx = fleft_shift(approx, (int) shift);

    return (approx);
}

/*--------------------------------------------------------------------------------*
 * rint_new()
 *
 * Round to the nearest integer with mid-point exception
 *---------------------------------------------------------------------------------*/

double rint_new(
    double x
)
{
    int a;

    /* middle value point test */
    if (ceil (x + 0.5) == floor (x + 0.5))
    {
        a = (int) ceil (x);

        if (a % 2 == 0)
        {
            return ceil (x);
        }
        else
        {
            return floor (x);
        }
    }
    else
    {
        return floor (x + 0.5);
    }
}


/*-------------------------------------------------------------------*
 * anint()
 *
 * Round to the nearest integer.
 *-------------------------------------------------------------------*/

double anint(
    double x
)
{
    return (x)>=0?(int)((x)+0.5):(int)((x)-0.5);
}

/*-------------------------------------------------------------------*
 * is_numeric_float()
 *
 * Returns 0 for all NaN and Inf values defined according to IEEE 754
 * floating point number's definition. Returns 1 for numeric values.
 *-------------------------------------------------------------------*/
short is_numeric_float(
    float x
)
{
    union float_int
    {
        float float_val;
        int int_val;
    } float_int;
    float_int.float_val=x;
    return ((float_int.int_val & 0x7f800000) != 0x7f800000);
}


