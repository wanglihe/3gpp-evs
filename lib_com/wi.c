/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <math.h>
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"


/*-------------------------------------------------------------------*
* Local constants
*
*
*-------------------------------------------------------------------*/

#define WARP_OS_RATE 8
#define LL 256
#define LL_OS (WARP_OS_RATE*LL)
#define OSLENGTH 12
#define CUTFREE_ABS_RANGE 6
#define CUTFREE_REL_RANGE 0.25
#define WI_THRESHLD 0.8
#define WI_SAMPLE_THLD 20
#define ERB_CBSIZE1 64
#define ERB_CBSIZE2 64
#define P_CBSIZE 64

#define _POLY1(x, c)    ((c)[0] * (x) + (c)[1])
#define _POLY2(x, c)    (_POLY1((x), (c)) * (x) + (c)[2])
#define _POLY3(x, c)    (_POLY2((x), (c)) * (x) + (c)[3])


/*-------------------------------------------------------------------*
* DTFS_new()
*
* DTFS structure initialization.
*-------------------------------------------------------------------*/

DTFS_STRUCTURE* DTFS_new(
    void
)
{
    short i ;
    DTFS_STRUCTURE* dtfs = NULL;

    dtfs = (DTFS_STRUCTURE *) calloc(1,sizeof(DTFS_STRUCTURE));
    dtfs->lag = 0 ;
    dtfs->nH=0;
    dtfs->nH_4kHz=0;
    dtfs->upper_cut_off_freq_of_interest=3300.0;
    dtfs->upper_cut_off_freq=4000.0;
    for(i=0; i<MAXLAG_WI; i++)
    {
        dtfs->a[i]=dtfs->b[i]=0.0 ;
    }
    return dtfs;     /* o: DTFS structure  */
}


/*-------------------------------------------------------------------*
* DTFS_copy()
*
* Copy from one DTFS STRUCTURE to another.
*-------------------------------------------------------------------*/

void DTFS_copy(
    DTFS_STRUCTURE *Xout,  /* o: DTFS structure  */
    DTFS_STRUCTURE Xinp    /* i: DTFS structure  */
)
{
    short k;
    for(k=0; k<MAXLAG_WI; k++)
    {
        Xout->a[k]=Xinp.a[k];
    }
    for(k=0; k<MAXLAG_WI; k++)
    {
        Xout->b[k]=Xinp.b[k];
    }
    Xout->lag=Xinp.lag;
    Xout->nH=Xinp.nH;
    Xout->nH_4kHz=Xinp.nH_4kHz;
    Xout->upper_cut_off_freq_of_interest=Xinp.upper_cut_off_freq_of_interest;
    Xout->upper_cut_off_freq=Xinp.upper_cut_off_freq;
}


/*-------------------------------------------------------------------*
* DTFS_sub()
*
* Difference of A and B coefficients in cartesian domain.
* Equivalent to time domain subtraction.
*-------------------------------------------------------------------*/

DTFS_STRUCTURE DTFS_sub(
    DTFS_STRUCTURE X1,  /* i: DTFS input 1 */
    DTFS_STRUCTURE X2   /* i: DTFS input 2 */
)
{
    DTFS_STRUCTURE tmp ;
    short i ;
    for(i=0; i<=X1.lag/2; i++)
    {
        tmp.a[i] = X1.a[i] ;
        tmp.b[i] = X1.b[i] ;
    }
    for(i=0; i<=X2.lag/2; i++)
    {
        tmp.a[i] -= X2.a[i] ;
        tmp.b[i] -= X2.b[i] ;
    }
    tmp.lag = max(X1.lag, X2.lag) ;
    tmp.nH=max(X1.nH,X2.nH);
    tmp.nH_4kHz=max(X1.nH_4kHz,X2.nH_4kHz);
    tmp.upper_cut_off_freq_of_interest=X1.upper_cut_off_freq_of_interest;
    tmp.upper_cut_off_freq=X1.upper_cut_off_freq;
    return tmp ;
}

/*-------------------------------------------------------------------*
* DTFS_fast_fs_inv()
*
* DTFS inverse.
*-------------------------------------------------------------------*/

void DTFS_fast_fs_inv(
    DTFS_STRUCTURE *X1_DTFS,  /* i : DTFS                     */
    float *out,      /* o : time domain output       */
    int   N          /* i : number of output samples */
)
{
    unsigned short i, M_2 = (unsigned short) min(X1_DTFS->lag>>1,X1_DTFS->nH), N_2 = (unsigned short) N>>1 ;
    float dbuf[256+1];     /* N can't be > 256 */

    if (N<X1_DTFS->lag)
    {
        N=X1_DTFS->lag;
    }

    /* Populate the dbuf array */
    dbuf[1] = X1_DTFS->a[0] ;
    dbuf[2] = 0.0 ;
    for (i=1; i<M_2; i++)
    {
        dbuf[2*i+1] =X1_DTFS->a[i] * N_2 ;
        dbuf[2*i+2] =X1_DTFS->b[i] * N_2 ;
    }

    if (N_2 != M_2)
    {
        dbuf[2*i+1] = X1_DTFS->a[i] * N_2 ;
        dbuf[2*i+2] = X1_DTFS->b[i] * N_2 ;
        i++;
    }

    /* Zero-padding in the frequency domain */
    for (   ; i<N_2; i++)
    {
        dbuf[2*i+1] = dbuf[2*i+2] = 0.0 ;
    }

    realft(dbuf, N_2, -1);

    for (i=1; i<=N; i++)
    {
        out[i-1] = dbuf[i]/N_2;
    }
}


/*-------------------------------------------------------------------*
* DTFS_freq_corr()
*
* Calculate correlation between two DTFS.
*-------------------------------------------------------------------*/

double DTFS_freq_corr(
    DTFS_STRUCTURE X1_DTFS,       /* i : X1 DTFS      */
    DTFS_STRUCTURE X2_DTFS,       /* i : X2 DTFS      */
    float lband,         /* i : low cutoff   */
    float hband          /* i : high cutoff  */
)
{
    short k ;
    double corr, fdiff, freq;
    double tmp;
    if(X1_DTFS.lag < X2_DTFS.lag)
    {
        DTFS_zeroPadd (X2_DTFS.lag,&X1_DTFS) ;
    }
    corr = freq = 0.0 ;
    fdiff = 12800.0/X2_DTFS.lag ;
    for(k=0; k<=min(X2_DTFS.lag>>1,X2_DTFS.nH_4kHz); k++, freq+=fdiff)
    {
        if(freq<hband &&  freq>=lband)
        {
            corr += (X1_DTFS.a[k]*X2_DTFS.a[k] + X1_DTFS.b[k]*X2_DTFS.b[k]);
        }
    }
    tmp = DTFS_getEngy_band (X1_DTFS,lband, hband)*DTFS_getEngy_band(X2_DTFS,lband, hband);
    if (tmp==0.0)
    {
        tmp=0.001;
    }
    return corr/sqrt(tmp) ;
}


/*-------------------------------------------------------------------*
* DTFS_alignment_weight()
*
* Estimate the shift to find the best match between the reference
* DTFS and the test DTFS.
*-------------------------------------------------------------------*/

float DTFS_alignment_weight(
    DTFS_STRUCTURE refX1_DTFS,      /* i : X1 the reference DTFS to keep fixed          */
    DTFS_STRUCTURE X2_DTFS,         /* i : X2 the test DTFS to shift to find best match */
    float Eshift,          /* i : Expected shift - coarse value                */
    const float *LPC1,           /* i : LPC to filter to find correlation in spch    */
    const float *LPC2            /* i : LPC to filter to find correlation in spch    */
)
{
    /* Eshift is w.r.t  X2 */
    short k ;
    float maxcorr, corr, Adiff, diff, tmp, tmp1, fshift, n, wcorr  ;
    float pwf = (float) 0.7, tmplpc[M+1] ;
    DTFS_STRUCTURE X1_DTFS;

    DTFS_copy (&X1_DTFS,refX1_DTFS);
    DTFS_adjustLag (&X1_DTFS,X2_DTFS.lag) ;
    DTFS_poleFilter (&X1_DTFS,LPC1,M+1) ;
    tmp=1.0;
    for(k=0,tmp=1.0; k<M+1; k++)
    {
        tmplpc[k]=LPC1[k]*(tmp*=pwf) ;
    }


    DTFS_zeroFilter (&X1_DTFS,tmplpc,M+1) ;
    DTFS_poleFilter (&X2_DTFS,LPC2,M+1) ;
    for(k=0,tmp=1.0; k<M+1; k++)
    {
        /* can be stored as a table */
        tmplpc[k]=LPC2[k]*(tmp*=pwf) ;
    }
    DTFS_zeroFilter (&X2_DTFS,tmplpc,M+1) ;
    maxcorr = (float) -HUGE_VAL;
    fshift = Eshift ;
    Adiff = max(6,0.15f*X2_DTFS.lag) ;
    if(X2_DTFS.lag < 60)
    {
        diff = 0.25;
    }
    else
    {
        diff = 0.5;
    }
    for(n=Eshift-Adiff; n<=Eshift+Adiff; n+=diff)
    {
        corr = tmp = 0.0f ;
        /* bit-exact optimization - PI2/X2_DTFS.lag should be counted as a single divide */
        tmp1 = (float) (PI2*n/X2_DTFS.lag) ;
        for(k=0; k<=min(X2_DTFS.lag>>1,X2_DTFS.nH_4kHz); k++, tmp+=tmp1)
        {
            /* Not counting math function cos and sin since they will be implemented as look-up tables */
            corr += (float) ((X1_DTFS.a[k]*X2_DTFS.a[k] + X1_DTFS.b[k]*X2_DTFS.b[k])* cos(tmp));
            corr += (float) ((X1_DTFS.b[k]*X2_DTFS.a[k] - X1_DTFS.a[k]*X2_DTFS.b[k])* sin(tmp));
        }
        wcorr=(float) (corr*(1.0f-0.01f* fabs(n-Eshift))) ;
        if( wcorr > maxcorr )
        {
            fshift = n ;
            maxcorr = wcorr ;
        }
    }
    return fshift ;  /* o : shift value to shift X2_DTFS by */
}


/*-------------------------------------------------------------------*
* DTFS_alignment_extract()
*
* Alignment for the best match between the reference DTFS and the test DTFS.
*-------------------------------------------------------------------*/

float DTFS_alignment_extract(
    DTFS_STRUCTURE refX1_DTFS,      /* i : X1 the reference DTFS to keep fixed          */
    DTFS_STRUCTURE X2_DTFS,         /* i : X2 the test DTFS to shift to find best match */
    float Eshift,          /* i : Expected shift - coarse value                */
    const float *LPC2            /* i : LPC to filter to find correlation in spch    */
)
{
    /* Eshift is w.r.t  X2 */
    short k ;
    float maxcorr, corr, Adiff, diff, tmp, tmp1, fshift, n ;
    float pwf = 0.7f, tmplpc[M+1] ;
    DTFS_STRUCTURE X1_DTFS;

    X1_DTFS = refX1_DTFS; /* copy into local copy */

    DTFS_adjustLag (&X1_DTFS,X2_DTFS.lag) ;

    DTFS_poleFilter (&X1_DTFS,LPC2,M+1) ;
    DTFS_poleFilter (&X2_DTFS,LPC2,M+1) ;

    for(k=0,tmp=1.0; k<M+1; k++)
    {
        tmplpc[k]=LPC2[k]*(tmp*=pwf) ;
    }
    DTFS_zeroFilter (&X1_DTFS,tmplpc,M+1) ;
    DTFS_zeroFilter (&X2_DTFS,tmplpc,M+1) ;

    maxcorr = (float) -HUGE_VAL;
    fshift = Eshift ;
    Adiff = max(4.0f,refX1_DTFS.lag/8);
    diff = 1.0f; /* Non-Fractional alignment */

    for(n=Eshift-Adiff; n<=Eshift+Adiff; n+=diff)
    {
        corr = tmp = 0.0f ;
        /* bit-exact optimization - PI2/X2_DTFS.lag should be counted as a single divide outside WI functions and passed in as input */
        tmp1 = (float) (PI2*n/X2_DTFS.lag) ;

        for(k=0; k<=min(X2_DTFS.lag>>1,X2_DTFS.nH_4kHz); k++, tmp+=tmp1)

        {
            corr += (float) ((X1_DTFS.a[k]*X2_DTFS.a[k] + X1_DTFS.b[k]*X2_DTFS.b[k])* cos(tmp));
            corr += (float) ((X1_DTFS.b[k]*X2_DTFS.a[k] - X1_DTFS.a[k]*X2_DTFS.b[k])* sin(tmp));
        }

        if( corr*(1.0f-0.01f*fabs(n-Eshift) ) > maxcorr )
        {
            fshift = n ;
            maxcorr = corr ;
        }
    }

    return fshift ;   /* o : shift value to shift X2 by   */
}



/*-------------------------------------------------------------------*
* DTFS_alignment_fine_new()
*
*  Shift value for DTFS finer alignment.
*-------------------------------------------------------------------*/

float DTFS_alignment_fine_new(
    DTFS_STRUCTURE X1_DTFS,    /* i : X1 the reference DTFS to keep fixed          */
    DTFS_STRUCTURE X2_DTFS,    /* i : X2 the test DTFS to shift to find best match */
    float Eshift      /* i : Expected shift - coarse value                */
)
{
    short k ;
    float maxcorr, corr, Adiff, diff, tmp, tmp1, fshift, n ;

    if (X1_DTFS.lag < X2_DTFS.lag)
    {
        DTFS_zeroPadd (X2_DTFS.lag,&X1_DTFS) ;
    }
    maxcorr = (float) -HUGE_VAL;
    fshift = Eshift ;
    Adiff = 20.0f ;
    diff =  1.0f;
    for(n=Eshift-Adiff+1; n<=Eshift+Adiff; n+=diff)
    {
        corr = tmp = 0.0 ;
        /* bit-exact optimization - PI2/X2_DTFS.lag should be counted as a single divide outside loops */
        tmp1 = (float) (PI2*n/X2_DTFS.lag) ;

        for(k=0; k<=min(X2_DTFS.lag>>1,X2_DTFS.nH); k++, tmp+=tmp1)
        {
            corr += (float) ((X1_DTFS.a[k]*X2_DTFS.a[k] + X1_DTFS.b[k]*X2_DTFS.b[k])* cos(tmp));
            corr += (float) ((X1_DTFS.b[k]*X2_DTFS.a[k] - X1_DTFS.a[k]*X2_DTFS.b[k])* sin(tmp));
        }
        if( corr*(1.0f-0.01f*fabs(n-Eshift) ) > maxcorr )
        {
            fshift = n ;
            maxcorr = corr ;
        }
    }
    return fshift ; /* o : shift value to shift X2 by    */
}


/*-------------------------------------------------------------------*
* DTFS_alignment_full()
*
*  Shift value for DTFS full alignment.
*-------------------------------------------------------------------*/

float DTFS_alignment_full(
    DTFS_STRUCTURE X1_DTFS,           /* i : reference DTFS     */
    DTFS_STRUCTURE X2_DTFS,           /* i : DTFS to shift      */
    int num_steps          /* i : resolution         */
)
{
    short k ;
    float maxcorr, corr, tmp, tmp1, fshift, n, diff ;

    if (X1_DTFS.lag < X2_DTFS.lag)
    {
        DTFS_zeroPadd (X2_DTFS.lag,&X1_DTFS) ;
    }

    maxcorr = (float) -HUGE_VAL;
    /* bit-exact optimization - 1/num_steps can be constant => should be counted as a multiply */
    diff = (float) X2_DTFS.lag / num_steps ;

    for(fshift=n=0.0; n<(float)X2_DTFS.lag; n+=diff)
    {
        corr = tmp = 0.0f ;
        tmp1 = (float) (PI2*n/X2_DTFS.lag) ;

        for(k=0; k<=min(X2_DTFS.lag>>1,X2_DTFS.nH_4kHz); k++, tmp+=tmp1)

        {
            corr += (float) ((X1_DTFS.a[k]*X2_DTFS.a[k] + X1_DTFS.b[k]*X2_DTFS.b[k])* cos(tmp));
            corr += (float) ((X1_DTFS.b[k]*X2_DTFS.a[k] - X1_DTFS.a[k]*X2_DTFS.b[k])* sin(tmp));
        }
        if(corr>maxcorr)
        {
            fshift = n ;
            maxcorr = corr ;
        }
    }
    return fshift ;
}


/*-------------------------------------------------------------------*
* DTFS_phaseShift()
*
*  Phase shift the DTFS coefficients.
*  ph is the amount of phase shift (between 0 and 2pi), where positive
*  value indicates a shift to right, and a negative ph value indicates a
*  left shift.
*-------------------------------------------------------------------*/

void DTFS_phaseShift(
    DTFS_STRUCTURE *X,      /* i/o: DTFS to shift    */
    float ph       /* i  : phase to shift   */
)
{
    short k ;
    float tmp, tmp2=0.0f ;
    for(k=0; k<=min(X->lag>>1,X->nH); k++, tmp2+=ph)
    {
        tmp = X->a[k] ;
        X->a[k] = (float) (tmp*cos(tmp2) - X->b[k]*sin(tmp2)) ;
        X->b[k] = (float) (tmp*sin(tmp2) + X->b[k]*cos(tmp2)) ;
    }
}


/*-------------------------------------------------------------------*
 * DTFS_zeroPadd()
 *
 *  Zero-pad the DTFS coefficients.
 *-------------------------------------------------------------------*/

void DTFS_zeroPadd(
    int   N,        /* i  : Target lag    */
    DTFS_STRUCTURE *X        /* i/o: DTFS          */
)
{
    int i ;
    float diff;
    if(N == X->lag)
    {
        return ;
    }
    for(i=(X->lag>>1)+1; i<=N>>1; i++)
    {
        X->a[i]=X->b[i]=0.0 ;
    }
    X->lag = N ;
    /* recompute nH for new lag */
    X->nH=(int)floor(X->upper_cut_off_freq/(12800.0/X->lag));
    diff = 12800.0f / X->lag ;
    if (X->upper_cut_off_freq-(diff*X->nH)>=diff)
    {
        X->nH++;
    }
}


/*-------------------------------------------------------------------*
* DTFS_to_fs()
*
* DTFS to fs conversion.
*-------------------------------------------------------------------*/

void DTFS_to_fs(
    const float *x,       /* i : time domain signal               */
    int   N,        /* i : Length of input vector           */
    DTFS_STRUCTURE *X,       /* o : DTFS structure with a, b, lag    */
    short Fs,       /* i : sampling rate                    */
    short FR_flag   /* i :  FR flag                         */
)
{
    short n;
    int nH, k, nH_band, nH_4kHz;
    float sum, tmp, diff;

    if (!FR_flag)
    {
        if (Fs==16000)
        {
            X->upper_cut_off_freq_of_interest=4000.0;
            X->upper_cut_off_freq=6400.0;
            X->Fs=INT_FS_12k8;
        }
        else if (Fs==8000)
        {
            X->upper_cut_off_freq_of_interest=3300.0;
            X->upper_cut_off_freq=4000.0;
            X->Fs=INT_FS_12k8;
        }
    }
    else
    {
        X->upper_cut_off_freq_of_interest=8000.0;
        X->upper_cut_off_freq=8000.0;
        X->Fs=16000.0;
    }

    X->lag = N ;
    nH_band=(int)floor(X->upper_cut_off_freq/(12800.0/X->lag));

    nH_4kHz=(int)floor(4000/(12800.0/X->lag));
    diff = 12800.0f / X->lag ;
    if (X->upper_cut_off_freq-(diff*nH_band)>=diff)
    {
        nH_band++;
    }
    if (4000-(diff*nH_4kHz)>=diff)
    {
        nH_4kHz++;
    }
    /* Number of harmonics excluding the ones at 0 and at pi */
    nH = (N-1) >> 1;
    /* The DC component */
    X->a[0] = 0.0 ;
    X->b[0] = 0.0 ;
    for( n=0; n<N; n++ )
    {
        X->a[0] += x[n] ;
    }
    X->a[0] /= N ;

    /* Strictly set the DC componet to zero */
    X->a[0] = 0.0 ;

    /* The harmonics excluding the one at pi */
    for( k=1; k<=nH; k++ )
    {
        X->a[k] = x[0] ;
        X->b[k] = 0.0 ;
        sum = tmp = (float) (PI2*k/N) ;
        for( n=1; n<N; n++, sum+=tmp )
        {
            X->a[k] += (float) (x[n] * cos(sum)) ;
            X->b[k] += (float) (x[n] * sin(sum)) ;
        }
        X->a[k] *= (2.0f/N) ;
        X->b[k] *= (2.0f/N) ;
    }

    /* The harmonic at 'pi' */
    if( N%2 == 0 )
    {
        X->a[k] = 0.0 ;
        tmp = 1.0 ;
        for( n=0; n<N; n++, tmp*=-1.0 )
        {
            X->a[k] += x[n] * tmp  ;
        }
        X->a[k] /= N ;
        X->b[k] = 0.0 ;
    }
    for(k=nH_band+1; k<=min((X->lag>>1),(MAXLAG_WI-1)); k++)
    {
        X->a[k]=0.0;
        X->b[k]=0.0;
    }
    X->nH=nH_band;
    X->nH_4kHz=nH_4kHz;
}


/*-------------------------------------------------------------------*
* DTFS_fs_inv()
*
*  DTFS Inverse
*-------------------------------------------------------------------*/

void DTFS_fs_inv(
    DTFS_STRUCTURE *X,      /* i : DTFS input                         */
    float *x,      /* o : time domain sig                    */
    int   N,       /* i : Output length                      */
    float ph0      /* i : phase shift applied to the output  */
)
{
    float phase, tmp ;
    short k, n ;

    for( n=0; n<N; n++ )
    {
        x[n] = X->a[0] ;
        tmp = phase = (float) (PI2*n/X->lag + ph0) ;
        for( k=1; k <=min(X->lag>>1,X->nH); k++, tmp+=phase )
        {
            x[n] += (float) (X->a[k]*cos(tmp) + X->b[k]*sin(tmp)) ;
        }
    }
}


/*-------------------------------------------------------------------*
* DTFS_transform()
*
*  DTFS transform.
*-------------------------------------------------------------------*/

void DTFS_transform(
    DTFS_STRUCTURE X,        /* i : Starting DTFS to use in WI    */
    DTFS_STRUCTURE X2,       /* i : Ending DTFS to use in WI      */
    const float *phase,   /* i : Phase contour                 */
    float *out,     /* o : Output time domain waveform   */
    int   N,        /* i : Number of samples to generate */
    short FR_flag   /* i : Flag to indicate called in FR context */
)
{
    short i, j,j1;
    float w, tmp ;
    float x1_256[256], x2_256[256];
    float sum1, sum2;
    short m, l1, k;
    int N1;
    float nrg_flag = 0;
    float x_r_fx[L_FRAME];
    float temp_w;

    DTFS_STRUCTURE *tmp1_dtfs=DTFS_new();
    DTFS_STRUCTURE *tmp2_dtfs=DTFS_new();
    DTFS_STRUCTURE *tmp3_dtfs=DTFS_new();
    DTFS_copy (tmp1_dtfs,X);
    DTFS_copy (tmp2_dtfs,X2);
    DTFS_fast_fs_inv (tmp1_dtfs,x1_256,256);
    DTFS_fast_fs_inv (tmp2_dtfs,x2_256,256);
    tmp = (float) (log(1.0 - WI_THRESHLD) / (X.lag - N)) ;
    for(i=0; i<N; i++)
    {
        if (FR_flag==0)
        {
            /* should not be counted inside the loop */
            if ( N - WI_SAMPLE_THLD > X.lag )
            {
                /* pre-computed and stored in a table */
                w = (float) (1.0 - exp (- (i+1) * tmp)) ;
            }
            else
            {
                /* can be a look-up table */
                w = (float) (i+1) / N ;
            }
        }
        else
        {
            if (nrg_flag)
            {
                w = (float) (i+1)/N;
            }
            else
            {
                if (N <= tmp2_dtfs->lag)
                {
                    N=tmp2_dtfs->lag+1;
                }

                N1 = N - tmp2_dtfs->lag;
                if ( i < N1 )
                {
                    w = (float) (i+1)/N1;
                }
                else w = 1.0;
            }
        }

        /* add sinc interpolation of two time domain waveforms at
        appropriate phase position */
        j = (LL_OS*10 + (int)rint_new (phase[i]*LL_OS/PI2)) % LL_OS;

        if (j<0)
        {
            j=0;
        }

        k=j%WARP_OS_RATE;
        l1=j/WARP_OS_RATE;

        set_f( x_r_fx, 0.0f, L_FRAME );

        temp_w =( 1 - w );

        for( j1=0; j1<12; j1++ )
        {
            m=(1000*LL+l1-OSLENGTH/2+j1)%LL;

            if( m < 0 )
            {
                m=0;
            }

            x_r_fx[m] =  x1_256[m]* temp_w +  x2_256[m] * w;
        }

        for(j1=0,sum1=sum2=0.0; j1<OSLENGTH; j1++)
        {
            /* mult or div by constants should be done once outside the loop */
            m=(1000*LL+l1-OSLENGTH/2+j1)%LL;

            if (m<0)
            {
                m=0;
            }

            sum1 += x_r_fx[m]*sinc[k][j1];
        }

        out[i]  = sum1;
    }



    free(tmp1_dtfs);
    free(tmp2_dtfs);
    free(tmp3_dtfs);
}


/*-------------------------------------------------------------------*
* DTFS_zeroFilter()
*
*  DTFS zero filter response.
*-------------------------------------------------------------------*/

void DTFS_zeroFilter(
    DTFS_STRUCTURE *X,      /* i/o: DTFS to zeroFilter inplace  */
    const float *LPC,    /* i  : LPCs                        */
    int   N        /* i  : LPC order                   */
)
{
    float tmp, tmp1, tmp2, sum1, sum2 ;
    short k, n ;
    tmp1 = (float) (PI2/X->lag) ;
    for( k=0 ; k<=min(X->lag>>1,X->nH) ; k++ )
    {
        tmp = tmp2 = k*tmp1 ;
        /* Calculate sum1 and sum2 */
        sum1 = 1.0;
        sum2 = 0.0;
        for( n=0 ; n<N ; n++, tmp2+=tmp )
        {
            sum1 += (float) (LPC[n] * cos(tmp2)) ;
            sum2 += (float) (LPC[n] * sin(tmp2)) ;
        }
        /* Calculate the circular convolution */
        tmp = X->a[k] ;
        X->a[k] = tmp * sum1 - X->b[k] * sum2 ;
        X->b[k] = X->b[k] * sum1 + tmp * sum2 ;
    }
}

/*-------------------------------------------------------------------*
* DTFS_poleFilter()
*
*  DTFS pole filter response.
*-------------------------------------------------------------------*/

void DTFS_poleFilter(
    DTFS_STRUCTURE *X,       /* i/o : DTFS to poleFilter inplace */
    const float *LPC,     /* i : LPCs                         */
    int   N         /* i : LPC order                    */
)
{
    float tmp, tmp1, tmp2, sum1, sum2 ;
    short k, n ;
    tmp1 = (float) (PI2/X->lag) ;
    for( k=0; k<=min(X->lag>>1,X->nH) ; k++ )
    {
        tmp = tmp2 = k*tmp1 ;
        /* Calculate sum1 and sum2 */
        sum1 = 1.0;
        sum2 = 0.0;
        for( n=0 ; n<N ; n++, tmp2+=tmp)
        {
            sum1 += (float) (LPC[n] * cos(tmp2)) ;
            sum2 += (float) (LPC[n] * sin(tmp2)) ;
        }
        /* Calculate the circular convolution */
        tmp = X->a[k] ;
        tmp2 = sum1 * sum1 + sum2 * sum2 ;
        X->a[k] = ( tmp * sum1 + X->b[k] * sum2 ) / tmp2 ;
        X->b[k] = ( -tmp * sum2 + X->b[k] * sum1 ) / tmp2 ;
    }
}


/*-------------------------------------------------------------------*
* DTFS_adjustLag()
*
*  Adjust DTFS lag based on a target lag.
*-------------------------------------------------------------------*/

void DTFS_adjustLag(
    DTFS_STRUCTURE *X_DTFS,      /* i/o : DTFS to adjust lag for */
    int N             /* i : Target lag               */
)
{
    int k ;
    float en,diff ;
    if (N==X_DTFS->lag)
    {
        return ;
    }

    if(N>X_DTFS->lag)
    {

        DTFS_zeroPadd (N, X_DTFS) ;
    }
    else
    {
        en = DTFS_getEngy (*X_DTFS) ;
        for(k=(N>>1)+1; k<=min(X_DTFS->lag>>1,X_DTFS->nH); k++)
        {
            X_DTFS->a[k] = 0.0 ;
            X_DTFS->b[k] = 0.0 ;
        }
        DTFS_setEngy (X_DTFS, en) ;
        X_DTFS->lag = N ;
        /* recompute nH for new lag */
        X_DTFS->nH=(int)floor(X_DTFS->upper_cut_off_freq/(12800.0/X_DTFS->lag));

        X_DTFS->nH_4kHz=(int)floor(4000.0/(12800.0/X_DTFS->lag));
        diff = 12800.0f / X_DTFS->lag ;
        if (X_DTFS->upper_cut_off_freq-(diff*X_DTFS->nH)>=diff)
        {
            X_DTFS->nH++;
        }
        if (4000.0-(diff*X_DTFS->nH_4kHz)>=diff)
        {
            X_DTFS->nH_4kHz++;
        }
    }
}


/*-------------------------------------------------------------------*
* DTFS_getEngy()
*
*  Get DTFS energy.
*-------------------------------------------------------------------*/

float DTFS_getEngy(
    DTFS_STRUCTURE X        /* i : DTFS to compute energy of */
)
{
    short k;
    float en=0.0 ;
    for(k=1; k<=min((X.lag-1)>>1,X.nH); k++)
    {
        en += X.a[k]*X.a[k] + X.b[k]*X.b[k] ;
    }
    en /= 2.0 ;
    en += X.a[0]*X.a[0] ;
    if(X.lag%2 == 0)
    {
        en += X.a[k]*X.a[k] + X.b[k]*X.b[k];
    }
    return en ;
}


/*-------------------------------------------------------------------*
 * DTFS_getEngy_band()
 *
 *  Get DTFS energy in the specified range from lband to hband.
 *-------------------------------------------------------------------*/

float DTFS_getEngy_band(
    DTFS_STRUCTURE X,          /* i : DTFS to compute energy of    */
    float lband,      /* i : low end of band of interest  */
    float hband       /* i : high end of band of interest */
)
{
    short k;
    float en=0.0f, freq, fdiff=12800.0f/X.lag ;
    for(freq=fdiff, k=1; k<=min((X.lag-1)>>1,X.nH_4kHz); k++, freq+=fdiff)
    {
        if(freq<hband &&  freq>=lband)
        {
            en += X.a[k]*X.a[k] + X.b[k]*X.b[k] ;
        }
    }
    en /= 2.0 ;
    if(lband == 0.0)
    {
        en += X.a[0]*X.a[0] ;
    }
    if((X.lag%2 == 0) &&  (hband == X.upper_cut_off_freq))
    {
        en += X.a[k]*X.a[k] + X.b[k]*X.b[k];
    }
    return en ;
}

/*-------------------------------------------------------------------*
 * DTFS_getEngy_band_wb()
 *
 *  Get DTFS energy in the specified range from lband to hband.
 *  This function is different to "DTFS_getEngy_band" as this can calculate
 *  lband, hband \in [1,6400] where "DTFS_getEngy_band" only upperlimited to
 *	 4Khz. Possibility: modify ""DTFS_getEngy_band"" and get rid of this
 *  function.
 *-------------------------------------------------------------------*/

float DTFS_getEngy_band_wb(
    DTFS_STRUCTURE X,          /* i : DTFS to compute energy of    */
    float lband,      /* i : low end of band of interest  */
    float hband       /* i : high end of band of interest */
)
{
    short k;
    float en=0.0f, freq, fdiff=12800.0f/X.lag ;
    for(freq=fdiff, k=1; k<=((X.lag-1)>>1); k++, freq+=fdiff)
    {
        if(freq<hband &&  freq>=lband)
        {
            en += X.a[k]*X.a[k] + X.b[k]*X.b[k] ;
        }
    }
    en /= 2.0 ;
    if(lband == 0.0)
    {
        en += X.a[0]*X.a[0] ;
    }
    if((X.lag%2 == 0) &&  (hband == X.upper_cut_off_freq))
    {
        en += X.a[k]*X.a[k] + X.b[k]*X.b[k];
    }
    return en ;
}


/*-------------------------------------------------------------------*
 * DTFS_setEngy()
 *
 *  Set DTFS energy.
 *-------------------------------------------------------------------*/

float DTFS_setEngy(
    DTFS_STRUCTURE *X_DTFS,      /* i/o : DTFS structure to set engy */
    float en2           /* i : Energy to set to             */
)
{
    short k ;
    float en1, tmp ;
    en1 = DTFS_getEngy (*X_DTFS);
    if(en1 ==  0.0)
    {
        return 0.0;
    }
    tmp = (float) sqrt(en2/en1) ;
    for(k=0; k<=min(X_DTFS->lag>>1,X_DTFS->nH); k++)
    {
        X_DTFS->a[k] *= tmp ;
        X_DTFS->b[k] *= tmp ;
    }
    return en1 ;
}


/*-------------------------------------------------------------------*
 * DTFS_car2pol()
 *
 *  DTFS cartesian to polar co-ordinates conversion.
 *-------------------------------------------------------------------*/

void DTFS_car2pol(
    DTFS_STRUCTURE *X       /* i/o : DTFS structure a, b, lag  */
)
{
    short k ;
    float tmp ;
    for( k=1 ; k<=min((X->lag-1)>>1,X->nH) ; k++ )
    {
        tmp = X->a[k] ;
        X->a[k] = (float) (0.5f * sqrt( tmp*tmp + X->b[k]*X->b[k] )) ;
        X->b[k] = (float) atan2(X->b[k], tmp) ;
    }
    if(X->lag%2==0)
    {
        tmp = X->a[k] ;
        X->a[k] = (float) sqrt( tmp*tmp + X->b[k]*X->b[k] ) ;
        X->b[k] = (float) atan2 (X->b[k], tmp) ;
    }
}


/*-------------------------------------------------------------------*
 * DTFS_pol2car()
 *
 *  DTFS polar to cartesian co-ordinates conversion.
 *-------------------------------------------------------------------*/

void DTFS_pol2car(
    DTFS_STRUCTURE *X       /* i/o : DTFS structure a, b, lag  */
)
{
    short k ;
    float tmp ;
    for( k=1 ; k<=min((X->lag-1)>>1,X->nH) ; k++ )
    {
        tmp = X->b[k] ;
        X->b[k] = (float) (2.0f * X->a[k] * sin(tmp)) ;
        X->a[k] = (float) (2.0f * X->a[k] * cos(tmp)) ;
    }
    if(X->lag%2==0)
    {
        tmp = X->b[k] ;
        X->b[k] = (float) (X->a[k] * sin(tmp)) ;
        X->a[k] = (float) (X->a[k] * cos(tmp)) ;
    }
}


/*-------------------------------------------------------------------*
* DTFS_setEngyHarm()
*
*  Adjust DTFS energy based on the input target energy.
*-------------------------------------------------------------------*/

float DTFS_setEngyHarm(
    float f1,             /* i  : lower band freq of input to control energy   */
    float f2,             /* i  : upper band freq of input to control energy   */
    float g1,             /* i  : lower band freq of output to control energy  */
    float g2,             /* i  : upper band freq of output to control energy  */
    float en2,            /* i  : Target Energy to set the DTFS to             */
    DTFS_STRUCTURE *X              /* i/o: DTFS to adjust the energy of                 */
)
{
    short k, count=0;
    float en1=0.0f, tmp, factor, diff=12800.0f/X->lag ;

    if(f1==0.0)
    {
        en1 += X->a[0]*X->a[0] ;
        count ++;
    }
    for(k=1, tmp=diff; k<=min((X->lag-1)>>1,X->nH); k++, tmp+=diff)
    {
        if( X->a[k] < EPSILON )
        {
            X->a[k]=0;
        }

        if(tmp>f1 &&  tmp<=f2)
        {
            en1 += X->a[k]*X->a[k] ;
            count ++;
        }
    }

    if (count <= 0.0)
    {
        count=1;
    }


    en1 /= count ;

    if (en2 < 0.0)
    {
        en2=0.0;
    }

    if(en1>0.0)
    {
        factor = (float) sqrt(en2/en1) ;
    }
    else
    {
        factor = 0.0f ;
    }
    if(g1==0.0)
    {
        X->a[k] *= factor ;
    }
    for(k=1, tmp=diff; k<=min((X->lag-1)>>1,X->nH); k++, tmp+=diff)
    {
        if(tmp>g1 &&  tmp<=g2)
        {
            X->a[k]*=factor ;
        }
    }
    return (float)(en1+1e-20);  /* o : Return Input RMS between f1 and f2 before scaling  */
}


/*-------------------------------------------------------------------*
 * cubicPhase()
 *
 * Compute coefficients of cubic phase function
 *-------------------------------------------------------------------*/
static
void cubicPhase(
    float ph1,        /* i  : phase offset       */
    float ph2,        /* i  : phase 2            */
    float L1,         /* i  : previous lag       */
    float L2,         /* i  : current lag        */
    int   N,          /* i  : input length       */
    float *phOut      /* o  : cubic phase output */
)
{
    float coef[4], f1, f2, c1, c2, factor ;
    short n ;
    double diff;

    N -= (int) L2;

    if (N<=0)
    {
        N=1;
    }

    /* Computation of the coefficients of the cubic phase function */
    f1 = (float) (PI2 / L1) ;
    f2 = (float) (PI2 / L2) ;
    ph1 = (float) fmod ((double)(ph1), PI2) ;
    ph2 = (float) fmod ((double)(ph2), PI2) ;
    coef[3] = ph1 ;
    coef[2] = f1 ;
    factor = (float) (anint (( ph1 - ph2 + 0.5*N*(f2+f1) ) / PI2 )) ;
    c1 = f2-f1 ;
    c2 = (float) (ph2 - ph1 - N*f1 + PI2 * factor) ;
    coef[0] = (N*c1 - 2*c2 ) / (N*N*N) ;
    coef[1] = (c1 - 3*N*N*coef[0]) / (2*N) ;
    /* Computation of the phase value at each sample point */
    phOut[0] = ph1 ;
    for(n=1; n<N; n++)
    {
        phOut[n] = _POLY3(n,coef) ;
    }
    diff = (float) (PI2/L2);
    for(   ; n<N+(int)L2; n++)
    {
        phOut[n] = (float) (phOut[n-1]+diff);
    }
}


/*-------------------------------------------------------------------*
 * DTFS_to_erb()
 *
 * DTFS to ERB conversion
 *-------------------------------------------------------------------*/

void DTFS_to_erb(
    DTFS_STRUCTURE X,        /* i : DTFS input   */
    float *out      /* o : ERB output       */
)
{
    short num_erb;
    unsigned short i, j, count[NUM_ERB_WB] ;
    float freq, diff ;

    const float *erb = NULL;
    num_erb=NUM_ERB_NB;
    if (X.upper_cut_off_freq==4000.0)
    {
        num_erb=NUM_ERB_NB;
        erb=&(erb_NB[0]);
    }
    else if (X.upper_cut_off_freq==6400.0)
    {
        num_erb=NUM_ERB_WB;
        erb=&(erb_WB[0]);
    }

    for(i=0; i<num_erb; i++)
    {
        count[i] = 0 ;
        out[i] = 0.0 ;
    }
    diff = 12800.0f / X.lag ;
    for(i=j=0, freq=0.0; i<=min(X.lag>>1,X.nH); i++, freq+=diff)
    {
        if (!(freq <= erb[num_erb]))
        {
            freq = erb[num_erb];
        }

        for(   ; j<num_erb; j++)
        {
            if(freq < erb[j+1])
            {
                if (X.a[i] < 0.0f)
                {
                    X.a[i]=0.0f;
                }

                out[j] += X.a[i] ;
                count[j]++ ;
                break ;
            }
        }
    }
    for(i=0; i<num_erb; i++)
    {
        if(count[i] > 1)
        {
            out[i] /= count[i] ;
        }
    }
}

/*-------------------------------------------------------------------*
 * erb_slot()
 *
 * Estimation of ERB slots.
 *-------------------------------------------------------------------*/
static
void erb_slot(
    int   lag,        /* i : input lag          */
    int   *out,       /* o : ERB slots          */
    float *mfreq,     /* i : ERB frequencies    */
    short num_erb     /* i : number of ERBs     */
)
{
    unsigned short i, j ;
    float freq, diff ;
    short upper_cut_off_freq;

    const float *erb=NULL;
    int nH_band;
    upper_cut_off_freq=4000;
    if (num_erb==NUM_ERB_NB)
    {
        upper_cut_off_freq=4000;
        erb=&(erb_NB[0]);
    }
    else if (num_erb==NUM_ERB_WB)
    {
        upper_cut_off_freq=6400;
        erb=&(erb_WB[0]);
    }
    nH_band=(int)floor(upper_cut_off_freq/(12800.0/lag));

    for(i=0; i<num_erb; i++)
    {
        out[i] = 0 ;
        mfreq[i] = 0.0 ;
    }
    diff = 12800.0f / lag ;
    if (upper_cut_off_freq-(diff*nH_band)>=diff)
    {
        nH_band++;
    }
    for(i=j=0, freq=0.0; i<=min(lag>>1,nH_band); i++, freq+=diff)
    {

        if (!(freq <= erb[num_erb]))
        {
            freq = erb[num_erb];
        }

        freq = min(freq, upper_cut_off_freq) ;

        for(   ; j<num_erb; j++)
        {
            if(freq < erb[j+1])
            {
                mfreq[j] += freq ;
                out[j]++ ;
                break ;
            }
        }
    }
    for(j=0; j<num_erb; j++)
    {
        if(out[j]>1)
        {
            mfreq[j]/=out[j] ;
        }
    }
}

/*-------------------------------------------------------------------*
 * DTFS_erb_inv()
 *
 * .DTFS after ERB inverse
 *-------------------------------------------------------------------*/

void DTFS_erb_inv(
    float *in,          /* i : ERB inpt                      */
    int   *slot,        /* i : ERB slots filled based on lag */
    float *mfreq,       /* i : erb frequence edges           */
    DTFS_STRUCTURE *X,           /* o : DTFS after erb-inv           */
    short num_erb       /* i : Number of ERB bands           */
)
{
    short i, j, m=0 ;
    float diff ;
    float freq, f[NUM_ERB_WB+2], amp[NUM_ERB_WB+2] ;
    short upper_cut_off_freq = 0;

    const float *erb = NULL;

    if (num_erb==NUM_ERB_NB)
    {
        upper_cut_off_freq=4000;
        erb=&(erb_NB[0]);
    }
    else if (num_erb==NUM_ERB_WB)
    {
        upper_cut_off_freq=6400;
        erb=&(erb_WB[0]);
    }

    f[m]=0.0;
    amp[m]=0.0;
    m++;
    for(i=0; i<num_erb; i++)
    {
        if(slot[i] != 0)
        {
            f[m]=mfreq[i];
            amp[m]=in[i];
            m++;
        }
    }
    f[m]=upper_cut_off_freq;
    amp[m]=0.0;
    m++;

    diff = 12800.0f / X->lag ;

    for(i=0, j=1, freq=0.0; i<=min(X->lag>>1,X->nH); i++, freq+=diff)
    {
        if (!(freq <= erb[num_erb]))
        {
            freq = erb[num_erb];
        }

        if (!(m<=num_erb+2))
        {
            m = num_erb+2;
        }

        if(freq>upper_cut_off_freq)
        {
            freq=upper_cut_off_freq;
        }
        for(   ; j<m; j++)
        {
            if(freq <= f[j])
            {
                X->a[i] = amp[j]*(freq-f[j-1]) + amp[j-1]*(f[j]-freq) ;
                if(f[j] != f[j-1])
                {
                    X->a[i] /= (f[j] - f[j-1]) ;
                }
                break ;
            }
        }

        X->a[0] = 0.0f;
    }
}


/*-------------------------------------------------------------------*
 * LPCPowSpect()
 *
 *  LPC power spectrum
 *-------------------------------------------------------------------*/
static
void LPCPowSpect(
    const float *freq,        /* i : ERB frequencies    */
    int   Nf,           /* i : Number of ERBs     */
    const float *LPC,         /* i : LPC coefficients   */
    int   Np,           /* i : Number of LPCs     */
    float *out          /* o : LPC power spectrum */
)
{
    float w, tmp, Re, Im ;
    short i, k ;

    for(k=0; k<Nf; k++)
    {
        Re = 1.0 ;
        Im = 0.0 ;
        /* Note that freq ranges between [0 UPPER_CUT_OFF_FREQ] */
        tmp = (float) (freq[k]/12800.0f * PI2) ;
        for(i=0, w=tmp; i<Np; i++, w+=tmp)
        {
            Re += (float) (LPC[i] * cos(w)) ;
            Im -= (float) (LPC[i] * sin(w)) ;
        }
        out[k] = 1.0f/(Re*Re + Im*Im) ;
    }
}


/*-------------------------------------------------------------------*
 * DTFS_getSpEngyFromResAmp()
 *
 *  Get speech energy from the DTFS
 *-------------------------------------------------------------------*/

float DTFS_getSpEngyFromResAmp(
    DTFS_STRUCTURE X,          /* i : DTFS                             */
    float lband,      /* i : Low band end to get energy from  */
    float hband,      /* i : High band end to get energy from */
    const float *curr_lsp   /* i : LPCs                             */
)
{
    short i, k;
    float w, tmp, Re, Im ;
    double en=0.0, freq, fdiff;
    fdiff=12800.0/X.lag ;
    if(hband == X.upper_cut_off_freq)
    {
        hband = 4001.0;
    }
    for(freq=0.0, k=0; k<=min(X.lag>>1,X.nH_4kHz); k++, freq+=fdiff)
    {
        if (X.a[k]<0.0)
        {
            X.a[k]=0.0;
        }

        if(freq<hband &&  freq>=lband)
        {
            Re = 1.0f ;
            Im = 0.0f ;
            tmp = (float) (PI2*freq/12800.0f) ;
            for(i=0, w=tmp; i<M+1; i++, w+=tmp)
            {
                Re += (float) (curr_lsp[i] * cos(w)) ;
                Im -= (float) (curr_lsp[i] * sin(w)) ;
            }
            if(k==0 ||  (X.lag%2==0 &&  k==X.lag>>1))
            {
                en += X.a[k]*X.a[k] / (Re*Re + Im*Im) ;
            }
            else
            {
                en += 2.0 * X.a[k]*X.a[k] / (Re*Re + Im*Im) ;
            }
        }
    }
    return ((float) en);
}


/*-------------------------------------------------------------------*
 * erb_diff()
 *
 * ERB difference
 *-------------------------------------------------------------------*/
static
void erb_diff(
    const float *prev_erb,       /* i  : previous ERB              */
    int   pl,              /* i  : previous lag              */
    const float *curr_erb,       /* i  : current ERB               */
    int   l,               /* i  : current lag               */
    const float *curr_lsp,       /* i  : current LSP coefficients  */
    float *out,            /* o  : ERB difference            */
    int   *index,          /* i  : ERB index                 */
    short num_erb          /* i  : Number of ERBs            */
)
{
    short i, j;
    int pslot[NUM_ERB_WB], cslot[NUM_ERB_WB], mmseindex ;
    float tmp, tmp1, t_prev_erb[NUM_ERB_WB], LPC[M+1], mfreq[NUM_ERB_WB], PowSpect[NUM_ERB_WB], mmse ;
    const float (*AmpCB1)[10] = 0;

    if (num_erb==NUM_ERB_NB)
    {
        AmpCB1=AmpCB1_NB;
    }
    else if (num_erb==NUM_ERB_WB)
    {
        AmpCB1=AmpCB1_WB;

    }

    erb_slot (l,cslot,mfreq,num_erb);
    erb_slot (pl,pslot,t_prev_erb,num_erb);
    for(i=0,tmp=1.0f; i<M+1; i++)
    {
        LPC[i]=curr_lsp[i]*(tmp*=0.78f) ;
    }
    LPCPowSpect (mfreq, num_erb, LPC, M+1, PowSpect);

    for(i=0; i<num_erb; i++)
    {
        t_prev_erb[i] = prev_erb[i] ;
    }

    if(pl>l)
    {
        tmp = t_prev_erb[0] ;
        for (i=0; i<num_erb; i++)
        {
            if (pslot[i] < 0)
            {
                pslot[i]=0;
            }

            if (pslot[i]!=0)
            {
                tmp = t_prev_erb[i] ;
            }
            else
            {
                t_prev_erb[i] = tmp ;
            }
        }
    }
    else if (l>pl)
    {
        tmp = t_prev_erb[num_erb-1] ;
        for(i=num_erb-1; i>=0; i--)
        {
            if (pslot[i]!=0)
            {
                tmp = t_prev_erb[i];
            }
            else
            {
                t_prev_erb[i] = tmp ;
            }
        }
    }

    for(i=0; i<num_erb; i++)
    {
        out[i] = curr_erb[i] - t_prev_erb[i] ;
    }

    /* First Band Amplitude Search */
    mmse=(float) HUGE_VAL;
    mmseindex=-1 ;
    for(j=0; j<ERB_CBSIZE1; j++)
    {
        tmp = 0.0;
        for(i=1; i<11; i++)
        {
            if(cslot[i]!=0)
            {
                float tmp1;
                if(AmpCB1[j][i-1]<-t_prev_erb[i])
                {
                    tmp1 = PowSpect[i] * SQR (curr_erb[i]) ;
                }
                else
                {
                    tmp1 =(float) (PowSpect[i] * SQR (out[i] - AmpCB1[j][i-1])) ;
                }
                if(AmpCB1[j][i-1] < out[i])
                {
                    tmp1 *= 0.9f ;
                }
                tmp += tmp1 ;
            }
        }

        if(tmp < mmse)
        {
            mmse = tmp ;
            mmseindex = j ;
        }
    }

    if (!(mmseindex<ERB_CBSIZE1 && mmseindex>=0))
    {
        mmseindex=0;
    }

    index[0] = mmseindex ;

    /* Second Band Amplitude Search */
    mmse = (float) HUGE_VAL;
    mmseindex=-1 ;
    for(j=0; j<ERB_CBSIZE2; j++)
    {
        tmp = 0.0;
        for(i=11; i<num_erb; i++)
        {
            if (num_erb==NUM_ERB_NB)
            {
                if(cslot[i]!=0)
                {
                    if(AmpCB2_NB[j][i-11]<-t_prev_erb[i])
                    {
                        tmp1 = PowSpect[i] * SQR (curr_erb[i]) ;
                    }
                    else
                    {
                        tmp1 = (float)(PowSpect[i] * SQR (out[i] - AmpCB2_NB[j][i-11])) ;
                    }

                    if(AmpCB2_NB[j][i-11] < out[i])
                    {
                        tmp1 *= 0.9f ;
                    }
                    tmp += tmp1 ;

                }
            }
            else if (num_erb==NUM_ERB_WB)
            {
                if(cslot[i]!=0)
                {
                    if(AmpCB2_WB[j][i-11]<-t_prev_erb[i])
                    {
                        tmp1 = PowSpect[i] * SQR (curr_erb[i]) ;
                    }
                    else
                    {
                        tmp1 = (float)(PowSpect[i] * SQR (out[i] - AmpCB2_WB[j][i-11])) ;
                    }

                    if(AmpCB2_WB[j][i-11] < out[i])
                    {
                        tmp1 *= 0.9f ;
                    }
                    tmp += tmp1 ;

                }
            }
        }

        if(tmp < mmse)
        {
            mmse = tmp ;
            mmseindex = j ;
        }
    }

    if (!(mmseindex<ERB_CBSIZE2 && mmseindex>=0))
    {
        mmseindex=0;
    }

    index[1] = mmseindex;

    return;
}

/*-------------------------------------------------------------------*
 * erb_add()
 *
 * Add current and past ERB
 *-------------------------------------------------------------------*/
static
void erb_add(
    float *curr_erb,       /* i/o:  current ERB    */
    int   l,               /* i  :  current lag    */
    const float *prev_erb,       /* i  :  previous ERB   */
    int   pl,              /* i  :  previous lag   */
    const int   *index,          /* i  :  ERB index      */
    short num_erb          /* i  :  number of ERBs */
)
{
    int i, pslot[NUM_ERB_WB], cslot[NUM_ERB_WB] ;
    float tmp, t_prev_erb[NUM_ERB_WB] ;
    const float (*AmpCB1)[10] = 0;

    if (num_erb==NUM_ERB_NB)
    {
        AmpCB1=AmpCB1_NB;
    }
    else if (num_erb==NUM_ERB_WB)
    {
        AmpCB1=AmpCB1_WB;
    }
    erb_slot (l, cslot, t_prev_erb, num_erb);
    erb_slot (pl, pslot, t_prev_erb, num_erb);

    for(i=0; i<num_erb; i++)
    {
        t_prev_erb[i] = prev_erb[i] ;
    }

    if(pl>l)
    {
        tmp = t_prev_erb[0] ;
        for(i=0; i<num_erb; i++)
        {
            if (!(pslot[i]>=0))
            {
                pslot[i] = 0;
            }

            if(pslot[i]!=0)
            {
                tmp = t_prev_erb[i] ;
            }
            else
            {
                t_prev_erb[i] = tmp ;
            }
        }
    }
    else if (l>pl)
    {
        tmp = t_prev_erb[num_erb-1] ;
        for(i=num_erb-1; i>=0; i--)
        {
            if(pslot[i]!=0)
            {
                tmp = t_prev_erb[i] ;
            }
            else
            {
                t_prev_erb[i] = tmp ;
            }
        }
    }

    for(i=1; i<11; i++)
    {
        if(cslot[i]!=0)
        {
            curr_erb[i] = (float) (AmpCB1[index[0]][i-1] + t_prev_erb[i]) ;
            curr_erb[i] = max(0.0f, curr_erb[i]);
        }
        else
        {
            curr_erb[i] = 0.0;
        }
    }
    for(i=11; i<(num_erb-2); i++)
    {
        if(cslot[i]!=0)
        {
            if (num_erb==NUM_ERB_NB)
            {
                curr_erb[i] = (float)(AmpCB2_NB[index[1]][i-11] + t_prev_erb[i]) ;
                curr_erb[i] = max(0.0f, curr_erb[i]);
            }
            else if (num_erb==NUM_ERB_WB)
            {
                curr_erb[i] = (float)(AmpCB2_WB[index[1]][i-11] + t_prev_erb[i]) ;
                curr_erb[i] = max(0.0f, curr_erb[i]);
            }
        }
        else
        {
            curr_erb[i] = 0.0f;
        }
    }
}

/*-------------------------------------------------------------------*
 * DTFS_quant_cw()
 *
 * DTFS quantization
 *-------------------------------------------------------------------*/

short DTFS_quant_cw(
    DTFS_STRUCTURE *X,           /* i/o: DTFS unquant inp, quant out  */
    int   pl,           /* i  : Previous lag                 */
    const float *curr_lpc,    /* i  : LPC                          */
    int   *POWER_IDX,   /* o  : Power index                  */
    int   *AMP_IDX,     /* o  : Amplitude index              */
    float *lastLgainE,  /* i/o: last frame lowband gain      */
    float *lastHgainE,  /* i/o: last frame highband gain     */
    float *lasterbE     /* i/o: last frame ERB vector        */
)
{
    short num_erb = 0;
    const float (*PowerCB)[2]=NULL;
    float G_CURR_ERB[NUM_ERB_WB];
    float G_A_POWER[2];
    float tmp, w[2], target1, target2, error, minerror ;
    float mfreq[NUM_ERB_WB], diff_erb[NUM_ERB_WB], curr_erb[NUM_ERB_WB] ;
    int j, slot[NUM_ERB_WB],bincount;
    short returnFlag = 1;
    float amperror;

    if (X->upper_cut_off_freq==4000.0)
    {
        num_erb=NUM_ERB_NB;
        PowerCB=PowerCB_NB;
    }
    else if (X->upper_cut_off_freq==6400.0)
    {
        num_erb=NUM_ERB_WB;
        PowerCB=PowerCB_WB;
    }

    /* Getting the Speech Domain Energy LOG Ratio */
    w[0] = (float) max(1E-10, log10(DTFS_getSpEngyFromResAmp (*X,0.0,1104.5,curr_lpc)));
    w[1] = (float) max(1E-10, log10(DTFS_getSpEngyFromResAmp (*X,1104.5,X->upper_cut_off_freq,curr_lpc)));
    tmp = w[0]+w[1] ;
    w[0] /= tmp;
    w[1] /= tmp;

    /* Power Quantization */
    G_A_POWER[0] = (float) log10(X->lag*DTFS_setEngyHarm (92.0,1104.5,0.0,1104.5,1.0,X));
    G_A_POWER[1] = (float) log10(X->lag*DTFS_setEngyHarm (1104.5,X->upper_cut_off_freq_of_interest,1104.5,X->upper_cut_off_freq,1.0,X));
    target1 = G_A_POWER[0]-*lastLgainE ;
    target2 = G_A_POWER[1]-*lastHgainE ;
    minerror=(float) HUGE_VAL ;
    *POWER_IDX=-1 ;
    for(j=0; j<P_CBSIZE; j++)
    {
        error = (float)(w[0]*fabs(target1 - PowerCB[j][0]) + w[1]*fabs(target2 - PowerCB[j][1])) ;
        if((target1>=PowerCB[j][0]) &&  (target2>=PowerCB[j][1]))
        {
            error*=0.8f ;
        }
        if(error < minerror)
        {
            minerror = error ;
            *POWER_IDX = j ;
        }
    }
    DTFS_to_erb (*X,curr_erb) ;

    for(j=0; j<num_erb; j++)
    {
        G_CURR_ERB[j] = curr_erb[j];
    }
    erb_slot (X->lag, slot, mfreq,num_erb) ;
    /* Amplitude Quantization */
    erb_diff (lasterbE, pl, curr_erb, X->lag, curr_lpc, diff_erb, AMP_IDX,num_erb) ;

    /* Amplitude Dequantization */
    erb_add (curr_erb, X->lag, lasterbE, pl, AMP_IDX,num_erb) ;
    curr_erb[0] = curr_erb[1] * 0.3f ;
    curr_erb[num_erb-2] = curr_erb[num_erb-3] * 0.3f;
    curr_erb[num_erb-1] = 0;
    /* Determine if the amplitude quantization is good enough */
    amperror=0.0;
    bincount=0;
    for(j=1; j<10; j++)
    {
        if(slot[j]!=0)
        {
            amperror += (float)(fabs(G_CURR_ERB[j]-curr_erb[j]));
            bincount++;
        }
    }
    amperror /= bincount ;

    if((amperror>0.47) &&  (target1>-0.4))
    {
        returnFlag = 0 ; /* Bumping up */
    }

    DTFS_erb_inv (curr_erb, slot, mfreq, X,num_erb) ;

    /* Back up the lasterbE memory after power normalization */
    DTFS_setEngyHarm (92.0,1104.5,0.0,1104.5,1.0,X) ;
    DTFS_setEngyHarm (1104.5,X->upper_cut_off_freq_of_interest,1104.5,X->upper_cut_off_freq,1.0,X) ;

    DTFS_to_erb(*X,lasterbE);

    /* Power Dequantization */
    *lastLgainE += (float) PowerCB[*POWER_IDX][0] ;
    *lastHgainE += (float) PowerCB[*POWER_IDX][1] ;
    target1 = (float) pow(10.0, (double)(*lastLgainE))/X->lag ;

    if (!(target1 >= 0.0))
    {
        target1 = 0;
    }

    DTFS_setEngyHarm (92.0f,1104.5f,0.0f,1104.5f,target1,X);
    target2 = (float) pow(10.0, (double)(*lastHgainE))/X->lag ;

    if (!(target2 >= 0.0))
    {
        target2 = 0;
    }

    DTFS_setEngyHarm (1104.5,X->upper_cut_off_freq_of_interest,1104.5,X->upper_cut_off_freq,target2,X);

    return returnFlag;   /* amp quant performance pass/fail   */
}

/*-------------------------------------------------------------------*
 * DTFS_dequant_cw()
 *
 * DTFS dequantization
 *-------------------------------------------------------------------*/

void DTFS_dequant_cw(
    int   pl,              /* i  : Previous lag                 */
    int   POWER_IDX,       /* i  : POWER index                  */
    const int   *AMP_IDX,        /* i  : Amp Shape index              */
    float *lastLgainD,     /* i/o: low band last gain           */
    float *lastHgainD,     /* i/o: high band last gain          */
    float *lasterbD,       /* i/o: last frame ERB vector        */
    DTFS_STRUCTURE *X,              /* o  : DTFS structure dequantized   */
    short num_erb          /* i  : Number of ERB bands          */
)
{
    float tmp, mfreq[NUM_ERB_WB], curr_erb[NUM_ERB_WB];
    int slot[NUM_ERB_WB] ;
    const float (*PowerCB)[2]=NULL;

    if (num_erb==NUM_ERB_NB)
    {
        PowerCB=PowerCB_NB;
    }
    else if (num_erb==NUM_ERB_WB)
    {
        PowerCB=PowerCB_WB;
    }

    /* Amplitude Dequantization */
    erb_add (curr_erb, X->lag, lasterbD, pl, AMP_IDX,num_erb) ;
    curr_erb[0] = curr_erb[1] * 0.3f ;
    curr_erb[num_erb-2] = curr_erb[num_erb-3] * 0.3f;
    curr_erb[num_erb-1] = 0;
    erb_slot(X->lag, slot, mfreq,num_erb) ;
    DTFS_erb_inv(curr_erb, slot, mfreq, X,num_erb) ;

    /* Back up the lasterbD memory after power normalization */
    DTFS_setEngyHarm(92.0,1104.5,0.0,1104.5,1.0,X) ;
    DTFS_setEngyHarm(1104.5,X->upper_cut_off_freq_of_interest,1104.5,X->upper_cut_off_freq,1.0,X) ;
    DTFS_to_erb(*X,lasterbD) ;

    /* Power Dequantization */
    *lastLgainD += (float) PowerCB[POWER_IDX][0] ;
    *lastHgainD += (float) PowerCB[POWER_IDX][1] ;
    tmp = (float) pow(10.0, (double)(*lastLgainD))/X->lag ;

    if (!(tmp >= 0.0))
    {
        tmp = 0.0;
    }

    DTFS_setEngyHarm (92.0,1104.5,0.0,1104.5,tmp,X) ;
    tmp = (float) pow(10.0, (double)(*lastHgainD))/X->lag ;

    if (!(tmp >= 0.0))
    {
        tmp = 0.0;
    }

    DTFS_setEngyHarm (1104.5,X->upper_cut_off_freq_of_interest,1104.5,X->upper_cut_off_freq,tmp,X) ;
}

/*-------------------------------------------------------------------*
 * WI_syn()
 *
 * Synthesis using waveform interpolation
 *-------------------------------------------------------------------*/

void WIsyn(
    DTFS_STRUCTURE PREVCW,            /* i  : Prev frame DTFS                                */
    DTFS_STRUCTURE *CURR_CW_DTFS,     /* i/o: Curr frame DTFS                                */
    const float *curr_lpc,         /* i  : LPC                                            */
    float *ph_offset,        /* i/o: Phase offset to line up at end of frame        */
    float *out,              /* o  : Waveform Interpolated time domain signal       */
    int   N,                 /* i  : Number of output samples to generate           */
    int FR_flag              /* i  : called for post-smoothing in FR                */
)
{
    DTFS_STRUCTURE *CURRCW_DTFS;
    unsigned short I=1, flag=0;
    float alignment, tmp, *phase;

    phase = (float *) calloc(N,sizeof(float));

    CURRCW_DTFS = DTFS_new();

    DTFS_copy (CURRCW_DTFS,*CURR_CW_DTFS);

    /* Calculating the expected alignment shift */
    alignment = (float) (*ph_offset /PI2 * PREVCW.lag) ;
    if(flag==1)
    {
        alignment *= I ;
    }
    /* Calculating the expected alignment shift */
    tmp = (float) fmod((N%((PREVCW.lag+CURRCW_DTFS->lag)>>1) + alignment), CURRCW_DTFS->lag) ;

    /* Compute the alignment shift */
    if (FR_flag==0)
    {
        alignment = DTFS_alignment_weight (PREVCW, *CURRCW_DTFS, tmp, curr_lpc, curr_lpc) ;
    }
    else /* FR case */
    {
        alignment = DTFS_alignment_full(PREVCW, *CURRCW_DTFS, CURRCW_DTFS->lag*2);
    }

    tmp = (float) (PI2*alignment/CURRCW_DTFS->lag) ;
    DTFS_phaseShift (CURRCW_DTFS,tmp) ;
    DTFS_phaseShift (CURR_CW_DTFS,(float) (PI2*alignment/CURR_CW_DTFS->lag)) ;

    /* Compute the cubic phase track and transform to 1-D signal */
    cubicPhase (*ph_offset, tmp, (float) PREVCW.lag , (float) CURRCW_DTFS->lag, N, phase) ;

    if (FR_flag==0)
    {
        DTFS_transform (PREVCW, *CURRCW_DTFS, phase, out, N,0) ;
    }
    else
    {
        DTFS_transform (PREVCW, *CURRCW_DTFS, phase, out, N,1) ;
    }

    /* Adjust the phase offset and wrap it between 0 and 2pi */
    if( flag == 2 )
    {
        tmp *= I ;
    }
    *ph_offset = (float) fmod ((double)(tmp), PI2) ;

    free(phase) ;
    free(CURRCW_DTFS);
}


/*-------------------------------------------------------------------*
 * ppp_extract_pitch_period()
 *
 *  Pitch period extraction
 *-------------------------------------------------------------------*/

short ppp_extract_pitch_period(
    const float *in,            /* i : input residual     */
    float *out,           /* o : output residual    */
    int   l,              /* i : lag                */
    short *out_of_bound   /* o : out of bound flag  */
)
{
    int i,j,k;
    int spike=0,range;
    float max1=0.0;
    const float *ptr=in+L_FRAME-l;
    float en1=0.0, en2=0.0, tmp;
    short spike_near_edge=0;
    float pos_max,neg_max;
    int spike_pos=0,spike_neg=0;
    float x;

    pos_max = (float) -HUGE_VAL;
    neg_max = 0.0;
    *out_of_bound = 0;

    for(i=0; i<l; i++)
    {
        if((x=(float) fabs(ptr[i]))>max1)
        {
            max1=x;
            spike=i;
        }
        en1 += ptr[i]*ptr[i] ;
    }

    if (ptr[spike]>0)
    {
        spike_pos=spike;

        /* search for neg spike around the pos spike */
        for(j=spike-10; j<spike+10; j++)
        {
            k=(j+l)%l;
            if (ptr[k]<neg_max)
            {
                neg_max=ptr[k];
                spike_neg=k;
            }
        }
    }
    else if (ptr[spike]<0)
    {
        spike_neg=spike;

        /* search for pos spike around the neg spike */
        for(j=spike-10; j<spike+10; j++)
        {
            k=(j+l)%l;
            if (ptr[k]>pos_max)
            {
                pos_max=ptr[k];
                spike_pos=k;
            }
        }
    }
    if (((l-1-max(spike_pos,spike_neg))<=2) ||(min(spike_pos,spike_neg)<=2))
    {
        *out_of_bound=1;
        return spike_near_edge;
    }
    range=(int)anint (max(CUTFREE_REL_RANGE*l,CUTFREE_ABS_RANGE));
    if((spike-range<0) ||  (spike+range>=l))
    {
        /* need to grab from one lag before
        ensure that there is no array bound read */
        if(L_FRAME -l -l < 0)
        {
            *out_of_bound=1;
            return spike_near_edge;
        }
        spike_near_edge=1;
    }
    if(spike-range<0)
    {
        for(i=0; i<l+(spike-range); i++)
        {
            out[i]=ptr[i];
        }

        /* Grab Remaining From One Lag Before */
        ptr-=l;
        for(; i<l; i++)
        {
            out[i]=ptr[i];
        }
    }
    else if (spike+range>=l)
    {
        for(i=0; i<spike-range; i++)
        {
            out[i]=ptr[i];
        }
        /* Grab Remaining From One Lag Before */
        if (ptr-l+i>=in)
        {
            for(ptr-=l; i<l; i++)
            {
                out[i]=ptr[i];
            }
        }
        else
        {
            for(; i<l; i++)
            {
                out[i]=ptr[i];
            }
        }
    }
    else
    {
        for(i=0; i<l; i++)
        {
            out[i]=ptr[i];
        }
    }

    /* Energy adjustment added to eliminate artifacts at the end of voicing */
    for(i=0; i<l; i++)
    {
        en2 += out[i]*out[i] ;
    }

    if(en1<en2)
    {
        tmp=(float) sqrt(en1/en2);
        for(i=0; i<l; i++)
        {
            out[i] *= tmp;
        }
    }

    return spike_near_edge;
}

/*-------------------------------------------------------------------*
 * DTFS_peaktoaverage()
 *
 * Estimate peak to average ratio in the DTFS
 *-------------------------------------------------------------------*/

void DTFS_peaktoaverage(
    DTFS_STRUCTURE X,          /* i : DTFS                  */
    float *pos,       /* o : positive peak to ave */
    float *neg        /* o : negative peak to ave */
)
{
    float tmp, time[PIT_MAX], sum=0.0, maxPosEn=0.0, maxNegEn=0.0;
    short i;

    DTFS_fs_inv (&X,time,X.lag,0.0);

    for(i=0; i<X.lag; i++)
    {
        tmp = SQR(time[i]);
        if(time[i] >= 0)
        {
            if(tmp>maxPosEn)
            {
                maxPosEn = tmp ;
            }
        }
        else
        {
            if(tmp>maxNegEn)
            {
                maxNegEn = tmp ;
            }
        }
        sum += tmp ;
    }


    if (sum==0.0)
    {
        *pos = 0.0f;
        *neg = 0.0f;
    }
    else
    {

        if (maxPosEn==0.0)
        {
            *pos=0.0;
        }
        else
        {
            *pos = (float) sqrt(maxPosEn*X.lag/sum);
        }

        if (maxPosEn==0.0)
        {
            *neg=0.0;
        }
        else
        {
            *neg = (float) sqrt(maxNegEn*X.lag/sum);
        }
    }
}

