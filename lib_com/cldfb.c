/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <math.h>
#include "stat_dec.h"
#include "prot.h"
#include "rom_com.h"
#include "string.h"
#include <assert.h>

#if defined __ICL
#define restrict __restrict
#else
#define restrict
#endif

#ifdef _MSC_VER
#pragma warning(disable : 4305) /* disable truncation from double to float warning (VC++)*/
#endif


/*-------------------------------------------------------------------*
 * local prototypes
 *--------------------------------------------------------------------*/

static void cldfb_init_proto_and_twiddles(HANDLE_CLDFB_FILTER_BANK hs);

static float GetEnergyCldfb( float *energyValuesSum, float *energyLookahead, float **realValues, float **imagValues,
                             int numberBands, int numberCols, HANDLE_TEC_ENC hTecEnc );

/*-------------------------------------------------------------------*
 * cplxMult()
 *
 * Conduct complex multiplication
 *--------------------------------------------------------------------*/
static void cplxMult(
    float *yr,   /* o  : real output */
    float *yi,   /* o  : imag output */
    float xr,    /* i  : real input 1*/
    float xi,    /* i  : imag input 1*/
    float cr,    /* i  : real input 1*/
    float ci     /* i  : imag input 1*/
)
{
    *yr = xr*cr - xi*ci;
    *yi = xr*ci + xi*cr;

    return;
}


/*-------------------------------------------------------------------*
 * cldfbAnalysis()
 *
 * Conduct multple overlap cmplex low delay MDCT
 *--------------------------------------------------------------------*/
void cldfbAnalysis(
    const float                 *timeIn,              /* i  : time buffer */
    float                       **realBuffer,         /* o  : real value buffer */
    float                       **imagBuffer,         /* o  : imag value buffer */
    int                          samplesToProcess,    /* i  : samples to process */
    HANDLE_CLDFB_FILTER_BANK     h_cldfb              /* i  : filterbank state */
)
{
    short i, k;
    short L2, M1, M2, M4;
    short no_col = h_cldfb->no_col;

    float r1, r2, rr12, ir12;
    float i1, i2, ri12, ii12;
    float rBuffer[2*CLDFB_NO_CHANNELS_MAX];
    float iBuffer[2*CLDFB_NO_CHANNELS_MAX];
    const float *rot_vctr_re;
    const float *rot_vctr_im;
    const float *ptr_pf;
    float *timeBuffer = h_cldfb->cldfb_state;
    int offset        = h_cldfb->p_filter_length - h_cldfb->no_channels;
    int frameSize     = h_cldfb->no_channels * h_cldfb->no_col;

    /* prepare input buffer */
    mvr2r( timeBuffer + frameSize, timeBuffer, offset );

    if( samplesToProcess > -1 )
    {
        mvr2r( timeIn, timeBuffer + offset, samplesToProcess );
        set_f( timeBuffer+offset+samplesToProcess, 0.0f, (frameSize-samplesToProcess) );
    }
    else
    {
        mvr2r( timeIn, timeBuffer + offset, frameSize );
    }

    /* only process needed cols */
    if( samplesToProcess > -1 )
    {
        no_col = min(no_col, (samplesToProcess + h_cldfb->no_channels - 1) / h_cldfb->no_channels);
    }

    M1  = h_cldfb->no_channels;
    M2 = M1 >> 1;
    M4 = M1 >> 2;
    L2 = M1 << 1;

    if( M2 & 1 )
    {
        M4 += 1;
    }

    rot_vctr_re = h_cldfb->rot_vec_ana_re;
    rot_vctr_im = h_cldfb->rot_vec_ana_im;

    ptr_pf = h_cldfb->p_filter;

    for( i = 0; i < no_col; i++ )
    {
        for (k=0; k < M4; k++ )
        {
            /* prototype filter */
            r1 = 0  - ptr_pf[(L2-M2-1-(2*k) + 0 * L2)] * timeBuffer[L2-M2-1-(2*k) + 0 * L2];
            r1 = r1 - ptr_pf[(L2-M2-1-(2*k) + 1 * L2)] * timeBuffer[L2-M2-1-(2*k) + 1 * L2];
            r1 = r1 - ptr_pf[(L2-M2-1-(2*k) + 2 * L2)] * timeBuffer[L2-M2-1-(2*k) + 2 * L2];
            r1 = r1 - ptr_pf[(L2-M2-1-(2*k) + 3 * L2)] * timeBuffer[L2-M2-1-(2*k) + 3 * L2];
            r1 = r1 - ptr_pf[(L2-M2-1-(2*k) + 4 * L2)] * timeBuffer[L2-M2-1-(2*k) + 4 * L2];

            r2 = 0  - ptr_pf[(L2-M2+(2*k) + 0 * L2)] * timeBuffer[L2-M2+(2*k) + 0 * L2];
            r2 = r2 - ptr_pf[(L2-M2+(2*k) + 1 * L2)] * timeBuffer[L2-M2+(2*k) + 1 * L2];
            r2 = r2 - ptr_pf[(L2-M2+(2*k) + 2 * L2)] * timeBuffer[L2-M2+(2*k) + 2 * L2];
            r2 = r2 - ptr_pf[(L2-M2+(2*k) + 3 * L2)] * timeBuffer[L2-M2+(2*k) + 3 * L2];
            r2 = r2 - ptr_pf[(L2-M2+(2*k) + 4 * L2)] * timeBuffer[L2-M2+(2*k) + 4 * L2];

            i1 = 0  - ptr_pf[(L2-3*M2+(2*k) + 0 * L2)] * timeBuffer[L2-3*M2+(2*k) + 0 * L2];
            i1 = i1 - ptr_pf[(L2-3*M2+(2*k) + 1 * L2)] * timeBuffer[L2-3*M2+(2*k) + 1 * L2];
            i1 = i1 - ptr_pf[(L2-3*M2+(2*k) + 2 * L2)] * timeBuffer[L2-3*M2+(2*k) + 2 * L2];
            i1 = i1 - ptr_pf[(L2-3*M2+(2*k) + 3 * L2)] * timeBuffer[L2-3*M2+(2*k) + 3 * L2];
            i1 = i1 - ptr_pf[(L2-3*M2+(2*k) + 4 * L2)] * timeBuffer[L2-3*M2+(2*k) + 4 * L2];

            i2 = 0  - ptr_pf[(L2-3*M2-1-(2*k) + 0 * L2)] * timeBuffer[L2-3*M2-1-(2*k) + 0 * L2];
            i2 = i2 - ptr_pf[(L2-3*M2-1-(2*k) + 1 * L2)] * timeBuffer[L2-3*M2-1-(2*k) + 1 * L2];
            i2 = i2 - ptr_pf[(L2-3*M2-1-(2*k) + 2 * L2)] * timeBuffer[L2-3*M2-1-(2*k) + 2 * L2];
            i2 = i2 - ptr_pf[(L2-3*M2-1-(2*k) + 3 * L2)] * timeBuffer[L2-3*M2-1-(2*k) + 3 * L2];
            i2 = i2 - ptr_pf[(L2-3*M2-1-(2*k) + 4 * L2)] * timeBuffer[L2-3*M2-1-(2*k) + 4 * L2];

            /* folding + pre modulation of DST IV */
            rr12 =  r1 - r2;
            ri12 = -i1 - i2;
            cplxMult(&rBuffer[2*k],&rBuffer[2*k+1],rr12,ri12,rot_vctr_re[k],rot_vctr_im[k]);

            /* folding + pre modulation of DCT IV */
            ir12 = r1 + r2;
            ii12 = i1 - i2;
            cplxMult(&iBuffer[2*k],&iBuffer[2*k+1],ir12,ii12,rot_vctr_re[k],rot_vctr_im[k]);
        }

        for (k=M4; k < M2; k++)
        {
            /* prototype filter */
            r1 = 0  - ptr_pf[(L2-M2-1-(2*k) + 0 * L2)] * timeBuffer[L2-M2-1-(2*k) + 0 * L2];
            r1 = r1 - ptr_pf[(L2-M2-1-(2*k) + 1 * L2)] * timeBuffer[L2-M2-1-(2*k) + 1 * L2];
            r1 = r1 - ptr_pf[(L2-M2-1-(2*k) + 2 * L2)] * timeBuffer[L2-M2-1-(2*k) + 2 * L2];
            r1 = r1 - ptr_pf[(L2-M2-1-(2*k) + 3 * L2)] * timeBuffer[L2-M2-1-(2*k) + 3 * L2];
            r1 = r1 - ptr_pf[(L2-M2-1-(2*k) + 4 * L2)] * timeBuffer[L2-M2-1-(2*k) + 4 * L2];

            r2 = 0  - ptr_pf[(L2-5*M2+(2*k) + 0 * L2)] * timeBuffer[L2-5*M2+(2*k) + 0 * L2];
            r2 = r2 - ptr_pf[(L2-5*M2+(2*k) + 1 * L2)] * timeBuffer[L2-5*M2+(2*k) + 1 * L2];
            r2 = r2 - ptr_pf[(L2-5*M2+(2*k) + 2 * L2)] * timeBuffer[L2-5*M2+(2*k) + 2 * L2];
            r2 = r2 - ptr_pf[(L2-5*M2+(2*k) + 3 * L2)] * timeBuffer[L2-5*M2+(2*k) + 3 * L2];
            r2 = r2 - ptr_pf[(L2-5*M2+(2*k) + 4 * L2)] * timeBuffer[L2-5*M2+(2*k) + 4 * L2];

            i1 = 0  - ptr_pf[(L2+M2-1-(2*k) + 0 * L2)] * timeBuffer[L2+M2-1-(2*k) + 0 * L2];
            i1 = i1 - ptr_pf[(L2+M2-1-(2*k) + 1 * L2)] * timeBuffer[L2+M2-1-(2*k) + 1 * L2];
            i1 = i1 - ptr_pf[(L2+M2-1-(2*k) + 2 * L2)] * timeBuffer[L2+M2-1-(2*k) + 2 * L2];
            i1 = i1 - ptr_pf[(L2+M2-1-(2*k) + 3 * L2)] * timeBuffer[L2+M2-1-(2*k) + 3 * L2];
            i1 = i1 - ptr_pf[(L2+M2-1-(2*k) + 4 * L2)] * timeBuffer[L2+M2-1-(2*k) + 4 * L2];

            i2 = 0  - ptr_pf[(L2-3*M2+(2*k) + 0 * L2)] * timeBuffer[L2-3*M2+(2*k) + 0 * L2];
            i2 = i2 - ptr_pf[(L2-3*M2+(2*k) + 1 * L2)] * timeBuffer[L2-3*M2+(2*k) + 1 * L2];
            i2 = i2 - ptr_pf[(L2-3*M2+(2*k) + 2 * L2)] * timeBuffer[L2-3*M2+(2*k) + 2 * L2];
            i2 = i2 - ptr_pf[(L2-3*M2+(2*k) + 3 * L2)] * timeBuffer[L2-3*M2+(2*k) + 3 * L2];
            i2 = i2 - ptr_pf[(L2-3*M2+(2*k) + 4 * L2)] * timeBuffer[L2-3*M2+(2*k) + 4 * L2];

            /* folding + pre modulation of DST IV */
            rr12 = r1 + r2;
            ri12 = i1 - i2;
            cplxMult(&rBuffer[2*k],&rBuffer[2*k+1],rr12,ri12,rot_vctr_re[k],rot_vctr_im[k]);

            /* folding + pre modulation of DCT IV */
            ir12 = r1 - r2;
            ii12 = i1 + i2;
            cplxMult(&iBuffer[2*k],&iBuffer[2*k+1],ir12,ii12,rot_vctr_re[k],rot_vctr_im[k]);
        }

        /* FFT of DST IV */
        fft_cldfb(rBuffer, M2);

        /* post modulation of DST IV */
        for (k=0; k < M2; k++)
        {
            cplxMult(&realBuffer[i][M1-1-(2*k)],&realBuffer[i][2*k],rBuffer[2*k],rBuffer[2*k+1],rot_vctr_re[k],rot_vctr_im[k]);
        }

        /* FFT of DCT IV */
        fft_cldfb(iBuffer, M2);

        /* post modulation of DCT IV */
        for (k=0; k < M2; k++)
        {
            /* do it inplace */
            cplxMult(&imagBuffer[i][2*k],&imagBuffer[i][M1-1-(2*k)],iBuffer[2*k],iBuffer[2*k+1],rot_vctr_re[k],rot_vctr_im[k]);
        }

        timeBuffer += L2 * 5;
        timeBuffer += h_cldfb->no_channels - h_cldfb->p_filter_length;
    }

    return;
}


/*-------------------------------------------------------------------*
 * cldfbSynthesis()
 *
 * Conduct inverse multple overlap cmplex low delay MDCT
 *--------------------------------------------------------------------*/
void cldfbSynthesis(
    float                    **realBuffer,         /* i  : real values */
    float                    **imagBuffer,         /* i  : imag values */
    float                     *timeOut,            /* o  : output time domain samples */
    int                        samplesToProcess,   /* i  : number of processed samples */
    HANDLE_CLDFB_FILTER_BANK   h_cldfb             /* i  : filter bank state */
)
{
    int i;
    int k;
    int L2;
    int M1;
    int M2;
    int M41;
    int M42;
    int Mz;

    float rBuffer[2*CLDFB_NO_CHANNELS_MAX];
    float iBuffer[2*CLDFB_NO_CHANNELS_MAX];
    const float *rot_vctr_re;
    const float *rot_vctr_im;
    float rr12, ir12;
    float ri12, ii12;

    float *synthesisBuffer;
    float new_samples[2*CLDFB_NO_CHANNELS_MAX];

    float *ptr_time_out;
    const float *p_filter;

    float accu0, accu1, accu2, accu3, accu4;
    int no_col = h_cldfb->no_col;

    M1  = h_cldfb->no_channels;
    L2  = M1 << 1;
    M2  = M1 >> 1;
    M41 = M2>>1;
    M42 = M2-M41;
    Mz  = M1 - h_cldfb->bandsToZero;

    /* only process needed cols */
    if(samplesToProcess > -1)
    {
        no_col = min(no_col, (samplesToProcess + h_cldfb->no_channels - 1) / h_cldfb->no_channels);
    }

    rot_vctr_re = h_cldfb->rot_vec_syn_re;
    rot_vctr_im = h_cldfb->rot_vec_syn_im;

    synthesisBuffer = h_cldfb->cldfb_state;
    p_filter = h_cldfb->p_filter;
    ptr_time_out = timeOut;

    mvr2r( synthesisBuffer, synthesisBuffer + M1 * h_cldfb->no_col, h_cldfb->p_filter_length );

    synthesisBuffer += M1 * h_cldfb->no_col;

    for (k=0; k < no_col; k++)
    {
        for (i=Mz; i < M1; i++)
        {
            realBuffer[k][i] = 0.0f;
            imagBuffer[k][i] = 0.0f;
        }

        for (i=0; i < M2; i++)
        {
            /* pre modulation of DST IV */
            cplxMult(&rBuffer[2*i], &rBuffer[2*i+1], realBuffer[k][2*i], realBuffer[k][M1-1-2*i], rot_vctr_re[i], rot_vctr_im[i]);

            /* pre modulation of DCT IV */
            cplxMult(&iBuffer[2*i], &iBuffer[2*i+1],-imagBuffer[k][2*i], imagBuffer[k][M1-1-2*i], rot_vctr_re[i], rot_vctr_im[i]);
        }

        /* FFT of DST IV */
        fft_cldfb(rBuffer, M2);

        /* FFT of DCT IV */
        fft_cldfb(iBuffer, M2);

        /* folding */
        for (i=0; i<M41; i++)
        {
            /* post modulation of DST IV */
            rr12 = rBuffer[M1-2-2*i]*rot_vctr_re[M2-1-i] - rBuffer[M1-1-2*i]*rot_vctr_im[M2-1-i];
            ri12 = rBuffer[M1-2-2*i]*rot_vctr_im[M2-1-i] + rBuffer[M1-1-2*i]*rot_vctr_re[M2-1-i];

            /* post modulation of DCT IV */
            ir12 = iBuffer[M1-2-2*i]*rot_vctr_re[M2-1-i] - iBuffer[M1-1-2*i]*rot_vctr_im[M2-1-i];
            ii12 = iBuffer[M1-2-2*i]*rot_vctr_im[M2-1-i] + iBuffer[M1-1-2*i]*rot_vctr_re[M2-1-i];

            new_samples[M1+M2+1+2*i] = -rr12 - ii12;
            new_samples[M2-2-2*i]   = -ri12 - ir12;

            new_samples[M1+M2-2-2*i] =  rr12 - ii12;
            new_samples[M2+1+2*i]   =  ir12 - ri12;
        }

        for (i=0; i<M42; i++)
        {
            /* post modulation of DST IV */
            rr12 = rBuffer[2*i]*rot_vctr_re[i] - rBuffer[2*i+1]*rot_vctr_im[i];
            ri12 = rBuffer[2*i]*rot_vctr_im[i] + rBuffer[2*i+1]*rot_vctr_re[i];

            /* post modulation of DCT IV */
            ir12 = iBuffer[2*i]*rot_vctr_re[i] - iBuffer[2*i+1]*rot_vctr_im[i];
            ii12 = iBuffer[2*i]*rot_vctr_im[i] + iBuffer[2*i+1]*rot_vctr_re[i];

            new_samples[M1+M2+2*i]   = ri12 +  ir12;
            new_samples[M2-1-2*i]   = rr12 +  ii12;

            new_samples[M1+M2-1-2*i] = ir12 - ri12;
            new_samples[M2+2*i]     = rr12 - ii12;
        }

        /* synthesis prototype filter */
        for (i=0; i < L2; i++)
        {
            accu0 = synthesisBuffer[0 * L2 + i] + p_filter[(0 * L2 + i)] * new_samples[L2 - 1 - i];
            accu1 = synthesisBuffer[1 * L2 + i] + p_filter[(1 * L2 + i)] * new_samples[L2 - 1 - i];
            accu2 = synthesisBuffer[2 * L2 + i] + p_filter[(2 * L2 + i)] * new_samples[L2 - 1 - i];
            accu3 = synthesisBuffer[3 * L2 + i] + p_filter[(3 * L2 + i)] * new_samples[L2 - 1 - i];
            accu4 = synthesisBuffer[4 * L2 + i] + p_filter[(4 * L2 + i)] * new_samples[L2 - 1 - i];

            synthesisBuffer[0 * L2 + i] = accu0;
            synthesisBuffer[1 * L2 + i] = accu1;
            synthesisBuffer[2 * L2 + i] = accu2;
            synthesisBuffer[3 * L2 + i] = accu3;
            synthesisBuffer[4 * L2 + i] = accu4;
        }

        for (i = 0; i < M1; i++)
        {
            ptr_time_out[M1 - 1 - i] = synthesisBuffer[4 * L2 + M1 + i];
        }

        ptr_time_out += M1;

        synthesisBuffer -= M1;

        set_f( synthesisBuffer, 0, M1 );
    }

    return;
}

/*-------------------------------------------------------------------*
 * configureClfdb()
 *
 * configures a CLDFB handle
 *--------------------------------------------------------------------*/
static void configureCldfb(
    HANDLE_CLDFB_FILTER_BANK h_cldfb,   /* i/o : filter bank handle */
    int samplerate                      /* i   : max samplerate to oeprate */
)
{
    short k;

    h_cldfb->no_col          = CLDFB_NO_COL_MAX;
    h_cldfb->bandsToZero     = 0;
    h_cldfb->nab             = 0;

    h_cldfb->no_channels     = samplerate * INV_CLDFB_BANDWIDTH + 0.5f;
    h_cldfb->p_filter_length = 10*h_cldfb->no_channels;

    cldfb_init_proto_and_twiddles (h_cldfb);

    h_cldfb->scale = 0.f;
    for ( k=0; k<h_cldfb->p_filter_length; k++ )
    {
        h_cldfb->scale += h_cldfb->p_filter[k] * h_cldfb->p_filter[k];
    }

    h_cldfb->scale *= (float)(6400/h_cldfb->no_channels);
    h_cldfb->scale = (float)sqrt( h_cldfb->scale );

    return;
}

/*-------------------------------------------------------------------*
 * openClfdb()
 *
 * open and configures a CLDFB handle
 *--------------------------------------------------------------------*/
int openCldfb(
    HANDLE_CLDFB_FILTER_BANK *h_cldfb,   /* i/o : filter bank handle */
    CLDFB_TYPE type,                     /* i   : analysis or synthesis */
    int samplerate                       /* i   : max samplerate to oeprate */
)
{
    HANDLE_CLDFB_FILTER_BANK hs;
    short buf_len;

    hs = (HANDLE_CLDFB_FILTER_BANK) calloc(1, sizeof (CLDFB_FILTER_BANK));
    if( hs == NULL )
    {
        return (1);
    }

    hs->type = type;

    configureCldfb (hs, samplerate);
    hs->memory = NULL;
    hs->memory_length = 0;

    if (type == CLDFB_ANALYSIS)
    {
        short timeOffset = hs->p_filter_length - hs->no_channels;
        buf_len = (hs->no_channels*hs->no_col+timeOffset);
    }
    else
    {
        buf_len = (hs->p_filter_length + hs->no_channels*hs->no_col);
    }

    hs->cldfb_state = (float *) calloc( buf_len, sizeof (float));
    if (hs->cldfb_state == NULL)
    {
        return (1);
    }

    set_f(hs->cldfb_state, 0.0f, buf_len);

    *h_cldfb = hs;

    return (0);
}


/*-------------------------------------------------------------------*
* resampleCldfb()
*
* Change sample rate of filter bank
*--------------------------------------------------------------------*/
void resampleCldfb(
    HANDLE_CLDFB_FILTER_BANK hs,
    int newSamplerate
)
{
    short timeOffset, newframeSize;

    /* keep old parameters before switching*/
    int timeOffsetold   = hs->p_filter_length - hs->no_channels;
    int old_no_channels = hs->no_channels;

    /* new settings */
    configureCldfb (hs, newSamplerate);
    timeOffset = hs->p_filter_length - hs->no_channels;
    newframeSize = hs->no_channels * hs->no_col;

    /*low complexity-resampling only stored previous samples that are needed for next frame modulation */
    lerp(hs->cldfb_state+(old_no_channels*hs->no_col),hs->cldfb_state+(old_no_channels*hs->no_col),timeOffset, timeOffsetold);
    mvr2r(hs->cldfb_state+(old_no_channels*hs->no_col),hs->cldfb_state+newframeSize,timeOffset);

    return;
}

/*-------------------------------------------------------------------*
* analysisCLDFBEncoder()
*
* Encoder CLDFB analysis + energy stage
*--------------------------------------------------------------------*/

void analysisCldfbEncoder(
    Encoder_State *st,                    /* i/o: encoder state structure                    */
    const float *timeIn,
    int samplesToProcess,
    float realBuffer[CLDFB_NO_COL_MAX][CLDFB_NO_CHANNELS_MAX],
    float imagBuffer[CLDFB_NO_COL_MAX][CLDFB_NO_CHANNELS_MAX],
    float *ppBuf_Ener
)
{
    short i;
    float *ppBuf_Real[CLDFB_NO_COL_MAX];
    float *ppBuf_Imag[CLDFB_NO_COL_MAX];

    for( i=0; i<CLDFB_NO_COL_MAX; i++ )
    {
        ppBuf_Real[i] = &realBuffer[i][0];
        ppBuf_Imag[i] = &imagBuffer[i][0];
    }

    cldfbAnalysis( timeIn,ppBuf_Real,ppBuf_Imag, samplesToProcess, st->cldfbAnaEnc );

    st->currEnergyHF = GetEnergyCldfb( ppBuf_Ener, &st->currEnergyLookAhead, ppBuf_Real, ppBuf_Imag,
                                       st->cldfbAnaEnc->no_channels, st->cldfbAnaEnc->no_col, &(st->tecEnc) );

    return;
}

/*-------------------------------------------------------------------*
* GetEnergyCldfb()
*
* Conduct energy from complex data
*--------------------------------------------------------------------*/
static float GetEnergyCldfb(
    float *energyValuesSum,/*!< the result of the operation */
    float *energyLookahead, /*!< the result in the core look-ahead slot */
    float **realValues,     /*!< the real part of the subsamples */
    float **imagValues,     /*!< the imaginary part of the subsamples */
    int   numberBands,      /*!< number of bands */
    int   numberCols,       /*!< number of subsamples */
    HANDLE_TEC_ENC hTecEnc
)
{
    short j, k;
    float energyValues[CLDFB_NO_COL_MAX][CLDFB_NO_CHANNELS_MAX];
    short numLookahead = 1;

    for (k = 0; k < numberCols; k++)
    {
        for (j = 0; j < numberBands; j++)
        {
            energyValues[k][j] = realValues[k][j] * realValues[k][j] +
                                 imagValues[k][j] * imagValues[k][j];
        }
    }

    if(numberBands >= freqTable[1])
    {
        float *tempEnergyValuesArry[CLDFB_NO_COL_MAX];

        assert(numberCols == CLDFB_NO_COL_MAX);
        for (j=0; j<numberCols; j++)
        {
            tempEnergyValuesArry[j] = &energyValues[j][0];
        }

        calcHiEnvLoBuff( numberCols, freqTable, 1, tempEnergyValuesArry, hTecEnc->loBuffer, hTecEnc->hiTempEnv );
    }

    for (j = 0; j < numberBands; j++)
    {
        energyValuesSum[j]=0;
        for (k = 0; k < CLDFB_NO_COL_MAX; k++)
        {
            energyValuesSum[j] +=  energyValues[k][j];
        }
    }

    if (numberBands > 20)
    {
        float energyHF = *energyLookahead; /* energy above 8 kHz */
        numberCols -= numLookahead;
        *energyLookahead = 6.1e-5f; /* process look-ahead region */

        for (j = 20; j < min(40, numberBands); j++)
        {
            energyHF += energyValuesSum[j];

            for (k = numberCols; k < CLDFB_NO_COL_MAX; k++)
            {
                energyHF -= energyValues[k][j];
                *energyLookahead += energyValues[k][j];
            }
        }

        return energyHF * OUTMAX_SQ_INV;
    }

    return 65535.0f;
}


/*-------------------------------------------------------------------*
* GetEnergyCldfb()
*
* Remove handle
*--------------------------------------------------------------------*/
void deleteCldfb(
    HANDLE_CLDFB_FILTER_BANK * h_cldfb
)
{
    HANDLE_CLDFB_FILTER_BANK hs = *h_cldfb;

    if (hs)
    {
        if (hs->cldfb_state)
        {
            free(hs->cldfb_state);
        }
        free(hs);
        *h_cldfb = NULL;
    }

    return;
}


/*-------------------------------------------------------------------*
* cldfb_init_proto_and_twiddles()
*
* Initializes rom pointer
*--------------------------------------------------------------------*/
static void cldfb_init_proto_and_twiddles(
    HANDLE_CLDFB_FILTER_BANK hs
)
{

    /*find appropriate set of rotVecs*/
    switch(hs->no_channels)
    {
    case 10:
        hs->rot_vec_ana_re = rot_vec_ana_re_L10;
        hs->rot_vec_ana_im = rot_vec_ana_im_L10;
        hs->rot_vec_syn_re = rot_vec_syn_re_L10;
        hs->rot_vec_syn_im = rot_vec_syn_im_L10;
        hs->p_filter       = CLDFB80_10;
        break;

    case 16:
        hs->rot_vec_ana_re = rot_vec_ana_re_L16;
        hs->rot_vec_ana_im = rot_vec_ana_im_L16;
        hs->rot_vec_syn_re = rot_vec_syn_re_L16;
        hs->rot_vec_syn_im = rot_vec_syn_im_L16;
        hs->p_filter       = CLDFB80_16;
        break;

    case 20:
        hs->rot_vec_ana_re = rot_vec_ana_re_L20;
        hs->rot_vec_ana_im = rot_vec_ana_im_L20;
        hs->rot_vec_syn_re = rot_vec_syn_re_L20;
        hs->rot_vec_syn_im = rot_vec_syn_im_L20;
        hs->p_filter       = CLDFB80_20;
        break;

    case 30:
        hs->rot_vec_ana_re = rot_vec_ana_re_L30;
        hs->rot_vec_ana_im = rot_vec_ana_im_L30;
        hs->rot_vec_syn_re = rot_vec_syn_re_L30;
        hs->rot_vec_syn_im = rot_vec_syn_im_L30;
        hs->p_filter       = CLDFB80_30;
        break;

    case 32:
        hs->rot_vec_ana_re = rot_vec_ana_re_L32;
        hs->rot_vec_ana_im = rot_vec_ana_im_L32;
        hs->rot_vec_syn_re = rot_vec_syn_re_L32;
        hs->rot_vec_syn_im = rot_vec_syn_im_L32;
        hs->p_filter       = CLDFB80_32;
        break;

    case 40:
        hs->rot_vec_ana_re = rot_vec_ana_re_L40;
        hs->rot_vec_ana_im = rot_vec_ana_im_L40;
        hs->rot_vec_syn_re = rot_vec_syn_re_L40;
        hs->rot_vec_syn_im = rot_vec_syn_im_L40;
        hs->p_filter       = CLDFB80_40;
        break;

    case 60:
        hs->rot_vec_ana_re = rot_vec_ana_re_L60;
        hs->rot_vec_ana_im = rot_vec_ana_im_L60;
        hs->rot_vec_syn_re = rot_vec_syn_re_L60;
        hs->rot_vec_syn_im = rot_vec_syn_im_L60;
        hs->p_filter       = CLDFB80_60;
        break;

    }

    return;
}


/*-------------------------------------------------------------------*
* cldfb_save_memory()
*
* Save the memory of filter; to be restored with cldfb_restore_memory()
*--------------------------------------------------------------------*/
int cldfb_save_memory(
    HANDLE_CLDFB_FILTER_BANK hs   /* i/o: cldfb handle */
)
{
    unsigned int offset = hs->p_filter_length - hs->no_channels;
    unsigned int frameSize = hs->no_channels * hs->no_col;

    if (hs->memory != NULL || hs->memory_length!=0)
    {
        /* memory already stored; Free memory first */
        return 1;
    }


    if (hs->type == CLDFB_ANALYSIS)
    {
        hs->memory_length = offset + frameSize;
    }
    else
    {
        hs->memory_length = hs->p_filter_length;
    }

    hs->memory = (float *) calloc( hs->memory_length, sizeof (float));
    if (hs->memory == NULL)
    {
        /* memory cannot be allocated */
        return (1);
    }

    /* save the memory */
    mvr2r (hs->cldfb_state, hs->memory, hs->memory_length);

    return 0;
}


/*-------------------------------------------------------------------*
* cldfb_restore_memory()
*
* Restores the memory of filter; memory to be save by cldfb_save_memory()
*--------------------------------------------------------------------*/
int cldfb_restore_memory(
    HANDLE_CLDFB_FILTER_BANK hs   /* i/o: cldfb handle */
)
{
    unsigned int offset = hs->p_filter_length - hs->no_channels;
    unsigned int frameSize = hs->no_channels * hs->no_col;
    unsigned int size;

    if (hs->memory == NULL || hs->memory_length==0)
    {
        /* memory not allocated */
        return 1;
    }

    if ( hs->type == CLDFB_ANALYSIS )
    {
        size = offset + frameSize;
    }
    else
    {
        size = hs->p_filter_length;
    }

    /* read the memory */
    mvr2r (hs->memory, hs->cldfb_state, hs->memory_length);

    /* adjust sample rate if it was changed in the meanwhile */
    if (hs->memory_length != size)
    {
        lerp(hs->cldfb_state, hs->cldfb_state, size, hs->memory_length);
    }

    hs->memory_length = 0;
    free(hs->memory);
    hs->memory = NULL;

    return 0;
}


/*-------------------------------------------------------------------*
* cldfb_reset_memory()
*
* Resets the memory of filter.
*--------------------------------------------------------------------*/
int cldfb_reset_memory(
    HANDLE_CLDFB_FILTER_BANK hs     /* i/o: cldfb handle */
)
{
    unsigned int offset = hs->p_filter_length - hs->no_channels;
    unsigned int frameSize = hs->no_channels * hs->no_col;
    int memory_length;

    if (hs->type == CLDFB_ANALYSIS)
    {
        memory_length = offset + frameSize;
    }
    else
    {
        memory_length = hs->p_filter_length;
    }

    /* save the memory */
    set_f (hs->cldfb_state, 0, memory_length);

    return 0;
}
