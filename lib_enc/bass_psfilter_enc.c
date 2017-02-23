/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "options.h"

#include "prot.h"
#include "rom_com.h"



/*---------------------------------------------------------------------*
 * bass_pf_enc()
 *
 * Low-frequency postfiltering, decoder parammeter estimation
 *---------------------------------------------------------------------*/

float bass_pf_enc(
    const float *orig,            /* (i) : 12.8kHz original signal                     */
    const float *syn,             /* (i) : 12.8kHz synthesis to postfilter             */
    const float pitch_buf[],      /* (i) : Pitch gain for all subframes (gainT_sf[16]) */
    const float gainT_sf[],       /* (i) : Pitch gain for all subframes (gainT_sf[16]) */
    const short l_frame,          /* (i) : frame length (should be multiple of l_subfr)*/
    const short l_subfr_in,       /* (i) : sub-frame length (80/64)                    */
    float mem_bpf[],        /* i/o : memory state [2*L_FILT16k]                  */
    float mem_error_bpf[],  /* i/o : memory state [2*L_FILT16k]                  */
    int *gain_factor_param, /* (o) : quantized gain factor                       */
    const short mode,             /* (i) : coding mode of adapt bpf                    */
    float *mem_deemph_err,  /* o  : Error deemphasis memory                      */
    float *lp_ener          /* o  : long_term error signal energy                */
)
{
    int i, j, sf, i_subfr, T, lg, l_subfr, l_filt;
    float d,n, snr, nrg1,nrg2, gain,nrg,tmp;
    float noise_buf[L_FILT16k+(2*L_SUBFR)], *noise, *noise_in;
    float error_buf[L_FILT16k+(2*L_SUBFR)], *error, *error_in;
    float cross_n_d, nrg_n;
    const float *pFilt;
    float ener2;


    if ( l_frame!=L_FRAME16k )
    {
        pFilt = filt_lp;
        l_filt = L_FILT;
    }
    else
    {
        pFilt = filt_lp_16kHz;
        l_filt = L_FILT16k;
    }

    noise = noise_buf + l_filt;
    noise_in = noise_buf + 2*l_filt;
    error = error_buf + l_filt;
    error_in = error_buf + 2*l_filt;

    sf = 0;
    snr=0.f;
    nrg_n=1e-6f;
    cross_n_d = 0.f;
    l_subfr = l_subfr_in;
    for (i_subfr=0; i_subfr<l_frame; i_subfr+=l_subfr, sf++)
    {
        T = (int)pitch_buf[sf];
        gain = gainT_sf[sf];

        if (gain > 1.0f) gain = 1.0f;
        if (gain < 0.0f) gain = 0.0f;

        lg = l_frame - T - i_subfr;
        if (lg < 0) lg = 0;
        if (lg > l_subfr) lg = l_subfr;


        if (gain > 0)
        {
            tmp = 0.01f;
            nrg=0.01f;
            for (i=0; i<lg; i++)
            {
                tmp += syn[i+i_subfr] * (0.5f*syn[i+i_subfr-T] + 0.5f*syn[i+i_subfr+T]);
                nrg += (0.5f*syn[i+i_subfr-T] + 0.5f*syn[i+i_subfr+T])*(0.5f*syn[i+i_subfr-T] + 0.5f*syn[i+i_subfr+T]);
            }
            for (i=lg; i<l_subfr; i++)
            {
                tmp += syn[i+i_subfr]*syn[i+i_subfr-T];
                nrg += syn[i+i_subfr-T]*syn[i+i_subfr-T];
            }
            gain=tmp/nrg;

            if(gain>1.0f)
            {
                gain=1.0f;
            }
            else if(gain<0.f)
            {
                gain=0.f;
            }

            ener2=0.01f;
            for (i=0; i<lg; i++)
            {
                error[i] = gain * (syn[i+i_subfr] - 0.5f*syn[i+i_subfr-T] - 0.5f*syn[i+i_subfr+T]);
                error[i] = error[i] + 0.9f* *mem_deemph_err;
                *mem_deemph_err = error[i];
                ener2 +=error[i]*error[i];
            }
            for (i=lg; i<l_subfr; i++)
            {
                error[i] = 0.5f*gain * (syn[i+i_subfr] - syn[i+i_subfr-T]);
                error[i] = error[i] + 0.9f* *mem_deemph_err;
                *mem_deemph_err = error[i];
                ener2 +=error[i]*error[i];
            }

            ener2=(float) (10.f*log10(ener2));
            *lp_ener = (float)(0.99f* *lp_ener + 0.01f*ener2);
            ener2=(float)pow(10.f,0.1f* *lp_ener);

            tmp=0.5f*tmp/(nrg+ener2);
            if(tmp>0.5f)
            {
                tmp=0.5f;
            }
            else if(tmp<0.f)
            {
                tmp=0.0f;
            }

            for (i=0; i<lg; i++)
            {
                noise_in[i] = tmp * (syn[i+i_subfr] - 0.5f*syn[i+i_subfr-T] - 0.5f*syn[i+i_subfr+T]);
                error_in[i] =(orig[i+i_subfr]-syn[i+i_subfr]);
            }
            for (i=lg; i<l_subfr; i++)
            {
                noise_in[i] = tmp * (syn[i+i_subfr] - syn[i+i_subfr-T]);
                noise_in[i] *= 0.5f;
                error_in[i] =(orig[i+i_subfr]-syn[i+i_subfr]);
            }
        }
        else
        {
            set_zero(noise_in, l_subfr);
            set_zero(error_in, l_subfr);
        }

        mvr2r(mem_bpf, noise_buf, 2*l_filt);
        mvr2r(noise_buf+l_subfr, mem_bpf, 2*l_filt);

        mvr2r(mem_error_bpf, error_buf, 2*l_filt);
        mvr2r(error_buf+l_subfr, mem_error_bpf, 2*l_filt);

        nrg1=1e-6f;
        nrg2=1e-6f;

        /* substract from voiced speech low-pass filtered noise */
        for (i=0; i<l_subfr; i++)
        {
            n = pFilt[0] * noise[i];
            d =  error[i];

            for(j=1; j<=l_filt; j++)
            {
                n += pFilt[j] * (noise[i-j] + noise[i+j]);
            }
            /*for optimal g*/
            nrg_n += n*n;
            cross_n_d += n*d;

            /*for evaluating SNR*/
            nrg1 += (d+n)*(d+n);
            nrg2 += d*d;
        }

        /*SegSNR*/
        snr += (float)log10( nrg2/nrg1 );
    }

    /*Compute and quantize optimal gain*/
    /* optimal gain = -<n,d>/<n,n> */
    if(mode==2)
    {
        *gain_factor_param = (int)(-2.f*(cross_n_d/nrg_n)+0.5f);
        if(*gain_factor_param>3)
        {
            *gain_factor_param=3;
        }
        else if(*gain_factor_param<0)
        {
            *gain_factor_param=0;
        }
        /*If optimal gain negatif or zero but snr still positif->gain=0.5f*/
        if(snr>0.f && *gain_factor_param==0)
        {
            *gain_factor_param=1;
        }
    }
    else
    {
        *gain_factor_param=2;
    }


    return(snr);

}
