/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"


/*-------------------------------------------------------------------*
 * modify_Fs()
 *
 * Function for resampling of signals
 *-------------------------------------------------------------------*/

short modify_Fs(           /* o  : length of output    */
    const float sigIn[],   /* i  : signal to decimate  */
    short lg,        /* i  : length of input + 0 delay signaling   */
    const int   fin,       /* i  : frequency of input  */
    float sigOut[],  /* o  : decimated signal    */
    const int   fout,      /* i  : frequency of output */
    float mem[]      /* i/o: filter memory       */
    ,int   nblp       /* i  : flag indicating if NB low-pass is applied */
)
{
    short i;
    short lg_out, fac_num, fac_den, filt_len, frac, mem_len;
    float num_den;
    short datastep, fracstep;
    float *sigIn_ptr;
    float signal_tab[3*L_FILT_MAX + L_FRAME48k], *signal, *signal_ana; /* 3* as 2* for memory and 1* for future prediction */
    float A[M+1], r[M+1], epsP[M+1], val;
    short mem_len_ana;
    short plus_sample_in;
    short j;
    float mu_preemph;
    float mem_preemph;
    const Resampling_cfg *cfg_ptr;
    short flag_low_order = 0;
    short filt_len_tmp;

    /*-------------------------------------------------------------------*
     * IIR filters for resampling to/from 8 kHz
     *-------------------------------------------------------------------*/

    /*-------------------------------------------------------------------*
     * Find the resampling configuration
     *-------------------------------------------------------------------*/

    if ( fin == fout )
    {
        /* just copy the signal and quit */
        for (i = 0; i < lg; i++)
        {
            sigOut[i] = sigIn[i];
        }

        return lg;
    }
    else
    {
        /* find the resampling configuration in the lookup table */
        for (cfg_ptr = &resampling_cfg_tbl[0]; (cfg_ptr->fin != 0) && !(cfg_ptr->fin == fin && cfg_ptr->fout == fout); cfg_ptr++)
        {
        }


        /* find config with NB 4kHz low-pass */
        if (nblp && (fin > 8000) && (fout == 12800))
        {
            flag_low_order = 1;
            for (cfg_ptr++; (cfg_ptr->fin != 0) && !(cfg_ptr->fin == fin && cfg_ptr->fout == fout); cfg_ptr++)
            {
            }
        }

        /* Retrieve and/or calculate the resampling parameters */
        fac_num = cfg_ptr->fac_num;
        fac_den = (short)((cfg_ptr->fin * fac_num) / cfg_ptr->fout);
        lg_out = (lg * fac_num) / fac_den;
        filt_len = cfg_ptr->filt_len;

        mem_len = 2*filt_len;
        plus_sample_in = 0; /* default, regular delay */
        frac = 0;

        if ( fin == 8000 && fout == 12800 )
        {
            plus_sample_in = 7;
            frac = 4;
        }

        signal = signal_tab+2*L_FILT_MAX + L_FRAME48k - mem_len - lg;
        signal_ana = signal;
        mem_len_ana = mem_len;

    }

    /*-------------------------------------------------------------------*
     * FIR filters for resampling to/from 12.8, 16, 32, 48 kHz
     *-------------------------------------------------------------------*/

    /* append filter memory */
    for (i=0; i<2*filt_len; i++)
    {
        signal[i] = mem[i];
    }

    for (i=0; i<lg; i++)
    {
        signal[i+(2*filt_len)] = sigIn[i];
    }

    if( plus_sample_in > 0 )
    {

        autocorr( signal_ana+mem_len_ana+lg-LEN_WIN_SSS, r, 1, LEN_WIN_SSS, wind_sss, 0, 0, 0 );

        mu_preemph = r[1] / r[0];
        mem_preemph = signal_ana[mem_len_ana+lg-LEN_WIN_SSS - 1];
        preemph( signal_ana+mem_len_ana+lg-LEN_WIN_SSS, mu_preemph, LEN_WIN_SSS, &mem_preemph );

        /* Autocorrelations */
        autocorr( signal_ana+mem_len_ana+lg-LEN_WIN_SSS, r, M, LEN_WIN_SSS, wind_sss, 0, 0, 0 );

        lag_wind( r, M, (float)fin, LAGW_STRONG );

        /* Levinson-Durbin */
        lev_dur( A, r, M, epsP );

        for (i=0; i<plus_sample_in; i++)
        {
            val = 0;
            for(j = 1; j <= M; j++)
            {
                val -= signal[i+lg+mem_len-j] * A[j];
            }
            signal[i+lg+mem_len] = val; /* AZ ringing padding */
        }

        mem_preemph = signal[mem_len+lg-LEN_WIN_SSS - 1];
        deemph( signal+mem_len+lg-LEN_WIN_SSS, mu_preemph, LEN_WIN_SSS+plus_sample_in, &mem_preemph );
    }

    /* interpolation */
    datastep = fac_den / fac_num;
    fracstep = fac_den - datastep * fac_num; /* equivalent to datastep = fac_den % fac_num */

    sigIn_ptr = signal + filt_len + plus_sample_in;

    filt_len_tmp = filt_len;
    if( flag_low_order )
    {
        filt_len_tmp = (filt_len + 1) >> 1;
    }

    for(i=0; i<lg_out; i++)
    {

        sigOut[i] = interpolation( sigIn_ptr, cfg_ptr->filter, frac, fac_num, filt_len_tmp );

        frac = frac + fracstep;
        if ( frac >= fac_num )
        {
            frac = frac - fac_num;
            sigIn_ptr++;
        }

        sigIn_ptr += datastep;
    }

    /* rescaling */
    if ( (fac_num > fac_den) == ((cfg_ptr->flags & RS_INV_FAC) != 0) )
    {
        num_den = (float)fac_num / fac_den;

        for( i=0; i<lg_out; i++ )
        {
            sigOut[i] *= num_den;
        }
    }

    /* update the filter memory */
    for (i=0; i<2*filt_len; i++)
    {
        mem[i] = signal[i+lg];
    }

    return lg_out;
}

short modify_Fs_intcub3m_sup(       /* o  : length of output    */
    const float sigIn[],          /* i  : signal to decimate with memory of 2 samples (indexes -2 & -1) */
    const short lg,               /* i  : length of input (suppose that lg is such that lg_out is integer, ex multiple of 5 in case of 16kHz to 12.8 kHz) */
    const int   fin,              /* i  : frequency of input  */
    float sigOut[],         /* o  : decimated signal    */
    const int   fout,             /* i  : frequency of output */
    short *delayout         /* o  : delay of output */
)
{
    short i, j, k, i1, i2, k1, k2, k3, kk, cind;
    short lg_out, fk1, k2d, k3d;
    float cc[4][4];
    float vv;
    const double (*cu)[3] = 0;

    /*-------------------------------------------------------------------*
     * Find the resampling configuration
     *-------------------------------------------------------------------*/

    /* check if fin and fout are the same */
    if ( fin == fout )
    {
        /* just copy the signal and quit */
        for (i = 0; i < lg; i++)
        {
            sigOut[i] = sigIn[i];
        }

        *delayout = 0;
        return lg;
    }
    else
    {
        /* length of the interpolated signal */
        lg_out = (short)(lg * fout / fin);

        /* cc[x][3]*s*s*s + cc[x][2]*s*s + cc[x][1]*s + cc[x][0]; indexes relatives of s : -1 0 1 2 */
        /* d : cc[x][0] = s[0] */
        /* b : cc[x][2] =(s[-1]+s[1])/2-s[0] */
        /* a : cc[x][3] = (s[-1]+s[2]-s[0]-s[1]-4*cc[x][2]) / 6 */
        /* c : cc[x][1] = s[1]-s[0]-cc[x][3]-cc[x][2] */

        /* coef inits using memory (indexes < 0) */
        /* cc[2][] : indexes -2 -1 0 1 */
        cc[2][0] = sigIn[-1]/3;
        cc[2][2] =(sigIn[-2]+sigIn[0])/2-sigIn[-1];
        cc[2][3] = (sigIn[-2]+sigIn[1]-sigIn[-1]-sigIn[0]-4*cc[2][2]) / 6;
        cc[2][1] = sigIn[0]-sigIn[-1]-cc[2][3]-cc[2][2];

        /* cc[3][] : indexes -1 0 1 2 */
        cc[3][0] = sigIn[0]/3;
        cc[3][2] =(sigIn[-1]+sigIn[1])/2-sigIn[0];
        cc[3][3] = (sigIn[-1]+sigIn[2]-sigIn[0]-sigIn[1]-4*cc[3][2]) / 6;
        cc[3][1] = sigIn[1]-sigIn[0]-cc[3][3]-cc[3][2];
        j = 0;

        if( fin == 12800 )
        {
            if( fout == 8000 )
            {
                cind = 0;
            }
            else if( fout == 16000 )
            {
                cind = 1;
            }
            else if(fout == 32000)
            {
                cind = 2;
            }
            else if(fout == 48000)
            {
                cind = 3;
            }
            else
            {
                printf("warning, output sampling frequency %d not implemented for input %d", fout, fin);
                return(-1);
            }
        }
        else if(fin == 16000)
        {
            if(fout == 12800)
            {
                cind = 4;
            }
            else if(fout == 32000)
            {
                cind = 5;
            }
            else if(fout == 48000)
            {
                cind = 6;
            }
            else
            {
                printf("warning, output sampling frequency %d not implemented for input %d", fout, fin);
                return(-1);
            }
        }
        else
        {
            printf("warning, input sampling frequency %d not implemented", fin);
            return(-1);
        }

        *delayout = ct2[cind][9];

        if( ct2[cind][12] == 15 )
        {
            cu = cu15;
        }

        if( ct2[cind][12] == 4 )
        {
            cu = cu4;
        }

        fk1 = 2*ct2[cind][12] - 2;
        k2d = fk1 / 2; /* shift of index in cu with respect to the next sample (ex 1.25 -> 0.25 ) */
        k3d = fk1 - 1; /* to compurte index in cu with respect to the last sample with - sign (ex 1.25 -> -0.75 ) */

        kk = 0;
        for(i = 0; i < lg-ct2[cind][11];)
        {
            sigOut[j++] = sigIn[i];
            for(k = 0; k < ct2[cind][10]; k++)
            {
                cc[kk][0] = sigIn[i+1]/3;
                cc[kk][2] =(sigIn[i]+sigIn[i+2])/2-sigIn[i+1];
                cc[kk][3] = (sigIn[i]+sigIn[i+3]-sigIn[i+1]-sigIn[i+2]-4*cc[kk][2]) / 6;
                cc[kk][1] = sigIn[i+2]-sigIn[i+1]-cc[kk][3]-cc[kk][2];
                i++;

                i2 = kk-2;
                i1 = kk-1;
                if( i1 < 0 )
                {
                    i1 += 4;
                }

                if( i2 < 0 )
                {
                    i2 += 4;
                }

                for(k1 = ct2[cind][k]; k1 < fk1; k1 += ct2[cind][8])
                {
                    k2 = k1 - k2d;
                    k3 = k3d - k1;
                    vv  = (float)( cu[k1][2]*cc[i2][3] + cu[k1][1]*cc[i2][2]  + cu[k1][0]*cc[i2][1] + cc[i2][0]);
                    vv += (float)( cu[k2][2]*cc[i1][3] + cu[k2][1]*cc[i1][2]  + cu[k2][0]*cc[i1][1] + cc[i1][0]);
                    vv += (float)(-cu[k3][2]*cc[kk][3] + cu[k3][1]*cc[kk][2]  - cu[k3][0]*cc[kk][1] + cc[kk][0]);
                    sigOut[j++] = vv;
                }

                kk++;
                if( kk == 4 )
                {
                    kk = 0;
                }
            }
        }

        sigOut[j++] = sigIn[i];

        for(k = 0; k < ct2[cind][11]-3; k++)
        {
            cc[kk][0] = sigIn[i+1]/3;
            cc[kk][2] =(sigIn[i]+sigIn[i+2])/2-sigIn[i+1];
            cc[kk][3] = (sigIn[i]+sigIn[i+3]-sigIn[i+1]-sigIn[i+2]-4*cc[kk][2]) / 6;
            cc[kk][1] = sigIn[i+2]-sigIn[i+1]-cc[kk][3]-cc[kk][2];
            i++;

            i2 = kk-2;
            i1 = kk-1;
            if( i1 < 0 )
            {
                i1 += 4;
            }

            if( i2 < 0 )
            {
                i2 += 4;
            }

            for(k1 = ct2[cind][k]; k1 < fk1; k1 += ct2[cind][8])
            {
                k2 = k1 - k2d;
                k3 = k3d - k1;
                vv  = (float)( cu[k1][2]*cc[i2][3] + cu[k1][1]*cc[i2][2]  + cu[k1][0]*cc[i2][1] + cc[i2][0]);
                vv += (float)( cu[k2][2]*cc[i1][3] + cu[k2][1]*cc[i1][2]  + cu[k2][0]*cc[i1][1] + cc[i1][0]);
                vv += (float)(-cu[k3][2]*cc[kk][3] + cu[k3][1]*cc[kk][2]  - cu[k3][0]*cc[kk][1] + cc[kk][0]);
                sigOut[j++] = vv;
            }

            kk++;

            if( kk == 4 )
            {
                kk = 0;
            }
        }

        kk--;
        if( kk == -1 )
        {
            kk = 3;
        }

        if( ct2[cind][10] == 1 )
        {
            sigOut[j++] = sigIn[i];
        }

        for(k1 = ct2[cind][k]; k1 < fk1; k1 += ct2[cind][8])
        {
            k2 = k1 - k2d;
            vv = (float)( cu[k2][2]*cc[kk][3] + cu[k2][1]*cc[kk][2]  + cu[k2][0]*cc[kk][1] + cc[kk][0]);
            sigOut[j++] = vv*3;
        }

        if( ct2[cind][10] < 3 )
        {
            sigOut[j++] = sigIn[i+1];
        }

        for( k1 = ct2[cind][k+1]; k1 < fk1; k1 += ct2[cind][8] )
        {
            vv = (float)( cu[k1][2]*cc[kk][3] + cu[k1][1]*cc[kk][2]  + cu[k1][0]*cc[kk][1] + cc[kk][0]);
            sigOut[j++] = vv*3;
        }

        if( ct2[cind][10] == 1 )
        {
            sigOut[j++] = sigIn[i+2];
        }
    }

    return lg_out;
}

/*-------------------------------------------------------------------*
 * Interpolate_allpass_steep()
 *
 * Interpolation by a factor 2
 *-------------------------------------------------------------------*/

void Interpolate_allpass_steep(
    const float *in,            /* i  : input array of size N       */
    float *state,            /* i/o: memory                      */
    const short N,                /* i  : number of input samples     */
    float *out            /* o  : output array of size 2*N    */
)
{
    short n, k;
    float temp[ALLPASSSECTIONS_STEEP - 1];

    /* upper allpass filter chain */
    for (k = 0; k < N; k++)
    {
        temp[0] = state[0] + AP2_STEEP[0] * in[k];
        state[0] = in[k] - AP2_STEEP[0] * temp[0];

        /* for better performance, unroll this loop */
        for (n = 1; n < ALLPASSSECTIONS_STEEP - 1; n++)
        {
            temp[n] = state[n] + AP2_STEEP[n] * temp[n - 1];
            state[n] = temp[n - 1] - AP2_STEEP[n] * temp[n];
        }

        out[2 * k + 1] = state[ALLPASSSECTIONS_STEEP - 1] + AP2_STEEP[ALLPASSSECTIONS_STEEP - 1] * temp[ALLPASSSECTIONS_STEEP - 2];
        state[ALLPASSSECTIONS_STEEP - 1] = temp[ALLPASSSECTIONS_STEEP - 2] - AP2_STEEP[ALLPASSSECTIONS_STEEP - 1] * out[2 * k + 1];
    }

    /* lower allpass filter chain */
    for (k = 0; k < N; k++)
    {
        temp[0] = state[ALLPASSSECTIONS_STEEP] + AP1_STEEP[0] * in[k];
        state[ALLPASSSECTIONS_STEEP] = in[k] - AP1_STEEP[0] * temp[0];

        /* for better performance, unroll this loop */
        for (n = 1; n < ALLPASSSECTIONS_STEEP - 1; n++)
        {
            temp[n] = state[ALLPASSSECTIONS_STEEP + n] + AP1_STEEP[n] * temp[n - 1];
            state[ALLPASSSECTIONS_STEEP + n] = temp[n - 1] - AP1_STEEP[n] * temp[n];
        }

        out[2 * k] = state[2 * ALLPASSSECTIONS_STEEP - 1] + AP1_STEEP[ALLPASSSECTIONS_STEEP - 1] * temp[ALLPASSSECTIONS_STEEP - 2];
        state[2 * ALLPASSSECTIONS_STEEP - 1] = temp[ALLPASSSECTIONS_STEEP - 2] - AP1_STEEP[ALLPASSSECTIONS_STEEP - 1] * out[2 * k];
    }

    return;
}

/*-------------------------------------------------------------------*
 * Decimate_allpass_steep()
 *
 * Decimation by a factor 2
 *-------------------------------------------------------------------*/

void Decimate_allpass_steep (
    const float *in,        /* i  : input array of size N                   */
    float *state,     /* i/o: memory                                  */
    const short N,          /* i  : number of input samples                 */
    float *out        /* o  : output array of size N/2                */
)
{
    short n, k;
    float temp[ALLPASSSECTIONS_STEEP];

    /* upper allpass filter chain */
    for (k = 0; k < N / 2; k++)
    {
        temp[0] = state[0] + AP1_STEEP[0] * in[2 * k];
        state[0] = in[2 * k] - AP1_STEEP[0] * temp[0];

        /* for better performance, unroll this loop */
        for (n = 1; n < ALLPASSSECTIONS_STEEP - 1; n++)
        {
            temp[n] = state[n] + AP1_STEEP[n] * temp[n - 1];
            if( fabs(temp[n]) < 1e-12 )
            {
                temp[n] = sign(temp[n])*1e-12;
            }
            state[n] = temp[n - 1] - AP1_STEEP[n] * temp[n];
        }

        out[k] = state[ALLPASSSECTIONS_STEEP - 1] + AP1_STEEP[ALLPASSSECTIONS_STEEP - 1] * temp[ALLPASSSECTIONS_STEEP - 2];
        state[ALLPASSSECTIONS_STEEP - 1] = temp[ALLPASSSECTIONS_STEEP - 2] - AP1_STEEP[ALLPASSSECTIONS_STEEP - 1] * out[k];
    }

    /* lower allpass filter chain */
    temp[0] = state[ALLPASSSECTIONS_STEEP] + AP2_STEEP[0] * state[2 * ALLPASSSECTIONS_STEEP];
    state[ALLPASSSECTIONS_STEEP] = state[2 * ALLPASSSECTIONS_STEEP] - AP2_STEEP[0] * temp[0];

    /* for better performance, unroll this loop */
    for (n = 1; n < ALLPASSSECTIONS_STEEP - 1; n++)
    {
        temp[n] = state[ALLPASSSECTIONS_STEEP + n] + AP2_STEEP[n] * temp[n - 1];
        if( fabs(temp[n]) < 1e-12 )
        {
            temp[n] = sign(temp[n])*1e-12;
        }
        state[ALLPASSSECTIONS_STEEP + n] = temp[n - 1] - AP2_STEEP[n] * temp[n];
    }

    temp[ALLPASSSECTIONS_STEEP - 1] = state[2 * ALLPASSSECTIONS_STEEP - 1] + AP2_STEEP[ALLPASSSECTIONS_STEEP - 1] *
                                      temp[ALLPASSSECTIONS_STEEP -2];

    state[2 * ALLPASSSECTIONS_STEEP - 1] = temp[ALLPASSSECTIONS_STEEP - 2] - AP2_STEEP[ALLPASSSECTIONS_STEEP - 1] *
                                           temp[ALLPASSSECTIONS_STEEP - 1];
    out[0] = (float)((out[0] + temp[ALLPASSSECTIONS_STEEP - 1]) * 0.5);

    for (k = 1; k < N / 2; k++)
    {
        temp[0] = state[ALLPASSSECTIONS_STEEP] + AP2_STEEP[0] * in[2 * k - 1];
        state[ALLPASSSECTIONS_STEEP] = in[2 * k - 1] - AP2_STEEP[0] * temp[0];

        /* for better performance, unroll this loop */
        for (n = 1; n < ALLPASSSECTIONS_STEEP - 1; n++)
        {
            temp[n] = state[ALLPASSSECTIONS_STEEP + n] + AP2_STEEP[n] * temp[n - 1];
            if( fabs(temp[n]) < 1e-12 )
            {
                temp[n] = sign(temp[n])*1e-12;
            }
            state[ALLPASSSECTIONS_STEEP + n] = temp[n - 1] - AP2_STEEP[n] * temp[n];
        }

        temp[ALLPASSSECTIONS_STEEP - 1] = state[2 * ALLPASSSECTIONS_STEEP - 1] + AP2_STEEP[ALLPASSSECTIONS_STEEP - 1] *
                                          temp[ALLPASSSECTIONS_STEEP - 2];
        state[2 * ALLPASSSECTIONS_STEEP - 1] = temp[ALLPASSSECTIONS_STEEP - 2] - AP2_STEEP[ALLPASSSECTIONS_STEEP - 1] *
                                               temp[ALLPASSSECTIONS_STEEP - 1];
        out[k] = (float)((out[k] + temp[ALLPASSSECTIONS_STEEP - 1]) * 0.5);
    }

    /* z^(-1) */
    state[2 * ALLPASSSECTIONS_STEEP] = in[N - 1];

    return;
}

/*-------------------------------------------------------------------*
 * interpolate_3_over_2_allpass()
 *
 * Interpolate 3/2 using allpass iir polyphase filter. Delay 4 samples @48k
 *-------------------------------------------------------------------*/

void interpolate_3_over_2_allpass(
    const float *input,           /* i  : input signal            */
    const short len,              /* i  : number of input samples */
    float *out,             /* o  : output signal           */
    float *mem,             /* i/o: memory                  */
    const float *filt_coeff       /* i  : filter coefficients     */
)
{
    short i, loop_len;
    float Vu[2], Vm[2], Vl[2];    /* Outputs of three cascaded allpass stages (upper, middle, and lower) */
    float out1_buff[L_FRAME32k*3];
    float * out1;
    float mem_temp;

    out1 = out1_buff;

    for (i = 0; i < len; i++ )
    {
        /* Upper branch */
        Vu[0] = mem[0] + filt_coeff[0] * ( input[i] - mem[1] );
        Vu[1] = mem[1] + filt_coeff[1] * ( Vu[0] - mem[2] );
        mem[3] = mem[2] + filt_coeff[2] * ( Vu[1] - mem[3] );

        mem[1] = Vu[0];
        mem[2] = Vu[1];
        *out1++ = mem[3];

        /* Middle branch */
        Vm[0] = mem[0] + filt_coeff[3] * (input[i]-mem[4]);
        Vm[1] = mem[4] + filt_coeff[4] * (Vm[0]-mem[5]);
        mem[6] = mem[5] + filt_coeff[5] * (Vm[1]-mem[6]);

        mem[4] = Vm[0];
        mem[5] = Vm[1];
        *out1++ = mem[6];

        /* Lower branch */
        Vl[0] = mem[0] + filt_coeff[6] * (input[i]-mem[7]);
        Vl[1] = mem[7] + filt_coeff[7] * (Vl[0]-mem[8]);
        mem[9] = mem[8] + filt_coeff[8] * (Vl[1]-mem[9]);

        mem[0] = input[i];
        mem[7] = Vl[0];
        mem[8] = Vl[1];
        *out1++ = mem[9];
    }

    loop_len = len*3/2;

    /*decimate by 2 and LPF*/
    for(i = 0; i < loop_len; i++)
    {
        mem_temp = out1_buff[2*i];
        out[i] = (((0.0473147f)*(mem_temp+mem[10]))+((-0.151521f)*(mem[11]+mem[14])));
        out[i] = (out[i]+((0.614152f)*(mem[12]+mem[13])));
        mem[10] = mem[11];
        mem[11] = mem[12];
        mem[12] = mem[13];
        mem[13] = mem[14];
        mem[14] = mem_temp;
    }
    return;
}

/*-------------------------------------------------------------------*
* decimate_2_over_3_allpass()
*
* Decimate 2/3 using allpass iir polyphase filter.
*-------------------------------------------------------------------*/

void decimate_2_over_3_allpass(
    const float *input,           /* i  : input signal            */
    const short len,              /* i  : number of input samples */
    float *out,             /* o  : output signal           */
    float *mem,             /* i/o: memory                  */
    const float *filt_coeff,      /* i  : filter coefficients     */
    const float *lp_num,
    const float *lp_den,
    float *lp_mem
)
{
    short i, loop_len;
    float Vu[2], Vm[2], Vl[2];    /* Outputs of three cascaded allpass stages (upper, middle, and lower) */
    float * out1;
    float *in;
    float out1_buff[L_FRAME48k*2];
    float tmp;

    /* Combine the 2nd order iir lpf with the decimation by 2 to improve the efficiency*/
    out1 = out1_buff;

    *out1++ = lp_num[0] * ( input[0] + lp_mem[0] ) - lp_den[2] * lp_mem[2];
    *out1++ = lp_num[1] * input[0] - lp_den[2] * lp_mem[1];

    for (i=1; i < len; i++)
    {
        tmp = lp_num[0] * ( input[i] + input[i-1] ) - lp_den[2] * out1[-2];
        *out1++ = tmp;
        tmp = lp_num[1] * input[i] - lp_den[2] * out1[-2];
        *out1++ = tmp;
    }
    lp_mem[0] = input[len-1];
    lp_mem[1] = out1[-1];
    lp_mem[2] = out1[-2];

    /* do the all pass polyphase filter with pi/3 cutoff */
    out1 = out;
    in = out1_buff;
    loop_len = (short) len*2/3;

    for (i = 0; i < loop_len; i++ )
    {
        /* Lower branch */
        Vl[0] = mem[8] + filt_coeff[6] * (*in - mem[9]);
        Vl[1] = mem[9] + filt_coeff[7] * (Vl[0] - mem[10]);
        mem[11] = mem[10] + filt_coeff[8] * (Vl[1] - mem[11]);

        mem[8] = *in++;
        mem[9] = Vl[0];
        mem[10] = Vl[1];
        *out1 = mem[11];

        /* Middle branch */
        Vm[0] = mem[4] + filt_coeff[3] * (*in - mem[5]);
        Vm[1] = mem[5] + filt_coeff[4] * (Vm[0]-mem[6]);
        mem[7] = mem[6] + filt_coeff[5] * (Vm[1]-mem[7]);

        mem[4] = *in++;
        mem[5] = Vm[0];
        mem[6] = Vm[1];
        *out1 += mem[7];

        /* Upper branch */
        Vu[0] = mem[0] + filt_coeff[0] * ( *in - mem[1] );
        Vu[1] = mem[1] + filt_coeff[1] * ( Vu[0] - mem[2] );
        mem[3] = mem[2] + filt_coeff[2] * ( Vu[1] - mem[3] );

        mem[0] = *in++;
        mem[1] = Vu[0];
        mem[2] = Vu[1];
        *out1++ += mem[3];
    }

    return;
}

/*-------------------------------------------------------------------*
  * interpolate_3_over_1_allpass()
  *
  * Interpolate 3/1 using allpass iir polyphase filter. Delay 4 samples @48k
  *-------------------------------------------------------------------*/

void interpolate_3_over_1_allpass(
    const float *input,           /* i  : input signal            */
    const short len,              /* i  : number of input samples */
    float *out,             /* o  : output signal           */
    float *mem,             /* i/o: memory                  */
    const float *filt_coeff       /* i  : filter coefficients     */
)
{
    short i;
    float Vu[2], Vm[2], Vl[2];    /* Outputs of three cascaded allpass stages (upper, middle, and lower) */
    float * out1;
    float mem_temp;
    out1 = &out[0];

    for (i = 0; i < len; i++ )
    {
        /* Upper branch */
        Vu[0] = mem[0] + filt_coeff[0] * ( input[i] - mem[1] );
        Vu[1] = mem[1] + filt_coeff[1] * ( Vu[0] - mem[2] );
        mem[3] = mem[2] + filt_coeff[2] * ( Vu[1] - mem[3] );

        mem[1] = Vu[0];
        mem[2] = Vu[1];
        *out1++ = mem[3];

        /* Middle branch */
        Vm[0] = mem[0] + filt_coeff[3] * (input[i]-mem[4]);
        Vm[1] = mem[4] + filt_coeff[4] * (Vm[0]-mem[5]);
        mem[6] = mem[5] + filt_coeff[5] * (Vm[1]-mem[6]);

        mem[4] = Vm[0];
        mem[5] = Vm[1];
        *out1++ = mem[6];

        /* Lower branch */
        Vl[0] = mem[0] + filt_coeff[6] * (input[i]-mem[7]);
        Vl[1] = mem[7] + filt_coeff[7] * (Vl[0]-mem[8]);
        mem[9] = mem[8] + filt_coeff[8] * (Vl[1]-mem[9]);

        mem[0] = input[i];
        mem[7] = Vl[0];
        mem[8] = Vl[1];
        *out1++ = mem[9];
    }

    /*LPF*/
    for(i = 0; i < len*3; i++)
    {
        mem_temp = out[i];
        out[i] = (((0.572769f)*(mem[12]+mem[11]))-((0.074005f)*(mem_temp+mem[10])));
        mem[10] = mem[11];
        mem[11] = mem[12];
        mem[12] = mem_temp;
    }

    return;
}


/*-------------------------------------------------------------------*
 * retro_interp4_5()
 *
 *
 *-------------------------------------------------------------------*/

void retro_interp4_5(
    const float *syn,
    float *pst_old_syn
)
{
    float *pf5, *pf4;
    short c;

    /* resample st->pst_old_syn in a reverse way to preserve time-alignment */
    pf4 = (float*) &pst_old_syn[58];
    pf5 = (float*) pst_old_syn;
    for (c=0; c<57; c++)
    {
        *pf5++ = pf4[0];
        *pf5++ = 0.2f * pf4[0] + 0.8f * pf4[1];
        *pf5++ = 0.4f * pf4[1] + 0.6f * pf4[2];
        *pf5++ = 0.6f * pf4[2] + 0.4f * pf4[3];
        *pf5++ = 0.8f * pf4[3] + 0.2f * pf4[4];
        pf4+=4;
    }
    *pf5++ = pf4[0];
    *pf5++ = 0.2f * pf4[0] + 0.8f * pf4[1];
    *pf5++ = 0.4f * pf4[1] + 0.6f * pf4[2];
    *pf5++ = 0.6f * pf4[2] + 0.4f * pf4[3];
    *pf5++ = 0.8f * pf4[3] + 0.2f * syn[0];
    /* all samples processed: NBPSF_PIT_MAX = 290 = (58*5) */

    return;
}


/*-------------------------------------------------------------------*
 * retro_interp5_4()
 *
 *
 *-------------------------------------------------------------------*/

void retro_interp5_4( float *pst_old_syn )
{
    float *pf5, *pf4;
    short c;

    /* resample st->pst_old_syn in a reverse way to preserve time-alignment */
    pf4 = (float*) &pst_old_syn[NBPSF_PIT_MAX-1];
    pf5 = pf4;
    for (c=0; c<58; c++)
    {
        *pf4-- = 0.75f * pf5[0] + 0.25f * pf5[-1];
        *pf4-- = 0.50f * pf5[-1] + 0.50f * pf5[-2];
        *pf4-- = 0.25f * pf5[-2] + 0.75f * pf5[-3];
        *pf4-- = pf5[-4];
        pf5-=5;
    }
    /* all samples processed: NBPSF_PIT_MAX = 290 = (58*5) */

    return;
}

