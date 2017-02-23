/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include "prot.h"



static void LpFilter2(float *x, float *y, int N, float *mem);
static float harmo(float *X, int n, float f);
static int Is_Periodic(float cov_max, int zp, float ener, float ener_mean,int pitch, int Framesize, float *mdctdata);
static float sig_tilt(float *s, int FrameSize);
static int zero_pass(float *s, int N);
static float dot(float *a, float *b, int N);
static int array_max_indx(float *s, int N);

/* Decoder side */
static void LpFilter2(float *x, float *y, int N, float *mem)
{
    int i;
    y[0] = 0.18f * mem[0] + 0.64f * mem[1] + 0.18f * x[0];
    y[1] = 0.18f * mem[1] + 0.64f * y[0] + 0.18f * x[1];
    for (i = 2; i < N; i++)
    {
        y[i] = 0.18f * y[i-2] + 0.64f * y[i-1] + 0.18f * x[i];
    }
}

static float harmo(float *X, int n, float f)
{
    int h, k, m = 8;
    float ener = 0, ener_harmo = 0;
    for (k = 1; k < m+1; k++)
    {
        h = (int)(k*f - 0.5f);
        if (k*f - h > 0.5f)
        {
            ener_harmo += X[h]*X[h] + X[h+1]*X[h+1];
        }
        else
        {
            ener_harmo += X[h]*X[h];
        }
    }
    for (k = 0; k < n; k++)
    {
        ener += X[k]*X[k];
    }
    return ener_harmo/(ener+EPSILON);
}

static int Is_Periodic(float cov_max, int zp, float ener, float ener_mean,int pitch, int Framesize, float *mdctdata)
{
    int flag =0;
    float f = 2.0f*Framesize/pitch;
    float harm;
    harm = harmo(mdctdata /*X*/, Framesize, f);
    if (ener < 50 || (ener < ener_mean - 8.0f && cov_max < 0.9f))
    {
        flag = 0;
    }
    else if (cov_max > 0.8f)
    {
        flag = 1;
    }
    else if (zp > 100)
    {
        flag = 0;
    }
    else if (ener < ener_mean - 6)
    {
        flag = 0;
    }
    else if (ener > ener_mean + 1 && cov_max > 0.6f)
    {
        flag = 1;
    }
    else if (harm < 0.7f)
    {
        flag = 0;
    }
    else
    {
        flag = 1;
    }
    return flag;
}

static int zero_pass(float *s, int N)
{
    int zp = 0, i;
    for (i = 1; i < N; i++)
    {
        if (s[i-1] * s[i] <= 0)
        {
            zp++;
        }
    }
    return zp;
}

static float dot(float *a, float *b, int N)
{
    float sum = 0;
    int i;
    for (i = 0; i < N; i++)
    {
        sum += (float)(a[i] * b[i]);
    }
    return sum;
}

static int array_max_indx(float *s, int N)
{
    int i, indx = 0;
    for (i = 0; i < N; i++)
    {
        if (s[i] > s[indx])
        {
            indx = i;
        }
    }
    return indx;
}

static float sig_tilt(float *s, int FrameSize)
{
    float tilt, enr1, enr2;
    int subFrameSize, shift = 2;
    float *p1, *p2;
    subFrameSize = FrameSize>>2;
    p1 = s+subFrameSize;
    p2 = s+subFrameSize-shift;
    enr1 = dot(p1, p2, FrameSize-subFrameSize);
    enr2 = dot(p1, p1, FrameSize-subFrameSize);
    tilt = enr1 / (enr2+EPSILON);
    return tilt;
}

static int pitch_search(float *s, /* lastPcmOut */
                        float *outx_new,
                        int   Framesize,
                        float *voicing,
                        int *zp,
                        float *ener,
                        float ener_mean,
                        float *mdct_data,
                        int curr_mode
                       )
{
    int   pitch = 0, t, i;
    float cov_max = 0, temp = 0, tmp, tilt, mdct_ener = 0, low_freq_rate;
    float s_LP[L_FRAME_MAX] = {0};
    float mem[2]= {0};
    short start_pos, end_pos;
    int  cov_size;
    int flag = 0, zp_current;
    int curr_frmsize = Framesize/curr_mode;
    float tmp_last    = 0; /* not needed, just to avoid compiler warning */
    float cov_max_tmp = 0;
    zp_current = zero_pass(outx_new, curr_frmsize);

    if (curr_mode==2)
    {
        zp_current = zp_current << 1;
    }
    if (Framesize <= 256)
    {
        if (zp_current > 70)
        {
            return 0;
        }
    }
    else
    {
        if (zp_current > 105)
        {
            return 0;
        }
    }
    if (curr_mode == 2)
    {
        mdct_data = mdct_data + (Framesize>>1);
    }

    for (i = 0; i < 30/curr_mode; i++)
    {
        mdct_ener += mdct_data[i]*mdct_data[i];
    }
    low_freq_rate = mdct_ener;
    for (i = 30/curr_mode; i < Framesize/curr_mode; i++)
    {
        mdct_ener += mdct_data[i]*mdct_data[i];
    }
    low_freq_rate /= (mdct_ener+EPSILON);
    if (curr_mode == 2)
    {
        mdct_data = mdct_data - (Framesize>>1);
    }
    if(low_freq_rate < 0.02f)
    {
        return 0;
    }

    LpFilter2(s, s_LP, Framesize, mem);
    tilt = sig_tilt(s_LP, Framesize);
    if (Framesize <= 320)
    {
        if (tilt < 0.5f)
        {
            return 0;
        }
    }
    else
    {
        if (tilt < 0.7f)
        {
            return 0;
        }
    }
    if (Framesize <= 320)
    {
        start_pos = (int)(Framesize *  34/256.0 + 0.5f);
        end_pos   = (int)(Framesize *     3/4.0 + 0.5f);
        for (t = start_pos; t < end_pos; t++)
        {
            cov_size =  Framesize-t;
            tmp = dot(s_LP, s_LP+t, cov_size)/cov_size;
            if (t > start_pos)                  /* don't use the first value */
            {
                if ( tmp > tmp_last)              /* store the current cov, if it is larger than the last one */
                {
                    cov_max_tmp = tmp;
                }
                else if (cov_max < cov_max_tmp)   /* otherwise */
                {
                    cov_max = cov_max_tmp;          /* use the last value cov, being a max */
                    pitch = t-1;                    /* and the last index as pitch */
                }
            }
            tmp_last = tmp;
        }
        temp = (float)(sqrt(dot(s_LP+pitch, s_LP+pitch, Framesize-pitch)) *
                       sqrt(dot(s_LP,s_LP,Framesize-pitch)));
        *voicing = cov_max * (Framesize - pitch) / (temp + EPSILON);
        {
            float temp2, voicing2;
            temp2 = (float)(sqrt(dot(s_LP+(pitch>>1), s_LP+(pitch>>1), Framesize-(pitch>>1))) *
                            sqrt(dot(s_LP,s_LP,Framesize-(pitch>>1))));
            voicing2 = dot(s_LP+(pitch>>1),s_LP,Framesize-(pitch>>1))/temp2;
            if (voicing2 > *voicing)
            {
                pitch = pitch>>1;
                *voicing = voicing2;
            }
        }
    }
    else
    {
        float s_tmp[L_FRAME_MAX];
        int Framesize_tmp = Framesize>>1;
        int i, pitch_tmp[3];
        for (i = 0; i < Framesize_tmp; i++)
        {
            s_tmp[i] = s_LP[2*i];
        }
        start_pos = (int)((34.0f*Framesize_tmp)/256 + 0.5f);
        end_pos  = (int)((Framesize_tmp>>1)*1.5f + 0.5f);
        cov_max = 0;
        pitch = 0;
        for (t = start_pos; t < end_pos; t++)
        {
            cov_size =  Framesize_tmp-t;
            tmp = dot(s_tmp, s_tmp+t, cov_size)/cov_size;
            if (t > start_pos)                  /* don't use the first value */
            {
                if ( tmp > tmp_last)              /* store the current cov, if it is larger than the last one */
                {
                    cov_max_tmp = tmp;
                }
                else if (cov_max < cov_max_tmp)   /* otherwise */
                {
                    cov_max = cov_max_tmp;          /* use the last value cov, being a max */
                    pitch = t-1;                    /* and the last index as pitch */
                }
            }
            tmp_last = tmp;
        }

        if (pitch > 0)
        {
            pitch_tmp[0] = max(2*pitch - 1,0);
            pitch_tmp[1] = 2*pitch;
            pitch_tmp[2] = 2*pitch + 1;
            cov_max = 0;
            pitch = 0;
            for (i = 0; i < 3; i++)
            {
                cov_size =  Framesize - pitch_tmp[i];
                temp = (float)(sqrt(dot(s_LP+pitch_tmp[i], s_LP+pitch_tmp[i], cov_size)) *
                               sqrt(dot(s_LP,s_LP,cov_size)));
                tmp = dot(s_LP, s_LP + pitch_tmp[i], cov_size)/(temp+EPSILON);
                if (tmp > cov_max)
                {
                    cov_max = tmp;
                    pitch = pitch_tmp[i];
                }
            }
            *voicing = cov_max;
        }
    }
    if (pitch > 0)
    {
        flag = Is_Periodic(*voicing, *zp, *ener, ener_mean, pitch, Framesize, mdct_data);
    }
    if (flag == 0 )
    {
        pitch = 0;
    }
    return pitch;
}

static  int OverlapAdd(float *pitch125_data, float *sbuf, int n, int pitch, int Bufsize)
{
    int pitch125 = (int)floor(0.5f+(1.25f*(float)pitch));
    int Loverlap = pitch125 - pitch;
    int n1 = min(Loverlap, Bufsize-n);
    int n2 = min(pitch125, Bufsize-n);
    int i;
    float tmp, dat;
    for (i = 0; i < n1; i++)
    {
        tmp = (float)i/(float)Loverlap;
        dat = sbuf[n+i];
        sbuf[n+i] = (float)((1.0 - tmp)*dat + tmp*pitch125_data[i]);
    }
    for (i = n1; i < n2; i++)
    {
        sbuf[n+i] = pitch125_data[i];
    }
    return (n+pitch);
}

static void add_noise (float      * const sbuf,
                       float      * const outx_new_n1,
                       float const* const noise_seg,
                       int          const Len,
                       float      * const gain,
                       float const* const gain_n,
                       int          const firstFrame)
{
    int i;

    if (!firstFrame)
    {
        sbuf[0] += *gain * (noise_seg[0] - 0.68f*(*outx_new_n1));
        *gain = 0.99f*(*gain) + 0.01f*(*gain_n);
    }
    for (i = 1; i < Len; i++)
    {
        sbuf[i] += *gain * (noise_seg[i] - 0.68f*noise_seg[i-1]);
        *gain = 0.99f*(*gain) + 0.01f*(*gain_n);
    }
    *outx_new_n1 = noise_seg[i-1];

}

static
int waveform_adj(
    float *overlapbuf,
    float *outdata2,
    float *outx_new,
    float *data_noise,
    float *outx_new_n1,
    float *nsapp_gain,
    float *nsapp_gain_n,
    int    Framesize,
    int    T_bfi,
    float  voicing,
    int    curr_mode,
    int    pitch)
{
    int i, zp1, zp2;
    float sbuf[L_FRAME_MAX] = {0};
    zp1 = zero_pass(outdata2, Framesize>>1);
    zp2 = zero_pass(outdata2+(Framesize>>1), Framesize>>1);

    /* judge if the pitch is usable */
    if( 4*max(zp1,1) < zp2 )
    {
        return 0;
    }

    /* adjust the pitch value */
    if (T_bfi && pitch <= Framesize>>1 && Framesize > 256 && curr_mode == 1)
    {
        int i1 = 0, i2 = 0;
        i1 = 1 + array_max_indx(outx_new, pitch);
        i2 = 1 + array_max_indx(outx_new+pitch, pitch);
        if ((float)(i2+pitch-i1)<(1.25f*pitch) && (1.25f*(i2+pitch-i1))>(float)pitch && (float)(i2+pitch-i1)<(float)(Framesize>>1))
        {
            pitch = i2+pitch-i1;
        }
    }

    {
        int pitch125 = 0, Loverlap = 0, n = 0, i;
        float pitch125_data[L_FRAME_MAX] = {0};

        pitch125 = (int)floor(0.5f+(1.25f*(float)pitch));
        Loverlap = pitch125 - pitch;
        memcpy(pitch125_data, outdata2+Framesize-pitch, sizeof(float)*pitch);
        memcpy(pitch125_data+pitch, outx_new, sizeof(float)*Loverlap);
        memcpy(sbuf, outx_new, sizeof(float)*Framesize);
        {
            float tmp[L_FRAME_MAX]= {0}, *p_tmp = tmp+1;
            memcpy(p_tmp, pitch125_data, sizeof(float)*pitch125);
            p_tmp[-1] = outdata2[Framesize-pitch-1];
            p_tmp[pitch125] = outx_new[Loverlap];
            for (i = 0; i < pitch125; i++)
            {
                pitch125_data[i] = 0.18f*p_tmp[i-1]+0.64f*p_tmp[i]+0.18f*p_tmp[i+1];
            }
        }
        while (n < Framesize)   /* periodical extension */
        {
            n = OverlapAdd(pitch125_data,sbuf,n,pitch,Framesize);
        }
        /* maximum pitch lag is 3/4 Framesize; pitch125_data is reused for
           temporary storage, since outdata2 (holding the pcm data of the
           last good frame) is still needed and overlapbuf overlaps outdata2 */
        mvr2r(&sbuf[Framesize/4], pitch125_data, (3*Framesize)/4);

        *nsapp_gain   = 0.0;
        *nsapp_gain_n = 1.0f - voicing/2;
        {
            int size = Framesize;
            /* use last good signal for noise generation */
            add_noise(sbuf, outx_new_n1, outdata2, size, nsapp_gain, nsapp_gain_n, 1);
            /* save current (noisy) output from IMDCT */
            mvr2r(outx_new, data_noise, size);
        }
        /* overlapbuf can now be filled with sbuf, needed for subsequently lost frames */
        mvr2r(pitch125_data, &overlapbuf[Framesize/4], (3*Framesize)/4);
    }
    for (i = 0; i < Framesize; i++)
    {
        outx_new[i] = sbuf[i];
    }
    return pitch;
}

void waveform_adj2(float *overlapbuf,
                   float *outx_new,
                   float *data_noise,
                   float *outx_new_n1,
                   float *nsapp_gain,
                   float *nsapp_gain_n,
                   float *recovery_gain,
                   float step_concealgain,
                   int    pitch,
                   int    Framesize,
                   int    delay,
                   int    bfi_cnt,
                   int    bfi)
{
    int i, n=0;
    float ratio;
    float sbuf[L_FRAME_MAX];
    if (pitch > 0)
    {
        while (n < Framesize)
        {
            /* periodical extension */
            for (i = 0; i < min(pitch, Framesize-n); i++)
            {
                sbuf[n+i] = overlapbuf[Framesize-pitch+i];
            }
            n += pitch;
        }
        for (i = 0; i < Framesize; i++)
        {
            overlapbuf[i] = sbuf[i];
        }
        {
            int    size      = Framesize;
            float* noise_ptr = data_noise;
            /* use last (noisy) output from IMDCT for noise generation */
            add_noise(sbuf, outx_new_n1, noise_ptr, size, nsapp_gain, nsapp_gain_n, 0);
            if (bfi)
            {
                /* save current (noisy) output from IMDCT */
                mvr2r(outx_new, noise_ptr, size);
            }
        }
        if (bfi_cnt == 4 || bfi == 0)
        {
            if (bfi == 0)
            {
                int gain_zero_start = 10000;
                /* overlap-and-add */
                if (step_concealgain > EPSILON)
                {
                    gain_zero_start = (int)min((float)L_FRAME48k, (*recovery_gain/step_concealgain)) + 1;
                }

                if (delay > 0)
                {
                    Framesize -= delay;
                }
                for (i = 0; i < min(gain_zero_start, Framesize); i++)
                {
                    ratio = (float)i/(float)Framesize;
                    outx_new[i] = (1-ratio)*sbuf[i]**recovery_gain+ratio*outx_new[i];
                    *recovery_gain -= step_concealgain;
                }
                for (i = gain_zero_start; i < Framesize; i++)
                {
                    ratio = (float)i/(float)Framesize;
                    outx_new[i] = ratio*outx_new[i];
                }
                if (*recovery_gain < 0.0f)
                {
                    *recovery_gain = 0.0f;
                }
            }
            else
            {
                /* overlap-and-add */
                for (i = 0; i < Framesize; i++)
                {
                    ratio = (float)i/(float)Framesize;
                    outx_new[i] = (1-ratio)*sbuf[i]+ratio*outx_new[i];
                }
            }
        }
        else
        {
            for (i = 0; i < Framesize; i++)
            {
                outx_new[i] = sbuf[i];
            }
        }
    }
    return;
}

void set_state(int *state, int num, int N)
{
    int i;
    for (i = 0; i < N-1; i++)
    {
        state[i] = state[i+1];
    }
    state[N-1] = num;
    return;
}

void concealment_init(int N, void *_plcInfo)
{
    T_PLCInfo *plcInfo = (T_PLCInfo*)_plcInfo;
    int i;

    plcInfo->FrameSize      = N;
    plcInfo->Pitch          = 0;
    plcInfo->T_bfi          = 0;
    plcInfo->outx_new_n1    = 0.0f;
    plcInfo->nsapp_gain     = 0.0f;
    plcInfo->nsapp_gain_n   = 0.0f;
    plcInfo->ener_mean      = 59.4260f;
    plcInfo->ener           = 0.0f;
    plcInfo->zp             = N;
    plcInfo->recovery_gain  = 0.0f;
    plcInfo->step_concealgain = 0.0f;
    plcInfo->concealment_method = TCX_NONTONAL;
    plcInfo->subframe       = 0;
    plcInfo->nbLostCmpt     = 0;
    plcInfo->seed = 21845;

    for (i = 0; i < TCX_TONALITY_INIT_CNT; i++)
    {
        plcInfo->TCX_Tonality[i] = 1;
    }
    for (i = TCX_TONALITY_INIT_CNT; i < DEC_STATE_LEN; i++)
    {
        plcInfo->TCX_Tonality[i] = 0;
    }
    for (i = 0; i < MAX_POST_LEN; i++)
    {
        plcInfo->Transient[i] = 1;
    }

    for (i = 0; i < L_FRAME_MAX; i++)
    {
        plcInfo->data_reci2[i] = 0;
    }
    return;
}

void concealment_decode(int curr_mode, float *invkoef, void *_plcInfo)
{
    T_PLCInfo *plcInfo = (T_PLCInfo*)_plcInfo;
    int i;
    int N = plcInfo->FrameSize;
    short *seed = &(plcInfo->seed);
    short sign;
    if (plcInfo->concealment_method == TCX_NONTONAL)
    {
        if (curr_mode == 1)
        {
            /* copy the data of the last frame */
            mvr2r(plcInfo->data_reci2, invkoef, N);
            /* sign randomization */
            for (i = 0; i < N; i++)
            {
                sign = ((own_random(seed)>>15)<<1)+1;
                invkoef[i] *= sign;
            }
        }
    }
    return;
}

void concealment_update(int bfi, int curr_mode, int tonality, float *invkoef, void *_plcInfo)
{
    T_PLCInfo *plcInfo  = (T_PLCInfo*)_plcInfo;
    float *data_reci2   = plcInfo->data_reci2;
    int   *tcx_tonality = plcInfo->TCX_Tonality;
    int    FrameSize    = plcInfo->FrameSize;
    int    subframe     = plcInfo->subframe;
    int    i;
    if (curr_mode == 1)
    {
        set_state(plcInfo->Transient, curr_mode, MAX_POST_LEN);
        for (i = 0; i < FrameSize; i++)
        {
            data_reci2[i] = invkoef[i];
        }
        if (!bfi)
        {
            set_state(tcx_tonality, tonality, DEC_STATE_LEN);
        }
    }
    else
    {
        if (subframe == 0)
        {
            set_state(plcInfo->Transient, curr_mode, MAX_POST_LEN);
            if (!bfi)
            {
                set_state(tcx_tonality, tonality, DEC_STATE_LEN);
            }
        }
        /* don't store the second subframe during frameloss; in
           pitch_search(), low_freq_rate is derived on the last good
           TCX-10 spectrum */
        if (!bfi || subframe == 0)
        {
            float *ptr = data_reci2+subframe;
            for (i = 0; i < FrameSize>>1; i++)
            {
                ptr[i] = invkoef[i];
            }
        }
    }
    return;
}

void concealment_update2(float *outx_new, void *_plcInfo, int FrameSize)
{
    T_PLCInfo *plcInfo = (T_PLCInfo*)_plcInfo;
    plcInfo->zp = zero_pass(outx_new, FrameSize);
    plcInfo->ener = dot(outx_new, outx_new, FrameSize)/FrameSize;
    plcInfo->ener = 10 * (float)log10(plcInfo->ener+EPSILON);
    if (plcInfo->zp < 100 && plcInfo->ener > 50)
    {
        plcInfo->ener_mean = 0.98f*plcInfo->ener_mean + 0.02f*plcInfo->ener;
    }
    return;
}

void concealment_signal_tuning(int bfi, int curr_mode, float *outx_new, void *_plcInfo, int nbLostCmpt, int pre_bfi,
                               float *OverlapBuf, int past_core_mode, float *outdata2, Decoder_State *st)
{
    T_PLCInfo *plcInfo = (T_PLCInfo*)_plcInfo;
    int    FrameSize   = plcInfo->FrameSize;
    float *data_reci2  = plcInfo->data_reci2;
    int    Pitch       = plcInfo->Pitch;
    float voicing      = 0;
    if (bfi)
    {

        if (st->enablePlcWaveadjust && plcInfo->concealment_method == TCX_NONTONAL)
        {
            if (nbLostCmpt == 1)
            {
                plcInfo->Pitch = pitch_search(outdata2,
                                              outx_new,
                                              FrameSize,
                                              &voicing,
                                              &plcInfo->zp,
                                              &plcInfo->ener,
                                              plcInfo->ener_mean,
                                              data_reci2,
                                              curr_mode
                                             );
                if (plcInfo->Pitch)    /* waveform adjustment for the first lost frame */
                {
                    plcInfo->Pitch = waveform_adj(OverlapBuf,
                                                  outdata2,
                                                  outx_new,
                                                  plcInfo->data_noise,
                                                  &plcInfo->outx_new_n1,
                                                  &plcInfo->nsapp_gain,
                                                  &plcInfo->nsapp_gain_n,
                                                  FrameSize,
                                                  plcInfo->T_bfi,
                                                  voicing,
                                                  curr_mode,
                                                  plcInfo->Pitch);
                }
            }
            else if (nbLostCmpt < 5)      /* waveform adjustment for the 2nd~4th lost frame */
            {
                waveform_adj2(OverlapBuf,
                              outx_new,
                              plcInfo->data_noise,
                              &plcInfo->outx_new_n1,
                              &plcInfo->nsapp_gain,
                              &plcInfo->nsapp_gain_n,
                              &plcInfo->recovery_gain,
                              plcInfo->step_concealgain,
                              Pitch,
                              FrameSize,
                              0,
                              nbLostCmpt,
                              bfi);
            }
        }
        plcInfo->T_bfi = 1;
    }
    else
    {
        if (pre_bfi &&
                past_core_mode != 0 &&
                st->last_total_brate >= HQ_48k &&
                st->last_codec_mode == MODE2)
        {
            if (plcInfo->concealment_method == TCX_NONTONAL)
            {
                if (plcInfo->nbLostCmpt < 4)   /* smoothing of the concealed signal with the good signal */
                {
                    waveform_adj2(OverlapBuf,
                                  outx_new,
                                  plcInfo->data_noise,
                                  &plcInfo->outx_new_n1,
                                  &plcInfo->nsapp_gain,
                                  &plcInfo->nsapp_gain_n,
                                  &plcInfo->recovery_gain,
                                  plcInfo->step_concealgain,
                                  Pitch,
                                  FrameSize,
                                  0,
                                  plcInfo->nbLostCmpt + 1,
                                  bfi);
                }
            }
        }
        else
        {
            plcInfo->T_bfi = 0;
        }
    }
    return;
}
