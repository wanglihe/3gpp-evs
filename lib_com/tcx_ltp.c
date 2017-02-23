/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "prot.h"
#include "rom_com.h"
#include "assert.h"


/*-------------------------------------------------------------------
 * Local constants
 *-------------------------------------------------------------------*/

#define ALPHA 0.85f


/*-------------------------------------------------------------------
 * tcx_ltp_get_lpc()
 *
 *
 *-------------------------------------------------------------------*/

static void tcx_ltp_get_lpc(
    float *input,
    int length,
    float *A,
    int lpcorder
)
{
    int i, j;
    float s, r[TCXLTP_LTP_ORDER+1];

    for (i = 0; i <= lpcorder; i++)
    {
        s = 0.0;

        for (j = 0; j < length-i; j++)
        {
            s += input[j]*input[j+i];
        }
        r[i] = s;
    }

    if (r[0] < 100.0f)
    {
        r[0] = 100.0f;
    }
    r[0] *= 1.0001f;

    lev_dur( A, r, lpcorder, NULL );

    return;
}


/*-------------------------------------------------------------------
 * tcx_ltp_get_zir()
 *
 *
 *-------------------------------------------------------------------*/

static void tcx_ltp_get_zir(
    float *zir,
    int length,
    float *synth_ltp,
    float *synth,
    float *A,
    int lpcorder,
    float gain,
    int pitch_int,
    int pitch_fr,
    int pitres,
    int filtIdx
)
{
    float buf[TCXLTP_LTP_ORDER], alpha, step;
    float *x0, *x1, s;
    float *y0, *y1, s2;
    const float *w0, *w1, *v0, *v1;
    int i, j, k, L;

    x0 = &synth_ltp[-pitch_int];
    x1 = x0 - 1;
    y0 = synth;
    y1 = y0 - 1;

    assert(filtIdx >= 0);

    w0 = &tcxLtpFilters[filtIdx].filt[pitch_fr];
    w1 = &tcxLtpFilters[filtIdx].filt[pitres - pitch_fr];
    v0 = &tcxLtpFilters[filtIdx].filt[0];
    v1 = &tcxLtpFilters[filtIdx].filt[pitres];
    L = tcxLtpFilters[filtIdx].length;

    for (j = 0; j < lpcorder; j++)
    {
        s = 0;
        s2 = 0;
        for (i = 0, k = 0; i < L; i++, k += pitres)
        {
            s += w0[k] * x0[i] + w1[k] * x1[-i];
            s2 += v0[k] * y0[i] + v1[k] * y1[-i];
        }
        s2 *= ALPHA;

        buf[j] = ( synth[j] - gain * s2 ) - ( synth_ltp[j] - gain * s );

        x0++;
        x1++;
        y0++;
        y1++;
    }

    set_f( zir, 0.0f, length );

    syn_filt( A, lpcorder, zir, zir, length, buf, 0);

    alpha = 1.f;
    step = 1.f/(float)(length/2);

    for ( j=length/2; j<length; j++ )
    {
        zir[j] *= alpha;
        alpha -= step;
    }

    return;
}


/*-------------------------------------------------------------------
 * predict_signal()
 *
 *
 *-------------------------------------------------------------------*/

void predict_signal(
    const float excI[],  /* i  : input excitation buffer  */
    float excO[],  /* o  : output excitation buffer */
    const short T0,      /* i  : integer pitch lag        */
    short frac,    /* i  : fraction of lag          */
    const short frac_max,/* i  : max fraction             */
    const short L_subfr  /* i  : subframe size            */
)
{
    short j;
    float s;
    const float *x0, *win;

    x0 = &excI[-T0-1];
    frac = -frac;

    if (frac < 0)
    {
        frac += frac_max;
        x0--;
    }

    if ( frac_max == 6 )
    {
        win = &inter6_2tcx2[frac][0];
    }
    else
    {
        win = &inter4_2tcx2[frac][0];
    }

    for (j=0; j<L_subfr; j++)
    {
        s = win[1]*x0[1] + win[2]*x0[2];
        excO[j] = s + win[0]*x0[0] + win[3]*x0[3];
        x0++;
    }

    return;
}


/*-------------------------------------------------------------------
 * tcx_ltp_synth_filter()
 *
 *
 *-------------------------------------------------------------------*/

static void tcx_ltp_synth_filter(
    float *synth_ltp,
    float *synth,
    int length,
    int pitch_int,
    int pitch_fr,
    float gain,
    int pitch_res,
    short filtIdx
)
{
    float *x0, *x1, s;
    float *y0, *y1, s2;
    const float *v0, *v1;
    const float *w0, *w1;
    int i, j, k, L;

    if ( gain > 0.f )
    {
        x0 = &synth_ltp[-pitch_int];
        x1 = x0 - 1;
        y0 = synth;
        y1 = y0 - 1;

        assert(filtIdx >= 0);

        w0 = &tcxLtpFilters[filtIdx].filt[pitch_fr];
        w1 = &tcxLtpFilters[filtIdx].filt[pitch_res - pitch_fr];
        v0 = &tcxLtpFilters[filtIdx].filt[0];
        v1 = &tcxLtpFilters[filtIdx].filt[pitch_res];

        L = tcxLtpFilters[filtIdx].length;


        for (j = 0; j < length; j++)
        {
            s = 0;
            s2 = 0;
            for (i = 0, k = 0; i < L; i++, k += pitch_res)
            {
                s += w0[k] * x0[i] + w1[k] * x1[-i];
                s2 += v0[k] * y0[i] + v1[k] * y1[-i];
            }

            s2 *= ALPHA;
            synth_ltp[j] = synth[j] - gain * s2 + gain * s;

            x0++;
            x1++;
            y0++;
            y1++;
        }
    }
    else
    {
        mvr2r( synth, synth_ltp, length );
    }

    return;
}


/*-------------------------------------------------------------------
 * tcx_ltp_synth_filter_zir()
 *
 *
 *-------------------------------------------------------------------*/

static void tcx_ltp_synth_filter_zir(
    float *synth_ltp,
    float *synth,
    int length,
    int pitch_int,
    int pitch_fr,
    float gain,
    int pitch_res,
    float *zir,
    short filtIdx
)
{
    float *x0, *x1, s;
    float *y0, *y1, s2;
    const float *v0, *v1;
    const float *w0, *w1;
    int i, j, k, L;

    x0 = &synth_ltp[-pitch_int];
    x1 = x0 - 1;
    y0 = synth;
    y1 = y0 - 1;

    assert(filtIdx >= 0);

    w0 = &tcxLtpFilters[filtIdx].filt[pitch_fr];
    w1 = &tcxLtpFilters[filtIdx].filt[pitch_res - pitch_fr];
    v0 = &tcxLtpFilters[filtIdx].filt[0];
    v1 = &tcxLtpFilters[filtIdx].filt[pitch_res];
    L = tcxLtpFilters[filtIdx].length;

    for (j = 0; j < length; j++)
    {
        s = 0;
        s2 = 0;

        for (i = 0, k = 0; i < L; i++, k += pitch_res)
        {
            s += w0[k] * x0[i] + w1[k] * x1[-i];
            s2 += v0[k] * y0[i] + v1[k] * y1[-i];
        }

        s2 *= ALPHA;

        synth_ltp[j] = ( synth[j] - gain * s2 + gain * s ) -zir[j];

        x0++;
        x1++;
        y0++;
        y1++;
    }

    return;
}


/*-------------------------------------------------------------------
 * tcx_ltp_synth_filter_fadein()
 *
 *
 *-------------------------------------------------------------------*/

static void tcx_ltp_synth_filter_fadein(
    float *synth_ltp,
    float *synth,
    int length,
    int pitch_int,
    int pitch_fr,
    float gain,
    int pitch_res,
    short filtIdx
)
{
    float *x0, *x1, s;
    float *y0, *y1, s2;
    const float *v0, *v1;
    const float *w0, *w1;
    int i, j, k, L;
    float alpha, step;

    if ( gain > 0.f )
    {
        x0 = &synth_ltp[-pitch_int];
        x1 = x0 - 1;
        y0 = synth;
        y1 = y0 - 1;

        assert(filtIdx >= 0);

        w0 = &tcxLtpFilters[filtIdx].filt[pitch_fr];
        w1 = &tcxLtpFilters[filtIdx].filt[pitch_res - pitch_fr];
        v0 = &tcxLtpFilters[filtIdx].filt[0];
        v1 = &tcxLtpFilters[filtIdx].filt[pitch_res];
        L = tcxLtpFilters[filtIdx].length;

        alpha = 0.f;
        step = 1.f/(float)(length);


        for (j = 0; j < length; j++)
        {
            s = 0;
            s2 = 0;

            for (i = 0, k = 0; i < L; i++, k += pitch_res)
            {
                s += w0[k] * x0[i] + w1[k] * x1[-i];
                s2 += v0[k] * y0[i] + v1[k] * y1[-i];
            }

            s2 *= ALPHA;

            synth_ltp[j] = synth[j] - alpha * gain * s2 + alpha * gain * s;

            alpha += step;

            x0++;
            x1++;
            y0++;
            y1++;
        }
    }
    else
    {
        mvr2r( synth, synth_ltp, length );
    }

    return;
}


/*-------------------------------------------------------------------
 * tcx_ltp_synth_filter_fadeout()
 *
 *
 *-------------------------------------------------------------------*/

static void tcx_ltp_synth_filter_fadeout(
    float *synth_ltp,
    float *synth,
    int length,
    int pitch_int,
    int pitch_fr,
    float gain,
    int pitch_res,
    short filtIdx
)
{
    float *x0, *x1, s;
    float *y0, *y1, s2;
    const float *v0, *v1;
    const float *w0, *w1;
    int i, j, k, L;
    float alpha, step;

    if ( gain > 0.f )
    {
        x0 = &synth_ltp[-pitch_int];
        x1 = x0 - 1;
        y0 = synth;
        y1 = y0 - 1;

        assert(filtIdx >= 0);

        w0 = &tcxLtpFilters[filtIdx].filt[pitch_fr];
        w1 = &tcxLtpFilters[filtIdx].filt[pitch_res - pitch_fr];
        v0 = &tcxLtpFilters[filtIdx].filt[0];
        v1 = &tcxLtpFilters[filtIdx].filt[pitch_res];
        L = tcxLtpFilters[filtIdx].length;

        alpha = 1.f;
        step = 1.f/(float)(length);

        for (j = 0; j < length; j++)
        {
            s = 0;
            s2 = 0;

            for (i = 0, k = 0; i < L; i++, k += pitch_res)
            {
                s += w0[k] * x0[i] + w1[k] * x1[-i];
                s2 += v0[k] * y0[i] + v1[k] * y1[-i];
            }

            s2 *= ALPHA;
            synth_ltp[j] = synth[j] - alpha * gain * s2 + alpha * gain * s;
            alpha -= step;

            x0++;
            x1++;
            y0++;
            y1++;
        }
    }
    else
    {
        mvr2r( synth, synth_ltp, length );
    }

    return;
}


/*-------------------------------------------------------------------
 * tcx_ltp_decode_params()
 *
 *
 *-------------------------------------------------------------------*/
int tcx_ltp_decode_params(
    int *ltp_param,
    int *pitch_int,
    int *pitch_fr,
    float *gain,
    int pitmin,
    int pitfr1,
    int pitfr2,
    int pitmax,
    int pitres
)
{
    int gainbits = 2;

    /* Decode Pitch and Gain */
    if ((ltp_param) && (ltp_param[0]))
    {
        if ( ltp_param[1] < ((pitfr2-pitmin)*pitres) )
        {
            *pitch_int = pitmin + (ltp_param[1]/pitres);
            *pitch_fr = ltp_param[1] - (*pitch_int - pitmin)*pitres;
        }
        else if (ltp_param[1] < ( (pitfr2-pitmin)*pitres + (pitfr1-pitfr2)*(pitres>>1)) )
        {
            *pitch_int = pitfr2 + ((ltp_param[1]-(pitfr2-pitmin)*pitres)/(pitres>>1));
            *pitch_fr = (ltp_param[1]-(pitfr2-pitmin)*pitres) - (*pitch_int - pitfr2)*(pitres>>1);
            *pitch_fr = *pitch_fr << 1;  /* was *= (pitres>>1); */
        }
        else
        {
            *pitch_int = ltp_param[1] + pitfr1 - ((pitfr2-pitmin)*pitres) - ((pitfr1-pitfr2)*(pitres>>1));
            *pitch_fr = 0;
        }
        *gain = (float)(ltp_param[2] + 1) * 0.625f/(float)(1<<gainbits);
        if(*pitch_int<PIT_MIN_SHORTER)
        {
            /*pitch out of range due to bit error */
            *pitch_int = PIT_MIN_SHORTER;
            return 1;
        }
        if(*pitch_int>PIT_MAX_MAX)
        {
            /*pitch out of range due to bit error */
            *pitch_int = PIT_MAX_MAX;
            return 1;
        }
    }
    else
    {
        *pitch_int = pitmax;
        *pitch_fr = 0;
        *gain = 0.0f;
    }
    return 0;
}


/*-------------------------------------------------------------------
 * tcx_ltp_post()
 *
 *
 *-------------------------------------------------------------------*/

void tcx_ltp_post(
    int tcxltp_on,
    short core,
    int L_frame,
    int L_frame_core,
    int delay,
    float *sig,
    float *tcx_buf,
    short tcx_buf_len,
    int bfi,
    int pitch_int,
    int pitch_fr,
    float gain,
    int *pitch_int_past,
    int *pitch_fr_past,
    float *gain_past,
    int *filtIdx_past,
    int pitres,
    int *pitres_past,
    float damping,
    int SideInfoOnly,
    float *mem_in,
    float *mem_out,
    int bitrate
)
{
    int tmp, L_transition, lpcorder, filtIdx;
    float gain2;
    float zir[L_FRAME_PLUS/4], A[TCXLTP_LTP_ORDER+1];
    float buf_in[TCXLTP_MAX_DELAY+L_FRAME48k+TCXLTP_MAX_DELAY], buf_out[2*L_FRAME48k];
    float *sig_in, *sig_out;


    filtIdx = 0;  /* just to avoid compilation warnings */

    /******** Init ********/

    /* Parameters */
    L_transition = L_frame/4;
    lpcorder = TCXLTP_LTP_ORDER;

    /* Input buffer */
    sig_in = buf_in + tcx_buf_len;
    mvr2r( mem_in, buf_in, tcx_buf_len );
    mvr2r( sig, buf_in+tcx_buf_len, L_frame );
    if ( core > ACELP_CORE )
    {
        mvr2r( tcx_buf, sig_in+L_frame, tcx_buf_len );
    }
    mvr2r( sig+L_frame-tcx_buf_len, mem_in, tcx_buf_len );

    /* Output buffer */
    sig_out = buf_out + L_frame;
    mvr2r( mem_out, buf_out, L_frame );

    /* TCX-LTP parameters: integer pitch, fractional pitch, gain */

    if ( !(SideInfoOnly || tcxltp_on) || core==ACELP_CORE )
    {
        /* No LTP */
        pitch_int = 0;
        pitch_fr = 0;
        gain = 0.f;
    }
    else if ( !bfi )
    {
        /* LTP and good frame */
        if (L_frame != L_frame_core)
        {
            tmp = pitch_int * pitres + pitch_fr;
            tmp = (tmp * L_frame + L_frame_core/2) / L_frame_core;
            pitch_int = tmp / pitres;
            pitch_fr = tmp % pitres;
        }

        if ( bitrate == HQ_48k && L_frame_core == L_FRAME16k )
        {
            gain *= 0.32f;
        }
        else if ( bitrate == HQ_48k && L_frame_core == 512 )
        {
            gain *= 0.40f;
        }
        else
        {
            gain *= 0.64f;
        }

    }
    else
    {
        /* PLC: [TCX: Fade-out]
         * PLC: LTP and bad frame (concealment) */

        pitch_int = *pitch_int_past;
        pitch_fr  = *pitch_fr_past;
        gain = *gain_past * damping;
        pitres = *pitres_past;
    }

    if ( SideInfoOnly )
    {
        gain = 0.f;
        if ( bfi )
        {
            *gain_past = 0.f;
        }
    }
    gain2 = gain;

    if (L_frame_core == L_FRAME)
    {
        switch ( L_frame )
        {
        case L_FRAME8k:
            filtIdx = 0;
            break;
        case L_FRAME16k:
            filtIdx = 1;
            break;
        case L_FRAME32k:
            filtIdx = 2;
            break;
        case L_FRAME48k:
            filtIdx = 3;
            break;
        default:
            assert(0);
            break;
        }
    }
    else if (L_frame_core == L_FRAME16k)
    {
        switch ( L_frame )
        {
        case L_FRAME8k:
            filtIdx = 4;
            break;
        case L_FRAME16k:
            filtIdx = 5;
            break;
        case L_FRAME32k:
            filtIdx = 6;
            break;
        case L_FRAME48k:
            filtIdx = 7;
            break;
        default:
            assert(0);
            break;
        }
    }
    else if (L_frame_core == 512)
    {
        switch ( L_frame )
        {
        case L_FRAME8k:
            filtIdx = 8;
            break;
        case L_FRAME16k:
            filtIdx = 9;
            break;
        case L_FRAME32k:
            filtIdx = 10;
            break;
        case L_FRAME48k:
            filtIdx = 11;
            break;
        default:
            assert(0);
            break;
        }
    }
    else
    {
        filtIdx = -1;
    }


    /******** Previous-frame part ********/
    tcx_ltp_synth_filter( sig_out, sig_in, delay, *pitch_int_past, *pitch_fr_past, *gain_past, *pitres_past, *filtIdx_past );

    /******** Transition part ********/
    if ( gain==0.f && *gain_past==0.f )
    {
        mvr2r( sig_in+delay, sig_out+delay, L_transition );
    }
    else if ( *gain_past==0.f )
    {
        tcx_ltp_synth_filter_fadein(  sig_out+delay, sig_in+delay, L_transition, pitch_int, pitch_fr, gain, pitres, filtIdx );
    }
    else if ( gain==0.f )
    {
        tcx_ltp_synth_filter_fadeout( sig_out+delay, sig_in+delay, L_transition, *pitch_int_past, *pitch_fr_past, *gain_past, *pitres_past, *filtIdx_past );
    }
    else if ( gain==*gain_past && pitch_int==*pitch_int_past && pitch_fr==*pitch_fr_past )
    {
        tcx_ltp_synth_filter( sig_out+delay, sig_in+delay, L_transition, pitch_int, pitch_fr, gain, pitres, filtIdx );
    }
    else
    {
        tcx_ltp_get_lpc( sig_out+delay-L_frame, L_frame, A, lpcorder );

        tcx_ltp_get_zir( zir, L_transition, sig_out+delay-lpcorder, sig_in+delay-lpcorder, A, lpcorder, gain, pitch_int, pitch_fr, pitres, filtIdx );

        tcx_ltp_synth_filter_zir( sig_out+delay, sig_in+delay, L_transition, pitch_int, pitch_fr, gain, pitres, zir, filtIdx );
    }

    /******** Current-frame part ********/
    tcx_ltp_synth_filter( sig_out+(delay+L_transition), sig_in+(delay+L_transition), L_frame-(delay+L_transition),
                          pitch_int, pitch_fr, gain, pitres, filtIdx );

    /******** Output ********/

    /* copy to output */
    mvr2r( sig_out, sig, L_frame );

    /* Update */
    *pitch_int_past = pitch_int;
    *pitch_fr_past = pitch_fr;
    *gain_past = gain2;
    *filtIdx_past = filtIdx;
    *pitres_past = pitres;
    mvr2r( sig_out, mem_out, L_frame );

    return;
}
