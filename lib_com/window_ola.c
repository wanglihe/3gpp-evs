/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"

/*--------------------------------------------------------------------------
* window_ola()
*
* Windowing, Overlap and Add
*--------------------------------------------------------------------------*/


void window_ola(
    const float *ImdctOut,          /* i  : input                       */
    float *auOut,             /* o  : output audio                */
    float *OldauOut,          /* i/o: audio from previous frame   */
    const short L,                  /* i  : length                      */
    const short right_mode,
    const short left_mode,          /* window overlap of current frame (0: full, 2: none, or 3: half) */
    const short use_bfi_win,        /* i  : use BFI windowing           */
    const short oldHqVoicing,       /* i  : previous HqVoicing          */
    float *oldgapsynth        /* i  : previous gapsynth           */
)
{
    short i, decimate, decay;
    short n ,n16, windecay48, windecay16;
    float win_right[R2_48];
    float win_int_left[R1_16];
    float win_left[R1_48];
    float win_int_right[R2_16];


    float SS2[L_FRAME48k-NS2SA(48000, N_ZERO_MDCT_NS)];
    float wret2[L_FRAME48k-NS2SA(48000, N_ZERO_MDCT_NS)];

    float *paout;


    n = (short)((float)L*N_ZERO_MDCT_NS/FRAME_SIZE_NS);
    n16 = (short)((float)L_FRAME16k*N_ZERO_MDCT_NS/FRAME_SIZE_NS);
    windecay48 = (short)(2*((float)L_FRAME48k*N_ZERO_MDCT_NS/FRAME_SIZE_NS))+R1_48;
    windecay16 = (short)(2*((float)L_FRAME16k*N_ZERO_MDCT_NS/FRAME_SIZE_NS))+R1_16;
    decimate = 1;
    decay = 0;

    tcx_get_windows_mode1(   left_mode,  right_mode, win_left, win_right, win_int_left, win_int_right, L);

    if ( L == L_FRAME32k || L == L_FRAME16k )
    {
        decimate = 3;
        decay = 1;
    }
    else if ( L == L_FRAME8k )
    {
        decimate = 6;
        decay = 2;
    }

    else if (L == 512)
    {
        windecay48 = (short)(2*((float)512*N_ZERO_MDCT_NS/FRAME_SIZE_NS))+R1_25;
        decimate = 1;
        decay = 0;
    }
    else if (L== 256)
    {
        windecay48 = (short)(2*((float)512*N_ZERO_MDCT_NS/FRAME_SIZE_NS))+R1_25;
        decimate = 2;
        decay = 0;
    }

    paout=auOut-n;

    if( use_bfi_win )
    {
        if ( L == L_FRAME32k )
        {
            for (i = 0; i < L/2; i+=2)
            {
                wret2[L/2-n + i +1] = win_left[(L_FRAME16k/2-i/2-1)*decimate+decay];
                wret2[L/2-n + i] =  win_int_left[(L_FRAME16k/2-i/2-1)];
            }

            for (i = n; i < L/2; i+=2)
            {
                wret2[ i -n] = win_left[(L_FRAME16k-i/2)*decimate-decay-1];
                wret2[ i +1-n] = win_int_left[L_FRAME16k-(i/2)-1];
            }
        }
        else
        {
            for (i = 0; i < L/2; i++)
            {
                wret2[i+L/2-n] =  win_left[(L/2-i-1)*decimate+decay];
            }

            for (i =n; i < L/2; i++)
            {
                wret2[i-n] =  win_left[(L-i)*decimate-decay-1];
            }
        }
        sinq(EVS_PI/(2*(L-n)), EVS_PI/(4*(L-n)), L-n, SS2);

        for (i = n; i < L/2; i++)
        {
            paout[i] = ImdctOut[L/2 + i];
        }

        for (i = 0; i < L/2; i++)
        {
            paout[L/2 + i] = -ImdctOut[L - 1 - i];
        }

        if( oldHqVoicing )
        {
            for (i=0 ; i < L-n; i++)
            {
                auOut[i] = auOut[i]*SS2[i]+ oldgapsynth[i+n]*(SS2[L-n-i-1]);
            }
        }
        else
        {
            for (i=0 ; i < L-n; i++)
            {
                auOut[i] =auOut[i]*SS2[i]+ OldauOut[i+n]*(SS2[L-n-i-1])/(wret2[i]+0.01f);
            }
        }
    }

    if ( L == L_FRAME32k )
    {
        if (use_bfi_win==0)
        {
            for (i = n; i < L/2; i+=2)
            {
                paout[i] = ImdctOut[L/2 + i] * win_right[(2*L_FRAME16k-(n16+(i-n)/2))*decimate-1-decay-windecay48]+OldauOut[i];
                paout[i+1] = ImdctOut[L/2 + i +1] * win_int_right[2*L_FRAME16k-(n16+(i-n)/2)-1-windecay16]+OldauOut[i+1];
            }

            for (i = 0; i < L/2-n; i+=2)
            {
                paout[L/2 + i +1] = -ImdctOut[L - 1 - (i+1)] * win_right[(3*L_FRAME16k/2-1-i/2)*decimate+decay-windecay48]+OldauOut[i+L/2+1];
                paout[L/2 + i ] = -ImdctOut[L - 1 - i] * win_int_right[(3*L_FRAME16k/2-1-i/2)-windecay16]+OldauOut[i+L/2];
            }

            for (i = L/2-n; i < L/2; i+=2)
            {
                paout[L/2 + i +1] = -ImdctOut[L - 1 - (i+1)]+OldauOut[i+L/2+1] ;
                paout[L/2 + i ] = -ImdctOut[L - 1 - i]+OldauOut[i+L/2];
            }
        }

        for (i = 0; i < L/2; i+=2)
        {
            OldauOut[L/2 + i +1] = -ImdctOut[i+1] * win_left[(L_FRAME16k/2-i/2-1)*decimate+decay];
            OldauOut[L/2 + i] = -ImdctOut[i] * win_int_left[(L_FRAME16k/2-i/2-1)];
        }


        for (i = n; i < L/2; i+=2)
        {
            OldauOut[ i] = -ImdctOut[L/2 - 1 - i]  *win_left[(L_FRAME16k-i/2)*decimate-decay-1];
            OldauOut[ i +1] = -ImdctOut[L/2 - 1 - (i +1)] * win_int_left[L_FRAME16k-(i/2)-1];
        }
    }
    else
    {
        if (use_bfi_win==0)
        {

            for (i = n; i < L/2; i++)
            {
                paout[i] = ImdctOut[L/2 + i] * win_right[(2*L-i)*decimate-1-decay-windecay48]+OldauOut[i];
            }

            for (i = 0; i < L/2-n; i++)
            {
                paout[L/2 + i] = -ImdctOut[L - 1 - i] * win_right[(3*L/2-1-i)*decimate+decay-windecay48]+OldauOut[i+L/2];
            }
            for (i = L/2-n; i < L/2; i++)
            {
                paout[L/2 + i] = -ImdctOut[L - 1 - i] + OldauOut[i+L/2];
            }
        }

        for (i = 0; i < L/2; i++)
        {
            OldauOut[L/2 + i] = -ImdctOut[i] * win_left[(L/2-i-1)*decimate+decay];
        }


        for (i = n; i < L/2; i++)
        {
            OldauOut[ i] = -ImdctOut[L/2 - 1 - i] * win_left[(L-i)*decimate-decay-1];
        }


    }

    for (i = 0; i < n; i++)
    {
        OldauOut[i] = -ImdctOut[L/2 - 1 - i];
    }
    for (i = 0; i < n; i++)
    {
        paout[L + i] = OldauOut[i];
    }

    return;
}

/*---------------------------------------------------------------------*
* core_switching_OLA()
*
* modify window after HQ core decoding
* Overlap ACELP and HQ
*---------------------------------------------------------------------*/

void core_switching_OLA(
    const float *mem_over_hp,       /* i  : upsampling filter memory            */
    const short last_L_frame,       /* i  : last L_frame lengthture             */
    const int   output_Fs,          /* i  : output sampling rate                */
    float *synth,             /* i/o: synthesized signal from HQ core     */
    const float *synth_subfr_out,   /* i  : synthesized signal from ACELP core  */
    float *synth_subfr_bwe,   /* i  : synthesized BWE from ACELP core     */
    const short output_frame,       /* i  : output frame length                 */
    const short bwidth              /* i  : output bandwidth                    */
)
{
    short i, L, Loverlapp, out_filt_length, filt_delay, decimate, decay;
    float tmp_buf_switch[SWITCH_MAX_GAP], tmp_buf_switch2[HQ_DELAY_COMP*HQ_DELTA_MAX+2];
    float delta;
    const float *win, *win_int;

    win = window_48kHz;
    win_int = window_8_16_32kHz;
    decimate = 1;
    decay = 0;

    if( output_frame == L_FRAME32k || output_frame == L_FRAME16k )
    {
        decimate = 3;
        decay = 1;
    }
    else if( output_frame == L_FRAME8k )
    {
        decimate = 6;
        decay = 2;
    }

    /* set multiplication factor according to the sampling rate */
    delta = 1;
    if( output_frame == L_FRAME16k )
    {
        delta = 2;
    }
    else if( output_frame == L_FRAME32k )
    {
        delta = 4;
    }
    else if( output_frame == L_FRAME48k )
    {
        delta = 6;
    }

    set_f( tmp_buf_switch, 0, SWITCH_MAX_GAP );
    set_f( tmp_buf_switch2, 0, HQ_DELAY_COMP*HQ_DELTA_MAX+2);

    Loverlapp = (short)(delta*SWITCH_OVERLAP_8k);

    mvr2r( synth_subfr_out,tmp_buf_switch, NS2SA(output_Fs, SWITCH_GAP_LENGTH_NS));  /* copy 6.25 ms subframe  */

    /* conversion from 12.8kHz to output_Fs */
    if( last_L_frame == L_FRAME)
    {
        /* resample filter memory */
        if( output_frame == L_FRAME8k )
        {
            mvr2r( synth_subfr_out + NS2SA(output_Fs, SWITCH_GAP_LENGTH_NS), tmp_buf_switch+NS2SA(output_Fs,SWITCH_GAP_LENGTH_NS),NS2SA(output_Fs, DELAY_CLDFB_NS)); /* copy subframe to tmp buffer */
        }
        else
        {
            out_filt_length = modify_Fs_intcub3m_sup( mem_over_hp+2, NS2SA(12800,DELAY_CLDFB_NS), 12800, tmp_buf_switch2, output_Fs , &filt_delay );

            for( i=0; i<filt_delay; i++ )
            {
                tmp_buf_switch2[out_filt_length+i-filt_delay] = tmp_buf_switch2[out_filt_length-1-filt_delay];
            }

            if ( (bwidth == NB && output_Fs>=16000) || (bwidth > NB && output_Fs>16000) )
            {
                /* mix cubic and CLDFB resampled buffers in case of resampling to higher frequency rates */
                for(i=0; i<out_filt_length; i++)
                {
                    float a = (float) i / out_filt_length;
                    float b = 1.f - a;

                    tmp_buf_switch[NS2SA(output_Fs, SWITCH_GAP_LENGTH_NS) + i] = a * tmp_buf_switch2[i] + b * synth_subfr_out[NS2SA(output_Fs, SWITCH_GAP_LENGTH_NS) + i];
                }
            }
            else
            {
                /* copy cubic resampled buffer (memory) */
                mvr2r( tmp_buf_switch2, tmp_buf_switch+NS2SA(output_Fs,SWITCH_GAP_LENGTH_NS), out_filt_length );
            }
        }
    }
    else
    {
        if( output_frame == L_FRAME16k )
        {
            /* no resampling */
            mvr2r( mem_over_hp+2, tmp_buf_switch+NS2SA(output_Fs,SWITCH_GAP_LENGTH_NS), NS2SA(output_Fs,DELAY_CLDFB_NS) );
        }
        else
        {
            if( output_frame == L_FRAME8k )
            {
                mvr2r( synth_subfr_out+NS2SA(output_Fs,SWITCH_GAP_LENGTH_NS), tmp_buf_switch+NS2SA(output_Fs,SWITCH_GAP_LENGTH_NS), NS2SA(output_Fs, DELAY_CLDFB_NS)); /* copy subframe to tmp buffer */
            }
            else
            {
                /* resample filter memory */
                out_filt_length = modify_Fs_intcub3m_sup( mem_over_hp+2, NS2SA(16000,DELAY_CLDFB_NS), 16000, tmp_buf_switch2, output_Fs, &filt_delay );

                for( i=0; i<filt_delay; i++ )
                {
                    tmp_buf_switch2[out_filt_length+i-filt_delay] = tmp_buf_switch2[out_filt_length-1-filt_delay];
                }

                if ( (bwidth == NB && output_Fs>=16000) || (bwidth > NB && output_Fs>16000) )
                {
                    /* mix cubic and CLDFB resampled buffers in case of resampling to higher frequency rates */
                    for(i=0; i<out_filt_length; i++)
                    {
                        float a = (float) i / out_filt_length;
                        float b = 1.f - a;

                        tmp_buf_switch[NS2SA(output_Fs, SWITCH_GAP_LENGTH_NS) + i] = a * tmp_buf_switch2[i] + b * synth_subfr_out[NS2SA(output_Fs, SWITCH_GAP_LENGTH_NS) + i];
                    }
                }
                else
                {
                    /* copy cubic resampled buffer (memory) */
                    mvr2r( tmp_buf_switch2, tmp_buf_switch+NS2SA(output_Fs,SWITCH_GAP_LENGTH_NS), out_filt_length );
                }
            }
        }
    }

    /* Windowing for overlapadd */
    L = NS2SA(output_Fs, SWITCH_GAP_LENGTH_NS + DELAY_CLDFB_NS); /* 6.25 ms gap + 1.25 ms resamp  */

    set_f( synth, 0, L - Loverlapp);

    if( output_frame == L_FRAME32k )
    {
        for( i=0; i<NS2SA(output_Fs,10000000.0f-N_ZERO_MDCT_NS); i+=2 )
        {
            synth[L-Loverlapp+i] /= win[(3*L_FRAME16k/2-1-i/2)*decimate+decay-(short)(2*((float)L_FRAME48k*N_ZERO_MDCT_NS/FRAME_SIZE_NS))];
            synth[L-Loverlapp+i+1] /= win_int[(3*L_FRAME16k/2-1-i/2)-(short)(2*((float)L_FRAME16k*N_ZERO_MDCT_NS/FRAME_SIZE_NS))];
        }
    }
    else
    {
        for( i=0 ; i<NS2SA(output_Fs,10000000.0f-N_ZERO_MDCT_NS); i++ )
        {
            synth[L-Loverlapp+i] /= win[(3*output_frame/2-1-i)*decimate+decay-(short)(2*((float)L_FRAME48k*N_ZERO_MDCT_NS/FRAME_SIZE_NS))];
        }
    }


    for( i=0; i<Loverlapp; i++ ) /* Windowing for overlapadd */
    {
        synth[i+(short)(L-Loverlapp)] *= (float)sin((i+1)*EVS_PI/(2*(Loverlapp+1)));
        synth_subfr_bwe[i+(short)(L-Loverlapp)-NS2SA(output_Fs,DELAY_CLDFB_NS)] *= (float)cos((i+1)*EVS_PI/(2*(Loverlapp+1)))*(float)cos((i+1)*EVS_PI/(2*(Loverlapp+1)));
        tmp_buf_switch[i+(short)(L-Loverlapp)] *= (float)cos((i+1)*EVS_PI/(2*(Loverlapp+1)))*(float)cos((i+1)*EVS_PI/(2*(Loverlapp+1)));
    }

    /* overlap-add ACELP (tmp_buf_switch) + HQ (synth) */
    for( i=0; i<NS2SA(output_Fs,DELAY_CLDFB_NS); i++ )
    {
        synth[i] += tmp_buf_switch[i];

    }
    for( i=NS2SA(output_Fs,DELAY_CLDFB_NS); i<L; i++ )
    {
        synth[i]+= tmp_buf_switch[i] + synth_subfr_bwe[i-NS2SA(output_Fs,DELAY_CLDFB_NS)];
    }

    return;
}


/*-------------------------------------------------------------------*
* sinq()
*
* Fast sinus generate sin(tmp*i+phi)
*-------------------------------------------------------------------*/

void sinq(
    const float tmp,    /* i  : sinus factor cos(tmp*i+phi)  */
    const float phi,    /* i  : sinus phase cos(tmp*i+phi)  */
    const short N,      /* i  : size of output */
    float x[]           /* o  : output vector  */
)
{
    float A;
    short i;

    x[0] = (float)sin(phi);
    x[1] = (float)sin(tmp+phi);
    x[2] = (float)sin(2*tmp+phi);

    if (fabs((float)(tmp)) > 0.0001f)
    {
        A = (x[2]+x[0])/x[1];
    }
    else
    {
        A = 0;
    }

    for (i=3; i<N; i++)
    {
        x[i] = A*x[i-1]-x[i-2];
    }

    return;
}

