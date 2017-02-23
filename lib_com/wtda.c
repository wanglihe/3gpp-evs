/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"
#include <assert.h>


/*--------------------------------------------------------------------------*
*  mvr2r_inv()
*
*
*--------------------------------------------------------------------------*/

static void mvr2r_inv(
    const float *in,
    float *out,
    short L,
    short decimate
)
{
    short i;

    in=in+(short) ((decimate-1)/2);
    out=out+L-1;

    for (i=0; i<L; i++)
    {
        *out=*in;
        in+=decimate;
        out--;
    }

    return;
}

/*--------------------------------------------------------------------------*
*  mvr2r_dec()
*
*
*--------------------------------------------------------------------------*/

static void mvr2r_dec(
    const float *in,
    float *out,
    short L,
    short decimate
)
{
    short i;
    in=in+(short) ((decimate-1)/2);
    for (i=0; i<L; i++)
    {
        *out=*in;
        in+=decimate;
        out++;
    }

    return;
}

/*--------------------------------------------------------------------------*
*  copy_win()
*
*
*--------------------------------------------------------------------------*/

static void copy_win(
    float *out_win,
    short nb_zero,
    const float *in_win,
    short win_lenght,
    short nb_one,
    short decimate
)
{
    if (decimate<0)
    {
        set_f(out_win,1,nb_one);
        mvr2r_inv(in_win,out_win+nb_one,win_lenght,-decimate);
        set_f(out_win+win_lenght+nb_one,0,nb_zero);
    }
    else
    {
        set_f(out_win,0,nb_zero);
        mvr2r_dec(in_win,out_win+nb_zero,win_lenght,decimate);
        set_f(out_win+nb_zero+win_lenght,1,nb_one);
    }

    return;
}

/*--------------------------------------------------------------------------*
*  tcx_get_windows_mode1()
*
*
*--------------------------------------------------------------------------*/

void tcx_get_windows_mode1(
    const short left_mode,          /* i: overlap mode of left window half          */
    const short right_mode,         /* i: overlap mode of right window half         */
    float  *left_win,               /* o: left overlap window                       */
    float  *right_win,              /* o: right overlap window                      */
    float  *left_win_int,           /* o: left overlap window                       */
    float  *right_win_int,          /* o: right overlap window                      */
    short const L
)
{
    /* Left part */
    if (left_mode == MIN_OVERLAP || left_mode == TRANSITION_OVERLAP)
    {
        if (L==256 || L == 512)
        {
            copy_win(left_win,R1_25-4*R2_25/7,small_overlap_25,R2_25/7,3*R2_25/7,1);
        }
        else
        {
            copy_win(left_win,R1_48-4*R2_48/7,small_overlap_48,R2_48/7,3*R2_48/7,1);
            copy_win(left_win_int,R1_16-4*R2_16/7,small_overlap_int,R2_16/7,3*R2_16/7,1);
        }
    }
    else if (left_mode == HALF_OVERLAP)
    {
        if (L==256 || L == 512)
        {
            copy_win(left_win,R1_25-5*R2_25/7,half_overlap_25,3*R2_25/7,2*R2_25/7,1);
        }
        else
        {
            copy_win(left_win,R1_48-5*R2_48/7,half_overlap_48,3*R2_48/7,2*R2_48/7,1);
            copy_win(left_win_int,R1_16-5*R2_16/7,half_overlap_int,3*R2_16/7,2*R2_16/7,1);
        }
    }
    else if (left_mode == ALDO_WINDOW)
    {
        /* ALDO */
        if (L==256 || L == 512)
        {
            mvr2r(window_256kHz,left_win,R1_25);
        }
        else
        {
            mvr2r(window_48kHz,left_win,R1_48);
            mvr2r(window_8_16_32kHz,left_win_int,R1_16);
        }
    }
    else
    {
        assert(!"Window not supported");
    }

    /* Right part */
    if (right_mode == MIN_OVERLAP || right_mode == TRANSITION_OVERLAP)
    {

        if (L==256 || L == 512)
        {
            copy_win(right_win,3*R2_25/7,small_overlap_25,R2_25/7,3*R2_25/7,-1);
        }
        else
        {
            copy_win(right_win,3*R2_48/7,small_overlap_48,R2_48/7,3*R2_48/7,-1);
            copy_win(right_win_int,3*R2_16/7,small_overlap_int,R2_16/7,3*R2_16/7,-1);
        }
    }
    else if (right_mode == HALF_OVERLAP)
    {

        if (L==256 || L == 512)
        {
            copy_win(right_win,2*R2_25/7,half_overlap_25,3*R2_25/7,2*R2_25/7,-1);
        }
        else
        {
            copy_win(right_win,2*R2_48/7,half_overlap_48,3*R2_48/7,2*R2_48/7,-1);
            copy_win(right_win_int,2*R2_16/7,half_overlap_int,3*R2_16/7,2*R2_16/7,-1);
        }
    }
    else if (right_mode == ALDO_WINDOW)
    {
        if (L==256 || L == 512)
        {
            mvr2r(window_256kHz+R1_25,right_win,R2_25);
        }
        else
        {
            mvr2r(window_48kHz+R1_48,right_win,R2_48);
            mvr2r(window_8_16_32kHz+R1_16,right_win_int,R2_16);
        }
    }
    else
    {
        assert(!"Window not supported");
    }

    return;
}


/*--------------------------------------------------------------------------*
*  wtda()
*
*  Windowing and time-domain aliasing
*--------------------------------------------------------------------------*/

void wtda(
    const float *new_audio,         /* i  : input audio                         */
    float *wtda_audio,        /* o  : windowed audio                      */
    float *old_wtda,          /* i/o: windowed audio from previous frame  */
    const short left_mode,
    const short right_mode,         /* window overlap of current frame (0: full, 2: none, or 3: half) */
    const short L                   /* i  : length                              */
)
{
    short i, decimate, decay;
    short n, windecay48, windecay16;
    const float *allsig_l, *allsig_r;
    float win_right[R2_48];
    float win_int_left[R1_16];
    float win_left[R1_48];
    float win_int_right[R2_16];

    tcx_get_windows_mode1(   left_mode,  right_mode, win_left, win_right, win_int_left,win_int_right, L);


    decimate = 1; /* L_FRAME 48k */
    decay = 0;
    windecay48 = (short)(2*((float)L_FRAME48k*N_ZERO_MDCT_NS/FRAME_SIZE_NS))+R1_48;


    if (L== L_FRAME32k || L== L_FRAME16k )
    {
        decimate = 3;
        decay = 1;

    }
    else if (L== L_FRAME8k)
    {
        decimate = 6;
        decay = 2;

    }


    n = (short)((float)L*N_ZERO_MDCT_NS/FRAME_SIZE_NS);

    windecay16 = (short)(2*((float)L_FRAME16k*N_ZERO_MDCT_NS/FRAME_SIZE_NS))+R1_16;

    /* algorithmic delay reduction */
    i = 0;

    if (old_wtda == NULL)
    {
        allsig_r = new_audio+n;
        allsig_l = new_audio+n-L;
    }
    else
    {
        allsig_r = new_audio+n;
        allsig_l = old_wtda+n;
    }


    if  ( L == L_FRAME32k )
    {
        {
            for (i=0; i<L/2-n; i+=2)
            {
                wtda_audio[i]=-allsig_r[L/2-i-1]*win_int_right[3*L_FRAME16k/2-i/2-1-windecay16]-allsig_r[L/2+i]*win_int_right[3*L_FRAME16k/2+i/2-windecay16];
                wtda_audio[i+1]=-allsig_r[L/2-(i+1)-1]*win_right[(3*L_FRAME16k/2-i/2-1)*decimate+decay-windecay48]-allsig_r[L/2+i+1]*win_right[(3*L_FRAME16k/2+1+i/2)*decimate-decay-1-windecay48];
            }

            for (i=L/2-n; i<L/2; i+=2)
            {
                wtda_audio[i]=-allsig_r[L/2-i-1];
                wtda_audio[i+1]=-allsig_r[L/2-(i+1)-1];
            }
            for (i=0; i<n; i+=2)
            {
                wtda_audio[i+L/2]=allsig_l[i]*win_left[(i/2)*decimate+decay]-new_audio[n-i-1];
                wtda_audio[i+L/2+1]=allsig_l[i+1]*win_int_left[i/2]-new_audio[n-(i+1)-1];
            }

            for (i=n; i<L/2; i+=2)
            {
                wtda_audio[i+L/2]=allsig_l[i]*win_left[(i/2)*decimate+decay]-allsig_l[L-i-1]*win_left[(L/2-i/2)*decimate-1-decay];
                wtda_audio[i+L/2+1]=allsig_l[i+1]*win_int_left[i/2]-allsig_l[L-(i+1)-1]*win_int_left[L/2-i/2-1];
            }

        }
    }
    else
    {
        {

            for (i=0; i<L/2-n; i++)
            {
                wtda_audio[i]=-allsig_r[L/2-i-1]*win_right[3*L/2*decimate-(i+1)*decimate+decay-windecay48]-allsig_r[L/2+i]*win_right[3*L/2*decimate-1+(i+1)*decimate-decay-windecay48];
            }

            for (i=L/2-n; i<L/2; i++)
            {
                wtda_audio[i]=-allsig_r[L/2-i-1];
            }

            for (i=0; i<n; i++)
            {
                wtda_audio[i+L/2]=allsig_l[i]*win_left[i*decimate+decay]-new_audio[n-i-1];
            }

            for (i=n; i<L/2; i++)
            {
                wtda_audio[i+L/2]=allsig_l[i]*win_left[i*decimate+decay]-allsig_l[L-i-1]*win_left[L*decimate-i*decimate-1-decay];
            }

        }

    }

    if (old_wtda != NULL)
    {
        mvr2r( new_audio, old_wtda, L );
    }

    return;
}
