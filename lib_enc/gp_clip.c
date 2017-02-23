/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "prot.h"
#include "options.h"
#include "prot.h"
#include "cnst.h"

/*-------------------------------------------------------------------*
* Local constants
*-------------------------------------------------------------------*/

#define DIST_ISF_MAX_IO    150.0f

#define DIST_ISF_MAX    120.0f
#define DIST_ISF_THRES  60
#define GAIN_PIT_THRES  0.9f
#define GAIN_PIT_MIN    0.6f

#define ALPHA1          0.98f
#define ALPHA4          0.99f
#define WINDOW_SIZE     50
#define THRESH_TYPE     0.85f
#define THRESH_VOICING  0.86f

/*-------------------------------------------------------------------*
* init_gp_clip
*
* Pitch Gain clipping initializations
*-------------------------------------------------------------------*/
void init_gp_clip(
    float mem[]  /* o:   memory of gain of pitch clipping algorithm */
)
{
    mem[0] = DIST_ISF_MAX;
    mem[1] = GAIN_PIT_MIN;
    mem[2] = 0.0f;              /* old energy of target (dB) */
    mem[3] = 0.0f;
    mem[4] = 0.0f;
    mem[5] = 0.8f;

    return;
}

/*-------------------------------------------------------------------*
 * Function gp_clip
 *
 * The gain needs to be limited (gain pitch < 1.0) when one of the
 * following cases occurs:
 * - a resonance on LPC filter (lp_disp < 60 Hz)  AND a good pitch
 *   prediction (lp_gp > 0.9)
 * - target energy drops by 6 dB AND a good pitch prediction (lp_gp>1.0)
 *-------------------------------------------------------------------*/

short gp_clip(
    const long  core_brate,  /* i  : core bitrate                               */
    const float *voicing,    /* i  : normalized correlations (from OL pitch)    */
    const short i_subfr,     /* i  : subframe index                             */
    const short coder_type,  /* i  : type of coder                              */
    const float xn[],        /* i  : target vector                              */
    float mem[]        /* i/o: memory of gain of pitch clipping algorithm */
)
{
    short clip;
    short i;
    float wener, tmp;

    clip = 0;
    if( (core_brate == ACELP_6k60) ||  (core_brate == ACELP_8k85) )
    {
        tmp = 0.9f + (0.1f * mem[0]/DIST_ISF_MAX_IO);	/* clipping is activated when filtered pitch gain > threshold (0.94 to 1) */
        if (mem[1] > tmp)
        {
            clip = 1;
        }
    }
    else if( (mem[0] < DIST_ISF_THRES) && (mem[1] > GAIN_PIT_THRES) )
    {
        clip = 1;
    }

    wener = 0.01f;
    for (i=0; i<L_SUBFR; i++)
    {
        wener += xn[i]*xn[i];
    }
    wener = 10.0f*(float)log10(wener);
    if ((wener < (mem[2]-6.0f)) && (mem[1] > 1.0f))
    {
        clip = 1;
    }

    mem[2] = wener;
    tmp = ALPHA1 * mem[4];

    if( coder_type == GENERIC || coder_type == TRANSITION || coder_type == INACTIVE )
    {
        tmp += (1-ALPHA1);
    }

    mem[4] = tmp;
    tmp = ALPHA4 * mem[5];
    if( i_subfr == 0 )
    {
        mem[5] = (1-ALPHA4) * voicing[0] + tmp;
    }
    else if( i_subfr == 2*L_SUBFR )
    {
        mem[5] = (1-ALPHA4) * voicing[1] + tmp;
    }
    if( mem[3] > WINDOW_SIZE )
    {
        if( ( mem[4] > THRESH_TYPE ) && ( mem[5] > THRESH_VOICING ) )
        {
            clip = 1;
        }
    }
    else
    {
        mem[3]++;
    }

    return (clip);
}

/*-------------------------------------------------------------------*
 * gp_clip_test_lsf()
 *
 * check the minimum distance of LSFs for pitch gain clipping flag
 *-------------------------------------------------------------------*/

void gp_clip_test_lsf(
    const long  core_brate,  /* i  : core bitrate                               */
    const float lsf[],        /* i  : LSF vector                                 */
    float mem[],        /* i/o: memory of gain of pitch clipping algorithm */
    const short Opt_AMR_WB    /* i  : flag indicating AMR-WB IO mode             */
)
{
    short i;
    short m;
    float dist, dist_min;

    dist_min = lsf[1] - lsf[0];

    if ( Opt_AMR_WB )
    {
        m = M-1;
    }
    else
    {
        m = M;
    }

    for (i=2; i<m; i++)
    {
        dist = lsf[i] - lsf[i-1];
        if (dist < dist_min)
        {
            dist_min = dist;
        }
    }

    dist = 0.8f*mem[0] + 0.2f*dist_min;

    if( (core_brate == ACELP_6k60) ||  (core_brate == ACELP_8k85) )
    {
        if (dist > DIST_ISF_MAX_IO)
        {
            dist = DIST_ISF_MAX_IO;
        }
    }
    else if (dist > DIST_ISF_MAX)
    {
        dist = DIST_ISF_MAX;
    }

    mem[0] = dist;

    return;
}

/*-------------------------------------------------------------------*
 * gp_clip_test_gain_pit()
 *
 * low-pass filtering of the pitch gain for pitch gain clipping flag
 *-------------------------------------------------------------------*/

void gp_clip_test_gain_pit(
    const long  core_brate,  /* i  : core bitrate                               */
    const float gain_pit,  /* i  :   gain of quantized pitch                  */
    float mem[]      /* i/o: memory of gain of pitch clipping algorithm */
)
{
    float gain;

    if( (core_brate == ACELP_6k60) || (core_brate == ACELP_8k85) )
    {
        gain = 0.98f*mem[1] + 0.02f*gain_pit;	  /* long term LTP gain average (>250ms) */
    }
    else
    {
        gain = 0.9f*mem[1] + 0.1f*gain_pit;
    }
    if (gain < GAIN_PIT_MIN)
    {
        gain = GAIN_PIT_MIN;
    }
    mem[1] = gain;

    return;
}
