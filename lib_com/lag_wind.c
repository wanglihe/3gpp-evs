/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <assert.h>
#include "options.h"
#include "prot.h"
#include "cnst.h"
#include "rom_com.h"

/*-------------------------------------------------------------------*
 * Local constants
 *-------------------------------------------------------------------*/

#define kLagWinThGain1      0.6f
#define kLagWinThGain2      0.3f


/*-------------------------------------------------------------*
 * procedure lag_wind()                                        *
 *           ~~~~~~~~~                                         *
 * lag windowing of the autocorrelations                       *
 *-------------------------------------------------------------*/

void lag_wind(
    float r[],           /* i/o: autocorrelations                       */
    const short m,       /* i  : order of LP filter                     */
    const int sr,        /* i  : sampling rate                          */
    const short strength /* i  : LAGW_WEAK, LAGW_MEDIUM, or LAGW_STRONG */
)
{
    short i;
    const float *wnd;

    assert(0 <= strength && strength <= NUM_LAGW_STRENGTHS);

    switch (sr)
    {
    case 8000:
        assert(m <= 16);
        assert(strength == LAGW_STRONG);
        wnd = lag_window_8k;
        break;
    case 12800:
        assert(m <= 16);
        wnd = lag_window_12k8[strength];
        break;
    case 16000:
        assert(m <= 16);
        wnd = lag_window_16k[strength];
        break;
    case 24000:
    case 25600:
        assert(m <= 16);
        wnd = lag_window_25k6[strength];
        break;
    case 32000:
        assert(m <= 16);
        wnd = lag_window_32k[strength];
        break;
    case 48000:
        assert(m <= 16);
        assert(strength == LAGW_STRONG);
        wnd = lag_window_48k;
        break;
    default:
        assert(!"Lag window not implemented for this sampling rate");
        return;
    }

    for( i=0; i<=m; ++i )
    {
        r[i] *= wnd[i];
    }

    return;
}

/*-------------------------------------------------------------*
 * procedure adapt_lag_wind()
 *
 *
 *-------------------------------------------------------------*/

void adapt_lag_wind(
    float r[],            /* i/o: autocorrelations                                       */
    int m,                /* i  : order of LP filter                                     */
    const int Top,        /* i  : open loop pitch lags from curr. frame (or NULL if n/a) */
    const float Tnc,      /* i  : open loop pitch gains from curr. frame (NULL if n/a)   */
    int sr                /* i  : sampling rate                                          */
)
{
    short strength;
    short pitch_lag;
    float pitch_gain;

    pitch_lag = (short)Top;
    pitch_gain = (float)Tnc;

    if (pitch_lag < 80)
    {
        if (pitch_gain > kLagWinThGain1)
        {
            strength = LAGW_STRONG;
        }
        else
        {
            strength = LAGW_MEDIUM;
        }
    }
    else if (pitch_lag < 160)
    {
        if (pitch_gain > kLagWinThGain2)
        {
            strength = LAGW_MEDIUM;
        }
        else
        {
            strength = LAGW_WEAK;
        }
    }
    else
    {
        strength = LAGW_WEAK;
    }

    lag_wind( r, m, sr, strength );

    return;
}
