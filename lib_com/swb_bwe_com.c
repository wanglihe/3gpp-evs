/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"

/*-------------------------------------------------------------------*
 * WB_BWE_gain_pred()
 *
 * predict WB frequency envelopes for 0b WB BWE
 *-------------------------------------------------------------------*/

short WB_BWE_gain_pred(
    float *WB_fenv,                 /* o  : WB frequency envelopes               */
    const float *core_dec_freq,           /* i  : Frequency domain core decoded signal */
    const short coder_type,               /* i  : coding type                          */
    short prev_coder_type,          /* i  : coding type of last frame            */
    float prev_WB_fenv,             /* i  : envelope for last frame              */
    float *voice_factors,           /* i  : voicing factors                      */
    const float pitch_buf[],              /* i  : pitch buffer                         */
    long  last_core_brate,          /* i  : previous frame core bitrate          */
    float last_wb_bwe_ener          /* i  : previous frame wb bwe signal energy  */
    ,short last_extl                 /* i  : extl. layer for last frame           */
    ,float tilt
)
{
    float enerL, alfa = 1.0f;
    short n_freq, mode;
    short ener_var_flag = 0;
    float voice_factor, pitch;
    short env_var_flag = 0;

    mode = NORMAL;

    enerL = EPSILON;
    for (n_freq = 128; n_freq<192; n_freq++)
    {
        enerL += core_dec_freq[n_freq] * core_dec_freq[n_freq];
    }
    WB_fenv[0] = EPSILON;
    for (n_freq = 192; n_freq<224; n_freq++)
    {
        WB_fenv[0] += core_dec_freq[n_freq] * core_dec_freq[n_freq];
    }

    WB_fenv[1] = EPSILON;
    for (n_freq = 224; n_freq<256; n_freq++)
    {
        WB_fenv[1] += core_dec_freq[n_freq] * core_dec_freq[n_freq];
    }

    voice_factor = sum_f( voice_factors, 4 );
    pitch = sum_f( pitch_buf, 4 ) + EPSILON;

    if(enerL < 16.0f*max(WB_fenv[0], WB_fenv[1]) && pitch < 308)
    {
        ener_var_flag = 1;
    }

    if(WB_fenv[0] > 2.0f*WB_fenv[1])
    {
        alfa = max(2.0f*WB_fenv[1]/WB_fenv[0], 0.1f);
        WB_fenv[0] *= alfa;
    }
    else if (2.0f*WB_fenv[0] < WB_fenv[1] && coder_type != UNVOICED)
    {
        alfa = max(2.0f*WB_fenv[0]/WB_fenv[1], 0.1f);
        WB_fenv[1] *= alfa;
    }

    WB_fenv[0] = (float)sqrt((WB_fenv[0]+WB_fenv[1])/64);

    if(coder_type != AUDIO && coder_type != UNVOICED && ener_var_flag == 0)
    {
        WB_fenv[0] *= 1.5f;
    }

    if( coder_type != TRANSITION && coder_type != AUDIO && coder_type != UNVOICED && sqrt(enerL) > 40.0f*WB_fenv[0] && alfa > 0.9f &&
            !(coder_type == prev_coder_type && WB_fenv[0] > prev_WB_fenv) )
    {
        WB_fenv[0] *= min((float)(0.025f*sqrt(enerL)/WB_fenv[0]), 4.0f);

        if( WB_fenv[0] > prev_WB_fenv )
        {
            WB_fenv[0] = 0.3f*WB_fenv[0] + 0.7f*prev_WB_fenv;
        }
    }

    alfa = min(1.5f, max(0.5f, 77.0f*voice_factor/pitch));
    if( sqrt(enerL) > 64.0f*alfa*WB_fenv[0] && 3.0f*WB_fenv[0]*WB_fenv[0] < sqrt(enerL) && prev_coder_type != UNVOICED )
    {
        env_var_flag = 1;
        WB_fenv[0] *= min((float)(0.015625f*sqrt(enerL)/WB_fenv[0]), 4.0f);

        if( WB_fenv[0] > prev_WB_fenv )
        {
            WB_fenv[0] = 0.3f*WB_fenv[0] + 0.7f*prev_WB_fenv;
        }
    }

    if( coder_type == UNVOICED || prev_coder_type == UNVOICED )
    {
        WB_fenv[0] *= 0.5f;
    }

    if( coder_type != AUDIO )
    {
        WB_fenv[0] /= max(1.2f * voice_factor, 1.0f);
        WB_fenv[0] *= min(2.0f, max(0.125f, pitch/400.0f));
    }

    if( last_core_brate > ACELP_8k00 && WB_fenv[0] > last_wb_bwe_ener )
    {
        WB_fenv[0] = 0.9f*last_wb_bwe_ener + 0.1f*WB_fenv[0];
    }

    if(last_extl != WB_BWE && (tilt < 8.f))
    {
        WB_fenv[0] *= min(0.5, 16.0f*tilt);
    }

    if(env_var_flag == 1)
    {
        WB_fenv[1] = 1.5f*WB_fenv[0];
        WB_fenv[0] *= 0.75f;
    }
    else
    {
        WB_fenv[1] = WB_fenv[0];
    }

    if(coder_type == UNVOICED || prev_coder_type == UNVOICED)
    {
        WB_fenv[1] *= 0.5f;
    }

    return (mode);
}

/*-------------------------------------------------------------------*
 * calc_normal_length()
 *
 *-------------------------------------------------------------------*/

void calc_normal_length(
    const short core,             /* i  : core                   */
    const float *sp,              /* i  : input signal           */
    const short mode,             /* i  : input mode             */
    const short extl,             /* i  : extension layer        */
    short *L_swb_norm,      /* o  : normalize length       */
    short *prev_L_swb_norm  /*i/o : last normalize length  */
)
{
    short i, n_freq, n_band, THRES;
    const float *pit;
    float peak, mean, mag;
    short L_swb_norm_trans, L_swb_norm_norm, L_swb_norm_harm, L_swb_norm_cur;
    short N;

    if( core == HQ_CORE || extl == SWB_BWE || extl == FB_BWE )
    {
        THRES = 8;
    }
    else
    {
        THRES = 4;
    }

    if( core == HQ_CORE && (mode == HQ_HARMONIC || mode == HQ_HVQ) )
    {
        N = 13;
    }
    else
    {
        N = 16;
    }

    n_band = 0;
    pit = sp;
    for(i = 0; i < N; i ++)
    {
        peak = 0.0f;
        mean = 0;

        for(n_freq = 0; n_freq < 16; n_freq ++)
        {
            mag = (float) fabs(*pit);
            if (mag > peak)
            {
                peak = mag;
            }
            mean += mag;
            pit ++;
        }

        if((15+THRES)*peak > THRES*mean && peak>10)
        {
            n_band += 1;
        }
    }

    if( core == ACELP_CORE )
    {
        L_swb_norm_trans = (short)(4 + 0.25f*n_band);
        L_swb_norm_norm = (short)(8 + 0.5f*n_band);
        L_swb_norm_harm = max((short)(32 + 2.0f*n_band), 24);

        if( mode == HARMONIC )
        {
            L_swb_norm_cur = L_swb_norm_harm;
        }
        else if( mode == NORMAL )
        {
            L_swb_norm_cur = L_swb_norm_norm;
        }
        else
        {
            L_swb_norm_cur = L_swb_norm_trans;
        }

        *L_swb_norm = (short)(0.5f*L_swb_norm_cur + 0.5f* (*prev_L_swb_norm));
        *prev_L_swb_norm = L_swb_norm_cur;
    }
    else
    {
        if( mode == HQ_HARMONIC || mode == HQ_HVQ )
        {
            L_swb_norm_cur = (short)( 32 + 2.5f*n_band);
        }
        else
        {
            L_swb_norm_cur = (short)(8 + 0.5f*n_band);
        }

        *L_swb_norm = (short)(0.1f*L_swb_norm_cur + 0.9f* (*prev_L_swb_norm) + 0.5f);
        *prev_L_swb_norm = L_swb_norm_cur;
    }

    return;
}
/*-------------------------------------------------------------------*
* calc_tilt_bwe()
*
* calculate tilt parameter
*-------------------------------------------------------------------*/

void calc_tilt_bwe(
    const float *sp,     /* i  : input signal    */
    float *tilt,   /* o  : signal tilt     */
    const short N        /* i  : signal length   */
)
{
    short i;
    float r0, r1;

    r0 = EPSILON;
    for(i=0; i<N; i++)
    {
        r0 += sp[i] * sp[i];
    }

    r1 = (float)fabs(sp[1] - sp[0]);
    for(i=2; i<N; i++)
    {
        if((sp[i] - sp[i-1]) * (sp[i-1] - sp[i-2]) < 0)
        {
            r1 += (float)fabs(sp[i] - sp[i-1]);
        }
    }

    *tilt = (float)(r1/sqrt(r0));

    return;
}

/*-------------------------------------------------------------------*
* calc_norm_envelop()
*
* calculate normalized parameter
*-------------------------------------------------------------------*/

void calc_norm_envelop(
    const float SWB_signal[],    /* i  : SWB spectrum            */
    float *envelope,       /* o  : normalized envelope     */
    const short L_swb_norm,      /* i  : length of envelope      */
    const short SWB_flength,     /* i  : Length of input/output  */
    const short st_offset        /* i  : offset                  */
)
{
    short i, lookback, env_index, n_freq, n_lag_now, n_lag;

    lookback = L_swb_norm/2;
    env_index = swb_bwe_subband[0]+st_offset;
    n_lag_now = L_swb_norm;
    for (n_freq = swb_bwe_trans_subband[0]+st_offset-lookback; n_freq<SWB_flength+st_offset-L_swb_norm; n_freq++)
    {
        /* Apply MA filter */
        envelope[env_index] = EPSILON;
        for (n_lag=0; n_lag<n_lag_now; n_lag++)
        {
            envelope[env_index] += (float)fabs(SWB_signal[n_freq+n_lag]);
        }
        env_index++;
    }

    for (n_freq = SWB_flength+st_offset-L_swb_norm, i = 0; n_freq<SWB_flength+st_offset-lookback; n_freq++, i++)
    {
        n_lag_now = L_swb_norm-i;
        /* Apply MA filter */
        envelope[env_index] = EPSILON;
        for (n_lag=0; n_lag<n_lag_now; n_lag++)
        {
            envelope[env_index] += (float)fabs(SWB_signal[n_freq+n_lag]);
        }
        env_index++;
    }

    return;
}


/*-------------------------------------------------------------------*
* calc_norm_envelop_lf()
*
* calc_envelope of low frequency spectrum
*-------------------------------------------------------------------*/
static
void calc_norm_envelop_lf(
    const float SWB_signal[],        /* i  : SWB spectrum                                    */
    float *envelope,                 /* o  : normalized envelope                             */
    short *L_swb_norm,               /* i/o: length of envelope                              */
    const short HQ_mode,             /* i  : HQ mode                                         */
    const short hq_generic_offset,   /* i  : frequency offset for representing hq swb bwe    */
    short *sfreq,              /* i  : starting frequency index                        */
    short *efreq)              /* i  : ending frequency index                          */
{
    short lookback, env_index, n_freq, n_lag_now, n_lag;

    *sfreq = 2;
    if ( hq_generic_offset == HQ_GENERIC_FOFFSET_24K4)
    {
        *efreq = 146;
        if ( HQ_mode == HQ_GEN_FB)
        {
            *efreq = 306;
        }
        if ( (328 - *efreq)*2 + 1 < *L_swb_norm)
        {
            *L_swb_norm = (328 - *efreq)*2 + 1;
        }
    }
    else
    {
        *efreq = 130;
        if (HQ_mode == HQ_GEN_FB)
        {
            *efreq = 290;
        }
        if ( (400 - *efreq)*2 + 1 < *L_swb_norm)
        {
            *L_swb_norm = (400 - *efreq)*2 + 1;
        }
    }

    lookback = *L_swb_norm/2;
    env_index = 0;
    n_lag_now = *L_swb_norm;
    for (n_freq = 0; n_freq<lookback; n_freq++)
    {
        envelope[env_index] = EPSILON;
        for (n_lag=0; n_lag<lookback+n_freq; n_lag++)
        {
            envelope[env_index] += (float)fabs(SWB_signal[n_lag]);
        }
        env_index++;
    }
    for (; n_freq<*efreq; n_freq++)
    {
        /* Apply MA filter */
        envelope[env_index] = EPSILON;
        for (n_lag=0; n_lag<n_lag_now; n_lag++)
        {
            envelope[env_index] += (float)fabs(SWB_signal[n_freq-lookback+n_lag]);
        }
        env_index++;
    }
    return;
}

/*-------------------------------------------------------------------*
 * WB_BWE_decoding()
 *
 * WB BWE decoder
 *-------------------------------------------------------------------*/

void WB_BWE_decoding(
    const float *core_dec_freq,    /* i  : Frequency domain core decoded signal */
    float *WB_fenv,          /* i  : WB frequency envelopes               */
    float *WB_signal,        /* o  : WB signal in MDCT domain             */
    const short WB_flength,        /* i  : Length of input/output               */
    const short mode,              /* i  : classification for WB signal         */
    const short last_extl,         /* i  : extl. layer for last frame           */
    float *prev_Energy,      /* i/o: energy for last frame                */
    float *prev_WB_fenv,     /* i/o: envelope for last frame              */
    short *prev_L_wb_norm,   /* i/o: length for last frame wb norm        */
    const short extl,              /* i  : extension layer                      */
    const short coder_type,        /* i  : coding type                          */
    const long  total_brate,       /* i  : core layer bitrate                   */
    short *Seed,             /* i/o: random generator seed                */
    short *prev_flag,        /* i/o: attenu flag of last frame            */
    short prev_coder_type    /* i  : coding type of last frame            */
)
{
    short n_freq, n_band;
    short i, L;
    float envelope[L_FRAME16k];
    float energy, wfenv[2], EnergyL;
    float *pit1;
    short L_wb_norm;
    float alfa, beta;
    short flag = 0;
    short core_type = 1;
    short signum[L_FRAME16k];
    float inv_L_wb_norm, weight;

    calc_normal_length( ACELP_CORE, core_dec_freq, mode, extl, &L_wb_norm, prev_L_wb_norm );

    set_f( WB_signal, 0, L_FRAME16k );

    /* copy excitation */
    if( coder_type != AUDIO && total_brate <= ACELP_8k00 )
    {
        core_type = 0;
    }

    if( core_type == 0 )
    {
        mvr2r(&core_dec_freq[160], &WB_signal[240], 80);
    }
    else
    {
        mvr2r(&core_dec_freq[80], &WB_signal[240], 80);
    }

    /* calculate envelope */
    calc_norm_envelop(WB_signal, envelope, L_wb_norm, WB_flength, 0);

    if( coder_type != UNVOICED && total_brate <= ACELP_8k00 )
    {
        inv_L_wb_norm = 1.0f/L_wb_norm;
        weight = (mode != HARMONIC) ? max(min(3.0f*inv_L_wb_norm, 0.5f), 0.25f) : 0.25f;
        for(n_freq = swb_bwe_subband[0]; n_freq<swb_bwe_subband[4]; n_freq++)
        {
            signum[n_freq] = 1;
            if (WB_signal[n_freq]<0)
            {
                signum[n_freq] = -1;
                WB_signal[n_freq] *= signum[n_freq];
            }

            WB_signal[n_freq] = WB_signal[n_freq] - 0.45f*envelope[n_freq]*inv_L_wb_norm;
            if(WB_signal[n_freq] > 0)
            {
                WB_signal[n_freq] *= (0.55f - weight);
            }
            WB_signal[n_freq] *= signum[n_freq];
        }
    }

    /* Normalize with envelope */
    for (n_freq = swb_bwe_subband[0]; n_freq<swb_bwe_subband[4]; n_freq++)
    {
        WB_signal[n_freq] /= envelope[n_freq];
    }

    if( mode == HARMONIC )
    {
        L = 4;
    }
    else
    {
        L = 1;
    }

    if( coder_type == UNVOICED )
    {
        for ( n_freq = swb_bwe_subband[0]; n_freq < swb_bwe_subband[4]; n_freq++ )
        {
            *Seed = (short)(12345*(*Seed) + 20101);
            WB_signal[n_freq] = (*Seed)/32768.0f;
        }
    }
    else
    {
        for( n_band = 0; n_band < 4; n_band += L )
        {
            energy = EPSILON;
            for (n_freq = swb_bwe_subband[n_band]; n_freq<swb_bwe_subband[n_band+L]; n_freq++)
            {
                energy += WB_signal[n_freq] * WB_signal[n_freq];
            }

            energy = (float)sqrt((swb_bwe_subband[n_band+L] - swb_bwe_subband[n_band])/energy);

            for (n_freq = swb_bwe_subband[n_band]; n_freq<swb_bwe_subband[n_band+L]; n_freq++)
            {
                WB_signal[n_freq] *= energy;
            }
        }
    }


    EnergyL = 0.0f;
    if( core_type == 1 )
    {
        if( prev_coder_type != AUDIO && total_brate <= ACELP_8k00 )
        {
            for(i=160; i<240; i++)
            {
                EnergyL += (float)fabs(core_dec_freq[i]);
            }
        }
        else
        {
            for(i=80; i<240; i++)
            {
                EnergyL += (float)fabs(core_dec_freq[i]);
            }
        }

        if(total_brate <= ACELP_8k00)
        {
            alfa = 0.8f;
            beta = 1.25f;
        }
        else
        {
            alfa = 0.5f;
            beta = 2.0f;
        }
    }
    else
    {
        if( prev_coder_type == AUDIO )
        {
            for(i=80; i<240; i++)
            {
                EnergyL += (float)fabs(core_dec_freq[i]);
            }
        }
        else
        {
            for(i=160; i<240; i++)
            {
                EnergyL += (float)fabs(core_dec_freq[i]);
            }
        }

        if( prev_coder_type == coder_type && WB_fenv[0] > prev_WB_fenv[0] )
        {
            alfa = 0.4f;
            beta = 2.5f;
        }
        else
        {
            alfa = 0.6f;
            beta = 1.67f;
        }

        if( coder_type == GENERIC || ((EnergyL > 0.5f*(*prev_Energy) && EnergyL < 2.0f*(*prev_Energy) && *prev_flag == 1)) )
        {
            WB_fenv[0] *= 0.5f;
            WB_fenv[1] *= 0.5f;
            flag = 1;
        }
    }
    if( (mode == HARMONIC && WB_fenv[1] < 0.25f*WB_fenv[0]) || mode == NORMAL )
    {
        if( last_extl == WB_BWE && ( (prev_coder_type == AUDIO && coder_type != AUDIO) || (prev_coder_type != AUDIO && coder_type == AUDIO)) && total_brate <= ACELP_8k00 )
        {
            if( WB_fenv[0] > prev_WB_fenv[0] )
            {
                wfenv[0] = 0.3f*WB_fenv[0] + 0.7f*prev_WB_fenv[0];
                wfenv[1] = 0.3f*WB_fenv[1] + 0.7f*prev_WB_fenv[1];
            }
            else
            {
                wfenv[0] = 0.5f*WB_fenv[0] + 0.5f*prev_WB_fenv[0];
                wfenv[1] = 0.4f*WB_fenv[1] + 0.4f*prev_WB_fenv[1];
            }
        }
        else if ( last_extl == WB_BWE && prev_WB_fenv[0]*EnergyL < WB_fenv[0]*(*prev_Energy) && WB_fenv[0] > prev_WB_fenv[0] && coder_type != AUDIO && coder_type != UNVOICED && total_brate <= ACELP_8k00)
        {
            wfenv[0] = 0.3f*WB_fenv[0] + 0.7f*prev_WB_fenv[0];
            wfenv[1] = 0.3f*WB_fenv[1] + 0.7f*prev_WB_fenv[1];
        }
        else if ( last_extl == WB_BWE && EnergyL > alfa*(*prev_Energy) && EnergyL < beta*(*prev_Energy) && prev_coder_type != UNVOICED )
        {
            wfenv[0] = 0.5f*(WB_fenv[0] + prev_WB_fenv[0]);
            wfenv[1] = 0.5f*(WB_fenv[1] + prev_WB_fenv[1]);
        }
        else
        {
            wfenv[0] = WB_fenv[0];
            wfenv[1] = WB_fenv[1];
        }
        for (n_freq = swb_bwe_subband[0]; n_freq<swb_bwe_subband[2]; n_freq++)
        {
            WB_signal[n_freq] *= wfenv[0];
        }

        for (n_freq = swb_bwe_subband[2]; n_freq<swb_bwe_subband[4]; n_freq++)
        {
            WB_signal[n_freq] *= wfenv[1];
        }

        prev_WB_fenv[0] = wfenv[0];
        prev_WB_fenv[1] = wfenv[1];
    }
    else
    {
        wfenv[0] = 0.5f*(WB_fenv[0]+WB_fenv[1]);

        if(last_extl == WB_BWE && EnergyL > 0.5f*(*prev_Energy) && EnergyL < 2.0f*(*prev_Energy))
        {
            wfenv[0] = 0.25f*wfenv[0] + 0.375f*(prev_WB_fenv[0] + prev_WB_fenv[1]);
        }
        for (n_freq = swb_bwe_subband[0]; n_freq<swb_bwe_subband[4]; n_freq++)
        {
            WB_signal[n_freq] *= wfenv[0];
        }

        prev_WB_fenv[0] = wfenv[0];
        prev_WB_fenv[1] = wfenv[0];
    }

    *prev_flag = flag;
    *prev_Energy = EnergyL;
    pit1 = &WB_signal[240];

    for(n_freq=0; n_freq<16; n_freq++)
    {
        *(pit1++) *= (0.2f+n_freq*0.05f);
    }

    if( core_type == 1 )
    {
        pit1 = &WB_signal[280];
        for(n_freq=0; n_freq<40; n_freq++)
        {
            *(pit1++) *= ( 1.0f-n_freq*0.02f);
        }
    }
    else
    {
        pit1 = &WB_signal[300];
        for(n_freq=0; n_freq<20; n_freq++)
        {
            *(pit1++) *= ( 1.0f-n_freq*0.04f);
        }
    }

    return;
}

/*-------------------------------------------------------------------*
 * SWB_BWE_decoding()
 *
 * SWB BWE decoder
 *-------------------------------------------------------------------*/

void SWB_BWE_decoding(
    const float *core_dec_freq,    /* i  : Frequency domain core decoded signal  */
    float *SWB_fenv,         /* i/o: SWB frequency envelopes               */
    float *SWB_signal,       /* o  : SWB signal in MDCT domain             */
    const short SWB_flength,       /* i  : Length of input/output                */
    const short mode,              /* i  : classification for SWB signal         */
    short *frica_flag,       /* o  : fricative signal flag                 */
    float *prev_Energy,      /* i/o: energy for last frame                 */
    float *prev_SWB_fenv,    /* i/o: envelope for last frame               */
    short *prev_L_swb_norm,  /* i/o: length for last frame wb norm         */
    const float tilt_nb,           /* i  : tilt of synthesis wb signal           */
    short *Seed,             /* i/o: random generator seed                 */
    const short st_offset,         /* i  : offset value due to different core    */
    float *prev_weight,      /* i/o: excitation weight value of last frame */
    const short extl               /* i  : extension layer                       */
    ,const short last_extl          /* i  : extension layer of last frame         */
)
{
    short n_freq, n_band, L, L_swb_norm;
    float *pit1;
    float envelope[L_FRAME32k];
    float fenvL, EnergyL, Energy, energy, weight, wfenv, factor;
    float mean, factor1, tmp1, tmp2, tmp3, tmp4, tmp_ener;
    short signum[L_FRAME32k];
    float inv_L_swb_norm;

    fenvL = 0.0f;
    EnergyL = EPSILON;
    for(n_freq=224+st_offset; n_freq<swb_bwe_trans_subband[0]+st_offset; n_freq++)
    {
        fenvL += core_dec_freq[n_freq] * core_dec_freq[n_freq];
    }

    for( n_freq = 16; n_freq < L_FRAME; n_freq++ )
    {
        EnergyL += core_dec_freq[n_freq] * core_dec_freq[n_freq];
    }
    fenvL = (float)sqrt(fenvL/16);
    EnergyL = (float)sqrt(EnergyL/240);
    if (fenvL > 8.0f*SWB_fenv[0])
    {
        fenvL = SWB_fenv[0];
    }
    calc_normal_length( ACELP_CORE, core_dec_freq, mode, extl, &L_swb_norm, prev_L_swb_norm );

    if( mode == TRANSIENT )
    {
        Energy = 0.0f;
        for(n_band = 0; n_band < SWB_FENV_TRANS; n_band++)
        {
            Energy += SWB_fenv[n_band]*SWB_fenv[n_band];
        }
        Energy /= SWB_FENV_TRANS;

        /* Reconstruct excitation from LF signal */
        mvr2r(&core_dec_freq[112], &SWB_signal[240+st_offset], 128);
        mvr2r(&core_dec_freq[112], &SWB_signal[368+st_offset], 128);
        mvr2r(&core_dec_freq[176], &SWB_signal[496+st_offset], 64);

        /* calculate envelope */
        calc_norm_envelop(SWB_signal, envelope, L_swb_norm, SWB_flength, st_offset);

        /* Normalize with envelope */
        for (n_freq = swb_bwe_trans_subband[0]+st_offset; n_freq<swb_bwe_trans_subband[SWB_FENV_TRANS]+st_offset; n_freq++)
        {
            SWB_signal[n_freq] /= envelope[n_freq];
        }

        for(n_band=0; n_band<SWB_FENV_TRANS; n_band++)
        {
            energy = EPSILON;
            for (n_freq = swb_bwe_trans_subband[n_band]+st_offset; n_freq<swb_bwe_trans_subband[n_band+1]+st_offset; n_freq++)
            {
                energy += SWB_signal[n_freq] * SWB_signal[n_freq];
            }

            tmp_ener = (float)(sqrt(swb_bwe_trans_subband_width[n_band] / energy) * SWB_fenv[n_band]);

            for (n_freq = swb_bwe_trans_subband[n_band]+st_offset; n_freq<swb_bwe_trans_subband[n_band+1]+st_offset; n_freq++)
            {
                SWB_signal[n_freq] *= tmp_ener;
            }
        }

        for( n_band = 0; n_band < 8; n_band++ )
        {
            prev_SWB_fenv[n_band] = SWB_fenv[n_band/4]*SWB_fenv[n_band/4];
        }

        for( n_band = 0; n_band < 6; n_band++ )
        {
            prev_SWB_fenv[8+n_band] = SWB_fenv[2+n_band/3]*SWB_fenv[2+n_band/3];
        }

        *prev_weight = 0.5f;
    }
    else
    {
        Energy = EPSILON;
        for(n_band = 0; n_band < SWB_FENV; n_band++)
        {
            Energy += SWB_fenv[n_band];
        }
        Energy /= SWB_FENV;
        if(last_extl != SWB_BWE && last_extl != FB_BWE)
        {
            if(16.0f*Energy < EnergyL && extl == FB_BWE)
            {
                for(n_band=0; n_band <SWB_FENV; n_band++)
                {
                    SWB_fenv[n_band] *= 0.2f;
                }
                fenvL *= 0.2f;
            }
            mvr2r(SWB_fenv, prev_SWB_fenv, SWB_FENV);
        }
        if( mode == HARMONIC )
        {
            mvr2r( core_dec_freq, &SWB_signal[240+st_offset], 240 );
            mvr2r( &core_dec_freq[128], &SWB_signal[480+st_offset], 80 );

            /* calculate envelope */
            calc_norm_envelop(SWB_signal, envelope, L_swb_norm, SWB_flength, st_offset);
        }
        else
        {
            if(mode == NOISE || ((Energy > EnergyL || (tilt_nb > 7 && Energy > 0.5f*EnergyL) || tilt_nb > 12) && Energy > 75 && fenvL > 25))
            {
                for (n_freq=swb_bwe_subband[0]+st_offset; n_freq<swb_bwe_subband[SWB_FENV]+st_offset; n_freq++)
                {
                    *Seed = (short)(12345*(*Seed) + 20101);
                    SWB_signal[n_freq] = (*Seed)/32768.0f;
                }
                if(mode != NOISE)
                {
                    *frica_flag = 1;
                }
            }
            else
            {
                /* modify SHB frequency envelopes when SHB spectrum is unflat */
                for(n_band=0; n_band <13; n_band++)
                {
                    if( SWB_fenv[n_band] * 0.9f > SWB_fenv[n_band+1] )
                    {
                        SWB_fenv[n_band+1] *= (0.8f+0.015f*n_band);
                    }
                    else if( SWB_fenv[n_band+1] * 0.9f > SWB_fenv[n_band] )
                    {
                        SWB_fenv[n_band] *= (0.8f+0.015f*n_band);
                    }
                }

                mvr2r(&core_dec_freq[112], &SWB_signal[240+st_offset], 128);
                mvr2r(&core_dec_freq[112], &SWB_signal[368+st_offset], 128);
                mvr2r(&core_dec_freq[176], &SWB_signal[496+st_offset], 64);

                tmp1 = (float)(fabs(SWB_signal[368+st_offset]) + fabs(SWB_signal[369+st_offset])) + EPSILON;
                tmp2 = (float)(fabs(SWB_signal[365+st_offset]) + fabs(SWB_signal[366+st_offset])) + EPSILON;
                pit1 = &SWB_signal[368+st_offset];

                tmp3 = tmp2/tmp1;
                if(tmp3 < 0.3)
                {
                    tmp3 = 0.3f;
                }

                while(tmp3 < 1)
                {
                    *pit1++ *= tmp3;
                    tmp3 += 0.1f;
                }

                pit1 = &SWB_signal[367+st_offset];
                tmp3 = tmp1/tmp2;

                if( tmp3 > 5 )
                {
                    tmp3 = 5;
                    while( tmp3 > 1 )
                    {
                        *pit1-- *= tmp3;
                        tmp3 -= 0.5f;
                    }
                }

                tmp1 = (float)(fabs(SWB_signal[496+st_offset]) + fabs(SWB_signal[497+st_offset])) + EPSILON;
                tmp2 = (float)(fabs(SWB_signal[492+st_offset]) + fabs(SWB_signal[493+st_offset]) + fabs(SWB_signal[494+st_offset]) + fabs(SWB_signal[495+st_offset])) + EPSILON;
                pit1 = &SWB_signal[496+st_offset];

                tmp3 = tmp2/tmp1;
                if( tmp3 < 0.3 )
                {
                    tmp3 = 0.3f;
                }

                while(tmp3 < 1)
                {
                    *pit1++ *= tmp3;
                    tmp3 += 0.1f;
                }

                pit1 = &SWB_signal[495+st_offset];

                tmp3 = tmp1/tmp2;
                tmp3 = 0.5f*tmp3;
                tmp4 = 0.05f*tmp3;

                while( tmp3 > 1 )
                {
                    *pit1-- *= tmp3;
                    tmp3 -= tmp4;
                }

                /* calculate envelope */
                calc_norm_envelop(SWB_signal, envelope, L_swb_norm, SWB_flength, st_offset);
            }
        }

        /* Normalize with envelope */
        if( *frica_flag == 0 && mode != NOISE )
        {
            L = swb_bwe_subband[0]+st_offset;
            inv_L_swb_norm = 1.0f/L_swb_norm;
            weight = (mode != HARMONIC) ? max(min(3.0f*inv_L_swb_norm, 0.5f), 0.2f) : 0.2f;
            weight = 0.4f*weight + 0.6f*(*prev_weight);
            for(n_freq = L; n_freq<swb_bwe_subband[SWB_FENV]+st_offset; n_freq++)
            {
                signum[n_freq] = 1;
                if (SWB_signal[n_freq]<0)
                {
                    signum[n_freq] = -1;
                    SWB_signal[n_freq] *= signum[n_freq];
                }

                SWB_signal[n_freq] = SWB_signal[n_freq] - envelope[n_freq]*inv_L_swb_norm;
                if(SWB_signal[n_freq] > 0)
                {
                    SWB_signal[n_freq] *= (1.2f - weight);
                }
                SWB_signal[n_freq] *= signum[n_freq];
            }

            for (n_freq = L; n_freq<swb_bwe_subband[SWB_FENV]+st_offset; n_freq++)
            {
                SWB_signal[n_freq] /= envelope[n_freq];
            }

            *prev_weight = weight;
        }
        else
        {
            *prev_weight = max(min(3.0f/L_swb_norm, 0.5f), 0.2f);
        }

        if( mode == HARMONIC )
        {
            pit1 = &SWB_signal[swb_bwe_subband[0]+st_offset];
            for(n_band=0; n_band<19; n_band++)
            {
                mean = 0;
                for(n_freq=0; n_freq<16; n_freq++)
                {
                    mean += (float) fabs(*pit1);
                    pit1++;
                }
                mean /= 16;
                pit1 -= 16;
                for(n_freq=0; n_freq<16; n_freq++)
                {
                    if(fabs(*pit1) < mean)
                    {
                        *pit1 *= 0.2f;
                    }
                    pit1++;
                }
            }
        }

        if( mode == HARMONIC )
        {
            L = 2;
        }
        else
        {
            L = 1;
        }

        for(n_band=0; n_band<SWB_FENV; n_band+=L)
        {
            energy = EPSILON;
            for (n_freq = swb_bwe_subband[n_band]+st_offset; n_freq<swb_bwe_subband[n_band+L]+st_offset; n_freq++)
            {
                energy += SWB_signal[n_freq] * SWB_signal[n_freq];
            }

            tmp_ener = (float)sqrt((swb_bwe_subband[n_band+L] - swb_bwe_subband[n_band])/energy);

            for(n_freq = swb_bwe_subband[n_band]+st_offset; n_freq<swb_bwe_subband[n_band+L]+st_offset; n_freq++)
            {
                SWB_signal[n_freq] *= tmp_ener;
            }
        }

        if( *prev_Energy > 1.25f*Energy && Energy > 0 )
        {
            weight = 0.5f*Energy/(*prev_Energy);
        }
        else
        {
            weight = 0.5f;
        }

        wfenv = weight*prev_SWB_fenv[0] + (1-weight)*SWB_fenv[0];
        factor = fenvL;
        factor1 = (wfenv - fenvL) * 0.125f;
        for (n_freq = swb_bwe_subband[0]+st_offset; n_freq < swb_bwe_subband[0]+8+st_offset; n_freq++)
        {
            SWB_signal[n_freq] *= factor;
            factor += factor1;
        }

        for(n_band = 0; n_band < 12; n_band++)
        {
            wfenv = weight*prev_SWB_fenv[n_band+1] + (1-weight)*SWB_fenv[n_band+1];
            factor = SWB_fenv[n_band];
            factor1 = (wfenv - SWB_fenv[n_band]) * smooth_factor[n_band];
            for (; n_freq < swb_bwe_sm_subband[n_band+1]+st_offset; n_freq++)
            {
                SWB_signal[n_freq] *= factor;
                factor += factor1;
            }
        }

        wfenv = weight*prev_SWB_fenv[13] + (1-weight)*SWB_fenv[13];
        factor = SWB_fenv[12];
        factor1 = (wfenv - SWB_fenv[12]) * smooth_factor[12];
        for ( ; n_freq < swb_bwe_sm_subband[13]+st_offset; n_freq++)
        {
            SWB_signal[n_freq] *= factor;
            factor += factor1;
        }

        for(n_band=13; n_band<SWB_FENV; n_band++)
        {
            wfenv = weight*prev_SWB_fenv[n_band] + (1-weight)*SWB_fenv[n_band];
            for ( ; n_freq<swb_bwe_subband[n_band+1]+st_offset; n_freq++)
            {
                SWB_signal[n_freq] *= wfenv;
            }
        }

        for(n_band = 0; n_band < SWB_FENV; n_band++)
        {
            prev_SWB_fenv[n_band] = SWB_fenv[n_band];
        }
    }

    pit1 = &SWB_signal[240+st_offset];
    for(n_freq=0; n_freq<4; n_freq++)
    {
        *(pit1++) *= 0.5f;
    }
    *prev_Energy = Energy;

    return;
}

/*-------------------------------------------------------------------*
 * time_envelop_shaping()
 *
 * time shaping of the SHB signal
 *-------------------------------------------------------------------*/

void time_envelop_shaping(
    float werr[],             /* i/o: SHB synthesis           */
    float SWB_tenv[],         /* i/o: frequency envelope      */
    const short L                 /* i  : frame length            */
)
{
    float *pit;
    float Energy;
    short i, j;
    float tmp_ener;

    pit = werr;
    for(i=0; i<SWB_TENV; i++)
    {
        Energy = EPSILON;
        for(j=0; j<L/4; j++)
        {
            Energy += *pit **pit;
            pit++;
        }
        Energy = (float)sqrt(4*Energy/L);

        if(SWB_tenv[i] < 2 && Energy < SWB_tenv[i])
        {
            SWB_tenv[i] = Energy;
        }

        pit -= L/4;
        tmp_ener = 1.0f / Energy;
        for (j = 0; j < L/4; j++)
        {
            *pit *= SWB_tenv[i] * tmp_ener;
            pit++;
        }
    }

    return;
}

/*-------------------------------------------------------------------*
 * time_reduce_pre_echo()
 *
 * Pre-echo reduction.
 *-------------------------------------------------------------------*/

void time_reduce_pre_echo(
    const float *synth,             /* i  : ACELP core synthesis    */
    float *error,             /* o  : SHB BWE synthesis       */
    float prev_td_energy,     /* o  : last td energy          */
    const short L                   /* i  : subframe length         */
)
{
    short i, j, pos = 0;
    float energy;
    float energyL[4];
    float tmp_ener;
    float *pit;
    float tmpi;

    for(i=0; i<4; i++)
    {
        energyL[i] = 0;
        for(j=0; j<L; j++)
        {
            energyL[i] += synth[L*i+j]*synth[L*i+j];
        }
        energyL[i] = (float)sqrt(energyL[i]/L);
    }

    for(i=0; i<3; i++)
    {
        if(energyL[i+1] > 1.8f*energyL[i] && energyL[i+1] > 50)
        {
            pos = i+1;
            break;
        }
    }

    if (pos > 0)
    {
        if(pos < 3)
        {
            pos ++;
        }

        energy = EPSILON;
        j = L*pos;
        for(i=0; i<j; i++)
        {
            energy += error[i]*error[i];
        }
        energy = (float)sqrt(energy/j);

        if(prev_td_energy < 0.2f*energy)
        {
            prev_td_energy = 0.2f*energy;
        }

        tmp_ener = prev_td_energy / energy;
        for (i = 0; i < j; i++)
        {
            error[i] *= tmp_ener;
        }

        energy = EPSILON;
        for(i=j; i<(j+L); i++)
        {
            energy += error[i]*error[i];
        }
        energy = (float)sqrt(energy/L);

        pit = &error[j];
        tmp_ener = prev_td_energy / energy;
        for (i = 0; i < L; i++)
        {
            tmpi = i/(1.0f*L);
            *pit++ *= ((1.0f - tmpi) * tmp_ener + tmpi);
        }
    }

    return;
}

/*-------------------------------------------------------------------*
 * hq_generic_hf_decoding()
 *
 *-------------------------------------------------------------------*/
void hq_generic_hf_decoding(
    const short HQ_mode,                     /* i  : HQ mode                                     */
    float *coeff_out1,                 /* i/o: BWE input & temporary buffer                */
    const float *hq_generic_fenv,            /* i  : SWB frequency envelopes                     */
    float *coeff_out,                  /* o  : SWB signal in MDCT domain                   */
    const short hq_generic_offset,           /* i  : frequency offset for representing hq swb bwe*/
    short *prev_L_swb_norm,            /* i/o: last normalize length                       */
    const short hq_generic_exc_clas,          /* i  : bwe excitation class                        */
    const short *R
)
{
    short i, n_freq, n_band, L_swb_norm;
    float fenvL, energy, wfenv, factor;
    float envelope[L_FRAME16k];
    float *pit1;
    float tmp1, tmp2, tmp3, tmp4;
    float mean_vector[20];
    short k;
    short nenv;
    short tenv;
    float rn_weight0;
    short blen,nband_lf,sfidx,efidx;
    short bwe_seed;
    short signum[L_FRAME16k];

    if ( hq_generic_offset <= HQ_GENERIC_FOFFSET_24K4 )
    {
        nenv = SWB_FENV;
    }
    else
    {
        nenv = SWB_FENV-2;
    }

    if ( HQ_mode == HQ_GEN_FB )
    {
        tenv = nenv + DIM_FB;
    }
    else
    {
        tenv = nenv;
    }

    fenvL = 0.0f;
    for( n_freq = HQ_GENERIC_ST_FREQ + hq_generic_offset; n_freq<swb_bwe_subband[0] + hq_generic_offset; n_freq++ )
    {
        fenvL += coeff_out1[n_freq] * coeff_out1[n_freq];
    }

    fenvL = (float)sqrt(fenvL/16);

    calc_normal_length( HQ_CORE, coeff_out1, HQ_GEN_SWB, -1, &L_swb_norm, prev_L_swb_norm );

    /* calculate envelope */
    calc_norm_envelop_lf( coeff_out1, envelope, &L_swb_norm, HQ_mode, hq_generic_offset, &sfidx, &efidx);

    blen = 16;
    if (hq_generic_exc_clas == HQ_GENERIC_EXC0 )
    {
        rn_weight0 = 0.8f;
    }
    else if (hq_generic_exc_clas == HQ_GENERIC_EXC1 )
    {
        rn_weight0 = 0.05f;
    }
    else
    {
        rn_weight0 = 0.2f;
    }

    nband_lf = (efidx-sfidx)/blen;
    for ( n_freq = sfidx; n_freq < efidx; n_freq++ )
    {
        if (coeff_out1[n_freq]<0)
        {
            signum[n_freq] = -1;
            coeff_out1[n_freq] *= signum[n_freq];
        }
        else
        {
            signum[n_freq] = 1;
        }
    }

    /* applying whitening */
    for ( n_freq = sfidx; n_freq < efidx; n_freq++ )
    {
        coeff_out1[n_freq] = coeff_out1[n_freq] / envelope[n_freq];
    }

    /* mean vector generation for controlling dynamic range */
    for( k =0 ; k < nband_lf; ++k )
    {
        energy = EPSILON;
        for ( i = k*blen+sfidx; i < (k+1)*blen+sfidx; ++i )
        {
            energy += coeff_out1[i];
        }
        mean_vector[k] = energy / blen;
    }

    /* dynamic range control */
    for( k =0 ; k < nband_lf; ++k )
    {
        for ( i = k*blen+sfidx; i < (k+1)*blen+sfidx; ++i )
        {
            coeff_out1[i] = coeff_out1[i] - (coeff_out1[i]-mean_vector[k])*rn_weight0;
        }
    }

    if (hq_generic_exc_clas == HQ_GENERIC_EXC0)
    {
        /* applying random sign */
        bwe_seed = R[0]*8+R[1]*4+R[2]*2+R[3];
        for ( n_freq = sfidx; n_freq < efidx; n_freq++ )
        {
            coeff_out1[n_freq] = coeff_out1[n_freq]*signum[n_freq]*(own_random(&bwe_seed)>0 ? 1.0f: -1.0f);
        }
    }
    else
    {
        for ( n_freq = sfidx; n_freq < efidx; n_freq++ )
        {
            coeff_out1[n_freq] =coeff_out1[n_freq]*signum[n_freq];
        }
    }

    /* normalizing modified low frequency spectrum */
    for( k =0 ; k < nband_lf; ++k )
    {
        energy = EPSILON;
        for ( i = k*blen+sfidx; i < (k+1)*blen+sfidx; ++i )
        {
            energy += coeff_out1[i] * coeff_out1[i];
        }
        energy = (float)sqrt(energy/blen);

        for ( i = k*blen+sfidx; i < (k+1)*blen+sfidx; ++i )
        {
            coeff_out1[i] /= energy;
        }
    }
    mvr2r(coeff_out1+HQ_GENERIC_OFFSET, &coeff_out[HQ_GENERIC_HIGH0+hq_generic_offset], HQ_GENERIC_LEN0);
    mvr2r(coeff_out1+HQ_GENERIC_OFFSET, &coeff_out[HQ_GENERIC_HIGH1+hq_generic_offset], HQ_GENERIC_LEN0);

    if ( hq_generic_offset <= HQ_GENERIC_FOFFSET_24K4 )
    {
        mvr2r( &coeff_out1[HQ_GENERIC_LOW0], &coeff_out[HQ_GENERIC_HIGH2+hq_generic_offset], HQ_GENERIC_END_FREQ - HQ_GENERIC_HIGH2 );
    }

    if ( HQ_mode == HQ_GEN_FB )
    {
        if ( hq_generic_offset <= HQ_GENERIC_FOFFSET_24K4 )
        {
            mvr2r(coeff_out1+HQ_GENERIC_LOW0 + HQ_GENERIC_END_FREQ - HQ_GENERIC_HIGH2, &coeff_out[fb_bwe_subband[0]], 160);
        }
        else
        {
            mvr2r(coeff_out1+HQ_GENERIC_OFFSET + HQ_GENERIC_LEN0, &coeff_out[fb_bwe_subband[0]], 160);
        }
    }
    tmp1 = EPSILON;
    tmp2 = EPSILON;
    for(i=0; i<5; ++i)
    {
        tmp1 += (float)fabs(coeff_out[HQ_GENERIC_HIGH1+hq_generic_offset+i]);
        tmp2 += (float)fabs(coeff_out[HQ_GENERIC_HIGH1-2+hq_generic_offset-i]);
    }

    pit1 = &coeff_out[HQ_GENERIC_HIGH1+hq_generic_offset];
    tmp3 = tmp2/tmp1;
    if( tmp3 < 0.3f )
    {
        tmp3 = 0.3f;
    }

    while( tmp3 < 1 )
    {
        *pit1++ *= tmp3;
        tmp3 += 0.1f;
    }

    pit1 = &coeff_out[HQ_GENERIC_HIGH1-1+hq_generic_offset];
    tmp3 = tmp1/tmp2;
    if( tmp3 > 5 )
    {
        tmp3 = 5;
        while( tmp3 > 1 )
        {
            *pit1-- *= tmp3;
            tmp3 -= 0.5f;
        }
    }

    if ( hq_generic_offset <= HQ_GENERIC_FOFFSET_24K4 )
    {
        tmp1 = (float)(fabs(coeff_out[HQ_GENERIC_HIGH2+hq_generic_offset]) + fabs(coeff_out[HQ_GENERIC_HIGH2+1+hq_generic_offset])) + EPSILON;
        tmp2 = (float)(fabs(coeff_out[HQ_GENERIC_HIGH2-4+hq_generic_offset]) + fabs(coeff_out[HQ_GENERIC_HIGH2-3+hq_generic_offset]) +
                       fabs(coeff_out[HQ_GENERIC_HIGH2-2+hq_generic_offset]) + fabs(coeff_out[HQ_GENERIC_HIGH2-1+hq_generic_offset])) + EPSILON;

        pit1 = &coeff_out[HQ_GENERIC_HIGH2+hq_generic_offset];
        tmp3 = tmp2/tmp1;
        if( tmp3 < 0.3f )
        {
            tmp3 = 0.3f;
        }

        while( tmp3 < 1 )
        {
            *pit1++ *= tmp3;
            tmp3 += 0.1f;
        }

        pit1 = &coeff_out[HQ_GENERIC_HIGH2-1+hq_generic_offset];
        tmp3 = tmp1/tmp2;
        tmp3 = 0.5f*tmp3;
        tmp4 = 0.05f*tmp3;
        while(tmp3 > 1)
        {
            *pit1-- *= tmp3;
            tmp3 -= tmp4;
        }
    }

    wfenv = hq_generic_fenv[0];
    for (n_freq = swb_bwe_subband[0]+hq_generic_offset,i=0; n_freq<swb_bwe_subband[0]+hq_generic_offset+8; n_freq++,i++)
    {
        factor = i*0.125f;
        coeff_out[n_freq] *= ((1-factor)*fenvL + factor*wfenv);
    }

    for(n_band=0; n_band<nenv-2; n_band++)
    {
        wfenv = hq_generic_fenv[n_band+1];
        for ( i=0; n_freq<swb_bwe_sm_subband[n_band+1]+hq_generic_offset; n_freq++, i++)
        {
            factor = i*smooth_factor[n_band];
            coeff_out[n_freq] *= ((1-factor)*hq_generic_fenv[n_band] + factor*wfenv);
        }
    }

    wfenv = hq_generic_fenv[nenv-1];
    for ( i=0; n_freq<swb_bwe_sm_subband[nenv-1]+hq_generic_offset; n_freq++, i++)
    {
        factor = i*smooth_factor[nenv-2];

        coeff_out[n_freq] *= ((1-factor)*hq_generic_fenv[nenv-2] + factor*wfenv);
    }

    if ( HQ_mode == HQ_GEN_SWB )
    {
        for(n_band=nenv-1; n_band<nenv; ++n_band)
        {
            wfenv = hq_generic_fenv[n_band];
            for ( ; n_freq<swb_bwe_subband[n_band+1]+hq_generic_offset; n_freq++)
            {
                coeff_out[n_freq] *= wfenv;
            }
        }
    }
    else
    {
        if ( hq_generic_fenv[nenv - 1] - hq_generic_fenv[nenv] > 15.f ||  hq_generic_fenv[nenv] < 5.f)
        {
            wfenv = hq_generic_fenv[nenv-1];
            for ( i=0; n_freq<fb_bwe_subband[0]; n_freq++, i++)
            {
                coeff_out[n_freq] *= wfenv;
            }

            for(n_band=0; n_band<DIM_FB ; n_band++)
            {
                wfenv = hq_generic_fenv[n_band + nenv];
                for ( i=0; n_freq<fb_bwe_subband[n_band+1]; n_freq++, i++)
                {
                    coeff_out[n_freq] *= wfenv;
                }
            }
        }
        else
        {
            for(n_band=0; n_band<DIM_FB ; n_band++)
            {
                wfenv = hq_generic_fenv[n_band + nenv - 1];

                for ( i=0; n_freq<fb_bwe_sm_subband[n_band]; n_freq++, i++)
                {
                    factor = i* fb_smooth_factor[n_band];
                    coeff_out[n_freq] *= ((1-factor)*hq_generic_fenv[n_band+nenv] + factor*wfenv);
                }
            }

            wfenv = hq_generic_fenv[tenv-1];

            for ( ; n_freq<fb_bwe_subband[DIM_FB]; n_freq++)
            {
                coeff_out[n_freq] *= wfenv;
            }
        }
    }
    return;
}


/*-------------------------------------------------------------------*
 * save_old_syn()
 *
 * Save and delay the ACELP core synthesis signal by
 * DELAY_FD_BWE_ENC_xxkx to be used by SWB BWE
 *-------------------------------------------------------------------*/

void save_old_syn(
    const short L_frame,        /* i  : frame length                */
    const float syn[],          /* i  : ACELP synthesis             */
    float old_syn[],      /* o  : old synthesis buffer        */
    float old_syn_mem[],  /* i/o: old synthesis buffer memory */
    const float preemph_fac,    /* i  : preemphasis factor          */
    float *mem_deemph     /* i/o: deemphasis filter memory    */
)
{
    short tmps;

    if( L_frame == L_FRAME )
    {
        tmps = NS2SA(12800, DELAY_FD_BWE_ENC_12k8_NS);
    }
    else
    {
        tmps = NS2SA(16000, DELAY_FD_BWE_ENC_16k_NS);
    }

    mvr2r( old_syn_mem, old_syn, tmps );
    mvr2r( syn, old_syn + tmps, L_frame - tmps );
    mvr2r( syn + L_frame - tmps, old_syn_mem, tmps );

    deemph( old_syn, preemph_fac, L_frame, mem_deemph );

    return;
}

