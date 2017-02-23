/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "prot.h"
#include "rom_enc.h"



static int comvad_hangover(
    float lt_snr_org,                 /*(i)original long time SNR*/
    float snr,                        /*(i) frequency domain SNR */
    float lf_snr,                     /*(i) long time frequency domain SNR calculated by l_speech_snr and l_silence_snr*/
    float snr_flux,                   /*(i) average tsnr*/
    int   bw_index,                   /*(i) band width index*/
    int   vad_flag,
    int   pre_res_hang_num,           /*(i) residual number of previous  hangover */
    int   continuous_speech_num2,     /*(i) number of continuous speech frames*/
    int   noisy_type                  /*(i) noisy type*/
);



int comvad_decision(T_CldfbVadState *st,
                    float snr,                       /*(i) frequency domain SNR */
                    float tsnr,                       /*(i) time domain SNR */
                    float snr_flux,                /*(i) average tsnr of several frames*/
                    float lt_snr,                  /*(i)long time SNR calculated by fg_energy and bg_energy*/
                    float lt_snr_org,               /*(i)original long time SNR*/
                    float lf_snr,                  /*(i) long time frequency domain
                                                   SNR calculated by l_speech_snr and l_silence_snr*/
                    float frame_energy,            /*(i) current frame energy */
                    int   music_backgound_f,       /*(i) background music flag*/
                    short *cldfb_addition,
                    short vada_flag
                   )
{
    int   speech_flag = st->speech_flag;
    int   fg_energy_count = st->fg_energy_count;
    int   bg_energy_count = st->bg_energy_count;
    float fg_energy = st->fg_energy;
    float bg_energy = st->bg_energy;
    int   l_speech_snr_count=st->l_speech_snr_count;
    int   vad_flag;
    float snr_thresh = 0.2f;
    float *ltd_stable_rate = st->ltd_stable_rate;
    float *sp_center = st->sp_center;
    int   frameloop = st->frameloop;
    int   bw_index = st->bw_index;
    int noisy_type = UNKNOWN_NOISE;
    short vadb_flag = 0;

    /*
    * ls_snr_org
    * -------NB---------------WB----------------SWB---------
    * l16  [13.5 17]       [13.5    16]       [13.1  15.2]
    * 126  [12   14]       [10.7  12.7]       [10.3  12.1]
    * l36  [9.5  12]       [8.5   10.5]       [7.8    9.2]
    * 15dB [3.5 4.5]
    * 20dB [4.5 5.5]
    */

    /*
    * lt_snr_org
    * -------NB---------------WB----------------SWB---------
    * l16  [-- --]       [--    --]       [--  --]
    * 126  [--   --]       [--  -]       [--  --]
    * l36  [--  --]       [--   --]       [-    --]
    * 15dB [1.2 1.9]
    * 20dB [2.1 2.7]                          [1.8  2.3]
    */

    if(st->lf_snr_smooth  > LS_MIN_SELENCE_SNR[bw_index - CLDFBVAD_NB_ID] && lt_snr_org > LT_MIN_SILENCE_SNR[bw_index - CLDFBVAD_NB_ID])

    {
        noisy_type = SILENCE;
    }
    snr_thresh = construct_snr_thresh( sp_center, snr_flux, lt_snr, lf_snr,
                                       st->continuous_speech_num, st->continuous_noise_num, st->fg_energy_est_start, bw_index);

    if(snr>(snr_thresh))
        vad_flag =1;
    else
        vad_flag =0;
    if(tsnr> 4.0)
    {
        vad_flag =1;
    }


    if(frameloop>25)
    {
        if(vad_flag == 1&&st->fg_energy_est_start==1)
        {
            if(fg_energy_count==512)
            {
                fg_energy = fg_energy*0.75f;
                fg_energy_count = 384;
            }
            if((frame_energy*bg_energy_count)>6*bg_energy)
            {
                fg_energy = fg_energy + frame_energy;
                fg_energy_count = fg_energy_count + 1;
            }
        }
    }
    if(music_backgound_f)
    {
        vad_flag =1;
    }
    if(vad_flag==1)
    {

        if(snr > st->l_silence_snr/st->l_silence_snr_count + 1.5)
        {
            if(l_speech_snr_count==512)
            {
                st->l_speech_snr = st->l_speech_snr*0.75f;
                l_speech_snr_count = 384;

                st->l_speech_snr += snr;
                l_speech_snr_count++;
            }
            else
            {
                st->l_speech_snr += snr;
                l_speech_snr_count++;
            }
        }
    }
    if(bw_index == CLDFBVAD_NB_ID)
    {
        if(snr_flux > 1.9+ lt_snr*0.28 )
        {
            vad_flag = 1;
        }
        if((snr_flux > 1.5)&&sp_center[3]>1.6&& lt_snr_org<3.5)
        {
            vad_flag = 1;
        }
        if((snr_flux > 1.2)&&sp_center[3]>1.9&& lt_snr_org<3.5)
        {
            vad_flag = 1;
        }
        if((snr_flux > 1.00f)&&sp_center[3]>3.2&& lt_snr_org<3.5)
        {
            vad_flag = 1;
        }
    }
    if(bw_index == CLDFBVAD_WB_ID)
    {
        if((snr_flux > 2.1+lt_snr*0.24) )
        {
            vad_flag = 1;
        }
        if((snr_flux > 1.6)&&sp_center[3]>2.5&& lt_snr_org<3.5)
        {
            vad_flag = 1;
        }
        if((snr_flux > 1.2)&&sp_center[3]>2.8&& lt_snr_org<3.5)
        {
            vad_flag = 1;
        }
        if((snr_flux > 1.0)&&sp_center[3]>4.5&& lt_snr_org<3.5)
        {
            vad_flag = 1;
        }
    }
    if(bw_index == CLDFBVAD_SWB_ID)
    {

        if((snr_flux > 2.1+ lt_snr*0.32) )
        {
            vad_flag = 1;
        }
        if((snr_flux > 1.68)&&sp_center[3]>2.76&& lt_snr_org<3.5)
        {
            vad_flag = 1;
        }
        if((snr_flux > 1.24)&&sp_center[3]>2.92&& lt_snr_org<3.5)
        {
            vad_flag = 1;
        }
        if((snr_flux > 1.10f)&&sp_center[3]>4.6&& lt_snr_org<3.5)
        {
            vad_flag = 1;
        }
    }
    if(st->fg_energy_est_start==0)
    {
        if(ltd_stable_rate[0]>0.08 && vad_flag == 1 && frame_energy>50)
        {
            st->fg_energy_est_start=1;
        }
    }

    /************************************************************************/
    /*   hangover                                                           */
    /************************************************************************/
    speech_flag = comvad_hangover(lt_snr_org, snr, lf_snr, snr_flux, bw_index, vad_flag, speech_flag,
                                  st->continuous_speech_num2, noisy_type);

    if(vad_flag==0&&speech_flag>0)
    {
        speech_flag--;
        vad_flag = 1;
    }

    vadb_flag = vad_flag;

    if(bw_index == CLDFBVAD_SWB_ID)
    {

        if(SILENCE == noisy_type && snr > 0.2 && vad_flag==0)
        {
            vad_flag = vada_flag;
        }
        else if(st->lf_snr_smooth < 10.5 || SILENCE!=noisy_type)
        {
            if((snr_flux > 2.0)
                    || (st->continuous_speech_num2 > 40 && (snr_flux > 1.8))
                    || music_backgound_f == 1)
            {
                vad_flag = vada_flag | vadb_flag;
            }
            /*only use for silence*/
            else if(noisy_type == SILENCE)
            {
                vad_flag = vada_flag;
            }
        }
    }
    else if(bw_index == CLDFBVAD_WB_ID)
    {

        if(SILENCE == noisy_type && snr > 0.2 && vad_flag==0)
        {
            vad_flag = vada_flag;
        }
        else
        {
            if(st->lf_snr_smooth < 10.5 || SILENCE!=noisy_type)
            {
                if(snr_flux > 2.2
                        || (st->continuous_speech_num2 > 40 && (snr_flux > 1.5))
                        || music_backgound_f == 1)
                {
                    vad_flag = vada_flag|vadb_flag;
                }
                else if(SILENCE == noisy_type)
                {
                    vad_flag = vada_flag;
                }

            }

        }

    }
    else
    {
        if(noisy_type == SILENCE  )
        {
            if(st->lf_snr_smooth > 12.5 && music_backgound_f == 0)
            {
                vad_flag = vada_flag;
            }
        }
        else
        {
            if((snr_flux > 2.0)
                    || (st->continuous_speech_num2 > 30 && (snr_flux > 1.5))
                    || music_backgound_f == 1)
            {
                vad_flag = vada_flag | vadb_flag;
            }
        }
    }
    if(vad_flag==0 )
    {
        if(st->l_silence_snr_count==512)
        {
            st->l_silence_snr = st->l_silence_snr*0.75f;
            st->l_silence_snr_count = 384;

            st->l_silence_snr += snr;
            st->l_silence_snr_count++;
        }
        else if(snr<0.8)
        {
            st->l_silence_snr += snr;
            st->l_silence_snr_count++;
        }
    }
    if((vad_flag + vada_flag) == 0)

    {
        if(bg_energy_count==512)
        {
            bg_energy = bg_energy*0.75f;
            bg_energy_count = 384;
        }

        if(tsnr<1.0)
        {
            bg_energy = bg_energy + frame_energy;
            bg_energy_count = bg_energy_count +1;
        }
    }

    st->lt_snr_org = lt_snr_org;
    st->speech_flag = speech_flag;

    st->fg_energy_count = fg_energy_count;
    st->bg_energy_count = bg_energy_count;
    st->fg_energy = fg_energy;
    st->bg_energy = bg_energy;
    st->l_speech_snr_count = l_speech_snr_count;


    st->vad_flag_for_bk_update = vad_flag;
    if(st->update_count < 12 && vadb_flag==1)
    {
        st->warm_hang_num = max(20, speech_flag);
    }


    if(vad_flag==0&&st->warm_hang_num>0)
    {
        st->warm_hang_num--;
        vad_flag = 1;
    }

    if(noisy_type == SILENCE
            && bw_index != CLDFBVAD_NB_ID)

    {
        *cldfb_addition = 2;
    }
    else
    {
        *cldfb_addition = 0;

        if(bw_index == CLDFBVAD_WB_ID)
        {
            *cldfb_addition = 3;
        }
        if(bw_index == CLDFBVAD_SWB_ID)
        {
            *cldfb_addition = 1;
        }
        if(st->bw_index == CLDFBVAD_NB_ID)
        {
            *cldfb_addition = 1;
        }
    }

    return vad_flag;
}

float construct_snr_thresh(float sp_center[],                 /*(i) spectral center*/
                           float snr_flux,                    /*(i) snr flux*/
                           float lt_snr,                      /*(i) long time time domain snr*/
                           float lf_snr,                      /*(i) long time frequency domain snr*/
                           int continuous_speech_num,         /*(i) number of continuous speech frames*/
                           int continuous_noise_num,          /*(i) number of continuous noise frames*/
                           int fg_energy_est_start,           /*(i) whether if estimated energy*/
                           int bw_index                       /*(i) band width index*/
                          )
{
    float test_l_snr=0.f;
    float snr_delta;
    float snr_thresh;
    float bw_snr;


    snr_delta = COMVAD_INIT_SNR_DELTA[bw_index];
    bw_snr = lt_snr;

    if(bw_index == CLDFBVAD_SWB_ID)
    {

        test_l_snr  = lt_snr;
        test_l_snr = test_l_snr*1.0f;

        if(sp_center[3]>2.80f)
        {
            snr_delta = snr_delta + 0.00f;
        }
        else if(sp_center[2]>2.6)
        {
            snr_delta = snr_delta + 0.03f;
        }
        else if(sp_center[2]>1.6)
        {
            snr_delta = snr_delta + 0.05f;
        }
        else if(sp_center[3]>1.4)
        {
            snr_delta = snr_delta + 0.10f;
        }
        else
        {
            snr_delta = snr_delta + 0.40f;
        }

        if(continuous_speech_num > 8&&fg_energy_est_start==1)
        {
            snr_delta  = snr_delta - 0.2f;
        }
        else if(continuous_noise_num > 12&&(snr_flux>0.6+lf_snr*0.1))
        {
            snr_delta = snr_delta + 0.1f;
        }
        else if(continuous_noise_num > 24)
        {
            snr_delta = snr_delta + 0.2f;
        }
        else if((continuous_noise_num > 4))
        {
            snr_delta = snr_delta + 0.1f;
        }

    }
    else if(bw_index == CLDFBVAD_WB_ID)
    {

        test_l_snr  = lt_snr;
        if(sp_center[3]>2.80f)
        {
            snr_delta = snr_delta + 0.00f;
        }
        else if(sp_center[2]>2.6)
        {
            snr_delta = snr_delta + 0.03f;
        }
        else if(sp_center[2]>1.6)
        {
            snr_delta = snr_delta + 0.05f;
        }
        else if(sp_center[3]>1.4)
        {
            snr_delta = snr_delta + 0.10f;
        }
        else
        {
            snr_delta = snr_delta + 0.30f;
        }

        if(continuous_speech_num > 8&&fg_energy_est_start==1)
        {
            snr_delta  = snr_delta - 0.1f;
        }
        else if(continuous_noise_num > 12&&(snr_flux>0.6+bw_snr*0.1))
        {
            snr_delta = snr_delta + 0.1f;
        }
        else if(continuous_noise_num > 24)
        {
            snr_delta   = snr_delta + 0.2f;
        }
        else if((continuous_noise_num > 4))
        {
            snr_delta = snr_delta + 0.1f;
        }
    }
    else if(bw_index == CLDFBVAD_NB_ID)
    {

        test_l_snr  = lt_snr;

        if(sp_center[3]>3.0)
        {
            snr_delta = snr_delta + 0.00f;
        }
        else if(sp_center[2]>2.6)
        {
            snr_delta = snr_delta + 0.02f;
        }
        else if(sp_center[2]>1.6)
        {
            snr_delta = snr_delta + 0.04f;
        }
        else if(sp_center[2]>1.46)
        {
            snr_delta = snr_delta + 0.10f;
        }
        else
        {
            snr_delta = snr_delta + 0.18f;
        }

        if(continuous_speech_num > 80&&fg_energy_est_start==1&&(sp_center[0]>1.4))
        {
            snr_delta  = snr_delta - 0.32f;
        }
        else if(continuous_speech_num > 8&&fg_energy_est_start==1&&(snr_flux>0.2+lf_snr*0.1))
        {
            snr_delta  = snr_delta - 0.1f;
        }
        else if(continuous_noise_num > 12&&(snr_flux>0.6+lf_snr*0.1))
        {
            snr_delta = snr_delta + 0.1f;
        }
        else if(continuous_noise_num > 24)
        {

            snr_delta = snr_delta + 0.2f;
        }
    }
    else
    {
        snr_delta = 1.0f;
    }

    snr_thresh  = snr_delta + test_l_snr  ;


    return snr_thresh;
}

static int comvad_hangover(float lt_snr_org,              /*(i)original long time SNR*/
                           float snr,                     /*(i) frequency domain SNR */
                           float lf_snr,                  /*(i) long time frequency domain
                                                   SNR calculated by l_speech_snr and l_silence_snr*/
                           float snr_flux,                /*(i) average tsnr*/
                           int   bw_index,                /*(i) band width index*/
                           int   vad_flag,
                           int   pre_res_hang_num,        /*(i) residual number of previous  hangover */
                           int   continuous_speech_num2,  /*(i) number of continuous speech frames*/
                           int   noisy_type               /*(i) noisy type*/
                          )

{
    int speech_flag = pre_res_hang_num;


    if(bw_index == CLDFBVAD_SWB_ID)
    {
        if(vad_flag)
        {
            if(lt_snr_org > 3.5f)
                speech_flag = 3;
            else
                speech_flag = 4;
            if((continuous_speech_num2 < 8)&& (lt_snr_org < 4.0f))
            {
                speech_flag = 8 - continuous_speech_num2;
            }
            else if((snr_flux > 0.8 )&&(continuous_speech_num2 > 24))
            {
                if(lt_snr_org > 3.6f)
                {
                    speech_flag = 3;
                }
                else if(lt_snr_org > 2.6f)
                {
                    speech_flag = 3;
                }
                else if(lt_snr_org > 1.6f)
                {
                    speech_flag = 4;
                }
                else
                {
                    speech_flag = 5;
                }
                speech_flag = speech_flag - 1;
            }
            if(continuous_speech_num2 < 120)
            {
                if(snr>1.5)
                {
                    speech_flag = 9;
                }
                else if(snr>1.0&&speech_flag<7)
                {
                    speech_flag = 7;
                }
                else if(speech_flag<3)
                {
                    speech_flag = 3;
                }
                if(speech_flag > 3)
                {
                    speech_flag -= 2;
                }
            }
            else
            {
                if(lt_snr_org > 3.6f)
                {
                    speech_flag = 1;
                }
                else if(lt_snr_org > 3.0f)
                {
                    speech_flag = 2;
                }
                else if(lt_snr_org > 2.5f)
                {
                    speech_flag = 3;
                }
                else if(lt_snr_org > 2.0f)
                {
                    speech_flag = 3;
                }
                else if(lt_snr_org > 1.5f)
                {
                    speech_flag = 4;
                }
                else
                {
                    speech_flag = 5;
                }
            }
            if(noisy_type==SILENCE)
            {
                speech_flag = 6;
            }
        }
    }
    else if(bw_index == CLDFBVAD_WB_ID)
    {
        if(vad_flag)
        {
            if(lt_snr_org > 3.5f)
                speech_flag = 1;
            else
                speech_flag = 2;
            if((continuous_speech_num2 < 8)&& (lt_snr_org < 4.0f))
            {
                speech_flag = 8 - continuous_speech_num2;
            }
            else if((snr_flux > 0.9 )&&(continuous_speech_num2 > 50))
            {
                if(lt_snr_org > 3.6f)
                {
                    speech_flag = 1;
                }
                else if(lt_snr_org > 2.6f)
                {
                    speech_flag = 5;
                }
                else if(lt_snr_org > 1.6f)
                {
                    speech_flag = 6;
                }
                else
                {
                    speech_flag = 7;
                }
                if(speech_flag > 1)
                {
                    speech_flag = speech_flag - 1;
                }

            }
            if(continuous_speech_num2 < 120)
            {
                if(snr>1.5)
                {
                    speech_flag = 6;
                }
                else if(snr>1.0&&speech_flag<5)
                {
                    speech_flag = 5;
                }
                else if(snr>0.8&&lt_snr_org < 2 &&speech_flag<4)
                {
                    speech_flag = 4;
                }
                else if(speech_flag<3)
                {
                    speech_flag = 3;
                }
            }
            else
            {
                if(lt_snr_org > 3.6f)
                {
                    speech_flag = 1;
                }
                else if(lt_snr_org > 3.0f)
                {
                    speech_flag = 2;
                }
                else if(lt_snr_org > 2.5f)
                {
                    speech_flag = 2;
                }
                else if(lt_snr_org > 2.0f)
                {
                    speech_flag = 3;
                }
                else
                {
                    speech_flag = 3;
                }
            }

            if(noisy_type==SILENCE)
            {
                speech_flag = 6;
            }
        }
    }
    else
    {
        if(vad_flag)
        {
            if(lt_snr_org > 3.5f)
                speech_flag = 3;
            else
                speech_flag = 4;
            if((continuous_speech_num2 < 8)&& (lt_snr_org < 4.0f))
            {
                speech_flag = 8 - continuous_speech_num2;
            }
            else if((snr_flux > 0.8 + lf_snr*0.1)&&(continuous_speech_num2 > 24))
            {
                if(lt_snr_org > 3.6f)
                {
                    speech_flag = 3;
                }
                else if(lt_snr_org > 2.6f)
                {
                    speech_flag = 8;
                }
                else if(lt_snr_org > 1.2f)
                {
                    speech_flag = 10;
                }
                else
                {
                    speech_flag = 12;
                }
                if(speech_flag > 2)
                {
                    speech_flag = speech_flag - 2;
                }

            }
            if(continuous_speech_num2 < 120)
            {
                if(snr>1.5)
                {
                    speech_flag = 10;
                }
                else if(snr>1.0&&speech_flag<7)
                {
                    speech_flag = 7;
                }
                else if(speech_flag<3&&continuous_speech_num2>12)
                {
                    speech_flag = 3;
                }
            }
            else
            {
                if(lt_snr_org > 3.6f)
                {
                    speech_flag = 2;
                }
                else if(lt_snr_org > 3.0f)
                {
                    speech_flag = 2;
                }
                else if(lt_snr_org > 2.5f)
                {
                    speech_flag = 3;
                }
                else if(lt_snr_org > 2.0f)
                {
                    speech_flag = 3;
                }
                else if(lt_snr_org > 1.5f)
                {
                    speech_flag = 4;
                }
                else
                {
                    speech_flag = 4;
                }
            }
            if(noisy_type==SILENCE)
            {
                speech_flag = 2;
            }
        }
    }

    if(vad_flag == 1)
    {
        if(noisy_type != SILENCE)
        {
            speech_flag--;
        }
        else
        {
            speech_flag = speech_flag - 3;
        }

        if(speech_flag < 0)
        {
            speech_flag = 0;
        }
    }

    return speech_flag;
}


