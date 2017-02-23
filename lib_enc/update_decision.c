/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "prot.h"
#include <stdio.h>
#include <math.h>



void bg_music_decision(T_CldfbVadState *st,
                       int   *music_backgound_f,   /*(i) background music flag*/
                       float frame_energy          /*(i) current frame energy 1*/
                      )
{
    int music_background_frame=0;
    float *sp_center = st->sp_center;
    float *ltd_stable_rate = st->ltd_stable_rate;
    float *sSFM = st->sfm;
    float *f_tonality_rate = st->f_tonality_rate;


    if((f_tonality_rate[1]>0.60)||(f_tonality_rate[0]>0.86))
    {
        if(ltd_stable_rate[0]<0.072&&sp_center[0]>1.2&&(sSFM[0]<0.76||sSFM[1]<0.88||sSFM[2]<0.96))
        {
            music_background_frame = 1;
        }
    }
    if(music_background_frame&&(4.6*st->fg_energy_count*frame_energy>st->fg_energy)&&(st->fg_energy_est_start==1))
    {
        st->music_background_rate  = st->music_background_rate*0.975f + 0.025f;
    }
    else if(music_background_frame)
    {
        st->music_background_rate  = st->music_background_rate*0.998f + 0.002f;
    }
    else
    {
        st->music_background_rate  = st->music_background_rate*0.997f;
    }
    if(st->music_background_rate>0.5)
    {
        *music_backgound_f = 1;
    }
    else
    {
        *music_backgound_f = 0;
    }
}

int update_decision(T_CldfbVadState *st,
                    float snr,                 /*(i) frequency domain SNR */
                    float tsnr,                /*(i) time domain SNR */
                    float frame_energy,        /*(i) current frame energy*/
                    float high_eng,            /*(i) current frame high frequency energy*/
                    int vad_flag,
                    int music_backgound_f      /*(i) background music flag*/
                   )
{
    float g_high_eng_sacle = 0.0f;
    float sp_center3_diff;
    int update_flag = 1;
    int tonality_flag = 0;
    int frameloop = st->frameloop;
    int bw = st->bw_index;
    float *sp_center = st->sp_center;
    float *ltd_stable_rate = st->ltd_stable_rate;
    float *sSFM = st->sfm;
    float *f_tonality_rate = st->f_tonality_rate;
    float tmpout = frame_energy - 4*st->frame_energy_smooth;


    g_high_eng_sacle = high_eng/(st->lt_bg_highf_eng + FLT_MIN);
    sp_center3_diff = sp_center[3] - st->lt_noise_sp_center3;
    if(frameloop>50)
    {

        if(ltd_stable_rate[0]>0.12)
        {
            update_flag=0;
        }

        if((bw == CLDFBVAD_WB_ID || bw == CLDFBVAD_SWB_ID) &&  st->frame_energy_smooth < 4*frame_energy)
        {
            if(g_high_eng_sacle > 3.0f && (sp_center3_diff > 0.4))
            {
                update_flag = 0;
            }
            if((sp_center[3]>2.8f)&&(ltd_stable_rate[0]>0.02f))
            {
                update_flag = 0;
            }
        }
    }



    if((f_tonality_rate[1] > 0.50) && (ltd_stable_rate[0]>0.1))
    {
        update_flag = 0;
    }
    if(sSFM[1] < 0.92 && sSFM[0] < 0.92 && sSFM[2] < 0.92)
    {
        update_flag=0;
    }
    if(sSFM[0] < 0.80 || sSFM[1] < 0.78||  sSFM[2] < 0.80)
    {
        update_flag=0;
    }

    if(frame_energy > 32*st->frame_energy_smooth)
    {
        update_flag = 0;
    }
    if((4.6*st->fg_energy_count*frame_energy>st->fg_energy)&&(st->fg_energy_est_start==1)&&(frame_energy>3))
    {
        update_flag = 0;
    }
    if((f_tonality_rate[1]>0.60)||(f_tonality_rate[0]>0.86))
    {
        update_flag   = 0;
        tonality_flag = 1;
    }


    if(tonality_flag)
    {
        st->tonality_rate3 = st->tonality_rate3*0.983f + 0.017f;
    }
    else
    {
        st->tonality_rate3 = st->tonality_rate3*0.983f;
    }

    if(st->tonality_rate3>0.5)
    {
        update_flag = 0;
    }


    if((sp_center[0] > 4.0f) && ltd_stable_rate[0]>0.04)
    {
        update_flag = 0;
    }
    if((f_tonality_rate[1] > 0.46) && ((sSFM[1] > 0.93)||(ltd_stable_rate[0]>0.09)))
    {
        update_flag = 0;

    }
    if((sSFM[1] < 0.93 && sSFM[0] < 0.92 && sSFM[2] < 0.97) && (f_tonality_rate[1] > 0.5))
    {
        update_flag = 0;
    }
    if((f_tonality_rate[1] > 0.43)&&(sSFM[0] < 0.95)&&(sp_center[1] > 1.94f))
    {
        update_flag = 0;
    }

    if(update_flag)
    {
        if(st->update_count < 1000)
        {
            st->update_count = st->update_count + 1;
        }

    }


    if(update_flag)
    {
        st->lt_noise_sp_center3 = 0.9f*st->lt_noise_sp_center3 + 0.1f*sp_center[3];
    }
    if((tmpout>0)&&(frameloop<100)&&(f_tonality_rate[1]<0.56)&&((sp_center[0]<1.36)||ltd_stable_rate[0]<0.03))
    {
        update_flag = 1;
    }
    if(snr<0.3 && tmpout < 0
            &&tsnr<1.2&&vad_flag==0&&f_tonality_rate[1]<0.5
            &&(music_backgound_f==0)&&ltd_stable_rate[3]<0.1)

    {
        update_flag = 1;
    }
    if(vad_flag && (snr > 1.0) && bw == CLDFBVAD_SWB_ID && tmpout > 0)
    {
        update_flag = 0;
    }

    if(vad_flag &&  (snr > 1.5) && bw != CLDFBVAD_SWB_ID && tmpout > 0)
    {
        update_flag = 0;
    }


    if(update_flag == 0)
    {
        st->update_num_with_snr = 0;
    }
    else
    {

        if(vad_flag && (snr > 3.0)  && st->update_num_with_snr < 10)
        {
            update_flag = 0;
            st->update_num_with_snr++;
        }
    }



    if(vad_flag==0||update_flag == 1)
    {
        float tmpp =(float)fabs(st->sp_center[2]-st->lt_noise_sp_center0);
        if(tmpp>2.5)
        {
            tmpp = 2.5f;
        }
        st->lt_noise_sp_center_diff_sum +=tmpp;
        st->lt_noise_sp_center_diff_counter++;
        if(st->lt_noise_sp_center_diff_counter==128)
        {
            st->lt_noise_sp_center_diff_sum = st->lt_noise_sp_center_diff_sum*0.75f;
            st->lt_noise_sp_center_diff_counter = 96;
        }
        if(fabs(sp_center[0]-st->lt_noise_sp_center0) >2.4 )
        {
            st->lt_noise_sp_center0 = 0.996f*st->lt_noise_sp_center0 + 0.004f*sp_center[0];
        }
        else if(fabs(sp_center[0]-st->lt_noise_sp_center0) >1.0 )
        {
            st->lt_noise_sp_center0 = 0.99f*st->lt_noise_sp_center0 + 0.01f*sp_center[0];
        }
        else
        {
            st->lt_noise_sp_center0 = 0.96f*st->lt_noise_sp_center0 + 0.04f*sp_center[0];
        }
    }

    if((fabs(sp_center[2]-st->lt_noise_sp_center0) > (6*(st->lt_noise_sp_center_diff_sum/st->lt_noise_sp_center_diff_counter)+0.3))&&st->frameloop>200)
    {
        update_flag = 0;
    }


    return update_flag;
}

