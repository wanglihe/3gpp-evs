/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "prot.h"
#include "cnst.h"
#include "rom_enc.h"



void background_update(T_CldfbVadState *st,
                       float frame_energy,		   /*(i) current frame energy*/
                       int   update_flag,          /*(i) current frame update flag*/
                       int   music_backgound_f     /*(i) background music flag*/
                       ,float snr
                      )
{
    int i;
    int SNR_sb_num;
    float *sb_bg_energy = st->sb_bg_energy;
    float *frame_sb_energy = st->frame_sb_energy;
    float *f_tonality_rate = st->f_tonality_rate;
    float *ltd_stable_rate = st->ltd_stable_rate;
    int   frameloop = st->frameloop;

    float t_bg_energy	= st->t_bg_energy;


    SNR_sb_num = ENERGY_BAND_NUM[st->bw_index - CLDFBVAD_NB_ID];


    frame_energy = frame_energy + 0.0001f;


    if((frameloop<60)&&(frameloop>5)&&(f_tonality_rate[0]<0.56)&&
            (f_tonality_rate[1]<0.5)&&ltd_stable_rate[1]<0.06&&snr<2.5f)
    {

        st->t_bg_energy_sum  += frame_energy;
        st->tbg_energy_count++;


        for(i=0; i<SNR_sb_num; i++)
        {
            sb_bg_energy[i] = sb_bg_energy[i]*0.90f + frame_sb_energy[i]*0.1f;
        }
    }
    if(update_flag==1&&frameloop>2&&music_backgound_f==0)
    {
        if(st->bg_update_count<16)
        {
            st->t_bg_energy_sum  += frame_energy;
            st->tbg_energy_count++;
            for(i=0; i<SNR_sb_num; i++)
            {
                sb_bg_energy[i] = sb_bg_energy[i]*0.96f + frame_sb_energy[i]*0.04f;
            }
            st->bg_update_count++;
        }
        else
        {
            float a = 0.94f;
            if((t_bg_energy<frame_energy)&&24*st->frame_energy_smooth<frame_energy)
            {
                for(i=0; i<SNR_sb_num; i++)
                {
                    sb_bg_energy[i] = (float)(sb_bg_energy[i]*0.999f+frame_sb_energy[i]*0.001);
                }
            }
            else if(12*t_bg_energy<frame_energy)
            {
                st->t_bg_energy_sum  += frame_energy;
                st->tbg_energy_count++;
                for(i=0; i<SNR_sb_num; i++)
                {
                    sb_bg_energy[i] = sb_bg_energy[i]*0.96f + frame_sb_energy[i]*0.04f;
                }
            }
            else
            {
                if(t_bg_energy>frame_energy)
                {
                    a = 0.95f;
                    st->t_bg_energy_sum  += frame_energy;
                    st->tbg_energy_count++;
                    for(i=0; i<SNR_sb_num; i++)
                    {
                        sb_bg_energy[i] = sb_bg_energy[i]*a + frame_sb_energy[i]*(1-a);
                    }
                }
                else
                {
                    a = 0.96f;
                    st->t_bg_energy_sum  += frame_energy;
                    st->tbg_energy_count++;
                    for(i=0; i<SNR_sb_num; i++)
                    {
                        sb_bg_energy[i] = sb_bg_energy[i]*a + frame_sb_energy[i]*(1-a);
                    }

                }
            }
        }
    }
    else
    {
        if((t_bg_energy>500*frame_energy)&&(sb_bg_energy[0] >10*frame_sb_energy[0]))
        {
            for(i=0; i<SNR_sb_num; i++)
            {
                sb_bg_energy[i] = sb_bg_energy[i]*0.96f + frame_sb_energy[i]*0.04f;
            }

        }
        else if(t_bg_energy>10*frame_energy)
        {
            for(i=0; i<SNR_sb_num; i++)
            {
                sb_bg_energy[i] = sb_bg_energy[i]*0.999f + frame_sb_energy[i]*0.001f;
            }
        }
    }

    if(st->t_bg_energy_sum>160*st->tbg_energy_count)
    {
        st->t_bg_energy_sum = 160.0f*st->tbg_energy_count;
    }

    if(music_backgound_f==1&&st->lt_snr_org<3.2&&t_bg_energy>1&&update_flag==0)
    {
        for(i=0; i<SNR_sb_num; i++)
        {
            sb_bg_energy[i] = sb_bg_energy[i]*0.98f+0.000001f;
        }

    }
    if(music_backgound_f==1&&frame_energy<5000*t_bg_energy)
    {
        for(i=0; i<SNR_sb_num; i++)
        {
            sb_bg_energy[i] = sb_bg_energy[i]*0.98f+0.000001f;
        }
    }

    if(st->tbg_energy_count==64)
    {
        st->tbg_energy_count = 48;
        st->t_bg_energy_sum = st->t_bg_energy_sum*0.75f;
    }

    st->t_bg_energy	 = st->t_bg_energy_sum/st->tbg_energy_count;



}

void est_energy(float sb_power[],                  /*(o) energy of sub-band divided uniformly*/
                float frame_sb_energy[],           /*(o) energy of sub-band divided non-uniformly*/
                float *p_frame_energy,             /*(o) frame energy 1*/
                float *p_frame_energy2,            /*(o) frame energy 2*/
                float *p_high_energy,              /*(o) high frequency energy*/
                int bw                             /*(i) band width*/
               )
{
    int i, j;
    float frame_energy2 = 0.0f;
    float high_energy = 0.0f;
    int band_num = BAND_NUM_TAB[bw];
    const float sb_power_scale[5] = {0.0f, 0.16f, 0.24f, 0.28f, 0.28f};

    const int *Nregion_index;
    int SNR_sb_num;


    for(i=0; i<band_num; i++)
    {
        if(i>0 && (i!=band_num-1))
        {
            frame_energy2+=sb_power[i];
        }
        if(i>5)
        {
            high_energy += sb_power[i];
        }
    }

    high_energy /= (32768.f * 32768.f);
    frame_energy2 /= (32768.f * 32768.f);

    Nregion_index = REGION_INDEX[bw-CLDFBVAD_NB_ID];
    SNR_sb_num = ENERGY_BAND_NUM[bw-CLDFBVAD_NB_ID];

    for(i=0; i<SNR_sb_num; i++)
    {
        frame_sb_energy[i] = 0;
        for(j=Nregion_index[i]; j<Nregion_index[i+1]; j++)
        {
            frame_sb_energy[i] += sb_power[j];
        }

        frame_sb_energy[i] /= (32768.f * 32768.f);
    }
    *p_high_energy = high_energy;
    *p_frame_energy2 = frame_energy2;
    *p_frame_energy = frame_energy2 + (sb_power_scale[bw]*sb_power[0] / (32768.f * 32768.f));
}
