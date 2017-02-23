/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <stdio.h>
#include "prot.h"
#include "rom_enc.h"






void SNR_calc(float frame_sb_energy[],  /*(i) energy of sub-band divided non-uniformly*/
              float sb_bg_energy[],     /*(i) sub-band background energy*/
              float t_bg_energy,        /*(i) time background energy of several frames*/
              float *snr,               /*(o) frequency domain SNR */
              float *tsnr,              /*(o) time domain SNR */
              float frame_energy,       /*(i) current frame energy */
              int bandwidth             /*(i) band width*/
             )
{
    int i;
    float snr_tmp,tmp;
    int SNR_sb_num;
    SNR_sb_num = ENERGY_BAND_NUM[bandwidth-CLDFBVAD_NB_ID];

    snr_tmp  = 0;
    for(i=0; i< SNR_sb_num; i++)
    {
        tmp  =  (frame_sb_energy[i]+0.0001f)/(sb_bg_energy[i]+0.0001f);
        tmp  =  (float)log10(tmp);
        if(tmp>-0.1)
        {
            tmp  =  tmp*3.3219f;
            snr_tmp += tmp;
        }
    }
    if (snr_tmp<0)
    {
        snr_tmp = 0;
    }
    *snr = snr_tmp/SNR_sb_num;
    tmp  = (frame_energy+0.0001f)/(t_bg_energy+0.0001f);
    tmp = (float)log10(tmp);
    *tsnr = tmp*3.3219f;
    if(bandwidth == CLDFBVAD_SWB_ID)
    {
        tmp  = (frame_energy)/(t_bg_energy+FLT_MIN);
        tmp = (float)log10(tmp+FLT_MIN);
        *tsnr = tmp*3.3219f;
    }

}

void calc_snr_flux(float tsnr,                  /*(i) time-domain SNR*/
                   float pre_snr[],             /*(io)time-domain SNR storage*/
                   float *snr_flux              /*(o) average tsnr*/
                  )
{
    int i;
    float snr_sum = 0.0f;


    if (tsnr < 2.6f&&tsnr>0)
        pre_snr[0] = tsnr;
    else if(tsnr<=0)
        pre_snr[0] = 0;
    else
        pre_snr[0] = 2.6f;

    snr_sum =0;
    for(i=0; i<PRE_SNR_NUM; i++)
    {
        snr_sum += pre_snr[i];
    }
    *snr_flux = snr_sum/PRE_SNR_NUM;
    for(i=PRE_SNR_NUM-1; i>0; i--)
    {
        pre_snr[i] = pre_snr[i-1];
    }

}

void calc_lt_snr(float *lt_snr_org,          /*(o) original long time SNR*/
                 float *lt_snr,              /*(o) long time SNR calculated by fg_energy and bg_energy*/
                 float fg_energy,            /*(i) foreground energy sum  */
                 int   fg_energy_count,      /*(i) number of the foreground energy frame */
                 float bg_energy,            /*(i) background energy sum  */
                 int   bg_energy_count,      /*(i) number of the background energy frame */
                 int   bw_index,             /*(i) band width index*/
                 float lt_noise_sp_center0   /*(i) long time noise spectral center by 0*/
                )
{
    float tmp_lt_noise_sp_center;
    float rtn_lt_snr;

    const float offset =  -0.00156247615814208984375f;


    tmp_lt_noise_sp_center = lt_noise_sp_center0-1.4f;
    if(tmp_lt_noise_sp_center>0.8)
    {
        tmp_lt_noise_sp_center = 0.8f;
    }
    if(tmp_lt_noise_sp_center<0)
    {
        tmp_lt_noise_sp_center = 0.0f;
    }
    rtn_lt_snr = (float)log10((fg_energy*bg_energy_count+FLT_MIN)/(bg_energy*fg_energy_count+FLT_MIN));
    *lt_snr_org = rtn_lt_snr;

    if(bg_energy_count<56||fg_energy_count<56)
    {
        rtn_lt_snr = 2.1f;
    }

    if(bw_index == CLDFBVAD_NB_ID)
    {
        rtn_lt_snr = (rtn_lt_snr-1.5f)*0.5f;
    }
    else if(bw_index == CLDFBVAD_WB_ID)
    {
        rtn_lt_snr = (rtn_lt_snr-1.5f)*0.50f;
    }
    else
    {
        rtn_lt_snr = (rtn_lt_snr-1.5f)*0.46f;

    }
    rtn_lt_snr = rtn_lt_snr + (rtn_lt_snr*0.4f+offset)*tmp_lt_noise_sp_center*0.4f;
    if(rtn_lt_snr<0)
    {
        rtn_lt_snr = 0.0f;
    }

    if(rtn_lt_snr>2.0)
    {
        rtn_lt_snr = 2.0f;
    }

    *lt_snr = rtn_lt_snr;


}

void calc_lf_snr(float *lf_snr_smooth,       /*(o) smoothed lf_snr*/
                 float *lf_snr,              /*(o) long time frequency domain
                                                       SNR calculated by l_speech_snr and l_silence_snr*/
                 float l_speech_snr,         /*(i) sum of active frames snr */
                 int   l_speech_snr_count,   /*(i) number of the active frame  */
                 float l_silence_snr,        /*(i) sum of the nonactive frames snr*/
                 int   l_silence_snr_count,  /*(i) number of the nonactive frame */
                 int   fg_energy_count,      /*(i) number of the foreground energy frame */
                 int   bg_energy_count,      /*(i) number of the background energy frame */
                 int   bw_index              /*(i) band width index*/
                )
{
    float l_snr;
    l_snr  = l_speech_snr/l_speech_snr_count - l_silence_snr/l_silence_snr_count;
    *lf_snr_smooth = *lf_snr_smooth*0.9f +  0.1f*l_snr;

    if(bg_energy_count<56||fg_energy_count<56)
    {
        l_snr  = 4.8f;
    }


    l_snr  = (l_snr - 3.0f)*0.12f;

    if(l_snr<0)
    {
        l_snr = 0;
    }

    if(l_snr>MAX_LF_SNR_TAB[bw_index])
    {
        l_snr = MAX_LF_SNR_TAB[bw_index];
    }

    *lf_snr = l_snr;
}

