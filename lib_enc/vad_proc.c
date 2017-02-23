/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>

#include "prot.h"
#include "rom_enc.h"


/*-------------------------------------------------------------------*
 * vad_init()
 *
 *
 *-------------------------------------------------------------------*/

short vad_init(
    T_CldfbVadState *st
)
{
    float sSFM[SFM_NUM]= {0.88f,0.92f,0.92f};
    short i = 0;

    if(st == NULL)
    {
        return -1;
    }

    st->frameloop = 0;
    st->lt_snr_org = 1.0f;
    st->lf_snr_smooth = 5.0f;
    st->l_silence_snr = 0.5f;
    st->l_speech_snr = 5.0f;
    st->l_silence_snr_count = 1;
    st->l_speech_snr_count = 1;
    st->fg_energy = 16*(3.0518e-5f);
    st->bg_energy = 16*(4.6566e-10f);
    st->fg_energy_count = 16;
    st->bg_energy_count = 16;
    st->tonality_rate3 = 0.46f;
    st->music_background_rate = 0.46f;
    st->lt_noise_sp_center_diff_sum = 0.4f;
    st->lt_noise_sp_center_diff_counter = 4;
    st->lt_noise_sp_center0 = 1.8f;
    st->lt_noise_sp_center3 = 2.0f;
    st->lt_bg_highf_eng = 2.0f;
    st->t_bg_energy = 0.01f;
    st->t_bg_energy_sum = 0.01f;
    st->tbg_energy_count = 1;
    st->bg_update_count = 0;
    st->frame_energy_smooth = 1.0f;
    st->fg_energy_est_start = 0;
    st->speech_flag = 0;
    st->continuous_noise_num = 0;
    st->continuous_speech_num = 0;
    st->continuous_speech_num2 = 0;
    st->update_num_with_snr = 0;                                 /* the number of the background update with SNR*/
    st->update_count = 0;
    st->warm_hang_num = 0;
    for(i = 0; i < PRE_SNR_NUM; i++)
    {
        st->pre_snr[i] = 0.0f;
    }

    for(i = 0; i < POWER_NUM; i++)
    {
        st->frames_power[i] = 0;
    }

    for(i = 0; i < SPEC_AMP_NUM; i++)
    {
        st->smooth_spec_amp[i] = 0;
    }

    for(i = 0; i < SFM_NUM; i++)
    {
        st->sfm[i] = sSFM[i];
    }

    for(i = 0; i < SP_CENTER_NUM; i++)
    {
        st->sp_center[i] = 1.2f;
    }

    for(i = 0; i < STABLE_NUM; i++)
    {
        st->ltd_stable_rate[i] = 0.07f;
    }

    for(i = 0; i < BG_ENG_NUM; i++)
    {
        st->sb_bg_energy[i] = 0.01f;
        st->frame_sb_energy[i] = 0.001f;
    }

    for(i = 0; i < TONA_NUM; i++)
    {
        st->f_tonality_rate[i] = 0.48f;
    }

    for(i = 0; i < PRE_SPEC_DIF_NUM; i++)
    {
        st->pre_spec_low_dif[i] = 1.0f;
    }


    return 0;
}


/*-------------------------------------------------------------------*
 * UpdateState()
 *
 *
 *-------------------------------------------------------------------*/

static void UpdateState(
    T_CldfbVadState *st,
    float frame_energy,         /*(i) current frame energy                  */
    float high_eng,             /*(i) current frame high frequency energy   */
    int   update_flag,          /*(i) current frame update flag             */
    int   music_backgound_f,    /*(i) background music flag                 */
    int   vad_flag
    ,float snr
)
{
    st->frame_energy_smooth = st->frame_energy_smooth*0.95f +  frame_energy*0.05f;

    if(vad_flag==0)
    {
        st->lt_bg_highf_eng = st->lt_bg_highf_eng*0.95f + high_eng*0.05f;
    }

    if(st->frameloop<1000)
    {
        st->frameloop++;
    }
    background_update( st, frame_energy, update_flag, music_backgound_f, snr);
    if( vad_flag == 0)
    {
        st->continuous_speech_num2 = 0;

        if(st->continuous_noise_num > 10)
        {
            st->continuous_speech_num = 0;
        }
        else if (st->continuous_speech_num > 9)
        {
            st->continuous_speech_num = 9;
        }

        st->continuous_noise_num++;

        if(st->continuous_noise_num > 2048)
        {
            st->continuous_noise_num = 2048;
        }
    }
    else
    {
        st->continuous_noise_num = 0;
        st->continuous_speech_num2++;
        st->continuous_speech_num++;

        if(st->continuous_speech_num > 2048)
        {
            st->continuous_speech_num = 2048;
        }

        if(st->continuous_speech_num2 > 2048)
        {
            st->continuous_speech_num2 = 2048;
        }
    }

    return;
}


/*-------------------------------------------------------------------*
 * vad_proc()
 *
 *
 *-------------------------------------------------------------------*/

short vad_proc(
    float realValues[16][60],       /* CLDFB real values        */
    float imagValues[16][60],       /* CLDFB imag values        */
    float *sb_power,                /* Energy of CLDFB data     */
    int numBands,                   /* number of input bands    */
    T_CldfbVadState *vad_st,
    short *cldfb_addition,
    short vada_flag
)
{
    float frame_energy,frame_energy2;
    float spec_amp[8*10]; /* 120 */

    float snr, tsnr;
    int   update_flag;
    int   vad_flag;
    int   music_backgound_f=0;
    float HB_Power=0;
    float snr_flux;
    float lt_snr;
    float lt_snr_org;
    float lf_snr;
    int bandwidth;

    if(numBands<20)
    {
        bandwidth = 1;
    }
    else if(numBands<40)
    {
        bandwidth = 2;
    }
    else
    {
        bandwidth = 3;
    }

    vad_st->bw_index = bandwidth;

    assert(numBands>=10);

    /* new optimized structure */
    est_energy(sb_power, vad_st->frame_sb_energy, &frame_energy, &frame_energy2, &HB_Power, bandwidth);

    subband_FFT(realValues,imagValues,spec_amp);

    spec_center(sb_power,vad_st->sp_center,bandwidth);

    ltd_stable(vad_st->frames_power,vad_st->ltd_stable_rate,frame_energy,vad_st->frameloop);

    spec_flatness( spec_amp, vad_st->smooth_spec_amp, vad_st->sfm );

    frame_spec_dif_cor_rate( spec_amp, vad_st->pre_spec_low_dif, vad_st->f_tonality_rate );

    bg_music_decision( vad_st, &music_backgound_f, frame_energy );

    SNR_calc( vad_st->frame_sb_energy, vad_st->sb_bg_energy, vad_st->t_bg_energy, &snr, &tsnr, frame_energy2, bandwidth );

    calc_snr_flux( tsnr, vad_st->pre_snr, &snr_flux );

    calc_lt_snr( &lt_snr_org, &lt_snr, vad_st->fg_energy, vad_st->fg_energy_count,
                 vad_st->bg_energy, vad_st->bg_energy_count, bandwidth, vad_st->lt_noise_sp_center0 );

    calc_lf_snr( &vad_st->lf_snr_smooth, &lf_snr, vad_st->l_speech_snr, vad_st->l_speech_snr_count,
                 vad_st->l_silence_snr, vad_st->l_silence_snr_count, vad_st->fg_energy_count,vad_st->bg_energy_count, bandwidth );

    vad_flag = comvad_decision( vad_st, snr, tsnr, snr_flux, lt_snr, lt_snr_org, lf_snr, frame_energy2, music_backgound_f, cldfb_addition, vada_flag );

    update_flag = update_decision( vad_st,snr,tsnr, frame_energy,HB_Power,vad_st->vad_flag_for_bk_update, music_backgound_f );
    UpdateState( vad_st, frame_energy2, HB_Power, update_flag,music_backgound_f, vad_st->vad_flag_for_bk_update,snr);


    return vad_flag;
}

