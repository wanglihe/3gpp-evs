/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <stdlib.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"
#include "basop_util.h"
#include "basop_proto_func.h"


/*-------------------------------------------------------------------*
 * para_pred_bws()
 *
 * predict SWB parameters for bandwidth switching
 *-------------------------------------------------------------------*/
static short para_pred_bws(
    Decoder_State *st,              /* i/o: decoder state structure       */
    float *signal_wb,       /* i  : wideband frequency signal     */
    float *SWB_fenv,        /* o  : frequency-domain BWE envelope */
    short coder_type        /* i  : coding type                   */
)
{
    short i, j, k;
    short mode;
    float *input_hi;
    float peak, mean[7], mag, min;
    float avrg1, avrg2;
    float att;

    mode = NORMAL;

    k = 0;
    input_hi = &signal_wb[SHARP_WIDTH];
    for(i = 0; i < 7; i ++)
    {
        peak = 0.0f;
        mean[i] = 0;
        for(j = 0; j < SHARP_WIDTH; j ++)
        {
            mag = (float) fabs(*input_hi);
            if (mag > peak)
            {
                peak = mag;
            }
            mean[i] += mag;
            input_hi ++;
        }

        if (peak*(SHARP_WIDTH+3.5f) > 4.5f*mean[i] && peak > 8.0f)
        {
            k += 1;
        }
    }

    avrg1 = 0.0f;
    avrg2 = 0.0f;
    for(i=1; i<4; i++)
    {
        avrg1 += mean[i];
        avrg2 += mean[i+3];
    }
    avrg1 /= 3;
    avrg2 /= 3;

    min = FLT_MAX;
    peak = 0.0f;
    for(i = 4; i < 7; i ++)
    {
        if(mean[i] > 2.0f*avrg2)
        {
            mean[i] *= 2*avrg2/mean[i];
        }
        if(mean[i] < min)
        {
            min = mean[i];
        }
        if(mean[i] > peak)
        {
            peak = mean[i];
        }
    }

    if(st->tilt_wb > 8)
    {
        min = min(st->tilt_wb/15.0f, 1.0f)*peak;
    }

    if( peak == 0 || min == 0 )
    {
        set_f( SWB_fenv, 0, SWB_FENV );
    }
    else
    {
        for(i = 0; i < SWB_FENV; i ++)
        {
            SWB_fenv[i] = min*mean[i/5+4]/(64*peak);
        }
    }

    for(j = 0, i = SWB_FENV/2; i < SWB_FENV; i ++)
    {
        SWB_fenv[i] *= (1.0f - (float)j++/SWB_FENV);
    }

    if(avrg1 > 8.0f*avrg2)
    {
        for(i = 0; i < SWB_FENV; i ++)
        {
            SWB_fenv[i] *= 0.5f;
        }
    }
    if( st->last_core != HQ_CORE && st->last_codec_mode == MODE1 &&
            (st->enerLH > 0.5f*st->prev_enerLH && st->enerLH < 2.0f*st->prev_enerLH) &&
            (st->enerLL > 0.5f*st->prev_enerLL && st->enerLL < 2.0f*st->prev_enerLL) )
    {
        for(i=0; i<SWB_FENV; i++)
        {
            if ( st->prev_coder_type != coder_type && SWB_fenv[i] > 2.0f*st->prev_SWB_fenv[i] )
            {
                SWB_fenv[i] = 0.1f * SWB_fenv[i] + 0.9f * st->prev_SWB_fenv[i];
            }
            else
            {
                SWB_fenv[i] = st->attenu1 * SWB_fenv[i] + (1.0f - st->attenu1) * st->prev_SWB_fenv[i];
            }
        }

        if( st->attenu1 < 0.9f )
        {
            st->attenu1 += 0.05f;
        }
    }
    else
    {
        if( st->core_brate != st->last_core_brate || (st->enerLH > 0.5f*st->prev_enerLH && st->enerLH < 2.0f*st->prev_enerLH) ||
                (st->enerLL > 0.5f*st->prev_enerLL && st->enerLL < 2.0f*st->prev_enerLL) )
        {
            for(i=0; i<SWB_FENV; i++)
            {
                if(SWB_fenv[i] > 2.0f*st->prev_SWB_fenv[i])
                {
                    SWB_fenv[i] = st->prev_SWB_fenv[i];
                }
            }
        }

        for(i=0; i<SWB_FENV; i++)
        {
            SWB_fenv[i] = 0.9f*SWB_fenv[i] + 0.1f*st->prev_SWB_fenv[i];
        }

        st->attenu1 = 0.1f;
    }

    if( k > 3 )
    {
        mode = HARMONIC;
    }

    att = ((float)N_WS2N_FRAMES - (float)st->bws_cnt) / (float)N_WS2N_FRAMES;
    if(st->L_frame == L_FRAME16k )
    {
        for( i = 0; i < 4; i++ )
        {
            SWB_fenv[i] *= att;
        }
    }

    for( i=4; i<SWB_FENV; i++ )
    {
        SWB_fenv[i] *= att;
    }

    return mode;
}

/*-------------------------------------------------------------------*
 * WB_BWE_gain_deq()
 *
 * Decoding of WB parameters
 *-------------------------------------------------------------------*/
static
short WB_BWE_gain_deq(
    Decoder_State *st,              /* i/o: decoder state structure   */
    float *WB_fenv
)
{
    short mode;
    short index;

    index = (short) get_next_indice( st, 5 );
    mode = (short) get_next_indice( st, 1 ) + 2;

    WB_fenv[0] = (float)pow(2, 0.5f*F_2_5[2*index]);
    WB_fenv[1] = (float)pow(2, 0.5f*F_2_5[2*index+1]);

    return (mode);
}

/*-------------------------------------------------------------------*
 * wb_bwe_dec()
 *
 * WB BWE decoder (only for 16kHz signals)
 *-------------------------------------------------------------------*/

void wb_bwe_dec(
    float *synth,            /* i/o: ACELP core synthesis/final synthesis    */
    float *hb_synth,         /* o  : SHB synthesis/final synthesis           */
    const short output_frame,      /* i  : frame length                            */
    Decoder_State *st,               /* i/o: decoder state structure                 */
    short coder_type,        /* i  : coding type                             */
    float *voice_factors,    /* i  : voicing factors                         */
    const float pitch_buf[]        /* i  : pitch buffer                            */
)
{
    float ysynth[L_FRAME48k];                        /* MDCT spectrum of core synthesis          */
    float yerror[L_FRAME48k];                        /* MDCT spectrum of error                   */
    float wtda_synth[2*L_FRAME48k];
    float WB_fenv[SWB_FENV];
    short mode;
    short i;

    /* MDCT of the core synthesis signal */
    wtda( synth, wtda_synth, st->old_wtda_swb, ALDO_WINDOW,ALDO_WINDOW, output_frame );
    direct_transform( wtda_synth, ysynth, 0, output_frame );

    if( !st->bfi )
    {
        if( st->total_brate == ACELP_13k20 )
        {
            /* de-quantization */
            mode = WB_BWE_gain_deq( st, WB_fenv );
            st->last_wb_bwe_ener = 0.5f*(WB_fenv[0] + WB_fenv[1]);
        }
        else
        {
            if( st->last_extl != WB_BWE )
            {
                st->prev_SWB_fenv[0] = 0.0f;
            }

            mode = WB_BWE_gain_pred( WB_fenv, ysynth, coder_type, st->prev_coder_type, st->prev_SWB_fenv[0], voice_factors, pitch_buf,
                                     st->last_core_brate, st->last_wb_bwe_ener, st->last_extl
                                     ,st->tilt_wb
                                   );
        }
    }
    else
    {
        /* FEC */
        mode = NORMAL;
        for( i=0; i<2; i++ )
        {
            WB_fenv[i] = 0.75f*st->prev_SWB_fenv[i];
        }
    }

    if( st->last_extl != WB_BWE || st->bfi )
    {
        mvr2r( WB_fenv, st->prev_SWB_fenv, 2 );
    }

    /* reconstruction of MDCT spectrum of the error signal */
    WB_BWE_decoding( ysynth, WB_fenv, yerror, L_FRAME16k, mode, st->last_extl,
                     &st->prev_Energy_wb, st->prev_SWB_fenv, &st->prev_L_swb_norm, st->extl,
                     coder_type, st->total_brate, &st->Seed, &st->prev_flag, st->prev_coder_type );

    if ( st->output_Fs == 32000)
    {
        set_f( &yerror[L_FRAME16k], 0, L_FRAME16k );
    }
    else if ( st->output_Fs == 48000 )
    {
        set_f( &yerror[L_FRAME16k], 0, L_FRAME32k );
    }

    inverse_transform( yerror, wtda_synth, 0, output_frame, -1 );

    window_ola( wtda_synth, hb_synth, st->mem_imdct,   output_frame,ALDO_WINDOW,ALDO_WINDOW, 0,0,0 );

    st->prev_mode = mode;

    return;
}

/*-------------------------------------------------------------------*
 * swb_bwe_gain_deq()
 *
 * Decoding of SWB parameters
 *-------------------------------------------------------------------*/

short swb_bwe_gain_deq(     /* o  : BWE class                       */
    Decoder_State *st,        /* i/o: decoder state structure   */
    const short core,       /* i  : core                            */
    float *SWB_tenv,  /* o  : time-domain BWE envelope        */
    float *SWB_fenv,  /* o  : frequency-domain BWE envelope   */
    const short hr_flag,    /* i  : high rate flag                  */
    const short hqswb_clas  /* i  : HQ BWE class                    */
)
{
    short index, mode, n_band;
    short indice[6];
    float quant_tmp[SWB_FENV/2], quant_tmp2[SWB_FENV/2];
    short nb_bits[6];
    short nenv;

    if ( hqswb_clas > 0)
    {
        mode = (short)get_next_indice( st, 1 );
        if (mode == 0)
        {
            mode = (short)get_next_indice( st, 1 );
        }
        else
        {
            mode = HQ_GENERIC_SP_EXC;
        }
    }
    else
    {
        mode = (short)get_next_indice( st, 2 );
    }

    if( mode == 1 && core == ACELP_CORE )
    {
        for( n_band = 0; n_band < SWB_TENV; n_band++ )
        {
            index = (short)get_next_indice( st, 4 );
            SWB_tenv[n_band] = (float)(1 << index);
        }

        indice[0] = (short)get_next_indice( st, 7 );
        indice[1] = (short)get_next_indice( st, 6 );

        for(n_band = 0; n_band < DIM_TR1; n_band++)
        {
            quant_tmp[2*n_band]= Env_TR_Cdbk1[indice[0]*DIM_TR1+n_band];
        }

        quant_tmp[1] = (quant_tmp[0]+quant_tmp[2])*0.5f+Env_TR_Cdbk2[indice[1]*DIM_TR2];
        quant_tmp[3] = quant_tmp[2]+Env_TR_Cdbk2[indice[1]*DIM_TR2+1];

        for(n_band = 0; n_band < SWB_FENV_TRANS; n_band++)
        {
            SWB_fenv[n_band] = (float)pow(10, 0.025f*(quant_tmp[n_band]+Mean_env_tr[n_band]));
        }

        /* in case of band-width switching, attenuate frame gain */
        if( st->bws_cnt1 > 0 )
        {
            for(n_band = 0; n_band < SWB_TENV; n_band++)
            {
                SWB_tenv[n_band] *= (float)st->bws_cnt1 / (float)N_NS2W_FRAMES;
            }

            for (n_band = 0; n_band < SWB_FENV_TRANS; n_band++)
            {
                SWB_fenv[n_band] *= (float)st->bws_cnt1 / (float)N_NS2W_FRAMES;
            }
        }
    }
    else
    {
        nb_bits[0] = 5;
        nb_bits[1] = 7;
        nb_bits[2] = 6;
        nb_bits[3] = 5;

        if ( hr_flag == 1 )
        {
            nb_bits[4] = 5;
            nenv = SWB_FENV - 2;
        }
        else
        {
            nb_bits[4] = 6;
            nenv = SWB_FENV;
        }

        for (n_band = 0; n_band < 5; n_band++)
        {
            indice[n_band] = (short) get_next_indice( st, nb_bits[n_band] );
        }

        if ( hqswb_clas == HQ_GEN_FB )
        {
            indice[n_band] = (short) get_next_indice( st, 5 );
        }

        mvr2r( &EnvCdbk11[indice[0] * DIM11], quant_tmp, DIM11 );
        mvr2r( &EnvCdbk1st[indice[1] * DIM1ST], quant_tmp2, DIM1ST );
        mvr2r( &EnvCdbk2nd[indice[2] * DIM2ND], quant_tmp2+DIM1ST, DIM2ND );

        for( n_band = 0; n_band < DIM11-1; n_band++ )
        {
            quant_tmp[n_band] += quant_tmp2[n_band];
            SWB_fenv[n_band*2] = quant_tmp[n_band];
        }

        if ( hr_flag == 1 )
        {
            quant_tmp[6] += quant_tmp2[6];
            SWB_fenv[11] = quant_tmp[6];


            mvr2r( &EnvCdbk3rd[indice[3] * DIM3RD], quant_tmp2, DIM3RD );
            mvr2r( &EnvCdbk3rd[indice[4] * DIM3RD], quant_tmp2+DIM3RD, DIM3RD );

            for(n_band = 0; n_band < 5; n_band++)
            {
                SWB_fenv[n_band*2+1] = ((quant_tmp[n_band]+quant_tmp[n_band+1])/2.f)+quant_tmp2[n_band+1];
            }

            SWB_fenv[0] += quant_tmp2[0];
        }
        else
        {
            quant_tmp[DIM11-1]+=quant_tmp2[DIM11-1];
            SWB_fenv[(DIM11-1)*2] = quant_tmp[DIM11-1];

            mvr2r( &EnvCdbk3rd[indice[3] * DIM3RD], quant_tmp2, DIM3RD );
            mvr2r( &EnvCdbk4th[indice[4] * DIM4TH], quant_tmp2+DIM3RD, DIM4TH );

            for( n_band = 0; n_band < DIM12-1; n_band++ )
            {
                SWB_fenv[n_band*2+1] = ((quant_tmp[n_band]+quant_tmp[n_band+1])/2.f)+quant_tmp2[n_band];
            }

            SWB_fenv[n_band*2+1] = quant_tmp[n_band]+quant_tmp2[n_band];
        }

        for( n_band = 0; n_band < nenv; n_band++ )
        {
            Word16 tmp,frac,exp;
            Word32 L_tmp;
            tmp = add((short)(SWB_fenv[n_band]*256),(short)(Mean_env[n_band]*256)); /*Q8 */

            L_tmp = L_mult(tmp, 21771); /* 0.166096 in Q17 -> Q26 */
            L_tmp = L_shr(L_tmp, 10); /* From Q26 to Q16 */
            frac = L_Extract_lc(L_tmp, &exp); /* Extract exponent of L_tmp */

            tmp = extract_l(Pow2(13, frac));/* Put 13 as exponent so that */
            /* output of Pow2() will be: */
            /* 16384 < Pow2() <= 32767 */
            exp = sub(exp, 13);
            tmp = shl(tmp, add(exp,1)); /*Q1 */
            SWB_fenv[n_band] = (float)tmp*0.5f;; /*Q1 */
        }

        if ( hqswb_clas == HQ_GEN_FB )
        {
            mvr2r( &EnvCdbkFB[indice[5] * DIM_FB], &SWB_fenv[nenv], DIM_FB );
            for( n_band = 0; n_band < DIM_FB; n_band++ )
            {
                Word16 tmp,frac,exp;
                Word32 L_tmp;

                tmp = add((short)(SWB_fenv[n_band + nenv]*128),(short)(Mean_env_fb[n_band]*128)); /*Q7 */
                L_tmp = L_mult(tmp, 21771); /* 0.166096 in Q17 -> Q25 */
                L_tmp = L_shr(L_tmp, 9); /* From Q25 to Q16 */
                frac = L_Extract_lc(L_tmp, &exp); /* Extract exponent of L_tmp */

                tmp = extract_l(Pow2(13, frac));/* Put 13 as exponent so that */
                /* output of Pow2() will be: */
                /* 16384 < Pow2() <= 32767 */
                exp = sub(exp, 13);
                tmp = shl(tmp, add(exp,1));
                move16();
                SWB_fenv[add(n_band,nenv)] = (float)tmp*0.5f;
            }
        }

    }

    return mode;
}

/*-------------------------------------------------------------------*
 * swb_bwe_dec()
 *
 * SWB BWE decoder (only for 32kHz signals)
 *-------------------------------------------------------------------*/

void swb_bwe_dec(
    Decoder_State *st,                /* i/o: decoder state structure                 */
    const float *synth,             /* i  : ACELP core synthesis/final synthesis    */
    float *hb_synth,          /* o  : SHB synthesis/final synthesis           */
    const short output_frame        /* i  : frame length                            */
    ,short coder_type               /* i  : coding type                             */
)
{
    short i, l_subfr;
    float ysynth[L_FRAME48k];
    float yerror[L_FRAME48k];
    float wtda_synth[2*L_FRAME48k];
    float SWB_tenv[SWB_TENV];
    float SWB_fenv[SWB_FENV];
    short L;
    short mode;
    short frica_flag = 0;
    float fb_ener_adjust = 0.0f;
    short j = 0;
    float ener_adjust_quan;
    short idxGain;

    /*---------------------------------------------------------------------*
     * SWB BWE decoding
     *---------------------------------------------------------------------*/

    /* windowing of the ACELP core synthesis */
    wtda( synth, wtda_synth, st->old_wtda_swb, ALDO_WINDOW, ALDO_WINDOW, output_frame );

    /* DCT of the ACELP core synthesis */
    direct_transform( wtda_synth, ysynth, 0, output_frame );

    if( !st->bfi )
    {
        if( st->bws_cnt > 0 )
        {
            /* estimate parameters */
            mode = para_pred_bws( st, ysynth, SWB_fenv, coder_type );
        }
        else
        {
            /* de-quantization */
            mode = swb_bwe_gain_deq( st, ACELP_CORE, SWB_tenv, SWB_fenv, 0, -1 );
        }

        L = mode == TRANSIENT ? SWB_FENV_TRANS : SWB_FENV;
        st->prev_ener_shb = 0.0f;
        for(i=0; i<L; i++)
        {
            st->prev_ener_shb += SWB_fenv[i];
        }
        st->prev_ener_shb /= L;
    }
    else
    {
        /* SHB FEC */
        if( st->prev_mode != TRANSIENT )
        {
            mode = st->prev_mode;
        }
        else
        {
            mode = NORMAL;
        }

        mvr2r( st->prev_SWB_fenv, SWB_fenv, SWB_FENV );
    }

    /* reconstruction of MDCT spectrum of the error signal */
    set_f( yerror, 0, output_frame );

    if (st->L_frame == L_FRAME16k )
    {
        SWB_BWE_decoding( ysynth, SWB_fenv, yerror, L_FRAME32k-80, mode, &frica_flag, &st->prev_Energy, st->prev_SWB_fenv,
                          &st->prev_L_swb_norm, st->tilt_wb, &st->Seed, 80, &st->prev_weight, st->extl
                          ,st->last_extl
                        );
    }
    else
    {
        SWB_BWE_decoding( ysynth, SWB_fenv, yerror, L_FRAME32k-80, mode, &frica_flag, &st->prev_Energy, st->prev_SWB_fenv,
                          &st->prev_L_swb_norm, st->tilt_wb, &st->Seed, 6, &st->prev_weight, st->extl
                          ,st->last_extl
                        );
    }

    if ( st->prev_frica_flag == 1 && frica_flag == 0 )
    {
        for( i = 0; i < L_SUBFR; i++ )
        {
            st->mem_imdct[i] *= 1.0f - i*0.015625f;
        }

        for( ; i < output_frame; i++ )
        {
            st->mem_imdct[i] = 0.0f;
        }
    }

    /* decode information */
    if ( st->extl == FB_BWE )
    {
        if( !st->bfi )
        {
            idxGain = (short)get_next_indice( st,NUM_BITS_FB_FRAMEGAIN);
            fb_ener_adjust = usdequant(idxGain, FB_GAIN_QLOW, FB_GAIN_QDELTA);
        }
        else if( st->bfi )
        {
            fb_ener_adjust = st->prev_fb_ener_adjust;
        }

        st->prev_fb_ener_adjust = fb_ener_adjust;
        if(mode == TRANSIENT)
        {
            ener_adjust_quan = fb_ener_adjust;
        }
        else
        {
            if(SWB_fenv[7] < 0.01f)
            {
                ener_adjust_quan = 0.0f;
            }
            else
            {
                ener_adjust_quan = min(SWB_fenv[13]/SWB_fenv[7], 4.0f);
            }
        }
        for( i = FB_BAND_BEGIN; i < FB_BAND_BEGIN+DE_OFFSET1; i++ )
        {
            yerror[i] = yerror[i-FB_BAND_WIDTH] * ((1.0f-j*FB_GAIN_QDELTA)*ener_adjust_quan + j*FB_GAIN_QDELTA*fb_ener_adjust);
            j++;
        }
        for(; i < FB_BAND_END; i++ )
        {
            yerror[i] = yerror[i-FB_BAND_WIDTH] * fb_ener_adjust;
        }
    }

    /* iDCT of the error signal */
    inverse_transform( yerror, wtda_synth, 0, output_frame, -1 );

    /* inverse windowing of the error signal */
    window_ola( wtda_synth, hb_synth, st->mem_imdct,   output_frame,ALDO_WINDOW,ALDO_WINDOW, 0,0,0 );
    l_subfr = output_frame/4;

    if( mode == TRANSIENT )
    {
        for(i = 0; i < SWB_TENV; i++)
        {
            SWB_tenv[i] *= 0.8f;
        }

        /* time envelope shaping when the current frame is TRANSIENT frame */
        time_envelop_shaping( hb_synth, SWB_tenv, output_frame );

        st->prev_td_energy = SWB_tenv[3];
    }
    else if( frica_flag == 1 && st->prev_frica_flag == 0 )
    {
        time_reduce_pre_echo( synth, hb_synth, st->prev_td_energy, l_subfr );
    }
    else
    {
        st->prev_td_energy = 0.0f;
        for( i=0; i<l_subfr; i++ )
        {
            st->prev_td_energy += hb_synth[3*l_subfr+i]*hb_synth[3*l_subfr+i];
        }
        st->prev_td_energy = (float)sqrt(st->prev_td_energy/l_subfr);
    }

    st->prev_frica_flag = frica_flag;
    st->prev_mode = mode;

    return;
}
