/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <stdlib.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"
#include "rom_dec.h"

/*-----------------------------------------------------------------*
 * Local functions
 *-----------------------------------------------------------------*/

static void dequantizeSHBparams( Decoder_State *st, const short extl, long  extl_brate, float *Q_lsf, float *Q_subgain,
                                 float *Q_framegrain, short *uv_flag, float *Q_shb_ener_sf, float *Q_shb_res_gshape, float *Q_mixFactors );
static void Dequant_lower_LSF( const short lsf_idx[], float lsf_q[] );
static void Map_higher_LSF(float lsf_q[], const float m, const float grid_in[]);
static void Dequant_mirror_point(const float lsf_q[], const short m_idx, float *m);

void InitSWBdecBuffer(
    Decoder_State *st  /* i/o: SHB decoder structure */
)
{
    set_f(st->old_bwe_exc, 0.0f, (PIT16k_MAX * 2));
    st->bwe_seed[0] = 23; /* 1; */
    st->bwe_seed[1] = 59; /* 10000; */
    set_f( st->old_bwe_exc_extended, 0.0f, NL_BUFF_OFFSET );
    st->bwe_non_lin_prev_scale = 0;
    st->last_voice_factor = 0.0f;

    set_f(st->genSHBsynth_Hilbert_Mem, 0.0f,HILBERT_MEM_SIZE);
    set_f(st->genSHBsynth_state_lsyn_filt_shb_local, 0.0f,2 * ALLPASSSECTIONS_STEEP);

    st->syn_dm_phase = 0;
    st->prev_fbbwe_ratio = 1.0f;
    st->prev_wb_bwe_frame_pow = 0.001f;
    st->prev_swb_bwe_frame_pow = 0.001f;
    st->prev_fb_ener_adjust = 0.0f;

    return;
}

void ResetSHBbuffer_Dec(
    Decoder_State *st   /* i/o: SHB encoder structure */
)
{
    short i;
    float f;
    float inc;

    if ( st->extl != WB_TBE )
    {
        f = 1.0f/22.0f;
        inc = 1.0f/22.0f;
    }
    else
    {
        f = 1.0f/6.0f;
        inc = 1.0f/6.0f;
    }

    /* states for the filters used in generating SHB excitation from WB excitation*/
    set_f(st->mem_csfilt,0,2);

    /* states for the filters used in generating SHB signal from SHB excitation*/
    set_f(st->state_syn_shbexc,0,L_SHB_LAHEAD);
    set_f( st->state_lpc_syn, 0, LPC_SHB_ORDER);
    if( sub(st->extl, FB_TBE) == 0 )
    {
        set_f( st->fb_state_lpc_syn, 0, LPC_SHB_ORDER );
        st->fb_tbe_demph = 0;
        fb_tbe_reset_synth( st->fbbwe_hpf_mem,&st->prev_fbbwe_ratio);
    }
    /* states for the filters used in generating SHB signal from SHB excitation in wideband*/
    set_f( st->mem_genSHBexc_filt_down_shb, 0.0f, (2*ALLPASSSECTIONS_STEEP+1) );
    set_f( st->mem_genSHBexc_filt_down_wb2, 0.0f, (2*ALLPASSSECTIONS_STEEP+1) );
    set_f( st->mem_genSHBexc_filt_down_wb3, 0.0f, (2*ALLPASSSECTIONS_STEEP+1) );

    set_f( st->state_lsyn_filt_shb,0, 2 * ALLPASSSECTIONS_STEEP );
    set_f( st->state_lsyn_filt_dwn_shb,0, 2 * ALLPASSSECTIONS_STEEP );
    set_f( st->mem_resamp_HB, 0, INTERP_3_1_MEM_LEN );

    /* States for the local synthesis filters */
    set_f(st->syn_overlap,0,L_SHB_LAHEAD);

    /* States for FEC */
    if ( st->extl != WB_TBE )
    {
        for (i=0; i<LPC_SHB_ORDER; i++)
        {
            st->lsp_prevfrm[i] = f;
            f += inc;
        }
    }
    else
    {
        for (i=0; i<LPC_SHB_ORDER_WB; i++)
        {
            st->lsp_prevfrm[i] = f;
            f += inc;
        }
    }

    st->GainFrame_prevfrm = 0.0f;
    st->GainAttn = 1.0;
    st->tbe_demph = 0.0f;
    st->tbe_premph = 0.0f;
    set_f(st->mem_stp_swb, 0, LPC_SHB_ORDER);
    st->gain_prec_swb = 1.0f;
    set_f( &st->GainShape_Delay[0], 0, NUM_SHB_SUBFR/2 );

    set_f(st->old_core_synth, 0, L_FRAME16k);
    set_f(st->old_tbe_synth, 0, L_SHB_TRANSITION_LENGTH);
    st->tilt_swb_fec = 0.0f;

    return;
}

/*-------------------------------------------------------------------*
 * wb_tbe_dec()
 *
 * WB TBE decoder, 6 - 8 kHz band decoding module
 *-------------------------------------------------------------------*/

void wb_tbe_dec(
    Decoder_State *st,              /* i/o: decoder state structure */
    const short coder_type,       /* i  : coding type */
    const float *bwe_exc_extended,/* i  : bandwidth extended excitation  */
    const float voice_factors[],  /* i  : voicing factors                */
    float *synth            /* o  : WB synthesis/final synthesis   */
)
{
    short i;
    float shaped_wb_excitation [ (L_FRAME16k + L_SHB_LAHEAD)/4 ];
    float exc4kWhtnd [L_FRAME16k / 4];
    float lsf_wb[LPC_SHB_ORDER_WB], lpc_wb[LPC_SHB_ORDER_WB + 1], GainShape[NUM_SHB_SUBFR], GainFrame ;
    float error[L_FRAME16k];
    float upsampled_synth[L_FRAME48k];
    float prev_pow, curr_pow, scale;
    float temp, curr_frame_pow;
    short j;
    float vf_modified[NB_SUBFR16k];
    short uv_flag = 0;

    if( st->bws_cnt == 0 )
    {

        if( !st->bfi )
        {
            if( st->use_partial_copy )
            {
                if(st->last_extl != WB_TBE)
                {
                    st->GainFrame_prevfrm = 0;
                    st->lsp_prevfrm[0] = 0.1f;
                    for (i=1; i<LPC_SHB_ORDER_LBR_WB; i++)
                    {
                        st->lsp_prevfrm[i] = st->lsp_prevfrm[i-i] + 0.1f;
                    }
                }

                mvr2r( st->lsp_prevfrm, lsf_wb, LPC_SHB_ORDER_LBR_WB );
                set_f( GainShape, RECIP_ROOT_EIGHT, NUM_SHB_SUBFR/2 );

                if( st->rf_frame_type == RF_NELP )
                {
                    /* Frame gain */
                    /* only four LSBs are valid */
                    st->rf_indx_tbeGainFr &= 0xF;
                    mvr2r( SHBCB_FrameGain16 + st->rf_indx_tbeGainFr, &GainFrame, 1 );
                    if( st->core == ACELP_CORE && st->last_core == ACELP_CORE && !st->prev_use_partial_copy && st->prev_coder_type == UNVOICED && GainFrame != st->GainFrame_prevfrm && st->last_extl == WB_TBE)
                    {
                        GainFrame = 0.2f*GainFrame + 0.8f*st->GainFrame_prevfrm;
                    }
                }
                else
                {
                    /* Frame gain */
                    temp = 0.0f;
                    switch (st->rf_indx_tbeGainFr)
                    {
                    case 0:
                        GainFrame = 0.5f;
                        if(st->GainFrame_prevfrm <= 1.25) temp = 0.8f;
                        break;
                    case 1:
                        GainFrame = 2.0f;
                        if(st->GainFrame_prevfrm > 1.25 && st->GainFrame_prevfrm <= 3) temp = 0.8f;
                        break;
                    case 2:
                        GainFrame = 4.0f;
                        if(st->GainFrame_prevfrm > 3 && st->GainFrame_prevfrm <= 6) temp = 0.8f;
                        break;
                    case 3:
                        GainFrame = 8.0f;
                        if(st->GainFrame_prevfrm > 6 && st->GainFrame_prevfrm <= 16) temp = 0.8f;
                        break;
                    default:
                        GainFrame = 1.0f;
                        fprintf(stderr, "RF WB-TBE gain bits not supported.");
                        break;
                    }

                    if(st->last_extl == WB_TBE)
                    {
                        GainFrame = (1 - temp)*GainFrame + temp*(st->GainFrame_prevfrm);
                    }

                    if (st->core == ACELP_CORE && st->last_core == ACELP_CORE)
                    {
                        if (!st->prev_use_partial_copy && st->last_coder_type == VOICED && st->rf_frame_type == RF_GENPRED
                                && st->prev_tilt_code_dec < 0.046f && st->prev_tilt_code_dec > 0.006f )
                        {
                            GainFrame *= 0.3f;
                        }
                    }
                }
            }
            else
            {
                /* de-quantization */
                dequantizeSHBparams( st, st->extl, st->extl_brate, lsf_wb, GainShape, &GainFrame, &uv_flag, 0, 0, 0 );
            }
        }
        else
        {
            if ( st->extl_brate == WB_TBE_0k35 )
            {
                mvr2r( st->lsp_prevfrm, lsf_wb, LPC_SHB_ORDER_LBR_WB );
            }
            else
            {
                mvr2r( st->lsp_prevfrm, lsf_wb, LPC_SHB_ORDER_WB );
            }
            set_f( GainShape, RECIP_ROOT_EIGHT, NUM_SHB_SUBFR/2 );
            st->GainAttn *= 0.85f;
            if( st->codec_mode == MODE1 )
            {
                GainFrame = st->GainAttn * st->GainFrame_prevfrm;
            }
            else
            {
                GainFrame =  st->GainFrame_prevfrm;
            }
        }

        if ( st->extl_brate == WB_TBE_0k35 )
        {
            /* convert LSPs back into LP coeffs */
            lsp2a( lpc_wb, lsf_wb, LPC_SHB_ORDER_LBR_WB );
            set_f( lpc_wb + LPC_SHB_ORDER_LBR_WB+1, 0.0f, (LPC_SHB_ORDER_WB - LPC_SHB_ORDER_LBR_WB) );
        }
        else
        {
            /* convert LSPs back into LP coeffs */
            lsp2a( lpc_wb, lsf_wb, LPC_SHB_ORDER_WB );
        }
        lpc_wb[0] = 1.0f;
        mvr2r( voice_factors, vf_modified, NB_SUBFR16k );
        if( coder_type == VOICED )
        {
            for (i = 1; i < NB_SUBFR; i++)
            {
                vf_modified[i] = 0.8f * voice_factors[i] + 0.2f * voice_factors[i-1];
            }

            if(st->L_frame != L_FRAME )
            {
                vf_modified[4] = 0.8f * voice_factors[4] + 0.2f * voice_factors[3];
            }
        }

        if(st->use_partial_copy && st->nelp_mode_dec)
        {
            set_f( vf_modified, 0.0f, NB_SUBFR16k );
        }

        /* From low band excitation, generate highband excitation */
        mvr2r( st->state_syn_shbexc, shaped_wb_excitation, L_SHB_LAHEAD / 4);

        GenShapedWBExcitation( shaped_wb_excitation+ L_SHB_LAHEAD / 4, lpc_wb, exc4kWhtnd, st->mem_csfilt, st->mem_genSHBexc_filt_down_shb,
                               st->mem_genSHBexc_filt_down_wb2, st->mem_genSHBexc_filt_down_wb3, st->state_lpc_syn, coder_type,
                               bwe_exc_extended, st->bwe_seed, vf_modified, uv_flag, st->igf );

        prev_pow = sum2_f( shaped_wb_excitation, L_SHB_LAHEAD/4 );
        curr_pow = sum2_f( shaped_wb_excitation + L_SHB_LAHEAD/4, L_SHB_LAHEAD/4 );

        if( voice_factors[0] > 0.75f )
        {
            curr_pow *= 0.25;
        }

        if( prev_pow == 0 )
        {
            scale = 0;
        }
        else
        {
            scale = sqrt(curr_pow / prev_pow);
        }
        for (i=0; i<L_SHB_LAHEAD/4 - 1; i++)
        {
            shaped_wb_excitation[i] *= scale;
        }
        scale = sqrt(scale);

        shaped_wb_excitation[L_SHB_LAHEAD/4 - 1] *= scale;

        /* Update SHB excitation */
        mvr2r( shaped_wb_excitation + L_FRAME16k/4, st->state_syn_shbexc, L_SHB_LAHEAD/4 );

        /* Adjust the subframe and frame gain of the synthesized shb signal */
        /* Scale the shaped excitation */
        ScaleShapedSHB( SHB_OVERLAP_LEN/2, shaped_wb_excitation, st->syn_overlap, GainShape, GainFrame, window_wb, subwin_wb );

        curr_frame_pow = sum2_f( shaped_wb_excitation, L_FRAME16k/4 ) + 0.001f;

        if( !st->bfi && (st->prev_bfi || st->prev_use_partial_copy) )
        {
            if( curr_frame_pow > 2.0f * st->prev_wb_bwe_frame_pow )
            {
                scale = root_a_over_b( st->prev_wb_bwe_frame_pow, curr_frame_pow );
                temp = (float) pow(scale, 0.125f);
            }
            else
            {
                scale = 1.0f;
                temp = 1.0f;
            }

            for( j=0; j<8; j++ )
            {
                GainShape[2*j] *= scale;
                GainShape[2*j+1] *= scale;
                for( i=0; i<L_FRAME16k/(4*8); i++ )
                {
                    shaped_wb_excitation[i + j*L_FRAME16k/(4*8)] *= scale;
                }
                scale /= temp;
            }
        }

        st->prev_wb_bwe_frame_pow = curr_frame_pow;

        /* generate 16kHz SHB signal (6 - 8 kHz) from 2kHz signal */
        GenWBSynth( shaped_wb_excitation, error , st->state_lsyn_filt_shb, st->state_lsyn_filt_dwn_shb );

        mvr2r( error + L_FRAME16k - L_SHB_TRANSITION_LENGTH, st->old_tbe_synth, L_SHB_TRANSITION_LENGTH );

        for ( i=0; i<L_FRAME16k; i++ )
        {
            synth[i] = 0.65f * error[i];
        }

        st->last_wb_bwe_ener = 0.0f;
        for( i=0; i<L_FRAME16k; i++ )
        {
            st->last_wb_bwe_ener += synth[i]*synth[i];
        }
        st->last_wb_bwe_ener = (float)sqrt(st->last_wb_bwe_ener/L_FRAME16k);

        if( st->output_Fs == 32000 )  /* 32kHz sampling rate, but only WB output - interpolate */
        {
            Interpolate_allpass_steep( synth, st->mem_resamp_HB, L_FRAME16k, upsampled_synth );
            mvr2r( upsampled_synth, synth, L_FRAME32k );
        }
        else if(st->output_Fs == 48000 )
        {
            interpolate_3_over_1_allpass( synth, L_FRAME16k, upsampled_synth, st->mem_resamp_HB, allpass_poles_3_ov_2 );
            mvr2r( upsampled_synth, synth, L_FRAME48k );
        }
    }
    else
    {
        for ( i = 0; i < LPC_SHB_ORDER_WB; i++ )
        {
            lsf_wb[i] = i/6;
        }
        GainFrame = 0;

        st->prev_wb_bwe_frame_pow = 0.001f;
    }

    /* Update previous frame parameters for FEC */
    mvr2r( lsf_wb, st->lsp_prevfrm, LPC_SHB_ORDER_WB);
    st->GainFrame_prevfrm = GainFrame;

    if( !st->bfi )
    {
        st->GainAttn = 1.0f;
    }

    return;
}

/*-------------------------------------------------------------------*
 * swb_tbe_dec()
 *
 * SWB TBE decoder, 6 - 14 kHz (or 7.5 - 15.5 kHz) band decoding module
 *-------------------------------------------------------------------*/

void swb_tbe_dec(
    Decoder_State *st,              /* i/o: decoder state structure */
    const short coder_type,       /* i  : coding type */
    const float *bwe_exc_extended,/* i  : bandwidth extended excitation  */
    const float voice_factors[],  /* i  : voicing factors                 */
    const float old_syn_12k8_16k[],/*i  : low band synthesis*/
    float *White_exc16k,    /* o  : shaped white excitation for the FB TBE */
    float *synth,           /* o  : SHB synthesis/final synthesis */
    float *pitch_buf
)
{
    short i, j;
    short stemp;
    float shaped_shb_excitation [ L_FRAME16k + L_SHB_LAHEAD ];
    float lsf_shb[LPC_SHB_ORDER], lpc_shb[LPC_SHB_ORDER + 1], GainShape[NUM_SHB_SUBFR], GainFrame;
    float error[L_FRAME32k];
    float ener;
    short is_fractive;
    float prev_pow, curr_pow, scale;
    float curr_frame_pow, temp;
    float GainShapeTemp[NUM_SHB_SUBFR/4], GainGrad0[3], GainGrad1[3], GainGradFEC[4];
    float vf_modified[NB_SUBFR16k];
    float f, inc;
    float GainFrame_prevfrm;
    float tilt_swb_fec;
    float prev_ener_ratio;
    float lsp_shb_1[LPC_SHB_ORDER], lsp_shb_2[LPC_SHB_ORDER], lsp_temp[LPC_SHB_ORDER];
    float lpc_shb_sf[4*(LPC_SHB_ORDER+1)];
    const float *ptr_lsp_interp_coef;
    float shb_ener_sf;
    float shb_res_gshape[NB_SUBFR16k];
    float mixFactors;
    short vind;
    float shb_res_dummy[L_FRAME16k];
    float shaped_shb_excitationTemp[L_FRAME16k];
    float ener_tmp[NUM_SHB_SUBGAINS];
    float GainShape_tmp[NUM_SHB_SUBGAINS];
    float pitch;
    short l_subframe;
    float formant_fac;
    float synth_scale;
    float lsf_diff[LPC_SHB_ORDER], w[LPC_SHB_ORDER];
    float refl[M];
    float tilt_para;

    /* initializations */
    GainFrame = 0.0f;
    mixFactors = 0.0f;
    shb_ener_sf = 0.0f;
    set_f( shaped_shb_excitationTemp, 0.0f, L_FRAME16k );
    st->shb_dtx_count = 0;
    is_fractive = 0;

    /* find tilt */
    calc_tilt_bwe( old_syn_12k8_16k, &tilt_swb_fec, L_FRAME);

    if( st->bfi && st->clas_dec != UNVOICED_CLAS )
    {
        tilt_swb_fec = st->tilt_swb_fec;
    }
    /* WB/SWB bandwidth switching */
    if( (st->tilt_wb > 5 && st->clas_dec == UNVOICED_CLAS) || st->tilt_wb > 10 )
    {
        if( (st->prev_fractive == 0
                && st->prev_enerLH < 2.0f*st->enerLH && st->prev_enerLH > 0.5f*st->enerLH
                && st->prev_enerLL < 2.0f*st->enerLL && st->prev_enerLL > 0.5f*st->enerLL )
                || (st->prev_fractive == 1 && st->prev_enerLH > 3.0f*st->enerLH)
                || (st->enerLL > 1.5f*st->enerLH && st->tilt_wb < 10.0f) )
        {
            is_fractive = 0;
        }
        else
        {
            is_fractive = 1;
        }
    }

    /* WB/SWB bandwidth switching */
    if( st->bws_cnt > 0 )
    {
        f = 1.0f/22.0f;
        inc = 1.0f/22.0f;

        if ( is_fractive == 1 )
        {
            mvr2r(lsf_tab, st->lsp_prevfrm, LPC_SHB_ORDER);
        }
        else
        {
            for (i=0; i<LPC_SHB_ORDER; i++)
            {
                st->lsp_prevfrm[i] = f;
                f += inc;
            }
        }
        if( (st->last_extl != SWB_TBE && st->last_extl != FB_TBE && !(st->prev_enerLH < 2.0f*st->enerLH && st->prev_enerLH > 0.5f*st->enerLH))
                || st->last_core != ACELP_CORE || (st->last_core == ACELP_CORE && labs(st->last_core_brate - st->core_brate) > 3600) || (is_fractive ^ st->prev_fractive) == 1 )
        {
            set_f( GainShape, 0.3536f, NUM_SHB_SUBFR );
        }
        else
        {
            st->prev_GainShape = (st->prev_GainShape > 0.3536f) ? 0.353f : st->prev_GainShape;
            set_f( GainShape, st->prev_GainShape, NUM_SHB_SUBFR );
        }

        /* this never happens */
        mvr2r( st->lsp_prevfrm, lsf_shb, LPC_SHB_ORDER );
        set_f( shb_res_gshape, 0.2f, NB_SUBFR16k );
    }
    else
    {
        if ( st->last_extl != SWB_TBE && st->last_extl != FB_TBE )
        {
            f = 1.0f/22.0f;
            inc = 1.0f/22.0f;
            for (i=0; i<LPC_SHB_ORDER; i++)
            {
                st->lsp_prevfrm[i] = f;
                f += inc;
            }
        }

        if( !st->bfi )
        {
            if( st->use_partial_copy )
            {
                if(st->last_extl != SWB_TBE)
                {
                    st->GainFrame_prevfrm = 0;
                    f = 1.0f/22.0f;
                    inc = 1.0f/22.0f;
                    for (i=0; i<LPC_SHB_ORDER; i++)
                    {
                        st->lsp_prevfrm[i] = f;
                        f += inc;
                    }
                }
                mvr2r( st->lsp_prevfrm, lsf_shb, LPC_SHB_ORDER );
                set_f( GainShape, RECIP_ROOT_EIGHT, NUM_SHB_SUBFR );

                if( st->rf_frame_type == RF_NELP )
                {
                    /* Frame gain */
                    GainFrame = usdequant(st->rf_indx_tbeGainFr, SHB_GAIN_QLOW, SHB_GAIN_QDELTA);
                    GainFrame = (float) pow(10.0, GainFrame);
                    if( st->core == ACELP_CORE && st->last_core == ACELP_CORE && !st->prev_use_partial_copy
                            && st->prev_coder_type == UNVOICED && GainFrame != st->GainFrame_prevfrm && st->next_coder_type != GENERIC && st->last_extl == SWB_TBE )
                    {
                        GainFrame = 0.2f*GainFrame + 0.8f*st->GainFrame_prevfrm;
                    }
                }
                else
                {
                    temp = 0.0f;
                    /* Frame gain */
                    switch (st->rf_indx_tbeGainFr)
                    {
                    case 0:
                        GainFrame = 0.5f;
                        if(st->GainFrame_prevfrm <= 1.25) temp = 0.8f;
                        break;
                    case 1:
                        GainFrame = 2.0f;
                        if(st->GainFrame_prevfrm > 1.25 && st->GainFrame_prevfrm <= 3) temp = 0.8f;
                        break;
                    case 2:
                        GainFrame = 4.0f;
                        if(st->GainFrame_prevfrm > 3 && st->GainFrame_prevfrm <= 6) temp = 0.8f;
                        break;
                    case 3:
                        GainFrame = 8.0f;
                        if(st->GainFrame_prevfrm > 6 && st->GainFrame_prevfrm <= 16) temp = 0.8f;
                        break;
                    default:
                        fprintf(stderr, "RF SWB-TBE gain bits not supported.");
                    }
                    if(st->last_extl == SWB_TBE)
                    {
                        GainFrame = (1 - temp)*GainFrame + temp*(st->GainFrame_prevfrm);
                    }

                    if( st->core == ACELP_CORE && st->last_core == ACELP_CORE )
                    {
                        if( !st->prev_use_partial_copy && st->last_coder_type == VOICED && st->rf_frame_type == RF_GENPRED && GainFrame > 8.0f && GainFrame < 11.67f )
                        {
                            GainFrame *= 0.3f;
                        }
                    }
                }
            }
            else
            {
                /* de-quantization */
                dequantizeSHBparams( st, st->extl, st->extl_brate, lsf_shb, GainShape, &GainFrame, &stemp, &shb_ener_sf, shb_res_gshape, &mixFactors );
            }
        }
        else
        {
            mvr2r( st->lsp_prevfrm, lsf_shb, LPC_SHB_ORDER );

            if( st->codec_mode == MODE1 )
            {
                /* the previous frame gainshape gradient and the gainshape gradient pattern for the current frame */
                for(j=0; j<3; j++)
                {
                    GainGrad0[j] = st->GainShape_Delay[j+1] - st->GainShape_Delay[j];
                    GainGrad1[j] = st->GainShape_Delay[j+5] - st->GainShape_Delay[j+4];
                    GainGradFEC[j+1] = GainGrad0[j]*0.4f + GainGrad1[j]*0.6f;
                }

                /* gradient for the first gainshape */
                if( ( GainGrad1[2] > 2 * GainGrad1[1] && GainGrad1[1] > 2 * GainGrad1[0] ) ||
                        ( GainGrad1[2] < 2 * GainGrad1[1] && GainGrad1[1] < 2 * GainGrad1[0] ) )
                {
                    GainGradFEC[0] = GainGrad1[1] * 0.1f + GainGrad1[2] * 0.9f;
                }
                else
                {
                    GainGradFEC[0] = GainGrad1[0] * 0.2f + GainGrad1[1] * 0.3f + GainGrad1[2] * 0.5f;
                }

                /* get the first gainshape template */
                if( (st->prev_coder_type == UNVOICED || st->last_good == UNVOICED_CLAS) && GainGradFEC[0] > 0 )
                {
                    GainShapeTemp[0] = st->GainShape_Delay[7] + GainGradFEC[0];
                }
                else if( GainGradFEC[0] > 0 )
                {
                    GainShapeTemp[0] = st->GainShape_Delay[7] + GainGradFEC[0] * 0.5f;
                }
                else
                {
                    GainShapeTemp[0] = st->GainShape_Delay[7];
                }

                /*Get the second the third and the fourth gainshape template*/
                if( ( GainGrad1[2] > 10.0f * GainGrad1[1] ) && GainGrad1[1] > 0 )
                {
                    for(i=1; i<NUM_SHB_SUBFR/4; i++)
                    {
                        GainShapeTemp[i] = GainShapeTemp[i-1] + GainGradFEC[i] * 0.8f;
                        GainShapeTemp[i] = max(GainShapeTemp[i], 0.01f);
                    }
                }
                else if( ( GainGrad1[2] > 10.0f * GainGrad1[1] ) && GainGrad1[1] < 0 )
                {
                    for( i=1; i<NUM_SHB_SUBFR/4; i++)
                    {
                        GainShapeTemp[i] = GainShapeTemp[i-1] + GainGradFEC[i] * 0.2f;
                        GainShapeTemp[i] = max(GainShapeTemp[i], 0.01f);
                    }
                }
                else
                {
                    for( i=1; i<NUM_SHB_SUBFR/4; i++)
                    {
                        GainShapeTemp[i] = GainShapeTemp[i-1] + GainGradFEC[i];
                        GainShapeTemp[i] = max(GainShapeTemp[i], 0.01f);
                    }
                }

                /* Get the gainshape and gain frame for the current frame*/
                if( (st->prev_coder_type == UNVOICED || st->last_good == UNVOICED_CLAS) && st->nbLostCmpt == 1 )
                {
                    for( i=0; i<NUM_SHB_SUBFR/4; i++)
                    {
                        for(j=0; j<4; j++)
                        {
                            GainShape[i * 4 + j ] = GainShapeTemp[i] * 1.2f;
                        }
                    }
                    st->GainAttn *= 0.95f;
                }
                else if( st->prev_coder_type == UNVOICED || st->last_good == UNVOICED_CLAS )
                {
                    for( i=0; i<NUM_SHB_SUBFR/4; i++)
                    {
                        for(j=0; j<4; j++)
                        {
                            GainShape[i * 4 + j ] = GainShapeTemp[i];
                        }
                    }
                    st->GainAttn *= 0.95f;
                }
                else if( st->nbLostCmpt > 1 )
                {
                    for( i=0; i<NUM_SHB_SUBFR/4; i++)
                    {
                        for(j=0; j<4; j++)
                        {
                            GainShape[i * 4 + j ] = GainShapeTemp[i] * 0.5f;
                        }
                    }
                    st->GainAttn *= 0.5f;
                }
                else
                {
                    for( i=0; i<NUM_SHB_SUBFR/4; i++)
                    {
                        for(j=0; j<4; j++)
                        {
                            GainShape[i * 4 + j] = GainShapeTemp[i];
                        }
                    }
                    st->GainAttn *= 0.85f;
                }

                GainFrame = st->GainAttn * st->GainFrame_prevfrm;
            }
            else
            {
                for( i=0; i<NUM_SHB_SUBFR/4; i++)
                {
                    for(j=0; j<4; j++)
                    {
                        GainShape[i * 4 + j] = st->cummulative_damping * st->GainShape_Delay[4+i];
                    }
                }

                if( tilt_swb_fec > 8 )
                {
                    if( st->nbLostCmpt == 1 )
                    {
                        GainFrame = 0.6* st->cummulative_damping * st->GainFrame_prevfrm;
                    }
                    else if( st->nbLostCmpt == 2 )
                    {
                        GainFrame = 0.35* st->cummulative_damping * st->GainFrame_prevfrm;
                    }
                    else
                    {
                        GainFrame = 0.2* st->cummulative_damping * st->GainFrame_prevfrm;
                    }
                }
                else
                {
                    GainFrame = st->GainFrame_prevfrm; /* gain locking */
                }
            }

            if( st->total_brate == ACELP_24k40 || st->total_brate == ACELP_32k )
            {
                if( st->codec_mode == MODE1 )
                {
                    scale = st->prev1_shb_ener_sf/sqrt((st->prev2_shb_ener_sf * st->prev3_shb_ener_sf) +0.0001);
                    scale = st->prev_res_shb_gshape * min(scale, 1.0f);
                    if(st->prev2_shb_ener_sf > 2.0f * st->prev1_shb_ener_sf || st->prev3_shb_ener_sf > 2.0f * st->prev2_shb_ener_sf)
                    {
                        shb_ener_sf = 0.5f * scale * st->prev1_shb_ener_sf;
                        if(st->nbLostCmpt > 1)
                        {
                            shb_ener_sf *= 0.5f;
                        }
                    }
                    else
                    {
                        shb_ener_sf = scale * scale * st->prev1_shb_ener_sf;
                    }
                }
                else
                {
                    if( st->prev2_shb_ener_sf > 2.0f * st->prev1_shb_ener_sf || st->prev3_shb_ener_sf > 2.0f * st->prev2_shb_ener_sf )
                    {
                        shb_ener_sf = 0.5f* st->cummulative_damping * st->prev1_shb_ener_sf;
                    }
                    else
                    {
                        shb_ener_sf = st->cummulative_damping * st->prev1_shb_ener_sf;
                    }
                }
            }

            shb_ener_sf = max(shb_ener_sf, 1.0f);
            mixFactors = st->prev_mixFactors;

            if( st->codec_mode == MODE2 )
            {
                set_f( shb_res_gshape, 1.0f, NB_SUBFR16k );
            }
            else
            {
                set_f( shb_res_gshape, 0.2f, NB_SUBFR16k );
            }
        }
    }

    /* get the gainshape delay */
    mvr2r( &st->GainShape_Delay[4], &st->GainShape_Delay[0], NUM_SHB_SUBFR/4 );
    for( i = 0; i<NUM_SHB_SUBFR/4; i++ )
    {
        st->GainShape_Delay[i+4] = GainShape[i*4];
    }

    mvr2r( voice_factors, vf_modified, NB_SUBFR16k );
    if( coder_type == VOICED || mean(voice_factors, 4) > 0.4f )
    {
        for( i = 1; i < NB_SUBFR; i++ )
        {
            vf_modified[i] = 0.8f * voice_factors[i] + 0.2f * voice_factors[i-1];
        }

        if(st->L_frame != L_FRAME )
        {
            vf_modified[4] = 0.8f * voice_factors[4] + 0.2f * voice_factors[3];
        }

    }

    if(st->use_partial_copy && st->nelp_mode_dec)
    {
        set_f( vf_modified, 0.0f, NB_SUBFR16k );
    }

    /* SHB LSF from current frame; and convert to LSP for interpolation */
    lsf2lsp( lsf_shb, lsp_shb_2, LPC_SHB_ORDER, 1.0f );

    if( st->last_extl == SWB_TBE || st->last_extl == FB_TBE )
    {
        /* SHB LSP values from prev. frame for interpolation */
        mvr2r( st->swb_lsp_prev_interp, lsp_shb_1, LPC_SHB_ORDER );
    }
    else
    {
        /* Use current frame's LSPs; in effect no interpolation */
        mvr2r( lsp_shb_2, lsp_shb_1, LPC_SHB_ORDER );
    }
    if( st->bws_cnt == 0 && st->bws_cnt1 == 0 && st->prev_use_partial_copy == 0 && st->use_partial_copy == 0)
    {
        lsf_diff[0] = lsf_diff[LPC_SHB_ORDER-1] = 0.5f;
        for(i=1; i<(LPC_SHB_ORDER-1); i++)
        {
            lsf_diff[i] = lsf_shb[i] - lsf_shb[i-1];
        }

        a2rc (st->cur_sub_Aq+1, refl, (short) M);
        tilt_para = 6.6956f * (1.0f + refl[0]) * (1.0f + refl[0]) - 3.8714f * (1.0f + refl[0]) + 1.3041f;
        if( st->last_extl != SWB_TBE && st->last_extl != FB_TBE )
        {
            for( i=1; i<LPC_SHB_ORDER-1; i++ )
            {
                st->prev_lsf_diff[i-1] = 0.5f*lsf_diff[i];
            }
        }

        if( st->total_brate <= ACELP_16k40 )
        {
            if(!(st->prev_tilt_para > 5.0f && (coder_type == TRANSITION || tilt_para < 1.0f)) && !(((st->prev_tilt_para < 3.0f && st->prev_coder_type >= VOICED)) && tilt_para > 5.0f))
            {
                for( i = 1; i < (LPC_SHB_ORDER-1); i++ )
                {
                    if( lsf_diff[i] < 0 || st->prev_lsf_diff[i-1] <= 0 ) /* safety check in case of bit errors */
                    {
                        w[i] = 0;
                        st->BER_detect = 1;
                    }
                    else
                    {
                        w[i] = (lsf_diff[i] < st->prev_lsf_diff[i-1]) ? min(max(0.8f*lsf_diff[i]/st->prev_lsf_diff[i-1], 0.5f), 1.0f) : min(max(0.8f*st->prev_lsf_diff[i-1]/lsf_diff[i], 0.5f), 1.0f);
                    }
                }
                w[0] = w[1];
                w[LPC_SHB_ORDER-1] = w[LPC_SHB_ORDER-2];

                for( i = 0; i < LPC_SHB_ORDER; i++ )
                {
                    lsp_temp[i] = lsp_shb_1[i]*(1.0f-w[i]) + lsp_shb_2[i]*w[i];
                }
            }
            else
            {
                mvr2r( lsp_shb_2, lsp_temp, LPC_SHB_ORDER );
            }

            /* convert from lsp to lsf */
            lsp2lsf( lsp_temp, lsf_shb, LPC_SHB_ORDER, 1.0f );
        }
        mvr2r( lsf_diff+1, st->prev_lsf_diff, LPC_SHB_ORDER-2 );
        st->prev_tilt_para = tilt_para;
    }
    else
    {
        mvr2r(lsp_shb_2, lsp_temp, LPC_SHB_ORDER);
    }

    if( st->total_brate == ACELP_24k40 || st->total_brate == ACELP_32k )
    {
        /* SHB LSP interpolation */
        ptr_lsp_interp_coef = interpol_frac_shb;
        for( j = 0; j < 4; j++ )
        {
            for( i = 0; i < LPC_SHB_ORDER; i++ )
            {
                lsp_temp[i] = lsp_shb_1[i]*(*ptr_lsp_interp_coef) + lsp_shb_2[i]*(*(ptr_lsp_interp_coef+1));
            }
            ptr_lsp_interp_coef += 2;

            /* convert from lsp to lsf */
            lsp2lsf( lsp_temp, lsp_temp, LPC_SHB_ORDER, 1.0f );

            /* convert lsf to lpc for SHB synthesis */
            lsp2a( lpc_shb_sf+j*(LPC_SHB_ORDER+1), lsp_temp, LPC_SHB_ORDER );
            lpc_shb_sf[j*(LPC_SHB_ORDER+1)] = 1.0f;
        }
    }

    /* Save the SWB LSP values from current frame for interpolation */
    mvr2r( lsp_shb_2, st->swb_lsp_prev_interp, LPC_SHB_ORDER );

    /* save the shb_ener and mixFactor values */
    st->prev3_shb_ener_sf = st->prev2_shb_ener_sf;
    st->prev2_shb_ener_sf = st->prev1_shb_ener_sf;
    st->prev1_shb_ener_sf = shb_ener_sf;
    st->prev_res_shb_gshape = shb_res_gshape[4];
    st->prev_mixFactors = mixFactors;

    /* SWB CNG/DTX - update memories */
    mvr2r( st->lsp_shb_prev, st->lsp_shb_prev_prev, LPC_SHB_ORDER );
    mvr2r( lsf_shb, st->lsp_shb_prev, LPC_SHB_ORDER );

    /* convert LSPs back into LP coeffs */
    lsp2a( lpc_shb, lsf_shb, LPC_SHB_ORDER );
    lpc_shb[0] = 1.0;
    vind = (short)(mixFactors*8.0f);

    /* Determine formant PF strength */
    formant_fac = swb_formant_fac( lpc_shb[1], &st->tilt_mem );
    if(st->total_brate > ACELP_32k)
    {
        for( j = 0; j < 4; j++ )
        {
            mvr2r(lpc_shb, &lpc_shb_sf[j*(LPC_SHB_ORDER+1)], LPC_SHB_ORDER+1);
        }
    }

    /* From low band excitation, generate highband excitation */
    mvr2r( st->state_syn_shbexc, shaped_shb_excitation, L_SHB_LAHEAD);
    GenShapedSHBExcitation( shaped_shb_excitation + L_SHB_LAHEAD, lpc_shb, White_exc16k, st->mem_csfilt, st->mem_genSHBexc_filt_down_shb, st->state_lpc_syn,
                            coder_type, bwe_exc_extended, st->bwe_seed, vf_modified, st->extl, &(st->tbe_demph), &(st->tbe_premph), lpc_shb_sf,
                            &shb_ener_sf, shb_res_gshape, shb_res_dummy, &vind, formant_fac, st->fb_state_lpc_syn, &(st->fb_tbe_demph), st->total_brate , st->prev_bfi);

    for( i=0; i<L_FRAME16k; i+=L_SUBFR16k )
    {
        /* TD BWE post-processing */
        PostShortTerm( &shaped_shb_excitation[L_SHB_LAHEAD+i], lpc_shb, &shaped_shb_excitationTemp[i], st->mem_stp_swb,
                       st->ptr_mem_stp_swb, &(st->gain_prec_swb), st->mem_zero_swb, formant_fac );
    }

    mvr2r( shaped_shb_excitationTemp, &shaped_shb_excitation[L_SHB_LAHEAD], L_FRAME16k );
    prev_pow = sum2_f( shaped_shb_excitation, L_SHB_LAHEAD + 10 );
    curr_pow = sum2_f( shaped_shb_excitation + L_SHB_LAHEAD + 10, L_SHB_LAHEAD + 10 );

    if( voice_factors[0] > 0.75f )
    {
        curr_pow *= 0.25;
    }

    if( prev_pow == 0 )
    {
        scale = 0;
    }
    else
    {
        scale = sqrt(curr_pow/prev_pow);
    }
    for( i=0; i<L_SHB_LAHEAD; i++ )
    {
        shaped_shb_excitation[i] *= scale;
    }
    for(   ; i<L_SHB_LAHEAD + 10 ; i++)
    {
        temp = (i-19)/10.0f;
        shaped_shb_excitation[i] *= (temp*1.0f + (1.0f-temp)*scale);
    }

    /* Update SHB excitation */
    mvr2r( shaped_shb_excitation + L_FRAME16k, st->state_syn_shbexc, L_SHB_LAHEAD );
    l_subframe = L_FRAME16k/NUM_SHB_SUBGAINS;
    ener = EPSILON;
    for(i=0; i<NUM_SHB_SUBGAINS; i++)
    {
        ener_tmp[i] = EPSILON;
        for(j=0; j<l_subframe; j++)
        {
            ener_tmp[i] += shaped_shb_excitation[i*l_subframe+j] * shaped_shb_excitation[i*l_subframe+j] * 0.0125f;
        }
        ener_tmp[i] = sqrt(ener_tmp[i]);
        ener += ener_tmp[i];
    }
    ener /= NUM_SHB_SUBGAINS;

    /* WB/SWB bandwidth switching */
    if( st->bws_cnt > 0 )
    {
        ener *= 0.35f;

        if( st->tilt_swb > 8 )
        {
            st->prev_fractive = 1;
        }

        if( is_fractive == 0 )
        {
            if( st->tilt_wb > 1.0 )
            {
                st->tilt_wb = 1.0f;
            }
            else if( st->tilt_wb < 0.5 )
            {
                st->tilt_wb = 0.5f;
            }

            if( st->prev_fractive == 1 && st->tilt_wb > 0.5 )
            {
                st->tilt_wb = 0.5f;
            }
        }
        else
        {
            if ( st->tilt_wb > 4)
            {
                if ( st->prev_fractive == 0)
                {
                    st->tilt_wb = 4;
                }
                else
                {
                    st->tilt_wb = 8;
                }
            }
            else
            {
                st->tilt_wb *= 2;
            }
        }

        if( ener != 0 )
        {
            if( ener*st->tilt_wb > st->enerLH )
            {
                st->tilt_wb = 0.5f*st->enerLH/ener;
            }
            else if( ener*st->tilt_wb < 0.05f*st->enerLH && is_fractive == 1 )
            {
                st->tilt_wb = 0.25f*st->enerLH/ener;
            }

            GainFrame_prevfrm = st->prev_ener_shb/ener;
        }
        else
        {
            GainFrame_prevfrm = 0;
        }

        if ( is_fractive == 1)
        {
            GainFrame = 8.0f*st->tilt_wb;
        }
        else
        {
            GainFrame = 2.0f*st->tilt_wb;
        }

        if ( (is_fractive & st->prev_fractive) == 1 && GainFrame > GainFrame_prevfrm )
        {
            GainFrame = 0.2f*GainFrame + 0.8f*GainFrame_prevfrm;
        }
        else
        {
            if ( (st->prev_enerLH < 2.0f*st->enerLH && st->prev_enerLH > 0.5f*st->enerLH )
                    && (st->prev_enerLL < 2.0f*st->enerLL && st->prev_enerLL > 0.5f*st->enerLL)
                    && (is_fractive ^ st->prev_fractive) == 0)
            {
                GainFrame = 0.5f*GainFrame + 0.5f*GainFrame_prevfrm;
            }
            else
            {
                if ( is_fractive == 0 && st->prev_fractive == 1 )
                {
                    GainFrame = (1.0f-0.1f*GainFrame)*GainFrame + 0.1f*GainFrame*GainFrame_prevfrm;
                }
                else
                {
                    GainFrame = 0.5f*GainFrame + 0.5f*GainFrame_prevfrm;
                }
            }
        }

        GainFrame *= ((float)N_WS2N_FRAMES - (float)st->bws_cnt) / (float)N_WS2N_FRAMES;
    }
    else
    {
        if( st->bws_cnt1 > 0 )
        {
            GainFrame *= (float)st->bws_cnt1 / (float)N_WS2N_FRAMES;
        }

        if( st->nbLostCmpt == 1 )
        {
            prev_ener_ratio = st->prev_ener_shb/ener;

            if( st->clas_dec != UNVOICED_CLAS && st->clas_dec != UNVOICED_TRANSITION &&st->tilt_swb_fec < 8.0 &&
                    ((st->enerLL > 0.5f*st->prev_enerLL && st->enerLL < 2.0f*st->prev_enerLL)|| (st->enerLH > 0.5f*st->prev_enerLH && st->enerLH < 2.0f*st->prev_enerLH)))
            {
                if( prev_ener_ratio > 4.0f * GainFrame )
                {
                    GainFrame = 0.4f * prev_ener_ratio + 0.6f * GainFrame;
                }
                else if( prev_ener_ratio > 2.0f * GainFrame )
                {
                    GainFrame = 0.8f * prev_ener_ratio + 0.2f * GainFrame;
                }
                else
                {
                    GainFrame = 0.2f * prev_ener_ratio + 0.8f * GainFrame;
                }

                if( tilt_swb_fec > st->tilt_swb_fec )
                {
                    GainFrame *=  st->tilt_swb_fec > 0 ? (min(5.0f,tilt_swb_fec/st->tilt_swb_fec)) : 1.0f;
                }

            }
            else if( (st->clas_dec != UNVOICED_CLAS || st->tilt_swb_fec > 8.0) && prev_ener_ratio > 4.0f * GainFrame &&
                     (st->enerLL > 0.5f*st->prev_enerLL ||st->enerLH > 0.5f*st->prev_enerLH) )
            {
                GainFrame = 0.2f * prev_ener_ratio + 0.8f * GainFrame;
            }
        }
        else if( st->nbLostCmpt > 1 )
        {
            prev_ener_ratio = st->prev_ener_shb/ener;
            if((prev_ener_ratio > 4.0 * GainFrame) && ((st->codec_mode == MODE1 && st->enerLL > st->prev_enerLL && st->enerLH > st->prev_enerLH) || st->codec_mode == MODE2))
            {
                if( tilt_swb_fec > 10.0f && st->tilt_swb_fec >10.0f )
                {
                    GainFrame =  min((prev_ener_ratio *0.8f + GainFrame * 0.2f),4.0f * GainFrame);
                }
                else
                {
                    GainFrame =  min((prev_ener_ratio *0.5f + GainFrame * 0.5f),4.0f * GainFrame);
                }
            }
            else if((prev_ener_ratio >  GainFrame) &&((st->codec_mode == MODE1 && st->enerLL > st->prev_enerLL && st->enerLH > st->prev_enerLH) || st->codec_mode == MODE2))
            {
                if( tilt_swb_fec > 10.0f && st->tilt_swb_fec >10.0f )
                {
                    GainFrame =  0.5f * prev_ener_ratio + 0.5f * GainFrame;
                }
                else
                {
                    GainFrame =  0.2f * prev_ener_ratio + 0.8f * GainFrame;
                }
            }
        }
    }

    st->prev_fractive = is_fractive;

    /* Adjust the subframe and frame gain of the synthesized shb signal */
    /* Scale the shaped excitation */
    if( st->L_frame == L_FRAME )
    {
        pitch = 0.25f*sum_f(pitch_buf, 4);
    }
    else
    {
        pitch = 0.2f*sum_f(pitch_buf, 5);
    }

    if( ((st->total_brate >= ACELP_24k40 && st->prev_coder_type == coder_type && coder_type != UNVOICED)
            || (st->total_brate <= ACELP_16k40 && (st->prev_coder_type == coder_type || (st->prev_coder_type == VOICED && coder_type == GENERIC) || (st->prev_coder_type == GENERIC && coder_type == VOICED))))
            && pitch > 70 && st->extl < FB_TBE)
    {
        for( i=0; i<NUM_SHB_SUBGAINS; i++ )
        {
            GainShape_tmp[i] = GainShape[i*4];
        }

        for( i=0; i<NUM_SHB_SUBGAINS; i++ )
        {
            if( ener_tmp[i]*GainShape_tmp[i] > st->prev_ener*st->prev_GainShape )
            {
                GainShape_tmp[i] = 0.5f*(st->prev_ener*st->prev_GainShape/ener_tmp[i] + GainShape_tmp[i]);
            }
            st->prev_ener = ener_tmp[i];
            st->prev_GainShape = GainShape_tmp[i];
        }
        for( i=0; i<NUM_SHB_SUBFR; i++ )
        {
            GainShape[i] = GainShape_tmp[i*NUM_SHB_SUBGAINS/NUM_SHB_SUBFR];
        }
    }

    ScaleShapedSHB( SHB_OVERLAP_LEN, shaped_shb_excitation, st->syn_overlap,GainShape,GainFrame,window_shb,subwin_shb );

    curr_frame_pow = sum2_f( shaped_shb_excitation, L_FRAME16k ) + 0.001f;

    if( !st->bfi && (st->prev_bfi || st->prev_use_partial_copy ) )
    {

        if( ( curr_frame_pow > 2.0f * st->prev_swb_bwe_frame_pow ) &&
                ( curr_frame_pow < 30.0f * st->prev_swb_bwe_frame_pow ) &&
                st->prev_coder_type == UNVOICED )
        {
            scale = sqrt( st->prev_swb_bwe_frame_pow/curr_frame_pow );
            if ((curr_frame_pow)==0) scale = 0;

            temp = (float)pow( scale, 0.125f );
        }
        else
        {
            scale = 1.0f;
            temp = 1.0f;
        }

        for( j=0; j<8; j++ )
        {
            GainShape[2*j] *= scale;
            GainShape[2*j+1] *= scale;
            for( i=0; i<L_FRAME16k/8; i++ )
            {
                shaped_shb_excitation[i + j*L_FRAME16k/8] *= scale;
            }

            scale /= temp;
        }
    }

    /* adjust the FEC frame energy */
    if( st->bfi )
    {
        scale = 1.0f;
        temp = 1.0f;
        if( st->nbLostCmpt == 1 )
        {
            if( curr_frame_pow > st->prev_swb_bwe_frame_pow &&
                    st->prev_coder_type != UNVOICED &&
                    st->last_good != UNVOICED_CLAS )
            {
                scale = sqrt( st->prev_swb_bwe_frame_pow/curr_frame_pow );
                if ((curr_frame_pow)==0) scale = 0;
                temp = (float) pow( scale, 0.125f );
            }
            else if( curr_frame_pow < 0.5f *st->prev_swb_bwe_frame_pow && st->nbLostCmpt == 1 &&
                     (st->enerLL > 0.5 * st->prev_enerLL || st->enerLH > 0.5 *st->prev_enerLH) &&
                     (st->prev_coder_type == UNVOICED || st->last_good == UNVOICED_CLAS || st->tilt_swb_fec > 5.0f) )
            {
                scale = sqrt(st->prev_swb_bwe_frame_pow / curr_frame_pow);
                if ((curr_frame_pow)==0) scale = 0;
                temp = (float) pow(scale, 0.125f);
            }
        }
        else if ( st->nbLostCmpt > 1 )
        {
            if( curr_frame_pow > st->prev_swb_bwe_frame_pow )
            {
                scale = sqrt( st->prev_swb_bwe_frame_pow / curr_frame_pow );
                if ((curr_frame_pow)==0) scale = 0;
                temp = (float) pow( scale, 0.125f );
            }
            else if( curr_frame_pow < 0.5f *st->prev_swb_bwe_frame_pow &&
                     (st->enerLL > 0.5 * st->prev_enerLL || st->enerLH > 0.5 *st->prev_enerLH) &&
                     (st->prev_coder_type == UNVOICED || st->last_good == UNVOICED_CLAS || st->tilt_swb_fec > 5.0f) )
            {
                scale = min(2.0f,sqrt(st->prev_swb_bwe_frame_pow/curr_frame_pow));
                if ((curr_frame_pow)==0) scale = 0;
                temp = (float) pow(scale, 0.125f);
            }
        }

        for( j=0; j<8; j++ )
        {
            GainShape[2 * j] *= scale;
            GainShape[2 * j + 1] *= scale;
            for( i=0; i<L_FRAME16k/8; i++ )
            {
                shaped_shb_excitation[i + j * L_FRAME16k/8] *= scale;
            }

            scale /= temp;
        }
    }

    st->prev_swb_bwe_frame_pow = curr_frame_pow;

    st->prev_ener_shb = EPSILON;
    for( i=0; i<L_FRAME16k; i++ )
    {
        st->prev_ener_shb += shaped_shb_excitation[i] * shaped_shb_excitation[i];
    }
    st->prev_ener_shb = (float)sqrt(st->prev_ener_shb/L_FRAME16k);

    for(i=0; i<SWB_FENV; i++)
    {
        st->prev_SWB_fenv[i] = (float)sqrt(curr_frame_pow/L_FRAME16k);
    }

    /* generate 32kHz SHB synthesis from 12.8(16)kHz signal */
    GenSHBSynth( shaped_shb_excitation, error, st->genSHBsynth_Hilbert_Mem,
                 st->genSHBsynth_state_lsyn_filt_shb_local, st->L_frame, &(st->syn_dm_phase) );

    mvr2r( error + L_FRAME32k - L_SHB_TRANSITION_LENGTH, st->old_tbe_synth, L_SHB_TRANSITION_LENGTH );

    /* resample SHB synthesis (if needed) and scale down */
    synth_scale = (st->codec_mode == MODE1) ? 0.9f : 1.f;

    if( st->output_Fs == 48000 )
    {
        if ( st->extl == FB_TBE)
        {
            for( i=0; i<L_FRAME16k; i++ )
            {
                White_exc16k[i] *= GainFrame * GainShape[NUM_SHB_SUBFR*i/L_FRAME16k];
            }
        }

        for( i=0; i<L_FRAME32k; i++ )
        {
            error[i] *= synth_scale;
        }

        interpolate_3_over_2_allpass( error, L_FRAME32k, synth, st->int_3_over_2_tbemem_dec, allpass_poles_3_ov_2 );
    }

    else if( st->output_Fs == 32000 )
    {
        for( i=0; i<L_FRAME32k; i++ )
        {
            synth[i] = synth_scale * error[i];
        }
    }
    else if( st->output_Fs == 16000 )
    {
        for( i=0; i<L_FRAME32k; i++ )
        {
            error[i] *= synth_scale;
        }

        Decimate_allpass_steep( error, st->mem_resamp_HB_32k, L_FRAME32k, synth );
    }

    /* Update previous frame parameters for FEC */
    mvr2r( lsf_shb, st->lsp_prevfrm, LPC_SHB_ORDER );
    if( st->codec_mode == MODE1 )
    {
        st->GainFrame_prevfrm = GainFrame;
        st->tilt_swb_fec = tilt_swb_fec;

        if( !st->bfi )
        {
            st->GainAttn = 1.0f;
        }
    }
    else
    {
        if( !st->bfi )
        {
            st->tilt_swb_fec = tilt_swb_fec;
            st->GainFrame_prevfrm = GainFrame; /* gain locking on lost frame */
            st->GainAttn = 1.0f;
        }
    }
    st->prev_ener = ener_tmp[NUM_SHB_SUBGAINS-1];
    st->prev_GainShape = GainShape[NUM_SHB_SUBFR-1];

    return;
}

/*-------------------------------------------------------------------*
 * Dequant_lower_LSF()
 *
 * Dequantized the lower LSFs
 *-------------------------------------------------------------------*/

static void Dequant_lower_LSF(
    const short lsf_idx[],             /* i  : LSF indices */
    float lsf_q[]                      /* o  : Quantized LSFs */
)
{
    short i;

    lsf_q[0] = lsf_q_cb[0][lsf_idx[0]];
    for (i = 1; i < NUM_Q_LSF; i++)
    {
        lsf_q[i] = lsf_q_cb[i][lsf_idx[i]] + lsf_q[i-1];
    }

    return;
}

/*-------------------------------------------------------------------*
 * Map_higher_LSF()
 *
 * Map the higher LSFs from the lower LSFs
 *-------------------------------------------------------------------*/

static void Map_higher_LSF(
    float lsf_q[],               /* i/o : Quantized lower LSFs */
    const float m,                     /* i   : Mirroring point */
    const float grid_in[]              /* i   : Input LSF smoohthing grid */
)
{
    float lsf_map[NUM_MAP_LSF];
    float grid[NUM_MAP_LSF];
    float last_q_lsf;
    float lsf_smooth[NUM_MAP_LSF];
    float offset;
    short i;
    float scale;

    for (i = 0; i < NUM_MAP_LSF; i++)
    {
        lsf_map[i] = 2*m - lsf_q[NUM_MAP_LSF - 1 - i];
    }

    if (m > MAX_LSF/2)
    {
        offset = lsf_map[0];
        scale = (MAX_LSF - m)/m;
        for (i = 0; i < NUM_MAP_LSF; i++)
        {
            lsf_map[i] = (lsf_map[i] - offset)*scale + offset;
        }
    }

    last_q_lsf = lsf_q[NUM_Q_LSF - 1];
    scale = MAX_LSF - last_q_lsf;

    for (i = 0; i < NUM_MAP_LSF; i++)
    {
        grid[i] = grid_in[i]*scale + last_q_lsf;
    }

    for (i = 0; i < NUM_MAP_LSF; i++)
    {
        lsf_smooth[i] = (1 - grid_smoothing[i])*lsf_map[i] + grid_smoothing[i]*grid[i];
    }

    for (i = 0; i < NUM_MAP_LSF; i++)
    {
        lsf_q[NUM_Q_LSF + i] = lsf_smooth[i];
    }

    return;
}

/*-------------------------------------------------------------------*
 * Map_higher_LSF()
 *
 * Map the higher LSFs from the lower LSFs
 *-------------------------------------------------------------------*/

static void Dequant_mirror_point(
    const float lsf_q[],               /* i/o : Quantized lower LSFs */
    const short m_idx,                 /* i   : Mirror point index */
    float *m                     /* i   : Mirroring point */
)
{
    *m = mirror_point_q_cb[m_idx] + lsf_q[NUM_Q_LSF - 1];

    return;
}

/*-------------------------------------------------------------------*
 * dequantizeSHBparams()
 *
 * Dequantize super highband spectral envolope, temporal gains and frame gain
 *-------------------------------------------------------------------*/

static void dequantizeSHBparams(
    Decoder_State *st,                /* i/o: decoder state structure   */
    const short extl,               /* i  : extension layer                         */
    long  extl_brate,         /* i  : extensiuon layer bitrate                */
    float *Q_lsf,             /* o  : SHB LSF from de-quantization            */
    float *Q_subgain,         /* o  : SHB subframe gains from de-quantization */
    float *Q_framegrain,      /* o  : SHB frame gain from de-quantization     */
    short *uv_flag,           /* o  : unvoiced flag                           */
    float *Q_shb_ener_sf,
    float *Q_shb_res_gshape,
    float *Q_mixFactors
)
{
    short i, j, idxLSF, idxSubGain, idxFrameGain;
    float Q_combined_gains[NUM_SHB_SUBFR/4];
    float lsf_q[LPC_SHB_ORDER];
    short lsf_idx[NUM_Q_LSF];
    short m_idx, grid_idx;
    float m;
    short idx_shb_fr_gain, idx_res_gs[5], idx_mixFac;

    /* LSFs */
    if( extl == WB_TBE )
    {
        if ( extl_brate == WB_TBE_0k35 )
        {
            idxFrameGain = st->gFrame_WB;
            idxLSF = st->lsf_WB;

            mvr2r( lbr_wb_bwe_lsfvq_cbook_2bit + idxLSF*LPC_SHB_ORDER_LBR_WB, Q_lsf, LPC_SHB_ORDER_LBR_WB );
            set_f( Q_subgain, RECIP_ROOT_EIGHT, NUM_SHB_SUBFR/2 );
            mvr2r( SHBCB_FrameGain16 + idxFrameGain, Q_framegrain, 1 );
        }
        else
        {
            /* read the information about UNVOICED frame */
            *uv_flag = (short)get_next_indice( st, 1 );

            idxSubGain = (short)get_next_indice( st, NUM_BITS_SHB_SUBGAINS );
            idxFrameGain = (short)get_next_indice( st, NUM_BITS_SHB_FrameGain );
            idxLSF = (short)get_next_indice( st, NUM_BITS_WB_LSF );

            mvr2r( wb_bwe_lsfvq_cbook_8bit + idxLSF*LPC_SHB_ORDER_WB, Q_lsf, LPC_SHB_ORDER_WB );
            mvr2r( HBCB_SubGain5bit + idxSubGain * NUM_SHB_SUBFR/4, Q_combined_gains, NUM_SHB_SUBFR/4 );

            for( i=0; i<NUM_SHB_SUBFR/4; i++ )
            {
                Q_combined_gains[i] = (float) pow(10.0f, Q_combined_gains[i] / 20.0f);
            }

            for( i=0; i<NUM_SHB_SUBFR/2; i+=2 )
            {
                Q_subgain[i] = Q_combined_gains[i/2];
                Q_subgain[i+1] = Q_combined_gains[i/2];
            }

            /* frame gain */
            mvr2r( SHBCB_FrameGain64 + idxFrameGain, Q_framegrain, 1 );
        }
    }
    else
    {
        if( st->codec_mode == MODE2 )
        {
            idxSubGain = st->idxSubGains;
            idxFrameGain = st->idxFrameGain;
        }
        else
        {
            idxSubGain = (short)get_next_indice( st, NUM_BITS_SHB_SUBGAINS );
            idxFrameGain = (short)get_next_indice( st, NUM_BITS_SHB_FRAMEGAIN );
        }

        if( st->total_brate == ACELP_24k40 || st->total_brate == ACELP_32k )
        {
            if( st->codec_mode == MODE2 )
            {
                idx_shb_fr_gain = st->idx_shb_fr_gain;
            }
            else
            {
                idx_shb_fr_gain = (short)get_next_indice( st, NUM_BITS_SHB_ENER_SF );
            }
            *Q_shb_ener_sf = usdequant(idx_shb_fr_gain, 0, 0.042f);
            *Q_shb_ener_sf = (float)pow(10.0, *Q_shb_ener_sf );

            for( i=0; i<5; i++ )
            {
                if( st->codec_mode == MODE2 )
                {
                    idx_res_gs[i] = st->idx_res_gs[i];
                }
                else
                {
                    idx_res_gs[i] = (short)get_next_indice( st, NUM_BITS_SHB_RES_GS );
                }
                Q_shb_res_gshape[i] = usdequant(idx_res_gs[i], 0.125f, 0.125f);
            }

            if( st->codec_mode == MODE2 )
            {
                idx_mixFac = st->idx_mixFac;
            }
            else
            {
                idx_mixFac = (short)get_next_indice( st, NUM_BITS_SHB_VF );
            }
            *Q_mixFactors = usdequant(idx_mixFac, 0.125f, 0.125f);
        }
        else
        {
            *Q_shb_ener_sf = 0;
            *Q_mixFactors = 0;
            set_f(Q_shb_res_gshape, 0, 5);
        }


        if( st->rf_flag == 0 && !((st->total_brate == ACELP_9k60) || ( (st->total_brate == 0) && ((st->last_total_brate == ACELP_9k60) || (st->last_total_brate == ACELP_13k20 && st->rf_flag_last)) )) )
        {
            /* LSFs */
            if( extl_brate == SWB_TBE_1k6 || extl_brate == FB_TBE_1k8 || extl_brate == SWB_TBE_2k8 || extl_brate == FB_TBE_3k0 )
            {
                for( i = 0; i < NUM_Q_LSF; i++ )
                {
                    if( st->codec_mode == MODE2 )
                    {
                        lsf_idx[i] = st->lsf_idx[i];
                    }
                    else
                    {
                        lsf_idx[i] = (short)get_next_indice(st, lsf_q_num_bits[i]);
                    }
                }
            }
            Dequant_lower_LSF( lsf_idx, lsf_q );
            if( st->codec_mode == MODE2 )
            {
                m_idx = st->m_idx;
                grid_idx = st->grid_idx;
            }
            else
            {
                m_idx = (short)get_next_indice( st, MIRROR_POINT_BITS );
                grid_idx = (short)get_next_indice( st, NUM_LSF_GRID_BITS );
            }

            Dequant_mirror_point(lsf_q, m_idx, &m);

            /* safety check in case of bit errors */
            if( m > MAX_LSF )
            {
                st->BER_detect = 1;
                m = MAX_LSF;
            }

            Map_higher_LSF(lsf_q, m, lsf_grid[grid_idx]);

            for( i = 0; i < LPC_SHB_ORDER; i++ )
            {
                /* safety check in case of bit errors */
                if ( lsf_q[LPC_SHB_ORDER - 1 - i] > MAX_LSF )
                {
                    st->BER_detect = 1;
                    lsf_q[LPC_SHB_ORDER - 1 - i] = MAX_LSF;
                }
                Q_lsf[i] = 0.5f - lsf_q[LPC_SHB_ORDER - 1 - i];
            }
        }
        else
        {
            set_s(lsf_idx, 0, 5);
            mvs2s(st->lsf_idx, lsf_idx, 5);
            grid_idx = 0;
            m_idx = 0;
            mvr2r( swb_tbe_lsfvq_cbook_8b + lsf_idx[0]*LPC_SHB_ORDER, Q_lsf, LPC_SHB_ORDER );
        }
        space_lsfs( Q_lsf, LPC_SHB_ORDER );

        /* Dequantize subgain indices */
        j =  idxSubGain*NUM_SHB_SUBGAINS;
        for ( i = 0; i < NUM_SHB_SUBGAINS; i++)
        {
            Q_subgain[i] = (float) pow(10.0, SHBCB_SubGain5bit[j++]);
        }

        for (i=NUM_SHB_SUBFR-1; i>=0; i--)
        {
            Q_subgain[i] = Q_subgain[i*NUM_SHB_SUBGAINS/NUM_SHB_SUBFR];
        }

        /* Frame gain */
        *Q_framegrain = usdequant(idxFrameGain, SHB_GAIN_QLOW, SHB_GAIN_QDELTA);
        *Q_framegrain = (float) pow(10.0, *Q_framegrain);
    }

    return;
}
/*-------------------------------------------------------------------*
 * fb_tbe_dec()
 *
 * FB TBE decoder, 14(resp. 15.5) - 20 kHz band decoding module
 *-------------------------------------------------------------------*/
void fb_tbe_dec(
    Decoder_State *st,                /* i/o: encoder state structure                 */
    const float fb_exc[],           /* i  : FB excitation from the SWB part         */
    float *hb_synth           /* o  : high-band synthesis                     */
)

{
    short i;
    float ratio = 0, fb_exc_energy = 0;
    float fb_synth[L_FRAME48k];

    /* decode FB slope information */
    if ( st->extl == FB_TBE && !st->bfi )
    {
        if( st->codec_mode == MODE2 )
        {
            i = st->idxGain;
        }
        else
        {
            i = (short)get_next_indice( st, 4 );
        }
        ratio = (float)(1 << i);
    }
    else if ( st->extl == FB_TBE && st->bfi )
    {
        ratio = st->prev_fbbwe_ratio;
    }

    fb_exc_energy = sum2_f(fb_exc,L_FRAME16k);

    /* FB TBE synthesis */
    synthesise_fb_high_band( fb_exc,fb_synth,fb_exc_energy,ratio, st->L_frame, st->bfi, &(st->prev_fbbwe_ratio), st->fbbwe_hpf_mem );

    /* add the fb_synth component to the hb_synth component */
    v_add( hb_synth, fb_synth, hb_synth, L_FRAME48k );

    return;
}


/*---------------------------------------------------------------------*
 * tbe_read_bitstream()
 *
 * Read TBE bitstream and populate the parameters for TD-BWE decoder.
 *---------------------------------------------------------------------*/
void tbe_read_bitstream(
    Decoder_State *st
)
{
    short i;

    if ( (st->rf_flag || st->total_brate == ACELP_9k60) && st->bwidth == WB)
    {
        /* WB LSF */
        st->lsf_WB = get_next_indice(st, NUM_BITS_LBR_WB_LSF);

        /* WB frame gain */
        st->gFrame_WB = get_next_indice(st, NUM_BITS_SHB_FrameGain_LBR_WB);
    }
    else if( st->total_brate >= ACELP_9k60 && st->total_brate <= ACELP_32k && (st->bwidth == SWB || st->bwidth == FB) )
    {
        if(st->rf_flag == 0 && st->total_brate > ACELP_9k60)
        {
            for (i=0; i<NUM_Q_LSF; i++)
            {
                st->lsf_idx[i] = get_next_indice(st, lsf_q_num_bits[i]);
            }
            st->m_idx = get_next_indice(st, MIRROR_POINT_BITS);
            st->grid_idx = get_next_indice(st, NUM_LSF_GRID_BITS);
        }
        else
        {
            st->lsf_idx[0] = get_next_indice(st, 8);
            st->m_idx = 0;
            st->grid_idx = 0;
        }
        /* shape gains */
        st->idxSubGains = get_next_indice(st, NUM_BITS_SHB_SUBGAINS);

        /* frame gain */
        st->idxFrameGain = get_next_indice(st, NUM_BITS_SHB_FRAMEGAIN);

        if( st->total_brate >= ACELP_24k40 )
        {
            /* sub frame energy*/
            st->idx_shb_fr_gain = get_next_indice(st, NUM_BITS_SHB_ENER_SF);

            /* gain shapes residual */
            for (i = 0; i < NB_SUBFR16k; i++)
            {
                st->idx_res_gs[i] = get_next_indice(st, NUM_BITS_SHB_RES_GS);
            }

            /* voicing factor */
            st->idx_mixFac = get_next_indice(st, NUM_BITS_SHB_VF);
        }

        if( st->tec_tfa == 1 )
        {
            st->tec_flag = get_next_indice(st, BITS_TEC);
            st->tfa_flag = get_next_indice(st, BITS_TFA);
            if( st->tfa_flag && st->tec_flag )
            {
                st->tec_flag = 2;
                st->tfa_flag = 0;
            }


        }
        else
        {
            st->tec_flag = 0;
            st->tfa_flag = 0;
        }
    }

    if( st->bwidth == FB )
    {
        st->idxGain = get_next_indice(st, 4);
    }

    return;
}


/*---------------------------------------------------------------------*
 * GenTransition()
 *
 * Generate a highband transition signal from the gain shape overlap
 * buffer to fill the gap caused by the delay alignment buffer when
 * switching from TBE to IGF
 *---------------------------------------------------------------------*/
void GenTransition(
    const float *input,                         /* i  : gain shape overlap buffer            */
    const float *old_hb_synth,                  /* i  : synthesized HB from previous frame   */
    short length,                         /* i  : targeted length of transition signal */
    float *output,                        /* o  : synthesized transitions signal       */
    float Hilbert_Mem[],                  /* i/o: memory                               */
    float state_lsyn_filt_shb_local[],    /* i/o: memory                               */
    short *syn_dm_phase,
    int   output_Fs,
    float *up_mem,
    int   rf_flag
    , int bitrate
)
{
    short i;
    float syn_overlap_32k[L_FRAME32k];

    /* upsample overlap snippet */
    Interpolate_allpass_steep( input, state_lsyn_filt_shb_local, SHB_OVERLAP_LEN, syn_overlap_32k );

    /* perform spectral flip and downmix with overlap snippet to match HB synth  */
    if (rf_flag || bitrate == ACELP_9k60)
    {
        flip_and_downmix_generic( syn_overlap_32k, syn_overlap_32k, 2*SHB_OVERLAP_LEN, Hilbert_Mem,
                                  Hilbert_Mem + HILBERT_ORDER1, Hilbert_Mem + (HILBERT_ORDER1+2*HILBERT_ORDER2), syn_dm_phase );
    }
    else
    {
        for(i = 0; i < 2*SHB_OVERLAP_LEN; i++)
        {
            syn_overlap_32k[i] = ((i%2)==0)?(-syn_overlap_32k[i]):(syn_overlap_32k[i]);
        }
    }

    /* cross fade of overlap snippet and mirrored HB synth from previous frame */
    for (i=0; i<2*L_SHB_LAHEAD; i++)
    {
        output[i] = window_shb_32k[i]*old_hb_synth[L_SHB_TRANSITION_LENGTH-1-i] + window_shb_32k[2*L_SHB_LAHEAD-1-i]*syn_overlap_32k[i];
    }

    /* fill transition signal with mirrored HB synth from previous frame to fully fill delay alignment buffer gap */
    for ( ; i < length; i++)
    {
        output[i] = old_hb_synth[L_SHB_TRANSITION_LENGTH-1-i];
    }

    if( output_Fs == 48000 )
    {
        interpolate_3_over_2_allpass( output, length, output, up_mem, allpass_poles_3_ov_2 );
    }

    return;
}


void GenTransition_WB(
    const float *input,                         /* i  : gain shape overlap buffer            */
    const float *old_hb_synth,                  /* i  : synthesized HB from previous frame   */
    short length,                         /* i  : targeted length of transition signal */
    float *output,                        /* o  : synthesized transitions signal       */
    float state_lsyn_filt_shb1[],
    float state_lsyn_filt_shb2[],
    int   output_Fs,
    float *up_mem
)
{
    short i;
    float speech_buf_16k1[L_FRAME16k], speech_buf_16k2[L_FRAME16k];
    float upsampled_synth[L_FRAME48k];

    /* upsample overlap snippet */
    Interpolate_allpass_steep( input, state_lsyn_filt_shb1, SHB_OVERLAP_LEN/2, speech_buf_16k1);
    Interpolate_allpass_steep( speech_buf_16k1, state_lsyn_filt_shb2, SHB_OVERLAP_LEN, speech_buf_16k2);


    /* perform spectral flip and downmix with overlap snippet to match HB synth  */
    for(i = 0; i < SHB_OVERLAP_LEN; i++)
    {
        speech_buf_16k2[i] = ((i%2)==0)?(-speech_buf_16k2[i]):(speech_buf_16k2[i]);
    }

    /* cross fade of overlap snippet and mirrored HB synth from previous frame */
    for (i=0; i<L_SHB_LAHEAD; i++)
    {
        output[i] = window_shb[i]*old_hb_synth[L_SHB_TRANSITION_LENGTH-1-i] + window_shb[L_SHB_LAHEAD-1-i]*speech_buf_16k2[i];
        output[i] *= 0.65f;
    }

    /* fill transition signal with mirrored HB synth from previous frame to fully fill delay alignment buffer gap */
    for ( ; i < length; i++)
    {
        output[i] = old_hb_synth[L_SHB_TRANSITION_LENGTH-1-i];
        output[i] *= 0.65f;
    }

    /* upsampling if necessary */
    if( output_Fs == 32000 )
    {
        Interpolate_allpass_steep( output, up_mem, L_FRAME16k, upsampled_synth );
        mvr2r( upsampled_synth, output, L_FRAME32k );
    }
    else if( output_Fs == 48000 )
    {
        interpolate_3_over_1_allpass( output, L_FRAME16k, upsampled_synth, up_mem, allpass_poles_3_ov_2 );
        mvr2r( upsampled_synth, output, L_FRAME48k );
    }

    return;
}


void TBEreset_dec(
    Decoder_State *st,                         /* i/o: decoder state structure                 */
    short  bandwidth                           /* i  : bandwidth mode                          */
)
{
    if( st->last_core != ACELP_CORE )
    {
        set_f( st->old_bwe_exc, 0.0f, PIT16k_MAX * 2 );
        st->bwe_non_lin_prev_scale = 0.f;
    }
    if( bandwidth == WB )
    {
        wb_tbe_extras_reset( st->mem_genSHBexc_filt_down_wb2, st->mem_genSHBexc_filt_down_wb3 );
        wb_tbe_extras_reset_synth( st->state_lsyn_filt_shb, st->state_lsyn_filt_dwn_shb, st->mem_resamp_HB );
        set_f( st->mem_genSHBexc_filt_down_shb, 0, 7 );
        set_f( st->state_lpc_syn, 0, 10 );
        set_f( st->state_syn_shbexc, 0, L_SHB_LAHEAD/4 );
        set_f( st->syn_overlap, 0, L_SHB_LAHEAD );
        set_f( st->mem_csfilt, 0, 2 );
    }
    else if( bandwidth == SWB || bandwidth == FB )
    {
        swb_tbe_reset( st->mem_csfilt, st->mem_genSHBexc_filt_down_shb, st->state_lpc_syn,
                       st->syn_overlap, st->state_syn_shbexc, &(st->tbe_demph), &(st->tbe_premph),
                       st->mem_stp_swb, &(st->gain_prec_swb) );

        set_f( st->GainShape_Delay, 0, NUM_SHB_SUBFR/2 );
        set_f( st->int_3_over_2_tbemem_dec, 0.f, INTERP_3_2_MEM_LEN );
        set_f( st->mem_resamp_HB_32k, 0, 2*ALLPASSSECTIONS_STEEP+1 );

        swb_tbe_reset_synth( st->genSHBsynth_Hilbert_Mem, st->genSHBsynth_state_lsyn_filt_shb_local );

        if ( bandwidth == FB )
        {
            st->prev_fb_ener_adjust = 0.0f;
            set_f(st->fb_state_lpc_syn, 0, LPC_SHB_ORDER);
            st->fb_tbe_demph = 0;
            fb_tbe_reset_synth( st->fbbwe_hpf_mem, &st->prev_fbbwe_ratio );
        }
    }

    return;
}
