/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"


/*-------------------------------------------------------------------*
 * reset_rf_indices()
 *
 * Initialization of oartial redundancy coding
 *-------------------------------------------------------------------*/

void reset_rf_indices(
    Encoder_State *st         /* i: state structure - contains partial RF indices */
)
{
    short i, j;

    st->rf_frame_type = 0; /* since this function is called every frame this will happen even for a SID frame, hence treating it as GSC frame, i.e no RF encoding */

    st->rf_mem_w0 = 0;
    set_f( st->rf_clip_var, 0, 6 );
    st->rf_tilt_code = 0;
    set_f( st->rf_mem_syn2, 0, M );
    set_f( st->rf_dispMem, 0, 8 );
    st->rf_gc_threshold = 0;
    set_f( st->rf_tilt_buf, 0, NB_SUBFR16k );

    st->rf_target_bits = 0;
    st->rf_target_bits_write = 0;
    st->rf_tcxltp_pitch_int_past = st->L_frame;
    st->rf_last_tns_active = 0;
    st->rf_second_last_tns_active = 0;
    st->rf_second_last_core= 0;

    for( i = 0; i < MAX_RF_FEC_OFFSET; i++ )
    {
        st->rf_indx_frametype[i] = RF_NO_DATA;
        st->rf_targetbits_buff[i] = 6;  /* rf_mode: 1, rf_frame_type: 3, and fec_offset: 2 */
        st->rf_indx_lsf[i][0] = 0;
        st->rf_indx_lsf[i][1] = 0;
        st->rf_indx_lsf[i][2] = 0;
        st->rf_indx_EsPred[i] = 0;
        st->rf_indx_nelp_fid[i] = 0;
        st->rf_indx_nelp_iG1[i] = 0;
        st->rf_indx_nelp_iG2[i][0] = 0;
        st->rf_indx_nelp_iG2[i][1] = 0;

        for( j = 0; j < NB_SUBFR16k; j++ )
        {
            st->rf_indx_ltfMode[i][j] = 0;
            st->rf_indx_pitch[i][j] = 0;
            st->rf_indx_fcb[i][j] = 0;
            st->rf_indx_gain[i][j] = 0;
        }

        st->rf_clas[i] = UNVOICED_CLAS;
        st->rf_gain_tcx[i] = 0;
        st->rf_tcxltp_param[i] = 0;

        st->rf_indx_tbeGainFr[i] = 0;
    }

    return;
}


/*-------------------------------------------------------------------*
 * coder_acelp_rf()
 *
 * Encode excitation signal (partial redundancy)
 *-------------------------------------------------------------------*/

void coder_acelp_rf(
    const short target_bits,    /* i:   target bits                 */
    const float speech[],       /* i  : speech[-M..lg]              */
    const short coder_type,     /* i  : coding type                 */
    const short rf_frame_type,  /* i  : rf_frame_type               */
    const float A[],            /* i  : coefficients 4xAz[M+1]      */
    const float Aq[],           /* i  : coefficients 4xAz_q[M+1]    */
    const float voicing[],      /* i  : open-loop LTP gain          */
    const short T_op[],         /* i  : open-loop LTP lag           */
    const float stab_fac,       /* i  : LP stability factor         */
    Encoder_State *st,            /* i/o : coder memory state         */
    ACELP_config *acelp_cfg,     /* i/o: configuration of the ACELP  */
    float *exc_rf,        /* i/o: pointer to RF excitation    */
    float *syn_rf         /* i/o: pointer to RF synthesis     */
)
{
    short i, i_subfr, nSubfr;
    int T0, T0_min, T0_min_frac, T0_max, T0_max_frac, T0_res;
    int T0_frac;
    float Es_pred_rf;
    float gain_pit, gain_code, voice_fac;
    float prev_gain_pit;
    ACELP_CbkCorr g_corr;
    float g_corr2[6];
    const float *p_A, *p_Aq;
    float code[L_SUBFR];
    float xn[L_SUBFR], cn[L_SUBFR], h1[L_SUBFR];
    float xn2[L_SUBFR], y1[L_SUBFR], y2[L_SUBFR];
    float res_save;
    float exc2[L_SUBFR];
    float exc_nelp[L_FRAME];
    float syn2[L_FRAME16k];
    float past_gcode, gain_inov;
    short clip_gain;
    float gain_code2;
    float code2[L_SUBFR];
    float y22[L_SUBFR];
    float lp_select;
    int *prm_rf;

    /*-----------------------------------------------------------------------*
     * Initialization
     *------------------------------------------------------------------------*/
    past_gcode = 0;
    gain_inov = 0;
    T0 = 0;
    T0_res = 0;
    T0_frac = 0;
    gain_pit = 0;
    gain_code = 0;
    voice_fac = 0;
    prev_gain_pit=0;
    Es_pred_rf = 0;
    set_f(code, 0.0f, L_SUBFR);

    /*-----------------------------------------------------------------------*
     * Configure ACELP partial copy                                          *
     *-----------------------------------------------------------------------*/

    BITS_ALLOC_config_acelp( target_bits, rf_frame_type, &(st->acelp_cfg_rf), 0, st->nb_subfr );

    /* Reset phase dispersion */
    if( st->last_core > 0 )
    {
        set_zero( st->rf_dispMem, 8 );
    }

    /*---------------------------------------------------------------*
     * Calculation of LP residual (filtering through A[z] filter)
     *---------------------------------------------------------------*/

    calc_residu( speech, exc_rf, Aq, st->L_frame );

    /*------------------------------------------------------------------------*
     * Find and quantize mean_ener_code for gain quantizer                    *
     *------------------------------------------------------------------------*/

    Es_pred_rf = 0;
    if( acelp_cfg->nrg_mode > 0 && rf_frame_type != RF_NELP )
    {
        Es_pred_enc( &Es_pred_rf, &st->rf_indx_EsPred[0], st->L_frame, L_SUBFR, exc_rf, voicing,
                     acelp_cfg->nrg_bits, acelp_cfg->nrg_mode>1 );
    }

    /*------------------------------------------------------------------*
     * ACELP subframe loop
     *------------------------------------------------------------------*/

    p_A = A;
    p_Aq = Aq;

    res_save = exc_rf[0];
    nSubfr = 0;

    for( i_subfr = 0; i_subfr < st->L_frame; i_subfr += L_SUBFR )
    {
        if( rf_frame_type != RF_NELP )
        {
            /* Restore exc[i_subfr] and save next exc[L_SUBFR+i_subfr] */
            exc_rf[i_subfr] = res_save;
            res_save = exc_rf[L_SUBFR + i_subfr];

            /*--------------------------------------------------------------------------*
             * Find target for pitch search (xn[]), target for innovation search (cn[]) *
             * and impulse response of the weighted synthesis filter (h1[]).            *
             *--------------------------------------------------------------------------*/

            find_targets( speech, &syn_rf[i_subfr-M], i_subfr, &(st->rf_mem_w0),
                          p_Aq, exc_rf, L_SUBFR, p_A, st->preemph_fac, xn, cn, h1 );
        }

        /* full frame nelp partial copy encoding */
        if( rf_frame_type == RF_NELP )
        {
            if( i_subfr == 0 )
            {
                nelp_encoder( st, exc_rf, exc_nelp
                              ,0
                            );
            }
            mvr2r( &exc_nelp[i_subfr], exc2, L_SUBFR );
            mvr2r( &exc_nelp[i_subfr], exc_rf, L_SUBFR );
        }
        else
        {
            /*-----------------------------------------------------------------*
            * Gain clipping test to avoid unstable synthesis on frame erasure
            * or in case of floating point encoder & fixed p. decoder
            *-----------------------------------------------------------------*/

            clip_gain = gp_clip( st->core_brate, voicing, i_subfr, coder_type, xn, st->rf_clip_var );

            /*-----------------------------------------------------------------*
            * - find unity gain pitch excitation (adaptive codebook entry)    *
            *   with fractional interpolation.                                *
            * - find filtered pitch exc. y1[]=exc[] convolved with h1[])      *
            * - compute pitch gain1                                           *
            *-----------------------------------------------------------------*/

            if( acelp_cfg->gains_mode[i_subfr/L_SUBFR] == 0 )
            {
                gain_pit = prev_gain_pit;
            }

            if ( acelp_cfg->ltp_bits != 0 )
            {
                prm_rf = &st->rf_indx_pitch[0][nSubfr];

                /* Adaptive Codebook (GC and VC) */
                Mode2_pit_encode( acelp_cfg->ltp_mode, i_subfr, &prm_rf, &exc_rf[i_subfr],
                                  T_op, &T0_min, &T0_min_frac, &T0_max, &T0_max_frac, &T0, &T0_frac, &T0_res, h1, xn,
                                  st->pit_min, st->pit_fr1, st->pit_fr1b, st->pit_fr2, st->pit_max, st->pit_res_max );

                /* find ACB excitation */
                if( T0_res == (st->pit_res_max>>1) )  /* st->pit_res_max is 4 for 12.8kHz core */
                {
                    pred_lt4( &exc_rf[i_subfr], &exc_rf[i_subfr], T0, T0_frac<<1, L_SUBFR+1, inter4_2, L_INTERPOL2, PIT_UP_SAMP);
                }
                else
                {
                    pred_lt4( &exc_rf[i_subfr], &exc_rf[i_subfr], T0, T0_frac, L_SUBFR+1, inter4_2, L_INTERPOL2, PIT_UP_SAMP);
                }


                /* filter adaptive codebook */
                lp_select =  lp_filt_exc_enc( MODE2, st->core_brate, 0, (acelp_cfg->gains_mode[i_subfr/L_SUBFR]>0)?(acelp_cfg->gains_mode[i_subfr/L_SUBFR]):(100),
                                              i_subfr, exc_rf, h1, xn, y1, xn2, L_SUBFR, st->L_frame, g_corr2, clip_gain, &(gain_pit), &(acelp_cfg->ltf_mode) );

                if( acelp_cfg->ltf_mode == NORMAL_OPERATION )
                {
                    st->rf_indx_ltfMode[0][nSubfr] = lp_select;
                }

                g_corr.y1y1 = g_corr2[0];
                g_corr.xy1 = -0.5f*(g_corr2[1]-0.01f)+0.01f;
            }
            else
            {
                gain_pit=0.f;
                g_corr.xy1=0.f;
                g_corr.y1y1=0.f;
                set_zero( y1, L_SUBFR );
                set_zero( exc_rf+i_subfr, L_SUBFR );
                T0 = L_SUBFR;
                T0_frac = 0;
                T0_res = 1;
            }

            /*----------------------------------------------------------------------*
             *                 Encode the algebraic innovation                      *
             *----------------------------------------------------------------------*/

            if( acelp_cfg->fixed_cdk_index[i_subfr/L_SUBFR] >= 0 )
            {
                prm_rf = &st->rf_indx_fcb[0][nSubfr];
                E_ACELP_innovative_codebook( exc_rf, T0, T0_frac, T0_res, gain_pit, st->rf_tilt_code,
                                             acelp_cfg->fixed_cdk_index[i_subfr/L_SUBFR],
                                             acelp_cfg->pre_emphasis,
                                             acelp_cfg->pitch_sharpening,
                                             acelp_cfg->phase_scrambling,
                                             acelp_cfg->formant_enh,
                                             acelp_cfg->formant_tilt,
                                             acelp_cfg->formant_enh_num,
                                             acelp_cfg->formant_enh_den,
                                             i_subfr, p_Aq, h1, xn, cn, y1, y2,
                                             st->acelp_autocorr, &prm_rf, code
                                             ,st->L_frame, st->last_L_frame, st->total_brate

                                           );
            }
            else
            {
                set_f( code, 0.0f, L_SUBFR );
                set_f( y2, 0.0f, L_SUBFR );
            }

            if( i_subfr < (st->L_frame - L_SUBFR) )
            {
                E_corr_xy2( xn, y1, y2, g_corr2, L_SUBFR );
                g_corr.y2y2 = 0.01F + g_corr2[2];
                g_corr.xy2 = 0.01F + -0.5f*g_corr2[3];
                g_corr.y1y2 = 0.01F + 0.5f*g_corr2[4];

                g_corr.xx = 0.01F + dotp( xn, xn, L_SUBFR );

                /*----------------------------------------------------------------------*
                 *                 Add Gaussian excitation                              *
                 *----------------------------------------------------------------------*/

                gain_code2 = 0.f;
                set_zero( code2, L_SUBFR );
                set_zero( y22, L_SUBFR );

                /*----------------------------------------------------------*
                 *  - Compute the fixed codebook gain                       *
                 *  - quantize fixed codebook gain                          *
                 *----------------------------------------------------------*/

                if( acelp_cfg->gains_mode[i_subfr/L_SUBFR] != 0 )
                {
                    prm_rf = &st->rf_indx_gain[0][nSubfr];
                    encode_acelp_gains( code, acelp_cfg->gains_mode[i_subfr/L_SUBFR], Es_pred_rf, clip_gain, &g_corr, &gain_pit, &gain_code,
                                        &prm_rf, &past_gcode, &gain_inov, L_SUBFR, code2, &gain_code2, st->flag_noisy_speech_snr );
                }

                gp_clip_test_gain_pit( st->core_brate, gain_pit, st->rf_clip_var );

                /*----------------------------------------------------------*
                 * - voice factor (for codebook tilt sharpening)            *
                 *----------------------------------------------------------*/

                st->rf_tilt_code = est_tilt( exc_rf+i_subfr, gain_pit, code, gain_code, &voice_fac, L_SUBFR, acelp_cfg->voice_tilt );

                /*-----------------------------------------------------------------*
                * Update memory of the weighting filter
                *-----------------------------------------------------------------*/

                st->rf_mem_w0 = xn[L_SUBFR-1] - gain_pit * y1[L_SUBFR-1] - gain_code * y2[L_SUBFR-1] - gain_code2*y22[L_SUBFR-1];

                /*-------------------------------------------------------*
                * - Find the total excitation.                          *
                *-------------------------------------------------------*/

                for( i = 0; i < L_SUBFR; i++ )
                {
                    exc2[i] = gain_pit * exc_rf[i+i_subfr];
                    exc2[i] +=  gain_code2 * code2[i];
                    exc_rf[i+i_subfr] = exc2[i] + gain_code * code[i];
                }

                /*---------------------------------------------------------*
                * Enhance the excitation                                  *
                *---------------------------------------------------------*/

                enhancer( MODE2, -1, acelp_cfg->fixed_cdk_index[i_subfr/L_SUBFR], 0, coder_type, st->L_frame, voice_fac, stab_fac,
                          past_gcode, gain_inov, &st->rf_gc_threshold, code, exc2, gain_pit, st->rf_dispMem );
            }
        }

        if( (i_subfr < (st->L_frame - L_SUBFR)) || (rf_frame_type != RF_NELP) )
        {
            /*----------------------------------------------------------*
             * - compute the synthesis speech                           *
             *----------------------------------------------------------*/

            syn_filt( p_Aq, M, exc2, &syn2[i_subfr], L_SUBFR, st->rf_mem_syn2, 1 );

            syn_filt( p_Aq, M, &exc_rf[i_subfr], &syn_rf[i_subfr], L_SUBFR, &syn_rf[i_subfr-M], 0 );

            /*----------------------------------------------------------*
             * Updates                                                  *
             *----------------------------------------------------------*/

            p_A += (M+1);
            p_Aq += (M+1);
            nSubfr++;

            /* copy current gain for next subframe use, in case there is no explicit encoding */
            prev_gain_pit = gain_pit;
        }

    } /* end of subframe loop */

    return;
}
