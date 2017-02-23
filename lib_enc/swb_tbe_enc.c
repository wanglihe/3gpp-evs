/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <stdlib.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"
#include "rom_enc.h"

/*-----------------------------------------------------------------*
 * Local functions
 *-----------------------------------------------------------------*/

static void return_M_Least( const float *inp, const short n_cols, const float *codebook, const short num_grp,
                            const float *weight, const short interNum, short *least );
static void singlevectortest_gain( const float *inp, const short dimen, const short cb_size, short *index,
                                   const float *weight, float *recon, const float *codebook );
static void determine_gain_weights( const float *gain, float *weights, const short dims );
static void QuantizeSHBsubgains( Encoder_State *st, float *subgains, const short extl );
static void QuantizeSHBframegain( Encoder_State *st, float *GainFrame, const short extl, long extl_brate, short *rf_gainFrame_ind );
static short closest_centroid( const float *data, const float *weights, const float *quantizer, const short centroids, const short length);
static void EstimateSHBFrameGain( const short length, const float *oriSHB,  const float *synSHB, float *subgain, float *GainFrame, const float *win_shb, const float *subwin_shb );
static void EstimateSHBGainShape( const short length, const float *oriSHB, const float *synSHB,  float *subgain, const float *subwin );
static float pow_off_pk(float a[], short len, short step);
static void Quant_lower_LSF( const float lsf[], float lsf_q[], short lsf_idx[] );
static short Quant_mirror_point(const float lsf[], const float lsf_q[], float *m);
static short Find_LSF_grid(const float lsf[], float lsf_q[], const float m);
static void gainFrSmooth_En(Encoder_State *st, float *shb_frame, const float *lpc_shb, const float *lsp_shb, float *MA_lsp_shb_spacing, short *frGainAttenuate, short *frGainSmoothEn );
static void Quant_BWE_LSF( Encoder_State *st, const float lsp_shb[], float Q_lsfs[] );
static void Quant_shb_ener_sf( Encoder_State *st, float *shb_ener_sf );
static void Quant_shb_res_gshape( Encoder_State *st, float *shb_res_gshape );


/*-------------------------------------------------------------------*
 * InitSWBencBuffer()
 *
 * Initialize SWB buffers
 *-------------------------------------------------------------------*/

void InitSWBencBuffer(
    Encoder_State *st  /* i/o: SHB encoder structure */
)
{
    set_f( st->old_speech_shb, 0.0f, L_LOOK_16k + L_SUBFR16k );
    set_f(st->old_bwe_exc, 0.0f, (PIT16k_MAX * 2));
    st->bwe_seed[0] = 23;
    st->bwe_seed[1] = 59;
    set_f( st->old_bwe_exc_extended, 0.0f, NL_BUFF_OFFSET );
    st->bwe_non_lin_prev_scale = 0;


    set_f( st->state_ana_filt_shb, 0.0f, (2*ALLPASSSECTIONS_STEEP+1) );

    set_f( st->elliptic_bpf_2_48k_mem[0], 0.0f, 4 );
    set_f( st->elliptic_bpf_2_48k_mem[1], 0.0f, 4 );
    set_f( st->elliptic_bpf_2_48k_mem[2], 0.0f, 4 );
    set_f( st->elliptic_bpf_2_48k_mem[3], 0.0f, 4 );
    st->prev_fb_energy = 0.0f;

    return;
}


/*-------------------------------------------------------------------*
 * ResetSHBbuffer_Enc()
 *
 *-------------------------------------------------------------------*/

void ResetSHBbuffer_Enc(
    Encoder_State *st      /* i/o: SHB encoder structure */
)
{
    /* states for the filters used in generating SHB excitation from WB excitation*/
    set_f( st->mem_genSHBexc_filt_down_shb,0, (2*ALLPASSSECTIONS_STEEP+1) );
    set_f( st->mem_csfilt, 0, 2 );

    /* states for the filters used in generating SHB signal from SHB excitation*/
    set_f( st->state_syn_shbexc, 0, L_SHB_LAHEAD );
    set_f( st->state_lpc_syn, 0, LPC_SHB_ORDER );
    if( sub(st->extl, FB_TBE) == 0 )
    {
        set_f( st->fb_state_lpc_syn, 0, LPC_SHB_ORDER );
        st->fb_tbe_demph = 0;
    }
    /* states for the filters used in generating WB signal from WB excitation*/
    set_f( st->decim_state1, 0.0f, (2*ALLPASSSECTIONS_STEEP+1) );
    set_f( st->decim_state2, 0.0f, (2*ALLPASSSECTIONS_STEEP+1) );
    set_f( st->mem_genSHBexc_filt_down_wb2, 0, (2*ALLPASSSECTIONS_STEEP+1) );
    set_f( st->mem_genSHBexc_filt_down_wb3, 0, (2*ALLPASSSECTIONS_STEEP+1) );


    /* overlap buffer used to Adjust SHB Frame Gain */
    set_f( st->mem_stp_swb, 0, LPC_SHB_ORDER );
    st->gain_prec_swb = 1.0f;
    set_f( st->syn_overlap,0,L_SHB_LAHEAD );
    st->tbe_demph = 0.0f;
    st->tbe_premph = 0.0f;

    return;
}


/*-------------------------------------------------------------------*
* wb_tbe_enc()
*
* WB TBE encoder, 6 - 8 kHz band encoding module
*-------------------------------------------------------------------*/

void wb_tbe_enc(
    Encoder_State *st,            /* i/o: encoder state structure             */
    const short coder_type,           /* i  : coding type                         */
    const float *hb_speech,           /* i  : HB target signal (6-8kHz) at 16kHz  */
    const float *bwe_exc_extended,    /* i  : bandwidth extended exciatation      */
    const float voice_factors[],      /* i  : voicing factors                     */
    const float pitch_buf[],          /* i  : pitch for each subframe             */
    const float voicing[]             /* i  : OL maximum normalized correlation   */
)
{
    short i, j, k, delay;
    float hb_old_speech[(L_LOOK_12k8 + L_SUBFR + L_FRAME) * 5/16];
    float shaped_wb_excitation [ (L_FRAME16k + L_SHB_LAHEAD)/4];
    float exc4kWhtnd [L_FRAME16k / 4];
    short ana_align_delay = - L_SHB_LAHEAD / 4 - 5;
    float GainFrame, GainShape[NUM_SHB_SUBFR];
    float lpc_wb[LPC_SHB_ORDER_WB+1];
    float lsp_wb[LPC_SHB_ORDER_WB], weights_lsp[LPC_SHB_ORDER_WB] = {1.0, 1.0};
    float *hb_new_speech, *hb_frame;
    float R[LPC_SHB_ORDER_WB+2], ervec[LPC_SHB_ORDER_WB+1];
    float prev_pow, curr_pow, scale;
    float p2m_in, p2m_out;
    short uv_flag;
    float pitBufAvg, voicingBufAvg;
    float vf_modified[NB_SUBFR16k];
    float temp_wb_fac, feedback;
    float lsp_spacing;
    float ervec_temp[LPC_SHB_ORDER_WB+1];
    float lsp_wb_temp[LPC_SHB_ORDER_WB], lpc_wb_temp[LPC_SHB_ORDER_WB+1];

    /* delay alignment */
    delay = (L_LOOK_12k8 + L_SUBFR) * 5/16;

    hb_new_speech = hb_old_speech + delay;
    hb_frame = hb_old_speech + L_SUBFR * 5/16 + ana_align_delay;

    mvr2r( st->old_speech_wb, hb_old_speech, delay );
    mvr2r( hb_speech, hb_new_speech, L_FRAME16k / 4 );
    mvr2r( hb_old_speech + L_FRAME16k / 4, st->old_speech_wb, delay );
    if( ( st->last_extl != WB_TBE && st->last_extl != WB_BWE ) &&
            ( st->clas == UNVOICED_CLAS || ( voicing[0] < 0.5f && voicing[1] < 0.5f && voicing[2] < 0.5f ) ) &&
            ( !st->igf ) )
    {
        /* In case of unvoiced signals after switching cores, back-propagate the target signal */
        mvr2r( hb_speech, hb_old_speech, delay );

        for( i = (L_LOOK_12k8 + L_SUBFR) * 5/16, j = k = 0; j < L_SUBFR16k; i--, j+=4, k++ )
        {
            hb_old_speech[i] *= ola_win_shb_switch_fold[j];
            hb_old_speech[i] += hb_speech[k] * ola_win_shb_switch_fold[L_SUBFR16k-4-j];
        }
    }

    autocorr( hb_old_speech, R, LPC_SHB_ORDER_WB+1, (NS2SA(INT_FS_12k8, 5000000L) + L_SUBFR + L_FRAME) * 5/16, win_lpc_hb_wb, 0, 1, 1 ); /* VE2QC: to be tuned later */
    lev_dur( lpc_wb_temp, R, LPC_SHB_ORDER_WB, ervec_temp );
    a2lsp( lsp_wb_temp, lpc_wb_temp, LPC_SHB_ORDER_WB );
    lsp_spacing = 0.5f;

    for( i = 0; i < LPC_SHB_ORDER_WB; i++ )
    {
        lsp_spacing = min(lsp_spacing, (float)( i == 0 ? lsp_wb_temp[0] : (lsp_wb_temp[i] - lsp_wb_temp[i -1])));
    }

    /* Spectral smoothing of autocorrelation coefficients */
    for (i = 0; i <= LPC_SHB_ORDER_WB; i++)
    {
        R[i] = R[i] * wac[i];
    }
    R[0] = max( R[0], 1.0e-8f );

    if ( st->rf_mode == 1 || st->extl_brate == WB_TBE_0k35 )
    {
        lev_dur( lpc_wb, R, LPC_SHB_ORDER_LBR_WB, ervec );

        /* Expand bandwidth of the LP coeffs */
        for (i = 0; i <= LPC_SHB_ORDER_LBR_WB; i++)
        {
            lpc_wb[i] *= lpc_weights[i];
        }

        /* convert into lsps and calculate weights */
        a2lsp(lsp_wb,lpc_wb, LPC_SHB_ORDER_LBR_WB);
        lsp_weights(lsp_wb, weights_lsp, LPC_SHB_ORDER_LBR_WB);

        /* Quantization of LSFs */
        i = closest_centroid( lsp_wb, weights_lsp, lbr_wb_bwe_lsfvq_cbook_2bit, 4, LPC_SHB_ORDER_LBR_WB );
        if( st->codec_mode == MODE2 )
        {
            st->lsf_WB = i;
        }
        else
        {
            push_indice( st, IND_SHB_LSF, i, NUM_BITS_LBR_WB_LSF );
        }

        mvr2r( lbr_wb_bwe_lsfvq_cbook_2bit + i*LPC_SHB_ORDER_LBR_WB, lsp_wb, LPC_SHB_ORDER_LBR_WB);

        lsp2a( lpc_wb, lsp_wb, LPC_SHB_ORDER_LBR_WB );
        set_f( lpc_wb + LPC_SHB_ORDER_LBR_WB+1, 0.0f, (LPC_SHB_ORDER_WB - LPC_SHB_ORDER_LBR_WB) );
    }
    else
    {
        lev_dur( lpc_wb, R, LPC_SHB_ORDER_WB, ervec );

        /* Expand bandwidth of the LP coeffs */
        for (i = 0; i <= LPC_SHB_ORDER_WB; i++)
        {
            lpc_wb[i] *= lpc_weights[i];
        }

        /* convert into lsps and calculate weights */
        a2lsp(lsp_wb,lpc_wb, LPC_SHB_ORDER_WB);
        lsp_weights(lsp_wb, weights_lsp, LPC_SHB_ORDER_WB);

        /* Quantization of LSFs */
        i = closest_centroid( lsp_wb, weights_lsp, wb_bwe_lsfvq_cbook_8bit, 256, LPC_SHB_ORDER_WB );
        push_indice( st, IND_SHB_LSF, i, NUM_BITS_WB_LSF );

        mvr2r( wb_bwe_lsfvq_cbook_8bit + i*LPC_SHB_ORDER_WB, lsp_wb, LPC_SHB_ORDER_WB);

        lsp2a( lpc_wb, lsp_wb, LPC_SHB_ORDER_WB );
    }

    uv_flag = 0;
    if ( st->extl_brate == WB_TBE_1k05 && st->coder_type_raw == UNVOICED )
    {
        uv_flag = 1;
    }

    mvr2r( voice_factors, vf_modified, NB_SUBFR16k );
    if( coder_type == VOICED )
    {
        for ( i = 1; i < NB_SUBFR; i++ )
        {
            vf_modified[i] = 0.8f * voice_factors[i] + 0.2f * voice_factors[i-1];
        }

        if(st->L_frame != L_FRAME )
        {
            vf_modified[4] = 0.8f * voice_factors[4] + 0.2f * voice_factors[3];
        }
    }

    /* From low band excitation, generate highband excitation */
    mvr2r( st->state_syn_shbexc, shaped_wb_excitation, L_SHB_LAHEAD/4 );
    GenShapedWBExcitation( shaped_wb_excitation + L_SHB_LAHEAD/4, lpc_wb, exc4kWhtnd, st->mem_csfilt, st->mem_genSHBexc_filt_down_shb,
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
        scale = sqrt(curr_pow/prev_pow);
    }

    for( i = 0; i < (L_SHB_LAHEAD/4 - 1); i++ )
    {
        shaped_wb_excitation[i] *= scale;
    }

    scale = sqrt( scale );

    shaped_wb_excitation[L_SHB_LAHEAD/4 - 1] *= scale;

    /* Update WB excitation */
    mvr2r( shaped_wb_excitation + L_FRAME16k/4, st->state_syn_shbexc, L_SHB_LAHEAD/4 );

    /* estimate the gain shape parameter */
    EstimateSHBGainShape( SHB_OVERLAP_LEN/2, hb_frame, shaped_wb_excitation, GainShape, subwin_wb );

    /* Gain frame adjustment factor */
    temp_wb_fac = (float)log( (GainShape[0]+0.00001f) / (st->prev_wb_GainShape+0.0001f) );
    feedback = temp_wb_fac * temp_wb_fac ;
    for (i = 1; i < NUM_SHB_SUBFR/4; i++)
    {
        temp_wb_fac = (float)log( (GainShape[i]+0.00001f) / (GainShape[i-1]+0.0001f) );
        feedback += (temp_wb_fac * temp_wb_fac);
    }
    feedback = 0.4f / (1 + 0.5f * feedback);

    temp_wb_fac = st->prev_wb_GainShape;
    for (i = 0; i < NUM_SHB_SUBFR/4; i++)
    {
        GainShape[i] = (1 - feedback) * GainShape[i] + feedback * temp_wb_fac;
        temp_wb_fac = GainShape[i];
    }
    st->prev_wb_GainShape = GainShape[NUM_SHB_SUBFR/4-1];

    /* Compute the power of gains away from the peak gain prior to quantization */
    p2m_in = pow_off_pk(GainShape, NUM_SHB_SUBFR / 4, 1);

    if ( st->extl_brate == WB_TBE_0k35 )
    {
        for (i=0; i<8; i++)
        {
            GainShape[i] = RECIP_ROOT_EIGHT;
        }
    }
    else
    {
        push_indice( st, IND_UV_FLAG, uv_flag, 1 );

        /* Quantization of the subframe gain parameter */
        QuantizeSHBsubgains( st, GainShape, st->extl );
    }

    /* Compute the power of gains away from the peak gain after quantization */
    p2m_out = pow_off_pk(GainShape, NUM_SHB_SUBFR/2, 2);

    /* Estimate the gain parameter */
    EstimateSHBFrameGain( SHB_OVERLAP_LEN/2, hb_frame, shaped_wb_excitation, GainShape, &GainFrame, window_wb, subwin_wb );

    /* If there's a big difference in the power of gains away from the peak gain   */
    /* due to poor quantization then suppress energy of the high band.             */
    if ( p2m_out > 2.0f * p2m_in )
    {
        float temp = 0;
        if(p2m_in >= 0 && p2m_out > 0)
        {
            temp = sqrt((2.0f * p2m_in) / p2m_out);
        }
        GainFrame *= temp;
    }

    pitBufAvg = 0.0025f * sum_f( pitch_buf, NB_SUBFR );
    voicingBufAvg = 0.333f * sum_f( voicing, 3 );
    if(voicingBufAvg <= 0.0f && (pitBufAvg != 0))
    {
        voicingBufAvg = pitBufAvg/1.001f;
    }
    else if(voicingBufAvg <= 0.0f)
    {
        voicingBufAvg = 1.0f;
    }

    GainFrame *= max(min((float)(pitBufAvg/voicingBufAvg), 1.0f), 0.7f);

    if( lsp_spacing < 0.01f )
    {
        GainFrame *= 0.65f;
    }

    /* Quantization of the frame gain parameter */
    if( st->igf && coder_type == VOICED )
    {
        GainFrame *= 0.5f;
    }
    else if( st->igf && (0.25f*sum_f(voice_factors, NB_SUBFR) > 0.35f) )
    {
        GainFrame *= 0.75f;
    }

    QuantizeSHBframegain( st, &GainFrame, st->extl, st->extl_brate, &st->RF_bwe_gainFr_ind );

    /* Adjust the subframe and frame gain of the synthesized SHB signal */
    /* Scale the shaped excitation*/
    ScaleShapedSHB( SHB_OVERLAP_LEN/2, shaped_wb_excitation, st->syn_overlap, GainShape, GainFrame, window_wb, subwin_wb );


    return;
}


/*-------------------------------------------------------------------*
 * swb_tbe_enc()
 *
 * SWB TBE encoder, 6 - 14 kHz (or 7.5 - 15.5 kHz) band encoding module
 *-------------------------------------------------------------------*/

void swb_tbe_enc(
    Encoder_State *st,               /* i/o: encoder state structure                */
    const short coder_type,        /* i  : coding type                            */
    float *shb_speech,       /* i  : SHB target signal (6-14kHz) at 16kHz   */
    const float *bwe_exc_extended, /* i  : bandwidth extended exciatation         */
    const float voice_factors[],   /* i  : voicing factors                        */
    float *White_exc16k,     /* o  : shaped white excitation for the FB TBE */
    const float voicing[],         /* i  : OL maximum normalized correlation      */
    const float pitch_buf[]        /* i  : pitch for each subframe                */
)
{
    short i, j, delay;
    float shb_old_speech[L_LOOK_16k + L_SUBFR16k + L_FRAME16k];
    float shaped_shb_excitation [L_FRAME16k + L_SHB_LAHEAD];
    float GainFrame, GainShape[NUM_SHB_SUBFR];
    float lpc_shb[LPC_SHB_ORDER+1];
    float weights_lsp[LPC_SHB_ORDER];
    float *shb_frame, *shb_new_speech;
    float lsf_shb_orig[LPC_SHB_ORDER];
    float lsp_shb_1[LPC_SHB_ORDER], lsp_shb_2[LPC_SHB_ORDER], lsp_temp[LPC_SHB_ORDER];
    float lpc_shb_sf[4*(LPC_SHB_ORDER+1)];
    const float *ptr_lsp_interp_coef;
    short tmp;
    float shb_ener_sf;
    float lsf_shb[LPC_SHB_ORDER];
    float shb_res[L_FRAME16k];
    float shb_res_gshape[NB_SUBFR16k], normFac;
    short vf_ind;
    float sd_uq_q, vf_modified[NB_SUBFR16k];
    float pitBufAvg, voicingBufAvg;
    float R[LPC_SHB_ORDER+2], ervec[LPC_SHB_ORDER+1];
    short ana_align_delay[2] = {-L_SHB_LAHEAD  - (NL_BUFF_OFFSET/2), -L_SHB_LAHEAD  - (NL_BUFF_OFFSET/2)};

    float prev_pow, curr_pow, scale;
    float p2m_in, p2m_out;
    short frGainAttenuate, frGainSmoothEn;
    float MA_lsp_shb_spacing;
    float temp_swb_fac, feedback;
    float shaped_shb_excitationTemp[L_FRAME16k];
    float lsf_diff[LPC_SHB_ORDER], w[LPC_SHB_ORDER];
    float refl[M];
    float tilt_para;
    float formant_fac;
    float temp;
    short stab_check = 1;

    /* initializations */
    set_f( shaped_shb_excitationTemp, 0.0f, L_FRAME16k );

    /* compensate for the delay in target generation and subframe LA */
    shb_frame = shb_old_speech + L_SUBFR16k + ana_align_delay[0];

    /* set up the speech buffers for TBE processing*/
    delay = L_LOOK_16k + L_SUBFR16k;
    shb_new_speech = shb_old_speech +  delay;
    mvr2r( st->old_speech_shb, shb_old_speech, delay );
    mvr2r( shb_speech, shb_new_speech, L_FRAME16k );
    mvr2r( shb_old_speech + L_FRAME16k, st->old_speech_shb, delay );

    autocorr( shb_old_speech,
              R,
              LPC_SHB_ORDER+1,
              NS2SA(INT_FS_16k, ACELP_LOOK_NS) + L_SUBFR16k + L_FRAME16k,
              win_lpc_shb, 0, 1, 1 );

    /* Spectral smoothing of autocorrelation coefficients */
    if(st->rf_mode || st->total_brate == ACELP_9k60)
    {
        for (i = 0; i <= LPC_SHB_ORDER; i++)
        {
            R[i] = R[i] * wac_swb[i];
        }
    }
    /* Set the autocorr[0] element to a non-negative value */
    R[0] = max( R[0], 1.0e-8f);

    lev_dur( lpc_shb, R, LPC_SHB_ORDER, ervec );
    {
        float enerG, lpc_shb1[M+1];

        /* extend the lpc_shb to a 16th order gain calc */
        set_f(lpc_shb1, 0, M+1);
        mvr2r(lpc_shb, lpc_shb1, LPC_SHB_ORDER + 1);

        /* estimate the LP gain */
        enerG = enr_1_Az(lpc_shb1, 2*L_SUBFR);

        /* if the LP gain is greater than a threshold, avoid saturation.
           The function 'is_numeric_float' used to check for infinity enerG */
        if(enerG > 32 || !(is_numeric_float(enerG)) )
        {
            set_f(lpc_shb, 0, LPC_SHB_ORDER+1);
            lev_dur( lpc_shb, R, 2, ervec );
        }
    }

    /* Expand bandwidth of the LP coeffs */
    if(st->rf_mode || st->total_brate == ACELP_9k60)
    {
        for( i = 0; i <= LPC_SHB_ORDER; i++ )
        {
            lpc_shb[i] *= lpc_weights[i];
        }
    }

    /* convert to LSFs */
    stab_check = a2lsp( lsf_shb, lpc_shb, LPC_SHB_ORDER );

    if( (st->last_extl != SWB_TBE && st->last_extl != FB_TBE) || st->ini_frame == 0 )
    {
        for( i=0; i < LPC_SHB_ORDER; i++ )
        {
            st->prev_lsp_shb[i] = i/20.0f;
        }
    }

    if( stab_check == 0 )
    {
        mvr2r( st->prev_lsp_shb, lsf_shb, LPC_SHB_ORDER );
    }

    mvr2r( lsf_shb, st->prev_lsp_shb, LPC_SHB_ORDER );

    mvr2r(lsf_shb, lsf_shb_orig, LPC_SHB_ORDER);

    gainFrSmooth_En( st, shb_frame, lpc_shb, lsf_shb, &MA_lsp_shb_spacing, &frGainAttenuate, &frGainSmoothEn);

    if(st->rf_mode || st->total_brate == ACELP_9k60)
    {
        lsp_weights( lsf_shb, weights_lsp, LPC_SHB_ORDER );

        /* to compensate for the 1.1* weighting done inside the function lsp_weights */
        weights_lsp[3]*=0.909091f;
        weights_lsp[4]*=0.909091f;

        /* 8-bit VQ, 10 dimension */
        i = closest_centroid(lsf_shb, weights_lsp, swb_tbe_lsfvq_cbook_8b, 256, LPC_SHB_ORDER);
        mvr2r(swb_tbe_lsfvq_cbook_8b + i * LPC_SHB_ORDER, lsf_shb, LPC_SHB_ORDER);

        set_s(st->lsf_idx, 0, NUM_Q_LSF);
        st->lsf_idx[0] = i;
    }
    else
    {
        /* Quantization of LSFs */
        Quant_BWE_LSF( st, lsf_shb, lsf_shb );
    }

    space_lsfs( lsf_shb, LPC_SHB_ORDER );

    /* voice factor adjustment and gainframe attenuation factor */
    sd_uq_q = 0;
    for( i = 0; i < LPC_SHB_ORDER; i++ )
    {
        /* Estimate the QD in lsps between UQ and Q*/
        sd_uq_q += (lsf_shb[i] - lsf_shb_orig[i]) * (lsf_shb[i] - lsf_shb_orig[i]);
    }

    mvr2r(voice_factors, vf_modified, NB_SUBFR16k);

    if( coder_type == VOICED || mean(voice_factors, 4) > 0.4f )
    {
        for( i = 1; i < NB_SUBFR; i++ )
        {
            vf_modified[i] = 0.8f * voice_factors[i] + 0.2f * voice_factors[i-1];
        }

        if( st->L_frame != L_FRAME )
        {
            vf_modified[4] = 0.8f * voice_factors[4] + 0.2f * voice_factors[3];
        }
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

    lsf_diff[0] = lsf_diff[LPC_SHB_ORDER-1] = 0.5f;
    for( i=1; i<(LPC_SHB_ORDER-1); i++ )
    {
        lsf_diff[i] = lsf_shb[i] - lsf_shb[i-1];
    }
    a2rc (st->cur_sub_Aq+1, refl, (short) M);

    /* LSP interpolation for 13.2 kbps and 16.4 kbps */
    tilt_para = 6.6956f * (1.0f + refl[0]) * (1.0f + refl[0]) - 3.8714f * (1.0f + refl[0]) + 1.3041f;
    if( st->last_extl != SWB_TBE )
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
                w[i] = (lsf_diff[i] < st->prev_lsf_diff[i-1]) ? min(max(0.8f*lsf_diff[i]/st->prev_lsf_diff[i-1], 0.5f), 1.0f) : min(max(0.8f*st->prev_lsf_diff[i-1]/lsf_diff[i], 0.5f), 1.0f);
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
            mvr2r(lsp_shb_2, lsp_temp, LPC_SHB_ORDER);
        }
    }

    mvr2r( lsf_diff+1, st->prev_lsf_diff, LPC_SHB_ORDER-2 );
    st->prev_tilt_para = tilt_para;

    if( st->total_brate == ACELP_24k40 || st->total_brate == ACELP_32k )
    {
        /* SHB LSP interpolation */
        ptr_lsp_interp_coef = interpol_frac_shb;
        for( j = 0; j < 4; j++ )
        {
            for( i = 0; i < LPC_SHB_ORDER; i++ )
            {
                lsp_temp[i] = lsp_shb_1[i]*(*ptr_lsp_interp_coef)
                              + lsp_shb_2[i]*(*(ptr_lsp_interp_coef+1));
            }
            ptr_lsp_interp_coef += 2;

            tmp = j*(LPC_SHB_ORDER+1);
            /* convert from lsp to lsf */
            lsp2lsf( lsp_temp, lpc_shb_sf+tmp, LPC_SHB_ORDER, 1.0f );
            /* convert lsf to lpc for SHB synthesis */
            lsp2a( lpc_shb_sf+tmp, lpc_shb_sf+tmp, LPC_SHB_ORDER );
            lpc_shb_sf[j*(LPC_SHB_ORDER+1)] = 1.0f;
        }

        /* -------- Calculate the SHB Energy --------  */
        shb_ener_sf = 0.003125f * sum2_f( shb_frame + L_SHB_LAHEAD, L_FRAME16k );
        Quant_shb_ener_sf( st, &shb_ener_sf );

        /* --------  calculate the residuals using the FOUR subframe LPCs --------  */

        set_f( shb_res, 0, L_FRAME16k );
        residu( lpc_shb_sf,                       LPC_SHB_ORDER, shb_frame + L_SHB_LAHEAD,       shb_res,       80 );
        residu( lpc_shb_sf + (LPC_SHB_ORDER+1),   LPC_SHB_ORDER, shb_frame + L_SHB_LAHEAD + 80,  shb_res + 80,  80 );
        residu( lpc_shb_sf + 2*(LPC_SHB_ORDER+1), LPC_SHB_ORDER, shb_frame + L_SHB_LAHEAD + 160, shb_res + 160, 80 );
        residu( lpc_shb_sf + 3*(LPC_SHB_ORDER+1), LPC_SHB_ORDER, shb_frame + L_SHB_LAHEAD + 240, shb_res + 240, 80 );

        set_f( shb_res_gshape, 0, NB_SUBFR16k );
        for( i = 0; i < NB_SUBFR16k; i++ )
        {
            shb_res_gshape[i] = sum2_f( shb_res+i*64, 64 );
        }

        maximum(shb_res_gshape, NB_SUBFR16k, &normFac);
        normFac = (float)1.0f/(0.0001f + normFac);
        for( i = 0; i < NB_SUBFR16k; i++ )
        {
            shb_res_gshape[i] = sqrt(shb_res_gshape[i]*normFac);
        }

        Quant_shb_res_gshape(st, shb_res_gshape);
    }

    /* Save the SWB LSP values from current frame for interpolation */
    mvr2r( lsp_shb_2, st->swb_lsp_prev_interp, LPC_SHB_ORDER );



    /* For 13.2 and 16.4kbps, convert LSPs back into LP coeffs */
    /* convert from lsp to lsf */
    lsp2lsf( lsp_temp, lpc_shb, LPC_SHB_ORDER, 1.0f );
    /* convert lsf to lpc for SHB synthesis */
    lsp2a( lpc_shb, lpc_shb, LPC_SHB_ORDER );
    lpc_shb[0] = 1.0f;

    /* Save the SWB LSP values from current frame for interpolation */
    mvr2r( lsp_shb_2, st->swb_lsp_prev_interp, LPC_SHB_ORDER );
    mvr2r( st->state_syn_shbexc, shaped_shb_excitation, L_SHB_LAHEAD );

    /* Determine formant PF strength */
    formant_fac = swb_formant_fac( lpc_shb[1], &st->tilt_mem );

    vf_ind = 20;
    GenShapedSHBExcitation( shaped_shb_excitation + L_SHB_LAHEAD, lpc_shb, White_exc16k, st->mem_csfilt,
                            st->mem_genSHBexc_filt_down_shb, st->state_lpc_syn, coder_type, bwe_exc_extended, st->bwe_seed,
                            vf_modified, st->extl, &(st->tbe_demph), &(st->tbe_premph), lpc_shb_sf, &shb_ener_sf, shb_res_gshape,
                            shb_res, &vf_ind, formant_fac, st->fb_state_lpc_syn,&(st->fb_tbe_demph), st->total_brate, 0 );

    if( st->total_brate == ACELP_24k40 || st->total_brate == ACELP_32k )
    {
        if(st->codec_mode == MODE2)
            st->idx_mixFac = vf_ind;
        else
            push_indice( st, IND_SHB_VF, vf_ind, NUM_BITS_SHB_VF);
    }

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
        scale = sqrt( curr_pow/ prev_pow );
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
    /* Estimate the gain-shape parameter */
    EstimateSHBGainShape( SHB_OVERLAP_LEN, shb_frame, shaped_shb_excitation, GainShape, subwin_shb );

    /* Gain shape BWS/high band low energy fix */
    if( st->cldfbHBLT < 1.0f )
    {
        /* There is not much HB past 10kHz; the temporal resolution is quite coarse, so reduce the dynamic range */
        for(i = 0; i < NUM_SHB_SUBGAINS; i++)
        {
            /* 0.35f +/- delta variation; delta = 0.1*(GS-0.35)*/
            GainShape[i] = 0.315f + 0.1f * GainShape[i];
        }
    }
    /* Gain frame adjustment factor */
    temp_swb_fac = (float)log( (GainShape[0]+0.00001f) / (st->prev_swb_GainShape+0.0001f) );
    feedback = temp_swb_fac * temp_swb_fac ;
    for (i = 1; i < NUM_SHB_SUBGAINS; i++)
    {
        temp_swb_fac = (float)log( (GainShape[i]+0.00001f) / (GainShape[i-1]+0.0001f) );
        feedback += (temp_swb_fac * temp_swb_fac);
    }
    feedback = 0.4f / (1 + 0.5f * feedback);

    if( frGainAttenuate == 1 || ( sum_f(voicing, 3) > 2.4f && sum_f(voice_factors, 4) > 0.8f ) )
    {
        temp_swb_fac = st->prev_swb_GainShape;
        for( i = 0; i < NUM_SHB_SUBGAINS; i++ )
        {
            GainShape[i] = (1 - feedback) * GainShape[i] + feedback * temp_swb_fac;
            temp_swb_fac = GainShape[i];
        }
    }
    st->prev_swb_GainShape = GainShape[3];

    /* Compute the power of gains away from the peak gain prior to quantization */
    p2m_in = pow_off_pk(GainShape, NUM_SHB_SUBGAINS, 1);

    /* Quantization of the gain shape parameter */
    QuantizeSHBsubgains( st, GainShape, st->extl );

    /* Compute the power of gains away from the peak gain after quantization */
    p2m_out = pow_off_pk(GainShape, NUM_SHB_SUBFR, 4);

    /* Estimate the gain parameter */
    EstimateSHBFrameGain( SHB_OVERLAP_LEN, shb_frame, shaped_shb_excitation, GainShape, &GainFrame, window_shb, subwin_shb );

    if( st->tec_tfa == 1 )
    {
        tfaCalcEnv( shb_frame, st->tfa_enr );
    }

    /* If there's a big difference in the power of gains away from the peak gain   */
    /* due to poor quantization then suppress energy of the high band.             */
    if( p2m_out > 2.0f * p2m_in )
    {
        float temp = 0;
        if(p2m_in >= 0 && p2m_out > 0)
        {
            temp = sqrt((2.0f * p2m_in)/p2m_out);
        }
        GainFrame *= temp;
    }

    if( frGainSmoothEn == 1 && st->prev_gainFr_SHB < GainFrame )
    {
        GainFrame = 0.5f * (st->prev_gainFr_SHB + GainFrame);
    }

    if( frGainAttenuate == 1 && MA_lsp_shb_spacing <= 0.0024f )
    {
        GainFrame = (float)pow( GainFrame, 0.8f );
    }
    else if( st->prev_frGainAtten == 1 && GainFrame > 3.0f * st->prev_gainFr_SHB )
    {
        GainFrame *= (0.8f + 0.5f*feedback);
    }
    st->prev_frGainAtten = frGainAttenuate;

    st->prev_gainFr_SHB = GainFrame;

    /* Gain attenuation when the SWB LSF quantization error is larger than a threshold */
    sd_uq_q = (sd_uq_q/0.0025f);
    if( st->L_frame == L_FRAME )
    {
        sd_uq_q = 1 - 0.2f* (sd_uq_q * sd_uq_q);
    }
    else
    {
        sd_uq_q = 1 - 0.1f* (sd_uq_q * sd_uq_q);
    }

    sd_uq_q =  max(min(sd_uq_q, 1.0f), 0.5f);

    pitBufAvg =  0.0025f * sum_f(pitch_buf, 4);
    voicingBufAvg = (sum_f(voice_factors, 4) > 0.6f) ? 0.333f : 0.1667f;
    voicingBufAvg = voicingBufAvg * sum_f(voicing, 3);

    if(voicingBufAvg <= 0.0f && (sd_uq_q*pitBufAvg != 0))
    {
        voicingBufAvg = sd_uq_q*pitBufAvg/1.001f;
    }
    else if(voicingBufAvg <= 0.0f)
    {
        voicingBufAvg = 1.0f;
    }

    /* Controlled gain evolution in SWB for stronger voiced segments */
    GainFrame *= max(min((float)(sd_uq_q*pitBufAvg/voicingBufAvg), 1.0f), 0.6f);
    if( st->L_frame == L_FRAME16k || st->rf_mode == 1 )
    {
        /* Compensate for energy increase mismatch due to memory-less synthesis*/
        GainFrame *= 0.85f;
    }

    /* Quantization of the frame gain parameter */
    QuantizeSHBframegain( st, &GainFrame, st->extl, 0, &st->RF_bwe_gainFr_ind );

    /* Adjust the subframe and frame gain of the synthesized SHB signal */
    /* Scale the shaped excitation */


    if( st->extl == FB_TBE)
    {
        for( i=0; i<L_FRAME16k; i++ )
        {
            White_exc16k[i] *= GainFrame * GainShape[NUM_SHB_SUBFR*i/L_FRAME16k];
        }
    }


    return;
}


/*-------------------------------------------------------------------*
 * EstimateSHBFrameGain()
 *
 * Estimate the overall gain factor needed to scale synthesized highband
 * to original highband signal level.
 *-------------------------------------------------------------------*/

static void EstimateSHBFrameGain(
    const short length,             /* i  : SHB overlap length          */
    const float *oriSHB,            /* i  : target original SHB frame   */
    const float *synSHB,            /* i  : shaped SHB excitation       */
    float *subgain,           /* o  : estimate of gain shape      */
    float *GainFrame,         /* o  : estimat of gain frame       */
    const float *win_shb,           /* i  : SHB window                  */
    const float *subwin_shb         /* i  : SHB subframe window         */
)
{
    const short *skip;
    short i, j, k, l_shb_lahead, l_frame;
    short join_length, num_join;
    float sig,mod_syn[L_FRAME16k+L_SHB_LAHEAD];
    float oriNrg = 0.0f, synNrg = 0.0f, sum_gain;
    float frame_gain;

    /* initilaization */
    l_frame = L_FRAME16k;
    l_shb_lahead = L_SHB_LAHEAD;
    skip = skip_bands_SWB_TBE;

    if( length == SHB_OVERLAP_LEN/2 )
    {
        skip = skip_bands_WB_TBE;
        l_frame = L_FRAME16k/4;
        l_shb_lahead = L_SHB_LAHEAD/4;
    }

    /* apply gain for each subframe, and store noise output signal using overlap-add*/
    set_f( mod_syn, 0, l_frame+l_shb_lahead );

    if ( length == SHB_OVERLAP_LEN/2 )
    {
        sum_gain = 0;
        for( k=0; k<length/2; k++ )
        {
            sum_gain = subwin_shb[2*k+2]*subgain[0];
            mod_syn[skip[0]+k] = synSHB[skip[0]+k] * sum_gain;
            mod_syn[skip[0]+k+length/2] = synSHB[skip[0]+k+length/2] * subgain[0];
        }

        for( i=1; i<NUM_SHB_SUBFR/2; i++ )
        {
            for( k=0; k<length; k++ )
            {
                sum_gain = subwin_shb[k+1]*subgain[i] + subwin_shb[length-k-1]*subgain[i-1];
                mod_syn[skip[i]+k] = synSHB[skip[i]+k] * sum_gain;
            }
        }

        for( k=0; k<length/2; k++ )
        {
            sum_gain = subwin_shb[length-2*k-2]*subgain[i-1];
            mod_syn[skip[i]+k] = synSHB[skip[i]+k] * sum_gain;
        }
    }
    else
    {
        num_join = NUM_SHB_SUBFR/NUM_SHB_SUBGAINS;
        join_length = num_join*length;
        for (k = 0, j = 0; k < length; k++)
        {
            mod_syn[j] = synSHB[j]*subwin_shb[k+1]*subgain[0];
            j++;
        }
        for (i = 0; i < NUM_SHB_SUBGAINS-1; i++)
        {
            for (k = 0; k < join_length - length; k++)
            {
                mod_syn[j] = synSHB[j]*subgain[i*num_join];
                j++;
            }

            for (k = 0; k < length; k++)
            {
                mod_syn[j] = synSHB[j]*(subwin_shb[length-k-1]*subgain[i*num_join] + subwin_shb[k+1]*subgain[(i+1)*num_join]);
                j++;
            }
        }
        for (k = 0; k < join_length - length; k++)
        {
            mod_syn[j] = synSHB[j]*subgain[(NUM_SHB_SUBGAINS-1)*num_join];
            j++;
        }
        for (k = 0; k < length; k++)
        {
            mod_syn[j] = synSHB[j]*subwin_shb[length-k-1]*subgain[(NUM_SHB_SUBGAINS-1)*num_join];
            j++;
        }
    }

    /* adjust frame energy */
    oriNrg = (float)1e-10;
    synNrg = (float)1e-10;

    for( i=0; i<l_shb_lahead; i++ )
    {
        sig = oriSHB[i]*win_shb[i];
        oriNrg += sig*sig;

        sig = mod_syn[i]*win_shb[i];
        synNrg += sig*sig;
    }

    for( ; i<l_frame; i++ )
    {
        oriNrg += oriSHB[i]*oriSHB[i];
        synNrg += mod_syn[i]*mod_syn[i];
    }

    for( ; i<l_frame+l_shb_lahead; i++ )
    {
        sig = oriSHB[i]*win_shb[l_frame+l_shb_lahead-1-i];
        oriNrg += sig*sig;

        sig = mod_syn[i]*win_shb[l_frame+l_shb_lahead-1-i];
        synNrg += sig*sig;
    }

    frame_gain = sqrt( oriNrg / synNrg );
    if ((synNrg)==0) frame_gain = 0;

    *GainFrame = (float)frame_gain;

    return;
}

/*-------------------------------------------------------------------*
 * pow_off_pk()
 *
 * Sums squares of SWB shape gain parameters away from peak values
 *-------------------------------------------------------------------*/

static float pow_off_pk(
    float a[],
    short len,
    short step
)
{
    short i, j = 0;
    float sum;

    sum = a[0] * a[0];
    for (j=0, i=1; i<len; i+=step)
    {
        sum += a[i] * a[i];
        if (a[i] > a[j])
        {
            j = i;
        }
    }

    sum -= a[j] * a[j];

    return (sum);
}


/*-------------------------------------------------------------------*
 * EstimateSHBGainShape()
 *
 * Estimate temporal gain parameters
 *-------------------------------------------------------------------*/

static void EstimateSHBGainShape(
    const short length,             /* i  : SHB overlap length          */
    const float *oriSHB,            /* i  : target original SHB frame   */
    const float *synSHB,            /* i  : shaped SHB excitation       */
    float *subgain,           /* o  : estimate of gain shape      */
    const float *subwin             /* i  : SHB subframe window         */
)
{
    const short *skip;
    short i, k;
    float sum_gain, oriNrg, synNrg, sig;
    short num_join, num_gains, join_length;

    float normFact;

    num_join = NUM_SHB_SUBFR/NUM_SHB_SUBGAINS;
    num_gains = NUM_SHB_SUBGAINS;
    skip = skip_bands_SWB_TBE;

    if( length == SHB_OVERLAP_LEN/2 )
    {
        num_gains = NUM_SHB_SUBFR/4;
        skip = skip_bands_WB_TBE;
    }

    /* calculate and normalize the subgain */
    sum_gain = (float)1e-10;
    oriNrg = 0.0f;
    synNrg = 0.0f;

    if( length == SHB_OVERLAP_LEN/2 )
    {
        for( i=0; i<NUM_SHB_SUBFR/2; i++ )
        {
            if( (i & 0x1) == 0 )
            {
                oriNrg = 1e-10f;
                synNrg = 1e-10f;
            }

            if( i == 0)
            {
                for( k=0; k<length/2; k++ )
                {
                    sig = oriSHB[skip[i]+k] * subwin[2*k+2];
                    oriNrg += sig*sig;

                    sig = synSHB[skip[i]+k] * subwin[2*k+2];
                    synNrg += sig*sig;
                }

                for( k=length/2; k<length; k++ )
                {
                    sig = oriSHB[skip[i]+k];
                    oriNrg += sig*sig;

                    sig = synSHB[skip[i]+k];
                    synNrg += sig*sig;
                }
            }
            else
            {
                for( k=0; k<length; k++ )
                {
                    sig = oriSHB[skip[i]+k] * subwin[k+1];
                    oriNrg += sig*sig;

                    sig = synSHB[skip[i]+k] * subwin[k+1];
                    synNrg += sig*sig;
                }
            }

            if( i == NUM_SHB_SUBFR/2 - 1 )
            {
                for( ; k<2*length - length/2; k++ )
                {
                    sig = oriSHB[skip[i]+k] * subwin[3*length-2*k-2];
                    oriNrg += sig*sig;

                    sig = synSHB[skip[i]+k] * subwin[3*length-2*k-2];
                    synNrg += sig*sig;
                }
            }
            else
            {
                for( ; k<2*length; k++ )
                {
                    sig = oriSHB[skip[i]+k] * subwin[2*length-k-1];
                    oriNrg += sig*sig;

                    sig = synSHB[skip[i]+k] * subwin[2*length-k-1];
                    synNrg += sig*sig;
                }
            }

            if ( (i & 0x1) == 1 )
            {
                subgain[i/2] = sqrt( oriNrg/synNrg );
                if ((synNrg)==0) subgain[i/2] = 0;

                sum_gain += subgain[i/2] * subgain[i/2];
            }
        }
    }
    else
    {
        join_length = num_join*length;

        for (i = 0; i < num_gains; i++)
        {
            oriNrg = 1e-10f;
            synNrg = 1e-10f;

            for (k = 0; k < length; k++)
            {
                sig = oriSHB[join_length*i + k] * subwin[k + 1];
                oriNrg += sig*sig;

                sig = synSHB[join_length*i + k] * subwin[k + 1];
                synNrg += sig*sig;
            }

            for (k = 0; k < (join_length - length); k++)
            {
                sig = oriSHB[length + join_length*i + k];
                oriNrg += sig*sig;

                sig = synSHB[length + join_length*i + k];
                synNrg += sig*sig;
            }

            for (k = 0; k < length; k++)
            {
                sig = oriSHB[join_length*(i+1) + k] * subwin[length - k - 1];
                oriNrg += sig*sig;

                sig = synSHB[join_length*(i+1) + k] * subwin[length - k - 1];
                synNrg += sig*sig;
            }
            subgain[i] = sqrt( oriNrg/synNrg );
            if ((synNrg)==0) subgain[i] = 0;
            sum_gain += subgain[i] * subgain[i];
        }
    }

    /* normalize the subgain */
    normFact = sqrt( 1.0f / sum_gain );
    if (sum_gain==0) normFact = 0;
    for (i = 0; i < num_gains; i++)
    {
        subgain[i] *= normFact;
    }



    return;
}

/*-------------------------------------------------------------------*
 * Quant_lower_LSF()
 *
 * Quantize the lower half of the LSF vector
 *-------------------------------------------------------------------*/

static void Quant_lower_LSF(
    const float lsf[],                  /* i  : Input LSFs             */
    float lsf_q[],                /* o  : Quantized LSFs         */
    short lsf_idx[]               /* o  : Quantized LSFs indices */
)
{
    short i;

    lsf_idx[0] = (short)squant(lsf[0], &lsf_q[0], lsf_q_cb[0], lsf_q_cb_size[0]);
    for (i = 1; i < NUM_Q_LSF; i++)
    {
        lsf_idx[i] = (short)squant(lsf[i] - lsf_q[i - 1], &lsf_q[i], lsf_q_cb[i], lsf_q_cb_size[i]);
        lsf_q[i] += lsf_q[i-1];
    }

    return;
}

/*-------------------------------------------------------------------*
 * Quant_mirror_point()
 *
 * Quantize the mirror point
 *-------------------------------------------------------------------*/

static short Quant_mirror_point(
    const float lsf[],                  /* i  : Input LSFs          */
    const float lsf_q[],
    float *m                      /* o  : Mirror point        */
)
{
    float m_diff;
    short m_idx;

    m_diff = 0.5f * (lsf[NUM_Q_LSF] - lsf_q[NUM_Q_LSF - 1]);
    m_idx = (short)squant(m_diff, m, mirror_point_q_cb, MIRROR_POINT_Q_CB_SIZE);

    *m += lsf_q[NUM_Q_LSF - 1];

    return m_idx;
}

/*-------------------------------------------------------------------*
 * Find_LSF_grid()
 *
 * Find the best grid for the LSFs
 *-------------------------------------------------------------------*/

static short Find_LSF_grid(
    const float lsf[],                  /* i  : Input LSFs             */
    float lsf_q[],                /* o  : Quantized LSFs         */
    const float m                       /* i  : Mirror point           */
)
{
    float lsf_map[NUM_MAP_LSF];
    float grid[NUM_LSF_GRIDS][NUM_MAP_LSF];
    float offset;
    float last_q_lsf;
    float lsf_t[NUM_MAP_LSF];
    float lsf_smooth[NUM_MAP_LSF];
    float D, D_best;
    short I_best = 0;
    short i, j;
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

    for (i = 0; i < NUM_LSF_GRIDS; i++)
    {
        for (j = 0; j < NUM_MAP_LSF; j++)
        {
            grid[i][j] = lsf_grid[i][j]*scale + last_q_lsf;
        }
    }

    D_best = QUANT_DIST_INIT;
    for (i = 0; i < NUM_LSF_GRIDS; i++)
    {
        D = EPSILON;
        for (j = 0; j < NUM_MAP_LSF; j++)
        {
            lsf_t[j] = (1 - grid_smoothing[j])*lsf_map[j] + grid_smoothing[j]*grid[i][j];
            D += (lsf_t[j] - lsf[NUM_Q_LSF + j])*(lsf_t[j] - lsf[NUM_Q_LSF + j]);
        }

        if (D < D_best)
        {
            mvr2r(lsf_t, lsf_smooth, NUM_MAP_LSF);
            D_best = D;
            I_best = i;
        }
    }

    for (i = 0; i < NUM_MAP_LSF; i++)
    {
        lsf_q[NUM_Q_LSF + i] = lsf_smooth[i];
    }

    return I_best;
}

/*-------------------------------------------------------------------*
* gainFrSmooth_En()
*
* Gain frame smoothing and attenuation control
*-------------------------------------------------------------------*/
static void gainFrSmooth_En(
    Encoder_State *st,
    float *shb_frame,
    const float *lpc_shb,
    const float *lsp_shb,
    float *MA_lsp_shb_spacing,
    short *frGainAttenuate,
    short *frGainSmoothEn
)
{
    float lsp_slow_evol_rate = 0, lsp_fast_evol_rate = 0, lsp_spacing = 0.5f;
    float temp_shb_frame[L_FRAME16k+L_SHB_LAHEAD];
    int i;

    for( i = 0; i < LPC_SHB_ORDER; i++ )
    {
        lsp_spacing = min(lsp_spacing, (float)( i == 0 ? lsp_shb[0] : (lsp_shb[i] - lsp_shb[i -1])));

        /* estimate the mean square error in lsps from current frame to past frames */
        lsp_slow_evol_rate += (lsp_shb[i] - st->lsp_shb_slow_interpl[i]) * (lsp_shb[i] - st->lsp_shb_slow_interpl[i]);
        lsp_fast_evol_rate += (lsp_shb[i] - st->lsp_shb_fast_interpl[i]) * (lsp_shb[i] - st->lsp_shb_fast_interpl[i]);

        /* update the slow and fast interpolation lsps for next frame */
        st->lsp_shb_slow_interpl[i] = 0.7f * st->lsp_shb_slow_interpl[i] + 0.3f * lsp_shb[i];
        st->lsp_shb_fast_interpl[i] = 0.3f * st->lsp_shb_fast_interpl[i] + 0.7f * lsp_shb[i];
    }

    if( st->last_extl != SWB_TBE && st->last_extl != FB_TBE && lsp_spacing < 0.008f )
    {
        st->lsp_shb_spacing[0] = lsp_spacing;
        st->lsp_shb_spacing[1] = lsp_spacing;
        st->lsp_shb_spacing[2] = lsp_spacing;
        st->prev_frGainAtten = 1;
    }

    *MA_lsp_shb_spacing = 0.1f*st->lsp_shb_spacing[0] + 0.2f*st->lsp_shb_spacing[1] + 0.3f*st->lsp_shb_spacing[2] + 0.4f*lsp_spacing;

    st->lsp_shb_spacing[0] = st->lsp_shb_spacing[1];
    st->lsp_shb_spacing[1] = st->lsp_shb_spacing[2];
    st->lsp_shb_spacing[2] = lsp_spacing;

    *frGainAttenuate = 0;
    *frGainSmoothEn = 0;

    if( (lsp_spacing < 0.008f && (*MA_lsp_shb_spacing < 0.005f || st->prev_frGainAtten == 1)) || lsp_spacing <= 0.0032f )
    {
        *frGainAttenuate = 1;
        mvr2r(shb_frame, temp_shb_frame, L_FRAME16k+L_SHB_LAHEAD);

        fir( temp_shb_frame, lpc_shb, shb_frame, st->shb_inv_filt_mem, L_FRAME16k+L_SHB_LAHEAD, LPC_SHB_ORDER, 1 );

        if( lsp_slow_evol_rate < 0.001f && lsp_fast_evol_rate < 0.001f )
        {
            *frGainSmoothEn = 1;
        }
    }
}

/*-------------------------------------------------------------------*
 * Quant_BWE_LSF()
 *
 * Quantize super highband spectral envolope
 *-------------------------------------------------------------------*/

static void Quant_BWE_LSF(
    Encoder_State *st,                    /* i/o: encoder state structure      */
    const float lsf_shb[],              /* i  : unquanitzed LSFs */
    float Q_lsfs[]                      /* o  : quanitzed LSFs   */
)
{
    float lsf[LPC_SHB_ORDER];
    float lsf_q[LPC_SHB_ORDER];
    short lsf_idx[NUM_Q_LSF];
    short i;
    short m_idx;
    float m;
    short grid_idx;

    for (i = 0; i < LPC_SHB_ORDER; i++)
    {
        lsf[i] = 0.5f - lsf_shb[LPC_SHB_ORDER - 1 - i];
    }

    Quant_lower_LSF( lsf, lsf_q, lsf_idx );

    for (i = 0; i < NUM_Q_LSF; i++)
    {
        if(st->codec_mode == MODE2)
        {
            st->lsf_idx[i] = lsf_idx[i];
        }
        else
        {
            push_indice( st, IND_SHB_LSF, lsf_idx[i], lsf_q_num_bits[i] );
        }
    }

    m_idx = Quant_mirror_point(lsf, lsf_q, &m);

    if(st->codec_mode == MODE2)
        st->m_idx = m_idx;
    else
        push_indice( st,IND_SHB_MIRROR, m_idx, MIRROR_POINT_BITS);

    grid_idx = Find_LSF_grid(lsf, lsf_q, m);

    if(st->codec_mode == MODE2)
        st->grid_idx = grid_idx;
    else
        push_indice( st,IND_SHB_GRID, grid_idx, NUM_LSF_GRID_BITS);

    for (i = 0; i < LPC_SHB_ORDER; i++)
    {
        Q_lsfs[i] = 0.5f - lsf_q[LPC_SHB_ORDER - 1 - i];
    }

    return;
}

/*-------------------------------------------------------------------*
 * closest_centroid()
 *
 * Determine a set of closest VQ centroids for a given input
 *-------------------------------------------------------------------*/

static short closest_centroid(
    const float *data,            /* i  : input data */
    const float *weights,         /* i  : weights */
    const float *quantizer,       /* i  : quantizer table */
    const short  centroids,        /* i  : number of centroids */
    const short   length)           /* i  : dimension of quantiser */
{
    short i,j, index;
    float tmp, werr, best_werr;

    index = 0;
    best_werr = 1.0E20f;

    for( i=0; i<centroids; i++ )
    {
        werr = 0.0f;
        for( j=0; j<length; j++ )
        {
            tmp = (float) *(data + j) - quantizer[i * length + j];
            werr += (float) (*(weights + j) * tmp * tmp);
            if(werr > best_werr) break;
        }

        if( werr < best_werr )
        {
            best_werr = werr;
            index = i;
        }
    }

    return index;
}

/*-------------------------------------------------------------------*
 * QuantizeSHBsubgains()
 *
 * Quantize super highband temporal gains
 *-------------------------------------------------------------------*/

static void QuantizeSHBsubgains(
    Encoder_State *st,                /* i/o: encoder state structure      */
    float subgains[],         /* i/o:  super highband temporal gains  */
    const short extl                /* i  : extension layer                 */
)
{
    short i, idxSubGain;
    float Unit_weights10[NUM_SHB_SUBFR];

    if( extl == WB_TBE )
    {
        set_f( Unit_weights10, 1.0f, (short)NUM_SHB_SUBFR/4 );

        for( i=0; i<NUM_SHB_SUBFR/4; i++ )
        {
            subgains[i+NUM_SHB_SUBFR/4] = 20.0f * (float) log10(subgains[i]);
        }

        idxSubGain = closest_centroid( subgains+NUM_SHB_SUBFR/4, Unit_weights10,
                                       HBCB_SubGain5bit, 1<<NUM_BITS_SHB_SUBGAINS,  NUM_SHB_SUBFR/4 );

        mvr2r( HBCB_SubGain5bit + idxSubGain * NUM_SHB_SUBFR/4, subgains, NUM_SHB_SUBFR/4 );

        push_indice( st, IND_SHB_SUBGAIN, idxSubGain, NUM_BITS_SHB_SUBGAINS );
        for( i=0; i<NUM_SHB_SUBFR/4; i++ )
        {
            subgains[i] = (float) pow(10.0f, subgains[i] / 20.0f);
        }

        for( i=NUM_SHB_SUBFR/2-1; i>=0; i-- )
        {
            subgains[i] = subgains[i/2];
        }
    }
    else
    {
        for(i = 0; i < NUM_SHB_SUBGAINS; i++)
        {
            subgains[i] = (float) log10(subgains[i] + 0.001f);
        }

        idxSubGain = (short)vquant(subgains, 0, subgains, SHBCB_SubGain5bit, NUM_SHB_SUBGAINS, 1<<NUM_BITS_SHB_SUBGAINS);

        for(i = 0; i < NUM_SHB_SUBGAINS; i++)
        {
            subgains[i] = (float) pow(10.0, subgains[i]);
        }

        for (i = NUM_SHB_SUBFR-1; i >= 0; i--)
        {
            subgains[i] = subgains[i*NUM_SHB_SUBGAINS/NUM_SHB_SUBFR];
        }

        if( st->codec_mode == MODE2 )
        {
            st->idxSubGains = idxSubGain;
        }
        else
        {
            push_indice(st, IND_SHB_SUBGAIN, idxSubGain, NUM_BITS_SHB_SUBGAINS);
        }
    }

    return;
}

/*-------------------------------------------------------------------*
 * Quant_shb_ener_sf()
 *
 * Quantize SHB subframe energies
 *-------------------------------------------------------------------*/

static void Quant_shb_ener_sf(
    Encoder_State *st,                /* i/o: encoder state structure      */
    float *shb_ener_sf         /* i/o:  super highband subframe energies  */
)
{
    short idxSubEner;

    float temp_shb_ener_sf, sum;

    sum = *shb_ener_sf;
    *shb_ener_sf = (float)log10(sum);

    temp_shb_ener_sf = 0;
    idxSubEner  = usquant(*shb_ener_sf, &temp_shb_ener_sf, 0, 0.042f, 1<<NUM_BITS_SHB_ENER_SF);
    *shb_ener_sf = (float)pow(10.0, temp_shb_ener_sf );

    if(st->codec_mode == MODE2)
        st->idx_shb_fr_gain = idxSubEner;
    else
        push_indice( st, IND_SHB_ENER_SF, idxSubEner, NUM_BITS_SHB_ENER_SF);

    return;
}

/*-------------------------------------------------------------------*
* Quant_shb_res_gshape()
*
* Quantize SHB gain shapes in residual domain
*-------------------------------------------------------------------*/

static void Quant_shb_res_gshape(
    Encoder_State *st,                     /* i/o: encoder state structure      */
    float shb_res_gshape[]         /* i/o:  super highband gain shapes */
)
{
    short i, idxSubGain[NB_SUBFR16k];

    for(i = 0; i < NB_SUBFR16k; i++)
    {
        idxSubGain[i] = usquant(shb_res_gshape[i], &shb_res_gshape[i], 0.125, 0.125f, 1<<NUM_BITS_SHB_RES_GS);
        if(st->codec_mode == MODE2)
            st->idx_res_gs[i] = idxSubGain[i];
        else
            push_indice( st, IND_SHB_RES_GS1+i, idxSubGain[i], NUM_BITS_SHB_RES_GS);
    }
}

/*-------------------------------------------------------------------*
 * QuantizeSHBframegains()
 *
 * Quantize super highband frame gain
 *-------------------------------------------------------------------*/

static void QuantizeSHBframegain(
    Encoder_State *st,                /* i/o: encoder state structure      */
    float *GainFrame,         /* i/o: Gain                            */
    const short extl,               /* i  : extension layer                 */
    long  extl_brate          /* i  : extension layer bitrate         */
    ,short *rf_gainFrame_ind
)
{
    short idxFrameGain;
    float Q_GainFrame;
    float Unit_weights1 = 1.0f;

    float GainFrameLog;

    if ( extl == WB_TBE )
    {
        determine_gain_weights( GainFrame, &(Unit_weights1), 1 );

        if( extl_brate == WB_TBE_0k35 )
        {
            singlevectortest_gain( GainFrame, 1, 1<<NUM_BITS_SHB_FrameGain_LBR_WB, &idxFrameGain,
                                   &(Unit_weights1), &Q_GainFrame, SHBCB_FrameGain16 );

            if ( Q_GainFrame > *GainFrame * 1.06f && idxFrameGain > 0 )  /* 1.06 = +0.5 dB */
            {
                idxFrameGain--;
                Q_GainFrame = SHBCB_FrameGain16[idxFrameGain];
            }

            st->gFrame_WB = idxFrameGain;
            *rf_gainFrame_ind = idxFrameGain;
        }
        else
        {
            singlevectortest_gain( GainFrame, 1, 1<<NUM_BITS_SHB_FrameGain, &idxFrameGain,
                                   &(Unit_weights1), &Q_GainFrame, SHBCB_FrameGain64 );

            push_indice( st, IND_SHB_FRAMEGAIN, idxFrameGain, NUM_BITS_SHB_FrameGain );
            *rf_gainFrame_ind = idxFrameGain;
        }
    }
    else
    {
        GainFrameLog = (float) log10(*GainFrame + 0.001f);
        idxFrameGain = (short)usquant(GainFrameLog, &Q_GainFrame, SHB_GAIN_QLOW, SHB_GAIN_QDELTA, 1<<NUM_BITS_SHB_FRAMEGAIN);

        while( Q_GainFrame > GainFrameLog + 0.495*SHB_GAIN_QDELTA  && idxFrameGain != 0)
        {
            idxFrameGain = idxFrameGain - 1;
            Q_GainFrame = idxFrameGain*SHB_GAIN_QDELTA + SHB_GAIN_QLOW;
        }

        Q_GainFrame = (float) pow(10.0, Q_GainFrame );
        if( st->codec_mode == MODE2 )
        {
            st->idxFrameGain = idxFrameGain;
        }
        else
        {
            push_indice( st, IND_SHB_FRAMEGAIN, idxFrameGain, NUM_BITS_SHB_FRAMEGAIN );
        }
        *rf_gainFrame_ind = idxFrameGain;
    }
    if( st->rf_mode )
    {
        /*Currently intended for SWB only. Modify for WB is needed later!*/
        if( st->rf_frame_type == RF_NELP )
        {
            *rf_gainFrame_ind = idxFrameGain; /* NELP Frame uses full 5 bits */
        }
        else /*RF_ALLPRED, RF_GENPRED, RF_NOPRED modes*/
        {
            if( *GainFrame <= 1.25 )     /* [0 to 1.25] range --> 0.5*/
                *rf_gainFrame_ind = 0;
            else if ( *GainFrame <= 3 )  /* (1.25 to 3] --> 2 */
                *rf_gainFrame_ind = 1;
            else if ( *GainFrame <= 6 )  /* (3 to 6] --> 4 */
                *rf_gainFrame_ind = 2;
            else                        /* (6 to Inf) --> 8 */
                *rf_gainFrame_ind = 3;
        }
    }

    *GainFrame = Q_GainFrame;

    return;
}



/*-------------------------------------------------------------------*
 * determine_gain_weights()
 *
 * Determine weights for gain quantization
 *-------------------------------------------------------------------*/

static void determine_gain_weights (
    const float *gain,        /* i  : Gain parameter    */
    float *weights,     /* o  : gain weights      */
    const short dims          /* i  : number of gains   */
)
{
    short j;

    for( j = 0; j < dims; j++ )
    {
        if( gain[j] > 1e-6 )
        {
            weights[j] = (float)(pow (fabs (gain[j]), -0.9f));
        }
        else
        {
            weights[j] = 501.187233628f;
        }
    }

    return;
}

/*-------------------------------------------------------------------*
 * singlevectortest_gain()
 *
 * VQ for coding superhigh band gain
 *-------------------------------------------------------------------*/

static void singlevectortest_gain (
    const float *inp,           /* i  : input gain vector               */
    const short dimen,          /* i  : dimension of the input vector   */
    const short cb_size,        /* i  : codebook size                   */
    short *index,         /* o  : quanitzation index              */
    const float *weight,        /* i  : Weights for the quanitzation    */
    float *recon,         /* o  : Reconstruction                  */
    const float *codebook       /* i  : Codebook                        */
)
{
    short k, interNum, flag;
    float meanU, meanQ;
    short least[4];

    interNum = 4;

    return_M_Least( inp, dimen, codebook, cb_size, weight, interNum, least );

    meanU = sum_f( inp, dimen );
    mvr2r( codebook + dimen*least[0], recon, dimen );

    index[0] = least[0];
    flag = 0;
    for( k = 0; k < interNum; k++ )
    {
        if( flag == 0 )
        {
            meanQ = sum_f( codebook + dimen*least[k], dimen );

            if( meanQ <= 1.1 * meanU )
            {
                flag = 1;
                mvr2r( codebook + dimen*least[k], recon, dimen );
                index[0] = least[k];
            }
        }
    }

    return;
}

static void return_M_Least(
    const float *inp,           /* i: input                                         */
    const short n_cols,         /* i: vector size                                   */
    const float *codebook,      /* i: codebook                                      */
    const short num_grp,        /* i: number of centroids                           */
    const float *weight,        /* i: gain weights                                  */
    const short interNum,       /* i: number on short list prior to 2nd stage search*/
    short *least          /* o: return value                                  */
)
{
    short i, k;
    float distance[1024], mindist, tmp;

    mindist = QUANT_DIST_INIT;
    for( i = 0; i < num_grp; i++ )
    {
        distance[i] = 0;
        for( k = 0; k < n_cols; k++ )
        {
            tmp = inp[k] - codebook[n_cols*i + k];
            distance[i] += weight[k] * tmp * tmp;
        }

        if( distance[i] < mindist )
        {
            mindist = distance[i];
            least[0] = i;
        }
    }

    distance[least[0]] = QUANT_DIST_INIT;

    for( k = 1; k < interNum; k++ )
    {
        mindist = QUANT_DIST_INIT;
        for( i = 0; i < num_grp; i++ )
        {
            if( distance[i] < mindist )
            {
                mindist = distance[i];
                least[k] = i;
            }
        }

        distance[least[k]] = QUANT_DIST_INIT;
    }

    return;
}


/*-------------------------------------------------------------------*
 * fb_tbe_reset_enc()
 *
 * Reset the extra parameters needed for  FB TBE encoding
 *-------------------------------------------------------------------*/

void fb_tbe_reset_enc(
    float elliptic_bpf_2_48k_mem[][4],
    float *prev_fb_energy
)
{
    set_f( elliptic_bpf_2_48k_mem[0], 0.0f, 4 );
    set_f( elliptic_bpf_2_48k_mem[1], 0.0f, 4 );
    set_f( elliptic_bpf_2_48k_mem[2], 0.0f, 4 );
    set_f( elliptic_bpf_2_48k_mem[3], 0.0f, 4 );
    *prev_fb_energy = 0.0f;

    return;
}

/*-------------------------------------------------------------------*
 * fb_tbe_enc()
 *
 * FB TBE encoder, 14(resp. 15.5) - 20 kHz band encoding module
 *-------------------------------------------------------------------*/

void fb_tbe_enc(
    Encoder_State *st,                /* i/o: encoder state structure                 */
    const float new_input[],        /* i  : input speech at 48 kHz sample rate      */
    const float fb_exc[]            /* i  : FB excitation from the SWB part         */
)
{
    float fb_exc_energy, ratio, temp2;
    float tmp_vec[L_FRAME48k];
    short idxGain;
    float input_fhb[L_FRAME48k];
    short Sample_Delay_HP;

    elliptic_bpf_48k_generic( new_input, tmp_vec, st->elliptic_bpf_2_48k_mem, full_band_bpf_2 );

    Sample_Delay_HP = NS2SA(48000, ACELP_LOOK_NS + DELAY_FD_BWE_ENC_12k8_NS + DELAY_FIR_RESAMPL_NS) - L_FRAME48k/2 ;

    if( st->last_extl != FB_TBE )
    {
        set_f( st->old_input_fhb, 0.0f, Sample_Delay_HP );
        set_f(tmp_vec, 0.0f, 320);
    }

    mvr2r( st->old_input_fhb, input_fhb, Sample_Delay_HP );
    mvr2r( tmp_vec, input_fhb + Sample_Delay_HP, L_FRAME48k-Sample_Delay_HP );
    mvr2r( tmp_vec + L_FRAME48k - Sample_Delay_HP, st->old_input_fhb, Sample_Delay_HP );

    /* Compute the energy of the Fullband component over 4kHz (16kHz to 20kHz) */
    temp2 = sum2_f( input_fhb, L_FRAME48k/2 ) + st->prev_fb_energy;
    st->prev_fb_energy = sum2_f( input_fhb + L_FRAME48k/2, L_FRAME48k/2 );
    fb_exc_energy = sum2_f( fb_exc, L_FRAME16k ) + EPSILON;
    ratio = (float) sqrt( temp2 / fb_exc_energy );
    idxGain = (short)( log2_f ((float)ratio) + 0.5f );
    idxGain = max( 0, min(15,idxGain) );
    ratio = (float)(1 << idxGain);

    if( st->codec_mode == MODE2 )
    {
        st->idxGain = idxGain;
    }
    else
    {
        push_indice( st, IND_FB_SLOPE, idxGain, 4 );
    }


    return;
}


/*---------------------------------------------------------------------*
 * tbe_write_bitstream()
 *
 * Write TBE bitstream.
 *---------------------------------------------------------------------*/
void tbe_write_bitstream(
    Encoder_State *st
)
{
    short i;

    if ( (st->rf_mode || st->total_brate == ACELP_9k60) && st->bwidth == WB)
    {
        /* WB LSF */
        push_next_indice( st, st->lsf_WB, NUM_BITS_LBR_WB_LSF );

        /* WB frame */
        push_next_indice( st, st->gFrame_WB, NUM_BITS_SHB_FrameGain_LBR_WB );
    }
    else if (st->total_brate >= ACELP_9k60 && st->total_brate <= ACELP_32k && (st->bwidth == SWB || st->bwidth == FB))
    {
        /* LSF coefficients */

        if(st->rf_mode || st->total_brate == ACELP_9k60)
        {
            push_next_indice( st, st->lsf_idx[0], 8 );
        }
        else
        {
            for (i = 0; i < NUM_Q_LSF; i++)
            {
                push_next_indice( st, st->lsf_idx[i], lsf_q_num_bits[i] );
            }

            /* LSF mirror points */
            push_next_indice( st, st->m_idx, MIRROR_POINT_BITS );

            /* LSF grid points */
            push_next_indice( st, st->grid_idx, NUM_LSF_GRID_BITS );
        }

        /* Gain shape */
        push_next_indice( st, st->idxSubGains, NUM_BITS_SHB_SUBGAINS );

        /* frame gain */
        push_next_indice( st, st->idxFrameGain, NUM_BITS_SHB_FRAMEGAIN );

        if (st->total_brate >= ACELP_24k40)
        {
            /* sub frame energy*/
            push_next_indice( st, st->idx_shb_fr_gain, NUM_BITS_SHB_ENER_SF );

            /* gain shapes residual */
            for (i = 0; i < NB_SUBFR16k; i++)
            {
                push_next_indice( st, st->idx_res_gs[i], NUM_BITS_SHB_RES_GS );
            }

            /* voicing factor */
            push_next_indice( st, st->idx_mixFac, NUM_BITS_SHB_VF );
        }

        if( st->tec_tfa == 1 )
        {
            push_next_indice( st, st->tec_flag, BITS_TEC );
            push_next_indice( st, st->tfa_flag, BITS_TFA );
        }
    }

    if (st->bwidth == FB)
    {
        push_next_indice( st, st->idxGain, 4 );
    }
}



void TBEreset_enc(
    Encoder_State *st,                        /* i/o: encoder state structure                 */
    short  bandwidth                  /* i  : bandwidth mode                          */
)
{
    if(st->last_core != ACELP_CORE)
    {
        set_f( st->old_bwe_exc, 0.0f, PIT16k_MAX * 2 );
        st->bwe_non_lin_prev_scale = 0.f;
    }

    if( bandwidth == WB )
    {
        wb_tbe_extras_reset( st->mem_genSHBexc_filt_down_wb2, st->mem_genSHBexc_filt_down_wb3 );
        set_f( st->mem_genSHBexc_filt_down_shb, 0, 7 );
        set_f( st->state_lpc_syn, 0, 10 );
        set_f( st->state_syn_shbexc, 0, L_SHB_LAHEAD/4 );
        set_f( st->syn_overlap, 0, L_SHB_LAHEAD );
        set_f( st->mem_csfilt, 0, 2 );
    }
    else if( bandwidth == SWB || bandwidth == FB )
    {
        set_f( st->state_ana_filt_shb, 0.0f, (2*ALLPASSSECTIONS_STEEP+1) );

        swb_tbe_reset( st->mem_csfilt, st->mem_genSHBexc_filt_down_shb, st->state_lpc_syn,
                       st->syn_overlap, st->state_syn_shbexc, &(st->tbe_demph), &(st->tbe_premph),
                       st->mem_stp_swb, &(st->gain_prec_swb) );


        if( bandwidth == FB )
        {
            set_f(st->fb_state_lpc_syn, 0, LPC_SHB_ORDER);
            st->fb_tbe_demph = 0;
            fb_tbe_reset_enc( st->elliptic_bpf_2_48k_mem, &st->prev_fb_energy );
        }

    }

    return;
}

