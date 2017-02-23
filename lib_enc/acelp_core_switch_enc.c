/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_enc.h"
#include "rom_com.h"
#include "prot.h"
#include "prot.h"

/*---------------------------------------------------------------------*
 * Local functions
 *---------------------------------------------------------------------*/

static void encod_gen_voic_core_switch( Encoder_State *st, LPD_state *mem,const short L_frame, const float inp[],
                                        const float Aq[], const float A[], const short coder_type, const short T_op[],
                                        const float voicing[], float *exc, const long core_bitrate );

static void bwe_switch_enc( Encoder_State *st, const float *old_input );

/*-------------------------------------------------------------------*
 * acelp_core_switch_enc()
 *
 * ACELP core encoder in the ACELP->HQ switching frame
 *--------------------------------------------------------------------*/

void acelp_core_switch_enc(
    Encoder_State *st,                    /* i/o: encoder state structure             */
    LPD_state    *mem,                    /* i/o: encoder memories                    */
    const float inp12k8[],              /* i  : input signal @12.8 kHz              */
    const float inp16k[],               /* i  : input signal @16 kHz                */
    const short T_op_orig[2],           /* i  : open-loop pitch values for quantiz. */
    const float voicing[3],             /* i  : Open-loop pitch gains               */
    const float A[NB_SUBFR16k*(M+1)]    /* i  : A(z) unquantized for the 4 subframes*/
)
{
    short i, T_op[2];
    float old_exc[L_EXC], *exc;         /* excitation signal buffer                 */
    const float *inp;
    long  cbrate;
    float Aq[2*(M+1)];

    /* initializations */
    exc = old_exc + L_EXC_MEM;          /* pointer to excitation signal in the current frame */
    mvr2r( mem->old_exc, old_exc, L_EXC_MEM );

    mvr2r( st->old_Aq_12_8, Aq, M+1 );
    mvr2r( st->old_Aq_12_8, Aq + (M+1), M+1 );

    T_op[0] = T_op_orig[0];
    T_op[1] = T_op_orig[1];

    /*----------------------------------------------------------------*
     * set switching frame bit-rate
     *----------------------------------------------------------------*/

    if( st->last_L_frame == L_FRAME )   /* ACELP@12k8 core */
    {
        inp = inp12k8;

        if( st->core_brate > ACELP_24k40 )
        {
            cbrate = ACELP_24k40;
        }
        else
        {
            cbrate = st->core_brate;
        }
    }
    else /* ACELP@16k core */
    {
        inp = inp16k;

        if( st->core_brate <= ACELP_8k00 )
        {
            cbrate = ACELP_8k00;
        }
        else if( st->core_brate <= ACELP_14k80 )
        {
            cbrate = ACELP_14k80;
        }
        else
        {
            cbrate = min( st->core_brate, ACELP_22k60 );
        }
    }

    if( st->last_L_frame != st->L_frame )
    {
        if( st->last_L_frame == L_FRAME )
        {
            T_op[0] = (short)(0.8f * T_op[0] + 0.5f);
            T_op[1] = (short)(0.8f * T_op[1] + 0.5f);
        }
        else
        {
            T_op[0] = (short)(1.25f * T_op[0] + 0.5f);
            T_op[1] = (short)(1.25f * T_op[1] + 0.5f);
        }
    }

    /*----------------------------------------------------------------*
     * Excitation encoding
     *----------------------------------------------------------------*/

    encod_gen_voic_core_switch( st, mem, st->last_L_frame, inp, Aq, A, GENERIC, T_op, voicing, exc, cbrate );


    /*----------------------------------------------------------------*
     * Manipulate ACELP subframe indices (move them to their proper place)
     *----------------------------------------------------------------*/

    for( i=0; i<20; i++ )
    {
        st->ind_list[IND_CORE_SWITCHING_CELP_SUBFRAME+i].value = st->ind_list[TAG_ACELP_SUBFR_LOOP_START+i].value;
        st->ind_list[IND_CORE_SWITCHING_CELP_SUBFRAME+i].nb_bits = st->ind_list[TAG_ACELP_SUBFR_LOOP_START+i].nb_bits;
        st->ind_list[TAG_ACELP_SUBFR_LOOP_START+i].nb_bits = -1;
    }

    /*----------------------------------------------------------------*
     * BWE encoding
     *----------------------------------------------------------------*/

    if( !( (st->last_L_frame == L_FRAME16k && inner_frame_tbl[st->bwidth]==L_FRAME16k ) || inner_frame_tbl[st->bwidth] == L_FRAME8k ) )
    {
        bwe_switch_enc( st, (const float *)st->old_input_signal );
    }

    return;
}


/*-------------------------------------------------------------------*
 * encod_gen_voic_core_switch()
 *
 * Encode excitation signal in ACELP->HQ switching frame
 *-------------------------------------------------------------------*/

static void encod_gen_voic_core_switch(
    Encoder_State *st,               /* i/o: state structure                  */
    LPD_state *mem,              /* i/o: encoder memories                 */
    const short L_frame,           /* i  : length of the frame              */
    const float inp[],             /* i  : input signal                     */
    const float Aq[],              /* i  : LP coefficients                  */
    const float A[],               /* i  : unquantized A(z) filter          */
    const short coder_type,        /* i  : coding type                      */
    const short T_op[],            /* i  : open loop pitch                  */
    const float voicing[],         /* i  : voicing                          */
    float *exc,              /* i/o: current non-enhanced excitation  */
    const long  core_bitrate       /* i  : switching frame bit-rate         */
)
{
    float res[L_SUBFR];            /* residual signal                       */
    float Ap[M+1];                 /* A(z) with spectral expansion          */
    float xn[L_SUBFR];             /* Target vector for pitch search        */
    float xn2[L_SUBFR];            /* Target vector for codebook search     */
    float cn[L_SUBFR];             /* Target vector in residual domain      */
    float h1[L_SUBFR+(M+1)];       /* Impulse response vector               */
    float code[L_SUBFR];           /* Fixed codebook excitation             */
    float y1[L_SUBFR];             /* Filtered adaptive excitation          */
    float y2[L_SUBFR];             /* Filtered algebraic excitation         */
    float gain_pit ;               /* Pitch gain                            */
    float voice_fac;               /* Voicing factor                        */
    float gain_code;               /* Gain of code                          */
    float gain_inov;               /* inovation gain                        */
    short i;                       /* tmp variables                         */
    short T0, T0_frac;             /* close loop integer pitch and fractional part */
    short T0_min, T0_max;          /* pitch variables                       */
    float pitch;                   /* floating pitch value                  */
    float g_corr[6];               /* ACELP correl, values + gain pitch     */
    short clip_gain;               /* ISF clip gain                         */
    short unbits;                  /* number of unused bits for  EVS_PI     */
    float norm_gain_code;
    short pitch_limit_flag;
    float tmpF;
    short lp_select, lp_flag;

    /*------------------------------------------------------------------*
     * Initializations
     *------------------------------------------------------------------*/

    unbits = 0;

    if( L_frame == L_FRAME )
    {
        T0_max = PIT_MAX;
        T0_min = PIT_MIN;
    }
    else /* L_frame == L_FRAME16k */
    {
        T0_max = PIT16k_MAX;
        T0_min = PIT16k_MIN;
    }

    /*------------------------------------------------------------------*
     * Calculation of LP residual (filtering through A[z] filter)
     *------------------------------------------------------------------*/

    residu( Aq, M, inp, res, L_SUBFR );

    /*------------------------------------------------------------------*
     * ACELP subframe loop
     *------------------------------------------------------------------*/

    mvr2r( res, exc, L_SUBFR );

    if( L_frame == L_FRAME16k )
    {
        weight_a( A, Ap, GAMMA16k, M ); /* Bandwidth expansion of A(z) filter coefficients */
        find_targets( inp, mem->mem_syn, 0, &mem->mem_w0, Aq, res, L_SUBFR, Ap, PREEMPH_FAC_16k, xn, cn, h1 );
    }
    else
    {
        weight_a( A, Ap, GAMMA1, M );   /* Bandwidth expansion of A(z) filter coefficients */
        find_targets( inp, mem->mem_syn, 0, &mem->mem_w0, Aq, res, L_SUBFR, Ap, PREEMPH_FAC, xn, cn, h1 );
    }


    /*----------------------------------------------------------------*
     * Close-loop pitch search and quantization
     * Adaptive exc. construction
     *----------------------------------------------------------------*/

    pitch = pit_encode( st, core_bitrate, 0, L_frame, coder_type, &pitch_limit_flag,
                        0, exc, L_SUBFR, T_op, &T0_min, &T0_max, &T0, &T0_frac, h1, xn );

    /*-----------------------------------------------------------------*
     * Find adaptive exitation
     *-----------------------------------------------------------------*/

    pred_lt4( exc, exc, T0, T0_frac, L_SUBFR+1, inter4_2, L_INTERPOL2, PIT_UP_SAMP );

    /*-----------------------------------------------------------------*
     * Gain clipping test to avoid unstable synthesis on frame erasure
     *   or in case of floating point encoder & fixed p. decoder
     *-----------------------------------------------------------------*/

    clip_gain = gp_clip( core_bitrate, voicing, 0, coder_type, xn, st->clip_var );

    /*-----------------------------------------------------------------*
     * LP filtering of the adaptive excitation, codebook target computation
     *-----------------------------------------------------------------*/

    lp_select = lp_filt_exc_enc( MODE1, core_bitrate, 0, coder_type, 0, exc, h1,xn, y1, xn2, L_SUBFR,
                                 L_frame, g_corr, clip_gain, &gain_pit, &lp_flag );

    if( lp_flag == NORMAL_OPERATION )
    {
        push_indice( st, IND_LP_FILT_SELECT, lp_select, 1 );
    }

    /*-----------------------------------------------------------------*
     * Innovation encoding
     *-----------------------------------------------------------------*/

    inov_encode( st, core_bitrate, 0, L_frame, L_frame, coder_type, st->bwidth, 0, 0, -1, Aq, gain_pit, cn, exc,
                 h1, mem->tilt_code, pitch, xn2, code, y2, &unbits );

    /*-----------------------------------------------------------------*
     * Gain encoding
     *-----------------------------------------------------------------*/

    if (L_frame == L_FRAME)
    {
        gain_enc_mless( st, core_bitrate, L_frame, TRANSITION, 0, -1, xn, y1, y2, code, st->old_Es_pred,
                        &gain_pit, &gain_code, &gain_inov, &norm_gain_code, g_corr, clip_gain );
    }
    else
    {
        gain_enc_mless( st, core_bitrate, L_frame, coder_type, 0, -1, xn, y1, y2, code, st->old_Es_pred,
                        &gain_pit, &gain_code, &gain_inov, &norm_gain_code, g_corr, clip_gain );
    }

    gp_clip_test_gain_pit( core_bitrate, gain_pit, st->clip_var);
    mem->tilt_code = est_tilt( exc, gain_pit, code, gain_code, &voice_fac, L_SUBFR, 0 );

    /*-----------------------------------------------------------------*
     * Construct adaptive part of the excitation
     *-----------------------------------------------------------------*/

    for( i = 0; i < L_SUBFR; i++ )
    {
        tmpF = gain_pit * exc[i];
        exc[i] = tmpF + gain_code * code[i];
    }

    /* write reserved bits */
    if( unbits )
    {
        push_indice( st, IND_UNUSED, 0, unbits );
    }

    /*-----------------------------------------------------------------*
     * long term prediction on the 2nd sub frame
     *-----------------------------------------------------------------*/

    pred_lt4( &exc[L_SUBFR], &exc[L_SUBFR], T0, T0_frac, L_SUBFR+1, inter4_2, L_INTERPOL2,PIT_UP_SAMP );

    for( i = 0; i < L_SUBFR; i++ )
    {
        exc[i+L_SUBFR] *= gain_pit;
    }

    return;
}


/*-------------------------------------------------------------------*
 * bwe_switch_enc()
 *
 * Encode BWE in ACELP->HQ switching frame
 *-------------------------------------------------------------------*/

static void bwe_switch_enc(
    Encoder_State *st,                /* i/o: encoder state structure */
    const float *new_speech         /* i  : original input signal   */
)
{
    short i, k, delta, Loverlapp, d1, d1m, maxd1, ind1, fdelay, gapsize;
    float accA, accB, min_sq_cross, min_corr, E1, E2, gain;
    float tmp_mem[2*L_FILT48k], tmp_mem2[2*L_FILT48k], hb_synth_tmp[NS2SA(48000, 10000000L)];
    const float *hp_filter;
    float synth_subfr_bwe[SWITCH_MAX_GAP];              /* synthesized bwe for core switching */
    short n, L, input_frame;

    input_frame = st->input_Fs/50;

    L = NS2SA(st->input_Fs,FRAME_SIZE_NS);
    n = ((float)L * N_ZERO_MDCT_NS/FRAME_SIZE_NS);

    /* set multiplication factor according to the sampling rate */
    hp_filter = hp16000_48000;
    fdelay = 48;
    if( st->input_Fs == 16000 )
    {
        delta = 1;
        if( st->last_L_frame == L_FRAME )
        {
            hp_filter = hp12800_16000;
            fdelay = 20;
        }
    }
    else if( st->input_Fs == 32000 )
    {
        delta = 2;
        if( st->last_L_frame == L_FRAME )
        {
            hp_filter = hp12800_32000;
            fdelay = 40;
        }
        else
        {
            hp_filter = hp16000_32000;
            fdelay = 32;
        }
    }
    else
    {
        delta = 3;
        if( st->last_L_frame == L_FRAME )
        {
            hp_filter = hp12800_48000;
            fdelay = 60;
        }
    }

    set_f( tmp_mem, 0, 2*L_FILT48k );
    set_f( tmp_mem2, 0, 2*L_FILT48k );

    Loverlapp = delta*SWITCH_OVERLAP_8k*2;
    gapsize = delta * (NS2SA(16000,SWITCH_GAP_LENGTH_NS));
    set_f( synth_subfr_bwe, 0, SWITCH_MAX_GAP );

    for( i=0; i<gapsize+fdelay; i++ )
    {
        /* target */
        synth_subfr_bwe[i] = new_speech[i+L/2+n+Loverlapp-gapsize];
    }

    for( i=0; i<fdelay; i++ )
    {
        /* put the 40 past samples into the memory */
        tmp_mem[i] = new_speech[i+L/2+n+Loverlapp-gapsize-fdelay];
    }

    /* HP filtered target */
    fir( synth_subfr_bwe, hp_filter, synth_subfr_bwe, tmp_mem, gapsize+fdelay, fdelay, 0 );
    mvr2r( synth_subfr_bwe+(short)(fdelay/2), synth_subfr_bwe, delta*(NS2SA(16000,SWITCH_GAP_LENGTH_NS))-(short)(fdelay/2) );

    /* codebook */
    fir( new_speech, hp_filter, hb_synth_tmp, tmp_mem2, input_frame>>1, fdelay, 1 );


    min_sq_cross = -1;
    min_corr = 0;
    d1m = 0;

    maxd1 = (short)(((input_frame>>1) - gapsize-fdelay)/delta);

    /* find delay */
    for( k = 0, d1 = 0; k < maxd1; d1 += delta, k++ )
    {
        accA = accB = 0;
        for( i = 0; i < gapsize; i += delta )
        {
            accA += hb_synth_tmp[d1+i+fdelay] * hb_synth_tmp[d1+i+fdelay];
            accB += hb_synth_tmp[d1+i+fdelay] * synth_subfr_bwe[i];
        }
        if( accB * accB * min_corr >= min_sq_cross *accA )
        {
            d1m = k;
            min_corr = accA;
            min_sq_cross = accB * accB;
        }
    }

    push_indice( st, IND_CORE_SWITCHING_AUDIO_DELAY, d1m, AUDIODELAYBITS );

    /* find gain */
    E1 = 0.0f;
    E2 = 1.0f;   /* to avoid /0  */

    for( i=0; i<gapsize; i++ )
    {
        E1 += synth_subfr_bwe[i] * synth_subfr_bwe[i];
        E2 += hb_synth_tmp[i+d1m*delta+fdelay] * hb_synth_tmp[i+d1m*delta+fdelay];
    }

    gain = (float)sqrt((float)(E1/E2));

    ind1 = usquant( gain, &gain, MINVALUEOFFIRSTGAIN, DELTAOFFIRSTGAIN, (1 << NOOFGAINBITS1) );
    push_indice( st, IND_CORE_SWITCHING_AUDIO_GAIN,  ind1, NOOFGAINBITS1 );

    return;
}
