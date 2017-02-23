/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <stdlib.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"

/*---------------------------------------------------------------------*
 * Local functions
 *---------------------------------------------------------------------*/

static void hp400_12k8(float signal[], const short lg,  float mem[]);

static void filt_6k_7k(float signal[], const short lg, float mem[]);

static void hf_synthesis( const long core_brate, const short output_frame, const float Aq[], const float exc[],
                          float synth[], float synth16k[], short *seed2, float *mem_hp400, float *mem_syn_hf,
                          float *mem_hf, float *delay_syn_hf, float *mem_hp_interp );

static void hf_synthesis_amr_wb( const long  core_brate, const short output_subfr, const float Ap[], float exc16k[],
                                 float synth_out[], float *mem_syn_hf, float *delay_syn_hf, float *mem_hp_interp, float p_r,
                                 float hf_gain_i, float til, float voice_factors, const float exc[] );

static void envelope( const long  core_brate, const float Aq[], float Ap[], float *r, float tilt0, float tilt,
                      float voice_factor, float *prev_r, float *voice_fac, float *unvoicing, float *unvoicing_sm, short *unvoicing_flag);

static void AdaptiveStartBand( short *start_band, const long  rate, const float *lsf, const float voicing_fac, const short clas, short *voicing_flag,
                               short *start_band_old, float *OptCrit_old );


/*-------------------------------------------------------------------*
 * hf_synth_init()
 *
 * hf synthesis filters initialization
 * - initialization of 400 Hz high pass filter
 * - initialization of band pass 6kHz to 7kHz FIR filter
 *-------------------------------------------------------------------*/

void hf_synth_init(
    float mem_hp400[],  /* o  : 400 Hz high pass filter memory initialization    */
    float mem_hf[]      /* o  : band pass 6kHz to 7kHz FIR filter initialization */
)
{
    set_f( mem_hp400, 0, 4 );

    set_f( mem_hf, 0, (L_FIR-1) );

    return;
}

/*-------------------------------------------------------------------*
 * hf_synth_amr_wb_init()
 *
 * hf synthesis filters initialization
 * - initialization of 1600 Hz low pass filter
 * - initialization of band pass 6kHz to 8kHz FIR filter for noise and line resampled signals
 *-------------------------------------------------------------------*/

void hf_synth_amr_wb_init(
    float *prev_r,                    /* o  : previous sub-frame gain                          */
    float *fmerit_w_sm,               /* o  : 1 sample memory fmerit_w param                   */
    float mem_syn_hf[],               /* o  : HF LPC synthesis filter initialization           */
    short *frame_count,               /* o  : frame counter initialization                     */
    float *ne_min,                    /* o  : minimum Noise gate - short-term energy initialization*/
    float *fmerit_m_sm,               /* o  : 1 sample memory fmerit_m param                   */
    float *voice_fac,                 /* o  : voice factor initialization                      */
    float *unvoicing,                 /* o  : unvoiced parameter                               */
    float *unvoicing_sm,              /* o  : smoothed unvoiced parameter                      */
    short *unvoicing_flag,            /* o  : unvoiced flag                                    */
    short *voicing_flag,              /* o  : voiced flag                                      */
    short *start_band_old,            /* o  : previous start point for copying frequency band  */
    float *OptCrit_old                /* o  : previous criterion value for deciding the start point */
)
{
    *prev_r = 0.0f;
    set_f( mem_syn_hf, 0.0f, M );
    *fmerit_w_sm = 0.0f;
    *frame_count = 0;
    *ne_min = -30.0f;
    *fmerit_m_sm = 0.0f;
    *voice_fac = 0.0f;
    *unvoicing = 0.0f;
    *unvoicing_sm = 1.0f;
    *unvoicing_flag = 0;
    *voicing_flag = 0;
    *start_band_old = 160;
    *OptCrit_old = 1.0f;

    return;
}

/*-------------------------------------------------------------------*
 * hf_synth_amr_wb_reset()
 *
 * reset of HF synthesis filters
 * - needed in switching scenarios
 *-------------------------------------------------------------------*/

void hf_synth_amr_wb_reset(
    short *seed2,                     /* i/o: random seed for HF noise gen                      */
    float mem_syn_hf[],               /* o  : HF synthesis memory                               */
    float mem_hp_interp[],            /* o  : interpol. memory                                  */
    float *prev_r,                    /* o  : previous sub-frame gain                           */
    float *fmerit_w_sm,               /* o  : 1 sample memory fmerit_w param                    */
    float delay_syn_hf[],             /* o  : HF synthesis memory                               */
    short *frame_count,               /* o  : frame counter memory                              */
    float *ne_min,                    /* o  : minimum Noise gate - short-term energy memory     */
    float *fmerit_m_sm,               /* o  : 1 sample memory fmerit_m param                    */
    float *voice_fac,                 /* o  : voice factor memory                               */
    float *unvoicing,                 /* o  : unvoiced parameter                               */
    float *unvoicing_sm,              /* o  : smoothed unvoiced parameter                      */
    short *unvoicing_flag,            /* o  : unvoiced flag                                    */
    short *voicing_flag,              /* o  : voiced flag                                      */
    short *start_band_old,            /* o  : previous start point for copying frequency band  */
    float *OptCrit_old                /* o  : previous criterion value for deciding the start point */
)
{
    short i;

    for( i=0; i<L_FRAME16k; i++ )
    {
        own_random( seed2 );
    }

    set_f( mem_syn_hf, 0.0f, M );
    set_f( delay_syn_hf, 0.0, NS2SA(16000,DELAY_CLDFB_NS) );
    set_f( mem_hp_interp, 0, INTERP_3_1_MEM_LEN );
    *prev_r = 0.0f;
    *fmerit_w_sm = 0.0f;
    *frame_count = 0;
    *ne_min = -30.0f;
    *fmerit_m_sm = 0.0f;
    *voice_fac = 0.0f;
    *unvoicing = 0.0f;
    *unvoicing_sm = 1.0f;
    *unvoicing_flag = 0;
    *voicing_flag = 0;
    *start_band_old = 160;
    *OptCrit_old = 1.0f;

    return;
}

/*-------------------------------------------------------------------*
 * hf_synth_amr_wb()
 *
 * HF synthesis in AMR-WB IO
 *-------------------------------------------------------------------*/

void hf_synth_amr_wb(
    const long  core_brate,                /* i  : core bitrate                                      */
    const short output_frame,              /* i  : output frame length                               */
    const float *Aq,                       /* i  : quantized Az                                      */
    const float *exc,                      /* i  : excitation at 12.8 kHz                            */
    float *synth,                    /* i/o: synthesis signal at 12.8k                         */
    float *mem_syn_hf,               /* i/o: HF synthesis memory                               */
    float *delay_syn_hf,             /* i/o: HF synthesis memory                               */
    float *prev_r,                   /* i/o  : previous sub-frame gain                         */
    float *fmerit_w_sm,              /* i/o: smoothed fmerit_w                                 */
    short *amr_io_class,             /* i  : signal class (determined by FEC algorithm)        */
    float *mem_hp_interp,            /* i/o: interpol. memory                                  */
    float *synth_out,                /* i/o: output signal at output Fs                        */
    float fmerit,                    /* i  : classify parameter from FEC                       */
    const short *hf_gain,                  /* i  : decoded HF gain                                   */
    const float *voice_factors,            /* i  : voicing factors                                   */
    const float pitch_buf[],               /* i  : pitch buffer                                      */
    const float ng_ener_ST,                /* i  : Noise gate - short-term energy                    */
    const float *lsf_new,                  /* i  : ISF vector                                        */
    short *frame_count,              /* i/o: frame counter                                     */
    float *ne_min,                   /* i/o: minimum Noise gate                                */
    float *fmerit_m_sm,              /* i/o: smoothed fmerit_m                                 */
    float *voice_facor_sm,           /* o  : voice factor memory                               */
    float *unvoicing,                /* o  : unvoiced parameter                                */
    float *unvoicing_sm,             /* o  : smoothed unvoiced parameter                       */
    short *unvoicing_flag,           /* o  : unvoiced flag                                     */
    short *voicing_flag,             /* o  : voiced flag                                       */
    short *start_band_old,           /* o  : previous start point for copying frequency band   */
    float *OptCrit_old               /* o  : previous criterion value for deciding the start point */
)
{
    const float *p_Aq;
    float *p_Ap;
    short i, j, i_subfr, output_subfr;
    float Ap[(M16k+1)*NB_SUBFR];
    float exc16k[L_FRAME16k], dct_exc[L_FRAME], dct_hb[L_FRAME16k];
    float ener, tmp, scale;
    float alpha, beta, sub_gain[NB_SUBFR];
    short core_type = 1;
    float pitch_var_cur, voice_fac, fmerit_m, fmerit_w;
    short start_band;
    float til[NB_SUBFR], til0[NB_SUBFR];
    float enr1, enr2;
    float *pt1, *pt2;
    float e_subfr1[NB_SUBFR], e_subfr2[NB_SUBFR], e_subfr3[NB_SUBFR];
    float HF_corr_gain[NB_SUBFR];
    float filt_weight[80];
    short filt_weight_coeff;
    float gamma;
    float hb_ener, g, hb_tonal[80], tonal_ener, hb_amb[80], inv_g;
    short fb, fn, signum[80];

    pt1 = (float *)synth+1;
    pt2 = (float *)synth;
    for( i=0; i<NB_SUBFR; i++ )
    {
        enr1 = 0.0f;
        enr2 = 0.0f;
        enr1 = dotp( pt2, pt2, L_SUBFR );
        enr2 = dotp( pt1, pt2, (L_SUBFR-1) );
        til[i] = enr2/(enr1+0.1f);
        til0[i] = til[i];
        pt1 += L_SUBFR;
        pt2 += L_SUBFR;
    }

    output_subfr = output_frame/NB_SUBFR;

    if( *amr_io_class != 7 )
    {
        core_type = 0;
    }

    /* modify LF parameters for excitation weighting or sub-frame gains calculating */
    pitch_var_cur = (float)(fabs(pitch_buf[0] - pitch_buf[1]) + fabs(pitch_buf[1] - pitch_buf[2]) + fabs(pitch_buf[2] - pitch_buf[3]));

    if( *frame_count > FRAME_COUNT_HF_SYNTH && *amr_io_class == UNVOICED_CLAS )
    {
        *frame_count = 0;
        *ne_min = -30.0f;
    }
    else
    {
        if ( *frame_count > 2*FRAME_COUNT_HF_SYNTH )
        {
            *frame_count = 2*FRAME_COUNT_HF_SYNTH;
        }
        else
        {
            (*frame_count)++;
        }

        if ( ng_ener_ST < *ne_min )
        {
            *ne_min = ng_ener_ST;
        }
    }

    voice_fac = 0.0f;
    for( i=0; i<NB_SUBFR; i++ )
    {
        voice_fac += voice_factors[i];
    }
    voice_fac *= 0.25f;

    fmerit_w = fmerit > 0.35f ? 0.35f : (fmerit < 0.15f ? 0.15f : fmerit);
    if ( core_type == 1 )
    {
        fmerit_w *= 0.5f;
    }

    fmerit_w *= (1.0f + voice_fac);
    *fmerit_w_sm = 0.9f*(*fmerit_w_sm) + 0.1f*fmerit_w;
    fmerit_w = *fmerit_w_sm;

    fmerit_m = (2.0f - (fmerit < 0.5f ? 1.0f : fmerit));
    *fmerit_m_sm = 0.5f*(*fmerit_m_sm) + 0.5f*fmerit_m;
    fmerit_m = *fmerit_m_sm;

    for (i=0; i<NB_SUBFR; i++)
    {
        if(pitch_var_cur < 10 && (til[i]) < 0)
        {
            til[i] = 0.2f;
        }

        til[i] = (1.0f - (til[i])) < 0.8f ? 0.8f : (1.0f - (til[i]));
        til[i] += (30.0f + *ne_min) * 0.007f;
        til[i] *= fmerit_m;
    }

    /* predict LPC coefficients and calculate sub-frame gains */
    p_Aq = Aq;
    p_Ap = Ap;
    for( i=0; i<NB_SUBFR; i++ )
    {
        envelope( core_brate, p_Aq, p_Ap, &sub_gain[i], til0[i], til[i], voice_factors[i], prev_r,
                  voice_facor_sm, unvoicing, unvoicing_sm, unvoicing_flag );

        p_Aq += (M+1);
        p_Ap += (M+1);
    }

    /* rate dependent adaptive start band */
    AdaptiveStartBand( &start_band, core_brate, lsf_new, voice_fac, *amr_io_class, voicing_flag, start_band_old, OptCrit_old );

    /* DCT transform of LF excitation */
    edct( exc, dct_exc, L_FRAME );
    set_f( dct_hb, 0.0f, L_FRAME16k );

    /* copy excitation from LF */
    for ( i=200; i<240; i++ )
    {
        dct_hb[i] = dct_exc[i];
    }

    hb_ener = 0.01f;
    for ( i=240; i<L_FRAME16k; i++ )
    {
        dct_hb[i] = dct_exc[i+start_band-240];
        signum[i-240] = 1;
        if (dct_hb[i]<0)
        {
            signum[i-240] = -1;
        }
        dct_hb[i] *= signum[i-240];
        hb_ener += dct_hb[i]*dct_hb[i];
    }

    fmerit_w *= (1.1f - start_band*0.00625f);
    alpha = (float)sqrt(fmerit_w);
    beta = 1.0f - fmerit_w;
    gamma=max(0.3f, min(1.0f, (1.05f-alpha*0.95f)));
    tonal_ener = 0.01f;
    for (i=0; i<8; i++)
    {
        fb = 0;
        fn = i+8;
        tmp = 0;
        for (j=0; j<fn; j++)
        {
            tmp += dct_hb[j+240];
        }
        hb_amb[i] = tmp/fn;
        hb_tonal[i] = dct_hb[i+240] - hb_amb[i];
        if ( hb_tonal[i] > 0 )
        {
            tonal_ener += hb_tonal[i]*hb_tonal[i];
        }
    }
    for (; i<L_SUBFR16k-8; i++)
    {
        fb = i-7;
        tmp = 0;
        for (j=fb; j<fb+15; j++)
        {
            tmp += dct_hb[j+240];
        }
        hb_amb[i] = tmp/15;
        hb_tonal[i] = dct_hb[i+240] - hb_amb[i];
        if ( hb_tonal[i] > 0 )
        {
            tonal_ener += hb_tonal[i]*hb_tonal[i];
        }
    }
    for (; i<L_SUBFR16k; i++)
    {
        fb = i-7;
        fn = L_SUBFR16k-i+7;
        tmp = 0;
        for (j=fb; j<L_SUBFR16k; j++)
        {
            tmp += dct_hb[j+240];
        }
        hb_amb[i] = tmp/fn;
        hb_tonal[i] = dct_hb[i+240] - hb_amb[i];
        if ( hb_tonal[i] > 0 )
        {
            tonal_ener += hb_tonal[i]*hb_tonal[i];
        }
    }
    g = beta * (hb_ener-tonal_ener)/(hb_ener - beta*tonal_ener);
    if (g<0.01f && g>-0.01f)
    {
        inv_g = sign(g)*100;
    }
    else
    {
        inv_g = 1/g;
    }
    ener = 0.01f;
    for (i=0; i<L_SUBFR16k; i++)
    {
        if (hb_tonal[i]>0)
        {
            hb_tonal[i] *= g;
        }
        hb_amb[i] *= inv_g;
        dct_hb[i+240] = hb_tonal[i] + hb_amb[i];
        dct_hb[i+240] *= signum[i];
        ener += dct_hb[i+240]*dct_hb[i+240];
    }
    scale = (float)(gamma*sqrt(hb_ener/ener));

    if( core_brate == ACELP_6k60 )
    {
        filt_weight_coeff=60;
    }
    else if( core_brate == ACELP_8k85 )
    {
        filt_weight_coeff=40;
    }
    else
    {
        filt_weight_coeff=20;
    }

    for (i=0; i<filt_weight_coeff; i++)
    {
        filt_weight[i]=(float)(-0.999/(filt_weight_coeff-1))*i+1.0f;
    }

    for( i=240; i<L_FRAME16k; i++ )
    {
        dct_hb[i] *= scale;

        if( core_brate < ACELP_23k85 && i>255 )
        {
            dct_hb[i] *= 0.59525f;
        }

        if (i >= 320-filt_weight_coeff)
        {
            dct_hb[i] *= filt_weight[i-320+filt_weight_coeff];
        }
    }

    for( i=200; i<256; i++ )
    {
        dct_hb[i] *= filt_hp[i-200];

        if( core_brate < ACELP_23k85 )
        {
            dct_hb[i] *= deem_tab[i-200];
        }
    }

    if ( core_brate == ACELP_23k85 )
    {
        for (i=0; i<NB_SUBFR; i++)
        {
            HF_corr_gain[i] = 2*HP_gain[hf_gain[i]];
        }
    }

    /* inverse DCT transform of HF excitation */
    set_f( exc16k, 0.0f, L_FRAME16k );
    edct( dct_hb, exc16k, L_FRAME16k );

    /* energy weighting in consecutive subframes */
    ener = sum2_f( exc, L_FRAME ) + 0.01f;
    tmp = sum2_f( exc16k, L_FRAME16k ) + 0.01f;

    for (i=0; i<NB_SUBFR; i++)
    {
        e_subfr1[i] = sum2_f( &exc[i*L_SUBFR], L_SUBFR ) + 0.01f;
        e_subfr2[i] = sum2_f( &exc16k[i*L_SUBFR16k], L_SUBFR16k ) + 0.01f;
        e_subfr3[i] = (e_subfr1[i] / ener) * tmp;

        for (j=i*L_SUBFR16k; j<(i+1)*L_SUBFR16k; j++)
        {
            exc16k[j] *= (float)sqrt(e_subfr3[i]/e_subfr2[i]);
        }
    }

    p_Ap = Ap;
    i = 0;
    for( i_subfr=0; i_subfr<L_FRAME16k; i_subfr+=L_SUBFR16k )
    {
        /* synthesis of the HF signal */
        hf_synthesis_amr_wb( core_brate, output_subfr, p_Ap, &exc16k[i_subfr], &synth_out[i * output_subfr],
                             mem_syn_hf, delay_syn_hf, mem_hp_interp, sub_gain[i], HF_corr_gain[i], til0[i], voice_factors[i], &exc[i*L_SUBFR] );

        p_Ap += (M+1);
        i++;
    }

    return;
}

/*-----------------------------------------------------------------------------------*
 * hf_synthesis_amr_wb()
 *
 * HF noise synthesis
 * - Generate HF noise between 6 and 8 kHz, mix it with signal obtained by linear resampling.
 * - Set energy of high band
 *-----------------------------------------------------------------------------------*/

static void hf_synthesis_amr_wb(
    const long  core_brate,                 /* i  : core bitrate                    */
    const short output_subfr,               /* i  : output sub-frame length         */
    const float Ap[],                       /* i  : quantized Aq                    */
    float exc16k[],                   /* i  : excitation at 12.8 kHz          */
    float synth_out[],                /* i/o: synthesis signal at output Fs   */
    float *mem_syn_hf,                /* i/o: HF synthesis memory             */
    float *delay_syn_hf,              /* i/o: HF synthesis memory             */
    float *mem_hp_interp,             /* i/o: interpol. memory                */
    float p_r,                        /* i  : sub-frame gain                  */
    float HF_corr_gain,               /* i  : HF gain index                   */
    float til0,
    float voice_factors,
    const float exc[]
)
{
    short i, delay;
    float HF_syn[L_SUBFR16k], upsampled_HF_syn[L_FRAME48k/NB_SUBFR];
    float temp_buffer[NS2SA(16000,DELAY_CLDFB_NS)];
    float ener, tmp, scale, exc2385[L_SUBFR16k];

    if( core_brate == ACELP_23k85 )
    {
        ener = (sum2_f( exc, L_SUBFR ) + 0.01f)/5;
        tmp = sum2_f( exc16k, L_SUBFR16k )+0.01f;
        scale = (float)sqrt(ener/tmp);

        for (i=0; i<L_SUBFR16k; i++)
        {
            exc2385[i] = exc16k[i] * scale * HF_corr_gain;
        }
    }

    for ( i=0; i<L_SUBFR16k; i++ )
    {
        exc16k[i] *= p_r;
    }

    if ( core_brate == ACELP_23k85 )
    {
        ener = (sum2_f( exc16k, L_SUBFR16k ) + 0.01f)*0.3f;
        tmp = sum2_f( exc2385, L_SUBFR16k ) + 0.01f;
        scale = (float)sqrt(ener/tmp);

        if ( scale > 1.0f || til0 < 0.0f )
        {
            mvr2r( exc2385, exc16k, L_SUBFR16k );
        }
        else
        {
            for ( i=0; i<L_SUBFR16k; i++ )
            {
                exc16k[i] = exc2385[i] * max( min(1.0f,(1-til0)*(1.6f - voice_factors)), scale );
            }
        }
    }

    syn_filt( Ap, M, exc16k, HF_syn, L_SUBFR16k, mem_syn_hf, 1 );

    /*-----------------------------------------------------------------*
     * Resample to output sampling rate
     * Synchronize LB and HB components (delay componsation)
     * Add synthesised high band to speech synthesis
     *-----------------------------------------------------------------*/

    /* compensate CLDFB resampling delay */
    delay = NS2SA(16000,DELAY_CLDFB_NS);
    mvr2r( HF_syn+L_SUBFR16k-delay, temp_buffer, delay );
    mvr2r( HF_syn, HF_syn+delay, L_SUBFR16k-delay );
    mvr2r( delay_syn_hf, HF_syn, delay );
    mvr2r( temp_buffer, delay_syn_hf, delay );

    /* interpolate the HF synthesis */
    if( output_subfr == L_FRAME48k/NB_SUBFR )   /* 48kHz sampled output */
    {
        interpolate_3_over_1_allpass( HF_syn, L_SUBFR16k, upsampled_HF_syn, mem_hp_interp, allpass_poles_3_ov_2 );
    }
    else if( output_subfr == L_FRAME32k/NB_SUBFR )  /* 32kHz sampled output */
    {
        Interpolate_allpass_steep( HF_syn, mem_hp_interp, L_SUBFR16k, upsampled_HF_syn );
    }
    else    /* 16kHz sampled output */
    {
        mvr2r( HF_syn, upsampled_HF_syn, L_SUBFR16k );
    }

    v_add( synth_out, upsampled_HF_syn, synth_out, output_subfr );

    return;
}


/*-----------------------------------------------------------------------------------*
 * EnhanceClass()
 *
 *
 *-----------------------------------------------------------------------------------*/

static short EnhanceClass(
    const float qq,
    const float pp,
    const float tilt0,           /* i  : spectrum tilt                */
    const float tilt,            /* i  : spectrum tilt                */
    const float voice_factor,    /* i  : voice factor                 */
    float *voice_fac,      /* i/o: smoothed voiced parameter    */
    float *unvoicing,      /* i/o: unvoiced parameter           */
    float *unvoicing_sm,   /* i/o: smoothed unvoiced parameter  */
    short *unvoicing_flag  /* i/o: unvoiced flag                */
)
{
    float unvoicing_tmp;

    /* Decide (*unvoicing_flag) to allow BWE enhancement when qq>pp */
    *voice_fac = 0.75f*(*voice_fac) + 0.25f*voice_factor;
    unvoicing_tmp = ((1.0f-tilt0)/2.0f) * (1-*voice_fac) * min(tilt/1.5f, 1.0f);
    *unvoicing = 0.5f*(*unvoicing) + 0.5f*unvoicing_tmp;

    if( *unvoicing_sm > *unvoicing )
    {
        *unvoicing_sm = 0.9f*(*unvoicing_sm) + 0.1f*(*unvoicing);
    }
    else
    {
        *unvoicing_sm = 0.99f*(*unvoicing_sm) + 0.01f*(*unvoicing);
    }

    if (*unvoicing - *unvoicing_sm > 0.1f)
    {
        *unvoicing_flag = 1;
    }

    if (*unvoicing - *unvoicing_sm < 0.05f)
    {
        *unvoicing_flag = 0;
    }

    return (*unvoicing_flag && qq>pp);
}

/*-----------------------------------------------------------------------------------*
 * envelope()
 *
 *
 *-----------------------------------------------------------------------------------*/

static void envelope(
    const long  core_brate,      /* i  : core bitrate                    */
    const float Aq[],            /* i  : de-quant. LPC coefficents       */
    float Ap[],            /* o  : extended LPC coefficents        */
    float *sub_gain,       /* o  : sub-frame gain                  */
    float tilt0,           /* i  : spectrum tilt                   */
    float tilt,            /* i  : spectrum tilt                   */
    float voice_factor,    /* i  : voice factor                    */
    float *prev_r,         /* i/o: previous sub-frame gain         */
    float *voice_fac,      /* i/o: smoothed voiced parameter       */
    float *unvoicing,      /* i/o: unvoiced parameter              */
    float *unvoicing_sm,   /* i/o: smoothed unvoiced parameter     */
    short *unvoicing_flag  /* i/o: unvoiced flag                   */
)
{
    float px, py, rx, ry, pp, rr;
    short i, Unvoicing_flag;
    float alpha;
    float est_level1, est_level2, qx, qy, qq, env_level[3];
    float As[3], k1, k2;

    /* LPC envelope weighting */
    if( core_brate == ACELP_6k60 )
    {
        weight_a( Aq, Ap, 0.9f, M );
    }
    else
    {
        weight_a( Aq, Ap, 0.6f, M );
    }

    /* LPC envelope level estimate */
    pp = 0.0f;
    px = 0.0f;
    py = 0.0f;
    rr = 0.0f;
    rx = 0.0f;
    ry = 0.0f;
    for ( i=0; i<17; i++ )
    {
        px += Ap[i]*exp_tab_p[i];
        py += Ap[i]*exp_tab_p[33-i];
        rx += Aq[i]*exp_tab_q[i];
        ry += Aq[i]*exp_tab_q[33-i];
    }

    pp = 1.0f/((float)sqrt(px*px+py*py));
    rr = 1.0f/((float)sqrt(rx*rx+ry*ry));

    for ( i=0; i<3; i++ )
    {
        As[i] = Aq[i];
    }

    if (As[2]==-1)
    {
        k2 = -0.6f;
        k1 = 0.99f;
        if (As[1]<0)
        {
            k1 = -k1;
        }
    }
    else
    {
        k1 = As[1]/(1+As[2]);
        k2 = As[2];
        if (k2>0.6f)
        {
            k2=0.6f;
        }
        if (k2<-0.6f)
        {
            k2=-0.6f;
        }
        if (k1>0.99f)
        {
            k1=0.99f;
        }
        if (k1<-0.99f)
        {
            k1=-0.99f;
        }
    }

    As[1]=(1+k2)*k1;
    As[2]=k2;

    qq = 0.0f;
    qx = 0.0f;
    qy = 0.0f;
    for ( i=0; i<3; i++ )
    {
        qx += As[i]*exp_tab_q[i];
        qy += As[i]*exp_tab_q[33-i];
    }

    qq = 1.0f/((float)sqrt(qx*qx+qy*qy));

    Unvoicing_flag = EnhanceClass( rr, pp, tilt0, tilt, voice_factor, voice_fac, unvoicing, unvoicing_sm, unvoicing_flag );

    alpha  = 0.0f;
    if ( Unvoicing_flag )
    {
        if (rr > (*prev_r))
        {
            rr = 0.5f*rr + 0.5f*(*prev_r);
        }

        *prev_r = rr;
        rr *= min(1.0f, tilt*(1.6f - voice_factor));
        qq *= max(1.0f, tilt*(1.6f - voice_factor));
        rr = min(rr,qq);
        rr= max(rr,pp);
    }
    else
    {
        if (rr < 1.0f && (*prev_r) < 1.0f)
        {
            alpha = (1-rr*rr);
        }

        rr = alpha*(*prev_r) + (1-alpha)*rr;
        *prev_r = rr;

        est_level1 = qq * min(1.0f, tilt*(1.6f - voice_factor));
        env_level[0] = pp;
        env_level[1] = qq;
        env_level[2] = rr;
        v_sort( env_level, 0, 2 );
        rr = env_level[0];
        est_level2 = rr * (1.0f + (float)fabs(tilt-1) * (1.6f - voice_factor));
        rr = min(est_level1, est_level2);
    }

    *sub_gain = min(5.0f, rr/pp);

    return;
}

/*---------------------------------------------------------------------*
 * AdaptiveStartBand()
 *
 * adaptively select the start band of bandwidth extension
 *---------------------------------------------------------------------*/

static void AdaptiveStartBand(
    short *start_band,    /* o  : start point of copied band                */
    const long  rate,           /* i  : core bitrate                              */
    const float *lsf,           /* i  : lsf frequency                             */
    const float voicing_fac,    /* i  : voicing factors                           */
    const short clas, 	        /* i  : signal class (determined by FEC algorithm)*/
    short *voicing_flag,
    short *start_band_old,
    float *OptCrit_old
)
{
    float lsf_diff[M], Crit, OptCrit = 1.0f, W;
    short i, pos, M2, voicing_flag_old;
    short tmp1, tmp2;

    /*voicing switching flag : to avoid switching start band frequently in VOICED or AUDIO area*/
    voicing_flag_old = *voicing_flag;
    if( voicing_fac > 0.4f || (voicing_fac > 0.3f && clas >= VOICED_CLAS) || clas == AUDIO_CLAS )
    {
        *voicing_flag = 1;
    }

    if( voicing_fac < 0.2f && clas < VOICED_CLAS )
    {
        *voicing_flag = 0;
    }

    /* rate adaptive start band */
    *start_band = 160;
    if( rate < ACELP_23k05 )
    {
        for(i=1; i<(M-1); i++)
        {
            lsf_diff[i] = lsf[i] - lsf[i-1];
        }

        W = SQR(1.0f*rate/ACELP_19k85) / 6000.0f;

        if (clas == AUDIO_CLAS)
        {
            W *= 0.75f;
        }

        pos = 2;
        M2 = M-2;
        if( *voicing_flag == 1 )
        {
            if( rate <= ACELP_8k85 )
            {
                M2 = M-8;
            }
            else if( rate <= ACELP_12k65 )
            {
                M2 = M-6;
            }
            else if( rate <= ACELP_15k85 )
            {
                M2 = M-4;
            }
        }

        for(i=2; i<M2; i++)
        {
            Crit = lsf_diff[i]*max(1.0f - lsf[i]*W, 0.001f);
            if( Crit <= OptCrit || i == 2 )
            {
                OptCrit = Crit;
                pos = i;
            }
        }
        /* *start_band = (short)( (0.5f*(lsf[pos]+lsf[pos-1])*40/1000.0f - 40) + 0.5f ); */
        /* emulate BASOP precision: */
        tmp1 = (short) ((0.02f * lsf[pos]  ) + 0.5f);
        tmp2 = (short) ((0.02f * lsf[pos-1]) + 0.5f);
        *start_band = tmp1 + tmp2 - 40;
        *start_band = min(max(*start_band, 40), 160);

        if ( voicing_flag_old != *voicing_flag || ( *voicing_flag == 0 && OptCrit < *OptCrit_old ) ||
                ( OptCrit < 0.7f * (*OptCrit_old) && *OptCrit_old > 64 ) )
        {
            *OptCrit_old = OptCrit;
            if ( abs((*start_band)-(*start_band_old))<20 && *voicing_flag==1 && voicing_flag_old==1 )
            {
                *start_band = *start_band_old;
            }
        }
        else
        {
            if (OptCrit<(*OptCrit_old) && (*voicing_flag)==1)
            {
                *OptCrit_old = OptCrit;
            }

            *start_band = *start_band_old;
        }

        if (clas == AUDIO_CLAS)
        {
            *start_band = min(*start_band, 120);
        }

        if( *start_band % 2 != 0 )
        {
            *start_band -= 1;
        }
    }

    *start_band_old = *start_band;

    return;
}


/*-------------------------------------------------------------------*
 * hf_synth_reset()
 *
 * Reset of HF synthesis filters (needed in switching scenarios)
 *-------------------------------------------------------------------*/

void hf_synth_reset(
    short *seed2,             /* i/o: random seed for HF noise gen    */
    float mem_hf[],           /* o  : HF band-pass filter memory      */
    float mem_syn_hf[],       /* o  : HF synthesis memory             */
    float mem_hp400[],        /* o  : memory of hp 400 Hz filter      */
    float mem_hp_interp[],    /* o  : interpol. memory                */
    float delay_syn_hf[]      /* o  : HF synthesis memory             */
)
{
    short i;

    for( i=0; i<L_FRAME16k; i++ )
    {
        own_random( seed2 );
    }

    set_f( mem_hf, 0.0f, (L_FIR-1) );
    set_f( mem_syn_hf, 0.0f, M );
    set_f( mem_hp400, 0.0f, 4 );
    set_f( delay_syn_hf, 0.0f, NS2SA(16000,DELAY_CLDFB_NS) );
    set_f( mem_hp_interp, 0, INTERP_3_1_MEM_LEN );

    return;
}


/*---------------------------------------------------------------------*
 * hf_synth()
 *
 * High frequency regeneration
 *---------------------------------------------------------------------*/

void hf_synth(
    const long  core_brate,    /* i  : core bitrate                   */
    const short output_frame,  /* i  : output frame length            */
    const float *Aq,           /* i  : quantized Az                   */
    const float *exc,          /* i  : excitation at 12.8 kHz         */
    float *synth,        /* i/o: 12.8kHz synthesis signal       */
    float *synth16k,     /* i/o: 16kHz synthesis signal         */
    short *seed2,        /* i/o: random seed for HF noise gen   */
    float *mem_hp400,    /* i/o: memory of hp 400 Hz filter     */
    float *mem_syn_hf,   /* i/o: HF synthesis memory            */
    float *mem_hf,       /* i/o: HF band-pass filter memory     */
    float *delay_syn_hf, /*i/o: HF synthesis memory             */
    float *mem_hp_interp /* i/o: interpol. memory               */
)
{
    const float *p_Aq;
    short i_subfr, output_subfr;

    output_subfr = output_frame/NB_SUBFR;

    p_Aq = Aq;
    for( i_subfr=0; i_subfr<L_FRAME; i_subfr+=L_SUBFR )
    {
        hf_synthesis( core_brate, output_subfr, p_Aq, &exc[i_subfr], &synth[i_subfr], &synth16k[i_subfr * output_subfr/L_SUBFR],
                      seed2, mem_hp400, mem_syn_hf, mem_hf, delay_syn_hf, mem_hp_interp );
        p_Aq  += (M+1);
    }

    return;
}


/*-----------------------------------------------------------------------------------*
 * hf_synthesis()
 *
 * HF noise synthesis
 * - Generate HF noise between 6 and 7 kHz.
 * - Set energy of noise according to synthesis tilt.
 *     tilt > 0.8 ==> - 14 dB (voiced)
 *     tilt   0.5 ==> - 6 dB  (voiced or noise)
 *     tilt < 0.0 ==>   0 dB  (noise)
 *-----------------------------------------------------------------------------------*/

static void hf_synthesis(
    const long  core_brate,                 /* i  : core bitrate                    */
    const short output_subfr,               /* i  : output sub-frame length         */
    const float Aq[],                       /* i  : quantized Aq                    */
    const float exc[],                      /* i  : excitation at 12.8 kHz          */
    float synth[],                    /* i/o: 12.8kHz synthesis signal        */
    float synth16k[],                 /* i/o: 16kHz synthesis signal          */
    short *seed2,                     /* i/o: random seed for HF noise gen    */
    float *mem_hp400,                 /* i/o: memory of hp 400 Hz filter      */
    float *mem_syn_hf,                /* i/o: HF synthesis memory             */
    float *mem_hf,                    /* i/o: HF band-pass filter memory      */
    float *delay_syn_hf,              /* i/o: HF synthesis memory             */
    float *mem_hp_interp              /* i/o: interpol. memory                */
)
{
    short i, delay;
    float tmp, ener, fac, scale;
    float HF_exc[L_SUBFR16k], HF_syn[L_SUBFR16k], upsampled_HF_syn[L_FRAME48k/NB_SUBFR];
    float temp_buffer[NS2SA(16000,DELAY_CLDFB_NS) - L_FILT16k];
    float Ap[M16k+1];

    /*-----------------------------------------------------------------*
     * generate white noise vector
     *-----------------------------------------------------------------*/

    for( i=0; i<L_SUBFR16k; i++ )
    {
        HF_exc[i] = (float)own_random(seed2);
    }

    /*-----------------------------------------------------------------*
     * calculate energy scaling factor so that white noise would have the
     * same energy as exc12k8
     * note: the scaling factor should be multiplied by L_SUBFR16k / L_SUBFR for
     * correctness but it will have only negligible impact so it has been left like this
     *-----------------------------------------------------------------*/

    ener = sum2_f( exc, L_SUBFR ) + 0.01f;
    tmp = sum2_f( HF_exc, L_SUBFR16k ) + 0.01f;
    scale = (float)(sqrt(ener/tmp));

    /*-----------------------------------------------------------------*
     * calculate energy scaling factor to respect tilt of synth12k8
     * (tilt: 1=voiced, -1=unvoiced)
     *-----------------------------------------------------------------*/

    hp400_12k8( synth, L_SUBFR, mem_hp400 );

    ener = 0.001f;
    tmp = 0.001f;

    for( i=1; i<L_SUBFR; i++ )
    {
        ener += synth[i] * synth[i];
        tmp += synth[i] * synth[i-1];
    }

    fac = tmp/ener;
    fac = (float)(1.0f - fac);

    if( core_brate == SID_2k40 || core_brate == FRAME_NO_DATA )
    {
        /* emphasize HF noise in CNG */
        fac *= 2.0f;
    }

    if( fac < 0.1f )
    {
        fac = 0.1f;
    }

    if( fac > 1.0f )
    {
        fac = 1.0f;
    }

    scale *= fac;
    /*-----------------------------------------------------------------*
     * modify HF excitation according to both calculated scaling factors
     *-----------------------------------------------------------------*/

    for( i=0; i<L_SUBFR16k; i++ )
    {
        HF_exc[i] *= scale;
    }

    /*-----------------------------------------------------------------*
     * high pass filtering (0.94ms of delay)
     *-----------------------------------------------------------------*/

    filt_6k_7k( HF_exc, L_SUBFR16k, mem_hf );

    /*-----------------------------------------------------------------*
     * synthesis of noise: 4.8kHz..5.6kHz --> 6kHz..7kHz
     *-----------------------------------------------------------------*/

    weight_a( Aq, Ap, 0.6f, M );
    syn_filt( Ap, M, HF_exc, HF_syn, L_SUBFR16k, mem_syn_hf, 1 );

    /*-----------------------------------------------------------------*
     * Add filtered HF noise to speech synthesis
     *-----------------------------------------------------------------*/

    /* delay by 5 samples @16kHz to compensate CLDFB resampling delay (20samples) and HP filtering delay (roughly 15 samples) */
    delay = NS2SA(16000,DELAY_CLDFB_NS) - 15;
    mvr2r( HF_syn+L_SUBFR16k-delay, temp_buffer, delay );
    mvr2r( HF_syn, HF_syn+delay, L_SUBFR16k-delay );
    mvr2r( delay_syn_hf, HF_syn, delay );
    mvr2r( temp_buffer, delay_syn_hf, delay );

    /* interpolate the HF synthesis */
    if( output_subfr == L_FRAME48k/NB_SUBFR )   /* 48kHz sampled output */
    {
        interpolate_3_over_1_allpass( HF_syn, L_SUBFR16k, upsampled_HF_syn, mem_hp_interp, allpass_poles_3_ov_2 );
    }
    else if( output_subfr == L_FRAME32k/NB_SUBFR )  /* 32kHz sampled output */
    {
        Interpolate_allpass_steep( HF_syn, mem_hp_interp, L_SUBFR16k, upsampled_HF_syn );
    }
    else    /* 16kHz sampled output */
    {
        mvr2r( HF_syn, upsampled_HF_syn, L_SUBFR16k );
    }

    v_add( synth16k, upsampled_HF_syn, synth16k, output_subfr );

    return;
}


/*-----------------------------------------------------------------------*
 * hp400_12k8()
 *
 * 2nd order Cheb2 high pass filter with cut off frequency at 400 Hz.
 * Optimized for fixed-point to get the following frequency response:
 *
 *  frequency  :   0Hz   100Hz  200Hz  300Hz  400Hz  630Hz  1.5kHz  3kHz
 *  dB loss  :   -infdB  -30dB  -20dB  -10dB  -3dB   +6dB    +1dB    0dB
 *
 * Algorithm  :
 *
 *  y[i] = b[0]*x[i] + b[1]*x[i-1] + b[2]*x[i-2]
 *                   + a[1]*y[i-1] + a[2]*y[i-2];
 *
 *  short b[3] = {3660, -7320,  3660};       in Q12
 *  short a[3] = {4096,  7320, -3540};       in Q12
 *
 *  float b[3] = {0.893554687, -1.787109375,  0.893554687};
 *  float a[3] = {1.000000000,  1.787109375, -0.864257812};
 *-----------------------------------------------------------------------*/

static void hp400_12k8(
    float signal[],  /* i/o: signal            */
    const short lg,        /* i  : lenght of signal  */
    float mem[]      /* i/o: filter memory [4] */
)
{
    short i;
    float x0, x1, x2;
    float yy0, yy1, y2;


    yy1 = mem[0];
    y2 = mem[1];
    x0 = mem[2];
    x1 = mem[3];

    for( i=0; i<lg; i++ )
    {
        x2 = x1;
        x1 = x0;
        x0 = signal[i];
        yy0 = yy1*a_hp400[1] + y2*a_hp400[2] + x0*b_hp400[0] + x1*b_hp400[1] + x2*b_hp400[2];

        signal[i] = yy0;
        y2 = yy1;
        yy1 = yy0;
    }

    mem[0] = yy1;
    mem[1] = y2;
    mem[2] = x0;
    mem[3] = x1;

    return;
}

/*-------------------------------------------------------------------*
 * filt_6k_7k()
 *
 * 15th order band pass 6kHz to 7kHz FIR filter
 *
 * frequency  :4kHz   5kHz  5.5kHz  6kHz  6.5kHz 7kHz  7.5kHz  8kHz
 * dB loss  : -60dB  -45dB  -13dB   -3dB   0dB   -3dB  -13dB  -45dB
 * (gain = 4.0)
 *-------------------------------------------------------------------*/

static void filt_6k_7k(
    float signal[],  /* i/o: signal        */
    const short lg,        /* i  : signal length */
    float mem[]      /* i/o: filter memory */
)
{
    short i, j;
    float s, x[L_FRAME48k/NB_SUBFR+(L_FIR-1)];

    for( i=0; i<(L_FIR-1); i++ )
    {
        x[i] = mem[i];
    }

    for(i=0; i<lg; i++)
    {
        x[i+(L_FIR-1)] = signal[i];
    }

    for(i=0; i<lg; i++)
    {
        s = 0.0;
        for(j=0; j<L_FIR; j++)
        {
            s += x[i+j] * fir_6k_7k[j];
        }

        /* gain of coef = 4.0 */
        signal[i] = (float)(s*0.25);
    }

    for( i=0; i<L_FIR-1; i++ )
    {
        mem[i] = x[i+lg];
    }

    return;
}
