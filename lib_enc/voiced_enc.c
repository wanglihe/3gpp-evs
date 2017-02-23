/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <stdlib.h>
#include "prot.h"
#include "cnst.h"
#include "rom_com.h"


/*-------------------------------------------------------------------*
 * Local functions
 *--------------------------------------------------------------------*/

static void SynthesisFilter( float *output, float *input, float *coef, float *memory, short order, short length );

/*-------------------------------------------------------------------*
 * ppp_voiced_encoder()
 *
 * voiced encoder for SC-VBR
 *--------------------------------------------------------------------*/

void ppp_voiced_encoder(
    Encoder_State *st,    /* i/o: state structure */
    float *in,    /* i  : residual signal */
    float *out,   /* o  : Quantized residual signal */
    short delay,  /* i  : open loop pitch */
    float *lpc1,  /* i  : prev frame de-emphasized LPC */
    float* lpc2,  /* i  : current frame de-emphasized LPC */
    float *exc,   /* i  : previous frame quantized excitation */
    float *pitch, /* o  : floating pitch values for each subframe */
    float vadsnr  /* i  : current frame SNR */
)
{
    short i;
    short spike_near_edge = 0;
    short flag;
    short delta_lag_E = 0, PPP_MODE_E, Q_delta_lag = 0;
    short out_of_bound = 0;
    float tmp, tmptmp, tmptmp1, res_enratio = 0, sp_enratio = 0;
    double tmp2;
    int pl, l;
    float interp_delay[3], temp_pl, temp_l;
    short upper_cut_off_freq_of_interest = 0, upper_cut_off_freq = 0;
    float pos_nq, neg_nq, pos_q, neg_q;
    float impzi[160];
    float impzo[160];
    float mem[10];
    float energy_impz = 0.0f;
    float pos_nq0, neg_nq0, tmpres;
    float sp_hb_enratio;
    float low_band_en;

    DTFS_STRUCTURE *CURRP_NQ = DTFS_new();
    DTFS_STRUCTURE *TMPDTFS = DTFS_new();
    DTFS_STRUCTURE *TMPDTFS2 = DTFS_new();
    DTFS_STRUCTURE *TMPDTFS3 = DTFS_new();
    DTFS_STRUCTURE *CURRP_Q_E = DTFS_new();
    DTFS_STRUCTURE *dtfs_temp = DTFS_new();

    if ( st->bwidth == WB )
    {
        upper_cut_off_freq_of_interest = 4000;
        upper_cut_off_freq = 6400;
    }
    else if ( st->bwidth == NB )
    {
        upper_cut_off_freq_of_interest = 3300;
        upper_cut_off_freq = 4000;
    }

    /* Initialization */
    if (st->firstTime_voicedenc)
    {
        st->firstTime_voicedenc = 0;
        st->dtfs_enc_lag = 0;
        st->dtfs_enc_nH = 0;
        st->dtfs_enc_nH_4kHz = 0;
        st->dtfs_enc_upper_cut_off_freq_of_interest = 3300.0;
        st->dtfs_enc_upper_cut_off_freq = 4000.0;

        set_f(st->dtfs_enc_a, 0, MAXLAG_WI);
        set_f(st->dtfs_enc_b, 0, MAXLAG_WI);
    }

    /* Figure out the PPP_MODE */
    if ( st->last_ppp_mode == 1 && !st->mode_QQF )
    {
        st->bump_up = 1;

        free(CURRP_NQ);
        free(TMPDTFS);
        free(TMPDTFS2);
        free(TMPDTFS3);
        free(CURRP_Q_E);
        free(dtfs_temp);

        return;
    }

    /* Use the aggresive bumpups if there are two consecutive Q frames */
    /* Aggresive bump upsare only used in the second Q frame */
    if ( st->last_ppp_mode == 1 )
    {
        st->rate_control = 0;
    }

    PPP_MODE_E = 'Q';
    pl = min(MAX_LAG_PIT, (int)rint_new(st->old_pitch_buf[(2*NB_SUBFR)-1]));
    l = min(MAX_LAG_PIT, (int)rint_new(delay));

    /* Bump up if the lag is out_fx of range */
    if ((l-pl)>13 || (l-pl)<-11  || (l < 19) || (pl < 19) )
    {
        st->bump_up = 1;

        free(CURRP_NQ);
        free(TMPDTFS);
        free(TMPDTFS2);
        free(TMPDTFS3);
        free(CURRP_Q_E);
        free(dtfs_temp);
        return;
    }

    if (st->last_ppp_mode!=1)
    {
        /* Obtain DTFS of last pl values of past excitation */
        DTFS_to_fs(exc-pl, pl, dtfs_temp, st->bwidth == WB ? (short)16000 : (short)8000, 0 );

    }

    if( st->last_coder_type_raw == UNVOICED )
    {
        pl = l; /* if prev frame was sil/uv */
    }

    /* Use the out array as a temp storage for currp */
    spike_near_edge = ppp_extract_pitch_period(in, out, l, &out_of_bound) ;

    if (out_of_bound == 1)
    {
        st->bump_up = 1;

        free(CURRP_NQ);
        free(TMPDTFS);
        free(TMPDTFS2);
        free(TMPDTFS3);
        free(CURRP_Q_E);
        free(dtfs_temp);

        return;
    }

    /* Get DTFS of current prototype */
    DTFS_to_fs(out, l, CURRP_NQ, st->bwidth == WB ? (short)16000 : (short)8000, 0);

    /* Ensure the extracted prototype is time-synchronous to the
     * last l samples of the frame. This proves to eliminate
     * some of the PPP-CELP transition problems.
     * Convert last samples into DTFS  */
    if (spike_near_edge != 0)
    {
        DTFS_to_fs(in+L_FRAME-l, l, TMPDTFS, st->bwidth == WB ? (short)16000 : (short)8000 ,0);

        tmp = DTFS_alignment_extract(*TMPDTFS, *CURRP_NQ, 0.0, lpc2) ; /* figure out how much to shift currp_nq to align with TMP */

        DTFS_phaseShift(CURRP_NQ, (float)(PI2*tmp/l)) ;
    }

    temp_pl = (float) pl;
    temp_l = (float) l;

    for(i = 0; i < NB_SUBFR; i++)
    {
        /* do the linear pitch interp to drive the nb_post_filt */
        Interpol_delay(interp_delay, &(temp_pl), &(temp_l), i, frac_4sf);
        pitch[i] = interp_delay[0];
    }



    /* Restoring PPP memories when the last frame is non-PPP */
    if ( st->last_ppp_mode != 1 )
    {
        st->ph_offset_E = 0.0 ;
        st->prev_cw_en = DTFS_getEngy(*dtfs_temp);

        /* Copy over dtfs_temp into TMPDTFS */
        DTFS_copy(TMPDTFS, *dtfs_temp);

        DTFS_car2pol(TMPDTFS);

        st->lastLgainE = (float) log10(TMPDTFS->lag*DTFS_setEngyHarm(92.0, 1104.5, 0.0, 1104.5, 1.0, TMPDTFS));
        st->lastHgainE = (float) log10(TMPDTFS->lag*DTFS_setEngyHarm(1104.5, upper_cut_off_freq_of_interest, 1104.5, upper_cut_off_freq, 1.0, TMPDTFS));

        DTFS_to_erb(*TMPDTFS, st->lasterbE);
    }
    else
    {
        /* Copy DTFS related parameters from 'st' to 'dtfs_temp' structure */
        dtfs_temp->lag = st->dtfs_enc_lag;
        dtfs_temp->nH = st->dtfs_enc_nH;
        dtfs_temp->nH_4kHz = st->dtfs_enc_nH_4kHz;
        dtfs_temp->upper_cut_off_freq_of_interest = st->dtfs_enc_upper_cut_off_freq_of_interest;
        dtfs_temp->upper_cut_off_freq = st->dtfs_enc_upper_cut_off_freq;

        mvr2r(st->dtfs_enc_a, dtfs_temp->a, MAXLAG_WI);
        mvr2r(st->dtfs_enc_b, dtfs_temp->b, MAXLAG_WI);
    }

    /*-----Open-loop Bump-Up-------- */

    /* Energy ratio calculation in residual and speech domain */
    /* Also, compute correlation between the previous and the */
    /* current prototype */
    res_enratio = DTFS_getEngy(*CURRP_NQ) / DTFS_getEngy(*dtfs_temp);

    /* Copy over CURRP_NQ into TMPDTFS */
    DTFS_copy(TMPDTFS, *CURRP_NQ);

    /* Copy over dtfs_temp into TMPDTFS2 */
    DTFS_copy(TMPDTFS2, *dtfs_temp);

    tmptmp = DTFS_alignment_full(*TMPDTFS2, *TMPDTFS, TMPDTFS->lag*2); /* align of prev_cw wrt curr_cw, new method */

    tmptmp1 = TMPDTFS->lag-tmptmp;
    tmp = tmptmp1;

    DTFS_phaseShift(TMPDTFS, (float)(-PI2*tmp/TMPDTFS->lag)); /* fixed bug , phase shift by tmp computed in TMP.lag domain (above) */
    tmpres = (float)(DTFS_freq_corr(*TMPDTFS, *TMPDTFS2, 100.0f, 3700.0f));

    DTFS_poleFilter(TMPDTFS, lpc2, M+1);

    DTFS_adjustLag(TMPDTFS2, TMPDTFS->lag); /* operate in CL domain */

    DTFS_poleFilter(TMPDTFS2, lpc1, M+1);

    tmp = (float)(DTFS_freq_corr(*TMPDTFS, *TMPDTFS2, 100.0f, 3700.0f));

    if ( DTFS_getEngy(*TMPDTFS2) > 0 )
    {
        sp_enratio = DTFS_getEngy(*TMPDTFS)/DTFS_getEngy(*TMPDTFS2);
    }
    else
    {
        sp_enratio = 0.0f;
    }

    if (PPP_MODE_E == 'Q')
    {
        /* Bump up if the lag is out of range */
        if (((l-pl)>13) ||((l-pl)<-11))
        {
            PPP_MODE_E = 'B' ;
        }
        else
        {
            delta_lag_E=(short) (l-pl);
        }

        /* Bump up if big change between the previous and the current CWs */
        if ( vadsnr < st->SNR_THLD )
        {
            if ( res_enratio > 5.0 && tmp < 0.65 )
            {
                PPP_MODE_E = 'B';
            }
        }
        else
        {
            if ( res_enratio > 3.0 && tmp < 1.2 )
            {
                PPP_MODE_E = 'B';
            }
        }
    }

    /* Rapid rampdown frame where time resolution is important */
    /* Not a suitable PPP frame -> Bump to CELP */
    if ( vadsnr < st->SNR_THLD )
    {
        if (res_enratio < 0.025)
        {
            st->bump_up = 1;

            free(CURRP_NQ);
            free(TMPDTFS);
            free(TMPDTFS2);
            free(TMPDTFS3);
            free(CURRP_Q_E);
            free(dtfs_temp);

            return;
        }
    }
    else
    {
        if ( res_enratio < 0.092f )
        {
            st->bump_up = 1;
        }
    }

    if ( min(res_enratio, sp_enratio) < 0.075f && tmp < -0.5f )
    {
        st->bump_up = 1;
    }

    /* Rapid rampup frame where time resolution is important */
    /* Not a suitable PPP frame -> Bump to CELP */
    if ( vadsnr < st->SNR_THLD )
    {
        if (res_enratio > 14.5)
        {
            st->bump_up = 1;

            free(CURRP_NQ);
            free(TMPDTFS);
            free(TMPDTFS2);
            free(TMPDTFS3);
            free(CURRP_Q_E);
            free(dtfs_temp);

            return;
        }
    }
    else
    {
        if (res_enratio > 7.0)
        {
            st->bump_up = 1;
        }
    }

    if ( st->bump_up == 1 )
    {

        free(CURRP_NQ);
        free(TMPDTFS);
        free(TMPDTFS2);
        free(TMPDTFS3);
        free(CURRP_Q_E);
        free(dtfs_temp);

        return;
    }

    /* Bump up when the previous frame is an unvoiced or a silent frame */
    if( st->last_coder_type_raw == UNVOICED )
    {
        st->bump_up=1;

        free(CURRP_NQ);
        free(TMPDTFS);
        free(TMPDTFS2);
        free(TMPDTFS3);
        free(CURRP_Q_E);
        free(dtfs_temp);

        return;
    }
    /* -----End Open-loop Bump-Up */

    /* PPP-WI Quantization */
    if (PPP_MODE_E == 'Q')
    {
        flag = 1;
        if (PPP_MODE_E == 'Q')
        {
            flag = ppp_quarter_encoder( st, CURRP_Q_E, TMPDTFS, dtfs_temp->lag, *CURRP_NQ, lpc2, &(st->lastLgainE), &(st->lastHgainE), &(st->lasterbE[0]), *dtfs_temp );
        }

        if (flag)
        {
            /* TMPDTFS : Target prototype: Amp Quantized + Phase Unquantized         */
            /* TMPDTFS2: Quantized prototype: Amp Quantized + Phase Quantized        */
            /* TMPDTFS3: Delta prototype: Diff betw. target and quant. in speech dom */

            /* ----- Closed-loop Bump-Up ---------- */
            DTFS_peaktoaverage(*TMPDTFS, &pos_nq, &neg_nq);
            DTFS_peaktoaverage(*CURRP_Q_E, &pos_q, &neg_q);

            /* Before we perform the peak-to-average ratio comparison, we have to */
            /* ensure that the energy is not decaying and also the pitch pulse */
            /* is clearly defined */

            /*  Usually triggers in the slow ramp down frames.	Does not fall under the test condition (res_enratio < 0.025) as
                both frames have little energy and the ratio is not very small. Not	suitable for PPP */

            if ( CURRP_Q_E->upper_cut_off_freq > 4000 )
            {
                /* Use this bumup only for WB signals */
                if ( DTFS_getEngy_band_wb(*CURRP_Q_E, 0.0, 2000.0) > 0 )
                {
                    sp_hb_enratio = DTFS_getEngy_band_wb(*CURRP_Q_E, 2000.0, 6400.0)/DTFS_getEngy_band_wb(*CURRP_Q_E, 0.0, 2000.0);
                }
                else
                {
                    sp_hb_enratio = 0;
                }

                low_band_en = (float)DTFS_getEngy_band_wb(*CURRP_Q_E, 0.0, 2000.0);

                if ( low_band_en < 25.0f && sp_hb_enratio < 1.6f )
                {
                    PPP_MODE_E = 'B';
                }
            }

            if ( vadsnr < st->SNR_THLD )
            {
                if ( DTFS_getEngy(*CURRP_NQ) > 0.8f * st->prev_cw_en && max(pos_nq, neg_nq) > 3.0f && st->rate_control )
                {
                    if ( pos_nq > neg_nq && pos_nq > 2.0f * pos_q )
                    {
                        PPP_MODE_E = 'B';
                    }

                    if ( pos_nq < neg_nq && neg_nq > 2.0f * neg_q )
                    {
                        PPP_MODE_E = 'B';
                    }
                }
            }
            else
            {
                if ((((DTFS_getEngy(*CURRP_NQ) >(st->prev_cw_en))&&(max(pos_nq,neg_nq)>3.5))&&(st->rate_control))||
                        (((DTFS_getEngy(*CURRP_NQ) >0.8*(st->prev_cw_en))&&(max(pos_nq,neg_nq)>3.0))&&(!st->rate_control)))
                {
                    if (((pos_nq > neg_nq) && (pos_nq > 2.5*pos_q)&&(st->rate_control))||
                            ((pos_nq > neg_nq) && (pos_nq > 2.0*pos_q)&&(!st->rate_control)))
                    {
                        PPP_MODE_E='B';
                    }

                    if ((((pos_nq < neg_nq) && (neg_nq > 2.5*neg_q))&&(st->rate_control))||
                            ((pos_nq < neg_nq) && (neg_nq > 2.0*neg_q)&&(!st->rate_control)))
                    {
                        PPP_MODE_E='B';
                    }
                }

                if ( st->rate_control )
                {
                    DTFS_peaktoaverage(*CURRP_NQ, &pos_nq0,&neg_nq0);

                    for (impzi[0]=1.0,i=1; i<160; i++)
                    {
                        impzi[i]=0.0;
                    }

                    for (i=0; i<160; i++)
                    {
                        impzo[i]=0.0;
                    }

                    for (i=0; i<10; i++)
                    {
                        mem[i]=0.0;
                    }

                    SynthesisFilter(&impzo[0],&impzi[0],lpc2,&mem[0],10,160);

                    for (i=0; i<160; i++)
                    {
                        energy_impz+=(impzo[i]*impzo[i]);
                    }

                    energy_impz = (float)(10*log10((float)energy_impz));

                    if ((DTFS_getEngy(*CURRP_Q_E) > st->prev_cw_en)&&(max(pos_q,neg_q)>3.5)&&energy_impz>15.0&&tmpres>0.7)
                    {
                        if ((pos_q > neg_q) && ((pos_q>3.0*pos_nq0) || ((pos_q > 1.5*pos_nq0) && (neg_q < 1.5*neg_q))))
                        {
                            PPP_MODE_E='B';
                        }

                        if ((pos_q <= neg_q) && ((neg_q>3.0*neg_nq0)|| ((neg_q > 1.5*neg_nq0) && (pos_q < 1.5*pos_q))))
                        {
                            PPP_MODE_E='B';
                        }
                    }
                }
            }

            DTFS_copy(TMPDTFS2,*CURRP_Q_E);

            DTFS_poleFilter(TMPDTFS,lpc2,M+1);
            DTFS_poleFilter(TMPDTFS2,lpc2,M+1);

            *TMPDTFS3 = DTFS_sub(*TMPDTFS,*TMPDTFS2);

            /* operate in ADR mode only the rate control is active. This adds some bumpups to improve the speech quality */
            tmp2 = DTFS_getEngy_band(*TMPDTFS, 1500.0, upper_cut_off_freq_of_interest)/DTFS_getEngy(*TMPDTFS);
            if (tmp2==0.0)
            {
                tmp2=0.001;
            }
            if ((tmp2> 0.05)&&(!st->rate_control))
            {
                tmp2=DTFS_getEngy_band(*TMPDTFS,1500.0,upper_cut_off_freq_of_interest)/DTFS_getEngy_band(*TMPDTFS3,1500.0,upper_cut_off_freq_of_interest);
                if (tmp2==0.0)
                {
                    tmp2=0.001;
                }
                if (10.0*log10(tmp2) < 0.1)
                {
                    if (res_enratio > 0.8)
                    {
                        PPP_MODE_E = 'B';
                    }
                }
            }

            /* To increase bump up, raise first threshold, lower second  */
            tmp = (float)(10.0*log10(DTFS_getEngy(*TMPDTFS)/DTFS_getEngy(*TMPDTFS3)));

            if ((tmp <= 0)&&(!st->rate_control))
            {
                PPP_MODE_E = 'B';
            }

            if ( vadsnr < st->SNR_THLD )
            {
                if ((( tmp < 3.05 && max(res_enratio,sp_enratio) > 0.8 )&&(st->rate_control))||
                        (( tmp < 2.8 && max(res_enratio,sp_enratio) > 0.65 )&&(!st->rate_control)))
                {
                    PPP_MODE_E = 'B';
                }
            }
            else
            {
                if ((( tmp < 2.4 && max(res_enratio,sp_enratio) > 0.94 ) &&(st->rate_control))||
                        (( tmp < 4.5 && max(res_enratio,sp_enratio) > 0.5 ) &&(!st->rate_control)))
                {
                    PPP_MODE_E = 'B';
                }
            }
            /* -----End closed-loop Bump-Up */
        }
        else
        {
            PPP_MODE_E = 'B' ; /* Amplitude quantization is failing */
        }
    }
    else
    {
    }

    if (PPP_MODE_E == 'B')
    {
        st->bump_up = 1;

        free(CURRP_NQ);
        free(TMPDTFS);
        free(TMPDTFS2);
        free(TMPDTFS3);
        free(CURRP_Q_E);
        free(dtfs_temp);

        return;
    }

    if ( st->Q_to_F )
    {
        st->patterncount += st->pattern_m;

        if (st->patterncount >= 1000)
        {
            st->patterncount -= 1000;
            PPP_MODE_E = 'B';
            st->bump_up = 1;

            free(CURRP_NQ);
            free(TMPDTFS);
            free(TMPDTFS2);
            free(TMPDTFS3);
            free(CURRP_Q_E);
            free(dtfs_temp);

            return;
        }
    }

    /* packetization of the delta lag in PPP */
    if (PPP_MODE_E == 'Q')
    {
        Q_delta_lag = delta_lag_E+11; /* to make it positive always */
        push_indice( st, IND_DELTALAG, Q_delta_lag, 5 );
    }

    WIsyn(*dtfs_temp, CURRP_Q_E, lpc2, &(st->ph_offset_E), out, L_FRAME, 0);

    DTFS_copy(dtfs_temp, *CURRP_Q_E);
    st->prev_cw_en = DTFS_getEngy(*CURRP_NQ);

    /* Copy DTFS related parameters from 'dtfs_temp' to 'st' structure */
    st->dtfs_enc_lag = dtfs_temp->lag;
    st->dtfs_enc_nH = dtfs_temp->nH;
    st->dtfs_enc_nH_4kHz = dtfs_temp->nH_4kHz;
    st->dtfs_enc_upper_cut_off_freq_of_interest = dtfs_temp->upper_cut_off_freq_of_interest;
    st->dtfs_enc_upper_cut_off_freq = dtfs_temp->upper_cut_off_freq;

    mvr2r(dtfs_temp->a, st->dtfs_enc_a, MAXLAG_WI);
    mvr2r(dtfs_temp->b, st->dtfs_enc_b, MAXLAG_WI);

    free(CURRP_NQ);
    free(TMPDTFS);
    free(TMPDTFS2);
    free(TMPDTFS3);
    free(CURRP_Q_E);
    free(dtfs_temp);

    return;
}

static void SynthesisFilter(
    float *output,
    float *input,
    float *coef,
    float *memory,
    short order,
    short length
)
{
    short i, j;
    float acc;

    /* IIR filter for each subframe */
    for (i=0; i<length; i++)
    {
        for ( j = order-1, acc = *input++; j>0; j-- )
        {
            acc -= coef[j] * memory[j];
            memory[j] = memory[j - 1];
        }

        acc -= coef[0] * memory[0];
        *output++ = acc;
        memory[0] = acc;
    }

    return;
}

