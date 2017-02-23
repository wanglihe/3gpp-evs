/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <stdlib.h>
#include "prot.h"
#include "cnst.h"
#include "rom_com.h"

/*-------------------------------------------------------------------*
 * ppp_voiced_decoder()
 *
 * Voiced decoder for SC-VBR
 *-------------------------------------------------------------------*/

void ppp_voiced_decoder(
    Decoder_State *st,           /* i/o: state structure */
    float *out,          /* o : residual signal */
    const float *lpc2,   /* i : current frame LPC */
    float *exc,          /* i : previous frame excitation */
    float *pitch         /* o : floating pitch values for each subframe */
    ,short bfi
)
{
    short k, delta_lag_D = 0;
    float upper_cut_off_freq_of_interest = 0, upper_cut_off_freq = 0;
    int pl, l;
    float interp_delay[3], temp_l, temp_pl, diff;

    DTFS_STRUCTURE *TMPDTFS = DTFS_new();
    DTFS_STRUCTURE *CURRP_Q_D = DTFS_new();

    DTFS_STRUCTURE *dtfs_temp = DTFS_new();

    if ( st->bwidth == NB )
    {
        upper_cut_off_freq_of_interest = 3300.0;
        upper_cut_off_freq = 4000.0;
    }
    else if ( st->bwidth == WB )
    {
        upper_cut_off_freq_of_interest = 4000.0;
        upper_cut_off_freq = 6400.0;
    }

    /* Initialization */
    if (st->firstTime_voiceddec)
    {
        st->firstTime_voiceddec=0;

        /* (st->PREV_CW_D) = DTFS_new();*/
        st->dtfs_dec_lag = 0;
        st->dtfs_dec_nH = 0;
        st->dtfs_dec_nH_4kHz = 0;
        st->dtfs_dec_upper_cut_off_freq_of_interest = 3300.0;
        st->dtfs_dec_upper_cut_off_freq = 4000.0;

        for(k = 0; k < MAXLAG_WI; k++)
        {
            st->dtfs_dec_a[k] = 0.0;
            st->dtfs_dec_b[k] = 0.0;
        }
    }
    pl = (int)min(rint_new(st->old_pitch_buf[(2*NB_SUBFR)-1]),MAX_LAG_PIT);
    delta_lag_D = (short) get_next_indice( st, 5 );

    l = min(MAX_LAG_PIT,pl+delta_lag_D-11);

    temp_pl = (float) pl;
    temp_l = (float) l;

    if (temp_pl != temp_l)
    {
        for(k=0; k<NB_SUBFR; k++)
        {
            /* do the linear pitch interp to drive the nb_post_filt */
            Interpol_delay(interp_delay,&(temp_pl), &(temp_l),k,frac_4sf);
            pitch[k] = min(MAX_LAG_PIT, max(19, interp_delay[0]));
        }
    }
    else
    {
        set_f(pitch, min(MAX_LAG_PIT, max(19, temp_l)), NB_SUBFR);
    }

    if (st->last_coder_type == UNVOICED)
    {
        pl = l; /* if prev frame was sil/uv*/
    }

    if (pl > (int)anint(1.85*l))
    {
        pl /= 2;
    }

    if (pl*2 <= PIT_MAX && pl <= (int)anint(0.54*l))
    {
        pl *= 2;
    }

    /* Restoring PPP memories when the last frame is non-PPP or full-rate PPP */
    if (st->last_ppp_mode_dec != 1)
    {
        DTFS_to_fs(exc-pl, pl, dtfs_temp, st->bwidth == WB ? (short)16000 : (short)8000 ,0);

        st->ph_offset_D = 0.0 ;

        /* Copy over PREV_CW_D into TMPDTFS */
        DTFS_copy(TMPDTFS, *dtfs_temp);

        DTFS_car2pol(TMPDTFS);

        st->lastLgainD = (float) log10(TMPDTFS->lag*DTFS_setEngyHarm(92.0,1104.5,0.0,1104.5,1.0,TMPDTFS));
        st->lastHgainD = (float) log10(TMPDTFS->lag*DTFS_setEngyHarm(1104.5,upper_cut_off_freq_of_interest,1104.5,upper_cut_off_freq,1.0,TMPDTFS));

        DTFS_to_erb(*TMPDTFS,st->lasterbD);
    }
    else
    {
        /* Copy DTFS related parameters from 'st' to 'dtfs_temp' structure */
        dtfs_temp->lag = st->dtfs_dec_lag;
        dtfs_temp->nH = st->dtfs_dec_nH;
        dtfs_temp->nH_4kHz = st->dtfs_dec_nH_4kHz;
        dtfs_temp->upper_cut_off_freq_of_interest = st->dtfs_dec_upper_cut_off_freq_of_interest;
        dtfs_temp->upper_cut_off_freq = st->dtfs_dec_upper_cut_off_freq;

        mvr2r(st->dtfs_dec_a, dtfs_temp->a, MAXLAG_WI);
        mvr2r(st->dtfs_dec_b, dtfs_temp->b, MAXLAG_WI);
    }

    CURRP_Q_D->lag = l;

    /* safety check in case of bit errors */
    if( CURRP_Q_D->lag <= 0 )
    {
        CURRP_Q_D->lag = 1;
        st->BER_detect = 1;
    }

    /* compute nH for lag */
    CURRP_Q_D->nH = (int)floor(upper_cut_off_freq/(INT_FS_12k8/CURRP_Q_D->lag));
    diff = INT_FS_12k8 / CURRP_Q_D->lag ;

    if (upper_cut_off_freq-(diff*CURRP_Q_D->nH)>=diff)
    {
        CURRP_Q_D->nH++;
    }

    CURRP_Q_D->nH_4kHz = (int)floor(4000.0/(INT_FS_12k8/CURRP_Q_D->lag));

    if (4000.0 - (diff*CURRP_Q_D->nH_4kHz) >= diff)
    {
        CURRP_Q_D->nH_4kHz++;
    }

    CURRP_Q_D->upper_cut_off_freq = dtfs_temp->upper_cut_off_freq;
    CURRP_Q_D->upper_cut_off_freq_of_interest = dtfs_temp->upper_cut_off_freq_of_interest;

    if ( bfi == 0 )
    {
        ppp_quarter_decoder( st, CURRP_Q_D, dtfs_temp->lag, &(st->lastLgainD), &(st->lastHgainD),
                             &(st->lasterbD[0]),
                             bfi,
                             *dtfs_temp );
    }

    WIsyn(*dtfs_temp, CURRP_Q_D, lpc2, &(st->ph_offset_D), out, (short) L_FRAME ,0 );

    DTFS_copy(dtfs_temp, *CURRP_Q_D);

    /* Copy DTFS related parameters from 'dtfs_temp' to 'st' structure */
    st->dtfs_dec_lag = dtfs_temp->lag;
    st->dtfs_dec_nH = dtfs_temp->nH;
    st->dtfs_dec_nH_4kHz = dtfs_temp->nH_4kHz;
    st->dtfs_dec_upper_cut_off_freq_of_interest = dtfs_temp->upper_cut_off_freq_of_interest;
    st->dtfs_dec_upper_cut_off_freq = dtfs_temp->upper_cut_off_freq;

    mvr2r(dtfs_temp->a, st->dtfs_dec_a, MAXLAG_WI);
    mvr2r(dtfs_temp->b, st->dtfs_dec_b, MAXLAG_WI);

    free(TMPDTFS);
    free(CURRP_Q_D);
    free(dtfs_temp);

    return;
}
