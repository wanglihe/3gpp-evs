/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <stdlib.h>
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"


/*-------------------------------------------------------------------*
 * ppp_quarter_encoder()
 *
 * PPP quarter encoder
 *--------------------------------------------------------------------*/

short ppp_quarter_encoder(
    Encoder_State *st,                        /* i/o: encoder state structure      */
    DTFS_STRUCTURE *CURRCW_Q,                  /* o  : Quantized (amp/phase) DTFS */
    DTFS_STRUCTURE *TARGETCW,                  /* o  : DTFS with quant phase but unquant Amp */
    int   prevCW_lag,                 /* i  : previous lag */
    DTFS_STRUCTURE vCURRCW_NQ,                 /* i  : Unquantized DTFS */
    const float *curr_lpc,                  /* i  : LPCS */
    float *lastLgainE,                /* i/o: last low band gain */
    float *lastHgainE,                /* i/o: last high band gain */
    float *lasterbE,                  /* i/o: last ERB vector */
    DTFS_STRUCTURE PREV_CW_E                   /* i  : past DTFS */
)
{
    DTFS_STRUCTURE *PREVDTFS;

    float tmp, temp_pl, temp_l;
    int l;
    short returnFlag = 1;
    int POWER_IDX;    /* Codebook index for the power quantization for PPP */
    int AMP_IDX[2];   /* Codebook index for the Amplitude quantization for PPP */
    float Erot = 0.0, z = 0.0;

    PREVDTFS = DTFS_new();

    DTFS_copy( CURRCW_Q, vCURRCW_NQ );
    DTFS_copy( PREVDTFS, PREV_CW_E );

    l = CURRCW_Q->lag;
    temp_l = (float) CURRCW_Q->lag;
    temp_pl = (float) prevCW_lag;

    DTFS_adjustLag( PREVDTFS, l );

    z = ((L_FRAME-temp_l)*(temp_l+temp_pl))/(2*temp_l*temp_pl);

    Erot = (float) (temp_l - rint_new(temp_l*(z - floor(z))));

    DTFS_phaseShift( PREVDTFS, (float)(PI2*Erot/CURRCW_Q->lag)) ;
    DTFS_car2pol(PREVDTFS);

    /* Amplitude Quantization */
    DTFS_car2pol(CURRCW_Q); /* at this point currCW_q=curr_nq */

    returnFlag = DTFS_quant_cw(CURRCW_Q,prevCW_lag, curr_lpc, &POWER_IDX, AMP_IDX, lastLgainE, lastHgainE, lasterbE);

    push_indice( st, IND_AMP0, AMP_IDX[0], 6 );
    push_indice( st, IND_AMP1, AMP_IDX[1], 6 );
    push_indice( st, IND_POWER, POWER_IDX, 6 );

    DTFS_copy(TARGETCW,*CURRCW_Q);

    /* Copying phase spectrum over */
    mvr2r(PREVDTFS->b, CURRCW_Q->b, (short)(CURRCW_Q->lag>>1)+1 );

    DTFS_pol2car(CURRCW_Q);
    DTFS_pol2car(TARGETCW);
    tmp = DTFS_alignment_fine_new(*TARGETCW,*CURRCW_Q, 0.0) ;

    if (((tmp+3)>7) || ((tmp+3)<0))
    {
        tmp = 0;
        returnFlag = 0;
    }

    DTFS_phaseShift( CURRCW_Q,(float)(PI2*tmp/CURRCW_Q->lag) );

    push_indice( st, IND_GLOBAL_ALIGNMENT, (short) (tmp+3), 3 );

    free( PREVDTFS );

    return returnFlag;
}

/*-------------------------------------------------------------------*
  * set_ppp_mode()
  *
  * Determine if the current frame should be coded by PPP or not
  * Impose PPP - CELP - CELP pattern
  *-------------------------------------------------------------------*/

void set_ppp_mode(
    Encoder_State *st,                  /* i/o: encoder state structure */
    short *coder_type,          /* i/o: coder type      */
    const short noisy_speech_HO,      /* i  : SC-VBR noisy speech HO flag */
    const short clean_speech_HO,      /* i  : SC-VBR clean speech HO flag */
    const short NB_speech_HO,         /* i  : SC-VBR NB speech HO flag */
    const short localVAD,
    const short localVAD_he,          /* i  : HE-SAD flag without hangover */
    short *vad_flag
    ,short T_op[]                /* i  : open loop pitch lag */
    ,short sp_aud_decision1      /* i  : Speech Audio Decision */
)
{
    if( *vad_flag == 1 &&
            ( noisy_speech_HO == 1 || clean_speech_HO == 1 || NB_speech_HO == 1 ) &&
            ( localVAD == 0 || localVAD_he == 0 ) )

    {
        *coder_type = UNVOICED;
    }

    if( *coder_type == INACTIVE && *vad_flag == 0 && st->last_nelp_mode == 1 ) /* avoid HO frame go to GSC */
    {
        *coder_type = UNVOICED;
    }

    /* force the coder to NELP mode during the first five frames */
    /* this will indicate the decoder that the coder is operating in the VBR mode */
    if ( st->ini_frame < 5 )
    {
        *coder_type = UNVOICED;
        *vad_flag = 1;
    }

    /* Pattern PPP-CELP-CELP (pppcountE holds number of consecutive PPP frames) */
    if ( *coder_type != VOICED || st->last_coder_type == TRANSITION )
    {
        /* ensure no transient to PPP transition */
        st->pppcountE = 0;
    }
    else
    {
        /* current mode is voiced */
        st->pppcountE++;

        if ( ( st->pppcountE == 1 && st->last_last_ppp_mode != 1 && !st->rate_control ) ||
                ( st->pppcountE == 1 && st->mode_QQF ) )
        {
            st->ppp_mode = 1;
            st->core_brate = PPP_NELP_2k80;
        }
        else if ( st->pppcountE == 2 )
        {
            if ( st->last_ppp_mode == 1 && !st->mode_QQF )
            {
                /* QFF mode */
                st->ppp_mode = 0;
            }
            else
            {
                /* QQF Mode */
                st->ppp_mode = 1;
                st->core_brate = PPP_NELP_2k80;
            }
        }
        else
        {
            st->ppp_mode = 0;
            st->pppcountE = 0;
        }
    }

    if ( st->ppp_mode == 0 && st->set_ppp_generic == 1 )
    {
        st->set_ppp_generic = 0;
        *coder_type = GENERIC;
    }

    if( st->last_core == HQ_CORE )
    {
        st->ppp_mode = 0;
        st->set_ppp_generic = 0;
        *coder_type = TRANSITION;
    }

    if(st->last_ppp_mode && !st->ppp_mode && sp_aud_decision1 && st->bwidth == NB && st->Opt_SC_VBR) /*if it were about to go from ppp->HQ*/
    {
        st->avoid_HQ_VBR_NB = 1;
        *coder_type = GENERIC;
    }

    if(st->last_nelp_mode && sp_aud_decision1 && st->bwidth == NB && st->Opt_SC_VBR) /*if it were about to go from nelp->HQ*/
    {
        st->avoid_HQ_VBR_NB = 1;
        *coder_type = GENERIC;
    }

    if ((( st->old_pitch_buf[(2*NB_SUBFR)-1] > PPP_LAG_THRLD) || (T_op[1]>PPP_LAG_THRLD) || !st->last_Opt_SC_VBR ) && (st->ppp_mode==1))
    {
        st->ppp_mode = 0;
        st->core_brate = ACELP_7k20;
    }

    return;
}
