/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "stat_enc.h"
#include "prot.h"


/*-------------------------------------------------------------------*
* coderLookAheadInnovation()
*
*
*-------------------------------------------------------------------*/

void coderLookAheadInnovation(
    const float A[],            /* input: coefficients NxAz[M+1]   */
    int *pT,                    /* out:   pitch                    */
    HANDLE_PLC_ENC_EVS st,      /* i/o:   coder memory state       */
    float *speechLookAhead,
    float *old_exc,
    int L_subfr,
    int L_frame
)
{
    int i;
    float *exc, exc_buf[L_EXC_MEM+2*L_SUBFR+8] = { 0.0f };
    int   T0=0;
    short prev_pitch;
    float ps, alp, max_ps;
    short subfr_len;
    short search_range;

    search_range = 9;

    /* Framing parameters */
    if( L_frame < L_FRAME16k )
    {
        subfr_len = (short)(1.75*L_subfr);
    }
    else
    {
        subfr_len = (short)(2*L_subfr);
    }

    /*------------------------------------------------------------------------*
     * Initialize buffers                                                     *
     *------------------------------------------------------------------------*/

    /* set excitation memory */
    exc = exc_buf + L_EXC_MEM + 8;
    mvr2r( old_exc, exc_buf, L_EXC_MEM+8);

    /*------------------------------------------------------------------------*
     * - Get residual signal and target at lookahead part.                    *
     *------------------------------------------------------------------------*/

    /* find LP residual signal for look-ahead part */
    getLookAheadResSig( speechLookAhead, A, exc, L_frame, L_subfr, M, 2 );

    /* Initialize excitation buffer */
    prev_pitch = st->T0_4th;
    /* find target signals */
    /* find best candidate of pitch lag */
    {
        max_ps = -1.0e10;
        T0 = st->T0_4th;
        for( i=-search_range; i<search_range; i++ )
        {
            if( prev_pitch+i>st->pit_max || prev_pitch+i<st->pit_min )
            {
                continue;
            }
            ps = dotp( exc, &exc[-(prev_pitch+i)], subfr_len );
            alp = dotp( &exc[-(prev_pitch+i)], &exc[-(prev_pitch+i)], subfr_len );
            ps /= (float)sqrt( alp + 1.0e-10 );
            if( max_ps < ps )
            {
                max_ps = ps;
                T0 = prev_pitch + i;
            }
        }
        if( max_ps < 0.0 )
        {
            T0 = st->T0_4th;
        }
    }

    pT[0] = T0;

    return;
}


/*-------------------------------------------------------------------*
* enc_prm_side_Info()
*
*
*-------------------------------------------------------------------*/

void enc_prm_side_Info(
    HANDLE_PLC_ENC_EVS hPlc_Ext,
    Encoder_State *st
)
{
    int diff_pitch;
    short bits_per_subfr, search_range;

    bits_per_subfr = 4;
    search_range = 8;

    if( hPlc_Ext->nBits > 1 )
    {
        push_next_indice(st, 1, 1);

        diff_pitch = hPlc_Ext->T0 - hPlc_Ext->T0_4th;

        if( ( diff_pitch > search_range-1 ) || ( diff_pitch < -search_range+1 ) )
        {
            diff_pitch = -8;
        }

        push_next_indice(st, (diff_pitch+search_range), bits_per_subfr);
    }
    else
    {
        push_next_indice(st, 0, 1);
    }

    return;
}

/*-------------------------------------------------------------------*
* encoderSideLossSimulation()
*
* Encoder side loss simulation
*-------------------------------------------------------------------*/

void encoderSideLossSimulation(
    Encoder_State *st,
    HANDLE_PLC_ENC_EVS hPlc_Ext,
    float *lsf_q,
    float stab_fac,
    int calcOnlylsf,
    int L_frame
)
{
    float lspLocal[M];
    float const* lsfBase;                      /* base for differential lsf coding */

    /*************************************************************
     * Decoder state could be stored with memcpy,
     * since Decoder_State does not contain pointer member.
     *************************************************************/

    /* Decoder State Update */
    lsf2lsp( lsf_q, lspLocal, M, st->sr_core );

    lsfBase = PlcGetlsfBase( st->lpcQuantization, st->narrowBand, st->sr_core );

    mvr2r( st->mem_MA, hPlc_Ext->mem_MA, M );

    /* lsf parameter processing for concealment */
    updatelsfForConcealment( hPlc_Ext, lsf_q);
    hPlc_Ext->stab_fac = stab_fac;

    /* Update Decoder State for the loss simulation at the next frame */
    mvr2r( lsf_q, hPlc_Ext->lsfold, M );
    mvr2r( lspLocal, hPlc_Ext->lspold, M );

    if (calcOnlylsf)
    {
        /* lsf concealment simulation */
        getConcealedlsf( hPlc_Ext, lsfBase, L_frame, st->clas);
        hPlc_Ext->T0 = hPlc_Ext->T0_4th;
    }
    else
    {
        float AqCon[(NB_SUBFR16k+1)*(M+1)];
        float *speechLookAhead;
        float old_exc[L_EXC_MEM+8];

        /*                 Initialize pointers here                  */
        mvr2r( hPlc_Ext->old_exc, old_exc, 8 );
        mvr2r( hPlc_Ext->LPDmem->old_exc, &old_exc[8], L_EXC_MEM );
        speechLookAhead = &( st->speech_enc_pe[L_frame] );

        /* lsf concealment simulation */
        getConcealedLP( hPlc_Ext, AqCon, lsfBase, st->sr_core, st->clas, L_frame);

        /* apply encoder side PLC simulation */
        hPlc_Ext->pit_min = st->pit_min;
        hPlc_Ext->pit_max = st->pit_max;
        coderLookAheadInnovation( AqCon, &(hPlc_Ext->T0), hPlc_Ext, speechLookAhead, old_exc, L_SUBFR, st->L_frame );

    }

    return;
}

/*-------------------------------------------------------------------*
* GplcTcxEncSetup()
*
*
*-------------------------------------------------------------------*/

void GplcTcxEncSetup(
    Encoder_State *st,
    HANDLE_PLC_ENC_EVS hPlc_Ext)
{
    hPlc_Ext->T0_4th = st->tcxltp_pitch_int;

    return;
}

/*-------------------------------------------------------------------*
* encSideSpecPowDiffuseDetector()
*
*
*-------------------------------------------------------------------*/

short encSideSpecPowDiffuseDetector(
    float *lsf_ref,
    float *lsf_con,
    int sr_core,
    float *prev_lsf4_mean,
    short sw
    ,short coder_type
)
{
    float lsf_mod[M];
    float dist1, dist2, cum_dist1, cum_dist2;
    float lsf4_mean;
    float th;
    float th_dif_lsf4_mean;
    short idx;
    int cnt_imprv, i;

    /* calculate the mean of the lowest 4 lsfs */
    lsf4_mean = 0;

    for(i = 0; i < 4; i++)
    {
        lsf4_mean += lsf_ref[i];
    }
    lsf4_mean /= 4.0f;

    if(sw)
    {
        mvr2r( lsf_con, lsf_mod, M );

        modify_lsf( lsf_mod, M, sr_core
                    , 1
                  );

        cum_dist1 = 0;
        cum_dist2 = 0;

        cnt_imprv = 0;


        for(i = 0; i < M; i++)
        {
            dist1 = (lsf_con[i] - lsf_ref[i]) * (lsf_con[i] - lsf_ref[i]);
            dist2 = (lsf_mod[i] - lsf_ref[i]) * (lsf_mod[i] - lsf_ref[i]);
            cum_dist1 += dist1;
            cum_dist2 += dist2;

            if(dist1 > dist2)
            {
                cnt_imprv++;
            }
        }

        th = 800;
        th_dif_lsf4_mean = 90;

        if(sr_core == 16000)
        {
            th *= 1.25;
            th_dif_lsf4_mean *= 1.25;
        }


        if(cum_dist1 > cum_dist2 * 1.15
                && lsf4_mean - *prev_lsf4_mean > th_dif_lsf4_mean
                && *prev_lsf4_mean < th
                && cnt_imprv > 2
                && coder_type == GENERIC
          )
        {
            idx = 1;
        }
        else
        {
            idx = 0;
        }
    }
    else
    {
        idx = 0;
    }

    /* update parameters */
    *prev_lsf4_mean = lsf4_mean;

    return idx;
}

/*-------------------------------------------------------------------*
* updateSpecPowDiffuseIdx()
*
*
*-------------------------------------------------------------------*/

void updateSpecPowDiffuseIdx(
    Encoder_State *st,
    const float gain_pitch_buf[],   /* i  : gain pitch values   */
    const float gain_code_buf[]     /* i  : gain pitch values   */
)
{
    float min_gp;
    int k;

    st->mean_gc[1] = gain_code_buf[0];
    min_gp = gain_pitch_buf[0];

    for(k = 1; k < 4; k++)
    {
        st->mean_gc[1] += gain_code_buf[k];

        if( gain_pitch_buf[k] < min_gp )
        {
            min_gp = gain_pitch_buf[k];
        }
    }

    if(st->mean_gc[1] / (st->mean_gc[0] + 1e-6) < 1.098 || min_gp > 0.82)
    {
        st->glr_idx [0]= 0;
    }
    st->mean_gc[0] = st->mean_gc[1];

    return;
}


/*-------------------------------------------------------------------*
* getConcealedlsf()
*
*
*-------------------------------------------------------------------*/

void getConcealedlsf(
    HANDLE_PLC_ENC_EVS memDecState,
    const float lsfBase[],
    int L_frame,
    int last_good
)
{
    float *lsf = memDecState->lsf_con;

    dlpc_bfi( L_frame, &lsf[0], memDecState->lsfold, last_good,
              1 /*assumes packet loss */, memDecState->mem_MA, memDecState->mem_AR, &(memDecState->stab_fac), memDecState->lsf_adaptive_mean,
              1, NULL, 0, NULL, NULL, lsfBase);

    return;
}
