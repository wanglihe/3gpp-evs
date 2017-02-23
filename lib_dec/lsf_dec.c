/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <stdlib.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"
#include "basop_proto_func.h"

/*---------------------------------------------------------------------*
 * Local functions
 *---------------------------------------------------------------------*/
static void lsf_mid_dec( Decoder_State *st, const float int_fs, float qlsp0[], float qlsp1[], short coder_type, float qlsp[],
                         const long core_brate, short ppp_mode, short nelp_mode, short prev_bfi, short *mid_lsf_int, short safety_net );

static void dqlsf_CNG( Decoder_State *st, float *lsf_q, unsigned int *p_offset_scale1, unsigned int *p_offset_scale2, short *p_no_scales );

/*---------------------------------------------------------------------*
 * lsf_dec()
 *
 * LSF decoder
 *---------------------------------------------------------------------*/

void lsf_dec(
    Decoder_State *st,          /* i/o: State structure                      */
    const short tc_subfr,       /* i  : TC subframe index                    */
    const short L_frame,        /* i  : length of the frame                  */
    const short coder_type,     /* i  : coding type                          */
    const short bwidth,         /* i  : input signal bandwidth               */
    float *Aq,            /* o  : quantized A(z) for 4 subframes       */
    short *LSF_Q_prediction,  /* o  : LSF prediction mode              */
    float *lsf_new,       /* o  : de-quantized LSF vector              */
    float *lsp_new,       /* o  : de-quantized LSP vector              */
    float *lsp_mid        /* o  : de-quantized mid-frame LSP vector    */
)
{
    short i, nBits = 0;
    float int_fs;
    float tmp_old[M+1], tmp_new[M+1], enr_old = 0.0f, enr_new = 0.0f;
    float lsf_diff = 0.0f;

    /* initialize */
    if( L_frame == L_FRAME )
    {
        int_fs = INT_FS_12k8;
    }
    else /* L_frame == L_FRAME16k */
    {
        int_fs = INT_FS_16k;
    }

    /* Find the number of bits for LSF quantization */
    if ( st->core_brate == SID_2k40 )
    {
        nBits = LSF_BITS_CNG;
    }
    else
    {
        if ( st->nelp_mode_dec == 0 && st->ppp_mode_dec == 0 )
        {
            nBits = LSF_bits_tbl[LSF_BIT_ALLOC_IDX(st->core_brate, coder_type)];
        }
        else if ( st->nelp_mode_dec == 1 )
        {
            if ( coder_type == UNVOICED )
            {
                if ( bwidth == NB )
                {
                    nBits = 32;

                }
                else
                {
                    nBits = 30;


                }
            }
        }
        else if ( st->ppp_mode_dec == 1 )
        {
            nBits = 26;

        }
    }

    /* LSF de-quantization */
    lsf_end_dec( st, coder_type, st->bwidth, nBits, lsf_new, st->mem_AR,st->mem_MA, int_fs, st->core_brate,
                 &st->offset_scale1[0][0], &st->offset_scale2[0][0], &st->offset_scale1_p[0][0], &st->offset_scale2_p[0][0],
                 &st->no_scales[0][0], &st->no_scales_p[0][0], &st->safety_net, NULL, LSF_Q_prediction, NULL );

    /* convert quantized LSFs to LSPs */
    lsf2lsp( lsf_new, lsp_new, M, int_fs );

    if( st->core_brate == SID_2k40 )
    {
        /* return if SID frame (conversion to A(z) done in the calling function) */
        return;
    }

    /*-------------------------------------------------------------------------------------*
     * FEC - update adaptive LSF mean vector
     *-------------------------------------------------------------------------------------*/

    for ( i=0; i<M; i++ )
    {
        st->lsf_adaptive_mean[i] = (st->lsfoldbfi1[i] + st->lsfoldbfi0[i] + lsf_new[i]) / 3;
    }

    if( st->prev_bfi && (coder_type == TRANSITION) && (tc_subfr == (L_frame-L_SUBFR)) )
    {
        lsf_diff = int_fs / (float)(2*(M+1));
        st->lsf_old[0] = lsf_diff;

        for ( i=1; i<M; i++ )
        {
            st->lsf_old[i] = st->lsf_old[i-1] + lsf_diff;
        }
        lsf2lsp( st->lsf_old, st->lsp_old, M, int_fs );
    }

    /*-------------------------------------------------------------------------------------*
     * Mid-frame LSF decoding
     * LSP interpolation and conversion of LSPs to A(z)
     *-------------------------------------------------------------------------------------*/
    if( st->rate_switching_reset )
    {
        /* extrapolation in case of unstable LSF convert */
        mvr2r( lsp_new, st->lsp_old, M );
        mvr2r( lsf_new, st->lsf_old, M );
    }

    lsf_mid_dec( st, int_fs, st->lsp_old, lsp_new, coder_type, lsp_mid, st->core_brate, st->ppp_mode_dec, st->nelp_mode_dec,
                 st->prev_bfi, &(st->mid_lsf_int), st->safety_net );

    if ( !( st->prev_bfi && (coder_type == TRANSITION) && (tc_subfr == (L_frame-L_SUBFR)) ) )
    {
        if (st->prev_bfi)
        {
            /* check, if LSP interpolation can be relaxed */
            lsp2a_stab( st->lsp_old, tmp_old, M);
            enr_old = enr_1_Az( tmp_old, 2*L_SUBFR );

            lsp2a_stab( lsp_new, tmp_new, M);
            enr_new = enr_1_Az( tmp_new, 2*L_SUBFR );
        }

        if (st->prev_bfi)
        {
            if ( enr_new/enr_old < 0.3f )
            {
                st->relax_prev_lsf_interp = -1;
                if ( st->clas_dec == UNVOICED_CLAS || st->clas_dec == SIN_ONSET || st->clas_dec == INACTIVE_CLAS || coder_type == GENERIC || coder_type == TRANSITION )
                {
                    st->relax_prev_lsf_interp = 1;
                }
            }
        }
    }

    if( st->last_core == HQ_CORE && st->core == ACELP_CORE )
    {
        /* update old LSPs/LSFs in case of HQ->ACELP core switching */
        mvr2r( lsp_mid, st->lsp_old, M );
        lsp2lsf( lsp_mid, st->lsf_old, M, int_fs );
    }

    /* LSP interpolation and conversion of LSPs to A(z) */
    int_lsp4( L_frame, st->lsp_old, lsp_mid, lsp_new, Aq, M, st->relax_prev_lsf_interp );

    /*------------------------------------------------------------------*
     * Check LSF stability (distance between old LSFs and current LSFs)
     *------------------------------------------------------------------*/

    st->stab_fac = lsf_stab( lsf_new, st->lsf_old, 0, st->L_frame);

    return;
}



/*------------------------------------------------------------------------------------------*
 * lsf_end_dec()
 *
 * De-quantize frame end LSF vector
 *------------------------------------------------------------------------------------------*/

void lsf_end_dec(
    Decoder_State *st,                /* i/o: decoder state structure                 */
    const short coder_type_org,     /* i  : coding type                             */
    const short bwidth,             /* i  : input signal bandwidth                  */
    const short nBits_in,           /* i  : number of bits used for ISF quantization*/
    float *qlsf,              /* o  : quantized LSFs in the cosine domain     */
    float *mem_AR,            /* i/o: quantizer memory for AR model           */
    float *mem_MA,            /* i/o: quantizer memory for MA model           */
    const float int_fs,             /* i  : sampling frequency                      */
    long  core_brate,         /* i  : Coding Bit Rate                         */
    unsigned int   *p_offset_scale1,
    unsigned int   *p_offset_scale2,
    unsigned int   *p_offset_scale1_p,
    unsigned int   *p_offset_scale2_p,
    short *p_no_scales,
    short *p_no_scales_p,
    short *safe_net,
    int   *lpc_param,
    short *LSF_Q_prediction,  /* o  : LSF prediction mode                     */
    int * nb_indices
)
{
    float pred0[M];               /* Prediction for the safety-net quantizer (usually mean)*/
    float pred1[M], pred2[M];     /* Prediction for the predictive quantizer*/
    short stages0;                /* Amount of stages used by safety-net quantizer*/
    short stages1;                /* Amount of stages used by predictive quantizer*/
    short levels0[MAX_VQ_STAGES]; /* Sizes of different codebook stages for safety-net quantizer*/
    short levels1[MAX_VQ_STAGES]; /* Sizes of different codebook stages for predictive quantizer*/
    short i;
    short TCQIdx[M/2+4];
    short bits0[MAX_VQ_STAGES], bits1[MAX_VQ_STAGES];
    int   cumleft;
    short lindice[MAX_VQ_STAGES+3]; /* Predictor selector needs 1 bit and the LVQ indice uses 3 shorts */
    short mode_lvq, mode_lvq_p;
    short safety_net, predmode, stages, *levels;
    const short *Bit_alloc1 = NULL, *bits;
    short num_bits;
    int * p_lpc_param;

    int nr_ind;
    short nBits;

    short coder_type;
    nBits = nBits_in;

    if( (coder_type_org == GENERIC) && (int_fs == INT_FS_16k) && (st->codec_mode == MODE1) )
    {
        coder_type = (short)get_next_indice( st, 1 );
        coder_type += 2;
        if( coder_type == GENERIC )
        {
            nBits--;
        }
    }
    else
    {
        coder_type = coder_type_org;
    }

    /*--------------------------------------------------------------------------------*
     * LSF de-quantization of SID frames
     *--------------------------------------------------------------------------------*/

    if ( core_brate == SID_2k40 )
    {
        dqlsf_CNG( st, qlsf, p_offset_scale1, p_offset_scale2, p_no_scales );
        v_sort( qlsf, 0, M-1);
        reorder_lsf( qlsf, MODE1_LSF_GAP, M, int_fs );

        return;
    }

    predmode = find_pred_mode(coder_type, bwidth, int_fs, &mode_lvq, &mode_lvq_p, st->total_brate);

    /*----------------------------------------------------------------*
     * Calculate number of stages and levels for each stage based on the allowed bit allocation
     * (subtract one bit for LSF predictor selection)
     *----------------------------------------------------------------*/

    lsf_allocate( nBits-(predmode>>1), mode_lvq, mode_lvq_p, &stages0, &stages1, levels0, levels1, bits0, bits1 );

    /*--------------------------------------------------------------------------*
     * Select safety_net or predictive mode
     *--------------------------------------------------------------------------*/

    p_lpc_param = lpc_param;

    nr_ind = 0;
    if( predmode == 0 )
    {
        safety_net = 1;
    }
    else if ( predmode == 1 )
    {
        safety_net = 0;
    }
    else
    {
        if( st->codec_mode == MODE2 )
        {
            nr_ind ++;
            /* read from param_lpc */
            safety_net = p_lpc_param[0];
            p_lpc_param++;
        }
        else
        {
            safety_net = (short)get_next_indice( st, 1 );
        }
    }

    *safe_net = safety_net;

    /*--------------------------------------------------------------------------*
     * Read indices from array
     *--------------------------------------------------------------------------*/

    if ( safety_net )
    {
        stages = stages0;
        levels = levels0;
        bits = bits0;
    }
    else
    {
        stages = stages1;
        levels = levels1;
        bits = bits1;
    }

    if( st->codec_mode == MODE2 )
    {
        /* VOICED_WB@16kHz */
        if ( int_fs == INT_FS_16k && coder_type == VOICED )
        {
            *nb_indices = 10;
            for(i=0; i<*nb_indices; i++)
            {
                TCQIdx[i] = (short)lpc_param[i];
            }
        }
        else
        {
            for ( i=0; i<stages-1; i++ )
            {
                num_bits = bits[i];
                lindice[i+1] = *p_lpc_param++;
                nr_ind ++;
            }

            cumleft = levels[stages-1];
            while ( cumleft > 0 )
            {
                if ( cumleft > LEN_INDICE )
                {
                    cumleft -= LEN_INDICE;
                    num_bits = LEN_INDICE;
                }
                else
                {
                    num_bits = (short)cumleft;
                    cumleft = 0;
                }

                lindice[i+1] = *p_lpc_param++;
                nr_ind++;
                i++;
            }
            *nb_indices = nr_ind;
        }
    }
    else
    {
        /* VOICED_WB@16kHz */
        if ( int_fs == INT_FS_16k && coder_type == VOICED )
        {
            Bit_alloc1 = &BC_TCVQ_BIT_ALLOC_40B[1];
            TCQIdx[0] = safety_net;
            for ( i=0; i<M/2+3; i++ )
            {
                TCQIdx[i+1] = (short)get_next_indice( st, Bit_alloc1[i] );
            }
        }
        else
        {
            for ( i=0; i<stages-1; i++ )
            {
                num_bits = bits[i];
                lindice[i+1] = (short)get_next_indice( st, num_bits );
            }

            cumleft = levels[stages-1];
            while ( cumleft > 0 )
            {
                if ( cumleft > LEN_INDICE )
                {
                    cumleft -= LEN_INDICE;
                    num_bits = LEN_INDICE;
                }
                else
                {
                    num_bits = (short)cumleft;
                    cumleft = 0;
                }

                lindice[i+1] = (short)get_next_indice( st, num_bits );
                i++;
            }
        }
    }

    if( st->reset_mem_AR == 1 )
    {
        for( i=0; i<M; i++ )
        {
            st->mem_AR[i] = ModeMeans[mode_lvq][i];
        }
        st->reset_mem_AR = 0;
    }

    /*------------------------------------------------------------------------------------------*
     * De-quantize LSF vector
     *------------------------------------------------------------------------------------------*/

    *LSF_Q_prediction = SAFETY_NET;

    /* VOICED_WB@16kHz */
    if( int_fs == INT_FS_16k && coder_type == VOICED )
    {
        /* BC-TCVQ decoder */
        safety_net = qlsf_ARSN_tcvq_Dec_16k ( qlsf, TCQIdx, nBits-1 );

        /* Update mem_MA */
        mvr2r( qlsf, mem_MA, M );

        if( safety_net )
        {
            mvr2r( ModeMeans[mode_lvq], pred0, M );
        }
        else
        {
            for(i = 0; i < M; i++)
            {
                pred0[i] = ModeMeans[mode_lvq][i] + Predictors[mode_lvq_p][i]*(mem_AR[i]-ModeMeans[mode_lvq][i]);
            }
            *LSF_Q_prediction = AUTO_REGRESSIVE;
        }
        v_add( qlsf, pred0, qlsf, M );
    }
    else
    {
        /* Safety-net */
        mvr2r( ModeMeans[mode_lvq], pred0, M );

        /* for mem_MA update */
        for (i=0; i<M; i++)
        {
            pred1[i] = pred0[i]+MU_MA*mem_MA[i];
        }

        if( safety_net )
        {
            /* LVQ */
            st->BER_detect = st->BER_detect |
                             vq_dec_lvq( 1, qlsf, &lindice[1], stages0, M, mode_lvq, levels0[stages0-1],
                                         p_offset_scale1, p_offset_scale2, p_offset_scale1_p, p_offset_scale2_p, p_no_scales, p_no_scales_p );

            v_add( qlsf, pred0, qlsf, M );
            v_sub( qlsf, pred1, mem_MA, M);
        }
        else
        {
            st->BER_detect = st->BER_detect |
                             vq_dec_lvq( 0, qlsf, &lindice[1], stages1, M, mode_lvq_p, levels1[stages1-1],
                                         p_offset_scale1, p_offset_scale2, p_offset_scale1_p, p_offset_scale2_p, p_no_scales, p_no_scales_p );

            if( predmode == 1 ) /* MA only */
            {
                mvr2r(qlsf, mem_MA, M);
                v_add( qlsf, pred1, qlsf, M );
                *LSF_Q_prediction = MOVING_AVERAGE;
            }
            else
            {
                /* AR  */
                for ( i=0; i<M; i++ )
                {
                    pred2[i] = pred0[i] + Predictors[mode_lvq_p][i] * (mem_AR[i] - pred0[i]);
                }
                v_add( qlsf, pred2, qlsf, M );
                v_sub( qlsf, pred1, mem_MA, M );
                *LSF_Q_prediction = AUTO_REGRESSIVE;
            }
        }
    }


    /*--------------------------------------------------------------------------*
     * Sort the quantized vector
     * Verify stability
     * Update AR-predictor memory
     *--------------------------------------------------------------------------*/

    /* Sort the quantized vector */
    v_sort( qlsf, 0, M-1 );

    /* Verify stability */
    reorder_lsf( qlsf, MODE1_LSF_GAP, M, int_fs );

    /* Update predictor memory */
    mvr2r( qlsf, mem_AR, M );

    st->mode_lvq = mode_lvq;


    return;
}


/*-------------------------------------------------------------------*
 * lsf_mid_dec()
 *
 * Decode mid-frame LSFs
 *-------------------------------------------------------------------*/

void lsf_mid_dec(
    Decoder_State *st,        /* i/o: decoder state structure             */
    const float int_fs,     /* i  : internal (ACELP) sampling frequency */
    float qlsp0[],    /* i  : quantized LSPs from frame beginning */
    float qlsp1[],    /* i  : quantized LSPs from frame end       */
    short coder_type, /* i  : Coder type                          */
    float qlsp[],     /* o  : quantized LSPs                      */
    const  long core_brate, /* i  : core bitrate                        */
    short ppp_mode,
    short nelp_mode,
    short prev_bfi,
    short *mid_lsf_int,
    short safety_net
)
{
    short j, idx;
    short nb_bits;
    float qlsf0[M], qlsf1[M], qlsf[M];
    const float *ratio = NULL;
    short bad_spacing;

    bad_spacing = 0;

    /* Convert LSPs to LSFs */
    lsp2lsf( qlsp0, qlsf0, M, int_fs );
    lsp2lsf( qlsp1, qlsf1, M, int_fs );

    /* Codebook selection */
    if( ppp_mode == 1 )
    {
        nb_bits = 1;
        ratio = &(tbl_mid_voi_wb_1b[0]);
    }
    else if( nelp_mode == 1 )
    {
        nb_bits = 4;
        ratio = &(tbl_mid_unv_wb_4b[0]);
    }
    else
    {
        nb_bits = mid_LSF_bits_tbl[LSF_BIT_ALLOC_IDX(core_brate, coder_type)];

        /* codebook selection */
        if( coder_type == VOICED )
        {
            switch ( nb_bits )
            {
            case 5:
            {
                ratio = tbl_mid_voi_wb_5b;
                break;
            }
            case 4:
            {
                ratio = tbl_mid_voi_wb_4b;
                break;
            }
            }
        }
        else if( coder_type == UNVOICED )
        {
            ratio = tbl_mid_unv_wb_5b;
        }
        else
        {
            /* GENERIC, TRANSITION, AUDIO and INACTIVE */
            switch ( nb_bits )
            {
            case 5:
            {
                ratio = tbl_mid_gen_wb_5b;
                break;
            }
            case 2:
            {
                ratio = tbl_mid_gen_wb_2b;
                break;
            }
            }
        }
    }

    /* Retrieve mid-frame LSF index */
    idx = (short)get_next_indice( st, nb_bits );

    /* Calculation of mid-LSF vector */
    for(j=0; j<M; j++)
    {
        qlsf[j] = (1.0f - ratio[idx*M+j]) * qlsf0[j] + ratio[idx*M+j] * qlsf1[j];
    }

    /* check for incorrect LSF ordering */
    if( *mid_lsf_int == 1 )
    {
        for (j=1; j<M; j++)
        {
            if ( qlsf[j] < qlsf[j-1] )
            {
                bad_spacing = 1;
                break;
            }
        }
    }

    /* Redo mid-LSF interpolation with 0.4 in case of LSF instability */
    if( prev_bfi || ( *mid_lsf_int == 1 && bad_spacing ) )
    {
        for(j=0; j<M; j++)
        {
            /* redo mid-LSF interpolation with 0.4 */
            qlsf[j] = 0.4f * qlsf0[j] + 0.6f * qlsf1[j];

            /* ensure correct ordering of LSF indices */
            if ( j > 0 && j < M && qlsf[j] < qlsf[j-1] + LSF_GAP_MID  )
            {
                qlsf[j] = qlsf[j-1] + LSF_GAP_MID;
            }

        }
    }
    else
    {
        /* otherwise, use regular LSF spacing and ordering as in the encoder */
        for (j=0; j<M; j++)
        {
            if ( j > 0 && j < M && qlsf[j] < qlsf[j-1] + LSF_GAP_MID )
            {
                qlsf[j] = qlsf[j-1] + LSF_GAP_MID;
            }

        }
    }

    if( prev_bfi )
    {
        /* continue redoing mid-LSF interpolation with 0.4 in order not to propagate the error */
        *mid_lsf_int = 1;
    }

    if( safety_net )
    {
        /* safety-net encountered -> stop redoing mid-LSF interpolation with 0.4 */
        *mid_lsf_int = 0;
    }

    reorder_lsf( qlsf, LSF_GAP_MID, M, int_fs );
    /* convert back to LSPs */
    lsf2lsp( qlsf, qlsp, M, int_fs );

    return;
}

/*----------------------------------------------------------------------------------------------*
 * dqlsf_CNG()
 *
 * LSF de-quantizer for SID frames (uses 28 bits, 4 for VQ, 24 for LVQ)
 *
 * Note:
 * LP-CNG LSF decoder does not need to know the sampling rate,
 * the sampling rate data is embedded inside the LSF coefficients
 * If the highest order LSF coefficient (lsf_q[M-1]) is smaller than 6350 then Fs=12.8kHz
 * If the highest order LSF coefficient (lsf_q[M-1]) is larger than 6350 then Fs=16kHz
 *----------------------------------------------------------------------------------------------*/

static void dqlsf_CNG(
    Decoder_State *st,                        /* i/o: decoder state structure   */
    float *lsf_q,                     /* o  : decoded LSFs */
    unsigned int *p_offset_scale1,           /* i  : offset for 1st LVQ subvector */
    unsigned int *p_offset_scale2,           /* i  : offset for second LVQ subvector */
    short *p_no_scales                /* i  : number of scales for LVQ struct */
)
{
    short indice[4];

    indice[0] = (short)get_next_indice( st, 4 );
    indice[1] = (short)get_next_indice( st, LEN_INDICE );
    indice[2] = (short)get_next_indice( st, LSF_BITS_CNG - 4 - LEN_INDICE );

    st->BER_detect = st->BER_detect |
                     deindex_lvq_cng( &indice[1], lsf_q, indice[0], LSF_BITS_CNG-4, p_offset_scale1, p_offset_scale2, p_no_scales );
    /* The sampling frequency of the LP-CNG frame can be determined by checking the value of the highest order LSF
       coefficient (last coefficient of lsf_q). If the last decoded LSF coefficient (lsf_q[15]) is larger than 6350
       the decoded frame is WB2 with sampling rate of 16 kHz, otherwise it is sampled at 12.8kHz and contains
       either NB or WB LSF data.  */

    v_add( lsf_q, &CNG_SN1[indice[0]*M], lsf_q, M );

    if( ((st->L_frame == L_FRAME16k)&&(lsf_q[M-1]<=WB_LIMIT_LSF)) || ((st->L_frame<L_FRAME16k)&&(lsf_q[M-1]>WB_LIMIT_LSF)) )
    {
        st->BER_detect = 1;
    }
    return;
}
