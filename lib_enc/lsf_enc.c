/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <stdlib.h>
#include "options.h"
#include "cnst.h"
#include "rom_enc.h"
#include "rom_com.h"
#include "prot.h"
#include "basop_proto_func.h"

/*-----------------------------------------------------------------*
 * Local constants
 *-----------------------------------------------------------------*/

#define MSVQ_MAXCNT      3000     /* was 300 */
#define LOWEMPH_LSFI     1.0f     /* bigger value means more weight for lower frequencies */

/*---------------------------------------------------------------------*
 * Local functions
 *---------------------------------------------------------------------*/

static void lsfq_CNG( Encoder_State *st, const float *lsf, const float *wghts, float *qlsf, unsigned int *p_offset_scale1, unsigned int * p_offset_scale2,
                      short * p_no_scales );

static void lsf_mid_enc( Encoder_State *st, const float int_fs, const float qisp0[], const float qisp1[], float isp[], const short coder_type,
                         const short bwidth, const long core_brate, float Bin_Ener[], short ppp_mode, short nelp_mode );

static float vq_lvq_lsf_enc( short pred_flag, short mode, float u[], short * levels, short stages, float w[], short Idx[], const float * lsf,
                             const float * pred, unsigned int  p_offset_scale1[][MAX_NO_SCALES+1], unsigned int  p_offset_scale2[][MAX_NO_SCALES+1],
                             short  p_no_scales[][2], float *resq, float * lsfq );

static float qlsf_ARSN_tcvq_Enc_16k ( const float *x, float *y, short *indice, const float *w, const short nBits, short safety_net );

/*-------------------------------------------------------------------*
 * lsf_enc()
 *
 * Quantization of LSF vector
 *-------------------------------------------------------------------*/

void lsf_enc(
    Encoder_State *st,               /* i/o: state structure                             */
    const short L_frame,           /* i  : length of the frame                         */
    const short coder_type,        /* i  : coding type                                 */
    float *lsf_new,          /* o  : quantized LSF vector                        */
    float *lsp_new,          /* i/o: LSP vector to quantize/quantized            */
    float *lsp_mid,          /* i/o: mid-frame LSP vector                        */
    float *Aq,               /* o  : quantized A(z) for 4 subframes              */
    float *stab_fac,         /* o  : LSF stability factor                        */
    const short Nb_ACELP_frames
)
{
    short nBits = 0;
    float int_fs;
    short force_sf = 0;
    float fec_lsf[M], stab;
    short i;

    /* initialize */
    if( L_frame == L_FRAME )
    {
        int_fs = INT_FS_12k8;
    }
    else /* L_frame == L_FRAME16k */
    {
        int_fs = INT_FS_16k;
    }

    /* convert LSPs to LSFs */
    lsp2lsf( lsp_new, lsf_new, M, int_fs );


    /* check resonance for pitch clipping algorithm */
    gp_clip_test_lsf( st->core_brate, lsf_new, st->clip_var, 0 );

    /*-------------------------------------------------------------------------------------*
     * Find the number of bits for LSF quantization
     *-------------------------------------------------------------------------------------*/

    if ( st->core_brate == SID_2k40 )
    {
        nBits = LSF_BITS_CNG;
    }
    else
    {
        if ( st->nelp_mode == 0 && st->ppp_mode == 0 )
        {
            nBits = LSF_bits_tbl[LSF_BIT_ALLOC_IDX(st->core_brate, coder_type)];
        }
        else if ( st->nelp_mode == 1 )
        {
            if ( st->bwidth == NB )
            {
                nBits = 32; /* Stole 1 bit for SID/NELP harmonization*/

            }
            else if ( st->bwidth == WB )
            {
                nBits = 30; /* Stole 1 bit for SID/NELP harmonization */

            }
        }
        else if ( st->ppp_mode == 1 )
        {

            /*The LSF bit scaling does not work here. */
            nBits = 26; /*Stole 1 bit for SID PPP harmonization*/

        }
    }

    if ( Nb_ACELP_frames < 3 )
    {
        /* first three ACELP frames after an HQ frame shall be processed only with safety-net quantizer */
        force_sf = 1;
    }

    if ( st->next_force_safety_net == 1 )
    {
        /* in case of unstable filter, choose safety-net to help FEC */
        force_sf = 1;
        st->next_force_safety_net = 0;
    }

    /*-------------------------------------------------------------------------------------*
     * LSF quantization
     *-------------------------------------------------------------------------------------*/


    lsf_end_enc( st, lsf_new, lsf_new, st->mem_AR, st->mem_MA, nBits, coder_type, st->bwidth, st->Bin_E, int_fs, st->core_brate,
                 &st->streaklimit, &st->pstreaklen, force_sf, 0, 0, NULL, NULL, NULL, st->coder_type_raw );

    /* convert quantized LSFs back to LSPs */
    lsf2lsp( lsf_new, lsp_new, M, int_fs );

    if ( st->last_core == HQ_CORE && st->core == ACELP_CORE )
    {
        /* don't use old LSF values if this is the first ACELP frame after HQ frames */
        mvr2r( lsf_new, st->lsf_old, M );
    }

    if ( st->core_brate == SID_2k40 )
    {
        /* return if SID frame (conversion to A(z) done in the calling function) */
        return;
    }

    /*-------------------------------------------------------------------------------------*
     * FEC - enforce safety-net in the next frame in case of unstable filter
     *-------------------------------------------------------------------------------------*/

    if( st->last_L_frame != st->L_frame)
    {
        /* FEC - in case of core switching, use old LSFs that have been smoothed with adaptive mean */
        mvr2r( st->lsf_old, st->lsfoldbfi1, M );
        mvr2r( st->lsf_old, st->lsfoldbfi0, M );
        mvr2r( st->lsf_old, st->lsf_adaptive_mean, M );
    }

    FEC_lsf_estim_enc( st, st->L_frame, fec_lsf );

    /* FEC - calculate LSF stability */
    stab = lsf_stab( lsf_new, fec_lsf, 0, st->L_frame);

    if (st->L_frame == L_FRAME16k && stab < STAB_FAC_LIMIT && coder_type == GENERIC )
    {
        st->next_force_safety_net = 1;
    }
    else
    {
        if (  stab < STAB_FAC_LIMIT
                && (st->clas == VOICED_CLAS
                    || (st->clas < VOICED_CLAS && coder_type == AUDIO) ) )
        {
            st->next_force_safety_net = 1;
        }
    }

    /* FEC - update adaptive LSF mean vector */
    for (i=0; i<M; i++)
    {
        st->lsf_adaptive_mean[i] = (st->lsfoldbfi1[i] + st->lsfoldbfi0[i] + lsf_new[i]) / 3;
    }

    /* FEC - update LSF memories */
    mvr2r( st->lsfoldbfi0, st->lsfoldbfi1, M );
    mvr2r( lsf_new, st->lsfoldbfi0, M );

    /*-------------------------------------------------------------------------------------*
     * Mid-frame LSF encoding
     * LSP interpolation and conversion of LSPs to A(z)
     *-------------------------------------------------------------------------------------*/

    if(st->rate_switching_reset)
    {
        /*extrapolation in case of unstable LSF convert*/
        mvr2r( lsp_new, st->lsp_old, M );
        mvr2r( lsf_new, st->lsf_old, M );
    }

    /* Mid-frame LSF encoding */
    lsf_mid_enc( st, int_fs, st->lsp_old, lsp_new, lsp_mid, coder_type, st->bwidth,
                 st->core_brate, st->Bin_E_old, st->ppp_mode, st->nelp_mode );


    if( st->last_core == HQ_CORE && st->core == ACELP_CORE )
    {
        /* don't use old LSP/LSF values if this is the first ACELP frame after HQ frames */
        mvr2r( lsp_mid, st->lsp_old, M );
        lsp2lsf( lsp_mid, st->lsf_old, M, int_fs );
    }

    /* LSP interpolation and conversion of LSPs to A(z) */
    int_lsp4( L_frame, st->lsp_old, lsp_mid, lsp_new, Aq, M, 0 );

    /*------------------------------------------------------------------*
     * Check LSF stability (distance between old LSFs and current LSFs)
     *------------------------------------------------------------------*/

    *stab_fac = lsf_stab( lsf_new, st->lsf_old, 0, st->L_frame);

    return;
}


/*-------------------------------------------------------------------*
 * lsfq_CNG()
 *
 * LSF quantizer for SID frames (uses 29 bits, 4 for VQ, 25 for LVQ)
 * Note:
 * The sampling frequency of the LP-CNG frame can be determined by checking the value of the highest order LSF
 * coefficient (last coefficient of lsf). If the last LSF coefficient (lsf[M-1]) is larger than 6350
 * the decoded frame is WB2 with sampling rate of 16 kHz, otherwise it is sampled at 12.8kHz and contains
 * either NB or WB LSF data.
 *-------------------------------------------------------------------*/

static void lsfq_CNG(
    Encoder_State *st,             /* i/o: encoder state structure      */
    const float *lsf,
    const float *wghts,
    float *qlsf_out,
    unsigned int *p_offset_scale1,
    unsigned int *p_offset_scale2,
    short *p_no_scales
)
{
    short i, j, idx_cv, idx_lvq[3];
    float min_dist, dist, dd[M], ddq[M];
    const float *p_cb;
    short first_cb, last_cb;
    int idx_lead_cng[2], idx_scale_cng[2];

    float qlsf[M];
    idx_cv = 0;

    /* quantize first stage with 4 bits */
    if ( lsf[M-1] > WB_LIMIT_LSF ) /* 16kHz sampled LSF vector*/
    {
        p_cb = &CNG_SN1[0];
        first_cb = 0;
        last_cb = 6;
    }
    else /* 12.8kHz sampled LSF vector*/
    {
        p_cb = &CNG_SN1[6*M];
        first_cb = 6;
        last_cb = M;
    }

    min_dist = 1.0e30f;
    for ( i = first_cb; i < last_cb; i++ )
    {
        dist = 0.0f;
        for( j=0; j<M; j++ )
        {
            dist += wghts[j] * (*p_cb) * (*p_cb -2*lsf[j]);
            p_cb++;
        }

        if ( dist < min_dist )
        {
            min_dist = dist;
            idx_cv = i;
        }
    }

    /* calculate difference */
    for( i = 0; i < M; i++ )
    {
        dd[i] = lsf[i] - CNG_SN1[idx_cv*M+i];
    }

    /* quantize the difference with LVQ */
    mslvq_cng( idx_cv, dd, qlsf, ddq, idx_lead_cng, idx_scale_cng, wghts, p_no_scales );

    index_lvq( ddq, idx_lead_cng, idx_scale_cng, START_CNG + idx_cv, idx_lvq,
               p_offset_scale1, p_offset_scale2, p_no_scales );

    v_add( qlsf, &CNG_SN1[idx_cv*M], qlsf, M );
    /* write the VQ index to the bitstream */
    push_indice( st, IND_ISF_0_0, idx_cv, 4 );

    /* write the LVQ index to the bitstream */
    push_indice( st, IND_ISF_0_1, idx_lvq[0], LEN_INDICE );
    push_indice( st, IND_ISF_0_1, idx_lvq[1], LSF_BITS_CNG - 4 - LEN_INDICE );

    mvr2r(qlsf, qlsf_out,M);
    return;
}

/*-------------------------------------------------------------------*
 * qlsf_Mode_Select()
 *
 * Mode selection for LSF quantizer
 *-------------------------------------------------------------------*/

static short qlsf_Mode_Select(
    const float *lsf,             /* i : LSF vector */
    float *w,               /* i : weighting vector */
    float *pred1,           /* i : prediction vector */
    float streaklimit       /* i : predictive streak limit */
    ,float op_loop_thr		/* i  : Open-loop Threshold    */
)
{
    short i;
    short safety_net;
    float pred_pow2[M], En;

    /* calculate the prediction residual */
    for( i = 0; i < M; i++ )
    {
        pred_pow2[i] = lsf[i] - pred1[i];
    }

    /* calculate its energy */
    v_mult( pred_pow2, pred_pow2, pred_pow2, M );
    En = dotp( pred_pow2, w, M );

    /* choose the mode */
    if( En > streaklimit * op_loop_thr )
    {
        /* Safety-net */
        safety_net = 1;
    }
    else
    {
        /* Predictive */
        safety_net = 0;
    }

    return safety_net;
}



/*-------------------------------------------------------------------*
  * lsf_end_enc()
  *
  * Quantization of LSF parameters
  *-------------------------------------------------------------------*/
void lsf_end_enc(
    Encoder_State *st,           /* i/o: encoder state structure                                */
    const float *lsf,          /* i  : LSF in the frequency domain (0..6400)                  */
    float *qlsf,         /* o  : quantized LSF                                          */
    float *mem_AR,       /* i/o: quantizer memory for AR model                          */
    float *mem_MA,       /* i/o: quantizer memory for MA model                          */
    const short nBits_in,      /* i  : number of bits to spend on ISF quantization            */
    const short coder_type_org,/* i  : coding type                                            */
    const short bwidth,        /* i  : input signal bandwidth                                 */
    float *Bin_Ener,     /* i  : FFT Bin energy 128 *2 sets                             */
    const float int_fs,        /* i  : sampling frequency                                     */
    long  core_brate,
    float *streaklimit,  /* i/o: Multiplier to prefer safety-net, more and more in long streaks of predictive frames*/
    short *pstreaklen,   /* i/o: Count of consecutive predictive frames */
    short force_sf,      /* i  : Force safety-net usage (if possible), due to filter instability, MDCT-core switching etc. */
    short rf_flag,
    short mode2_flag,
    int * lpc_param,
    short * no_indices,
    short * bits_param_lpc,
    short coder_type_raw
)
{
    short i;
    short Idx0[MAX_VQ_STAGES+3];      /* Optimal codebook indices for safety-net quantizer                 */
    short Idx1[MAX_VQ_STAGES+3];      /* Optimal codebook indices for predictive quantizer                 */
    short indice[MAX_VQ_STAGES+3];    /* Temp. array of indice for vector de-quantizer                     */
    short mode_lvq = 0, mode_lvq_p = 0;
    short bits0[MAX_VQ_STAGES], bits1[MAX_VQ_STAGES];
    const short *Bit_alloc1 = NULL;
    float Err[2];                     /* Quantization error for safety-net(0) and predictive(1) quantizers */
    float Tmp [M];                    /* Temporary target vector (mean and prediction removed)             */
    float pred0[M];                   /* Prediction for the safety-net quantizer (usually mean)            */
    float pred1[M];                   /* Prediction for the predictive quantizer                           */
    float pred2[M];                   /* Prediction for the predictive quantizer                           */
    float wghts[M];                   /* Weighting used for quantizer (currently GSM based)                */
    short stages0;                    /* Amount of stages used by safety-net quantizer                     */
    short stages1;                    /* Amount of stages used by predictive quantizer                     */
    short levels0[MAX_VQ_STAGES];     /* Sizes of different codebook stages for safety-net quantizer       */
    short levels1[MAX_VQ_STAGES];     /* Sizes of different codebook stages for predictive quantizer       */
    short predmode;                   /* 0: safety-net only, 1: predictive only, 2: best of the two        */
    short safety_net, cumleft, num_bits;
    short *Idx, stages, *bits;
    float Tmp1[M], Tmp2[M];           /* Temporary target vectors for MA and AR quantizers respectively     */
    float abs_threshold;              /* Abslute Error value that is considered as "good enough" (in practice close to SD of 1.0dB)*/
    float lsfq[M*2], resq[M*2];

    short coder_type;
    short nBits;
    short TCQIdx0[M+2];               /* Optimal codebook indices for VQ-TCQ quantizer                 */
    short *TCQIdx;


    nBits = nBits_in;

    if((coder_type_org == GENERIC) && (int_fs == INT_FS_16k) && (rf_flag == 0) && (mode2_flag == 0))
    {

        if (coder_type_raw == VOICED)
        {
            coder_type = VOICED; /* Reflect Inactive mode */
        }
        else
        {
            nBits--; /* This is for real Generic*/
            coder_type = coder_type_org;
        }
    }
    else
    {
        coder_type = coder_type_org;
    }
    /*--------------------------------------------------------------------------------*
    * Calculate the number of stages and levels for each stage based on allowed bit budget
    * Set absolute threshold for codebook-type decision logic
    *--------------------------------------------------------------------------------*/
    if ( bwidth == NB )
    {
        abs_threshold = SFNETLOWLIMIT_NB;
    }
    else
    {
        abs_threshold = SFNETLOWLIMIT_WB;
    }
    Unified_weighting(&Bin_Ener[L_FFT/2], lsf, wghts, bwidth==NB, coder_type==UNVOICED, (int)(int_fs),M);

    /*--------------------------------------------------------------------------------*
    * LSF quantization of SID frames
    *--------------------------------------------------------------------------------*/

    if ( core_brate == SID_2k40 )
    {
        lsfq_CNG( st, lsf,  wghts, qlsf, &st->offset_scale1[0][0], &st->offset_scale2[0][0], &st->no_scales[0][0] );
        v_sort( qlsf, 0, M-1 );
        reorder_lsf( qlsf, MODE1_LSF_GAP, M, int_fs );

        return;
    }

    predmode = find_pred_mode(coder_type, bwidth, int_fs, &mode_lvq, &mode_lvq_p,st->total_brate);

    /*----------------------------------------------------------------*
    * Calculate number of stages and levels for each stage based on the allowed bit allocation
    * (subtract one bit for LSF predictor selection)
    *----------------------------------------------------------------*/

    lsf_allocate( nBits-(predmode>>1), mode_lvq, mode_lvq_p, &stages0, &stages1, levels0, levels1, bits0, bits1 );

    /*--------------------------------------------------------------------------------*
    * LSF quantization of all other frames but SID frames
    * Select safety-net or predictive mode
    *--------------------------------------------------------------------------------*/

    Err[0] = FLT_MAX;
    Err[1] = FLT_MAX;
    /* for mem_MA update */
    for (i=0; i<M; i++)
    {
        pred1[i] = ModeMeans[mode_lvq][i]+MU_MA*mem_MA[i];
    }

    if ( predmode == 0 ) /* Safety-net only */
    {
        /* Subtract only mean */
        mvr2r(ModeMeans[mode_lvq], pred0,M);
        v_sub(lsf, pred0, Tmp, M);

        /* LVQ quantization (safety-net only) */
        Err[0] = vq_lvq_lsf_enc(0, mode_lvq, Tmp, levels0, stages0,wghts, Idx0, lsf, pred0,
                                st->offset_scale1,st->offset_scale2, st->no_scales, resq, lsfq);
        safety_net = 1;
        *pstreaklen = 0; /* Streak is ended with safety-net */
    }
    else if (predmode ==1) /* only MA prediction */
    {
        v_sub(lsf, pred1, Tmp1, M);
        Err[1] = vq_lvq_lsf_enc(2, mode_lvq_p, Tmp1, levels1, stages1, wghts, Idx1, lsf, pred1,
                                st->offset_scale1_p,st->offset_scale2_p,st->no_scales_p,resq, lsfq );

        safety_net = 0;
    }
    else /* Switched Safety-Net/AR prediction */
    {
        /* Subtract mean and AR prediction */
        mvr2r(ModeMeans[mode_lvq], pred0, M);
        /* subtract only mean */
        v_sub(lsf, pred0, Tmp, M);
        for (i = 0; i < M; i++)
        {
            /* subtract mean and AR prediction */
            pred2[i] = pred0[i] + Predictors[mode_lvq_p][i]*(mem_AR[i]-pred0[i]);
            Tmp2[i] = lsf[i] - pred2[i];
        }
        /* Adaptive scaling factor (multiplier) is updated in order to reduce the amount of consecutive predictive frames in
           case of possible frame erasure. AR-predictive usage for VOICED mode is allowed to be higher than other modes. */
        if ( ((*pstreaklen > (STREAKLEN+3))&&(coder_type==VOICED)) || ((*pstreaklen > (STREAKLEN)) &&(coder_type!=VOICED)))
        {
            /* update the adaptive scaling factor to become smaller with increasing number of concecutive predictive frames. */
            *streaklimit *= STREAKMULT;
        }

        if ( *pstreaklen == 0 )
        {
            /* reset the adaptive scaling factor */
            *streaklimit = 1.0f;
        }


        /* VOICED_WB@16kHz */
        if ( int_fs == INT_FS_16k && coder_type == VOICED )
        {
            /* select safety_net or predictive in open loop*/
            safety_net = qlsf_Mode_Select( lsf, wghts, pred2, *streaklimit, OP_LOOP_THR_HVO );

            if ( force_sf == 1 )
            {
                safety_net = 1;
            }

            if ( safety_net )
            {
                /* Safety-net - BC-TCQ quantization : SN */
                Err[0] = qlsf_ARSN_tcvq_Enc_16k( Tmp, lsfq, TCQIdx0, wghts, nBits-1, safety_net );
                *pstreaklen = 0;
            }
            else
            {
                /* predictive - BC-TCQ quantization : AR */
                /* For consistency Err[1] contains predictive error */
                Err[1] = qlsf_ARSN_tcvq_Enc_16k( Tmp2, lsfq, TCQIdx0, wghts, nBits-1, safety_net );
                (*pstreaklen)++;
            }
        }
        /* all other frames (not VOICED@16kHz) */
        else
        {
            /* safety-net */
            Err[0] = vq_lvq_lsf_enc(0, mode_lvq, Tmp, levels0, stages0, wghts, Idx0, lsf, pred0,
                                    st->offset_scale1,st->offset_scale2,st->no_scales, resq, lsfq);
            /* Predictive quantizer is calculated only if it can be selected, this saves computation */
            if (!force_sf || Err[0] > abs_threshold)
            {
                Err[1] = vq_lvq_lsf_enc(2, mode_lvq_p, Tmp2, levels1, stages1, wghts, Idx1, lsf, pred2,
                                        st->offset_scale1_p, st->offset_scale2_p, st->no_scales_p, &resq[M], &lsfq[M]);
            }
            /* Select whether to use safety-net or predictive LSF quantizer. The decision is based on following:
               if the non-predictive (safety-net) quantization error (Err[0]) is low enough it is selected
               or if the predictively quantized error (Err[1]) is by at least adaptive margin smaller than non-predictive quantizer.
               or if the in case of frame erasure the resulting concealed predictive LSF would be unstable safety-net is selected */
            if ( force_sf || Err[0]*(*streaklimit) < PREFERSFNET * Err[1] || Err[0] < abs_threshold )
            {
                safety_net = 1;
                *pstreaklen = 0; /* Reset the consecutive predictive frame counter */
            }
            else
            {
                safety_net = 0;
                (*pstreaklen)++; /* Increase the consecutive predictive frame counter by one */
            }
        }
    }

    /*--------------------------------------------------------------------------*
     * Write indices to array
     *--------------------------------------------------------------------------*/

    if (mode2_flag == 0)
    {
        /* write coder_type bit for VOICED@16kHz or GENERIC@16kHz */
        if((coder_type_org == GENERIC)
                &&(int_fs == INT_FS_16k)
          )
        {
            /* VOICED =2 and GENERIC=3, so "coder_type-2" means VOICED =0 and GENERIC=1*/
            push_indice( st, IND_LSF_PREDICTOR_SELECT_BIT, coder_type-2, 1 );
        }
        /* write predictor selection bit */
        if ( predmode == 2 )
        {
            push_indice( st, IND_LSF_PREDICTOR_SELECT_BIT, safety_net, 1 );
        }

        if ( (coder_type == VOICED) && (int_fs == INT_FS_16k) )
        {
            /* BC-TCVQ (only for VOICED@16kHz) */
            TCQIdx = &TCQIdx0[1];
            Bit_alloc1 = &BC_TCVQ_BIT_ALLOC_40B[1];
            for( i=0; i<M/2+3; i++ )
            {
                push_indice( st, IND_LSF, TCQIdx[i], Bit_alloc1[i]);
            }
        }
        else
        {
            cumleft = nBits;
            if( predmode == 2 )
            {
                /* subtract predictor selection bit */
                cumleft = nBits - 1;
            }

            if ( safety_net )
            {
                stages = stages0;
                Idx = Idx0;
                bits = bits0;
            }
            else
            {
                stages = stages1;
                Idx = Idx1;
                bits = bits1;
            }

            for ( i=0; i<stages-1; i++ )
            {
                indice[i] = Idx[i];
                num_bits = bits[i];
                cumleft -= num_bits;

                push_indice( st, IND_LSF, indice[i], num_bits );
            }

            while ( cumleft > 0 )
            {
                indice[i] = Idx[i];

                if ( cumleft > LEN_INDICE )
                {
                    num_bits = LEN_INDICE;
                }
                else
                {
                    num_bits = cumleft;
                }

                cumleft -= num_bits;

                push_indice( st, IND_LSF, indice[i], num_bits );

                i++;
            }
        }
    }
    else
    {
        if ( (coder_type == VOICED) && (int_fs == INT_FS_16k) )
        {
            /* BC-TCVQ (only for VOICED@16kHz) */
            /* Number of quantization indices */
            *no_indices = 10;
            for(i=0; i<*no_indices; i++)
            {
                lpc_param[i] = (int)TCQIdx0[i];
                bits_param_lpc[i] = BC_TCVQ_BIT_ALLOC_40B[i];
            }
        }
        else
        {
            /* Number of quantization indices */

            /* there are 31 bits */
            if (safety_net==1)
            {
                Idx = Idx0;
                *no_indices = stages0+1;
                for(i=0; i<stages0; i++)
                {
                    lpc_param[i] = (int)(Idx[i]);
                    indice[i] = Idx[i];
                    bits_param_lpc[i] = bits0[i];
                }
                lpc_param[stages0] = (int)(Idx[stages0]);
                indice[stages0] = Idx[stages0];

                bits_param_lpc[stages0-1] = LEN_INDICE;
                bits_param_lpc[stages0] = bits0[stages0-1]-LEN_INDICE;

            }
            else
            {
                *no_indices = stages1+1;
                Idx = Idx1;
                for(i=0; i<stages1; i++)
                {
                    lpc_param[i] = (int)(Idx[i]);
                    indice[i] = Idx[i];
                    bits_param_lpc[i] = bits1[i];
                }
                lpc_param[stages1] = (int)(Idx[stages1]);
                indice[stages1] = Idx[stages1];

                bits_param_lpc[stages1-1] = LEN_INDICE;
                bits_param_lpc[stages1] = bits1[stages1-1]-LEN_INDICE;
            }
            if (predmode==2)
            {
                for (i=*no_indices; i>0; i--)
                {
                    lpc_param[i] = lpc_param[i-1];
                    bits_param_lpc[i] = bits_param_lpc[i-1];
                }
                lpc_param[0] = safety_net; /* put the safety net info on the last param */
                bits_param_lpc[0] = 1;
                *no_indices = *no_indices+1;
            }
        }
    }


    /*--------------------------------------------------------------------------*
     *  De-quantize encoded LSF vector
     *--------------------------------------------------------------------------*/

    if ( safety_net )
    {
        /* Safety-net */
        if ( coder_type == VOICED && int_fs == INT_FS_16k )
        {
            /* BC-TCQ */
            mvr2r( lsfq, mem_MA, M );
            v_add( lsfq, pred0, qlsf, M );
        }
        else
        {
            vq_dec_lvq( 1, qlsf, &indice[0], stages0, M, mode_lvq, levels0[stages0-1],
                        &st->offset_scale1[0][0], &st->offset_scale2[0][0], &st->offset_scale1_p[0][0], &st->offset_scale2_p[0][0],
                        &st->no_scales[0][0], &st->no_scales_p[0][0] );

            v_add(qlsf, pred0, qlsf, M );
            v_sub(qlsf, pred1,mem_MA, M);
        }
    }
    else
    {
        if ( (coder_type == VOICED) && (int_fs == INT_FS_16k) )
        {
            /* BC-TCVQ */
            mvr2r( lsfq, mem_MA, M );
            v_add( lsfq, pred2, qlsf, M );
        }
        else
        {
            /* LVQ */
            vq_dec_lvq( 0, qlsf, &indice[0], stages1, M, mode_lvq_p, levels1[stages1-1],
                        &st->offset_scale1[0][0], &st->offset_scale2[0][0], &st->offset_scale1_p[0][0],
                        &st->offset_scale2_p[0][0], &st->no_scales[0][0], &st->no_scales_p[0][0] );
            if (predmode==1)
            {
                mvr2r(qlsf, mem_MA, M);
                v_add( qlsf, pred1, qlsf, M );
            }
            else
            {
                v_add( qlsf, pred2, qlsf, M );
                v_sub(qlsf, pred1, mem_MA, M);
            }
        }
    }

    /*--------------------------------------------------------------------------*
     * Sort the quantized vector
     * Verify stability
     * Update AR-predictor memory
     *--------------------------------------------------------------------------*/

    /* Sort the quantized vector to ascending order */
    v_sort( qlsf, 0, M-1 );

    /* Verify stability by adding minimum separation */
    reorder_lsf( qlsf, MODE1_LSF_GAP, M, int_fs );

    /* Update predictor memories */
    mvr2r( qlsf, mem_AR, M );


    return;
}



/*-------------------------------------------------------------------*
 * first_VQstages()
 *
 *
 *-------------------------------------------------------------------*/

static void first_VQstages(
    const float *  const *cb,
    float u[],              /* i  : vector to be encoded (prediction and mean removed)  */
    short *levels,          /* i  : number of levels in each stage                      */
    short stagesVQ,         /* i  : number of stages                                    */
    float w[],              /* i  : weights                                             */
    short N,                /* i  : vector dimension                                    */
    short max_inner,        /* i  : maximum number of swaps in inner loop               */
    short indices_VQstage[]
)
{
    float ftmp, resid_buf[2*LSFMBEST*M], *resid[2], dist_buf[2*LSFMBEST], *dist[2];
    float en, tmp, Tmp[M], *pTmp;
    short *pTmp_short, idx_buf[2*LSFMBEST*MAX_VQ_STAGES], parents[LSFMBEST], counter=0, j,
                                                                             m, s,c, c2, p_max, *indices[2];
    short maxC = LSFMBEST;

    /*float dd[16];*/
    const float *cb_stage, *cbp;

    /* Set pointers to previous (parent) and current node (parent node is indexed [0], current node is indexed [1]) */
    indices[0] = idx_buf;
    indices[1] = idx_buf + maxC*stagesVQ;
    resid[0] = resid_buf;
    resid[1] = resid_buf + maxC*N;
    dist[0] = dist_buf;
    dist[1] = dist_buf + maxC;

    set_s( idx_buf, 0, (const short)(2*stagesVQ*maxC) );
    set_s( parents, 0, maxC ) ;

    /* Set up inital distance vector */
    for( tmp = 0, j=0; j<N; j++ )
    {
        tmp += u[j]*u[j]*w[j] ;
    }
    set_f( dist[1], tmp, maxC ) ;

    /* Set up inital error (residual) vectors */
    pTmp = resid[1];
    for( c = 0; c < maxC; c++ )
    {
        mvr2r( u, pTmp, N );
        pTmp += N;
    }

    /*----------------------------------------------------------------*
     * LSF quantization
     *----------------------------------------------------------------*/

    /* Loop over all stages */
    for ( m = 1, s = 0; s < stagesVQ; s++ )
    {
        /* set codebook pointer to point to first stage */
        cbp = cb[s];

        /* save pointer to the beginning of the current stage */
        cb_stage = cbp;

        /* swap pointers to parent and current nodes */
        pTmp_short = indices[0];
        indices[0] = indices[1];
        indices[1] = pTmp_short;

        pTmp = resid[0];
        resid[0] = resid[1];
        resid[1] = pTmp;

        pTmp = dist[0];
        dist[0] = dist[1];
        dist[1] = pTmp;

        /* p_max points to maximum distortion node (worst of best) */
        p_max = 0;

        /* set distortions to a large value */
        set_f( dist[1], 99e10f, maxC );

        for ( j = 0; j < levels[s]; j++ )
        {
            /* compute weighted codebook element and its energy */
            for ( c2 = 0; c2 < N; c2++ )
            {
                Tmp[c2] = w[c2] * cbp[c2];
            }

            en = cbp[0] * Tmp[0];
            for ( c2 = 1; c2 < N; c2++ )
            {
                en += cbp[c2] * Tmp[c2];
            }
            cbp += N ;

            /* iterate over all parent nodes */
            for( c = 0; c < m; c++ )
            {
                pTmp = &resid[0][c*N];
                tmp = pTmp[0] * Tmp[0];
                for ( c2=1; c2<N; c2++ )
                {
                    tmp += pTmp[c2] * Tmp[c2];
                }
                tmp = en - 2.0f * tmp;
                tmp += dist[0][c];

                if ( tmp <= dist[1][p_max] )
                {
                    /* replace worst */
                    dist[1][p_max] = tmp;
                    indices[1][p_max*stagesVQ+s] = j;
                    parents[p_max] = c;

                    /* limit number of times inner loop is entered */
                    if ( counter < max_inner )
                    {
                        counter++;
                        if ( counter < max_inner )
                        {
                            /* find new worst */
                            p_max = maximum(dist[1], maxC, &ftmp);
                        }
                        else
                        {
                            /* find minimum distortion */
                            p_max = minimum(dist[1], maxC, &ftmp);
                        }
                    }
                }
            }
        }

        /*------------------------------------------------------------*
         * Compute error vectors for each node
         *------------------------------------------------------------*/

        for ( c = 0; c < maxC; c++ )
        {
            /* subtract codebook entry from the residual vector of the parent node */
            pTmp = resid[1]+c*N ;
            mvr2r( resid[0]+parents[c]*N, pTmp, N );
            v_sub( pTmp, cb_stage+(indices[1][c*stagesVQ+s])*N, pTmp, N );

            /* get indices that were used for parent node */
            mvs2s( indices[0]+parents[c]*stagesVQ, indices[1]+c*stagesVQ, s );
        }

        m = maxC;
    }

    mvs2s(indices[1],indices_VQstage,maxC*stagesVQ );

    return;
}

/*---------------------------------------------------------------------------
 * vq_enc_lsf_lvq()
 *
 *  Multi-stage VQ encoder for LSF quantization. Trained codebooks are used in initial stages
 *  and lattice-VQ quantization is applied on residual vector in other stages.
 *
 * Note:
 *    Compared to normal multistage VQ resulting LSF vector is reordered to ascending order before
 *    weighted error calculation (spectral distortion) at the final stage.
 *
 * Returns:
 *    Weighted error
 *--------------------------------------------------------------------------*/

static float vq_lvq_lsf_enc(
    short pred_flag,
    short mode,
    float u[],
    short * levels,
    short stages,
    float w[],
    short Idx[],
    const float * lsf,
    const float * pred,
    unsigned int  p_offset_scale1[][MAX_NO_SCALES+1],
    unsigned int  p_offset_scale2[][MAX_NO_SCALES+1],
    short  p_no_scales[][2],
    float *resq,
    float * lsfq
)
{
    int i;
    const float *const *cb, *cb_stage;
    float cand[LSFMBEST][M];
    int maxC=LSFMBEST, stagesVQ;
    short  mode_glb, j, indices_firstVQ[LSFMBEST*MAX_VQ_STAGES], c2;
    float dd[M];
    float quant[LSFMBEST][M], diff[M];
    float lat_cv[LSFMBEST][M], e[LSFMBEST], ftmp, tmp;
    int idx_lead[LSFMBEST][2], idx_scale[LSFMBEST][2];

    stagesVQ = stages-1;
    /* Codebook selection */
    if (pred_flag==0) /* safety net*/
    {
        cb = &Quantizers[CB[mode]];
        mode_glb = offset_lvq_modes_SN[mode] + offset_in_lvq_mode_SN[mode][levels[stagesVQ] - min_lat_bits_SN[mode]];
    }
    else /*  pred */
    {
        cb = &Quantizers_p[CB_p[mode]];
        mode_glb = offset_lvq_modes_pred[mode] + offset_in_lvq_mode_pred[mode][levels[stagesVQ] - min_lat_bits_pred[mode]];
    }
    if (stagesVQ>0)
    {
        /* first VQ stages */
        first_VQstages( cb, u, levels, stagesVQ, w, M, MSVQ_MAXCNT, indices_firstVQ );
    }


    for ( i=0; i<maxC; i++ )
    {
        mvr2r( pred, cand[i], M );
        for ( j=0; j<stagesVQ; j++ )
        {
            Idx[j] = indices_firstVQ[i*stagesVQ+j];
        }

        for ( j=0; j<stagesVQ; j++ )
        {
            cb_stage = cb[j];
            v_add( cand[i], cb_stage+Idx[j]*M, cand[i], M );
        }


        /* LVQ quantization */
        v_sub( lsf, cand[i], dd, M );
        mslvq(dd, quant[i], lat_cv[i], idx_lead[i], idx_scale[i], w, mode, mode_glb, pred_flag, p_no_scales);
        v_add( cand[i], quant[i], cand[i], M );

        /* arrange the LSF candidate vector prior to selection to an ascending order*/
        v_sort(cand[i],0,M-1);

        /* calculate the spectral distortion using weighted MSE of sorted LSF vector*/
        v_sub( cand[i], lsf, diff, M );
        v_mult( diff, diff, diff, M );
        e[i] = dotp( diff, w, M );
    }

    /* find the optimal candidate */
    c2 = minimum( e, maxC, &ftmp );
    set_f(resq,0.0f, M);
    for ( j=0; j<stagesVQ; j++ )
    {
        Idx[j] = indices_firstVQ[c2*stagesVQ+j];
        cb_stage = cb[j];
        v_add(resq,cb_stage+Idx[j]*M, resq, M);
    }
    v_add(resq, quant[c2],resq, M); /* quantized prediction residual */
    /*for(j=0;j<M;j++)
      {
      resq_ind[i] = LSFM(resq[i]);
      }
    */
    mvr2r(cand[c2], lsfq, M);
    index_lvq( lat_cv[c2], idx_lead[c2], idx_scale[c2], mode_glb, &Idx[stagesVQ],
               &p_offset_scale1[0][0], &p_offset_scale2[0][0], &p_no_scales[0][0] );

    tmp = e[c2];

    return tmp;
}



/*---------------------------------------------------------------------------
 * BcTcvq_1st()
 *
 *
 *--------------------------------------------------------------------------*/

static void BcTcvq_1st(
    float x[][2],
    const float CB[][128][2],
    short s[][16],
    short c[][16],
    float cDist[][16],
    float Q[][16][2],
    float W[][2]
)
{
    short state, prev_state;
    short index, bestCode;
    float dist, minDist;

    for(state = 0; state < NUM_STATE; state +=2)
    {
        prev_state = NTRANS[0][state];
        index	   = NTRANS[2][state];

        minDist  = (x[0][0] - CB[0][index][0]) * (x[0][0] - CB[0][index][0]) * W[0][0];
        minDist += (x[0][1] - CB[0][index][1]) * (x[0][1] - CB[0][index][1]) * W[0][1];
        bestCode = index;

        for(index = index+8; index < 128; index += 8)
        {
            dist  = (x[0][0] - CB[0][index][0]) * (x[0][0] - CB[0][index][0]) * W[0][0];
            dist += (x[0][1] - CB[0][index][1]) * (x[0][1] - CB[0][index][1]) * W[0][1];

            if(dist < minDist)
            {
                minDist  = dist;
                bestCode = index;
            }
        }

        /* Update */
        s[0][state]		= prev_state;
        c[0][state]		= bestCode;
        cDist[0][state] = minDist;
        Q[0][state][0]  = CB[0][bestCode][0];
        Q[0][state][1]  = CB[0][bestCode][1];
    }

    return;
}

/*---------------------------------------------------------------------------
 * BcTcvq_2nd()
 *
 *
 *--------------------------------------------------------------------------*/

static void BcTcvq_2nd(
    float x[][2],
    const float CB[][128][2],
    short s[][16],
    short c[][16],
    float cDist[][16],
    float Q[][16][2],
    float W[][2],
    const float itc[][2][2]
)
{
    short state, prev_state;
    short index, bestCode;
    float dist, minDist;
    float pred[N_DIM], target[N_DIM];

    for(state = 0; state < NUM_STATE; state++)
    {
        prev_state = NTRANS[0][state];
        index	   = NTRANS[2][state];

        /* Prediction */
        pred[0]   = itc[0][0][0] * Q[0][prev_state][0] + itc[0][0][1] * Q[0][prev_state][1];
        pred[1]   = itc[0][1][0] * Q[0][prev_state][0] + itc[0][1][1] * Q[0][prev_state][1];
        target[0] = x[1][0] - pred[0];
        target[1] = x[1][1] - pred[1];

        minDist  = (target[0] - CB[1][index][0]) * (target[0] - CB[1][index][0]) * W[1][0];
        minDist += (target[1] - CB[1][index][1]) * (target[1] - CB[1][index][1]) * W[1][1];
        bestCode = index;

        for(index = index + 8; index < 128; index += 8)
        {
            dist  = (target[0] - CB[1][index][0]) * (target[0] - CB[1][index][0]) * W[1][0];
            dist += (target[1] - CB[1][index][1]) * (target[1] - CB[1][index][1]) * W[1][1];

            if(dist < minDist)
            {
                minDist  = dist;
                bestCode = index;
            }
        }

        /* Update */
        s[1][state]		= prev_state;
        c[1][state]		= bestCode;
        cDist[1][state] = cDist[0][prev_state] + minDist;
        Q[1][state][0]  = CB[1][bestCode][0] + pred[0];
        Q[1][state][1]  = CB[1][bestCode][1] + pred[1];
    }

    return;
}


/*---------------------------------------------------------------------------
 * BcTcvq_SubBlock()
 *
 *
 *--------------------------------------------------------------------------*/

static void BcTcvq_SubBlock(
    float x[][2],
    const float CB[][64][2],
    short s[][16],
    short c[][16],
    float cDist[][16],
    float Q[][16][2],
    short stage,
    float W[][2],
    const float itc[][2][2]
)
{
    short stage1, stage2, state, prev_state, branch;
    short index, bestCode, brCode[N_DIM];
    float dist, minDist, brDist[N_DIM];
    float pred[N_DIM], target[N_DIM], brQuant[N_DIM][N_DIM];

    stage1 = stage-1;
    stage2 = stage-2;

    for(state = 0; state < NUM_STATE; state ++)
    {

        /* 1st brarnch search */
        prev_state = NTRANS[0][state];
        index	   = NTRANS[2][state];

        /* Prediction */
        pred[0]   = itc[stage1][0][0] * Q[stage1][prev_state][0] + itc[stage1][0][1] * Q[stage1][prev_state][1];
        pred[1]   = itc[stage1][1][0] * Q[stage1][prev_state][0] + itc[stage1][1][1] * Q[stage1][prev_state][1];
        target[0] = x[stage][0] - pred[0];
        target[1] = x[stage][1] - pred[1];

        minDist  = (target[0] - CB[stage2][index][0]) * (target[0] - CB[stage2][index][0]) * W[stage][0];
        minDist += (target[1] - CB[stage2][index][1]) * (target[1] - CB[stage2][index][1]) * W[stage][1];
        bestCode = index;

        for(index = index + 8; index < 64; index += 8)
        {
            dist  = (target[0] - CB[stage2][index][0]) * (target[0] - CB[stage2][index][0]) * W[stage][0];
            dist += (target[1] - CB[stage2][index][1]) * (target[1] - CB[stage2][index][1]) * W[stage][1];

            if(dist < minDist)
            {
                minDist  = dist;
                bestCode = index;
            }
        }

        brCode[0]	  = bestCode;
        brDist[0]	  = cDist[stage1][prev_state] + minDist;
        brQuant[0][0] = CB[stage2][bestCode][0] + pred[0];
        brQuant[0][1] = CB[stage2][bestCode][1] + pred[1];

        /* 2nd branch search */
        prev_state = NTRANS[1][state];
        index	   = NTRANS[3][state];

        /* Prediction */
        pred[0]   = itc[stage1][0][0] * Q[stage1][prev_state][0] + itc[stage1][0][1] * Q[stage1][prev_state][1];
        pred[1]   = itc[stage1][1][0] * Q[stage1][prev_state][0] + itc[stage1][1][1] * Q[stage1][prev_state][1];
        target[0] = x[stage][0] - pred[0];
        target[1] = x[stage][1] - pred[1];

        minDist  = (target[0] - CB[stage2][index][0]) * (target[0] - CB[stage2][index][0]) * W[stage][0];
        minDist += (target[1] - CB[stage2][index][1]) * (target[1] - CB[stage2][index][1]) * W[stage][1];
        bestCode = index;

        for(index = index + 8; index < 64; index += 8)
        {

            dist  = (target[0] - CB[stage2][index][0]) * (target[0] - CB[stage2][index][0]) * W[stage][0];
            dist += (target[1] - CB[stage2][index][1]) * (target[1] - CB[stage2][index][1]) * W[stage][1];

            if(dist < minDist)
            {
                minDist  = dist;
                bestCode = index;
            }
        }

        brCode[1]	  = bestCode;
        brDist[1]	  = cDist[stage1][prev_state] + minDist;
        brQuant[1][0] = CB[stage2][bestCode][0] + pred[0];
        brQuant[1][1] = CB[stage2][bestCode][1] + pred[1];

        /* Select Best branch */
        branch = 1;

        if(brDist[0] <= brDist[1])
        {
            branch = 0;
        }

        /* Update */
        s[stage][state]		= NTRANS[branch][state];
        c[stage][state]		= brCode[branch];
        cDist[stage][state] = brDist[branch];
        Q[stage][state][0]  = brQuant[branch][0];
        Q[stage][state][1]  = brQuant[branch][1];
    }

    return;
}


/*---------------------------------------------------------------------------
 * BcTcvq_FixSearch()
 *
 *
 *--------------------------------------------------------------------------*/

static float BcTcvq_FixSearch(
    float x[][2],
    const float CB[][32][2],
    short c[][4],
    float Q[][16][2],
    const short FixBranch[][4][4],
    short stage,
    short inis,
    short fins,
    short *prev_state,
    float W[][2],
    const float itc[][2][2]
)
{
    short stage1, stage4, branch;
    short index, bestCode;
    float dist, minDist;
    float pred[N_DIM], target[N_DIM];

    stage1 = stage-1;
    stage4 = stage-4;

    branch = FixBranch[inis>>2][fins][stage4];
    index  = NTRANS2[branch+2][*prev_state];

    /* Prediction */
    pred[0]   = itc[stage1][0][0] * Q[stage1][*prev_state][0] + itc[stage1][0][1] * Q[stage1][*prev_state][1];
    pred[1]   = itc[stage1][1][0] * Q[stage1][*prev_state][0] + itc[stage1][1][1] * Q[stage1][*prev_state][1];
    target[0] = x[stage][0] - pred[0];
    target[1] = x[stage][1] - pred[1];

    minDist  = (target[0] - CB[stage4][index][0]) * (target[0] - CB[stage4][index][0]) * W[stage][0];
    minDist += (target[1] - CB[stage4][index][1]) * (target[1] - CB[stage4][index][1]) * W[stage][1];
    bestCode = index;

    for(index = index + 8; index < 32; index += 8)
    {
        dist  = (target[0] - CB[stage4][index][0]) * (target[0] - CB[stage4][index][0]) * W[stage][0];
        dist += (target[1] - CB[stage4][index][1]) * (target[1] - CB[stage4][index][1]) * W[stage][1];

        if(dist < minDist)
        {
            minDist  = dist;
            bestCode = index;
        }
    }

    /* Update */
    *prev_state				 = NTRANS2[branch][*prev_state];
    c[fins][stage4]			 = bestCode;
    Q[stage][*prev_state][0] = CB[stage4][bestCode][0] + pred[0];
    Q[stage][*prev_state][1] = CB[stage4][bestCode][1] + pred[1];

    return minDist;

}

/*---------------------------------------------------------------------------
 * optimalPath()
 *
 *
 *--------------------------------------------------------------------------*/

static short optimalPath(
    float cDist[][16],
    float blockDist[],
    short blockCodeword[][4],
    short bestCodeword[],
    short codeWord[][16],
    short bestState[],
    short preState[][16]
)
{
    short stage, state;
    float opDist[NUM_STATE];
    float minDist;
    short fBlock;
    short prev_state;

    for(state = 0; state < NUM_STATE; state++)
    {
        opDist[state] = cDist[3][state] + blockDist[state];
    }

    minDist = opDist[0];
    fBlock  = 0;

    for(state = 1; state < NUM_STATE; state++)
    {
        if(opDist[state] < minDist)
        {
            minDist = opDist[state];
            fBlock  = state;
        }
    }

    prev_state = bestState[4] = fBlock;

    for(stage = N_STAGE_VQ - 5; stage >= 0; stage --)
    {
        bestCodeword[stage] = codeWord[stage][prev_state];
        bestState[stage]	= preState[stage][prev_state];
        prev_state		    = bestState[stage];
    }

    for(stage = 0; stage < 4; stage ++)
    {
        bestCodeword[stage + 4] = blockCodeword[fBlock][stage];
    }

    return fBlock;
}

/*---------------------------------------------------------------------------
 * quantEnc()
 *
 *
 *--------------------------------------------------------------------------*/

static void quantEnc(
    float *y,
    short c[],
    const float CB_SUB1[][128][2],
    const float CB_SUB2[][64][2],
    const float CB_SUB3[][32][2],
    const float itc[][2][2]
)
{
    short i,j;
    short stage;
    float pred[N_DIM], Y[8][2];

    /* stage #1 */
    Y[0][0]  = CB_SUB1[0][c[0]][0];
    Y[0][1]  = CB_SUB1[0][c[0]][1];

    /* stage #2 */
    pred[0]		= itc[0][0][0] * Y[0][0] + itc[0][0][1]*Y[0][1];
    pred[1]		= itc[0][1][0] * Y[0][0] + itc[0][1][1]*Y[0][1];
    Y[1][0] = CB_SUB1[1][c[1]][0] + pred[0];
    Y[1][1] = CB_SUB1[1][c[1]][1] + pred[1];

    /* stage #3 - #4 */
    for(stage = 2; stage < N_STAGE_VQ-4; stage ++)
    {
        pred[0]		= itc[stage-1][0][0] * Y[stage-1][0] + itc[stage-1][0][1]*Y[stage-1][1];
        pred[1]		= itc[stage-1][1][0] * Y[stage-1][0] + itc[stage-1][1][1]*Y[stage-1][1];

        Y[stage][0] = CB_SUB2[stage-2][c[stage]][0] + pred[0];
        Y[stage][1] = CB_SUB2[stage-2][c[stage]][1] + pred[1];
    }

    /* stage #5 - #8 */
    for(stage = N_STAGE_VQ-4; stage < N_STAGE_VQ; stage ++)
    {
        pred[0]		= itc[stage-1][0][0] * Y[stage-1][0] + itc[stage-1][0][1]*Y[stage-1][1];
        pred[1]		= itc[stage-1][1][0] * Y[stage-1][0] + itc[stage-1][1][1]*Y[stage-1][1];

        Y[stage][0] = CB_SUB3[stage-4][c[stage]][0] + pred[0];
        Y[stage][1] = CB_SUB3[stage-4][c[stage]][1] + pred[1];
    }

    /* Transform Vector to Scalar */
    for(i = 0; i < N_STAGE_VQ; i++)
    {
        for(j = 0; j < N_DIM; j++)
        {
            y[i*N_DIM+j] = Y[i][j];
        }
    }

    return;
}

/*---------------------------------------------------------------------------
 * buildCode()
 *
 *
 *--------------------------------------------------------------------------*/

static void buildCode(
    short *ind,
    short fins,
    short c[],
    short s[]
)
{
    short stage;
    short BrIndex[4];

    set_s(BrIndex, 0, (N_STAGE_VQ - 4));


    for(stage = N_STAGE_VQ - 4; stage >= 1; stage--)
    {
        if(s[stage] > 7)
        {
            BrIndex[stage-1] = 1;
        }
    }
    ind[0] = fins;

    /* stage #1 - #2 */
    for(stage = 0; stage < 2; stage++)
    {
        ind[stage+1]  = BrIndex[stage]      << 4;
        ind[stage+1] += c[stage] >> 3;
    }

    /* stage #3 - #4 */
    for(stage = 2; stage < N_STAGE_VQ - 4; stage++)
    {
        ind[stage+1]  = BrIndex[stage]      << 3;
        ind[stage+1] += c[stage] >> 3;
    }

    /* Stage #5 - #8 */
    for(stage = N_STAGE_VQ-4; stage < N_STAGE_VQ; stage++)
    {
        ind[stage+1]  = c[stage] >> 3;
    }

    return;
}

/*---------------------------------------------------------------------------
 * BcTcvq()
 *
 *
 *--------------------------------------------------------------------------*/

static void BcTcvq(
    short snFlag,
    const float *x,
    float *y,
    const float *weight,
    short *ind
)
{
    float X[N_STAGE_VQ][N_DIM], W[N_STAGE_VQ][N_DIM];

    /* Count Variable */
    short i,j;

    /* TCVQ Structure */
    short stage, state, prev_state;
    short preState[N_STAGE_VQ][NUM_STATE];
    short codeWord[N_STAGE_VQ][NUM_STATE];
    float acumDist[N_STAGE_VQ-4][NUM_STATE];
    short inis, fins, ptr_fins;
    short fBlock;
    short fState[NUM_STATE];
    short fCodeword[4][4];
    short iniBlock[NUM_STATE];
    short blockCodeword[NUM_STATE][4];

    /* Prediction variable */
    float quant[N_STAGE_VQ][NUM_STATE][N_DIM];

    /* Distortion variable */
    float minDist;
    float fDist;
    float blockDist[NUM_STATE];

    /* Decoding variable */
    short bestCodeword[N_STAGE_VQ];
    short bestState[N_STAGE_VQ];

    /* Code Share variable */
    const float (*TCVQ_CB_SUB1)[128][2], (*TCVQ_CB_SUB2)[64][2], (*TCVQ_CB_SUB3)[32][2];
    const float (*IntraCoeff)[2][2];

    /* Memoryless Module */
    if(snFlag)
    {
        TCVQ_CB_SUB1 = SN_TCVQ_CB_SUB1;
        TCVQ_CB_SUB2 = SN_TCVQ_CB_SUB2;
        TCVQ_CB_SUB3 = SN_TCVQ_CB_SUB3;
        IntraCoeff	 = SN_IntraCoeff;
    }
    else /* Memory Module */
    {
        TCVQ_CB_SUB1 = AR_TCVQ_CB_SUB1;
        TCVQ_CB_SUB2 = AR_TCVQ_CB_SUB2;
        TCVQ_CB_SUB3 = AR_TCVQ_CB_SUB3;
        IntraCoeff	 = AR_IntraCoeff;
    }

    /* Transform Scalar to Vector */
    for(i = 0; i < N_STAGE_VQ; i++)
    {
        for(j = 0; j < N_DIM; j++)
        {
            X[i][j] = x[(N_DIM*i) + j];
            W[i][j] = weight[(N_DIM*i) + j];
        }
    }

    /* Initialzie */
    for(i=0; i<N_STAGE_VQ-4; i++)
    {
        for(j=0; j<NUM_STATE; j++)
        {
            acumDist[i][j]=0.f;
        }
    }

    /* BcTcvq Search */
    /* stage #1 */
    BcTcvq_1st(X, TCVQ_CB_SUB1, preState, codeWord, acumDist, quant, W);

    /* stage #2 */
    BcTcvq_2nd(X, TCVQ_CB_SUB1, preState, codeWord, acumDist, quant, W, IntraCoeff);

    /* stage #3 - #4 */
    for(stage = 2; stage < N_STAGE_VQ - 4; stage++)
    {
        BcTcvq_SubBlock(X, TCVQ_CB_SUB2, preState, codeWord, acumDist, quant, stage, W, IntraCoeff);
    }

    /* Search initial state at each block */
    for(state = 0; state < NUM_STATE; state++)
    {
        prev_state = state;

        for(stage = N_STAGE_VQ - 5; stage >= 0; stage--)
        {
            prev_state = preState[stage][prev_state];
        }
        iniBlock[state] = prev_state;
    }

    /* stage #5 - #8 */
    for(state = 0; state < NUM_STATE; state++)
    {
        inis	  = iniBlock[state];
        ptr_fins  = inis >> 2;

        minDist = FLT_MAX;

        for(i = 0; i < 4; i++)
        {
            fins  = ptr_fins*4 + i;

            prev_state = state;
            fDist = BcTcvq_FixSearch(X, TCVQ_CB_SUB3, fCodeword, quant, FixBranch, N_STAGE_VQ-4, inis, i, &prev_state, W, IntraCoeff);

            for(stage = N_STAGE_VQ-3; stage < N_STAGE_VQ; stage++)
            {
                fDist += BcTcvq_FixSearch(X, TCVQ_CB_SUB3, fCodeword, quant, FixBranch, stage, inis, i, &prev_state, W, IntraCoeff);
            }

            if(fDist < minDist)
            {
                minDist			 = fDist;
                blockDist[state] = minDist;
                fState[state]	 = fins;

                for(j = 0; j < 4; j++)
                {
                    blockCodeword[state][j] = fCodeword[i][j];
                }
            }
        }
    }

    /* Select optimal path */
    fBlock = optimalPath(acumDist, blockDist, blockCodeword, bestCodeword, codeWord, bestState, preState);

    /* Select Quantized Value */
    quantEnc(y, bestCodeword, TCVQ_CB_SUB1, TCVQ_CB_SUB2, TCVQ_CB_SUB3, IntraCoeff);

    /* Buid Code for Decoder */
    buildCode(ind, fState[fBlock], bestCodeword, bestState);

    return;
}


/*---------------------------------------------------------------------------
 * SVQ_2d()
 *
 *
 *--------------------------------------------------------------------------*/

static short SVQ_2d(
    float *x,
    float *y,
    const float *W,
    const float CB[][8],
    short Size
)
{
    short i, j;
    short index = 0;
    float distortion;
    float temp;

    temp = FLT_MAX;

    for(i = 0; i < Size; i++)
    {
        distortion = 0.0;
        for(j = 0; j < 8; j++)
        {
            distortion += (x[j] - CB[i][j])*(x[j] - CB[i][j])*W[j];
        }

        if(distortion < temp)
        {
            temp  = distortion;
            index = i;
        }
    }

    for(i = 0; i < M/2; i++)
    {
        y[i] = CB[index][i];
    }

    return index;
}


/*---------------------------------------------------------------------------
 * qlsf_ARSN_tcvq_Enc_16k()
 *
 * Predictive BC-TCQ encoder for LSF quantization
 *--------------------------------------------------------------------------*/

static float qlsf_ARSN_tcvq_Enc_16k (
    const float *x,           /* i  : Vector to be encoded    */
    float *y,           /* o  : Quantized LSF vector    */
    short *indice,      /* o  : Indices                 */
    const float *w,           /* i  : LSF Weights             */
    const short nBits,        /* i  : number of bits          */
    short safety_net    /* i  : safety_net flag         */
)
{
    short i;
    float temp_f;
    float yy[M];
    float error_svq[M], error_svq_q[M];
    float x_q[M];
    if(safety_net == 1)
    {
        indice[0] = 1;
        BcTcvq(1, x, x_q, w, &indice[1]);

        if(nBits>30)
        {
            /* SVQ  */
            for(i = 0; i < M; i++)
            {
                error_svq[i] = (x[i] - x_q[i]) * scale_inv_ARSN[i];
            }

            /* 5bits 1st Split VQ for Residual*/
            indice[10] = SVQ_2d(error_svq, error_svq_q, w, AR_SVQ_CB1, 32);
            /* 4bits 2nd Split VQ for Residual*/
            indice[11] = SVQ_2d(&error_svq[8], &error_svq_q[8], &w[8], AR_SVQ_CB2, 16 );

            for(i = 0; i < M; i++)
            {
                x_q[i] = x_q[i] + (error_svq_q[i] * scale_ARSN[i]);
            }
        }
    }
    else
    {
        indice[0] = 0;
        BcTcvq(0, x, x_q, w, &indice[1]);

        if(nBits>30)
        {
            /* SVQ */
            for(i = 0; i < M; i++)
            {
                error_svq[i] = x[i]   - x_q[i];
            }

            /* 5bits 1st Split VQ for Residual*/
            indice[10] = SVQ_2d(error_svq, error_svq_q, w, AR_SVQ_CB1, 32 );
            /* 4bits 2nd Split VQ for Residual*/
            indice[11] = SVQ_2d(&error_svq[8], &error_svq_q[8], &w[8], AR_SVQ_CB2, 16 );

            for(i = 0; i < M; i++)
            {
                x_q[i] = x_q[i] + error_svq_q[i];
            }
        }
    }

    v_sub( x_q, x, yy, M );
    v_mult( yy, yy, yy, M );
    temp_f = dotp( yy, w, M );

    /* Recover the quantized LSF */
    mvr2r(x_q, y, M);

    return temp_f;
}

/*-------------------------------------------------------------------*
 * lsf_mid_enc()
 *
 * Mid-frame LSF quantization
 --------------------------------------------------------------------*/

static void lsf_mid_enc(
    Encoder_State *st,          /* i/o: encoder state structure             */
    const float int_fs,         /* i  : internal (ACELP) sampling frequency */
    const float qlsp0[],        /* i  : quantized LSPs from frame beginning */
    const float qlsp1[],        /* i  : quantized LSPs from frame end       */
    float lsp[],                /* i/o: mid-frame LSP                       */
    const short coder_type,     /* i  : coding type                         */
    const short bwidth,         /* i  : input signal bandwidth              */
    const long  core_brate,     /* i  : core bitrate                        */
    float Bin_Ener[],           /* i  : per bin log energy spectrum         */
    short ppp_mode,
    short nelp_mode
)
{
    float ftemp, lsf[M], qlsf[M], qlsf1[M], qlsf0[M], wghts[M], err, err_min;
    short j, k, idx, nb_bits = 0, size = 0;
    const float *ratio = 0;

    /* convert LSPs to LSFs */
    lsp2lsf( lsp, lsf, M, int_fs );
    lsp2lsf( qlsp0, qlsf0, M, int_fs );
    lsp2lsf( qlsp1, qlsf1, M, int_fs );

    Unified_weighting( Bin_Ener, lsf, wghts, bwidth==NB, coder_type==UNVOICED, (int)(int_fs), M );

    /* codebook selection, number of bits, size of the codebook */
    if ( ppp_mode == 0 && nelp_mode == 0 )
    {
        nb_bits = mid_LSF_bits_tbl[LSF_BIT_ALLOC_IDX(core_brate, coder_type)];

        /* codebook selection */
        if ( coder_type == VOICED )
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
        else if ( coder_type == UNVOICED )
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

        size = (short) pow2[nb_bits];
    }
    else if ( ppp_mode == 1 )
    {
        ratio = tbl_mid_voi_wb_1b;
        nb_bits = 1;
        size = 2;
    }
    else if ( nelp_mode == 1 )
    {
        ratio = tbl_mid_unv_wb_4b;
        nb_bits = 4;
        size = 16;
    }

    /* loop over codevectors */
    err_min = 1e30f;
    idx = 0;
    for ( k = 0; k < size; k++ )
    {
        err = 0;

        for (j=0; j<M; j++)
        {
            qlsf[j] = (1.0f - ratio[k*M+j]) * qlsf0[j] + ratio[k*M+j] * qlsf1[j];

            if ( j > 0 && j < M && qlsf[j] < qlsf[j-1] + LSF_GAP_MID )
            {
                qlsf[j] = qlsf[j-1] + LSF_GAP_MID;
            }

            ftemp = lsf[j] - qlsf[j];
            err +=  wghts[j] * ftemp * ftemp;
        }

        if ( err < err_min )
        {
            err_min = err;
            idx = k;
        }
    }

    /* calculate the quantized LSF vector */
    for ( j = 0; j < M; j++ )
    {
        qlsf[j] = (1.0f - ratio[idx*M+j]) * qlsf0[j] + ratio[idx*M+j] * qlsf1[j];

        if ( j > 0 && j < M && qlsf[j] < qlsf[j-1] + LSF_GAP_MID )
        {
            qlsf[j] = qlsf[j-1] + LSF_GAP_MID;
        }
    }

    reorder_lsf( qlsf, LSF_GAP_MID, M, int_fs );

    /* convert LSFs back to LSPs */
    lsf2lsp( qlsf, lsp, M, int_fs );

    push_indice( st, IND_MID_FRAME_LSF_INDEX, idx, nb_bits );

    return;
}
