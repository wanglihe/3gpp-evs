/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"
#include "rom_enc.h"
#include "basop_proto_func.h"


#define kMaxC    8

/*--------------------------------------------------------------------------*
 * msvq_enc()
 *
 * MSVQ encoder
 *--------------------------------------------------------------------------*/

void msvq_enc
(
    const float *const *cb,    /* i  : Codebook (indexed cb[*stages][levels][p])     */
    const int dims[],   /* i  : Dimension of each codebook stage (NULL: full dim.)   */
    const int offs[],   /* i  : Starting dimension of each codebook stage (NULL: 0)  */
    float u[],          /* i  : Vector to be encoded (prediction and mean removed)   */
    const int *levels,  /* i  : Number of levels in each stage                       */
    int maxC,           /* i  : Tree search size (number of candidates kept from     */
    /*      one stage to the next == M-best)                     */
    int stages,         /* i  : Number of stages                                     */
    float w[],          /* i  : Weights                                              */
    int N,              /* i  : Vector dimension                                     */
    int maxN,           /* i  : Codebook dimension                                   */
    int Idx[]           /* o  : Indices                                              */
)
{
    float *resid[2], *dist[2];
    float en, tmp, *pTmp,*p1;
    const float *cb_stage, *cbp, *p2;
    int *indices[2], j, m, s,c, c2, p_max, i;
    float resid_buf[2*LSFMBEST_MAX*M_MAX], dist_buf[2*LSFMBEST_MAX], Tmp[M_MAX];
    int idx_buf[2*LSFMBEST_MAX*MAX_VQ_STAGES_USED], parents[LSFMBEST_MAX];
    int n, maxn, start;

    /*----------------------------------------------------------------*
     * Allocate memory for previous (parent) and current nodes.
     *   Parent node is indexed [0], current node is indexed [1].
     *----------------------------------------------------------------*/

    indices[0] = idx_buf;
    indices[1] = idx_buf + maxC*stages;
    set_i( idx_buf, 0, 2*stages*maxC );

    resid[0] = resid_buf;
    pTmp = resid_buf + maxC*N;
    resid[1] = pTmp;

    dist[0] = dist_buf;
    dist[1] = dist_buf + maxC;
    set_i( parents, 0, maxC);

    /* Set up initial distance vector */
    for(tmp=0.0, j=0; j<N; j++)
    {
        tmp += u[j]*u[j]*w[j] ;
    }
    set_f( dist[1], tmp, maxC);

    /* Set up initial error (residual) vectors */
    for(c=0; c<maxC; c++)
    {
        for (i = 0; i < N; i++)
        {
            *pTmp++ = u[i];
        }
    }

    /* Loop over all stages */
    for (m=1, s=0; s<stages; s++)
    {
        /* codebook pointer is set to point to first stage */
        cbp = cb[s];
        /* Save pointer to beginning of current stage */
        cb_stage = cbp;

        /* Set up pointers to parent and current nodes */
        swap (indices[0], indices[1], int*);
        swap (resid[0], resid[1], float*);
        swap (dist[0], dist[1], float*);

        /* p_max points to maximum distortion node (worst of best) */
        p_max = 0;

        if (dims)
        {
            n = dims[s];
            maxn = dims[s];
        }
        else
        {
            n = N;
            maxn = maxN;
        }

        if (offs)
        {
            start = offs[s];
        }
        else
        {
            start = 0;
        }

        set_zero(Tmp, start);
        set_zero(Tmp + start + n, N - (start + n));
        /* Set distortions to a large value */
        for (j=0; j < maxC; j++)
        {
            dist[1][j] = FLT_MAX;
        }
        if (!s) /* means: m==1 */
        {
            /* This loop is identical to the one below, except, that the inner
               loop over c=0..m is hardcoded to c=0, since m=1.                 */
            /* dist[0][0] */
            for (j=0; j<levels[s]; j++)
            {
                en = 0.0f;
                /* w,Tmp */
                /* Compute weighted codebook element and its energy */
                for (c2=0; c2 < n; c2++)
                {
                    Tmp[start+c2] = w[start+c2] * cbp[c2];
                    en += cbp[c2] * Tmp[start+c2];
                }
                cbp += maxn;                                              /*  pointer is incremented */

                pTmp = &resid[0][0];
                /* Tmp */
                tmp = (*pTmp++) * Tmp[0];
                for ( c2=1; c2<N; c2++ )
                {
                    tmp += (*pTmp++) * Tmp[c2];
                }
                tmp = en - 2.0f * tmp;
                tmp += dist[0][0];
                if (tmp < dist[1][p_max])
                {
                    /* Replace worst */
                    dist[1][p_max] = tmp;
                    indices[1][p_max*stages] = j;
                    parents[p_max] = 0;

                    p_max = 0;
                    for (c2=1; c2 < maxC; c2++)
                    {
                        if (dist[1][c2] > dist[1][p_max])
                        {
                            p_max = c2;
                        }
                    }
                } /* if (tmp <= dist[1][p_max]) */
            } /* for (j=0; j<levels[s]; j++) */
        }
        else
        {
            /* dist[0][0] */
            for (j=0; j<levels[s]; j++)
            {
                en = 0.0f;
                /* w,Tmp */
                /* Compute weighted codebook element and its energy */
                for (c2=0; c2 < n; c2++)
                {
                    Tmp[start+c2] = w[start+c2] * cbp[c2];
                    en += cbp[c2] * Tmp[start+c2];
                }
                cbp += maxn;                                              /*  pointer is incremented */

                /* dist[0][0] */
                pTmp = &resid[0][0];
                /* Iterate over all parent nodes */
                for(c=0; c<m; c++)
                {
                    /* Tmp[0] */
                    tmp = (*pTmp++) * Tmp[0];
                    for ( c2=1; c2<N; c2++ )
                    {
                        tmp += (*pTmp++) * Tmp[c2];
                    }
                    tmp = en - 2.0f * tmp;
                    tmp += dist[0][c];
                    if (tmp < dist[1][p_max])
                    {
                        /* Replace worst */
                        dist[1][p_max] = tmp;
                        indices[1][p_max*stages+s] = j;
                        parents[p_max] = c;

                        p_max = 0;
                        for (c2=1; c2 < maxC; c2++)
                        {
                            if (dist[1][c2] > dist[1][p_max])
                            {
                                p_max = c2;
                            }
                        }
                    } /* if (tmp <= dist[1][p_max]) */
                } /* for(c=0; c<m; c++) */
            } /* for (j=0; j<levels[s]; j++) */
        }

        /*------------------------------------------------------------*
         * Compute error vectors for each node
         *------------------------------------------------------------*/
        /* parents[0] */
        pTmp = resid[1];
        for (c=0; c<maxC; c++)
        {
            /* Subtract codebook entry from residual vector of parent node and multiply with scale factor */
            p1 = resid[0]+parents[c]*N;
            p2 = cb_stage+(indices[1][c*stages+s])*maxn;
            mvr2r(p1, pTmp, start);
            for (j=0; j<n; j++)
            {
                pTmp[start+j] = (p1[start+j] - p2[j]);
            }
            mvr2r(p1+start+n, pTmp+start+n, N-(start+n));
            pTmp += N;                                                /*  pointer is incremented */

            /* Get indices that were used for parent node */
            mvi2i(indices[0]+parents[c]*stages, indices[1]+c*stages, s);
        } /* for (c=0; c<maxC; c++) */

        m = maxC;
    } /* for (m=1, s=0; s<stages; s++) */

    /* Find the optimum candidate */
    c2 = minimum (dist[1], maxC, 0);
    mvi2i (indices[1]+c2*stages, Idx, stages);

    return;
}

/*--------------------------------------------------------------------------*
 * lsf_msvq_ma_encprm()
 *
 *
 *--------------------------------------------------------------------------*/

int lsf_msvq_ma_encprm(
    Encoder_State * st,
    int *param_lpc,
    int core,
    int acelp_mode,
    int acelp_midLpc,
    short * bits_param_lpc,
    short no_indices
)
{
    int i, nbits_lpc;
    int bits_midlpc = 5;

    nbits_lpc = 0;
    for (i=0; i<no_indices; i++)
    {
        push_next_indice(st, *param_lpc, bits_param_lpc[i]);
        param_lpc++;
        nbits_lpc += bits_param_lpc[i];
    }
    if ( acelp_mode!=VOICED )
    {
        if ( core==0 && acelp_midLpc)
        {
            push_next_indice(st, *param_lpc, bits_midlpc);
            nbits_lpc += bits_midlpc;
        }
    }

    return nbits_lpc;
}

/*--------------------------------------------------------------------------*
 * lsf_msvq_ma_encprm()
 *
 *
 *--------------------------------------------------------------------------*/

void midlsf_enc(
    float qlsf0[],
    float qlsf1[],
    const float lsf[],
    short *idx,
    int N,
    float *Bin_Ener,
    int narrowBand,
    int sr_core,
    int coder_type
)
{
    float pred[M], wghts[M], err, err_min, tmp;
    int   NS=0, j, k;
    const float *ratio=NULL;

    /* Select codebook */
    if ( coder_type == UNVOICED )
    {
        ratio = tbl_mid_unv_wb_5b;
    }
    else
    {
        ratio = tbl_mid_gen_wb_5b;
    }
    NS = 32;

    /* Weights */
    Unified_weighting(Bin_Ener, lsf, wghts, narrowBand, coder_type==UNVOICED, sr_core,M);
    err_min = FLT_MAX;
    *idx = 0;
    for ( k = 0; k < NS; k++ )
    {
        err = 0;

        for (j=0; j<N; j++)
        {
            pred[j] = (1.0f - ratio[k*N+j]) * qlsf0[j] + ratio[k*N+j] * qlsf1[j];

            if ( j > 0 && j < N && pred[j] < pred[j-1] + LSF_GAP_MID )
            {
                pred[j] = pred[j-1] + LSF_GAP_MID;
            }

            tmp = lsf[j] - pred[j];
            err +=  wghts[j] * tmp * tmp;
        }

        if ( err < err_min )
        {
            err_min = err;
            *idx = k;
        }
    }

    return;
}


/*--------------------------------------------------------------------------*
 * Q_lsf_tcxlpc()
 *
 *
 *--------------------------------------------------------------------------*/

int Q_lsf_tcxlpc(           /* o  : number of indices                   */
    /* const */ float lsf[],  /* (I) original lsf                         */
    float lsf_q[],            /* (O) quantized lsf                        */
    Word16 lsp_q_ind[],       /* (O) quantized lsp (w/o MA prediction)    */
    int indices[],            /* (O) VQ indices                           */
    int narrowband,           /* (I) narrowband flag                      */
    int cdk,                  /* (I) codebook selector                    */
    float mem_MA[],           /* (I) MA memory                            */
    int coder_type,           /* (I) acelp extended mode                  */
    float *Bin_Ener           /* (I) Spectrum energy                      */
)
{
    float weights[M];
    float pred[M16k];
    int i;
    int NumIndices;
    const float *means;
    Word16 lsf_q_ind[M16k];
    float lsf_rem[M];
    float lsf_rem_q[M];
    Word16 lsf_rem_q_ind[M];

    Unified_weighting( &Bin_Ener[L_FFT/2], lsf, weights, narrowband, coder_type==UNVOICED, 12800, M );

    NumIndices = 0;

    /* Put disabled flag */
    indices[NumIndices++] = 0;

    /* Inter-frame prediction */
    means = lsf_means[narrowband];
    for (i=0; i<M; ++i)
    {
        pred[i] = means[i] + MU_MA * mem_MA[i];
    }

    /* Subtract prediction */
    for (i=0; i<M; ++i)
    {
        lsf[i] -= pred[i];
    }

    msvq_enc( lsf_codebook[narrowband][cdk], lsf_dims, lsf_offs, lsf, lsf_numlevels,
              kMaxC, TCXLPC_NUMSTAGES, weights, M, M, indices + NumIndices );
    msvq_dec( lsf_codebook[narrowband][cdk], lsf_dims, lsf_offs, TCXLPC_NUMSTAGES, M,
              M, indices + NumIndices, lsf_q, lsf_q_ind );

    NumIndices += TCXLPC_NUMSTAGES;

    /* Update flag */
    indices[0] = lsf_ind_is_active(lsf_q_ind, lsf_means[narrowband], narrowband, cdk);

    /* Get residual vector */
    for (i=0; i<M; ++i)
    {
        lsf_rem[i] = (pred[i] + lsf[i]) - (lsf_means[narrowband][i] + lsf_q_ind[i]/(float)(2.0f*1.28f));
    }
    /* Quantize using extra stage(s) */
    msvq_enc( lsf_ind_codebook[narrowband][cdk], lsf_ind_dims, lsf_ind_offs, lsf_rem, lsf_ind_numlevels,
              kMaxC, TCXLPC_IND_NUMSTAGES, weights, M, M, indices + NumIndices );

    /* Only add contribution if flag is enabled */
    if (indices[0])
    {
        /* Decode */
        msvq_dec( lsf_ind_codebook[narrowband][cdk], lsf_ind_dims, lsf_ind_offs, TCXLPC_IND_NUMSTAGES, M, M,
                  indices + NumIndices, lsf_rem_q, lsf_rem_q_ind );
        NumIndices += TCXLPC_IND_NUMSTAGES;

        /* Add to MA-removed vector */
        for (i=0; i<M; ++i)
        {
            lsf_q_ind[i] = add(lsf_q_ind[i], lsf_rem_q_ind[i]);
        }
    }

    /* Add inter-frame prediction */
    for (i=0; i<M; ++i)
    {
        lsf_q[i] += pred[i];
        lsf[i] += pred[i];
    }

    reorder_lsf(lsf_q, TCXLPC_LSF_GAP, M, INT_FS_12k8);
    for (i=0; i<M; ++i)
    {
        lsf_q_ind[i] = add(lsf_q_ind[i], LSFM(lsf_means[narrowband][i]));
        move16();
    }

    basop_reorder_lsf(lsf_q_ind, LSF_GAP_VAL(TCXLPC_LSF_GAP), M, INT_FS_FX);
    if( lsp_q_ind )
    {
        basop_lsf2lsp(lsf_q_ind, lsp_q_ind);
    }

    return NumIndices;
}


/*--------------------------------------------------------------------------*
 * enc_lsf_tcxlpc()
 *
 *
 *--------------------------------------------------------------------------*/

int enc_lsf_tcxlpc(          /* Returns: number of bits written */
    int **indices,             /* (I) Ptr to VQ indices */
    Encoder_State *st          /* (I/O) Encoder state   */
)
{
    int i, NumBits;
    int flag;

    /* Read flag */
    flag = (*indices)[0];
    ++*indices;

    NumBits = TCXLPC_NUMBITS;
    for (i=0; i<TCXLPC_NUMSTAGES; ++i)
    {
        push_next_indice(st, **indices, lsf_numbits[i]);
        ++*indices;
    }
    if (flag)
    {
        NumBits += TCXLPC_IND_NUMBITS;
        for (i=0; i<TCXLPC_IND_NUMSTAGES; ++i)
        {
            push_next_indice(st, **indices, lsf_ind_numbits[i]);
            ++*indices;
        }
    }
    return NumBits;
}


/*--------------------------------------------------------------------------*
 * lsf_bctcvq_encprm()
 *
 *
 *--------------------------------------------------------------------------*/

int lsf_bctcvq_encprm(
    Encoder_State *st,
    int *param_lpc,
    short * bits_param_lpc,
    short no_indices
)
{
    short i, nbits_lpc;

    nbits_lpc = 0;

    for (i=0; i<no_indices; i++)
    {
        push_next_indice(st, *param_lpc, bits_param_lpc[i]);
        param_lpc++;
        nbits_lpc += bits_param_lpc[i];
    }

    return nbits_lpc;
}
