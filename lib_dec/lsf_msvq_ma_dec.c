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
#include "basop_proto_func.h"

/*---------------------------------------------------------------------*
 * lsf_msvq_ma_decprm()
 *
 *
 *---------------------------------------------------------------------*/

int lsf_msvq_ma_decprm(
    Decoder_State *st,
    int *param_lpc,
    int core,
    int acelp_mode,
    int acelp_midLpc,
    int narrowBand,
    int sr_core
)
{
    int i, nbits_lpc;
    int  bits_midlpc=5;
    short bits0[MAX_VQ_STAGES], bits1[MAX_VQ_STAGES], stages0, stages1, stages,
          levels0[MAX_VQ_STAGES], levels1[MAX_VQ_STAGES], * bits;
    short predmode, mode_lvq, mode_lvq_p,  safety_net;


    if( (sr_core==INT_FS_16k)&&(acelp_mode==UNVOICED) )
    {
        predmode = find_pred_mode(GENERIC, 1-narrowBand/*st->bwidth*/, sr_core,
                                  &mode_lvq, &mode_lvq_p, st->total_brate);
    }
    else
    {
        if (core == TCX_20_CORE)
        {
            predmode = find_pred_mode(AUDIO, 1-narrowBand/*st->bwidth*/, sr_core,
                                      &mode_lvq, &mode_lvq_p, st->total_brate );
        }
        else
        {
            predmode = find_pred_mode(acelp_mode, 1-narrowBand/*st->bwidth*/, sr_core,
                                      &mode_lvq, &mode_lvq_p, st->total_brate );
        }
    }

    lsf_allocate( 31-(predmode>>1), mode_lvq, mode_lvq_p, &stages0, &stages1, levels0, levels1, bits0, bits1 );

    nbits_lpc = 0;

    if (predmode == 2)
    {
        /* there is choice between SN and AR prediction */
        safety_net = get_next_indice(st, 1);

        if (safety_net==1)
        {
            stages = stages0;
            bits = bits0;
        }
        else
        {
            stages = stages1;
            bits = bits1;
        }
        *param_lpc = safety_net;
        param_lpc++;
        nbits_lpc++;

    }
    else
    {
        stages = stages1;
        bits = bits1;
    }

    for (i=0; i<stages-1; i++)
    {
        *param_lpc = get_next_indice(st, bits[i]);
        param_lpc++;
        nbits_lpc += bits[i];

    }
    *param_lpc = get_next_indice(st, LEN_INDICE);
    param_lpc++;
    nbits_lpc += LEN_INDICE;


    *param_lpc = get_next_indice(st, bits[i]-LEN_INDICE);
    param_lpc++;
    nbits_lpc += bits[i]-LEN_INDICE;


    if ( acelp_mode!=VOICED && core==0 && acelp_midLpc)
    {

        *param_lpc = get_next_indice(st, bits_midlpc);
        nbits_lpc += bits_midlpc;
    }

    return nbits_lpc;
}


/*---------------------------------------------------------------------*
 * lsf_bctcvq_decprm()
 *
 *
 *---------------------------------------------------------------------*/

int lsf_bctcvq_decprm(
    Decoder_State * st,
    int *param_lpc
)
{
    int i, nbits_lpc;
    int num_par;
    const short *bits1;

    num_par = 10;
    bits1 = BC_TCVQ_BIT_ALLOC_40B;

    nbits_lpc = 0;

    for (i=0; i<num_par; i++)
    {
        *param_lpc = get_next_indice(st, bits1[i]);
        param_lpc++;
        nbits_lpc += bits1[i];
    }

    return nbits_lpc;
}


/*---------------------------------------------------------------------*
 * D_lsf_tcxlpc()
 *
 *
 *---------------------------------------------------------------------*/

int D_lsf_tcxlpc(            /* (O) number of indices */
    const int indices[],       /* (I) VQ indices        */
    float lsf_q[],             /* (O) quantized lsf     */
    Word16 lsp_q_ind[],        /* (O) quantized lsp (w/o MA prediction) */
    int narrowband,            /* (I) narrowband flag   */
    int cdk,                   /* (I) codebook selector */
    float mem_MA[]             /* (I) MA memory         */
)
{
    int i;
    int NumIndices;
    float pred[M16k];
    const float *means;
    Word16 lsf_q_ind[M16k];
    float lsf_rem_q[M];
    Word16 lsf_rem_q_ind[M];

    NumIndices = 1;

    msvq_dec( lsf_codebook[narrowband][cdk], lsf_dims, lsf_offs, TCXLPC_NUMSTAGES, M, M, indices + NumIndices, lsf_q, lsf_q_ind );

    NumIndices += TCXLPC_NUMSTAGES;

    if (indices[0])
    {
        /* Only add contribution if flag is enabled */

        msvq_dec( lsf_ind_codebook[narrowband][cdk], lsf_ind_dims, lsf_ind_offs, TCXLPC_IND_NUMSTAGES, M, M, indices + NumIndices, lsf_rem_q, lsf_rem_q_ind );
        NumIndices += TCXLPC_IND_NUMSTAGES;

        /* Add to MA-removed vector */
        for (i=0; i<M; ++i)
        {
            lsf_q_ind[i] = add(lsf_q_ind[i], lsf_rem_q_ind[i]);
        }
    }

    /* Inter-frame prediction */
    means = lsf_means[narrowband];

    for (i=0; i<M; ++i)
    {
        pred[i] = means[i] + MU_MA * mem_MA[i];
    }

    /* Add prediction */
    for (i=0; i<M; ++i)
    {
        lsf_q[i] += pred[i];
        lsf_q_ind[i] = add(lsf_q_ind[i], LSFM(means[i]));
    }
    reorder_lsf(lsf_q, TCXLPC_LSF_GAP, M, INT_FS_12k8);

    basop_reorder_lsf(lsf_q_ind, LSF_GAP_VAL(TCXLPC_LSF_GAP), M, INT_FS_FX);

    if (lsp_q_ind)
    {
        basop_lsf2lsp(lsf_q_ind, lsp_q_ind);
    }

    return NumIndices;
}


/*---------------------------------------------------------------------*
 * dec_lsf_tcxlpc()
 *
 *
 *---------------------------------------------------------------------*/

int dec_lsf_tcxlpc(          /* (O) number of bits read */
    Decoder_State *st,         /* (I/O) Decoder state   */
    int **indices,             /* (O) Ptr to VQ indices */
    int narrowband,            /* (I) narrowband flag   */
    int cdk                    /* (I) codebook selector */
)
{
    int i, start_bit_pos;
    float lsf_q_ignored[M];
    Word16 lsf_q_ind[M];
    int *flag;

    flag = *indices; /* Save pointer */
    *flag = 0; /* Set flag to disabled */
    ++*indices;

    start_bit_pos = st->next_bit_pos;

    for (i=0; i<TCXLPC_NUMSTAGES; ++i)
    {
        **indices = get_next_indice(st, lsf_numbits[i]);
        ++*indices;
    }

    /* Decode independent lsf */
    msvq_dec( lsf_codebook[narrowband][cdk], lsf_dims, lsf_offs, TCXLPC_NUMSTAGES, M, M, flag+1, lsf_q_ignored, lsf_q_ind );

    /* Update flag */
    *flag = lsf_ind_is_active(lsf_q_ind, lsf_means[narrowband], narrowband, cdk);

    if (*flag)
    {
        for (i=0; i<TCXLPC_IND_NUMSTAGES; ++i)
        {
            **indices = get_next_indice(st, lsf_ind_numbits[i]);
            ++*indices;
        }
    }

    return st->next_bit_pos - start_bit_pos;
}

