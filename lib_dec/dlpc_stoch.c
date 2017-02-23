/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "options.h"
#include "prot.h"
#include "rom_com.h"
#include "basop_proto_func.h"

/*------------------------------------------------------------------*
* lpc_unquantize()
*
*
*------------------------------------------------------------------*/

void lpc_unquantize(
    Decoder_State * st,
    float *lsfold,
    float *lspold,
    float *lsf,
    float *lsp,
    int lpcQuantization,
    int *param_lpc,
    int numlpc,
    int core,
    float *mem_MA,
    float *lspmid,
    float *lsfmid,
    short coder_type,
    int acelp_midLpc,
    int narrow_band,
    short *seed_acelp,
    int sr_core,
    short *mid_lsf_int,
    short prev_bfi,
    short *LSF_Q_prediction,  /* o  : LSF prediction mode                     */
    short *safety_net
)
{
    int nb_indices=0, k;

    mvr2r(lsfold, &lsf[0], M);
    mvr2r(lspold, &lsp[0], M);

    if( lpcQuantization == 0 )
    {
        nb_indices = dlpc_avq(param_lpc, &lsf[M], numlpc, sr_core );
        for ( k=0; k<numlpc; k++ )
        {
            lsf2lsp(&lsf[(k+1)*M], &lsp[(k+1)*M], M, sr_core);
        }
    }
    else if( lpcQuantization == 1 )
    {
        if( sr_core == INT_FS_16k && coder_type == UNVOICED )
        {
            lsf_end_dec( st, GENERIC, 1-narrow_band /* st->bwidth */ , ENDLSF_NBITS, &lsf[M], st->mem_AR, mem_MA, sr_core, st->core_brate,
                         &st->offset_scale1[0][0], &st->offset_scale2[0][0], &st->offset_scale1_p[0][0], &st->offset_scale2_p[0][0],
                         &st->no_scales[0][0], &st->no_scales_p[0][0], &st->safety_net, param_lpc, LSF_Q_prediction, &nb_indices);
        }
        else
        {
            if (st->core == TCX_20_CORE)
            {
                lsf_end_dec( st, AUDIO, 1-narrow_band /* st->bwidth */ , ENDLSF_NBITS, &lsf[M], st->mem_AR, mem_MA, sr_core, st->core_brate,
                             &st->offset_scale1[0][0], &st->offset_scale2[0][0], &st->offset_scale1_p[0][0], &st->offset_scale2_p[0][0],
                             &st->no_scales[0][0], &st->no_scales_p[0][0], &st->safety_net, param_lpc, LSF_Q_prediction, &nb_indices);
            }
            else
            {
                lsf_end_dec( st, coder_type, 1-narrow_band /* st->bwidth */ , 31, &lsf[M], st->mem_AR, mem_MA, sr_core, st->core_brate,
                             &st->offset_scale1[0][0], &st->offset_scale2[0][0], &st->offset_scale1_p[0][0], &st->offset_scale2_p[0][0],
                             &st->no_scales[0][0], &st->no_scales_p[0][0], &st->safety_net, param_lpc, LSF_Q_prediction, &nb_indices);
            }
        }

        lsf2lsp(&lsf[M], &lsp[M], M, sr_core);
    }
    else
    {
        assert(0);
    }

    *seed_acelp=0;
    for(k=nb_indices-1; k>=0; k--)
    {
        /* rightshift before *seed_acelp+param_lpc[i] to avoid overflows*/
        *seed_acelp=(short)((((*seed_acelp)>>1)+param_lpc[k]) * 31821L + 13849L);
    }

    /* Decoded mid-frame lsf */
    if( lpcQuantization && acelp_midLpc && core == ACELP_CORE && st->rate_switching_reset == 0 )
    {
        midlsf_dec( &lsf[0], &lsf[M], (short)param_lpc[nb_indices], lsfmid, M, coder_type, mid_lsf_int, prev_bfi, *safety_net );

        reorder_lsf( lsfmid, LSF_GAP_MID, M, sr_core );
        lsf2lsp( lsfmid, lspmid, M, sr_core );
    }


    return;
}
