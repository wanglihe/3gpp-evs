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
#include "stl.h"
#include "basop_proto_func.h"

/*---------------------------------------------------------------------*
 * midlsf_dec()
 *
 *
 *---------------------------------------------------------------------*/

void midlsf_dec(
    float qlsf0[],
    float qlsf1[],
    short idx,
    float qlsf[],
    int N,
    int coder_type,
    short *mid_lsf_int,
    short prev_bfi,
    short safety_net
)
{
    const float *ratio=NULL;
    int   j;
    short bad_spacing = 0;
    /* Select codebook */
    if ( coder_type == UNVOICED )
    {
        ratio = tbl_mid_unv_wb_5b;
    }
    else
    {
        ratio = tbl_mid_gen_wb_5b;
    }

    for (j=0; j<N; j++)
    {
        qlsf[j] = (1.0f - ratio[idx*N+j]) * qlsf0[j] + ratio[idx*N+j] * qlsf1[j];
    }

    if(mid_lsf_int != NULL) /*at the decoder*/
    {
        /* check for incorrect LSF ordering */
        if ( *mid_lsf_int == 1 )
        {
            for (j=1; j<N; j++)
            {
                if ( qlsf[j] < qlsf[j-1] )
                {
                    bad_spacing = 1;
                    break;
                }
            }
        }
        /* Redo mid-LSF interpolation with 0.4 in case of LSF instability */
        if ( prev_bfi || ( *mid_lsf_int == 1 && bad_spacing ) )
        {
            for (j=0; j<N; j++)
            {
                /* redo mid-LSF interpolation with 0.4 */
                qlsf[j] = 0.4f * qlsf0[j] + 0.6f * qlsf1[j];

                /* ensure correct ordering of LSF indices */
                if ( j > 0 && j < N && qlsf[j] < qlsf[j-1] + LSF_GAP_MID  )
                {
                    qlsf[j] = qlsf[j-1] + LSF_GAP_MID;
                }

            }
        }
        else
        {
            /* otherwise, use regular LSF spacing and ordering as in the encoder */
            for (j=0; j<N; j++)
            {
                if ( j > 0 && j < N && qlsf[j] < qlsf[j-1] + LSF_GAP_MID )
                {
                    qlsf[j] = qlsf[j-1] + LSF_GAP_MID;
                }

            }
        }
        if ( prev_bfi )
        {
            /* continue redoing mid-LSF interpolation with 0.4 in order not to propagate the error */
            *mid_lsf_int = 1;
        }
        if ( safety_net )
        {
            /* safety-net encountered -> stop redoing mid-LSF interpolation with 0.4 */
            *mid_lsf_int = 0;
        }
    }
    else
    {
        /* use regular LSF spacing */
        for (j=0; j<N; j++)
        {
            if ( j > 0 && j < N && qlsf[j] < qlsf[j-1] + LSF_GAP_MID )
            {
                qlsf[j] = qlsf[j-1] + LSF_GAP_MID;
            }
        }
    }

    return;
}


/*---------------------------------------------------------------------*
 * lsf_ind_is_active()
 *
 *
 *---------------------------------------------------------------------*/

int lsf_ind_is_active(
    const Word16 lsf_q_ind[],
    const float means[],
    int narrowband,
    int cdk
)
{
    Word16 lsf[2], min_distance;

    lsf[0] = add(lsf_q_ind[0], LSFM(means[0]));
    move16();
    lsf[1] = add(lsf_q_ind[1], LSFM(means[1]));
    move16();

    min_distance = lsf[0];
    min_distance = s_min(min_distance, sub(lsf[1], lsf[0]));

    assert(narrowband == 0 || narrowband == 1);
    assert(cdk == 0 || cdk == 1);

    return sub(min_distance, min_distance_thr[narrowband][cdk]) < 0;
}
