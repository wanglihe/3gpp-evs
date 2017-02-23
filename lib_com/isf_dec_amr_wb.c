/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"

/*---------------------------------------------------------------------*
* isf_dec_amr_wb()
*
* Decoding of ISF parameters in AMR-WB IO mode
*---------------------------------------------------------------------*/

void isf_dec_amr_wb(
    Decoder_State *st,        /* i/o: State structure                             */
    float *Aq,          /* o  : quantized A(z) for 4 subframes              */
    float *isf_new,     /* o  : de-quantized ISF vector                     */
    float *isp_new      /* o  : de-quantized ISP vector                     */
)
{
    short i;
    short indice[7];

    set_s( indice, -1, 7 );

    /*---------------------------------*
     * ISF de-quantization of SID frames
     *---------------------------------*/

    if ( st->core_brate == SID_1k75 )
    {
        indice[0] = (short)get_next_indice( st, 6 );
        indice[1] = (short)get_next_indice( st, 6 );
        indice[2] = (short)get_next_indice( st, 6 );
        indice[3] = (short)get_next_indice( st, 5 );
        indice[4] = (short)get_next_indice( st, 5 );

        disf_ns_28b( indice, isf_new );

        reorder_isf( isf_new, ISF_GAP, M, INT_FS_12k8 );

        isf2isp( isf_new, isp_new, M, INT_FS_12k8 );

        /* return if SID frame (conversion to A(z) done in the calling function) */
        return;
    }

    /*-----------------------------------------------------------------*
     * ISF de-quantization of all other frames
     *-----------------------------------------------------------------*/

    if( st->core_brate == ACELP_6k60 )
    {
        indice[0] = (short)get_next_indice( st, 8 );
        indice[1] = (short)get_next_indice( st, 8 );
        indice[2] = (short)get_next_indice( st, 7 );
        indice[3] = (short)get_next_indice( st, 7 );
        indice[4] = (short)get_next_indice( st, 6 );

        disf_2s_36b( indice, isf_new, st->mem_AR, st->mem_MA );
    }
    else
    {
        indice[0] = (short)get_next_indice( st, 8 );
        indice[1] = (short)get_next_indice( st, 8 );
        indice[2] = (short)get_next_indice( st, 6 );
        indice[3] = (short)get_next_indice( st, 7 );
        indice[4] = (short)get_next_indice( st, 7 );
        indice[5] = (short)get_next_indice( st, 5 );
        indice[6] = (short)get_next_indice( st, 5 );

        disf_2s_46b( indice, isf_new, st->mem_AR, st->mem_MA);
    }

    reorder_isf( isf_new, ISF_GAP, M, INT_FS_12k8 );

    /* convert quantized ISFs to ISPs */
    isf2isp( isf_new, isp_new, M, INT_FS_12k8 );

    /*-------------------------------------------------------------------------------------*
     * FEC - update adaptive mean ISF vector
     *-------------------------------------------------------------------------------------*/

    for ( i=0; i<M; i++ )
    {
        st->lsf_adaptive_mean[i] = (st->lsfoldbfi1[i] + st->lsfoldbfi0[i] + isf_new[i]) / 3;
    }

    /*-------------------------------------------------------------------------------------*
     * ISP interpolation
     * A(z) calculation
     *-------------------------------------------------------------------------------------*/

    if( st->rate_switching_reset )
    {
        /*extrapolation instead of interpolation*/
        mvr2r( isp_new, st->lsp_old, M );
        mvr2r( isf_new, st->lsf_old, M );
    }

    /* ISP interpolation and A(z) calculation */
    int_lsp( L_FRAME, st->lsp_old, isp_new, Aq, M, interpol_isp_amr_wb, 1 );

    /*------------------------------------------------------------------*
     * Check ISF stability : distance between old ISF and current ISF
     *------------------------------------------------------------------*/

    st->stab_fac = lsf_stab( isf_new, st->lsf_old, 1, st->L_frame );

    return;
}

/*-------------------------------------------------------------------*
 * disf_ns_28b()
 *
 * ISF de-quantizer for SID_1k75 frames (only for AMR-WB IO mode)
 *-------------------------------------------------------------------*/

void disf_ns_28b(
    short *indice,      /* i  : quantized indices (use indice[0] = -1 in the decoder) */
    float *isf_q        /* o  : ISF in the frequency domain (0..6400) */
)
{
    short i;

    for (i = 0; i < 2; i++)
    {
        isf_q[i] =  dico1_ns_28b[indice[0]*2+i];
    }

    for (i = 0; i < 3; i++)
    {
        isf_q[i+2] =  dico2_ns_28b[indice[1]*3+i];
        isf_q[i+5] =  dico3_ns_28b[indice[2]*3+i];
    }

    for (i = 0; i < 4; i++)
    {
        isf_q[i+8] =  dico4_ns_28b[indice[3]*4+i];
        isf_q[i+12] =  dico5_ns_28b[indice[4]*4+i];
    }

    for (i=0; i<M; i++)
    {
        isf_q[i] += mean_isf_noise_amr_wb[i];
    }

    return;
}

/*-------------------------------------------------------------------*
 * disf_2s_46b()
 *
 * ISF de-quantizer for 46b. codebooks (only for AMR-WB IO mode)
 *-------------------------------------------------------------------*/

void disf_2s_46b(
    short *indice,    /* i  : quantized indices (use indice[0] = -1 in the decoder) */
    float *isf_q,     /* o  : quantized ISFs in the cosine domain */
    float *mem_AR,    /* o  : quantizer memory for AR model       */
    float *mem_MA     /* i/o: quantizer memory for MA model       */
)
{
    short i;

    for (i = 0; i < 9; i++)
    {
        isf_q[i] =  dico1_isf[indice[0]*9+i];
    }

    for (i = 0; i < 7; i++)
    {
        isf_q[i+9] =  dico2_isf[indice[1]*7+i];
    }

    for (i = 0; i < 3; i++)
    {
        isf_q[i] +=  dico21_isf_46b[indice[2]*3+i];
        isf_q[i+3] +=  dico22_isf_46b[indice[3]*3+i];
        isf_q[i+6] +=  dico23_isf_46b[indice[4]*3+i];
        isf_q[i+9] +=  dico24_isf_46b[indice[5]*3+i];
    }

    for (i = 0; i < 4; i++)
    {
        isf_q[i+12] +=  dico25_isf_46b[indice[6]*4+i];
    }

    for (i = 0; i < M; i++)
    {
        mem_AR[i] = (float)(isf_q[i] + MU_MA * mem_MA[i]);  /* Update with quantized ISF vector for AR model */
        mem_MA[i] = isf_q[i];                               /* Update with quantized prediction error for MA model */
        isf_q[i] = mem_AR[i] + mean_isf_amr_wb[i];          /* Quantized ISFs */
    }

    return;
}


/*-------------------------------------------------------------------*
 * disf_2s_36b()
 *
 * ISF de-quantizer for 36b. codebooks (only for AMR-WB IO mode)
 *-------------------------------------------------------------------*/

void disf_2s_36b(
    short *indice,    /* i  : quantized indices (use indice[0] = -1 in the decoder) */
    float *isf_q,     /* o  : quantized ISFs in the cosine domain */
    float *mem_AR,    /* i/o: quantizer memory for AR model       */
    float *mem_MA     /* i/o: quantizer memory for MA model       */
)
{
    short i;
    const float *pt_dico1;

    pt_dico1 = dico1_isf;    /* Pointer of the 1st stage, 1st plit */

    for (i = 0; i < 9; i++)
    {
        isf_q[i] =  pt_dico1[indice[0]*9+i];
    }

    for (i = 0; i < 7; i++)
    {
        isf_q[i+9] =  dico2_isf[indice[1]*7+i];
    }

    for (i = 0; i < 5; i++)
    {
        isf_q[i] +=  dico21_isf_36b[indice[2]*5+i];
    }

    for (i = 0; i < 4; i++)
    {
        isf_q[i+5] +=  dico22_isf_36b[indice[3]*4+i];
    }

    for (i = 0; i < 7; i++)
    {
        isf_q[i+9] +=  dico23_isf_36b[indice[4]*7+i];
    }

    for (i = 0; i < M; i++)
    {
        mem_AR[i] = (float)(isf_q[i] + MU_MA * mem_MA[i]); /* Update with quantized ISF vector for AR model */
        mem_MA[i] = isf_q[i];                              /* Update with quantized prediction error for MA model */
        isf_q[i] = mem_AR[i] + mean_isf_amr_wb[i];         /* Quantized ISFs */
    }

    return;
}
