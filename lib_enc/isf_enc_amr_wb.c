/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_enc.h"
#include "rom_com.h"
#include "prot.h"


/*-----------------------------------------------------------------*
 * Local constants
 *-----------------------------------------------------------------*/

#define DICO1_NS_19b     16       /* codebook dimensions for SID ISF quantizers */
#define DICO2_NS_19b     16
#define DICO3_NS_19b     16
#define DICO4_NS_19b     8
#define DICO5_NS_19b     16

#define DICO1_NS_28b     64
#define DICO2_NS_28b     64
#define DICO3_NS_28b     64
#define DICO4_NS_28b     32
#define DICO5_NS_28b     32

#define N_SURV_MAX       4        /* maximum number of survivors */

/*---------------------------------------------------------------------*
 * Local functions
 *---------------------------------------------------------------------*/

static void qisf_ns_28b( Encoder_State *st, float *isf );
static void qisf_2s_46b( Encoder_State *st, float *isf, short nb_surv, float *mem_AR, float *mem_MA );
static void qisf_2s_36b( Encoder_State *st, float *isf, short nb_surv, float *mem_AR, float *mem_MA );
static void VQ_stage1(const float *x, const float *dico, const short dim, const short dico_size, short *index, const short surv);
static short sub_VQ(float *x, const float *dico, const short dim, const short dico_size, float *distance);

/*-------------------------------------------------------------------*
 * isf_enc_amr_wb()
 *
 * Quantization of ISF parameters in AMR-WB IO mode
 *-------------------------------------------------------------------*/

void isf_enc_amr_wb(
    Encoder_State *st,          /* i/o: state structure                             */
    float *isf_new,     /* i/o  : quantized ISF vector                      */
    float *isp_new,     /* i/o: ISP vector to quantize/quantized            */
    float *Aq,          /* o  : quantized A(z) for 4 subframes              */
    float *stab_fac     /* o  : ISF stability factor                        */
)
{

    /*---------------------------------*
     * ISF quantization of SID frames
     *---------------------------------*/

    if ( st->core_brate == SID_1k75 )
    {
        qisf_ns_28b( st, isf_new );

        reorder_isf( isf_new, ISF_GAP, M, INT_FS_12k8 );

        isf2isp( isf_new, isp_new, M, INT_FS_12k8 );

        /* return if SID frame (conversion to A(z) done in the calling function) */
        return;
    }

    /* check resonance for pitch clipping algorithm */
    gp_clip_test_lsf( st->core_brate, isf_new, st->clip_var, 1 );

    /*---------------------------------------*
     * ISF quantization of all other frames
     *---------------------------------------*/

    if ( st->core_brate == ACELP_6k60 )
    {
        qisf_2s_36b( st, isf_new, 4, st->mem_AR, st->mem_MA);
    }
    else if( st->core_brate >= ACELP_8k85 )
    {
        qisf_2s_46b( st, isf_new, 4, st->mem_AR, st->mem_MA);
    }

    reorder_isf( isf_new, ISF_GAP, M, INT_FS_12k8 );

    /* convert quantized ISFs back to ISPs */
    isf2isp( isf_new, isp_new, M, INT_FS_12k8 );

    /*------------------------------------------------------------------*
     * ISP interpolation
     * A(z) calculation
     *------------------------------------------------------------------*/

    if( st->rate_switching_reset )
    {
        mvr2r( isf_new, st->lsf_old, M );
        mvr2r( isp_new, st->lsp_old, M );
    }

    int_lsp( L_FRAME, st->lsp_old, isp_new, Aq, M, interpol_isp_amr_wb, 1 );

    /*------------------------------------------------------------------*
     * Calculate ISF stability (distance between old ISF and current ISF)
     *------------------------------------------------------------------*/

    *stab_fac = lsf_stab( isf_new, st->lsf_old, 1, st->L_frame );


    return;
}

/*-------------------------------------------------------------------*
* qisf_ns_28b()
*
* ISF quantizer for SID frames (only in AMR-WB IO mode)
*-------------------------------------------------------------------*/

static void qisf_ns_28b(
    Encoder_State *st,          /* i/o: encoder state structure             */
    float *isf          /* i/o: unquantized/quantized ISF vector    */
)
{
    short i, indice[5];
    float tmp;

    for (i=0; i<M; i++)
    {
        isf[i] -= mean_isf_noise_amr_wb[i];
    }

    indice[0] = sub_VQ(&isf[0], dico1_ns_28b, 2, DICO1_NS_28b, &tmp);
    indice[1] = sub_VQ(&isf[2], dico2_ns_28b, 3, DICO2_NS_28b, &tmp);
    indice[2] = sub_VQ(&isf[5], dico3_ns_28b, 3, DICO3_NS_28b, &tmp);
    indice[3] = sub_VQ(&isf[8], dico4_ns_28b, 4, DICO4_NS_28b, &tmp);
    indice[4] = sub_VQ(&isf[12], dico5_ns_28b+4, 4, DICO5_NS_28b-1, &tmp) + 1;   /* First vector has a problem -> do not allow   */

    /* write indices to array */
    push_indice( st, IND_ISF_0_0, indice[0], 6 );
    push_indice( st, IND_ISF_0_1, indice[1], 6 );
    push_indice( st, IND_ISF_0_2, indice[2], 6 );
    push_indice( st, IND_ISF_0_3, indice[3], 5 );
    push_indice( st, IND_ISF_0_4, indice[4], 5 );

    /* decoding the ISFs */
    disf_ns_28b( indice, isf );

    return;
}


/*---------------------------------------------------------------------*
 * qisf_2s_36b()
 *
 * ISF quantizer for AMR-WB 6k60 frames
 *
 * The ISF vector is quantized using two-stage MA-prediction VQ with split-by-2
 * in 1st stage and split-by-3 in the second stage.
 *---------------------------------------------------------------------*/

static void qisf_2s_36b(
    Encoder_State *st,       /* i/o: encoder state structure              */
    float *isf,      /* i/o: unquantized/quantized ISF vector     */
    short nb_surv,   /* i  : number of survivors (1, 2, 3 or 4)   */
    float *mem_AR,   /* o  : quantizer memory for AR model        */
    float *mem_MA    /* i/o: quantizer memory for MA model        */
)
{
    short i, k, indice[5], tmp_ind[2];
    short surv1[N_SURV_MAX];     /* indices of survivors from 1st stage */
    float temp, min_err, distance, isf2[M];

    /*------------------------------------------------------------------------*
     * Subtract mean
     *------------------------------------------------------------------------*/

    for (i=0; i<M; i++)
    {
        isf[i] -= mean_isf_amr_wb[i] + MU_MA * mem_MA[i];
    }

    /*------------------------------------------------------------------------*
     * Quantize ISFs 0 - 8
     *------------------------------------------------------------------------*/

    VQ_stage1(&isf[0], dico1_isf, 9, SIZE_BK1, surv1, nb_surv);

    distance = 1.0e30f;
    if(nb_surv > N_SURV_MAX)
    {
        nb_surv = N_SURV_MAX;
    }

    for (k=0; k<nb_surv; k++)
    {
        for (i = 0; i < 9; i++)
        {
            isf2[i] = isf[i] - dico1_isf[i+surv1[k]*9];
        }

        tmp_ind[0] = sub_VQ(&isf2[0], dico21_isf_36b, 5, SIZE_BK21_36b, &min_err);
        temp = min_err;

        tmp_ind[1] = sub_VQ(&isf2[5], dico22_isf_36b, 4, SIZE_BK22_36b, &min_err);
        temp += min_err;

        if (temp < distance)
        {
            distance = temp;
            indice[0] = surv1[k];
            for (i=0; i<2; i++)
            {
                indice[i+2] = tmp_ind[i];
            }
        }
    }

    /*------------------------------------------------------------------------*
     * Quantize ISFs 9 - 15
     *------------------------------------------------------------------------*/

    VQ_stage1(&isf[9], dico2_isf, 7, SIZE_BK2, surv1, nb_surv);

    distance = 1.0e30f;
    for (k=0; k<nb_surv; k++)
    {
        for (i = 0; i < 7; i++)
        {
            isf2[9+i] = isf[9+i] - dico2_isf[i+surv1[k]*7];
        }

        tmp_ind[0] = sub_VQ(&isf2[9], dico23_isf_36b, 3, SIZE_BK23_36b, &min_err);
        temp = min_err;
        if (temp < distance)
        {
            distance = temp;
            indice[1] = surv1[k];
            indice[4] = tmp_ind[0];
        }
    }

    /*------------------------------------------------------------------------*
     * write indices to array
     *------------------------------------------------------------------------*/

    push_indice( st, IND_ISF_0_0, indice[0], 8 );
    push_indice( st, IND_ISF_0_1, indice[1], 8 );
    push_indice( st, IND_ISF_1_0, indice[2], 7 );
    push_indice( st, IND_ISF_1_1, indice[3], 7 );
    push_indice( st, IND_ISF_1_2, indice[4], 6 );

    /*------------------------------------------------------------------------*
     * decoding the ISFs
     *------------------------------------------------------------------------*/

    disf_2s_36b( indice, isf, mem_AR, mem_MA );

    return;
}


/*-------------------------------------------------------------------*
 * qisf_2s_46b()
 *
 * ISF quantizer for all other AMR-WB frames
 *
 * The ISF vector is quantized using two-stage VQ with split-by-2
 * in 1st stage and split-by-5 in the second stage.
 *-------------------------------------------------------------------*/

static void qisf_2s_46b(
    Encoder_State *st,       /* i/o: encoder state structure              */
    float *isf,      /* i/o: unquantized/quantized ISF vector     */
    short nb_surv,   /* i  : number of survivors (1, 2, 3 or 4)   */
    float *mem_AR,   /* o  : quantizer memory for AR model        */
    float *mem_MA    /* i/o: quantizer memory for MA model        */
)
{
    short i, k, indice[7], tmp_ind[5];
    short surv1[N_SURV_MAX];     /* indices of survivors from 1st stage */
    float temp, min_err, distance, isf2[M];


    /*------------------------------------------------------------------------*
     * Subtract mean
     *------------------------------------------------------------------------*/

    for (i=0; i<M; i++)
    {
        isf[i] -= mean_isf_amr_wb[i] + MU_MA * mem_MA[i];
    }

    /*------------------------------------------------------------------------*
     * Quantize ISFs 0 - 8
     *------------------------------------------------------------------------*/

    VQ_stage1(&isf[0], dico1_isf, 9, SIZE_BK1, surv1, nb_surv);

    distance = 1.0e30f;
    if(nb_surv > N_SURV_MAX)
    {
        nb_surv = N_SURV_MAX;
    }

    for (k=0; k<nb_surv; k++)
    {
        for (i = 0; i < 9; i++)
        {
            isf2[i] = isf[i] - dico1_isf[i+surv1[k]*9];
        }

        tmp_ind[0] = sub_VQ(&isf2[0], dico21_isf_46b, 3, SIZE_BK21, &min_err);
        temp = min_err;
        tmp_ind[1] = sub_VQ(&isf2[3], dico22_isf_46b, 3, SIZE_BK22, &min_err);
        temp += min_err;
        tmp_ind[2] = sub_VQ(&isf2[6], dico23_isf_46b, 3, SIZE_BK23, &min_err);
        temp += min_err;
        if (temp < distance)
        {
            distance = temp;
            indice[0] = surv1[k];
            for (i=0; i<3; i++)
            {
                indice[i+2]=tmp_ind[i];
            }
        }
    }

    /*------------------------------------------------------------------------*
     * Quantize ISFs 9 - 15
     *------------------------------------------------------------------------*/

    VQ_stage1(&isf[9], dico2_isf, 7, SIZE_BK2, surv1, nb_surv);

    distance = 1.0e30f;
    for (k=0; k<nb_surv; k++)
    {
        for (i = 0; i < 7; i++)
        {
            isf2[9+i] = isf[9+i] - dico2_isf[i+surv1[k]*7];
        }
        tmp_ind[0] = sub_VQ(&isf2[9], dico24_isf_46b, 3, SIZE_BK24, &min_err);
        temp = min_err;

        tmp_ind[1] = sub_VQ(&isf2[12], dico25_isf_46b, 4, SIZE_BK25, &min_err);
        temp += min_err;

        if (temp < distance)
        {

            distance = temp;
            indice[1] = surv1[k];
            for (i=0; i<2; i++)
            {
                indice[i+5]=tmp_ind[i];
            }
        }
    }

    /*------------------------------------------------------------------------*
     * write indices to array
     *------------------------------------------------------------------------*/

    push_indice( st, IND_ISF_0_0, indice[0], 8 );
    push_indice( st, IND_ISF_0_1, indice[1], 8 );
    push_indice( st, IND_ISF_1_0, indice[2], 6 );
    push_indice( st, IND_ISF_1_1, indice[3], 7 );
    push_indice( st, IND_ISF_1_2, indice[4], 7 );
    push_indice( st, IND_ISF_1_3, indice[5], 5 );
    push_indice( st, IND_ISF_1_4, indice[6], 5 );

    /*------------------------------------------------------------------------*
     * decoding the ISFs
     *------------------------------------------------------------------------*/

    disf_2s_46b( indice, isf, mem_AR, mem_MA );

    return;
}

/*-------------------------------------------------------------------*
 * VQ_stage1()
 *
 * 1st stage of ISF quantization
 *-------------------------------------------------------------------*/

static void VQ_stage1(
    const float *x,         /* i  : ISF vector                         */
    const float *dico,      /* i  : ISF codebook                       */
    const short dim,        /* i  : codebook dimension                 */
    const short dico_size,  /* i  : codebook size                      */
    short *index,     /* o  : indices of best vector candidates  */
    const short surv        /* i  : nb of surviving best candidates    */
)
{
    float dist_min[N_SURV_MAX];
    float dist, temp;
    const float *p_dico;
    short   i, j, k, l;


    for (i=0; i<surv; i++)
    {
        dist_min[i] = 1.0e30f;
        index[i] = i;
    }

    p_dico = dico;

    for (i = 0; i < dico_size; i++)
    {
        dist = 0.0;
        for (j = 0; j < dim; j++)
        {
            temp = x[j] - *p_dico++;
            dist += temp * temp;
        }

        for (k=0; k<surv; k++)
        {
            if (dist < dist_min[k])
            {
                for (l=surv-1; l>k; l--)
                {
                    dist_min[l] = dist_min[l-1];
                    index[l] = index[l-1];
                }
                dist_min[k] = dist;
                index[k] = i;
                break;
            }
        }
    }
    return;
}

/*-------------------------------------------------------------------*
 * sub_VQ()
 *
 * Quantization of a subvector in Split-VQ of ISFs
 *-------------------------------------------------------------------*/

static short sub_VQ(        /* o  : selected codebook vector index      */
    float *x,         /* i/o: ISF vector                          */
    const float *dico,      /* i  : ISF codebook                        */
    const short dim,        /* i  : codebook dimension                  */
    const short dico_size,  /* i  : codebook size                       */
    float *distance   /* o  : quantization error (min. distance)  */
)
{
    float dist_min, dist, temp;
    const float *p_dico;
    short   i, j, index;


    dist_min = 1.0e30f;
    p_dico = dico;

    index = 0;
    for (i = 0; i < dico_size; i++)
    {
        dist = 0.0f;
        for (j = 0; j < dim; j++)
        {
            temp = x[j] - *p_dico++;
            dist += temp * temp;
        }
        if (dist < dist_min)
        {
            dist_min = dist;
            index = i;
        }
    }

    *distance = dist_min;

    /* Reading the selected vector */
    p_dico = &dico[index * dim];
    for (j = 0; j < dim; j++)
    {
        x[j] = *p_dico++;
    }
    return index;
}
