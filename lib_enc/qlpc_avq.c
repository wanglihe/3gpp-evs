/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <stdio.h>
#include "prot.h"


/*-------------------------------------------------------------------*
* qlpc_avq()
*
*
*--------------------------------------------------------------------*/

void qlpc_avq(
    const float *lsf,        /* (i) Input LSF vectors              */
    const float *lsfmid,
    float *lsf_q,            /* (o) Quantized LFS vectors          */
    float *lsfmid_q,
    int *index,              /* (o) Quantization indices           */
    int *nb_indices,         /* (o) Number of quantization indices */
    int *nbbits,             /* (o) Number of quantization bits    */
    int core,
    float sr_core
)
{
    int i;
    float lsfmid_q0[M];
    int *tmp_index, indxt[256], nbits, nbt, nit;
    float dummy[M];

    /* Init */
    tmp_index = &index[0];
    *nb_indices = 0;

    /* Quantize end_frame LPC */
    for (i=0; i<M; i++)
    {
        lsf_q[i] = 0.0f;
    }

    tmp_index[0] = vlpc_1st_cod(lsf, lsf_q, sr_core, dummy );

    nbt = vlpc_2st_cod( lsf, lsf_q, &tmp_index[1], 0, sr_core );

    nit = 1 + 2 + index[1] + index[2];
    tmp_index += nit;
    *nb_indices += nit;
    nbbits[0] = 8 + nbt;

    *tmp_index = 0;

    /* Quantize mid-frame LPC */

    if( core == TCX_10_CORE )
    {
        tmp_index++;
        *nb_indices +=1;

        /* LPC2: Abs? */
        for (i=0; i<M; i++)
        {
            lsfmid_q[i] = 0.0f;
        }

        tmp_index[0] = vlpc_1st_cod(lsfmid, lsfmid_q, sr_core, dummy );

        nbits = vlpc_2st_cod(lsfmid, lsfmid_q, &tmp_index[1], 0, sr_core);

        nbt = 8 + nbits;
        nit = 1 + 2 + tmp_index[1] + tmp_index[2];

        /* LPC2: RelR? */
        for (i=0; i<M; i++)
        {
            lsfmid_q0[i] = lsf_q[i];
        }

        nbits = vlpc_2st_cod(lsfmid, lsfmid_q0, indxt, 3, sr_core);

        if (nbits < nbt)
        {
            nbt = nbits;
            nit = 2 + indxt[0] + indxt[1];
            tmp_index[-1] = 1;

            for (i=0; i<M; i++)
            {
                lsfmid_q[i] = lsfmid_q0[i];
            }

            for (i=0; i<nit; i++)
            {
                tmp_index[i] = indxt[i];
            }
        }

        tmp_index += nit;
        *nb_indices += nit;
        nbbits[1] = 1 + nbt;
    }

    return;
}


/*-------------------------------------------------------------------*
* unary_code()
*
*
*--------------------------------------------------------------------*/

static int unary_code(int ind, Encoder_State *st)
{
    int nb_bits;

    nb_bits = 1;

    /* Index bits */
    ind -= 1;

    while (ind >= 16)
    {
        push_next_indice(st, 0xffffU, 16);
        nb_bits += 16;
        ind -= 16;
    }
    if (ind > 0)
    {
        push_next_indice(st, (1U<<ind)-1, ind);
        nb_bits += ind;
    }

    /* Stop bit */
    push_next_indice(st, 0, 1);

    return(nb_bits);
}


/*-------------------------------------------------------------------*
* unpack4bits()
*
*
*--------------------------------------------------------------------*/

static int unpack4bits(int nbits, const int *prm, Encoder_State *st)
{
    int i;

    if (nbits == 0)
    {
        push_next_indice(st, 0, 0);
        i = 1;
    }
    else
    {
        i=0;

        while (nbits > 4)
        {
            push_next_indice(st, prm[i], 4);
            nbits -= 4;
            i++;
        }
        push_next_indice(st, prm[i], nbits);
        i++;
    }

    return(i);
}


/*-------------------------------------------------------------------*
* encode_lpc_avq()
*
*
*--------------------------------------------------------------------*/

int encode_lpc_avq(
    Encoder_State *st,
    int numlpc,
    int *param_lpc,
    int mode
)
{
    int k,j;
    int q_type, nb_ind;
    int i,qn1,qn2,nb,avqBits,st1=0;
    int nb_bits;

    j = 0;
    nb_bits = 0;

    for (k=0; k<numlpc; k++)
    {
        /* Retrieve quantizer type */
        if (k==0)
        {
            q_type = 0;
        }
        else
        {
            q_type = param_lpc[j++];
        }

        /* Determine number of AVQ indices */
        nb_ind = 0;

        if (q_type==0)
        {
            st1 = param_lpc[j++];
        }
        qn1 = param_lpc[j++];
        qn2 = param_lpc[j++];
        nb_ind = qn1 + qn2;

        if ( k==0 || (k==1 && mode!=1) )
        {
            /* Encode quantizer type */
            if (k==0)
            {
                nb = 0;
            }
            else
            {
                nb = 1;
                push_next_indice(st, q_type, nb);
            }
            nb_bits += nb;

            /* Encode quantizer data */

            if (q_type==0)
            {
                /* Absolute quantizer with 1st stage stochastic codebook */
                push_next_indice(st, st1, 8);
                nb_bits += 8;
            }

            /* 2 bits to specify Q2,Q3,Q4,ext */
            nb_bits += 4;
            i = qn1-2;

            if ((i<0) || (i>3))
            {
                i = 3;
            }
            push_next_indice(st, i, 2);

            i = qn2-2;

            if ((i<0) || (i>3))
            {
                i = 3;
            }
            push_next_indice(st, i, 2);

            /* Unary code for abs and rel LPC0/LPC2 */
            /* Q5 = 0, Q6=10, Q0=110, Q7=1110, ... */
            nb = qn1;

            if (nb > 6)
            {
                nb -= 3;
            }
            else if (nb > 4)
            {
                nb -= 4;
            }
            else if (nb == 0)
            {
                nb = 3;
            }
            else
            {
                nb = 0;
            }

            if (nb > 0)
            {
                unary_code(nb, st);
            }
            nb_bits += nb;

            nb = qn2;

            if (nb > 6)
            {
                nb -= 3;
            }
            else if (nb > 4)
            {
                nb -= 4;
            }
            else if (nb == 0)
            {
                nb = 3;
            }
            else
            {
                nb = 0;
            }

            if (nb > 0)
            {
                unary_code(nb, st);
            }
            nb_bits += nb;

            avqBits = 4*qn1;
            unpack4bits(avqBits, &param_lpc[j], st);
            j += qn1;
            nb_bits += avqBits;

            avqBits = 4*qn2;
            unpack4bits(avqBits, &param_lpc[j], st);
            j += qn2;
            nb_bits += avqBits;
        }
        else
        {
            j += nb_ind;
        }
    }

    return(nb_bits);
}
