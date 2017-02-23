/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <assert.h>
#include "prot.h"
#include <stdio.h>
#include <stdlib.h>


/*------------------------------------------------------------------*
* dlpc_avq()
*
* Variable bit-rate multiple LPC un-quantizer
*------------------------------------------------------------------*/

int dlpc_avq(
    int *index,       /* (i)   Quantization indices                       */
    float *LSF_Q,     /* (o)   Quantized LSF vectors                      */
    int numlpc,       /* (i) Number of sets of lpc */
    float sr_core
)
{
    int i, nbi, last;
    int *p_index, q_type;

    /* Last LPC index */

    if ( numlpc==1 )
    {
        last = 0;
    }
    else
    {
        last = M;
    }

    p_index = index;

    /* Decode last LPC */
    for (i=0; i<M; i++)
    {
        LSF_Q[last+i] = 0.0f;
    }

    vlpc_1st_dec(p_index[0], &LSF_Q[last], sr_core );
    p_index++;
    vlpc_2st_dec(&LSF_Q[last], &p_index[0], 0, sr_core );
    nbi = 2 + p_index[0] + p_index[1];
    p_index += nbi;

    /* Decode intermediate LPC (512 framing) */

    if ( numlpc==2 )
    {
        q_type = p_index[0];
        p_index++;

        if (q_type == 0)
        {

            for (i=0; i<M; i++)
            {
                LSF_Q[i] = 0.0f;
            }
            vlpc_1st_dec(p_index[0], &LSF_Q[0], sr_core );
            p_index++;
            vlpc_2st_dec(&LSF_Q[0], &p_index[0], 0, sr_core );
        }
        else if (q_type == 1)
        {
            for (i=0; i<M; i++)
            {
                LSF_Q[i] = LSF_Q[M+i];
            }

            vlpc_2st_dec(&LSF_Q[0], &p_index[0], 3, sr_core );

        }

        nbi = 2 + p_index[0] + p_index[1];
        p_index += nbi;
    }

    return p_index-index;
}

/*------------------------------------------------------------------*
* unary_decode()
*
*
*------------------------------------------------------------------*/

static int unary_decode(
    Decoder_State *st,
    int *ind
)
{
    int start_bit_pos;

    start_bit_pos = st->next_bit_pos;

    /* Index bits */
    *ind = 0;

    while (get_next_indice_1(st) && !st->BER_detect)
    {
        *ind += 1;
    }

    if (*ind != 0)
    {
        *ind += 1;
    }

    return st->next_bit_pos - start_bit_pos;

}


/*------------------------------------------------------------------*
* pack4bits()
*
*
*------------------------------------------------------------------*/

static int pack4bits(
    int nbits,
    Decoder_State *st,
    int *prm
)
{
    int i;

    i=0;

    while (nbits > 4)
    {
        prm[i] = get_next_indice(st, 4);
        nbits -= 4;
        i++;
    }
    prm[i] = get_next_indice(st, nbits);
    i++;

    return(i);
}


/*------------------------------------------------------------------*
* decode_lpc_avq()
*
*
*------------------------------------------------------------------*/

int decode_lpc_avq(
    Decoder_State *st,
    int numlpc,
    int *param_lpc
)
{
    int k,j;
    int nb, qn1, qn2, avqBits, q_type;
    int start_bit_pos;


    j = 0;
    start_bit_pos = st->next_bit_pos;


    for (k=0; k<numlpc; k++)
    {
        /* Decode quantizer type */

        if (k==0)
        {
            q_type = 0;
            nb = 0;
        }
        else
        {
            nb = 1;
            q_type = get_next_indice(st, nb);
            param_lpc[j++] = q_type;
        }

        /* Decode quantization indices */

        if (q_type==0)
        {
            /* Absolute quantizer with 1st stage stochastic codebook */
            param_lpc[j++] = get_next_indice(st, 8);
        }

        /* 2 bits to specify Q2,Q3,Q4,ext */
        qn1 = 2 + get_next_indice(st, 2);
        qn2 = 2 + get_next_indice(st, 2);

        /* Unary code */
        /* Q5 = 0, Q6=10, Q0=110, Q7=1110, ... */

        if (qn1 > 4)
        {
            nb = unary_decode(st, &qn1);

            if (nb == 1) qn1 += 5;
            else if (nb == 2) qn1 += 4;
            else if (nb == 3) qn1 = 0;
            else qn1 += 3;
        }

        if (qn2 > 4)
        {
            nb = unary_decode(st, &qn2);

            if (nb == 1) qn2 += 5;
            else if (nb == 2) qn2 += 4;
            else if (nb == 3) qn2 = 0;
            else qn2 += 3;
        }

        /* check for potential bit errors */
        if( qn1 > NB_SPHERE || qn2 > NB_SPHERE )
        {
            qn1 = 0;
            qn2 = 0;
            st->BER_detect = 1;
        }

        param_lpc[j] = qn1;
        j++;
        param_lpc[j] = qn2;
        j++;

        /* Decode Split-by-2 algebraic VQ */
        avqBits = 4*qn1;
        pack4bits(avqBits, st, &param_lpc[j]);
        j += qn1;

        avqBits = 4*qn2;
        pack4bits(avqBits, st, &param_lpc[j]);
        j += qn2;
    }

    return st->next_bit_pos - start_bit_pos;
}
