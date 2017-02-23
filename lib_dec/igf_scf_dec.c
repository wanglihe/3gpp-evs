/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <math.h>
#include "prot.h"
#include "options.h"
#include "stat_dec.h"


/*---------------------------------------------------------------------*
 * IGFSCFDecoderOpen()
 *
 * initialization of an instance of this module, pass a ptr to a hPublicData
 *---------------------------------------------------------------------*/

void IGFSCFDecoderOpen(
    IGFSCFDEC_INSTANCE_HANDLE                 hPublicData,      /**< inout: handle to public data */
    int                                       scfCountLongBlock,
    int                                       bitRate,
    int                                       mode,
    int                                       rf_mode             /**< in: flag to signal the RF mode */
)
{

    hPublicData->scfCountLongBlock  = scfCountLongBlock;
    hPublicData->t = 0; /* protect against the invalid request of starting decoding with a dependent block */

    IGFCommonFuncsIGFGetCFTables( bitRate, mode, rf_mode, &hPublicData->cf_se00, &hPublicData->cf_se01, &hPublicData->cf_off_se01,
                                  &hPublicData->cf_se02, &hPublicData->cf_off_se02, &hPublicData->cf_se10, &hPublicData->cf_off_se10,
                                  &hPublicData->cf_se11, &hPublicData->cf_off_se11 );

    return;
}

/*---------------------------------------------------------------------*
 * quant_ctx()
 *
 *
 *---------------------------------------------------------------------*/

static int quant_ctx(int ctx
                    )
{
    /*
      ctx ... -5 -4 -3 -2 -1 0 1 2 3 4 5 ...
    Q(ctx)... -3 -3 -3 -2 -1 0 1 2 3 3 3 ...
    */
    if (abs(ctx) <= 3)
    {
        return ctx;
    }
    else if (ctx > 3)
    {
        return 3;
    }
    else
    {
        /* ctx < -3 */
        return -3;
    }

}

/*---------------------------------------------------------------------*
 * arith_decode_bits()
 *
 *
 *---------------------------------------------------------------------*/

static int arith_decode_bits(
    IGFSCFDEC_INSTANCE_HANDLE  hPrivateData,   /**< instance handle */
    Decoder_State             *st,             /**< in: pointer to decoder state */
    int nBits                                  /**< number of bits to decode */
)
{
    int i;
    int x;
    int bit;

    x = 0;
    for (i = 0; i < nBits; ++i)
    {
        /* decode one bit using the new raw AC function */
        ari_decode_14bits_bit_ext(st, &bit, &hPrivateData->acState);
        x = (x << 1) | bit;
    }

    return x;
}


/*---------------------------------------------------------------------*
 * arith_decode_residual()
 *
 *
 *---------------------------------------------------------------------*/

static int arith_decode_residual(
    IGFSCFDEC_INSTANCE_HANDLE  hPrivateData,         /**< instance handle */
    Decoder_State             *st,                   /**< in: pointer to decoder state */
    const unsigned short* cumulativeFrequencyTable,  /**< cumulative frequency table to be used */
    int tableOffset                                  /**< offset used to align the table */
)
{
    int val;
    int x;
    int extra;
    int extra_tmp;

    /* decode one of the IGF_SYMBOLS_IN_TABLE == 27 alphabet symbols using the new raw AC function */
    ari_decode_14bits_s27_ext(st, &val, &hPrivateData->acState, cumulativeFrequencyTable);

    /* meaning of the values of val: */
    /* esc_{0} IGF_MIN_ENC_SEPARATE ... IGF_MAX_ENC_SEPARATE esc_{IGF_SYMBOLS_IN_TABLE - 1} */
    if ((val != 0) && (val != IGF_SYMBOLS_IN_TABLE - 1))
    {
        x = (val - 1) + IGF_MIN_ENC_SEPARATE;


        x -= tableOffset;
        return x;
    }

    /* decode one of the tails of the distribution */
    /* decode extra with 4 bits */
    extra = arith_decode_bits(hPrivateData, st, 4);
    if (extra == 15)
    {
        /* escape code 15 to indicate extra >= 15 */
        /* decode addtional extra with 6 bits */
        extra_tmp = arith_decode_bits(hPrivateData, st, 6);
        if (extra_tmp == 63)
        {
            /* escape code 63 to indicate extra_tmp >= 63 */
            /* decode safety extra with 7 bits */
            extra_tmp = 63 + arith_decode_bits(hPrivateData, st, 7);
        }
        extra = 15 + extra_tmp;
    }

    if (val == 0)
    {
        /* escape code 0 to indicate x <= IGF_MIN_ENC_SEPARATE - 1 */
        x = (IGF_MIN_ENC_SEPARATE - 1) - extra;
    }
    else
    {
        /* escape code (IGF_SYMBOLS_IN_TABLE - 1) to indicate x >= IGF_MAX_ENC_SEPARATE + 1 */
        x = (IGF_MAX_ENC_SEPARATE + 1) + extra;
    }

    x -= tableOffset;

    return x;
}


/*---------------------------------------------------------------------*
 * arith_decode_flush()
 *
 *
 *---------------------------------------------------------------------*/

static void arith_decode_flush(
    Decoder_State              *st                        /**< in: pointer to decoder state */
)
{
    get_next_indice_tmp(st, -14); /* return back the least significant 14 bits to the bitstream */

    return;
}


/*---------------------------------------------------------------------*
 * decode_sfe_vector()
 *
 *
 *---------------------------------------------------------------------*/

static void decode_sfe_vector(
    IGFSCFDEC_INSTANCE_HANDLE   hPrivateData,             /**< instance handle */
    Decoder_State              *st,                       /**< in: pointer to decoder state */
    int                         t,                        /**< frame counter reset to 0 at each independent frame */
    int                        *prev_x,                   /**< previous vector */
    int                        *x,                        /**< current vector to decode */
    int                         length                    /**< number of elements to decode */
)
{
    /*
       f
       ^
       |  d a x
       |    c b
       |      e  --> t
    */
    int f;
    int pred;
    int ctx;
    int ctx_f;
    int ctx_t;


    for (f = 0; f < length; f++)
    {
        if (t == 0)
        {
            if (f == 0)
            {
                /* decode one of the IGF_SYMBOLS_IN_TABLE == 27 alphabet symbols using the new raw AC function */
                ari_decode_14bits_s27_ext(st, &pred, &hPrivateData->acState, hPrivateData->cf_se00);
                x[f] = pred << 2;
                x[f] += arith_decode_bits(hPrivateData, st, 2); /* LSBs as 2 bit raw */
            }
            else if (f == 1)
            {
                pred = x[f - 1]; /* pred = b */
                x[f] = pred + arith_decode_residual( hPrivateData, st, hPrivateData->cf_se01, hPrivateData->cf_off_se01 );
            }
            else
            {
                /* f >= 2 */
                pred = x[f - 1]; /* pred = b */
                ctx = quant_ctx(x[f - 1] - x[f - 2]); /* Q(b - e) */
                x[f] = pred + arith_decode_residual( hPrivateData, st, &hPrivateData->cf_se02[(IGF_SYMBOLS_IN_TABLE + 1) * (IGF_CTX_OFFSET + ctx)],
                                                     hPrivateData->cf_off_se02[IGF_CTX_OFFSET + ctx] );
            }
        }
        else
        {
            /* t == 1 */
            if (f == 0)
            {
                pred = prev_x[f]; /* pred = a */
                x[f] = pred + arith_decode_residual( hPrivateData, st, hPrivateData->cf_se10, hPrivateData->cf_off_se10 );
            }
            else
            {
                /* (t == 1) && (f >= 1) */
                pred = prev_x[f] + x[f - 1] - prev_x[f - 1]; /* pred = a + b - c */
                ctx_f = quant_ctx(prev_x[f] - prev_x[f - 1]); /* Q(a - c) */
                ctx_t = quant_ctx(x[f - 1] - prev_x[f - 1]); /* Q(b - c) */
                x[f] = pred + arith_decode_residual(hPrivateData, st,
                                                    &hPrivateData->cf_se11[(IGF_SYMBOLS_IN_TABLE + 1) * IGF_CTX_COUNT * (IGF_CTX_OFFSET + ctx_t) + (IGF_SYMBOLS_IN_TABLE + 1) * (IGF_CTX_OFFSET + ctx_f)],
                                                    hPrivateData->cf_off_se11[IGF_CTX_COUNT * (IGF_CTX_OFFSET + ctx_t) + (IGF_CTX_OFFSET + ctx_f)]);
            }
        }
        if (x[f] < 0)
        {
            x[f] = 0;
            st->BER_detect = 1;
        }
        if (x[f] > 91)
        {
            x[f] = 91;
            st->BER_detect = 1;
        }
    }

    return;
}


/*---------------------------------------------------------------------*
 * IGFSCFDecoderReset()
 *
 * resets the internal decoder memory (context memory)
 *---------------------------------------------------------------------*/

void IGFSCFDecoderReset(
    IGFSCFDEC_INSTANCE_HANDLE           hPublicData       /**< inout: handle to public data or NULL in case there was no instance created */
)
{

    /* reset of coder */
    hPublicData->t       = 0;
    /* we do not need to fill hPublicData->prev with zeros, because when t = 0 no previous information is used */

    return;
}


/*---------------------------------------------------------------------*
 * IGFSCFDecoderDecode()
 *
 * main IGF decoder function
 *---------------------------------------------------------------------*/

void IGFSCFDecoderDecode(
    IGFSCFDEC_INSTANCE_HANDLE           hPublicData,       /**< inout: handle to public data or NULL in case there was no instance created */
    Decoder_State                      *st,                /**< inout: pointer to decoder state */
    int                                *sfe,               /**< out: ptr to an array which will contain the decoded quantized coefficients */
    int                                 indepFlag          /**< in: if  1 on input the encoder will be forced to reset,
                                                                              if  0 on input the encodder will be forced to encode without a reset */
)
{


    /* insert data */
    hPublicData->bitsRead = st->next_bit_pos;
    ari_start_decoding_14bits(st, &hPublicData->acState); /* start AC decoding */

    /* check if coder needs a reset and do it if neccessary */
    if (indepFlag)
    {
        /* reset of coder */
        IGFSCFDecoderReset(hPublicData);
    }

    decode_sfe_vector( hPublicData, st, hPublicData->t, hPublicData->prev, sfe, hPublicData->scfCountLongBlock );

    arith_decode_flush(st); /* finish decoding */


    /* advance history */
    mvi2i(sfe, hPublicData->prev, hPublicData->scfCountLongBlock);
    hPublicData->t++;

    hPublicData->bitsRead = st->next_bit_pos - hPublicData->bitsRead;

    return;
}
