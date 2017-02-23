/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <assert.h>
#include "prot.h"
#include "options.h"
#include "stat_enc.h"
#include "stat_com.h"
#include "cnst.h"

/*---------------------------------------------------------------------*
 * IGFSCFEncoderOpen()
 *
 * initialization of an instance of this module, pass a ptr to a hPublicData
 *---------------------------------------------------------------------*/

void IGFSCFEncoderOpen(
    IGFSCFENC_INSTANCE_HANDLE                  hPublicData,       /**< inout: handle to public data */
    int                                        scfCountLongBlock,
    int                                        bitRate,
    int                                        mode,
    int                                        rf_mode             /**< in: flag to signal the RF mode */
)
{

    hPublicData->ptrBitIndex       = 0;
    hPublicData->bitCount          = 0;
    hPublicData->Tsave             = 0;
    hPublicData->contex_saved      = 0;
    hPublicData->acState.low       = 0;
    hPublicData->acState.high      = 0;
    hPublicData->acState.vobf      = 0;
    set_i(hPublicData->prev, 0, 64);
    set_i(hPublicData->prevSave, 0, 64);

    hPublicData->scfCountLongBlock  = scfCountLongBlock;
    hPublicData->t = 0; /* protect against the invalid request of starting encoding with a dependent block */

    IGFCommonFuncsIGFGetCFTables( bitRate, mode, rf_mode, &hPublicData->cf_se00, &hPublicData->cf_se01, &hPublicData->cf_off_se01,
                                  &hPublicData->cf_se02, &hPublicData->cf_off_se02, &hPublicData->cf_se10,
                                  &hPublicData->cf_off_se10, &hPublicData->cf_se11, &hPublicData->cf_off_se11 );

    return;
}


/*---------------------------------------------------------------------*
 * quant_ctx()
 *
 *
 *---------------------------------------------------------------------*/

static int quant_ctx(
    int ctx
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
    else     /* ctx < -3 */
    {
        return -3;
    }
}


/*---------------------------------------------------------------------*
 * arith_encode_bits()
 *
 *
 *---------------------------------------------------------------------*/

static void arith_encode_bits(
    IGFSCFENC_INSTANCE_HANDLE  hPrivateData,    /**< instance handle */
    int                       *ptr,             /**< pointer to expanded bit buffer, one bit in each int */
    int                        x,               /**< value to encode */
    int                        nBits            /**< number of bits to encode */
)
{
    int i;
    int bit;

    for (i = nBits - 1; i >= 0; --i)
    {
        bit = (x >> i) & 1;
        hPrivateData->ptrBitIndex = ari_encode_14bits_sign(ptr, hPrivateData->ptrBitIndex, 32767, /* disable the bit count limitation */ &hPrivateData->acState, bit );
    }

    return;
}


/*---------------------------------------------------------------------*
 * arith_encode_residual()
 *
 *
 *---------------------------------------------------------------------*/

static void arith_encode_residual(
    IGFSCFENC_INSTANCE_HANDLE  hPrivateData,              /**< instance handle */
    int                       *ptr,                       /**< pointer to expanded bit buffer, one bit in each int */
    int                        x,                         /**< prediction residual to encode */
    const unsigned short      *cumulativeFrequencyTable,  /**< cumulative frequency table to be used */
    int                        tableOffset                /**< offset used to align the table */
)
{
    int extra;

    x += tableOffset;
    if ((x >= IGF_MIN_ENC_SEPARATE) && (x <= IGF_MAX_ENC_SEPARATE))
    {
        /* encode one of the IGF_SYMBOLS_IN_TABLE == 27 alphabet symbols using the new raw AC function */
        hPrivateData->ptrBitIndex = ari_encode_14bits_ext( ptr, hPrivateData->ptrBitIndex, &hPrivateData->acState, (x - IGF_MIN_ENC_SEPARATE) + 1, cumulativeFrequencyTable );

        return;
    }
    else if (x < IGF_MIN_ENC_SEPARATE)
    {
        /* send escape code 0 to indicate x <= IGF_MIN_ENC_SEPARATE - 1 */
        extra = (IGF_MIN_ENC_SEPARATE - 1) - x;
        hPrivateData->ptrBitIndex = ari_encode_14bits_ext( ptr, hPrivateData->ptrBitIndex, &hPrivateData->acState, 0, cumulativeFrequencyTable );
    }
    else
    {
        /* x > IGF_MAX_ENC_SEPARATE */
        /* send escape code (IGF_SYMBOLS_IN_TABLE - 1) to indicate x >= IGF_MAX_ENC_SEPARATE + 1 */
        extra = x - (IGF_MAX_ENC_SEPARATE + 1);
        hPrivateData->ptrBitIndex = ari_encode_14bits_ext( ptr, hPrivateData->ptrBitIndex, &hPrivateData->acState, IGF_SYMBOLS_IN_TABLE - 1, cumulativeFrequencyTable );
    }

    /* encode one of the tails of the distribution */
    if (extra < 15)
    {
        /* encode extra with 4 bits if extra < 15 */
        arith_encode_bits(hPrivateData, ptr, extra, 4);
    }
    else
    {
        /* extra >= 15 */
        /* send escape code 15 to indicate extra >= 15 */
        arith_encode_bits(hPrivateData, ptr, 15, 4);
        extra -= 15;

        if (extra < 63)
        {
            /* encode additional extra with 6 bits */
            arith_encode_bits(hPrivateData, ptr, extra, 6);
        }
        else     /* extra >= 63 */
        {
            arith_encode_bits(hPrivateData, ptr, 63, 6);
            extra -= 63;
            /* encode safety extra with 7 bits */
            arith_encode_bits(hPrivateData, ptr, extra, 7);
        }
    }

    return;
}


/*---------------------------------------------------------------------*
 * encode_sfe_vector()
 *
 *
 *---------------------------------------------------------------------*/

static void encode_sfe_vector(
    IGFSCFENC_INSTANCE_HANDLE           hPrivateData,     /**< instance handle */
    int                                *ptr,              /**< pointer to expanded bit buffer, one bit in each int */
    int                                 t,                /**< frame counter reset to 0 at each independent frame */
    int                                *prev_x,           /**< previous vector */
    int                                *x,                /**< current vector to encode */
    int                                 length,           /**< number of elements to encode */
    int                                 do_real_encoding  /**< whether the real encoding is needed, otherwise only the number of bits is used */
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


    (void) do_real_encoding; /* avoid compiler warning: unreferenced formal parameter */

    for (f = 0; f < length; f++)
    {
        if (t == 0)
        {
            if (f == 0)
            {
                /* encode one of the IGF_SYMBOLS_IN_TABLE == 27 alphabet symbols using the new raw AC function */
                hPrivateData->ptrBitIndex = ari_encode_14bits_ext( ptr, hPrivateData->ptrBitIndex, &hPrivateData->acState, x[f] >> 2, hPrivateData->cf_se00 );
                arith_encode_bits(hPrivateData, ptr, x[f] & 3, 2); /* LSBs as 2 bit raw */
            }
            else if (f == 1)
            {
                pred = x[f - 1]; /* pred = b */
                arith_encode_residual( hPrivateData, ptr, x[f] - pred, hPrivateData->cf_se01, hPrivateData->cf_off_se01 );
            }
            else
            {
                /* f >= 2 */
                pred = x[f - 1]; /* pred = b */
                ctx = quant_ctx(x[f - 1] - x[f - 2]); /* Q(b - e) */
                arith_encode_residual( hPrivateData, ptr, x[f] - pred, &hPrivateData->cf_se02[(IGF_SYMBOLS_IN_TABLE + 1) * (IGF_CTX_OFFSET + ctx)],
                                       hPrivateData->cf_off_se02[IGF_CTX_OFFSET + ctx] );
            }
        }
        else
        {
            /* t == 1 */
            if (f == 0)
            {
                pred = prev_x[f]; /* pred = a */
                arith_encode_residual( hPrivateData, ptr, x[f] - pred, hPrivateData->cf_se10, hPrivateData->cf_off_se10 );
            }
            else
            {
                /* (t == 1) && (f >= 1) */
                pred = prev_x[f] + x[f - 1] - prev_x[f - 1]; /* pred = a + b - c */
                ctx_f = quant_ctx(prev_x[f] - prev_x[f - 1]); /* Q(a - c) */
                ctx_t = quant_ctx(x[f - 1] - prev_x[f - 1]); /* Q(b - c) */
                arith_encode_residual( hPrivateData, ptr, x[f] - pred,
                                       &hPrivateData->cf_se11[(IGF_SYMBOLS_IN_TABLE + 1) * IGF_CTX_COUNT * (IGF_CTX_OFFSET + ctx_t) + (IGF_SYMBOLS_IN_TABLE + 1) * (IGF_CTX_OFFSET + ctx_f)],
                                       hPrivateData->cf_off_se11[IGF_CTX_COUNT * (IGF_CTX_OFFSET + ctx_t) + (IGF_CTX_OFFSET + ctx_f)]);
            }
        }
    }

    return;
}


/*---------------------------------------------------------------------*
 * IGFSCFEncoderReset()
 *
 * Reset of Arith enc context memory
 *---------------------------------------------------------------------*/

int IGFSCFEncoderReset(
    IGFSCFENC_INSTANCE_HANDLE           hPublicData
)
{


    hPublicData->t       = 0;
    /* we do not need to fill hPublicData->prev with zeros, because when t = 0 no previous information is used */

    return 0;
}


/*---------------------------------------------------------------------*
 * IGFSCFEncoderEncode()
 *
 * main IGF encoder function
 *---------------------------------------------------------------------*/

int IGFSCFEncoderEncode(
    IGFSCFENC_INSTANCE_HANDLE           hPublicData,       /**< inout: handle to public data or NULL in case there was no instance created */
    Encoder_State                      *st,                /**< inout: pointer to decoder state */
    int                                 bitCount,          /**< in: offset to the first bit in bitbuffer which should be readed by iisArithDecoderDecode function */
    int                                *sfe,               /**< in: ptr to an array which contain quantized scalefactor energies */
    int                                 indepFlag,         /**< in: if  1 on input the encoder will be forced to reset,
                                                                  if  0 on input the encodder will be forced to encode without a reset */
    int                                 doRealEncoding     /**< in: whether the real encoding is needed, otherwise only the number of bits is used */
)
{
    int ptr[BITBUFSIZE]; /* temporary expanded bit buffer, one bit in each int */
    int i;


    /* insert data: */
    hPublicData->ptrBitIndex = 0;
    hPublicData->bitCount                  = bitCount;
    ari_start_encoding_14bits(&hPublicData->acState); /* start AC encoding */

    /* check if coder needs a reset and do it if necessary */
    if (indepFlag)
    {
        IGFSCFEncoderReset( hPublicData );
    }

    encode_sfe_vector( hPublicData, ptr, hPublicData->t, hPublicData->prev, sfe, hPublicData->scfCountLongBlock, doRealEncoding );

    hPublicData->ptrBitIndex = ari_done_encoding_14bits( ptr, hPublicData->ptrBitIndex, &hPublicData->acState ); /* finish AC encoding */
    hPublicData->bitCount = hPublicData->bitCount + hPublicData->ptrBitIndex;

    /* advancing history: */
    mvi2i(sfe, hPublicData->prev, hPublicData->scfCountLongBlock);
    hPublicData->t++;

    /* copy the bits from the temporary bit buffer, if doRealEncoding is enabled */
    if (doRealEncoding != 0)
    {
        for (i = 0; i < hPublicData->ptrBitIndex; ++i)
        {
            push_next_indice(st, ptr[i], 1);
        }
    }

    /* return next bit offset in the stream */
    return hPublicData->bitCount;
}


/*---------------------------------------------------------------------*
 * IGFSCFEncoderSaveContextState()
 *
 * for a closed loop enc, the ArithEncoder needs to memorize the context
 *---------------------------------------------------------------------*/

void IGFSCFEncoderSaveContextState(
    IGFSCFENC_INSTANCE_HANDLE           hPublicData       /**< inout: handle to public data or NULL in case there was no instance created */
)
{

    hPublicData->Tsave                      = hPublicData->t;

    mvi2i(hPublicData->prev, hPublicData->prevSave, hPublicData->scfCountLongBlock);


    return;
}

/*---------------------------------------------------------------------*
 * IGFSCFEncoderRestoreContextState()
 *
 * for a closed loop enc, the ArithEncoder needs to memorize the context
 *---------------------------------------------------------------------*/
void IGFSCFEncoderRestoreContextState(
    IGFSCFENC_INSTANCE_HANDLE           hPublicData       /**< inout: handle to public data or NULL in case there was no instance created */
)
{


    hPublicData->t                      = hPublicData->Tsave;

    mvi2i(hPublicData->prevSave, hPublicData->prev, hPublicData->scfCountLongBlock);


    return;
}
