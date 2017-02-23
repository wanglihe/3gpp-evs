/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <math.h>
#include "options.h"
#include "prot.h"
#include "rom_com.h"

/*------------------------------------------------------------------------*
 * decode_envelope_indices()
 *
 * Decode envelope indices
 *------------------------------------------------------------------------*/

short decode_envelope_indices(  /* o  : Number of bits                    */
    Decoder_State *st,            /* i/o: decoder state structure           */
    const short start_norm,     /* i  : starting band index               */
    const short num_sfm,        /* i  : Number of subbands                */
    const short numnrmibits,    /* i  : Bitrate of fall-back coding mode  */
    short *difidx,        /* o  : Diff indices/encoded diff indices */
    const short flag_HQ2,       /* i  : indicator of HQ2 core             */
    const short is_transient    /* i  : transient flag                    */
)
{
    short hcode_l;
    short i;
    short LCmode;

    if( flag_HQ2 == LOW_RATE_HQ_CORE || flag_HQ2 == LOW_RATE_HQ_CORE_TRAN)
    {
        LCmode = (short)get_next_indice( st, BITS_DE_HMODE);
        difidx[start_norm] = (short)get_next_indice( st, BITS_DE_FCOMP);
    }
    else
    {
        LCmode = (short)get_next_indice( st, 2 );
        difidx[start_norm] = (short)get_next_indice( st, NORM0_BITS );
    }

    if(is_transient && flag_HQ2 == LOW_RATE_HQ_CORE_TRAN)
    {
        hcode_l = 0;
        if(LCmode == 1 )
        {
            hdecnrm_tran(st, num_sfm, &difidx[start_norm + 1] );
            for( i = start_norm + 1; i < start_norm + num_sfm; i++ )
            {
                hcode_l += huffsizn_tran[difidx[i]];
            }
        }
        else
        {
            hdecnrm_context(st,num_sfm, &difidx[start_norm], &hcode_l);
        }
    }
    else
    {
        hcode_l = 0;
        if ( LCmode == 0 )
        {
            hdecnrm_context( st, num_sfm, &difidx[start_norm], &hcode_l );
        }
        else if( LCmode == 1 )
        {
            hdecnrm_resize( st, num_sfm, &difidx[start_norm + 1] );

            for( i = start_norm + 1; i < start_norm + num_sfm; i++ )
            {
                hcode_l += resize_huffsizn[difidx[i]];
            }

            for( i = start_norm + 2; i< start_norm + num_sfm; i++ )
            {
                if( difidx[i-1]>17 )
                {
                    difidx[i] = difidx[i] - min(difidx[i-1]-17,3);
                }

                if( difidx[i-1]<13 )
                {
                    difidx[i] = difidx[i] - max(difidx[i-1]-13,-3);
                }
            }
        }
        else if ( LCmode == 2 )
        {
            hdecnrm( st, num_sfm, &difidx[start_norm + 1] );
            for( i = start_norm + 1; i < start_norm + num_sfm; i++ )
            {
                hcode_l += huffsizn[difidx[i]];
            }
        }
        else
        {
            for( i = start_norm + 1; i < start_norm + num_sfm; i++ )
            {
                difidx[i] = (short)get_next_indice( st, NORMI_BITS );

            }
            hcode_l = numnrmibits;
        }
    }

    return hcode_l;
}

/*------------------------------------------------------------------------*
 * dequantize_norms()
 *
 * De-quantization of norms
 *------------------------------------------------------------------------*/

void dequantize_norms(         /* o  : Number of bits                    */
    Decoder_State *st,            /* i/o: decoder state structure           */
    const short start_norm,     /* i  : First SDE encoded norm            */
    const short num_sfm,        /* i  : Number of norms                   */
    const short is_transient,   /* i  : Transient flag                    */
    short *ynrm,                /* o  : Decoded norm indices              */
    short *normqlg2             /* o  : Log2 of decoded norms             */
)
{
    short i,j,k;
    short idxbuf[NB_SFM];

    /* First sub-frame */
    normqlg2[start_norm] = dicnlg2[ynrm[start_norm]];

    /* Other sub-frames */
    if( is_transient )
    {
        /* Recover quantization indices and quantized norms */
        idxbuf[0] = ynrm[0];
        for( i = 1; i < num_sfm; i++ )
        {
            idxbuf[i] = ynrm[i] + idxbuf[i-1] - 15;
            /* safety check in case of bit errors */
            if ( idxbuf[i] < 0 || idxbuf[i] > 39 )
            {
                idxbuf[i] = 39;
                st->BER_detect = 1;
            }
        }

        recovernorm( idxbuf, ynrm, normqlg2, num_sfm );
    }
    else
    {
        for (i = start_norm + 1; i < start_norm + num_sfm; i++)
        {
            j = i - 1;
            k = ynrm[j] - 15;
            ynrm[i] = ynrm[i] + k;
            /* safety check in case of bit errors */
            if ( ynrm[i] < 0 || ynrm[i] > 39 )
            {
                ynrm[i] = 39;
                st->BER_detect = 1;
            }
            normqlg2[i] = dicnlg2[ynrm[i]];
        }
    }

    return;
}
