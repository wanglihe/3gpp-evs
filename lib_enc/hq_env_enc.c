/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <math.h>
#include "options.h"
#include "prot.h"
#include "rom_com.h"
#include "rom_enc.h"

/*--------------------------------------------------------------------------------------*
 * encode_envelope_indices()
 *
 * Encode envelope indices
 *--------------------------------------------------------------------------------------*/

short encode_envelope_indices(  /* o  : Number of bits if flag_pack=0,0 if flag_pack=1  */
    Encoder_State *st,            /* i/o: encoder state structure                         */
    const short num_sfm,        /* i  : Number of subbands                              */
    const short numnrmibits,    /* i  : Bitrate of fall-back coding mode                */
    short *difidx,        /* i/o: Diff indices/encoded diff indices               */
    short *LCmode,        /* o  : Coding mode if flag_pack=0, i : if flag_pack=1  */
    const short flag_pack,      /* i  : indicator of packing or estimating bits         */
    const short flag_HQ2,       /* i  : indicator of HQ2 core                           */
    const short is_transient    /* i  : transient flag                                  */
)
{
    short bits;
    short prevj;
    short hcode_l;
    short i,j;
    short difidx_flag;
    short index_max, index_min, index_rad;
    short difidx_org[NB_SFM];       /* lenght of this buffer is max(BANDS_MAX,NB_SFM) */
    short m, r;
    short v, k;

    set_s( difidx_org, 0, NB_SFM );
    difidx_flag = 0;

    /*------------------------------------------------------------------*
     * Check Huffman encoding for QNorm indices
     *------------------------------------------------------------------*/

    /* LC mode index is changed to synchronize LR-MDCT signaling    */
    /* LC mode 0 = Context based coding                             */
    /* LC mode 1 = resized huffman coding                           */
    /* LC mode 2 = normal Huffman Coding                            */
    /* LC mode 3 = bit packing                                      */
    if ( flag_pack == 0 )
    {
        if(is_transient && flag_HQ2 == LOW_RATE_HQ_CORE_TRAN)
        {
            bits = 0;
            index_max = 0;
            index_min = 31;
            for( i = 0; i< num_sfm; i++ )
            {
                if( difidx[i] > index_max )
                {
                    index_max = difidx[i];
                }
                if( difidx[i] < index_min )
                {
                    index_min = difidx[i];
                }
            }
            if(index_min > 10 && index_max < 22)
            {
                for( i = 1; i < num_sfm; i++ )
                {
                    j = difidx[i];
                    bits += huffsizn_tran[j];
                }
            }
            hcode_l= 0;
            *LCmode = 0;
            prevj = difidx[0] + OFFSET_NORM;
            /* LC mode 0 = Context based coding        */
            for( i = 1; i < num_sfm; i++ )
            {
                j = difidx[i];
                if( prevj>HTH_NORM )
                {
                    /* above */
                    hcode_l += huffsizn_n[31-j];
                }
                else
                {
                    if( prevj<LTH_NORM )
                    {
                        /* less */
                        hcode_l += huffsizn_n[j];
                    }
                    else
                    {
                        /* equal */
                        hcode_l += huffsizn_e[j];
                    }
                }
                prevj = j;
            }
            if( hcode_l >= bits && bits !=0)
            {
                /* LC mode 1 Transient Huffman Coding   */
                *LCmode = 1;
                hcode_l = bits;
            }
        }
        else
        {
            /* Check bits if LC mode == 3 -> Check bits if LC mode == 0 */
            hcode_l= 0;
            prevj = difidx[0] + OFFSET_NORM;
            for( i = 1; i < num_sfm; i++ )
            {
                j = difidx[i];
                if( prevj>HTH_NORM )
                {
                    /* above */
                    hcode_l += huffsizn_n[31-j];
                }
                else
                {
                    if( prevj<LTH_NORM )
                    {
                        /* less */
                        hcode_l += huffsizn_n[j];
                    }
                    else
                    {
                        /* equal */
                        hcode_l += huffsizn_e[j];
                    }
                }
                prevj = j;
            }

            *LCmode = 0;

            /* LR-MDCT core doesn't have coding mode 2 and 3 */
            if( flag_HQ2 == NORMAL_HQ_CORE )
            {
                /* Check bits if LC mode == 1 -> Check bits if LC mode == 2 */
                bits = 0;
                for( i = 1; i < num_sfm; i++ )
                {
                    j = difidx[i];
                    bits += huffsizn[j];
                }

                /*------------------------------------------------------------------------------*
                 * comparing bit expenses of coding mode 2 with that of the optimal coding mode
                 *------------------------------------------------------------------------------*/

                if( hcode_l > bits )
                {
                    *LCmode = 2;
                    hcode_l = bits;
                }
            }

            /* Check bits if LC mode == 2 -> Check bits if LC mode == 1  */
            bits = 0;
            index_max = 0;
            index_min = 31;
            for( i = 1; i < num_sfm; i++ )
            {
                difidx_org[i] = difidx[i];
            }

            difidx_flag = 0;
            for( i = 2; i < num_sfm; i++ )
            {
                if( difidx_org[i-1] > 17 )
                {
                    difidx[i] = difidx_org[i] + min((difidx_org[i-1]-17),3);
                    if( difidx[i] > 31 )
                    {
                        difidx_flag = 1;
                        break;
                    }
                }

                if( difidx_org[i-1] < 13 )
                {
                    difidx[i] = difidx_org[i] + max((difidx_org[i-1]-13),-3);
                    if( difidx[i] < 0 )
                    {
                        difidx_flag = 1;
                        break;
                    }
                }
            }

            index_rad = 0;
            if( difidx_flag != 1 )
            {
                for( i = 1; i< num_sfm; i++ )
                {
                    if( difidx[i] > index_max )
                    {
                        index_max = difidx[i];
                    }

                    if( difidx[i] < index_min )
                    {
                        index_min = difidx[i];
                    }
                }

                index_rad = max((15 - index_min),(index_max - 15));

                if( index_rad <= HUFF_THR )
                {
                    for( i = 1; i < num_sfm; i++ )
                    {
                        j = difidx[i];
                        bits += resize_huffsizn[j];
                    }

                    /*------------------------------------------------------------------*
                     * comparing bit expenses of coding mode 1 with that of coding mode 0
                     *------------------------------------------------------------------*/

                    if( hcode_l > bits )
                    {
                        hcode_l = bits;
                        *LCmode = 1;
                    }
                }
            }

            /* LR-MDCT core doesn't have coding mode 2 and 3 */
            if( flag_HQ2 == NORMAL_HQ_CORE )
            {
                /*------------------------------------------------------------------------------*
                 * comparing bit expenses of coding mode 3 with that of the optimal coding mode
                 *------------------------------------------------------------------------------*/

                if( hcode_l >= numnrmibits )
                {
                    hcode_l = numnrmibits;
                    *LCmode = 3;
                }
            }

            if( (*LCmode != 1 && flag_HQ2 == NORMAL_HQ_CORE ) || flag_HQ2 == LOW_RATE_HQ_CORE )
            {
                for(i = 2; i< num_sfm; i++)
                {
                    difidx[i] = difidx_org[i];
                }
            }
        }
    }
    else
    {
        if(flag_HQ2 == LOW_RATE_HQ_CORE_TRAN || flag_HQ2 == LOW_RATE_HQ_CORE)
        {
            push_indice( st, IND_HQ2_DENG_HMODE, *LCmode, BITS_DE_HMODE);
            push_indice( st, IND_HQ2_DIFF_ENERGY, difidx[0], BITS_DE_FCOMP);
        }
        else
        {
            push_indice( st, IND_LC_MODE, *LCmode, 2 );
            push_indice( st, IND_YNRM, difidx[0], NORM0_BITS );
        }

        if(is_transient && flag_HQ2 == LOW_RATE_HQ_CORE_TRAN)
        {
            hcode_l = 0;
            if ( *LCmode == 1 )
            {
                /* LC mode 0 Transient Huffman Coding   */
                for( i = 1; i < num_sfm; i++ )
                {
                    j = difidx[i];
                    m = huffnorm_tran[j];
                    r = huffsizn_tran[j];
                    v = 0;

                    /* Bit reverse */
                    for( k = 0; k < r; k++ )
                    {
                        v <<= 1;
                        v |= m & 1;
                        m >>= 1;
                    }

                    push_indice(st,IND_HQ2_DIFF_ENERGY, v, r);
                }
            }
            else
            {
                /* LC mode 1 context based Coding   */
                prevj = difidx[0] + OFFSET_NORM;
                for( i = 1; i < num_sfm; i++ )
                {
                    j = difidx[i];

                    if( prevj > HTH_NORM )
                    {
                        /* above */
                        r = huffsizn_n[31-j];
                        m = huffnorm_n[31-j];
                    }
                    else
                    {
                        if( prevj<LTH_NORM )
                        {
                            /* less */
                            r = huffsizn_n[j];
                            m = huffnorm_n[j];
                        }
                        else
                        {
                            /* equal */
                            r = huffsizn_e[j];
                            m = huffnorm_e[j];
                        }
                    }
                    push_indice(st,IND_HQ2_DIFF_ENERGY, m, r);
                    prevj = j;
                }
            }
        }
        else
        {
            hcode_l = 0;
            if ( *LCmode == 0 )
            {
                /* LC mode 3 -> LC mode 0 */
                prevj = difidx[0] + OFFSET_NORM;
                for( i = 1; i < num_sfm; i++ )
                {
                    j = difidx[i];

                    if( prevj > HTH_NORM )
                    {
                        /* above */
                        r = huffsizn_n[31-j];
                        m = huffnorm_n[31-j];
                    }
                    else
                    {
                        if( prevj<LTH_NORM )
                        {
                            /* less */
                            r = huffsizn_n[j];
                            m = huffnorm_n[j];
                        }
                        else
                        {
                            /* equal */
                            r = huffsizn_e[j];
                            m = huffnorm_e[j];
                        }
                    }

                    if( flag_HQ2 == LOW_RATE_HQ_CORE )
                    {
                        push_indice( st,IND_HQ2_DIFF_ENERGY, m, r);
                    }
                    else
                    {
                        push_indice( st, IND_YNRM, m, r );
                    }

                    prevj = j;
                }
            }
            else if( *LCmode == 1 )
            {
                if ( flag_HQ2 == 1 )
                {
                    index_max = 0;
                    index_min = 31;
                    for(i = 1; i< num_sfm; i++)
                    {
                        difidx_org[i] = difidx[i];
                    }

                    for(i = 2; i< num_sfm; i++)
                    {
                        if(difidx_org[i-1] > 17)
                        {
                            difidx[i] = difidx_org[i] + min((difidx_org[i-1]-17),3);
                            if(difidx[i] > 31)
                            {
                                difidx_flag = 1;
                                break;
                            }
                        }

                        if(difidx_org[i-1] < 13)
                        {
                            difidx[i] = difidx_org[i] + max((difidx_org[i-1]-13),-3);
                            if(difidx[i] < 0)
                            {
                                difidx_flag = 1;
                                break;
                            }
                        }
                    }

                    if( difidx_flag != 1 )
                    {
                        for(i = 1; i< num_sfm; i++)
                        {
                            if(difidx[i]>index_max)
                            {
                                index_max = difidx[i];
                            }

                            if(difidx[i]<index_min)
                            {
                                index_min = difidx[i];
                            }
                        }

                        index_rad = max((15 - index_min),(index_max - 15));

                        if(index_rad <= HUFF_THR)
                        {
                            for (i = 1; i < num_sfm; i++)
                            {
                                j = difidx[i];
                            }
                        }
                    }
                }

                /* LC mode 2 -> LC mode 1 */
                for( i = 1; i < num_sfm; i++ )
                {
                    j = difidx[i];

                    m = resize_huffnorm[j];
                    r = resize_huffsizn[j];
                    v = 0;

                    /* Bit reverse */
                    for( k = 0; k < r; k++ )
                    {
                        v <<= 1;
                        v |= m & 1;
                        m >>= 1;
                    }

                    if ( flag_HQ2 )
                    {
                        push_indice( st,IND_HQ2_DIFF_ENERGY, v, r);
                    }
                    else
                    {
                        push_indice( st, IND_YNRM, v, r );
                    }
                }
            }
            else if( *LCmode == 2 )
            {
                /* LC mode 1 -> LC mode 2 */
                for( i = 1; i < num_sfm; i++ )
                {
                    j = difidx[i];

                    m = huffnorm[j];
                    r = huffsizn[j];

                    push_indice( st, IND_YNRM, m, r );
                }
            }
            else
            {
                for( i = 1; i < num_sfm; i++ )
                {
                    push_indice( st, IND_YNRM, difidx[i], NORMI_BITS );
                }
            }
        }
    }

    return hcode_l;
}

/*--------------------------------------------------------------------------*
 * diff_envelope_coding()
 *
 * Create differential code of norm indices
 *--------------------------------------------------------------------------*/

void diff_envelope_coding(
    const short is_transient,       /* i  : transient indicator                 */
    const short num_env_bands,      /* i  : number of envelope bands to code    */
    const short start_norm,         /* i  : start of envelope coding            */
    short *ynrm,              /* i/o: quantization indices for norms      */
    short *normqlg2,          /* i/o: quantized norms                     */
    short *difidx             /* o  : differential code                   */
)
{
    short i;
    short idxbuf[NB_SFM];
    short normbuf[NB_SFM];

    /* Differential coding for indices of quantized norms */
    if( is_transient )
    {
        /* Reorder quantization indices and quantized norms */
        reordernorm( ynrm, normqlg2, idxbuf, normbuf, num_env_bands );
        diffcod( num_env_bands, idxbuf, &difidx[1] );
        difidx[0] = idxbuf[0];
        recovernorm( idxbuf, ynrm, normqlg2, num_env_bands );
    }
    else
    {
        diffcod( num_env_bands, &ynrm[start_norm], &difidx[1] );
        difidx[0] = ynrm[start_norm];

        for( i = start_norm; i < start_norm + num_env_bands; i++ )
        {
            normqlg2[i] = dicnlg2[ynrm[i]];
        }
    }

    return;
}
