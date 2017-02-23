/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"
#include "rom_dec.h"

/*--------------------------------------------------------------------------*/
/*  Function  decode_huff_context                                           */
/*  ~~~~~~~~~~~~~~~~~                                                       */
/*                                                                          */
/*  Context based Huffman decoding for indices of quantized norms           */
/*--------------------------------------------------------------------------*/
/*  const short *hufftab,   (i)    Huffman table                            */
/*  short       *rbits      (i/o)  the number of read bits                  */
/*--------------------------------------------------------------------------*/
static
short decode_huff_context(
    Decoder_State *st,              /* i/o: decoder state structure   */
    const short *hufftab,
    short *rbits
)
{
    while( *hufftab > 0)
    {
        *rbits+=(*hufftab & 0xf);
        hufftab += (*hufftab >> 4) + get_next_indice( st, *hufftab & 0xf );
    }

    return (-*hufftab);
}


/*--------------------------------------------------------------------------*/
/*  Function  hdecnrm                                                       */
/*  ~~~~~~~~~~~~~~~~~                                                       */
/*                                                                          */
/*  Huffman decoding for indices of quantized norms                         */
/*--------------------------------------------------------------------------*/
/*  short       N           (i)    number of norms                          */
/*  short       *index      (o)    indices of quantized norms               */
/*--------------------------------------------------------------------------*/

void hdecnrm(
    Decoder_State *st,              /* i/o: decoder state structure   */
    const short N,
    short *index
)
{
    short i, j, k, n, m;
    short temp;
    short *pidx;

    pidx  = index;

    m = N - 1;
    for (i=0; i<m; i++)
    {
        j = 0;
        k = 0;
        if ( get_next_indice_1( st ) )
        {
            j = 1;
        }

        if ( get_next_indice_1( st ) )
        {
            k = 1;
        }
        n = j * 2 + k;
        j = j * 4;
        temp = 16 + n - j;

        if ( get_next_indice_1( st ) )
        {
            temp = 12 + n + j;

            if ( get_next_indice_1( st ) )
            {
                j = 0;
                if ( get_next_indice_1( st ) )
                {
                    j = 1;
                }

                temp = 8 + n;

                if (j!=0)
                {
                    temp += 12;
                }

                if ( get_next_indice_1( st ) )
                {
                    temp = n;

                    if ( get_next_indice_1( st ) )
                    {
                        temp = n + 4;
                    }

                    if (j!=0)
                    {
                        temp += 24;
                    }
                }
            }
        }

        *pidx++ = temp;
    }

    return;
}

/*--------------------------------------------------------------------------
 *  huff_dec()
 *
 *  Huffman decoding
 *--------------------------------------------------------------------------*/

void huff_dec(
    Decoder_State *st,                /* i/o: decoder state structure                         */
    const short N,                  /* i  : Number of codewords to decode                   */
    const short buffer_len,         /* i  : Number of bits to read                          */
    const short num_lengths,        /* i  : Number of different huffman codeword lengths    */
    const short *thres,             /* i  : Threshold of first codeword of each length      */
    const short *offset,            /* i  : Offset for first codeword                       */
    const short *huff_tab,          /* i  : Huffman table order by codeword lengths         */
    short *index              /* o  : Decoded index                                   */
)
{
    short i, j, k;
    unsigned short val;
    short last_bits = buffer_len;

    val = 0;
    j = 0;
    for (i = 0; i < N; i++)
    {
        last_bits = buffer_len - j;
        val <<= last_bits;
        val &= (1<<buffer_len) - 1; /* 0xFFF; */
        val |= (short)get_next_indice( st, last_bits );

        /* Find codeword length */
        j = num_lengths - 1;
        while (val < thres[j])
        {
            j--;
        }
        k = (val - thres[j]) >> j;
        *index++ = huff_tab[offset[j] + k];
    }

    /* Put back unused bits */
    st->next_bit_pos -= j;

    return;
}


/*--------------------------------------------------------------------------
 *  hdecnrm_context()
 *
 *  Huffman decoding for indices of quantized norms
 *--------------------------------------------------------------------------*/

void hdecnrm_context(
    Decoder_State *st,              /* i/o: decoder state structure   */
    const short N,                /* i  : number of norms           */
    short *index,           /* o  : indices of quantized norms */
    short *n_length         /* o  : decoded stream length     */
)
{
    short i, prevj;

    prevj = index[0] + OFFSET_NORM;
    for( i=1; i < N; i++)
    {
        if( prevj > HTH_NORM )
        {
            /* above */
            index[i] = decode_huff_context( st, hntable, n_length);
            index[i] = 31 - index[i];
        }
        else
        {
            if( prevj < LTH_NORM )
            {
                /* less */
                index[i] = decode_huff_context( st, hntable, n_length );
            }
            else
            {
                /* equal */
                index[i] = decode_huff_context( st, hetable, n_length );
            }
        }
        prevj = index[i];
    }

    return;
}

void hdecnrm_resize(
    Decoder_State *st,              /* i/o: decoder state structure   */
    const short N,
    short *index
)
{
    short i, j, k, m;
    short temp;
    short *pidx;

    pidx  = index;

    m = N - 1;
    for (i=0; i<m; i++)
    {
        j = 0;
        k = 0;

        for( j = 0; j < 11; j++)
        {
            if ( get_next_indice_1( st ) )
            {
                k++;
            }
            else
            {
                break;
            }
        }

        if(k == 11)
        {
            temp = 25;
        }
        else if (k == 10)
        {
            temp = 5;
        }
        else if (k == 9)
        {
            temp = 6;
        }
        else
        {
            if ( get_next_indice_1( st ) )
            {
                temp = 16 + k;
            }
            else
            {
                temp = 15 - k;
            }

        }

        *pidx++ = temp;
    }

    return;
}


/*--------------------------------------------------------------------------
 * hdecnrm_trans()
 *
 * Huffman decoding for indices of quantized norms
 *--------------------------------------------------------------------------*/

void hdecnrm_tran(
    Decoder_State *st,              /* i/o: decoder state structure   */
    const short N,                /* i  : number of norms           */
    short *index            /* o  : indices of quantized norms */
)
{
    short i, j, k, n, m;
    short temp;
    short *pidx;
    short l;

    pidx  = index;

    m = N - 1;
    for (i=0; i<m; i++)
    {
        j = 0;
        k = 0;
        if ( get_next_indice_1(st ) )
        {
            j=1;
        }

        if ( get_next_indice_1(st ) )
        {
            k=1;
        }

        n = k * 2 + j;
        l = k * 4;
        if((j==0 && k==0) || (j==1 && k==0) || (j==1 && k==1))
        {
            temp = 15 + l - n;
        }
        else
        {
            if ( get_next_indice_1(st ) )
            {
                temp = 15+n-l;
            }
            else
            {
                temp = 15+l-n;
                if ( get_next_indice_1(st ) )
                {
                    for(k=0; k<3;)
                    {
                        if(get_next_indice_1(st ))
                        {
                            k++;
                        }
                        else
                        {
                            break;
                        }
                    }

                    if(k==0 || k==3)
                    {
                        temp-=5;
                        if(k==3)
                        {
                            temp--;
                        }
                    }
                    else if(k==1)
                    {
                        temp++;
                    }
                    else
                    {
                        temp +=2;
                        if ( get_next_indice_1(st ) )
                        {
                            temp ++;
                            if ( get_next_indice_1(st ) )
                            {
                                temp++;
                            }
                        }
                    }
                }
            }
        }

        *pidx++ = temp;
    }

    return;
}


