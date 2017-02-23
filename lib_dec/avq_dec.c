/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "prot.h"
#include "options.h"
#include "cnst.h"


/*-----------------------------------------------------------------*
 * AVQ_demuxdec()
 *
 * Demultiplex and decode subvectors using
 * split algebraic vector dequantizer based on RE8 latice.
 *-----------------------------------------------------------------*/

void AVQ_demuxdec(
    Decoder_State *st,        /* i/o: decoder state structure         */
    int   xriq[],     /* o:   decoded subvectors [0..8*Nsv-1] */
    short *nb_bits,   /* i/o: number of allocated bits        */
    const short Nsv,        /* i:   number of subvectors            */
    short nq_out[]    /* i/o: AVQ nq index                    */
)
{
    short i,j, bits, order_v;
    long I[NSV_MAX];
    int nq[NSV_MAX], *kv, code[8];

    kv = xriq;        /* reuse vector to save memory */
    bits = *nb_bits;

    for(i = 0; i <NSV_MAX; i++)
    {
        I[i] = -1;
    }

    for( i=0; i<Nsv; i++ )
    {
        nq[i] = 0;      /* initialization and also forced if the budget is exceeded */

        if( bits > 8 )
        {
            /* read the unary code including the stop bit for nq[i] */
            nq[i] = -1;
            do
            {
                (nq[i])++;

                if (5 * nq[i] + 4 == bits)
                {
                    break;
                }
            }
            while ( get_next_indice_1( st ) );

            if( 5*nq[i]+4 == bits ) /* check the overflow */
            {
                bits++;     /* overflow stop bit */
            }

            /* check for potential bit errors */
            if( nq[i] > NB_SPHERE )
            {
                st->BER_detect = 1;
                set_i( xriq, 0, Nsv*8 );
                set_s( nq_out, 0, Nsv );
                *nb_bits = 0;

                return;
            }

            bits -= (short)nq[i];
            bits--;         /* count the stop bit */

            if( nq[i] > 0 )
            {
                nq[i]++;
            }

            /* read codebook indices (rank I and event. Voronoi index kv) */
            if( nq[i] == 0 ) /* Q0 */
            {
                /* nothing to read */
            }
            else if( nq[i] < 5 )    /* Q2, Q3, Q4 */
            {
                I[i] = get_next_indice( st, (short)(4*nq[i]) );
                bits -= (short)(4*nq[i]);
            }
            else if( nq[i]%2 == 0 )    /* Q4 + Voronoi extensions r=1,2,3,... */
            {
                I[i] = get_next_indice( st, 4*4 );
                bits -= 4*4;
                order_v = (short)(nq[i]/2) - 2;

                for( j=0; j<8; j++ )
                {
                    kv[i*8+j] = get_next_indice( st, order_v );
                }
                bits -= 8*order_v;
            }
            else /* Q3 + Voronoi extensions r=1,2,3,... */
            {
                I[i] = get_next_indice( st, 4*3 );
                bits -= 4*3;
                order_v = (short)(nq[i]/2) - 1;

                for( j=0; j<8; j++ )
                {
                    kv[i*8+j] = get_next_indice( st, order_v );
                }
                bits -= 8*order_v;
            }
        }
    } /* for */

    /* decode all subvectors */
    for( i=0; i<Nsv; i++ )
    {
        /* multi-rate RE8 decoder */
        re8_dec( nq[i], I[i], &kv[8*i], code );

        /* write decoded RE8 vector to decoded subvector #i */
        for( j=0; j<8; j++ )
        {
            xriq[i*8+j] = (short)code[j];
        }
    }

    *nb_bits = bits;

    for( i=0; i<Nsv; i++ )
    {
        nq_out[i] =(short) nq[i];
    }

    return;
}


/*-----------------------------------------------------------------*
 * AVQ_dec_lpc()
 *
 * Demultiplex and decode subvectors for LPC dequantization
 * using split algebraic vector dequantizer
 *-----------------------------------------------------------------*/

void AVQ_dec_lpc(
    const int   indx[],     /* i  : index[] (4 bits per words)      */
    int   nvecq[],    /* o  : vector quantized                */
    const short Nsv         /* i  : number of subvectors (lg=Nsv*8) */
)
{
    int    i, l, n, nq, nk, pos, ival, c[8], kv[8];
    long   I;

    /* last index word */
    pos = Nsv-1;

    for( l=0; l<Nsv; l++ )
    {
        pos += indx[l];
    }

    /* decode all subvectors */
    for( l=Nsv-1; l>=0; l-- )
    {
        nq = indx[l];        /* quantizer number (0,2,3..n) */

        nk = 0;
        n = nq;

        if( nq > 4 )
        {
            nk = (nq-3)>>1;
            n = nq - nk*2;
        }

        /* read n groups of 4-bit for Voronoi index (k[]) */
        for( i=0; i<8; i++)
        {
            kv[i] = 0;
        }

        while( nk-- > 0 )
        {
            ival = (indx[pos--] & 0x0F);
            ival <<= 4;
            ival += (indx[pos--] & 0x0F);

            for( i=7; i>=0; i-- )
            {
                kv[i] <<= 1;
                kv[i] += (ival & 0x01);
                ival >>= 1;
            }
        }

        /* read n groups of 4-bit for base codebook index (I) */
        I = 0;
        while( n-- > 0 )
        {
            I <<= 4;
            I += (indx[pos--] & 0x0F);
        }

        /* multi-rate RE8 decoder */
        re8_dec( nq, I, kv, c );

        /* write decoded RE8 vector */
        for( i=0; i<8; i++ )
        {
            nvecq[(l*8)+i] = c[i];
        }
    }

    return;
}
