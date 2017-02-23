/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "prot.h"
#include "rom_com.h"

/*-------------------------------------------------------------------*
 * Function AVQ_cod()                                                *
 *                                                                   *
 * Split algebraic vector quantizer (AVQ) based on RE8 latice        *
 *-------------------------------------------------------------------*/

float AVQ_cod(              /* o:   comfort noise gain factor        */
    const float xri[],      /* i:   vector to quantize               */
    int   xriq[],     /* o:   quantized normalized vector (assuming the bit budget is enough) */
    const short nb_bits,    /* i:   number of allocated bits         */
    const short Nsv         /* i:   number of subvectors (lg=Nsv*8)  */
)
{
    short  i, j, iter;
    int    c[8];
    float  gain_inv, x1[8], ener, tmp, nbits, nbits_max, fac, offset;
    float  ebits[NSV_MAX];

    /* find energy of each subvector in log domain (scaled for bits estimation) */
    for( i=0; i<Nsv; i++ )
    {
        ener = 2.0f; /* to set ebits >= 0 */
        for( j=0; j<8; j++ )
        {
            x1[j] = xri[i*8+j];
            ener += x1[j]*x1[j];
        }

        /* estimated bit consumption when gain=1 */
        ebits[i] = 5.0f * FAC_LOG2 * (float)log10( ener*0.5f );
    }

    /* estimate gain according to number of bits allowed */
    fac = 128.0f;      /* start at the middle (offset range = 0 to 255.75) */
    offset = 0.0f;
    nbits_max = 0.95f * ((float)(nb_bits - Nsv));

    /* tree search with 10 iterations : offset with step of 0.25 bits (0.3 dB) */
    for( iter=0; iter<10; iter++ )
    {
        offset += fac;
        /* calculate the required number of bits */
        nbits = 0.0;
        for( i=0; i<Nsv; i++ )
        {
            tmp = ebits[i] - offset;
            if( tmp < 0.0 )
            {
                tmp = 0.0;
            }
            nbits += tmp;
        }
        /* decrease gain when no overflow occurs */
        if( nbits <= nbits_max )
        {
            offset -= fac;
        }
        fac *= 0.5;
    }

    /* estimated gain (when offset=0, estimated gain=1) */
    gain_inv = 1.0f / (float)pow(10.0f, (float) (offset / (2.0f*5.0f*FAC_LOG2)) );

    /* quantize all subvector using estimated gain */
    for( i=0; i<Nsv; i++ )
    {
        for( j=0; j<8; j++ )
        {
            x1[j] = xri[i*8+j] * gain_inv;
        }

        re8_PPV( x1, c );
        for( j=0; j<8; j++ )
        {
            xriq[i*8+j] = c[j];
        }
    }

    fac = 0;

    /* round bit allocations and save */
    for( i=0; i<Nsv; i++ )
    {
        xriq[(Nsv*8)+i] = (int) floor( ebits[i]*128.0f );
    }

    return( fac );

}


/*-----------------------------------------------------------------*
 * AVQ_encmux()
 *
 * Encode subvectors and write indexes into the bitstream
 *-----------------------------------------------------------------*/

void AVQ_encmux(
    Encoder_State *st,        /* i/o: encoder state structure      */
    const short extl,       /* i  : extension layer                                 */
    int   xriq[],     /* i/o: rounded subvectors [0..8*Nsv-1] followed
                                        by rounded bit allocations [8*Nsv..8*Nsv+Nsv-1] */
    short *nb_bits,   /* i/o: number of allocated bits                        */
    const short Nsv,        /* i:   number of subvectors                            */
    short nq_out[]    /* o  : AVQ nq index                                    */
)
{
    short i, j=0, bits, pos, pos_max, overflow;
    short sort_idx[NSV_MAX];
    int   *t, nq[NSV_MAX], kv[NSV_MAX*8];
    long  I[NSV_MAX];
    short nq_ind, i_ind, kv_ind;

    if( extl == SWB_BWE_HIGHRATE || extl == FB_BWE_HIGHRATE )
    {
        nq_ind = IND_NQ2;
        i_ind = IND_I2;
        kv_ind = IND_KV2;
    }
    else
    {
        nq_ind = IND_NQ;
        i_ind = IND_I;
        kv_ind = IND_KV;
    }

    for(i = 0; i <NSV_MAX; i++)
    {
        I[i] = -1;
    }

    /*-----------------------------------------------------------------
     * Encode subvectors and fix possible overflows in total bit budget,
     * i.e. find for each subvector a codebook index nq (nq=0,2,3,4,...,NSV_MAX),
     * a base codebook index (I), and a Voronoi index (kv)
     *-----------------------------------------------------------------*/

    /* sort subvectors by estimated bit allocations in decreasing order */
    t = kv;   /* reuse vector to save memory */
    for( i=0; i<Nsv; i++ )
    {
        t[i] = xriq[8*Nsv+i];
    }

    for( i=0; i<Nsv; i++ )
    {
        bits = (short)t[0];
        pos = 0;
        for( j=1; j<Nsv; j++ )
        {
            if( t[j] > bits )
            {
                bits = (short)t[j];
                pos = j;
            }
        }
        sort_idx[i] = pos;
        t[pos] = -1;
    }

    /* compute multi-rate indices and avoid bit budget overflow */
    pos_max = 0;
    bits = 0;
    for( i=0; i<Nsv; i++ )
    {
        /* find vector to quantize (criteria: nb of estimated bits) */
        pos = sort_idx[i];

        /* compute multi-rate index of rounded subvector (nq,I,kv[]) */
        re8_cod( &xriq[pos*8], &nq[pos], &I[pos], &kv[8*pos] );

        if( nq[pos] > 0 )
        {
            j = pos_max;
            if( pos > j )
            {
                j = pos;
            }

            /* compute (number of bits -1) to describe Q #nq */
            if( nq[pos] >= 2 )
            {
                overflow = (short)(nq[pos]*5-1);
            }
            else
            {
                overflow = 0;
            }

            /* check for overflow and compute number of bits-1 (n) */
            if( (bits+overflow+j) > *nb_bits )
            {
                /* if budget overflow */
                for( j=pos*8; j<(pos*8)+8; j++ )
                {
                    xriq[j] = 0;
                }
                nq[pos] = 0; /* force Q0 */
            }
            else
            {
                bits += overflow;
                pos_max = j; /* update index of the last described subvector */
            }
        }
    }

    /* write indexes to the bitstream */
    /* ============================== */

    bits = *nb_bits;
    overflow = 0;

    for( i=0; i<Nsv; i++ )
    {
        if( 5*nq[i]-1 == bits ) /* check the overflow */
        {
            overflow = 1;
        }

        if( bits > 8 )
        {
            /* write the unary code for nq[i] */
            j = (short)(nq[i] - 1);
            if ( nq[i] > 0 )
            {
                /* write the unary code */
                while ( j > 16 )
                {
                    push_indice( st, nq_ind, 65535, 16 );
                    bits -= 16;
                    j -= 16;
                }

                if ( j > 0 )
                {
                    push_indice( st, nq_ind, (1<<j)-1, j );
                    bits -= j;
                }
            }

            if ( !overflow )
            {
                /* write the stop bit */
                push_indice( st, nq_ind, 0, 1 );
                bits--;
            }

            /* write codebook indices (rank I and event. Voronoi index kv) */
            if( nq[i] == 0 )    /* Q0 */
            {
                /* nothing to write */
            }
            else if( nq[i] < 5 )    /* Q2, Q3, Q4 */
            {
                push_indice( st, i_ind, I[i], (short)(4*nq[i]) );
                bits -= (short)(4*nq[i]);
            }
            else if( nq[i]%2 == 0 )    /* Q4 + Voronoi extensions r=1,2,3,... */
            {
                push_indice( st, i_ind, I[i], 4*4 );
                bits -= 4*4;
                pos = (short)(nq[i]/2 - 2);  /* Voronoi order determination */
                for( j=0; j<8; j++ )
                {
                    push_indice( st, kv_ind, kv[i*8+j], pos );
                }

                bits -= 8*pos;
            }
            else    /* Q3 + Voronoi extensions r=1,2,3,... */
            {
                push_indice( st, i_ind, I[i], 4*3 );
                bits -= 4*3;

                pos = (short)(nq[i]/2 - 1);  /* Voronoi order determination */
                for( j=0; j<8; j++ )
                {
                    push_indice( st, kv_ind, kv[i*8+j], pos );
                }

                bits -= 8*pos;
            }
        }
    } /* for */

    *nb_bits = bits;

    for( i=0; i<Nsv; i++ )
    {
        nq_out[i] = (short) nq[i];
    }

    return;
}


/*-------------------------------------------------------------------*
 * Function AVQ_cod_lpc()                                            *
 *                                                                   *
 * Split algebraic vector quantizer (AVQ) for LPC quantization       *
 *-------------------------------------------------------------------*/

void AVQ_cod_lpc(
    const float nvec[],     /* i:   vector to quantize              */
    int   nvecq[],    /* o:   quantized normalized vector (assuming the bit budget is enough) */
    int   *indx,      /* o:   index[] (4 bits per words)      */
    const short Nsv         /* i:   number of subvectors (lg=Nsv*8) */
)
{
    int    i, l, n, nq, nk, pos, ival, c[8], kv[8];
    float  x1[8];
    long   I;

    /* quantize all subvector using estimated gain */
    pos = Nsv;

    for( l=0; l<Nsv; l++ )
    {
        for( i=0; i<8; i++ )
        {
            x1[i] = nvec[l*8+i];
        }

        re8_PPV( x1, c );

        re8_cod( c, &nq, &I, kv );

        for( i=0; i<8; i++ )
        {
            nvecq[l*8+i] = c[i];
        }

        indx[l] = nq;      /* index[0..Nsv-1] = quantizer number (0,2,3,4...) */

        nk = 0;
        n = nq;

        if( nq > 4 )
        {
            nk = (nq-3)>>1;
            n = nq - nk*2;
        }

        /* write n groups of 4-bit for base codebook index (I) */
        while( n-- > 0 )
        {
            indx[pos++] = (I & 0x0F);
            I >>= 4;
        }

        /* write n groups of 4-bit for Voronoi index (k[]) */
        while( nk-- > 0 )
        {
            ival = 0;

            for( i=0; i<8; i++ )
            {
                ival <<= 1;
                ival += (kv[i] & 0x01);
                kv[i] >>= 1;
            }
            indx[pos++] = (ival & 0x0F);
            ival >>= 4;
            indx[pos++] = (ival & 0x0F);
        }
    }

    return;
}
