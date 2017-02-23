/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"
#include "basop_util.h"
#include "basop_proto_func.h"

/*---------------------------------------------------------------------*
 * tcq_core_LR_enc()
 *
 * Main Generic Audio Decoder Routine for LR-MDCT
 *---------------------------------------------------------------------*/

void tcq_core_LR_dec(
    Decoder_State *st,
    int	  *inp_vector,
    const short bit_budget,
    const short BANDS,
    const short *band_start,
    const short *band_width,
    Word32 *Rk_fx,
    int   *npulses,
    short *k_sort,
    const short *p2a_flags,
    const short p2a_bands,
    const short *last_bitalloc,
    const short input_frame,
    const short adjustFlag,
    const short *is_transient
)
{
    short i, j, k;
    float Rk_sort[NB_SFM];
    short flag_wbnb = 0;
    short USQ_TCQ[NB_SFM];              /* TCQ is selected by default*/
    short nb_bytes, pulsesnum, nz;
    int positions[L_FRAME32k];
    short k_num[2];
    ARCODEC ardec, *pardec;
    BITSTREAM bs, *pbs;

    int nzbands = 0;
    int lsbtcq_bits = TCQ_AMP;
    int tcq_arbits  = 2;

    /* LSB TCQ variables*/
    short dpath[280];
    short bcount = 0;
    float mbuffer[560];

    Word32 leftbits = 0;
    Word32 sepbits = 0;
    Word32 divider = 0;

    /*Word32 Rk_fx[NB_SFM];*/      /* Q16 */
    Word32 Rk_sort_fx[NB_SFM]; /* Q16 */
    Word32 bsub_fx = 0;

    Word16 nzb = 0;
    Word32 delta_fx;
    Word32 surplus_fx;
    Word32 bit_surplus_fx[2];

    /* initialization */
    set_s(dpath, 0, 280);
    set_f(mbuffer, 0.f, 560);
    set_f( Rk_sort, 0.f, NB_SFM );
    set_s( USQ_TCQ, 0, NB_SFM );
    set_i( positions, 0, L_FRAME32k );

    if( input_frame <= L_FRAME16k  && adjustFlag == 0 && *is_transient == 0 )
    {
        flag_wbnb = 1;
        lsbtcq_bits = 0;
        tcq_arbits  = 0;
    }

    pardec = &ardec;
    pbs = &bs;
    pbs->curPos = 7;
    pbs->numbits = 0;
    pbs->numByte = 0;

    /* Bits distribution analysis*/
    for( i = 0; i < BANDS; i++ )
    {
        if( L_sub( ar_div(Rk_fx[i], band_width[i]), 49152) >= 0)
        {
            /* USQ used for high importance bands*/
            USQ_TCQ[i] = 1;
        }
        else
        {
            /* TCQ used for usual bands*/
            USQ_TCQ[i] = 0;
        }
        if( Rk_fx[i] > 0.0f )
        {
            nzbands++;
        }
    }

    for( j = 0; j < BANDS; j++ )
    {
        if( Rk_fx[j] > 0.0f )
        {
            nzb++;
        }
    }

    bsub_fx = L_shl(L_add(tcq_arbits, lsbtcq_bits), 16);
    IF( bsub_fx > 0)
    {
        bsub_fx = L_add( bsub_fx,  2048);
    }
    for( j = BANDS - 1; j >= 0; j-- )
    {
        if( Rk_fx[j] > 0 )
        {
            Rk_fx[j] = L_sub(Rk_fx[j], ar_div(bsub_fx, nzb));
            if( Rk_fx[j] < 0)
            {
                bsub_fx = L_sub(bsub_fx, L_add(ar_div(bsub_fx, nzb), Rk_fx[j]));
                Rk_fx[j] = 0;
            }
            else
            {
                bsub_fx = L_sub(bsub_fx, ar_div(bsub_fx, nzb));
            }
            nzb = sub(nzb, 1);
        }
    }

    srt_vec_ind_fx( Rk_fx, Rk_sort_fx, k_sort, BANDS);

    /*read the bits*/
    nb_bytes = bit_budget >> 3;
    k = bit_budget - (nb_bytes << 3);
    for( i = 0; i < nb_bytes; i++ )
    {
        pbs->buf[i] = (unsigned char)get_next_indice(st, 8);
    }

    if( k > 0 )
    {
        pbs->buf[nb_bytes] = (unsigned char)get_next_indice(st, (short)k);
        pbs->buf[nb_bytes] <<= (8 - k);
        i++;
        nb_bytes++;
    }
    /* set two more bytes, which are used to flush the arithmetic coder, to 0
       -> this avoids reading of uninitialized memory */
    nb_bytes = min(nb_bytes + 2, MAX_SIZEBUF_PBITSTREAM);
    for( ; i < nb_bytes; i++ )
    {
        pbs->buf[i] = 0;
    }

    pbs->maxBytes = nb_bytes;

    ar_decoder_start( pardec, pbs );

    delta_fx = 0;
    surplus_fx = 0;

    if( input_frame <= L_FRAME16k  && adjustFlag == 0 && *is_transient == 0 )
    {
        surplus_fx = -131072;

        bit_allocation_second_fx( Rk_fx, Rk_sort_fx, BANDS, band_width,
                                  k_sort, k_num, p2a_flags, p2a_bands, last_bitalloc, input_frame);

        nzbands = 0;
        for ( j = 0; j < BANDS; j++ )
        {
            if ( sub(j, k_num[0]) == 0 || sub(j, k_num[1]) == 0)
            {
                sepbits = L_add( sepbits, Rk_fx[k_sort[j]]);
            }
            else
            {
                leftbits = L_add( leftbits, Rk_fx[k_sort[j]]);
                if( Rk_fx[k_sort[j]] > 0 )
                {
                    nzbands = add(nzbands, 1);
                }
            }
        }

        for( k = 0; k < BANDS; k++ )
        {
            if( k != k_num[0] && k != k_num[1])
            {
                if (Rk_fx[k_sort[k]] > 0 && USQ_TCQ[k_sort[k]] == 0)
                {
                    /* When number of bits per band is less than
                       arithmetic bits overhead, this band is not encoded.
                       It may happens when the actual number of bits per
                       band is near same to estimated number of bits, for
                       most bands (very unprobable but possible) */
                    if( L_add( Rk_fx[k_sort[k]], delta_fx) < 0 )
                    {
                        pulsesnum = 0;
                        for( i = 0; i < band_width[k_sort[k]]; i++ )
                        {
                            inp_vector[band_start[k_sort[k]] + i] = 0;
                        }
                        if( surplus_fx != 0 )
                        {
                            surplus_fx = L_add( Rk_fx[k_sort[k]], surplus_fx);
                            surplus_fx = L_add( delta_fx, surplus_fx);
                        }
                    }
                    else
                    {
                        /*get number of pulses */
                        pulsesnum = GetScale_fx( band_width[k_sort[k]],
                                                 L_add( Rk_fx[k_sort[k]], delta_fx),
                                                 &surplus_fx );

                        leftbits = L_sub( leftbits, L_add( Rk_fx[k_sort[k]], delta_fx) );

                        decode_position_ari_fx( pardec, band_width[k_sort[k]], pulsesnum, &nz, &positions[band_start[k_sort[k]]] );
                        decode_mangitude_tcq_fx( pardec, band_width[k_sort[k]], pulsesnum, nz, &positions[band_start[k_sort[k]]], &inp_vector[band_start[k_sort[k]]], &surplus_fx );
                        decode_signs_fx( pardec, band_width[k_sort[k]], &inp_vector[band_start[k_sort[k]]] );

                    }
                    nzbands--;
                }
                else if (Rk_fx[k_sort[k]] > 0 && USQ_TCQ[k_sort[k]] == 1)
                {
                    /* When number of bits per band is less than
                       arithmetic bits overhead, this band is not encoded.
                       It may happens when the actual number of bits per
                       band is near same to estimated number of bits, for
                       most bands (very unprobable but possible) */
                    if( L_add( Rk_fx[k_sort[k]], delta_fx) < 0 )
                    {
                        pulsesnum = 0;
                        for( i = 0; i < band_width[k_sort[k]]; i++ )
                        {
                            inp_vector[band_start[k_sort[k]] + i] = 0;
                        }
                        if( surplus_fx != 0 )
                        {
                            surplus_fx = L_add( Rk_fx[k_sort[k]], surplus_fx);
                            surplus_fx = L_add( delta_fx, surplus_fx);
                        }
                    }
                    else
                    {

                        pulsesnum = GetScale_fx(band_width[k_sort[k]],
                                                L_add( Rk_fx[k_sort[k]], delta_fx),
                                                &surplus_fx);

                        leftbits = L_sub( leftbits, L_add( Rk_fx[k_sort[k]], delta_fx) );

                        decode_position_ari_fx( pardec, band_width[k_sort[k]], pulsesnum, &nz, &positions[band_start[k_sort[k]]] );
                        decode_magnitude_usq_fx( pardec, band_width[k_sort[k]], pulsesnum, nz, &positions[band_start[k_sort[k]]], &inp_vector[band_start[k_sort[k]]] );
                        decode_signs_fx( pardec, band_width[k_sort[k]], &inp_vector[band_start[k_sort[k]]] );

                    }
                    nzbands--;
                }
                else
                {
                    pulsesnum = 0;
                    for( i = 0; i < band_width[k_sort[k]]; i++ )
                    {
                        inp_vector[band_start[k_sort[k]] + i] = 0;
                    }
                }

                npulses[k_sort[k]] = pulsesnum;

                if( Rk_fx[k_sort[k]] > 0 && surplus_fx < 0 )
                {
                    IF( nzbands <= 1 )
                    {
                        divider = 0;
                    }
                    ELSE
                    {
                        divider = 2;
                    }

                    IF( L_add( L_add( surplus_fx, sepbits), ar_div( leftbits, divider ) ) < 0 )
                    {
                        /* Overflow possible => start to distribute negative surplus */
                        delta_fx = ar_div( surplus_fx + sepbits, nzbands);
                    }
                    else
                    {
                        delta_fx = 0;
                    }
                    surplus_fx = L_sub(surplus_fx, delta_fx);
                }
                else
                {
                    delta_fx = 0;
                }
            }
        }

        if (( L_sub(surplus_fx,524288) > 0 && sub(input_frame,L_FRAME8k) == 0 ) || ( L_sub(surplus_fx,786432) > 0 && sub(input_frame,L_FRAME16k) == 0 ))
        {
            bit_surplus_fx[0] = Mult_32_16(surplus_fx,24576);/* Q16 */
            bit_surplus_fx[1] = Mult_32_16(surplus_fx,8192);/* Q16 */
        }
        else
        {
            bit_surplus_fx[0] = surplus_fx;
            bit_surplus_fx[1] = 0;
        }

        for( k = 0; k < BANDS; k++ )
        {
            for( j = 0; j < 2; j++ )
            {
                if( k == k_num[j] )
                {
                    Rk_fx[k_sort[k]] = L_add(Rk_fx[k_sort[k]],bit_surplus_fx[j]);
                    if( Rk_fx[k_sort[k]] > 0 && USQ_TCQ[k_sort[k]] == 0 )
                    {
                        /* get number of pulses */
                        pulsesnum = GetScale_fx( band_width[k_sort[k]], Rk_fx[k_sort[k]], &surplus_fx );

                        decode_position_ari_fx( pardec, band_width[k_sort[k]], pulsesnum, &nz, &positions[band_start[k_sort[k]]] );
                        /* decode tcq magniitude and update the surplus bits. */
                        decode_mangitude_tcq_fx( pardec, band_width[k_sort[k]], pulsesnum, nz, &positions[band_start[k_sort[k]]], &inp_vector[band_start[k_sort[k]]], &surplus_fx );
                        decode_signs_fx( pardec, band_width[k_sort[k]], &inp_vector[band_start[k_sort[k]]] );
                    }
                    else if( Rk_fx[k_sort[k]] > 0 && USQ_TCQ[k_sort[k]] == 1 )
                    {
                        pulsesnum = GetScale_fx( band_width[k_sort[k]], Rk_fx[k_sort[k]], &surplus_fx );

                        decode_position_ari_fx( pardec, band_width[k_sort[k]], pulsesnum, &nz, &positions[band_start[k_sort[k]]] );
                        /* decode usq magnitude and don't need to update surplus bits */
                        decode_magnitude_usq_fx( pardec, band_width[k_sort[k]], pulsesnum, nz, &positions[band_start[k_sort[k]]], &inp_vector[band_start[k_sort[k]]] );
                        decode_signs_fx( pardec, band_width[k_sort[k]], &inp_vector[band_start[k_sort[k]]] );
                    }
                    else
                    {
                        pulsesnum = 0;
                        for ( i = 0; i < band_width[k_sort[k]]; i++ )
                        {
                            inp_vector[band_start[k_sort[k]] + i] = 0;
                        }
                    }
                    npulses[k_sort[k]] = pulsesnum;
                }
            }
        }
    }
    else
    {
        for( k = 0; k < BANDS; k++ )
        {
            if( Rk_fx[k_sort[k]] > 0 )
            {
                pulsesnum = GetScale_fx(band_width[k_sort[k]], Rk_fx[k_sort[k]] + delta_fx, &surplus_fx);

                decode_position_ari_fx( pardec, band_width[k_sort[k]], pulsesnum, &nz, &positions[band_start[k_sort[k]]] );

                /*decode usq magnitude and don't need to update surplus bits*/
                decode_magnitude_usq_fx( pardec, band_width[k_sort[k]], pulsesnum, nz, &positions[band_start[k_sort[k]]], &inp_vector[band_start[k_sort[k]]] );
                decode_signs_fx( pardec, band_width[k_sort[k]], &inp_vector[band_start[k_sort[k]]] );

                nzbands = sub(nzbands, 1);
            }
            else
            {
                pulsesnum = 0;
                for( i = 0; i < band_width[k_sort[k]]; i++ )
                {
                    inp_vector[band_start[k_sort[k]] + i] = 0;
                }
            }

            npulses[k_sort[k]] = pulsesnum;

            /* surplus distribution */
            if ( surplus_fx > 0 && nzbands > 0 )
            {
                delta_fx = ar_div(surplus_fx, nzbands);
                surplus_fx = L_sub(surplus_fx, delta_fx);
            }
        }
    }
    /* Load TCQ path from bitstream */
    LoadTCQdata( pardec, dpath, lsbtcq_bits);

    TCQLSBdec( dpath, mbuffer, 2*lsbtcq_bits );

    ar_decoder_done( pardec );

    /* Restore TCQ */
    if( !flag_wbnb )
    {
        for( k = 0; k < BANDS; k++)
        {
            if( Rk_fx[k_sort[k]] > 0 )
            {
                RestoreTCQdec( &inp_vector[ band_start[ k_sort[ k]]], band_width[k_sort[k]], &bcount, mbuffer );
            }
        }
    }
    else
    {
        for( k = 0; k < BANDS; k++)
        {
            if( Rk_fx[k_sort[k]] > 0 && k != k_num[0] && k != k_num[1] )
            {
                RestoreTCQdec( &inp_vector[ band_start[ k_sort[ k]]], band_width[k_sort[k]], &bcount, mbuffer );
            }
        }
        for( k = 0; k < BANDS; k++)
        {
            if( Rk_fx[k_sort[k]] > 0 && (k == k_num[0] || k == k_num[1]) )
            {
                RestoreTCQdec( &inp_vector[ band_start[ k_sort[ k]]], band_width[k_sort[k]], &bcount, mbuffer );
            }
        }
    }



    return;
}
