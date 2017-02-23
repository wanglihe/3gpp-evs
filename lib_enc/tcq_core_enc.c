/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "options.h"
#include "rom_com.h"
#include "prot.h"
#include "cnst.h"
#include "basop_util.h"
#include "basop_proto_func.h"


/*---------------------------------------------------------------------*
 * tcq_core_LR_enc()
 *
 * Main Generic Audio Encoder Routine for LR-MDCT
 *---------------------------------------------------------------------*/

void tcq_core_LR_enc(
    Encoder_State *st,
    int   inp_vector[],
    const float coefs_norm[],
    float coefs_quant[],
    const short bit_budget,     /* number of bits */
    const short BANDS,
    const short *sfm_start,
    const short *sfm_end,
    const short *sfmsize,
    Word32 *Rk_fx,
    int   *npulses,
    short *k_sort,
    const short *p2a_flags,
    const short p2a_bands,
    const short *last_bitalloc,
    const short input_frame,
    const short adjustFlag,
    const short is_transient
)
{
    short i, j, k, size, nb_bytes;
    int nzp;

    float gain;
    float step_scale[NB_SFM];

    short pos_index[NB_SFM];
    float Rk_sort[NB_SFM];
    int USQ_TCQ[NB_SFM];              /* TCQ is selected by default*/
    float coefs_norm_dec[L_FRAME32k]; /* New output buffer (TCQ+USQ)*/

    float pulses, crosscorr, selfcorr;
    int savedstates[TCQ_MAX_BAND_SIZE];
    ARCODEC arenc, *parenc;
    BITSTREAM bs, *pbs;
    short k_num[2];

    int flag_wbnb = 0;
    int lsbtcq_bits = TCQ_AMP;
    int tcq_arbits  = 2;
    int nzbands = 0;
    short bcount = 0;
    float abuffer[560];
    float mbuffer[560];
    float sbuffer[560];
    short dpath[280];
    /*Word32 Rk_fx[NB_SFM];*/      /* Q16 */
    Word32 Rk_sort_fx[NB_SFM]; /* Q16 */
    Word32 bsub_fx = 0;
    Word32 est_frame_bits_fx;

    Word16 nzb = 0;
    Word32 delta_fx;
    Word32 surplus_fx;
    Word32 bit_surplus_fx[2];

    Word32 leftbits = 0;
    Word32 sepbits = 0;
    Word32 divider = 0;

    set_s(dpath, 0, 280);
    set_f(abuffer, 0.f, 560);
    set_f(mbuffer, 0.f, 560);

    set_f(sbuffer, FLT_MAX, 560);

    /* initialization */
    set_f( Rk_sort, 0.f, NB_SFM );
    set_i( USQ_TCQ, 0, NB_SFM );
    set_f( coefs_norm_dec, 0.f, L_FRAME32k );
    InitLSBTCQ(&bcount);

    if( input_frame <= L_FRAME16k  && adjustFlag == 0 && is_transient == 0 )
    {
        flag_wbnb = 1;
        lsbtcq_bits = 0;
        tcq_arbits  = 0;
    }

    parenc = &arenc;
    pbs = &bs;

    pbs->curPos = 7;
    pbs->numbits = 0;
    pbs->numByte = 0;
    memset(pbs->buf, 0, MAX_SIZEBUF_PBITSTREAM);
    ar_encoder_start(parenc, pbs, bit_budget);

    /* TCQ Index initialize */
    memset( pos_index, 0, sizeof(short)*NB_SFM );

    /* Bits distribution analysis */
    for( i = 0; i < BANDS; i++ )
    {
        if ( L_sub(ar_div(Rk_fx[i], sfmsize[i]), 49152) >= 0 )
        {
            /* USQ used for high importance bands*/
            USQ_TCQ[i] = 1;
        }
        else
        {
            /* TCQ used for usual bands */
            USQ_TCQ[i] = 0;
        }
        if( Rk_fx[i] > 0 )
        {
            nzbands++;
        }
    }

    for( j = 0; j < BANDS; j++ )
    {
        if( Rk_fx[j] > 0 )
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
    /* Quantize spectral band shapes using TCQ */
    /* Select ISC */
    set_f( coefs_quant, 0.0, sfm_end[BANDS-1]+1 );

    mvr2r( coefs_norm, coefs_quant,sfm_end[BANDS-1]+1 );

    delta_fx = 0;
    est_frame_bits_fx = 0;
    if( input_frame <= L_FRAME16k  && adjustFlag == 0 && is_transient == 0 )
    {
        surplus_fx = -131072;
        bit_allocation_second_fx( Rk_fx, Rk_sort_fx, BANDS, sfmsize, k_sort, k_num, p2a_flags, p2a_bands, last_bitalloc, input_frame );

        nzbands = 0;
        for ( j = 0; j < BANDS; j++ )
        {
            if( sub(j, k_num[0]) != 0 && sub(j, k_num[1]) != 0)
            {
                leftbits = L_add( leftbits, Rk_fx[k_sort[j]]);
                if( Rk_fx[k_sort[j]] > 0 )
                {
                    nzbands = add(nzbands, 1);
                }
            }
            else
            {
                sepbits = L_add( sepbits, Rk_fx[k_sort[j]]);
            }
        }

        /* Separate the position information from the input signal(coefs_norm) */
        /* Gather the NZ coefficients*/
        for( k = 0; k < BANDS; k++) /* Loop through non-zero blocks  */
        {
            if( k != k_num[0] && k != k_num[1])
            {
                if( Rk_fx[k_sort[k]] > 0 && USQ_TCQ[k_sort[k]] == 0 ) /* Then have non-zero block AND WILL BE ENCODED BY TCQ */
                {
                    /* Encode Position Info, NZ Info, Signs */
                    size = sfmsize[k_sort[k]];
                    /* Determine scale step, ISC and TCQ quantizer */
                    GetISCScale( &coefs_quant[sfm_start[k_sort[k]]], size, L_add( Rk_fx[k_sort[k]], delta_fx),
                                 &coefs_norm_dec[ sfm_start[ k_sort[k]]], &step_scale[k_sort[k]], &surplus_fx, &pulses, savedstates, 0, &nzp, 0, 0, 0, 0 );
                    leftbits = L_sub( leftbits, L_add( Rk_fx[k_sort[k]], delta_fx) );

                    npulses[ k_sort[k]] = (int)pulses;

                    encode_position_ari_fx( parenc, &coefs_norm_dec[ sfm_start[ k_sort[k]]], size, &est_frame_bits_fx );
                    encode_magnitude_tcq_fx( parenc, &coefs_norm_dec[ sfm_start[ k_sort[k]]], size, npulses[k_sort[k]], nzp, savedstates, &est_frame_bits_fx );
                    encode_signs_fx( parenc, &coefs_norm_dec[ sfm_start[ k_sort[k]]], size, nzp, &est_frame_bits_fx );
                    nzbands--;

                }
                /* Have USQ coded band */
                else if( Rk_fx[k_sort[k]] > 0 && sub(USQ_TCQ[k_sort[k]], 1) == 0 )
                {
                    size = sfmsize[k_sort[k]];
                    GetISCScale( &coefs_quant[ sfm_start[ k_sort[k]]], size, L_add( Rk_fx[k_sort[k]], delta_fx),
                                 &coefs_norm_dec[ sfm_start[ k_sort[k]]], &step_scale[k_sort[k]], &surplus_fx, &pulses, savedstates, 1, &nzp, 0, 0, 0, 0 );
                    leftbits = L_sub( leftbits, L_add( Rk_fx[k_sort[k]], delta_fx) );

                    npulses[ k_sort[ k]] = (int)pulses;

                    encode_position_ari_fx( parenc, &coefs_norm_dec[sfm_start[k_sort[k]]], size, &est_frame_bits_fx );
                    encode_magnitude_usq_fx( parenc, &coefs_norm_dec[sfm_start[k_sort[k]]], size, npulses[k_sort[k]], nzp, &est_frame_bits_fx );
                    encode_signs_fx( parenc, &coefs_norm_dec[ sfm_start[ k_sort[ k]]], size, nzp, &est_frame_bits_fx );

                    nzbands--;
                }
                else /* Then have  zero block  */
                {
                    npulses[ k_sort[ k]] = 0;
                    size = sfmsize[k_sort[k]];
                }

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

                    IF( L_add( L_add( surplus_fx, sepbits), ar_div( leftbits, divider) ) < 0 )
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
            bit_surplus_fx[0] = Mult_32_16(surplus_fx,24576); /* Q16 */
            bit_surplus_fx[1] = Mult_32_16(surplus_fx,8192);  /* Q16 */
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

                    if( Rk_fx[k_sort[k]] > 0 && USQ_TCQ[k_sort[k]] == 0 ) /* Then have non-zero block AND WILL BE ENCODED BY TCQ */
                    {
                        /* Encode Position Info, NZ Info, Signs */
                        size = sfmsize[k_sort[k]];

                        /* Determine scale step, ISC and TCQ quantizer */
                        GetISCScale( &coefs_quant[sfm_start[k_sort[k]]], size, Rk_fx[k_sort[k]], &coefs_norm_dec[ sfm_start[ k_sort[k]]],
                                     &step_scale[k_sort[k]], &surplus_fx, &pulses, savedstates, 0, &nzp, 0, 0, 0, 0 );

                        npulses[ k_sort[k]] = pulses;

                        encode_position_ari_fx( parenc, &coefs_norm_dec[ sfm_start[ k_sort[k]]], size, &est_frame_bits_fx );
                        encode_magnitude_tcq_fx( parenc, &coefs_norm_dec[ sfm_start[ k_sort[k]]], size, npulses[k_sort[k]], nzp, savedstates, &est_frame_bits_fx );
                        encode_signs_fx( parenc, &coefs_norm_dec[ sfm_start[ k_sort[k]]], size, nzp, &est_frame_bits_fx );
                    }
                    /* Have USQ coded band */
                    else if( Rk_fx[k_sort[k]] > 0 && sub(USQ_TCQ[k_sort[k]], 1) == 0 )
                    {
                        size = sfmsize[k_sort[k]];

                        GetISCScale( &coefs_quant[ sfm_start[ k_sort[k]]], size, Rk_fx[k_sort[k]], &coefs_norm_dec[ sfm_start[ k_sort[k]]],
                                     &step_scale[k_sort[k]], &surplus_fx, &pulses, savedstates, 1, &nzp, 0, 0, 0, 0 );

                        npulses[ k_sort[k]] = pulses;

                        encode_position_ari_fx( parenc, &coefs_norm_dec[sfm_start[k_sort[k]]], size, &est_frame_bits_fx );
                        encode_magnitude_usq_fx( parenc, &coefs_norm_dec[sfm_start[k_sort[k]]], size, npulses[k_sort[k]], nzp, &est_frame_bits_fx );
                        encode_signs_fx( parenc, &coefs_norm_dec[ sfm_start[ k_sort[k]]], size, nzp, &est_frame_bits_fx );
                    }
                    else /* Then have  zero block */
                    {
                        npulses[ k_sort[k]] = 0;
                        size = sfmsize[k_sort[k]];
                    }
                }
            }
        }
    }
    else
    {
        surplus_fx = 0;

        /* Separate the position information from the input signal(coefs_norm) */
        /* Gather the NZ coefficients*/
        for( k = 0; k < BANDS; k++) /* Loop through non-zero blocks  */
        {
            if( Rk_fx[k_sort[k]] > 0 )
            {
                size = sfmsize[k_sort[k]];
                GetISCScale( &coefs_quant[ sfm_start[ k_sort[k]]], size, L_add( Rk_fx[k_sort[k]], delta_fx), &coefs_norm_dec[ sfm_start[ k_sort[k]]], &step_scale[k_sort[k]], &surplus_fx, &pulses, savedstates, 1, &nzp, &bcount, abuffer, mbuffer, sbuffer);

                npulses[ k_sort[k]] = pulses;
                encode_position_ari_fx(parenc, &coefs_norm_dec[sfm_start[k_sort[k]]], size, &est_frame_bits_fx);
                encode_magnitude_usq_fx(parenc, &coefs_norm_dec[sfm_start[k_sort[k]]], size, npulses[k_sort[k]], nzp, &est_frame_bits_fx);
                encode_signs_fx(parenc, &coefs_norm_dec[ sfm_start[ k_sort[k]]], size, nzp, &est_frame_bits_fx);

                /* nzbands--;  */
                nzbands = sub(nzbands, 1);
            }
            else /* Then have zero block  */
            {
                npulses[ k_sort[k]] = 0;
                size = sfmsize[k_sort[k]];
            }

            /* Surplus distribution */
            if( surplus_fx > 0 && nzbands > 0 )
            {
                delta_fx = ar_div(surplus_fx, nzbands);
                surplus_fx = L_sub(surplus_fx, delta_fx);
            }
        }
    }

    TCQLSB( bcount, abuffer, mbuffer, sbuffer, dpath );

    /* Save TCQ path to bitstream */
    SaveTCQdata( parenc, dpath, lsbtcq_bits );

    /* Add tcq sequence to decoding buffer */
    InitLSBTCQ( &bcount );

    ar_encoder_done( parenc );

    /* Loop through non-zero blocks   */
    if( !flag_wbnb )
    {
        for( k = 0; k < BANDS; k++ )
        {
            if( Rk_fx[k_sort[k]] > 0 )
            {
                size = sfmsize[k_sort[k]];
                RestoreTCQ( &coefs_norm_dec[ sfm_start[ k_sort[k]]], size, &bcount, mbuffer );
            }
        }
    }

    nb_bytes = bit_budget >> 3;
    j = bit_budget - (nb_bytes << 3);
    for( i = 0; i < nb_bytes; i++ )
    {
        push_indice(st, IND_HQ2_SUBBAND_TCQ, pbs->buf[i], 8);
    }
    if( j > 0 )
    {
        push_indice(st, IND_HQ2_SUBBAND_TCQ, (pbs->buf[nb_bytes] >> (8 - j)), j);
    }

    /* Clear decoding buffer */
    set_f( coefs_quant, 0.0, sfm_end[BANDS-1]+1 );
    /* New analysis of decoded frame */
    for( i = 0; i < BANDS; i++ )
    {
        if( Rk_fx[ k_sort[i]] > 0 )
        {
            gain = 0.0f;

            crosscorr = 0.0f;
            selfcorr = EPSILON;
            for( j = 0; j < sfmsize[k_sort[i]]; j++ )
            {
                crosscorr += (coefs_norm[sfm_start[k_sort[i]]+j] * coefs_norm_dec[sfm_start[k_sort[i]]+j]);
                selfcorr += (coefs_norm_dec[sfm_start[k_sort[i]]+j] * coefs_norm_dec[sfm_start[k_sort[i]]+j]);
            }

            gain = crosscorr / selfcorr;

            if( gain == 0 )
            {
                gain = 1e-10f;
            }

            /* Use optimal gain */
            for( j = 0; j < sfmsize[k_sort[i]]; j++ )
            {
                inp_vector[sfm_start[k_sort[i]]+j] = round_f( ((1.0f/QTCQ)*coefs_norm_dec[sfm_start[k_sort[i]]+j]) );
                coefs_quant[sfm_start[k_sort[i]]+j] = gain*coefs_norm_dec[sfm_start[k_sort[i]]+j];
            }
        }
    }


    return;

}
