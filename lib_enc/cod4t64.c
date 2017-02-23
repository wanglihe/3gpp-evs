/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_enc.h"
#include "prot.h"
#include "rom_com.h"


/*---------------------------------------------------------------------*
 * Local functions
 *---------------------------------------------------------------------*/

static short quant_1p_N1( const short pos, const short N );
static short quant_2p_2N1( const short pos1, const short pos2, const short N );
static short quant_3p_3N1( const short pos1, const short pos2, const short pos3, const short N );
static long quant_4p_4N( const short pos[], const short N );
static long quant_5p_5N( const short pos[], const short N );
static long quant_6p_6N_2( const short pos[], const short N );
static int pre_process( const float v[], short pos_vector[], int pos_vector_num[], int *pulse_pos_num);
static int fcb_encode_position( short pos_vector[], int n, int pos_num, int flag );
static int  fcb_encode_class( int *buffer, int pulse_num, int pos_num );
static int fcb_encode_PI( const float v[], int pulse_num );

/*---------------------------------------------------------------------*
 * ACELP_4t64()
 *
 * 20, 36, 44, 52, 64, 72, 88 bits algebraic codebook.
 * 4 tracks x 16 positions per track = 64 samples.
 *
 * 20 bits --> 4 pulses in a frame of 64 samples.
 * 36 bits --> 8 pulses in a frame of 64 samples.
 * 44 bits 13 + 9 + 13 + 9 --> 10 pulses in a frame of 64 samples.
 * 52 bits 13 + 13 + 13 + 13 --> 12 pulses in a frame of 64 samples.
 * 64 bits 2 + 2 + 2 + 2 + 14 + 14 + 14 + 14 -->
 *                               16 pulses in a frame of 64 samples.
 * 72 bits 10 + 2 + 10 + 2 + 10 + 14 + 10 + 14 -->
 *                               18 pulses in a frame of 64 samples.
 * 88 bits 11 + 11 + 11 + 11 + 11 + 11 + 11 + 11 -->
 *                               24 pulses in a frame of 64 samples.
 * All pulses can have two (2) possible amplitudes: +1 or -1.
 * Each pulse can have sixteen (16) possible positions.
 *---------------------------------------------------------------------*/

short acelp_4t64(
    Encoder_State *st,        /* i/o: encoder state structure                       */
    float dn[],         /* i  : corr. between target and h[].                 */
    const float cn[],         /* i  : residual after long term prediction           */
    const float H[],          /* i  : impulse response of weighted synthesis filter */
    float R[],          /* i  : autocorrelation values                        */
    const short acelpautoc,   /* i  : autocorrealtion flag                          */
    float code[],       /* o  : algebraic (fixed) codebook excitation         */
    float y[],          /* o  : filtered fixed codebook excitation            */
    short nbbits,       /* i  : number of bits per codebook                   */
    const short cmpl_flag,    /* i  : coomplexity reduction flag                    */
    const short Opt_AMR_WB      /* i  : flag indicating AMR-WB IO mode              */
)
{
    short i, k, index, track;
    long L_index;

    short ind[NPMAXPT*NB_TRACK_FCB_4T+32];
    short saved_bits = 0;
    PulseConfig config;
    int indexing_indices[6], wordcnt, bitcnt;


    /*-----------------------------------------------------------------*
     * Configuration
     *-----------------------------------------------------------------*/


    switch (nbbits)
    {
    case 20:          /* EVS/AMR-WB pulse indexing: 20 bits, 4 pulses, 4 tracks  */
        config.nbiter = 4;    /* 4x12x16=768 loop                                        */
        config.alp = 2.0f;
        config.nb_pulse = 4;
        config.fixedpulses = 0;
        config.nbpos[0] = 4;
        config.nbpos[1] = 8;
        break;

    case 28:          /* EVS pulse indexing: 28 bits, 6 pulses, 4 tracks     */
        config.nbiter = 4;    /* 4x20x16=1280 loops                                  */
        config.alp = 1.0f;    /* coeff for sign setting                              */
        config.nb_pulse = 6;
        config.fixedpulses = 0;
        config.nbpos[0] = 6;
        config.nbpos[1] = 6;
        config.nbpos[2] = 8;
        break;

    case 36:          /* EVS/AMR-WB pulse indexing: 36 bits, 8 pulses, 4 tracks  */
        config.nbiter = 4;    /* 4x20x16=1280 loops                                      */
        config.alp = 1.0f;    /* coeff for sign setting                                  */
        config.nb_pulse = 8;
        config.fixedpulses = 2;
        config.nbpos[0] = 4;
        config.nbpos[1] = 8;
        config.nbpos[2] = 8;
        break;

    case 43:          /* EVS pulse indexing:    43 bits, 10 pulses, 4 tracks */
    case 44:          /* AMR-WB pulse indexing: 44 bits, 10 pulses, 4 tracks */
        config.nbiter = 4;    /* 4x26x16=1664 loops                                  */
        config.alp = 1.0f;
        config.nb_pulse = 10;
        config.fixedpulses = 2;
        config.nbpos[0] = 4;
        config.nbpos[1] = 6;
        config.nbpos[2] = 8;
        config.nbpos[3] = 8;
        break;

    case 50:          /* EVS pulse indexing:    50 bits, 12 pulses, 4 tracks */
    case 52:          /* AMR-WB pulse indexing: 52 bits, 12 pulses, 4 tracks */
        config.nbiter = 4;    /* 4x26x16=1664 loops                                  */
        config.alp = 1.0f;
        config.nb_pulse = 12;
        config.fixedpulses = 4;
        config.nbpos[0] = 4;
        config.nbpos[1] = 6;
        config.nbpos[2] = 8;
        config.nbpos[3] = 8;
        break;

    case 62:          /* EVS pulse indexing:    62 bits, 16 pulses, 4 tracks */
    case 64:          /* AMR-WB pulse indexing: 64 bits, 16 pulses, 4 tracks */
        config.nbiter = 3;    /* 3x36x16=1728 loops                                  */
        config.alp = 0.8F;
        config.nb_pulse = 16;
        config.fixedpulses = 4;
        config.nbpos[0] = 4;
        config.nbpos[1] = 4;
        config.nbpos[2] = 6;
        config.nbpos[3] = 6;
        config.nbpos[4] = 8;
        config.nbpos[5] = 8;
        break;

    case 72:          /* AMR-WB pulse indexing: 72 bits, 18 pulses, 4 tracks */
        config.nbiter = 3;    /* 3x35x16=1680 loops                                  */
        config.alp = 0.75F;
        config.nb_pulse = 18;
        config.fixedpulses = 4;
        config.nbpos[0] = 2;
        config.nbpos[1] = 3;
        config.nbpos[2] = 4;
        config.nbpos[3] = 5;
        config.nbpos[4] = 6;
        config.nbpos[5] = 7;
        config.nbpos[6] = 8;
        break;

    case 88:          /* AMR-WB pulse indexing: 88 bits, 24 pulses, 4 tracks */
        config.nbiter = 2;    /* 2x53x16=1696 loop                                   */
        config.alp = 0.5f;
        config.nb_pulse = 24;
        config.fixedpulses = 4;
        config.nbpos[0] = 2;
        config.nbpos[1] = 2;
        config.nbpos[2] = 3;
        config.nbpos[3] = 4;
        config.nbpos[4] = 5;
        config.nbpos[5] = 6;
        config.nbpos[6] = 7;
        config.nbpos[7] = 8;
        config.nbpos[8] = 8;
        config.nbpos[9] = 8;
        break;

    case 87:          /* EVS pulse indexing:   87 bits, 26 pulses, 4 tracks  */
        config.nbiter = 1;
        config.alp = 0.5F;
        config.nb_pulse = 26;
        config.fixedpulses = 4;
        config.nbpos[0] = 4;
        config.nbpos[1] = 6;
        config.nbpos[2] = 6;
        config.nbpos[3] = 8;
        config.nbpos[4] = 8;
        config.nbpos[5] = 8;
        config.nbpos[6] = 8;
        config.nbpos[7] = 8;
        config.nbpos[8] = 8;
        config.nbpos[9] = 8;
        config.nbpos[10] = 8;
        break;
    }

    /* reduce the number of iterations as a compromise between the performance and complexity */
    if( cmpl_flag > 0 )
    {
        config.nbiter = cmpl_flag;
    }

    config.codetrackpos = TRACKPOS_FIXED_FIRST;
    config.bits = nbbits;

    /*-----------------------------------------------------------------*
     * Search
     *-----------------------------------------------------------------*/

    if( acelpautoc )
    {
        E_ACELP_4tsearchx( dn, cn, R, code, &config, ind );

        /* Generate weighted code */
        set_f( y, 0.0f, L_SUBFR );
        for( i=0; i<L_SUBFR; i++ )
        {
            /* Code is sparse, so check which samples are non-zero */
            if( code[i] != 0 )
            {
                for( k=0; k<L_SUBFR-i; k++ )
                {
                    y[i+k] += code[i] * H[k];
                }
            }
        }
    }
    else
    {
        E_ACELP_4tsearch( dn, cn, H, code, &config, ind, y );
    }

    /*-----------------------------------------------------------------*
     * Indexing
     *-----------------------------------------------------------------*/

    if( !Opt_AMR_WB )
    {
        /* EVS pulse indexing */

        saved_bits = E_ACELP_indexing( code, config, NB_TRACK_FCB_4T, indexing_indices );

        saved_bits = 0;

        wordcnt = nbbits >> 4;
        bitcnt = nbbits & 15;
        for ( i = 0; i < wordcnt; i++ )
        {
            push_indice( st, IND_ALG_CDBK_4T64, indexing_indices[i], 16 );
        }
        if ( bitcnt )
        {
            push_indice( st, IND_ALG_CDBK_4T64, indexing_indices[i], bitcnt );
        }

    }
    else
    {
        /* AMR-WB pulse indexing */

        if (nbbits == 20)
        {
            for (track = 0; track < NB_TRACK_FCB_4T; track++)
            {
                k = track * NPMAXPT;
                index = quant_1p_N1(ind[k], 4);
                push_indice( st, IND_ALG_CDBK_4T64, index, 5 );
            }
        }
        else if (nbbits == 36)
        {
            for (track = 0; track < NB_TRACK_FCB_4T; track++)
            {
                k = track * NPMAXPT;
                index = quant_2p_2N1(ind[k], ind[k+1], 4);
                push_indice( st, IND_ALG_CDBK_4T64, index, 9 );
            }
        }
        else if (nbbits == 44)
        {
            for (track = 0; track < (NB_TRACK_FCB_4T - 2); track++)
            {
                k = track * NPMAXPT;
                index = quant_3p_3N1(ind[k], ind[k+1], ind[k+2], 4);
                push_indice( st, IND_ALG_CDBK_4T64, index, 13 );
            }

            for (track = 2; track < NB_TRACK_FCB_4T; track++)
            {
                k = track * NPMAXPT;
                index = quant_2p_2N1(ind[k], ind[k+1], 4);
                push_indice( st, IND_ALG_CDBK_4T64, index, 9 );
            }
        }
        else if (nbbits == 52)
        {
            for (track = 0; track < NB_TRACK_FCB_4T; track++)
            {
                k = track*NPMAXPT;
                index = quant_3p_3N1(ind[k], ind[k+1], ind[k+2], 4);
                push_indice( st, IND_ALG_CDBK_4T64, index, 13 );
            }
        }
        else if (nbbits == 64)
        {
            for (track = 0; track < NB_TRACK_FCB_4T; track++)
            {
                k = track * NPMAXPT;
                L_index = quant_4p_4N(&ind[k], 4);
                index = ((L_index >> 14) & 3);
                push_indice( st, IND_ALG_CDBK_4T64_1, index, 2 );
            }

            for (track = 0; track < NB_TRACK_FCB_4T; track++)
            {
                k = track * NPMAXPT;
                L_index = quant_4p_4N(&ind[k], 4);
                index = (L_index & 0x3FFF);
                push_indice( st, IND_ALG_CDBK_4T64_2, index, 14 );
            }
        }
        else if (nbbits == 72)
        {
            for (track=0; track< (NB_TRACK_FCB_4T - 2); track++)
            {
                k = track * NPMAXPT;
                L_index = quant_5p_5N(&ind[k], 4);
                index = ((L_index >> 10) & 0x03FF);
                push_indice( st, IND_ALG_CDBK_4T64_1, index, 10 );
            }

            for (track = 2; track < NB_TRACK_FCB_4T; track++)
            {
                k = track * NPMAXPT;
                L_index = quant_4p_4N(&ind[k], 4);
                index = ((L_index >> 14) & 3);
                push_indice( st, IND_ALG_CDBK_4T64_1, index, 2 );
            }

            for (track=0; track< (NB_TRACK_FCB_4T - 2); track++)
            {
                k = track * NPMAXPT;
                L_index = quant_5p_5N(&ind[k], 4);
                index = (L_index & 0x03FF);
                push_indice( st, IND_ALG_CDBK_4T64_2, index, 10 );
            }

            for (track = 2; track < NB_TRACK_FCB_4T; track++)
            {
                k = track * NPMAXPT;
                L_index = quant_4p_4N(&ind[k], 4);
                index = (L_index & 0x3FFF);
                push_indice( st, IND_ALG_CDBK_4T64_2, index, 14 );
            }
        }
        else if (nbbits == 88)
        {
            for (track = 0; track < NB_TRACK_FCB_4T; track++)
            {
                k = track * NPMAXPT;
                L_index = quant_6p_6N_2(&ind[k], 4);
                index = ((L_index >> 11) & 0x07FF);
                push_indice( st, IND_ALG_CDBK_4T64_1, index, 11 );
            }

            for (track = 0; track < NB_TRACK_FCB_4T; track++)
            {
                k = track * NPMAXPT;
                L_index = quant_6p_6N_2(&ind[k], 4);
                index = (L_index & 0x07FF);
                push_indice( st, IND_ALG_CDBK_4T64_2, index, 11 );
            }
        }
    }

    return saved_bits;
}


/*---------------------------------------------------------------------*
 * Quantization of 1 pulse with N+1 bits:                              *
 *---------------------------------------------------------------------*/
static short quant_1p_N1(  /* o:   return N+1 bits             */
    const short pos,       /* i:   position of the pulse       */
    const short N          /* i:   number of bits for position */
)
{
    short mask, index;

    mask = ((1<<N)-1);

    index = (pos & mask);

    if ((pos & NB_POS_FCB_4T) != 0)
    {
        index += 1 << N;
    }

    return index;
}


/*---------------------------------------------------------------------*
 * Quantization of 2 pulses with 2*N+1 bits:                           *
 *---------------------------------------------------------------------*/
static short quant_2p_2N1( /* o:   return (2*N)+1 bits         */
    const short pos1,      /* i:   position of the pulse 1     */
    const short pos2,      /* i:   position of the pulse 2     */
    const short N          /* i:   number of bits for position */
)
{
    short mask, index;

    mask = ((1<<N)-1);

    /*-----------------------------------------------------------------*
     * sign of 1st pulse == sign of 2nd pulse
     *-----------------------------------------------------------------*/

    if (((pos2 ^ pos1) & NB_POS_FCB_4T) == 0)
    {
        if ((pos1 - pos2) <= 0)
        {
            index = ((pos1 & mask) << N) + (pos2 & mask);
        }
        else
        {
            index = ((pos2 & mask) << N) + (pos1 & mask);
        }
        if ((pos1 & NB_POS_FCB_4T) != 0)
        {
            index += 1 << (2*N);
        }
    }
    else
    {

        /*-----------------------------------------------------------------*
         * sign of 1st pulse != sign of 2nd pulse
         *-----------------------------------------------------------------*/
        if (((pos1 & mask) - (pos2 & mask)) <= 0)
        {
            index = ((pos2 & mask) << N) + (pos1 & mask);
            if ((pos2 & NB_POS_FCB_4T) != 0)
            {
                index += 1 << (2*N);
            }
        }
        else
        {
            index = ((pos1 & mask) << N) + (pos2 & mask);
            if ((pos1 & NB_POS_FCB_4T) != 0)
            {
                index += 1 << (2*N);
            }
        }
    }

    return index;
}

/*---------------------------------------------------------------------*
 * Quantization of 3 pulses with 3*N+1 bits:                           *
 *---------------------------------------------------------------------*/
static short quant_3p_3N1( /* o:   return (3*N)+1 bits         */
    const short pos1,      /* i:   position of the pulse 1     */
    const short pos2,      /* i:   position of the pulse 2     */
    const short pos3,      /* i:   position of the pulse 3     */
    const short N          /* i:   number of bits for position */
)
{
    short index, nb_pos;

    nb_pos = (1 << (N-1));

    /* Quantization of 3 pulses with 3*N+1 bits */
    if (((pos1 ^ pos2) & nb_pos) == 0)
    {
        index = quant_2p_2N1(pos1, pos2, (N - 1));
        index += (pos1 & nb_pos) << N;
        index += quant_1p_N1(pos3, N) << (2 * N);
    }
    else if (((pos1 ^ pos3) & nb_pos) == 0)
    {
        index = quant_2p_2N1(pos1, pos3, (N - 1));
        index += (pos1 & nb_pos) << N;
        index += quant_1p_N1(pos2, N) << (2 * N);
    }
    else
    {
        index = quant_2p_2N1(pos2, pos3, (N - 1));
        index += (pos2 & nb_pos) << N;
        index += quant_1p_N1(pos1, N) << (2 * N);
    }

    return index;
}

/*---------------------------------------------------------------------*
 * Quantization of 4 pulses with 4*N+1 bits:                           *
 *---------------------------------------------------------------------*/
static long quant_4p_4N1( /* o:   return (4*N)+1 bits         */
    const short pos1,      /* i:   position of the pulse 1     */
    const short pos2,      /* i:   position of the pulse 2     */
    const short pos3,      /* i:   position of the pulse 3     */
    const short pos4,      /* i:   position of the pulse 4     */
    const short N          /* i:   number of bits for position */
)
{
    long index, nb_pos;

    nb_pos = (1 << (N-1));

    /* Quantization of 4 pulses with 4*N+1 bits */
    if (((pos1 ^ pos2) & nb_pos) == 0)
    {
        index = quant_2p_2N1(pos1, pos2, (N - 1));
        index += (pos1 & nb_pos) << N;
        index += quant_2p_2N1(pos3, pos4, N) << (2 * N);
    }
    else if (((pos1 ^ pos3) & nb_pos) == 0)
    {
        index = quant_2p_2N1(pos1, pos3, (N - 1));
        index += (pos1 & nb_pos) << N;
        index += quant_2p_2N1(pos2, pos4, N) << (2 * N);
    }
    else
    {
        index = quant_2p_2N1(pos2, pos3, (N - 1));
        index += (pos2 & nb_pos) << N;
        index += quant_2p_2N1(pos1, pos4, N) << (2 * N);
    }

    return (index);
}

/*---------------------------------------------------------------------*
 * Quantization of 4 pulses with 4*N bits:                             *
 *---------------------------------------------------------------------*/
static long quant_4p_4N(  /* o:   return 4*N bits             */
    const short pos[],     /* i:   position of the pulse 1..4  */
    const short N          /* i:   number of bits for position */
)
{
    short i, j, k, n_1;
    short posA[4], posB[4];
    long nb_pos, index = 0;

    n_1 = N - 1;
    nb_pos = (1 << n_1);

    i = 0;
    j = 0;
    for (k = 0; k < 4; k++)
    {
        if ((pos[k] & nb_pos) == 0)
        {
            posA[i++] = pos[k];
        }
        else
        {
            posB[j++] = pos[k];
        }
    }

    switch (i)
    {
    case 0:
        index = 1 << ((4 * N) - 3);
        index += quant_4p_4N1(posB[0], posB[1], posB[2], posB[3], n_1);
        break;
    case 1:
        index = quant_1p_N1(posA[0], n_1) << (( 3 * n_1) + 1);
        index += quant_3p_3N1(posB[0], posB[1], posB[2], n_1);
        break;
    case 2:
        index = quant_2p_2N1(posA[0], posA[1], n_1) << (( 2 * n_1) + 1);
        index += quant_2p_2N1(posB[0], posB[1], n_1);
        break;
    case 3:
        index = quant_3p_3N1(posA[0], posA[1], posA[2], n_1) << N;
        index += quant_1p_N1(posB[0], n_1);
        break;
    case 4:
        index = quant_4p_4N1(posA[0], posA[1], posA[2], posA[3], n_1);
        break;
    }
    index += (i & 3) << ((4 * N) - 2);

    return (index);
}

/*---------------------------------------------------------------------*
 * Quantization of 5 pulses with 5*N bits:                             *
 *---------------------------------------------------------------------*/
static long quant_5p_5N(  /* o:   return 5*N bits             */
    const short pos[],     /* i:   position of the pulse 1..5  */
    const short N          /* i:   number of bits for position */
)
{
    short i, j, k, n_1, nb_pos;
    short posA[5], posB[5];
    long index = 0;

    n_1 = N-1;
    nb_pos = (1 << n_1);

    i = 0;
    j = 0;
    for (k = 0; k < 5; k++)
    {
        if ((pos[k] & nb_pos) == 0)
        {
            posA[i++] = pos[k];
        }
        else
        {
            posB[j++] = pos[k];
        }
    }
    switch (i)
    {
    case 0:
        index = 1 << ((5 * N) - 1);
        index += quant_3p_3N1(posB[0], posB[1], posB[2], n_1) << ((2 * N) + 1);
        index += quant_2p_2N1(posB[3], posB[4], N);
        break;
    case 1:
        index = 1 << ((5 * N) - 1);
        index += quant_3p_3N1(posB[0], posB[1], posB[2], n_1) << ((2 * N) + 1);
        index += quant_2p_2N1(posB[3], posA[0], N);
        break;
    case 2:
        index = 1 << ((5 * N) - 1);
        index += quant_3p_3N1(posB[0], posB[1], posB[2], n_1) << ((2 * N) + 1);
        index += quant_2p_2N1(posA[0], posA[1], N);
        break;
    case 3:
        index = quant_3p_3N1(posA[0], posA[1], posA[2], n_1) << ((2 * N) + 1);
        index += quant_2p_2N1(posB[0], posB[1], N);
        break;
    case 4:
        index = quant_3p_3N1(posA[0], posA[1], posA[2], n_1) << ((2 * N) + 1);
        index += quant_2p_2N1(posA[3], posB[0], N);
        break;
    case 5:
        index = quant_3p_3N1(posA[0], posA[1], posA[2], n_1) << ((2 * N) + 1);
        index += quant_2p_2N1(posA[3], posA[4], N);
        break;
    }

    return (index);
}

/*---------------------------------------------------------------------*
 * Quantization of 6 pulses with 6*N-2 bits:                           *
 *---------------------------------------------------------------------*/
static long quant_6p_6N_2(/* o:   return 6*N-2 bits           */
    const short pos[],     /* i:   position of the pulse 1..6  */
    const short N          /* i:   number of bits for position */
)
{
    short i, j, k, n_1;
    short posA[6], posB[6];
    long nb_pos, index = 0;

    n_1 = N - 1;
    nb_pos = 1 << n_1;

    i = 0;
    j = 0;
    for (k = 0; k < 6; k++)
    {
        if ((pos[k] & nb_pos) == 0)
        {
            posA[i++] = pos[k];
        }
        else
        {
            posB[j++] = pos[k];
        }
    }
    switch (i)
    {
    case 0:
        index = 1 << ((6 * N) - 5);
        index += quant_5p_5N(posB, n_1) << N;
        index += quant_1p_N1(posB[5], n_1);
        break;
    case 1:
        index = 1 << ((6 * N) - 5);
        index += quant_5p_5N(posB, n_1) << N;
        index += quant_1p_N1(posA[0], n_1);
        break;
    case 2:
        index = 1 << ((6 * N) - 5);
        index += quant_4p_4N(posB, n_1) << ((2 * n_1) + 1);
        index += quant_2p_2N1(posA[0], posA[1], n_1);
        break;
    case 3:
        index = quant_3p_3N1(posA[0], posA[1], posA[2], n_1) << ((3 * n_1) + 1);
        index += quant_3p_3N1(posB[0], posB[1], posB[2], n_1);
        break;
    case 4:
        i = 2;
        index = quant_4p_4N(posA, n_1) << ((2 * n_1) + 1);
        index += quant_2p_2N1(posB[0], posB[1], n_1);
        break;
    case 5:
        i = 1;
        index = quant_5p_5N(posA, n_1) << N;
        index += quant_1p_N1(posB[0], n_1);
        break;
    case 6:
        i = 0;
        index = quant_5p_5N(posA, n_1) << N;
        index += quant_1p_N1(posA[5], n_1);
        break;
    }

    index += (i & 3) << ((6 * N) - 4);

    return (index);
}

/*---------------------------------------------------------------------*
*order the pulse position                                             *
*---------------------------------------------------------------------*/

static int pre_process( /* o:   return sign value of pulse on a track              */
    const float v[],        /* i:   the pulse vector                                   */
    short pos_vector[],     /* o:   position of the pulse on a track                   */
    int pos_vector_num[],   /* o:   the pulse number on the position which have pulse  */
    int *pulse_pos_num      /* i:   the number of position which have pulse            */
)
{
    int  j,k;
    int  sign;

    sign = 0;
    j = 0;
    for (k=0; k<64; k+=4)
    {
        if (v[k])
        {
            pos_vector[j] = k>>2;
            pos_vector_num[j] = fabs(v[k]);
            if (v[k]>0)
                sign = sign << 1;
            else
                sign = ( sign << 1 ) + 1;
            j++;
        }
    }
    *pulse_pos_num = j;

    return sign;
}

/*---------------------------------------------------------------------*
 *encode the position                                                  *
 *---------------------------------------------------------------------*/

static int fcb_encode_position( /* o:   return index of the positions which have pulse*/
    short pos_vector[],         /* i:   position of the pulse on a track              */
    int n,
    int pos_num,                /* i:   the number of position which have pulse   */
    int flag
)
{
    int i;
    int mmm1;
    int temp2;
    mmm1 = PI_select_table[n][pos_num] - 1;
    temp2 = pos_num;

    if (flag)        /* no decrease */
    {
        for (i=0; i<pos_num; i++)
        {
            mmm1 -= PI_select_table[n-pos_vector[i]-1][temp2--];
        }
    }
    else
    {
        for (i=0; i<pos_num; i++)
        {
            mmm1 -= PI_select_table[n-pos_vector[i]-1][temp2--];
            n--;
        }
    }

    return mmm1;
}

/*---------------------------------------------------------------------*
 *encode class for 3p 4p 5p 6p/track                                   *
 *---------------------------------------------------------------------*/

static int fcb_encode_cl(/* o:   class index of the pulse on a track       */
    int buffer[],      /* i:   pulses on a track                         */
    int pulse_num,     /* i:   pulses number on a track                  */
    int pos_num        /* i:   number of the position which have pulse   */
)
{
    int  i,k;
    int  temp1,temp2;
    temp1 = pos_num + pulse_num - 1;
    temp2 = pulse_num;
    k = PI_select_table[temp1][pulse_num] - 1;
    temp1 --;
    for (i=0; i<pulse_num; i++)
    {
        k -= PI_select_table[temp1-buffer[i]][temp2--];
        temp1--;
    }

    return k;
}

/*---------------------------------------------------------------------*
 *encode the class and compute class offset                            *
 *---------------------------------------------------------------------*/

static int fcb_encode_class(/* o:   class offset        */
    int sector_6p_num[],   /* i:   position which have pulse on a track             */
    int pulse_num,         /* i:   pulse number on a track                          */
    int pulse_pos_num      /* i:   number of position which have pulse on a track   */
)
{
    int i,j,k;
    int mn9_offet;
    int vector_class[6];
    int *vector_class_ptr;
    mn9_offet = 0;
    if ( pulse_pos_num < pulse_num )
    {
        vector_class_ptr = vector_class;
        for (i=0; i<pulse_pos_num; i++)
        {
            for (j=0; j<(sector_6p_num[i]-1); j++)
            {
                *vector_class_ptr++ = i ;
            }
        }
        k = fcb_encode_cl(vector_class,pulse_num-pulse_pos_num,pulse_pos_num);
        mn9_offet = PI_factor[pulse_pos_num] * k ;
    }

    return mn9_offet;
}


/*---------------------------------------------------------------------*
 *encode fcb pulse index                                               *
 *---------------------------------------------------------------------*/

static int fcb_encode_PI( /* o:   return index of the  pulse on a track */
    const float v[],    /* i:   the pulse vector                      */
    int pulse_num       /* i:   number of the pulse on a track        */
)
{
    short vector_p[7];
    int  pulse_pos_num;
    int  vector_p_num[7];
    int  code_index;
    int  sign;

    /*order the pulse position*/
    sign = pre_process(v, vector_p, vector_p_num, &pulse_pos_num);

    /*encode the position*/
    code_index  = fcb_encode_position(vector_p,16,pulse_pos_num,1);

    /*encode the class and compute class offset*/
    code_index += fcb_encode_class(vector_p_num,pulse_num,pulse_pos_num);

    code_index = PI_offset[pulse_num][pulse_num + 1 - pulse_pos_num]  +  ( code_index << pulse_pos_num ) + sign;

    return code_index;
}



/*--------------------------------------------------------------------------*
* E_ACELP_code43bit
*
* Fixed bit-length arithmetic coding of pulses
* v - (input) pulse vector
* s - (output) encoded state
* n - (output) range of possible states (0...n-1)
* p - (output) number of pulses found
* len - (input) length of pulse vector
* trackstep - (input) step between tracks
*--------------------------------------------------------------------------*/

short E_ACELP_code43bit(const float code[], long unsigned *ps, int *p, unsigned short idxs[])
{
    int j,k,track;
    int ind[32];

    int tmp;
    int joint_index;
    int joint_offset = 3611648;        /*offset for 3 pulses per track*/
    short saved_bits = 0;

    for (track = 0; track< 2; track++)
    {
        k = track * NPMAXPT;
        ps[track] = fcb_encode_PI(code+track,3);
        p[track] = 3;
    }

    for (track = 2; track < NB_TRACK_FCB_4T; track++)
    {
        j = track * NPMAXPT;
        for (k=track; k<64; k+=4)
        {
            if (code[k])
            {
                tmp = k>>2;
                if (code[k]<0)
                {
                    tmp += 16;
                }
                if (fabs(code[k])>1)
                {
                    ind[j] = tmp;
                    ind[j+1] = tmp;
                    break;
                }
                else
                {
                    ind[j] = tmp;
                    j++;
                }
            }
        }
        k = track * NPMAXPT;
        ps[track] = quant_2p_2N1(ind[k], ind[k+1], 4);
        p[track] = 2;
    }
    joint_index = ps[0]*5472 + ps[1];
    if (joint_index >= joint_offset)
    {
        joint_index += joint_offset;
    }
    else
    {
        saved_bits += 1;
    }

    idxs[0] = ( ( ps[2] << 9 ) + ps[3] ) & 0xffff;
    idxs[1] = ( ( joint_index << 2 ) + ( ps[2] >> 7 ) ) & 0xffff;
    idxs[2] = joint_index >> 14;

    return saved_bits;
}

