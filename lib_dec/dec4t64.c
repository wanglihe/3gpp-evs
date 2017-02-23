/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"


/*-------------------------------------------------------------------*
 * Local functions
 *-------------------------------------------------------------------*/

static void add_pulses(const short pos[], const short nb_pulse, const short track, float code[]);
static void dec_1p_N1(const long index, const short N, const short offset, short pos[]);
static void dec_2p_2N1(const long index, const short N, const short offset, short pos[]);
static void dec_3p_3N1(const long index, const short N, const short offset, short pos[]);
static void dec_4p_4N1(const long index, const short N, const short offset, short pos[]);
static void dec_4p_4N(const long index, const short N, const short offset, short pos[]);
static void dec_5p_5N(const long index, const short N, const short offset, short pos[]);
static void dec_6p_6N2(const long index, const short N, const short offset, short pos[]);
static int fcb_decode_PI( int code_index, short sector_6p[], int pulse_num );

/*----------------------------------------------------------------------------------*
 * dec_acelp_4t64()
 *
 * 20, 36       bits algebraic codebook decoder.
 * 4 tracks x 16 positions per track = 64 samples.
 *
 * 20 bits --> 4 pulses in a frame of 64 samples.
 * 36 bits --> 8 pulses in a frame of 64 samples.
 *
 * All pulses can have two (2) possible amplitudes: +1 or -1.
 * Each pulse can have sixteen (16) possible positions.
 *
 * See cod4t64.c for more details of the algebraic code.
 *----------------------------------------------------------------------------------*/

void dec_acelp_4t64(
    Decoder_State *st,            /* i/o: decoder state structure   */
    short nbbits,         /* i  : number of bits per codebook           */
    float code[],         /* o  : algebraic (fixed) codebook excitation */
    const short Opt_AMR_WB                  /* i  : flag indicating AMR-WB IO mode                */
)
{
    short i, k, pos[7];
    long L_index;
    long ind1[NB_TRACK_FCB_4T], ind2[NB_TRACK_FCB_4T];
    PulseConfig config;
    int indexing_indices[6], wordcnt, bitcnt;

    if ( !Opt_AMR_WB )
    {
        switch (nbbits)
        {
        case 20:
            config.nb_pulse = 4;
            break;

        case 28:
            config.nb_pulse = 6;
            break;

        case 36:
            config.nb_pulse = 8;
            break;

        case 43:
            config.nb_pulse = 10;
            break;

        case 50:
            config.nb_pulse = 12;
            break;

        case 62:
            config.nb_pulse = 16;
            break;


        case 87:
            config.nb_pulse = 26;
            break;
        }

        config.bits = nbbits;
        config.codetrackpos = TRACKPOS_FIXED_FIRST;


        wordcnt = nbbits >> 4;
        bitcnt = nbbits & 15;
        for ( i = 0; i < wordcnt; i++ )
        {
            indexing_indices[i] = get_next_indice( st, 16 );
        }
        if ( bitcnt )
        {
            indexing_indices[i] = get_next_indice( st, bitcnt );
        }

        D_ACELP_indexing( code, config, NB_TRACK_FCB_4T, indexing_indices, &st->BER_detect );
    }
    else
    {
        for (i=0; i<L_SUBFR; i++)
        {
            code[i] = 0.0f;
        }

        if (nbbits == 20)
        {
            for (k=0; k<NB_TRACK_FCB_4T; k++)
            {
                L_index = get_next_indice( st, 5 );
                dec_1p_N1(L_index, 4, 0, pos);
                add_pulses(pos, 1, k, code);
            }
        }
        else if (nbbits == 36)
        {
            for (k=0; k<NB_TRACK_FCB_4T; k++)
            {
                L_index = get_next_indice( st, 9 );
                dec_2p_2N1(L_index, 4, 0, pos);
                add_pulses(pos, 2, k, code);
            }
        }
        else if (nbbits == 44)    /* AMR-WB pulse indexing */
        {
            for(k = 0; k < NB_TRACK_FCB_4T - 2; k++)
            {
                L_index = get_next_indice( st, 13 );
                dec_3p_3N1(L_index, 4, 0, pos);
                add_pulses(pos, 3, k, code);
            }

            for(k = 2; k < NB_TRACK_FCB_4T; k++)
            {
                L_index = get_next_indice( st, 9 );
                dec_2p_2N1(L_index, 4, 0, pos);
                add_pulses(pos, 2, k, code);
            }
        }
        else if (nbbits == 52)    /* AMR-WB pulse indexing */
        {
            for(k = 0; k < NB_TRACK_FCB_4T; k++)
            {
                L_index = get_next_indice( st, 13 );
                dec_3p_3N1(L_index, 4, 0, pos );
                add_pulses(pos, 3, k, code );
            }
        }
        else if (nbbits == 64)    /* AMR-WB pulse indexing */
        {
            for(k = 0; k < NB_TRACK_FCB_4T; k++)
            {
                ind1[k] = get_next_indice( st, 2 );
            }

            for (k = 0; k < NB_TRACK_FCB_4T; k++)
            {
                ind2[k] = get_next_indice( st, 14 );
            }

            for (k = 0; k < NB_TRACK_FCB_4T; k++)
            {
                L_index = ((ind1[k]<<14) + ind2[k]);
                dec_4p_4N(L_index, 4, 0, pos);
                add_pulses(pos, 4, k, code);
            }
        }
        else if (nbbits == 72)
        {
            for(k = 0; k < NB_TRACK_FCB_4T - 2; k++)
            {
                ind1[k] = get_next_indice( st, 10 );
            }

            for(k = 2; k < NB_TRACK_FCB_4T; k++)
            {
                ind1[k] = get_next_indice( st, 2 );
            }

            for(k = 0; k < NB_TRACK_FCB_4T - 2; k++)
            {
                ind2[k] = get_next_indice( st, 10 );
            }

            for(k = 2; k < NB_TRACK_FCB_4T; k++)
            {
                ind2[k] = get_next_indice( st, 14 );
            }

            for(k = 0; k < NB_TRACK_FCB_4T - 2; k++)
            {
                L_index = ((ind1[k]<<10) + ind2[k]);
                dec_5p_5N(L_index, 4, 0, pos);
                add_pulses(pos, 5, k, code);
            }

            for(k = 2; k < NB_TRACK_FCB_4T; k++)
            {
                L_index = ((ind1[k]<<14) + ind2[k]);
                dec_4p_4N(L_index, 4, 0, pos);
                add_pulses(pos, 4, k, code);
            }

        }
        else if (nbbits == 88)
        {
            for(k = 0; k < NB_TRACK_FCB_4T; k++)
            {
                ind1[k] = get_next_indice( st, 11 );
            }

            for(k = 0; k < NB_TRACK_FCB_4T; k++)
            {
                ind2[k] = get_next_indice( st, 11 );
            }

            for(k = 0; k < NB_TRACK_FCB_4T; k++)
            {
                L_index = ((ind1[k]<<11) + ind2[k]);
                dec_6p_6N2(L_index, 4, 0, pos);
                add_pulses(pos, 6, k, code);
            }
        }
    }

    return;
}

/*-------------------------------------------------------*
 * add_pulses()
 *
 * Add decoded pulses to the codeword
 *-------------------------------------------------------*/

static void add_pulses(
    const short pos[],     /* i:   pulse position     */
    const short nb_pulse,  /* i:   nb. of pulses      */
    const short track,     /* i:   no. of the tracks  */
    float code[]     /* i/o: decoded codevector */
)
{
    short i, k;

    for (k=0; k<nb_pulse; k++)
    {
        i = ((pos[k] & (NB_POS_FCB_4T-1))*NB_TRACK_FCB_4T) + track;
        if ((pos[k] & NB_POS_FCB_4T) == 0)
        {
            code[i] += 1.0f;
        }
        else
        {
            code[i] -= 1.0f;
        }
    }

    return;
}

/*-------------------------------------------------------*
 * dec_1p_N1()
 *
 * Decode 1 pulse with N+1 bits
 *-------------------------------------------------------*/

void dec_1p_N1(
    const  long index,   /* i:   quantization index    */
    const short N,       /* i:   nb. of bits           */
    const short offset,  /* i:   pulse position offset */
    short pos[]    /* o:   pulse position        */
)
{
    short i, pos1;
    long mask;

    mask = ((1<<N)-1);

    pos1 = (short)(index & mask) + offset;

    i = (short)(index >> N) & 1;

    if (i == 1)
    {
        pos1 += NB_POS_FCB_4T;
    }

    pos[0] = pos1;

    return;
}

/*-------------------------------------------------------*
 * dec_2p_2N1()
 *
 * Decode 2 pulses with 2*N+1 bits:
 *-------------------------------------------------------*/

void dec_2p_2N1(
    const long  index,   /* i:   quantization index    */
    const short N,       /* i:   nb. of bits           */
    const short offset,  /* i:   pulse position offset */
    short pos[]    /* o:   pulse position        */
)
{
    short i, pos1, pos2;
    long mask;

    mask = ((1<<N)-1);

    pos1 = (short)((index >> N) & mask) + offset;
    i = (short)(index >> (2*N)) & 1;
    pos2 = (short)(index & mask) + offset;
    if ((pos2 - pos1) < 0)
    {
        if (i == 1)
        {
            pos1 += NB_POS_FCB_4T;
        }
        else
        {
            pos2 += NB_POS_FCB_4T;
        }
    }
    else
    {
        if (i == 1)
        {
            pos1 += NB_POS_FCB_4T;
            pos2 += NB_POS_FCB_4T;
        }
    }

    pos[0] = pos1;
    pos[1] = pos2;

    return;
}

/*-------------------------------------------------------*
 * dec_3p_3N1()
 *
 * Decode 3 pulses with 3*N+1 bits:
 *-------------------------------------------------------*/

static void dec_3p_3N1(
    const long  index,   /* i:   quantization index    */
    const short N,       /* i:   nb. of bits           */
    const short offset,  /* i  : pulse position offset */
    short pos[]    /* o  : pulse position        */
)
{
    short j;
    long idx, mask;

    mask = ((1 << ((2 * N) - 1)) - 1);
    idx = index & mask;
    j = offset;

    if(((index >> ((2 * N) - 1)) & 1) == 1)
    {
        j += (1 << (N - 1));
    }

    dec_2p_2N1(idx, N - 1, j, pos);
    mask = ((1 << (N + 1)) - 1);
    idx = (index >> (2 * N)) & mask;
    dec_1p_N1(idx, N, offset, pos + 2);

    return;
}

/*-------------------------------------------------------*
 * dec_4p_4N1()
 *
 * Decode 4 pulses with 4*N+1 bits:
 *-------------------------------------------------------*/

static void dec_4p_4N1(
    const long index,     /* i  : quantization index    */
    const short N,        /* i  : nb. of bits           */
    const short offset,   /* i  : pulse position offset */
    short pos[]     /* o  : pulse position        */
)
{
    short j;
    long mask, idx;

    mask = ((1 << ((2 * N) - 1)) - 1);
    idx = index & mask;
    j = offset;

    if(((index >> ((2 * N) - 1)) & 1) == 1)
    {
        j += (1 << (N - 1));
    }

    dec_2p_2N1(idx, N - 1, j, pos);
    mask = ((1 << ((2 * N) + 1)) - 1);
    idx = (index >> (2 * N)) & mask;
    dec_2p_2N1(idx, N, offset, pos + 2);

    return;
}

/*-------------------------------------------------------*
 * dec_4p_4N()
 *
 * Decode 4 pulses with 4*N bits:
 *-------------------------------------------------------*/

static void dec_4p_4N(
    const long index,   /* i  : quantization index    */
    const short N,      /* i  : nb. of bits           */
    const short offset, /* i  : pulse position offset */
    short pos[]   /* o  : pulse position        */
)
{
    short j, n_1;

    n_1 = N - 1;
    j = offset + (1 << n_1);
    switch((index >> ((4 * N) - 2)) & 3)
    {
    case 0:
    {
        if(((index >> ((4 * n_1) + 1)) & 1) == 0)
        {
            dec_4p_4N1(index, n_1, offset, pos);
        }
        else
        {
            dec_4p_4N1(index, n_1, j, pos);
        }
        break;
    }
    case 1:
    {
        dec_1p_N1((index >> ((3 * n_1) + 1)), n_1, offset, pos);
        dec_3p_3N1(index, n_1, j, pos + 1);
        break;
    }
    case 2:
    {
        dec_2p_2N1((index >> ((2 * n_1) + 1)), n_1, offset, pos);
        dec_2p_2N1(index, n_1, j, pos + 2);
        break;
    }
    case 3:
    {
        dec_3p_3N1((index >> (n_1 + 1)), n_1, offset, pos);
        dec_1p_N1(index, n_1, j, pos + 3);
        break;
    }
    }

    return;
}

/*-------------------------------------------------------*
 * dec_5p_5N()
 *
 * Decode 5 pulses with 5*N bits:
 *-------------------------------------------------------*/

static void dec_5p_5N(
    const long index,    /* i  : quantization index    */
    const short N,       /* i  : nb. of bits           */
    const short offset,  /* i  : pulse position offset */
    short pos[]    /* o  : pulse position        */
)
{
    short  j, n_1;
    long  idx;

    n_1 = N - 1;
    j = offset + (1 << n_1);
    idx = (index >> ((2 * N) + 1));

    if(((index >> ((5 * N) - 1)) & 1) == 0)
    {
        dec_3p_3N1(idx, n_1, offset, pos);
        dec_2p_2N1(index, N, offset, pos + 3);
    }
    else
    {
        dec_3p_3N1(idx, n_1, j, pos);
        dec_2p_2N1(index, N, offset, pos + 3);
    }

    return;
}

/*-------------------------------------------------------*
 * dec_6p_6N2()
 *
 * Decode 6 pulses with 6*N+2 bits:
 *-------------------------------------------------------*/

static void dec_6p_6N2(
    const long index,     /* i  : quantization index    */
    const short N,        /* i  : nb. of bits           */
    const short offset,   /* i  : pulse position offset */
    short pos[]     /* o  : pulse position        */
)
{
    short j, n_1, offsetA, offsetB;

    n_1 = N - 1;
    j = offset + (1 << n_1);
    offsetA = offsetB = j;
    if(((index >> ((6 * N) - 5)) & 1) == 0)
    {
        offsetA = offset;
    }
    else
    {
        offsetB = offset;
    }

    switch((index >> ((6 * N) - 4)) & 3)
    {
    case 0:
    {
        dec_5p_5N(index >> N, n_1, offsetA, pos);
        dec_1p_N1(index, n_1, offsetA, pos + 5);
        break;
    }

    case 1:
    {
        dec_5p_5N(index >> N, n_1, offsetA, pos);
        dec_1p_N1(index, n_1, offsetB, pos + 5);
        break;
    }

    case 2:
    {
        dec_4p_4N(index >> ((2 * n_1) + 1), n_1, offsetA, pos);
        dec_2p_2N1(index, n_1, offsetB, pos + 4);
        break;
    }

    case 3:
    {
        dec_3p_3N1(index >> ((3 * n_1) + 1), n_1, offset, pos);
        dec_3p_3N1(index, n_1, j, pos + 3);
        break;
    }
    }

    return;
}

/*---------------------------------------------------------------------*
 * fcb_decode_class_all_p()
 *
 * Get the position number and the pulse number on every position
 *---------------------------------------------------------------------*/

static int fcb_decode_class_all_p(   /* o:   The index of pulse positions                          */
    int *code_index,                 /* i:   fcb index information                              */
    short sector_6p_num[],           /* o:   Number of pulses for each position                 */
    int pulse_num,                   /* i:   Number of pulses on a track.                       */
    int *pos_num                     /* o:   Number of positions which have pulses on the track.*/
)
{
    int i,j,k;
    int mn9;
    int pulse_pos_num;
    for (i=1; i<=pulse_num; i++)
    {
        if ((*code_index) < PI_offset[pulse_num][i])
        {
            break;
        }
    }

    (*code_index) -= PI_offset[pulse_num][i-1];

    pulse_pos_num = pulse_num - i + 2;
    j = (*code_index)>>pulse_pos_num;

    k = j/PI_select_table[16][pulse_pos_num];
    mn9 = j - k * PI_select_table[16][pulse_pos_num];
    if ( ( pulse_pos_num < pulse_num ) && ( pulse_pos_num > 1 ) )
    {
        for (i=0; i<pulse_pos_num; i++)
        {
            sector_6p_num[i] = 1;
        }
        sector_6p_num[k] ++;
    }
    else
    {
        if ( pulse_pos_num == 1 )
        {
            sector_6p_num[0] = (short)pulse_num;
        }
        else
        {
            for (i=0; i<pulse_num; i++)
            {
                sector_6p_num[i] = 1;
            }
        }
    }

    *pos_num = pulse_pos_num;

    return mn9;
}

/*---------------------------------------------------------------------*
 * fcb_decode_position()
 *
 * Decode the pulse position not same to the others
 *---------------------------------------------------------------------*/

static void fcb_decode_position(
    int index,                /* i:   position index information        */
    short pos_vector[],        /* o:   the positon vector                */
    int pos_num                /* i:   number of positions               */
)
{
    int  i,k,l;
    int  temp;

    k = index ;
    l = 0;
    temp = pos_num;
    for (i=0; i<pos_num-1; i++)
    {
        k = PI_select_table[16-l][temp] - k ;

        for (; PI_select_table[16 - l][temp] >= k; l+=2);

        if (k > PI_select_table[17-l][temp])
        {
            l--;
        }

        k = PI_select_table[17-l][temp--] - k ;
        pos_vector[i] = (short)(l-1);
    }
    pos_vector[i] = (short)(l+k);

    return;
}

/*---------------------------------------------------------------------*
 * fcb_decode_PI()
 *
 * decode fcb pulse index
 *---------------------------------------------------------------------*/

static int fcb_decode_PI(  /* o:   return pulse position number   */
    int code_index,        /* i:   fcb index information          */
    short sector_6p[],     /* o:   decoded pulse position         */
    int pulse_num          /* i:   pulse number for the track     */
)
{
    int i,l;
    int mn9;
    int pulse_pos_num;
    short sector_6p_temp[7], sector_6p_num_temp[7];
    short *sector_6p_ptr0;
    short *sector_6p_ptr1;

    /*get the position number and the pulse number on every position */
    mn9 = fcb_decode_class_all_p(&code_index, sector_6p_num_temp, pulse_num, &pulse_pos_num);

    /* rebuild the vector*/
    /* decode the pulse position not same to the others*/
    fcb_decode_position(mn9, sector_6p_temp, pulse_pos_num);
    for (i=pulse_pos_num-1; i>=0; i--)
    {
        sector_6p_temp[i] +=  ((code_index&0x1)<<4) ;
        code_index = code_index>>1;
    }

    /* decode the pulse position maybe some pulse position same to other pulse */
    sector_6p_ptr0 = &sector_6p[pulse_num];
    sector_6p_ptr1 = &sector_6p_temp[pulse_pos_num];
    for (i=0; i<pulse_pos_num; i++)
    {
        sector_6p_ptr1 -- ;
        for (l=0; l<sector_6p_num_temp[pulse_pos_num-1-i]; l++)
        {
            *--sector_6p_ptr0 = *sector_6p_ptr1 ;
        }
    }

    return pulse_pos_num;
}

/*---------------------------------------------------------------------*
 * Read FCB index                                                      *
 *---------------------------------------------------------------------*/

void D_ACELP_decode_43bit(unsigned short idxs[], float code[], int *pulsestrack)
{
    int ps[8];
    short pos[7];
    int joint_index;
    int joint_offset = 3611648;        /*offset for 3 pulses per track*/

    set_f( code, 0.0f, L_SUBFR );

    ps[3] = idxs[0] & 0x1ff;
    ps[2] = ( ( idxs[1] & 3 ) << 7) + ( idxs[0] >> 9 );
    joint_index = ( ( idxs[2] << 16 ) + idxs[1] ) >> 2;

    if ( joint_index >= joint_offset)
    {
        joint_index = joint_index - joint_offset;
    }

    ps[0] = joint_index/5472;
    ps[1] = joint_index - ps[0]*5472;
    fcb_decode_PI(ps[0], pos, 3);
    add_pulses(pos, pulsestrack[0], 0, code);
    fcb_decode_PI(ps[1], pos, 3);
    add_pulses(pos, pulsestrack[1], 1, code);

    dec_2p_2N1(ps[2], 4, 0, pos);
    add_pulses(pos, pulsestrack[2], 2, code);
    dec_2p_2N1(ps[3], 4, 0, pos);
    add_pulses(pos, pulsestrack[3], 3, code);

    return;
}

