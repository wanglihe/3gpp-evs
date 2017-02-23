/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <assert.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "stat_enc.h"
#include "stat_dec.h"
#include "rom_com.h"
#include "mime.h"




/*-------------------------------------------------------------------*
* pack_bit()
*
* insert a bit into packed octet
*-------------------------------------------------------------------*/
static void pack_bit(
    const Word16 bit,    /* i:   bit to be packed */
    UWord8 **pt,         /* i/o: pointer to octet array into which bit will be placed */
    UWord8 *omask        /* i/o: output mask to indicate where in the octet the bit is to be written */
)
{
    if (*omask == 0x80)
    {
        **pt = 0;
    }
    if (bit != 0)
    {
        **pt = **pt | *omask;
    }
    *omask >>= 1;
    if (*omask == 0)
    {
        *omask = 0x80;
        (*pt)++;
    }

    return;
}

/*-------------------------------------------------------------------*
* unpack_bit()
*
* unpack a bit from packed octet
*-------------------------------------------------------------------*/
static Word16 unpack_bit(
    UWord8 **pt,         /* i/o: pointer to octet array from which bit will be read */
    UWord8 *mask         /* i/o: mask to indicate the bit in the octet */
)
{
    Word16 bit;

    bit = (**pt & *mask) != 0;
    *mask >>= 1;
    if (*mask == 0)
    {
        *mask = 0x80;
        (*pt)++;
    }

    return bit;
}

/*-------------------------------------------------------------------*
* rate2AMRWB_IOmode()
*
* lookup AMRWB IO mode
*-------------------------------------------------------------------*/

static Word16 rate2AMRWB_IOmode(
    Word32 rate                    /* i: bit rate */
)
{
    switch ( rate )
    {
    /* EVS AMR-WB IO modes */
    case SID_1k75      :
        return AMRWB_IO_SID;
    case ACELP_6k60    :
        return AMRWB_IO_6600;
    case ACELP_8k85    :
        return AMRWB_IO_8850;
    case ACELP_12k65   :
        return AMRWB_IO_1265;
    case ACELP_14k25   :
        return AMRWB_IO_1425;
    case ACELP_15k85   :
        return AMRWB_IO_1585;
    case ACELP_18k25   :
        return AMRWB_IO_1825;
    case ACELP_19k85   :
        return AMRWB_IO_1985;
    case ACELP_23k05   :
        return AMRWB_IO_2305;
    case ACELP_23k85   :
        return AMRWB_IO_2385;
    default:
        return -1;
    }
}

/*-------------------------------------------------------------------*
* rate2EVSmode()
*
* lookup EVS mode
*-------------------------------------------------------------------*/
static Word16 rate2EVSmode(
    Word32 rate                    /* i: bit rate */
)
{
    switch ( rate )
    {
    /* EVS Primary modes */
    case FRAME_NO_DATA :
        return NO_DATA;
    case SID_2k40      :
        return PRIMARY_SID;
    case PPP_NELP_2k80 :
        return PRIMARY_2800;
    case ACELP_7k20    :
        return PRIMARY_7200;
    case ACELP_8k00    :
        return PRIMARY_8000;
    case ACELP_9k60    :
        return PRIMARY_9600;
    case ACELP_13k20   :
        return PRIMARY_13200;
    case ACELP_16k40   :
        return PRIMARY_16400;
    case ACELP_24k40   :
        return PRIMARY_24400;
    case ACELP_32k     :
        return PRIMARY_32000;
    case ACELP_48k     :
        return PRIMARY_48000;
    case ACELP_64k     :
        return PRIMARY_64000;
    case HQ_96k        :
        return PRIMARY_96000;
    case HQ_128k       :
        return PRIMARY_128000;
    default            :
        return rate2AMRWB_IOmode(rate);
    }
}

/*-------------------------------------------------------------------*
 * push_indice()
 *
 * Push a new indice into the buffer
 *-------------------------------------------------------------------*/

void push_indice(
    Encoder_State *st,            /* i/o: encoder state structure */
    short id,            /* i  : ID of the indice */
    unsigned short value,         /* i  : value of the quantized indice */
    short nb_bits        /* i  : number of bits used to quantize the indice */
)
{
    short i;


    if ( st->last_ind == id )
    {
        /* indice with the same name as the previous one */
        i = st->next_ind;
    }
    else
    {
        /* new indice - find an empty slot in the list */
        i = id;
        while (st->ind_list[i].nb_bits != -1)
        {
            i++;
        }
    }

    /* store the new indice in the list */
    st->ind_list[i].value = value;
    st->ind_list[i].nb_bits = nb_bits;

    /* updates */
    st->next_ind = i + 1;
    st->last_ind = id;
    st->nb_bits_tot += nb_bits;

    return;
}

/*-------------------------------------------------------------------*
 * push_next_indice()
 *
 * Push a new indice into the buffer at the next position
 *-------------------------------------------------------------------*/

void push_next_indice(
    Encoder_State *st,           /* i/o: encoder state structure */
    unsigned short value,        /* i  : value of the quantized indice */
    short nb_bits       /* i  : number of bits used to quantize the indice */
)
{

    /* store the values in the list */
    st->ind_list[st->next_ind].value   = value;
    st->ind_list[st->next_ind].nb_bits = nb_bits;
    st->next_ind++;

    /* update the total number of bits already written */
    st->nb_bits_tot += nb_bits;

    return;
}

/*-------------------------------------------------------------------*
 * push_next_bits()
 * Push a bit buffer into the buffer at the next position
 *-------------------------------------------------------------------*/

void push_next_bits(
    Encoder_State *st,         /* i/o: encoder state structure */
    int bits[],      /* i  : bit buffer to pack, sequence of single bits */
    short nb_bits      /* i  : number of bits to pack */
)
{
    unsigned short code;
    int i, nb_bits_m15;
    Indice *ptr;

    ptr = &st->ind_list[st->next_ind];
    nb_bits_m15 = nb_bits - 15;

    for (i=0; i<nb_bits_m15; i += 16)
    {
        code = (unsigned short)((bits[i] << 15) | ((bits[i+1] << 14) | ((bits[i+2] << 13) | ((bits[i+3] << 12) |
                                ((bits[i+4] << 11) | ((bits[i+5] << 10) | ((bits[i+6] << 9) | ((bits[i+7] << 8) |
                                        ((bits[i+8] << 7) | ((bits[i+9] << 6) | ((bits[i+10] << 5) | ((bits[i+11] << 4) |
                                                ((bits[i+12] << 3) | ((bits[i+13] << 2) | ((bits[i+14] << 1) | bits[i+15])))))))))))))));

        ptr->value   = code;
        ptr->nb_bits = 16;
        ++ptr;
    }
    for (; i<nb_bits; ++i)
    {
        ptr->value   = bits[i];
        ptr->nb_bits = 1;
        ++ptr;
    }
    st->next_ind = (int)(ptr - st->ind_list);
    st->nb_bits_tot = st->nb_bits_tot + nb_bits;
}

/*-------------------------------------------------------------------*
 * get_next_indice()
 *
 * Get the next indice from the buffer
 *-------------------------------------------------------------------*/

unsigned short get_next_indice(   /* o  : value of the indice */
    Decoder_State *st,                    /* i/o: decoder state structure */
    short nb_bits                 /* i  : number of bits that were used to quantize the indice */
)
{
    unsigned short value;
    short i;

    assert(nb_bits <= 16);

    /* detect corrupted bitstream */
    if( st->next_bit_pos + nb_bits > st->total_brate/50 )
    {
        st->BER_detect = 1;
        return(0);
    }

    value = 0;
    for (i = 0; i < nb_bits; i++)
    {
        value <<= 1;
        value += st->bit_stream[st->next_bit_pos+i];
    }

    /* update the position in the bitstream */
    st->next_bit_pos += nb_bits;

    return value;
}

/*-------------------------------------------------------------------*
 * get_next_indice_1()
 *
 * Get the next 1-bit indice from the buffer
 *-------------------------------------------------------------------*/

unsigned short get_next_indice_1(     /* o  : value of the indice */
    Decoder_State *st                         /* i/o: decoder state structure */
)
{

    /* detect corrupted bitstream */
    if( ( st->next_bit_pos + 1 > st->total_brate/50 && st->codec_mode == MODE1 ) ||
            ( (st->next_bit_pos + 1 > st->total_brate/50 + (2*8) ) && st->codec_mode == MODE2 ) /* add two zero bytes for arithmetic coder flush */
      )
    {
        st->BER_detect = 1;
        return(0);
    }

    return st->bit_stream[st->next_bit_pos++];
}

/*-------------------------------------------------------------------*
 * get_next_indice_tmp()
 *
 * update the total number of bits and the position in the bitstream
 *-------------------------------------------------------------------*/

void get_next_indice_tmp(
    Decoder_State *st,                     /* o  : decoder state structure */
    short nb_bits                  /* i  : number of bits that were used to quantize the indice */
)
{
    /* update the position in the bitstream */
    st->next_bit_pos += nb_bits;

}

/*-------------------------------------------------------------------*
 * get_indice()
 *
 * Get indice at specific position in the buffer
 *-------------------------------------------------------------------*/

unsigned short get_indice(    /* o  : value of the indice */
    Decoder_State *st,        /* i/o: decoder state structure */
    short pos,        /* i  : absolute position in the bitstream (update after the read) */
    short nb_bits     /* i  : number of bits that were used to quantize the indice */
)
{
    unsigned short value;
    int i;

    assert(nb_bits <= 16);

    /* detect corrupted bitstream */
    if( pos + nb_bits > st->total_brate/50 )
    {
        st->BER_detect = 1;
        return(0);
    }

    value = 0;
    for (i = 0; i < nb_bits; i++)
    {
        value <<= 1;
        value += st->bit_stream[pos+i];
    }

    return value;
}

/*-------------------------------------------------------------------*
 * get_indice_1()
 *
 * Get a 1-bit indice at specific position in the buffer
 *-------------------------------------------------------------------*/

unsigned short get_indice_1(  /* o  : value of the indice */
    Decoder_State *st,        /* i/o: decoder state structure */
    short pos         /* i  : absolute position in the bitstream (update after the read) */
)
{
    /* detect corrupted bitstream */
    if( pos+1 > st->total_brate/50 )
    {
        st->BER_detect = 1;
        return(0);
    }

    return st->bit_stream[pos];
}

/*-------------------------------------------------------------------*
 * reset_indices_enc()
 *
 * Reset the buffer of encoder indices
 *-------------------------------------------------------------------*/

void reset_indices_enc(
    Encoder_State *st
)
{
    short i;

    st->nb_bits_tot = 0;
    st->next_ind = 0;
    st->last_ind = -1;

    for (i=0; i<MAX_NUM_INDICES; i++)
    {
        st->ind_list[i].nb_bits = -1;
    }

    return;
}

/*-------------------------------------------------------------------*
 * reset_indices_dec()
 *
 * Reset the buffer of decoder indices
 *-------------------------------------------------------------------*/

void reset_indices_dec(
    Decoder_State *st
)
{
    st->next_bit_pos = 0;

    return;
}

/*-------------------------------------------------------------------*
* write_indices()
*
* Write the buffer of indices to a file
*-------------------------------------------------------------------*/

void write_indices(
    Encoder_State *st,  /* i/o: encoder state structure                                     */
    FILE *file,         /* i  : output bitstream file                                       */
    UWord8 *pFrame,     /* i: byte array with bit packet and byte aligned coded speech data */
    Word16 pFrame_size  /* i: size of the binary encoded access unit [bits]                 */
)
{
    short i, k;
    unsigned short stream[2+MAX_BITS_PER_FRAME], *pt_stream;
    int mask;
    short value, nb_bits;
    UWord8 header;

    if( st->bitstreamformat == G192 )
    {
        /*-----------------------------------------------------------------*
         * Encode Sync Header and Frame Length
         *-----------------------------------------------------------------*/

        pt_stream = stream;
        for (i=0; i<(2 + MAX_BITS_PER_FRAME); ++i)
        {
            stream[i] = 0;
        }
        *pt_stream++ = SYNC_GOOD_FRAME;
        *pt_stream++ = st->nb_bits_tot;

        /*----------------------------------------------------------------*
         * Bitstream packing (conversion of individual indices into a serial stream)
         * Writing the serial stream into file
         * Clearing of indices
         *----------------------------------------------------------------*/

        for (i=0; i<MAX_NUM_INDICES; i++)
        {
            value = st->ind_list[i].value;
            nb_bits = st->ind_list[i].nb_bits;
            if (nb_bits != -1)
            {
                /* mask from MSB to LSB */
                mask = 1 << (nb_bits - 1);

                /* write bit by bit */
                for (k=0; k < nb_bits; k++)
                {
                    if ( value & mask )
                    {
                        *pt_stream++ = G192_BIN1;
                    }
                    else
                    {
                        *pt_stream++ = G192_BIN0;
                    }

                    mask >>= 1;
                }
            }
        }

    }
    else
    {
        /* Create and write ToC header */
        /*  qbit always  set to  1 on encoder side  for AMRWBIO , no qbit in use for EVS,  but set to 0(bad)  */
        header = (UWord8)(st->Opt_AMR_WB << 5 | st->Opt_AMR_WB << 4 | rate2EVSmode(st->nb_bits_tot * 50));
        fwrite( &header, sizeof(UWord8), 1, file );
        /* Write speech bits */
        fwrite( pFrame, sizeof(UWord8), (pFrame_size + 7)>>3, file );
    }
    /* Clearing of indices */
    for (i=0; i<MAX_NUM_INDICES; i++)
    {
        st->ind_list[i].nb_bits = -1;
    }

    if( st->bitstreamformat == G192 )
    {
        /* write the serial stream into file */
        fwrite( stream, sizeof(unsigned short), 2+stream[1], file );
    }


    /* reset index pointers */
    st->nb_bits_tot = 0;
    st->next_ind = 0;
    st->last_ind = -1;

    return;
}

/*-------------------------------------------------------------------*
 * indices_to_serial()
 *
 * pack indices into serialized payload format
 *-------------------------------------------------------------------*/

void indices_to_serial(
    const Encoder_State *st,                  /* i: encoder state structure */
    UWord8 *pFrame,           /* o: byte array with bit packet and byte aligned coded speech data */
    Word16 *pFrame_size       /* o: size of the binary encoded access unit [bits] */
)
{
    Word16 i, k, j;
    Word16 cmi = 0, core_mode=0;
    Word32  mask;
    Word16 amrwb_bits[(ACELP_23k85 / 50)];
    UWord8 omask= 0x80;
    UWord8 *pt_pFrame=pFrame;

    if ( st->Opt_AMR_WB )
    {
        cmi = rate2EVSmode(st->total_brate);
        core_mode = rate2EVSmode(st->nb_bits_tot * 50);

        j=0;
        for (i=0; i<MAX_NUM_INDICES; i++)
        {
            if (st->ind_list[i].nb_bits != -1)
            {
                /* mask from MSB to LSB */
                mask = 1 << (st->ind_list[i].nb_bits - 1);

                /* temporarily save bit */
                for (k=0; k < st->ind_list[i].nb_bits; k++)
                {
                    amrwb_bits[j++] = (st->ind_list[i].value & mask) > 0;
                    mask >>= 1;
                }
            }
        }
    }

    *pFrame_size = st->nb_bits_tot;

    /*----------------------------------------------------------------*
    * Bitstream packing (conversion of individual indices into a serial stream)
    *----------------------------------------------------------------*/
    j=0;
    for (i=0; i<MAX_NUM_INDICES; i++)
    {
        if (st->ind_list[i].nb_bits != -1)
        {
            /* mask from MSB to LSB */
            mask = 1 << (st->ind_list[i].nb_bits - 1);

            /* write bit by bit */
            for (k=0; k < st->ind_list[i].nb_bits; k++)
            {
                if (st->Opt_AMR_WB )
                {
                    pack_bit(amrwb_bits[sort_ptr[core_mode][j++]], &pt_pFrame, &omask);
                }
                else
                {
                    pack_bit(st->ind_list[i].value & mask, &pt_pFrame, &omask);
                    j++;
                }
                mask >>= 1;
            }
        }
    }

    if ( st->Opt_AMR_WB && core_mode == AMRWB_IO_SID)    /* SID frame */
    {
        /* insert STI bit and CMI */
        pack_bit(1, &pt_pFrame, &omask);
        for (mask=0x08; mask>0; mask >>= 1)
        {
            pack_bit(cmi & mask, &pt_pFrame, &omask);
        }
    }

    return;
}

/*-------------------------------------------------------------------*
 * indices_to_serial_generic()
 *
 * pack indices into serialized payload format
 *-------------------------------------------------------------------*/

void indices_to_serial_generic(
    const Indice *ind_list,         /* i: indices list */
    const Word16  num_indices,      /* i: number of indices to write */
    UWord8 *pFrame,           /* o: byte array with bit packet and byte aligned coded speech data */
    Word16 *pFrame_size       /* i/o: number of bits in the binary encoded access unit [bits] */
)
{
    Word16 i, k, j;
    Word32 mask;
    UWord8 omask;
    UWord8 *pt_pFrame = pFrame;
    Word16 nb_bits_tot = 0;

    omask = (0x80 >> (*pFrame_size & 0x7));
    pt_pFrame += *pFrame_size >> 3;

    /*----------------------------------------------------------------*
    * Bitstream packing (conversion of individual indices into a serial stream)
    *----------------------------------------------------------------*/
    j=0;
    for (i=0; i<num_indices; i++)
    {
        if (ind_list[i].nb_bits != -1)
        {
            /* mask from MSB to LSB */
            mask = 1 << (ind_list[i].nb_bits - 1);

            /* write bit by bit */
            for (k=0; k < ind_list[i].nb_bits; k++)
            {
                pack_bit(ind_list[i].value & mask, &pt_pFrame, &omask);
                j++;
                mask >>= 1;
            }
            nb_bits_tot += ind_list[i].nb_bits;
        }
    }

    *pFrame_size += nb_bits_tot;

    return;
}

/*-------------------------------------------------------------------*
 * decoder_selectCodec()
 *
 *
 *-------------------------------------------------------------------*/

static void decoder_selectCodec(
    Decoder_State *st,            /* i/o: decoder state structure                */
    const long total_brate,    /* i  : total bitrate                          */
    const short bit0           /* i  : first bit                              */
)
{
    /* set the AMR-WB IO flag */
    if( total_brate == SID_1k75 ||
            total_brate == ACELP_6k60  || total_brate == ACELP_8k85  || total_brate == ACELP_12k65 ||
            total_brate == ACELP_14k25 || total_brate == ACELP_15k85 || total_brate == ACELP_18k25 ||
            total_brate == ACELP_19k85 || total_brate == ACELP_23k05 || total_brate == ACELP_23k85 )
    {
        st->Opt_AMR_WB = 1;
    }
    else if ( total_brate != FRAME_NO_DATA )
    {
        st->Opt_AMR_WB = 0;
    }

    if ( st->Opt_AMR_WB )
    {
        st->codec_mode = MODE1;
    }
    else
    {
        switch ( total_brate )
        {
        case 0:
            st->codec_mode = st->last_codec_mode;
            break;
        case 2400:
            st->codec_mode = st->last_codec_mode;
            break;
        case 2800:
            st->codec_mode = MODE1;
            break;
        case 7200:
            st->codec_mode = MODE1;
            break;
        case 8000:
            st->codec_mode = MODE1;
            break;
        case 9600:
            st->codec_mode = MODE2;
            break;
        case 13200:
            st->codec_mode = MODE1;
            break;
        case 16400:
            st->codec_mode = MODE2;
            break;
        case 24400:
            st->codec_mode = MODE2;
            break;
        case 32000:
            st->codec_mode = MODE1;
            break;
        case 48000:
            st->codec_mode = MODE2;
            break;
        case 64000:
            st->codec_mode = MODE1;
            break;
        case 96000:
            st->codec_mode = MODE2;
            break;
        case 128000:
            st->codec_mode = MODE2;
            break;
        default :   /* validate that total_brate (derived from RTP packet or a file header) is one of the defined bit rates  */
            st->codec_mode = st->last_codec_mode;
            st->bfi = 1;
            break;
        }
    }

    if( st->ini_frame == 0 )
    {
        if( st->codec_mode == -1 )
        {
            st->codec_mode = MODE1;
        }
        st->last_codec_mode = st->codec_mode;
    }

    /* set SID/CNG type */
    if( total_brate == SID_2k40 )
    {
        if( bit0 == G192_BIN0 )
        {
            st->cng_type = LP_CNG;

            /* force MODE1 when selecting LP_CNG */
            st->codec_mode = MODE1;
        }
        else
        {
            st->cng_type = FD_CNG;
            if ( st->last_codec_mode == MODE2 && st->last_total_brate == ACELP_13k20 )
            {
                st->codec_mode = MODE1;
            }
        }
        st->last_cng_type = st->cng_type;     /* CNG type switching at the first correctly received SID frame */
    }


    return;
}

/*-------------------------------------------------------------------*
 * dec_prm_core()
 *
 *
 *-------------------------------------------------------------------*/

void dec_prm_core(
    Decoder_State *st
)
{
    int n, frame_size_index = -1;

    st->core = -1;

    if (st->total_brate == FRAME_NO_DATA)
    {
        st->m_frame_type = ZERO_FRAME;
    }
    else if (st->total_brate == SID_2k40)
    {
        st->m_frame_type = SID_FRAME;
    }
    else
    {
        st->m_frame_type = ACTIVE_FRAME;
        for (n=0; n<FRAME_SIZE_NB; ++n)
        {
            if (FrameSizeConfig[n].frame_bits == st->total_brate/50)
            {
                frame_size_index =  n;
                break;
            }
        }

        /* Get bandwidth mode */
        st->bwidth = get_next_indice(st, FrameSizeConfig[frame_size_index].bandwidth_bits);
        st->bwidth += FrameSizeConfig[frame_size_index].bandwidth_min;
        if (st->bwidth > FB)
        {
            st->bwidth = FB;
            st->BER_detect = 1;
        }

        if (st->bwidth > SWB && st->total_brate < ACELP_16k40)
        {
            st->bwidth = SWB;
            st->BER_detect = 1;
        }
        /* Skip reserved bit */
        get_next_indice_tmp(st, FrameSizeConfig[frame_size_index].reserved_bits);

        if (get_next_indice_1(st)) /* TCX */
        {
            if (get_next_indice_1(st))
            {
                st->core = HQ_CORE;
            }
            else
            {
                st->core = TCX_20_CORE;
            }
        }
        else /* ACELP */
        {
            st->core = ACELP_CORE;
        }
    }

    return;
}

/*-----------------------------------------------------------------*
 * decision_matrix_core_dec()
 *
 * Read core mode signalling bits from the bitstream
 * Set st->core, and st->bwidth if signalled together with the core.
 *-----------------------------------------------------------------*/

static void decision_matrix_core_dec(
    Decoder_State *st                 /* i/o: decoder state structure                   */
)
{
    short start_idx;
    long ind;
    short nBits;

    assert(st->bfi != 1);

    st->core = -1;
    st->bwidth = -1;

    if ( st->total_brate == FRAME_NO_DATA || st->total_brate == SID_2k40 )
    {
        st->core = ACELP_CORE;
    }
    /* SC-VBR */
    else if ( st->total_brate == PPP_NELP_2k80 )
    {
        st->core = ACELP_CORE;
        return;
    }

    /*---------------------------------------------------------------------*
     * ACELP/HQ core selection
     *---------------------------------------------------------------------*/

    if ( st->total_brate < ACELP_24k40 )
    {
        st->core = ACELP_CORE;
    }
    else if ( st->total_brate >= ACELP_24k40 && st->total_brate <= ACELP_64k )
    {
        /* read the ACELP/HQ core selection bit */
        st->core = ((short) get_next_indice( st, 1 ))*HQ_CORE;
    }
    else
    {
        st->core = HQ_CORE;
    }

    /*-----------------------------------------------------------------*
     * Read ACELP signalling bits from the bitstream
     *-----------------------------------------------------------------*/

    if ( st->core == ACELP_CORE )
    {
        /* find the section in the ACELP signalling table corresponding to bitrate */
        start_idx = 0;
        while ( acelp_sig_tbl[start_idx] != st->total_brate )
        {
            start_idx++;
        }

        /* skip the bitrate */
        start_idx += 1;

        /* retrieve the number of bits */
        nBits = (short) acelp_sig_tbl[start_idx++];

        /* retrieve the signalling indice */
        ind = acelp_sig_tbl[start_idx + get_next_indice( st, nBits )];
        st->bwidth = (ind >> 3) & 0x7;

        /* convert signalling indice into signalling information */
        if ( (ind & 0x7) == LR_MDCT )
        {
            st->core = HQ_CORE;
        }
    }

    /*-----------------------------------------------------------------*
     * Read HQ signalling bits from the bitstream
     * Set HQ core type
     *-----------------------------------------------------------------*/

    if ( st->core == HQ_CORE )
    {
        /* read the HQ/TCX core switching flag */
        if ( get_next_indice( st, 1 ) )
        {
            st->core = TCX_20_CORE;
        }

        /* For TCX: read/set band-width (needed for different I/O sampling rate support) */
        if( st->core == TCX_20_CORE && st->total_brate > ACELP_16k40 )
        {
            ind = get_next_indice( st, 2 );

            if( ind == 0 )
            {
                st->bwidth = NB;
            }
            else if( ind == 1 )
            {
                st->bwidth = WB;
            }
            else if( ind == 2 )
            {
                st->bwidth = SWB;
            }
            else
            {
                st->bwidth = FB;
            }
        }
    }

    return;
}

/*-------------------------------------------------------------------*
 * mdct_switching_dec()
 *
 * Set up MDCT core switching if indicated in the bit stream
 *-------------------------------------------------------------------*/

static void mdct_switching_dec(
    Decoder_State *st           /* i/o: decoder state structure     */
)
{
    if (st->Opt_AMR_WB)
    {
        return;
    }

    if (st->total_brate == ACELP_13k20 || st->total_brate == ACELP_32k)
    {
        st->mdct_sw_enable = MODE1;
    }
    else if (ACELP_16k40 <= st->total_brate && st->total_brate <= ACELP_24k40)
    {
        st->mdct_sw_enable = MODE2;
    }

    if ( st->codec_mode == MODE1 && st->mdct_sw_enable == MODE1 )
    {
        /* Read ahead core mode signaling */
        short next_bit_pos_save = st->next_bit_pos;
        short core_save         = st->core;
        short bwidth_save       = st->bwidth;

        decision_matrix_core_dec(st); /* sets st->core */

        if (st->core == TCX_20_CORE)
        {
            /* Trigger TCX */
            st->codec_mode = MODE2;
            st->mdct_sw = MODE1;
        }
        else
        {
            /* Rewind bitstream */
            st->next_bit_pos = next_bit_pos_save;
            if (st->bfi)
            {
                st->core   = core_save;
                st->bwidth = bwidth_save;
            }
        }
    }
    else if (st->codec_mode == MODE2 && st->mdct_sw_enable == MODE2)
    {
        /* Read ahead core mode signaling */
        short next_bit_pos_save = st->next_bit_pos;
        short core_save         = st->core;
        short bwidth_save       = st->bwidth;

        dec_prm_core(st); /* sets st->core */

        if (st->core == HQ_CORE)
        {
            /* Trigger HQ_CORE */
            st->codec_mode = MODE1;
            st->mdct_sw = MODE2;
        }
        else
        {
            /* Rewind bitstream */
            st->next_bit_pos = next_bit_pos_save;
            if (st->bfi)
            {
                st->core   = core_save;
            }

            /* always reset bwidth, to not interfere with BER logic */
            st->bwidth = bwidth_save;
        }
    }

    return;
}



/*-------------------------------------------------------------------*
 * read_indices()
 *
 * Read indices from an ITU-T G.192 bitstream to the buffer
 * Simulate packet losses by inserting frame erasures
 *-------------------------------------------------------------------*/

short read_indices(                   /* o  : 1 = reading OK, 0 = problem            */
    Decoder_State *st,                /* i/o: decoder state structure                */
    FILE *file,              /* i  : bitstream file                         */
    const short rew_flag            /* i  : rewind flag (rewind file after reading)*/
)
{
    short k;
    unsigned short utmp, stream[2+MAX_BITS_PER_FRAME], *pt_stream;
    unsigned short *bit_stream_ptr;
    unsigned short num_bits;
    long total_brate;
    short curr_ft_good_sp, curr_ft_bad_sp;
    short g192_sid_first,sid_upd_bad, sid_update;
    short speech_bad, speech_lost;
    unsigned short num_bits_read;

    st->bfi = 0;
    st->BER_detect = 0;
    st->mdct_sw_enable = 0;
    st->mdct_sw = 0;
    reset_indices_dec( st );

    /* read the Sync Header field from the bitstream */
    /* in case rew_flag is set, read until first good frame is encountered */
    do
    {
        /* read the Sync header */
        if ( fread( &utmp, sizeof(unsigned short), 1, file ) != 1 )
        {
            if( ferror( file ) )
            {
                /* error during reading */
                fprintf(stderr, "\nError reading the bitstream !");
                exit(-1);
            }
            else
            {
                /* end of file reached */
                return 0;
            }
        }

        /* set the BFI indicator according the value of Sync Header */
        if ( utmp == SYNC_BAD_FRAME )
        {
            st->bfi = 1;
        }



        else
        {
            st->bfi = 0;
        }

        /* read the Frame Length field from the bitstream */
        if ( fread( &num_bits, sizeof(unsigned short), 1, file ) != 1 )
        {
            if( ferror( file ) )
            {
                /* error during reading */
                fprintf(stderr, "\nError reading the bitstream !");
                exit(-1);
            }
            else
            {
                /* end of file reached */
                return 0;
            }
        }

        /* convert the frame length to total bitrate */
        total_brate = (long)(num_bits * 50);

        /* read ITU-T G.192 serial stream of indices from file to the local buffer */
        /* Validate that the G.192 length is within the defined  bit rate range
        to not allow writing past the end of the "stream" buffer  */
        if( num_bits > MAX_BITS_PER_FRAME )
        {
            fprintf(stderr, "\nError, too large G.192 frame (size(%d))! Exiting ! \n", num_bits);
            exit(-1);
        }

        /*  verify that a  valid  num bits value  is present in the G.192 file */
        /*  only AMRWB or EVS bit rates  or 0(NO DATA)  are  allowed  in G.192 file frame reading  */
        if( rate2EVSmode(total_brate) < 0 ) /* negative value means that a valid rate was not found */
        {
            fprintf(stderr, "\nError, illegal bit rate (%ld) in  the  G.192 frame ! Exiting ! \n", total_brate);
            exit(-1);
        }
        pt_stream = stream;
        num_bits_read = (unsigned short) fread( pt_stream, sizeof(unsigned short), num_bits, file );
        if( num_bits_read != num_bits )
        {
            fprintf(stderr, "\nError, invalid number of bits read ! Exiting ! \n");
            exit(-1);
        }

    }
    while ( rew_flag && (st->bfi || total_brate < 2800) );

    /* G.192 RX DTX handler*/
    if( !rew_flag )
    {
        /* handle SID_FIRST, SID_BAD, SPEECH_LOST,  NO_DATA as properly  as possible for the ITU-T  G.192 format  */
        /* (total_brate, bfi , st_CNG)   =  rx_handler(received frame type, [previous frame type],  past CNG state, past core) */
        curr_ft_good_sp = 0;
        curr_ft_bad_sp  = 0;

        if( total_brate > SID_2k40 )
        {
            if( st->bfi == 0 )
            {
                curr_ft_good_sp = 1;
            }
            else
            {
                curr_ft_bad_sp = 1;
            }
        }

        sid_update = 0;
        sid_upd_bad = 0;
        if( total_brate == SID_1k75 || total_brate == SID_2k40 )
        {
            if( st->bfi == 0 )
            {
                sid_update = 1;
            }
            else
            {
                sid_upd_bad = 1;   /* this frame type may happen in ETSI/3GPP CS cases ,  a corrupt sid frames  */
            }
        }

        /* AMRWB  26.173 G.192  file reader (read_serial)  does not declare/use SID_BAD ft,
            it declares every bad synch marked frame initially as  a lost_speech frame,
           and then the RXDTX handler CNG state decides the decoding mode CNG/SPEECH.
           While In the AMRWB ETSI/3GPP format eid a CRC error in a detected  SID_UPDATE frames triggers SID_BAD.

            Here we inhibit use of the SID-length info, even though it is available in the G.192 file format after STL/EID-XOR .
         */
        if ( sid_upd_bad )
        {
            sid_upd_bad     = 0;
            total_brate     = FRAME_NO_DATA ;         /* treat SID_BAD  as a  stolen signaling frame --> SPEECH LOST */

        }

        g192_sid_first = 0;
        if( st->core == AMR_WB_CORE && st->prev_ft_speech && total_brate == FRAME_NO_DATA && st->bfi == 0 )
        {
            g192_sid_first = 1;   /*  SID_FIRST detected for previous AMRWB/AMRWBIO  active frames only  */
            /*
            It is not possible to perfectly simulate rate switching conditions EVS->AMRWBIO  where:
            the very first SID_FIRST detection is based on a past EVS active frame
            and  a  good length 0  "SID_FIRST"(NO_DATA)   frame is sent in AMRWBIO,
            , due to the one frame state memory in the AMRWB legacy  G.192 SID_FIRST encoding
            */
        }

        speech_bad = 0;
        if( total_brate > SID_2k40 && st->bfi != 0 ) /*   CS-type of CRC failure   frame */
        {
            speech_bad = 1;        /* initial ft assumption, CNG_state decides what to do */
        }

        speech_lost = 0;
        if( total_brate == 0 && st->bfi != 0 ) /*  unsent  NO_DATA or stolen NO_DATA/signaling  frame  */
        {
            speech_lost = 1;       /* initial ft assumption, CNG_state decides what to do */
        }

        /* Do not allow decoder to enter CNG-synthesis for  any instantly  received  GOOD+LENGTH==0  frame
           as this frame was never transmitted, one  can not know it is good and has a a length of zero ) */

        if( st->CNG != 0 )
        {
            /* We were in CNG synthesis  */
            if( curr_ft_good_sp != 0  )
            {
                /* only a good speech frame makes you leave CNG synthesis */
                st->CNG = 0;
            }
        }
        else
        {
            /* We were in SPEECH synthesis  */
            /* only a received/detected SID frame can make the decoder enter into CNG synthsis  */
            if( g192_sid_first || sid_update || sid_upd_bad )
            {
                st->CNG = 1;
            }
        }

        /* set bfi, total_brate pair  for proper decoding  */
        /*  handle the  G.192   _simulated_ untransmitted NO_DATA frame,  setting  for decoder  SPEECH synthesis  */
        if ( (st->CNG==0)  &&  (total_brate==0  && st->bfi == 0 ) )
        {
            st->bfi= 1;       /*  SPEECH PLC code will now become active as in a real system */
            /* total_brate= 0  */
        }

        /* handle bad/lost speech frame(and CS bad sid frame) in the decoders CNG synthesis settings pair (total_brate, bfi) */
        if( ((st->CNG != 0) && ( (speech_bad != 0) || (speech_lost != 0) ))  || /* SP_BAD or SPEECH_LOST)   --> stay in CNG */
                ( sid_upd_bad != 0 ))                                               /* SID_UPD_BAD              --> start CNG */
        {
            st->bfi = 0;         /* bfi=0 needed to activate CNG code */
            total_brate = 0;
        }
        /* update for next frame's G.192 file format's  odd SID_FIRST detection (primarily for AMRWBIO)  */
        st->prev_ft_speech = ((curr_ft_good_sp != 0) || (curr_ft_bad_sp != 0));

        /*   st->total brate= total_brate ;   updated in a good frame below */
    } /* rew_flag */

    if ( st->bfi == 0 && !rew_flag )
    {
        /* select Mode 1 or Mode 2 */
        decoder_selectCodec( st, total_brate, *pt_stream );
    }

    /* in case rew_flag is set, rewind the file and return */
    /* (used in io_enc() to print out info about technologies and to initialize the codec) */
    if ( rew_flag )
    {
        rewind( file );
        st->total_brate = total_brate;
        return 1;
    }

    /* GOOD frame */
    if ( st->bfi == 0  )
    {
        /* GOOD frame - convert ITU-T G.192 words to short values */
        bit_stream_ptr = st->bit_stream;

        for ( k = 0; k < (short)num_bits; ++k )
        {
            *bit_stream_ptr++ = (*pt_stream++ == G192_BIN1);
        }

        /* add two zero bytes for arithmetic coder flush */
        for (k=0; k<8*2; ++k)
        {
            *bit_stream_ptr++ = 0;
        }
        /* a change of the total bitrate should not be
           known to the decoder, if the received frame was
           lost */
        st->total_brate = total_brate;

        mdct_switching_dec(st);
    }

    return 1;
}



/*-------------------------------------------------------------------*
 * read_indices_mime()
 *
 * Read indices from MIME formatted bitstream to the buffer
 *-------------------------------------------------------------------*/

Word16 read_indices_mime(                /* o  : 1 = reading OK, 0 = problem            */
    Decoder_State *st,                   /* i/o: decoder state structure                */
    FILE *file,                          /* i  : bitstream file                         */
    Word16 rew_flag                      /* i  : rewind flag (rewind file after reading)*/
)
{
    Word16 k, isAMRWB_IOmode, cmi, core_mode = -1, qbit,sti;
    UWord8 header;
    UWord8 pFrame[(MAX_BITS_PER_FRAME + 7) >> 3];
    UWord8 mask= 0x80, *pt_pFrame=pFrame;
    UWord16 *bit_stream_ptr;
    Word16 num_bits;
    Word32 total_brate;
    Word16 curr_ft_good_sp;
    Word16 amrwb_sid_first, sid_upd_bad, sid_update;
    Word16 speech_bad, speech_lost;
    Word16 no_data;
    Word16 num_bytes_read;

    st->BER_detect = 0;
    st->bfi = 0;
    st->mdct_sw_enable = 0;
    st->mdct_sw = 0;
    reset_indices_dec( st );

    /* read the FT Header field from the bitstream */
    /* read the FT header */
    if ( fread( &header, sizeof(UWord8), 1, file ) != 1 )
    {
        if( ferror( file ) )
        {
            /* error during reading */
            fprintf(stderr, "\nError reading the bitstream !");
            exit(-1);
        }
        else
        {
            /* end of file reached */
            return 0;
        }
    }


    /* init local RXDTX flags */
    curr_ft_good_sp = 0;
    speech_lost = 0;
    speech_bad = 0;

    sid_update = 0;
    sid_upd_bad = 0;
    sti = -1;
    amrwb_sid_first = 0;  /* derived from sti  SID_FIRST indicator in AMRWB payload */
    no_data = 0;

    if( st->amrwb_rfc4867_flag != 0 )
    {
        /*   RFC 4867
        5.3 ....
        Each stored speech frame starts with a one-octet frame header with
        the following format:
        0 1 2 3 4 5 6 7
        +-+-+-+-+-+-+-+-+
        |P| FT    |Q|P|P|
        +-+-+-+-+-+-+-+-+
        The FT field and the Q bit are defined in the same way as in
        Section 4.3.2. The P bits are padding and MUST be set to 0, and MUST be ignored. */

        isAMRWB_IOmode   = 1;
        qbit             = (header>>2)&0x01 ;         /* b2 bit       (b7 is the F bit ) */
        st->bfi = !qbit;
        core_mode  = ((header>>3) & 0x0F);     /*  b6..b3      */
        total_brate = AMRWB_IOmode2rate[core_mode];   /* get the frame length from the header */
    }
    else
    {
        /*0 1 2 3 4 5 6 7   MS-bit ---> LS-bit
         +-+-+-+-+-+-+-+-+
         |H|F|E|x| brate |
         +-+-+-+-+-+-+-+-+
          where :
            "E|x|  brate "  is the 6 bit "FT" -field
             x is unused    if E=0, (should be 0 )
             x is the q-bit if E=1, q==1(good), Q==0(bad, maybe bit errors in payload )
             H,F  always   0 in RTP format.
        */
        isAMRWB_IOmode = (header & 0x20) > 0;   /* get EVS mode-from header */ /*    b2   */
        core_mode      = (header & 0x0F);        /* b4,b5,b6,b7 */

        if( isAMRWB_IOmode )
        {
            qbit = (header & 0x10) > 0;      /* get Q bit,    valid for IO rates */ /* b3 */
            total_brate = AMRWB_IOmode2rate[core_mode];
        }
        else
        {
            qbit = 1;  /* assume good q_bit for the unused EVS-mode bit,    complete ToC validity checked later */
            total_brate = PRIMARYmode2rate[ core_mode ];
        }
        st->bfi = !qbit;
    }

    /* set up RX-DTX-handler input */
    if(   core_mode == 14  )
    {
        /* SP_LOST */
        speech_lost=1;
    }
    if ( core_mode  == 15)
    {
        /*  NO_DATA unsent CNG frame OR  any frame marked or injected  as no_data  by e.g a signaling layer or dejitter buffer */
        no_data=1;
    }

    num_bits = (Word16)(total_brate/50);
    if( total_brate < 0 )
    {
        /* validate that total_brate (derived from RTP packet or a file header) is one of the defined bit rates  */
        fprintf(stderr, "\n  Error illegal total bit rate (= %d) in MIME ToC header \n",     total_brate );
        num_bits = -1;
    }

    /* Check correctness of ToC headers  */
    if( st->amrwb_rfc4867_flag == 0 )
    {
        /* EVS ToC header (FT field(b2-b7), H bit (b0),    F bit (b1)  ,  (EVS-modebit(b2)=0  unused(Qbit)(b3)==0)   */
        if ( (isAMRWB_IOmode == 0) &&  ((num_bits < 0)  ||  ((header & 0x80) > 0) || ((header & 0x40) > 0)  || (header & 0x30) != 0x00 )  )
        {
            /* incorrect FT header */
            fprintf(stderr, "\nError in EVS  FT ToC header(%02x) ! ",header);
            exit(-1);
        }
        else if( (isAMRWB_IOmode != 0) && ( (num_bits < 0) ||  ((header & 0x80) > 0) || ((header & 0x40) > 0) )  )  /* AMRWBIO */
        {
            /* incorrect IO FT header */
            fprintf(stderr, "\nError in EVS(AMRWBIO)  FT ToC header(%02x) ! ",header);
            exit(-1);
        }
    }
    else
    {
        /* legacy AMRWB ToC, is only using  Padding bits which MUST be ignored */
        if ( num_bits < 0  )
        {
            /* incorrect FT header */
            fprintf(stderr, "\nError in AMRWB RFC4867  Toc(FT)  header(%02x) !", header);
            exit(-1);
        }
    }



    /* read serial stream of indices from file to the local buffer */
    num_bytes_read = (Word16) fread( pFrame, sizeof(UWord8), (num_bits + 7)>>3, file );
    if( num_bytes_read != (num_bits + 7)>>3 )
    {
        fprintf(stderr, "\nError, invalid number of bytes read ! Exiting ! \n");
        exit(-1);
    }




    /* in case rew_flag is set, rewind the file and return */
    if ( rew_flag )
    {

        st->total_brate = total_brate;
        /* select MODE1 or MODE2 */
        if( st->bfi == 0 && speech_lost == 0 && no_data == 0 )
        {
            decoder_selectCodec( st, total_brate, unpack_bit(&pt_pFrame,&mask) ? G192_BIN1 : G192_BIN0);
        }
        return 1;
    }

    /* unpack speech data */
    bit_stream_ptr = st->bit_stream;
    for( k=0; k<num_bits; k++ )
    {
        if( isAMRWB_IOmode )
        {
            st->bit_stream[sort_ptr[core_mode][k]] = unpack_bit(&pt_pFrame,&mask);
            bit_stream_ptr++;
        }
        else
        {
            *bit_stream_ptr++ = unpack_bit(&pt_pFrame,&mask);
        }
    }

    /* unpack auxiliary bits */
    /* Note: these cmi bits are unpacked for demo purposes;  they are actually not needed  */
    if( isAMRWB_IOmode && total_brate == SID_1k75 )
    {
        sti = unpack_bit(&pt_pFrame,&mask);
        cmi  = unpack_bit(&pt_pFrame,&mask) << 3;
        cmi |= unpack_bit(&pt_pFrame,&mask) << 2;
        cmi |= unpack_bit(&pt_pFrame,&mask) << 1;
        cmi |= unpack_bit(&pt_pFrame,&mask);

        if( sti == 0 )
        {
            total_brate = 0;     /* signal received SID_FIRST as a good frame with no bits */
            for(k=0; k<35; k++)
            {
                st->bfi  |= st->bit_stream[k] ; /* partity check of 35 zeroes,  any single 1 gives BFI */
            }
        }

    }

    /*add two zero bytes for arithmetic coder flush*/
    for( k=0; k< 2*8; ++k )
    {
        *bit_stream_ptr++ = 0;
    }
    /* MIME RX_DTX handler */
    if( !rew_flag )
    {
        /* keep st->CNG , st_bfi and total_brate  updated  for proper synthesis in DTX and FER  */
        if( total_brate > SID_2k40 )
        {
            if( st->bfi == 0 )   /* so  far derived from q bit in AMRWB/AMRWBIO cases   */
            {
                curr_ft_good_sp = 1;
            }
        }

        /* handle q_bit and  lost_sp  clash ,  assume worst case  */
        if( speech_lost != 0 )  /*  overrides  a good q_bit */
        {
            curr_ft_good_sp = 0;
            st->bfi      = 1;     /* override  qbit */
        }

        /* now_bfi_fx has been set based on q_bit and ToC fields */

        /* SID_UPDATE check */
        if( total_brate == SID_1k75 || total_brate == SID_2k40 )
        {
            if( st->bfi == 0 )
            {
                /* typically from q bit  */
                sid_update = 1;
            }
            else
            {
                sid_upd_bad = 1;  /* may happen in saving from e.g. a CS-connection */
            }
        }

        if( isAMRWB_IOmode && total_brate == 0 && sti == 0 )
        {
            if( st->bfi )
            {
                sid_upd_bad = 1;          /*  corrupt sid_first, signaled as bad sid  */
            }
            else
            {
                amrwb_sid_first =  1;     /* 1-sti  */
            }
        }

        if( total_brate > SID_2k40 && st->bfi )  /* typically from q bit  */
        {
            speech_bad = 1;    /* initial assumption,   CNG synt state decides what to actually do */
        }
        /* all frame types decoded */

        /*    update CNG synthesis state */
        /*    Decoder can only  enter CNG-synthesis  for  CNG frame types (sid_upd,  sid_bad, sid_first) */
        if( st->CNG != 0 )
        {
            /* We were in CNG synthesis  */
            if( curr_ft_good_sp != 0  )
            {
                /* only a good speech frame makes decoder leave CNG synthesis */
                st->CNG = 0;
            }
        }
        else
        {
            /*   We were in SPEECH synthesis  */
            /*   only a received SID frame can make the decoder enter into CNG synthesis  */
            if( amrwb_sid_first || sid_update || sid_upd_bad )
            {
                st->CNG = 1;
            }
        }

        /* Now modify bfi flag for the  decoder's  SPEECH/CNG synthesis logic  */
        /*   in SPEECH synthesis, make sure to activate speech plc for a received no_data frame,
             bo_data frames may be injected by the network or by the dejitter buffer   */
        /*   modify bfi_flag to stay/move to the into the correct decoder PLC section  */
        if ( (st->CNG == 0)  &&  ( no_data != 0 )  )
        {
            /*  treat no_data received in speech synthesis as  SP_LOST frames, SPEECH PLC code will now become active */
            st->bfi = 1;
            /* total_brate= 0;    always zero for no_data */
        }

        /* in CNG  */
        /* handle bad speech frame(and bad sid frame) in the decoders CNG synthesis settings pair (total_brate, bfi)  */
        if( ( st->CNG != 0 && ( speech_bad || speech_lost || no_data ))  || /* SP_BAD or SPEECH_LOST)   --> stay in CNG */
                sid_upd_bad )                                                 /* SID_UPD_BAD               --> start/stay  CNG   */
        {
            st->bfi     = 0;   /* mark as good to not start speech PLC */
            total_brate = 0;   /* this zeroing needed  for  speech_bad, sid_bad frames */
        }
    }

    /*  now  bfi, total_brate are set by RX-DTX handler::
    bfi==0, total_brate!=0    cng or speech pending  bitrate
    bfi==0, total_brate==0    cng will continue or start(sid_first, sid_bad)
    bfi==1, total_brate!=0    speech plc
    bfi==1, total_brate==0 ,  speech plc
    */
    /*  handle available AMRWB/AMRWBIO MIME header ToC rate-info at startup   */
    if(  ( st->bfi != 0 && rew_flag == 0 && st->ini_frame == 0) &&    /*  ini_frame can not be used when rewflag is 1  */
            ( (st->amrwb_rfc4867_flag != 0)  || (st->amrwb_rfc4867_flag == 0 && isAMRWB_IOmode != 0  )) )  /*AMRWB ToC */
    {
        Word32 init_rate;

        init_rate = total_brate;       /* default , may have been modified from original ToC value */
        if (speech_lost !=0 || no_data != 0 )
        {
            init_rate =  ACELP_12k65;  /* make sure the decoder starts up in a selected AMRWB mode */
        }
        else if( speech_bad !=0  )
        {
            init_rate   =   AMRWB_IOmode2rate[core_mode];   /* read from from ToC */
        }
        st->total_brate  = init_rate; /* not updated on bfi as  decoderSelectCodec is not called below */
        st->core_brate   = init_rate;
    }
    if( st->bfi == 0 )
    {
        /* select MODE1 or MODE2 in  MIME */

        decoder_selectCodec( st, total_brate, *st->bit_stream ? G192_BIN1 : G192_BIN0);

        /* a change of the total bitrate should not be known to the decoder, if the received frame was truly lost */
        st->total_brate  = total_brate;
        mdct_switching_dec(st);
    }
    /* else{ bfi stay in past synthesis mode(SP,CNG) } */

    return 1;
}



/*-------------------------------------------------------------------*
* get_rfFrameType()
*
* Extract the rf frame type
*-------------------------------------------------------------------*/

static void get_rfFrameType(
    Decoder_State *st,              /* i  : decoder state structure */
    short *rf_frame_type    /* o  : RF frame type           */
)
{
    short num_bits;
    num_bits = st->total_brate/50;

    if( st->rf_flag == 1)
    {
        /* the last three bits in a packet is the RF frame type */
        *rf_frame_type = get_indice( st, num_bits - 3, 3 );
    }
    else
    {
        *rf_frame_type = 0;
    }

    return;
}

/*-------------------------------------------------------------------*
* get_rfFlag()
*
* Check if rf flag is present in the bitstream
*-------------------------------------------------------------------*/

static void get_rfFlag(
    Decoder_State *st,              /* i: decoder state structure    */
    short *rf_flag,         /* o  : check for the RF flag    */
    short *nBits,
    long *ind
)
{
    short start_idx, nBits_tmp;
    long ind_tmp;

    /* Init */
    *rf_flag = 0;

    /* check for rf_flag in the packet and extract the rf_frame_type and rf_fec_offset */
    if( st->total_brate == ACELP_13k20 && (st->bfi == FRAMEMODE_NORMAL || st->bfi == FRAMEMODE_FUTURE) )
    {
        /* find the section in the ACELP signalling table corresponding to bitrate */
        start_idx = 0;
        while ( acelp_sig_tbl[start_idx] != st->total_brate )
        {
            start_idx++;
            assert((start_idx < MAX_ACELP_SIG) && "ERROR: start_idx larger than acelp_sig_tbl[].\n");
        }

        /* skip the bitrate */
        start_idx += 1;

        /* retrieve the number of bits */
        nBits_tmp = (short) acelp_sig_tbl[start_idx++];

        /* retrieve the signalling indice */
        ind_tmp = acelp_sig_tbl[start_idx + get_indice( st, 0, nBits_tmp )];

        /* convert signalling indice into RF flag. */
        *rf_flag = (ind_tmp >> 7) & 0x1;

        if( ind )
        {
            *ind = ind_tmp;
        }

        if( nBits )
        {
            *nBits = nBits_tmp;
        }
    }

    return;
}

/*-------------------------------------------------------------------*
* get_rf_fec_offset()
*
* Extract the FEC offset
*-------------------------------------------------------------------*/

static void get_rf_fec_offset(
    Decoder_State *st,                      /* i  : decoder state structure       */
    short *rf_fec_offset            /* o  : RF fec offset                 */
)
{
    short num_bits, tmp;
    num_bits = st->total_brate/50;

    if( st->rf_flag == 1)
    {
        /* the two bits before the rf frame type contains the fec offset  */
        tmp = get_indice( st, num_bits - 5, 2 );

        if (tmp == 0)
        {
            *rf_fec_offset = 2;
        }
        else
        {
            *rf_fec_offset = 2*tmp + 1;
        }

    }
    else
    {
        *rf_fec_offset = 0;
    }

    return;
}

/*-------------------------------------------------------------------*
* get_rfTargetBits()
*
* Return the number of RF target bits
*-------------------------------------------------------------------*/

static void get_rfTargetBits(
    short rf_frame_type,           /* i  : RF frame type                 */
    short *rf_target_bits          /* o  : Number of RF target bits      */
)
{
    /* Number of RF bits for different RF coder types */

    switch (rf_frame_type)
    {
    case RF_NO_DATA:
        *rf_target_bits = 5;
        break;
    case RF_TCXFD:
        *rf_target_bits = 27;
        break;
    case RF_TCXTD1:
        *rf_target_bits = 16;
        break;
    case RF_TCXTD2:
        *rf_target_bits = 16;
        break;
    case RF_ALLPRED:
        /* Es_pred bits 3 bits, LTF: 1, pitch: 8,5,5,5, FCB: 0, gain: 7,0,7,0, Diff GFr: 4*/
        *rf_target_bits = 63;
        break;
    case RF_NOPRED:
        /* Es_pred bits 3 bits, LTF: 0, pitch: 0, FCB: 7,7,7,7, gain: 6,0,6,0, Diff GFr: 2*/
        *rf_target_bits = 66;
        break;
    case RF_GENPRED:
        /* Es_pred bits 3 bits, LTF: 1, pitch: 8,0,8,0, FCB: 6,7,5,5, gain: 5,0,5,0, Diff GFr: 0*/
        *rf_target_bits = 70;
        break;
    case RF_NELP:
        /* gain: 19, Diff GFr: 5 */
        *rf_target_bits =  45;
        break;
    }

    return;
}

/*-------------------------------------------------------------------*
* berCheck()
*
* Check for bit errors in channel aware signalling.
*-------------------------------------------------------------------*/

static void berCheck(
    Decoder_State *st,          /* i/o: decoder state structure     */
    short *coder_type           /* i/o: coder type                  */
)
{
    /* In case of RF flag = 1, and valid RF packet with primary and partial copy */
    if( st->bwidth == NB || st->bwidth == FB || *coder_type >= TRANSITION )
    {
        if( st->use_partial_copy == 1 )
        {
            st->use_partial_copy = 0;
        }

        st->bfi = 1;
        st->bwidth = st->last_bwidth;
        st->BER_detect = 1;
        *coder_type = GENERIC;
    }

    return;
}

/*-------------------------------------------------------------------*
 * getPartialCopyInfo()
 *
 * Check if the frame includes a partial copy for channel aware processing.
 *-------------------------------------------------------------------*/

void getPartialCopyInfo(
    Decoder_State *st,              /* i/o: decoder state structure       */
    short *coder_type,
    short *sharpFlag
)
{
    short nBits = 0;
    long ind = 0;

    /* check the rf flag in the packet */
    get_rfFlag( st, &(st->rf_flag), &nBits , &ind);

    /* get rf frame type info */
    get_rfFrameType( st, &(st->rf_frame_type) );

    /* Get the FEC offset info */
    get_rf_fec_offset( st, &(st->rf_fec_offset) );

    /* reset number of target bits in case of rate switching */
    st->rf_target_bits = 0;

    /* Get the number of bits used for RF*/
    if( st->rf_flag == 1 )
    {
        *coder_type = ind & 0x7;
        st->bwidth = (ind >> 3) & 0x7;
        *sharpFlag = (ind >> 6) & 0x1;
        st->codec_mode = MODE2;
        get_rfTargetBits( st->rf_frame_type, &(st->rf_target_bits) );

        if( st->bfi == FRAMEMODE_FUTURE )
        {
            st->use_partial_copy = 1;
            /* now set the frame mode to normal mode */
            if(st->rf_frame_type >= RF_TCXFD && st->rf_frame_type <= RF_TCXTD2)
            {
                st->bfi = 1;
                st->core = 1;
            }
            else
            {
                st->bfi = FRAMEMODE_NORMAL;
                st->core = 0;
            }
        }

        /* check for bit errors */
        berCheck( st, coder_type );

        get_next_indice_tmp(st, nBits);
    }

    return;
}

/*-------------------------------------------------------------------*
 * get_NextCoderType()
 *
 * Extract the coder type of next frame
 *-------------------------------------------------------------------*/

void get_NextCoderType(
    unsigned char *bitsteam,            /* i : bitstream            */
    short *next_coder_type              /* o : next coder type      */
)
{
    short k;
    short start_idx;
    char bit_stream[ACELP_13k20/50];
    long tmp;
    short nBits_tmp;


    for( k = 0; k < ACELP_13k20/50; k++ )
    {
        bit_stream[k] = (bitsteam[k / 8] >> (7 - (k % 8))) & 0x1;
    }
    start_idx = 0;
    while ( acelp_sig_tbl[start_idx] != ACELP_13k20 )
    {
        start_idx++;
        assert((start_idx < MAX_ACELP_SIG) && "ERROR: start_idx larger than acelp_sig_tbl[].\n");
    }

    /* skip the bitrate */
    start_idx += 1;

    tmp = 0;
    nBits_tmp = (short) acelp_sig_tbl[start_idx++];
    for (k = 0; k < nBits_tmp; k++)
    {
        tmp <<= 1;
        tmp += bit_stream[k];
    }

    /* retrieve the signalling indice */
    *next_coder_type = acelp_sig_tbl[start_idx + tmp] & 0x7;

    return;
}

/*-------------------------------------------------------------------*
 * read_indices_from_djb()
 *
 * Read indices from the de-jitter buffer payload (works also for AMR-WB IO mode)
 *-------------------------------------------------------------------*/

void read_indices_from_djb(
    Decoder_State *st,                      /* i/o: decoder state structure       */
    unsigned char *pt_stream,               /* i  : bitstream file                */
    int num_bits,                 /* i  : input frame length in bits    */
    short partialframe,             /* i  : partial frame information     */
    short next_coder_type           /* i  : next coder type information   */
)
{
    int k;
    unsigned short *bit_stream_ptr;
    long total_brate;
    short bit0;

    st->BER_detect = 0;
    bit0 = 0;
    /* There is no FRAME_NO_DATA or BAD frame indicator in RTP, frames are just missing.
     * In case of comfort noise handle missing frame as FRAME_NO_DATA, otherwise use PLC. */
    if(num_bits != 0)
    {
        st->bfi = 0;
        bit0 = ((*pt_stream & 0x80) != 0) ? G192_BIN1 : G192_BIN0;
    }
    else if(st->total_brate == SID_1k75 || st->total_brate == SID_2k40 ||
            st->total_brate == FRAME_NO_DATA)
    {
        st->bfi = 0;
    }
    else
    {
        st->bfi = 1;
    }
    if( partialframe || st->prev_use_partial_copy)
    {
        st->next_coder_type = next_coder_type;
    }
    else
    {
        st->next_coder_type = INACTIVE;
    }

    st->mdct_sw_enable = 0;
    st->mdct_sw = 0;
    reset_indices_dec( st );
    total_brate = num_bits * 50;

    if(partialframe == 1)
    {
        st->bfi = 2;
    }
    if ( st->bfi != 1 )
    {
        /* select Mode 1 or Mode 2 */
        decoder_selectCodec( st, total_brate, bit0 );

        /* convert bitstream from compact bytes to short values and store it in decoder state */
        bit_stream_ptr = st->bit_stream;
        for( k = 0; k < num_bits; k++ )
        {
            *bit_stream_ptr++ = (pt_stream[k / 8] >> (7 - (k % 8))) & 0x1;
        }
        /* add two zero bytes for arithmetic coder flush */
        for( k=0; k < 8*2; ++k )
        {
            *bit_stream_ptr++ = 0;
        }

        /* a change of the total bitrate should not be known to the decoder, if the received frame was lost */
        st->total_brate = total_brate;

        mdct_switching_dec(st);
    }

    return;
}

/*-------------------------------------------------------------------*
 * get_indice_preview()
 *
 * Indices preview to parse for the presence of partial copy
 *-------------------------------------------------------------------*/
static unsigned short get_indice_preview(
    unsigned char *bitstream,
    int bitstreamSize,
    short pos,
    short nb_bits
)
{
    unsigned short value;
    int i;
    unsigned short bitstreamShort[MAX_BITS_PER_FRAME+16];
    unsigned short *bitstreamShortPtr;

    /* convert bitstream from compact bytes to short values */
    bitstreamShortPtr = bitstreamShort;
    for( i = 0; i < bitstreamSize; i++ )
    {
        *bitstreamShortPtr++ = (bitstream[i / 8] >> (7 - (i % 8))) & 0x1;
    }

    assert(nb_bits <= 16);
    value = 0;
    for (i = 0; i < nb_bits; i++)
    {
        value <<= 1;
        value += bitstreamShort[pos+i];
    }
    return value;
}

/*-------------------------------------------------------------------*
 * evs_dec_previewFrame()
 *
 * Signalling index preview
 *-------------------------------------------------------------------*/
void evs_dec_previewFrame(
    unsigned char *bitstream,
    int bitstreamSize,
    short *partialCopyFrameType,
    short *partialCopyOffset
)
{
    long total_brate;
    short start_idx, nBits;
    long ind;
    short rf_flag;

    rf_flag = 0;
    *partialCopyFrameType = 0;
    *partialCopyOffset = 0;
    total_brate = bitstreamSize * 50;

    if( total_brate == ACELP_13k20 )
    {
        /* find the section in the ACELP signalling table corresponding to bitrate */
        start_idx = 0;
        while ( acelp_sig_tbl[start_idx] != total_brate )
        {
            start_idx++;
            assert((start_idx < MAX_ACELP_SIG) && "ERROR: start_idx larger than acelp_sig_tbl[].\n");
        }

        /* skip the bitrate */
        start_idx += 1;
        /* retrieve the number of bits */
        nBits = (short) acelp_sig_tbl[start_idx++];

        /* retrieve the signalling indice */
        ind = acelp_sig_tbl[start_idx + get_indice_preview( bitstream, bitstreamSize, 0, nBits )];

        /* convert signalling indice into RF flag. */
        rf_flag = (ind >> 7) & 0x1;
        if(rf_flag != 0)
        {
            /* read the fec offset at which the partial copy is received */
            ind = get_indice_preview( bitstream, bitstreamSize, (bitstreamSize-5), 2 );
            if(ind== 0) *partialCopyOffset = 2;
            else if(ind == 1) *partialCopyOffset = 3;
            else if(ind == 2) *partialCopyOffset = 5;
            else if(ind == 3) *partialCopyOffset = 7;
            /* the last three bits in a packet is the RF frame type */
            *partialCopyFrameType = get_indice_preview( bitstream, bitstreamSize, bitstreamSize - 3, 3 );
        }
    }

    return;
}

