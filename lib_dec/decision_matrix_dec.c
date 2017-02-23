/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <assert.h>
#include "options.h"
#include "stat_dec.h"
#include "rom_com.h"
#include "prot.h"

/*-----------------------------------------------------------------*
 * decision_matrix_dec()
 *
 * ACELP/HQ core selection
 * Read ACELP signalling bits from the bitstream
 * Set extension layers
 *-----------------------------------------------------------------*/

void decision_matrix_dec(
    Decoder_State *st,                /* i/o: decoder state structure                   */
    short *coder_type,        /* o  : coder type                                */
    short *sharpFlag,         /* o  : formant sharpening flag                   */
    short *hq_core_type,      /* o  : HQ core type                              */
    short *core_switching_flag/* o  : ACELP->HQ switching frame flag            */
)
{
    short start_idx;
    short ppp_nelp_mode;
    long ind;
    short nBits;
    short tmp;
    st->core = -1;
    st->core_brate = 0;
    st->extl = -1;
    st->extl_brate = 0;
    st->ppp_mode_dec = 0;
    st->nelp_mode_dec = 0;
    st->igf = 0;

    if( st->total_brate > ACELP_8k00 )
    {
        st->vbr_hw_BWE_disable_dec = 0;
    }

    if( st->mdct_sw == MODE2 )
    {
        st->core = HQ_CORE;
    }
    else
    {
        if ( st->total_brate == FRAME_NO_DATA || st->total_brate == SID_2k40 )
        {
            st->core = ACELP_CORE;
            st->core_brate = st->total_brate;

            if( st->total_brate != FRAME_NO_DATA )
            {
                st->cng_type = get_next_indice( st, 1 );

                if( st->cng_type == LP_CNG )
                {
                    st->L_frame = L_FRAME;

                    if( get_next_indice( st, 1 ) == 1 )
                    {
                        st->L_frame = L_FRAME16k;
                    }
                }
                else
                {
                    st->bwidth = get_next_indice(st, 2);

                    if( get_next_indice(st, 1) == 0 )
                    {
                        st->L_frame = L_FRAME;
                    }
                    else
                    {
                        st->L_frame = L_FRAME16k;
                    }
                }
            }

            if( st->output_Fs >= 32000 && st->bwidth >= SWB )
            {
                st->extl = SWB_CNG;
            }
            if ( st->total_brate == FRAME_NO_DATA && st->prev_bfi && !st->bfi && st->L_frame > L_FRAME16k)
            {
                st->L_frame = st->last_CNG_L_frame;
            }

            return;
        }

        /* SC-VBR */
        else if ( st->total_brate == PPP_NELP_2k80 )
        {
            st->core = ACELP_CORE;
            st->core_brate = PPP_NELP_2k80;
            st->L_frame = L_FRAME;
            st->fscale = sr2fscale((int)INT_FS_12k8);

            if( st->ini_frame == 0 )
            {
                /* avoid switching of internal ACELP Fs in the very first frame */
                st->last_L_frame = st->L_frame;
                st->last_core = st->core;
                st->last_core_brate = st->core_brate;
                st->last_extl = st->extl;
            }

            st->vbr_hw_BWE_disable_dec = 1;

            get_next_indice( st, 1 );

            ppp_nelp_mode = (short) get_next_indice( st, 2 );

            /* 0 - PPP_NB, 1 - PPP_WB, 2 - NELP_NB, 3 - NELP_WB */
            if ( ppp_nelp_mode == 0 )
            {
                st->ppp_mode_dec = 1;
                *coder_type = VOICED;
                st->bwidth = NB;
            }
            else if ( ppp_nelp_mode == 1 )
            {
                st->ppp_mode_dec = 1;
                *coder_type = VOICED;
                st->bwidth = WB;
            }
            else if ( ppp_nelp_mode == 2 )
            {
                st->nelp_mode_dec = 1;
                *coder_type = UNVOICED;
                st->bwidth = NB;
            }
            else if ( ppp_nelp_mode == 3 )
            {
                st->nelp_mode_dec = 1;
                *coder_type = UNVOICED;
                st->bwidth = WB;
            }

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
            tmp = (short) get_next_indice( st, 1 );

            if( tmp == 0 )
            {
                st->core = ACELP_CORE;
            }
            else
            {
                st->core = HQ_CORE;
            }
        }
    }

    /*-----------------------------------------------------------------*
     * Read ACELP signalling bits from the bitstream
     *-----------------------------------------------------------------*/

    if( st->core == ACELP_CORE )
    {
        /* find the section in the ACELP signalling table corresponding to bitrate */
        start_idx = 0;
        while ( acelp_sig_tbl[start_idx] != st->total_brate )
        {
            start_idx++;
            if( start_idx >= MAX_ACELP_SIG )
            {
                st->BER_detect = 1;
                start_idx--;
                break;
            }
        }

        /* skip the bitrate */
        start_idx += 1;

        /* retrieve the number of bits */
        nBits = (short) acelp_sig_tbl[start_idx++];

        start_idx += get_next_indice( st, nBits );
        if( start_idx >= MAX_ACELP_SIG )
        {
            ind = 0;
            st->BER_detect = 1;
        }
        else
        {
            /* retrieve the signalling indice */
            ind = acelp_sig_tbl[start_idx];

            /* convert signalling indice into signalling information */
            *coder_type = ind & 0x7;
            if ( *coder_type == LR_MDCT )
            {
                st->core = HQ_CORE;
                st->bwidth = (ind >> 3) & 0x7;
            }
            else
            {
                st->bwidth = (ind >> 3) & 0x7;
                *sharpFlag = (ind >> 6) & 0x1;
            }
        }

        /* detect corrupted signalling (due to bit errors) */
        if( ( st->BER_detect ) ||
                ( ind >= 1<<7 ) ||
                ( st->total_brate <= ACELP_13k20 && st->bwidth == FB ) ||
                ( st->total_brate >= ACELP_32k && st->bwidth == NB ) ||
                ( st->total_brate >= ACELP_32k && !(*coder_type == GENERIC || *coder_type == TRANSITION || *coder_type == INACTIVE ) ) ||
                ( st->total_brate < ACELP_13k20 && st->bwidth != NB && *coder_type == LR_MDCT ) ||
                ( st->total_brate >= ACELP_13k20 && *coder_type == UNVOICED ) ||
                ( st->total_brate >= ACELP_13k20 && *coder_type == AUDIO && st->bwidth == NB )
          )
        {
            st->BER_detect = 0;
            st->bfi = 1;
            if( st->ini_frame == 0 )
            {
                st->core = ACELP_CORE;
                st->L_frame = L_FRAME;
                st->last_core = st->core;
                st->last_core_brate = st->core_brate;
            }
            else if( st->last_total_brate == -1 ) /* can happen in case of BER when no good frame was received before */
            {
                *coder_type = st->last_coder_type;
                st->bwidth = st->last_bwidth;
                st->total_brate = st->last_total_brate_ber;
                if( st->last_core == AMR_WB_CORE )
                {
                    st->core = ACELP_CORE;
                    st->codec_mode = MODE1;
                }
                else if( st->last_core_bfi == TCX_20_CORE || st->last_core_bfi == TCX_10_CORE )
                {
                    st->core = st->last_core_bfi;
                    st->codec_mode = MODE2;
                }
                else
                {
                    st->core = st->last_core;
                    st->codec_mode = MODE1;
                }
                st->core_brate = st->last_core_brate;
                st->extl = st->last_extl;
                st->extl_brate = st->total_brate - st->core_brate;
            }
            else
            {
                *coder_type = st->last_coder_type;
                st->bwidth = st->last_bwidth;
                st->total_brate = st->last_total_brate;
                if( st->last_core == AMR_WB_CORE )
                {
                    st->core = ACELP_CORE;
                    st->codec_mode = MODE1;
                }
                else if( st->last_core == TCX_20_CORE || st->last_core == TCX_10_CORE )
                {
                    st->core = st->last_core;
                    st->codec_mode = MODE2;
                }
                else
                {
                    st->core = st->last_core;
                    st->codec_mode = MODE1;
                }
                st->core_brate = st->last_core_brate;
                st->extl = st->last_extl;
                st->extl_brate = st->total_brate - st->core_brate;
            }
            return;
        }
    }

    /*-----------------------------------------------------------------*
     * Set extension layers
     *-----------------------------------------------------------------*/

    if ( st->core == ACELP_CORE && st->bwidth == WB && st->total_brate < ACELP_9k60 )
    {
        if ( st->vbr_hw_BWE_disable_dec == 0 )
        {
            st->extl = WB_BWE;
        }
    }
    else if ( st->core == ACELP_CORE && st->bwidth == WB && st->total_brate >= ACELP_9k60 && st->total_brate <= ACELP_16k40 )
    {
        /* read the WB TBE/BWE selection bit */
        if ( get_next_indice( st, 1 ) )
        {
            st->extl = WB_BWE;
            st->extl_brate = WB_BWE_0k35;
        }
        else
        {
            st->extl = WB_TBE;
            st->extl_brate = WB_TBE_1k05;
        }
    }
    else if ( st->core == ACELP_CORE && (st->bwidth == SWB || st->bwidth == FB) && st->total_brate >= ACELP_13k20 )
    {
        if (st->total_brate >= ACELP_48k)
        {
            st->extl = SWB_BWE_HIGHRATE;
            if( st->bwidth == FB )
            {
                st->extl = FB_BWE_HIGHRATE;
            }

            st->extl_brate = SWB_BWE_16k;
        }

        /* read the SWB TBE/BWE selection bit */
        else
        {
            tmp = get_next_indice( st, 1 );

            if( tmp )
            {
                st->extl = SWB_BWE;
                st->extl_brate = SWB_BWE_1k6;
            }
            else
            {
                st->extl = SWB_TBE;
                st->extl_brate = SWB_TBE_1k6;
                if( st->total_brate >= ACELP_24k40 )
                {
                    st->extl_brate = SWB_TBE_2k8;
                }
            }
        }

        /* set FB TBE and FB BWE extension layers */
        if ( st->bwidth == FB && st->total_brate >= ACELP_24k40 )
        {
            if ( st->extl == SWB_BWE )
            {
                st->extl = FB_BWE;
                st->extl_brate = FB_BWE_1k8;
            }
            else if ( st->extl == SWB_TBE )
            {
                st->extl = FB_TBE;
                st->extl_brate = FB_TBE_3k0;
            }
        }
    }

    /* set core bitrate */
    st->core_brate = st->total_brate - st->extl_brate;

    /*-----------------------------------------------------------------*
     * Read HQ signalling bits from the bitstream
     * Set HQ core type
     *-----------------------------------------------------------------*/


    if ( st->core == HQ_CORE )
    {
        if( st->mdct_sw != MODE2 )
        {
            /* skip the HQ/TCX core switching flag */
            get_next_indice_tmp( st, 1 );
        }

        /* read ACELP->HQ core switching flag */
        *core_switching_flag = (short)get_next_indice( st, 1 );

        if( *core_switching_flag == 1 )
        {
            st->last_L_frame_ori = st->last_L_frame;

            /* read ACELP L_frame info */
            if( get_next_indice( st, 1 ) == 0 )
            {
                st->last_L_frame = L_FRAME;
            }
            else
            {
                st->last_L_frame = L_FRAME16k;
            }
        }

        if( st->mdct_sw != MODE2 )
        {
            /* read/set band-width (needed for different I/O sampling rate support) */
            if( st->total_brate > ACELP_16k40 )
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

        /* detect bit errors in signalling */
        if( ( st->total_brate >= ACELP_24k40 && st->bwidth == NB ) ||
                ( st->core == HQ_CORE && st->total_brate <= LRMDCT_CROSSOVER_POINT && st->bwidth == FB)
          )
        {
            st->bfi = 1;

            st->core_brate = st->total_brate;
            st->extl = -1;
            st->extl_brate = 0;
            if( st->last_core == AMR_WB_CORE )
            {
                st->core = ACELP_CORE;
                st->L_frame = L_FRAME;
                st->codec_mode = MODE1;
                st->last_L_frame = L_FRAME;

                if( st->total_brate >= ACELP_16k40 )
                {
                    st->total_brate = ACELP_13k20;
                    st->core_brate = st->total_brate;
                }
            }
        }

        /* set HQ core type */
        *hq_core_type = NORMAL_HQ_CORE;
        if( (st->bwidth == SWB || st->bwidth == WB) && st->total_brate <= LRMDCT_CROSSOVER_POINT )
        {
            *hq_core_type = LOW_RATE_HQ_CORE;
        }
        else if( st->bwidth == NB )
        {
            *hq_core_type = LOW_RATE_HQ_CORE;
        }
    }

    /*-----------------------------------------------------------------*
     * Set ACELP frame length
     *-----------------------------------------------------------------*/

    if( st->core_brate == FRAME_NO_DATA )
    {
        /* prevent "L_frame" changes in CNG segments */
        st->L_frame = st->last_L_frame;
    }
    else if ( st->core_brate == SID_2k40 && st->bwidth == WB && st->first_CNG && st->act_cnt2 < MIN_ACT_CNG_UPD )
    {
        /* prevent "L_frame" changes in SID frame after short segment of active frames */
        st->L_frame = st->last_CNG_L_frame;
    }
    else if ( ( st->core_brate == SID_2k40 && st->total_brate >= ACELP_9k60 && st->bwidth == WB ) ||
              ( st->total_brate > ACELP_24k40 && st->total_brate < HQ_96k ) || ( st->total_brate == ACELP_24k40 && st->bwidth >= WB ) )
    {
        st->L_frame = L_FRAME16k;
    }
    else
    {
        st->L_frame = L_FRAME;
    }

    if ( st->L_frame == L_FRAME16k )
    {
        st->nb_subfr = NB_SUBFR16k;
    }
    else
    {
        st->nb_subfr = NB_SUBFR;
    }

    if( st->output_Fs == 8000 )
    {
        st->extl = -1;
    }
    else if( st->output_Fs == 16000 && st->L_frame == L_FRAME16k )
    {
        st->extl = -1;
        st->extl_brate = 0;
    }

    if ( st->ini_frame == 0 )
    {
        /* avoid switching of internal ACELP Fs in the very first frame */
        st->last_L_frame = st->L_frame;
        st->last_core = st->core;
        st->last_core_brate = st->core_brate;
        st->last_extl = st->extl;
    }

    return;
}
