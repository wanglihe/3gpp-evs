/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <assert.h>
#include "options.h"
#include "prot.h"
#include "stat_enc.h"
#include "stat_dec.h"
#include "rom_com.h"


/*-----------------------------------------------------------------*
 * decision_matrix_enc()
 *
 * Select operating point (combination of technologies) based on input signal properties and command-line parameters:
 *
 *             7.20        8.00        9.60        13.20        16.40         24.40            32               48               64               96      128
 *  Mode       1           1           2           1            2             2                1                2                1                2       2
 * ----------------------------------------------------------------------------------------------------------------------------------------------------------------
 *  NB
 *  speech     ACELP@12k8  ACELP@12k8  ACELP@12k8  ACELP@12k8
 *  audio      LR MDCT     LR MDCT     TCX         LR MDCT
 *  inactive   GSC@12k8    GSC@12k8    TCX         GSC@12k8
 * ----------------------------------------------------------------------------------------------------------------------------------------------------------------
 *  WB
 *  speech     ACELP@12k8  ACELP@12k8  ACELP@12k8  ACELP@12k8   ACELP@16k     ACELP@16k        ACELP@16k        TCX              ACELP@16k        TCX     TCX
 *             +0b WB BWE  +0b WB BWE  +TD WB BWE  +TD WB BWE
 *  audio      GSC@12k8    GSC@12k8    TCX         LR MDCT      TCX           TCX              HQ               TCX              HQ               TCX     TCX
 *             +0b WB BWE  +0b WB BWE  +IGF
 *  inactive   GSC@12k8    GSC@12k8    TCX         GSC@12k8     TCX           TCX              AVQ@16k          TCX              AVQ@16k          TCX     TCX
 *             +0b WB BWE  +0b WB BWE  +IGF        +FD WB BWE
 * ----------------------------------------------------------------------------------------------------------------------------------------------------------------
 *  SWB
 *  speech                                         ACELP@12k8   ACELP@16k     ACELP@16k        ACELP@16k        TCX              ACELP@16k        TCX     TCX
 *                                                 +TD SWB BWE  +TD SWB BWE   +TD SWB BWE      +TD SWB BWE      +IGF             +HR SWB BWE
 *  audio                                          LR MDCT/GSC  TCX           TCX              HQ               TCX              HQ               TCX     TCX
 *                                                 +FD SWB BWE  +IGF          +IGF             +FD SWB BWE      +IGF
 *  inactive                                       GSC@12k8     TCX           TCX              AVQ@16k          TCX              AVQ@16k          TCX     TCX
 *                                                 +FD SWB BWE  +IGF          +IGF             +FD SWB BWE      +IGF             +HR SWB BWE
 * ----------------------------------------------------------------------------------------------------------------------------------------------------------------
 *  FB
 *  speech                                                      ACELP@16k     ACELP@16k        ACELP@16k        TCX              ACELP@16k        TCX     TCX
 *                                                              +TD FB BWE    +TD FB BWE       +TD FB BWE       +IGF             +HR FB BWE
 *  audio                                                       TCX           TCX              HQ               TCX              HQ               TCX     TCX
 *                                                              +IGF          +IGF             +FD FB BWE       +IGF
 *  inactive                                                    TCX           TCX              AVQ@16k          TCX              AVQ@16k          TCX     TCX
 *                                                              +IGF          +IGF             +FD FB BWE       +IGF             +HR FB BWE
 * -----------------------------------------------------------------------------------------------------------------------------------------------------------------
 *
 * Note: the GSC technology is part of the ACELP core as AUDIO coder_type (it is used also at 13.2 for SWB unvoiced noisy speech)
 * Note2: FB processing is optional and is activated via "-band FB" option on the encoder command line
 * Note3: NB (0-4kHz), WB (0-8kHz), SWB (0-16kHz), FB (0-20kHz)
 *
 * Signalling of modes (x marks a mode that must be signalled in the bitstream)
 *
 *                     7.20           8.00           9.6            13.2           16.4           24.4           32             48             64
 *                     NB WB SWB FB   NB WB SWB FB   NB WB SWB FB   NB WB SWB FB   NB WB SWB FB   NB WB SWB FB   NB WB SWB FB   NB WB SWB FB   NB WB SWB FB
 * GC, 12k8            x  x           x  x           x  x           x  x   x       x  x   x       x
 * UC, 12k8            x  x           x  x           x  x
 * VC, 12k8            x  x           x  x           x  x           x  x   x       x  x   x       x
 * TC, 12k8            x  x           x  x           x  x           x  x   x       x  x   x       x
 * GC, 16k                                                                                           x   x   x      x   x   x      x   x   x      x   x   x
 * TC, 16k                                                                                           x   x   x      x   x   x      x   x   x      x   x   x
 * AC(GSC)             x  x           x  x           x  x           x  x   x       x  x   x       x
 * IC                  x  x           x  x           x  x           x  x   x       x  x   x       x  x   x   x      x   x   x      x   x   x      x   x   x
 *
 * GC, 12k8, FS        x  x           x  x           x  x           x  x   x       x  x   x       x
 * GC, 16k, FS                                                                                       x   x   x      x   x   x      x   x   x          x   x
 * VC, 12k8, FS                                                     x  x   x       x  x   x       x
 * TC, 12k8, FS                                                                                   x
 * TC, 16k, FS                                                                                       x   x   x      x   x   x      x   x   x          x   x
 *
 * LR MDCT             x              x              x              x  x   x   x   x  x   x   x
 *
 *-----------------------------------------------------------------*/

void decision_matrix_enc(
    Encoder_State *st,                /* i  : encoder state structure                   */
    const short sp_aud_decision1,   /* i  : 1st stage speech/music classification     */
    const short sp_aud_decision2,   /* i  : 2nd stage speech/music classification     */
    const short coder_type,         /* i  : coder type                                */
    const short vad_flag,
    short *hq_core_type       /* o  : HQ core type                              */
)
{
    /* initialization */
    st->core = -1;
    st->extl = -1;
    st->extl_brate = 0;
    *hq_core_type = -1;
    st->igf = 0;

    /* SID and FRAME_NO_DATA frames */
    if( st->Opt_DTX_ON && (st->core_brate == SID_2k40 || st->core_brate == FRAME_NO_DATA ) )
    {
        st->core = ACELP_CORE;

        if( st->input_Fs >= 32000 && st->bwidth >= SWB )
        {
            st->extl = SWB_CNG;
        }

        st->rf_mode = 0;

        return;
    }

    st->core_brate = 0;

    /* SC-VBR */
    if ( st->Opt_SC_VBR )
    {
        /* SC-VBR */
        st->core = ACELP_CORE;
        st->core_brate = ACELP_7k20;
        st->total_brate = ACELP_7k20;

        if ( st->ppp_mode == 1 )
        {
            /* PPP mode */
            st->core_brate = PPP_NELP_2k80;
        }
        else if( ((coder_type == UNVOICED || coder_type == TRANSITION) && !sp_aud_decision1) || st->bwidth != NB )
        {

            if( coder_type == UNVOICED && vad_flag == 1 && ( ( st->last_bwidth >= SWB  && st->last_Opt_SC_VBR ) || st->last_bwidth < SWB )
                    && (st->last_core!=HQ_CORE || st->bwidth != NB) )
            {
                /* NELP mode */
                st->nelp_mode = 1;
                st->core_brate = PPP_NELP_2k80;
            }
            else if ( coder_type == TRANSITION || ( coder_type == UNVOICED && st->nelp_mode != 1 ) ||
                      ( ( coder_type == AUDIO || coder_type == INACTIVE ) && st->bwidth != NB ) )
            {
                /* silence portions */
                st->core_brate = ACELP_8k00;
                st->total_brate = ACELP_8k00;
            }
        }

        return;
    }

    /*---------------------------------------------------------------------*
     * NB
     *---------------------------------------------------------------------*/

    else if ( st->bwidth == NB )
    {
        st->core = ACELP_CORE;

        if( st->total_brate >= HQCORE_NB_MIN_RATE && sp_aud_decision1 == 1 )
        {
            st->core = HQ_CORE;
        }
    }

    /*---------------------------------------------------------------------*
     * WB
     *---------------------------------------------------------------------*/

    else if ( st->bwidth == WB )
    {
        st->core = ACELP_CORE;

        if( ( st->total_brate >= HQCORE_WB_MIN_RATE && sp_aud_decision1 == 1 ) || st->total_brate >= HQ_96k )
        {
            st->core = HQ_CORE;
        }
        else
        {
            if ( st->bwidth == WB && st->total_brate < ACELP_9k60 )
            {
                st->extl = WB_BWE;
            }
            else if ( st->bwidth == WB && st->total_brate >= ACELP_9k60 && st->total_brate <= ACELP_16k40 )
            {
                /* Note: WB BWE is used exceptionally at 13.2 kbps if GSC is selected instead of LR-MDCT */
                if ( sp_aud_decision1 == 1 || coder_type == INACTIVE || ( sp_aud_decision1 == 0 && sp_aud_decision2 == 1 ) )
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
        }
    }

    /*---------------------------------------------------------------------*
     * SWB and FB
     *---------------------------------------------------------------------*/

    else if ( st->bwidth == SWB || st->bwidth == FB )
    {
        if( ( st->total_brate >= HQCORE_SWB_MIN_RATE && sp_aud_decision1 == 1 ) || st->total_brate >= HQ_96k )
        {
            st->core = HQ_CORE;
        }
        else
        {
            st->core = ACELP_CORE;

            if ( st->total_brate >= ACELP_13k20 && st->total_brate < ACELP_48k )
            {
                /* Note: SWB BWE is not used in case of GSC noisy speech */
                /* Note: SWB BWE is used exceptionally at 13.2 kbps if GSC is selected instead of LR-MDCT */
                if ( ( sp_aud_decision1 == 1 || coder_type == INACTIVE || ( sp_aud_decision1 == 0 && sp_aud_decision2 == 1 ) ) && !st->GSC_noisy_speech )
                {
                    st->extl = SWB_BWE;
                    st->extl_brate = SWB_BWE_1k6;

                    if ( st->bwidth == FB && st->total_brate >= ACELP_24k40 )
                    {
                        st->extl = FB_BWE;
                        st->extl_brate = FB_BWE_1k8;
                    }
                }
                else
                {
                    st->extl = SWB_TBE;
                    st->extl_brate = SWB_TBE_1k6;

                    if( st->total_brate >= ACELP_24k40 )
                    {
                        st->extl_brate = SWB_TBE_2k8;
                    }

                    if ( st->bwidth == FB && st->total_brate >= ACELP_24k40 )
                    {
                        st->extl = FB_TBE;
                        st->extl_brate = FB_TBE_3k0;
                    }
                }
            }
            else if ( st->total_brate >= ACELP_48k )
            {
                st->extl = SWB_BWE_HIGHRATE;
                st->extl_brate = SWB_BWE_16k;

                if( st->bwidth == FB )
                {
                    st->extl = FB_BWE_HIGHRATE;
                }
            }
        }
    }

    /*-----------------------------------------------------------------*
     * Set HQ core type
     *-----------------------------------------------------------------*/

    if( st->core == HQ_CORE )
    {
        *hq_core_type = NORMAL_HQ_CORE;

        if( (st->bwidth == SWB || st->bwidth == WB) && st->total_brate <= LRMDCT_CROSSOVER_POINT )
        {
            /* note that FB (bit-rate >= 24400bps) is always coded with NORMAL_HQ_CORE */
            *hq_core_type = LOW_RATE_HQ_CORE;
        }
        else if( st->bwidth == NB )
        {
            *hq_core_type = LOW_RATE_HQ_CORE;
        }
    }

    /* set core bitrate */
    st->core_brate = st->total_brate - st->extl_brate;

    if ( st->ini_frame == 0 )
    {
        /* avoid switching in the very first frame */
        st->last_core = st->core;
        st->last_core_brate = st->core_brate;
        st->last_extl = st->extl;
    }

    return;
}

/*---------------------------------------------------------------------*
 * signalling_mode1_tcx20_enc()
 *
 * write MODE1 TCX20 signalling information into the bit-stream
 *---------------------------------------------------------------------*/

short signalling_mode1_tcx20_enc(
    Encoder_State *st,                /* i/o: encoder state structure   */
    const short push                /* i  : flag to push indice       */
)
{
    short num_bits;
    short nBits, idx, start_idx;

    assert(st->core == TCX_20_CORE);

    num_bits = 0;

    /* Use ACELP signaling for LR MDCT */
    if ( st->total_brate <= ACELP_16k40 )
    {
        /* find the section in the ACELP signalling table corresponding to bitrate */
        idx = 0;
        while ( acelp_sig_tbl[idx] != st->total_brate )
        {
            idx++;
        }

        /* retrieve the number of bits for signalling */
        nBits = (short) acelp_sig_tbl[++idx];

        /* retrieve the signalling index */
        start_idx = ++idx;
        while( acelp_sig_tbl[idx] != SIG2IND(LR_MDCT, st->bwidth, 0, 0) )
        {
            idx++;
        }

        num_bits += nBits;
        if( push )
        {
            push_indice( st, IND_ACELP_SIGNALLING, idx - start_idx, nBits );
        }

        /* HQ/TCX core switching flag */
        ++num_bits;
        if( push )
        {
            push_indice( st, IND_MDCT_CORE, 1, 1 );
        }
    }
    else
    {
        if( st->core_brate <= ACELP_64k )
        {
            /* write ACELP/HQ core indication flag */
            ++num_bits;
            if( push )
            {
                push_indice( st, IND_CORE, 1, 1 );
            }
        }

        /* HQ/TCX core switching flag */
        ++num_bits;
        if( push )
        {
            push_indice( st, IND_MDCT_CORE, 1, 1 );
        }

        num_bits += 2;
        if( push )
        {
            /* write band-width (needed for different I/O sampling rate support) */
            if( st->bwidth == NB )
            {
                push_indice( st, IND_HQ_BWIDTH, 0, 2 );
            }
            else if( st->bwidth == WB )
            {
                push_indice( st, IND_HQ_BWIDTH, 1, 2 );
            }
            else if( st->bwidth == SWB )
            {
                push_indice( st, IND_HQ_BWIDTH, 2, 2 );
            }
            else  /* st->bwidth == FB */
            {
                push_indice( st, IND_HQ_BWIDTH, 3, 2 );
            }
        }
    }

    return num_bits;
}


/*---------------------------------------------------------------------*
 * signalling_enc()
 *
 * write signalling information into the bit-stream
 *---------------------------------------------------------------------*/

void signalling_enc(
    Encoder_State *st,                /* i  : encoder state structure   */
    const short coder_type,         /* i  : coder type                */
    const short sharpFlag           /* i  : formant sharpening flag   */
)
{
    short nBits, idx, start_idx;

    if( st->mdct_sw == MODE2 )
    {

        assert(!st->tcxonly);
        assert(st->core == HQ_CORE);

        push_next_indice(st, 1, 1); /* TCX */
        push_next_indice(st, 1, 1); /* HQ_CORE */

        /* write ACELP->HQ core switching flag */
        if( st->last_core == ACELP_CORE || st->last_core == AMR_WB_CORE )
        {
            push_indice( st, IND_HQ_SWITCHING_FLG, 1, 1 );

            /* write ACELP L_frame info */
            if( st->last_L_frame == L_FRAME )
            {
                push_indice( st, IND_LAST_L_FRAME, 0, 1 );
            }
            else
            {
                push_indice( st, IND_LAST_L_FRAME, 1, 1 );
            }
        }
        else
        {
            push_indice( st, IND_HQ_SWITCHING_FLG, 0, 1 );
        }

        return;
    }

    if( st->core == ACELP_CORE )
    {
        if( st->ppp_mode == 1 || st->nelp_mode == 1 )
        {
            /* 1 bit to distinguish between 2.8kbps PPP/NELP frame and SID frame */
            push_indice( st, IND_CORE, 0, 1 );

            /* SC-VBR: 0 - PPP_NB, 1 - PPP_WB, 2 - NELP_NB, 3 - NELP_WB */
            if ( coder_type == VOICED && st->bwidth == NB && st->ppp_mode == 1 )
            {
                push_indice( st, IND_PPP_NELP_MODE, 0, 2 );
            }
            else if ( coder_type == VOICED && st->bwidth != NB && st->ppp_mode == 1 )
            {
                push_indice( st, IND_PPP_NELP_MODE, 1, 2 );
            }
            else if ( coder_type == UNVOICED && st->bwidth == NB && st->nelp_mode == 1 )
            {
                push_indice( st, IND_PPP_NELP_MODE, 2, 2);
            }
            else if ( coder_type == UNVOICED && st->bwidth != NB && st->nelp_mode == 1 )
            {
                push_indice( st, IND_PPP_NELP_MODE, 3, 2 );
            }
        }
        else if( st->core_brate != SID_2k40 && st->core_brate != FRAME_NO_DATA )
        {
            /* write the ACELP/HQ core selection bit */
            if (st->total_brate >= ACELP_24k40 )
            {
                push_indice( st, IND_CORE, 0, 1 );
            }

            /* find the section in the ACELP signalling table corresponding to bitrate */
            idx = 0;
            while ( acelp_sig_tbl[idx] != st->total_brate )
            {
                idx++;
            }

            /* retrieve the number of bits for signalling */
            nBits = (short) acelp_sig_tbl[++idx];

            /* retrieve the signalling index */
            start_idx = ++idx;
            while( acelp_sig_tbl[idx] != SIG2IND(coder_type, st->bwidth, sharpFlag, st->rf_mode) )
            {
                idx++;
            }

            push_indice( st, IND_ACELP_SIGNALLING, idx - start_idx, nBits );
        }

        /* write extension layer flag to distinguish between TBE (0) and BWE (1) */
        if( st->extl_brate > 0 )
        {
            if( st->extl == WB_TBE || st->extl == SWB_TBE || st->extl == FB_TBE )
            {
                push_indice( st, IND_BWE_FLAG, 0, 1 );
            }
            else if( st->extl == WB_BWE || st->extl == SWB_BWE || st->extl == FB_BWE )
            {
                push_indice( st, IND_BWE_FLAG, 1, 1 );
            }
        }
    }
    else  /* HQ core */
    {
        /* write ACELP->HQ core switching flag */
        if( st->last_core == ACELP_CORE || st->last_core == AMR_WB_CORE )
        {
            push_indice( st, IND_HQ_SWITCHING_FLG, 1, 1 );

            /* write ACELP L_frame info */
            if( st->last_L_frame == L_FRAME )
            {
                push_indice( st, IND_LAST_L_FRAME, 0, 1 );
            }
            else
            {
                push_indice( st, IND_LAST_L_FRAME, 1, 1 );
            }
        }
        else
        {
            push_indice( st, IND_HQ_SWITCHING_FLG, 0, 1 );
        }

        /* HQ/TCX core switching flag */
        push_indice( st, IND_MDCT_CORE, 0, 1 );

        /* Use ACELP signaling for LR MDCT */
        if ( st->total_brate <= ACELP_16k40 )
        {
            /* find the section in the ACELP signalling table corresponding to bitrate */
            idx = 0;
            while ( acelp_sig_tbl[idx] != st->total_brate )
            {
                idx++;
            }

            /* retrieve the number of bits for signalling */
            nBits = (short) acelp_sig_tbl[++idx];

            /* retrieve the signalling index */
            start_idx = ++idx;
            while( acelp_sig_tbl[idx] != SIG2IND(LR_MDCT, st->bwidth, 0, 0) )
            {
                idx++;
            }

            push_indice( st, IND_ACELP_SIGNALLING, idx - start_idx, nBits );
        }
        else
        {
            if( st->core_brate <= ACELP_64k )
            {
                /* write ACELP/HQ core indication flag */
                push_indice( st, IND_CORE, 1, 1 );
            }

            /* write band-width (needed for different I/O sampling rate support) */
            if( st->bwidth == NB )
            {
                push_indice( st, IND_HQ_BWIDTH, 0, 2 );
            }
            else if( st->bwidth == WB )
            {
                push_indice( st, IND_HQ_BWIDTH, 1, 2 );
            }
            else if( st->bwidth == SWB )
            {
                push_indice( st, IND_HQ_BWIDTH, 2, 2 );
            }
            else  /* st->bwidth == FB */
            {
                push_indice( st, IND_HQ_BWIDTH, 3, 2 );
            }
        }
    }

    return;
}

/*---------------------------------------------------------------------*
 * signalling_enc_rf()
 *
 * write channel-aware signalling information into the bit-stream
 *---------------------------------------------------------------------*/

void signalling_enc_rf(
    Encoder_State *st               /* i/o: encoder state structure   */
)
{
    short i, tmp_rf, sfr;

    /* write partial copy into bitstream */
    if( st->rf_mode == 1 )
    {
        enc_prm_rf( st, st->rf_indx_frametype[st->rf_fec_offset], st->rf_fec_offset );
        st->rf_indx_tbeGainFr[0] = st->RF_bwe_gainFr_ind;
    }

    /* Shift the RF indices such that the partial copy associated with
       (n-fec_offset)th frame is included in the bitstream in nth frame. */
    tmp_rf = st->rf_fec_offset;
    for(i = tmp_rf; i >= 0 ; i--)
    {
        /* rf frame type */
        st->rf_indx_frametype[i+1] = st->rf_indx_frametype[i];

        /* rf target bits buffer */
        st->rf_targetbits_buff[i+1] = st->rf_targetbits_buff[i];

        /* lsf indx */
        st->rf_indx_lsf[i+1][0] = st->rf_indx_lsf[i][0];
        st->rf_indx_lsf[i+1][1] = st->rf_indx_lsf[i][1];
        st->rf_indx_lsf[i+1][2] = st->rf_indx_lsf[i][2];

        /* ES pred energy */
        st->rf_indx_EsPred[i+1] = st->rf_indx_EsPred[i];

        /* LTF mode, sfr params: pitch, fcb and gain */
        for(sfr = 0; sfr < st->nb_subfr; sfr++)
        {
            st->rf_indx_ltfMode[i+1][sfr] = st->rf_indx_ltfMode[i][sfr];
            st->rf_indx_pitch[i+1][sfr] = st->rf_indx_pitch[i][sfr];
            st->rf_indx_fcb[i+1][sfr] = st->rf_indx_fcb[i][sfr];
            st->rf_indx_gain[i+1][sfr] = st->rf_indx_gain[i][sfr];
        }

        /* shift the nelp indices */
        st->rf_indx_nelp_iG1[i+1] = st->rf_indx_nelp_iG1[i];
        st->rf_indx_nelp_iG2[i+1][0] = st->rf_indx_nelp_iG2[i][0];
        st->rf_indx_nelp_iG2[i+1][1] = st->rf_indx_nelp_iG2[i][1];
        st->rf_indx_nelp_fid[i+1] = st->rf_indx_nelp_fid[i];

        /* tbe gain Fr shift */
        st->rf_indx_tbeGainFr[i+1] = st->rf_indx_tbeGainFr[i];
        st->rf_clas[i+1] = st->rf_clas[i];
        st->rf_gain_tcx[i+1] = st->rf_gain_tcx[i];
        st->rf_tcxltp_param[i+1] = st->rf_tcxltp_param[i];
    }

    return;
}
