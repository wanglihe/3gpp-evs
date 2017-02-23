/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <limits.h>
#include <string.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"


/*-------------------------------------------------------------------*
 * enc_prm_hm()
 *
 *
 *-------------------------------------------------------------------*/

static void enc_prm_hm(
    int *prm_hm,
    Encoder_State *st,
    int L_frame
)
{
    /* Disable HM for non-GC,VC modes */
    if( st->tcx_cfg.coder_type != VOICED && st->tcx_cfg.coder_type != GENERIC )
    {
        return;
    }

    /* Flag */
    push_next_indice(st, prm_hm[0], 1);

    if( prm_hm[0] )
    {
        /* Periodicy index */
        EncodeIndex( L_frame >= 256, prm_hm[1], st );

        if( st->tcx_cfg.coder_type == VOICED )
        {
            /* Gain index */
            push_next_indice(st, prm_hm[2], kTcxHmNumGainBits);
        }
    }

    return;
}

/*-----------------------------------------------------------------*
 * Function  enc_prm_rf()                                          *
 * ~~~~~~~~~~~~~~~~~~~~~~                                          *
 *                                                                 *
 * encode RF parameters for ACELP and TCX partial copy             *
 *-----------------------------------------------------------------*/

void enc_prm_rf( Encoder_State *st,
                 const short rf_frame_type,
                 const short fec_offset
               )
{
    short sfr, nb_subfr, n, index;
    short ltp_mode, ltf_mode, gains_mode;

    nb_subfr = st->nb_subfr;

    /* partial copy bitstream writing */
    if ( rf_frame_type >= RF_TCXFD && rf_frame_type <= RF_TCXTD2)
    {
        /* TCX frames partial copy write */
        if(rf_frame_type == RF_TCXFD)
        {
            push_next_indice(st, st->rf_indx_lsf[fec_offset][0], lsf_numbits[0]);  /* VQ 1 */
            push_next_indice(st, st->rf_indx_lsf[fec_offset][1], lsf_numbits[1]);  /* VQ 2 */
            push_next_indice(st, st->rf_indx_lsf[fec_offset][2], lsf_numbits[2]);  /* VQ 3 */
        }

        /* classification */
        if( st->rf_clas[fec_offset] == UNVOICED_CLAS )
        {
            index = 0;
        }
        else if( (st->rf_clas[fec_offset] == VOICED_TRANSITION) || (st->rf_clas[fec_offset] == UNVOICED_TRANSITION) )
        {
            index = 1;
        }
        else if( st->rf_clas[fec_offset] == VOICED_CLAS )
        {
            index = 2;
        }
        else
        {
            index = 3;
        }
        push_next_indice(st, index, 2);

        if(rf_frame_type == RF_TCXFD)
        {
            /* TCX global gain  = 7 bits */
            push_next_indice(st, st->rf_gain_tcx[fec_offset], 7);
        }
        else
        {
            /* LTP data */
            if ( (rf_frame_type == RF_TCXTD1 || rf_frame_type == RF_TCXTD2) && st->tcxltp )
            {
                push_next_indice(st, st->rf_tcxltp_param[fec_offset], 9);
            }
        }
    }
    else if( rf_frame_type == 7 )   /* NELP bitstream writing */
    {
        /* LSF indices */
        push_next_indice(st, st->rf_indx_lsf[fec_offset][0], 8);  /* VQ 1 */
        push_next_indice(st, st->rf_indx_lsf[fec_offset][1], 8);  /* VQ 2 */

        /* NELP gain indices */
        push_next_indice( st, st->rf_indx_nelp_iG1[fec_offset], 5 );
        push_next_indice( st, st->rf_indx_nelp_iG2[fec_offset][0], 6 );
        push_next_indice( st, st->rf_indx_nelp_iG2[fec_offset][1], 6 );

        /* NELP filter selection index */
        push_next_indice( st, st->rf_indx_nelp_fid[fec_offset], 2 );

        /* tbe gainFr */
        push_next_indice( st, st->rf_indx_tbeGainFr[fec_offset], 5 );
    }
    else if ( rf_frame_type >= 4 ) /* rf_frame_type ALL_PRED: 4, NO_PRED: 5, GEN_PRED: 6 */
    {
        /* LSF indices */
        push_next_indice(st, st->rf_indx_lsf[fec_offset][0], 8);  /* VQ 1 */
        push_next_indice(st, st->rf_indx_lsf[fec_offset][1], 8);  /* VQ 2 */

        /* ES pred */
        push_next_indice(st, st->rf_indx_EsPred[fec_offset], 3);

        ltp_mode = ACELP_LTP_MODE[1][1][rf_frame_type];
        ltf_mode = ACELP_LTF_MODE[1][1][rf_frame_type];
        gains_mode = ACELP_GAINS_MODE[1][1][rf_frame_type];

        /* Subframe parameters */
        for( sfr = 0; sfr < nb_subfr; sfr++ )
        {
            /* Pitch lag (5, or 8 bits) */
            n = ACELP_LTP_BITS_SFR[ltp_mode][sfr];
            if (n != 0)
            {
                push_next_indice(st, st->rf_indx_pitch[fec_offset][sfr], n);
            }

            /* Adaptive codebook filtering (1 bit) */
            if( ltf_mode == 2 )
            {
                push_next_indice(st, st->rf_indx_ltfMode[fec_offset][sfr], 1);
            }

            /*Innovative codebook*/
            if( (rf_frame_type == RF_NOPRED)
                    || ( rf_frame_type == RF_GENPRED && (sfr == 0 || sfr == 2)) )
            {
                push_next_indice(st, st->rf_indx_fcb[fec_offset][sfr], 7);
            }

            /* Gains (5b, 6b or 7b / subfr) */
            if( sfr == 0 || sfr == 2 )
            {
                n = ACELP_GAINS_BITS[gains_mode];
                push_next_indice(st, st->rf_indx_gain[fec_offset][sfr], n);
            }
        }
        /* tbe gainFr */
        push_next_indice( st, st->rf_indx_tbeGainFr[fec_offset], 2 );
    }

    /***************/
    /*IMPORTANT: The last three bits are always the rf_frame_type in the bitstream (for both acelp and tcx partial copy);
                 the rf_frame_type indicates the length of the partial copy payload at the decoder.
                 The 2 bits before the rf_frame_type contains the fec_offset */

    /***************/
    /* write FEC offset just before the rf_frame_type */
    if(fec_offset == 2 )
    {
        push_next_indice(st, 0, 2);
    }
    else if(fec_offset == 3 || fec_offset == 5 || fec_offset == 7)
    {
        push_next_indice(st, (fec_offset - 1)/2, 2);
    }

    /* write RF frame type last in the bitstream */
    push_next_indice(st, rf_frame_type, 3);

}

/*-----------------------------------------------------------------*
 * Function  enc_prm()                                             *
 * ~~~~~~~~~~~~~~~~~~~~~~                                          *
 *                                                                 *
 * encode parameters according to selected mode including          *
 * the FAC parameters when transition occurs.                      *
 *-----------------------------------------------------------------*/

void enc_prm(
    const short coder_type,   /* i  : coding type                     */
    int param[],        /* i  : parameters                      */
    int param_lpc[],    /* i  : LPC parameters                  */
    Encoder_State *st,        /* i/o: quantization Analysis values    */
    const short L_frame,      /* i  : frame length                    */
    CONTEXT_HM_CONFIG hm_cfg[],
    short * bits_param_lpc,
    short no_param_lpc
)
{
    short j, k, n, sfr, core, last_core;
    int *prm;
    short nbits_start, total_nbbits, nbits_header, nbits_lpc, nbits_tcx;
    short lg, nb_subfr;
    int lgFB;
    int nTnsParams, nTnsBits;
    int ix, j_old, wordcnt, bitcnt;
    short hm_size;
    short numlpc, index;
    short flag_ctx_hm;

    int idx;
    int start_idx;
    short nBits;

    /*--------------------------------------------------------------------------------*
     * initialization
     *--------------------------------------------------------------------------------*/

    nbits_lpc = 0;
    nbits_tcx = 0;

    flag_ctx_hm = 0;

    /* Useful parameters */
    nb_subfr = st->nb_subfr;
    core = st->core;
    last_core = st->last_core;

    /* Initialize pointers */
    prm = param;

    /* Init counters */
    j = 0;
    nbits_start = st->nb_bits_tot;

    /*--------------------------------------------------------------------------------*
     * HEADER
     *--------------------------------------------------------------------------------*/

    if (st->mdct_sw == MODE1)
    {
        /* Adjust st->bits_frame_core not to subtract MODE2 bandwidth signaling */
        st->bits_frame_core += FrameSizeConfig[st->frame_size_index].bandwidth_bits;

        /* Write MODE1 core mode signaling */
        signalling_mode1_tcx20_enc(st, 1);
    }

    /* EVS header */
    /* Modes (ACE_GC, ACE_UC, TCX20, TCX10...) */
    if ( st->tcxonly )
    {
        push_next_indice(st, core == TCX_10_CORE, 1);
        {
            if( st->clas == UNVOICED_CLAS )
            {
                index = 0;
            }
            else if( (st->clas == VOICED_TRANSITION) || (st->clas == UNVOICED_TRANSITION) )
            {
                index = 1;
            }
            else if( st->clas == VOICED_CLAS )
            {
                index = 2;
            }
            else
            {
                index = 3;
            }
            push_next_indice(st, index, 2);
        }
    }
    else
    {
        if ( core == ACELP_CORE )
        {
            /* write the RF signalling information */
            if( st->rf_mode ==1 )
            {
                /* find the section in the ACELP signalling table corresponding to bitrate */
                idx = 0;
                while ( acelp_sig_tbl[idx] != st->total_brate  )  /* total bitrate is kept at 13.2kbps */
                {
                    idx++;
                }

                /* retrieve the number of bits for signalling */
                nBits = (short) acelp_sig_tbl[++idx];

                /* retrieve the signalling index */
                start_idx = ++idx;
                while( acelp_sig_tbl[idx] != SIG2IND(coder_type, st->bwidth, st->sharpFlag, st->rf_mode) )
                {
                    idx++;
                }
                push_next_indice( st, idx - start_idx, nBits);
                push_next_indice( st, 0, 1); /* Indicate to the decoder that the core is ACELP*/
                nbits_start = 3;
            }
            else
            {
                push_next_indice( st, coder_type, 3 );
            }
        }
        else
        {
            if (st->mdct_sw == MODE1)
            {
                /* 2 bits instead of 3 as TCX is already signaled */
                push_next_indice(st, st->tcx_cfg.coder_type, 2 );
            }
            else
            {
                if (st->mdct_sw_enable == MODE2)
                {
                    push_next_indice(st, 1, 1); /* TCX */
                    push_next_indice(st, 0, 1); /* not HQ_CORE */
                    push_next_indice(st, st->tcx_cfg.coder_type, 2);
                }
                else
                {
                    /*write the RF signalling information*/
                    if( st->rf_mode ==1 )
                    {
                        /* find the section in the ACELP signalling table corresponding to bitrate */
                        idx = 0;
                        while ( acelp_sig_tbl[idx] != st->total_brate )
                        {
                            idx++;
                        }

                        /* retrieve the number of bits for signalling */
                        nBits = (short) acelp_sig_tbl[++idx];

                        if(st->tcx_cfg.coder_type == VOICED || st->tcx_cfg.coder_type == GENERIC || st->tcx_cfg.coder_type == TRANSITION)
                        {
                            st->sharpFlag=1;
                        }
                        else
                        {
                            st->sharpFlag=0;
                        }

                        /* retrieve the signalling index */
                        start_idx = ++idx;
                        while( acelp_sig_tbl[idx] != SIG2IND(st->tcx_cfg.coder_type, st->bwidth, st->sharpFlag, st->rf_mode) )
                        {
                            idx++;
                        }
                        push_next_indice( st, idx - start_idx, nBits);
                        push_next_indice( st, 1, 1); /* Indicate to the decoder that the core is TCX*/
                        nbits_start = 3;
                    }
                    else
                    {
                        push_next_indice(st, ACELP_MODE_MAX + st->tcx_cfg.coder_type, 3 );
                    }
                }
            }
        }
    }

    /* Encode previous mode for error concealment */
    if( !( core == ACELP_CORE && st->tcx_cfg.lfacNext<=0 ) )
    {
        push_next_indice(st, ((last_core!=ACELP_CORE) || (core==TCX_10_CORE)), 1);
    }

    /* write TCX overlap mode (1 bit: full, 2 bits: half or no overlap) */
    if( core != ACELP_CORE )
    {
        int overlap_code;
        assert(st->tcx_cfg.tcx_curr_overlap_mode != NOT_SUPPORTED && st->tcx_cfg.tcx_curr_overlap_mode <= ALDO_WINDOW && st->tcx_cfg.tcx_curr_overlap_mode >= FULL_OVERLAP);  /*1 is not allowed!*/
        if (st->tcx_cfg.tcx_curr_overlap_mode == MIN_OVERLAP)
        {
            nbits_tcx = 2;
            overlap_code = 2;
        }
        else if (st->tcx_cfg.tcx_curr_overlap_mode == HALF_OVERLAP)
        {
            nbits_tcx = 2;
            overlap_code = 3;
        }
        else
        {
            nbits_tcx = 1;
            overlap_code = 0;
        }
        push_next_indice(st, overlap_code, nbits_tcx);
    }

    if( st->plcExt.enableGplc )
    {
        /* encode side information. */
        enc_prm_side_Info( &st->plcExt, st );
    }

    if( st->glr )
    {
        if( core != ACELP_CORE || coder_type == INACTIVE ||
                (st->last_core == ACELP_CORE && st->last_coder_type_raw == INACTIVE) ||
                st->glr_reset )
        {
            st->glr_idx[0] = 0;
        }

        if( core == ACELP_CORE )
        {
            push_next_indice(st, st->glr_idx[0], G_LPC_RECOVERY_BITS);
        }
    }

    st->glr_reset = 0;

    nbits_header = st->nb_bits_tot - nbits_start;

    /*--------------------------------------------------------------------------------*
     * LPC PARAMETERS
     *--------------------------------------------------------------------------------*/

    if( st->enableTcxLpc && core != ACELP_CORE )
    {
        /* Encode the indices */
        nbits_lpc = enc_lsf_tcxlpc(&param_lpc, st);
    }
    else
    {
        if( st->lpcQuantization == 0 )
        {
            /* LPC quantizer */
            if(core == TCX_20_CORE)
            {
                numlpc = 1;
            }
            else
            {
                numlpc = 2;
            }

            nbits_lpc = encode_lpc_avq(st, numlpc, param_lpc, core);
        }
        else if( st->lpcQuantization == 1 )
        {
            if(st->sr_core == 16000 && coder_type == VOICED
                    && core==ACELP_CORE
              )
            {
                nbits_lpc = lsf_bctcvq_encprm(st, param_lpc, bits_param_lpc, no_param_lpc);
            }
            else
            {
                nbits_lpc = lsf_msvq_ma_encprm(st, param_lpc, core, coder_type, st->acelp_cfg.midLpc, bits_param_lpc, no_param_lpc );
            }
        }
        else
        {
            assert(0);
        }
    }


    /*--------------------------------------------------------------------------------*
     * PRINT BIT ALLOCATION
     *--------------------------------------------------------------------------------*/



    /*--------------------------------------------------------------------------------*
     * ACELP
     *--------------------------------------------------------------------------------*/
    if( core == ACELP_CORE )
    {
        /* Adaptive BPF (2 bits)*/
        n = ACELP_BPF_BITS[st->acelp_cfg.bpf_mode];

        if( n!=0 )
        {
            push_next_indice(st, st->bpf_gain_param, n);
        }

        /* Mean energy (2 or 3 bits) */
        n = ACELP_NRG_BITS[st->acelp_cfg.nrg_mode];

        if( n!=0 )
        {
            push_next_indice(st, prm[j++], n);
        }

        /* Subframe parameters */
        for( sfr=0; sfr<nb_subfr; sfr++ )
        {
            /* Pitch lag (4, 5, 6, 8 or 9 bits) */
            n = ACELP_LTP_BITS_SFR[st->acelp_cfg.ltp_mode][sfr];

            if (n!=0)
            {
                push_next_indice(st, prm[j++], n);
            }

            /* Adaptive codebook filtering (1 bit) */
            if( st->acelp_cfg.ltf_mode == 2 )
            {
                push_next_indice(st, prm[j++], 1);
            }

            /*Innovative codebook*/

            j_old = j;
            if ((st->acelp_cfg.fixed_cdk_index[sfr] >= ACELP_FIXED_CDK_NB) || (st->acelp_cfg.fixed_cdk_index[sfr] < 0))
            {
                fprintf(stderr,"ACELP bits allocation: wrong fixed cdk bit allocation\n");
                assert(0);
            }

            wordcnt = ACELP_FIXED_CDK_BITS(st->acelp_cfg.fixed_cdk_index[sfr]) >> 4;
            bitcnt = ACELP_FIXED_CDK_BITS(st->acelp_cfg.fixed_cdk_index[sfr]) & 15;

            for (ix = 0; ix < wordcnt; ix++)
            {
                push_next_indice(st, prm[j++], 16);
            }

            if (bitcnt)
            {
                push_next_indice(st, prm[j++], bitcnt);
            }

            j = j_old+8;

            /* Gains (5b, 6b or 7b / subfr) */
            n=ACELP_GAINS_BITS[st->acelp_cfg.gains_mode[sfr]];
            push_next_indice(st, prm[j++], n);
        }

    }


    /*--------------------------------------------------------------------------------*
     * TCX20
     *--------------------------------------------------------------------------------*/

    if( core == TCX_20_CORE )
    {
        if( st->enablePlcWaveadjust )
        {
            short * index = NULL;
            index = &st->Tonal_SideInfo;
            push_next_indice(st, *index, 1);
        }

        /* TCX Gain = 7 bits */
        push_next_indice(st, prm[j++], 7);

        /* TCX Noise Filling = NBITS_NOISE_FILL_LEVEL bits */
        push_next_indice(st, prm[j++], NBITS_NOISE_FILL_LEVEL);

        /* LTP data */
        if ( st->tcxltp || (st->sr_core > 25600) ) /* PLC pitch info for HB */
        {
            if ( prm[j] )
            {
                push_next_indice(st, 1       , 1);
                push_next_indice(st, prm[j+1], 9);
                push_next_indice(st, prm[j+2], 2);
            }
            else
            {
                push_next_indice(st, 0, 1);
            }
        }

        j += 3;

        /* TCX spectral data */
        lg = L_frame;
        lgFB = st->tcx_cfg.tcx_coded_lines;
        if ( last_core == ACELP_CORE )
        {
            /* ACE->TCX transition */
            lg += st->tcx_cfg.tcx_offset;
            lgFB += lgFB >> 2;
            if(st->tcx_cfg.lfacNext<0)
            {
                lg -= st->tcx_cfg.lfacNext;
            }
        }

        /* TNS data */
        nTnsParams = 0;
        nTnsBits = 0;

        if (st->tcx_cfg.fIsTNSAllowed)
        {
            WriteTnsData(st->tcx_cfg.pCurrentTnsConfig, prm+j, &nTnsParams, st, &nTnsBits);
            j += nTnsParams;
        }

        hm_size = (int)(2.0f*st->tcx_cfg.bandwidth*(float)lg);

        if (st->tcx_lpc_shaped_ari && last_core != ACELP_CORE)
        {
            enc_prm_hm( &prm[j], st, hm_size );
        }

        /*Context HM flag*/
        if ( st->tcx_cfg.ctx_hm && (last_core != ACELP_CORE) )
        {
            push_next_indice(st, prm[j], 1);

            if (prm[j])
            {
                EncodeIndex(hm_size >= 256, prm[j+1], st);

                flag_ctx_hm=1;
            }
        }
        j += NPRM_CTX_HM;
        /* IGF data */
        if (st->igf)
        {
            st->hIGFEnc.infoTotalBitsPerFrameWritten = 0;
            IGFEncWriteBitstream( &st->hIGFEnc, st, &st->hIGFEnc.infoTotalBitsPerFrameWritten,
                                  (st->last_core == ACELP_CORE)?IGF_GRID_LB_TRAN:IGF_GRID_LB_NORM, 1 );
        }

        total_nbbits = st->nb_bits_tot - nbits_start;
        if(st->rf_mode)
        {
            total_nbbits += st->rf_target_bits_write;
        }
        nbits_tcx = st->bits_frame_core - total_nbbits;

        if (st->tcx_lpc_shaped_ari)
        {
            push_next_bits(st, &prm[++j], nbits_tcx);
            j += nbits_tcx;
        }
        else
        {
            ACcontextMapping_encode2_no_mem_s17_LC( st, prm+j, lgFB, prm[j-1], /* lastnz */
                                                    nbits_tcx, NPRM_RESQ * st->tcx_cfg.resq, flag_ctx_hm?hm_cfg:NULL );
        }

    }


    /*--------------------------------------------------------------------------------*
     * TCX10
     *--------------------------------------------------------------------------------*/

    if( core == TCX_10_CORE )
    {
        int nbits_igf = 0;
        if (st->igf)
        {
            nbits_igf = IGFEncWriteConcatenatedBitstream( &st->hIGFEnc, st );
        }

        for (k = 0; k < 2; k++)
        {
            flag_ctx_hm = 0;

            prm = param + (k*NPRM_DIV);
            j = 0;

            nbits_tcx = total_nbbits = st->nb_bits_tot - nbits_start;

            if(st->enablePlcWaveadjust && k)
            {
                short * index = NULL;
                index = &st->Tonal_SideInfo;
                push_next_indice(st, *index, 1);
            }
            /* TCX Gain = 7 bits */
            push_next_indice(st, prm[j++], 7);

            /* TCX Noise Filling = NBITS_NOISE_FILL_LEVEL bits */
            push_next_indice(st, prm[j++], NBITS_NOISE_FILL_LEVEL);

            /* LTP data */
            if ( (k == 0) && (st->tcxltp || (st->sr_core > 25600) ) ) /* PLC pitch info for HB */
            {
                if ( prm[j] )
                {
                    push_next_indice(st, 1       , 1);
                    push_next_indice(st, prm[j+1], 9);
                    push_next_indice(st, prm[j+2], 2);
                }
                else
                {
                    push_next_indice(st, 0       , 1);
                }
            }
            j += 3;

            /* TCX spectral data */
            lg = L_frame >> 1;
            lgFB = st->tcx_cfg.tcx_coded_lines >> 1;

            if ( k==0 && last_core == ACELP_CORE )
            {
                /* ACE->TCX transition */
                lg += st->tcx_cfg.tcx_offset;
                lgFB += lgFB >> 1;

                if(st->tcx_cfg.lfacNext<0)
                {
                    lg -= st->tcx_cfg.lfacNext;
                }
            }

            /* TNS data */
            nTnsParams = 0;
            nTnsBits = 0;

            if( st->tcx_cfg.fIsTNSAllowed )
            {
                SetTnsConfig(&st->tcx_cfg, 0, (last_core == ACELP_CORE) && (k == 0));
                WriteTnsData(st->tcx_cfg.pCurrentTnsConfig, prm+j, &nTnsParams, st, &nTnsBits);
                j += nTnsParams;
            }

            hm_size = (int)(2.0f*st->tcx_cfg.bandwidth*(float)lgFB);

            /*Context HM flag*/
            if ( st->tcx_cfg.ctx_hm && !(last_core == ACELP_CORE && k == 0) )
            {
                push_next_indice(st, prm[j], 1);

                if (prm[j])
                {
                    EncodeIndex(hm_size >= 256, prm[j+1], st);
                    flag_ctx_hm=1;
                }
            }
            j += NPRM_CTX_HM;

            total_nbbits = st->nb_bits_tot - nbits_start;
            nbits_tcx = (st->bits_frame_core - nbits_header - nbits_lpc - nbits_igf + 1 - k) / 2 - (total_nbbits - nbits_tcx);

            ACcontextMapping_encode2_no_mem_s17_LC( st, prm+j, lgFB, prm[j-1], /* lastnz */nbits_tcx,
                                                    NPRM_RESQ * st->tcx_cfg.resq, flag_ctx_hm?&hm_cfg[k]:NULL );

        } /* k, window index */
    }


    /*--------------------------------------------------------------------------------*
     * END
     *--------------------------------------------------------------------------------*/

    total_nbbits = st->nb_bits_tot - nbits_start;



    return;
}

