/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"

/*-----------------------------------------------------------------*
 * dec_prm_hm()
 *
 *
 *-----------------------------------------------------------------*/

static void dec_prm_hm(
    Decoder_State *st,
    int *prm_hm,
    const short L_frame
)
{
    /* Disable HM for non-GC,VC modes */
    if( st->tcx_cfg.coder_type != VOICED && st->tcx_cfg.coder_type != GENERIC )
    {
        prm_hm[0] = 0;
        return;
    }

    prm_hm[1] = -1;
    prm_hm[2] = 0;

    /* Flag */
    prm_hm[0] = get_next_indice(st, 1);

    if (prm_hm[0])
    {
        /* Periodicity index */
        DecodeIndex( st, L_frame >= 256, &prm_hm[1] );

        /* Gain index */
        if( st->tcx_cfg.coder_type == VOICED )
        {
            prm_hm[2] = get_next_indice(st, kTcxHmNumGainBits );
        }
    }

    return;
}


/*-----------------------------------------------------------------*
 * Function  dec_prm()                                             *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~                                   *
 *
 * SQ is used for TCX modes
 *
 * decode parameters according to selected mode                    *
 *-----------------------------------------------------------------*/

void dec_prm(
    short *core,            /* o : current frame mode               */
    short *last_core,       /* o : last frame mode                  */
    short *coder_type,      /* o : coder type                       */
    int param[],            /* o : decoded parameters               */
    int param_lpc[],        /* o : LPC parameters                   */
    Word16 *total_nbbits,   /* i/o : number of bits / decoded bits  */
    Decoder_State *st,      /* i/o:  decoder memory state           */
    const int L_frame,
    short *bitsRead
)
{
    short j, k, n, sfr;
    int *prm;
    int lg;
    int lgFB;
    int start_bit_pos, bits_common;
    int acelp_target_bits=-1;
    int tmp;
    int nTnsParams;
    int nTnsBits;
    int nb_subfr;
    int nbits_tcx;
    int ix, j_old, wordcnt, bitcnt;
    int hm_size;
    int *prms;
    CONTEXT_HM_CONFIG hm_cfg;
    int indexBuffer[N_MAX+1];
    short flag_ctx_hm;
    short ind;
    short ltp_mode, gains_mode;
    int prm_ltp[LTPSIZE];

    /*--------------------------------------------------------------------------------*
     * INIT
     *--------------------------------------------------------------------------------*/

    hm_cfg.indexBuffer = indexBuffer;

    if (st->mdct_sw == MODE1)
    {
        start_bit_pos = 0; /* count from frame start */

        /* Adjust st->bits_frame_core not to subtract MODE2 bandwidth signaling */
        for (n=0; n<FRAME_SIZE_NB; n++)
        {
            if (FrameSizeConfig[n].frame_bits == st->total_brate/50)
            {
                st->bits_frame_core += FrameSizeConfig[n].bandwidth_bits;
                break;
            }
        }
    }
    else
    {
        if( st->rf_flag == 1 )
        {
            /*Inherent adjustment to accommodate the compact packing used in the RF mode*/
            start_bit_pos = st->next_bit_pos - 2;
        }
        else
        {
            start_bit_pos = st->next_bit_pos;
        }
    }

    /* Framing parameters */
    nb_subfr = st->nb_subfr;

    /* Initialize pointers */
    prm = param;

    /* Init counters */
    j = 0;

    /* Init LTP data */
    st->tcx_hm_LtpPitchLag = -1;
    st->tcxltp_gain = 0.0f;


    /*--------------------------------------------------------------------------------*
     * EVS HEADER
     *--------------------------------------------------------------------------------*/

    /* Modes (ACE_GC, ACE_UC, TCX20, TCX10...) */
    if ( st->tcxonly )
    {
        tmp = get_next_indice(st, 1);
        *core = tmp+1;

        ind = get_next_indice(st, 2);
        st->clas_dec = ONSET;
        if( ind == 0 )
        {
            st->clas_dec = UNVOICED_CLAS;
        }
        else if( ind == 1 )
        {
            if( st->last_good >= VOICED_TRANSITION )
            {
                st->clas_dec = VOICED_TRANSITION;
            }
            else
            {
                st->clas_dec = UNVOICED_TRANSITION;
            }
        }
        else if( ind == 2 )
        {
            st->clas_dec = VOICED_CLAS;
        }
        *coder_type = INACTIVE;
        st->VAD = 0;
    }
    else
    {
        if (st->mdct_sw == MODE1)
        {
            /* 2 bits instead of 3 as TCX is already signaled */
            *core = TCX_20_CORE;
            st->tcx_cfg.coder_type = get_next_indice(st, 2);
            *coder_type = st->tcx_cfg.coder_type;
        }
        else
        {
            if (st->mdct_sw_enable == MODE2)
            {
                if (get_next_indice_1(st)) /* TCX */
                {
                    tmp = get_next_indice(st, 3);
                    assert(!(tmp & 4) || !"HQ_CORE encountered in dec_prm");
                    *core = TCX_20_CORE;
                    st->tcx_cfg.coder_type = tmp;
                    *coder_type = st->tcx_cfg.coder_type;
                }
                else /* ACELP */
                {
                    *core = ACELP_CORE;
                    *coder_type = get_next_indice(st, 2);
                }
            }
            else
            {
                if(st->rf_flag == 1)
                {
                    if( !( st->use_partial_copy ) )
                    {
                        tmp = get_next_indice(st, 1);
                        if(tmp == 0)
                        {
                            *core = ACELP_CORE;
                        }
                        else
                        {
                            *core = TCX_20_CORE;
                            st->tcx_cfg.coder_type = *coder_type;
                        }
                    }
                }
                else
                {
                    tmp = get_next_indice(st, 3);
                    if( tmp < ACELP_MODE_MAX )
                    {
                        *core = ACELP_CORE;
                        *coder_type = tmp;
                    }
                    else
                    {
                        *core = TCX_20_CORE;
                        st->tcx_cfg.coder_type = tmp-ACELP_MODE_MAX;
                        *coder_type = st->tcx_cfg.coder_type;
                    }
                }
            }
        }

        if( st->igf && *core == ACELP_CORE )
        {
            st->bits_frame_core -= get_tbe_bits(st->total_brate, st->bwidth, st->rf_flag );
        }

        if( st->rf_flag )
        {
            st->bits_frame_core -= (st->rf_target_bits+1); /* +1 as flag-bit not considered in rf_target_bits */
        }

        /* Inactive frame detection on non-DTX mode */
        if( *coder_type == INACTIVE )
        {
            st->VAD = 0;
        }
        else
        {
            st->VAD = 1;
        }
    }

    /*Core extended mode mapping for correct PLC classification*/
    st->core_ext_mode=*coder_type;
    if( *coder_type == INACTIVE )
    {
        st->core_ext_mode = UNVOICED;
    }

    /* Decode previous mode for error concealment */
    if( !(*core==0 && st->tcx_cfg.lfacNext<=0) && !st->use_partial_copy )
    {
        tmp = get_next_indice(st, 1);
        *last_core = tmp;

        /*for TCX 10 force last_core to be TCX since ACELP as previous core is forbidden*/
        if(*core==TCX_10_CORE)
        {
            *last_core = 1;
        }
        if((st->prev_bfi == 0) && (st->last_core != *last_core))
        {
            st->BER_detect = 1;
        }
    }

    if(st->rf_flag && st->use_partial_copy && !st->tcxonly)
    {
        st->bits_frame_core = st->rf_target_bits;

        /* offset the indices to read the acelp partial copy */
        get_next_indice_tmp(st, start_bit_pos + st->total_brate/50 - st->rf_target_bits - 3 - st->next_bit_pos);
    }

    if( !st->use_partial_copy )
    {
        int overlap_code;

        /* Set the last overlap mode based on the previous and current frame type and coded overlap mode */
        if ((*last_core == ACELP_CORE) || (*last_core == AMR_WB_CORE))
        {
            st->tcx_cfg.tcx_last_overlap_mode = TRANSITION_OVERLAP;
        }
        else if ((*core == TCX_10_CORE) && (st->tcx_cfg.tcx_curr_overlap_mode == ALDO_WINDOW))
        {
            st->tcx_cfg.tcx_last_overlap_mode = FULL_OVERLAP;
        }
        else if ((*core != TCX_10_CORE) && (st->tcx_cfg.tcx_curr_overlap_mode == FULL_OVERLAP))
        {
            st->tcx_cfg.tcx_last_overlap_mode = ALDO_WINDOW;
        }
        else
        {
            st->tcx_cfg.tcx_last_overlap_mode = st->tcx_cfg.tcx_curr_overlap_mode;
        }

        /* Set the current overlap mode based on the current frame type and coded overlap mode */
        st->tcx_cfg.tcx_curr_overlap_mode = ALDO_WINDOW;

        if( *core != ACELP_CORE )
        {
            overlap_code = 0;
            if (get_next_indice(st, 1))
            {
                overlap_code = 2 + get_next_indice(st, 1);
            }
            assert(MIN_OVERLAP == 2 && HALF_OVERLAP == 3);
            st->tcx_cfg.tcx_curr_overlap_mode = overlap_code;

            /*TCX10 : always symmetric windows*/
            if ((*core == TCX_20_CORE) && (overlap_code == 0) &&
                    (*last_core != ACELP_CORE) && (*last_core != AMR_WB_CORE))
            {
                st->tcx_cfg.tcx_curr_overlap_mode = ALDO_WINDOW;
            }
        }

        if( st->enableGplc )
            /* SIDE INFO. DECODING */
        {
            short pitchDiff;
            short bits_per_subfr, search_range;
            bits_per_subfr = 4;
            search_range = 8;
            st->flagGuidedAcelp = get_next_indice(st, 1);

            pitchDiff = 0;
            if( st->flagGuidedAcelp )
            {
                pitchDiff = get_next_indice(st, bits_per_subfr);
                st->guidedT0 = (pitchDiff - search_range);
            }
            if( pitchDiff==0 && st->flagGuidedAcelp )
            {
                st->flagGuidedAcelp = 0;
            }
        }
        else
        {
            st->flagGuidedAcelp = 0;
        }

        if( st->dec_glr )
        {
            if( *core == ACELP_CORE )
            {
                st->dec_glr_idx = get_next_indice(st, G_LPC_RECOVERY_BITS);
            }
            else
            {
                st->dec_glr_idx = -1;
            }
        }
    }

    /*--------------------------------------------------------------------------------*
     * LPC PARAMETERS
     *--------------------------------------------------------------------------------*/

    /*Initialization of LPC Mid flag*/
    if( (st->lpcQuantization == 1 && *coder_type == VOICED) || (st->use_partial_copy))
    {
        (&(st->acelp_cfg))->midLpc = 0;
    }
    else
    {
        (&(st->acelp_cfg))->midLpc = st->acelp_cfg.midLpc_enable;
    }

    if( st->use_partial_copy == 0 )
    {
        /* Number of sets of LPC parameters (does not include mid-lpc) */
        if ( st->tcxonly==0 || *core < TCX_10_CORE )
        {
            st->numlpc = 1;
        }
        else
        {
            st->numlpc = 2;
        }

        /* Decode LPC parameters */
        if( st->enableTcxLpc && *core != ACELP_CORE )
        {
            int tcx_lpc_cdk;
            tcx_lpc_cdk = tcxlpc_get_cdk( *coder_type );
            dec_lsf_tcxlpc( st, &param_lpc, st->narrowBand, tcx_lpc_cdk );
        }
        else
        {
            if( st->lpcQuantization == 0 )
            {
                decode_lpc_avq( st, st->numlpc, param_lpc );
            }
            else if( st->lpcQuantization == 1 )
            {
                if(st->sr_core == 16000 && *coder_type == VOICED && *core == ACELP_CORE)
                {
                    lsf_bctcvq_decprm(st, param_lpc);
                }
                else
                {
                    lsf_msvq_ma_decprm( st, param_lpc, *core, *coder_type, st->acelp_cfg.midLpc, st->narrowBand, st->sr_core );
                }
            }
            else
            {
                assert(0);
            }
        }
    }
    else
    {
        st->numlpc = 1;

        if( st->rf_frame_type == RF_TCXFD )
        {
            param_lpc[0] = 0;
            param_lpc[1] = get_next_indice(st, lsf_numbits[0]);  /* VQ 1 */
            param_lpc[2] = get_next_indice(st, lsf_numbits[1]);  /* VQ 2 */
            param_lpc[3] = get_next_indice(st, lsf_numbits[2]);  /* VQ 3 */
        }
        else if( st->rf_frame_type >= RF_ALLPRED && st->rf_frame_type <= RF_NELP )
        {
            /* LSF indices */
            param_lpc[0] = get_next_indice(st, 8);  /* VQ 1 */
            param_lpc[1] = get_next_indice(st, 8);  /* VQ 2 */
        }
    }

    bits_common = st->next_bit_pos - start_bit_pos;

    /*--------------------------------------------------------------------------------*
     * ACELP
     *--------------------------------------------------------------------------------*/

    if( *core == ACELP_CORE && st->use_partial_copy == 0 )
    {
        /* Target Bits */

        /* needed in decoder to read the bitstream */
        acelp_target_bits = st->bits_frame_core - bits_common;

        /*Configure ACELP*/
        BITS_ALLOC_config_acelp( acelp_target_bits, *coder_type, &(st->acelp_cfg), st->narrowBand, st->nb_subfr );


        /* Adaptive BPF (2 bits)*/
        n = ACELP_BPF_BITS[st->acelp_cfg.bpf_mode];
        if( n!=0 )
        {
            st->bpf_gain_param = get_next_indice(st, n);
        }
        else
        {
            st->bpf_gain_param=(st->acelp_cfg.bpf_mode)*2;
        }

        /* Mean energy (2 or 3 bits) */
        n = ACELP_NRG_BITS[st->acelp_cfg.nrg_mode];
        if( n!=0 )
        {
            prm[j++] = get_next_indice(st, n);
        }

        /* Subframe parameters */
        for( sfr=0; sfr<nb_subfr; sfr++ )
        {
            /* Pitch lag (4, 5, 6, 8 or 9 bits) */
            n = ACELP_LTP_BITS_SFR[st->acelp_cfg.ltp_mode][sfr];

            if( n!=0 )
            {
                prm[j] = get_next_indice(st, n);
                j++;
            }

            /* Adaptive codebook filtering (1 bit) */
            if( st->acelp_cfg.ltf_mode == 2 )
            {
                prm[j] = get_next_indice(st, 1);
                j++;
            }

            /* Innovative codebook */
            {
                /* Decode pulse positions. */
                j_old = j;
                wordcnt = ACELP_FIXED_CDK_BITS(st->acelp_cfg.fixed_cdk_index[sfr])/16;
                bitcnt = ACELP_FIXED_CDK_BITS(st->acelp_cfg.fixed_cdk_index[sfr]) & 15;

                /* sanity check - can happen in case of bit errors */
                if ((st->acelp_cfg.fixed_cdk_index[sfr] >= ACELP_FIXED_CDK_NB) || (st->acelp_cfg.fixed_cdk_index[sfr] < 0))
                {
                    st->acelp_cfg.fixed_cdk_index[sfr] = 0;
                    st->BER_detect = 1;
                }

                for (ix = 0; ix < wordcnt; ix++)
                {
                    prm[j] = get_next_indice(st, 16);
                    j++;
                }

                if( bitcnt )
                {
                    prm[j] = get_next_indice(st, bitcnt);
                }

                j = j_old + 8;
            }

            /* Gains (5b, 6b or 7b / subfr) */
            n = ACELP_GAINS_BITS[st->acelp_cfg.gains_mode[sfr]];
            prm[j++] = get_next_indice(st, n);
        }/*end of subfr loop*/
    }
    else if ( st->rf_frame_type >= RF_ALLPRED  && st->use_partial_copy )
    {
        BITS_ALLOC_config_acelp( st->rf_target_bits,         /* target bits ranges from 56 to 72 depending on rf_type */
                                 st->rf_frame_type,          /* already offset by 4 to parse the config elements for partial copy */
                                 &(st->acelp_cfg),           /* acelp_cfg_rf*/
                                 0,                          /* is narrowBand */
                                 st->nb_subfr );

        /* rf_frame_type NELP: 7 */
        if(st->rf_frame_type == RF_NELP)
        {
            /* NELP gain indices */
            st->rf_indx_nelp_iG1 = get_next_indice( st, 5 );
            st->rf_indx_nelp_iG2[0] = get_next_indice( st, 6 );
            st->rf_indx_nelp_iG2[1] = get_next_indice( st, 6 );

            /* NELP filter selection index */
            st->rf_indx_nelp_fid = get_next_indice( st, 2 );

            /* tbe gainFr */
            st->rf_indx_tbeGainFr = get_next_indice( st, 5 );
        }
        else
        {
            /* rf_frame_type ALL_PRED: 4, NO_PRED: 5, GEN_PRED: 6*/
            /* ES pred */
            prm[j++] = get_next_indice(st, 3);

            ltp_mode = ACELP_LTP_MODE[1][1][st->rf_frame_type];
            gains_mode = ACELP_GAINS_MODE[1][1][st->rf_frame_type];

            /* Subframe parameters */
            for( sfr = 0; sfr < nb_subfr; sfr++ )
            {
                /* Pitch lag (5, or 8 bits) */
                n = ACELP_LTP_BITS_SFR[ltp_mode][sfr];
                if (n != 0)
                {
                    prm[j++] = get_next_indice(st, n);
                }


                /*Innovative codebook*/
                if( (st->rf_frame_type == RF_NOPRED)
                        || ( st->rf_frame_type == RF_GENPRED && (sfr == 0 || sfr == 2)) )
                {
                    /* NOTE: FCB actual bits need to be backed up as well */
                    /*n = ACELP_FIXED_CDK_BITS(st->rf_indx_fcb[fec_offset][sfr]) & 15;*/
                    prm[j] = get_next_indice(st, 7);
                    j = j + 8;
                }

                /* Gains (5b, 6b or 7b / subfr) */
                if( sfr == 0 || sfr == 2)
                {
                    n = ACELP_GAINS_BITS[gains_mode];
                    prm[j++] = get_next_indice(st, n);
                }
            }
            st->rf_indx_tbeGainFr = get_next_indice( st, 2 );
        }
    }

    /*--------------------------------------------------------------------------------*
     * TCX20
     *--------------------------------------------------------------------------------*/

    if( *core == TCX_20_CORE && st->use_partial_copy == 0 )
    {
        flag_ctx_hm = 0;

        if( st->enablePlcWaveadjust )
        {
            ind = get_next_indice(st, 1);
            st->tonality_flag = ind;
        }

        /* TCX Gain = 7 bits */
        prm[j++] = get_next_indice(st, 7);

        /* TCX Noise Filling = NBITS_NOISE_FILL_LEVEL bits */
        prm[j++] = get_next_indice(st, NBITS_NOISE_FILL_LEVEL);

        /* LTP data */
        if ( st->tcxltp || st->sr_core > 25600 )  /* PLC pitch info for HB */
        {

            prm[j] = get_next_indice(st, 1);
            if( prm[j] )
            {
                prm[j+1] = get_next_indice(st, 9);
                prm[j+2] = get_next_indice(st, 2);
            }
            st->BER_detect = st->BER_detect |
                             tcx_ltp_decode_params(
                                 &prm[j], &(st->tcxltp_pitch_int), &(st->tcxltp_pitch_fr), &(st->tcxltp_gain),
                                 st->pit_min, st->pit_fr1, st->pit_fr2, st->pit_max, st->pit_res_max );
            st->tcxltp_last_gain_unmodified = st->tcxltp_gain;
            st->tcx_hm_LtpPitchLag = ((!st->tcxonly) && (st->tcxltp_pitch_int < L_frame)
                                      ? (((2 * st->L_frame * st->pit_res_max) << kLtpHmFractionalResolution) / (st->tcxltp_pitch_int * st->pit_res_max + st->tcxltp_pitch_fr))
                                      : -1);
        }

        j += 3;

        /* TCX spectral data */
        lg = L_frame;
        lgFB = st->tcx_cfg.tcx_coded_lines;

        if( *last_core == ACELP_CORE )
        {
            /* ACE->TCX transition */
            lg += st->tcx_cfg.tcx_offset;
            lgFB += lgFB >> 2;

            if( st->tcx_cfg.lfacNext < 0 )
            {
                lg -= st->tcx_cfg.lfacNext;
            }
        }

        /* TNS data */
        nTnsParams = 0;
        nTnsBits = 0;

        if( st->tcx_cfg.fIsTNSAllowed )
        {
            SetTnsConfig(&st->tcx_cfg, 1, *last_core == ACELP_CORE);

            ReadTnsData(st->tcx_cfg.pCurrentTnsConfig, st, &nTnsBits, prm+j, &nTnsParams);

            j += nTnsParams;
        }

        hm_size = (int)(2.0f*st->TcxBandwidth*(float)lg);

        if( st->tcx_lpc_shaped_ari && *last_core != ACELP_CORE)
        {
            dec_prm_hm(st, &prm[j], hm_size );
        }

        nbits_tcx = ( st->bits_frame_core -  (st->next_bit_pos - start_bit_pos) );

        /*Context HM flag*/
        if ( st->tcx_cfg.ctx_hm && (*last_core != ACELP_CORE) )
        {
            prm[j] = get_next_indice(st, 1);
            nbits_tcx--;

            if( prm[j] )
            {
                int NumIndexBits = DecodeIndex( st, hm_size >= 256, prm+j+1 );

                flag_ctx_hm=1;

                ConfigureContextHm(
                    lgFB,
                    nbits_tcx,
                    *(prm+j+1),
                    st->tcx_hm_LtpPitchLag,
                    &hm_cfg);

                nbits_tcx-=NumIndexBits;
            }
        }
        j += NPRM_CTX_HM;
        /* read IGF payload */
        if (st->igf)
        {
            n = st->next_bit_pos;
            IGFDecReadLevel( &st->hIGFDec, st, (st->last_core == ACELP_CORE)?IGF_GRID_LB_TRAN:IGF_GRID_LB_NORM, 1 );

            IGFDecReadData( &st->hIGFDec, st, (st->last_core == ACELP_CORE)?IGF_GRID_LB_TRAN:IGF_GRID_LB_NORM, 1 );

            nbits_tcx  -= (st->next_bit_pos - n);
        }

        nbits_tcx = ( st->bits_frame_core -  (st->next_bit_pos - start_bit_pos) );
        if( st->tcx_lpc_shaped_ari )
        {
            prm[j++] = nbits_tcx;  /* store length of buffer */
            prms = &prm[j];
            for (ix = 0; ix < nbits_tcx; ix++)
            {
                prms[ix] = get_next_indice_1(st);
            }
            for (ix = 0; ix < 32; ix++)
            {
                prms[ix+nbits_tcx] = 1;
            }
            j += nbits_tcx;
        }
        else
        {

            st->resQBits[0] = ACcontextMapping_decode2_no_mem_s17_LC( st, prm+j,lgFB,
                              nbits_tcx, NPRM_RESQ*st->tcx_cfg.resq, flag_ctx_hm ? &hm_cfg : NULL );
            j += lg;
        }
    }
    if( st->rf_frame_type >= RF_TCXFD && st->rf_frame_type <= RF_TCXTD2 && st->use_partial_copy == 1 )
    {
        /* classification */
        ind = get_next_indice(st, 2);
        st->clas_dec = ONSET;
        if( ind == 0 )
        {
            st->clas_dec = UNVOICED_CLAS;
        }
        else if( ind == 1 )
        {
            if( st->last_good >= VOICED_TRANSITION )
            {
                st->clas_dec = VOICED_TRANSITION;
            }
            else
            {
                st->clas_dec = UNVOICED_TRANSITION;
            }
        }
        else if( ind == 2 )
        {
            st->clas_dec = VOICED_CLAS;
        }

        if( st->rf_frame_type == RF_TCXFD )
        {
            /* TCX Gain = 7 bits */
            st->old_gaintcx_bfi = get_next_indice(st, 7);
        }
        else
        {
            /* LTP data */
            if( st->tcxltp )
            {
                if( st->rf_frame_type == RF_TCXTD2 || st->rf_frame_type == RF_TCXTD1 )
                {
                    prm_ltp[0] = 1; /* LTP active*/
                    prm_ltp[1] = get_next_indice(st, 9);
                    prm_ltp[2] = 3; /* max ampl. quantizer output (2bits), anyway not used later*/

                    if(!st->prev_bfi)
                    {
                        st->BER_detect = st->BER_detect |
                                         tcx_ltp_decode_params(&prm_ltp[0], &(st->tcxltp_pitch_int), &(st->tcxltp_pitch_fr), &(st->tcxltp_gain),
                                                               st->pit_min, st->pit_fr1, st->pit_fr2, st->pit_max, st->pit_res_max );

                        st->tcxltp_last_gain_unmodified = st->tcxltp_gain;
                    }
                }
            }
        }
    }


    /*--------------------------------------------------------------------------------*
     * TCX10
     *--------------------------------------------------------------------------------*/

    if( *core == TCX_10_CORE )
    {
        int tcxltp_prm_0 = 0;
        int tcxltp_prm_1 = 0;
        int tcxltp_prm_2 = 0;
        int nbits_igf = 0;
        /* read IGF payload */
        if (st->igf)
        {
            for (k = 0; k < 2; k++)
            {
                n = st->next_bit_pos;

                IGFDecReadLevel( &st->hIGFDec, st, IGF_GRID_LB_SHORT, k == 0 ? 1 : 0 );

                IGFDecReadData( &st->hIGFDec, st, IGF_GRID_LB_SHORT, k == 0 ? 1 : 0 );

                IGFDecStoreTCX10SubFrameData( &st->hIGFDec, k );

                nbits_igf  += st->next_bit_pos - n;
            }
        }

        for( k = 0; k < 2; k++ )
        {
            flag_ctx_hm = 0;

            prm = param + (k*DEC_NPRM_DIV);
            j = 0;
            nbits_tcx = st->next_bit_pos - start_bit_pos;

            if(st->enablePlcWaveadjust && k)
            {
                ind = get_next_indice(st, 1);
                st->tonality_flag = ind;
            }
            /* TCX Gain = 7 bits */
            prm[j++] = get_next_indice(st, 7);

            /* TCX Noise Filling = NBITS_NOISE_FILL_LEVEL bits */
            prm[j++] = get_next_indice(st, NBITS_NOISE_FILL_LEVEL);

            /* LTP data */
            if( (k == 0) && ( st->tcxltp || (st->sr_core > 25600) ) )      /* PLC pitch info for HB */
            {
                prm[j] = get_next_indice(st, 1);

                if( prm[j] )
                {
                    prm[j+1] = get_next_indice(st, 9);
                    prm[j+2] = get_next_indice(st, 2);

                    tcxltp_prm_0 = prm[j];
                    tcxltp_prm_1 = prm[j+1];
                    tcxltp_prm_2 = prm[j+2];
                }

                st->BER_detect = st->BER_detect |
                                 tcx_ltp_decode_params(
                                     &prm[j], &(st->tcxltp_pitch_int), &(st->tcxltp_pitch_fr), &(st->tcxltp_gain),
                                     st->pit_min, st->pit_fr1, st->pit_fr2, st->pit_max, st->pit_res_max );

                st->tcxltp_last_gain_unmodified = st->tcxltp_gain;
                st->tcx_hm_LtpPitchLag = -1;
                j += 3;
            }
            else
            {
                prm[j++] = tcxltp_prm_0;
                prm[j++] = tcxltp_prm_1;
                prm[j++] = tcxltp_prm_2;
            }


            /* TCX spectral data */
            lgFB = st->tcx_cfg.tcx_coded_lines >> 1;;

            if( k==0 && *last_core == ACELP_CORE )
            {
                /* ACE->TCX transition */
                lgFB += lgFB >> 1;
            }

            /* TNS data */
            nTnsParams = 0;
            nTnsBits = 0;

            if( st->tcx_cfg.fIsTNSAllowed )
            {
                if((*last_core == ACELP_CORE) && (k == 0))
                {
                    st->BER_detect = 1;
                    *last_core = 1;
                }
                SetTnsConfig(&st->tcx_cfg, 0, 0 /*(*last_core == ACELP_CORE) && (k == 0)*/); /* lcm : tnsConfig[0][1] does not exist !!! */

                ReadTnsData(st->tcx_cfg.pCurrentTnsConfig, st, &nTnsBits, prm+j, &nTnsParams);

                j += nTnsParams;
            }

            hm_size = (int)(2.0f*st->TcxBandwidth*(float)lgFB);

            /*compute target bits*/
            nbits_tcx = (( st->bits_frame_core - bits_common - nbits_igf + 1 - k)/2) - ( (st->next_bit_pos - start_bit_pos) - nbits_tcx);

            /*Context HM flag*/
            if( st->tcx_cfg.ctx_hm && !(*last_core == ACELP_CORE && k == 0) )
            {
                prm[j] = get_next_indice(st, 1);
                nbits_tcx--;

                if( prm[j] )
                {
                    /* Read PeriodicityIndex */
                    int NumIndexBits = DecodeIndex( st, hm_size >= 256, prm+j+1 );

                    flag_ctx_hm = 1;

                    ConfigureContextHm(
                        lgFB,
                        nbits_tcx,
                        *(prm+j+1),
                        -1,
                        &hm_cfg);

                    nbits_tcx -= NumIndexBits;
                }
            }
            j += NPRM_CTX_HM;

            st->resQBits[k] = ACcontextMapping_decode2_no_mem_s17_LC( st, prm+j,
                              lgFB,
                              nbits_tcx, NPRM_RESQ*st->tcx_cfg.resq, flag_ctx_hm ? &hm_cfg : NULL );
            j += lgFB;

        } /* k, window index */
    }


    if(!st->use_partial_copy)
    {
        if (*total_nbbits-bitsRead[0] < (st->next_bit_pos - start_bit_pos))
        {
            st->BER_detect = 1;
            st->next_bit_pos = start_bit_pos + *total_nbbits - bitsRead[0];
        }
        bitsRead[0] = st->next_bit_pos - start_bit_pos;
    }
    return;
}
