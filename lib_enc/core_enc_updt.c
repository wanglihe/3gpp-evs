/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "options.h"
#include "prot.h"
#include "rom_com.h"
#include "cnst.h"


/*-------------------------------------------------------------------*
* core_encode_update()
*
* Common updates of buffers
*-------------------------------------------------------------------*/

void core_encode_update(
    Encoder_State *st           /* i/o: Encoder state structure     */
)
{
    short n;

    /* Update Input Signal Buffers */
    n = st->encoderPastSamples_enc + st->encoderLookahead_enc;

    mvr2r( st->buf_speech_enc_pe + st->L_frame, st->buf_speech_enc_pe, n );
    mvr2r( st->buf_speech_enc + st->L_frame, st->buf_speech_enc, n );

    if( !st->tcxonly )
    {
        n = st->L_frame + st->L_frame/4;
        mvr2r( st->buf_wspeech_enc + st->L_frame, st->buf_wspeech_enc, n );
    }

    if( st->core == ACELP_CORE || st->core == AMR_WB_CORE || st->core_brate == SID_2k40 || st->core_brate == FRAME_NO_DATA )
    {
        mvr2r( st->buf_speech_enc + st->L_frame, st->buf_speech_ltp + st->L_frame, st->L_frame );
    }

    n = st->encoderPastSamples_enc + st->encoderLookahead_enc;
    mvr2r( st->buf_speech_ltp + st->L_frame, st->buf_speech_ltp, n );
    mvr2r( st->buf_synth + st->L_frame, st->buf_synth, st->L_frame+L_SUBFR );

    if( (st->core_brate <= SID_2k40 && st->cng_type == FD_CNG) || (st->tcxonly && st->codec_mode == MODE2) )
    {
        /* reset LP memories */
        set_zero( st->mem_MA, M );
        mvr2r( GEWB_Ave, st->mem_AR, M );
    }


    return;
}

/*-------------------------------------------------------------------*
* core_encode_update_cng()
*
* Common updates in case of CNG
*-------------------------------------------------------------------*/

void core_encode_update_cng(
    Encoder_State *st,
    float *timeDomainBuffer,
    float *A,
    const float Aw[]      /* i  : weighted A(z) unquant. for subframes*/
)
{
    float lsp[M], lsf[M];
    float *synth, synth_buf[M+1+L_FRAME_PLUS+L_FRAME_PLUS/2], wsyn[L_FRAME_PLUS];
    float tmp;
    float enr;
    float att;
    float enr_index;

    /* LPC -> LSP/lsp */
    a2lsp_stab( A, lsp, st->lsp_old );

    /* LSP/lsp -> LSF/lsf */
    if( st->L_frame == L_FRAME16k )
    {
        lsp2lsf( lsp, lsf, M, INT_FS_16k );
    }
    else
    {
        lsp2lsf( lsp, lsf, M, INT_FS_12k8 );
    }

    /* Update synth memory */
    synth = synth_buf + (1+M);
    mvr2r( st->LPDmem.syn, synth_buf, 1+M );
    mvr2r( timeDomainBuffer, synth, st->L_frame);
    mvr2r( synth+st->L_frame-(1+M), st->LPDmem.syn, 1+M );
    mvr2r( synth, st->synth, st->L_frame);

    /* Update ZIR */
    set_zero( synth+st->L_frame, st->L_frame/2 );
    syn_filt( A, M,synth+st->L_frame, synth+st->L_frame, st->L_frame/2, &synth[st->L_frame-M], 0 );
    mvr2r( synth+st->L_frame-(st->L_frame/2), st->LPDmem.Txnq, st->L_frame/2 );

    /* Update pe-synth memory */
    tmp = synth[-(1+M)];
    preemph( synth-M, st->preemph_fac, M+st->L_frame, &tmp );
    mvr2r( synth+st->L_frame-M, st->LPDmem.mem_syn, M );
    mvr2r( synth+st->L_frame-M, st->LPDmem.mem_syn2, M );

    /* Update excitation memory */
    mvr2r( st->LPDmem.old_exc+st->L_frame, st->LPDmem.old_exc, max(L_EXC_MEM-st->L_frame,0));
    residu( A, M,synth, st->LPDmem.old_exc+max(L_EXC_MEM-st->L_frame,0), st->L_frame );
    if ( st->core_brate == SID_2k40 )
    {
        enr = dotp( st->LPDmem.old_exc+max(L_EXC_MEM-st->L_frame,0), st->LPDmem.old_exc+max(L_EXC_MEM-st->L_frame,0), st->L_frame ) / st->L_frame;
        enr = (float)log10( enr + 0.1f ) / (float)log10( 2.0f );

        /* decrease the energy in case of WB input */
        if( st->bwidth != NB )
        {
            if( st->bwidth == WB )
            {
                if( st->CNG_mode >= 0 )
                {
                    /* Bitrate adapted attenuation */
                    att = ENR_ATT[st->CNG_mode];
                }
                else
                {
                    /* Use least attenuation for higher bitrates */
                    att = ENR_ATT[4];
                }
            }
            else
            {
                att = 1.5f;
            }

            enr -= att;
        }

        enr_index = (short)( (enr + 2.0f) * STEP_SID );
        if( enr_index > 127 )
        {
            enr_index = 127;
        }

        if( enr_index < 0 )
        {
            enr_index = 0;
        }
        st->old_enr_index = enr_index;
    }

    /* Update weighted synthesis memory */
    calc_residu( synth, wsyn, Aw, st->L_frame );
    tmp = st->wspeech_enc[-1]-st->LPDmem.mem_w0;
    deemph( wsyn, st->preemph_fac, st->L_frame, &tmp );
    st->LPDmem.mem_w0 = st->wspeech_enc[st->L_frame-1]-wsyn[st->L_frame-1];

    /* Update LPC-related memories */
    mvr2r( lsp, st->lsp_old, M );
    mvr2r( lsf, st->lsf_old, M );
    st->envWeighted = 0;
    mvr2r( A, st->old_Aq_12_8, M+1 );
    st->old_Es_pred = 0;

    /* Reset acelp memories */
    set_zero( st->dispMem, 8 );
    st->LPDmem.tilt_code = TILT_CODE;
    st->LPDmem.gc_threshold = 0.0f;

    /* Update ace/tcx mode */
    st->core = ACELP_CORE;

    /* Reset TCX overlap */
    st->tcx_cfg.tcx_curr_overlap_mode = st->tcx_cfg.tcx_last_overlap_mode = ALDO_WINDOW;

    if( st->first_CNG == 0 )
    {
        mvr2r( st->lsp_old, st->lspCNG, M );
    }

    return;
}
