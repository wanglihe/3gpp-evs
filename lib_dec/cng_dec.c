/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"

/*---------------------------------------------------------------------*
 * Local function
 *---------------------------------------------------------------------*/

static void shb_CNG_decod( Decoder_State *st, const float *synth, float *shb_synth, const short sid_bw );


/*---------------------------------------------------------------------*
 * CNG_dec()
 *
 * Decoding of CNG parameters
 *---------------------------------------------------------------------*/

void CNG_dec(
    Decoder_State *st,               /* i/o: State structure                          */
    const short L_frame,           /* i  : length of the frame                      */
    float Aq[],              /* o  : LP coefficients                          */
    const long  core_brate,        /* i  : core bitrate                             */
    float *lsp_new,          /* i/o: current frame LSPs                       */
    float *lsf_new,          /* i/o: current frame LSFs                       */
    short *allow_cn_step     /* o  : allow CN step                            */
    ,short *sid_bw            /* i  : 0-NB/WB, 1-SWB SID                       */
    ,float *q_env
)
{
    float step;
    short i, enr_index;
    short num_bits;
    float enr, weights;
    short m = 0;
    short ptr, j, k;
    short m1;
    float tmp[HO_HIST_SIZE*M];
    short burst_ho_cnt = 0;
    short ll, s_ptr;
    float lsf_tmp[M];
    float C[M];
    float max[2];
    short max_idx[2];
    float ftmp;
    float lsp_tmp[M];
    float dev;
    float max_dev;
    float dist;
    short env_idx[2];
    float enr1;
    float env[NUM_ENV_CNG];
    float tmp_env[HO_HIST_SIZE*NUM_ENV_CNG];
    short LSF_Q_prediction;  /* o  : LSF prediction mode - just temporary variable in CNG                */

    /*-----------------------------------------------------------------*
     * Decode CNG spectral envelope (only in SID frame)
     *-----------------------------------------------------------------*/

    if ( core_brate == SID_1k75 || core_brate == SID_2k40 )
    {
        /* de-quantize the LSF vector */
        if ( st->Opt_AMR_WB )
        {
            isf_dec_amr_wb( st, Aq, lsf_new, lsp_new );
        }
        else
        {
            lsf_dec( st, 0, L_frame, INACTIVE, -1, Aq, &LSF_Q_prediction, lsf_new, lsp_new, 0 );
        }
    }
    else
    {
        /* Use old LSP vector */
        mvr2r( st->lsp_old, lsp_new, M );
        mvr2r( st->lsf_old, lsf_new, M );
    }

    /* Initialize the CNG spectral envelope in case of the very first CNG frame */
    if( st->first_CNG == 0 )

    {
        mvr2r( st->lsp_old, st->lspCNG, M );
    }

    /*-----------------------------------------------------------------*
     * Decode residual signal energy
     *-----------------------------------------------------------------*/

    *allow_cn_step = 0;
    if( core_brate == SID_1k75 || core_brate == SID_2k40 )
    {
        if( core_brate == SID_2k40 )
        {
            step  = STEP_SID;
        }
        else
        {
            step  = STEP_AMR_WB_SID;
        }

        /* intialize the energy quantization parameters */
        if( !st->Opt_AMR_WB )
        {
            num_bits = 7;
        }
        else
        {
            num_bits = 6;
        }

        /* decode the energy index */
        enr_index = (short) get_next_indice( st, num_bits );

        if ( st->last_core_brate > SID_2k40 &&
                st->first_CNG != 0 &&
                st->old_enr_index >= 0 &&
                enr_index > st->old_enr_index + 1
           )
        {
            *allow_cn_step = 1;
        }
        st->old_enr_index = enr_index;
        if ( enr_index == 0 )
        {
            enr_index = -5;
        }
        st->Enew = enr_index / step - 2.0f;

        /* find the new energy value */
        st->Enew = (float)( pow( 2.0, st->Enew ) );

        /* decode SID type */
        if( core_brate == SID_2k40 )
        {
            burst_ho_cnt   = (short) get_next_indice( st, 3 );

            *sid_bw = (short) get_next_indice( st, 1 );
            if ( *sid_bw == 0 )
            {
                env_idx[0] = (short) get_next_indice( st, 6 );

                /* get quantized res_env_details */
                for ( i=0; i<NUM_ENV_CNG; i++ )
                {
                    q_env[i] = CNG_details_codebook[env_idx[0]][i];
                }
            }
        }

        /* Reset CNG history if CNG frame length is changed */
        if ( st->bwidth == WB && st->first_CNG && st->L_frame != st->last_CNG_L_frame )
        {
            st->ho_hist_size = 0;
        }
    }

    /*---------------------------------------------------------------------*
     * CNG spectral envelope update
     * Find A(z) coefficients
     *---------------------------------------------------------------------*/

    if( st->last_core_brate == FRAME_NO_DATA || st->last_core_brate == SID_1k75 || st->last_core_brate == SID_2k40 )
    {
        /* Reset hangover counter if not first SID period */
        if( st->core_brate > FRAME_NO_DATA )
        {
            st->num_ho = 0;
        }

        /* Update LSPs if last SID energy not outlier or insufficient number of hangover frames */
        if( st->num_ho < 3 || st->Enew < 1.5f * st->lp_ener )
        {
            for( i=0; i<M; i++ )
            {
                /* AR low-pass filter  */
                st->lspCNG[i] = CNG_ISF_FACT * st->lspCNG[i] + (1-CNG_ISF_FACT) * lsp_new[i];
            }
        }
    }
    else
    {
        /* Update CNG_mode, if allowed */
        if( ( st->Opt_AMR_WB || st->bwidth == WB ) && ( !st->first_CNG || st->act_cnt2 >= MIN_ACT_CNG_UPD ) )
        {
            if( st->last_active_brate > ACELP_16k40 )
            {
                st->CNG_mode = -1;
            }
            else if( st->last_active_brate > ACELP_13k20 )
            {
                st->CNG_mode = 4;
            }
            else if( st->last_active_brate > ACELP_9k60 )
            {
                st-> CNG_mode = 3;
            }
            else if( st->last_active_brate > ACELP_8k00 )
            {
                st->CNG_mode = 2;
            }
            else if( st->last_active_brate > ACELP_7k20 )
            {
                st->CNG_mode = 1;
            }
            else
            {
                st->CNG_mode = 0;
            }
        }

        /* If first sid after active burst update LSF history from circ buffer */
        burst_ho_cnt = min(burst_ho_cnt, st->ho_circ_size);
        st->act_cnt = 0;
        s_ptr = st->ho_circ_ptr-burst_ho_cnt+1;
        if( s_ptr < 0 )
        {
            s_ptr += st->ho_circ_size;
        }

        for( ll= (short) burst_ho_cnt ; ll > 0 ; ll-- )
        {
            if( ++(st->ho_hist_ptr) == HO_HIST_SIZE )
            {
                st->ho_hist_ptr = 0;
            }
            /* Conversion between 12.8k and 16k LSPs */
            if( ( L_frame == L_FRAME16k && st->ho_16k_lsp[s_ptr] == 0 ) || ( L_frame == L_FRAME && st->ho_16k_lsp[s_ptr] == 1 ) )
            {
                lsp_convert_poly( &st->ho_lsp_circ[s_ptr*M], L_frame, 0 );
            }

            /* update circular buffers */
            mvr2r(&(st->ho_lsp_circ[s_ptr*M]), &(st->ho_lsp_hist[st->ho_hist_ptr*M]), M );
            mvr2r(&(st->ho_ener_circ[s_ptr]), &(st->ho_ener_hist[st->ho_hist_ptr]), 1 );
            st->ho_sid_bw = ( st->ho_sid_bw & 0x3fffffffL ) << 1;
            mvr2r(&(st->ho_env_circ[s_ptr*NUM_ENV_CNG]), &(st->ho_env_hist[st->ho_hist_ptr*NUM_ENV_CNG]), NUM_ENV_CNG );

            st->ho_hist_size++;
            if (st->ho_hist_size > HO_HIST_SIZE)
            {
                st->ho_hist_size = HO_HIST_SIZE;
            }

            s_ptr++;

            if( s_ptr == st->ho_circ_size )
            {
                s_ptr = 0;
            }
        }
        if( st->ho_hist_size > 0)
        {
            *allow_cn_step |= ( st->ho_ener_hist[st->ho_hist_ptr] > 4.0f * st->lp_ener );
        }

        if ( !*allow_cn_step && st->ho_hist_size > 0 )
        {
            ptr = st->ho_hist_ptr;
            mvr2r( &(st->ho_lsp_hist[ptr*M]), tmp, M );
            m1 = 0;
            if( (st->ho_sid_bw & 0x1L) == 0 )
            {
                mvr2r( &st->ho_env_hist[ptr*NUM_ENV_CNG], tmp_env, NUM_ENV_CNG );
                m1 = 1;
            }
            enr = W_DTX_HO[0] * st->ho_ener_hist[ptr];
            weights = W_DTX_HO[0];
            m = 1;
            for( k=1; k<st->ho_hist_size; k++ )
            {
                ptr--;
                if( ptr < 0 )
                {
                    ptr = HO_HIST_SIZE - 1;
                }

                if ( st->ho_ener_hist[ptr] <  st->ho_ener_hist[st->ho_hist_ptr] * BUF_H_NRG &&
                        st->ho_ener_hist[ptr] > st->ho_ener_hist[st->ho_hist_ptr] * BUF_L_NRG )
                {
                    enr += W_DTX_HO[k] * st->ho_ener_hist[ptr];
                    weights += W_DTX_HO[k];
                    mvr2r( &st->ho_lsp_hist[ptr*M], &tmp[m*M], M );
                    if( (st->ho_sid_bw & (0x1L << k)) == 0 )
                    {
                        mvr2r( &st->ho_env_hist[ptr*NUM_ENV_CNG], &tmp_env[m1*NUM_ENV_CNG], NUM_ENV_CNG );
                        m1++;
                    }
                    m++;
                }
            }

            enr /= weights;
            st->lp_ener = enr;

            set_f( max, 0.0f, 2 );
            set_s( max_idx, 0, 2 );

            for( i=0; i<m; i++ )
            {
                if (st->L_frame == L_FRAME )
                {
                    lsp2lsf( &tmp[i*M], lsf_tmp, M, INT_FS_12k8 );
                    ftmp = 6400.0f / (M+1);
                    C[i] = (6400.0f - lsf_tmp[M-1] - ftmp) * (6400.0f - lsf_tmp[M-1] - ftmp);
                }
                else
                {
                    lsp2lsf( &tmp[i*M], lsf_tmp, M, INT_FS_16k );
                    ftmp = 8000.0f / (M+1);
                    C[i] = (8000.0f - lsf_tmp[M-1] - ftmp) * (8000.0f - lsf_tmp[M-1] - ftmp);
                }

                C[i] += (lsf_tmp[0] - ftmp) * (lsf_tmp[0] - ftmp);

                for ( j=0; j<M-1; j++ )
                {
                    C[i] += (lsf_tmp[j+1] - lsf_tmp[j] - ftmp) * (lsf_tmp[j+1] - lsf_tmp[j] - ftmp);
                }

                C[i] *= 0.0588235f; /* 1/M+1 */

                if ( C[i] > max[0] )
                {
                    max[1] = max[0];
                    max_idx[1] = max_idx[0];
                    max[0] = C[i];
                    max_idx[0] = i;
                }
                else if ( C[i] > max[1] )
                {
                    max[1] = C[i];
                    max_idx[1] = i;
                }
            }

            if ( m == 1 )
            {
                mvr2r( tmp, lsp_tmp, M );
            }
            else if ( m < 4 )
            {
                for ( i=0; i<M; i++ )
                {
                    lsp_tmp[i] = 0.0f;
                    for ( j=0; j<m; j++ )
                    {
                        lsp_tmp[i] += tmp[j*M+i];
                    }

                    lsp_tmp[i] -= tmp[max_idx[0]*M+i];
                    lsp_tmp[i] /= (float)(m - 1);
                }
            }
            else
            {
                for ( i=0; i<M; i++ )
                {
                    lsp_tmp[i] = 0.0f;
                    for ( j=0; j<m; j++ )
                    {
                        lsp_tmp[i] += tmp[j*M+i];
                    }

                    lsp_tmp[i] -= (tmp[max_idx[0]*M+i] + tmp[max_idx[1]*M+i]);
                    lsp_tmp[i] /= (float)(m - 2);
                }
            }

            dist = 0.0f;
            max_dev = 0.0f;
            for ( i=0; i<M; i++ )
            {
                dev = (float)fabs(lsp_tmp[i] - lsp_new[i]);
                dist += dev;
                if ( dev > max_dev )
                {
                    max_dev = dev;
                }
            }

            if ( dist > 0.4f || max_dev > 0.1f )
            {
                for( i=0; i<M; i++ )
                {
                    st->lspCNG[i] = lsp_tmp[i];
                }
            }
            else
            {
                for( i=0; i<M; i++ )
                {
                    /* AR low-pass filter  */
                    st->lspCNG[i] = 0.8f * lsp_tmp[i] + (1-0.8f) * lsp_new[i];
                }
            }
            if( m1 > 0 )
            {
                for ( i=0; i<NUM_ENV_CNG; i++ )
                {
                    env[i] = 0;
                    for ( j=0; j<m1; j++ )
                    {
                        env[i] += tmp_env[j*NUM_ENV_CNG+i];
                    }

                    env[i] /= (float)m1;
                    env[i] = env[i] - 2*st->lp_ener;
                }
                mvr2r(env, st->lp_env, NUM_ENV_CNG);
            }
        }
        else
        {
            mvr2r( lsp_new, st->lspCNG, M );  /* use newly analyzed parameters */
        }
    }

    if( core_brate == SID_1k75 || core_brate == SID_2k40 )
    {
        /* Update hangover memory during CNG */
        if ( !*allow_cn_step && (st->Enew < 1.5f * st->lp_ener) )
        {
            /* update the pointer to circular buffer of old LSP vectors */
            if( ++(st->ho_hist_ptr) == HO_HIST_SIZE )
            {
                st->ho_hist_ptr = 0;
            }

            /* update the circular buffer of old LSP vectors with the new LSP vector */
            mvr2r( lsp_new, &(st->ho_lsp_hist[(st->ho_hist_ptr)*M]), M );

            /* update the hangover energy buffer */
            st->ho_ener_hist[st->ho_hist_ptr] = st->Enew;
            if ( core_brate == SID_2k40 && *sid_bw == 0 )
            {
                enr1 = (float)log10( st->Enew*L_frame + 0.1f ) / (float)log10( 2.0f );
                for ( i=0; i<NUM_ENV_CNG; i++ )
                {
                    /* get quantized envelope */
                    env[i] = pow(2.0f,(enr1 - q_env[i])) + 2*st->Enew;
                }
                st->ho_sid_bw = ( st->ho_sid_bw & 0x3fffffffL ) << 1;
                mvr2r( env, &(st->ho_env_hist[(st->ho_hist_ptr)*NUM_ENV_CNG]), NUM_ENV_CNG );
            }
            else if( *sid_bw != 0 )
            {
                st->ho_sid_bw = ( st->ho_sid_bw & 0x3fffffffL ) << 1;
                st->ho_sid_bw |= 0x1L;
            }
            if(++(st->ho_hist_size) > HO_HIST_SIZE)
            {
                st->ho_hist_size = HO_HIST_SIZE;
            }
        }

        /* Update the frame length memory */
        st->last_CNG_L_frame = st->L_frame;

        if( core_brate != SID_1k75 )
        {
            st->num_ho = m;
        }
    }

    if ( st->Opt_AMR_WB )
    {
        isp2a( st->lspCNG, Aq, M );
    }
    else
    {
        lsp2a_stab( st->lspCNG, Aq, M );
    }

    for( i=1; i<L_frame/L_SUBFR; i++ )
    {
        mvr2r( Aq, &Aq[i*(M+1)], M+1 );
    }

    return;
}


/*---------------------------------------------------------------------*
 * swb_CNG_dec()
 *
 * Comfort noise generation for SHB signal
 *---------------------------------------------------------------------*/

void swb_CNG_dec(
    Decoder_State *st,               /* i/o: State structure                          */
    const float *synth,            /* i  : ACELP core synthesis at 32kHz            */
    float *shb_synth,        /* o  : high-band CNG synthesis                  */
    const short sid_bw             /* i  : 0-NB/WB, 1-SWB SID                       */
)
{
    if ( st->core_brate == FRAME_NO_DATA || st->core_brate == SID_2k40 )
    {
        /* SHB SID decoding and CNG */
        if ( st->cng_type == LP_CNG && st->extl == SWB_CNG )
        {
            shb_CNG_decod( st, synth, shb_synth, sid_bw );
        }
        st->last_vad = 0;
        st->burst_cnt = 0;
    }
    else
    {
        st->last_vad = 1;
        st->burst_cnt++;
        if ( st->burst_cnt > 10 )
        {
            st->burst_cnt = 0;
        }
    }

    return;
}

/*---------------------------------------------------------------------*
 * shb_CNG_decod()
 *
 * Main routine of SHB SID decoding and CNG
 *---------------------------------------------------------------------*/

static void shb_CNG_decod(
    Decoder_State *st,               /* i/o: State structure                          */
    const float *synth,            /* i  : ACELP core synthesis at 32kHz            */
    float *shb_synth,         /* o  : high-band CNG synthesis                  */
    const short sid_bw             /* i  : 0-NB/WB, 1-SWB SID                       */
)
{
    short i;
    short idx_ener;
    float shb_lpcCNG[LPC_SHB_ORDER+1];
    float shb_lspCNG[LPC_SHB_ORDER];
    float excTmp[L_FRAME16k];
    float excSHB[L_FRAME16k];
    float ener_excSHB;
    float wb_ener;
    float gain;
    float shb_syn16k[L_FRAME16k];
    float ftmp;
    float step;
    float interp;
    float ener;
    short allow_cn_step=0;

    if( !st->bfi )
    {
        if ( st->core_brate == SID_2k40 && sid_bw == 1 )
        {
            idx_ener = (short) get_next_indice( st, 4 );

            if ( idx_ener == 0 )
            {
                idx_ener = -15;
            }

            /* de-quantization of SHB CNG parameters */
            st->last_shb_cng_ener = ((float)idx_ener/0.9f - 6.0f) * (float)log10(2.0f) * 10.0f;
        }
    }

    /* SHB spectrum estimation */
    interp = ((float)st->shb_dtx_count)/32;
    interp = min( interp, 1.0f );
    for ( i=0; i<LPC_SHB_ORDER; i++ )
    {
        shb_lspCNG[i] = interp * st->lsp_shb_prev[i];
        shb_lspCNG[i] += (1 - interp) * st->lsp_shb_prev_prev[i];
    }

    if( st->shb_dtx_count < 1000 )
    {
        st->shb_dtx_count++;
    }

    lsp2a ( shb_lpcCNG, shb_lspCNG, LPC_SHB_ORDER );
    shb_lpcCNG[0] = 1.0f;

    /* SHB energy estimation */
    wb_ener = 0.001f;
    for ( i=0; i<L_FRAME32k; i++ )
    {
        wb_ener += synth[i] * synth[i];
    }

    wb_ener /= L_FRAME32k;
    wb_ener = 10 * (float)log10(wb_ener);

    if ( st->first_CNG == 0 )
    {
        st->wb_cng_ener = wb_ener;
    }
    if ( fabs(wb_ener - st->wb_cng_ener) > 12.0f )
    {
        allow_cn_step = 1;
    }

    if ( allow_cn_step == 1 )
    {
        st->wb_cng_ener = wb_ener;
    }
    else
    {
        ftmp = wb_ener - st->wb_cng_ener;
        st->wb_cng_ener += 0.9f * ftmp;
    }
    if ( st->core_brate == SID_2k40 && sid_bw == 1 && !st->bfi )
    {
        st->last_wb_cng_ener = st->wb_cng_ener;

        if ( st->first_CNG == 0 )
        {
            st->shb_cng_ener = st->last_shb_cng_ener;
        }
    }

    gain = st->wb_cng_ener - st->last_wb_cng_ener;
    if ( gain > 15 )
    {
        gain = 15;
    }
    step = gain + st->last_shb_cng_ener - st->shb_cng_ener;
    if ( allow_cn_step == 1 || st->last_core_brate > SID_2k40 )
    {
        st->shb_cng_ener += step;
    }
    else
    {
        st->shb_cng_ener += 0.25f * step;
    }
    /* generate white noise excitation */
    for ( i=0; i<L_FRAME16k; i++ )
    {
        excTmp[i] = (float)own_random( &st->swb_cng_seed );
    }

    /* synthesis filtering */
    syn_filt( shb_lpcCNG, LPC_SHB_ORDER, excTmp, excSHB, L_FRAME16k, st->state_lpc_syn , 1 );

    /* synthesis signal gain shaping */
    ener_excSHB = 0.001f;
    for ( i=0; i<L_FRAME16k; i++ )
    {
        ener_excSHB += excSHB[i] * excSHB[i];
    }

    if ( st->last_vad == 1 )
    {
        if ( st->burst_cnt > 3 && st->last_core != HQ_CORE )
        {
            st->trans_cnt = 5;
        }
        else
        {
            st->trans_cnt = 0;
        }
    }

    if ( st->trans_cnt > 0 )
    {
        i = (short)((float)st->trans_cnt / 15.0f * 255);
        ener = st->shb_cng_ener + sin_table256[i] * (st->last_shb_ener - st->shb_cng_ener);
        st->trans_cnt--;
    }
    else
    {
        ener = st->shb_cng_ener;
    }

    gain = (float)sqrt( pow(10, 0.1f*ener)*L_FRAME16k/ener_excSHB );

    for ( i=0; i<L_FRAME16k; i++ )
    {
        shb_syn16k[i] = gain * excSHB[i];
    }

    /* generate 32kHz SHB signal (12.8 - 14.4kHz) from 12.8kHz signal */
    GenSHBSynth( shb_syn16k, shb_synth, st->genSHBsynth_Hilbert_Mem, st->genSHBsynth_state_lsyn_filt_shb_local, st->L_frame, &(st->syn_dm_phase) );

    if ( st->output_Fs == 48000 )
    {
        interpolate_3_over_2_allpass( shb_synth, L_FRAME32k, shb_synth, st->interpol_3_2_cng_dec, allpass_poles_3_ov_2 );
    }
    ResetSHBbuffer_Dec( st );
    return;
}
