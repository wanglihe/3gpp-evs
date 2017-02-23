/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_enc.h"
#include "rom_com.h"
#include "prot.h"

/*---------------------------------------------------------------------*
 * Local constants
 *---------------------------------------------------------------------*/

#define MAX_DELTA_CNG             1

/*---------------------------------------------------------------------*
 * Local functions
 *---------------------------------------------------------------------*/

static short shb_DTX( Encoder_State *st, const float *shb_speech, const float *syn_12k8_16k );

static void shb_CNG_encod( Encoder_State *st, const short update
                         );

/*---------------------------------------------------------------------*
 * CNG_enc()
 *
 * Confort noise generation for the coder
 *---------------------------------------------------------------------*/

void CNG_enc(
    Encoder_State *st,             /* i/o: State structure                                 */
    const short L_frame,         /* i  : length of the frame                             */
    float Aq[],            /* o  : LP coefficients                                 */
    const float *speech,         /* i  : pointer to current frame input speech buffer    */
    float enr,             /* i  : residual energy from Levinson-Durbin            */
    float *lsp_new,        /* i/o: current frame ISPs                              */
    float *lsf_new,        /* i/o: current frame ISFs                              */
    short *allow_cn_step,  /* o  : allow CN step                                   */
    short burst_ho_cnt,    /* i  : hangover frames at end of speech burst          */
    float *q_env,
    short *sid_bw,
    float *exc_mem2
)
{
    short enr_index, i;
    float step, res[L_FRAME16k];
    short maxl, num_bits;
    short j, k, ptr;
    short m1;
    float weights;
    float sp_enr;
    short m = 0;
    float tmp[HO_HIST_SIZE*M];
    short ll, s_ptr;
    float att=1.0f;
    float lsf_tmp[M];
    float C[M];
    float max[2];
    short max_idx[2];
    float ftmp;
    float lsp_tmp[M];
    float dev;
    float max_dev;
    float dist;
    short max_idx1[2]= {0,0};
    float fft_io[L_FRAME16k];
    float *ptR,*ptI;
    float enr1=0;
    float env[NUM_ENV_CNG];
    float min1;
    short min1_idx;
    float d;
    float res1[L_FRAME16k];
    float tmp_env[HO_HIST_SIZE*NUM_ENV_CNG];
    short force_cn_step=0;

    /* calculate input energy */
    sp_enr = (float) log10( sum2_f( speech, L_frame )/L_frame + 0.1f )/ (float)log10(2.0f);

    if (sp_enr < 0.0f)
    {
        sp_enr = 0.0f;
    }

    if ( st->first_CNG == 0 || st->old_enr_index < 0 )
    {
        st->lp_sp_enr = sp_enr;
    }
    else
    {
        if ( st->last_core_brate > SID_2k40 && burst_ho_cnt > 0 && st->lp_sp_enr < 6.0f && (sp_enr - st->lp_sp_enr) > 4.0f && sp_enr > 6.0f )
        {
            st->lp_sp_enr = sp_enr;
            force_cn_step = 1;
        }
        else
        {
            st->lp_sp_enr  = 0.1f * sp_enr  + 0.9f * st->lp_sp_enr;
        }
    }

    /* update the pointer to circular buffer of old LSP vectors */
    if( ++(st->cng_hist_ptr) == DTX_HIST_SIZE )
    {
        st->cng_hist_ptr = 0;
    }

    /* update the circular buffer of old LSP vectors with the new LSP vector */
    mvr2r( lsp_new, &(st->cng_lsp_hist[(st->cng_hist_ptr)*M]), M );

    /*-----------------------------------------------------------------*
     * Find CNG spectral envelope
     * Find LSP median
     *-----------------------------------------------------------------*/

    if( (st->core_brate == SID_2k40 || st->core_brate == SID_1k75) && st->cng_cnt >= (st->cng_hist_size-1) )
    {
        set_f( max, 0.0f, 2 );
        set_s( max_idx, 0, 2 );

        for( i=0; i<st->cng_hist_size; i++ )
        {
            if (st->L_frame == L_FRAME )
            {
                lsp2lsf( &st->cng_lsp_hist[i*M], lsf_tmp, M, INT_FS_12k8 );
                ftmp = 6400.0f / (M+1);
                C[i] = (6400.0f - lsf_tmp[M-1] - ftmp) * (6400.0f - lsf_tmp[M-1] - ftmp);
            }
            else
            {
                lsp2lsf( &st->cng_lsp_hist[i*M], lsf_tmp, M, INT_FS_16k );
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

        for ( i=0; i<M; i++ )
        {
            lsp_new[i] = 0.0f;
            for ( j=0; j<st->cng_hist_size; j++ )
            {
                lsp_new[i] += st->cng_lsp_hist[j*M+i];
            }

            lsp_new[i] -= (st->cng_lsp_hist[max_idx[0]*M+i] + st->cng_lsp_hist[max_idx[1]*M+i]);
            lsp_new[i] /= (float)(st->cng_hist_size - 2);
        }
        max_idx1[0] = max_idx[0];
        max_idx1[1] = max_idx[1];
    }

    /*-----------------------------------------------------------------*
     * Quantize CNG spectral envelope (only in SID frame)
     * Quantize the LSF vector
     *-----------------------------------------------------------------*/

    *allow_cn_step = ((st->cng_cnt == 0) &&
                      (st->lp_sp_enr > 6.0f) &&
                      ((st->lp_sp_enr + 4.0f) < sp_enr) &&
                      (st->first_CNG != 0 ) &&
                      (st->old_enr_index >= 0) &&
                      (st->last_core_brate > SID_2k40)) ||
                     force_cn_step;

    if( st->core_brate == SID_2k40 || st->core_brate == SID_1k75 )
    {
        /* LSF quantization */
        if ( st->Opt_AMR_WB )
        {
            isf_enc_amr_wb( st, lsf_new, lsp_new, 0, 0 );
        }
        else
        {
            lsf_enc( st, L_frame, INACTIVE, lsf_new, lsp_new, 0, 0, 0, 100 );
        }

        /* Reset CNG history if CNG frame length is changed */
        if ( st->bwidth == WB && st->first_CNG && st->L_frame != st->last_CNG_L_frame )
        {
            st->ho_hist_size = 0;
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
        /* Update CNG_mode if allowed */
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
                st->CNG_mode = 3;
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

        /* If first SID after active burst update LSF history from circ buffer */
        st->act_cnt = 0;
        s_ptr = st->ho_circ_ptr-burst_ho_cnt+1;
        if( s_ptr < 0 )
        {
            s_ptr += st->ho_circ_size;
        }

        for( ll = burst_ho_cnt; ll > 0; ll-- )
        {
            if( ++(st->ho_hist_ptr) == HO_HIST_SIZE )
            {
                st->ho_hist_ptr = 0;
            }
            /* Conversion between 12.8k and 16k LSPs */
            if( L_frame == L_FRAME && st->ho_16k_lsp[s_ptr] == 1 )
            {
                /* Conversion from 16k LPSs to 12k8 */
                lsp_convert_poly( &(st->ho_lsp_circ[s_ptr*M]), L_frame, 0 );
            }
            else if ( L_frame == L_FRAME16k && st->ho_16k_lsp[s_ptr] == 0 )
            {
                /* 16k LSPs already converted and stored, just copy to the other buffer */
                mvr2r(&(st->ho_lsp_circ2[s_ptr*M]), &(st->ho_lsp_circ[s_ptr*M]), M );
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
        if ( burst_ho_cnt > 0)
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
                mvr2r(tmp, lsp_tmp, M);
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

    /*-----------------------------------------------------------------*
     * Find residual signal
     * Calculate residual signal energy per sample
     *-----------------------------------------------------------------*/

    /* calculate the residual signal */
    residu( Aq, M, speech, res, L_frame );

    mvr2r(res, res1, L_frame);
    if( st->bwidth != NB )
    {
        if( st->bwidth == WB && st->CNG_mode >= 0 )
        {
            ftmp = HO_ATT[st->CNG_mode];
        }
        else
        {
            ftmp = 0.6f;
        }

        att = ftmp/6.0f;
        att = 1.0f/(1 + att * 8);

        if ( att < ftmp )
        {
            att = ftmp;
        }

        for( i = 0; i < st->L_frame; i++ )
        {
            res1[i] *= att;
        }
    }

    /* calculate the spectrum of residual signal */
    mvr2r(res1, fft_io, st->L_frame);

    if ( st->L_frame == L_FRAME16k )
    {
        modify_Fs( fft_io, L_FRAME16k, 16000, fft_io, 12800, exc_mem2, 0 );
    }

    fft_rel(fft_io, L_FFT, LOG2_L_FFT);
    ptR = &fft_io[1];
    ptI = &fft_io[L_FFT-1];
    for ( i=0; i<NUM_ENV_CNG; i++ )
    {
        env[i] = 2.0f*(*ptR **ptR + *ptI **ptI)/L_FFT;
        ptR++;
        ptI--;
    }

    mvr2r( env, &(st->cng_res_env[(st->cng_hist_ptr)*NUM_ENV_CNG]), NUM_ENV_CNG );
    /* calculate the residual signal energy */
    enr = dotp( res, res, L_frame ) / L_frame;

    /* convert log2 of residual signal energy */
    enr = (float)log10( enr + 0.1f ) / (float)log10( 2.0f );

    /* update the circular buffer of old energies */
    st->cng_ener_hist[st->cng_hist_ptr] = enr;

    /*-----------------------------------------------------------------*
     * Quantize residual signal energy (only in SID frame)
     *-----------------------------------------------------------------*/

    if( st->core_brate == SID_2k40 || st->core_brate == SID_1k75 )
    {
        if( st->cng_cnt >= st->cng_hist_size - 1 )
        {
            /* average the envelope except outliers */
            for ( i=0; i<NUM_ENV_CNG; i++ )
            {
                for ( j=0; j<st->cng_hist_size; j++ )
                {
                    env[i] += st->cng_res_env[j*NUM_ENV_CNG+i];
                }

                env[i] -= (st->cng_res_env[max_idx1[0]*NUM_ENV_CNG+i] + st->cng_res_env[max_idx1[1]*NUM_ENV_CNG+i]);
                env[i] /= (float)(st->cng_hist_size - 2);
            }
            /* compute average excitation energy */
            enr = 0;
            weights = 0;
            ptr = st->cng_hist_ptr;

            for( k=0; k<st->cng_hist_size; k++ )
            {
                enr += W_HIST[k] * st->cng_ener_hist[ptr--];
                if( ptr < 0 )
                {
                    ptr = DTX_HIST_SIZE - 1;
                }

                weights += W_HIST[k];
            }

            /* normalize the average value */
            enr /= weights;
        }

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

        /* intialize the energy quantization parameters */
        if( !st->Opt_AMR_WB )
        {
            step = STEP_SID;
            maxl = 127;
            num_bits = 7;
        }
        else
        {
            step = STEP_AMR_WB_SID;
            maxl = 63;
            num_bits = 6;
        }

        /* calculate the energy quantization index */
        enr_index = (short)( (enr + 2.0f) * step );

        /* limit the energy quantization index */
        if( enr_index > maxl )
        {
            enr_index = maxl;
        }

        if( enr_index < 0 )
        {
            enr_index = 0;
        }

        /* allow only slow energy increase */
        if( st->first_CNG && enr_index > st->old_enr_index + MAX_DELTA_CNG )
        {
            if( *allow_cn_step == 1 )
            {
                enr_index = st->old_enr_index + (short) (0.85f*(enr_index - st->old_enr_index)) ;
            }
            else
            {
                enr_index = st->old_enr_index + MAX_DELTA_CNG;
            }
        }
        st->old_enr_index = enr_index;

        push_indice( st, IND_ENERGY, enr_index, num_bits );
        if ( enr_index == 0 )
        {
            enr_index = -5;
        }
        /* find the quatized energy */
        st->Enew = (float)enr_index / step - 2.0f;
        st->Enew = (float)( pow( 2.0f, st->Enew ) );
        if ( st->core_brate == SID_2k40 )
        {
            enr1 = (float)log10( st->Enew*L_frame + 0.1f ) / (float)log10( 2.0f );
            for ( i=0; i<NUM_ENV_CNG; i++ )
            {
                env[i] -= 2 * st->Enew;

                if ( env[i] < 0.0f )
                {
                    env[i] = 0.1f;
                }

                env[i] = (float)log10( env[i] + 0.1f ) / (float)log10( 2.0f );
                env[i] -= att;

                if ( env[i] < 0 )
                {
                    env[i] = 0;
                }

                env[i] = enr1 - env[i];
            }

            /* codebook search */
            min1 = 9999.0f;;
            min1_idx = 0;

            for ( i=0; i<64; i++ )
            {
                d = 0.0f;
                for ( j=0; j<NUM_ENV_CNG; j++ )
                {
                    d += (env[j] - CNG_details_codebook[i][j]) * (env[j] - CNG_details_codebook[i][j]);
                }

                if ( d < min1 )
                {
                    min1 = d;
                    min1_idx = i;
                }
            }
            push_indice( st, IND_CNG_ENV1, min1_idx, 6 );
            /* get quantized res_env_details */
            for ( i=0; i<NUM_ENV_CNG; i++ )
            {
                q_env[i] = CNG_details_codebook[min1_idx][i];
            }
        }
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
            if ( st->core_brate == SID_2k40 )
            {
                for ( i=0; i<NUM_ENV_CNG; i++ )
                {
                    /* get quantized envelope */
                    env[i] = pow(2.0f,(enr1 - q_env[i])) + 2*st->Enew;
                }
                mvr2r( env, &(st->ho_env_hist[(st->ho_hist_ptr)*NUM_ENV_CNG]), NUM_ENV_CNG );
            }
            if(++(st->ho_hist_size) > HO_HIST_SIZE)
            {
                st->ho_hist_size = HO_HIST_SIZE;
            }
        }
    }

    /* dithering bit for AMR-WB IO mode is always set to 0 */
    if( st->core_brate == SID_1k75 )
    {
        push_indice( st, IND_DITHERING, 0, 1 );
    }

    if ( st->core_brate == SID_2k40 )
    {
        push_indice( st, IND_ACELP_16KHZ, st->L_frame == L_FRAME16k ? 1 : 0, 1 );
    }

    if ( st->core_brate == SID_2k40 )
    {
        /*  transmit ho_cnt for use at decoder side as  CNG synthesis assistance   */
        if( st->burst_ho_cnt > (HO_HIST_SIZE-1) )
        {
            push_indice( st, IND_CNG_HO, (HO_HIST_SIZE-1), 3 );  /* send max allowed value , limited to 7 */
        }
        else
        {
            push_indice( st, IND_CNG_HO, st->burst_ho_cnt, 3 );   /* send actual value */
        }
        st->num_ho = m;
        push_indice( st, IND_SID_TYPE, 0, 1 );

        if ( st->input_Fs < 32000 )
        {
            push_indice( st, IND_SID_BW, 0, 1 );
            *sid_bw = 0;
        }
    }

    /*-----------------------------------------------------------------*
     * Updates
     *-----------------------------------------------------------------*/

    /* update the SID frames counter */
    if( st->core_brate == SID_2k40 || st->core_brate == SID_1k75 )
    {
        st->cng_cnt = 0;
        st->cng_hist_ptr = -1;

        /* update frame length memory */
        st->last_CNG_L_frame = st->L_frame;
    }
    else
    {
        st->cng_cnt++;
    }



    return;
}

/*---------------------------------------------------------------------*
 * swb_CNG_enc()
 *
 * SWB DTX/CNG encoding
 *---------------------------------------------------------------------*/

void swb_CNG_enc(
    Encoder_State *st,             /* i/o: State structure                                 */
    const float *shb_speech,     /* i  : SHB target signal (6-14kHz) at 16kHz            */
    const float *syn_12k8_16k    /* i  : ACELP core synthesis at 12.8kHz or 16kHz        */
)
{
    short shb_SID_updt;

    if( st->core_brate == SID_2k40 || st->core_brate == FRAME_NO_DATA )
    {
        if ( st->cng_type == LP_CNG )
        {
            /* decide if SHB SID encoding or not */
            shb_SID_updt = shb_DTX( st, shb_speech, syn_12k8_16k );

            /* SHB CNG encoding */
            shb_CNG_encod( st, shb_SID_updt );
        }
        st->last_vad = 0;
    }
    else
    {
        st->last_vad = 1;
    }

    return;
}

/*---------------------------------------------------------------------*
 * shb_CNG_encod()
 *
 * SID parameters encoding for SHB signal
 *---------------------------------------------------------------------*/

static void shb_CNG_encod(
    Encoder_State *st,             /* i/o: State structure                                 */
    const short update           /* i  : SID update flag                                 */
)
{
    short idx_ener = 0;

    if ( update == 1 )
    {
        /* SHB energy quantization */
        idx_ener = (short)(0.9f * (0.1f*st->mov_shb_cng_ener/(float)log10(2.0f) + 6.0f) + 0.5f );
        if ( st->bwidth < SWB )
        {
            idx_ener = 0;
        }

        if ( idx_ener > 15 )
        {
            idx_ener = 15;
        }
        else if ( idx_ener < 0 )
        {
            idx_ener = 0;
        }

        push_indice( st, IND_SHB_CNG_GAIN, idx_ener, 4 );
        push_indice( st, IND_SID_BW, 1, 1 );
        st->nb_bits_tot = st->nb_bits_tot - st->ind_list[IND_CNG_ENV1].nb_bits;
        st->ind_list[IND_CNG_ENV1].nb_bits = -1;
        push_indice( st, IND_UNUSED, 0, 2 );
        st->ho_sid_bw = ( st->ho_sid_bw & 0x3fffffffL ) << 1;
        st->ho_sid_bw |= 0x1L;
    }
    else if ( st->core_brate == SID_2k40 )
    {
        st->ho_sid_bw = ( st->ho_sid_bw & 0x3fffffffL ) << 1;
        push_indice( st, IND_SID_BW, 0, 1 );
    }


    return;
}

/*---------------------------------------------------------------------*
 * shb_DTX()
 *
 * Decide if encoding SHB SID or not
 *---------------------------------------------------------------------*/

static short shb_DTX(
    Encoder_State *st,             /* i/o: State structure                                 */
    const float *shb_speech,     /* i  : SHB target signal (6-14kHz) at 16kHz            */
    const float *syn_12k8_16k    /* i  : ACELP core synthesis at 12.8kHz or 16kHz        */
)
{
    short i;
    short update;
    float shb_old_speech[(L_LOOK_12k8 + L_SUBFR + L_FRAME) * 5/4];
    float *shb_new_speech;
    float wb_ener;
    float shb_ener;
    float log_wb_ener;
    float log_shb_ener;
    float ftmp;
    short allow_cn_step=0;

    shb_new_speech = shb_old_speech + (L_LOOK_12k8 + L_SUBFR) * 5/4;
    mvr2r( st->old_speech_shb, shb_old_speech, (L_LOOK_12k8 + L_SUBFR) * 5/4 );
    mvr2r( shb_speech, shb_new_speech, L_FRAME16k );
    mvr2r( shb_old_speech + L_FRAME16k, st->old_speech_shb, (L_LOOK_12k8 + L_SUBFR) * 5/4 );

    shb_ener = 0;
    for ( i=0; i<L_FRAME16k; i++ )
    {
        shb_ener += shb_old_speech[i] * shb_old_speech[i];
    }
    shb_ener /= L_FRAME16k;

    wb_ener = sum2_f( syn_12k8_16k, st->L_frame) + 0.001f;
    wb_ener = wb_ener/st->L_frame;

    log_wb_ener = 10 * (float)log10(wb_ener);
    log_shb_ener = 10 * (float)log10(shb_ener) - 6.5f;

    if ( st->first_CNG == 0 )
    {
        st->mov_wb_cng_ener = log_wb_ener;
        st->mov_shb_cng_ener = log_shb_ener;
        st->last_wb_cng_ener = log_wb_ener;
        st->last_shb_cng_ener = log_shb_ener;
    }
    if ( fabs(log_wb_ener - st->mov_wb_cng_ener) > 12.0f )
    {
        allow_cn_step = 1;
    }

    if ( allow_cn_step == 1 )
    {
        st->mov_wb_cng_ener = log_wb_ener;
        st->mov_shb_cng_ener = log_shb_ener;
    }
    else
    {
        ftmp = log_wb_ener - st->mov_wb_cng_ener;

        st->mov_wb_cng_ener += 0.9f * ftmp;

        ftmp = log_shb_ener - st->mov_shb_cng_ener;

        st->mov_shb_cng_ener += 0.25f * ftmp;
    }
    st->shb_NO_DATA_cnt++;

    update = 0;
    if ( st->core_brate == SID_2k40 )
    {
        if ( st->first_CNG == 0 )
        {
            update = 1;
        }
        else if ( st->shb_cng_ini_cnt > 0 )
        {
            st->shb_cng_ini_cnt--;
            update = 1;
        }
        else if ( st->last_vad == 1 )
        {
            update = 1;
        }
        else if ( st->shb_NO_DATA_cnt >= 100 )
        {
            update = 1;
        }
        else if ( fabs((st->mov_wb_cng_ener - st->mov_shb_cng_ener) - (st->last_wb_cng_ener - st->last_shb_cng_ener)) > 3.0f )
        {
            update = 1;
        }
        else if ( (st->bwidth >=SWB && st->last_SID_bwidth < SWB) || (st->bwidth <SWB && st->last_SID_bwidth >= SWB) )
        {
            update = 1;
        }

        st->last_SID_bwidth = st->bwidth;
    }

    if ( update == 1 )
    {
        st->last_wb_cng_ener = st->mov_wb_cng_ener;
        st->last_shb_cng_ener = st->mov_shb_cng_ener;
        st->shb_NO_DATA_cnt = 0;
    }

    return (update);
}
