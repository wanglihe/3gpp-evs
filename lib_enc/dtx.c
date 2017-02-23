/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <assert.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"

/*-------------------------------------------------------------------*
 * Local constants
 *-------------------------------------------------------------------*/

#define ALPHA_ENER_SLOW   0.99f   /* Slow adaptation (noise up, speech down) */
#define ALPHA_ENER_FAST   0.90f   /* Fast adaptation (noise down, speech up) */
#define MIN_CNT           50      /* Minimum frame number before SID interval adaptation */

#define SNR_H             51.0f   /* Estimated SNR and corresponding SID interval        */
/* 51dB corresponds to 25dB SNR before noise supressor */
#define SNR_L             36.0f
#define INT_H             50
#define INT_L             8

#define LTE_VAR           -4.0f

#define CNG_TYPE_HO       20      /* hangover for switching between CNG types */


/*-------------------------------------------------------------------*
 * Local functions
 *-------------------------------------------------------------------*/

static void update_SID_cnt( Encoder_State *st );


/*-------------------------------------------------------------------*
 * dtx()
 *
 * Discontinuous transmission operation
 *-------------------------------------------------------------------*/

void dtx(
    Encoder_State *st,                  /* i/o: encoder state structure                  */
    const short vad,                  /* i  : vad flag                                 */
    const float speech[]              /* i  : Pointer to the speech frame              */
)
{
    float alpha;

    /* Initialization */
    if( st->ini_frame == 0 )
    {
        st->active_cnt = CNG_TYPE_HO;
        if( st->codec_mode == MODE1 || st->Opt_AMR_WB )
        {
            st->cng_type = LP_CNG;
        }
        else
        {
            st->cng_type = FD_CNG;
        }
    }

    if( st->Opt_DTX_ON && vad == 0 &&
            st->ini_frame > 2 &&                     /* CNG coding starts after 3 frames */
            st->fd_cng_reset_flag == 0 &&
            st->last_core != AMR_WB_CORE &&
            st->Opt_AMR_WB == 0 )
    {
        if ( st->last_core_brate > SID_2k40 &&
                st->last_total_brate_cng != -1 &&
                st->last_total_brate_cng != st->total_brate &&
                (st->last_total_brate_cng <=  ACELP_24k40 || st->lp_noise < 15)
           )
        {
            st->total_brate = st->last_total_brate_cng;
            if( !(st->total_brate == ACELP_7k20 && st->Opt_SC_VBR) )
            {
                st->Opt_SC_VBR = 0;
            }

            st->Opt_RF_ON = 0;
            if( st->rf_mode && st->rf_fec_offset > 0 && st->total_brate == ACELP_13k20 && st->bwidth != NB )
            {
                st->Opt_RF_ON = 1;
            }
            st->rf_mode = st->Opt_RF_ON;
            st->bwidth = st->last_bwidth_cng;
            st->codec_mode = st->last_codec_mode_cng;
        }
        if ( st->last_core_brate <= SID_2k40 &&
                st->last_total_brate != st->total_brate &&
                ( st->last_total_brate <= ACELP_24k40 || st->lp_noise < 15 ) )
        {
            st->total_brate = st->last_total_brate;
            if( !(st->total_brate == ACELP_7k20 && st->Opt_SC_VBR) )
            {
                st->Opt_SC_VBR = 0;
            }

            st->Opt_RF_ON = 0;
            if( st->rf_mode && st->rf_fec_offset > 0 && st->total_brate == ACELP_13k20 && st->bwidth != NB )
            {
                st->Opt_RF_ON = 1;
            }
            st->rf_mode = st->Opt_RF_ON;
            st->bwidth = st->last_bwidth;
            switch ( st->total_brate )
            {
            case 5900:
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
            }
        }
    }

    /*------------------------------------------------------------------------*
     * Select SID or FRAME_NO_DATA frame if DTX is enabled
     *------------------------------------------------------------------------*/

    if( st->Opt_DTX_ON && vad == 0 &&
            st->ini_frame > 2 &&                     /* CNG coding starts after 3 frames */
            ( st->total_brate <= ACELP_24k40 || st->lp_noise < 15 ) && /* at higher bitrates, DTX kicks in only when the level of background noise is low */
            st->fd_cng_reset_flag == 0 )
    {
        /* reset counter */
        st->active_cnt = 0;

        if ( st->Opt_AMR_WB )
        {
            st->last_total_brate_cng = -1;
        }
        else
        {
            st->last_total_brate_cng = st->total_brate;
            st->last_bwidth_cng = st->bwidth;
            st->last_codec_mode_cng = st->codec_mode;
        }

        if( st->cnt_SID == 0 )
        {
            /* this will be a SID frame */
            if ( st->Opt_AMR_WB )
            {
                st->core_brate = SID_1k75;
            }
            else
            {
                st->core_brate = SID_2k40;
            }
        }
        else
        {
            /* this will be a no data frame */
            st->core_brate = FRAME_NO_DATA;
        }

        if( st->core_brate == FRAME_NO_DATA && st->last_core != ACELP_CORE && !st->Opt_AMR_WB )
        {
            /* force SID frame when switching from HQ core or AMR-WB IO mode into inactive frame in ACELP core when DTX is on */
            st->core_brate = SID_2k40;
        }

        if( st->cng_type == FD_CNG && st->total_brate <= ACELP_24k40 )     /* at highest bit-rates, use exclusively LP_CNG */
        {
            if ( st->total_brate == ACELP_9k60 || st->total_brate == ACELP_16k40 || st->total_brate == ACELP_24k40 )
            {
                st->codec_mode = MODE2;
            }
        }
        else
        {
            st->cng_type = LP_CNG;
            if ( st->codec_mode == MODE2 )
            {
                st->lp_cng_mode2 = 1;
            }
            st->codec_mode = MODE1;
        }
    }

    /*------------------------------------------------------------------------*
     * Reset counters when in active frame (neither SID nor FRAME_NO_DATA frame)
     *------------------------------------------------------------------------*/

    if( st->core_brate != SID_2k40 && st->core_brate != SID_1k75 && st->core_brate != 0 )
    {
        st->cnt_SID = 0;

        /* change SID update rate */
        /* first SID update is only 8 (3 in AMR-WB IO mode) frames after the active speech end */
        if( !st->Opt_AMR_WB )
        {
            st->max_SID = FIXED_SID_RATE;
        }
        else
        {
            st->max_SID = 3;
        }

        if ( st->max_SID > st->interval_SID )
        {
            st->max_SID = st->interval_SID;
        }

        /* reset the counter of CNG frames for averaging */
        st->cng_cnt = 0;

        if( st->active_cnt >= CNG_TYPE_HO && !st->Opt_AMR_WB )
        {
            if ( st->cng_type == LP_CNG && ( (st->input_bwidth == NB && st->bckr_tilt_lt > 9.f) || (st->input_bwidth > NB && st->bckr_tilt_lt > 45.f) ) )
            {
                st->cng_type = FD_CNG;
            }
            else if ( st->cng_type == FD_CNG && ( (st->input_bwidth == NB && st->bckr_tilt_lt < 2.f) || (st->input_bwidth > NB && st->bckr_tilt_lt < 10.f) ) )
            {
                st->cng_type = LP_CNG;
            }
            st->last_total_brate_cng = -1;
        }
        else if( st->Opt_AMR_WB )
        {
            st->cng_type = LP_CNG;
        }

        st->active_cnt++;
        st->active_cnt = min(st->active_cnt,200);
    }

    /*------------------------------------------------------------------------*
     * Update speech and background noise long-term energy
     *------------------------------------------------------------------------*/

    st->frame_ener = 0.0f;

    if ( st->Opt_DTX_ON )
    {
        st->frame_ener = sum2_f( speech, L_FRAME );

        /* Active speech (voiced) */
        if ( st->clas == VOICED_CLAS )
        {
            alpha = ALPHA_ENER_SLOW;
            if ( st->frame_ener > st->lt_ener_voiced )
            {
                alpha = ALPHA_ENER_FAST;
            }

            st->lt_ener_voiced = alpha * st->lt_ener_voiced + (1.0f-alpha) * st->frame_ener;

            (st->VarDTX_cnt_voiced)++;

            if (st->VarDTX_cnt_voiced > MIN_CNT)
            {
                st->VarDTX_cnt_voiced = MIN_CNT;
            }
        }

        /* Background noise */
        else if( !st->Opt_AMR_WB )
        {
            alpha = ALPHA_ENER_SLOW;
            if (st->frame_ener < st->lt_ener_noise)
            {
                alpha = ALPHA_ENER_FAST;
            }

            st->lt_ener_noise = alpha * st->lt_ener_noise + (1.0f-alpha) * st->frame_ener;

            (st->VarDTX_cnt_noise)++;

            if ( st->VarDTX_cnt_noise > MIN_CNT)
            {
                st->VarDTX_cnt_noise = MIN_CNT;
            }
        }
    }

    /* Update of the SID counter */
    update_SID_cnt( st );

    /* Update encoded bandwidth */
    if(  st->Opt_DTX_ON && (st->core_brate == SID_2k40 || st->core_brate == FRAME_NO_DATA ) )
    {
        st->bwidth = st->last_bwidth;
        if( st->last_core_brate > SID_2k40 && st->last_total_brate_cng != -1 )
        {
            st->bwidth = st->last_bwidth_cng;
        }

        if( st->Opt_RF_ON && st->total_brate == ACELP_13k20 && st->bwidth == NB )
        {
            st->codec_mode = MODE1;
            reset_rf_indices(st);
            st->Opt_RF_ON = 0;
            st->rf_mode = 0;
        }

        if ( st->codec_mode == MODE2 )
        {
            short n, bits_frame_nominal, tmpBandwidthMin;

            bits_frame_nominal = st->total_brate / 50;
            for( n=0; n<FRAME_SIZE_NB; n++ )
            {
                if( FrameSizeConfig[n].frame_bits == bits_frame_nominal )
                {
                    break;
                }
            }
            if( n == FRAME_SIZE_NB )
            {
                assert(!"Bitrate not supported: not part of EVS");
            }

            tmpBandwidthMin = FrameSizeConfig[n].bandwidth_min;

            if( st->rf_mode )
            {
                tmpBandwidthMin = WB;
            }

            st->bwidth = max(min(st->bwidth, FrameSizeConfig[n].bandwidth_max), tmpBandwidthMin);
        }

    }

    return;
}


/*---------------------------------------------------------------------*
 * update_SID_cnt()
 *
 * Update of the SID counter
 *---------------------------------------------------------------------*/

static void update_SID_cnt(
    Encoder_State *st       /* i/o: State structure                     */
)
{
    float EstimatedSNR, delta;

    if( st->core_brate == SID_2k40 || st->core_brate == SID_1k75 || st->core_brate == FRAME_NO_DATA )
    {
        /* Adapt the SID interval */
        if ( st->var_SID_rate_flag && st->VarDTX_cnt_voiced == MIN_CNT && st->VarDTX_cnt_noise == MIN_CNT )
        {
            EstimatedSNR = 10.0f * (float)log10( (0.01f + st->lt_ener_voiced) / (0.01f + st->lt_ener_noise) );
            if ( EstimatedSNR > SNR_H )
            {
                st->interval_SID = INT_H;
            }
            else if ( EstimatedSNR < SNR_L )
            {
                st->interval_SID = INT_L;
            }
            else
            {
                st->interval_SID = INT_L + (short)((INT_H - INT_L) * (EstimatedSNR - SNR_L)/(SNR_H - SNR_L));
            }

            if( !st->Opt_AMR_WB || st->max_SID != 3 )
            {
                st->max_SID = st->interval_SID;      /* change SID update rate */
            }
        }

        if( st->Opt_DTX_ON == 1 && st->cnt_SID != 0 )
        {
            /* Send SID frame only if long-term energy variation is above threshold */
            delta = 10.0f * (float)log10((0.01f + st->lt_ener_noise)/(0.01f + st->lt_ener_last_SID));
            if ( delta < LTE_VAR && st->VarDTX_cnt_voiced == MIN_CNT && st->VarDTX_cnt_noise == MIN_CNT )
            {
                /* Send SID frame, and reset lt_ener_noise */
                st->lt_ener_noise = st->frame_ener;
            }
        }
        else
        {
            /* If SID frame was sent, update long-term energy */
            st->lt_ener_last_SID = st->lt_ener_noise;
        }

        (st->cnt_SID)++;

        if( st->var_SID_rate_flag )
        {
            if( st->Opt_AMR_WB && st->max_SID == 3 && st->cnt_SID == 3 )
            {
                /* set the size of CNG history buffer for averaging to 3 frames */
                st->cng_hist_size = 3;
            }
            else if ( st->max_SID != 3 && st->cnt_SID == DTX_HIST_SIZE )
            {
                /* set the size of CNG history buffer for averaging to DTX_HIST_SIZE frames */
                /* be sure that DTX_HIST_SIZE >= INT_L */
                st->cng_hist_size = DTX_HIST_SIZE;
            }
        }

        if( !st->var_SID_rate_flag && st->interval_SID > 1 )
        {
            /* set the size of CNG history buffer for averaging to interval_SID frames */
            st->cng_hist_size = st->interval_SID;
            if (st->cng_hist_size > DTX_HIST_SIZE)
            {
                st->cng_hist_size = DTX_HIST_SIZE;
            }
        }

        if( st->cnt_SID >= st->max_SID )
        {
            /* adaptive SID update interval */
            st->max_SID = st->interval_SID;
            st->cnt_SID = 0;
        }
    }

    return;
}


/*-------------------------------------------------------------------*
 * dtx_hangover_control()
 *
 *
 *-------------------------------------------------------------------*/

void dtx_hangover_control(
    Encoder_State *st,            /* i/o: encoder state structure  */
    const float lsp_new[M]      /* i  : current frame LSPs       */
)
{
    short ptr;
    short i, j, m;
    float tmp_lsp[max(DTX_HIST_SIZE,HO_HIST_SIZE)*M];
    float tmp_enr[max(DTX_HIST_SIZE,HO_HIST_SIZE)];
    float tmp[max(DTX_HIST_SIZE,HO_HIST_SIZE)*M];
    float enr_new;
    float weights;
    float enr_est, lsp_est[M];
    float Dlsp, Denr;
    float lsf_tmp[M];
    float C[M];
    float max[2];
    short max_idx[2];
    float ftmp;
    float Dlsp_n2e, Denr_n2e;

    /* get current frame exc energy in log2 */
    enr_new = (float)(log10(st->ho_ener_circ[st->ho_circ_ptr]) / log10(2.0f));

    if( enr_new < 0.0f )
    {
        enr_new = 0.0f;
    }

    /* get energies and lsps of hangover frames  */
    ptr = st->ho_circ_ptr - (st->burst_ho_cnt-1);
    if( ptr < 0 )
    {
        ptr += st->ho_circ_size;
    }

    for( i=0; i<st->burst_ho_cnt-1; i++ )
    {
        mvr2r( &(st->ho_lsp_circ[ptr*M]), &(tmp_lsp[i*M]), M );
        tmp_enr[i] = st->ho_ener_circ[ptr];

        ptr++;
        if( ptr == st->ho_circ_size )
        {
            ptr = 0;
        }
    }

    /* get estimated CNG energy and lsps assuming terminate hangover at current frame */
    ptr = st->burst_ho_cnt - 2;
    enr_est = W_DTX_HO[0] * tmp_enr[ptr];
    weights = W_DTX_HO[0];
    mvr2r( &(tmp_lsp[ptr*M]), tmp, M );
    m = 1;

    for( i = 1; i < st->burst_ho_cnt-2; i++ )
    {
        if ( tmp_enr[ptr-i] < tmp_enr[ptr] * BUF_H_NRG && tmp_enr[ptr-i] > tmp_enr[ptr] * BUF_L_NRG )
        {
            enr_est += W_DTX_HO[i] * tmp_enr[ptr-i];
            weights += W_DTX_HO[i];
            mvr2r( &tmp_lsp[(ptr-i)*M], &tmp[m*M], M );
            m++;
        }
    }

    enr_est /= weights;

    if( enr_est < 1.0f )
    {
        enr_est = 1.0f;
    }

    Denr_n2e = (float)fabs(enr_new - log10(enr_est) / log10(2.0f));

    if( m < 3 )
    {
        enr_est = 0.8f*enr_est + (1-0.8f)*st->ho_ener_circ[st->ho_circ_ptr];
    }
    else
    {
        enr_est = 0.95f*enr_est + (1-0.95f)*st->ho_ener_circ[st->ho_circ_ptr];
    }

    enr_est = (float)(log10(enr_est) / log10(2.0f));

    if( enr_est < 0.0f )
    {
        enr_est = 0.0f;
    }

    set_f( max, 0.0f, 2 );
    set_s( max_idx, 0, 2 );

    for( i=0; i<m; i++ )
    {
        if( st->L_frame == L_FRAME )
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

        C[i] *= 0.0588235f; /* 0.0588235f = 1/(M+1) */

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

    if( m == 1 )
    {
        mvr2r( tmp, lsp_est, M );
    }
    else if( m < 4 )
    {
        for( i=0; i<M; i++ )
        {
            lsp_est[i] = 0.0f;
            for( j=0; j<m; j++ )
            {
                lsp_est[i] += tmp[j*M+i];
            }

            lsp_est[i] -= tmp[max_idx[0]*M+i];
            lsp_est[i] /= (float)(m - 1);
        }
    }
    else
    {
        for( i=0; i<M; i++ )
        {
            lsp_est[i] = 0.0f;
            for( j=0; j<m; j++ )
            {
                lsp_est[i] += tmp[j*M+i];
            }

            lsp_est[i] -= (tmp[max_idx[0]*M+i] + tmp[max_idx[1]*M+i]);
            lsp_est[i] /= (float)(m - 2);
        }
    }

    Dlsp_n2e = 0.0f;
    for( i=0; i<M; i++ )
    {
        Dlsp_n2e += (float)fabs(lsp_new[i] - lsp_est[i]);
        lsp_est[i] = 0.8f * lsp_est[i] + (1-0.8f) * lsp_new[i];
    }

    /* get deviation of CNG parameters between newly estimated and current state memory */
    Dlsp = 0.0f;
    max[0] = 0.0f;

    for( i=0; i<M; i++ )
    {
        Dlsp += (float)fabs(st->lspCNG[i] - lsp_est[i]);
        if ( fabs(st->lspCNG[i] - lsp_est[i]) > max[0] )
        {
            max[0] = (float)fabs(st->lspCNG[i] - lsp_est[i]);
        }
    }
    Denr = (float)fabs((log10(st->lp_ener+0.1f)/log10(2.0f)) - enr_est);

    /* make decision if DTX hangover can be terminated */
    st->hangover_terminate_flag = 0;

    if (( Dlsp < 0.4f && Denr < 1.4f && max[0] < 0.1f && Dlsp_n2e < 0.4f && Denr_n2e < 1.2f && st->Opt_SC_VBR)||
            ( Dlsp < 0.4f && Denr < 0.8f && max[0] < 0.1f && Dlsp_n2e < 0.4f && Denr_n2e < 0.8f && !st->Opt_SC_VBR))
    {
        st->hangover_terminate_flag = 1;
    }

    return;
}
