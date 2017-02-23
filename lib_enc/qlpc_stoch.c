/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"
#include "rom_enc.h"
#include "basop_proto_func.h"

/*-------------------------------------------------------------------*
* lpc_quantization()
*
*
*--------------------------------------------------------------------*/

void lpc_quantization(
    Encoder_State * st,
    int core,
    int lpcQuantization,
    float lsf_old[],
    const float lsp[],
    const float lspmid[],
    float lsp_q[],
    float lsf_q[],
    float lspmid_q[],
    float mem_MA[],
    float mem_AR[],
    int narrowBand,
    short coder_type,
    int acelp_midLpc,
    int param_lpc[],
    int nbits_lpc[],
    short *seed_acelp,
    int sr_core,
    float *Bin_Ener,
    float *Bin_Ener_old,
    short * bits_param_lpc,
    short * no_param_lpc
)
{
    int nb_indices=0;
    float lsfmid_q[M];
    short lsfmid_idx;
    int i;
    float lsf[M], lsfmid[M];
    short force_sf;

    float fec_lsf[M], stab;


    /****** High-rate LPC quantizer *******/
    if( lpcQuantization == 0 )
    {
        lsp2lsf( lsp, lsf, M, sr_core );

        if( core == TCX_10_CORE )
        {
            lsp2lsf( lspmid, lsfmid, M, sr_core );
        }

        qlpc_avq( lsf, lsfmid, lsf_q, lsfmid_q, param_lpc, &nb_indices, nbits_lpc, core, sr_core );

        lsf2lsp( lsf_q, lsp_q, M, sr_core );

        if( core == TCX_10_CORE )
        {
            lsf2lsp( lsfmid_q, lspmid_q, M, sr_core );
        }

        assert(nb_indices<=NPRM_LPC_NEW);
    }

    /****** Low-rate LPC quantizer *******/

    else if( lpcQuantization == 1 )
    {

        lsp2lsf( lsp, lsf, M, sr_core );

        force_sf = 0;
        /*Force safety net when possible in case of transitions*/
        if( st->tc_cnt >= 1 || st->last_core_brate <= SID_2k40 || st->next_force_safety_net )
        {
            force_sf = 1;
            st->next_force_safety_net = 0;
        }

        if( st->next_force_safety_net == 1 && st->Opt_RF_ON )
        {
            force_sf = 1;
            st->next_force_safety_net = 0;
        }

        if( sr_core == INT_FS_16k && coder_type == UNVOICED )
        {
            lsf_end_enc( st, lsf, lsf_q, mem_AR, mem_MA, ENDLSF_NBITS, GENERIC, st->bwidth, Bin_Ener, sr_core, st->core_brate,
                         &st->streaklimit, &st->pstreaklen, force_sf, 0, 1, param_lpc, no_param_lpc, bits_param_lpc, GENERIC );

            nb_indices = (int)(*no_param_lpc);
        }
        else
        {
            lsf_end_enc( st, lsf, lsf_q, mem_AR, mem_MA, ENDLSF_NBITS, coder_type, st->bwidth, Bin_Ener, sr_core, st->core_brate,
                         &st->streaklimit, &st->pstreaklen, force_sf, 0, 1, param_lpc, no_param_lpc, bits_param_lpc, coder_type );

            nb_indices = (int)(*no_param_lpc);
        }

        FEC_lsf_estim_enc( st, st->L_frame, fec_lsf );

        /* FEC - calculate LSF stability */
        stab = lsf_stab( lsf_q, fec_lsf, 0, st->L_frame);

        if (stab < ( STAB_FAC_LIMIT + 0.2 ) && ( coder_type == VOICED  || coder_type == GENERIC) && st->Opt_RF_ON )
        {
            st->next_force_safety_net = 1;
        }

        lsf2lsp( lsf_q, lsp_q, M, sr_core );

        *nbits_lpc = ENDLSF_NBITS;
    }
    else
    {
        assert(0);
    }

    *seed_acelp = 0;
    for( i=nb_indices-1; i>=0; i-- )
    {
        /* rightshift before *seed_acelp+param_lpc[i] to avoid overflows*/
        *seed_acelp=(short)((((*seed_acelp)>>1)+param_lpc[i]) * 31821L + 13849L);
    }

    /* Mid-frame LPC quantization */
    if( lpcQuantization && acelp_midLpc )
    {
        if( st->rate_switching_reset == 0 )
        {
            lsp2lsf( lspmid, lsfmid, M, sr_core );
            midlsf_enc( lsf_old, lsf_q, lsfmid, &lsfmid_idx, M, Bin_Ener_old, narrowBand, sr_core, coder_type );
            param_lpc[nb_indices++] = (int)lsfmid_idx;

            midlsf_dec( lsf_old, lsf_q, lsfmid_idx, lsfmid_q, M, coder_type, NULL, 0, 1 );

            reorder_lsf( lsfmid_q, LSF_GAP_MID, M, sr_core );
            lsf2lsp( lsfmid_q, lspmid_q, M, sr_core );
        }
        else
        {
            param_lpc[nb_indices++] = 0;
        }
    }

    return;
}



/*-------------------------------------------------------------------*
 * Unified_weighting()
 *
 * LSF weighting
 *-------------------------------------------------------------------*/

void Unified_weighting(
    float Bin_Ener_128[], /* i  : FFT Bin energy 128 bins in two sets     */
    const float lsf[],          /* i  : LSF vector                              */
    float w[],            /* o  : LP weighting filter (numerator)         */
    const short narrowBand,     /* i  : flag for Narrowband                     */
    const short unvoiced,       /* i  : flag for Unvoiced frame                 */
    const short sr_core,        /* i  : sampling rate of core-coder             */
    const int   order           /* i  : LP order                                */
)
{
    short i;
    float nf;
    const Word16 (*ptr_lsf_fit_model)[M];
    int norm_lsf[M];
    float tmp, min;
    float w_fft[M];
    float Bin_Ener_160[160];
    float *Bin_Ener;
    const float *Freq_w_Table;
    short last_bin;
    /*float compen;*/


    /*Config. weighting*/
    if( narrowBand )
    {
        ptr_lsf_fit_model = lsf_unified_fit_model_nb;
        nf = 6400.f;

        last_bin = 127;
        Bin_Ener = Bin_Ener_128;
    }
    else if( sr_core == 12800 )
    {
        ptr_lsf_fit_model = lsf_unified_fit_model_wb;
        nf = 6400.f;

        last_bin = 127;
        Bin_Ener = Bin_Ener_128;
    }
    else
    {
        ptr_lsf_fit_model = lsf_unified_fit_model_wbhb;
        nf = 8000.f;

        /* Fill the missing part (128~159) of the bin energy */
        last_bin = 159;

        mvr2r( Bin_Ener_128, Bin_Ener_160, L_FFT/2 );

        /* Find average bin energy (32 Energy) */
        tmp = 0.f;
        for( i=95; i<127; i++ )
        {
            tmp += Bin_Ener_160[i];
        }

        tmp = tmp/32.f;
        for( i=127; i<160; i++ )
        {
            Bin_Ener_160[i] = tmp;
        }

        Bin_Ener = Bin_Ener_160;
    }

    /* 1) FFT weights*/
    if( unvoiced )
    {
        Freq_w_Table = Freq_Weight_UV;
    }
    else
    {
        Freq_w_Table = Freq_Weight_Com;
    }

    /* Use Envelope */
    min = 1000.0f;
    for( i=0; i<M; i++ )
    {
        norm_lsf[i] = (int)((lsf[i]/50.f) + 0.5f);

        if( norm_lsf[i] == 0 )
        {
            w_fft[i] = Bin_Ener[1];
        }
        else if( norm_lsf[i] == last_bin )
        {
            w_fft[i] = Bin_Ener[last_bin-1];
        }
        else
        {
            tmp = max(Bin_Ener[norm_lsf[i]], Bin_Ener[norm_lsf[i]-1]);
            w_fft[i] = max(Bin_Ener[norm_lsf[i]+1], tmp);
        }

        if( w_fft[i] < MIN_LOG_60dB )
        {
            w_fft[i] = (float)MIN_LOG_VAL_60dB;
        }
        else
        {
            w_fft[i] = (float)(10.0 * log10(w_fft[i]));
        }

        if ( w_fft[i] < min )
        {
            min = w_fft[i];
        }
    }

    for( i=0; i<M; i++ )
    {
        w_fft[i] = (float)(sqrt(w_fft[i]-min)+2.f);
        w_fft[i] *= Freq_w_Table[norm_lsf[i]];
    }

    /* 2) IHM weights*/
    w[0] = 1.f/( lsf[0] - 0 ) + 1.f/( lsf[1] - lsf[0] );
    for( i=1; i<order-1; i++ )
    {
        w[i] = 1.f/( lsf[i] - lsf[i-1] ) + 1.f/( lsf[i+1] - lsf[i] );
    }
    w[order-1] = 1.f/( lsf[order-1] - lsf[order-2] ) +
                 1.f/( nf - lsf[order-1] );

    /* 3) Fitting model combining the two weights*/
    for( i=0; i<order; i++ )
    {
        w[i] *= (nf/EVS_PI);
        w[i] = ((float)(ptr_lsf_fit_model[0][i])/(1<<10))
               + w[i]*((float)(ptr_lsf_fit_model[1][i])/(1<<15))
               + w[i]*w[i]*((float)(ptr_lsf_fit_model[2][i])/(1<<18))
               + w_fft[i]*((float)(ptr_lsf_fit_model[3][i])/(1<<12));

        if( w[i] < 1.f/(i+1) )
        {
            w[i] = 1.f/(i+1);
        }
    }

    return;
}

