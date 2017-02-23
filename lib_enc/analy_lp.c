/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"

/*-------------------------------------------------------------------*
 * analy_lp()
 *
 * Perform LP analysis
 *
 * - autocorrelations + lag windowing
 * - Levinson-Durbin algorithm to find A(z)
 * - convert A(z) to LSPs
 * - find interpolated LSPs and convert back to A(z) for all subframes
 * - update LSPs for the next frame
 *-------------------------------------------------------------------*/

void analy_lp(
    const float speech[],     /* i  : pointer to the speech frame                         */
    const short L_frame,      /* i  : length of the frame                                 */
    const short L_look,       /* i  : look-ahead                                          */
    float *ener,        /* o  : residual energy from Levinson-Durbin                */
    float A[],          /* o  : A(z) filter coefficients                            */
    float epsP[],       /* o  : LP analysis residual energies for each iteration    */
    float lsp_new[],    /* o  : current frame LSPs                                  */
    float lsp_mid[],    /* o  : current mid-frame LSPs                              */
    float lsp_old[],    /* i/o: previous frame unquantized LSPs                     */
    const short Top[2],       /* i  : open loop pitch lag                                 */
    const float Tnc[2],       /* i  : open loop pitch gain                                */
    const float sr_core       /* i  : internal sampling rate                              */
)
{
    short i_subfr, wind_length, half_frame;
    float r[M+1], *lsp;
    const float *wind, *pt;
    short half_frame_idx;

    if( L_frame == L_FRAME )
    {
        wind_length = L_LP;
        wind = LP_assym_window;
    }
    else  /* L_frame == L_FRAME16k */
    {
        wind_length = L_LP_16k;
        wind = LP_assym_window_16k;
    }
    lsp = lsp_mid;
    half_frame = L_frame>>1;

    half_frame_idx = 0;

    for( i_subfr = half_frame; i_subfr <= L_frame; i_subfr = i_subfr + half_frame )
    {
        pt = speech + i_subfr + L_look - wind_length;

        /* Autocorrelations */
        autocorr( pt, r, M, wind_length, wind, 0, 0, 0 );

        /* Lag windowing */
        adapt_lag_wind( r, M, Top[half_frame_idx], Tnc[half_frame_idx], sr_core );
        ++half_frame_idx;

        /* Levinson-Durbin */
        lev_dur( A, r, M, epsP );

        /* Conversion of A(z) to LSPs */
        a2lsp_stab( A, lsp, lsp_old );

        lsp = lsp_new;
    }

    /* LSP interpolation */
    int_lsp4( L_frame, lsp_old, lsp_mid, lsp_new, A, M, 0 );

    /* updates */
    mvr2r( lsp_new, lsp_old, M );

    *ener = epsP[M];

    return;
}


/*-------------------------------------------------------------------*
 * analy_lp_AMR_WB()
 *
 * Perform LP analysis for AMR-WB IO mode
 *
 * - autocorrelations + lag windowing
 * - Levinson-Durbin algorithm to find A(z)
 * - convert A(z) to ISPs
 * - find interpolated ISPs and convert back to A(z) for all subframes
 * - update ISPs for the next frame
 *-------------------------------------------------------------------*/

void analy_lp_AMR_WB(
    const float speech[],    /* i  : pointer to the speech frame                      */
    float *ener,       /* o  : residual energy from Levinson-Durbin             */
    float A[],         /* o  : A(z) filter coefficients                         */
    float epsP[],      /* o  : LP analysis residual energies for each iteration */
    float isp_new[],   /* o  : current frame ISPs                               */
    float isp_old[],   /* i/o: previous frame unquantized ISPs                  */
    float isf_new[],   /* o  : current frame ISFs                               */
    const int   Top,         /* i  : open loop pitch lag                              */
    const float Tnc          /* i  : open loop pitch gain                             */
)
{
    short wind_length;
    float r[M+1];
    const float *wind;

    /* Initialization */
    wind_length = L_LP_AMR_WB;
    wind = hamcos_window;

    /* Autocorrelations */
    autocorr( speech - L_SUBFR, r, M, wind_length, wind, 0, 0, 0 );

    /* Lag windowing */
    adapt_lag_wind( r, M, Top, Tnc, 12800 );

    /* Levinson-Durbin  */
    lev_dur( A, r, M, epsP );

    a2isf( A, isf_new, stable_ISF, M);
    isf2isp( isf_new, isp_new, M, INT_FS_12k8 );

    /* ISP interpolation */
    int_lsp( L_FRAME, isp_old, isp_new, A, M, interpol_isp_amr_wb, 1 );

    *ener = epsP[M];

    /* updates */
    mvr2r( isp_new, isp_old, M );

    return;
}

