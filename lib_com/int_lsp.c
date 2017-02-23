/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"

/*---------------------------------------------------------------------*
 * int_lsp()
 *
 * Find the interpolated LSP parameters for all subframes
 *---------------------------------------------------------------------*/

void int_lsp(
    const short L_frame,     /* i  : length of the frame               */
    const float lsp_old[],   /* i  : LSPs from past frame              */
    const float lsp_new[],   /* i  : LSPs from present frame           */
    float *Aq,         /* o  : LP coefficients in both subframes */
    const short m,           /* i  : order of LP filter                */
    const float *int_coeffs, /* i  : interpolation coefficients        */
    const short Opt_AMR_WB   /* i  : flag indicating AMR-WB IO mode    */
)
{
    float lsp[M], fnew, fold;
    short i, k;
    const float *pt_int_coeffs;

    if( L_frame == L_FRAME )
    {
        pt_int_coeffs = int_coeffs;
    }
    else /* L_frame == L_FRAME16k */
    {
        pt_int_coeffs = interpol_frac_16k;
    }

    for( k=0; k<L_frame/L_SUBFR; k++ )
    {
        fnew = pt_int_coeffs[k];
        fold = (float)(1.0f - fnew);

        for (i = 0; i < m; i++)
        {
            lsp[i] = lsp_old[i] * fold + lsp_new[i] * fnew;
        }

        if ( Opt_AMR_WB )
        {
            isp2a( lsp, Aq, m );
        }
        else
        {
            lsp2a_stab(lsp, Aq, m);
        }

        Aq += (m+1);
    }

    return;
}

/*---------------------------------------------------------------------*
 * int_lsp4()
 *
 * Interpolate LSPs find the A[z] parameters for all subframes by interpolating between
 * old end-frame LSPs, current mid-frame LSPs and current end-frame LSPs
 *---------------------------------------------------------------------*/

void int_lsp4(
    const short L_frame,    /* i  : length of the frame               */
    const float lsp_old[],  /* i  : LSPs from past frame              */
    const float lsp_mid[],  /* i  : LSPs from mid-frame               */
    const float lsp_new[],  /* i  : LSPs from present frame           */
    float *Aq,        /* o  : LP coefficients in both subframes */
    const short m,          /* i  : order of LP filter                */
    short relax_prev_lsf_interp /* i  : relax prev frame lsf interp after erasure */
)
{
    float lsp[M];
    short i;
    short j;
    const float *pt_int_coeffs;

    if( L_frame == L_FRAME )
    {
        if ( relax_prev_lsf_interp == 1 )
        {
            pt_int_coeffs = interpol_frac_mid_relaxprev_12k8;
        }
        else if ( relax_prev_lsf_interp == 2 )
        {
            pt_int_coeffs = interpol_frac_mid_FEC;
        }
        else if ( relax_prev_lsf_interp == -1 )
        {
            pt_int_coeffs = interpol_frac_mid_relaxprev_pred_12k8;
        }
        else
        {
            pt_int_coeffs = interpol_frac_mid;
        }
    }
    else /* L_frame == L_FRAME16k */
    {
        if ( relax_prev_lsf_interp == 1 )
        {
            pt_int_coeffs = interpol_frac_mid_relaxprev_16k;
        }
        else if ( relax_prev_lsf_interp == 2 )
        {
            pt_int_coeffs = interpol_frac_mid_16k_FEC;
        }
        else if ( relax_prev_lsf_interp == -1 )
        {
            pt_int_coeffs = interpol_frac_mid_relaxprev_pred_16k;
        }
        else
        {
            pt_int_coeffs = interpol_frac_mid_16k;
        }
    }

    for( j=0; j<L_frame/L_SUBFR; j++ )
    {
        for( i=0; i<m; i++ )
        {
            lsp[i] = lsp_old[i]*(*pt_int_coeffs) + lsp_mid[i]*(*(pt_int_coeffs+1)) + lsp_new[i]*(*(pt_int_coeffs+2));
        }
        pt_int_coeffs += 3;

        lsp2a_stab( lsp, Aq, m );

        Aq += (m+1);
    }

    return;
}
