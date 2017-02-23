/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"

/*-------------------------------------------------------------------*
 * FEC_lsf_estim()
 *
 * - LSP calculation
 * - A(z) calculation
 *-------------------------------------------------------------------*/

void FEC_lsf2lsp_interp(
    Decoder_State *st,              /* i/o: Decoder static memory                        */
    const short L_frame,          /* i  : length of the frame                          */
    float *Aq,              /* o  : calculated A(z) for 4 subframes              */
    float *lsf,             /* o  : estimated LSF vector                         */
    float *lsp              /* o  : estimated LSP vector                         */
)
{
    /* convert LSFs to LSPs */
    if ( st->Opt_AMR_WB )
    {
        isf2isp( lsf, lsp, M, INT_FS_12k8 );
    }
    else
    {
        if( L_frame == L_FRAME )
        {
            lsf2lsp( lsf, lsp, M, INT_FS_12k8 );
        }
        else /* L_frame == L_FRAME16k */
        {
            lsf2lsp( lsf, lsp, M, INT_FS_16k );
        }
    }

    /*----------------------------------------------------------------------*
     * Interpolate LSP vector and find A(z)
     *----------------------------------------------------------------------*/

    if ( st->Opt_AMR_WB )
    {
        int_lsp( L_frame, st->lsp_old, lsp, Aq, M, interpol_isp_amr_wb, 1 );
    }
    else
    {
        int_lsp( L_frame, st->lsp_old, lsp, Aq, M, interpol_frac_12k8, 0 );
    }

    return;
}
