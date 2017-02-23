/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "prot.h"


/*---------------------------------------------------------*
 * stat_noise_uv_dec()
 *
 * Modifies excitation signal in UC mode when the noise is stationary
 *---------------------------------------------------------*/

void stat_noise_uv_dec(
    Decoder_State *st,          /* i/o: decoder static memory                */
    const short coder_type,   /* i  : coding type                          */
    float *lsp_new,     /* i  : end-frame LSP vector                 */
    float *lsp_mid,     /* i  : mid-frame LSP vector                 */
    float *Aq,          /* o  : A(z) quantized for the 4 subframes   */
    float *exc2         /* i/o: excitation buffer                    */
)
{
    short i;
    float ftmp, noisiness = 0;

    /*-----------------------------------------------------------------*
     * Decode the VAD flag
     *-----------------------------------------------------------------*/

    if( coder_type == UNVOICED || ( coder_type == INACTIVE && st->core_brate <= ACELP_9k60 ) )
    {
        /* read the noisiness parameter */
        noisiness = (float)get_next_indice( st, 5 );
    }

    /*-----------------------------------------------------------------*
     * Update long-term energies for FEC
     * Update LSP vector for CNG
     *-----------------------------------------------------------------*/

    if( coder_type == INACTIVE )
    {
        if( st->unv_cnt > 20 )
        {
            ftmp = st->lp_gainc * st->lp_gainc;
            st->lp_ener = 0.7f * st->lp_ener + 0.3f * ftmp;
            for( i=0; i<M; i++ )
            {
                st->lspCNG[i] = (float)(0.9f * st->lspCNG[i] + 0.1f * lsp_new[i]);
            }
        }
        else
        {
            st->unv_cnt++;
        }
    }
    else
    {
        st->unv_cnt = 0;
    }

    /*-----------------------------------------------------------------*
     * Modify the excitation signal
     *-----------------------------------------------------------------*/

    if ( !st->Opt_AMR_WB )
    {
        stat_noise_uv_mod( coder_type, noisiness, st->lsp_old, lsp_new, lsp_mid, Aq, exc2, 0, &st->ge_sm,
                           &st->uv_count, &st->act_count, st->lspold_s, &st->noimix_seed, &st->min_alpha,
                           &st->exc_pe, st->core_brate, st->bwidth );
    }

    return;
}
