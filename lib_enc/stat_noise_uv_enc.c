/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "prot.h"
#include "rom_com.h"

/*-----------------------------------------------------------------*
 * stat_noise_uv_enc()
 *
 * Modifies excitation signal when the noise is stationary
 *-----------------------------------------------------------------*/

void stat_noise_uv_enc(
    Encoder_State *st,           /* i/o: state structure                      */
    const short coder_type,    /* i  : coding type                          */
    const float *epsP,         /* i  : LP prediction errors                 */
    float *isp_new,      /* i  : immittance spectral pairs at 4th sfr */
    float *isp_mid,      /* i  : immittance spectral pairs at 2nd sfr */
    float *Aq,           /* i  : A(z) quantized for the 4 subframes   */
    float *exc2          /* i/o: excitation buffer                    */
)
{
    short noisiness = 0;

    /*-----------------------------------------------------------------*
     * Calculate and write the noisiness parameter
     *-----------------------------------------------------------------*/

    if ( coder_type == UNVOICED || ( coder_type == INACTIVE && st->core_brate <= ACELP_9k60 ) )
    {
        if ( st->bwidth != NB )
        {
            /* WB case */
            noisiness = (short)(((epsP[2] / epsP[16]) - 1) * 2.0f * 32);
        }
        else if ( coder_type == INACTIVE && st->bwidth == NB )
        {
            /* NB GSC case */
            noisiness = (short)(((epsP[2] / epsP[16]) - 1) * 0.25f * 32);
        }
        else
        {
            /* NB case */
            noisiness = (short)(((epsP[2] / epsP[16]) - 1) * 0.5f * 32);
        }

        if( noisiness < 0 )
        {
            noisiness = 0;
        }

        if( noisiness > 31 )
        {
            noisiness = 31;
        }

        push_indice( st, IND_NOISINESS, noisiness, 5 );
    }

    /*-----------------------------------------------------------------*
     * Modify the stationary noise excitation signal
     *-----------------------------------------------------------------*/

    stat_noise_uv_mod( coder_type, noisiness, st->lsp_old, isp_new, isp_mid, Aq, exc2, 0, &st->ge_sm,
                       &st->uv_count, &st->act_count, st->lspold_s, &st->noimix_seed, &st->min_alpha, &st->exc_pe,
                       st->core_brate, st->bwidth );

    return;
}
