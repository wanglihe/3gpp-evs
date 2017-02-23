/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"

/*-------------------------------------------------------------------*
 * Local constants
 *-------------------------------------------------------------------*/

#define ALPHA0  0.5f
#define BETA0  (1.0f-ALPHA0)

/*-------------------------------------------------------------------*
 * inact_switch_ematch()
 *
 * Apply energy matching when swithcing to INACTIVE frame coded by the GSC technology
 *-------------------------------------------------------------------*/

void inact_switch_ematch(
    float exc2[],               /* i/o: CELP/GSC excitation buffer   */
    float dct_exc_tmp[],        /* i  : GSC excitation in DCT domain */
    float lt_ener_per_band[],   /* i/o: Long term energy per band    */
    const short coder_type,           /* i  : Coding mode                  */
    const short L_frame,              /* i  : Frame lenght                 */
    const long  core_brate,           /* i  : Core bit rate                */
    const short bfi                   /* i  : frame lost indicator         */
    ,const short last_core,            /* i  : Last core used               */
    const short last_codec_mode       /* i  : Last codec mode              */
)
{
    float Ener_per_bd[MBANDS_GN];
    float ftmp;
    float *pt_exc;
    short j, i;


    /*--------------------------------------------------------------------------
     * average energy per band
     *--------------------------------------------------------------------------*/

    if( coder_type == AUDIO && bfi == 0)
    {
        Ener_per_band_comp( dct_exc_tmp, Ener_per_bd, MBANDS_GN, 1 );

        /* reset long-term energy per band */
        for( i = 0; i < MBANDS_GN; i++ )
        {
            lt_ener_per_band[i] = Ener_per_bd[i];
        }
    }
    else if( coder_type == VOICED || coder_type == GENERIC || coder_type == TRANSITION || last_core != ACELP_CORE || last_codec_mode != MODE1 )
    {
        /* Find spectrum and energy per band for GC and VC frames */
        edct( exc2, dct_exc_tmp, L_frame );
        Ener_per_band_comp( dct_exc_tmp, Ener_per_bd, MBANDS_GN, 1 );

        /* reset long-term energy per band */
        for(i = 0; i < MBANDS_GN; i++)
        {
            lt_ener_per_band[i] = Ener_per_bd[i];
        }
    }
    else if( coder_type == INACTIVE && core_brate <= ACELP_24k40 )
    {
        /* Find spectrum and energy per band for inactive frames */
        edct( exc2, dct_exc_tmp, L_frame );
        Ener_per_band_comp( dct_exc_tmp, Ener_per_bd, MBANDS_GN, 1 );

        /* More agressive smoothing in the first 50 frames */
        pt_exc = dct_exc_tmp;
        for( i = 0; i < MBANDS_GN; i++ )
        {
            /* Compute smoothing gain to apply with gain limitation */
            lt_ener_per_band[i] = ALPHA0*lt_ener_per_band[i] + BETA0*Ener_per_bd[i];

            ftmp = lt_ener_per_band[i] - Ener_per_bd[i];
            ftmp = (float)pow(10, ftmp);

            if( i < 2 )
            {
                for (j = 0; j < 8; j ++)
                {
                    *pt_exc  *= ftmp;
                    pt_exc++;
                }
            }
            else
            {
                for (j = 0; j < 16; j ++)
                {
                    *pt_exc  *= ftmp;
                    pt_exc++;
                }
            }
        }

        /* Going back to time */
        edct( dct_exc_tmp, exc2, L_frame );
    }

    return;
}
