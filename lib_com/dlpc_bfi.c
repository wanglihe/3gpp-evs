/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "options.h"
#include "prot.h"
#include "rom_com.h"


/*---------------------------------------------------------------------*
* routine:   dlpc_bfi()
*
*
*---------------------------------------------------------------------*/

void dlpc_bfi(
    int L_frame,
    float *lsf_q,            /* o  : quantized lsfs                         */
    const float *lsfold,           /* i  : past quantized lsf                     */
    const short last_good,         /* i  : last good received frame               */
    const short nbLostCmpt,        /* i  : counter of consecutive bad frames      */
    float mem_MA[],          /* i/o: quantizer memory for MA model          */
    float *mem_AR,           /* i/o: quantizer memory for AR model          */
    float *stab_fac,         /* i  : lsf stability factor                   */
    float *lsf_adaptive_mean,/* i  : lsf adaptive mean, updated when BFI==0 */
    int   numlpc,            /* i  : Number of division per superframe      */
    float lsf_cng[],
    int   plcBackgroundNoiseUpdated,
    float *lsf_q_cng,        /* o  : quantized lsfs of background noise      */
    float *old_lsf_q_cng,    /* o  : old quantized lsfs for background noise */
    const float lsfBase[]          /* i  : base for differential lsf coding        */
)
{
    lsf_dec_bfi( MODE2, &lsf_q[0], lsfold, lsf_adaptive_mean, lsfBase, mem_MA, mem_AR, *stab_fac, 0, L_frame, last_good,
                 nbLostCmpt, plcBackgroundNoiseUpdated, lsf_q_cng, lsf_cng, old_lsf_q_cng, 0, 0, 0 );

    if( numlpc == 2 )
    {
        /* Decode the second LPC */
        lsf_dec_bfi( MODE2, &lsf_q[M], &lsf_q[0], lsf_adaptive_mean, lsfBase, mem_MA, mem_AR, *stab_fac, 0, L_frame, last_good,
                     nbLostCmpt+1, plcBackgroundNoiseUpdated, lsf_q_cng, lsf_cng, old_lsf_q_cng, 0, 0, 0 );
    }

    return;
}
