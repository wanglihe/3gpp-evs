/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "cnst.h"
#include "options.h"
#include "prot.h"
#include "rom_com.h"


/*-------------------------------------------------------------------*
 * Interpol_delay()
 *
 * Interpolation of pitch lag
 *--------------------------------------------------------------------*/

void Interpol_delay(
    float *out,            /* o : pitch interpolation output   */
    float *last,           /* i : last frame pitch lag         */
    float *current,        /* i : current frame pitch lag      */
    short SubNum,          /* i : subframe number              */
    const float *frac            /* i : interpolation constant       */
)
{
    out[0] = (1.0f - frac[SubNum]) **last + frac[SubNum] **current;
    out[1] = (1.0f - frac[SubNum + 1]) **last + frac[SubNum + 1] **current;
    out[2] = (1.0f - frac[SubNum + 2]) **last + frac[SubNum + 2] **current;

    return;
}


/*-------------------------------------------------------------------*
 * deemph_lpc()
 *
 * De-emphasis of LP coefficients
 * convolve LPC with [1 -PREEMPH_FAC] to de-emphasise LPC
 *--------------------------------------------------------------------*/

void deemph_lpc(
    float *p_Aq_curr,           /* i : LP coefficients current frame                 */
    float *p_Aq_old,            /* i : LP coefficients previous frame                */
    float *LPC_de_curr,         /* o : De-emphasized LP coefficients current frame   */
    float *LPC_de_old           /* o : De-emphasized LP coefficients previous frame  */
    ,short deemph_old

)
{
    short k;
    float b[M+2];
    float a[2] = {-PREEMPH_FAC,1.0};

    b[0] = 1.0;
    for(k = 0; k < M; k++)
    {
        b[k+1] = p_Aq_curr[k];
    }
    b[M+1] = 0.0;

    for(k = 0; k <= M; k++)
    {
        LPC_de_curr[k] = a[0]*b[k] + a[1]*b[k+1];
    }

    if ( deemph_old == 1)
    {
        /* ignoring the 1st value which is 1.0 in this case */
        b[0] = 1.0;
        for(k = 0; k < M; k++)
        {
            b[k+1] = p_Aq_old[k+1];
        }
        b[M+1] = 0.0;

        for(k = 0; k <= M; k++)
        {
            LPC_de_old[k] = a[0]*b[k] + a[1]*b[k+1];
        }
    }

    return;
}
