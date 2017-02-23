/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdio.h>
#include <assert.h>
#include "prot.h"
#include "cnst.h"





/*----------------------------------------------------------------------------------*
* Damping_fact()
*
*  Estimate damping factor
*----------------------------------------------------------------------------------*/

float Damping_fact( const short coder_type, int nbLostCmpt, short last_good, float stab_fac, float *lp_gainp, int core)
{
    float alpha, gain;

    if(core == ACELP_CORE)
    {
        alpha = ALPHA_VT;/* rapid convergence to 0 */
        if( ( coder_type == UNVOICED) && (nbLostCmpt <= 3))       /* Clear unvoiced last good frame   */
        {
            alpha = ALPHA_UU;
        }
        else if( last_good == UNVOICED_CLAS )
        {
            if( nbLostCmpt == 1 )
            {
                /* If stable, do not decrease the energy, pitch gain = 0 */
                alpha = stab_fac * (1.0f - 2.0f*ALPHA_U) + 2.0f*ALPHA_U;  /* [0.8, 1.0] */
            }
            else if (nbLostCmpt == 2 )
            {
                alpha = ALPHA_U*1.5f;   /* 0.6 */
            }
            else
            {
                alpha = ALPHA_U;   /* 0.4 go rapidly to CNG gain, pitch gain = 0 */
            }
        }
        else if( last_good == UNVOICED_TRANSITION )
        {
            alpha = ALPHA_UT;
        }
        else if( (last_good == ONSET) && (nbLostCmpt <= 3) && (coder_type == GENERIC))
        {
            alpha = 0.8f;
        }
        else if( ( (last_good == VOICED_CLAS) || (last_good == ONSET) ) && (nbLostCmpt <= 3) )
        {
            alpha = ALPHA_V;      /* constant for the first 3 erased frames */
        }
        if(last_good >= VOICED_CLAS)
        {
            if( nbLostCmpt == 1 ) /* if first erased frame in a block, reset harmonic gain */
            {
                gain = (float)sqrt( *lp_gainp );  /* move pitch gain towards 1 for voiced to remove energy fluctuations */

                if( gain > 0.98f )
                {
                    gain = 0.98f;
                }
                else if( gain < 0.85f )
                {
                    gain = 0.85f;
                }
                alpha *= gain;
            }
            else if( nbLostCmpt == 2 )
            {
                alpha = (0.6f+0.35f*stab_fac)**lp_gainp;
            }
            else
            {
                *lp_gainp *= (0.7f+0.2f*stab_fac);
                alpha = *lp_gainp;
            }
        }
    }
    else
    {
        if (nbLostCmpt < 2)
        {
            alpha = (0.7f+0.3f*stab_fac);
        }
        else if (nbLostCmpt == 2)
        {
            alpha = (0.45f+0.4f*stab_fac);
        }
        else
        {
            alpha = 0.35f+0.4f*stab_fac;
        }
    }
    return alpha;
}

