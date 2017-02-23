/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"

/*---------------------------------------------------------------------*
* routine:   lsf_dec_bfi()
*
* Estimate the lsfs in case of FER
* Bad frame, all active speech coders
*---------------------------------------------------------------------*/

void lsf_dec_bfi(
    short codec_mode,                 /* i: : codec mode: MODE1 | MODE2              */
    float *lsf,                       /* o  : estimated LSF vector                   */
    const float *lsfold,              /* i  : past quantized lsf                     */
    float *lsf_adaptive_mean,         /* i  : lsf adaptive mean, updated when BFI==0 */
    const float lsfBase[],            /* i  : base for differential lsf coding       */
    float *mem_MA,                    /* i/o: quantizer memory for MA model          */
    float *mem_AR,                    /* o  : quantizer memory for AR model          */
    float stab_fac,                   /* i  : lsf stability factor                   */
    short last_coder_type,            /* i  : last coder type                        */
    short L_frame,                    /* i  : frame length                           */
    const short last_good,            /* i  : last good received frame               */
    const short nbLostCmpt,           /* i  : counter of consecutive bad frames      */
    int   plcBackgroundNoiseUpdated,  /* i  : background noise alreadyupdated?       */
    float *lsf_q_cng,                 /* o  : quantized lsfs of background noise     */
    float *lsf_cng,                   /* i  : long term target for fading to bg noise*/
    float *old_lsf_q_cng,             /* i  : old quantized lsfs for background noise*/
    short Last_GSC_pit_band_idx,      /* i  : AC mode (GSC) - Last pitch band index  */
    short Opt_AMR_WB,                 /* i  : IO flag                                */
    const short MODE1_bwidth          /* i  : coded bandwidth                        */
)
{
    short i;
    float alpha, tmp;
    float lsf_mean[M];
    const float* pt_meansForFading;
    const float* pt_meansForMemUpdate;
    float beta;

    /* init vectors */
    if (codec_mode == MODE1)
    {
        pt_meansForMemUpdate = &lsf_mean[0];

        if( Opt_AMR_WB )
        {
            pt_meansForFading = mean_isf_amr_wb;
        }
        else
        {
            if( L_frame == L_FRAME )
            {
                if( MODE1_bwidth == NB )
                {
                    pt_meansForFading = GENB_Ave;
                }
                else
                {
                    pt_meansForFading = GEWB_Ave;
                }
            }
            else
            {
                pt_meansForFading = GEWB2_Ave;
            }
        }
    }
    else
    {
        pt_meansForFading = pt_meansForMemUpdate = lsfBase;
        if (lsf_cng != NULL && plcBackgroundNoiseUpdated)
        {
            pt_meansForFading = lsf_cng;
        }
    }

    /*----------------------------------------------------------------------*
     * Initialize the alpha factor
     *----------------------------------------------------------------------*/

    if( nbLostCmpt <= 3 )
    {
        if(last_coder_type == UNVOICED)
        {
            /* clearly unvoiced */
            alpha = ALPHA_UU;
        }
        else if( last_coder_type == AUDIO || last_good == INACTIVE_CLAS )
        {
            if( Last_GSC_pit_band_idx > 0 && nbLostCmpt > 1 )
            {
                alpha = 0.8f;
            }
            else
            {
                alpha = 0.995f;
            }
        }
        else if(last_good == UNVOICED_CLAS )
        {
            if( nbLostCmpt <= 1 )
            {
                /* if stable, do not flatten the spectrum in the 1st erased frame  */
                alpha = stab_fac * (1.0f - 2.0f*ALPHA_U) + 2.0f*ALPHA_U;  /* [0.8, 1.0] */
            }
            else if( nbLostCmpt == 2 )
            {
                alpha = ALPHA_U*1.5f;   /* 0.6 */
            }
            else
            {
                /* go rapidly to CNG spectrum  */
                alpha = ALPHA_U;
            }
        }
        else if(last_good == UNVOICED_TRANSITION )
        {
            alpha = ALPHA_UT;
        }
        else if(last_good == VOICED_CLAS || last_good == ONSET )
        {
            /* clearly voiced -  mild convergence to the CNG spectrum for the first 3 erased frames */
            alpha = ALPHA_V;
        }
        else if(last_good == SIN_ONSET)
        {
            alpha = ALPHA_S;
        }
        else
        {
            alpha = ALPHA_VT;
        }
    }
    else
    {
        alpha = 1.f/nbLostCmpt;
    }

    if( codec_mode == MODE1 )
    {
        beta = BETA_FEC;
    }
    else
    {
        if( plcBackgroundNoiseUpdated )
        {
            beta = 0.f;
        }
        else
        {
            beta = 0.25f;
        }
    }

    for( i=0; i<M; i++ )
    {
        lsf_mean[i] = beta * pt_meansForFading[i] + (1-beta)*lsf_adaptive_mean[i];
        lsf[i] = alpha * lsfold[i] + (1.0f - alpha) * lsf_mean[i];   /*   towards the CNG spectral envelope   */

        if( lsf_q_cng != NULL )
        {
            lsf_q_cng[i] = max(alpha,0.8f) * old_lsf_q_cng[i] + (1.0f - max(alpha,0.8f)) * pt_meansForFading[i];
        }
    }

    /*-----------------------------------------------------------------*
     * - estimate past quantized residual to be used in next frame
     * - Check A(z) filter stability through lsf ordering
     *-----------------------------------------------------------------*/

    if ( Opt_AMR_WB )
    {
        reorder_isf( lsf, ISF_GAP, M, INT_FS_12k8 );
    }
    else
    {
        if( L_frame == L_FRAME )
        {
            reorder_lsf( lsf, codec_mode==MODE1?MODE1_LSF_GAP:LSF_GAP, M, INT_FS_12k8 );

            if( lsf_q_cng != NULL )
            {
                reorder_lsf(lsf_q_cng, LSF_GAP, M, INT_FS_12k8);
            }
        }
        else  /* L_frame > L_FRAME */
        {
            reorder_lsf( lsf, MODE1_LSF_GAP, M, L_frame * 50 );
            if(lsf_q_cng!=NULL)
            {
                reorder_lsf(lsf_q_cng, MODE1_LSF_GAP, M, L_frame * 50);
            }
        }
    }
    /* update the AR memory to be used in the next frame */
    mvr2r( lsf, mem_AR, M );

    for(i=0; i<M; i++)
    {
        tmp = lsf[i] - pt_meansForMemUpdate[i];
        mem_MA[i] = (float)(tmp - MU_MA * mem_MA[i]); /* Update with quantized prediction error for MA model */
        mem_MA[i] *= 0.5f;                            /* Attenuate the MA Q memory                           */
    }

    return;
}

/*---------------------------------------------------------------------*
* routine:   PlcGetlsfBase()
*
*
*---------------------------------------------------------------------*/

float const * PlcGetlsfBase(
    int const lpcQuantization,
    int const narrowBand,
    int const sr_core
)
{
    if( lpcQuantization == 0 )
    {
        /* high rates, return value is never used; the correct value changes
           dynamically and is not available during PLC; therefore, the setting
           is kept as before (without the define PLC_FIX_XSF_HANDLING); the
           correct value would be isf[m] as returned by lpc_unquantize()
           during normal decoding */

        if( sr_core == 32000 )
        {
            return means_swb_cleanspeech_lsf32k0;
        }
        else if( sr_core == 25600 )
        {
            return means_swb_cleanspeech_lsf25k6;
        }
        else
        {
            return means_wb_cleanspeech_lsf16k0;
        }
    }

    /* lpcQuntization == 1 is left */

    if( sr_core == 16000 )
    {
        return GEWB2_Ave;
    }

    /* sr_core == 12.8k is left */

    if( narrowBand == 0 )
    {
        return GEWB_Ave;
    }

    /* narrowBand == 1 is left */
    return GENB_Ave;
}
