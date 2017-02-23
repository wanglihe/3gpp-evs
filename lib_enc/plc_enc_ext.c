/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include "prot.h"
#include "typedef.h"
#include "stat_enc.h"
#include "cnst.h"
#include "rom_com.h"


#define NBITS_GACELP              5

/*-------------------------------------------------------------------*
* open_PLC_ENC_EVS()
*
*
*-------------------------------------------------------------------*/

void open_PLC_ENC_EVS(
    HANDLE_PLC_ENC_EVS hPlcExt,
    int sampleRate
)
{
    short itr;


    hPlcExt->enableGplc = 0;
    hPlcExt->calcOnlylsf = 1;
    hPlcExt->nBits = NBITS_GACELP;
    hPlcExt->stab_fac = 0;
    set_f(hPlcExt->lsfoldbfi0, 0.0f, M);
    set_f(hPlcExt->lsfoldbfi1, 0.0f, M);
    set_f(hPlcExt->lsf_adaptive_mean, 0.0f, M);
    set_f(hPlcExt->old_exc, 0.0f, 8);
    set_f(hPlcExt->mem_MA, 0.0f, M);
    set_f(hPlcExt->lsfold, 0.0f, M);
    set_f(hPlcExt->lspold, 0.0f, M);
    if( sampleRate == 12800 )
    {
        hPlcExt->T0_4th = L_SUBFR;
        hPlcExt->T0 = L_SUBFR;
        for( itr=0; itr<M; itr++ )
        {
            hPlcExt->lsf_con[itr] = lsf_init[itr];
            hPlcExt->last_lsf_ref[itr] = lsf_init[itr];
            hPlcExt->last_lsf_con[itr] = lsf_init[itr];
        }
    }
    else
    {
        hPlcExt->T0_4th = L_SUBFR;
        hPlcExt->T0 = L_SUBFR;
        for( itr=0; itr<M; itr++ )
        {
            hPlcExt->lsf_con[itr] = lsf_init[itr] * 1.25;
            hPlcExt->last_lsf_ref[itr] = lsf_init[itr] * 1.25;
            hPlcExt->last_lsf_con[itr] = lsf_init[itr] * 1.25;
        }
    }


    return;
}


/*-------------------------------------------------------------------*
* gPLC_encInfo()
*
* Function to extract and write guided information
*-------------------------------------------------------------------*/

void gPLC_encInfo(
    HANDLE_PLC_ENC_EVS self,
    const int brate,
    const int bwidth,
    const short last_clas,
    const int coder_type
)
{
    if (self)
    {
        self->calcOnlylsf = 1;
        if(bwidth >= WB  && brate == ACELP_24k40 )
        {
            self->enableGplc = 1;

            if( ( last_clas == VOICED_CLAS || last_clas == ONSET )
                    && (coder_type == VOICED || coder_type == GENERIC) )
            {
                self->nBits = NBITS_GACELP;
                self->calcOnlylsf = 0;
            }
            else
            {
                self->nBits = 1;
            }
        }
        else
        {
            self->enableGplc = 0;
            self->nBits = 0;
        }

    }

    return;
}
