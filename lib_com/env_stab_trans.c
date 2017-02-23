/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"
#include <math.h>


/*--------------------------------------------------------------------------*
 * env_stab_transient_detect()
 *
 * Transient detector for envelope stability measure
 *--------------------------------------------------------------------------*/

void env_stab_transient_detect(
    const short is_transient,         /* i:   Transient flag                                        */
    const short length,               /* i  : Length of spectrum (32 or 48 kHz)                     */
    const short norm[],               /* i  : quantization indices for norms                        */
    short *no_att_hangover,     /* i/o: Frame counter for attenuation hangover                */
    float *energy_lt,           /* i/o: Long-term energy measure for transient detection      */
    const short HQ_mode,              /* i  : HQ coding mode                                        */
    const short bin_th,               /* i  : HVQ cross-over frequency bin                          */
    const float *coeff                /* i  : Coded spectral coefficients                           */
)
{
    float d_max;
    float e_frame;
    short blk;
    short i;
    float E_sub[4];
    float delta_e_sub;
    short norm_ind;

    short num_subframes = 4;
    short bands_per_subframe = 9;

    if( HQ_mode == HQ_HVQ )
    {
        e_frame = 0.0f;

        for (i = 0; i < bin_th; i++)
        {
            e_frame += coeff[i]*coeff[i];
        }

        e_frame = (float)sqrt(e_frame / bin_th);

        if (e_frame > ENERGY_TH)
        {
            *energy_lt = ENERGY_LT_BETA*(*energy_lt) + (1-ENERGY_LT_BETA)*e_frame;
        }

        if (*no_att_hangover > 0)
        {
            (*no_att_hangover)--;
        }
    }
    else
    {
        d_max = 0.0f;
        e_frame = 0.0f;
        if (is_transient && length == L_FRAME32k)
        {
            /* Measure subframe energies */
            for (blk = 0; blk < num_subframes; blk++)
            {
                E_sub[blk] = 0.0f;
                for (i=0; i<bands_per_subframe; i++)
                {
                    norm_ind = subf_norm_groups[blk][i];
                    E_sub[blk] += dicn[norm[norm_ind]];
                }
                E_sub[blk] = E_sub[blk] / bands_per_subframe;
                e_frame +=  E_sub[blk];
            }
            /* Test for transient */
            if (e_frame > ENERGY_TH * num_subframes)
            {
                for (blk = 0; blk < num_subframes-1; blk++)
                {
                    delta_e_sub = (E_sub[blk+1]-E_sub[blk]) / *energy_lt;
                    if (delta_e_sub > d_max)
                    {
                        d_max = delta_e_sub;
                    }
                }
            }
        }
        else
        {
            /* Update long-term energy measure */
            e_frame = 0.0f;

            for (i = 0; i < SFM_N_ENV_STAB; i++)
            {
                e_frame += dicn[norm[i]];
            }

            e_frame = e_frame / SFM_N_ENV_STAB;

            if (e_frame > ENERGY_TH)
            {
                *energy_lt = ENERGY_LT_BETA*(*energy_lt) + (1-ENERGY_LT_BETA)*e_frame;
            }
        }

        /* Add hang-over for conservative application of stability-dependent attenuation */
        if(d_max > DELTA_TH)
        {
            *no_att_hangover = ATT_LIM_HANGOVER;
        }
        else if (*no_att_hangover > 0)
        {
            (*no_att_hangover)--;
        }

    }

    return;
}
