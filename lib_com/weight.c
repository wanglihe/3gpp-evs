/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"

/*--------------------------------------------------------------------------*
 * sfm2mqb()
 *
 * Map sub-vectors to pbands
 *--------------------------------------------------------------------------*/

static void sfm2mqb(
    short spe[],        /* i  : sub-vectors   */
    short spe2q[],      /* o  : pbands        */
    const short nb_sfm        /* i  : number of norms */
)
{
    short tmp, i;

    /* short groups */
    spe2q[0]  = spe[0] + 3;
    spe2q[1]  = spe[1] + 3;
    spe2q[2]  = spe[2] + 3;
    spe2q[3]  = spe[3] + 3;
    spe2q[4]  = spe[4] + 3;
    spe2q[5]  = spe[5] + 3;
    spe2q[6]  = spe[6] + 3;
    spe2q[7]  = spe[7] + 3;
    spe2q[8]  = spe[8] + 3;
    spe2q[9]  = spe[9] + 3;

    spe2q[10]  = ((spe[10] + spe[11]) >> 1) + 4;
    spe2q[11]  = ((spe[12] + spe[13]) >> 1) + 4;
    spe2q[12]  = ((spe[14] + spe[15]) >> 1) + 4;

    spe2q[13]  = ((spe[16] + spe[17]) >> 1) + 5;
    spe2q[14]  = ((spe[18] + spe[19]) >> 1) + 5;

    tmp = 0;
    for (i=20; i < 24; i++)
    {
        tmp += spe[i];
    }
    spe2q[15] = (short)(((int)tmp * 8192L) >> 15) + 6;

    tmp = 0;
    for (i=24; i < 27; i++)
    {
        tmp += spe[i];
    }
    spe2q[16] = (short)(((int)tmp * 10923L) >> 15) + 6;

    if (nb_sfm > 27)
    {
        tmp = 0;
        for (i=27; i < 30; i++)
        {
            tmp += spe[i];
        }
        spe2q[17] = (short)(((int)tmp * 10923L) >> 15) + 6;

        if (nb_sfm > 30)
        {
            tmp = 0;
            for (i=30; i < 35; i++)
            {
                tmp += spe[i];
            }
            spe2q[18] = (short)(((int)tmp * 6553L) >> 15) + 7;

            tmp = 0;
            for (i=35; i < 44; i++)
            {
                tmp += spe[i];
            }
            spe2q[19] = (short)(((int)tmp * 3641L) >> 15) + 8;
        }
    }

    return;
}

/*--------------------------------------------------------------------------*
 * mqb2sfm()
 *
 * Map pbands to sub-vectors
 *--------------------------------------------------------------------------*/

static void mqb2sfm(
    short spe2q[],      /* i  : pbands         */
    short spe[],        /* o  : sub-vectors    */
    const short lnb_sfm       /* i  : number of norms */
)
{
    short i;

    spe[0]  = spe2q[0];
    spe[1]  = spe2q[1];
    spe[2]  = spe2q[2];
    spe[3]  = spe2q[3];
    spe[4]  = spe2q[4];
    spe[5]  = spe2q[5];
    spe[6]  = spe2q[6];
    spe[7]  = spe2q[7];
    spe[8]  = spe2q[8];
    spe[9]  = spe2q[9];

    spe[10] = spe2q[10];
    spe[11] = spe2q[10];

    spe[12] = spe2q[11];
    spe[13] = spe2q[11];

    spe[14] = spe2q[12];
    spe[15] = spe2q[12];

    spe[16] = spe2q[13];
    spe[17] = spe2q[13];

    spe[18] = spe2q[14];
    spe[19] = spe2q[14];
    for (i=20; i < 24; i++)
    {
        spe[i] = spe2q[15];
    }

    for (i=24; i < 27; i++)
    {
        spe[i] = spe2q[16];
    }

    if (lnb_sfm>SFM_N_STA_8k)
    {
        for (i=27; i < 30; i++)
        {
            spe[i] = spe2q[17];
        }

        if (lnb_sfm>SFM_N_STA_10k)
        {
            for (i=30; i < 35; i++)
            {
                spe[i] = spe2q[18];
            }

            for (i=35; i < 44; i++)
            {
                spe[i] = spe2q[19];
            }
        }
    }

    return;
}

/*--------------------------------------------------------------------------*
 * map_quant_weight()
 *
 * Calculate the quantization weights
 *--------------------------------------------------------------------------*/

void map_quant_weight(
    const short normqlg2[],            /* i  : quantized norms   */
    short wnorm[],               /* o  : weighted norm     */
    const short is_transient           /* i  : transient flag    */
)
{
    short sfm;
    short tmp16;
    short spe2q[NUM_MAP_BANDS];
    short spe[NB_SFM];

    short spe2q_max;
    short spe2q_min;
    short norm_max;
    short shift;
    short sum;
    short k;
    short lnb_sfm,num_map_bands;

    if (is_transient)
    {
        lnb_sfm = NB_SFM;
        num_map_bands = NUM_MAP_BANDS;

        for (sfm = 0; sfm < lnb_sfm; sfm+=4)
        {
            sum = 0;
            for (k=0; k < 4; k++)
            {
                sum = sum + normqlg2[sfm+k];
            }
            sum = sum >> 2;
            for (k=0; k < 4; k++)
            {
                spe[sfm +k] = sum;
            }
        }
    }
    else
    {
        lnb_sfm = NB_SFM;
        num_map_bands = NUM_MAP_BANDS;


        for (sfm = 0; sfm < lnb_sfm; sfm++)
        {
            spe[sfm] = normqlg2[sfm];
        }
    }

    sfm2mqb(spe, spe2q, lnb_sfm);

    for (sfm = 0; sfm < num_map_bands; sfm++)
    {
        spe2q[sfm] = spe2q[sfm] - 10;
    }

    /* spectral smoothing */
    for (sfm = 1; sfm < num_map_bands; sfm++)
    {
        tmp16 = spe2q[sfm-1] - 4;
        if (spe2q[sfm]<tmp16)
        {
            spe2q[sfm] = tmp16;
        }
    }

    for (sfm = num_map_bands-2 ; sfm >= 0 ; sfm--)
    {
        tmp16 = spe2q[sfm+1] - 8;
        if (spe2q[sfm]<tmp16)
        {
            spe2q[sfm] = tmp16;
        }
    }

    for (sfm = 0; sfm < num_map_bands ; sfm++)
    {
        if (spe2q[sfm]<a_map[sfm])
        {
            spe2q[sfm] = a_map[sfm];
        }
    }

    /* Saturate by the Absolute Threshold of Hearing */
    spe2q_max = MIN16B;
    spe2q_min = MAX16B;

    for (sfm = 0; sfm < num_map_bands ; sfm++)
    {
        spe2q[sfm] = sfm_width[sfm] - spe2q[sfm];

        if (spe2q_max<spe2q[sfm])
        {
            spe2q_max = spe2q[sfm];
        }

        if (spe2q_min>spe2q[sfm])
        {
            spe2q_min = spe2q[sfm];
        }
    }

    for (sfm = 0; sfm < num_map_bands; sfm++)
    {
        spe2q[sfm] = spe2q[sfm] - spe2q_min;
    }

    spe2q_max = spe2q_max - spe2q_min;

    if (spe2q_max == 0)
    {
        norm_max = 0;
    }
    else
    {
        if (spe2q_max < 0)
        {
            spe2q_max = ~spe2q_max;
        }
        for (norm_max=0; spe2q_max<0x4000; norm_max++)
        {
            spe2q_max <<= 1;
        }
    }

    shift = norm_max - 13;
    for (sfm = 0; sfm < num_map_bands ; sfm++)
    {
        if (shift<0)
        {
            spe2q[sfm] = spe2q[sfm] >> (-shift);
        }
        else
        {
            spe2q[sfm] = spe2q[sfm] << shift;
        }
    }

    mqb2sfm(spe2q,spe,lnb_sfm);

    if (is_transient)
    {
        for (sfm = 0; sfm < lnb_sfm; sfm+=4)
        {
            sum = 0;
            for (k=0; k < 4; k++)
            {
                sum = sum + spe[sfm+k];
            }

            sum = sum >> 2;

            for (k=0; k < 4; k++)
            {
                spe[sfm +k] = sum;
            }
        }
    }

    /* modify the norms for bit-allocation */
    for (sfm = 0; sfm < lnb_sfm ; sfm++)
    {
        wnorm[sfm] = spe[sfm] + normqlg2[sfm];
    }

    return;
}
