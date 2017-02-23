/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"
#include "prot.h"        /* Function prototypes                      */

/*--------------------------------------------------------------------------
 * logqnorm()
 *
 * Log quantization for norms of sub-vectors
 *--------------------------------------------------------------------------*/

void logqnorm(
    const float *x,             /* i  : coefficient vector */
    short *k,             /* o  : index */
    const short L,              /* i  : codebook length */
    const short N,              /* i  : sub-vector size */
    const float *thren          /* i  : quantization thresholds */
)
{
    short i, j, j1, j2;
    float temp, power;


    temp = 0.0;
    for (i=0; i<N; i++)
    {
        temp += (x[i] * x[i]);
    }

    /* sqrt will be done later */
    temp *= inv_tbl[N];

    if ( thren[0]*thren[0] - temp <= 0 )
    {
        *k = 0;
    }
    else if (thren[L-2]*thren[L-2] - temp > 0)
    {
        *k = L - 1;
    }
    else
    {
        power = (float)sqrt(temp);

        j1 = 0;
        j2 = L - 1;
        while ((j2-j1)>1)
        {
            j = (j1 + j2) >> 1;
            if ( power >= thren[j] )
            {
                j2 = j;
            }
            else
            {
                j1 = j;
            }
        }

        *k = j2;
    }

    return;
}


/*--------------------------------------------------------------------------
 * logqnorm_2()
 *
 *
 *--------------------------------------------------------------------------*/

void logqnorm_2(
    const float *env_fl,        /* o  : index */
    const short L,              /* i  : codebook length */
    const short n_env_band,     /* i  : sub-vector size */
    const short nb_sfm,         /* i  : sub-vector size */
    short *ynrm,
    short *normqlg2,
    const float *thren          /* i  : quantization thresholds */
)
{
    short i, j, j1, j2;
    float temp, power;

    for(i=n_env_band; i<nb_sfm; i++)
    {
        temp = env_fl[i-n_env_band];
        if ( thren[0] - temp <= 0 )
        {
            *ynrm = 0;
        }
        else if (thren[L-2] - temp > 0)
        {
            *ynrm = L - 1;
        }
        else
        {
            power = temp;
            j1 = 0;
            j2 = L - 1;
            while ((j2-j1)>1)
            {
                j = (j1 + j2) >> 1;
                if ( power >= thren[j] )
                {
                    j2 = j;
                }
                else
                {
                    j1 = j;
                }
            }
            *ynrm = j2;
        }
        *normqlg2 = dicnlg2[*ynrm];
        normqlg2++;
        ynrm++;
    }

    return;
}

/*--------------------------------------------------------------------------
 *  calc_norm()
 *
 *  Calculate the norms for the spectral envelope
 *--------------------------------------------------------------------------*/

void calc_norm(
    const float *x,                         /* i  : Input vector.                       */
    short *norm,                      /* o  : Quantization indices for norms      */
    short *normlg,                    /* o  : Quantized norms in log2             */
    const short start_band,                 /* i  : Indice of band to start coding      */
    const short num_bands,                  /* i  : Number of bands                     */
    const short *band_len,                  /* i  : Length of bands                     */
    const short *band_start                 /* i  : Start of bands                      */
)
{
    short nrm;
    short band;

    set_s(norm, 0, start_band);

    logqnorm(&x[band_start[start_band]], &nrm, 32, band_len[start_band], thren);
    norm[start_band] = nrm;
    normlg[start_band] = dicnlg2[nrm];

    for (band = start_band + 1; band < start_band + num_bands; band++)
    {
        logqnorm(&x[band_start[band]], &nrm, 40, band_len[band], thren );
        norm[band] = nrm;
        normlg[band] = dicnlg2[nrm];
    }

    return;
}
