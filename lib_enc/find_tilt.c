/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"

/*---------------------------------------------------------------------*
 * Local constants
 *---------------------------------------------------------------------*/

#define FACT    3.0f                            /* background noise energy estimation adjusting factor - to maintain the ADR about the same */
#define TH_COR  0.6f                            /* Minimum correlation for per bin processing                    */
#define TH_D    50.0f                           /* Difference limit between nearest harmonic and a frequency bin */
#define TH_PIT  (INT_FS_12k8 / (2.0f * TH_D))        /* Maximum pitch for per bin processing                          */

/*-------------------------------------------------------------------*
 * find_tilt()
 *
 * Find LF/HF energy ratio
 *-------------------------------------------------------------------*/

void find_tilt(
    const float fr_bands[],  /* i  : energy in frequency bands                  */
    const float bckr[],      /* i  : per band background noise energy estimate  */
    float ee[2],       /* o  : lf/hf E ration for present frame           */
    const short pitch[3],    /* i  : open loop pitch values for 3 half-frames   */
    const float voicing[3],  /* i  : normalized correlation for 3 half-frames   */
    const float *lf_E,       /* i  : per bin energy  for low frequencies        */
    const float corr_shift,  /* i  : normalized correlation correction          */
    const short bwidth,      /* i  : input signal bandwidth                     */
    const short max_band,    /* i  : maximum critical band                      */
    float hp_E[],      /* o  : energy in HF                               */
    const short codec_mode,  /* i  : Mode 1 or 2                                */
    float *bckr_tilt_lt,               /* i/o: lf/hf E ratio of background noise          */
    short Opt_vbr_mode
)
{
    float lp_bckr, hp_bckr, lp_E, freq, f0, f1, f2, mean_voi, bin;
    const float *pt_bands, *pt_bckr, *tmp_E, *hf_bands, *pt_E;
    short cnt, i, nb_bands;
    float th_pit;

    /*-----------------------------------------------------------------*
     * Initializations
     *-----------------------------------------------------------------*/

    th_pit = TH_PIT;

    if( bwidth != NB )
    {
        /* WB processing */
        bin = BIN;                  /* First useful frequency bin ~ 50 Hz     */
        pt_bands = fr_bands;
        tmp_E = lf_E;
        pt_bckr = bckr;
        nb_bands = 10;
    }
    else
    {
        /* NB processing */
        bin = 3.0f*BIN;             /* first useful frequency bin ~ 150 Hz    */
        pt_bands = fr_bands+1;      /* exlcude 1st critical band              */
        tmp_E = lf_E + 2;           /* start at the 3rd bin (150 Hz)          */
        pt_bckr = bckr+1;           /* exclude 1st critical band              */
        nb_bands = 9;               /* nb. of "low" frequency bands taken into account in NB processing      */
    }

    /*-----------------------------------------------------------------*
     * Find spectrum tilt
     *-----------------------------------------------------------------*/

    pt_E = tmp_E;                   /* pointer at the 1st useful element of the per-bin energy vector  */
    hf_bands = fr_bands;

    /* bckr + voicing */
    lp_bckr = mean( pt_bckr, nb_bands );                              /* estimated noise E in first critical bands, up to 1270 Hz */
    hp_bckr = 0.5f * (bckr[max_band-1] + bckr[max_band]);             /* estimated noise E in last 2 critical bands */
    *bckr_tilt_lt = 0.9f * *bckr_tilt_lt + 0.1f * lp_bckr / hp_bckr;

    if ( codec_mode == MODE2
            || Opt_vbr_mode
       )
    {
        lp_bckr *= FACT;
        hp_bckr *= FACT;
    }

    mean_voi = 0.5f * (voicing[1] + voicing[2]) + corr_shift;
    f0 = INT_FS_12k8 / pitch[2];

    for( i=0; i<2; i++ )
    {
        hp_E[i] = 0.5f * (hf_bands[max_band-1] + hf_bands[max_band]) - hp_bckr; /* average E in last 2 critical bands */

        if ( !Opt_vbr_mode)
        {
            if( hp_E[i] < E_MIN )
            {
                /* to avoid division by 0 */
                hp_E[i] = E_MIN;
            }

        }
        else
        {
            if( hp_E[i] < 1.0f )
            {
                /* to avoid division by 0 */
                hp_E[i] = 1.0f;
            }
        }


        if( (mean_voi > TH_COR) && (pitch[2] < th_pit) )
        {
            /* high-pitched voiced frames */
            freq = bin;                           /* 1st useful frequency bin */
            cnt = 0;
            lp_E = 0.0f;
            f1 = 1.5f * f0;                       /* middle between 2 harmonics */
            f2 = f0;

            while( freq <= 1270.0f )              /* end frequency of 10th critical band */
            {
                /*pt_E*/
                while( freq <= f1 )
                {
                    if( fabs(freq-f2) < TH_D )    /* include only bins sufficiently close to harmonics */
                    {
                        lp_E += *pt_E;
                        cnt++;
                    }
                    freq += BIN;
                    pt_E++;
                }
                f1 += f0;
                f2 += f0;                         /* next harmonic */
            }

            lp_E = lp_E / (float)cnt - lp_bckr;
            pt_E = tmp_E + VOIC_BINS;             /* update for next half-frame */
        }
        else
        {
            /* other than high-pitched voiced frames */
            lp_E = mean( pt_bands, nb_bands ) - lp_bckr; /* average E in first critical bands, up to 1270 Hz */
        }

        if ( !Opt_vbr_mode)
        {
            if( lp_E < E_MIN )
            {
                /* avoid negative E due to noise subtraction */
                lp_E = E_MIN;
            }
        }
        else
        {

            if( lp_E < 0.0f )
            {
                /* avoid negative E due to noise subtraction */
                lp_E = 0.0f;
            }
        }

        /* calculate the tilt (LF/HF ratio) */
        ee[i] = lp_E / hp_E[i];

        if( bwidth == NB )                        /* for NB input, compensate for the missing bands */
        {
            ee[i] *= 6.0f;
        }

        pt_bands += NB_BANDS;                     /* update pointers for the next half-frame */
        hf_bands += NB_BANDS;
    }

    return;
}
