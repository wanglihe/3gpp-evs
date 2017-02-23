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

#define BIN_4000       80    /* The frequency bin corresponding to 4kHz */
#define MAX_SNR1       45.0f
#define INV_MAX_SNR    (1.0f / (MAX_SNR1-1.0f)) /* Max. SNR considered for noise subtraction in voiced segments */
#define MAX_SNR_SNR1   (MAX_SNR1 * INV_MAX_SNR) /* 45 * (1 / (MAX_SNR1-1)) */
#define MAX_BANDEXC    20
#define TH_0_MAX       (1.5f*3.125f)
#define TH_1_MAX       (1.5f*2.8125f)
#define TH_2_MAX       (1.5f*2.1875f)
#define TH_3_MAX       (1.5f*1.875f)

#define TH_UP          0.15625f
#define TH_DW          0.15625f
#define NB_TH3_MIN     30
#define NB_TH1_MIN     30

#define MAX_THR        0.92f

/*------------------------------------------------------------------------*
 * stab_est()
 *
 * Signal stability estimation based on energy variation
 *------------------------------------------------------------------------*/

short stab_est(
    float etot,                /* i   : Total energy of the current frame   */
    float *lt_diff_etot,       /* i/o : Long term total energy variation    */
    float *mem_etot,           /* i/o : Total energy memory                 */
    short *nb_thr_3,           /* i/o : Number of consecutives frames of level 3 */
    short *nb_thr_1,           /* i/o : Number of consecutives frames of level 1 */
    float *thresh,             /* i/o : Detection thresold                 */
    short *last_music_flag,    /* i/o : Previous music detection ouptut    */
    short vad_flag
)
{
    short i, music_flag2;
    float mean_diff;
    float ftmp_c, fcorr, dev;

    /*------------------------------------------------------------------------*
     * Find mean of the past MAX_LT frames energy variation
     *------------------------------------------------------------------------*/

    mean_diff = 0.0f;
    for (i = 1; i<MAX_LT; i++)
    {
        mean_diff += lt_diff_etot[i-1] * INV_MAX_LT; /* divide by MAX_LT */
        lt_diff_etot[i-1] = lt_diff_etot[i];
    }

    mean_diff += lt_diff_etot[i-1] * INV_MAX_LT; /* divide by MAX_LT */

    /*------------------------------------------------------------------------*
     * Find statistical deviation of the energy variation history
     * against the last 15 frames
     *------------------------------------------------------------------------*/

    fcorr = 0.0f;
    for (i = MAX_LT-15; i<MAX_LT; i++)
    {
        ftmp_c = lt_diff_etot[i] - mean_diff;
        fcorr += ftmp_c * ftmp_c;
    }

    lt_diff_etot[i-1] = etot - *mem_etot;
    *mem_etot = etot;

    /*------------------------------------------------------------------------*
     * Compute statistical deviation
     *------------------------------------------------------------------------*/

    dev = (float)sqrt(fcorr / (MAX_LT-15));

    /*------------------------------------------------------------------------*
     * State machine to decide level of inter-harmonic noise reduction and
     * the first bins where this inter-harmonic noise reduction will be applied
     * (only if this frame is GOOD or if we are already far from NB_BFI_THR)
     * (if music_flag2 is 0, the spectral modification is deactivated, otherwise, it is activated)
     *------------------------------------------------------------------------*/

    music_flag2 = 0;

    /*--------------------------------------------------------------------*
     * statistical deviation < Thresh3 and last signal category type >= 3
     * (last category was "tonal" and the new one is "very tonal")
     *--------------------------------------------------------------------*/

    if ( dev < thresh[3] && *last_music_flag >= 3 )
    {
        music_flag2 = 4;
        *nb_thr_3 += 1;
        *nb_thr_1 = 0;
    }

    /*--------------------------------------------------------------------*
     * statistical deviation < Thresh2 and last signal category type >= 2
     * (last category was "moderatly tonal" and the new one is a "tonal" )
     *--------------------------------------------------------------------*/

    else if ( dev < thresh[2] && *last_music_flag >= 2 )
    {
        music_flag2 = 3;
        *nb_thr_3 += 1;
        *nb_thr_1 = 0;
    }

    /*--------------------------------------------------------------------*
     * statistical deviation < Thresh1 and last signal category type  >= 1
     * (last category was "slightly tonal" and the new one is a "moderatly tonal")
     *--------------------------------------------------------------------*/

    else if ( dev < thresh[1] && *last_music_flag >= 1 )
    {
        music_flag2 = 2;
    }

    /*--------------------------------------------------------------------*
     * statistical deviation < Thresh0
     * (last category was "not tonal" and the new one is "slightly tonal")
     *--------------------------------------------------------------------*/

    else if ( dev < thresh[0] )
    {
        music_flag2 = 1;
    }

    /*--------------------------------------------------------------------*
     * statistical deviation > Thresh0
     * (Statistical deviation is high: the new tonal category is not tonal)
     *--------------------------------------------------------------------*/

    else
    {
        *nb_thr_1 += 1;
        *nb_thr_3 = 0;
    }

    /*------------------------------------------------------------------------*
     * Update the thresholds
     *------------------------------------------------------------------------*/

    if (*nb_thr_3 > NB_TH3_MIN)
    {
        /* the number of consecutive categories type 3 or 4 (most tonal and tonal) */
        /* is greater than 30 frames ->increase the deviations thresholds to allow more variation */
        thresh[0] += TH_UP;
        thresh[1] += TH_UP;
        thresh[2] += TH_UP;
        thresh[3] += TH_UP;
    }
    else if (*nb_thr_1 > NB_TH1_MIN)
    {
        /* the number of consecutive categories type 0 (non tonal frames) */
        /* is greater than 30 frames -> decrease the deviations thresholds to allow less variation */
        thresh[0] -= TH_DW;
        thresh[1] -= TH_DW;
        thresh[2] -= TH_DW;
        thresh[3] -= TH_DW;
    }

    /* limitation of the threshold (this local macro stores the highest of the two and it also counts the # of operations) */
    set_max( &thresh[0], TH_0_MIN2 );
    set_max( &thresh[1], TH_1_MIN2 );
    set_max( &thresh[2], TH_2_MIN2 );

    set_min( &thresh[0], TH_0_MAX );
    set_min( &thresh[1], TH_1_MAX );
    set_min( &thresh[2], TH_2_MAX );

    set_max( &thresh[3], TH_3_MIN2 );
    set_min( &thresh[3], TH_3_MAX );
    /*------------------------------------------------------------------------*
     * Final updates
     *------------------------------------------------------------------------*/

    *last_music_flag = music_flag2;
    if ( vad_flag == 0 )
    {
        /* overwrite decision in unvoiced frames */
        music_flag2 = 0;
    }

    return music_flag2;
}

