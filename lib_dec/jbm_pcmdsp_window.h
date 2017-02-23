/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

/*! @file jbm_jbm_pcmdsp_window.h Window functions. */

#ifndef JBM_PCMDSP_WINDOW_H
#define JBM_PCMDSP_WINDOW_H JBM_PCMDSP_WINDOW_H

/* local headers */
#include "jbm_types.h"

/*! Generates a Hann window (cos-shaped) of length n.
 *  Roughly:
 *
 *                       1    __
 *                           /  \
 *                       0 _/    \_
 *                         <------>
 *                            n
 */
void hannWindow(uint16_t n, Float * w);

/** Overlap/Add of two signal with a given window. */
/** @param[in] fadeOut signal to fade out
 *  @param[in] fadeIn signal to fade in
 *  @param[in] out buffer to store the output signal
 *  @param[in] n number of samples
 *  @param[in] nChannels number of channels
 *  @param[in] fadeOutWin window for fade out
 *  @param[in] fadeInWin window for fade in */
void overlapAdd(const int16_t *fadeOut, const int16_t *fadeIn, int16_t *out,
                uint16_t n, uint16_t nChannels, const float *fadeOutWin, const float *fadeInWin);

#endif /* JBM_PCMDSP_WINDOW_H */
