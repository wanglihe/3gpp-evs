/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

/*! @file jbm_pcmdsp_similarityestimation.h Algorithms for correlation and similarity estimation. */

#ifndef JBM_PCMDSP_SIMILARITYESTIMATION_H
#define JBM_PCMDSP_SIMILARITYESTIMATION_H JBM_PCMDSP_SIMILARITYESTIMATION_H

/* local headers */
#include "jbm_types.h"

/*
********************************************************************************
*
*     Function        : cross_correlation_self
*     Tables          : <none>
*     Compile Defines : <none>
*     Return          : (Float) cross correlation coefficient
*     Information     : Calculate cross correlation coefficient for template
*                       segment.
*                       The returned value is signal-energy dependant.
*
*                       Used formula:
*
*                       corr_len-1
*                       ----
*                       \
*                       /    (j+x)*(j+y)
*                       ----
*                       j=0
*
*
*     23-JUL-04  S.Doehla        initial version
*
********************************************************************************
*/
Float cross_correlation_self(const int16_t * signal,
                             uint16_t x,
                             uint16_t y,
                             uint16_t corr_len);

/*
********************************************************************************
*
*     Function        : cross_correlation_subsampled_self
*     Tables          : <none>
*     Compile Defines : <none>
*     Return          : (Float) cross correlation coefficient
*     Information     : Calculate cross correlation coefficient for template
*                       segment.
*                       The returned value is signal-energy dependant.
*
*                       Used formula:
*
*                       corr_len-1
*                       ----
*                       \
*                       /    (j+x)*(j+y)
*                       ----
*                       j=0
*
*
*     23-JUL-04  S.Doehla        initial version
*
********************************************************************************
*/
Float cross_correlation_subsampled_self(const int16_t * signal,
                                        uint16_t x,
                                        uint16_t y,
                                        uint16_t corr_len,
                                        uint16_t subsampling);

/*
********************************************************************************
*
*     Function        : normalized_cross_correlation_self
*     Tables          : <none>
*     Compile Defines : <none>
*     Return          : (Float) normalized cross correlation coefficient
*     Information     : Calculate normalized cross correlation coefficient
*                       for template segment.
*                       The returned value is signal-energy independant.
*                       This means, no matter how loud your signal is, equal
*                       signals will return 1.0, cross-phased signals -1.0.
*
*                       Complexity is very high due to many floating point
*                       operations and using squared root!
*
*                       This function fills parameter energy with the common
*                       energy of signal x and signal y. This might be useful
*                       for silence detection.
*
*                       Used formula:
*
*                       corr_len-1
*                       ----
*                       \         (j+x)*(j+y)
*                        \    __________________
*                        /      --------------
*                       /     -/ (j+x)�+(j+y)�
*                       ----
*                       j=0
*
*
*     23-JUL-04  S.Doehla        initial version
*
********************************************************************************
*/
Float normalized_cross_correlation_self(const int16_t * signal,
                                        uint16_t x,
                                        uint16_t y,
                                        uint16_t corr_len,
                                        uint16_t subsampling,
                                        Float * energy);

/* Splits the signal into segments and checks if all of them have very low energy. */
bool_t isSilence(const int16_t * signal, uint32_t len, uint32_t segments);

#endif /* JBM_PCMDSP_SIMILARITYESTIMATION_H */
