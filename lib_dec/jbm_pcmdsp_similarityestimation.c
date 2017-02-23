/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

/*! @file jbm_pcmdsp_similarityestimation.c Algorithms for correlation and similarity estimation. */

/* system headers */
#include <stdlib.h>
#include <math.h>
#include "options.h"
/* local headers */
#include "jbm_pcmdsp_similarityestimation.h"


/* Calculates cross correlation coefficient for template segment. */
Float cross_correlation_self(const int16_t * signal,
                             uint16_t x, uint16_t y, uint16_t corr_len)
{
    Float c_c;
    int j;

    c_c = 0.0f;
    for (j = 0; j < corr_len; j++)
    {
        c_c += ((Float) signal[j + x] * (Float) signal[j + y]);
    }
    return c_c;
}

/* Calculates cross correlation coefficient for template segment. */
Float cross_correlation_subsampled_self(const int16_t * signal,
                                        uint16_t x, uint16_t y, uint16_t corr_len, uint16_t subsampling)
{
    Float c_c;
    int j;

    c_c = 0.0f;
    for (j = 0; j < corr_len; j += subsampling)
    {
        c_c += ((Float) signal[j + x] * (Float) signal[j + y]);
    }
    return c_c;
}

/* Calculates normalized cross correlation coefficient for template segment. */
Float normalized_cross_correlation_self(const int16_t * signal,
                                        uint16_t x, uint16_t y, uint16_t corr_len,
                                        uint16_t subsampling, Float * energy)
{
    Float32 c_c;
    Float32 energy_xy, energy_x, energy_y;
    uint16_t j;
    const int16_t *signal_a, *signal_b;

    c_c = 0.0f;
    energy_x = 0.0f;
    energy_y = 0.0f;
    signal_a = &signal[x];
    signal_b = &signal[y];
    for (j = 0; j < corr_len; j += subsampling)
    {
        c_c += ((Float32) signal_a[j] * (Float32) signal_b[j]);
        energy_x += ((Float32) signal_a[j]) * ((Float32) signal_a[j]);
        energy_y += ((Float32) signal_b[j]) * ((Float32) signal_b[j]);
    }
    energy_xy = (Float32)sqrt((Float32)energy_x * (Float32)energy_y);
    if(energy_xy < 1.0f)
    {
        energy_xy = 1.0f;   /* conceal silent frames */
    }

    c_c = c_c / energy_xy;
    *energy = energy_xy;
    return c_c;
}

/* Splits the signal into segments and checks if all of them have very low energy. */
bool_t isSilence(const int16_t * signal, uint32_t len, uint32_t segments)
{
    uint32_t i, samplesPerSegment;
    Float energy;

    energy = 0;
    samplesPerSegment = len / segments;
    for(i = 0; i < len; i++)
    {
        energy += (signal[i] / 32768.f) * (signal[i] / 32768.f);
        if( ( i != 0U && i % samplesPerSegment == 0U ) || i + 1 == len )
        {
            /* check energy of current segment */
            energy = 10 * (Float)log10( energy / samplesPerSegment );
            if( energy > -65 )
            {
                return false;
            }
            energy = 0;
        }
    }
    return true;
}

