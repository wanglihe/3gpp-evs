/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <memory.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include "typedef.h"
#include "prot.h"
#include "rom_enc.h"
#include "rom_com.h"




/*
 * E_GAIN_norm_corr_interpolate
 *
 * Parameters:
 *    x           I: input vector
 *    frac        I: fraction (-4..+3)
 *
 * Function:
 *    Interpolating the normalized correlation
 *
 * Returns:
 *    interpolated value
 */
static Float32 E_GAIN_norm_corr_interpolate(Float32 *x, Word32 frac)
{
    Float32 s, *x1, *x2;
    const Float32 *c1, *c2;

    if (frac < 0)
    {
        frac += 4;
        x--;
    }

    x1 = &x[0];
    x2 = &x[1];
    c1 = &E_ROM_inter4_1[frac];
    c2 = &E_ROM_inter4_1[4 - frac];
    s = x1[0] * c1[0] + x2[0] * c2[0];
    s += x1[-1] * c1[4] + x2[1] * c2[4];
    s += x1[-2] * c1[8] + x2[2] * c2[8];
    s += x1[-3] * c1[12] + x2[3] * c2[12];
    return s;
}

static Float32 E_GAIN_norm_corr_interpolate6(Float32 *x, Word32 frac)
{
    Float32 s, *x1, *x2;
    const Float32 *c1, *c2;

    if (frac < 0)
    {
        frac += 6;
        x--;
    }

    x1 = &x[0];
    x2 = &x[1];
    c1 = &E_ROM_inter6_1[frac];
    c2 = &E_ROM_inter6_1[6 - frac];
    s = x1[0] * c1[0] + x2[0] * c2[0];
    s += x1[-1] * c1[6] + x2[1] * c2[6];
    s += x1[-2] * c1[12] + x2[2] * c2[12];
    s += x1[-3] * c1[18] + x2[3] * c2[18];
    return s;
}

/*
 * E_GAIN_closed_loop_search
 *
 * Parameters:
 *    exc            I: excitation buffer
 *    xn             I: target signal
 *    h              I: weighted synthesis filter impulse response
 *    dn             I: residual domain target signal
 *    t0_min         I: minimum value in the searched range
 *    t0_max         I: maximum value in the searched range
 *    pit_frac       O: chosen fraction
 *    i_subfr        I: flag to first subframe
 *    t0_fr2         I: minimum value for resolution 1/2
 *    t0_fr1         I: minimum value for resolution 1
 *
 * Function:
 *    Find the closed loop pitch period with 1/4 subsample resolution.
 *
 * Returns:
 *    chosen integer pitch lag
 */
Word32 E_GAIN_closed_loop_search(Float32 exc[],
                                 Float32 xn[], Float32 h[],
                                 Word32 t0_min, Word32 t0_min_frac, Word32 t0_max, Word32 t0_max_frac, Word32 t0_min_max_res, Word32 *pit_frac, Word32 *pit_res, Word32 pit_res_max,
                                 Word32 i_subfr, Word32 pit_min, Word32 pit_fr2, Word32 pit_fr1, Word32 L_subfr)
{
    Float32 corr_v[32 + 2 * L_INTERPOL1 + 1];
    Float32 cor_max, max, temp;
    Float32 *corr;
    Word32 i, fraction, frac1, frac2, step;
    Word32 t0, t_min, t_max;
    (void)t0_min_frac;
    (void)t0_max_frac;
    (void)t0_min_max_res;

    /* Find interval to compute normalized correlation */
    if (t0_min_frac>0)
    {
        t0_min++;
    }
    t_min = t0_min - L_INTERPOL1;
    t_max = t0_max + L_INTERPOL1;

    /* allocate memory to normalized correlation vector */
    corr = &corr_v[-t_min];      /* corr[t_min..t_max] */

    /* Compute normalized correlation between target and filtered excitation */
    norm_corr(exc, xn, h, t_min, t_max, corr, L_subfr);

    /*  find integer pitch */
    max = corr[t0_min];
    t0  = t0_min;
    for(i = t0_min + 1; i <= t0_max; i++)
    {
        if( corr[i] >= max)
        {
            max = corr[i];
            t0 = i;
        }
    }



    /* If first subframe and t0 >= pit_fr1, do not search fractionnal pitch */
    if((i_subfr == 0) & (t0 >= pit_fr1))
    {
        *pit_frac = 0;
        *pit_res = 1;
        return(t0);
    }


    /*
     * Search fractionnal pitch
     * Test the fractions around t0 and choose the one which maximizes
     * the interpolated normalized correlation.
     */

    if ( t0_min_max_res == (pit_res_max>>1) )
    {
        t0_min_frac = t0_min_frac << 1;
        t0_max_frac = t0_max_frac << 1;
    }

    step = 1;
    frac1 = -(pit_res_max-1);
    frac2 = pit_res_max-1;
    if (((i_subfr == 0) & (t0 >= pit_fr2)) | (pit_fr2 <= pit_min))
    {
        step = 2;
        frac1 = -(pit_res_max-2);
        frac2 = pit_res_max-2;
    }

    if ( (t0 == t0_min) && (t0_min_frac==0) )
    {
        frac1 = t0_min_frac;
    }
    else if ( (t0 == t0_min) && (frac1+pit_res_max<t0_min_frac) )
    {
        frac1 = t0_min_frac-pit_res_max;
    }
    if (t0 == t0_max)
    {
        frac2 = t0_max_frac;
    }
    assert(frac1<=0 && frac2>=0 && frac2>frac1);
    if (pit_res_max == 6)
    {
        cor_max = E_GAIN_norm_corr_interpolate6(&corr[t0], frac1);
        fraction = frac1;
        for (i = (frac1 + step); i <= frac2; i += step)
        {
            ;
            temp = E_GAIN_norm_corr_interpolate6(&corr[t0], i);
            if (temp > cor_max)
            {
                cor_max = temp;
                fraction = i;
            }

        }
    }
    else
    {
        cor_max = E_GAIN_norm_corr_interpolate(&corr[t0], frac1);
        fraction = frac1;
        for (i = (frac1 + step); i <= frac2; i += step)
        {
            ;
            temp = E_GAIN_norm_corr_interpolate(&corr[t0], i);
            if (temp > cor_max)
            {
                cor_max = temp;
                fraction = i;
            }

        }
    }

    /* limit the fraction value */
    if (fraction < 0)
    {
        fraction += pit_res_max;
        t0 -= 1;
    }
    if (((i_subfr == 0) & (t0 >= pit_fr2)) | (pit_fr2 <= pit_min))
    {
        *pit_res = pit_res_max>>1;
        *pit_frac = fraction>>1;
    }
    else
    {
        *pit_res = pit_res_max;
        *pit_frac = fraction;
    }
    return (t0);
}



