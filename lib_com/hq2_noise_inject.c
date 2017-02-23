/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "prot.h"

/*--------------------------------------------------------------------------*
 * hq2_noise_inject()
 *
 * HQ2 noise injection for WB signals
 *--------------------------------------------------------------------------*/

void hq2_noise_inject(
    float y2hat[],
    const short band_start[],
    const short band_end[],
    const short band_width[],
    float Ep[],
    float Rk[],
    const int   npulses[],
    short ni_seed,
    const short bands,
    const short ni_start_band,
    const short bw_low,
    const short bw_high,
    const float enerL,
    const float enerH,
    float last_ni_gain[],
    float last_env[],
    short *last_max_pos_pulse,
    short *p2a_flags,
    short p2a_bands,
    const short hqswb_clas,
    const short bwidth,
    const long  bwe_br
)
{
    short i, j, k, ni_end_band, satur = 0, count[BANDS_MAX], max_pos_pulse, pos;
    float ni_gain[BANDS_MAX], pd[BANDS_MAX], rand, peak[BANDS_MAX], fac, env[BANDS_MAX];
    short sb = bands;

    if( (hqswb_clas == HQ_HARMONIC || hqswb_clas ==  HQ_NORMAL) && (bwe_br == HQ_16k40 || bwe_br == HQ_13k20) && bwidth == SWB )
    {
        sb = (bwe_br == HQ_16k40) ? 19 : 17;
    }

    /* calculate the envelopes/ the decoded peak coeff./number of the decoded coeff./ the last subbands of the bit-allocated/saturation of bit-allocation */
    ni_end_band = bands;
    max_pos_pulse = bands;
    for (k = ni_start_band; k < ni_end_band; k++)
    {
        pd[k] = (float) Rk[k] / band_width[k];
        env[k] = (float)sqrt(Ep[k]/band_width[k]);

        peak[k] = 0.0f;
        count[k] = 0;
        if(npulses[k] != 0)
        {
            for (i = band_start[k]; i <= band_end[k]; i++)
            {
                Ep[k] -= y2hat[i] * y2hat[i];
                if( fabs(y2hat[i]) > peak[k])
                {
                    peak[k] = (float)fabs(y2hat[i]);
                }

                if(y2hat[i] != 0)
                {
                    count[k]++;
                }
            }

            max_pos_pulse = k;
            Ep[k] = Ep[k] > 0 ? (float)sqrt(Ep[k]/band_width[k]) : 0.0f;
        }
        else
        {
            Ep[k] = env[k];
        }
    }

    for (k = ni_start_band; k < ni_end_band; k++)
    {
        /* calculate the noise gain */
        satur = (pd[k] >= 0.8f) ? 1 : 0;
        if (satur == 0 && Ep[k] > 0)
        {
            if(npulses[k] != 0)
            {
                if( bwidth == SWB )
                {
                    if(hqswb_clas != HQ_TRANSIENT)
                    {
                        if(hqswb_clas ==  HQ_HARMONIC)
                        {
                            fac = (k <= sb) ? 6.0f*(1.5f - pd[k])*env[k]*Ep[k]/(peak[k]*peak[k]) : 1.5f*Ep[k]/peak[k];
                        }
                        else
                        {
                            fac = (k <= sb) ? 5.0f*(1.5f - pd[k])*Ep[k]/(peak[k]) : 4.0f*Ep[k]/peak[k];
                        }
                    }
                    else
                    {
                        fac = 1.1f;
                    }
                }
                else
                {
                    fac = 20.0f*min(1.0f, (1.5f - pd[k]))*env[k]*Ep[k]/(peak[k]*peak[k]);
                    if(k > 1 && k < ni_end_band-1)
                    {
                        if(count[k+1] == 0 && env[k] > 0.5f*env[k+1] && env[k] < 2.0f*env[k-1])
                        {
                            fac = 1.5f*env[k+1]/env[k];
                        }
                        else if(count[k-1] == 0 && peak[k] > 2.0f*env[k])
                        {
                            fac = env[k-1]/env[k];
                        }
                    }

                    if(k >= ni_end_band - p2a_bands && bwidth == WB)
                    {
                        if(bw_low*enerH > bw_high*enerL && peak[k] < 2.0f*env[k])
                        {
                            fac *= (2.0f - Ep[k]/peak[k]);
                        }

                        if( p2a_flags[k] == 0 && fac > 0 )
                        {
                            fac *= min(1.25f*env[k]/(Ep[k]*fac), 1.5f);
                        }
                    }
                }
            }
            else
            {
                fac = (hqswb_clas == HQ_HARMONIC && bwidth == SWB) ? 0.8f : 1.1f;
            }

            ni_gain[k] = fac * Ep[k];
        }
        else
        {
            ni_gain[k] = 0.0;
        }

        /* smooth the noise gain between the current frame and the previous frame */
        pos = bwidth == SWB ? ni_end_band-1 : max(max_pos_pulse, *last_max_pos_pulse);

        if(k <= pos)
        {
            if(k > 0 && k < ni_end_band-1)
            {
                if( (env[k] > 0.5f*last_env[k] && env[k] < 2.0f*last_env[k]) ||
                        (((env[k]+env[k-1]+env[k+1]) > 0.5f*(last_env[k]+last_env[k-1]+last_env[k+1])) &&
                         ((env[k]+env[k-1]+env[k+1]) < 2.0f*(last_env[k]+last_env[k-1]+last_env[k+1]))) )
                {
                    if(ni_gain[k] > last_ni_gain[k])
                    {
                        ni_gain[k] = 0.2f*ni_gain[k] + 0.8f*last_ni_gain[k];
                    }
                    else
                    {
                        ni_gain[k] = 0.6f*ni_gain[k] + 0.4f*last_ni_gain[k];
                    }
                }
            }
            else if (k == ni_end_band-1)
            {
                if( (env[k] > 0.5f*last_env[k] && env[k] < 2.0f*last_env[k]) ||
                        ((env[k]+env[k-1]) > 0.5f*(last_env[k]+last_env[k-1]) &&
                         (env[k]+env[k-1]) < 2.0f*(last_env[k]+last_env[k-1])) )
                {
                    if(ni_gain[k] > last_ni_gain[k])
                    {
                        ni_gain[k] = 0.2f*ni_gain[k] + 0.8f*last_ni_gain[k];
                    }
                    else
                    {
                        ni_gain[k] = 0.6f*ni_gain[k] + 0.4f*last_ni_gain[k];
                    }
                }
            }
        }

        /* inject noise into the non-decoded coeffs */
        if(k >= ni_end_band - p2a_bands && p2a_flags[k] == 0 && bwidth != SWB)
        {
            for (i = band_start[k]; i <= band_end[k]; i++)
            {
                if (y2hat[i] != 0)
                {
                    y2hat[i] *= 0.8f;
                }
            }
        }

        if(k == max_pos_pulse && k < bands - p2a_bands && satur != 1 && bwidth != SWB)
        {
            j = 0;
            if(ni_gain[k] < 0.01f)
            {
                fac = 0.0f;
            }
            else
            {
                fac = max(Ep[k]/ni_gain[k], 1.0f);
            }
            for (i = band_start[k]; i <= band_end[k]; i++)
            {
                if (y2hat[i] == 0)
                {
                    rand = own_random (&ni_seed) / 32768.0f;
                    y2hat[i] += (fac - ((fac-1.0f)*j)/band_width[k])*ni_gain[k] * rand;
                }

                j++;
            }
        }
        else
        {
            for (i = band_start[k]; i <= band_end[k]; i++)
            {
                if (y2hat[i] == 0)
                {
                    rand = own_random (&ni_seed) / 32768.0f;
                    y2hat[i] += ni_gain[k] * rand;
                }
            }
        }
    }

    mvr2r(env, last_env, ni_end_band);
    mvr2r(ni_gain, last_ni_gain, ni_end_band);
    *last_max_pos_pulse = max_pos_pulse;

    return;
}

