/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "prot.h"


/*-----------------------------------------------------------------*
 * Local constants
 *-----------------------------------------------------------------*/

#define THR_CORR_MAX    60      /* upper threshold of multi-harm. correlation */
#define THR_CORR_MIN    49      /* lower threshold of multi-harm. correlation */
#define THR_CORR_STEP   0.2f    /* step for the threshold of multi-harm. correlation */

/*---------------------------------------------------------------------*
 * multi_harm()
 *
 * Perform multi-harmonic analysis, information used for UV and VAD decision
 *---------------------------------------------------------------------*/

short multi_harm(                  /* o  : frame multi-harmonicity (1-harmonic, 0-not)   */
    const float Bin_E[],           /* i  : log-energy spectrum of the current frame      */
    float old_S[],           /* i/o: prev. log-energy spectrum w. subtracted floor */
    float cor_map_LT[],      /* i/o: LT correlation map                            */
    float *multi_harm_limit, /* i/o: multi harminic threshold                      */
    const long  total_brate,       /* i  : total bitrate             */
    const short bwidth,            /* i  : input signal bandwidth    */
    short *cor_strong_limit, /* i/o: HF correlation indicator  */
    float *st_mean_avr_dyn,  /* i/o: long term average dynamic */
    float *st_last_sw_dyn,   /* i/o: last dynamic              */
    float *cor_map_sum,
    float *sp_floor          /* o: noise floor estimate        */
)
{
    short i, j, k, L, stemp, N_mins, ind_mins[L_FFT/4], *pt_mins, harm;
    float ftemp, ftemp2, flor, step, corx2, cory2, corxy, cor, cor_map_LT_sum, cor_strong, S[L_FFT/2];
    float mean_dyn;

    /*------------------------------------------------------------------*
     * initialization
     *------------------------------------------------------------------*/

    if( bwidth == NB )
    {
        /* length of the useful part of the spectrum (up to 4kHz) */
        L = 76;
    }
    else
    {
        /* length of the useful part of the spectrum (up to 6.4kHz) */
        L = L_FFT/2;
    }

    mvr2r( Bin_E, S, L );

    /*------------------------------------------------------------------*
     * searching of spectral maxima and minima
     *------------------------------------------------------------------*/

    pt_mins = ind_mins;

    /* index of the first minimum */
    if (Bin_E[0] < Bin_E[1])

    {
        *pt_mins++ = 0;
    }

    for (i=1; i<L-1; i++)
    {
        /* minimum found */
        if ( Bin_E[i] < Bin_E[i-1] && Bin_E[i] < Bin_E[i+1] )
        {
            *pt_mins++ = i;
        }
    }

    /* index of the last minimum */
    if (Bin_E[L-1] < Bin_E[L-2])
    {
        *pt_mins++ = L-1;
    }

    /* total number of minimas found */
    N_mins = (short)(pt_mins - ind_mins - 1);

    /*------------------------------------------------------------------*
     * calculation of the spectral floor
     * subtraction of the spectral floor
     *------------------------------------------------------------------*/

    set_f(S, 0, L);
    if (N_mins > 0)
    {
        *sp_floor = 0;
        for (i=0; i<N_mins; ++i)
        {
            *sp_floor += Bin_E[ind_mins[i]];
        }
        *sp_floor /= (float)N_mins;
        set_f( S, 0, ind_mins[0]);
        set_f( &S[ind_mins[N_mins]], 0, L - ind_mins[N_mins]);

        pt_mins = ind_mins;
        flor = 0;
        step = 0;

        for (i=ind_mins[0]; i<ind_mins[N_mins]; i++)
        {
            /* we are at the end of the next minimum */
            if (i == *pt_mins)
            {
                pt_mins++;
                flor = Bin_E[i];

                /* calculate the new step */
                step = (Bin_E[*pt_mins] - Bin_E[i]) / (*pt_mins-i);
            }

            /* subtract the floor */
            if (Bin_E[i] > flor)
            {
                S[i] = Bin_E[i] - flor;
            }
            else
            {
                S[i] = 0;
            }

            /* update the floor */
            flor += step;
        }
    }
    else
    {
        *sp_floor = Bin_E[0];
    }
    *sp_floor *= 1.0f/log(10);

    /* calculate the maximum dynamic per band */
    mean_dyn = mean( &S[L-40], 40 );
    mean_dyn = 0.6f * *st_mean_avr_dyn + 0.4f * mean_dyn;

    if ( mean_dyn < 9.6f && *cor_strong_limit != 0 )
    {
        *cor_strong_limit = 0;
        *st_last_sw_dyn = mean_dyn;
    }
    else if ((mean_dyn - *st_last_sw_dyn) > 4.5f)
    {
        *cor_strong_limit = 1;
    }

    if( total_brate < ACELP_9k60 || total_brate > ACELP_16k40 )
    {
        *cor_strong_limit = 1;
    }

    *st_mean_avr_dyn = mean_dyn;

    /*------------------------------------------------------------------*
     * calculation of the correlation map
     *------------------------------------------------------------------*/

    if (N_mins > 0)
    {
        corx2 = 0;
        corxy = 0;
        stemp = ind_mins[0];
        ftemp = old_S[stemp];
        cory2 = ftemp * ftemp;
        k = 1;
        for (i = stemp+1; i <= ind_mins[N_mins]; i++)
        {
            if (i == ind_mins[k])
            {
                /* include the last peak point (new minimum) to the corr. sum */
                ftemp = old_S[i];
                cory2 += ftemp * ftemp;

                /* calculation of the norm. peak correlation */
                if ((corx2 == 0) || (cory2 == 0))
                {
                    cor = 0;
                }
                else
                {
                    cor = corxy * corxy / (corx2 * cory2);
                }

                /* save the norm. peak correlation in the correlation map */
                for (j=ind_mins[k-1]; j<ind_mins[k]; j++)
                {
                    old_S[j] = S[j];
                    S[j] = cor;
                }

                corx2 = 0;
                cory2 = 0;
                corxy = 0;

                k++;
            }

            ftemp = S[i];
            ftemp2 = old_S[i];
            corx2 += ftemp * ftemp;
            cory2 += ftemp2 * ftemp2;
            corxy += ftemp * ftemp2;
        }

        mvr2r( S, old_S, ind_mins[0]);
        mvr2r( &S[ind_mins[N_mins]], &old_S[ind_mins[N_mins]], L - ind_mins[N_mins]);
    }

    /*------------------------------------------------------------------*
     * updating of the long-term correlation map
     * summation of the long-term correlation map
     *------------------------------------------------------------------*/

    cor_strong = 0;
    *cor_map_sum = 0;

    for (i=0; i<L; i++)
    {
        *cor_map_sum += S[i];
        cor_map_LT[i] = M_ALPHA * cor_map_LT[i] + (1-M_ALPHA) * S[i];
        if (cor_map_LT[i] > 0.95f)
        {
            cor_strong = 1;
        }
    }


    /* summation of the LT correlation map */
    cor_map_LT_sum = sum_f(cor_map_LT, L);

    if ( bwidth == NB )
    {
        cor_map_LT_sum *= 1.53f;
        *cor_map_sum *= 1.53f;
    }


    /* final decision about multi-harmonicity */
    if ( (cor_map_LT_sum > *multi_harm_limit) || (cor_strong == 1) )
    {
        harm = 1;
    }
    else
    {
        harm = 0;
    }

    /*------------------------------------------------------------------*
     * updating of the decision threshold
     *------------------------------------------------------------------*/

    if (cor_map_LT_sum > THR_CORR)
    {
        *multi_harm_limit -= THR_CORR_STEP;
    }
    else
    {
        *multi_harm_limit += THR_CORR_STEP;
    }

    if (*multi_harm_limit > THR_CORR_MAX)
    {
        *multi_harm_limit = THR_CORR_MAX;
    }

    if (*multi_harm_limit < THR_CORR_MIN)
    {
        *multi_harm_limit = THR_CORR_MIN;
    }

    if (N_mins <= 0)
    {
        set_f(old_S, 0, L);
    }

    return harm;
}
