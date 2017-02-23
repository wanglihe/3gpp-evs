/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_enc.h"
#include "prot.h"

/*-------------------------------------------------------------------*
 * Local constants
 *-------------------------------------------------------------------*/

#define MAX_DELTA    16                   /* half-length of the delta search      */
#define COR_BUF_LEN  (L_INTERPOL1*2 + MAX_DELTA*2 + 1)

/*-------------------------------------------------------------------*
 * pitch_ol2()
 *
 * Open-loop pitch precision improvement with 1/4 resolution
 * The pitch is searched in the interval <pitch_ol-delta, pitch_ol+delta),
 * i.e. the value pitch_ol + delta is not a part of the interval
 *-------------------------------------------------------------------*/

void pitch_ol2(
    const short pit_min,        /* i  : pit_min value                                    */
    const short pitch_ol,       /* i  : pitch to be improved                             */
    float *pitch_fr,      /* o  : adjusted 1/4 fractional pitch                    */
    float *voicing_fr,    /* o  : adjusted 1/4 fractional voicing                  */
    const short pos,            /* i  : position in frame where to calculate the improv. */
    const float *wsp,           /* i  : weighted speech for current frame and look-ahead */
    const short delta           /* i  : delta for pitch search                           */
)
{
    short i, t, t0, t1, step, fraction, t0_min, t0_max, t_min, t_max;
    float temp, cor_max, enr_wsp, enr_old, cor[COR_BUF_LEN], *pt_cor, wsp_fr[L_SUBFR];
    const float *pt_wsp;

    t0_min = pitch_ol - delta;
    t0_max = pitch_ol + delta - 1;

    if( t0_min < pit_min )
    {
        t0_min = pit_min;
    }
    t_min = t0_min - L_INTERPOL1;

    if( t0_max > PIT_MAX )
    {
        t0_max = PIT_MAX;
    }
    t_max = t0_max + L_INTERPOL1;

    pt_wsp = wsp + pos;
    pt_cor = cor;
    for ( t=t_min; t<=t_max; t++ )
    {
        *pt_cor++ = dotp( pt_wsp, pt_wsp-t, L_SUBFR );
    }

    pt_cor = cor + L_INTERPOL1;
    cor_max = *pt_cor++;
    t1 = t0_min;
    for ( t=t0_min+1; t<=t0_max; t++ )
    {
        if ( *pt_cor > cor_max  )
        {
            cor_max = *pt_cor;
            t1 = t;
        }
        pt_cor++;
    }

    /*------------------------------------------------------------------*
     * Search fractional pitch with 1/4 subsample resolution.
     * search the fractions around t0 and choose the one which maximizes
     * the interpolated normalized correlation.
     *-----------------------------------------------------------------*/

    pt_cor = cor + L_INTERPOL1 - t0_min;
    t0 = t1;

    step = 1;                /* 1/4 subsample resolution */
    fraction = 1;

    if (t0 == t0_min)        /* Limit case */
    {
        fraction = 0;
        cor_max = interpolation( &pt_cor[t0], E_ROM_inter4_1, fraction, PIT_UP_SAMP, 4 );
    }
    else                     /* Process negative fractions */
    {
        t0--;
        cor_max = interpolation( &pt_cor[t0], E_ROM_inter4_1, fraction, PIT_UP_SAMP, 4 );
        for ( i=(fraction+step); i<=3; i=i+step )
        {
            temp = interpolation( &pt_cor[t0], E_ROM_inter4_1, i, PIT_UP_SAMP, 4 );
            if (temp > cor_max)
            {
                cor_max = temp;
                fraction = i;
            }
        }
    }

    for ( i=0; i<=3; i=i+step )     /* Process positive fractions */
    {
        temp = interpolation( &pt_cor[t1], E_ROM_inter4_1, i, PIT_UP_SAMP, 4 );
        if (temp > cor_max)
        {
            cor_max = temp;
            fraction = i;
            t0 = t1;
        }
    }

    *pitch_fr = t0 + (float)fraction / 4.0f;
    pred_lt4( pt_wsp, wsp_fr, t0, fraction, L_SUBFR, E_ROM_inter4_1, 4, PIT_UP_SAMP);

    enr_wsp = dotp( pt_wsp, pt_wsp, L_SUBFR ) + 0.01f;
    enr_old = dotp( wsp_fr, wsp_fr, L_SUBFR ) + 0.01f;
    *voicing_fr = cor_max * inv_sqrt(enr_wsp * enr_old);

    return;
}


/*-------------------------------------------------------------------*
 * StableHighPitchDetect()
 *
 * Very short stable pitch detection
 *-------------------------------------------------------------------*/

void StableHighPitchDetect(
    short *flag_spitch,       /* o  : flag to indicate very short stable pitch */
    short pitch[],            /* i/o: OL pitch buffer                         */
    const float voicing[],          /* i  : OL pitch gains                          */
    const float Bin_E[],            /* i  : per bin log energy spectrum             */
    const float wsp[],              /* i  : weighted speech                         */
    const short localVAD,
    float *voicing_sm,        /* i/o: smoothed open-loop pitch gains          */
    float *voicing0_sm,       /* i/o: smoothed high pitch gains               */
    float *LF_EnergyRatio_sm, /* i/o: smoothed [0, 300Hz] relative peak energy*/
    short *predecision_flag,  /* i/o: predecision flag                        */
    float *diff_sm,           /* i/o: smoothed pitch frequency difference     */
    float *energy_sm          /* i/o: smoothed energy around pitch frequency  */
)
{
    short i, pitch_freq_point, pit_min_up;
    short T, Tp, pit_min;

    float voicing_m;
    float energy0, energy1, ratio, cor_max, diff, sum_energy;
    const float *pt_wsp;

    voicing_m = mean( voicing, 3 );
    *voicing_sm = 0.75f * (*voicing_sm) + 0.25f * voicing_m;


    /* initial short pitch possibility pre-decision */
    pitch_freq_point = (short)(L_FFT/pitch[1] + 0.5f );
    diff = 0.0f;
    sum_energy = 0.0f;

    for( i=1; i<2*pitch_freq_point; i++ )
    {
        diff += (Bin_E[pitch_freq_point] - Bin_E[i]);
        sum_energy += Bin_E[i];
    }
    sum_energy /= (2*pitch_freq_point-1);

    *diff_sm = 0.2f * diff  + 0.8f * *diff_sm;
    *energy_sm = 0.2f * sum_energy + 0.8f * *energy_sm;
    diff /= sum_energy;

    if( *diff_sm < -10 && *energy_sm < 38.5 && diff < -0.8 )
    {
        *predecision_flag = 1;
    }

    if( *diff_sm > 10 && *energy_sm > 83 && diff > 0.5 )
    {
        *predecision_flag = 0;
    }

    /* short pitch possiblity pre-decision */
    maximum(Bin_E, 7, &energy0);
    maximum(Bin_E+8, 7, &energy1);
    ratio = max(energy1-energy0,0);
    ratio *= max(voicing_m,0);

    *LF_EnergyRatio_sm = (15*(*LF_EnergyRatio_sm) + ratio)/16;

    if( *LF_EnergyRatio_sm > 35 || ratio > 50 )
    {
        *predecision_flag = 1;
    }

    if( *LF_EnergyRatio_sm < 16 )
    {
        *predecision_flag = 0;
    }

    /* short pitch candidate detection */
    Tp = pitch[1];
    cor_max = 0;

    pt_wsp = wsp + 3*L_SUBFR;
    pit_min = PIT_MIN_DOUBLEEXTEND;
    pit_min_up = PIT_MIN;

    for( T=pit_min; T<=pit_min_up; T++ )
    {
        energy1 = dotp( pt_wsp, pt_wsp-T, L_SUBFR );

        if( energy1 > cor_max || T == pit_min )
        {
            cor_max = energy1;
            Tp = T;
        }
    }

    energy0 = dotp( pt_wsp, pt_wsp, L_SUBFR ) + 0.01f;
    energy1 = dotp( pt_wsp-Tp, pt_wsp-Tp, L_SUBFR ) + 0.01f;
    cor_max *= inv_sqrt( energy0*energy1 );
    *voicing0_sm = 0.75f*(*voicing0_sm) + 0.25f*cor_max;

    /* final short pitch correction */
    *flag_spitch = 0;
    if( localVAD && *predecision_flag && *voicing0_sm > 0.65f && *voicing0_sm > 0.7f * (*voicing_sm) )
    {
        *flag_spitch = 1;

        pitch[0] = Tp;
        pitch[1] = Tp;
        pitch[2] = Tp;
    }

    return;
}

/*-------------------------------------------------------------------*
 * pitchDoubling_det()
 * Multiple pitch doubling detector
 *
 *-------------------------------------------------------------------*/

void pitchDoubling_det(
    float *wspeech,
    short *pitch_ol,
    float *pitch_fr,
    float *voicing_fr
)
{
    float new_op_fr[2];
    float new_voicing[2];
    short new_Top[2];
    int m, T;

    /*save initial values*/
    new_Top[0]=pitch_ol[0];
    new_Top[1]=pitch_ol[1];
    for(m=2; m<5; m++)
    {
        T = pitch_ol[0]/m;
        if(T>=PIT_MIN_12k8)
        {
            pitch_ol2( PIT_MIN_SHORTER, T, &new_op_fr[0], &new_voicing[0], 0, wspeech, 2 );
            pitch_ol2( PIT_MIN_SHORTER, T, &new_op_fr[1], &new_voicing[1], L_SUBFR, wspeech, 2 );

            if((new_voicing[0]+new_voicing[1])>(voicing_fr[0]+voicing_fr[1]))
            {
                new_Top[0]=T;
                pitch_fr[0]=new_op_fr[0];
                pitch_fr[1]=new_op_fr[1];
                voicing_fr[0]=new_voicing[0];
                voicing_fr[1]=new_voicing[1];
            }
        }

        T = pitch_ol[1]/m;
        if(T>=PIT_MIN_12k8)
        {
            pitch_ol2( PIT_MIN_SHORTER, T, &new_op_fr[0], &new_voicing[0], 2*L_SUBFR, wspeech, 2 );
            pitch_ol2( PIT_MIN_SHORTER, T, &new_op_fr[1], &new_voicing[1], 3*L_SUBFR, wspeech, 2 );

            if((new_voicing[0]+new_voicing[1])>(voicing_fr[2]+voicing_fr[3]))
            {
                new_Top[1]=T;
                pitch_fr[2]=new_op_fr[0];
                pitch_fr[3]=new_op_fr[1];
                voicing_fr[2]=new_voicing[0];
                voicing_fr[3]=new_voicing[1];
            }
        }
    }
    pitch_ol[0]=new_Top[0];
    pitch_ol[1]=new_Top[1];

    return;
}
