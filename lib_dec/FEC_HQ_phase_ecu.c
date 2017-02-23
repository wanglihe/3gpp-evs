/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "rom_dec.h"
#include "rom_com.h"
#include "cnst.h"
#include "prot.h"
#include "assert.h"


/*---------------------------------------------------------------------*
 * Local constants
 *---------------------------------------------------------------------*/

#define FEC_MAX                               512
#define FEC_NB_PULSE_MAX                      20
#define FEC_FFT_MAX_SIZE                      512
#define FEC_DCIM_FILT_SIZE_MAX                60

#define PHASE_DITH                            (PI2)

#define DELTA_CORR                            6                       /* Range for phase correction around peak */
#define THRESH_TR_dB                          10.0f
#define THRESH_TR_LIN                         (float)pow(10.0f,THRESH_TR_dB/10.0f)
#define THRESH_TR_LIN_INV                     (float)pow(10.0f,-THRESH_TR_dB/10.0f)
#define MAX_INCREASE_GRPOW                    0.0f                    /* maximum amplification in case of transients */
#define MAX_INCREASE_GRPOW_LIN                (float)pow(10.0f,MAX_INCREASE_GRPOW/10.0f)

#define PHASE_DITH_SCALE                      (float)pow(2.0,-16.0)   /* for scaling random short values  to  +/- pi */

#define BURST_PHDITH_THRESH                   (4-1)                   /* speech start phase dither with <burst_phdith_thresh> losses in a row */
#define BURST_PHDITH_RAMPUP_LEN               2                       /* speech ramp up degree of phase dither over a length of <burst_phdith_rampup_len> frames */
#define BURST_ATT_THRESH                      (3-1)                   /* speech start attenuate with <burst_att_thresh> losses in a row */
#define ATT_PER_FRAME                         4                       /* speech attenuation in dB */
#define BETA_MUTE_THR                         10                      /* time threshold to start beta-noise attenuation */
#define BETA_MUTE_FAC                         0.5f                    /* attenuation factor per additional bad frame */

#define LGW32k                                7
#define LGW16k                                6
#define LGW48k                                LGW32k+1                /* Use the same frequency groups as for SWB + 1 */

#define L_TRANA_LOG32k                        8
#define L_TRANA_LOG16k                        7

#define PFIND_SENS                            0.97f                   /* peakfinder sensitivity */

/*---------------------------------------------------------------------*
 * Local function
 *---------------------------------------------------------------------*/

short rand_phase( const short seed, float *sin_F, float *cos_F );


/*-------------------------------------------------------------------*
 * mult_rev2()
 *
 * Multiplication of two vectors second vector is multiplied in reverse order
 *-------------------------------------------------------------------*/

static void mult_rev2(
    const float x1[],   /* i  : Input vector 1                                   */
    const float x2[],   /* i  : Input vector 2                                   */
    float y[],    /* o  : Output vector that contains vector 1 .* vector 2 */
    const short N       /* i  : Vector length                                    */
)
{
    short i,j ;

    for (i=0,j=N-1 ; i<N ; i++,j--)
    {
        y[i] = x1[i] * x2[j] ;
    }

    return;
}


/*-------------------------------------------------------------------*
 * fft_spec2()
 *
 * Square magnitude of fft spectrum
 *-------------------------------------------------------------------*/

static void fft_spec2(
    float x[],          /* i/o : Input vector: complex spectrum -> square magnitude spectrum  */
    const short N       /* i   : Vector lenght                                                */
)
{
    short i, j;

    for (i = 1, j = N-1; i < N/2; i++, j--)
    {
        x[i] = x[i]*x[i] + x[j]*x[j];
    }

    x[0] *= x[0];
    x[N/2] *= x[N/2];

    return;
}

/*------------------------------------------------------------------*
 * rand_phase()
 *
 * randomized phase in form of sin and cos components
 *------------------------------------------------------------------*/

short rand_phase(
    const short seed,
    float *sin_F,
    float *cos_F
)
{
    const float *sincos = sincos_t_ext + 128;
    short seed2 = seed;
    own_random(&seed2);

    if (seed2 & 0x40)
    {
        *sin_F = sincos[seed2 >> 8];
    }
    else
    {
        *sin_F = -sincos[seed2 >> 8];
    }

    if (seed2 & 0x80)
    {
        *cos_F = sincos[-(seed2 >> 8)];
    }
    else
    {
        *cos_F = -sincos[-(seed2 >> 8)];
    }

    return seed2;
}

/*------------------------------------------------------------------*
 * trans_ana()
 *
 * Transient analysis
 *------------------------------------------------------------------*/

static void trans_ana(
    const float *xfp,                /* i  : Input signal                                         */
    float *mag_chg,            /* i/o: Magnitude modification                               */
    float *ph_dith,            /* i/o: Phase dither                                         */
    float *mag_chg_1st,        /* i/o: per band magnitude modifier for transients           */
    const short output_frame,        /* i  : Frame length                                         */
    const short time_offs,           /* i  : Time offset                                          */
    const float est_mus_content,     /* i  : 0.0=speech_like ... 1.0=Music    (==st->env_stab )   */
    const short last_fec,            /* i  : signal that previous frame was concealed with fec_alg*/
    float *alpha,                    /* o  : Magnitude modification factors for fade to average   */
    float *beta,                     /* o  : Magnitude modification factors for fade to average   */
    float *beta_mute,                /* o  : Factor for long-term mute                            */
    float Xavg[LGW_MAX]              /* o  : Frequency group average gain to fade to              */
)
{
    const float *w_hamm;
    float grp_pow_chg, att_val, att_degree;
    float xfp_left[L_TRANA48k], xfp_right[L_TRANA48k];
    float gr_pow_left[LGW_MAX], gr_pow_right[LGW_MAX];
    const float *xfp_;
    short Ltrana, Ltrana_2, Lprot, LtranaLog = 0, Lgw, k,  burst_len;
    short att_always[LGW_MAX];             /* fixed attenuation per frequency group if set to 1*/
    short burst_phdith_thresh = BURST_PHDITH_THRESH;     /* speech settings */
    short burst_att_thresh = BURST_ATT_THRESH ;
    float att_per_frame = ATT_PER_FRAME;
    short burst_phdith_rampup_len = BURST_PHDITH_RAMPUP_LEN;
    short tr_dec[LGW_MAX];

    /* check burst error */
    burst_len = time_offs/output_frame + 1;

    set_s(att_always,0,LGW_MAX);
    *ph_dith = 0.0f;

    /* softly shift attenuation just a bit later for estimated "stable" music_content */
    burst_phdith_thresh = BURST_PHDITH_THRESH+(short)(est_mus_content*1.0f + 0.5f);
    burst_att_thresh = BURST_ATT_THRESH + (short)(est_mus_content*1.0f + 0.5f);
    att_per_frame  = ATT_PER_FRAME - (short)(est_mus_content*1.0f + 0.5f);     /* only slighty less att for music */
    att_per_frame *=0.1f;

    if ( burst_len > burst_phdith_thresh )
    {
        /* increase degree of dither */
        *ph_dith = PHASE_DITH * min(1.0f,((float)burst_len - (float)burst_phdith_thresh)/(float)burst_phdith_rampup_len);
    }

    att_degree = 0;
    if (burst_len > burst_att_thresh)
    {
        set_s(att_always,1,LGW_MAX);

        /* increase degree of attenuation */
        if( burst_len - burst_att_thresh <= PH_ECU_MUTE_START )
        {
            att_degree = (float)(burst_len - burst_att_thresh) * att_per_frame;
        }
        else
        {
            att_degree = (float) PH_ECU_MUTE_START * att_per_frame + (burst_len - burst_att_thresh - PH_ECU_MUTE_START) * 6.0206f;
        }
    }

    Lprot = ( 2* output_frame * 4)/ 5;   /* 4/5==1024/1280, keep mult within short */
    Ltrana = Lprot/QUOT_LPR_LTR;
    Ltrana_2 = Ltrana/2;

    if (output_frame == L_FRAME48k)
    {
        w_hamm = w_hamm48k_2;
        Lgw = LGW48k;
    }
    else if (output_frame == L_FRAME32k)
    {
        w_hamm = w_hamm32k_2;
        LtranaLog = L_TRANA_LOG32k;
        Lgw = LGW32k;
    }
    else
    {
        w_hamm = w_hamm16k_2;
        LtranaLog = L_TRANA_LOG16k;
        Lgw = LGW16k;
    }

    if (burst_len <= 1 || (burst_len == 2 && last_fec) )
    {
        set_f(alpha,1.0f,LGW_MAX);
        set_f(beta,0.0f,LGW_MAX);
        *beta_mute = BETA_MUTE_FAC_INI;

        /* apply hamming window */
        v_mult( xfp, w_hamm, xfp_left, Ltrana_2 );
        mult_rev2( xfp+Ltrana_2, w_hamm, xfp_left+Ltrana_2, Ltrana_2 );

        xfp_ = xfp + Lprot - Ltrana;
        v_mult( xfp_, w_hamm, xfp_right, Ltrana_2 );
        mult_rev2( xfp_+Ltrana_2, w_hamm, xfp_right+Ltrana_2, Ltrana_2 );

        /* spectrum */
        if (output_frame == L_FRAME48k)
        {
            fft3(xfp_left,xfp_left,Ltrana);
            fft3(xfp_right,xfp_right,Ltrana);
        }
        else
        {
            fft_rel(xfp_left,Ltrana,LtranaLog);
            fft_rel(xfp_right,Ltrana,LtranaLog);
        }

        /* square representation */
        fft_spec2(xfp_left,Ltrana);
        fft_spec2(xfp_right,Ltrana);

        /* band powers in frequency groups
        exclude bin at 0 and at EVS_PI from calculation */
        xfp_left[Ltrana_2] = 0.0f;
        xfp_right[Ltrana_2] = 0.0f;
    }

    for ( k = 0; k < Lgw; k++ )
    {
        if (burst_len <= 1 || (burst_len == 2 && last_fec) )
        {
            gr_pow_left[k]  = sum_f(xfp_left +gw[k],gw[k+1]-gw[k]);
            gr_pow_right[k] = sum_f(xfp_right+gw[k],gw[k+1]-gw[k]);

            /* check if transient in any of the bands */
            gr_pow_left[k] += FLT_MIN;       /* otherwise div by zero may occur */
            gr_pow_right[k] += FLT_MIN;

            Xavg[k] = sqrt(0.5f*(gr_pow_left[k]+gr_pow_right[k])/(float)(gw[k+1]-gw[k]));

            grp_pow_chg = gr_pow_right[k] / gr_pow_left[k];

            /* dither phase in case of transient */
            /* separate transition detection and application of forced burst dithering */
            tr_dec[k] = (grp_pow_chg > THRESH_TR_LIN) || (grp_pow_chg < THRESH_TR_LIN_INV);

            /* magnitude modification */
            if ( tr_dec[k] || att_always[k])
            {
                att_val = min(MAX_INCREASE_GRPOW_LIN,grp_pow_chg);
                att_val = (float)sqrt(att_val);
                mag_chg_1st[k] = att_val;
                mag_chg[k] = att_val;
            }
            else
            {
                mag_chg_1st[k] = 1.0f;
                mag_chg[k] = 1.0f;
            }
        }
        else
        {
            if( burst_len < OFF_FRAMES_LIMIT )
            {
                mag_chg[k] = mag_chg_1st[k] * (float)pow(10.0,-att_degree/20.0);
            }
            else
            {
                mag_chg[k] = 0;
            }
            if (burst_len > BETA_MUTE_THR)
            {
                *beta_mute *= BETA_MUTE_FAC;
            }
            alpha[k] = mag_chg[k];
            beta[k] = sqrt(1.0f - SQR(alpha[k])) **beta_mute;
            if (k>=LGW32k-1)
            {
                beta[k] *=0.1f;
            }
            else if (k>=LGW16k-1)
            {
                beta[k] *=0.5f;
            }
        }
    }

    return;
}

/*------------------------------------------------------------------*
 * peakfinder()
 *
 * Peak-picking algorithm
 *------------------------------------------------------------------*/

static void peakfinder(
    const float *x0,    /* i : vector from which the maxima will be found                     */
    const short len0,   /* i : length of input vector                                         */
    short *plocs,       /* o : the indicies of the identified peaks in x0                     */
    short *cInd,        /* o : number of identified peaks                                     */
    const float sel     /* i : The amount above surrounding data for a peak to be identified  */
)
{
    float minMag, tempMag, leftMin;
    float dx0[L_PROT48k_2], x[L_PROT48k_2+1], peakMag[MAX_PLOCS];
    short k, i, len, tempLoc = 0, foundPeak, ii, xInd;
    short *ind, indarr[L_PROT48k_2+1], peakLoc[MAX_PLOCS];

    ind = indarr;

    /* Find derivative */
    v_sub(x0+1,x0,dx0,len0-1);

    /* This is so we find the first of repeated values */
    for (i=0; i<len0-1; i++)
    {
        if (dx0[i] == 0.0f)
        {
            dx0[i] = -1.0e-12f;
        }
    }

    /* Find where the derivative changes sign
       Include endpoints in potential peaks and valleys */
    k=0;
    x[k] = x0[0];
    ind[k++] = 0;
    for ( i = 1; i < len0-1; i++)
    {
        if (dx0[i-1] * dx0[i] < 0)
        {
            ind[k] = i;
            x[k++] = x0[i];
        }
    }
    ind[k] = len0-1;
    x[k++] = x0[len0-1];

    /* x only has the peaks, valleys, and endpoints */
    len = k;
    minimum(x,len,&minMag);

    if (len > 2)
    {
        /* Set initial parameters for loop */
        tempMag = minMag;
        foundPeak = 0;
        leftMin = minMag;

        /* Deal with first point a little differently since tacked it on
           Calculate the sign of the derivative since we taked the first point
           on it does not necessarily alternate like the rest. */

        /* The first point is larger or equal to the second */
        if (x[0] >= x[1])
        {
            ii = -1;
            if (x[1] >= x[2]) /* x[1] is not extremum -> overwrite with x[0] */
            {
                x[1] = x[0];
                ind[1] = ind[0];
                ind++;
                len--;
            }
        }
        else /* First point is smaller than the second */
        {
            ii = 0;
            if (x[1] < x[2]) /* x[1] is not extremum -> overwrite with x[0] */
            {
                x[1] = x[0];
                ind[1] = ind[0];
                ind++;
                len--;
            }
        }

        *cInd = 0;

        /* Loop through extrema which should be peaks and then valleys */
        while (ii < len-1)
        {
            ii++;   /* This is a peak */

            /*Reset peak finding if we had a peak and the next peak is bigger
              than the last or the left min was small enough to reset.*/
            if (foundPeak)
            {
                tempMag = minMag;
                foundPeak = 0;
            }

            /* Make sure we don't iterate past the length of our vector */
            if (ii == len-1)
            {
                break;  /* We assign the last point differently out of the loop */
            }

            /* Found new peak that was larger than temp mag and selectivity larger
               than the minimum to its left. */
            if ((x[ii] > tempMag) && (x[ii] > leftMin + sel))
            {
                tempLoc = ii;
                tempMag = x[ii];
            }

            ii++; /* Move onto the valley */

            /* Come down at least sel from peak */
            if (!foundPeak && (tempMag > sel + x[ii]))
            {
                foundPeak = 1;                  /* We have found a peak */
                leftMin = x[ii];
                peakLoc[*cInd] = tempLoc;       /* Add peak to index */
                peakMag[*cInd] = tempMag;
                (*cInd)++;
            }
            else if (x[ii] < leftMin) /* New left minimum */
            {
                leftMin = x[ii];
            }
        }

        /* Check end point */
        if (x[len-1] > tempMag && x[len-1] > leftMin + sel)
        {
            peakLoc[*cInd] = len-1;
            peakMag[*cInd] = x[len-1];
            (*cInd)++;
        }
        else if (!foundPeak && tempMag > minMag) /* Check if we still need to add the last point */
        {
            peakLoc[*cInd] = tempLoc;
            peakMag[*cInd] = tempMag;
            (*cInd)++;
        }

        /* Create output */
        for (i=0; i<*cInd; i++)
        {
            plocs[i] = ind[peakLoc[i]];
        }
    }
    else /* This is a monotone function where an endpoint is the only peak */
    {
        xInd = (x[0] > x[1]) ? 0 : 1;
        peakMag[0] = x[xInd];
        if (peakMag[0] > minMag + sel)
        {
            plocs[0] = ind[xInd];
            *cInd=1;
        }
        else
        {
            *cInd=0;
        }
    }

    return;
}

/*-------------------------------------------------------------------*
 * imax()
 *
 * Get interpolated maximum position
 *-------------------------------------------------------------------*/

static float imax(const float *y)
{
    float posi, y1, y2,y3,y3_y1,y2i;
    /* Seek the extrema of the parabola P(x) defined by 3 consecutive points so that P([-1 0 1]) = [y1 y2 y3] */
    float ftmp_den1,ftmp_den2;
    y1 = y[0];
    y2 = y[1];
    y3 = y[2];
    y3_y1 = y3-y1;
    ftmp_den1 =   (y1+y3-2*y2);
    ftmp_den2 =   (4*y2 - 2*y1 - 2*y3);

    if( ftmp_den2==0.0f  || ftmp_den1==0.0f )
    {
        return (0.0f);  /* early exit with left-most value */
    }

    y2i = -0.125f * SQR(y3_y1) /(ftmp_den1) + y2;
    /* their corresponding normalized locations */
    posi = y3_y1/(ftmp_den2);

    /* Interpolated maxima if locations are not within [-1,1], calculated extrema are ignored */
    if (posi >= 1.0f  || posi <= -1.0f)
    {
        posi = y3 > y1 ? 1.0f : -1.0f;
    }
    else
    {
        if (y1 >= y2i)
        {
            posi = (y1 > y3) ? -1.0f : 1.0f;
        }
        else if (y3 >= y2i)
        {
            posi = 1.0f;
        }
    }

    return posi + 1.0f;
}


/*-------------------------------------------------------------------*
 * spec_ana()
 *
 * Spectral analysis
 *-------------------------------------------------------------------*/

static void spec_ana(
    const float *prevsynth,             /* i  : Input signal                                    */
    short *plocs,                       /* o  : The indicies of the identified peaks            */
    float *plocsi,                      /* o  : Interpolated positions of the identified peaks  */
    short *num_plocs,                   /* o  : Number of identified peaks                      */
    float *X_sav,                       /* o  : Stored fft spectrum                             */
    const short output_frame,           /* i  : Frame length                                    */
    const short bwidth                  /* i  : Encoded bandwidth index                          */
)
{
    short i, Lprot, LprotLog2=0, hamm_len2=0, Lprot2_1, m;
    const float *w_hamm = NULL;
    float xfp[L_PROT48k];
    float Xmax, Xmin, sel;
    short stop_band_start;
    short stop_band_length;

    Lprot = 2*output_frame * L_PROT32k / 1280;
    Lprot2_1 = Lprot/2+1;

    if (output_frame == L_FRAME48k)
    {
        w_hamm = w_hamm_sana48k_2;
        hamm_len2 = L_PROT_HAMM_LEN2_48k;
    }
    else if (output_frame == L_FRAME32k)
    {
        w_hamm = w_hamm_sana32k_2;
        hamm_len2 = L_PROT_HAMM_LEN2_32k;
        LprotLog2 = 10;
    }
    else
    {
        w_hamm = w_hamm_sana16k_2;
        hamm_len2 = L_PROT_HAMM_LEN2_16k;
        LprotLog2 = 9;
    }

    /* Apply hamming-rect window */
    v_mult(prevsynth,w_hamm,xfp,hamm_len2);
    mvr2r(prevsynth+hamm_len2,xfp+hamm_len2,Lprot - 2*hamm_len2);
    mult_rev2(prevsynth + Lprot - hamm_len2,w_hamm,xfp  + Lprot - hamm_len2,hamm_len2);

    /* Spectrum */
    if (output_frame == L_FRAME48k)
    {
        fft3(xfp,xfp,Lprot);
    }
    else
    {
        fft_rel(xfp,Lprot,LprotLog2);
    }

    /* Apply zeroing of non-coded FFT spectrum */
    if (output_frame > inner_frame_tbl[bwidth])
    {
        stop_band_start = 128 << bwidth;
        stop_band_length = Lprot - (stop_band_start << 1);
        stop_band_start = stop_band_start + 1;
        set_f( xfp + stop_band_start, 0, stop_band_length );
    }

    mvr2r(xfp,X_sav,Lprot);

    /* Magnitude representation */
    fft_spec2(xfp, Lprot);

    for (i=0; i<Lprot2_1 ; i++)
    {
        xfp[i] = (float) sqrt((double)xfp[i]);
    }

    /* Find maxima */
    maximum(xfp,Lprot2_1,&Xmax);
    minimum(xfp,Lprot2_1,&Xmin);
    sel = (Xmax-Xmin)*(1.0f-PFIND_SENS);
    peakfinder(xfp,Lprot2_1, plocs, num_plocs,sel); /* NB peak at xfp[0] and xfp Lprot2_1-1 may occur */

    /* Refine peaks */
    for (m = 0; m < *num_plocs; m++)
    {
        if (plocs[m] == 0)
        {
            plocsi[m] = plocs[m] + imax(&xfp[plocs[m]]);
        }
        else if (plocs[m] == Lprot/2)
        {
            plocsi[m] = plocs[m] - 2 + imax(&xfp[plocs[m]-2]);
        }
        else
        {
            plocsi[m] = plocs[m] - 1 + imax(&xfp[plocs[m]-1]);
        }
    }

    return;
}

/*-------------------------------------------------------------------*
 * subst_spec()
 *
 * Substitution spectrum calculation
 *-------------------------------------------------------------------*/

static void subst_spec(
    const short *plocs,                 /* i   : The indicies of the identified peaks               */
    const float *plocsi,                /* i   : Interpolated positions of the identified peaks     */
    short *num_plocs,             /* i/o : Number of identified peaks                         */
    const short time_offs,              /* i   : Time offset                                        */
    float *X,                     /* i/o : FFT spectrum                                       */
    const float *mag_chg,               /* i   : Magnitude modification                             */
    const float ph_dith,                /* i   : Phase dither                                       */
    const short *is_trans,              /* i   : Transient flags                                    */
    const short output_frame,           /* i   : Frame length                                       */
    short *seed,                  /* i/o : Random seed                                        */
    const float *alpha,                 /* i   : Magnitude modification factors for fade to average */
    const float *beta,                  /* i   : Magnitude modification factors for fade to average */
    float beta_mute,                    /* i   : Factor for long-term mute                          */
    const float Xavg[LGW_MAX]           /* i   : Frequency group averages to fade to                */

)
{
    const float *sincos;
    short Xph_short;
    float corr_phase[MAX_PLOCS], Xph;
    float Lprot_1, cos_F, sin_F, tmp;
    short Lprot, Lecu, m, i, e, im_ind, delta_corr_up, delta_corr_dn, delta_tmp;
    float mag_chg_local;   /* for peak attenuation in burst */
    short k;

    sincos = sincos_t_ext + 128;
    Lprot = (short)(L_PROT32k*output_frame/640);
    Lprot_1 = 1.0f/Lprot;
    Lecu = output_frame*2;

    /* Correction phase of the identified peaks */
    if (is_trans[0] || is_trans[1])
    {
        *num_plocs = 0;
    }
    else
    {
        tmp = PI2 * (Lecu-(Lecu-Lprot)/2 + NS2SA(output_frame*50,PH_ECU_ALDO_OLP2_NS-PH_ECU_LOOKAHEAD_NS) - output_frame/2 + time_offs) * Lprot_1;

        for (m = 0; m < *num_plocs; m++)
        {
            corr_phase[m] = plocsi[m]*tmp;
        }
    }

    i = 1;
    k = 0;
    im_ind = Lprot-1;
    for (m = 0; m < *num_plocs; m++)
    {
        delta_corr_dn = DELTA_CORR;
        delta_corr_up = DELTA_CORR;

        if (m > 0)
        {
            delta_tmp = (plocs[m] - plocs[m-1] - 1)/2;
            if (delta_tmp < DELTA_CORR)
            {
                delta_corr_dn = delta_tmp;
            }
        }

        if (m < *num_plocs-1 )
        {
            delta_tmp = (plocs[m+1] - plocs[m] - 1)/2;
            if ( delta_tmp < DELTA_CORR )
            {
                delta_corr_up = delta_tmp;
            }
        }

        /* Input Xph */
        while ( i < plocs[m]-delta_corr_dn )
        {
            *seed = own_random(seed);

            if (*seed & 0x40)
            {
                sin_F = sincos[*seed >> 8];
            }
            else
            {
                sin_F = -sincos[*seed >> 8];
            }

            if (*seed & 0x80)
            {
                cos_F = sincos[-(*seed >> 8)];
            }
            else
            {
                cos_F = -sincos[-(*seed >> 8)];
            }

            tmp =       (X[i] * cos_F - X[im_ind] * sin_F);
            X[im_ind] = (X[i] * sin_F + X[im_ind] * cos_F);
            if (alpha[k] < 1.0f)
            {
                *seed = rand_phase(*seed,&sin_F,&cos_F);
                X[i] = alpha[k] * tmp            + beta[k] * Xavg[k] * cos_F;
                X[im_ind] = alpha[k] * X[im_ind] + beta[k] * Xavg[k] * sin_F;
            }
            else
            {
                X[i] = mag_chg[k] * tmp;
                X[im_ind] *= mag_chg[k];
            }
            i++;
            im_ind--;
            if (i >= gwlpr[k+1])
            {
                k++;
            }
        }

        e = plocs[m] + delta_corr_up;
        if (e > Lprot/2 - 1)
        {
            e = Lprot/2 - 1;
        }

        Xph = corr_phase[m];
        Xph_short = (short)(((long)(Xph * 512/EVS_PI)) % 32768) & 0x03ff;
        if (Xph_short >= 512)
        {
            sin_F = -sincos_t_ext[Xph_short-512];
            if (Xph_short < 768)
            {
                cos_F = -sincos_t_ext[Xph_short-512+256];
            }
            else
            {
                cos_F = sincos_t_ext[-Xph_short+1024+256];
            }
        }
        else
        {
            sin_F = sincos_t_ext[Xph_short];
            if (Xph_short < 256)
            {
                cos_F = sincos_t_ext[Xph_short+256];
            }
            else
            {
                cos_F = -sincos_t_ext[-Xph_short+256+512];
            }
        }

        while ( i <= e )
        {
            mag_chg_local = mag_chg[k];

            if( ph_dith != 0.0f )
            {
                /* Call phase randomization only when needed */
                Xph = corr_phase[m];
                *seed = own_random(seed);
                Xph +=  *seed  * ph_dith * PHASE_DITH_SCALE;  /* where ph_dith is  0..2PI,  or -2PI (in transient), bin phase scaling factor from trans_ana */

                if( ph_dith > 0.0f )
                {
                    /* up to 6 dB additional att of peaks in non_transient longer bursts, (when peak phase is randomized) */
                    /* 0.5~= sqrt((float)pow(10.0,-6/10.0));    ph_dith= 0..2pi,--> scale=1.0 ...5 */
                    mag_chg_local *= 0.5f + (1.0f-(1.0f/PHASE_DITH)*ph_dith)*0.5f;
                }

                Xph_short = (short)((int)(Xph * 512/EVS_PI) & 0x03ff);

                if (Xph_short >= 512)
                {
                    sin_F = -sincos_t_ext[Xph_short-512];
                    if (Xph_short < 768)
                    {
                        cos_F = -sincos_t_ext[Xph_short-512+256];
                    }
                    else
                    {
                        cos_F = sincos_t_ext[-Xph_short+1024+256];
                    }
                }
                else
                {
                    sin_F = sincos_t_ext[Xph_short];
                    if (Xph_short < 256)
                    {
                        cos_F = sincos_t_ext[Xph_short+256];
                    }
                    else
                    {
                        cos_F = -sincos_t_ext[-Xph_short+256+512];
                    }
                }
            }

            tmp =       (X[i] * cos_F - X[im_ind] * sin_F);
            X[im_ind] = (X[i] * sin_F + X[im_ind] * cos_F);
            if (alpha[k] < 1.0f)
            {
                float alpha_local = mag_chg_local;
                float beta_local = beta_mute * sqrt(1.0f - SQR(alpha_local));
                if (k>=LGW32k-1)
                {
                    beta_local *=0.1f;
                }
                else if (k>=LGW16k-1)
                {
                    beta_local *=0.5f;
                }

                *seed = rand_phase(*seed,&sin_F,&cos_F);
                X[i] = alpha_local * tmp            + beta_local * Xavg[k] * cos_F;
                X[im_ind] = alpha_local * X[im_ind] + beta_local * Xavg[k] * sin_F;
            }
            else
            {
                X[i] = mag_chg_local * tmp;
                X[im_ind] *= mag_chg_local;
            }

            i++;
            im_ind--;
            if (i >= gwlpr[k+1])
            {
                k++;
            }
        }
    }

    while ( i < Lprot/2 )
    {
        *seed = own_random(seed);

        if (*seed & 0x40)
        {
            sin_F = sincos[*seed >> 8];
        }
        else
        {
            sin_F = -sincos[*seed >> 8];
        }

        if (*seed & 0x80)
        {
            cos_F = sincos[-(*seed >> 8)];
        }
        else
        {
            cos_F = -sincos[-(*seed >> 8)];
        }


        tmp =       (X[i] * cos_F - X[im_ind] * sin_F);
        X[im_ind] = (X[i] * sin_F + X[im_ind] * cos_F);
        if (alpha[k] < 1.0f)
        {
            *seed = rand_phase(*seed,&sin_F,&cos_F);
            X[i] = alpha[k] * tmp              + beta[k] * Xavg[k] * cos_F;
            X[im_ind] = alpha[k] * X[im_ind] + beta[k] * Xavg[k] * sin_F;
            im_ind--;
        }
        else
        {
            X[i] = mag_chg[k] * tmp;
            X[im_ind--] *= mag_chg[k];
        }
        i++;

        if (i >= gwlpr[k+1])
        {
            k++;
        }
    }

    return;
}

/*--------------------------------------------------------------------------
 *  rec_wtda()
 *
 *  Windowing and TDA of reconstructed frame
 *--------------------------------------------------------------------------*/

static void rec_wtda(
    float *X,                           /* i  : FFT spectrum                          */
    float *ecu_rec,                     /* o  : Reconstructed frame in tda domain     */
    const short output_frame,           /* i  : Frame length                          */
    const short Lprot,                  /* i  : Prototype frame length                */
    const float fs
)
{
    float old_wtda[L_FRAME48k],xinit[L_FRAME48k],xsubst_[2*L_FRAME48k];
    short timesh;
    short xf_len;
    short i;
    float *p_ecu;
    float g;
    float tbl_delta;

    /* extract reconstructed frame with aldo window */
    set_f(xsubst_,0.0f,output_frame-Lprot/2);
    mvr2r(X,xsubst_+output_frame-Lprot/2,Lprot);
    set_f(xsubst_+output_frame+Lprot/2,0.0f,output_frame-Lprot/2);

    /* Smoothen onset of ECU frame */
    xf_len = (short)((float)output_frame*N_ZERO_MDCT_NS/FRAME_SIZE_NS) - (output_frame - Lprot/2);
    p_ecu = xsubst_ + (output_frame - Lprot/2);
    tbl_delta = 64.0/xf_len; /* 64 samples = 1/4 cycle in sincos_t */
    for ( i = 0; i < xf_len; i++, p_ecu++ )
    {
        g = sincos_t[((short)(i * tbl_delta))];
        g *= g;
        *p_ecu = g * (*p_ecu);
    }

    timesh = NS2SA(fs, 10000000L - PH_ECU_ALDO_OLP2_NS);
    set_f( xinit, 0, timesh );
    mvr2r( xsubst_, xinit+timesh, output_frame-timesh );

    set_f( old_wtda,0,output_frame);

    /* warm up old_wtda */
    wtda( xinit, ecu_rec, old_wtda, ALDO_WINDOW, ALDO_WINDOW, output_frame );
    wtda( xsubst_+output_frame-timesh, ecu_rec, old_wtda, ALDO_WINDOW, ALDO_WINDOW, output_frame );

    return;
}


/*--------------------------------------------------------------------------
 *  rec_frame()
 *
 *  Frame reconstruction
 *--------------------------------------------------------------------------*/

static void rec_frame(
    float *X,                           /* i  : FFT spectrum                          */
    float *ecu_rec,                     /* o  : Reconstructed frame in tda domain     */
    const short output_frame            /* i  : Frame length                          */
)
{
    short Lprot, LprotLog2=0;
    float fs;

    fs = (float)output_frame*50.0f;
    Lprot = 2*output_frame * L_PROT32k / 1280;

    if (output_frame == L_FRAME48k)
    {
        LprotLog2 = 9;
    }
    else if (output_frame == L_FRAME32k)
    {
        LprotLog2 = 10;
    }
    else
    {
        LprotLog2 = 9;
    }

    /* extend spectrum and IDFT */
    if (output_frame == L_FRAME48k)
    {
        ifft3( X, X, Lprot );
    }
    else
    {
        ifft_rel( X, Lprot, LprotLog2 );
    }

    rec_wtda( X, ecu_rec, output_frame, Lprot, fs );

    return;
}


/*--------------------------------------------------------------------------
 *  fir_dwn()
 *
 *
 *--------------------------------------------------------------------------*/

static void fir_dwn(
    const float x[],          /* i  : input vector                              */
    const float h[],          /* i  : impulse response of the FIR filter        */
    float y[],          /* o  : output vector (result of filtering)       */
    const short L,            /* i  : input vector size                         */
    const short K,            /* i  : order of the FIR filter (K+1 coefs.)      */
    const short decimation    /* i  : decimation                                */
)
{
    float s;
    short i, j, k, mmax;

    k = 0;

    /* do the filtering */
    for ( i = K/2; i < L; i += decimation )
    {
        s = x[i] * h[0];

        if ( i < K )
        {
            mmax = i;
        }
        else
        {
            mmax = K;
        }

        for ( j = 1; j <= mmax; j++ )
        {
            s += h[j] * x[i-j];
        }

        y[k] = s;
        k++;
    }

    for ( ; i < L+K/2; i += decimation )
    {
        s = 0;

        for (j = i-L+1; j <= K; j++)
        {
            s += h[j] * x[i-j];
        }

        y[k] = s;
        k++;
    }

    return;
}


/*--------------------------------------------------------------------------
 *  fec_ecu_pitch()
 *
 *
 *--------------------------------------------------------------------------*/

static void fec_ecu_pitch(
    const float *prevsynth,
    float *prevsynth_LP,
    const short L,
    short *N,
    float *min_corr,
    short *decimatefator,
    const short HqVoicing
)
{
    short i,filt_size;
    float accA,accB,accC,Ryy;
    short delay_ind,k,cb_start,cb_end,tmp_short,Lon20;
    float Asr_LP[FEC_DCIM_FILT_SIZE_MAX+1];

    switch(L)
    {
    case L_FRAME48k:
        *decimatefator=6;
        filt_size=60;
        mvr2r(Asr_LP48,Asr_LP,filt_size+1);
        break;
    case L_FRAME32k:
        *decimatefator=4;
        filt_size=40;
        mvr2r(Asr_LP32,Asr_LP,filt_size+1);
        break;
    case L_FRAME16k:
        *decimatefator=2;
        filt_size=20;
        mvr2r(Asr_LP16,Asr_LP,filt_size+1);
        break;
    default:
        *decimatefator=2;
        filt_size=40;
        mvr2r(Asr_LP16,Asr_LP,filt_size+1);
        break;
    }

    /* We need to inverse the ALDO window */

    /* Resampling to work at 8Khz */
    fir_dwn( prevsynth,Asr_LP,prevsynth_LP,2*L,filt_size,*decimatefator ); /* resampling without delay */

    Lon20=(short) ((L/20)/ *decimatefator);

    /* Correlation analysis */
    *min_corr=0;
    accC=0;
    for (k = 0; k < 6*Lon20; k++)
    {
        accC+=prevsynth_LP[34*Lon20+k]*prevsynth_LP[34*Lon20+k];
    }

    if (HqVoicing==1)
    {
        cb_start=0;
        cb_end=33*Lon20;
    }
    else
    {
        cb_start=0;
        cb_end=28*Lon20;
    }

    tmp_short=34*Lon20;
    accB=0;
    delay_ind=cb_start;
    for (i=cb_start; i<cb_end; i++) /* cb_end = 35 let 6 ms min of loop size */
    {
        accA=0;
        if (i==cb_start)
        {
            accB=0;
            for (k = 0; k < 6*Lon20; k++)
            {
                accA+=prevsynth_LP[i+k]*prevsynth_LP[tmp_short+k];
                accB+=prevsynth_LP[i+k]*prevsynth_LP[i+k];
            }
        }
        else
        {
            accB=accB-prevsynth_LP[i-1]*prevsynth_LP[i-1]+prevsynth_LP[i+6*Lon20-1]*prevsynth_LP[i+6*Lon20-1];
            for (k = 0; k < 6*Lon20; k++)
            {
                accA+=prevsynth_LP[i+k]*prevsynth_LP[tmp_short+k];
            }
        }

        /* tests to avoid division by zero */
        if( accB <= 0 )
        {
            accB = 1.0f;
        }

        if( accC <= 0 )
        {
            accC = 1.0f;
        }

        if( accB*accC <= 0 )
        {
            Ryy = accA;
        }
        else
        {
            Ryy = accA/(float)sqrt((accB*accC));
        }

        if( Ryy  > *min_corr  )
        {
            *min_corr = Ryy;
            delay_ind = i;
        }

        if ( HqVoicing == 0 && *min_corr > 0.95f )
        {
            break;
        }
    }

    *N = 40*Lon20-delay_ind-6*Lon20;

    return;
}


/*--------------------------------------------------------------------------
 *  fec_ecu_dft()
 *
 *
 *--------------------------------------------------------------------------*/

static void fec_ecu_dft(
    const float *prevsynth_LP, /* i: */
    const short N,             /* i: */
    float *Tfr,                /* o: */
    float *Tfi,                /* o: */
    float *sum_Tf_abs,         /* o: */
    float *Tf_abs,             /* o: */
    short *Nfft                /* o: */
)
{
    float target[2*L_FRAME8k],tmp;
    short i, Lon20,tmp_short,N_LP,k;

    Lon20 = (short)160/20;

    for (i=0; i<N; i++)
    {
        target[i] = prevsynth_LP[2*160-3*Lon20-N+i];
    }

    /* DFT  */
    *sum_Tf_abs = 0;
    *Nfft = (short)pow(2, (short)ceil(log(N)/log(2)));
    tmp =((float)N-1.0f)/((float)*Nfft-1.0f);

    set_f(Tfr,0.0f,*Nfft);
    set_f(Tfi,0.0f,*Nfft);
    Tfr[0] = target[0];
    Tfr[*Nfft-1] = target[N-1];
    for (i=1; i<*Nfft-1; i++) /* interpolation for FFT */
    {
        tmp_short = (short)floor(i*tmp);
        Tfr[i] = target[tmp_short]+((float)i*tmp-((float)tmp_short))*(target[tmp_short+1]-target[tmp_short]);
    }

    DoRTFTn(Tfr,Tfi,*Nfft);
    N_LP=(short)ceil(*Nfft/2);

    for (k=0; k<N_LP; k++)
    {
        Tf_abs[k]=(float)sqrt(Tfr[k]*Tfr[k]+Tfi[k]*Tfi[k]);
        *sum_Tf_abs+=Tf_abs[k];
    }

    return;
}

/*--------------------------------------------------------------------------*
 * singenerator()
 *
 * fast cosinus generator Amp*cos(2*pi*freq+phi)
 *--------------------------------------------------------------------------*/


static
void singenerator(
    const short L,      /* i  : size of output */
    const float cosfreq, /* i  : cosine of 1-sample dephasing at the given frequency */
    const float sinfreq, /* i  : sine   of 1-sample dephasing at the given frequency */
    const float a_re,    /* i  : real part of complex spectral coefficient at the given frequency */
    const float a_im,    /* i  : imag part of complex spectral coefficient at the given frequency */
    float xx[]           /* o  : output vector */
)
{

    float *ptr, L_C0, L_S0, L_C1, L_S1;
    float C0, S0, C1, S1;
    short i;

    L_C0 = a_re;
    S0 = a_im;

    ptr = xx;

    *ptr = *ptr + L_C0;
    ptr++;

    for (i=0; i<L/2-1; i++)
    {
        C0 = L_C0;
        L_C1 = C0*cosfreq;
        L_C1 = L_C1 - S0*sinfreq;
        L_S1 = C0 * sinfreq;
        S1 = L_S1 + S0 * cosfreq;
        *ptr = *ptr + L_C1;
        ptr++;

        C1 = L_C1;
        L_C0 = C1 * cosfreq;
        L_C0 = L_C0 - S1 * sinfreq;
        L_S0 = C1 * sinfreq;
        S0 = L_S0 + S1 * cosfreq;
        *ptr = *ptr + L_C0;
        ptr++;
    }

    C0 = L_C0;
    L_C1 = C0 * cosfreq;
    L_C1 = L_C1 - S0 * sinfreq;
    *ptr = *ptr + L_C1;
    ptr++;

    return;
}


/*--------------------------------------------------------------------------
 *  sinusoidal_synthesis()
 *
 *
 *--------------------------------------------------------------------------*/

static void sinusoidal_synthesis(
    const float *Tfr,            /* i: */
    const float *Tfi,            /* i: */
    float *Tf_abs,               /* i: */
    const short N,               /* i: */
    const short L,               /* i: */
    const short decimate_factor, /* i: */
    const short Nfft,            /* i: */
    const float sum_Tf_abs,      /* i: */
    float *synthesis,            /* o: */
    const short HqVoicing        /* i: */
)
{
    short i,k,nb_pulses,indmax = 0,nb_pulses_final;
    short pulses[FEC_MAX/2];
    float a_re[FEC_NB_PULSE_MAX], a_im[FEC_NB_PULSE_MAX], cosfreq, sinfreq;
    float freq[FEC_NB_PULSE_MAX], tmp;
    float mmax,cumsum;
    short Lon20=8;

    /* peak selection  */
    short PL,cpt;
    float old,new_s;
    short* p_pulses;
    short glued;

    p_pulses=pulses;
    nb_pulses=0;
    new_s=Tf_abs[1];
    glued=1;
    cpt=0;
    old=0;

    PL=0;
    if (N>Lon20*10 || HqVoicing )
    {
        PL=1;
    }
    while(cpt<=N/2-1-2)
    {
        if(Tf_abs[cpt]>old && Tf_abs[cpt]>new_s)
        {
            glued=cpt;

            for (i=glued; i<cpt+PL+1; i++)
            {
                *p_pulses++=i;
                nb_pulses++;
            }
            old=Tf_abs[cpt+PL];
            new_s=Tf_abs[cpt+2+PL];
            cpt=cpt+PL+1;
            glued=1;
        }
        else
        {
            old=Tf_abs[cpt];
            new_s=Tf_abs[cpt+2];
            cpt++;
            glued=0;
        }
    }


    nb_pulses_final=0;

    /* peak selection : keep the more energetics (max 20) */
    tmp=1.0f/(float)(Nfft/2);
    cumsum=0;
    for ( i=0; i<min(FEC_NB_PULSE_MAX,nb_pulses); i++)
    {
        mmax=0;
        for ( k=0; k<nb_pulses; k++)
        {
            if(Tf_abs[pulses[k]]>mmax)
            {
                mmax=Tf_abs[pulses[k]];
                indmax=pulses[k];
            }
        }
        cumsum += Tf_abs[indmax];

        if  (HqVoicing || cumsum<sum_Tf_abs*0.7f)
        {
            a_re[i] = Tfr[indmax]*tmp;
            a_im[i] = Tfi[indmax]*tmp;
            freq[i]=(float)indmax*2/((float)N);
            nb_pulses_final++;
            Tf_abs[indmax]=-1;
        }
        else
        {
            a_re[i] = Tfr[indmax]*tmp;
            a_im[i] = Tfi[indmax]*tmp;
            freq[i]=(float)indmax*2/((float)N);
            nb_pulses_final++;
            break;
        }
    }

    nb_pulses = nb_pulses_final;

    /* sinusoidal synthesis */
    set_f( synthesis, 0.0f, 40*L/20 );

    if ( HqVoicing )
    {
        k = 40*L/20;
    }
    else
    {
        k = 40*L/20;
    }

    for ( i=0; i<nb_pulses; i++)
    {
        cosfreq = (float) cos(EVS_PI*freq[i]/(float)decimate_factor);
        sinfreq = (float) sin(EVS_PI*freq[i]/(float)decimate_factor);
        singenerator(k, cosfreq, sinfreq, a_re[i], a_im[i], synthesis);
    }

    return;
}

/*--------------------------------------------------------------------------
 *  fec_noise_filling()
 *
 *
 *--------------------------------------------------------------------------*/

static void fec_noise_filling(
    const float *prevsynth,
    float *synthesis,
    const short L,
    const short N,
    const short HqVoicing,
    float *gapsynth,
    short *ni_seed_forfec
)
{
    float SS[L_FRAME48k/2],tmp;
    short Rnd_N_noise;
    short k,kk,i;
    short N_noise;
    float noisevect[34*L_FRAME48k/20];

    mvr2r( prevsynth+2*L-3*L/20-N,noisevect,N);

    /* Noise addition on full band  */
    /* residual  */
    if (N<L)
    {
        N_noise=(N/2);
        for (k=0; k<N; k++)
        {
            noisevect[k]=(noisevect[k]-synthesis[k]);
        }
    }
    else
    {
        N_noise=L/2;
        for (k=0; k<L; k++)
        {
            noisevect[k]=(noisevect[N-L+k]-synthesis[N-L+k]);
        }
    }

    if (HqVoicing)
    {
        for (i=0; i<N; i++)
        {
            noisevect[i]*=0.25f;
        }
    }

    kk=0;
    k=0;
    Rnd_N_noise=N_noise;

    while (k<2*L)
    {
        if (kk==0)
        {
            tmp=((own_random(ni_seed_forfec)/ (32768.0f)*0.2f)+0.5f);
            kk=1;
        }
        else
        {
            tmp=((own_random(ni_seed_forfec)/ (32768.0f)*0.3f)+0.7f);
            kk=0;
        }

        Rnd_N_noise=(short)((float)N_noise*tmp);

        sinq((const float) EVS_PI/(2.0f*(float)Rnd_N_noise),(const float) EVS_PI/(4.0f*(float)Rnd_N_noise),(const short)  Rnd_N_noise,SS);

        for (i=0; i<Rnd_N_noise; i++)
        {
            if (k<2*L)
            {
                synthesis[k]+=(noisevect[N_noise-Rnd_N_noise+i]*SS[i]+noisevect[N_noise+i]*SS[Rnd_N_noise-i-1]);/* *noisefact; */
            }
            k++;
        }
    }

    kk = 7*L/20;

    /* overlappadd with the ms of valid mdct of the last frame */
    sinq( EVS_PI/(6.0f*L/20.0f), EVS_PI/(12.0f*L/20.0f), 3*L/20, SS );

    for( k=0; k<3*L/20; k++ )
    {
        tmp = SS[k]*SS[k];
        synthesis[k] = prevsynth[k+37*L/20]*(1-tmp)+synthesis[k]*tmp;
    }

    mvr2r( synthesis, synthesis+kk, 2*L-kk );
    mvr2r( synthesis+L, gapsynth, L );

    mvr2r( prevsynth+2*L-3*L/20-kk, synthesis, kk );

    return;
}


/*--------------------------------------------------------------------------
 *  fec_alg()
 *
 *
 *--------------------------------------------------------------------------*/

static void fec_alg(
    const float *prevsynth,
    const float *prevsynth_LP,
    float *ecu_rec,
    const short output_frame,
    const short N,
    const short decimatefactor,
    const short HqVoicing,
    float *gapsynth,
    short *ni_seed_forfec
)
{
    short Nfft;
    float sum_Tf_abs;
    float Tfr[FEC_FFT_MAX_SIZE];
    float Tfi[FEC_FFT_MAX_SIZE];
    float Tf_abs[FEC_FFT_MAX_SIZE/2];
    float synthesis[2*L_FRAME48k];
    short n;

    fec_ecu_dft( prevsynth_LP, N, Tfr, Tfi, &sum_Tf_abs, Tf_abs, &Nfft );

    sinusoidal_synthesis( Tfr, Tfi, Tf_abs, N, output_frame, decimatefactor, Nfft, sum_Tf_abs, synthesis, HqVoicing );

    fec_noise_filling(prevsynth,synthesis,output_frame,N*decimatefactor,HqVoicing,gapsynth, ni_seed_forfec );

    n = (short)((float)output_frame*N_ZERO_MDCT_NS/FRAME_SIZE_NS);
    wtda( synthesis+(output_frame-n), ecu_rec, NULL, ALDO_WINDOW, ALDO_WINDOW, output_frame );

    return;
}


/*--------------------------------------------------------------------------
 *  hq_phase_ecu()
 *
 *  Main routine for HQ phase ECU
 *--------------------------------------------------------------------------*/

static void hq_phase_ecu(
    const float *prevsynth,            /* i  : buffer of previously synthesized signal   */
    float *ecu_rec,              /* o  : reconstructed frame in tda domain         */
    short *time_offs,            /* i/o: Sample offset for consecutive frame losses*/
    float *X_sav,                /* i/o: Stored spectrum of prototype frame        */
    short *num_p,                /* i/o: Number of identified peaks                */
    short *plocs,                /* i/o: Peak locations                            */
    float *plocsi,               /* i/o: Interpolated peak locations               */
    const float env_stab,              /* i  : Envelope stability parameter              */
    short *last_fec,             /* i/o: Flag for usage of pitch dependent ECU     */
    const short prev_bfi,              /* i  : indicating burst frame error              */
    const short old_is_transient[2],   /* i  : flags indicating previous transient frames*/
    float *mag_chg_1st,          /* i/o: per band magnitude modifier for transients*/
    float Xavg[LGW_MAX],         /* i/o: Frequency group average gain to fade to   */
    float *beta_mute,            /* o  : Factor for long-term mute                 */
    const short bwidth,                /* i  : Encoded bandwidth                         */
    const short output_frame           /* i  : frame length                              */
)
{
    short Lprot;
    float mag_chg[LGW_MAX], ph_dith, X[L_PROT48k];
    short seed;
    float alpha[LGW_MAX], beta[LGW_MAX];

    Lprot = (2*output_frame * 4)/ 5;

    if ( !prev_bfi || ( prev_bfi && *last_fec && (*time_offs == output_frame)) )
    {
        if( prev_bfi && *last_fec )
        {
            *time_offs += 0;
        }
        else
        {
            *time_offs = 0;
        }

        trans_ana( prevsynth + 2*output_frame-Lprot - *time_offs /* /2 */, mag_chg, &ph_dith, mag_chg_1st, output_frame,
                   *time_offs, env_stab, *last_fec, alpha, beta, beta_mute, Xavg /* 1.0 stable-music,  0.0 speech-like */ );
        spec_ana( prevsynth+2*output_frame-Lprot - *time_offs, plocs, plocsi, num_p, X_sav, output_frame, bwidth );

        if( prev_bfi && *last_fec )
        {
            *time_offs += output_frame;
        }
    }
    else
    {
        *time_offs += output_frame;
        if(*time_offs <= 0)
        {
            /* detect wrap around  of st->time_offs */
            *time_offs = 32767 ;   /* continued muting will ensure that the now fixed seeds are not creating tones */
        }

        trans_ana( prevsynth + 2*output_frame-Lprot, mag_chg, &ph_dith, mag_chg_1st,
                   output_frame, *time_offs, env_stab, 0, alpha, beta, beta_mute, Xavg );  /* 1.0 stable-music,  0.0 speech-like */
    }

    mvr2r( X_sav, X, Lprot );

    /* seed for own_rand2 */
    seed = *time_offs;
    if (*num_p > 0)
    {
        seed += plocs[*num_p-1];
    }

    subst_spec( plocs, plocsi, num_p, *time_offs, X, mag_chg, ph_dith, old_is_transient, output_frame, &seed, alpha, beta, *beta_mute, Xavg );

    /* reconstructed frame in tda domain */
    rec_frame( X, ecu_rec, output_frame );

    return;
}


/*--------------------------------------------------------------------------
 *  hq_ecu()
 *
 *  Main routine for HQ ECU
 *--------------------------------------------------------------------------*/

void hq_ecu(
    const float *prevsynth,            /* i  : buffer of previously synthesized signal   */
    float *ecu_rec,              /* o  : reconstructed frame in tda domain         */
    short *time_offs,            /* i/o: Sample offset for consecutive frame losses*/
    float *X_sav,                /* i/o: Stored spectrum of prototype frame        */
    short *num_p,                /* i/o: Number of identified peaks                */
    short *plocs,                /* i/o: Peak locations                            */
    float *plocsi,               /* i/o: Interpolated peak locations               */
    const float env_stab,              /* i  : Envelope stability parameter              */
    short *last_fec,             /* i/o: Flag for usage of pitch dependent ECU     */
    const short ph_ecu_HqVoicing,      /* i  : HQ Voicing flag                           */
    short *ph_ecu_active,        /* i  : Phase ECU active flag                     */
    float *gapsynth,             /* o  : Gap synthesis                             */
    const short prev_bfi,              /* i  : indicating burst frame error              */
    const short old_is_transient[2],   /* i  : flags indicating previous transient frames*/
    float *mag_chg_1st,          /* i/o: per band magnitude modifier for transients*/
    float Xavg[LGW_MAX],         /* i/o: Frequency group average gain to fade to   */
    float *beta_mute,            /* o   : Factor for long-term mute                */
    const short output_frame,          /* i  : frame length                              */
    Decoder_State *st            /* i/o: decoder state structure                   */
)
{
    short N;
    float corr = 0.0f;
    short decimatefactor;
    float prevsynth_LP[2*L_FRAME8k];

    /* find pitch and R value */
    if (!(output_frame < L_FRAME16k) )
    {
        fec_ecu_pitch( prevsynth+NS2SA(output_frame*50,ACELP_LOOK_NS/2-PH_ECU_LOOKAHEAD_NS), prevsynth_LP, output_frame, &N, &corr, &decimatefactor, ph_ecu_HqVoicing );
    }
    else
    {
        corr = 0.0f;
        decimatefactor=4;
        N=output_frame/4;
    }

    if ((st->total_brate >= 48000 && ( output_frame >= L_FRAME16k && !prev_bfi && (!old_is_transient[0] || old_is_transient[1] ) &&
                                       (ph_ecu_HqVoicing || ( ((st->env_stab_plc > 0.5) && ( corr < 0.6)) || (st->env_stab_plc < 0.5 && (corr > 0.85 )))))) ||
            (st->total_brate < 48000  && ( ( ph_ecu_HqVoicing || corr > 0.85 ) && !prev_bfi && (!old_is_transient[0] || old_is_transient[1]) )) )
    {
        fec_alg( prevsynth+NS2SA(output_frame*50,ACELP_LOOK_NS/2-PH_ECU_LOOKAHEAD_NS), prevsynth_LP, ecu_rec, output_frame, N, decimatefactor, ph_ecu_HqVoicing, gapsynth, &st->ni_seed_forfec );
        *last_fec = 1;
        *ph_ecu_active = 0;
        *time_offs = output_frame;
    }
    else
    {
        hq_phase_ecu( prevsynth, ecu_rec, time_offs, X_sav, num_p, plocs, plocsi, env_stab, last_fec,
                      prev_bfi, old_is_transient, mag_chg_1st, Xavg, beta_mute, st->bwidth, output_frame );

        *last_fec = 0;
        *ph_ecu_active = 1;
    }

    return;
}
