/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "prot.h"
#include "rom_enc.h"
#include "rom_com.h"




static float interpolate_corr(          /* o  : interpolated value   */
    const float *x,       /* i  : input vector         */
    const short frac,     /* i  : fraction of lag      */
    const short frac_max  /* i  : max fraction         */
)
{
    short i;
    float s;
    const float *c1, *c2, *x1, *x2, *win;


    if ( frac_max == 6 )
    {
        win = E_ROM_inter6_1;
    }
    else
    {
        win = E_ROM_inter4_1;
    }

    x1 = &x[0];
    x2 = &x[1];
    c1 = &win[frac];
    c2 = &win[frac_max-frac];
    s = 0.0f;
    for(i=0; i<4; i++, c1+=frac_max, c2+=frac_max)
    {
        s+= (*x1--) * (*c1) + (*x2++) * (*c2);
    }


    return s;
}

static void tcx_ltp_pitch_search(
    int pitch_ol,
    int *pitch_int,
    int *pitch_fr,
    int *index,
    float *norm_corr,
    const short len,
    float *wsp,
    int pitmin,
    int pitfr1,
    int pitfr2,
    int pitmax,
    int pitres
)
{
    short i, t, t0, t1, step, fraction, t0_min, t0_max, t_min, t_max, delta;
    float temp, cor_max, cor[256], *pt_cor;


    if ( pitres == 6 )
    {
        delta = 8;
    }
    else
    {
        delta = 16;
    }

    t0_min = (short)pitch_ol - (delta>>1);
    t0_max = (short)pitch_ol + (delta>>1) - 1;

    if( t0_min < pitmin )
    {
        t0_min = pitmin;
        t0_max = t0_min + delta - 1;
    }
    if( t0_max > pitmax )
    {
        t0_max = pitmax;
        t0_min = t0_max - delta + 1;
    }
    t_min = t0_min - L_INTERPOL1;
    t_max = t0_max + L_INTERPOL1;

    pt_cor = cor;
    for ( t=t_min; t<=t_max; t++ )
    {
        *pt_cor++ = dotp( wsp, wsp-t, len );
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
    temp = dotp( wsp, wsp, len ) * dotp( wsp-t1, wsp-t1, len );
    *norm_corr = cor_max / (float)sqrt( temp + 0.1f );
    if ( t1 >= pitfr1 )
    {
        *pitch_int = t1;
        *pitch_fr = 0;
        *index = t1 - pitfr1 + ((pitfr2-pitmin)*pitres) + ((pitfr1-pitfr2)*(pitres>>1));
        return;
    }

    /*------------------------------------------------------------------*
     * Search fractional pitch with 1/4 subsample resolution.
     * search the fractions around t0 and choose the one which maximizes
     * the interpolated normalized correlation.
     *-----------------------------------------------------------------*/

    pt_cor = cor + L_INTERPOL1 - t0_min;
    t0 = t1;
    if ( t0 >= pitfr2 )
    {
        step = 2;
        fraction = 2;
    }
    else
    {
        step = 1;
        fraction = 1;
    }

    if (t0 == t0_min)        /* Limit case */
    {
        fraction = 0;
        cor_max = interpolate_corr( &pt_cor[t0], fraction, pitres );
    }
    else                     /* Process negative fractions */
    {
        t0--;
        cor_max = interpolate_corr( &pt_cor[t0], fraction, pitres );
        for ( i=(fraction+step); i<=pitres-1; i=i+step )
        {
            temp = interpolate_corr( &pt_cor[t0], i, pitres );
            if (temp > cor_max)
            {
                cor_max = temp;
                fraction = i;
            }
        }
    }

    for ( i=0; i<=pitres-1; i=i+step )     /* Process positive fractions */
    {
        temp = interpolate_corr( &pt_cor[t1], i, pitres );
        if (temp > cor_max)
        {
            cor_max = temp;
            fraction = i;
            t0 = t1;
        }
    }
    *pitch_int = t0;
    *pitch_fr = fraction;
    if ( t0 >= pitfr2 )
    {
        *index = t0*(pitres>>1) + (fraction>>1) - (pitfr2*(pitres>>1)) + ((pitfr2-pitmin)*pitres);
    }
    else
    {
        *index = t0*pitres + fraction - (pitmin*pitres);
    }


}

static void tcx_ltp_find_gain( float *speech, float *pred_speech, int L_frame, float *gain, int *gain_index )
{
    int gainbits = 2;


    /* Find gain */
    *gain = get_gain( speech, pred_speech, L_frame, NULL);

    /* Quantize gain */
    if (*gain >= 0.875f)
    {
        *gain_index = 3;  /* 1.00/2 */
    }
    else if (*gain >= 0.625f)
    {
        *gain_index = 2;  /* 0.75/2 */
    }
    else if (*gain >= 0.375f)
    {
        *gain_index = 1;  /* 0.50/2 */
    }
    else if (*gain >= 0.125f)
    {
        *gain_index = 0;  /* 0.25/2 */
    }
    else
    {
        *gain_index = -1; /* escape */
    }
    /* Dequantize gain */
    *gain = (float)(*gain_index + 1) * 0.625f/(float)(1<<gainbits);

}

void tcx_ltp_encode( int tcxltp_on,
                     int tcxOnly,
                     int tcxMode,
                     int L_frame,
                     int L_subfr,
                     float *speech,
                     float *speech_ltp,
                     float *wsp,
                     int Top,
                     int *ltp_param,
                     int *ltp_bits,
                     int *pitch_int,
                     int *pitch_fr,
                     float *gain,
                     int *pitch_int_past,
                     int *pitch_fr_past,
                     float *gain_past,
                     float *norm_corr_past,
                     int last_core,
                     int pitmin,
                     int pitfr1,
                     int pitfr2,
                     int pitmax,
                     int pitres,
                     struct TransientDetection const * pTransientDetection,
                     int SideInfoOnly,
                     float *A,
                     int lpcorder
                   )
{
    int n;
    float norm_corr=0.0f;
    float pred_speech[L_FRAME_PLUS];
    float tempFlatness;
    float maxEnergyChange;
    float buf_zir[M+L_FRAME_PLUS/4], *zir;
    float r[M+1], Aest[M+1];
    float alpha, step;

#define L_frameTCX  L_frame


    /* Reset memory if past frame is acelp */
    if ( last_core == ACELP_CORE )
    {
        *pitch_int_past = L_frameTCX;
        *pitch_fr_past = 0;
        *gain_past = 0.f;
    }

    /* By default, LTP is off */
    ltp_param[0] = 0;

    if ( tcxltp_on || SideInfoOnly )
    {
        /* Find pitch lag */
        tcx_ltp_pitch_search( Top, pitch_int, pitch_fr, &ltp_param[1], &norm_corr, L_frame, wsp, pitmin, pitfr1, pitfr2, pitmax, pitres );

        tempFlatness = GetTCXAvgTemporalFlatnessMeasure(pTransientDetection, NSUBBLOCKS, 1+min(NSUBBLOCKS, (int)ceil(0.5f+NSUBBLOCKS*(1.0f*(*pitch_int)/L_frame))));

        maxEnergyChange = GetTCXMaxenergyChange(pTransientDetection,
                                                tcxMode == TCX_10,
                                                NSUBBLOCKS, 1+min(NSUBBLOCKS, (int)ceil(0.5f+NSUBBLOCKS*(float)(*pitch_int)/(float)L_frame)));

        /* Switch LTP on */
        if ( ( tcxOnly == 0 && tcxMode == TCX_20 && norm_corr **norm_corr_past > 0.25f && tempFlatness < 3.5f ) ||
                ( tcxOnly == 1 && tcxMode == TCX_10 && max(norm_corr, *norm_corr_past) > 0.5f && maxEnergyChange < 3.5f ) ||
                /* Use LTP for lower correlation when pitch lag is big, L_frame*(1.2f-norm_corr) < *pitch_int <=> norm_corr > 1.2f-*pitch_int/L_frame */
                ( tcxOnly == 1 && norm_corr > 0.44f && L_frameTCX*(1.2f-norm_corr) < *pitch_int) ||
                ( tcxOnly == 1 && tcxMode == TCX_20 && norm_corr > 0.44f && (tempFlatness < 6.0f || (tempFlatness <  7.0f && maxEnergyChange < 22.0f)) ) )
        {
            {
                ltp_param[0] = 1;
            }
        }
    }
    if ( ltp_param[0] )
    {
        /* Find predicted signal */
        predict_signal( speech, pred_speech, *pitch_int, *pitch_fr, pitres, L_frameTCX );

        /* Find gain */
        tcx_ltp_find_gain( speech, pred_speech, L_frameTCX, gain, &ltp_param[2] );

        /* Total number of bits for LTP */
        if (ltp_param[2] + 1)   /* gain > 0 */
        {
            *ltp_bits = 12;
        }
        else   /* gain <= 0 -> turn off LTP */
        {
            ltp_param[0] = 0;
        }
    }
    if (!ltp_param[0])
    {
        /* No LTP -> set everything to zero */
        *pitch_int = L_frameTCX;
        *pitch_fr = 0;
        ltp_param[1] = 0;
        set_zero( pred_speech, L_frameTCX );
        *gain = 0.f;
        ltp_param[2] = 0;
        if ( tcxltp_on || SideInfoOnly)
        {
            *ltp_bits = 1;
        }
        else
        {
            *ltp_bits = 0;
        }
    }

    if (SideInfoOnly)
    {
        *gain = 0.f;
    }

    if ( *gain_past==0.f && *gain==0.f )
    {
        mvr2r( speech, speech_ltp, L_subfr );
    }
    else if ( *gain_past==0.f )
    {
        alpha = 0.f;
        step = 1.f/(float)(L_subfr);
        for ( n=0; n<L_subfr; n++ )
        {
            speech_ltp[n] = speech[n] - alpha **gain * pred_speech[n];
            alpha += step;
        }
    }
    else
    {
        if ( A==NULL )
        {
            int i,j;
            float s;
            for (i = 0; i <= lpcorder; i++)
            {
                s = 0.0;
                for (j = 0; j < L_frameTCX-i; j++)
                {
                    s += speech[j-L_frameTCX]*speech[j+i-L_frameTCX];
                }
                r[i] = s;
            }
            if (r[0] < 100.0f)
            {
                r[0] = 100.0f;
            }
            r[0] *= 1.0001f;
            lev_dur( Aest, r, lpcorder, NULL );
            A = Aest;
        }
        if ( *gain>0.f )
        {
            predict_signal( speech-lpcorder, buf_zir, *pitch_int, *pitch_fr, pitres, lpcorder );
        }
        else
        {
            set_f( buf_zir, 0.0f, lpcorder );
        }
        for ( n=0; n<lpcorder; n++ )
        {
            buf_zir[n] = speech_ltp[n-lpcorder] - speech[n-lpcorder] + *gain * buf_zir[n];
        }
        zir = buf_zir + lpcorder;
        set_f( zir, 0.0f, L_subfr );
        syn_filt( A, lpcorder,zir, zir, L_subfr, buf_zir, 0);
        alpha = 1.f;
        step = 1.f/(float)(L_subfr/2);
        for ( n=L_subfr/2; n<L_subfr; n++ )
        {
            zir[n] *= alpha;
            alpha -= step;
        }
        for ( n=0; n<L_subfr; n++ )
        {
            speech_ltp[n] = ( speech[n] - *gain * pred_speech[n] ) + zir[n];
        }
    }

    if ( SideInfoOnly || *gain == 0.0f)
    {
        for ( n=L_subfr; n<L_frameTCX; n++ )
        {
            speech_ltp[n] = speech[n];
        }
    }
    else
    {
        for ( n=L_subfr; n<L_frameTCX; n++ )
        {
            speech_ltp[n] = speech[n] - *gain * pred_speech[n];
        }
    }

    /* Update */
    *pitch_int_past = *pitch_int;
    *pitch_fr_past = *pitch_fr;
    *gain_past = *gain;
    *norm_corr_past = norm_corr;

}

