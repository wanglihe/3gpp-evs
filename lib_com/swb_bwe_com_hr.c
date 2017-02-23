/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "prot.h"

/*-------------------------------------------------------------------*
 * swb_hr_nonzero_subband_noise_fill()
 *
 * SWB BWE HR noise filling of zero subbands
 *-------------------------------------------------------------------*/

static void swb_hr_nonzero_subband_noise_fill(
    const float tilt_wb,               /* i  : tilt of wideband signal      */
    float *t_audio,              /* i/o: mdct spectrum                */
    short *bwe_highrate_seed,    /* i/o: seed of random noise         */
    const short N,                     /* i  : length of subband            */
    const short Nsv                    /* i  : number of subband            */
)
{
    short i, j;
    float *ptr;
    float min_bwe, max_bwe, tmpF;

    if( tilt_wb > 5.0f )
    {
        for( i=0; i<Nsv; i++ )
        {
            min_bwe = FLT_MAX;
            max_bwe = FLT_MIN;
            for( j=0; j<N; j++ )
            {
                tmpF = (float) fabs( t_audio[i*N + j] );

                if( tmpF > 0 )
                {
                    min_bwe = min(tmpF, min_bwe);
                    max_bwe = max(tmpF, max_bwe);
                }
            }

            if( max_bwe == min_bwe && min_bwe > 1.0f )
            {
                min_bwe *= 0.5f;
            }

            ptr = &t_audio[i*N];
            for( j=0; j<N; j++ )
            {
                if( *ptr == 0 )
                {
                    if( own_random( bwe_highrate_seed ) > 0 )
                    {
                        *ptr = 0.5f*min_bwe;
                    }
                    else
                    {
                        *ptr = -0.5f*min_bwe;
                    }
                }
                ptr++;
            }
        }
    }

    return;
}

/*-------------------------------------------------------------------*
 * swb_hr_noise_fill()
 *
 * SWB BWE HR noise filling
 *-------------------------------------------------------------------*/

void swb_hr_noise_fill(
    const short is_transient,         /* i  : transient flag          */
    const short spect_start,          /* i  : spectrum start point    */
    const short spect_end,            /* i  : spectrum end point      */
    const float tilt_wb,              /* i  : tilt of wideband signal */
    const float pitch,                /* i  : pitch value             */
    const short nq[],                 /* i  : AVQ nq index            */
    short Nsv,                  /* i  : number of subband       */
    short *bwe_highrate_seed,   /* i/o: seed of random noise    */
    float *t_audio              /* i/o: mdct spectrum           */
)
{
    short i, j, k;
    float alpha, beta;
    short pos_start, pos_end, incr;

    if( is_transient )
    {
        for( i=0; i<Nsv; i++ )
        {
            if( nq[i] == 0 )
            {
                for( j=0; j<WIDTH_BAND; j++ )
                {
                    t_audio[i*WIDTH_BAND + j] = own_random( bwe_highrate_seed )/65536.0f;
                }
            }
        }
    }
    else
    {
        Nsv = (spect_end - spect_start)/WIDTH_BAND;
        alpha = 0.25f;
        if( tilt_wb > 5.0f )
        {
            beta = 0.25f;
        }
        else
        {
            if( 100 > 0.25f*pitch )
            {
                beta = 0.25f;
            }
            else
            {
                beta = 100/pitch;
            }
        }

        i = 0;
        if( nq[i] == 0 )
        {
            i = 1;
            while( i < Nsv && nq[i] == 0 )
            {
                i++;
            }

            pos_start = i;


            while(  i< Nsv && nq[i] != 0 )
            {
                i++;
            }

            pos_end = i - 1;

            if( pos_end - pos_start > pos_start )
            {
                pos_end = 2*pos_start - 1;
            }

            incr = pos_end;

            for( j = pos_start-1; j >= 0; j-- )
            {
                for( k=0; k<WIDTH_BAND; k++ )
                {
                    t_audio[j*WIDTH_BAND + k] = alpha*t_audio[incr*WIDTH_BAND + k] + beta*own_random( bwe_highrate_seed )/32768.0f;
                }

                incr = (incr < pos_start) ? pos_end : incr - 1;
            }
        }

        while( i < Nsv )
        {
            if( nq[i] == 0 )
            {
                pos_start = i;
                pos_end = i;

                for( j=pos_start; j<Nsv; j++ )
                {
                    if( nq[j] != 0 )
                    {
                        i = j;
                        pos_end = j;
                        break;
                    }
                }

                if( pos_start == pos_end )
                {
                    i = Nsv;
                    pos_end = Nsv;
                }

                incr = (pos_start-1);

                for( j = pos_end-1; j >= pos_start; j-- )
                {
                    for( k=0; k<WIDTH_BAND; k++ )
                    {
                        t_audio[j*WIDTH_BAND + k] = alpha*t_audio[incr*WIDTH_BAND + k] + beta*own_random( bwe_highrate_seed )/32768.0f;
                    }
                    incr = (incr == 0) ? (pos_start-1) : incr - 1;
                }
            }
            else
            {
                i++;
            }
        }
    }

    swb_hr_nonzero_subband_noise_fill( tilt_wb, t_audio, bwe_highrate_seed, WIDTH_BAND, Nsv );

    return;
}


/*-------------------------------------------------------------------*
 * td_postprocess()
 *
 * post processing in time domain for td/fd switching
 *-------------------------------------------------------------------*/

float td_postprocess(           /* o  : gain                        */
    float hb_synth[],     /* i/o: high-band synthesis         */
    const short input_frame,    /* i  : frame length                */
    const short last_extl       /* i  : last extension layer        */
)
{
    short i;
    short pos, ind1, ind2;
    float max_samp, gain, tmpF;
    float alpha, beta, step;

    max_samp = (float)fabs(hb_synth[0]);
    pos = 0;
    for( i = 1; i < input_frame; i++ )
    {
        if( fabs(hb_synth[i]) > max_samp )
        {
            max_samp = (float)fabs(hb_synth[i]);
            pos = i;
        }
    }

    if( pos < 160 )
    {
        gain = sum2_f( hb_synth + input_frame - 80, 80 ) + EPSILON;
    }
    else
    {
        gain = sum2_f( hb_synth, 80 ) + EPSILON;
    }

    gain = (float)sqrt(gain/80);

    ind1 = max(0, pos-40);
    ind2 = min( input_frame, pos+40 );
    tmpF = sum2_f( hb_synth + ind1, ind2 - ind1 ) + EPSILON;
    tmpF = (float)sqrt( tmpF/(ind2-ind1) );

    gain = min( 1.0f, gain/tmpF );

    if( last_extl == SWB_BWE || last_extl == FB_BWE )
    {
        for( i = ind1; i < input_frame; i++ )
        {
            hb_synth[i] *= gain;
        }
    }
    else
    {
        for( i = ind1; i < ind2; i++ )
        {
            hb_synth[i] *= gain;
        }

        if( ind2 != input_frame )
        {
            step = 0.0f;
            alpha = (gain > 0.5f) ? 1.0f : 0.5f;
            beta = (alpha - gain)/(input_frame - ind2);

            for( i = ind2; i < input_frame; i++ )
            {
                hb_synth[i] *= (gain + step);
                step += beta;
            }
        }
    }

    return( gain );
}
