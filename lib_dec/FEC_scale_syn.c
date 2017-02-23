/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"

#define SCLSYN_LAMBDA           0.3f

/*-------------------------------------------------------------------*
 * FEC_scale_syn()
 *
 * Smooth speech energy evolution when recovering after erasure(s)
 *-------------------------------------------------------------------*/

void FEC_scale_syn(
    const short L_frame,          /* i  : length of the frame                     */
    short clas,             /* i/o: frame classification                    */
    const short last_good,        /* i:   last good frame classification          */
    float *synth,           /* i/o: synthesized speech at Fs = 12k8 Hz      */
    const float *pitch,           /* i:   pitch values for each subframe          */
    float enr_old,          /* i:   energy at the end of previous frame     */
    float enr_q,            /* i:   transmitted energy for current frame    */
    const short coder_type,       /* i:   coder type                              */
    const short LSF_Q_prediction, /* i  : LSF prediction mode                     */
    short *scaling_flag,    /* i/o: flag to indicate energy control of syn  */
    float *lp_ener_FEC_av,  /* i/o: averaged voiced signal energy           */
    float *lp_ener_FEC_max, /* i/o: averaged voiced signal energy           */
    const short bfi,              /* i:   current  frame BFI                      */
    const long total_brate,       /* i:   total bitrate                           */
    const short prev_bfi,         /* i:   previous frame BFI                      */
    const long last_core_brate,   /* i:   previous frame core bitrate             */
    float *exc,             /* i/o: excitation signal without enhancement   */
    float *exc2,            /* i/o: excitation signal with enhancement      */
    const float Aq[],             /* i:   LP filter coefs                         */
    float *old_enr_LP,      /* i/o: LP filter E of last good voiced frame   */
    const float *mem_tmp,         /* i:   temp. initial synthesis filter states   */
    float *mem_syn,         /* o:   initial synthesis filter states         */
    int   avoid_lpc_burst_on_recovery, /* i  : if true the excitation energy is limited if LP has big gain */
    short  force_scaling     /* i: force scaling                             */
)
{
    short i;
    float enr1, enr2, gain1, gain2;
    float scaling, ener_max, enr2_av, enr2_max;
    float  enr_LP;
    float h1[L_FRAME/2], tilt, rr0, rr1, pitch_dist, mean_pitch;
    short k;

    gain2 = 0.0f;
    gain1 = 0.0f;
    scaling = 1.0f;
    enr_LP = 0.0f;

    /*-----------------------------------------------------------------*
     * Find the synthesis filter impulse response on voiced
     *-----------------------------------------------------------------*/

    if( clas >= VOICED_TRANSITION && clas < INACTIVE_CLAS )
    {
        if( L_frame == L_FRAME )
        {
            enr_LP = enr_1_Az( Aq+(NB_SUBFR-1)*(M+1), L_SUBFR );
        }
        else  /* L_frame == L_FRAME16k */
        {
            enr_LP = enr_1_Az( Aq+(NB_SUBFR16k-1)*(M+1), L_SUBFR );
        }
    }

    /*-----------------------------------------------------------------*
     * Define when to scale the synthesis
     *-----------------------------------------------------------------*/

    if( bfi )
    {
        *scaling_flag = 1;				    /* Always check synthesis on bad frames */
    }
    else if( prev_bfi )
    {
        if( ( LSF_Q_prediction == AUTO_REGRESSIVE ) || ( LSF_Q_prediction == MOVING_AVERAGE ) )
        {
            *scaling_flag = 2;				/* Decoded LSFs affected  */
        }
        else if( coder_type != TRANSITION )
        {
            *scaling_flag = 1;				/* SN, but not TC mode - LSF still affected by the interpolation */
        }
        else
        {
            *scaling_flag = 0;				/* LSF still possibly affected due to interpolation */
        }
        scaling = 1.5f;
    }
    else
    {
        if( (LSF_Q_prediction == AUTO_REGRESSIVE) && (*scaling_flag == 2) )
        {
            *scaling_flag = 2;				/* Continue with energy control till the end of AR prediction */
        }
        else if( *scaling_flag > 0 )
        {
            (*scaling_flag)--;					/* If scaling flag was equal to 2, add one control frame to account for the LSF interpolation */
        }
        scaling = 2.0f;
    }

    /*-----------------------------------------------------------------*
    * Find the energy/gain at the end of the frame
    *-----------------------------------------------------------------*/

    fer_energy( L_frame, clas, synth, pitch[(L_frame>>6)-1], &enr2, L_frame );

    if( bfi || (total_brate == ACELP_7k20) || (total_brate == ACELP_8k00) )
    {
        /* previous frame erased and no TC frame */
        if( *scaling_flag > 0 )
        {
            enr2 += 0.01f;

            if( bfi )				/* In all bad frames, limit the gain to 1  */
            {
                gain2 = (float)sqrt( enr_old / enr2 );
                if( gain2 > 1.0f )
                {
                    gain2 = 1.0f;
                }

                /* find the energy/gain at the beginning of the frame */
                fer_energy( L_frame, clas, synth, pitch[0], &enr1, 0 );

                enr1 += 0.1f;
                gain1 = (float)sqrt( enr_old / enr1 );
                if( gain1 > 1.0f )
                {
                    gain1 = 1.0f;
                }
            }
            else					/* good frame  */
            {
                if( enr_q == 0 )          /* If E info (FEC protection bits) is not available in the bitstream */
                {
                    enr_q = enr2;

                    set_f( h1, 0.0f, L_FRAME/2 );
                    h1[0] = 1.0f;
                    syn_filt( Aq+(3*(M+1)), M, h1, h1, L_FRAME/2, h1+(M+1), 0 );
                    rr0 = dotp( h1, h1, L_FRAME/2-1 ) + 0.001f;
                    rr1 = dotp( h1, h1+1, L_FRAME/2-1 );
                    tilt = rr1 / rr0;

                    pitch_dist = 0.0f;
                    mean_pitch = pitch[0];
                    for( k=0; k<(NB_SUBFR - 1); k++ )
                    {
                        pitch_dist += abs((short)(pitch[k+1]+0.5f)-(short)(pitch[k]+0.5f));
                        mean_pitch += pitch[k+1];
                    }
                    pitch_dist /= (float)(NB_SUBFR-1);
                    mean_pitch /= (float)(NB_SUBFR);

                    if( ( tilt > 0.7f ) &&											/* HF resonnant filter */
                            ( (pitch_dist > 8.0f) || (mean_pitch < PIT_MIN) ) &&        /* pitch unstable or very short      */
                            ( (prev_bfi) || ( (coder_type == GENERIC) && (LSF_Q_prediction == AUTO_REGRESSIVE) ) ) )
                    {
                        if( enr_q > scaling * enr_old )
                        {
                            enr_q = scaling * enr_old;
                        }
                    }
                    else
                    {
                        if( (clas <= VOICED_TRANSITION) || (clas >= INACTIVE_CLAS) )
                        {
                            ener_max = *lp_ener_FEC_av;
                        }
                        else
                        {
                            ener_max = *lp_ener_FEC_max;
                        }

                        if( enr_old > ener_max )
                        {
                            ener_max = enr_old;
                        }
                        if( enr_q > scaling * ener_max )
                        {
                            enr_q = scaling * ener_max;
                        }
                    }
                }

                gain2 = (float)sqrt( enr_q / enr2 );


                /*-----------------------------------------------------------------*
                * Find the energy/gain at the beginning of the frame to ensure smooth transition after erasure(s)
                *-----------------------------------------------------------------*/

                if( ( (last_good >= VOICED_TRANSITION && last_good < INACTIVE_CLAS && (clas == UNVOICED_CLAS || clas == INACTIVE_CLAS)) ||
                        last_core_brate == SID_1k75 || last_core_brate == SID_2k40 || last_core_brate == FRAME_NO_DATA ) && prev_bfi )
                {
                    /* voiced -> unvoiced signal transition */
                    /* CNG -> active signal transition */
                    gain1 = gain2;
                }
                else
                {
                    /* find the energy at the beginning of the frame */
                    fer_energy( L_frame, clas, synth, pitch[0], &enr1, 0 );

                    enr1 += 0.1f;
                    gain1 = (float)sqrt( enr_old / enr1 );
                    if( gain1 > 1.2f )
                    {
                        /* prevent clipping */
                        gain1 = 1.2f;
                    }

                    /* prevent amplifying the unvoiced or inactive part of the frame in case an offset is followed by an onset */
                    if( clas == ONSET && gain1 > gain2 && prev_bfi )
                    {
                        gain1 = gain2;
                    }
                }

                enr2 = enr_q;     /* Set the end frame energy to the scaled energy, to be used in the lp_ener_FEC  */
            }

            /*------------------------------------------------------------------------------*
            * Smooth the energy evolution by exponentially evolving from gain1 to gain2
            *------------------------------------------------------------------------------*/

            gain2 *= ( 1.0f - AGC );
            for( i=0; i<L_frame; i++ )
            {
                gain1 = gain1 * AGC + gain2;
                exc[i] *= gain1;
                exc2[i] *= gain1;
            }

            /* smoothing is done in excitation domain, so redo synthesis */
            mvr2r( mem_tmp, mem_syn, M );
            syn_12k8( L_frame, Aq, exc2, synth, mem_syn, 1 );
        }
    }
    else
    {
        /* previous frame erased and no TC frame */
        if( prev_bfi && coder_type != TRANSITION )
        {
            enr2 += 0.01f;
            if( enr_q == 0 )          /* If E info (FEC protection bits) is not available in the bitstream */
            {
                enr_q = enr2;

                set_f( h1, 0.0f, L_FRAME/2 );
                h1[0] = 1.0f;
                syn_filt( Aq+(3*(M+1)), M, h1, h1, L_FRAME/2, h1+(M+1), 0 );
                rr0 = dotp( h1, h1, L_FRAME/2-1 ) + 0.001f;
                rr1 = dotp( h1, h1+1, L_FRAME/2-1 );
                tilt = rr1 / rr0;

                if( ( ( (total_brate == ACELP_13k20) || (total_brate == ACELP_12k85) || (total_brate == ACELP_12k15) || (total_brate == ACELP_11k60) ||
                        (total_brate == ACELP_9k60) ) &&
                        ( tilt > 0.7f ) &&											    /* HF resonnant filter */
                        ( (clas == UNVOICED_CLAS) || (clas == INACTIVE_CLAS) ) ) )	    /* unvoiced classification */
                {
                    if( enr_q > scaling * enr_old )
                    {
                        enr_q = scaling * enr_old;
                    }
                }
                else if( last_good >= VOICED_TRANSITION && last_good < INACTIVE_CLAS && clas >= VOICED_TRANSITION && clas < INACTIVE_CLAS )
                {
                    /* voiced -> voiced recovery */
                    if( (*old_enr_LP != 0.0f) && (enr_LP > 2 * *old_enr_LP) )
                    {
                        enr_q /= enr_LP;
                        enr_q *= 2* *old_enr_LP;
                    }
                    else if (avoid_lpc_burst_on_recovery && enr_LP > 20.0f)
                    {
                        enr_q *= sqrt(20.0f / enr_LP);
                    }
                }

                if( (last_good >= VOICED_TRANSITION && last_good < INACTIVE_CLAS && clas >= VOICED_TRANSITION && clas < INACTIVE_CLAS )
                        || force_scaling)
                {
                    if( enr_q > enr_old )
                    {
                        enr_q = (1 - SCLSYN_LAMBDA) * enr_old + SCLSYN_LAMBDA * enr_q;
                    }
                }
            }

            gain2 = (float)sqrt( enr_q / enr2 );

            /* do not allow E increase if enr_q index == 0 (lower end Q clipping) */
            if( enr_q < 1.1f )
            {
                if( gain2 > 1.0f )
                {
                    gain2 = 1.0f;
                }
            }
            else
            {
                if( gain2 > 1.2f )
                {
                    gain2 = 1.2f;
                }
            }

            /*-----------------------------------------------------------------*
            * Find the energy/gain at the beginning of the frame to ensure smooth transition after erasure(s)
            *-----------------------------------------------------------------*/

            if( clas == SIN_ONSET )
            {
                /* allow only slow increase */
                gain1 = 0.5f * gain2;
            }
            else if( (last_good >= VOICED_TRANSITION && last_good < INACTIVE_CLAS && (clas == UNVOICED_CLAS || clas == INACTIVE_CLAS)) ||
                     last_core_brate == SID_1k75 || last_core_brate == SID_2k40 || last_core_brate == FRAME_NO_DATA )
            {
                /* voiced -> unvoiced signal transition */
                /* CNG -> active signal transition */
                gain1 = gain2;
            }
            else
            {
                /* find the energy at the beginning of the frame */
                fer_energy( L_frame, clas, synth, pitch[0], &enr1, 0 );

                enr1 += 0.1f;
                gain1 = (float)sqrt( enr_old / enr1 );
                if( gain1 > 1.2f )
                {
                    /* prevent clipping */
                    gain1 = 1.2f;
                }
                if (avoid_lpc_burst_on_recovery && (enr_LP > 20.0f) && (enr_LP <= 2 * *old_enr_LP) && (gain1 > 1.0f))
                {
                    gain1 = 1.0f;
                }

                /* prevent amplifying the unvoiced or inactive part of the frame in case an offset is followed by an onset */
                if( clas == ONSET && gain1 > gain2 )
                {
                    gain1 = gain2;
                }
            }

            /*-----------------------------------------------------------------*
            * Smooth the energy evolution by exponentially evolving from gain1 to gain2
            *-----------------------------------------------------------------*/

            gain2 *= ( 1.0f - AGC );
            for( i=0; i<L_frame; i++ )
            {
                gain1 = gain1 * AGC + gain2;
                exc[i] *= gain1;
                exc2[i] *= gain1;
            }

            /* smoothing is done in excitation domain, so redo synthesis */
            mvr2r( mem_tmp, mem_syn, M );
            syn_12k8( L_frame, Aq, exc2, synth, mem_syn, 1 );
        }
    }

    /*-----------------------------------------------------------------*
    * Update low-pass filtered energy for voiced frames
    *-----------------------------------------------------------------*/

    if( !bfi && (clas >= VOICED_TRANSITION && clas < INACTIVE_CLAS) )
    {
        if( clas == VOICED_TRANSITION )
        {
            enr2_av = enr2;
            fer_energy( L_frame, VOICED_CLAS, synth, pitch[(L_frame>>6)-1], &enr2_max, L_frame );
        }
        else
        {
            enr2_max = enr2;
            fer_energy( L_frame, UNVOICED_CLAS, synth, pitch[(L_frame>>6)-1], &enr2_av, L_frame );
        }

        *lp_ener_FEC_av = 0.05f * enr2_av + 0.95f * *lp_ener_FEC_av;
        *lp_ener_FEC_max = 0.05f * enr2_max + 0.95f * *lp_ener_FEC_max;
    }

    /*-----------------------------------------------------------------*
     * Update the LP filter energy for voiced frames
     *-----------------------------------------------------------------*/

    if( clas >= VOICED_TRANSITION && clas < INACTIVE_CLAS )
    {
        *old_enr_LP = enr_LP;
    }


    return;
}
