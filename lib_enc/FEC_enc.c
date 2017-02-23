/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"

/*-------------------------------------------------------------------*
 * FEC_encode()
 *
 * Encoder supplementary information for FEC
 *-------------------------------------------------------------------*/

void FEC_encode(
    Encoder_State *st,               /* i/o: encoder state structure                         */
    const float *synth,            /* i  : pointer to synthesized speech for E computation */
    const short coder_type,        /* i  : type of coder                                   */
    short clas,              /* i  : signal clas for current frame                   */
    const float *fpit,             /* i  : close loop fractional pitch buffer              */
    const float *res,              /* i  : LP residual signal frame                        */
    short *last_pulse_pos,   /* i/o: Position of the last pulse                      */
    const short L_frame,           /* i  : Frame length                                    */
    const long  total_brate,       /* i  : total codec bitrate                             */
    const long  core_brate         /* i  : total codec bitrate                             */
)
{
    short tmpS, index;
    short maxi, sign = 0, tmp_FER_pitch;
    float enr_q;

    tmpS = 0;
    enr_q = 1.0f;

    if( coder_type > UNVOICED && coder_type < AUDIO && core_brate >= ACELP_11k60 )
    {
        /*-----------------------------------------------------------------*
         * encode signal class (not needed for VC mode since it is clearly voiced) (2 bits)
         *-----------------------------------------------------------------*/

        if ( coder_type != VOICED )
        {
            /* encode signal clas with 2 bits */
            if( clas == UNVOICED_CLAS )
            {
                index = 0;
            }
            else if( clas == VOICED_TRANSITION || clas == UNVOICED_TRANSITION )
            {
                index = 1;
            }
            else if( clas == VOICED_CLAS )
            {
                index = 2;
            }
            else
            {
                index = 3;
            }

            push_indice( st, IND_FEC_CLAS, index, FEC_BITS_CLS);
        }

        /*-----------------------------------------------------------------*
         * Encode frame energy (5 bits)
         *-----------------------------------------------------------------*/

        if( total_brate >= ACELP_16k40 && coder_type != TRANSITION )    /* GENERIC and VOICED frames */
        {
            /* frame energy (maximum energy per pitch period for voiced frames or mean energy per sample over 2nd halframe for unvoiced frames) */
            fer_energy( L_frame, clas, synth, fpit[(L_frame>>6)-1], &enr_q, L_frame );

            /* linearly quantize the energy in the range 0 : FEC_ENR_STEP : 96 dB */
            tmpS = (short)( 10.0 * log10( enr_q + 0.001f ) / FEC_ENR_STEP );

            if( tmpS > FEC_ENR_QLIMIT )
            {
                tmpS = FEC_ENR_QLIMIT;
            }

            if( tmpS < 0 )
            {
                tmpS = 0;
            }

            push_indice( st, IND_FEC_ENR, tmpS, FEC_BITS_ENR );
        }

        /*-----------------------------------------------------------------*
         * Encode last glottal pulse position (8 bits)
         *-----------------------------------------------------------------*/

        if( total_brate >= ACELP_32k && coder_type != TRANSITION )   /* GENERIC frames */
        {
            /* retrieve the last glottal pulse position of the previous frame */
            /* use the current pitch information to scale or not the quantization */
            tmp_FER_pitch = (short)(fpit[0]); /* take the 1st subframe pitch, since it is easier to retieve it on decoder side */


            sign = 0;
            maxi = *last_pulse_pos;
            if ( maxi < 0 )
            {
                sign = 1;
                maxi = -maxi;
            }

            if ( tmp_FER_pitch >= 128 )
            {
                maxi /= 2;
            }

            if ( maxi > 127 )
            {
                /* better not use the glottal pulse position at all instead of using a wrong pulse */
                /* can happen only with pitch > 254 and max pit = 289 and should happen very rarely */
                maxi = 0;
            }

            if( sign == 1 )
            {
                maxi += 128; /* use 8 bits (MSB represents the sign of the pulse) */
            }

            push_indice( st, IND_FEC_POS, maxi, FEC_BITS_POS );
        }

        /* find the glottal pulse position of the current frame (could be sent as extra FEC info in the next frame) */
        maxi = 0;
        if( clas >= VOICED_CLAS && total_brate >= ACELP_24k40 )
        {
            maxi = findpulse( L_frame, res, (short)(fpit[(L_frame>>6)-1]), 0, &sign );
            if ( sign == 1 )
            {
                maxi = -maxi;
            }
        }

        *last_pulse_pos = maxi;
    }
    else
    {
        *last_pulse_pos = 0;
    }

    return;
}


/*-------------------------------------------------------------------*
 * FEC_lsf_estim_enc()
 *
 * Simulates LSF estimation in case of FEC in the encoder ( only one frame erasure is considered )
 * The estimated LSF vector is then used to check LSF stability and may invoke safety-net usage in the next frame
 *-------------------------------------------------------------------*/

void FEC_lsf_estim_enc(
    Encoder_State *st,              /* i  : Encoder static memory                        */
    const short L_frame,          /* i  : length of the frame                          */
    float *lsf              /* o  : estimated LSF vector                         */
)
{
    short i;
    float alpha, lsf_mean[M];

    if( L_frame == L_FRAME )
    {
        mvr2r( UVWB_Ave, lsf_mean, M );
    }
    else
    {
        mvr2r( GEWB2_Ave, lsf_mean, M );
    }

    /*----------------------------------------------------------------------*
     * Initialize the alpha factor
     *----------------------------------------------------------------------*/

    if( st->last_coder_type == UNVOICED )
    {
        /* clearly unvoiced */
        alpha = ALPHA_UU;
    }
    else if( st->last_coder_type == AUDIO || st->clas == INACTIVE_CLAS )
    {
        alpha = 0.995f;
    }
    else if( st->clas == UNVOICED_CLAS )
    {
        /* if stable, do not flatten the spectrum in the first erased frame  */
        alpha = st->stab_fac * (1.0f - 2.0f*ALPHA_U) + 2.0f*ALPHA_U;
    }
    else if( st->clas == UNVOICED_TRANSITION )
    {
        alpha = ALPHA_UT;
    }
    else if( st->clas == VOICED_CLAS || st->clas == ONSET )
    {
        /* clearly voiced -  mild convergence to the CNG spectrum for the first three erased frames */
        alpha = ALPHA_V;
    }
    else if( st->clas == SIN_ONSET )
    {
        alpha = ALPHA_S;
    }
    else
    {
        /* long erasures and onsets - rapid convergence to the CNG spectrum */
        alpha = ALPHA_VT;
    }

    /*----------------------------------------------------------------------*
     * Extrapolate LSF vector
     *----------------------------------------------------------------------*/

    /* extrapolate the old LSF vector */
    for (i=0; i<M; i++)
    {
        /* calculate mean LSF vector */
        lsf_mean[i] = BETA_FEC * lsf_mean[i] + (1-BETA_FEC) * st->lsf_adaptive_mean[i];

        /* move old LSF vector towards the mean LSF vector */
        lsf[i] = alpha * st->lsf_old[i] + (1.0f - alpha) * lsf_mean[i];
    }

    /* check LSF stability through LSF ordering */
    if( L_frame == L_FRAME )
    {
        reorder_lsf( lsf, MODE1_LSF_GAP, M, INT_FS_12k8 );
    }
    else  /* L_frame == L_FRAME16k */
    {
        reorder_lsf( lsf, MODE1_LSF_GAP, M, INT_FS_16k );
    }

    return;
}
