/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <math.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"

/*---------------------------------------------------------------------*
 * Local constants
 *---------------------------------------------------------------------*/

#define K_COR_ENC           2.857f
#define C_COR_ENC           -1.286f
#define K_EE_ENC            0.04167f
#define C_EE_ENC            0.0f
#define K_ZC_ENC            -0.04f
#define C_ZC_ENC            2.4f
#define K_RELE_ENC          0.05f
#define C_RELE_ENC          0.45f
#define K_PC_ENC            -0.07143f
#define C_PC_ENC            1.857f
#define K_SNR_ENC           0.1111f
#define C_SNR_ENC           -0.3333f

/*-------------------------------------------------------------------*
 * signal_clas()
 *
 * Classification state machine for FEC
 * Coder type modification
 *-------------------------------------------------------------------*/

short signal_clas(               /* o  : classification for current frames                 */
    Encoder_State *st,                   /* i/o: encoder state structure                           */
    short *coder_type,           /* i/o: coder type                                        */
    const float voicing[3],            /* i  : normalized correlation for 3 half-frames          */
    const float *speech,               /* i  : pointer to speech signal for E computation        */
    const short localVAD,              /* i  : vad without hangover                              */
    const short pit[3],                /* i  : open loop pitch values for 3 half-frames          */
    const float *ee,                   /* i  : lf/hf E ration for 2 half-frames                  */
    const float relE,                  /* i  : frame relative E to the long term average         */
    const short L_look,                /* i  : look-ahead                                        */
    short *uc_clas               /* o  : flag for unvoiced class used in sp/mus classifier */
)
{
    float mean_voi2, mean_ee2, tmp;
    float een, corn, zcn, relEn, pcn, fmerit1;
    short i, clas, pc, zc;
    short unmod_coder_type;

    /*----------------------------------------------------------------*
     * Calculate average voicing
     * Calculate average spectral tilt
     * Calculate zero-crossing rate
     * Calculate pitch stability
     *----------------------------------------------------------------*/

    /* average voicing on second half-frame and look-ahead  */
    mean_voi2 = 0.5f * (voicing[1] + voicing[2]);

    /* average spectral tilt in dB */
    tmp = ee[0] * ee[1];
    if( tmp < 1.0f )
    {
        tmp = 1.0f;
    }

    mean_ee2 = 0.5f * 20.0f * (float)log10( tmp );

    /* compute zero crossing rate */
    zc = 0;
    for( i = L_look; i < L_FRAME + L_look; i++ )
    {
        if( speech[i] <= 0.0f && speech[i-1] > 0.0f )
        {
            zc++;
        }
    }

    /* compute pitch stability */
    pc = (short)(abs(pit[1] - pit[0]) + abs(pit[2] - pit[1]));

    /*-----------------------------------------------------------------*
     * Transform parameters to the range <0:1>
     * Compute the merit function
     *-----------------------------------------------------------------*/

    een = K_EE_ENC * mean_ee2 + C_EE_ENC;
    if( een > 1.0f )
    {
        een = 1.0f;
    }
    else if( een < 0.0f )
    {
        een = 0.0f;
    }

    corn = K_COR_ENC * mean_voi2 + C_COR_ENC;
    if( corn > 1.0f )
    {
        corn = 1.0f;
    }
    else if( corn < 0.0f )
    {
        corn = 0.0f;
    }

    zcn = K_ZC_ENC * zc + C_ZC_ENC;
    if( zcn > 1.0f )
    {
        zcn = 1.0f;
    }
    else if( zcn < 0.0f )
    {
        zcn = 0.0f;
    }

    relEn = K_RELE_ENC * relE + C_RELE_ENC;
    if( relEn > 1.0f )
    {
        relEn = 1.0f;
    }
    else if( relEn < 0.5f )
    {
        relEn = 0.5f;
    }

    pcn = K_PC_ENC * pc + C_PC_ENC;
    if( pcn > 1.0f )
    {
        pcn = 1.0f;
    }
    else if( pcn < 0.0f )
    {
        pcn = 0.0f;
    }

    fmerit1 = (1.0f/6.0f) * (een + 2.0f*corn + zcn + relEn + pcn);

    /*-----------------------------------------------------------------*
     * FEC classification
     * Onset classification
     *-----------------------------------------------------------------*/


    /* FEC classification */
    if( localVAD == 0 || *coder_type == UNVOICED || relE < -6.0f )
    {
        clas = UNVOICED_CLAS;
    }
    else
    {
        switch( st->last_clas )
        {
        case VOICED_CLAS:
        case ONSET:
        case VOICED_TRANSITION:
            if( fmerit1 < 0.49f )
            {
                clas = UNVOICED_CLAS;
            }
            else if( fmerit1 < 0.66f )
            {
                clas = VOICED_TRANSITION;
            }
            else
            {
                clas = VOICED_CLAS;
            }

            break;

        case UNVOICED_CLAS:
        case UNVOICED_TRANSITION:
            if( fmerit1 > 0.63f )
            {
                clas = ONSET;
            }
            else if( fmerit1 > 0.585f )
            {
                clas = UNVOICED_TRANSITION;
            }
            else
            {
                clas = UNVOICED_CLAS;
            }

            break;

        default:
            clas = UNVOICED_CLAS;

            break;
        }
    }

    /* set flag for unvoiced class, it will be used in sp/mus classifier */
    *uc_clas = clas;
    if( ( ( *coder_type == UNVOICED ) ||
            ( st->input_bwidth != NB && fmerit1 < 0.41f && st->mold_corr > 0.65f ) ||              /* WB case */
            ( st->input_bwidth == NB && fmerit1 * 0.88f < 0.41f  && st->mold_corr > 0.55f ) ) &&   /* NB case */
            relE > -15.0f && st->lt_dec_thres < 1.5f )
    {
        *uc_clas = UNVOICED_CLAS;
    }

    /* Onset classification */

    /* tc_cnt == -1: frame after TC frame in continuous block of GC/VC frames */
    /* tc_cnt ==  0: UC frame */
    /* tc_cnt ==  1: onset/transition frame, coded by GC mode */
    /* tc_cnt ==  2: frame after onset/transition frame, coded by TC mode */

    if( clas == UNVOICED_CLAS )
    {
        st->tc_cnt = 0;
    }

    if( clas >= VOICED_TRANSITION && st->tc_cnt >= 0 )
    {
        st->tc_cnt += 1;
    }

    if( st->tc_cnt > 2 )
    {
        st->tc_cnt = -1;
    }

    if ( st->codec_mode == MODE1 )
    {
        /*---------------------------------------------------------------------*
         * Coder type modification
         *
         * Prevent UC mode in certain conditions
         * Prevent VC mode in certain conditions
         * Select TC mode in appropriate frames
         *---------------------------------------------------------------------*/

        /* At higher rates, use GC coding instead of UC coding to improve quality */
        if( st->total_brate > ACELP_9k60 && *coder_type == UNVOICED )
        {
            *coder_type = GENERIC;
        }

        /* Prevent UC coding on mixed content at 9.6 kb/s */
        if( st->total_brate == ACELP_9k60 && *coder_type == UNVOICED && st->audio_frame_cnt != 0 )
        {
            *coder_type = GENERIC;
        }

        unmod_coder_type = *coder_type;

        /* Enforce GC mode on inactive signal (this can be later overwritten to INACTIVE) */
        if( localVAD == 0 && ( (*coder_type == UNVOICED && (!st->Opt_SC_VBR ||
                                ( st->Opt_SC_VBR && st->vbr_generic_ho == 0 && st->last_coder_type > UNVOICED ))  )
                               || *coder_type == TRANSITION || *coder_type == VOICED )
          )
        {
            *coder_type = GENERIC;
        }
        if( *coder_type == GENERIC && unmod_coder_type == UNVOICED && st->Opt_SC_VBR )
        {
            st->vbr_generic_ho = 1;
        }

        if ( *coder_type > UNVOICED && st->Opt_SC_VBR )
        {
            st->vbr_generic_ho = 0;
        }

        if( localVAD == 0 && *coder_type == UNVOICED )
        {
            st->last_7k2_coder_type = GENERIC;
        }
        else
        {
            st->last_7k2_coder_type = *coder_type;
        }

        /* Select TC mode for appropriate frames which is in general MODE1_VOICED_TRANSITION, VOICED_CLAS or MODE1_ONSET frames following UNVOICED_CLAS frames */
        if( localVAD != 0 && st->tc_cnt >= 1 )
        {
            if ( st->tc_cnt == 1 )
            {
                /* onset/transition frame is always coded using GC mode */
                *coder_type = GENERIC;
            }
            else
            {
                /* frame after onset/transition frame is coded by TC mode */
                *coder_type = TRANSITION;
            }
        }

        /* At higher rates and with 16kHz core, allow only GC and TC mode */
        if( st->total_brate >= ACELP_24k40 && *coder_type != GENERIC && *coder_type != TRANSITION )
        {
            *coder_type = GENERIC;
        }

        /* Patch for certain low-level signals for which the gain quantizer sometimes goes out of its dynamic range */
        if( *coder_type == VOICED && st->input_bwidth == NB && relE < -10.0f && st->total_brate <= ACELP_8k00 )
        {
            *coder_type = GENERIC;
        }
    }

    return clas;
}
