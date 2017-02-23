/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "cnst.h"
#include "options.h"
#include "prot.h"
#include "rom_com.h"
#include "stat_dec.h"

/*-------------------------------------------------------------------*
 * Local constants
 *-------------------------------------------------------------------*/

#define K_COR_DEC           0.8547f         /* <-0.29, 0.88>        */
#define C_COR_DEC           0.2479f

#define K_TILT_DEC          0.8333f         /* <-0.35, 0.85>        */
#define C_TILT_DEC          0.2917f

#define K_ZC_DEC            -0.04f          /* <63, 38>             */
#define C_ZC_DEC            2.52f

#define K_ENR_DEC           0.04f           /* <-14, 11>            */
#define C_ENR_DEC           0.56f

#define K_PC_DEC            -0.0357f        /* <45, 17>             */
#define C_PC_DEC            1.6071f

/*-------------------------------------------------------------------*
 * Local functions
 *-------------------------------------------------------------------*/

static float calculate_zero_crossings( const float *synth, const short L_frame );

static float calculate_pitch_synchr_norm_correlation( const float *pitch, const float *synth, const short L_frame, const short L_subfr );

static float calculate_spectral_tilt( const float *synth, const short nbSubfr, const short L_subfr );

static short calculate_classification_result( short const last_good, float fmerit1, float const ener, short const codec_mode );

static void FEC_ClassifierCore( const float *synth, const float *pitch, short const L_frame, const short codec_mode, short *clas, float *lp_speech, float *frame_ener,
                                float *ener, float *enern, float *class_para, const float LTP_Gain, const int narrowBand, const SIGNAL_CLASSIFIER_MODE mode );

static short FEC_dec_class( Decoder_State *st, const long core_brate, const short coder_type, float *enr_q, const short last_good );

static void FEC_classificationMusic( const short coder_type, short *decision_hyst, short *clas );


/*------------------------------------------------------------------------*
 * calculate_zero_crossings()
 *
 *
 *------------------------------------------------------------------------*/

static float calculate_zero_crossings(
    const float *synth,
    const short L_frame
)
{
    short i;
    float zc_frame;

    zc_frame = 0.f;

    for ( i = 0; i <  L_frame; ++i )
    {
        if ( (synth[i] <= 0.0f) && (synth[i-1] > 0.0f) )
        {
            zc_frame += 1.0f;
        }
    }

    /* Renormalization for 12.8kHz core*/
    zc_frame *= 256.0f / (float)L_frame;

    return zc_frame;
}

/*------------------------------------------------------------------------*
 * calculate_pitch_synchr_norm_correlation()
 *
 *
 *------------------------------------------------------------------------*/

static float calculate_pitch_synchr_norm_correlation(
    const float *pitch,
    const float *synth,
    const short L_frame,
    const short L_subfr
)
{
    short T0, pos;
    short j, i;
    float cor_max[16] = {0.f};
    float enr1t, enr2t, voicing;

    T0 = (short) pitch[3];

    if ( T0 > (3 * (L_subfr / 2) ))
    {
        T0 = (short) ( 0.5f * ( pitch[2] + pitch[3]) + 0.5f );
    }

    pos = L_frame;
    j = 0;

    while ( pos > (L_frame/L_subfr==4?3*L_subfr:4*L_subfr) )
    {
        pos -= T0;

        cor_max[j] = dotp( &synth[pos], &synth[pos-T0], T0);
        enr1t = dotp( &synth[pos-T0], &synth[pos-T0], T0);
        enr2t = dotp( &synth[pos], &synth[pos], T0);

        cor_max[j] *= inv_sqrt( enr1t * enr2t + 0.1f );


        if ( (pos - T0) <  (L_frame - L_subfr) )
        {
            T0 = (short) ( 0.5f * ( pitch[2] + pitch[3] ) + 0.5f );
        }
        j++;
    }

    voicing = cor_max[0];

    for ( i = 1; i < j; ++i )
    {
        voicing += cor_max[i];
    }

    if( j > 1 )
    {
        voicing /= j;
    }

    return voicing;
}


/*------------------------------------------------------------------------*
 * calculate_spectral_tilt()
 *
 *
 *------------------------------------------------------------------------*/

static float calculate_spectral_tilt(
    const float *synth,
    const short nbSubfr,
    const short L_subfr
)
{
    short i;
    float tilt;

    float enr1t = 0.0f;
    float enr2t = 0.0f;
    float const* pt1  = synth + L_subfr;
    float const* pt2  = synth + L_subfr - 1;

    for ( i = 0; i < nbSubfr-1; ++i )
    {
        enr1t += dotp( pt1, pt1,  L_subfr );
        enr2t += dotp( pt1, pt2,  L_subfr );

        pt1 += L_subfr;
        pt2 += L_subfr;
    }

    tilt = enr2t / ( enr1t + 0.1f );

    return tilt;
}

/*------------------------------------------------------------------------*
 * calculate_classification_result()
 *
 *
 *------------------------------------------------------------------------*/
static short calculate_classification_result(
    short const last_good,
    float fmerit1,
    float const ener,
    short const codec_mode
)
{
    short result = UNVOICED_CLAS;

    switch ( last_good )
    {
    case VOICED_CLAS:
    case ONSET:
    case SIN_ONSET:
    case VOICED_TRANSITION:

        if ( fmerit1 < 0.39f )
        {
            result = UNVOICED_CLAS;
        }
        else if ( fmerit1 < 0.63f && (codec_mode == MODE2 || ener < -15.0f) )
        {
            result = VOICED_TRANSITION;
        }
        else
        {
            result = VOICED_CLAS;
        }
        break;

    case UNVOICED_CLAS:
    case UNVOICED_TRANSITION:
    case INACTIVE_CLAS: /* relevant for MODE2 only */

        if ( fmerit1 > 0.56f )
        {
            result = ONSET;
        }
        else if ( fmerit1 > 0.45f )
        {
            result = UNVOICED_TRANSITION;
        }
        else
        {
            result = UNVOICED_CLAS;
        }
        break;
    }

    return result;
}


/*------------------------------------------------------------------------*
 * FEC_ClassifierCore()
 *
 *
 *------------------------------------------------------------------------*/

static void FEC_ClassifierCore(
    const float *synth,
    const float *pitch,             /* i  : pitch values for each subframe          */
    short const L_frame,            /* i  : length of the frame                     */
    const short codec_mode,         /* i  : codec mode                              */
    short *clas,              /* i/o: frame classification                    */
    float *lp_speech,         /* i/o: long term active speech energy average  */
    float *frame_ener,        /* o  :                                         */
    float *ener,              /* o  :                                         */
    float *enern,             /* o  :                                         */
    float *class_para,        /* o  : classification para. fmerit1            */
    const float LTP_Gain,           /* i  :                                         */
    const int   narrowBand,         /* i  :                                         */
    const SIGNAL_CLASSIFIER_MODE mode    /* i  :                                    */
)
{
    float zc_frame;
    float tiltn, corn, zcn, pcn, fmerit1, tilt;
    float voicing;
    float pc = 0.0f;

    /*------------------------------------------------------------------------*
     * Compute the zero crossing rate for all subframes
     *------------------------------------------------------------------------*/

    zc_frame = calculate_zero_crossings(synth, L_frame);

    /*------------------------------------------------------------------------*
     * Compute the normalized correlation pitch-synch at the end of the frame
     *------------------------------------------------------------------------*/

    voicing = calculate_pitch_synchr_norm_correlation(pitch, synth, L_frame, L_SUBFR);

    /*------------------------------------------------------------------------*
     * Compute pitch coherence
     *------------------------------------------------------------------------*/

    pc = 0;       /* just to remove MSVC warnings */
    if (codec_mode == MODE1 || !(LTP_Gain != -1.0f && mode == CLASSIFIER_TCX) )
    {
        pc = (float) fabs( pitch[3] + pitch[2] - pitch[1] - pitch[0] ) * 256.0f / (float)L_frame;
    }

    /*------------------------------------------------------------------------*
     * Compute spectral tilt
     *------------------------------------------------------------------------*/

    tilt = calculate_spectral_tilt(synth, L_frame/L_SUBFR, L_SUBFR );

    /*------------------------------------------------------------------------*
     * Compute pitch-synchronous energy at the frame end
     *------------------------------------------------------------------------*/

    *ener = frame_energy( L_frame, pitch, synth, *lp_speech, frame_ener );

    /*------------------------------------------------------------------------*
     * transform parameters between 0 & 1
     * find unique merit function
     *------------------------------------------------------------------------*/

    *enern = K_ENR_DEC  **ener    + C_ENR_DEC;
    tiltn  = K_TILT_DEC * tilt     + C_TILT_DEC;
    corn   = K_COR_DEC  * voicing  + C_COR_DEC;
    zcn    = K_ZC_DEC   * zc_frame + C_ZC_DEC;

    if ( codec_mode == MODE2 && LTP_Gain != -1.0f && mode == CLASSIFIER_TCX )
    {
        pcn   = LTP_Gain * C_PC_DEC;
    }
    else
    {
        pcn   = K_PC_DEC * pc       + C_PC_DEC;
    }

    if( pcn > 1.0f )
    {
        pcn = 1.0f;
    }
    else if( pcn < 0.0f )
    {
        pcn = 0.0f;
    }

    fmerit1 = (1.0f / 6.0f ) * ( tiltn + 2.0f * corn + zcn + pcn + *enern );

    if( codec_mode == MODE2 && narrowBand )
    {
        fmerit1 *= 0.9f;
    }

    /*------------------------------------------------------------------------*
     * frame classification
     *------------------------------------------------------------------------*/

    if( codec_mode == MODE1 )
    {
        *class_para = fmerit1;
    }
    *clas = calculate_classification_result(*clas, fmerit1, *ener, codec_mode);

    return;
}


/*------------------------------------------------------------------------*
 * FEC_clas_estim()
 *
 * Estimation of frame class, if not available in the bitstream
 *------------------------------------------------------------------------*/

void FEC_clas_estim(
    const float *syn,
    const float *pitch,               /* i  : pitch values for each subframe                 */
    const short L_frame,              /* i  : length of the frame                            */
    const short coder_type,           /* i  : coder type                                     */
    const short codec_mode,           /* i  : codec mode                                     */
    float *mem_syn_clas_estim,  /* i/o: memory of the synthesis signal for frame class estimation */
    short *clas,                /* i/o: frame classification                           */
    float *lp_speech,           /* i/o: long term active speech energy average         */
    long  const bitrate,              /* i  : Decoded bitrate                                */
    short const Opt_AMR_WB,           /* i  : flag indicating AMR-WB IO mode                 */
    short *decision_hyst,       /* i/o: hysteresis of the music/speech decision        */
    short *locattack,           /* i/o: detection of attack (mainly to localized speech burst) */
    short *UV_cnt,              /* i/o: number of consecutives frames classified as UV */
    float *LT_UV_cnt,           /* i/o: long term consecutives frames classified as UV */
    float *Last_ener,           /* i/o: last_energy frame                              */
    short *amr_io_class,        /* i/o: classification for AMR-WB IO mode              */
    float *lt_diff_etot,        /* i/o: long-term total energy variation               */
    float *class_para,          /* o  : classification para. fmerit1                   */
    const float LTP_Gain,             /* i  :                                                */
    const int   narrowBand,           /* i  :                                                */
    const SIGNAL_CLASSIFIER_MODE mode,/* i  :                                                */
    const int   bfi,                  /* i  :                                                */
    const float preemph_fac,          /* i  :                                                */
    const int   tcxonly,              /* i  :                                                */
    const int   last_core_brate       /* i  : last core bitrate                              */
)
{
    float diff_ener;
    float mean_diff;
    short i;
    float ftmp_c, fcorr, dev;
    float frame_ener = 0.0f;
    float ener       = 0.0f;
    float enern      = 0.0f;
    float tmp;
    float old_synth[L_SYN_CLAS_ESTIM], *synth;

    /* After Rate Switching st->last_core is reset to 0. Check for last_core_brate is needed */
    if((last_core_brate == SID_1k75 || last_core_brate == ACELP_6k60 || last_core_brate == ACELP_8k85
            || last_core_brate == ACELP_12k65 || last_core_brate == ACELP_14k25 || last_core_brate == ACELP_15k85
            || last_core_brate == ACELP_18k25 || last_core_brate == ACELP_19k85 || last_core_brate == ACELP_23k05
            || last_core_brate == ACELP_23k85) && !Opt_AMR_WB && codec_mode == MODE2
            && L_frame > L_FRAME)
    {
        int oldLenClasBuff, newLenClasBuff;
        oldLenClasBuff = L_SYN_MEM_CLAS_ESTIM * L_FRAME/L_frame;
        newLenClasBuff = L_SYN_MEM_CLAS_ESTIM;
        lerp( &mem_syn_clas_estim[L_SYN_MEM_CLAS_ESTIM-oldLenClasBuff], &mem_syn_clas_estim[L_SYN_MEM_CLAS_ESTIM-newLenClasBuff], newLenClasBuff, oldLenClasBuff );
    }

    synth = old_synth + L_SYN_MEM_CLAS_ESTIM;                     /* Set pointer to current frame */
    mvr2r( mem_syn_clas_estim, old_synth, L_SYN_MEM_CLAS_ESTIM ); /* Copy old synthesis into local buffer */
    mvr2r( syn, synth, L_frame );                                 /* Copy current synthesis into local buffer */

    /* TCX outputs non-pre-speech */
    if( codec_mode == MODE2 && mode == CLASSIFIER_TCX )
    {
        tmp = syn[-1];
        preemph(synth, preemph_fac, L_frame, &tmp);
    }

    /* Do the classification only
       - MODE1: when the class is not transmitted in the bitstream
       - MODE2: on good frames (classifier is also called for bfi=1) */
    if ((codec_mode == MODE1 && ( bitrate < ACELP_11k60 || coder_type <= UNVOICED || Opt_AMR_WB)) ||
            (codec_mode == MODE2 && bfi!=1 && !tcxonly ))
    {

        /*------------------------------------------------------------------------*
          * classification decision depending on coder_type information
         *------------------------------------------------------------------------*/

        if( coder_type == UNVOICED )
        {
            *clas = UNVOICED_CLAS;
        }
        else if( coder_type == VOICED )
        {
            *clas = VOICED_CLAS;
        }
        else if( coder_type == INACTIVE && !Opt_AMR_WB)
        {
            *clas = INACTIVE_CLAS;
        }
        else
        {
            /*------------------------------------------------------------------------*
             * GC, TC and AC frames
             *------------------------------------------------------------------------*/

            FEC_ClassifierCore( synth, pitch, L_frame, codec_mode, clas, lp_speech, &frame_ener, &ener,
                                &enern, class_para, LTP_Gain, narrowBand, mode );
        }

        /*------------------------------------------------------------------------*
         * Overwrite classification decision in case of music
         *------------------------------------------------------------------------*/

        if( codec_mode == MODE1 )
        {
            FEC_classificationMusic( coder_type, decision_hyst, clas );
        }

        /*---------------------------------------------------------------------------------*
         * Measure energy on active voice frames (to improve FEC performance)
         *---------------------------------------------------------------------------------*/

        if ( *clas == VOICED_CLAS )
        {
            if ( ( codec_mode == MODE2 && coder_type == VOICED ) ||
                    ( codec_mode == MODE1 && ( Opt_AMR_WB || (coder_type != GENERIC && coder_type != TRANSITION ) ) ) )
            {
                enern = frame_energy( L_frame, pitch, synth, *lp_speech, &frame_ener );
            }

            *lp_speech = 0.99f * *lp_speech + 0.01f * frame_ener;
        }

        if( codec_mode == MODE1 )
        {
            /*---------------------------------------------------------------------------------*
             * Overwrite classification decision to UNVOICED_CLAS in case of INACTIVE frame
             *---------------------------------------------------------------------------------*/

            if( coder_type == INACTIVE && *clas != INACTIVE_CLAS )
            {
                *clas = UNVOICED_CLAS;
            }

            /*---------------------------------------------------------------------------------*
             * Classification refinement to improve noise coding (only in AMR-WB IO mode)
             *---------------------------------------------------------------------------------*/

            if( Opt_AMR_WB )
            {
                *locattack = 0;


                if( *clas == UNVOICED_CLAS
                        && coder_type != INACTIVE
                  )
                {
                    /* unvoiced signal but not silence */
                    if( *lp_speech <= 40 )
                    {
                        *UV_cnt = 16;
                    }
                    else
                    {
                        *UV_cnt -= 8;
                    }
                }
                else if (
                    coder_type != INACTIVE
                )
                {
                    /* neither unvoiced nor clean silence */
                    (*UV_cnt)++;
                }

                if( *UV_cnt > 300 )
                {
                    *UV_cnt = 300;
                }
                else if( *UV_cnt < 0 )
                {
                    *UV_cnt = 0;
                }

                /* Update long-term average number of frames between UNVOICED */
                if( coder_type == INACTIVE )
                {
                    /* tends to speech if no activity */
                    *LT_UV_cnt = 0.95f * *LT_UV_cnt;

                    if ( *UV_cnt > 125 )
                    {
                        *UV_cnt = 125;
                    }
                }
                else
                {
                    *LT_UV_cnt = 0.9f * *LT_UV_cnt + 0.1f * *UV_cnt;
                }

                /*-----------------------------------------------------------------------------*
                 * Compute frame energy difference
                 * IF long-term average is high and energy difference is relatively low
                 *   classification is overwritten to AUDIO
                 * IF energy difference > 6.0dB
                 *   consider an attack
                 *-----------------------------------------------------------------------------*/

                diff_ener = ener - *Last_ener;
                *Last_ener = ener;
                *amr_io_class = *clas;

                if( *LT_UV_cnt > LT_UV_THR && diff_ener < 12.0f )
                {
                    *amr_io_class = AUDIO_CLAS;
                }

                if( (diff_ener > 6.0f && *clas == AUDIO_CLAS) || (diff_ener > 9.0f) )
                {
                    *locattack = 1;
                }

                /*------------------------------------------------------------------------*
                 * Compute statistical deviation of long-term energy variation and
                 * overwrite classification, if needed
                 *------------------------------------------------------------------------*/

                if( coder_type != INACTIVE )
                {
                    /* calculate mean energy variation of past MAX_LT frames */
                    mean_diff = 0.0f;
                    for (i = 1; i<MAX_LT; i++)
                    {
                        mean_diff += lt_diff_etot[i-1] * INV_MAX_LT;
                        lt_diff_etot[i-1] = lt_diff_etot[i];
                    }
                    mean_diff += lt_diff_etot[i-1] * INV_MAX_LT;

                    /* find statistical deviation of the energy variation history against the last 15 frames */
                    fcorr = 0.0f;
                    for (i = MAX_LT-15; i<MAX_LT; i++)
                    {
                        ftmp_c = lt_diff_etot[i] - mean_diff;
                        fcorr += ftmp_c * ftmp_c;
                    }
                    lt_diff_etot[i-1] = diff_ener;

                    /* compute statistical deviation */
                    dev = (float)sqrt(fcorr / (MAX_LT-15));

                    /* overwrite classification, if needed */
                    if ( *amr_io_class == AUDIO_CLAS && dev > 5.0f )
                    {
                        *amr_io_class = *clas;
                        *UV_cnt =  (short)(80 + *UV_cnt*0.2f);
                    }
                }
            }
        }
    }

    /* update the memory of synthesis for frame class estimation */
    mvr2r( old_synth + L_frame, mem_syn_clas_estim, L_SYN_MEM_CLAS_ESTIM );

    return;
}


/*------------------------------------------------------------------------*
 * FEC_dec_class()
 *
 * Decode class and energy information
 *------------------------------------------------------------------------*/

static short FEC_dec_class(
    Decoder_State *st,            /* i/o: decoder state structure             */
    const long  bitrate,        /* i  : core bitrate                        */
    const short coder_type,     /* i  : coder type                          */
    float *enr_q,         /* i  : decoded energy                      */
    const short last_good       /* i  : Last good FEC classification        */
)
{
    short clas, tmpS;

    clas = ONSET;

    if( coder_type != VOICED )
    {
        /* decode the class */
        tmpS = (short)get_next_indice( st, FEC_BITS_CLS );

        if( tmpS == 0 )
        {
            clas = UNVOICED_CLAS;
        }
        else if( tmpS == 1 )
        {
            if( last_good >= VOICED_TRANSITION )
            {
                clas = VOICED_TRANSITION;
            }
            else
            {
                clas = UNVOICED_TRANSITION;
            }
        }
        else if( tmpS == 2 )
        {
            clas = VOICED_CLAS;
        }
    }
    else
    {
        clas = VOICED_CLAS;
    }

    /* decode the energy */
    if( bitrate > ACELP_13k20 && coder_type != TRANSITION && coder_type < AUDIO )
    {
        tmpS = (short)get_next_indice( st, FEC_BITS_ENR );

        /* convert from logarithmic to linear domain (the range is 0 : 1.55 : 96 dB) */
        *enr_q = (float)pow( 10.0f, ( (float)tmpS * FEC_ENR_STEP ) / 10.0f );
    }

    return clas;
}


/*------------------------------------------------------------------------*
 * FEC_classificationMusic()
 *
 *
 *------------------------------------------------------------------------*/

static void FEC_classificationMusic(
    const short  coder_type,      /* i  : coder type                              */
    short *decision_hyst,   /* i/o: hysteresis of the music/speech decision */
    short *clas             /* i/o: frame classification                    */
)
{
    if( coder_type == AUDIO )
    {
        (*decision_hyst) += 4;
    }
    else
    {
        (*decision_hyst)--;
    }

    if( coder_type == INACTIVE )
    {
        *decision_hyst -= 10;
    }

    if( *decision_hyst > 200 )
    {
        *decision_hyst = 200;
    }
    else if( *decision_hyst < 0 )
    {
        *decision_hyst = 0;
    }

    if( *decision_hyst > 16 && *clas < VOICED_CLAS && coder_type == AUDIO )
    {
        *clas = VOICED_CLAS;
    }

    return;
}


/*------------------------------------------------------------------------*
 * FEC_pos_dec()
 *
 * Decode class, energy and last glottal pulse position at higher bitrates
 * ( last glottal pulse position is encoded only in GENERIC frames )
 *------------------------------------------------------------------------*/

short FEC_pos_dec(
    Decoder_State *st,               /* i/o: decoder state structure                 */
    const short coder_type,        /* i  : coder type                              */
    const short last_good,         /* i  : last good classification                */
    short *last_pulse_pos,   /* o  : last glotal pulse position in the lost ACB */
    short *clas,             /* o  : decoded classification                  */
    float *enr_q,            /* o  : decoded energy                          */
    const long  core_brate         /* i  : decoded bitrate                         */
)
{
    short pitch_index, T0, T0_frac, T0_min, T0_max;
    short bit_pos_pitch_index, nBits;

    T0 = 0;
    if( coder_type > UNVOICED )
    {
        /* decode the clas and energy information */
        if( coder_type < AUDIO )
        {
            *clas = FEC_dec_class( st, core_brate, coder_type, enr_q, last_good );

            if( coder_type == GENERIC && *clas == VOICED_CLAS && ( last_good <= UNVOICED_CLAS || last_good == INACTIVE_CLAS) )
            {
                *clas = SIN_ONSET;
            }
        }

        if( coder_type == GENERIC && core_brate > ACELP_24k40 )
        {
            nBits = 0;
            if( coder_type != AUDIO )
            {
                /* find the number of bits */
                nBits = ACB_bits_16kHz_tbl[BIT_ALLOC_IDX_16KHZ(core_brate, coder_type, 0, 0)];
            }

            bit_pos_pitch_index = 71; /* 64 kbps WB,SWB and FB*/

            if( core_brate <= ACELP_32k )
            {
                bit_pos_pitch_index = 72; /* 32 kbp, WB*/
                if (st->bwidth > WB)
                {
                    bit_pos_pitch_index = 73; /* 32 kbp, SWB, FB*/
                }
            }
            /* retrieve the pitch index */
            pitch_index = (short)get_indice( st, bit_pos_pitch_index, nBits );

            /* decode pitch period */
            T0_min = PIT_MIN;
            T0_max = PIT_MAX;
            pit16k_Q_dec( pitch_index, 10, 1, &T0, &T0_frac, &T0_min, &T0_max, &st->BER_detect );

            /* decode last pulse position */
            *last_pulse_pos = (short)get_next_indice( st, FEC_BITS_POS );

            /* respect the sign */
            if (*last_pulse_pos >= 128)
            {
                *last_pulse_pos = -(*last_pulse_pos & 0x7F);
            }

            if ( T0 >= 128 )
            {
                *last_pulse_pos *= 2;
            }

            if( st->BER_detect )
            {
                *last_pulse_pos = 0;
            }
        }
    }

    return T0;
}
