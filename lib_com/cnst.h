/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#ifndef CNST_H
#define CNST_H

#include "options.h"

/*----------------------------------------------------------------------------------*
 * General constants
 *----------------------------------------------------------------------------------*/

#define MODE1                                 1
#define MODE2                                 2

#define EVS_PI                                3.14159265358979323846264338327950288f

#define PI2                                   (2*EVS_PI)
#define RANDOM_INITSEED                       21845     /* Seed for random generators */
#ifndef FLT_MIN
#define FLT_MIN                               (1.175494351e-38F)
#endif
#ifndef FLT_MAX
#define FLT_MAX                               (3.402823466e+38F)
#endif
#define TRUE                                  1
#define FALSE                                 0


#define MAX_FRAME_COUNTER                     200
#define MAX_BITS_PER_FRAME                    2560

#define ENC                                   0         /* Index for "encoder" */
#define DEC                                   1         /* Index for "decoder" */

#define NB                                    0         /* Indicator of 4 kHz bandwidth */
#define WB                                    1         /* Indicator of 8 kHz bandwidth */
#define SWB                                   2         /* Indicator of 14 kHz bandwidth */
#define FB                                    3         /* Indicator of 20 kHz bandwidth */

/* Conversion of bandwidth string into numerical constant */
#define CONV_BWIDTH(bw)                       ( !strcmp(bw, "NB") ? NB : !strcmp(bw, "WB") ? WB : !strcmp(bw, "SWB") ? SWB : !strcmp(bw, "FB") ? FB : -1)

#define L_FRAME48k                            960       /* Frame size in samples at 48kHz */
#define L_FRAME32k                            640       /* Frame size in samples at 32kHz */
#define L_FRAME16k                            320       /* Frame size in samples at 16kHz */
#define L_FRAME8k                             160       /* Frame size in samples at 8kHz */

/* Conversion of ns to samples for a given sampling frequency */
#define NS2SA(fs,x)                           (short)((((long)(fs)/100L) * ((x)/100L)) / 100000L)

#define SYNC_GOOD_FRAME                       (unsigned short) 0x6B21         /* synchronization word of a "good" frame */
#define SYNC_BAD_FRAME                        (unsigned short) 0x6B20         /* synchronization word of a "bad" frame */
#define G192_BIN0                             (unsigned short) 0x007F         /* binary "0" according to ITU-T G.192 */
#define G192_BIN1                             (unsigned short) 0x0081         /* binary "1" according to ITU-T G.192 */

#define ACTIVE_FRAME                          0xFF
#define SID_FRAME                             0xFA
#define ZERO_FRAME                            0xF0
#define FRAME_SIZE_NB                         13

#define RATE_MODE_MAX                         2  /* Number of rate mode */
#define BANDWIDTH_MODE_MAX                    2  /* Number of different bandwidth (NB/WB-FB) */

#define MIN_LOG_60dB                          0.000001f
#define MIN_LOG_VAL_60dB                      -60.0f

#define INV_LOG_2                             1.442695040888963f /* 1/log(2) */


/*----------------------------------------------------------------------------------*
 * Layers
 *----------------------------------------------------------------------------------*/

#define ACELP_CORE                            0         /* ACELP core layer                                         */
#define TCX_20_CORE                           1         /* TCX 20ms core layer                                      */
#define TCX_10_CORE                           2         /* TCX 10ms core layer                                      */
#define HQ_CORE                               3         /* HQ core layer                                            */
#define AMR_WB_CORE                           4         /* AMR-WB IO core                                           */


#define WB_TBE                                5         /* WB TBE layer (16/32/48kHz signals)                       */
#define WB_BWE                                6         /* WB BWE layer optimized for music (16/32/48kHz signals)   */

#define SWB_CNG                               7         /* SWB CNG layer (32/48kHz signals)                         */
#define SWB_TBE                               8         /* SWB TBE layer optimized for speech (32/48kHz signals)    */
#define SWB_BWE                               9         /* SWB BWE layer optimized for music (32/48kHz signals)     */
#define SWB_BWE_HIGHRATE                      10         /* SWB BWE layer optimized for highrate speech (32/48kHz)   */

#define FB_TBE                                11         /* FB TBE layer (48kHz signals)                             */
#define FB_BWE                                12        /* FB BWE layer optimized for music (48kHz)                 */
#define FB_BWE_HIGHRATE                       13        /* FB BWE layer optimized for highrate speech (48kHz)       */

#define IGF_BWE                               14        /* IGF layer for music (16.4 and 24.4kbps), 32kHz signals */

#define LP_CNG                                0         /* LP-based CNG in DTX operation */
#define FD_CNG                                1         /* FD-based CNG in DTX operation */

/*----------------------------------------------------------------------------------*
 * Bitrates
 *----------------------------------------------------------------------------------*/

#define FRAME_NO_DATA                         0         /* Frame with no data */
#define SID_1k75                              1750      /* SID at 1.75 kbps               (used only in AMR-WB IO mode   */
#define SID_2k40                              2400      /* SID at 2.40 kbps                                              */
#define PPP_NELP_2k80                         2800      /* PPP and NELP at 2.80 kbps      (used only for SC-VBR)         */
#define ACELP_5k90                            5900      /* ACELP core layer at average bitrate of 5.90 kbps (used only in SC-VBR mode)      */
#define ACELP_6k60                            6600      /* ACELP core layer at 6.60  kbps (used only in AMR-WB IO mode)  */
#define ACELP_7k20                            7200      /* ACELP core layer at 7.20  kbps                                */
#define ACELP_8k00                            8000      /* ACELP core layer at 8     kbps                                */
#define ACELP_8k85                            8850      /* ACELP core layer at 8.85  kbps (used only in AMR-WB IO mode)  */
#define ACELP_9k60                            9600      /* ACELP core layer at 9.60  kbps                                */
#define ACELP_11k60                           11600     /* ACELP core layer at 11.60 kbps (used for SWB TBE)             */
#define ACELP_12k15                           12150     /* ACELP core layer at 12.15 kbps (used for WB TBE)              */
#define ACELP_12k65                           12650     /* ACELP core layer at 12.65 kbps (used only in AMR-WB IO mode)  */
#define ACELP_12k85                           12850     /* ACELP core layer at 12.85 kbps (used for WB BWE)              */
#define ACELP_13k20                           13200     /* ACELP core layer at 13.20 kbps                                */
#define ACELP_14k25                           14250     /* ACELP core layer at 14.25 kbps (used only in AMR-WB IO mode)  */
#define ACELP_14k80                           14800     /* ACELP core layer at 14.80 kbps (used only in core switching)  */
#define ACELP_15k85                           15850     /* ACELP core layer at 15.85 kbps (used only in AMR-WB IO mode)  */
#define ACELP_16k40                           16400     /* ACELP core layer at 16.40 kbps                                */
#define ACELP_18k25                           18250     /* ACELP core layer at 18.25 kbps (used only in AMR-WB IO mode)  */
#define ACELP_19k85                           19850     /* ACELP core layer at 19.85 kbps (used only in AMR-WB IO mode)  */
#define ACELP_22k60                           22600     /* ACELP core layer at 22.60 kbps (used only in core switching)  */
#define ACELP_23k05                           23050     /* ACELP core layer at 23.05 kbps (used only in AMR-WB IO mode)  */
#define ACELP_23k85                           23850     /* ACELP core layer at 23.85 kbps (used only in AMR-WB IO mode)  */
#define ACELP_24k40                           24400     /* ACELP core layer at 24.40 kbps                                */
#define ACELP_29k00                           29000     /* ACELP core layer at 29.00 kbps (used for FB + SWB TBE)        */
#define ACELP_29k20                           29200     /* ACELP core layer at 29.20 kbps (used for SWB TBE)             */
#define ACELP_30k20                           30200     /* ACELP core layer at 30.20 kbps (used for FB + SWB BWE)        */
#define ACELP_30k40                           30400     /* ACELP core layer at 30.40 kbps (used for SWB BWE)             */
#define ACELP_32k                             32000     /* ACELP core layer at 32    kbps                                */
#define ACELP_48k                             48000     /* ACELP core layer at 48    kbps                                */
#define ACELP_64k                             64000     /* ACELP core layer at 64    kbps                                */

#define HQ_16k40                              16400     /* HQ core at 16.4 kbps   */
#define HQ_13k20                              13200     /* HQ core at 13.2 kbps */
#define HQ_24k40                              24400     /* HQ core at 24.4 kbps */
#define HQ_32k                                32000     /* HQ core at 32 kbps */
#define HQ_48k                                48000     /* HQ core at 48 kbps */
#define HQ_64k                                64000     /* HQ core at 64 kbps */
#define HQ_96k                                96000     /* HQ core at 96 kbps */
#define HQ_128k                               128000    /* HQ core at 128 kbps */

#define WB_TBE_0k35                           350       /* WB TBE layer (used only at 9.6 kbps on top of ACELP@12k8 core for 16kHz signals) */
#define WB_BWE_0k35                           350       /* WB BWE layer (used only on top of ACELP@12k8 core for 16kHz signals) */
#define WB_TBE_1k05                           1050      /* WB TBE layer (used only on top of ACELP@12k8 core for 16kHz signals) */
#define SWB_TBE_1k6                           1600      /* SWB TBE layer */
#define SWB_BWE_1k6                           1600      /* SWB BWE layer */
#define FB_TBE_1k8                            1800      /* SWB+FB TBE layer (used only for 48kHz signals) */
#define FB_BWE_1k8                            1800      /* SWB+FB BWE layer (used only for 48kHz signals) */
#define SWB_TBE_2k8                           2800      /* SWB TBE layer @32kbps */
#define FB_TBE_3k0                            3000      /* SWB+FB TBE layer @32kbps (used only for 48kHz signals) */
#define SWB_BWE_16k                           16000     /* SWB BWE layer for highrate SWB speech */

#define SIZE_BRATE_TBL                        11

#define BRATE2IDX(brate)                      ( brate == ACELP_7k20  ?  0: \
                                                brate == ACELP_8k00  ?  1 : \
                                                brate == ACELP_11k60 ?  2 : \
                                                brate == ACELP_12k15 ?  3 : \
                                                brate == ACELP_12k85 ?  4 : \
                                                brate == ACELP_13k20 ?  5 : \
                                                brate == ACELP_14k80 ?  6 : \
                                                brate == ACELP_16k40 ?  7 : \
                                                brate == ACELP_22k60 ?  8 : \
                                                brate == ACELP_24k40 ?  9 : \
                                                brate == ACELP_29k00 ? 10 : \
                                                brate == ACELP_29k20 ? 11 : \
                                                brate == ACELP_30k20 ? 12 : \
                                                brate == ACELP_30k40 ? 13 : \
                                                brate == ACELP_32k   ? 14 : \
                                                brate == ACELP_48k   ? 15 : \
                                                brate == ACELP_64k   ? 16 : \
                                                brate == HQ_96k      ? 17 : \
                                                brate == HQ_128k     ? 18 : -1 )

#define BRATE2IDX16k( brate )                 ( brate == ACELP_8k00  ? 0 :  \
                                                brate == ACELP_14k80 || brate == ACELP_16k40? 1 : \
                                                brate == ACELP_22k60 ? 2 : \
                                                brate == ACELP_24k40 ? 3 : \
                                                brate == ACELP_29k00 ? 4 : \
                                                brate == ACELP_29k20 ? 5 : \
                                                brate == ACELP_30k20 ? 6 : \
                                                brate == ACELP_30k40 ? 7 : \
                                                brate == ACELP_32k   ? 8 : \
                                                brate == ACELP_48k   ? 9 : \
                                                brate == ACELP_64k   ? 10: -1 )

/* Combine parameters into a single index (used to retrieve number of bits from bit allocation tables) */
#define LSF_BIT_ALLOC_IDX(brate, ctype)       ( 6*BRATE2IDX(brate) + (ctype) )

#define BIT_ALLOC_IDX(brate, ctype, sfrm, tc) \
                                              ( ( sfrm != -1 ? NB_SUBFR : 1 ) * \
                                              ( ( tc == -1 ? 4 : 10 ) * BRATE2IDX(brate) + (ctype == INACTIVE ? GENERIC : ctype) - 1 + (tc == -1 ? 0 : tc) )  + \
                                                ( sfrm != -1 ? sfrm/L_SUBFR : 0 ) )

#define BIT_ALLOC_IDX_16KHZ(brate, ctype, sfrm, tc) \
                                              ( ( sfrm > -1 ? NB_SUBFR16k : 1 ) * \
                                              ( ( tc == -1 ? 3 : 7 ) * BRATE2IDX16k(brate) + (ctype == TRANSITION ? 2 : (ctype == GENERIC ? 1 :0) ) + (tc == -1 ? 0 : tc) ) + \
                                                ( sfrm != -1 ? sfrm/L_SUBFR : 0 ) )

/* Combine coder_type, bandwidth, formant sharpening flag, and channel-aware flag into one indice */
#define SIG2IND(ctype, bw, sf, ca_rf)         ( ctype | (bw << 3) | (sf << 6) | (ca_rf << 7) )

#define MAX_ACELP_SIG 100

/*----------------------------------------------------------------------------------*
 * Bitstream indices
 *----------------------------------------------------------------------------------*/
#define MAX_PVQ_PUSH_IND           320     /* Maximum number of (fwd+reverse) calls to push_indices for the PVQ_range encoder */
enum
{
    IND_CORE,
    IND_PPP_NELP_MODE,
    IND_SID_TYPE,
    IND_ACELP_16KHZ,
    IND_ACELP_SIGNALLING,
    IND_MDCT_CORE,
    IND_BWE_FLAG,
    IND_HQ_SWITCHING_FLG,
    IND_LAST_L_FRAME,
    IND_VAD_FLAG,
    IND_HQ_BWIDTH,
    IND_TC_SUBFR,
    IND_LSF_PREDICTOR_SELECT_BIT          = IND_TC_SUBFR + 4,
    IND_LSF,
    IND_MID_FRAME_LSF_INDEX                = IND_LSF + 17,

    IND_ISF_0_0,
    IND_ISF_0_1,
    IND_ISF_0_2,
    IND_ISF_0_3,
    IND_ISF_0_4,
    IND_ISF_1_0,
    IND_ISF_1_1,
    IND_ISF_1_2,
    IND_ISF_1_3,
    IND_ISF_1_4,

    IND_GSC_ATTACK,
    IND_GSC_SWB_SPEECH,
    IND_NOISE_LEVEL,
    IND_HF_NOISE,
    IND_PIT_CONTR_IDX,
    IND_FEC_CLAS,
    IND_FEC_ENR,
    IND_FEC_POS,
    IND_ES_PRED,
    IND_HARM_FLAG_ACELP,
    /* ------------- Loop for alg. codebook indices at 24.4 kbps (special case) -------------- */
    TAG_ALG_CDBK_4T64_24KBIT_START,
    IND_ALG_CDBK_4T64_1_24KBIT             = TAG_ALG_CDBK_4T64_24KBIT_START,
    IND_ALG_CDBK_4T64_2_24KBIT             = TAG_ALG_CDBK_4T64_24KBIT_START,
    TAG_ALG_CDBK_4T64_24KBIT_END           = TAG_ALG_CDBK_4T64_24KBIT_START + 40,
    /* ------------------------------------------------ */

    /* ------------- ACELP subframe loop -------------- */
    TAG_ACELP_SUBFR_LOOP_START,
    IND_PITCH                              = TAG_ACELP_SUBFR_LOOP_START,
    IND_LP_FILT_SELECT                     = TAG_ACELP_SUBFR_LOOP_START,
    IND_ALG_CDBK_1T64                      = TAG_ACELP_SUBFR_LOOP_START,
    IND_ALG_CDBK_2T32                      = TAG_ACELP_SUBFR_LOOP_START,
    IND_ALG_CDBK_4T64                      = TAG_ACELP_SUBFR_LOOP_START,
    IND_ALG_CDBK_4T64_1                    = TAG_ACELP_SUBFR_LOOP_START,
    IND_ALG_CDBK_4T64_2                    = TAG_ACELP_SUBFR_LOOP_START,
    IND_ALG_CDBK_4T64_1BIT                 = TAG_ACELP_SUBFR_LOOP_START,
    IND_GAUS_CDBK_INDEX                    = TAG_ACELP_SUBFR_LOOP_START,
    IND_TILT_FACTOR                        = TAG_ACELP_SUBFR_LOOP_START,
    IND_GAIN                               = TAG_ACELP_SUBFR_LOOP_START,
    IND_GAIN_CODE                          = TAG_ACELP_SUBFR_LOOP_START,
    IND_TC_IMP_SHAPE                       = TAG_ACELP_SUBFR_LOOP_START,
    IND_TC_IMP_POS                         = TAG_ACELP_SUBFR_LOOP_START,
    IND_TC_IMP_SIGN                        = TAG_ACELP_SUBFR_LOOP_START,
    IND_TC_IMP_GAIN                        = TAG_ACELP_SUBFR_LOOP_START,
    IND_GAIN_PIT                           = TAG_ACELP_SUBFR_LOOP_START,
    IND_PIT_IDX                            = TAG_ACELP_SUBFR_LOOP_START,
    IND_AVQ_GAIN                           = TAG_ACELP_SUBFR_LOOP_START,
    IND_I                                  = TAG_ACELP_SUBFR_LOOP_START,
    IND_KV                                 = TAG_ACELP_SUBFR_LOOP_START,
    IND_NQ                                 = TAG_ACELP_SUBFR_LOOP_START,
    IND_HF_GAIN_MODIFICATION               = TAG_ACELP_SUBFR_LOOP_START,
    TAG_ACELP_SUBFR_LOOP_END               = TAG_ACELP_SUBFR_LOOP_START + 300,
    /* ------------------------------------------------ */

    IND_MEAN_GAIN2,
    IND_Y_GAIN_TMP                         = IND_MEAN_GAIN2 + 32,
    IND_Y_GAIN_HF                          = IND_Y_GAIN_TMP + 32,
    IND_HQ_VOICING_FLAG,
    IND_HQ_SWB_CLAS,
    IND_NF_IDX,
    IND_LC_MODE,
    IND_YNRM,
    IND_HQ_SWB_EXC_SP_CLAS                 = IND_YNRM + 44,
    IND_HQ_SWB_EXC_CLAS                    = IND_HQ_SWB_EXC_SP_CLAS,
    IND_SWB_FENV_HQ                        = IND_HQ_SWB_EXC_CLAS,
    IND_FB_FENV_HQ                         = IND_SWB_FENV_HQ + 5,
    IND_DELTA_ENV_HQ                       = IND_FB_FENV_HQ + 5,
    IND_HVQ_BWE_NL,
    IND_NUM_PEAKS                          = IND_HVQ_BWE_NL + 2,
    IND_POS_IDX,
    IND_FLAGN                              = IND_POS_IDX + 280,
    IND_PG_IDX,
    IND_HVQ_PEAKS                          = IND_PG_IDX + 27,
    IND_HVQ_NF_GAIN                        = IND_HVQ_PEAKS + 54,
    IND_HQ2_SWB_CLAS                       = IND_HVQ_NF_GAIN + 2,
    IND_HQ2_DENG_MODE,
    IND_HQ2_DENG_8SMODE,
    IND_HQ2_DENG_8SMODE_N0,
    IND_HQ2_DENG_8SMODE_N1,
    IND_HQ2_DENG_8SPOS,
    IND_HQ2_DENG_8SDEPTH,
    IND_HQ2_DENG_HMODE,
    IND_HQ2_DIFF_ENERGY,
    IND_HQ2_P2A_FLAGS                      = IND_HQ2_DIFF_ENERGY + 100,
    IND_HQ2_LAST_BA_MAX_BAND               = IND_HQ2_P2A_FLAGS + 60,
    IND_RC_START                           = IND_HQ2_LAST_BA_MAX_BAND + 2,
    IND_RC_END                             = IND_RC_START + MAX_PVQ_PUSH_IND ,
    IND_HVQ_PVQ_GAIN                       = IND_RC_END+1,
    IND_NOISINESS                          = IND_HVQ_PVQ_GAIN + 8,
    IND_ENERGY,
    IND_CNG_HO,
    IND_SID_BW,
    IND_CNG_ENV1,
    IND_WB_FENV,
    IND_WB_CLASS,
    IND_IG1,
    IND_IG2A,
    IND_IG2B,
    IND_NELP_FID,
    IND_DELTALAG,
    IND_POWER,
    IND_AMP0,
    IND_AMP1,
    IND_GLOBAL_ALIGNMENT,
    IND_PVQ_FINE_GAIN,
    IND_UV_FLAG,
    IND_SHB_SUBGAIN                        = IND_PVQ_FINE_GAIN + 44,
    IND_SHB_FRAMEGAIN,
    IND_SHB_ENER_SF,
    IND_SHB_RES_GS1,
    IND_SHB_RES_GS2,
    IND_SHB_RES_GS3,
    IND_SHB_RES_GS4,
    IND_SHB_RES_GS5,
    IND_SHB_VF,
    IND_SHB_LSF,
    IND_SHB_MIRROR                         = IND_SHB_LSF + 5,
    IND_SHB_GRID,
    IND_SWB_CLASS,
    IND_SWB_TENV,
    IND_SWB_FENV                           = IND_SWB_TENV + 4,
    IND_SHB_CNG_GAIN                       = IND_SWB_FENV + 4,
    IND_DITHERING,
    IND_FB_SLOPE,

    IND_HQ2_SPT_SHORTEN,
    IND_HQ2_SUBBAND_TCQ,
    IND_HQ2_SUBBAND_GAIN                   = IND_HQ2_SUBBAND_TCQ + 100,
    IND_HQ2_DUMMY                          = IND_HQ2_SUBBAND_GAIN + 20,

    IND_LAGINDICES,
    IND_NOISEG,
    IND_AUDIO_GAIN,
    IND_AUDIO_DELAY,

    /* ------------- HR SWB BWE loop -------------- */
    TAG_HR_BWE_LOOP_START                  = IND_AUDIO_DELAY + 4,
    IND_HR_IS_TRANSIENT                    = TAG_HR_BWE_LOOP_START,
    IND_HR_GAIN                            = TAG_HR_BWE_LOOP_START,
    IND_HR_ENVELOPE                        = TAG_HR_BWE_LOOP_START,
    IND_HR_HF_GAIN                         = TAG_HR_BWE_LOOP_START,
    IND_I2                                 = TAG_HR_BWE_LOOP_START,
    IND_KV2                                = TAG_HR_BWE_LOOP_START,
    IND_NQ2                                = TAG_HR_BWE_LOOP_START,
    TAG_HR_BWE_LOOP_END                    = TAG_HR_BWE_LOOP_START + 200,
    /* ------------------------------------------------ */

    IND_CORE_SWITCHING_CELP_SUBFRAME,
    IND_CORE_SWITCHING_AUDIO_DELAY         = IND_CORE_SWITCHING_CELP_SUBFRAME + 20,
    IND_CORE_SWITCHING_AUDIO_GAIN,

    IND_UNUSED,
    MAX_NUM_INDICES                        = IND_UNUSED + 127
};

/*----------------------------------------------------------------------------------*
 * Delays
 *----------------------------------------------------------------------------------*/

#define FRAME_SIZE_NS                         20000000L

#define ACELP_LOOK_NS                         8750000L
#define DELAY_FIR_RESAMPL_NS                  937500L
#define DELAY_CLDFB_NS                        1250000L

#define DELAY_SWB_TBE_12k8_NS                 1250000L
#define DELAY_SWB_TBE_16k_NS                  1125000L
#define MAX_DELAY_TBE_NS                      1312500L
#define DELAY_BWE_TOTAL_NS                    2312500L
#define DELAY_FD_BWE_ENC_12k8_NS              (DELAY_BWE_TOTAL_NS - (MAX_DELAY_TBE_NS - DELAY_SWB_TBE_12k8_NS))
#define DELAY_FD_BWE_ENC_16k_NS               (DELAY_BWE_TOTAL_NS - (MAX_DELAY_TBE_NS - DELAY_SWB_TBE_16k_NS))
#define DELAY_FD_BWE_ENC_NS                   2250000L

#define L_LOOK_12k8                           NS2SA(INT_FS_12k8, ACELP_LOOK_NS)     /* look-ahead length at 12.8kHz */
#define L_LOOK_16k                            NS2SA(INT_FS_16k, ACELP_LOOK_NS)      /* look-ahead length at 16kHz   */

/* core switching constants @16kHz */
#define SWITCH_GAP_LENGTH_NS                  6250000L          /* lenght of ACELP->HQ switching gap in ms  */
#define HQ_DELAY_COMP                         NS2SA(8000, DELAY_CLDFB_NS)
#define HQ_DELTA_MAX                          6                 /* maximum multiplication factor (==48kHz/8kHz) for core switching modules */

#define N_ZERO_MDCT_NS                        5625000L          /* number of zeros in ms for MDCT */
#define NL_BUFF_OFFSET                        12

#define N_WS2N_FRAMES                         40        /* number of frames for attenuation during the band-width switching */
#define N_NS2W_FRAMES                         20        /* number of frames for attenuation during the band-width switching */

/*----------------------------------------------------------------------------------*
 * Coder types (only for ACELP core when not running in AMR-WB IO mode)
 *----------------------------------------------------------------------------------*/

#define INACTIVE                              0         /* inactive      */
#define UNVOICED                              1         /* unvoiced      */
#define VOICED                                2         /* purely voiced */
#define GENERIC                               3         /* generic       */
#define TRANSITION                            4         /* transition    */
#define AUDIO                                 5         /* audio (GSC)   */
#define LR_MDCT                               6         /* low-rate MDCT core */

/*--------------------------------------------------*
 * Partial copy frame types (only for ACELP core )
 *--------------------------------------------------*/

#define ACELP_MODE_MAX                        4
#define RF_MODE_MAX                           4

/* TCX partial copy frame types */
#define RF_NO_DATA                            0
#define RF_TCXFD                              1
#define RF_TCXTD1                             2
#define RF_TCXTD2                             3

/* ACELP partial copy frame types */
#define RF_ALLPRED                            ACELP_MODE_MAX
#define RF_NOPRED                             ACELP_MODE_MAX + 1
#define RF_GENPRED                            ACELP_MODE_MAX + 2
#define RF_NELP                               ACELP_MODE_MAX + 3


/*--------------------------------------------------------------*
 * Frame length constants in mode 2
 *---------------------------------------------------------------*/

#define ACELP_TCX_TRANS_NS                    1250000     /* Duration of the ACELP->TCX overlap - 1.25 ms */
#define L_FRAME_MAX                           960         /* Max 20ms frame size @48kHz              */
#define L_FRAME_PLUS                          1200        /* Max frame size (long TCX frame)            */
#define L_MDCT_OVLP_MAX                       NS2SA(48000,ACELP_LOOK_NS) /* = Max mdct overlap */
#define N_TCX10_MAX                           480         /* Max size of TCX10 MDCT spectrum */
#define BITS_TEC                              1           /* number of bits for TEC */
#define BITS_TFA                              1           /* number of bits for TTF */
#define N_TEC_TFA_SUBFR                       16          /* number of subframes of TEC/TFA */
#define L_TEC_TFA_SUBFR16k                    (L_FRAME16k/N_TEC_TFA_SUBFR)       /* TEC/TFA subframe size @ 16kHz*/
#define MAX_TEC_SMOOTHING_DEG                 6           /* max degree of smoothing for TEC */
#define N_MAX                                 1200        /* Max size of MDCT spectrum = 25ms @ 48kHz */
#define N_MAX_TCX                             1000        /* Max size of TCX/IGF coded lines */
#define IGF_START_MN                          164         /* MDCT lines not used by IGF*/
#define IGF_START_MX                          800         /* max. MDCT lines used by IGF*/

#define NUM_DCT_LENGTH                        24

#define NB_DIV                                2           /* number of division (frame) per 20ms frame  */
#define L_MDCT_HALF_OVLP_MAX                  (L_MDCT_OVLP_MAX/2)   /* Size of the MDCT half overlap @ 48 kHz */
#define L_MDCT_MIN_OVLP_MAX                   60                    /* Size of the MDCT minimal overlap @ 48 kHz - 1.25ms */
#define L_MDCT_TRANS_OVLP_MAX                 NS2SA(48000, ACELP_TCX_TRANS_NS) /* Size of the ACELP->MDCT transition overlap - 1.25ms */

#define L_NEXT_MAX_16k                        NS2SA(16000, ACELP_LOOK_NS) /* 140 */       /* maximum encoder lookahead at 16kHz        */
#define L_NEXT_MAX_32k                        NS2SA(32000, ACELP_LOOK_NS) /* 280 */       /* maximum encoder lookahead at 32kHz        */
#define L_PAST_MAX_32k                        360         /* maximum encoder past samples at 32kHz     */

/*----------------------------------------------------------------------------------*
 * ACELP core constants
 *----------------------------------------------------------------------------------*/

#define SAFETY_NET                            0
#define MOVING_AVERAGE                        1
#define AUTO_REGRESSIVE                       2

#define INT_FS_12k8                           12800.0f                          /* internal sampling frequency                */
#define M                                     16                                /* order of the LP filter @ 12.8kHz           */
#define L_FRAME                               256                               /* frame size at 12.8kHz                      */
#define NB_SUBFR                              4                                 /* number of subframes per frame              */
#define L_SUBFR                               (L_FRAME/NB_SUBFR)                /* subframe size                              */

#define L_INP_MEM                             (L_LOOK_16k + ((L_LP_16k - (NS2SA(INT_FS_16k, ACELP_LOOK_NS) + L_SUBFR16k/2)) - 3*L_SUBFR16k/2))            /* length of memory of input signal, given by the Look-Ahead + the past memory (max needed for the LP window at 16 kHz) */
#define L_INP_12k8                            (L_INP_MEM + L_FRAME)             /* length of input signal buffer @12.8kHz */
#define L_INP                                 (L_INP_MEM + L_FRAME32k)          /* length of input signal buffer @32kHz */

#define L_EXC_MEM                             L_FRAME16k                        /* length of memory of excitation signal @16kHz     */
#define L_EXC_MEM_12k8                        (PIT_MAX + L_INTERPOL)            /* length of memory of excitation signal @12.8kHz   */
#define L_EXC_MEM_DEC                         (3*L_FRAME16k/2)                  /*Half-frame needed for PLC in case of TCX->ACELP*/
#define L_EXC                                 (L_EXC_MEM + L_FRAME16k + 1)      /* length of encoder excitation signal buffer @16kHz*/
#define L_EXC_DEC                             (L_EXC_MEM_DEC + L_FRAME16k + 1 + L_SUBFR)  /* length of decoder excitation signal buffer @16kHz*/
#define L_SYN_MEM                             NS2SA(48000,DELAY_CLDFB_NS)       /* synthesis memory length, 1.25ms @ 48kHz          */
#define L_SYN                                 (L_SYN_MEM + L_FRAME16k)          /* length of synthesis signal buffer @16kHz         */
#define L_WSP_MEM                             (PIT_MAX + L_INTERPOL)            /* length of memory for weighted input signal @12.8kHz*/
#define L_WSP                                 (L_WSP_MEM + L_FRAME + L_LOOK_12k8) /* length of weighted input signal buffer @12.8kHz*/

#define OLD_SYNTH_SIZE_DEC                    (2*L_FRAME_MAX)                   /* decoder past synthesis; needed for LTP, PLC and rate switching*/
#define OLD_SYNTH_INTERNAL_DEC                (2*L_FRAME32k)                    /* decoder past synthesis @ internal sampling rate; needed for LTP, PLC and rate switching*/
#define OLD_SYNTH_SIZE_ENC                    L_FRAME32k+L_FRAME32k/4           /* encoder synth memory */
#define OLD_EXC_SIZE_DEC                      (3*L_FRAME_MAX/2+2*L_FIR_FER2)    /*old excitation needed for decoder for PLC*/

#define TILT_CODE                             0.3f                              /* ACELP code preemphasis factor               */

#define L_SUBFR16k                            (L_FRAME16k/NB_SUBFR)             /* subframe size at 16kHz   */
#define L_HALFR16k                            (2*L_SUBFR16k)                    /* half-frame size at 16kHz */

#define L_INTERPOL2                           16        /* Length of filter for interpolation         */
#define L_INTERPOL                            (L_INTERPOL2+1) /* Length of filter for interpolation   */
#define TILT_FAC                              0.68f     /* tilt factor (denominator)                  */
#define M16k                                  20        /* order of the LP filter @ 16kHz             */
#define PIT_SHARP                             0.85f     /* pitch sharpening factor */
#define PIT_UP_SAMP                           4         /* upsampling factor for 1/4 interpolation filter */
#define PIT_L_INTERPOL2                       16
#define PIT_FIR_SIZE2                         (PIT_UP_SAMP*PIT_L_INTERPOL2+1)
#define PIT_UP_SAMP6                          6
#define PIT_L_INTERPOL6_2                     17
#define PIT_FIR_SIZE6_2                       (PIT_UP_SAMP6*PIT_L_INTERPOL6_2+1)
#define E_MIN                                 0.0035f   /* minimum allowable energy */
#define STEP_DELTA                            0.0625f   /* quantization step for tilt compensation of gaussian cb. excitation */
#define GAMMA_EV                              0.92f     /* weighting factor for core synthesis error weighting */
#define FORMANT_SHARPENING_NOISE_THRESHOLD    21.0f     /* lp_noise level above which formant sharpening is deactivated */
#define BWD_N_BINS_MAX                        13
#define LP_NOISE_THRESH                       20.f

#define L_FILT_UP8k                           24        /* Resampling - delay of filter for  8 kHz output signals (at 12.8 kHz sampling rate) */
#define LEN_WIN_SSS                           120
#define L_FILT                                12        /* Resampling - delay of the resampling low-pass filter @12.8kHz                  */
#define L_FILT16k                             15        /* Resampling - delay of filter for 16 kHz input signals (at 16kHz sampling rate) */
#define L_FILT32k                             30        /* Resampling - delay of filter for 32 kHz input signals (at 32kHz sampling rate) */
#define L_FILT48k                             45        /* Resampling - delay of filter for 48 kHz input signals (at 48kHz sampling rate) */
#define L_FILT_UP16k                          12        /* Resampling - delay of filter for 16 kHz output signals (at 12.8 kHz sampling rate) */
#define L_FILT_UP32k                          12        /* Resampling - delay of filter for 32 kHz output signals (at 12.8 kHz sampling rate) */
#define L_FILT_UP48k                          12        /* Resampling - delay of filter for 48 kHz output signals (at 12.8 kHz sampling rate) */
#define L_FILT_MAX                            L_FILT48k /* Resampling - maximum length of all filters - for memories */
#define RS_INV_FAC                            0x8000    /* Resampling - flag needed in rom_com and modif_fs to allow pre-scaled and non pre-scaled filters */

#define CLDFB_NO_CHANNELS_MAX                 60        /* CLDFB resampling - max number of CLDFB channels */
#define CLDFB_NO_COL_MAX                      16        /* CLDFB resampling - max number of CLDFB col. */
#define CLDFB_NO_COL_MAX_SWITCH               6         /* CLDFB resampling - max number of CLDFB col. for switching */
#define CLDFB_NO_COL_MAX_SWITCH_BFI           8         /* CLDFB resampling - max number of CLDFB col. for switching, BFI */
#define INV_CLDFB_BANDWIDTH                   (1.f/800.f)

typedef enum
{
    CLDFB_ANALYSIS,
    CLDFB_SYNTHESIS
} CLDFB_TYPE;

#define L_FFT                                 256       /* Spectral analysis - length of the FFT */
#define LOG2_L_FFT                            8         /* Spectral analysis - log2 of L_FFT */

#define BIN                                   (INT_FS_12k8/L_FFT)/* Spectral analysis - Width of one frequency bin in Hz */
#define NB_BANDS                              20        /* Spectral analysis - number of frequency bands */
#define VOIC_BINS                             74        /* Spectral analysis - max number of frequency bins considered as voiced (related to VOIC_BAND and L_FFT) */
#define VOIC_BAND                             17        /* Spectral analysis - number of critical bands considered as voiced (related to VOIC_BINS) */
#define VOIC_BINS_8k                          115       /* Spectral analysis - max number of frequency bins considered as voiced in NB (related to VOIC_BAND_8k and L_FFT) */
#define VOIC_BAND_8k                          17        /* Spectral analysis - number of critical bands considered as voiced in NB (related to VOIC_BINS_8k) */

#define M_ALPHA                               0.9f      /* Multi-harm analysis - forgetting factor of LT correlation map */
#define M_GAMMA                               0.99f     /* Multi-harm analysis - forgetting factor of active speech decision predictor */
#define THR_CORR                              56        /* Multi-harm analysis - starting threshold of multi-harm. correlation */

#define L_LP                                  320       /* LP analysis - LP window size */
#define L_LP_16k                              400       /* LP analysis @16kHz - LP window size for 16kHz */
#define L_LP_AMR_WB                           384       /* LP analysis - windows size (only for AMR-WB IO mode) */
#define GRID50_POINTS                         51        /* LP analysis - half-number of points to evaluate Chebyshev polynomials used in the LP coefs. conversion */
#define GRID40_POINTS                         41        /* LP analysis - half-number of points to evaluate Chebyshev polynomials used in the LP coefs. conversion */
#define GRID100_POINTS                        100       /* LP analysis - number of points to evaluate Chebyshev polynomials */

#define PIT_MIN                               34        /* OL pitch analysis - Minimum pitch lag       */
#define PIT_MAX                               231       /* OL pitch analysis - Maximum pitch lag                          */
#define PIT_MIN_EXTEND                        20        /* OL pitch analysis - Minimum pitch lag of extended range     */
#define PIT_MIN_DOUBLEEXTEND                  17        /* OL pitch analysis - Minimum pitch lag of double-extended range     */
#define OPL_DECIM                             2         /* OL pitch analysis - decimation factor */
#define L_INTERPOL1                           4         /* OL pitch analysis - interval to compute normalized correlation */
#define FIR_SIZE1                             (PIT_UP_SAMP*L_INTERPOL1+1) /* OL pitch analysis - total length of the 1/4 interpolation filter */

#define PIT_MIN_SHORTER                       29        /* OL pitch analysis - minimum for wider pitch */

#define PIT_MIN_12k8                          29        /* Minimum pitch lag with resolution 1/4      */
#define PIT_FR2_12k8                          121       /* Minimum pitch lag with resolution 1/2      */
#define PIT_FR1_12k8                          154       /* Minimum pitch lag with resolution 1        */
#define PIT_MAX_12k8                          231       /* Maximum pitch lag                          */
#define PIT_FR1_8b_12k8                       82        /* Minimum pitch lag with resolution 1 for low bit-rate pitch delay codings*/
#define PIT_MIN_16k                           36
#define PIT_FR2_16k                           36
#define PIT_FR1_16k                           165
#define PIT_FR1_8b_16k                        165
#define PIT_MIN_25k6                          58
#define PIT_FR2_25k6                          58
#define PIT_FR1_25k6                          164
#define PIT_MAX_25k6                          463
#define PIT_FR1_8b_25k6                       164
#define PIT_MIN_32k                           72
#define PIT_FR2_32k                           72
#define PIT_FR1_32k                           75
#define PIT_MAX_32k                           577
#define PIT_FR1_8b_32k                        75
#define PIT_MAX_MAX                           PIT_MAX_32k

#define PIT_FR1_8b                            92        /* Pitch encoding - Minimum pitch lag with resolution 1        */
#define PIT_FR2_9b                            128       /* Pitch encoding - Minimum pitch lag with resolution 1/2      */
#define PIT_FR1_9b                            160       /* Pitch encoding - Minimum pitch lag with resolution 1        */
#define PIT_FR1_EXTEND_8b                     64        /* Pitch encoding - Minimum pitch lag with resolution 1 of extended range       */
#define PIT_FR2_EXTEND_9b                     116       /* Pitch encoding - Minimum pitch lag with resolution 1/2 of extended range     */
#define PIT_FR1_EXTEND_9b                     128       /* Pitch encoding - Minimum pitch lag with resolution 1 of extended range       */
#define PIT_FR1_DOUBLEEXTEND_8b               58        /* Pitch encoding - Minimum pitch lag with resolution 1 of double-extended range       */
#define PIT_FR2_DOUBLEEXTEND_9b               112       /* Pitch encoding - Minimum pitch lag with resolution 1/2 of double-extended range     */
#define PIT_FR1_DOUBLEEXTEND_9b               124       /* Pitch encoding - Minimum pitch lag with resolution 1 of double-extended range       */

#define LOW_PASS                              0         /* LP filtering - flag for low-pass filtering of the excitation */
#define FULL_BAND                             1         /* LP filtering - flag for no low-pass filtering of the excitation */
#define NORMAL_OPERATION                      2         /* LP filtering - flag for selecting the best of the two above */

#define NB_TRACK_FCB_2T                       2         /* Algebraic codebook - number of tracks in algebraic fixed codebook search with 2 tracks */
#define NB_POS_FCB_2T                         32        /* Algebraic codebook - number of positions in algebraic fixed codebook search with 2 tracks */
#define NB_TRACK_FCB_4T                       4         /* Algebraic codebook - number of tracks in algebraic fixed codebook search with 4 tracks */
#define NB_POS_FCB_4T                         16        /* Algebraic codebook - number of positions in algebraic fixed codebook search with 4 tracks */
#define NB_PULSE_MAX                          36
#define NPMAXPT                               ((NB_PULSE_MAX+NB_TRACK_FCB_4T-1)/NB_TRACK_FCB_4T)
#define MAX_IDX_LEN                           9

#define GAIN_PRED_ORDER                       4         /* Gain quantization - prediction order for gain quantizer (only for AMR-WB IO mode) */
#define MEAN_ENER                             30        /* Gain quantization - average innovation energy */

#define DTX_HIST_SIZE                         8         /* CNG & DTX - number of last signal frames used for CNG averaging */
#define CNG_ISF_FACT                          0.9f      /* CNG & DTX - CNG spectral envelope smoothing factor */
#define STEP_AMR_WB_SID                       2.625f    /* CNG & DTX - CNG energy quantization step */
#define HO_HIST_SIZE                          8         /* CNG & DTX - maximal number of hangover frames used for averaging */
#define NUM_ENV_CNG                           20
#define BUF_L_NRG                             0.7f      /* CNG & DTX - lower threshold factor for hangover updates */
#define BUF_H_NRG                             1.03f     /* CNG & DTX - higher threshold factor for hangover updates */

#define BUF_DEC_RATE                          25        /* CNG & DTX - buffer size decrease rate for active frames */
#define STEP_SID                              5.25f     /* CNG & DTX - CNG energy quantization step */

#define MIN_ACT_CNG_UPD                       20        /* DTX - Minimum number of consecutive active frames for CNG mode update */
#define FIXED_SID_RATE                        8         /* DTX SID rate */

#define TOTALNOISE_HIST_SIZE                  4

#define UNKNOWN_NOISE                         0        /* unknown noisy type  */
#define SILENCE                               1        /* speech with high SNR  */
#define CLDFBVAD_NB_ID                        1
#define CLDFBVAD_WB_ID                        2
#define CLDFBVAD_SWB_ID                       3
#define CLDFBVAD_FB_ID                        4
#define SP_CENTER_NUM                         4        /* number of spectral centroid  */
#define TONA_NUM                              3        /* number of tonal  */
#define PRE_SNR_NUM                           32       /* number of snr to calculate average SNR of all sub-bands  */
#define POWER_NUM                             56       /* number of energy of several frames*/
#define PRE_SPEC_DIF_NUM                      56       /* number of energy of several frames*/

#define MAX_SUBBAND_NUM                       12       /* max number of sub-band divided non-uniformly*/
#define BG_ENG_NUM                            MAX_SUBBAND_NUM  /* number of energy of sub-band divided non-uniformly*/
#define MIN_AMP_ID                            5
#define MAX_AMP_ID                            64
#define SPEC_AMP_NUM                          (MAX_AMP_ID-MIN_AMP_ID+1)
#define STABLE_NUM                            4        /* number of time-domain stable rate*/
#define SFM_NUM                               3        /* number of spectral flatness  */


#define START_NG                              5         /* Stationary noise UV modification */
#define FULL_NG                               10        /* Stationary noise UV modification */
#define ISP_SMOOTHING_QUANT_A1                0.9f      /* Stationary noise UV modification */

#define KP559016994                           0.55901699f  /* EDCT & EMDCT constants */
#define KP951056516                           0.95105652f  /* EDCT & EMDCT constants */
#define KP587785252                           0.58778525f  /* EDCT & EMDCT constants */
#define KP866025403                           0.86602540f  /* EDCT & EMDCT constants */
#define KP250000000                           0.25000000f  /* EDCT & EMDCT constants */

#define FEC_BITS_CLS                          2          /* FEC - number of bits for clas information */
#define FEC_BITS_ENR                          5          /* FEC - number of bits for energy information */
#define FEC_ENR_STEP                          (96.0f/(1<<FEC_BITS_ENR))
#define FEC_ENR_QLIMIT                        ((1<<FEC_BITS_ENR)-1)
#define FEC_BITS_POS                          8          /* FEC - number of bits for glottal pulse position */
#define L_SYN_MEM_CLAS_ESTIM                  (2*PIT16k_MAX - L_FRAME16k) /* FEC - memory of the synthesis signal for frame class estimation */
#define L_SYN_CLAS_ESTIM                      (L_SYN_MEM_CLAS_ESTIM + L_FRAME16k)   /* FEC - length of the synthesis signal buffer for frame class estimation */

#define FRAME_COUNT_HF_SYNTH                  100        /* Threshold of frame counter in HF synthesis */

#define AUDIO_COUNTER_INI                     200        /* Counter initialization      */
#define AUDIO_COUNTER_STEP                    10         /* Counter increase on each audio frame      */
#define AUDIO_COUNTER_MAX                     1000       /* Counter saturation      */

#define BWD_TOTAL_WIDTH                       320        /* BWD width */
#define BWD_COUNT_MAX                         100        /* maximum value of BWD counter              */

#define PREEMPH_FAC                           0.68f      /* preemphasis factor at 12.8kHz                */

#define PREEMPH_FAC_16k                       0.72f


#define PREEMPH_FAC_SWB                       0.9f       /* preemphasis factor for super wide band       */


#define GAMMA1                                0.92f      /* weighting factor (numerator) default:0.92    */
#define GAMMA16k                              0.94f

#define FORMANT_SHARPENING_G1                 0.75f      /* Formant sharpening numerator weighting at 12.8kHz   */
#define FORMANT_SHARPENING_G2                 0.9f       /* Formant sharpening denominator weighting at 12.8kHz */
#define FORMANT_SHARPENING_G1_16k             0.8f       /* Formant sharpening numerator weighting at 16kHz     */
#define FORMANT_SHARPENING_G2_16k             0.92f      /* Formant sharpening denominator weighting at 16kHz   */

#define FSCALE_DENOM                          512

#define ACELP_FIXED_CDK_NB                    41
#define ACELP_FIXED_CDK_BITS(n)               PulseConfTable[n].bits

#define L_FIR                                 31

enum TRACKPOS
{
    TRACKPOS_FIXED_FIRST        =0,         /* Fill tracks from left */
    TRACKPOS_FIXED_EVEN         =1,         /* Even tracks */
    TRACKPOS_FIXED_FIRST_TWO    =2,         /* One track of 32 */
    TRACKPOS_FIXED_TWO          =3,         /* Two tracks of 32 instead of four of 16 */
    TRACKPOS_FREE_ONE           =4,         /* One freely moving track with extra pulse */
    TRACKPOS_FREE_THREE         =6,         /* Three freely moving tracks with single extra pulses */
    TRACKPOS_GRADIENT           =7
};

enum
{
    LAGW_WEAK,    /* weak lag window            */
    LAGW_MEDIUM,  /* medium strength lag window */
    LAGW_STRONG,  /* strong lag window          */

    NUM_LAGW_STRENGTHS
};

/*----------------------------------------------------------------------------------*
 * ACELP@16kHz core constants
 *----------------------------------------------------------------------------------*/

#define NB_SUBFR16k                           5         /* number of subframes per frame @16kHz             */
#define INT_FS_16k                            16000.0f  /* CELP core internal sampling frequency @16kHz     */

#define PIT16k_MIN                            42        /* Minimum pitch lag @16kHz                         */
#define PIT16k_MAX                            289       /* Maximum pitch lag @16kHz                         */
#define PIT16k_FR2_TC0_2SUBFR                 83        /* Minimum pitch lag with resolution 1/2 @16kHz for TC02, 2nd subframe */
#define PIT16k_MIN_EXTEND                     21        /* Minimum pitch lag of extended range @16kHz       */
#define PIT16k_FR2_EXTEND_9b                  88        /* Minimum 9 bit pitch lag with resolution 1/2 of extended range @16kHz     */
#define PIT16k_FR1_EXTEND_9b                  130       /* Minimum 9 bit pitch lag with resolution 1 of extended range @16kHz       */
#define PIT16k_FR2_EXTEND_10b                 264       /* Minimum 10 bit pitch lag with resolution 1/2 of extended range @16kHz    */

#define WIDTH_BAND                            8         /* sub-band width in AVQ coding                             */
#define G_AVQ_MIN                             0.80f     /* lower limit for gain Q in higher-rate ACELP contribution */
#define G_AVQ_MAX                             96.0f     /* upper limit for gain Q in higher-rate ACELP contribution */
#define FAC_PRE_AVQ                           0.3f      /* preemhasis factor in ACELP pre-quantizer                 */

#define G_AVQ_MIN_INACT                       0.70f     /* lower limit for gain Q in higher-rate ACELP contribution, inactive segments */
#define G_AVQ_MAX_INACT                       4.1f      /* upper limit for gain Q in higher-rate ACELP contribution, inactive segments */
#define G_AVQ_MIN_INACT_48k                   0.35f     /* lower limit for gain Q in higher-rate ACELP contribution, inactive segments, 48 kbit/s */
#define G_AVQ_MAX_INACT_48k                   2.8f      /* upper limit for gain Q in higher-rate ACELP contribution, inactive segments, 48 kbit/s */
#define G_AVQ_MIN_INACT_64k                   0.25f     /* lower limit for gain Q in higher-rate ACELP contribution, inactive segments, 64 kbit/s */
#define G_AVQ_MAX_INACT_64k                   1.5f      /* upper limit for gain Q in higher-rate ACELP contribution, inactive segments, 64 kbit/s */
#define G_AVQ_DELTA_INACT_48k                 (G_AVQ_MAX_INACT_48k - G_AVQ_MIN_INACT_48k) / ((1 << G_AVQ_BITS) - 1)
#define G_AVQ_DELTA_INACT_64k                 (G_AVQ_MAX_INACT_64k - G_AVQ_MIN_INACT_64k) / ((1 << G_AVQ_BITS) - 1)
#define G_AVQ_BITS                            6         /* number of bits to quantize the AVQ gain in higher-rate ACELP contribtuion   */
#define G_AVQ_DELTA                           (G_AVQ_MAX - G_AVQ_MIN) / ((1 << G_AVQ_BITS) - 1)
#define G_AVQ_DELTA_INACT                     (G_AVQ_MAX_INACT - G_AVQ_MIN_INACT) / ((1 << G_AVQ_BITS) - 1)

#define G_PITCH_MIN                           0.00f     /* SQ of gains: pitch gain lower limit */
#define G_PITCH_MAX                           1.22f     /* SQ of gains: pitch gain upper limit */
#define G_CODE_MIN                            0.02f     /* SQ of gains: code gain lower limit */
#define G_CODE_MAX                            5.00f     /* SQ of gains: code gain upper limit */

#define G_PITCH_MIN_TC192                     0.1f
#define G_PITCH_MAX_TC192                     0.95f
#define G_CODE_MIN_TC192                      0.6f
#define G_CODE_MAX_TC192                      41.0f

/*--------------------------------------------------------------*
 * TCX constants
 *---------------------------------------------------------------*/

#define NOISE_FILL_RANGES                     1
#define NBITS_NOISE_FILL_LEVEL                3   /* Number of bits used for coding noise filling level for each range */
#define MIN_NOISE_FILLING_HOLE                8
#define HOLE_SIZE_FROM_LTP(gain)              (4+(int)(2.0f*gain*(4.0f/0.625f)))
#define FDNS_NPTS                             64
#define AVG_TCX20_LSF_BITS                    40
#define AVG_TCX10_LSF_BITS                    59
#define LTPSIZE                               3
#define TCXLTP_DELAY_NS                       250000
#define TCXLTP_MAX_DELAY                      NS2SA(48000,TCXLTP_DELAY_NS)
#define TCXLTP_LTP_ORDER                      24
#define TCX_RES_Q_BITS_GAIN                   3

/* Use arithmetic coder with LPC shaping also at high (i.e. TCX-only) bitrates */
#define LPC_SHAPED_ARI_MAX_RATE               ACELP_9k60
#define N_MAX_ARI                             800

/*----------------------------------------------------------------------------------*
 * TNS constants
 *----------------------------------------------------------------------------------*/

#define R1_48                                 690
#define R2_48                                 420
#define R1_16                                 230
#define R2_16                                 140
#define R1_25                                 368
#define R2_25                                 224
#define TNS_MAX_NUM_OF_FILTERS                2                                                     /* TNS maximum number of filters                                    */
#define TNS_MAX_FILTER_ORDER                  8                                                     /* TNS maximum filter order                                         */
#define ITF_MAX_FILTER_ORDER                  16                                                    /* ITF maximum filter order                                         */
#define NPRM_TNS                              (2+TNS_MAX_NUM_OF_FILTERS*(3+TNS_MAX_FILTER_ORDER))   /* TNS total number of quantized parameters                         */
#define NPRM_RESQ                             100                                                   /* Maximum number of parameter for residual Q in TCX                */
#define NPRM_CTX_HM                           3                                                     /* Number of Parameters for Context HM : flag+index*/
#define NPRM_DIV                              (2+NPRM_TNS+N_MAX/2+NPRM_RESQ+NPRM_CTX_HM)            /* Total number of quantized parameter in 1 division                */
#define DEC_NPRM_DIV                          NPRM_DIV                                              /* Total number of quantized parameter in 1 division (decoder side) */
#define NPRM_LPC_NEW                          2*(4+2*NB_SPHERE)                                     /* LPC total number of quantized parameters                         */

#define BITBUFSIZE (128000/50)

#define TNS_COEF_RES                          4     /* Bit resolution of the coefficients. */
#define INDEX_SHIFT                           (1 << (TNS_COEF_RES-1)) /* For shifting the index so that index 0 points to 0. */

/*----------------------------------------------------------------------------------*
 * LSF quantization constants
 *----------------------------------------------------------------------------------*/

#define GENERIC_MA_LIMIT                      ACELP_9k60 /* crossover limit, for Generic WB mode, < use SN/AR, >= use MA-predictor */

#define NC16k                                 (M16k/2)
#define NO_ITER                               4           /* number of iterations for tracking the root */
#define SPC                                   0.0234952f
#define SPC_plus                              SPC * 1.001f
#define ALPHA_SQ                              ((0.5f / PI2) * (0.5f / PI2))

#define NC                                    M/2
#define LSF_GAP                               50.0f
#define LSF_BITS_CNG                          29

#define MU_MA                                 (1.0f/3.0f) /* original prediction factor (only for AMR-WB IO mode) */
#define ISF_GAP                               50        /* Minimum ISF separation for end-frame ISFs (only in AMR-WB IO mode)  */
#define LSF_GAP_MID                           80.0f     /* Minimum LSF separation for mid-frame LSFs  */
#define MODE1_LSF_GAP                               70.0f     /* Minimum LSF separation for end-frame ISFs  */
#define PREFERSFNET                           1.05
#define SFNETLOWLIMIT_WB                      35000     /* new sampling rate dependent thresholds used in LSF codebook decision logic, WB case */
#define SFNETLOWLIMIT_NB                      38000     /* new sampling rate dependent thresholds used in LSF codebook decision logic, NB case */
#define LSFMBEST                              2         /* number of survivors from one stage to another */
#define STREAKLEN                             3         /* Allow this many predictive frames, before starting limiting */
#define STREAKMULT                            0.8f      /* Exponential limiting multiplier */

#define LSFMBEST_MAX                          16

#define TCXLPC_NUMSTAGES                      3
#define TCXLPC_NUMBITS                        13
#define TCXLPC_IND_NUMSTAGES                  1
#define TCXLPC_IND_NUMBITS                    2
#define TCXLPC_LSF_GAP                        80.0f

#define MAX_VQ_STAGES                         4
#define MAX_VQ_STAGES_USED                    9         /* this is the maximum number of stages currently used and changing this will affect the memory allocated
                                                          MAX_VQ_STAGES is also used as offset for addressing some arrays, so this should NOT be changed*/
#define MIDLSF_NBITS                          5
#define ENDLSF_NBITS                          31

#define LEN_INDICE                            15
#define LATTICE_DIM                           8
#define NO_LEADERS                            49
#define MAX_NO_BR_LVQ                         28
#define MAX_NO_SCALES                         3
#define MAX_NO_VALS                           4
#define WB_LIMIT_LSF                          6350
#define CNG_LVQ_MODES                         16

#define MAX_NO_MODES                         128
#define START_CNG                            112
#define MAX_NO_MODES_p                       145
#define NO_CODING_MODES                        6
#define LVQ_COD_MODES                         18

/* BC-TCVQ */
#define N_STAGE_VQ                            8
#define N_DIM                                 2
#define NUM_SUBSET                            8
#define OP_LOOP_THR_HVO                       3784536.3f /* 80% : Open-loop Threshold  */
#define NUM_STATE                             16        /* BC-TCQ - Number of state of the Trellis */
#define N_STAGE                               16        /* BC-TCQ - Smaple number in a frame */

#define SIZE_BK1                              256
#define SIZE_BK2                              256
#define SIZE_BK21                             64
#define SIZE_BK22                             128
#define SIZE_BK23                             128
#define SIZE_BK24                             32
#define SIZE_BK25                             32
#define SIZE_BK21_36b                         128
#define SIZE_BK22_36b                         128
#define SIZE_BK23_36b                         64

#define NB_QUA_GAIN5B                         32     /* Number of quantization level        */
#define NB_QUA_GAIN6B                         64     /* Number of quantization level        */
#define NB_QUA_GAIN7B                         128    /* Number of quantization level        */

/*----------------------------------------------------------------------------------*
 * Transient detection
 *----------------------------------------------------------------------------------*/

#define NSUBBLOCKS                            8              /* Number of subblocks per frame, one transient per a sub-block can be found */
#define MAX_TD_DELAY                          2*NSUBBLOCKS   /* Maximum allowed delay (in number of subblocks) of the transient detection, affects required memory */

#define NO_TCX                                0
#define TCX_20                                1
#define TCX_10                                2
#define TCX_5                                 3

#define TRANSITION_OVERLAP                    (-2)
#define RECTANGULAR_OVERLAP                   (-1)
#define FULL_OVERLAP                          0
#define NOT_SUPPORTED                         1
#define MIN_OVERLAP                           2
#define HALF_OVERLAP                          3
#define ALDO_WINDOW                           4

#define SWITCH_OVERLAP_8k                     15        /* == NS2SA(8000, SWITCH_GAP_LENGTH_NS) - NS2SA(8000, 10000000.0f - N_ZERO_MDCT_NS) */
#define SWITCH_GAP_LENGTH_8k                  50

/*----------------------------------------------------------------------------------*
 * FEC constants
 *----------------------------------------------------------------------------------*/

#define UNVOICED_CLAS                         0         /* Unvoiced, silence, noise, voiced offset   */
#define UNVOICED_TRANSITION                   1         /* Transition from unvoiced to voiced components - possible onset, but too small */
#define VOICED_TRANSITION                     2         /* Transition from voiced - still voiced, but with very weak voiced characteristics     */
#define VOICED_CLAS                           3         /* Voiced frame, previous frame was also voiced or ONSET                   */
#define ONSET                                 4         /* Voiced onset sufficiently well built to follow with a voiced concealments    */
#define SIN_ONSET                             5         /* Artificial harmonic+noise onset (used only in decoder)                       */
#define INACTIVE_CLAS                         6         /* Inactive frame (used only in decoder)                  */
#define AUDIO_CLAS                            7         /* Audio frame (used only in AMR-WB IO mode) */

#define BETA_FEC                              0.75f     /* FEC - weighting factor for LSF estimation in FER */
#define STAB_FAC_LIMIT                        0.25f     /* FEC - limit at which safety net is forced for next frame */

#define MODE1_L_FIR_FER                       5         /* FEC - impulse response length for low- and high-pass filters in FEC */
#define L_FIR_FER                             3         /* impulse response length for low- & high-pass filters in FER concealment  */
#define L_FIR_FER2                            11        /* new filter tuning: 11*/
#define MAX_UPD_CNT                           5         /* FEC - maximum number of frames since last pitch update */
#define ALPHA_S                               0.6f      /* FEC - damping factor for SIN_ONSET frames */
#define ALPHA_V                               1.0f      /* FEC - damping factor for VOICED_CLAS frames */
#define ALPHA_VT                              0.4f      /* FEC - damping factor for VOICED_TRANSITION frames */
#define ALPHA_UT                              0.8f      /* FEC - damping factor for UNVOICED_TRANSITION frames */
#define ALPHA_U                               0.4f      /* FEC - damping factor for UNVOICED_CLAS frames */
#define ALPHA_UU                              1.0f      /* FEC - damping factor for UNVOICED_CLAS frames */

#define AGC                                   0.98f

#define PLC_MIN_CNG_LEV                       0.01f  /* minimum background level */
#define PLC_MIN_STAT_BUFF_SIZE                50     /* buffer size for minimum statistics */
#define PLC_MIN_CNG_LEV                       0.01f
#define G_LPC_RECOVERY_BITS                   1

/*----------------------------------------------------------------------------------*
 * Transition mode (TC) constants
 *----------------------------------------------------------------------------------*/

/* Conversion of tc_subfr to index */
#define TC_SUBFR2IDX(x)                     ( x == 0   ? 0 : \
                                              x == 1   ? 0 : \
                                              x == 2   ? 1 : \
                                              x == 3   ? 2 : \
                                              x == 4   ? 3 : \
                                              x == 64  ? 4 : \
                                              x == 128 ? 5 : \
                                              x == 192 ? 6 : \
                                              x == 256 ? 7 : 0 )

#define TC_SUBFR2IDX_16KHZ(x)               ( x == 0   ? 0 : \
                                              x == 64  ? 1 : \
                                              x == 128 ? 2 : \
                                              x == 192 ? 3 : \
                                              x == 256 ? 4 : 0 )

#define L_IMPULSE                             17        /* TC - length of one prototype impulse                                                          */
#define L_IMPULSE2                            8         /* TC - half-length of one prototype impulse == floor(L_IMPULSE/2)                               */
#define NUM_IMPULSE                           8         /* TC - number of prototype impulses                                                             */
#define N_GAIN_CODE_TC                        8         /* TC - number of levels for gain_code quantization for subrames without glot. impulse(s) -      */
#define N_GAIN_TC                             8         /* TC - number of levels for gain_trans quantization                                             */
/* TC - attention: DO NOT CHANGE the following constants - needed for correct bit-allocations        */
#define TC_0_0                                1         /* TC - subframe ID for TC: first glottal impulse in the 1st subframe, second in the 1st subframe    */
#define TC_0_64                               2         /* TC - subframe ID for TC: first glottal impulse in the 1st subframe, second in the 2nd subframe    */
#define TC_0_128                              3         /* TC - subframe ID for TC: first glottal impulse in the 1st subframe, second in the 3rd subframe    */
#define TC_0_192                              4         /* TC - subframe ID for TC: first glottal impulse in the 1st subframe, second in the 4th subframe    */

/*----------------------------------------------------------------------------------*
 * AVQ constants
 *----------------------------------------------------------------------------------*/

#define NB_LDQ3                               9         /* RE8 constants */
#define NB_SPHERE                             32
#define NB_LEADER                             36
#define NB_LDQ4                               27
#define FAC_LOG2                              3.321928095f

#define NSV_MAX                               34        /* maximal number of sub-vectors used by the AVQ */

/*----------------------------------------------------------------------------------*
 * Arithmetic coder
 *----------------------------------------------------------------------------------*/

#define A_THRES_SHIFT                         2
#define A_THRES                               (1<<A_THRES_SHIFT)
#define VAL_ESC                               16
#define NBITS_CONTEXT                         8
#define NBITS_RATEQ                           2

#define cbitsnew                              16
#define stat_bitsnew                          14

#define ari_q4new                            (((long)1<<cbitsnew)-1)
#define ari_q1new                            (ari_q4new/4+1)
#define ari_q2new                            (2*ari_q1new)
#define ari_q3new                            (3*ari_q1new)

#define kLtpHmFractionalResolution            7
#define kLtpHmFlag                            (1 << 8)

#define MAX_LENGTH                            L_FRAME_MAX

/*----------------------------------------------------------------------------------*
 * TCQ constants
 *----------------------------------------------------------------------------------*/

#define MAX_PULSES                            560

#define NORMAL_HQ_CORE                        0         /* Signal use of Normal HQ core */
#define LOW_RATE_HQ_CORE                      1         /* Signal use of Low Rate MDCT core */
#define LOW_RATE_HQ_CORE_TRAN                 2         /* Signal use of Low Rate MDCT core Tran SWB */
#define NORM_MDCT_FACTOR                      L_FRAME8k /* Normalize Low Rate MDCT coefficients to this frame size */
#define BANDS_MAX                             (4*8)
#define MAX_GQLEVS                            32        /* Max fine gain levels */
#define BITS_DE_CMODE                         1
#define BITS_DE_HMODE                         1
#define BITS_DE_8SMODE                        1
#define MAXIMUM_ENERGY_LOWBRATE               255
#define MINIMUM_ENERGY_LOWBRATE               -256
#define BITS_ABS_ENG                          7
#define ABS_ENG_OFFSET                        64
#define BITS_MAX_DEPTH                        3
#define BITS_DE_8SMODE_N0                     1
#define BITS_DE_8SMODE_N1                     1
#define BITS_DE_8SPOS                         5
#define BITS_DE_FCOMP                         5
#define BITS_DE_LSB                           1
#define DE_OFFSET0                            46
#define DE_OFFSET1                            32
#define DE_LIMIT                              64
#define LRMDCT_BE_OFFSET                      15
#define LRMDCT_BE_LIMIT                       31

#define HQCORE_NB_MIN_RATE                    7200     /* NB LR MDCT coding down to this bit rate */
#define HQCORE_WB_MIN_RATE                    13200    /* WB LR MDCT coding down to this bit rate */
#define HQCORE_SWB_MIN_RATE                   13200    /* SWB LR MDCT coding down to this bit rate */

#define LRMDCT_CROSSOVER_POINT                16400    /* Use LR MDCT core at this rate and below */

#define HTH_NORM                              17
#define LTH_NORM                              13
#define OFFSET_NORM                            3


/*----------------------------------------------------------------------------------*
 * SWB TBE constants
 *----------------------------------------------------------------------------------*/

#define STEPSNUM                              4                    /* Number of steps in a2lsp routine for SHB LPC */
#define ALLPASSSECTIONS_STEEP                 3                    /* Size of all pass filters for interpolation and decimation by a factor of 2 */
#define INTERP_3_1_MEM_LEN                    13                   /* Size of all pass filters for interpolation and decimation by a factor of 3:1 */
#define INTERP_3_2_MEM_LEN                    15
#define L_SHB_LAHEAD                          20                   /* Size of lookahead for SHB */
#define NUM_SHB_SUBFR                         16
#define LPC_SHB_ORDER                         10
#define LPC_WHTN_ORDER                        4                    /* Order of whitening filter for SHB excitation */
#define SHB_OVERLAP_LEN                       (L_FRAME16k-L_SHB_LAHEAD)/(NUM_SHB_SUBFR-1)
#define QUANT_DIST_INIT                       (10000000000.0f)     /* Quantiser search distance initialisation */
#define HIBND_ACB_L_FAC                       5/2                  /* SHB Interpolation Factor */
#define NUM_HILBERTS                          2
#define HILBERT_ORDER1                        5
#define HILBERT_ORDER2                        4
#define HILBERT_MEM_SIZE                      (HILBERT_ORDER1 + (2*HILBERT_ORDER2) + (2*HILBERT_ORDER2))

#define L_SHB_TRANSITION_LENGTH               (2*NS2SA(48000, DELAY_BWE_TOTAL_NS))

#define NUM_BITS_SHB_SubGain                  6
#define NUM_BITS_SHB_FrameGain                6

#define NUM_BITS_SHB_FrameGain_LBR_WB         4
#define RECIP_ROOT_EIGHT                      0.3535534f          /* 1.0 / sqrt(8.0) - constant Gain Shape over TD BWE subframes */

#define LPC_SHB_ORDER_WB                      6
#define LPC_WHTN_ORDER_WB                     2                    /* Order of whitening filter for WB excitation */
#define NUM_BITS_WB_LSF                       8
#define LPC_SHB_ORDER_LBR_WB                  4
#define NUM_BITS_LBR_WB_LSF                   2

#define COMP_FIL_ORDER                        19
#define MAX_BIQ_N                             L_FRAME32k

#define NUM_SHB_SUBGAINS                      4                    /* Number of subframe gains */
#define NUM_BITS_SHB_SUBGAINS                 5                    /* Number of bits for subframe gains for SWB */
#define NUM_BITS_SHB_FRAMEGAIN                5                    /* Number of bits for framegain for SWB */
#define NUM_BITS_SHB_ENER_SF                  6
#define NUM_BITS_SHB_RES_GS                   3
#define NUM_BITS_SHB_VF                       3
#define NUM_BITS_SHB_SUBGAINS_RF              5                    /* Number of bits for subframe gains for SWB in RF */
#define SHB_GAIN_QLOW                         -1.0f                /* SHB gain lowest scalar quantizer value */
#define SHB_GAIN_QDELTA                       0.15f                /* SHB gain scalar quantizer stepsize */
#define NUM_Q_LSF                             5                    /* Number of quantized LSFs */
#define MIRROR_POINT_BITS                     2                    /* Number of bits used to quantize mirror point */
#define MIRROR_POINT_Q_CB_SIZE                4                    /* Size of codebook used to quantize mirror point */
#define MAX_LSF                               0.5f                 /* Maximum value of the LSFs */
#define NUM_MAP_LSF                           5                    /* Number of mapped LSFs */
#define NUM_LSF_GRIDS                         4                    /* Number of LSF grids */
#define NUM_LSF_GRID_BITS                     2                    /* Number of bits used for the LSF grids */

#define VF_0th_PARAM                          0.34f
#define VF_1st_PARAM                          0.5f
#define VF_2nd_PARAM                          (VF_1st_PARAM - VF_0th_PARAM)

#define GAMMA0                                0.65f                /* Mean value of gamma1/gamma2 for formant PF   */
#define GAMMA_SHARP                           0.15f                /* Largest sharpening for gamma1/gamma2 (0.83/0.67)*/
#define SWB_NOISE_MIX_FAC                     0.15f                /* Noise mixing adjustment factor for active PF */
#define SWB_TILT_LOW                          1.0f                 /* Lower threshold for PF tilt adaptation */
#define SWB_TILT_HIGH                         2.0f                 /* Higher threshold for PF tilt adaptation */
#define SWB_TILT_DELTA                        (1.0f/(SWB_TILT_HIGH-SWB_TILT_LOW)) /* Inclination between thresholds */
#define GAMMA3_PLUS_WB                        0.65f                /* WB post-filter  */
#define GAMMA3_MINUS_WB                       0.85f                /* WB post-filter  */
#define AGC_FAC_WB                            0.85f                /* WB post-filter - gain adjustment factor */
#define AGC_FAC1_WB                           (1.0f-AGC_FAC_WB)    /* WB post-filter - gain adjustment factor complement */

/*----------------------------------------------------------------------------------*
 * SWB BWE constants
 *----------------------------------------------------------------------------------*/

#define INV_L_SUBFR16k                        0.0125f
#define SWB_L_SUBFR                           160
#define FB_L_SUBFR                            240
#define SWB_FENV                              14
#define FB_GAIN_QLOW                          0.0f
#define FB_GAIN_QDELTA                        0.03125f
#define FB_MAX_GAIN_VAR                       0.5f


#define NUM_BITS_FB_FRAMEGAIN                 4   /* Number of bits for framegain for FB */
#define FB_BAND_BEGIN                         620
#define FB_BAND_END                           800
#define FB_BAND_WIDTH                         180
#define N_CAND                                2
#define N_CB11                                32   /* 5bits */
#define N_CB1ST                               128  /* 7bits */
#define N_CB2ND                               64   /* 6bits */
#define N_CB3RD                               32   /* 5bits */
#define N_CB4TH                               64   /* 6bits */
#define DIM1ST                                3
#define DIM2ND                                4
#define DIM3RD                                3
#define DIM4TH                                4
#define DIM11                                 (DIM1ST+DIM2ND)
#define DIM12                                 (DIM3RD+DIM4TH)
#define N_CAND_TR                             3
#define N_CB_TR1                              128
#define N_CB_TR2                              64
#define DIM_TR1                               2
#define DIM_TR2                               2
#define SWB_FENV_TRANS                        4
#define SWB_TENV                              4
#define NUM_SHARP                             9
#define SHARP_WIDTH                           32

#define HARMONIC                              3
#define NORMAL                                2
#define TRANSIENT                             1
#define NOISE                                 0

/*----------------------------------------------------------------------------------*
 * HR SWB BWE constants
 *----------------------------------------------------------------------------------*/

#define NSV_OVERLAP                           2             /* number of sub-bands overlaping with lower-band (0-8kHz) */   /* note that NSV_MAX >= END_FREQ_BWE_FULL/(8*50) + NSV_OVERLAP ! */
#define N_BANDS_BWE_HR                        4             /* number of frequency bands in non-transient frame */
#define N_BANDS_TRANS_BWE_HR                  2             /* number of frequency bands in transient frame */
#define END_FREQ_BWE                          14400         /* maximum frequency coded by AVQ */
#define END_FREQ_BWE_FULL                     16000         /* maximum frequency coded by HR SWB BWE */
#define END_FREQ_BWE_FULL_FB                  20000         /* maximum frequency coded by HR FB BWE */

#define NBITS_GLOB_GAIN_BWE_HR                5             /* number of bits of the global gain quantizer */
#define MIN_GLOB_GAIN_BWE_HR                  3             /* minimum value of the global gain quantizer */
#define MAX_GLOB_GAIN_BWE_HR                  500           /* maximum value of the global gain quantizer */

#define NBITS_ENVELOPE_BWE_HR1                6             /* number of bits for envelope VQ - first two subbands in non-transient frame */
#define NBITS_ENVELOPE_BWE_HR2                5             /* number of bits for envelope VQ - second two subbands in non-transient frame */
#define NBITS_ENVELOPE_BWE_HR_TR              4             /* number of bits for envelope VQ - two subbands in transient frame */
#define NUM_ENVLOPE_CODE_HR1                  64            /* dimension of envelope VQ - first two subbands in non-transient frame */
#define NUM_ENVLOPE_CODE_HR2                  32            /* dimension of envelope VQ - second two subbands in non-transient frame */
#define NUM_ENVLOPE_CODE_HR_TR                16            /* dimension of envelope VQ - two subbands in transient frame */
#define NUM_ENVLOPE_CODE_HR_TR2               8             /* dimension of envelope VQ - two subbands in transient frame */

#define NUM_NONTRANS_START_FREQ_COEF          (L_FRAME32k/2 - NSV_OVERLAP*WIDTH_BAND)                   /* start frequency coefficient (==7.6kHz) in non-transient frame */
#define NUM_NONTRANS_END_FREQ_COEF            (L_FRAME32k*END_FREQ_BWE/END_FREQ_BWE_FULL)               /* end frequency coefficient (==14.4kHz) in non-transient frame */
#define NUM_TRANS_START_FREQ_COEF             (NUM_NONTRANS_START_FREQ_COEF/NUM_TIME_SWITCHING_BLOCKS)  /* start frequency coefficient (==7.6kHz) in transient frame */
#define NUM_TRANS_END_FREQ_COEF               (NUM_NONTRANS_END_FREQ_COEF/NUM_TIME_SWITCHING_BLOCKS)    /* end frequency coefficient (==14.4kHz) in transient frame */
#define NUM_TRANS_END_FREQ_COEF_EFF           140
#define WIDTH_NONTRANS_FREQ_COEF              ((NUM_NONTRANS_END_FREQ_COEF - NUM_NONTRANS_START_FREQ_COEF)/N_BANDS_BWE_HR)  /* number of coefficients per band in non-transient frame */
#define WIDTH_TRANS_FREQ_COEF                 ((NUM_TRANS_END_FREQ_COEF - NUM_TRANS_START_FREQ_COEF)/N_BANDS_TRANS_BWE_HR)  /* number of coefficients per band in transient frame */

#define NBITS_THRESH_BWE_HR                   400           /* BWE HR number of bits threshold */

#define NBITS_HF_GAIN_BWE_HR                  2             /* number of bits for HF (noncoded) energy estimation */
#define BWE_HR_TRANS_EN_LIMIT1                0.1f          /* HF (noncoded) energy equalization limit 1, transient frames */
#define BWE_HR_TRANS_EN_LIMIT2                0.3f          /* HF (noncoded) energy equalization limit 2, transient frames */
#define BWE_HR_TRANS_EN_LIMIT3                0.5f          /* HF (noncoded) energy equalization limit 3, transient frames */
#define BWE_HR_NONTRANS_EN_LIMIT1             0.5f          /* HF (noncoded) energy equalization limit 1, non-transient frames */
#define BWE_HR_NONTRANS_EN_LIMIT2             1.2f          /* HF (noncoded) energy equalization limit 2, non-transient frames */
#define BWE_HR_NONTRANS_EN_LIMIT3             0.8f          /* HF (noncoded) energy equalization limit 3, non-transient frames */

/*----------------------------------------------------------------------------------*
 * FD CNG
 *----------------------------------------------------------------------------------*/

#define OUTMAX_INV                            0.000030517578125f  /* 1/2^15 */
#define OUTMAX_SQ                             1073741824.f  /* 2^30 */
#define OUTMAX_SQ_INV                         0.00000000093132257461547852f  /* 1/2^30 */
#define DELTA                                 (1e-20f)

#define CLDFB_SCALING                           (1.5f)

#define FFTLEN                                640
#define FFTLEN2                               (FFTLEN/2)
#define CORECLDFBLEN                          20
#define TOTCLDFBLEN                           40
#define FFTCLDFBLEN                           (FFTLEN2+TOTCLDFBLEN-CORECLDFBLEN)
#define PERIODOGLEN                           (FFTLEN2-2)
#define NPART                                 24
#define NPARTCLDFB                            10
#define NPART_SHAPING                         62

#define MSSUBFRLEN                            12
#define MSNUMSUBFR                            6
#define MSBUFLEN                              5
#define MSALPHACORALPHA                       0.7f
#define MSALPHACORMAX                         0.3f
#define MSALPHAMAX                            0.96f
#define MSALPHAHATMIN                         0.05f /* It is used for all bands except the first one to get a stable bass */
#define MSQEQINVMAX                           (1.f/5.f)
#define MSAV                                  2.12f
#define MSBETAMAX                             0.8f
#define MSSNREXP                              (-0.02f/0.064f)

#define NB_LAST_BAND_SCALE                    0.8f
#define SWB_13k2_LAST_BAND_SCALE              0.8f

#define CNG_LOG_SCALING                       512.f /*2^9*/

#define M_MAX                                 32
#define N_GAIN_MIN                            4
#define N_GAIN_MAX                            17

#define CHEAP_NORM_SIZE                       161

#define CNA_MAX_BRATE                         ACELP_13k20

/*----------------------------------------------------------------------------------*
 * Bass post-filter constants
 *----------------------------------------------------------------------------------*/

#define NBPSF_PIT_MAX                         (PIT16k_MAX+1)  /* maximum pitch value for bass post-filter */
#define L_TRACK_HIST                          10

/*----------------------------------------------------------------------------------*
 * NB post-filter constants
 *----------------------------------------------------------------------------------*/

#define THRESCRIT                             0.5f      /* NB post-filter - threshold LT pst switch off */
#define AGC_FAC                               0.9875f   /* NB post-filter - gain adjustment factor */
#define AGC_FAC1                              (1.0f-AGC_FAC) /* NB post-filter - gain adjustment factor complement */
#define LONG_H_ST                             20        /* NB post-filter - impulse response length */
#define POST_G1                               0.75f     /* NB post-filter - denominator weighting factor 12kbps */
#define POST_G2                               0.7f      /* NB post-filter - numerator weighting factor 12kbps */
#define GAMMA1_PST                            0.7f      /* denominator weighting factor     */
#define GAMMA2_PST                            0.55f     /* numerator  weighting factor      */
#define GAMMA3_PLUS                           0.2f      /* NB post-filter - tilt weighting factor when k1>0 */
#define GAMMA3_MINUS                          0.9f      /* NB post-filter - tilt weighting factor when k1<0 */
#define F_UP_PST                              8         /* NB post-filter - resolution for fractionnal delay */
#define LH2_S                                 4         /* NB post-filter - length of INT16 interp. subfilters */
#define LH2_L                                 16        /* NB post-filter - length of long interp. subfilters  */
#define MIN_GPLT                              (1.0f/1.5f) /* NB post-filter - LT gain minimum */
#define LH_UP_S                               (LH2_S/2)
#define LH_UP_L                               (LH2_L/2)
#define LH2_L_P1                              (LH2_L + 1)
#define DECMEM_RES2                           (PIT16k_MAX + 2 + LH_UP_L)
#define SIZ_RES2                              (DECMEM_RES2 + L_SUBFR)
#define SIZ_Y_UP                              ((F_UP_PST-1) * (L_SUBFR+1))
#define SIZ_TAB_HUP_L                         ((F_UP_PST-1) * LH2_L)
#define SIZ_TAB_HUP_S                         ((F_UP_PST-1) * LH2_S)
#define POST_G1_MIN                           0.65f
#define POST_G2_MIN                           0.55f
#define POST_G1_NOIS                          0.15f
#define POST_G2_NOIS                          0.10f
#define BG1                                   (-0.01f)
#define BG2                                   (-0.05f)
#define CG1                                   0.9f
#define CG2                                   1.45f
#define C_LP_NOISE                            (0.1f/4.0f)
#define K_LP_NOISE                            15.0f
#define LP_NOISE_THR                          25.0f

/*----------------------------------------------------------------------------------*
 * Stability estimation
 *----------------------------------------------------------------------------------*/

#define NB_BFI_THR                            2         /* threshold for counter of last bad frames */
#define MAX_LT                                40
#define INV_MAX_LT                            (1.0f/MAX_LT)

#define TH_0_MIN                              2.5f
#define TH_1_MIN                              1.875f
#define TH_2_MIN                              1.5625f
#define TH_3_MIN                              1.3125f

/*----------------------------------------------------------------------------------*
 * Speech/music classifier constants
 *----------------------------------------------------------------------------------*/

#define N_FEATURES                            12       /* number of features */
#define N_MIXTURES                            6        /* number of mixtures */
#define M_LSP_SPMUS                           6        /* number of LSPs used in speech/music classifier */
#define NB_BANDS_SPMUS                        15
#define START_BAND_SPMUS                      2
#define N_OLD_BIN_E                           42       /* == (L_FFT/2-2)/3 */

#define LOWEST_FBIN                           3        /* lowest frequency bin for feature vector preparation */
#define HIGHEST_FBIN                          70       /* highest frequency bin for feature vector preparation */
#define HANG_LEN_INIT                         8        /* number of frames for hang-over (causes delay of decision) */
#define HANG_LEN                              8
#define BUF_LEN                               60
#define L_OVR                                 8

/*----------------------------------------------------------------------------------*
 * LD music post-filter constants
 *----------------------------------------------------------------------------------*/

#define TH_0_MIN2                             1.875f
#define TH_1_MIN2                             1.25f
#define TH_2_MIN2                             0.9375f
#define TH_3_MIN2                             0.625f

#define DCT_L_POST                            640
#define OFFSET2                               192

#define VOIC_BINS_HR                          640
#define BIN_16kdct                            (6400.0f/DCT_L_POST)
#define NB_LIMIT_BAND                         16
#define MBANDS_GN_LD                          20       /* number of bands for gain coding in the postfilter */

/*----------------------------------------------------------------------------------*
 * AC mode (GSC) constants
 *----------------------------------------------------------------------------------*/

#define NOISE_LEVEL_SP0                       8
#define NOISE_LEVEL_SP1a                      9
#define NOISE_LEVEL_SP1                       10
#define NOISE_LEVEL_SP2                       12
#define NOISE_LEVEL_SP3                       14

#define MAX_DYNAMIC                           82
#define MIN_DYNAMIC                           50
#define DYNAMIC_RANGE                         (MAX_DYNAMIC-MIN_DYNAMIC)
#define MAX_GSC_NF_BITS                       3
#define GSC_NF_STEPS                          (1 << MAX_GSC_NF_BITS)

#define CRIT_NOIS_BAND                        23

#define SSF                                   32        /* Sub-subframe length for energy estimation in UC decision */
#define NB_SSF                                (L_FRAME / SSF) /* number of sub-subframes per frame */

#define MBANDS_GN                             16       /* Number of band for gain coding in GSC */
#define BAND1k2                               3

#define MBANDS_LOC                            (MBANDS_GN-1)
#define BIN_SIZE                              25.0f
#define SWNB_SUBFR                            1

#define VAR_COR_LEN                           10

#define CFREQ_BITRATE                         ACELP_11k60

#define LT_UV_THR                             100
#define LT_UV_THRMID                          70

#define PIT_EXC_L_SUBFR                       L_FRAME
#define LOCAL_CT                              VOICED

/*----------------------------------------------------------------------------------*
 * Core switching constants
 *----------------------------------------------------------------------------------*/

#define SWITCH_MAX_GAP                        360  /* 6.25 + 1.25 of filter mem max == NS2SA(48000, SWITCH_GAP_LENGTH_NS+DELAY_CLDFB_NS) */

/*----------------------------------------------------------------------------------*
 * HQ core constants
 *----------------------------------------------------------------------------------*/

#define HQ_NORMAL                             0
#define HQ_TRANSIENT                          1
#define HQ_HARMONIC                           2
#define HQ_HVQ                                3
#define HQ_GEN_SWB                            4
#define HQ_GEN_FB                             5

#define PREECHO_SMOOTH_LEN                    20
#define INV_PREECHO_SMOOTH_LENP1              (1 / (PREECHO_SMOOTH_LEN + 1.0));

#define MAX16B                                32767
#define MIN16B                                (-32768)

#define EPSILON                               0.000000000000001f

#define MAX_SEGMENT_LENGTH                    480
#define NUM_TIME_SWITCHING_BLOCKS             4
#define NUM_MAP_BANDS                         20
#define NUM_MAP_BANDS_HQ_24k4                 17
#define NUM_MAP_BANDS_HQ_32k                  18
#define FREQ_LENGTH                           800

#define STOP_BAND                             800

#define SFM_G1                                16
#define SFM_G1G2                              24
#define SFM_N_NB                              18
#define SFM_N_WB                              26
#define SFM_N_STA_8k                          27
#define SFM_N_STA_10k                         30
#define SFM_N_ENV_STAB                        SFM_N_STA_8k  /* Number of bands for env_stab stability measure */
#define SFM_N_ENV_STAB_WB                     SFM_N_WB      /* Number of bands for env_stab stability measure used in HQPLC decision for WB signals */
#define SFM_N_HARMONIC                        39
#define SFM_N                                 36

#define L_HQ_WB_BWE                           20            /* == band_end_wb[SFM_N_WB-1] - (band_start_wb[SFM_N_WB-1]+12) */
#define N_INTL_GRP_16                         2             /* Number of interleaving band groups at 16kHz samplerate */
#define N_INTL_GRP_32                         2             /* Number of interleaving band groups at 32kHz samplerate */
#define N_INTL_GRP_48                         3             /* Number of interleaving band groups at 48kHz samplerate */
#define SFM_N_SWB                             39
#define SFM_N_HARM                            31
#define SFM_N_HARM_FB                         33
#define NB_SFM                                44
#define NB_SFM_MAX                            58
#define WID_G1                                8
#define WID_G2                                16
#define WID_G3                                24
#define WID_GX                                32
#define NUMC_N                                544
#define HQ_MAX_BAND_LEN                       96            /* Largest bandwidth in HQ mode (band_len_harm[32]) */
#define HVQ_PVQ_BUF_LEN                       (HVQ_PVQ_COEFS*(MAX_PVQ_BANDS-1) + HQ_MAX_BAND_LEN) /* 24*7+96 = 216 */

#define QBIT_MAX2                             9

#define FLAGN_BITS                            1
#define GAIN0_BITS                            5
#define GAINI_BITS                            5

#define FLAGS_BITS                            2
#define FLAGS_BITS_FB                         3
#define NORM0_BITS                            5
#define NORMI_BITS                            5
#define NUMNRMIBITS_SWB_STA_8k                5*(SFM_N_STA_8k-1)
#define NUMNRMIBITS_SWB_STA_10k               5*(SFM_N_STA_10k-1)
#define NUMNRMIBITS_SWB_HARMONIC              185
#define NUMNRMIBITS_SWB                       190
#define NUMNRMIBITS                           215
#define NUMNRMIBITS_WB                        125

#define NOHUFCODE                             0
#define HUFCODE                               1
#define HUFF_THR                              10
#define NOSUPERPOSITION                       40

#define MAXVALUEOFFIRSTGAIN                   2.5f
#define MINVALUEOFFIRSTGAIN                   -2.5f
#define NOOFGAINBITS1                         6

#define AUDIODELAYBITS                        6
#define DELTAOFFIRSTGAIN                      (float)(MAXVALUEOFFIRSTGAIN - MINVALUEOFFIRSTGAIN) / (float)((1 << NOOFGAINBITS1) - 1)

#define MAX_D1M_16k                           ((L_FRAME16k>>1) - NS2SA(16000,SWITCH_GAP_LENGTH_NS) - 16)
#define MAX_D1M_12k8                          ((L_FRAME16k>>1) - NS2SA(16000,SWITCH_GAP_LENGTH_NS) - 20)

#define MAX_P_ATT                             40   /* Maximum number of pulses for gain attenuation factor */
#define NB_G                                  4    /* Number of band groups */
#define MAX_GAIN_BITS                         5    /* Maximum number of gain bits */

#define ENV_ADJ_START                         6    /* Number of consecutive bands for which the attenuation is maximum */
#define ENV_ADJ_INCL                          5    /* Inclination for mapping between attenuation region width and attenuation limit */

#define ENV_SMOOTH_FAC                        0.1f  /* Smoothing factor for envelope stability measure  */
#define L_STAB_TBL                            10    /* Number of elements in stability transition table */
#define M_STAB_TBL                            2.571757f    /* Mid point where the transition table is mirrored */
#define D_STAB_TBL                            0.103138f    /* Stability measure step size in transition table  */
#define NUM_ENV_STAB_PLC_STATES               2         /* Number of states of markov model */

#define ATT_LIM_HANGOVER                      150       /* Number of hangover frames for disabling stability dependent attenuation */
#define DELTA_TH                              5.0f      /* Delta energy threshold for transient detection for envelope stability */
#define ENERGY_TH                             100.0f    /* Energy threshold for transient detection */
#define ENERGY_LT_BETA                        0.93f     /* Smoothing factor for long-term energy measure */

#define START_EXC                             60
#define L_HARMONIC_EXC                        202

#define HQ_GENERIC_OFFSET                     2
#define HQ_GENERIC_END_FREQ                   560
#define HQ_GENERIC_END_FREQ_14P2KHZ           568
#define HQ_GENERIC_END_FREQ_16P0KHZ           640

#define HQ_GENERIC_FOFFSET_24K4               80
#define HQ_GENERIC_FOFFSET_32K                144
#define HQ_GENERIC_SWB_NBITS                  31
#define HQ_GENERIC_SWB_NBITS2                 30
#define HQ_GENERIC_FB_NBITS                   5

#define HQ_GENERIC_ST_FREQ                    224
#define HQ_GENERIC_LOW0                       80
#define HQ_GENERIC_HIGH0                      240
#define HQ_GENERIC_HIGH1                      368
#define HQ_GENERIC_HIGH2                      496
#define HQ_GENERIC_LEN0                       128
#define HQ_GENERIC_NVQIDX                     6

#define HQ_GENERIC_EXC0                       0
#define HQ_GENERIC_EXC1                       1
#define HQ_GENERIC_SP_EXC                     2

#define LF_EMP_FAC                            1.2f

#define DIM_FB                                3
#define HQ_FB_FENV                            SWB_FENV + DIM_FB
#define N_CB_FB                               32

#define HVQ_THRES_BIN_24k                     224
#define HVQ_THRES_SFM_24k                     22
#define HVQ_THRES_BIN_32k                     320
#define HVQ_THRES_SFM_32k                     25
#define HVQ_MIN_PEAKS                         2
#define HVQ_MAX_PEAKS_32k                     23
#define HVQ_MAX_PEAKS_24k                     17
#define HVQ_MAX_PEAKS_24k_CLAS                20    /* Limit for HVQ mode */
#define HVQ_MAX_PEAKS                         HVQ_MAX_PEAKS_32k + 1
#define HVQ_NUM_SFM_24k                       (SFM_N_HARMONIC - 1 - HVQ_THRES_SFM_24k)
#define HVQ_NUM_SFM_32k                       (SFM_N_HARMONIC - 1 - HVQ_THRES_SFM_32k)
#define HVQ_E_PEAK_SMOOTH_FAC                 (0.3f)

#define HVQ_MAX_RATE                          32000

#define NUMNRMIBITS_SWB_HVQ_24k               35
#define NUMNRMIBITS_SWB_HVQ_32k               25

#define MAX_PVQ_BANDS                         8
#define HVQ_MAX_PVQ_WORDS                     ((HVQ_MAX_RATE/50)/16 + MAX_PVQ_BANDS)
#define HVQ_MAX_POS_WORDS                     40
#define HVQ_PVQ_COEFS                         24
#define HVQ_BAND_MIN_PULSES                   2
#define HVQ_BAND_MAX_BITS_24k                 80
#define HVQ_BAND_MAX_BITS_32k                 95
#define HVQ_NEW_BAND_BIT_THR                  30

#define HVQ_NF_GROUPS                         2
#define HVQ_NF_WEIGHT1                        0.9578f       /* HVQ Classifier - Noise floor estimate weight 1 */
#define HVQ_NF_WEIGHT2                        0.6472f       /* HVQ Classifier - Noise floor estimate weight 2 */
#define HVQ_PE_WEIGHT1                        0.42237f      /* HVQ Classifier - Peak envelope estimate weight 1 */
#define HVQ_PE_WEIGHT2                        0.80285f      /* HVQ Classifier - Peak envelope estimate weight 2 */
#define HVQ_THR_POW                           0.88f         /* HVQ Classifier power factor for threshold calc */
#define HVQ_SHARP_THRES                       9             /* HVQ Classifier - Sharpness threshold */

#define HVQ_PA_FAC                            0.7071f       /* HVQ Classifier peak allocation factor */
#define HVQ_PA_PEAKS_SHARP1                   9             /* HVQ Classifier - Maximum number of peaks for band with high sharpness */
#define HVQ_PA_PEAKS_SHARP2                   3             /* HVQ Classifier - Maximum number of peaks for band with medium sharpness */
#define HVQ_PA_PEAKS_SHARP3                   2             /* HVQ Classifier - Maximum number of peaks for band with low sharpness */
#define HVQ_PA_SHARP_THRES2                   16.0f         /* HVQ Classifier - Sharpness threshold for band with medium sharpness */
#define HVQ_PA_SHARP_THRES3                   12.0f         /* HVQ Classifier - Sharpness threshold for band with low sharpness */

#define HVQ_BW                                32            /* HVQ Classifier subband bandwidth */
#define HVQ_NSUB_32k                          10
#define HVQ_NSUB_24k                          7             /* HVQ Classifier number of subbands */

#define HVQ_BWE_NOISE_BANDS                   2             /* Number of BWE noise bands */
#define HVQ_BWE_WEIGHT1                       0.95f
#define HVQ_BWE_WEIGHT2                       0.2f
#define HVQ_NFPE_FACTOR                       6.4f
#define HVQ_LB_NFPE_FACTOR                    3.2f

#define HVQ_VQ_DIM                            5             /* HVQ peak VQ dimension */
#define HVQ_PVQ_GAIN_BITS                     5             /* Number of bits to encode PVQ gains in HVQ */
#define HVQ_NUM_CLASS                         4             /* Number of codebook classes */
#define HVQ_CB_SIZE                           256

#define NUM_PG_HUFFLEN                        9             /* Number of Huffman codewords for peak gains */
#define MAX_PG_HUFFLEN                        12            /* Length of the longest codeword for peak gain Huffman coding */

#define HVQ_CP_HUFF_OFFSET                    3             /* HVQ Code Pos - Delta offset */
#define HVQ_CP_HUFF_MAX                       51            /* HVQ Code Pos - Maximum delta for huffman coding */
#define HVQ_CP_HUFF_MAX_CODE                  13            /* HVQ Code Pos - Size of largest code word */
#define HVQ_CP_HUFF_NUM_LEN                   11            /* HVQ Code Pos - Number of different huffman lengths */
#define HVQ_CP_L2_MAX                         64            /* HVQ Code Pos - Layer 2 maximum size */
#define HVQ_CP_L1_LEN                         5             /* HVQ Code Pos - Layer 1 block size */
#define HVQ_CP_MAP_LEN                        8             /* HVQ Code Pos - Mapping table size */
#define HVQ_CP_MAP_IDX_LEN                    3             /* HVQ Code Pos - Mapping index size */
#define HVQ_CP_DELTA                          0             /* HVQ Code Pos - Use Delta coding */
#define HVQ_CP_SPARSE                         1             /* HVQ Code Pos - Use Sparse coding */

#define MAX_SPLITS                            10            /* Maximum number of PVQ band splits */
#define QUANTAQ3OFFSET                        1
enum QuantaMode { NEAREST = 0, CONS };
#define DS_INDEX_LINEAR_END                   21
#define PYR_OFFSET                            1
#define RCF_INIT_SHIFT                        14
#define THR_ADD_SPLIT                         7             /* Threshold for using additional split */
#define PVQ_MAX_BAND_SIZE                     64            /* Maxiumum supported band size for PVQ search */
#define MIN_BAND_SIZE                         1             /* Minimum supported band size for PVQ search */
#define RC_BITS_RESERVED                      1
#define MAX_PVQ_BITS_PER_COEFFICIENT          80            /* Maximum bits per coefficient allocated per PVQ band. Q3. */
#define MAX_SRT_LEN                           NB_SFM_MAX    /* Maximum length of input for srt_vec_ind() */

/*  index_pvq constants  */
#define KMAX                                512
#define KMAX_NON_DIRECT                     96                /* max K  for non-direct indexing  recursion rows  */
#define ODD_DIV_SIZE                        48                /* ind0=1/1  ind1 =1/3  ...  ind47=1/95 */


/* TCQ */
#define TCQ_MAX_BAND_SIZE                   120      /* Maxiumum supported band size for TCQ+USQ search */
#define STATES                              8
#define MAX_AR_FREQ                         16383
#define AR_BITS                             16
#define STATES_LSB                          4
#define TCQ_LSB_SIZE                        24
#define TCQ_AMP                             10
#define QTCQ                                (0.2f)

#define AR_TOP                              ( ( 1 << AR_BITS ) - 1 )
#define AR_FIRST                            ( AR_TOP / 4 + 1 )
#define AR_HALF                             ( 2 * AR_FIRST )
#define AR_THIRD                            ( 3 * AR_FIRST )

#define MAX_SIZEBUF_PBITSTREAM              1024

/*----------------------------------------------------------------------------------*
 * SWB BWE for LR MDCT core
 *----------------------------------------------------------------------------------*/

#define G1_RANGE                              4
#define G1G2_RANGE                            15
#define GRP_SB                                4      /*Maximum subband groups*/
#define THR1                                  4         /* Bit allocation threshold value */
#define THR2                                  5         /* Bit allocation threshold value */
#define THR3                                  6         /* Bit allocation threshold value */

#define NB_SWB_SUBBANDS                       4         /* maximum number of subbands in normal2 subband coding */
#define SWB_SB_LEN0_12KBPS                    55/* length of subband number X in lowest bit rate operation */
#define SWB_SB_LEN1_12KBPS                    68
#define SWB_SB_LEN2_12KBPS                    84
#define SWB_SB_LEN3_12KBPS                    105
#define SWB_HIGHBAND_12KBPS                   (SWB_SB_LEN0_12KBPS+SWB_SB_LEN1_12KBPS+SWB_SB_LEN2_12KBPS+SWB_SB_LEN3_12KBPS)
#define SWB_LOWBAND_12KBPS                    (HQ_GENERIC_END_FREQ_14P2KHZ - SWB_HIGHBAND_12KBPS)
#define SWB_HIGHBAND_MAX                      SWB_HIGHBAND_12KBPS
#define SWB_LOWBAND_MAX                       SWB_LOWBAND_12KBPS

#define SWB_SB_OFF0_12KBPS                    0         /* subband offsets are based on the subband lengths */
#define SWB_SB_OFF1_12KBPS                    (SWB_SB_OFF0_12KBPS + SWB_SB_LEN0_12KBPS)
#define SWB_SB_OFF2_12KBPS                    (SWB_SB_OFF1_12KBPS + SWB_SB_LEN1_12KBPS)
#define SWB_SB_OFF3_12KBPS                    (SWB_SB_OFF2_12KBPS + SWB_SB_LEN2_12KBPS)
#define SWB_SB_OFF4_12KBPS                    (SWB_SB_OFF3_12KBPS + SWB_SB_LEN3_12KBPS)

/* 16.4 kbps */
#define SWB_SB_LEN0_16KBPS                    59/* length of subband number X in lowest bit rate operation */
#define SWB_SB_LEN1_16KBPS                    74
#define SWB_SB_LEN2_16KBPS                    92
#define SWB_SB_LEN3_16KBPS                    115
#define SWB_HIGHBAND_16KBPS                   (SWB_SB_LEN0_16KBPS+SWB_SB_LEN1_16KBPS+SWB_SB_LEN2_16KBPS+SWB_SB_LEN3_16KBPS)
#define SWB_LOWBAND_16KBPS                    (HQ_GENERIC_END_FREQ_16P0KHZ - SWB_HIGHBAND_16KBPS)

#define SWB_SB_OFF0_16KBPS                    0         /* subband offsets are based on the subband lengths */
#define SWB_SB_OFF1_16KBPS                    (SWB_SB_OFF0_16KBPS + SWB_SB_LEN0_16KBPS)
#define SWB_SB_OFF2_16KBPS                    (SWB_SB_OFF1_16KBPS + SWB_SB_LEN1_16KBPS)
#define SWB_SB_OFF3_16KBPS                    (SWB_SB_OFF2_16KBPS + SWB_SB_LEN2_16KBPS)
#define SWB_SB_OFF4_16KBPS                    (SWB_SB_OFF3_16KBPS + SWB_SB_LEN3_16KBPS)

/* SpectrumSmoothing */
#define L_SB                                  12       /* subband length for SpectrumSmoothing */

/* SpectrumSmoothing for NSS */
#define L_SB_NSS                              8
#define L_SB_NSS_HALF                         (L_SB_NSS/2)
#define NUM_SUBBAND_SMOOTH_MAX                (SWB_HIGHBAND_12KBPS/L_SB_NSS+1)
#define MA_LEN                                7

/* Harmonic mode */
#define NB_SWB_SUBBANDS_HAR_SEARCH_SB         2        /* search number of subbands in harmonic subband coding */
#define NB_SWB_SUBBANDS_HAR                   4        /* maximum number of subbands in harmonic subband coding */
#define N_NBIGGEST_PULSEARCH                  18
#define N_NBIGGEST_SEARCH_LRG_B               32


/* 13.2 kbps */
#define SWB_SB_BW_LEN0_12KBPS_HAR             56  /* Group 1 length for BWE */
#define SWB_SB_BW_LEN1_12KBPS_HAR             100  /* Group 2 Length for BWE */
#define SWB_SB_BW_LEN2_12KBPS_HAR             SWB_SB_BW_LEN1_12KBPS_HAR
#define SWB_SB_BW_LEN3_12KBPS_HAR             SWB_SB_BW_LEN0_12KBPS_HAR

/* 16.4 kbps */
#define SWB_SB_BW_LEN0_16KBPS_HAR             60  /* Group 1 length for BWE */
#define SWB_SB_BW_LEN1_16KBPS_HAR             110  /* Group 2 Length for BWE */
#define SWB_SB_BW_LEN2_16KBPS_HAR             SWB_SB_BW_LEN1_16KBPS_HAR
#define SWB_SB_BW_LEN3_16KBPS_HAR             SWB_SB_BW_LEN0_16KBPS_HAR

#define SWB_SB_OFF0_SUB5_12KBPS_HAR           0     /* subband offsets are based on the subband lengths */
#define SWB_SB_OFF1_SUB5_12KBPS_HAR           (SWB_SB_OFF0_SUB5_12KBPS_HAR + SWB_SB_BW_LEN0_12KBPS_HAR)
#define SWB_SB_OFF2_SUB5_12KBPS_HAR           (SWB_SB_OFF1_SUB5_12KBPS_HAR + SWB_SB_BW_LEN1_12KBPS_HAR)
#define SWB_SB_OFF3_SUB5_12KBPS_HAR           (SWB_SB_OFF2_SUB5_12KBPS_HAR + SWB_SB_BW_LEN2_12KBPS_HAR)

#define SWB_SB_OFF0_SUB5_16KBPS_HAR           0     /* subband offsets are based on the subband lengths */
#define SWB_SB_OFF1_SUB5_16KBPS_HAR           (SWB_SB_OFF0_SUB5_16KBPS_HAR + SWB_SB_BW_LEN0_16KBPS_HAR)
#define SWB_SB_OFF2_SUB5_16KBPS_HAR           (SWB_SB_OFF1_SUB5_16KBPS_HAR + SWB_SB_BW_LEN1_16KBPS_HAR)
#define SWB_SB_OFF3_SUB5_16KBPS_HAR           (SWB_SB_OFF2_SUB5_16KBPS_HAR + SWB_SB_BW_LEN2_16KBPS_HAR)

#define LR_BLK_LEN                            16
#define LR_HLF_PK_BLK_LEN                     8
#define LR_LOWBAND_DIF_PK_LEN                 10
#define SWB_HAR_RAN1                          80
#define SWB_HAR_RAN2                          140
#define SWB_HAR_RAN3                          200
#define SPT_SHORTEN_SBNUM                     4

/* LRMDCT fix precision */
#define SWB_BWE_LR_Qs                         12
#define SWB_BWE_LR_Qbe                        14
#define SWB_BWE_LR_QRk                        16


/*----------------------------------------------------------------------------------*
 * FEC for HQ core
 *----------------------------------------------------------------------------------*/

#define MAX_PGF                               7
#define MAX_ROW                               2

#define MAX_SB_NB                             3

#define NELP_LP_ORDER                         8
#define NUM_NELP_GAINS                        10

#define MINIMUM_RATE_TO_ENCODE_VOICING_FLAG   45000
#define FRAC_BWE_SMOOTH                       2.0f    /*  >= 1 */
#define FRAMECTTOSTART_MDCT                   3

/*----------------------------------------------------------------------------------*
 * Channel aware mode (FEC)
 *----------------------------------------------------------------------------------*/

#define FEC_OFFSET                            3
#define MAX_RF_FEC_OFFSET                     9


/*----------------------------------------------------------------------------------*
 * HQ FEC
 *----------------------------------------------------------------------------------*/

#define POST_HQ_DELAY_NS                      DELAY_BWE_TOTAL_NS                    /* delay of post processing after core HQ coding */
#define PH_ECU_ALDO_OLP2_NS                   (ACELP_LOOK_NS/2)                     /* half length of ALDO WINDOW overlap */
#define PH_ECU_LOOKAHEAD_NS                   (11*ACELP_LOOK_NS/(7*2))              /* Number of nanoseconds look-ahead ahead from the end of the past synthesized frame */
#define PH_ECU_MEM_NS                         ((L_PROT48k/48 - 20)*1000000-PH_ECU_LOOKAHEAD_NS) /* Number of nanoseconds memory for Phase ECU before the old_synthFB_fx pointer */

#define N_LEAD_NB                             70                                    /* (N_LEAD_MDCT*(L_FRAME8k/20)) */
#define N_ZERO_NB                             45                                    /* (N_ZERO_MDCT*(L_FRAME8k/20)) */
#define N_LEAD_O_NB                           90                                    /* (20.f-N_LEAD_MDCT)*(L_FRAME8k/20) */
#define N_ZERO_O_NB                           35                                    /* (10.f-N_ZERO_MDCT)*(L_FRAME8k/20) */
#define N_Z_L_NB                              115                                   /* (N_Z_L_MDCT*(float)(L/20)) = N_ZERO_NB + N_LEAD_NB*/
#define N_Z_L_O_NB                            205                                   /* (N_Z_L_O_MDCT*(float)(L/20)) = N_ZERO_NB + N_LEAD_NB + N_LEAD_O_NB */

#define L_PROT32k                             1024                                  /* HQ phase ECU prototype frame length */
#define MAX_PLOCS                             L_PROT48k/4+1                         /* maximum number of spectral peaks to be searched */
#define QUOT_LPR_LTR                          4
#define LGW_MAX                               9                                     /* maximum number frequency group widths */
#define BETA_MUTE_FAC_INI                     0.5f                                  /* initial noise attenuation factor */
#define L_TRANA32k                            (L_PROT32k/QUOT_LPR_LTR)              /* transient analysis frame length */
#define L_TRANA16k                            (L_TRANA32k/2)
#define L_TRANA8k                             (L_TRANA32k/4)
#define L_PROT_HAMM_LEN2_48k                  NS2SA(48000,6000000L)
#define L_PROT_HAMM_LEN2_32k                  NS2SA(32000,6000000L)
#define L_PROT_HAMM_LEN2_16k                  NS2SA(16000,6000000L)
#define L_PROT48k                             L_PROT32k * 3/2                       /* HQ phase ECU prototype frame length */
#define L_PROT48k_2                           L_PROT48k/2
#define L_TRANA48k                            (L_PROT48k/QUOT_LPR_LTR)              /* transient analysis frame length */
#define PH_ECU_SPEC_SIZE                      L_PROT48k
#define T_SIN_PI_2                            (PH_ECU_SPEC_SIZE/4)
#define HQ_FEC_SIGN_SFM                       16
#define OFF_FRAMES_LIMIT                      30                                    /* HQ phase ECU, burst length for muting to zero */
#define PH_ECU_MUTE_START                     15                                    /* HQ phase ECU, burst length to start steep muting */

#define SCALE_DOWN_3dB                        0.7071f
#define MAX_TILT                              0.f
#define ED_THRES                              1.0f

#define ED_THRES_12P                          0.032209f
#define ED_THRES_50P                          0.159063f
#define ED_THRES_90P                          0.532669
#define MAXDELAY_FEC                          224

#define RANDOM_START                          1
#define HQ_FEC_SIGN_THRES                     6
#define HQ_FEC_SIGN_THRES_TRANS               3
#define HQ_FEC_BAND_SIZE                      4


/*--------------------------------------------------------------*
 * Tonal MDCT PLC
 *---------------------------------------------------------------*/

#define MAX_NUMBER_OF_IDX                     30
#define GROUP_LENGTH                          7
#define MAX_PEAKS_FROM_PITCH                  10
#define LAST_HARMONIC_POS_TO_CHECK            128       /* 128 because we check harmonics only up to 3.2 kHz */
#define ALLOWED_SIDE_LOBE_FLUCTUATION         3.0f      /*  4.8 dB */
#define LEVEL_ABOVE_ENVELOPE                  7.59f     /*  8.8 dB */
#define UNREACHABLE_THRESHOLD                 16.0f     /*   12 dB Increase of LEVEL_ABOVE_ENVELOPE so that the threshold is not reached */
#define SMALL_THRESHOLD                       1.10f     /* 0.41 dB Increase of LEVEL_ABOVE_ENVELOPE for the peak detection at a definitive peak in the estimated spectrum */
#define BIG_THRESHOLD                         1.5f      /* 1.76 dB Increase of LEVEL_ABOVE_ENVELOPE for the peak detection at a probable peak in the estimated spectrum */

#define kSmallerLagsTargetBitsThreshold       150
#define kCtxHmOlRSThr                         2.6f


#define kTcxHmNumGainBits                     2         /* Number of bits for the gain index */
#define kTcxHmParabolaHalfWidth               4         /* Parabola half width */
#define kLtpHmGainThr                         0.46f     /* Use the LTP pitch lag in the harmonic model? */

#define LOWRATE_TCXLPC_MAX_BR                 ACELP_9k60

/*--------------------------------------------------------------*
 * Waveform-adjustment MDCT PLC
 *---------------------------------------------------------------*/

#define DEC_STATE_LEN                         10
#define MAX_POST_LEN                          3
#define TCX_TONALITY_INIT_CNT                 7

#define TCX_NONTONAL                          0
#define TCX_TONAL                             1

/*---------------------------------------------------------------*
 * IGF                                                           *
 *---------------------------------------------------------------*/

#define IGF_MAX_TILES                         5
#define IGF_MAX_GRANULE_LEN                   1200
#define IGF_TRANS_FAK                         2
#define IGF_MAX_SFB                           23
#define IGF_NOF_GRIDS                         3
#define IGF_MAX_SUBFRAMES                     2

#define IGF_MODE_WB                           1
#define IGF_MODE_SWB                          2
#define IGF_MODE_FB                           3

#define IGF_BITRATE_WB_9600         0
#define IGF_BITRATE_RF_WB_13200     1
#define IGF_BITRATE_SWB_9600        2
#define IGF_BITRATE_SWB_13200       3
#define IGF_BITRATE_RF_SWB_13200    4
#define IGF_BITRATE_SWB_16400       5
#define IGF_BITRATE_SWB_24400       6
#define IGF_BITRATE_SWB_32000       7
#define IGF_BITRATE_SWB_48000       8
#define IGF_BITRATE_FB_16400        9
#define IGF_BITRATE_FB_24400        10
#define IGF_BITRATE_FB_32000        11
#define IGF_BITRATE_FB_48000        12
#define IGF_BITRATE_FB_96000        13
#define IGF_BITRATE_FB_128000       14
#define IGF_BITRATE_UNKNOWN         15

#define IGF_WHITENING_OFF                     0
#define IGF_WHITENING_MID                     1
#define IGF_WHITENING_STRONG                  2

#define IGF_GRID_LB_NORM                      0
#define IGF_GRID_LB_TRAN                      1
#define IGF_GRID_LB_SHORT                     2

/* constants for IGFSCFDecoder and IGFSCFEncoder */
#define IGF_CTX_OFFSET                        3                                                           /* offset added to make the context values nonnegative, for array indexing */
#define IGF_CTX_COUNT                         (2 * IGF_CTX_OFFSET + 1)                                    /* number of contexts for the AC statistical model */
#define IGF_MIN_ENC_SEPARATE                  -12                                                          /* minimum residual value coded separately, without escape coding */
#define IGF_MAX_ENC_SEPARATE                  +12                                                          /* maximum residual value coded separately, without escape coding */
#define IGF_SYMBOLS_IN_TABLE                  (1 + (IGF_MAX_ENC_SEPARATE - IGF_MIN_ENC_SEPARATE + 1) + 1) /* alphabet size */

/*----------------------------------------------------------------------------------*
 * SC-VBR
 *----------------------------------------------------------------------------------*/

#define UVG1_CBSIZE                           32        /* NELP unvoiced gain-1 codebook size */
#define UVG2_CBSIZE                           64        /* NELP unvoiced gain-2 codebook size */

/* PPP constants */
#define NUM_ERB_WB                            24        /* Number of ERB bands in wideband */
#define NUM_ERB_NB                            22        /* Number of ERB bands in narrowband */

#define VBR_ADR_MAX_TARGET                    6.15f     /* max target ADR for VBR. This rate is used in the closed loop rate control */
#define PPP_LAG_THRLD                         180       /* max lag allowed for PPP coding */

#define MAXLAG_WI                             (PPP_LAG_THRLD/2 + 12)  /* Maximum lag used in waveform interpolation */
#define MAX_LAG_PIT                           (PPP_LAG_THRLD + 21) /* Max possible pitch lag after adding delta lag */

/*----------------------------------------------------------------------------------*
 * JBM
 *----------------------------------------------------------------------------------*/

#define MAX_JBM_SLOTS                         100 /* every primary copy and partial copy stored in JBM needs one slot */
#define MAX_AU_SIZE                           (128000/50/8) /* max frame size in bytes */

/*----------------------------------------------------------------------------------*
 * TEC/TFA
 *----------------------------------------------------------------------------------*/
#define DELAY_TEMP_ENV_BUFF_TEC               9
#define EXT_DELAY_HI_TEMP_ENV                 2


/*----------------------------------------------------------------------------------*
 * BASOP ROM Tables
 *----------------------------------------------------------------------------------*/

#define LD_INT_TAB_LEN    120
#define INV_TABLE_SIZE    256
#define SQRT_TABLE_SIZE   256


/*----------------------------------------------------------------------------------*
 * Decoder modes
 *----------------------------------------------------------------------------------*/


enum
{
    PRIMARY_2800,
    PRIMARY_7200,
    PRIMARY_8000,
    PRIMARY_9600,
    PRIMARY_13200,
    PRIMARY_16400,
    PRIMARY_24400,
    PRIMARY_32000,
    PRIMARY_48000,
    PRIMARY_64000,
    PRIMARY_96000,
    PRIMARY_128000,
    PRIMARY_SID,
    PRIMARY_FUT1,
    SPEECH_LOST,
    NO_DATA
};

enum
{
    AMRWB_IO_6600,
    AMRWB_IO_8850,
    AMRWB_IO_1265,
    AMRWB_IO_1425,
    AMRWB_IO_1585,
    AMRWB_IO_1825,
    AMRWB_IO_1985,
    AMRWB_IO_2305,
    AMRWB_IO_2385,
    AMRWB_IO_SID/*,
    AMRWB_IO_FUT1,
    AMRWB_IO_FUT2,
    AMRWB_IO_FUT3,
    AMRWB_IO_FUT4,
    SPEECH_LOST,
    NO_DATA */
};

enum
{
    G192,
    MIME
};

#endif /* CNST_H */
