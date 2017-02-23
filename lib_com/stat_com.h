/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/


#ifndef STAT_COM_H
#define STAT_COM_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "options.h"
#include "typedef.h"
#include "cnst.h"


/*----------------------------------------------------------------------------------*
 * Declaration of structures
 *----------------------------------------------------------------------------------*/

/* Forward declaration of Decoder_State */
struct Decoder_State;

/* NB postfilter static variables */
typedef struct
{
    int   on;                     /* On/off flag           */
    short reset;                  /* reset flag           */
    float mem_pf_in[L_SUBFR];     /* Input memory          */
    float mem_stp[L_SUBFR];       /*  1/A(gamma1) memory   */
    float mem_res2[DECMEM_RES2];  /* A(gamma2) residual    */
    float mem_zero[M];            /* null memory to compute i.r. of A(gamma2)/A(gamma1) */
    float gain_prec;              /* for gain adjustment   */
} PFSTAT;

typedef struct
{
    float a[MAXLAG_WI];
    float b[MAXLAG_WI];
    int lag;
    int nH;
    int nH_4kHz;
    float upper_cut_off_freq_of_interest;
    float upper_cut_off_freq;
    float Fs;
} DTFS_STRUCTURE;

typedef struct
{
    short lead_sign_ind;
    unsigned int index, size;
    short dim, k_val;
} PvqEntry;

typedef struct
{
    unsigned char buf[MAX_SIZEBUF_PBITSTREAM];
    signed char curPos;
    unsigned int numByte;
    unsigned int numbits;
    unsigned int maxBytes;
} BITSTREAM, *PBITSTREAM;

typedef struct
{
    PBITSTREAM	bsInst;

    unsigned int	low;
    unsigned int	high;

    unsigned int	value;
    int		bits_to_follow;

    int		num_bits;
    int     max_bits;

} ARCODEC, *PARCODEC;

/*---------------------------------------------------------------*
 * Encoder/Decoder Static RAM                                    *
 *---------------------------------------------------------------*/
typedef struct
{
    /*Coding mode info*/
    short mode_index;

    /*LPC info*/
    short midLpc;             /*Flag for using or not mid LPC for the current frame*/
    short midLpc_enable;      /*Flag enabling or not the mid LPC for the current mode*/

    /*ICB flags*/
    short pre_emphasis;       /*Flag for triggering pre-emphasis*/
    short pitch_sharpening;   /*Flag for triggering pitch sharpening*/
    short phase_scrambling;   /*Flag for triggering phase scrambling*/
    short formant_enh;        /*Flag for triggering formant enhancement*/
    short formant_tilt;       /*Flag for triggering formant tile*/

    short voice_tilt;         /*Flag for triggering new voice factor tilt*/

    float formant_enh_num;
    float formant_enh_den;

    short bpf_mode;

    short nrg_mode;
    short nrg_bits;

    short ltp_mode;
    short ltp_bits;

    short ltf_mode;
    short ltf_bits;

    short gains_mode[NB_SUBFR16k];

    int fixed_cdk_index[NB_SUBFR16k];


} ACELP_config;

/* tns_base.h */
/** TNS configuration.
  * Use InitTnsConfiguration to initialize it.
  */
typedef struct
{
    unsigned char maxOrder;

    /** Maximum number of filters. */
    unsigned char nMaxFilters;

    /** Parameters for each TNS filter */
    struct TnsParameters const * pTnsParameters;

    /** Lower borders for each filter.
      * Upper border for the first filter is nsbBorders-1.
      * Upper borders for other filters is the lower border of previous filter.
      */
    short iFilterBorders[TNS_MAX_NUM_OF_FILTERS+1];

} STnsConfig;

/** TNS filter.
  * Parameters that define a TNS filter.
  */
typedef struct
{
    /** Number of subbands covered by the filter. */
    int spectrumLength;
    /** Filter order. */
    int order;
    /** Quantized filter coefficients. */
    int coefIndex[TNS_MAX_FILTER_ORDER];
    /** Prediction gain. The ratio of a signal and TNS residual energy. */
    float predictionGain;
    /** Average squared filter coefficient. */
    float avgSqrCoef;
} STnsFilter;

/** TNS data.
  * TNS data describing all active filters.
  */
typedef struct
{
    /** Number of active filters. */
    int nFilters;
    /** Active filters. */
    STnsFilter filter[TNS_MAX_NUM_OF_FILTERS];
} STnsData;

typedef enum
{
    TNS_NO_ERROR = 0,
    TNS_FATAL_ERROR
} TNS_ERROR;

typedef struct
{
    /* TCX mdct window */
    float tcx_mdct_window[L_MDCT_OVLP_MAX];   /* Sine window for OL decision and DTX transitions*/
    float tcx_aldo_window_1[L_FRAME_MAX];     /* ALDO window long slope */
    float tcx_aldo_window_2[L_MDCT_OVLP_MAX]; /* ALDO window short slope */
    float *tcx_aldo_window_1_trunc;           /* ALDO window truncated long slope */

    short last_aldo;
    float tcx_mdct_window_half[L_MDCT_HALF_OVLP_MAX];
    float tcx_mdct_window_minimum[L_MDCT_MIN_OVLP_MAX];
    float tcx_mdct_window_trans[L_MDCT_TRANS_OVLP_MAX]; /* transition window for ACELP->TCX */

    int tcx5Size;  /* Size of the TCX5 spectrum. Always 5ms. */

    short tcx_curr_overlap_mode; /* window overlap of current frame (0: full, 2: none, or 3: half) */
    short tcx_last_overlap_mode; /* previous window overlap, i.e. left-side overlap of current frame */

    int tcx_mdct_window_length; /* length of overlap-add region*/
    int tcx_mdct_window_half_length; /* length of the "half" overlap */
    int tcx_mdct_window_min_length; /* length of the "minimum" overlap */
    int tcx_mdct_window_trans_length; /* length of the ACELP->TCX overlap */

    int tcx_mdct_window_delay; /*length of window delay*/
    int tcx_offset; /* Offset of the folding point relative to the end of the previous synthetised frame */
    int tcx_mdct_window_length_old; /*for keeping old value for sample rate switching */
    float tcx_mdct_windowFB[L_MDCT_OVLP_MAX];
    float tcx_aldo_window_1_FB[L_FRAME_MAX];      /* ALDO window long slope */
    float tcx_aldo_window_2_FB[L_MDCT_OVLP_MAX];  /* ALDO window short slope */
    float *tcx_aldo_window_1_FB_trunc;            /* ALDO window truncated long slope */

    float tcx_mdct_window_halfFB[L_MDCT_HALF_OVLP_MAX];
    float tcx_mdct_window_minimumFB[L_MDCT_MIN_OVLP_MAX];
    float tcx_mdct_window_transFB[L_MDCT_TRANS_OVLP_MAX]; /* transition window for ACELP->TCX */

    int tcx5SizeFB;  /* Size of the TCX5 spectrum. Always 5ms. */

    int tcx_mdct_window_lengthFB; /*length of window, length of overlap-add region*/
    int tcx_mdct_window_half_lengthFB; /*length of the "half" overlap window*/
    int tcx_mdct_window_min_lengthFB; /*length of the "minimum" overlap window*/
    int tcx_mdct_window_trans_lengthFB; /* length of the ACELP->TCX overlap */

    int tcx_mdct_window_delayFB; /*length of window delay*/
    int tcx_offsetFB;

    int tcx_coded_lines;  /* max number of coded lines, depending on bandwidth mode */

    /* FAC window */
    int lfacNext;
    int lfacNextFB;

    /* TNS */
    int fIsTNSAllowed;
    STnsConfig tnsConfig[2][2];
    STnsConfig const * pCurrentTnsConfig;

    /*Quantization*/
    float sq_rounding; /*set the sq deadzone (no deadzone=0.5f)*/
    short tcxRateLoopOpt;

    /*Bandwidth*/
    float preemph_fac;     /*preemphasis factor*/
    float bandwidth;

    /* Context HM - Residual Quantization*/
    short ctx_hm;    /*Flag for enabling Context HM*/
    short resq;      /*Flag for enabling Residual Quantization*/
    int coder_type;  /*GC,VC,UC*/

    float na_scale;

    float SFM2;

} TCX_config;

/* prot.h */
typedef struct
{
    int bits;            /* bits per subframe */
    int nbiter;          /* number of iterations */
    float alp;           /* initial energy of all fixed pulses */
    int nb_pulse;        /* number of pulses */
    int fixedpulses;     /* number of pulses whose position is determined from correlation and not by iteration */
    int nbpos[13];       /* number of positions tried in the pair-wise search */
    enum TRACKPOS codetrackpos;         /* ordering of tracks -mode */
} PulseConfig;

/* ari.h */
typedef struct
{
    int	low,high,vobf;
} Tastat;

/* FD-based CNG setup */
typedef struct
{
    int fftlen;                  /* FFT length */
    int stopFFTbin;              /* Number of FFT bins to be actually processed */
    int numPartitions;           /* Number of partitions */
    const int* sidPartitions;    /* Upper boundaries for grouping the (sub)bands into partitions when transmitting SID frames (define as NULL pointer to avoid grouping) */

    int numShapingPartitions;    /* Number of partitions */
    const int* shapingPartitions;/* Upper boundaries for grouping the (sub)bands into partitions for shaping at the decoder (define as NULL pointer to avoid grouping) */
} FD_CNG_SETUP;


/* Scale setup */
typedef struct
{
    int bwmode;

    int bitrateFrom;
    int bitrateTo;

    float scale;

} SCALE_SETUP;


/* Arrays and variables common to encoder and decoder */
typedef struct
{
    FD_CNG_SETUP FdCngSetup;

    int     numSlots;       /* Number of time slots in CLDFB matrix */
    int     regularStopBand;/* Number of CLDFB bands to be considered */

    int     numCoreBands;   /* Number of core bands to be decomposed into FFT subbands */
    int     stopBand;       /* Total number of (sub)bands to be considered */
    int     startBand;      /* First (sub)band to be considered */
    int     stopFFTbin;     /* Total number of FFT subbands */
    int     frameSize;      /* Frame size in samples */
    int     fftlen;         /* FFT length used for the decomposition */

    float  timeDomainBuffer[L_FRAME16k];
    float  fftBuffer[FFTLEN];
    float  olapBufferAna[FFTLEN];
    float  olapBufferSynth[FFTLEN];
    float  olapBufferSynth2[FFTLEN];
    const float * olapWinAna;
    const float * olapWinSyn;
    const float * fftSineTab;

    float   msM_win;
    float   msM_subwin;

    short   msFrCnt_init_counter;          /* Frame counter at initialization */
    short   msFrCnt_init_thresh;
    float   init_old;
    int     msFrCnt;        /* Frame counter */
    float   msAlphaCor[2];  /* Correction factor (smoothed) */
    float   msSlope[2];
    float   msQeqInvAv[2];
    int     msMinBufferPtr;

    float   msPsdSum[2];
    float   msPeriodogSum[2];

    short   offsetflag;

    float   periodog[PERIODOGLEN];     /* Periodogram */
    float   cngNoiseLevel[FFTCLDFBLEN];  /* Noise level applied for the CNG in each (sub)band */
    short   seed;           /* Seed memory (for random function) */

    int     npart;          /* Number of partitions */
    int     midband[NPART];        /* Central band of each partition */
    int     nFFTpart;       /* Number of hybrid spectral partitions */
    int     part[NPART];           /* Partition upper boundaries (band indices starting from 0) */
    float   psize[NPART];          /* Partition sizes */
    float   psize_inv[NPART];      /* Inverse of partition sizes */
    float   FFTscalingFactor;      /* Squared ratio between core signal analysis FFT and noise estimator FFT */
    float   scalingFactor;
    int     nCLDFBpart;              /* Number of CLDFB spectral partitions */
    int     CLDFBpart[NPARTCLDFB];     /* CLDFB Partition upper boundaries (band indices starting from 0 above the core coder bands) */
    float   CLDFBpsize_inv[NPARTCLDFB];/* Inverse of CLDFB partition sizes */

    short   inactive_frame_counter;
    short   sid_frame_counter;
    short   active_frame_counter;

    float   sidNoiseEst[NPART];    /* Transmitted noise level */

    int     frame_type_previous;

    float   A_cng[17];
    float   exc_cng[L_FRAME16k];

    int     CngBitrate;
    short   CngBandwidth;

    short   flag_noisy_speech;
    float   likelihood_noisy_speech;

}
FD_CNG_COM;
typedef FD_CNG_COM *HANDLE_FD_CNG_COM;

/* parameter_bitmapping.h */
/** Function that gets specific value from p.
  * @param p Pointer to a variable that can also be structure or array.
  * @param index Index of a variable when p is an array, otherwise 0.
  * @param pValue Pointer to the value.
  * @return Substructure associated with this value or NULL if there is none.
  */
typedef void const * (* TGetParamValue)(void const * p, int index, int * pValue);

/** Function that puts specific value to p.
  * @param p Pointer to a variable that can also be structure or array.
  * @param index Index of a variable when p is an array, otherwise 0.
  * @param value The value.
  * @return Substructure associated with this value or NULL if there is none.
  */
typedef void * (* TSetParamValue)(void * p, int index, int value);

/** Function that return required number of bits for a value when it is coded.
  * @param value The value.
  * @param index Index of a variable when it is an element of an array, otherwise 0.
  * @return Number of bits required to code the value.
  */
typedef int (* TGetNumberOfBits)(int value, int index);

/** Function that encodes a value.
  * @param value The value.
  * @param index Index of a variable when it is an element of an array, otherwise 0.
  * @return Coded value.
  */
typedef int (* TEncodeValue)(int value, int index);

/** Function that decodes a value.
  * @param bitstream Bitstream containing a coded value.
  * @param index Index of a variable when it is an element of an array, otherwise 0.
  * @param pValue A pointer where the decoded value should be stored.
  * @return Number of bits read from the bitstream.
  */
typedef int (* TDecodeValue)(struct Decoder_State *st, int index, int * pValue);

/** Structure that defines mapping between a parameter and a bistream. */
typedef struct ParamBitMap
{
    /** Number of bits in a bitstream required for the parameter.
      * If nBits is equal to 0 then GetNumberOfBits is used.
      */
    int nBits;
    /** Function to get the number of bits required for a value of this parameter.
      * If nBits != 0 it is not used and can be set to NULL.
      * If fZeroAllowed == 0 then GetNumberOfBits must take care of this.
      */
    TGetNumberOfBits GetNumberOfBits;
    /** If fZeroAllowed is 0 then the value can be zero.
      * If the value can't be zero then value-1 is stored in a bitstream.
      * If EncodeValue is not equal to NULL, then the encode/decode function
      * must take care of this flag - there is no additional processing in parameter bitmapping.
      * If EncodeValue is equal to NULL, then the encode/decode function takes care of this.
      */
    int fZeroAllowed;
    /** Function to get the value of this parameter.
      * The function returns a pointer to be used in functions in pSubParamBitMap.
      * If the function returns NULL then the same pointer as for the current
      * parameter is to be used.
      * The function should not do any additional processing if fZeroAllowed == 0,
      * but just set the value as it is.
      */
    TGetParamValue GetParamValue;
    /** Function to set the value of this parameter.
      * The function returns a pointer to be used in functions in pSubParamBitMap.
      * If the function returns NULL then the same pointer as for the current
      * parameter is to be used.
      * The function should not do any additional processing if fZeroAllowed == 0,
      * but just set the value as it is.
      */
    TSetParamValue SetParamValue;

    /** Function to encode a value of this parameter.
      * When it is equal to NULL, fixed-width coding is used.
      * If fZeroAllowed == 0 then EncodeValue must take care of this.
      */
    TEncodeValue EncodeValue;

    /** Function to decode a value of this parameter.
      * When it is equal to NULL, fixed-width coding is used.
      * If fZeroAllowed == 0 then DecodeValue must take care of this.
      */
    TDecodeValue DecodeValue;

    /** Pointer to the map for substructure.
      * The number of structures is determined by this parameter's value.
      * NULL means that there is no substructure.
      */
    struct ParamsBitMap const * pSubParamBitMap;
} ParamBitMap;

/** Structure that defines mapping between parameters and a bistream. */
typedef struct ParamsBitMap
{
    /** Number of parameters in params. */
    int nParams;
    /** Definition of the mapping for each parameter. */
    ParamBitMap params[10];
} ParamsBitMap;

/* tns_tables.h */
struct TnsParameters
{
    /* Parameters for each TNS filter */
    int startLineFrequency;   /* Starting lower frequency of the TNS filter [20..16000] */
    int nSubdivisions;        /* Number of spectrum subdivisions in which the filter operates [1..8) */
    float minPredictionGain;  /* Minimum prediction gain required to turn on the TNS filter */
    float minAvgSqrCoef;      /* Minimum average square of coefficients required to turn on the TNS filter */
};

/**********************************************/
/* Helper structures for hufmann table coding */
/**********************************************/

typedef struct
{
    unsigned char value;
    unsigned short int code;
    unsigned char nBits;
} Coding;


/* Scale TCX setup */
typedef struct
{
    int bwmode;

    int bitrateFrom;
    int bitrateTo;

    float scale;

} SCALE_TCX_SETUP;


/* glob_con.h */
typedef struct
{
    unsigned short frame_bits;          /*Bits per frame*/
    unsigned short frame_net_bits;      /*Net Bits per frame*/
    unsigned char  transmission_bits;    /*max=1*/
    unsigned char  transmission_mode[2]; /*SID,VBR/CBR*/
    unsigned char  bandwidth_bits;      /*max=2*/
    unsigned char  bandwidth_min;       /*first valid bandwidth (NB,WB,SWB,FB)*/
    unsigned char  bandwidth_max;       /*last valid bandwidth (NB,WB,SWB,FB)*/
    unsigned char  reserved_bits;       /*max=1*/
} FrameSizeParams;

typedef struct
{
    int no_channels;
    int no_col;
    int p_filter_length;
    CLDFB_TYPE type;

    const float *p_filter;

    /* rotation vectors */
    const float *rot_vec_ana_re;
    const float *rot_vec_ana_im;
    const float *rot_vec_syn_re;
    const float *rot_vec_syn_im;

    /* memory helper states */
    float *memory;
    unsigned int memory_length;

    /* main filter state */
    float *cldfb_state;

    /* other parameters */
    int bandsToZero;                  /* bands not synthesized */
    int nab;                          /* number of active bands */
    float scale;                      /* scaling of frequency domain */
}
CLDFB_FILTER_BANK;

typedef CLDFB_FILTER_BANK *HANDLE_CLDFB_FILTER_BANK;

typedef struct
{
    float pGainTemp[CLDFB_NO_COL_MAX];
    float loBuffer[CLDFB_NO_COL_MAX + MAX_TEC_SMOOTHING_DEG];

} TEMPORAL_ENVELOPE_CODING_DECODER;
typedef TEMPORAL_ENVELOPE_CODING_DECODER* HANDLE_TEC_DEC;

typedef struct
{
    float loBuffer[CLDFB_NO_COL_MAX + MAX_TEC_SMOOTHING_DEG + DELAY_TEMP_ENV_BUFF_TEC];
    float loTempEnv[CLDFB_NO_COL_MAX];
    float loTempEnv_ns[CLDFB_NO_COL_MAX];
    float hiTempEnv[CLDFB_NO_COL_MAX + DELAY_TEMP_ENV_BUFF_TEC + EXT_DELAY_HI_TEMP_ENV];
    int tranFlag;
    int corrFlag;
} TEMPORAL_ENVELOPE_CODING_ENCODER;
typedef TEMPORAL_ENVELOPE_CODING_ENCODER* HANDLE_TEC_ENC;


typedef enum
{
    FRAME_0 = 0,
    FRAME_2 =  40,
    FRAME_2_4 = 48,
    FRAME_4 = 80,
    FRAME_5_6 = 112,
    FRAME_7_2 = 144,
    FRAME_8 = 160,
    FRAME_9_6 = 192,
    FRAME_13_2 = 264,
    FRAME_16_4 = 328,
    FRAME_24_4 = 488,
    FRAME_32 = 640,
    FRAME_48 = 960,
    FRAME_64 = 1280,
    FRAME_96 = 1920,
    FRAME_128 = 2560
} FRAME_SIZE;

/*---------------------------------------------------------------*
 * IGF                                                           *
 *---------------------------------------------------------------*/

typedef struct igf_grid_struct
{
    int                   swb_offset[IGF_MAX_SFB];
    int                   swb_offset_len;
    int                   startFrequency;
    int                   stopFrequency;
    int                   startLine;
    int                   stopLine;
    int                   startSfb;
    int                   stopSfb;
    int                   sfbWrap[IGF_MAX_TILES+1];
    int                   sbWrap[IGF_MAX_TILES];
    int                   nTiles;
    int                   minSrcSubband;
    int                   minSrcFrequency;
    int                   tile[IGF_MAX_TILES];
    int                   infoIsRefined;
    int                   infoGranuleLen;
    float                 infoTransFac;
    float                 whiteningThreshold[2][IGF_MAX_TILES];
    float                 gFactor;
    float                 fFactor;
    float                 lFactor;
} IGF_GRID, *H_IGF_GRID;

typedef struct IGF_INFO_struct
{
    short                 nfSeed;
    int                   sampleRate;
    int                   frameLength;
    int                   maxHopsize;
    IGF_GRID              grid[IGF_NOF_GRIDS];
    short                 bitRateIndex;
} IGF_INFO, *H_IGF_INFO;


typedef struct
{
    int                 *indexBuffer;
    int                 *peakIndices, *holeIndices;
    int                 numPeakIndices, numHoleIndices;
} CONTEXT_HM_CONFIG;

#endif
