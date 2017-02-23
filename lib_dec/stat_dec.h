/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#ifndef STAT_DEC_H
#define STAT_DEC_H

#include "options.h"
#include "stat_com.h"          /* Common structures                      */
#include "cnst.h"

typedef enum _DEC_MODE
{
    DEC_NO_FRAM_LOSS         = 0x0,
    DEC_CONCEALMENT_EXT      = 0x1
} DEC_MODE;

typedef enum
{
    FRAMEMODE_NORMAL        = 0x0,   /**< frame available            */
    FRAMEMODE_MISSING       = 0x1,   /**< frame missing => conceal   */
    FRAMEMODE_FUTURE        = 0x2
} frameMode;


/* Arrays and variables specific to decoder */
typedef struct
{
    HANDLE_FD_CNG_COM hFdCngCom;

    float   msPeriodog[NPART_SHAPING]; /* Periodogram */
    float   msBminWin[NPART_SHAPING];
    float   msBminSubWin[NPART_SHAPING];
    float   msPsd[NPART_SHAPING];          /* Power Spectral Density estimate (i.e., smoothed periodogram) */
    float   msAlpha[NPART_SHAPING];        /* Optimal smoothing parameter */
    float   msMinBuf[MSNUMSUBFR*NPART_SHAPING];       /* Buffer of minima */
    float   msCurrentMinOut[NPART_SHAPING];
    float   msCurrentMin[NPART_SHAPING];
    float   msCurrentMinSubWindow[NPART_SHAPING];
    int     msLocalMinFlag[NPART_SHAPING];
    int     msNewMinFlag[NPART_SHAPING];
    float   msPsdFirstMoment[NPART_SHAPING];
    float   msPsdSecondMoment[NPART_SHAPING];
    float   msNoiseFloor[NPART_SHAPING];   /* Estimated noise floor */
    float   msNoiseEst[NPART_SHAPING];     /* Estimated noise level */
    float   msLogPeriodog[NPART_SHAPING]; /* Periodogram */
    float   msLogNoiseEst[NPART_SHAPING]; /* Estimated noise level */
    int     npart_shaping;          /* Number of partitions */
    int     nFFTpart_shaping;       /* Number of hybrid spectral partitions */
    int     part_shaping[NPART_SHAPING];           /* Partition upper boundaries (band indices starting from 0) */
    int     midband_shaping[NPART_SHAPING];        /* Central band of each partition */
    float   psize_shaping[NPART_SHAPING];          /* Partition sizes */
    float   psize_inv_shaping[NPART_SHAPING];      /* Inverse of partition sizes */
    float   bandNoiseShape[FFTLEN2];/* CNG spectral shape computed at the decoder */
    float   partNoiseShape[NPART];/* CNG spectral shape computed at the decoder */

    short   flag_dtx_mode;
    float   lp_speech;
    float   lp_noise;

    float   msPeriodogBuf[MSBUFLEN*NPART_SHAPING];
    int     msPeriodogBufPtr;

}
FD_CNG_DEC;
typedef FD_CNG_DEC *HANDLE_FD_CNG_DEC;

typedef struct
{
    int FrameSize;

    int Pitch;
    int T_bfi;

    int Transient[MAX_POST_LEN];
    int TCX_Tonality[DEC_STATE_LEN];

    float outx_new_n1;
    float nsapp_gain;
    float nsapp_gain_n;
    float data_reci2[L_FRAME_MAX];
    float data_noise[L_FRAME_MAX];
    float ener_mean;
    float ener;
    int zp;
    float recovery_gain;
    float step_concealgain;

    int concealment_method;
    int subframe;
    int nbLostCmpt;

    short seed;

} T_PLCInfo;


/*---------------------------------------------------------------*
 * Structures for Tonal MDCT PLC                                 *
 *---------------------------------------------------------------*/

typedef enum
{
    TONALMDCTCONCEAL_OK = 0,

    __error_codes_start = -100,

    TONALMDCTCONCEAL_NSAMPLES_LARGER_THAN_MAXBLOCKSIZE,
    TONALMDCTCONCEAL_INVALIDPOINTER,
    TONALMDCTCONCEAL_UNEXPECTED_ERROR,

    __error_codes_end
} TONALMDCTCONCEAL_ERROR;

typedef struct
{
    unsigned int nSamples;
    unsigned int nSamplesCore;
    Float32 * spectralData;
    float * scaleFactors;
    int blockIsValid;
    int blockIsConcealed;
    int tonalConcealmentActive;
} blockData;

typedef struct
{
    unsigned int numIndexes;
    unsigned short int indexOfTonalPeak[MAX_NUMBER_OF_IDX];
    unsigned short int lowerIndex[MAX_NUMBER_OF_IDX];
    unsigned short int upperIndex[MAX_NUMBER_OF_IDX];
    Float32 phaseDiff[MAX_NUMBER_OF_IDX]; /* This one can be stored with 16 bits in range 0..2*PI */
    Float32 phase_currentFramePredicted[MAX_NUMBER_OF_IDX*GROUP_LENGTH]; /* This one can be stored with 16 bits in range 0..2*PI, but the code has to be adapted to use moduo(2*PI) after adding */
} TonalComponentsInfo;

typedef void (*ApplyScaleFactorsPointer)(float x[], int lg, float const scaleFactors[]);

struct tonalmdctconceal
{
    TCX_config * tcx_cfg;
    void * pMDSTData;
    unsigned int nSamples;
    unsigned int nSamplesCore;
    unsigned int nNonZeroSamples;
    unsigned int nScaleFactors;

    float lastPitchLag;

    blockData lastBlockData;
    blockData secondLastBlockData;

    Float32 scaleFactorsBuffers[2][FDNS_NPTS]; /* Contains also global gain. If it can not be stored in 16 bits with global gain included, then store global gain separately. */
    Float32 spectralDataBuffers[2][L_FRAME_MAX]; /* 16 bits is enough, because it is stored before applying scale factors. Take care that power spectrum is also stored here. */
    Float32 timeDataBuffer[(3*L_FRAME_MAX)/2];
    Float32 * lastPcmOut;
    Float32 * secondLastPcmOut;
    float * secondLastPowerSpectrum;

    float nFramesLost;

    TonalComponentsInfo * pTCI;
};

typedef struct tonalmdctconceal* TonalMDCTConcealPtr;

typedef enum SIGNAL_CLASSIFER_MODE
{
    CLASSIFIER_ACELP,
    CLASSIFIER_TCX
} SIGNAL_CLASSIFIER_MODE;

/*---------------------------------------------------------------*
 * Structures for IGF decoder                                    *
 *---------------------------------------------------------------*/

/* IGFSCFDecoder.h */
typedef struct
{
    int               bitsRead;  /* after a call bitsRead contains the number of bits consumed by the decoder */
    int               prev[64];  /* no more than 64 SCFs for the IGF energy envelope of one block */
    int               scfCountLongBlock;
    int               t;
    int               bitrate;
    const unsigned short   *cf_se00;
    const unsigned short   *cf_se01;
    short                   cf_off_se01;
    const unsigned short   *cf_se02;
    const short            *cf_off_se02;
    const unsigned short   *cf_se10;
    short                   cf_off_se10;
    const unsigned short   *cf_se11;
    const short            *cf_off_se11;
    Tastat            acState;
} IGFSCFDEC_INSTANCE, *IGFSCFDEC_INSTANCE_HANDLE;

/* IGFDec.h */
typedef struct igfdec_private_data_struct
{

    IGF_INFO                  igfInfo;
    /* envelope reconstruction: */
    float                     igf_sN[IGF_MAX_SFB];                     /* only with short blocks as static needed */
    float                     igf_pN[IGF_MAX_SFB];                     /* only with short blocks as static needed */
    int                       igf_curr[IGF_MAX_SFB];                   /* current igf energies */
    int                       igf_prev[IGF_MAX_SFB];                   /* needed for concealment or indepflag==0 */
    int                       igf_curr_subframe[IGF_MAX_SUBFRAMES][IGF_TRANS_FAK][IGF_MAX_SFB];    /* current igf energies per subframe*/
    int                       igf_prev_subframe[IGF_MAX_SUBFRAMES][IGF_MAX_SFB];                   /* needed for concealment or indepflag==0 */
    int                       igf_flatteningTrigger_subframe[IGF_MAX_SUBFRAMES];

    /* spectral whitening: */
    float                     pSpecFlat[IGF_START_MX];
    int                       currWhiteningLevel[IGF_MAX_TILES];
    int                       prevWhiteningLevel[IGF_MAX_TILES];                                /* needed for concealment */
    int                       currWhiteningLevel_subframe[IGF_MAX_SUBFRAMES][IGF_MAX_TILES];
    int                       prevWhiteningLevel_subframe[IGF_MAX_SUBFRAMES][IGF_MAX_TILES];    /* needed for concealment */

    float                     totalNoiseNrg;
    int                       n_noise_bands;

    float                     totalNoiseNrg_off;
    int                       n_noise_bands_off;

    /* IGF SCF decoding: */
    IGFSCFDEC_INSTANCE        hArithSCFdec;

    /* concealment: */
    int                       frameLossCounter;

} IGFDEC_PRIVATE_DATA,*IGF_DEC_PRIVATE_DATA_HANDLE;

typedef struct igfdec_instance_struct
{
    int             isIGFActive;
    int             infoIGFAllZero;
    int             infoIGFStopLine;
    int             infoIGFStartLine;
    int             infoIGFStopFreq;
    int             infoIGFStartFreq;
    unsigned char   infoTCXNoise[IGF_START_MX];
    int             flag_sparse[N_MAX_TCX-IGF_START_MN];
    float           virtualSpec[N_MAX_TCX-IGF_START_MN];
    int             flatteningTrigger;
    IGFDEC_PRIVATE_DATA igfData;
} IGFDEC_INSTANCE, *IGF_DEC_INSTANCE_HANDLE;

/*----------------------------------------------------------------------------------*
 *
 * Main decoder structure
 *
 *----------------------------------------------------------------------------------*/

typedef struct Decoder_State
{

    /*----------------------------------------------------------------------------------*
     * Common parameters
     *----------------------------------------------------------------------------------*/

    short codec_mode;                                   /* Mode 1 or 2 */
    short mdct_sw_enable;                               /* MDCT switching enable flag */
    short mdct_sw;                                      /* MDCT switching indicator */
    short last_codec_mode;                              /* last used codec mode */

    unsigned short bit_stream[MAX_BITS_PER_FRAME+16];
    short next_bit_pos;                                 /* position of the next bit to be read from the bitstream */
    short bitstreamformat;                              /* Bitstream format flag (G.192/MIME) */
    short amrwb_rfc4867_flag;                           /* MIME from rfc4867 is used */
    short BER_detect;                                   /* flag to signal detected bit error in the bitstream */
    int   output_Fs;                                    /* output sampling rate */
    long  total_brate;                                  /* total bitrate in kbps of the codec */
    long  last_total_brate;                             /* last total bitrate in kbps of the codec */
    long  last_total_brate_ber;                         /* last total bitrate in kbps of the codec - used only when first frame is lost and BER is detected afterwards */
    short core;                                         /* core (ACELP_CORE, TCX_20_CORE, TCX_10_CORE, HQ_CORE, AMR_WB_CORE) */
    long  core_brate;                                   /* core bitrate */
    long  last_core_brate;                              /* previous frame core bitrate */
    short extl;                                         /* extension layer */
    short last_extl;                                    /* previous extension layer */
    long  extl_brate;                                   /* extension layer bitrate */
    short L_frame;                                      /* ACELP core internal frame length */
    short bwidth;                                       /* encoded signal bandwidth */
    short Opt_AMR_WB;                                   /* flag indicating AMR-WB IO mode */
    short Opt_VOIP;                                     /* flag indicating VOIP mode with JBM */
    short ini_frame;                                    /* initialization frames counter */
    Word16 CNG;                                      /* RXDTX handler: CNG=1, nonCNG=0 */
    Word16 prev_ft_speech;                           /* RXDTX handler: previous frametype flag for  G.192 format AMRWB SID_FIRST detection */

    /*----------------------------------------------------------------------------------*
     * ACELP core parameters
     *----------------------------------------------------------------------------------*/

    float old_exc[L_EXC_MEM_DEC];                       /* old excitation */
    float old_excFB[L_FRAME48k];                        /* old excitation FB */
    float lsp_old[M];                                   /* old LSP vector at the end of the frame */
    float lsf_old[M];                                   /* old LSF vector at the end of the frame */
    unsigned int offset_scale1[MAX_NO_MODES+1][MAX_NO_SCALES+1]; /* offsets for LSF LVQ structure 1st 8-dim subvector*/
    unsigned int offset_scale2[MAX_NO_MODES+1][MAX_NO_SCALES+1]; /* offsets for LSF LVQ structure 2nd 8-dim subvector*/
    unsigned int offset_scale1_p[MAX_NO_MODES_p+1][MAX_NO_SCALES+1]; /* offsets for LSF LVQ structure, pred. case, 1st 8-dim subvector*/
    unsigned int offset_scale2_p[MAX_NO_MODES_p+1][MAX_NO_SCALES+1]; /* offsets for LSF LVQ structure, pred. case, 2nd 8-dim subvector*/
    short no_scales[MAX_NO_MODES][2];                     /* LSF LVQ structure */
    short no_scales_p[MAX_NO_MODES_p][2];                 /* LSF LVQ structure */
    float tilt_code;                                    /* tilt of code */
    float mem_syn1[M];                                  /* synthesis filter memory (for core switching and FD BWE) */
    float mem_syn2[M];                                  /* synthesis filter memory */
    float mem_syn3[M];
    float mem_deemph;                                   /* deemphasis filter memory */
    float mem_hp20_out[4];                              /* HP filter memory for synthesis */
    float mem_AR[M];                                    /* AR memory of LSF quantizer (past quantized LSFs without mean) */
    float mem_MA[M];                                    /* MA memory of LSF quantizer (past quantized residual) */
    float stab_fac;                                     /* LSF stability factor */
    float stab_fac_smooth;                              /* low-pass filtered stability factor */
    short last_coder_type;                              /* previous coder type */
    float agc_mem2[2];                                   /* memory of AGC for saturation control */
    float past_qua_en[GAIN_PRED_ORDER];                 /* gain quantization memory (used also in AMR-WB IO mode) */
    short mid_lsf_int;
    short safety_net;

    short seed_tcx ;                                    /* AC mode (GSC) - seed for noise fill*/
    short GSC_noisy_speech;                             /* AC mode (GSC) - flag to indicate GSC osn SWB noisy speech */
    short Last_GSC_noisy_speech_flag;                   /* AC mode (GSC) - mem of the past flag to indicate GSC osn SWB noisy speech */
    short cor_strong_limit;                             /* AC mode (GSC) - Indicator about high spectral correlation per band */
    float old_y_gain[MBANDS_GN];                        /* AC mode (GSC) - AR mem for low rate gain quantization  */
    short noise_lev;                                    /* AC mode (GSC) - noise level                             */
    float lt_ener_per_band[MBANDS_GN];
    float Last_frame_ener;                              /* AC mode (GSC) - last frame energy  */
    float Last_GSC_spectrum[L_FRAME];                   /* AC mode (GSC) - Last good GSC spectrum */
    short Last_GSC_pit_band_idx;                        /* AC mode (GSC) - Last pitch band index */
    float last_exc_dct_in[L_FRAME];                     /* AC mode (GSC) - previous excitation                        */
    float last_ener;                                    /* AC mode (GSC) - previous energy                           */
    short last_bitallocation_band[6];                   /* AC mode (GSC) - previous bit allocation of each band      */

    float gc_threshold;                                 /* Noise enhancer - threshold for gain_code */
    float dispMem[8];                                   /* Noise enhancer - phase dispersion algorithm memory */

    float prev_r;                                       /* HF BWE - previous sub-frame gain                */
    float fmerit_w_sm;                                  /* HF BWE - fmerit parameter memory */
    short frame_count;                                  /* HF BWE - frame count             */
    float ne_min;                                       /* HF BWE - minimum Noise gate - short-term energy */
    float fmerit_m_sm;                                  /* HF BWE - memory of fmerit_m param               */
    float voice_fac_amr_wb_hf;                          /* HF BWE - voice factor                           */
    float unvoicing;                                    /* HF BWE - unvoiced parameter                     */
    float unvoicing_sm;                                 /* HF BWE - smoothed unvoiced parameter            */
    short unvoicing_flag;                               /* HF BWE - unvoiced flag                          */
    short voicing_flag;                                 /* HF BWE - voiced flag                            */
    short start_band_old;                               /* HF BWE - previous start point for copying frequency band */
    float OptCrit_old;                                  /* HF BWE - previous criterion value for deciding the start point */

    short seed2;                                        /* HF (6-7kHz) BWE - seed for random signal generator */
    float mem_hp400[4];                                 /* HF (6-7kHz) BWE - hp400 filter memory */
    float mem_hf[(L_FIR-1)];                            /* HF (6-7kHz) BWE - band-pass filter memory */
    float mem_syn_hf[M];                                /* HF (6-7kHz) BWE - synthesis filter memory */
    float delay_syn_hf[NS2SA(16000,DELAY_CLDFB_NS)];      /* HF (6-7kHz) BWE - To synchronise BWE content with postfiltered synthesis */
    float mem_hp_interp[INTERP_3_1_MEM_LEN];            /* HF (6-7 kHz) BWE - interp. memory */

    short unv_cnt;                                      /* Stationary noise UV modification - unvoiced frame counter */
    short uv_count;                                     /* Stationary noise UV modification - unvoiced counter */
    short act_count;                                    /* Stationary noise UV modification - activation counter */
    float ge_sm;                                        /* Stationary noise UV modification - smoothed excitation gain */
    float lspold_s[M];                                  /* Stationary noise UV modification - old LSP vector */
    short noimix_seed;                                  /* Stationary noise UV modification - mixture seed */
    float min_alpha;                                    /* Stationary noise UV modification - minimum alpha */
    float exc_pe;                                       /* Stationary noise UV modification - memory of the preemphasis filter */

    short bfi;                                          /* FEC - bad frame indicator */
    short prev_bfi;                                     /* FEC - previous bad frame indicator */
    short prev_old_bfi;                                 /* FEC - previous old bad frame indicator */
    short seed;                                         /* FEC - seed for random generator for excitation */
    float lp_ener_bfi;                                  /* FEC - long-term active-signal average energy */
    short last_good;                                    /* FEC - clas of last good received */
    float lp_gainp;                                     /* FEC - low-pass filtered pitch gain */
    float lp_gainc;                                     /* FEC - low-pass filtered code gain */
    float lp_ener;                                      /* FEC - low-pass filtered energy */
    float enr_old;                                      /* FEC - energy of the concealed frame */
    float bfi_pitch;                                    /* FEC - pitch for FEC */
    short bfi_pitch_frame;                              /* FEC - frame length when pitch for FEC is saved */
    float old_pitch_buf[2*NB_SUBFR16k+2];               /* FEC - buffer of old subframe pitch values */
    short upd_cnt;                                      /* FEC - counter of frames since last update */
    short scaling_flag;                                 /* FEC - flag to indicate energy control of syn  */
    float lp_ener_FEC_av;                               /* FEC - averaged voiced signal energy           */
    float lp_ener_FEC_max;                              /* FEC - averaged voiced signal energy           */
    float old_enr_LP;                                   /* FEC - LP filter gain */
    short prev_nbLostCmpt;                              /* FEC - compt for number of consecutive lost frame at the previous frame*/
    short mode_lvq;                                     /* FEC - index for LSF mean vector */
    float lsfoldbfi0[M];                                /* FEC - LSF vector of the previous frame */
    float lsfoldbfi1[M];                                /* FEC - LSF vector of the past previous frame */
    float lsf_adaptive_mean[M];                         /* FEC - adaptive mean LSF vector for FEC */
    short decision_hyst;                                /* FEC - hysteresis of the music/speech decision */
    float old_exc2[L_EXC_MEM];                          /* FEC - old excitation2 used in fast recovery */
    float old_syn2[L_EXC_MEM];                          /* FEC - old syn speech used in fast recovery */
    short relax_prev_lsf_interp;
    float mem_syn_clas_estim[L_SYN_MEM_CLAS_ESTIM];     /* FEC - memory of the synthesis signal for frame class estimation */
    float tilt_swb_fec;                                 /* FEC - SWB TBE TILT */

    short cng_seed;                                     /* DTX/CNG - seed for white noise random generator */
    float lspCNG[M];                                    /* DTX/CNG - LP filtered ISPs */
    short first_CNG;                                    /* DTX/CNG - first CNG frame flag */
    float Enew;                                         /* DTX/CNG - decoded residual energy */
    short old_enr_index;                                /* DTX/CNG - index of last encoded CNG energy */
    short cng_ener_seed;                                /* DTX/CNG - seed for random generator for variation of excitation energy */
    short cng_ener_seed1;
    short last_allow_cn_step;
    short ho_hist_size;                                 /* DTX/CNG - size of DTX hangover history buffer for averaging, <0,HO_HIST_SIZE> */
    short ho_hist_ptr;                                  /* DTX/CNG - pointer for averaging buffers */
    long  ho_sid_bw;                                    /* DTX/CNG - SID bandwidth flags */
    float ho_lsp_hist[HO_HIST_SIZE*M];                  /* DTX/CNG - old LSP buffer for averaging */
    float ho_ener_hist[HO_HIST_SIZE];                   /* DTX/CNG - energy buffer for averaging */
    float ho_env_hist[HO_HIST_SIZE*NUM_ENV_CNG];
    short act_cnt;                                      /* DTX/CNG - counter of active frames */
    short ho_circ_size;                                 /* DTX/CNG - size of DTX hangover history buffer for averaging, <0,HO_HIST_SIZE> */
    short ho_circ_ptr;                                  /* DTX/CNG - pointer for averaging buffers */
    float ho_lsp_circ[HO_HIST_SIZE*M];                  /* DTX/CNG - old LSP buffer for averaging */
    float ho_ener_circ[HO_HIST_SIZE];                   /* DTX/CNG - energy buffer for averaging */
    float ho_env_circ[HO_HIST_SIZE*NUM_ENV_CNG];
    short num_ho;                                       /* DTX/CNG - number of selected hangover frames */
    short ho_16k_lsp[HO_HIST_SIZE];                     /* DTX/CNG - 16k LSPs flags */
    short CNG_mode;                                     /* DTX/CNG - mode for DTX configuration */
    long last_active_brate;                             /* DTX/CNG - last active frame bitrate used for CNG_mode control */
    short last_CNG_L_frame;                             /* DTX/CNG - last CNG frame length */
    short act_cnt2;                                     /* DTX/CNG - counter of active frames for CNG_mode switching */
    short cng_type;                                     /* DTX/CNG - flag indicating LP or CLDFB based SID/CNG */
    short last_cng_type;                                /* DTX/CNG - flag indicating last frame LP or CLDFB based SID/CNG */
    float old_env[20];
    float lp_env[20];
    float exc_mem[24];
    float exc_mem1[30];

    short bpf_off;                                      /* Bass post-filter - do not use BPF when this flag is set to 1 */
    float pst_old_syn[NBPSF_PIT_MAX];                   /* Bass post-filter - old synthesis buffer 1 */
    float pst_mem_deemp_err;                            /* Bass post-filter - filter memory of noise LP filter */
    float pst_lp_ener;                                  /* Bass post-filter - long-term energy */
    short Track_on_hist[L_TRACK_HIST];                  /* Bass post-filter - History of half frame usage */
    short vibrato_hist[L_TRACK_HIST];                   /* Bass post-filter - History of frames declared as vibrato */
    float psf_att;                                      /* Bass post-filter - post filter attenuation factor */
    float mem_mean_pit[L_TRACK_HIST];                   /* Bass post-filter - average pitch memory */

    HANDLE_CLDFB_FILTER_BANK cldfbAna;                  /* main analysis filter bank handle */
    HANDLE_CLDFB_FILTER_BANK cldfbBPF;                  /* BPF analysis filter bank handle */
    HANDLE_CLDFB_FILTER_BANK cldfbSyn;                  /* main synthesis  filter bank handle */

    short last_active_bandsToZero_bwdec;
    short flag_NB_bwddec;
    short last_flag_filter_NB;
    float perc_bwddec;
    int active_frame_cnt_bwddec;
    short flag_buffer[20];
    int total_frame_cnt_bwddec;
    float avg_nrg_LT;
    float ng_ener_ST;                                   /* Noise gate - short-term energy */

    short last_L_frame;                                 /* ACELP@16kHz - last value of st->L_frame */
    float mem_preemp_preQ;                              /* ACELP@16kHz - prequantizer preemhasis memory */
    short last_nq_preQ;                                 /* ACELP@16kHz - AVQ subquantizer number of the last sub-band of the last subframe  */
    short use_acelp_preq;                               /* ACELP@16kHz - flag of prequantizer usage */

    /* Improvement of unvoiced and audio signals in AMR-WB IO mode */
    short UV_cnt;                                       /* number of consecutives frames classified as UV */
    float LT_UV_cnt;                                    /* long-term consecutives frames classified as UV */
    float Last_ener;                                    /* last_energy frame                              */
    float lt_diff_etot[MAX_LT];                         /* stability estimation - long-term total energy variation */
    float old_Aq[NB_SUBFR*(M+1)];                       /* old LPC filter coefficient                     */
    float lt_voice_fac;                                 /* average voice factor over 4 sub-frames         */

    /* NB and formant post-filter */
    PFSTAT pfstat;                                      /* NB and formant post-filter states    */
    float psf_lp_noise;                                 /* NB post-filter - long-term noise     */

    /*----------------------------------------------------------------------------------*
     * SC-VBR
     *----------------------------------------------------------------------------------*/

    /* PPP decoder variables */
    short last_ppp_mode_dec;
    short ppp_mode_dec;
    short last_nelp_mode_dec;
    short nelp_mode_dec;
    int firstTime_voiceddec;

    /* DTFS variables */
    float dtfs_dec_a[MAXLAG_WI];
    float dtfs_dec_b[MAXLAG_WI];
    int   dtfs_dec_lag;
    int   dtfs_dec_nH;
    int   dtfs_dec_nH_4kHz;
    float dtfs_dec_upper_cut_off_freq_of_interest;
    float dtfs_dec_upper_cut_off_freq;
    float ph_offset_D;
    float lastLgainD;                                   /* previous gain value for the low band */
    float lastHgainD;                                   /* previous gain value for the high band */
    float lasterbD[NUM_ERB_WB];                         /* previous amplitude spectrum (ERB) */

    /* NELP decoder variables */
    float bp1_filt_mem_nb_dec[14];
    float bp1_filt_mem_wb_dec[8];
    float shape1_filt_mem_dec[20];
    float shape2_filt_mem_dec[20];
    float shape3_filt_mem_dec[20];

    short nelp_dec_seed;
    float FadeScale;
    float prev_gain_pit_dec;
    float prev_tilt_code_dec;
    short vbr_hw_BWE_disable_dec;
    short last_vbr_hw_BWE_disable_dec;

    /*----------------------------------------------------------------------------------*
    * channel-aware mode
    *----------------------------------------------------------------------------------*/

    float tilt_code_dec[NB_SUBFR16k];

    short rf_frame_type;
    short use_partial_copy;
    short prev_use_partial_copy;
    short rf_flag;
    short rf_flag_last;

    short rf_fec_offset;
    short next_coder_type;
    short prev_rf_frame_type;
    short rf_target_bits;

    short rf_indx_nelp_fid;
    short rf_indx_nelp_iG1;
    short rf_indx_nelp_iG2[2];
    short rf_indx_tbeGainFr;

    /*----------------------------------------------------------------------------------*
     * HR SWB BWE parameters
     *----------------------------------------------------------------------------------*/

    short bwe_highrate_seed;
    float t_audio_prev[2*END_FREQ_BWE_FULL_FB/50 - NUM_NONTRANS_START_FREQ_COEF];
    short old_is_transient_hr_bwe;
    float mem_EnergyLT;


    /*----------------------------------------------------------------------------------*
     * HQ core parameters
     *----------------------------------------------------------------------------------*/

    float synth_history[L_FRAME48k+OLD_SYNTH_SIZE_DEC+NS2SA(48000, PH_ECU_LOOKAHEAD_NS)];            /* unified synthesis memory */
    float *old_synthFB;
    float old_out[L_FRAME48k];                          /* HQ core - previous synthesis for OLA */

    float old_outLB[L_FRAME32k];
    float old_coeffs[L_FRAME8k];                        /* HQ core - old coefficients (for FEC) */
    float oldIMDCTout[L_FRAME8k/2];
    float prev_oldauOut[L_FRAME8k];
    float diff_energy;
    short stat_mode_out;
    short stat_mode_old;
    short phase_mat_flag;
    short phase_mat_next;
    short old_Min_ind;
    float old_auOut_2fr[L_FRAME8k*2];
    short old_is_transient[3];                          /* HQ core - previous transient flag (for FEC)  */
    float old_out_pha[2][N_LEAD_NB];                    /* FEC for HQ Core, 0-phase matching old_out, 1-overlapping original old_out and  phase matching old_out*/
    short old_bfi_cnt;                                  /* HQ core - # of bfi until previous frame(for FEC)  */
    float ynrm_values[MAX_SB_NB][MAX_PGF];
    float r_p_values[MAX_SB_NB][MAX_ROW];
    float Norm_gain[SFM_N_NB];
    short HQ_FEC_seed;
    float energy_MA_Curr[2];

    short last_core;
    short prev_last_core;
    short last_hq_core_type;
    short last_L_frame_ori;
    float previoussynth[L_FRAME48k];
    float old_synth_sw[NS2SA(48000,FRAME_SIZE_NS-ACELP_LOOK_NS-DELAY_BWE_TOTAL_NS)];
    float delay_buf_out[HQ_DELTA_MAX*HQ_DELAY_COMP];
    short mem_norm[SFM_N_ENV_STAB];
    float mem_env_delta;
    short no_att_hangover;
    float energy_lt;
    short hq_generic_seed;
    float prev_noise_level[2];
    short prev_hqswb_clas;
    short prev_R;                                       /* the table of bit allocation of last frame  */
    float prev_coeff_out[L_HQ_WB_BWE];                  /* the highest coefficients of last frame     */
    short prev_SWB_peak_pos[SPT_SHORTEN_SBNUM];
    float old_Aq_12_8[M+1];                             /* old Aq[] for core switching */
    float old_Es_pred;                                  /* old Es_pred for core switching */

    short HqVoicing;
    float fer_samples[L_FRAME48k];
    float prev_normq[SFM_N_WB];                         /* previous norms                             */
    float prev_env[SFM_N_WB];                           /* previous noise envelopes                   */

    float last_ni_gain[BANDS_MAX];
    float last_env[BANDS_MAX];
    short last_max_pos_pulse;

    /* pre-echo reduction */
    float memfilt_lb;
    float mean_prev_hb;
    float smoothmem;
    float mean_prev;
    float mean_prev_nc;
    float wmold_hb;
    short prevflag;
    short pastpre;
    short prev_frm_hfe2;
    short prev_stab_hfe2;
    float prev_ni_ratio;
    float prev_En_sb[NB_SWB_SUBBANDS];

    /* PVQ range coder state */
    unsigned int rc_low;
    unsigned int rc_range;
    unsigned int rc_help;
    short rc_num_bits;
    short rc_offset;
    short rc_end;

    /*----------------------------------------------------------------------------------*
     * TBE parameters
     *----------------------------------------------------------------------------------*/

    /* states for the filters used in generating SHB excitation from WB excitation */
    float state_lpc_syn[LPC_SHB_ORDER];
    float mem_csfilt [2];

    /* states for the filters used in generating SHB signal from SHB excitation*/
    float state_syn_shbexc[L_SHB_LAHEAD];
    float syn_overlap[L_SHB_LAHEAD];                    /* overlap buffer used to Adjust SHB Frame Gain*/

    /* previous frame parameters for frame error concealment */
    float lsp_prevfrm[ LPC_SHB_ORDER];
    float GainFrame_prevfrm;
    float GainShape_Delay[NUM_SHB_SUBFR/2];
    float GainAttn;

    float old_bwe_exc[PIT16k_MAX * 2];       /* old excitation */
    short bwe_seed[2];
    float bwe_non_lin_prev_scale;
    float old_bwe_exc_extended[NL_BUFF_OFFSET];
    float last_voice_factor;

    float genSHBsynth_Hilbert_Mem[HILBERT_MEM_SIZE];

    float mem_genSHBexc_filt_down_shb[(2*ALLPASSSECTIONS_STEEP+1)];
    float mem_genSHBexc_filt_down_wb2[(2*ALLPASSSECTIONS_STEEP+1)];
    float mem_genSHBexc_filt_down_wb3[(2*ALLPASSSECTIONS_STEEP+1)];
    float genSHBsynth_state_lsyn_filt_shb_local[ 2 * ALLPASSSECTIONS_STEEP ];
    float state_lsyn_filt_shb[ 2 * ALLPASSSECTIONS_STEEP];
    float state_lsyn_filt_dwn_shb[ 2 * ALLPASSSECTIONS_STEEP];
    float mem_resamp_HB[INTERP_3_1_MEM_LEN];
    float mem_resamp_HB_32k[2*ALLPASSSECTIONS_STEEP+1];
    float prev_synth_buffer[NS2SA(48000,DELAY_BWE_TOTAL_NS - DELAY_CLDFB_NS)];
    float hb_prev_synth_buffer[NS2SA(48000, DELAY_BWE_TOTAL_NS)];
    short old_bwe_delay;

    short syn_dm_phase;
    float fbbwe_hpf_mem[4][4];
    float prev_wb_bwe_frame_pow;
    float prev_swb_bwe_frame_pow;
    float prev_ener;
    float prev_GainShape;
    float fb_state_lpc_syn[LPC_SHB_ORDER];
    float fb_tbe_demph;
    float prev_fbbwe_ratio;

    /* WB/SWB bandwidth switching */
    float tilt_wb;
    float tilt_swb;
    float prev_ener_shb;
    float enerLH;
    float prev_enerLH;
    float enerLL;
    float prev_enerLL;
    short prev_fractive;
    short prev_bws_cnt;
    short bws_cnt;
    short bws_cnt1;
    float attenu1;
    short last_inner_frame;
    short last_bwidth;
    float prev_weight1;
    float t_audio_q[L_FRAME];
    float tbe_demph;
    float tbe_premph;
    float mem_stp_swb[LPC_SHB_ORDER];
    float *ptr_mem_stp_swb;
    float gain_prec_swb;
    float mem_zero_swb[LPC_SHB_ORDER];

    float swb_lsp_prev_interp[LPC_SHB_ORDER];
    float prev1_shb_ener_sf, prev2_shb_ener_sf, prev3_shb_ener_sf, prev_res_shb_gshape, prev_mixFactors;
    float tilt_mem;                                     /* Formant factor adaptation tilt smoothing memory */
    float prev_lsf_diff[LPC_SHB_ORDER-2];
    float prev_tilt_para;
    float cur_sub_Aq[M+1];

    /* quantized data */
    short lsf_idx[NUM_Q_LSF];
    short m_idx;
    short grid_idx;
    short idxSubGains;
    short idxFrameGain;
    short idx_shb_fr_gain;
    short idx_res_gs[NB_SUBFR16k];
    short idx_mixFac;

    short lsf_WB;
    short gFrame_WB;

    short idxGain;

    float old_core_synth[L_FRAME16k];
    float old_tbe_synth[L_SHB_TRANSITION_LENGTH];

    float int_3_over_2_tbemem_dec[INTERP_3_2_MEM_LEN];
    float interpol_3_2_cng_dec[INTERP_3_2_MEM_LEN];

    /*----------------------------------------------------------------------------------*
     * SWB BWE parameters
     *----------------------------------------------------------------------------------*/

    float old_wtda_swb[L_FRAME48k];
    float old_syn_12k8_16k[NS2SA(16000, DELAY_FD_BWE_ENC_NS)];
    float mem_deemph_old_syn;
    short prev_mode;
    float prev_SWB_fenv[SWB_FENV];
    float prev_Energy;
    float prev_Energy_wb;
    short prev_L_swb_norm;
    short Seed;
    short prev_frica_flag;
    float mem_imdct[L_FRAME48k];
    float prev_td_energy;
    float prev_weight;
    short prev_coder_type;
    short prev_flag;
    float last_wb_bwe_ener;
    float prev_fb_ener_adjust;

    /*----------------------------------------------------------------------------------*
     * SWB DTX/CNG parameters
     *----------------------------------------------------------------------------------*/

    float shb_cng_ener;
    float wb_cng_ener;
    float last_wb_cng_ener;
    float last_shb_cng_ener;
    short swb_cng_seed;
    float lsp_shb_prev_prev[LPC_SHB_ORDER];
    float lsp_shb_prev[LPC_SHB_ORDER];
    short shb_dtx_count;
    short last_vad;
    short trans_cnt;
    short burst_cnt;
    float last_shb_ener;

    /*----------------------------------------------------------------------------------*
     * HQ FEC
     *----------------------------------------------------------------------------------*/
    float *prev_good_synth;
    short prev_sign_switch[HQ_FEC_SIGN_SFM];
    short prev_sign_switch_2[HQ_FEC_SIGN_SFM];

    /* HQ PHASE ECU internal state */
    short time_offs;
    float X_sav[PH_ECU_SPEC_SIZE];
    short num_p;
    short plocs[MAX_PLOCS];
    float plocsi[MAX_PLOCS];
    float env_stab;
    short mem_norm_hqfec[SFM_N_ENV_STAB];
    float mem_env_delta_hqfec;
    float env_stab_plc;
    float env_stab_state_p[NUM_ENV_STAB_PLC_STATES];
    short envstabplc_hocnt;

    float mag_chg_1st[LGW_MAX];          /* i/o: per band magnitude modifier for transients*/
    float Xavg[LGW_MAX];                 /* Frequency group average gain to fade to   */
    float beta_mute;                     /* Factor for long-term mute            */

    short last_fec;
    short ph_ecu_HqVoicing;
    short oldHqVoicing;
    float oldgapsynth[L_FRAME48k];
    short ph_ecu_active;      /* Set to 1 if Phase ECU was used in last bad frame; Set to 2 if TCX TD PLC was used */
    short ni_seed_forfec;
    short ber_occured_in_pvq;    /* flag for BER detection from PVQ routines */

    /*----------------------------------------------------------------------------------*
     * LD music post-filter
     *----------------------------------------------------------------------------------*/

    float LDm_mem_etot;                               /* LD music post-filter - total energy memory  */
    short LDm_last_music_flag;                        /* LD music post-filter - last music flag */
    short LDm_nb_thr_1;                               /* LD music post-filter - number of consecutive frames of level 1 */
    short LDm_nb_thr_3;
    float dct_post_old_exc[DCT_L_POST-OFFSET2];
    float LDm_thres[4];                               /* LD music post-filter - Classification threshold */
    float LDm_lt_diff_etot[MAX_LT];                   /* LD music post-filter - long-term total energy variation */
    float LDm_enh_lp_gbin[VOIC_BINS_HR];              /* LD music post-filter - smoothed suppression gain, per bin FFT */
    float LDm_enh_lf_EO[VOIC_BINS_HR];                /* LD music post-filter - old per bin E for previous half frame */
    float LDm_enh_min_ns_gain;                        /* LD music post-filter - minimum suppression gain */
    float LDm_bckr_noise[MBANDS_GN_LD];               /* LD music post-filter - background noise estimation per critical band */
    float filt_lfE[DCT_L_POST];
    short last_nonfull_music;


    /*ACELP config*/
    short force_lpd_reset;
    ACELP_config acelp_cfg;       /*configuration set for each frame*/

    ACELP_config acelp_cfg_rf; /* configuration for RF frame */

    /*TCX config*/
    TCX_config tcx_cfg;
    int L_frameTCX;

    /*dec_prm.c*/
    int bits_frame;               /*bit per frame overall */
    int bits_frame_core;          /*bit per frame for the core*/
    int narrowBand;

    int last_is_cng;

    float old_syn_Overl[L_FRAME32k/2];

    float syn_Overl_TDAC[L_FRAME32k/2];
    float syn_Overl_TDACFB[L_FRAME_MAX/2];

    float syn_Overl[L_FRAME32k/2];
    float syn_OverlFB[L_FRAME_MAX/2];

    float *acelp_zir;
    float old_synth[OLD_SYNTH_INTERNAL_DEC];/* synthesis memory                 */
    int   old_synth_len;
    int   old_synth_lenFB;
    float syn[M+1];

    /* bass_pf.c */
    int   bpf_gain_param;                      /* bass post-filter gain factor parameter (0->noBpf)*/

    int   L_frame_past;
    int   L_frameTCX_past;

    float lsfold_uw[M];                  /* old lsf (unweighted) */
    float lspold_uw[M];                  /* old lsp (unweighted) */
    short seed_tcx_plc;                  /* seed memory (for random function in TCX PLC) */
    float past_gpit;                     /* past gain of pitch (for frame recovery) */
    float past_gcode;                    /* past energy (!) of code  (for frame recovery) */
    float lsf_cng[M];                    /* lsf coefficients used for CNG generation (long term) */
    float lspold_cng[M];                 /* lsp coefficients used for CNG generation (long term) */
    float lsp_q_cng[M];                  /* lsp coefficients used for CNG generation (short term interpolated) */
    float old_lsp_q_cng[M];              /* lsp coefficients used for CNG generation (short term interpolated) */
    float lsf_q_cng[M];                  /* lsf coefficients used for CNG generation (short term interpolated) */
    float old_lsf_q_cng[M];              /* lsf: old quantized lsfs for background noise */
    float Aq_cng[(NB_SUBFR16k+1)*(M+1)]; /* LPC coefficients derived from CNG estimate  */
    float mem_syn_unv_back[M];           /* filter memory for unvoiced synth */
    int   plcBackgroundNoiseUpdated;     /* flag: Is background noise estimate updated? */
    float last_gain_syn_deemph;
    float last_concealed_gain_syn_deemph;

    int enableTcxLpc; /* global toggle for the TCX LPC quantizer */
    int envWeighted;  /* are is{p,f}_old[] weighted or not? */

    /* variables for framing */
    int nb_subfr;

    int fscale;
    int fscale_old;
    int sr_core;

    int pit_min;
    int pit_fr1;
    int pit_fr1b;
    int pit_fr2;
    int pit_max;
    int pit_res_max;
    int pit_res_max_past;

    int pit_min_TCX;
    int pit_max_TCX;

    /*Preemphasis factor*/
    float preemph_fac;
    float gamma;

    /*for AMR-WB like 6.4 to 7 kHz upsampling and noise filling*/
    float mem_Aq[NB_SUBFR16k*(M+1)];

    /* Error concealment */
    int   last_core_bfi;                    /* PLC - mode in previous frame                                 */
    int   nbLostCmpt;                       /* PLC - compt for number of consecutive lost frame             */
    int   noise_filling_index;              /* PLC - last decoded noise filling index                       */
    float old_fpitch;                       /* PLC - last pitch of previous frame (as transmitted)          */
    float old_fpitchFB;                     /* PLC - last pitch of previous FB frame (depends on output sr) */
    short clas_dec;                         /* PLC - frame class at the decoder                             */
    float mem_pitch_gain[2*NB_SUBFR16k+2];  /* PLC - Pitch gain memory                                      */
    short plc_use_future_lag;               /* PLC - flag indicating if info (pitch lag / pitch gain) about
                                                     future frame is usable */
    int   prev_widow_left_rect;
    float CngLevelBackgroundTrace_bfi;       /* PLC - long term gain estimate for background level, used for PLC fade out */
    /* state variables for the minimum statistics used for PLC */
    float NoiseLevelMemory_bfi[PLC_MIN_STAT_BUFF_SIZE];
    int   NoiseLevelIndex_bfi;
    int   CurrLevelIndex_bfi;
    float LastFrameLevel_bfi;
    float old_gaintcx_bfi;
    float cummulative_damping;
    float cngTDLevel;
    float conceal_eof_gain;
    float damping;
    float gainHelper;
    float stepCompensate;
    int reset_mem_AR;
    short tcxConceal_recalc_exc;

    short rate_switching_init;

    /* LPC quantization */
    int lpcQuantization;
    int numlpc;

    /* Bandwidth */
    float TcxBandwidth;

    float voice_fac;

    int tcxonly;

    /*TCX resisual Q*/
    int resQBits[NB_DIV];       /* number of bits read for the residual Quantization in TCX*/

    int last_ctx_hm_enabled;

    /* TCX-LTP */
    int tcxltp;                 /* TCX LTP enabled for operation mode */
    float tcxltp_gain;
    int tcxltp_pitch_int;
    int tcxltp_pitch_fr;

    float tcxltp_mem_in[TCXLTP_MAX_DELAY];
    float tcxltp_mem_out[L_FRAME48k];
    int tcxltp_pitch_int_post_prev;
    int tcxltp_pitch_fr_post_prev;
    float tcxltp_gain_post_prev;
    int tcxltp_filt_idx_prev;

    struct tonalmdctconceal tonalMDCTconceal;
    int tonal_mdct_plc_active;
    int last_tns_active;
    int second_last_tns_active;
    float cummulative_damping_tcx;
    int second_last_core;
    float tcxltp_second_last_pitch;
    float tcxltp_third_last_pitch;
    float tcxltp_last_gain_unmodified;

    float FBTCXdelayBuf[111]; /* 2.3125ms at 48kHz -> 111 samples */

    /* parameters for switching */
    float mem_syn_r[L_SYN_MEM];         /*LPC synthesis memory needed for rate switching*/
    short rate_switching_reset;

    float bpf_noise_buf[L_FRAME16k];
    float *p_bpf_noise_buf;

    int   enableGplc;
    int   flagGuidedAcelp;
    int   T0_4th;
    int   guidedT0;

    short enablePlcWaveadjust;
    short tonality_flag;
    T_PLCInfo plcInfo;

    short VAD;
    short flag_cna;
    short last_flag_cna;

    float lp_noise;

    short seed_acelp;
    int core_ext_mode; /*GC,VC,UC,TC: core extended mode used for PLC or Acelp-external modules.*/

    short dec_glr;
    short dec_glr_idx;

    short tcx_hm_LtpPitchLag;
    short tcx_lpc_shaped_ari;

    DEC_MODE         m_decodeMode;
    unsigned char    m_frame_type;                       /*ZERO_FRAME/SID_FRAME/ACTIVE_FRAME*/
    unsigned char    m_old_frame_type;                   /*ZERO_FRAME/SID_FRAME/ACTIVE_FRAME*/

    /*Frequency-domain-based CNG*/
    HANDLE_FD_CNG_DEC            hFdCngDec;

    IGFDEC_INSTANCE  hIGFDec;
    short igf;

    short tec_tfa;
    short tec_flag;
    short tfa_flag;
    TEMPORAL_ENVELOPE_CODING_DECODER tecDec;

    short old_ppp_mode;
    float old_hb_synth[L_FRAME48k];
    short con_tcx;
    short last_con_tcx;

    short writeFECoffset;


} Decoder_State;

#endif

