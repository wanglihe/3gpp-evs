/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#ifndef STAT_ENC_H
#define STAT_ENC_H

#include "options.h"
#include "stat_com.h"          /* Common structures                      */
#include "cnst.h"


/*------------------------------------------------------------------------------------------*
 * Indice
 *------------------------------------------------------------------------------------------*/

typedef struct
{
    unsigned short value;        /* value of the quantized indice */
    short nb_bits;      /* number of bits used for the quantization of the indice */
} Indice;

typedef struct
{
    /* signal memory */
    float syn[1+M];               /* Synthesis memory (non-pe)                            */

    /* ACELP memories*/
    float old_exc[L_EXC_MEM];     /* ACELP exc memory (Aq)                                */
    float mem_w0;
    float mem_syn[M];             /* ACELP synthesis memory (pe) before post-proc         */
    float mem_syn2[M];            /* ACELP synthesis memory (pe) after post-proc          */
    float mem_syn_r[L_SYN_MEM];   /* ACELP synthesis memory for 1.25ms                    */
    float mem_syn3[M];

    float tilt_code;
    float gc_threshold;           /* Noise enhancer - threshold for gain_code             */

    /* TCX memory */
    float Txnq[L_FRAME32k/2+64];/* Q target (overlap or ACELP+ZIR, use Aq)              */
    float *acelp_zir;
    float tcx_target_bits_fac;


} LPD_state;

typedef struct PLC_ENC_EVS
{

    int nBits;                  /* number of bits */

    int enableGplc;
    int T0_4th;
    int T0;
    int calcOnlylsf;
    int pit_min;
    int pit_max;

    float mem_MA[M];
    float mem_AR[M];

    float lsfold[M];            /* old lsf (frequency domain)               */
    float lspold[M];            /* old lsp (immittance spectral pairs)      */

    float lsfoldbfi0[M];        /* Previous frame lsf                       */
    float lsfoldbfi1[M];        /* Past previous frame lsf                  */
    float lsf_adaptive_mean[M]; /*  Mean lsf for bfi cases                  */
    float stab_fac;

    LPD_state *LPDmem;
    float old_exc[8];           /* ACELP exc memory (Aq)                     */

    float lsf_con[M];
    float last_lsf_ref[M];
    float last_lsf_con[M];
} PLC_ENC_EVS;

/* Arrays and variables specific to encoder */
typedef struct
{
    HANDLE_FD_CNG_COM hFdCngCom;

    float   msPeriodog[NPART]; /* Periodogram */
    float   msBminWin[NPART];
    float   msBminSubWin[NPART];
    float   msPsd[NPART];          /* Power Spectral Density estimate (i.e., smoothed periodogram) */
    float   msAlpha[NPART];        /* Optimal smoothing parameter */
    float   msMinBuf[MSNUMSUBFR*NPART];       /* Buffer of minima */
    float   msCurrentMinOut[NPART];
    float   msCurrentMin[NPART];
    float   msCurrentMinSubWindow[NPART];
    int     msLocalMinFlag[NPART];
    int     msNewMinFlag[NPART];
    float   msPsdFirstMoment[NPART];
    float   msPsdSecondMoment[NPART];
    float   msNoiseFloor[NPART];   /* Estimated noise floor */
    float   msNoiseEst[NPART];     /* Estimated noise level */
    float   energy_ho[NPART];
    float   msNoiseEst_old[NPART];
    float   msLogPeriodog[NPART]; /* Periodogram */
    float   msLogNoiseEst[NPART]; /* Estimated noise level */

    float   msPeriodogBuf[MSBUFLEN*NPART];
    int     msPeriodogBufPtr;

    int     stopFFTbinDec;
    int     startBandDec;
    int     stopBandDec;
    int     npartDec;
    int     midbandDec[NPART];
    int     nFFTpartDec;
    int     partDec[NPART];
}
FD_CNG_ENC;
typedef FD_CNG_ENC *HANDLE_FD_CNG_ENC;

/* transient_detection.h */
/** Delay buffer.
 * Used to buffer input samples and to define the subblock size of a transient detector.
 */
typedef struct
{
    /** Subblock size of a transient detector that uses this delay buffer. */
    int   nSubblockSize;
    /** Delay buffer */
    float buffer[L_FRAME_MAX/NSUBBLOCKS];
    /** Size of the delay buffer in use. Maximum delay from all users of this buffer. */
    int   nDelay;
} DelayBuffer;

/** Subblock energies.
 * Holds subblock energies and recursively accumulated energies.
 * Also buffers the energies.
 */
typedef struct
{
    /** Delay buffer. */
    DelayBuffer * pDelayBuffer;
    /** Subblock energies with a delay buffering. */
    float subblockNrg[NSUBBLOCKS+MAX_TD_DELAY];
    /** Recursively accumulated subblock energies with a delay buffering.
        At index i the value corresponds to the accumulated subblock energy up to i-1,
        including block i-1 and without block i. */
    float accSubblockNrg[NSUBBLOCKS+MAX_TD_DELAY+1];
    /** subblockNrgChange[i] = max(subblockNrg[i]/subblockNrg[i-1], subblockNrg[i-1]/subblockNrg[i]) */
    float subblockNrgChange[NSUBBLOCKS+MAX_TD_DELAY];
    /** Size of the delay buffer in use, as number of subblocks. Maximum delay from all users of this buffer. */
    int   nDelay;
    /* Delay of the input (modulo pDelayBuffer->nSubblockSize), nPartialDelay <= pDelayBuffer->nDelay. */
    int   nPartialDelay;

    /** Decay factor for the recursive accumulation */
    float facAccSubblockNrg;

    /** High-pass filter states (delay line) */
    float firState1;
    float firState2;

} SubblockEnergies;

struct TransientDetector;

/** Attack detection function.
 * Definition of a function used for checking the presence of an attack, given subblock energies, accumulated subblock energies and a threshold.
 * @param pSubblockNrg Subblock energies.
 * @param pAccSubblockNrg Recursively accumulated subblock energies.
 At index i the value corresponds to the accumulated subblock energy up to i-1,
 including block i-1 and without block i.
 * @param nSubblocks Number of subblocks available (those with an index >= 0). Subblocks from 0 to NSUBBLOCKS-1 correspond to the current frame.
 * @param nPastSubblocks Number of past subblocks available (those with a negative index).
 * @param attackRatioThreshold Attack ratio threshold.
 * @param pbIsAttackPresent Pointer to an output variable that will be set to TRUE if an attack is found, otherwise set to FALSE.
 * @param pAttackIndex Pointer to an output variable that will hold an attack position.
 */
typedef void (* TCheckSubblocksForAttack)(float const * pSubblockNrg, float const * pAccSubblockNrg, int nSubblocks, int nPastSubblocks, float attackRatioThreshold, int * pbIsAttackPresent, int * pAttackIndex);

/** Transient detector.
 */
typedef struct TransientDetector
{
    /** Subblock energies used in this transient detector. */
    SubblockEnergies * pSubblockEnergies;
    /* Delay of the transient detector in number of subblocks, nDelay <= pSubblockEnergies->nDelay. */
    int   nDelay;
    /** Number of subblocks to check for transients. */
    int   nSubblocksToCheck;
    /** Function for checking a presence of an attack. */
    TCheckSubblocksForAttack CheckSubblocksForAttack;
    /** Attack ratio threshold. */
    float attackRatioThreshold;
    /** True when an attack was detected. */
    int   bIsAttackPresent;
    /** The index of an attack. */
    int   attackIndex;

} TransientDetector;

/** Transient detection.
  * Holds all transient detectors and buffers used by them.
  */
typedef struct TransientDetection
{
    /** Transient detector. */
    TransientDetector transientDetector;
    /** Delay buffer used by the transient detectors. */
    DelayBuffer delayBuffer;
    /** Subblock energies used by the transient detector. */
    SubblockEnergies subblockEnergies;
} TransientDetection;


typedef struct
{
    int bw_index;                                      /* index of band width */

    /* feature */
    float sp_center[SP_CENTER_NUM];                    /* spectral center*/
    float ltd_stable_rate[STABLE_NUM];                 /* time-domain stable rate*/
    float sfm[SFM_NUM];                                /* spectral flatness*/
    float f_tonality_rate[TONA_NUM];                   /* tonality rate*/
    float frame_sb_energy[BG_ENG_NUM];                 /* energy of sub-band divided non-uniformly*/
    float frames_power[POWER_NUM];                     /* energy of several frames*/
    float pre_spec_low_dif[PRE_SPEC_DIF_NUM];          /* low frequency spectral different*/
    float t_bg_energy;                                 /* time background energy of several frames*/
    float t_bg_energy_sum;                             /* number of time background energy*/
    int   tbg_energy_count;                            /* sum of time background energy of several frames*/
    int   bg_update_count;                             /* time of background update*/
    float frame_energy_smooth;                         /* smoothed energy of several frames*/


    /************************************************************************/
    /* history  parameters                                                  */
    /************************************************************************/
    float smooth_spec_amp[SPEC_AMP_NUM];               /* smoothed spectral amplitude*/
    float sb_bg_energy[BG_ENG_NUM];                    /* sub-band background energy*/
    float pre_snr[PRE_SNR_NUM];                        /* previous time SNR*/
    float lt_snr_org;                                  /* original long time SNR*/
    float lf_snr_smooth;                               /* smoothed lf_snr*/
    float l_silence_snr;                               /* sum of snr's of non active frames*/
    float l_speech_snr;                                /* sum of snr's of active frames*/
    int   l_silence_snr_count;                         /* number of non active frames*/
    int   l_speech_snr_count;                          /* number of active frames*/
    float fg_energy ;                                  /* foreground energy sum  */
    float bg_energy ;                                  /* background energy sum  */
    int   fg_energy_count ;                            /* number of the foreground energy frame */
    int   bg_energy_count ;                            /* number of the background energy frame */
    int   fg_energy_est_start;                         /* flag by that indicate whether if estimate energy*/
    int   speech_flag;                                 /* residual number of hangover 1 */
    int   continuous_noise_num;                        /* time of continuous noise frames*/
    int   continuous_speech_num;                       /* time of continuous speech frames*/
    int   continuous_speech_num2;                      /* time 2 of continuous speech frames*/
    int   frameloop;                                   /* number of frame*/
    float tonality_rate3 ;                             /* tonality rate*/
    float music_background_rate;                       /* music background rate*/
    float lt_noise_sp_center_diff_sum;                 /* different sum of long time noise sp_center*/
    float lt_noise_sp_center_diff_counter;             /* number of the member lt_noise_sp_center_diff_sum*/
    float lt_noise_sp_center0 ;                        /* long time noise sp_center0*/
    float lt_noise_sp_center3 ;                        /* long time noise sp_center3*/
    float lt_bg_highf_eng;                             /* average  of long time high frequency energy*/
    int  update_num_with_snr;                          /* the number of the background update with SNR*/
    int  update_count;
    short  warm_hang_num;                                    /* the number of hangover for warm up*/
    short  vad_flag_for_bk_update;
} T_CldfbVadState;


typedef struct igfscfenc_public_data_struct
{
    int               ptrBitIndex;
    int               bitCount;
    int               prev[64];  /* no more than 64 SCFs for the IGF energy envelope of one block */
    int               prevSave[64];
    int               scfCountLongBlock;
    int               t;
    int               Tsave;
    int               contex_saved;
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
} IGFSCFENC_INSTANCE, *IGFSCFENC_INSTANCE_HANDLE;


typedef struct igf_enc_private_data_struct
{
    IGF_INFO                    igfInfo;
    int                         igfScfQuantized[IGF_MAX_SFB];
    IGFSCFENC_INSTANCE          hIGFSCFArithEnc;


    int                         igfCurrWhiteningLevel[IGF_MAX_TILES];
    int                         igfPrevWhiteningLevel[IGF_MAX_TILES];

    float                       prevSFM_FIR[IGF_MAX_TILES];
    float                       prevSFM_IIR[IGF_MAX_TILES];
    int                         wasTransient;

    UWord8                      igfBitstream[BITBUFSIZE/8];
    short                       igfBitstreamBits;


} IGF_ENC_PRIVATE_DATA, *IGF_ENC_PRIVATE_DATA_HANDLE;

/**********************************************************************/ /**
Definition of public data for this module
**************************************************************************/
typedef struct igf_enc_instance_struct
{
    IGF_ENC_PRIVATE_DATA  igfData;
    int                   infoSamplingRate;
    int                   infoStartFrequency;
    int                   infoStopFrequency;
    int                   infoStartLine;
    int                   infoStopLine;
    int                   infoTotalBitsWritten;
    int                   infoTotalBitsPerFrameWritten;
    int                   infoFrameCount;
    int                   flatteningTrigger;
    float                 spec_be_igf[N_MAX_TCX-IGF_START_MN];
    float                 tns_predictionGain;
} IGF_ENC_INSTANCE, *IGF_ENC_INSTANCE_HANDLE;

/*----------------------------------------------------------------------------------*
 *
 * Main encoder structure
 *
 *----------------------------------------------------------------------------------*/

typedef struct
{

    /*----------------------------------------------------------------------------------*
     * Common parameters
     *----------------------------------------------------------------------------------*/

    short codec_mode;                                   /* Mode1 or Mode2 */
    short last_codec_mode;                              /* previous frame Mode 1 or 2 */
    short last_codec_mode_cng;                          /* previous inactive frame Mode 1 or 2 */
    short mdct_sw_enable;                               /* MDCT switching enable flag */
    short mdct_sw;                                      /* MDCT switching indicator */
    float prev_hi_ener;
    short prev_hi_sparse;
    float clas_sec_old;
    short clas_final_old;
    float last_gain1;
    float last_gain2;

    short nb_bits_tot;                                  /* total number of bits already written */
    Indice *ind_list;                                   /* list of indices */
    short bitstreamformat;                              /* Bitstream format flag (G.192/MIME) */
    short next_ind;                                     /* pointer to the next empty slot in the list of indices */
    short last_ind;                                     /* last written indice */

    int   input_Fs;                                     /* input signal sampling frequency in Hz */
    long  total_brate;                                  /* total bitrate in kbps of the codec */
    long  last_total_brate;                             /* last frame's total bitrate in kbps of the codec */
    long  last_total_brate_cng;                         /* last inactive frame's total bitrate in kbps of the codec */
    short core;                                         /* core (ACELP_CORE, TCX_20_CORE, TCX_10_CORE, HQ_CORE, AMR_WB_CORE) */
    long  core_brate;                                   /* core bitrate */
    long  last_core_brate;                              /* previous frame core bitrate */
    short extl;                                         /* extension layer */
    short last_extl;                                    /* previous extension layer */
    long  extl_brate;                                   /* extension layer bitrate */
    short input_bwidth;                                 /* input signal bandwidth */
    short bwidth;                                       /* encoded bandwidth NB, WB, SWB or FB */
    short max_bwidth;                                   /* maximum encoded bandwidth */
    short last_input_bwidth;                            /* input signal bandwidth in the previous frame */
    short last_bwidth;                                  /* coded bandwidth in the previous frame */
    short last_bwidth_cng;                              /* coded bandwidth in the previous inactive frame */
    short L_frame;                                      /* ACELP core internal frame length */
    short Opt_AMR_WB;                                   /* flag indicating AMR-WB IO mode */
    short Opt_DTX_ON;                                   /* flag indicating DTX operation */
    short cng_type;                                     /* flag indicating LP or CLDFB based SID/CNG */
    short active_cnt;                                   /* counter of active frames */
    short Opt_SC_VBR;                                   /* flag indicating SC-VBR mode */
    short last_Opt_SC_VBR;                              /* flag indicating prev frame's SC-VBR mode */
    short lp_cng_mode2;
    /*----------------------------------------------------------------------------------*
     * ACELP core parameters
     *----------------------------------------------------------------------------------*/

    LPD_state LPDmem;                                  /* ACELP memories */

    short clas;                                         /* current frame clas */
    short last_clas;                                    /* previous frame signal classification */
    float Bin_E[L_FFT];                                 /* per bin energy of two frames */
    float Bin_E_old[L_FFT/2];                           /* per bin energy of old 2nd frames */
    float lsp_old1[M];                                  /* old unquantized LSP vector at the end of the frame at 12k8 */
    float lsf_old1[M];                                  /* old unquantized LSF vector at the end of the frame at 12k8 */
    float lsp_old[M];                                   /* old LSP vector at the end of the frame */
    float lsf_old[M];                                   /* old LSF vector at the end of the frame */
    float lsp_old16k[M];                                /* old LSP vector at the end of the frame @16kHz */
    float lspold_enc[M];                                /* old lsp (immittance spectral pairs) */
    short pstreaklen;                                   /* LSF quantizer */
    float streaklimit;                                  /* LSF quantizer */
    unsigned int offset_scale1[MAX_NO_MODES+1][MAX_NO_SCALES+1]; /* offsets for LSF LVQ structure 1st 8-dim subvector*/
    unsigned int offset_scale2[MAX_NO_MODES+1][MAX_NO_SCALES+1]; /* offsets for LSF LVQ structure 2nd 8-dim subvector*/
    unsigned int offset_scale1_p[MAX_NO_MODES_p+1][MAX_NO_SCALES+1]; /* offsets for LSF LVQ structure, pred. case, 1st 8-dim subvector*/
    unsigned int offset_scale2_p[MAX_NO_MODES_p+1][MAX_NO_SCALES+1]; /* offsets for LSF LVQ structure, pred. case, 2nd 8-dim subvector*/
    short no_scales[MAX_NO_MODES][2];                   /* LSF LVQ structure */
    short no_scales_p[MAX_NO_MODES_p][2];               /* LSF LVQ structure */
    float stab_fac;                                     /* LSF stability factor */
    float mem_decim[2*L_FILT_MAX];                      /* decimation filter memory */
    float mem_deemph;                                   /* deemphasis filter memory */
    float mem_preemph;                                  /* preemphasis filter memory */
    float mem_hp20_in[4];                               /* HP filter memory */
    float old_inp_12k8[L_INP_MEM];                      /* memory of input signal at 12.8kHz */
    float old_wsp[L_WSP_MEM];                           /* old weighted signal vector */
    float old_wsp2[(L_WSP_MEM - L_INTERPOL)/OPL_DECIM]; /* old decimated weighted signal vector */
    float fr_bands1[NB_BANDS];            /* spectrum per critical bands of the previous frame  */
    float fr_bands2[NB_BANDS];            /* spectrum per critical bands 2 frames ago  */
    float mem_wsp;                                      /* weighted signal vector memory */
    float mem_decim2[3];                                /* weighted signal decimation filter memory */
    float mem_syn1[M];                                  /* synthesis filter memory (for core switching and FD BWE) */

    float clip_var[6];                                  /* pitch gain clipping memory */
    float past_qua_en[4];                        /* gain quantization memory (used also in AMR-WB IO mode) */
    float mem_AR[M];                                    /* AR memory of LSF quantizer (past quantized LSFs without mean) */
    float mem_MA[M];                             /* MA memory of LSF quantizer (past quantized residual) (used also in AMR-WB IO mode) */
    short mCb1;                                         /* LSF quantizer - counter of stationary frames after a transition frame */
    short coder_type_raw;                               /* raw coder_type (before UNVOICED is lost) */
    short last_coder_type_raw;                          /* raw last_coder_type (coming from the sigal classification) */
    short last_coder_type;                              /* previous coding type */
    short ini_frame;                                    /* initialization frames counter */
    float old_thres;                                    /* normalized correlation weighting in open-loop pitch */
    float old_corr;                                     /* normalized correlation in previous frame (mean value) */
    short old_pitch;                                    /* previous pitch for open-loop pitch search */
    short delta_pit;                                    /* open-loop pitch extrapolation correction */
    float ee_old;                                       /* previous frame low/high frequency energy ratio */
    short min_band;                                     /* minimum critical band of useful bandwidth */
    short max_band;                                     /* maximum critical band of useful bandwidth */
    short tc_cnt;                                       /* TC frame counter */
    short audio_frame_cnt;                              /* Counter of relative presence of audio frames      */
    float old_dE1;                                      /* Maximum energy increase in previous frame */
    short old_ind_deltaMax;                             /* Index of the sub-subframe of maximum energy in previous frame */
    float old_enr_ssf[2*NB_SSF];                        /* Maxima of energies per sub-subframes of previous frame */
    short spike_hyst;                                   /* Hysteresis to prevent UC after sharp energy spike */
    short music_hysteresis;                             /* Counter of frames after AUDIO coding mode to prevent UC */
    short last_harm_flag_acelp;                         /* harmonicity flag for ACELP @32kbps rate */
    float old_Aq_12_8[M+1];                             /* old Aq[] for core switching */
    float old_Es_pred;                                  /* old Es_pred for core switching */
    short high_stable_cor;

    short seed_tcx;                                     /* AC mode (GSC) - seed for noise fill*/
    short cor_strong_limit;                             /* AC mode (GSC) - Indicator about high spectral correlation per band */
    short GSC_noisy_speech;                             /* AC mode (GSC) - flag to indicate GSC on SWB noisy speech */
    short mem_last_pit_band;                            /* AC mode (GSC) - memory of the last band where pitch contribution was significant */
    float mem_w0_tmp;
    float mem_syn_tmp[M];
    float var_cor_t[VAR_COR_LEN];
    float mid_dyn;                                      /* AC mode (GSC) - signal dynamic                         */
    short noise_lev;                                    /* AC mode (GSC) - noise level                             */
    short past_dyn_dec;                                 /* AC mode (GSC) - Past noise level decision */
    float Last_frame_ener;                              /* AC mode (GSC) - Last frame energy  */
    short pit_exc_hangover;                              /* AC mode (GSC) - Hangover for the time contribution switching*/

    float last_exc_dct_in[L_FRAME];                     /* AC mode (GSC) - previous excitation                       */
    float last_ener;                                    /* AC mode (GSC) - previous energy                           */
    short last_bitallocation_band[6];                   /* AC mode (GSC) - previous bit allocation of each band      */


    short last_act_dec_hang;                            /* Speech/music classifier - decision from the last active frame */
    float past_PS[HIGHEST_FBIN-LOWEST_FBIN];
    float past_ps_diff;
    float past_epsP2;
    short inact_cnt;
    float wdrop;
    float wdlp_0_95_sp;
    short sp_mus_state;
    short past_dec[HANG_LEN-1];                         /* Speech/music classifier - buffer of past binary decisions */
    float past_dlp[HANG_LEN-1];                         /* Speech/music classifier - buffer of past non-binary decisions */
    float last_lsp[M_LSP_SPMUS];
    float last_cor_map_sum;
    float last_non_sta;
    float past_log_enr[NB_BANDS_SPMUS];                 /* Speech/music classifier - last average per-band log energy used for non_staX */
    float gsc_thres[4];                                 /* Speech/music classifier - classification threshold */
    float gsc_lt_diff_etot[MAX_LT];                     /* Speech/music classifier - long-term total energy variation */
    float gsc_mem_etot;                                 /* Speech/music classifier - total energy memory  */
    short gsc_last_music_flag;                          /* Speech/music classifier - last music flag */
    short gsc_nb_thr_1;                                 /* Speech/music classifier - number of consecutives frames of level 1 */
    short gsc_nb_thr_3;                                 /* Speech/music classifier - number of consecutives frames of level 3 */
    float mold_corr;
    float lt_gpitch;
    float mean_avr_dyn;                                 /* Speech/music classifier - long term average dynamic */
    float last_sw_dyn;                                  /* Speech/music classifier - last dynamic              */
    float lt_dec_thres;                                 /* Speech/music classifier - Long term speech/music thresold values */
    float ener_RAT;                                     /* Speech/music classifier - LF/to total energy ratio */

    /* speech/music classifier improvement parameters */
    float old_Bin_E[3*N_OLD_BIN_E];
    float buf_flux[BUF_LEN];
    float buf_pkh[BUF_LEN];
    float buf_epsP_tilt[BUF_LEN];
    float buf_cor_map_sum[BUF_LEN];
    float buf_Ntonal[BUF_LEN];
    float buf_Ntonal2[BUF_LEN];
    float buf_Ntonal_lf[BUF_LEN];
    float buf_dlp[10];
    short onset_cnt;
    float buf_etot[4];
    short attack_hangover;
    float dec_mov;
    float dec_mov1;
    float mov_log_max_spl;
    float old_lt_diff[2];
    short UV_cnt1;
    float LT_UV_cnt1;

    short lt_music_hangover;
    float tonality2_buf[HANG_LEN_INIT];
    float tonality3_buf[HANG_LEN_INIT];
    float LPCErr_buf[HANG_LEN_INIT];
    short lt_music_state;
    short lt_speech_state;
    short lt_speech_hangover;
    float lpe_buf[HANG_LEN_INIT];
    float voicing_buf[HANG_LEN_INIT];
    short gsc_hangover;
    float sparse_buf[HANG_LEN_INIT];
    float hf_spar_buf[HANG_LEN_INIT];
    float LT_sparse;
    short gsc_cnt;

    short Last_pulse_pos;                               /* FEC - last position of the first glottal pulse in the frame */
    float lsfoldbfi0[M];                                /* FEC - LSF vector of the previous frame */
    float lsfoldbfi1[M];                                /* FEC - LSF vector of the past previous frame */
    float lsf_adaptive_mean[M];                         /* FEC - adaptive mean LSF vector for FEC */
    short next_force_safety_net;                        /* FEC - flag to force safety net in next frame */
    float old_S[L_FFT/2];                               /* Tonal detector - prev. log-energy spectrum with subtracted floor */
    float cor_map[L_FFT / 2];                           /* Tonal detector - LT correlation map */
    float noise_char;                                   /* Tonal detector - LT noise character */
    float ave_enr2[NB_BANDS];                           /* Tonal detector - LT average E per crit. band (for non_sta2) */
    float act_pred;                                     /* Tonal detector - prediction of speech activity from 0 to 1 (0-inactive, 1-active) */
    float multi_harm_limit;                             /* Tonal detector - adaptive threshold */
    float enrO[NB_BANDS];                               /* Noise estimator - previous energy per critical band */
    float bckr[NB_BANDS];                               /* Noise estimator - background noise estimation per critical band */
    float ave_enr[NB_BANDS];                            /* Noise estimator - long-term average energy per critical band */
    short pitO;                                         /* Noise estimator - previous open-loop pitch value */
    short aEn;                                          /* Noise estimator - noise estimator adaptation flag */
    float totalNoise;                                   /* Noise estimator - total noise energy */
    short first_noise_updt;                             /* Noise estimator - flag used to determine if the first noise update frame */
    short harm_cor_cnt;                                 /* Noise estimator - 1st memory counter of harm or correlation frame */
    short bg_cnt;                                       /* Noise estimator - pause length counter */
    float prim_act_quick;                               /* Noise estimator - primary activity quick */
    float prim_act_slow;                                /* Noise estimator - primary activity slow  */
    float prim_act;                                     /* Noise estimator - primary activity slow rise quick fall */
    float prim_act_quick_he;                            /* Noise estimator - primary activity quick */
    float prim_act_slow_he;                             /* Noise estimator - primary activity slow  */
    float prim_act_he;                                  /* Noise estimator - primary activity slow rise quick fall */
    float Etot_l;                                       /* Noise estimator - Track energy from below  */
    float Etot_h;                                       /* Noise estimator - Track energy from above  */
    float Etot_l_lp;                                    /* Noise estimator - Smoothed low energy      */
    float Etot_last;                                    /* Noise estimator - Energy of last frame     */
    float Etot_lp;                                      /* Noise estimator - Filtered input energy   */
    float lt_tn_track;
    float lt_tn_dist;
    float lt_Ellp_dist;
    float lt_haco_ev;
    short low_tn_track_cnt;
    float epsP_0_2_lp;
    float epsP_0_2_ad_lp;
    float epsP_2_16_lp;
    float epsP_2_16_lp2;
    float epsP_2_16_dlp_lp;
    float epsP_2_16_dlp_lp2;
    float lt_aEn_zero;

    float Etot_st_est;                                  /* Noise estimation - short term estimate of E{ Etot } */
    float Etot_sq_st_est;                               /* Noise estimation - short term estimate of E{ Etot^2 } */


    short nb_active_frames;
    short hangover_cnt;
    float lp_speech;
    float Etot_v_h2;
    float sign_dyn_lp;
    short nb_active_frames_he;
    short hangover_cnt_he;
    short nb_active_frames_HE_SAD;
    long  vad_flag_reg_H;
    long  vad_flag_reg_L;
    long  vad_prim_reg;

    short vad_flag_cnt_50;
    short vad_prim_cnt_16;

    short hangover_cnt_dtx;
    short hangover_cnt_music;

    float bcg_flux;
    short soft_hangover;
    short voiced_burst;
    short bcg_flux_init;
    float voicing_old;
    short nb_active_frames_he1;
    short hangover_cnt_he1;

    float bckr_tilt_lt;

    short var_SID_rate_flag;                            /* CNG and DTX - flag for variable SID rate */
    float lp_ener;                                      /* CNG and DTX - low-pass filtered energy for CNG */
    short cng_seed;                                     /* CNG and DTX - seed for white noise random generator */
    float lspCNG[M];                                    /* CNG and DTX - LP filtered ISPs */
    short first_CNG;                                    /* CNG and DTX - first CNG frame flag */
    float lp_noise;                                     /* CNG and DTX - LP filterend total noise estimation */
    short cnt_SID;                                      /* CNG and DTX - counter of SID update for the interop. mode or dtx, if enabled */
    short max_SID;                                      /* CNG and DTX - max allowed number of CNG FRAME_NO_DATA frames */
    short interval_SID;                                 /* CNG and DTX - interval of SID update, default 8 */
    short old_enr_index;                                /* CNG and DTX - index of last encoded CNG energy */
    float Enew;                                         /* CNG and DTX - CNG target residual energy */
    short VarDTX_cnt_voiced;                            /* CNG and DTX - counter for variable DTX activation (speech) */
    float lt_ener_voiced;                               /* CNG and DTX - long-term energy of signal (measured on voiced parts) */
    short VarDTX_cnt_noise;                             /* CNG and DTX - counter for variable DTX activation (noise) */
    float lt_ener_noise;                                /* CNG and DTX - long-term energy of background noise */
    float lt_ener_last_SID;                             /* CNG and DTX - long-term energy of last SID frame */
    short cng_hist_size;                                /* CNG and DTX - size of CNG history buffer for averaging, <0,DTX_HIST_SIZE> */
    short cng_hist_ptr;                                 /* CNG and DTX - pointer for averaging buffers */
    float cng_lsp_hist[DTX_HIST_SIZE*M];                /* CNG and DTX - old LSP buffer for averaging */
    float cng_ener_hist[DTX_HIST_SIZE];                 /* CNG and DTX - log energy buffer for averaging */
    short cng_cnt;                                      /* CNG and DTX - counter of CNG frames for averaging */
    short cng_ener_seed;                                /* CNG and DTX - seed for random generator for variation of excitation energy */
    short cng_ener_seed1;
    float frame_ener;
    float lp_sp_enr;
    short last_allow_cn_step;
    short ho_hist_size;                                 /* CNG and DTX - size of DTX hangover history buffer for averaging, <0,HO_HIST_SIZE> */
    short ho_hist_ptr;                                  /* CNG and DTX - pointer for averaging buffers */
    long  ho_sid_bw;                                    /* CNG and DTX - SID bandwidth flags */
    float ho_lsp_hist[HO_HIST_SIZE*M];                  /* CNG and DTX - old LSP buffer for averaging */
    float ho_ener_hist[HO_HIST_SIZE];                   /* CNG and DTX - energy buffer for averaging */
    float ho_env_hist[HO_HIST_SIZE*NUM_ENV_CNG];
    short act_cnt;                                      /* CNG and DTX - counter of active frames */
    short ho_circ_size;                                 /* CNG and DTX - size of DTX hangover history buffer for averaging, <0,HO_HIST_SIZE> */
    short ho_circ_ptr;                                  /* CNG and DTX - pointer for averaging buffers */
    float ho_lsp_circ[HO_HIST_SIZE*M];                  /* CNG and DTX - old LSP buffer for averaging */
    float ho_ener_circ[HO_HIST_SIZE];                   /* CNG and DTX - energy buffer for averaging */
    float ho_env_circ[HO_HIST_SIZE*NUM_ENV_CNG];
    short burst_ho_cnt;                                 /* CNG and DTX - counter of hangover frames at end of active burst */
    short cng_buf_cnt;                                  /* CNG and DTX - Counter of buffered CNG parameters */
    float cng_exc2_buf[HO_HIST_SIZE*L_FFT];             /* CNG and DTX - exc2 buffer for storing */
    long  cng_brate_buf[HO_HIST_SIZE];                  /* CNG and DTX - buffer for storing last_active_brate */

    short CNG_mode;                                     /* CNG and DTX - mode for DTX configuration */
    long last_active_brate;                             /* CNG and DTX - last active frame bitrate used for CNG_mode control */
    short ho_16k_lsp[HO_HIST_SIZE];                     /* CNG and DTX - 16k LSPs flags */
    short last_CNG_L_frame;                             /* CNG and DTX - last CNG frame length */
    short act_cnt2;                                     /* CNG and DTX - counter of active frames for CNG_mode switching */
    float ho_lsp_circ2[HO_HIST_SIZE*M];                 /* CNG and DTX - second buffer of LSPs */
    short num_ho;                                       /* CNG and DTX - number of selected hangover frames */
    short hangover_terminate_flag;                      /* CNG and DTX - flag indicating whether to early terminate DTX hangover */
    float old_env[NUM_ENV_CNG];
    float lp_env[NUM_ENV_CNG];
    float cng_res_env[NUM_ENV_CNG*HO_HIST_SIZE];
    float exc_mem[24];
    float exc_mem1[30];
    float exc_mem2[30];

    float dispMem[8];                                   /* Noise enhancer - phase dispersion algorithm memory */

    short uv_count;                                     /* Stationary noise UV modification - unvoiced counter */
    short act_count;                                    /* Stationary noise UV modification - activation counter */
    float ge_sm;                                        /* Stationary noise UV modification - smoothed excitation gain */
    float lspold_s[M];                                  /* Stationary noise UV modification - old LSP vector */
    short noimix_seed;                                  /* Stationary noise UV modification - mixture seed */
    float min_alpha;                                    /* Stationary noise UV modification - minimum alpha */
    float exc_pe;                                       /* Stationary noise UV modification - memory of the preemphasis filter */

    short last_L_frame;                                 /* ACELP@16kHz - last L_frame value */
    float mem_decim16k[2*L_FILT_MAX];                   /* ACELP@16kHz - decimation filter memory @16kHz */
    float mem_preemph16k;                               /* ACELP@16kHz - preemphasis filter memory @16kHz */
    float old_inp_16k[L_INP_MEM];                       /* ACELP@16kHz - memory of input signal @16 kHz */
    float mem_deemp_preQ;                               /* ACELP@16kHz - prequantizer deemhasis memory */
    float mem_preemp_preQ;                              /* ACELP@16kHz - prequantizer preemhasis memory */
    short last_nq_preQ;                                 /* ACELP@16kHz - AVQ subquantizer number of the last sub-band of the last subframe  */
    short use_acelp_preq;                               /* ACELP@16kHz - flag of prequantizer usage */

    short bpf_off;                                      /* Bass post-filter - do not use BPF when this flag is set to 1 */
    float old_pitch_buf[2*NB_SUBFR16k];                 /* Bass post-filter - buffer of old subframe pitch values */
    float pst_mem_deemp_err;                            /* Bass post-filter - filter memory of noise LP filter */
    float pst_lp_ener;                                  /* Bass post-filter - long-term energy */
    float lps;
    float lpm;

    /* stable short pitch detection */
    float voicing0_sm;
    float voicing_sm;
    float LF_EnergyRatio_sm;
    short predecision_flag;
    float diff_sm;
    float energy_sm;


    /*----------------------------------------------------------------------------------*
     * HF WB BWE for AMR-WB IO mode at 23.85 kbps
     *----------------------------------------------------------------------------------*/

    float gain_alpha;
    float mem_hf2_enc[L_FIR-1];
    float mem_hp400_enc[4];
    float mem_hf_enc[L_FIR-1];
    float mem_syn_hf_enc[M];
    short seed2_enc;

    /*----------------------------------------------------------------------------------*
     * CLDFB analysis
     *----------------------------------------------------------------------------------*/

    HANDLE_CLDFB_FILTER_BANK cldfbAnaEnc;          /* main analysis filter bank handle */
    HANDLE_CLDFB_FILTER_BANK cldfbSynTd;


    /*----------------------------------------------------------------------------------*
     * FD CNG handle
     *----------------------------------------------------------------------------------*/

    HANDLE_FD_CNG_ENC hFdCngEnc;
    short fd_cng_reset_flag;
    float last_totalNoise;
    float totalNoise_increase_hist[TOTALNOISE_HIST_SIZE];
    short totalNoise_increase_len;

    /*----------------------------------------------------------------------------------*
     * SC-VBR parameters
     *----------------------------------------------------------------------------------*/

    float vadsnr;
    float vadnoise;

    /* NELP variables */
    float shape1_filt_mem[20];
    float shape2_filt_mem[20];
    float shape3_filt_mem[20];
    float txlpf1_filt1_mem[20];
    float txlpf1_filt2_mem[20];
    float txhpf1_filt1_mem[20];
    float txhpf1_filt2_mem[20];
    float bp1_filt_mem_wb[8];
    float nelp_lp_fit_mem[NELP_LP_ORDER*2];
    float bp1_filt_mem_nb[14];
    short nelp_enc_seed;
    float nelp_gain_mem;
    short last_nelp_mode;
    short nelp_mode;

    /* PPP variables */
    short pppcountE;
    short bump_up;
    short last_ppp_mode;
    short last_last_ppp_mode;
    short ppp_mode;
    float prev_ppp_gain_pit;
    float prev_tilt_code;

    /* voiced encoder variables */
    int   firstTime_voicedenc;

    /* DTFS variables */
    float dtfs_enc_a[MAXLAG_WI];
    float dtfs_enc_b[MAXLAG_WI];
    int   dtfs_enc_lag;
    int   dtfs_enc_nH;
    int   dtfs_enc_nH_4kHz;
    float dtfs_enc_upper_cut_off_freq_of_interest;
    float dtfs_enc_upper_cut_off_freq;

    float prev_cw_en;
    float ph_offset_E;
    float lastLgainE;                                   /* Previous gain value for the low band */
    float lastHgainE;                                   /* Previous gain value for the high band */
    float lasterbE[NUM_ERB_WB];                         /* Previous Amplitude spectrum (ERB) */

    short mode_QQF;
    short rate_control;
    float SNR_THLD;
    short Q_to_F;
    short pattern_m;
    short patterncount;
    short Last_Resort;
    short numactive;                                    /* keep the count of the frames inside current 600 frame block */
    float sum_of_rates;                                 /* sum of the rates of past 600 active frames */
    float global_avr_rate;                              /* global rate up to current time. recorded a (rate in kbps) * 6000 */
    int global_frame_cnt;                               /* 600 active frame block count. Used to update the global rate */
    short set_ppp_generic;
    short avoid_HQ_VBR_NB;

    /*----------------------------------------------------------------------------------*
     * HQ core parameters
     *----------------------------------------------------------------------------------*/

    float input_buff[L_FRAME48k+L_FRAME48k+NS2SA(48000, DELAY_FIR_RESAMPL_NS)];
    float * input;
    float * old_input_signal;

    float old_hpfilt_in;
    float old_hpfilt_out;
    float EnergyLT;
    float Energy_Old;
    short TransientHangOver;
    float old_out[L_FRAME32k];             /* buffer for OLA; at the encoder, the maximum length is L_FRAME32k (corresponds to maximum internal L_frame length) */
    short last_core;
    short hq_generic_speech_class;
    short mode_count;                      /* HQ_HARMONIC mode count                     */
    short mode_count1;                     /* HQ_NORMAL mode count                       */
    short Nb_ACELP_frames;
    short prev_Npeaks;                     /* number of peaks in previous frame */
    short prev_peaks[HVQ_MAX_PEAKS];       /* indices of the peaks in previous frame */
    short hvq_hangover;
    short prev_hqswb_clas;
    short prev_SWB_peak_pos[SPT_SHORTEN_SBNUM];

    /* speech/music classification */
    short last_vad_spa;
    short lt_old_mode[3];
    float lt_voicing;
    float lt_corr;
    float lt_tonality;
    short lt_corr_pitch[3];
    short lt_hangover;
    float lowrate_pitchGain;

    short prev_frm_index[NB_SWB_SUBBANDS_HAR_SEARCH_SB];
    short prev_frm_hfe2;
    short prev_stab_hfe2;
    float prev_ni_ratio;
    float prev_En_sb[NB_SWB_SUBBANDS];
    short last_bitalloc_max_band[2];
    float last_ni_gain[BANDS_MAX];
    float last_env[BANDS_MAX];
    short last_max_pos_pulse;

    /* PVQ range coder state */
    unsigned int rc_low;
    unsigned int rc_range;
    short rc_cache;
    short rc_carry;
    short rc_carry_count;
    short rc_num_bits;
    short rc_tot_bits;
    short rc_offset;
    /*----------------------------------------------------------------------------------*
     * TBE parameters
     *----------------------------------------------------------------------------------*/

    float old_speech_shb[L_LOOK_16k + L_SUBFR16k];          /* Buffer memories */
    float old_speech_wb[(L_LOOK_12k8 + L_SUBFR) * 5/16];    /* Buffer memories */
    float old_input_fhb[NS2SA(48000, ACELP_LOOK_NS + DELAY_FD_BWE_ENC_12k8_NS + DELAY_FIR_RESAMPL_NS) - L_FRAME48k/2];
    float prev_lsp_shb[LPC_SHB_ORDER];
    float state_ana_filt_shb[ 2 * ALLPASSSECTIONS_STEEP + 1 ]; /* states for the analysis filters */
    float cldfbHBLT;
    /* states for the filters used in generating SHB excitation from WB excitation*/
    float mem_csfilt[2];

    /* states for the filters used in generating SHB signal from SHB excitation*/
    float state_syn_shbexc[L_SHB_LAHEAD];
    float state_lpc_syn[LPC_SHB_ORDER];
    float old_bwe_exc[PIT16k_MAX * 2];                  /* old excitation */
    short bwe_seed[2];
    float bwe_non_lin_prev_scale;
    float old_bwe_exc_extended[NL_BUFF_OFFSET];
    float syn_overlap[L_SHB_LAHEAD];                /* overlap buffer used to Adjust SHB Frame Gain*/
    float decim_state1[(2*ALLPASSSECTIONS_STEEP+1)];
    float decim_state2[(2*ALLPASSSECTIONS_STEEP+1)];
    float mem_genSHBexc_filt_down_wb2[(2*ALLPASSSECTIONS_STEEP+1)];
    float mem_genSHBexc_filt_down_wb3[(2*ALLPASSSECTIONS_STEEP+1)];
    float mem_genSHBexc_filt_down_shb[(2*ALLPASSSECTIONS_STEEP+1)];

    float elliptic_bpf_2_48k_mem[4][4];
    float prev_fb_energy;
    float prev_gainFr_SHB;
    float lsp_shb_slow_interpl[LPC_SHB_ORDER];
    float lsp_shb_fast_interpl[LPC_SHB_ORDER];
    float shb_inv_filt_mem[LPC_SHB_ORDER];
    float lsp_shb_spacing[3];
    float prev_swb_GainShape;
    short prev_frGainAtten;

    short spectral_tilt_reset;
    short consec_inactive;
    float ra_deltasum;
    short trigger_SID;
    float running_avg;
    float snr_sum_vad;

    float prev_wb_GainShape;
    float swb_lsp_prev_interp[LPC_SHB_ORDER];
    float fb_state_lpc_syn[LPC_SHB_ORDER];
    float fb_tbe_demph;
    float tilt_mem;


    short prev_coder_type;
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
    float dec_2_over_3_mem[12];
    float dec_2_over_3_mem_lp[6];


    /*----------------------------------------------------------------------------------*
     * SWB BWE parameters
     *----------------------------------------------------------------------------------*/

    float new_input_hp[NS2SA(16000, ACELP_LOOK_NS + DELAY_FD_BWE_ENC_NS + DELAY_FIR_RESAMPL_NS - DELAY_CLDFB_NS)];
    float old_input[NS2SA(48000, DELAY_FD_BWE_ENC_NS + DELAY_FIR_RESAMPL_NS)];
    float old_input_wb[NS2SA(16000, DELAY_FD_BWE_ENC_NS)];
    float old_input_lp[NS2SA(16000, ACELP_LOOK_NS + DELAY_FD_BWE_ENC_NS)];
    float old_syn_12k8_16k[NS2SA(16000, DELAY_FD_BWE_ENC_NS)];
    float old_fdbwe_speech[L_FRAME48k];
    float mem_deemph_old_syn;
    short prev_mode;
    float old_wtda_swb[L_FRAME48k];
    short prev_L_swb_norm1;
    float prev_global_gain;
    short modeCount;
    float EnergyLF;
    float tbe_demph;
    float tbe_premph;
    float mem_stp_swb[LPC_SHB_ORDER];
    float *ptr_mem_stp_swb;
    float gain_prec_swb;
    float mem_zero_swb[LPC_SHB_ORDER];

    /*----------------------------------------------------------------------------------*
     * WB, SWB and FB bandwidth detector
     *----------------------------------------------------------------------------------*/

    float lt_mean_NB;
    float lt_mean_WB;
    float lt_mean_SWB;
    short count_WB;
    short count_SWB;
    short count_FB;

    /*----------------------------------------------------------------------------------*
     * SWB DTX/CNG parameters
     *----------------------------------------------------------------------------------*/

    short last_vad;
    float last_wb_cng_ener;
    float last_shb_cng_ener;
    float mov_wb_cng_ener;
    float mov_shb_cng_ener;
    short shb_cng_ini_cnt;
    short last_SID_bwidth;
    short shb_NO_DATA_cnt;

    /*----------------------------------------------------------------------------------*
     * Channel-aware mode
     *----------------------------------------------------------------------------------*/

    short rf_mode;                    /* flag to signal the RF mode */
    short rf_mode_last;
    short Opt_RF_ON;
    short rf_frame_type;

    short rf_target_bits_write;
    short rf_fec_offset;
    short rf_targetbits_buff[MAX_RF_FEC_OFFSET];
    int rf_indx_frametype[MAX_RF_FEC_OFFSET];

    float rf_mem_w0;
    float rf_clip_var[6];
    float rf_tilt_code;
    float rf_mem_syn2[M];
    float rf_dispMem[8];
    float rf_gc_threshold;

    short rf_target_bits;
    float rf_tilt_buf[NB_SUBFR16k];

    short rf_indx_lsf[MAX_RF_FEC_OFFSET][3];
    int rf_indx_pitch[MAX_RF_FEC_OFFSET][NB_SUBFR16k];
    int rf_indx_fcb[MAX_RF_FEC_OFFSET][NB_SUBFR16k];
    int rf_indx_gain[MAX_RF_FEC_OFFSET][NB_SUBFR16k];
    int rf_indx_EsPred[MAX_RF_FEC_OFFSET];
    int rf_indx_ltfMode[MAX_RF_FEC_OFFSET][NB_SUBFR16k];

    short rf_indx_nelp_fid[MAX_RF_FEC_OFFSET];
    short rf_indx_nelp_iG1[MAX_RF_FEC_OFFSET];
    short rf_indx_nelp_iG2[MAX_RF_FEC_OFFSET][2];

    short rf_indx_tbeGainFr[MAX_RF_FEC_OFFSET];

    int rf_tcxltp_pitch_int_past;
    int rf_last_tns_active;
    int rf_second_last_tns_active;
    int rf_second_last_core;
    short rf_clas[MAX_RF_FEC_OFFSET];
    int rf_gain_tcx[MAX_RF_FEC_OFFSET];
    int rf_tcxltp_param[MAX_RF_FEC_OFFSET];

    short rf_fec_indicator;
    short RF_bwe_gainFr_ind;

    /*----------------------------------------------------------------------------------*
     * Local synthesis parameters
     *----------------------------------------------------------------------------------*/


    /*----------------------------------------------------------------------------------*
     *
     * Mode2
     *
     *----------------------------------------------------------------------------------*/

    int frame_size_index;         /*0-FRAME_SIZE_NB-1: index determining the frame size*/
    int bits_frame_nominal;       /*avg bits per frame on active frame*/
    int bits_frame;               /*bits per frame overall */
    int bits_frame_core;          /*bits per frame for the core*/
    int narrowBand;

    /*ACELP config*/
    ACELP_config acelp_cfg;  /*configuration set for each frame*/

    ACELP_config acelp_cfg_rf; /* configuration for RF frame */

    /*TCX config*/
    TCX_config tcx_cfg;
    int L_frameTCX;

    /* cod_main.c */
    float mem_preemph_enc;               /* speech preemph filter memory (at encoder-sampling-rate) */

    /* Signal Buffers and Pointers at encoder-sampling-rate */
    float *speech_enc;
    float *speech_enc_pe;
    float *new_speech_enc;
    float *new_speech_enc_pe;
    float buf_speech_enc[L_PAST_MAX_32k+L_FRAME32k+L_NEXT_MAX_32k];
    float buf_speech_enc_pe[L_PAST_MAX_32k+L_FRAME32k+L_NEXT_MAX_32k];
    float *wspeech_enc;
    float buf_wspeech_enc[L_FRAME16k+L_SUBFR+L_FRAME16k+L_NEXT_MAX_16k];
    float *synth;
    float buf_synth[OLD_SYNTH_SIZE_ENC+L_FRAME32k]; /*can be reduced to PIT_MAX_MAX+L_FRAME_MAX if no rate switching*/
    float *speech_TCX;
    float *new_speech_TCX;

    /* Core Signal Analysis Outputs */
    float *spectrum[2]; /* MDCT output for a short block */
    float spectrum_long[N_MAX]; /* MDCT output for a long block. Points to spectrum */
    float noiseTiltFactor; /* compensation for LPC tilt in noise filling  */
    int noiseLevelMemory;  /* counter of consecutive low TCX noise levels */
    STnsData tnsData[2];
    int fUseTns[2];

    int enableTcxLpc; /* global toggle for the TCX LPC quantizer */
    int envWeighted;  /* are is{p,f}_old_q[] weighted or not? */

    int acelpEnabled; /* Flag indicating if ACELP can be used */
    int tcx10Enabled; /* Flag indicating if TCX 10 can be used */
    int tcx20Enabled; /* Flag indicating if TCX 20 can be used */

    short tcxMode; /* Chosen TCX mode for this frame */

    float mem_wsp_enc;                       /* wsp vector memory */

    int  nb_bits_header_ace;          /* number of bits for the header */
    int  nb_bits_header_tcx;          /* number of bits for the header */

    /* restrict the possible in EVS: 0 base 10 = d.c.b.a base 2*/
    /* a = 0/1 : ACELP on/off*/
    /* b = 0/1 : TCX20 on/off*/
    /* c = 0/1 : TCX40 on/off*/
    /* d = 0/1 : TCX80 on/off*/
    short restrictedMode;




    /* Framing */
    int nb_subfr;

    float preemph_fac;   /*Preemphasis factor*/
    float gamma;

    TransientDetection transientDetection;
    int transient_info[3];

    int   acelpFramesCount;
    float prevTempFlatness;

    float prevEnergyHF;
    float currEnergyHF;
    float currEnergyLookAhead;

    int lpcQuantization;

    int encoderLookahead_enc;
    int encoderPastSamples_enc;
    int encoderLookahead_FB;

    /* pitch_ol for adaptive lag window */
    int old_pitch_la;                  /* past open loop pitch lag from look-ahead  */
    float old_voicing_la;              /* past open loop pitch gain from look-ahead */


    int tcxonly;

    short flag_noisy_speech_snr;          /* encoder detector for noisy speech */

    int fscale;
    int sr_core;
    int acelp_autocorr;                   /* Optimize acelp in 0 covariance or 1 correlation domain */

    int pit_min;
    int pit_fr1;
    int pit_fr1b;
    int pit_fr2;
    int pit_max;
    int pit_res_max;

    /* for FAC */
    int L_frame_past;

    int memQuantZeros[L_FRAME_PLUS];      /* Quantization deadzone flags */

    /*Adaptive BPF*/
    int bpf_gain_param;
    float mem_bpf[2*L_FILT16k];
    float mem_error_bpf[2*L_FILT16k];

    /* TCX-LTP */
    int tcxltp;
    int tcxltp_pitch_int;
    int tcxltp_pitch_fr;
    float tcxltp_gain;
    int tcxltp_pitch_int_past;
    int tcxltp_pitch_fr_past;
    float tcxltp_gain_past;
    float tcxltp_norm_corr_past;
    float buf_speech_ltp[L_PAST_MAX_32k+L_FRAME32k+L_NEXT_MAX_32k];
    float *speech_ltp;
    float *new_speech_ltp;
    short tcxltp_filt_idx;
    int tcxltp_bits;
    int tcxltp_param[LTPSIZE];

    float measuredBwRatio; /* measured bw; used for TCX noise-filling */
    int nmStartLine; /* Starting line for the noise measurement */

    short glr;
    short glr_idx[2];
    float mean_gc[2];
    float prev_lsf4_mean;
    short glr_reset;
    int last_sr_core;
    float last_stab_fac;

    /*for rate switching*/
    short rate_switching_reset; /*Rate switching flag requiring a reset of memories at least partially */
    short rate_switching_reset_16kHz;

    short enablePlcWaveadjust;
    short Tonal_SideInfo;

    short seed_acelp;

    short tcx_lpc_shaped_ari;

    PLC_ENC_EVS plcExt;
    T_CldfbVadState vad_st;

    IGF_ENC_INSTANCE hIGFEnc;
    short igf;

    short tec_tfa;
    TEMPORAL_ENVELOPE_CODING_ENCODER tecEnc;
    short tec_flag;
    short tfa_flag;
    float tfa_enr[N_TEC_TFA_SUBFR];

    short vbr_generic_ho;
    float last_7k2_coder_type;

    short sharpFlag;


    short Local_VAD;

} Encoder_State;


/*---------------------------------------------------------------*
 * Encoder Static RAM                                            *
 *---------------------------------------------------------------*/

typedef struct GainItemStr
{
    float nmrValue;
    short gainIndex;
} GainItem;

/* Structure for storing correlations between ACELP codebook components and target */
typedef struct
{
    float xx;   /* energy of target x */
    float y1y1; /* energy of adaptive cbk contribution y1 */
    float y2y2; /* energy of fixed cbk contribution y2 */
    float xy1;  /* correlation of x and y1 */
    float xy2;  /* correlation of x and y2 */
    float y1y2; /* correlation of y1 and y2 */
} ACELP_CbkCorr;

struct PLC_ENC_EVS;
typedef struct PLC_ENC_EVS * HANDLE_PLC_ENC_EVS;


#endif
