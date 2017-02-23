/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#ifndef ROM_COM_H
#define ROM_COM_H

#include <stdio.h>
#include "options.h"
#include "stat_enc.h"
#include "stat_dec.h"
#include "stl.h"
#include "basop_util.h"

typedef struct
{
    int fin;       /* input frequency                   */
    int fout;      /* output frequency                  */
    short fac_num;   /* numerator of resampling factor    */
    const float *filter;   /* resampling filter coefficients    */
    short filt_len;  /* number of filter coeff.           */
    unsigned short flags;     /* flags from config. table          */
} Resampling_cfg;

typedef struct
{
    short bands;
    short bw;
    const short *band_width;
    Word32 L_qint;
    Word16 eref_fx;
    Word16 bit_alloc_weight_fx;
    short gqlevs;
    short Ngq;
    short p2a_bands;
    float p2a_th;
    float pd_thresh;
    float ld_slope;
    float ni_coef;
    float ni_pd_th;
} Xcore_Config;

/*-----------------------------------------------------------------*
 * Table of bitrates
 *-----------------------------------------------------------------*/

extern const long brate_tbl[SIZE_BRATE_TBL];
extern const long acelp_sig_tbl[MAX_ACELP_SIG];

/*-----------------------------------------------------------------*
 * Bit-allocation tables
 *-----------------------------------------------------------------*/

extern const short LSF_bits_tbl[];                      /* Bit allocation table for end-frame ISF quantizer */
extern const short mid_LSF_bits_tbl[];                  /* Bit allocation table for mid-frame ISF quantizer */
extern const short Es_pred_bits_tbl[];                  /* Bit allocation table for scaled innovation energy prediction */
extern const short gain_bits_tbl[];                     /* Bit allocation table for gain quantizer */
extern const short ACB_bits_tbl[];                      /* Bit allocation table for adaptive codebook (pitch) */
extern const short FCB_bits_tbl[];                      /* Bit allocation table for algebraic (fixed) codebook (innovation) */
extern const short reserved_bits_tbl[];                 /* Bit allocation table for reseved bits */

extern const short ACB_bits_16kHz_tbl[];                /* Bit allocation table for adaptive codebook (pitch) @16kHz */
extern const short FCB_bits_16kHz_tbl[];                /* Bit allocation table for algebraic (fixed) codebook (innovation) @16kHz */
extern const short gain_bits_16kHz_tbl[];               /* Bit allocation table for gain quantizer @16kHz */
extern const short AVQ_bits_16kHz_tbl[];                /* Bit allocation table for AVQ bits @16kHz ACELP, active segments */

extern const long unsigned pulsestostates[17][9];       /* Number of states for any combination of pulses in any combination of vector length */

extern const unsigned char ACELP_NRG_MODE[RATE_MODE_MAX][BANDWIDTH_MODE_MAX][ACELP_MODE_MAX+RF_MODE_MAX];
extern const unsigned char ACELP_NRG_BITS[3];

extern const unsigned char ACELP_LTP_MODE[RATE_MODE_MAX][BANDWIDTH_MODE_MAX][ACELP_MODE_MAX+RF_MODE_MAX];
extern const unsigned char ACELP_LTP_BITS[RATE_MODE_MAX][BANDWIDTH_MODE_MAX][ACELP_MODE_MAX+RF_MODE_MAX];
extern const unsigned char ACELP_LTP_BITS_SFR[8+RF_MODE_MAX][5];

extern const unsigned char ACELP_LTF_MODE[RATE_MODE_MAX][BANDWIDTH_MODE_MAX][ACELP_MODE_MAX+RF_MODE_MAX];
extern const unsigned char ACELP_LTF_BITS[4];

extern const unsigned char ACELP_GAINS_MODE[RATE_MODE_MAX][BANDWIDTH_MODE_MAX][ACELP_MODE_MAX+RF_MODE_MAX];
extern const unsigned char ACELP_GAINS_BITS[10];

extern const unsigned char ACELP_BPF_MODE[RATE_MODE_MAX][BANDWIDTH_MODE_MAX][ACELP_MODE_MAX+RF_MODE_MAX];
extern const unsigned char ACELP_BPF_BITS[3];


/*----------------------------------------------------------------------------------*
 * Pre-processing
 *----------------------------------------------------------------------------------*/

extern const float inv_tbl[];                           /* Table of 1/x values */

extern const Resampling_cfg resampling_cfg_tbl[];       /* table of resampling configurations */
extern const FrameSizeParams FrameSizeConfig[FRAME_SIZE_NB];

extern const float h_high[];                            /* HP filter for filtering random part of excitation in FEC */
extern const float crit_bands[];                        /* Table of critical bands */
extern const int   pow2[];                              /* Table with power of 2 values */
extern const float sincos_t[];                          /* FFT - sinus and cosinus tables */
extern const float sincos_t_ext[];
extern const float sincos_t_rad3[];
extern const short fft256_read_indexes[];               /* FFT */
extern const float inter4_2[];                          /* 1/4 resolution interpolation filter */
extern const float LP_assym_window[];                   /* Assymetric window for LP analysis @12.8kHz */
extern const float LP_assym_window_16k[];               /* Assymetric window for LP analysis @16kHz   */
extern const float hamcos_window[];                     /* Hamming-Cosinus window */
extern const float grid50[];                            /* Table of grid points for evaluating Chebyshev polynomials */
extern const float grid40[];                            /* Table of grid points for evaluating Chebyshev polynomials */
extern const float grid100[];                           /* Table of 100 grid points for evaluating Chebyshev polynomials */

extern const float wind_sss[LEN_WIN_SSS];               /* window for modify_sf ana */
extern const float filter5_39s320_120[];                /* LP FIR filter for 8kHz signal resampling */

extern const float lag_window_8k[17];
extern const float lag_window_12k8[][17];
extern const float lag_window_16k[][17];
extern const float lag_window_25k6[][17];
extern const float lag_window_32k[][17];
extern const float lag_window_48k[17];

extern const float interpol_frac_12k8[];                /* LPC interpolation coefficients */
extern const float interpol_isp_amr_wb[];               /* LPC interpolation coefficients for AMR-WB interoperable mode */
extern const float interpol_frac_16k[];                 /* LPC interpolation coefficients @ 16kHz */
extern const float interpol_frac_mid[];                 /* LPC interpolation coefficients with mid-ISFs */
extern const float interpol_frac_mid_16k[];             /* LPC interpolation coefficients with mid-ISFs @ 16kHz */
extern const float interpol_frac_mid_relaxprev_12k8[];	/* LPC interpolation coefficients with mid-ISFs @ 16kHz - relaxed prev frame interp */
extern const float interpol_frac_mid_FEC[];             /* LPC interpolation coefficients with mid-ISFs - FEC */
extern const float interpol_frac_mid_relaxprev_16k[];   /* LPC interpolation coefficients with mid-ISFs @ 16kHz - relaxed prev frame interp */
extern const float interpol_frac_mid_16k_FEC[];         /* LPC interpolation coefficients with mid-ISFs @ 16kHz - FEC */
extern const float interpol_frac_mid_relaxprev_pred_12k8[];
extern const float interpol_frac_mid_relaxprev_pred_16k[];

extern const float inter6_2[PIT_FIR_SIZE6_2];
extern const Float32 inter4_2tcx2[4][4];
extern const Float32 inter6_2tcx2[6][4];
typedef struct TCX_LTP_FILTER
{
    const float *filt;
    int length;
} TCX_LTP_FILTER;
extern const TCX_LTP_FILTER tcxLtpFilters[12];

extern const float gain_qua_mless_7b[];                 /* Gain quantization - gain quantization table */
extern const float gain_qua_mless_6b[];                 /* Gain quantization - gain quantization table */
extern const float gain_qua_mless_5b[];                 /* Gain quantization - gain quantization table */
extern const float pred_gain[];                  /* Gain quantization - MA predicition coefficients for gain quantizer */
extern const float t_qua_gain6b[];                      /* Gain quantization - gain quantization table for AMR-WB interoperable mode */
extern const float t_qua_gain7b[];                      /* Gain quantization - gain quantization table for AMR-WB interoperable mode */
extern const float Es_pred_qua_5b[];                    /* Gain quantization - quantization table for scaled innovation energy prediciton */
extern const float Es_pred_qua_4b[];                    /* Gain quantization - quantization table for scaled innovation energy prediciton */
extern const float Es_pred_qua_3b[];                    /* Gain quantization - quantization table for scaled innovation energy prediciton */
extern const float Es_pred_qua_4b_no_ltp[];             /* Gain quantization - quantization table for scaled innovation energy prediciton */

extern const float b_1sfr[];                            /* Gain quantization - gain estimation constants for gain quantizer at 7.2 and 8.0 kbps */
extern const float b_2sfr[];                            /* Gain quantization - gain estimation constants for gain quantizer at 7.2 and 8.0 kbps */
extern const float b_3sfr[];                            /* Gain quantization - gain estimation constants for gain quantizer at 7.2 and 8.0 kbps */
extern const float b_4sfr[];                            /* Gain quantization - gain estimation constants for gain quantizer at 7.2 and 8.0 kbps */

extern const float gp_gamma_1sfr_8b[];                  /* Gain quantization - gain quantization table for gain quantizer at 7.2 and 8.0 kbps */
extern const float gp_gamma_2sfr_8b[];                  /* Gain quantization - gain quantization table for gain quantizer at 7.2 and 8.0 kbps */
extern const float gp_gamma_3sfr_8b[];                  /* Gain quantization - gain quantization table for gain quantizer at 7.2 and 8.0 kbps */
extern const float gp_gamma_4sfr_8b[];                  /* Gain quantization - gain quantization table for gain quantizer at 7.2 and 8.0 kbps */

extern const float gp_gamma_1sfr_7b[];                  /* Gain quantization - gain quantization table for gain quantizer at 7.2 and 8.0 kbps */
extern const float gp_gamma_2sfr_7b[];                  /* Gain quantization - gain quantization table for gain quantizer at 7.2 and 8.0 kbps */
extern const float gp_gamma_3sfr_7b[];                  /* Gain quantization - gain quantization table for gain quantizer at 7.2 and 8.0 kbps */
extern const float gp_gamma_4sfr_7b[];                  /* Gain quantization - gain quantization table for gain quantizer at 7.2 and 8.0 kbps */

extern const float gp_gamma_1sfr_7b[];                  /* Gain quantization - gain quantization table for gain quantizer at 7.2 and 8.0 kbps */
extern const float gp_gamma_2sfr_7b[];                  /* Gain quantization - gain quantization table for gain quantizer at 7.2 and 8.0 kbps */
extern const float gp_gamma_3sfr_7b[];                  /* Gain quantization - gain quantization table for gain quantizer at 7.2 and 8.0 kbps */
extern const float gp_gamma_4sfr_7b[];                  /* Gain quantization - gain quantization table for gain quantizer at 7.2 and 8.0 kbps */

extern const float gp_gamma_1sfr_6b[];                  /* Gain quantization - gain quantization table for gain quantizer at 7.2 and 8.0 kbps */
extern const float gp_gamma_2sfr_6b[];                  /* Gain quantization - gain quantization table for gain quantizer at 7.2 and 8.0 kbps */
extern const float gp_gamma_3sfr_6b[];                  /* Gain quantization - gain quantization table for gain quantizer at 7.2 and 8.0 kbps */
extern const float gp_gamma_4sfr_6b[];                  /* Gain quantization - gain quantization table for gain quantizer at 7.2 and 8.0 kbps */

extern const short E_ROM_qua_gain5b_const[];
extern const short E_ROM_qua_gain6b_const[];
extern const short E_ROM_qua_gain7b_const[];

extern const float gain_qua_mless[];

extern const float tbl_gain_code_tc[];                  /* TC - code gain quantization table */
extern const float tbl_gain_trans_tc[];                 /* TC - gain quantization table for g_trans */
extern const float glottal_cdbk[];                      /* TC - table of prototype glottal impulses */

extern const int   PI_select_table[23][8];              /* selection table for Pulse indexing          */
extern const int   PI_offset[8][8];                     /* offset table for Pulse indexing             */
extern const short PI_factor[];                         /* EVS_PI factor table for Pulse indexing          */

extern const float deem_tab[];                          /* HF BWE - de-emphasis coefficients       */
extern const float filt_hp[];
extern const float exp_tab_p[];                         /* HF BWE - Table of values exp(-j*w*i)   */
extern const float exp_tab_q[];                         /* HF BWE - Table of values exp(-j*w*i)   */
extern const float HP_gain[];                           /* HF BWE - quantization table for 23.85 */
extern const float fir_6k_8k[];                         /* HF BWE - band-pass filter coefficients */

extern const float b_hp400[];                           /* HF (6-7kHz) BWE - 400Hz HP filter coefficients */
extern const float a_hp400[];                           /* HF (6-7kHz) BWE - 400Hz HP filter coefficients */
extern const float fir_6k_7k[];                         /* HF (6-7kHz) BWE - 6.0 - 7.0 kHz BP filter coefficients */

extern const float low_H[];                             /* Enhacer - 2.0 - 6.4 kHz impulse response with phase dispersion */
extern const float mid_H[];                             /* Enhancer - 3.2 - 6.4 kHz impulse response with phase dispersion */

extern const float filt_lp[1 + L_FILT];
extern const float filt_lp_16kHz[1+L_FILT16k];

extern const float tab_hup_l[SIZ_TAB_HUP_L];            /* NB post-filter */
extern const float tab_hup_s[SIZ_TAB_HUP_S];            /* NB post-filter */

extern const float edct_table_80[];                     /* EDCT */
extern const float edct_table_120[];                    /* EDCT */
extern const float edct_table_320[];                    /* EDCT */
extern const float edct_table_480[];                    /* EDCT */
extern const float edct_table_128[];                    /* EDCT */
extern const float edct_table_160[];                    /* EDCT */
extern const float edct_table_40[];                     /* EDCT */
extern const float edct_table_20[];                     /* EDCT */
extern const float edct_table_64[];
extern const float edct_table_100[];                    /* EDCT */
extern const float edct_table_200[];
extern const float edct_table_240[];
extern const float edct_table_256[];
extern const float edct_table_400[];
extern const float edct_table_600[];                    /* EDCT */

extern const short crit_bins[];                         /* (used only in AMR-WB IO mode) */
extern const float crit_bins_corr[CRIT_NOIS_BAND];
extern const float crit_bands_loc[];                    /* (used only in AMR-WB IO mode) */
extern const float mfreq_loc_LD[];                      /* LD music post-filter */
extern const short mfreq_bindiv_LD[];
extern const float post_dct_wind[OFFSET2];
extern const float MAX_SNR_SNR1_tab[];
extern const float INV_MAX_SNR_tab[];
extern const float sc_qnoise[];

extern const float W_DTX_HO[HO_HIST_SIZE];
extern const float ENR_ATT[5];
extern const float HO_ATT[5];

extern const short hq_swb_bwe_nb_bits[];

/*----------------------------------------------------------------------------------*
 * ISF quantization (AMR-WB IO mode)
 *----------------------------------------------------------------------------------*/

extern const float mean_isf_amr_wb[M];                  /* Mean ISF vector (only in AMR-WB IO mode) */
extern const float mean_isf_noise_amr_wb[];             /* Mean ISF vector for SID frame (only in AMR-WB IO mode) */
extern const float gaus_dico[];                         /* Gaussian codebook */
extern const float gaus_dico_swb[];                     /* Gaussian codebook for SWB TBE */

extern const float dico1_isf[];                         /* ISF codebook - common 1st stage, 1st split (only in AMR-WB IO mode) */
extern const float dico2_isf[];                         /* ISF codebook - common 1st stage, 2nd split (only in AMR-WB IO mode) */

extern const float dico21_isf_46b[];                    /* ISF codebook - 46b, 2nd stage, 1st split (only in AMR-WB IO mode) */
extern const float dico22_isf_46b[];                    /* ISF codebook - 46b, 2nd stage, 2st split (only in AMR-WB IO mode) */
extern const float dico23_isf_46b[];                    /* ISF codebook - 46b, 2nd stage, 3rd split (only in AMR-WB IO mode) */
extern const float dico24_isf_46b[];                    /* ISF codebook - 46b, 2nd stage, 4th split (only in AMR-WB IO mode) */
extern const float dico25_isf_46b[];                    /* ISF codebook - 46b, 2nd stage, 5th split (only in AMR-WB IO mode) */

extern const float dico21_isf_36b[];                    /* ISF codebook - 36b, 2nd stage, 1st split (only in AMR-WB IO mode) */
extern const float dico22_isf_36b[];                    /* ISF codebook - 36b, 2nd stage, 2nd split (only in AMR-WB IO mode) */
extern const float dico23_isf_36b[];                    /* ISF codebook - 36b, 2nd stage, 3rd split (only in AMR-WB IO mode) */

extern const float dico1_ns_28b[];                      /* ISF codebook for SID frames - 28b, 1st split */
extern const float dico2_ns_28b[];                      /* ISF codebook for SID frames - 28b, 2nd spilt */
extern const float dico3_ns_28b[];                      /* ISF codebook for SID frames - 28b, 3rd spilt */
extern const float dico4_ns_28b[];                      /* ISF codebook for SID frames - 28b, 4th spilt */
extern const float dico5_ns_28b[];                      /* ISF codebook for SID frames - 28b, 5th spilt */

extern const float dico1_cng_ev[];
extern const float dico2_cng_ev[];
extern const float dico3_cng_ev[];
extern const float dico4_cng_ev[];
extern const float dico5_cng_ev[];

/*----------------------------------------------------------------------------------*
 * LSF quantization - MSVQ tables
 *----------------------------------------------------------------------------------*/

extern const float stable_ISP[];
extern const float stable_LSP[];
extern const float stable_ISF[];

extern const float UVWB_Ave[];
extern const float UVNB_Ave[];
extern const float SVWB_Ave[];
extern const float SVNB_Ave[];
extern const float IAWB_Ave[];
extern const float IANB_Ave[];
extern const float GENB_Ave[];
extern const float GEWB_Ave[];
extern const float GEWB2_Ave[];
extern const float TRWB_Ave[];
extern const float TRWB2_Ave[];
extern const float means_wb_cleanspeech_lsf16k0[];
extern const float means_swb_cleanspeech_lsf25k6[];
extern const float means_swb_cleanspeech_lsf32k0[];
extern const float ModeMean12[];

extern const float Predictor0[];
extern const float Predictor1[];
extern const float Predictor2[];
extern const float Predictor3[];
extern const float Predictor4[];
extern const float Predictor5[];
extern const float Predictor6[];
extern const float Predictor7[];
extern const float Predictor8[];
extern const float CNG_SN1[];

extern const short CB[];
extern const short CB_p[];
extern const float *const ModeMeans[];
extern const float *const Predictors[];
extern const short CBsizes[];
extern const short CBbits[];

extern const short CBbits_p[];
extern const float vals[NO_LEADERS][MAX_NO_VALS];
extern const int   no_vals[NO_LEADERS];
extern const int   no_vals_ind[NO_LEADERS][MAX_NO_VALS];
extern const int   C[LATTICE_DIM+1][LATTICE_DIM+1];
extern const short BitsVQ[];
extern const short BitsVQ_p[];
extern const int   no_lead[][MAX_NO_SCALES*2];
extern const int   no_lead_p[][MAX_NO_SCALES*2];
extern const float sigma[][16];
extern const float sigma_p[][16];
extern const float inv_sigma[][16];
extern const float inv_sigma_p[][16];
extern const float scales[][MAX_NO_SCALES*2];
extern const float scales_p[][MAX_NO_SCALES*2];
extern const short predmode_tab[][6];
extern const float pl[];
extern const int   pi0[];
extern const unsigned int table_no_cv[];
extern const int   pl_par[];
extern const float *const Quantizers[];
extern const float *const Quantizers_p[];
extern const short cng_sort[];
extern const short perm[][4];
extern const short min_lat_bits_SN[];
extern const short min_lat_bits_pred[];
extern const short offset_in_lvq_mode_SN[][21];
extern const short offset_in_lvq_mode_pred[][32];
extern const short offset_lvq_modes_SN[];
extern const short offset_lvq_modes_pred[];

/*-----------------------------------------------------------------*
 * LSF quantization - BC-TCVQ tables
 *-----------------------------------------------------------------*/

extern const short FixBranch[4][4][N_STAGE_VQ - 4];

extern const short BC_TCVQ_BIT_ALLOC_40B[];

extern const int   NTRANS[4][16];
extern const int   NTRANS2[4][16];

extern const float AR_IntraCoeff[N_STAGE_VQ-1][2][2];
extern const float SN_IntraCoeff[N_STAGE_VQ-1][2][2];

extern const float scale_ARSN[];
extern const float scale_inv_ARSN[];

extern const float AR_TCVQ_CB_SUB1[2][128][2];
extern const float AR_TCVQ_CB_SUB2[2][64][2];
extern const float AR_TCVQ_CB_SUB3[4][32][2];

extern const float SN_TCVQ_CB_SUB1[2][128][2];
extern const float SN_TCVQ_CB_SUB2[2][64][2];
extern const float SN_TCVQ_CB_SUB3[4][32][2];

extern const float AR_SVQ_CB1[32][8];
extern const float AR_SVQ_CB2[16][8];


extern const short uniform_model[];

/*-----------------------------------------------------------------*
 * LSF quantization - mid-frame quantization tables
 *-----------------------------------------------------------------*/

extern const float tbl_mid_gen_wb_2b[];
extern const float tbl_mid_gen_wb_5b[];

extern const float tbl_mid_unv_wb_4b[];
extern const float tbl_mid_unv_wb_5b[];

extern const float tbl_mid_voi_wb_1b[];
extern const float tbl_mid_voi_wb_4b[];
extern const float tbl_mid_voi_wb_5b[];

/*-----------------------------------------------------------------*
 * LSF quantization - Mode 2 quantization tables
 *-----------------------------------------------------------------*/

extern const float lsf_init[16];
extern const Float32 means_wb_31bits_ma_lsf[16];
extern const Float32 means_nb_31bits_ma_lsf[16];
extern const float *lsf_means[2];
extern const float *const lsf_codebook[2][2][TCXLPC_NUMSTAGES];
extern const int lsf_numbits[TCXLPC_NUMSTAGES];
extern const int lsf_dims[TCXLPC_NUMSTAGES];
extern const int lsf_offs[TCXLPC_NUMSTAGES];
extern const float lsf_q_diff_cb_8b_rf[];
extern const float lsf_cdk_nb_gc_stg1[];
extern const float lsf_cdk_nb_gc_stg2[];
extern const float lsf_cdk_nb_gc_stg3[];
extern const float lsf_ind_cdk_nb_gc_stg4[];
extern const float lsf_cdk_nb_vc_stg1[];
extern const float lsf_cdk_nb_vc_stg2[];
extern const float lsf_cdk_nb_vc_stg3[];
extern const float lsf_ind_cdk_nb_vc_stg4[];
extern const float lsf_cdk_wb_gc_stg1[];
extern const float lsf_cdk_wb_gc_stg2[];
extern const float lsf_cdk_wb_gc_stg3[];
extern const float lsf_ind_cdk_wb_gc_stg4[];
extern const float lsf_cdk_wb_vc_stg1[];
extern const float lsf_cdk_wb_vc_stg2[];
extern const float lsf_cdk_wb_vc_stg3[];
extern const float lsf_ind_cdk_wb_vc_stg4[];

extern const float *const lsf_ind_codebook[2][2][TCXLPC_IND_NUMSTAGES];
extern const int lsf_ind_numbits[TCXLPC_IND_NUMSTAGES];
extern const int lsf_ind_dims[TCXLPC_IND_NUMSTAGES];
extern const int lsf_ind_offs[TCXLPC_IND_NUMSTAGES];
extern const Word16 min_distance_thr[2][2];

typedef float lsp_unw_triplet[3];
extern const lsp_unw_triplet p16_gamma0_92to1[16];
extern const lsp_unw_triplet p16_gamma0_94to1[16];


/*------------------------------------------------------------------------------*
 * AVQ - RE8 tables
 *------------------------------------------------------------------------------*/

extern const int select_table22[][9];
extern const int vals_a[][4];                           /* value of leader element */
extern const int vals_q[][4];                           /* code parameter for every leader */
extern const unsigned int Is[];                         /* codebook start address for every leader */
extern const int AA3[];                                 /* A3 - Number of the absolute leaders in codebook Q3 */
extern const int AA4[];                                 /* A4 - Number of the absolute leaders in codebook Q4 */
extern const unsigned int II3[];                        /* I3 - Cardinality offsets for absolute leaders in Q3 */
extern const unsigned int II4[];                        /* I4 - Cardinality offset for absolute leaders in Q4 */
extern const int Da_pos[];                              /* Position of the first absolute leader on a spherical shell (or sphere) */
extern const int Da_nb[];                               /* Number of absolute leaders on a spherical shell */
extern const int Da_id[];                               /* identification code of an absolute leader */
extern const int Da_nq[];                               /* Codebook number for each absolute leader */

/*------------------------------------------------------------------------------*
 * SWB TBE tables
 *------------------------------------------------------------------------------*/

extern const short skip_bands_SWB_TBE[];                /* bands for SWB TBE quantisation */
extern const short skip_bands_WB_TBE[];                 /* bands for WB TBE quantisation  */

extern const float interpol_frac_shb[];

extern const float AP1_STEEP[];                         /* All pass filter coeffs for interpolation and decimation by a factor of 2 */
extern const float AP2_STEEP[];                         /* All pass filter coeffs for interpolation and decimation by a factor of 2 */
extern const float STEPS[];                             /* Granuality in conversion from lpc to lsp */

extern const float cos_fb_exc[];
extern const float recip_order[];

extern const float win_lpc_shb[];                       /* Window for calculating SHB LPC coeffs */
extern const float win_lpc_hb_wb[];
extern const float ola_win_shb_switch_fold[];
extern const float win_flatten[];                       /* Window for calculating whitening filter for SHB excitation */
extern const float win_flatten_4k[];                    /* Window for calculating whitening filter for WB excitation */
extern const float window_shb[];                        /* Overlap add window for SHB excitation used in anal and synth */
extern const float window_shb_32k[];                    /* Upsampled overlap add window for SHB excitation used transition generation */
extern const float subwin_shb[];                        /* Short overlap add window for SHB excitation used in anal and synth  */
extern const float window_wb[];
extern const float subwin_wb[];                         /* Short overlap add window for SHB excitation used in anal and synth  */

extern const float Hilbert_coeffs[4*NUM_HILBERTS][HILBERT_ORDER1+1];

extern const float wac[];
extern const float wac_swb[];

extern const float wb_bwe_lsfvq_cbook_8bit[];
extern const float lbr_wb_bwe_lsfvq_cbook_2bit[];
extern const float swb_tbe_lsfvq_cbook_8b[];
extern const float SHBCB_SubGain5bit[];                 /* 5 bit Quantizer table for SHB gain shapes */
extern const float HBCB_SubGain5bit[];                  /* 5-bit TD WB BWE temporal shaping codebook */
extern const float SHBCB_FrameGain64[];                 /* 6 bit Quantizer table for SHB overall gain */
extern const float SHBCB_FrameGain16[];

extern const float full_band_bpf_1[][5];
extern const float full_band_bpf_2[][5];
extern const float full_band_bpf_3[][5];

extern const float lsf_q_cb_4b[];                       /* 4 bit differential scalar quantizer table for TD SWB BWE LSFs 1 and 2*/
extern const float lsf_q_cb_3b[];                       /* 3 bit differential scalar quantizer table for TD SWB BWE LSFs 3, 4 and 5*/
extern const float *const lsf_q_cb[];                   /* Codebook array for each LSF */
extern const short lsf_q_cb_size[];                     /* Size of each element of the above */
extern const short lsf_q_num_bits[];                    /* Size of each element of the above, in bits */
extern const float mirror_point_q_cb[];                 /* LSF mirroring point codebook */
extern const float lsf_grid[4][5];                      /* LSF mirroring adjustment grid */
extern const float grid_smoothing[];                    /* LSF mirroring smoothing table */

extern const float overlap_coefs[];                     /* HR SWB BWE - overlap coefficients */
extern const float overlap_coefs_48kHz[];               /* HR SWB BWE - overlap coefficients @48kHz */
extern const float swb_hr_env_code1[];                  /* HR SWB BWE - envelope Q table - first two subabnds in non-transient frames */
extern const float swb_hr_env_code2[];                  /* HR SWB BWE - envelope Q table - second two subabnds in non-transient frames*/
extern const float swb_hr_env_code3[];                  /* HR SWB BWE - envelope Q table - two subabnds in transient frames */

extern const float allpass_poles_3_ov_2[9];
extern const float decimate_3_ov_2_lowpass_num[3];
extern const float decimate_3_ov_2_lowpass_den[3];


/*------------------------------------------------------------------------------*
 * WB BWE tables
 *------------------------------------------------------------------------------*/

extern const float F_2_5[64];

/*------------------------------------------------------------------------------*
 * SWB BWE tables
 *------------------------------------------------------------------------------*/

extern const short swb_bwe_trans_subband[];
extern const short swb_bwe_trans_subband_width[];
extern const short swb_bwe_subband[];
extern const float swb_inv_bwe_subband_width[];
extern const short swb_bwe_sm_subband[];
extern const float smooth_factor[];
extern const short fb_bwe_subband[];
extern const float fb_inv_bwe_subband_width[];
extern const short fb_bwe_sm_subband[];
extern const float fb_smooth_factor[];
extern const float EnvCdbk11 [];
extern const float EnvCdbk1st [];
extern const float EnvCdbk2nd [];
extern const float EnvCdbk3rd [];
extern const float EnvCdbk4th [];
extern const float EnvCdbkFB[];
extern const float Env_TR_Cdbk1 [];
extern const float Env_TR_Cdbk2 [];
extern const float w_NOR[];
extern const float Mean_env[];
extern const float Mean_env_fb[];
extern const float Mean_env_tr[];

/*------------------------------------------------------------------------------*
 * ACEPL/HQ core switching tables
 *------------------------------------------------------------------------------*/

extern const float hp12800_32000[];
extern const float hp16000_32000[];
extern const float hp12800_48000[];
extern const float hp16000_48000[];
extern const float hp12800_16000[];


extern const double cu15[28][3];
extern const double cu4[6][3];
extern const short ct2[7][13];

/*------------------------------------------------------------------------------*
 * HQ core tables
 *------------------------------------------------------------------------------*/

extern const float window_48kHz[];
extern const float window_256kHz[];
extern const float half_overlap_25[];
extern const float half_overlap_48[];
extern const float half_overlap_int[];
extern const float small_overlap_25[];
extern const float small_overlap_48[];
extern const float small_overlap_int[];
extern const float window_8_16_32kHz[];

extern const float short_window_48kHz[];
extern const float short_window_32kHz[];
extern const float short_window_16kHz[];
extern const float short_window_8kHz[];

extern const float wscw16q15[];
extern const float wscw16q15_8[];
extern const float wscw16q15_16[];
extern const float wscw16q15_32[];

/* Band structure */
extern const short band_len[];
extern const short band_start[];
extern const short band_end[];
extern const short band_len_wb[];
extern const short band_start_wb[];
extern const short band_end_wb[];
extern const short band_len_harm[];
extern const short band_start_harm[];
extern const short band_end_harm[];
extern const float rat[SFM_N_WB];

extern const short intl_bw_16[N_INTL_GRP_16];
extern const short intl_bw_32[N_INTL_GRP_32];
extern const short intl_bw_48[N_INTL_GRP_48];
extern const short intl_cnt_16[N_INTL_GRP_16];
extern const short intl_cnt_32[N_INTL_GRP_32];
extern const short intl_cnt_48[N_INTL_GRP_48];
extern const short norm_order_48[NB_SFM];
extern const short norm_order_32[SFM_N_SWB];
extern const short norm_order_16[SFM_N_WB];

extern const float dicn_pg[45];
extern const short expPkEnrg_tbl[45];
extern const int   manPkEnrg_tbl[45];
extern const int   E_max5_tbl[40];

extern const float thren_pg[44];

extern const float dicn[40];
extern const float dicn_inv[40];
extern const float thren[39];
extern const short dicnlg2[40];
extern const short huffnorm[32];
extern const short huffsizn[32];
extern const short huffcoef[60];
extern const short pgain_huffnorm[32];
extern const short pgain_huffsizn[32];

extern const short resize_huffnorm[32];
extern const short resize_huffsizn[32];

extern const short huffnorm_tran[32];
extern const short huffsizn_tran[32];

extern const short sfm_width[20];
extern const short a_map[20];
extern const short subf_norm_groups[4][11];

extern const Word32 SQRT_DIM_fx[];

/* HQ inner_frame signallisation table */
extern const short inner_frame_tbl[];

/* NB short win: 7200/8000/9600, 13200/16400 */
extern const short band_width_40_4_6_0_0_0[4];
extern const short band_width_40_5_6_0_0_0[5];

/* NB long win: 7200, 8000, 9600, 13200, 16400 */
extern const short band_width_160_18_6_4_0_0[18];
extern const short band_width_160_17_6_3_0_0[17];
extern const short band_width_160_14_6_3_0_0[14];
extern const short band_width_160_13_6_2_0_0[13];

/* WB short win: 13200/16400 */
extern const short band_width_80_7_6_0_0_0[7];

/* WB long win: 13200, 16400 */
extern const short band_width_320_18_6_3_0_0[18];
extern const short band_width_320_20_6_3_0_0[20];

/* SWB short win: 13200, 16400 */
extern const short band_width_142_8_8_0_0_0[8];
extern const short band_width_160_8_8_0_0_0[8];

/* SWB long win: 13200, 16400 */
extern const short band_width_568_22_6_2_0_0[22];
extern const short band_width_640_24_6_4_0_0[24];

/* LR-MDCT configuration tables */
extern const Xcore_Config xcore_config_8kHz_007200bps_long;
extern const Xcore_Config xcore_config_8kHz_008000bps_long;
extern const Xcore_Config xcore_config_8kHz_013200bps_long;
extern const Xcore_Config xcore_config_8kHz_016400bps_long;

extern const Xcore_Config xcore_config_8kHz_007200bps_short;
extern const Xcore_Config xcore_config_8kHz_008000bps_short;
extern const Xcore_Config xcore_config_8kHz_013200bps_short;
extern const Xcore_Config xcore_config_8kHz_016400bps_short;

extern const Xcore_Config xcore_config_16kHz_013200bps_long;
extern const Xcore_Config xcore_config_16kHz_016400bps_long;

extern const Xcore_Config xcore_config_16kHz_013200bps_short;
extern const Xcore_Config xcore_config_16kHz_016400bps_short;

extern const Xcore_Config xcore_config_32kHz_013200bps_long;
extern const Xcore_Config xcore_config_32kHz_016400bps_long;

extern const Xcore_Config xcore_config_32kHz_013200bps_short;
extern const Xcore_Config xcore_config_32kHz_016400bps_short;


extern const short Nb[NB_SFM];
extern const short LNb[ NB_SFM];

extern const Word32 pow_getbitsfrompulses_fx[16];
extern const Word32 table_logcum_fx[563];
extern const Word16 DDP_fx[4];
extern const float DDP[4];


extern const short   step_tcq[8][STATES];
extern const short   denc[8][STATES];
extern const short   ddec[8][STATES];

extern const short step_LSB[STATES_LSB][2];
extern const short denc_LSB[STATES_LSB][2];
extern const short dqnt_LSB[STATES_LSB][4];
extern const short dstep_LSB[4][2];
extern const short ddec_LSB[4][2];

extern const short nextstate[STATES][2];

extern const short fine_gain_bits[];
extern const float *const finegain[];

extern const unsigned char hBitsMinus1_N01[2];
extern const unsigned char hBitsMinus1_N02[65];
extern const unsigned char hBitsMinus1_N03[65];
extern const unsigned char hBitsMinus1_N04[65];
extern const unsigned char hBitsMinus1_N05[54];
extern const unsigned char hBitsMinus1_N06[42];
extern const unsigned char hBitsMinus1_N07[34];
extern const unsigned char hBitsMinus1_N08[29];
extern const unsigned char hBitsMinus1_N09[25];
extern const unsigned char hBitsMinus1_N10[22];
extern const unsigned char hBitsMinus1_N11[19];
extern const unsigned char hBitsMinus1_N12[17];
extern const unsigned char hBitsMinus1_N13[16];
extern const unsigned char hBitsMinus1_N14[14];
extern const unsigned char hBitsMinus1_N15[13];
extern const unsigned char hBitsMinus1_N16[13];
extern const unsigned char hBitsMinus1_N17[12];
extern const unsigned char hBitsMinus1_N18[12];
extern const unsigned char hBitsMinus1_N19[11];
extern const unsigned char hBitsMinus1_N20[11];
extern const unsigned char hBitsMinus1_N21[10];
extern const unsigned char hBitsMinus1_N22[10];
extern const unsigned char hBitsMinus1_N23[10];
extern const unsigned char hBitsMinus1_N24[10];
extern const unsigned char hBitsMinus1_N25[9];
extern const unsigned char hBitsMinus1_N26[9];
extern const unsigned char hBitsMinus1_N27[9];
extern const unsigned char hBitsMinus1_N28[9];
extern const unsigned char hBitsMinus1_N29[9];
extern const unsigned char hBitsMinus1_N30[8];
extern const unsigned char hBitsMinus1_N31[8];
extern const unsigned char hBitsMinus1_N32[8];
extern const unsigned char hBitsMinus1_N33[8];
extern const unsigned char hBitsMinus1_N34[8];
extern const unsigned char hBitsMinus1_N35[8];
extern const unsigned char hBitsMinus1_N36[8];
extern const unsigned char hBitsMinus1_N37[8];
extern const unsigned char hBitsMinus1_N38[8];
extern const unsigned char hBitsMinus1_N39[8];
extern const unsigned char hBitsMinus1_N40[8];
extern const unsigned char hBitsMinus1_N41[7];
extern const unsigned char hBitsMinus1_N42[7];
extern const unsigned char hBitsMinus1_N43[7];
extern const unsigned char hBitsMinus1_N44[7];
extern const unsigned char hBitsMinus1_N45[7];
extern const unsigned char hBitsMinus1_N46[7];
extern const unsigned char hBitsMinus1_N47[7];
extern const unsigned char hBitsMinus1_N48[7];
extern const unsigned char hBitsMinus1_N49[7];
extern const unsigned char hBitsMinus1_N50[7];
extern const unsigned char hBitsMinus1_N51[7];
extern const unsigned char hBitsMinus1_N52[7];
extern const unsigned char hBitsMinus1_N53[7];
extern const unsigned char hBitsMinus1_N54[7];
extern const unsigned char hBitsMinus1_N55[7];
extern const unsigned char hBitsMinus1_N56[7];
extern const unsigned char hBitsMinus1_N57[7];
extern const unsigned char hBitsMinus1_N58[7];
extern const unsigned char hBitsMinus1_N59[7];
extern const unsigned char hBitsMinus1_N60[7];
extern const unsigned char hBitsMinus1_N61[6];
extern const unsigned char hBitsMinus1_N62[6];
extern const unsigned char hBitsMinus1_N63[6];
extern const unsigned char hBitsMinus1_N64[6];
extern const unsigned char* const hBitsN[65];
extern const short dsHighDiracsTab[43];
extern const unsigned int intLimCDivInvDQ31[68];
extern const unsigned char obtainEnergyQuantizerDensity_f[57];
extern const short lim_neg_inv_tbl_fx[11];
extern const short fg_inv_tbl_fx[13] ;
/* functions and tables for pvq_indexing */
extern const unsigned int exactdivodd[ODD_DIV_SIZE] ;

extern const float gain_att[];
extern const float att_step[];
extern const float gain_qlow[];
extern const short gain_cb_size[];
extern const float stab_trans[];
extern const float env_stab_tp[2][2];
/*----------------------------------------------------------------------------------*
 * SWB BWE for LR MDCT core
 *----------------------------------------------------------------------------------*/
extern const float gain_table[NB_SWB_SUBBANDS];
/* HQ_NORMAL mode */
extern const short bits_lagIndices[NB_SWB_SUBBANDS];
extern const short subband_offsets_12KBPS[NB_SWB_SUBBANDS];
extern const short subband_offsets_16KBPS[NB_SWB_SUBBANDS];
extern const short subband_search_offsets[NB_SWB_SUBBANDS];
extern const short bw_SPT_tbl[2][SPT_SHORTEN_SBNUM];

/* HQ_HARMONIC mode */
extern const short bits_lagIndices_mode0_Har[NB_SWB_SUBBANDS_HAR_SEARCH_SB];
extern const short subband_offsets_sub5_13p2kbps_Har[NB_SWB_SUBBANDS_HAR];
extern const short subband_search_offsets_13p2kbps_Har[NB_SWB_SUBBANDS_HAR_SEARCH_SB];
extern const short subband_offsets_sub5_16p4kbps_Har[NB_SWB_SUBBANDS_HAR];
extern const short subband_search_offsets_16p4kbps_Har[NB_SWB_SUBBANDS_HAR_SEARCH_SB] ;

/*----------------------------------------------------------------------------------*
 * SC-VBR
 *----------------------------------------------------------------------------------*/

/* NELP filter coefficients */
extern const float bp1_num_coef_wb[5];
extern const float bp1_den_coef_wb[5];
extern const float shape1_num_coef[11];
extern const float shape1_den_coef[11];
extern const float shape2_num_coef[11];
extern const float shape2_den_coef[11];
extern const float shape3_num_coef[11];
extern const float shape3_den_coef[11];
extern const float txlpf1_num_coef[11];
extern const float txlpf1_den_coef[11];
extern const float txhpf1_num_coef[11];
extern const float txhpf1_den_coef[11];
extern const float bp1_num_coef_nb_fx_order7[8];
extern const float bp1_den_coef_nb_fx_order7[8];
extern const float num_nelp_lp[NELP_LP_ORDER+1];
extern const float den_nelp_lp[NELP_LP_ORDER+1];

extern const float UVG1CB_WB[UVG1_CBSIZE][2];
extern const float UVG2CB1_WB[UVG2_CBSIZE][5];
extern const float UVG2CB2_WB[UVG2_CBSIZE][5];

extern const float UVG1CB_NB[UVG1_CBSIZE][2];
extern const float UVG2CB1_NB[UVG2_CBSIZE][5];
extern const float UVG2CB2_NB[UVG2_CBSIZE][5];

extern const float frac_4sf[NB_SUBFR + 2];

extern const float erb_WB[NUM_ERB_WB + 1];
extern const float erb_NB[NUM_ERB_NB + 1];

extern const float AmpCB1_WB[64][10];
extern const float AmpCB2_WB[64][NUM_ERB_WB - 11];

extern const float AmpCB1_NB[64][10];
extern const float AmpCB2_NB[64][NUM_ERB_NB - 11];

extern const float PowerCB_WB[64][2];
extern const float PowerCB_NB[64][2];
extern const float sinc[8][12];

extern const float hvq_thr_adj[5];
extern const float hvq_peak_cb[1024];
extern const float hvq_class_c[16];
extern const short hvq_cb_search_overlap24k[17];
extern const short hvq_cb_search_overlap32k[21];

extern const short hvq_pg_huff_offset[NUM_PG_HUFFLEN];
extern const short hvq_pg_huff_thres[NUM_PG_HUFFLEN];
extern const short hvq_pg_huff_tab[32];

extern const short hvq_cp_huff_len[52];
extern const short hvq_cp_huff_val[52];
extern const short hvq_cp_layer1_map5[HVQ_CP_MAP_LEN];

extern const short hvq_cp_huff_thres[HVQ_CP_HUFF_NUM_LEN];
extern const short hvq_cp_huff_offset[HVQ_CP_HUFF_NUM_LEN];
extern const short hvq_cp_huff_tab[52];

extern const short pow2_angle_res[8];

/*------------------------------------------------------------------------------*
 * GSC mode
 *------------------------------------------------------------------------------*/

extern const float sin_table256[];

extern const short gsc_sfm_start[];
extern const short gsc_sfm_end[];
extern const short gsc_sfm_size[];

extern const float sm_table[];

extern const float mfreq_loc[];
extern const short mfreq_bindiv_loc[];

extern const short mfreq_bindiv[];
extern const float mean_gp[];
extern const float dic_gp[];

extern const float Gain_mean[];
extern const float Gain_meanHR[];
extern const float Gain_mean_dic[];
extern const float Gain_mean_dicHR[];

extern const float YG_mean16[];
extern const float YG_dicMR_1[];
extern const float YG_dicMR_2[];
extern const float YG_dicMR_3[];
extern const float YG_dicMR_4[];

extern const float mean_m[];
extern const float mean_gain_dic[];

extern const float YGain_mean_LR[];
extern const float YGain_dic1_LR[];
extern const float YGain_dic2_LR[];
extern const float YGain_dic3_LR[];
extern const float Gain_dic2_NBHR[];
extern const float Gain_dic3_NBHR[];

extern const short GSC_freq_bits[];
extern const short GSC_freq_DL0_bits[];

extern const float Gain_meanNB[];
extern const float Gain_mean_dicNB[];
extern const float Mean_dic_NB[];
extern const float Gain_dic1_NB[];
extern const float Gain_dic2_NB[];
extern const float Gain_dic3_NB[];


/*------------------------------------------------------------------------------*
 * FFT transform
 *------------------------------------------------------------------------------*/

extern const short Odx_fft64[64];
extern const float w_fft64[32];
extern const short Ip_fft64[6];
extern const short Odx_fft32_15[32];
extern const float w_fft32[16];
extern const short Ip_fft32[6];
extern const short Odx_fft32_5[32];
extern const short Odx_fft16[16];
extern const float w_fft16[8];
extern const short Ip_fft16[6];
extern const float w_fft8[8];
extern const short Ip_fft8[6];
extern const short Idx_dortft80[80];
extern const short Idx_dortft120[120];
extern const short Idx_dortft160[160];
extern const short Idx_dortft320[320];
extern const short Idx_dortft480[480];
extern const short Ip_fft128[10];
extern const float w_fft128[64];
extern const short Ip_fft256[10];
extern const float w_fft256[128];
extern const short Ip_fft512[18];
extern const float w_fft512[256];
extern const short Idx_dortft40[40];
extern const short Odx_fft8_5[8];
extern const short ip_edct2_64[6];
extern const float w_edct2_64[80];
extern const short Idx_dortft20[20];
extern const short Odx_fft4_5[4];
extern const float w_fft4[2];
extern const short Ip_fft4[6];

/*----------------------------------------------------------------------------------*
 * FEC for HQ core
 *----------------------------------------------------------------------------------*/

extern const float Asr_LP32[41];
extern const float Asr_LP16[21];
extern const float Asr_LP48[61];

extern const short Num_bands_NB[];
extern const float SmoothingWin_NB875[];
extern const float SmoothingWin_NB2[];

/*----------------------------------------------------------------------------------*
 * CLDFB
 *----------------------------------------------------------------------------------*/

extern const int freqTable[2];

extern const float CLDFB80_10[100];
extern const float CLDFB80_16[160];
extern const float CLDFB80_20[200];
extern const float CLDFB80_30[300];
extern const float CLDFB80_32[320];
extern const float CLDFB80_40[400];
extern const float CLDFB80_60[600];

extern const float *const CLDFB80_COEF[];

extern const float rot_vec_ana_re_L10[5];
extern const float rot_vec_ana_im_L10[5];
extern const float rot_vec_ana_re_L16[8];
extern const float rot_vec_ana_im_L16[8];
extern const float rot_vec_ana_re_L20[10];
extern const float rot_vec_ana_im_L20[10];
extern const float rot_vec_ana_re_L30[15];
extern const float rot_vec_ana_im_L30[15];
extern const float rot_vec_ana_re_L32[16];
extern const float rot_vec_ana_im_L32[16];
extern const float rot_vec_ana_re_L40[20];
extern const float rot_vec_ana_im_L40[20];
extern const float rot_vec_ana_re_L60[30];
extern const float rot_vec_ana_im_L60[30];

extern const float rot_vec_syn_re_L10[5];
extern const float rot_vec_syn_im_L10[5];
extern const float rot_vec_syn_re_L16[8];
extern const float rot_vec_syn_im_L16[8];
extern const float rot_vec_syn_re_L20[10];
extern const float rot_vec_syn_im_L20[10];
extern const float rot_vec_syn_re_L30[15];
extern const float rot_vec_syn_im_L30[15];
extern const float rot_vec_syn_re_L32[16];
extern const float rot_vec_syn_im_L32[16];
extern const float rot_vec_syn_re_L40[20];
extern const float rot_vec_syn_im_L40[20];
extern const float rot_vec_syn_re_L60[30];
extern const float rot_vec_syn_im_L60[30];

extern const float bpf_weights_16[CLDFB_NO_COL_MAX];

extern const float CNG_details_codebook[64][NUM_ENV_CNG];
/*----------------------------------------------------------------------------------*
 * FD CNG
 *----------------------------------------------------------------------------------*/

extern const int d_array[18];
extern const float m_array[18];
extern const float msQeqInvAv_thresh[3];
extern const float msNoiseSlopeMax[4];

extern const SCALE_SETUP scaleTable[20];
extern const SCALE_SETUP scaleTable_cn_only[18];
extern const float scaleTable_cn_only_amrwbio[3][2];

extern const int sidparts_encoder_noise_est[24];


extern const FD_CNG_SETUP FdCngSetup_nb;
extern const FD_CNG_SETUP FdCngSetup_wb1;
extern const FD_CNG_SETUP FdCngSetup_wb2;
extern const FD_CNG_SETUP FdCngSetup_wb3;
extern const FD_CNG_SETUP FdCngSetup_swb1;
extern const FD_CNG_SETUP FdCngSetup_swb2;

extern const short maxN_37bits;
extern const short maxC_37bits;
extern const short stages_37bits;
extern const int levels_37bits[6];
extern const short bits_37bits[6];
extern const float * const cdk_37bits[];

extern const float fftSineTab640[321];

extern const float olapWinAna512[512];
extern const float olapWinAna640[640];

extern const float olapWinSyn256[256];
extern const float olapWinSyn320[320];


/*----------------------------------------------------------------------------------*
 * TCX
 *----------------------------------------------------------------------------------*/

extern float const gain_corr_fac[];
extern float const gain_corr_inv_fac[];

extern const SCALE_TCX_SETUP scaleTcxTable[13];


/*----------------------------------------------------------------------------------*
 * Arithmetic coder
 *----------------------------------------------------------------------------------*/

extern const unsigned char ari_lookup_s17_LC[4096];
extern const unsigned short ari_pk_s17_LC_ext[64][18];

extern const int NumRatioBits[2][17];
extern const float Ratios_WB_2[32];
extern const float Ratios_WB_3[32];
extern const float Ratios_WB_4[32];
extern const float Ratios_WB_5[32];
extern const float Ratios_WB_6[32];
extern const float Ratios_WB_7[32];
extern const float Ratios_WB_8[16];
extern const float Ratios_WB_9[16];
extern const float Ratios_WB_10[16];
extern const float Ratios_WB_11[16];
extern const float Ratios_WB_12[16];
extern const float Ratios_WB_13[16];
extern const float Ratios_WB_14[16];
extern const float Ratios_WB_15[16];
extern const float Ratios_WB_16[4];
extern const float Ratios_WB_17[4];
extern const float Ratios_WB_18[4];
extern const float Ratios_NB_2[32];
extern const float Ratios_NB_3[16];
extern const float Ratios_NB_4[16];
extern const float Ratios_NB_5[16];
extern const float Ratios_NB_6[16];
extern const float Ratios_NB_7[16];
extern const float Ratios_NB_8[16];
extern const float Ratios_NB_9[8];
extern const float Ratios_NB_10[8];
extern const float Ratios_NB_11[8];
extern const float Ratios_NB_12[8];
extern const float Ratios_NB_13[4];
extern const float Ratios_NB_14[4];
extern const float Ratios_NB_15[4];
extern const float Ratios_NB_16[4];
extern const float Ratios_NB_17[4];
extern const float Ratios_NB_18[4];
extern const float *const Ratios[2][17];

extern const Word16 qGains[2][1 << kTcxHmNumGainBits];

/*----------------------------------------------------------------------------------*
 * Innovative codebook
 *----------------------------------------------------------------------------------*/

extern const PulseConfig PulseConfTable[];

/*----------------------------------------------------------------------------------*
 * TNS
 *----------------------------------------------------------------------------------*/

extern struct TnsParameters const tnsParametersIGF32kHz_LowBR[1];
extern struct TnsParameters const tnsParameters32kHz[2];
extern struct TnsParameters const tnsParameters32kHz_grouped[2];
extern struct TnsParameters const tnsParameters16kHz[1];
extern struct TnsParameters const tnsParameters16kHz_grouped[2];
extern struct TnsParameters const tnsParameters48kHz_grouped[2];

extern float const tnsAcfWindow[TNS_MAX_FILTER_ORDER];

extern ParamsBitMap const tnsEnabledSWBTCX20BitMap;
extern ParamsBitMap const tnsEnabledSWBTCX10BitMap;
extern ParamsBitMap const tnsEnabledWBTCX20BitMap;

extern const Coding codesTnsAllCoeffs[];

extern const Coding codesTnsCoeff0SWBTCX20[];
extern const Coding codesTnsCoeff0SWBTCX10[];
extern const Coding codesTnsCoeff1SWBTCX20[];
extern const Coding codesTnsCoeff1SWBTCX10[];
extern const Coding codesTnsCoeff2SWBTCX20[];
extern const Coding codesTnsCoeff2SWBTCX10[];
extern const Coding codesTnsCoeff3SWBTCX20[];
extern const Coding codesTnsCoeff3SWBTCX10[];
extern const Coding codesTnsCoeff4SWBTCX20[];
extern const Coding codesTnsCoeff4SWBTCX10[];
extern const Coding codesTnsCoeff5[];
extern const Coding codesTnsCoeff6[];

extern const Coding codesTnsCoeff0WBTCX20[];
extern const Coding codesTnsCoeff1WBTCX20[];
extern const Coding codesTnsCoeff2WBTCX20[];
extern const Coding codesTnsCoeff3WB[];

extern const Coding codesTnsCoeff456[];
extern const Coding codesTnsCoeff7[];

extern int const nTnsCoeffCodes;

extern const Coding * const codesTnsCoeffSWBTCX20[];
extern const Coding * const codesTnsCoeffSWBTCX10[];
extern const Coding * const codesTnsCoeffWBTCX20[];
extern int const nTnsCoeffTables;

extern const Coding codesTnsOrderTCX20[];
extern const Coding codesTnsOrderTCX10[];
extern const Coding codesTnsOrder[];
extern int const nTnsOrderCodes;

extern float const tnsCoeff4[16];

/**********************************************************************/ /**
igf settings structure for each bitrate mode
**************************************************************************/
typedef struct igf_mode_type
{
    int     sampleRate;
    int     frameLength;
    int     bitRate;
    int     igfMinFq;
    float   transFac;
    int     maxHopsize;
} IGF_MODE,*H_IGF_MODE;

extern const IGF_MODE igfMode[15];
extern const int swb_offset_LB_new[15][IGF_MAX_SFB];
extern const float igf_whitening_TH[15][2][IGF_MAX_TILES];
extern const short cf_off_se01_tab[9];
extern const short cf_off_se10_tab;
extern const short cf_off_se02_tab[9][IGF_CTX_COUNT];
extern const short cf_off_se11_tab[IGF_CTX_COUNT][IGF_CTX_COUNT];
extern const unsigned short cf_se00_tab[IGF_SYMBOLS_IN_TABLE + 1];
extern const unsigned short cf_se01_tab[9][IGF_SYMBOLS_IN_TABLE + 1];
extern const unsigned short cf_se02_tab[9][IGF_CTX_COUNT][IGF_SYMBOLS_IN_TABLE + 1];
extern const unsigned short cf_se10_tab[IGF_SYMBOLS_IN_TABLE + 1];
extern const unsigned short cf_se11_tab[IGF_CTX_COUNT][IGF_CTX_COUNT][IGF_SYMBOLS_IN_TABLE + 1];

extern const int bwMode2fs[4];

extern const float normReciprocal[CHEAP_NORM_SIZE];
extern const float *w_a[7];

extern const PWord16 SineTable512_fx[];
extern const Word16 ldCoeff[7];
extern const UWord32 exp2_tab_long[32];
extern const UWord32 exp2w_tab_long[32];
extern const UWord32 exp2x_tab_long[32];
extern const Word16 invTable[INV_TABLE_SIZE+1];
extern const Word16 sqrtTable[SQRT_TABLE_SIZE+1];
extern const Word16 invSqrtTable[SQRT_TABLE_SIZE+1];

/* ACELP pulse coding */
extern const int low_len[10];
extern const int low_mask[10];
extern const int indx_fact[10];
extern const int index_len[3];
extern const int index_mask[3];

#endif

