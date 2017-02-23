/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#ifndef ROM_ENC_H
#define ROM_ENC_H

#include <stdio.h>
#include "options.h"
#include "stat_enc.h"
#include "cnst.h"

/*----------------------------------------------------------------------------------*
 * General tables
 *----------------------------------------------------------------------------------*/

extern const float sqrt_han_window[];                   /* Half of the square root hanning window */
extern const short tipos[];                             /* Starting points for pulse position search in Algebraic innovation codebook */
extern const Float32 E_ROM_inter4_1[PIT_UP_SAMP * L_INTERPOL1 + 1];
extern const Float32 E_ROM_inter6_1[PIT_UP_SAMP6 * L_INTERPOL1 + 1];
extern const Float32 E_ROM_inter4_1[PIT_UP_SAMP * L_INTERPOL1 + 1];
extern const Float32 E_ROM_inter6_1[PIT_UP_SAMP6 * L_INTERPOL1 + 1];
extern const float W_HIST[DTX_HIST_SIZE];               /* CNG & DTX - table for calculation of average excitation energy */

extern const short bwd_start_bin[];
extern const short bwd_end_bin[];

extern const float h_fir[];                                   /* 2nd order fir filter for wsp, decimation by 2 */

extern const float preemphCompensation[];

/*----------------------------------------------------------------------------------*
 * VAD tables
 *----------------------------------------------------------------------------------*/
extern const short hangover_hd_tbl[3];
extern const short hangover_sf_tbl[6];
/*----------------------------------------------------------------------------------*
 * Open-loop pitch search tables
 *----------------------------------------------------------------------------------*/

extern const short nb_sect_12k8[];
extern const short nb_subsect_12k8[];
extern const short len_12k8[];
extern const short len1_12k8[];
extern const short sublen_12k8[];
extern const short sublen1_12k8[];
extern const short pit_max_12k8[];
extern const short sec_length_12k8[];
extern const short sec_length1_12k8[];

/*----------------------------------------------------------------------------------*
 * LSF quantizer
 *----------------------------------------------------------------------------------*/

extern const int lsf_numlevels[TCXLPC_NUMSTAGES];
extern const int lsf_ind_numlevels[TCXLPC_IND_NUMSTAGES];

extern const short lsf_unified_fit_model_nb[4][16];
extern const short lsf_unified_fit_model_wb[4][16];
extern const short lsf_unified_fit_model_wbhb[4][16];

extern const float Freq_Weight_Com[160];
extern const float Freq_Weight_UV[160];

/*----------------------------------------------------------------------------------*
 * Speech/music classification
 *----------------------------------------------------------------------------------*/

extern const float w[HANG_LEN][HANG_LEN];

extern const float m_speech[];
extern const float invV_speech[];
extern const float lvm_speech[];

extern const float m_music[];
extern const float invV_music[];
extern const float lvm_music[];

extern const float m_noise[];
extern const float invV_noise[];
extern const float lvm_noise[];

extern const float SF[];
extern const float SF_8k[];

/*----------------------------------------------------------------------------------*
 * SWB TBE
 *----------------------------------------------------------------------------------*/

extern const float lpc_weights[];

/*----------------------------------------------------------------------------------*
 * WB, SWB and FB bandwidth detector
 *----------------------------------------------------------------------------------*/

extern const float hann_window_320[];

/*----------------------------------------------------------------------------------*
 * Huffman coding
 *----------------------------------------------------------------------------------*/

extern const short huffsizn_e[32];
extern const short huffsizn_n[32];

extern const short huffnorm_e[32];
extern const short huffnorm_n[32];
extern const short hessize[8];
extern const short hescode[8];

/*----------------------------------------------------------------------------------*
 * VAD
 *----------------------------------------------------------------------------------*/

extern const short BAND_NUM_TAB[5];

extern const float M_inr[16];
extern const float M_ini[16];
extern const float M_r[8];
extern const float M_i[8];
extern const float M_Wr[16];
extern const float M_Wi[16];

extern const short SP_CENTER_BAND_NUM_TAB[5];
extern const float DELTA1[5];
extern const float DELTA2[5];

extern const int NREGION_INDEX_NB[9];
extern const int NREGION_INDEX_WB[13];
extern const int NREGION_INDEX_SWB[16];
extern const int NREGION_INDEX_FB[16];
extern const int ENERGY_BAND_NUM[4];
extern const int *const REGION_INDEX[4];
extern const float MAX_LF_SNR_TAB[4];

extern const float COMVAD_INIT_SNR_DELTA[5];
extern const float LS_MIN_SELENCE_SNR[3];
extern const float LT_MIN_SILENCE_SNR[3];

/*----------------------------------------------------------------------------------*
 * Starting line for the noise measurement in TCX.
 *----------------------------------------------------------------------------------*/
extern const int startLineWB[11];
extern const int startLineSWB[9];


#endif

