/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#ifndef BASOP_PROTO_FUNC_H
#define BASOP_PROTO_FUNC_H

#include "stl.h"
#include "basop_util.h"


/* tcx_lpc_cdk.h */
#define LSF_GAP_VAL(x) (Word16)((x)*2.0f*1.28f)
#define LSFM(x) FL2WORD16_SCALE(x*1.28, 15-1)  /* 14Q1*1.28 */

/* cnst.h */
#define GAMMA1_INV        17809    /* weighting factor (numerator) default:0.92 (1Q14format)     */
#define GAMMA16k_INV      17430    /* weighting factor (numerator) default:0.94 (1Q14format)     */
#define FS_2              16384    /* isf max value (Use in reorder_fx.c) */
#define INT_FS_FX         12800    /* internal sampling frequency         */

void basop_lsp2a_stab(const Word16 *lsp, Word16 *a);
void basop_lsf2lsp(const Word16 lsf[], Word16 lsp[]);
void basop_weight_a(const Word16 *a, Word16 *ap, const Word16 gamma);
void basop_weight_a_inv(const Word16 *a, Word16 *ap, const Word16 inv_gamma);
void basop_E_LPC_a_add_tilt(const Word16 *a, Word16 *ap, Word16 gamma);
void basop_reorder_lsf(Word16 *lsf, const Word16 min_dist, const Word16 n, const Word32 fs);
void basop_E_LPC_f_lsp_a_conversion(const Word16 *lsp, Word16 *a,  const Word16 m);

/* tcx_utils.c */
void basop_lpc2mdct(Word16 *lpcCoeffs, Word16 lpcOrder,
                    Word16 *mdct_gains, Word16 *mdct_gains_exp,
                    Word16 *mdct_inv_gains, Word16 *mdct_inv_gains_exp);

void basop_PsychAdaptLowFreqDeemph(Word32 x[], const Word16 lpcGains[], const Word16 lpcGains_e[], Word16 lf_deemph_factors[]);
void basop_mdct_noiseShaping_interp(Word32 x[], Word16 lg, Word16 gains[], Word16 gains_exp[]);


#endif
