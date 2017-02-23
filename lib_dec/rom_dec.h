/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#ifndef ROM_DEC_H
#define ROM_DEC_H

#include <stdio.h>
#include "options.h"
#include "cnst.h"

extern const float h_low[];                      /* LP filter for filtering periodic part of excitation in artificial onset construction after FEC */

extern const int mult_avq_tab[];
extern const int shift_avq_tab[];

extern const short hntable[55];
extern const short hetable[57];
extern const short hestable[15];

extern const float lsf_tab[LPC_SHB_ORDER];

extern const short gw[LGW_MAX];
extern const short gwlpr[LGW_MAX];
extern const float w_hamm32k_2[L_TRANA32k/2];
extern const float w_hamm16k_2[L_TRANA16k/2];
extern const float w_hamm8k_2[L_TRANA8k/2];
extern const float w_hamm_sana32k_2[L_PROT_HAMM_LEN2_32k];
extern const float w_hamm_sana16k_2[L_PROT_HAMM_LEN2_16k];
extern const float w_hamm48k_2[L_TRANA48k/2];
extern const float w_hamm_sana48k_2[L_PROT_HAMM_LEN2_48k];

extern const float h_high3_32[L_FIR_FER2];
extern const float h_high3_16[L_FIR_FER2];


#endif

