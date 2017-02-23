/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "cnst.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "options.h"
#include "prot.h"
#include "rom_com.h"



void tcx_decoder_memory_update(
    float *xn_buf,          /* i/o: mdct output buffer used also as temporary buffer */
    float *synthout,        /* o: synth                                              */
    int L_frame_glob,       /* i: global frame length                                */
    float *A,               /* i: Quantized LPC coefficients                         */
    Decoder_State *st,      /* i/o : decoder memory state                            */
    float *syn              /* o: st->syn                                            */
)
{
    float tmp;
    float *synth;
    float buf[1+M+L_FRAME_PLUS];
    float preemph_f = st->preemph_fac;

    /*TCX must be aligned with ACELP*/
    assert(st->tcx_cfg.lfacNext<=0);

    /* Output synth */
    mvr2r(xn_buf, synthout, L_frame_glob);

    /* Update synth */
    synth = buf + 1 + M;
    mvr2r(syn, buf, 1+M);
    mvr2r(xn_buf, synth, L_frame_glob);
    mvr2r(synth+L_frame_glob-M-1, syn, 1+M);

    /* Emphasis of synth -> synth_pe */
    tmp = synth[-M-1];
    preemph(synth-M, preemph_f, M+L_frame_glob, &tmp);

    mvr2r(synth+L_frame_glob-M, st->mem_syn2, M);
    mvr2r( synth+L_frame_glob-L_SYN_MEM, st->mem_syn_r, L_SYN_MEM);

    if ( !st->tcxonly || L_frame_glob == L_FRAME16k)
    {
        /* Update excitation */
        assert(L_frame_glob < L_EXC_MEM_DEC);
        mvr2r(st->old_exc+(L_frame_glob), st->old_exc, L_EXC_MEM_DEC-(L_frame_glob) );
        residu(A, M, synth, st->old_exc+L_EXC_MEM_DEC-(L_frame_glob), (L_frame_glob));
    }

    /* Update old_Aq */
    mvr2r(A, st->old_Aq_12_8, M+1);
}

/*---------------------------------------------------------------
 * Residual Quantization
 *--------------------------------------------------------------*/

/* Returns: number of bits used (including "bits") */
int tcx_ari_res_invQ_spec(
    float x_Q[],          /* i/o: quantized spectrum               */
    int L_frame,          /* i: number of lines                    */
    const int prm[],      /* i: bit-stream                         */
    int target_bits,      /* i: number of bits available           */
    int bits,             /* i: number of bits used so far         */
    float deadzone,       /* i: quantizer deadzone                 */
    const float x_fac[]   /* i: spectrum post-quantization factors */
)
{
    int i, j, num_zeros;
    int zeros[L_FRAME_PLUS];
    float fac_m, fac_p, sign;


    /* Limit the number of residual bits */
    target_bits = min(target_bits, NPRM_RESQ);


    /* Requantize the spectrum line-by-line */
    fac_m = deadzone * 0.5f;
    num_zeros = 0;
    for (i=0; i<L_frame; ++i)
    {
        if (bits >= target_bits)
        {
            /* no bits left */
            break;
        }
        if (x_Q[i] != 0)
        {
            if (x_Q[i] > 0)
            {
                sign = x_fac[i];
            }
            else
            {
                sign = -x_fac[i];
            }

            x_Q[i] += sign*(prm[bits++] * 0.5f - fac_m);
        }
        else
        {
            zeros[num_zeros++] = i;
        }
    }

    /* Requantize zeroed-lines of the spectrum */
    fac_p = (1.0f - deadzone)*0.33f*2.0f;
    --target_bits; /* reserve 1 bit for the check below */
    for (j=0; j<num_zeros; ++j)
    {
        if (bits >= target_bits)
        {
            /* 1 or 0 bits left */
            break;
        }

        i = zeros[j];

        if (prm[bits++] != 0)
        {
            x_Q[i] = (2*prm[bits++]-1) * fac_p * x_fac[i];
        }
    }


    return bits;
}

int tcx_res_invQ_gain(
    float *gain_tcx,
    const int *prm,
    int resQBits
)
{
    int bits;


    /*Refine the gain quantization*/
    for(bits=0; bits<min(resQBits,TCX_RES_Q_BITS_GAIN); bits++)
    {
        if(prm[bits]==0)
        {
            *gain_tcx=(*gain_tcx)*gain_corr_inv_fac[bits];
        }
        else
        {
            *gain_tcx=(*gain_tcx)*gain_corr_fac[bits];
        }
    }

    return(bits);
}

int tcx_res_invQ_spec(
    float *x,
    int L_frame,
    const int *prm,
    int resQBits,
    int bits,
    float sq_round,
    const float lf_deemph_factors[]
)
{
    int i;
    float fac_m, fac_p, thres;


    /* Limit the number of residual bits */
    resQBits = min(resQBits, NPRM_RESQ);

    /* Requantize the spectrum line-by-line */
    fac_p = 0.5f - sq_round * 0.5f;
    fac_m = sq_round * 0.5f;
    if (!lf_deemph_factors)
    {
        for (i = 0; (i < L_frame) && (bits < resQBits); i++)
        {
            /* bits < resQBits */
            if (x[i] != 0.0f)
            {
                if (prm[bits++] == 0)
                {
                    x[i] -= (x[i] > 0.0f) ? fac_m : fac_p;
                }
                else
                {
                    x[i] += (x[i] > 0.0f) ? fac_p : fac_m;
                }
            }
        }
        resQBits--; /* Quantize zeroed lines of the spectrum */
        for (i = 0; (i < L_frame) && (bits < resQBits); i++)
        {
            if (x[i] == 0.0f)
            {
                if (prm[bits++] != 0)
                {
                    x[i] = (prm[bits++] * 2.64f - 1.32f) * fac_p;
                }
            }
        }

        return bits;
    }
    for(i=0; i<L_frame; i++)
    {
        if(bits>=resQBits)
        {
            break;
        }
        if(x[i]!=0.f && lf_deemph_factors[i] > 0.5f)
        {
            if(prm[bits++]==0)
            {
                x[i]-=(x[i]>0)?fac_m*lf_deemph_factors[i]:fac_p*lf_deemph_factors[i];
            }
            else
            {
                x[i]+=(x[i]>0)?fac_p*lf_deemph_factors[i]:fac_m*lf_deemph_factors[i];
            }
        }
    }

    /*Quantize zeroed-line of the spectrum*/
    for(i=0; i<L_frame; i++)
    {
        if(bits>=(resQBits-1))  /*need at least two bits*/
        {
            break;
        }
        if(x[i]==0.f && lf_deemph_factors[i] > 0.5f)
        {
            if(prm[bits++]==1)
            {
                thres=(1-sq_round)*0.66f*lf_deemph_factors[i];
                x[i]=(prm[bits++]*2.f-1.f)*thres;
            }
        }
    }

    return bits;
}
