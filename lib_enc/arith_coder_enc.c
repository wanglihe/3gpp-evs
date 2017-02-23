/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"
#include "basop_util.h"
#include "basop_proto_func.h"


/*-------------------------------------------------------------------*
 * Local constants
 *-------------------------------------------------------------------*/

#define LOG2_E  1.44269504089f

#define kMaxNumHeapElems 10

typedef struct HeapElem
{
    float mScore; /* Sort key                    */
    int   mIndex; /* Original index              */
} HeapElem;

typedef struct Heap
{
    HeapElem mElem[2*kMaxNumHeapElems+1];
    int      mSize;
} Heap;


/*-------------------------------------------------------------------*
 * tcx_arith_estimate_scale()
 *
 *
 *-------------------------------------------------------------------*/

static float tcx_arith_estimate_scale(  /* o: estimated SQ scale            */
    const float abs_spectrum[],           /* i: absolute MDCT coefficients    */
    int L_frame,                          /* i: number of spectral lines      */
    const Word16 envelope[],              /* i: scaled envelope (Q15-e)       */
    Word16 envelope_e                     /* i: scaled envelope exponent (Q0) */
)
{
    float scale, tmp;
    int k;

    /* compute normalised standard deviation and determine approximate scale */
    scale = 0.01f;
    for (k = 0; k < L_frame; k++)
    {
        tmp = abs_spectrum[k] * envelope[k];
        scale += tmp * tmp;
    }
    tmp = (float)(1 << (15-envelope_e));
    scale = (float)sqrt((L_frame * tmp*tmp*4.0f) / scale);


    return scale;
}


/*-------------------------------------------------------------------*
 * MinHeapify_i()
 *
 *
 *-------------------------------------------------------------------*/

static void MinHeapify_i(Heap *H, int i)
{
    int left, right, largest;
    HeapElem T;

    left    = 2*i + 1;
    right   = left + 1;
    largest = i;

    if (H->mElem[left].mScore < H->mElem[largest].mScore)
    {
        largest = left;
    }
    if (H->mElem[right].mScore < H->mElem[largest].mScore)
    {
        largest = right;
    }
    while (largest != i)
    {
        T.mIndex = H->mElem[i].mIndex;
        T.mScore = H->mElem[i].mScore;

        H->mElem[i].mIndex = H->mElem[largest].mIndex;
        H->mElem[i].mScore = H->mElem[largest].mScore;

        H->mElem[largest].mIndex = T.mIndex;
        H->mElem[largest].mScore = T.mScore;

        i = largest;

        left    = 2*i + 1;
        right   = left + 1;

        if (H->mElem[left].mScore < H->mElem[largest].mScore)
        {
            largest = left;
        }
        if (H->mElem[right].mScore < H->mElem[largest].mScore)
        {
            largest = right;
        }
    }

    return;
}


/*-------------------------------------------------------------------*
 * tcx_arith_find_max_scale()
 *
 *
 *-------------------------------------------------------------------*/

static float tcx_arith_find_max_scale(
    const float abs_spectrum[],    /* i: absolute MDCT coefficients    */
    int L_frame,                   /* i: number of spectral lines      */
    const Word16 envelope[],       /* i: scaled envelope (Q15-e)       */
    Word16 envelope_e,             /* i: scaled envelope exponent (Q0) */
    const Word16 exps[],           /* i: expfp(-(int)envelope[]/2)     */
    float deadzone                 /* i: deadzone (0.5f = no deadzone) */
)
{
    int i, k, q;
    float p, scale;
    Heap heap;
    Word16 tmpi1, tmpi2;
    float envelope_scale;
    const float limit = -9.70406052784f; /* = ln(1/16384): log of smallest allowed probability */

    /* Find the top most offending lines according to probability estimates */
    heap.mSize = kMaxNumHeapElems;
    heap.mElem[0].mScore = 0;     /* mal: just to silnce the compiler */

    for (i=0; i<kMaxNumHeapElems; ++i)
    {
        heap.mElem[i].mIndex = 0;
        heap.mElem[i].mScore = 0;
    }
    for (; i<2*kMaxNumHeapElems+1; ++i)
    {
        heap.mElem[i].mScore = FLT_MAX;
    }
    for (k=0; k<L_frame; ++k)
    {
        p = envelope[k] * abs_spectrum[k];
        if (p > heap.mElem[0].mScore)
        {
            heap.mElem[0].mScore = p;
            heap.mElem[0].mIndex = k;
            MinHeapify_i(&heap, 0);
        }
    }

    /* Make sure the scale is limited so that the offending lines don't cause probability underflow. */
    /* Also limit scale to avoiding saturation of the gain quantizer */
    scale = 1.0f/(float)sqrt(L_frame*0.5f);
    envelope_scale = -(float)pow(2, envelope_e-16);
    for (i=0; i<heap.mSize; ++i)
    {
        k = heap.mElem[i].mIndex;

        /* Get approximate maximum allowed magnitude */
        q = (int)ceil(((limit - log(1.0f - (exps[k]/32768.0) * (exps[k]/32768.0))) / (envelope[k]*envelope_scale) - 1) / 2.0f);

        /* Refinement: get the exact q */
        powfp_odd2(exps[k], q, &tmpi1, &tmpi2);
        if (tmpi1 - tmpi2 >= 2)
        {
            /* q may be too low */
            powfp_odd2(exps[k], q+1, &tmpi1, &tmpi2);
            while (tmpi1 - tmpi2 >= 2)
            {
                ++q;
                powfp_odd2(exps[k], q+1, &tmpi1, &tmpi2);
            }
        }
        else
        {
            /* q is too high */
            --q;
            powfp_odd2(exps[k], q, &tmpi1, &tmpi2);
            while (tmpi1 - tmpi2 < 2)
            {
                --q;
                powfp_odd2(exps[k], q, &tmpi1, &tmpi2);
            }
        }

        /* Find the largest scale so that the quantized magnitude is at most q */
        p = (q+0.99f-deadzone)/(abs_spectrum[k] + 0.000001f);
        assert((int)(abs_spectrum[k] * p + deadzone) <= q);
        scale = min(scale, p);
    }


    return scale;
}


/*-------------------------------------------------------------------*
 * tcx_arith_find_kMax()
 *
 *
 *-------------------------------------------------------------------*/

static int tcx_arith_find_kMax(  /* o: index of highest freq. nonzero line (-1 if all zeros) */
    const float abs_spectrum[],    /* i: absolute MDCT coefficients    */
    int L_frame,                   /* i: number of spectral lines      */
    float scale,                   /* i: scalar quantizer scale        */
    float deadzone,                /* i: deadzone (0.5f = no deadzone) */
    const int deadzone_flags[]     /* i: line-wise deadzone control    */
)
{
    int kMax;


    kMax = L_frame - 1;
    while ((kMax >= 0) && (abs_spectrum[kMax] * scale < (1.0f - deadzone) + deadzone * deadzone_flags[kMax]))
    {
        kMax--;
    }

    return kMax;
}


/*-------------------------------------------------------------------*
 * tcx_arith_rateloop()
 *
 *
 *-------------------------------------------------------------------*/

static float tcx_arith_rateloop( /* o: best scale                       */
    const float abs_spectrum[],    /* i: absolute MDCT coefficients       */
    int L_frame,                   /* i: number of spectral lines         */
    const Word16 envelope[],       /* i: scaled envelope (Q15-e)          */
    Word16 envelope_e,             /* i: scaled envelope exponent (Q0)    */
    const Word16 exps[],           /* i: expfp(-(int)envelope[]/2)        */
    int target_bits,               /* i: target bit budget                */
    float deadzone,                /* i: deadzone (0.5f = no deadzone)    */
    const int deadzone_flags[],    /* i: line-wise deadzone control       */
    float *target_bits_fac         /* i/o: scale estimator compensation   */
)
{
    int k, kMax, q;
    float s, adjust;
    float fixed_bits[2][N_MAX_ARI];
    float estimator_undershoot;
    float max_complexity;
    int   iter;           /* rate loop iteration counter               */
    float scale;          /* SQ scale factor to try next               */
    float scale_best;     /* best SQ scale factor                      */
    float scale_max;      /* maximum allowable scale factor            */
    float lob;            /* lower bound of SQ scale factor            */
    float hib;            /* upper bound of SQ scale factor            */
    int   flag;           /* 1:bit surplus, -1:bit deficit, 0:unknown  */
    float complexity;     /* cumulative rate loop complexity           */
    float bits;           /* number of bits (approximate)              */
    float envelope_scale;


    scale = tcx_arith_estimate_scale(abs_spectrum, L_frame, envelope, envelope_e);
    scale *= *target_bits_fac;

    scale_max = tcx_arith_find_max_scale(abs_spectrum, L_frame, envelope, envelope_e, exps, deadzone);
    if (scale > scale_max) scale = scale_max;

    scale_best     = scale;
    lob            = 0.0f;
    hib            = 0.0f;
    flag           = 0;
    complexity     = 0;
    bits           = 0;
    max_complexity = 96.0f * L_frame;
    iter           = 0;
    envelope_scale = (float)pow(2, envelope_e-15);

    estimator_undershoot = 0;
    /* Precalculate fixed bit costs */
    for (k=0; k<L_frame; ++k)
    {
        s = envelope[k] * envelope_scale;

        fixed_bits[0][k] = -log2_f(1-exps[k]/32768.0f);
        fixed_bits[1][k] = 1-s*0.5f*LOG2_E - log2_f(1-(exps[k]/32768.0f)*(exps[k]/32768.0f));
    }
    while (complexity + 48 + L_frame * 11 < max_complexity)
    {
        kMax = tcx_arith_find_kMax( abs_spectrum, L_frame, scale, deadzone, deadzone_flags );
        complexity += 16 + (L_frame - kMax) * 5 + (kMax + 1) * 2;

        bits = estimator_undershoot * kMax + 1;

        for (k=0; k<=kMax; ++k)
        {
            s = envelope[k] * envelope_scale;
            q = (int)(abs_spectrum[k] * scale + deadzone);
            bits += fixed_bits[min(1,q)][k];
            bits += s*q*LOG2_E;
        }
        complexity += 32 + 6*kMax;
        if (iter == 0)
        {
            /* First rate loop iteration */
            if (scale < scale_max)
            {
                /* Only update in non-degenerate case */
                /* Update estimator temporal compensation factor */
                *target_bits_fac *= target_bits / (float)bits;
                if (*target_bits_fac > 1.25f) *target_bits_fac = 1.25f;
                if (*target_bits_fac < 0.75f) *target_bits_fac = 0.75f;
            }
        }
        if (bits <= target_bits)
        {
            /* Bits leftover => scale is too small */
            if (flag <= 0 || scale >= scale_best)
            {
                scale_best = scale;
                flag       = 1;
            }

            lob = scale;
            if (hib > 0)
            {
                /* Bisection search */
                scale = (lob + hib)*0.5f;
            }
            else
            {
                /* Initial scale adaptation */
                adjust = 1.25f * target_bits / (float)bits;
                if (adjust > 2.0f) adjust = 2.0f;
                scale *= adjust;
                if (scale > scale_max) scale = scale_max;
            }
        }
        else
        {
            /* Ran out of bits => scale is too large */
            hib = scale;
            if (lob > 0)
            {
                /* Bisection search */
                scale = (lob + hib)*0.5f;
            }
            else
            {
                /* Initial scale adaptation */
                adjust = 0.8f * target_bits / (float)bits;
                if (adjust < 0.5f) adjust = 0.5f;
                scale *= adjust;
            }
            if (flag <= 0)
            {
                scale_best = scale;
                flag = 0;
            }
        }
        ++iter;
    }


    return scale_best;
}


/*-------------------------------------------------------------------*
 * tcx_arith_encode()
 *
 *
 *-------------------------------------------------------------------*/

static int tcx_arith_encode(     /* o: number of bits consumed */
    int q_abs_spectrum[],          /* i/o: scalar quantized absolute spectrum     */
    const int signs[],             /* i: signs                                    */
    int kMax,                      /* i: number of nonzero spectral lines to code */
    int L_frame,                   /* i: nominal number of spectral lines         */
    const Word16 exps[],           /* i: expfp(-(int)envelope[]/2)                */
    int target_bits,               /* i: target bit budget                        */
    int prm[]                      /* o: bit-stream                               */
)
{
    Tastat as, as_lastgood;
    int bp, bp_lastgood;
    int k;
    int kEncoded;
    Word16 tmpi1, tmpi2;

    /* Final coding */
    ari_start_encoding_14bits(&as);
    ari_copy_states(&as, &as_lastgood);
    bp = bp_lastgood = 0;
    kEncoded = kMax;
    for (k=0; k<=kMax; ++k)
    {
        if (q_abs_spectrum[k] == 0)
        {
            assert(exps[k] >= 2);
            bp = ari_encode_14bits_range(prm, bp, target_bits, &as, exps[k]>>1, 16384);
        }
        else
        {
            /* q_abs_spectrum[k] != 0 */
            powfp_odd2(exps[k], q_abs_spectrum[k], &tmpi1, &tmpi2);
            while (tmpi1 < tmpi2 + 2)
            {
                --q_abs_spectrum[k];
                powfp_odd2(exps[k], q_abs_spectrum[k], &tmpi1, &tmpi2);
            }
            bp = ari_encode_14bits_range(prm, bp, target_bits, &as, tmpi2>>1, tmpi1>>1);
            bp = ari_encode_14bits_sign(prm, bp, target_bits, &as, signs[k]);
        }
        /* Check bit budget status */
        if (as.high <= as.low)
        {
            /* no bits left */
            /* printf("\noverflow at %d\n\n", k); */
            if (q_abs_spectrum[k] > 1)   /* Lower magnitude is still > 0 */
            {
                /* Restore state */
                ari_copy_states(&as_lastgood, &as);
                bp = bp_lastgood;

                /* Quantize to lower magnitude */
                --q_abs_spectrum[k];

                /* Retry encoding */
                powfp_odd2(exps[k], q_abs_spectrum[k], &tmpi1, &tmpi2);
                bp = ari_encode_14bits_range(prm, bp, target_bits, &as, tmpi2>>1, tmpi1>>1);
                bp = ari_encode_14bits_sign(prm, bp, target_bits, &as, signs[k]);
                if (as.high > as.low)   /* Success */
                {
                    ari_copy_states(&as, &as_lastgood);
                    bp_lastgood = bp;
                    kEncoded = k;
                    for (++k; k <= kMax; k++)
                    {
                        q_abs_spectrum[k] = 0;
                    }
                    break;
                }
            }
            ari_copy_states(&as_lastgood, &as);
            bp = bp_lastgood;
            kEncoded = k-1;
            for (; k <= kMax; k++)
            {
                q_abs_spectrum[k] = 0;
            }
            break;
        }
        else
        {
            ari_copy_states(&as, &as_lastgood);
            bp_lastgood = bp;
        }
    }

    /* Send zeros until L_frame */
    for (k=kEncoded+1, kEncoded=L_frame-1; k<L_frame; ++k)
    {
        assert(exps[k] >= 2);
        bp = ari_encode_14bits_range(prm, bp, target_bits, &as, exps[k]>>1, 16384);
        /* Check bit budget status */
        if (as.high <= as.low)
        {
            /* no bits left */
            ari_copy_states(&as_lastgood, &as);
            bp = bp_lastgood;
            kEncoded = k-1;
            break;
        }
        else
        {
            ari_copy_states(&as, &as_lastgood);
            bp_lastgood = bp;
        }
    }

    if (kEncoded == L_frame-1)
    {
        /* RESQ bits possibly available */
        /* Limit target bits to actually needed bits */
        bp = ari_done_cbr_encoding_14bits(prm, bp, bp + 16 + as.vobf, &as);
    }
    else
    {
        bp = ari_done_cbr_encoding_14bits(prm, bp, target_bits, &as);
    }

    return bp;
}


/*-------------------------------------------------------------------*
 * tcx_arith_encode_envelope()
 *
 *
 *-------------------------------------------------------------------*/

void tcx_arith_encode_envelope(
    float spectrum[],                       /* i/o: MDCT coefficients           */
    int signs[],                            /* o: signs (spectrum[.]<0)         */
    int L_frame,                            /* i: frame or MDCT length          */
    int L_spec,                             /* i: length w/o BW limitation      */
    Encoder_State *st,                   /* i/o: coder state                 */
    const Word16 A_ind[],                   /* i: quantised LPC coefficients    */
    int target_bits,                        /* i: number of available bits      */
    int prm[],                              /* o: bitstream parameters          */
    int use_hm,                             /* i: use HM in current frame?      */
    int prm_hm[],                           /* o: HM parameter area             */
    short tcxltp_pitch,                     /* i: TCX LTP pitch in FD, -1 if n/a*/
    int *arith_bits,                        /* o: bits used for ari. coding     */
    int *signaling_bits,                    /* o: bits used for signaling       */
    int low_complexity                      /* i: low-complexity flag           */
)
{
    Word16 tmp;
    Word32 env[N_MAX_ARI];             /* unscaled envelope (Q16) */
    Word16 *envelope; /* scaled envelope (Q15-e) */
    Word16 envelope_e;
    Word16 exponents[N_MAX_ARI]; /* Q15 */
    int L_spec_core;
    int *q_spectrum;
    TCX_config *tcx_cfg;
    float scale;
    int k, kMax;
    float deadzone;
    const int *deadzone_flags;
    float gamma_w, gamma_uw;
    int hm_bits;

    assert(L_spec<=N_MAX_ARI);

    tcx_cfg        = &st->tcx_cfg;
    deadzone       = tcx_cfg->sq_rounding;
    deadzone_flags = st->memQuantZeros;
    *signaling_bits = 0;
    assert(st->enableTcxLpc);
    gamma_w  = 1.0f;
    gamma_uw = 1.0f/st->gamma;
    tcx_arith_render_envelope( A_ind, L_frame, L_spec, FL2WORD16(tcx_cfg->preemph_fac), FL2WORD16(gamma_w), FL2WORD16(0.5f*gamma_uw), env );

    for (k=0; k<L_spec; ++k)
    {
        if(spectrum[k] < 0)
        {
            spectrum[k]=-spectrum[k];
            signs[k]=1;
        }
        else
        {
            signs[k]=0;
        }
    }

    if (use_hm)
    {
        tcx_hm_analyse( spectrum, L_spec, env, target_bits, tcx_cfg->coder_type, prm_hm, tcxltp_pitch, st->tcxltp_gain, &hm_bits );

        target_bits     -= hm_bits;
        *signaling_bits += hm_bits;
    }
    else
    {
        prm_hm[0] = 0;  /* just to be sure */
        hm_bits   = 0;
    }

    L_spec_core = L_spec;
    if (st->igf)
    {
        L_spec_core = min(L_spec_core, st->hIGFEnc.infoStartLine);
    }
    envelope = (Word16*)env;
    tcx_arith_scale_envelope( L_spec, L_spec_core, env, target_bits, low_complexity, envelope, &envelope_e );

    tmp = sub(envelope_e, 1);
    FOR (k = 0; k < L_spec; k++)
    {
        exponents[k] = expfp(negate(envelope[k]), tmp);
    }
    scale = tcx_arith_rateloop( spectrum, L_spec, envelope, envelope_e, exponents, target_bits, deadzone, deadzone_flags, &st->LPDmem.tcx_target_bits_fac );

    /* Final quantization */
    kMax = tcx_arith_find_kMax( spectrum, L_spec, scale, deadzone, deadzone_flags );

    q_spectrum     = (int*)env; /* Reuse buffer */
    for (k=0; k<=kMax; ++k)
    {
        /* quantise using dead-zone */
        q_spectrum[k] = (int)(spectrum[k] * scale + deadzone);
    }

    /* Final encoding */
    *arith_bits = tcx_arith_encode( q_spectrum, signs, kMax, L_spec, exponents, target_bits, prm );

    /* Multiply back the signs */
    for (k=0; k<=kMax; ++k)
    {
        spectrum[k] = (float)(q_spectrum[k] * (1-2*signs[k]));
    }
    for (; k<max(L_frame, L_spec); ++k)
    {
        spectrum[k] = 0;
    }

    return;
}
