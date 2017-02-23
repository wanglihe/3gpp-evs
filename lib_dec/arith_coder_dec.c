/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <assert.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"
#include "basop_util.h"
#include "basop_proto_func.h"


/*-------------------------------------------------------*
 * tcx_arith_decode()
 *
 *
 *-------------------------------------------------------*/

static int tcx_arith_decode(     /* o: number of bits consumed                */
    int L_frame,                   /* i: number of spectral lines               */
    const Word16 envelope[],       /* i: scaled envelope (Q15-e)                */
    Word16 envelope_e,             /* i: scaled envelope exponent (Q0)          */
    int target_bits,               /* i: target bit budget                      */
    const int prm[],               /* i: bit-stream                             */
    float q_spectrum[]             /* o: scalar quantized spectrum              */
)
{
    int bp, k, q, s;
    Tastat as;
    unsigned int exp_k;
    Word16 tmp;


    bp = ari_start_decoding_14bits_prm(prm, 0, &as);
    tmp = sub(envelope_e, 1);

    for (k=0; k<L_frame; ++k)
    {
        if( envelope[k] == 0 )  /* safety check in case of bit errors */
        {
            set_zero(q_spectrum, L_frame);
            return -1;
        }
        else
        {
            exp_k = expfp(negate(envelope[k]), tmp);
        }
        /* decode line magnitude */
        bp = ari_decode_14bits_pow(prm, bp, target_bits, &q, &as, exp_k);
        if (q)
        {
            /* line is non-zero, decode sign */
            bp = ari_decode_14bits_sign(prm, bp, target_bits, &s, &as);
            q_spectrum[k] = (float)q * (3 - 2*s);
        }
        else
        {
            /* line is zero, no sign needed */
            q_spectrum[k] = 0.0f;
        }

        if (as.high <= as.low)
        {
            if( bp < target_bits )  /* safety check in case of bit errors */
            {
                bp = -1;
            }
            break;    /* no bits left, so exit loop */
        }
    }
    for (; k < L_frame; k++)
    {
        q_spectrum[k] = 0.0f;
    }

    return bp;
}


/*-------------------------------------------------------*
 * tcx_arith_decode_envelope()
 *
 *
 *-------------------------------------------------------*/

void tcx_arith_decode_envelope(
    float q_spectrum[],               /* o: quantised MDCT coefficients   */
    int L_frame,                      /* i: frame or MDCT length          */
    int L_spec,                       /* i: length w/o BW limitation      */
    Decoder_State *st,                /* i/o: coder state                 */
    const short coder_type,           /* i  : coder type                  */
    const Word16 A_ind[],             /* i: quantised LPC coefficients    */
    float tcxltp_gain,                /* i: TCX LTP gain                  */
    int target_bits,                  /* i: number of available bits      */
    const int prm[],                  /* i: bitstream parameters          */
    int use_hm,                       /* i: use HM in current frame?      */
    const int prm_hm[],               /* i: HM parameter area             */
    short tcxltp_pitch,               /* i: TCX LTP pitch in FD, -1 if n/a*/
    int *arith_bits,                  /* o: bits used for ari. coding     */
    int *signaling_bits,              /* o: bits used for signaling       */
    int low_complexity                /* i: low-complexity flag           */
)
{
    Word32 env[N_MAX_ARI];             /* unscaled envelope (Q16) */
    Word16 *envelope; /* scaled envelope (Q15-e) */
    Word16 envelope_e;
    int L_spec_core;
    TCX_config *tcx_cfg;
    int k;
    float gamma_w, gamma_uw;
    int hm_bits;

    if( L_spec > N_MAX_ARI || (target_bits > (ACELP_13k20/50)) || (target_bits <= 0) )
    {
        /* this could happen in case of bit errors */
        st->BER_detect = 1;
        L_spec = N_MAX_ARI;
        *signaling_bits = 0;
        *arith_bits = 0;
        set_zero( q_spectrum, L_frame );

        return;
    }

    tcx_cfg = &st->tcx_cfg;
    *signaling_bits = 0;

    assert(st->enableTcxLpc);
    gamma_w  = 1.0f;
    gamma_uw = 1.0f/st->gamma;

    tcx_arith_render_envelope( A_ind, L_frame, L_spec, FL2WORD16(tcx_cfg->preemph_fac), FL2WORD16(gamma_w), FL2WORD16(0.5f*gamma_uw), env );

    if (use_hm)
    {
        if (prm_hm[0])
        {
            tcx_hm_decode( L_spec, env, target_bits, coder_type, prm_hm, tcxltp_pitch, tcxltp_gain, &hm_bits );

            if (hm_bits < 0)
            {
                st->BER_detect = 1;
                *signaling_bits = 0;
                *arith_bits = 0;
                set_zero( q_spectrum, L_frame );

                return;
            }
        }
        else
        {
            hm_bits = 1;
        }
        *signaling_bits += hm_bits;
    }

    L_spec_core = L_spec;
    if (st->igf)
    {
        L_spec_core = min(L_spec_core, st->hIGFDec.infoIGFStartLine);
    }
    envelope = (Word16*)env;
    tcx_arith_scale_envelope( L_spec, L_spec_core, env, target_bits, low_complexity, envelope, &envelope_e );

    *arith_bits = tcx_arith_decode( L_spec, envelope, envelope_e, target_bits, prm, q_spectrum );

    /* safety check in case of bit errors */
    if( *arith_bits < 0 )
    {
        st->BER_detect = 1;
        set_zero( q_spectrum, L_frame );
    }

    for (k=L_spec; k<L_frame; ++k)
    {
        q_spectrum[k] = 0;
    }

    return;
}
