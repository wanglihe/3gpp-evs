/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_enc.h"
#include "rom_com.h"

/*--------------------------------------------------------------------------*
 * hvq_enc()
 *
 * Harmonic VQ encoder
 *--------------------------------------------------------------------------*/

short hvq_enc(                      /* o  : Consumed bits                   */
    Encoder_State *st,                /* i/o: encoder state structure         */
    const long  brate,              /* i  : Total bit rate                  */
    const short hvq_bits,           /* i  : HVQ bit budget                  */
    const short Npeaks,             /* i  : Number of peaks                 */
    const short *ynrm,              /* i  : Envelope coefficients           */
    short *R,                 /* i/o: Bit allocation/updated bit allocation */
    short *peaks,             /* i  : Peak pos. / Encoded peak pos.   */
    float *nf_gains,          /* i/o: Noise fill gains / Quant. nf gains */
    float *noise_level,       /* o  : Quantized noise level           */
    const float *pe_gains,          /* i  : Peak gains                      */
    const float *coefs,             /* i  : spectrum coefficients           */
    float *coefs_out          /* o  : encoded spectrum coefficients   */
)
{
    short bin_th,q,j,i;
    short nf_cnt;
    short q_noise_level_idx[HVQ_BWE_NOISE_BANDS];
    float q_noise_level[HVQ_BWE_NOISE_BANDS];
    float d, nf, nf_mean, pe, pe_mean, nfpe;
    short bits_used;
    float lb_nfpe;

    bits_used = 0;

    bin_th = HVQ_THRES_BIN_32k;
    if( brate == HQ_24k40 )
    {
        bin_th = HVQ_THRES_BIN_24k;
    }

    nf = 800;
    pe = 800;
    q = bin_th;

    /* Find HB noise level */
    for( i = 0; i < HVQ_BWE_NOISE_BANDS; i++ )
    {
        nf_cnt = 0;
        nf_mean = EPSILON;
        pe_mean = EPSILON;
        for( j = 0; j < (L_FRAME32k-bin_th)/HVQ_BWE_NOISE_BANDS; j++, q++ )
        {
            d = (float)fabs(coefs[q]);

            if( d > pe )
            {
                pe = HVQ_BWE_WEIGHT2 * pe + (1 - HVQ_BWE_WEIGHT2) * d;
            }
            else
            {
                pe = HVQ_BWE_WEIGHT1 * pe + (1 - HVQ_BWE_WEIGHT1) * d;
                if( d > nf )
                {
                    nf = HVQ_BWE_WEIGHT1 * nf + (1 - HVQ_BWE_WEIGHT1) * d;
                }
                else
                {
                    nf = HVQ_BWE_WEIGHT2 * nf + (1 - HVQ_BWE_WEIGHT2) * d;
                }
                nf_mean += nf;
                nf_cnt++;
            }

            pe_mean += pe;
        }

        if( nf_cnt > 0 )
        {
            nf_mean /= nf_cnt;
        }
        pe_mean /= (L_FRAME32k-bin_th)/HVQ_BWE_NOISE_BANDS;

        nfpe = HVQ_NFPE_FACTOR*nf_mean/pe_mean;
        noise_level[i] = nfpe * nfpe * nfpe;

        q_noise_level_idx[i] = usquant( noise_level[i], &q_noise_level[i], 0.0f, 0.1f, 4 );
        push_indice( st, IND_HVQ_BWE_NL, q_noise_level_idx[i], 2 );
        bits_used += 2;

        noise_level[i] = q_noise_level[i];
    }

    for( i = 0; i < HVQ_NF_GROUPS; i ++ )
    {
        lb_nfpe = HVQ_LB_NFPE_FACTOR*nf_gains[i]/pe_gains[i];
        lb_nfpe = lb_nfpe * lb_nfpe * lb_nfpe;

        if( lb_nfpe > 0.5f )
        {
            lb_nfpe = 0.5f;
        }
        nf_gains[i] *= 2*lb_nfpe;
    }

    bits_used += peak_vq_enc( st, coefs, coefs_out, (short)brate, hvq_bits - bits_used,
                              Npeaks, ynrm, R, peaks, &nf_gains[0] );

    return bits_used;
}
