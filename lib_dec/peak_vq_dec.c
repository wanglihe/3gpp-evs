/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"

/*------------------------------------------------------------------------*
 * Local functions
 *------------------------------------------------------------------------*/

static void dequant_peaks( Decoder_State *st, float *vect_out, const float *peak_gain );

static short hvq_dec_pos( Decoder_State *st, short *pos_vec, const short length, const short num_peaks );
static short sparse_dec_pos(Decoder_State *st, short *out, const short length );

/*--------------------------------------------------------------------------
 * hvq_dec()
 *
 * HVQ decoder
 *--------------------------------------------------------------------------*/

void hvq_dec(
    Decoder_State *st,                /* i/o: decoder state structure       */
    const short num_bits,           /* i : Number of available bits       */
    const long  core_brate,         /* i : Core bit-rate                  */
    const short *ynrm,              /* i : Envelope coefficients          */
    short *R,                 /* i/o: Bit allocation/updated bit allocation */
    float *noise_level,       /* o : Noise level                    */
    short *peak_idx,          /* o : Peak position vector           */
    short *Npeaks,            /* o : Total number of peaks          */
    float *coefsq_norm,       /* o : Output vector                  */
    const short core
)
{
    short i;
    short bits;
    short noise_level_idx;

    bits = num_bits;

    for( i = 0; i < HVQ_BWE_NOISE_BANDS; i++ )
    {
        noise_level_idx = (short) get_next_indice( st, 2 );
        noise_level[i] = usdequant( noise_level_idx, 0.0f, 0.1f );

        bits -= 2;
    }

    peak_vq_dec( st, coefsq_norm, (short)core_brate, bits, ynrm, R, peak_idx,
                 Npeaks, core );
}

/*--------------------------------------------------------------------------
 * peak_vq_dec()
 *
 * Vector de-quantization of MDCT peaks
 *--------------------------------------------------------------------------*/

void peak_vq_dec(
    Decoder_State *st,                    /* i/o: decoder state structure     */
    float *coefs_out,             /* o  : Output coefficient vector   */
    const short brate,                  /* i  : Core bitrate                */
    const short num_bits,               /* i  : Number of bits for HVQ      */
    const short *ynrm,                  /* i  : Envelope coefficients       */
    short *R,                     /* i/o: Bit allocation/updated bit allocation */
    short *vq_peak_idx,           /* o  : Peak position vector        */
    short *Npeaks,                /* o  : Number of peaks             */
    const short core                    /* i  : Core type                   */
)
{
    short vq_peaks, i, j, k, FlagN, hcode_l, diff;
    short bin_th, max_peaks, pvq_bands;
    short nf_seed = RANDOM_INITSEED;
    short nf_gains_idx[HVQ_NF_GROUPS], pgain_difidx[HVQ_MAX_PEAKS_32k], pvq_norm[MAX_PVQ_BANDS];
    short gain_bits_array[MAX_PVQ_BANDS];
    short pos_bits;
    float nf_gains[HVQ_NF_GROUPS], peak_gains[HVQ_MAX_PEAKS_32k];
    int   manE_peak, manPkEnrg;  /* (man, exp) representation ported from BASOP for interoperability */
    short expE_peak, expPkEnrg;
    float pvq_vector[HVQ_PVQ_BUF_LEN];
    short res_vec[HVQ_THRES_BIN_32k];
    short k_sort[HVQ_MAX_PVQ_WORDS];
    short pvq_inp_vector[HVQ_PVQ_BUF_LEN];
    short npulses[MAX_PVQ_BANDS];
    short pvq_bits, Rk[MAX_PVQ_BANDS];
    float fg_pred[NB_SFM_MAX];
    short Rk_f[MAX_PVQ_BANDS]; /* Q3 */
    short sel_bnds[HVQ_NUM_SFM_24k];
    short n_sel_bnds;
    short hvq_band_end[MAX_PVQ_BANDS];
    short hvq_band_start[MAX_PVQ_BANDS];
    short hvq_band_width[MAX_PVQ_BANDS];
    short n;
    short s;
    float normq;

    set_s( gain_bits_array, 0, MAX_PVQ_BANDS );
    set_f( pvq_vector, 0.0f, HVQ_PVQ_BUF_LEN );
    set_s( npulses, 0, MAX_PVQ_BANDS );
    set_s( pvq_inp_vector, 0, HVQ_PVQ_BUF_LEN );

    /* Set bitrate dependent variables */
    if (brate == HQ_24k40)
    {
        max_peaks = HVQ_MAX_PEAKS_24k;
        bin_th = HVQ_THRES_BIN_24k;
    }
    else
    {
        max_peaks = HVQ_MAX_PEAKS_32k;
        bin_th = HVQ_THRES_BIN_32k;
    }

    /* Get number of peaks */
    vq_peaks = (short) get_next_indice( st, 5 );
    vq_peaks = max_peaks - vq_peaks;
    *Npeaks = vq_peaks;
    diff = 5;

    /* safety check in case of bit errors */
    if( *Npeaks < HVQ_MIN_PEAKS )
    {
        st->BER_detect = 1;
        vq_peaks = HVQ_MIN_PEAKS;
        *Npeaks = HVQ_MIN_PEAKS;
    }

    /* De-quantize peak positions */
    for (i = 0; i < bin_th; i++)
    {
        res_vec[i] = 0;
    }

    /* Unpack HVQ codewords */
    pos_bits = hvq_dec_pos(st, res_vec, bin_th, vq_peaks);
    diff += pos_bits;

    for (i = 0, j = 0; i < bin_th && j < vq_peaks; i++) /* safety check in case of bit errors */
    {
        if (res_vec[i])
        {
            vq_peak_idx[j++] = i;
        }
    }

    /* safety check in case of bit errors */
    if( j < vq_peaks )
    {
        st->BER_detect = 1;
        vq_peaks = j - 1;
        *Npeaks = j - 1;
    }

    /* Huffman or differential coding */
    FlagN = (short) get_next_indice( st, 1 );

    /* De-quantize peak gains */
    pgain_difidx[0] = (short) get_next_indice( st, GAIN0_BITS );

    /* safety check in case of bit errors */
    if( pgain_difidx[0] > 44 )
    {
        st->BER_detect = 1;
        pgain_difidx[0] = 44;
    }
    peak_gains[0] = dicn_pg[pgain_difidx[0]]*sign((float) res_vec[vq_peak_idx[0]]);

    hcode_l = 0;
    if(FlagN)
    {
        huff_dec( st, vq_peaks-1, MAX_PG_HUFFLEN, NUM_PG_HUFFLEN, hvq_pg_huff_thres, hvq_pg_huff_offset, hvq_pg_huff_tab, &pgain_difidx[1] );
        for (i = 1; i < vq_peaks; i++)
        {
            hcode_l += pgain_huffsizn[pgain_difidx[i]];
        }
    }
    else
    {
        for (i = 1; i < vq_peaks; i++)
        {
            pgain_difidx[i] = (short) get_next_indice( st, GAINI_BITS );
            hcode_l += GAINI_BITS;
        }
    }

    for (i = 1; i < vq_peaks; i++)
    {
        pgain_difidx[i] += pgain_difidx[i - 1] - 15;

        /* safety check in case of bit errors */
        if( pgain_difidx[i] > 44 || pgain_difidx[i] < 0)
        {
            st->BER_detect = 1;
            pgain_difidx[i] = 44;
        }

        peak_gains[i] = dicn_pg[pgain_difidx[i]]*sign((float) res_vec[vq_peak_idx[i]]);
    }

    /* Scale up peak gains and accumulate peak energy */
    /* Simulating BASOP code for interoperability     */
    manE_peak = 0;
    expE_peak = 32;
    for (i = 0; i < vq_peaks; i++)
    {
        peak_gains[i] *= 4.0f;
        manPkEnrg = manPkEnrg_tbl[pgain_difidx[i]];
        expPkEnrg = expPkEnrg_tbl[pgain_difidx[i]];
        floating_point_add(&manE_peak, &expE_peak, manPkEnrg, expPkEnrg);
    }

    /* Number of bits used for peak gain quantization */
    diff += FLAGN_BITS + GAIN0_BITS + hcode_l;

    /* De-quantize peaks */
    for (i = 0; i < vq_peaks; i++)
    {
        dequant_peaks( st, &coefs_out[vq_peak_idx[i]-2], &peak_gains[i] );

        diff += 9;
    }

    for (i = 0; i < HVQ_NF_GROUPS; i++)
    {
        nf_gains_idx[i] = (short) get_next_indice( st, 5 );
        nf_gains[i] = 0.5f*dicn[nf_gains_idx[i]];
        diff += 5;
    }

    pvq_bits = num_bits - diff;

    /* Calculate number of PVQ bands to code and assign bits */
    pvq_bands = hvq_pvq_bitalloc(pvq_bits, brate, st->bwidth, ynrm, manE_peak, expE_peak, Rk, R, sel_bnds, &n_sel_bnds );

    /* safety check in case of bit errors */
    if (pvq_bands == 0)
    {
        st->BER_detect = 1;
    }

    pvq_bits -= HVQ_PVQ_GAIN_BITS*pvq_bands;

    /* Get band limits for concatenated PVQ target */
    hvq_concat_bands( pvq_bands, sel_bnds, n_sel_bnds, hvq_band_start, hvq_band_width, hvq_band_end );

    s = 0;
    for (k = 0; k < pvq_bands; k++)
    {
        k_sort[k] = k;
        Rk_f[k] = Rk[k] * 8;
    }

    pvq_decode_frame( st, pvq_vector, npulses, pvq_inp_vector, hvq_band_start,
                      hvq_band_end, hvq_band_width, pvq_bands, Rk_f, pvq_bits, core );

    fine_gain_pred( hvq_band_start, hvq_band_end, hvq_band_width, k_sort, npulses, NULL, NULL,
                    pvq_bands, pvq_vector, pvq_inp_vector, fg_pred, core );

    fine_gain_dec( st, k_sort, pvq_bands, gain_bits_array, fg_pred );

    apply_gain(k_sort, hvq_band_start, hvq_band_end, pvq_bands, fg_pred, pvq_vector );

    i = 0;
    n = 0;
    s = 0;
    for (k = 0; k < pvq_bands; k++)
    {
        pvq_norm[k] = (short) get_next_indice( st, HVQ_PVQ_GAIN_BITS );
        pvq_norm[k] += 8;

        diff += HVQ_PVQ_GAIN_BITS;

        j = 0;
        if( k >= pvq_bands - n_sel_bnds)
        {
            i = band_start_harm[sel_bnds[s++]];
        }
        while (j < hvq_band_width[k])
        {
            normq = dicn[pvq_norm[k]];
            if (coefs_out[i] == 0)
            {
                coefs_out[i] = pvq_vector[n] * normq;
                j++;
                n++;
            }
            i++;
        }
    }
    /* Noise fill unqantized coeffs with one gain per group */
    for (i = 0; i < HVQ_NF_GROUPS; i++)
    {
        for (j = i*(bin_th/HVQ_NF_GROUPS); j < (i+1)*(bin_th/HVQ_NF_GROUPS); j++)
        {
            if (coefs_out[j] == 0)
            {
                coefs_out[j] = ((float)own_random(&nf_seed)/MAX16B)*nf_gains[i];
            }
        }
    }

    return;
}

/*--------------------------------------------------------------------------
 * dequant_peaks()
 *
 * Reads codebook vector and scales peak
 *--------------------------------------------------------------------------*/

static void dequant_peaks(
    Decoder_State *st,                      /* i/o: decoder state structure */
    float *vect_out,                /* o  : Quantized vector        */
    const float *peak_gain                /* i  : Peak gain               */
)
{
    float xq[4];
    const float *tmp;
    short i, hvq_cb_rev;
    short cb_idx;

    hvq_cb_rev = (short) get_next_indice( st, 1 );
    cb_idx = (short) get_next_indice( st, 8 );

    if( hvq_cb_rev )
    {
        tmp = &hvq_peak_cb[cb_idx*4+3];
        for (i = 0; i < 4; i++)
        {
            xq[i] = tmp[-i];
        }
    }
    else
    {
        mvr2r(&hvq_peak_cb[cb_idx*4], xq, 4);
    }
    if(vect_out[0] == 0)
    {
        vect_out[0] = xq[0] **peak_gain;
        vect_out[1] = xq[1] **peak_gain;
    }
    else
    {
        if( fabs(peak_gain[-1]) <= fabs(*peak_gain) )
        {
            vect_out[0] = xq[0] **peak_gain;
            vect_out[1] = xq[1] **peak_gain;
        }
        else
        {
            if(vect_out[1] == 0 || fabs(peak_gain[-1]) <= fabs(*peak_gain))
            {
                vect_out[1] = xq[1] **peak_gain;
            }
        }
    }
    vect_out[2] = *peak_gain;
    vect_out[3] = xq[2] **peak_gain;
    vect_out[4] = xq[3] **peak_gain;

    return;
}

/*--------------------------------------------------------------------------
 * hvq_dec_pos()
 *
 * HVQ decode peak positions
 *--------------------------------------------------------------------------*/

static short hvq_dec_pos(
    Decoder_State *st,                      /* i/o: decoder state structure   */
    short *pos_vec,
    const short length,
    const short num_peaks
)
{
    short peak_idx[HVQ_MAX_PEAKS_32k];
    short delta[HVQ_MAX_PEAKS_32k];
    short sign_vec[HVQ_MAX_PEAKS_32k];

    short mode;
    short num_bits;
    short i, j;

    num_bits = 0;
    set_s(pos_vec, 0, length);

    mode = (short)get_next_indice(st, 1);
    num_bits += 1;

    if (mode == HVQ_CP_DELTA)
    {
        huff_dec(st, num_peaks, HVQ_CP_HUFF_MAX_CODE, HVQ_CP_HUFF_NUM_LEN, hvq_cp_huff_thres, hvq_cp_huff_offset, hvq_cp_huff_tab, delta);

        for (i = 0; i < num_peaks; i++)
        {
            num_bits += hvq_cp_huff_len[delta[i]];
        }

        peak_idx[0] = delta[0] - HVQ_CP_HUFF_OFFSET;
        /* safety check in case of bit errors */
        if (peak_idx[0] < 2)
        {
            peak_idx[0] = 2;
            st->BER_detect = 1;
        }
        for (i = 1; i < num_peaks; i++)
        {
            peak_idx[i] = delta[i] + peak_idx[i-1] + HVQ_CP_HUFF_OFFSET;
            /* safety check in case of bit errors */
            if (peak_idx[i] >= HVQ_THRES_BIN_32k)
            {
                peak_idx[i] = HVQ_THRES_BIN_32k - 1;
                st->BER_detect = 1;
            }
        }

        for (i = 0; i < num_peaks; i++)
        {
            pos_vec[peak_idx[i]] = 1;
        }
    }
    else
    {
        num_bits += sparse_dec_pos(st, pos_vec, length);
    }

    for (i = 0; i < num_peaks; i++)
    {
        sign_vec[i] = (get_next_indice_1(st) == 0) ? -1 : 1;
    }
    num_bits += num_peaks;

    for (i = 0, j = 0; i < length && j < num_peaks; i++) /* safety check in case of bit errors */
    {
        if (pos_vec[i])
        {
            pos_vec[i] *= sign_vec[j++];
        }
    }

    return num_bits;
}

/*--------------------------------------------------------------------------
 * sparse_dec_pos()
 *
 * Sparse decode positions
 *--------------------------------------------------------------------------*/

static short sparse_dec_pos(
    Decoder_State *st,                      /* i/o: decoder state structure   */
    short *out,
    const short length
)
{
    short layer2[HVQ_CP_L2_MAX];
    short layer_length;
    short i, j;
    short bits;
    short idx, val;

    set_s(layer2, 0, HVQ_CP_L2_MAX);
    set_s(out, 0, length);
    bits = 0;

    layer_length = (short)((float)length/HVQ_CP_L1_LEN + 0.5);

    for (i = 0; i < layer_length; i++)
    {
        layer2[i] = (short)get_next_indice_1(st);
    }
    bits += layer_length;

    for (j = 0; j < layer_length; j++)
    {
        if (layer2[j])
        {
            idx = (short)get_next_indice(st, HVQ_CP_MAP_IDX_LEN);
            bits += HVQ_CP_MAP_IDX_LEN;

            val = hvq_cp_layer1_map5[idx];

            /* safety check in case of bit errors */
            if ( j == 0 && val > 4 ) /* out[0] and out[1] are invalid positions */
            {
                st->BER_detect = 1;
                val = 4;
            }
            for (i = min((j+1)*HVQ_CP_L1_LEN, length)-1; i >= j*HVQ_CP_L1_LEN; i--)
            {
                out[i] = val&1;
                val >>= 1;
            }
        }
    }

    return bits;
}
