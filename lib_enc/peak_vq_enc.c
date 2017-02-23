/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <math.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"

/*--------------------------------------------------------------------------
 * Local functions
 *--------------------------------------------------------------------------*/

void quant_peaks( Encoder_State *st, const float *vect_in, float *vect_out, const float *peak_gain, short *vq_idx ,const short overlap,
                  const short brate, const short Npeaks );

static short hvq_code_pos( Encoder_State *st, const short *inp, const short length, const short num_peaks );
static short sparse_code_pos( const short *inp, const short length, short *result );

/*--------------------------------------------------------------------------
 * peak_vq_enc()
 *
 * Vector Quantization of MDCT peaks
 *--------------------------------------------------------------------------*/

short peak_vq_enc(
    Encoder_State *st,                    /* i/o: encoder state structure     */
    const float *coefs,                 /* i  : Input coefficient vector    */
    float *coefs_out,             /* o  : Quantized output vector     */
    const short brate,                  /* i  : Core bitrate                */
    const short num_bits,               /* i  : Number of bits for HVQ      */
    const short vq_peaks,               /* i  : Number of identified peaks  */
    const short *ynrm,                  /* i  : Envelope coefficients       */
    short *R,                     /* i/o: Bit allocation/updated bit allocation */
    short *vq_peak_idx,           /* i  : Peak index vector           */
    float *nf_gains               /* i  : Estimated noise floor gains */
)
{
    short pos_bits;
    float normq;
    float pgain_q[HVQ_MAX_PEAKS_32k];
    float peak_gains[HVQ_MAX_PEAKS_32k];
    float coefs_pvq[HVQ_PVQ_BUF_LEN];
    float pvq_vector[HVQ_PVQ_BUF_LEN];
    float *pPvqVectorBandStart;
    float fg_pred[NB_SFM_MAX];
    short i, j, k, m, r, pvq_bands, num_overlap_bins;
    short hcode_l, FlagN, low_peak_bin, vq_cb_idx, max_peaks, bin_th;
    short bits = 0;
    short q_nf_gain_idx[HVQ_NF_GROUPS];
    short nf_seed = RANDOM_INITSEED;
    short pgain_cb_idx[HVQ_MAX_PEAKS], pgain_difidx[HVQ_MAX_PEAKS];
    short pvq_norm[MAX_PVQ_BANDS];
    short pvq_bits, bit_budget;
    short pos_vec[HVQ_THRES_BIN_32k];
    short npulses[MAX_PVQ_BANDS];
    short pvq_inp_vector[HVQ_PVQ_BUF_LEN];
    short k_sort[MAX_PVQ_BANDS];
    short Rk[MAX_PVQ_BANDS];
    short Rk_f[MAX_PVQ_BANDS]; /*Q3*/
    float gopt[NB_SFM];
    short sel_bnds[HVQ_NUM_SFM_24k];
    short n_sel_bnds;
    int   manE_peak, manPkEnrg;
    short expE_peak, expPkEnrg;
    short hvq_band_end[MAX_PVQ_BANDS];
    short hvq_band_start[MAX_PVQ_BANDS];
    short hvq_band_width[MAX_PVQ_BANDS];
    short n;
    short s;
    set_f( coefs_pvq, 0.0f, HVQ_PVQ_BUF_LEN );
    set_f( pvq_vector, 0.0f, HVQ_PVQ_BUF_LEN );
    set_s( npulses, 0, MAX_PVQ_BANDS );

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

    for (i = 0; i < bin_th; i++)
    {
        pos_vec[i] = 0;
    }

    /* Quantize noise floor gains */
    for (i = 0; i < HVQ_NF_GROUPS; i++)
    {
        logqnorm(&nf_gains[i], &q_nf_gain_idx[i], 32, 1, &thren[0]);
        nf_gains[i] = 0.5f*dicn[q_nf_gain_idx[i]];
        push_indice( st, IND_HVQ_NF_GAIN , q_nf_gain_idx[i], 5);
        bits += 5;
    }

    /* Signal number of peaks */
    i = max_peaks - vq_peaks;
    push_indice( st, IND_NUM_PEAKS, i, 5);
    bits += 5;

    /* Identify position of first peak and arrange peak gains by position */
    low_peak_bin = bin_th;
    for (i = 0; i < vq_peaks; i++)
    {
        if (vq_peak_idx[i] < low_peak_bin)
        {
            low_peak_bin = vq_peak_idx[i];
        }
        pos_vec[vq_peak_idx[i]] = (short) sign((float) coefs[vq_peak_idx[i]]);
    }

    for (i = 0, j = 0; i < bin_th; i++)
    {
        if(pos_vec[i] != 0)
        {
            peak_gains[j] = (float) fabs(coefs[i]);
            vq_peak_idx[j] = i;
            j++;
        }
    }

    /* Scale down peak gains */
    for (i = 0; i < vq_peaks; i++)
    {
        peak_gains[i] *= 0.25f;
    }

    /* Quantize peak gains */
    logqnorm(&peak_gains[0],&pgain_cb_idx[0],32,1,&thren_pg[0]);
    for (i = 1; i < vq_peaks; i++)
    {
        logqnorm(&peak_gains[i],&pgain_cb_idx[i],45,1,&thren_pg[0]);
    }

    /* Code quantized peak gain indices */
    diffcod(vq_peaks, pgain_cb_idx, &pgain_difidx[1]);
    for(i = 0; i < vq_peaks; i++)
    {
        pgain_q[i] = dicn_pg[pgain_cb_idx[i]];
    }
    pgain_difidx[0] = pgain_cb_idx[0];

    /* Scale up peak gains and accumulate peak energy */
    manE_peak = 0;
    expE_peak = 32;
    for (i = 0; i < vq_peaks; i++)
    {
        pgain_q[i] *= 4.0f;
        manPkEnrg = manPkEnrg_tbl[pgain_cb_idx[i]];
        expPkEnrg = expPkEnrg_tbl[pgain_cb_idx[i]];
        floating_point_add(&manE_peak, &expE_peak, manPkEnrg, expPkEnrg);
    }

    /* Huffman coding */
    hcode_l = 0;
    for (i = 1; i < vq_peaks; i++)
    {
        hcode_l += pgain_huffsizn[pgain_difidx[i]];
    }

    FlagN = HUFCODE;

    if (hcode_l >= GAINI_BITS * (vq_peaks - 1))
    {
        hcode_l = GAINI_BITS * (vq_peaks - 1);
        FlagN = NOHUFCODE;
    }

    push_indice( st, IND_FLAGN, FlagN, 1);
    push_indice( st, IND_PG_IDX, pgain_difidx[0], GAIN0_BITS);

    if (FlagN)
    {
        for (i = 1; i < vq_peaks; i++)
        {
            j = pgain_difidx[i];
            m = pgain_huffnorm[j];
            r = pgain_huffsizn[j];

            push_indice( st, IND_PG_IDX, m, r );
        }
    }
    else
    {
        for (i = 1; i < vq_peaks; i++)
        {
            push_indice( st, IND_PG_IDX, pgain_difidx[i], GAINI_BITS );
        }
    }

    /* Number of bits used for peak gain quantization */
    bits += FLAGN_BITS + GAIN0_BITS + hcode_l;

    /* Add sign for peak shape normalization */
    for (i = 0; i < vq_peaks; i++)
    {
        peak_gains[i] = pos_vec[vq_peak_idx[i]] * pgain_q[i];
    }

    /* Quantize peak shapes */
    for (i = 0; i < vq_peaks-1; i++)
    {
        num_overlap_bins = 5-(vq_peak_idx[i+1]-vq_peak_idx[i]);
        quant_peaks( st, &coefs[vq_peak_idx[i]-2], &coefs_out[vq_peak_idx[i]-2], &peak_gains[i], &vq_cb_idx, num_overlap_bins, brate, vq_peaks );
        push_indice( st, IND_HVQ_PEAKS, vq_cb_idx, 8 );
        bits += 9;
    }

    quant_peaks( st, &coefs[vq_peak_idx[i]-2], &coefs_out[vq_peak_idx[i]-2], &peak_gains[i], &vq_cb_idx, 0, brate, vq_peaks );
    push_indice( st, IND_HVQ_PEAKS, vq_cb_idx, 8 );
    bits += 9;

    pos_bits = hvq_code_pos( st, pos_vec, bin_th, vq_peaks );

    bits += pos_bits;
    bit_budget = num_bits - bits;

    /* Calculate number of PVQ bands to code and assign bits */
    pvq_bands = hvq_pvq_bitalloc(bit_budget, brate, st->bwidth, ynrm, manE_peak, expE_peak, Rk, R, sel_bnds, &n_sel_bnds );

    /* Get band limits for concatenated PVQ target */
    hvq_concat_bands( pvq_bands, sel_bnds, n_sel_bnds, hvq_band_start, hvq_band_width, hvq_band_end );

    /* Quantize PVQ bands */
    i = 0;
    n = 0;
    s = 0;
    for (k = 0; k < pvq_bands; k++)
    {
        if( k >= pvq_bands - n_sel_bnds)
        {
            i = band_start_harm[sel_bnds[s]];
            s++;
        }
        k_sort[k] = (int) k;
        j = 0;
        pPvqVectorBandStart = &pvq_vector[n];
        while ( j < hvq_band_width[k] )
        {
            if (coefs_out[i] == 0)
            {
                pvq_vector[n] = coefs[i];
                j++;
                n++;
            }
            i++;
        }
        logqnorm( pPvqVectorBandStart, &pvq_norm[k], 40, hvq_band_width[k], &thren[0] );
    }

    normalizecoefs( pvq_vector, pvq_norm, pvq_bands, hvq_band_start, hvq_band_end );

    bit_budget -= HVQ_PVQ_GAIN_BITS*pvq_bands;
    for (k = 0; k < pvq_bands; k++)
    {
        Rk_f[k] = Rk[k] * 8;
    }

    pvq_bits = bit_budget;
    set_s( npulses, 0, MAX_PVQ_BANDS );

    pvq_encode_frame( st, pvq_vector, coefs_pvq, gopt, npulses, pvq_inp_vector, hvq_band_start,
                      hvq_band_end, hvq_band_width, pvq_bands, Rk_f, pvq_bits, HQ_CORE );

    for(i=0; i<pvq_bands; i++)
    {
        k_sort[i] = i;
    }


    fine_gain_pred( hvq_band_start, hvq_band_end, hvq_band_width, k_sort, npulses, NULL, NULL,
                    pvq_bands, coefs_pvq, pvq_inp_vector, fg_pred, HQ_CORE );

    i = 0;
    n = 0;
    s = 0;
    for (k = 0; k < pvq_bands; k++)
    {
        normq = dicn[pvq_norm[k]] * (gopt[k] / fg_pred[k]);

        logqnorm(&normq, &pvq_norm[k], 40, 1, &thren[0]);
        pvq_norm[k] -= 8;
        if (pvq_norm[k] < 0)
        {
            pvq_norm[k] = 0;
        }

        push_indice( st, IND_HVQ_PVQ_GAIN, pvq_norm[k], HVQ_PVQ_GAIN_BITS);
        pvq_bits += HVQ_PVQ_GAIN_BITS;

        pvq_norm[k] += 8;
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
                coefs_out[i] = coefs_pvq[n]*fg_pred[k];
                coefs_out[i] = coefs_out[i]*normq;
                j++;
                n++;
            }
            i++;
        }

    }

    bits += pvq_bits;

    /* Noise fill unquantized coeffs with one gain per group */
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

    return bits;
}

/*--------------------------------------------------------------------------
 * quant_peaks()
 *
 * Applies VQ on input vector
 *--------------------------------------------------------------------------*/

void quant_peaks(
    Encoder_State *st,              /* i/o: encoder state structure         */
    const float *vect_in,         /* i  : Target vector                   */
    float *vect_out,        /* i/o: Quantized vector                */
    const float *peak_gain,       /* i  : Peak gain vector                */
    short *vq_idx,          /* o  : Codebook index                  */
    const short overlap,          /* i  : Overlap indicator               */
    const short brate,            /* i  : Core bitrate                    */
    const short Npeaks            /* i  : Number of peaks                 */
)
{
    float x[4];
    float xq[4];
    short weights[4];
    short i, cb_class, search_overlap;

    set_s(weights,1,4);

    x[0] = vect_in[0] / (*peak_gain);
    x[1] = vect_in[1] / (*peak_gain);
    x[2] = vect_in[3] / (*peak_gain);
    x[3] = vect_in[4] / (*peak_gain);

    if(vect_out[0] != 0 )
    {
        if( fabs(peak_gain[-1]) > fabs(*peak_gain))
        {
            weights[0] = 0;

            if (vect_out[1] != 0)
            {
                weights[1] = 0;
            }
        }
    }
    if(overlap > 0)
    {
        if( fabs(peak_gain[1]) > fabs(*peak_gain) )
        {
            for(i = 3; i> 3-overlap; i--)
            {
                weights[i] = 0;
            }
        }
    }


    /* Classify */
    cb_class = (short) w_vquant(x, 0, weights, 0, hvq_class_c, HVQ_VQ_DIM-1, HVQ_NUM_CLASS, 0);
    search_overlap = (brate == HQ_24k40) ? hvq_cb_search_overlap24k[HVQ_MAX_PEAKS_24k-Npeaks] : hvq_cb_search_overlap32k[HVQ_MAX_PEAKS_32k-Npeaks];

    /* Quantize */
    if( cb_class == 0 )
    {
        *vq_idx = (short) w_vquant(x, 0, weights, xq, hvq_peak_cb, 4, HVQ_CB_SIZE/2+search_overlap, 0);
        push_indice( st, IND_HVQ_PEAKS, 0, 1 );
    }
    else if( cb_class == 1 )
    {
        *vq_idx = (short) w_vquant(x, 0, weights, xq, &hvq_peak_cb[HVQ_CB_SIZE*2-search_overlap*4], 4, HVQ_CB_SIZE/2+search_overlap, 0);
        *vq_idx += HVQ_CB_SIZE/2-search_overlap;
        push_indice( st, IND_HVQ_PEAKS, 0, 1 );
    }
    else if( cb_class == 2 )
    {
        *vq_idx = (short) w_vquant(x, 0, weights, xq, &hvq_peak_cb[HVQ_CB_SIZE*2-search_overlap*4], 4, HVQ_CB_SIZE/2+search_overlap, 1);
        *vq_idx += HVQ_CB_SIZE/2-search_overlap;
        push_indice( st, IND_HVQ_PEAKS, 1, 1 );
    }
    else
    {
        *vq_idx = (short) w_vquant(x, 0, weights, xq, hvq_peak_cb, 4, HVQ_CB_SIZE/2+search_overlap, 1);
        push_indice( st, IND_HVQ_PEAKS, 1, 1 );
    }

    vect_out[0] = weights[0] * (xq[0] * (*peak_gain)) + (weights[0]^1)*vect_out[0];
    vect_out[1] = weights[1] * (xq[1] * (*peak_gain)) + (weights[1]^1)*vect_out[1];
    vect_out[2] = *peak_gain;
    vect_out[3] = weights[2] * (xq[2] * (*peak_gain)) + (weights[2]^1)*vect_out[3];
    vect_out[4] = weights[3] * (xq[3] * (*peak_gain)) + (weights[3]^1)*vect_out[4];

    return;
}

/*--------------------------------------------------------------------------
 * code_pos()
 *
 * Code pulse positions
 *--------------------------------------------------------------------------*/

static short sparse_code_pos(
    const short *inp,
    const short length,
    short *result
)
{
    short layer2[HVQ_CP_L2_MAX];
    short layer_length;
    short i,j;
    short val, idx;
    short bits = 0;
    short mask;

    set_s(layer2, 0, HVQ_CP_L2_MAX);

    layer_length = (short)((float)length/HVQ_CP_L1_LEN + 0.5);

    for (j = 0; j < layer_length; j++)
    {
        for (i = j*HVQ_CP_L1_LEN; i < min((j+1)*HVQ_CP_L1_LEN, length); i++)
        {
            if (inp[i])
            {
                layer2[j] = 1;
                break;
            }
        }
    }

    for (i = 0; i < layer_length; i++)
    {
        result[i] = layer2[i];
    }
    bits += layer_length;

    for (j = 0; j < layer_length; j++)
    {
        if (layer2[j])
        {
            val = 0;
            for (i = j*HVQ_CP_L1_LEN; i < min((j+1)*HVQ_CP_L1_LEN, length); i++)
            {
                val <<= 1;
                val |= inp[i];
            }

            for (idx = 0; idx < HVQ_CP_MAP_LEN; idx++)
            {
                if (hvq_cp_layer1_map5[idx] == val)
                {
                    break;
                }
            }

            mask = 1<<(HVQ_CP_MAP_IDX_LEN - 1);
            for (i = 0; i < HVQ_CP_MAP_IDX_LEN; i++)
            {
                result[bits++] = (idx&mask) >> (HVQ_CP_MAP_IDX_LEN - 1 - i);
                mask >>= 1;
            }
        }
    }

    return bits;
}

/*--------------------------------------------------------------------------
 * hvq_code_pos()
 *
 * Code pulse positions
 *--------------------------------------------------------------------------*/

static short hvq_code_pos(
    Encoder_State *st,              /* i/o: encoder state structure      */
    const short *inp,
    const short length,
    const short num_peaks
)
{
    short sparse_result[4*HVQ_THRES_BIN_32k/HVQ_CP_L1_LEN];
    short delta[HVQ_MAX_PEAKS_32k];
    short peak_idx[HVQ_MAX_PEAKS_32k];
    short inp_abs[HVQ_THRES_BIN_32k];
    short inp_sign[HVQ_MAX_PEAKS_32k];

    short i, j;
    short bits;
    short delta_max;
    short delta_bits, sparse_bits;

    bits = 0;

    /* Extract sorted peak index vector and sign vector */
    for (i = 0, j = 0; i < length; i++)
    {
        inp_abs[i] = (short)abs(inp[i]);
        if (inp[i])
        {
            peak_idx[j] = i;
            inp_sign[j++] = inp[i];
        }
    }

    /* Calculate delta */
    delta[0] = peak_idx[0] + HVQ_CP_HUFF_OFFSET;
    delta_max = delta[0];
    for (i = 1; i < num_peaks; i++)
    {
        delta[i] = peak_idx[i] - peak_idx[i-1] - HVQ_CP_HUFF_OFFSET;
        if (delta_max < delta[i])
        {
            delta_max = delta[i];
        }
    }

    /* Calculate bits needed for huffman coding of deltas */
    delta_bits = -1;
    if (delta_max <= HVQ_CP_HUFF_MAX)
    {
        delta_bits = 0;
        for (i = 0; i < num_peaks; i++)
        {
            delta_bits += hvq_cp_huff_len[delta[i]];
        }
    }

    /* Calculate bits neeed for sparse coding */
    sparse_bits = sparse_code_pos(inp_abs, length, sparse_result);

    /* Decide which coding mode to use */
    if (delta_bits > sparse_bits || delta_bits < 0)
    {
        push_indice(st, IND_POS_IDX, HVQ_CP_SPARSE, 1);

        for (i = 0; i < sparse_bits; i++)
        {
            push_indice(st, IND_POS_IDX, sparse_result[i], 1);
        }
        bits += sparse_bits + 1;
    }
    else
    {
        push_indice(st, IND_POS_IDX, HVQ_CP_DELTA, 1);

        for (i = 0; i < num_peaks; i++)
        {
            j = delta[i];
            push_indice(st, IND_POS_IDX, hvq_cp_huff_val[j], hvq_cp_huff_len[j]);
        }
        bits += delta_bits + 1;
    }

    /* Send sign */
    for (i = 0; i < num_peaks; i++)
    {
        push_indice(st, IND_POS_IDX, (inp_sign[i] < 0 ? 0 : 1), 1);
    }
    bits += num_peaks;

    return bits;
}

