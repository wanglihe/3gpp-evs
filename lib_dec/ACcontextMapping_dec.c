/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"

/*-------------------------------------------------------------------*
* ACcontextMapping_decode2_no_mem_s17_LC()
*
* Arithmetic decoder
*-------------------------------------------------------------------*/

long ACcontextMapping_decode2_no_mem_s17_LC(    /* o: resQBits              */
    Decoder_State *st,        /* i/o: decoder state                           */
    int *x,                   /* o: decoded spectrum                          */
    long nt,                  /* i: size of spectrum                          */
    int nbbits,               /* i: bit budget                                */
    int resQMaxBits,          /* i: residual coding maximum bits              */
    CONTEXT_HM_CONFIG *hm_cfg /* i: context-based harmonic model configuration*/
)
{
    Tastat as;
    int start_bit_pos, lsbs_bit_pos, overflow_bit_pos;
    int a, b, t, a1, b1, a1_i, b1_i, k;
    int lev, pki, esc_nb;
    int rateFlag;
    int r, lastnz, n;
    int resQBits;
    int rest_bits, rest_bits_overflow;
    int nt_half;
    int c[2], *ctx;
    int p1, p2;
    int ii[2], idx1, idx2, idx;
    int numPeakIndicesOrig, numHoleIndices;
    int nbbits_m2;


    set_i(x, 0, nt);

    /* Rate flag */
    if (nbbits>400)
    {
        rateFlag = 2 << NBITS_CONTEXT;
    }
    else
    {
        rateFlag = 0;
    }

    /*Decode number of ntuples*/
    start_bit_pos = st->next_bit_pos;
    lsbs_bit_pos = start_bit_pos + nbbits - 1;

    n = 0;
    k = 1;
    nt_half = nt>>1;

    while (k<nt_half)
    {
        ++n;
        k = k<<1;
    }
    n = get_next_indice(st, n) + 1;

    /* Init */
    c[0] = c[1] = 0;

    t = 0;

    lastnz = n << 1;

    if (lastnz > nt)
    {
        st->BER_detect = 1;
        return 0;
    }

    if (hm_cfg)
    {
        /* mapped domain */
        numPeakIndicesOrig = hm_cfg->numPeakIndices;
        hm_cfg->numPeakIndices = min(hm_cfg->numPeakIndices, lastnz );
        numHoleIndices = lastnz - hm_cfg->numPeakIndices;

        /* Mark hole indices beyond lastnz as pruned */
        for (k=numHoleIndices; k<hm_cfg->numHoleIndices; ++k)
        {
            hm_cfg->holeIndices[k] = hm_cfg->holeIndices[k] + nt;
        }

        ii[0] = numPeakIndicesOrig;
        ii[1] = 0;

        p1 = p2 = 0;    /* to avoid compilation warnings */
    }
    else
    {
        /* unmapped domain */
        ii[0] = 0;
        p1 = p2 = 0;
    }

    /* Start Decoding */

    ari_start_decoding_14bits(st, &as);
    overflow_bit_pos = st->next_bit_pos;

    nbbits_m2 = nbbits + cbitsnew - 2;
    rest_bits_overflow = rest_bits = -nbbits_m2;

    /* Main Loop through the 2-tuples */
    for (k=0; k<lastnz; k+=2)
    {
        if (hm_cfg)
        {
            a1_i = get_next_coeff_mapped(ii, &p1, &idx1, hm_cfg);
            b1_i = get_next_coeff_mapped(ii, &p2, &idx2, hm_cfg);
        }
        else
        {
            a1_i = get_next_coeff_unmapped(ii, &idx1);
            b1_i = get_next_coeff_unmapped(ii, &idx2);
        }

        idx = min(idx1, idx2);

        /* Get context */
        ctx = &c[p1 | p2];

        t = *ctx + rateFlag;
        t -= ((nt_half - idx) >> (sizeof(int)*8-1)) << NBITS_CONTEXT;

        esc_nb = 0;
        r = 0;
        if (t < 0 || lsbs_bit_pos <= 0)
        {
            st->BER_detect = 1;
            return 0;
        }

        a = b = 0;
        /* MSBs decoding */
        for (lev = 0; lev < 15 && lsbs_bit_pos > 0; ++lev)
        {
            esc_nb = min(lev, 3);
            pki = ari_lookup_s17_LC[t + (esc_nb << (NBITS_CONTEXT + NBITS_RATEQ))];
            ari_decode_14bits_s17_ext(st,&r,&as,ari_pk_s17_LC_ext[pki]);

            if (r<VAL_ESC)
            {
                break;
            }

            /* LSBs decoding */
            a += get_indice_1(st, lsbs_bit_pos--) << lev;
            b += get_indice_1(st, lsbs_bit_pos--) << lev;
        }

        if ( (lsbs_bit_pos < -1 && r >= VAL_ESC) || (lev > 14) )
        {
            x[a1_i] = 0;
            x[b1_i] = 0;
            st->BER_detect = 1;
            return 0;
        }

        /* MSBs contributions */
        b1 = r>>2;
        a1 = r&0x3;
        a += a1 << lev;
        b += b1 << lev;

        /* lsbs bits sign bits */
        rest_bits += 2*lev;

        rest_bits += min(a, 1);
        rest_bits += min(b, 1);

        /* Dectect overflow */

        if (st->next_bit_pos-start_bit_pos+rest_bits>0)
        {
            /* Roll back bit-stream position to overflow_bit_pos */
            get_next_indice_tmp(st, overflow_bit_pos - st->next_bit_pos);
            rest_bits = rest_bits_overflow;
            x[a1_i] = 0;
            x[b1_i] = 0;
            break;
        }

        overflow_bit_pos = st->next_bit_pos;
        rest_bits_overflow = rest_bits;

        /* Store decoded data */
        x[a1_i] = a;
        x[b1_i] = b;

        /* Update context for next 2-tuple */
        if (p1 == p2)
        {
            /* peak-peak or hole-hole context */
            lev = esc_nb - 1;

            if (lev <= 0)
            {
                t = 1 + (a1 + b1)*(lev+2);
            }
            else
            {
                t = 13 + lev;
            }

            *ctx = (*ctx & 0xf) * 16 + t;
        }
        else
        {
            /* mixed context */

            if (idx1 & 1)
            {
                /* update first context */
                c[p1] = update_mixed_context(c[p1], a);
            }

            if (idx2 & 1)
            {
                /* update second context */
                c[p2] = update_mixed_context(c[p2], b);
            }
        }
    }

    /* Total number of decoded AC bits */
    get_next_indice_tmp(st, -(cbitsnew - 2));

    /* detect overflow */

    if(k!=lastnz)
    {
        rest_bits += nbbits_m2;
        /* Set bit-stream position to (start_bit_pos+nbbits-rest_bits) */
        get_next_indice_tmp(st, (start_bit_pos+nbbits-rest_bits)-st->next_bit_pos);
    }

    /* Decode signs */
    if (hm_cfg)
    {
        n = nt;
    }
    else
    {
        n = lastnz;
    }

    for (k = 0; k < n; k++)
    {
        if (x[k] > 0)
        {
            x[k] *= 1-2*get_next_indice_1(st);
        }
    }

    /*Decode Residual Q*/
    resQBits = min(resQMaxBits, lsbs_bit_pos+1 - st->next_bit_pos);

    for (k=0; k<resQBits; ++k)
    {
        x[nt+k] = get_indice_1(st, lsbs_bit_pos - k);
    }

    /* Set bit-stream pointer to end of buffer */
    get_next_indice_tmp(st, (start_bit_pos+nbbits)-st->next_bit_pos);


    return resQBits;
}
