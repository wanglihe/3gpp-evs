/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include <assert.h>
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"


/*-------------------------------------------------------------------*
* ACcontextMapping_encode2_no_mem_s17_LC()
*
* Arithmetic encoder
*-------------------------------------------------------------------*/

void ACcontextMapping_encode2_no_mem_s17_LC(
    Encoder_State *st,
    int *x,
    int nt,
    int lastnz,
    int nbbits,
    int resQMaxBits,
    CONTEXT_HM_CONFIG *hm_cfg)
{
    int ptr[BITBUFSIZE];
    Tastat as, as_overflow;
    int bp, bp_overflow;
    int a1, b1, a1_i, b1_i, k;
    int t, pki, lev1;
    int rateFlag;
    int value;
    int nbbits_ntuples, nbbits_lsbs, nbbits_signs,nbbits_signs_overflow, nbbits_lsbs_overflow, flag_overflow;
    int *lsbs_bits;
    int nt_half;
    int c[2], *ctx;
    int p1, p2;
    int ii[2], idx1, idx2, idx;
    int numPeakIndicesOrig, numHoleIndices;
    int signs[N_MAX];
    int nbbits_m2;

    a1 = 0;       /* to avoid compilation warnings */
    b1 = 0;       /* to avoid compilation warnings */

    /* Rate flag */
    if (nbbits > 400)
    {
        rateFlag = 2 << NBITS_CONTEXT;
    }
    else
    {
        rateFlag = 0;
    }

    /* Init */
    nt_half = nt >> 1;
    c[0] = c[1] = 0;

    /* Bits for encoding the number of encoded tuples */
    nbbits_ntuples = 0;
    k = 1;

    while (k<nt/2)
    {
        nbbits_ntuples++;
        k = k<<1;
    }

    t = 0;
    bp = nbbits_ntuples;
    nbbits_signs = 0;
    nbbits_lsbs = 0;
    nbbits_m2 = nbbits - 2;
    flag_overflow = 0;

    if (hm_cfg)
    {
        /* mapped domain */
        numPeakIndicesOrig = hm_cfg->numPeakIndices;
        hm_cfg->numPeakIndices = min(hm_cfg->numPeakIndices, lastnz);
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

        /* Find last non-zero tuple */
        /* ensure termination of while loop by dummy value */
        a1 = x[0];
        x[0] = 1;

        while (x[lastnz-1] == 0 && x[lastnz-2] == 0)
        {
            lastnz -= 2;
        }
        x[0] = a1;
    }

    lsbs_bits=ptr+nbbits-1;

    /*Start Encoding*/
    ari_start_encoding_14bits(&as);

    /*Main Loop through the 2-tuples*/
    b1_i = -1;

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
        t += (nt_half >= idx) ? 0 : (1 << NBITS_CONTEXT);

        /* Init current 2-tuple encoding */

        if (flag_overflow != 0)
        {
            x[a1_i] = 0;
            x[b1_i] = 0;
        }

        a1 = abs(x[a1_i]);
        b1 = abs(x[b1_i]);

        lev1 = -1;

        /*Copy states*/
        ari_copy_states(&as, &as_overflow);
        bp_overflow = bp;
        nbbits_signs_overflow = nbbits_signs;
        nbbits_lsbs_overflow = nbbits_lsbs;

        /*Signs encoding*/

        if (a1 > 0)
        {
            signs[nbbits_signs++] = ((unsigned int)x[a1_i] >> (sizeof(unsigned int)*8-1));
        }

        if (b1 > 0)
        {
            signs[nbbits_signs++] = ((unsigned int)x[b1_i] >> (sizeof(unsigned int)*8-1));
        }

        /* MSBs coding */
        while (max(a1, b1) >= A_THRES)
        {
            pki = ari_lookup_s17_LC[t + ((lev1+1) << (NBITS_CONTEXT + NBITS_RATEQ))];
            bp = ari_encode_14bits_ext(ptr,bp,&as,VAL_ESC,ari_pk_s17_LC_ext[pki]);

            *lsbs_bits-- = a1 & 1;
            *lsbs_bits-- = b1 & 1;

            /* LSBs bit counting */
            nbbits_lsbs += 2;

            a1 >>= 1;
            b1 >>= 1;

            lev1 = min(lev1+1, 2);
        }

        pki = ari_lookup_s17_LC[t + ((lev1+1) << (NBITS_CONTEXT + NBITS_RATEQ))];
        bp = ari_encode_14bits_ext(ptr,bp,&as,a1+A_THRES*b1,ari_pk_s17_LC_ext[pki]);


        /* Check bit budget */
        if (bp+as.vobf+nbbits_signs+nbbits_lsbs > nbbits_m2)
        {
            ari_copy_states(&as_overflow, &as);
            bp = bp_overflow;

            if (!flag_overflow)
            {
                nbbits_signs = nbbits_signs_overflow;
                nbbits_lsbs = nbbits_lsbs_overflow;

                if (hm_cfg)
                {
                    flag_overflow = 1;

                    /* Code from now only zeros */
                    x[a1_i] = 0;
                    x[b1_i] = 0;
                    lev1 = -1;

                    pki = ari_lookup_s17_LC[t];
                    bp = ari_encode_14bits_ext(ptr,bp,&as,0,ari_pk_s17_LC_ext[pki]);

                    if (bp+as.vobf+nbbits_signs+nbbits_lsbs>nbbits_m2)
                    {
                        ari_copy_states(&as_overflow, &as);
                        bp = bp_overflow;
                        break;
                    }
                }
                else
                {
                    break;
                }
            }
            else
            {
                break;
            }
        }

        /* Update context for next 2-tuple */
        if (p1 == p2)
        {
            /* peak-peak or hole-hole context */

            if (lev1 <= 0)
            {
                t = 1 + (a1 + b1)*(lev1+2);
            }
            else
            {
                t = 13 + lev1;
            }

            *ctx = (*ctx & 0xf) * 16 + t;
        }
        else
        {
            /* mixed context */

            if (idx1 & 1)
            {
                /* update first context */
                c[p1] = update_mixed_context(c[p1], abs(x[a1_i]));
            }

            if (idx2 & 1)
            {
                /* update second context */
                c[p2] = update_mixed_context(c[p2], abs(x[b1_i]));
            }
        }
    } /*end of the 2-tuples loop*/


    /* End arithmetic coder, overflow management */
    bp = ari_done_encoding_14bits(ptr,bp,&as);

    /*Overflow is detected*/

    if (k!=lastnz)
    {

        if (hm_cfg)
        {
            /*Fill with zero to be sure that decoder finish at the same position the MSB decoding*/

            for(; bp<nbbits-(nbbits_signs+nbbits_lsbs);)
            {
                ptr[bp++] = 0;
            }
        }
        else
        {
            lastnz = k;
        }
    }

    /* Push number of encoded tuples */
    value = (lastnz>>1)-1;
    push_next_indice(st, value, nbbits_ntuples);

    /* Push arithmetic coded bits */
    push_next_bits(st, &ptr[nbbits_ntuples], bp - nbbits_ntuples);

    /* Push sign bits */
    push_next_bits(st, signs, nbbits_signs);
    bp += nbbits_signs;

    /*write residual Quantization bits*/

    for(k=0; k<min(nbbits-bp-nbbits_lsbs,resQMaxBits); k++)
    {
        ptr[nbbits-1-nbbits_lsbs-k]=x[nt+k];
    }
    /* Write filler bits */
    for (; k<nbbits-bp-nbbits_lsbs; ++k)
    {
        ptr[nbbits-1-nbbits_lsbs-k]=0;
    }

    /* Check for debugging */
    assert( bp+k <= nbbits );

    /* Push the rest of the buffer */
    push_next_bits(st, &ptr[bp], nbbits - bp);

    /* return (bp+nbbits_lsbs);*/ /*return only for debug plot*/
    return;
}

/*-------------------------------------------------------------------*
* find_last_nz_pair()
*
*
*-------------------------------------------------------------------*/

static int find_last_nz_pair(
    const int x[],
    int length,
    const CONTEXT_HM_CONFIG *hm_cfg
)
{
    int last_nz, i;
    const int *tmp;

    last_nz = 2;

    if (hm_cfg)
    {
        /* mapped kernel */
        tmp = hm_cfg->indexBuffer;

        for (i=length; i >= 4; i-=2)
        {

            if (x[tmp[i-2]] != 0 || x[tmp[i-1]] != 0)
            {
                last_nz = i;
                break;
            }
        }
    }
    else
    {
        /* unmapped kernel */

        for (i=length; i >= 4; i-=2)
        {

            if (x[i-2] != 0 || x[i-1] != 0)
            {
                last_nz = i;
                break;
            }
        }
    }

    return last_nz;
}



/*-------------------------------------------------------------------*
* ACcontextMapping_encode2_estimate_no_mem_s17_LC()
*
*
*-------------------------------------------------------------------*/

int ACcontextMapping_encode2_estimate_no_mem_s17_LC(
    const int *x,
    int nt,
    int *lastnz_out,
    int *nEncoded,
    int target,
    int *stop,
    CONTEXT_HM_CONFIG *hm_cfg
)
{
    int a1, b1, a1_i, b1_i;
    int k, t, pki, lev1;
    int lastnz, lastnz2;
    int rateFlag;
    int nbits_old, nbits;
    int stop2;
    int proba;
    short nlz;
    const unsigned short *cum_freq;
    long symbol;
    const unsigned char *lookup;
    int nt_half;
    int c[2], *ctx;
    int p1, p2;
    int ii[2], idx1, idx2, idx;
    int numPeakIndicesOrig = 0, numHoleIndices = 0; /* only to avoid compiler warning */

    /* Rate flag */
    if (target > 400)
    {
        rateFlag = 2 << NBITS_CONTEXT;
    }
    else
    {
        rateFlag = 0;
    }

    /* 2 bits = arithmetic coder initialization interval = 1 bits for rounding last proba + 1 bit?*/
    nbits = 2;
    /*proba coded on 14bits -> proba=1*/
    proba = 16384;

    /* Init */
    nt_half = nt >> 1;
    stop2 = 0;
    c[0] = c[1] = 0;

    /* bits to encode lastnz */
    k = 1;

    while (k<nt/2)
    {
        nbits++;
        k = k<<1;
        /* check while condition */
    }
    nbits_old = nbits;

    nbits -= target;

    /* Find last non-zero tuple in the mapped domain signal */
    lastnz = find_last_nz_pair(x, nt, hm_cfg);

    /* At least one tuple is coded */
    lastnz2 = 2;

    if (hm_cfg)
    {
        /* mapped domain */
        numPeakIndicesOrig = hm_cfg->numPeakIndices;
        hm_cfg->numPeakIndices = min(hm_cfg->numPeakIndices, lastnz);
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
        t += (nt_half >= idx) ? 0 : (1 << NBITS_CONTEXT);

        /* Init current 2-tuple encoding */
        a1 = abs(x[a1_i]);
        b1 = abs(x[b1_i]);
        lev1 = -(1 << (NBITS_CONTEXT+NBITS_RATEQ));

        /* Signs Bits */
        nbits += min(a1, 1);
        nbits += min(b1, 1);

        /* pre-compute address of ari_pk_s17_LC_ext[0][Val_esc] to avoid doing it multiple times inside the loop */
        lookup = &ari_lookup_s17_LC[t] + (1 << (NBITS_CONTEXT+NBITS_RATEQ));

        /* check while condition */
        /* MSBs coding */
        while (max(a1, b1) >= A_THRES)
        {
            pki = lookup[lev1];
            cum_freq = ari_pk_s17_LC_ext[pki] + VAL_ESC;
            /*p1*p2=proba on 28 bits: p=0.5->power(2,27)*/
            proba *= *(cum_freq);
            /*Number of leading zero computed in one cycle=norm_l() in BASOP*/
            nlz=2;
            while(proba<134217728)    /*power(2,27)*/
            {
                nlz++;
                proba=proba<<1;
            }
            nbits+=nlz;
            /*addition added as shift not done in norm_l(): real shift = 14-nlz*/
            proba>>=14; /*proba is rounded down on 14 bits ->automatic over-estimation of bit consumption*/

            (a1) >>= 1;
            (b1) >>= 1;

            lev1 = min(lev1 + (1 << (NBITS_CONTEXT+NBITS_RATEQ)), 2 << (NBITS_CONTEXT+NBITS_RATEQ));
            /* check while condition */
        }
        pki = lookup[lev1];
        symbol = a1 + A_THRES*b1;
        cum_freq = ari_pk_s17_LC_ext[pki] + symbol;
        /*p1*p2=proba on 28 bits: p=0.5->power(2,27)*/
        proba *= (cum_freq[0] - cum_freq[1]);
        /*Number of leading zero computed in one cycle=norm_l() in BASOP*/
        nlz=0;
        while(proba<134217728)  /*power(2,27)*/
        {
            nlz++;
            proba=proba<<1;
        }

        nbits+=nlz;
        proba>>=14; /*proba is rounded down on 14 bits ->automatic over-estimation of bit consumption*/

        /* Should we truncate? */
        if ( nbits>0 )
        {
            stop2 = 1;

            if (*stop)
            {
                break;
            }
        }
        else
        {
            if (hm_cfg || k==0 || x[a1_i] || x[b1_i])
            {
                nbits_old = nbits+target;
                lastnz2 = b1_i+1;
            }
        }

        /* Update context for next 2-tuple */
        if (p1 == p2)   /* peak-peak or hole-hole context */
        {
            lev1 >>= NBITS_CONTEXT+NBITS_RATEQ;

            if (lev1 <= 0)
            {
                t = 1 + (a1 + b1)*(lev1+2);
            }
            else
            {
                t = 13 + lev1;
            }

            *ctx = (*ctx & 0xf) * 16 + t;
        }
        else
        {
            /* mixed context */

            if (idx1 & 1)
            {
                /* update first context */
                c[p1] = update_mixed_context(c[p1], abs(x[a1_i]));
            }

            if (idx2 & 1)
            {
                /* update second context */
                c[p2] = update_mixed_context(c[p2], abs(x[b1_i]));
            }
        }
    } /*end of the 2-tuples loop*/

    nbits += target;

    /* Output */

    if (*stop)
    {
        nbits = nbits_old;
    }

    if (stop2)
    {
        stop2 = nbits;
    }
    *nEncoded = lastnz2;
    *stop = stop2;
    *lastnz_out = lastnz;

    if (hm_cfg)
    {
        /* Restore hole indices beyond lastnz */
        for (k=numHoleIndices; k<hm_cfg->numHoleIndices; ++k)
        {
            hm_cfg->holeIndices[k] = hm_cfg->holeIndices[k] - nt;
        }
        hm_cfg->numPeakIndices = numPeakIndicesOrig;
    }

    return nbits_old;
}
