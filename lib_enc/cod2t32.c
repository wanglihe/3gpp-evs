/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "prot.h"


/*---------------------------------------------------------------------*
 * Local constants
 *---------------------------------------------------------------------*/

#define STEP      2
#define MSIZE     1024


/*----------------------------------------------------------------------------------
 * Function  acelp_2t32()
 *
 * 12 bits algebraic codebook.
 * 2 tracks x 32 positions per track = 64 samples.
 *
 * 12 bits --> 2 pulses in a frame of 64 samples.
 *
 * All pulses can have two (2) possible amplitudes: +1 or -1.
 * Each pulse can have 32 possible positions.
 *----------------------------------------------------------------------------------*/

void acelp_2t32(
    Encoder_State *st,              /* i/o: encoder state structure                       */
    const float dn[],             /* i  : corr. between target and h[].                 */
    const float h[],              /* i  : impulse response of weighted synthesis filter */
    float code[],           /* o  : algebraic (fixed) codebook excitation         */
    float y[]               /* o  : filtered fixed codebook excitation            */
)
{
    short i, j, k, i0, i1, ix, iy, pos, pos2, index;
    float psk, ps1, ps2, alpk, alp1, alp2, sq;
    float pol[L_SUBFR], dn_p[L_SUBFR], r0;
    short ii,jj;
    float s, cor, sign0, sign1;
    float *p0, *p1, *p2;
    const float *ptr_h1, *ptr_h2, *ptr_hf;
    float rrixix[NB_TRACK_FCB_2T][NB_POS_FCB_2T];
    float rrixiy[MSIZE];


    /*----------------------------------------------------------------*
     * Compute rrixix[][] needed for the codebook search.
     *----------------------------------------------------------------*/

    /* Init pointers to last position of rrixix[] */
    p0 = &rrixix[0][NB_POS_FCB_2T - 1];
    p1 = &rrixix[1][NB_POS_FCB_2T - 1];

    ptr_h1 = h;
    cor = 0.0f;
    for (i = 0; i < NB_POS_FCB_2T; i++)
    {
        cor += *ptr_h1 **ptr_h1;
        ptr_h1++;
        *p1-- = cor;
        cor += *ptr_h1 **ptr_h1;
        ptr_h1++;
        *p0-- = cor;
    }

    p0 = rrixix[0];
    p1 = rrixix[1];

    for (i = 0; i < NB_POS_FCB_2T; i++)
    {
        *p0 = 0.5f * (*p0);
        p0++;
        *p1 = 0.5f * (*p1);
        p1++;
    }

    /*------------------------------------------------------------*
     * Compute rrixiy[][] needed for the codebook search.
     *------------------------------------------------------------*/

    pos = MSIZE - 1;
    pos2 = MSIZE - 2;
    ptr_hf = h + 1;

    for (k = 0; k < NB_POS_FCB_2T; k++)
    {
        /* Init pointers to last position of diagonals */
        p1 = &rrixiy[pos];
        p0 = &rrixiy[pos2];

        cor = 0.0f;
        ptr_h1 = h;
        ptr_h2 = ptr_hf;

        for (i = k+1; i < NB_POS_FCB_2T; i++)
        {
            cor += *ptr_h1++ **ptr_h2++;
            *p1 = cor;

            cor += *ptr_h1++ **ptr_h2++;
            *p0 = cor;

            p1 -= (NB_POS_FCB_2T + 1);
            p0 -= (NB_POS_FCB_2T + 1);
        }

        cor += *ptr_h1++ **ptr_h2;
        *p1 = cor;

        pos -= NB_POS_FCB_2T;
        pos2--;
        ptr_hf += STEP;
    }

    /*----------------------------------------------------------------*
     * computing reference vector and pre-selection of polarities
     *----------------------------------------------------------------*/

    for(i=0; i<L_SUBFR; i++)
    {
        /* FIR high-pass filtering */
        if (i==0)
        {
            r0=dn[i]-dn[i+1]*0.35f;
        }
        else if (i==L_SUBFR-1)
        {
            r0=-dn[i-1]*0.35f+dn[i];
        }
        else
        {
            r0=-dn[i-1]*0.35f+dn[i]-dn[i+1]*0.35f;
        }

        /* pre-selection of polarities */
        if (r0>=0.0f)
        {
            pol[i]= 1.0f;
        }
        else
        {
            pol[i]=-1.0f;
        }

        /* including polarities into dn[] */
        dn_p[i]=dn[i]*pol[i];
    }

    /*----------------------------------------------------------------*
     * compute denominator ( multiplied by polarity )
     *----------------------------------------------------------------*/

    k=0;
    ii=0;
    for(i=0; i<NB_POS_FCB_2T; i++)
    {
        jj=1;
        for(j=0; j<NB_POS_FCB_2T; j++)
        {
            rrixiy[k+j]*=pol[ii]*pol[jj];
            jj+=2;
        }
        ii+=2;
        k+=NB_POS_FCB_2T;
    }

    /*----------------------------------------------------------------*
     * search 2 pulses
     * All combinaisons are tested:
     * 32 pos x 32 pos x 2 signs = 2048 tests
     *----------------------------------------------------------------*/

    p0 = rrixix[0];
    p1 = rrixix[1];
    p2 = rrixiy;

    psk = -1;
    alpk = 1;
    ix = 0;
    iy = 1;
    for (i0 = 0; i0 < L_SUBFR; i0 += STEP)
    {
        ps1 = dn_p[i0];
        alp1 = *p0++;
        pos = -1;
        for (i1 = 1; i1 < L_SUBFR; i1 += STEP)
        {
            ps2 = ps1 + dn_p[i1];
            alp2 = alp1 + *p1++ + *p2++;
            sq = ps2 * ps2;
            s = alpk * sq - psk * alp2;
            if (s > 0)
            {
                psk = sq;
                alpk = alp2;
                pos = i1;
            }
        }
        p1 -= NB_POS_FCB_2T;
        if (pos >= 0)
        {
            ix = i0;
            iy = pos;
        }
    }

    i0 = ix/STEP;
    i1 = iy/STEP;
    sign0 = pol[ix];
    sign1 = pol[iy];


    /*-------------------------------------------------------------------*
     * Build the codeword, the filtered codeword and index of codevector.
     *-------------------------------------------------------------------*/

    set_f( code, 0.0f, L_SUBFR );

    code[ix] = sign0;
    code[iy] = sign1;
    index = (i0<<6) + i1;

    if(sign0 < 0.0f)
    {
        index += 0x800;
    }

    if(sign1 < 0.0f)
    {
        index += 0x20;
    }

    set_f( y, 0.0f, L_SUBFR );
    for(i=ix; i<L_SUBFR; i++)
    {
        y[i] = (sign0*h[i-ix]);
    }

    for(i=iy; i<L_SUBFR; i++)
    {
        y[i] += (sign1*h[i-iy]);
    }

    /* write index to array of indices */
    push_indice( st, IND_ALG_CDBK_2T32, index, 12 );

    return;
}

/*----------------------------------------------------------------------------------
 * acelp_1t64()
 *
 * 7 bits algebraic codebook.
 * 1 track x 64 positions per track = 64 samples.
 *
 * The pulse can have 64 possible positions and two (2) possible amplitudes: +1 or -1.
 *----------------------------------------------------------------------------------*/

void acelp_1t64(
    Encoder_State *st,              /* i/o: encoder state structure                       */
    const float dn[],             /* i  : corr. between target and h[].                 */
    const float h[],              /* i  : impulse response of weighted synthesis filter */
    float code[],           /* o  : algebraic (fixed) codebook excitation         */
    float y[]               /* o  : filtered fixed codebook excitation            */
)
{
    short i, pos, sgn, index;
    float tmp;

    /*-------------------------------------------------------------------*
     * Find position and sign of maximum impulse.
     *-------------------------------------------------------------------*/

    pos = emaximum( dn, L_SUBFR, &tmp );
    sgn = (short)sign( dn[pos] );

    /*-------------------------------------------------------------------*
     * Build the codeword, the filtered codeword and index of codevector.
     *-------------------------------------------------------------------*/

    set_f( code, 0.0f, L_SUBFR );
    code[pos] = sgn;

    set_f( y, 0.0f, L_SUBFR );

    for( i=pos; i<L_SUBFR; i++ )
    {
        y[i] = (sgn*h[i-pos]);
    }

    index = pos;

    if( sgn > 0 )
    {
        index += L_SUBFR;
    }

    push_indice( st, IND_ALG_CDBK_1T64, index, 7 );

    return;
}
