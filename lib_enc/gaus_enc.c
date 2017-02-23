/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"

/*-------------------------------------------------------------------*
 * Local constants
 *-------------------------------------------------------------------*/

#define NMAX          8  /* Control of the routine's complexity */
#define FAC_DELTA  16.0f

/*---------------------------------------------------------------------*
 * Prototypes
 *---------------------------------------------------------------------*/

static short cod_2pos( const short ind1,const short ind2,const float sign1,
                       const float sign2,const short n );

static void gauss2v( Encoder_State *st, const float h[], const float xn[], const float dn[], float code[],
                     float y1[], float *gain, const short lg, const short nb_bits );

/*-------------------------------------------------------------------*
 * Gaus_encode
 *
 * Encoder UnVoiced excitation coding using Gaussian codebooks
 * - ACELP quantized Gaussian excitation
 * - gain quantization
 * - Total excitation for UnVoiced coders
 * - Updates
 *-------------------------------------------------------------------*/

float gaus_encode(
    Encoder_State *st,                  /* i/o: encoder state structure      */
    const short i_subfr,              /* i  : subframe index                             */
    const float *h1,                  /* i  : weighted filter input response             */
    const float *xn,                  /* i  : target vector                              */
    float *exc,                 /* o  : pointer to excitation signal frame         */
    float *mem_w0,              /* o  : weighting filter denominator memory        */
    float *gp_clip_mem,         /* o  : memory of gain of pitch clipping algorithm */
    float *tilt_code,           /* o  : synthesis excitation spectrum tilt         */
    float *code,                /* o  : algebraic excitation                       */
    float *gain_code,           /* o  : Code gain.                                 */
    float *y2,                  /* o  : zero-memory filtered adaptive excitation   */
    float *gain_inov,           /* o  : innovation gain                            */
    float *voice_fac,           /* o  : voicing factor                             */
    float *gain_pit,            /* o  : adaptive excitation gain                   */
    float *norm_gain_code,      /* o  : normalized innovative cb. gain             */
    const long  core_brate            /* i  : core bitrate                               */
)
{
    short i = 0, nb_bits, idx;
    float dn[L_SUBFR];              /* Correlation between xn and h1     */

    /*----------------------------------------------------------------*
     *  Encode gaussian excitation
     *----------------------------------------------------------------*/

    corr_xh( xn, dn, h1, L_SUBFR ); /* Correlation between target xn[] and impulse response h1[] */

    nb_bits = FCB_bits_tbl[BIT_ALLOC_IDX(core_brate, UNVOICED, i_subfr, TC_SUBFR2IDX(-1))];

    gauss2v( st, h1, xn, dn, code, y2, gain_code, L_SUBFR, nb_bits>>1 );

    /*----------------------------------------------------------------*
     *  Encode gaussian gain
     *----------------------------------------------------------------*/

    /* codeword energy computation */
    *gain_inov = 1.0f / (float)sqrt( ( dotp( code, code, L_SUBFR ) + 0.01f ) / L_SUBFR );

    nb_bits = gain_bits_tbl[BIT_ALLOC_IDX(core_brate, UNVOICED, i_subfr, TC_SUBFR2IDX(-1))];

    idx = gain_enc_gaus( gain_code, nb_bits, -30.0f, 190.0f );
    push_indice( st, IND_GAIN, idx, nb_bits );

    /*-----------------------------------------------------------------*
     * Total excitation for Unvoiced coders
     *-----------------------------------------------------------------*/

    for (i = 0; i < L_SUBFR;  i++)
    {
        exc[i+i_subfr] = *gain_code * code[i];
    }

    /*-----------------------------------------------------------------*
     * Updates: last value of new target is stored in mem_w0
     *-----------------------------------------------------------------*/

    *mem_w0 = xn[L_SUBFR-1] - *gain_code * y2[L_SUBFR-1];

    init_gp_clip(gp_clip_mem);          /* reset pitch clipping parameters  */
    *gain_pit = 0.0f;
    *tilt_code = 0.0f;                  /* purely unvoiced  */
    *voice_fac = -1.0f;                 /* purely unvoiced  */

    *norm_gain_code = *gain_code / *gain_inov;

    return L_SUBFR;
}

/*-------------------------------------------------------------------*
 * gauss2v()
 *
 * encoder of Gaussian Codebook for unvoiced
 * consisting of addition of 2 Gaussian vectors
 *
 * One Gaussian vector of 192 values vectors delayed by 2
 *-------------------------------------------------------------------*/

static void gauss2v(
    Encoder_State *st,        /* i/o: encoder state structure      */
    const float h[],        /* i  : weighted LP filter impulse response     */
    const float xn[],       /* i  : target signal                           */
    const float dn[],       /* i  : backward filtered target                */
    float code[],     /* o  : gaussian excitation                     */
    float y11[],      /* o  : zero-memory filtered gauss. excitation  */
    float *gain,      /* o  : excitation gain                         */
    const short lg,         /* i  : subframe size                           */
    const short nb_bits     /* i  : nb ob bits per track (max 6)            */
)
{
    short i, j, ind1, ind2, idx;
    short nvec, step;
    float cor, cora, cor2, cor2w, eneri, enerw;
    float *pt1, *pt2;
    float max[NMAX+1], *pos[NMAX+1], sign[NMAX+1];
    float ener[NMAX+1], corr[NMAX+1], ener1;
    float dico2[L_SUBFR*NMAX];
    float c1, c0;
    float gxx, gcc;
    float gaus_dico2[190];
    float hg[190];
    float delta;
    short index_delta;


    /*-----------------------------------------------------------------*
     * Encode the tilt of gaussian excitation
     *-----------------------------------------------------------------*/

    c0=0.0f;                              /* Compute spectral tilt of target */
    c1=0.0f;
    for (i=1; i<L_SUBFR; i++)
    {
        c0 += xn[i]*xn[i];
        c1 += xn[i]*xn[i-1];
    }
    if(c0 < FLT_MIN)
    {
        gxx = 0.f;
    }
    else
    {
        gxx = c1/c0;
    }

    set_f(hg, 0, 190);                    /* Compute spectral tilt of filtered codebook */
    mvr2r(h, hg, L_SUBFR);
    conv(gaus_dico, hg, gaus_dico2, 190);

    c0=0.0f;
    c1=0.0f;
    for (i=1; i<190; i++)
    {
        c0 += gaus_dico2[i]*gaus_dico2[i];
        c1 += gaus_dico2[i]*gaus_dico2[i-1];
    }
    gcc = c1/c0;
    delta = (1-gcc*gxx) / (2*gcc+gxx);    /* Compute and quantize spectral tilt modification factor (3b) */

    index_delta = (short)(FAC_DELTA * delta);
    if (index_delta < 0)
    {
        index_delta = 0;
    }
    if (index_delta > 7)
    {
        index_delta = 7;
    }

    delta = STEP_DELTA * (float)index_delta;
    if( delta > 0.0f )                    /* Adapt spectral tilt of initial codebook */
    {
        gaus_dico2[0] = gaus_dico[0];
        for (i=1; i<190; i++)
        {
            gaus_dico2[i] = (gaus_dico[i] - delta*gaus_dico[i-1])/(1+delta*delta);
        }
    }
    else
    {
        for (i=0; i<190; i++)
        {
            gaus_dico2[i] = gaus_dico[i];
        }
    }

    /*-----------------------------------------------------------------*
     * Codebook search initializations
     *-----------------------------------------------------------------*/

    ind1 = 0;
    ind2 = 0;

    nvec = 1 << nb_bits;
    step = 0x80 >> nb_bits;

    /*-----------------------------------------------------------------*
     * dot product between dn and gaussian codevectors,
     * keep NMAX best vectors
     *-----------------------------------------------------------------*/

    set_f(max, 0, NMAX+1);
    set_f(sign, 0, NMAX+1);

    for (i = 0; i < NMAX+1; i++)
    {
        pos[i] = (float *)gaus_dico2;
    }

    pt1 = (float *)gaus_dico2;

    for (i=0; i < nvec; i++, pt1+=step)
    {
        cor = dotp(pt1, dn, lg);
        cora = (float)fabs(cor);
        j = NMAX-1;
        do
        {
            if(cora >= max[j])
            {
                max[j+1]  = max[j];
                pos[j+1]  = pos[j];
                sign[j+1] = sign[j];
                max[j]    = cora;
                pos[j]    = pt1;
                sign[j]   = cor;
            }

            j--;
        }
        while (j >= 0);
    }

    /*-----------------------------------------------------------------*
     * filter selected vectors
     * put sign
     * compute energy
     *-----------------------------------------------------------------*/

    pt1 = dico2;
    for (i=0; i<NMAX; i++, pt1+=L_SUBFR)
    {
        conv(pos[i], h, pt1, lg);

        /* put sign and compute energy */
        if (sign[i] < 0.0f)
        {
            for(j=0; j<lg; j++)
            {
                pt1[j] = -pt1[j]; /*Store into dico2*/
            }
        }
        ener[i] = dotp(pt1, pt1, lg);
        corr[i] = dotp(pt1, xn,  lg);   /* must be equal to sign[i] */
    }

    /*-----------------------------------------------------------------*
     * try all combinations of NMAX best vectors
     *-----------------------------------------------------------------*/

    pt1 = dico2;

    /* Initial values for search algorithm */
    enerw = 1.0f;
    cor2w = -1.0f;

    for(i=0; i< NMAX; i++, pt1+=L_SUBFR)
    {
        pt2 = pt1;
        for(j=i; j<NMAX; j++, pt2+=L_SUBFR)
        {
            cor = corr[i] + corr[j];
            eneri = ener[i] + ener[j] + 2.0f * dotp(pt1, pt2, lg);
            cor2  = cor*cor;
            if(cor2 * enerw > cor2w * eneri)
            {
                cor2w = cor2;
                enerw = eneri;
                ind1 = i;
                ind2 = j;
            }
        }
    }

    /*-----------------------------------------------------------------*
     * Compute zero-memory filtered gauss. excitation y
     *-----------------------------------------------------------------*/

    pt1 = dico2 + ind1*L_SUBFR;
    pt2 = dico2 + ind2*L_SUBFR;
    for(i=0; i<lg; i++)
    {
        y11[i] = pt1[i] + pt2[i];
    }

    /*-----------------------------------------------------------------*
     * Signs of vectors
     *-----------------------------------------------------------------*/

    if(sign[ind1] >= 0.0f)
    {
        sign[ind1] = 1.0f;
    }
    else
    {
        sign[ind1] = -1.0f;
    }

    if(sign[ind2] >= 0.0f)
    {
        sign[ind2] = 1.0f;
    }
    else
    {
        sign[ind2] = -1.0f;
    }

    /*-----------------------------------------------------------------*
     * Compute code
     *-----------------------------------------------------------------*/

    pt1 = pos[ind1];
    pt2 = pos[ind2];
    for(i=0; i<lg; i++)
    {
        code[i] = pt1[i]*sign[ind1] + pt2[i]*sign[ind2];
    }
    cor = corr[ind1] + corr[ind2];

    /*-----------------------------------------------------------------*
     * Compute index
     *-----------------------------------------------------------------*/

    i = (short)((pos[ind1] - gaus_dico2) / step);  /* Division by step can be replaced by shift */
    j = (short)((pos[ind2] - gaus_dico2) / step);  /* Division by step can be replaced by shift */

    idx = cod_2pos(i, j, sign[ind1], sign[ind2], nvec);

    push_indice( st, IND_GAUS_CDBK_INDEX, idx, 2*nb_bits+1 );
    push_indice( st, IND_TILT_FACTOR, index_delta, 3 );

    /*-----------------------------------------------------------------*
     * Find quantized gain
     *-----------------------------------------------------------------*/

    *gain = cor / enerw;
    ener1 = dotp( xn, xn, lg );
    ener1 = (float)sqrt( ener1 / enerw );  /* Minimize ebergy error */
    *gain = *gain * 0.6f + ener1 * 0.4f;

    return;
}

/*---------------------------------------------------------------------*
 * Put selected codevector positions and signs into quantization index
 *---------------------------------------------------------------------*/

static short cod_2pos(  /* o  : codebook quantization index   */
    const short ind1,   /* i  : index of 1st gaussian vector  */
    const short ind2,   /* i  : index of 2nd gaussian vector  */
    const float sign1,  /* i  : sign of 1st gaussian vector   */
    const float sign2,  /* i  : sign of 2nd gaussian vector   */
    const short n       /* i  : nb. of codebook vectors       */
)
{
    short i1, i2, index, s1, s2;


    s1=1;
    if(sign1 > 0.0f)
    {
        s1=0;
    }

    s2=1;
    if(sign2 > 0.0f)
    {
        s2=0;
    }

    if(s1 == s2)
    {
        if (ind1<=ind2)
        {
            i1 = ind1;
            i2 = ind2;
        }
        else
        {
            i1 = ind2;
            i2 = ind1;
        }
    }
    else
    {
        if (ind1>ind2)
        {
            i1 = ind1;
            i2 = ind2;
        }
        else
        {
            i1 = ind2;
            i2 = ind1;
            s1 = s2;
        }
    }
    index = i1*n + i2;
    index = (index<<1) + s1;

    return index;
}
