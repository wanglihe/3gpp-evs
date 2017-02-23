/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"

/*---------------------------------------------------------------------*
 * Local functions
 *---------------------------------------------------------------------*/

static void gaus_dec2v( Decoder_State *st, float *code, const short lg, const short nb_bits );
static void dec_2pos( short index, short *ind1, short *ind2, float *sign1, float *sign2, const short n );

/*---------------------------------------------------------------------*
 * gaus_dec()
 *
 * no adaptive excitation constructed
 * - decode the codebook indices,
 * - find the excitation
 *---------------------------------------------------------------------*/

void gaus_dec(
    Decoder_State *st,               /* i/o: decoder state structure   */
    const long core_brate,        /* i  : core bitrate                                */
    const short i_subfr,           /* i:   subframe index                              */
    float *code,             /* o:   unvoiced excitation                         */
    float *norm_gain_code,   /* o:   gain of the normalized gaussian excitation  */
    float *lp_gainp,         /* i/o: LP-filtered pitch gain (FEC)                */
    float *lp_gainc,         /* i/o: LP-filtered code gain (FEC)                 */
    float *gain_inov,        /* o:   unscaled innovation gain                    */
    float *tilt_code,        /* o:   synthesis excitation spectrum tilt          */
    float *voice_fac,        /* o:   estimated voicing factor                    */
    float *gain_pit,         /* o:   pitch gain                                  */
    float *pt_pitch,         /* o:   floating pitch buffer                       */
    float *exc,              /* o:   excitation signal frame                     */
    float *gain_code,        /* o:   gain of the gaussian excitation             */
    float *exc2              /* o  : Scaled excitation signal frame              */
)
{
    short i, idx, nb_bits;

    /*----------------------------------------------------------------*
     * Decode Gaussian excitation
     *----------------------------------------------------------------*/

    nb_bits = FCB_bits_tbl[BIT_ALLOC_IDX(core_brate, UNVOICED, i_subfr, TC_SUBFR2IDX(-1))];

    gaus_dec2v( st, code, L_SUBFR, nb_bits>>1 );

    /*-----------------------------------------------------------------*
     * Decode gain of Gaussian excitation and normalized Gaussian excitation
     *-----------------------------------------------------------------*/

    *gain_inov = 1.0f / (float)sqrt( ( dotp( code, code, L_SUBFR ) + 0.01f ) / L_SUBFR );

    nb_bits = gain_bits_tbl[BIT_ALLOC_IDX(core_brate, UNVOICED, i_subfr, TC_SUBFR2IDX(-1))];
    idx = (short)get_next_indice( st, nb_bits );

    /* safety check in case of bit errors */
    if( idx > 78 )
    {
        idx = 78;
        st->BER_detect = 1;
    }

    *gain_code = gain_dec_gaus( idx, nb_bits, -30.0f, 190.0f, *gain_inov, norm_gain_code );

    /* update LP filtered gains for the case of frame erasures */
    lp_gain_updt( i_subfr, 0.0f, *norm_gain_code, lp_gainp, lp_gainc, L_FRAME );   /* supposes that gain_dec_gaus() is used for ACELP@12k8 only */

    /*-----------------------------------------------------------------*
     * Updates
     *-----------------------------------------------------------------*/

    *tilt_code = 0.0f;
    *voice_fac = -1.0f;           /* only unvoiced              */
    *gain_pit = 0.0f;             /* needed for BASS postfitler */
    *pt_pitch = (float)L_SUBFR;   /* floating pitch buffer      */

    /*-----------------------------------------------------------------*
     * Construct scaled excitation
     *-----------------------------------------------------------------*/

    set_f(&exc2[i_subfr], 0, L_SUBFR);
    for (i = 0; i < L_SUBFR;  i++)
    {
        exc[i+i_subfr] = *gain_code * code[i];
    }

    return;
}

/*-----------------------------------------------------*
 * gaus_dec2v()
 *
 * decoder of Gaussian Codebook for unvoiced
 * consisting of addition of 2 Gaussian vectors
 *
 * One Gaussian vector of 190 values
 *-----------------------------------------------------*/

static void gaus_dec2v(
    Decoder_State *st,        /* i/o: decoder state structure   */
    float *code,      /* o:   decoded gaussian codevector */
    const short lg,         /* i:   codevector length           */
    const short nb_bits     /* i:   nb ob bits per track (max 6)*/
)
{
    short i, ind1, ind2, idx, index_delta;
    short nvec, step;
    float sign1, sign2;
    float *pt1, *pt2;
    float gaus_dico2[190];
    float delta;

    nvec = 1 << nb_bits;
    step = 0x80 >> nb_bits;

    idx = (short)get_next_indice( st, 2*nb_bits+1 );
    index_delta = (short)get_next_indice( st, 3 );

    dec_2pos( idx, &ind1, &ind2, &sign1, &sign2, nvec );

    delta = STEP_DELTA * (float)(index_delta);
    if( delta > 0.0f )
    {
        gaus_dico2[0] = gaus_dico[0];
        for (i=1; i<190; i++)
        {
            gaus_dico2[i] = (gaus_dico[i] - delta*gaus_dico[i-1])/(1 + delta*delta);
        }
    }
    else
    {
        for (i=0; i<190; i++)
        {
            gaus_dico2[i] = gaus_dico[i];
        }
    }

    pt1 = &gaus_dico2[ind1 * step];
    pt2 = &gaus_dico2[ind2 * step];

    for(i=0; i<lg; i++)
    {
        code[i] = pt1[i]*sign1 + pt2[i]*sign2;
    }

    return;
}

/*-----------------------------------------------------*
 * dec_2pos()
 *
 * Decode the codevectors positions and signs
 *-----------------------------------------------------*/

static void dec_2pos(
    short index,   /* i:   quantization index    */
    short *ind1,   /* o:   1st vector index      */
    short *ind2,   /* o:   2nd vector index      */
    float *sign1,  /* o:   1st vector sign       */
    float *sign2,  /* o:   2nd vector sign       */
    const short n        /* i:   nb. of vectors in cb. */
)
{
    short i;

    i = index & 1;
    if(i == 0)
    {
        *sign1 = 1.0f;
    }
    else
    {
        *sign1 = -1.0f;
    }

    index = index>>1;

    *ind1 = index/n;
    *ind2 = index -(*ind1*n);
    if( *ind1 > *ind2)
    {
        *sign2 = -*sign1;
    }
    else
    {
        *sign2 = *sign1;
    }

    return;
}





/*-----------------------------------------------------*
 * gaus_L2_dec :
 *
 * decoder of Gaussian Codebook for unvoiced as Layer 2
 *
 * One Gaussian vector
 *-----------------------------------------------------*/

void gaus_L2_dec(
    float *code,            /* o:   decoded gaussian codevector */
    float tilt_code,
    const float *Aq,
    float formant_enh_num,
    short *seed_acelp       /*i/o: random seed */
)
{
    short i;

    /*Generate white gaussian noise using central limit theorem method (N only 4 as E_util_random is not purely uniform)*/
    for( i=0; i<L_SUBFR; i++ )
    {
        code[i]=(float)(own_random(seed_acelp))/(1<<15);
        code[i]+=(float)(own_random(seed_acelp))/(1<<15);
        code[i]+=(float)(own_random(seed_acelp))/(1<<15);
    }

    /*Shape the gaussian excitation*/
    cb_shape( 1, 0, 0, 1, 0, formant_enh_num, FORMANT_SHARPENING_G2, Aq, code, tilt_code, 0 );

    return;
}

