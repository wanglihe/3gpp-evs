/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"

/*-------------------------------------------------------------------
 * decod_ppp()
 *
 * PPP decoder
 *-------------------------------------------------------------------*/

void decod_ppp(
    Decoder_State *st,           /* i/o: state structure */
    const float Aq[],          /* i  : 12k8 Lp coefficient */
    float *pitch_buf,    /* i/o: floating pitch values for each subframe */
    float *exc,          /* i/o: current non-enhanced excitation */
    float *exc2,         /* i/o: current enhanced excitation */
    float *voice_factors, /* o  : voicing factors */
    float *bwe_exc       /* o  : excitation for SWB TBE */
    , float *gain_buf
    , short bfi
)
{
    short k;
    float p_Aq_old[M+1], excQ_ppp[L_FRAME], p_Aq_curr[M], LPC_de_old[M+1];
    float LPC_de_curr[M+1], pitch[NB_SUBFR];

    /* call voiced decoder at this point */
    for( k=0; k<M; k++)
    {
        p_Aq_curr[k] = Aq[k+(3*(M+1))+1];
    }

    lsp2a_stab( st->lsp_old, p_Aq_old, M );

    deemph_lpc( p_Aq_curr, p_Aq_old, LPC_de_curr, LPC_de_old
                ,0
              );

    /* last frame-end lpc and curr frame-end lpc */
    ppp_voiced_decoder( st, excQ_ppp, LPC_de_curr, exc, pitch
                        ,bfi
                      );

    st->tilt_code = st->prev_tilt_code_dec;

    mvr2r( excQ_ppp, exc, L_FRAME );
    mvr2r( exc, exc2, L_FRAME );

    st->dispMem[0] = 2;
    st->dispMem[2] = st->prev_gain_pit_dec;

    for(k=3; k<7; k++)
    {
        st->dispMem[k] = st->dispMem[k-1];
    }

    mvr2r( pitch, pitch_buf, NB_SUBFR );

    interp_code_5over2( exc2, bwe_exc, L_FRAME );
    set_f( voice_factors, 0.0f, NB_SUBFR16k );
    set_f(gain_buf,0,NB_SUBFR16k);

    return;
}
