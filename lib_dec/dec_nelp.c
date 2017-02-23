/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "prot.h"

/*-------------------------------------------------------------------*
* decod_nelp()
*
* Decode unvoiced NELP
*-------------------------------------------------------------------*/

void decod_nelp(
    Decoder_State *st,              /* i/o: decoder static memory                   */
    const short coder_type,       /* i : coding type                              */
    float *tmp_noise,       /* o : long term temporary noise energy         */
    float *pitch_buf,       /* o : floating pitch values for each subframe  */
    float *exc,             /* o : adapt. excitation exc                    */
    float *exc2,            /* o : adapt. excitation/total exc              */
    float *voice_factors,   /* o  : voicing factors                         */
    float *bwe_exc,         /* o  : excitation for SWB TBE                  */
    const short bfi,              /* i : bad frame indicator                      */
    float *gain_buf
)
{
    short i;
    float exc_nelp[L_FRAME];

    *tmp_noise = 0;

    nelp_decoder( st,exc_nelp,exc,bfi,coder_type, gain_buf );

    mvr2r( exc_nelp, exc, L_FRAME );
    mvr2r( exc_nelp, exc2, L_FRAME );

    st->tilt_code = 0.0f; /* purely unvoiced */
    set_f( st->tilt_code_dec, 0, NB_SUBFR16k );
    st->prev_tilt_code_dec = st->tilt_code;

    st->dispMem[0] = 0;
    st->prev_gain_pit_dec = 0.0;
    st->dispMem[2] = st->prev_gain_pit_dec;

    for(i=3; i<7; i++)
    {
        st->dispMem[i] = st->dispMem[i-1];
    }

    set_f(pitch_buf, L_SUBFR, NB_SUBFR);
    interp_code_5over2( exc2, bwe_exc, L_FRAME );

    set_f( voice_factors, 0.0f, NB_SUBFR16k );

    return;
}
