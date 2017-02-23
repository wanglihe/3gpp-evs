/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"

/*----------------------------------------------------------------------------------*
 * CNG_reset_dec()
 *
 * Reset decoder static variables in case of CNG frame
 *----------------------------------------------------------------------------------*/

void CNG_reset_dec(
    Decoder_State *st,          /* i/o: decoder state structure            */
    float *pitch_buf,     /* o  : floating pitch for each subframe   */
    float *voice_factors  /* o  : voicing factors                    */
)
{
    mvr2r( UVWB_Ave, st->mem_AR, M );
    set_f(st->mem_MA, 0, M );
    set_f( st->dispMem, 0, 8 );
    st->tilt_code = 0.0f;
    st->gc_threshold = 0.0f;

    /* last good received frame for FEC in ACELP */
    st->clas_dec = UNVOICED_CLAS;
    st->last_good = UNVOICED_CLAS;

    /* LP-filtered pitch gain set to 0 */
    st->lp_gainp = 0.0f;

    /* convert CNG energy into CNG gain for ACELP FEC */
    st->lp_gainc = (float)sqrt( st->lp_ener );

    /* reset the pitch buffer in case of FRAME_NO_DATA or SID frames */
    if(st->L_frame == L_FRAME )
    {
        set_f( pitch_buf, (float)L_SUBFR, NB_SUBFR );
    }
    else  /* st->L_frame == L_FRAME16k */
    {
        set_f( pitch_buf, (float)L_SUBFR16k, NB_SUBFR16k );
    }

    set_f( voice_factors, 1.0, NB_SUBFR16k );

    /* deactivate bass post-filter */
    st->bpf_off = 1;

    /* Reset active frame counter */
    st->act_cnt2 = 0;

    return;
}
