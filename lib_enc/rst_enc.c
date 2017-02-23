/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"


/*-------------------------------------------------------------------*
 * CNG_reset_enc()
 *
 * Reset encoder static variables after a CNG frame
 *-------------------------------------------------------------------*/

void CNG_reset_enc(
    Encoder_State *st,             /* i/o: encoder state structure          */
    LPD_state *mem,            /* i/o: encoder memories                 */
    float *pitch_buf,      /* o  : floating pitch for each subframe */
    float *voice_factors   /* o  : voicing factors                  */
    ,short VBR_cng_reset_flag
)
{
    init_gp_clip(st->clip_var);
    mvr2r( UVWB_Ave, st->mem_AR, M );
    set_f(st->mem_MA, 0, M );
    mem->mem_w0 = 0.0f;
    mem->tilt_code = 0.0f;
    mem->gc_threshold = 0.0f;
    if( VBR_cng_reset_flag )
    {
        set_f( mem->mem_syn, 0, M );
    }
    set_f( st->dispMem, 0, 8 );

    /* last good received frame for FEC in ACELP */
    st->clas = UNVOICED_CLAS;
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

    /* Reset active frame counter */
    st->act_cnt2 = 0;

    /* deactivate bass post-filter */
    st->bpf_off = 1;

    return;
}
