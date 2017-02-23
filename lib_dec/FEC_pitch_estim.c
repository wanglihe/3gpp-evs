/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/


#include <math.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"


/*------------------------------------------------------------------------*
 * FEC_pitch_estim()
 *
 * Estimation of pitch for FEC
 *------------------------------------------------------------------------*/

void FEC_pitch_estim(
    const short Opt_AMR_WB,                 /* i  : flag indicating AMR-WB IO mode          */
    const short last_core,                  /* i  : last core                               */
    const short L_frame,                    /* i  : length of the frame                     */
    const short clas,                       /* i  : current frame classification            */
    const short last_good,                  /* i  : last good clas information              */
    const float pitch_buf[],                /* i  : Floating pitch   for each subframe      */
    const float old_pitch_buf[],            /* i  : buffer of old subframe pitch values     */
    float *bfi_pitch,                 /* i/o: update of the estimated pitch for FEC   */
    short *bfi_pitch_frame,           /* o  : frame length when pitch was updated     */
    short *upd_cnt,                   /* i/o: update counter                          */
    const short coder_type                  /* i  : coder_type                              */
)
{
    if( last_core == HQ_CORE )
    {
        *bfi_pitch = pitch_buf[(L_frame/L_SUBFR)-1];
        *bfi_pitch_frame = L_frame;
        *upd_cnt = MAX_UPD_CNT;
    }

    if( (clas == VOICED_CLAS && last_good >= VOICED_TRANSITION) || (Opt_AMR_WB && clas == VOICED_CLAS) )
    {
        /* update pitch for FEC if pitch is coherent */
        if( ( (pitch_buf[3] < 1.4f * pitch_buf[1]) && (pitch_buf[3] > 0.7f * pitch_buf[1]) &&
                (pitch_buf[1] < 1.4f * old_pitch_buf[2*NB_SUBFR-1]) && (pitch_buf[1] > 0.7f * old_pitch_buf[2*NB_SUBFR-1]) &&
                (L_frame == L_FRAME) ) ||
                ( (pitch_buf[3] < 1.4f * pitch_buf[1]) && (pitch_buf[3] > 0.7f * pitch_buf[1]) &&
                  (pitch_buf[1] < 1.4f * old_pitch_buf[2*NB_SUBFR16k-1]) && (pitch_buf[1] > 0.7f * old_pitch_buf[2*NB_SUBFR16k-1]) &&
                  (L_frame == L_FRAME16k) ) || ( coder_type == TRANSITION ) )
        {
            *bfi_pitch = pitch_buf[(L_frame/L_SUBFR)-1];
            *bfi_pitch_frame = L_frame;
            *upd_cnt = 0;
        }
    }

    return;
}
