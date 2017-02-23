/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "prot.h"

/*------------------------------------------------------------------
 * syn_12k8()
 *
 * perform the synthesis filtering 1/A(z)
 *------------------------------------------------------------------*/

void syn_12k8(
    const short L_frame,  /* i  : length of the frame                        */
    const float *Aq,      /* i  : LP filter coefficients                     */
    const float *exc,     /* i  : input signal                               */
    float *synth,   /* o  : output signal                              */
    float *mem,     /* i/o: initial filter states                      */
    const short update_m  /* i  : update memory flag: 0 --> no memory update */
)                         /*                          1 --> update of memory */
{
    const float *p_Aq;
    short i_subfr;

    p_Aq = Aq;       /* pointer to interpolated LPC parameters */
    for ( i_subfr=0; i_subfr<L_frame; i_subfr+=L_SUBFR )
    {
        syn_filt( p_Aq, M, &exc[i_subfr], &synth[i_subfr], L_SUBFR, mem, update_m );
        p_Aq  += (M+1);    /* interpolated LPC parameters for next subframe */
    }

    return;
}
