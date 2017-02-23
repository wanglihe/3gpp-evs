/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "prot.h"

/*-------------------------------------------------------------------*
 * find_targets()
 *
 * Find the target vectors for excitaiton search:
 *-------------------------------------------------------------------*/

void find_targets(
    const float *speech,     /* i  : pointer to the speech frame                      */
    const float *mem_syn,    /* i  : memory of the synthesis filter                   */
    const short i_subfr,     /* i  : subframe index                                   */
    float *mem_w0,     /* i/o: weighting filter denominator memory              */
    const float *p_Aq,       /* i  : interpolated quantized A(z) filter               */
    const float *res,        /* i  : residual signal                                  */
    const short L_subfr,     /* i  : length of vectors for gain quantization          */
    const float *Ap,         /* i  : unquantized A(z) filter with bandwidth expansion */
    const float tilt_fac,    /* i  : tilt factor                                      */
    float *xn,         /* o  : Close-loop Pitch search target vector            */
    float *cn,         /* o  : target vector in residual domain                 */
    float *h1          /* o  : impulse response of weighted synthesis filter    */
)
{
    short i;
    float error[M+5*L_SUBFR];       /* error of quantization */
    float tmp_fl[L_SUBFR+M], tmp;

    /*------------------------------------------------------------------------*
     * Find the target vector for excitation search:
     *
     *             |------|  res[n]
     * speech[n]---| A(z) |--------
     *             |------|       |   |--------| error[n]  |------|
     *                   zero -- (-)--| 1/A(z) |-----------| W(z) |-- target
     *                   exc          |--------|           |------|
     *
     * Instead of subtracting the zero-input response of filters from
     * the weighted input speech, the above configuration is used to
     * compute the target vector.
     *-----------------------------------------------------------------------*/

    for ( i=0; i<M; i++ )
    {
        error[i] = speech[i+i_subfr-M] - mem_syn[i];
    }

    syn_filt( p_Aq, M, &res[i_subfr], error+M, L_subfr, error, 0 );
    residu( Ap, M, error+M, xn, L_subfr );
    deemph( xn, tilt_fac, L_subfr, mem_w0 );

    /*-----------------------------------------------------------------*
     * Find target in residual domain (cn[]) for innovation search
     *--------------------------------------------------------------*/
    if( cn!=NULL )
    {
        /* first half: xn[] --> cn[] */
        set_f( tmp_fl, 0, M );
        mvr2r( xn, tmp_fl+M, L_subfr/2 );
        tmp = 0.0f;

        preemph( tmp_fl+M, tilt_fac, L_subfr/2, &tmp );
        syn_filt( Ap, M, tmp_fl+M, tmp_fl+M, L_subfr/2, tmp_fl, 0 );
        residu( p_Aq, M, tmp_fl+M, cn, L_subfr/2 );

        /* second half: res[] --> cn[] (approximated and faster) */
        mvr2r( &res[i_subfr+(L_subfr/2)], cn+(L_subfr/2), L_subfr/2 );
    }

    /*-----------------------------------------------------------------*
     * Compute impulse response h1[] of the weighted synthesis filter
     *-----------------------------------------------------------------*/

    set_f( h1, 0, L_subfr );
    mvr2r( Ap, h1, M+1 );
    syn_filt( p_Aq, M, h1, h1, L_subfr, h1+(M+1), 0 );
    tmp = 0.0f;
    deemph( h1, tilt_fac, L_subfr, &tmp );

    return;
}
