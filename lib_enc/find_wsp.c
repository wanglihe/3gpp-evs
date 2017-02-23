/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "prot.h"

/*-------------------------------------------------------------------*
 * find_wsp()
 *
 * Compute weighted speech used in open-loop pitch search
 *-------------------------------------------------------------------*/

void find_wsp(
    const short L_frame,     /* i  : length of the frame                   */
    const short L_subfr,     /* i  : length of subframe                    */
    const short nb_subfr,    /* i  : number of subframes                   */
    const float *A,          /* i  : A(z) filter coefficients              */
    float *Aw,         /* o  : weighted A(z) filter coefficients     */
    const float *speech,     /* i  : pointer to the denoised speech frame  */
    const float tilt_fact,   /* i  : tilt factor                           */
    float *wsp,        /* o  : poitnter to the weighted speech frame */
    float *mem_wsp,    /* i/o: W(Z) denominator memory               */
    const float gamma,       /* i  : weighting factor                      */
    const short L_look       /* i  : look-ahead                            */
)
{
    float *p_Aw, tmp;
    short i_subfr;


    /*-----------------------------------------------------------------*
     *  Compute weighted A(z) unquantized for subframes
     *-----------------------------------------------------------------*/

    weight_a_subfr( nb_subfr, A, Aw, gamma, M );

    /*-----------------------------------------------------------------*
     *  Compute weighted speech for all subframes
     *-----------------------------------------------------------------*/

    p_Aw = Aw;
    for (i_subfr=0; i_subfr<L_frame; i_subfr += L_subfr )
    {
        residu( p_Aw, M, &speech[i_subfr], &wsp[i_subfr], L_subfr );
        p_Aw += (M+1);
    }
    p_Aw -= (M+1);

    /*-----------------------------------------------------------------*
     *  Weighted speech computation is extended on look-ahead
     *-----------------------------------------------------------------*/

    deemph( wsp, tilt_fact, L_frame, mem_wsp );
    residu( p_Aw, M, &speech[L_frame], &wsp[L_frame], L_look );
    tmp = *mem_wsp;
    deemph( &wsp[L_frame], tilt_fact, L_look, &tmp );

    return;
}
