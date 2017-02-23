/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"

/*-------------------------------------------------------------------*
 * encod_nelp()
 *
 * Encode Unvoiced frames in SC-VBR
 *-------------------------------------------------------------------*/

void encod_nelp(
    Encoder_State *st,              /* i/o: state structure                         */
    LPD_state *mem,             /* i/o: encoder memories                        */
    const float *speech,          /* i  : input speech                            */
    const float Aw[],             /* i  : weighted A(z) unquantized for subframes */
    const float *Aq,              /* i  : 12k8 Lp coefficient                     */
    float *res,             /* o  : residual signal                         */
    float *synth,           /* o  : core synthesis                          */
    float *tmp_noise,       /* o  : long-term noise energy                  */
    float *exc,             /* i/o: current non-enhanced excitation         */
    float *exc2,            /* i/o: current enhanced excitation             */
    float *pitch_buf,       /* o  : floating pitch values for each subframe */
    float *voice_factors,   /* o  : voicing factors                         */
    float *bwe_exc          /* o  : excitation for SWB TBE                  */
)
{
    float xn[L_SUBFR];      /* Target vector for pitch search               */
    float h1[L_SUBFR];      /* Impulse response vector                      */
    float exc_nelp[L_FRAME];
    const float *p_Aw,*p_Aq;/* pointer to LP filter coeff. vector           */
    short i_subfr, j;


    short reduce_gains = 0;

    if ( st->bwidth == NB && st->input_Fs >= 16000)
    {
        if (st->last_nelp_mode == 0)
        {
            set_f( st->nelp_lp_fit_mem, 0, NELP_LP_ORDER*2 );
        }
        polezero_filter( res, res, L_FRAME, num_nelp_lp, den_nelp_lp, NELP_LP_ORDER, st->nelp_lp_fit_mem );  /*16-Q of filter coeff*/
    }

    p_Aw = Aw;
    p_Aq = Aq;

    for (i_subfr=0; i_subfr<L_FRAME; i_subfr+=L_SUBFR)
    {
        /*----------------------------------------------------------------*
         * - Find the excitation search target "xn" and innovation
         * target in residual domain "cn"
         * - Compute impulse response, h1[], of weighted synthesis filter
         *----------------------------------------------------------------*/
        mvr2r( &res[i_subfr], &exc[i_subfr], L_SUBFR );

        find_targets( speech, mem->mem_syn, i_subfr, &mem->mem_w0, p_Aq,res, L_SUBFR, p_Aw, TILT_FAC, xn, NULL, h1 );

        if (i_subfr == 0)
        {
            if ( (st->Local_VAD == 1 ) && ( st->bwidth == NB) )
            {
                reduce_gains = 1;
            }

            nelp_encoder(st, res, exc_nelp, reduce_gains);
        }
        *tmp_noise = 0;

        /*-----------------------------------------------------------------*
         * Synthesize speech to update mem_syn[].
         * Update A(z) filters
         *-----------------------------------------------------------------*/

        syn_filt( p_Aq, M, &exc_nelp[i_subfr], &synth[i_subfr], L_SUBFR, mem->mem_syn, 1 );

        p_Aw += (M+1);
        p_Aq += (M+1);
        *pitch_buf = L_SUBFR;
        pitch_buf++;
    }

    mvr2r( exc_nelp, exc, L_FRAME );

    /*-----------------------------------------------------------------*
     * Updates: last value of new target is stored in mem_w0
     *-----------------------------------------------------------------*/

    mem->mem_w0 = xn[L_SUBFR-1] - (exc[L_FRAME-1]);
    mem->tilt_code = 0.0f; /* purely unvoiced */
    st->prev_tilt_code = mem->tilt_code;

    mvr2r(exc, exc2, L_FRAME);

    st->prev_ppp_gain_pit = 0.0;
    st->dispMem[0] = 0;
    st->dispMem[2] = st->prev_ppp_gain_pit;

    for( j=3; j<7; j++ )
    {
        st->dispMem[j] = st->dispMem[j-1];
    }

    interp_code_5over2( exc2, bwe_exc, L_FRAME );
    set_f( voice_factors, 0.0f, NB_SUBFR16k );

    return;
}
