/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "options.h"
#include "prot.h"
#include "rom_com.h"


/*-------------------------------------------------------------------*
 * enc_acelp_tcx_main()
 *
 * encoder function for coding ACELP/TCX
 *--------------------------------------------------------------------*/

void enc_acelp_tcx_main(
    const float new_samples[],          /* i  : new samples                         */
    Encoder_State *st,                    /* i/o: encoder state structure             */
    const short coder_type,             /* i  : coding type                         */
    const short pitch[3],               /* i  : open-loop pitch values for quantiz. */
    const float voicing[3],             /* i  : open-loop pitch gains               */
    float Aw[NB_SUBFR16k*(M+1)],  /* i  : weighted A(z) unquant. for subframes*/
    const float lsp_new[M],             /* i  : LSPs at the end of the frame        */
    const float lsp_mid[M],             /* i  : LSPs at the middle of the frame     */
    HANDLE_FD_CNG_ENC hFdCngEnc,        /* i/o: FD CNG handle                       */
    float bwe_exc_extended[],     /* i/o: bandwidth extended excitation       */
    float *voice_factors,         /* o  : voicing factors                     */
    float pitch_buf[],            /* o  : floating pitch for each subframe    */
    short vad_hover_flag
)
{

    float old_bwe_exc[(PIT16k_MAX + (L_FRAME16k + 1) + L_SUBFR16k) * 2]; /* excitation buffer */
    float *ptr_bwe_exc;              /* pointer to BWE excitation signal in the current frame */

    ptr_bwe_exc = old_bwe_exc + PIT16k_MAX * 2;

    if (st->last_core == ACELP_CORE)
    {
        set_f( old_bwe_exc + PIT16k_MAX * 2, 0.f, ((L_FRAME16k + 1) + L_SUBFR16k) * 2 );
        mvr2r( st->old_bwe_exc, old_bwe_exc, PIT16k_MAX * 2 );
    }
    else
    {
        set_f( old_bwe_exc, 0.f, ((L_FRAME16k + 1) + L_SUBFR16k + PIT16k_MAX) * 2 );
    }


    /* Guided ACELP PLC */
    gPLC_encInfo( &st->plcExt, st->total_brate, st->bwidth, st->clas, coder_type );

    if( st->core_brate != FRAME_NO_DATA && st->core_brate != SID_2k40 )
    {
        /* Run Core Coder */
        if( st->tcxonly == 0 )
        {
            core_encode_openloop( st, coder_type, pitch, voicing, Aw, lsp_new, lsp_mid, pitch_buf, voice_factors, ptr_bwe_exc, vad_hover_flag );
        }
        else
        {
            core_encode_twodiv( new_samples, st, coder_type, pitch, voicing, Aw );
        }

        /* Apply non linearity to the SHB excitation */
        if( st->core == ACELP_CORE && st->igf )
        {
            non_linearity( ptr_bwe_exc, bwe_exc_extended, st->old_bwe_exc_extended, L_FRAME32k, &st->bwe_non_lin_prev_scale, coder_type, voice_factors, st->L_frame );

            /* update the old_BWE_exc memory */
            mvr2r( &old_bwe_exc[L_FRAME32k], st->old_bwe_exc, PIT16k_MAX * 2 );
        }
        else
        {
            set_f( st->old_bwe_exc_extended, 0, NL_BUFF_OFFSET );
            set_f( st->old_bwe_exc, 0, PIT16k_MAX * 2 );   /* reset old non_linear exc during igf frames */
            st->bwe_non_lin_prev_scale = 0.0f;
        }
    }
    else
    {
        /* Run SID Coder */
        if( st->core_brate == SID_2k40 )
        {
            FdCng_encodeSID( hFdCngEnc, st, st->preemph_fac );
        }

        /* Generate Comfort Noise */
        generate_comfort_noise_enc( st );

        /* Update Core Encoder */
        core_encode_update_cng( st, hFdCngEnc->hFdCngCom->timeDomainBuffer, hFdCngEnc->hFdCngCom->A_cng, Aw );
    }

    /* coreSwitching update of Mode 1 parameters in the last frame */
    st->last_coder_type = coder_type;


    return;
}
