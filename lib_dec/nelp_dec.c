/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <math.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"

/*-------------------------------------------------------------------*
 * nelp_decoder()
 *
 * NELP decoder
 *-------------------------------------------------------------------*/

void nelp_decoder(
    Decoder_State *st,           /* i/o: decoder static memory      */
    float *exc_nelp,     /* o  : adapt. excitation/total exc*/
    float *exc,          /* o  : adapt. excitation exc      */
    short bfi,           /* i  : frame error rate           */
    const short coder_type,    /* i  : coding type                */
    float *gain_buf
)
{
    int i, fid = 0;
    float ptr[L_FRAME], filtRes[L_FRAME], gain_fac;
    float Gains[10],Gain,E3,E2,R;
    float ptr_tmp[L_FRAME];
    int iG1, iG2[2];

    if ((st->last_nelp_mode_dec == 1) && (st->bwidth != st->last_bwidth))
    {
        st->last_nelp_mode_dec = 0;
    }

    if ( coder_type == UNVOICED && st->bwidth == NB )
    {
        if ( st->last_nelp_mode_dec != 1 )
        {
            set_f( st->bp1_filt_mem_nb_dec, 0, 7*2 );
        }
    }
    else if ( coder_type == UNVOICED && (st->bwidth == WB || st->bwidth == SWB) )
    {
        if ( st->last_nelp_mode_dec != 1 )
        {
            set_f( st->bp1_filt_mem_wb_dec, 0, 4*2 );
        }
    }

    if (st->last_nelp_mode_dec != 1)
    {
        set_f( st->shape1_filt_mem_dec, 0, 20 );
        set_f( st->shape2_filt_mem_dec, 0, 20 );
        set_f( st->shape3_filt_mem_dec, 0, 20 );
    }

    if (bfi == 0)
    {
        if(st->rf_frame_type == RF_NELP && st->use_partial_copy)
        {
            iG1 = st->rf_indx_nelp_iG1;
            iG2[0] = st->rf_indx_nelp_iG2[0];
            iG2[1] = st->rf_indx_nelp_iG2[1];
        }
        else
        {
            /* Do Unvoiced/NELP Decoding */
            iG1 = get_next_indice( st, 5 );
            iG2[0] = get_next_indice( st, 6 );
            iG2[1] = get_next_indice( st, 6 );
        }

        if ( coder_type == UNVOICED && (st->bwidth == WB || st->bwidth == SWB) )
        {
            if( st->rf_frame_type == RF_NELP && st->use_partial_copy )
            {
                fid = st->rf_indx_nelp_fid;
            }
            else
            {
                fid = get_next_indice( st, 2 );
            }
        }

        dequantize_uvg( iG1, iG2, Gains, st->bwidth );
    }
    else
    {
        for (i=1,Gain=0.001f; i<=L_SUBFR; i++)
        {
            Gain += SQR(exc[-i]);
        }

        Gain = (float) (sqrt(Gain/L_SUBFR));
        Gain *= 0.8f;/* Some scale down of energy since it is an erasure */

        set_f(Gains, Gain, 10);
    }

    if ( coder_type == UNVOICED && (st->bwidth == WB || st->bwidth == SWB) )
    {
        gain_fac = 1.16f;
    }
    else
    {
        gain_fac = 1.37f;
    }

    generate_nelp_excitation( &(st->nelp_dec_seed), Gains, ptr, gain_fac );

    if ( coder_type == UNVOICED && (st->bwidth == WB || st->bwidth == SWB) )
    {
        polezero_filter( ptr, ptr_tmp, L_FRAME, bp1_num_coef_wb, bp1_den_coef_wb, 4, st->bp1_filt_mem_wb_dec );
        mvr2r( ptr_tmp, ptr, L_FRAME );
    }

    if ( coder_type == UNVOICED && st->bwidth == NB )
    {
        polezero_filter( ptr, ptr_tmp, L_FRAME, bp1_num_coef_nb_fx_order7, bp1_den_coef_nb_fx_order7, 7, st->bp1_filt_mem_nb_dec );
        mvr2r(ptr_tmp,ptr,L_FRAME);
    }

    for( i=0, E3=0.001f; i<L_FRAME; i++ )
    {
        E3 += SQR(ptr[i]);
    }

    if ( coder_type == UNVOICED && (st->bwidth == WB || st->bwidth == SWB) )
    {
        polezero_filter( ptr, ptr_tmp, L_FRAME, shape1_num_coef, shape1_den_coef, 10, st->shape1_filt_mem_dec );
        mvr2r( ptr_tmp, ptr, L_FRAME );

        switch(fid)
        {
        case 1:
            /* Update other filter memory */
            polezero_filter( ptr, filtRes, L_FRAME, shape3_num_coef, shape3_den_coef, 10, st->shape3_filt_mem_dec );

            /* filter the residual to desired shape */
            polezero_filter( ptr, ptr_tmp, L_FRAME, shape2_num_coef, shape2_den_coef, 10, st->shape2_filt_mem_dec );
            mvr2r( ptr_tmp, ptr, L_FRAME );

            break;
        case 2:
            /* Update other filter memory */
            polezero_filter( ptr, filtRes, L_FRAME, shape2_num_coef, shape2_den_coef, 10, st->shape2_filt_mem_dec );

            /* filter the residual to desired shape */
            polezero_filter( ptr, ptr_tmp, L_FRAME, shape3_num_coef, shape3_den_coef, 10, st->shape3_filt_mem_dec );

            mvr2r( ptr_tmp, ptr, L_FRAME );

            break;
        default:
            /* Update other filter memory */
            polezero_filter( ptr, filtRes, L_FRAME, shape2_num_coef, shape2_den_coef, 10, st->shape2_filt_mem_dec );

            polezero_filter( ptr, filtRes, L_FRAME, shape3_num_coef, shape3_den_coef, 10, st->shape3_filt_mem_dec );

            break;
        }

        for (i=0, E2=0.001f; i<L_FRAME; i++)
        {
            E2 += SQR(ptr[i]);
        }

        R = (float)sqrt(E3/E2);

        for (i=0; i<L_FRAME; i++)
        {
            ptr[i] *= R;
        }
    }

    mvr2r( ptr, exc_nelp, L_FRAME );
    set_f( gain_buf, 0.f, NB_SUBFR16k );

    return;
}
