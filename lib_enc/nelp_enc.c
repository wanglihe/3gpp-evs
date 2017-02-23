/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"

/*-------------------------------------------------------------------*
 * quantize_uvg()
 *
 * Quantize unvoiced gains
 *--------------------------------------------------------------------*/

static void quantize_uvg(
    float *G,             /* i : unvoiced gain     */
    int *iG1,             /* i : gain 1 index      */
    int *iG2,             /* i : gain 2 index      */
    float *quantG,        /* o : quantized gain    */
    short bwidth
)
{
    float G1[2], G2[10];
    int i, j, k;
    float mse, mmse;
    const float (*UVG1CB)[2] = NULL;
    const float (*UVG2CB1)[5] = NULL;
    const float (*UVG2CB2)[5] = NULL;

    if ( bwidth == NB )
    {
        UVG1CB = UVG1CB_NB;
        UVG2CB1 = UVG2CB1_NB;
        UVG2CB2 = UVG2CB2_NB;
    }
    else if ( bwidth == WB || bwidth == SWB )
    {
        UVG1CB = UVG1CB_WB;
        UVG2CB1 = UVG2CB1_WB;
        UVG2CB2 = UVG2CB2_WB;
    }

    for (i=0; i<2; i++)
    {
        G1[i] = 0;
        for (j=0; j<5; j++)
        {
            G1[i] += SQR(G[i*5+j]);
        }
        G1[i] = (float) log10(sqrt(G1[i]/5));
    }

    mmse = (float) 1e30;
    *iG1 = 0;

    for (i=0; i< UVG1_CBSIZE; i++)
    {
        mse = SQR(G1[0]-(UVG1CB[i][0]))+SQR(G1[1]-(UVG1CB[i][1]));
        if (mse < mmse)
        {
            *iG1 = i;
            mmse = mse;
        }
    }

    G1[0] = (float) pow(10.0,(UVG1CB[*iG1][0]));
    G1[1] = (float) pow(10.0,(UVG1CB[*iG1][1]));

    for (i=0; i<2; i++)
    {
        for (j=0; j<5; j++)
        {
            G2[i*5+j] = G[i*5+j]/G1[i];
        }
    }

    for (i=0; i<2; i++)
    {
        mmse = (float) 1e30;
        iG2[i] = 0;
        for (j=0; j<UVG2_CBSIZE; j++)
        {
            mse = 0;
            for (k=0; k<5; k++)
            {
                if (i == 0)
                {
                    mse += SQR(G2[i*5+k]-UVG2CB1[j][k]);
                }
                else if (i == 1)
                {
                    mse += SQR(G2[i*5+k]-UVG2CB2[j][k]);
                }
            }

            if (mse < mmse)
            {
                mmse = mse;
                iG2[i] = j;
            }
        }
    }

    mvr2r(G, G2, 10);

    dequantize_uvg( *iG1, iG2, quantG, bwidth );

    return;
}


/*-------------------------------------------------------------------*
 * nelp_encoder()
 *
 * NELP encoder
 *--------------------------------------------------------------------*/

void nelp_encoder(
    Encoder_State *st,         /* i/o: encoder state                      */
    float *in,         /* i  : residual signal                    */
    float *exc         /* o  : NELP quantized excitation signal   */
    , short reduce_gains
)
{
    int i,j;
    float *ptr = exc;
    int lag = 25; /* to cover 25*9 + 31 */
    float Gains[10], gain_fac;
    int iG1, iG2[2], fid;
    float fdbck, var_dB, tmp;
    float E1 = 0, E2, E3, R, EL1 = 0, EH1 = 0, EL2, EH2, RL, RH;
    float filtRes[L_FRAME];
    float ptr_tmp[L_FRAME];
    short rf_flag;

    rf_flag = st->rf_mode;

    if ((st->last_nelp_mode == 1) && (st->bwidth != st->last_bwidth))
    {
        st->last_nelp_mode = 0;
    }

    if ( st->bwidth == NB )
    {
        if (st->last_nelp_mode != 1)
        {
            set_f(st->bp1_filt_mem_nb, 0, 7*2);
        }
    }
    else if ( st->bwidth == WB || st->bwidth == SWB )
    {
        if (st->last_nelp_mode != 1)
        {
            set_f(st->bp1_filt_mem_wb, 0, 4*2);
        }
    }

    if (st->last_nelp_mode != 1)
    {
        if ( st->bwidth == WB || st->bwidth == SWB )
        {
            set_f(st->shape1_filt_mem, 0, 20);
            set_f(st->shape2_filt_mem, 0, 20);
            set_f(st->shape3_filt_mem, 0, 20);
            set_f(st->txlpf1_filt1_mem, 0, 20);
            set_f(st->txlpf1_filt2_mem, 0, 20);
            set_f(st->txhpf1_filt1_mem, 0, 20);
            set_f(st->txhpf1_filt2_mem, 0, 20);
        }
    }

    /* Start Unvoiced/NELP Processing */
    if ( st->bwidth == WB || st->bwidth == SWB )
    {
        for (i=0, E1=0.001f; i<L_FRAME; i++)
        {
            E1 += SQR(in[i]);
        }

        polezero_filter(in, filtRes, L_FRAME, txlpf1_num_coef, txlpf1_den_coef, 10, st->txlpf1_filt1_mem);

        for (i=0, EL1=0.001f; i<L_FRAME; i++)
        {
            EL1 += SQR(filtRes[i]);
        }

        polezero_filter(in, filtRes, L_FRAME, txhpf1_num_coef, txhpf1_den_coef, 10, st->txhpf1_filt1_mem);

        for (i=0, EH1=0.001f; i<L_FRAME; i++)
        {
            EH1 += SQR(filtRes[i]);
        }
    }

    for (i=0; i<9; i++)
    {
        for (j = i*lag, Gains[i] = 0.001f; j<(i+1)*lag; j++)
        {
            Gains[i] += SQR(in[j]);
        }

        Gains[i] = (float) sqrt(Gains[i]/lag);
    }

    for (j = i*lag, Gains[i] = 0.001f; j<L_FRAME; j++)
    {
        Gains[i] += SQR(in[j]);
    }

    Gains[i] = (float) sqrt(Gains[i]/(L_FRAME-(lag*i)));

    if (reduce_gains==1)
    {
        for (i=0; i<10; i++)
        {
            Gains[i] = ( Gains[i]*0.6 );
        }
    }

    if (st->last_nelp_mode != 1) /* if prev frame was not NELP then init mem*/
    {
        st->nelp_gain_mem = Gains[0];
    }

    tmp = (float) (20.0 * (log10 (Gains[0]) - log10 (st->nelp_gain_mem) ) );
    var_dB = tmp * tmp;
    for (i = 1; i < 10; i++)
    {
        tmp = (float) (20.0 * (log10 (Gains[i]) - log10 (Gains[i - 1])));
        var_dB += tmp * tmp;
    }

    if (st->last_nelp_mode!=1)
    {
        var_dB *= 0.111f;
    }
    else
    {
        var_dB *= 0.1f;
    }

    fdbck = (float) (0.82f / (1.0f + exp (0.25f * (var_dB - 20.0f))));

    for (i = 0; i < 10; i++)
    {
        Gains[i] = (float)((1.0f - fdbck) * Gains[i] + fdbck * st->nelp_gain_mem);
        st->nelp_gain_mem = Gains[i];
    }

    quantize_uvg(Gains, &iG1, iG2, Gains, (short) st->bwidth );

    if( rf_flag )
    {
        st->rf_indx_nelp_iG1[0] = (short) iG1;
        st->rf_indx_nelp_iG2[0][0] = (short) iG2[0];
        st->rf_indx_nelp_iG2[0][1] = (short) iG2[1];
    }
    else
    {
        push_indice( st, IND_IG1, iG1, 5 );
        push_indice( st, IND_IG2A, iG2[0], 6 );
        push_indice( st, IND_IG2B, iG2[1], 6 );
    }

    if ( st->bwidth == WB || st->bwidth == SWB )
    {
        gain_fac = 1.16f;
    }
    else
    {
        gain_fac = 1.37f;
    }

    generate_nelp_excitation( &(st->nelp_enc_seed), Gains, ptr, gain_fac );

    if ( st->bwidth == WB || st->bwidth == SWB )
    {
        polezero_filter(ptr,ptr_tmp,L_FRAME,bp1_num_coef_wb,bp1_den_coef_wb,4,st->bp1_filt_mem_wb);
        mvr2r(ptr_tmp,ptr,L_FRAME);
    }
    else if ( st->bwidth == NB )
    {
        polezero_filter(ptr,ptr_tmp,L_FRAME,bp1_num_coef_nb_fx_order7,bp1_den_coef_nb_fx_order7,7,st->bp1_filt_mem_nb);
        mvr2r(ptr_tmp,ptr,L_FRAME);
    }

    for (i=0, E3=0.001f; i<L_FRAME; i++)
    {
        E3 += SQR(ptr[i]);
    }

    if ( st->bwidth == WB || st->bwidth == SWB )
    {
        polezero_filter(ptr,ptr_tmp,L_FRAME,shape1_num_coef,shape1_den_coef,10, st->shape1_filt_mem);
        mvr2r(ptr_tmp,ptr,L_FRAME);

        for (i=0, E2=0.001f; i<L_FRAME; i++)
        {
            E2 += SQR(ptr[i]);
        }

        R = (float) sqrt(E1/E2);

        for (i=0; i<L_FRAME; i++)
        {
            filtRes[i] = R*ptr[i];
        }

        polezero_filter(filtRes,ptr_tmp,L_FRAME,txlpf1_num_coef,txlpf1_den_coef,10,st->txlpf1_filt2_mem);
        mvr2r(ptr_tmp,filtRes,L_FRAME);

        for (i=0, EL2=0.001f; i<L_FRAME; i++)
        {
            EL2 += SQR(filtRes[i]);
        }

        for (i=0; i<L_FRAME; i++)
        {
            filtRes[i] = R*ptr[i];
        }

        polezero_filter(filtRes,ptr_tmp,L_FRAME,txhpf1_num_coef,txhpf1_den_coef,10,st->txhpf1_filt2_mem);
        mvr2r(ptr_tmp,filtRes,L_FRAME);

        for (i=0, EH2=0.001f; i<L_FRAME; i++)
        {
            EH2 += SQR(filtRes[i]);
        }

        RL = (float) 10.0f * (float)log10(EL1/EL2);
        RH = (float) 10.0f * (float)log10(EH1/EH2);

        fid = 0;
        if (RL < -3)
        {
            fid = 1;
        }
        else if (RH < -3)
        {
            fid = 2;
        }

        if(rf_flag==0)
        {

            switch(fid)
            {
            case 1:
                /* Update other filter memory */
                polezero_filter(ptr,filtRes,L_FRAME,shape3_num_coef,shape3_den_coef,10,st->shape3_filt_mem);

                /* filter the residual to desired shape */
                polezero_filter(ptr,ptr_tmp,L_FRAME,shape2_num_coef,shape2_den_coef,10,st->shape2_filt_mem);

                mvr2r(ptr_tmp,ptr,L_FRAME);

                break;
            case 2:
                /* Update other filter memory */
                polezero_filter(ptr,filtRes,L_FRAME,shape2_num_coef,shape2_den_coef,10,st->shape2_filt_mem);

                /* filter the residual to desired shape */
                polezero_filter(ptr,ptr_tmp,L_FRAME,shape3_num_coef,shape3_den_coef,10,st->shape3_filt_mem);

                mvr2r(ptr_tmp,ptr,L_FRAME);

                break;
            default:
                /* Update other filter memory */
                polezero_filter(ptr,filtRes,L_FRAME,shape2_num_coef,shape2_den_coef,10,st->shape2_filt_mem);
                polezero_filter(ptr,filtRes,L_FRAME,shape3_num_coef,shape3_den_coef,10,st->shape3_filt_mem);

                break;
            }

            for (i=0, E2=0.001f; i<L_FRAME; i++)
            {
                E2 += SQR(ptr[i]);
            }

            R = (float) sqrt(E3/E2);
            for (i=0; i<L_FRAME; i++)
            {
                ptr[i] *= R;
            }
        }

        if(rf_flag)
        {
            st->rf_indx_nelp_fid[0] = (short)fid;
        }
        else
        {
            push_indice( st, IND_NELP_FID, fid, 2 );
        }
    }

    if(rf_flag==0)
    {
        for (i = 0; i < L_FRAME; i++)
        {
            exc[i] = ptr[i];
        }
    }
    return;
}
