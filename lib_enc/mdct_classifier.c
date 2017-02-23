/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "prot.h"

/*--------------------------------------------------------------------------*
 * mdct_classifier()
 *
 * MDCT signal classifier for HQ_CORE/TCX_20_CORE
 *--------------------------------------------------------------------------*/

#define SMOOTH_FILT_COEFF 0.8f
#define THRESH_UP 1.6f
#define THRESH_DOWN 1.1f
#define HQ_LOCAL 3
#define TCX_LOCAL 1

short mdct_classifier(      /* o  : HQ_CORE/TCX_20_CORE decision            */
    const float *Y,         /* i  : FFT spectrum from fft_rel               */
    Encoder_State *st,      /* i/o: Encoder state variable                  */
    short vadflag,
    float *cldfbBuf_Ener    /* i  : CLDFB energies                          */
)
{
    short c;
    float X[129];
    short k;
    float y[2];
    float nf;
    float pe;
    short np;
    float max_cand;
    short max_i;
    float p_energy;
    float n_energy;
    short d_acc;
    short pos_last;
    float clas_sec;
    short clas_final;
    short i;
    float gain1, gain2, gain3, gain11, gain4;
    float peak_l, peak_h, avrg_l, avrg_h, peak_H1, avrg_H1, peak_H2, avrg_H2;
    short condition1, condition2;
    short condition3, condition4;

    for (k = 0; k < 127; k++)
    {
        X[1 + k] = Y[1 + k] * Y[1 + k] + Y[255 - k] * Y[255 - k];
    }

    for (k = 0; k < 2; k++)
    {
        y[k] = Y[k << 7] * Y[k << 7];
    }

    for (k = 0; k < 2; k++)
    {
        X[k << 7] = y[k];
    }

    nf = X[0];
    pe = X[0];
    np = 0;
    max_cand = -1.0F;
    max_i = 0;
    p_energy = 0.0F;
    n_energy = 0.0F;
    d_acc = 0;
    pos_last = -1;

    for (k = 0; k < 128; k++)
    {
        if (X[k + 1] > nf)
        {
            nf = 0.9578F * nf + 0.0422F * X[k + 1];
        }
        else
        {
            nf = 0.6472F * nf + 0.3528F * X[k + 1];
        }

        if (X[k + 1] > pe)
        {
            pe = 0.42237F * pe + 0.57763F * X[k + 1];
        }
        else
        {
            pe = 0.80285F * pe + 0.19715F * X[k + 1];
        }

        if (X[k + 1] > pe * 0.64F)
        {
            if (X[k + 1] > max_cand)
            {
                max_cand = X[k + 1];
                max_i = (short)(2 + k);
            }
        }
        else
        {
            if (max_i > 0)
            {
                if (np > 0)
                {
                    d_acc = (short)((short)(d_acc + max_i) - pos_last);
                }
                np++;
                pos_last = max_i;
            }

            max_cand = -1.0F;
            max_i = 0;
        }

        p_energy += pe * pe;
        n_energy += nf * nf;
    }

    if (np > 1)
    {
        nf = (float)d_acc / ((float)np - 1.0F);
    }
    else
    {
        nf = 0.0F;
    }
    gain1 = 0.0f;
    gain2 = 0.0f;
    gain3 = 0.0f;
    for(i=0; i<8; i++)
    {
        gain1 += cldfbBuf_Ener[i]/8;
        gain2 += cldfbBuf_Ener[8+i]/8;
        gain3 += cldfbBuf_Ener[16+i]/8;
    }

    gain11 = 8*(gain1 - cldfbBuf_Ener[0]/8)/7;
    gain4 = 0.0f;
    for(i=0; i<12; i++)
    {
        gain4 += cldfbBuf_Ener[12+i]/12;
    }


    peak_H1 = cldfbBuf_Ener[25];
    avrg_H1 = cldfbBuf_Ener[25];
    for(i=1; i<5; i++)
    {
        if(cldfbBuf_Ener[25+i] > peak_H1)
        {
            peak_H1 = cldfbBuf_Ener[25+i];
        }
        avrg_H1 += cldfbBuf_Ener[25+i];
    }

    peak_H2 = cldfbBuf_Ener[20];
    avrg_H2 = cldfbBuf_Ener[20];
    for(i=1; i<5; i++)
    {
        if(cldfbBuf_Ener[20+i] > peak_H2)
        {
            peak_H2 = cldfbBuf_Ener[20+i];
        }
        avrg_H2 += cldfbBuf_Ener[20+i];
    }


    peak_l = 0.0f;
    avrg_l = EPSILON;
    peak_h = 0.0f;
    avrg_h = EPSILON;
    for(i=0; i<32; i++)
    {
        avrg_l += X[20+i];
        avrg_h += X[96+i];
        if(X[20+i] > peak_l)
        {
            peak_l = X[20+i];
        }
        if(X[96+i] > peak_h)
        {
            peak_h = X[96+i];
        }
    }

    condition1 = nf > 12.0F;
    condition2 = p_energy - n_energy * 147.87276f > 0;
    condition3 = gain3 > 1.2f*gain2 ||
                 (gain3 >= 0.8f*gain2 && 5*peak_H1 > 2.0f*avrg_H1) ||
                 (2.6f*peak_l*avrg_h < peak_h *avrg_l || peak_l*avrg_h > 2.6f*peak_h *avrg_l);
    condition4 = (gain4 > 0.8f*gain11 && 2.56f*peak_l*avrg_h > peak_h *avrg_l && peak_l*avrg_h < 5.12f*peak_h *avrg_l) ||
                 (gain4 > 0.3f*gain11 && 32*peak_h < 1.5f*avrg_h && 5*peak_H2 < 1.5f*avrg_H2) ||
                 (2.56f*peak_l*avrg_h < peak_h *avrg_l && 32*peak_h > 1.5f*avrg_h) || (peak_l*avrg_h > 2.56f*peak_h *avrg_l && 32*peak_h < 1.5f*avrg_h);

    if ( ( st->total_brate == HQ_32k && ((!condition1 && condition2) || (condition1 && !condition2)|| condition3) )
            || (st->total_brate == HQ_24k40 && condition4))
    {
        c = HQ_LOCAL;
    }
    else
    {
        c = TCX_LOCAL;
    }

    /* Smooth decision from instantaneous decision*/
    clas_sec = (SMOOTH_FILT_COEFF * st->clas_sec_old) + ((1-SMOOTH_FILT_COEFF) * c);

    /* Do thresholding with hysteresis */
    if((st->clas_final_old == HQ_CORE || st->clas_final_old == TCX_20_CORE)
            && ((st->last_gain1 > 0.5f*gain1 && st->last_gain1 < 2.0f*gain1) && (st->last_gain2 > 0.5f*gain2 && st->last_gain2 < 2.0f*gain2)))
    {
        clas_final = st->clas_final_old;
    }
    else if (clas_sec > st->clas_sec_old && clas_sec > THRESH_UP) /* Going up? */
    {
        clas_final = HQ_CORE;
    }
    else if (clas_sec < THRESH_DOWN)/* Going down */
    {
        clas_final = TCX_20_CORE;
    }
    else
    {
        clas_final = st->clas_final_old;
    }

    /* Prevent the usage of HQ_CORE on noisy-speech or inactive */
    if ( st->mdct_sw_enable==MODE2 && (st->flag_noisy_speech_snr == 1 || vadflag == 0 ) && clas_final == HQ_CORE )
    {
        clas_final = TCX_20_CORE;
    }


    /* Memory update */
    st->clas_sec_old = clas_sec;
    st->clas_final_old = clas_final;
    st->last_gain1 = gain1;
    st->last_gain2 = gain2;

    return clas_final;
}


