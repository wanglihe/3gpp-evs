/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <assert.h>
#include "options.h"
#include "prot.h"
#include "cnst.h"
#include "rom_com.h"


/*-------------------------------------------------------------------------*
 * Local functions                                                         *
 *-------------------------------------------------------------------------*/

static int gain_enc_uv( const float *code, int lcode, float *gain_pit, float *gain_code,
                        ACELP_CbkCorr *coeff, float *past_gcode, float *gain_inov );

static int gain_enc_gacelp_uv( const float *code, const float *code2, int lcode, const float mean_ener, float *gain_pit, float *gain_code,
                               float *gain_code2, ACELP_CbkCorr *coeff, float *past_gcode, float *gain_inov, short noisy_speech_flag );


/*-------------------------------------------------------------------------*
 * procedure q_gain2_plus                                                  *
 * ~~~~~~~~~~~~~~~~~~~~~~                                                  *
 * Quantization of pitch and codebook gains.                               *
 * The following routines is Q_gains updated for AMR_WB_PLUS.              *
 * MA prediction is removed and MEAN_ENER is now quantized with 2 bits and *
 * transmitted once every ACELP frame to the gains decoder.                *
 * The pitch gain and the code gain are vector quantized and the           *
 * mean-squared weighted error criterion is used in the quantizer search.  *
 *-------------------------------------------------------------------------*/

void encode_acelp_gains(
    float *code,
    int gains_mode,
    float mean_ener_code,
    short clip_gain,
    ACELP_CbkCorr *g_corr,
    float *gain_pit,
    float *gain_code,
    int **pt_indice,
    float *past_gcode,
    float *gain_inov,
    int L_subfr,
    float *code2,
    float *gain_code2,
    short noisy_speech_flag
)
{
    int index = 0;

    if (((gains_mode>0) && (gains_mode <4)))
    {
        /* Memory-less gain coding */
        index = Mode2_gain_enc_mless( code, L_subfr, gain_pit, gain_code, g_corr, mean_ener_code,
                                      clip_gain, past_gcode, gain_inov,gains_mode-1 );
    }
    else if (gains_mode == 6)
    {
        /* UV gains quantizer (6bits/subfr) */
        index = gain_enc_uv(code, L_subfr, gain_pit, gain_code, g_corr, past_gcode, gain_inov );
    }
    else if (gains_mode == 7)
    {
        /* GACELP_UV gains quantizer (7=5-2bits/subfr) */
        index = gain_enc_gacelp_uv(code, code2, L_subfr, mean_ener_code, gain_pit, gain_code,gain_code2, g_corr, past_gcode, gain_inov, noisy_speech_flag );
    }
    else
    {
        fprintf(stderr, "invalid gains coding for acelp!\n");
        assert(0);
    }

    **pt_indice = index;
    (*pt_indice)++;

    return;
}

/*---------------------------------------------------------------------*
 * procedure Mode2_gain_enc_mless
 * Quantization of pitch and codebook gains.
 * - an initial predicted gain, gcode0, is first determined based on
 *   the predicted scaled innovation energy
 * - the correction  factor gamma = g_code / gcode0 is then vector quantized
 *   along with gain_pit
 * - the mean-squared weighted error criterion is used for the quantizer search
 *---------------------------------------------------------------------*/

int Mode2_gain_enc_mless(
    const float *code,     /* i    : algebraic excitation                                            */
    int lcode,             /* (i)  : Subframe size                                                   */
    float *gain_pit,       /* o    : quantized pitch gain                                            */
    float *gain_code,      /* o    : quantized codebook gain                                         */
    ACELP_CbkCorr *pcoeff, /* i/o  : correlations <y1,y1>, -2<xn,y 1>,<y2,y2>, -2<xn,y2> and 2<y1,y2>*/
    float mean_ener,       /* (i)  : mean_ener defined in open-loop (3 bits)                         */
    const short clip_gain, /* i    : gain pitch clipping flag (1 = clipping)                         */
    float *past_gcode,     /* (i/o): past gain of code                                               */
    float *gain_inov,      /* (o)  : unscaled innovation gain                                        */
    const short coder_type /* (i)  : coder type                                                      */
)
{
    short index, i, size,size_clip;
    const Word16 *p, *t_qua_gain;
    float dist, dist_min, g_pitch, g_code, gcode0, ener_code;
    ACELP_CbkCorr coeff;

    /*-----------------------------------------------------------------*
     * - calculate the unscaled innovation energy
     * - calculate the predicted gain code
     *-----------------------------------------------------------------*/

    if(coder_type == 0)
    {
        *gain_inov = 1.0f/ (float)sqrt( ( dotp( code, code, lcode ) + 0.01f ) / lcode );
        ener_code = 10 * (float)log10( ( dotp( code, code, lcode ) + 0.01f ) / lcode );
        gcode0 = (float) pow(10, 0.05 * (mean_ener - ener_code));
    }
    else
    {
        ener_code = 0.01F;
        for(i=0; i<lcode; i++)
        {
            ener_code += code[i] * code[i];
        }

        *gain_inov = (float)sqrt((float)lcode / ener_code);

        ener_code = (float)(-10.0 * log10((float)lcode / ener_code));
        gcode0 = mean_ener - ener_code;
        gcode0 = (float)pow(10.0,gcode0/20.0); /* predicted gain */
    }

    /*-----------------------------------------------------------------*
     * gain quantization initializations
     * - find the initial quantization pitch index
     * - set the gains searching range
     *-----------------------------------------------------------------*/

    if( coder_type == 0)
    {
        t_qua_gain = E_ROM_qua_gain5b_const;
        size_clip=9;
        size = NB_QUA_GAIN5B;
    }
    else if(coder_type == 1)
    {
        t_qua_gain = E_ROM_qua_gain6b_const;
        size_clip=6;
        size = NB_QUA_GAIN6B;                            /* searching range of the gain quantizer */
    }
    else
    {
        t_qua_gain = E_ROM_qua_gain7b_const;
        size_clip=21;
        size = NB_QUA_GAIN7B;
    }

    if( clip_gain == 1 )
    {
        size -= size_clip;                       /* limit pitch gain  to 1.0 */
    }

    coeff = *pcoeff;
    coeff.xy1 *= -2.0f;
    coeff.xy2 *= -2.0f;
    coeff.y1y2 *= 2.0f;

    /*-----------------------------------------------------------------*
     * search for the best quantizer
     *-----------------------------------------------------------------*/

    p=t_qua_gain;
    dist_min = FLT_MAX;
    index = 0;

    for (i = 0; i<size; i++)
    {
        g_pitch = (float)((*p++))/(1<<14);         /* pitch gain */
        g_code = gcode0 * (float)((*p++))/(1<<11);           /* codebook gain */

        dist = g_pitch*g_pitch * coeff.y1y1
               + g_pitch         * coeff.xy1
               + g_code*g_code   * coeff.y2y2
               + g_code          * coeff.xy2
               + g_pitch*g_code  * coeff.y1y2;

        if ( dist < dist_min )
        {
            dist_min = dist;
            index = i;
        }
    }

    *gain_pit  = (float)(t_qua_gain[index*2])/(1<<14);
    *gain_code = (float)(t_qua_gain[index*2+1])/(1<<11) * gcode0;

    *past_gcode = *gain_code / *gain_inov;

    return index;
}

/*---------------------------------------------------------------------*
* procedure gain_enc_uv
* Quantization of pitch and codebook gains.
* - an initial predicted gain, gcode0, is first determined based on
*   the predicted scaled innovation energy
* - the correction  factor gamma = g_code / gcode0 is then vector quantized
*   along with gain_pit
* - the mean-squared weighted error criterion is used for the quantizer search
*---------------------------------------------------------------------*/

static int gain_enc_uv(
    const float *code,         /* i    : algebraic excitation                                            */
    int lcode,         /* (i)  : Subframe size                                                   */
    float *gain_pit,     /* o    : quantized pitch gain                                            */
    float *gain_code,    /* o    : quantized codebook gain                                         */
    ACELP_CbkCorr *coeff,        /* i/o  : correlations <y1,y1>, -2<xn,y1>,<y2,y2>, -2<xn,y2> and 2<y1,y2> */
    float *past_gcode,   /* (i/o): past gain of code                                               */
    float *gain_inov     /* (o)  : unscaled innovation gain                                        */
)
{
    short index;
    float g_code,g_code_corr,g_code_nrg;

    /*-----------------------------------------------------------------*
     * - calculate the unscaled innovation energy
     *-----------------------------------------------------------------*/

    *gain_inov = 1.0f/ (float)sqrt( ( dotp( code, code, lcode ) + 0.01f ) / lcode );

    g_code_corr = coeff->xy2/(coeff->y2y2*(*gain_inov));     /*Correlation based*/

    (void) g_code_nrg;
    g_code=g_code_corr;

    /*90dB max of gain for 2^15 amplitude code, */
    if( g_code>0.000001f )
    {
        index=(int)(((20.f*log10(g_code)+30.f)/1.9f)+0.5f);

        if( index>63 )
        {
            index=63;
        }
        else if(index<0)
        {
            index=0;
        }
    }
    else
    {
        index=0;
    }

    *gain_code= (float) pow(10.f,(((index*1.9f)-30.f)/20.f));
    *past_gcode=*gain_code;  /*unscaled gain*/
    *gain_code *= *gain_inov; /*scaled gain*/
    *gain_pit=0.f;

    return index;
}

/*---------------------------------------------------------------------*
* procedure gain_enc_gacelp_uv
* Quantization of pitch and codebook gains.
* - an initial predicted gain, gcode0, is first determined based on
*   the predicted scaled innovation energy
* - the correction  factor gamma = g_code / gcode0 is then vector quantized
*   along with gain_pit
* - the mean-squared weighted error criterion is used for the quantizer search
*---------------------------------------------------------------------*/

static int gain_enc_gacelp_uv(
    const float *code,         /* i    : algebraic excitation                                            */
    const float *code2,        /* i    : gaussian excitation                                             */
    int lcode,         /* (i)  : Subframe size                                                   */
    const float mean_ener,     /* i : quantized mean energy of the frame                                 */
    float *gain_pit,     /* o    : quantized pitch gain                                            */
    float *gain_code,    /* o    : quantized codebook gain                                         */
    float *gain_code2,    /* o    : quantized codebook gain                                        */
    ACELP_CbkCorr *coeff,        /* i/o  : correlations <y1,y1>, -2<xn,y1>,<y2,y2>, -2<xn,y2> and 2<y1,y2> */
    float *past_gcode,   /* (i/o): past gain of code                                               */
    float *gain_inov    /* (o)  : unscaled innovation gain                                        */
    ,short noisy_speech_flag /* (i) : noisy speech flag                                            */
)
{
    short index;
    float gcode,gcode2,pred_nrg_frame;
    float norm_code2;
    short i;
    float c, c_index2,c_first;
    short index2;

    /*-----------------------------------------------------------------*
     * - calculate the unscaled innovation energy
     *-----------------------------------------------------------------*/

    *gain_inov = 1.0f/ (float)sqrt( ( dotp( code, code, lcode ) + 0.01f ) / lcode );
    pred_nrg_frame = (float)pow(10.0,mean_ener/20.0);
    gcode=pred_nrg_frame*(*gain_inov);
    norm_code2 = 1.0f/ (float)sqrt( ( dotp( code2, code2, lcode ) + 0.01f ) / lcode );
    gcode2=pred_nrg_frame*(norm_code2);

    /*-----------------------------------------------------------------*
     * search for the best quantizer
     *-----------------------------------------------------------------*/

    *gain_code = coeff->xy2/(coeff->y2y2*gcode);

    if(*gain_code>0.000001f)
    {
        index=(int)(((20.f*log10(*gain_code)+20.f)/1.25f)+0.5f);

        if(index>31)
        {
            index=31;
        }
        else if(index<0)
        {
            index=0;
        }
    }
    else
    {
        index=0;
    }

    *gain_code= (float) pow(10.f,(((index*1.25f)-20.f)/20.f));
    *gain_code *= gcode;

    if( noisy_speech_flag )
    {
        c_first=0.8f*coeff->xx-(*gain_code)*(*gain_code)*coeff->y2y2;
    }
    else
    {
        c_first=coeff->xx-(*gain_code)*(*gain_code)*coeff->y2y2;
    }
    index2=0;
    *gain_code2 = (float) (index2*0.25f+0.25f)*(*gain_code*(gcode2/gcode));

    c_index2=c_first-(*gain_code2)*(*gain_code2)*coeff->y1y1-2*(*gain_code)*(*gain_code2)*coeff->y1y2;

    for(i=1; i<4; i++)
    {
        *gain_code2 = (float) (i*0.25f+0.25f)*(*gain_code*(gcode2/gcode));

        c=c_first-(*gain_code2)*(*gain_code2)*coeff->y1y1-2*(*gain_code)*(*gain_code2)*coeff->y1y2;

        if(fabs(c)<fabs(c_index2))
        {
            index2=i;
            c_index2=c;
        }
    }

    *gain_code2 = (float) (index2*0.25f+0.25f)*(*gain_code*(gcode2/gcode));
    index+=index2*32;

    *gain_pit  = 0.f;
    *past_gcode=*gain_code/ *gain_inov;  /*unscaled gain*/

    return index;
}
