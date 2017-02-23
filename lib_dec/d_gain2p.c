/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include "options.h"
#include "typedef.h"
#include "prot.h"
#include "cnst.h"
#include "rom_com.h"


/*--------------------------------------------------------------------------*
* Mode2_gain_dec_mless
*
* Decoding of pitch and codebook gains without updating long term energies
*-------------------------------------------------------------------------*/

static void Mode2_gain_dec_mless(
    int index,          /* (i)  : index of quantizer                      */
    float code[],       /* (i)  : Innovative code vector                  */
    int lcode,          /* (i)  : Subframe size                           */
    float *gain_pit,    /* (o)  : Quantized pitch gain                    */
    float *gain_code,   /* (o)  : Quantized codebook gain                 */
    float mean_ener,    /* (i)  : mean_ener defined in open-loop (2 bits) */
    float *past_gpit,   /* (i/o): past gain of pitch                      */
    float *past_gcode,  /* (i/o): past energy of code                     */
    float *gain_inov,   /* (o)  : un-scaled innovation gain               */
    short coder_type    /* (i)  : coder type for number of bits           */
)
{
    short i;
    float ener_code,gcode0;
    const Word16 *t_qua_gain;

    ener_code = 0.0f;

    if( coder_type == 0 )
    {
        *gain_inov = 1.0f / (float)sqrt( ( dotp( code, code, lcode ) + 0.01f ) / lcode);
    }
    else
    {
        ener_code = 0.01f;

        for(i=0; i<lcode; i++)
        {
            ener_code += code[i] * code[i];
        }
        *gain_inov = (float)sqrt((float)lcode / ener_code);
    }

    /*-----------------------------------------------------------------*
     * Select the gains quantization table
     *-----------------------------------------------------------------*/

    if( coder_type == 0 )
    {
        t_qua_gain = E_ROM_qua_gain5b_const;
    }
    else if(coder_type == 1)
    {
        t_qua_gain = E_ROM_qua_gain6b_const;
    }
    else
    {
        t_qua_gain = E_ROM_qua_gain7b_const;
    }

    /*-----------------------------------------------------------------*
     * decode pitch gain
     *-----------------------------------------------------------------*/

    *gain_pit = (float)(t_qua_gain[index*2])/(1<<14);

    /*-----------------------------------------------------------------*
     * calculate the predicted gain code
     *-----------------------------------------------------------------*/

    if( coder_type == 0 )
    {
        ener_code = 10 * (float)log10( ( dotp( code, code, lcode ) + 0.01f ) / lcode );
        gcode0 = (float) pow(10, 0.05*(mean_ener - ener_code));
    }
    else
    {
        ener_code = (float)(-10.0 * log10((float)lcode / ener_code));
        gcode0 = mean_ener - ener_code;
        gcode0 = (float)pow(10.0,gcode0/20.0);   /* predicted gain */
    }

    /*-----------------------------------------------------------------*
     * decode normalized codebook gain
     *-----------------------------------------------------------------*/

    *gain_code = (float)(t_qua_gain[index*2+1])/(1<<11) * gcode0;

    *past_gpit = *gain_pit;
    *past_gcode = *gain_code / *gain_inov;

    return;
}


/*---------------------------------------------------------------------*
 * gain_dec_uv
 *
 * Decoding of pitch and codebook gains for Unvoiced mode
 *---------------------------------------------------------------------*/

static void gain_dec_uv(
    int index,            /* i/o: Quantization index vector             */
    float *code,          /* i  : algebraic code excitation             */
    int lcode,            /* i  : Subframe size                         */
    float *gain_pit,      /* o  : Quantized pitch gain                  */
    float *gain_code,     /* o  : Quantized codebook gain               */
    float *past_gpit,     /* i/o: past gain of pitch                    */
    float *past_gcode,    /* i/o: past energy of code                   */
    float *gain_inov      /* o  : unscaled innovation gain              */
)
{
    /*-----------------------------------------------------------------*
     * Innovation energy (without gain)
     *-----------------------------------------------------------------*/

    *gain_inov = 1.0f / (float)sqrt( ( dotp( code, code, lcode ) + 0.01f ) / lcode );

    /*-----------------------------------------------------------------*
    * Decode pitch gain
    *-----------------------------------------------------------------*/
    *gain_pit = 0.0f;

    /*-----------------------------------------------------------------*
     * Decode codebook gain
     *-----------------------------------------------------------------*/

    *gain_code= (float)pow(10.f,(((index*1.9f)-30.f)/20.f));

    /*-----------------------------------------------------------------*
     * past gains for error concealment
     *-----------------------------------------------------------------*/

    *past_gpit = *gain_pit;
    *past_gcode = *gain_code;
    *gain_code *= *gain_inov;

    return;
}


/*---------------------------------------------------------------------*
 * gain_dec_gacelp_uv
 *
 * Decoding of pitch and codebook gains for Unvoiced mode
 *---------------------------------------------------------------------*/
static void gain_dec_gacelp_uv(
    int index,              /* i/o: Quantization index vector             */
    float *code,            /* i  : algebraic code excitation             */
    float *code2,           /* i  : algebraic code excitation             */
    float mean_ener,        /* i : mean energy                            */
    int lcode,              /* i  : Subframe size                         */
    float *gain_pit,        /* o  : Quantized pitch gain                  */
    float *gain_code,       /* o  : Quantized codebook gain               */
    float *gain_code2,      /* o  : Quantized codebook gain               */
    float *past_gpit,       /* i/o: past gain of pitch                    */
    float *past_gcode,      /* i/o: past energy of code                   */
    float *gain_inov        /* o  : unscaled innovation gain              */
)
{
    float pred_nrg_frame,norm_code2;
    float gcode, gcode2;
    short index2;

    pred_nrg_frame = (float)pow(10.0,mean_ener/20.0);

    /*-----------------------------------------------------------------*
     * Prediction gains
     *-----------------------------------------------------------------*/

    *gain_inov = 1.0f / (float)sqrt( ( dotp( code, code, lcode ) + 0.01f ) / lcode );
    gcode=pred_nrg_frame*(*gain_inov);

    norm_code2 = 1.0f/ (float)sqrt( ( dotp( code2, code2, lcode ) + 0.01f ) / lcode );
    gcode2=pred_nrg_frame*(norm_code2);

    /*-----------------------------------------------------------------*
     * Decode pitch gain
     *-----------------------------------------------------------------*/

    *gain_pit = 0.0f;
    *past_gpit = *gain_pit;

    /*-----------------------------------------------------------------*
     * past gains for error concealment
     *-----------------------------------------------------------------*/

    index2=index>>5;
    index=index&0x1F;

    *gain_code= (float)pow(10.f,(((index*1.25f)-20.f)/20.f))*gcode;
    *gain_code2 = (float) (index2*0.25f+0.25f)*(*gain_code*(gcode2/gcode));

    *past_gcode=*gain_code/ *gain_inov;  /*unscaled gain*/

    return;
}


/*---------------------------------------------------------------------*
 * decode_acelp_gains
 *
 *
 *---------------------------------------------------------------------*/

void decode_acelp_gains(
    float *code,
    int gains_mode,
    float mean_ener_code,
    float *gain_pit,
    float *gain_code,
    int **pt_indice,
    float *past_gpit,
    float *past_gcode,
    float *gain_inov,
    int L_subfr,
    float *code2,
    float *gain_code2
)
{
    int index = 0;

    index = **pt_indice;
    (*pt_indice)++;

    if (((gains_mode > 0) && (gains_mode < 4)))
    {
        /* EVS gains quantizer (5bits/subfr) */
        Mode2_gain_dec_mless(index, code, L_subfr, gain_pit, gain_code, mean_ener_code, past_gpit, past_gcode, gain_inov, gains_mode-1 );
    }
    else if(gains_mode == 6)
    {
        /* UV gains quantizer (6bits/subfr) */
        gain_dec_uv(index, code, L_subfr, gain_pit, gain_code, past_gpit, past_gcode, gain_inov );
    }
    else if(gains_mode == 7)
    {
        /* GACELP_UV gains quantizer (7=5-2bits/subfr) */
        gain_dec_gacelp_uv(index, code, code2, mean_ener_code, L_subfr, gain_pit, gain_code, gain_code2, past_gpit, past_gcode, gain_inov );
    }
    else
    {
        fprintf(stderr, "invalid gains coding for acelp!\n");
        assert(0);
    }

    return;
}
