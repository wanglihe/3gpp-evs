/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "options.h"
#include "rom_com.h"
#include "prot.h"

/*-------------------------------------------------------------------*
* Local function
*--------------------------------------------------------------------*/

static int BITS_ALLOC_adjust_acelp_fixed_cdk( int bits_frame, int *fixed_cdk_index, int nb_subfr );


/*-------------------------------------------------------------------*
* BITS_ALLOC_init_config_acelp()
*
* initial configuration for Mode 2 ACELP
*--------------------------------------------------------------------*/

void BITS_ALLOC_init_config_acelp(
    int bit_rate,
    int narrowBand,
    int nb_subfr,
    ACELP_config *acelp_cfg        /*o:  configuration structure of ACELP*/
)
{
    short rate_mode_index;

    if( bit_rate <= ACELP_9k60 )
    {
        rate_mode_index=0;
    }
    else
    {
        rate_mode_index=1;
    }

    acelp_cfg->mode_index=rate_mode_index;

    /*LPC: midLpc should be swithced off?*/
    acelp_cfg->midLpc_enable = 1;

    /*ACELP ICB config*/
    if( (rate_mode_index==0) || (narrowBand==1) )
    {
        acelp_cfg->pre_emphasis = 1;
        acelp_cfg->formant_enh = 1;
        acelp_cfg->formant_enh_num = FORMANT_SHARPENING_G1;
        acelp_cfg->formant_enh_den = FORMANT_SHARPENING_G2;
        acelp_cfg->formant_tilt = 0;
        acelp_cfg->voice_tilt = 0;
    }
    else
    {
        acelp_cfg->pre_emphasis = 0;
        acelp_cfg->formant_enh = 1;
        acelp_cfg->formant_enh_num = FORMANT_SHARPENING_G1;
        acelp_cfg->formant_enh_den = FORMANT_SHARPENING_G2;
        acelp_cfg->formant_tilt = 1;
        acelp_cfg->voice_tilt = 1;
    }

    /*Wide band @ 16kHz*/
    if ( nb_subfr == NB_SUBFR16k )
    {
        acelp_cfg->pre_emphasis = 1;
        acelp_cfg->formant_enh = 1;
        acelp_cfg->formant_enh_num = FORMANT_SHARPENING_G1_16k;
        acelp_cfg->formant_enh_den = FORMANT_SHARPENING_G2_16k;
        acelp_cfg->formant_tilt = 0;
        acelp_cfg->voice_tilt = 2;
    }

    return;
}

/*-------------------------------------------------------------------*
* BITS_ALLOC_config_acelp()
*
* configure all Mode 2 ACELP modes and allocate the bits
*--------------------------------------------------------------------*/

int BITS_ALLOC_config_acelp(
    const int   bits_frame,         /* i  : remaining bit budget for the frame  */
    const short coder_type,         /* i  : acelp extended mode index           */
    ACELP_config *acelp_cfg,         /* i/o: configuration structure of ACELP    */
    const short narrowBand,         /* i  : narrowband flag                     */
    const short nb_subfr            /* i  : number of subframes                 */
)
{
    short mode_index;
    short band_index;
    short i;
    short remaining_bits, bits;

    /*Sanity check*/

    mode_index = acelp_cfg->mode_index;
    band_index = (narrowBand==0);
    bits = 0;

    if ( band_index==0 )
    {
        if(coder_type == INACTIVE)
        {
            acelp_cfg->formant_enh = 0;
        }
        else
        {
            acelp_cfg->formant_enh = 1;
        }
    }

    if( band_index==1 && nb_subfr == NB_SUBFR )
    {

        if( coder_type == INACTIVE)
        {
            acelp_cfg->pre_emphasis = 0;
            acelp_cfg->formant_enh = 0;
            acelp_cfg->formant_enh_num = FORMANT_SHARPENING_G1_16k;
            acelp_cfg->formant_tilt = 1;
            acelp_cfg->voice_tilt = 1;
        }
        else
        {
            acelp_cfg->pre_emphasis = 1;
            acelp_cfg->formant_enh = 1;
            acelp_cfg->formant_enh_num = FORMANT_SHARPENING_G1;
            acelp_cfg->formant_tilt = 0;
            acelp_cfg->voice_tilt = 0;
        }
    }

    if( coder_type == UNVOICED )
    {
        if( ACELP_GAINS_MODE[mode_index][band_index][coder_type]==6 )
        {
            acelp_cfg->pitch_sharpening = 0;
            acelp_cfg->phase_scrambling = 1;
        }
        else
        {
            acelp_cfg->pitch_sharpening = 0;
            acelp_cfg->phase_scrambling = 0;
        }
    }
    else
    {
        acelp_cfg->pitch_sharpening = 1;
        acelp_cfg->phase_scrambling = 0;
    }

    if( coder_type > ACELP_MODE_MAX )
    {
        /* keep pitch sharpening for RF_ALLPRED mode */
        acelp_cfg->pitch_sharpening = 0;
        acelp_cfg->phase_scrambling = 0;
    }

    /*Allocate bits and different modes*/
    acelp_cfg->bpf_mode=ACELP_BPF_MODE[mode_index][band_index][coder_type];
    bits+=ACELP_BPF_BITS[acelp_cfg->bpf_mode];

    acelp_cfg->nrg_mode=ACELP_NRG_MODE[mode_index][band_index][coder_type];
    acelp_cfg->nrg_bits=ACELP_NRG_BITS[acelp_cfg->nrg_mode];
    bits+=acelp_cfg->nrg_bits;

    acelp_cfg->ltp_mode=ACELP_LTP_MODE[mode_index][band_index][coder_type];
    acelp_cfg->ltp_bits=0;
    acelp_cfg->ltf_mode=ACELP_LTF_MODE[mode_index][band_index][coder_type];
    acelp_cfg->ltf_bits=ACELP_LTF_BITS[acelp_cfg->ltf_mode];

    if( nb_subfr == NB_SUBFR16k && acelp_cfg->ltf_bits == 4 )
    {
        acelp_cfg->ltf_bits++;
    }
    bits+=acelp_cfg->ltf_bits;


    for ( i=0; i<nb_subfr; i++ )
    {
        acelp_cfg->gains_mode[i] = ACELP_GAINS_MODE[mode_index][band_index][coder_type];

        /* skip subframe 1, 3 gain encoding, and use from subframe 0, and 3, respectively */
        if(coder_type >= ACELP_MODE_MAX && (i == 1 || i == 3))
        {
            acelp_cfg->gains_mode[i] = 0;
        }

        bits += ACELP_GAINS_BITS[acelp_cfg->gains_mode[i]];
        bits += ACELP_LTP_BITS_SFR[acelp_cfg->ltp_mode][i];
        acelp_cfg->ltp_bits += ACELP_LTP_BITS_SFR[acelp_cfg->ltp_mode][i];
    }

    /*Innovation*/
    if( bits_frame < bits )
    {
        printf("\nWarning: bits per frame too low\n");
        return -1;
    }

    if( coder_type == RF_ALLPRED )
    {
        set_i(acelp_cfg->fixed_cdk_index, -1, nb_subfr);
    }
    else if ( coder_type == RF_GENPRED )
    {
        acelp_cfg->fixed_cdk_index[0] = 0;  /* 7 bits */
        acelp_cfg->fixed_cdk_index[1] = -1;
        acelp_cfg->fixed_cdk_index[2] = 0;  /* 7 bits */
        acelp_cfg->fixed_cdk_index[3] = -1;
        acelp_cfg->fixed_cdk_index[4] = -1;
        bits += 14;
    }
    else if( coder_type == RF_NOPRED )
    {
        set_i(acelp_cfg->fixed_cdk_index, 0, nb_subfr);
        bits += 28;
    }
    else
    {
        bits += BITS_ALLOC_adjust_acelp_fixed_cdk(bits_frame - bits, acelp_cfg->fixed_cdk_index, nb_subfr );
    }

    remaining_bits = bits_frame-bits;

    /*Sanity check*/
    if( remaining_bits < 0 )
    {
        bits = -1;
    }


    return( bits );
}

/*-------------------------------------------------------------------*
* BITS_ALLOC_adjust_acelp_fixed_cdk()
*
*
*--------------------------------------------------------------------*/

static int BITS_ALLOC_adjust_acelp_fixed_cdk(
    int bits_frame, /*i: bit budget*/
    int *fixed_cdk_index,
    int nb_subfr
)
{
    int bits_subframe2;
    int sfr, k, bitsused, bits_currsubframe;

    bits_subframe2 = bits_frame;

    if( bits_subframe2 < ACELP_FIXED_CDK_BITS(0)*nb_subfr )
    {
        return(bits_frame+1 ); /* Not enough bits for lowest mode. -> trigger alarm*/
    }

    /* search cdk-index for first subframe */
    for (k=0; k<ACELP_FIXED_CDK_NB-1; k++)
    {
        if (ACELP_FIXED_CDK_BITS(k)*nb_subfr > bits_subframe2)
        {
            k--;    /* previous mode did not exceed bit-budget */
            break;
        }
    }

    if( ACELP_FIXED_CDK_BITS(k)*nb_subfr > bits_subframe2 )
    {
        k--;    /* previous mode did not exceed bit-budget */
    }
    fixed_cdk_index[0] = k;
    bitsused = ACELP_FIXED_CDK_BITS(k);

    for (sfr=1; sfr < nb_subfr; sfr++)
    {
        bits_currsubframe = (sfr*bits_subframe2 + bits_subframe2) - bitsused*nb_subfr;

        /* try increasing mode while below threshold */
        while ( (k < ACELP_FIXED_CDK_NB-1) && (ACELP_FIXED_CDK_BITS(k+1)*nb_subfr <= bits_currsubframe) )
        {
            k++;
        }

        /* try decreasing mode until below threshold */
        while (ACELP_FIXED_CDK_BITS(k)*nb_subfr > bits_currsubframe)
        {
            k--;
            if ( k == 0 )
            {
                break;
            }
        }

        /* store mode */
        fixed_cdk_index[sfr] = k;
        bitsused += ACELP_FIXED_CDK_BITS(k);
    }

    return bitsused;
}
