/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"


/*-------------------------------------------------------------------*
 * Local constants
 *-------------------------------------------------------------------*/

#define BAND3k          15
#define BIN_1k2         48
#define BAND_2k         12
#define BAND_0k4        4
#define ALP             0.7f
#define MALP            (1.0f-ALP)
#define ALPMY           (0.86f)
#define ALPY            1.5f

/*-------------------------------------------------------------------*
 * Local functions
 *-------------------------------------------------------------------*/

static void NoiseFill( float *exc_diffQ, short *seed_tcx, const short Mbands_gn );
static void Ener_per_band( const float exc_diff[], float y_gain4[] );
static void Apply_gain( float exc_diffQ[], float Ener_per_bd_iQ[], float Ener_per_bd_yQ[]);
static void normalize_spec( float fac_up, float fy_norm[], const short L_frame );
static void gs_dec_amr_wb( const long core_brate, short *seed_tcx, const float dct_in[],
                           float dct_out[], const float pitch[], const float voice_fac
                           ,const short clas, const short coder_type
                         );


/*-------------------------------------------------------------------*
 * NoiseFill()
 *
 * noise fill function for unvoiced/inactive frames (used only in AMR-WB IO mode)
 *-------------------------------------------------------------------*/

static void NoiseFill(
    float *exc_diffQ,             /* i/o: Noise per band                          */
    short *seed_tcx,              /* i  : Seed for noise                          */
    const short Mbands_gn               /* i  : number of bands                         */
)
{
    short StartBin, NB_Qbins, i_band, k;
    StartBin = 0;
    NB_Qbins  = 0;

    for( i_band = 0; i_band < Mbands_gn; i_band++ )
    {
        StartBin += NB_Qbins;
        NB_Qbins = crit_bins[i_band];
        for( k=StartBin; k<NB_Qbins + StartBin; k++ )
        {
            exc_diffQ[k] += 0.75f*((float)own_random(seed_tcx)/32768.0f);
        }
    }

    return;
}

/*-------------------------------------------------------------------*
 * Ener_per_band()
 *
 * Computed the energy per band (used only in AMR-WB IO mode)
 *-------------------------------------------------------------------*/

static void Ener_per_band(
    const float exc_diff[],    /* i  : target signal                     */
    float y_gain4[]      /* o  : Energy per band to quantize       */
)
{
    float etmp;
    const float *pt;
    short i,j;

    pt  = exc_diff;
    for(j = 0; j < CRIT_NOIS_BAND; j++)
    {
        etmp = 0.01f;
        for(i = 0; i < crit_bins[j]; i++)
        {
            etmp += (*pt **pt);
            pt++;
        }
        etmp = max( etmp, .01f );
        y_gain4[j] = (float)sqrt(etmp);
    }

    return;
}


/*-------------------------------------------------------------------*
 * Apply_gain()
 *
 * Rescaling of the modified excitation vector (used only in AMR-WB IO mode)
 *-------------------------------------------------------------------*/

static void Apply_gain(
    float exc_diffQ[],       /* i/o: Quantized excitation              */
    float Ener_per_bd_iQ[],  /* o  : Target ener per band              */
    float Ener_per_bd_yQ[]   /* o  : Ener per band for norm vector     */
)
{
    short i, i_band;
    short StartBin, NB_Qbins;
    float y_gain;

    /*------------------------------------------------------------------
     * For all the bands
     * Find the energy ratio between modified vector and original vector
     *------------------------------------------------------------------*/

    StartBin = 0;
    NB_Qbins  = 0;
    for( i_band = 0; i_band < CRIT_NOIS_BAND; i_band++ )
    {
        StartBin += NB_Qbins;
        NB_Qbins = crit_bins[i_band];
        y_gain = (float)(Ener_per_bd_iQ[i_band]/Ener_per_bd_yQ[i_band]);

        /*------------------------------------------------------------------
         * For bands below 400 Hz or for unvoiced/inactive frames
         *   only apply the energy ratio
         *------------------------------------------------------------------*/

        for( i = StartBin; i < NB_Qbins + StartBin; i++ )
        {
            exc_diffQ[i] *= y_gain;
        }
    }

    return;
}

/*-------------------------------------------------------------------*
 * normalize_spec()
 *
 * Spectrum normalization (zeroed of bins below a certain threshold) (used only in AMR-WB IO mode)
 *-------------------------------------------------------------------*/

static void normalize_spec(
    float fac_up,           /* i  : Core bitrate                  */
    float fy_norm[],        /* i/o: Frequency quantized parameter */
    const short L_frame           /* i  : Section length                */
)
{
    float max_val;
    short idx, j;

    idx = emaximum(fy_norm, L_frame, &max_val);
    max_val = (float)fabs(fac_up/fy_norm[idx]);

    for(j = 0; j < L_frame; j++)
    {
        fy_norm[j] *= max_val;
    }

    return;
}

/*-------------------------------------------------------------------*
  * gs_dec_amr_wb()
  *
  * Modification of decoded excitation vector depending of the content type (used only in AMR-WB IO mode)
  *-------------------------------------------------------------------*/

static void gs_dec_amr_wb(
    const long core_brate,        /* i  : bitrate allocated to the core       */
    short *seed_tcx,        /* i/o: seed used for noise generation      */
    const float dct_in[],         /* i  : cdt of residual signal              */
    float dct_out[],        /* i/o: dct of pitch only excitation        */
    const float pitch[],          /* i  : pitch buffer                        */
    const float voice_fac         /* i  : gain pitch                          */
    ,const short clas,            /* i  : signal frame class                  */
    const short coder_type       /* i  : coder type                          */
)
{
    short i, mDiff_len;
    float exc_diffQ[L_FRAME16k];
    short j;
    float etmp14;
    float ftmp, ftmp1;
    float Ener_per_bd_iQ[CRIT_NOIS_BAND];
    float Ener_per_bd_yQ[CRIT_NOIS_BAND];

    /*--------------------------------------------------------------------------------------*
     * compute the energy per band for the decoded excitation (in frequency domain)
     *--------------------------------------------------------------------------------------*/

    Ener_per_band( dct_in, Ener_per_bd_iQ );

    /*--------------------------------------------------------------------------------------*
     * adjust quantization noise for the low level to compensate for the poor 6 bit gainQ
     *--------------------------------------------------------------------------------------*/

    if( core_brate < ACELP_12k65)
    {
        ftmp = 0;
        for(i = 0; i < CRIT_NOIS_BAND; i++)
        {
            ftmp  = max(Ener_per_bd_iQ[i], ftmp);
        }

        if( (coder_type == INACTIVE || clas == VOICED_TRANSITION) &&  ftmp < 20.0f )
        {
            for(i = 0; i < CRIT_NOIS_BAND; i++)
            {
                Ener_per_bd_iQ[i] *= crit_bins_corr[i];
            }
        }
    }
    /*--------------------------------------------------------------------------------------*
     * Find the length of the temporal contribution, with a minimum contribution of 1.2kHz
     *--------------------------------------------------------------------------------------*/

    minimum( pitch, 4, &etmp14 );
    etmp14 = 12800.0f/etmp14;
    etmp14 *= 8.0f;

    if( core_brate >= ACELP_12k65 )
    {
        etmp14 *=2;
    }

    mDiff_len = (short)(etmp14+0.5f);
    etmp14 = 32768.0f;
    j = 0;
    for(i = 0; i < CRIT_NOIS_BAND; i++)
    {
        if( fabs(crit_bands_loc[i] - mDiff_len) < etmp14 )
        {
            etmp14 = (float)fabs(crit_bands_loc[i] - mDiff_len);
            j += crit_bins[i];
        }
    }

    mDiff_len = j;
    if( mDiff_len < BIN_1k2 )
    {
        mDiff_len = BIN_1k2;
    }
    mvr2r(dct_in, exc_diffQ, mDiff_len);
    set_f(exc_diffQ+mDiff_len, 0, L_FRAME-mDiff_len);

    /*--------------------------------------------------------------------------------------*
     * normalization of the spectrum and noise fill
     *--------------------------------------------------------------------------------------*/

    normalize_spec(4, exc_diffQ , mDiff_len);
    NoiseFill( exc_diffQ, seed_tcx, CRIT_NOIS_BAND);


    /*--------------------------------------------------------------------------------------*
     * Recompute energy per band of the modified excitation vector (in frequency domain)
     *--------------------------------------------------------------------------------------*/

    Ener_per_band( exc_diffQ, Ener_per_bd_yQ );

    /*--------------------------------------------------------------------------------------*
     * Compute tilt factor and amplify HF accordingly
     *--------------------------------------------------------------------------------------*/

    ftmp = (0.5f * (1.0f - voice_fac));  /* 1=unvoiced, 0=voiced */

    for(i = 240; i < L_FRAME; i++)
    {
        ftmp1 = ftmp*(0.067f * i - 15.0f);
        ftmp1 = max(1.0f, ftmp1);
        exc_diffQ[i] *=  ftmp1;
    }

    /*--------------------------------------------------------------------------------------*
     * Match the energy of the modified excitation vector to the decoded excitation
     *--------------------------------------------------------------------------------------*/

    Apply_gain( exc_diffQ, Ener_per_bd_iQ, Ener_per_bd_yQ);

    /*--------------------------------------------------------------------------------------*
     * Copy to the output vector
     *--------------------------------------------------------------------------------------*/

    mvr2r( exc_diffQ, dct_out, L_FRAME );

    return;
}

/*-------------------------------------------------------------------*
 * improv_amr_wb_gs()
 *
 * Modify the decoded excitation to increase quality of
 * unvoiced and audio signals (used only in AMR-WB IO mode)
 *-------------------------------------------------------------------*/

void improv_amr_wb_gs(
    const short clas,                             /* i  : signal frame class                  */
    const short coder_type,                       /* i  : coder type                          */
    const long core_brate,                        /* i  : bitrate allocated to the core       */
    short *seed_tcx,                        /* i/o: Seed used for noise generation      */
    float *old_Aq,                          /* i/o: old LPC filter coefficient          */
    float *mem_syn2,                        /* i/o: synthesis memory                    */
    const float lt_voice_fac,                     /* i/o: long term voice factor              */
    const short locattack,                        /* i  : Flag for a detected attack          */
    float *Aq,                              /* i/o: Decoded LP filter coefficient       */
    float *exc2,                            /* i/o: Decoded complete excitation         */
    float *mem_tmp,                         /* i/o: synthesis temporary memory          */
    float *syn,                             /* i/o: Decoded synthesis to be updated     */
    const float *pitch_buf,                       /* i  : Decoded pitch buffer                */
    const float Last_ener                         /* i  : Last energy                         */
    ,const short rate_switching_reset              /* i  : rate switching reset flag           */
    ,const short last_coder_type                    /* i  : Last coder_type */
)
{
    short i;
    float dct_exc_in[L_FRAME], dct_exc_out[L_FRAME];
    float Aq_orig[NB_SUBFR*(M+1)], enr_LP_old, enr_LP_new;

    /*------------------------------------------------------------*
     * Condition to enter the section on excitation modification
     *------------------------------------------------------------*/

    /* Enter the modification for all inactive frames and also for unvoiced frames if bit rate is below 8k85 */
    if( ( locattack == 0 && core_brate <= ACELP_12k65) &&
            ( (core_brate < ACELP_8k85 && clas != AUDIO_CLAS && (clas == UNVOICED_CLAS || clas == VOICED_TRANSITION)) || coder_type == INACTIVE ) )
    {
        /*------------------------------------------------------------*
         * two differents paths:
         *   unvoiced or inactive
         *   generic audio sound
         * LP filter smoothing for inactive parts
         *------------------------------------------------------------*/
        *seed_tcx = (short)((short)(pitch_buf[0]*64.0f)*(short)(pitch_buf[3]*64.0f));

        if( coder_type == INACTIVE && Last_ener > -3.0f && last_coder_type == UNVOICED && rate_switching_reset == 0 )
        {
            mvr2r( Aq, Aq_orig, NB_SUBFR*(M+1) );

            for(i =0; i < NB_SUBFR * (M+1); i++)
            {
                Aq[i] = ALP*old_Aq[i] + MALP*Aq[i];
            }

            /* check the smoothed LP filter stability */
            enr_LP_old = enr_1_Az( old_Aq, L_SUBFR );
            for(i =0; i < NB_SUBFR; i++)
            {
                enr_LP_new = enr_1_Az( Aq+i*(M+1), L_SUBFR );

                if( enr_LP_new > 128 * enr_LP_old )
                {
                    /* filter is unstable, do not modify the excitation */
                    mvr2r( Aq_orig, Aq, NB_SUBFR*(M+1) );
                    return;
                }

                enr_LP_old = enr_LP_new;
            }
        }

        /*------------------------------------------------------------*
         * Find frequency representation of the excitation
         * Do the excitation modification according to the content
         * Go back to time domain -> Overwrite excitation
         *------------------------------------------------------------*/

        edct( exc2, dct_exc_in, L_FRAME );
        set_f( exc2, 0, L_FRAME );
        set_f( dct_exc_out, 0, L_FRAME );

        gs_dec_amr_wb( core_brate, seed_tcx, dct_exc_in, dct_exc_out, pitch_buf, lt_voice_fac, clas, coder_type );

        edct( dct_exc_out, exc2, L_FRAME );

        /*------------------------------------------------------------*
         * Redo core synthesis at 12k8 Hz with the modified excitation
         *------------------------------------------------------------*/

        mvr2r( mem_tmp, mem_syn2, M );
        syn_12k8( L_FRAME, Aq, exc2, syn, mem_syn2, 1 );
    }

    return;
}
