/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"

/*-------------------------------------------------------------------*
 * Local constants
 *-------------------------------------------------------------------*/

#define L_ENR            (NB_SSF+2)
#define VOI_THRLD        0.2f

/*-------------------------------------------------------------------*
 * Local functions
 *-------------------------------------------------------------------*/

static float find_ener_decrease( const short ind_deltaMax, const float *pt_enr_ssf );

/*-------------------------------------------------------------------*
 * find_ener_decrease()
 *
 * Find maximum energy ratio between short sub-subframes in case
 * energy is trailing off after a spike
 *-------------------------------------------------------------------*/

static float find_ener_decrease(        /* o  : maximum energy ratio                              */
    const short ind_deltaMax,           /* i  : index of the beginning of maximum energy search   */
    const float *pt_enr_ssf             /* i  : Pointer to the energy buffer                      */
)
{
    short i, j, end, flag;
    float maxEnr, minEnr, dE2;

    dE2 = 0.0f;
    j = ind_deltaMax+2;
    end = j+L_ENR;
    maxEnr = pt_enr_ssf[j];
    j++;
    flag = 0;
    for( i=j; i<end; i++ )
    {
        if( pt_enr_ssf[i] > maxEnr && flag == 0 )
        {
            maxEnr = pt_enr_ssf[i];
            j++;
        }
        else
        {
            flag = 1;
        }
    }

    minEnr = maxEnr;
    for( i=j; i<end; i++ )
    {
        if( pt_enr_ssf[i] < minEnr )
        {
            minEnr = pt_enr_ssf[i];
        }
    }

    dE2 = maxEnr / (minEnr + 1.0e5f);

    return dE2;
}


/*-------------------------------------------------------------------*
 * find_uv()
 *
 * Decision about coder type
 *-------------------------------------------------------------------*/

short find_uv(                              /* o  : coding type                                   */
    Encoder_State *st,                        /* i/o: encoder state structure                       */
    const float *pitch_fr,                  /* i  : pointer to adjusted fractional pitch (4 val.) */
    const float *voicing_fr,                /* i  : refined correlation for each subframes        */
    const float *voicing,                   /* i  : correlation for 3 half-frames                 */
    const float *speech,                    /* i  : pointer to speech signal for E computation    */
    const short localVAD,                   /* i  : vad without hangover                          */
    const float *ee,                        /* i  : lf/hf Energy ratio for present frame          */
    const float corr_shift,                 /* i  : normalized correlation correction in noise    */
    const float relE,                       /* i  : relative frame energy                         */
    const float Etot,                       /* i  : total energy                                  */
    const float hp_E[],                     /* i  : energy in HF                                  */
    short *flag_spitch,                        /* i/o: flag to indicate very short stable pitch and high correlation */
    float voicing_sm                  /* i/o: smoothed open-loop pitch gains                */
    , const short last_core_orig              /* i  : original last core                            */
)
{
    const float *pt_speech;
    short i, coder_type, ind_deltaMax, tmp_offset_flag, nb_cond, flag_low_relE;
    float fac, mean_voi3, mean_ee, relE_thres;
    float enr_ssf[4*NB_SSF+2];
    float dE1, *pt_enr_ssf, *pt_enr_ssf1, dE2, dE3, dE2_th, ee0_th, ee1_th, voi_th;
    float mean_voi3_offset;
    float voicing_m, dpit1, dpit2, dpit3;

    /*-----------------------------------------------------------------*
     * Detect sudden energy increases to catch voice and music attacks (dE1)
     *
     * - Find maximum energy per short sub-subframe
     *   two sub-subframe sets are used, shifted by half the sub-subframe length
     * - Find maximum energy increase (ratio) between adjacent sub-subframes
     *-----------------------------------------------------------------*/

    /* find maximum energy per sub-subframe  */
    pt_speech = speech - SSF;
    pt_enr_ssf = enr_ssf + 2*NB_SSF;
    for( i=0; i < 2*(NB_SSF+1); i++ )
    {
        emaximum( pt_speech, SSF, pt_enr_ssf );
        pt_speech += (SSF/2);
        pt_enr_ssf++;
    }

    dE1 = 0.0f;
    ind_deltaMax = 0;
    pt_enr_ssf = enr_ssf + 2*NB_SSF;
    pt_enr_ssf1 = pt_enr_ssf + 2;

    /* test on energy increase between adjacent sub-subframes */
    for( i=0 ; i < 2*NB_SSF ; i++ )
    {
        fac = *pt_enr_ssf1 / (*pt_enr_ssf + 1.0f);
        if( fac > dE1 )
        {
            dE1 = fac;
            ind_deltaMax = i;
        }

        pt_enr_ssf++;
        pt_enr_ssf1++;
    }

    /*-----------------------------------------------------------------*
     * Average spectral tilt
     * Average voicing (normalized correlation)
     *-----------------------------------------------------------------*/

    mean_ee = 1.0f/3.0f * (st->ee_old + ee[0] + ee[1]);
    mean_voi3 = 1.0f/3.0f * (voicing[0] + voicing[1] + voicing[2]);

    /*-----------------------------------------------------------------*
     * Total frame energy difference (dE3)
     *-----------------------------------------------------------------*/

    dE3 = Etot - st->Etot_last;

    /*-----------------------------------------------------------------*
     * Energy decrease after spike (dE2)
     *-----------------------------------------------------------------*/

    /* set different thresholds and conditions for NB and WB input */
    if ( st->input_bwidth == NB )
    {
        dE2_th = 21.0f;
        nb_cond = ( mean_voi3 + corr_shift ) < 0.68f;
    }
    else
    {
        dE2_th = 30.0f;
        nb_cond = 1;        /* no additional condition for WB input */
    }

    /* calcualte maximum energy decrease */
    dE2 = 0.0f;
    pt_enr_ssf = enr_ssf + 2*NB_SSF;

    if( dE1 > 30.0f && nb_cond )
    {
        if( 2 * NB_SSF - ind_deltaMax < L_ENR )
        {
            st->old_ind_deltaMax = ind_deltaMax;
            mvr2r( pt_enr_ssf, st->old_enr_ssf, 2*NB_SSF );
        }
        else
        {
            st->old_ind_deltaMax = -1;
            dE2 = find_ener_decrease( ind_deltaMax, pt_enr_ssf );
            if( dE2 > dE2_th )
            {
                st->spike_hyst = 0;
            }
        }
    }
    else
    {
        if( st->old_ind_deltaMax >= 0 )
        {
            mvr2r( st->old_enr_ssf, enr_ssf, 2*NB_SSF );
            dE2 = find_ener_decrease( st->old_ind_deltaMax, enr_ssf );
            if( dE2 > dE2_th )
            {
                st->spike_hyst = 1;
            }
        }
        st->old_ind_deltaMax = -1;
    }

    /*-----------------------------------------------------------------*
     * Detection of voiced offsets (tmp_offset_flag)
     *-----------------------------------------------------------------*/

    tmp_offset_flag = 1;

    if ( st->input_bwidth != NB )
    {
        ee0_th = 2.4f;
        voi_th = 0.74f;
    }
    else
    {
        ee0_th = 9.8f;
        voi_th = 0.76f;
    }

    if( ( st->last_coder_type_raw == UNVOICED ) ||             /* previous frame was unvoiced  */
            ( ( ee[0] < ee0_th ) && ( hp_E[0] > (float)E_MIN ) &&  /* energy is concentrated in high frequencies provided that some energy is present in HF. The cast to (float) is needed for Linux i686 (gcc version 4.7.2), otherwise the criterion hp_E[0] > E_MIN holds true if E_MIN was assigned to hp_E[0] before */
              ( voicing[0] + corr_shift < voi_th ) ) )             /* normalized correlation is low */
    {
        tmp_offset_flag = 0;
    }

    /*-----------------------------------------------------------------*
     * Decision about UC
     *-----------------------------------------------------------------*/

    /* SC-VBR - set additional parameters and thresholds for SC-VBR */
    mean_voi3_offset = 0.0f;
    flag_low_relE = 0;
    ee1_th = 9.5f;
    if ( st->Opt_SC_VBR )
    {
        ee1_th = 8.5f;

        /* SC-VBR - determine the threshold on relative energy as a function of lp_noise */
        if ( st->input_bwidth != NB )
        {
            if ( st->Last_Resort == 0 )
            {
                relE_thres = 0.650f * st->lp_noise - 33.5f;
            }
            else
            {
                relE_thres = 0.700f * st->lp_noise - 33.5f;
            }
        }
        else
        {
            relE_thres = 0.60f * st->lp_noise - 28.2f;

        }

        if( relE_thres < -25.0f )
        {
            relE_thres = -25.0f;
        }

        /* SC-VBR = set flag on low relative energy */
        if ( relE < relE_thres )
        {
            flag_low_relE = 1;
        }

        /* SC-VBR - correction of voicing threshold for NB inputs (important only in noisy conditions) */
        if ( st->input_bwidth == NB && st->vadnoise < 20.0f )
        {
            mean_voi3_offset = 0.05f;
        }
    }

    /* make decision whether frame is unvoiced */
    coder_type = GENERIC;
    if ( st->input_bwidth == NB )
    {
        if( ( ( mean_voi3 + corr_shift < 0.68f + mean_voi3_offset ) &&         /* normalized correlation low  */
                ( ( voicing[2] + corr_shift ) < 0.79f ) &&      /* normalized correlation low on look-ahead - onset detection */
                ( ee[0] < 10.0f ) && ( hp_E[0] > (float)E_MIN ) &&     /* energy concentrated in high frequencies provided that some energy is present in HF...  */
                ( ee[1] < ee1_th ) && ( hp_E[1] > (float)E_MIN ) &&    /* ... biased towards look-ahead to detect onsets. The cast to (float) is needed for Linux i686 (gcc version 4.7.2), otherwise the criterion hp_E[0] > E_MIN holds true if E_MIN was assigned to hp_E[] before */
                ( tmp_offset_flag == 0 ) &&                     /* take care of voiced offsets */
                ( st->music_hysteresis == 0 ) &&                /*  ... and in segment after AUDIO frames   */
                ( dE1 <= 29.0f ) &&                             /* avoid on sharp energy spikes  */
                ( st->old_dE1 <= 29.0f ) &&                     /*   + one frame hysteresis   */
                ( st->spike_hyst < 0 ) ) ||                     /* avoid after sharp energy spikes followed by decay (e.g. castanets) */
                flag_low_relE )                                 /* low relative frame energy (only for SC-VBR) */
        {
            coder_type = UNVOICED;
        }
    }
    else
    {
        if( ( ( mean_voi3 + corr_shift < 0.695f + mean_voi3_offset ) &&        /* normalized correlation low  */
                ( ee[0] < 6.2f ) && ( hp_E[0] > (float)E_MIN ) &&      /* energy concentrated in high frequencies provided that some energy is present in HF */
                ( ee[1] < 6.2f ) && ( hp_E[1] > (float)E_MIN ) &&      /* ... biased towards look-ahead to detect onsets. The cast to (float) is needed for Linux i686 (gcc version 4.7.2), otherwise the criterion hp_E[0] > E_MIN holds true if E_MIN was assigned to hp_E[] before */
                ( tmp_offset_flag == 0 ) &&                     /* take care of voiced offsets */
                ( st->music_hysteresis == 0 ) &&                /*  ... and in segment after AUDIO frames   */
                ( dE1 <= 30.0f ) &&                             /* avoid on sharp energy spikes  */
                ( st->old_dE1 <= 30.0f ) &&                     /*   + one frame hysteresis   */
                ( st->spike_hyst < 0 ) ) ||                     /* avoid after sharp energy spikes followed by decay (e.g. castanets) */
                ( flag_low_relE
                  &&  st->old_dE1 <= 30.0f
                ))                                 /* low relative frame energy (only for SC-VBR) */
        {
            coder_type = UNVOICED;
        }
    }

    /*-----------------------------------------------------------------*
     * Decision about VC
     *-----------------------------------------------------------------*/

    st->set_ppp_generic = 0;

    if( localVAD == 1 && coder_type == GENERIC && last_core_orig != AMR_WB_CORE)
    {
        if(  ( voicing_fr[0] > 0.605f ) &&                    /* normalized correlation high in 1st sf.  */
                ( voicing_fr[1] > 0.605f ) &&                    /* normalized correlation high in 2st sf.  */
                ( voicing_fr[2] > 0.605f ) &&                    /* normalized correlation high in 3st sf.  */
                ( voicing_fr[3] > 0.605f ) &&                    /* normalized correlation high in 4st sf.  */
                ( mean_ee > 4.0f ) &&                            /* energy concentrated in low frequencies  */
                ( fabs( pitch_fr[1] - pitch_fr[0] ) < 3.0f ) &&    /* small OL pitch difference in 1st sf.    */
                ( fabs( pitch_fr[2] - pitch_fr[1] ) < 3.0f ) &&    /* small OL pitch difference in 2nd sf.    */
                ( fabs( pitch_fr[3] - pitch_fr[2] ) < 3.0f ) )     /* small OL pitch difference in 3rd sf.    */
        {
            coder_type = VOICED;
        }
        else if ( st->Opt_SC_VBR && st->input_bwidth == NB && st->vadnoise < 20 )
        {
            if( ( voicing_fr[0] > 0.25f) &&                   /* normalized correlation high in 1st sf. */
                    ( voicing_fr[1] > 0.25f ) &&                  /* normalized correlation high in 2st sf. */
                    ( voicing_fr[2] > 0.25f ) &&                  /* normalized correlation high in 3st sf. */
                    ( voicing_fr[3] > 0.25f ) &&                  /* normalized correlation high in 4st sf. */
                    ( mean_ee > 1.0f ) &&                         /* energy concentrated in low frequencies (used 1.0 for WB) */
                    ( fabs( pitch_fr[1] - pitch_fr[0] ) < 5.0f ) && /* small OL pitch difference in 1st sf.    */
                    ( fabs( pitch_fr[2] - pitch_fr[1] ) < 5.0f ) && /* small OL pitch difference in 2nd sf.    */
                    ( fabs( pitch_fr[3] - pitch_fr[2] ) < 5.0f ) )  /* small OL pitch difference in 3rd sf.    */
            {
                st->set_ppp_generic = 1;
                coder_type = VOICED;
            }
        }
        /* set VOICED mode for frames with very stable pitch and high correlation
           and avoid to switch to AUDIO/MUSIC later                              */
        voicing_m = mean( voicing_fr, NB_SUBFR );

        dpit1 = (float)fabs( pitch_fr[0] - pitch_fr[1] );
        dpit2 = (float)fabs( pitch_fr[1] - pitch_fr[2] );
        dpit3 = (float)fabs( pitch_fr[2] - pitch_fr[3] );

        if( *flag_spitch || ( dpit1 <= 3.0f && dpit2 <= 3.0f && dpit3 <= 3.0f &&
                              voicing_m > 0.95f && voicing_sm > 0.97f ) )
        {
            coder_type = VOICED;
            *flag_spitch = 1;  /*to avoid switch to AUDIO/MUSIC later*/
        }
    }

    /*-----------------------------------------------------------------*
     * Channel-aware mode - set RF mode and total bitrate
     *-----------------------------------------------------------------*/

    st->rf_mode = st->Opt_RF_ON;

    if ( coder_type == GENERIC )
    {
        if( ( voicing_fr[0] < VOI_THRLD ) &&          /* normalized correlation high in 1st sf.  */
                ( voicing_fr[1] < VOI_THRLD ) &&          /* normalized correlation high in 2st sf.  */
                ( voicing_fr[2] < VOI_THRLD ) &&          /* normalized correlation high in 3st sf.  */
                ( voicing_fr[3] < VOI_THRLD ) &&          /* normalized correlation high in 4st sf.  */
                ( st->vadnoise > 25.0f ) )                /* when speech is clean */
        {
            st->rf_mode = 0;
            /* Current frame cannot be compressed to pack the partial redundancy;*/

            if( st->rf_mode != st->Opt_RF_ON )
            {
                core_coder_mode_switch( st, st->bwidth, st->total_brate );
            }
        }
    }

    /*-----------------------------------------------------------------*
     * Updates
     *-----------------------------------------------------------------*/

    /* update spike hysteresis parameters */
    if( st->spike_hyst >= 0 && st->spike_hyst < 2 )
    {
        st->spike_hyst++;
    }

    /* reset spike hysteresis */
    if( ( st->spike_hyst > 1 ) &&
            ( dE3 > 5.0f ||                                                 /* energy increases               */
              ( relE > -13.0f && ( mean_voi3 + corr_shift > 0.695f) ) ) )   /* normalized correlation is high */
    {
        st->spike_hyst = -1;
    }

    /* update tilt parameters */
    st->ee_old = ee[1];
    st->old_dE1 = dE1;

    /* save the raw coder_type for various modules later in the codec (the reason is that e.g. UNVOICED is not used (rewritten) at higher rates) */
    st->coder_type_raw = coder_type;

    return coder_type;
}
