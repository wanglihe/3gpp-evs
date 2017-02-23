/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <assert.h>
#include "options.h"
#include "prot.h"
#include "cnst.h"
#include "stat_com.h"


/* PLC: [Common: Fade-out]
 * PLC: and for PLC fade out */

void minimumStatistics(
    float*      noiseLevelMemory,
    int*        noiseLevelIndex,
    int*        currLevelIndex,
    float*      noiseEstimate,
    float*      lastFrameLevel,
    float       currentFrameLevel,
    float const minLev,
    int   const buffSize
)
{
    float aOpt;
    float f;
    int p;
    int i;

    if (currentFrameLevel < minLev)
    {
        currentFrameLevel = minLev;
    }
    /* compute optimal factor aOpt for recursive smoothing of frame minima */
    if (*lastFrameLevel >= *noiseEstimate)
    {
        aOpt = *noiseEstimate / *lastFrameLevel;
    }
    else
    {
        aOpt = *lastFrameLevel / *noiseEstimate;
    }
    aOpt *= aOpt;
    *lastFrameLevel = currentFrameLevel;
    /* recursively compute smoothed frame minima using optimal factor aOpt */
    f = currentFrameLevel * (1.0f - aOpt);
    f += aOpt * noiseLevelMemory[(*currLevelIndex ? *currLevelIndex : buffSize)-1];
    /* if current frame min is a new local min, set index to current index */
    p = *noiseLevelIndex;
    if (noiseLevelMemory[p] >= f)
    {
        noiseLevelMemory[*currLevelIndex] = f;
        p = *currLevelIndex;
    }
    else
    {
        noiseLevelMemory[*currLevelIndex] = f;
        /* current min is not a new min, so check if min must be re-searched */
        if (p != *currLevelIndex)
        {
            f = noiseLevelMemory[p];   /* min is still in memory, so return it */
        }
        else
        {
            /* p == currLevelIndex; min was removed from memory, re-search min */
            for (i = *currLevelIndex + 1; i < buffSize; i++)
            {
                if (f >= noiseLevelMemory[i])
                {
                    f = noiseLevelMemory[i];
                    p = i;
                }
            }
            for (i = 0; i <= *currLevelIndex; i++)
            {
                if (f >= noiseLevelMemory[i])
                {
                    f = noiseLevelMemory[i];
                    p = i;
                }
            }
        }
    }
    /* update local-minimum-value index and current circular-buffer index */
    *noiseLevelIndex = p;
    p = *currLevelIndex + 1;
    *currLevelIndex = (p == buffSize) ? 0 : p;

    *noiseEstimate = f;

    return;
}


/*----------------------------------------------------------------------*
 * PLC: [ACELP: Fade-out]
 * PLC: getLevelSynDeemph: derives on frame or subframe basis the level
 *      of LPC synthesis and deeemphasis based on the given input
 *----------------------------------------------------------------------*/
float getLevelSynDeemph(float const h1Init[],     /* i: input value or vector to be processed */
                        float const A[],          /* i: LPC coefficients                      */
                        int   const lenLpcExc,    /* i: length of the LPC excitation buffer   */
                        float const preemph_fac,  /* i: preemphasis factor                    */
                        int   const numLoops)     /* i: number of loops                       */
{
    float levelSynDeemphSub;
    float levelSynDeemph = 0;
    float h1[L_FRAME_PLUS/4];
    float mem[M];
    float tmp = 0;
    int   loop;

    for (loop = 0; loop  < numLoops; loop++)
    {
        set_zero(h1, lenLpcExc);
        set_zero(mem, M);

        h1[0] = *h1Init;

        syn_filt(A, M, h1, h1, lenLpcExc, mem, 0);
        deemph(h1, preemph_fac, lenLpcExc, &tmp);
        A += (M+1);

        /* gain introduced by synthesis+deemphasis */
        levelSynDeemphSub = (float)sqrt(dotp( h1, h1, lenLpcExc));

        /* mean of the above across all subframes */
        levelSynDeemph += (1.0f/(float)numLoops) * levelSynDeemphSub;
    }
    return levelSynDeemph;
}

void genPlcFiltBWAdap(int   const sr_core,     /* i: core sampling rate                                         */
                      float*      lpFiltAdapt, /* o: filter coefficients for filtering codebooks in case of flc */
                      int   const type,        /* i: type of filter, either 0 : lowpass or 1 : highpass         */
                      float const alpha        /* i: fade out factor [0 1) used decrease filter tilt            */
                     )
{
    float a;
    switch (sr_core)
    {
    case 16000 :
        a = 0.4000f;
        break;
    default    :
        a = 0.2813f; /*sr_core = 12800*/
        break;
    }
    switch (type)
    {
    case 0 :
        *lpFiltAdapt++ =   a/(2.f*a+1.f);
        *lpFiltAdapt++ = 1.f/(2.f*a+1.f);
        *lpFiltAdapt   =   a/(2.f*a+1.f);
        break;
    case 1 :
        a *= alpha;
        *lpFiltAdapt++ =  -a/(2.f*a+1.f);
        *lpFiltAdapt++ = 1.f/(2.f*a+1.f);
        *lpFiltAdapt   =  -a/(2.f*a+1.f);
        break;
    default   :
        fprintf(stderr,"PLC: Filter type neither lowpass nor highpass.\n");
        assert(0);
        break;
    }

}


/*-----------------------------------------------------------------*
 * PLC: [ACELP: general]
 * PLC: high pass filtering
 *-----------------------------------------------------------------*/
void highPassFiltering(const short last_good,     /* i:   last classification type                           */
                       const int   L_buffer,      /* i:   buffer length                                      */
                       float       exc2[],        /* i/o: unvoiced excitation before the high pass filtering */
                       const float hp_filt[],     /* i:   high pass filter coefficients                      */
                       const int   l_fir_fer)     /* i:   high pass filter length                            */
{
    int   i;

    if( last_good > UNVOICED_TRANSITION )
    {
        for( i=0 ; i< L_buffer; i++ )
        {
            exc2[i] = dotp(&exc2[i], hp_filt, l_fir_fer);
        }
    }
}


/*----------------------------------------------------------------------------------*
 * PLC: [Common: mode decision]
 * PLC: Decide which Concealment to use. Update pitch lags if needed
 *----------------------------------------------------------------------------------*/
int GetPLCModeDecision(Decoder_State *st   /* i/o:    decoder memory state pointer */
                      )
{
    int core;
    int numIndices = 0;
    if( st->flagGuidedAcelp == 1 )
    {
        /* update mem_lag according to info available on future frame */
        st->old_pitch_buf[2*st->nb_subfr] = (float)st->guidedT0;
        st->old_pitch_buf[2*st->nb_subfr+1] = (float)st->guidedT0;
        st->mem_pitch_gain[0] = st->mem_pitch_gain[1] = 1.f;
    }
    if(( st->last_core > ACELP_CORE && st->tcxltp_last_gain_unmodified != 0 ) || ( st->flagGuidedAcelp == 1 ))
    {
        /* no updates needed here, because already updated in last good frame */
        st->plc_use_future_lag = 1;
    }
    else
    {
        st->plc_use_future_lag = 0;
    }
    if (st->last_core == -1)
    {
        if (st->Opt_AMR_WB)
        {
            core = 0;
        }
        else
        {
            core = 1;
        }
        st->last_core = ACELP_CORE;
        st->tonal_mdct_plc_active = 0;
    }
    else
    {
        core = 0;
        if (st->nbLostCmpt > 1)
        {
            core = st->last_core_bfi;
        }

        /* no FD TCX PLC after a TCX transition frame: the appropriate framing is not implemented */
        if (st->nbLostCmpt == 1)
        {
            st->tonal_mdct_plc_active = 0;
            if ( !(st->rf_flag && st->use_partial_copy && (st->rf_frame_type == RF_TCXTD1 || st->rf_frame_type == RF_TCXTD2)))
            {
                if ((st->last_core == TCX_20_CORE)
                        && (st->second_last_core == TCX_20_CORE)
                        && ((st->old_fpitch <= 0.5f*st->L_frame) || (st->tcxltp_last_gain_unmodified <= 0.4f))
                        /* it is fine to call the detection even if no ltp information
                           is available, meaning that st->old_fpitch ==
                           st->tcxltp_second_last_pitch == st->L_frame */
                        && (st->old_fpitch == st->tcxltp_second_last_pitch)
                        && !st->last_tns_active
                        && !st->second_last_tns_active)
                {

                    TonalMDCTConceal_Detect(&st->tonalMDCTconceal,
                                            (st->tcxltp_last_gain_unmodified > 0) ? st->old_fpitch : 0,
                                            &numIndices);
                    if ((numIndices > 10)
                            || ((numIndices > 5)
                                && (fabs(st->tcxltp_third_last_pitch-st->tcxltp_second_last_pitch) < 0.5f)
                               )
                            || ((numIndices > 0) && ((st->last_good <= UNVOICED_TRANSITION) || (st->tcxltp_last_gain_unmodified <= 0.4f))
                                && (fabs(st->tcxltp_third_last_pitch-st->tcxltp_second_last_pitch) < 0.5f)
                               ))
                    {
                        core = 1;
                        st->tonal_mdct_plc_active = 1;
                    }
                    else if (st->last_good <= UNVOICED_TRANSITION || st->tcxltp_last_gain_unmodified <= 0.4f)
                    {
                        core = 1;
                    }
                }
                else if (st->last_core != ACELP_CORE)
                {
                    if (st->last_good <= UNVOICED_TRANSITION || st->tcxltp_last_gain_unmodified <= 0.4f)
                    {
                        core = st->last_core;
                    }
                }
            }
        }
    }
    return core;
}
