/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_enc.h"
#include "rom_com.h"
#include "prot.h"
#include <assert.h>

/*-----------------------------------------------------------------*
 * Local constants
 *-----------------------------------------------------------------*/

#define MDCT_SW_SIG_LINE_THR       2.85f  /* Significant spectral line threshold above Etot (dB)                               */
#define MDCT_SW_SIG_PEAK_THR       36.0f  /* Significant peak threshold below Etot (dB)                                        */
#define MDCT_SW_HI_SPARSE_THR      0.25f  /* Max. ratio of significant spectral lines for the spectrum to be considered sparse */
#define MDCT_SW_HI_ENER_LO_THR     7.5f   /* Hi band low energy threshold (dB)                                                 */
#define MDCT_SW_1_VOICING_THR      0.9f   /* Voicing threshold                                                                 */
#define MDCT_SW_1_HI_ENER_LO_THR   12.5f  /* Hi band low energy threshold (dB)                                                 */
#define MDCT_SW_1_SIG_HI_LEVEL_THR 28.0f  /* High signal level threshold above noise floor (dB)                                */
#define MDCT_SW_1_SIG_LO_LEVEL_THR 22.5f  /* Low signal level threshold above noise floor (dB)                                 */
#define MDCT_SW_1_COR_THR          80.0f  /* Threshold on cor_map_sum to indicate strongly tonal signal                        */
#define MDCT_SW_1_SPARSENESS_THR   0.65f  /* Threshold on spectrum sparseness                                                  */

#define MDCT_SW_2_VOICING_THR      0.6f   /* Voicing threshold                                                                 */
#define MDCT_SW_2_HI_ENER_LO_THR   9.5f   /* Hi band low energy threshold (dB)                                                 */
#define MDCT_SW_2_SIG_HI_LEVEL_THR 19.0f  /* High signal level threshold above noise floor (dB)                                */
#define MDCT_SW_2_SIG_LO_LEVEL_THR 23.5f  /* Low signal level threshold above noise floor (dB)                                 */
#define MDCT_SW_2_COR_THR          62.5f  /* Threshold on cor_map_sum to indicate strongly tonal signal                        */
#define MDCT_SW_2_SPARSENESS_THR   0.4f   /* Threshold on spectrum sparseness                                                  */

#define MDCT_SW_HYST_FAC           0.8f   /* Hysteresis tolerance factor                                                       */

#define LOG_10                     2.30258509299f

/*--------------------------------------------------------------------------*
 * get_sparseness()
 *
 *
 *--------------------------------------------------------------------------*/

static float get_sparseness(
    const float Bin_E[],
    short n,
    float thr
)
{
    short num_max, i;

    thr *= log(10); /* Convert to 10*log() domain from 10*log10() domain */

    thr = max(thr, 3.0f); /* Set an absolute minimum for close to silent signals */

    num_max = 0;
    for (i=1; i<n-1; ++i)
    {
        if (Bin_E[i] > max(max(Bin_E[i-1], Bin_E[i+1]), thr))
        {
            ++num_max;
        }
    }

    return 1.0f - num_max / (float)((n-2)/2);
}

/*--------------------------------------------------------------------------*
 * MDCT_selector()
 *
 *
 *--------------------------------------------------------------------------*/

void MDCT_selector(
    Encoder_State *st,        /* i/o: Encoder State           */
    float sp_floor,           /* i  : Noise floor estimate    */
    float Etot,               /* i  : Total energy            */
    float cor_map_sum,        /* i  : harmonicity factor      */
    const float voicing[],    /* i  : voicing factors         */
    const float enerBuffer[], /* i  : CLDFB buffers           */
    short vadflag
)
{
    if( st->mdct_sw_enable == MODE1 || st->mdct_sw_enable == MODE2 )
    {
        float hi_ener, frame_voicing, sparseness;
        int peak_count;
        short prefer_tcx, prefer_hq_core, switching_point, hi_sparse, sparse;
        int lob_cldfb, hib_cldfb, lob_fft, hib_fft;
        int i;
        float sig_lo_level_thr, sig_hi_level_thr, cor_thr, voicing_thr, sparseness_thr, hi_ener_lo_thr;
        int last_core;

        if (st->bwidth == NB)
        {
            lob_cldfb = 3200/400;
            hib_cldfb = 4000/400;
            lob_fft = (L_FFT/2)/2; /* 3.2 KHz */
            hib_fft = (40*(L_FFT/2))/64; /* 4.0 KHz */
        }
        else if (st->bwidth == WB)
        {
            lob_cldfb = 4800/400;
            hib_cldfb = 8000/400;
            lob_fft = 3*L_FFT/2/4; /* 4.8 KHz */
            hib_fft = L_FFT/2; /* 6.4 KHz (should be 8 KHz) */
        }
        else
        {
            lob_cldfb = 6400/400;
            hib_cldfb = 16000/400;
            if (st->bwidth == FB)
            {
                hib_cldfb = 24000/400;
            }
            lob_fft = L_FFT/2; /* 6.4 KHz */
            hib_fft = L_FFT/2; /* 6.4 KHz (should be 8 KHz) */
        }

        /* st->last_core is reset to TCX_20_CORE in init_acelp() => fix it here */
        last_core = st->last_core;
        if( st->last_codec_mode == MODE1 && last_core == TCX_20_CORE )
        {
            last_core = HQ_CORE;
        }

        /* Voicing */
        frame_voicing = (voicing[0]+voicing[1])*0.5f;

        /* Spectral sparseness */
        sparseness = get_sparseness(st->Bin_E, lob_fft, Etot - MDCT_SW_SIG_PEAK_THR);

        /* Hi band energy */
        hi_ener = (float)log10(mean(&enerBuffer[lob_cldfb], hib_cldfb-lob_cldfb)+0.0001f);

        /* Hi band sparseness */
        if (st->bwidth >= SWB)
        {
            /* For SWB, assume hi band sparseness based on 4.8 KHz-6.4 KHz band */
            lob_fft = 3*L_FFT/2/4; /* 4.8 KHz */
        }

        peak_count = 0;
        for (i=lob_fft; i<hib_fft; ++i)
        {
            if (st->Bin_E[i] >= Etot + MDCT_SW_SIG_LINE_THR*LOG_10)
            {
                ++peak_count;
            }
        }
        hi_sparse = peak_count <= anint((hib_fft - lob_fft) * MDCT_SW_HI_SPARSE_THR);
        sparse    = peak_count <= anint((hib_fft - lob_fft) * MDCT_SW_HI_SPARSE_THR/MDCT_SW_HYST_FAC);

        /* Hysteresis */
        if (st->prev_hi_sparse > 0 && sparse > 0 && min(min(voicing[0], voicing[1]), voicing[2]) >= MDCT_SW_1_VOICING_THR)
        {
            hi_sparse = 1;
        }

        /* Allowed switching point? */
        switching_point = (last_core != HQ_CORE && last_core != TCX_20_CORE) || /* previous core was non-MDCT */
                          (st->prev_hi_ener <= MDCT_SW_HI_ENER_LO_THR || hi_ener <= MDCT_SW_HI_ENER_LO_THR) || /* hi band is close to silent */
                          (last_core == HQ_CORE && (st->mdct_sw_enable == MODE1 || (hi_sparse > 0 && st->prev_hi_sparse >= 0 && st->prev_hi_sparse <= 1))) || /* HQ_CORE and hi band became sparse */
                          (last_core == TCX_20_CORE && (hi_sparse == 0 && st->prev_hi_sparse > 0)); /* TCX and hi band became dense */

        if (st->mdct_sw_enable == MODE1)
        {
            sig_lo_level_thr = MDCT_SW_1_SIG_LO_LEVEL_THR;
            sig_hi_level_thr = MDCT_SW_1_SIG_HI_LEVEL_THR;
            cor_thr          = MDCT_SW_1_COR_THR;
            voicing_thr      = MDCT_SW_1_VOICING_THR;
            sparseness_thr   = MDCT_SW_1_SPARSENESS_THR;
            hi_ener_lo_thr   = MDCT_SW_1_HI_ENER_LO_THR;
        }
        else
        {
            /* st->mdct_sw_enable == MODE2 */
            sig_lo_level_thr = MDCT_SW_2_SIG_LO_LEVEL_THR;
            sig_hi_level_thr = MDCT_SW_2_SIG_HI_LEVEL_THR;
            cor_thr          = MDCT_SW_2_COR_THR;
            voicing_thr      = MDCT_SW_2_VOICING_THR;
            sparseness_thr   = MDCT_SW_2_SPARSENESS_THR;
            hi_ener_lo_thr   = MDCT_SW_2_HI_ENER_LO_THR;
        }

        prefer_tcx = (Etot - sp_floor >= sig_hi_level_thr) && /* noise floor is low */
                     (cor_map_sum >= cor_thr || frame_voicing >= voicing_thr || sparseness >= sparseness_thr) && /* strong tonal components */
                     (hi_ener <= hi_ener_lo_thr || hi_sparse > 0); /* high freqs have low energy or are sparse */

        prefer_hq_core = (Etot - sp_floor < sig_lo_level_thr) || /* noise floor is very high */
                         (cor_map_sum < cor_thr*MDCT_SW_HYST_FAC && frame_voicing < voicing_thr*MDCT_SW_HYST_FAC && sparseness < sparseness_thr*MDCT_SW_HYST_FAC) || /* too weak tonal components */
                         (st->mdct_sw_enable == MODE1 && !prefer_tcx && st->transientDetection.transientDetector.bIsAttackPresent == 1);

        /* Prefer HQ_CORE on transients */
        if ( st->mdct_sw_enable == MODE2 && st->transientDetection.transientDetector.bIsAttackPresent == 1 )
        {
            prefer_tcx = 0;
            prefer_hq_core = 1;
        }

        if (switching_point && (prefer_tcx || prefer_hq_core))
        {
            if (prefer_tcx)
            {
                st->core = TCX_20_CORE;
            }
            else /* prefer_hq_core */
            {
                st->core = HQ_CORE;
            }
        }
        else if (last_core == HQ_CORE || last_core == TCX_20_CORE)
        {
            st->core = last_core;
        }

        /* Prevent the usage of HQ_CORE on noisy-speech or inactive */
        if ( st->mdct_sw_enable==MODE2 && st->core == HQ_CORE && (st->flag_noisy_speech_snr == 1 || vadflag == 0 ) )
        {
            st->core = TCX_20_CORE;
        }


        /* Update memories */
        if (hi_sparse <= 0)
        {
            st->prev_hi_sparse = hi_sparse;
        }
        else
        {
            st->prev_hi_sparse += hi_sparse;
            if (st->prev_hi_sparse >= 2)
            {
                st->prev_hi_sparse = 2;
            }
        }
        st->prev_hi_ener = hi_ener;
    }

    return;
}

/*--------------------------------------------------------------------------*
 * MDCT_selector_reset()
 *
 *
 *--------------------------------------------------------------------------*/

void MDCT_selector_reset(
    Encoder_State *st          /* i/o: Encoder State */
)
{
    st->prev_hi_ener = 0;
    st->prev_hi_sparse = -1;

    return;
}
