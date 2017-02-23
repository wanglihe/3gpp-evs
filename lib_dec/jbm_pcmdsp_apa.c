/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

/*! @file jbm_pcmdsp_apa.c Adaptive Playout for Audio (apa). */

/* system headers */
#include <assert.h>
#include <math.h>
#include <stdlib.h>   /* malloc(), free() */
#include <stdio.h>
#include "options.h"
/* instrumentation */
/* local headers */
#include "jbm_pcmdsp_apa.h"
#include "jbm_pcmdsp_similarityestimation.h"
#include "jbm_pcmdsp_window.h"
#include "jbm_types.h"
#include "cnst.h"


/*
********************************************************************************
*                         LOCAL DATA DEFINITIONS AND PROTOTYPES
********************************************************************************
*/

/* maximum number of segments/iterations in extend_frm() */
#define MAXN 10

/* definition of state struct */
struct apa_state_t
{
    /* output buffer */
    int16_t  buf_out[APA_BUF];
    uint16_t l_buf_out;

    /* Hann window */
    Float win[APA_BUF];
    uint16_t l_halfwin;

    /* sampling rate [Hz] */
    uint16_t rate;

    /* length of a segment [samples] */
    uint16_t l_seg;

    /* length of a frame [samples] */
    uint16_t l_frm;

    /* total number of processed input samples since apa_reset() */
    uint32_t l_in_total;

    /* sum of inserted/removed samples since last apa_set_scale() */
    int32_t  diffSinceSetScale;
    /* number of input frames since last apa_set_scale() */
    uint32_t nFramesSinceSetScale;

    /* current and previous  scaling ratio [%] */
    uint16_t scale;

    /* minimum pitch length [samples] */
    uint16_t p_min;

    /* search length [samples] */
    uint16_t l_search;

    uint16_t wss;                   /* waveform subsampling */
    uint16_t css;                   /* correlation subsampling */

    Float targetQuality;
    uint16_t qualityred;           /* quality reduction threshold */
    uint16_t qualityrise;          /* quality rising for adaptive quality thresholds */

    uint16_t last_pitch;           /* last pitch/sync position */
    uint16_t bad_frame_count;      /* # frames before quality threshold is lowered */
    uint16_t good_frame_count;     /* # scaled frames */

    uint16_t num_channels;         /* number of input/output channels */
};


/* prototypes for local functions */

/** Converts the correlation energy to dB. */
float apa_corrEnergy2dB(float energy, uint16_t corr_len);

/** Increases the calculated quality of signals with low energy. */
float apa_getQualityIncreaseForLowEnergy(float energydB);

static bool_t logarithmic_search(const apa_state_t * ps,
                                 const int16_t * signal,
                                 int16_t s_start,
                                 uint16_t inlen,
                                 uint16_t offset,
                                 uint16_t fixed_pos,
                                 uint16_t corr_len,
                                 uint16_t wss,
                                 uint16_t css,
                                 int16_t * synchpos);

static bool_t find_synch (apa_state_t * ps,
                          const int16_t * in,
                          uint16_t l_in,
                          int16_t s_start,
                          uint16_t s_len,
                          int16_t fixed_pos,
                          uint16_t corr_len,
                          uint16_t offset,
                          Float * energy,
                          Float * quality,
                          int16_t * synch_pos);

static bool_t copy_frm (apa_state_t * ps,
                        const int16_t frm_in[],
                        int16_t frm_out[],
                        uint16_t * l_frm_out);

static bool_t shrink_frm (apa_state_t * ps,
                          const int16_t frm_in[],
                          uint16_t maxScaling,
                          int16_t frm_out[],
                          uint16_t * l_frm_out);

static bool_t extend_frm (apa_state_t * ps,
                          const int16_t frm_in[],
                          int16_t frm_out[],
                          uint16_t * l_frm_out);

/*
********************************************************************************
*                         PUBLIC PROGRAM CODE
********************************************************************************
*/

/* Allocates memory for state struct and initializes elements. */
uint8_t apa_init (apa_state_t ** pps)
{
    apa_state_t *ps = NULL;

    /* make sure pointer is valid */
    if(!pps)
    {
        return 1;
    }

    /* allocate state struct */
    ps = (apa_state_t *) malloc (sizeof (apa_state_t));
    if(!ps)
    {
        return 2;
    }
    apa_reset (ps);
    *pps = ps;
    return 0;
}

/* Sets state variables to initial value. */
void apa_reset (apa_state_t * ps)
{
    /* init state struct */
    ps->l_buf_out = 0;
    ps->l_halfwin = 0;
    ps->rate = 0;
    ps->l_seg = 0;
    ps->l_frm = 0;
    ps->l_in_total = 0;
    ps->diffSinceSetScale = 0;
    ps->nFramesSinceSetScale = 0;
    ps->scale = 100;
    ps->p_min = 0;
    ps->l_search = 0;
    ps->wss = 1;
    ps->css = 1;
    ps->targetQuality = 0.0f;
    ps->qualityred = 0;
    ps->qualityrise = 0;
    ps->last_pitch = 0;
    ps->bad_frame_count = 0;
    ps->good_frame_count = 0;
    ps->num_channels = 0;
}

/* Sets the audio configuration. */
bool_t apa_set_rate (apa_state_t * ps, uint16_t rate, uint16_t num_channels)
{
    /* make sure pointer is valid */
    if (ps == (apa_state_t *) NULL)
    {
        return 1;
    }

    /* check range */
    if ((rate < APA_MIN_RATE) || (rate > APA_MAX_RATE))
    {
        return 1;
    }

    /* reset state struct */
    apa_reset (ps);

    /* copy rate to state struct */
    ps->rate = rate;

    /* set number of channels */
    ps->num_channels = num_channels;

    /*
     * several other parameters depend on the sampling rate
     * and are set below. Some "magic numbers" are used here
     * which are based on typical values of a "pitch" in
     * human voice. The pitch length is the period of the
     * base frequency and is usually assumed to be 40-240
     * samples at 16 kHz.
     */

    /* set segment size */
    /* in the order of a pitch, set to 160 samples at 16 kHz */
    /* used for windowing and as the correlation length, i.e., */
    /* the size of the template segment. */
    ps->l_seg = (uint16_t) ((Float) ps->rate / 100.0f) * ps->num_channels;

    /* init Hann window */
    /* Note: l_win < APA_BUF is required, which is assured */
    /* because APA_MAX_RATE/100 = l_win = 441 < 2048 = APA_BUF */
    /* Length of Hann window should be independent of
     * number of channels - same window applied to all channels */
    ps->l_halfwin = ps->rate / 100;
    hannWindow (ps->l_halfwin*2, ps->win);

    /* set frame size */
    /* set to 320 samples at 16 kHz */
    ps->l_frm = (uint16_t) ((Float) ps->rate / 50.0f) * ps->num_channels;

    /* set minimum pitch */
    /* set to 40 samples at 16 kHz */
    /* (defines min change in number of samples, i.e., abs(l_in-l_out) >= p_min) */
    ps->p_min = (uint16_t) ((Float) ps->rate / 400.0f) * ps->num_channels;

    /* set search length */
    /* must cover one pitch, set to 200 samples at 16 kHz */
    /* (the resulting maximum pitch is then p_min+l_search = 240 samples at 16 kHz) */
    ps->l_search = (uint16_t) ((Float) ps->rate / 80.0f) * ps->num_channels;

    return 0;
}

/* Set scaling. */
bool_t apa_set_scale (apa_state_t * ps, uint16_t scale)
{
    /* make sure pointer is valid */
    if (ps == (apa_state_t *) NULL)
    {
        return 1;
    }

    /* check range */
    if ((scale < APA_MIN_SCALE) || (scale > APA_MAX_SCALE))
    {
        return 1;
    }

    /* do nothing if same scale is set multiple times */
    /* (otherwise scale control is confused) */
    if (ps->scale == scale)
    {
        return 0;
    }

    /* copy to state struct */
    ps->scale = scale;

    /* reset scaling statistics */
    ps->diffSinceSetScale = 0;
    ps->nFramesSinceSetScale = 0;

    return 0;
}

/*
********************************************************************************
*
*     Function        : apa_set_quality
*     Tables          : <none>
*     Compile Defines : <none>
*     Return          : 0 on success, 1 on failure
*     Information     : Set quality thresholds.
*
*                       quality is lower limit for minimum quality
*                       Range is [-2;2] - where positive values allow
*                       only pasting with same phase information
*                       Negative values would yield cross phased pasting
*
*                       qualityred allows dynamic lowering of lower quality
*                       bound - this gives better results for rhythmic signals
*                       Range is [0;20], meaning 0.1 lowering*qualityred
*
*                       undocumented: qualityrise (same as qualityred - other
*                       direction)
*
********************************************************************************
*/
bool_t apa_set_quality(
    apa_state_t *ps,
    float        quality,
    uint16_t     qualityred,
    uint16_t     qualityrise)
{
    assert(ps != (apa_state_t *) NULL);
    assert(-2.0f <= quality && quality <= 3.1f);
    assert(qualityred > 0 && qualityred <= 20);
    assert(qualityrise > 0 && qualityrise <= 20);

    ps->targetQuality = quality;
    ps->qualityred = qualityred;
    ps->qualityrise = qualityrise;
    ps->bad_frame_count = 0;
    ps->good_frame_count = 0;
    return 0;
}

/*
********************************************************************************
*
*     Function        : apa_set_complexity_options
*     Tables          : <none>
*     Compile Defines : <none>
*     Return          : 0 on success, 1 on failure
*     Information     : Set complexity options
*                       Waveform subsampling computes the correlation function
*                         for certain positions only
*                       Correlation function subsampling computes the maxima
*                         for certain positions only
*
********************************************************************************
*/
bool_t apa_set_complexity_options (apa_state_t * ps, uint16_t wss, uint16_t css)
{
    /* make sure pointer is valid */
    if (ps == (apa_state_t *) NULL)
        return 1;
    if (wss == 0 || wss > 1000)
        return 1;
    if (css == 0 || css > 1000)
        return 1;
    ps->wss = wss;
    ps->css = css;
    return 0;
}

/*
********************************************************************************
*
*     Function        : apa_exit
*     Tables          : <none>
*     Compile Defines : <none>
*     Return          : 0 on success, 1 on failure
*     Information     : The memory used for storing the state is freed.
*                       The state struct pointer is set to NULL.
*
********************************************************************************
*/
bool_t apa_exit (apa_state_t ** pps)
{
    /* ignore NULL pointer input */
    if (*pps == (apa_state_t *) NULL)
    {
        return 0;
    }

    /* deallocate state struct */
    free (*pps);
    /* set pointer to NULL */
    *pps = NULL;
    return 0;
}

/*
********************************************************************************
*
*     Function        : apa_exec
*     Tables          : <none>
*     Compile Defines : <none>
*     Return          : 0 on success, 1 on failure
*     Information     : Execute adaptive playout for audio, i.e., audio scaling.
*                       Will take l_in input samples from a_in[] and
*                       try to extend/shrink the amount of samples according
*                       to the last scaling set by using apa_set_scale().
*                       The actual amount of samples after scaling may vary
*                       and is given in l_out. The scaled audio samples
*                       are contained in a_out[]. Note that the scaling is
*                       achieved only in average. The input buffer must be
*                       filled with 20ms audio. The output buffer must be
*                       allocated externally and must be at least of size
*                       APA_BUF.
*                       Scaling can only be performed when a sampling rate
*                       is specified using apa_set_rate(). Otherwise,
*                       an error is returned.
*
*                       The amount of scaling is achieved by controlling the
*                       frequency of scaling. Note that the exact amount of
*                       scaling is signal dependent and is an integer
*                       multiple of a pitch. Hence, when we want to achieve
*                       a scaling of e.g. 110% then the APA module will typically
*                       forward several frames without any modification and
*                       then scale one frame by a higher amount, e.g. 143%.
*
********************************************************************************
*/
uint8_t apa_exec (apa_state_t * ps, /* i/o: state struct */
                  const int16_t a_in[], /* i:   input samples */
                  uint16_t l_in, /* i:   number of input samples */
                  uint16_t maxScaling, /* i: allowed number of inserted/removed samples */
                  int16_t a_out[], /* o:   output samples */
                  uint16_t * l_out /* o:   number of output samples */)
{
    uint16_t i;
    int16_t frm_in[APA_BUF];
    uint16_t l_frm_out;
    int16_t l_rem;
    int32_t dl_scaled, dl_copied, l_frm_out_target;
    int32_t expScaling, actScaling;
    uint32_t statsResetThreshold, statsResetShift;

    statsResetThreshold = 1637;
    statsResetShift     = 2;

    /* make sure no invalid output is used */
    *l_out = 0;
    l_frm_out = 0;

    /* make sure pointer is valid */
    if (ps == (apa_state_t *) NULL)
    {
        return 1;
    }
    /* check available rate */
    if (ps->rate == 0)
    {
        return 2;
    }
    /* check size of input */
    if( l_in != ps->l_frm )
    {
        return 3;
    }

    /* get target length */
    if(ps->scale > 100)
    {
        expScaling = (int32_t)((ps->l_frm * (ps->scale - 100.0f) / 100.0f) * (ps->nFramesSinceSetScale + 1) + 0.5f);
    }
    else if(ps->scale < 100)
    {
        expScaling = (int32_t)((ps->l_frm * (ps->scale - 100.0f) / 100.0f) * (ps->nFramesSinceSetScale + 1) - 0.5f);
    }
    else
    {
        expScaling = 0;
    }
    actScaling = ps->diffSinceSetScale - ps->l_frm;
    l_frm_out_target = expScaling - actScaling;

    /* Wait until we have l_frm outputs samples */
    /* (required to search for correlation in the past). */
    /* If we don't have enough samples, simply copy input to output */
    if (ps->l_buf_out < ps->l_frm)
    {
        for (i = 0; i < ps->l_frm; i++)
        {
            a_out[i] = a_in[i];
        }
        l_frm_out = ps->l_frm;
    }
    else
    {
        int16_t * buf_out_ptr = &(ps->buf_out[ps->l_buf_out - ps->l_frm]);
        int16_t * frm_in_ptr = &(frm_in[ps->l_frm]);

        /* fill input frame */
        /* 1st input frame: previous output samples */
        for (i = 0; i < ps->l_frm; i++)
        {
            frm_in[i] = buf_out_ptr[i];
        }
        /* 2nd input frame: new input samples */
        for (i = 0; i < ps->l_frm; i++)
        {
            frm_in_ptr[i] = a_in[i];
        }
        /* no scaling */
        if (ps->scale == 100)
        {
            copy_frm (ps, frm_in, a_out, &l_frm_out);
        }
        /* shrink */
        else if (ps->scale < 100)
        {
            shrink_frm (ps, frm_in, maxScaling, a_out, &l_frm_out);
        }
        /* extend */
        else
        {
            extend_frm (ps, frm_in, a_out, &l_frm_out);
        }
        /* control the amount/frequency of scaling */
        if( l_frm_out != ps->l_frm )
        {
            if( maxScaling != 0U &&
                    abs( (int16_t)(ps->l_frm - l_frm_out) ) > maxScaling )
            {
                /* maxScaling exceeded -> discard scaled frame */
                copy_frm (ps, frm_in, a_out, &l_frm_out);
            }
            else if( abs(l_frm_out_target) > ps->l_frm )   /* ignore small difference */
            {
                dl_copied = l_frm_out_target - (int32_t) ps->l_frm;
                dl_scaled = l_frm_out_target - (int32_t) l_frm_out;
                /* discard scaled frame if copied frame is closer to target length */
                if( abs(dl_copied) < abs(dl_scaled) )
                {
                    copy_frm (ps, frm_in, a_out, &l_frm_out);
                }
            }
        }
    }

    /* copy output to internal buffer */
    /* avoid buffer overflow: */
    /* discard old samples; always keep at least most recent l_frm samples */
    if ((ps->l_buf_out + l_frm_out) > APA_BUF)
    {
        int16_t * buf_out_ptr1 = ps->buf_out;
        int16_t * buf_out_ptr2;

        l_rem = (ps->l_frm - l_frm_out);
        if (l_rem < 0)
        {
            l_rem = 0;
        }
        buf_out_ptr2 = &(ps->buf_out[ps->l_buf_out - l_rem]);
        for (i = 0; i < l_rem; i++)
        {
            buf_out_ptr1[i] = buf_out_ptr2[i];
        }
        ps->l_buf_out = l_rem;
    }
    /* append new output samples */
    if ((ps->l_buf_out + l_frm_out) > APA_BUF)
    {
        return 5;
    }
    {
        int16_t * buf_out_ptr = &(ps->buf_out[ps->l_buf_out]);
        for (i = 0; i < l_frm_out; i++)
        {
            buf_out_ptr[i] = a_out[i];
        }
    }
    ps->l_buf_out += l_frm_out;

    *l_out = l_frm_out;
    /* update time */
    ps->l_in_total += ps->l_frm;

    if( abs(ps->diffSinceSetScale) < (0x7FFFFF - (l_frm_out - ps->l_frm)) &&
            ps->nFramesSinceSetScale < statsResetThreshold )
    {
        ps->diffSinceSetScale += l_frm_out - ps->l_frm;
        ++ps->nFramesSinceSetScale;
    }
    else   /* scale statistics down to avoid overflow */
    {
        ps->diffSinceSetScale    >>= statsResetShift;
        ps->nFramesSinceSetScale >>= statsResetShift;
    }

    return 0;
}


/*
********************************************************************************
*                         LOCAL PROGRAM CODE
********************************************************************************
*/


/*
********************************************************************************
*
*     Function        : get_scaling_quality
*     Tables          : <none>
*     Compile Defines : <none>
*     Return          : 0 on success, 1 on failure
*     Information     : Uses pitch, half pitch, three halves and double pitch
*                       to evaluate the quality of the scaled frame by checking
*                       periodicity.
*                       Silence can be detected as additional feature. This must
*                       be set in global struct apa_state.
*
*                       If search length is very narrow then use fewer points
*                       to evaluate periodicity and silence.
*
*                       Computationally not very efficient by using normalized
*                       cross-correlation: Using sqrt() for energy calculation
*                       adds complexity.
*
*     03-AUG-04  S.Doehla        initial version
*
********************************************************************************
*/
static void get_scaling_quality(const apa_state_t * ps,
                                const int16_t * signal,
                                uint16_t s_len,
                                uint16_t offset,
                                uint16_t corr_len,
                                uint16_t pitch,
                                Float * energydB,
                                Float * quality)
{
    Float maxEnergy = 0.0f;
    Float qualityOfMaxEnergy = 0.0f;  /* we measure the quality for all channels and select the one with highest energy */

    Float half_pitch_cn = 0.0f;
    Float pitch_cn = 0.0f;
    Float three_halves_pitch_cn = 0.0f;
    Float double_pitch_cn = 0.0f;

    Float pitch_energy = 0.0f;
    Float half_pitch_energy = 0.0f;
    Float three_halves_pitch_energy = 0.0f;
    Float double_pitch_energy = 0.0f;

    uint16_t i = 0;

    for(i=0; i<ps->num_channels; i++)
    {
        Float energy;
        offset = 0;

        pitch_cn = normalized_cross_correlation_self(signal, pitch+offset, offset, corr_len,
                   ps->num_channels * 2, &pitch_energy);
        if(pitch_cn > 0.0f)
        {
            /* calculate correlation for double pitch */
            if( 2*pitch + offset + corr_len <= s_len )
            {
                double_pitch_cn = normalized_cross_correlation_self(signal, 2*pitch + offset,
                                  offset, corr_len, ps->num_channels * 2, &double_pitch_energy);
            }
            else
            {
                double_pitch_cn     = pitch_cn;
                double_pitch_energy = pitch_energy;
            }
            /* calculate correlation for three/half pitch */
            if( (3*pitch)/2 + offset + corr_len <= s_len )
            {
                three_halves_pitch_cn = normalized_cross_correlation_self(signal, (3*pitch)/2 + offset,
                                        offset, corr_len, ps->num_channels * 2, &three_halves_pitch_energy);
            }
            else
            {
                three_halves_pitch_cn     = pitch_cn;
                three_halves_pitch_energy = pitch_energy;
            }
            /* calculate correlation for half pitch */
            if( pitch/2 + offset + corr_len <= s_len )
            {
                half_pitch_cn = normalized_cross_correlation_self(signal, pitch/2 + offset,
                                offset, corr_len, ps->num_channels * 2, &half_pitch_energy);
            }
            else
            {
                half_pitch_cn     = pitch_cn;
                half_pitch_energy = pitch_energy;
            }

            /* combine correlation results */
            *quality = (half_pitch_cn * three_halves_pitch_cn) + (pitch_cn * double_pitch_cn);
            energy = pitch_energy + half_pitch_energy + three_halves_pitch_energy + double_pitch_energy;
        }
        else
        {
            *quality = pitch_cn; /* value is negative, thus pass it */
            energy = pitch_energy;
        }

        /* update the quality by the quality of the signal with the highest energy */
        if(energy > maxEnergy)
        {
            qualityOfMaxEnergy = *quality;
            maxEnergy = energy;
        }

        /* go to next channel */
        ++signal;
    }
    *quality = qualityOfMaxEnergy;

    /* increase calculated quality of signals with low energy */
    *energydB = apa_corrEnergy2dB(maxEnergy, corr_len);
    *quality += apa_getQualityIncreaseForLowEnergy(*energydB);
}

/* Converts the correlation energy to dB. */
float apa_corrEnergy2dB(float energy, uint16_t corr_len)
{
    float energydB = 10.0f * (float)log10(energy / (32768.0f * 32768.0f * corr_len * 4.0f));
    return energydB;
}

/* Increases the calculated quality of signals with low energy. */
float apa_getQualityIncreaseForLowEnergy(float energydB)
{
    const float qualIncreaseMinEnergy = -65;
    const float qualIncreaseMaxEnergy = -40;
    float qualIncForLowEnergy = 0;
    if(energydB < qualIncreaseMaxEnergy)
    {
        qualIncForLowEnergy = energydB;
        if(qualIncForLowEnergy < qualIncreaseMinEnergy)
        {
            qualIncForLowEnergy = qualIncreaseMinEnergy;
        }
        if(qualIncForLowEnergy > qualIncreaseMaxEnergy)
        {
            qualIncForLowEnergy = qualIncreaseMaxEnergy;
        }
        qualIncForLowEnergy = (qualIncForLowEnergy - qualIncreaseMaxEnergy) /
                              (qualIncreaseMinEnergy - qualIncreaseMaxEnergy) * 2;
        assert(qualIncForLowEnergy >= 0 && qualIncForLowEnergy <= 2);
    }
    return qualIncForLowEnergy;
}

/*
********************************************************************************
*
*     Function        : logarithmic_search
*     Tables          : <none>
*     Compile Defines : <none>
*     Return          : 0 on success, 1 on failure
*     Information     : Search for best match of a template segment using
*                       hierarchical search method:
*                       Parameter css is used for sampling every css'd correlation
*                       value. The area around the best match so far is used for
*                       further correlation value with half css-value until css=1.
*                       Search area length is always half previous search length.
*                       Parameter wss is passed to the correlation computation
*                       If the search area passes the boundaries, the search
*                       window is reduced so that it's entirely inside the
*                       boundaries.
*
********************************************************************************
*/
static bool_t logarithmic_search(const apa_state_t * ps,
                                 const int16_t * signal,
                                 int16_t s_start,
                                 uint16_t inlen,
                                 uint16_t offset,
                                 uint16_t fixed_pos,
                                 uint16_t corr_len,
                                 uint16_t wss,
                                 uint16_t css,
                                 int16_t * synchpos)
{
    int i;
    Float coeff;
    Float coeff_max;
    int16_t s_start_old = 0;
    uint16_t s_len_old = 0;
    do
    {
        coeff_max = -FLT_MAX; /* will always be overwritten with result of first correlation */
        for (i = s_start; i < s_start+inlen; i += css)
        {
            if((wss == 1) && (ps->num_channels == 1))
            {
                coeff = cross_correlation_self(signal, i+offset, fixed_pos+offset, corr_len);
            }
            else
            {
                coeff = cross_correlation_subsampled_self( signal, i+offset, fixed_pos+offset,
                        corr_len, wss*ps->num_channels);
            }

            /* update max corr */
            if( ps->scale < 100 )
            {
                /* shrinking: prefer greater synchpos for equal coeff */
                if (coeff >= coeff_max)
                {
                    coeff_max = coeff;
                    *synchpos = i;
                }
            }
            else
            {
                /* extending: prefer smaller synchpos for equal coeff */
                if (coeff > coeff_max)
                {
                    coeff_max = coeff;
                    *synchpos = i;
                }
            }
        }
        /* backup old search range */
        s_start_old = s_start;
        s_len_old = inlen;

        css = css/2;
        inlen = inlen/2;
        s_start = *synchpos - inlen/2;
        if(s_start < s_start_old)
        {
            s_start = s_start_old;
        }
        if((s_start+inlen) > (s_start_old+s_len_old))
        {
            inlen = s_start_old - s_start + s_len_old;
        }
    }
    while (css > 2);

    return 0;
}


/*
********************************************************************************
*
*     Function        : find_synch
*     Tables          : <none>
*     Compile Defines : <none>
*     Return          : 0 on success, 1 on failure
*     Information     : Find the best match of an template segment within
*                       a search region by similarity measures.
*
*                       Typical example:
*
*         0         10        20        30        40        50        60
*  in[] = abcdefghijk_abcdefghijk_abcdefghijk_abcdEFGHIJk_abcdefghijk_a
*  l_in = 61
*  offset = 30                          |
*  s_start = -20    <-------------------|
*  s_len   = 15     <------------->     |
*  search range:    ***************     |
*  fixed_pos = 10                       |--------->
*  corr_len = 6                         |         <---->
*  template segment:                    |         ******
*  synch_pos: -14         <-------------|
*
*                       All positions are given relative to offset. The
*                       search region starts at offset+s_start and ends
*                       at offset+s_start+s_len. The template segment
*                       starts at offset+fixed_pos and ends at
*                       offset+fixed_pos+corr_len. For correlation, the
*                       template segment (EFGHIJ) is matched against the
*                       segment in the search region, e.g., against (k_abcd)
*                       in the first search position. The search position
*                       with the best match (-14: EFGHIJ <-> efghij) is
*                       returned.
*
*      19-JUN-03  N.Faerber       initial version
*      23-APR-04  S.Doehla        added subsampling
*
********************************************************************************
*/
static bool_t find_synch (apa_state_t * ps,
                          const int16_t * in,
                          uint16_t l_in,
                          int16_t s_start,
                          uint16_t s_len,
                          int16_t fixed_pos,
                          uint16_t corr_len,
                          uint16_t offset,
                          Float * energy,
                          Float * quality,
                          int16_t * synch_pos)
{
    assert( (corr_len - 1 + s_start + s_len - 1 + offset) < l_in );
    assert( (corr_len - 1 + fixed_pos + offset) < l_in);

    /* pass last pitch to search function as prediction value */
    *synch_pos = ps->last_pitch;

    logarithmic_search(ps,
                       in,
                       s_start,
                       s_len,
                       offset,
                       fixed_pos,
                       corr_len,
                       ps->wss,
                       ps->css * ps->num_channels,
                       synch_pos);
    /* assert synch_pos is cleanly divisible by number of channels */
    assert( *synch_pos % ps->num_channels == 0 );

    *quality = 0;
    get_scaling_quality(ps, in, l_in, offset, corr_len,
                        abs(fixed_pos - *synch_pos), energy, quality);
    ps->last_pitch = *synch_pos;
    return 0;
}


/*
********************************************************************************
*
*     Function        : copy_frm
*     Tables          : <none>
*     Compile Defines : <none>
*     Return          : 0 on success, 1 on failure
*     Information     : Copy an audio.
*
*                       The frame size is fixed to ps->l_frm. The input data
*                       is stored in frm_in[], where the first ps->l_frm samples
*                       shall include the previous output frame and the second
*                       ps->l_frm samples shall contain the current input frame.
*                       The output frame is stored in frm_out[] and contains
*                       l_frm_out = ps->l_frm.
*
*                       The first ps->l_frm input samples are not used by
*                       this function and are only provided for a consistent
*                       function call with shrink_frm() and extend_frm().
*
********************************************************************************
*/
static bool_t copy_frm (apa_state_t * ps,
                        const int16_t frm_in[], int16_t frm_out[], uint16_t * l_frm_out)
{
    uint16_t i;

    /* only 2nd input frame is used */
    frm_in += ps->l_frm;

    /* copy frame */
    for (i = 0; i < ps->l_frm; i++)
    {
        frm_out[i] = frm_in[i];
    }

    /* set output length */
    *l_frm_out = ps->l_frm;

    return 0;
}

/*
********************************************************************************
*
*     Function        : shrink_frm
*     Tables          : <none>
*     Compile Defines : <none>
*     Return          : 0 on success, 1 on failure
*     Information     : Shrink the length of an audio frame using the WSOLA
*                       algorithm.
*
*                       The frame size is fixed to ps->l_frm. The input data
*                       is stored in frm_in[], where the first ps->l_frm samples
*                       shall include the previous output frame and the second
*                       ps->l_frm samples shall contain the current input frame.
*                       The output frame is stored in frm_out[] and contains
*                       l_frm_out samples. The amount of shrinking is signal
*                       dependent.
*
*                       The first ps->l_frm input samples are not used by
*                       this function and are only provided for a consistent
*                       function call with extend_frm().
*
********************************************************************************
*/
static bool_t shrink_frm (apa_state_t * ps,
                          const int16_t frm_in[], uint16_t maxScaling, int16_t frm_out[], uint16_t * l_frm_out)
{
    bool_t findSynchResult = 0;
    int16_t xtract, l_rem, s_start, s_end;
    uint16_t i;
    uint16_t over;
    Float energy, quality = 0.0f;
    uint16_t l_frm;
    uint16_t l_seg;

    l_frm = ps->l_frm;
    l_seg = ps->l_seg;

    /* only 2nd input frame is used */
    frm_in += l_frm;

    /* set search range */
    s_start = (ps->p_min / ps->num_channels) * ps->num_channels;
    s_end = s_start + ps->l_search;
    if ((s_end + l_seg) >= l_frm)
    {
        s_end = (l_frm - l_seg);
    }

    /* calculate overlap position */
    if( isSilence( frm_in, l_seg, 10 ) )
    {
        /* maximum scaling */
        energy  = -65;
        quality = 5;
        if( maxScaling != 0U && s_end > maxScaling + 1 )
        {
            xtract = maxScaling;
        }
        else
        {
            /* set to last valid element (i.e. element[len - 1] but note for stereo last element is last pair of samples) */
            xtract = s_end - ps->num_channels;
        }
    }
    else
    {
        /* find synch */
        findSynchResult = find_synch (ps, frm_in, l_frm,
                                      s_start, (uint16_t) (s_end - s_start), 0,
                                      l_seg, 0, &energy, &quality, &xtract);
    }
    /* assert synch_pos is cleanly divisible by number of channels */
    assert( xtract % ps->num_channels == 0 );

    /* set frame overlappable - reset if necessary */
    over = 1;

    /* test whether frame has sufficient quality */
    if(quality < (ps->targetQuality - (ps->bad_frame_count * 0.1f) + (ps->good_frame_count * 0.2f)))
    {
        /* not sufficient */
        over = 0;
        if( ps->bad_frame_count < ps->qualityred )
        {
            ++ps->bad_frame_count;
        }
        if( ps->good_frame_count > 0U )
        {
            --ps->good_frame_count;
        }
    }
    else
    {
        /* sufficient quality */
        if( ps->bad_frame_count > 0U )
        {
            --ps->bad_frame_count;
        }
        if( ps->good_frame_count < ps->qualityrise )
        {
            ++ps->good_frame_count;
        }
    }

    /* Calculate output data */
    if(over && xtract)
    {
        if (findSynchResult == 1)
        {
            return 1;
        }
        overlapAdd(frm_in, frm_in + xtract, frm_out, l_seg, ps->num_channels,
                   ps->win + ps->l_halfwin, ps->win);
    }
    else
    {
        xtract = 0;
        for (i = 0; i < l_seg; i++)
        {
            frm_out[i] = frm_in[i];
        }
    }

    /* append remaining samples */
    l_rem = l_frm - xtract - l_seg;
    for (i = 0; i < l_rem; i++)
    {
        frm_out[l_seg + i] = frm_in[l_frm - l_rem + i];
    }

    /* set output length */
    *l_frm_out = l_seg + l_rem;

    return 0;
}

/*
********************************************************************************
*
*     Function        : extend_frm
*     Tables          : <none>
*     Compile Defines : <none>
*     Return          : 0 on success, 1 on failure
*     Information     : Extend the length of an audio frame using the WSOLA
*                       algorithm.
*
*                       The frame size is fixed to ps->l_frm. The input data
*                       is stored in frm_in[], where the first ps->l_frm samples
*                       shall include the previous output frame and the second
*                       ps->l_frm samples shall contain the current input frame.
*                       The output frame is stored in frm_out[] and contains
*                       l_frm_out samples. The amount of extension is signal
*                       dependent.
*
********************************************************************************
*/
static bool_t extend_frm (apa_state_t * ps,
                          const int16_t frm_in[], int16_t frm_out[], uint16_t * l_frm_out)
{
    bool_t findSynchResult = 0;
    uint16_t l_frm_out_target;
    uint16_t n, i;
    int16_t N;
    int16_t s[MAXN + 2], s_max, s_min;
    int16_t xtract[MAXN + 2], sync_start, s_end;
    uint16_t over[MAXN + 2];
    int16_t l_rem;
    int16_t s_start = 0;
    Float energy, quality = 0.0f;
    uint16_t l_frm, l_seg;
    const int16_t *fadeOut, *fadeIn;
    int16_t *out;

    l_frm = ps->l_frm;
    l_seg = ps->l_seg;

    /* number of segments/iterations */
    l_frm_out_target = (uint16_t) ((Float) l_frm * 1.5f);
    N = (l_frm_out_target / l_seg) - 1;
    if (N < 1)
    {
        N = 1;
    }
    if (N > MAXN)
    {
        return 1;
    }
    /* calculate equally spaced search regions */
    /* s[n] are given relative to 2nd frame and point to the start of */
    /* the search region. The first segment (n=1) will not be moved. */
    /* Hence, the iterations will start with n=2. */
    s_min = -(ps->l_search) - (ps->p_min);
    /* (make sure not to exceed array dimension) */
    if (l_frm + s_min < 0)
    {
        s_min = -(l_frm);
    }
    s_max = l_frm - 2 * l_seg - ps->l_search;
    if (s_max < s_min)
    {
        N = 1;
    }
    /* for just one segment start at s_min */
    if (N == 1)
    {
        s[2] = s_min;
    }
    /* else, spread linear in between s_min and s_max */
    /* (including s_min and s_max) */
    else
    {
        for (n = 2; n <= (N + 1); n++)
        {
            s[n] = s_min + ((s_max - s_min) * (n - 2)) / (N - 1);
        }
    }

    /*
     *  Planning Phase
     */

    xtract[1] = -(l_seg); /* make sync_start=0 in 1st iteration */
    n = 2;
    {
        /* define synch segment (to be correlated with search region) */
        sync_start = xtract[n - 1] + l_seg;
        over[n] = 1; /* will be reset if overlap is not required */
        /* check end of search region: should be at least p_min */
        /* samples on the left of synch_start */
        if ((s[n] + ps->l_search) < (sync_start - (ps->p_min)))
        {
            s_start = s[n];
            s_end = s_start + ps->l_search;
        }
        else
        {
            /* shrink search region to enforce minimum shift */
            s_end = sync_start - (ps->p_min);
            if (s[n] + ps->l_search < sync_start)
            {
                s_start = s[n]; /* just do it with normal start position */
            }
            else if (n == (N + 1))   /* move search region left for last segment */
            {
                s_start = s_end - (ps->l_search - ps->p_min);
            }
            else
            {
                over[n] = 0; /* don't search/overlap (just copy down) */
            }
        }

        if (over[n])
        {
            /* calculate overlap position */
            if( isSilence( frm_in, l_seg, 10 ) )
            {
                /* maximum scaling */
                energy  = -65;
                quality = 5;
                xtract[n] = s_start + ps->num_channels;
            }
            else
            {
                /* find synch */
                findSynchResult = find_synch( ps, frm_in, 2 * l_frm,
                                              s_start, s_end - s_start, sync_start,
                                              l_seg, l_frm, &energy, &quality, &xtract[n] );
            }
            /* assert synch_pos is cleanly divisible by number of channels */
            assert( xtract[n] % ps->num_channels == 0 );

            /* test for sufficient quality */
            if(quality < (ps->targetQuality - (ps->bad_frame_count * 0.1f) + (ps->good_frame_count * 0.2f)))
            {
                /* not sufficient */
                over[n] = 0;
                xtract[n] = sync_start;
                if( ps->bad_frame_count < ps->qualityred )
                {
                    ++ps->bad_frame_count;
                }
                if( ps->good_frame_count > 0U )
                {
                    --ps->good_frame_count;
                }
            }
            else
            {
                /* sufficient quality */
                if( ps->bad_frame_count > 0U )
                {
                    --ps->bad_frame_count;
                }
                if( ps->good_frame_count < ps->qualityrise )
                {
                    ++ps->good_frame_count;
                }
            }
            if (findSynchResult)
            {
                return 1;
            }
        }
        else
        {
            xtract[n] = sync_start;
        }

    }

    /* Calculate output data */
    for (n = 2; n <= N; n++)
    {
        if (over[n] && xtract[n-1] + l_seg != xtract[n])
        {
            /* mix 2nd half of previous segment with 1st half of current segment */
            fadeOut = frm_in + l_frm + xtract[n - 1] + l_seg;
            fadeIn  = frm_in + l_frm + xtract[n];
            out     = frm_out + (n - 2) * l_seg;
            overlapAdd(fadeOut, fadeIn, out, l_seg, ps->num_channels,
                       ps->win + ps->l_halfwin, ps->win);
        }
        else
        {
            /* just copy down 1st half of current segment (= 2nd half of previous segment) */
            int16_t * frm_out_ptr;
            const int16_t * frm_in_ptr;
            frm_out_ptr = &(frm_out[(n - 2) * l_seg]);
            frm_in_ptr  = &(frm_in[l_frm + xtract[n]]);
            for (i = 0; i < l_seg; i++)
            {
                frm_out_ptr[i] = frm_in_ptr[i];
            }
        }
    }

    /* append remaining samples */
    l_rem = l_frm - (xtract[N] + l_seg);
    for (i = 0; i < l_rem; i++)
    {
        frm_out[(N - 1) * l_seg + i] = frm_in[2 * l_frm - l_rem + i];
    }

    /* set output length */
    *l_frm_out = (N - 1) * l_seg + l_rem;

    return 0;
}

