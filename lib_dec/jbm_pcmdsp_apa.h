/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

/*! @file jbm_pcmdsp_apa.h Adaptive Playout for Audio (apa). */

#ifndef JBM_PCMDSP_APA_H
#define JBM_PCMDSP_APA_H JBM_PCMDSP_APA_H

#include "jbm_types.h"

/*
********************************************************************************
*                         DEFINITION OF CONSTANTS
********************************************************************************
*/

/* size of IO buffers (a_in[], a_out[]) for apa_exec() */
#define APA_BUF 4096*3

/* min/max sampling rate [Hz] */
#define APA_MIN_RATE 1000
#define APA_MAX_RATE 48000

/* min/max  scaling [%] */
#define APA_MIN_SCALE 50
#define APA_MAX_SCALE 150

#define APA_SM_SURROUND           1
#define APA_SM_LOGARITHMIC        2
#define APA_SM_FULLSUBSAMPLED     3

#define APA_SIM_CCF                11
#define APA_SIM_NCCF               12
#define APA_SIM_AMDF               13
#define APA_SIM_SSE                14

/*
********************************************************************************
*                         DEFINITION OF DATA TYPES
********************************************************************************
*/

struct apa_state_t;
typedef struct apa_state_t apa_state_t;
/*! handle for APA */
typedef struct apa_state_t* PCMDSP_APA_HANDLE;


/*
********************************************************************************
*                         DECLARATION OF PROTOTYPES
********************************************************************************
*/

/*! Allocates memory for state struct and initializes elements.
 *  @return 0 on success, 1 on failure */
uint8_t apa_init(apa_state_t **s);

/*! Sets state variables to initial value. */
void apa_reset(apa_state_t *s);

/*! Sets the audio configuration.
 *  Must be called once before processing can start.
 *  If called again during processing it will reset the state struct!
 *  Typical sample rates: 8000, 16000, 22050, 44100. Must be in range [APA_MIN_RATE,APA_MAX_RATE].
 *  Will also set a number of other state variables that depend on the sampling rate.
 *  @param[in,out] ps state
 *  @param[in] rate sample rate [Hz]
 *  @param[in] num_channels number of channels
 *  @return 0 on success, 1 on failure */
bool_t apa_set_rate(
    apa_state_t *ps,
    uint16_t      rate,
    uint16_t      num_channels);

/*! Set scaling.
 *  The scale is given in % and will be valid until changed again.
 *  Must be in range [APA_MIN_SCALE,APA_MAX_SCALE].
 *  @return 0 on success, 1 on failure */
bool_t apa_set_scale(apa_state_t *s, uint16_t scale);

bool_t apa_set_complexity_options(
    apa_state_t *s,
    uint16_t     wss,
    uint16_t     css);

bool_t apa_set_quality(
    apa_state_t *s,
    float       quality,
    uint16_t     qualityred,
    uint16_t     qualityrise);

bool_t apa_exit(
    apa_state_t **s);

uint8_t apa_exec(
    apa_state_t   *s,
    const int16_t a_in[],
    uint16_t      l_in,
    uint16_t      maxScaling,
    int16_t       a_out[],
    uint16_t     *l_out);

#endif /* JBM_PCMDSP_APA_H */
