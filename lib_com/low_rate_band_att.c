/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <math.h>
#include "options.h"
#include "rom_com.h"
#include "prot.h"


/*--------------------------------------------------------------------------*
 * fine_gain_pred()
 *
 * Fine gain prediction
 *--------------------------------------------------------------------------*/

void fine_gain_pred(
    const short *sfm_start,                /* i  : Sub band start indices    */
    const short *sfm_end,                  /* i  : Sub band end indices      */
    const short *sfm_size,                 /* i  : Sub band bandwidths       */
    const short *i_sort,                   /* i  : Energy sorting indices    */
    const short *K,                        /* i  : Number of pulses per band */
    const short *maxpulse,                 /* i  : Maximum pulse per band    */
    const short *R,                        /* i  : Bits per sub band (Q3)    */
    const short num_sfm,                   /* i  : Number of sub bands       */
    float *xq,                       /* i/o: Quantized vector /quantized vector with finegain adj */
    short *y,                        /* i/o: Quantized vector (int)    */
    float *fg_pred,                  /* o  : Predicted fine gains      */
    const short core                       /* i  : Core                      */
)
{
    short i, band;
    float gp;
    float xx;
    float accuracy;
    short k, bw;
    float att;

    for( band = 0; band < num_sfm; band++)
    {

        k  = K[i_sort[band]];
        if(k > 0)
        {
            bw = sfm_size[i_sort[band]];
            xx = 0;
            for(i = sfm_start[i_sort[band]]; i < sfm_end[i_sort[band]]; i++)
            {
                xx += xq[i] * xq[i];
            }

            if ( xx > 0)
            {
                /* Normalize synthesis to RMS=1.0 */
                gp = (float) sqrt(bw / xx);

                if (core == HQ_CORE && R != NULL && R[i_sort[band]] <= 256)
                {
                    accuracy = ((float)k/(float)bw)*maxpulse[i_sort[band]];
                    att = 1.0f - 0.05f / accuracy;
                    att = max( 0.840896f, att); /* Limit attenuation to norm quantizer error, 2^-0.25 */
                    gp *= att;
                }

                fg_pred[band] = gp;
            }
            else
            {
                fg_pred[band] = 0;
            }
        }
        else
        {
            fg_pred[band] = 0;
            for(i = sfm_start[i_sort[band]]; i < sfm_end[i_sort[band]]; i++)
            {
                y[i] = 0;
                xq[i] = 0;
            }
        }
    }
    return;
}

/*--------------------------------------------------------------------------*
 * get_max_pulses()
 *
 * Find the maximum pulse height (in unit pulses) in each band
 *--------------------------------------------------------------------------*/

void get_max_pulses(
    const short *band_start,               /* i  : Sub band start indices    */
    const short *band_end,                 /* i  : Sub band end indices      */
    const short *k_sort,                   /* i  : Indices for sorting by energy */
    const short *npulses,                  /* i  : Pulses per sub band       */
    const short  BANDS,                    /* i  : Number of bands           */
    short *inp_vector,               /* i/o: Encoded shape vectors (int)*/
    short *maxpulse                  /* o  : Maximum pulse height per band */
)
{
    short i, k;
    int npul;
    int maxp;

    for (k = 0; k < BANDS; k++)
    {
        npul = npulses[k_sort[k]];
        maxp = 0;
        if (npul > 0)
        {
            for (i = band_start[k_sort[k]]; i < band_end[k_sort[k]]; i++)
            {
                if (abs(inp_vector[i]) > maxp)
                {
                    maxp = abs(inp_vector[i]);
                }
            }
        }
        maxpulse[k_sort[k]] = maxp;
    }

    return;
}

/*--------------------------------------------------------------------------*
 * fine_gain_dec()
 *
 * Fine gain decoder. Decodes fine gain adjustments and applies correction to
 * predicted fine gains
 *--------------------------------------------------------------------------*/

void fine_gain_dec
(
    Decoder_State *st,          /* i/o: Decoder state struct                 */
    const short *ord,           /* i  : Indices for energy order             */
    const short num_sfm,        /* i  : Number of bands                      */
    const short *gain_bits,     /* i  : Gain adjustment bits per sub band    */
    float *fg_pred        /* i/o: Predicted gains / Corrected gains    */
)
{
    short band;
    short gbits;
    short idx;
    float gain_dbq;

    for ( band = 0; band < num_sfm; band++)
    {
        gbits = gain_bits[ord[band]];
        if ( fg_pred[band] != 0 && gbits > 0 )
        {
            idx = (short)get_next_indice( st, (short)gbits );
            gain_dbq = finegain[gbits-1][idx];

            /* Update prediced gain with quantized correction */
            fg_pred[band] *= (float)pow(10, gain_dbq * 0.05f);
        }
    }

    return;
}


