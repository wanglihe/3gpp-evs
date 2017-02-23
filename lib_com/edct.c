/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"

static float const * get_edct_table(short length)
{
    float const * edct_table = NULL;
    switch (length)
    {
    case 1200:
        edct_table = edct_table_600;
        break;
    case 960:
        edct_table = edct_table_480;
        break;
    case 640:
        edct_table = edct_table_320;
        break;
    case 320:
        edct_table = edct_table_160;
        break;
    case 256:
        edct_table = edct_table_128;
        break;
    case 240:
        edct_table = edct_table_120;
        break;
    case 200:
        edct_table = edct_table_100;
        break;
    case 160:
        edct_table = edct_table_80;
        break;
    case 40:
        edct_table = edct_table_20;
        break;
    case 800:
        edct_table = edct_table_400;
        break;
    case 512:
        edct_table = edct_table_256;
        break;
    case 480:
        edct_table = edct_table_240;
        break;
    case 400:
        edct_table = edct_table_200;
        break;
    case 128:
        edct_table = edct_table_64;
        break;
    case 80:
        edct_table = edct_table_40;
        break;
    default:
        fprintf(stderr, "edct/edst(): length is not in table!");
        exit(-1);
        break;
    }
    return edct_table;
}

/*-----------------------------------------------------------------*
 * edct()
 *
 * DCT transform
 *-----------------------------------------------------------------*/

void edct(
    const float *x,         /* i  : input signal        */
    float *y,         /* o  : output transform    */
    short length      /* i  : length              */
)
{
    short i;
    float re[L_FRAME_PLUS/2];
    float im[L_FRAME_PLUS/2];
    const float *edct_table = 0;
    float local;

    edct_table = get_edct_table(length);

    /* Twiddling and Pre-rotate */
    for (i = 0; i < length/2; i++)
    {
        re[i] = x[2*i] * edct_table[i] + x[length-1-2*i] * edct_table[length/2-1-i];
        im[i] = x[length-1-2*i] * edct_table[i] - x[2*i] * edct_table[length/2-1-i];
    }

    DoFFT(re, im, length/2);

    local = (0.75f * EVS_PI)/ length;

    for (i = 0; i < length/2; i++)
    {
        float tmp;
        tmp = re[i] - im[i] * local;  /* 3*pi/(4*length) is a constant */
        im[i] = im[i] + re[i] * local;
        re[i] = tmp;
    }

    /* Post-rotate and obtain the output data */
    for (i = 0; i < length/2; i++)
    {
        y[2*i] = re[i] * edct_table[i] + im[i] * edct_table[length/2-1-i];
        y[length-1-2*i] = re[i] * edct_table[length/2-1-i] - im[i] * edct_table[i];
    }

    return;
}

/*-------------------------------------------------------------------------*
 * FUNCTION : edst()
 *
 * PURPOSE : DST_IV transform
 *-------------------------------------------------------------------------*/
void edst(
    const float *x,         /* i  : input signal        */
    float *y,         /* o  : output transform    */
    short length      /* i  : length              */
)
{
    short i;
    float re[L_FRAME_PLUS/2];
    float im[L_FRAME_PLUS/2];
    const float *edct_table = 0;
    float local;

    edct_table = get_edct_table(length);

    /* Twiddling and Pre-rotate */
    for (i = 0; i < length/2; i++)
    {
        re[i] = x[length-1-2*i] * edct_table[i] + x[2*i] *edct_table[length/2-1-i];
        im[i] = x[2*i] * edct_table[i] - x[length-1-2*i] * edct_table[length/2-1-i];
    }

    DoFFT(re, im, length/2);

    local = (0.75f * EVS_PI)/ length;

    for (i = 0; i < length/2; i++)
    {
        float tmp;
        tmp = re[i] - im[i] * local;  /* 3*pi/(4*length) is a constant */
        im[i] = im[i] + re[i] * local;
        re[i] = tmp;
    }

    /* Post-rotate and obtain the output data */
    for (i = 0; i < length/2; i++)
    {
        y[2*i] = re[i] * edct_table[i] + im[i] * edct_table[length/2-1-i];
        y[length-1-2*i] = im[i] * edct_table[i] - re[i] * edct_table[length/2-1-i];
    }

    return;
}

/*-----------------------------------------------------------------*
 * iedct_short()
 *
 * Inverse EDCT for short frames
 *-----------------------------------------------------------------*/

void iedct_short(
    const float *in,                /* i  : input vector     */
    float *out,               /* o  : output vector    */
    const short segment_length      /* i  : length           */
)
{
    float alias[MAX_SEGMENT_LENGTH];
    short i;

    edct(in, alias, segment_length/2);

    for (i = 0; i < segment_length/4; i++)
    {
        out[i] = alias[segment_length/4 + i];
        out[segment_length/4 + i] = -alias[segment_length/2 - 1 - i];
        out[segment_length/2 + i] = -alias[segment_length/4 - 1 - i];
        out[3*segment_length/4 + i] = -alias[i];
    }

    return;
}
