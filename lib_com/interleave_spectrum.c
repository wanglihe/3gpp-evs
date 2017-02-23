/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"

/*--------------------------------------------------------------------------*
 * interleave_spectrum()
 *
 * Interleave the spectrum
 *--------------------------------------------------------------------------*/

void interleave_spectrum(
    float *coefs,     /* i/o: input and output coefficients                 */
    const short length      /* i  : length of spectrum                            */
)
{
    short i, j, k;
    float *p1, *p2, *p3, *p4;
    float *p_out;
    float coefs_out[STOP_BAND];
    short sublen;
    short grps;
    const short *bw;
    const short *cnt;

    if ( length == L_FRAME48k )
    {
        bw = intl_bw_48;
        cnt = intl_cnt_48;
        grps = N_INTL_GRP_48;
    }
    else if( length == L_FRAME32k )
    {
        bw = intl_bw_32;
        cnt = intl_cnt_32;
        grps = N_INTL_GRP_32;
    }
    else /* length == L_FRAME16k */
    {
        bw = intl_bw_16;
        cnt = intl_cnt_16;
        grps = N_INTL_GRP_16;
    }

    sublen = length/4;
    p1 = coefs;
    p2 = coefs + sublen;
    p3 = coefs + sublen*2;
    p4 = coefs + sublen*3;
    p_out = coefs_out;

    for (i = 0; i < grps; i++)
    {
        for (j = 0; j < cnt[i]; j++)
        {
            for (k = 0; k < bw[i]; k++)
            {
                *p_out++ = *p1++;
            }

            for (k = 0; k < bw[i]; k++)
            {
                *p_out++ = *p2++;
            }

            for (k = 0; k < bw[i]; k++)
            {
                *p_out++ = *p3++;
            }

            for (k = 0; k < bw[i]; k++)
            {
                *p_out++ = *p4++;
            }
        }
    }

    /* For FB the interleaved spectrum is 800 samples */
    mvr2r(coefs_out, coefs, (short)(p_out - coefs_out));

    return;
}


/*--------------------------------------------------------------------------*
 * de_interleave_spectrum()
 *
 * Deinterleave the spectrum
 *--------------------------------------------------------------------------*/

void de_interleave_spectrum(
    float *coefs,         /* i/o: input and output coefficients  */
    short length          /* i  : length of spectrum             */
)
{
    short   i, j, k;
    float   *p1, *p2, *p3, *p4;
    float   *p_in;
    float   coefs_out[L_FRAME48k];
    short   sublen;
    short   grps;
    const short   *bw;
    const short   *cnt;

    if ( length == L_FRAME48k )
    {
        bw = intl_bw_48;
        cnt = intl_cnt_48;
        grps = N_INTL_GRP_48;
    }
    else if( length == L_FRAME32k )
    {
        bw = intl_bw_32;
        cnt = intl_cnt_32;
        grps = N_INTL_GRP_32;
    }
    else /* length == L_FRAME16k */
    {
        bw = intl_bw_16;
        cnt = intl_cnt_16;
        grps = N_INTL_GRP_16;
    }

    set_f(coefs_out, 0, L_FRAME48k);
    sublen = length/4;
    p1 = coefs_out;
    p2 = coefs_out + sublen;
    p3 = coefs_out + sublen*2;
    p4 = coefs_out + sublen*3;
    p_in = coefs;

    for (i = 0; i < grps; i++)
    {
        for (j = 0; j < cnt[i]; j++)
        {
            for (k = 0; k < bw[i]; k++)
            {
                *p1++ = *p_in++;
            }
            for (k = 0; k < bw[i]; k++)
            {
                *p2++ = *p_in++;
            }
            for (k = 0; k < bw[i]; k++)
            {
                *p3++ = *p_in++;
            }
            for (k = 0; k < bw[i]; k++)
            {
                *p4++ = *p_in++;
            }
        }
    }

    mvr2r(coefs_out, coefs, length);

    return;
}
