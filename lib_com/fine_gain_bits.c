/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "rom_com.h"
#include "prot.h"


/*--------------------------------------------------------------------------
 * subband_gain_bits()
 *
 * HQ core encoder
 *--------------------------------------------------------------------------*/

void subband_gain_bits(
    const short *Rk,            /* i  : bit allocation per band (Q3)*/
    const short N,              /* i  : number of bands         */
    short *bits,          /* o  : gain bits per band      */
    const short *sfmsize        /* i  : Size of bands           */
)
{
    short i,b,tot;
    short bps;

    tot = 0;

    for ( i = 0; i < N; i++ )
    {
        bps = (short)(((int) Rk[i]*(int)fg_inv_tbl_fx[sfmsize[i]>>3])>>18) ;
        if (((sfmsize[i]*(bps+1)) << 3) - Rk[i] == 0)
        {
            /* correct approx. division result, to obtain exact integer division output */
            bps++;
        }
        bps = min(7, bps);



        b = fine_gain_bits[bps];

        bits[i] = b;
        tot += b;
    }

    if ( tot == 0)
    {
        /* If no gain bits were assigned, use one bit anyway for potential PVQ overage */
        bits[0] = 1;
    }

    return;
}

/*--------------------------------------------------------------------------*
 * assign_gain_bits()
 *
 * Assign gain adjustment bits and update bit budget
 *--------------------------------------------------------------------------*/

short assign_gain_bits(              /* o  : Number of assigned gain bits      */
    const short core,                /* i  : HQ core                           */
    const short BANDS,               /* i  : Number of bands                   */
    const short *band_width,         /* i  : Sub band bandwidth                */
    short *Rk,                 /* i/o: Bit allocation/Adjusted bit alloc. (Q3)*/
    short *gain_bits_array,    /* o  : Assigned gain bits                */
    short *Rcalc               /* o  : Bit budget for shape quantizer (Q3)*/
)
{
    short subband_cnt;
    short gain_bits_tot;
    short i;

    /* Allocate gain bits for every subband used, based on bit rate and bandwidth */
    if( core == HQ_CORE )
    {
        subband_gain_bits(Rk, BANDS, gain_bits_array, band_width);
    }
    else
    {
        set_s( gain_bits_array, 0, BANDS );
    }

    /* Re-adjust bit budget for gain quantization */
    subband_cnt = 0;
    gain_bits_tot = 0;
    *Rcalc = 0.0f;
    for (i = 0; i < BANDS; i++)
    {
        if (Rk[i] > 0)
        {
            subband_cnt++;
            Rk[i] -= gain_bits_array[i] * 8;
            gain_bits_tot += gain_bits_array[i];
            *Rcalc += Rk[i];
        }
    }

    return gain_bits_tot;
}
