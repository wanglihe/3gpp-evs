/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"

/*-------------------------------------------------------------------*
 * rc_get_bits2()
 *
 *  Get number of bits needed to finalize range coder
 *-------------------------------------------------------------------*/

short rc_get_bits2(                 /* o: Number of bits needed         */
    const short N,                  /* i: Number of bits currently used */
    const unsigned int range        /* i: Range of range coder          */
)
{
    short tmp;
    tmp =   N + 2 + norm_ul(range);  /* aligned to BASOP */
    return tmp;
}


void rangeCoderFinalizationFBits(
    short Brc,
    unsigned int INTrc,
    short *FBits
)
{
    unsigned int   Bq15, UL_tmp;
    unsigned short Bq15ui16, B, E;
    short k;

    B =  30 - norm_ul(INTrc);  /* aligned to BASOP */


#define RCF_INIT_SHIFTp1        (RCF_INIT_SHIFT+1)
#define RCF_FINALIZE_LIMIT      ((1L << 16) - 1)
    *FBits = (Brc + 32) * 8;

    Bq15 = 0;
    k = B - RCF_INIT_SHIFT;
    if (k >= 0 )
    {
        Bq15 = INTrc >> ( k ); /* single op */
    }

    E = 2;
    for( k = 1; k < 4; k++)
    {
        Bq15ui16 = Bq15 >> (E & 1);
        UL_tmp   = (unsigned int)Bq15ui16;
        Bq15     = (UL_tmp * Bq15ui16) >> RCF_INIT_SHIFTp1 ;
        E = 2*B;
        if( Bq15 > RCF_FINALIZE_LIMIT )
        {
            E++;
        }
        B = E;
    }
    *FBits -= B;
#undef RCF_INIT_SHIFTp1
#undef RCF_FINALIZE_LIMIT

    return;
}

