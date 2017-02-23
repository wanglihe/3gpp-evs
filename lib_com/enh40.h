/*
  ===========================================================================
   File: ENH40.H                                         v.2.3 - 30.Nov.2009
  ===========================================================================

            ITU-T  STL  BASIC OPERATORS

            40-BIT ARITHMETIC OPERATORS

   History:
   07 Nov 04   v2.0     Incorporation of new 32-bit / 40-bit / control
                        operators for the ITU-T Standard Tool Library as
                        described in Geneva, 20-30 January 2004 WP 3/16 Q10/16
                        TD 11 document and subsequent discussions on the
                        wp3audio@yahoogroups.com email reflector.
   March 06   v2.1      Changed to improve portability.

  ============================================================================
*/


#ifndef _ENH40_H
#define _ENH40_H


#include "stl.h"


#ifdef _MSC_VER
#define MAX_40 (0x0000007fffffffff)
#define MIN_40 (0xffffff8000000000)
#endif /* ifdef _MSC_VER */



#define L40_OVERFLOW_OCCURED(  L40_var1) (Overflow = 1, exit(1), L40_var1)
#define L40_UNDERFLOW_OCCURED( L40_var1) (Overflow = 1, exit(2), L40_var1)



/*****************************************************************************
*
*  Prototypes for enhanced 40 bit arithmetic operators
*
*****************************************************************************/
Word40 L40_shr(   Word40 L40_var1, Word16 var2);
Word40 L40_shr_r( Word40 L40_var1, Word16 var2);
Word40 L40_shl(   Word40 L40_var1, Word16 var2);
Word40 L40_shl_r( Word40 L40_var1, Word16 var2);

static __inline Word40 L40_mult(      Word16 var1,   Word16 var2);

static __inline Word40 L40_mac(      Word40 L40_var1, Word16 var1,   Word16 var2);
static __inline Word16 mac_r40(      Word40 L40_var1, Word16 var1,   Word16 var2);

static __inline Word40 L40_msu(      Word40 L40_var1, Word16 var1,   Word16 var2);
static __inline Word16 msu_r40(      Word40 L40_var1, Word16 var1,   Word16 var2);


void Mpy_32_16_ss( Word32 L_var1, Word16 var2,   Word32 *L_varout_h, UWord16 *varout_l);
void Mpy_32_32_ss( Word32 L_var1, Word32 L_var2, Word32 *L_varout_h, UWord32 *L_varout_l);


Word40 L40_lshl(  Word40 L40_var1, Word16 var2);
Word40 L40_lshr(  Word40 L40_var1, Word16 var2);

static __inline Word40 L40_set(      Word40 L40_var1);
static __inline UWord16 Extract40_H( Word40 L40_var1);
static __inline UWord16 Extract40_L( Word40 L40_var1);
static __inline UWord32 L_Extract40(  Word40 L40_var1);

static __inline Word40 L40_deposit_h( Word16 var1);
static __inline Word40 L40_deposit_l( Word16 var1);
static __inline Word40 L40_deposit32( Word32 L_var1);

static __inline Word40 L40_round(     Word40 L40_var1);
static __inline Word16 round40(       Word40 L40_var1);


Word40 L40_add( Word40 L40_var1, Word40 L40_var2);
Word40 L40_sub( Word40 L40_var1, Word40 L40_var2);
Word40 L40_abs( Word40 L40_var1);
Word40 L40_negate( Word40 L40_var1);
Word40 L40_max( Word40 L40_var1, Word40 L40_var2);
Word40 L40_min( Word40 L40_var1, Word40 L40_var2);
Word32 L_saturate40( Word40 L40_var1);
Word16 norm_L40(     Word40 L40_var1);



/*#ifdef _MSC_VER*/
static __inline Word40 L40_set( Word40 L40_var1)
{
    Word40 L40_var_out;

#if defined(_MSC_VER) && (_MSC_VER <= 1200)
    L40_var_out =  L40_var1 & 0x000000ffffffffff;

    if( L40_var1 & 0x8000000000)
        L40_var_out = L40_var_out | 0xffffff0000000000;
#else
    L40_var_out =  L40_var1 & 0x000000ffffffffffLL;

    if( L40_var1 & 0x8000000000LL)
        L40_var_out = L40_var_out | 0xffffff0000000000LL;
#endif


    return( L40_var_out);
}
/*#endif*/ /* ifdef _MSC_VER */



static __inline UWord16 Extract40_H( Word40 L40_var1)
{
    UWord16 var_out;

    var_out = ( UWord16)( L40_var1 >> 16);


    return( var_out);
}


static __inline UWord16 Extract40_L( Word40 L40_var1)
{
    UWord16 var_out;

    var_out = ( UWord16)( L40_var1);


    return( var_out);
}


static __inline UWord32 L_Extract40( Word40 L40_var1)
{
    UWord32 L_var_out;

    L_var_out = ( UWord32) L40_var1;


    return(L_var_out);
}


static __inline Word40 L40_deposit_h( Word16 var1)
{
    Word40 L40_var_out;

    L40_var_out = (( Word40) var1) << 16;

#if defined(_MSC_VER) && (_MSC_VER <= 1200)
    if( var1 & 0x8000)
    {
        L40_var_out = L40_set( L40_var_out | 0xff00000000);
#else
    if( var1 & 0x8000)
    {
        L40_var_out = L40_set( L40_var_out | 0xff00000000LL);
#endif
    }


    return( L40_var_out);
}


static __inline Word40 L40_deposit_l( Word16 var1)
{
    Word40 L40_var_out;

    L40_var_out = var1;

#if defined(_MSC_VER) && (_MSC_VER <= 1200)
    if( var1 & 0x8000)
    {
        L40_var_out = L40_set( L40_var_out | 0xffffff0000);
#else
    if( var1 & 0x8000)
    {
        L40_var_out = L40_set( L40_var_out | 0xffffff0000LL);
#endif
    }


    return( L40_var_out);
}


static __inline Word40 L40_deposit32( Word32 L_var1)
{
    Word40 L40_var_out;

    L40_var_out = ( Word40) L_var1;

#if defined(_MSC_VER) && (_MSC_VER <= 1200)
    if( L_var1 & 0x80000000)
    {
        L40_var_out = L40_set( L40_var_out | 0xff00000000);
#else
    if( L_var1 & 0x80000000)
    {
        L40_var_out = L40_set( L40_var_out | 0xff00000000LL);
#endif
    }


    return( L40_var_out);
}








static __inline Word40 L40_round( Word40 L40_var1)
{
    Word40 L40_var_out;
    Word40 L40_constant;

#if defined(_MSC_VER) && (_MSC_VER <= 1200)
    L40_constant = L40_set( 0xffffff0000);
#else
    L40_constant = L40_set( 0xffffff0000LL);
#endif

    L40_var_out = L40_add( 0x8000, L40_var1);
    L40_var_out = L40_var_out & L40_constant;


    return( L40_var_out);
}


static __inline Word16 round40( Word40 L40_var1)
{
    Word16 var_out;

    var_out = extract_h( L_saturate40( L40_round( L40_var1)));


    return( var_out);
}


static __inline Word40 L40_mult( Word16 var1, Word16 var2)
{
    Word32 L_var_out;
    Word40 L40_var_out;

    L_var_out = ( Word32) var1 * ( Word32) var2;
    L40_var_out = ( Word40) L_var_out;

    /* Below line can not overflow, so we can use << instead of L40_shl. */
    L40_var_out = L40_var_out << 1;


    return( L40_var_out);
}












static __inline Word40 L40_mac( Word40 L40_var1, Word16 var2, Word16 var3)
{
    Word40 L40_var_out;

    L40_var_out = L40_mult( var2, var3);
    L40_var_out = L40_add( L40_var1, L40_var_out);


    return( L40_var_out);
}






static __inline Word16 mac_r40( Word40 L40_var1, Word16 var2, Word16 var3)
{
    Word40 L40_var_out;
    Word16 var_out;

    L40_var_out = L40_mac( L40_var1, var2, var3);
    var_out = round40( L40_var_out);


    return( var_out);
}






static __inline Word40 L40_msu( Word40 L40_var1, Word16 var2, Word16 var3)
{
    Word40 L40_var_out;

    L40_var_out = L40_mult( var2, var3);
    L40_var_out = L40_sub( L40_var1, L40_var_out);


    return( L40_var_out);
}






static __inline Word16 msu_r40( Word40 L40_var1, Word16 var2, Word16 var3)
{
    Word40 L40_var_out;
    Word16 var_out;

    L40_var_out = L40_msu( L40_var1, var2, var3);
    var_out = round40( L40_var_out);


    return( var_out);
}










































#endif /*_ENH40_H*/


/* end of file */


