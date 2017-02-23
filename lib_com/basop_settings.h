/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#ifndef __BASOP_SETTINGS_H
#define __BASOP_SETTINGS_H

#include "stl.h"
#include "basop_mpy.h"

#define _LONG                long
#define _SHORT               short
#ifdef _WIN32
#define _INT64               __int64
#else
#define _INT64               long long
#endif

#define WORD32_BITS         32
#define MAXVAL_WORD32       ((signed)0x7FFFFFFF)
#define MINVAL_WORD32       ((signed)0x80000000)
#define WORD32_FIX_SCALE    ((_INT64)(1)<<(WORD32_BITS-1))

#define WORD16_BITS         16
#define MAXVAL_WORD16       (((signed)0x7FFFFFFF)>>16)
#define MINVAL_WORD16       (((signed)0x80000000)>>16)
#define WORD16_FIX_SCALE    ((_INT64)(1)<<(WORD16_BITS-1))

/*!
  \def  Macro converts a float < 1 to Word32 fixed point with saturation and rounding
*/
#define FL2WORD32(val)                                                                                                     \
(Word32)( ( (val) >= 0) ?                                                                                                               \
((( (double)(val) * (WORD32_FIX_SCALE) + 0.5 ) >= (double)(MAXVAL_WORD32) ) ? (_LONG)(MAXVAL_WORD32) : (_LONG)( (double)(val) * (double)(WORD32_FIX_SCALE) + 0.5)) : \
((( (double)(val) * (WORD32_FIX_SCALE) - 0.5) <=  (double)(MINVAL_WORD32) ) ? (_LONG)(MINVAL_WORD32) : (_LONG)( (double)(val) * (double)(WORD32_FIX_SCALE) - 0.5)) )

/*!
  \def   Macro converts a float < 1 to Word16 fixed point with saturation and rounding
*/
#define FL2WORD16(val)                                                                                                     \
(Word16)( ( (val) >= 0) ?                                                                                                               \
((( (double)(val) * (WORD16_FIX_SCALE) + 0.5 ) >= (double)(MAXVAL_WORD16) ) ? (_LONG)(MAXVAL_WORD16) : (_LONG)( (double)(val) * (double)(WORD16_FIX_SCALE) + 0.5)) : \
((( (double)(val) * (WORD16_FIX_SCALE) - 0.5) <=  (double)(MINVAL_WORD16) ) ? (_LONG)(MINVAL_WORD16) : (_LONG)( (double)(val) * (double)(WORD16_FIX_SCALE) - 0.5)) )

/*!
  \def   Macro converts a Word32 fixed point to Word16 fixed point <1 with saturation
*/
#define WORD322WORD16(val)                                                                           \
 ( ( ((((val) >> (WORD32_BITS-WORD16_BITS-1)) + 1) > (((_LONG)1<<WORD16_BITS)-1)) && ((_LONG)(val) > 0) ) ? \
 (Word16)(_SHORT)(((_LONG)1<<(WORD16_BITS-1))-1):(Word16)(_SHORT)((((val) >> (WORD32_BITS-WORD16_BITS-1)) + 1) >> 1) )

/*!
  \def   Macro converts a Word32 fixed point < 1 to float shifts result left by scale
*/
#define WORD322FL_SCALE(x,scale) ( ((float)((_LONG)(x))) / (((_INT64)1<<(WORD32_BITS-1 - (scale)))) )

/*!
  \def   Macro converts a float < 1 to Word32 fixed point with saturation and rounding, shifts result right by scale
*/
/* Note: Both x and scale must be constants at compile time, scale must be in range -31..31 */
#define FL2WORD32_SCALE(x,scale) FL2WORD32((double)(x) *(((_INT64)1<<(WORD32_BITS-1 - (scale)))) / ((_INT64)1<<(WORD32_BITS-1)))

/*!
  \def   Macro converts a Word16 fixed point < 1 to float shifts result left by scale
*/
#define WORD162FL_SCALE(x,scale) ( ((float)((_LONG)(x))) / (((_INT64)1<<(WORD16_BITS-1 - (scale)))) )

/*!
  \def   Macro converts a float < 1 to Word16 fixed point with saturation and rounding, shifts result right by scale
*/
/* Note: At compile time, x must be a float constant and scale must be an integer constant in range -15..15 */
#define FL2WORD16_SCALE(x,scale) FL2WORD16((float)(x) *(((_INT64)1<<(WORD16_BITS-1 - (scale)))) / ((_INT64)1<<(WORD16_BITS-1)))


/* Word16 Packed Type */
typedef struct
{
    struct
    {
        Word16 re;
        Word16 im;
    } v;
} PWord16;

#endif /* __BASOP_SETTINGS_H */
