/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#ifndef __BASTYPES_H
#define __BASTYPES_H

typedef unsigned char  BYTE;
typedef unsigned short WORD;
#if defined(__alpha__) || defined(__alpha) || defined(__sgi)
typedef unsigned int   DWORD; /* long is 64 bits on these machines */
#else
typedef unsigned long  DWORD;
#endif

typedef int            BOOL;
typedef signed   int   INT;
typedef signed long    LONG;
typedef unsigned long  ULONG;
typedef unsigned int   UINT;
typedef float          FLOAT;
typedef double         DOUBLE;
typedef unsigned char  UCHAR;
typedef char           CHAR;

/* from uld_types.h: */
typedef short          SHORT;
typedef unsigned short USHORT;
typedef long int       LINT;
typedef unsigned long int ULINT;

#ifndef TRUE
#define TRUE  1
#endif

#ifndef FALSE
#define FALSE 0
#endif

#ifndef NULL
#define NULL 0
#endif

#define INVALID_HANDLE NULL

#endif
