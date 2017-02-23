/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

/*! @file jbm_types.h Data types used for JBM. */

#ifndef JBM_TYPES_H
#define JBM_TYPES_H JBM_TYPES_H

/**************************
* internally used types
***************************/
#ifndef _WIN32
#include <inttypes.h> /* part of C99 */
#else

#ifndef _UINT8_T
#define _UINT8_T
typedef unsigned char      uint8_t;
#endif /* _UINT8_T */

#ifndef _INT8_T
#define _INT8_T
typedef signed char        int8_t;
#endif /* _INT8_T */

#ifndef _UINT16_T
#define _UINT16_T
typedef unsigned short     uint16_t;
#endif /* _UINT16_T */

#ifndef _INT16_T
#define _INT16_T
typedef signed short       int16_t;
#endif /* _INT16_T */

#ifndef _UINT32_T
#define _UINT32_T
typedef unsigned int       uint32_t;
#endif /* _UINT32_T */

#ifndef _INT32_T
#define _INT32_T
typedef signed int         int32_t;
#endif /* _INT32_T */

#ifndef _UINT64_T
#define _UINT64_T
typedef unsigned __int64   uint64_t;
#endif /* _UINT64_T */

#ifndef _INT64_T
#define _INT64_T
typedef signed __int64     int64_t;
#endif /* _INT64_T */

#endif
#ifndef _BOOL_T
#define _BOOL_T
typedef unsigned char      bool_t;
#endif /* _BOOL_T */

#define true  (1)
#define false (0)

typedef float              Float;

#include "typedef.h"

#endif /* JBM_TYPES_H */
