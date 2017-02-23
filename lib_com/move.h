/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#ifndef _MOVE_H
#define _MOVE_H

/* BASOP -> FLC brigde: data move counting */

#include "stl.h"

static __inline void move16( void)
{
}

static __inline void move32( void)
{
}

static __inline void test( void)
{
}

static __inline void logic16( void)
{
}

static __inline void logic32( void)
{
}


/*-------- legacy ----------*/
#define data_move()          move16()
#define L_data_move()        move32()
#define data_move_external() move16()
#define compare_zero()       test()
/*-------- end legacy ----------*/

#define cast16 move16

#endif /* _MOVE_H */
