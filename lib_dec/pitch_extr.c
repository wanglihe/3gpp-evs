/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include "cnst.h"
#include "prot.h"
#include "options.h"
#include "basop_util.h"




/*-----------------------------------------------------------------*
 * Pitch prediction for frame erasure using linear fitting         *
 *-----------------------------------------------------------------*/
void pitch_pred_linear_fit(
    const short nbLostCmpt,         /* i:   bfi counter                                  */
    const short last_good,          /* i:   last classification type                     */
    float *old_pitch_buf,      /* i:   pitch lag buffer                             */
    float *old_fpitch,         /* i:                                                */
    float *T0_out,             /* o:   estimated close loop pitch                   */
    int pit_min,            /* i:   Minimum pitch lag                            */
    int pit_max,            /* i:   Maximum pitch lag                            */
    float *mem_pitch_gain,     /* i:   pitch gain [0] is the most recent subfr gain */
    int  limitation,
    short  plc_use_future_lag, /* i:   */
    short *extrapolationFailed,/* o:   flag if extrap decides not to change the pitch  */
    int  nb_subfr            /* i:   number of ACELP subframes                     */
)
{
    float T0 = 0;
    float mdy, dy[5], ftmp;
    short lcor = 5;
    short imax, i;
    float pg[8], ml[8]; /* local buffer for pitch gain and mem_lag*/
    short no_subfr_pred;
    float mem_lag[2*NB_SUBFR16k+2];


    if( nb_subfr == 4 )
    {
        for (i=0; i<2*NB_SUBFR+2; i++)
        {
            mem_lag[i] = old_pitch_buf[2*NB_SUBFR+1 - i];
        }
    }
    else  /* L_frame == L_FRAME16k */
    {
        for (i=0; i<2*NB_SUBFR16k+2; i++)
        {
            mem_lag[i] = old_pitch_buf[2*NB_SUBFR16k+1 - i];
        }
    }

    if( (int) *old_fpitch > pit_max )
    {
        *extrapolationFailed = 1;
        *T0_out = pit_max;
        printf("\n WARNING: (int)*old_fpitch > pit_max : old_fpitch = %f, pit_max = %i \n\n",*old_fpitch, pit_max );
        return;
    }

    if (nbLostCmpt == 1 && last_good >= UNVOICED_TRANSITION && last_good < ONSET)
    {
        if (plc_use_future_lag)
        {
            no_subfr_pred = 2;
        }
        else
        {
            no_subfr_pred = 4;
        }

        /* copy to local buffers, depending on availability of info about future subframes */
        mvr2r(mem_pitch_gain+no_subfr_pred-2,pg,8);
        mvr2r(mem_lag+no_subfr_pred-2,ml,8);

        mdy = 0.0f;
        for (i = (lcor-1); i >= 0; i--)
        {
            dy[i] = (ml[i] - ml[i+1]);
            mdy += dy[i];
        }

        /*---------------------------------------------------*
         * remove maximum variation
         *---------------------------------------------------*/
        ftmp = (float)fabs(dy[0]);
        imax = 0;
        for (i = 1; i <lcor; i++)
        {
            if (ftmp < (float)fabs(dy[i]))
            {
                ftmp = (float)fabs(dy[i]);
                imax = i;
            }
        }
        if ((fabs(dy[imax]) < 0.15f* *old_fpitch) && ((limitation == 1) || (fabs(dy[imax]) < fabs(mdy))))
        {
            Word16 pg_fx[5];
            Word32 ml_fx[5];
            Word32 pit, a, b, pita, pitb;
            Word16 sum0;
            Word16 const timeWeight[5] = {20480/*1.25f Q14*/, 18432/*1.125f Q14*/, 16384/*1.f Q14*/, 14336/*0.875f Q14*/, 12288/*.75f Q14*/}; /*Q14*/
            Word16 a1, a2, a3, a4, a5, tmpa, tmpb, b1, b2, b3, b4, b5;
            Word16 a_e, b_e, sum0_q;

            /* convert gains and lags to fixed precision */
            for( i=0; i<lcor; i++ )
            {
                pg_fx[i] = (short)(pg[i] / pow(2.f, -15 + 1)); /* Q14 */
                ml_fx[i] = (long)(ml[i] / pow(2.f, -31 + 15)); /* Q16 */
            }

            FOR( i=0; i<lcor; i++ )
            {
                pg_fx[i] = mult(mult(pg_fx[i] , pg_fx[i]) , timeWeight[i]); /*Q12 'til pg[lcor-1], Q14 'til pg[8]*/  move16();
            }

            /* Linear prediction (estimation) of pitch */
            /* sum0=(pg[1]+4*pg[2]+9*pg[3]+16*pg[4])*pg[0]+(pg[2]+4*pg[3]+9*pg[4])*pg[1]+(pg[3]+4*pg[4])*pg[2]+pg[4]*pg[3];*/
            {
                Word32 t1, t2, t3, t4, t5, t6, t7;
                Word16 e1, e2, e3, e4, e5, e6, e7;
                t1 = L_mult0(pg_fx[4], pg_fx[3]); /*Q24*/ /* t1 = pg[4]*pg[3] */
                e1 = 7;
                t2 = L_add(L_deposit_l(pg_fx[3]), L_shl(L_deposit_l(pg_fx[4]), 2)); /*Q12*/
                e2 = norm_l(t2);
                t2 = L_shl(t2, e2); /*Q12,-e2*/
                t2 = Mpy_32_16_1(t2, pg_fx[2]); /*Q9,-e2*/ /* t2 = (pg[3]+4*pg[4])*pg[2] */
                e2 = sub(22, e2);
                t3 = L_add(L_deposit_l(pg_fx[2]), L_add(L_shl(L_deposit_l(pg_fx[3]), 2), L_add(L_shl(L_deposit_l(pg_fx[4]), 3), L_deposit_l(pg_fx[4])))); /*Q12*/
                e3 = norm_l(t3);
                t3 = L_shl(t3, e3); /*Q12,-e3*/
                t3 = Mpy_32_16_1(t3, pg_fx[1]); /*Q9,-e3*/ /* t3 = (pg[2]+4*pg[3]+9*pg[4])*pg[1] */
                e3 = sub(22, e3);
                t4 = L_add(pg_fx[1], L_add(L_shl(L_deposit_l(pg_fx[2]), 2), L_add(L_add(L_shl(L_deposit_l(pg_fx[3]), 3), L_deposit_l(pg_fx[3])), L_shl(L_deposit_l(pg_fx[4]), 4)))); /*Q12*/
                e4 = norm_l(t4);
                t4 = L_shl(t4, e4); /*Q12,-e4*/
                t4 = Mpy_32_16_1(t4, pg_fx[0]); /*Q9,-e4*/ /* t4 = (pg[1]+4*pg[2]+9*pg[3]+16*pg[4])*pg[0] */
                e4 = sub(22, e4);
                t5 = BASOP_Util_Add_Mant32Exp(t1, e1, t2, e2, &e5);
                t6 = BASOP_Util_Add_Mant32Exp(t3, e3, t4, e4, &e6);
                t7 = BASOP_Util_Add_Mant32Exp(t5, e5, t6, e6, &e7); /*Q31,e7*/
                sum0_q = norm_l(t7);
                sum0 = round_fx(L_shl(t7, sum0_q)); /*Q15,e7-sum0_q*/
                sum0_q = add(15, sub(sum0_q, e7)); /* sum0 is now Qsum0_q*/
            }

            pit = 0;
            move16();
            IF( sum0 != 0 )
            {
                /* Shift to the right, changing Q as long as no precision is lost */
                WHILE( s_and(sum0, 1) == 0 )
                {
                    sum0 = shr(sum0, 1);
                    sum0_q = sub(sum0_q, 1);
                }

                /* float:
                a=-(
                  (  3*pg[1]+4*pg[2]+3*pg[3])*pg[0]                                           *//*a1*//*
      *ml[0]   +(
     ( 2*pg[2]+2*pg[3])*pg[1]-4*pg[1]*pg[0]                                   *//*a2*//*
      )*ml[1]  +(
     - 8*pg[2]*pg[0]-3*pg[2]*pg[1]+pg[3]*pg[2]                                *//*a3*//*
      )*ml[2]  +(
     -12*pg[3]*pg[0]-6*pg[3]*pg[1]-2*pg[3]*pg[2]                              *//*a4*//*
      )*ml[3]  +(
     -16*pg[4]*pg[0]   -9*pg[4]*pg[1]  -4*pg[4]*pg[2]   -pg[4]*pg[3]          *//*a5*//*
      )*ml[4] ) /sum0;
*/

                /*magic numbers: Q11 if not DIRECTLY marked otherwise*/
                a5 = mac_r(L_mac(L_mac(L_mult(mult_r(-32768,pg_fx[0])  /*Q8*/,pg_fx[4])/*Q5+16*/,  mult_r(-9*2048,pg_fx[1])/*Q8*/ , pg_fx[4]/*Q12*/ )/*Q5+16*/ ,mult_r(-4*2048,pg_fx[2])/*Q8*/, pg_fx[4]/*Q12*/)/*Q5+16*/,mult_r(pg_fx[4],-4096/*Q12->Q9*/),mult_r(pg_fx[3],16384/*Q12->Q11*/))/*Q5*/;
                a4 = mac_r(L_mac(L_mult(      mult_r(-12*2048,pg_fx[0])/*Q8*/,pg_fx[3] /*Q12*/)/*Q5+16*/,mult_r(-6*2048,pg_fx[1])/*Q8*/,pg_fx[3]/*Q12*/)/*Q5+16*/,mult_r(-2*2048,pg_fx[2])/*Q8*/,pg_fx[3]/*Q12*/)/*Q5*/;
                a3 = mac_r(L_mac(L_mult(      mult_r(-8*2048,pg_fx[0]) /*Q8*/,pg_fx[2]),mult_r(-3*2048,pg_fx[1])/*Q8*/,pg_fx[2]),mult_r(pg_fx[2],4096/*Q12->Q9*/),mult_r(pg_fx[3],16384/*12->Q11*/));/*Q5*/
                a2 = mac_r(L_mac(L_mult(      mult_r(2*2048,pg_fx[1])  /*Q8*/,pg_fx[2])/*Q5+16*/,mult_r(2*2048,pg_fx[1])/*Q8*/,pg_fx[3])/*Q5+16*/,mult_r(-4*2048,pg_fx[0])/*Q8*/,pg_fx[1]/*Q12*/)/*Q5*/;
                a1 = mac_r(L_mac(L_mult(      mult_r(3*2048,pg_fx[0])  /*Q8*/,pg_fx[1])/*Q5+16*/,mult_r(4*2048,pg_fx[0])/*Q8*/,pg_fx[2]/*Q12*/)/*Q5+16*/,mult_r(3*2048,pg_fx[0])/*Q8*/,pg_fx[3]/*Q12*/)/*Q5*/;

                a = L_mac(L_mac(L_mac(L_mac(L_mult(a1
                                                   , round_fx(L_shl(ml_fx[0],4)))/*Q4*/
                                            , round_fx(L_shl(ml_fx[1],4)) /*Q4*/, a2)
                                      , round_fx(L_shl(ml_fx[2],4)) /*Q4*/, a3)
                                , round_fx(L_shl(ml_fx[3],4)) /*Q4*/, a4)
                          , round_fx(L_shl(ml_fx[4],4)) /*Q4*/, a5);    /*Q-6+16 = Q10*/

                a_e = norm_l(a);
                a = L_shl(a, a_e);

                a1 =  BASOP_Util_Divide3216_Scale( L_negate(a),     /* Numerator  */ /*scalefactor 21*/
                                                   sum0,            /* Denominator*/  /*scalefactor 10*/
                                                   &tmpa);          /* scalefactor for result */

                /* Float:
                b=((   pg[1]+2*pg[2]+3*pg[3]+4*pg[4])*pg[0]                        *//*b1*//*
        *ml[0]   +
       (( pg[2]+2*pg[3]+3*pg[4])*pg[1]-pg[1]*pg[0])                 *//*b2*//*
       *ml[1]   +
       ( -2*pg[2]*pg[0]-pg[2]*pg[1]+(pg[3]+2*pg[4])*pg[2])          *//*b3*//*
     *ml[2]   +
     ( -3*pg[3]*pg[0]-2*pg[3]*pg[1]-pg[3]*pg[2]+pg[4]*pg[3])      *//*b4*//*
     *ml[3]   +
     ( -4*pg[4]*pg[0]-3*pg[4]*pg[1]-2*pg[4]*pg[2]-pg[4]*pg[3])    *//*b5*//*
     *ml[4]   )/sum0;                                                                                */

                /*magic numbers in Q13 if not DIRECTLY marked otherwise*/
                b1 = mac_r(L_mac(L_mac(L_mult(mult_r(pg_fx[1],pg_fx[0]),32768/4)/*Q7+16*/,mult_r(2*8192,pg_fx[0])/*Q10*/,pg_fx[2]/*Q12*/)/*Q7+16*/,mult_r(3*8192,pg_fx[0])/*Q10*/,pg_fx[3]/*Q12*/)/*Q7+16*/,  /*mult_r(4*8192,pg_fx[0])*/   pg_fx[0]/*Q10*/,pg_fx[4]/*Q12*/)/*Q7*/;
                b2 = mac_r(L_mac(L_mac(L_mult(mult_r(pg_fx[2],pg_fx[1]),32768/4)/*Q7+16*/,mult_r(2*8192,pg_fx[1]),pg_fx[3]),mult_r(3*8192,pg_fx[1]),pg_fx[4])/*Q7+16*/,mult_r(pg_fx[1],-32768/2/*Q12->Q12*/),mult_r(pg_fx[0],32768/2/*Q12->Q10*/))/*Q7*/;
                b3 = mac_r(L_mac(L_mac(L_mult(mult_r(-2*8192,pg_fx[0]),pg_fx[2])/*Q7+16*/,mult_r(pg_fx[2],-32768/2),mult_r(pg_fx[1],32768/2)),mult_r(pg_fx[3],32768/2),mult_r(pg_fx[2],32768/2))/*Q5+16*/,mult_r(2*8192,pg_fx[2]),pg_fx[4])/*Q7*/;
                b4 = mac_r(L_mac(L_mac(L_mult(mult_r(-3*8192,pg_fx[0]),pg_fx[3]),mult_r(-2*8192,pg_fx[1]),pg_fx[3]),mult_r(-32768/2,pg_fx[3]),mult_r(32768/2,pg_fx[2])),mult_r(32768/2,pg_fx[4]),mult_r(32768/2,pg_fx[3]));/*Q7*/
                b5 = mac_r(L_mac(L_mac(L_mult(mult_r(-32768/*(-4*8192)*/,pg_fx[0]),pg_fx[4]),mult_r(-3*8192,pg_fx[1]),pg_fx[4]),mult_r(-2*8192,pg_fx[2]),pg_fx[4]),mult_r(-32768/2,pg_fx[4]),mult_r(32768/2,pg_fx[3]))/*Q7*/;

                b = L_mac(L_mac(L_mac(L_mac(L_mult(b1
                                                   , round_fx(L_shl(ml_fx[0],4)))/*Q4*/
                                            , round_fx(L_shl(ml_fx[1],4)) /*Q4*/, b2)
                                      , round_fx(L_shl(ml_fx[2],4)) /*Q4*/, b3)
                                , round_fx(L_shl(ml_fx[3],4)) /*Q4*/, b4)
                          , round_fx(L_shl(ml_fx[4],4)) /*Q4*/, b5);    /*Q-4+16 = Q12*/
                /*predict pitch for 4th future subframe*/

                b_e = norm_l(b);
                b = L_shl(b, b_e);

                b1 =  BASOP_Util_Divide3216_Scale(b,       /* Numerator  */ /*scalefactor 19*/
                                                  sum0,                                /* Denominator*/  /*scalefactor 10*/
                                                  &tmpb);                              /* scalefactor for result*/

                /*pit = a + b * ((float)no_subfr_pred + (float)nb_subfr);*/
                pita = L_shl(  L_deposit_l(a1),add(add(sum0_q, 16-10+1),sub(tmpa, a_e)))/*Q16*/;
                pitb = L_shl_r(L_mult(b1/*Q15*/,add(no_subfr_pred,nb_subfr)/*Q0*/ ),add(add(sum0_q, 16-12),sub(tmpb, b_e)));
                pit = L_add(  pita , pitb ); /*Q16*/

                /* convert pitch back to float precision */
                T0 = (float)pit*pow(2.f, -31+15);

                /*limit pitch to allowed range*/
                if( T0 > pit_max )
                {
                    T0 = (float)pit_max;
                }

                if( T0 < pit_min )
                {
                    T0 = (float)pit_min;
                }

                *extrapolationFailed = 0;
            }
            else
            {
                T0 = 0;
                *extrapolationFailed = 1;
            }
        }
        else
        {
            T0 = 0;
            *extrapolationFailed = 1;
        }
    }
    else
    {
        T0 = *old_fpitch;
        *extrapolationFailed = 1;
    }

    *T0_out = T0;

    return;
}

void get_subframe_pitch(
    int nSubframes,         /* i:   number of subframes                             */
    float  pitchStart,         /* i:   starting pitch lag (in subframe -1)             */
    float  pitchEnd,           /* i:   ending pitch lag (in subframe nSubframes-1)     */
    float *pitchBuf            /* o:   interpolated pitch lag per subframe             */
)
{
    int i;
    float pitchDelta;

    assert((nSubframes > 0) && (pitchBuf != NULL) && (pitchStart >= 0) && (pitchEnd > 0));

    pitchDelta = (pitchEnd - pitchStart)/nSubframes;
    pitchBuf[0] = pitchStart+pitchDelta;
    for (i = 1; i < nSubframes; i++)
    {
        pitchBuf[i] = pitchBuf[i-1] + pitchDelta;
    }

    return;
}
