/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_dec.h"
#include "prot.h"


/*-------------------------------------------------------------------------
 * FEC_synchro_exc()
 *
 * Perform resynchronisation of the last glottal pulse in voiced frame lost
 *------------------------------------------------------------------------*/

static short FEC_synchro_exc(     /* o  : do_WI flag                                       */
    const short L_frame,          /* i  : length of the frame                               */
    float *exc,                   /* i/o: exc vector to modify                              */
    const short desire_puls_pos,  /* i  : Pulse position send by the encoder                */
    const short true_puls_pos,    /* i  : Present pulse location                            */
    const short Old_pitch         /* i  : Pitch use to create temporary adaptive codebook   */
)
{
    float exc_tmp[L_FRAME16k+L_SUBFR], min_energy, *pt_exc,*pt_exc1, ftmp, fact;
    short i, j, point_to_remove, point_to_add = -1, nb_min;
    short min_pos[L_FRAME16k/PIT_MIN_DOUBLEEXTEND], points_by_pos[L_FRAME16k/PIT_MIN_DOUBLEEXTEND], total_point, tmp_len;
    short *pt_pos, pos, start_search, tmp16;
    short remaining_len;

    /* Init */
    for( i=0; i<L_FRAME16k/PIT_MIN_DOUBLEEXTEND; i++ )
    {
        min_pos[i] = 10000;
        points_by_pos[i] = 0;
    }

    /* Find number of point to remove and number of minimum */
    point_to_remove = true_puls_pos - desire_puls_pos;  /* if it is negative it means remove point else it means add point */

    pos = L_frame - true_puls_pos;
    nb_min = pos/Old_pitch;

    /* if Old pitch < 128, must have at least 2 min */
    if (Old_pitch <= 128 && nb_min < 2)
    {
        nb_min = 2;
    }

    /* Must have at least 1 min */
    if (nb_min == 0)
    {
        /* Safety check, should be very rare but still needed*/
        nb_min = 1;
    }

    pt_exc = exc + pos;

    /* Find start_search for minimum energy search*/
    start_search = -3*Old_pitch/4;

    if (start_search  + pos < 0)
    {
        /* This case is rare but still need to be taken care*/
        start_search   = -pos;

        if(abs(start_search) < Old_pitch/8)
        {
            /* it's not safe to remove/add point inside 1/8 of the pulse position */
            return 0;
        }
    }

    /* Find min energy in the first pitch section */
    min_energy = 65536*65536.0f;

    /* --------------------------------------------------------------------
     * The minimum energy regions are determined by the computing the energy
     * using a sliding 5-sample window. The minimum energy position is set
     * at the middle of the window at which the energy is at minimum
     * --------------------------------------------------------------------*/

    ftmp = (pt_exc[start_search]*pt_exc[start_search]);
    ftmp += (pt_exc[start_search+1]* pt_exc[start_search+1]);
    ftmp += (pt_exc[start_search+2]* pt_exc[start_search+2]);
    ftmp += (pt_exc[start_search+3]* pt_exc[start_search+3]);
    ftmp += (pt_exc[start_search+4]* pt_exc[start_search+4]);

    if (ftmp < min_energy)
    {
        min_energy = ftmp;
        min_pos[0] = (pos + start_search +2);
    }

    for (i = start_search; i < -5; i++)
    {
        ftmp  -= pt_exc[i]*pt_exc[i];
        ftmp  += pt_exc[i+5]*pt_exc[i+5];
        if (ftmp < min_energy)
        {
            min_energy = ftmp;
            min_pos[0] = pos + i + 2;
        }
    }

    for (j = 1; j < nb_min; j++)
    {
        min_pos[j] = min_pos[j-1]-Old_pitch;
        /* If the first minimum is in the past, forget this minimum */
        if (min_pos[j] < 0)
        {
            /* Safety check */
            min_pos[j] = -10000;
            nb_min = nb_min-1;
        }
    }

    /* safety-measure against not properly initialized min_pos[] */
    if( min_energy >= 65536*65536.0f )
    {
        return 0;
    }

    /*--------------------------------------------------------------------
     * Determine the number of samples to be added or removed at each pitch
     * cycle whereby less samples are added/removed at the beginning and
     * more towards the end of the frame
     * --------------------------------------------------------------------*/

    if(nb_min == 1 || abs(point_to_remove) == 1)
    {
        nb_min = 1;
        points_by_pos[0] = (short)abs(point_to_remove);
    }
    else
    {
        /* First position */
        fact = (float)abs(point_to_remove) / (nb_min*nb_min);

        total_point = (short)(fact+0.5);
        points_by_pos[0] = total_point;

        for (i = 2; i <= nb_min; i++)
        {
            points_by_pos[i-1] = (short)(fact*(i*i) - total_point+0.5) ;
            total_point += points_by_pos[i-1];

            /* ensure a constant increase */
            if (points_by_pos[i-1] < points_by_pos[i-2])
            {
                tmp16 = points_by_pos[i-2];
                points_by_pos[i-2] = points_by_pos[i-1];
                points_by_pos[i-1] = tmp16;
            }
        }
    }

    /* --------------------------------------------------------------------
     * Sample deletion or insertion is performed in minimum energy regions.
     * At the end of this section the last maximum pulse in the concealed
     * excitation is forced to align to the actual maximum pulse position
     * at the end of the frame which is transmitted in the future frame.
     * --------------------------------------------------------------------*/

    if ( point_to_remove > 0 )
    {
        point_to_add  =point_to_remove;
    }

    pt_exc = exc_tmp;
    pt_exc1 = exc;

    i = 0;
    pt_pos = min_pos + nb_min - 1;
    if ( point_to_add > 0 ) /* add some points */
    {
        remaining_len = L_frame;

        for (i = 0; i < nb_min; i++)
        {
            /* Copy section */
            if(i == 0)
            {
                /* Compute len to copy */
                tmp_len = *pt_pos;
            }
            else
            {
                /* Compute len to copy */
                tmp_len = *pt_pos - *(pt_pos+1) - points_by_pos[i-1];
            }

            mvr2r( pt_exc1, pt_exc, tmp_len);
            remaining_len -= tmp_len;
            pt_exc1 += tmp_len;
            pt_exc  += tmp_len;

            /* add some points */
            ftmp = -(*pt_exc1/20);
            for (j=0; j <points_by_pos[i];  j++)
            {
                *pt_exc++ = ftmp;
                ftmp = -ftmp;
            }

            remaining_len -= points_by_pos[i];
            pt_pos--;
        }

        /* Copy remaining length */
        mvr2r( pt_exc1, pt_exc, remaining_len);

        /* Update true excitation buffer*/
        mvr2r( exc_tmp, exc, L_frame );
    }
    else /* Remove points */
    {
        remaining_len = L_frame;

        for(i = 0; i < nb_min; i++)
        {
            if (i == 0)
            {
                /* Compute len to copy */
                tmp_len = *pt_pos;
            }
            else
            {
                /* Compute len to copy */
                tmp_len = *pt_pos - *(pt_pos+1) - points_by_pos[i-1];
            }

            mvr2r( pt_exc1, pt_exc, tmp_len);
            remaining_len -= tmp_len;
            pt_exc1+= tmp_len;
            pt_exc += tmp_len;

            /* Remove points */
            for (j=0; j <points_by_pos[i];  j++)
            {
                pt_exc1++ ;
            }
            pt_pos--;
        }

        /* Copy remaining length */
        mvr2r( pt_exc1, pt_exc, remaining_len);

        /* Update true excitation buffer*/
        mvr2r( exc_tmp, exc, L_frame );
    }

    return 1;
}


/*---------------------------------------------------------------------*
 * FEC_SinOnset()
 *
 * Create an artificial onset when it is lost
 *---------------------------------------------------------------------*/

void FEC_SinOnset (
    float *exc,        /* i/o : exc vector to modify                                           */
    short puls_pos,    /* i   : last pulse position desired                                    */
    short T0,
    float enr_q,       /* i   : energy provide by the encoder                                  */
    float *Aq,         /* i   : A(z) filter                                                    */
    const short L_frame      /* i   : frame length                                                   */
)
{
    short P0, onset_len, sign = 0, i, len, L_subfr;
    float h1[L_SUBFR16k], mem[M], exc_tmp[L_FRAME16k + MODE1_L_FIR_FER];
    float *pt_end, *pt_exc, enr_LP, gain;

    L_subfr = L_SUBFR;
    if( L_frame == L_FRAME16k )
    {
        L_subfr = L_SUBFR16k;
    }

    if (T0 < 2*L_subfr)
    {
        onset_len = 2*L_subfr;
    }
    else
    {
        onset_len = T0;
    }


    P0 = puls_pos;

    if (P0 < 0)
    {
        sign = 1;
        P0 = -P0;
    }

    if ( P0 > PIT_MAX && L_frame == L_FRAME )
    {
        P0 = PIT_MAX; /* Should never be the case, however... */
    }
    else if ( P0 > PIT16k_MAX && L_frame == L_FRAME16k )
    {
        P0 = PIT16k_MAX; /* Should never be the case, however... */
    }
    set_f( exc_tmp, 0, L_frame+ MODE1_L_FIR_FER);                      /* Reset excitation vector */

    /*--------------------------------------------*
     * Find LP filter impulse response energy
     *--------------------------------------------*/

    set_f( h1, 0, L_subfr );                          /* Find the impulse response */
    set_f( mem, 0, M );
    h1[0] = 1.0f;
    syn_filt(Aq, M, h1, h1, L_subfr, mem, 0);
    enr_LP = dotp( h1, h1, L_subfr ) + 0.01f;  /* Find the impulse response energy */

    /*------------------------------------------------------------------------------------------*
     * Construct the harmonic part as a train of low-pass filtered pulses
     *------------------------------------------------------------------------------------------*/

    pt_exc = exc_tmp + L_frame- 1 - MODE1_L_FIR_FER/2 - P0;              /* beginning of the 1st pulse */
    pt_end = exc_tmp + onset_len;

    len = (short)(pt_exc - pt_end);
    if( len > MODE1_L_FIR_FER )
    {
        len  = MODE1_L_FIR_FER;
    }
    if(!sign)
    {
        for( i=0; i< len; i++ )
        {
            /* The filter response would have E=1 in full band. */
            pt_exc[i] += h_low[i];      /* As it is lp filtered, the E is somewhat lower */
        }
    }
    else
    {
        for( i=0; i< len; i++ )
        {
            /* The filter response would have E=1 in full band. */
            pt_exc[i] -= h_low[i];      /* As it is lp filtered, the E is somewhat lower */
        }
    }

    gain = (float)sqrt(1.5*enr_q / enr_LP);    /* divide by LP filter E, scale by transmitted E */
    gain *= 0.96f;
    for( i=0; i< L_frame; i++ )
    {
        exc_tmp[i] *= gain;
    }

    mvr2r(&exc_tmp[L_frame - L_EXC_MEM], exc, L_EXC_MEM);

    return;
}

short FEC_enhACB(
    const short L_frame,                   /* i   : frame length                                                              */
    const short last_L_frame,              /* i   : frame length of last frame                                                */
    float *exc_io,                   /* i/o : adaptive codebook memory                                                  */
    const short new_pit,                   /* i   : decoded first frame pitch                                                 */
    const short puls_pos,                  /* i   : decoder position of the last glottal pulses decoded in the previous frame */
    const float bfi_pitch                  /* i   : pitch used for concealment                                                */
)
{
    short Tc, P0, sign, pit_search;
    short Tlist[10], Terr, diff_pit, dist_Plast;
    float ftmp;
    float exc[L_FRAME16k + L_SUBFR];
    short do_WI = 1;

    set_f(exc, 0.0f, L_FRAME16k - L_EXC_MEM);
    set_f(exc+L_FRAME16k, 0.0f, L_SUBFR);
    mvr2r(exc_io, exc + L_FRAME16k - L_EXC_MEM, L_EXC_MEM);

    Tc = (short)bfi_pitch;
    mvr2r(exc + L_FRAME16k - Tc, exc + L_FRAME16k, L_SUBFR);

    /*------------------------------------------------------------
     * Decode phase information transmitted in the bitstream
     * (The position of the absolute maximum glottal pulse from
     * the end of the frame and its sign)
     *------------------------------------------------------------*/

    P0 = puls_pos;
    sign = 0;
    if (P0 < 0)
    {
        sign = 1;
        P0 = -P0;
    }

    if( L_frame == L_FRAME )
    {
        if( P0>PIT_MAX )
        {
            P0 = PIT_MAX; /* Should never be the case, however... */
        }
    }
    else  /* L_frame == L_FRAME16k */
    {
        if( P0>PIT16k_MAX )
        {
            P0 = PIT16k_MAX; /* Should never be the case, however... */
        }
    }

    /*----------------------------------------------------------------------------------
     * Find the position of the first the maximum(minimum) lp_filtered pulse
     * <----- Mem --->|<--------------------- L_frame ------------>|<----- L_SUBFR --->|
     *                |<-------pit_search---->                     |                   |
     *----------------------------------------------------------------------------------*/

    pit_search = Tc;

    Tlist[0] = findpulse( L_frame, exc+L_frame-pit_search , pit_search,DEC, &sign);

    Terr = (short) abs(pit_search-Tlist[0]-P0);

    dist_Plast = Tc-Tlist[0];

    Tlist[1] = findpulse( L_frame, exc+L_frame-pit_search , pit_search+L_SUBFR,DEC, &sign);

    if(Terr > abs(Tlist[1]-Tc + P0))
    {
        dist_Plast = Tc-Tlist[1];
        Terr = (short)abs(Tlist[1]-Tc + P0);
    }

    diff_pit = (short)abs(new_pit - Tc);
    ftmp = (float) (int)((float)L_frame/(float)Tc+0.5);

    if ( Terr <= ftmp * diff_pit
            && Terr != 0       /* If Terr = 0, no resynchronization required */
            && Terr < L_SUBFR) /* prevent catastrophy search */
    {
        /* perform excitation resynchronization here */
        do_WI = FEC_synchro_exc( L_frame, exc, P0, dist_Plast, Tc );
        mvr2r( exc + L_FRAME16k - L_EXC_MEM, exc_io, L_EXC_MEM );
    }
    else
    {
        do_WI = 0;
    }

    if( last_L_frame != L_FRAME16k )
    {
        /* This is essential in case of bit rate switching + FEC*/
        do_WI = 0;
    }

    return do_WI;
}

