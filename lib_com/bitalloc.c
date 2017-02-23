/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"
#include "basop_util.h"
#include "basop_proto_func.h"


/*--------------------------------------------------------------------------
 * bitalloc()
 *
 * Adaptive bit allocation for 20kHz audio codec
 *--------------------------------------------------------------------------*/

void bitalloc (
    short *y,                /* i  : reordered norm of sub-vectors                 */
    short *idx,              /* i  : reordered sub-vector indices                  */
    short sum,               /* i  : number of available bits                      */
    short N,                 /* i  : number of norms                               */
    short K,                 /* i  : maximum number of bits per dimension          */
    short *r,                /* o  : bit-allacation vector                         */
    const short *sfmsize,          /* i  : band length                                   */
    const short hqswb_clas         /* i  : signal classification flag                    */
)
{
    short i, j, k, n, m, v, im;
    short diff, temp;
    short fac;
    short ii;
    short SFM_thr = SFM_G1G2;

    N -= 1;

    if ( hqswb_clas == HQ_HARMONIC )
    {
        SFM_thr = 22;
    }

    fac = 3;
    K -= 2;
    im = 1;
    diff = sum;
    n = sum >> 3;
    for ( i=0; i<n; i++ )
    {
        k = 0;
        temp = y[0];
        for (m=1; m<im; m++)
        {
            if ( temp < y[m] )
            {
                temp = y[m];
                k = m;
            }
        }

        if ( temp < y[m] )
        {
            k = m;
            if ( im < N)
            {
                im++;
            }
        }

        j = idx[k];
        if ( sum >= sfmsize[j] && r[j] < K )
        {
            y[k] -= fac;
            r[j]++;
            if ( r[j] >= K )
            {
                y[k] = MIN16B;
            }
            sum -= sfmsize[j];
        }
        else
        {
            y[k] = MIN16B;
            k++;
            if ( k == im && im < N )
            {
                im++;
            }
        }

        if ( sum < WID_G1 || diff == sum )
        {
            break;
        }

        diff = sum;
        v = N - 1;

        if (k > v)
        {
            for ( ii=0; ii<=N; ii++ )
            {
                if ( y[ii] > MIN16B )
                {
                    if( ii < N )
                    {
                        im = ii + 1;
                    }

                    break;
                }
            }
        }
    }

    if ( sum >= WID_G2)
    {
        for (i=0; i<=N; i++)
        {
            j = idx[i];
            if ( j >= SFM_G1 && j < SFM_thr && r[j] == 0 )
            {
                r[j] = 1;
                sum -= WID_G2;
                if ( sum < WID_G2 )
                {
                    break;
                }
            }
        }
    }

    if ( sum >= WID_G2 )
    {
        for (i=0; i<=N; i++)
        {
            j = idx[i];
            if ( j >= SFM_G1 && j < SFM_thr && r[j] == 1 )
            {
                r[j] = 2;
                sum -= WID_G2;
                if ( sum < WID_G2 )
                {
                    break;
                }
            }
        }
    }

    if ( sum >= WID_G1 )
    {
        for (i=0; i<=N; i++)
        {
            j = idx[i];
            if ( j < SFM_G1 && r[j] == 0 )
            {
                r[j] = 1;
                sum -= WID_G1;
                if ( sum < WID_G1 )
                {
                    break;
                }
            }
        }
    }

    if ( sum >= WID_G1 )
    {
        for (i=0; i<=N; i++)
        {
            j = idx[i];
            if ( j < SFM_G1 && r[j] == 1 )
            {
                r[j] = 2;
                sum -= WID_G1;
                if ( sum < WID_G1 )
                {
                    break;
                }
            }
        }
    }

    return;
}


/*-------------------------------------------------------------------*
 * BitAllocF()
 *
 * Fractional bit allocation
 *-------------------------------------------------------------------*/

short BitAllocF (
    short *y,                /* i  : norm of sub-vectors                           */
    long  bit_rate,          /* i  : bitrate                                       */
    short B,                 /* i  : number of available bits                      */
    short N,                 /* i  : number of sub-vectors                         */
    short *R,                /* o  : bit-allocation indicator                      */
    short *Rsubband,         /* o  : sub-band bit-allocation vector (Q3)           */
    const short hqswb_clas,  /* i  : hq swb class                                  */
    const short num_env_bands   /* i  : Number sub bands to be encoded for HQ_SWB_BWE  */
)
{
    Word16 fac;
    Word16 i, n, Nmin, Bits, bs, low_rate = 0;

    Word16 m_fx;
    Word32 t_fx, B_fx;
    Word32 L_tmp1, L_tmp2;
    Word16 tmp, exp1, exp2;
    Word32 Rsubband_w32_fx[NB_SFM];                                             /* Q15  */
    Word16 B_w16_fx;

    set_i( Rsubband_w32_fx, 0, NB_SFM);

    fac = 3;
    if (L_sub(bit_rate, 32000) < 0)
    {
        bs = 2;
    }
    else
    {
        bs = 3;
    }
    low_rate = 1;

    Nmin = N;
    if ( sub(Nmin,SFM_N) > 0)
    {
        Nmin = SFM_N;
    }

    /* Initial bits distribution */
    if (sub(hqswb_clas , HQ_GEN_SWB) == 0 || sub(hqswb_clas , HQ_GEN_FB) == 0)
    {
        /* Initial bits distribution */
        L_tmp1 = 0;
        m_fx = 0;
        for ( i = 0; i < num_env_bands ; i++)
        {
            L_tmp1 = L_mac0(L_tmp1, Nb[i], y[i]);
        }
        L_tmp1 = L_msu0(L_tmp1, fac, B);

        t_fx = 0;
        n = 0;
        tmp = add(band_end[num_env_bands-1], shl(band_end[num_env_bands-1], 1));
        exp1 = norm_s(tmp);
        tmp = div_s(16384, shl(tmp, exp1));/*15 + 14 - exp1*/
        exp2 = norm_s(tmp);
        tmp = shl(tmp, exp2);
        exp1 = add(29, sub(exp2, exp1));

        for ( i = 0; i < N; i++)
        {
            L_tmp2 = L_sub(L_mult0(y[i], band_end[num_env_bands-1]), L_tmp1);
            Rsubband_w32_fx[i] = L_mult0(extract_l(L_tmp2), Nb[i]);
            move32();/*Q0*/
            if (Rsubband_w32_fx[i] > 0)
            {
                n = add(n,Nb[i]);
                Rsubband_w32_fx[i] = Mpy_32_16(Rsubband_w32_fx[i], tmp);
                move32();/*exp1 - 15*/
                Rsubband_w32_fx[i] = L_shl(Rsubband_w32_fx[i], sub(30, exp1));/*Q15*/

                t_fx = L_add(t_fx, Rsubband_w32_fx[i]);/*Q0*/
            }
            else
            {
                Rsubband_w32_fx[i] = 0;
                move32();
            }
        }
    }
    else
    {
        /* Initial bits distribution */
        L_tmp1 = 0;
        m_fx = 0;
        for ( i = 0; i < N ; i++)
        {
            L_tmp1 = L_mac0(L_tmp1, Nb[i], y[i]);
        }
        L_tmp1 = L_msu0(L_tmp1, fac, B);


        t_fx = 0;
        n = 0;
        tmp = add(band_end[N-1], shl(band_end[N-1], 1));
        exp1 = norm_s(tmp);
        tmp = div_s(16384, shl(tmp, exp1));/*15 + 14 - exp1*/
        exp2 = norm_s(tmp);
        tmp = shl(tmp, exp2);
        exp1 = add(29, sub(exp2, exp1));
        for ( i = 0; i < N; i++)
        {
            L_tmp2 = L_sub(L_mult0(y[i], band_end[N-1]), L_tmp1);
            Rsubband_w32_fx[i] = L_mult0(extract_l(L_tmp2), Nb[i]);
            move32();/*Q0*/
            if (Rsubband_w32_fx[i] > 0)
            {
                n = add(n,Nb[i]);
                Rsubband_w32_fx[i] = Mpy_32_16(Rsubband_w32_fx[i], tmp);
                move32();/*exp1 - 15*/
                Rsubband_w32_fx[i] = L_shl(Rsubband_w32_fx[i], sub(30, exp1));/*Q15*/

                t_fx = L_add(t_fx, Rsubband_w32_fx[i]);/*Q0*/
            }
            else
            {
                Rsubband_w32_fx[i] = 0;
                move32();
            }
        }
    }

    /* Distribute the remaining bits to subbands with non-zero bits */
    B_fx = L_shl(B, 15);
    WHILE (L_sub(L_shr(L_add(t_fx, 16384), 15) , B) != 0)
    {
        L_tmp1 = L_sub(t_fx, B_fx);
        exp1 = sub(norm_l(L_tmp1), 1);
        exp2 = norm_s(n);
        tmp = div_s(extract_h(L_shl(L_tmp1, exp1)), shl(n, exp2));/*15 + 15 + exp1 - 16 - exp2*/
        m_fx = shl(tmp, sub(exp2, exp1));/*Q14*/

        t_fx = 0;
        n = 0;
        for ( i = 0; i < N; i++)
        {
            if (Rsubband_w32_fx[i] > 0)
            {
                Rsubband_w32_fx[i] = L_msu(Rsubband_w32_fx[i], m_fx, Nb[i]);
                move32();

                if (Rsubband_w32_fx[i] > 0)
                {
                    n = add(n,Nb[i]);

                    t_fx = L_add(t_fx, Rsubband_w32_fx[i]);
                }
                else
                {
                    Rsubband_w32_fx[i] = 0;
                    move32();
                }
            }
        }
    }
    Bits = B;

    /* Impose bit-constraints to subbands with less than minimum bits*/
    t_fx = 0;
    n = 0;
    for ( i = 0; i < N; i++)
    {
        if (Rsubband_w32_fx[i] > 0)
        {
            test();
            test();
            if ((L_sub(Rsubband_w32_fx[i] , L_shl(add(bs, LNb[i]), 15)) <0) && (sub(low_rate,1) == 0))
            {
                Rsubband_w32_fx[i] = 0;
                move32();
            }
            else if ( L_sub(Rsubband_w32_fx[i] , L_shl(Nb[i], 15)) <=0)
            {
                B = sub(B,Nb[i]);
                Rsubband_w32_fx[i] = L_shl(Nb[i], 15);
                move32();
            }
            else
            {
                n = add(n,Nb[i]);
                t_fx = L_add(t_fx, Rsubband_w32_fx[i]);
            }
        }
    }

    /* Distribute the remaining bits to subbands with more than 1-bit per sample */
    WHILE (L_sub(L_shr(L_add(t_fx, 16384), 15) ,B) != 0)
    {
        L_tmp1 = L_sub(t_fx, L_shl(B, 15));
        L_tmp2 = L_abs(L_tmp1);

        if( n>0 )
        {
            exp1 = sub(norm_l(L_tmp2), 1);
            exp2 = norm_s(n);
            tmp = div_s(extract_h(L_shl(L_tmp2, exp1)), shl(n, exp2));/*15 + 15 + exp1 - 16 - exp2*/
            m_fx = shl(tmp, sub(exp2, exp1));/*Q14*/
            if (L_tmp1 < 0)
            {
                m_fx = negate(m_fx);
            }

            t_fx = 0;
            n = 0;
            for( i = 0; i < N; i++)
            {
                if (L_sub(Rsubband_w32_fx[i] , L_shl(Nb[i], 15)) > 0)
                {
                    Rsubband_w32_fx[i] = L_msu(Rsubband_w32_fx[i], m_fx, Nb[i]);
                    if (L_sub(Rsubband_w32_fx[i] ,L_shl(Nb[i], 15)) > 0)
                    {
                        n   = add(n,Nb[i]);

                        t_fx = L_add(t_fx, Rsubband_w32_fx[i]);
                    }
                    else
                    {
                        B = sub(B,Nb[i]);

                        Rsubband_w32_fx[i] = L_shl(Nb[i], 15);
                        move32();
                    }
                }
            }
        }
        /*In case no subband has enough bits more than 1-bit per sample, take bits off the higher subbands */
        if (t_fx == 0)
        {
            for ( i = N-1; i >= 0; i--)
            {
                if (Rsubband_w32_fx[i] > 0)
                {
                    B = add( B, Nb[i] );
                    Rsubband_w32_fx[i] = 0;
                    move32();
                    if ( B >= 0)
                    {
                        BREAK;
                    }
                }
            }
            BREAK;
        }
    }

    /* fine redistribution of over-allocated or under-allocated bits */
    tmp = 0;
    for ( i = 0; i < N; i++)
    {
        Rsubband[i] = extract_l(L_shr(Rsubband_w32_fx[i], 12));
        tmp = add(tmp, Rsubband[i]);
    }

    B = Bits;
    B_w16_fx = shl(B, 3);
    if (sub(tmp ,B_w16_fx)>0)
    {
        tmp = sub(tmp, B_w16_fx);
        for ( i = 0; i < N; i++)
        {
            if (sub(Rsubband[i], add(shl(Nb[i], 3), tmp)) >= 0)
            {
                Rsubband[i] = sub(Rsubband[i], tmp);
                BREAK;
            }
        }
    }
    else
    {
        tmp = sub(tmp, B_w16_fx);
        for ( i = 0; i < N; i++)
        {
            if (Rsubband[i] > 0)
            {
                Rsubband[i] = sub(Rsubband[i], tmp);
                BREAK;
            }
        }
    }

    /* Calculate total used bits and initialize R to be used for Noise Filling */
    tmp = 0;
    for ( i = 0; i < N; i++)
    {
        tmp = add(tmp, Rsubband[i]);
        R[i] = shr(Rsubband[i], 3);
    }

    return shr(tmp, 3);
}

/*-------------------------------------------------------------------*
 * Bit_group()
 *
 * bit allocation in group
 *-------------------------------------------------------------------*/
static void Bit_group_fx (
    Word16 *y,                 /* i  : norm of sub-band                              Q0*/
    Word16 start_band,         /* i  : start band indices                            Q0*/
    Word16 end_band,           /* i  : end band indices                              Q0*/
    Word16 Bits,               /* i  : number of allocation bits in group            Q0*/
    Word16 thr,                /* i  : smallest bit number for allocation in group   Q0*/
    Word32 *Rsubband_fx,          /* o  : bit allocation of sub-band                    Q21*/
    Word16 *fac_fx                /* i  : weight factor for norm of sub-band            Q13*/
)
{
    Word16 i, j, k, m, y_index[16], index[16], bit_band, band_num, norm_sum;
    Word16 tmp,exp;
    Word16 factor_fx;
    Word32 R_temp_fx[16], R_sum_fx = 0, R_sum_org_fx = 0, Bits_avg_fx = 0;
    Word32 L_tmp;
    UWord32 lo;

    /* initialization for bit allocation in one group*/
    tmp = 6554;
    move16();  /*Q15  1/5    */
    IF(sub(thr,5) == 0)
    {
        tmp = 6554;
        move16();  /*Q15  1/5    */
    }
    IF(sub(thr,6) == 0)
    {
        tmp = 5462;
        move16();/*Q15  1/6 */
    }
    IF(sub(thr,7) == 0)
    {
        tmp = 4682;
        move16();/*Q15  1/7 */
    }
    bit_band = mult(tmp, Bits); /*0+15-15=0, Q0 */
    band_num = sub(end_band,start_band);

    FOR( i = 0; i < band_num; i++ )
    {
        y_index[i] = y[add(i,start_band)];
        move16();
        index[i] = i;
        move16();
    }

    /* Rearrange norm vector in decreasing order */
    reordvct(y_index, band_num, index);
    /* norm vector modification */

    factor_fx = div_s(1, band_num);/*Q15 */
    IF ( sub(thr,5) > 0 )
    {
        FOR ( i = 0; i < band_num; i++ )
        {
            L_tmp = L_mult(i,factor_fx);/*Q16 */
            tmp = extract_h(L_shl(L_tmp, 13)); /*Q13 */
            tmp = sub(fac_fx[1],tmp);/*Q13 */
            L_tmp = L_mult(y_index[i],tmp);/*Q14 */
            y_index[i] = extract_h(L_shl(L_tmp, 2));/*Q0 */
        }
    }
    ELSE
    {
        FOR ( i = 0; i < band_num; i++ )
        {
            L_tmp = L_mult(i,factor_fx);/*Q16 */
            tmp = extract_h(L_shl(L_tmp, 13)); /*Q13 */
            tmp = sub(fac_fx[0],tmp);/*Q13 */
            L_tmp = L_mult(y_index[i],tmp);/*Q14 */
            y_index[i] = extract_h(L_shl(L_tmp, 2));/*Q0 */
        }
    }

    /* bit allocation based on modified norm */
    L_tmp = L_mult(band_num,24576);/*Q16 */
    tmp = extract_h(L_shl(L_tmp,7));/*Q7 */
    IF ( sub(shl(bit_band,7),tmp) >= 0 )
    {
        FOR ( j = 0; j < band_num; j++)
        {
            IF ( y_index[j] < 0 )
            {
                y_index[j] = 0;
                move16();
            }
            R_temp_fx[j] = 2097152;
            move16();/*Q21 = 1     move16(); */
        }

        i = sub(band_num,1);
        norm_sum = 0;/*Q0 */
        FOR (k = 0; k <= i; k++)
        {
            norm_sum = add(norm_sum,y_index[k]);
        }

        FOR (j = 0; j < band_num; j++)
        {
            IF(norm_sum == 0)
            {
                FOR (k = 0; k <= i; k++)
                {
                    R_temp_fx[k] = 0;
                    move32();/*Q21 */
                }
            }
            ELSE
            {
                exp = norm_s(norm_sum);
                tmp = shl(norm_sum, exp);/*Q(exp) */
                tmp = div_s(16384,tmp); /*Q(15+14-exp) */
                Bits_avg_fx = L_mult(tmp, Bits);/*Q(30-exp) */

                FOR (k = 0; k <= i; k++)
                {
                    L_tmp = L_shl(L_deposit_l(y_index[k]),24);
                    Mpy_32_32_ss(Bits_avg_fx,L_tmp,&L_tmp,&lo);

                    R_temp_fx[k] = L_shl(L_tmp,sub(exp,2));
                    move32();/*Q21 */
                }
            }

            L_tmp = L_shl(L_deposit_l(thr),21);/*Q21 */
            IF ( L_sub(R_temp_fx[i],L_tmp) < 0 )
            {
                R_temp_fx[i] = 0;
                move32();
                norm_sum = sub(norm_sum,y_index[i]);
                i--;
            }
            ELSE
            {
                BREAK;
            }
        }
    }
    ELSE
    {
        FOR ( j = 0; j < bit_band; j++ )
        {
            IF ( y_index[j] < 0 )
            {
                y_index[j] = 0;
                move16();
            }
            R_temp_fx[j] = 2097152;
            move32();/*Q21 = 1 */
        }

        FOR ( j = bit_band; j < band_num; j++ )
        {
            R_temp_fx[j] = 0;
            move32();
        }

        norm_sum = 0;
        FOR (k = 0; k < bit_band; k++)
        {
            norm_sum = add(norm_sum,y_index[k]);
        }

        i = bit_band;
        FOR (j = 0; j < bit_band; j++)
        {
            IF(norm_sum == 0)
            {
                FOR (k = 0; k < i; k++)
                {
                    R_temp_fx[k] = 0;
                    move32();/*Q21                    */
                }
            }
            ELSE
            {
                exp = norm_s(norm_sum);
                tmp = shl(norm_sum, exp);/*Q(exp) */
                tmp = div_s(16384,tmp); /*Q(15+14-exp) */
                Bits_avg_fx = L_mult(tmp, Bits);/*Q(30-exp) */
                FOR (k = 0; k < i; k++)
                {
                    L_tmp = L_shl(L_deposit_l(y_index[k]),24);
                    Mpy_32_32_ss(Bits_avg_fx,L_tmp,&L_tmp,&lo);
                    R_temp_fx[k] = L_shl(L_tmp,sub(exp,2));
                    move32();/*Q21 */
                }
            }
            R_sum_fx = 0;
            L_tmp = L_shl(L_deposit_l(thr),21);/*Q21 */
            FOR (k = 0; k < i; k++)
            {
                IF (L_sub(R_temp_fx[k],L_tmp) < 0)
                {
                    FOR(m = k; m < i; m++)
                    {
                        norm_sum = sub(norm_sum,y_index[m]);
                        R_temp_fx[m] = 0;
                        move32();/*Q21 */
                    }
                    i = k;
                    BREAK;
                }
                ELSE
                {
                    R_sum_fx = L_add(R_sum_fx,R_temp_fx[k]);
                }
            }
            IF (L_sub(R_sum_fx,R_sum_org_fx) == 0)
            {
                BREAK;
            }

            R_sum_org_fx = R_sum_fx;
        }
    }

    /*  index comeback */
    FOR ( k = 0 ; k < band_num; k++ )
    {
        j = index[k];
        move16();
        Rsubband_fx[add(j,start_band)] = R_temp_fx[k];
        move32();
    }

    return;

}

/*-------------------------------------------------------------------*
 * BitAllocWB()
 *
 * WB bit allocation
 *-------------------------------------------------------------------*/

short BitAllocWB (
    short *y,                /* i  : norm of sub-vectors                           */
    short B,                 /* i  : number of available bits                      */
    short N,                 /* i  : number of sub-vectors                         */
    short *R,                /* o  : bit-allocation indicator                      */
    short *Rsubband          /* o  : sub-band bit-allocation vector (Q3)           */
)
{
    Word16 t_fx;
    Word16 i, j, k, B1, B2, B3, B_saved;
    Word16 Rsum_fx, Rsum_sub_fx[3];
    Word32 Ravg_sub_32_fx[3], R_diff_32_fx[2];
    Word16 factor_fx[2];/*Q13 */
    Word16 BANDS;
    Word16 tmp,exp;
    Word32 L_tmp,L_tmp1;
    Word32 Rsubband_buf[NB_SFM];
    UWord16 lo;

    BANDS = N;
    move16();
    IF( sub(BANDS,SFM_N) > 0)
    {
        BANDS = SFM_N;
        move16();
    }
    /* Init Rsubband to non-zero values for bands to be allocated bits */
    FOR (k = 0; k < BANDS; k++)
    {
        Rsubband_buf[k] = 2097152;
        move32();/*Q21 */
    }
    /* Calculate the norm sum and average of sub-band */
    Rsum_sub_fx[0] = 0;
    FOR ( j = 0; j < SFM_G1; j++ )
    {
        IF ( y[j] > 0 )
        {
            Rsum_sub_fx[0] = add(Rsum_sub_fx[0],y[j]);
            move16();/*Q0 */
        }
    }
    Ravg_sub_32_fx[0] = L_mult(Rsum_sub_fx[0], 2048);
    move32();/*Q16  0+15+1   , q15 1/16 =2048 */

    Rsum_sub_fx[1] = 0;
    move32();
    FOR ( j = SFM_G1; j < SFM_G1G2; j++ )
    {
        IF ( y[j] > 0 )
        {
            Rsum_sub_fx[1] = add(Rsum_sub_fx[1],y[j]);
            move16();/*Q0 */
        }
    }
    Ravg_sub_32_fx[1] = L_mult(Rsum_sub_fx[1], 4096); /*16  0+15+1   , q15 1/8 =4096 */

    Rsum_sub_fx[2] = 0;
    move16();
    FOR ( j = SFM_G1G2; j < BANDS; j++ )
    {
        IF ( y[j] > 0 )
        {
            Rsum_sub_fx[2] = add(Rsum_sub_fx[2],y[j]);
            move16();/*Q0 */
        }
    }
    tmp = div_s(1, BANDS-SFM_G1G2); /*Q15 */
    Ravg_sub_32_fx[2] = L_mult(Rsum_sub_fx[2], tmp);
    move32();/*Q16 */

    /* Bit allocation for every group */
    tmp = add(Rsum_sub_fx[0],Rsum_sub_fx[1]);
    Rsum_fx = add(tmp,Rsum_sub_fx[2]);/*Q0     */

    factor_fx[0] = 16384;/*Q13     move16(); */
    factor_fx[1] = 24576;/*Q13     move16(); */
    {
        R_diff_32_fx[0] = L_sub(Ravg_sub_32_fx[0], Ravg_sub_32_fx[1]);
        move32();/*Q16 */
        R_diff_32_fx[1] = L_sub(Ravg_sub_32_fx[1], Ravg_sub_32_fx[2]);
        move32();/*Q16 */

        IF ( L_sub(R_diff_32_fx[0],393216) < 0 && L_sub(R_diff_32_fx[1],245760) < 0 )
        {
            IF(Rsum_fx == 0)
            {
                B1 = 0;
                move16();
                B2 = 0;
                move16();
                B3 = 0;
                move16();
            }
            ELSE
            {
                exp = norm_s(Rsum_fx);
                tmp = shl(Rsum_fx,exp);/*Q(exp) */
                tmp = div_s(16384,tmp);/*Q(15+14-exp) */
                L_tmp1 = L_mult(B,Rsum_sub_fx[0]);/*Q1 */
                Mpy_32_16_ss(L_tmp1,tmp,&L_tmp,&lo);
                B1 = extract_h(L_shl(L_tmp,add(exp,1)));/*Q0 */
                IF(L_sub(L_tmp1,L_mult(B1,Rsum_fx)) > 0 && L_sub(L_tmp1,L_mult(add(B1,1),Rsum_fx)) >= 0)
                {
                    B1 = add(B1,1);
                }
                L_tmp1 = L_mult(B,Rsum_sub_fx[1]);/*Q1 */
                Mpy_32_16_ss(L_tmp1,tmp,&L_tmp,&lo);
                B2 = extract_h(L_shl(L_tmp,add(exp,1)));/*Q0 */
                IF(L_sub(L_tmp1,L_mult(B2,Rsum_fx)) > 0 && L_sub(L_tmp1,L_mult(add(B2,1),Rsum_fx)) >= 0)
                {
                    B2 = add(B2,1);
                }
                L_tmp1 = L_mult(B,Rsum_sub_fx[2]);/*Q1 */
                Mpy_32_16_ss(L_tmp1,tmp,&L_tmp,&lo);
                B3 = extract_h(L_shl(L_tmp,add(exp,1)));/*Q0 */
                IF(L_sub(L_tmp1,L_mult(B3,Rsum_fx)) > 0 && L_sub(L_tmp1,L_mult(add(B3,1),Rsum_fx)) >= 0)
                {
                    B3 = add(B3,1);
                }
            }
            IF ( L_sub(Ravg_sub_32_fx[2],786432) > 0 )
            {
                B_saved = 0;
                move16();
                IF ( sub(B1,288) > 0 )
                {
                    B_saved = sub(B1,288);
                    B1 = 288;
                    move16();
                }

                IF ( sub(B2,256) > 0 )
                {
                    tmp = sub(B2,256);
                    B_saved = add(B_saved,tmp);
                    B2 = 256;
                    move16();
                }

                IF ( sub(B3,96) > 0 )
                {
                    tmp = sub(B3,96);
                    B_saved = add(B_saved,tmp);
                    B3 = 96;
                    move16();
                }

                IF ( B_saved > 0 )
                {
                    IF ( sub(B1,288) == 0 )
                    {
                        tmp = shr(B_saved,1);
                        B2 = add(B2,tmp);
                        tmp = sub(B,B1);
                        B3 = sub(tmp,B2);
                    }
                    ELSE
                    {
                        tmp = shr(B_saved,1);
                        B1 = add(B1,tmp);
                        IF ( sub(B2,256) == 0 )
                        {
                            tmp = sub(B,B1);
                            B3 = sub(tmp,B2);
                        }
                        ELSE
                        {
                            tmp = sub(B,B1);
                            B2 = sub(tmp,B3);
                        }
                    }
                }
            }

            factor_fx[0] = 16384;
            move16();/*Q13 */
            factor_fx[1] = 12288;
            move16();/*Q13 */
        }
        ELSE
        {
            IF(Rsum_fx == 0)
            {
                B1 = 0;
                move16();
                B2 = 0;
                move16();
                B3 = B;
                move16();
            }
            ELSE
            {
                exp = norm_s(Rsum_fx);
                tmp = shl(Rsum_fx,exp);/*Q(exp) */
                tmp = div_s(16384,tmp);/*Q(15+14-exp) */
                L_tmp1 = L_mult(B,Rsum_sub_fx[0]);/*Q1 */
                Mpy_32_16_ss(L_tmp1,tmp,&L_tmp,&lo);
                B1 = extract_h(L_shl(L_tmp,add(exp,1)));/*Q0 */
                IF(L_sub(L_tmp1,L_mult(B1,Rsum_fx)) > 0 && L_sub(L_tmp1,L_mult(add(B1,1),Rsum_fx)) >= 0)
                {
                    B1 = add(B1,1);
                }
                Mpy_32_16_ss(1975684956,shl(B,5),&L_tmp1,&lo);
                Mpy_32_16_ss(L_tmp1,shl(Rsum_sub_fx[1],7),&L_tmp1,&lo);
                Mpy_32_16_ss(L_tmp1,tmp,&L_tmp,&lo);
                B2 = extract_h(L_shl(L_tmp,sub(exp,11)));/*Q0 */
                IF(L_sub(L_tmp1,L_shl(L_mult(B2,Rsum_fx),12)) > 0 && L_sub(L_add(L_tmp1,2),L_shl(L_mult(add(B2,1),Rsum_fx),12)) >= 0)
                {
                    B2 = add(B2,1);
                }
                tmp = sub(B,B1);
                B3 = sub(tmp,B2);
            }
        }
    }

    IF ( sub(Rsum_sub_fx[2],3) < 0 )
    {
        B2 = add(B2,B3);
        B3 = 0;
        move16();
    }

    /* Bit allocation in group */
    Bit_group_fx( y, 0, SFM_G1, B1, 5, Rsubband_buf, factor_fx);
    Bit_group_fx( y, SFM_G1, SFM_G1G2, B2, 6, Rsubband_buf, factor_fx);
    Bit_group_fx( y, SFM_G1G2, BANDS, B3, 7, Rsubband_buf, factor_fx);
    FOR (i = 0; i < BANDS; i++)
    {
        Rsubband[i] = extract_l(L_shr(Rsubband_buf[i], 18));
        move16();
    }

    /* Calcuate total used bits and initialize R to be used for Noise Filling */
    L_tmp = 0;
    move32();
    FOR( i = 0; i < N; i++)
    {
        L_tmp = L_add(L_tmp,Rsubband_buf[i]);/*Q21 */
        R[i] = extract_h(L_shr(Rsubband_buf[i],5));/*Q0 */
    }
    t_fx = extract_h(L_shr(L_tmp, 5)); /*Q0 */

    return (Word16)t_fx;
}
