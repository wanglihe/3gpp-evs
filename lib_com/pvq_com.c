
/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"
#include "stl.h"

/*-------------------------------------------------------------------*
 * Local definitions
 *-------------------------------------------------------------------*/


typedef unsigned long long ui64_t;


short  shrtCDivSignedApprox(
    const short  num,
    const short  den
)
{
    Word16 pool_part;
    pool_part   =  extract_h( L_mult( negate(abs_s(num)), lim_neg_inv_tbl_fx[den] ));
    /* neg_in always, positive out always, so that positive truncation(rounding) towards zero is used */
    if ( num < 0 )
    {
        pool_part  = negate(pool_part);  /* make negative,  one op */
    }
    return pool_part;
}

int intLimCDivPos( int  NUM,
                   short DEN
                 )
{
    int result;

    result = (int )(((ui64_t) (NUM << 1)* intLimCDivInvDQ31[DEN]) >> 32) ;

    return result;
}

int   intLimCDivSigned( int   NUM,
                        short DEN
                      )
{
    int result;
    int tmp_num;
    tmp_num  = max(NUM,-NUM);
    /*    result   = (int)sign(NUM)*intLimCDivPos(tmp_num, DEN); */
    result =  (int)sign(NUM)*(int)(((ui64_t) (tmp_num << 1)* intLimCDivInvDQ31[DEN]) >> 32) ;
    return result;
}




static void  nearProjQ15(
    short  x,
    short *result
)
{
    int  NearProjC[4] = {0x00003a77, 0xffff9c52, 0x00000d57, 0x00007e5f};
    int  P12, P3;

    P12 = ((int)x * NearProjC[0] >> 16) + NearProjC[1];
    P3 = ((int)x * P12 >> 14) + NearProjC[2];
    *result =(short) (((int)x * P3 >> 15) + NearProjC[3]);

    return;
}


void obtainEnergyQuantizerDensity(
    short L_in,
    short R_in,
    short *Density
)
{
    int  Rnrg;
    short den;
    int  L, R;
    L=L_in;
    R=R_in;
    den = (2*L -1) ;
    if( den<=67)
    {
        Rnrg =  (int)(((ui64_t) (R << 1)* intLimCDivInvDQ31[den]) >> 32) ;
    }
    else
    {
        Rnrg = (R / den);
    }
    Rnrg += 28;

    if ((R - 96 < 56) && (R - 96 < Rnrg))
    {
        Rnrg = R - 96;
    }
    else if (56 < Rnrg)
    {
        Rnrg = 56;
    }

    Rnrg     = max(0, Rnrg);  /* _f table set to 1 for low Rnrg */
    *Density = (short)obtainEnergyQuantizerDensity_f[Rnrg];
    return ;
}


void   dsDirac2Dirac(
    short dsDiracIndex,
    short *diracs
)
{
    short dsHighIndex;
    dsHighIndex = dsDiracIndex - DS_INDEX_LINEAR_END - 1;
    *diracs     = dsDiracIndex;
    if (dsDiracIndex > DS_INDEX_LINEAR_END)
    {
        *diracs = dsHighDiracsTab[dsHighIndex];
    }
    return ;
}


void dsDiracPerQuanta(short  td,
                      short  t_quanta,
                      short  dsm,
                      const unsigned char* const *frQuanta,
                      short *DsIdx
                     )
{
    const unsigned char *sv;
    short nsv;
    short t_quanta_o;
    short dsIndex;
    short i;

    sv = frQuanta[td];
    nsv = sv[0];

    t_quanta_o = t_quanta - QUANTAQ3OFFSET;
    if ( t_quanta_o >= sv[nsv] )
    {
        *DsIdx = nsv;
        return;
    }
    if ( t_quanta_o <= sv[1] )
    {
        *DsIdx = 1;
        return ;
    }
    dsIndex = 1 << frQuanta[0][td];
    ;
    if (t_quanta_o > sv[nsv >> 1])
    {
        dsIndex = nsv - dsIndex;  /*single op*/
    }
    for( i = frQuanta[0][td]-1; i >= 0; i--)
    {
        /*
        tmp_diff = (t_quanta_o > sv[dsIndex] );
        dsIndex += ((tmp_diff<<1) - 1) << i;  */ /* clang warns for (neg)<<i  */
        dsIndex += shl(sub(shl(lshr(sub(sv[dsIndex], t_quanta_o), 15), 1), 1), i);
    }
    dsIndex += (t_quanta_o > sv[dsIndex]);
    dsIndex -= (dsIndex > 1);

    if (dsm == CONS)
    {
        *DsIdx = dsIndex;
        return ;
    }
    *DsIdx =  dsIndex + ((sv[dsIndex + 1] + sv[dsIndex]) < (t_quanta_o << 1)) ;
    return;
}

void   QuantaPerDsDirac(
    short td,
    short  dsDiracIndex,
    const unsigned char* const*  dimFrQuanta,
    short *Quanta
)
{
    if ( !dsDiracIndex )
    {
        *Quanta = 0;
        return;
    }
    *Quanta = ((short) dimFrQuanta[td][dsDiracIndex]) + QUANTAQ3OFFSET ;
    return ;
}



void conservativeL1Norm(
    short L,
    short Qvec,
    short Fcons,
    short Qavail,
    short Qreserv,
    short Dspec,
    short *Dvec,        /*o*/
    short *Qspare,      /*o*/
    short *Qreservplus, /*o*/
    short *Dspecplus    /*o*/
)
{
    short Minit, Mprime;
    short Qtestminus;
    const unsigned char *frQuantaL;

    frQuantaL    = hBitsN[L];

    *Qreservplus = Qreserv + Qvec - QUANTAQ3OFFSET;

    dsDiracPerQuanta(L, Qvec, Fcons, hBitsN, &Minit);
    Mprime       = Minit;
    do
    {
        Qtestminus = frQuantaL[Mprime];
        *Qspare    = Qavail - Qtestminus;
        Mprime     = Mprime - 1;
    }
    while ( Mprime >= 0 && *Qspare < QUANTAQ3OFFSET);

    if (Mprime < 0)
    {
        *Qspare = Qavail + QUANTAQ3OFFSET;
    }

    dsDirac2Dirac(Mprime + 1,Dvec);

    *Dspecplus    = Dspec + *Dvec;

    *Qreservplus -= frQuantaL[Minit];

    *Qspare      -= QUANTAQ3OFFSET;

    return;
}




void bandBitsAdjustment(
    short Brc,
    unsigned int  INTrc,
    short Bavail,
    short Nbands,
    short D,
    short L,
    short Bband,
    short Breserv,
    short *Bband_adj,  /* o */
    short *Brem,       /* o */
    short *Breservplus /* o */
)
{
    int  Btemp;
    short Bff ;
    short Dtmp, tmp_short1 ,tmp_short2;

    rangeCoderFinalizationFBits(Brc, INTrc,&Bff);

    if (D < Nbands)
    {
        Dtmp  =  min(D, 3); /* avoid branch cost */
        Btemp = intLimCDivSigned( Breserv-Bff, Dtmp);

        *Breservplus = Bband + Breserv;
    }
    else
    {
        Btemp        = 0;
        *Breservplus = Bband + Bff;
    }

    *Bband_adj = Bband;
    tmp_short1 = L * 80;
    if ( tmp_short1 < Bband )
    {
        *Bband_adj = tmp_short1;
    }

    tmp_short1 = Bavail - Bff;
    tmp_short2 = *Bband_adj + Btemp ;

    *Bband_adj = tmp_short2;
    if ( tmp_short1 < tmp_short2 )
    {
        *Bband_adj = tmp_short1 ;
    }
    if (*Bband_adj < 0)
    {
        *Bband_adj = 0;
    }
    *Brem = Bavail - Bff;
}



/*-------------------------------------------------------------------*
 * local_norm_s()
 *
 *
 *-------------------------------------------------------------------*/
static
short local_norm_s (       /* o : n shifts needed to normalize */
    short short_var      /* i : signed 16 bit variable */
)
{
    short short_res;

    if (short_var == -1 )
    {
        return (16-1);
    }
    else
    {
        if (short_var == 0)
        {
            return 0;
        }

        else
        {
            if (short_var < 0)
            {
                short_var = ~short_var;
            }

            for (short_res = 0; short_var < 0x4000; short_res++)
            {
                short_var <<= 1;
            }
        }
    }

    return (short_res);
}




static short  Ratio_base2Q11(        /* o : ,  Q11 */
    short opp,                       /* i :    Q15 */
    short near                       /* i :    Q15 */
)
{
    Word16 mc, nc, ms, ns, d, z;
    Word16 result;
    Word32 acc;

    ns = local_norm_s(opp );    /* exponent*/
    nc = local_norm_s(near );    /* exponent */

    ms  = opp << ns;   /* mantissa */
    mc  = near<< nc;   /* mantissa */

    acc = L_mac(538500224L, mc, -2776);  /* a0*mc + a1, acc(Q27), a0(Q11), a1(Q27)*/
    z = mac_r(acc, ms, -2776);           /* z in Q11, a0 in Q11 */
    d = sub(ms, mc);                     /* d in Q15 */
    z = mult_r(z, d);                    /* z in Q11 */

    result = add(z, shl(sub(nc, ns), 11));
    return result;
}

void    Ratio_rQ3(
    short opp,
    short near,
    short *result
)
{
    short  tmp = 128;
    tmp  +=  Ratio_base2Q11(opp, near);
    *result = (tmp >> 8);
    return ;
}





void densityAngle2RmsProjDec(
    short  D,         /* i: density */
    short  indexphi,  /* i: decoded  index from AR dec */

    short *oppQ15,    /* o: */
    short *nearQ15,   /* o: */
    short *oppRatioQ3 /* o: */
)
{
    short phiQ14q;
    int   L_2xPhiQ14q;
    short oppTail, nearTail;
    phiQ14q = intLimCDivPos( labs(indexphi)<<13,D>>1);
    if ( (0xFFFe&D)  == 0 )  /* even  -> positive,  odd -> 0 */
    {
        phiQ14q = 1 << 13; /* one op */
    }
#define  C_TOP             ((1 << 14) - (1 << 6))
#define  C_BOT             (1 << 6)
#define  C_SHRT_MAXABS     (1L << 15)
#define  C_SHRT_MAX_POS    ((1 << 15) - 1)
#define  C_Q14_MAX         (1<<14)

    oppTail  = phiQ14q >= C_TOP;
    nearTail = phiQ14q < C_BOT;
    if (!(oppTail || nearTail))
    {
        L_2xPhiQ14q = 2*(int)phiQ14q;
        nearProjQ15( (short)(C_SHRT_MAXABS - L_2xPhiQ14q )  , oppQ15 );
        nearProjQ15(phiQ14q << 1, nearQ15 );
        Ratio_rQ3(*oppQ15, *nearQ15,oppRatioQ3 );
    }
    else
    {
        *oppRatioQ3 = (1 - 2*nearTail)*C_Q14_MAX ;
        *oppQ15     = oppTail * C_SHRT_MAX_POS;
        *nearQ15    = nearTail * C_SHRT_MAX_POS;
    }
#undef   C_TOP
#undef   C_BOT
#undef   C_SHRT_MAXABS
#undef   C_SHRT_MAX_POS
#undef   C_Q14_MAX
    return;
}

void densityAngle2RmsProjEnc
(
    short D,
    short phiQ14uq,
    short *indexphi,
    short *oppQ15,
    short *nearQ15,
    short *oppRatioQ3
)
{

#define C_MAX13 (1L << 13)
    *indexphi = ((int)D * phiQ14uq + C_MAX13 ) >> 14;
#undef C_MAX13

    if( (0xFFFE&D) == 0 )
    {
        *indexphi = -1; /* one op */
    }
    densityAngle2RmsProjDec(D, *indexphi, oppQ15, nearQ15, oppRatioQ3);
    return;
}

void NearOppSplitAdjustment(
    short qband,
    short qzero,
    const short Qac,
    const unsigned int  INTac,
    short qglobal,
    short FlagCons,
    short Np,
    short Nhead,
    short Ntail,
    short Nnear,
    short Nopp,
    short oppRQ3,
    short *qnear,
    short *qopp,
    short *qglobalupd
)
{
    short qac, qskew,  qboth, QIb;
    int  QIa, L_QIb, qnum;
    short qmin, qavg, Midx;

    rangeCoderFinalizationFBits(Qac, INTac, &qac);
    qboth = qband - (qac - qzero);
    qskew = 0; /* skew calc code  */
    if(Nhead > 1)
    {
        qavg = intLimCDivSigned(qboth, Np);
        dsDiracPerQuanta(Ntail, qavg, FlagCons,  hBitsN, &Midx );
        QuantaPerDsDirac( Nhead, Midx, hBitsN, &qmin);

        qskew = qavg - qmin;
        qskew = max(0, qskew);
    } /* end of skew calc code*/

    QIa = intLimCDivPos(Nopp, Nnear) + 1;
    qnum = qband + qzero - qac - qskew - Nopp * (int)oppRQ3;
    L_QIb = 0;
    if( qnum > 0 )
    {
        L_QIb = intLimCDivPos(qnum, QIa);
    }
    QIb    = (short) min(32767L,L_QIb); /* saturate, L_QIb >= 0  */
    *qnear = QIb;
    if(  QIb > qboth )
    {
        *qnear = qboth; /*single op*/
    }
    *qopp       =  qboth   - *qnear;
    *qglobalupd =  qglobal - (qac - qzero);
    return;
}

/*--------------------------------------------------------------------------*
 * apply_gain()
 *
 * Apply gain
 *--------------------------------------------------------------------------*/

void apply_gain
(
    const short *ord,                       /* i  : Indices for energy order                     */
    const short *band_start,                /* i  : Sub band start indices                       */
    const short *band_end,                  /* i  : Sub band end indices                         */
    const short num_sfm,                    /* i  : Number of bands                              */
    const float *gains,                     /* i  : Band gain vector                             */
    float *xq                         /* i/o: Float synthesis / Gain adjusted synth        */
)
{
    short band,i;
    float g;

    for ( band = 0; band < num_sfm; band++)
    {
        g = gains[ord[band]];

        for( i = band_start[band]; i < band_end[band]; i++)
        {
            xq[i] *= g;
        }
    }

    return;
}


/*--------------------------------------------------------------------------*
 * fine_gain_quant()
 *
 * Fine gain quantization
 *--------------------------------------------------------------------------*/

void fine_gain_quant
(
    Encoder_State *st,
    const short *ord,                       /* i  : Indices for energy order                     */
    const short num_sfm,                    /* i  : Number of bands                              */
    const short *gain_bits,                 /* i  : Gain adjustment bits per sub band            */
    float *fg_pred,                   /* i/o: Predicted gains / Corrected gains            */
    const float *gopt                       /* i  : Optimal gains                                */
)
{
    short band;
    short gbits;
    short idx;
    float gain_db,gain_dbq;
    float err;

    for ( band = 0; band < num_sfm; band++)
    {
        gbits = gain_bits[ord[band]];
        if ( fg_pred[band] != 0 && gbits > 0 )
        {
            err = gopt[band] / fg_pred[band];
            gain_db = 20.0f*(float)log10(err);
            idx = (short) squant(gain_db, &gain_dbq, finegain[gbits-1], gain_cb_size[gbits-1]);
            push_indice( st, IND_PVQ_FINE_GAIN, idx, gbits );

            /* Update prediced gain with quantized correction */
            fg_pred[band] *= (float)pow(10, gain_dbq * 0.05f);
        }
    }

    return;
}

/*-------------------------------------------------------------------*
 * srt_vec_ind()
 *
 * sort vector and save sorting indices
 *-------------------------------------------------------------------*/

void srt_vec_ind (
    const short *linear,      /* i: linear input */
    short *srt,         /* o: sorted output*/
    short *I,           /* o: index for sorted output  */
    short length        /* i: length of vector */
)
{
    float linear_f[MAX_SRT_LEN];
    float srt_f[MAX_SRT_LEN];

    mvs2r(linear, linear_f, length);
    srt_vec_ind_f(linear_f, srt_f, I, length);
    mvr2s(srt_f, srt, length);

    return;
}

/*-------------------------------------------------------------------*
 * srt_vec_ind_f()
 *
 * sort vector and save sorting indices, using float input values
 *-------------------------------------------------------------------*/

void srt_vec_ind_f(
    const float *linear,      /* i: linear input */
    float *srt,         /* o: sorted output*/
    short *I,           /* o: index for sorted output  */
    short length        /* i: length of vector */
)
{
    short pos,npos;
    short idxMem;
    float valMem;

    /*initialize */
    for (pos = 0; pos < length; pos++)
    {
        I[pos] = pos;
    }

    mvr2r(linear, srt,length);

    /* now iterate */
    for (pos = 0; pos < (length - 1); pos++)
    {
        for (npos = (pos + 1); npos < length;  npos++)
        {
            if (srt[npos] < srt[pos])
            {
                idxMem    = I[pos];
                I[pos]    = I[npos];
                I[npos]   = idxMem;

                valMem    = srt[pos];
                srt[pos]  = srt[npos];
                srt[npos] = valMem;
            }
        }
    }

    return;
}


unsigned int UMult_32_32(unsigned int UL_var1, unsigned int UL_var2)
{
    ui64_t tmp;

    tmp = (ui64_t)UL_var1 * (ui64_t)UL_var2;
    return (unsigned int)(tmp >> 32);
}

/*-------------------------------------------------------------------*
 *  UL_div
 *
 *  Calculate UL_num/UL_den. UL_num assumed to be Q31, UL_den assumed
 *  to be Q32, then result is in Q32.
 *-------------------------------------------------------------------*/
static unsigned int UL_div(const unsigned int UL_num, const unsigned int UL_den)
{
    unsigned int UL_e, UL_Q;
    unsigned int UL_msb;
    short i;

    UL_e = 0xffffffff - UL_den;
    UL_Q = UL_num;

    for (i = 0; i < 5; i++)
    {
        UL_msb = UMult_32_32(UL_Q, UL_e);
        UL_Q = UL_Q + UL_msb;
        UL_e = UMult_32_32(UL_e, UL_e);
    }

    return UL_Q;
}
/*-------------------------------------------------------------------*
 *  UL_inverse
 *
 *  Calculate inverse of UL_val. Output in Q_exp.
 *-------------------------------------------------------------------*/
unsigned int UL_inverse(const unsigned int UL_val, short *exp)
{
    unsigned int UL_tmp;


    *exp = norm_ul(UL_val);  /* aligned to BASOP */
    UL_tmp = UL_val << (*exp);   /* Q32*/

    *exp = 32 + 31 - *exp;

    return UL_div(0x80000000, UL_tmp);
}

/*-----------------------------------------------------------------------------
 * ratio()
 *
 * Divide the numerator by the denominator.
 *----------------------------------------------------------------------------*/
Word16 ratio(const Word32 numer, const Word32 denom, Word16 *expo)
{
    Word16 expNumer, expDenom;
    Word16 manNumer, manDenom;
    Word16 quotient;

    expDenom = norm_l(denom);                     /* exponent */
    manDenom = extract_h(L_shl(denom, expDenom)); /* mantissa */
    expNumer = norm_l(numer);                     /* exponent */
    manNumer = extract_h(L_shl(numer, expNumer)); /* mantissa */
    manNumer = shr(manNumer, 1); /* Ensure the numerator < the denominator */
    quotient = div_s(manNumer, manDenom); /* in Q14 */

    *expo = sub(expNumer, expDenom);

    return quotient; /* Q14 */
}

/*-----------------------------------------------------------------------------
 * atan2_fx():
 *
 * Approximates arctan piecewise with various 4th to 5th order least square fit
 * polynomials for input in 5 segments:
 *   - 0.0 to 1.0
 *   - 1.0 to 2.0
 *   - 2.0 to 4.0
 *   - 4.0 to 8.0
 *   - 8.0 to infinity
 *---------------------------------------------------------------------------*/
Word16 atan2_fx(  /* o: Angle between 0 and EVS_PI/2 radian (Q14) */
    const Word32 y,  /* i: Argument must be positive (Q15) */
    const Word32 x   /* i: Q15 */
)
{
    Word32 acc, arg;
    Word16 man, expo, reciprocal;
    Word16 angle, w, z;

    IF (L_sub(x, 0) == 0)
    {
        return 25736; /* EVS_PI/2 in Q14 */
    }
    man = ratio(y, x, &expo); /*  man in Q14 */
    expo = sub(expo, (15 - 14)); /*  Now, man is considered in Q15 */
    arg = L_shr((Word32)man, expo);

    IF (L_shr(arg, 3+15) != 0)
    /*===============================*
     *      8.0 <= x < infinity      *
     *===============================*/
    {
        /* atan(x) = EVS_PI/2 - 1/x + 1/(3x^3) - 1/(5x^5) + ...
         *         ~ EVS_PI/2 - 1/x, for x >= 8.
         */
        expo = norm_l(arg);
        man = extract_h(L_shl(arg, expo));
        reciprocal = div_s(0x3fff, man);
        expo = sub(15 + 1, expo);
        reciprocal = shr(reciprocal, expo);   /*  Q14*/
        angle = sub(25736, reciprocal);       /* Q14   (EVS_PI/2 - 1/x) */

        /* For 8.0 <= x < 10.0, 1/(5x^5) is not completely negligible.
         * For more accurate result, add very small correction term.
         */
        IF (L_sub(L_shr(arg, 15), 10L) < 0)
        {
            angle = add(angle, 8); /* Add tiny correction term. */
        }
    }
    ELSE IF (L_shr(arg, 2+15) != 0)
    /*==========================*
     *      4.0 <= x < 8.0      *
     *==========================*/
    {
        /* interval: [3.999, 8.001]
         * atan(x) ~ (((a0*x     +   a1)*x   + a2)*x   + a3)*x   + a4
         *         = (((a0*8*y   +   a1)*8*y + a2)*8*y + a3)*8*y + a4   Substitute 8*y -> x
         *         = (((a0*8^3*y + a1*8^2)*y + a2*8)*y + a3)*8*y + a4
         *         = (((    c0*y +     c1)*y +   c2)*y + c3)*8*y + c4,
         *  where y = x/8
         *   and a0 = -1.28820869667651e-04, a1 = 3.88263533346295e-03,
         *       a2 = -4.64216306484597e-02, a3 = 2.75986060068931e-01,
         *       a4 = 7.49208077809799e-01.
         */
        w = extract_l(L_shr(arg, 3));              /* Q15  y = x/8 */
        acc = 533625337L;                          /* Q31  c1 = a1*8^2 */       move32();
        z = mac_r(acc, w, -2161);                  /* Q15  c0 = a0*8^3 */
        acc = -797517542L;                         /* Q31  c2 = a2*8 */         move32();
        z = mac_r(acc, w, z);                      /* Q15 */
        acc = 592675551L;                          /* Q31  c3 = a3 */           move32();
        z = mac_r(acc, w, z);                      /* z (in:Q15, out:Q12) */
        acc = 201114012L;                          /* Q28  c4 = a4 */           move32();
        acc = L_mac(acc, w, z);                    /* Q28 */
        angle = extract_l(L_shr(acc, (28 - 14)));  /* Q14 result of atan(x), where 4 <= x < 8 */
    }
    ELSE IF (L_shr(arg, 1+15) != 0)
    /*==========================*
     *      2.0 <= x < 4.0      *
     *==========================*/
    {
        /* interval: [1.999, 4.001]
         * atan(x) ~ (((a0*x    + a1)*x   +   a2)*x   + a3)*x   + a4
         *         = (((a0*4*y  + a1)*4*y +   a2)*4*y + a3)*4*y + a4   Substitute 4*y -> x
         *         = (((a0*16*y + a1*4)*y +   a2)*4*y + a3)*4*y + a4
         *         = (((a0*32*y + a1*8)*y + a2*2)*2*y + a3)*4*y + a4
         *         = (((   c0*y +   c1)*y +   c2)*2*y + c3)*4*y + c4,
         *  where y = x/4
         *   and a0 = -0.00262378195660943, a1 = 0.04089687039888652,
         *       a2 = -0.25631148958325911, a3 = 0.81685854627399479,
         *       a4 = 0.21358070563097167
         * */
        w = extract_l(L_shr(arg, 2));              /* Q15  y = x/4 */
        acc = 702602883L;                          /* Q31  c1 = a1*8 */         move32();
        z = mac_r(acc, w, -2751);                  /* Q15  c0 = a0*32 */
        acc = -1100849465L;                        /* Q31  c2 = a2*2 */         move32();
        z = mac_r(acc, w, z);                      /* z (in:Q15, out:Q14) */
        acc = 877095185L;                          /* Q30  c3 = a3 */           move32();
        z = mac_r(acc, w, z);                      /* z (in:Q14, out:Q12) */
        acc = 57332634L;                           /* Q28  c4 = a4 */           move32();
        acc = L_mac(acc, w, z);                    /* Q28 */
        angle = extract_l(L_shr(acc, (28 - 14)));  /* Q14  result of atan(x) where 2 <= x < 4 */
    }
    ELSE IF (L_shr(arg, 15) != 0)
    /*==========================*
     *      1.0 <= x < 2.0      *
     *==========================*/
    {
        /* interval: [0.999, 2.001]
         * atan(x) ~ (((a0*x   +  1)*x   + a2)*x   +   a3)*x   + a4
         *         = (((a0*2*y + a1)*2*y + a2)*2*y +   a3)*2*y + a4    Substitute 2*y -> x
         *         = (((a0*4*y + a1*2)*y + a2)*2*y +   a3)*2*y + a4
         *         = (((a0*4*y + a1*2)*y + a2)*y   + a3/2)*4*y + a4
         *         = (((  c0*y +   c1)*y + c2)*y   +   c3)*4*y + c4,
         *  where y = x/2
         *   and a0 = -0.0160706457245251, a1 = 0.1527106504065224,
         *       a2 = -0.6123208404800871, a3 = 1.3307896976322915,
         *       a4 = -0.0697089375247448
         */
        w = extract_l(L_shr(arg, 1));              /* Q15  y= x/2 */
        acc = 655887249L;                          /* Q31  c1 = a1*2 */         move32();
        z = mac_r(acc, w, -2106);                  /* Q15  c0 = a0*4 */
        acc = -1314948992L;                        /* Q31  c2 = a2 */           move32();
        z = mac_r(acc, w, z);
        acc = 1428924557L;                         /* Q31  c3 = a3/2 */         move32();
        z = mac_r(acc, w, z);                      /* z (in:Q15, out:Q13) */
        acc = -37424701L;                          /* Q29  c4 = a4 */           move32();
        acc = L_mac(acc, w, z);                    /* Q29 */
        angle = extract_l(L_shr(acc, (29 - 14)));  /* Q14  result of atan(x) where 1 <= x < 2 */
    }
    ELSE
    /*==========================*
     *      0.0 <= x < 1.0      *
     *==========================*/
    {
        /* interval: [-0.001, 1.001]
         * atan(x) ~ ((((a0*x   +   a1)*x   + a2)*x + a3)*x + a4)*x + a5
         *         = ((((a0*2*x + a1*2)*x/2 + a2)*x + a3)*x + a4)*x + a5
         *         = ((((  c0*x +   c1)*x/2 + c2)*x + c3)*x + c4)*x + c5
         *  where
         *    a0 = -5.41182677118661e-02, a1 = 2.76690449232515e-01,
         *    a2 = -4.63358392562492e-01, a3 = 2.87188466598566e-02,
         *    a4 =  9.97438122814383e-01, a5 = 5.36158556179092e-05.
         */
        w = extract_l(arg);                         /* Q15 */
        acc = 1188376431L;                          /* Q31  c1 = a1*2 */        move32();
        z = mac_r(acc, w, -3547);                   /* Q15  c0 = a0*2 */
        acc = -995054571L;                          /* Q31  c2 = a2 */          move32();
        z = extract_h(L_mac0(acc, w, z));           /* Q15  non-fractional mode multiply */
        acc = 61673254L;                            /* Q31  c3 = a3 */          move32();
        z = mac_r(acc, w, z);
        acc = 2141982059L;                          /* Q31  c4 = a4 */          move32();
        z = mac_r(acc, w, z);
        acc = 115139L;                              /* Q31  c5 = a5 */          move32();
        acc = L_mac(acc, w, z);                     /* Q31 */
        angle = extract_l(L_shr(acc, 31 - 14));     /* Q14  result of atan(x), where 0 <= x < 1 */
    }
    return angle;  /* Q14 between 0 and EVS_PI/2 radian. */
}
