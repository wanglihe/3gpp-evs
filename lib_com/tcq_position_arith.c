/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "options.h"
#include "prot.h"
#include "rom_com.h"
#include "basop_util.h"
#include "basop_proto_func.h"

void bitstream_save_bit(PBITSTREAM pBS, int bit);
unsigned int bitstream_load_bit(PBITSTREAM pBS);
void bitstream_rollback(PBITSTREAM pBS, int numBits);


static int ar_make_model( const short *freq, short *model, int len );
static int ar_decode( PARCODEC arInst, const short *model );
static void ar_encode( PARCODEC arInst, const short *model, int symbol );
static void ar_encode_uniform( PARCODEC arInst, unsigned int data, int bits );


/* 32x16 multiply: */
Word32 Mult_32_16(Word32 a, Word16 b)
{
    Word32 result;
    UWord16 lo;
    /* use Mpy_32_16_ss(): */
    Mpy_32_16_ss(a, b, &result, &lo);
    return result;
}

/* 32x32 multiply: */
Word32 Mult_32_32(Word32 a, Word32 b)
{
    Word32 result;
    UWord32 lo;
    /* use Mpy_32_32_ss(): */
    Mpy_32_32_ss(a, b, &result, &lo);
    return result;
}

static
void set16_fx(
    Word16 y[],  /* i/o: Vector to set                       */
    const Word16 a,    /* i  : Value to set the vector to          */
    const Word16 N     /* i  : Lenght of the vector                */
)
{
    Word16 i;

    for (i=0 ; i<N ; i++)
    {
        y[i] = a;
    }

    return;
}

static
void set32_fx(
    Word32 y[],  /* i/o: Vector to set                       */
    const Word32 a,    /* i  : Value to set the vector to          */
    const Word16 N     /* i  : Lenght of the vector                */
)
{
    Word16 i;
    for (i=0 ; i<N ; i++)
    {
        y[i] = a;
    }

    return;
}

Word32 ar_div(Word32 num, Word32 denum)
{
    Word16 exp1, exp2, exp, i;
    Word32 varout;
    Word32 sign;

    sign = L_and(L_xor(num, denum), 0x80000000);

    num   = L_abs(num);
    denum = L_abs(denum);

    if(L_sub(num, denum) < 0 || denum == 0)
    {
        return 0;
    }
    else if(L_sub(num, denum) == 0)
    {
        return 1;
    }
    else
    {
        exp1 = norm_l(num);
        exp2 = norm_l(denum);
        exp = sub(exp2, exp1);
        denum = L_shl(denum, exp);
        exp = add(exp, 1);
        varout = 0;
        for (i = 0; i < exp; i++)
        {
            num = L_sub(num, denum);
            varout = L_shl(varout, 1);
            if (num >= 0)
            {
                num = L_shl(num, 1);
                varout = L_add(varout, 1);
            }
            else
            {
                num = L_add(num, denum);
                num = L_shl(num, 1);
            }
        }
    }

    if( sign != 0)
    {
        varout = L_negate(varout);
    }

    return varout;
}

void srt_vec_ind_fx (
    const Word32 *linear,      /* linear input */
    Word32 *srt,         /* sorted output*/
    Word16 *I,           /* index for sorted output  */
    Word16 length
)
{
    Word16  pos,npos;
    Word16 idxMem;
    Word32  valMem;

    /*initilize */
    for (pos = 0; pos < length; pos++)
    {
        I[pos] = pos;
    }

    for( pos = 0; pos < length; pos++ )
    {
        srt[pos] = linear[pos];
    }

    /* now iterate */
    for (pos = 0; pos < (length - 1); pos++)
    {
        for (npos = (pos + 1); npos < length;  npos++)
        {
            if (L_sub(srt[npos], srt[pos]) < 0)
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

static Word32 GetBitsFromPulses_fx(Word16 m, Word16 n)
{
    Word16 i, integer_fx, temp_fx1, temp_fx2, exp1, exp2;
    Word32 temp32;
    Word32 frac_fx32;
    Word32 logCoeff_fx;
    Word16 exp = 0;
    Word32 mantissa_fx = 0;

    if (m == 0)
    {
        return 0;
    }

    for (i = 0; i < min(m, n); i++)
    {
        logCoeff_fx = L_add(L_shl(i + 1, 16), L_sub(table_logcum_fx[n+1], L_add(table_logcum_fx[i + 2], table_logcum_fx[n - i])));
        logCoeff_fx = L_add(logCoeff_fx, L_sub(table_logcum_fx[m], L_add(table_logcum_fx[i + 1], table_logcum_fx[m - i])));/*Q16 */

        integer_fx = extract_h(logCoeff_fx);/*Q0 */
        frac_fx32 = L_sub(logCoeff_fx, L_shl(integer_fx, 16));/*Q16 */

        /*ln2, 0.987, ln2 * ln2, 0.977 */
        /*temp1 = (int) (frac / 0.0625); */
        /*temp2 = frac - (float)temp1 * 0.0625f; */

        /* frac = pow(2.0, temp1 * 0.0625) * (1 + 0.693 * temp2 + 0.480 * temp2 * temp2 * 0.5);*/
        /*frac = pow_getbitsfrompulses[temp1] * (1 + 0.693f * temp2 + 0.480f * temp2 * temp2 * 0.5f); */

        temp_fx1 = extract_h(L_shl(frac_fx32, 4));
        temp_fx2 = extract_l(L_and(frac_fx32, 0xfff));/*Q16	 */

        frac_fx32 =L_add(L_mult(temp_fx2, 22708), Mult_32_16(L_mult0(temp_fx2, temp_fx2), 7864));/*Q32 */
        frac_fx32 = L_add(0x40000000, L_shr(frac_fx32, 2));/*30 */

        exp1 = norm_l(pow_getbitsfrompulses_fx[temp_fx1]);
        exp2 = norm_l(frac_fx32);
        frac_fx32 = Mult_32_32(L_shl(pow_getbitsfrompulses_fx[temp_fx1], exp1), L_shl(frac_fx32, exp2));/*21 + exp1 + 30 + exp2 - 31 */
        frac_fx32 = L_shr(frac_fx32, exp1 + exp2) + 1;/*20 */

        if (sub(exp, integer_fx) < 0)
        {
            mantissa_fx = L_shr(mantissa_fx, sub(integer_fx, exp));
            mantissa_fx = L_add(mantissa_fx, frac_fx32);

            exp = integer_fx;
        }
        else
        {
            mantissa_fx = L_add(mantissa_fx, L_shr(frac_fx32, sub(exp, integer_fx)));
        }
        if (L_sub(mantissa_fx, 0x200000) >= 0)
        {
            exp++;

            mantissa_fx = L_shr(mantissa_fx, 1);
        }
    }

    mantissa_fx = L_shl(mantissa_fx, 2);/*22 */
    temp_fx1 = extract_h(mantissa_fx);
    temp32 = L_shl(L_sub(mantissa_fx, L_deposit_h(temp_fx1)), 15);/*31 */
    exp1 = sub(norm_l(temp32), 1);
    temp32 = ar_div(L_shl(temp32, exp1), temp_fx1); /*31 + exp1 */
    temp32 = L_shr(temp32, exp1 + 1);/*30 */

    frac_fx32 = L_sub(0x40000000, L_shr(temp32, 1));/*30 */
    frac_fx32 = Mult_32_32(frac_fx32, temp32);/*29 */
    frac_fx32 = L_shr(frac_fx32, 13);/*16 */
    exp1 = norm_l(temp_fx1);
    temp_fx1 = Log2_norm_lc(L_shl(temp_fx1, exp1));/*15 */
    frac_fx32 = frac_fx32 + Mult_32_32(frac_fx32, 950680361); /* frac_fx32 *= 1/ln(2) */
    return L_add(L_deposit_h(exp), L_add(L_shl(temp_fx1, 1), frac_fx32));
}

void decode_position_ari_fx(PARCODEC pardec, Word16 size, Word16 npulses, Word16* nz, Word32* position)
{
    Word16 i, nzp;
    Word16 mode_num_nz[TCQ_MAX_BAND_SIZE];
    Word16 prob[TCQ_MAX_BAND_SIZE];

    Word32 btcq_fx, pnzp_fx;
    Word16 integer, frac;

    Word32 cp, scp, fxone, fxp1;
    Word16 stpos = 0, pos, ovrflag, temppos, storepos;

    fxone = 32768;
    fxp1  = 512*32768;
    temppos = 0;
    storepos = 0;
    ovrflag = 0;


    set16_fx( mode_num_nz, 0, TCQ_MAX_BAND_SIZE );
    set16_fx( prob, 0, TCQ_MAX_BAND_SIZE );

    for (i = 0; i < size; i++)
    {
        position[i] = 0;
    }

    if (L_sub(npulses, 1) > 0)
    {
        btcq_fx = GetBitsFromPulses_fx(npulses, size);
        for (i = 0; i < L_min(npulses, size); i++)
        {
            /*calculate the probability of #nz */

            pnzp_fx = L_sub(L_deposit_h(add(i, 1)), btcq_fx);
            pnzp_fx = L_add(pnzp_fx, L_add(L_sub(table_logcum_fx[size + 1], L_add(table_logcum_fx[i + 2], table_logcum_fx[size - i])),
                                           L_sub(table_logcum_fx[npulses], L_add(table_logcum_fx[i + 1], table_logcum_fx[npulses - i]))));
            pnzp_fx = L_add(pnzp_fx, 917498);/*16 */
            integer = extract_h(pnzp_fx);
            frac = extract_l(L_shr(L_sub(pnzp_fx, L_deposit_h(integer)), 1));/*15 */
            prob[i] = extract_h(L_shl(Pow2(integer, frac), 16));/*0 */
            if (prob[i] == 0)
            {
                prob[i] = 1;
            }
        }

        ar_make_model(prob, mode_num_nz, min(npulses, size));
        *nz = add(1, ar_decode(pardec, mode_num_nz));/*get #nz */
        nzp = *nz;

        if( nzp == 1 )
        {
            mode_num_nz[0] = MAX_AR_FREQ;
            for (i = 0; i < size; i ++)
            {
                mode_num_nz[i+1] = round_fx( L_shr( L_deposit_h( div_l( L_deposit_h(size - i - 1), size)), 1) );
            }

            position[ ar_decode(pardec, mode_num_nz) ] = 1;
        }
        else
        {
            mode_num_nz[0] = MAX_AR_FREQ;

            for( ; nzp > 0; nzp-- )
            {
                scp = fxp1;
                temppos  = 0;
                storepos = 0;

                for( i = stpos; i < size; i++)
                {
                    ovrflag = 0;

                    if( nzp == (size - i) )
                    {
                        cp = 0;
                    }
                    else
                    {
                        cp = L_sub( fxone, div_l( L_deposit_h(nzp), (size - i)) );
                    }
                    scp = Mult_32_16( scp, extract_l(cp) );
                    mode_num_nz[i+1-storepos-stpos] = round_fx( L_shl( scp, 6) );

                    if( (mode_num_nz[i+1-storepos-stpos] == 0 && scp > 0) || mode_num_nz[i-storepos-stpos] == mode_num_nz[i+1-storepos-stpos] )
                    {
                        ovrflag = 1;
                        temppos = ar_decode(pardec, mode_num_nz);
                        storepos += temppos;
                        scp = fxp1;

                        if( temppos == i-stpos) /* esc transmitted */
                        {
                            i--;
                        }
                        else
                        {
                            break;
                        }
                    }
                }
                if( !ovrflag )
                {
                    pos = ar_decode(pardec, mode_num_nz) + storepos;
                }
                else
                {
                    pos = storepos;
                }

                position[ stpos + pos] = 1;
                stpos += pos + 1;
            }
        }

    }
    else if (L_sub(npulses, 1) == 0)
    {
        *nz = npulses;
        nzp = *nz;
        mode_num_nz[0] = MAX_AR_FREQ;
        for (i = 0; i < size; i ++)
        {
            mode_num_nz[i+1] = round_fx( L_shr( L_deposit_h( div_l( L_deposit_h(size - i - 1), size)), 1) );
        }

        position[ ar_decode(pardec, mode_num_nz) ] = 1;
    }
    else
    {
        *nz = 0;
    }

    return;
}
void decode_magnitude_usq_fx(ARCODEC* pardec, Word16 size, Word16 npulses, Word16 nzpos, Word32* positions, Word32* out)
{
    Word16 i, magnp, magnzp;
    Word16 magns[TCQ_MAX_BAND_SIZE], magncout = 0;

    Word16 storemagn, ovrflag, pos, tempmagn = 0, mmodel[MAX_PULSES+2];
    Word32 cp, scp, fxone, fxp1;

    fxone = 32768;
    fxp1  = 512*32768;
    ovrflag = 0;


    set16_fx( magns, 1, TCQ_MAX_BAND_SIZE );
    if (sub(nzpos, npulses) == 0)
    {
        for (i = 0; i < size; i++)
        {
            out[i] = positions[i];
        }
        return;
    }
    else if (sub(nzpos, 1) == 0)
    {
        for (i = 0; i < size; i++)
        {
            if ( positions[i] != 0 )
            {
                out[i] = npulses;
                return;
            }
        }
    }

    magnzp = sub(nzpos, 1);
    magnp = sub(npulses, 1);

    magncout = 0;

    set32_fx( out, 0, size );
    set16_fx( mmodel, 0, MAX_PULSES+2 );

    mmodel[0] = MAX_AR_FREQ;
    magncout = 0;
    for( pos = 0; pos < size; pos++)
    {
        scp = fxp1;
        if( positions[pos] != 0)
        {
            storemagn = 0;

            for( i = 0; i < magnp; i++)
            {
                ovrflag = 0;

                if( magnzp == (magnp-i) )
                {
                    cp = 0;
                }
                else
                {
                    cp = L_sub( fxone, div_l( L_deposit_h(magnzp), magnp-i) );
                }

                if( cp == fxone )
                {
                    break;
                }

                scp = Mult_32_16( scp, extract_l(cp) );
                mmodel[i+1-storemagn] = round_fx( L_shl( scp, 6) );

                if( (mmodel[i+1- storemagn] == 0 && scp > 0) || mmodel[i- storemagn] == mmodel[i+1- storemagn] )
                {
                    mmodel[i+1-storemagn] = 0;
                    /* read data */
                    tempmagn = ar_decode( pardec, mmodel );
                    storemagn += tempmagn;

                    if( tempmagn < i )
                    {
                        /* just magnitude */
                        ovrflag = 1;
                        break;
                    }
                    else
                    {
                        /* esc code */
                        scp = fxp1;
                        i--;
                    }
                }
            }

            if( ovrflag )
            {
                out[magncout] = storemagn + 1;
            }
            else
            {
                out[magncout] = ar_decode( pardec, mmodel ) + storemagn + 1;
            }
            magnp -= out[magncout];
            magnzp--;
            magncout++;

            if (magnzp == 0) /* last magnitude generation */
            {
                for( pos = pos+1; pos < size; pos++)
                {
                    if( positions[pos] != 0)
                    {
                        out[magncout] = magnp + 1;
                        return;
                    }
                    else
                    {
                        out[magncout] = 0;
                        magncout++;
                    }
                }
            }
            else if(magnzp == magnp) /* rest magnitudes generation */
            {
                for( pos = pos+1; pos < size; pos++)
                {
                    out[magncout] = positions[pos];
                    magncout++;
                }
                return;
            }
        }
        else
        {
            out[magncout] = 0;
            magncout++;
        }
    }

    return;
}

static Word16 quantize_fx( Word16 val, Word16 D)
{
    Word16 qval4_fx;
    Word16 retval_fx;

    qval4_fx = shr(abs_s(add(val, 512)), 12);
    retval_fx = add(shl(qval4_fx, 2), DDP_fx[D]);
    /* 2nd zero check */
    if (D == 0)
    {
        if (sub(abs_s(sub(shl(abs_s(retval_fx), 10), abs_s(val))), abs_s(val)) > 0)
        {
            retval_fx = 0;
        }
    }
    return retval_fx;
}

void decode_mangitude_tcq_fx(ARCODEC* pardec, Word16 size, Word16 npulses, Word16 nzpos, Word32* positions, Word32* out, Word32* surplus_fx)
{
    Word32 tcq_bits_fx, bits_fx/*, surplus_fx*/;
    Word16 prob0_fx, prob1_fx, num, denum, quantum1_fx, quantum2_fx;
    Word16 exp, exp1, exp2, tmp16;
    Word32 tmp32;

    Word16 i, j, symbol, st;
    Word16 leftp  = npulses;/*pulsesnum; */
    Word16 leftnz = nzpos; /*nzpos; */
    Word16 magn_mode[3] = {MAX_AR_FREQ, 0, 0};

    bits_fx = 0;
    tcq_bits_fx = L_sub(table_logcum_fx[npulses], L_add(table_logcum_fx[nzpos], table_logcum_fx[npulses - (nzpos - 1)]));

    if (sub(nzpos, npulses) == 0)
    {
        for (i = 0; i < size; i++)
        {
            out[i] = positions[i];
        }

        return;
    }
    else if (sub(nzpos, 1) == 0)
    {
        for (i = 0; i < size; i++)
        {
            if ( positions[i] != 0 )
            {
                out[i] = npulses;
                return;
            }
        }
    }
    st = 0;
    for (i = 0; i < size && leftnz > 1; i++)
    {
        out[i] = positions[i];
        if (positions[i] != 0)
        {
            /*generate the trellis path */
            symbol = 0;
            for (j = 0; j < leftp; j++)
            {
                num = sub(leftnz, 1);
                denum = sub(leftp, add(j, 1));
                if (sub(num, denum) >= 0)
                {
                    prob1_fx = MAX_16;
                    prob0_fx = 0;
                }
                else
                {
                    exp1 = sub(norm_s(num), 1);
                    exp2 = norm_s(denum);
                    prob1_fx = div_s(shl(num, exp1), shl(denum, exp2));/*15 + exp1 - exp2 */
                    exp = 15 + exp1 - exp2;
                    prob1_fx = shl(prob1_fx, sub(15, exp));
                    prob0_fx = sub(MAX_16, prob1_fx);
                }
                if (L_sub(sub(leftp, j), leftnz) == 0)
                {
                    symbol = add(j, 1);
                    break;
                }

                quantum1_fx = quantize_fx(shl(add(j, 1), 10), ddec[st][0]);
                quantum2_fx = quantize_fx(shl(add(j, 1), 10), ddec[st][1]);

                if (sub(quantum1_fx, add(j, 1)) != 0 && sub(quantum2_fx, add(j, 1)) != 0)
                {
                    prob0_fx = MAX_16;
                    move16();
                    prob1_fx = 0;
                    move16();
                }
                if (sub(prob0_fx, MAX_16) == 0 || sub(prob1_fx, MAX_16) == 0)
                {
                    symbol = add(j, 1);
                    continue;
                }

                /*magn_mode[1] = (short)(prob1 * MAX_AR_FREQ); */
                magn_mode[1] = mult(prob1_fx, MAX_AR_FREQ);

                if (ar_decode(pardec, magn_mode))
                {
                    exp1 = norm_s(prob1_fx);
                    tmp32 = L_deposit_h(shl(prob1_fx, exp1));/*exp1 + 15 + 16 */
                    tmp16 = Log2_norm_lc(tmp32);/*15 */
                    bits_fx = L_sub(bits_fx, L_sub(tmp16, L_shl(add(exp1, 1), 15)));/*15 */

                    symbol = add(j, 1);
                    break;
                }
                else
                {
                    exp1 = norm_s(prob0_fx);
                    tmp32 = L_deposit_h(shl(prob0_fx, exp1));/*exp1 + 15 + 16 */
                    tmp16 = Log2_norm_lc(tmp32);/*15 */
                    bits_fx = L_sub(bits_fx, L_sub(tmp16, L_shl(add(exp1, 1), 15)));/*15 */
                }
            }
            out[i] = symbol;
            /*leftp -= symbol; */
            leftp = sub(leftp, symbol);
            leftnz--;
        }

        quantum1_fx = quantize_fx(out[i], ddec[st][0]);
        quantum2_fx = quantize_fx(out[i], ddec[st][1]);

        /*generate the next state */
        if (sub(quantum1_fx, out[i]) == 0)
        {
            st = nextstate[st][0];
        }
        else
        {
            st = nextstate[st][1];
        }
    }

    /*generate the magnitudes */
    for (; i < size; i++)
    {
        out[i] = 0;
        if (positions[i] != 0)
        {
            out[i] = add(sub(leftp, leftnz), 1);
        }

    }

    if (sub(nzpos, npulses) != 0 && sub(nzpos, 1) > 0)
    {
        /*update the surplus */
        *surplus_fx = L_add(*surplus_fx, L_sub(tcq_bits_fx, L_shl(bits_fx, 1)));
    }

    return;
}

Word16 GetScale_fx( Word16 blen, Word32 bits_fx/*Q16*/, Word32 *surplus_fx/*Q16*/)
{
    Word16 pulses = MAX_PULSES, p_est, exp, exp1, exp2, magicnum;
    Word32 t, a, b, ab, estbits_fx = 0;

    magicnum = 24773; /*Q17: 0.188992013101951f; */

    t = L_shr( L_mult( magicnum, blen), 2);
    exp = norm_l(t);
    a = L_shl( 14 - exp, 15) + Log2_norm_lc( L_shl( t, exp ) );

    exp1 = sub( norm_l(bits_fx), 1);
    exp2 = norm_s( blen - 1 );
    b = L_shr( L_deposit_l( div_l( L_shl( bits_fx, exp1), shl( blen - 1, exp2) ) ), exp1-exp2 );

    ab = L_add( a, b);

    p_est = extract_h(L_shl(Pow2( extract_l( L_shr(ab,15) ), ab&0x7FFF ),16));  /* enable saturation of pow2 result  */

    pulses = min( p_est, MAX_PULSES );

    for( ; pulses >= 0; pulses--)
    {
        estbits_fx = GetBitsFromPulses_fx( pulses, blen);
        if( L_sub( bits_fx, estbits_fx) >= 0)
        {
            break;
        }
    }

    if ( surplus_fx != 0 )
    {
        *surplus_fx = L_add(*surplus_fx, L_sub(bits_fx, estbits_fx));
    }

    return pulses;
}

void decode_signs_fx(ARCODEC* pardec, Word16 size, Word32* out)
{
    Word16 i;

    for ( i = 0; i < size; i++)
    {
        if ( out[i] != 0  )
        {
            out[i] = ( ar_decode( pardec, uniform_model ) > 0) ? out[i] : -out[i];
        }
    }

    return;
}

Word32 encode_position_ari_fx(PARCODEC parenc, float* quants, Word16 size, Word32* est_bits_frame_fx)
{
    Word16 i;
    Word16 nz = 0, pulses = 0;
    Word16 prob[TCQ_MAX_BAND_SIZE];
    Word16 model_num_nz[TCQ_MAX_BAND_SIZE];
    float *cur_quants = quants;
    Word16 integer, frac;
    Word32 /*est_bits_frame_fx, */btcq_fx = 0, bits_fx = 0, pnzp_fx;
    Word32 cp, scp, fxone, fxp1;
    Word16 pos;

    fxone = 32768;
    fxp1  = 512*32768;

    set16_fx( prob, 0, TCQ_MAX_BAND_SIZE );
    set16_fx( model_num_nz, 0, TCQ_MAX_BAND_SIZE );

    for (i = 0; i < size; i ++)
    {
        pulses = add(pulses, abs_s((short)cur_quants[i]));
        if (cur_quants[i] != 0)
        {
            nz++;
        }
    }

    btcq_fx = GetBitsFromPulses_fx(pulses, size);
    /* Estimate TCQ bits */
    bits_fx = L_sub(table_logcum_fx[size + 1], L_add(table_logcum_fx[nz + 1], table_logcum_fx[size - nz + 1]));
    bits_fx = L_add(bits_fx, L_sub(btcq_fx, L_sub(table_logcum_fx[size + 1], L_add(table_logcum_fx[nz + 1], table_logcum_fx[size - nz + 1]))));
    bits_fx = L_sub(bits_fx, L_sub(table_logcum_fx[pulses], L_add(table_logcum_fx[nz], table_logcum_fx[pulses - (nz - 1)])));

    bits_fx = L_sub(bits_fx, nz);
    *est_bits_frame_fx = L_add(*est_bits_frame_fx, bits_fx);

    /*caculate the #nz probability */
    for (i = 0; i < min(pulses, size); i++)
    {
        pnzp_fx = L_sub(L_deposit_h(add(i, 1)), btcq_fx);

        pnzp_fx = L_add(pnzp_fx, L_add(L_sub(table_logcum_fx[size + 1], L_add(table_logcum_fx[i + 2], table_logcum_fx[size - i])),
                                       L_sub(table_logcum_fx[pulses], L_add(table_logcum_fx[i + 1], table_logcum_fx[pulses - i]))));

        pnzp_fx = L_add(pnzp_fx, 917498);/*16 */
        integer = extract_h(pnzp_fx);
        frac = extract_l(L_shr(L_sub(pnzp_fx, L_deposit_h(integer)), 1));/*15 */
        prob[i] = extract_h(L_shl(Pow2(integer, frac), 16));/*0 */

        /*zero probability will incur problems in ar_make_model() */
        if (prob[i] == 0)
        {
            prob[i] = 1;
        }
    }

    ar_make_model(prob, model_num_nz, min(pulses, size));

    if (sub(nz, 1) > 0)
    {
        ar_encode( parenc, model_num_nz, nz - 1);/*encode #nz */

        scp = fxp1;
        pos = 0;
        for( i = 0; i < size && nz > 0; i++)
        {
            if( nz == (size - i) )
            {
                cp = 0;
            }
            else
            {
                cp = L_sub( fxone, div_l( L_deposit_h(nz), (size - i)) );
            }
            scp = Mult_32_16( scp, extract_l(cp) );
            model_num_nz[pos+1] = round_fx( L_shl( scp, 6) );

            if( (model_num_nz[pos+1] == 0 && scp > 0) || model_num_nz[pos] == model_num_nz[pos+1] )
            {
                model_num_nz[pos+1] = 0;
                ar_encode( parenc, model_num_nz, pos );
                i--;
                scp = fxp1;
                pos = 0;
                continue;
            }

            if( cur_quants[i] != 0 )
            {
                ar_encode( parenc, model_num_nz, pos );
                pos = 0;
                scp = fxp1;
                nz--;
            }
            else
            {
                pos++;
            }
        }
    }
    else if (sub(nz, 1) == 0)
    {
        if (sub(pulses, 1) > 0)
        {
            /*temp -= log2_f((float)(model_num_nz[nz-1] - model_num_nz[nz]) / MAX_AR_FREQ); */
            ar_encode(parenc, model_num_nz, 0);/*encode #nz */
        }

        pos = 0;
        for( i = 0; i < size; i++)
        {
            model_num_nz[i+1] = round_fx( L_shr( L_deposit_h( div_l( L_deposit_h(size - i - 1), size)), 1) );

            if( cur_quants[i] != 0 )
            {
                pos = i;
            }
        }
        ar_encode( parenc, model_num_nz, pos ); /* encode pos */
    }
    return bits_fx;
}

Word32 encode_magnitude_tcq_fx(ARCODEC* parenc, float* magn_fx, Word16 size, Word16 npulses, Word16 nzpos, Word32* savedstates, Word32* est_frame_bits_fx)
{
    Word32 tcq_bits_fx, bits_fx/*, est_frame_bits_fx*/;
    Word16 prob0_fx, prob1_fx, num, denum, quantum1_fx, quantum2_fx;
    Word16 exp, exp1, exp2;

    Word16 i, j;
    Word32 st;
    Word16 magn_mode[3] = {MAX_AR_FREQ, 0, 0};

    Word16 leftp  = npulses;/*pulsesnum; */
    Word16 leftnz = nzpos;/*nzpos; */

    bits_fx = 0;

    tcq_bits_fx = L_sub(table_logcum_fx[npulses], L_add(table_logcum_fx[nzpos], table_logcum_fx[npulses - (nzpos - 1)]));

    *est_frame_bits_fx = L_add(*est_frame_bits_fx, tcq_bits_fx);

    if (sub(nzpos, npulses) == 0 || sub(nzpos, 1) == 0)
    {
        return bits_fx;
    }

    st = 0;
    for ( i = 0; i < size && leftnz > 1; i++)
    {
        st = savedstates[i];
        if (magn_fx[i] != 0)
        {
            for ( j = 0; j < leftp; j++)
            {
                /*calculate the two path probs point to next two states */
                num = sub(leftnz, 1);
                denum = sub(leftp, add(j, 0x1));
                if (sub(num, denum) >= 0)
                {
                    prob1_fx = MAX_16;
                    prob0_fx = 0;
                }
                else
                {
                    exp1 = sub(norm_s(num), 1);
                    exp2 = norm_s(denum);
                    prob1_fx = div_s(shl(num, exp1), shl(denum, exp2));/*15 + exp1 - exp2 */
                    exp = 15 + exp1 - exp2;
                    prob1_fx = shl(prob1_fx, sub(15, exp));
                    prob0_fx = sub(MAX_16, prob1_fx);
                }

                quantum1_fx = quantize_fx(shl(add(j, 1), 10), ddec[st][0]);
                quantum2_fx = quantize_fx(shl(add(j, 1), 10), ddec[st][1]);

                if (sub(quantum1_fx, add(j, 1)) != 0 && sub(quantum2_fx, add(j, 1)) != 0)
                {
                    prob0_fx = MAX_16;
                    prob1_fx = 0;
                }
                if (sub(prob0_fx, MAX_16) == 0 || sub(prob1_fx, MAX_16) == 0)
                {
                    continue;
                }

                magn_mode[1] = mult(prob1_fx, MAX_AR_FREQ);
                if (sub(j, sub(abs_s( (short)magn_fx[i] ), 1)) < 0)
                {
                    ar_encode(parenc, magn_mode, 0);
                }
                else
                {
                    if (sub(leftp, j) > leftnz)
                    {
                        ar_encode(parenc, magn_mode, 1);
                    }
                    break;
                }
            }

            leftnz--;
            leftp = sub(leftp, abs_s( (short)magn_fx[i] ));
        }
    }

    return bits_fx;
}

Word32 encode_signs_fx(ARCODEC* parenc, float* magn, Word16 size, Word16 npos, Word32* est_frame_bits_fx)
{
    Word32 i, sign;

    *est_frame_bits_fx = L_add(*est_frame_bits_fx, L_deposit_h(npos));
    for (i = 0; i < size; i++)
    {
        if (magn[i] != 0)
        {
            sign = (magn[i] > 0) ? 1 : 0;
            ar_encode_uniform(parenc, sign, 1);
        }
    }

    return L_deposit_h(npos);
}

Word32 encode_magnitude_usq_fx(ARCODEC* parenc, float* magn_fx, Word16 size, Word16 npulses, Word16 nzpos, Word32* est_frame_bits_fx)
{
    Word16 i, j, k, magnp, magnzp;
    Word16 magn_position[MAX_PULSES];
    Word32 /*est_frame_bits_fx, */bits_fx;

    Word16 pos, model_m[MAX_PULSES + 2];
    Word32 fxone, fxp1, cp, scp;

    fxone = 32768;
    fxp1  = 512*32768;

    /*estimate fac bits */
    bits_fx = L_sub(table_logcum_fx[npulses], L_add(table_logcum_fx[nzpos], table_logcum_fx[npulses - nzpos + 1]));

    *est_frame_bits_fx = L_add(*est_frame_bits_fx, bits_fx);

    if (sub(npulses, nzpos) == 0 || sub(nzpos, 1) == 0)
    {
        return bits_fx;
    }
    magnp = sub(npulses, 1);
    magnzp = sub(nzpos, 1);

    /*generate the binary sequences of magnitudes */
    k = 0;
    for (i = 0; i < size; i++)
    {
        if (magn_fx[i] != 0)
        {
            for (j = 0; j < abs_s((short)magn_fx[i]) - 1; j++)
            {
                magn_position[k++] = 0;
            }
            magn_position[k++] = 1;
        }
    }

    set16_fx( model_m, 0, MAX_PULSES + 2);
    scp = fxp1;
    model_m[0] = MAX_AR_FREQ;
    pos = 0;
    for( i = 0; i < npulses-1 && magnzp > 0; i++ )
    {
        if( magnzp == magnp )
        {
            cp = 0;
        }
        else
        {
            cp = L_sub( fxone, div_l( L_deposit_h(magnzp), magnp) );
        }
        scp = Mult_32_16( scp, extract_l(cp) );
        model_m[pos+1] = round_fx( L_shl( scp, 6) );

        if( (model_m[pos+1] == 0 && scp > 0) || model_m[pos] == model_m[pos+1] )
        {
            model_m[pos+1] = 0;

            ar_encode( parenc, model_m, pos );
            pos = 0;
            i--;
            scp = fxp1;
            continue;
        }

        if( magn_position[i] != 0 )
        {
            ar_encode( parenc, model_m, pos );
            pos = 0;
            magnzp--;
            scp = fxp1;
        }
        else
        {
            pos++;
        }

        magnp--;
    }
    return bits_fx;
}


static void transmission_bits( PARCODEC arInst, int bit )
{
    bitstream_save_bit( arInst->bsInst, bit );
    arInst->num_bits++;
    bit = !bit;

    for( ; arInst->bits_to_follow > 0 && arInst->num_bits < arInst->max_bits; arInst->bits_to_follow --)
    {
        bitstream_save_bit( arInst->bsInst, bit );
        arInst->num_bits++;
    }

    return;
}

void ar_encoder_start( PARCODEC arInst, PBITSTREAM bsInst, int max_bits
                     )
{
    arInst->bsInst = bsInst;

    arInst->low	= 0;
    arInst->high = AR_TOP;
    arInst->bits_to_follow	= 0;

    arInst->num_bits = 0;
    arInst->max_bits = max_bits;
}


static void ar_encode( PARCODEC arInst, const short *model, int symbol )
{
    unsigned int range, high, low;

    high = arInst->high;
    low = arInst->low;

    symbol ++;
    range  = high - low + 1;

    high  = low + ( range * model[symbol - 1]) / model[0] - 1;
    low   = low + ( range * model[symbol] ) / model[0];

    for( ; ; )
    {
        if( high < AR_HALF )
        {
            transmission_bits( arInst, 0 );
        }
        else
        {
            if( low >= AR_HALF )
            {
                transmission_bits( arInst, 1 );

                low  -= AR_HALF;
                high -= AR_HALF;
            }
            else
            {
                if( low >= AR_FIRST && high < AR_THIRD )
                {
                    arInst->bits_to_follow ++;

                    low  -= AR_FIRST;
                    high -= AR_FIRST;
                }
                else
                {
                    break;
                }
            }
        }

        low  = low << 1;
        high = ( high << 1 ) + 1;
    }

    arInst->high = high;
    arInst->low = low;

    return;
}


static void ar_encode_uniform( PARCODEC arInst, unsigned int data, int bits )
{
    int i;

    for( i = 0 ; i < bits ; i ++ )
    {
        ar_encode( arInst, uniform_model, data & 0x1 );
        data >>= 1;
    }

    return;
}


void ar_encoder_done( PARCODEC arInst )
{
    arInst->bits_to_follow ++;
    transmission_bits( arInst, arInst->low >= AR_FIRST );

    return;
}


void ar_decoder_start( PARCODEC arInst, PBITSTREAM bsInst )
{
    int i;

    arInst->bsInst	= bsInst;

    arInst->low		= 0;
    arInst->high	= AR_TOP;
    arInst->value	= 0;

    for( i = 0; i < AR_BITS ; i ++ )
    {
        arInst->value = ( arInst->value << 1 ) + bitstream_load_bit( arInst->bsInst );
    }

    return;
}


static int ar_decode( PARCODEC arInst, const short *model )
{
    unsigned int range;
    short cum;
    int   symbol;

    range = (unsigned int)( arInst->high - arInst->low ) + 1;
    cum   = (short)( ( ( (unsigned int)( arInst->value - arInst->low ) + 1 ) * model[0] - 1 ) / range );

    for( symbol = 1 ; model[symbol] > cum ; symbol ++ );

    arInst->high	= arInst->low + ( range * model[symbol - 1] ) / model[0] - 1;
    arInst->low		= arInst->low + ( range * model[symbol] ) / model[0];

    for( ; ; )
    {
        if( arInst->high >= AR_HALF )
        {
            if( arInst->low >= AR_HALF)
            {
                arInst->value -= AR_HALF;
                arInst->low   -= AR_HALF;
                arInst->high  -= AR_HALF;
            }
            else
            {
                if( arInst->low >= AR_FIRST && arInst->high < AR_THIRD )
                {
                    arInst->value -= AR_FIRST;
                    arInst->low   -= AR_FIRST;
                    arInst->high  -= AR_FIRST;
                }
                else
                {
                    break;
                }
            }
        }
        arInst->low  <<= 1;
        arInst->high   = ( arInst->high  << 1 ) + 1;
        arInst->value  = ( arInst->value << 1 ) + bitstream_load_bit( arInst->bsInst );
    }

    return (symbol - 1);
}

void ar_decoder_done( PARCODEC arInst )
{
    bitstream_rollback( arInst->bsInst, AR_BITS - 2 );

    return;
}


static int ar_make_model( const short *freq, short *model, int len )
{
    short dist;
    unsigned int sum = 0;
    unsigned int cum = 0;
    int i;

    for( i = 0 ; i < len ; i ++ )
    {
        sum += freq[i];
    }

    if( sum == 0 )
    {
        return 0;
    }

    for( i = len ; ; i -- )
    {
        model[i] = (short)( ( cum * MAX_AR_FREQ ) / sum );

        if( !i )
        {
            break;
        }

        cum += freq[i - 1];
    }


    for( i = 0 ; i < len - 1 ; i ++ )
    {
        dist = model[i] - model[i + 1];

        if( dist <= 0  )
        {
            model[i + 1] += dist - 1;
        }
    }

    for( i = len ; i ; i -- )
    {
        dist = model[i - 1] - model[i];

        if( dist <= 0  )
        {
            model[i - 1] -= dist - 1;
        }
    }

    return (model[0] > model[1]);
}

void bitstream_save_bit(PBITSTREAM pBS, int bit)
{
    unsigned char cur;

    cur = pBS->buf[pBS->numByte];


    cur = (unsigned char)(cur | (bit << pBS->curPos--));
    pBS->buf[pBS->numByte] = cur;
    pBS->numbits++;

    if (pBS->curPos < 0)
    {
        pBS->curPos = 7;
        pBS->numByte++;
    }

    return;
}


unsigned int bitstream_load_bit(PBITSTREAM pBS)
{
    unsigned int bit;
    signed char * curPos;

    /* safety check in case of bit errors */
    if( pBS->numByte >= pBS->maxBytes )
    {
        return 0;
    }

    curPos = &pBS->curPos;
    bit = (( pBS->buf[pBS->numByte] >> (*curPos)--) & 0x00000001);

    if (*curPos < 0)
    {
        pBS->numByte++;
        *curPos = 7;
    }

    return bit;
}


void bitstream_rollback(PBITSTREAM pBS, int numBits)
{

    while (numBits > 0)
    {
        pBS->curPos++;
        pBS->numbits--;
        if (pBS->curPos == 8)
        {
            pBS->curPos = 0;
            pBS->numByte--;
        }
        numBits--;
    }

    return;
}


static float quantize( float val, int D)
{
    int qval4;
    float retval;/* = qval4*4.0f + DD[D]; */

    qval4 = (int)fabs( (val + 0.5)/4.0 );
    retval = qval4*4.0f + DDP[D];

    /* 2nd zero check */
    if( D == 0 )
    {
        if( fabs( fabs(retval) - fabs( val ) ) > fabs( val ) )
        {
            retval = 0;
        }
    }

    return retval;
}


static void TCQnew( float *v, float scale, int length, float *vout, int pulses, int *pulsesout, int* nzposout, int *savedstates, int * lasttrellislevel, int terminate)
{
    short i, st, dminpos, position, pulsesnum, nzpos = 0;
    float dmin, quantum1, quantum2, curdist1, curdist2, newdist1, newdist2, signq;

    float metric[STATES][TCQ_MAX_BAND_SIZE];
    short path[STATES][TCQ_MAX_BAND_SIZE];
    short quant[STATES][TCQ_MAX_BAND_SIZE];
    short pused[STATES][TCQ_MAX_BAND_SIZE];

    set_f( *metric, 0.0f, STATES*TCQ_MAX_BAND_SIZE );
    set_s( *path, 0, STATES*TCQ_MAX_BAND_SIZE );
    set_s( *quant, 0, STATES*TCQ_MAX_BAND_SIZE );
    set_s( *pused, 0, STATES*TCQ_MAX_BAND_SIZE );

    /* Initialize metric */
    for( st = 1; st < STATES; st++)
    {
        metric[st][0] = 1000;
    }
    for( st = 2; st < STATES; st++)
    {
        metric[st][1] = 1000;
    }
    for( st = 4; st < STATES; st++)
    {
        metric[st][2] = 1000;
    }

    /* Viterbi for input sequence */
    for( i = 0; i < length; i++) /* cycle over symbols */
    {
        for( st = 0; st < STATES; st++) /* cycle over conditions */
        {
            curdist1 = metric[ step_tcq[st][0]][i];
            curdist2 = metric[ step_tcq[st][1]][i];

            /* step 1 */
            quantum1 = quantize( v[i]*scale, denc[st][0]);
            newdist1 = (float)( quantum1 - fabs(v[i])*scale );
            newdist1 *= newdist1;

            if( quantum1 + pused[ step_tcq[st][0]][i] > pulses && terminate)
            {
                newdist1 = 10000.0f; /* pulses check */
            }

            /* step 2 */
            quantum2 = quantize( v[i]*scale, denc[st][1]);
            newdist2 = (float)( quantum2 - fabs(v[i])*scale );
            newdist2 *= newdist2;

            if( quantum2 + pused[ step_tcq[st][1]][i] > pulses && terminate)
            {
                newdist2 = 10000.0f; /* pulses check */
            }

            /* decision */
            if( curdist1 + newdist1 < curdist2 + newdist2 )
            {
                path[st][i+1] = step_tcq[st][0];
                metric[st][i+1] = curdist1 + newdist1;
                quant[st][i+1] = (int)quantize( v[i]*scale, denc[st][0]);
                pused[st][i+1] = (int)(pused[ step_tcq[st][0]][i] + abs( quant[st][i+1] ));
            }
            else
            {
                path[st][i+1] = step_tcq[st][1];
                metric[st][i+1] = curdist2 + newdist2;
                quant[st][i+1] = (int)quantize( v[i]*scale, denc[st][1]);
                pused[st][i+1] = (int)(pused[ step_tcq[st][1]][i] + abs( quant[st][i+1] ));
            }
        }
    }

    /* Find path with minimal metric */
    dminpos = 0;
    dmin    = metric[ dminpos][ length];
    for( i = 1; i < STATES; i++)
    {
        if( (dmin > metric[ i][ length] && pused[i][ length] == pulses) ||
                ( pused[dminpos][ length] != pulses && pused[i][ length] == pulses) )
        {
            dmin = metric[ i][ length];
            dminpos = i;
        }
    }

    /* Trace back to get output */
    pulsesnum = 0;
    position = dminpos;

    for( i = length; i > 0; i--)
    {
        signq = ( v[i-1] > 0.0f )?(1.0f):(-1.0f);
        vout[i-1] = signq*quant[position][i];

        position  = path[position][i];
        savedstates[i-1] = position;

        /* calculate output pulses number & nz */
        pulsesnum += (int)fabs(vout[i-1]);/*quant[position][i]; */
        if( fabs(vout[i-1]) > 0.0f )
        {
            if( nzpos == 0 ) *lasttrellislevel = i;

            nzpos++;
        }
    }

    if( pulsesout != 0 )
    {
        *pulsesout = pulsesnum;
    }
    if( nzposout != 0 )
    {
        *nzposout  = nzpos;
    }

    return;
}


float GetISCScale( float *quants, int size, Word32 bits_fx, float *magn, float *qscale, Word32 *surplus_fx, float *pulses, int* savedstates,
                   int noTCQ, int *nzpout, short *bcount, float *abuffer, float *mbuffer, float *sbuffer)
{
    float scale, m, t, actualt, magnbits = 0.0f;
    int pulsesnum, pos, terminate, leftp, leftnz, trellislevel, st;
    int i, j, nzpos, direction, pulsescurr, nzposcurr, lasttrellislevel;
    float dist[TCQ_MAX_BAND_SIZE];
    float aquants[TCQ_MAX_BAND_SIZE];
    float dmin, prob0, prob1, quantum1, quantum2;
    float sx2 = 0, sy2 = 0, sxy = 0, g;
    int   pn = 0;
    float pt = 0.f;
    int diff;
    int sign, m_int;
    int flag_g1;

    set_f( dist, 0.0f, TCQ_MAX_BAND_SIZE );
    set_f( aquants, 0.0f, TCQ_MAX_BAND_SIZE );

    if( bits_fx < 0 )
    {
        pulsesnum =  0;

        if( surplus_fx != 0 )
        {
            *surplus_fx = L_add( *surplus_fx, bits_fx);
        }
    }
    else
    {
        pulsesnum =  GetScale_fx( size, bits_fx, surplus_fx);
    }
    *nzpout = 0;

    if( pulses != 0 )
    {
        *pulses = (float)pulsesnum;
    }

    if( pulsesnum > 0 )
    {
        /* Initial quantization */
        for( i = 0, m = 0; i < size; i++)
        {
            aquants[ i] = (float)fabs( quants[i] );
            m += aquants[i];
        }

        scale = (pulsesnum + EPSILON) / (m + EPSILON);

        for( i = 0, t = 0.0f; i < size; i++)
        {
            magn[i] = (float)((int)( 0.5f + aquants[ i] * scale));

            t += magn[i];
        }

        /* Pulses redistribution */
        while( t != pulsesnum )
        {
            pn = 0;
            pt = 0.f;


            for( i = 0, nzpos = 0; i < size; i++)
            {
                if( magn[i] > 0.0f )
                {
                    pn += (int)magn[i];
                    pt += aquants[i];
                }
            }

            direction = ( pulsesnum - t > 0 )?(1):(-1);

            /* new alg */
            {
                for( i = 0; i < size; i++)
                {
                    sxy += aquants[ i]*magn[ i];
                    sx2 += aquants[ i]*aquants[ i];
                    sy2 += magn[ i]*magn[ i];
                }
                for( i = 0; i < size; i++)
                {
                    if( magn[i] > 0.0f )
                    {
                        g = (pt)/(pn + direction + EPSILON);
                    }
                    else
                    {
                        g = (pt + aquants[ i])/(pn + direction + EPSILON);
                    }

                    dist[i] = sx2 - 2.0f*(sxy + direction*aquants[i])*g + g*g*(sy2 + 2.0f*magn[i]*direction + 1.0f);
                }
            }

            {
                pos = 0;
                dmin = dist[0];

                /* find min */
                for( i = 1; i < size; i++)
                {
                    if( dmin > dist[i] )
                    {
                        pos = i;
                        dmin = dist[i];
                    }
                }

                magn[pos] += direction;
                t += direction;
            }
        }

        /* calculate actual nz positions */
        actualt = 0.0f;
        for( i = 0, nzpos = 0; i < size; i++)
        {
            if( magn[i] > 0.0f )
            {
                if (quants[i] < 0)
                {
                    magn[i] *= -1;
                }

                actualt += aquants[i];
                nzpos++;
            }
        }

        /* calculate scale */
        if( actualt > 0)
        {
            scale = pulsesnum / actualt;
        }
        else
        {
            scale = FLT_MAX;
        }
        *qscale = scale;
        *nzpout = nzpos;

        if( (nzpos != pulsesnum && nzpos > 1) && noTCQ == 0 )
        {
            terminate = 1;
            TCQnew( quants, scale, size, magn, pulsesnum, &pulsescurr, &nzposcurr, savedstates, &lasttrellislevel, terminate);

            if( pulsesnum > pulsescurr )
            {
                scale *= 1.1f;
            }

            if( pulsesnum < pulsescurr )
            {
                scale *= 0.9f;
            }
            if( pulsesnum > pulsescurr )
            {
                diff = pulsesnum - pulsescurr;

                for( i = size-1; i >=0; i--)
                {
                    if( fabs(magn[i]) > 0 )
                    {
                        sign = (magn[i]>0)?(1) : (-1);
                        magn[i] = (float)(sign*( fabs(magn[i]) + diff));

                        break;
                    }
                }
            }
            else if( pulsesnum < pulsescurr )
            {
                diff = pulsescurr - pulsesnum;

                for( i = size-1; i >=0 && diff > 0; i--)
                {
                    if( fabs(magn[i]) > 0 )
                    {
                        sign = (magn[i]>0)?(1) : (-1);
                        m_int = (int)fabs(magn[i]);

                        if( diff < m_int )
                        {
                            magn[i] = (float)(sign*( fabs(magn[i]) - diff));
                            break;
                        }
                        else
                        {
                            diff = diff - m_int;
                            magn[i] = 0;
                            nzposcurr--;
                        }
                    }
                }
            }

            pulsescurr = pulsesnum;

            /* Magnitudes coding */
            {
                leftp  = pulsescurr;/*pulsesnum; */
                leftnz = nzposcurr; /*nzpos; */
                trellislevel = 0;

                for( i = 0; i < size && leftnz > 1; i++)
                {
                    if( magn[i] != 0.0f )
                    {
                        for( j = 0; j < leftp; j++)
                        {
                            prob1 = (leftnz - 1.0f) / (leftp - j - 1.0f);
                            prob0 = 1.0f - prob1;

                            st = savedstates[ trellislevel ];
                            quantum1 = (float)quantize( (float)(j + 1), ddec[st][ 0 ]);
                            quantum2 = (float)quantize( (float)(j + 1), ddec[st][ 1 ]);

                            if( quantum1 != (j+1) && quantum2 != (j+1) )
                            {
                                /* this magnitude is not possible so set probabilities */
                                prob0 = 1.0f;
                                prob1 = 0.0f;
                            }

                            if( j < fabs(magn[i]) - 1 )
                            {
                                magnbits -= log2_f( prob0 );
                            }
                            else
                            {
                                magnbits -= log2_f( prob1 );
                                break;
                            }
                        }

                        leftnz--;
                        leftp -= (int)fabs( magn[i] );
                    }

                    trellislevel++;
                }

                /* Update actual occured surplus */
                *nzpout = nzposcurr;
            }
        }

        if( *nzpout > 1 && bcount != 0 )
        {
            flag_g1 = 0;

            for( i = 0; i < size; i++)
            {
                if( fabs(magn[i]) > 1.0f )
                {
                    flag_g1 = 1;
                }
            }
            /* prepare vector for TCQ */
            for( i = 0; i < size && flag_g1 && *bcount < 2*TCQ_AMP; i++)
            {
                if( fabs(magn[i]) > 0.0f )
                {
                    abuffer[*bcount] = quants[i];
                    mbuffer[*bcount] = magn[i];
                    sbuffer[*bcount] = scale;

                    (*bcount)++;

                }
            }
        }

        if( actualt > 0)
        {
            *qscale = pulsesnum / actualt;
        }
        else
        {
            *qscale = FLT_MAX;
        }
    }

    return magnbits;
}


void InitLSBTCQ(short *bcount)
{
    *bcount = 0;

    return;
}


void TCQLSB(short bcount, float *abuffer, float *mbuffer, float *sbuffer, short *dpath)

{
    short i, st, dminpos, position;
    float q = QTCQ;
    float dmin, curdist1, curdist2, newdist1, newdist2;
    float metric[STATES_LSB][TCQ_LSB_SIZE];
    short path[STATES_LSB][TCQ_LSB_SIZE];
    short quant[STATES_LSB][TCQ_LSB_SIZE];
    short dquant[STATES_LSB][TCQ_LSB_SIZE];
    short qout[TCQ_LSB_SIZE];
    float q1, q2, s1, s2, a1, a2, sign1, sign2;
    float dbuffer[MAX_PULSES];

    set_f( *metric, 0.0f, STATES_LSB*TCQ_LSB_SIZE );
    set_s( *path, 0, STATES_LSB*TCQ_LSB_SIZE );
    set_s( *quant, 0, STATES_LSB*TCQ_LSB_SIZE );
    set_s( *dquant, 0, STATES_LSB*TCQ_LSB_SIZE );
    set_s( qout, 0, TCQ_LSB_SIZE );
    set_f( dbuffer, 0.0f, MAX_PULSES );

    metric[1][0] = 16777216.0f;
    metric[2][0] = 16777216.0f;
    metric[3][0] = 16777216.0f;

    for( i = 0; i < 2*TCQ_AMP; i+=2 )
    {
        q1 = mbuffer[i];
        q2 = mbuffer[i + 1];

        s1 = sbuffer[i];
        s2 = sbuffer[i + 1];

        a1 = abuffer[i];
        a2 = abuffer[i + 1];

        /* cycle over conditions */
        for( st = 0; st < 4; st++)
        {
            curdist1 = metric[ step_LSB[st][0] ][i/2];
            curdist2 = metric[ step_LSB[st][1] ][i/2];

            /* step 1 */
            sign1 = (denc_LSB[st][0] & 0x1)?(q):(-q);
            sign2 = (denc_LSB[st][0] & 0x2)?(q):(-q);
            newdist1 = (a1 - (q1+sign1)/s1)*(a1 - (q1+sign1)/s1) + (a2 - (q2+sign2)/s2)*(a2 - (q2+sign2)/s2);

            /* step 2 */
            sign1 = (denc_LSB[st][1] & 0x1)?(q):(-q);
            sign2 = (denc_LSB[st][1] & 0x2)?(q):(-q);
            newdist2 = (a1 - (q1+sign1)/s1)*(a1 - (q1+sign1)/s1) + (a2 - (q2+sign2)/s2)*(a2 - (q2+sign2)/s2);

            /* decision */
            if( curdist1 + newdist1 < curdist2 + newdist2 )
            {
                path[st][i/2+1] = step_LSB[st][0];
                metric[st][i/2+1] = curdist1 + newdist1;
                quant[st][i/2+1] = 0;
                dquant[st][i/2+1] = dqnt_LSB[ step_LSB[st][0] ][st];
            }
            else
            {
                path[st][i/2+1] = step_LSB[st][1];
                metric[st][i/2+1] = curdist2 + newdist2;
                quant[st][i/2+1] = 1;
                dquant[st][i/2+1] = dqnt_LSB[ step_LSB[st][0] ][st];
            }
        }

    }

    /* Find path with minimal metric */
    dminpos = 0;
    dmin    = metric[ dminpos][ i/2];
    for( st = 1; st < 4; st++)
    {
        if( dmin > metric[ st][ i/2] )
        {
            dmin = metric[ st][ i/2];
            dminpos = st;
        }
    }

    /* Trace back to get output */
    position = dminpos;

    for( ; i >= 0; i-=2)
    {
        qout[i/2] = quant[position][ i/2+1 ];
        dpath[i/2] = dquant[position][ i/2+1 ];

        dbuffer[i]     = (denc_LSB[position][qout[i/2]] & 0x1)?(q):(-q);
        dbuffer[i + 1] = (denc_LSB[position][qout[i/2]] & 0x2)?(q):(-q);

        position  = path[position][i/2+1];
    }

    /* add decoded sequence to quanta */
    for( i = 0; i < bcount; i++ )
    {
        mbuffer[i] += dbuffer[i];
    }

    return;
}

void TCQLSBdec( short *dpath, float *mbuffer, short bcount )
{
    float q = QTCQ;
    int i, state = 0;

    for( i = 0; i < bcount/2; i++)
    {
        mbuffer[2*i]     = ( ddec_LSB[state][dpath[i]] & 0x1)?(q):(-q);
        mbuffer[2*i + 1] = ( ddec_LSB[state][dpath[i]] & 0x2)?(q):(-q);

        state  = dstep_LSB[state][dpath[i]];
    }

    return;
}

void SaveTCQdata( PARCODEC arInst, short *dpath, short bcount)
{
    int i;
    for( i = 0; i < bcount; i++)
    {
        ar_encode_uniform( arInst, dpath[i], 1);
    }

    return;
}

void LoadTCQdata( PARCODEC arInst, short *dpath, short bcount)
{
    int i;
    for( i = 0; i < bcount; i++)
    {
        dpath[i] = ar_decode( arInst, uniform_model );
    }

    return;
}

void RestoreTCQdec( int * magn, int size, short *bcount, float *mbuffer)
{
    int i, nzpos = 0, flag_g1 = 0;

    /* calculate actual nz positions */
    for( i = 0, nzpos = 0; i < size; i++)
    {
        if( magn[i] != 0 )
        {
            nzpos++;
            if( abs(magn[i]) > 1 )
            {
                flag_g1 = 1;
            }
            magn[i] *= (1.0f/QTCQ);
        }
    }

    if( nzpos > 1)
    {
        for( i = 0; i < size && flag_g1 && *bcount < 2*TCQ_AMP; i++)
        {
            if( magn[i] != 0 )
            {
                mbuffer[*bcount] = magn[i] + (1.0f/QTCQ)*mbuffer[*bcount];
                magn[i] = round_f( mbuffer[*bcount] );
                (*bcount)++;
            }
        }
    }

    return;
}


void RestoreTCQ( float * magn, int size, short *bcount, float *mbuffer)
{
    int i, nzpos = 0, flag_g1 = 0;

    /* calculate actual nz positions */
    for( i = 0, nzpos = 0; i < size; i++)
    {
        if( magn[i] != 0.0f )
        {
            nzpos++;
            if( fabs(magn[i]) > 1.0f )
            {
                flag_g1 = 1;
            }
        }
    }

    if( nzpos > 1)
    {
        for( i = 0; i < size && flag_g1 && *bcount < 2*TCQ_AMP; i++)
        {
            if( fabs(magn[i]) > 0.0f )
            {
                magn[i] = mbuffer[*bcount];
                (*bcount)++;
            }
        }
    }

    return;
}

