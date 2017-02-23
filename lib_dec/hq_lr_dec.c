/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "rom_dec.h"
#include "prot.h"
#include "stl.h"
#include "basop_util.h"

/*--------------------------------------------------------------------------*
 * Local functions
 *--------------------------------------------------------------------------*/

static short p2a_threshold_dequant( Decoder_State *st, short *p2a_flags, const short bands, const short p2a_bands );

static void mdct_spectrum_fine_gain_dec( Decoder_State *st, float y2[], const short band_start[], const short band_end[],
        const short k_sort[], const short bands, const Word32 L_qint,
        const short Ngq, const short gqlevs, const short gqbits );

static float band_energy_dequant( Decoder_State *st, float band_energy[], const short bands, const Word32 L_qint, const Word16 eref_fx, const short is_transient );

static void spt_shorten_domain_set_dec( Decoder_State *st, const short p2a_flags[], const short new_band_start[],
                                        const short new_band_end[], const short new_band_width[], const short bands,
                                        short band_start[], short band_end[], short band_width[], short *bit_budget );

/*-------------------------------------------------------------------*
 * hq_lr_dec()
 *
 * HQ low rate decoding routine
 *-------------------------------------------------------------------*/

void hq_lr_dec(
    Decoder_State *st,            /* i/o: decoder state structure         */
    float yout[],         /* o  : transform-domain output coefs.  */
    const short inner_frame,    /* i  : inner frame length              */
    short num_bits,       /* i  : number of available bits        */
    short *is_transient   /* o  : transient flag                  */
)
{
    short i, k1, pbits, p2a_flags[BANDS_MAX], bit_budget, bands, length, gqlevs, gqbits, Ngq, p2a_bands, ni_seed;
    short band_start[BANDS_MAX], band_end[BANDS_MAX], band_width[BANDS_MAX];
    short k_sort[BANDS_MAX];
    int npulses[BANDS_MAX], inp_vector[L_FRAME48k];
    float ni_coef, ni_pd_th, pd_thresh, ld_slope;
    float ebits, Rk[BANDS_MAX], band_energy[BANDS_MAX], y2[L_FRAME48k], p2a_th;
    Word32 Rk_fx[BANDS_MAX];
    Word32 L_qint;              /* Q29 */
    Word16 eref_fx;             /* Q10 */
    Word16 bit_alloc_weight_fx; /* Q13 */
    float y2_ni[L_FRAME48k],y2_org[L_FRAME48k];
    short hqswb_clas = 0;
    short lowlength, highlength, har_bands = 0;
    float m[L_FRAME32k];
    float Ep[BANDS_MAX], enerH = 0.0f, enerL = 0.0f;
    short lowband, highband, bw_low = 0, bw_high = 20;
    float Ep_tmp[BANDS_MAX];
    float band_energy_tmp[BANDS_MAX];
    short last_bitalloc_max_band[2];
    long bwe_br;
    short trans_bit, p2a_flags_tmp[BANDS_MAX];
    short adjustFlag = 0;
    short prev_SWB_peak_pos_tmp[SPT_SHORTEN_SBNUM];
    int   j, k;
    short flag_spt;
    short org_band_start[SPT_SHORTEN_SBNUM];
    short org_band_end[SPT_SHORTEN_SBNUM];
    short org_band_width[SPT_SHORTEN_SBNUM];
    short new_band_start[SPT_SHORTEN_SBNUM];
    short new_band_end[SPT_SHORTEN_SBNUM];
    short new_band_width[SPT_SHORTEN_SBNUM];
    Word32 L_tmp,L_tmp2,L_tmp3;
    Word16 exp,tmp,exp2,tmp1,tmp2,tmp3,alpha_fx,frac1;
    Word32 enerH_fx;
    Word32 enerL_fx;
    Word32 Ep_fx[BANDS_MAX];
    Word32 Ep_avrg_fx, Ep_vari_fx;
    Word32 Ep_avrgL_fx;
    Word32 Ep_peak_fx;
    Word32 Ep_tmp_fx[BANDS_MAX];
    Word16 gama_fx;/*Q15 0.85f; */
    Word16 beta_fx;/*Q14 1.05f; */
    Word32 L_band_energy[BANDS_MAX],L_band_energy_tmp[BANDS_MAX];
    UWord16 lo;
    Word16 Q_band_energy;

    set_s(last_bitalloc_max_band, 0, 2);
    set_f( y2, 0.0f, L_FRAME48k );
    set_i( inp_vector, 0, inner_frame );
    flag_spt = 0;
    set_s(prev_SWB_peak_pos_tmp, 0, SPT_SHORTEN_SBNUM);
    bwe_br = st->core_brate;

    if( st->bwidth == SWB && ( bwe_br == HQ_16k40 || bwe_br == HQ_13k20 ) )
    {
        hqswb_clas = (short)get_next_indice( st,2);
        num_bits -= 2;

        *is_transient = 0;
        if ( hqswb_clas == HQ_TRANSIENT )
        {
            *is_transient = 1;
        }
    }
    else
    {
        /* decode transient flag */
        *is_transient = (short)get_next_indice( st, 1 );
        num_bits--;
    }


    /* Configure decoder for different bandwidths, bit rates, etc. */
    hq2_core_configure( inner_frame, num_bits, *is_transient, &bands, &length, band_width, band_start, band_end, &L_qint, &eref_fx,
                        &bit_alloc_weight_fx, &gqlevs, &Ngq, &p2a_bands, &p2a_th, &pd_thresh, &ld_slope, &ni_coef, &ni_pd_th, bwe_br );

    highlength = band_end[bands-1];
    har_bands = bands;

    if( st->bwidth == SWB  && *is_transient == 0 && (bwe_br == HQ_16k40 || bwe_br == HQ_13k20) )
    {
        /* reserve bits for HQ_NORMAL and HQ_HARMONIC */
        if( hqswb_clas == HQ_NORMAL || hqswb_clas==HQ_HARMONIC)
        {
            num_bits -= (short)get_usebit_npswb( hqswb_clas );
        }

        if( hqswb_clas == HQ_NORMAL )
        {
            flag_spt = 1;
        }
    }

    if( (bwe_br == HQ_16k40 || bwe_br == HQ_13k20) && st->bwidth == SWB )
    {
        if( st->prev_hqswb_clas != HQ_NORMAL )
        {
            j = 0;
            for(k=bands-SPT_SHORTEN_SBNUM; k<bands; k++)
            {
                st->prev_SWB_peak_pos[j] = 0;
                j++;
            }
        }
    }

    /* Spectral energy calculation/quantization */
    ebits = band_energy_dequant( st, band_energy, bands, L_qint, eref_fx, *is_transient );

    /* simple check: band_energy is too large, Abnormal Situation of bit errors */
    for( k=0; k<bands; k++ )
    {
        /* Max: 45.0(737279,Q14) at 32kHz, highest band, Min: -6.600037(-108135,Q14) at 8kHz(NB),8kbps, is_transient (-6.7f is safty-threshold) */
        if( band_energy[k] > 45.0f || band_energy[k] < -6.7f )
        {
            st->BER_detect = 1;
            set_f( yout, 0, inner_frame );
            return;
        }
    }

    for (i = 0; i < bands; i++)
    {
        L_band_energy[i] = (Word32)(band_energy[i] * pow(2.0f, SWB_BWE_LR_Qbe));
    }

    /* First pass bit budget for TCQ of spectral band information */
    gqbits = (short int) log2_f ((float) gqlevs);
    bit_budget = num_bits - (short) ceil (ebits) - Ngq * gqbits;

    pbits = 0;
    if( st->bwidth == SWB && (bwe_br == HQ_16k40 || bwe_br == HQ_13k20) )
    {
        if ( hqswb_clas == HQ_HARMONIC)
        {
            set_s( p2a_flags, 1, har_bands );
        }
        else
        {
            /* High band tonality detector based on per band peak-to-average ratio */
            pbits = p2a_threshold_dequant( st, p2a_flags, bands, p2a_bands );
            bit_budget -= pbits;

            if( hqswb_clas == HQ_NORMAL )
            {
                return_bits_normal2( &bit_budget, p2a_flags, bands, bits_lagIndices );
            }
        }
    }
    else
    {
        /* High band tonality detector based on per band peak-to-average ratio */
        pbits = p2a_threshold_dequant( st, p2a_flags, bands, p2a_bands );
        bit_budget -= pbits;
    }

    if( flag_spt == 1 )
    {
        /* initalize the desired parameters for SPT */
        spt_shorten_domain_band_save( bands, band_start, band_end, band_width, org_band_start, org_band_end, org_band_width );
        spt_shorten_domain_pre( band_start, band_end, st->prev_SWB_peak_pos, bands, bwe_br, new_band_start, new_band_end, new_band_width );
        spt_shorten_domain_set_dec( st, p2a_flags, new_band_start, new_band_end, new_band_width, bands, band_start, band_end, band_width, &bit_budget );
    }

    /* safety check in case of bit errors */
    if( bit_budget < 2 )
    {
        st->BER_detect = 1;
        set_f( yout, 0, inner_frame );
        return;
    }

    /* Estimate number of bits per sub-band */
    Q_band_energy = SWB_BWE_LR_Qbe;
    FOR(i = 0; i < bands; i++)
    {
        L_tmp = L_shl(L_band_energy[i],sub(16,Q_band_energy));/*Q16 */

        frac1 = L_Extract_lc(L_tmp, &exp); /* Extract exponent of L_tmp */
        L_tmp =  Pow2(30, frac1);
        exp = sub(exp, 30);
        Ep_fx[i] = L_shl(L_tmp , sub(exp,6)); /* Q -6 */
        Ep[i] = (float)(Ep_fx[i]/pow(2.0,-6));
    }

    FOR( i = 0; i < bands; i++ )
    {
        L_tmp2 = Ep_fx[i];
        L_tmp = L_max(1, L_tmp2);
        exp = norm_l(L_tmp);
        tmp = extract_h(L_shl(L_tmp, exp));

        L_tmp3 = (Word32)band_width[i];
        exp2 = norm_l(L_tmp3);
        tmp2 = extract_h(L_shl(L_tmp3, exp2));

        exp2 = sub(exp, exp2); /* Denormalize and substract */

        tmp3 = sub(tmp2, tmp);
        IF (tmp3 > 0)
        {
            tmp2 = shr(tmp2, 1);
        }
        IF (tmp3 > 0)
        {
            exp2 = add(exp2, 1);
        }
        tmp = div_s(tmp2, tmp);
        L_tmp = L_deposit_h(tmp);
        L_tmp = Isqrt_lc1(L_tmp, &exp2);
        move32();/*Q(31-exp2) */
        Ep_tmp_fx[i] = L_shr(L_tmp,sub(15,exp2));/*Q13 */
        Ep_tmp[i] = (float)(Ep_tmp_fx[i]/pow(2.0,13));
    }

    if ( *is_transient == 0 && inner_frame == L_FRAME8k && st->core_brate <= ACELP_13k20)
    {
        /* decode the last p2a_bands-1 subbands bit-allocation index of the previous frame */
        j = 0;
        for(i = 0; i < 2; i++)
        {
            last_bitalloc_max_band[i] = (short)get_next_indice( st, 1 );
        }

        lowband = 6;
        move16();
        trans_bit = 2;
        move16();
        bit_budget =sub(bit_budget,trans_bit);
        gama_fx = 27852; /*Q15 0.85f;*/
        beta_fx = 17203;
        move16();/*Q14 1.05f; */
        set_s( &p2a_flags_tmp[bands-trans_bit], 0, 2 );

        IF( st->core_brate == ACELP_13k20 )
        {
            beta_fx = 13107;
            move16();/*14 1.25f; */
            gama_fx = 31130;
            move16();/*0.95f; */
            mvs2s(&p2a_flags[sub(bands,trans_bit)], &p2a_flags_tmp[sub(bands,trans_bit)], trans_bit);
        }

        /* calculate the the low band/high band energy and the variance/avrage of the envelopes */
        Ep_vari_fx = 0;
        move32();
        Ep_avrg_fx = 0;
        move32();
        Ep_avrgL_fx = 0;
        move32();
        Ep_peak_fx = 0;
        move32();
        FOR( i = 0; i < bands; i++ )
        {
            IF( sub(i,lowband) >= 0)
            {
                Ep_vari_fx = L_add(Ep_vari_fx,L_abs(L_sub(Ep_tmp_fx[i],Ep_tmp_fx[sub(i,1)])));/*Q15 */
                Ep_avrg_fx = L_add(Ep_avrg_fx,Ep_tmp_fx[i]);/*Q15 */

            }
            ELSE
            {
                Ep_avrgL_fx = L_add(Ep_avrgL_fx,Ep_tmp_fx[i]);/*Q15 */
                IF(L_sub(Ep_tmp_fx[i],Ep_peak_fx) > 0)
                {
                    Ep_peak_fx = Ep_tmp_fx[i];
                    move32();/*Q15 */
                }
            }
        }
        /* modify the last p2a_bands subbands band_energies */
        k = (int)bands;
        mvi2i( L_band_energy,L_band_energy_tmp,k); /*Q_band_energy */
        Mpy_32_16_ss(Ep_peak_fx,24576,&L_tmp,&lo);
        Mpy_32_16_ss(Ep_peak_fx,shl(sub(bands,lowband),9),&L_tmp2,&lo);
        Mpy_32_16_ss(Ep_avrg_fx,1126,&L_tmp3,&lo);

        IF(( (L_sub(L_tmp, L_shr(Ep_avrgL_fx,1)) < 0  && st->core_brate == ACELP_13k20 ) || st->core_brate < ACELP_13k20 )&&
           L_sub(L_tmp2, L_tmp3) < 0 && L_sub(L_tmp2, L_shr(Ep_avrg_fx,7)) > 0)
        {
            FOR(i = lowband; i < bands; i++)
            {
                Mpy_32_16_ss(Ep_avrg_fx,24576,&L_tmp,&lo);
                IF(L_sub(L_shr(Ep_tmp_fx[i],1), L_tmp) < 0)
                {
                    Mpy_32_16_ss(Ep_peak_fx,sub(bands,lowband),&L_tmp,&lo);
                    tmp = extract_h(L_shl(L_tmp,14));/*Q-4 */
                    IF(tmp != 0)
                    {
                        exp = norm_s(tmp);
                        tmp = shl(tmp,exp);/*Q(exp) */
                        tmp = div_s(16384,tmp);/*Q(15+14-exp=29-exp) */
                        exp = sub(29,exp);
                    }
                    ELSE
                    {
                        tmp = 0x7fff;
                        move16();
                        exp = 0;
                        move16();
                    }
                    Mpy_32_16_ss(Ep_avrg_fx,tmp,&L_tmp,&lo);
                    L_tmp = L_shl(L_tmp,sub(13,exp));/*Q(13+exp-15 +13-exp +4 = 15) */
                    L_tmp2 = L_add(L_tmp,13107); /*15 */
                    tmp2 = extract_l(L_min(L_max(L_tmp2,16384),gama_fx)); /*15 = 15 */
                    Mpy_32_16_ss(L_band_energy_tmp[i],tmp2,&L_band_energy_tmp[i],&lo);
                }
            }
        }
        ELSE
        {
            FOR(i = sub(bands,trans_bit); i < bands; i++)
            {
                alpha_fx = 16384;
                move16();/*Q14 */
                IF( sub(p2a_flags_tmp[i],1) == 0)
                {
                    Mpy_32_16_ss(Ep_tmp_fx[i],sub(bands,lowband),&L_tmp,&lo);
                    tmp = extract_h(L_shl(L_tmp,14));/*Q-4 */
                    IF(tmp != 0)
                    {
                        exp = norm_s(tmp);
                        tmp = shl(tmp,exp);/*Q(exp) */
                        tmp = div_s(16384,tmp);/*Q(15+14-exp=29-exp) */
                        exp = sub(29,exp);
                    }
                    ELSE
                    {
                        tmp = 0x7fff;
                        move16();
                        exp = 0;
                        move16();
                    }
                    Mpy_32_16_ss(Ep_vari_fx,3277,&L_tmp,&lo);
                    Mpy_32_16_ss(L_tmp,tmp,&L_tmp,&lo);
                    L_tmp = L_shl(L_tmp,sub(12,exp));/*Q(13+exp-15 +12-exp +4 = 14) */

                    tmp2 = extract_h(Ep_avrg_fx);/*Q13-16=-3 */
                    IF(tmp2 != 0)
                    {
                        exp = norm_s(tmp2);
                        tmp2 = shl(tmp2,exp);/*Q(exp) */
                        tmp2 = div_s(16384,tmp2);/*Q(15+14-exp=29-exp) */
                        exp = sub(29,exp);
                    }
                    ELSE
                    {
                        /*when the divisor is zero, happens rarely*/
                        tmp2 = 0x7fff;
                        move16();
                        exp = 0;
                        move16();
                    }
                    Mpy_32_16_ss(Ep_vari_fx,6554,&L_tmp2,&lo);
                    Mpy_32_16_ss(L_tmp2,tmp2,&L_tmp2,&lo);
                    L_tmp2 = L_shl(L_tmp2,sub(13,exp));/*Q(13+exp-15 +13-exp +3 = 14) */
                    L_tmp=L_min(L_tmp,L_tmp2);/*14 */
                    tmp=extract_l(L_min(L_tmp,13107));/*14 */
                    alpha_fx =add(16384,tmp);

                }
                IF(sub(last_bitalloc_max_band[j++], 1) == 0)
                {
                    Mpy_32_16_ss(Ep_tmp_fx[i],sub(bands,lowband),&L_tmp,&lo);
                    tmp = extract_h(L_shl(L_tmp,14));/*Q-2 */
                    IF(tmp != 0)
                    {
                        exp = norm_s(tmp);
                        tmp = shl(tmp,exp);/*Q(exp) */
                        tmp = div_s(16384,tmp);/*Q(15+14-exp=29-exp) */
                        exp = sub(29,exp);
                    }
                    ELSE
                    {
                        tmp = 0x7fff;
                        move16();
                        exp = 0;
                        move16();
                    }
                    Mpy_32_16_ss(Ep_avrg_fx,tmp,&L_tmp,&lo);
                    L_tmp = L_shl(L_tmp,sub(14,exp));/*Q(13+exp-15 +14-exp+2 = 14) */
                    L_tmp =L_max(L_tmp,16384); /*14 */
                    tmp=extract_l(L_min(L_tmp,beta_fx)); /*14 */
                    alpha_fx=shl(mult(alpha_fx,tmp),1);/*14+14-15 +1=14 */
                }
                ELSE
                {
                    tmp2 = extract_h(Ep_avrg_fx);/*13 -16 =-3 */
                    IF(tmp2 != 0)
                    {
                        exp = norm_s(tmp2);
                        tmp2 = shl(tmp2,exp);/*Q(exp) */
                        tmp2 = div_s(16384,tmp2);/*Q(15+14-exp=29-exp) */
                        exp = sub(29,exp);
                    }
                    ELSE
                    {
                        /*when the divisor is zero, happens rarely*/
                        tmp2 = 0x7fff;
                        move16();
                        exp = 0;
                        move16();
                    }
                    Mpy_32_16_ss(Ep_tmp_fx[i],tmp2,&L_tmp,&lo);
                    L_tmp = L_shl(L_tmp,sub(19,exp));/*Q(13+exp-15 +19-exp +3 = 20) */
                    Mpy_32_16_ss(L_tmp,shl(sub(bands,lowband),9),&L_tmp,&lo);
                    L_tmp =L_max(L_tmp,13926); /*14 */
                    tmp2 =extract_l(L_min(L_tmp,16384)); /*14 */
                    alpha_fx=shl(mult(alpha_fx,tmp2),1);/*14+14-15+1 =14 */
                }
                Mpy_32_16_ss(L_band_energy_tmp[i],alpha_fx,&L_tmp,&lo);
                L_band_energy_tmp[i] = L_shl(L_tmp,1);/*Q(Q_band_energy+14-15 +1= Q_band_energy) */
            }
        }
        lowband = 3;
        move16();
        Ep_avrg_fx = 0;
        move32();
        Ep_avrgL_fx = 0;
        move32();
        Ep_peak_fx = 0;
        move32();
        FOR(i = 0; i < bands; i++)
        {
            IF(sub(i,lowband) >=0 )
            {
                Ep_avrg_fx = L_add(Ep_avrg_fx,Ep_tmp_fx[i]);/*Q15 */
            }
            ELSE
            {
                Ep_avrgL_fx = L_add(Ep_avrgL_fx,L_shr(Ep_tmp_fx[i],1));/*Q12 */
                IF(L_sub(Ep_tmp_fx[i],Ep_peak_fx) > 0)
                {
                    Ep_peak_fx = Ep_tmp_fx[i];
                    move32();/*Q13 */
                }
            }
        }
        Mpy_32_16_ss(Ep_peak_fx,28262,&L_tmp,&lo);
        Mpy_32_16_ss(Ep_avrgL_fx,24576,&L_tmp2,&lo);
        IF( L_sub(L_shr(Ep_avrg_fx,2), L_tmp2) > 0 && L_sub(L_shr(Ep_avrg_fx,4), L_tmp2) < 0 &&  L_sub(L_tmp, Ep_avrgL_fx)>0)
        {
            adjustFlag = 1;
            move16();
            FOR (i = 0; i < lowband; i++)
            {
                tmp = extract_h(Ep_avrgL_fx);/*Q-4 */
                IF(tmp != 0)
                {
                    exp = norm_s(tmp);
                    tmp = shl(tmp,exp);/*Q(exp) */
                    tmp = div_s(16384,tmp);/*Q(15+14-exp=29-exp) */
                    exp = sub(29,exp);
                }
                ELSE
                {
                    tmp = 0x7fff;
                    move16();
                    exp = 0;
                    move16();
                }
                Mpy_32_16_ss(Ep_peak_fx,tmp,&L_tmp,&lo);
                Mpy_32_16_ss(L_tmp,lowband,&L_tmp,&lo);
                Mpy_32_16_ss(L_tmp,18842,&L_tmp,&lo);
                L_tmp = L_shl(L_tmp,sub(27,exp));/*Q14 0.5 */
                tmp2=extract_l(L_min(L_tmp,19661));/*14 */
                Mpy_32_16_ss(L_band_energy_tmp[i],tmp2,&L_tmp,&lo);
                L_band_energy_tmp[i] = L_shl(L_tmp,1); /*Q_band_energy  */
            }
        }
        for (i = 0; i < bands; i++)
        {
            band_energy_tmp[i] = (float)(L_band_energy_tmp[i]/pow(2.0f, SWB_BWE_LR_Qbe));
        }

        hq2_bit_alloc( band_energy_tmp, bands, Rk_fx, &bit_budget, p2a_flags, bit_alloc_weight_fx,
                       band_width, num_bits, hqswb_clas, st->bwidth, *is_transient );
    }
    else if( *is_transient == 0 && inner_frame == L_FRAME16k )
    {
        bit_budget = sub(bit_budget,2);/* bits in high bands to indicate the last 2 subbands is allocated bits or not */
        /* decode the last p2a_bands-1 subbands bit-allocation index of the previous frame */
        for(i = 0; i < 2; i++)
        {
            last_bitalloc_max_band[i] = (short)get_next_indice( st, 1 );
        }
        FOR( i = 0; i < bands; i++ )
        {
            Ep_tmp_fx[i] = L_shl(Ep_tmp_fx[i],2);
        }
        IF( st->core_brate == ACELP_13k20 )
        {
            lowband = 8;
            move16();
            highband = 15;
            move16();
            bw_low = sub(band_start[highband],band_start[lowband]);
            bw_high = sub(add(band_end[sub(bands,1)],1),band_start[highband]);
        }
        ELSE
        {
            lowband = 8;
            move16();
            highband = 16;
            move16();
            bw_low = sub(band_start[highband],band_start[lowband]);
            bw_high = sub(add(band_end[sub(bands,1)],1),band_start[highband]);
        }
        /* calculate the the low band/high band energy and the variance/avrage of the envelopes */
        enerL_fx = 0;
        move32();
        enerH_fx = 0;
        move32();
        Ep_vari_fx = 0;
        move32();
        Ep_avrg_fx = 0;
        move32();
        FOR( i = 0; i < bands; i++ )
        {
            IF( sub(i,lowband) >= 0 && add(sub(i,bands),p2a_bands) < 0)
            {
                Ep_vari_fx = L_add(Ep_vari_fx,L_abs(L_sub(Ep_tmp_fx[i],Ep_tmp_fx[sub(i,1)])));/*Q15 */
                Ep_avrg_fx = L_add(Ep_avrg_fx,Ep_tmp_fx[i]);/*Q15 */
            }

            IF(sub(i,highband) >= 0)
            {
                enerH_fx  = L_add(enerH_fx,L_shl(Ep_fx[i],2));/*Q0 */
            }
            ELSE IF(sub(i,lowband) >= 0)
            {
                enerL_fx  = L_add(enerL_fx,L_shl(Ep_fx[i],2));/*Q0 */
            }
        }

        enerL = (float)(enerL_fx/pow(2.0,-4));
        enerH = (float)(enerH_fx/pow(2.0,-4));
        /* modify the last p2a_bands subbands band_energies */
        k = (int)bands;
        mvi2i( L_band_energy,L_band_energy_tmp,k); /*Q_band_energy */

        L_tmp = L_max(enerH_fx,enerL_fx);
        tmp = s_max(bw_low,bw_high);
        i = norm_l(L_tmp);
        j = norm_s(tmp);
        Mpy_32_16_ss(L_shl(enerH_fx,i),shl(bw_low,j),&L_tmp,&lo);
        Mpy_32_16_ss(L_shl(enerL_fx,i),shl(bw_high,j),&L_tmp2,&lo);
        L_tmp2 = L_sub(L_tmp,L_tmp2);

        FOR( i = sub(bands,p2a_bands); i < bands; i++ )
        {
            IF( sub(p2a_flags[i],1) == 0 || L_tmp2 > 0 )
            {
                tmp = sub(bands,p2a_bands);
                tmp = sub(tmp,lowband);/*Q0 */

                tmp1 = extract_h(L_shl(Ep_avrg_fx,1));/*Q0 */
                IF(tmp1 != 0)
                {
                    exp = norm_s(tmp1);
                    tmp1 = shl(tmp1,exp);/*Q(exp) */
                    tmp1 = div_s(16384,tmp1);/*Q(15+14-exp = 29-exp) */
                    exp = sub(29,exp);
                }
                ELSE
                {
                    /*when the divisor is zero, happens rarely*/
                    tmp1 = 0x7fff;
                    move16();
                    exp = 0;
                    move16();
                }
                Mpy_32_16_ss(Ep_tmp_fx[i],tmp1,&L_tmp,&lo);
                Mpy_32_16_ss(L_tmp,tmp,&L_tmp,&lo);
                Mpy_32_16_ss(L_tmp,16384,&L_tmp,&lo);
                L_tmp = L_shl(L_tmp,sub(32,exp));/*Q15 */
                tmp = extract_l(L_min(L_tmp,6554));/*Q15 */
                Mpy_32_16_ss(Ep_vari_fx,tmp1,&L_tmp,&lo);
                Mpy_32_16_ss(L_tmp,tmp,&L_tmp,&lo);
                L_tmp = L_shl(L_tmp,sub(15,exp));/*Q15 */
                tmp = extract_l(L_shr(L_min(L_tmp,13107),1));/*Q14 */
                alpha_fx = add(tmp,16384);/*Q14 */
            }
            ELSE
            {
                alpha_fx = 16384;
                move16();/*Q14 */
            }

            IF(add(sub(i,bands),p2a_bands) > 0)
            {
                tmp = sub(bands, p2a_bands);
                IF(sub(last_bitalloc_max_band[sub(i, add(tmp, 1))], 1) == 0)
                {
                    tmp = sub(tmp,lowband);
                    Mpy_32_16_ss(Ep_tmp_fx[i],tmp,&L_tmp,&lo);
                    tmp = extract_h(L_shl(L_tmp,16));/*Q0 */
                    IF(tmp != 0)
                    {
                        exp = norm_s(tmp);
                        tmp = shl(tmp,exp);/*Q(exp) */
                        tmp = div_s(16384,tmp);/*Q(15+14-exp=29-exp) */
                        exp = sub(29,exp);
                    }
                    ELSE
                    {
                        tmp = 0x7fff;
                        move16();
                        exp = 0;
                        move16();
                    }
                    Mpy_32_16_ss(Ep_avrg_fx,tmp,&L_tmp,&lo);
                    L_tmp = L_shl(L_tmp,sub(14,exp));/*Q14 */
                    tmp = extract_l(L_min(L_max(L_tmp,16384),20480));/*Q14 */
                    L_tmp = L_mult(alpha_fx,tmp);/*Q(14+14+1=29) */
                    alpha_fx = extract_l(L_shr(L_tmp,15)); /*Q14                    */
                }
                ELSE
                {
                    tmp = sub(tmp,lowband);

                    tmp1 = extract_h(L_shl(Ep_avrg_fx,1));/*Q0 */
                    IF(tmp1 != 0)
                    {
                        exp = norm_s(tmp1);
                        tmp1 = shl(tmp1,exp);/*Q(exp) */
                        tmp1 = div_s(16384,tmp1);/*Q(15+14-exp=29-exp) */
                        exp = sub(29,exp);
                    }
                    ELSE
                    {
                        /*when the divisor is zero, happens rarely*/
                        tmp1 = 0x7fff;
                        move16();
                        exp = 0;
                        move16();
                    }
                    Mpy_32_16_ss(Ep_tmp_fx[i],tmp1,&L_tmp,&lo);
                    Mpy_32_16_ss(L_tmp,tmp,&L_tmp,&lo);
                    L_tmp = L_shl(L_tmp,sub(29,exp));/*Q14 */
                    tmp = extract_l(L_min(L_max(L_tmp,13926),16384));/*Q14 */
                    L_tmp = L_mult(alpha_fx,tmp);/*Q(14+14+1=29) */
                    alpha_fx = extract_l(L_shr(L_tmp,15)); /*Q14  */
                }
            }
            Mpy_32_16_ss(L_band_energy_tmp[i],alpha_fx,&L_tmp,&lo);
            L_band_energy_tmp[i] = L_shl(L_tmp,1);/*Q Q_band_energy */
        }
        lowband = 6;
        move16();
        Ep_avrg_fx = 0;
        move32();
        Ep_avrgL_fx = 0;
        move32();
        Ep_peak_fx = 0;
        move32();
        FOR(i = 0; i < bands; i++)
        {
            IF(sub(i,lowband) >= 0)
            {
                Ep_avrg_fx = L_add(Ep_avrg_fx,Ep_tmp_fx[i]);/*Q15 */
            }
            ELSE
            {
                Ep_avrgL_fx = L_add(Ep_avrgL_fx,Ep_tmp_fx[i]);/*Q15 */
                IF(L_sub(Ep_tmp_fx[i],Ep_peak_fx) > 0)
                {
                    Ep_peak_fx = Ep_tmp_fx[i];
                    move32();/*Q15 */
                }
            }
        }

        Mpy_32_16_ss(Ep_peak_fx,24576,&L_tmp,&lo);
        Mpy_32_16_ss(Ep_peak_fx,19661,&L_tmp2,&lo);
        Mpy_32_16_ss(Ep_avrgL_fx,24576,&L_tmp3,&lo);

        IF( (L_sub(L_shr(Ep_avrgL_fx,1), Ep_avrg_fx)>0 && L_sub(L_tmp,L_shr(Ep_avrgL_fx,2)) > 0  && L_sub(L_shr(Ep_avrgL_fx,1),L_tmp2) < 0 ) ||
            (L_sub(L_shr(Ep_avrg_fx,1), Ep_avrgL_fx)>0 && L_sub(L_shr(Ep_avrg_fx,3),L_tmp3) < 0 && L_sub(L_tmp,L_shr(Ep_avrgL_fx,2)) > 0 ) )
        {
            adjustFlag = 1;
            move16();
            FOR (i = 0; i < lowband; i++)
            {
                tmp = extract_h(L_shl(Ep_avrgL_fx,1));/*Q0 */
                IF(tmp != 0)
                {
                    exp = norm_s(tmp);
                    tmp = shl(tmp,exp);/*Q(exp) */
                    tmp = div_s(16384,tmp);/*Q(15+14-exp=29-exp) */
                    exp = sub(29,exp);
                }
                ELSE
                {
                    tmp = 0x7fff;
                    move16();
                    exp = 0;
                    move16();
                }
                Mpy_32_16_ss(Ep_peak_fx,tmp,&L_tmp,&lo);
                Mpy_32_16_ss(L_tmp,lowband,&L_tmp,&lo);
                L_tmp = L_shl(L_tmp,sub(28,exp));/*Q14 0.5 */
                tmp = extract_l(L_min(L_tmp,19661));/* Q14 */
                Mpy_32_16_ss(L_band_energy_tmp[i],tmp,&L_tmp,&lo);
                L_band_energy_tmp[i] = L_shl(L_tmp,1); /*Q_band_energy  */
            }
        }
        for (i = 0; i < bands; i++)
        {
            band_energy_tmp[i] = (float)(L_band_energy_tmp[i]/pow(2.0f, SWB_BWE_LR_Qbe));
        }

        hq2_bit_alloc( band_energy_tmp, bands, Rk_fx, &bit_budget, p2a_flags, bit_alloc_weight_fx,
                       band_width, num_bits, hqswb_clas, st->bwidth,*is_transient );
    }
    else if( st->bwidth == SWB && hqswb_clas == HQ_HARMONIC && (bwe_br == HQ_16k40 || bwe_br == HQ_13k20) )
    {
        hq2_bit_alloc_har( band_energy, bit_budget, bands, Rk_fx, p2a_bands,bwe_br,p2a_flags,band_width);
    }
    else
    {
        hq2_bit_alloc( band_energy, bands, Rk_fx, &bit_budget, p2a_flags, bit_alloc_weight_fx,
                       band_width, num_bits, hqswb_clas, st->bwidth, *is_transient );
    }

    tcq_core_LR_dec( st, inp_vector, bit_budget, bands, band_start, band_width, Rk_fx, npulses, k_sort,
                     p2a_flags, p2a_bands, last_bitalloc_max_band, inner_frame, adjustFlag, is_transient );

    /* Prepare floating Rk for next modules */
    for( k = 0; k < bands; k++)
    {
        Rk[k] = WORD322FL_SCALE( Rk_fx[k], SWB_BWE_LR_QRk - 1);
    }

    /* Denormalize the coded MDCT spectrum */
    mdct_spectrum_denorm( inp_vector, y2, band_start, band_end, band_width, band_energy, npulses, bands, ld_slope, pd_thresh );

    /* Apply fine gain to denormalized coded spectrum */
    mdct_spectrum_fine_gain_dec( st, y2, band_start, band_end, k_sort, bands, L_qint, Ngq, gqlevs, gqbits );

    /*restore the band information */
    if( flag_spt == 1 )
    {
        spt_shorten_domain_band_restore( bands, band_start, band_end, band_width, org_band_start, org_band_end, org_band_width );
    }

    mvr2r( y2, y2_org, L_FRAME32k );

    /* Inject noise into components having relatively low pulse energy per band */
    ni_seed = npulses[0] + npulses[1] + npulses[2] + npulses[3];

    for(i=0; i<band_end[bands-1]+1; i++) y2_ni[i] = y2[i];
    hq2_noise_inject( y2_ni, band_start, band_end, band_width, Ep, Rk, npulses, ni_seed, bands, 0, bw_low, bw_high, enerL, enerH,
                      st->last_ni_gain, st->last_env, &st->last_max_pos_pulse, p2a_flags, p2a_bands, hqswb_clas, st->bwidth, bwe_br );

    if( st->bwidth == SWB && (bwe_br == HQ_16k40 || bwe_br == HQ_13k20) )
    {
        if( hqswb_clas == HQ_NORMAL || hqswb_clas == HQ_HARMONIC )
        {
            preset_hq2_swb( hqswb_clas, band_end, &har_bands, p2a_bands,length, bands, &lowlength, &highlength, m );

            /*Gap filling for the core coder*/
            swb_bwe_dec_lr( st, y2, m, bwe_br , bands, band_start, band_end, band_energy, p2a_flags, hqswb_clas, lowlength,
                            highlength, har_bands, &st->prev_frm_hfe2, &st->prev_stab_hfe2, band_width,y2_ni, &ni_seed );

            post_hq2_swb( m, lowlength, highlength, hqswb_clas, har_bands, bands, p2a_flags, band_start, band_end, y2, npulses );

            if( hqswb_clas == HQ_NORMAL )
            {
                spt_swb_peakpos_tmp_save( y2, bands, band_start, band_end, prev_SWB_peak_pos_tmp );
                for( k=0; k<SPT_SHORTEN_SBNUM; k++)
                {
                    if( p2a_flags[bands-SPT_SHORTEN_SBNUM+k] == 0 || npulses[bands-SPT_SHORTEN_SBNUM+k] == 0 )
                    {
                        prev_SWB_peak_pos_tmp[k]=0;
                    }
                }
            }

            mvr2r( y2_ni, y2, lowlength );
        }
        else
        {
            mvr2r( y2_ni, y2, band_end[bands-1]+1 ); /* HQ_TRANSIENT */
        }
    }
    else
    {
        mvr2r( y2_ni, y2, band_end[bands-1]+1 );/* NB, WB */
    }


    /* bandwidth switching */
    if( !(st->last_inner_frame >= L_FRAME16k && st->bws_cnt > 0) )
    {
        k1 = *is_transient ? bands - 2 : bands - 6;
        st->prev_ener_shb = 0.0f;

        for( i = k1; i < bands; i++ )
        {
            st->prev_ener_shb += Ep_tmp[i]/(bands-k1);
        }
    }

    if( st->last_inner_frame >= L_FRAME32k )
    {
        set_f(st->prev_SWB_fenv, st->prev_ener_shb, SWB_FENV);
    }

    updat_prev_frm( y2, yout, bwe_br, length, inner_frame, bands, st->bwidth, *is_transient, hqswb_clas, &st->prev_hqswb_clas,
                    st->prev_SWB_peak_pos, prev_SWB_peak_pos_tmp, &st->prev_frm_hfe2, &st->prev_stab_hfe2, st->bws_cnt );

    return;
}


/*------------------------------------------------------------------------------------
 * small_symbol_dec_tran()
 *
 * Huffman decoding of differential energies
 *--------------------------------------------------------------------------------------*/
static short small_symbol_dec_tran(
    Decoder_State *st,            /* i/o: decoder state structure   */
    int   *qbidx,                 /* o  : output of dequantized differential energy */
    const short bands,            /* i  : number of bands */
    const short is_transient      /* i  : transient flag  */
)
{
    short i,  bits;
    short difidx[BANDS_MAX];

    /* Decoding differential energies*/
    bits = decode_envelope_indices(st, 0, bands, 0, difidx, LOW_RATE_HQ_CORE_TRAN ,is_transient);
    bits += BITS_DE_FCOMP;

    /* counting 1 bit for band_energy_huff_coding_mode */
    bits += BITS_DE_HMODE;

    /* converting to original values */
    for( i=0; i<bands; i++ )
    {
        qbidx[i] = difidx[i]-LRMDCT_BE_OFFSET;
    }

    return( bits );
}


/*--------------------------------------------------------------------------
 * small_symbol_dec()
 *
 * Huffman decoding of differential energies (MSB and LSB)
 *--------------------------------------------------------------------------*/
static short small_symbol_dec( /* o  : bits                                        */
    Decoder_State *st,            /* i/o: decoder state structure   */
    int   *qbidx,         /* o  : output of dequantized differential energy   */
    const short bands           /* i  : number of bands                             */
    ,const short is_transient
)
{
    short i, LSB, bits;
    short difidx[BANDS_MAX];

    /* Decoding MSB bits */
    bits = decode_envelope_indices( st, 0, bands, 0, difidx, LOW_RATE_HQ_CORE, is_transient );
    bits += BITS_DE_FCOMP;

    /* counting 1 bit for band_energy_huff_coding_mode */
    bits += BITS_DE_HMODE;

    /* Decoding LSB bit packing */
    for( i=0; i<bands; ++i )
    {
        LSB = (short)get_next_indice( st, BITS_DE_LSB );
        difidx[i] = (difidx[i]<<1)|LSB;
    }

    /* counting bits for LSB */
    bits += bands;

    /* converting to original values */
    for( i=1; i<bands; ++i )
    {
        qbidx[i] = difidx[i]-DE_OFFSET1;
    }

    qbidx[0] = difidx[0]-DE_OFFSET0;

    return( bits );
}


static
short decode_huff_8s(
    Decoder_State *st,
    const short *hufftab,
    short *rbits
)
{
    short bit;

    while( *hufftab > 0)
    {
        *rbits+=(*hufftab & 0xf);
        bit = (short)get_next_indice( st, *hufftab & 0xf );
        hufftab += (*hufftab >> 4) + bit;
    }

    return (-*hufftab);
}

static short large_symbol_dec( /* o  : bits                                        */
    Decoder_State *st,                    /* i/o: decoder state structure                     */
    int   *qbidx,                 /* o  : output of dequantized differential energy   */
    const short bands                   /* i  : number of bands                             */
)
{
    short i, bits;
    short LSB[BANDS_MAX];
    short basic_shift,cntbits,ns2mode;
    short pos_outlyer;
    short ns2mode0,ns2mode1;

    cntbits = BITS_DE_8SMODE;
    ns2mode = (short)get_next_indice (st, BITS_DE_8SMODE);

    if (ns2mode == 0 )
    {
        ns2mode0 = (short)get_next_indice (st, BITS_DE_8SMODE_N0);
        ns2mode1 = (short)get_next_indice (st, BITS_DE_8SMODE_N1);
        cntbits += BITS_DE_8SMODE_N0+BITS_DE_8SMODE_N1;

        if (ns2mode0 == 0)
        {
            if (ns2mode1 == 1)
            {
                pos_outlyer = (short)get_next_indice (st, BITS_DE_8SPOS);
                cntbits+=BITS_DE_8SPOS;
                qbidx[pos_outlyer] = ((short)get_next_indice (st, BITS_ABS_ENG) - ABS_ENG_OFFSET);
                cntbits+=BITS_ABS_ENG;
            }
            else
            {
                pos_outlyer = -1;
            }

            for( i=0; i<pos_outlyer; ++i )
            {
                bits = 0;
                qbidx[i] = (decode_huff_8s(st, hestable, &bits) - 4);
                cntbits += bits;
            }

            for( i=pos_outlyer+1; i<bands; ++i )
            {
                bits = 0;
                qbidx[i] = (decode_huff_8s(st, hestable, &bits) - 4);
                cntbits += bits;
            }
        }
        else
        {
            if (ns2mode1 == 1)
            {
                pos_outlyer = (short)get_next_indice (st, BITS_DE_8SPOS);
                cntbits+=BITS_DE_8SPOS;
                qbidx[0] = ((short)get_next_indice (st, BITS_ABS_ENG) - ABS_ENG_OFFSET);
                cntbits+=BITS_ABS_ENG;
                qbidx[pos_outlyer] = ((short)get_next_indice (st, BITS_ABS_ENG) - ABS_ENG_OFFSET);
                cntbits+=BITS_ABS_ENG;
            }
            else
            {
                pos_outlyer = 0;
                qbidx[0] = ((short)get_next_indice (st, BITS_ABS_ENG) - ABS_ENG_OFFSET);
                cntbits+=BITS_ABS_ENG;
            }

            for( i=1; i<pos_outlyer; ++i )
            {
                bits = 0;
                qbidx[i] = (decode_huff_8s(st, hestable, &bits) - 4);
                cntbits += bits;
            }

            for( i=pos_outlyer+1; i<bands; ++i )
            {
                bits = 0;
                qbidx[i] = (decode_huff_8s(st, hestable, &bits) - 4);
                cntbits += bits;
            }
        }
    }
    else
    {
        basic_shift = (short)get_next_indice (st, BITS_MAX_DEPTH);
        cntbits += BITS_MAX_DEPTH;

        for(i=0; i<bands; ++i)
        {
            bits = 0;
            qbidx[i] = (decode_huff_8s(st, hestable, &bits) - 4);
            qbidx[i] = qbidx[i]<<basic_shift;
            cntbits += bits;
        }

        for(i=0; i<bands; ++i)
        {
            LSB[0] = (short)get_next_indice (st, basic_shift);
            qbidx[i] += LSB[0];
            cntbits += basic_shift;
        }
    }

    return cntbits;  /* xx bits for diff. energies + 1 bit for LC coding mode */
}


/*--------------------------------------------------------------------------*
 * band_energy_dequant()
 *
 *
 *--------------------------------------------------------------------------*/

static float band_energy_dequant(
    Decoder_State *st,              /* i/o: decoder state structure   */
    float band_energy[],
    const short bands,
    const Word32 L_qint,  /* Q29 */
    const Word16 eref_fx, /* Q10 */
    const short is_transient
)
{
    short k;
    short deng_cmode;
    int   bq0;
    int   bq1[BANDS_MAX], bq2[BANDS_MAX];
    float deng_bits;
    Word32 L_band_energy[BANDS_MAX]; /* SWB_BWE_LR_Qbe */
    Word16 exp_normd;
    Word16 rev_qint_fx;
    Word16 Qrev_qint;

    /* parsing energy difference coding mode */
    deng_cmode = (short) get_next_indice( st, BITS_DE_CMODE );

    if ( deng_cmode == 0 )
    {
        deng_bits = (float)large_symbol_dec( st, bq2, bands );

        /* counting 1 bit for deng coding mode */
        deng_bits += BITS_DE_CMODE;
    }
    else
    {
        if(is_transient)
        {
            deng_bits = (float)small_symbol_dec_tran(st, bq2, bands , is_transient);
        }
        else
        {
            deng_bits = (float)small_symbol_dec( st, bq2, bands, is_transient );
        }

        /* counting 1 bit for deng coding mode */
        deng_bits += BITS_DE_CMODE;
    }

    /* Get the reference energy */
    exp_normd = norm_l(L_qint);
    rev_qint_fx = div_s(0x4000, round_fx(L_shl(L_qint, exp_normd))); /* Q14-(29+exp_normd-16)+15 */
    Qrev_qint = sub(14-(29-16)+15, exp_normd);

    bq0 = (int)round_fx(L_shl(L_mult(eref_fx, rev_qint_fx), sub(5, Qrev_qint))); /* 16-(10+Qrev_qint+1) */

    /* Reconstruct quantized spectrum */
    bq1[0] = bq2[0] + bq0;
    for (k = 1; k < bands; k++)
    {
        bq1[k] = bq2[k] + bq1[k - 1];
    }

    for (k = 0; k < bands; k++)
    {
        L_band_energy[k] = L_mls(L_qint, (short)bq1[k]);
        move32();/* 29+0-15 -> Qbe(Q14) */
        band_energy[k] = (float)(L_band_energy[k]/pow(2.0f, SWB_BWE_LR_Qbe));
    }

    if( is_transient )
    {
        reverse_transient_frame_energies( band_energy, bands );
    }

    return( deng_bits );
}


/*--------------------------------------------------------------------------*
 * p2a_threshold_dequant()
 *
 *
 *--------------------------------------------------------------------------*/

static short p2a_threshold_dequant(
    Decoder_State *st,              /* i/o: decoder state structure   */
    short *p2a_flags,
    const short bands,
    const short p2a_bands
)
{
    short j, k;

    for( k = 0; k < bands - p2a_bands; k++ )
    {
        p2a_flags[k] = 1;
    }

    j = 0;
    for( k = bands - p2a_bands; k < bands; k++ )
    {
        p2a_flags[k] = (short) get_next_indice( st, 1 );
        j++;
    }

    return( j );
}


/*--------------------------------------------------------------------------*
 * mdct_spectrum_fine_gain_dec()
 *
 *
 *--------------------------------------------------------------------------*/

static void mdct_spectrum_fine_gain_dec(
    Decoder_State *st,              /* i/o: decoder state structure   */
    float y2[],
    const short band_start[],
    const short band_end[],
    const short k_sort[],
    const short bands,
    const Word32 L_qint,
    const short Ngq,
    const short gqlevs,
    const short gqbits
)
{
    short i, k, imin;
    float gamma;
    float gain_table[MAX_GQLEVS];
    Word16 exp_normn, exp_normd;
    Word16 delta_fx, Qdelta;
    Word32 L_delta, L_q;
    Word32 L_temp;
    Word16 gain_table_fx[MAX_GQLEVS];
    Word16 Qgt;
    Word16 temp_lo_fx, temp_hi_fx;

    /* Fine gain quantization on only the most significant energy bands */
    /*delta = qint / gqlevs; */
    exp_normn = norm_l(L_qint);
    exp_normn = sub(exp_normn, 1);
    exp_normd = norm_s(gqlevs);
    delta_fx = div_l(L_shl(L_qint, exp_normn), shl(gqlevs, exp_normd));
    Qdelta = add(sub(exp_normn, exp_normd), 28); /* 29+exp_normn-(exp_normd)-1; */
    L_delta = L_shl(L_deposit_h(delta_fx), sub(13, Qdelta));
    /*q = (-qint + delta) / 2.0f; */
    L_q = L_shr(L_sub(L_delta, L_qint), 1);

    FOR (i = 0; i < gqlevs; i++)
    {
        /*gain_table[i] = (float) pow (2.0f, q * 0.5f); */
        L_temp = L_shr(L_shr(L_q, 1), sub(29, 16));
        temp_lo_fx = L_Extract_lc(L_temp, &temp_hi_fx);
        Qgt = sub(14, temp_hi_fx);
        gain_table_fx[i] = extract_l(Pow2(14, temp_lo_fx));        /* Qgt */

        /*q += delta; */
        L_q = L_add(L_q, L_delta);
        gain_table_fx[i] = shl(gain_table_fx[i], sub(14, Qgt)); /* Qgt -> Q14 */
        gain_table[i] = (float)(gain_table_fx[i]/pow(2.0f, 14));
    }

    for( k = bands - Ngq; k < bands; k++ )
    {
        imin = (short)get_next_indice( st, gqbits );
        gamma = gain_table[imin];

        for( i = band_start[k_sort[k]]; i <= band_end[k_sort[k]]; i++ )
        {
            y2[i] *= gamma;
        }
    }

    return;
}

/*--------------------------------------------------------------------------*
 * spt_shorten_domain_set_dec()
 *
 * update the shorten band information based on p2a analysis
 *--------------------------------------------------------------------------*/

static void spt_shorten_domain_set_dec(
    Decoder_State *st,                 /* i:   encoder state structure             */
    const short p2a_flags[],         /* i:   p2a anlysis information             */
    const short new_band_start[],    /* i:   new band start position             */
    const short new_band_end[],      /* i:   new band end position               */
    const short new_band_width[],    /* i:   new subband band width              */
    const short bands,               /* i:   total number of subbands            */
    short band_start[],        /* o:   band start position                 */
    short band_end[],          /* o:   band end position                   */
    short band_width[],        /* o:   sub band band width                 */
    short *bit_budget          /* i/o: bit budget                          */
)
{
    int   j,k;
    short kpos;
    short spt_shorten_flag[SPT_SHORTEN_SBNUM];

    kpos = 0;
    j = 0;
    for( k = bands - SPT_SHORTEN_SBNUM; k < bands; k++ )
    {
        spt_shorten_flag[j] = 0;
        if(p2a_flags[k] == 1)
        {
            spt_shorten_flag[j] = (short)get_next_indice (st, 1 );
            *bit_budget -= 1;
            if(spt_shorten_flag[j] == 1)
            {
                band_start[k] = new_band_start[j];
                band_end[k]   = new_band_end[j];
                band_width[k] = new_band_width[j];
            }
        }

        kpos++;
        j++;
    }

    return;
}
