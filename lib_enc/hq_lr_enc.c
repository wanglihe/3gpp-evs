/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "rom_enc.h"
#include "rom_com.h"
#include "prot.h"
#include "stl.h"
#include "basop_util.h"

/*--------------------------------------------------------------------------*
 * Local functions
 *--------------------------------------------------------------------------*/

static short small_symbol_enc( Encoder_State *st, const int *qbidx, const short bands, short *hLCmode, const short flag_pack, const short is_transient );

static short small_symbol_enc_tran( Encoder_State *st, const int *qbidx, const short bands, short *hLCmode, const short flag_pack ,const short is_transient );

static float band_energy_quant( Encoder_State *st, const float *t_audio, const short band_start[], const short band_end[], float band_energy[],
                                const short bands, const Word32 L_qint, const Word16 eref_fx, const short is_transient );

static short p2a_threshold_quant( Encoder_State *st, const float *t_audio, const short band_start[], const short band_end[], const short band_width[],
                                  const short bands, const short p2a_bands, const float p2a_th, short *p2a_flags );

static void mdct_spectrum_fine_gain_enc( Encoder_State *st, const float ybuf[], float y2[], const short band_start[], const short band_end[],
        const short k_sort[], const short bands, const Word32 L_qint, const short Ngq, const short gqlevs, const short gqbits );

/*--------------------------------------------------------------------------*
 * spt_shorten_domain_set()
 *
 * Track the spectral peak based on peak -avg analysis
 *--------------------------------------------------------------------------*/

static void spt_shorten_domain_set(
    Encoder_State *st,                 /* i:   encoder state structure             */
    const float t_audio[],           /* i:   input spectrum                      */
    const short p2a_flags[],         /* i:   p2a anlysis information             */
    const short new_band_start[],    /* i:   new band start position             */
    const short new_band_end[],      /* i:   new band end position               */
    const short new_band_width[],    /* i:   new subband band width              */
    const short bands,               /* i:   total number of subbands            */
    short       band_start[],        /* i/o: band start position                 */
    short       band_end[],          /* i/o: band end position                   */
    short       band_width[],        /* i:   sub band band width                 */
    short       *bit_budget          /* i/o: bit budget                          */
)
{
    short   i, j, k;
    short kpos;
    float max_y2;
    short max_y2_pos;
    short spt_shorten_flag[SPT_SHORTEN_SBNUM];

    kpos = 0;
    j = 0;
    for( k=bands-SPT_SHORTEN_SBNUM; k<bands; k++ )
    {
        if( p2a_flags[k] == 1 )
        {
            spt_shorten_flag[j] = 0;
            if ( st->prev_SWB_peak_pos[kpos] != 0)
            {
                max_y2 = 0.0f;
                max_y2_pos = 0;
                for( i=band_start[k]; i<=band_end[k]; i++ )
                {
                    if( max_y2 < fabs(t_audio[i]) )
                    {
                        max_y2 =(float) fabs(t_audio[i]);
                        max_y2_pos = i;
                    }
                }
                if( max_y2_pos >= new_band_start[j] && max_y2_pos <= new_band_end[j] )
                {
                    band_start[k] = new_band_start[j];
                    band_end[k]   = new_band_end[j];
                    band_width[k] = new_band_width[j];
                    spt_shorten_flag[j] = 1;
                }
            }
            push_indice( st, IND_HQ2_SPT_SHORTEN, spt_shorten_flag[j], 1 );
            *bit_budget -= 1;
        }

        kpos++;
        j++;
    }

    return;
}

/*--------------------------------------------------------------------------*
 * hq_lr_enc()
 *
 * HQ low rate encoding routine
 *--------------------------------------------------------------------------*/

void hq_lr_enc(
    Encoder_State *st,             /* i/o: encoder state structure             */
    float t_audio[],      /* i/o: transform-domain coefs.             */
    const short inner_frame,    /* i  : inner frame length                  */
    short *num_bits,      /* i/o: number of available bits            */
    const short is_transient    /* i  : transient flag                      */
)
{
    short i, k1, k2;
    short bit_budget, pbits;
    short bands, length, ni_seed, gqlevs, gqbits, Ngq, p2a_bands;
    short p2a_flags[BANDS_MAX];
    short band_start[BANDS_MAX], band_end[BANDS_MAX], band_width[BANDS_MAX];
    float band_energy[BANDS_MAX], Rk[BANDS_MAX];
    Word32 Rk_fx[BANDS_MAX];
    float ebits;
    float p2a_th, ni_coef, ni_pd_th, pd_thresh, ld_slope;
    Word32 L_qint;              /* Q29 */
    Word16 eref_fx;             /* Q10 */
    Word16 bit_alloc_weight_fx; /* Q13 */
    short k_sort[BANDS_MAX];
    int npulses[BANDS_MAX];
    int inp_vector[L_FRAME48k];
    float y2[L_FRAME48k];
    float y2_ni[L_FRAME48k];
    short hqswb_clas;
    short lowlength;
    short highlength;
    float m[L_FRAME32k];
    short har_bands;
    float Ep[BANDS_MAX], enerH = 0.0f, enerL = 0.0f;
    short lowband, highband, bw_low = 0, bw_high = 20;
    float band_energy_tmp[BANDS_MAX];
    long bwe_br;
    short trans_bit, p2a_flags_tmp[BANDS_MAX];
    short adjustFlag = 0;
    short prev_SWB_peak_pos_tmp[SPT_SHORTEN_SBNUM];
    int  k,j;
    short flag_spt;
    short org_band_start[SPT_SHORTEN_SBNUM];
    short org_band_end[SPT_SHORTEN_SBNUM];
    short org_band_width[SPT_SHORTEN_SBNUM];
    short new_band_start[SPT_SHORTEN_SBNUM];
    short new_band_end[SPT_SHORTEN_SBNUM];
    short new_band_width[SPT_SHORTEN_SBNUM];
    short bws_cnt=0;
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
    set_f( y2, 0.0f, L_FRAME48k );
    set_i( inp_vector, 0, inner_frame );
    flag_spt = 0;
    set_s(prev_SWB_peak_pos_tmp, 0, SPT_SHORTEN_SBNUM);

    bwe_br = st->core_brate;
    hqswb_clas = HQ_NORMAL;
    if( st->bwidth == SWB && (bwe_br == HQ_16k40 || bwe_br == HQ_13k20) )
    {
        if ( is_transient == 1 )
        {
            hqswb_clas = HQ_TRANSIENT;
        }
        else
        {
            /* classification of HQ_HARMONIC and HQ_NORMAL frames for SWB BWE */
            hqswb_clas = peak_avrg_ratio( st->total_brate, t_audio, NUMC_N, &st->mode_count, &st->mode_count1);
        }

        /* write the classification information into the bitstream */
        push_indice( st, IND_HQ2_SWB_CLAS, hqswb_clas, 2 );
        (*num_bits) -= 2;

        if( hqswb_clas == HQ_NORMAL )
        {
            flag_spt = 1;
        }
    }
    else
    {
        /* write the transient bit into the bitstream */
        push_indice( st, IND_HQ2_SWB_CLAS, is_transient, 1 );

        /* subtract one bit for the transient flag */
        (*num_bits)--;
    }

    /* Configure encoder for different bandwidths, bit rates, etc. */
    hq2_core_configure( inner_frame, *num_bits, is_transient, &bands, &length, band_width, band_start, band_end, &L_qint, &eref_fx,
                        &bit_alloc_weight_fx, &gqlevs, &Ngq, &p2a_bands, &p2a_th, &pd_thresh, &ld_slope, &ni_coef, &ni_pd_th, bwe_br );

    highlength = band_end[bands-1];
    har_bands = bands;

    if( st->bwidth == SWB && is_transient == 0 && (bwe_br == HQ_16k40 || bwe_br == HQ_13k20) )
    {
        /* reserve bits for HQ_NORMAL and HQ_HARMONIC modes */
        if( hqswb_clas == HQ_NORMAL || hqswb_clas == HQ_HARMONIC )
        {
            (*num_bits) -= (short)get_usebit_npswb(hqswb_clas);
        }
    }

    if( (bwe_br == HQ_16k40 || bwe_br == HQ_13k20) && st->bwidth == SWB )
    {
        if( st->prev_hqswb_clas != HQ_NORMAL )
        {
            j = 0;
            for( k=bands-SPT_SHORTEN_SBNUM; k<bands; k++ )
            {
                st->prev_SWB_peak_pos[j] = 0;
                j++;
            }
        }
    }

    /* Check if input frame is larger than coded bandwidth */
    if( inner_frame > length && is_transient )
    {
        /* If so, collapse transient frame (4 short transforms) to remove uncoded coefficients */
        for( i = 1; i < NUM_TIME_SWITCHING_BLOCKS; i++ )
        {
            k1 = i*length/NUM_TIME_SWITCHING_BLOCKS;
            k2 = i*inner_frame/NUM_TIME_SWITCHING_BLOCKS;

            mvr2r( &t_audio[k2], &t_audio[k1], length/NUM_TIME_SWITCHING_BLOCKS );
        }
    }

    /* Spectral energy calculation/quantization */
    ebits = band_energy_quant( st, t_audio, band_start, band_end, band_energy, bands, L_qint, eref_fx, is_transient );

    for( i = 0; i < bands; i++)
    {
        L_band_energy[i] = (Word32)(band_energy[i] * pow(2.0f, SWB_BWE_LR_Qbe));
    }

    /* First pass bit budget for TCQ of spectral band information */
    gqbits = (short int) log2_f ((float) gqlevs);
    bit_budget = (*num_bits) - (short) ceil (ebits) - Ngq * gqbits;

    pbits = 0;
    if( st->bwidth == SWB && (bwe_br == HQ_16k40 || bwe_br == HQ_13k20) )
    {
        if( hqswb_clas == HQ_HARMONIC )
        {
            set_s( p2a_flags, 1, har_bands );
        }
        else
        {
            /* High band tonality detector based on per band peak-to-average ratio */
            pbits = p2a_threshold_quant( st, t_audio, band_start, band_end, band_width, bands, p2a_bands, p2a_th, p2a_flags );
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
        pbits = p2a_threshold_quant( st, t_audio, band_start, band_end, band_width, bands, p2a_bands, p2a_th, p2a_flags );
        bit_budget -= pbits;
    }

    if( flag_spt == 1 )
    {
        /* initalize the desired parameters for SPT */
        spt_shorten_domain_band_save( bands, band_start, band_end, band_width, org_band_start, org_band_end, org_band_width );
        spt_shorten_domain_pre( band_start, band_end, st->prev_SWB_peak_pos, bands, bwe_br, new_band_start, new_band_end, new_band_width );
        spt_shorten_domain_set( st, t_audio, p2a_flags, new_band_start, new_band_end, new_band_width, bands, band_start, band_end, band_width, &bit_budget );
    }

    /* Estimate number of bits per band */

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
    }


    if ( is_transient == 0 && inner_frame == L_FRAME8k && st->core_brate <= ACELP_13k20 )
    {
        lowband = 6;
        move16();
        trans_bit = 2;
        move16();
        bit_budget =sub(bit_budget,trans_bit);
        gama_fx = 27852; /*Q15 0.85f; */
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
                        /*when the divisor is zero, happens rarely*/
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
            j = 0;
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
                        /*when the divisor is zero, happens rarely*/
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
                IF(sub(st->last_bitalloc_max_band[j++], 1) == 0)
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
                    /*when the divisor is zero, happens rarely*/
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

        hq2_bit_alloc( band_energy_tmp, bands, Rk_fx, &bit_budget, p2a_flags, bit_alloc_weight_fx, band_width,
                       *num_bits, hqswb_clas, st->bwidth, is_transient );

        /* encode the last p2a_bands-1 subbands bit-allocation index of the previous frame */
        for(i = 0; i < 2; i++)
        {
            push_indice ( st, IND_HQ2_LAST_BA_MAX_BAND, st->last_bitalloc_max_band[i], 1 );
        }
    }
    else if( is_transient == 0 && inner_frame == L_FRAME16k )
    {
        bit_budget = sub(bit_budget,2);/* bits in high bands to indicate the last 2 subbands is allocated bits or not */
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
                IF(sub(st->last_bitalloc_max_band[sub(i, add(tmp, 1))], 1) == 0)
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
                    alpha_fx = extract_l(L_shr(L_tmp,15)); /*Q14*/
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
                tmp = extract_l(L_min(L_tmp,19661));/*Q14 */
                Mpy_32_16_ss(L_band_energy_tmp[i],tmp,&L_tmp,&lo);
                L_band_energy_tmp[i] = L_shl(L_tmp,1); /*Q_band_energy  */
            }
        }
        for (i = 0; i < bands; i++)
        {
            band_energy_tmp[i] = (float)(L_band_energy_tmp[i]/pow(2.0f, SWB_BWE_LR_Qbe));
        }

        hq2_bit_alloc( band_energy_tmp, bands, Rk_fx, &bit_budget, p2a_flags, bit_alloc_weight_fx, band_width,
                       *num_bits, hqswb_clas, st->bwidth, is_transient );

        /* encode the last p2a_bands-1 subbands bit-allocation index of the previous frame */
        for(i = 0; i < 2; i++)
        {
            push_indice( st, IND_HQ2_LAST_BA_MAX_BAND, st->last_bitalloc_max_band[i], 1 );
        }
    }
    else if( st->bwidth == SWB && hqswb_clas == HQ_HARMONIC && (bwe_br == HQ_16k40 || bwe_br == HQ_13k20) )
    {
        /* bit allocation for harmonic mode */
        hq2_bit_alloc_har( band_energy, bit_budget, bands, Rk_fx, p2a_bands,bwe_br, p2a_flags, band_width);
    }
    else
    {

        /* estimate number of bits per band */
        hq2_bit_alloc( band_energy, bands, Rk_fx, &bit_budget, p2a_flags, bit_alloc_weight_fx, band_width,
                       *num_bits, hqswb_clas, st->bwidth, is_transient );
    }

    tcq_core_LR_enc( st, inp_vector, t_audio, y2, bit_budget, bands, band_start, band_end, band_width, Rk_fx,
                     npulses, k_sort,  p2a_flags, p2a_bands, st->last_bitalloc_max_band, inner_frame, adjustFlag, is_transient );

    if((inner_frame == L_FRAME8k && st->core_brate <= ACELP_13k20) || inner_frame == L_FRAME16k)
    {
        j = 0;
        for(i = 2; i > 0; i--)
        {
            if(npulses[bands-i] > 0)
            {
                st->last_bitalloc_max_band[j] = 1;
            }
            else
            {
                st->last_bitalloc_max_band[j] = 0;
            }
            j++;
        }
    }

    /* Prepare floating Rk for next modules */
    for( k = 0; k < bands; k++)
    {
        Rk[k] = WORD322FL_SCALE( Rk_fx[k], SWB_BWE_LR_QRk - 1);
    }

    /* Denormalize the coded MDCT spectrum */
    mdct_spectrum_denorm( inp_vector, y2, band_start, band_end, band_width, band_energy, npulses, bands, ld_slope, pd_thresh );

    /* Apply fine gain quantization to denormalized coded spectrum */
    mdct_spectrum_fine_gain_enc( st, t_audio, y2, band_start, band_end, k_sort, bands, L_qint, Ngq, gqlevs, gqbits );

    /* reStore the subband information*/
    if(flag_spt == 1)
    {
        spt_shorten_domain_band_restore(bands, band_start, band_end, band_width, org_band_start, org_band_end, org_band_width);
    }

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

            swb_bwe_enc_lr( st, y2, t_audio, m , bwe_br, bands, band_start, band_end, band_energy, p2a_flags,
                            hqswb_clas, lowlength, highlength, st->prev_frm_index, har_bands, &st->prev_frm_hfe2, &st->prev_stab_hfe2,band_width, y2_ni, &ni_seed );

            post_hq2_swb( m, lowlength, highlength, hqswb_clas, har_bands, bands, p2a_flags, band_start, band_end, y2, npulses );

            if( hqswb_clas == HQ_NORMAL )
            {
                spt_swb_peakpos_tmp_save(y2, bands, band_start, band_end, prev_SWB_peak_pos_tmp);
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
            mvr2r( y2_ni, y2, band_end[bands-1]+1 );            /* HQ_TRANSIENT */
        }
    }
    else
    {
        mvr2r(y2_ni, y2, band_end[bands-1]+1);        /* NB, WB */
    }


    updat_prev_frm( y2, t_audio, bwe_br, length, inner_frame, bands, st->bwidth, is_transient, hqswb_clas, &st->prev_hqswb_clas,
                    st->prev_SWB_peak_pos, prev_SWB_peak_pos_tmp, &st->prev_frm_hfe2, &st->prev_stab_hfe2, bws_cnt );

    if( st->bwidth != SWB )
    {
        /* reset HQ classifier memories */
        st->mode_count = 0;
        st->mode_count1 = 0;
    }
    if( hqswb_clas != HQ_HARMONIC && (bwe_br == HQ_16k40 || bwe_br == HQ_13k20) && st->bwidth == SWB )
    {
        st->prev_frm_index[0] = -1;
        st->prev_frm_index[1] = -1;
    }
    /* update number of unused bits */
    *num_bits = 0;

    st->hvq_hangover = 0;

    return;
}

/*--------------------------------------------------------------------------*
 * small_symbol_enc_tran()
 *
 * Huffman encoding of differential energies, estimating or packing bits
 * if flag_pack = 0, LC mode info. is output else LC mode info. is input
 * if flag_pack = 0, estimatng else packing bits
 *--------------------------------------------------------------------------*/
static short small_symbol_enc_tran(  /* o  : bits                                     */
    Encoder_State *st,               /* i/o: encoder state structure                  */
    const int   *qbidx,              /* i  : input of dequantized differential energy */
    const short BANDS,               /* i  : number of bands                          */
    short *hLCmode,            /* i/o: LC mode info                             */
    const short flag_pack,           /* i  : indicator of packing or estimating bits  */
    const short is_transient
)
{
    short i, bits;
    short difidx[BANDS_MAX];

    for( i=0; i<BANDS; i++ )
    {
        difidx[i] = (short)(qbidx[i]+LRMDCT_BE_OFFSET);
    }

    for( i=0; i<BANDS; ++i )
    {
        if (difidx[i]>LRMDCT_BE_LIMIT || difidx[i]<0)
        {
            /* Huffman cannot encode this vector */
            return -1;
        }
    }

    /* Preparing lossless coding input */
    if ( flag_pack == 0 )
    {
        /* estimating # of bits */
        bits = encode_envelope_indices(st, BANDS, -1, difidx, hLCmode, flag_pack, LOW_RATE_HQ_CORE_TRAN , is_transient );
        bits += BITS_DE_FCOMP;   /* xx bits for diff. energies + BITS_DE_FCOMP bits for first energies */
    }
    else
    {
        bits = 0;
        encode_envelope_indices(st, BANDS, -1, difidx, hLCmode, flag_pack, LOW_RATE_HQ_CORE_TRAN ,is_transient );
    }

    return bits + BITS_DE_HMODE;  /* xx bits for diff. energies + 1 bit for LC coding mode */
}


/*--------------------------------------------------------------------------*
 * small_symbol_enc()
 *
 * Huffman encoding of differential energies, estimating or packing bits
 * if flag_pack = 0, LC mode info. is output else LC mode info. is input
 * if flag_pack = 0, estimatng else packing bits
 *--------------------------------------------------------------------------*/

static short small_symbol_enc(  /* o  : bits                                     */
    Encoder_State *st,            /* i/o: encoder state structure                  */
    const int   *qbidx,         /* i  : input of dequantized differential energy */
    const short BANDS,          /* i  : number of bands                          */
    short *hLCmode,       /* i/o: LC mode info                             */
    const short flag_pack,      /* i  : indicator of packing or estimating bits  */
    const short is_transient
)
{
    short i, bits;
    short difidx[BANDS_MAX], LSB[BANDS_MAX];

    /* Preparing lossless coding input */
    difidx[0] = (short)(qbidx[0]+DE_OFFSET0);

    for( i=1; i<BANDS; ++i )
    {
        difidx[i] = (short)(qbidx[i]+DE_OFFSET1);
    }

    for( i=0; i<BANDS; ++i )
    {
        if (difidx[i]>=DE_LIMIT || difidx[i]<0)
        {
            /* Huffman cannot encode this vector */
            return -1;
        }
    }

    /* splitting MSB and LSB */
    for( i=0; i<BANDS; ++i )
    {
        LSB[i] = difidx[i]&1;
        difidx[i] >>= 1;
    }

    /* Preparing lossless coding input */
    if ( flag_pack == 0 )
    {
        /* estimating # of bits */
        /* Encoding MSB bits */
        bits = encode_envelope_indices( st, BANDS, -1, difidx, hLCmode, flag_pack, LOW_RATE_HQ_CORE, is_transient );
        bits += BITS_DE_FCOMP;   /* xx bits for diff. energies + BITS_DE_FCOMP bits for first energies */

        /* Encoding LSB bit packing */
        bits += BANDS;
    }
    else
    {
        /* Encoding MSB bits */
        bits = 0;
        encode_envelope_indices( st, BANDS, -1, difidx, hLCmode, flag_pack, LOW_RATE_HQ_CORE, is_transient );

        /* Encoding LSB bit packing */
        for( i=0; i<BANDS; ++i )
        {
            push_indice( st,IND_HQ2_DIFF_ENERGY, LSB[i], BITS_DE_LSB);
        }
    }

    return bits + BITS_DE_HMODE;  /* xx bits for diff. energies + 1 bit for LC coding mode */
}


/*--------------------------------------------------------------------------*
 * large_symbol_enc()
 *
 *
 *--------------------------------------------------------------------------*/

static short large_symbol_enc(  /* o  : bits                                     */
    Encoder_State *st,          /* i  : encoder state structure                  */
    int   *qbidx,       /* i  : input of dequantized differential energy */
    const short BANDS,        /* i  : number of bands                          */
    short *hLCmode0,    /* i/o: LC mode info                             */
    const short flag_pack     /* i  : indicator of packing or estimating bits  */
)
{
    short i, bits;
    short LSB1[BANDS_MAX];
    short min_q=513,max_q=-1,offset0;
    short min_bits,min_bits_pos;
    short tdifidx0[BANDS_MAX], tdifidx1[BANDS_MAX];
    short basic_shift;
    short bitsmode0,bitsmode1;
    short lsbdepth1;
    short cnt_outlyer,pos_outlyer,cnt_outlyer0;

    cnt_outlyer0 = 0;
    cnt_outlyer = 0;
    bitsmode0 = 0;
    bitsmode1 = 0;
    pos_outlyer = 0;
    lsbdepth1 = 0;

    if ( flag_pack == 0 || ( flag_pack == 1 && *hLCmode0 == 0) )
    {
        if (qbidx[0]>ABS_ENG_OFFSET-1 || qbidx[0]<-ABS_ENG_OFFSET)
        {
            cnt_outlyer0 = 2;
        }
        else if (qbidx[0]>3 || qbidx[0]<-4)
        {
            cnt_outlyer0 = 1;
        }
        else
        {
            cnt_outlyer0 = 0;
        }

        cnt_outlyer=0;
        pos_outlyer = -1;
        for( i=1; i<BANDS; ++i )
        {
            if (qbidx[i]>3 || qbidx[i]<-4)
            {
                cnt_outlyer++;
                pos_outlyer = i;
            }

            if (qbidx[i]>ABS_ENG_OFFSET-1 || qbidx[i]<-ABS_ENG_OFFSET)
            {
                cnt_outlyer++;
            }
        }

        if (cnt_outlyer0 == 0 && cnt_outlyer<=1)
        {
            bitsmode0 = BITS_DE_8SMODE + BITS_DE_8SMODE_N0 + BITS_DE_8SMODE_N1;
            if (cnt_outlyer == 1)
            {
                /* 01 */
                bitsmode0 += BITS_DE_8SPOS + BITS_ABS_ENG;
            }

            for( i=0; i<pos_outlyer; ++i )
            {
                tdifidx0[i] = (short)(qbidx[i]);
                bitsmode0 += hessize[tdifidx0[i]+4];
            }

            for( i=pos_outlyer+1; i<BANDS; ++i )
            {
                tdifidx0[i] = (short)(qbidx[i]);
                bitsmode0 += hessize[tdifidx0[i]+4];
            }
        }
        else if (cnt_outlyer0 == 1 && cnt_outlyer<=1)
        {
            bitsmode0 = BITS_DE_8SMODE + BITS_DE_8SMODE_N0 + BITS_DE_8SMODE_N1;
            tdifidx0[0] = (short)(qbidx[0]);
            bitsmode0 += BITS_ABS_ENG;
            if (cnt_outlyer == 1)
            {
                /* 11 */
                bitsmode0 += BITS_DE_8SPOS + BITS_ABS_ENG;
            }
            else
            {
                pos_outlyer = 0;
            }

            for( i=1; i<pos_outlyer; ++i )
            {
                tdifidx0[i] = (short)(qbidx[i]);
                bitsmode0 += hessize[tdifidx0[i]+4];
            }

            for( i=pos_outlyer+1; i<BANDS; ++i )
            {
                tdifidx0[i] = (short)(qbidx[i]);
                bitsmode0 += hessize[tdifidx0[i]+4];
            }
        }
        else
        {
            bitsmode0 = 20000;
        }
    }

    if (flag_pack==0 || (flag_pack==1 && *hLCmode0==1) )
    {
        /* components 0 range : -256~255 */
        max_q = MINIMUM_ENERGY_LOWBRATE;
        min_q = MAXIMUM_ENERGY_LOWBRATE;
        for( i=0; i<BANDS; ++i )
        {
            if (qbidx[i]>max_q)
            {
                max_q = (short)qbidx[i];
            }

            if (qbidx[i]<min_q)
            {
                min_q = (short)qbidx[i];
            }
        }

        /* Counting bits for transmitting all components using same method */
        for(i=0;; ++i)
        {
            if (max_q <= ((2<<(i+1))-1) && min_q >= -(2<<(i+1)))
            {
                break;
            }
        }
        basic_shift = i;

        min_bits=1000;
        min_bits_pos = basic_shift;
        for(offset0=basic_shift; offset0<basic_shift+3; offset0++)
        {
            max_q = MINIMUM_ENERGY_LOWBRATE;
            min_q = MAXIMUM_ENERGY_LOWBRATE;

            bitsmode1 = BITS_DE_8SMODE + BITS_MAX_DEPTH;
            for(i=0; i<BANDS; ++i)
            {
                bitsmode1 += (hessize[((short)(qbidx[i])>>offset0)+4]+(offset0));
            }

            if (min_bits>bitsmode1)
            {
                min_bits_pos = offset0;
                min_bits = bitsmode1;
            }
        }

        bitsmode1 = min_bits;
        lsbdepth1 = min_bits_pos;

        for( i=0; i<BANDS; ++i )
        {
            LSB1[i] = (short)(qbidx[i])&((1<<lsbdepth1)-1);
            tdifidx1[i] = (short)(qbidx[i])>>lsbdepth1;
        }
    }

    /* Preparing lossless coding input */
    if ( flag_pack == 0 )
    {
        /* estimating # of bits */
        /* Encoding MSB bits */
        if (bitsmode0<bitsmode1)
        {
            bits = bitsmode0;
            *hLCmode0 = 0;
        }
        else
        {
            bits = bitsmode1;
            *hLCmode0 = 1;
        }
    }
    else
    {
        /* Encoding MSB bits */
        if (*hLCmode0==0)
        {
            push_indice(st, IND_HQ2_DENG_8SMODE, 0, BITS_DE_8SMODE);
            bits = BITS_DE_8SMODE;
            if (cnt_outlyer0 == 0)
            {
                push_indice(st, IND_HQ2_DENG_8SMODE_N0, 0, BITS_DE_8SMODE_N0);
                bits += BITS_DE_8SMODE_N0;
                if (cnt_outlyer == 1)
                {
                    /* 01 */
                    push_indice(st, IND_HQ2_DENG_8SMODE_N1, 1, BITS_DE_8SMODE_N1);
                    bits+=BITS_DE_8SMODE_N1;
                    push_indice(st, IND_HQ2_DENG_8SPOS, pos_outlyer, BITS_DE_8SPOS);
                    bits+=BITS_DE_8SPOS;
                    push_indice(st, IND_HQ2_DIFF_ENERGY, qbidx[pos_outlyer]+ABS_ENG_OFFSET, BITS_ABS_ENG);
                    bits+=BITS_ABS_ENG;
                }
                else
                {
                    /* 00 */
                    push_indice(st, IND_HQ2_DENG_8SMODE_N1, 0, BITS_DE_8SMODE_N1);
                    bits+=BITS_DE_8SMODE_N1;
                }

                for( i=0; i<pos_outlyer; ++i )
                {
                    push_indice(st, IND_HQ2_DIFF_ENERGY, hescode[tdifidx0[i]+4], hessize[tdifidx0[i]+4]);
                    bitsmode0 += hessize[tdifidx0[i]+4];
                }

                for( i=pos_outlyer+1; i<BANDS; ++i )
                {
                    push_indice(st, IND_HQ2_DIFF_ENERGY, hescode[tdifidx0[i]+4], hessize[tdifidx0[i]+4]);
                    bitsmode0 += hessize[tdifidx0[i]+4];
                }
            }
            else if (cnt_outlyer0 == 1)
            {
                push_indice(st, IND_HQ2_DENG_8SMODE_N0, 1, BITS_DE_8SMODE_N0);
                bits += BITS_DE_8SMODE_N0;
                if (cnt_outlyer == 1)
                {
                    push_indice(st, IND_HQ2_DENG_8SMODE_N1, 1, BITS_DE_8SMODE_N1);
                    bits+=BITS_DE_8SMODE_N1;
                    push_indice(st, IND_HQ2_DENG_8SPOS, pos_outlyer, BITS_DE_8SPOS);
                    bits+=BITS_DE_8SPOS;
                    push_indice(st, IND_HQ2_DIFF_ENERGY, qbidx[0]+ABS_ENG_OFFSET, BITS_ABS_ENG);
                    bits += BITS_ABS_ENG;
                    push_indice(st, IND_HQ2_DIFF_ENERGY, qbidx[pos_outlyer]+ABS_ENG_OFFSET, BITS_ABS_ENG);
                    bits+=BITS_ABS_ENG;
                }
                else
                {
                    push_indice(st, IND_HQ2_DENG_8SMODE_N1, 0, BITS_DE_8SMODE_N1);
                    bits+=BITS_DE_8SMODE_N1;
                    push_indice(st, IND_HQ2_DIFF_ENERGY, qbidx[0]+ABS_ENG_OFFSET, BITS_ABS_ENG);
                    bits += BITS_ABS_ENG;
                }

                for( i=1; i<pos_outlyer; ++i )
                {
                    push_indice(st, IND_HQ2_DIFF_ENERGY, hescode[tdifidx0[i]+4], hessize[tdifidx0[i]+4]);
                    bits += hessize[tdifidx0[i]+4];
                }

                for( i=pos_outlyer+1; i<BANDS; ++i )
                {
                    push_indice(st, IND_HQ2_DIFF_ENERGY, hescode[tdifidx0[i]+4], hessize[tdifidx0[i]+4]);
                    bits += hessize[tdifidx0[i]+4];
                }
            }
        }
        else
        {
            bits = BITS_DE_8SMODE + BITS_MAX_DEPTH;
            push_indice(st, IND_HQ2_DENG_8SMODE, 1, BITS_DE_8SMODE);
            push_indice(st, IND_HQ2_DENG_8SDEPTH, lsbdepth1, BITS_MAX_DEPTH);

            for(i=0; i<BANDS; ++i)
            {
                push_indice(st, IND_HQ2_DIFF_ENERGY, hescode[tdifidx1[i]+4], hessize[tdifidx1[i]+4]);
                bits += hessize[tdifidx1[i]+4];
            }

            if (lsbdepth1 > 0)
            {
                for(i=0; i<BANDS; ++i)
                {
                    push_indice(st, IND_HQ2_DIFF_ENERGY, LSB1[i], lsbdepth1);
                }
                bits += BANDS * lsbdepth1;
            }
        }
    }

    return bits;  /* xx bits for diff. energies + 1 bit for LC coding mode */
}

/*-------------------------------------------------------------------*
 * band_energy_quant()
 *
 *
 *-------------------------------------------------------------------*/

static float band_energy_quant(
    Encoder_State *st,             /* i/o: encoder state structure      */
    const float *t_audio,
    const short band_start[],
    const short band_end[],
    float band_energy[],
    const short bands,
    const Word32 L_qint,  /* Q29 */
    const Word16 eref_fx, /* Q10 */
    const short is_transient
)
{
    short i, k;
    float E;
    short ebits;
    short hLCmode0,hLCmode1,deng_bits;
    int bq1_temp[BANDS_MAX],bq2_temp[BANDS_MAX];
    int bq0;
    int bq1[BANDS_MAX];
    int bq2[BANDS_MAX];
    short deng_cmode = 0;
    short hbits;
    Word32 L_tmp;
    Word32 L_band_energy[BANDS_MAX];
    Word16 exp_normd;
    Word16 rev_qint_fx;
    Word16 Qrev_qint;

    /* Calculate the band energies */
    for (k = 0; k < bands; k++)
    {
        E = 0.0f;
        for (i = band_start[k]; i <= band_end[k]; i++)
        {
            E += t_audio[i] * t_audio[i];
        }

        band_energy[k] = (float) log2_f (E + 0.18e-1f);
    }

    if (is_transient)
    {
        reverse_transient_frame_energies( band_energy, bands );
    }

    /* Quantize the reference and band energies */
    for(k=0; k<bands; k++) L_band_energy[k] = (Word32)(band_energy[k]*pow(2.0f, SWB_BWE_LR_Qbe));

    exp_normd = norm_l(L_qint);
    rev_qint_fx = div_s(0x4000, round_fx(L_shl(L_qint, exp_normd))); /* Q14-(29+exp_normd-16)+15 */
    Qrev_qint = sub(14-(29-16)+15, exp_normd);

    bq0 = (int)round_fx(L_shl(L_mult(eref_fx, rev_qint_fx), sub(5, Qrev_qint))); /* 16-(10+Qrev_qint+1) */
    FOR (k = 0; k < bands; k++)
    {
        /*bq1[k] = round_f (band_energy[k] / qint); */
        L_tmp = L_mls(L_band_energy[k], rev_qint_fx); /* Q14+Qrev_qint-15 */
        bq1[k] = (int)round_fx( L_shl(L_tmp, sub(17, Qrev_qint)) ); /* 16-(14+Qrev_qint-15) */
    }

    if(is_transient)
    {
        mvi2i( bq1, bq1_temp, bands );

        /* Calculate the differential energies */
        diffcod_lrmdct(bands,bq0,bq1_temp,bq2_temp,is_transient);
    }

    /* Calculate the differential energies */
    bq2[0] = bq1[0] - bq0;
    for (k = 1; k < bands; k++)
    {
        bq2[k] = bq1[k] - bq1[k - 1];
    }

    /* Modifying qbidx to be located in the range -256~255 */
    for( i=0; i<bands; ++i )
    {
        if (bq2[i]>MAXIMUM_ENERGY_LOWBRATE)
        {
            bq2[i] = MAXIMUM_ENERGY_LOWBRATE;
        }
        if (bq2[i]<MINIMUM_ENERGY_LOWBRATE)
        {
            bq2[i] = MINIMUM_ENERGY_LOWBRATE;
        }
    }

    /* Get number of bits by Huffman0 coding */
    ebits = large_symbol_enc( st, bq2, bands, &hLCmode0, 0 );

    if(is_transient)
    {
        /* Get number of bits by Huffman coding */
        hbits = small_symbol_enc_tran(st, bq2_temp, bands, &hLCmode1, 0 ,is_transient);
    }
    else
    {
        /* Get number of bits by Huffman coding */
        hbits = small_symbol_enc( st, bq2, bands, &hLCmode1, 0, is_transient );
    }

    /* comparing used bits */
    if ( ebits < hbits || hbits == -1 )
    {
        deng_cmode = 0;
        push_indice ( st, IND_HQ2_DENG_MODE, deng_cmode, BITS_DE_CMODE );
        large_symbol_enc( st, bq2, bands, &hLCmode0, 1 );
        deng_bits = ebits+BITS_DE_CMODE;
    }
    else
    {
        /* setting energy difference coding mode and storing it */
        deng_cmode = 1;
        push_indice( st, IND_HQ2_DENG_MODE, deng_cmode, BITS_DE_CMODE );

        deng_bits = hbits + BITS_DE_CMODE;

        /* packing indice */
        if(is_transient)
        {
            mvi2i( bq2_temp, bq2, bands );

            small_symbol_enc_tran(st, bq2, bands, &hLCmode1, 1 ,is_transient);
        }
        else
        {
            small_symbol_enc( st, bq2, bands, &hLCmode1, 1, is_transient );
        }
    }

    /* Reconstruct quantized spectrum */
    bq1[0] = bq2[0] + bq0;
    for (k = 1; k < bands; k++)
    {
        bq1[k] = bq2[k] + bq1[k - 1];
    }

    for (k = 0; k < bands; k++)
    {
        L_band_energy[k] = L_mls(L_qint, (Word16)bq1[k]);
        move32();/* 29+0-15 -> Qbe(Q14) */
        band_energy[k] = (float)(L_band_energy[k]/pow(2.0f, SWB_BWE_LR_Qbe));
    }

    if( is_transient )
    {
        reverse_transient_frame_energies( band_energy, bands );
    }

    return( deng_bits );
}


/*-------------------------------------------------------------------*
 * p2a_threshold_quant()
 *
 *
 *-------------------------------------------------------------------*/

static short p2a_threshold_quant(
    Encoder_State *st,             /* i/o: encoder state structure      */
    const float *t_audio,
    const short band_start[],
    const short band_end[],
    const short band_width[],
    const short bands,
    const short p2a_bands,
    const float p2a_th,
    short *p2a_flags
)
{
    short i, j, k;
    float p, a, e;
    float p2a;

    set_s( p2a_flags, 1, bands );

    j = 0;
    for (k = bands - p2a_bands; k < bands; k++)
    {
        a = 0.0f;
        p = 0.0f;
        for (i = band_start[k]; i <= band_end[k]; i++)
        {
            e = t_audio[i] * t_audio[i];
            if (e > p)
            {
                p = e;
            }
            a += e;
        }

        if (a > 0.0f)
        {
            a /= band_width[k];
            p2a = 10.0f * (float) log10 (p / a);

            if (p2a <= p2a_th)
            {
                p2a_flags[k] = 0;
            }
        }

        push_indice( st, IND_HQ2_P2A_FLAGS, p2a_flags[k], 1 );
        j++;
    }

    return( j );
}


/*-------------------------------------------------------------------*
 * mdct_spectrum_fine_gain_enc()
 *
 *
 *-------------------------------------------------------------------*/

static void mdct_spectrum_fine_gain_enc(
    Encoder_State *st,             /* i/o: encoder state structure      */
    const float ybuf[],
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
    float Eyy, Exy, gamma;
    float dmin, d;
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

    for (k = bands - Ngq; k < bands; k++)
    {
        Eyy = 0.0f;
        Exy = 0.0f;
        for (i = band_start[k_sort[k]]; i <= band_end[k_sort[k]]; i++)
        {
            Eyy += y2[i] * y2[i];
            Exy += ybuf[i] * y2[i];
        }

        if (Eyy > 0.0f && Exy > 0.0f)
        {
            gamma = Exy / Eyy;
            dmin = FLT_MAX;
            imin = -1;
            for (i = 0; i < gqlevs; i++)
            {
                d = (float) fabs (gamma - gain_table[i]);
                if (d < dmin)
                {
                    dmin = d;
                    imin = i;
                }
            }

            gamma = gain_table[imin];

            for (i = band_start[k_sort[k]]; i <= band_end[k_sort[k]]; i++)
            {
                y2[i] *= gamma;
            }
        }
        else
        {
            imin = 0;
        }

        push_indice( st, IND_HQ2_SUBBAND_GAIN, imin, gqbits );
    }

    return;
}
