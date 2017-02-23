
/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <math.h>
#include "options.h"
#include "rom_com.h"
#include "prot.h"
#include "prot.h"
#include "stl.h"

/*-------------------------------------------------------------------*
* Local functions
*--------------------------------------------------------------------*/

static short calc_pvq_splits( Encoder_State *st, const short band_bits, const short sfmsize, const float *y, short *bits );
static void obtainEnergyParameter( float Enear, float Eopp, short *param);

static void densityIndexSymbolEncode( Encoder_State *st, short density , short  r_dim, short l_dim , short index_phi);

/*-------------------------------------------------------------------*
* pvq_encode_band()
*
* Encode band with PVQ
*--------------------------------------------------------------------*/

static void pvq_encode_band(
    Encoder_State *st,
    const float *coefs_norm,
    short *pulse_vector,
    short *npulses,
    float *coefs_quant,
    const short sfmsize,
    const short band_bits,
    short *bits_left,
    const short strict_bits
)
{
    short K_val;
    short j, Np;
    float enr, E_part[MAX_SPLITS+1];
    short part_start[MAX_SPLITS+1], dim_part[MAX_SPLITS+1], bits_part[MAX_SPLITS+1];
    short pool_tot, pool_part, dim_parts;
    float g_part[MAX_SPLITS];
    short g_part_s[MAX_SPLITS];
    short sg_part[MAX_SPLITS+1];
    short idx_sort[MAX_SPLITS+1];
    short js, band_bits_tot, split_bit;

    Np = calc_pvq_splits(st, band_bits, sfmsize, coefs_norm, &split_bit);
    band_bits_tot = band_bits - split_bit;

    enr = 0.0f;
    for(j = 0; j<sfmsize; j++)
    {
        enr += coefs_norm[j]*coefs_norm[j];
    }

    dim_parts = intLimCDivPos( (int)sfmsize,Np) ;
    set_s(dim_part,dim_parts,Np-1);
    dim_part[Np-1] = sfmsize-dim_parts*(Np-1);

    part_start[0] = 0;
    for(j = 1; j<Np; j++)
    {
        part_start[j] = part_start[j-1] + dim_part[j-1];
    }

    set_s( g_part_s, -32768, Np );
    if(Np > 1)
    {
        encode_energies( st, coefs_norm, Np, dim_part, E_part, bits_part, g_part_s, band_bits_tot, bits_left, enr, sfmsize, strict_bits );
    }
    else
    {
        bits_part[0] = band_bits_tot;
    }

    pool_tot = 0;
    pool_part = 0;

    for (j = 0; j < Np; j++)
    {
        g_part[j] = -((float)g_part_s[j])/32768;
        g_part_s[j] = -g_part_s[j];
    }

    srt_vec_ind(g_part_s,sg_part,idx_sort,Np);
    for(j = 0; j < Np; j++)
    {
        js = idx_sort[Np-1-j];


        pool_part     = shrtCDivSignedApprox(pool_tot, Np-j);
        bits_part[js] = max(0, min(bits_part[js]+pool_part, 256));

        conservativeL1Norm(dim_part[js],bits_part[js], strict_bits, *bits_left, pool_tot , *npulses, /* inputs */
                           &K_val, bits_left, &pool_tot, npulses);                                   /* outputs */
        if( K_val >= 1 )
        {
            pvq_encode(st, coefs_norm + part_start[js], pulse_vector + part_start[js],
                       coefs_quant + part_start[js], K_val, dim_part[js], g_part[js]);
        }
        else
        {
            set_f(coefs_quant + part_start[js],0.0f,dim_part[js]);
            set_s(pulse_vector + part_start[js],0,dim_part[js]);
        }
    }

    return;
}

/*-------------------------------------------------------------------*
* pvq_encode_frame()
*
*
*--------------------------------------------------------------------*/

void pvq_encode_frame(
    Encoder_State *st,
    const float *coefs_norm,       /* i  : normalized coefficients to encode */
    float *coefs_quant,      /* o  : quantized coefficients */
    float *gopt,             /* o  : optimal shape gains */
    short *npulses,          /* o  : number of pulses per band */
    short *pulse_vector,     /* o  : non-normalized pulse shapes */
    const short *sfm_start,        /* i  : indices of first coefficients in the bands */
    const short *sfm_end,          /* i  : indices of last coefficients in the bands */
    const short *sfmsize,          /* i  : band sizes */
    const short nb_sfm,            /* i  : total number of bands */
    const short *R,                /* i  : bitallocation per band (Q3)*/
    const short pvq_bits,          /* i  : number of bits avaiable */
    const short core               /* i  : core */
)
{
    short i, j;
    short band_bits, bits_left;
    short bit_pool = 0;
    short coded_bands, bands_to_code;
    short curr_bits;
    short R_sort[NB_SFM]; /*Q3*/
    short is, i_sort[NB_SFM];
    short strict_bits;

    rc_enc_init(st, pvq_bits);
    curr_bits = (pvq_bits - RC_BITS_RESERVED)<<3;

    bands_to_code = 0;
    for (i = 0; i < nb_sfm; i++)
    {
        if (R[i] > 0)
        {
            bands_to_code++;
        }
    }

    if (core == ACELP_CORE)
    {
        strict_bits = 1;
        srt_vec_ind (R, R_sort, i_sort, nb_sfm);
    }
    else
    {
        strict_bits = 0;
        for(i=0; i<nb_sfm; i++)
        {
            i_sort[i] = i;
        }
    }

    coded_bands = 0;
    for (i = 0; i < nb_sfm; i++)
    {
        is = i_sort[i];
        gopt[is] = 0;
        if(R[is] > 0)
        {
            bandBitsAdjustment(st->rc_num_bits, st->rc_range, curr_bits, bands_to_code, bands_to_code-coded_bands, sfmsize[is] ,R[is], bit_pool, /* inputs  */
                               &band_bits, &bits_left, &bit_pool);                                                                               /* outputs */
            pvq_encode_band( st, &coefs_norm[sfm_start[is]], &pulse_vector[sfm_start[is]],
                             &npulses[is], &coefs_quant[sfm_start[is]], sfmsize[is], band_bits,
                             &bits_left, strict_bits);
            gopt[is] = dotp(coefs_quant+sfm_start[is], coefs_norm+sfm_start[is], sfmsize[is]) /
                       (dotp(coefs_quant+sfm_start[is], coefs_quant+sfm_start[is], sfmsize[is]) + 1e-15f);
            if (gopt[is] == 0.0f)
            {
                gopt[is] = 1e-10f;
            }
            /* Updates */
            coded_bands++;
        }
        else
        {
            for (j = sfm_start[is]; j < sfm_end[is]; j++)
            {
                coefs_quant[j] = 0.0f;
                pulse_vector[j] = 0;
            }
        }
    }

    rc_enc_finish(st);

    return;
}

/*---------------------------------------------------------------------*
 * pvq_core_enc()
 *
 * Main Generic Audio Encoder Routine
 *---------------------------------------------------------------------*/

short pvq_core_enc (
    Encoder_State *st,
    float coefs_norm[],
    float coefs_quant[],
    short bits_tot,                           /* total number of bits */
    short nb_sfm,
    const short *sfm_start,
    const short *sfm_end,
    const short *sfmsize,
    short *R,
    short *Rs,
    short *npulses,
    short *maxpulse,
    const short core
)
{
    short i;
    short  R_upd; /*Q3*/
    short  ord[NB_SFM_MAX];
    float  fg_pred[NB_SFM_MAX];
    short  pvq_bits;
    short  pulse_vector[L_FRAME48k];
    float gopt[NB_SFM];
    short gain_bits_array[NB_SFM];
    short gain_bits_tot;

    R_upd = bits_tot * 8;
    gain_bits_tot = assign_gain_bits( core, nb_sfm, sfmsize, R, gain_bits_array, &R_upd );
    pvq_bits = R_upd >> 3;
    pvq_encode_frame( st, coefs_norm, coefs_quant, gopt, npulses, pulse_vector,
                      sfm_start, sfm_end, sfmsize, nb_sfm, R, pvq_bits, core );
    bits_tot = pvq_bits + gain_bits_tot;



    if( Rs != NULL )
    {
        for(i=0; i<nb_sfm; i++)
        {
            Rs[i]   = Rs[i] * (npulses[i] > 0); /* Update Rs in case no pulses were assigned */
        }
    }

    for(i=0; i<nb_sfm; i++)
    {
        ord[i] = i;
        R[i]   = R[i] * (npulses[i] > 0); /* Update in case no pulses were assigned */
    }

    get_max_pulses( sfm_start, sfm_end, ord, npulses, nb_sfm, pulse_vector, maxpulse );

    fine_gain_pred( sfm_start, sfm_end, sfmsize, ord, npulses, maxpulse, R, nb_sfm,
                    coefs_quant, pulse_vector, fg_pred, core);

    fine_gain_quant(st, ord, nb_sfm, gain_bits_array, fg_pred, gopt);

    apply_gain(ord, sfm_start, sfm_end, nb_sfm, fg_pred, coefs_quant);

    return bits_tot;
}


/*-------------------------------------------------------------------*
* encode_energies()
*
*
*--------------------------------------------------------------------*/

void encode_energies(
    Encoder_State *st,
    const float *coefs,
    short Np,
    short *dim_part,
    float *E_part,
    short *bits_part,
    short *g_part,
    short qband,

    short *bits_left,
    float enr,
    short dim,
    const short strict_bits
)
{
    short i, j, l_Np, r_Np;
    short l_bits, r_bits, l_dim, r_dim;
    float l_enr, r_enr;
    short il, ir;
    short oppRQ3, qzero;
    short density;
    short phi;
    short index_phi=-1;

    l_Np = Np>>1;
    r_Np = Np-l_Np;

    l_enr = 0.0f;
    l_bits = 0;
    l_dim = 0;
    for(i=0; i<l_Np; i++)
    {
        l_dim += dim_part[i];
    }
    for(j = 0; j<l_dim; j++)
    {
        l_enr += coefs[j]*coefs[j];
    }
    r_enr = enr - l_enr;
    r_dim = dim - l_dim;

    obtainEnergyQuantizerDensity( dim, qband, &density);
    obtainEnergyParameter( l_enr, r_enr, &phi);
    rangeCoderFinalizationFBits(st->rc_num_bits, st->rc_range, &qzero);
    densityAngle2RmsProjEnc(density, phi , &index_phi, &ir, &il, &oppRQ3);
    densityIndexSymbolEncode(st, density, r_dim, l_dim, index_phi);

    for(i = 0; i<l_Np; i++)
    {
        g_part[i] = ((int)g_part[i] * il + 16384) >> 15;
    }

    for(i = l_Np; i<Np; i++)
    {
        g_part[i] = ((int)g_part[i] * ir + 16384) >> 15;
    }
    NearOppSplitAdjustment( qband, qzero, st->rc_num_bits, st->rc_range,  *bits_left,
                            strict_bits, Np, dim_part[0], dim_part[Np-1],
                            l_dim, r_dim, oppRQ3,
                            &l_bits, &r_bits, bits_left);
    if(l_Np > 1)
    {
        encode_energies( st, coefs, l_Np, dim_part, E_part, bits_part, g_part, l_bits, bits_left, l_enr, l_dim, strict_bits );
    }
    else
    {
        E_part[0] = l_enr;
        bits_part[0] = l_bits;
    }
    if(r_Np > 1)
    {
        encode_energies( st, &coefs[l_dim], r_Np, &dim_part[l_Np], &E_part[l_Np], &bits_part[l_Np], &g_part[l_Np], r_bits, bits_left, r_enr, r_dim, strict_bits );
    }
    else
    {
        E_part[1] = r_enr;
        bits_part[1] = r_bits;
    }
    return;
}

static void densityIndexSymbolEncode(Encoder_State *st,
                                     short density ,
                                     short opp_sz,
                                     short near_sz ,
                                     short index_phi)
{
    short angle, c ;
    int sym_freq = 1, cum_freq, tot;
    short densityPlOne , densitySubC;
    short densitySubIndex,nearFlag;

    if ((0xFFFE&density)  != 0  )  /* even */
    {
        angle = atan2_fx(SQRT_DIM_fx[opp_sz], SQRT_DIM_fx[near_sz]);
        angle = shl(angle, 1);
        angle = mult_r(angle, 20861);
        c = mult_r(density, angle);

        densityPlOne    = 1 + density;
        densitySubC     = density - c;

        tot      = 2 + density*densityPlOne;           /* c==0, c==density*/
        sym_freq = 1 + 2*index_phi ;                   /* c==density */
        cum_freq = index_phi*index_phi;                /* c==density*/

        if (c == 0)
        {
            sym_freq = 2*densityPlOne  - sym_freq;
            cum_freq = 2*index_phi*(densityPlOne) - cum_freq;
        }
        else if( densitySubC != 0  ) /* c  n.eq. density */
        {
            densitySubIndex = density - index_phi;
            nearFlag        = (index_phi <= c);

            tot = densityPlOne + density*c*(densitySubC);
            sym_freq = nearFlag ? ( 1 + 2*index_phi*densitySubC ) : ( 1 + 2*densitySubIndex*c );
            cum_freq = nearFlag ? ( index_phi*((index_phi-1)*densitySubC + 1)) :
                       ( tot - densityPlOne - densitySubIndex*(densitySubIndex+1)*c + index_phi);
        }
        /* else keep values for c==density*/
        rc_encode(st, cum_freq, sym_freq, tot);
    }

    return;
}


/*--------------------------------------------------------------------------*
 * calc_pvq_splits()
 *
 * Calculate the number of segments needed
 *--------------------------------------------------------------------------*/

static short calc_pvq_splits(          /* o  : Number of segments           */
    Encoder_State *st,                 /* i/o: Encoder state                */
    const short band_bits,             /* i  : Band bit rate                */
    const short sfmsize,               /* i  : Band width                   */
    const float *y,                    /* i  : Target vector                */
    short *bits                  /* o  : Consumed bits                */
)
{
    short Np;
    short Npart;
    short i,j,k;
    float E[MAX_SPLITS];
    float Emean;
    float tmp;
    float max_dev;
    Np = (short)(intLimCDivPos( (int)band_bits, 67)>>2);
    if  (band_bits - 268*Np != 0 || Np == 0 ) /* L_msu */
    {
        Np++;            /* ceil */
    }
    *bits = 0;

    if ( Np < MAX_SPLITS && (band_bits - (8*sfmsize * THR_ADD_SPLIT ) > 0))
    {
        Npart  =   (short) intLimCDivPos( (int)sfmsize, Np);
        *bits = 8;
        Emean = 0;
        k = 0;
        for ( i = 0; i < Np; i++)
        {
            E[i] = EPSILON;
            for (j = 0; j < Npart; j++, k++)
            {
                E[i] += y[k]*y[k];
            }
            E[i]   =  30 - norm_l(max((int)E[i], 1)); /* a 0 input integer  Ei yields a zeroed log2(Ei) out */
            Emean += E[i];
        }
        Emean /= Np;

        max_dev = -1;
        for ( i = 0; i < Np; i++)
        {
            tmp = (float)fabs(E[i] - Emean);
            if ( tmp > max_dev )
            {
                max_dev = tmp;
            }
        }

        if ( max_dev > (32 - band_bits/(8*Np)) )
        {
            rc_enc_bits(st, 1, 1);
            Np += 1;
        }
        else
        {
            rc_enc_bits(st, 0, 1);
        }
    }

    Np = max(Np, (short)(ceil((float)sfmsize/PVQ_MAX_BAND_SIZE)));
    Np = min(MAX_SPLITS, Np);
    Np = min((short)floor((float)sfmsize/MIN_BAND_SIZE), Np);
    return Np;
}


#define EPSILON_obtainEnergyParameter 0.00373f /**/
#define C_obtainEnergyParameter ((float)(1 << 15) / EVS_PI + EPSILON_obtainEnergyParameter)

void obtainEnergyParameter(
    float Enear,
    float Eopp,
    short *param
)
{
    if (!Eopp)
    {
        *param=0;
        return;
    }
    if (!Enear)
    {
        *param = 1 << 14;
        return;
    }
    *param = (short) (C_obtainEnergyParameter  * atan2( (float)sqrt(Eopp) , (float) sqrt(Enear)) + 0.5f);
    return;
}
#undef C_obtainEnergyParameter
#undef EPSILON_obtainEnergyParameter


