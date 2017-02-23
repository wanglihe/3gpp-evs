/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"

/*-----------------------------------------------------------------*
 * Local functions
 *-----------------------------------------------------------------*/

static void make_offset_scale(int j, const unsigned int tab_no_cv[], const int no_ld[], short no_scl, unsigned int offset_scale[][MAX_NO_SCALES+1]);
static void init_offset( unsigned int offset_scale1[][MAX_NO_SCALES+1], unsigned int offset_scale2[][MAX_NO_SCALES+1], unsigned int offset_scale1_p[][MAX_NO_SCALES+1],
                         unsigned int offset_scale2_p[][MAX_NO_SCALES+1], short no_scales[][2], short no_scales_p[][2]);
static void decode_comb( int index, float *cv, int idx_lead );
static void decode_sign_pc1( float *c, int idx_sign, int parity );
static void put_value( float *cv, int *p, float val, int dim, int no_new_val );
static void decode_leaders( int index, int idx_lead, float *cv );
static void idx2c( int n, int *p, int k, int val );
static void divide_64_32( short *x, unsigned int y, unsigned int *result, unsigned int *rem );
static short decode_indexes( short *index, int no_bits, const float *p_scales, short *p_no_scales,
                             unsigned int *p_offset_scale1, unsigned int *p_offset_scale2, float *x_lvq, short mode_glb );

/*-----------------------------------------------------------------*
 * permute()
 * used in CNG-LP coding
 *-----------------------------------------------------------------*/

void permute(
    float *pTmp1,         /* i/o: vector whose components are to be permuted */
    const short *perm     /* i  : permutation info (indexes that should be interchanged), max two perms */
)
{
    int p1, p2;
    float tmp;

    p1 = perm[0];
    p2 = perm[1];
    tmp = pTmp1[p1];
    pTmp1[p1] = pTmp1[p2];
    pTmp1[p2] = tmp;
    p1 = perm[2];

    if ( p1 > -1 )
    {
        p2 = perm[3];
        tmp = pTmp1[p1];
        pTmp1[p1] = pTmp1[p2];
        pTmp1[p2] = tmp;
    }

    return;
}

/*-----------------------------------------------------------------*
 * init_lvq()
 *
 *-----------------------------------------------------------------*/

void init_lvq(
    unsigned int offset_scale1[][MAX_NO_SCALES+1],
    unsigned int offset_scale2[][MAX_NO_SCALES+1],
    unsigned int offset_scale1_p[][MAX_NO_SCALES+1],
    unsigned int offset_scale2_p[][MAX_NO_SCALES+1],
    short no_scales[][2],
    short no_scales_p[][2]
)
{
    short i, j;
    /* safety-net mode */
    for(i=0; i<MAX_NO_MODES; i++)
    {
        j=0;
        while ((j<MAX_NO_SCALES) && (no_lead[i][j] >0 ))
        {
            j++;
        }
        no_scales[i][0] = j;
        j = MAX_NO_SCALES;
        while ((j<MAX_NO_SCALES*2) && (no_lead[i][j] >0 ))
        {
            j++;
        }
        no_scales[i][1] = j-MAX_NO_SCALES;
    }
    /* predictive mode */
    for(i=0; i<MAX_NO_MODES_p; i++)
    {
        j=0;
        while ((j<MAX_NO_SCALES) && (no_lead_p[i][j] >0 ))
        {
            j++;
        }
        no_scales_p[i][0] = j;
        j = MAX_NO_SCALES;
        while ((j<MAX_NO_SCALES*2) && (no_lead_p[i][j] >0 ))
        {
            j++;
        }
        no_scales_p[i][1] = j-MAX_NO_SCALES;
    }
    /* index offsets for each truncation */
    init_offset( offset_scale1, offset_scale2, offset_scale1_p, offset_scale2_p, no_scales, no_scales_p );
}

/*-----------------------------------------------------------------*
 * make_offset_scale()
 *
 *-----------------------------------------------------------------*/
static
void make_offset_scale(
    int j,
    const unsigned int tab_no_cv[],
    const int no_ld[],
    short no_scl,
    unsigned int offset_scale[][MAX_NO_SCALES+1]
)
{
    int i;

    offset_scale[j][0] = 1;
    for( i=1; i<=no_scl; i++ )
    {
        offset_scale[j][i] = offset_scale[j][i-1] + tab_no_cv[no_ld[i-1]];
    }

    return;
}

/*-----------------------------------------------------------------*
 * init_offset()
 *
 *-----------------------------------------------------------------*/
static
void init_offset(
    unsigned int offset_scale1[][MAX_NO_SCALES+1],
    unsigned int offset_scale2[][MAX_NO_SCALES+1],
    unsigned int offset_scale1_p[][MAX_NO_SCALES+1],
    unsigned int offset_scale2_p[][MAX_NO_SCALES+1],
    short no_scales[][2],
    short no_scales_p[][2]
)
{
    int j;
    /* safety-net */
    for( j=0; j<MAX_NO_MODES; j++ )
    {
        make_offset_scale( j, table_no_cv, no_lead[j], no_scales[j][0], offset_scale1 );
        make_offset_scale( j, table_no_cv, &no_lead[j][MAX_NO_SCALES], no_scales[j][1], offset_scale2 );
    }

    /* predictive modes AR and MA */
    for( j=0; j<MAX_NO_MODES_p; j++ )
    {
        make_offset_scale(j, table_no_cv, no_lead_p[j], no_scales_p[j][0], offset_scale1_p);
        make_offset_scale(j, table_no_cv, &no_lead_p[j][MAX_NO_SCALES], no_scales_p[j][1], offset_scale2_p);
    }

    offset_scale1[MAX_NO_MODES][0] = 1;
    offset_scale2[MAX_NO_MODES][0] = 1;
    offset_scale1_p[MAX_NO_MODES_p][0] = 1;
    offset_scale2_p[MAX_NO_MODES_p][0] = 1;

    return;
}


/*-----------------------------------------------------------------*
 * decode_indexes()
 *
 *-----------------------------------------------------------------*/

static short decode_indexes(
    short * index,
    int no_bits,
    const float * p_scales,
    short * p_no_scales,
    unsigned int * p_offset_scale1,
    unsigned int * p_offset_scale2,
    float * x_lvq,
    short mode_glb
)
{
    unsigned int index1=0, index2, i, idx_scale;
    short len_scales = MAX_NO_SCALES*2, no_modes;
    float scale;

    no_modes = MAX_NO_SCALES+1;

    if (no_bits <= 2*LEN_INDICE) /* the third short is not used */
    {
        index[2] = 0;
        if ( no_bits <= LEN_INDICE )
        {
            index[1] =0;
        }
    }

    /* safety check in case of bit errors */
    for( i = 0; i<3; i++ )
    {
        if( index[i] < 0 )
        {
            set_f( x_lvq, 0.0f, 2*LATTICE_DIM );
            index[i] = 0;
            return 1;
        }
    }

    /* first subvector */
    if ( p_offset_scale2[mode_glb*no_modes + p_no_scales[mode_glb*2+1]] > 0 )
    {
        divide_64_32( index, (unsigned int)p_offset_scale2[mode_glb*no_modes+p_no_scales[mode_glb*2+1]], &index1, &index2 );
    }
    else
    {
        index1 = (unsigned int)(index[0]); /* this is for very low bitrates, so there is no loss in truncation */
        index2 = 0;
    }

    if ( index1 == 0 )
    {
        for( i=0; i<LATTICE_DIM; i++ )
        {
            x_lvq[i] = 0.0;
        }
    }
    else
    {
        if( index1 >= p_offset_scale1[mode_glb*no_modes+p_no_scales[mode_glb*2]] )
        {
            /* safety check in case of bit errors */
            set_f( x_lvq, 0.0f, 2*LATTICE_DIM );
            return 1;
        }

        /* find idx_scale */
        i = 1;
        while( (short)i <= p_no_scales[mode_glb*2] && index1 >= p_offset_scale1[mode_glb*no_modes +i] )
        {
            i++;
        }

        idx_scale = i-1;
        index1 -= p_offset_scale1[mode_glb*no_modes+idx_scale];

        /* find idx_leader */
        i = 1;
        while( index1 >= table_no_cv[i] )
        {
            i++;
        }
        decode_comb(index1-table_no_cv[i-1], x_lvq, i-1 );
        scale = p_scales[mode_glb*len_scales+idx_scale];
        for( i=0; i<LATTICE_DIM; i++ )
        {
            x_lvq[i] *= scale;
        }
    }

    /* second subvector */
    if ( index2 == 0 )
    {
        for( i=LATTICE_DIM; i<2*LATTICE_DIM; i++ )
        {
            x_lvq[i] = 0.0;
        }
    }
    else
    {

        /* find the index for the scale/truncation */
        i = 1;
        while( index2 >= p_offset_scale2[mode_glb*no_modes+i] )
        {
            i++;
        }

        idx_scale = i-1;
        index2 -= p_offset_scale2[mode_glb*no_modes+idx_scale];
        /* find the index of the leader vector */
        i = 1;
        while ( index2 >= table_no_cv[i] )
        {
            i++;
        }

        decode_comb( index2-table_no_cv[i-1], &x_lvq[LATTICE_DIM], i-1 );

        scale = p_scales[mode_glb*len_scales+MAX_NO_SCALES+idx_scale];
        for( i=LATTICE_DIM; i<2*LATTICE_DIM; i++ )
        {
            x_lvq[i] *= scale;
        }
    }

    return 0;
}

/*-----------------------------------------------------------------*
 * deindex_lvq()
 *
 *-----------------------------------------------------------------*/

short deindex_lvq(
    short *index,                  /* i  : index to be decoded, as an array of 3 short */
    float *x_lvq,                  /* o  : decoded codevector                          */
    short mode,                    /* i  : LVQ  coding mode (select scales & no_lead ), or idx_cv */
    short sf_flag,                 /* i  : safety net flag                 */
    short no_bits,                 /* i  : number of bits for lattice      */
    unsigned int *p_offset_scale1, /* i  : offset for first subvector      */
    unsigned int *p_offset_scale2, /* i  : offset for the second subvector */
    short *p_no_scales
)
{
    int i;
    const float * p_scales;
    short mode_glb;
    short ber_flag;

    if ( sf_flag == 1 )
    {
        mode_glb = offset_lvq_modes_SN[mode] + offset_in_lvq_mode_SN[mode][no_bits-min_lat_bits_SN[mode]];
        p_scales = &scales[0][0];
    }
    else
    {
        mode_glb = offset_lvq_modes_pred[mode] + offset_in_lvq_mode_pred[mode][no_bits-min_lat_bits_pred[mode]];
        p_scales = &scales_p[0][0];
    }

    /* decode the lattice index into the lattice codevectors for the two subvectors */
    ber_flag = decode_indexes( index, no_bits, p_scales, p_no_scales, p_offset_scale1, p_offset_scale2, x_lvq, mode_glb );

    if ( sf_flag == 1 )
    {
        for( i=0; i<2*LATTICE_DIM; i++ )
        {
            x_lvq[i] *= sigma[mode][i];
        }
    }
    else
    {
        for( i=0; i<2*LATTICE_DIM; i++ )
        {
            x_lvq[i] *= sigma_p[mode][i];
        }
    }

    return ber_flag;
}

/*------------------------------------------------------------------------------------------------------------*
 * deindex_lvq_cng()
 *
 * Note:
 * The sampling frequency for the LVQ CNG decoder frame can be determined by checking the fully decoded
 * value of the highest order LSF coefficient. Thus sampling rate information, nor extra codebooks are
 * not needed for deindex_lvq_cng(), since it is embedded inside the LSF codebooks.
 *------------------------------------------------------------------------------------------------------------*/

short deindex_lvq_cng(
    short *index,                   /* i  : index to be decoded, as an array of 3 short */
    float *x_lvq,                   /* o  : decoded codevector */
    short idx_cv,                   /* i  : relative mode_lvq, wrt START_CNG */
    int no_bits,                    /* i  : number of bits for lattice */
    unsigned int * p_offset_scale1,
    unsigned int * p_offset_scale2,
    short * p_no_scales
)
{
    int i;
    const float *p_scales;
    short mode_glb, mode;
    short ber_flag;

    mode_glb = START_CNG + idx_cv;
    mode = LVQ_COD_MODES + idx_cv;

    p_scales = &scales[0][0];
    ber_flag = decode_indexes( index, no_bits, p_scales, p_no_scales, p_offset_scale1, p_offset_scale2, x_lvq, mode_glb );

    for( i=0; i<2*LATTICE_DIM; i++ )
    {
        x_lvq[i] *= sigma[mode][i];
    }

    /* check if permutting needed */
    if ( cng_sort[idx_cv] )
    {
        permute( x_lvq, perm[idx_cv] );
    }

    return ber_flag;
}



/*-----------------------------------------------------------------*
 * decode_comb()
 *
 * combinatorial deindexing of a codevector including the signs
*
 *-----------------------------------------------------------------*/
static void decode_comb(
    int index,          /* i  : index to be decoded */
    float *cv,          /* o  : decoded codevector  */
    int idx_lead        /* i  : leader class index  */
)
{
    int idx_sign;

    idx_sign = (int)(index/pi0[idx_lead]);
    index -= idx_sign*pi0[idx_lead];
    decode_leaders((int)index, idx_lead, cv);
    decode_sign_pc1(cv, idx_sign, pl_par[idx_lead]);

    return;
}

/*-----------------------------------------------------------------*
 * decode_leaders()
 *
 * decode index of a codevector from the leader class idx_lead
*-----------------------------------------------------------------*/

static void decode_leaders(
    int index,          /* i  : index to be decoded */
    int idx_lead,       /* i  : leader class index  */
    float *cv           /* o  : decoded codevector  */
)
{
    int i, no_vals_loc, no_vals_last, p[LATTICE_DIM], index1, dim_loc, n_crt;
    float val_crt;

    no_vals_loc = no_vals[idx_lead];
    val_crt = vals[idx_lead][no_vals_loc-1];
    no_vals_last = no_vals_ind[idx_lead][no_vals_loc-1];

    for( i=0; i<no_vals_last; i++ )
    {
        cv[i] = val_crt;
    }

    val_crt = 1;
    dim_loc = no_vals_last;

    switch ( no_vals_loc )
    {
    case 1:
        break;
    case 2:
        idx2c(LATTICE_DIM, p, no_vals_ind[idx_lead][0], index);
        put_value(cv, p, vals[idx_lead][0], no_vals_last, no_vals_ind[idx_lead][0]);
        break;
    case 4:
        dim_loc += no_vals_ind[idx_lead][2];
        n_crt = no_vals_ind[idx_lead][2];
        index1 = (int)index/C[dim_loc][n_crt];
        index -= index1*C[dim_loc][n_crt];
        idx2c(dim_loc, p, n_crt, index);
        put_value(cv, p, vals[idx_lead][2], no_vals_last, no_vals_ind[idx_lead][2]);
        index = index1;
    case 3:
        dim_loc +=  no_vals_ind[idx_lead][1];
        n_crt = no_vals_ind[idx_lead][1];
        index1 = (int)index/C[dim_loc][n_crt];
        index -= index1*C[dim_loc][n_crt];
        idx2c(dim_loc, p, n_crt, index);
        put_value(cv, p, vals[idx_lead][1], dim_loc-n_crt, n_crt);
        idx2c(LATTICE_DIM, p, no_vals_ind[idx_lead][0], index1);
        put_value(cv, p, vals[idx_lead][0], dim_loc, no_vals_ind[idx_lead][0]);
        break;
    }

    return;
}

/*-----------------------------------------------------------------*
 * put_value()
 *
 * inserts no_new_val values of val in the codevector cv at the positions given by the array p
 *-----------------------------------------------------------------*/

static void put_value(
    float *cv,                /* i  : input codevector */
    int *p,                   /* i  : array with positions */
    float val,                /* i  : value to be inserted */
    int dim,                  /* i  : vector dimension  */
    int no_new_val            /* i  : number of values to be inserted */
)
{
    float cv_out[LATTICE_DIM];
    int i, occ[LATTICE_DIM], cnt;

    for( i=0; i<dim+no_new_val; i++ )
    {
        occ[i] = 0;
    }

    for( i=0; i<no_new_val; i++ )
    {
        cv_out[p[i]] = val;
        occ[p[i]]  = 1;
    }

    cnt = 0;
    for( i=0; i<dim+no_new_val; i++ )
    {
        if (occ[i] == 0)
        {
            cv_out[i] = cv[cnt++];
        }
    }

    for( i=0; i<dim+no_new_val; i++ )
    {
        cv[i] = cv_out[i];
    }

    return;
}

/*-----------------------------------------------------------------*
 * idx2c()
 *
 * decode index of binomial combinations, find the positions of k components out of n total components
 *-----------------------------------------------------------------*/

static void idx2c(
    int n,                     /* i  : total number of positions (components)  */
    int *p,                    /* o  : array with positions of the k components */
    int k,                     /* i  : number of components whose position is to be determined */
    int val                    /* i  : index to be decoded */
)
{
    int i, skip, pos, k1;

    skip = 0;
    pos = 0;
    k1 = k-1;
    while( skip+C[n-pos-1][k1] -1 < val )
    {
        skip += C[n-pos-1][k1];
        pos++;
    }

    p[0] = pos;
    n -= pos+1;
    val -= skip;
    if ( k == 1 )
    {
        return;
    }

    idx2c( n, p+1, k1, val );

    /* pos+1 */
    for( i=1; i<k; i++ )
    {
        p[i] += pos+1;
    }

    return;
}

/*-----------------------------------------------------------------*
 * decode_sign_pc1()
 *
 *-----------------------------------------------------------------*/

void decode_sign_pc1(
    float *c,       /* o  : decoded codevector */
    int idx_sign,   /* i  : sign index */
    int parity      /* i  : parity flag (+1/-1/0) */
)
{
    int i, len = LATTICE_DIM, cnt_neg = 1;

    if ( parity )
    {
        len -= 1;
    }

    for( i=0; i<len; i++ )
    {
        if (c[i] > 0)
        {
            if (idx_sign % 2)
            {
                c[i] = -c[i];
                cnt_neg = -cnt_neg;
            }
            idx_sign >>= 1;
        }
    }

    if ( len < LATTICE_DIM )
    {
        if (cnt_neg != parity)
        {
            c[len] = -c[len];
        }
    }

    return;
}

/*-----------------------------------------------------------------*
 * extract_low()
 *
 * (!!!!! function for int64 !!!!)
 *-----------------------------------------------------------------*/

static unsigned int extract_low(
    unsigned int x
)
{
    return (x&(0xffff));
}

/*-----------------------------------------------------------------*
 * extract_high()
 *
 * (!!!!! function for int64 !!!!)
 *-----------------------------------------------------------------*/

static unsigned int extract_high(
    unsigned int x
)
{
    return (x>>16);
}

/*-----------------------------------------------------------------*
 * multiply32_32_64()
 *
 * (!!!!! function for int64 !!!!)
 *-----------------------------------------------------------------*/

void multiply32_32_64(
    unsigned int x,
    unsigned int y,
    unsigned int *res
)
{
    unsigned int tmp, x_tmp[2], y_tmp[2];
    unsigned int high = 0;

    x_tmp[0] = extract_low(x); /* lowest 16 bits */
    x_tmp[1] = extract_high(x);
    y_tmp[0] = extract_low(y);
    y_tmp[1] = extract_high(y);
    tmp = x_tmp[0]*y_tmp[0];
    high = extract_high(tmp);
    res[0] = extract_low(tmp);
    tmp = x_tmp[1]*y_tmp[0]+ x_tmp[0]*y_tmp[1] + high; /* x and y are not using all 32 bits */
    high = extract_high(tmp);
    res[0] += (extract_low(tmp)<<16);
    tmp = x_tmp[1]*y_tmp[1]+high;
    res[1] = tmp;

    return;
}

/*-----------------------------------------------------------------*
 * get_no_bits()
 *
 * (!!!!! function for int64 !!!!)
 *-----------------------------------------------------------------*/

static int get_no_bits(
    unsigned int x
)
{
    int nb = 0;

    while( x > 0 )
    {
        x >>= 1;
        nb++;
    }

    return nb;
}

/*-----------------------------------------------------------------*
 * divide_64_32()
 *
 * (!!!!! function for int64 !!!!)
 *-----------------------------------------------------------------*/

static void divide_64_32(
    short *xs,              /* i  : denominator as array of two int32 */
    unsigned int y,         /* i  : nominator on 32 bits */
    unsigned int *result,   /* o  : integer division result on 32 bits */
    unsigned int *rem       /* o  : integer division reminder on 32 bits */
)
{
    int nb_x1;
    unsigned int r, q, q1, x_tmp, x[2];

    x[0] = (((unsigned int)xs[2]&((1<<2)-1))<<(LEN_INDICE*2)) + (xs[1]<<LEN_INDICE) + xs[0];
    x[1] = xs[2]>>2;

    /* find number of bits of x[0] and x[1] */
    nb_x1 = get_no_bits(x[1]);

    /* take the first 32 bits */
    if ( nb_x1 > 0 )
    {
        x_tmp = (x[1]<<(32-nb_x1)) + (x[0]>>nb_x1);
        q = (unsigned int)(x_tmp/y+0.5);
        r = x_tmp-q*y; /* this is the first reminder */
        r = (r<<nb_x1)+(x[0]&((1<<nb_x1) - 1));


        q1 = (unsigned int)(r/y +0.5);
        *result = (q<<nb_x1) + q1;
        *rem = r - q1*y;
    }
    else
    {
        x_tmp = x[0];
        q = ((unsigned int)(x_tmp/y+0.5));
        *result = q;
        *rem = x_tmp-q*y;
    }

    return;
}
