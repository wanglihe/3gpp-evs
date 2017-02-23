/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "prot.h"
#include "rom_com.h"
#include "cnst.h"

/*-----------------------------------------------------------------*
 * Local functions
 *-----------------------------------------------------------------*/

static float quantize_data( float *data, const float *w_in, float *qin, float *cv_out, int *idx_lead, int *idx_scale, const float *sigma,
                            const float * inv_sigma,
                            const float *scales, short no_scales, const int *no_lead );
static float q_data( float *pTmp1, const float *w1, float *quant, float *cv_out, int *idx_lead, int *idx_scale,
                     const float * p_inv_sigma,
                     const float *p_sigma, const float *p_scales, short *p_no_scales, const int *p_no_lead );
static void prepare_data( float *xsort, int *sign,
                          float *data, float *w, const float *w_in,
                          const float * sigma,
                          const float * inv_sigma,
                          int *p_sig );
static float calculate_min_dist( float cv_pot[LATTICE_DIM], int no_scales, const float *scale, const float *w, int *p_best_scale,
                                 int *p_best_idx, const int *no_leaders, int sig, int * indx);
static int find_pos( float *c, int len,    float arg, int *p );
static void take_out_val( float *v, float *v_out, float val, int len );
static int index_leaders( float *cv, int idx_lead, int dim );
static int c2idx( int n, int *p, int k );
static int encode_sign_pc1( int parity, float *cv);
static int encode_comb( float *cv, int idx_lead );
static void sort_desc_ind( float *s, int len, int *ind );


/*-----------------------------------------------------------------*
 * mslvq()
 *
 * Encodes the LSF residual
 *-----------------------------------------------------------------*/
float mslvq (
    float *pTmp,               /* i  : M-dimensional input vector */
    float *quant,              /* o  : quantized vector */
    float *cv_out,             /* o  : corresponding 8-dim lattice codevectors (without the scaling) */
    int   *idx_lead,           /* o  : leader index for each 8-dim subvector  */
    int   *idx_scale,          /* o  : scale index for each subvector */
    float *w,                  /* i  : weights for LSF quantization */
    short mode,                /* i  : number indicating the coding type (V/UV/G...)*/
    short mode_glb,            /* i  : LVQ coding mode */
    int   pred_flag,            /* i  : prediction flag (0: safety net, 1,2 - predictive  )*/
    short  no_scales[][2]
)
{
    float dist;
    const float * p_scales, *p_sigma
    ,*p_inv_sigma
    ;
    const int  *p_no_lead;
    short * p_no_scales;

    dist = 0.0f;
    p_no_scales = no_scales[mode_glb];
    if ( pred_flag == 0 )
    {
        p_sigma = sigma[mode];
        /* inverse sigma is precomputed to save complexity */
        p_inv_sigma = inv_sigma[mode];
        p_scales = scales[mode_glb];
        p_no_lead = no_lead[mode_glb];

    }
    else
    {
        p_sigma = sigma_p[mode];
        /* inverse sigma is precomputed to save complexity */
        p_inv_sigma = inv_sigma_p[mode];
        p_scales = scales_p[mode_glb];
        p_no_lead = no_lead_p[mode_glb];
    }
    /* first subvector */
    dist += quantize_data( pTmp, w, quant, cv_out, idx_lead, idx_scale,
                           p_sigma,
                           p_inv_sigma,
                           p_scales, p_no_scales[0], p_no_lead );
    /* second subvector */
    dist += quantize_data( pTmp+LATTICE_DIM, w+LATTICE_DIM, quant+LATTICE_DIM,
                           cv_out+LATTICE_DIM, &idx_lead[1], &idx_scale[1],
                           p_sigma+LATTICE_DIM,
                           p_inv_sigma+LATTICE_DIM,
                           p_scales+MAX_NO_SCALES,
                           p_no_scales[1], p_no_lead+MAX_NO_SCALES );

    return dist;
}
/*-----------------------------------------------------------------*
 * q_data()
 *
 * (used for LSF quantization in CNG)
 *-----------------------------------------------------------------*/

static float q_data (
    float *pTmp1,
    const float *w1,
    float *quant,
    float *cv_out,
    int   *idx_lead,
    int   *idx_scale,
    const float *p_sigma,
    const float *p_inv_sigma,
    const float *p_scales,
    short *p_no_scales,
    const int *p_no_lead
)
{
    float dist = 0.0f;
    /* first subvector */
    dist += quantize_data( pTmp1, w1, quant, cv_out, idx_lead, idx_scale, p_sigma,
                           p_inv_sigma, p_scales, p_no_scales[0], p_no_lead );
    /* second subvector */
    dist += quantize_data( pTmp1+LATTICE_DIM, w1+LATTICE_DIM, quant+LATTICE_DIM, cv_out+LATTICE_DIM, &idx_lead[1], &idx_scale[1],
                           p_sigma+LATTICE_DIM,p_inv_sigma+LATTICE_DIM,
                           p_scales+MAX_NO_SCALES, p_no_scales[1], p_no_lead+MAX_NO_SCALES );

    return dist;
}

/*-----------------------------------------------------------------*
 * mslvq_cng()
 *
 * Encodes the LSF residual in SID frames
 *-----------------------------------------------------------------*/

float mslvq_cng (
    short idx_cv,              /* i  : index of cv from previous stage */
    float *pTmp,               /* i  : 16 dimensional input vector */
    float *quant,              /* o  : quantized vector */
    float *cv_out,             /* o  : corresponding 8-dim lattice codevectors (without the scaling) */
    int   *idx_lead,           /* o  : leader index for each 8-dim subvector  */
    int   *idx_scale,          /* o  : scale index for each subvector */
    const float *w,            /* i  : weights for LSF quantization */
    short * no_scales
)
{
    float dist;
    const float *p_scales, *p_sigma
    , *p_inv_sigma
    ;
    const int *p_no_lead;
    short *p_no_scales;
    short mode_glb, mode, i;
    float pTmp1[M], w1[M];

    dist = 0.0f;
    mode = (short)LVQ_COD_MODES + idx_cv;

    /* for CNG there is only one bitrate but several quantizer structures, depending on the previous VQ stage */
    mode_glb = START_CNG + idx_cv;

    p_sigma = sigma[mode];
    p_inv_sigma = inv_sigma[mode];
    p_scales = scales[mode_glb];
    p_no_lead = no_lead[mode_glb];
    p_no_scales = &no_scales[mode_glb*2];
    /* check if LSF component permutation is needed or not */
    if ( cng_sort[idx_cv] )
    {
        /* change order in subvecs */
        for( i=0; i<M; i++ )
        {
            pTmp1[i] = pTmp[i];
            w1[i] = w[i];
        }
        /* sorting the quantizer input and the corresponding weights according to the specified permutations */
        permute(pTmp1, perm[idx_cv]);
        permute(w1, perm[idx_cv]);

        dist = q_data( pTmp1, w1, quant, cv_out, idx_lead, idx_scale, p_sigma,
                       p_inv_sigma, p_scales, p_no_scales, p_no_lead );
        /* permute back */
        permute( quant, perm[idx_cv] );
    }
    else
    {
        dist = q_data( pTmp, w, quant, cv_out, idx_lead, idx_scale, p_sigma,
                       p_inv_sigma,
                       p_scales, p_no_scales, p_no_lead );
    }

    return dist;
}
/*-----------------------------------------------------------------*
 * prepare_data()
 *
 *-----------------------------------------------------------------*/

static void prepare_data (
    float *xsort,
    int   *sign,
    float *data,
    float *w,
    const float *w_in,
    const float *sigma,
    const float * inv_sigma,
    int   *p_sig
)
{
    int   j, sig;
    float s
    ,inv_s
    ;

    /* scale data */
    for( j=0; j<LATTICE_DIM; j++ )
    {
        s = sigma[j];
        inv_s = inv_sigma[j];
        xsort[j] = data[j]*inv_s;
        w[j] = w_in[j]*(s*s);
    }

    sig = 1;
    for( j=0; j<LATTICE_DIM; j++ )
    {
        if ( xsort[j] < 0 )
        {
            sign[j] = -1;
            sig = -sig;
            xsort[j] = -xsort[j];
        }
        else
        {
            sign[j] = 1;
        }
    }
    *p_sig = sig;

    return;
}

/*-----------------------------------------------------------------*
 * calculate_min_dist()
 *
 *-----------------------------------------------------------------*/
static float calculate_min_dist(
    float cv_pot[LATTICE_DIM],
    int   no_scales,
    const float *scale,
    const float *w,
    int   *p_best_scale,
    int   *p_best_idx,
    const int *no_leaders,
    int sig,
    int * indx
)
{
    int k, l,j, best_scale = -1, best_idx = -1;
    float s, s2, tmp_dist, min_dist, wx[LATTICE_DIM], wind[LATTICE_DIM];
    float sum1[NO_LEADERS], sum2[NO_LEADERS];
    const float *pl_crt;
    float p;

    /* compare first with the origin */
    min_dist = 0.0f;
    for(j=0; j<LATTICE_DIM; j++)
    {
        /* sorting the weight based on the ordering indx[] of the input vector */
        wind[j] = w[indx[j]];
        wx[j] = 2.0f*wind[j]*cv_pot[j];
    }
    s = scale[0];
    s2 =  s*s;
    pl_crt = &pl[0];

    for(j=0; j<no_leaders[0]; j++)
    {
        sum1[j] = 0;
        sum2[j] = 0;
        l = 0;
        while(l<LATTICE_DIM-1)
        {
            p = *pl_crt;
            if (p)
            {
                sum1[j] += wx[l]*p;
                sum2[j] += wind[l]*p*p;
                pl_crt++;
                l++;
            }
            else
            {
                pl_crt += LATTICE_DIM-l;
                l = LATTICE_DIM;
            }
        }
        if ((l-LATTICE_DIM+1)==0)
        {
            p = *pl_crt;
            /* if it went up to 7th position */
            if (pl_par[j])
            {
                if (sig != pl_par[j])
                {
                    sum1[j] -= wx[l]*p;
                    sum2[j] += wind[l]*p*p;
                    pl_crt++;
                }
                else
                {
                    sum1[j] += wx[l]*p;
                    sum2[j] += wind[l]*p*p;
                    pl_crt++;
                }
            }
            else
            {
                sum1[j] += wx[l]*p;
                sum2[j] += wind[l]*p*p;
                pl_crt++;
            }
        }
        /* distance between the potential codevector and the input calculated in ordered space */
        tmp_dist = s2*sum2[j] - s*sum1[j];
        if (tmp_dist <min_dist)
        {
            min_dist = tmp_dist;
            best_scale = 0;
            best_idx = j;
        }
    }
    tmp_dist = min_dist + 1.0f;
    for (k=1; k<no_scales; k++)
    {
        s = scale[k];
        s2 = s*s;
        for (j=0; j<no_leaders[k]; j++)
        {
            /* distance between the potential codevector and the input calculated in ordered space */
            tmp_dist = s2*sum2[j] - s*sum1[j];
            if (tmp_dist <min_dist)
            {
                min_dist = tmp_dist;
                best_scale = k;
                best_idx = j;
            }
        }
    }
    *p_best_scale = best_scale;
    *p_best_idx = best_idx;
    return min_dist;
}
/*-----------------------------------------------------------------*
 * quantize_data()
 *
 *-----------------------------------------------------------------*/

static float quantize_data(
    float *data,             /* i  : residual LSF data to quantize */
    const float *w_in,         /* i  : weights                       */
    float *qin,              /* o  : quantized output (scaled)     */
    float *cv_out,           /* o  : codevectors   */
    int   *idx_lead,         /* o  : leader indexes for each subvector */
    int   *idx_scale,        /* o  : scale indexes for each subvector */
    const float *sigma,      /* i  : standard deviation             */
    const float * inv_sigma, /* i  : inverse of standard deviation   */
    const float *scale,      /* i  : scales for each truncation     */
    short no_scales,         /* i  : number of truncation for each subvector */
    const int *no_leaders    /* i  : number of leader vectors for each truncation of each subvector */
)
{
    int j
    ;
    float
    w[LATTICE_DIM],  min_dist = 0;
    int best_idx = 0, best_scale = -1;
    float s;
    float cv_pot[LATTICE_DIM];
    int indx[LATTICE_DIM];
    int  sig, sign[LATTICE_DIM]
    ;
    int smallest;
    int id[LATTICE_DIM];
    if (no_scales>0)
    {
        prepare_data( cv_pot, sign, data, w, w_in, sigma, inv_sigma, &sig );
        /* sorting of the input vector based on its absolute values; indx: permutation corresponding to the sorting */
        sort_desc_ind(cv_pot, LATTICE_DIM, indx);
        smallest = indx[LATTICE_DIM-1];
        min_dist = calculate_min_dist(cv_pot, no_scales, scale, w, &best_scale, &best_idx, no_leaders, sig, indx);
        if (best_scale > -1)
        {
            for(j=0; j<LATTICE_DIM; j++)
            {
                id[indx[j]] = j;
            }
            for(j=0; j<LATTICE_DIM; j++)
            {
                cv_out[j] = sign[j]*pl[best_idx*LATTICE_DIM+id[j]];
            }
            if (pl_par[best_idx])
            {
                if ( sig - pl_par[best_idx] != 0 )
                {
                    cv_out[smallest] = -cv_out[smallest];
                }
            }
            s = scale[best_scale];
            for( j=0; j<LATTICE_DIM; j++ )
            {
                qin[j] = cv_out[j]*s*sigma[j];
            }
            *idx_lead = best_idx;
            *idx_scale = best_scale;
        }
        else
        {
            for( j=0; j<LATTICE_DIM; j++ )
            {
                qin[j] = 0;
            }

            *idx_lead = best_idx;
            *idx_scale = best_scale;
        }
    }
    else
    {
        *idx_lead = 0;
        *idx_scale = -1;
        for( j=0; j<LATTICE_DIM; j++ )
        {
            cv_out[j] = 0;
            qin[j] = 0;
        }
    }

    return min_dist;
}

/*-----------------------------------------------------------------*
 * sort_desc_ind()
 *
 * sorts in descending order and computes indices in the sorted vector
 *-----------------------------------------------------------------*/

static void sort_desc_ind(
    float *s,         /* i/o: vector to be sorted */
    int   len,        /* i  : vector length */
    int   *ind        /* o  : array of indices */
)
{
    int i, k, sorted, a;
    float t;

    for ( i=0 ; i<len ; i++ )
    {
        ind[i] = i;
    }
    sorted = 0;
    for ( k=len-1 ; k && !sorted ; k-- )
    {
        sorted = 1;
        for ( i=0 ; i < k ; i++ )
        {
            if ( s[i] < s[i+1] )
            {
                sorted = 0;
                t = s[i];
                s[i] = s[i+1];
                s[i+1] = t;
                a = ind[i];
                ind[i] = ind[i+1];
                ind[i+1] = a;
            }
        }
    }

    return;
}


/*-----------------------------------------------------------------*
 * index_lvq()
 *
 * sorts in descending order and computes indices in the sorted vector
 *-----------------------------------------------------------------*/

void index_lvq (
    float *quant,           /* i  : codevector to be indexed (2 8-dim subvectors)*/
    int   *idx_lead,        /* i  : leader class index for each subvector */
    int   *idx_scale,       /* i  : scale index for each subvector */
    int   mode,             /* i  : integer signalling the quantizer structure for the current bitrate */
    short *index,           /* o  : encoded index (represented on 3 short each with 15 bits ) */
    unsigned int * p_offset_scale1,
    unsigned int * p_offset_scale2,
    short * p_no_scales
)
{
    unsigned int index1, index2, tmp, idx[2];
    int len_offset;

    len_offset = MAX_NO_SCALES+1;

    index1 = 0;

    /* for first subvector */
    if ( idx_scale[0] > -1 )
    {
        index1 = encode_comb(quant, idx_lead[0]) + table_no_cv[idx_lead[0]] + p_offset_scale1[mode*len_offset +idx_scale[0]];
    }

    /* for second subvector */
    index2 = 0;
    if ( idx_scale[1] > -1 )
    {
        index2 = encode_comb(&quant[LATTICE_DIM], idx_lead[1]) + table_no_cv[idx_lead[1]] + p_offset_scale2[mode*len_offset+idx_scale[1]];
    }

    multiply32_32_64(index1, (unsigned int)(p_offset_scale2[mode*len_offset+p_no_scales[mode*2+1]]), idx);

    tmp = idx[0] + index2;
    if ( tmp < idx[0] || tmp < index2 )
    {
        idx[1] +=1;
    }

    idx[0] = tmp;

    /* convert to 3 short */
    index[0] = ((idx[0])&(0xffff>>1));
    index[1] = ((idx[0])>>15)&(0xffff>>1);
    index[2] = ((idx[0])>>30) + (((idx[1])<<2)&(0xffff>>1));

    return;
}



/*-----------------------------------------------------------------*
 * encode_comb()
 *
 * creates an index for the lattice codevector
 *-----------------------------------------------------------------*/

static int encode_comb(        /* o  : index of the absolute valued codevector */
    float *cv,                 /* i  : codevector to be indexed */
    int   idx_lead             /* i  : leader class index, to know the values */
)
{
    int idx_sign, idx_ld_class;

    idx_sign = encode_sign_pc1( pl_par[idx_lead], cv );
    idx_ld_class = index_leaders( cv, idx_lead,  LATTICE_DIM );

    return idx_sign * pi0[idx_lead] + idx_ld_class;
}


/*-----------------------------------------------------------------*
 * index_leaders()
 *
 * gives the index in a class of leaders without considering the sign yet
 *-----------------------------------------------------------------*/

static int index_leaders(     /* o  : index */
    float *cv,                /* i  : codevector to be indexed */
    int   idx_lead,           /* i  : leader class index */
    int   dim                 /* i  : vector dimension */
)
{
    int index, i, no_vals_loc, nr, p[LATTICE_DIM], dim_loc;
    float cv_copy[LATTICE_DIM], val_crt;

    no_vals_loc = no_vals[idx_lead];

    if ( no_vals_loc == 1 )
    {
        return 0;
    }

    for( i=0; i<LATTICE_DIM; i++ )
    {
        cv_copy[i] = (float)fabs(cv[i]);
    }

    val_crt = vals[idx_lead][0];
    nr = find_pos(cv_copy, dim, val_crt, p);
    index = c2idx(LATTICE_DIM, p, nr);

    if ( no_vals_loc == 2 )
    {
        return index;
    }

    take_out_val( cv_copy, cv_copy, val_crt, dim );
    dim_loc = dim-no_vals_ind[idx_lead][0];
    index *= C[dim_loc][no_vals_ind[idx_lead][1]];
    val_crt = vals[idx_lead][1];
    nr = find_pos( cv_copy, dim_loc, val_crt, p );
    index += c2idx( dim_loc, p, nr );

    if ( no_vals_loc == 3 )
    {
        return index;
    }

    take_out_val(cv_copy, cv_copy, val_crt, dim_loc);
    dim_loc = dim_loc-no_vals_ind[idx_lead][1];
    index *= C[dim_loc][no_vals_ind[idx_lead][2]];
    val_crt = vals[idx_lead][2];
    nr = find_pos(cv_copy, dim_loc, val_crt, p);
    index += c2idx(dim_loc, p, nr);
    /* maximum 4 values */

    return index;
}

/*-----------------------------------------------------------------*
 * find_pos()
 *
 * Finds the positions in vector c for which the vector components are equal to 'arg'.
 * It returns the number of such positions and their values in the array 'p'.
 *-----------------------------------------------------------------*/

int find_pos(         /* o  : number of positions */
    float *c,         /* i  : input vector */
    int   len,        /* i  : input vector dim */
    float arg,        /* i  : argument to be compared with */
    int   *p          /* o  : vector of positions */
)
{
    int i=0, j=0;

    /* how many (j) and which (p) positions are in the relation pred(arg,c[i]) */
    for( i=0; i<len; i++)
    {
        if( (int)(arg*10) == (int)(c[i]*10))
        {
            p[j++]=i;
        }
    }

    return j;
}

/*-----------------------------------------------------------------*
 * encode_sign_pc1()
 *
 * Creates an index for signs of the significant codevector components
 * Gives the index of the signs - binary representation where negative sign stands for 1
 * and positive sign stands for 1.
 *-----------------------------------------------------------------*/

static int encode_sign_pc1(    /* o  : index of signs */
    int   parity,              /* i  : parity of the leader class to which the codevector belongs */
    float *cv                  /* i  : input codevector */
)
{
    int idx_sign, cnt, i, len=LATTICE_DIM;

    idx_sign = 0;
    cnt = 0;

    if ( parity )
    {
        len -= 1;
    }

    for( i=0; i<len; i++ )
    {
        if ( cv[i] < 0 )
        {
            idx_sign += (1<<cnt);
            cnt++;
        }

        if ( cv[i] > 0 )
        {
            cnt ++;
        }
    }

    return idx_sign;
}

/*-----------------------------------------------------------------*
 * take_out_val()
 *
 * removes the value val from the vector v
 *-----------------------------------------------------------------*/

static void take_out_val(
    float *v,            /* i  : input vector */
    float *v_out,        /* o  : output vector without the value val*/
    float val,           /* i  : value to be removed */
    int   len            /* i  : input vector length */
)
{
    int i, cnt;

    cnt = 0;

    for( i=0; i<len; i++ )
    {
        if ((int)(v[i]*10)!= (int)(val*10))
        {
            v_out[cnt++] = v[i];
        }
    }

    return;
}


/*-----------------------------------------------------------------*
 * c2idx()
 *
 *-----------------------------------------------------------------*/

int c2idx(
    int n,
    int *p,
    int k
)
{
    int skip, i, p0;

    if ( k == 1 )
    {
        return p[0];
    }
    else
    {
        skip = 0;
        for( i=1; i<=p[0]; i++ )
        {
            skip += C[n-i][k-1];
        }

        p0 = p[0];
        for( i=1; i<k; i++)
        {
            p[i] -= p0+1;
        }

        return skip + c2idx( n-p0-1, p+1, k-1 );
    }
}
