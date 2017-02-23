/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"


/*-------------------------------------------------------------------*
 * Local functions
 *-------------------------------------------------------------------*/

static int re8_identify_absolute_leader(int y[]);
static void re8_coord(int y[], int k[]);

/*----------------------------------------------------------------*
 * re8_vor()
 *
 * Multi-rate RE8 indexing by Voronoi extension
 *----------------------------------------------------------------*/

void re8_vor(
    int y[],  /* i  : point in RE8 (8-dimensional integer vector)      */
    int *n,   /* o  : codebook number n=0,2,3,4,... (scalar integer)   */
    int k[],  /* o  : Voronoi index (integer vector of dimension 8) used only if n>4 */
    int c[],  /* o  : codevector in Q0, Q2, Q3, or Q4 if n<=4, y=c     */
    int *ka   /* o  : identifier of absolute leader (needed to index c)*/
)
{
    int   i, r, m, v[8], c_tmp[8], k_mod[8], k_tmp[8], iter, ka_tmp, n_tmp, mask;
    float sphere;


    /*----------------------------------------------------------------*
     * verify if y is in Q0, Q2, Q3 or Q4
     *   (a fast search is used here:
     *    the codebooks Q0, Q2, Q3 or Q4 are specified in terms of RE8 absolute leaders
     *    (see FORinstance Xie and Adoul's paper in ICASSP 96)
     *    - a unique code identifying the absolute leader related to y is computed
     *      in re8_identify_absolute_leader()
     *      this code is searched FORin a pre-defined list which specifies Q0, Q2, Q3 or Q4)
     *      the absolute leader is identified by ka
     *    - a translation table maps ka to the codebook number n)
     *----------------------------------------------------------------*/

    *ka = re8_identify_absolute_leader( y );

    /*----------------------------------------------------------------*
     * compute codebook number n of Qn (by table look-up)
     *   at this stage, n=0,2,3,4 or out=100
     *----------------------------------------------------------------*/

    *n = Da_nq[*ka];

    /*----------------------------------------------------------------*
     * decompose y into :
     *     (if n<=4:)
     *     y = c        where c is in Q0, Q2, Q3 or Q4
     *   or
     *     (if n>4:)
     *     y = m c + v  where c is in Q3 or Q4, v is a Voronoi codevector
     *                        m=2^r (r integer >=2)
     *
     *   in the latter case (if n>4), as a side-product, compute the (Voronoi) index k[] of v
     *   and replace n by n = n' + 2r where n' = 3 or 4 (c is in Qn') and r is defined above
     *----------------------------------------------------------------*/

    if( *n <= 4 )
    {
        for( i=0; i<8; i++ )
        {
            c[i] = y[i];
        }
    }
    else
    {
        /*------------------------------------------------------------*
         * initialize r and m=2^r based on || y ||^2/8
         *------------------------------------------------------------*/

        sphere = 0.0;
        for( i=0; i<8; i++ )
        {
            sphere += (float)y[i]*(float)y[i];
        }
        sphere *= 0.125;
        r = 1;
        sphere *= 0.25;                         /* not counted, merge 0.125*0.25 */
        while( sphere > 11.0 )
        {
            r++;
            sphere *= 0.25;
        }
        /*------------------------------------------------------------*
         * compute the coordinates of y in the RE8 basis
         *------------------------------------------------------------*/

        re8_coord( y, k_mod );

        /*------------------------------------------------------------*
         * compute m and the mask needed for modulo m (for Voronoi coding)
         *------------------------------------------------------------*/

        m = 1<<r;
        mask = m-1; /* 0x0..011...1 */

        /*------------------------------------------------------------*
         * find the minimal value of r (or equivalently of m) in 2 iterations
         *------------------------------------------------------------*/

        for( iter=0; iter<2; iter++ )
        {
            /*--------------------------------------------------------*
             * compute v such that y is in m RE_8 +v (by Voronoi coding)
             *--------------------------------------------------------*/

            for( i=0 ; i<8; i++ )
            {
                k_tmp[i] = k_mod[i] & mask;
            }

            re8_k2y( k_tmp, m, v );

            /*--------------------------------------------------------*
             * compute c = (y-v)/m
             * (y is in RE8, c is also in RE8 by definition of v)
             *--------------------------------------------------------*/

            for( i=0; i<8; i++ )
            {
                c_tmp[i] = (y[i]-v[i]) / m;
                /* M IS A FACTOR OF 2 */
            }

            /*--------------------------------------------------------*
             *  verify if c_tmp is in Q2, Q3 or Q4
             *--------------------------------------------------------*/

            ka_tmp = re8_identify_absolute_leader( c_tmp );

            /*--------------------------------------------------------*
             * at this stage, n_tmp=2,3,4 or out = 100 -- n=0 is not possible
             *--------------------------------------------------------*/

            n_tmp = Da_nq[ka_tmp];

            if( n_tmp > 4 )
            {
                /*--------------------------------------------------------*
                 * if c is not in Q2, Q3, or Q4 (i.e. n_tmp>4), use m = 2^(r+1) instead of 2^r
                 *--------------------------------------------------------*/

                r++;
                m = m<<1;
                mask = ((mask<<1)+1);            /* mask = m-1; <- this is less complex */
            }
            else
            {
                /*--------------------------------------------------------*
                 * c is in Q2, Q3, or Q4 -> the decomposition of y as y = m c + v is valid
                 *
                 * since Q2 is a subset of Q3, indicate n=3 instead of n=2 (this is because
                 * for n>4, n=n'+2r with n'=3 or 4, so n'=2 is not valid)
                 *--------------------------------------------------------*/

                if( n_tmp < 3 )
                {
                    n_tmp = 3;
                }

                /*--------------------------------------------------------*
                 * save current values into ka, n, k and c
                 *--------------------------------------------------------*/

                *ka = ka_tmp;
                *n = n_tmp + 2*r;
                for (i=0; i<8; i++)
                {
                    k[i] = k_tmp[i];
                    c[i] = c_tmp[i];
                }

                /*--------------------------------------------------------*
                 * try  m = 2^(r-1) instead of 2^r to be sure that m is minimal
                 *--------------------------------------------------------*/

                r--;
                m = m>>1;
                mask = mask>>1;
            }
        }
    }

    return;
}

/*-----------------------------------------------------------------------*
 * re8_identify_absolute_leader()
 *
 * Identify the absolute leader related to y using a pre-defined table which
 * specifies the codebooks Q0, Q2, Q3 and Q4
 -----------------------------------------------------------------------*/

static int re8_identify_absolute_leader(    /* o : integer indicating if y if in Q0, Q2, Q3 or Q4 (or if y is an outlier)   */
    int y[]                                 /* i : point in RE8 (8-dimensional integer vector)                              */
)
{
    int i,s,id,C[8],nb,pos,ka;
    long tmp;


    /*-----------------------------------------------------------------------*
     * compute the RE8 shell number s = (y1^2+...+y8^2)/8 and C=(y1^2, ..., y8^2)
     *-----------------------------------------------------------------------*/

    s = 0;
    for( i=0; i<8; i++ )
    {
        C[i] = y[i]*y[i];
        s += C[i];
    }
    s >>= 3;

    /*-----------------------------------------------------------------------*
     * compute the index 0 <= ka <= NB_LEADER+1 which identifies an absolute leader of Q0, Q2, Q3 or Q4
     *
     * by default, ka=index of last element of the table (to indicate an outlier)
     *-----------------------------------------------------------------------*/

    ka = NB_LEADER+1;
    if( s == 0 )
    {
        /*-------------------------------------------------------------------*
         * if s=0, y=0 i.e. y is in Q0 -> ka=index of element indicating Q0
         *-------------------------------------------------------------------*/

        ka = NB_LEADER;
    }
    else
    {
        /*-------------------------------------------------------------------*
         * the maximal value of s for y in  Q0, Q2, Q3 or Q4 is NB_SPHERE
         *   if s> NB_SPHERE, y is an outlier (the value of ka is set correctly)
         *-------------------------------------------------------------------*/

        if( s <= NB_SPHERE )
        {
            /*---------------------------------------------------------------*
             * compute the unique identifier id of the absolute leader related to y:
             * s = (y1^4 + ... + y8^4)/8
             *---------------------------------------------------------------*/

            tmp = 0;
            for( i=0; i<8; i++ )
            {
                tmp += C[i]*C[i];
            }
            id = (int)(tmp >> 3);

            /*---------------------------------------------------------------*
             * search for id in table Da_id
             * (containing all possible values of id if y is in Q2, Q3 or Q4)
             * this search is focused based on the shell number s so that
             * only the id's related to the shell of number s are checked
             *---------------------------------------------------------------*/

            nb = Da_nb[s-1];        /* get the number of absolute leaders used on the shell of number s */
            pos = Da_pos[s-1];      /* get the position of the first absolute leader of shell s in Da_id */
            for( i=0; i<nb; i++ )
            {
                if( id == Da_id[pos] )
                {
                    ka = pos; /* get ka */
                    break;
                }
                pos++;
            }
        }
    }

    return( ka );
}

/*-------------------------------------------------------------------------
 * re8_coord()
 *
 * COMPUTATION OF RE8 COORDINATES
 -----------------------------------------------------------------------*/

static void re8_coord(
    int *y,  /* i  : 8-dimensional point y[0..7] in RE8 */
    int *k   /* o  : coordinates k[0..7]                */
)
{
    int i, tmp, sum;


    /*---------------------------------------------------------------*
     * compute k = y M^-1
     *   M = 1/4 [ 1          ]
     *           [-1  2       ]
     *           [ |    \     ]
     *           [-1       2  ]
     *           [ 5 -2 _ -2 4]
     *
     *---------------------------------------------------------------*/

    k[7]=y[7];
    tmp = y[7];
    sum = 5*y[7];
    for( i=6; i>=1; i-- )
    {
        /* apply factor 2/4 from M^-1 */
        k[i] = (y[i]-tmp)>>1;
        sum -= y[i];
    }

    /* apply factor 1/4 from M^-1 */
    k[0] = (y[0]+sum)>>2;

    return;
}


/*-------------------------------------------------------------------------
 * re8_k2y()
 *
 * Voronoi indexing (index decoding) k -> y
 -----------------------------------------------------------------------*/

void re8_k2y(
    const int *k,  /* i  : Voronoi index k[0..7]                                   */
    const int m,   /* i  : Voronoi modulo (m = 2^r = 1<<r, where r is integer >=2) */
    int *y   /* o  : 8-dimensional point y[0..7] in RE8                      */
)
{
    int i, v[8], tmp, sum, *ptr1, *ptr2;
    float z[8];


    /*---------------------------------------------------------------*
     * compute y = k M and z=(y-a)/m, where
     *   M = [4        ]
     *       [2 2      ]
     *       [|   \    ]
     *       [2     2  ]
     *       [1 1 _ 1 1]
     *   a=(2,0,...,0)
     *---------------------------------------------------------------*/

    for (i=0; i<8; i++)
    {
        y[i] = k[7];
    }

    z[7] = (float)y[7]/m;
    sum = 0;
    for( i=6; i>=1; i-- )
    {
        tmp   = 2*k[i];
        sum  += tmp;
        y[i] += tmp;
        z[i] = (float)y[i]/m;   /* m is a factor of 2*/
    }
    y[0] += (4*k[0] + sum);
    z[0] = (float)(y[0]-2)/m;

    /*---------------------------------------------------------------*
     * find nearest neighbor v of z in infinite RE8
     *---------------------------------------------------------------*/

    re8_PPV( z, v );

    /*---------------------------------------------------------------*
     * compute y -= m v
     *---------------------------------------------------------------*/

    ptr1 = y;
    ptr2 = v;
    for( i=0; i<8; i++ )
    {
        *ptr1++ -= m **ptr2++;
    }

    return;
}
