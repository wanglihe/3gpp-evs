/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"


/*-------------------------------------------------------------------*
 *  LOCAL DEFINITIONS
 *-------------------------------------------------------------------*/

#define N_OPT           5           /* complexity setting,  direct functional calculation  limit              */
#define TABLE_LIM_OPT   96          /* odd divisor table recursion limit setting,  due to dim 6               */

/* local typedefs for optimized pvq indexing, used only in this c-file to vectorize common function calls for c-code clarity */
typedef  void          (*VEC2INDFUNCM) (const short* ,  short* , unsigned int*, unsigned int*);
typedef  unsigned int  (*NFUNCM)       (short);
typedef  unsigned int  (*H_FUNCM)      ( unsigned int );
typedef  void          (*IND2VECFUNCM) ( short, short, unsigned int, short* ) ;
typedef  void          (*NDIM_FUNCM)   ( short, short, unsigned int, short* );

typedef unsigned long long ui64;  /*  Basop Mpy32_32_uu simulation */

/* local constants for indexing functions  c-code clarity */
#ifndef ONE_U

#define ZERO        0
#define ONE         1
#define ONE_U       1u
#define TWO         2
#define SIGNBIT      0x80000000u

#define MAXBIT      0x40000000L  /* if Word32 */
#define MINNEG      0xffffffffL

#define UDIVBY3     2863311531U


/*-------------------------------------------------------------------*
 *  local_norm_l_opt
 *
 *  rewritten version of STL norm_l for "int" fixed point normalization
 *  in floating point c-code.
 *-------------------------------------------------------------------*/
static
short local_norm_l_opt (  /* o  : shifts needed for normalization     */
    int l32var          /* i  : signed  32bit input value           */
)
{
    short l32res;


    if (l32var == (int)MINNEG )
    {
        return (32-ONE);
    }
    else
    {
        if (l32var == ZERO)
        {
            return ZERO;
        }
        else
        {
            if (l32var < ZERO)
            {
                l32var = ~l32var;
            }

            for (l32res = ZERO; l32var < (int)MAXBIT; l32res++)
            {
                l32var <<= ONE;
            }
        }
    }

    return (l32res);
}



/*-------------------------------------------------------------------*
 * floor_sqrt_exact()
 *   returns  x = floor(sqrt(input));     where (x*x) <= input
 *-------------------------------------------------------------------*/
unsigned int floor_sqrt_exact(  /* o  : floor(sqrt(input))                   */
    unsigned int input          /* i  : unsigned input  [0.. UINT_MAX/4]     */
)
{
    double _tmp;
    unsigned int x;
    if (input == ZERO)
    {
        return ZERO;
    }

    _tmp  = (double) input ;
    x     = (unsigned int)(sqrt( _tmp ));   /* floor is a part of the cast  */
    return x;
}

/*-------------------------------------------------------------------*
 * f_odd_exact_div_opt()
 *
 * find   1/(den1*2+1) * ( num1p*num2p - num3) ,
 *        if the result is known a priori to be exactly a 32bit unsigned integer
 *-------------------------------------------------------------------*/
static
unsigned int f_odd_exact_div_opt(  /* o  :  see Eq.        */
    unsigned int num1p,        /* i  :  see Eq.        */
    unsigned int num2p,        /* i  :  see Eq.        */
    unsigned int num3,         /* i  :  see Eq.        */
    unsigned int den1          /* i  :  see Eq.        */
)
{
    unsigned int tmp;
    tmp = exactdivodd[den1] * (num1p*num2p - num3);

    return (tmp);
}

/*---------------------------------------------------------------------------*
 * f_even_exact_div_opt()
 *
 *    returns  (num1p*num2p - num3 )/ den1
 *     if the result is known a priori to be exactly a 32bit unsigned integer
 *--------------------------------------------------------------------------*/
static
unsigned int f_even_exact_div_opt(  /* o  :  see Eq.        */
    unsigned int num1p,             /* i  :  see Eq.      */
    unsigned int num2p,             /* i  :  see Eq.   range should be larger than num1p   */
    unsigned int num3,              /* i  :  see Eq.                   */
    int den1               /* i  :  see Eq.                   */
)
{
    unsigned int tmp1, oddfactor, UL_tmp;
    int den1_m1;
    short even_sh;
    unsigned int UL_tmp_h;
    ui64         ULL_tmp;

    den1_m1 = den1 - ONE; /* remove top bit */
    even_sh =  (31) - local_norm_l_opt((den1_m1)^den1);  /*  NB STL operation defined for  signed positive 32 bit variable  */
    oddfactor = exactdivodd[den1_m1>>even_sh];
    even_sh  -= ONE;

    ULL_tmp  =  (ui64)num1p*(ui64)num2p;            /* use STL Mpy_32_32_uu functionality */
    UL_tmp_h =  (unsigned int)(ULL_tmp>>32);        /* high output from basicop */
    UL_tmp   =  (unsigned int) (ULL_tmp);           /* low output  from basicop */

    if(num3 > UL_tmp)
    {
        UL_tmp_h = UL_tmp_h - ONE_U;
    }
    UL_tmp = (UL_tmp - num3);                             /* can and should underflow */
    UL_tmp = (UL_tmp_h<<(32-even_sh))|(UL_tmp>>even_sh);  /* bitwise OR */

    /* use tabled modular multiplicative inverse for the odd part division */
    tmp1 = UL_tmp*oddfactor;


    return tmp1;
}

/*-------------------------------------------------------------------*
 * a_three()
 *-------------------------------------------------------------------*/
static
unsigned int a_three(           /* o:  offset for dim 3 */
    unsigned int k_val          /* i:  nb unit pulses   */
)                               /* k_val may be higher than  16 bit signed   */
{
    if( k_val )
    {
        return (ONE_U + k_val*((k_val - ONE_U) << ONE));
    }
    else
    {
        return ZERO;
    }
}

/*-------------------------------------------------------------------*
 * a_four()
 *-------------------------------------------------------------------*/
static
unsigned int a_four(       /* o:  offset for dim 4 */
    unsigned int k_val     /* i:  nb unit pulses   */
)
{
    if(k_val)
    {
        return  UDIVBY3*((k_val<<ONE)*(4 + ((k_val<<ONE) - 3)*k_val)  -   3); ;
    }
    return ZERO;
}

/*-------------------------------------------------------------------*
 * a_five()
 *-------------------------------------------------------------------*/
static
unsigned int a_five(            /* o:  offset for dim 5 */
    unsigned int k_val          /* i:  nb unit pulses   */
)
{
    unsigned int offset;

    if(k_val==0)
    {
        offset=0;
    }
    else if(k_val==1)
    {
        offset=1;
    }
    else
    {
        offset = (ONE_U + ((((((k_val - TWO)*k_val + 5)*k_val - 4)*k_val)*(UDIVBY3))<<ONE));
    }
    return offset;
}

/*-------------------------------------------------------------------*
 * direct_msize()
 *-------------------------------------------------------------------*/
static
unsigned int direct_msize(short dim_in, short k_in)
{
    unsigned int msize,k,ksq;

    msize = ONE_U; /* k==0 or dim==1 , and base */
    if(k_in)
    {
        k = (unsigned int)k_in;
        ksq=k*k;
        switch (dim_in)
        {
        case (5):
            /* k*k = 238*238 < 2^16  */
            msize += ( ((ksq*(5  +  ksq))* UDIVBY3 )<<ONE );
            break;
        case (4 ):
            msize  =   ((k*(ksq + 2))*UDIVBY3)<<TWO;
            break;
        case ( 3 ):
            msize  +=  ((ksq)<<ONE) ;
            break;
        case ( 2 ):
            msize  = k<<ONE;
            break;
        }
    }

    return msize;
}


/* update h_mem[0.. k_val_in+1] ,  with starting offsets for A+U recursion */
static
void initOffsets( short dim_in, unsigned int* h_mem, short k_val_in)
{
    unsigned int k_val_curr, k_val_prev;
    short k_val;
    h_mem[0]          =  ZERO;                                  /*    % A(=>0,k=0)      */
    h_mem[1]          =  ONE_U;                                 /*    % A(*,k=1)        */

    if(dim_in==2)
    {
        for( k_val = TWO; k_val <= k_val_in; k_val++)
        {
            h_mem[k_val] = (unsigned int)((k_val<<ONE) - ONE);    /* A(2, 2 .. k ) */
        }
        h_mem[k_val] =   (unsigned int)(k_val_in);             /* U(2,k+1) */
    }
    else  if(dim_in==3)
    {
        k_val_prev        =   ONE;
        for(k_val_curr = TWO; k_val_curr <= (unsigned int)k_val_in; k_val_curr++)
        {
            h_mem[k_val_curr]=  (ONE_U + k_val_curr*(k_val_prev<<ONE) );
            k_val_prev = k_val_curr;
        }
        h_mem[k_val_curr] = k_val_curr* k_val_prev;
    }
    else if(dim_in==4)
    {
        for(k_val_curr = TWO; k_val_curr <= (unsigned int)k_val_in; k_val_curr++)
        {
            h_mem[k_val_curr] = UDIVBY3*((k_val_curr<<ONE)*(4 + ((k_val_curr<<ONE) - 3)*k_val_curr)  -   3);
        }
        h_mem[k_val_curr] = (UDIVBY3*((k_val_curr<<ONE)*(4 + ((k_val_curr<<ONE) - 3)*k_val_curr)  -   3))>>1 ;
    }
    return;
}


/*-------------------------------------------------------------------*
 * a_fwd()
 *
 *  create offsets for A(n,k)  from lower A(n-1,k)
 *-------------------------------------------------------------------*/
static
void a_fwd(
    unsigned int   *a_in,           /* i/o: offsets   */
    short n_items          /* i  :  items, k's  */
)
{
    unsigned int a_1;
    short i,i_prev;
    unsigned int a_in0 ;         /* i  :  start column value   */

    a_in0 = ONE_U;
    i_prev=ZERO;
    for(i=ONE; i <= n_items; i++) /*basic A   fwd row  recursion */
    {
        a_1           = a_in0  +  a_in[i_prev] + a_in[i] ;   /* a_in addressed in at least two locations */
        a_in[i_prev]  = a_in0;
        a_in0         = a_1;
        i_prev        = i;          /* no real need to count as it is a ptr update */
    }
    a_in[i_prev] = a_in0;

    return;
}

/*-------------------------------------------------------------------*
 * a_bwd()
 *
 *  create offsets for A(n,k)  from higher A(n+1,k)
 *-------------------------------------------------------------------*/
static
void a_bwd(
    unsigned int   *a_in,      /* i/o: offsets   */
    short n_items     /* i:  n_items  */
)
{
    unsigned int a_1;
    unsigned int a_in0;
    short i;
    short i_prev;

    a_in0  = ZERO;
    i_prev = ZERO;
    for(i=ONE; i<=n_items; i++) /*basic A   reverse row  recursion */
    {
        a_1          = a_in[i] - a_in0 - a_in[i_prev];
        a_in[i_prev] = a_in0;
        a_in0        = a_1;
        i_prev       = i;
    }
    a_in[i_prev] = a_in0;
    return;
}

static
unsigned int  direct_row_A2U_rec_calc(short dim_in , short k_val_in, unsigned int a_km2, unsigned int a_km1)
{

    /*  U(n,k)    =  (A(n,k-2)-1)/2     +     ((2*n-1)*A(n,k-1) - A(n,k-2) )/2*(k-1)           */
    /*  U(n,k) = floor(A(n,k-2)/2) + (n*A(n,k-1) - floor(A(n,k-1)/2) - floor(A(n,k-2)/2) +1 )/(k-1) */
    /*  U(n,k) = floor(A(n,k-2)/2) + (n*A(n,k-1) - (floor(A(n,k-1)/2) + floor(A(n,k-2)/2) +1) ) /(k-1) */

    unsigned int divisor, km2_size, result;

    divisor     = (unsigned int)(k_val_in-ONE);
    km2_size = (a_km1>>ONE) + (a_km2>>ONE) + ONE_U  ;

    if(divisor&ONE_U)
    {
        /* odd */
        result = ( (a_km2>>ONE ) + f_odd_exact_div_opt((unsigned int)(dim_in), a_km1, km2_size , divisor>>ONE) );
    }
    else
    {
        /* even divisor */
        result = ( (a_km2>>ONE) + f_even_exact_div_opt((unsigned int)dim_in,a_km1, km2_size, divisor ) );
    }
    return result;
}

static
void a_u_fwd(
    unsigned int   *a_u_in,
    short k_val_in,
    short mem_size_m1
)
{
    unsigned int u_kp1_prev;
    unsigned int a_k_prev ;

    /* mem_size_m1 =  1 + k_val_in   */
    u_kp1_prev = a_u_in[mem_size_m1]; /* previous  n  U (n,k+1) value*/
    a_k_prev   = a_u_in[k_val_in];    /* previous  n  A(n,k) value*/

    a_fwd(&a_u_in[ONE], k_val_in);  /* a_u_in[k==ZERO] = zero if n>0 */

    /*      low dynamic last offset entry mixed recursion */
    /*      used for size calculation  */
    /*      U(n,k+1) = 1 + U(n-1,k+1) + U(n-1,k)        + U(n,k)                            */
    /*      U(n,k+1) = 1 + U(n-1,k+1) + (A(n-1,k)-1)/2  + (A(n,k)-1)/2                      */
    /*                  Note, A(n,k) always odd for k>0 , subtracted one always shifted out */
    /*                  assert(a_k_prev>0, a_k-curr>0)                                      */

    a_u_in[mem_size_m1] = ONE_U  +  u_kp1_prev   +  (a_k_prev>>ONE)  + (a_u_in[k_val_in] >>ONE);

    return;
}

/*-------------------------------------------------------------------*
 * nm_h_prep_opt()
 *
 * find  and return  N_MPVQ(n,k) and also offsets A(n, 0  to  k ) and  U(n,k+1).
 *-------------------------------------------------------------------*/
static
unsigned int nm_h_prep_opt(      /* o:  msize for dim     */
    short dim_in,       /* i:  dimension         */
    short k_val_in,     /* i:  nb unit pulses    */
    unsigned int   *h            /* o:  A/U offsets array */
)
{
    short mem_size_m1, k_val ;
    short dim_tmp, d_start;
    unsigned int h_saveA, h_saveB;   /*  registers   for alternating  A(n,k-1), A(n,k-2)*/
    unsigned int numDsub1;           /*  k_val_curr, k_val_prev*/;

    mem_size_m1 = k_val_in + ONE;

    if( k_val_in > TABLE_LIM_OPT )
    {
        if( dim_in >= 3 )
        {
            d_start = 3;
        }
        else
        {
            d_start = 2;
        }
        initOffsets(d_start, h, k_val_in);

        for(dim_tmp = d_start; dim_tmp < dim_in; dim_tmp++)
        {
            a_u_fwd(h, k_val_in, mem_size_m1);
        }
    }
    else
    {
        h[ZERO] = ZERO;
        h[ONE]  = ONE_U;
        numDsub1=(unsigned int) ((dim_in << ONE) - ONE);
        h[TWO] = numDsub1;

        /* interleaved odd even calls */
        h_saveA   =   numDsub1 ;
        h_saveB   =   ONE_U;
        for (k_val = 3; k_val < (mem_size_m1); k_val++ )
        {
            /* A(n,k)  =  A(n,k-2) + ((2*n-1)*A(n,k-1)-A(n,k-2)) /(k-1)  */
            /* first odd  k, even divisor */
            h_saveB   += f_even_exact_div_opt(numDsub1, h_saveA, h_saveB, k_val - ONE);
            h[k_val] = h_saveB;

            k_val++; /* next even k, odd divisor */
            if( k_val >= (mem_size_m1))
            {
                break;
            }
            h_saveA      += f_odd_exact_div_opt(numDsub1, h_saveB, h_saveA, (k_val - ONE)>>ONE);
            h[k_val]      = h_saveA;
        }
        /*  always do the last (k+1) recursion based  on  U(n,k+1) = func( A(n-2,k+1), A(n-1,k+1) )   */
        h[mem_size_m1] = direct_row_A2U_rec_calc(dim_in, mem_size_m1 , h[mem_size_m1-2], h[k_val_in]);
    }

    /*  N_MPVQ(n,k) = 1 + U(n,k+1) + U(n,k) =  1 + U(n,k+1) + (A(n,k)-1)/2 ;   */ /* A(n,k) always odd */
    return ( ONE + h[mem_size_m1]  +  (h[k_val_in]>>ONE) );
}


/*
 find_amp_split_offset_func_mem()
   find first offset  in range 0..k_val_in  that is less than ind_in
   using a tree search with direct function calls or memory iteration
*/
static
short find_amp_split_offset_func_mem(
    unsigned int   ind_in,
    short k_val_in,
    H_FUNCM h_func_ptr,          /* i: offset function pointer  */
    unsigned int   *h_mem ,
    short k_test,              /* o:  k_value  */
    unsigned int   *tmp_offset          /* o:  offset found  */
)
{
    short  not_ready, low,high ;

    low  = 0;
    high = k_val_in;
    /* split over  A(n,k)= h_mem(k), or use direct function  */
    not_ready = ONE ;
    while(not_ready)
    {
        k_test  = (low+high)>>ONE;              /*% split range in half */

        if(h_mem)
        {
            *tmp_offset = h_mem[k_test];                         /* memory search  */
        }
        else
        {
            *tmp_offset = (*h_func_ptr)((unsigned int)k_test);   /* function search. NB only line difference to the memory search*/
        }

        if(ind_in > *tmp_offset )
        {
            if(k_test  < high)
            {
                low =  1 + k_test ;
            }
            else
            {
                not_ready=0;
            }
        }
        else
        {
            if (*tmp_offset  == ind_in )
            {
                not_ready=0;
            }
            else
            {
                high  = k_test  - 1;
            }
        }
    }
    return k_test;
}


/*
  get_lead_sign()
  update  index and return leading sign
*/
static
short get_lead_sign(unsigned int *ind)
{
    short leading_sign;

    if( (*ind)&ONE_U)       /*  leading sign stored in LSB  */
    {
        leading_sign =  -1;
    }
    else
    {
        leading_sign = 1;
    }
    (*ind) = (*ind)>>ONE;

    return leading_sign;
}


/*-------------------------------------------------------------------*
 * mind2vec_one()
 *-------------------------------------------------------------------*/
static
void mind2vec_one(
    short k_val_in,       /* i:  nb unit pulses   */
    short leading_sign,   /* i: leading sign */
    unsigned int   ind,            /* i:  index            */
    short *vec_out        /* o:  pulse train      */
)
{
    ind = 0;    /* to avoid compiler warings */

    vec_out[ind] = (leading_sign*k_val_in);  /*  NB  input k_val_in can be zero  */
}

/*-------------------------------------------------------------------*
 * mind2vec_two()
 *-------------------------------------------------------------------*/
static
void mind2vec_two(
    short k_val_in,        /* i:  nb unit pulses   */
    short leading_sign,    /* i: leading sign */
    unsigned int   ind_in,          /* i:  index            */
    short *vec_out         /* o:  pulse train      */
)
{
    unsigned int ind_tmp;
    short val1;

    if(k_val_in > 0)    /* k_val_in check   */
    {
        if (ind_in==0)
        {
            vec_out[0] = (leading_sign*k_val_in);
            vec_out[1] = 0;
        }
        else if (ind_in == ( (unsigned int)(k_val_in<<ONE) - ONE_U) )
        {
            vec_out[0] = 0;
            vec_out[1] = (leading_sign*k_val_in);
        }
        else
        {
            ind_tmp       =  ind_in - ONE_U;
            val1          =  (short)(ONE_U + (ind_tmp>>ONE));
            vec_out[0]    =  leading_sign*(k_val_in - val1) ;

            if(ind_tmp&ONE_U)
            {
                vec_out[1] =  -val1 ;
            }
            else
            {
                vec_out[1] =  val1;
            }
        }
    }
}

static
short setval_update_sign(
    short k_delta,
    short k_max_local,
    short *leading_sign,
    unsigned int   *ind_in,
    short *vec_out
)
{
    if(k_delta != 0 )
    {
        *vec_out       = (*leading_sign)*k_delta;
        *leading_sign  =  get_lead_sign(ind_in);
        k_max_local    = k_max_local-k_delta ;
    }
    return k_max_local;
}

/*-------------------------------------------------------------------*
 * mind2vec_three()
 *-------------------------------------------------------------------*/
static
void mind2vec_three(
    short k_max_local,     /* i:  nb unit pulses   */
    short leading_sign,    /* i: leading sign */
    unsigned int   ind_in,          /* i:  index            */
    short *vec_out         /* o:  pulse train      */
)
{
    short k_delta ;
    unsigned int   acc_val;

    /*
        use direct calculation of first amplitude
        (to find amplitudes faster than using split or linear iteration)
    */
    if(ind_in==0)
    {
        vec_out[0] = leading_sign*k_max_local;
    }
    else
    {
        acc_val = ((ONE_U  + floor_sqrt_exact((ind_in<<ONE) - ONE_U))>>ONE );   /* in BASOP use approximation + search for exact sqrt )*/

        k_delta = k_max_local - (short)acc_val;
        ind_in -= a_three(acc_val);                  /* remove amplitude offset A(n,k_acc) */

        k_max_local = setval_update_sign( k_delta, k_max_local, &leading_sign,&ind_in,vec_out);

        mind2vec_two( k_max_local, leading_sign, ind_in ,&vec_out[1] );
    }
    return;
}

/*-------------------------------------------------------------------*
 * mind2vec_direct ,
    general function for direct decoding using direct functions
    (no  memory recursion)
 *-------------------------------------------------------------------*/
static
void mind2vec_direct(
    short k_max_local,         /* i:  nb unit pulses   */
    short leading_sign,        /* i: leading sign  */
    unsigned int   ind,                 /* i:  index            */
    H_FUNCM h_func_ptr,          /* i : offset function */
    NDIM_FUNCM nd_func_ptr,         /* i : next dimension function  */
    short *vec_out             /* o:  pulse train      */
)
{
    short    k_delta, k_test=0;
    unsigned int   tmp_offset;

    if(ind==0)
    {
        vec_out[0]               =  leading_sign*k_max_local;
    }
    else
    {
        k_test = find_amp_split_offset_func_mem(ind,k_max_local, h_func_ptr , NULL, k_test, &tmp_offset);

        k_delta  =  k_max_local - k_test;
        ind      =  ind - tmp_offset;                           /* %  remove amplitude offset A(n,k_acc) */
        k_max_local = setval_update_sign( k_delta, k_max_local, &leading_sign, &ind,vec_out);
        (*nd_func_ptr)( k_max_local, leading_sign, ind , &vec_out[1] );
    }
    return;
}

/*-------------------------------------------------------------------*
 * mind2vec_four()
 *-------------------------------------------------------------------*/
static
void mind2vec_four(
    short k_val_in,        /* i:  nb unit pulses   */
    short leading_sign,    /* i: leading sign */
    unsigned int   ind_in,          /* i:  index            */
    short *vec_out         /* o:  pulse train      */
)
{
    mind2vec_direct(k_val_in,leading_sign, ind_in, a_four, mind2vec_three, vec_out);
    return;
}

/*-------------------------------------------------------------------*
 * mind2vec_five()
 *-------------------------------------------------------------------*/
static
void mind2vec_five(
    short k_val_in ,       /* i:  nb unit pulses   */
    short leading_sign,    /* i: leading sign  */
    unsigned int   ind_in,          /* i:  index            */
    short *vec_out         /* o:  pulse train      */
)
{
    mind2vec_direct(k_val_in,leading_sign, ind_in, a_five, mind2vec_four, vec_out);
    return;
}


/*-------------------------------------------------------------------*
 * mind2vec()
 *-------------------------------------------------------------------*/
static
void mind2vec(
    short dim_in,              /* i:  dimension        */
    short k_max_local,         /* i:  nb unit pulses   */
    short leading_sign,        /* i: leading sign  */
    unsigned int   ind,                 /* i:  index            */
    short *vec_out,            /* o:  pulse train      */
    unsigned int   *h_in                /* i:  offset vector   A=1+2U  */
)
{
    short pos;
    short k_acc, k_delta;
    unsigned int tmp_offset;

    k_acc = k_max_local;

    pos = ZERO;
    while (pos < dim_in)        /* first to last position decoding */
    {
        if(ind == 0)
        {
            vec_out[pos] = leading_sign*k_max_local;
            break;              /* "fast" recursion exit */
        }
        else
        {
            {
                /* linear magnitude search */
                k_acc       = k_max_local;
                tmp_offset  = h_in[k_acc];
                while(tmp_offset > ind)
                {
                    k_acc        = k_acc - 1;
                    tmp_offset   = h_in[k_acc];
                }
            }
            k_delta   = k_max_local - k_acc; /* amplitude decoding */
        }
        ind          = ind - tmp_offset;     /* remove amplitude index offset A(n,k_acc) */

        k_max_local = setval_update_sign( k_delta, k_max_local, &leading_sign, &ind, &vec_out[pos]);

        /* move from  A(n,kmax) to A(n-1, k_max_local), */
        a_bwd( h_in,k_max_local + 1 );  /* [0..k_max_local], no need to calculate U(n,k_max_local+1) in index decoding   */
        pos    = pos + 1;
    }
    return;
}


/*-------------------------------------------------------------------*
 * vec2mind_one()
 *-------------------------------------------------------------------*/
static
void  vec2mind_one(
    const short *vec_in,             /* i  : PVQ abs pulse train      */
    short *k_val_out_ptr ,     /* o  : number of unit pulses    */
    unsigned int   *next_sign_ind,      /* i/o: next sign ind            */
    unsigned int   *ind                 /* o  : MPVQ index               */
)
{
    *k_val_out_ptr = -1;         /* just to avoid compiler warnings */

    *ind  =   ZERO;
    /* *k_val_out_ptr =  (short)abs(*vec_in);  */      /* dim==1, function not called recursively */
    *next_sign_ind = (unsigned int)(*vec_in < ZERO);   /* single sign always pushed out of MPVQ */
    return;
}

/*-------------------------------------------------------------------*
* vec2mind_two()
*-------------------------------------------------------------------*/
static
void vec2mind_two(
    const short *vec_in,             /* i : PVQ  pulse train          */
    short *k_val_out_ptr,      /* o : number of unit pulses    */
    unsigned int   *next_sign_ind,      /* i/o: next sign ind */
    unsigned int   *ind                 /* o: MPVQ index */
)
{
    unsigned int lead_sign_ind;
    short abs0,abs1,abs01;

    abs0 = (short) abs(vec_in[0]);
    abs1 = (short) abs(vec_in[1]);

    abs01          = abs0+abs1;
    *k_val_out_ptr = abs01;

    if(abs01==0)            /* zeroes  can happen in a recursive encoding call */
    {
        *next_sign_ind = SIGNBIT;
        *ind = ZERO;
    }
    else
    {
        *next_sign_ind=0;
        if(abs1 == 0)
        {
            *ind           = ZERO;
        }
        else if(abs0 == 0)
        {
            *ind            =  (unsigned int)(abs1<<ONE) - ONE_U;
            *next_sign_ind = 1;
        }
        else
        {
            lead_sign_ind   =  (unsigned int)(vec_in[1]<0); /*% sign always shifted to first pos */
            *ind            = ONE_U + ((unsigned int)(abs1-1)<<ONE)  +  lead_sign_ind;

        }
        *next_sign_ind= (unsigned int)(vec_in[*next_sign_ind]<0);
    }
    return;
}


/*-------------------------------------------------------------------*
 *  vec2mind_gen345( vec_in kval,   next_dim_func , offset_func,....)
 *   generic call saves  PROM ,
 *-------------------------------------------------------------------*/
static
void vec2mind_gen345(
    const short *vec_in,                 /* i : PVQ abs pulse train          */
    short *k_val_out_ptr,          /* o : number of unit pulses    */
    unsigned int   *next_sign_ind ,         /* i/o: next sign ind */
    unsigned int   *index ,                 /* o: MPVQ index */
    VEC2INDFUNCM vec2indfunc_ptr,         /* i: */
    H_FUNCM a_func_ptr               /*i:  pffset function  */
)
{
    short  tmp_val;

    tmp_val = vec_in[ZERO];
    (*vec2indfunc_ptr)(&vec_in[ONE], k_val_out_ptr, next_sign_ind, index);

    if( ((*next_sign_ind&SIGNBIT)==0 && tmp_val!=0)!= ZERO )
    {
        *index   = (*index<<ONE)  +  *next_sign_ind;
    }

    if(tmp_val != 0)       /* push sign */
    {
        *next_sign_ind  =  (unsigned int)(tmp_val<0);
    }

    *index +=   (*a_func_ptr)((unsigned int)*k_val_out_ptr);
    if( tmp_val > 0 )
    {
        *k_val_out_ptr += tmp_val;
    }
    else
    {
        *k_val_out_ptr -= tmp_val;
    }
    return ;
}


/*-------------------------------------------------------------------*
 * vec2mind_three()
 *-------------------------------------------------------------------*/
static
void vec2mind_three(
    const short *vec_in,             /* i : PVQ   pulse train          */
    short *k_val_out_ptr,      /* o : number of unit pulses    */
    unsigned int   *next_sign_ind,      /* i/o: next sign ind */
    unsigned int   *index               /* o: MPVQ index */
)
{


    vec2mind_gen345(vec_in,k_val_out_ptr, next_sign_ind, index, vec2mind_two, a_three);

    return ;
}

/*-------------------------------------------------------------------*
 * vec2mind_four()
 *-------------------------------------------------------------------*/
static
void vec2mind_four(
    const short *vec_in,                 /* i : PVQ  pulse train          */
    short *k_val_out_ptr,          /* o : number of unit pulses    */
    unsigned int   *next_sign_ind,          /* i/o: next sign ind */
    unsigned int   *index                   /* o: MPVQ index */
)
{

    vec2mind_gen345(vec_in,k_val_out_ptr, next_sign_ind, index, vec2mind_three, a_four);

    return ;
}

/*-------------------------------------------------------------------*
 * vec2mind_five()
 *-------------------------------------------------------------------*/
static
void vec2mind_five(
    const short *vec_in,                 /* i : PVQ abs pulse train          */
    short *k_val_out_ptr,          /* o : number of unit pulses    */
    unsigned int   *next_sign_ind,          /* i/o: next sign ind */
    unsigned int   *index                   /* o: MPVQ index */
)
{

    vec2mind_gen345(vec_in,k_val_out_ptr, next_sign_ind, index, vec2mind_four, a_five);

    return ;
}



/*-------------------------------------------------------------------*
 * vec2mind()
 *-------------------------------------------------------------------*/
static
void vec2mind(
    short dim_in,              /* i :  dim                       */
    short k_val_in,            /* i :  number of unit pulses     */
    const short *vec_in,             /* i :  PVQ pulse train           */
    unsigned int   *next_sign_ind,      /* o :  pushed leading sign       */
    unsigned int   *index  ,            /* o :  MPVQ index                */
    unsigned int   *N_MPVQ_ptr,         /* o :  size(N_MPVQ(dim,K_val_in))*/
    unsigned int   *h_mem               /* o :  offsets                   */
)
{
    short  pos,  mem_size_m1 ;
    short  k_val_acc ;
    short  tmp_val;


    mem_size_m1  = k_val_in + ONE;

    *next_sign_ind = SIGNBIT; /*  % should always be 0 or 1 out,  highest bit set signals no sign found yet*/


    pos       =  dim_in - 2;                        /*  % address 2nd last sample */
    vec2mind_two(&vec_in[pos],&k_val_acc,next_sign_ind ,index);
    initOffsets( 3, h_mem, k_val_in) ;


    for (pos--; pos>=0; pos--)
    {
        /*
        %   Check if the  leading sign 'bit' is to be added
        */
        tmp_val = vec_in[pos];
        if(  ((*next_sign_ind&SIGNBIT)==0 && (tmp_val != 0)) )
        {
            *index = (*index<<ONE) + (*next_sign_ind);
        }

        /* % push sign fwd,  for encoding in the next non_zero position  */
        if(tmp_val != 0)
        {
            *next_sign_ind  = (unsigned int)(tmp_val<0);
        }

        /*%  add  indexing offset up to this reverse (r_l) accumulated  amplitude point */
        *index += h_mem[k_val_acc];            /* % k_val_acc==0 ==>0 */

        k_val_acc += (short)abs(tmp_val);/*% now increase acc k value for next N */


        if(pos)
        {
            a_u_fwd(h_mem, k_val_in ,mem_size_m1);
            /*%  update A(n,k=1:k_val_in) and U(n,k_val_in+1) , NB here  (k_val_in + 2 elements always has to be updated */
        }
    }

    /*   size is needed for the subseqent arithmetic encoding/transmission of the index.        */
    /* use relation N_MPVQ(n,K) =  1 +  (A(n, K)-1)/2 + U(n, 1 + K)                        */
    /* =            N_MPVQ(n,K) =  1 +  (A(n, K)>>1)  + U(n, 1 + K) ,  as A(n,K) is odd)   */

    *N_MPVQ_ptr  =   ONE_U   +  (h_mem[k_val_acc]>>ONE)  +  h_mem[ mem_size_m1 ]  ; /*  total size size */

    return;
}

/*--------------------------------------------------------------------------*
 * mpvq_encode_vec()
 *
 * returns struct with  leading sign index, MPVQ-index ,  dim   and  N_MPVQ
 *-------------------------------------------------------------------------*/

PvqEntry mpvq_encode_vec(       /* o : leading_sign_index, index, size, k_val        */
    const short *vec_in,        /* i : signed pulse train        */
    short dim_in,         /* i : dimension                 */
    short k_val_local     /* i : nb unit pulses            */
)
{
    PvqEntry result;
    unsigned int h_mem[1+KMAX_NON_DIRECT+1];    /* allocate max offset memory  for dim 6    */
    /* OPT: actually only 1+k_val_in+1 needed ) */
    unsigned int lead_sign_ind;

    VEC2INDFUNCM  vec2mind_f[1+N_OPT] = { (VEC2INDFUNCM)NULL, vec2mind_one, vec2mind_two, vec2mind_three, vec2mind_four, vec2mind_five };
    /* VEC2INDFUNCM can be a static global  */


    result.k_val =  k_val_local;
    result.dim   =  dim_in;

    /* NB, k_val_local  may be changed in some sub encoding routines */
    if( dim_in > N_OPT)     /* use the generic dimension  function */
    {
        vec2mind(dim_in, k_val_local, vec_in, &lead_sign_ind,  &result.index, &result.size,  h_mem);
    }
    else                    /* if (dim_in<=N_OPT) */
    {
        (vec2mind_f[dim_in])(vec_in, &k_val_local, &lead_sign_ind, &result.index);
        result.size = direct_msize(dim_in, k_val_local);    /* k_val_local not used for  dim==1 */
    }
    result.lead_sign_ind=(short)lead_sign_ind;

    return result;
}

/*-------------------------------------------------------------------*
 * get_size_mpvq_calc_offset()
 *
 *  unsigned int h_mem[1 + KMAX +1 ];
 *   example using fixed size of offset vector input help variable
 *-------------------------------------------------------------------*/

PvqEntry get_size_mpvq_calc_offset(     /* o : size, dim, k_val        */
    short dim_in,              /* i : dimension                */
    short k_val_in,            /* i : nb unit pulses           */
    unsigned int* h_mem                 /* o : offsets                  */
)
{
    PvqEntry entry;

    entry.dim   = dim_in;
    entry.k_val = k_val_in;
    entry.index         = 0U; /* avoid gcc warning in struct passing */
    entry.lead_sign_ind = 0;  /* avoid gcc warning in struct passing */
    if(dim_in > N_OPT )      /* non-direct solutions,  use A+U relation */
    {
        entry.size = nm_h_prep_opt(entry.dim, entry.k_val, h_mem);
    }
    else
    {
        entry.size =  direct_msize(dim_in, entry.k_val);  /* ues equations,  h_mem is not used */
    }


    return entry;
}

/*-------------------------------------------------------------------*
 * mpvq_decode_vec()
 *-------------------------------------------------------------------*/

void mpvq_decode_vec(               /* o :  void                        */
    const PvqEntry *entry,          /* i :  sign_ind, index, dim, k_val */
    unsigned int   *h_mem,          /* i :  A/U offsets                 */
    short *vec_out         /* o :  pulse train                 */
)
{
    short i, leading_sign;
    IND2VECFUNCM  mind2vec_f[N_OPT+1] = { (IND2VECFUNCM)NULL, mind2vec_one, mind2vec_two, mind2vec_three, mind2vec_four, mind2vec_five };
    /*IND2VECFUNCM vector can be static global  */


    for(i=0; i<entry->dim; i++)
    {
        vec_out[i]=ZERO;  /* set all of short output vector to zero */
    }

    leading_sign    =  1;
    if(entry->lead_sign_ind)
    {
        leading_sign = -1;
    }

    if(entry->k_val != 0)
    {
        if(entry->dim > N_OPT )  /* N_OPT  */
        {
            /*  generic */
            mind2vec(entry->dim, entry->k_val, leading_sign, entry->index, vec_out, h_mem);
        }
        else
        {
            /*  specialized functions */
            (mind2vec_f[entry->dim])(entry->k_val, leading_sign, entry->index, vec_out);
        }
    }
    return;
}

#ifdef ONE_U
#undef ZERO
#undef ONE
#undef ONE_U
#undef TWO
#undef MAXBIT
#undef MINNEG
#undef SIGNBIT
#undef UDIVBY3
#endif



#endif
