/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <assert.h>
#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"

/*---------------------------------------------------------------------*
 * Local constants
 *---------------------------------------------------------------------*/

#define FORMAT_POST_FILT_G1 0.75f           /*0.75f*/ /*denominator 0.9,0.75,0.15,0.9*/


/*--------------------------------------------------------------------------
 * Local functions
 *--------------------------------------------------------------------------*/

static void Dec_postfilt( const short L_subfr, PFSTAT *pfstat, const int t0, const float *signal_ptr, const float *coeff,
                          float *sig_out, const float gamma1, const float gamma2, const float gain_factor, const short disable_hpf );

static void pst_ltp( const int t0, const float *ptr_sig_in, float *ptr_sig_pst0, float gain_factor, const short L_subfr );

static void search_del( const int t0, const float *ptr_sig_in, int *ltpdel, int *phase, float *num_gltp, float *den_gltp,
                        float *y_up, int *off_yup, const short L_subfr );

static void filt_plt( const float *s_in, const float *s_ltp, float *s_out, const float gain_plt, const short L_subfr );

static void compute_ltp_l( const float *s_in, const int ltpdel, const int phase, float *y_up, float *num, float *den, const short L_subfr );

static int select_ltp( const float num1, const float den1, const float num2, const float den2 );

static void modify_pst_param( const float psf_lp_noise, float *g1, float *g2, const short coder_type, float *gain_factor );

static void Dec_formant_postfilt( PFSTAT *pfstat, const float *signal_ptr, const float *coeff, float *sig_out,
                                  const float gamma1, const float gamma2, const short l_subfr );

/*--------------------------------------------------------------------------*
 *  Function  Init_post_filter()
 *
 *  Post-filter initialization
 *--------------------------------------------------------------------------*/

void Init_post_filter(
    PFSTAT * pfstat         /* i  : post-filter state memories  */
)
{
    /* It is off by default */
    pfstat->on = 0;

    /* Reset */
    pfstat->reset = 0;

    /* Initialize arrays and pointers */
    set_zero( pfstat->mem_pf_in, L_SUBFR );

    /* res2 =  A(gamma2) residual */
    set_zero( pfstat->mem_res2, DECMEM_RES2 );

    /* 1/A(gamma1) memory */
    set_zero( pfstat->mem_stp, L_SUBFR );

    /* null memory to compute i.r. of A(gamma2)/A(gamma1) */
    set_zero( pfstat->mem_zero, M );

    /* for gain adjustment */
    pfstat->gain_prec = 1.0f;

    return;
}


/*--------------------------------------------------------------------------
 * nb_post_filt()
 *
 * Main routine to perform post filtering of NB signals
 *--------------------------------------------------------------------------*/

void nb_post_filt(
    const short L_frame,       /* i  : frame length                       */
    const short L_subfr,       /* i  : sub-frame length                   */
    PFSTAT *pfstat,       /* i/o: Post filter related memories       */
    float *psf_lp_noise, /* i/o: long term noise energy             */
    const float tmp_noise,     /* i  : noise energy                       */
    float *synth,        /* i/o: synthesis                          */
    const float *Aq,           /* i  : LP filter coefficient              */
    const float *pitch_buf,    /* i  : floating pitch for each subframe   */
    const short coder_type,    /* i  : coder_type                         */
    const short BER_detect,    /* i  : BER detect flag                    */
    const short disable_hpf    /* i  : flag to disabled HPF               */
)
{
    short t0_first, i, j;
    const float *p_Aq;
    float *pf_in, post_G1, post_G2, gain_factor;
    float pf_in_buffer[M+L_FRAME16k];

    if(!BER_detect)
    {
        /* update long-term background noise energy during inactive frames */
        if( coder_type == INACTIVE )
        {
            *psf_lp_noise = 0.95f * *psf_lp_noise + 0.05f * tmp_noise;
        }
    }

    /* set post-filter input */
    modify_pst_param( *psf_lp_noise, &post_G1, &post_G2, coder_type, &gain_factor );


    if( pfstat->reset )
    {
        set_zero( pfstat->mem_res2, DECMEM_RES2 );
        mvr2r( &synth[L_frame-L_SYN_MEM], pfstat->mem_pf_in, L_SYN_MEM );
        mvr2r( &synth[L_frame-L_SYN_MEM], pfstat->mem_stp, L_SYN_MEM );
        pfstat->gain_prec = 1.0f;
        pfstat->reset = 0;
        return;
    }

    pf_in = &pf_in_buffer[M];
    mvr2r( pfstat->mem_pf_in+L_SYN_MEM-M, &pf_in[-M], M );
    mvr2r( synth, pf_in, L_frame );
    mvr2r( &synth[L_frame - L_SYN_MEM], pfstat->mem_pf_in, L_SYN_MEM );

    /* deactivation of the post filter in case of AUDIO because it causes problems to singing sequences */
    if( coder_type == AUDIO )
    {
        post_G1 = 1.f;
        post_G2 = 1.f;
        gain_factor = 1.f;
    }

    /* run the post filter */
    p_Aq = Aq;
    for( i=0, j=0; i<L_frame; i+=L_subfr, j++ )
    {
        t0_first = (short)(pitch_buf[j] + 0.5f);

        Dec_postfilt( L_subfr, pfstat, t0_first, &pf_in[i], p_Aq, &synth[i], post_G1, post_G2, gain_factor, disable_hpf );

        p_Aq += (M+1);
    }

    return;
}


/*--------------------------------------------------------------------------
 *  formant_post_filt:
 *
 *  WB and SWB formant post-filtering
 *--------------------------------------------------------------------------*/

void formant_post_filt(
    PFSTAT *pfstat,        /* i/o: Post filter related memories    */
    float *synth_in,      /* i  : 12k8 synthesis                  */
    const float *Aq,            /* i  : LP filter coefficient           */
    float *synth_out,     /* i/o: input signal                    */
    const short L_frame,        /* i  : frame length                    */
    const short L_subfr,        /* i  : sub-frame length                */
    const float lp_noise,       /* i  : background noise energy         */
    const long  rate,           /* i  : bit-rate                        */
    const short off_flag        /* i  : off flag                        */
)
{
    short i_subfr;
    const float *p_Aq;
    float post_G1, post_G2;

    /*default parameter for noisy speech and high bit-rates*/
    if ( L_frame == L_FRAME )
    {
        post_G2 = 0.7f;
        if( lp_noise < LP_NOISE_THRESH )    /* Clean speech */
        {
            if( rate < ACELP_13k20 )   /*Low rates*/
            {
                post_G1 = 0.8f;
            }
            else if( rate < ACELP_24k40 )
            {
                post_G1 = 0.75f;
            }
            else
            {
                post_G1 = 0.72f;
            }
        }
        else    /*Noisy speech*/
        {
            if( rate < ACELP_15k85 ) /*Low rates*/
            {
                post_G1 = 0.75f;
            }
            else    /*High rates*/
            {
                post_G1 = 0.7f;
            }
        }
    }
    else
    {
        post_G2 = 0.76f;
        if ( lp_noise >= LP_NOISE_THRESH )
        {
            post_G1 = 0.76f;
        }
        else if( rate == ACELP_13k20 )
        {
            post_G1 = 0.82f;
        }
        else if( rate == ACELP_16k40 )
        {
            post_G1 = 0.80f;
        }
        else if( rate == ACELP_24k40 || rate == ACELP_32k )
        {
            post_G1 = 0.78f;
        }
        else
        {
            post_G1 = 0.76f;
        }
    }

    /* Switch off post-filter*/
    if( off_flag )
    {
        post_G1 = post_G2;
    }

    /* Reset post filter */
    if( pfstat->reset )
    {
        pfstat->reset = 0;
        mvr2r( &synth_in[L_frame-L_SYN_MEM], pfstat->mem_pf_in, L_SYN_MEM);
        mvr2r( &synth_in[L_frame-L_SYN_MEM], pfstat->mem_stp, L_SYN_MEM );
        pfstat->gain_prec = 1.f;
        mvr2r( synth_in,synth_out, L_frame );

        return;
    }

    /* input memory*/
    mvr2r( pfstat->mem_pf_in, synth_in-L_SYN_MEM, L_SYN_MEM);
    mvr2r( &synth_in[L_frame-L_SYN_MEM], pfstat->mem_pf_in, L_SYN_MEM);

    /* run the post filter */
    p_Aq = Aq;
    for( i_subfr=0; i_subfr<L_frame; i_subfr+=L_subfr)
    {
        Dec_formant_postfilt( pfstat, &synth_in[i_subfr], p_Aq, &synth_out[i_subfr],post_G1, post_G2, L_subfr );

        p_Aq += (M+1);
    }

    return;
}


/*----------------------------------------------------------------------------
 * Dec_postfilt()
 *
 * Adaptive postfilter main function
 *   Short-term postfilter :
 *     Hst(z) = Hst0(z) Hst1(z)
 *     Hst0(z) = 1/g0 A(gamma2)(z) / A(gamma1)(z)
 *     if {hi} = i.r. filter A(gamma2)/A(gamma1) (truncated)
 *     g0 = SUM(|hi|) if > 1
 *     g0 = 1. else
 *     Hst1(z) = 1/(1 - |mu|) (1 + mu z-1)
 *     with mu = k1 * gamma3
 *     k1 = 1st parcor calculated on {hi}
 *     gamma3 = gamma3_minus if k1<0, gamma3_plus if k1>0
 *   Long-term postfilter :
 *     harmonic postfilter :   H0(z) = gl * (1 + b * z-p)
 *       b = gamma_g * gain_ltp
 *       gl = 1 / 1 + b
 *     computation of delay p on A(gamma2)(z) s(z)
 *     sub optimal search
 *       1. search around 1st subframe delay (3 integer values)
 *       2. search around best integer with fract. delays (1/8)
 *----------------------------------------------------------------------------*/

static void Dec_postfilt(
    const short L_subfr,      /* i  : sub-frame length                          */
    PFSTAT *pfstat,      /* i/o: Post filter related memories              */
    const int   t0,           /* i  : pitch delay given by coder                */
    const float *signal_ptr,  /* i  : input signal (pointer to current subframe */
    const float *coeff,       /* i  : LPC coefficients for current subframe     */
    float *sig_out,     /* o  : postfiltered output                       */
    const float gamma1,       /* i  : short term postfilt. den. weighting factor*/
    const float gamma2,       /* i  : short term postfilt. num. weighting factor*/
    const float gain_factor,  /* i  : Gain Factor                               */
    const short disable_hpf   /* i  : flag to disable HPF                       */
)
{
    float apond1[M+1];         /* s.t. denominator coeff. */
    float apond2[LONG_H_ST];
    float sig_ltp[L_SUBFR+1];  /* H0 output signal */
    float res2[SIZ_RES2];
    float *sig_ltp_ptr;
    float *res2_ptr;
    float *ptr_mem_stp;
    float parcor0;

    /* Init pointers and restore memories */
    res2_ptr = res2 + DECMEM_RES2;
    ptr_mem_stp = pfstat->mem_stp + L_SYN_MEM - 1;
    mvr2r( pfstat->mem_res2, res2, DECMEM_RES2 );

    /* Compute weighted LPC coefficients */
    weight_a( coeff, apond1, gamma1, M );
    weight_a( coeff, apond2, gamma2, M );
    set_f( &apond2[M+1], 0, LONG_H_ST-(M+1) );

    /* Compute A(gamma2) residual */
    residu( apond2, M, signal_ptr, res2_ptr, L_subfr );

    /* Harmonic filtering */
    sig_ltp_ptr = sig_ltp + 1;

    if( !disable_hpf )
    {
        pst_ltp( t0, res2_ptr, sig_ltp_ptr, gain_factor, L_subfr );
    }
    else
    {
        mvr2r( res2_ptr, sig_ltp_ptr, L_subfr );
    }

    /* Save last output of 1/A(gamma1)  */
    /* (from preceding subframe)        */
    sig_ltp[0] = *ptr_mem_stp;

    /* Controls short term pst filter gain and compute parcor0   */
    calc_st_filt( apond2, apond1, &parcor0, sig_ltp_ptr, pfstat->mem_zero, L_subfr, -1 );

    syn_filt( apond1, M,sig_ltp_ptr, sig_ltp_ptr, L_subfr, pfstat->mem_stp+L_SYN_MEM-M, 0 );
    mvr2r( sig_ltp_ptr+L_SUBFR-L_SYN_MEM, pfstat->mem_stp, L_SYN_MEM );

    /* Tilt filtering */
    filt_mu( sig_ltp, sig_out, parcor0, L_subfr, -1 );

    /* Gain control */
    scale_st( signal_ptr, sig_out, &(pfstat->gain_prec), L_subfr, -1 );

    /* Update for next subframe */
    mvr2r( &res2[L_subfr], pfstat->mem_res2, DECMEM_RES2 );

    return;
}



/*----------------------------------------------------------------------------
 * Dec_formant_postfilt
 *
 * Post - adaptive postfilter main function
 *   Short term postfilter :
 *     Hst(z) = Hst0(z) Hst1(z)
 *     Hst0(z) = 1/g0 A(gamma2)(z) / A(gamma1)(z)
 *     if {hi} = i.r. filter A(gamma2)/A(gamma1) (truncated)
 *     g0 = SUM(|hi|) if > 1
 *     g0 = 1. else
 *     Hst1(z) = 1/(1 - |mu|) (1 + mu z-1)
 *     with mu = k1 * gamma3
 *     k1 = 1st parcor calculated on {hi}
 *     gamma3 = gamma3_minus if k1<0, gamma3_plus if k1>0
 *----------------------------------------------------------------------------*/

static void Dec_formant_postfilt(
    PFSTAT *pfstat,        /* i/o: states strucure                           */
    const float *signal_ptr,    /* i  : input signal (pointer to current subframe */
    const float *coeff,         /* i  : LPC coefficients for current subframe     */
    float *sig_out,       /* o  : postfiltered output                       */
    const float gamma1,         /* i  : short term postfilt. den. weighting factor*/
    const float gamma2,         /* i  : short term postfilt. num. weighting factor*/
    const short l_subfr         /* i  : subframe length                           */
)
{
    /* Local variables and arrays */
    float apond1[M+1];         /* s.t. denominator coeff. */
    float apond2[LONG_H_ST];
    float res2[L_SUBFR];
    float resynth[L_SUBFR+1];
    float parcor0;

    /* Compute weighted LPC coefficients */
    weight_a( coeff, apond1, gamma1, M );
    weight_a( coeff, apond2, gamma2, M );

    set_zero( &apond2[M+1], LONG_H_ST-(M+1) );

    /* Compute A(gamma2) residual */
    residu( apond2, M, signal_ptr, res2, l_subfr );

    /* Controls short term pst filter gain and compute parcor0 */
    calc_st_filt( apond2, apond1, &parcor0, res2, pfstat->mem_zero, l_subfr, -1 );

    /* 1/A(gamma1) filtering, mem_stp is updated */
    resynth[0] = *(pfstat->mem_stp + L_SYN_MEM - 1);

    syn_filt( apond1, M,res2, &(resynth[1]), l_subfr, pfstat->mem_stp+L_SYN_MEM-M, 0 );

    mvr2r( &(resynth[1])+l_subfr-L_SYN_MEM, pfstat->mem_stp, L_SYN_MEM);

    /* Tilt filtering */
    filt_mu( resynth, sig_out, parcor0, l_subfr, -1 );

    /* Gain control */
    scale_st( signal_ptr, sig_out, &pfstat->gain_prec, l_subfr, -1 );

    return;
}


/*----------------------------------------------------------------------------
 * pst_ltp()
 *
 * Perform harmonic postfilter
 *----------------------------------------------------------------------------*/

static void pst_ltp(
    const   int t0,               /* i  : pitch delay given by coder          */
    const float *ptr_sig_in,      /* i  : postfilter input filter (residu2)   */
    float *ptr_sig_pst0,    /* o  : harmonic postfilter output          */
    float gain_factor,      /* i  : gain factor                         */
    const short L_subfr           /* i  : sub-frame length                    */
)
{
    int ltpdel, phase;
    float num_gltp, den_gltp;
    float num2_gltp, den2_gltp;
    float gain_plt;
    float y_up[SIZ_Y_UP];
    const float *ptr_y_up;
    int off_yup;

    /* Suboptimal delay search */
    search_del( t0, ptr_sig_in, &ltpdel, &phase, &num_gltp, &den_gltp, y_up, &off_yup, L_subfr );

    if( num_gltp == 0.0f )
    {
        mvr2r( ptr_sig_in, ptr_sig_pst0, L_subfr );
    }
    else
    {
        if( phase == 0 )
        {
            ptr_y_up = ptr_sig_in - ltpdel;
        }
        else
        {
            /* filtering with long filter */
            compute_ltp_l( ptr_sig_in, ltpdel, phase, ptr_sig_pst0, &num2_gltp, &den2_gltp, L_subfr );

            if( select_ltp( num_gltp, den_gltp, num2_gltp, den2_gltp ) == 1 )
            {
                /* select short filter */
                ptr_y_up = y_up + ((phase - 1) * (L_subfr+1) + off_yup);
            }
            else
            {
                /* select long filter */
                num_gltp = num2_gltp;
                den_gltp = den2_gltp;
                ptr_y_up = ptr_sig_pst0;
            }
        }

        if( num_gltp >= den_gltp )
        {
            /* beta bounded to 1 */
            gain_plt = MIN_GPLT;
        }
        else
        {
            gain_plt = den_gltp / (den_gltp + ((float)0.5) * num_gltp);
        }

        /* decrease gain in noisy condition */
        gain_plt += ((1.0f-gain_plt)*gain_factor);

        /* filtering by H0(z) = harmonic filter */
        filt_plt( ptr_sig_in, ptr_y_up, ptr_sig_pst0, gain_plt, L_subfr );
    }

    return;
}

/*----------------------------------------------------------------------------
 * search_del()
 *
 * Computes best (shortest) integer LTP delay + fine search
 *---------------------------------------------------------------------------*/

static void search_del(
    const   int t0,            /* i  : pitch delay given by coder       */
    const float *ptr_sig_in,   /* i  : input signal (with delay line)   */
    int *ltpdel,       /* o  : delay = *ltpdel - *phase / f_up  */
    int *phase,        /* o  : phase                            */
    float *num_gltp,     /* o  : numerator of LTP gain            */
    float *den_gltp,     /* o  : denominator of LTP gain          */
    float *y_up,         /* o  : LT delayed signal if fract. delay*/
    int *off_yup,      /* o  : offset in y_up                   */
    const short L_subfr        /* i  : sub-frame length                 */
)
{
    const float *ptr_h;
    float tab_den0[F_UP_PST - 1], tab_den1[F_UP_PST - 1];
    float *ptr_den0, *ptr_den1;
    const float *ptr_sig_past, *ptr_sig_past0;
    const float *ptr1;
    int i, n, ioff, i_max;
    float ener, num, numsq, den0, den1;
    float den_int, num_int;
    float den_max, num_max, numsq_max;
    int phi_max;
    int lambda, phi;
    float temp0, temp1;
    float *ptr_y_up;

    /*-------------------------------------
     * Computes energy of current signal
     *-------------------------------------*/

    ener = 0.0f;
    for (i = 0; i < L_subfr; i++)
    {
        ener += ptr_sig_in[i] * ptr_sig_in[i];
    }

    if (ener < 0.1f)
    {
        *num_gltp = 0.0f;
        *den_gltp = 1.0f;
        *ltpdel = 0;
        *phase = 0;

        return;
    }

    /*-------------------------------------
     * Selects best of 3 integer delays
     * Maximum of 3 numerators around t0
     *-------------------------------------*/

    lambda = t0 - 1;
    ptr_sig_past = ptr_sig_in - lambda;
    num_int = -1.0e30f;
    i_max = 0;

    for (i = 0; i < 3; i++)
    {
        num = 0.0f;
        for (n = 0; n < L_subfr; n++)
        {
            num += ptr_sig_in[n] * ptr_sig_past[n];
        }
        if (num > num_int)
        {
            i_max = i;
            num_int = num;
        }
        ptr_sig_past--;
    }

    if (num_int <= 0.0f)
    {
        *num_gltp = 0.0f;
        *den_gltp = 1.0f;
        *ltpdel = 0;
        *phase = 0;

        return;
    }

    /* Calculates denominator for i_max */
    lambda += i_max;
    ptr_sig_past = ptr_sig_in - lambda;
    den_int = (float) 0.;
    for (n = 0; n < L_subfr; n++)
    {
        den_int += ptr_sig_past[n] * ptr_sig_past[n];
    }

    if (den_int < (float) 0.1)
    {
        *num_gltp = (float) 0.;
        *den_gltp = (float) 1.;
        *ltpdel = 0;
        *phase = 0;
        return;
    }

    /*----------------------------------
     * Select best phase around lambda
     * Compute y_up & denominators
     *----------------------------------*/

    ptr_y_up = y_up;
    den_max = den_int;
    ptr_den0 = tab_den0;
    ptr_den1 = tab_den1;
    ptr_h = tab_hup_s;
    ptr_sig_past0 = ptr_sig_in + LH_UP_S - 1 - lambda;        /* points on lambda_max+1 */

    /* loop on phase  */
    for (phi = 1; phi < F_UP_PST; phi++)
    {
        /* Computes criterion for (lambda+1) - phi/F_UP_PST     */
        /* and lambda - phi/F_UP_PST                            */
        ptr_sig_past = ptr_sig_past0;
        /* computes y_up[n] */
        for (n = 0; n <= L_subfr; n++)
        {
            ptr1 = ptr_sig_past++;
            temp0 = (float) 0.;
            for (i = 0; i < LH2_S; i++)
            {
                temp0 += ptr_h[i] * ptr1[-i];
            }
            ptr_y_up[n] = temp0;
        }

        /* compute den0 (lambda+1) and den1 (lambda) */
        /* part common to den0 and den1 */
        temp0 = (float) 0.;
        for (n = 1; n < L_subfr; n++)
        {
            temp0 += ptr_y_up[n] * ptr_y_up[n];
        }

        /* den0 */
        den0 = temp0 + ptr_y_up[0] * ptr_y_up[0];
        *ptr_den0++ = den0;

        /* den1 */
        den1 = temp0 + ptr_y_up[L_subfr] * ptr_y_up[L_subfr];
        *ptr_den1++ = den1;
        if (fabs (ptr_y_up[0]) > fabs (ptr_y_up[L_subfr]))
        {
            if (den0 > den_max)
            {
                den_max = den0;
            }
        }
        else
        {
            if (den1 > den_max)
            {
                den_max = den1;
            }
        }
        ptr_y_up += (L_subfr+1);
        ptr_h += LH2_S;
    }
    if (den_max < 0.1f)
    {
        *num_gltp = 0.0f;
        *den_gltp = 1.0f;
        *ltpdel = 0;
        *phase = 0;
        return;
    }

    /* Computation of the numerators                */
    /* and selection of best num*num/den            */
    /* for non null phases                          */

    /* Initialize with null phase */
    num_max = num_int;
    den_max = den_int;
    numsq_max = num_max * num_max;
    phi_max = 0;
    ioff = 1;

    ptr_den0 = tab_den0;
    ptr_den1 = tab_den1;
    ptr_y_up = y_up;

    /* if den_max = 0 : will be selected and declared unvoiced */
    /* if num!=0 & den=0 : will be selected and declared unvoiced */
    /* degenerated seldom cases, switch off LT is OK */

    /* Loop on phase */
    for (phi = 1; phi < F_UP_PST; phi++)
    {

        /* computes num for lambda+1 - phi/F_UP_PST */
        num = 0.0f;
        for (n = 0; n < L_subfr; n++)
        {
            num += ptr_sig_in[n] * ptr_y_up[n];
        }
        if (num < 0.0f)
        {
            num = 0.0f;
        }
        numsq = num * num;

        /* selection if num/sqrt(den0) max */
        den0 = *ptr_den0++;
        temp0 = numsq * den_max;
        temp1 = numsq_max * den0;
        if (temp0 > temp1)
        {
            num_max = num;
            numsq_max = numsq;
            den_max = den0;
            ioff = 0;
            phi_max = phi;
        }

        /* computes num for lambda_max - phi/F_UP_PST */
        ptr_y_up++;
        num = (float) 0.;
        for (n = 0; n < L_subfr; n++)
        {
            num += ptr_sig_in[n] * ptr_y_up[n];
        }
        if (num < (float) 0.)
        {
            num = (float) 0.;
        }
        numsq = num * num;

        /* selection if num/sqrt(den1) max */
        den1 = *ptr_den1++;
        temp0 = numsq * den_max;
        temp1 = numsq_max * den1;
        if (temp0 > temp1)
        {
            num_max = num;
            numsq_max = numsq;
            den_max = den1;
            ioff = 1;
            phi_max = phi;
        }
        ptr_y_up += L_subfr;
    }

    /*---------------------------------------------------
     * test if normalized crit0[iopt] > THRESHCRIT
     *--------------------------------------------------*/

    if ((num_max == 0.0f) || (den_max <= 0.1f))
    {
        *num_gltp = 0.0f;
        *den_gltp = 1.0f;
        *ltpdel = 0;
        *phase = 0;
        return;
    }

    /* comparison num * num             */
    /* with ener * den x THRESCRIT      */
    temp1 = den_max * ener * THRESCRIT;
    if (numsq_max >= temp1)
    {
        *ltpdel = lambda + 1 - ioff;
        *off_yup = ioff;
        *phase = phi_max;
        *num_gltp = num_max;
        *den_gltp = den_max;
    }
    else
    {
        *num_gltp = 0.0f;
        *den_gltp = 1.0f;
        *ltpdel = 0;
        *phase = 0;
    }

    return;
}

/*----------------------------------------------------------------------------
 * filt_plt()
 *
 * Perform long term postfilter
 *----------------------------------------------------------------------------*/

static void filt_plt(
    const float *s_in,      /* i  : input signal with past      */
    const float *s_ltp,     /* i  : filtered signal with gain 1 */
    float *s_out,     /* o  : output signal               */
    const float gain_plt,   /* i  : filter gain                 */
    const short L_subfr     /* i  : the length of subframe      */
)
{
    int n;
    float gain_plt_1;

    gain_plt_1 = (float) 1. - gain_plt;

    for (n = 0; n < L_subfr; n++)
    {
        s_out[n] = gain_plt * s_in[n] + gain_plt_1 * s_ltp[n];
    }

    return;
}

/*----------------------------------------------------------------------------
 * compute_ltp_l()
 *
 * compute delayed signal, num & den of gain for fractional delay
 * with long interpolation filter
 *----------------------------------------------------------------------------*/

static void compute_ltp_l(
    const float *s_in,    /* i  : input signal with past  */
    const   int ltpdel,   /* i  : delay factor            */
    const   int phase,    /* i  : phase factor            */
    float *y_up,    /* o  : delayed signal          */
    float *num,     /* o  : numerator of LTP gain   */
    float *den,     /* o  : denominator of LTP gain */
    const short L_subfr   /* i  : the length of subframe  */
)
{
    const float *ptr_h;
    int n, i;
    const float *ptr2;
    float temp;

    /* Filtering with long filter */
    ptr_h = tab_hup_l + (phase - 1) * LH2_L;
    ptr2 = s_in - ltpdel + LH_UP_L;

    /* Compute y_up */
    for (n = 0; n < L_subfr; n++)
    {
        temp = 0.0f;
        for (i = 0; i < LH2_L; i++)
        {
            temp += ptr_h[i] **ptr2--;
        }
        y_up[n] = temp;
        ptr2 += LH2_L_P1;
    }

    /* Compute num */
    *num = 0.0f;
    for (n = 0; n < L_subfr; n++)
    {
        *num += y_up[n] * s_in[n];
    }

    if (*num < 0.0f)
    {
        *num = 0.0f;
    }

    /* Compute den */
    *den = 0.0f;
    for (n = 0; n < L_subfr; n++)
    {
        *den += y_up[n] * y_up[n];
    }

    return;
}

/*----------------------------------------------------------------------------
 *  select_ltp()
 *
 *  selects best of (gain1, gain2)
 *  with gain1 = num1 * 2** sh_num1 / den1 * 2** sh_den1
 *  and  gain2 = num2 * 2** sh_num2 / den2 * 2** sh_den2
 *----------------------------------------------------------------------------*/

static int select_ltp( /* o  : 1 = 1st gain, 2 = 2nd gain */
    const float num1,  /* i  : numerator of gain1   */
    const float den1,  /* i  : denominator of gain1 */
    const float num2,  /* i  : numerator of gain2   */
    const float den2   /* i  : denominator of gain2 */
)
{
    if (den2 == (float) 0.)
    {
        return (1);
    }

    if (num2 * num2 * den1 > num1 * num1 * den2)
    {
        return (2);
    }
    else
    {
        return (1);
    }
}


/*------------------------------------------------------------------------------------
 * modify_pst_param()
 *
 * Modify gamma1 and gamma2 values in function of the long-term noise level
 *-----------------------------------------------------------------------------------*/

static void modify_pst_param(
    const float psf_lp_noise, /* i  : Long term noise energy                */
    float *g1,          /* o  : Gamma1 used in post filter            */
    float *g2,          /* o  : Gamma2 used in post filter            */
    const short coder_type,   /* i  : coder type                            */
    float *gain_factor  /* o  : Gain factor applied in post filtering */
)
{
    float ftmp;

    if( coder_type != INACTIVE && psf_lp_noise < LP_NOISE_THR )
    {
        ftmp = psf_lp_noise*BG1 + CG1;
        if( ftmp > POST_G1 )
        {
            ftmp  = POST_G1;
        }
        else if( ftmp < POST_G1_MIN )
        {
            ftmp  = POST_G1_MIN;
        }
        *g1 = ftmp;

        ftmp = psf_lp_noise*BG2 + CG2;
        if( ftmp > POST_G2 )
        {
            ftmp  = POST_G2;
        }
        else if( ftmp < POST_G2_MIN )
        {
            ftmp  = POST_G2_MIN;
        }
        *g2 = ftmp;
    }
    else
    {
        *g1 = POST_G1_NOIS;
        *g2 = POST_G2_NOIS;
    }

    /* Set gain_factor of the harmonic filtering */
    ftmp = (psf_lp_noise - K_LP_NOISE) * C_LP_NOISE;

    if( ftmp >= 0.25f )
    {
        /* the noise is really high */
        *gain_factor = 0.25f;
    }
    else if ( ftmp < 0 )
    {
        *gain_factor = 0.0f;
    }
    else
    {
        *gain_factor = ftmp;
    }

    return;
}
