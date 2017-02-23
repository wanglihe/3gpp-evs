/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <assert.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"

/*---------------------------------------------------------------------*
 * Local constants
 *---------------------------------------------------------------------*/

#define A2       0.2f
#define GAIN_VAR 0.000011f

/*-------------------------------------------------------*
 * CNG_exc()
 *
 * Comfort noise generation routine
 *-------------------------------------------------------*/

void CNG_exc(
    const long  core_brate,       /* i  : core bitrate                          */
    const short L_frame,          /* i  : length of the frame                   */
    float *Enew,            /* i/o: decoded SID energy                    */
    short *seed,            /* i/o: random generator seed                 */
    float exc[],            /* o  : current non-enhanced excitation       */
    float exc2[],           /* o  : current enhanced excitation           */
    float *lp_ener,         /* i/o: LP filtered E                         */
    const long  last_core_brate,  /* i  : previous frame core bitrate           */
    short *first_CNG,       /* i/o: first CNG frame flag for energy init. */
    short *cng_ener_seed,   /* i/o: random generator seed for CNG energy  */
    float bwe_exc[],        /* o  : excitation for SWB TBE                */
    const short allow_cn_step,    /* i  : allow CN step                         */
    short *last_allow_cn_step,  /* i/o: last allow step                   */
    const short num_ho,           /* i  : number of selected hangover frames    */
    float q_env[],
    float *lp_env,
    float *old_env,
    float *exc_mem,
    float *exc_mem1,
    short *sid_bw,
    short *cng_ener_seed1,
    float exc3[],
    short Opt_AMR_WB
)
{
    float enr;
    short i;
    float ener_lp;
    short i_subfr;
    short pit_max;
    float ftmp;
    float *ptR,*ptI;
    float fft_io[L_FRAME16k];
    float itmp[129];
    float env[NUM_ENV_CNG];
    float enr1;
    float denv[NUM_ENV_CNG];

    /*------------------------------------------------------------------*
     * Initializations
     *------------------------------------------------------------------*/

    if( L_frame == L_FRAME )
    {
        pit_max = PIT_MAX;
    }
    else /* L_frame == L_FRAME16k */
    {
        pit_max = PIT16k_MAX;
    }

    /*---------------------------------------------------------------------*
     * Initialization of CNG energy for the first CNG frame
     *---------------------------------------------------------------------*/

    if( *first_CNG == 0 )
    {
        if( core_brate == FRAME_NO_DATA )
        {
            /* needed only in decoder when the very first SID frame was erased and this frame is FRAME_NO_DATA frame */
            *Enew = dotp( exc-pit_max, exc-pit_max, pit_max ) / pit_max;
        }

        *lp_ener = *Enew;
    }

    /*---------------------------------------------------------------------*
     * Update CNG energy
     *---------------------------------------------------------------------*/

    if( last_core_brate != SID_1k75 && last_core_brate != FRAME_NO_DATA && last_core_brate != SID_2k40 )
    {
        /* Partially reset CNG energy after active speech period */
        if ( allow_cn_step == 0 && *last_allow_cn_step == 0 )
        {
            if( num_ho < 3 || *Enew < 1.5f * *lp_ener )
            {
                *lp_ener = 0.8f * *lp_ener + 0.2f * *Enew;
            }
            else
            {
                *lp_ener = 0.95f * *lp_ener + 0.05f * *Enew;
            }
        }
        else
        {
            *lp_ener = *Enew;
            *last_allow_cn_step = 0;
        }
    }
    else
    {
        /* normal CNG update */
        if ( *last_allow_cn_step == 0 )
        {
            *lp_ener = (float)(A2 **Enew + (1-A2) **lp_ener);
        }
        else
        {
            if ( core_brate == SID_1k75 || core_brate == SID_2k40 )
            {
                *last_allow_cn_step = 0;
            }

            *lp_ener = *Enew;

        }
    }

    if ( allow_cn_step == 1 )
    {
        *last_allow_cn_step = 1;
    }

    /*---------------------------------------------------------------------*
     * Generate white noise vector
     *---------------------------------------------------------------------*/

    for ( i=0; i<L_frame; i++ )
    {
        exc2[i] = (float)own_random( seed );
    }

    /*------------------------------------------------------------*
     * Insert random variation for excitation energy
     *  (random variation is scaled according to *lp_ener value)
     *------------------------------------------------------------*/

    for ( i_subfr=0; i_subfr<L_frame; i_subfr += L_SUBFR )
    {
        ener_lp = ( (own_random( cng_ener_seed ) * (*lp_ener) ) * GAIN_VAR ) + (*lp_ener);

        if( ener_lp < 0.0f)
        {
            ener_lp = 0.01f;
        }
        enr = dotp( &exc2[i_subfr], &exc2[i_subfr], L_SUBFR ) + 0.01f;

        enr = (float)sqrt( ener_lp * L_SUBFR / enr );

        for ( i=0; i<L_SUBFR; i++ )
        {
            exc2[i_subfr+i] *= enr;
        }
    }

    if ( Opt_AMR_WB != 1 )
    {
        mvr2r( exc2, exc3, L_FRAME16k);

        enr1 = (float)log10( *Enew*L_frame + 0.1f ) / (float)log10( 2.0f );

        if ( core_brate == SID_2k40 )
        {
            if ( *sid_bw == 0 )
            {
                for ( i=0; i<NUM_ENV_CNG; i++ )
                {
                    /* get quantized envelope */
                    env[i] = pow(2.0f,(enr1 - q_env[i]));
                }
            }

            /* initialize CNG envelope */
            if( *first_CNG == 0 && *sid_bw == 0 )
            {
                mvr2r(env, lp_env, NUM_ENV_CNG);
            }

            if ( *sid_bw == 0 )
            {
                mvr2r(env, old_env, NUM_ENV_CNG);
            }
        }

        for ( i=0; i<NUM_ENV_CNG; i++ )
        {
            /* get AR low-passed envelope */
            lp_env[i] = 0.9f*lp_env[i] + (1-0.9f)*old_env[i];
        }

        /* calculate the spectrum of random excitation signal */
        mvr2r(exc2, fft_io, L_frame);

        if ( L_frame == L_FRAME16k )
        {
            modify_Fs( fft_io, L_FRAME16k, 16000, fft_io, 12800, exc_mem1, 0 );
        }

        fft_rel(fft_io, L_FFT, LOG2_L_FFT);
        ptR = &fft_io[1];
        ptI = &fft_io[L_FFT-1];
        for ( i=0; i<NUM_ENV_CNG; i++ )
        {
            env[i] = 2.0f*(*ptR **ptR + *ptI **ptI)/L_FFT;
            ptR++;
            ptI--;
        }

        for ( i=0; i<NUM_ENV_CNG; i++ )
        {
            denv[i] = lp_env[i] + 2*(*lp_ener) - env[i];

            if ( denv[i] < 0 )
            {
                denv[i] = 0;
            }
        }
        set_f(itmp, 0.0f, NUM_ENV_CNG);

        set_f(fft_io, 0.0f, L_FFT);
        ptR = &fft_io[1];
        ptI = &fft_io[L_FFT-1];
        for (i=0; i<NUM_ENV_CNG; i++)
        {
            *ptR = own_random( cng_ener_seed1 );
            *ptI = own_random( cng_ener_seed1 );

            env[i] = 2.0f*(*ptR **ptR + *ptI **ptI)/L_FFT;
            ptR++;
            ptI--;
        }

        for ( i=0; i<NUM_ENV_CNG; i++ )
        {
            itmp[i] += own_random( cng_ener_seed1 )*denv[i]*0.000011f + denv[i];

            if (itmp[i] < 0.0f)
            {
                itmp[i] = 0;
            }
        }

        ptR = &fft_io[1];
        ptI = &fft_io[L_FFT-1];
        for ( i=0; i<NUM_ENV_CNG; i++ )
        {
            *ptR *= sqrt(itmp[i]/env[i]);
            *ptI *= sqrt(itmp[i]/env[i]);

            ptR++;
            ptI--;
        }

        ifft_rel(fft_io, L_FFT, LOG2_L_FFT);

        if ( L_frame == L_FRAME16k )
        {
            modify_Fs( fft_io, L_FFT, 12800, fft_io, 16000, exc_mem, 0 );
        }

        enr1 = dotp( fft_io, fft_io, L_frame ) / L_frame;

        /* add time domain randomization */
        for ( i_subfr=0; i_subfr<L_frame; i_subfr += L_SUBFR )
        {
            enr = dotp( &fft_io[i_subfr], &fft_io[i_subfr], L_SUBFR ) + 0.01f;
            ener_lp = ( (own_random( cng_ener_seed1 ) * (enr1) ) * 0.000011f ) + (enr1);
            ener_lp *= L_SUBFR;
            enr = (float)sqrt( ener_lp / enr );

            if( last_core_brate != SID_2k40 && last_core_brate != SID_1k75 && last_core_brate != FRAME_NO_DATA && core_brate == SID_2k40 )
            {
                if ( enr > 1 )
                {
                    enr = 1;
                }
            }

            for ( i=0; i<L_SUBFR; i++ )
            {
                fft_io[i_subfr+i] *= enr;
            }
        }

        for ( i=0; i<L_frame; i++ )
        {
            fft_io[i] = 0.75f*fft_io[i] + exc2[i];
        }

        enr = (dotp( fft_io, fft_io, L_frame ) / L_frame) + 0.01f;
        enr = (*lp_ener)/enr;

        if ( enr > 1 )
        {
            enr = 1;
        }

        ftmp = sqrt(enr);
        for (i=0; i<L_frame; i++)
        {
            fft_io[i] *= ftmp;
        }

        mvr2r( fft_io, exc2, L_frame );
    }
    if ( Opt_AMR_WB != 1 )
    {
        mvr2r( exc3, exc, L_frame );
    }
    else
    {
        mvr2r( exc2, exc, L_frame );
    }

    if( L_frame == L_FRAME )
    {
        interp_code_5over2( exc2, bwe_exc, L_frame );
    }
    else
    {
        interp_code_4over2( exc2, bwe_exc, L_frame );
    }

    return;
}


/*-------------------------------------------------------*
 * cng_params_postupd
 *
 * Post-update of CNG parameters
 *-------------------------------------------------------*/
void cng_params_postupd(
    const short ho_circ_ptr,           /* i  : pointer for CNG averaging buffers   */
    short *cng_buf_cnt,          /* i/o: counter for CNG store buffers       */
    const float *const cng_exc2_buf,   /* i  : Excitation buffer                   */
    const long  *const cng_brate_buf,  /* i  : bit rate buffer                     */
    float ho_env_circ[]          /* i/o: Envelope buffer                     */
)
{
    short i, j;
    const float *exc2;
    float fft_io[L_FFT];
    float sp[129];
    float *ptR,*ptI;
    float env[NUM_ENV_CNG];
    short CNG_mode;
    short ptr;
    float att;
    long  last_active_brate;

    ptr = ho_circ_ptr - *cng_buf_cnt + 1;
    if( ptr < 0 )
    {
        ptr += HO_HIST_SIZE;
    }

    for( j = 0; j < *cng_buf_cnt; j++ )
    {
        exc2  = &cng_exc2_buf[ptr*L_FFT];
        last_active_brate = cng_brate_buf[ptr];

        /* calculate the spectrum of residual signal */
        /* calculate the spectrum of residual signal */
        mvr2r(exc2, fft_io, L_FFT);

        fft_rel(fft_io, L_FFT, LOG2_L_FFT);

        ptR = &fft_io[1];
        ptI = &fft_io[L_FFT-1];
        for (i=0; i<NUM_ENV_CNG; i++)
        {
            sp[i] = 2.0f*(*ptR **ptR + *ptI **ptI)/L_FFT;
            ptR++;
            ptI--;
        }

        mvr2r(sp,env,NUM_ENV_CNG);
        if( last_active_brate > ACELP_13k20 )
        {
            CNG_mode = 4;
        }
        else if( last_active_brate > ACELP_9k60 )
        {
            CNG_mode = 3;
        }
        else if( last_active_brate > ACELP_8k00 )
        {
            CNG_mode = 2;
        }
        else if( last_active_brate > ACELP_7k20 )
        {
            CNG_mode = 1;
        }
        else
        {
            CNG_mode = 0;
        }

        att = 1/pow(2,ENR_ATT[CNG_mode]);

        for ( i=0; i<NUM_ENV_CNG; i++ )
        {
            env[i] *= att;
        }

        /* update the circular buffer of old residual envelope */
        mvr2r( env, &(ho_env_circ[(ho_circ_ptr)*NUM_ENV_CNG]), NUM_ENV_CNG );

        ptr++;
        if(ptr == HO_HIST_SIZE)
        {
            ptr = 0;
        }
    }

    *cng_buf_cnt = 0;

    return;

}



/*-------------------------------------------------------*
 * cng_params_upd()
 *
 * update CNG parameters
 *-------------------------------------------------------*/

void cng_params_upd(
    const float lsp_new[],        /* i  : LSP parameters                                      */
    const float exc2[],           /* i  : current enhanced excitation                         */
    const short L_frame,          /* i  : frame length                                        */
    short *ho_circ_ptr,     /* i/o: pointer for CNG averaging buffers                   */
    float ho_ener_circ[],   /* o  : energy buffer for CNG averaging                     */
    short *ho_circ_size,    /* i/o: size of DTX hangover history buffer for averaging   */
    float ho_lsp_circ[],    /* o  : old LSP buffer for CNG averaging                    */
    const short enc_dec_flag,     /* i  : Flag indicating encoder or decoder (ENC,DEC)    */
    float ho_env_circ[],    /* i/o: Envelope buffer                                 */
    short *cng_buf_cnt,     /* i/o: Counter of postponed FFT-processing instances   */
    float cng_exc2_buf[],   /* i/o: Excitation buffer                               */
    long  cng_brate_buf[],  /* i/o: last_active_brate buffer                        */
    const long  last_active_brate /* i  : Last active bit rate                            */
)
{
    float enr;
    float fft_io[L_FRAME16k];
    float sp[129];
    float *ptR,*ptI;
    float env[NUM_ENV_CNG];
    short i;
    short CNG_mode;
    float att;

    /* update the pointer to circular buffer of old LSP vectors */
    (*ho_circ_ptr)++;
    if( *ho_circ_ptr == HO_HIST_SIZE )
    {
        *ho_circ_ptr = 0;
    }

    /* update the circular buffer of old LSP vectors with the new LSP vector */
    mvr2r( lsp_new, &(ho_lsp_circ[(*ho_circ_ptr)*M]), M );

    /* calculate the residual signal energy */
    enr = dotp( exc2, exc2, L_frame ) / L_frame;

    /* update the circular buffer of old energies */
    ho_ener_circ[*ho_circ_ptr] = enr;

    if( enc_dec_flag == ENC)
    {
        /* Store residual signal for postponed FFT-processing*/
        (*cng_buf_cnt)++;
        if( *cng_buf_cnt > HO_HIST_SIZE)
        {
            *cng_buf_cnt = HO_HIST_SIZE;
        }
        mvr2r( exc2, &(cng_exc2_buf[(*ho_circ_ptr)*L_FFT]), L_FFT );
        cng_brate_buf[*ho_circ_ptr] = last_active_brate;
    }
    else
    {

        /* calculate the spectrum of residual signal */
        mvr2r(exc2, fft_io, L_frame);

        fft_rel(fft_io, L_FFT, LOG2_L_FFT);

        ptR = &fft_io[1];
        ptI = &fft_io[L_FFT-1];
        for (i=0; i<NUM_ENV_CNG; i++)
        {
            sp[i] = 2.0f*(*ptR **ptR + *ptI **ptI)/L_FFT;
            ptR++;
            ptI--;
        }

        mvr2r(sp,env,NUM_ENV_CNG);
        if( last_active_brate > ACELP_13k20 )
        {
            CNG_mode = 4;
        }
        else if( last_active_brate > ACELP_9k60 )
        {
            CNG_mode = 3;
        }
        else if( last_active_brate > ACELP_8k00 )
        {
            CNG_mode = 2;
        }
        else if( last_active_brate > ACELP_7k20 )
        {
            CNG_mode = 1;
        }
        else
        {
            CNG_mode = 0;
        }

        att = 1/pow(2,ENR_ATT[CNG_mode]);

        for ( i=0; i<NUM_ENV_CNG; i++ )
        {
            env[i] *= att;
        }

        /* update the circular buffer of old residual envelope */
        mvr2r( env, &(ho_env_circ[(*ho_circ_ptr)*NUM_ENV_CNG]), NUM_ENV_CNG );

    }

    (*ho_circ_size)++;
    if( *ho_circ_size > HO_HIST_SIZE )
    {
        *ho_circ_size = HO_HIST_SIZE;
    }

    return;
}
