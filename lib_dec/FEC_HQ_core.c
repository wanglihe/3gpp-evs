/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_dec.h"
#include "rom_com.h"
#include "prot.h"

/*---------------------------------------------------------------------*
 * Local prototypes
 *---------------------------------------------------------------------*/

static short FEC_phase_matching( Decoder_State *st, float *ImdctOut, float *auOut, float *OldauOut, float OldauOut_pha[2][N_LEAD_NB] );

static void FEC_phase_matching_nextgood( const float *ImdctOut, float *auOut, float *OldauOut,
        float OldauOut_pha[2][N_LEAD_NB], float mean_en_high );

static void FEC_phase_matching_burst( const float *ImdctOut, float *auOut, float *OldauOut,
                                      float OldauOut_pha[2][N_LEAD_NB], float *prev_oldauOut );

static void Repetition_smoothing_nextgood( const float *ImdctOut, float *auOut, float *OldImdctOut, float *OldauOut,
        short cur_data_use_flag, short overlap_time );

static int Repetition_smoothing( const float *ImdctOut, float *auOut, float *OldImdctOut, float *OldauOut,
                                 const short L, float *prev_oldauOut, short overlap_time );

static void Windowing_1st_NB( float *ImdctOutWin, const float *ImdctOut, const float *win, const float *smoothingWin, short smoothing_flag );

static void Windowing_2nd_NB( float *ImdctOutWin, const float *ImdctOut, const float *win );

static void common_overlapping( float *auOut, float *ImdctOutWin, float *OldauOut, short end1, short offset1,
                                short start2, short end2, short offset_i2, short offset2 );

static void Smoothing_vector_NB( const float OldauOutnoWin[], const float ImdctOutWin[], const float SmoothingWin[],
                                 float auOut[], const short ol_size );

static void Smoothing_vector_scaledown_NB( const float OldauOutnoWin[], const float ImdctOutWin[], const float SmoothingWin[],
        float auOut[], const short ol_size );

static void Scaledown( float x[], float y[], float scale_v, const short N );


/*--------------------------------------------------------------------------*
 * Regression_Anal()
 *
 *
 *--------------------------------------------------------------------------*/

static void Regression_Anal(
    const float *values,    /* i : Previous values                  */
    float *r_p,       /* o : Output r[a b] array : y=ax+b     */
    const short num_pgf     /* i : Number of previous good frame    */
)
{
    short i;
    float aindex[MAX_PGF+1], b_p[MAX_PGF+1];

    /* Initialize */
    for(i=0; i<num_pgf+1; i++)
    {
        aindex[i]=0.f;
        b_p[i]=0.f;
    }

    /* [aindex[0] aindex[1]][r[0]]=[b[0]]*/
    /* [aindex[1] aindex[2]][r[1]] [b[1]]*/
    /* r[0] is the y-intercept(initial value). r[1] is slope*/

    /* r[0] = (b[0]a[2]-a[1]b[1])/(a[0]a[2]-a[1]a[1])
       r[1] = (b[0]a[1]-a[0]b[1])/(a[1]a[1]-a[0]a[2]) */

    aindex[0]=num_pgf;
    for(i=1; i<num_pgf+1; i++)
    {
        aindex[1] += (float)i;
        aindex[2] += ((float)i*(float)i);
    }
    /* Calculate b[] */
    for(i=0; i<num_pgf; i++)
    {
        b_p[0] += (values[i]);
        b_p[1] += ((float)(num_pgf-i)*values[i]);
    }
    r_p[0] = (b_p[0]*aindex[2] - aindex[1]*b_p[1]) / (aindex[0]*aindex[2] - aindex[1]*aindex[1]);
    r_p[1] = (b_p[0]*aindex[1] - aindex[0]*b_p[1]) / (aindex[1]*aindex[1] - aindex[0]*aindex[2]);

    return;
}


/*--------------------------------------------------------------------------*
 * FEC_scaling()
 *
 *
 *--------------------------------------------------------------------------*/
static void FEC_scaling(
    float *old_coeffs,            /* i/o  : Pointer to old MDCT coeffs.                    */
    float *t_audio_q,             /* o  : MDCT coeffs. (for synthesis)                     */
    float *Norm_gain,             /* i  : Gain for Norm of each band                       */
    short *HQ_FEC_seed,           /* i/o  : Seed for Ransom number Generator               */
    short nb_sfm,                 /* i  : Number of sub-band                               */
    const short *start_band,
    const short *end_band
)
{
    short i, j;

    for (i = 0; i < RANDOM_START; i++)
    {
        for (j = start_band[i]; j < end_band[i]; j++)
        {
            t_audio_q[j] = Norm_gain[i] * old_coeffs[j];
        }
    }

    for (i = RANDOM_START; i < nb_sfm; i++)
    {
        for (j = start_band[i]; j < end_band[i]; j++)
        {
            if((float)own_random( HQ_FEC_seed )<0.f)
            {
                t_audio_q[j] = Norm_gain[i] * (-old_coeffs[j]);
            }
            else
            {
                t_audio_q[j] = Norm_gain[i] * old_coeffs[j];
            }
        }
    }
    return;
}


/*--------------------------------------------------------------------------*
 * HQ_FEC_processing()
 *
 *
 *--------------------------------------------------------------------------*/

void HQ_FEC_processing(
    Decoder_State *st,            /* i/o: decoder state structure                          */
    float *t_audio_q,             /* o  : MDCT coeffs. (for synthesis)                     */
    short is_transient,           /* i  : Old flag for transient                           */
    float ynrm_values[][MAX_PGF], /* i  : Old average Norm values for each group of bands  */
    float r_p_values[][MAX_ROW],  /* i  : Computed y-intercept and slope by Regression     */
    short num_Sb,                 /* i  : Number of sub-band group                         */
    short nb_sfm,                 /* i  : Number of sub-band                               */
    short *Num_bands_p,           /* i  : Number of coeffs. for each sub-band              */
    short output_frame,           /* i  : Frame size                                       */
    const short *sfm_start,             /* i  : Start of bands                                   */
    const short *sfm_end                /* i  : End of bands                                     */
)
{
    short i, j, k;
    float energy_diff;
    short mute_start,num_pgf;
    short Num_sb_bwe;
    float tmp, norm_p[MAX_SB_NB];
    float *norm_values;
    float *r_p;
    short sfm;

    /* Make sure 'energy_MA_Curr[0]' is not zero to prevent a 'div by 0' error. */
    if (st->energy_MA_Curr[0] < 1.0f)
    {
        st->energy_MA_Curr[0] = 1.0f; /*It seems to be in Q0 the FxP*/
    }
    /* Decide the start frame number for adaptive muting */
    /* Normalized energy difference between the current frame and the moving average */
    energy_diff = (float)fabs((st->energy_MA_Curr[1] - st->energy_MA_Curr[0])/st->energy_MA_Curr[0]);

    if((energy_diff<ED_THRES) && (is_transient==0)) /* First erasure frame */
    {
        mute_start = 5;
    }
    else
    {
        mute_start = 2;
    }

    if( st->prev_old_bfi == 1 && st->nbLostCmpt == 1 && output_frame == L_FRAME8k )
    {
        st->nbLostCmpt++;
    }

    /* Frequency-domain FEC */
    if( st->nbLostCmpt == 1 ) /* First erasure frame */
    {
        if ( is_transient == 0 )
        {
            if(energy_diff < ED_THRES)
            {
                for (i=0; i < output_frame; i++)
                {
                    t_audio_q[i] = st->old_coeffs[i];
                }
            }
            else
            {
                for (i=0; i < output_frame; i++)
                {
                    st->old_coeffs[i] *= SCALE_DOWN_3dB;
                    t_audio_q[i] = st->old_coeffs[i];
                }
            }

            /* Sign prediction in 4-dim bands up to 1.6 kHz*/
            if (st->old_is_transient[1] == 0)
            {
                if (st->old_is_transient[2] == 0)
                {
                    for (sfm = 0; sfm < HQ_FEC_SIGN_SFM; sfm++)
                    {
                        if (st->prev_sign_switch[sfm] >= HQ_FEC_SIGN_THRES)
                        {
                            for(i = 0; i < HQ_FEC_BAND_SIZE; i++)
                            {
                                t_audio_q[i+sfm*HQ_FEC_BAND_SIZE]*= -1.0f;
                            }
                        }
                    }
                }
                else
                {
                    for (sfm = 0; sfm < HQ_FEC_SIGN_SFM; sfm++)
                    {
                        if (st->prev_sign_switch[sfm] >= HQ_FEC_SIGN_THRES_TRANS)
                        {
                            for(i = 0; i < HQ_FEC_BAND_SIZE; i++)
                            {
                                t_audio_q[i+sfm*HQ_FEC_BAND_SIZE]*= -1.0f;
                            }
                        }
                    }
                }
            }
            else
            {
                for (i = RANDOM_START*8; i < output_frame; i++)
                {
                    if( (float)own_random( &st->HQ_FEC_seed ) < 0.0f )
                    {
                        t_audio_q[i] *= -1.0f;
                    }
                }
            }
        }
        else
        {
            if( st->old_is_transient[1] )      /* hangover */
            {
                for (i=0; i < output_frame; i++)
                {
                    st->old_coeffs[i] *= SCALE_DOWN_3dB;
                    t_audio_q[i] = st->old_coeffs[i];
                }
            }
            else
            {
                for (i = 0; i < RANDOM_START*8; i++)
                {
                    st->old_coeffs[i] *= SCALE_DOWN_3dB;
                    t_audio_q[i] = st->old_coeffs[i];
                }

                for (i = RANDOM_START*8; i < output_frame; i++)
                {
                    st->old_coeffs[i] *= SCALE_DOWN_3dB;
                    if( (float)own_random( &st->HQ_FEC_seed ) < 0.0f )
                    {
                        t_audio_q[i] = (-st->old_coeffs[i]);
                    }
                    else
                    {
                        t_audio_q[i] = st->old_coeffs[i];
                    }
                }
            }
        }
    }
    else /* st->nbLostCmpt > 1 */
    {
        if( energy_diff < ED_THRES && is_transient == 0 )
        {
            num_pgf = 4;
        }
        else
        {
            num_pgf = 2;
        }

        Num_sb_bwe = num_Sb;
        if( st->nbLostCmpt == 2 )
        {
            for( i=0; i<Num_sb_bwe; i++ )
            {
                norm_values = &ynrm_values[i][0];
                r_p = &r_p_values[i][0];
                Regression_Anal(norm_values, r_p, num_pgf);
            }
        }

        /* Fade-out Norm by the result of Regression */
        if( st->nbLostCmpt >= mute_start )
        {
            /* Scaling */
            for ( i=0; i < output_frame; i++ )
            {
                st->old_coeffs[i] *= SCALE_DOWN_3dB;
            }
        }

        k=0;
        for( i=0; i<Num_sb_bwe; i++)
        {
            norm_values = &ynrm_values[i][0];
            r_p =  &r_p_values[i][0];

            /* Predict the average energy of each sub-band using Regression */
            /* Linear Regression */
            if( r_p[1] > MAX_TILT )
            {
                r_p[1] = MAX_TILT;
                norm_p[i] = (float)(norm_values[0]+r_p[1]*(float)(st->nbLostCmpt-1));
            }
            else
            {
                norm_p[i] = (float)(r_p[0]+r_p[1]*(float)(st->nbLostCmpt-1+num_pgf));
            }

            if( norm_values[0] != 0.0f && norm_p[i] > 0.0f ) /* Avoid negative value of the predicted energy */
            {
                tmp = norm_p[i]/norm_values[0]; /* Pred_new / Old */

                if( tmp > 1.0f )
                {
                    tmp = 1.f;
                }

                for( j=0; j<Num_bands_p[i]; j++ )
                {
                    st->Norm_gain[k++] = tmp;
                }
            }
            else
            {
                /* Scale down the last gain with the fixed gain(-3dB) */
                for(j=0; j<Num_bands_p[i]; j++)
                {
                    st->Norm_gain[k++] *= SCALE_DOWN_3dB;
                }
            }
        }

        /* Scaling for core band */
        FEC_scaling( st->old_coeffs, t_audio_q, st->Norm_gain, &st->HQ_FEC_seed, nb_sfm, sfm_start, sfm_end );
    }

    return;
}


/*--------------------------------------------------------------------------*
 * HQ_FEC_Mem_update()
 *
 *
 *--------------------------------------------------------------------------*/

void HQ_FEC_Mem_update(
    Decoder_State *st,               /* i/o: decoder state structure            */
    float *t_audio_q,
    float *normq,
    short *ynrm,
    short *Num_bands_p,
    short is_transient,
    short hqswb_clas,
    short c_switching_flag,
    short nb_sfm,
    short num_Sb,
    float *mean_en_high,
    short hq_core_type,              /* i : normal or low-rate MDCT(HQ) core */
    short output_frame
)
{
    short i, j, k;
    float tmp;
    float tmp_energy = 0;
    float *norm_values;
    short offset;
    short Min_ind;
    int Min_value;
    float Max_coeff;
    short Max_ind;
    short stat_mode_curr;
    float en_high[MAX_SB_NB];

    if (is_transient)
    {
        set_s(st->prev_sign_switch_2, 0, HQ_FEC_SIGN_SFM);
        set_s(st->prev_sign_switch, 0, HQ_FEC_SIGN_SFM);
    }
    else
    {
        for (j = 0; j < HQ_FEC_SIGN_SFM; j++)
        {
            st->prev_sign_switch[j] = st->prev_sign_switch_2[j];
            st->prev_sign_switch_2[j] = 0;

            for(i = 0; i < HQ_FEC_BAND_SIZE; i++)
            {
                tmp = st->old_coeffs[i+j*HQ_FEC_BAND_SIZE]*t_audio_q[i+j*HQ_FEC_BAND_SIZE];
                if (tmp < 0)
                {
                    st->prev_sign_switch[j]++;
                    st->prev_sign_switch_2[j]++;
                }
            }
        }
    }

    if( output_frame == L_FRAME8k )
    {
        /* if LR MDCT core is used, recalculate norms from decoded MDCT spectrum (using code from hq_hr_enc()) */
        if( ( hqswb_clas == HQ_HVQ ) || ( hq_core_type == LOW_RATE_HQ_CORE ) )
        {
            /* First group */
            logqnorm(t_audio_q, ynrm, 32, WID_G1, thren);
            j = ynrm[0];
            offset = WID_G1;

            for( i=1; i<SFM_G1; i++ )
            {
                logqnorm( &t_audio_q[offset], &ynrm[i], 40, WID_G1, thren );
                offset += WID_G1;
            }

            /* Second group */
            for( i=SFM_G1; i<SFM_G1+2; i++ )
            {
                logqnorm( &t_audio_q[offset], &ynrm[i], 40, WID_G2 , thren);
                offset += WID_G2;
            }

        }

        /* Memory update for the LGF log2 Norm */
        for (i = 0; i < nb_sfm; i++)
        {
            normq[i] = dicn[ynrm[i]];
        }

        k=0;
        for(i=0; i<num_Sb; i++)
        {
            norm_values = &st->ynrm_values[i][0];
            mvr2r( norm_values, &norm_values[1], MAX_PGF-1 );

            tmp = 0.f;
            for(j=0; j<Num_bands_p[i]; j++)
            {
                tmp += (float)normq[k++];
            }
            norm_values[0]=tmp/(float)Num_bands_p[i];
            tmp_energy += tmp;
        }

        if((c_switching_flag)||((st->last_core == ACELP_CORE)&&(st->core == HQ_CORE)))
        {
            for(i=0; i<MAX_SB_NB; i++)
            {
                for(j=1; j<MAX_PGF; j++)
                {
                    st->ynrm_values[i][j]=st->ynrm_values[i][0];
                }
            }
        }
        set_f(st->Norm_gain, 1.f, SFM_N_NB);

        /* st->energy_MA_Curr[1]=Energy of the current frame */
        st->energy_MA_Curr[1] = tmp_energy/(float)nb_sfm;

        /* Moving Average */
        st->energy_MA_Curr[0] = 0.8f*st->energy_MA_Curr[0]+0.2f*st->energy_MA_Curr[1];

        st->diff_energy = (float)fabs((st->energy_MA_Curr[1] - st->energy_MA_Curr[0])/st->energy_MA_Curr[0]);
        /* Classify the stationary mode : 12%  */
        if((st->diff_energy<ED_THRES_12P))
        {
            stat_mode_curr = 1;
        }
        else
        {
            stat_mode_curr = 0;
        }

        /* Apply Hysteresis to prevent frequent mode changing */
        if(st->stat_mode_old == stat_mode_curr)
        {
            st->stat_mode_out = stat_mode_curr;
        }

        st->stat_mode_old = stat_mode_curr;


        /* Find max. band index (Minimum value means maximum energy) */
        Min_ind=0;
        Min_value = 100;
        for( i=0; i<num_Sb; i++ )
        {
            if( Min_value>ynrm[i] )
            {
                Min_value = ynrm[i];
                Min_ind = i;
            }
        }

        /* Find max. coeff in band 0 */
        Max_ind = 0;
        if(Min_ind==0)
        {
            Max_coeff = 0.f;
            for(i=0; i<8; i++)
            {
                tmp = (float)fabs(t_audio_q[i]);
                if(Max_coeff<tmp)
                {
                    Max_coeff = tmp;
                    Max_ind = i;
                }
            }
        }

        /* Find energy difference from band 16 */
        k=1;

        for(i=k; i<num_Sb; i++)
        {
            en_high[i] = 0.f;
            for(j=0; j<2; j++)
            {
                en_high[i] += 0.5f*st->ynrm_values[i][j+1];
            }
        }

        *mean_en_high = 0.f;
        for(i=k; i<num_Sb; i++)
        {
            *mean_en_high += (float)(en_high[i]/st->ynrm_values[i][0]);
        }
        *mean_en_high /= (float)(num_Sb-k);

        if ((Min_ind<5) && ( abs(Min_ind - st->old_Min_ind)< 2) &&(st->diff_energy<ED_THRES_90P)&&(!st->bfi) && (!st->prev_bfi)&&(!st->prev_old_bfi)
                &&(!is_transient)&&(!st->old_is_transient[1]) && (st->prev_last_core==HQ_CORE) && (st->last_core==HQ_CORE))
        {
            if((Min_ind==0)&&(Max_ind<3))
            {
                st->phase_mat_flag = 0;
            }
            else
            {
                st->phase_mat_flag = 1;
            }
        }
        else
        {
            st->phase_mat_flag = 0;
        }

        st->old_Min_ind = Min_ind;
    }

    for (i=0; i < L_FRAME8k; i++)
    {
        st->old_coeffs[i] = t_audio_q[i];
    }

    st->old_is_transient[2] = st->old_is_transient[1];
    st->old_is_transient[1] = st->old_is_transient[0];
    st->old_is_transient[0] = is_transient;

    return;
}

/*--------------------------------------------------------------------------*
 * find_best_delay()
 *
 *
 *--------------------------------------------------------------------------*/

static short find_best_delay(
    float *mu_o,
    float *in,
    short mind1,
    short maxd1,
    short lin,
    short delta,
    short *false_flag

)
{
    short i, d1, k;
    short d1m = 0;
    float min_sq_cross, min_corr;
    float accA, accB;
    float Rxy[MAXDELAY_FEC], Ryy[MAXDELAY_FEC];

    for( k = 0, d1 = mind1; k < (maxd1-mind1)/delta; d1 += delta, k++ )
    {
        accA = accB = 0;
        for( i = 0; i < lin; i += delta )
        {
            accA += mu_o[d1+i] * mu_o[d1+i];
            accB += mu_o[d1+i] * in[i];
        }

        Rxy[k] = accB;
        Ryy[k] = accA;
    }

    /*  Obtain the best delay values */
    min_sq_cross = -FLT_MAX;
    min_corr = 0;

    for( d1 = 0; d1 < (maxd1-mind1)/delta; d1++ )
    {
        if( Rxy[d1] * min_corr >= min_sq_cross * Ryy[d1] )
        {
            d1m = d1;
            min_corr = Ryy[d1];
            min_sq_cross = Rxy[d1];
        }
    }
    d1m *= delta;


    accA =  min_sq_cross/min_corr;
    if(accA < 0.5 || accA > 1.5)
    {
        *false_flag = 1;
    }
    else
    {
        *false_flag = 0;
    }

    return d1m;
}

/*--------------------------------------------------------------------------*
 * Search_Max_Corr()
 *
 *
 *--------------------------------------------------------------------------*/

static short Search_Max_Corr(
    float *mu_o,                    /* i : *old_auOut_2fr       */
    short old_Min_ind,              /* i : * Old Minimum index  */
    const short L                   /* i : L/2                  */
)
{
    short pos;
    short pos2,delta2;
    short lin, delta;
    short mind1,maxd1;
    float *in;
    short false_flag;

    short min_d1, max_d1;

    if(old_Min_ind == 0)
    {
        lin = 8*L/20;               /* Basic size of the block for phase matching */

        mind1 = 0;                  /* min value of delay d1 to search for */
        maxd1 = 12*L/20;            /* max value of delay d1 to search for */

        in = mu_o + 2*L - lin;

        /* generate correlation */
        delta  = 2;
        delta2 = 1;

        pos = find_best_delay(mu_o, in, mind1, maxd1, lin, delta, &false_flag);

        if(false_flag)
        {
            return 0;
        }
        min_d1 = max(mind1,mind1+pos-delta+1);
        max_d1 = min(maxd1,mind1+pos+delta);
        pos2 = find_best_delay(mu_o, in, min_d1, max_d1, lin, delta2, &false_flag);

        if(mind1 > (mind1+pos-delta+1))
        {
            pos = pos2;
        }
        else
        {
            pos = pos + pos2 - delta + 1 ;
        }
        pos = pos+lin+mind1;
    }
    else
    {
        lin = 6*L/20;

        mind1 =  9*L/20;                /* min value of delay d1 to search for */
        maxd1 = 14*L/20;                /* max value of delay d1 to search for */

        in = mu_o + 2*L - lin;

        /* generate correlation */
        delta  = 2;
        delta2 = 1;

        pos = find_best_delay(mu_o, in, mind1, maxd1, lin, delta, &false_flag);

        if(false_flag)
        {
            return 0;
        }

        min_d1 = max(mind1, mind1+pos-delta+1);
        max_d1 = min(maxd1, mind1+pos+delta);
        pos2 = find_best_delay(mu_o, in, min_d1, max_d1, lin, delta2, &false_flag);

        if(mind1 > (mind1+pos-delta+1) )
        {
            pos = pos2;
        }
        else
        {
            pos = pos + pos2 - delta + 1 ;
        }

        pos = pos+lin+mind1;
    }

    return pos;
}

/*--------------------------------------------------------------------------*
 * FEC_phase_matching()
 *
 *
 *--------------------------------------------------------------------------*/

static short FEC_phase_matching(
    Decoder_State *st,                      /* i  : Decoder State                           */
    float *ImdctOut,                        /* i  : input                                   */
    float *auOut,                           /* o  : output audio                            */
    float *OldauOut,
    float OldauOut_pha[2][N_LEAD_NB]
)
{
    short i;
    float ImdctOutWin[2*L_FRAME8k];
    short pos, remain;
    short ol_size;
    float OldauOutnoWin[L_FRAME8k];
    short L_overlap, L;
    float OldauOut2[L_FRAME8k];
    float pow1=0, pow2=0;
    float win_NB[L_FRAME8k + 25];
    float SmoothingWin_NB3[24];

    L = L_FRAME8k;

    for(i=0; i<3*L/20; i++)
    {
        SmoothingWin_NB3[i] = SmoothingWin_NB875[i*3];
    }

    for(i=0; i<L + 25; i++)
    {
        win_NB[i] = window_48kHz[i*6+3];
    }

    set_f(ImdctOutWin, 0.0, 2*L);

    /* OLA */
    ol_size = 2*L/20;
    pos = Search_Max_Corr(st->old_auOut_2fr, st->old_Min_ind, L);

    if(pos == 0)
    {
        return 1;
    }

    /* Repetition */
    remain = L+N_Z_L_NB - ((2*L)-pos);
    mvr2r(&st->old_auOut_2fr[pos], &ImdctOutWin[N_ZERO_NB], (2*L)-pos);

    /* OldauOut without windowing */
    for(i = N_ZERO_NB; i < L/2; i++)
    {
        OldauOutnoWin[i-N_ZERO_NB] = -st->oldIMDCTout[L/2 - 1 - i];
    }
    for(i = 0; i < L/2; i++)
    {
        OldauOutnoWin[i+N_ZERO_O_NB] = -st->oldIMDCTout[i];
    }

    mvr2r(OldauOutnoWin, &ImdctOutWin[N_ZERO_NB+(2*L)-pos], remain);

    for(i = 0; i < L; i++)
    {
        pow1 += (float)fabs(st->old_auOut_2fr[L+i]);
        pow2 += (float)fabs(ImdctOutWin[N_ZERO_NB+i]);
    }
    if(pow2!=0)
    {
        pow1 /= pow2;
        for(i = N_ZERO_NB; i < 2*L; i++)
        {
            ImdctOutWin[i] *= pow1;
        }
    }
    Smoothing_vector_NB(OldauOutnoWin, &ImdctOutWin[N_ZERO_NB], SmoothingWin_NB2, auOut, ol_size);

    for (i = 0; i < L/2; i++)
    {
        ImdctOutWin[3*L/2 + i] *= win_NB[L/2-i-1];
    }

    for (i = N_ZERO_NB; i < L/2; i++)
    {
        ImdctOutWin[L + i] *= win_NB[(L-1-i)];
    }
    mvr2r(&ImdctOutWin[N_Z_L_O_NB], &OldauOut_pha[0][0], N_LEAD_NB);
    mvr2r(&ImdctOutWin[ol_size+N_ZERO_NB], &auOut[ol_size], N_Z_L_NB-ol_size);
    mvr2r(&ImdctOutWin[L], &auOut[N_Z_L_NB], N_ZERO_NB);
    mvr2r(&ImdctOutWin[L], OldauOut, L);

    for(i = 0; i < L/2; i++)
    {
        OldauOut2[i] = -ImdctOut[L/2 - 1 - i];
        OldauOut2[L/2+i] = -ImdctOut[i];
    }

    L_overlap = 3*L/20;
    Smoothing_vector_NB(&ImdctOutWin[N_Z_L_O_NB], &OldauOut2[N_ZERO_NB], SmoothingWin_NB3, &OldauOut_pha[1][0], L_overlap);

    for(i=L_overlap; i<N_LEAD_NB; i++)
    {
        OldauOut_pha[1][i] = OldauOut2[i+N_ZERO_NB];
    }

    return 0;
}

/*--------------------------------------------------------------------------*
 * FEC_phase_matching_nextgood()
 *
 *
 *--------------------------------------------------------------------------*/

void FEC_phase_matching_nextgood(
    const float *ImdctOut,                /* i  : input                           */
    float *auOut,                         /* o  : output audio                    */
    float *OldauOut,                      /* i/o: audio from previous frame       */
    float OldauOut_pha[2][N_LEAD_NB],
    float mean_en_high
)
{
    short i;
    float ImdctOutWin[2*L_FRAME8k];
    short L_overlap, L;
    short oldout_pha_idx;
    float *OldOut_pha;
    float win_NB[L_FRAME8k + 25];

    L = L_FRAME8k;
    for(i=0; i<L + 25; i++)
    {
        win_NB[i] = window_48kHz[i*6+3];
    }

    if((mean_en_high>2.f)||(mean_en_high<0.5f))
    {
        oldout_pha_idx = 1;
    }
    else
    {
        oldout_pha_idx = 0;
    }

    /* Overlapping with next good frame : Overlapping to remove the discontinuity */
    L_overlap = N_LEAD_NB;
    OldOut_pha = OldauOut_pha[oldout_pha_idx];
    for (i = 0; i < N_LEAD_NB; i++)
    {
        OldOut_pha[i] *= SmoothingWin_NB875[L_overlap-i-1];
    }

    if(oldout_pha_idx == 1)
    {
        /* Use phase matching and overlapping with the Oldauout*/
        /* Windowing */
        Windowing_1st_NB(ImdctOutWin, ImdctOut, win_NB, NULL, 0);
        Windowing_2nd_NB(ImdctOutWin, ImdctOut, win_NB);
    }
    else
    {
        /* Only use phase matching */
        /* Windowing */
        Windowing_1st_NB(ImdctOutWin, ImdctOut, win_NB, SmoothingWin_NB875, 1);
        Windowing_2nd_NB(ImdctOutWin, ImdctOut, win_NB);
    }

    common_overlapping(auOut, ImdctOutWin, OldOut_pha, N_LEAD_NB, 0, N_LEAD_NB, L, N_ZERO_NB, 0);
    mvr2r(&ImdctOutWin[L], OldauOut, L);

    return;
}

/*--------------------------------------------------------------------------*
 * FEC_phase_matching_burst()
 *
 *
 *--------------------------------------------------------------------------*/

static void FEC_phase_matching_burst(
    const float *ImdctOut,              /* i  : input                           */
    float *auOut,                       /* o  : output audio                    */
    float *OldauOut,                    /* i/o: audio from previous frame       */
    float OldauOut_pha[2][N_LEAD_NB],
    float *prev_oldauOut                /* i : OldauOut from previous frame     */
)
{
    short i;
    short L_overlap;
    float OldauOut2[L_FRAME8k];
    float ImdctOutWin[2*L_FRAME8k];
    short L;
    float win_NB[L_FRAME8k + 25];
    float SmoothingWin_NB3[24];

    L = L_FRAME8k;

    for(i=0; i<3*L/20; i++)
    {
        SmoothingWin_NB3[i] = SmoothingWin_NB875[i*3];
    }

    for(i=0; i<L + 25; i++)
    {
        win_NB[i] = window_48kHz[i*6+3];
    }

    /* Windowing */
    Windowing_1st_NB(ImdctOutWin, ImdctOut, win_NB, NULL, 0);
    Windowing_2nd_NB(ImdctOutWin, ImdctOut, win_NB);

    /* Repetition with old frame to reserve energy */
    common_overlapping(auOut, ImdctOutWin, prev_oldauOut, N_Z_L_NB, 0, N_Z_L_NB, L, N_ZERO_NB, 0);

    /* data transition from OldauOut to auOut using smoothing win*/
    Smoothing_vector_NB(OldauOut_pha[0], auOut, SmoothingWin_NB875, auOut, N_LEAD_NB);

    /* Update the OldauOut array for next overlapping */
    mvr2r(&ImdctOutWin[N_Z_L_O_NB], &OldauOut_pha[0][0], N_LEAD_NB);
    mvr2r(&ImdctOutWin[L], OldauOut, L);
    Scaledown(prev_oldauOut, prev_oldauOut, SCALE_DOWN_3dB, L);

    for(i = 0; i < L/2; i++)
    {
        OldauOut2[i] = -ImdctOut[L/2 - 1 - i];
        OldauOut2[L/2+i] = -ImdctOut[i];
    }

    L_overlap = 3*L/20;
    Smoothing_vector_NB(&ImdctOutWin[N_Z_L_O_NB], &OldauOut2[N_ZERO_NB], SmoothingWin_NB3, &OldauOut_pha[1][0], L_overlap);

    for(i=L_overlap; i<N_LEAD_NB; i++)
    {
        OldauOut_pha[1][i] = OldauOut2[i+N_ZERO_NB];
    }

    return;
}


/*--------------------------------------------------------------------------*
 * Repetition_smoothing_nextgood()
 *
 *
 *--------------------------------------------------------------------------*/

static void Repetition_smoothing_nextgood(
    const float *ImdctOut,          /* i  : input                           */
    float *auOut,                   /* o  : output audio                    */
    float *OldImdctOut,             /* i  : input                           */
    float *OldauOut,                /* i/o: audio from previous frame       */
    short cur_data_use_flag,        /* i  : current imdct data use flag     */
    short overlap_time
)
{
    short i;
    float ImdctOutWin[2*L_FRAME8k];
    float win_NB[L_FRAME8k + 25];
    short L_overlap;
    short ol_size;
    short L;

    L = L_FRAME8k;

    for(i=0; i<L_FRAME8k + 25; i++)
    {
        win_NB[i] = window_48kHz[i*6+3];
    }

    for (i = N_ZERO_NB; i < L/2; i++)
    {
        OldauOut[i-N_ZERO_NB] = -OldImdctOut[L/2 - 1 - i];
    }
    for(i = 0; i < L/2; i++)
    {
        OldauOut[i+N_ZERO_O_NB] = -OldImdctOut[i];
    }

    /* Overlapping with next good frame : Overlapping to remove the discontinuity */
    if(cur_data_use_flag)
    {
        ol_size = N_LEAD_NB;

        for (i = N_ZERO_NB; i < L/2; i++)
        {
            ImdctOutWin[i+L] = -ImdctOut[L/2 - 1 - i];
        }
        for(i = 0; i < L/2; i++)
        {
            ImdctOutWin[i+3*L/2] = -ImdctOut[i];
        }

        /*a = (float)(1./(float)(ol_size)); y = ax */
        Smoothing_vector_scaledown_NB(OldauOut, &ImdctOutWin[N_Z_L_O_NB], SmoothingWin_NB875, OldauOut, ol_size);

        /* Scale down the overlapped signal */
        Scaledown(&ImdctOutWin[ol_size+N_Z_L_O_NB], &OldauOut[ol_size], SCALE_DOWN_3dB, N_Z_L_NB-ol_size);
    }

    L_overlap = overlap_time;
    for (i = 0; i < L_overlap; i++)
    {
        OldauOut[i] *= SmoothingWin_NB875[L_overlap-i-1];
    }
    for(i=L_overlap; i < L; i++)
    {
        OldauOut[i] = 0.f;
    }

    /* Windowing */
    Windowing_1st_NB(ImdctOutWin, ImdctOut, win_NB, SmoothingWin_NB875, 1);
    Windowing_2nd_NB(ImdctOutWin, ImdctOut, win_NB);

    v_add(&ImdctOutWin[N_ZERO_NB], OldauOut, auOut, L);
    mvr2r(&ImdctOutWin[L], OldauOut, L);

    return;
}

/*--------------------------------------------------------------------------*
 * Repetition_smoothing()
 *
 *
 *--------------------------------------------------------------------------*/

static int Repetition_smoothing(
    const float *ImdctOut,        /* i  : input                           */
    float *auOut,                 /* o  : output audio                    */
    float *OldImdctOut,           /* i  : input                           */
    float *OldauOut,              /* i/o: audio from previous frame       */
    const short L,                /* i  : length                          */
    float *prev_oldauOut,         /* i : OldauOut from previous frame     */
    short overlap_time            /* i : overlap time                     */
)
{
    short i;
    float ImdctOutWin[2*L_FRAME8k];
    float pow1 = 0.f;
    float pow2 = 0.f;
    float OldauOutnoWin[L_FRAME8k];
    float win_NB[L_FRAME8k + 25];

    for(i=0; i<L_FRAME8k + 25; i++)
    {
        win_NB[i] = window_48kHz[i*6+3];
    }

    /* Windowing */
    Windowing_1st_NB(ImdctOutWin, ImdctOut, win_NB, NULL, 0);
    Windowing_2nd_NB(ImdctOutWin, ImdctOut, win_NB);

    /* Repetition with old frame to reserve energy */
    common_overlapping( auOut, ImdctOutWin, prev_oldauOut, N_Z_L_NB, 0, N_Z_L_NB, L, N_ZERO_NB, 0);

    /* OldauOut without windowing */
    for (i = N_ZERO_NB; i < L/2; i++)
    {
        OldauOutnoWin[i-N_ZERO_NB] = -OldImdctOut[L/2 - 1 - i];
    }
    for(i = 0; i < L/2; i++)
    {
        OldauOutnoWin[i+N_ZERO_O_NB] = -OldImdctOut[i];
    }

    /* data transition from OldauOut to auOut using smoothing win*/
    Smoothing_vector_NB(OldauOutnoWin, auOut, SmoothingWin_NB875, auOut, overlap_time);

    pow1 = sum2_f(&auOut[1*L/20], 4*L/20);
    pow2 = sum2_f(&auOut[N_LEAD_NB], 4*L/20);


    if(pow2 > pow1*3)
    {
        return 1;
    }

    /* Update the OldauOut array for next overlapping */
    mvr2r(&ImdctOutWin[L], OldauOut, L);
    Scaledown(prev_oldauOut, prev_oldauOut, SCALE_DOWN_3dB, L);

    return 0;
}


/*--------------------------------------------------------------------------*
 * Windowing_1st()
 *
 *
 *--------------------------------------------------------------------------*/

static void Windowing_1st_NB(
    float *ImdctOutWin,                      /* o : Output                                  */
    const float *ImdctOut,                   /* i : Input                                   */
    const float *win,                        /* i : Window                                  */
    const float *smoothingWin,               /* i : Smoothing Window                        */
    short smoothing_flag                     /* i : 1=Smoothing window, 0=Original window   */
)
{
    short i;
    short L;

    L = L_FRAME8k;
    if( smoothing_flag == 0 )
    {
        for (i = N_ZERO_NB; i < L/2; i++)
        {
            ImdctOutWin[i] = ImdctOut[L/2 + i] * win[(2*L-1-i)-N_LEAD_O_NB];
        }

        for (i = 0; i < N_ZERO_O_NB; i++)
        {
            ImdctOutWin[L/2 + i] = -ImdctOut[L - 1 - i] * win[(3*L/2-1-i)-N_LEAD_O_NB];
            ImdctOutWin[3*L/2 + i] = -ImdctOut[i] * win[(L/2-i-1)];
        }
    }
    else
    {
        for (i = N_ZERO_NB; i < L/2; i++)
        {
            ImdctOutWin[i] = ImdctOut[L/2 + i] * smoothingWin[(i-N_ZERO_NB)]; /*win[(2*L-i)*decimate-1-decay-14*L_FRAME48k/20];*/
        }

        for (i = 0; i < N_ZERO_O_NB; i++)
        {
            ImdctOutWin[L/2 + i] = -ImdctOut[L - 1 - i] * smoothingWin[(i+N_ZERO_O_NB)];   /*win[(3*L/2-1-i)*decimate+decay-L_FRAME48k*14/20];*/
            ImdctOutWin[3*L/2 + i] = -ImdctOut[i] * win[(L/2-i-1)];
        }
    }

    return;
}

/*--------------------------------------------------------------------------*
 * Windowing_2nd()
 *
 *
 *--------------------------------------------------------------------------*/

static void Windowing_2nd_NB(
    float *ImdctOutWin,         /* o : Output                   */
    const float *ImdctOut,      /* i : Input                    */
    const float *win           /* i : Window                   */
)
{
    short i;
    short L;

    L = L_FRAME8k;
    for (i = N_ZERO_O_NB; i < L/2; i++)
    {
        ImdctOutWin[L/2 + i] = -ImdctOut[L - 1 - i];
        ImdctOutWin[3*L/2 + i] = -ImdctOut[i] * win[L/2-i-1];
    }

    for (i = 0; i < N_ZERO_NB; i++)
    {
        ImdctOutWin[L + i] = -ImdctOut[L/2 - 1 - i];
    }

    for (i = N_ZERO_NB; i < L/2; i++)
    {
        ImdctOutWin[L + i] = -ImdctOut[L/2 - 1 - i] * win[L - 1 - i];
    }

    return;
}

/*--------------------------------------------------------------------------*
 * common_overlapping()
 *
 *
 *--------------------------------------------------------------------------*/

static void common_overlapping(
    float *auOut,           /* i : Input                                   */
    float *ImdctOutWin,     /* o : Output                                  */
    float *OldauOut,        /* i : Window                                  */
    short end1,             /* i : Decay                                   */
    short offset1,
    short start2,
    short end2,
    short offset_i2,
    short offset2
)
{
    short i;

    /* Common Overlapping */
    for (i=0 ; i < end1; i++)
    {
        auOut[i] = ImdctOutWin[i+N_ZERO_NB] + OldauOut[i+offset1];
    }
    for (i=start2 ; i < end2; i++)
    {
        auOut[i+offset2]  =  ImdctOutWin[i+offset_i2];
    }

    return;
}


/*--------------------------------------------------------------------------*
 * Smoothing_vector()
 *
 *
 *--------------------------------------------------------------------------*/

static void Smoothing_vector_NB(
    const float OldauOutnoWin[],  /* i  : Input vector 1                                   */
    const float ImdctOutWin[],    /* i  : Input vector 2                                   */
    const float SmoothingWin[],   /* i  : Smoothing window                                 */
    float auOut[],          /* o  : Output vector that contains vector 1 .* vector 2 */
    const short ol_size           /* i  : Overlap size                                     */
)
{
    short i;
    float weight;

    for (i=0 ; i < ol_size; i++)
    {
        weight = SmoothingWin[i];
        auOut[i] = (OldauOutnoWin[i]*(1.f-weight)) + (ImdctOutWin[i]*weight);
    }

    return;
}

/*--------------------------------------------------------------------------*
 * Smoothing_vector_scaledown()
 *
 *
 *--------------------------------------------------------------------------*/

static void Smoothing_vector_scaledown_NB(
    const float OldauOutnoWin[],  /* i  : Input vector 1                                   */
    const float ImdctOutWin[],    /* i  : Input vector 2                                   */
    const float SmoothingWin[],   /* i  : Smoothing window                                 */
    float auOut[],          /* o  : Output vector that contains vector 1 .* vector 2 */
    const short ol_size           /* i  : Overlap size                                     */
)
{
    short i;
    float weight;

    for (i=0 ; i < ol_size; i++)
    {
        weight = SmoothingWin[i];
        auOut[i] = (OldauOutnoWin[i]*(1.f-weight)) + (ImdctOutWin[i]*SCALE_DOWN_3dB*weight);
    }

    return;
}

/*--------------------------------------------------------------------------*
 * Scaledown()
 *
 *
 *--------------------------------------------------------------------------*/

static void Scaledown(
    float x[],          /* i  : Input vector                                  */
    float y[],          /* o  : Output vector that contains vector 1 .* vector 2 */
    float scale_v,
    const short N             /* i  : Overlap size                                    */
)
{
    short i;

    for(i=0; i<N; i++)
    {
        y[i] = x[i]*scale_v;
    }

    return;
}

/*--------------------------------------------------------------------------*
 * time_domain_FEC_HQ()
 *
 *
 *--------------------------------------------------------------------------*/

void time_domain_FEC_HQ(
    Decoder_State *st,          /* i : Decoder State                        */
    float *wtda_audio,          /* i : input                                */
    float *out,                 /* o : output audio                         */
    float mean_en_high,         /* i : transient flag                       */
    const short output_frame
)
{
    if( (st->nbLostCmpt==1)&&(st->phase_mat_flag==1)&&(st->phase_mat_next==0) )
    {
        if (FEC_phase_matching(st, wtda_audio, out, st->old_out, st->old_out_pha) )
        {
            window_ola( wtda_audio, out, st->old_out, output_frame, st->tcx_cfg.tcx_last_overlap_mode, st->tcx_cfg.tcx_curr_overlap_mode, st->prev_bfi, st->oldHqVoicing, st->oldgapsynth );
            st->phase_mat_next = 0;
        }
        else
        {
            st->phase_mat_next = 1;
        }
    }
    else if((st->prev_bfi == 1)&&(st->bfi == 0) &&(st->phase_mat_next == 1))
    {
        FEC_phase_matching_nextgood( wtda_audio, out, st->old_out, st->old_out_pha, mean_en_high);
        st->phase_mat_next = 0;
    }
    else if((st->prev_bfi == 1)&&(st->bfi == 1) &&(st->phase_mat_next == 1))
    {
        FEC_phase_matching_burst( wtda_audio, out, st->old_out, st->old_out_pha, st->prev_oldauOut);
        st->phase_mat_next = 1;
    }
    else
    {
        if(st->bfi == 0 && st->prev_bfi == 1)
        {
            if((st->stat_mode_out==1) || (st->diff_energy<ED_THRES_50P) )/* st->diff_energy<ED_THRES_L1)*/
            {
                Repetition_smoothing_nextgood( wtda_audio, out, st->oldIMDCTout, st->old_out, st->old_bfi_cnt> 1 ? 1:0, N_LEAD_NB);
            }
            else if(st->old_bfi_cnt > 1)
            {
                Next_good_after_burst_erasures( wtda_audio, out, st->old_out, N_LEAD_NB );
            }
            else
            {
                window_ola( wtda_audio, out, st->old_out, output_frame, st->tcx_cfg.tcx_last_overlap_mode, st->tcx_cfg.tcx_curr_overlap_mode, st->prev_bfi, st->oldHqVoicing, st->oldgapsynth );
            }
        }
        else /* if(st->bfi == 1) */
        {
            if( (st->stat_mode_out==1) || (st->diff_energy<ED_THRES_50P) )
            {
                if( Repetition_smoothing( wtda_audio, out, st->oldIMDCTout, st->old_out, output_frame, st->prev_oldauOut, N_LEAD_NB) )
                {
                    window_ola( wtda_audio, out, st->old_out, output_frame, st->tcx_cfg.tcx_last_overlap_mode, st->tcx_cfg.tcx_curr_overlap_mode, st->prev_bfi, st->oldHqVoicing, st->oldgapsynth );
                }
            }
            else
            {
                window_ola( wtda_audio, out, st->old_out, output_frame, st->tcx_cfg.tcx_last_overlap_mode, st->tcx_cfg.tcx_curr_overlap_mode, st->prev_bfi, st->oldHqVoicing, st->oldgapsynth );
            }
        }

        st->phase_mat_next = 0;
    }

    return;
}


/*--------------------------------------------------------------------------*/
/*  Function  Next_good_after_burst_erasures                                */
/*  ~~~~~~~~~~~~~~~~~~~~                                                    */
/*                                                                          */
/*  Windowing, Overlap and Add                                              */
/*--------------------------------------------------------------------------*/
/*  float     ImdctOut[]  (i)   input                                       */
/*  short     auOut[]     (o)   output audio                                */
/*  float     OldauOut[]  (i/o) audio from previous frame                   */
/*  short     ol_size     (i)   overlap size                                */
/*  short     L           (i)   length                                      */
/*--------------------------------------------------------------------------*/

void Next_good_after_burst_erasures(
    const float *ImdctOut,
    float *auOut,
    float *OldauOut,
    const short ol_size
)
{
    float ImdctOutWin[2*L_FRAME8k];
    short i, L;
    float win_NB[L_FRAME8k + 25];

    L = L_FRAME8k;
    for(i=0; i<L + 25; i++)
    {
        win_NB[i] = window_48kHz[i*6+3];
    }

    /* Windowing */
    Windowing_1st_NB(ImdctOutWin, ImdctOut, win_NB, NULL, 0);
    Windowing_2nd_NB(ImdctOutWin, ImdctOut, win_NB);

    /* Overlapping with next good frame : Overlapping to remove the discontinuity */
    Smoothing_vector_scaledown_NB(&OldauOut[N_ZERO_NB], &ImdctOutWin[N_Z_L_O_NB], SmoothingWin_NB875, &OldauOut[N_ZERO_NB], ol_size);

    /* Scale down the overlapped signal */
    Scaledown(&ImdctOutWin[ol_size+N_Z_L_O_NB], &OldauOut[ol_size+N_ZERO_NB], SCALE_DOWN_3dB, N_Z_L_NB-ol_size);

    /* Common Overlapping */
    common_overlapping(auOut, ImdctOutWin, OldauOut, N_Z_L_NB, N_ZERO_NB, 0, N_ZERO_NB, L, N_Z_L_NB);
    mvr2r( &ImdctOutWin[L], OldauOut, L );

    return;
}
