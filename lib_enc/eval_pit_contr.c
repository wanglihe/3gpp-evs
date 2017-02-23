/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"


/*-------------------------------------------------------------------*
 * Local constantes
 *-------------------------------------------------------------------*/

#define NB_VOIC         13
#define DIV_NB_VOIC     (1.0f/NB_VOIC)
#define ALPA            0.95f
#define ALPAM1          (1.0f-ALPA)
#define BETA            (ALPAM1/2.0f)
#define AFREQ_THR       2
#define GPIT_THR        0.4f
#define HANGOVER_DELAY  2

/*-------------------------------------------------------------------*
 * Pit_exc_contribution_len()
 *
 * Determine up to which band the pit contribution is significant
 *-------------------------------------------------------------------*/

short Pit_exc_contribution_len(   /* o  : bin where pitch contribution is significant */
    Encoder_State *st,              /* i/o: state structure                                  */
    const float *dct_res,         /* i  : DCT of residual                 */
    float *dct_pitex,       /* i/o: DCT of pitch contribution       */
    float *pitch_buf,       /* i/o: Pitch per subframe              */
    short *hangover,        /* i  : hangover for the time contribution switching */
    const short coder_type        /* i  : coding type */
)
{
    float corr_dct_pit[MBANDS_LOC], corr_tmp;
    float av_corr, min_corr, ftmp;
    short freq, i, j;
    short last_pit_band, pit_contr_idx, last_pit_bin;
    float ener_res;
    float ener_pit;
    float low_pit, F1st_harm, F8th_harm;
    float corr_dct_pit_tmp[MBANDS_LOC];
    short time_flg = 0;
    short Len, max_len;
    short tmp_dec;
    short Mbands_loc = MBANDS_LOC-2;

    if(st->L_frame == L_FRAME16k )
    {
        Mbands_loc = MBANDS_LOC;
    }

    minimum( pitch_buf, st->L_frame >> 6, &low_pit );

    F1st_harm = 12800.0f/low_pit;
    F8th_harm = 8.0f*F1st_harm;

    freq = 0;
    for (i = 0; i <Mbands_loc; i++)
    {
        corr_tmp = 0.0f;
        ener_res = 0.1f;
        ener_pit = 0.1f;

        for (j = 0; j < mfreq_bindiv_loc[i]; j++)
        {
            corr_tmp+=dct_res[j+freq]*dct_pitex[j+freq];
            ener_res+=dct_res[j+freq]*dct_res[j+freq];
            ener_pit+=dct_pitex[j+freq]*dct_pitex[j+freq];
        }

        corr_dct_pit[i] = (float)(corr_tmp / sqrt(ener_res*ener_pit));
        freq += mfreq_bindiv_loc[i];
    }

    /* Smooth the inter-correlation value and skip the last band for the average (since last band is almost always 0)*/
    corr_dct_pit_tmp[0] = ALPA*corr_dct_pit[0]+ALPAM1*corr_dct_pit[1];
    if( corr_dct_pit_tmp[0] < 0.5f )
    {
        corr_dct_pit_tmp[0] = 0.5f;
    }
    corr_dct_pit_tmp[0] = (corr_dct_pit_tmp[0]-0.5f) *2.0f ;

    for (i = 1; i <Mbands_loc-1; i++)
    {
        corr_dct_pit_tmp[i] = ALPA*corr_dct_pit[i]+BETA*corr_dct_pit[i+1]+BETA*corr_dct_pit[i-1];
        if(corr_dct_pit_tmp[i] < 0.5f)
        {
            corr_dct_pit_tmp[i] = 0.5f;
        }
        corr_dct_pit_tmp[i] = (corr_dct_pit_tmp[i]-0.5f) *2.0f ;
    }

    corr_dct_pit_tmp[i] = ALPA*corr_dct_pit[i]+ALPAM1*corr_dct_pit[i-1];

    if( corr_dct_pit_tmp[i] < 0.5f )
    {
        corr_dct_pit_tmp[i] = 0.5f;
    }
    corr_dct_pit_tmp[i] = (corr_dct_pit_tmp[i] - 0.5f) * 2.0f ;

    for (i = 0; i <Mbands_loc; i++)
    {
        corr_dct_pit[i] = corr_dct_pit_tmp[i];
    }

    av_corr = DIV_NB_VOIC * corr_dct_pit[0];
    for (i = 1; i <NB_VOIC; i++)
    {
        av_corr += DIV_NB_VOIC*corr_dct_pit[i];
    }

    /* Find the cut-off freq similarly to HSX */
    last_pit_band = 0;

    av_corr *= 6400;

    if( st->core_brate < ACELP_9k60 )
    {
        /* Correlation really poor at low rate, time domain still valide */
        av_corr *= 2.0;
    }

    min_corr =(float)fabs(mfreq_loc[0] - av_corr);

    for (i = 1; i <Mbands_loc; i++)
    {
        ftmp = (float)fabs(mfreq_loc[i] - av_corr);
        if( ftmp < min_corr )
        {
            last_pit_band = i;
            min_corr = ftmp;
        }
    }

    if( F8th_harm > mfreq_loc[last_pit_band] )
    {
        do
        {
            last_pit_band++;
        }
        while( F8th_harm >=  mfreq_loc[last_pit_band] );
    }

    if( last_pit_band > 7+BAND1k2 && (st->core_brate < CFREQ_BITRATE || st->bwidth == NB) )
    {
        last_pit_band = 7+BAND1k2;
    }
    else if ( last_pit_band > 10+BAND1k2 && st->core_brate >= CFREQ_BITRATE )
    {
        last_pit_band = 10+BAND1k2;
    }

    time_flg = 0;
    if( (st->mem_last_pit_band > 0 && st->old_corr > 0.5f && st->mold_corr > 0.5f && st->lt_gpitch >= 1.5f*GPIT_THR )
            || (last_pit_band > 6)
            || (last_pit_band >= 4 && st->lt_gpitch >= 1.5f*GPIT_THR && st->old_corr > 0.7f)
            || (last_pit_band > BAND1k2  &&  st->mold_corr > 0.80f && st->lt_gpitch >= GPIT_THR) )
    {
        tmp_dec = 1;
    }
    else
    {
        tmp_dec = 0;
    }

    /* Different past and current decision */
    if ( (st->mem_last_pit_band == 0 && tmp_dec == 1) || (st->mem_last_pit_band > 0 && tmp_dec == 0) )
    {
        if( *hangover == 0 )
        {
            time_flg = tmp_dec;
            *hangover = HANGOVER_DELAY;
        }
        else
        {
            time_flg = 0;
            if( st->mem_last_pit_band > 0 )
            {
                time_flg = 1;
            }

            (*hangover) -= 1;
            if( *hangover < 0 )
            {
                *hangover = 0;
            }
        }
    }
    else
    {
        time_flg = tmp_dec;
        *hangover = HANGOVER_DELAY;
    }

    /* Decicison on final lenght of time contribution */
    pit_contr_idx = 0;
    if( time_flg == 1 || coder_type != INACTIVE || st->GSC_noisy_speech )
    {
        if( st->core_brate  <ACELP_9k60 && low_pit < 64 )
        {
            last_pit_band = 9+BAND1k2;
            if(st->bwidth == NB)
            {
                last_pit_band = 7+BAND1k2;
            }
        }
        else if(st->core_brate  < ACELP_9k60 && low_pit < 128)
        {
            last_pit_band = 5+BAND1k2;
        }
        else if(st->core_brate  < ACELP_9k60 )
        {
            last_pit_band = 3+BAND1k2;
        }
        else if( last_pit_band < BAND1k2+1 )
        {
            last_pit_band = BAND1k2+1;
        }

        last_pit_bin = (short)(mfreq_loc[last_pit_band]/BIN_SIZE);

        st->bpf_off = 0;

        max_len = st->L_frame-last_pit_bin;
        if( st->bwidth == NB )
        {
            max_len = 160-last_pit_bin;
        }

        Len = 80;
        if( max_len < 80 )
        {
            Len = max_len;
        }

        if(st->core_brate == ACELP_8k00 && st->bwidth != NB )
        {
            for (i=0; i < max_len; i++)
            {
                dct_pitex[i+last_pit_bin] = 0.0f;
            }
        }
        else
        {
            for (i = 0; i < Len; i++)
            {
                dct_pitex[i+last_pit_bin] *= sm_table[i];
            }

            for (; i < max_len; i++)
            {
                dct_pitex[i+last_pit_bin] = 0.0f;
            }
        }

        st->mem_last_pit_band = last_pit_band;
        pit_contr_idx = last_pit_band-BAND1k2;
    }
    else
    {
        set_f(dct_pitex, 0.0f, st->L_frame);
        st->bpf_off = 1;
        last_pit_bin = 0;
        last_pit_band = 0;
        pit_contr_idx = 0;
        st->mem_last_pit_band = 0;
        set_f( pitch_buf, (float)L_SUBFR, NB_SUBFR );

        /* pitch contribution useless - delete all previously written indices belonging to pitch contribution */
        for( i = TAG_ACELP_SUBFR_LOOP_START; i < TAG_ACELP_SUBFR_LOOP_END; i++ )
        {
            if( st->ind_list[i].nb_bits != -1 )
            {
                st->nb_bits_tot -= st->ind_list[i].nb_bits;
                st->ind_list[i].nb_bits = -1;
            }
        }

        if( st->ind_list[IND_ES_PRED].nb_bits != -1 )
        {
            st->nb_bits_tot -= st->ind_list[IND_ES_PRED].nb_bits;
            st->ind_list[IND_ES_PRED].nb_bits = -1;
        }
    }

    if( st->core_brate < CFREQ_BITRATE )
    {
        if( st->core_brate < ACELP_9k60 )
        {
            if( pit_contr_idx > 0 )
            {
                pit_contr_idx = 1;
            }

            if( coder_type == INACTIVE )
            {
                push_indice( st, IND_PIT_CONTR_IDX, pit_contr_idx, 1 );
            }
        }
        else
        {
            push_indice( st, IND_PIT_CONTR_IDX, pit_contr_idx, 3 );
        }
    }
    else
    {
        push_indice( st, IND_PIT_CONTR_IDX, pit_contr_idx, 4 );
    }

    return last_pit_bin;
}
