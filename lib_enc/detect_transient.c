/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "prot.h"
#include "cnst.h"
#include "rom_com.h"


/*--------------------------------------------------------------------------
 * hp_filter()
 *
 * High pass filter
 *--------------------------------------------------------------------------*/

static void hp_filter(
    const float *x,         /* i  : input signal            */
    float *y,         /* o  : output signal           */
    float *oldx,      /* i/o: previous filter input   */
    float *oldy,      /* i/o: previous filter output  */
    const short L           /* i  : length (32 or 48Hz)     */
)
{
    short i;

    y[0] = 0.4931f * *oldy + 0.7466f*(x[0] - *oldx);

    for(i = 1; i < L; i++)
    {
        y[i] = 0.4931f*y[i-1] + 0.7466f*(x[i] - x[i-1]);
    }

    *oldx = x[L - 1];
    *oldy = y[L - 1];

    return;
}

/*--------------------------------------------------------------------------
 * detect_transient()
 *
 * Detect if the signal is a transient
 *--------------------------------------------------------------------------*/

short detect_transient(                /* o  : transient flag           */
    const float *in,                 /* i  : input signal             */
    Encoder_State *st,                 /* i/o: Encoder state structure  */
    const short L,                   /* i  : length (32 or 48kHz)     */
    const short coder_type           /* i  : coder type               */
)
{
    float Energy;
    float EnergyLT;
    short i, blk;
    short IsTransient;
    float out_filt[L_FRAME48k];
    short position = 0;
    float thr;
    float Thres = 0.f;
    float Energy_in[5];
    float E_low, E_high;
    float E_in = 0.0f, E_out = 0.0f;

    IsTransient = 0;

    if( st->last_extl != st->extl )
    {
        st->TransientHangOver = 0;
        st->old_hpfilt_in = 0;
        st->old_hpfilt_out = 0;
        st->Energy_Old = 0;
    }

    /* High-pass filter */
    hp_filter( in, out_filt, &(st->old_hpfilt_in), &(st->old_hpfilt_out), L );

    /* Long-term energy */
    if( st->last_extl != st->extl || (st->last_extl == st->extl && st->last_core != st->core) || st->last_codec_mode == MODE2 )
    {
        EnergyLT = EPSILON;
        for ( i = 0; i < L/4; i++ )
        {
            EnergyLT += out_filt[i] * out_filt[i];
        }
    }
    else
    {
        EnergyLT = st->EnergyLT;
    }

    if(L == L_FRAME8k)
    {
        Energy_in[0] = st->Energy_Old;
        E_in = 0;
        E_out = 0;

        /* Compute block energy */
        for ( blk = 0; blk < 4; blk++ )
        {
            Energy = EPSILON;
            Energy_in[blk+1] = EPSILON;

            for ( i = 0; i < L/4; i++ )
            {
                Energy += out_filt[i + blk*(L/4)] * out_filt[i + blk*(L/4)];
                Energy_in[blk+1] += in[i + blk*(L/4)] * in[i + blk*(L/4)];
            }

            E_in += Energy_in[blk+1];
            E_out += Energy;

            Thres = 15.f;

            if(Energy > 6.0f * EnergyLT)
            {
                IsTransient = 1;
                position = blk;
            }

            EnergyLT = 0.75f*EnergyLT + 0.25f*Energy;
        }
    }
    else
    {
        /* Compute block energy */
        for ( blk = 0; blk < 4; blk++ )
        {
            Energy = EPSILON;
            for ( i = 0; i < L/4; i++ )
            {
                Energy += out_filt[i + blk*(L/4)] * out_filt[i + blk*(L/4)];
            }

            if( st->extl == SWB_BWE || st->extl == FB_BWE )
            {
                if( ( Energy > 13.5f * EnergyLT ) || ( Energy > 10.0f * EnergyLT && coder_type == INACTIVE ) )
                {
                    IsTransient = 1;
                    position = blk;
                }
            }
            else
            {
                if( st->total_brate <= HQ_16k40 && st->bwidth == SWB )
                {
                    thr = 13.5f;
                }
                else
                {
                    thr = 6.0f;
                }

                if( Energy > thr * EnergyLT )
                {
                    IsTransient = 1;
                    position = blk;
                }
            }

            EnergyLT = 0.75f*EnergyLT + 0.25f*Energy;
        }
    }

    st->EnergyLT = EnergyLT;

    if( ( st->last_extl != SWB_BWE && st->last_extl != SWB_TBE && st->extl == SWB_BWE ) ||
            ( st->last_extl != FB_BWE && st->last_extl != FB_TBE && st->extl == FB_BWE ) )
    {
        IsTransient = 0;
    }

    if( IsTransient && L == L_FRAME8k )
    {
        E_low = 0.f;
        blk = 0;
        for(i=0; i<position+1; i++)
        {
            E_low += Energy_in[i];
            blk++;
        }
        E_low /= (float)blk;

        E_high = 0.f;
        blk = 0;
        for(i=position+1; i<5; i++)
        {
            E_high += Energy_in[i];
            blk++;
        }
        E_high /= (float)blk;

        if( ((E_high/E_low<2.0f)&&(E_high/E_low>0.7f)) && ((E_in/E_out)>Thres) )
        {
            IsTransient = 0;
        }
    }

    if ( IsTransient )
    {
        if(L == L_FRAME8k)
        {
            if( position == 3 )
            {
                st->TransientHangOver = 1;
            }
        }
        else
        {
            st->TransientHangOver = 1;
        }
    }
    else
    {
        if( st->TransientHangOver )
        {
            st->TransientHangOver = 0;
            IsTransient = 1;
        }
    }

    st->Energy_Old = Energy_in[4];

    return IsTransient;
}
