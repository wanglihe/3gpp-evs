/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <memory.h>
#include <assert.h>
#include "typedef.h"
#include "prot.h"
#include "rom_com.h"


/*---------------------------------------------------------------------*
* Local functions
*---------------------------------------------------------------------*/

static void D_ACELP_decode_arithtrack( float v[], long unsigned s, int p, int trackstep, int tracklen );

/*---------------------------------------------------------------------*
* Function D_ACELP_indexing()
*
*---------------------------------------------------------------------*/

void D_ACELP_indexing(
    Float32 code[],
    PulseConfig config,
    int num_tracks,
    int index[]
    ,short *BER_detect
)
{
    int track, pulses, k, pulsestrack[NB_TRACK_FCB_4T];
    long unsigned s;
    long unsigned index_n[NB_TRACK_FCB_4T];
    unsigned short trackpos, idxs[MAX_IDX_LEN];
    int restpulses, wordcnt;

    assert(num_tracks == NB_TRACK_FCB_4T);

    wordcnt = (config.bits + 15) >> 4;          /* ceil(bits/16) */

    /* check if some tracks have more pulses */
    restpulses = config.nb_pulse & (num_tracks-1);

    /* cast to short */
    for (k=0; k<wordcnt; k++)
    {
        idxs[k] = (unsigned short)index[k];
    }

    if (restpulses)
    {
        /* check if we need to code track positions */
        switch (config.codetrackpos)
        {
        case TRACKPOS_FREE_THREE:
            /* Find track with less pulses */
            trackpos = idxs[0] & 3;
            longshiftright(idxs,2,idxs,wordcnt,wordcnt);

            /* set number of pulses per track */
            set_i( pulsestrack,(config.nb_pulse>>2)+1,4);
            pulsestrack[trackpos]--;    /* this one has less pulses */
            break;
        case TRACKPOS_FREE_ONE:
            /* Find track with more pulses */
            trackpos = idxs[0] & 3;
            longshiftright(idxs,2,idxs,wordcnt,wordcnt);

            /* set number of pulses per track */
            set_i( pulsestrack,(config.nb_pulse>>2),4);
            pulsestrack[trackpos]++;    /* this one has more pulses */
            break;
        case TRACKPOS_FIXED_EVEN:
            /* Pulses on even tracks */
            pulsestrack[0] = (config.nb_pulse+1) >> 1;
            pulsestrack[1] = 0;
            pulsestrack[2] = config.nb_pulse >> 1;
            pulsestrack[3] = 0;
            break;
        case TRACKPOS_FIXED_FIRST:
            /* set number of pulses per track */
            set_i( pulsestrack,config.nb_pulse/num_tracks,4);
            for (k=0; k<restpulses; k++)
            {
                pulsestrack[k]++;
            }
            break;
        case TRACKPOS_FIXED_TWO:
            /* 1100, 0110, 0011, 1001 */
            /* Find track with less pulses */
            trackpos = idxs[0] & 3;
            longshiftright(idxs,2,idxs,wordcnt,wordcnt);

            /* set number of pulses per track */
            set_i( pulsestrack,(config.nb_pulse>>2),4);
            pulsestrack[trackpos]++;
            trackpos++;
            trackpos &= 3;
            pulsestrack[trackpos]++;
            break;
        default:
            assert(0);
            break;
        }
    }
    else
    {
        /* set number of pulses per track */
        set_i( pulsestrack,(config.nb_pulse/num_tracks),num_tracks);
    }

    if (config.bits == 43)
    {
        D_ACELP_decode_43bit(idxs, code, pulsestrack);
    }
    else
    {

        fcb_pulse_track_joint_decode(idxs, wordcnt, index_n, pulsestrack, num_tracks);

        for (track=num_tracks-1; track >=1; track--)
        {
            pulses = pulsestrack[track];

            if (pulses)
            {
                s = index_n[track];

                /* decode state to actual pulse positions on track */
                /*D_ACELP_decode_arithtrack_old(code+track, s, pulses, 4);                    */
                D_ACELP_decode_arithtrack(code+track, s, pulses, num_tracks, 16);
            }
            else
            {
                /* track is empty */
                for (k=track; k < 16*num_tracks; k+=num_tracks)
                {
                    code[k] = 0.0f;
                }
            }
        }
        s = index_n[0];
        pulses = pulsestrack[0];
        /* safety check in case of bit errors */
        if (s >= pulsestostates[16][pulses-1])
        {
            set_f( code, 0.0f, L_SUBFR );
            *BER_detect = 1;
            return;
        }
        if (pulses)
        {
            D_ACELP_decode_arithtrack(code, s, pulses, num_tracks, 16);
        }
        else
        {
            /* track is empty */
            for (k=0; k < 16*num_tracks; k+=num_tracks)
            {
                code[k] = 0.0f;
            }
        }
    }

    return;
}


static void D_ACELP_decode_arithtrack(
    float v[],
    long unsigned s,
    int p,
    int trackstep,
    int tracklen
)
{
    int k;

    for (k=(tracklen)-1; k>= 0; k--)
    {
        v[k*trackstep] = 0.0f;          /* default: there is no pulse here */
        while((p) && (s >= pulsestostates[k][p-1] ))
        {
            s -= pulsestostates[k][p-1];
            if (v[k*trackstep])
            {
                /* there is a pulse here already = sign is known */
                if (v[k*trackstep] > 0.0f)
                {
                    v[k*trackstep]++; /* place one more pulse here */
                }
                else
                {
                    v[k*trackstep]--; /* place one more pulse here */
                }
            }
            else
            {
                /* this is the first pulse here -> determine sign */
                if (s & 1)
                {
                    v[k*trackstep] = -1.0f; /* place a negative pulse here */
                }
                else
                {
                    v[k*trackstep] = +1.0f; /* place a negative pulse here */
                }
                s >>= 1;
            }
            p--;                /* one pulse placed, so one less left */
        }
    }

    return;
}


void fcb_pulse_track_joint_decode(unsigned short *idxs, int wordcnt, long unsigned *index_n, int *pulse_num, int track_num)
{
    int hi_to_low[10] = { 0, 0, 0, 3, 9, 5, 3, 1, 8, 8};

    unsigned long long index;
    int indx_tmp,indx_flag,indx_flag_1;
    int track,track_num1,pulse_num0,pulse_num1;
    int div_tmp;
    int indx_flag_2;

    indx_flag=0;
    indx_flag_1=0;
    indx_flag_2 = 0;
    for (track=0; track < track_num; track++)
    {
        indx_flag += (pulse_num[track]>>2);
        indx_flag_1 += (pulse_num[track]>>1);
        indx_flag_2 += (pulse_num[track]>>3);
    }

    if (indx_flag >= track_num)
    {
        hi_to_low[4] = 9;
    }
    else
    {
        hi_to_low[4] = 1;
    }

    if (indx_flag_2 >= 1)
    {
        hi_to_low[7] = 9;
    }
    else
    {
        hi_to_low[7] = 1;
    }
    if (indx_flag_1>=track_num)
    {
        if (indx_flag>=track_num)
        {
            index = 0;
            if (indx_flag_2 >= 1)
            {
                for (track=(wordcnt-1); track >= 6; track--)
                {
                    index = ( index << 16 ) + idxs[track] ;
                }
                index_n[3] = ( ((unsigned int)idxs[5]) << 8 ) + ( ( idxs[4] >> 8 ) & 0xff );
                index_n[2] = ( ( ((unsigned int)idxs[4]) << 16 ) + idxs[3] ) & 0xffffffUL;
                index_n[1] = ( ((unsigned int)idxs[2]) << 8 ) + ( ( idxs[1] >> 8 ) & 0xff );
                index_n[0] = ( ( ((unsigned int) idxs[1]) << 16 ) + idxs[0] ) & 0xffffffUL;
            }
            else
            {
                for (track=(wordcnt-1); track >= track_num; track--)
                {
                    index = ( index << 16 ) + idxs[track] ;
                }
                for (track=0; track < track_num; track++)
                {
                    index_n[track] = idxs[track];
                }
            }
        }
        else
        {
            index = 0;
            for (track=(wordcnt-1); track >= 2; track--)
            {
                index = ( index << 16 ) + idxs[track] ;
            }

            index_n[3] = idxs[1] & 0xff;
            index_n[2] = idxs[1] >> 8;
            index_n[1] = idxs[0] & 0xff;
            index_n[0] = idxs[0] >> 8;
        }

        track_num1 = track_num - 1;
        pulse_num1 = pulse_num[track_num1];
        index = ( index << hi_to_low[pulse_num1] ) + ( index_n[track_num1] >> low_len[pulse_num1] );
        for (track=(track_num-1); track > 0; track--)
        {
            track_num1 = track - 1;
            pulse_num0 = pulse_num[track_num1];
            pulse_num1 = pulse_num[track];
            index = ( index << hi_to_low[pulse_num0] ) + ( index_n[track_num1] >> low_len[pulse_num0] );

            div_tmp = index / indx_fact[pulse_num1];
            indx_tmp = index - div_tmp * indx_fact[pulse_num1];
            index_n[track] = ( index_n[track] & low_mask[pulse_num1]) + ( indx_tmp << low_len[pulse_num1] );
            index = div_tmp;
        }
        pulse_num1 = pulse_num[0];
        index_n[0] = ( index_n[0] & low_mask[pulse_num1]) + ( index << low_len[pulse_num1] );
    }
    else
    {
        index = 0;
        for (track=(wordcnt-1); track >= 0; track--)
        {
            index = ( index << 16 ) + idxs[track];
        }
        for (track=3; track > 0; track--)
        {
            pulse_num1 = pulse_num[track];
            index_n[track] =  index & index_mask[pulse_num1];
            index = index >> index_len[pulse_num1];
        }
        index_n[0] = index;
    }

    return;
}
