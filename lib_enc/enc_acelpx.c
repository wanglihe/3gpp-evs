/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <memory.h>
#include <assert.h>
#include "typedef.h"
#include "options.h"
#include "prot.h"
#include "rom_enc.h"



/* Iterations: nb_pos_ix*16 */
static
void E_ACELP_2pulse_searchx(Word32 nb_pos_ix, UWord8 track_x,
                            UWord8 track_y, Float32 *R, Float32 *ps, Float32 *alp,
                            Word16 *ix, Word16 *iy, Float32 dn[],
                            Word16 *dn2, Float32 cor[], Float32 sign[])
{
    Word32 i,y, x_save=0, y_save=0;
    Word16 x, *pos_x;
    Float32 ps0, alp0, alp1, ps1, alp2, ps2, sq, s, sqk, alpk, *pR, sgnx, *pRx, *pRy, sign_x, sign_y;
    /* x_save=y_save=0 */
    /* eight dn2 max positions per track */
    pos_x = &dn2[track_x << 3];
    /* save these to limit memory searches */
    ps0 = *ps;
    alp0 = *alp + 2.0f*R[0];

    sqk = -1.0F;
    alpk = 1.0F;

    x = pos_x[0];
    sgnx = sign[track_y];
    if (sign[x] < 0)
    {
        sgnx = -sgnx;
    }
    if ( (alp0 + (cor[x] * sign[x]) + (cor[track_y] * sign[track_y]) + (R[track_y-x] * sgnx)) < 0.0F )
    {
        sqk = 1.0F;
    }


    /* loop track 1 */
    for (i=0; i<nb_pos_ix; i++)
    {
        x = pos_x[i];
        sgnx = sign[x];
        /* dn[x] has only nb_pos_ix positions saved */
        ps1 = ps0 + dn[x];
        alp1 = alp0 + 2*sgnx*cor[x];
        pR = R-x;

        for (y = track_y; y < L_SUBFR; y += 4)
        {
            ps2 = ps1 + dn[y];
            alp2 = alp1 + 2.0f*sign[y]*(cor[y] + sgnx*pR[y]);
            sq = ps2 * ps2;

            s = (alpk * sq) - (sqk * alp2);
            if (s > 0.0F)
            {
                sqk = sq;
                alpk = alp2;
                y_save = y;
                x_save = x;
            }
        }
    }
    /* Update numerator */
    *ps = ps0 + dn[x_save] + dn[y_save];
    /* Update denominator */
    *alp = alpk;

    /* Update product of autocorrelation and already fixed pulses. with the
     * two newly found ones */
    pRx = R-x_save;
    pRy = R-y_save;
    sign_x = sign[x_save];
    sign_y = sign[y_save];
    for (i=0; i<L_SUBFR; i++)
    {
        cor[i] += pRx[i]*sign_x + pRy[i]*sign_y;
    }


    *ix = x_save;
    *iy = y_save;
    if (((x_save & 3) != track_x) || ((y_save & 3) != track_y))
    {
        /* sanity check */
        assert(0);
    }
    return;
}


static
void E_ACELP_1pulse_searchx(UWord8 track_x,
                            UWord8 track_y, Float32 *R, Float32 *ps, Float32 *alp,
                            Word16 *ix, Float32 dn[],
                            Float32 cor[], Float32 sign[])
{
    Word32 x, x_save = 0;
    Float32 ps0, alp0;
    Float32 ps1, sq, sqk;
    Float32 alp1, alpk;
    Float32 s;

    /* save these to limit memory searches */
    ps0 = *ps;
    alp0 = *alp + R[0];
    sqk = -1.0F;
    alpk = 1.0F;

    if ( (alp0 + (cor[track_x] * sign[track_x])) < 0 )
    {
        sqk = 1.0F;
    }

    x_save = track_x;
    for (x = track_x; x < L_SUBFR; x += 4)
    {
        ps1 = ps0 + dn[x];
        alp1 = alp0 + 2*sign[x]*cor[x];
        sq = ps1 * ps1;
        s = (alpk * sq) - (sqk * alp1);
        if (s > 0.0F)
        {
            sqk = sq;
            alpk = alp1;
            x_save = x;
        }
    }
    if(track_y!=track_x)
    {
        for (x = track_y; x < L_SUBFR; x += 4)
        {
            ps1 = ps0 + dn[x];
            alp1 = alp0 + 2*sign[x]*cor[x];
            sq = ps1 * ps1;
            s = (alpk * sq) - (sqk * alp1);
            if (s > 0.0F)
            {
                sqk = sq;
                alpk = alp1;
                x_save = x;
            }
        }
    }

    *ps = ps0 + dn[x_save];
    *alp = alpk;
    *ix = x_save;
    return;
}


/* Autocorrelation method for searching pulse positions effectively
 * Algorithm is identical to traditional covariance method. */
void E_ACELP_4tsearchx(
    Float32 dn[],
    const Float32 cn[],
    Float32 Rw[],
    float code[],
    PulseConfig *config
    , Word16 ind[]
)
{
    Float32 sign[L_SUBFR], vec[L_SUBFR];
    Float32 cor[L_SUBFR];
    Float32 R_buf[2*L_SUBFR-1], *R;
    Float32 dn2[L_SUBFR];
    Float32 psk = 0.0F, ps2k, ps, ps2, alpk, alp = 0.0F;
    Word16 codvec[NB_PULSE_MAX];
    Word16 pos_max[4];
    Word16 dn2_pos[8 * 4];
    UWord8 ipos[NB_PULSE_MAX];
    Float32 *p0;
    Word32 i, j, k, l, st, pos = 0, index, track;
    UWord8 iPulse;
    Float32 val;
    Float32 s;
    UWord8 restpulses;


    alp = config->alp;
    for (k=0; k<config->nb_pulse; k++)
    {
        codvec[k] = (k & 3);
    }

    set_f( cor, 0.0f, L_SUBFR);

    /* Set up autocorrelation vector */
    R = R_buf+L_SUBFR-1;
    R[0] = Rw[0];
    for (k=1; k<L_SUBFR; k++)
    {
        R[k] = R[-k] = Rw[k];
    }



    /*
     * Find sign for each pulse position.
     */

    acelp_pulsesign(cn, dn, dn2, sign, vec, alp);

    /*
     * Select the most important 8 position per track according to dn2[].
     */
    acelp_findcandidates(dn2, dn2_pos, pos_max, L_SUBFR, NB_TRACK_FCB_4T);



    /*
     * Deep first search:
     * ------------------
     * 20 bits (4p):  4 iter x ((4x16)+(8x16))              = 768 tests
     * 36 bits (8p):  4 iter x ((1x1)+(4x16)+(8x16)+(8x16)) = 1280 tests
     * 52 bits (12p): 3 iter x ((1x1)+(1x1)+(4x16)+(6x16)
     *                                      +(8x16)+(8x16)) = 1248 tests
     * 64 bits (16p): 2 iter x ((1x1)+(1x1)+(4x16)+(6x16)
     *                        +(6x16)+(8x16)+(8x16)+(8x16)) = 1280 tests
     */
    ps2k = -1.0;
    alpk = 1000.0;
    /*Number of iterations*/
    for (k = 0; k < config->nbiter; k++)
    {
        /* copy search order from hash-table */
        assert(config->nb_pulse+(k*4) <= 40);
        for ( l=0; l<config->nb_pulse; l++ )
        {
            ipos[l] = tipos[(k * 4) + l];
        }

        /* if all tracks do not have equal number of pulses */
        restpulses = config->nb_pulse & 3;
        if (restpulses)
        {
            switch (config->codetrackpos)
            {
            case TRACKPOS_FIXED_FIRST:  /* fixed track positions, starting from left */
                /* add tracks from left */
                for (iPulse=0; iPulse < restpulses; iPulse++)
                {
                    ipos[config->nb_pulse-restpulses+iPulse] = iPulse;
                }
                /* Put the same track on the next position, because the 1-pulse search
                 * will access it to determine if this could be in any track. */
                ipos[config->nb_pulse] = ipos[config->nb_pulse-1];
                break;
            case TRACKPOS_FIXED_EVEN:  /* fixed track positions, odd tracks */
                /* odd tracks, switch order for every iteration */
                ipos[config->nb_pulse-restpulses] = (k<<1) & 2;    /* 0 for even k, 2 for odd*/
                ipos[config->nb_pulse-restpulses+1] = ipos[config->nb_pulse-restpulses] ^ 2; /* 2 for even k, 0 for odd*/
                break;
            case TRACKPOS_FIXED_TWO:  /* two tracks instead of four */
                /* Put the next track on the next position, because the 1-pulse search
                 * will access it to determine if this could be in any track. */
                ipos[config->nb_pulse] = (ipos[config->nb_pulse-1]+1) & 3;
                break;
            default:              /* one or three free track positions */
                /* copy an extra position from table - 1pulse search will access this */
                ipos[config->nb_pulse] = tipos[(k*4)+config->nb_pulse];
                break;
            }
        }
        if (config->fixedpulses == 0)/* 1100, 11, 1110, 1111, 2211 */
        {
            pos = 0;
            ps = 0.0F;
            alp = 0.0F;
            memset(cor, 0, L_SUBFR * sizeof(Float32));
        }
        else if (config->fixedpulses == 2) /* 2222 and 3322 */
        {
            /* --- first stage: fix 2 pulses --- */
            /* index to first non-fixed position */
            pos = 2;

            /* set fixed positions */
            ind[0] = pos_max[ipos[0]];
            ind[1] = pos_max[ipos[1]];

            /* correlation of fixed part with residual */
            ps = dn[ind[0]] + dn[ind[1]];

            /* multiplication of autocorrelation with signed fixed pulses */
            /* first pulse */
            p0 = R-ind[0];
            if (sign[ind[0]]>0)
            {
                for (i=0; i<L_SUBFR; i++)
                {
                    cor[i] = *p0;
                    p0++;
                }
            }
            else
            {
                for (i=0; i<L_SUBFR; i++)
                {
                    cor[i] = -*p0;
                    p0++;
                }
            }
            /* second pulse */
            p0 = R-ind[1];
            if (sign[ind[1]]>0)
            {
                for (i=0; i<L_SUBFR; i++)
                {
                    cor[i] += *p0;
                    p0++;
                }
            }
            else
            {
                for (i=0; i<L_SUBFR; i++)
                {
                    cor[i] -= *p0;
                    p0++;
                }
            }

            /* normalisation contribution of fixed part */
            alp = sign[ind[0]]*cor[ind[0]] + sign[ind[1]]*cor[ind[1]];


        }
        else /* if (config->fixedpulses == 4) */ /* 3333 and above */
        {
            /* first stage: fix 4 pulses */
            pos = 4;

            ind[0] = pos_max[ipos[0]];
            ind[1] = pos_max[ipos[1]];
            ind[2] = pos_max[ipos[2]];
            ind[3] = pos_max[ipos[3]];

            /* correlation of fixed part with residual */
            ps = dn[ind[0]] + dn[ind[1]] + dn[ind[2]] + dn[ind[3]];

            /* multiplication of autocorrelation with signed fixed pulses */
            /* first pulse */
            p0 = R-ind[0];
            if (sign[ind[0]]>0)
            {
                for (i=0; i<L_SUBFR; i++)
                {
                    cor[i] = *p0;
                    p0++;
                }
            }
            else
            {
                for (i=0; i<L_SUBFR; i++)
                {
                    cor[i] = -*p0;
                    p0++;
                }
            }
            /* pulses 1..3 */
            for (j=1; j<4; j++)
            {
                p0 = R-ind[j];
                if (sign[ind[j]]>0)
                {
                    for (i=0; i<L_SUBFR; i++)
                    {
                        cor[i] += *p0;
                        p0++;
                    }
                }
                else
                {
                    for (i=0; i<L_SUBFR; i++)
                    {
                        cor[i] -= *p0;
                        p0++;
                    }
                }
            }


            /* normalisation contribution of fixed part */
            alp = sign[ind[0]]*cor[ind[0]]
                  + sign[ind[1]]*cor[ind[1]]
                  + sign[ind[2]]*cor[ind[2]]
                  + sign[ind[3]]*cor[ind[3]];

        }

        /* other stages of 2 pulses */
        for (j = pos, st = 0; j < config->nb_pulse; j += 2, st++)
        {
            if ((config->nb_pulse-j) >= 2)
            {
                /*pair-wise search*/

                /*
                 * Calculate correlation of all possible positions
                 * of the next 2 pulses with previous fixed pulses.
                 * Each pulse can have 16 possible positions.
                 */

                E_ACELP_2pulse_searchx(config->nbpos[st], ipos[j], ipos[j + 1], R, &ps, &alp,&ind[j], &ind[j+1], dn, dn2_pos, cor, sign);

            }
            else
            {
                /*single pulse search*/
                E_ACELP_1pulse_searchx(ipos[j], ipos[j + 1], R, &ps, &alp,&ind[j], dn, cor, sign);
            }
        }

        /* memorise the best codevector */
        ps2 = ps * ps;
        s = (alpk * ps2) - (ps2k * alp);
        if (s > 0.0F)
        {
            ps2k = ps2;
            psk = ps;
            alpk = alp;
            memcpy(codvec, ind, config->nb_pulse * sizeof(Word16));
        }

    }

    /*
     * Store weighted energy of code, build the codeword and index of codevector.
     */

    memset(code, 0,  L_SUBFR * sizeof(float));
    memset(ind, 0xffffffff, NPMAXPT * 4 * sizeof(Word16));
    for (k = 0; k < config->nb_pulse; k++)
    {
        i = codvec[k];  /* read pulse position  */
        val = sign[i];  /* read sign            */

        index = i / 4;  /* pos of pulse (0..15) */
        track = i % 4;
        if (val*psk > 0)
        {
            code[i] += 1.0f;
            codvec[k] += (2 * L_SUBFR);
        }
        else
        {
            code[i] -= 1.0f;
            index += 16;
        }

        i = track * NPMAXPT;
        while (ind[i] >= 0)
        {
            i++;
        }

        ind[i] = index;

    }

    return;

}
