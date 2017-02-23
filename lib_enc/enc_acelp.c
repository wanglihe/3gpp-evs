/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <memory.h>
#include <assert.h>
#include "typedef.h"
#include "options.h"
#include "prot.h"
#include "rom_com.h"
#include "rom_enc.h"



/*---------------------------------------------------------------------*
* Local functions
*---------------------------------------------------------------------*/

static void E_ACELP_codearithp(const float v[], long unsigned *n, long unsigned *ps, int *p, int trackstep, int tracklen);

/*---------------------------------------------------------------------*
* Local constants
*---------------------------------------------------------------------*/

#define NB_MAX          8

/*
 * E_ACELP_h_vec_corrx
 *
 * Parameters:
 *    h              I: scaled impulse response
 *    vec            I: vector to correlate with h[]
 *    track          I: track to use
 *    sign           I: sign vector
 *    rrixix         I: correlation of h[x] with h[x]
 *    cor            O: result of correlation (16 elements)
 *
 * Function:
 *    Calculate the correlations of h[] with vec[] for the specified track
 *
 * Returns:
 *    void
 */
static void acelp_h_vec_corr1(
    Float32 h[],
    Float32 vec[],
    UWord8 track,
    Float32 sign[],
    Float32 (*rrixix)[16],
    Float32 cor[],
    Word16 dn2_pos[],
    Word32 nb_pulse
)
{
    Word16 i, j;
    Word32 dn;
    Word16 *dn2;
    Float32 *p0;
    Float32 s;

    dn2 = &dn2_pos[track * 8];
    p0 = rrixix[track];
    for (i = 0; i < nb_pulse; i++)
    {
        dn = dn2[i];
        s = 0.0F;
        /* L_SUBFR-dn */
        /* vec[dn] */
        for (j = 0; j < (L_SUBFR - dn); j++)
        {
            s += h[j] * vec[dn + j];
        }

        cor[dn >> 2] = sign[dn] * s + p0[dn >> 2];
    }
    return;
}

static void acelp_h_vec_corr2(
    Float32 h[],
    Float32 vec[],
    UWord8 track,
    Float32 sign[],
    Float32 (*rrixix)[16],
    Float32 cor[])
{
    Word32 i, j;
    Float32 *p0;
    Float32 s;


    p0 = rrixix[track];
    /* sign[track] */
    for (i = 0; i < 16; i++)
    {
        s = 0.0F;
        /* h[0], vec[track] */
        /* L_SUBFR-track */
        for (j = 0; j < L_SUBFR - track; j++)
        {
            s += h[j] * vec[track + j];
        }

        cor[i] = s * sign[track] + p0[i];
        track += 4;
    }
    return;
}


/*
 * acelp_2pulse_search
 *
 * Parameters:
 *    nb_pos_ix      I: nb of pos for pulse 1 (1..8)
 *    track_x        I: track of pulse 1
 *    track_y        I: track of pulse 2
 *    ps           I/O: correlation of all fixed pulses
 *    alp          I/O: energy of all fixed pulses
 *    ix             O: position of pulse 1
 *    iy             O: position of pulse 2
 *    dn             I: corr. between target and h[]
 *    dn2            I: vector of selected positions
 *    cor_x          I: corr. of pulse 1 with fixed pulses
 *    cor_y          I: corr. of pulse 2 with fixed pulses
 *    rrixiy         I: corr. of pulse 1 with pulse 2
 *
 * Function:
 *    Find the best positions of 2 pulses in a subframe
 *
 * Returns:
 *    void
 */
static void acelp_2pulse_search(  Word32 nb_pos_ix, UWord8 track_x,
                                  UWord8 track_y, Float32 *ps, Float32 *alp,
                                  Word16 *ix, Word16 *iy, Float32 dn[],
                                  Word16 *dn2, Float32 cor_x[],
                                  Float32 cor_y[], Float32 (*rrixiy)[256])
{
    Word16 x, x2, y, x_save = 0, y_save = 0, i, *pos_x;
    Float32 ps0, alp0;
    Float32 ps1, ps2, sq, sqk;
    Float32 alp1, alp2, alpk;
    Float32 *p1, *p2;
    Float32 s;

    /* x_save=y_save=0 */
    /* eight dn2 max positions per track */
    pos_x = &dn2[track_x << 3];
    /* save these to limit memory searches */
    ps0 = *ps;
    alp0 = *alp;

    alpk = 1.0F;
    sqk = -1.0F;
    x2 = pos_x[0] >> 2;
    if( (alp0 + cor_x[x2] + cor_y[0] + rrixiy[track_x][x2 << 4]) < 0 )
    {
        sqk = 1.0F;
    }

    /* loop track 1 */
    for (i = 0; i < nb_pos_ix; i++)
    {
        x = pos_x[i];
        x2 = x >> 2;
        /* dn[x] has only nb_pos_ix positions saved */
        ps1 = ps0 + dn[x];
        alp1 = alp0 + cor_x[x2];
        p1 = cor_y;
        p2 = &rrixiy[track_x][x2 << 4];
        for (y = track_y; y < L_SUBFR; y += 4)
        {
            ps2 = ps1 + dn[y];
            alp2 = alp1 + (*p1++) + (*p2++);

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

    *ps = ps0 + dn[x_save] + dn[y_save];
    *alp = alpk;
    *ix = x_save;
    *iy = y_save;
    return;
}


/*
 * E_ACELP_1pulse_search
 *
 * Parameters:
 *    track_x        I: track of pulse 1
 *    track_y        I: track of pulse 2
 *    ps           I/O: correlation of all fixed pulses
 *    alp          I/O: energy of all fixed pulses
 *    ix             O: position of pulse 1
 *    dn             I: corr. between target and h[]
 *    cor_x          I: corr. of pulse 1 with fixed pulses
 *    cor_y          I: corr. of pulse 2 with fixed pulses
 *
 * Function:
 *    Find the best positions of 1 pulse in a subframe
 *
 * Returns:
 *    void
 */
static void E_ACELP_1pulse_search(UWord8 track_x,
                                  UWord8 track_y,
                                  Float32 *ps,
                                  Float32 *alp,
                                  Word16 *ix,
                                  Float32 dn[],
                                  Float32 cor_x[],
                                  Float32 cor_y[])
{
    Word16 x, x_save = 0;
    Float32 ps0, alp0;
    Float32 ps1, sq, sqk;
    Float32 alp1, alpk;
    Float32 s;

    /* save these to limit memory searches */
    ps0 = *ps;
    alp0 = *alp;
    alpk = 1.0F;
    sqk = -1.0F;

    if ( (alp0 + cor_x[(track_x >> 2)]) < 0)
    {
        sqk = 1.0F;
    }
    for (x = track_x; x < L_SUBFR; x += 4)
    {
        ps1 = ps0 + dn[x];
        alp1 = alp0 + cor_x[x>>2];
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
            alp1 = alp0 + cor_y[x>>2];
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


/*
 * acelp_pulsesign
 *
 * Parameters:
 *    cn          I: residual after Word32 term prediction
 *    dn          I: corr. between target and h[].
 *    dn2         O: dn2[] = mix of dn[] and cn[]
 *    sign        O: sign of pulse
 *    vec         O: negative sign of pulse
 *
 * Function:
 *    Determine sign of each pulse position, store them in "sign"
 *    and change dn to all positive.
 *    Subframe size = L_SUBFR
 * Returns:
 *    void
 */
void acelp_pulsesign(
    const float cn[],
    float dn[],
    float dn2[],
    float sign[],
    float vec[],
    float alp
)
{
    int i;
    float val;
    float s, cor;

    /* calculate energy for normalization of cn[] and dn[] */
    val = (cn[0] * cn[0]) + 1.0F;
    cor = (dn[0] * dn[0]) + 1.0F;
    for (i = 1; i < L_SUBFR; i ++)
    {
        val += (cn[i] * cn[i]);
        cor += (dn[i] * dn[i]);
    }

    s = (float)sqrt(cor / val);
    for (i = 0; i < L_SUBFR; i++)
    {
        cor = (s * cn[i]) + (alp * dn[i]);
        if (cor >= 0.0F)
        {
            sign[i] = 1.0F;
            vec[i] = -1.0F;
            dn2[i] = cor;     /* dn2[] = mix of dn[] and cn[]   */
        }
        else
        {
            sign[i] = -1.0F;
            vec[i] = 1.0F;
            dn[i] = -dn[i];   /* modify dn[] according to the fixed sign */
            dn2[i] = -cor;    /* dn2[] = mix of dn[] and cn[]            */
        }
    }
    return;
}

void acelp_findcandidates(float dn2[], short dn2_pos[], short pos_max[], int L_subfr, int tracks)
{
    int i,k,j;
    float *ps_ptr;
    /* &pos_max[0], &dn2_pos[0] */
    for (i = 0; i < tracks; i++)
    {
        for (k = 0; k < NB_MAX; k++)
        {
            ps_ptr = &dn2[i];
            for (j = i+tracks; j < L_subfr; j += tracks)
            {
                if (dn2[j] > *ps_ptr)
                {
                    ps_ptr = &dn2[j];
                }
            }
            *ps_ptr = (float)k - NB_MAX;    /* dn2 < 0 when position is selected */
            dn2_pos[i * 8 + k] = ps_ptr - dn2;
        }
        pos_max[i] = dn2_pos[i * 8];
    }
}

static void acelp_hbuf(float *h_buf, float **h, float **h_inv, const float *H)
{
    int i;
    *h = h_buf + L_SUBFR;
    *h_inv = h_buf + (3*L_SUBFR);
    for (i=0; i<L_SUBFR; i++)
    {
        (*h)[-1-i] = 0.0f;
        (*h_inv)[-1-i] = 0.0f;
        (*h)[i] = H[i];
        (*h_inv)[i] = -H[i];
    }
}


static void E_ACELP_corrmatrix(Float32 h[], Float32 sign[], Float32 vec[], Float32 rrixix[4][16], Float32 rrixiy[4][256])
{
    Float32 *p0, *p1, *p2, *p3, *psign0, *psign1, *psign2, *psign3;
    Float32 *ptr_h1, *ptr_h2, *ptr_hf;
    Float32 cor;
    int i,k,pos;
    /*
    * Compute rrixix[][] needed for the codebook search.
    */

    /* storage order --> i3i3, i2i2, i1i1, i0i0 */

    /* Init pointers to last position of rrixix[] */
    p0 = &rrixix[0][16 - 1];
    p1 = &rrixix[1][16 - 1];
    p2 = &rrixix[2][16 - 1];
    p3 = &rrixix[3][16 - 1];

    ptr_h1 = h;
    cor    = 0.0F;
    for(i = 0; i < 16; i++)
    {
        cor += (*ptr_h1) * (*ptr_h1);
        ptr_h1++;
        *p3-- = cor * 0.5F;
        cor += (*ptr_h1) * (*ptr_h1);
        ptr_h1++;
        *p2-- = cor * 0.5F;
        cor += (*ptr_h1) * (*ptr_h1);
        ptr_h1++;
        *p1-- = cor * 0.5F;
        cor += (*ptr_h1) * (*ptr_h1);
        ptr_h1++;
        *p0-- = cor * 0.5F;
    }


    /*
     * Compute rrixiy[][] needed for the codebook search.
     */

    /* storage order --> i2i3, i1i2, i0i1, i3i0 */

    pos = 256 - 1;
    ptr_hf = h + 1;
    for(k = 0; k < 16; k++)
    {

        p3 = &rrixiy[2][pos];
        p2 = &rrixiy[1][pos];
        p1 = &rrixiy[0][pos];
        p0 = &rrixiy[3][pos - 16];

        cor = 0.0F;
        ptr_h1 = h;
        ptr_h2 = ptr_hf;
        for(i = k; i < 15; i++)
        {
            cor += (*ptr_h1) * (*ptr_h2);
            ptr_h1++;
            ptr_h2++;
            *p3 = cor;
            cor += (*ptr_h1) * (*ptr_h2);
            ptr_h1++;
            ptr_h2++;
            *p2 = cor;
            cor += (*ptr_h1) * (*ptr_h2);
            ptr_h1++;
            ptr_h2++;
            *p1 = cor;
            cor += (*ptr_h1) * (*ptr_h2);
            ptr_h1++;
            ptr_h2++;
            *p0 = cor;

            p3 -= (16 + 1);
            p2 -= (16 + 1);
            p1 -= (16 + 1);
            p0 -= (16 + 1);
        }

        cor += (*ptr_h1) * (*ptr_h2);
        ptr_h1++;
        ptr_h2++;
        *p3 = cor;
        cor += (*ptr_h1) * (*ptr_h2);
        ptr_h1++;
        ptr_h2++;
        *p2 = cor;
        cor += (*ptr_h1) * (*ptr_h2);
        ptr_h1++;
        ptr_h2++;
        *p1 = cor;

        pos -= 16;
        ptr_hf += 4;
    }

    /* storage order --> i3i0, i2i3, i1i2, i0i1 */

    pos = 256 - 1;
    ptr_hf = h + 3;
    for(k = 0; k < 16; k++)
    {

        p3 = &rrixiy[3][pos];
        p2 = &rrixiy[2][pos - 1];
        p1 = &rrixiy[1][pos - 1];
        p0 = &rrixiy[0][pos - 1];

        cor = 0.0F;
        ptr_h1 = h;
        ptr_h2 = ptr_hf;
        for(i= k + 1; i < 16; i++ )
        {
            cor += (*ptr_h1) * (*ptr_h2);
            ptr_h1++;
            ptr_h2++;
            *p3 = cor;
            cor += (*ptr_h1) * (*ptr_h2);
            ptr_h1++;
            ptr_h2++;
            *p2 = cor;
            cor += (*ptr_h1) * (*ptr_h2);
            ptr_h1++;
            ptr_h2++;
            *p1 = cor;
            cor += (*ptr_h1) * (*ptr_h2);
            ptr_h1++;
            ptr_h2++;
            *p0 = cor;

            p3 -= (16 + 1);
            p2 -= (16 + 1);
            p1 -= (16 + 1);
            p0 -= (16 + 1);
        }

        cor += (*ptr_h1) * (*ptr_h2);
        *p3 = cor;

        pos--;
        ptr_hf += 4;
    }

    /*
     * Modification of rrixiy[][] to take signs into account.
     */

    p0 = &rrixiy[0][0];
    /* speed-up: 11% */
    p1 = &rrixiy[1][0];
    p2 = &rrixiy[2][0];
    p3 = &rrixiy[3][0];
    for(i = 0; i < L_SUBFR; i += 4)
    {
        if (sign[i+0] < 0.0F) psign0 = &vec[1];
        else                  psign0 = &sign[1];
        if (sign[i+1] < 0.0F) psign1 = &vec[2];
        else                  psign1 = &sign[2];
        if (sign[i+2] < 0.0F) psign2 = &vec[3];
        else                  psign2 = &sign[3];
        if (sign[i+3] < 0.0F) psign3 = &vec[0];
        else                  psign3 = &sign[0];
        p0[0]  = p0[0]  * psign0[ 0];
        p0[1]  = p0[1]  * psign0[ 4];
        p0[2]  = p0[2]  * psign0[ 8];
        p0[3]  = p0[3]  * psign0[12];
        p0[4]  = p0[4]  * psign0[16];
        p0[5]  = p0[5]  * psign0[20];
        p0[6]  = p0[6]  * psign0[24];
        p0[7]  = p0[7]  * psign0[28];
        p0[8]  = p0[8]  * psign0[32];
        p0[9]  = p0[9]  * psign0[36];
        p0[10] = p0[10] * psign0[40];
        p0[11] = p0[11] * psign0[44];
        p0[12] = p0[12] * psign0[48];
        p0[13] = p0[13] * psign0[52];
        p0[14] = p0[14] * psign0[56];
        p0[15] = p0[15] * psign0[60];
        p0 += 16;

        p1[0]  = p1[0]  * psign1[ 0];
        p1[1]  = p1[1]  * psign1[ 4];
        p1[2]  = p1[2]  * psign1[ 8];
        p1[3]  = p1[3]  * psign1[12];
        p1[4]  = p1[4]  * psign1[16];
        p1[5]  = p1[5]  * psign1[20];
        p1[6]  = p1[6]  * psign1[24];
        p1[7]  = p1[7]  * psign1[28];
        p1[8]  = p1[8]  * psign1[32];
        p1[9]  = p1[9]  * psign1[36];
        p1[10] = p1[10] * psign1[40];
        p1[11] = p1[11] * psign1[44];
        p1[12] = p1[12] * psign1[48];
        p1[13] = p1[13] * psign1[52];
        p1[14] = p1[14] * psign1[56];
        p1[15] = p1[15] * psign1[60];
        p1 += 16;

        p2[0]  = p2[0]  * psign2[ 0];
        p2[1]  = p2[1]  * psign2[ 4];
        p2[2]  = p2[2]  * psign2[ 8];
        p2[3]  = p2[3]  * psign2[12];
        p2[4]  = p2[4]  * psign2[16];
        p2[5]  = p2[5]  * psign2[20];
        p2[6]  = p2[6]  * psign2[24];
        p2[7]  = p2[7]  * psign2[28];
        p2[8]  = p2[8]  * psign2[32];
        p2[9]  = p2[9]  * psign2[36];
        p2[10] = p2[10] * psign2[40];
        p2[11] = p2[11] * psign2[44];
        p2[12] = p2[12] * psign2[48];
        p2[13] = p2[13] * psign2[52];
        p2[14] = p2[14] * psign2[56];
        p2[15] = p2[15] * psign2[60];
        p2 += 16;

        p3[0]  = p3[0]  * psign3[ 0];
        p3[1]  = p3[1]  * psign3[ 4];
        p3[2]  = p3[2]  * psign3[ 8];
        p3[3]  = p3[3]  * psign3[12];
        p3[4]  = p3[4]  * psign3[16];
        p3[5]  = p3[5]  * psign3[20];
        p3[6]  = p3[6]  * psign3[24];
        p3[7]  = p3[7]  * psign3[28];
        p3[8]  = p3[8]  * psign3[32];
        p3[9]  = p3[9]  * psign3[36];
        p3[10] = p3[10] * psign3[40];
        p3[11] = p3[11] * psign3[44];
        p3[12] = p3[12] * psign3[48];
        p3[13] = p3[13] * psign3[52];
        p3[14] = p3[14] * psign3[56];
        p3[15] = p3[15] * psign3[60];
        p3 += 16;
    }
}

void E_ACELP_4tsearch(Float32 dn[], const Float32 cn[], const Float32 H[], float code[],
                      PulseConfig *config, Word16 ind[], Float32 y[])
{
    Float32 sign[L_SUBFR], vec[L_SUBFR];
    Float32 cor_x[16], cor_y[16], h_buf[4 * L_SUBFR];
    Float32 rrixix[4][16];
    Float32 rrixiy[4][256];
    Float32 dn2[L_SUBFR];
    Float32 psk, ps, alpk, alp = 0.0F;
    Word16 codvec[NB_PULSE_MAX];
    Word16 pos_max[4];
    Word16 dn2_pos[8 * 4];
    UWord8 ipos[NB_PULSE_MAX];
    Float32 *p0, *p1, *p2, *p3;
    Float32 *h, *h_inv;
    Word32 i, j, k, l, st, pos = 0, index, track;
    UWord8 iPulse;
    Float32 val;
    Float32 s;
    UWord8 restpulses;


    alp = config->alp;                    /* initial value for energy of all fixed pulses */
    /* set_i( (int *)codvec, 0, config->nb_pulse);                    */

    for (k=0; k<config->nb_pulse; k++) codvec[k] = 0;

    /*
     * Find sign for each pulse position.
     */

    acelp_pulsesign(cn, dn, dn2, sign, vec, alp);

    /*
     * Select the most important 8 position per track according to dn2[].
     */

    acelp_findcandidates(dn2, dn2_pos, pos_max, L_SUBFR, NB_TRACK_FCB_4T);

    /*
     * Compute h_inv[i].
     */

    acelp_hbuf( h_buf, &h, &h_inv, H );

    /*
     * Compute correlation matrices needed for the codebook search.
     */

    E_ACELP_corrmatrix(h, sign, vec, rrixix, rrixiy);


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
    psk = -1.0;
    alpk = 1.0;
    /*Number of iterations*/
    for (k = 0; k < config->nbiter; k++)
    {
        /* copy search order from hash-table */
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
            memset(vec, 0, L_SUBFR * sizeof(Float32));
        }
        else if (config->fixedpulses == 2) /* 2222 and 3322 */
        {
            /* first stage: fix 2 pulses */
            pos = 2;

            ind[0] = pos_max[ipos[0]];
            ind[1] = pos_max[ipos[1]];
            ps = dn[ind[0]] + dn[ind[1]];


            /*ind[1]>>2 and ind[0]>>2 and save*/
            /* ipos[1] and ipos[0]    and save*/
            alp = rrixix[ipos[0]][ind[0] >> 2] + rrixix[ipos[1]][ind[1] >> 2] +
                  rrixiy[ipos[0]][((ind[0] >> 2) << 4) + (ind[1] >> 2)];

            if (sign[ind[0]] < 0.0)
            {
                p0 = h_inv - ind[0];
            }
            else
            {
                p0 = h - ind[0];
            }
            if (sign[ind[1]] < 0.0)
            {
                p1 = h_inv - ind[1];
            }
            else
            {
                p1 = h - ind[1];
            }
            /*ptx = &vec   p1 and p0 already initialize*/
            vec[0] = p0[0] + p1[0];
            vec[1] = p0[1] + p1[1];
            vec[2] = p0[2] + p1[2];
            vec[3] = p0[3] + p1[3];
            for (i = 4; i < L_SUBFR; i += 6)
            {
                vec[i] = p0[i] + p1[i];
                vec[i + 1] = p0[i + 1] + p1[i + 1];
                vec[i + 2] = p0[i + 2] + p1[i + 2];
                vec[i + 3] = p0[i + 3] + p1[i + 3];
                vec[i + 4] = p0[i + 4] + p1[i + 4];
                vec[i + 5] = p0[i + 5] + p1[i + 5];
            }
        }
        else /* 3333 and above */
        {
            /* first stage: fix 4 pulses */
            pos = 4;

            ind[0] = pos_max[ipos[0]];
            ind[1] = pos_max[ipos[1]];
            ind[2] = pos_max[ipos[2]];
            ind[3] = pos_max[ipos[3]];
            ps = dn[ind[0]] + dn[ind[1]] + dn[ind[2]] + dn[ind[3]];

            p0 = h - ind[0];
            if (sign[ind[0]] < 0.0)
            {
                p0 = h_inv - ind[0];
            }

            p1 = h - ind[1];
            if (sign[ind[1]] < 0.0)
            {
                p1 = h_inv - ind[1];
            }

            p2 = h - ind[2];
            if (sign[ind[2]] < 0.0)
            {
                p2 = h_inv - ind[2];
            }

            p3 = h - ind[3];
            if (sign[ind[3]] < 0.0)
            {
                p3 = h_inv - ind[3];
            }
            /* pt =&vec; others already defined*/
            vec[0] = p0[0] + p1[0] + p2[0] + p3[0];
            for (i = 1; i < L_SUBFR; i += 3)
            {
                vec[i] = p0[i] + p1[i] + p2[i] + p3[i];
                vec[i + 1] = p0[i + 1] + p1[i + 1] + p2[i + 1] + p3[i + 1];
                vec[i + 2] = p0[i + 2] + p1[i + 2] + p2[i + 2] + p3[i + 2];
            }

            alp = 0.0F;
            alp += vec[0] * vec[0] + vec[1] * vec[1];
            alp += vec[2] * vec[2] + vec[3] * vec[3];

            for (i = 4; i < L_SUBFR; i += 6)
            {
                alp += vec[i] * vec[i];
                alp += vec[i + 1] * vec[i + 1];
                alp += vec[i + 2] * vec[i + 2];
                alp += vec[i + 3] * vec[i + 3];
                alp += vec[i + 4] * vec[i + 4];
                alp += vec[i + 5] * vec[i + 5];
            }

            alp *= 0.5F;
        }

        /* other stages of 2 pulses */
        for (j = pos, st = 0; j < config->nb_pulse; j += 2, st++)
        {
            if ((config->nb_pulse-j) >= 2)  /*pair-wise search*/
            {

                /*
                 * Calculate correlation of all possible positions
                 * of the next 2 pulses with previous fixed pulses.
                 * Each pulse can have 16 possible positions.
                 */
                acelp_h_vec_corr1(h, vec, ipos[j], sign, rrixix, cor_x, dn2_pos,
                                  config->nbpos[st]);
                acelp_h_vec_corr2(h, vec, ipos[j + 1], sign, rrixix, cor_y);

                /*
                 * Find best positions of 2 pulses.
                 */
                acelp_2pulse_search(config->nbpos[st], ipos[j], ipos[j + 1], &ps, &alp,
                                    &ind[j], &ind[j+1], dn, dn2_pos, cor_x, cor_y, rrixiy);

            }
            else  /*single pulse search*/
            {
                acelp_h_vec_corr2(h, vec, ipos[j], sign, rrixix, cor_x);
                acelp_h_vec_corr2(h, vec, ipos[j + 1], sign, rrixix, cor_y);
                E_ACELP_1pulse_search(ipos[j], ipos[j+1], &ps, &alp,
                                      &ind[j], dn, cor_x, cor_y);
            }
            if( j < (config->nb_pulse - 2) )
            {
                p0 = h - ind[j];
                if (sign[ind[j]] < 0.0)
                {
                    p0 = h_inv - ind[j];
                }

                p1 = h - ind[j + 1];
                if (sign[ind[j + 1]] < 0.0)
                {
                    p1 = h_inv - ind[j + 1];
                }

                vec[0] += p0[0] + p1[0];
                vec[1] += p0[1] + p1[1];
                vec[2] += p0[2] + p1[2];
                vec[3] += p0[3] + p1[3];
                for (i = 4; i < L_SUBFR; i += 6)
                {
                    vec[i] += p0[i] + p1[i];
                    vec[i + 1] += p0[i + 1] + p1[i + 1];
                    vec[i + 2] += p0[i + 2] + p1[i + 2];
                    vec[i + 3] += p0[i + 3] + p1[i + 3];
                    vec[i + 4] += p0[i + 4] + p1[i + 4];
                    vec[i + 5] += p0[i + 5] + p1[i + 5];

                }
            }
        }

        /* memorise the best codevector */

        ps = ps * ps;
        s = (alpk * ps) - (psk * alp);

        if (psk < 0)
        {
            s = 1.0F;
        }
        if (s > 0.0F)
        {
            psk = ps;
            alpk = alp;
            memcpy(codvec, ind, config->nb_pulse * sizeof(Word16));
        }
    }

    /*
     * Build the codeword, the filtered codeword and index of codevector, as well as store weighted correlations.
     */

    memset(code, 0,  L_SUBFR * sizeof(float));
    memset(y, 0,  L_SUBFR * sizeof(float));
    memset(ind, 0xffffffff, NPMAXPT * 4 * sizeof(Word16));
    for (k = 0; k < config->nb_pulse; k++)
    {
        i = codvec[k];  /* read pulse position  */
        val = sign[i];  /* read sign            */

        index = i / 4;  /* pos of pulse (0..15) */
        track = i % 4;
        if (val > 0)
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

        p0 = h_inv - codvec[k];
        for(i=0; i<L_SUBFR; i++)
        {
            y[i] += *p0++;
        }
    }
}


/*
 * E_ACELP_4t
 *
 * Parameters:
 *    dn          I: corr. between target and h[].
 *    cn          I: residual after Word32 term prediction
 *    H           I: impulse response of weighted synthesis filter (Q12)
 *    code        O: algebraic (fixed) codebook excitation (Q9)
 *    y           O: filtered fixed codebook excitation (Q9)
 *    nbbits      I: 20, 36, 44, 52, 64, 72 or 88 bits
 *    mode        I: speech mode
 *    _index      O: index
 *
 * Function:
 *    20, 36, 44, 52, 64, 72, 88 bits algebraic codebook.
 *    4 tracks x 16 positions per track = 64 samples.
 *
 *    20 bits 5 + 5 + 5 + 5 --> 4 pulses in a frame of 64 samples.
 *    36 bits 9 + 9 + 9 + 9 --> 8 pulses in a frame of 64 samples.
 *    44 bits 13 + 9 + 13 + 9 --> 10 pulses in a frame of 64 samples.
 *    52 bits 13 + 13 + 13 + 13 --> 12 pulses in a frame of 64 samples.
 *    64 bits 2 + 2 + 2 + 2 + 14 + 14 + 14 + 14 -->
 *                                  16 pulses in a frame of 64 samples.
 *    72 bits 10 + 2 + 10 + 2 + 10 + 14 + 10 + 14 -->
 *                                  18 pulses in a frame of 64 samples.
 *    88 bits 11 + 11 + 11 + 11 + 11 + 11 + 11 + 11 -->
 *                                  24 pulses in a frame of 64 samples.
 *
 *    All pulses can have two (2) possible amplitudes: +1 or -1.
 *    Each pulse can sixteen (16) possible positions.
 *
 * Returns:
 *    void
 */
static void E_ACELP_4t(
    Float32 dn[],
    Float32 cn[],
    Float32 H[],
    Float32 R[],
    int acelpautoc,
    float code[],
    int cdk_index,
    int _index[]
    ,const short L_frame,
    const short last_L_frame,
    const long total_brate,
    const short i_subfr
)
{
    PulseConfig config;
    Word16 ind[NPMAXPT*4];
    Float32 y[L_SUBFR];

    config = PulseConfTable[cdk_index];

    if( L_frame != last_L_frame && total_brate == ACELP_24k40 && i_subfr < 5*L_SUBFR )
    {
        (config.nbiter)--;
        config.nbiter = max(config.nbiter,1);
    }

    if (acelpautoc & 0x01)
    {
        E_ACELP_4tsearchx(dn, cn, R, code, &config, ind);
    }
    else
    {
        E_ACELP_4tsearch(dn, cn, H, code, &config, ind, y);
    }

    E_ACELP_indexing(code, config, NB_TRACK_FCB_4T, _index);

    return;
}

short E_ACELP_indexing(
    Float32 code[],
    PulseConfig config,
    int num_tracks,
    int prm[]
)
{
    unsigned short track;
    int p[NB_TRACK_FCB_4T], wordcnt;
    int k;
    unsigned short idxs[MAX_IDX_LEN], maxppos;
    unsigned long s[NB_TRACK_FCB_4T], n[NB_TRACK_FCB_4T];
    int maxp;
    short saved_bits;

    assert(num_tracks == NB_TRACK_FCB_4T);

    saved_bits = 0;

    /*
     * Code state of pulses of all tracks
     * */
    wordcnt = (config.bits + 15) >> 4;     /* ceil(bits/16) */
    for (k=0; k<wordcnt; k++)
    {
        idxs[k] = 0;
    }
    if (config.bits == 43)    /* EVS pulse indexing */
    {
        saved_bits = E_ACELP_code43bit(code,s, p, idxs);
    }
    else
    {

        for (track=0; track < num_tracks; track++)
        {
            /* Code track of length 2^4 where step between tracks is 4. */
            E_ACELP_codearithp(code+track,n+track,s+track, p+track, num_tracks, 16);
        }
        fcb_pulse_track_joint( idxs, wordcnt, s, p, num_tracks );

    }

    /* check if we need to code track positions */
    switch(config.codetrackpos)
    {
    case TRACKPOS_FIXED_TWO:
        /* Code position of consecutive tracks with single extra pulses */

        /* Find track with one pulse less. */
        if (p[0] == p[1])
        {
            /* Either 1100 or 0011 */
            if (p[1]>p[2])
            {
                track = 0;    /* 1100 */
            }
            else
            {
                track = 2;      /* 0011 */
            }
        }
        else
        {
            /* Either 0110 or 1001 */
            if (p[0] < p[1])
            {
                track = 1;    /* 0110 */
            }
            else
            {
                track = 3;    /* 1001 */
            }
        }
        /* Multiply by number of possible states (=shift by two) and
         * add actual state. */
        longshiftleft(idxs,2,idxs,wordcnt);
        longadd(idxs, &track, wordcnt, 1);
        break;
    case TRACKPOS_FREE_THREE:
        /* Code position of track with one pulse less than others */

        /* Find track with one pulse less. */
        maxp = p[0];
        maxppos = 0;
        for (track=1; track < 4; track++)
        {
            if (p[track] < maxp)
            {
                maxppos = track;
                break;
            }
        }
        /* Multiply by number of possible states (=shift by two) and
         * add actual state. */
        longshiftleft(idxs,2,idxs,wordcnt);
        longadd(idxs, &maxppos, wordcnt, 1);
        break;
    case TRACKPOS_FREE_ONE:
        /* Code position of track with one pulse more than others */

        /* Find track with one pulse more. */
        maxp = p[0];
        maxppos = 0;
        for (track=1; track < 4; track++)
        {
            if (p[track] > maxp)
            {
                maxppos = track;
                break;
            }
        }
        /* Multiply by number of possible states (=shift by two) and
         * add actual state. */
        longshiftleft(idxs, 2, idxs, wordcnt);
        longadd(idxs, &maxppos, wordcnt, 1);
        break;
    case TRACKPOS_FIXED_EVEN:
    case TRACKPOS_FIXED_FIRST:
        break;
    default:
        printf("Codebook mode not implemented.");
        assert(0);   /* mode not yet implemented*/
        break;
    }

    /* cast to output buffer */
    for (k=0; k<wordcnt; k++)
    {
        prm[k] = idxs[k];
    }


    return (saved_bits);
}


/*--------------------------------------------------------------------------*
 * E_ACELP_innovative_codebook
 *
 * Find innovative codebook.
 *--------------------------------------------------------------------------*/
void E_ACELP_innovative_codebook(
    float *exc,        /* i  : pointer to the excitation frame                  */
    int T0,          /* i  : integer pitch lag                                */
    int T0_frac,     /* i  : fraction of lag                                  */
    int T0_res,      /* i  : pitch resolution                                 */
    float gain_pit,    /* i  : adaptive codebook gain                           */
    float tilt_code,   /* i  : tilt factor                                      */
    int mode,        /* i  : innovative codebook mode                         */
    int pre_emphasis,/* i  : use pre_emphasis                                 */
    int pitch_sharpening, /* i  : use pitch sharpening                        */
    int phase_scrambling, /* i  : use phase scrambling                        */
    int formant_enh, /* i  : use formant enhancement                          */
    int formant_tilt,/* i  : use tilt of formant enhancement                  */
    float formant_enh_num, /* i  : formant enhancement numerator weighting factor*/
    float formant_enh_den, /* i  : formant enhancement denominator weighting factor*/
    const short i_subfr,     /* i  : subframe index                                   */
    const float *Aq,         /* i  : quantized LPC coefficients                       */
    float *h1,         /* i  : impulse response of weighted synthesis filter    */
    float *xn,         /* i  : Close-loop Pitch search target vector            */
    float *cn,         /* i  : Innovative codebook search target vector         */
    float *y1,         /* i  : zero-memory filtered adaptive excitation         */
    float *y2,         /* o  : zero-memory filtered algebraic excitation        */
    int acelpautoc,  /* i  : autocorrelation mode enabled                     */
    int **pt_indice, /* i/o: quantization indices pointer                     */
    float *code        /* o  : innovative codebook                              */
    ,const short L_frame,     /* i  : length of the frame                              */
    const short last_L_frame,/* i  : length of the last frame                         */
    const long total_brate   /* i  : total bit-rate                                   */
)
{
    float xn2[L_SUBFR], cn2[L_SUBFR], dn[L_SUBFR], h2[L_SUBFR];
    float Rw2[L_SUBFR];
    int i,k;
    float pitch;

    pitch = (float)T0+(float)T0_frac/(float)T0_res;

    /* Update target vector for ACELP codebook search */
    updt_tar(xn, xn2, y1, gain_pit, L_SUBFR);

    /* Include fixed-gain pitch contribution into impulse resp. h1[] */
    mvr2r( h1, h2, L_SUBFR);
    cb_shape( pre_emphasis, pitch_sharpening, phase_scrambling, formant_enh, formant_tilt,
              formant_enh_num, formant_enh_den, Aq, h2, tilt_code, pitch );

    /* Correlation between target xn2[] and impulse response h1[] */
    if (acelpautoc & 0x01)
    {
        corr_xh(h2, Rw2, h2, L_SUBFR);
        for (k=0; k<L_SUBFR; k++)
        {
            cn2[k] = xn2[k];
            for (i=0; i<k; i++)
            {
                cn2[k]-=cn2[i]*h2[k-i];
            }
        }

        E_ACELP_toeplitz_mul( Rw2, cn2, dn );
    }
    else
    {
        updt_tar(cn, cn2, &exc[i_subfr], gain_pit, L_SUBFR);
        corr_xh(xn2, dn, h2, L_SUBFR);
    }

    /* Innovative codebook search */
    if (mode < ACELP_FIXED_CDK_NB)
    {
        E_ACELP_4t( dn, cn2, h2, Rw2, acelpautoc, code, mode, *pt_indice, L_frame, last_L_frame, total_brate, i_subfr );
    }
    else
    {
        assert(0);
    }
    *pt_indice += 8;


    /* Generate weighted code */
    set_f( y2, 0.0f, L_SUBFR);
    for (i=0; i<L_SUBFR; i++)
    {
        /* Code is sparse, so check which samples are non-zero */
        if (code[i] != 0)
        {
            for (k=0; k<L_SUBFR-i; k++)
            {
                y2[i+k] += code[i]*h2[k];
            }
        }
    }

    /*-------------------------------------------------------*
     * - Add the fixed-gain pitch contribution to code[].    *
     *-------------------------------------------------------*/
    cb_shape( pre_emphasis, pitch_sharpening, phase_scrambling, formant_enh, formant_tilt,
              formant_enh_num, formant_enh_den, Aq, code, tilt_code, pitch );

    return;

}



/*--------------------------------------------------------------------------*
 * E_ACELP_codearithp
 *
 * Fixed bit-length arithmetic coding of pulses
 * v - (input) pulse vector
 * s - (output) encoded state
 * n - (output) range of possible states (0...n-1)
 * p - (output) number of pulses found
 * len - (input) length of pulse vector
 * trackstep - (input) step between tracks
 *--------------------------------------------------------------------------*/
static void E_ACELP_codearithp(const float v[], long unsigned *n, long unsigned *ps, int *p, int trackstep, int tracklen)
{
    int k, h, t, pos[9], sig[9], posno, tmp, L_subfr;
    long unsigned s;
    posno = 0;
    L_subfr = trackstep*tracklen;
    for (k=t=0; k < L_subfr; k+=trackstep, t++)
    {
        tmp = (v[k]>0?1:-1);        /* sign */
        for (h=0; h < v[k]*tmp; h++)
        {
            pos[posno] = t;
            sig[posno] = tmp;
            posno++;
            if (posno > 9)
            {
                assert(0);
            }
        }
    }
    *p = posno;

    s = 0;
    for (k=0; k < posno; k++)
    {
        /* check if next position is the same as this one */
        if ((k==posno-1) || (pos[k] != pos[k+1]))
        {
            /* next position is not the same (or we are at the last position)
             * -> save sign */
            s <<= 1;
            if (sig[k] < 0)
            {
                s++;
            }
        }
        s += pulsestostates[pos[k]][k];
    }
    *ps = s;
    if (posno)
    {
        *n = pulsestostates[tracklen][posno-1];
    }
    else
    {
        *n = 0;
    }

}

void fcb_pulse_track_joint( unsigned short *idxs, int wordcnt, long unsigned *index_n, int *pulse_num, int track_num )
{
    int hi_to_low[10] = { 0, 0, 0, 3, 9, 5, 3, 1, 8, 8};

    unsigned long long index;
    unsigned int index_mask;
    int indx_tmp,indx_flag,indx_flag_1;
    int track,track_num1,pulse_num0,pulse_num1;
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

    if (indx_flag_2 >= 1)
    {
        hi_to_low[7] = 9;
        index_mask = 0xffffff;
    }
    else
    {
        hi_to_low[7] = 1;
        if (indx_flag>=track_num)
        {
            hi_to_low[4]=9;
            index_mask = 0xffff;
        }
        else
        {
            hi_to_low[4]=1;
            index_mask = 0xff;
        }
    }

    if (indx_flag_1>=track_num)
    {
        indx_tmp = 0;
        index = index_n[0] >> low_len[pulse_num[0]];
        for (track=1; track < track_num; track++)
        {
            pulse_num0 = pulse_num[track-1];
            pulse_num1 = pulse_num[track];
            indx_tmp = index_n[track] >> low_len[pulse_num1];
            index = index * indx_fact[pulse_num1] + indx_tmp;

            index_n[track-1] = ( index_n[track-1] & low_mask[pulse_num0] ) + ( ( index << low_len[pulse_num0] )& index_mask ) ;
            index = index >> hi_to_low[pulse_num0];
        }
        track_num1 = track_num - 1;
        pulse_num1 = pulse_num[track_num1];
        index_n[track_num1] = ( ( index_n[track_num1] & low_mask[pulse_num1] ) + ( index << low_len[pulse_num1] ) ) & index_mask;
        index = index >> hi_to_low[pulse_num1];
        if (indx_flag>=track_num)
        {
            if (indx_flag_2 >= 1)
            {
                idxs[0] = index_n[0] & 0xffff;
                idxs[1] = ( ( index_n[1] << 8 ) + ( index_n[0] >> 16) ) & 0xffff;
                idxs[2] = index_n[1] >> 8;
                idxs[3] = index_n[2] & 0xffff;
                idxs[4] = ( ( index_n[3] << 8 ) + ( index_n[2] >> 16) ) & 0xffff;
                idxs[5] = index_n[3] >> 8;
                for (track=6; track < wordcnt; track++)
                {
                    idxs[track] = index&0xffff;
                    index = index >> 16;
                }
            }
            else
            {
                for (track=0; track < track_num; track++)
                {
                    idxs[track] = index_n[track];
                }
                for (track=track_num; track < wordcnt; track++)
                {
                    idxs[track] = index&0xffff;
                    index = index >> 16;
                }
            }
        }
        else
        {
            idxs[0] = ( index_n[0] << 8 ) + index_n[1];
            idxs[1] = ( index_n[2] << 8 ) + index_n[3];
            for (track=2; track < wordcnt; track++)
            {
                idxs[track] = index&0xffff;
                index = index >> 16;
            }
        }
    }
    else
    {
        index = index_n[0];
        for (track=1; track < 4; track++)
        {
            pulse_num1 = pulse_num[track];
            index = ( index << index_len[pulse_num1] ) + index_n[track];
        }
        for (track=0; track < wordcnt; track++)
        {
            idxs[track] = index&0xffff;
            index = index >> 16;
        }
    }

    return;
}
