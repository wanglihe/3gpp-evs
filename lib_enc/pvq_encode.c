/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "prot.h"
#include "rom_com.h"
#include "prot.h"
#include "math.h"

static void pyramidSearch(const float  *s, short L, short  Ptot, float A,  short *ztak, float *stak ) ;

void pvq_encode(
    Encoder_State *st,
    const float *x,         /* i:   vector to quantize               */
    short *y,               /* o:   quantized vector (non-scaled int)*/
    float *xq,              /* o:   quantized vector (scaled float)  */
    const short pulses,     /* i:   number of allocated pulses       */
    const short dim,        /* i:   Length of vector                 */
    const float gain        /* i:   Gain                             */
)
{
    PvqEntry entry;

    pyramidSearch(x, dim, pulses, gain, y, xq );

    entry = mpvq_encode_vec(y, dim, pulses);
    if(dim != 1)
    {
        rc_enc_bits(st, (unsigned int)entry.lead_sign_ind, 1);
        rc_enc_uniform(st, entry.index, entry.size);
    }
    else
    {
        rc_enc_bits(st, (unsigned int)entry.lead_sign_ind, 1);
    }
    return;

}

static float L1norm(
    const float *s,
    float *abso_s,
    short *sign_s,
    short L
)
{
    int i;
    float ftmp, r = 0.0f;

    for( i=0; i < L; i++ )
    {
        ftmp      = s[i];
        sign_s[i] = (short)sign(ftmp);
        abso_s[i] = fabs(ftmp);
        r +=  abso_s[i];
    }
    return r;
}


void pyramidSearch(const float  *s,
                   short L,
                   short  Ptot,
                   float A,
                   short *ztak,
                   float *stak
                  )
{
    short i, high_pulse_dens;
    short z[PVQ_MAX_BAND_SIZE];
    float sL1;
    float energy,xcorr;
    short P;
    float energy1,energy2;
    float xcorr1sq, xcorr2sq ;
    float zscale, sscale;
    float abso_s[PVQ_MAX_BAND_SIZE];
    short sign_s[PVQ_MAX_BAND_SIZE];
    short n = 0;

    high_pulse_dens = Ptot > (short )(0.501892089843750f * L);
    sL1 = L1norm(s, abso_s,sign_s, L );
    if (high_pulse_dens  &&  sL1 > 0 && A>0)
    {
        P      =     Ptot - PYR_OFFSET;
        zscale =  P/sL1;
        for(i=0; i < L; i++)
        {
            z[i] = floor(zscale * abso_s[i]) ;
        }
    }
    else
    {
        for( i=0; i < L; i++)
        {
            z[i] = 0;
        }
    }

    energy = 0.0f;
    if (sL1 > 0 && A > 0)
    {
        xcorr = 0.0f;
        P = 0;
        for( i = 0; i < L; i++)
        {
            xcorr  += abso_s[i] * z[i];
            energy += z[i] * z[i];
            P += z[i];
        }
        energy = 0.5f*energy;
        while( P < Ptot)
        {
            energy  += 0.5f;
            n = 0;
            xcorr1sq   = xcorr + abso_s[0];
            xcorr1sq  *= xcorr1sq;
            energy1  = energy + z[0];
            for( i = 1; i < L; i++)
            {
                xcorr2sq   = xcorr + abso_s[i] ;
                xcorr2sq   *= xcorr2sq;
                energy2  = energy + z[i] ;

                if ((xcorr1sq*energy2) < (xcorr2sq*energy1))
                {
                    n = i;
                    xcorr1sq = xcorr2sq;
                    energy1  = energy2;
                }
            }
            z[n]++;
            xcorr  += abso_s[n];
            energy = energy1;
            P++;
        }
    }
    else
    {
        if (L > 1)
        {
            z[0] = (short)(Ptot*0.5f);
            z[L-1] = - (Ptot - z[0]);
            energy  = 0.5f*(z[0] * z[0] + z[L-1] * z[L-1]);
        }
        else
        {
            z[0] = Ptot;
            energy = 0.5f*(z[0] * z[0]);
        }
    }

    sscale = (float) A/(float)sqrt(2.0f * energy);
    for(i = 0; i < L; i++)
    {
        ztak[i] = z[i] * sign_s[i];
        stak[i] = sscale * ztak[i] ;
    }
    return;
}
