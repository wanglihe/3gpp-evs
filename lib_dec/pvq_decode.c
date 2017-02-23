/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "prot.h"
#include "rom_com.h"
#include "math.h"


/*-------------------------------------------------------------------*
 * Function pvq_decode()                                             *
 *                                                                   *
 * PVQ subvector decoding algorithm                                  *
 *-------------------------------------------------------------------*/

void pvq_decode(
    Decoder_State *st,
    float *xq,        /* o:   decoded vector (scaled float)    */
    short *y,         /* o:   decoded vector (non-scaled short)*/
    const short k_val,      /* i:   number of allocated pulses       */
    const short dim,        /* i:   Length of vector                 */
    const float gain        /* i:   Gain                             */
)
{
    short i;
    float gain_fac;
    float yy;
    short output[PVQ_MAX_BAND_SIZE];       /* short interface as in STL-FIP */
    unsigned int h_mem[1+KMAX_NON_DIRECT+1];              /* allocate max offset memory  for dim 6 */
    PvqEntry entry;

    entry = get_size_mpvq_calc_offset(dim, k_val, h_mem); /* get size & prepare H(adaptive table for entry.size=N_MPVQ(dim,k_val) */

    if( dim != 1)
    {
        entry.lead_sign_ind  = (short)rc_dec_bits(st, 1);
        entry.index = rc_dec_uniform(st, entry.size);    /* NB so far no PVQ-size wc  is exactly 2^32-1 */
        /* safety check in case of bit errors */
        if( entry.index >= entry.size || st->ber_occured_in_pvq != 0 )
        {
            st->ber_occured_in_pvq = 1;
            st->BER_detect = 1;
            entry.index = 0;   /* a zero index will essentially disable PVQ index decompostion complexity */
        }

    }
    else
    {
        entry.lead_sign_ind = (short)rc_dec_bits(st, 1);    /* always a single sign bit */
        entry.index         = 0;
    }
    mpvq_decode_vec(&entry, h_mem, output);

    for(i=0; i<dim ; i++)
    {
        y[i]= (int)output[i];
    }
    /* Find decoded vector energy */
    yy = 0;
    for( i = 0; i < dim; i++ )
    {
        yy += y[i]*y[i];
    }

    /* Apply scaling , always energy in yy */
    gain_fac = gain * 1.0f/(float)sqrt(yy);
    for( i = 0; i < dim; i++)
    {
        xq[i] = y[i] * gain_fac;
    }

    return;
}
