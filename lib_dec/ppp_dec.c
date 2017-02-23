/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <stdlib.h>
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"

/*-------------------------------------------------------------------
 * ppp_quarter_decoder()
 *
 * PPP quarter decoder
 *-------------------------------------------------------------------*/

void ppp_quarter_decoder(
    Decoder_State *st,               /* i/o: decoder state structure   */
    DTFS_STRUCTURE *CURRCW_Q_DTFS,    /* i/o: Current CW DTFS */
    int   prevCW_lag,        /* i  : Previous lag */
    float *lastLgainD,       /* i/o: Last gain lowband */
    float *lastHgainD,       /* i/o: Last gain highwband */
    float *lasterbD,         /* i/o: Last ERB vector */
    short bfi,					/* i  : FER flag */
    DTFS_STRUCTURE PREV_CW_D          /* i  : Previous DTFS */
)
{
    DTFS_STRUCTURE *PREVDTFS = DTFS_new();

    float tmp, temp_pl = (float) prevCW_lag, temp_l = (float) CURRCW_Q_DTFS->lag;
    int l = CURRCW_Q_DTFS->lag;
    int POWER_IDX,AMP_IDX[2];
    float Erot = 0.0, z = 0.0;
    short num_erb = 24;

    if ( CURRCW_Q_DTFS->upper_cut_off_freq == 4000.0 )
    {
        num_erb = 22;
    }
    else if ( CURRCW_Q_DTFS->upper_cut_off_freq == 6400.0 )
    {
        num_erb = 24;
    }

    DTFS_copy(PREVDTFS, PREV_CW_D);
    if( bfi == 0 )
    {
        POWER_IDX = get_next_indice( st, 6 );
        AMP_IDX[0] = get_next_indice( st, 6 );
        AMP_IDX[1] = get_next_indice( st, 6 );

        /* Amplitude Dequantization */
        DTFS_dequant_cw(prevCW_lag,POWER_IDX,AMP_IDX, lastLgainD, lastHgainD, lasterbD,CURRCW_Q_DTFS,num_erb);
    }

    /* Copying phase spectrum over */
    DTFS_adjustLag(PREVDTFS,l);

    z=((L_FRAME-temp_l)*(temp_l+temp_pl))/(2*temp_l*temp_pl);

    Erot=(float) (temp_l - rint_new(temp_l*(z - floor(z))));

    DTFS_phaseShift(PREVDTFS,(float)(PI2*Erot/CURRCW_Q_DTFS->lag)) ;
    DTFS_car2pol(PREVDTFS);

    mvr2r(PREVDTFS->b, CURRCW_Q_DTFS->b, (short)(CURRCW_Q_DTFS->lag>>1)+1);

    DTFS_pol2car(CURRCW_Q_DTFS);

    tmp = (float) get_next_indice( st, 3 );
    DTFS_phaseShift(CURRCW_Q_DTFS,(float)(PI2*(tmp-3)/CURRCW_Q_DTFS->lag)) ;

    free( PREVDTFS );

    return;
}
