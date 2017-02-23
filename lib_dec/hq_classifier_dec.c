/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"

/*--------------------------------------------------------------------------*
 * hq_classifier_dec()
 *
 * HQ mode selector (decision_matrix)
 *--------------------------------------------------------------------------*/

short hq_classifier_dec(            /* o  : Consumed bits                   */
    Decoder_State *st,                /* i/o: decoder state structure         */
    const long  core_brate,         /* i  : Core bit rate                   */
    const short length,             /* i  : Frame length                    */
    short *is_transient,      /* o  : Transient flag                  */
    short *hqswb_clas         /* o  : HQ class                        */
)
{
    short bits;

    if( length >= L_FRAME32k && core_brate <= HQ_32k )
    {
        *hqswb_clas = (short)get_next_indice( st, 2 );
        if ( length == L_FRAME48k )
        {
            if ( *hqswb_clas == 0 )
            {
                *hqswb_clas = HQ_GEN_FB;
            }
        }
    }
    else
    {
        *hqswb_clas = (short)get_next_indice( st, 1 );
    }

    *is_transient = 0;
    if( *hqswb_clas == HQ_TRANSIENT )
    {
        *is_transient = 1;
    }

    if ( *hqswb_clas == HQ_NORMAL && length == L_FRAME32k && core_brate <= HQ_32k)
    {
        *hqswb_clas = HQ_GEN_SWB;
    }

    if( length >= L_FRAME32k && core_brate <= HQ_32k )
    {
        bits = 2;
    }
    else
    {
        bits = 1;
    }

    return bits;
}
