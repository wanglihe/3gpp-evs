/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "options.h"
#include "prot.h"
#include "cnst.h"
#include "rom_com.h"
#include "rom_enc.h"


/*---------------------------------------------------------------------------

   function name: SetModeIndex
   description:   function for configuring the codec between frames - currently bit rate and bandwidth only
                  must not be called while a frame is encoded - hence mutexes must be used in MT environments

   format:        BANDWIDTH*16 + BITRATE (mode index)

  ---------------------------------------------------------------------------*/

void SetModeIndex(
    Encoder_State *st,
    const long    total_brate,
    const short   bwidth)
{

    /* Reconfigure the core coder */
    if( (st->last_total_brate != total_brate) || (st->last_bwidth != bwidth) || (st->last_codec_mode == MODE1) || (st->rf_mode_last != st->rf_mode) )
    {
        core_coder_mode_switch( st, st->bwidth, total_brate );
    }


    return;
}
