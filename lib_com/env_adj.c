/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"


/*--------------------------------------------------------------------------*
 * env_adj()
 *
 * Adjust the band energies of noise-fill and low resolution bands
 *--------------------------------------------------------------------------*/

void env_adj (
    const short *pulses,            /* i  : number of pulses per band           */
    const short length,             /* i  : length of spectrum                  */
    const short last_sfm,           /* i  : index of the last band              */
    float *adj,               /* o  : adjustment factors for the envelope */
    const float env_stab,           /* i  : Envelope stability parameter        */
    const short *sfmsize            /* i  : Band widths                         */
)
{
    short i, j, group;
    int npul;
    short att_state;
    short start, len;
    float tmp;
    float gain_adj;
    short idx;

    att_state = 0;
    len  = 0;
    start = 0;

    /* Find attenuation levels */
    for( i = 0; i <= last_sfm ; i++ )
    {
        group = (sfmsize[i] >> 3) - 1;
        npul  = pulses[i];

        if( length == L_FRAME32k )
        {
            if( npul == 0 )
            {
                /* Noise filled band */
                if ( group <= 1 )
                {
                    if ( i > 0 && pulses[i-1] != 0 && pulses[i+1] != 0 )
                    {
                        adj[i] = 0.36f;
                    }
                    else if ( i > 0 && ( pulses[i-1] == 0 || pulses[i+1] == 0) )
                    {
                        adj[i] = 0.54f;
                    }
                    else
                    {
                        adj[i] = 0.72f;
                    }
                }
                else if (i < last_sfm)
                {
                    if ( pulses[i-1] != 0 && pulses[i+1] != 0 )
                    {
                        adj[i] = 0.54f;
                    }
                    else
                    {
                        adj[i] = 0.72f;
                    }
                }
                else
                {
                    adj[i] = 0.72f;
                }

                if( att_state == 0 )
                {
                    start = i;
                }

                len++;
                att_state = 1;
            }
            else
            {
                adj[i] = 1.0f;
                if(att_state == 1) /* End of attenuation region found */
                {
                    tmp = min(1, max(0, len-ENV_ADJ_START)*(1.0f/ENV_ADJ_INCL));
                    for( j = start; j < i ; j++ )
                    {
                        adj[j] = max(tmp + (1-tmp)*adj[j],env_stab);
                    }
                    len = 0;
                    att_state = 0;
                }
            }
        }
        /* length == L_FRAME16k */
        else
        {
            /* Calculate low accuracy band attenuation */
            gain_adj = 1.0f;
            if( npul > 0 && npul < MAX_P_ATT )
            {
                idx = (short)(npul * att_step[group] + 0.5f) - 1;
                if( idx < MAX_P_ATT )
                {
                    gain_adj = gain_att[idx];
                }
            }
            adj[i] = gain_adj;
        }
    }

    /* Check if the sequence ended with an attenuation region */
    if( att_state == 1 )
    {
        tmp = min(1, max(0, len-ENV_ADJ_START)*(1.0f/ENV_ADJ_INCL));

        for( j = start; j < i ; j++ )
        {
            adj[j] = max(tmp + (1-tmp)*adj[j],env_stab);
        }
    }

    return;
}
