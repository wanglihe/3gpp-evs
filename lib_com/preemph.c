/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "prot.h"


/*-------------------------------------------------------------*
 * preemph()
 *
 * Preemphasis: filtering through 1 - mu z^-1
 *-------------------------------------------------------------*/

void preemph(
    float *signal,  /* i/o: signal             */
    const float mu,       /* i  : preemphasis factor */
    const short L,        /* i  : vector size        */
    float *mem      /* i/o: memory (x[-1])     */
)
{
    short i;
    float temp;

    temp = signal[L-1];
    for (i=L-1; i>0; i--)
    {
        signal[i] = signal[i] - mu*signal[i-1];
    }

    signal[0] = signal[0] - mu*(*mem);
    *mem = temp;

    return;
}
