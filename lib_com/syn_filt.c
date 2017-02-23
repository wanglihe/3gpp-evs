/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"

/*------------------------------------------------------------------*
 * syn_filt()
 *
 * perform the synthesis filtering 1/A(z)
 *------------------------------------------------------------------*/

void syn_filt(
    const float a[],      /* i  : LP filter coefficients                     */
    const short m,        /* i  : order of LP filter                         */
    const float x[],      /* i  : input signal                               */
    float y[],      /* o  : output signal                              */
    const short l,        /* i  : size of filtering                          */
    float mem[],    /* i/o: initial filter states                      */
    const short update_m  /* i  : update memory flag: 0 --> no memory update */
)                         /*                          1 --> update of memory */
{
    short i, j;
#if !defined(TCXLTP_LTP_ORDER)
    float buf[L_FRAME48k + L_FRAME48k/2 + M];    /* temporary synthesis buffer */
#else
    float buf[L_FRAME48k + L_FRAME48k/2 + TCXLTP_LTP_ORDER];    /* temporary synthesis buffer */
#endif
    float s, *yy;

    yy = &buf[0];

    /*------------------------------------------------------------------*
     * copy initial filter states into synthesis buffer and do synthesis
     *------------------------------------------------------------------*/

    for (i = 0; i < m; i++)
    {
        *yy++ = mem[i];
    }

    /*-----------------------------------------------------------------------*
     * Do the filtering
     *-----------------------------------------------------------------------*/

    for (i = 0; i < l; i++)
    {
        s = x[i];
        for (j = 1; j <= m; j++)
        {
            s -= a[j]*yy[i-j];
        }

        yy[i] = s;
        y[i] = s;
    }

    /*------------------------------------------------------------------*
     * Update memory if required
     *------------------------------------------------------------------*/

    if (update_m)
    {
        for (i = 0; i < m; i++)
        {
            mem[i] = yy[l-m+i];
        }
    }

    return;
}


/*-------------------------------------------------------------------*
 * synth_mem_updt2()
 *
 * Update of synthesis filter memories in case of ACELP@12k8 <-> ACELP@16k switching
 *--------------------------------------------------------------------*/

void synth_mem_updt2(
    const short L_frame,        /* i  : frame length                            */
    const short last_L_frame,   /* i  : frame length                            */
    float old_exc[],      /* i/o: excitation buffer                       */
    float mem_syn_r[],    /* i/o: synthesis filter memory                 */
    float mem_syn2[],     /* o  : synthesis filter memory for find_target */
    float mem_syn[],      /* o  : synthesis filter memory for find_target */
    const short dec             /* i  : flag for decoder indication             */
)
{
    short mem_syn_r_size_old, mem_syn_r_size_new;

    /* Residual and update old_exc */
    if( dec == DEC )
    {
        lerp( old_exc+L_EXC_MEM_DEC-(last_L_frame+last_L_frame/2), old_exc+L_EXC_MEM_DEC-(L_frame+L_frame/2), L_frame+L_frame/2, last_L_frame+last_L_frame/2 );
    }
    else
    {
        lerp( old_exc+L_EXC_MEM-last_L_frame, old_exc+L_EXC_MEM-L_frame, L_frame, last_L_frame );
    }

    /*Resamp memory*/
    /*Size of LPC syn memory*/
    mem_syn_r_size_old = (short)(1.25*last_L_frame/20.f);
    mem_syn_r_size_new = (short)(1.25*L_frame/20.f);

    lerp( mem_syn_r + L_SYN_MEM - mem_syn_r_size_old, mem_syn_r + L_SYN_MEM - mem_syn_r_size_new, mem_syn_r_size_new, mem_syn_r_size_old );

    mvr2r( mem_syn_r+L_SYN_MEM-M, mem_syn2, M );

    if( mem_syn != NULL )
    {
        mvr2r( mem_syn2, mem_syn, M );
    }

    return;
}

