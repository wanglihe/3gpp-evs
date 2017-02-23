/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "prot.h"
#include "rom_com.h"


/*-------------------------------------------------------------------*
 * Local function
 *-------------------------------------------------------------------*/

static void reorder_lsfs( float *lsf, float min_dist, const short n, int sr_core );


/*-------------------------------------------------------------------*
 * getLookAheadResSig()
 *
 *
 *-------------------------------------------------------------------*/

void getLookAheadResSig(
    float *speechLookAhead,
    const float *A,
    float *res,
    int L_frame,
    int L_subfr,
    int m,
    int numSubFrame
)
{
    const float *p_A;
    short i_subfr;
    short subfr_len[2] = { 0 };

    subfr_len[0] = L_subfr;
    subfr_len[1] = L_frame < L_FRAME16k ? (short)(0.75*L_subfr) : L_subfr;

    p_A = A;
    for (i_subfr=0; i_subfr<numSubFrame*L_subfr; i_subfr+=L_subfr)
    {
        residu( p_A, m, &speechLookAhead[i_subfr], &res[i_subfr], subfr_len[i_subfr/L_subfr] );
        p_A += (m+1);
    }

    return;
}

/*-------------------------------------------------------------------*
 * updatelsfForConcealment()
 *
 *
 *-------------------------------------------------------------------*/

void updatelsfForConcealment(
    HANDLE_PLC_ENC_EVS decState,
    float *lsf
)
{
    short i;

    for (i=0; i<M; i++)
    {
        decState->lsf_adaptive_mean[i] = ( decState->lsfoldbfi1[i]+ decState->lsfoldbfi0[i]+lsf[i])/3;
        decState->lsfoldbfi1[i] = decState->lsfoldbfi0[i];
        decState->lsfoldbfi0[i] = lsf[i];
    }

    return;
}


/*-------------------------------------------------------------------*
 * getConcealedLP()
 *
 *
 *-------------------------------------------------------------------*/

void getConcealedLP(
    HANDLE_PLC_ENC_EVS memDecState,
    float *AqCon,
    const float lsfBase[],
    const float sr_core,
    int last_good,
    int L_frame
)
{
    float *lsf = memDecState->lsf_con;
    float lsp[(NB_DIV+1)*M];
    short k;

    dlpc_bfi( L_frame, &lsf[0], memDecState->lsfold, last_good,
              1 /* assumes packet loss */ , memDecState->mem_MA, memDecState->mem_AR, &(memDecState->stab_fac), memDecState->lsf_adaptive_mean,
              1, NULL, 0, NULL, NULL, lsfBase);

    mvr2r( memDecState->lspold, lsp, M );

    for ( k=0; k<1; k++ )
    {
        lsf2lsp(&lsf[k*M], &lsp[(k+1)*M], M, sr_core);
    }

    int_lsp( L_FRAME, &lsp[0], &lsp[M], AqCon, M, interpol_frac_12k8, 0 );

    return;
}

/*-------------------------------------------------------------------*
 * RecLpcSpecPowDiffuseLc()
 *
 *
 *-------------------------------------------------------------------*/

void RecLpcSpecPowDiffuseLc(
    float *lspq,
    float *lsp_old,
    float *lsfq,
    Decoder_State *st
    , int reset_q
)
{
    const float *means;
    float lsf_old[M];
    short i;


    means = PlcGetlsfBase ( st->lpcQuantization, st->narrowBand, st->sr_core );

    mvr2r( st->lsf_old, lsf_old, M );

    modify_lsf( lsf_old, M, st->sr_core
                , reset_q
              );

    lsf2lsp(lsf_old, lsp_old, M, st->sr_core);

    if (reset_q)
    {
        for (i=0; i<M; i++)
        {
            lsfq[i] = st->mem_MA[i] + means[i];
        }

        v_sort( lsfq, 0, M - 1 );

        reorder_lsfs( lsfq, LSF_GAP, M, st->sr_core);

        lsf2lsp(lsfq, lspq, M, st->sr_core);
    }
    else
    {
        modify_lsf( lsfq, M, st->sr_core, reset_q );

        lsf2lsp(lsfq, lspq, M, st->sr_core);
    }


    return;
}


/*-------------------------------------------------------------------*
 * modify_lsf()
 *
 *
 *-------------------------------------------------------------------*/

void modify_lsf(
    float *lsf,
    const short n,
    const int sr_core
    , int reset_q
)
{
    short i, k;
    float gap;
    float th;

    th = 1900;

    if (reset_q==0)
    {
        th = 800;
    }

    if(sr_core == 16000)
    {
        th *= 1.25;
    }


    i = 1;

    while(lsf[i] < th && i < n)
    {
        i++;
    }
    gap = lsf[i - 1] / i;

    for(k = 0; k < i - 1; k++)
    {
        lsf[k] = gap * (k + 1);
    }

    return;
}


/*-------------------------------------------------------------------*
 * reorder_lsfs()
 *
 *
 *-------------------------------------------------------------------*/

static void reorder_lsfs(
    float *lsf,      /* i/o: vector of lsfs in the frequency domain (0..0.5)*/
    float min_dist,  /* i  : minimum required distance */
    const short n,         /* i  : LPC order                 */
    int sr_core
)
{
    short i;
    float lsf_min;
    float lsf_max;
    float fac;
    float th1, th2;

    th1 = 1000.0f;
    th2 = 1900.0f;


    if(sr_core == 16000)
    {
        min_dist *= 1.25;
        th1 *= 1.25;
        th2 *= 1.25;
    }

    /*-----------------------------------------------------------------*
     * Verify the LSF ordering and minimum GAP
     *-----------------------------------------------------------------*/
    fac = 3.0;

    lsf_min = min_dist * fac;
    for (i = 0; i < n; i++)
    {

        if (lsf[i] > th1)
        {
            fac = 2.0;
        }
        else
        {
            if (lsf[i] > 1900.0)
            {
                fac = 1.0;
            }
        }

        if (lsf[i] < lsf_min)
        {
            lsf[i] = lsf_min;
        }

        lsf_min = lsf[i] + min_dist * fac;
    }

    /*------------------------------------------------------------------------------------------*
     * Reverify the LSF ordering and minimum GAP in the reverse order (security)
     *------------------------------------------------------------------------------------------*/

    lsf_max = (float)(sr_core)/2.0f - min_dist * fac;

    if( lsf[n-1] > lsf_max )        /* If danger of unstable filter in case of resonance in HF */
    {
        for (i = n-1; i >=0; i--)   /* Reverify the minimum LSF gap in the reverse sens */
        {

            if (lsf[i] <= th2)
            {
                fac = 2.0;
            }
            else
            {
                if (lsf[i] <= th1)
                {
                    fac = 3.0;
                }
            }

            if (lsf[i] > lsf_max)
            {
                lsf[i] = lsf_max;
            }

            lsf_max = lsf[i] - min_dist * fac;
        }
    }

    return;
}
