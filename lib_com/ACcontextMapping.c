/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "prot.h"


/*-------------------------------------------------------------------*
* get_next_coeff_mapped()
*
*
*-------------------------------------------------------------------*/

int get_next_coeff_mapped(  /* o  : index of next coefficient */
    int ii[2],                /* i/o: coefficient indexes       */
    int *pp,                  /* o  : peak(1)/hole(0) indicator */
    int *idx,                 /* o  : index in unmapped domain  */
    CONTEXT_HM_CONFIG *hm_cfg /* i  : HM configuration          */
)
{
    unsigned int p;

    p = (ii[1] - hm_cfg->numPeakIndices) & (hm_cfg->indexBuffer[ii[1]] - hm_cfg->indexBuffer[ii[0]]);
    p >>= sizeof(p)*8-1;
    *pp = p;
    *idx = ii[p];
    ii[p]++;

    return hm_cfg->indexBuffer[*idx];
}


/*-------------------------------------------------------------------*
* get_next_coeff_unmapped()
*
*
*-------------------------------------------------------------------*/

int get_next_coeff_unmapped(/* o  : index of next coefficient */
    int *ii,                  /* i/o: coefficient indexes       */
    int *idx                  /* o  : index in unmapped domain  */
)
{
    *idx = *ii;
    (*ii)++;

    return *idx;
}

/*-------------------------------------------------------------------*
* update_mixed_context()
*
*
*-------------------------------------------------------------------*/

int update_mixed_context(
    int ctx,
    int a
)
{
    int t;

    t = 1-13 + (a & ~1)*((a>>2)+1);

    if (t > 0)
    {
        t = min((a>>3), 2);
    }

    return (ctx & 0xf) * 16 + t + 13;
}
