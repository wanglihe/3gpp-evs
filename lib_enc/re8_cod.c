/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "prot.h"


/*------------------------------------------------------------------------
 * re8_cod()
 *
 * MULTI-RATE INDEXING OF A POINT y in THE LATTICE RE8 (INDEX COMPUTATION)
 * note: the index I is defined as a 32-bit word, but only
 *       16 bits are required (long can be replaced by unsigned integer)
 *--------------------------------------------------------------------------*/

void re8_cod(
    int x[],  /* i  : point in RE8 (8-dimensional integer vector)                         */
    int *n,   /* i  : codebook number (*n is an integer defined in {0,2,3,4,..,n_max})    */
    long *I,  /* o  : index of c (pointer to unsigned 16-bit word)                        */
    int k[]   /* o  : index of v (8-dimensional vector of binary indices) = Voronoi index */
)
{
    int ka, c[8];

    /*----------------------------------------------------------------------
     * decompose x as x = 2^r c + v, where r is an integer >=0, c is an element
     *  of Q0, Q2, Q3 or Q4, and v is an element of a Voronoi code in RE8
     *  (if r=0, x=c)
     *  this decomposition produces as a side-product the index k[] of v
     *  and the identifier ka of the absolute leader related to c
     *
     *  the index of y is split into 2 parts :
     *  - the index I of c
     *  - the index k[] of v
     ----------------------------------------------------------------------*/

    re8_vor( x, n, k, c, &ka );

    /* compute the index I (only if c is in Q2, Q3 or Q4) */
    if( *n > 0 )
    {
        re8_compute_base_index( c, ka, I );
    }

    return;
}
