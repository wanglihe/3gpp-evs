/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "prot.h"

/*-------------------------------------------------------------------*
 * Function prototypes
 *-------------------------------------------------------------------*/
void re8_decode_base_index(int n, long I, int *x);

/*--------------------------------------------------------------------------
  re8_dec(n, I, k, y)
  MULTI-RATE INDEXING OF A POINT y in THE LATTICE RE8 (INDEX DECODING)
  (i) n: codebook number (*n is an integer defined in {0,2,3,4,..,n_max})
  (i) I: index of c (pointer to unsigned 16-bit word)
  (i) k: index of v (8-dimensional vector of binary indices) = Voronoi index
  (o) y: point in RE8 (8-dimensional integer vector)
  note: the index I is defined as a 32-bit word, but only
  16 bits are required (long can be replaced by unsigned integer)
  --------------------------------------------------------------------------
 */
void re8_dec(int n, long I, int k[], int y[])
{
    int i, m, v[8];


    /*------------------------------------------------------------------------*
     * decode the sub-indices I and kv[] according to the codebook number n:
     *  if n=0,2,3,4, decode I (no Voronoi extension)
     *  if n>4, Voronoi extension is used, decode I and kv[]
     *------------------------------------------------------------------------*/
    if (n <= 4)
    {
        re8_decode_base_index(n, I, y);
    }
    else
    {
        /*--------------------------------------------------------------------*
         * compute the Voronoi modulo m = 2^r where r is extension order
         *--------------------------------------------------------------------*/
        m = 1;
        while (n > 4)
        {
            m *= 2;
            n -= 2;
        }

        /*--------------------------------------------------------------------*
         * decode base codebook index I into c (c is an element of Q3 or Q4)
         *  [here c is stored in y to save memory]
         *--------------------------------------------------------------------*/

        re8_decode_base_index(n, I, y);

        /*--------------------------------------------------------------------*
         * decode Voronoi index k[] into v
         *--------------------------------------------------------------------*/
        re8_k2y(k, m, v);

        /*--------------------------------------------------------------------*
         * reconstruct y as y = m c + v (with m=2^r, r integer >=1)
         *--------------------------------------------------------------------*/
        for (i=0; i<8; i++)
        {
            y[i] = m*y[i] + v[i];
        }
    }

    return;
}
