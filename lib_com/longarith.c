/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdio.h>
#include <assert.h>
#include "prot.h"




/**
 * \brief          long addition: a[] = a[] + b[]
 *                 Addition of unsigned vectors a[] and b[] without saturation.
 *                 All words of b[] are added with carry to the corresponding words in a.
 *                 An assertion failure occurs, if lenb exceeds lena or if overflow occurs.
 *
 * \param          unsigned short a[]
 *                 Input/Output: vector of the length lena
 * \param          unsigned short b[]
 *                 Input: vector of the length lenb
 * \param          int lena
 *                 Input: length of vector a[] in units of 'unsigned short'
 *                 Note:  lena must be greater/equal lenb
 * \param          int lenb
 *                 Input: length of vector b[] in units of 'unsigned short'
 *
 * \return         void
 */

void longadd(unsigned short a[], unsigned short b[], int lena, int lenb)
{
    int h;
    long carry = 0;

    assert(lena >= lenb);
    for (h=0; h < lenb; h++)
    {
        carry += ((unsigned long)a[h]) + ((unsigned long)b[h]);
        a[h]  = (unsigned short) carry;
        carry = carry >> 16;
    }
    for (; h < lena; h++)
    {
        carry = ((unsigned long)a[h]) + carry;
        a[h] = (unsigned short) carry;
        carry = carry >> 16;
    }

    assert(carry == 0);   /* carry != 0 indicates addition overflow */
    return;
}

/**
 * \brief          long shift right: d[] = a[] >> b
 *                 Logical shift right of unsigned vectors a[] by b bit-positions.
 *                 Vector d[] is filled with leading zeroes, where lend exceeds lena.
 *                 It is allowed to overlay d[] with a[].
 *
 * \param          unsigned short a[]
 *                 Input: vector of the length lena
 * \param          int b
 *                 Input: number of bit positions to shift right
 * \param          unsigned short d[]
 *                 Output: vector of the length lend
 *                 Note:   vector d[] can be overlaid with vector a[], i.e. a[] = a[] >> b
 * \param          int lena
 *                 Input: length of vector a[] in units of 'unsigned short'
 * \param          int lend
 *                 Input: length of vector d[] in units of 'unsigned short'
 *
 * \return         void
 */

void longshiftright(unsigned short a[], int b, unsigned short d[], int lena, int lend)
{
    int intb, fracb, fracb_u, k;


    intb = b >> 4;

    a += intb;
    lena -= intb;

    fracb = b & 0xF;
    if (fracb)
    {
        fracb_u = 16-fracb;
        for (k=0; k < lena-1; k++)
        {
            d[k] = ((a[k] >> fracb) | (a[k+1] << fracb_u)) & 0xFFFF;
        }
        d[k] = (a[k] >> fracb);
        k++;
    }
    else
    {
        for (k=0; k < lena; k++)
        {
            d[k] = a[k];
        }
    }
    /* fill remaining upper bits with zero */
    for (; k < lend; k++)
    {
        d[k] = 0;
    }
    return;
}

/**
 * \brief          long shift left: d[] = a[] << b
 *                 Logical shift left of unsigned vectors a[] by b bit-positions.
 *                 It is allowed to overlay d[] with a[].
 *
 * \param          unsigned short a[]
 *                 Input: vector of the length len
 * \param          int b
 *                 Input: number of bit positions to shift left
 * \param          unsigned short d[]
 *                 Output: vector of the length len
 *                 Note:   vector d[] can be overlaid with vector a[], i.e. a[] = a[] << b
 * \param          int len
 *                 Input: length of vector a[] and d[] in units of 'unsigned short'
 *
 * \return         void
 */

void longshiftleft(unsigned short a[], int b, unsigned short d[], int len)
{
    int intb;      /* integer part of b */
    int fracb;     /* shift left value for all upper words a[k] */
    int fracb_l;   /* shift right value for all lower words a[k-1] */
    int k = len - 1;


    intb  = b >> 4;
    fracb = b & 0xF;

    if (fracb)
    {
        fracb_l = 16-fracb;
        for (; k >intb; k--)
        {
            d[k] = (a[k-intb] << fracb) | (a[k-intb-1] >> fracb_l);
        }
        d[k] = (a[k-intb] << fracb);
        k--;
    }
    else
    {
        for (; k >= intb; k--)
        {
            d[k] = a[k-intb];
        }
    }
    for ( ; k >= 0; k--)
    {
        d[k] = 0;
    }
    return;
}
