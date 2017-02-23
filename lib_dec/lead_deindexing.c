/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"
#include "rom_dec.h"

/*-------------------------------------------------------------------*
 * Local function prototype
 *-------------------------------------------------------------------*/
static void fcb_decode_pos(int index, int pos_vector[], int pulse_num, int pos_num);

/*-------------------------------------------------------------------*
 * re8_decode_base_index
 *
 * Decode RE8 base index
 *-------------------------------------------------------------------*/

void re8_decode_base_index(
    int n,    /* i  : codebook number (*n is an integer defined in {0,2,3,4,..,n_max}) */
    long I,   /* i  : index of c (pointer to unsigned 16-bit word)                     */
    int *x    /* o  : point in RE8 (8-dimensional integer vector)                      */
)
{
    int i,j,k1,l,m,m1,m2;
    int setor_8p_temp[8],setor_8p_temp_1[8],setor_8p_temp_2[8]= {0};
    int sign_8p;
    int code_level;
    const int *a1,*a2;

    int ka;
    int offset;
    int code_index;
    int element_a10,element_a11,element_a12;


    element_a11 = 0,
    element_a12 = 0;
    if (n < 2)
    {
        for (i=0; i<8; i++)
        {
            x[i]=0;
        }
    }
    else
    {
        if ( I > 65519L )
        {
            I = 0;
        }
        /*-------------------------------------------------------------------*
         * search for the identifier ka of the absolute leader (table-lookup)
         * Q2 is a subset of Q3 - the two cases are considered in the same branch
         *-------------------------------------------------------------------*/
        if ( n <= 3 )
        {
            for (i=1; i<NB_LDQ3; i++)
            {
                if (I<(long)II3[i])
                    break;
            }
            ka = AA3[i-1];
        }
        else
        {
            for (i=1; i<NB_LDQ4; i++)
            {
                if (I<(long)II4[i])
                    break;
            }
            ka = AA4[i-1];
        }

        /*-------------------------------------------------------*
         * decode
         *-------------------------------------------------------*/

        a1 = vals_a[ka];
        a2 = vals_q[ka];
        k1 = a2[0];
        code_level = a2[1];

        offset = Is[ka];
        code_index = I - offset;

        sign_8p = code_index & ( (1<<k1) -1 );

        code_index = code_index >> k1 ;

        m  = 0;
        m1 = 0;
        m2 = 0;

        element_a10 = a1[0];
        switch (code_level)
        {
        case 4:

            m2 = 1; /*a2[4];*/
            i = code_index & 1;
            if( i == 0)
                setor_8p_temp_2[0] = 0;
            else
                setor_8p_temp_2[0] = 1;
            code_index = code_index >> 1;

        case 3:

            m  = a2[2];
            m1 = a2[3];
            l = select_table22[m1][m];
            j = ( code_index * mult_avq_tab[l] ) >> shift_avq_tab[l];
            code_index = code_index - j * l;
            fcb_decode_pos(code_index,setor_8p_temp_1,m,m1);
            code_index = j;
            element_a12 = a1[2];

        case 2:

            m  = a2[2];
            fcb_decode_pos(code_index,setor_8p_temp,8,m);
            element_a11 = a1[1];
        }

        for (i=0; i<8; i++)
        {
            x[i] = element_a10;
        }

        for (i=0; i<m; i++)
        {
            x[setor_8p_temp[i]] = element_a11;
        }

        for (i=0; i<m1; i++)
        {
            x[setor_8p_temp[setor_8p_temp_1[i]]] = element_a12;
        }

        for (i=0; i<m2; i++)
        {
            x[setor_8p_temp[setor_8p_temp_1[setor_8p_temp_2[0]]]] = 6;
        }

        /*--------------------------------------------------------------------*
         * add the sign of all elemnt ( except the last one in some case )
         *--------------------------------------------------------------------*/
        m1 = k1-1;

        for (i=0; i<8; i++)
        {
            if ( x[i] != 0)
            {
                if ( ( sign_8p >> m1 ) & 1 )
                {
                    x[i] *= -1;
                }
                m1 --;
            }
        }

        /*--------------------------------------------------------------------*
         * recover the sign of last element if needed
         *--------------------------------------------------------------------*/
        if ( k1 == 7 )
        {
            m1 = 0;

            for (i=0; i<8; i++)
            {
                m1 += x[i] ;
            }
            if ( m1 & 3 )
            {
                x[7] *= -1;
            }
        }
    }

    return;
}


/*-------------------------------------------------------------------*
 * fcb_decode_pos
 *
 * base function for decoding position index
 *-------------------------------------------------------------------*/

void fcb_decode_pos(
    int index,        /* i  : Index to decoder    */
    int pos_vector[], /* o  : Position vector     */
    int pulse_num,    /* i  : Number of pulses    */
    int pos_num       /* i  : Number of positions */
)
{
    int i,k,l;
    int temp1,temp2;

    const int *select_table23;
    const int *select_table24;

    k = index;
    l = 0;
    temp1 = pos_num;
    temp2 = pulse_num+1;

    for( i=0; i<pos_num-1; i++ )
    {
        select_table23 = select_table22[temp1];
        select_table24 = &select_table23[pulse_num-l];

        k = *select_table24 - k;

        while ( k <= (*select_table24--) )
        {
            l++;
        }
        k = select_table23[temp2-l] - k;

        pos_vector[i] = l - 1;
        temp1--;
    }

    pos_vector[i] = l+k;

    return;
}
