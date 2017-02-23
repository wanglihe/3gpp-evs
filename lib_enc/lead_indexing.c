/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <math.h>
#include "options.h"
#include "rom_com.h"
#include "prot.h"


/*-------------------------------------------------------------------*
 * Local function prototypes
 *-------------------------------------------------------------------*/
static int fcb_encode_pos(int pos_vector[], int pulse_num, int pos_num);

/*-------------------------------------------------------------------*
 * re8_compute_base_index:
 *
 * Compute base index for RE8
 *-------------------------------------------------------------------*/

void re8_compute_base_index(
    const int   *x,   /* i  : Elemen of Q2, Q3 or Q4                         */
    const int   ka,   /* i  : Identifier of the absolute leader related to x */
    long  *I    /* o  : index                                          */
)
{
    int i, j, k1,m;
    int setor_8p[8], setor_8p_temp[8];
    int sign_8p;
    int code_index, code_level, code_area;
    int offset;
    const int *a1,*a2;


    a1 = vals_a[ka];
    a2 = vals_q[ka];

    /* the sign process */

    sign_8p = 0;
    m = 0;
    code_index = 0;
    k1 = a2[0];

    if ((a2[1]==2) && (a1[0] ^ 1) && (ka != 5))
    {
        for (i=0; i<8; i++)
        {
            if (x[i] != 0)
            {
                sign_8p = sign_8p * 2;
                setor_8p_temp[m] = i;
                m++;
            }
            if (x[i] < 0)
            {
                sign_8p += 1;
            }
        }

        code_index = fcb_encode_pos(setor_8p_temp,8,m);
        code_index = (code_index << k1) + sign_8p;

        offset = Is[ka];


        *I = offset + code_index;
    }
    else
    {
        for (i=0; i<8; i++)
        {
            setor_8p[i] = abs(x[i]);

            if (x[i] < 0)
            {
                sign_8p = sign_8p * 2 + 1;
                m++;
            }
            if (x[i] > 0)
            {
                sign_8p = sign_8p * 2;
                m++;
            }
        }

        if ( k1 != m )
        {
            sign_8p = sign_8p >> 1;
        }

        /* code level by level */

        code_level = a2[1];
        code_area = 8;

        if ( a2[2] != 1 )
        {
            for (j=0; j<code_level-1; j++)
            {
                m = 0;

                for (i=0; i<code_area; i++)
                {
                    if (setor_8p[i] != a1[j])
                    {
                        setor_8p_temp[m] = i;
                        setor_8p[m] = setor_8p[i];
                        m++;
                    }
                }
                code_index *= select_table22[m][code_area];

                code_index += fcb_encode_pos(setor_8p_temp,code_area,m);

                code_area = m;
            }
        }
        else
        {
            for (i=0; i<code_area; i++)
            {
                if (setor_8p[i] == a1[1])
                {
                    code_index += i;
                }
            }
        }

        code_index = ( code_index << k1 ) + sign_8p;

        offset = Is[ka];

        *I = offset + code_index;
    }

    return;
}

/*-------------------------------------------------------------------*
 * fcb_encode_pos:
 *
 * Base function to compute base index for RE8
 *-------------------------------------------------------------------*/

static int fcb_encode_pos(      /* o  : Code index              */
    int pos_vector[],       /* i  : Position vectort        */
    int pulse_num,          /* i  : Pulse number            */
    int pos_num             /* i  : Position number         */
)
{
    int i,j;
    int code_index;
    int temp,temp1;
    const int *select_table23;


    temp = pulse_num-1;

    select_table23 = select_table22[pos_num];

    code_index = select_table23[pulse_num] - select_table23[pulse_num-pos_vector[0]];

    for (i=0,j=1; i<pos_num-1; i++,j++)
    {
        temp1 = pos_num-j;

        select_table23 = select_table22[temp1];

        code_index += select_table23[temp-pos_vector[i]] - select_table23[pulse_num-pos_vector[j]] ;
    }

    return code_index;
}
