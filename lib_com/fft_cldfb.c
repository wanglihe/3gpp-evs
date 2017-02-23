/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <string.h>
#include <assert.h>
#include "prot.h"



#if defined __ICL
#define restrict __restrict
#else
#define restrict
#endif

#ifdef _MSC_VER
#pragma warning(disable : 4305) /* disable truncation from double to float warning (VC++)*/
#endif

static void fft8(float *vec);
static void fft10(float *vec);
static void fft16(float *vec);
static void fft20(float *vec);
static void fft30(float *vec);
static void fft5s(float *x, int stride);


static const float INV_SQRT2    = 7.071067811865475e-1;
static const float COS_PI_DIV8  = 9.238795325112867e-1;
static const float COS_3PI_DIV8 = 3.826834323650898e-1;
static const float SQRT2PLUS1   = 2.414213562373095;
static const float SQRT2MINUS1  = 4.142135623730952e-1;


#ifdef _MSC_VER
#pragma warning(default : 4305)
#endif

/*******************************************************************************
 Functionname:  fft8
 *******************************************************************************

 Description:   8-point FFT. Complex-valued input takes 52 real additions
                and 4 real multiplications.

 Arguments:     vec - pointer to data (interleaved real / imaginary parts)

 Return:        none

*******************************************************************************/
static void fft8(float * restrict vec)
{
    float temp1[16];
    float temp2[16];


    /* Pre-additions */
    temp1[0]  = vec[0] + vec[8];
    temp1[2]  = vec[0] - vec[8];
    temp1[1]  = vec[1] + vec[9];
    temp1[3]  = vec[1] - vec[9];
    temp1[4]  = vec[2] + vec[10];
    temp1[6]  = vec[2] - vec[10];
    temp1[5]  = vec[3] + vec[11];
    temp1[7]  = vec[3] - vec[11];
    temp1[8]  = vec[4] + vec[12];
    temp1[10] = vec[4] - vec[12];
    temp1[9]  = vec[5] + vec[13];
    temp1[11] = vec[5] - vec[13];
    temp1[12] = vec[6] + vec[14];
    temp1[14] = vec[6] - vec[14];
    temp1[13] = vec[7] + vec[15];
    temp1[15] = vec[7] - vec[15];

    /* Pre-additions and core multiplications */
    temp2[0]  =  temp1[0] + temp1[8];
    temp2[4]  =  temp1[0] - temp1[8];
    temp2[1]  =  temp1[1] + temp1[9];
    temp2[5]  =  temp1[1] - temp1[9];
    temp2[8]  =  temp1[2] - temp1[11];
    temp2[10] =  temp1[2] + temp1[11];
    temp2[9]  =  temp1[3] + temp1[10];
    temp2[11] =  temp1[3] - temp1[10];
    temp2[2]  =  temp1[4] + temp1[12];
    temp2[7]  =  temp1[4] - temp1[12];
    temp2[3]  =  temp1[5] + temp1[13];
    temp2[6]  =  temp1[13]- temp1[5];

    temp1[1]  =  temp1[6] + temp1[14];
    temp1[2]  =  temp1[6] - temp1[14];
    temp1[0]  =  temp1[7] + temp1[15];
    temp1[3]  =  temp1[7] - temp1[15];

    temp2[12] = (temp1[0] + temp1[2]) *  INV_SQRT2;
    temp2[14] = (temp1[0] - temp1[2]) *  INV_SQRT2;
    temp2[13] = (temp1[3] - temp1[1]) *  INV_SQRT2;
    temp2[15] = (temp1[1] + temp1[3]) * -INV_SQRT2;

    /* Post-additions */
    vec[0]  = temp2[0] + temp2[2];
    vec[8]  = temp2[0] - temp2[2];
    vec[1]  = temp2[1] + temp2[3];
    vec[9]  = temp2[1] - temp2[3];
    vec[4]  = temp2[4] - temp2[6];
    vec[12] = temp2[4] + temp2[6];
    vec[5]  = temp2[5] - temp2[7];
    vec[13] = temp2[5] + temp2[7];
    vec[6]  = temp2[8] + temp2[14];
    vec[14] = temp2[8] - temp2[14];
    vec[7]  = temp2[9] + temp2[15];
    vec[15] = temp2[9] - temp2[15];
    vec[2]  = temp2[10]+ temp2[12];
    vec[10] = temp2[10]- temp2[12];
    vec[3]  = temp2[11]+ temp2[13];
    vec[11] = temp2[11]- temp2[13];

}



/*******************************************************************************
 Functionname:  fft16
 *******************************************************************************

 Description:   16-point FFT. Complex-valued input takes 144 real additions and
                24 real multiplications.

 Arguments:     vec - pointer to data (interleaved real / imaginary parts)

 Return:        none

*******************************************************************************/
/* fast implementation, completely unrolled and inlined */
static void fft16(float * restrict vec)
{
    float temp10, temp11, temp12, temp13, temp14, temp15, temp16, temp17,
          temp18, temp19, temp110, temp111, temp112, temp113, temp114, temp115;
    float temp20, temp21, temp22, temp23, temp24, temp25, temp26, temp27,
          temp28, temp29, temp210, temp211, temp212, temp213, temp214, temp215;
    float vec0, vec1, vec2, vec3, vec4, vec5, vec6, vec7,
          vec8, vec9, vec10, vec11, vec12, vec13, vec14, vec15;


    /* even */
    vec0 = vec[0] + vec[16];
    vec1 = vec[1] + vec[17];
    vec2 = vec[2] + vec[18];
    vec3 = vec[3] + vec[19];
    vec4 = vec[4] + vec[20];
    vec5 = vec[5] + vec[21];
    vec6 = vec[6] + vec[22];
    vec7 = vec[7] + vec[23];
    vec8 = vec[8] + vec[24];
    vec9 = vec[9] + vec[25];
    vec10 = vec[10] + vec[26];
    vec11 = vec[11] + vec[27];
    vec12 = vec[12] + vec[28];
    vec13 = vec[13] + vec[29];
    vec14 = vec[14] + vec[30];
    vec15 = vec[15] + vec[31];

    /* Pre-additions */
    temp10  = vec0 + vec8;
    temp12  = vec0 - vec8;
    temp11  = vec1 + vec9;
    temp13  = vec1 - vec9;
    temp14  = vec2 + vec10;
    temp16  = vec2 - vec10;
    temp15  = vec3 + vec11;
    temp17  = vec3 - vec11;
    temp18  = vec4 + vec12;
    temp110 = vec4 - vec12;
    temp19  = vec5 + vec13;
    temp111 = vec5 - vec13;
    temp112 = vec6 + vec14;
    temp114 = vec6 - vec14;
    temp113 = vec7 + vec15;
    temp115 = vec7 - vec15;

    /* Pre-additions and core multiplications */
    temp20  =  temp10 + temp18;
    temp24  =  temp10 - temp18;
    temp21  =  temp11 + temp19;
    temp25  =  temp11 - temp19;
    temp28  =  temp12 - temp111;
    temp210 =  temp12 + temp111;
    temp29  =  temp13 + temp110;
    temp211 =  temp13 - temp110;
    temp22  =  temp14 + temp112;
    temp27  =  temp14 - temp112;
    temp23  =  temp15 + temp113;
    temp26  =  temp113- temp15;

    temp11  =  temp16 + temp114;
    temp12  =  temp16 - temp114;
    temp10  =  temp17 + temp115;
    temp13  =  temp17 - temp115;

    temp212 = (temp10 + temp12) *  INV_SQRT2;
    temp214 = (temp10 - temp12) *  INV_SQRT2;
    temp213 = (temp13 - temp11) *  INV_SQRT2;
    temp215 = (temp11 + temp13) * -INV_SQRT2;



    /* odd */
    vec0 = vec[0] - vec[16];
    vec1 = vec[1] - vec[17];
    vec2 = vec[2] - vec[18];
    vec3 = vec[3] - vec[19];
    vec4 = vec[4] - vec[20];
    vec5 = vec[5] - vec[21];
    vec6 = vec[6] - vec[22];
    vec7 = vec[7] - vec[23];
    vec8 = vec[8] - vec[24];
    vec9 = vec[9] - vec[25];
    vec10 = vec[10] - vec[26];
    vec11 = vec[11] - vec[27];
    vec12 = vec[12] - vec[28];
    vec13 = vec[13] - vec[29];
    vec14 = vec[14] - vec[30];
    vec15 = vec[15] - vec[31];

    /* Pre-additions and core multiplications */
    temp19  = (vec2 + vec14) * -COS_3PI_DIV8;
    temp110 = (vec2 - vec14) * COS_PI_DIV8;
    temp18  = (vec3 + vec15) * COS_3PI_DIV8;
    temp111 = (vec3 - vec15) *  COS_PI_DIV8;
    temp15  = (vec4 + vec12) * -INV_SQRT2;
    temp16  = (vec4 - vec12) * INV_SQRT2;
    temp14  = (vec5 + vec13) * INV_SQRT2;
    temp17  = (vec5 - vec13) *  INV_SQRT2;
    temp113 = (vec6 + vec10) * -COS_PI_DIV8;
    temp114 = (vec6 - vec10) * COS_3PI_DIV8;
    temp112 = (vec7 + vec11) * COS_PI_DIV8;
    temp115 = (vec7 - vec11) *  COS_3PI_DIV8;

    /* Core multiplications */
    vec2 = temp18  * SQRT2PLUS1  - temp112 * SQRT2MINUS1;
    vec3 = temp19  * SQRT2PLUS1  - temp113 * SQRT2MINUS1;
    vec4 = temp110 * SQRT2MINUS1 - temp114 * SQRT2PLUS1;
    vec5 = temp111 * SQRT2MINUS1 - temp115 * SQRT2PLUS1;

    /* Post-additions */
    temp18 += temp112;
    temp19 += temp113;
    temp110+= temp114;
    temp111+= temp115;

    vec6   = vec0  + temp14;
    vec10  = vec0  - temp14;
    vec7   = vec1  + temp15;
    vec11  = vec1  - temp15;

    vec12  = temp16 - vec9;
    vec14  = temp16 + vec9;
    vec13  = vec8  + temp17;
    vec15  = vec8  - temp17;

    temp10  = vec6  - vec14;
    temp12  = vec6  + vec14;
    temp11  = vec7  + vec15;
    temp13  = vec7  - vec15;
    temp14  = vec10 + vec12;
    temp16  = vec10 - vec12;
    temp15  = vec11 + vec13;
    temp17  = vec11 - vec13;

    vec10  = temp18 + temp110;
    temp110 = temp18 - temp110;
    vec11  = temp19 + temp111;
    temp111 = temp19 - temp111;

    temp112 = vec2  + vec4;
    temp114 = vec2  - vec4;
    temp113 = vec3  + vec5;
    temp115 = vec3  - vec5;


    /* Post-additions */
    *vec++ = temp20 + temp22;
    *vec++ = temp21 + temp23;
    *vec++ = temp12 + vec10;
    *vec++ = temp13 + vec11;
    *vec++ = temp210+ temp212;
    *vec++ = temp211+ temp213;
    *vec++ = temp10 + temp112;
    *vec++ = temp11 + temp113;
    *vec++ = temp24 - temp26;
    *vec++ = temp25 - temp27;
    *vec++ = temp16 + temp114;
    *vec++ = temp17 + temp115;
    *vec++ = temp28 + temp214;
    *vec++ = temp29 + temp215;
    *vec++ = temp14 + temp110;
    *vec++ = temp15 + temp111;
    *vec++ = temp20 - temp22;
    *vec++ = temp21 - temp23;
    *vec++ = temp12 - vec10;
    *vec++ = temp13 - vec11;
    *vec++ = temp210- temp212;
    *vec++ = temp211- temp213;
    *vec++ = temp10 - temp112;
    *vec++ = temp11 - temp113;
    *vec++ = temp24 + temp26;
    *vec++ = temp25 + temp27;
    *vec++ = temp16 - temp114;
    *vec++ = temp17 - temp115;
    *vec++ = temp28 - temp214;
    *vec++ = temp29 - temp215;
    *vec++ = temp14 - temp110;
    *vec++ = temp15 - temp111;

}


/*******************************************************************************
 Functionname:  fft15
 *******************************************************************************

 Description:   15-point FFT. Complex-valued input takes 176 real additions
                and 34 real multiplications.

 Arguments:     vec - pointer to data (interleaved real / imaginary parts)

 Return:        none

*******************************************************************************/
static void fft15(float * restrict vec)
{

    float r0, r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13, r14, r15, r16, r17;
    float i0, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, i12, i13, i14, i15, i16, i17;
    float tmp0, tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9,
          tmp10, tmp11, tmp12, tmp13, tmp14, tmp15, tmp16, tmp17, tmp18, tmp19,
          tmp20, tmp21, tmp22, tmp23, tmp24, tmp25, tmp26, tmp27, tmp28, tmp29;


    /* Pre-additions real part */
    r1  = vec[2]  + vec[8];
    r2  = vec[2]  - vec[8];
    r3  = vec[4]  + vec[16];
    r4  = vec[4]  - vec[16];
    r5  = vec[6]  + vec[24];
    r6  = vec[6]  - vec[24];
    r7  = vec[10] + vec[20];
    r8  = vec[10] - vec[20];
    r9  = vec[12] + vec[18];
    r10 = vec[12] - vec[18];
    r11 = vec[14] + vec[26];
    r12 = vec[14] - vec[26];
    r13 = vec[22] + vec[28];
    r14 = vec[22] - vec[28];

    tmp2  = r1  + r3;
    tmp4  = r1  - r3;
    tmp6  = r2  + r14;
    tmp8  = r2  - r14;
    tmp10 = r4  + r12;
    tmp12 = r4  - r12;
    tmp14 = r5  + r9;
    tmp16 = r5  - r9;
    tmp18 = r11 + r13;
    tmp20 = r11 - r13;

    /* Pre-additions imaginary part */
    i1  = vec[3]  + vec[9];
    i2  = vec[3]  - vec[9];
    i3  = vec[5]  + vec[17];
    i4  = vec[5]  - vec[17];
    i5  = vec[7]  + vec[25];
    i6  = vec[7]  - vec[25];
    i7  = vec[11] + vec[21];
    i8  = vec[11] - vec[21];
    i9  = vec[13] + vec[19];
    i10 = vec[13] - vec[19];
    i11 = vec[15] + vec[27];
    i12 = vec[15] - vec[27];
    i13 = vec[23] + vec[29];
    i14 = vec[23] - vec[29];

    tmp3  = i1  + i3;
    tmp5  = i1  - i3;
    tmp7  = i2  + i14;
    tmp9  = i2  - i14;
    tmp11 = i4  + i12;
    tmp13 = i4  - i12;
    tmp15 = i5  + i9;
    tmp17 = i5  - i9;
    tmp19 = i11 + i13;
    tmp21 = i11 - i13;


    /* Pre-additions and core multiplications */
    tmp28=  tmp4  + tmp20;
    tmp29=  tmp5  + tmp21;
    r4   =  tmp2  + tmp18;
    i4   =  tmp3  + tmp19;
    r3   = (r4    + tmp14) * -1.25f;
    i3   = (i4    + tmp15) * -1.25f;
    r2   = (tmp29 - i8)    * -8.660254037844387e-1f;
    i2   = (tmp28 - r8)    *  8.660254037844387e-1f;
    r1   =  r4    + r7;
    i1   =  i4    + i7;
    r0   =  r1    + vec[0]   +  tmp14;
    i0   =  i1    + vec[1]   +  tmp15;
    r7   =  tmp4  - tmp20;
    i7   =  tmp5  - tmp21;
    r8   = (tmp3  - tmp19) * -4.841229182759272e-1f;
    i8   = (tmp2  - tmp18) *  4.841229182759272e-1f;
    tmp0 =  tmp6  + r10;
    tmp1 =  tmp7  + i10;
    tmp2 =  r6    - tmp10;
    tmp3 =  i6    - tmp11;
    r10  =  tmp7  * -2.308262652881440f;
    i10  =  tmp6  *  2.308262652881440f;
    r11  =  tmp8  *  1.332676064001459f;
    i11  =  tmp9  *  1.332676064001459f;
    r6   = (r7    - tmp16) *  5.590169943749475e-1f;
    i6   = (i7    - tmp17) *  5.590169943749475e-1f;
    r12  = (tmp1  + tmp3)  *  5.877852522924733e-1f;
    i12  = (tmp0  + tmp2)  * -5.877852522924733e-1f;
    r13  = (tmp7  - tmp11) * -8.816778784387098e-1f;
    i13  = (tmp6  - tmp10) *  8.816778784387098e-1f;
    r14  = (tmp8  + tmp12) *  5.090369604551274e-1f;
    i14  = (tmp9  + tmp13) *  5.090369604551274e-1f;
    r16  =  tmp11 *  5.449068960040204e-1f;
    i16  =  tmp10 * -5.449068960040204e-1f;
    r17  =  tmp12 *  3.146021430912046e-1f;
    i17  =  tmp13 *  3.146021430912046e-1f;

    r4 *=  1.875f;
    i4 *=  1.875f;
    r1 *= -1.5f;
    i1 *= -1.5f;
    r7 *= -8.385254915624212e-1f;
    i7 *= -8.385254915624212e-1f;
    r5  = tmp29 *  1.082531754730548f;
    i5  = tmp28 * -1.082531754730548f;
    r9  = tmp1  *  1.538841768587627f;
    i9  = tmp0  * -1.538841768587627f;
    r15 = tmp3  *  3.632712640026803e-1f;
    i15 = tmp2  * -3.632712640026803e-1f;


    /* Post-additions real part */
    tmp2  = r0  + r1;
    tmp4  = r3  + r6;
    tmp6  = r3  - r6;
    tmp8  = r4  + r5;
    tmp10 = r4  - r5;
    tmp12 = r7  + r8;
    tmp14 = r7  - r8;
    tmp16 = r13 + r16;
    tmp18 = r14 + r17;
    tmp20 = r10 - r13;
    tmp22 = r11 - r14;
    tmp24 = r12 + r15;
    tmp26 = r12 - r9;

    r1  = tmp2  + r2;
    r2  = tmp2  - r2;
    r3  = tmp4  + tmp26;
    r4  = tmp4  - tmp26;
    r5  = tmp6  + tmp24;
    r6  = tmp6  - tmp24;
    r7  = tmp16 + tmp18;
    r8  = tmp16 - tmp18;
    r9  = tmp20 - tmp22;
    r10 = tmp20 + tmp22;
    r11 = r1    + tmp8;
    r12 = r2    + tmp10;
    r13 = r11   - tmp12;
    r14 = r12   - tmp14;
    r15 = r12   + tmp14;
    r16 = r11   + tmp12;

    /* Post-additions imaginary part */
    tmp3  = i0  + i1;
    tmp5  = i3  + i6;
    tmp7  = i3  - i6;
    tmp9  = i4  + i5;
    tmp11 = i4  - i5;
    tmp13 = i7  + i8;
    tmp15 = i7  - i8;
    tmp17 = i13 + i16;
    tmp19 = i14 + i17;
    tmp21 = i10 - i13;
    tmp23 = i11 - i14;
    tmp25 = i12 + i15;
    tmp27 = i12 - i9;

    i1  = tmp3  + i2;
    i2  = tmp3  - i2;
    i3  = tmp5  + tmp27;
    i4  = tmp5  - tmp27;
    i5  = tmp7  + tmp25;
    i6  = tmp7  - tmp25;
    i7  = tmp17 + tmp19;
    i8  = tmp17 - tmp19;
    i9  = tmp21 - tmp23;
    i10 = tmp21 + tmp23;
    i11 = i1    + tmp9;
    i12 = i2    + tmp11;
    i13 = i11   - tmp13;
    i14 = i12   - tmp15;
    i15 = i12   + tmp15;
    i16 = i11   + tmp13;

    *vec++ = r0;
    *vec++ = i0;
    *vec++ = r13 + r5 + r7;
    *vec++ = i13 + i5 + i7;
    *vec++ = r15 + r3 - r9;
    *vec++ = i15 + i3 - i9;
    *vec++ = r0  + r4;
    *vec++ = i0  + i4;
    *vec++ = r13 + r6 - r7;
    *vec++ = i13 + i6 - i7;
    *vec++ = r2;
    *vec++ = i2;
    *vec++ = r0  + r5;
    *vec++ = i0  + i5;
    *vec++ = r16 + r3 - r10;
    *vec++ = i16 + i3 - i10;
    *vec++ = r15 + r4 + r9;
    *vec++ = i15 + i4 + i9;
    *vec++ = r0  + r6;
    *vec++ = i0  + i6;
    *vec++ = r1;
    *vec++ = i1;
    *vec++ = r14 + r5 + r8;
    *vec++ = i14 + i5 + i8;
    *vec++ = r0  + r3;
    *vec++ = i0  + i3;
    *vec++ = r16 + r4 + r10;
    *vec++ = i16 + i4 + i10;
    *vec++ = r14 + r6 - r8;
    *vec++ = i14 + i6 - i8;


}

/*******************************************************************************
 Functionname:  fft5s
 *******************************************************************************

 Description:   5-point FFT.

 Arguments:     x      - pointer to input data (interleaved real / imaginary parts)
                stride - stride for input data

 Return:        none

*******************************************************************************/
static const float C51 =  0.9510565162951535f;
static const float C52 = -1.5388417685876270f;
static const float C53 = -0.3632712640026803f;
static const float C54 =  0.5590169943749475f;
static const float C55 = -1.25f;

static void fft5s(float *x, int stride)
{
    float r1,r2,r3,r4;
    float s1,s2,s3,s4;
    float t;
    /* real part */
    r1     = x[1*stride] + x[4*stride];
    r4     = x[1*stride] - x[4*stride];
    r3     = x[2*stride] + x[3*stride];
    r2     = x[2*stride] - x[3*stride];
    t      = (r1-r3) * C54;
    r1     = r1 + r3;
    x[0] = x[0] + r1;
    r1     = x[0] + (r1 * C55);
    r3     = r1 - t;
    r1     = r1 + t;
    t      = (r4 + r2) * C51;
    r4     = t + (r4 * C52);
    r2     = t + (r2 * C53);

    /* imaginary part */
    s1     = x[1*stride+1] + x[4*stride+1];
    s4     = x[1*stride+1] - x[4*stride+1];
    s3     = x[2*stride+1] + x[3*stride+1];
    s2     = x[2*stride+1] - x[3*stride+1];
    t      = (s1 - s3) * C54;
    s1     = s1 + s3;
    x[1] = x[1] + s1;
    s1     = x[1] + (s1 * C55);
    s3     = s1 - t;
    s1     = s1 + t;
    t      = (s4 + s2) * C51;
    s4     = t + (s4 * C52);
    s2     = t + (s2 * C53);

    /* combination */
    x[1*stride] = r1 + s2;
    x[4*stride] = r1 - s2;
    x[2*stride] = r3 - s4;
    x[3*stride] = r3 + s4;

    x[1*stride+1] = s1 - r2;
    x[4*stride+1] = s1 + r2;
    x[2*stride+1] = s3 + r4;
    x[3*stride+1] = s3 - r4;
}


/**
 * \brief    Function performs a complex 10-point FFT
 *           The FFT is performed inplace. The result of the FFT
 *           is scaled by SCALEFACTOR10 bits.
 *
 *           WOPS FLC version:                    1093 cycles
 *           WOPS with 32x16 bit multiplications:  196 cycles
 *
 * \param    [i/o] re    real input / output
 * \param    [i/o] im    imag input / output
 * \param    [i  ] s     stride real and imag input / output
 *
 * \return   void
 */
static void fft10(float * restrict vec)
{
    float t;
    float r1,r2,r3,r4;
    float s1,s2,s3,s4;
    float y00,y01,y02,y03,y04,y05,y06,y07,y08,y09;
    float y10,y11,y12,y13,y14,y15,y16,y17,y18,y19;

    /* 2 fft5 stages */

    /* real part */
    r1  = vec[12] + vec[8];
    r4  = vec[12] - vec[8];
    r3  = vec[4] + vec[16];
    r2  = vec[4] - vec[16];
    t   = (r1 - r3) * C54;
    r1  = r1 + r3;
    y00 = vec[0] + r1;
    r1  = y00 + (r1 * C55);
    r3  = r1 - t;
    r1  = r1 + t;
    t   = (r4 + r2) * C51;
    r4  = t + (r4 * C52);
    r2  = t + (r2 * C53);

    /* imaginary part */
    s1  = vec[13] + vec[9];
    s4  = vec[13] - vec[9];
    s3  = vec[5] + vec[17];
    s2  = vec[5] - vec[17];
    t   = (s1 - s3) * C54;
    s1  = s1 + s3;
    y01 = vec[1] + s1;
    s1  = y01 + (s1 * C55);
    s3  = s1 - t;
    s1  = s1 + t;
    t   = (s4 + s2) * C51;
    s4  = t + (s4 * C52);
    s2  = t + (s2 * C53);

    /* combination */
    y04 = r1 + s2;
    y16 = r1 - s2;
    y08 = r3 - s4;
    y12 = r3 + s4;
    y05 = s1 - r2;
    y17 = s1 + r2;
    y09 = s3 + r4;
    y13 = s3 - r4;

    /* real part */
    r1  = vec[2] + vec[18];
    r4  = vec[2] - vec[18];
    r3  = vec[14] + vec[6];
    r2  = vec[14] - vec[6];
    t   = (r1 - r3) * C54;
    r1  = r1 + r3;
    y02 = vec[10] + r1;
    r1  = y02 + (r1 * C55);
    r3  = r1 - t;
    r1  = r1 + t;
    t   = (r4 + r2) * C51;
    r4  = t + (r4 * C52);
    r2  = t + (r2 * C53);

    /* imaginary part */
    s1  = vec[3] + vec[19];
    s4  = vec[3] - vec[19];
    s3  = vec[15] + vec[7];
    s2  = vec[15] - vec[7];
    t   = (s1 - s3) * C54;
    s1  = s1 + s3;
    y03 = vec[11] + s1;
    s1  = y03 + (s1 * C55);
    s3  = s1 - t;
    s1  = s1 + t;
    t   = (s4 + s2) * C51;
    s4  = t + (s4 * C52);
    s2  = t + (s2 * C53);

    /* combination */
    y06 = r1 + s2;
    y18 = r1 - s2;
    y10 = r3 - s4;
    y14 = r3 + s4;
    y07 = s1 - r2;
    y19 = s1 + r2;
    y11 = s3 + r4;
    y15 = s3 - r4;

    /* 5 fft2 stages */
    vec[0]  = y00 + y02;
    vec[1]  = y01 + y03;
    vec[2]  = y12 - y14;
    vec[3]  = y13 - y15;
    vec[4]  = y04 + y06;
    vec[5]  = y05 + y07;
    vec[6]  = y16 - y18;
    vec[7]  = y17 - y19;
    vec[8]  = y08 + y10;
    vec[9]  = y09 + y11;
    vec[10] = y00 - y02;
    vec[11] = y01 - y03;
    vec[12] = y12 + y14;
    vec[13] = y13 + y15;
    vec[14] = y04 - y06;
    vec[15] = y05 - y07;
    vec[16] = y16 + y18;
    vec[17] = y17 + y19;
    vec[18] = y08 - y10;
    vec[19] = y09 - y11;
}

/**
 * \brief    Function performs a complex 20-point FFT
 *           The FFT is performed inplace. The result of the FFT
 *           is scaled by SCALEFACTOR20 bits.
 *
 *           WOPS FLC version:                    1509 cycles
 *           WOPS with 32x16 bit multiplications:  432 cycles
 *
 * \param    [i/o] re    real input / output
 * \param    [i/o] im    imag input / output
 * \param    [i  ] s     stride real and imag input / output
 *
 * \return   void
 */
static void fft20(float *signal)
{
    const int s = 2;
    float *re = signal, *im = signal+1;
    float r1,r2,r3,r4;
    float s1,s2,s3,s4;
    float x0,x1,x2,x3,x4;
    float t,t0,t1,t2,t3,t4,t5,t6,t7;
    float y00,y01,y02,y03,y04,y05,y06,y07,y08,y09;
    float y10,y11,y12,y13,y14,y15,y16,y17,y18,y19;
    float y20,y21,y22,y23,y24,y25,y26,y27,y28,y29;
    float y30,y31,y32,y33,y34,y35,y36,y37,y38,y39;

    /*  */

    /* 1. FFT5 stage */

    /* real part */
    x0  = re[s* 0];
    x1  = re[s*16];
    x2  = re[s*12];
    x3  = re[s* 8];
    x4  = re[s* 4];
    r1  = x1 + x4;
    r4  = x1 - x4;
    r3  = x2 + x3;
    r2  = x2 - x3;
    t   = (r1 - r3) * C54;
    r1  = r1 + r3;
    y00 = x0 + r1;
    r1  = y00 + (r1 * C55);
    r3  = r1 - t;
    r1  = r1 + t;
    t   = (r4 + r2) * C51;
    r4  = t  + (r4 * C52);
    r2  = t + (r2 * C53);

    /* imaginary part */
    x0  = im[s* 0];
    x1  = im[s*16];
    x2  = im[s*12];
    x3  = im[s* 8];
    x4  = im[s* 4];
    s1  = x1 + x4;
    s4  = x1 - x4;
    s3  = x2 + x3;
    s2  = x2 - x3;
    t   = (s1 - s3) * C54;
    s1  = (s1 + s3);
    y01 = (x0 + s1);
    s1  = y01 + (s1 * C55);
    s3  = (s1 - t);
    s1  = (s1 + t);
    t   = (s4 + s2) * C51;
    s4  = t + (s4 * C52);
    s2  = t + (s2 * C53);

    /* combination */
    y08 = (r1 + s2);
    y32 = (r1 - s2);
    y16 = (r3 - s4);
    y24 = (r3 + s4);

    y09 = (s1 - r2);
    y33 = (s1 + r2);
    y17 = (s3 + r4);
    y25 = (s3 - r4);

    /* 2. FFT5 stage */

    /* real part */
    x0  = re[s* 5];
    x1  = re[s* 1];
    x2  = re[s*17];
    x3  = re[s*13];
    x4  = re[s* 9];
    r1  = (x1 + x4);
    r4  = (x1 - x4);
    r3  = (x2 + x3);
    r2  = (x2 - x3);
    t   = (r1 - r3) * C54;
    r1  = (r1 + r3);
    y02 = (x0 + r1);
    r1  = y02 + (r1 * C55);
    r3  = (r1 - t);
    r1  = (r1 + t);
    t   = (r4 + r2) * C51;
    r4  = t + (r4 * C52);
    r2  = t + (r2 * C53);

    /* imaginary part */
    x0  = im[s* 5];
    x1  = im[s* 1];
    x2  = im[s*17];
    x3  = im[s*13];
    x4  = im[s* 9];
    s1  = (x1 + x4);
    s4  = (x1 - x4);
    s3  = (x2 + x3);
    s2  = (x2 - x3);
    t   = (s1 - s3) * C54;
    s1  = (s1 + s3);
    y03 = (x0 + s1);
    s1  = y03 + (s1 * C55);
    s3  = (s1 - t);
    s1  = (s1 + t);
    t   = (s4 + s2) * C51;
    s4  = t + (s4 * C52);
    s2  = t + (s2 * C53);

    /* combination */
    y10 = (r1 + s2);
    y34 = (r1 - s2);
    y18 = (r3 - s4);
    y26 = (r3 + s4);

    y11 = (s1 - r2);
    y35 = (s1 + r2);
    y19 = (s3 + r4);
    y27 = (s3 - r4);

    /* 3. FFT5 stage */

    /* real part */
    x0  = re[s*10];
    x1  = re[s* 6];
    x2  = re[s* 2];
    x3  = re[s*18];
    x4  = re[s*14];
    r1  = (x1 + x4);
    r4  = (x1 - x4);
    r3  = (x2 + x3);
    r2  = (x2 - x3);
    t   = (r1 - r3) * C54;
    r1  = (r1 + r3);
    y04 = (x0 + r1);
    r1  = y04 + (r1 * C55);
    r3  = (r1 - t);
    r1  = (r1 + t);
    t   = (r4 + r2) * C51;
    r4  = t + (r4 * C52);
    r2  = t + (r2 * C53);

    /* imaginary part */
    x0  = im[s*10];
    x1  = im[s* 6];
    x2  = im[s* 2];
    x3  = im[s*18];
    x4  = im[s*14];
    s1  = (x1 + x4);
    s4  = (x1 - x4);
    s3  = (x2 + x3);
    s2  = (x2 - x3);
    t   = (s1 - s3) * C54;
    s1  = (s1 + s3);
    y05 = (x0 + s1);
    s1  = y05 + (s1 * C55);
    s3  = (s1 - t);
    s1  = (s1 + t);
    t   = (s4 + s2) * C51;
    s4  = t + (s4 * C52);
    s2  = t + (s2 * C53);

    /* combination */
    y12 = (r1 + s2);
    y36 = (r1 - s2);
    y20 = (r3 - s4);
    y28 = (r3 + s4);

    y13 = (s1 - r2);
    y37 = (s1 + r2);
    y21 = (s3 + r4);
    y29 = (s3 - r4);

    /* 4. FFT5 stage */

    /* real part */
    x0  = re[s*15];
    x1  = re[s*11];
    x2  = re[s* 7];
    x3  = re[s* 3];
    x4  = re[s*19];
    r1  = (x1 + x4);
    r4  = (x1 - x4);
    r3  = (x2 + x3);
    r2  = (x2 - x3);
    t   = (r1 - r3) * C54;
    r1  = (r1 + r3);
    y06 = (x0 + r1);
    r1  = y06 + (r1 * C55);
    r3  = (r1 - t);
    r1  = (r1 + t);
    t   = (r4 + r2) * C51;
    r4  = t + (r4 * C52);
    r2  = t + (r2 * C53);

    /* imaginary part */
    x0  = im[s*15];
    x1  = im[s*11];
    x2  = im[s* 7];
    x3  = im[s* 3];
    x4  = im[s*19];
    s1  = (x1 + x4);
    s4  = (x1 - x4);
    s3  = (x2 + x3);
    s2  = (x2 - x3);
    t   = (s1 - s3) * C54;
    s1  = (s1 + s3);
    y07 = (x0 + s1);
    s1  = y07 + (s1 * C55);
    s3  = (s1 - t);
    s1  = (s1 + t);
    t   = (s4 + s2) * C51;
    s4  = t + (s4 * C52);
    s2  = t + (s2 * C53);

    /* combination */
    y14 = (r1 + s2);
    y38 = (r1 - s2);
    y22 = (r3 - s4);
    y30 = (r3 + s4);

    y15 = (s1 - r2);
    y39 = (s1 + r2);
    y23 = (s3 + r4);
    y31 = (s3 - r4);


    /* 1. FFT4 stage */

    /* Pre-additions */
    t0 = (y00 + y04);
    t2 = (y00 - y04);
    t1 = (y01 + y05);
    t3 = (y01 - y05);
    t4 = (y02 + y06);
    t7 = (y02 - y06);
    t5 = (y07 + y03);
    t6 = (y07 - y03);

    /* Post-additions */
    re[s* 0] = (t0 + t4);
    im[s* 0] = (t1 + t5);
    re[s* 5] = (t2 - t6);
    im[s* 5] = (t3 - t7);
    re[s*10] = (t0 - t4);
    im[s*10] = (t1 - t5);
    re[s*15] = (t2 + t6);
    im[s*15] = (t3 + t7);

    /* 2. FFT4 stage */

    /* Pre-additions */
    t0 = (y08 + y12);
    t2 = (y08 - y12);
    t1 = (y09 + y13);
    t3 = (y09 - y13);
    t4 = (y10 + y14);
    t7 = (y10 - y14);
    t5 = (y15 + y11);
    t6 = (y15 - y11);

    /* Post-additions */
    re[s* 4] = (t0 + t4);
    im[s* 4] = (t1 + t5);
    re[s* 9] = (t2 - t6);
    im[s* 9] = (t3 - t7);
    re[s*14] = (t0 - t4);
    im[s*14] = (t1 - t5);
    re[s*19] = (t2 + t6);
    im[s*19] = (t3 + t7);


    /* 3. FFT4 stage */

    /* Pre-additions */
    t0 = (y16 + y20);
    t2 = (y16 - y20);
    t1 = (y17 + y21);
    t3 = (y17 - y21);
    t4 = (y18 + y22);
    t7 = (y18 - y22);
    t5 = (y23 + y19);
    t6 = (y23 - y19);

    /* Post-additions */
    re[s* 8] = (t0 + t4);
    im[s* 8] = (t1 + t5);
    re[s*13] = (t2 - t6);
    im[s*13] = (t3 - t7);
    re[s*18] = (t0 - t4);
    im[s*18] = (t1 - t5);
    re[s* 3] = (t2 + t6);
    im[s* 3] = (t3 + t7);

    /* 4. FFT4 stage */

    /* Pre-additions */
    t0 = (y24 + y28);
    t2 = (y24 - y28);
    t1 = (y25 + y29);
    t3 = (y25 - y29);
    t4 = (y26 + y30);
    t7 = (y26 - y30);
    t5 = (y31 + y27);
    t6 = (y31 - y27);

    /* Post-additions */
    re[s*12] = (t0 + t4);
    im[s*12] = (t1 + t5);
    re[s*17] = (t2 - t6);
    im[s*17] = (t3 - t7);
    re[s* 2] = (t0 - t4);
    im[s* 2] = (t1 - t5);
    re[s* 7] = (t2 + t6);
    im[s* 7] = (t3 + t7);

    /* 5. FFT4 stage */

    /* Pre-additions */
    t0 = (y32 + y36);
    t2 = (y32 - y36);
    t1 = (y33 + y37);
    t3 = (y33 - y37);
    t4 = (y34 + y38);
    t7 = (y34 - y38);
    t5 = (y39 + y35);
    t6 = (y39 - y35);

    /* Post-additions */
    re[s*16] = (t0 + t4);
    im[s*16] = (t1 + t5);
    re[s* 1] = (t2 - t6);
    im[s* 1] = (t3 - t7);
    re[s* 6] = (t0 - t4);
    im[s* 6] = (t1 - t5);
    re[s*11] = (t2 + t6);
    im[s*11] = (t3 + t7);

    /*  */
    /*  */
}

/*******************************************************************************
 Functionname:  fft30
 *******************************************************************************

 Description:   30-point FFT.

 Arguments:     in - pointer to data (interleaved real / imaginary parts)

 Return:        none

*******************************************************************************/

static void fft30(float * restrict in)
{
    int i;
    float temp[60];
    float * temp_l  = temp;
    float * temp_lu = temp + 2*8;
    float * temp_h  = temp + 2*15;
    float * temp_hu = temp + 2*15 + 2*8;
    float *in_l     = in + 2*0;
    float *in_h     = in + 2*15;
    for(i=0; i<7; i++)
    {
        *temp_l++ = *in_l++;
        *temp_l++ = *in_l++;
        *temp_h++ = *in_h++;
        *temp_h++ = *in_h++;
        *temp_l++ = *in_h++;
        *temp_l++ = *in_h++;
        *temp_h++ = *in_l++;
        *temp_h++ = *in_l++;
    }
    *temp_l++ = *in_l++;
    *temp_l++ = *in_l++;
    *temp_h++ = *in_h++;
    *temp_h++ = *in_h++;
    temp_l = temp;
    temp_h = temp + 30;
    fft15(temp_l);
    fft15(temp_h);

    in_l = in + 2*0;
    in_h = in + 2*15;
    for(i=0; i<7; i++)
    {
        *in_l++ = *temp_l   + *temp_h;
        *in_h++ = *temp_l++ - *temp_h++;
        *in_l++ = *temp_l   + *temp_h;
        *in_h++ = *temp_l++ - *temp_h++;

        *in_h++ = *temp_lu   + *temp_hu;
        *in_l++ = *temp_lu++ - *temp_hu++;
        *in_h++ = *temp_lu   + *temp_hu;
        *in_l++ = *temp_lu++ - *temp_hu++;
    }
    *in_l++ = *temp_l   + *temp_h;
    *in_h++ = *temp_l++ - *temp_h++;
    *in_l++ = *temp_l   + *temp_h;
    *in_h++ = *temp_l++ - *temp_h++;
}

/*-------------------------------------------------------------------*
 * fft_cldfb()
 *
 * Interface functions FFT subroutines
 *--------------------------------------------------------------------*/
void fft_cldfb (
    float *data,       /* i/o: input/output vector */
    int size           /* size of fft operation */
)
{

    switch(size)
    {
    case 5:
        fft5s(data,2);
        break;
    case 8:
        fft8(data);
        break;
    case 10:
        fft10(data);
        break;
    case 16:
        fft16(data);
        break;
    case 20:
        fft20(data);
        break;
    case 30:
        fft30(data);
        break;

    default:
        assert(0);
        break;
    }
}

