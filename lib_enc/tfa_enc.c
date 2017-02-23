/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <assert.h>
#include "prot.h"



void tfaCalcEnv(const float* shb_speech, float* enr)
{
    int i, j, k;
    for (i=0, k=0; i<N_TEC_TFA_SUBFR; i++)
    {
        enr[i] = 1e-12f;
        for(j=0; j<L_TEC_TFA_SUBFR16k; j++) /* XX/2 since Fs = 16kHz */
        {
            enr[i] += shb_speech[k] * shb_speech[k];
            k++;
        }
    }
}

short tfaEnc_TBE(float* enr, short last_core, float* voicing, float* pitch_buf)
{
    int i;
    float m_g, m_a; /* m_g: geometrical mean, m_a: arithmetical mean */
    float flatness = 0;
    float voicing_sum;
    short tfa_flag;

    float pitch_buf_sum;
    const float m_a_bottom = 10000.0 / N_TEC_TFA_SUBFR;

    tfa_flag = 0;

    m_a = 1e-12f;
    m_g = 0.0;
    for (i=0; i<N_TEC_TFA_SUBFR; i++)
    {
        m_a += enr[i];
        m_g += log10(enr[i]);
    }
    m_a /= N_TEC_TFA_SUBFR;
    m_g /= N_TEC_TFA_SUBFR;
    m_g = pow(10.0, m_g);

    flatness = m_g / m_a;

    voicing_sum = voicing[0] + voicing[1];

    pitch_buf_sum = pitch_buf[0] + pitch_buf[1] + pitch_buf[2] + pitch_buf[3];

    if ((flatness > 0.70 && pitch_buf_sum > 4.0 * 110 && voicing_sum > 2.0 * 0.70) ||
            (last_core == TCX_20_CORE && flatness > 0.50 && voicing_sum < 2.0 * 0.70))
    {
        tfa_flag = 1;
    }
    /* energy lower limit */
    if(m_a < m_a_bottom)
    {
        tfa_flag = 0;
    }

    return tfa_flag;
}
