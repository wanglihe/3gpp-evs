/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include <math.h>
#include "options.h"
#include "cnst.h"
#include "prot.h"
#include "rom_com.h"


/*-------------------------------------------------------------------*
 * hf_parinitiz()
 *
 *
 *-------------------------------------------------------------------*/

void hf_parinitiz(
    const long  total_brate,
    const short hqswb_clas,
    short lowlength,
    short highlength,
    short wBands[],
    const short **subband_search_offset,
    const short **subband_offsets,
    short *nBands,
    short *nBands_search,
    short *swb_lowband,
    short *swb_highband
)
{
    *swb_lowband = lowlength;
    *swb_highband = highlength;

    if( hqswb_clas == HQ_HARMONIC )
    {
        /* Mode dependent initializations (performed every frame in case mode-switching implemented) */
        *nBands = NB_SWB_SUBBANDS_HAR;
        *nBands_search = NB_SWB_SUBBANDS_HAR_SEARCH_SB;

        if ( total_brate == HQ_13k20 )
        {
            wBands[0] = SWB_SB_BW_LEN0_12KBPS_HAR;
            wBands[1] = SWB_SB_BW_LEN1_12KBPS_HAR;
            wBands[2] = SWB_SB_BW_LEN2_12KBPS_HAR;
            wBands[3] = SWB_SB_BW_LEN3_12KBPS_HAR;
            *subband_offsets = subband_offsets_sub5_13p2kbps_Har;
            *subband_search_offset = subband_search_offsets_13p2kbps_Har;
        }
        else
        {
            wBands[0] = SWB_SB_BW_LEN0_16KBPS_HAR;
            wBands[1] = SWB_SB_BW_LEN1_16KBPS_HAR;
            wBands[2] = SWB_SB_BW_LEN2_16KBPS_HAR;
            wBands[3] = SWB_SB_BW_LEN3_16KBPS_HAR;
            *subband_offsets = subband_offsets_sub5_16p4kbps_Har;
            *subband_search_offset = subband_search_offsets_16p4kbps_Har;
        }
    }
    else
    {
        /* Mode-dependent initializations (performed every frame in case mode-switching implemented) */
        *nBands = NB_SWB_SUBBANDS;
        *nBands_search = NB_SWB_SUBBANDS;

        if ( total_brate == HQ_13k20 )
        {
            wBands[0] = SWB_SB_LEN0_12KBPS;
            wBands[1] = SWB_SB_LEN1_12KBPS;
            wBands[2] = SWB_SB_LEN2_12KBPS;
            wBands[3] = SWB_SB_LEN3_12KBPS;
            *subband_offsets = subband_offsets_12KBPS;
        }
        else
        {
            wBands[0] = SWB_SB_LEN0_16KBPS;
            wBands[1] = SWB_SB_LEN1_16KBPS;
            wBands[2] = SWB_SB_LEN2_16KBPS;
            wBands[3] = SWB_SB_LEN3_16KBPS;
            *subband_offsets = subband_offsets_16KBPS;
        }
    }
    return;
}

/*-------------------------------------------------------------------*
 * GetPredictedSignal()
 *
 * Routine for calculating the predicted signal
 *-------------------------------------------------------------------*/
void GetPredictedSignal(
    const float *predBuf,   /* i  : prediction buffer        */
    float *outBuf,    /* o  : output buffer            */
    const short lag,        /* i  : prediction buffer offset */
    const short fLen,       /* i  : length of loop (output)  */
    const float gain        /* i  : gain to be applied       */
)
{
    short i;

    predBuf += lag;

    for(i = 0; i < fLen; i++)
    {
        *outBuf++ = *predBuf++ * gain;
    }

    return;
}
/*-------------------------------------------------------------------*
 * est_freq_har_decis()
 *
 * Harmonic frequency decision matrix
 *-------------------------------------------------------------------*/
static
void est_freq_har_decis(
    short *har_freq_est1,       /* o: harmonic analysis 1*/
    short *har_freq_est2,       /* o: harmonic analysis 2*/
    short sharp,                /* i: pka-avg for group 1 */
    short sharp1,               /* i: pka-avg for group 2 */
    short hfe_est_countk1,      /* i: group pks count 1*/
    short hfe_est_countk2,      /* i: group pks count 2*/
    short k,                    /* i: group count */
    short k1,
    short k2,
    short *prev_frm_hfe2        /* i: harmonic estimation  */
)
{
    short temp_hfe2 = 0;

    if( k != 0)
    {
        *har_freq_est1 = (int) sharp /k;
    }
    if(k1 > 1)
    {
        *har_freq_est2 = (int) sharp1 /k1;
    }
    else if( (k1 < 2 && (k2 != 0 || k>1)))
    {
        *har_freq_est2 = (int) *har_freq_est1;
    }
    else
    {
        if((hfe_est_countk1 != 0 || hfe_est_countk2 != 0) && (k1 == 0 && k2 == 0))
        {
            *har_freq_est2 = (*har_freq_est1);
        }
        else
        {
            *har_freq_est2 = 2*(*har_freq_est1);
        }
    }

    /* Consider Estimation Error upto 200Hz */
    if( *prev_frm_hfe2 != 0 && ( (short)abs(*prev_frm_hfe2 - *har_freq_est2) < 10 || abs(*prev_frm_hfe2-2*(*har_freq_est2))<10))
    {
        *har_freq_est2 = *prev_frm_hfe2;
    }
    else if(*prev_frm_hfe2 != 0 && abs(*har_freq_est2-2*(*prev_frm_hfe2))<10)
    {
        *har_freq_est2 =2*( *prev_frm_hfe2);
    }
    else
    {
        temp_hfe2 = (short)(*prev_frm_hfe2+*har_freq_est2)/2;

        if(abs(temp_hfe2 -*prev_frm_hfe2) < 2)
        {
            temp_hfe2 =*prev_frm_hfe2;
            *har_freq_est2 = temp_hfe2;
        }
    }
    if( *har_freq_est2 <*har_freq_est1 && (k>1 && k1<2))
    {
        *har_freq_est2 = *har_freq_est1;
    }
    return;
}

/*--------------------------------------------------------------------------*
 * har_est()
 *
 * Harmonic Structure analysis using LF spectrum
 *--------------------------------------------------------------------------*/

short har_est(
    float spectra[],                    /* i  : coded spectrum                 */
    short N,                            /* i  : length of the desired spectrum */
    short *har_freq_est1,               /* i/o: Estimation harmonics 1         */
    short *har_freq_est2,               /* o  : Estimation harmonics 2         */
    short *flag_dis,                    /* i/o: flag for BWE reconstruction    */
    short *prev_frm_hfe2,               /* i/o: Estimated harmonic update      */
    const short subband_search_offset[],/* i  : Subband Search range           */
    const short sbWidth[],              /* i  : Subband Search range           */
    short *prev_stab_hfe2               /* i/o: Estimated harmonic position    */
)
{
    float peak;
    short i, j, q, k , k1, k2;
    float input_abs[L_FRAME32k],blk_peak[30];
    short blk_end,blk_st;
    short peak_pos,blk_peak_pos[30],diff_peak_pos[30], sharp, sharp1, sharp2;
    int min_har_pos;
    short blk_peak_pos_te[30];
    float blk_peak_te[30];
    short temp;
    short hfe_est_countk,hfe_est_countk1,hfe_est_countk2;
    short r1, r2, r3;
    short start_pos;
    float blk_peak_max;
    short blk_peak_pos_max;

    short nlags, ct_hfsb2,sum_diff=0;
    short blk_peak_pos_hfsb2[30],diff_peak_pos_hfsb2[30];
    short rem_hfe2 ,q_diffpos_hfe2 = 0,diff_posmax_hfe2, q_diffpos_prevhfe2;
    short blk_peak_max_idx,blk_peak_pos_max_diff,diff_peak_pos_te[30];
    rem_hfe2 = 0;
    q_diffpos_hfe2 = 0;
    diff_posmax_hfe2 = 0;
    q_diffpos_prevhfe2 = 0;

    set_s(blk_peak_pos_te,0,30);
    set_s(diff_peak_pos,0,30);

    r1 = SWB_HAR_RAN1;
    r2 = SWB_HAR_RAN2;
    r3 = SWB_HAR_RAN3;
    start_pos =r1;

    /* Copy the abs values of LF spectrum*/
    for ( i = start_pos; i < N; i++)
    {
        input_abs[i] = (float) fabs(spectra[i]);
    }

    blk_end = (short) N/LR_BLK_LEN;
    blk_st = (short) start_pos/LR_BLK_LEN;

    if( (float)N/(LR_BLK_LEN) - blk_end > 0.0f)
    {
        blk_end++;
    }

    /* initialization of over buffer for fractional point */
    for(i=N; i<(blk_end*LR_BLK_LEN); i++)
    {
        input_abs[i] = 0.0f;
    }

    q = start_pos;

    /* Block Processing, to detect the spectral peaks*/
    for(i = blk_st; i < blk_end; i++)
    {
        peak = 0.0f;
        peak_pos = 0;

        for(j = 0; j < LR_BLK_LEN; j ++, q ++)
        {
            if (input_abs[q] > peak)
            {
                peak = input_abs[q];
                peak_pos = q;
            }

            if( i > blk_st && input_abs[q] != 0 && input_abs[q] == peak && (peak_pos - blk_peak_pos[i-1]) < LR_HLF_PK_BLK_LEN )
            {
                peak = input_abs[q];
                peak_pos = q;
            }
        }

        blk_peak[i] = peak;
        blk_peak_pos[i] = peak_pos;
    }

    for(i = blk_st; i < blk_end; i++)
    {
        if(i > blk_st)
        {
            if(blk_peak_pos[i] != 0 && blk_peak_pos[i-1] != 0)
            {
                if((blk_peak_pos[i] - blk_peak_pos[i-1]) < LR_LOWBAND_DIF_PK_LEN)
                {
                    if(blk_peak[i] > blk_peak[i-1])
                    {
                        blk_peak[i-1] = 0.0f;
                        blk_peak_pos[i-1] = 0;
                    }
                    else
                    {
                        blk_peak[i] = blk_peak[i-1];
                        blk_peak_pos[i] = blk_peak_pos[i-1];
                        blk_peak[i-1] = 0.0f;
                        blk_peak_pos[i-1] = 0;
                    }
                }
            }
        }
    }

    /* peak counts in each group */
    j = 0;
    hfe_est_countk = 0;
    hfe_est_countk1 = 0;
    hfe_est_countk2 = 0;
    for(i = blk_st; i < blk_end; i++)
    {
        if(blk_peak_pos[i] != 0 )
        {
            blk_peak_pos_te[j] = blk_peak_pos[i];
            if( blk_peak_pos[i] <r2 )
            {
                hfe_est_countk++;
            }
            else if( blk_peak_pos[i] <r3 )
            {
                hfe_est_countk1++;
            }
            else
            {
                hfe_est_countk2++;
            }

            blk_peak_te[j] = blk_peak[i];
            j++;
        }
    }

    min_har_pos = SWB_HAR_RAN1;
    temp = 0;
    blk_peak_max = blk_peak_te[0];
    blk_peak_pos_max = blk_peak_pos_te[0] ;
    blk_peak_max_idx = 0;
    for(i = 1; i < j; i++)
    {
        diff_peak_pos[i-1] = blk_peak_pos_te[i] - blk_peak_pos_te[i-1];
        if( diff_peak_pos[i-1] <= min_har_pos )
        {
            min_har_pos = diff_peak_pos[i-1];
        }

        if(blk_peak_te[i-1] > blk_peak_max)
        {
            blk_peak_max = blk_peak_te[i-1];
            blk_peak_pos_max = blk_peak_pos_te[i-1] ;
            blk_peak_max_idx = i-1;
        }

        temp++;
    }
    blk_peak_pos_max_diff = diff_peak_pos[blk_peak_max_idx];

    /* Decision for BWE reconstruction */
    if((hfe_est_countk < 2 && hfe_est_countk1 <2 && hfe_est_countk2 <2 ) || min_har_pos >= SWB_HAR_RAN1)
    {
        *flag_dis = 0;
        if((hfe_est_countk == 1 && hfe_est_countk1 == 1)&& (hfe_est_countk2 == 1 || hfe_est_countk2 == 0))
        {
            *flag_dis = 1;
        }
    }
    for(i=0; i<temp; i++)
    {
        if(blk_peak_pos_max_diff+LR_LOWBAND_DIF_PK_LEN < diff_peak_pos[i])
        {
            diff_peak_pos[i] = 0;
        }
    }
    mvs2s(diff_peak_pos,diff_peak_pos_te,temp);
    set_s(diff_peak_pos,-1,temp);
    j=0;
    for(i=0; i<temp; i++)
    {
        if(diff_peak_pos_te[i] != 0)
        {
            diff_peak_pos[j] = diff_peak_pos_te[i];
            j++;
        }
    }
    temp = j;
    /* harmonic estimation analysis to perform BWE Reconstruction */
    if( *flag_dis )
    {
        sharp = 0;
        k = 0;
        k1 = 0;
        sharp1 = 0;
        sharp2 = 0;
        k2 = 0;

        for(i=0, q=1; i<temp; i++, q++)
        {
            if((diff_peak_pos[i] <= (min_har_pos+LR_LOWBAND_DIF_PK_LEN)) && diff_peak_pos[i] > 0)
            {
                sharp += diff_peak_pos[i];
                k++;
            }
            else if((diff_peak_pos[i] <= (min_har_pos+2*LR_LOWBAND_DIF_PK_LEN)) && diff_peak_pos[i] > 0)
            {
                sharp1 += diff_peak_pos[i];
                k1++;
            }
            else if(diff_peak_pos[i] > 0)
            {
                sharp2 += diff_peak_pos[i];
                k2++;
            }
        }

        est_freq_har_decis(har_freq_est1,har_freq_est2,sharp,sharp1,hfe_est_countk1,hfe_est_countk2,k,k1,k2,prev_frm_hfe2);

        blk_peak_pos_max = blk_peak_pos_te[temp-1];

        if((*prev_stab_hfe2) > 0 && (*prev_frm_hfe2) > 0 && *prev_stab_hfe2 < N)
        {
            rem_hfe2 = *har_freq_est2%(*prev_frm_hfe2);
            diff_posmax_hfe2 = (short) abs(blk_peak_pos_max-*prev_stab_hfe2);
            if(rem_hfe2 == 0)
            {
                if(diff_posmax_hfe2 < 9 || *har_freq_est2 == 0)
                {
                    blk_peak_pos_max = *prev_stab_hfe2;
                }
                else
                {
                    q_diffpos_hfe2 = diff_posmax_hfe2/(*har_freq_est2);
                    q_diffpos_prevhfe2 = diff_posmax_hfe2/(*prev_frm_hfe2);
                    if(q_diffpos_hfe2 < 10 || q_diffpos_prevhfe2 < 10)
                    {
                        blk_peak_pos_max = *prev_stab_hfe2;
                    }
                    else
                    {
                        *prev_stab_hfe2 = blk_peak_pos_max;
                    }
                }
            }
            else
            {
                *prev_stab_hfe2 = blk_peak_pos_max;
            }
        }
        else
        {
            *prev_stab_hfe2 = blk_peak_pos_max;
        }

        if(*har_freq_est1==0  ||  *har_freq_est2 ==0)
        {
            *flag_dis =0;
        }
    }

    if(*flag_dis == 0)
    {
        if(*prev_frm_hfe2 != 0)
        {
            *har_freq_est2 = *prev_frm_hfe2;
        }
        else
        {
            nlags = (short)pow(2, bits_lagIndices_mode0_Har[0]);
            ct_hfsb2 = 0;
            for(i = 0; i < j; i++)
            {
                if(blk_peak_pos_te[i] >=(subband_search_offset[0]- nlags/2)
                        && blk_peak_pos_te[i] <(subband_search_offset[0] + sbWidth[0]+ nlags/2))
                {
                    blk_peak_pos_hfsb2[ct_hfsb2] = blk_peak_pos_te[i];
                    ct_hfsb2++;
                }
            }

            if(ct_hfsb2 >1)
            {
                for(i=1; i<ct_hfsb2; i++)
                {
                    diff_peak_pos_hfsb2[i-1] = blk_peak_pos_hfsb2[i] - blk_peak_pos_hfsb2[i-1];
                    sum_diff +=diff_peak_pos_hfsb2[i-1];
                }
                *har_freq_est2 = sum_diff/ct_hfsb2;
            }
            else
            {
                *har_freq_est2 = (short) min_har_pos;
            }
        }
    }


    return blk_peak_pos_max;
}

void genhf_noise(
    float noise_flr[],                /* i  : smoothed non tonal                           */
    float xSynth_har[],               /* o  : hf non tonal components                      */
    float *predBuf,                   /* i  : smoothed tonal compone                       */
    short bands,                      /* i  : total number of subbands in a frame          */
    short harmonic_band,              /* i  : Number of LF harmonic frames                 */
    short har_freq_est2,              /* i  : harmonic signal parameter                    */
    short pos_max_hfe2,               /* i  : last pulse in core coder                     */
    short *pul_res,                   /* o  : pulse resolution                             */
    GainItem pk_sf[],                 /* o  : representative region                        */
    const short fLenLow,                    /* i  : low frequency length                         */
    const short fLenHigh,                   /* i  : high frequency length                        */
    const short sbWidth[],                  /* i  : bandwidth for high bands                     */
    const short lagIndices[],               /* i  : correlation indices for most representative  */
    const short subband_offsets[],          /* i  : band offsets for HF reconstruction           */
    const short subband_search_offset[]     /* i  : most representative regions offsets in LF    */
)
{
    short k,j,ii,st_pos,dst_pos;
    short nlags[NB_SWB_SUBBANDS_HAR_SEARCH_SB];
    float tmpbuf[L_FRAME32k];
    short hfband_end[NB_SWB_SUBBANDS];
    short rem_hfe, temp_last_peakpos,i,l,pos,res;
    float hf_pulse_peaks[160],pulse_peak_sb[320];
    short st_last_peakpos;

    set_f(tmpbuf,0.0f,L_FRAME32k);
    for(k=0; k<3; k++)
    {
        hfband_end[k] = fLenLow+subband_offsets[k+1];
    }
    hfband_end[3] = fLenLow+fLenHigh;

    rem_hfe = (short)(fLenLow-pos_max_hfe2-1)/har_freq_est2;
    st_last_peakpos = pos_max_hfe2+(rem_hfe*har_freq_est2);
    temp_last_peakpos = st_last_peakpos;
    i=0;

    for(k=0; k<2; k++)
    {
        nlags[k] = (short)pow(2, bits_lagIndices_mode0_Har[k]);

        l=0;
        while(st_last_peakpos < (fLenLow+subband_offsets[k]))
        {
            st_last_peakpos +=har_freq_est2;
        }
        st_last_peakpos-=har_freq_est2;

        if(k==0)
        {
            st_pos=subband_search_offset[k]-nlags[k]/2+lagIndices[k];

            /*Copy the LF Smoothed Noise to the HF*/
            for(j=0; j<sbWidth[k]; j++)
            {
                xSynth_har[j]= noise_flr[st_pos+j];
                tmpbuf[j] = xSynth_har[j];
                if(predBuf[st_pos+j] != 0.0f)
                {
                    hf_pulse_peaks[l] = predBuf[st_pos+j];
                    l++;
                }
            }
        }
        else
        {
            st_pos=subband_search_offset[k]+nlags[k]/2 -lagIndices[k];
            dst_pos = st_pos-sbWidth[k];
            ii = sbWidth[k-1];
            /*Copy the LF Smoothed Noise floor to the HF*/
            for(j= st_pos; j>(dst_pos) && ii<(sbWidth[k]+sbWidth[k-1]); j--)
            {
                xSynth_har[ii] = noise_flr[j];
                tmpbuf[ii] = xSynth_har[ii];
                if(predBuf[j] != 0.0f)
                {
                    hf_pulse_peaks[l] = predBuf[j];
                    l++;
                }
                ii++;
            }
        }
        pos =0;
        for(j = 0; j< l; j++)
        {
            st_last_peakpos+=har_freq_est2;
            if(st_last_peakpos< hfband_end[k])
            {
                pk_sf[k*8+pos].nmrValue = hf_pulse_peaks[j];
                pk_sf[k*8+pos].gainIndex = st_last_peakpos-fLenLow;
                pul_res[k]++;
                pulse_peak_sb[i] = hf_pulse_peaks[j];
                i++;
                pos++;
            }
        }
        st_last_peakpos = temp_last_peakpos;
    }
    res=i-1;
    l=1;
    ii = hfband_end[k-1]-fLenLow-1;
    for(; k<(bands-harmonic_band); k++)
    {
        for(j=hfband_end[k-1]-fLenLow; j<(hfband_end[k]-fLenLow); j++)
        {
            xSynth_har[j] = tmpbuf[ii];
            tmpbuf[j] = xSynth_har[j];
            ii--;
        }
        pos=0;
        while(st_last_peakpos < hfband_end[k-1])
        {
            st_last_peakpos +=har_freq_est2;
        }
        while(st_last_peakpos < hfband_end[k] && pul_res[k]<pul_res[2-l] && l<=2)
        {
            pk_sf[k*8+pos].nmrValue = pulse_peak_sb[res];
            pk_sf[k*8+pos].gainIndex = st_last_peakpos-fLenLow;
            pul_res[k]++;
            res--;
            pos++;
            st_last_peakpos+=har_freq_est2;
        }
        l++;
    }

    return;
}

/*-------------------------------------------------------------------*
 * SmoothSpec()
 *
 * Smoothes specified samples using moving average method. The number
 * of points in the average is given by 'span'. Note that current
 * implementation does not accept 'span' to be smaller than 'fLen'.
 *-------------------------------------------------------------------*/
static
void SmoothSpec(
    const float *inBuf,         /* i  : input               */
    float       *outBuf,        /* o  : output              */
    const short fLen,           /* i  : length              */
    short       span            /* i  : averaging length    */
)
{
    short i, span1, nItems;
    float sum, ispan;
    const float *oldPtr, *newPtr;

    /* not accepted */
    if( span > fLen )
    {
        mvr2r( inBuf, outBuf, fLen );
        return;
    }

    /* span must be odd */
    if( (span & 0x1) == 0 )
    {
        span--;
    }

    span1 = span >> 1;

    /* first sample */
    sum = *inBuf;
    *outBuf++ = sum;

    oldPtr = inBuf;
    newPtr = inBuf + 2;

    /* handle start */
    inBuf++;
    sum += *inBuf;

    for( i = 1, nItems = 3; i < span1; i++, inBuf++ )
    {
        sum += *newPtr++;
        *outBuf++ = sum / nItems;

        sum += *newPtr++;
        nItems += 2;
    }

    ispan = 1.0f / span;
    inBuf++;
    i++;
    sum += *newPtr++;
    *outBuf++ = sum * ispan;

    /* moving average */
    for( ; i < fLen - span1; i++, inBuf++ )
    {
        sum += *newPtr++;
        sum -= *oldPtr++;
        *outBuf++ = sum * ispan;
    }

    /* handle end */
    nItems = span - 2;
    sum -= *oldPtr++;

    for( ; i < fLen - 1; i++, inBuf++ )
    {
        sum -= *oldPtr++;
        *outBuf++ = sum / nItems;
        nItems -= 2;
        sum -= *oldPtr++;
    }

    /* last sample */
    *outBuf = *inBuf;

    return;
}

/*-------------------------------------------------------------------*
 * SpectrumSmoothing()
 *
 * Smoothing of the low-frequency envelope
 *-------------------------------------------------------------------*/

void SpectrumSmoothing(
    float *inBuf,               /* i  : input     */
    float *outBuf,              /* o  : output    */
    const short fLen,           /* i  : length    */
    const float th_cut          /* i  : threshold of cut */
)
{
    short i, k;
    float inBuf_pss[L_FRAME32k];
    float outBuf_pss[L_FRAME32k];
    float max_val[L_FRAME32k/L_SB];
    float max_val_norm;
    float inBuf_abs;
    short j;
    short num_subband_smooth;
    short m, n;
    short cnt_zero_cont;
    short n_list[BANDS_MAX];
    short reset_flag;

    num_subband_smooth = (short) fLen/L_SB;
    if( (float)fLen / L_SB - num_subband_smooth > 0.0f )
    {
        num_subband_smooth++;
    }

    for( i=0; i<fLen; i++ )
    {
        inBuf_pss[i] = inBuf[i];
        outBuf_pss[i] = 0.0f;
    }

    for( i=fLen; i<fLen + (num_subband_smooth * L_SB - fLen); i++ )
    {
        inBuf_pss[i] = 0.0f;
        outBuf_pss[i] = 0.0f;
    }

    j = 0;
    for ( i=0; i<num_subband_smooth; i++ )
    {
        max_val[i] = 0;
        for( k=0; k<L_SB; k++ )
        {
            inBuf_abs = (float)fabs(inBuf_pss[j]);
            if(max_val[i] < inBuf_abs )
            {
                max_val[i] = inBuf_abs;
            }

            j++;
        }
    }

    /* convert to maximum amplitude frequency log scale envelope */
    j = 0;
    for ( i=0; i<num_subband_smooth; i++ )
    {
        max_val_norm = 0.0f;
        if( max_val[i] != 0.0f )
        {
            max_val_norm = 10.0f / max_val[i];
        }

        for( k = 0; k < L_SB; k++ )
        {
            if(inBuf_pss[j] == 0.0f)
            {
                outBuf_pss[j] = 0.0f;
            }
            else if( fabs(inBuf_pss[j]) < max_val[i] )
            {
                outBuf_pss[j] = inBuf_pss[j] * max_val_norm;
            }
            else
            {
                /* CLIP , for avoiding computational difference */
                outBuf_pss[j] = 10.0f;
                if( inBuf_pss[j] < 0.0f )
                {
                    outBuf_pss[j] = -10.0f;
                }
            }
            j++;
        }
    }
    k = 0;
    m = 0;
    n = 0;
    reset_flag = 0;
    n_list[0] = 0;
    for(j=0; j<num_subband_smooth; j++)
    {
        cnt_zero_cont = 0;
        for(i=0; i<L_SB; i++)
        {
            if(outBuf_pss[k] == 0.0f)
            {
                cnt_zero_cont++;
            }
            else
            {
                cnt_zero_cont=0;
            }
            k++;
        }

        if(cnt_zero_cont != 0)
        {
            if(j > subband_search_offsets[0]/L_SB && reset_flag == 0)
            {
                n = 0;
                reset_flag = 1;
            }
            n_list[n] = j;
            n++;
        }

        if(reset_flag == 1 && n == 1)
        {
            m = 0;
        }

        if(cnt_zero_cont > 3*L_SB/4)
        {
            for(i=0; i<L_SB; i++)
            {
                if(outBuf_pss[k-L_SB+i] == 0.0f)
                {
                    outBuf_pss[k-L_SB+i] = outBuf_pss[n_list[m]*L_SB+i]*0.5f;
                }
            }
            m++;
        }
    }

    for(i=0; i<fLen; i++)
    {
        outBuf[i] = 0.0f;
        if(fabs(outBuf_pss[i]) > th_cut)
        {
            outBuf[i] = outBuf_pss[i];
        }
    }

    return;
}

/*-------------------------------------------------------------------*
 * Get20Log10Spec()
 *
 * Calculates 20*log10() for the specified samples. Input and output buffers can be the same.
 *-------------------------------------------------------------------*/
static
void Get20Log10Spec(
    const float *inBuf,         /* i  : input           */
    float *outBuf,        /* o  : output          */
    const short fLen            /* i  : loop length     */
)
{
    short i;

    for( i = 0; i < fLen; i++, inBuf++ )
    {
        *outBuf++ = (float) (20.0f * log10(fabs(*inBuf + 1.0)));
    }

    return;
}

void convert_lagIndices_pls2smp(
    short lagIndices_in[],
    short nBands_search,
    short lagIndices_out[],
    const float sspectra[],
    const short sbWidth[],
    const short fLenLow
)
{
    int sb;
    int i, cnt;

    for( sb = 0; sb < nBands_search; sb++ )
    {
        cnt = 0;
        i = 0;
        while( cnt <= lagIndices_in[sb] )
        {
            if( sspectra[subband_search_offsets[sb]+i] != 0.0f )
            {
                cnt++;
            }

            i++;

            if( subband_search_offsets[sb] + i + sbWidth[sb] >= fLenLow )
            {
                /* over fLenLow, no need for more search */
                break;
            }
        }

        lagIndices_out[sb] = (short)i-1+subband_search_offsets[sb];
    }

    return;
}

int get_usebit_npswb(
    short hqswb_clas
)
{
    short i;
    short bits;
    short up_lmt;
    const short *bits_req;

    up_lmt = 0;
    bits_req = bits_lagIndices;
    bits =0;
    if(hqswb_clas == HQ_NORMAL)
    {
        up_lmt = NB_SWB_SUBBANDS;
        bits_req = bits_lagIndices;
    }
    else if ( hqswb_clas == HQ_HARMONIC )
    {
        up_lmt = NB_SWB_SUBBANDS_HAR_SEARCH_SB;
        bits_req =bits_lagIndices_mode0_Har;
        bits = 2; /*noise gain*/
    }
    for( i = 0; i < up_lmt; i++ )
    {
        bits += bits_req[i];
    }

    return bits;
}

void SpectrumSmoothing_nss(
    float *inBuf,
    float *outBuf,
    int fLen
)
{
    int i, k;
    float inBufw[L_FRAME32k+L_SB_NSS];
    float outBufw[L_FRAME32k+L_SB_NSS];
    float temp_sum_1[NUM_SUBBAND_SMOOTH_MAX];
    float temp_sum_2[NUM_SUBBAND_SMOOTH_MAX];
    float temp_sum_3[NUM_SUBBAND_SMOOTH_MAX];
    float temp_sum_log[NUM_SUBBAND_SMOOTH_MAX];
    float temp_sum_smooth[NUM_SUBBAND_SMOOTH_MAX];
    float temp_sum_div[NUM_SUBBAND_SMOOTH_MAX];
    short num_subband_smooth;
    float smr, avg_val, r0;
    float clip_cof, avg_val2, thre;
    float thre_min;
    float max_peak;

    /* calculate subband number */
    num_subband_smooth = (short)(fLen/L_SB_NSS);
    if( fLen/(L_SB_NSS + 0.0f) - num_subband_smooth > 0.0f )
    {
        num_subband_smooth++;
    }

    /* buffer copy for fractional point */
    for( i = 0; i < fLen; i++ )
    {
        inBufw[i] = inBuf[i];
        outBufw[i] = 0.0f;
    }

    /* initialization of over buffer for fractional point */
    for( i = fLen; i < fLen + L_SB_NSS; i++ )
    {
        inBufw[i] = 0.0f;
        outBufw[i] = 0.0f;
    }

    avg_val = 0.0f;
    for( i = 0; i < fLen; i++ )
    {
        r0 = (float)fabs(inBufw[i]);
        avg_val += r0;
    }
    avg_val /= (float)fLen;

    max_peak = 0.0f;
    for( i = 0; i < fLen; i++ )
    {
        r0 = (float)fabs(inBufw[i]);
        if( max_peak < r0 )
        {
            max_peak = r0;
        }
    }

    smr = 10.0f * (float)log10( max_peak/(avg_val + 1.0e-20) + 1.0e-20 );

    for( i = 0; i < num_subband_smooth; i++ )
    {
        temp_sum_1[i] = 0.0;
        temp_sum_2[i] = 0.0;

        for( k = 0; k < L_SB_NSS_HALF; k++ )
        {
            temp_sum_1[i] += (float)(fabs(inBufw[k + L_SB_NSS * i]));
        }

        for( k = L_SB_NSS_HALF; k < L_SB_NSS; k++ )
        {
            temp_sum_2[i] += (float)(fabs(inBufw[k + L_SB_NSS * i]));
        }

        temp_sum_1[i] *= 0.25f;
        temp_sum_2[i] *= 0.25f;
        temp_sum_3[i] = temp_sum_1[i] * temp_sum_2[i];

        if( temp_sum_3[i] == 0.0f )
        {
            temp_sum_3[i] = temp_sum_1[i] + temp_sum_2[i];
        }
    }

    Get20Log10Spec( temp_sum_3, temp_sum_log, num_subband_smooth );

    for( i = 0; i < num_subband_smooth; i++ )
    {
        temp_sum_log[i] *= 0.5f;
    }

    SmoothSpec( temp_sum_log, temp_sum_smooth, num_subband_smooth, MA_LEN );

    for (i = 0; i < num_subband_smooth; i++)
    {
        temp_sum_div[i] = (float)pow(10, -1.0f * temp_sum_smooth[i]/20);
    }

    for (i = 0; i < num_subband_smooth; i++)
    {
        for (k = 0; k < L_SB_NSS; k++)
        {
            outBufw[k + L_SB_NSS * i] = inBufw[k + L_SB_NSS * i] * temp_sum_div[i];
        }
    }

    avg_val2 = 0.0f;
    for( i = 0; i < fLen; i++ )
    {
        r0 = (float)fabs(outBufw[i]);
        avg_val2 += r0;
    }
    avg_val2 /= (float)fLen;

    clip_cof = smr - 16.0f;
    if( clip_cof < 0.0f )
    {
        clip_cof = 0.0f;
    }
    clip_cof += 2.5f;

    thre = avg_val2 * clip_cof;
    thre_min = avg_val2 * 0.25f;
    for(i = 0; i < fLen; i++)
    {
        if( fabs(outBufw[i]) > thre )
        {
            if( outBufw[i] < 0.0f )
            {
                outBufw[i] = -1.0f*thre;
            }
            else
            {
                outBufw[i] = thre;
            }
        }

        if( fabs(outBufw[i]) < thre_min )
        {
            outBufw[i] = 0.0f;
        }
    }

    for(i = 0; i < fLen; i++)
    {
        outBuf[i] = outBufw[i];
    }

    return;
}

static float get_sigma(
    float x_abs[],
    float avg,
    int length
)
{
    int i;
    float d;
    float sigma;

    d = 0;
    for( i=0; i<length; i++ )
    {
        d += x_abs[i] * x_abs[i];
    }

    d /= (length-1);
    d -= avg*avg;

    sigma = (float)sqrt(d);

    return sigma;
}

/*--------------------------------------------------------------------------*
 * FindNBiggest2_simple()
 *
 * Finds N biggest components from input
 * Maximum value allowed for nIdx is currently 140 and the maximum value of n is currently 20
 *--------------------------------------------------------------------------*/

void FindNBiggest2_simple(
    const float *inBuf,          /* i  : input buffer (searched)                     */
    GainItem *g,              /* o  : N biggest components found                  */
    const short nIdx,            /* i  : search length                               */
    short *n,              /* i  : number of components searched (N biggest)   */
    short N_NBIGGESTSEARCH
)
{
    short j;

    float abs_in[400];
    float avg_in;
    float max_in;
    float thr;
    short   peak_cnt;
    float sigma;

    max_in = 0;
    avg_in = 0;
    for (j = 0; j < nIdx; j++)
    {
        abs_in[j] = (float)fabs(inBuf[j]);

        if( max_in < abs_in[j] )
        {
            max_in = abs_in[j];
        }

        avg_in += abs_in[j];
    }

    avg_in /= (float)nIdx;

    peak_cnt = 0;
    if( max_in <= 0.0001f )
    {
        for (j = 0; j < N_NBIGGESTSEARCH; j++)
        {
            g[peak_cnt].nmrValue = 0.0f;
            g[peak_cnt].gainIndex = j;
            peak_cnt++;
        }
    }

    sigma = get_sigma( abs_in, avg_in, nIdx );
    thr = avg_in + sigma * 1.15f;

    if( peak_cnt < N_NBIGGESTSEARCH )
    {
        for (j = 0; j < nIdx; j++)
        {
            if( abs_in[j] > thr )
            {
                g[peak_cnt].nmrValue = abs_in[j];
                g[peak_cnt].gainIndex = j;
                abs_in[j] = 0.0f;
                peak_cnt++;
            }

            if( peak_cnt == N_NBIGGESTSEARCH )
            {
                break;
            }
        }
    }

    thr *= (0.3f / N_NBIGGESTSEARCH) * peak_cnt + 0.7f;

    if( peak_cnt < N_NBIGGESTSEARCH )
    {
        for (j = 0; j < nIdx; j++)
        {
            if( abs_in[j] > thr )
            {
                g[peak_cnt].nmrValue = abs_in[j];
                g[peak_cnt].gainIndex = j;
                abs_in[j] = 0.0f;
                peak_cnt++;
            }

            if( peak_cnt == N_NBIGGESTSEARCH )
            {
                break;
            }
        }
    }

    thr *= (0.6f / N_NBIGGESTSEARCH) * peak_cnt + 0.3f;
    if( peak_cnt < N_NBIGGESTSEARCH )
    {
        for (j = 0; j < nIdx; j++)
        {
            if( abs_in[j] > thr )
            {
                g[peak_cnt].nmrValue = abs_in[j];
                g[peak_cnt].gainIndex = j;
                abs_in[j] = 0.0f;
                peak_cnt++;
            }

            if( peak_cnt == N_NBIGGESTSEARCH )
            {
                break;
            }
        }
    }

    *n = peak_cnt;

    return;
}

/*--------------------------------------------------------------------------*
 * spectrumsmooth_noiseton()
 * Spectrum normalization for the the core coder
 *--------------------------------------------------------------------------*/
float spectrumsmooth_noiseton(
    float spectra[],          /* i  : core coder                                  */
    const float spectra_ni[],       /* i  : core coder with sparse filling              */
    float sspectra[],         /* o  : Smoothed tonal information from core coder  */
    float sspectra_diff[],    /* o  : non tonal infomration for gap filling       */
    float sspectra_ni[],      /* o  : smoothed core coder                         */
    const short fLenLow,            /* i  : low frequency boundaries                    */
    short *ni_seed            /* io : random seed                                 */
)
{
    float spectra_diff[L_FRAME32k];
    float ni_ratio,ss_min,cut_sig_th,cut_ni_th;
    short i,pcnt,sign;
    float  spectra_rm[L_FRAME32k];
    float  cut_input=0.1f;
    float  rand_a[L_FRAME32k];

    /* pre-prepare random array for float-fix interoperability */
    for(i=0; i<fLenLow; i++)
    {
        rand_a[i] = own_random(ni_seed)/32768.0f;
    }

    /*Get the pulse resolution for the core coder*/
    pcnt = 0;
    for(i=0; i<fLenLow; i++)
    {
        if(spectra[i] != 0.0)
        {
            pcnt++;
        }
    }

    ni_ratio = 4.0f*(pcnt)/(fLenLow+0.0f);
    ni_ratio = min(0.9f, ni_ratio);

    ss_min = ni_ratio*10.0f;
    cut_sig_th = ss_min/4.0f;
    cut_sig_th = max(0.95f, cut_sig_th);
    /*core coder normalization for gap filling*/
    for(i=0; i<fLenLow; i++)
    {
        spectra_rm[i] = 0.0f;
        if( fabs(spectra[i]) > cut_input )
        {
            spectra_rm[i] = spectra[i];
        }
    }
    SpectrumSmoothing( spectra_rm, sspectra, fLenLow, cut_sig_th);
    /*Extract noise informaton from the core coder*/
    mvr2r(sspectra, sspectra_ni, fLenLow);
    for(i=0; i<fLenLow; i++)
    {
        spectra_diff[i] = spectra_ni[i] - spectra[i];
    }
    cut_ni_th = 0.0f;
    /*normalize sparse filled components*/
    for(i=0; i<fLenLow; i++)
    {
        spectra_rm[i] = 0.0f;
        if( fabs(spectra_diff[i]) > cut_input )
        {
            spectra_rm[i] = spectra_diff[i];
        }
    }
    SpectrumSmoothing( spectra_rm, sspectra_diff, fLenLow, cut_ni_th);
    /*Normalized corecoder for Gap filling */
    for(i=0; i<fLenLow; i++)
    {
        sign = 1;
        if(sspectra[i] < 0)
        {
            sign=-1;
        }
        if(fabs(sspectra[i]) > ss_min )
        {
            sspectra[i] =sign*((10-ss_min)/10.0f*(float)fabs(sspectra[i])+ss_min);
        }
        if(sspectra[i] != 0.0)
        {
            sspectra_ni[i] = sspectra[i];
        }
        else
        {
            sspectra_ni[i] = sspectra_diff[i]*ni_ratio;
        }
        if( sspectra_ni[i] == 0.0f)
        {
            sspectra_ni[i] = 0.5f*10.0f*ni_ratio* rand_a[i];
        }
    }
    return(ss_min);
}

/*--------------------------------------------------------------------------*
 * noiseinj_hf()
 * level adjustments for the missing bands in the core coder
 *--------------------------------------------------------------------------*/

void noiseinj_hf(
    float xSynth_har[],             /* o  : gap filled information            */
    float th_g[],                   /* i  : level adjustment information      */
    float band_energy[],            /* i  : subband energies                  */
    float *prev_En_sb,              /* i/o: band Energies                     */
    const short p2a_flags[],              /* i  : Missing bands in the core coder   */
    short BANDS,                    /* i  : total bands                       */
    short band_start[],             /* i  : band start indices                */
    short band_end[],               /* i  : band end indices                  */
    const short fLenLow                   /* i  : low frequency bandwidth           */
)
{

    float *p_En,ni_scale,*p_Enn_sm_sb,En[NB_SWB_SUBBANDS],Enn_sm_sb[NB_SWB_SUBBANDS];
    short k,i;

    short map_pulse_t[L_FRAME32k];
    short map_pulse[L_FRAME32k];

    set_s(map_pulse_t, 0, band_end[BANDS-1]+1);
    set_s(map_pulse, 0, band_end[BANDS-1]+1);
    /*level adjust the missing bands in the core coder */
    p_En = En;
    for(k=BANDS-NB_SWB_SUBBANDS; k<BANDS; k++)
    {
        *p_En = 0.0f;
        if(p2a_flags[k] == 0)
        {
            for(i=band_start[k]; i<=band_end[k]; i++)
            {
                if( fabs(xSynth_har[i-fLenLow]) <= th_g[k-(BANDS-NB_SWB_SUBBANDS)] )
                {
                    *p_En += xSynth_har[i-fLenLow]*xSynth_har[i-fLenLow];
                }
                else
                {
                    map_pulse_t[i] = 1;
                }
            }
            *p_En = (float)sqrt(*p_En);
        }
        p_En++;
    }

    p_En = En;
    p_Enn_sm_sb = Enn_sm_sb;
    for(k=BANDS-NB_SWB_SUBBANDS; k<BANDS; k++)
    {
        *p_Enn_sm_sb = prev_En_sb[k-(BANDS-NB_SWB_SUBBANDS)];
        if(p2a_flags[k] == 0)
        {
            if( prev_En_sb[k-(BANDS-NB_SWB_SUBBANDS)] < 0.8f*band_energy[k] )
            {
                *p_Enn_sm_sb = (0.15f*(*p_En)) + (0.85f*prev_En_sb[k-(BANDS-NB_SWB_SUBBANDS)]);
            }
            else
            {
                *p_Enn_sm_sb = (0.8f*(*p_En)) + (0.2f*prev_En_sb[k-(BANDS-NB_SWB_SUBBANDS)]);
            }
        }

        p_Enn_sm_sb++;
        p_En++;
    }

    p_En = En;
    p_Enn_sm_sb = Enn_sm_sb;
    map_pulse[fLenLow] = (map_pulse_t[fLenLow] | map_pulse_t[fLenLow+1]);
    for(i=fLenLow+1; i<band_end[BANDS-1]; i++)
    {
        map_pulse[i] = ( map_pulse_t[i-1] | map_pulse_t[i] | map_pulse_t[i+1] );
    }
    map_pulse[i] = (map_pulse_t[i-1] | map_pulse_t[i]);
    for(k=BANDS-NB_SWB_SUBBANDS; k<BANDS; k++)
    {
        if(p2a_flags[k] == 0 && *p_En != 0.0f)
        {
            ni_scale = sqrt( (*p_Enn_sm_sb)/(*p_En) );
            ni_scale = min(1.25f, ni_scale);
            ni_scale = max(0.75f, ni_scale);
            ni_scale *= 0.8f;
            for(i=band_start[k]; i<=band_end[k]; i++)
            {
                if( fabs(xSynth_har[i-fLenLow]) <= th_g[k-(BANDS-NB_SWB_SUBBANDS)] )
                {
                    if(map_pulse[i] == 0)
                    {
                        xSynth_har[i-fLenLow] *= ni_scale;
                    }
                }
            }
            prev_En_sb[k-(BANDS-NB_SWB_SUBBANDS)] = *p_Enn_sm_sb;
        }
        p_Enn_sm_sb++;
        p_En++;
    }
    return;
}

/*--------------------------------------------------------------------------*
 * noise_extr_corcod()
 * Spectrum normalization for the core coder
 *--------------------------------------------------------------------------*/

void noise_extr_corcod(
    float spectra[],                  /* i  : core coder                                  */
    const float spectra_ni[],               /* i  : core coder with sparse filling              */
    float sspectra[],                 /* o  : Smoothed tonal information from core coder  */
    float sspectra_diff[],            /* o  : non tonal infomration for gap filling       */
    float sspectra_ni[],              /* o  : smoothed core coder                         */
    const short fLenLow,                    /* i  : low frequency bands width                   */
    short prev_hqswb_clas,            /* i  : classification information                  */
    float *prev_ni_ratio              /* i  : noise parameter                             */
)
{
    short i,pulse_num;
    float spectra_diff[L_FRAME32k];
    float ni_ratio, ni_ratio_cur,br_adj;

    /*Spectrum Smoothing for tonal signals*/
    SpectrumSmoothing_nss( spectra, sspectra, fLenLow );

    mvr2r(sspectra, sspectra_ni, fLenLow);
    /*noise extraction*/
    for(i=0; i<fLenLow; i++)
    {
        spectra_diff[i] = spectra_ni[i] - spectra[i];
    }
    SpectrumSmoothing_nss( spectra_diff, sspectra_diff, fLenLow );
    /*Smoothing the noise components*/
    br_adj = 0.9f;
    pulse_num=0;
    for(i=0; i<fLenLow; i++)
    {
        if(spectra[i] != 0.0f)
        {
            pulse_num++;
        }
    }
    ni_ratio_cur = 0.0f;
    if(pulse_num != 0)
    {
        ni_ratio_cur = (fLenLow-pulse_num)/(fLenLow+0.0f);
        ni_ratio_cur *= br_adj;
    }
    if( prev_hqswb_clas == HQ_HARMONIC )
    {
        if(ni_ratio_cur >(*prev_ni_ratio))
        {
            ni_ratio = 0.8f*ni_ratio_cur +(*prev_ni_ratio)*0.2f;
        }
        else
        {
            ni_ratio = 0.6f*ni_ratio_cur +(*prev_ni_ratio)*0.4f;
        }
    }
    else
    {
        ni_ratio = 0.7f*ni_ratio_cur;
    }
    *prev_ni_ratio =ni_ratio;

    for(i=0; i<fLenLow; i++)
    {
        sspectra_diff[i] *= ni_ratio;
        sspectra_ni[i] = sspectra[i] + sspectra_diff[i];
    }
    return;
}

/*--------------------------------------------------------------------------*
 * ton_ene_est()
 * band energies for missing bands in the core coder
 *--------------------------------------------------------------------------*/

void ton_ene_est(
    float xSynth_har[],                  /* i  : buffer with non tonal compoents    */
    float be_tonal[],                    /* o  : tonal energy of the missing bands  */
    float band_energy[],                 /* i  : subband energies                   */
    short band_start[],                  /* i  : subband start indices              */
    short band_end[],                    /* i  : subband end indices                */
    short band_width[],                  /* i  : subband widths                     */
    const short fLenLow,                       /* i  : low frequency width                */
    const short fLenHigh,                      /* i  : High frequency width               */
    short bands,                         /* i  : total subbands                     */
    short har_bands,                     /* i  : total number of harmonics bands    */
    float ni_lvl,                        /* i  : noise enve for the hf bands        */
    GainItem pk_sf[],                    /* i  : subband widths                     */
    short *pul_res                       /* i  : tonal resolution                   */
)
{
    short sb_ton_loc[SWB_HAR_RAN1];
    float sb_ton[SWB_HAR_RAN1],peak[NB_SWB_SUBBANDS];
    short pos,count_pos_st,count_pos_end;
    short pul_res_bnd[NB_SWB_SUBBANDS];
    short k,i,j;
    float E=0.0f,E_r=0.0f,fac,ni_gain[NB_SWB_SUBBANDS],avg_pe[NB_SWB_SUBBANDS];

    set_s(sb_ton_loc, -1, SWB_HAR_RAN1);
    set_f(ni_gain,0.0f,NB_SWB_SUBBANDS);
    set_f(avg_pe,0.0f,NB_SWB_SUBBANDS);
    set_f(sb_ton, 0.0f,NB_SWB_SUBBANDS);
    set_f(peak, 0.0f,NB_SWB_SUBBANDS);
    /*non tonal adjustments*/
    for(i=0; i<fLenHigh; i++)
    {
        xSynth_har[i]*=ni_lvl;
    }

    pos=0;
    for(k=0; k<bands-har_bands; k++)
    {
        for(j=0; j<pul_res[k]; j++)
        {
            sb_ton_loc[pos] = pk_sf[k*8+j].gainIndex;
            sb_ton[pos] = pk_sf[k*8+j].nmrValue;
            pos++;
        }
    }
    k=0;
    pos=0;
    do
    {
        count_pos_st = pos;
        while(sb_ton_loc[pos] <=(band_end[k+har_bands]-fLenLow) && sb_ton_loc[pos]>=0 )
        {
            pos++;
        }
        count_pos_end = pos;
        pul_res_bnd[k] = count_pos_end-count_pos_st;
        if(pul_res_bnd[k] >0)
        {
            peak[k]=(float) fabs(sb_ton[count_pos_st]);
        }
        k++;
    }
    while(k<NB_SWB_SUBBANDS);

    k=0;
    /*energy calculation for tonal components*/
    for(i=har_bands; i<bands; i++)
    {
        E = sum2_f(&xSynth_har[band_start[i]-fLenLow],band_width[i]); /*noise energy*/
        E_r = (float) E/(float)pow(2.0f,band_energy[i]);

        if(E_r <0.06f)
        {
            avg_pe[k] = (float) sqrt(pow(2.0f,band_energy[i])/band_width[i]);
            fac = 0.6f;
            if( pul_res_bnd[k] != 0 )
            {
                fac = ((float)sqrt(E/band_width[i])/peak[k]);
            }
            ni_gain[k] = fac*avg_pe[k];

            ni_gain[k] = max(((ni_gain[k]*ni_gain[k]*E_r) >= 0.12f) ? 0.05f*ni_gain[k] : 1.0f*ni_gain[k], 1.4f);

            for(j=band_start[i]; j<=band_end[i]; j++)
            {
                xSynth_har[j-fLenLow] *=ni_gain[k];
            }
            E = sum2_f(&xSynth_har[band_start[i]-fLenLow],band_width[i]); /*noise energy*/
        }
        k++;
        be_tonal[i] = (float) pow (2.0f, band_energy[i])- E; /*tonal energy*/

        if(be_tonal[i] <0.0f)
        {
            E=0;
            for(j=(band_start[i]-fLenLow); j<=(band_end[i]-fLenLow); j++)
            {
                xSynth_har[j]*=0.25f;
                E+=xSynth_har[j]*xSynth_har[j];
            }
            be_tonal[i] = (float) pow (2.0f, band_energy[i])- E; /*tonal energy*/
        }
    }
    return;
}

/*--------------------------------------------------------------------------*
 * Gettonl_scalfact()
 * Gap filling for the core coder
 *--------------------------------------------------------------------------*/
void Gettonl_scalfact (
    float *outBuf,                   /* o  : synthesized spectrum                        */
    const float *codbuf,                   /* i  : core coder                                  */
    const short fLenLow,                   /* i  : lowband length                              */
    const short fLenHigh,                  /* i  : highband length                             */
    short harmonic_band,             /* i  : total number of Low frequency bands         */
    short bands,                     /* i  : total number of subbands in a frame         */
    float *band_energy,              /* i  : band energy of each subband                 */
    short *band_start,               /* i  : subband start indices                       */
    short *band_end,                 /* i  : subband end indices                         */
    const short p2aflags[],                /* i  : missing bands in the core coder             */
    float be_tonal[],                /* i  : tonal energy                                */
    GainItem *pk_sf,                 /* i  : toanl information for Sparse filling        */
    short *pul_res_pk                /* i  : pulse resolution information                */
)
{
    short k,i,band_pos;
    short sb_ton_loc[SWB_HAR_RAN1];
    short pos_tmp;
    float sb_ton[SWB_HAR_RAN1],est_ton_ene[NB_SWB_SUBBANDS],ton_sf;
    float step,enrd_r = 0.9f;
    float band_sf[SWB_HAR_RAN1];
    short pos,count_pos_st,count_pos_end,j;

    set_f(est_ton_ene,0.0f,NB_SWB_SUBBANDS);
    set_s(sb_ton_loc, -1,SWB_HAR_RAN1);
    /* Get the tonal information for sparse filling  */
    pos=0;
    for(k=0; k<bands-harmonic_band; k++)
    {
        for(j=0; j<pul_res_pk[k]; j++)
        {
            sb_ton_loc[pos] = pk_sf[k*8+j].gainIndex;
            sb_ton[pos] = pk_sf[k*8+j].nmrValue;
            pos++;
        }
    }
    k=0;
    pos=0;
    pos_tmp=0;

    do
    {
        band_pos = k+harmonic_band;
        count_pos_st = pos;
        while(sb_ton_loc[pos] <=(band_end[band_pos]-fLenLow) && sb_ton_loc[pos]>=0 )
        {
            pos++;
        }
        count_pos_end = pos;
        for(i=count_pos_st; i<count_pos_end; i++)
        {
            est_ton_ene[k] +=(sb_ton[i]*sb_ton[i]);
        }
        if(est_ton_ene[k] <=0.0f)
        {
            est_ton_ene[k] = 0.01f;
        }
        ton_sf = 0.0f;
        if(be_tonal[band_pos] > 0.0f)
        {
            ton_sf= (float) sqrt(be_tonal[band_pos]/est_ton_ene[k]);
        }
        for(i=count_pos_st; i<count_pos_end; i++)
        {
            band_sf[pos_tmp] = ton_sf;
            pos_tmp++;
        }
        k++;
    }
    while(k<NB_SWB_SUBBANDS);
    /* Gap filling for the core coder  */
    step = 1.0f / (0.077f * fLenHigh);
    pos_tmp=0;
    for(k=0; k<bands-harmonic_band; k++)
    {
        band_pos = k+harmonic_band;
        if( be_tonal[band_pos] > 0.0f )
        {
            enrd_r *= (float)sqrt(be_tonal[band_pos]/pow(2.0f,band_energy[band_pos]));
        }
        else
        {
            enrd_r = 0.0f;
        }
        enrd_r -=step;
        if(p2aflags[band_pos] == 1)
        {
            for(i= band_start[band_pos]; i<=band_end[band_pos]; i++)
            {
                outBuf[i-fLenLow] = codbuf[i];
            }
        }
        else
        {
            pos =0;
            pos+=pos_tmp;
            for(j=0; j<pul_res_pk[k]; j++)
            {
                outBuf[pk_sf[k*8+j].gainIndex] = pk_sf[k*8+j].nmrValue*band_sf[pos]*enrd_r;
                pos++;
            }
        }
        pos_tmp+=pul_res_pk[k];
    }

    return;

}


/*--------------------------------------------------------------------------*
 * return_bits_normal2()
 *
 *
 *--------------------------------------------------------------------------*/

void return_bits_normal2(
    short *bit_budget,
    const short p2a_flags[],
    const short bands,
    const short bits_lagIndices[]
)
{
    int i;

    for( i=0 ; i < NB_SWB_SUBBANDS; i++ )
    {
        if( p2a_flags[bands-NB_SWB_SUBBANDS+i] == 1 )
        {
            *bit_budget += bits_lagIndices[i];
        }
    }

    return;
}

/*--------------------------------------------------------------------------*
 * preset_hq2_swb()
 *
 *
 *--------------------------------------------------------------------------*/

void preset_hq2_swb
(
    const short hqswb_clas,
    const short band_end[],
    short *har_bands,
    short  p2a_bands,
    const short length,
    const short bands,
    short *lowlength,
    short *highlength,
    float m[]
)
{
    if( hqswb_clas == HQ_HARMONIC )
    {
        *har_bands = bands - p2a_bands + 1;
        *lowlength = band_end[*har_bands-1] + 1;
    }
    else
    {
        *lowlength = band_end[bands-NB_SWB_SUBBANDS-1] + 1;
    }

    *highlength = (length-*lowlength);
    set_f( m, 0.0f, length );

    return;
}


/*--------------------------------------------------------------------------*
 * post_hq2_swb()
 *
 *
 *--------------------------------------------------------------------------*/

void post_hq2_swb(
    const float m[],
    const short lowlength,
    const short highlength,
    const short hqswb_clas,
    const short har_bands,
    const short bands,
    const short p2a_flags[],
    const short band_start[],
    const short band_end[],
    float y2[],
    int npulses[]
)
{
    int  i, k;

    /* copy the scratch buffer to the output */
    mvr2r( &m[lowlength], &y2[lowlength], highlength );

    if( hqswb_clas == HQ_HARMONIC )
    {
        k = har_bands;
    }
    else
    {
        k = bands-NB_SWB_SUBBANDS;
    }

    for( ; k<bands; k++ )
    {
        if( p2a_flags[k] == 0 && npulses[k] == 0 )
        {
            for( i = band_start[k]; i <= band_end[k]; i++ )
            {
                if( y2[i] != 0.0f )
                {
                    npulses[k]++;
                }
            }
        }
    }

    return;
}


/*--------------------------------------------------------------------------*
 * GetSynthesizedSpecThinOut()
 *
 * Synthesize the spectrum in generic subband coding
 *--------------------------------------------------------------------------*/
void GetSynthesizedSpecThinOut(
    const float *predBuf,           /* i  : prediction buffer (i.e., lowband)   */
    float *outBuf,            /* o  : synthesized spectrum                */
    const short nBands,             /* i  : number of subbands calculated       */
    const short *sbWidth,           /* i  : subband lengths                     */
    const short *lagIndices,        /* i  : lowband index for each subband      */
    const float *lagGains,          /* i  : first gain for each subband         */
    const short predBufLen          /* i  : lowband length                      */
)
{
    short sb;
    short fLen, lag;
    float *ptr_in_outBuf;

    ptr_in_outBuf = outBuf;

    for( sb = 0; sb < nBands; sb++ )
    {
        fLen = sbWidth[sb];
        lag = lagIndices[sb];

        if( lag + fLen > predBufLen )
        {
            /* should never happen */
            lag = predBufLen - fLen;
        }

        GetPredictedSignal( predBuf, outBuf, lag, fLen, lagGains[sb]);
        outBuf += fLen;
    }

    outBuf = ptr_in_outBuf;

    return;
}


/*--------------------------------------------------------------------------*
 * GetlagGains()
 *
 *
 *--------------------------------------------------------------------------*/

void GetlagGains(
    const float *predBuf,              /* i: predictve buffer                */
    const float *band_energy,          /* i: band Energies                   */
    const short nBands,                /* i: high frequency bands            */
    const short *sbWidth,              /* i: high frequency band resolution  */
    const short *lagIndices,           /* i: correlation indices             */
    const short predBufLen,            /* i: predictive buffer length        */
    float *lagGains              /* o: lag gains                       */
)
{
    short i;
    short sb, fLen, lag;
    float outBuf[L_FRAME32k];
    float lagEnergy;

    /* Get the gain information for the missing bands*/
    for( sb = 0; sb < nBands; sb++ )
    {
        fLen = sbWidth[sb];
        lag = lagIndices[sb];

        if(lag + fLen > predBufLen)
        {
            /* should never happen */
            lag = predBufLen - fLen;
        }

        GetPredictedSignal( predBuf, outBuf, lag, fLen, 1.0 );

        lagEnergy = 0.0f;
        for( i=0; i<fLen; i++ )
        {
            lagEnergy += outBuf[i] * outBuf[i];
        }

        if( lagEnergy != 0.0f )
        {
            lagGains[sb] = (float)sqrt( pow(2.0f, band_energy[sb]) / lagEnergy );
        }
        else
        {
            lagGains[sb] = (float)sqrt( pow(2.0f, band_energy[sb]) / (lagEnergy+0.001f) );
        }
    }

    return;
}

/*--------------------------------------------------------------------------*
 * updat_prev_frm()
 *
 *
 *--------------------------------------------------------------------------*/
void updat_prev_frm(
    float y2[],                     /* i/o: core coder buffer                 */
    float t_audio[],                /*   o: core coder buffer                 */
    long bwe_br,                    /*   i: core bitrate                      */
    short length,                   /*   i: frame length coded bw             */
    const short inner_frame,              /*   i: input frame length                */
    short bands,                    /*   i: sub band resolution               */
    short bwidth,                   /*   i: NB/WB/SWB indicator               */
    const short is_transient,             /*   i: signal class information          */
    short hqswb_clas,               /*   i: signal class information          */
    short *prev_hqswb_clas,         /*   o: update signal class information   */
    short prev_SWB_peak_pos[],      /*   o: update core coder last coded peaks*/
    short prev_SWB_peak_pos_tmp[],  /*   o: update core coder last coded peaks*/
    short *prev_frm_hfe2,           /*   o: update harmonics                  */
    short *prev_stab_hfe2,          /*   o: update harmonics                  */
    short bws_cnt                   /*   i: band width detector               */
)
{
    short i,k,k1,k2,j;

    /* Copy the coded MDCT coefficient to the output buffer */
    if ( !is_transient )
    {
        /* Copy the scratch buffer to the output */
        mvr2r( y2, t_audio, length );

        /* If the input frame is larger than coded bandwidth, zero out uncoded MDCT coefficients */
        if ( inner_frame > length )
        {
            set_f( t_audio + length, 0.0f, inner_frame - length );
        }
    }
    else /* transient frame */
    {
        if( inner_frame == length || bws_cnt > 0)
        {
            /* Copy the scratch buffer to the output */
            mvr2r( y2, t_audio, length );
        }
        else
        {
            /* un-collapse transient frame and interleave zeros */
            for( i = 0; i < NUM_TIME_SWITCHING_BLOCKS; i++ )
            {
                k1 = i*length/NUM_TIME_SWITCHING_BLOCKS;
                k2 = i*inner_frame/NUM_TIME_SWITCHING_BLOCKS;

                mvr2r( y2 + k1, t_audio + k2, length/NUM_TIME_SWITCHING_BLOCKS );
                set_f( t_audio + k2 + length/NUM_TIME_SWITCHING_BLOCKS, 0.0f, (inner_frame-length)/NUM_TIME_SWITCHING_BLOCKS );
            }
        }
    }

    /* update */
    if( (bwe_br == HQ_16k40 || bwe_br == HQ_13k20) && bwidth == SWB )
    {
        *prev_hqswb_clas = hqswb_clas;
        if( hqswb_clas != HQ_HARMONIC )
        {
            *prev_frm_hfe2 = 0;
            *prev_stab_hfe2 = 0;
        }
    }
    else
    {
        *prev_hqswb_clas = is_transient;
    }

    if((bwe_br == HQ_16k40 ||bwe_br == HQ_13k20) && bwidth == SWB && hqswb_clas == HQ_NORMAL )
    {
        j = 0;
        for(k=bands-SPT_SHORTEN_SBNUM; k<bands; k++)
        {
            prev_SWB_peak_pos[j] = prev_SWB_peak_pos_tmp[j];
            j++;
        }
    }

    return;
}
