/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <math.h>
#include <stdlib.h>
#include "options.h"
#include "cnst.h"
#include "rom_com.h"
#include "prot.h"

/*--------------------------------------------------------------------------*
 * Local constants
 *--------------------------------------------------------------------------*/

#define NUMSF       8
#define LOG2_NUMSF       3
#define INV_NUMSF (float)0.125

/*--------------------------------------------------------------------------*
 * preecho_sb()
 *
 * Time-domain sub-band based pre-echo reduction
 *--------------------------------------------------------------------------*/

void preecho_sb(
    const long  brate,              /* i  : core bit-rate                                           */
    const float wtda_audio[],       /* i  : imdct signal                                            */
    float *rec_sig,           /* i  : reconstructed signal, output of the imdct transform     */
    const short framelength,        /* i  : frame length                                            */
    float *memfilt_lb,        /* i/o: memory                                                  */
    float *mean_prev_hb,      /* i/o: memory                                                  */
    float *smoothmem,         /* i/o: memory                                                  */
    float *mean_prev,         /* i/o: memory                                                  */
    float *mean_prev_nc,      /* i/o: memory                                                  */
    float *wmold_hb,          /* i/o: memory                                                  */
    short *prevflag,          /* i/o: flag                                                    */
    short *pastpre,           /* i/o: flag                                                    */
    const short bwidth
)
{
    short i, j, len3xLp20;
    float es_mdct[9];               /* 0..3 (0..7): energy of the 4 (8) subframes, 4..5: (8..10) energy of the future subframes */
    float es_mdct_hb[9];            /* 0..3 (0..7): energy of the 4 (8) subframes, 4..5: (8..10) energy of the future subframes */
    short zcr[9];                   /* 0..3 (0..7): zero crossing of the 4 (8) subframes, 4..5: (8..10) zero crossing of the future subframes */
    short maxnzcr[8], cntnzcr;      /* max number of samples without zero crossing  */
    float plus_es_mdct[64];         /* 8*8 subsubframes */
    float imdct_mem[L_FRAME48k];    /* memory of the imdct transform, used in the next frame   */

    short maxind, stind, stind_hb, cnt2, cnt5, adv, advmem;
    short ind2, ind3, pluslim;
    float max_es;
    float *ptr;
    float min_g[13], g, gt[13];
    float min_g_hb[13], gt_hb[13];
    float preechogain[960+PREECHO_SMOOTH_LEN];
    float preechogain_hb[960];
    short subframelength, subsubframelength;
    float maxcrit, savehalfe, savehalfe_hb;
    float rec_sig_lb[960], rec_sig_hb[960]; /* 960 max frame length at 48 kHz */
    float eshbmean2, eshbmean3, sxyhb2, sxylb3, sxyhb3;
    float *fptr1, *fptr2, *fptr3, *fptr4, *fptr5, *fptr6;
    short *sptr1, *sptr2;
    float wmold;
    float lim16, lim32;
    short limzcr, limmaxnzcr;
    float max_es_hb;
    short num_subsubframes, log2_num_subsubframes;
    float ftmp, fattnext;
    float sum_plus_es, mean_plus_es[65];
    float last2, last2_hb;
    float max_plus_es_mdct;
    short nb_flag, smooth_len;
    short firstnzcr;
    float es_mdct_half[9];
    float es_mdct_quart[9];
    double invsmoothlenp1; /*double just to have bit-exactness with the previous version, otherwise rounding difference of +-1*/
    short subframelength_s2, subframelength_s34;
    float oldgain, oldgain_hb;

    if( brate <= HQ_32k )
    {
        nb_flag = 0;
        if( bwidth == NB )
        {
            nb_flag = 1;
        }

        limzcr = 16;
        smooth_len = 4;
        invsmoothlenp1 = 0.2;

        if( nb_flag == 1 )
        {
            limzcr = 10;
            smooth_len = PREECHO_SMOOTH_LEN;
            invsmoothlenp1 = INV_PREECHO_SMOOTH_LENP1;
        }

        limmaxnzcr = framelength/24;
        num_subsubframes = 8;
        log2_num_subsubframes = 3;

        if( framelength == L_FRAME8k )
        {
            num_subsubframes = 4;
            log2_num_subsubframes = 2;
        }
        len3xLp20 = framelength/2-(short)((float)framelength*N_ZERO_MDCT_NS/FRAME_SIZE_NS);
        for( i = 0; i < len3xLp20; i++ )
        {
            imdct_mem[i] = -wtda_audio[len3xLp20 - 1 - i];
        }

        for( i = 0; i < framelength/2; i++ )
        {
            imdct_mem[len3xLp20 + i] = -wtda_audio[i];
        }


        subframelength = framelength >> LOG2_NUMSF;
        subsubframelength = subframelength >> log2_num_subsubframes;
        wmold = *smoothmem;
        subframelength_s2 = subframelength/2;
        subframelength_s34 = subframelength*3/4;

        cntnzcr = -1;

        lim16 = 0.1f;
        lim32 = 0.01f;
        savehalfe = 0.0f;
        savehalfe_hb = 0.0f;

        if( *pastpre == 0 )
        {
            /* if past frame mean energies are not known (no preecho_sb in the past frame), limit max attenuation to 1*/
            lim16 = (float)1;
            lim32 = (float)1;
        }

        *pastpre = 2;
        fptr1 = rec_sig_lb;
        fptr2 = rec_sig;
        fptr3 = rec_sig + 1;
        fptr4 = rec_sig + 2;
        *fptr1 = (float)(0.25*(*memfilt_lb + *fptr3) + 0.5 * *fptr2);
        fptr1++;

        for(j = 1; j < framelength - 1; j++)
        {
            *fptr1 = (float)(0.25*(*fptr2+*fptr4)+0.5* *fptr3);
            fptr1++;
            fptr2++;
            fptr3++;
            fptr4++;
        }

        *fptr1 = (float)(0.25*(*fptr2)+0.5* *fptr3);
        fptr1 = rec_sig_lb;
        fptr2 = rec_sig;
        fptr3 = rec_sig_hb;

        for(j = 0; j < framelength; j++)
        {
            *fptr3 = *fptr2  - *fptr1;
            fptr1++;
            fptr2++;
            fptr3++;
        }

        fptr2--;
        *memfilt_lb = *fptr2;

        /* energy of low bands 8 present and 1 future sub-frames */
        fptr1 = es_mdct;
        fptr5 = es_mdct_half;
        fptr6 = es_mdct_quart;
        fptr4 = es_mdct_hb;
        fptr2 = rec_sig;
        fptr3 = rec_sig_hb;
        sptr1 = zcr;
        *sptr1 = 0;
        sptr2 = maxnzcr;
        for (j = 0; j < NUMSF; j++)                 /* 8 present subframes */
        {
            *fptr1 = 100 + *fptr2 **fptr2;
            *fptr4 = 100 + *fptr3 **fptr3;

            *sptr2 = 0;
            firstnzcr = 1;
            if( j == 0 )
            {
                firstnzcr = 0;
            }

            fptr2++;
            fptr3++;

            for (i = 1; i < subframelength; i++)
            {
                if( i == subframelength_s2 )
                {
                    *fptr5 = *fptr1;
                }

                if( i == subframelength_s34 )
                {
                    *fptr6 = *fptr1;
                }
                *fptr1 += *fptr2 **fptr2;
                *fptr4 += *fptr3 **fptr3;
                cntnzcr++;
                if( *fptr2 **(fptr2-1) < 0 )
                {
                    (*sptr1)++;
                    if(cntnzcr > *sptr2)
                    {
                        *sptr2 = cntnzcr;
                    }

                    if( (firstnzcr > 0) && (cntnzcr > maxnzcr[j-1]) )
                    {
                        maxnzcr[j-1] = cntnzcr;
                    }

                    firstnzcr = 0;
                    cntnzcr = -1;
                }
                fptr2++;
                fptr3++;
            }

            if(cntnzcr > *sptr2)
            {
                *sptr2 = cntnzcr;
            }
            fptr4++;
            sptr1++;
            sptr2++;

            if( (firstnzcr > 0) && (cntnzcr > maxnzcr[j-1]) )
            {
                maxnzcr[j-1] = cntnzcr;
            }

            *sptr1 = 0;
            if((j < NUMSF-1) && ((*fptr2 **(fptr2-1) < 0))) /* zcr between 2 subframes */
            {
                (*sptr1)++;  /* counts for the nexte subframe */
                cntnzcr = -1;
            }

            if( *fptr5 < (*fptr1)/2 )
            {
                *fptr5 = 2*(*fptr1 - *fptr5);
            }
            else
            {
                *fptr5 = *fptr1;
            }
            fptr1++;
            fptr5++;
            fptr6++;
        }

        fptr2 = imdct_mem;
        j = NUMSF;
        *fptr1 = 100 + *fptr2 **fptr2;
        *sptr1 = 0;
        fptr2++;
        for (i = 1; i < len3xLp20; i++)  /* one future subframe but 140 samples (not 80) (enough with ALDO window) */
        {
            *fptr1 += *fptr2 **fptr2;

            if(*fptr2 **(fptr2-1) < 0)
            {
                (*sptr1)++;
            }

            fptr2++;
        }

        fptr2 = imdct_mem;
        fptr3 = imdct_mem + 1;
        fptr4 = imdct_mem + 2;
        ftmp = (float)(-0.25*(rec_sig[framelength-1] + *fptr3) + 0.5 * *fptr2);
        es_mdct_hb[NUMSF] =  100 + ftmp * ftmp;

        for(j = 1; j < len3xLp20 - 1; j++)
        {
            ftmp = (float)(-0.25*(*fptr2+*fptr4)+0.5* *fptr3);
            es_mdct_hb[NUMSF] +=  ftmp * ftmp;
            fptr2++;
            fptr3++;
            fptr4++;
        }

        ftmp = (float)(-0.25*(*fptr2)+0.5* *fptr3);
        es_mdct_hb[NUMSF] +=  ftmp * ftmp;
        max_es_hb = es_mdct_hb[0];                        /* for memorising the max energy */

        max_es = es_mdct[0];                        /* for memorising the max energy */
        maxind = 0;
        for (i = 1; i <= NUMSF; i++)
        {
            if (es_mdct_hb[i] >= max_es_hb)               /* '='  to handle the first window*/
            {
                max_es_hb = es_mdct_hb[i];                /* max energy low band, 8 present and 1 future subframes */
            }

            if (es_mdct[i] >= max_es)               /* '='  to handle the first window*/
            {
                max_es = es_mdct[i];                /* max energy low band, 8 present and 1 future subframes */
                maxind = i;
            }
        }

        cnt2 = cnt5 = 0;
        if( *prevflag != 0 || max_es < subframelength * 10000 )
        {
            maxind = 0;
        }

        if( max_es < 4 * *mean_prev )
        {
            maxind = 0;
        }
        *prevflag = 0;

        for (i = 0; i < maxind; i++)                /* only subbands before max energy subband are handled */
        {
            g = 1;                                  /* default gain */
            min_g[i] = 1;
            min_g_hb[i] = 1;
            if( (es_mdct_half[i] > 100*(*mean_prev_nc+500000))  || /* less then 20% energy in 3/4 of the subframe -> starting onset in the last quarter */
                    ((es_mdct_half[i] > 10*(*mean_prev_nc+500000)) &&
                     ((zcr[i] < limzcr) || (es_mdct_quart[i] < es_mdct[i]/6)) )) /* already an offset, plosif, do not touch */
            {
                maxind = i;                         /* no preecho reduction after the first subframe with gain 1 */
                *prevflag = 1;
                for(j = i-1; j >= 0; j--)
                {
                    if (es_mdct[j] > es_mdct[i] / 2)
                    {
                        maxind = j;
                    }
                }
            }
            else
            {
                if (es_mdct[i] < max_es / 16)
                {
                    g = lim16;
                    cnt5++;

                    if (es_mdct[i] < max_es / 32)
                    {
                        g = lim32;
                        cnt2++;
                    }

                    min_g[i] = (float)sqrt((double)(*mean_prev / es_mdct[i]));              /*limitation of attenuation gain */
                    min_g_hb[i] = (float)sqrt((double)(*mean_prev_hb / es_mdct_hb[i]));     /*limitation of attenuation gain */
                    if( (zcr[i] < limzcr/2) || (maxnzcr[i] > limmaxnzcr) )
                    {
                        if(min_g[i] < 1) /* *mean_prev < es_mdct[i]) */
                        {
                            *mean_prev = es_mdct[i];
                        }
                        min_g[i] = 1;               /* not noise-like, do not touch the amplitude, but may do in HB*/
                    }

                    if(min_g[i] > 1)                /* in fix point will be min_g = min(min_g, 1), weight 1 */
                    {
                        min_g[i] = (float)1;
                    }

                    if(min_g_hb[i] > 1)             /* in fix point will be min_g = min(min_g, 1), weight 1 */
                    {
                        min_g_hb[i] = (float)1;
                    }
                }
                else
                {
                    if( i > 0 && maxind < NUMSF )
                    {
                        *prevflag = 1;
                    }
                    maxind = i;                     /* no preecho reduction after the first subframe with gain 1*/
                }
            }
            gt[i] = g;
            gt_hb[i] = g;
        }

        for ( i = maxind; i <= NUMSF; i++ )         /* also for the first memory subframe */
        {
            gt[i] = 1;
            min_g[i] = 1;
            gt_hb[i] = 1;
            min_g_hb[i] = 1;

        }

        ind2 = 0;
        for( i = 0; i < NUMSF; i++ )
        {
            if( gt[i] < 1 )                         /*gt not yet limited by min_g*/
            {
                ind2 = i+1;                         /* first subframe with gain = 1 after last gain < 1 --> frame with the attack*/
            }
        }

        if( (wmold > 0.5) && ((cnt2+cnt5) < 2) )    /* mini either 1 cnt2 (and so also cnt5) or 2 cnt5 */
        {
            /* maxind = 0; false alarm, no echo reduction */
            ind2 = 0;
        }

        fptr3 = gt;
        fptr4 = gt_hb;
        fptr5 = min_g;
        fptr6 = min_g_hb;

        for (i = 0; i < ind2; i++) /* only subbands before max energy subband are handled*/
        {
            if(*fptr3 < *fptr5)                     /* in fix point will be g = max(g, min_g), weight 1 */
            {
                *fptr3 = *fptr5;
            }

            if(*fptr4 < *fptr6)                     /* in fix point will be g = max(g, min_g), weight 1 */
            {
                *fptr4 = *fptr6;
            }

            fptr1 = preechogain + i*subframelength+smooth_len;
            fptr2 = preechogain_hb + i*subframelength;
            for(j = 0; j < subframelength; j++)
            {
                *fptr1 = *fptr3;
                *fptr2 = *fptr4;
                fptr1++;
                fptr2++;
            }

            fptr3++;
            fptr4++;
            fptr5++;
            fptr6++;
        }

        max_plus_es_mdct = 0;
        adv = smooth_len;                                    /* samples needed to have near 1 gain after smoothing at the beggining of the attack subframe*/
        advmem = adv;

        if( ind2 > 0 || wmold < 1 || *wmold_hb < 1 )
        {
            ptr = imdct_mem;
            pluslim = num_subsubframes; /* if ind2 == NUMSF */

            if( ind2 < NUMSF )
            {
                ptr =  rec_sig + subframelength*ind2;
                pluslim = (NUMSF-ind2)*num_subsubframes;
            }

            maxcrit = *mean_prev_nc;
            if( ind2 == 0 )
            {
                sum_plus_es = *mean_prev_nc; /* 8 times mean sususb enenrgy (=maxcrit)*/
                pluslim = num_subsubframes;
                oldgain = wmold;
                oldgain_hb = *wmold_hb;
            }
            else /* ind2 > 0*/
            {
                sum_plus_es = es_mdct[ind2-1]; /* 8 times mean sususb enenrgy (=maxcrit)*/
                oldgain = gt[ind2-1];
                oldgain_hb = gt_hb[ind2-1];
                maxcrit = es_mdct[ind2-1]*gt[ind2-1]*gt[ind2-1];              /* /1 (iso /8) : 8 times of the pevious subframe mean*/

                if( (max_es/80 > maxcrit) && (zcr[ind2] > limzcr) )
                {
                    maxcrit = max_es/80;                /* still 10 times smaller then mean max_es*/
                }
            }

            for (j = 0; j < pluslim; j++)           /* 8  sub-subframes */
            {
                plus_es_mdct[j] = 100;
                for (i = 0; i < subsubframelength; i++)
                {

                    ftmp = *ptr **ptr;
                    plus_es_mdct[j] += ftmp;
                    if( plus_es_mdct[j] > max_plus_es_mdct )
                    {
                        max_plus_es_mdct = plus_es_mdct[j];
                    }
                    ptr++;
                }

                sum_plus_es += 4 * plus_es_mdct[j];
                mean_plus_es[j+1] = sum_plus_es / (j+2); /* 4/j could be tabulated in fixed point */
                if( mean_plus_es[j+1] < maxcrit)
                {
                    mean_plus_es[j+1] = maxcrit;
                }
            }
            mean_plus_es[0] = plus_es_mdct[0];
            if( mean_plus_es[0] < maxcrit)
            {
                mean_plus_es[0] = maxcrit;
            }
            mean_plus_es[pluslim] = -1;

            j = 0;
            while((plus_es_mdct[j] < mean_plus_es[j]) && (plus_es_mdct[j] < max_plus_es_mdct/8 ))
            {
                j++;
            }

            adv -= j*subsubframelength;

            if( ind2 < NUMSF )                    /* onset not in future frame */
            {
                fptr1 = preechogain + ind2*subframelength+smooth_len;
                fptr2 = preechogain_hb + ind2*subframelength;

                for (i = 0; i < j*subsubframelength; i++)
                {
                    *fptr1 = oldgain;        /*keep the gain of the previous subframe*/
                    *fptr2 = oldgain_hb;     /*keep the gain of the previous subframe*/
                    fptr1++;
                    fptr2++;
                }
            }
        }

        if( ind2 > 0 )
        {
            /* check increasing energy of preecho by regression last 3 subframes (if possible) */
            ind3 = ind2 + (j >> log2_num_subsubframes);  /* return (with rounding) to subframe basis */
            if( ind3 > 1 )
            {
                /* case of 3 points is simply */
                eshbmean2 = es_mdct_hb[ind3-1] + es_mdct_hb[ind3-2];
                sxyhb2 = (es_mdct_hb[ind3 - 1] - es_mdct_hb[ind3 - 2]); /* / eshbmean2 * 2; 04042013:  division not needed, only sign of sxyhb2 is used*/

                if( ind3 > 2 )
                {
                    eshbmean3 = (eshbmean2 + es_mdct_hb[ind3-3]) / 3;
                    sxylb3 = (es_mdct[ind3 - 1] - es_mdct[ind3 - 3]); /* /eslbmean3 / 2;           /2 for 3 points regression calc; 04042013:  division not needed, only sign of sxylb3 is used*/
                    sxyhb3 = (es_mdct_hb[ind3 - 1] - es_mdct_hb[ind3 - 3])/eshbmean3 / 2;
                    if ((sxyhb3 < (float)0.2) || (sxylb3 < 0))
                    {
                        ind2 = 0;
                        adv = advmem;
                    }
                }
                else
                {
                    if (sxyhb2 < (float)0.0)
                    {
                        ind2 = 0;
                        adv = advmem; /* 04042013: small bug corection*/
                    }
                }

                eshbmean3 = (eshbmean2 + es_mdct_hb[ind3]) / 3; /*verif toward future subsubframe*/
                sxyhb3 = (es_mdct_hb[ind3] - es_mdct_hb[ind3 - 2])/eshbmean3 / 2;
                if (sxyhb3 < (float)0.2)
                {
                    ind2 = 0;
                    adv = advmem;
                }
            }
        }

        stind = ind2*subframelength-adv;
        stind_hb = stind + advmem;
        if( stind < 0 )
        {
            stind = 0;
        }

        if( stind_hb < 0 )
        {
            stind_hb = 0;
        }

        fptr1 = preechogain + stind + smooth_len;
        fptr2 = preechogain_hb + stind_hb;

        for (i = stind + smooth_len; i < framelength; i++)       /* rest of the gains, without 4 (PREECHO_SMOOTH_LEN) 1 for fadeout */
        {
            *(fptr1++) = 1;
        }

        for (i = stind_hb; i < framelength; i++)    /* rest of the gains*/
        {
            *(fptr2++) = 1;
        }

        for (i = 0; i < smooth_len; i++)
        {
            preechogain[i] =  *smoothmem;
        }

        fattnext = 1;
        if( stind > framelength )
        {
            fattnext = gt[ind2-1];
        }

        for (i = 0; i < smooth_len; i++)
        {
            preechogain[framelength + i] = fattnext;
        }

        fptr1 = preechogain;
        for (i = 0; i < framelength; i++)
        {
            fptr2 = fptr1;
            for(j = 1; j <= smooth_len; j++)
            {
                *fptr1 += *(++fptr2);
            }

            *fptr1 *= (float)invsmoothlenp1;
            fptr1++;
        }

        *smoothmem = fattnext;
        *wmold_hb = preechogain_hb[framelength-1];

        /* apply gain */
        fptr1 = preechogain;
        fptr2 = preechogain_hb;
        fptr3 = rec_sig;
        fptr4 = rec_sig_lb;
        fptr5 = rec_sig_hb;

        for (i = 0; i < framelength; i++)
        {
            *fptr3 = *fptr4 **fptr1 + *fptr5 **fptr2;
            fptr1++;
            fptr2++;
            fptr3++;
            fptr4++;
            fptr5++;
        }

        *mean_prev_nc = es_mdct[0];                 /* compute mean not corrected by the actual gains*/

        for (i = 1; i < NUMSF; i++)                 /* all present subbands */
        {
            if( i == NUMSF/2 )
            {
                savehalfe = *mean_prev_nc;
            }
            *mean_prev_nc += es_mdct[i];
        }

        if( savehalfe < *mean_prev_nc/2 )
        {
            *mean_prev_nc = 2*(*mean_prev_nc - savehalfe);
        }
        *mean_prev_nc = *mean_prev_nc * INV_NUMSF; /* >> LOG2_NUMSF in fixpoint */

        for( i = 0; i < ind2; i++ )                 /* only subbands before max energy subband are handled*/
        {
            es_mdct[i] = es_mdct[i]*gt[i]*gt[i];
            es_mdct_hb[i] = es_mdct_hb[i]*gt_hb[i]*gt_hb[i];
        }

        *mean_prev = es_mdct[0];                    /* compute mean used in next frame to limit gain*/
        *mean_prev_hb = es_mdct_hb[0];              /* compute mean used in next frame to limit gain*/

        for (i = 1; i < NUMSF; i++)                 /* all present subbands */
        {
            if(i == NUMSF/2)
            {
                savehalfe = *mean_prev;
                savehalfe_hb = *mean_prev_hb;
            }

            *mean_prev += es_mdct[i];
            *mean_prev_hb += es_mdct_hb[i];
        }

        if( savehalfe < *mean_prev/2 )
        {
            *mean_prev = 2*(*mean_prev - savehalfe);
        }

        if( savehalfe_hb < *mean_prev_hb/2 )
        {
            *mean_prev_hb = 2*(*mean_prev_hb - savehalfe_hb);
        }

        *mean_prev = *mean_prev * INV_NUMSF; /* >> LOG2_NUMSF in fixpoint */
        *mean_prev_hb = *mean_prev_hb * INV_NUMSF; /* >> LOG2_NUMSF in fixpoint */
        last2 = (es_mdct[NUMSF-1] + es_mdct[NUMSF-2])/2;
        last2_hb = (es_mdct_hb[NUMSF-1] + es_mdct_hb[NUMSF-2])/2;

        if( last2 > *mean_prev )
        {
            *mean_prev = last2;
        }

        if( last2_hb > *mean_prev_hb )
        {
            *mean_prev_hb = last2_hb;
        }
    }

    return;
}


/*--------------------------------------------------------------------------*
 * inverse_transform()
 *
 * Inverse transform from the DCT domain to time domain
 *--------------------------------------------------------------------------*/

void inverse_transform(
    const float *in_mdct,           /* i  : input MDCT vector               */
    float *out,               /* o  : output vector                   */
    const short is_transient,       /* i  : transient flag                  */
    const short L,                  /* i  : output frame length             */
    const short L_inner             /* i  : length of the transform         */
)
{
    float out_alias[L_FRAME48k];
    float alias[MAX_SEGMENT_LENGTH];
    const float *in_segment;
    float *out_segment;
    float tmp;
    short ta, seg;
    short segment_length;
    float in_mdct_modif[L_FRAME48k];
    float *in_segment_modif;
    const float *win;

    segment_length = L/2;

    if (is_transient)
    {
        if (L == L_FRAME48k)
        {
            win = short_window_48kHz;
        }
        else if (L == L_FRAME32k)
        {
            win = short_window_32kHz;
        }
        else if( L == L_FRAME16k )
        {
            win = short_window_16kHz;
        }
        else  /* L == L_FRAME8k */
        {
            win = short_window_8kHz;
        }

        set_f( out_alias, 0.0f, L );

        in_segment = in_mdct;
        in_segment_modif = in_mdct_modif;

        if( L == L_inner )
        {
            mvr2r( in_mdct, in_mdct_modif, L );
        }
        else if( L > L_inner )
        {
            for( seg = 0; seg < NUM_TIME_SWITCHING_BLOCKS; seg++ )
            {
                for( ta = 0; ta < L_inner/NUM_TIME_SWITCHING_BLOCKS; ta++ )
                {
                    *in_segment_modif++ = *in_segment++;
                }

                for( ta = 0; ta < (L-L_inner)/NUM_TIME_SWITCHING_BLOCKS; ta++ )
                {
                    *in_segment_modif++ = 0;
                }
            }
        }
        else /* L < L_inner */
        {
            for( seg = 0; seg < NUM_TIME_SWITCHING_BLOCKS; seg++ )
            {
                for( ta = 0; ta < segment_length/2; ta++ )
                {
                    *in_segment_modif++ = *in_segment++;
                }
                in_segment += (L_inner-L)/NUM_TIME_SWITCHING_BLOCKS;
            }
        }

        out_segment = out_alias - segment_length/4;
        in_segment = in_mdct_modif;

        iedct_short( in_segment, alias, segment_length );

        for( ta = segment_length/4; ta < segment_length/2; ta++ )
        {
            out_segment[ta] = alias[ta];
        }

        for( ta = segment_length/2; ta < segment_length; ta++ )
        {
            out_segment[ta] = alias[ta] * win[ta];
        }

        out_segment = out_segment + segment_length/2;
        in_segment = in_segment + segment_length/2;

        for( seg = 1; seg <  NUM_TIME_SWITCHING_BLOCKS-1; seg++ )
        {
            iedct_short( in_segment, alias, segment_length );

            for( ta = 0; ta < segment_length; ta++ )
            {
                out_segment[ta] = out_segment[ta] + alias[ta] * win[ta];
            }

            in_segment = in_segment + segment_length/2;
            out_segment = out_segment + segment_length/2;
        }

        iedct_short( in_segment, alias, segment_length );

        for (ta = 0; ta < segment_length/2; ta++)
        {
            out_segment[ta] = out_segment[ta] + alias[ta] * win[ta];
        }

        for (ta = segment_length/2; ta < 3*segment_length/4; ta++)
        {
            out_segment[ta] = alias[ta];
        }

        for (ta = 0; ta < L/2; ta++)
        {
            tmp = out_alias[ta];
            out[ta] = out_alias[L-1-ta];
            out[L-1-ta] = tmp;
        }
    }
    else
    {
        edct( in_mdct, out, L );
    }

    return;
}

