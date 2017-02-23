/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#include <stdlib.h>
#include "prot.h"
#include "cnst.h"
#include "rom_com.h"
#include "rom_enc.h"

/*---------------------------------------------------------------------*
 * Local constants
 *---------------------------------------------------------------------*/

#define PIT_MIN2     20      /* pit_min for pitch tracking                                                  */
#define PIT_MIN_1    44      /* pitch tracking */
#define PIT_MIN2_1   24
#define THR_relE     -11.0f
#define THRES0       1.17f   /* Threshold to favor smaller pitch lags                                       */
#define DELTA0       2.0f    /* initial range of multiples search                                           */
#define STEP         1.0f    /* increment in range of multiples search                                      */
#define THRES1       0.4f    /* Threshold to favor pitch lags coherence for neighbours                      */
#define DELTA_COH    14      /* Maximum pitch lags difference for neighbours to be considered as coherent   */
#define THRES3       0.7f    /* Threshold to favor pitch lags coherence with previous frames                */
#define CORR_TH0     0.4f    /* Noise threshold for reinforcement with past frame pitch                     */
#define CORR_TH1     0.5f    /* Noise threshold for reinforcement with neighbourhood pitch correlations     */
#define LEN_X        ((PIT_MAX/OPL_DECIM)-(PIT_MIN2/OPL_DECIM)+1)  /* Correlation buffer length             */
#define COH_FAC 1.4f         /* Factor for measuring the pitch coherence                                    */
#define NSECT 4
#define NHFR  3
#define L_FIR_PO  5
#define L_MEM  (L_FIR_PO-2)

/*-----------------------------------------------------------------*
 * Local function prototypes
 *-----------------------------------------------------------------*/

static void pitch_neighbour( const short sect0, const short pitch_tmp[], short pitch[NHFR][2*NSECT],
                             const float corr_tmp[], float corr[3][2*NSECT], const float thres1[2*NHFR], const short ind_tmp[2*NHFR] );

static void find_mult( float *fac, const short pitch0, const short pitch1, const short pit_max0,
                       float *corr, short *old_pitch, float *old_corr, float delta, const float step);

static short pitch_coherence( const short pitch0, const short pitch1, const float fac_max, const short diff_max );

static void lp_decim2( const float x[], float y[], const short l, float *mem );


/*-----------------------------------------------------------------*
 * pitch_ol_init
 *
 * Open loop pitch variable initialization
 *-----------------------------------------------------------------*/

void pitch_ol_init(
    float *old_thres,   /* o  : threshold for reinforcement of past pitch influence */
    short *old_pitch,   /* o  : pitch  of the 1st half-frame of previous frame      */
    short *delta_pit,   /* o  : pitch evolution extrapolation                       */
    float *old_corr     /* o  : correlation                                         */
)
{
    *old_thres = 0.0f;
    *old_pitch = 0;
    *delta_pit = 0;
    *old_corr  = 0.0f;

    return;
}

/*-----------------------------------------------------------------*
 * pitch_ol()
 *
 * Compute the open loop pitch lag.
 *
 * The pitch lag search is divided in three sections of two sets.
 * Each section cannot have a pitch multiple.
 * We find a maximum for each section.
 * We compare the maxima of each section.
 *
 * As there is a margin between section overlaps, especially for
 * longer delays, this section selection is more robust for not
 * to find multiples in the same section when pitch evolves rapidly
 *
 * For each section, the length of the vectors to correlate is
 * greater than or equal to the longest pitch delay
 *------------------------------------------------------------------*/

void pitch_ol(
    short pitch[3],       /* o  : open loop pitch lag for each half-frame                        */
    float voicing[3],     /* o  : maximum normalized correlation for each half-frame             */
    short *old_pitch,     /* i/o: OL pitch of the 2nd half-frame of the last frame               */
    float *old_corr,      /* i/o: correlation                                                    */
    float corr_shift,     /* i  : normalized correlation correction                              */
    float *old_thres,     /* i/o: maximum correlation weighting with respect to past frame pitch */
    short *delta_pit,     /* i/o: old pitch extrapolation correction (added to old pitch)        */
    float *st_old_wsp2,   /* i/o: weighted speech memory                                         */
    const float *wsp,           /* i  : weighted speech for current frame and look-ahead               */
    float mem_decim2[3],  /* i/o: wsp decimation filter memory                                   */
    const float relE,           /* i  : relative frame energy                                          */
    const short L_look,         /* i  : look-ahead                                                     */
    const short last_class,     /* i  : frame classification of last frame                             */
    const short bwidth,         /* i  : bandwidth                                                      */
    const short Opt_SC_VBR      /* i  : SC-VBR flag                                                    */
)
{
    float old_wsp2[(L_WSP-L_INTERPOL)/OPL_DECIM], *wsp2;
    float tmp_mem[3], scale1[2*DELTA_COH-1];
    float scaled_buf[2*LEN_X + 3*(DELTA_COH-1)];
    float cor_temp, cor_buf[2*LEN_X], *pt1, *pt2, *pt3, *pt4, *pt5 = 0, *pt6 = 0, *pt_cor0, *pt_cor1, *pt_cor2, *pt_cor3, *pt_cor4;
    float thres1[2*NHFR];
    short diff, cnt, ind, ind1, offset, offset1, offset_la = 0, offset_la1 = 0, coh_flag, coh_flag1;

    short i, j, k, m, pit_min, pit_min1, sect0, subsect0, old_tmp, old_tmp1, len_x, len_x1, len_temp, len_temp1;
    short pitchX[NHFR][2*NSECT], pitch_tmp[2*NHFR], ind_tmp[2*NHFR], tmp_buf[NHFR+1];
    float enr, enr1 = 0.0f, enr_norm[NSECT], enr_norm1[NSECT], fac;
    float scaledX[NHFR][2*NSECT], corX[NHFR][2*NSECT], cor_tmp[2*NHFR], cor_mean;

    const short *nb_sect, *nb_subsect, *len, *len1, *sublen, *sublen1, *pit_max, *sec_length, *sec_length1;
    short ind_corX, ind1_corX;
    float pit_min_coding;

    /*--------------------------------------------------------------*
     * Initialization
     *--------------------------------------------------------------*/

    nb_sect = nb_sect_12k8;
    nb_subsect = nb_subsect_12k8;

    len = len_12k8;
    len1 = len1_12k8;
    sublen = sublen_12k8;
    sublen1 = sublen1_12k8;
    pit_max = pit_max_12k8;
    sec_length = sec_length_12k8;
    sec_length1 = sec_length1_12k8;

    if( last_class < VOICED_TRANSITION && bwidth != NB )
    {
        /* reset last pitch reinforcement in case of unvoiced or transitions: it avoids some pitch doublings */
        *old_thres = 0.0f;
    }

    pit_min_coding = PIT_MIN_EXTEND;
    if( ( bwidth != NB && *old_pitch > PIT_MIN ) ||
            ( bwidth == NB && (*old_pitch > PIT_MIN2_1 || *old_thres < 0.1) ) )
    {
        pit_min = PIT_MIN/OPL_DECIM;
        pit_min1= PIT_MIN_1/OPL_DECIM;
        subsect0 = 2;
        sect0 = 1;
    }
    else
    {
        pit_min = PIT_MIN2/OPL_DECIM;
        pit_min1= PIT_MIN2_1/OPL_DECIM;
        subsect0 = 0 ;
        sect0 = 0;
    }

    len_x = PIT_MAX/OPL_DECIM - pit_min + 1;
    len_x1= PIT_MAX/OPL_DECIM - pit_min1 + 1;

    /*--------------------------------------------------------------*
     * Find decimated weighted speech
     * Update wsp buffer with the memory
     * decimation of wsp[] to search pitch in LF and to reduce complexity
     * Extend the decimation of wsp to the end of the speech buffer
     * Update wsp memory
     *--------------------------------------------------------------*/

    mvr2r( st_old_wsp2, old_wsp2, (L_WSP_MEM-L_INTERPOL)/OPL_DECIM );
    wsp2 = old_wsp2 + ((L_WSP_MEM-L_INTERPOL)/OPL_DECIM);
    lp_decim2( wsp, wsp2, L_FRAME, mem_decim2 );

    mvr2r( mem_decim2, tmp_mem, 3 );
    lp_decim2( &wsp[L_FRAME], &wsp2[L_FRAME/OPL_DECIM], L_look, tmp_mem );

    mvr2r( &old_wsp2[L_FRAME/OPL_DECIM], st_old_wsp2, (L_WSP_MEM-L_INTERPOL)/OPL_DECIM );

    /*-----------------------------------------------------------------*
     * attenuate the correlation correction factor due to noise
     * reset correlation buffer outside the useful range
     * Find the scaling functions for immediate neigbours and
     * further ones
     *-----------------------------------------------------------------*/

    corr_shift *= 0.5f;
    set_f( scaled_buf, 0, DELTA_COH-1 );
    set_f( scaled_buf + (DELTA_COH-1) + len_x, 0, DELTA_COH-1 );
    set_f( scaled_buf + 2*(DELTA_COH-1) + len_x + len_x1, 0, DELTA_COH-1 );

    pt1 = scale1 + DELTA_COH-1;
    pt2 = pt1;
    for( i=0 ; i< DELTA_COH ; i++ )
    {
        *pt1 = -(*old_thres)/DELTA_COH * i + *old_thres + 1.0f;
        *pt2-- = *pt1++;
    }

    /*-----------------------------------------------------------------*
     * Estimate the new pitch by extrapolating the old pitch value
     * for 2 half-frames
     *-----------------------------------------------------------------*/

    old_tmp = *old_pitch + *delta_pit;

    if( old_tmp > PIT_MAX/OPL_DECIM )
    {
        old_tmp = PIT_MAX/OPL_DECIM;
    }

    if( old_tmp < pit_min )
    {
        old_tmp = pit_min;
    }

    old_tmp1 = old_tmp + *delta_pit;
    if( old_tmp1 > PIT_MAX/OPL_DECIM )
    {
        old_tmp1 = PIT_MAX/OPL_DECIM;
    }

    if( old_tmp1 < pit_min )
    {
        old_tmp1 = pit_min;
    }


    /*-----------------------------------------------------------------*
     * Loop for all three half-frames (current frame + look-ahead)
     *-----------------------------------------------------------------*/
    pt_cor0 = scaled_buf + DELTA_COH-1;

    pt_cor2 = pt_cor0 - pit_min + old_tmp;
    pt_cor4 = pt_cor0 - pit_min1 + old_tmp + (DELTA_COH-1) + len_x;

    for( i=0 ; i<NHFR ; i++ )
    {
        pt1 = wsp2 + i*2*(L_SUBFR/OPL_DECIM);
        pt2 = pt1 - pit_min;
        enr = 0.01f;
        pt_cor1 = pt_cor0;
        pt4 = pt1 - pit_min1;
        pt_cor3 = pt_cor0 + (DELTA_COH-1) + len_x;

        /*-----------------------------------------------------------------*
         * First two half-frames (corresponding to current frame)
         *-----------------------------------------------------------------*/
        if( i < NHFR-1 )
        {
            pt3 = pt1;
            pt5 = pt1;

            for( j = sect0; j < nb_sect[i]; j++ )    /* loop for each section */
            {
                /*-----------------------------------------------------------------*
                 * Find fixed vector energy
                 *-----------------------------------------------------------------*/

                /* 1st set */
                k = (short)(pt1 - pt3);
                for( k = k+len[j]; k > 0; k-- )
                {
                    enr += *pt3 **pt3;
                    pt3++;
                }
                enr_norm[j] = enr;

                /* Reduce complexity (length of 'enr1' section is equal or larger than 'enr') */
                pt5 = pt3;
                enr1 = enr;

                /* 2nd set */
                k = (short)(pt1 - pt5);
                for( k = k+len1[j]; k > 0; k-- )
                {
                    enr1 += *pt5 **pt5 ;
                    pt5++ ;
                }
                enr_norm1[j] = enr1 ;
            }

            /*-----------------------------------------------------------------*
             * Find correlation for the non-overlapping pitch lag values
             *-----------------------------------------------------------------*/
            k = (short)(pt2 - pt1);
            for( k = k + pit_max[subsect0]; k >= 0; k-- )
            {
                *pt_cor1++ = dotp( pt1, pt2--, sublen[0] ) ;
            }

            /*-----------------------------------------------------------------*
             * For each subsection, find the correlation
             *-----------------------------------------------------------------*/
            for( j = subsect0; j < nb_subsect[i]; j++ )
            {
                len_temp = sublen[j];
                len_temp1= sublen1[j];

                k = (short)(pt2-pt1);
                if (len_temp < len_temp1)
                {
                    for ( k = k + pit_max[j+1]; k >= 0; k-- )
                    {
                        cor_temp = pt1[0] * pt2[0];
                        for ( m = 1; m < len_temp; m++ )
                        {
                            cor_temp += pt1[m] * pt2[m];
                        }
                        *pt_cor1++ = cor_temp;
                        for ( ; m < len_temp1; m++ )
                        {
                            cor_temp += pt1[m] * pt2[m];
                        }
                        *pt_cor3++ = cor_temp;
                        pt2--;
                    }
                }
                else
                {
                    for ( k = k + pit_max[j+1]; k >= 0; k-- )
                    {
                        cor_temp = pt1[0] * pt2[0];
                        for ( m = 1; m < len_temp1; m++ )
                        {
                            cor_temp += pt1[m] * pt2[m];
                        }
                        *pt_cor3++ = cor_temp;
                        for ( ; m < len_temp; m++ )
                        {
                            cor_temp += pt1[m] * pt2[m];
                        }
                        *pt_cor1++ = cor_temp;
                        pt2--;
                    }
                }
            }
        }

        /*-----------------------------------------------------------------*
         * Third half-frame (look-ahead)
         *-----------------------------------------------------------------*/

        else
        {
            /*-----------------------------------------------------------------*
             * For each section in both sets, find fixed vector energy
             *-----------------------------------------------------------------*/

            pt6 = pt1 + L_look/OPL_DECIM - 1;
            pt3 = pt6;
            pt5 = pt6;

            for( j = sect0; j < nb_sect[i]; j++ )    /* loop for each section */
            {
                /* 1st set */
                k = (short)(pt3 - pt6);
                for( k = k + len[j]; k > 0; k-- )
                {
                    enr += *pt3 **pt3;
                    pt3--;
                }

                enr_norm[j] = enr;

                pt5 = pt3;
                enr1 = enr;

                /* 2nd set */
                k = (short)(pt5 - pt6);
                for( k = k + len1[j]; k > 0; k-- )
                {
                    enr1 += *pt5 **pt5 ;
                    pt5-- ;
                }

                enr_norm1[j] = enr1 ;
            }

            /* Set pointers */
            if( sect0 == 0 )
            {
                pt2 = pt6 - pit_min;
                k = 2;
            }
            else
            {
                pt2 = pt6 - pit_max[1] - 1;
                k = pit_max[2] - pit_max[1];
            }

            /*-----------------------------------------------------------------*
             * Find correlation for the non-overlapping pitch lag values
             *-----------------------------------------------------------------*/

            for( ; k > 0; k-- )
            {
                *pt_cor1 = 0;
                for( m=0; m<sublen[0]; m++ )
                {
                    *pt_cor1 += pt6[-m] * pt2[-m];
                }
                pt_cor1++;
                pt2--;
            }

            /*-----------------------------------------------------------------*
             * For each subsection, find the correlation (overlapping pitch lag values)
             *-----------------------------------------------------------------*/

            for( j = subsect0; j < nb_subsect[i]; j++ )
            {
                len_temp = sublen[j];
                len_temp1= sublen1[j];

                k = pit_max[j+1] - pit_max[j];
                if (len_temp < len_temp1)
                {
                    for ( ; k > 0; k-- )
                    {
                        cor_temp = pt6[0] * pt2[0];
                        for ( m = 1; m < len_temp; m++ )
                        {
                            cor_temp += pt6[-m] * pt2[-m];
                        }
                        *pt_cor1++ = cor_temp;
                        for ( ; m < len_temp1; m++ )
                        {
                            cor_temp += pt6[-m] * pt2[-m];
                        }
                        *pt_cor3++ = cor_temp;
                        pt2--;
                    }
                }
                else
                {
                    for ( ; k > 0; k-- )
                    {
                        cor_temp = pt6[0] * pt2[0];
                        for ( m = 1; m < len_temp1; m++ )
                        {
                            cor_temp += pt6[-m] * pt2[-m];
                        }
                        *pt_cor3++ = cor_temp;
                        for ( ; m < len_temp; m++ )
                        {
                            cor_temp += pt6[-m] * pt2[-m];
                        }
                        *pt_cor1++ = cor_temp;
                        pt2--;
                    }
                }
            }
        }

        /* Save unscaled correlation vector  */
        mvr2r( pt_cor0, cor_buf, len_x );
        mvr2r( pt_cor0 + (DELTA_COH-1) + len_x, cor_buf + len_x, len_x1 );

        /*-----------------------------------------------------------------*
         * scale correlation function in the neighbourhood of
         * the extrapolated pitch
         *-----------------------------------------------------------------*/
        pt_cor1 = pt_cor2 - (DELTA_COH-1);
        pt_cor3 = pt_cor4 - (DELTA_COH-1) ;
        pt2 = scale1;

        for( k=0 ; k < 2*DELTA_COH-1 ; k++ )
        {
            *pt_cor1++ *= (*pt2);
            *pt_cor3++ *= (*pt2++);
        }

        /* update for next half-frame */
        pt_cor2 = pt_cor0 - pit_min + old_tmp1;
        pt_cor4 = pt_cor0 - pit_min1 + old_tmp1 + (DELTA_COH-1) + len_x;
        /*-----------------------------------------------------------------*
         * For each section, find maximum correlation and compute
         * normalized correlation
         *-----------------------------------------------------------------*/

        pt_cor1 = pt_cor0;
        offset = 0;
        pt_cor3 = pt_cor0 + (DELTA_COH-1) + len_x;
        offset1 = 0;
        for( j=sect0; j < nb_sect[i] ; j++ )       /* loop for each section */
        {
            /* 1st set */
            if( i==2 )
            {
                offset_la = L_look/OPL_DECIM - 1 - (len[j]-1);
            }
            else
            {
                offset_la = 0;
            }

            /* 2nd set */
            if( i==2 )
            {
                offset_la1 = L_look/OPL_DECIM - 1 - (len1[j]-1);
            }
            else
            {
                offset_la1 = 0;
            }

            /* 1st set of candidates */
            ind = maximum( pt_cor1, sec_length[j], 0 ) + offset;
            pitchX[i][j] = ind + pit_min;
            pt2 = pt1 - pitchX[i][j] + offset_la;  /* selected moving vector  */
            enr1 = dotp( pt2, pt2, len[j] ) + 0.01f;
            enr1 = inv_sqrt( enr_norm[j] * enr1 ); /* 1/sqrt(energy)    */
            corX[i][j] = cor_buf[ind] * enr1;      /* find best normalized correlation per section  */
            scaledX[i][j] = pt_cor0[ind] * enr1;   /* find best scaled normalized correlation per section  */

            pt_cor1 += sec_length[j];
            offset = offset + sec_length[j];

            /* 2nd set of candidates */
            ind1 = maximum( pt_cor3, sec_length1[j], 0 ) + offset1;
            pitchX[i][j+NSECT] = ind1 + pit_min1;
            pt4 = pt1 - pitchX[i][j+NSECT] + offset_la1;   /* selected moving vector  */
            enr1 = dotp( pt4, pt4, len1[j] ) + 0.01f;
            enr1 = inv_sqrt( enr_norm1[j] * enr1 );        /* 1/sqrt(energy) */
            corX[i][j+NSECT] = cor_buf[ind1+len_x] * enr1; /* find best normalized correlation per section */
            scaledX[i][j+NSECT] = pt_cor0[ind1+(DELTA_COH-1)+len_x] * enr1; /* find best scaled normalized correlation per section */

            pt_cor3 += sec_length1[j];
            offset1 = offset1 + sec_length1[j];
        }
    }
    /*-----------------------------------------------------------------*
     * Favor a smaller delay if it happens that it has its multiple
     * in the longer-delay sections  (harmonics check)
     *-----------------------------------------------------------------*/

    for( i=0; i<2; i++ )          /* loop for the 2 half-frames */
    {
        fac = THRES0;
        find_mult( &fac, pitchX[i][2], pitchX[i][3], pit_max[7], &scaledX[i][2], old_pitch, old_corr, DELTA0, STEP );      /* Multiples in longest-delay section */
        find_mult( &fac, pitchX[i][1], pitchX[i][2], pit_max[5], &scaledX[i][1], old_pitch, old_corr, DELTA0, STEP );      /* Multiples in 3rd section */

        if((sect0==0) && ((pitchX[i][0]*2) >= pit_min_coding))
        {
            find_mult( &fac, pitchX[i][0], pitchX[i][1], pit_max[3], &scaledX[i][0],
                       old_pitch, old_corr, DELTA0, STEP );                                                                /* Multiples in 2nd section */
        }
        fac = THRES0;
        find_mult( &fac, pitchX[i][NSECT+2], pitchX[i][NSECT+3], pit_max[7], &scaledX[i][NSECT+2], old_pitch, old_corr, DELTA0, STEP ) ; /* Multiples in longest-delay section */
        find_mult( &fac, pitchX[i][NSECT+1], pitchX[i][NSECT+2], pit_max[6], &scaledX[i][NSECT+1], old_pitch, old_corr, DELTA0, STEP ) ; /* Multiples in 3rd section */

        if( (sect0 == 0) && (pitchX[i][NSECT+0]*2 >= pit_min_coding) )
        {
            find_mult( &fac, pitchX[i][NSECT+0], pitchX[i][NSECT+1], pit_max[4], &scaledX[i][NSECT+0], old_pitch, old_corr, DELTA0, STEP ); /* Multiples in 2nd section */
        }
    }

    fac = THRES0;                 /*  the look-ahead  */
    find_mult( &fac, pitchX[i][2], pitchX[i][3], pit_max[7], &scaledX[i][2], old_pitch, old_corr, 2.0f, 2.0f );      /* Multiples in longest-delay section */
    find_mult( &fac, pitchX[i][1], pitchX[i][2], pit_max[5], &scaledX[i][1], old_pitch, old_corr, DELTA0, STEP );    /* Multiples in 3rd section */

    if( (sect0 == 0) && (pitchX[i][0]*2 >= pit_min_coding) )
    {
        find_mult( &fac, pitchX[i][0], pitchX[i][1], pit_max[3], &scaledX[i][0], old_pitch, old_corr, DELTA0, STEP );  /* Multiples in 2nd section */
    }

    fac = THRES0;                 /*  the look-ahead  */
    find_mult( &fac, pitchX[i][NSECT+2], pitchX[i][NSECT+3], pit_max[7], &scaledX[i][NSECT+2], old_pitch, old_corr, 2.0f, 2.0f ) ;   /* Multiples in longest-delay section */
    find_mult( &fac, pitchX[i][NSECT+1], pitchX[i][NSECT+2], pit_max[6], &scaledX[i][NSECT+1], old_pitch, old_corr, DELTA0, STEP ) ; /* Multiples in 3rd section */

    if( (sect0 == 0) && (pitchX[i][NSECT+0]*2 >= pit_min_coding))
    {
        find_mult( &fac, pitchX[i][NSECT+0], pitchX[i][NSECT+1], pit_max[4], &scaledX[i][NSECT+0], old_pitch, old_corr, DELTA0, STEP );                                                              /* Multiples in 2nd section */
    }

    /*-----------------------------------------------------------------*
     * Do 1st estimate for pitch values
     * Adjust the normalized correlation using estimated noise level
     * Compute the maximum scaling for the neighbour correlation
     * reinforcement
     *-----------------------------------------------------------------*/

    for( i=0; i<NHFR; i++ )
    {
        ind = maximum( scaledX[i]+sect0, (short)(NSECT-sect0), 0 );
        ind = ind + sect0;
        ind_tmp[i] = ind;
        pitch_tmp[i] = pitchX[i][ind];
        cor_tmp[i] = corX[i][ind];
        cor_tmp[i] += corr_shift;

        if( cor_tmp[i] > 1.0f )
        {
            cor_tmp[i] = 1.0f;
        }
        thres1[i] = THRES1 * cor_tmp[i];    /* Higher is the neighbour's correlation, higher is the weighting  */

        /* 2nd set of pitch candidates */
        ind1 = maximum( scaledX[i]+sect0+NSECT, (short)(NSECT-sect0), 0 );
        ind1 += (sect0 + NSECT);
        ind_tmp[i+NHFR] = ind1;
        pitch_tmp[i+NHFR] = pitchX[i][ind1];
        cor_tmp[i+NHFR] = corX[i][ind1];
        cor_tmp[i+NHFR] += corr_shift;

        if( cor_tmp[i+NHFR] > 1.0f )
        {
            cor_tmp[i+NHFR] = 1.0f;
        }
        thres1[i+NHFR] = THRES1 * cor_tmp[i+NHFR]; /* Higher is the neighbour's correlation, higher is the weighting */
    }
    /*-----------------------------------------------------------------*
     * Take into account previous and next pitch values of the present
     * frame and look-ahead. Choose the pitch lags and normalize
     * correlations for each half-frame & look-ahead
     *-----------------------------------------------------------------*/

    pitch_neighbour( sect0, pitch_tmp, pitchX, cor_tmp, scaledX, thres1, ind_tmp );
    for( i=0; i<NHFR; i++ )
    {
        ind = maximum( scaledX[i]+sect0, (short)(NSECT-sect0), 0 );
        ind = ind + sect0;
        ind_corX = maximum( corX[i]+sect0, (short)(NSECT-sect0), 0 );
        ind_corX = ind_corX + sect0;

        ind1 = maximum( scaledX[i]+sect0+NSECT, (short)(NSECT-sect0), 0 );
        ind1 += (sect0 + NSECT);
        ind1_corX = maximum( corX[i]+sect0+NSECT, (short)(NSECT-sect0), 0 );
        ind1_corX += (sect0+NSECT);

        if( scaledX[i][ind1] > scaledX[i][ind] )
        {
            ind = ind1;
        }

        if( Opt_SC_VBR && corX[i][ind1_corX] > corX[i][ind_corX] )
        {
            ind_corX = ind1_corX;
        }

        if( Opt_SC_VBR && (pitchX[i][ind]*0.4 < pitchX[i][ind_corX]) && (pitchX[i][ind]*0.6 > pitchX[i][ind_corX]) && (corX[i][ind_corX]>=0.9) ) /* && (pitchX[i][ind]>50)) */
        {
            pitch[i] = pitchX[i][ind_corX];
            voicing[i] = corX[i][ind_corX];
        }
        else
        {
            pitch[i] = pitchX[i][ind];
            voicing[i] = corX[i][ind];
        }
    }

    /*-----------------------------------------------------------------*
     * Increase the threshold for correlation reinforcement with
     * the past if correlation is high and pitch is stable
     *-----------------------------------------------------------------*/

    cor_mean = 0.5f * (voicing[0] + voicing[1]) + corr_shift;
    if( cor_mean > 1.0f )
    {
        cor_mean = 1.0f;
    }

    /* pitch unstable in present frame or from previous frame or normalized correlation too low  */
    coh_flag  = pitch_coherence( (short)pitch[0], (short)pitch[1], COH_FAC, DELTA_COH );
    coh_flag1 = pitch_coherence( (short)pitch[0], (short)*old_pitch, COH_FAC, DELTA_COH );
    if( ( coh_flag == 0 ) || ( coh_flag1 == 0 ) || ( cor_mean < CORR_TH0 ) || (relE < THR_relE) )
    {
        *old_thres = 0.0f;                          /* Reset the threshold          */
    }
    else
    {
        *old_thres += (0.16f * cor_mean);           /* The threshold increase is directly dependent on normalized correlation  */
    }

    if( *old_thres > THRES3 )
    {
        *old_thres = THRES3;
    }

    if( voicing[1] > voicing[0] )
    {
        *old_corr = voicing[1];
    }
    else
    {
        *old_corr = cor_mean;
    }

    /*-----------------------------------------------------------------*
     * Extrapolate the pitch value for the next frame by estimating
     * the pitch evolution. This value is added to the old_pitch
     * in the next frame and is then used when the normalized
     * correlation is reinforced by the past estimate
     *-----------------------------------------------------------------*/

    tmp_buf[0] = *old_pitch;
    for( i=0; i<NHFR; i++ )
    {
        tmp_buf[i+1] = pitch[i];
    }

    *delta_pit = 0;
    cnt = 0;
    for( i=0 ; i<NHFR ; i++ )
    {
        diff = tmp_buf[i+1] - tmp_buf[i];
        coh_flag = pitch_coherence( (short)tmp_buf[i], (short)tmp_buf[i+1], COH_FAC, DELTA_COH );
        if( coh_flag != 0 )
        {
            *delta_pit = *delta_pit + diff;
            cnt++;
        }
    }

    if( cnt == 2 )
    {
        *delta_pit /= 2;
    }

    if( cnt == 3 )
    {
        *delta_pit /= 3;
    }

    /*-----------------------------------------------------------------*
     * update old pitch, upsample pitch
     *-----------------------------------------------------------------*/

    *old_pitch = pitch[1];

    for( i=0; i<NHFR; i++ )
    {
        pitch[i] *= OPL_DECIM;
    }


    return;
}

/*-----------------------------------------------------------------*
 * find_mult
 *
 * Verifies whether max pitch delays in higher sections have multiples
 * in lower sections
 *-----------------------------------------------------------------*/

static void find_mult(
    float *fac,       /* i/o: correlation scaling factor                           */
    const short pitch0,     /* i  : pitch of max correlation in the c section            */
    const short pitch1,     /* i  : pitch of max correlation in the longer-delay section */
    const short pit_max0,   /* i  : max pitch delay in the longer-delay section          */
    float *corr,      /* i/o: max correlation in the shorter-delay section         */
    short *old_pitch, /* i  : pitch from previous frame                            */
    float *old_corr,  /* i  : max correlation from previous frame                  */
    float delta,      /* i  : initial multiples search range                       */
    const float step        /* i  : increment in range of multiples search               */
)
{
    short pit_min;

    pit_min = 2 * pitch0;         /* double the shorter-delay section pitch */
    while( pit_min <= pit_max0 + (short)delta )            /* check for section boundary */
    {
        if( abs( pit_min - pitch1 ) <= (short)delta )       /* if multiple in the allowed range */
        {
            if ( *old_corr < 0.6f || (float)pitch0 > (float)*old_pitch * 0.4f )
            {
                /* reinforce the normalized correlation */
                *corr *= *fac;
            }
            *fac *= THRES0;
        }
        pit_min = pit_min + pitch0;  /* next multiple */
        delta += step;               /* add the incertitude to the allowed range */
    }
    return;
}

/*---------------------------------------------------------------------------*
 * pitch_neighbour
 *
 * Verifies if the maximum correlation pitch lag is coherent with neighbour
 * values
 *---------------------------------------------------------------------------*/
static void pitch_neighbour(
    const short sect0,                /* i  : indicates whether section 0 (below PIT_MIN) is used     */
    const short pitch_tmp[],          /* i  : estimated pitch values for each half-frame & look-ahead */
    short pitch[NHFR][2*NSECT], /* i  : tested pitch values for each half-frame & look-ahead    */
    const float corr_tmp[],           /* i  : raw normalized correlation (before different scalings)  */
    float corr[NHFR][2*NSECT],  /* i/o: normalized correlation for each half-frame & look-ahead */
    const float thres1[2*NHFR],       /* i  : maximum scaling for the immediate neighbours            */
    const short ind_tmp[2*NHFR]       /* i  : maximum section indices                                 */
)
{
    short delta, i, j, k, K, coh_flag;

    for( k=sect0 ; k<NSECT ; k++ )      /* loop for each section */
    {
        if( k == (NSECT-1) )
        {
            K = 2;                      /* the number of tests depends on the section */
        }
        else
        {
            K = 3;
        }
        for( i=0; i<K; i++ )            /* loop for the 2 half-frames and the look-ahead  */
        {
            /* Compare pitch values of the present frame */
            for( j=0; j<K; j++ )        /* Verify pitch coherence with neighbours (including past pitch) */
            {
                if( j != i )            /* Exclude itself */
                {
                    if( corr_tmp[j] >= CORR_TH1 )       /* reinforcement can happen only if the correlation is high enough */
                    {
                        delta = (short)abs( pitch[i][k] - pitch_tmp[j] );       /* Find difference of pitch values  */
                        coh_flag = pitch_coherence( (short)pitch[i][k], (short)pitch_tmp[j], COH_FAC, DELTA_COH );
                        if( coh_flag != 0 )
                        {
                            /* Favour section-wise stability */
                            if ( ind_tmp[j] == k )
                            {
                                corr[i][k] *= ( -thres1[j]/DELTA_COH * delta + thres1[j]+1.0f );    /* Favour closer values */
                            }
                            else
                            {
                                corr[i][k] *= ( -thres1[j]/DELTA_COH * 0.625f * delta + thres1[j] * 0.625f +1.0f ) ;
                            }
                        }
                    }
                }
            }
        }
    }
    /*---------------------*
     * 2nd set of sections
     *---------------------*/
    for( k=sect0 ; k<NSECT ; k++ )      /* loop for each section */
    {
        if( k == (NSECT-1) )
        {
            K = 2;                      /* the number of tests depends on the section */
        }
        else
        {
            K = 3;
        }
        for( i=0; i<K; i++ )            /* loop for the 2 half-frames and the look-ahead */
        {
            /* Compare pitch values of the present frame */
            for( j=0; j<K; j++ )        /* Verify pitch coherence with neighbours (including past pitch) */
            {
                if( j != i )            /* Exclude itself */
                {
                    if( corr_tmp[j+NHFR] >= CORR_TH1 )           /* reinforcement can happen only if the correlation is high enough */
                    {
                        delta = (short)abs( pitch[i][NSECT+k] - pitch_tmp[j+NHFR] );       /* Find difference of pitch values */
                        coh_flag = pitch_coherence( (short)pitch[i][NSECT+k], (short)pitch_tmp[j+NHFR], COH_FAC, DELTA_COH );
                        if( coh_flag != 0 )
                        {
                            /* Favour section-wise stability */
                            if ( ind_tmp[j+NHFR] == (NSECT + k) )
                            {
                                corr[i][NSECT+k] *= ( -thres1[j+NHFR]/DELTA_COH * delta + thres1[j+NHFR]+1.0f );       /* Favour closer values */
                            }
                            else
                            {
                                corr[i][NSECT+k] *= ( -thres1[j+NHFR]/DELTA_COH * 0.625f * delta + thres1[j+NHFR] * 0.625f +1.0f ) ;
                            }
                        }
                    }
                }
            }
        }
    }
    return;
}

/*-----------------------------------------------------------------*
 * pitch_coherence
 *
 * Verify if pitch evolution is smooth
 *-----------------------------------------------------------------*/

static short pitch_coherence(
    const short pitch0,   /* i  : first pitch to compare         */
    const short pitch1,   /* i  : 2nd pitch to compare           */
    const float fac_max,  /* i  : max ratio of both values       */
    const short diff_max  /* i  : max difference of both values  */
)
{
    short smaller, larger;
    if( pitch1 < pitch0 )     /* Finds smaller and larger of 2 short values */
    {
        smaller = pitch1;
        larger = pitch0;
    }
    else
    {
        smaller = pitch0;
        larger = pitch1;
    }
    if( ((float)larger < fac_max * (float)smaller) && ( (larger-smaller) < diff_max) )
    {
        return 1;
    }
    else
    {
        return 0;
    }
}


/*-----------------------------------------------------------------*
 * lp_decim2:
 *
 * Decimate a vector by 2 with 2nd order fir filter.
 *-----------------------------------------------------------------*/
static void lp_decim2(
    const float x[],  /* i  : signal to process         */
    float y[],  /* o  : processed signals         */
    const short l,    /* i  : size of filtering         */
    float *mem  /* i/o: memory (size=3)           */
)
{
    const float *p_h;
    float temp, *p_x, x_buf[L_FRAME_PLUS+L_MEM];
    short i, j, k;

    /* copy initial filter states into buffer */
    p_x = x_buf;
    for (i=0; i<L_MEM; i++)
    {
        *p_x++ = mem[i];
    }
    for (i=0; i<l; i++)
    {
        *p_x++ = x[i];
    }
    for (i=0; i<L_MEM; i++)
    {
        mem[i] = x[l-L_MEM+i];
    }
    for (i=0, j=0; i<l; i+=2, j++)
    {
        p_x = &x_buf[i];
        p_h = h_fir;
        temp = 0.0f;

        for (k=0; k<L_FIR_PO; k++)
        {
            temp += *p_x++ **p_h++;
        }
        y[j] = temp;
    }
    return;
}

