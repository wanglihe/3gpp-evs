/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

#ifndef PROT_H
#define PROT_H

#include <stdio.h>
#include <stdarg.h>
#include "options.h"
#include "typedef.h"
#include "stat_enc.h"
#include "stat_dec.h"
#include "stat_com.h"
#include "cnst.h"
#include "stl.h"


/*----------------------------------------------------------------------------------*
 * Prototypes of global macros
 *----------------------------------------------------------------------------------*/

#ifndef min
#define min(x,y)                            ((x)<(y)?(x):(y))
#endif

#ifndef max
#define max(x,y)                            ((x)>(y)?(x):(y))
#endif

#define log_base_2(x)                       ((double)log((double)(x))*1.4426950408889634074f)
#define round_f(x)                          (((x)>0)?(int)((x)+0.5f):(-(int)((-x)+0.5f)))

#ifndef ABSVAL
#define ABSVAL(a)                           ((a)>=0?(a):(-(a)))
#endif

#ifndef SQR
#define SQR(a)                              ((a)*(a))
#endif

#ifndef SWAP
#define SWAP(a,b)                           {tempr=(a); (a)=(b); (b)=tempr;}
#endif

#ifndef swap
#define swap(x,y,type)                      {type u__p; u__p=x; x=y; y=u__p;}
#endif

#define set_max(a, b)         { if ((b) > *a) { *a = (b); } }   /* If the first argument is already the highes or lowest, nothing is done. */
#define set_min(a, b)         { if ((b) < *a) { *a = (b); } }   /* Otherwise, the 2nd arg is stored at the address of the first arg. */

static __inline Word16 L_Extract_lc(const Word32 L_32, Word16 *p_hi)
{
    *p_hi = extract_h(L_32);
    return lshr(extract_l(L_32), 1);
}

/*----------------------------------------------------------------------------------*
 * MODE1 prototypes
 *----------------------------------------------------------------------------------*/

float inv_sqrt(                             /* o  : inverse square root of input value */
    const float x                           /* i  : input value                        */
);

short own_random(                           /* o  : output random value */
    short *seed                             /* i/o: random seed         */
);

float sign(                                 /* o  : sign of x (+1/-1) */
    const float x                           /* i  : input value of x  */
);

float log2_f(                               /* o  : logarithm2 of x    */
    const float x                           /* i  : input value of x   */
);

short norm_ul(
    unsigned int UL_var1
);
short sum_s(                                /* o  : sum of all vector elements            */
    const short *vec,                       /* i  : input vector                          */
    const short lvec                        /* i  : length of input vector                */
);

float sum_f(                                /* o  : sum of all vector elements            */
    const float *vec,                       /* i  : input vector                          */
    const short lvec                        /* i  : length of input vector                */
);

float sum2_f(                               /* o  : sum of all squared vector elements    */
    const float *vec,                       /* i  : input vector                          */
    const short lvec                        /* i  : length of input vector                */
);

void set_c(
    char y[],                        /* i/o: Vector to set                       */
    const char a,                          /* i  : Value to set the vector to          */
    const short N                           /* i  : Lenght of the vector                */
);

void set_s(
    short y[],                        /* i/o: Vector to set                       */
    const short a,                          /* i  : Value to set the vector to          */
    const short N                           /* i  : Lenght of the vector                */
);

void set_i(
    int y[],                        /* i/o: Vector to set                       */
    const int a,                          /* i  : Value to set the vector to          */
    const short N                           /* i  : Lenght of the vector                */
);

void set_f(
    float y[],                        /* i/o: Vector to set                       */
    const float a,                          /* i  : Value to set the vector to          */
    const short N                           /* i  : Lenght of the vector                */
);

void set_zero(
    float *vec,                       /* o  : input vector         */
    int lvec                        /* i  : length of the vector */
);

void mvr2r(
    const float x[],                        /* i  : input vector  */
    float y[],                        /* o  : output vector */
    const short n                           /* i  : vector size   */
);

void mvs2s(
    const short x[],                        /* i  : input vector  */
    short y[],                        /* o  : output vector */
    const short n                           /* i  : vector size   */
);

unsigned int mvr2s(
    const float x[],                        /* i  : input vector  */
    short y[],                        /* o  : output vector */
    const short n                           /* i  : vector size   */
);

void mvs2r(
    const short x[],                        /* i  : input vector  */
    float y[],                        /* o  : output vector */
    const short n                           /* i  : vector size   */
);

void mvi2i(
    const int x[],                        /* i  : input vector  */
    int y[],                        /* o  : output vector */
    const int n                           /* i  : vector size   */
);
void AGC_dec(
    float x[],
    float mem[],
    const short n
);

short maximum(                              /* o  : index of the maximum value in the input vector */
    const float *vec,                       /* i  : input vector                                   */
    const short lvec,                       /* i  : length of input vector                         */
    float *max                        /* o  : maximum value in the input vector              */
);

short minimum(                              /* o  : index of the minimum value in the input vector */
    const float *vec,                       /* i  : input vector                                   */
    const short lvec,                       /* i  : length of input vector                         */
    float *min                        /* o  : minimum value in the input vector              */
);

short emaximum(                             /* o  : return index with max energy value in vector   */
    const float *vec,                       /* i  : input vector                                   */
    const short lvec,                       /* i  : length of input vector                         */
    float *ener_max                   /* o  : maximum energy value                           */
);

float mean(                                 /* o  : vector mean                            */
    const float *vec,                       /* i  : input vector                           */
    const short lvec                        /* i  : length of input vector                 */
);

float dotp(                                 /* o  : dot product of x[] and y[]    */
    const float  x[],                       /* i  : vector x[]                    */
    const float  y[],                       /* i  : vector y[]                    */
    const short  n                          /* i  : vector length                 */
);

void conv(
    const float x[],                        /* i  : input vector                              */
    const float h[],                        /* i  : impulse response (or second input vector) */
    float y[],                        /* o  : output vetor (result of convolution)      */
    const short L                           /* i  : vector size                               */
);

void fir(
    const float x[],                        /* i  : input vector                              */
    const float h[],                        /* i  : impulse response of the FIR filter        */
    float y[],                        /* o  : output vector (result of filtering)       */
    float mem[],                      /* i/o: memory of the input signal (M samples)    */
    const short L,                          /* i  : input vector size                         */
    const short K,                          /* i  : order of the FIR filter (M+1 coefs.)      */
    const short upd                         /* i  : 1 = update the memory, 0 = not            */
);

void v_add(
    const float x1[],                       /* i  : Input vector 1                       */
    const float x2[],                       /* i  : Input vector 2                       */
    float y[],                        /* o  : Output vector that contains vector 1 + vector 2  */
    const short N                           /* i  : Vector lenght                                    */
);

void v_sub(
    const float x1[],                       /* i  : Input vector 1                                   */
    const float x2[],                       /* i  : Input vector 2                                   */
    float y[],                        /* o  : Output vector that contains vector 1 - vector 2  */
    const short N                           /* i  : Vector lenght                                    */
);

void v_mult(
    const float x1[],                       /* i  : Input vector 1                                   */
    const float x2[],                       /* i  : Input vector 2                                   */
    float y[],                        /* o  : Output vector that contains vector 1 .* vector 2 */
    const short N                           /* i  : Vector lenght                                    */
);

void v_multc(
    const float x[],                        /* i  : Input vector                                     */
    const float c,                          /* i  : Constant                                         */
    float y[],                        /* o  : Output vector that contains c*x                  */
    const short N                           /* i  : Vector lenght                                    */
);

int squant(                                 /* o: index of the winning codeword   */
    const float x,                          /* i: scalar value to quantize        */
    float *xq,                        /* o: quantized value                 */
    const float cb[],                       /* i: codebook                        */
    const int   cbsize                      /* i: codebook size                   */
);

int vquant(                                 /* o: index of the winning codevector */
    float x[],                        /* i: vector to quantize              */
    const float x_mean[],                   /* i: vector mean to subtract (0 if none) */
    float xq[],                       /* o: quantized vector                */
    const float cb[],                       /* i: codebook                        */
    const int   dim,                        /* i: dimension of codebook vectors   */
    const int   cbsize                      /* i: codebook size                   */
);

int w_vquant(                               /* o: index of the winning codevector */
    float x[],                      /* i: vector to quantize              */
    const float x_mean[],                 /* i: vector mean to subtract (0 if none) */
    const short weights[],                /* i: error weights                   */
    float xq[],                     /* o: quantized vector                */
    const float cb[],                     /* i: codebook                        */
    const int   dim,                      /* i: dimension of codebook vectors   */
    const int   cbsize,                   /* i: codebook size                   */
    const short reverse                   /* i: reverse codebook vectors        */
);

short usquant(                              /* o: index of the winning codeword   */
    const float x,                          /* i: scalar value to quantize        */
    float *xq,                        /* o: quantized value                 */
    const float qlow,                       /* i: lowest codebook entry (index 0) */
    const float delta,                      /* i: quantization step               */
    const short cbsize                      /* i: codebook size                   */
);

float usdequant(                            /* o: dequanzited gain                */
    const int idx,                          /* i: quantizer index                 */
    const float qlow,                       /* i: lowest codebook entry (index 0) */
    const float delta                       /* i: quantization step               */
);

void v_sort(
    float *r,                         /* i/o: Vector to be sorted in place */
    const short lo,                         /* i  : Low limit of sorting range   */
    const short up                          /* i  : High limit of sorting range  */
);

float var(                                  /* o: variance of vector                    */
    const float *x,                          /* i: input vector                          */
    const int len                            /* i: length of inputvector                 */
);

float std_dev(                              /* o: standard deviation                    */
    const float *x,                         /* i: input vector                          */
    const int len                           /* i: length of the input vector            */
);

float dot_product_mat(                      /* o  : the dot product x'*A*x        */
    const float *x,                         /* i  : vector x                      */
    const float *A,                         /* i  : matrix A                      */
    const short  m                          /* i  : vector length                 */
);

float root_a(
    float a
);

float root_a_over_b(
    float a,
    float b
);

void polezero_filter (
    const float *in,                        /* i  : input vector                              */
    float *out,                       /* o  : output vector                             */
    const short N,                          /* i  : input vector size                         */
    const float *b,                         /* i  : numerator coefficients                    */
    const float *a,                         /* i  : denominator coefficients                  */
    const short order,                      /* i  : filter order                              */
    float *mem                        /* i/o: filter memory                             */
);

double rint_new(
    double x                          /* i/o: Round to the nearest integer with mid point exception */
);

double anint(
    double x                          /* i/o: Round to the nearest integer */
);
short is_numeric_float(                     /* o : Output either 1 if Numeric, 0 if NaN or Inf */
    float x                                 /* i : Input value which is checked if numeric or not */
);

void push_indice(
    Encoder_State *st,                      /* i/o: encoder state structure */
    short id,                      /* i  : ID of the indice */
    unsigned short value,                   /* i  : value of the quantized indice */
    short nb_bits                  /* i  : number of bits used to quantize the indice */
);

void push_next_indice(
    Encoder_State *st,                      /* i/o: encoder state structure */
    unsigned short value,                   /* i  : value of the quantized indice */
    short nb_bits                  /* i  : number of bits used to quantize the indice */
);

void push_next_bits(
    Encoder_State *st,                     /* i/o: encoder state structure */
    int bits[],                 /* i  : bit buffer to pack, sequence of single bits */
    short nb_bits                 /* i  : number of bits to pack */
);

unsigned short get_next_indice(             /* o  : value of the indice */
    Decoder_State *st,                        /* i/o: decoder state structure */
    short nb_bits                     /* i  : number of bits that were used to quantize the indice */
);

unsigned short get_next_indice_1(           /* o  : value of the indice */
    Decoder_State *st                         /* i/o: decoder state structure */
);

void get_next_indice_tmp(
    Decoder_State *st,                      /* o  : decoder state structure */
    short nb_bits                           /* i  : number of bits that were used to quantize the indice */
);

unsigned short get_indice(                  /* o  : value of the indice */
    Decoder_State *st,                        /* i/o: decoder state structure */
    short pos,                        /* i  : absolute position in the bitstream */
    short nb_bits                     /* i  : number of bits that were used to quantize the indice */
);

unsigned short get_indice_1(                /* o  : value of the indice */
    Decoder_State *st,                        /* i/o: decoder state structure */
    short pos                         /* i  : absolute position in the bitstream */
);

void reset_indices_enc(
    Encoder_State *st                         /* i/o: encoder state structure */
);

void reset_indices_dec(
    Decoder_State *st                         /* i/o: decoder state structure */
);

void write_indices(
    Encoder_State *st,                        /* i/o: encoder state structure */
    FILE *file                       /* i  : output bitstream file                     */
    , UWord8 *pFrame,                    /* i  : byte array with bit packet and byte aligned coded speech data */
    Word16 pFrame_size                 /* i  : size of the binary encoded access unit [bits] */
);

short read_indices(                         /* o  : 1 = OK, 0 = something wrong            */
    Decoder_State *st,                        /* i/o: decoder state structure */
    FILE *file,                      /* i  : bitstream file                         */
    const short rew_flag                    /* i  : rewind flag (rewind file after reading) */
);

Word16 read_indices_mime(                   /* o  : 1 = reading OK, 0 = problem            */
    Decoder_State *st,                      /* i/o: decoder state structure                */
    FILE *file,                             /* i  : bitstream file                         */
    Word16 rew_flag                         /* i  : rewind flag (rewind file after reading)*/
);

void indices_to_serial(
    const Encoder_State *st,                /* i: encoder state structure */
    UWord8 *pFrame,                         /* o: byte array with bit packet and byte aligned coded speech data */
    Word16 *pFrame_size
);

void indices_to_serial_generic(
    const Indice *ind_list,                /* i: indices list */
    const Word16  num_indices,             /* i: number of indices to write */
    UWord8 *pFrame,                        /* o: byte array with bit packet and byte aligned coded speech data */
    Word16 *pFrame_size                    /* o: size of the binary encoded access unit [bits] */
);

void evs_dec_previewFrame(
    unsigned char *bitstream,               /* i  : bitstream pointer */
    int bitstreamSize,                      /* i  : bitstream size    */
    short *partialCopyFrameType,            /* o  : frame type of the partial copy */
    short *partialCopyOffset                /* o  : offset of the partial copy relative to the primary copy */
);

void read_indices_from_djb(
    Decoder_State *st,                        /* i/o: decoder state structure                */
    unsigned char *pt_stream,                 /* i  : bitstream file                         */
    int nbits                       /* i  : number of bits                         */
    ,short partialframe              /* i  : partial frame information     */
    ,short next_coder_type           /* i  : next coder type information     */
);

void getPartialCopyInfo(
    Decoder_State *st,                        /* i  : decoder state structure       */
    short *coder_type,
    short *sharpFlag
);

void get_NextCoderType(
    unsigned char *bitsteam,                  /* i  : bitstream            */
    short *next_coder_type            /* o  : next coder type      */
);

int print_disclaimer (
    FILE *fPtr
);

void autocorr(
    const float *x,                         /* i  : input signal               */
    float *r,                         /* o  : autocorrelations vector    */
    const short m,                          /* i  : order of LP filter         */
    const short len,                        /* i  : window size                */
    const float *wind,                      /* i  : window                     */
    const short rev_flag,                   /* i  : flag to reverse window     */
    const short sym_flag,                   /* i  : symmetric window flag      */
    const short no_thr                      /* i  : flag to avoid thresholding */
);

short lev_dur(                              /* o:   energy of prediction error   */
    float *a,                         /* o:   LP coefficients (a[0] = 1.0) */
    const float *r,                         /* i:   vector of autocorrelations   */
    const short m,                          /* i:   order of LP filter           */
    float epsP[]                      /* o:   prediction error energy      */
);

float get_delay(                            /* o  : delay value in ms                         */
    const short what_delay,                 /* i  : what delay? (ENC or DEC)                  */
    const int   io_fs                       /* i  : input/output sampling frequency           */
);

void decision_matrix_enc(
    Encoder_State *st,                        /* i/o: encoder state structure                   */
    const short sp_aud_decision1,           /* i  : 1st stage speech/music classification     */
    const short sp_aud_decision2,           /* i  : 2nd stage speech/music classification     */
    const short coder_type,                 /* i  : coder type                                */
    const short vad_flag,
    short *hq_core_type               /* o  : HQ core type                              */
);

void signalling_enc(
    Encoder_State *st,                        /* i  : encoder state structure                   */
    const short coder_type,                 /* i  : coder type                                */
    const short sharpFlag                   /* i  : formant sharpening flag                   */
);

short signalling_mode1_tcx20_enc(
    Encoder_State *st,                        /* i/o: encoder state structure                   */
    const short push                        /* i  : flag to push indice                       */
);

void decision_matrix_dec(
    Decoder_State *st,                        /* i/o: decoder state structure                   */
    short *coder_type,                /* o  : coder type                                */
    short *sharpFlag,                 /* o  : formant sharpening flag                   */
    short *hq_core_type,              /* o  : HQ core type                              */
    short *core_switching_flag        /* o  : ACELP->HQ switching frame flag            */
);

float lsf_stab(                             /* o  : LP filter stability                         */
    const float *lsf,                       /* i  : LSF vector                                  */
    const float *lsfold,                    /* i  : old LSF vector                              */
    const short Opt_AMR_WB,                 /* i  : flag indicating AMR-WB IO mode              */
    const short L_frame                     /* i  : frame lenght                                */
);

void hf_synth_amr_wb_init(
    float *prev_r,                    /* o  : 1 sample memory for preemphasis/deemphasis*/
    float *fmerit_w_sm,               /* o  : 1 sample memory fmerit_w param*/
    float mem_syn_hf[],               /* o  : HF LPC synthesis filter initialization           */
    short *frame_count,               /* o  : frame counter initialization                     */
    float *ne_min,                    /* o  : minimum Noise gate - short-term energy initialization*/
    float *fmerit_m_sm,               /* o  : 1 sample memory fmerit_m param                   */
    float *voice_fac,                 /* o  : voice factor initialization                      */
    float *unvoicing,                 /* o  : unvoiced parameter                               */
    float *unvoicing_sm,              /* o  : smoothed unvoiced parameter                      */
    short *unvoicing_flag,            /* o  : unvoiced flag                                    */
    short *voicing_flag,              /* o  : voiced flag                                      */
    short *start_band_old,            /* o  : previous start point for copying frequency band  */
    float *OptCrit_old                /* o  : previous criterion value for deciding the start point */
);

void hf_synth_amr_wb_reset(
    short *seed2,                     /* i/o: random seed for HF noise gen                     */
    float mem_syn_hf[],               /* o  : HF synthesis memory                              */
    float mem_hp_interp[],            /* o  : interpol. memory                                 */
    float *prev_r,                    /* o  : 1 sample memory for deemphasis                   */
    float *fmerit_w_sm,               /* o  : 1 sample memory fmerit_w param                   */
    float delay_syn_hf[],             /* o  : HF synthesis memory                               */
    short *frame_count,               /* o  : frame counter memory                             */
    float *ne_min,                    /* o  : minimum Noise gate - short-term energy memory    */
    float *fmerit_m_sm,               /* o  : 1 sample memory fmerit_m param                   */
    float *voice_fac,                 /* o  : voice factor memory                              */
    float *unvoicing,                 /* o  : unvoiced parameter                               */
    float *unvoicing_sm,              /* o  : smoothed unvoiced parameter                      */
    short *unvoicing_flag,            /* o  : unvoiced flag                                    */
    short *voicing_flag,              /* o  : voiced flag                                      */
    short *start_band_old,            /* o  : previous start point for copying frequency band  */
    float *OptCrit_old                /* o  : previous criterion value for deciding the start point */
);

void hf_synth_amr_wb(
    const long  core_brate,                 /* i  : core bitrate                                      */
    const short output_frame,               /* i  : output frame length                               */
    const float *Aq,                        /* i  : quantized Az                                      */
    const float *exc,                       /* i  : excitation at 12.8 kHz                            */
    float *synth,                     /* i/o: synthesis signal at 12.8 kHz                      */
    float *mem_syn_hf,                /* i/o: HF synthesis memory                               */
    float *delay_syn_hf,              /* i/o: HF synthesis memory                               */
    float *prev_r,                    /* i/o: preemphasis/deemphasis filter memory              */
    float *fmerit_w_sm,               /* i/o: smoothed fmerit_w                                 */
    short *amr_io_class,              /* i  : signal class (determined by FEC algorithm)        */
    float *mem_hp_interp,             /* i/o: interpol. memory                                  */
    float *synth_out,                 /* i/o: synthesis signal at output Fs                     */
    float fmerit,                     /* i  : classify parameter from FEC                       */
    const short *hf_gain,                   /* i  : decoded HF gain                                   */
    const float *voice_factors,             /* i  : voicing factors                                   */
    const float pitch_buf[],                /* i  : pitch buffer                                      */
    const float ng_ener_ST,                 /* i  : Noise gate - short-term energy                    */
    const float *lsf_new,                   /* i  : ISF vector                                        */
    short *frame_count,               /* i/o: frame counter                                     */
    float *ne_min,                    /* i/o: minimum Noise gate                                */
    float *fmerit_m_sm,               /* i/o: smoothed fmerit_m                                 */
    float *voice_facor_sm,            /* i/o: voice factor memory                               */
    float *unvoicing,                 /* i/o: unvoiced parameter                                */
    float *unvoicing_sm,              /* i/o: smoothed unvoiced parameter                       */
    short *unvoicing_flag,            /* i/o: unvoiced flag                                     */
    short *voicing_flag,              /* i/o: voiced flag                                       */
    short *start_band_old,            /* i/o: previous start point for copying frequency band   */
    float *OptCrit_old                /* i/o: previous criterion value for deciding the start point */
);

void hf_cod_init(
    float *mem_hp400_enc,             /* o: memory of hp 400 Hz filter   */
    float *mem_hf1_enc,               /* o: HF band-pass filter memory   */
    float *mem_syn_hf_enc,            /* o: HF synthesis memory          */
    float *mem_hf2_enc,               /* o: HF band-pass filter memory   */
    float *gain_alpha                 /* o: smoothing gain for transitions between active and inactive frames */
);

void hf_cod(
    const long  core_brate,                 /* i  : core bitrate                 */
    const float *speech16k,                 /* i  : original speech at 16 kHz    */
    const float Aq[],                       /* i  : quantized Aq                 */
    const float exc[],                      /* i  : excitation at 12.8 kHz       */
    float synth[],                    /* i  : 12.8kHz synthesis signal     */
    short *seed2_enc,                 /* i/o: random seed for HF noise gen */
    float *mem_hp400_enc,             /* i/o: memory of hp 400 Hz filter   */
    float *mem_syn_hf_enc,            /* i/o: HF synthesis memory          */
    float *mem_hf1_enc,               /* i/o: HF band-pass filter memory   */
    float *mem_hf2_enc,               /* i/o: HF band-pass filter memory   */
    const short *dtxHangoverCount,
    float *gain_alpha,                /* i/o: smoothing gain for transitions between active and inactive frames */
    short *hf_gain                    /* o  :  HF gain to be transmitted to decoder */
);

void hf_synth_init(
    float mem_hp400[],                /* o:   400 Hz high pass filter memory initialization     */
    float mem_hf[]                    /* o:   band pass 6kHz to 7kHz FIR filter initialization  */
);

void hf_synth_reset(
    short *seed2,                     /* i/o: random seed for HF noise gen                */
    float mem_hf[],                   /* o  : HF band-pass filter memory                  */
    float mem_syn_hf[],               /* o  : HF synthesis memory                         */
    float mem_hp400[],                /* o  : memory of hp 400 Hz filter                  */
    float mem_hp_interp[],            /* o  : interpol. memory                            */
    float delay_syn_hf[]              /* o  : HF synthesis memory                         */
);

void hf_synth(
    const long  core_brate,                 /* i  : core bitrate                 */
    const short output_frame,               /* i  : output frame length          */
    const float *Aq,                        /* i  : quantized Az                 */
    const float *exc,                       /* i  : excitation at 12.8 kHz       */
    float *synth,                     /* i/o: 12.8kHz synthesis signal     */
    float *synth16k,                  /* i/o: 16kHz synthesis signal       */
    short *seed2,                     /* i/o: random seed for HF noise gen */
    float *mem_hp400,                 /* i/o: memory of hp 400 Hz filter   */
    float *mem_syn_hf,                /* i/o: HF synthesis memory          */
    float *mem_hf,                    /* i/o: HF band-pass filter memory   */
    float *delay_syn_hf,              /* i/o: HF synthesis memory          */
    float *mem_hp_interp              /* i/o: interpol. memory             */
);

short lsp_convert_poly(
    float  w[],                       /* i/o: LSP or ISP parameters          */
    const short L_frame,                    /* i  : flag for up or down conversion */
    const short Opt_AMRWB                   /* i  : flag for the AMR-WB IO mode    */
);

short findpulse(                            /* o  : pulse position        */
    const short L_frame,                    /* i  : length of the frame   */
    const float res[],                      /* i  : residual signal       */
    const short T0,                         /* i  : integer pitch         */
    const short enc_dec,                    /* i  : flag enc/dec, 0 - enc, 1 - dec          */
    short *sign                       /* i/o: sign of the maximum */
);

void fft_rel(
    float x[],                        /* i/o: input/output vector    */
    const short n,                          /* i  : vector length          */
    const short m                           /* i  : log2 of vector length  */
);

void ifft_rel(
    float io[],                       /* i/o: input/output vector   */
    const short n,                          /* i  : vector length         */
    const short m                           /* i  : log2 of vector length */
);

void preemph(
    float *signal,                    /* i/o: signal             */
    const float mu,                         /* i  : preemphasis factor */
    const short L,                          /* i  : vector size        */
    float *mem                        /* i/o: memory (x[-1])     */
);

void cb_shape(
    const short preemphFlag,                /* i  : flag for pre-emphasis                       */
    const short pitchFlag,                  /* i  : flag for pitch sharpening                   */
    const short scramblingFlag,             /* i  : flag for phase scrambling                   */
    const short formantFlag,                /* i  : flag for formant sharpening                 */
    const short formantTiltFlag,            /* i  : flag for formant tilt                       */
    const float g1,                         /* i  : formant sharpening numerator weighting      */
    const float g2,                         /* i  : formant sharpening denominator weighting    */
    const float *p_Aq,                      /* i  : LP filter coefficients                      */
    float *code,                      /* i/o: signal to shape                             */
    const float tilt_code,                  /* i  : tilt of code                                */
    const float pt_pitch                    /* i  : pointer to current subframe fractional pitch*/
);

void isp2a(
    const float *isp,                       /* i  : ISP vector (in the cosine domain)       */
    float *a,                         /* o  : LP filter coefficients                  */
    const short m                           /* i  : order of LP analysis                    */
);

void isp2isf(
    const float isp[],                      /* i  : isp[m] (range: -1<=val<1)               */
    float isf[],                      /* o  : isf[m] normalized (range: 0<=val<=fs/2) */
    const short m,                          /* i  : LPC order                               */
    const float fs                          /* i  : sampling frequency                      */
);

void isf2isp(
    const float isf[],                      /* i  : isf[m] normalized (range: 0<=val<=fs/2) */
    float isp[],                      /* o  : isp[m] (range: -1<=val<1)               */
    const short m,                          /* i  : LPC order                               */
    const float fs                          /* i  : sampling frequency                      */
);

void reorder_isf(
    float *isf,                       /* i/o: vector of isfs in the frequency domain (0..0.5)*/
    const float min_dist,                   /* i  : minimum required distance               */
    const short n,                          /* i  : LPC order                               */
    const float fs                          /* i  : sampling frequency                      */
);

void lsp2lsf(
    const float lsp[],                      /* i  : isp[m] (range: -1<=val<1)               */
    float lsf[],                      /* o  : isf[m] normalized (range: 0<=val<=fs/2) */
    const short m,                          /* i  : LPC order                               */
    const float fs                          /* i  : sampling frequency                      */
);

void lsf2lsp(
    const float lsf[],                      /* i  : isf[m] normalized (range: 0<=val<=fs/2) */
    float lsp[],                            /* o  : isp[m] (range: -1<=val<1)               */
    const short m,                          /* i  : LPC order                               */
    const float fs                          /* i  : sampling frequency                      */
);

void lsp2isp(
    const float *lsp,                       /* i  : LSP vector                              */
    float *isp,                       /* o  : ISP filter coefficients                 */
    float *stable_isp,                /* i/o: ISP filter coefficients                 */
    const short m                           /* i  : order of LP analysis                    */
);

void isp2lsp(
    const float *isp,                       /* i  : LSP vector                              */
    float *lsp,                       /* o  : ISP filter coefficients                 */
    float *stable_lsp,                /* i/o: stable LSP filter coefficients          */
    const short m                           /* i  : order of LP analysis                    */
);

void reorder_lsf(
    float *lsf,                             /* i/o: vector of lsfs in the frequency domain (0..0.5)*/
    const float min_dist,                   /* i  : minimum required distance               */
    const short n,                          /* i  : LPC order                               */
    const float fs                          /* i  : sampling frequency                      */
);

void CNG_exc(
    const long  core_brate,                 /* i  : core bitrate                            */
    const short L_frame,                    /* i  : length of the frame                     */
    float *Enew,                      /* i/o: decoded SID energy                      */
    short *seed,                      /* i/o: random generator seed                   */
    float exc[],                      /* o  : current non-enhanced excitation         */
    float exc2[],                     /* o  : current enhanced excitation             */
    float *lp_ener,                   /* i/o: LP filtered E                           */
    const long  last_core_brate,            /* i  : previous frame core bitrate             */
    short *first_CNG,                 /* i/o: first CNG frame flag for energy init.   */
    short *cng_ener_seed,             /* i/o: random generator seed for CNG energy    */
    float bwe_exc[],                  /* o  : excitation for SWB TBE                  */
    const short allow_cn_step,              /* i  : allow CN step                           */
    short *last_allow_cn_step,        /* i/o: last CN_step                            */
    const short num_ho,                     /* i  : number of selected hangover frames      */
    float q_env[],
    float *lp_env,
    float *old_env,
    float *exc_mem,
    float *exc_mem1,
    short *sid_bw,
    short *cng_ener_seed1,
    float exc3[],
    short Opt_AMR_WB
);

void cng_params_upd(
    const float lsp_new[],                  /* i  : LSP aprameters                                      */
    const float exc2[],                     /* i  : current enhanced excitation                         */
    const short L_frame,                    /* i  : frame length                                        */
    short *ho_circ_ptr,               /* i/o: pointer for CNG averaging buffers                   */
    float ho_ener_circ[],             /* o  : energy buffer for CNG averaging                     */
    short *ho_circ_size,              /* i/o: size of DTX hangover history buffer for averaging   */
    float ho_lsp_circ[],              /* o  : old LSP buffer for CNG averaging                    */
    const short enc_dec_flag,               /* i  : Flag indicating encoder or decoder (ENC,DEC)        */
    float ho_env_circ[],              /* i/o: Envelope buffer                                     */
    short *cng_buf_cnt,               /* i/o: Counter of postponed FFT-processing instances       */
    float cng_exc2_buf[],             /* i/o: Excitation buffer                                   */
    long  cng_brate_buf[],            /* i/o: last_active_brate buffer                            */
    const long  last_active_brate           /* i  : Last active bit rate                                */
);

void cng_params_postupd(
    const short ho_circ_ptr,                /* i  : pointer for CNG averaging buffers   */
    short *cng_buf_cnt,               /* i/o: counter for CNG store buffers       */
    const float *const cng_exc2_buf,        /* i  : Excitation buffer                   */
    const long  *const cng_brate_buf,       /* i  : bit rate buffer                     */
    float ho_env_circ[]               /* i/o: Envelope buffer                     */
);

void disf_ns_28b(
    short *indice,                    /* i  : quantized indices (use indice[0] = -1 in the decoder) */
    float *isf_q                      /* o  : ISF in the frequency domain (0..6400) */
);

void limit_T0(
    const short L_frame,                    /* i  : length of the frame                                  */
    const short delta,                      /* i  : Half the close-loop searched interval                */
    const short pit_flag,                   /* i  : selecting absolute(0) or delta(1) pitch quantization */
    const short limit_flag,                 /* i  : flag for limits (0=restrained, 1=extended)           */
    const short T0,                         /* i  : rough pitch estimate around which the search is done */
    const short T0_frac,                    /* i  : pitch estimate fractional part                       */
    short *T0_min,                    /* o  : lower pitch limit                                    */
    short *T0_max                     /* o  : higher pitch limit                                   */
);

float interpolation(                        /* o  : interpolated value              */
    const float *x,                         /* i  : input vector                    */
    const float *win,                       /* i  : interpolation window            */
    const short frac,                       /* i  : fraction                        */
    const short up_samp,                    /* i  : upsampling factor               */
    const short nb_coef                     /* i  : nb of filter coef               */
);

void deemph(
    float *signal,                    /* i/o: signal                          */
    const float mu,                         /* i  : deemphasis factor               */
    const short L,                          /* i  : vector size                     */
    float *mem                        /* i/o: memory (y[-1])                  */
);

float est_tilt(                             /* o  : tilt of the code                */
    const float *adpt_exc,                  /* i  : adaptive excitation vector      */
    const float gain_pit,                   /* i  : adaptive gain                   */
    const float *fixe_exc,                  /* i  : algebraic exctitation vector    */
    const float gain_code,                  /* i  : algebraic code gain             */
    float *voice_fac,                 /* o  : voicing factor                  */
    const short L_subfr,                    /* i : subframe size                    */
    const short flag_tilt                   /* i : flag for special tilt            */
);

void weight_a(
    const float *a,                         /* i  : LP filter coefficients          */
    float *ap,                        /* o  : weighted LP filter coefficients */
    const float gamma,                      /* i  : weighting factor                */
    const short m                           /* i  : order of LP filter              */
);

void weight_a_subfr(
    const short nb_subfr,               /* i  : number of subframes             */
    const float *a,                     /* i  : LP filter coefficients          */
    float *ap,                    /* o  : weighted LP filter coefficients */
    const float gamma,                  /* i  : weighting factor                */
    const short m                       /* i  : order of LP filter              */
);

void syn_12k8(
    const short L_frame,                    /* i  : length of the frame                         */
    const float *Aq,                        /* i  : LP filter coefficients                      */
    const float *exc,                       /* i  : input signal                                */
    float *synth,                     /* o  : output signal                               */
    float *mem,                       /* i/o: initial filter states                       */
    const short update_m                    /* i  : update memory flag: 0 --> no memory update  */
);                                          /*                          1 --> update of memory  */

void syn_filt(
    const float a[],                        /* i  : LP filter coefficients                     */
    const short m,                          /* i  : order of LP filter                         */
    const float x[],                        /* i  : input signal                               */
    float y[],                        /* o  : output signal                              */
    const short l,                          /* i  : size of filtering                          */
    float mem[],                      /* i/o: initial filter states                      */
    const short update_m                    /* i  : update memory flag: 0 --> no memory update */
);                                          /*                          1 --> update of memory */
void synth_mem_updt2(
    const short L_frame,                    /* i  : frame length                                */
    const short last_L_frame,               /* i  : frame length                                */
    float old_exc[],                        /* i/o  : excitation buffer                           */
    float mem_syn_r[],                      /* i/o: synthesis filter memory                     */
    float mem_syn2[],                       /* o  : synthesis filter memory for find_target     */
    float mem_syn[],                        /* o  : synthesis filter memory for find_target     */
    const short dec                         /* i: flag for decoder indication           */
);

void int_lsp(
    const short L_frame,                    /* i  : length of the frame               */
    const float lsp_old[],                  /* i  : LSPs from past frame              */
    const float lsp_new[],                  /* i  : LSPs from present frame           */
    float *Aq,                        /* o  : LP coefficients in both subframes */
    const short m,                          /* i  : order of LP filter                */
    const float *int_coeffs,                /* i  : interpolation coefficients        */
    const short Opt_AMR_WB                  /* i  : flag indicating AMR-WB IO mode    */
);

void int_lsp4(
    const short L_frame,                    /* i  : length of the frame               */
    const float lsp_old[],                  /* i  : previous end-frame LSPs           */
    const float lsp_mid[],                  /* i  : current mid-frame LSPs            */
    const float lsp_new[],                  /* i  : current end-frame LSPs            */
    float *Aq,                        /* o  : LP coefficients in both subframes */
    const short m,                          /* i  : order of LP filter                */
    short relax_prev_lsf_interp       /* i  : relax prev frame lsf interp after erasure */
);

short modify_Fs(                            /* o  : length of output               */
    const float sigIn[],                    /* i  : signal to decimate             */
    short lg,                         /* i  : length of input                */
    const int   fin,                        /* i  : frequency of input             */
    float sigOut[],                   /* o  : decimated signal               */
    const int   fout,                       /* i  : frequency of output            */
    float mem[]                       /* i/o: filter memory                  */
    ,int   nblp                        /* i  : flag indicating if NB low-pass is applied */
);

void pred_lt4(
    const float excI[],                     /* i  : input excitation buffer     */
    float excO[],                     /* o  : output excitation buffer    */
    const short T0,                         /* i:   integer pitch lag           */
    short       frac,                       /* i:   fraction of lag             */
    const short L_subfr,                    /* i:   subframe size               */
    const float *win,                       /* i:   interpolation window        */
    const short nb_coef,                    /* i  : nb of filter coef           */
    const short up_sample                   /* i  : up_sample                   */
);

void pred_lt4_tc(
    float       exc[],                      /* i:   excitation buffer            */
    const short T0,                         /* i:   integer pitch lag            */
    short frac,                       /* i:   fraction of lag              */
    const float *win,                       /* i  : interpolation window         */
    const short imp_pos,                    /* i:   glottal impulse position     */
    const short i_subfr                     /* i:   subframe index               */
);

void residu(
    const float *a,                         /* i  : LP filter coefficients                  */
    const short m,                          /* i  : order of LP filter                      */
    const float *x,                         /* i  : input signal (usually speech)           */
    float *y,                         /* o  : output signal (usually residual)        */
    const short l                           /* i  : size of filtering                       */
);

void calc_residu(
    const float *speech,                    /* i  : weighted speech signal                  */
    float *res,                       /* o  : residual signal                         */
    const float *p_Aq,                      /* i  : quantized LP filter coefficients        */
    const short L_frame                     /* i  : size of frame                           */
);

float enr_1_Az(                             /* o  : impulse response energy                 */
    const float Aq[],                       /* i  : LP filter coefs                         */
    const short len                         /* i  : impulse response length                 */
);


void Es_pred_enc(
    float *Es_pred,                   /* o  : predicited scaled innovation energy     */
    int *Es_pred_indice,              /* o  : indice corresponding to above parameter */
    const short L_frame,                    /* i  : length of the frame                     */
    const short L_subfr,                    /* i  : length of the subframe                  */
    const float *res,                       /* i  : residual signal                         */
    const float *voicing,                   /* i  : normal. correlattion in three 1/2frames */
    const short nb_bits,                    /* i  : allocated number of bits                */
    const short no_ltp                      /* i  : no_ltp flag                             */
);

void init_lvq(
    unsigned int offset_scale1[][MAX_NO_SCALES+1],
    unsigned int offset_scale2[][MAX_NO_SCALES+1],
    unsigned int offset_scale1_p[][MAX_NO_SCALES+1],
    unsigned int offset_scale2_p[][MAX_NO_SCALES+1],
    short no_scales[][2],
    short no_scales_p[][2]
);

float mslvq (
    float *pTmp,                            /* i  : M-dimensional input vector */
    float *quant,                           /* o  : quantized vector */
    float *cv_out,                          /* o  : corresponding 8-dim lattice codevectors (without the scaling) */
    int   *idx_lead,                        /* o  : leader index for each 8-dim subvector  */
    int   *idx_scale,                       /* o  : scale index for each subvector */
    float *w,                               /* i  : weights for LSF quantization */
    short mode,                             /* i  : number indicating the coding type (V/UV/G...)*/
    short mode_glb,                         /* i  : LVQ coding mode */
    int   pred_flag,                        /* i  : prediction flag (0: safety net, 1 - predictive )*/
    short  no_scales[][2]
);

void permute(
    float *pTmp1,                           /* (i/o): vector whose components are to be permuted */
    const short *perm                             /* (i)  : permutation info (indexes that should be interchanged), max two perms */
);

float mslvq_cng(
    short idx_cv,                           /* (i): index of cv from previous stage */
    float *pTmp,                            /* (i): 16 dimensional input vector */
    float *quant,                           /* (o): quantized vector */
    float *cv_out,                          /* (o): corresponding 8-dim lattice codevectors (without the scaling) */
    int   *idx_lead,                        /* (o): leader index for each 8-dim subvector  */
    int   *idx_scale,                       /* (o): scale index for each subvector */
    const float *w,                               /* (i): weights for LSF quantization */
    short * no_scales
);

short deindex_lvq_cng(
    short *index,                           /* i  : index to be decoded, as an array of 3 short */
    float *x_lvq,                           /* o  : decoded codevector */
    short idx_cv,                           /* i  : relative mode_lvq, wrt START_CNG */
    int no_bits,                            /* i  : number of bits for lattice */
    unsigned int * p_offset_scale1,
    unsigned int * p_offset_scale2,
    short * p_no_scales
);

void multiply32_32_64(
    unsigned int x,                         /* (i): operand 1 */
    unsigned int y,                         /* (i): operand 2 */
    unsigned int *res                       /* (o): result as array of two uint32 */
);

short deindex_lvq(
    short *index,                           /* (i): index to be decoded, as an array of 3 short */
    float *x_lvq,                           /* (o): decoded codevector */
    short mode,                             /* (i): LVQ  coding mode (select scales & no_lead ), or idx_cv */
    short sf_flag,                          /* (i): safety net flag  */
    short no_bits,                          /* (i): number of bits for lattice */
    unsigned int *p_offset_scale1,          /* i  : offset for first subvector */
    unsigned int *p_offset_scale2,          /* i  : offset for the second subvector */
    short * p_no_scales
);
short vq_dec_lvq (
    short sf_flag,                          /* i  : safety net flag                           */
    float x[],                              /* o  : Decoded vector                            */
    short indices[],                        /* i  : Indices                                   */
    short stages,                           /* i  : Number of stages                          */
    short N,                                /* i  : Vector dimension                          */
    short mode,                             /* i  : lvq coding type                           */
    short no_bits,                          /* i  : no. bits for lattice                      */
    unsigned int * p_offset_scale1,
    unsigned int * p_offset_scale2,
    unsigned int * p_offset_scale1_p,
    unsigned int * p_offset_scale2_p,
    short * p_no_scales,
    short * p_no_scales_p
);

void index_lvq (
    float *quant,                           /* i  : codevector to be indexed (2 8-dim subvectors)*/
    int   *idx_lead,                        /* i  : leader class index for each subvector */
    int   *idx_scale,                       /* i  : scale index for each subvector */
    int   mode,                             /* i  : integer signalling the quantizer structure for the current bitrate */
    short *index,                           /* o  : encoded index (represented on 3 short each with 15 bits ) */
    unsigned int * p_offset_scale1,
    unsigned int * p_offset_scale2,
    short * p_no_scales
);

short qlsf_ARSN_tcvq_Dec_16k (
    float *y,                         /* o  : Quantized LSF vector    */
    short *indice,                    /* i  : Indices                 */
    const short nBits                       /* i  : number of bits          */
);

int lsf_bctcvq_encprm(
    Encoder_State *st,
    int *param_lpc,
    short *bits_param_lpc,
    short no_indices
);

int lsf_bctcvq_decprm(
    Decoder_State *st,
    int *param_lpc
);

void lsf_allocate(
    const short nBits,                      /* i  : Number of bits to use for quantization      */
    const short framemode,                  /* i  : ISF quantizer mode                          */
    const short framemode_p,                /* i  : ISF quantizer mode                          */
    short *stages0,                   /* o  : Number of stages for safety-net quantizer   */
    short *stages1,                   /* o  : Number of stages for predictive quantizer   */
    short levels0[],                  /* o  : Number of vectors for each stage for SFNET  */
    short levels1[],                  /* o  : Number of vectors for each stage for pred   */
    short bits0[],                    /* o  : Number of bits for each stage safety net    */
    short bits1[]                     /* o  : Number of bits for each stage predictive    */
);

void disf_2s_36b(
    short *indice,                    /* i  : quantized indices (use indice[0] = -1 in the decoder) */
    float *isf_q,                     /* o  : quantized ISFs in the cosine domain */
    float *mem_AR,                    /* i/o: quantizer memory for AR model       */
    float *mem_MA                     /* i/o: quantizer memory for MA model       */
);

void disf_2s_46b(
    short *indice,                    /* i  : quantized indices (use indice[0] = -1 in the decoder) */
    float *isf_q,                     /* o  : quantized ISFs in the cosine domain */
    float *mem_AR,                    /* o  : quantizer memory for AR model       */
    float *mem_MA                     /* i/o: quantizer memory for MA model       */
);

void re8_k2y(
    const int *k,                           /* i  : Voronoi index k[0..7]                                   */
    const int m,                            /* i  : Voronoi modulo (m = 2^r = 1<<r, where r is integer >=2) */
    int *y                            /* o  : 8-dimensional point y[0..7] in RE8                      */
);

void re8_PPV(
    float x[],                        /* i  : point in R^8                                */
    int   y[]                         /* o  : point in RE8 (8-dimensional integer vector) */
);

void enhancer(
    const short codec_mode,                 /* i  : flag indicating Codec Mode              */
    const long  core_brate,                 /* i  : core bitrate                            */
    const short cbk_index,                  /* i  :                                         */
    const short Opt_AMR_WB,                 /* i  : flag indicating AMR-WB IO mode          */
    const short coder_type,                 /* i  : coding type                             */
    const short L_frame,                    /* i  : frame size                              */
    const float voice_fac,                  /* i  : subframe voicing estimation             */
    const float stab_fac,                   /* i  : LP filter stablility measure            */
    const float norm_gain_code,             /* i  : normalized innovative cb. gain          */
    const float gain_inov,                  /* i  : gain of the unscaled innovation         */
    float *gc_threshold,              /* i/o: code threshold                          */
    float *code,                      /* i/o: innovation                              */
    float *exc2,                      /* i/o: adapt. excitation/total exc.            */
    const float gain_pit,                   /* i  : Quantized pitch gain                    */
    float *dispMem                    /* i/o: Phase dispersion algorithm memory       */
);

void phase_dispersion(
    const float gain_code,                  /* i  : gain of code                            */
    const float gain_pit,                   /* i  : gain of pitch                           */
    float code[],                     /* i/o: code vector                             */
    const short mode,                       /* i  : level, 0=hi, 1=lo, 2=off                */
    float disp_mem[]                  /* i/o: static memory (size = 8)                */
);

void re8_vor(
    int y[],                          /* i  : point in RE8 (8-dimensional integer vector)      */
    int *n,                           /* o  : codebook number n=0,2,3,4,... (scalar integer)   */
    int k[],                          /* o  : Voronoi index (integer vector of dimension 8) used only if n>4 */
    int c[],                          /* o  : codevector in Q0, Q2, Q3, or Q4 if n<=4, y=c     */
    int *ka                           /* o  : identifier of absolute leader (needed to index c)*/
);

void edct(
    const float *x,                         /* i  : input signal        */
    float *y,                         /* o  : output transform    */
    short length                      /* i  : length              */
);

void edst(
    const float *x,                         /* i  : input signal        */
    float *y,                         /* o  : output transform    */
    short length                      /* i  : length              */
);

void iedct_short(
    const float *in,                        /* i  : input vector        */
    float *out,                       /* o  : output vector       */
    const short segment_length              /* i  : length              */
);


void DoRTFT480(
    float *x,                         /* i/o : real part of input and output data       */
    float *y                          /* i/o : imaginary part of input and output data  */
);

void DoRTFT320(
    float *x,                         /* i/o : real part of input and output data       */
    float *y                          /* i/o : imaginary part of input and output data  */
);

void DoRTFT160(
    float *x,                         /* i/o : real part of input and output data       */
    float *y                          /* i/o : imaginary part of input and output data  */
);

void DoRTFT128(
    float *x,                         /* i/o : real part of input and output data       */
    float *y                          /* i/o : imaginary part of input and output data  */
);

void DoRTFT120(
    float *x,                         /* i/o : real part of input and output data       */
    float *y                          /* i/o : imaginary part of input and output data  */
);

void DoRTFT80(
    float *x,                         /* i/o : real part of input and output data       */
    float *y                          /* i/o : imaginary part of input and output data  */
);

void DoRTFT20(
    float *x,                         /* i/o : real part of input and output data       */
    float *y                          /* i/o : imaginary part of input and output data  */
);

void DoRTFT40(
    float *x,                         /* i/o : real part of input and output data       */
    float *y                          /* i/o : imaginary part of input and output data  */
);

void DoRTFTn(
    float *x,                               /* i/o : real part of input and output data       */
    float *y,                               /* i/o : imaginary part of input and output data  */
    const short n                           /* i   : size of the FFT n=(2^k) up to 1024 */
);

void BASOP_cfft(
    Word32 *re,                             /* i/o: real part                                 */
    Word32 *im,                             /* i/o: imag part                                 */
    Word16 s,                               /* i  : stride real and imag part                 */
    Word16 *scale                           /* i  : scalefactor                               */
);

void sinq(
    const float tmp,                        /* i  : sinus factor cos(tmp*i+phi)  */
    const float phi,                        /* i  : sinus phase cos(tmp*i+phi)  */
    const short N,                          /* i  : size of output */
    float x[]                               /* o  : output vector  */
);

void edct2(
    short n,
    short isgn,
    float *in,
    float *a,
    const short *ip,
    const float *w
);

void stat_noise_uv_mod(
    const short coder_type,                 /* i  : coding type                          */
    float noisiness,                  /* i  : noisiness parameter                  */
    const float *isp_old,                   /* i  : old ISP vector at 4th sfr            */
    const float *isp_new,                   /* i  : ISP vector at 4th sfr                */
    const float *isp_mid,                   /* i  : ISP vector at 2nd sfr                */
    float *Aq,                        /* o  : A(z) quantized for the 4 subframes   */
    float *exc2,                      /* o  : excitation buffer                    */
    const short bfi,                        /* i  : bad frame indicator                  */
    float *ge_sm,                     /* i/o: ????                                 */
    short *uv_count,                  /* i/o: unvoiced counter                     */
    short *act_count,                 /* i/o: activation counter                   */
    float lspold_s[],                 /* i/o: old LSP                              */
    short *noimix_seed,               /* i/o: mixture seed                         */
    float *st_min_alpha,              /* i/o: minimum alpha                        */
    float *exc_pe,                    /* i/o: scale Q_stat_noise                   */
    const long  bitrate,                    /* i  : core bitrate                         */
    const short bwidth                      /* i  : Bandwidth                            */
);

void SWB_TBE_describe_envelope(             /* Function to represent energies in MDCT frequency bands */
    const float *yos,                       /* i  : MDCT coefficients of weighted original */
    float *SWB_env_gain               /* i/o: energy of SWB envelope                 */
);

void dct2(
    const float in[],                       /* i  : time domain input       */
    float out[]                       /* o  : transform domain output */
);

void idct2(
    const float in[],                       /* i  : transform domain input  */
    float out[]                       /* o  : time domain output      */
);

void pre_echo_att(
    float *Last_frame_ener,           /* i/o: Energy of the last frame         */
    float *exc,                       /* i/o: Excitation of the current frame  */
    const short attack_flag,                /* i  : flag signalling attack encoded by AC mode (GSC) */
    const short last_coder_type             /* i  : Last coder type                  */
);

void hq_swb_harmonic_calc_norm_envelop(
    float *SWB_signal,                      /* i  : input signal                */
    float *envelope,                        /* o  : output envelope             */
    int   L_swb_norm,                       /* i  : length of normaliztion      */
    int   SWB_flength                       /* i  : length of input signal      */
);

void limit_band_noise_level_calc(
    short *wnorm,                     /* i  : reordered norm of sub-vectors           */
    short *limit,                     /* o  : highest band of bit allocation          */
    long core_brate,                 /* o  : bit rate                                */
    float *noise_level                /* o  : noise level                             */
);

short peak_avrg_ratio(
    const long total_brate,
    const float *input_hi,                  /* i  : input signal               */
    const short N,                          /* i  : number of coefficients     */
    short *mode_count,                /* i/o: HQ_HARMONIC mode count     */
    short *mode_count1                /* i/o: HQ_NORMAL mode count       */
);

short build_nf_codebook(                    /* o  : Number of coefficients in nf codebook   */
    const short flag_32K_env_ho,            /* i  : Envelope attenuation hangover flag      */
    const float *coeff,                     /* i  : Coded spectral coefficients             */
    const short *sfm_start,                 /* i  : Subband start indices                   */
    const short *sfmsize,                   /* i  : Subband widths                          */
    const short *sfm_end,                   /* i  : Subband end indices                     */
    const short nb_sfm,                     /* i  : Number of subbands                      */
    const short *R,                         /* i  : Per-band bit allocation                 */
    float *CodeBook,                  /* o  : Noise-fill codebook                     */
    float *CodeBook_mod               /* o  : Densified noise-fill codebook           */
);

void apply_noisefill_HQ(
    const short *R,                         /* i  : bit allocation                          */
    const short length,                     /* i  : input frame length                      */
    const short flag_32K_env_ho,            /* i  : envelope stability hangover flag        */
    const long  core_brate,                 /* i  : core bit rate                           */
    const short last_sfm,                   /* i  : last coded subband                      */
    const float *CodeBook,                  /* i  : Noise-fill codebook                     */
    const float *CodeBook_mod,              /* i  : Densified noise-fill codebook           */
    const short cb_size,                    /* i  : Codebook length                         */
    const short *sfm_start,                 /* i  : Subband start coefficient               */
    const short *sfm_end,                   /* i  : Subband end coefficient                 */
    const short *sfmsize,                   /* i  : Subband band width                      */
    float *coeff                      /* i/o: coded/noisefilled spectrum              */
);

void harm_bwe_fine(
    const short *R,                         /* i  : bit allocation                          */
    const short last_sfm,                   /* i  : last coded subband                      */
    const short high_sfm,                   /* i  : higher transition band to BWE           */
    const short num_sfm,                    /* i  : total number of bands                   */
    const short *norm,                      /* i  : quantization indices for norms          */
    const short *sfm_start,                 /* i  : Subband start coefficient               */
    const short *sfm_end,                   /* i  : Subband end coefficient                 */
    short *prev_L_swb_norm,           /* i/o: last normalize length                   */
    float *coeff,                     /* i/o: coded/noisefilled normalized spectrum   */
    float *coeff_out,                 /* o  : coded/noisefilled spectrum              */
    float *coeff_fine                 /* o  : BWE fine structure                      */
);

void hvq_bwe_fine(
    const short last_sfm,                   /* i  : last coded subband                      */
    const short num_sfm,                    /* i  : total number of bands                   */
    const short *sfm_end,                   /* i  : Subband end coefficient                 */
    const short *peak_idx,                  /* i  : Peak index                              */
    const short Npeaks,                     /* i  : Number of peaks                         */
    short *peak_pos,                  /* i/o: Peak positions                          */
    short *prev_L_swb_norm,           /* i/o: last normalize length                   */
    float *coeff,                     /* i/o: coded/noisefilled normalized spectrum   */
    short *bwe_peaks,                 /* o  : Positions of peaks in BWE               */
    float *coeff_fine                 /* o  : HVQ BWE fine structure                  */
);

void hq_fold_bwe(
    const short last_sfm,                   /* i  : last coded subband              */
    const short *sfm_end,                   /* i  : Subband end coefficient         */
    const short num_sfm,                    /* i  : Number of subbands              */
    float *coeff                      /* i/o: coded/noisefilled normalized spectrum */
);

void apply_nf_gain(
    const short nf_idx,
    const short last_sfm,                   /* i  : last coded subband              */
    const short *R,                         /* i  : bit allocation                  */
    const short *sfm_start,                 /* i  : Subband start coefficient       */
    const short *sfm_end,                   /* i  : Subband end coefficient         */
    float *coeff                      /* i/o: coded/noisefilled normalized spectrum */

);

void hq_generic_fine(
    float *coeff,                     /* i  : coded/noisefilled normalized spectrum */
    const short last_sfm,                   /* i  : Last coded band                 */
    const short *sfm_start,                 /* i  : Subband start coefficient       */
    const short *sfm_end,                   /* i  : Subband end coefficient         */
    short *bwe_seed,                  /* i/o: random seed for generating BWE input */
    float *coeff_out1                 /* o  : HQ GENERIC input                */
);

void harm_bwe(
    const float *coeff_fine,                /* i  : fine structure for BWE                */
    const float *coeff,                     /* i  : coded/noisefilled normalized spectrum */
    const short num_sfm,                    /* i  : Number of subbands              */
    const short *sfm_start,                 /* i  : Subband start coefficient       */
    const short *sfm_end,                   /* i  : Subband end coefficient         */
    const short last_sfm,                   /* i  : last coded subband              */
    const short high_sfm,                   /* i  :                                 */
    const short *R,                         /* i  : bit allocation                  */
    const short prev_hq_mode,               /* i  : previous hq mode                */
    short *norm,                      /* i/o: quantization indices for norms  */
    float *noise_level,               /* i/o: noise levels for harmonic modes */
    float *prev_noise_level,          /* i/o: noise factor in previous frame  */
    short *bwe_seed,                  /* i/o: random seed for generating BWE input */
    float *coeff_out                  /* o  : coded/noisefilled spectrum      */
);

void hvq_bwe(
    const float *coeff,                     /* i  : coded/noisefilled spectrum              */
    const float *coeff_fine,                /* i  : BWE fine structure                      */
    const short *sfm_start,                 /* i  : Subband start coefficient               */
    const short *sfm_end,                   /* i  : Subband end coefficient                 */
    const short *sfm_len,                   /* i  : Subband length                          */
    const short last_sfm,                   /* i  : last coded subband                      */
    const short prev_hq_mode,               /* i  : previous hq mode                        */
    const short *bwe_peaks,                 /* i  : HVQ bwe peaks                           */
    const short bin_th,                     /* i  : HVQ transition bin                      */
    const short num_sfm,                    /* i  : Number of bands                         */
    const long  core_brate,                 /* i  : Core bit-rate                           */
    const short *R,                         /* i  : Bit allocation                          */
    short *norm,                      /* i/o: quantization indices for norms          */
    float *noise_level,               /* i/o: noise levels for harmonic modes         */
    float *prev_noise_level,          /* i/o: noise factor in previous frame          */
    short *bwe_seed,                  /* i/o: random seed for generating BWE input    */
    float *coeff_out                  /* o  : coded/noisefilled spectrum              */
);

void hvq_concat_bands
(
    const short pvq_bands,                  /* i  : Number of bands in concatenated PVQ target  */
    const short *sel_bnds,                  /* i  : Array of selected high bands                */
    const short n_sel_bnds,                 /* i  : Number of selected high bands               */
    short *hvq_band_start,            /* i  : Band start indices                          */
    short *hvq_band_width,            /* i  : Band widths                                 */
    short *hvq_band_end               /* i  : Band end indices                            */
);

void hq_generic_bwe(
    const short HQ_mode,                    /* i  : HQ mode                                     */
    float *coeff_out1,                /* i/o: BWE input & temporary buffer                */
    const float *hq_generic_fenv,           /* i  : SWB frequency envelopes                     */
    float *coeff_out,                 /* o  : SWB signal in MDCT domain                   */
    const short hq_generic_offset,          /* i  : frequency offset for representing hq generic*/
    short *prev_L_swb_norm,           /* i/o: last normalize length                       */
    const short hq_generic_exc_clas,        /* i  : hf excitation class                         */
    const short *sfm_end,                   /* i  : End of bands                                */
    const short num_sfm,
    const short num_env_bands,
    const short *R
);

void logqnorm_2(
    const float *env_fl,                    /* o  : index */
    const short L,                          /* i  : codebook length */
    const short n_env_band,                 /* i  : sub-vector size */
    const short nb_sfm,                     /* i  : sub-vector size */
    short *ynrm,
    short *normqlg2,
    const float *thren                      /* i  : quantization thresholds */
);

void map_hq_generic_fenv_norm(
    const short hqswb_clas,
    const float *hq_generic_fenv,
    short *ynrm,
    short *normqlg2,
    const short num_env_bands,
    const short nb_sfm,
    const short hq_generic_offset
);

short calc_nor_delta_hf(
    Encoder_State *st,
    const float *t_audio,
    short *ynrm,
    short *Rsubband,
    const short num_env_bands,
    const short nb_sfm,
    const short *sfmsize,
    const short *sfm_start,
    const short core_sfm
);

short get_nor_delta_hf(
    Decoder_State *st,
    short *ynrm,
    short *Rsubband,
    const short num_env_bands,
    const short nb_sfm,
    const short core_sfm
);

void hq_wb_nf_bwe(
    const float *coeff,                     /* i  : coded/noisefilled normal. spectrum  */
    const short is_transient,
    const short prev_bfi,                   /* i  : previous bad frame indicator        */
    const float *normq_v,
    const short num_sfm,                    /* i  : Number of subbands                  */
    const short *sfm_start,                 /* i  : Subband start coefficient           */
    const short *sfm_end,                   /* i  : Subband end coefficient             */
    const short *sfmsize,                   /* i  : Subband band width                  */
    const short last_sfm,                   /* i  : last coded subband                  */
    const short *R,                         /* i  : bit allocation                      */
    const short prev_is_transient,          /* i  : previous transient flag             */
    float *prev_normq,                /* i/o: previous norms                      */
    float *prev_env,                  /* i/o: previous noise envelopes            */
    short *bwe_seed,                  /* i/o: random seed for generating BWE input*/
    float *prev_coeff_out,            /* i/o: decoded spectrum in previous frame  */
    short *prev_R,                    /* i/o: previous frame bit allocation info. */
    float *coeff_out                  /* o  : coded/noisefilled spectrum          */
);

short encode_envelope_indices(              /* o  : Number of bits                          */
    Encoder_State *st,                        /* i/o: encoder state structure                 */
    const short num_sfm,                    /* i  : Number of subbands                      */
    const short numnrmibits,                /* i  : Bitrate of fall-back coding mode        */
    short *difidx,                    /* i/o: Diff indices/encoded diff indices       */
    short *LCmode,                    /* o  : Coding mode                             */
    const short flag_pack,                  /* i  : indicator of packing or estimating bits */
    const short flag_HQ2                    /* i  : indicator of HQ2 core                   */
    ,const short is_transient
);

void diff_envelope_coding(
    const short is_transient,               /* i  : transient indicator                 */
    const short num_env_bands,              /* i  : number of envelope bands to code    */
    const short start_norm,                 /* i  : start of envelope coding            */
    short *ynrm,                      /* i/o: quantization indices for norms      */
    short *normqlg2,                  /* i/o: quantized norms                     */
    short *difidx                     /* o  : differential code                   */
);

short decode_envelope_indices(              /* o  : Number of bits                      */
    Decoder_State *st,                        /* i/o: decoder state structure             */
    const short start_norm,                 /* i  : First SDE encoded norm              */
    const short num_sfm,                    /* i  : Number of norms                     */
    const short numnrmibits,                /* i  : Bitrate of fall-back coding mode    */
    short *ynrm,                      /* o  : Decoded norm indices                */
    const short flag_HQ2                    /* i  : indicator of HQ2 core               */
    ,const short is_transient
);

void dequantize_norms(                      /* o  : Number of bits                    */
    Decoder_State *st,                        /* i/o: decoder state structure             */
    const short start_norm,                 /* i  : First SDE encoded norm            */
    const short num_sfm,                    /* i  : Number of norms                   */
    const short is_transient,               /* i  : Transient flag                    */
    short *ynrm,                      /* o  : Decoded norm indices              */
    short *normqlg2                   /* o  : Log2 of decoded norms             */
);

void hq_configure(
    const short length,                     /* i  : Frame length                      */
    const short hqswb_clas,                 /* i  : HQ SWB class                      */
    const long  brate,                      /* i  : Codec bitrate                     */
    short *num_sfm,                   /* o  : Total number of subbands          */
    short *nb_sfm,                    /* o  : Total number of coded bands       */
    short *start_norm,                /* o  : First norm to be SDE encoded      */
    short *num_sde_norm,              /* o  : Number of norms for SDE encoding  */
    short *numnrmibits,               /* o  : Number of bits in fall-back norm encoding   */
    short *hq_generic_offset,         /* o  : Freq offset for HQ GENERIC        */
    short const **sfmsize,            /* o  : Subband bandwidths                */
    short const **sfm_start,          /* o  : Subband start coefficients        */
    short const **sfm_end             /* o  : Subband end coefficients          */
);

short hvq_enc(                              /* o  : Consumed bits                   */
    Encoder_State *st,                        /* i/o: encoder state structure         */
    const long  brate,                      /* i  : Total bit rate                  */
    const short hvq_bits,                   /* i  : HVQ bit budget                  */
    const short Npeaks,                     /* i  : Number of peaks                 */
    const short *ynrm,                      /* i  : Envelope coefficients           */
    short *R,                         /* i/o: Bit allocation/updated bit allocation */
    short *peaks,                     /* i  : Peak pos. / Encoded peak pos.   */
    float *nf_gains,                  /* i/o: Noise fill gains / Quant. nf gains */
    float *noise_level,               /* o  : Quantized noise level           */
    const float *pe_gains,                  /* i  : Peak gains                      */
    const float *coefs,                     /* i  : spectrum coefficients           */
    float *coefs_out                  /* o  : encoded spectrum coefficients   */
);

short hq_classifier_enc(                    /* o  : Consumed bits                   */
    Encoder_State *st,                      /* i/o: encoder state structure         */
    const short length,                     /* i  : Frame length                    */
    const float *coefs,                     /* i  : Spectral coefficients           */
    const short is_transient,               /* i  : Transient flag                  */
    short *Npeaks,                    /* o  : Number of identified peaks      */
    short *peaks,                     /* o  : Peak indices                    */
    float *pe_gains,                  /* o  : Peak gains                      */
    float *nf_gains,                  /* o  : Noise-fill gains                */
    short *hqswb_clas                 /* o  : HQ class                        */
);

short hq_classifier_dec(                    /* o  : Consumed bits                   */
    Decoder_State *st,                        /* i/o: decoder state structure         */
    const long  core_brate,                 /* i  : Core bit rate                   */
    const short length,                     /* i  : Frame length                    */
    short *is_transient,              /* o  : Transient flag                  */
    short *hqswb_clas                 /* o  : HQ class                        */
);


void hq_bit_allocation(
    const long core_brate,                  /* i  : Core bit-rate                      */
    const short length,                     /* i  : Frame length                       */
    const short hqswb_clas,                 /* i  : HQ class                           */
    short *num_bits,                  /* i/o: Remaining bit budget               */
    const short *normqlg2,                  /* i  : Quantized norms                    */
    const short nb_sfm,                     /* i  : Number sub bands to be encoded     */
    const short *sfmsize,                   /* i  : Sub band bandwidths                */
    float *noise_level,               /* o  : HVQ noise level                    */
    short *R,                         /* o  : Bit allocation per sub band        */
    short *Rsubband,                  /* o  : Fractional bit allocation (Q3)     */
    short *sum,                       /* o  : Sum of allocated shape bits        */
    short *core_sfm,                  /* o  : Last coded band in core            */
    const short num_env_bands
);

void enforce_zero_for_min_envelope(
    const short hqswb_clas,                 /* i  : HQ coding mode                     */
    const short *ynrm,                      /* i  : Envelope indices                   */
    float *coefsq,                    /* i/o: Quantized spectrum/zeroed spectrum */
    short  nb_sfm,                    /* i  : Number of coded sub bands          */
    const short *sfm_start,                 /* i  : Sub band start indices             */
    const short *sfm_end                    /* i  : Sub band end indices               */
);

short assign_gain_bits(                     /* o  : Number of assigned gain bits      */
    const short core,                       /* i  : HQ core                           */
    const short BANDS,                      /* i  : Number of bands                   */
    const short *band_width,                /* i  : Sub band bandwidth                */
    short *Rk,                        /* i/o: Bit allocation/Adjusted bit alloc. (Q3)*/
    short *gain_bits_array,           /* o  : Assigned gain bits                */
    short *Rcalc                      /* o  : Bit budget for shape quantizer (Q3)*/
);

void apply_envelope(
    const float *coeff,                 /* i/o: Coded/noisefilled normalized spectrum   */
    const short *norm,                  /* i  : Envelope                                */
    const float *norm_adj,              /* i  : Envelope adjustment                     */
    const short num_sfm,                /* i  : Total number of bands                   */
    const short last_sfm,               /* i  : Last coded band                         */
    const short HQ_mode,                /* i  : HQ mode                                 */
    const short length,                 /* i  : Frame length                            */
    const short *sfm_start,             /* i  : Sub band start indices                  */
    const short *sfm_end,               /* i  : Sub band end indices                    */
    float *normq_v,               /* o  : Envelope with adjustment                */
    float *coeff_out,             /* o  : coded/noisefilled spectrum              */
    float *coeff_out1             /* o  : noisefilled spectrum for HQ SWB BWE     */
);

PvqEntry mpvq_encode_vec(               /* o : Leading_sign_index, index, size, k_val       */
    const short* vec_in,                /* i : Signed pulse train                           */
    short dim_in,                 /* i : Dimension                                    */
    short k_val_local             /* i : Num unit pulses                              */
);

PvqEntry get_size_mpvq_calc_offset(      /* o : Size, dim, k_val            */
    short dim_in,                        /* i : Dimension                   */
    short k_val_in,                      /* i : Num unit pulses             */
    unsigned int* h_mem                  /* o : Offsets                     */
);

void mpvq_decode_vec(
    const PvqEntry* entry,              /* i :  Sign_ind, index, dim, k_val */
    unsigned int* h_mem,                /* i :  A/U offsets                 */
    short* vec_out                      /* o :  Pulse train                 */
);

unsigned int UMult_32_32(
    unsigned int UL_var1,
    unsigned int UL_var2
);

unsigned int UL_inverse(
    const unsigned int UL_val,
    short *exp
);

Word16 ratio(
    const Word32 numer,
    const Word32 denom,
    Word16 *expo
);

Word16 atan2_fx(                        /* o: Angle between 0 and EVS_PI/2 radian (Q14) */
    const Word32 y,                     /* i: Argument must be positive (Q15)       */
    const Word32 x                      /* i: Q15                                   */
);

void encode_energies(
    Encoder_State *st,
    const float *coefs,
    short Np,
    short *dim_part,
    float *E_part,
    short *bits_part,
    short *g_part,
    short qband,
    short *bits_left,
    float enr,
    short n,
    const short strict_bits
);

void decode_energies(
    Decoder_State *st,
    short Np,
    short *dim_part,
    short *bits_part,
    short *g_part,
    short qband,
    short *bits_left,
    short n,
    const short strict_bits
);

void pvq_encode_frame(
    Encoder_State *st,
    const float *coefs_norm,               /* i  : normalized coefficients to encode */
    float *coefs_quant,              /* o  : quantized coefficients */
    float *gopt,                     /* o  : optimal shape gains */
    short *npulses,                  /* o  : number of pulses per band */
    short *pulse_vector,             /* o  : non-normalized pulse shapes */
    const short *sfm_start,                /* i  : indices of first coefficients in the bands */
    const short *sfm_end,                  /* i  : indices of last coefficients in the bands */
    const short *sfmsize,                  /* i  : band sizes */
    const short nb_sfm,                    /* i  : total number of bands */
    const short *R,                        /* i  : bitallocation per band (Q3)*/
    const short pvq_bits,                  /* i  : number of bits avaiable */
    const short core                       /* i  : core */
);

void pvq_decode_frame(
    Decoder_State *st,
    float *coefs_quant,              /* o  : quantized coefficients */
    short *npulses,                  /* o  : number of pulses per band */
    short *pulse_vector,             /* o  : non-normalized pulse shapes */
    const short *sfm_start,                /* i  : indices of first coefficients in the bands */
    const short *sfm_end,                  /* i  : indices of last coefficients in the bands */
    const short *sfmsize,                  /* i  : band sizes */
    const short nb_sfm,                    /* i  : total number of bands */
    const short *R,                        /* i  : bitallocation per band (Q3) */
    const short pvq_bits,                  /* i  : number of bits avaiable */
    const short core                       /* i  : core */
);



void srt_vec_ind (
    const short *linear,                   /* linear input */
    short *srt,                      /* sorted output*/
    short *I,                        /* index for sorted output  */
    short length
);
void srt_vec_ind_f(
    const float *linear,                   /* linear input */
    float *srt,                      /* sorted output*/
    short *I,                        /* index for sorted output  */
    short length                     /* length of vector */
);


unsigned int floor_sqrt_exact(
    unsigned int input
);




void fine_gain_quant(
    Encoder_State *st,
    const short *ord,                       /* i  : Indices for energy order                     */
    const short num_sfm,                    /* i  : Number of bands                              */
    const short *gain_bits,                 /* i  : Gain adjustment bits per sub band            */
    float *fg_pred,                   /* i/o: Predicted gains / Corrected gains            */
    const float *gopt                       /* i  : Optimal gains                                */
);

void apply_gain(
    const short *ord,                       /* i  : Indices for energy order                     */
    const short *band_start,                /* i  : Sub band start indices                       */
    const short *band_end,                  /* i  : Sub band end indices                         */
    const short num_sfm,                    /* i  : Number of bands                              */
    const float *gains,                     /* i  : Band gain vector                             */
    float *xq                         /* i/o: Float synthesis / Gain adjusted synth        */
);

void fine_gain_pred(
    const short *sfm_start,                /* i  : Sub band start indices    */
    const short *sfm_end,                  /* i  : Sub band end indices      */
    const short *sfm_size,                 /* i  : Sub band bandwidths       */
    const short *i_sort,                   /* i  : Energy sorting indices    */
    const short *K,                        /* i  : Number of pulses per band */
    const short *maxpulse,                 /* i  : Maximum pulse per band    */
    const short *R,                        /* i  : Bits per sub band (Q3)    */
    const short num_sfm,                   /* i  : Number of sub bands       */
    float *xq,                       /* i/o: Quantized vector /quantized vector with finegain adj */
    short *y,                        /* i/o: Quantized vector (int)    */
    float *fg_pred,                  /* o  : Predicted fine gains      */
    const short core                       /* i  : Core                      */
);

void fine_gain_dec(
    Decoder_State *st,
    const short *ord,                       /* i  : Indices for energy order                     */
    const short num_sfm,                    /* i  : Number of bands                              */
    const short *gain_bits,                 /* i  : Gain adjustment bits per sub band            */
    float *fg_pred                    /* i/o: Predicted gains / Corrected gains            */
);

void get_max_pulses(
    const short *band_start,                /* i  : Sub band start indices    */
    const short *band_end,                  /* i  : Sub band end indices      */
    const short *k_sort,                    /* i  : Indices for sorting by energy */
    const short *npulses,                   /* i  : Pulses per sub band       */
    const short  BANDS,                     /* i  : Number of bands           */
    short *inp_vector,                /* i/o: Encoded shape vectors (int)*/
    short *maxpulse                   /* o  : Maximum pulse height per band */
);

Word32 ar_div(
    Word32 num,
    Word32 denum
);

void ar_encoder_start(
    PARCODEC arInst,
    PBITSTREAM bsInst,
    int max_bits
);

void ar_decoder_start(
    PARCODEC arInst,
    PBITSTREAM bsInst
);

void ar_encoder_done(
    PARCODEC arInst
);

void ar_decoder_done(
    PARCODEC arInst
);

float GetISCScale(
    float *quants,
    int size,
    Word32 bits_fx,
    float *magn,
    float *qscale,
    Word32 *surplus_fx,
    float *pulses,
    int* savedstates,
    int noTCQ,
    int *nzpout,
    short *bcount,
    float *abuffer,
    float *mbuffer,
    float *sbuffer
);

Word32 Mult_32_16(
    Word32 a,
    Word16 b
);

Word32 Mult_32_32(
    Word32 a,
    Word32 b
);

void decode_position_ari_fx(
    PARCODEC pardec,
    Word16 size,
    Word16 npulses,
    Word16 *nz,
    Word32 *position
);

void decode_magnitude_usq_fx(
    ARCODEC *pardec,
    Word16 size,
    Word16 npulses,
    Word16 nzpos,
    Word32 *positions,
    Word32 *out
);

void decode_mangitude_tcq_fx(
    ARCODEC *pardec,
    Word16 size,
    Word16 npulses,
    Word16 nzpos,
    Word32 *positions,
    Word32 *out,
    Word32 *surplus_fx
);

void decode_signs_fx(
    ARCODEC *pardec,
    Word16 size,
    Word32 *out
);

void srt_vec_ind_fx(
    const Word32 *linear,
    Word32 *srt,
    Word16 *I,
    Word16 length
);

Word16 GetScale_fx(
    Word16 blen,
    Word32 bits_fx/*Q16*/,
    Word32 *surplus_fx/*Q16*/
);

void bit_allocation_second_fx(
    Word32 *Rk,
    Word32 *Rk_sort,
    Word16  BANDS,
    const Word16 *band_width,
    Word16 *k_sort,
    Word16 *k_num,
    const Word16 *p2a_flags,
    const Word16  p2a_bands,
    const Word16 *last_bitalloc,
    const Word16  input_frame
);

Word32 encode_position_ari_fx(
    PARCODEC parenc,
    float *quants,
    Word16 size,
    Word32 *est_bits_frame_fx
);

Word32 encode_magnitude_tcq_fx(
    ARCODEC *parenc,
    float *magn_fx,
    Word16 size,
    Word16 npulses,
    Word16 nzpos,
    Word32 *savedstates,
    Word32 *est_frame_bits_fx
);

Word32 encode_signs_fx(
    ARCODEC *parenc,
    float *magn,
    Word16 size,
    Word16 npos,
    Word32 *est_frame_bits_fx
);

Word32 encode_magnitude_usq_fx(
    ARCODEC *parenc,
    float *magn_fx,
    Word16 size,
    Word16 npulses,
    Word16 nzpos,
    Word32 *est_frame_bits_fx
);

void tcq_core_LR_enc(
    Encoder_State *st,
    int   inp_vector[],
    const float coefs_norm[],
    float coefs_quant[],
    const short bit_budget,  /* number of bits */
    const short nb_sfm,
    const short *sfm_start,
    const short *sfm_end,
    const short *sfmsize,
    Word32 *Rk_fx,
    int   *npulses,
    short *k_sort,
    const short *p2a_flags,
    const short p2a_bands,
    const short *last_bitalloc,
    const short input_frame,
    const short adjustFlag,
    const short is_transient
);

void tcq_core_LR_dec(
    Decoder_State *st,
    int   *inp_vector,
    const short bit_budget,
    const short bands,
    const short *band_start,
    const short *band_width,
    Word32 *Rk_fx,
    int   *npulses,
    short *k_sort,
    const short *p2a_flags,
    const short p2a_bands,
    const short *last_bitalloc,
    const short input_frame,
    const short adjustFlag,
    const short *is_transient
);

void InitLSBTCQ(
    short *bcount
);

void TCQLSB(
    short bcount,
    float *abuffer,
    float *mbuffer,
    float *sbuffer,
    short *dpath
);

void RestoreTCQ(
    float * magn,
    int size,
    short *bcount,
    float *mbuffer
);

void SaveTCQdata(
    PARCODEC arInst,
    short *dpath,
    short bcount
);

void LoadTCQdata(
    PARCODEC arInst,
    short *dpath,
    short bcount
);

void RestoreTCQdec(
    int * magn,
    int size,
    short *bcount,
    float *mbuffer
);

void TCQLSBdec(
    short *dpath,
    float *mbuffer,
    short bcount
);

void bit_allocation_second_fx(
    Word32 *Rk,
    Word32 *Rk_sort,
    Word16  BANDS,
    const Word16 *band_width,
    Word16 *k_sort,
    Word16 *k_num,
    const Word16 *p2a_flags,
    const Word16  p2a_bands,
    const Word16 *last_bitalloc,
    const Word16  input_frame
);

void io_ini_enc(
    const int   argc,                       /* i  : command line arguments number             */
    char  *argv[],                    /* i  : command line arguments                    */
    FILE  **f_input,                  /* o  : input signal file                         */
    FILE  **f_stream,                 /* o  : output bitstream file                     */
    FILE  **f_rate,                   /* o  : bitrate switching profile (0 if N/A)      */
    FILE  **f_bwidth,                 /* o  : bandwidth switching profile (0 if N/A)    */
    FILE  **f_rf,                     /* o  : channel aware configuration file          */
    short  *quietMode,                /* o  : limit printouts                           */
    short  *noDelayCmp,               /* o  : turn off delay compensation               */
    Encoder_State *st                         /* o  : state structure                           */
);

void read_next_rfparam(
    short *rf_fec_offset,                  /* o: rf offset                         */
    short *rf_fec_indicator,               /* o: rf FEC indicator                  */
    FILE *f_rf                            /* i: file pointer to read parameters   */
);

void read_next_brate(
    long  *total_brate,                     /* i/o: total bitrate                             */
    const int last_total_brate,             /* i  : last total bitrate                        */
    FILE  *f_rate,                          /* i  : bitrate switching profile (0 if N/A)      */
    int   input_Fs,                         /* i  : input sampling frequency                  */
    short *Opt_AMR_WB,                      /* i  : flag indicating AMR-WB IO mode            */
    short *Opt_SC_VBR,                      /* i/o: SC-VBR flag                               */
    short *codec_mode                       /* i/o: Mode 1 or 2                               */
);

void read_next_bwidth(
    short *max_bwidth,                      /* i/o: maximum encoded bandwidth                 */
    FILE  *f_bwidth,                        /* i  : bandwidth switching profile (0 if N/A)    */
    long  *bwidth_profile_cnt,              /* i/o: counter of frames for bandwidth switching profile file */
    int    input_Fs                         /* i  : input sampling rate                       */
);


void init_encoder(
    Encoder_State *st                         /* i/o: state structure   */
);

void destroy_encoder(
    Encoder_State *st                         /* i/o: state structure   */
);

void evs_enc(
    Encoder_State *st,                        /* i/o: state structure             */
    const short *data,                      /* i  : input signal                */
    const short n_samples                   /* i  : number of input samples     */
);

void amr_wb_enc(
    Encoder_State *st,                        /* i/o: encoder state structure     */
    const short *data,                      /* i  : input signal                */
    const short n_samples                   /* i  : number of input samples     */
);

void pre_proc(
    Encoder_State *st,                        /* i/o: encoder state structure                  */
    const short input_frame,                /* i  : frame length                             */
    const float signal_in[],                /* i  : new samples                              */
    float old_inp_12k8[],             /* i/o: buffer of old input signal               */
    float old_inp_16k[],              /* i/o: buffer of old input signal @ 16kHz       */
    float **inp,                      /* o  : ptr. to inp. signal in the current frame */
    short *sp_aud_decision1,          /* o  : 1st stage speech/music classification    */
    short *sp_aud_decision2,          /* o  : 2nd stage speech/music classification    */
    float fr_bands[2*NB_BANDS],       /* o  : energy in frequency bands                */
    short *vad_flag,
    short *localVAD,
    float *Etot,                      /* o  : total energy, correlation shift          */
    float *ener,                      /* o  : residual energy from Levinson-Durbin     */
    short pitch[3],                   /* o  : open-loop pitch values for quantiz.      */
    float voicing[3],                 /* o  : OL maximum normalized correlation        */
    float A[NB_SUBFR16k*(M+1)],       /* o  : A(z) unquantized for the 4 subframes     */
    float Aw[NB_SUBFR16k*(M+1)],      /* o  : weighted A(z) unquantized for subframes  */
    float epsP[M+1],                  /* o  : LP prediction errors                     */
    float lsp_new[M],                 /* o  : LSPs at the end of the frame             */
    float lsp_mid[M],                 /* o  : LSPs in the middle of the frame          */
    short *coder_type,                /* o  : coder type                               */
    short *sharpFlag,                 /* o  : formant sharpening flag                  */
    short *vad_hover_flag,
    short *attack_flag,               /* o  : flag signalling attack encoded by AC mode (GSC)    */
    float *new_inp_resamp16k,         /* o  : new input signal @16kHz, non pre-emphasised, used by the WB TBE/BWE */
    short *Voicing_flag,              /* o  : voicing flag for HQ FEC                  */
    float realBuffer[CLDFB_NO_COL_MAX][CLDFB_NO_CHANNELS_MAX], /* i/o : real buffer   */
    float imagBuffer[CLDFB_NO_COL_MAX][CLDFB_NO_CHANNELS_MAX], /* i/o : imag buffer   */
    short *hq_core_type               /* o  : HQ core type                            */
);

short mdct_classifier(                      /* o: Class decision, 0 = Mode 1, 1 = Mode 2    */
    const float *Y,                         /* i: FFT spectrum from fft_rel                 */
    Encoder_State *st,                      /* i/o: Encoder state variable                  */
    short vadflag
    ,float *cldfbBuf_Ener
);

void MDCT_selector(
    Encoder_State *st,                        /* i/o: Encoder State                       */
    float sp_floor,                         /* i  : Noise floor estimate                */
    float Etot,                             /* i  : Total energy                        */
    float cor_map_sum,                      /* i  : harmonicity factor                  */
    const float voicing[],                  /* i  : voicing factors                     */
    const float enerBuffer[],               /* i  : CLDFB buffers                       */
    short vadflag
);

void MDCT_selector_reset(
    Encoder_State *st                        /* i/o: Encoder State                       */
);

void acelp_core_enc(
    Encoder_State *st,                        /* i/o: encoder state structure             */
    LPD_state *mem,                       /* i/o: acelp memories                      */
    const float inp[],                      /* i  : input signal of the current frame   */
    const short vad_flag,
    const float ener,                       /* i  : residual energy from Levinson-Durbin*/
    const short pitch[3],                   /* i  : open-loop pitch values for quantiz. */
    const float voicing[3],                 /* i  : Open-loop pitch gains               */
    float A[NB_SUBFR16k*(M+1)],       /* i  : A(z) unquantized for the 4 subframes*/
    float Aw[NB_SUBFR16k*(M+1)],      /* i  : weighted A(z) unquant. for subframes*/
    const float epsP[M+1],                  /* i  : LP prediction errors                */
    float lsp_new[M],                 /* i  : LSPs at the end of the frame        */
    float lsp_mid[M],                 /* i  : LSPs in the middle of the frame     */
    short coder_type,                 /* i  : coding type                         */
    const short sharpFlag,                  /* i  : formant sharpening flag             */
    short vad_hover_flag,
    const short attack_flag,                /* i  : flag signalling attack encoded by AC mode (GSC) */
    float bwe_exc_extended[],         /* i/o: bandwidth extended excitation       */
    float *voice_factors,             /* o  : voicing factors                     */
    float old_syn_12k8_16k[],         /* o  : ACELP core synthesis at 12.8kHz or 16kHz to be used by SWB BWE */
    float pitch_buf[NB_SUBFR16k],     /* o  : floating pitch for each subframe    */
    short *unbits                     /* o  : number of unused bits               */
);

void acelp_core_switch_dec_bfi(
    Decoder_State *st,                        /* i/o: decoder state structure                 */
    float synth_out[],                /* o  : synthesis                               */
    const short coder_type                  /* i  : coder type                              */
);

void acelp_core_switch_enc(
    Encoder_State *st,                        /* i/o: encoder state structure                 */
    LPD_state *mem,                       /* i/o: encoder memories                        */
    const float inp12k8[],                  /* i  : input signal @12.8 kHz                  */
    const float inp16k[],                   /* i  : input signal @16 kHz                    */
    const short T_op[2],                    /* i  : open-loop pitch values for quantiz.     */
    const float voicing[3],                 /* i  : Open-loop pitch gains                   */
    const float A[NB_SUBFR16k*(M+1)]        /* i  : A(z) unquantized for the 4 subframes    */
);

short modify_Fs_intcub3m_sup(               /* o  : length of output                        */
    const float sigIn[],                    /* i  : signal to decimate with memory of 2 samples (indexes -2 & -1) */
    const short lg,                         /* i  : length of input                         */
    const int   fin,                        /* i  : frequency of input                      */
    float sigOut[],                   /* o  : decimated signal                        */
    const int   fout,                       /* i  : frequency of output                     */
    short *delayout                   /* o  : delay of output                         */
);

void core_switching_OLA(
    const float *mem_over_hp,               /* i  : upsampling filter memory                */
    const short last_L_frame,               /* i  : last L_frame lengthture                 */
    const int   output_Fs,                  /* i  : output sampling rate                    */
    float *synth,                     /* i/o: synthesized signal from HQ core         */
    const float *synth_subfr_out,           /* i  : synthesized signal from ACELP core      */
    float *synth_subfr_bwe,           /* i  : synthesized BWE from ACELP core         */
    const short output_frame,               /* i  : output frame length                     */
    const short bwidth                      /* i  : output bandwidth                        */
);

void retro_interp4_5(
    const float *syn,
    float *pst_old_syn
);

void retro_interp5_4(
    float *pst_old_syn
);

void core_switching_hq_prepare_enc(
    Encoder_State *st,                        /* i/o: encoder state structure                 */
    short *num_bits,                  /* i/o: bit budget update                       */
    const short input_frame,                /* i  : input frame length                      */
    float *wtda_audio,
    const float *audio
);

void core_switching_hq_prepare_dec(
    Decoder_State *st,                        /* i/o: encoder state structure                 */
    short *num_bits,                  /* i/o: bit budget update                       */
    const short input_frame                 /* i  : input frame length                      */
);

void acelp_core_switch_dec(
    Decoder_State *st,                        /* i/o: decoder structure                       */
    float *synth_subfr_out,           /* o  : synthesized ACELP subframe              */
    float *tmp_synth_bwe,             /* o  : synthesized ACELP subframe BWE          */
    const short output_frame,               /* i  : input frame length                      */
    const short core_switching_flag,        /* i  : core switching flag                     */
    float *mem_synth                  /* o  : synthesis to overlap                    */
);

void space_lsfs(
    float *lsfs,
    const short order
);

void lsp2a(
    float *pc_in,
    float *freq,
    const short order
);

void lsp_weights(
    const float *lsps,
    float *weight,
    const short order
);

void a2lsp_stab(
    const float *a,                         /* i:   LP filter coefficients              */
    float *lsp,                       /* o:   Line spectral pairs                 */
    const float *old_lsp                    /* i:   LSP vector from past frame          */
);

void lsp2a_stab(
    const float *lsp,                       /* i  : LSF vector (in the cosine domain)   */
    float *a,                         /* o  : LP filter coefficients              */
    const short m                           /* i  : order of LP analysis                */
);

void isf2lsf(
    const float *isf,                       /* i  : ISF vector                          */
    float *lsf,                       /* o  : LSF vector                          */
    float *stable_lsp,                /* i/o: LSF vector                          */
    const short m,                          /* i  : order of LP analysis                */
    const float int_fs
);

void lsf2isf(
    const float *lsf,                       /* i  : LSF vector                          */
    float *isf,                       /* o  : ISF vector                          */
    float *stable_isp,                /* i/o: ISP vector                          */
    const short m,                          /* i  : order of LP analysis                */
    const float int_fs
);

short a2lsp(
    float *freq,                      /* o  : LSP vector                              */
    const float *a,                         /* i  : predictor coefficients                  */
    const short order                       /* i  : order of LP analysis                    */
);

void ResetSHBbuffer_Enc(
    Encoder_State *st                       /* i/o: encoder state structure                 */
);

void ResetSHBbuffer_Dec(
    Decoder_State *st                       /* i/o: decoder state structure                 */
);

void calc_st_filt(
    const float *apond2,                    /* i  : coefficients of numerator               */
    const float *apond1,                    /* i  : coefficients of denominator             */
    float *parcor0,                   /* o  : 1st parcor calcul. on composed filter   */
    float *sig_ltp_ptr,               /* i/o: input of 1/A(gamma1) : scaled by 1/g0   */
    float *mem_zero,                  /* i/o: All zero memory                         */
    const short L_subfr,                    /* i  : the length of subframe                  */
    const short extl                        /* i  : extension layer info                    */
);

void scale_st(
    const float *sig_in,                    /* i  : postfilter input signal                 */
    float *sig_out,                   /* i/o: postfilter output signal                */
    float *gain_prec,                 /* i/o: last value of gain for subframe         */
    const short L_subfr,                    /* i  : the length of subframe                  */
    const short extl                        /* i  : extension layer info                    */
);

void filt_mu(
    const float *sig_in,                    /* i  : signal (beginning at sample -1)         */
    float *sig_out,                   /* o  : output signal                           */
    const float parcor0,                    /* i  : parcor0 (mu = parcor0 * gamma3)         */
    const short L_subfr,                    /* i  : the length of subframe                  */
    const short extl                        /* i  : extension layer info                    */
);

void PostShortTerm(
    float *sig_in,                    /* i  : input signal (ptr. to current subframe  */
    float *lpccoeff,                  /* i  : LPC coefficients for current subframe   */
    float *sig_out,                   /* o  : postfiltered output                     */
    float *mem_stp,                   /* i/o: postfilter memory                       */
    float *ptr_mem_stp,               /* i/o: pointer to postfilter memory            */
    float *ptr_gain_prec,             /* i/o: for gain adjustment                     */
    float *mem_zero,                  /* i/o: null memory to compute h_st             */
    const float formant_fac                 /* i  : Strength of post-filter [0,1]           */
);

float swb_formant_fac(                      /* o  : Formant filter strength [0,1]           */
    const float lpc_shb2,                   /* i  : 2nd HB LPC coefficient                  */
    float *tilt_mem                   /* i/o: Tilt smoothing memory                   */
);

void GenShapedSHBExcitation(
    float *excSHB,                    /* o  : synthesized shaped shb exctiation       */
    const float *lpc_shb,                   /* i  : lpc coefficients                        */
    float *exc16kWhtnd,               /* o  : whitened synthesized shb excitation     */
    float *mem_csfilt,                /* i/o: memory                                  */
    float *mem_genSHBexc_filt_down_shb,/* i/o: memory                                 */
    float *state_lpc_syn,             /* i/o: memory                                  */
    const short coder_type,                 /* i  : coding type                             */
    const float *bwe_exc_extended,          /* i  : bandwidth extended excitation           */
    short bwe_seed[],                 /* i/o: random number generator seed            */
    float voice_factors[],            /* i  : voicing factor                          */
    const short extl,                       /* i  : extension layer                         */
    float *tbe_demph,                 /* i/o: de-emphasis memory                      */
    float *tbe_premph,                /* i/o: pre-emphasis memory                     */
    float *lpc_shb_sf,                /* i  : LP coefficients                         */
    float *shb_ener_sf,               /* i  : SHB subframe energies                   */
    float *shb_res_gshape,            /* i  : SHB LP residual gain shape              */
    float *shb_res,                   /* i  : SHB residual used in encoder only       */
    short *vf_ind,                    /* i/o: Mixing factor index                     */
    const float formant_fac,                /* i   : Formant sharpening factor [0..1]       */
    float fb_state_lpc_syn[],         /* i/o: memory                                  */
    float *fb_tbe_demph,              /* i/o: fb de-emphasis memory                   */
    const long bitrate,                     /* i  : overall bitrate                         */
    const short prev_bfi                    /* i  : previous frame was lost flag            */
);

void Estimate_mix_factors(
    const float *shb_res,                   /* i  : SHB LP residual */
    const float *exc16kWhtnd,               /* i  : SHB transformed low band excitation */
    const float *White_exc16k,              /* i  : Modulated envelope shaped white noise  */
    const float pow1,                       /* i  : SHB exc. power for normalization */
    const float pow2,                       /* i  : White noise excitation for normalization */
    float *vf_modified,               /* o  : Estimated voice factors */
    short *vf_ind                     /* o  : voice factors VQ index */
);

void GenSHBSynth(
    const float *shb_target_speech,         /* i  : input synthesized speech                */
    float *shb_syn_speech_32k,        /* o  : output highband component               */
    float Hilbert_Mem[],              /* i/o: memory                                  */
    float state_lsyn_filt_shb_local[],/* i/o: memory                                  */
    const short L_frame,                    /* i  : ACELP Frame length                      */
    short *syn_dm_phase
);

void ScaleShapedSHB(
    const short length,                     /* i  : SHB overlap length                      */
    float *synSHB,                    /* i/o: synthesized shb signal                  */
    float *overlap,                   /* i/o: buffer for overlap-add                  */
    const float *subgain,                   /* i  : subframe gain                           */
    const float frame_gain,                 /* i  : frame gain                              */
    const float *win,                       /* i  : window                                  */
    const float *subwin                     /* i  : subframes window                        */
);

void Interpolate_allpass_steep (
    const float *in,                        /* i  : input array of size N                   */
    float *state,                     /* i/o: memory                                  */
    const short N,                          /* i  : number of input samples                 */
    float *out                        /* o  : output array of size 2*N                */
);

void Decimate_allpass_steep (
    const float *in,                        /* i  : input array of size N                   */
    float *state,                     /* i/o: memory                                  */
    const short N,                          /* i  : number of input samples                 */
    float *out                        /* o  : output array of size N/2                */
);

void  interpolate_3_over_2_allpass(
    const float *input,                   /* i  : input signal                            */
    const short len,                      /* i  : number of input samples                 */
    float *out,                     /* o  : output signal                           */
    float *mem,                     /* i/o: memory                                  */
    const float *filt_coeff               /* i  : filter coefficients                     */
);

void decimate_2_over_3_allpass(
    const float *input,                   /* i  : input signal                            */
    const short len,                      /* i  : number of input samples                 */
    float *out,                     /* o  : output signal                           */
    float *mem,                     /* i/o: memory                                  */
    const float *filt_coeff,              /* i  : filter coefficients                     */
    const float *lp_num,
    const float *lp_den,
    float *lp_mem
);

void interpolate_3_over_1_allpass(
    const float *input,                   /* i  : input signal                            */
    const short len,                      /* i  : number of input samples                 */
    float *out,                     /* o  : output signal                           */
    float *mem,                     /* i/o: memory                                  */
    const float *filt_coeff               /* i  : filter coefficients                     */
);

void InitSWBencBuffer(
    Encoder_State *st                         /* i/o: encoder state structure                 */
);

void swb_tbe_enc(
    Encoder_State *st,                        /* i/o: encoder state structure                 */
    const short coder_type,                 /* i  : coding type                             */
    float *new_speech,                /* i  : original input signal                   */
    const float *bwe_exc_extended,          /* i  : bandwidth extended exciatation          */
    const float voice_factors[],            /* i  : voicing factors                         */
    float *White_exc16k,              /* o  : shaped white excitation for the FB TBE  */
    const float voicing[],                  /* i  : OL maximum normalized correlation       */
    const float pitch_buf[]                 /* i  : pitch for each subframe                 */
);

void InitSWBdecBuffer(
    Decoder_State *swb_dnc                    /* i/o: SHB decoder structure                   */
);

void swb_tbe_dec(
    Decoder_State *st,                        /* i/o: decoder state structure                 */
    const short coder_type,                 /* i  : coding type */
    const float *bwe_exc_extended,          /* i  : bandwidth extended exciatation          */
    const float voice_factors[],            /* i  : voicing factors                         */
    const float old_syn_12k8_16k[],
    float *White_exc16k,              /* o  : shaped white excitation for the FB TBE  */
    float *synth,                     /* i/o: ACELP core synthesis/final synthesis    */
    float *pitch_buf
);

void flip_and_downmix_generic(
    float input[],                    /* i  : input spectrum                          */
    float output[],                   /* o  : output  spectrum                        */
    const short length,                     /* i  : length of spectra                       */
    float mem1_ext[HILBERT_ORDER1],   /* i/o: Hilbert filter memory                   */
    float mem2_ext[2*HILBERT_ORDER2], /* i/o: memory                                  */
    float mem3_ext[2*HILBERT_ORDER2], /* i/o: memory                                  */
    short *phase_state                /* i/o: Phase state in case frequency isn't multiple of 50 Hz */
);

void non_linearity(
    const float input[],                    /* i  : input signal                            */
    float output[],                   /* i  : output signal                           */
    float old_bwe_exc_extended[],     /* i/o: memory bugffer                          */
    const short length,                     /* i  : input length                            */
    float *prev_scale                 /* i/o: memory                                  */
    ,short  coder_type,                 /* i  : Coder Type                              */
    float   *voice_factors,             /* i  : Voice Factors                           */
    const short  L_frame			        /* i  : ACELP frame length                      */
);

void interp_code_5over2(
    const float inp_code[],                /* i  :  input vector                           */
    float interp_code[],             /* o  :  output vector                          */
    const short inp_length                 /* i  :  length of the input vector             */
);

void interp_code_4over2(
    const float inp_code[],                 /* i  :  input vector                           */
    float interp_code[],              /* o  :  output vector                          */
    const short inp_length                  /* i  :  length of the input vector             */
);

void flip_spectrum_and_decimby4(
    const float input[],                    /* i  : input spectrum                          */
    float output[],                   /* o  : output  spectrum                        */
    const short length,                     /* i  : vector length                           */
    float mem1[],                     /* i/o: memory                                  */
    float mem2[],                     /* i/o: memory                                  */
    const short ramp_flag                   /* i  : flag to trigger slow ramp-up of output  */
);

void GenShapedWBExcitation(
    float *excSHB,                    /* o   : synthesized shaped shb exctiation      */
    const float *lpc_shb,                   /* i   : lpc coefficients                       */
    float *exc4kWhtnd,                /* o   : whitened synthesized shb excitation    */
    float *mem_csfilt,                /* i/o : memory                                 */
    float *mem_genSHBexc_filt_down1,  /* i/o : memory                                 */
    float *mem_genSHBexc_filt_down2,  /* i/o : memory                                 */
    float *mem_genSHBexc_filt_down3,  /* i/o : memory                                 */
    float *state_lpc_syn,             /* i/o : memory                                 */
    const short coder_type,                 /* i   : coding type                            */
    const float *bwe_exc_extended,          /* i   : bandwidth extended exciatation         */
    short bwe_seed[],                 /* i/o : random number generator seed           */
    const float voice_factors[],            /* i   : voicing factor                         */
    const short uv_flag,                    /* i   : unvoiced flag                          */
    const short igf_flag
);

void GenWBSynth(
    const float *input_synspeech,           /* i  : input synthesized speech                */
    float *shb_syn_speech_16k,        /* o  : output highband compnent                */
    float *state_lsyn_filt_shb1,      /* i/o: memory                                  */
    float *state_lsyn_filt_shb2       /* i/o: memory                                  */
);

void wb_tbe_enc(
    Encoder_State *st,                         /* i/o: encoder state structure                */
    const short coder_type,                  /* i  : coding type                            */
    const float *hb_speech,                  /* i  : HB target signal (6-8kHz) at 16kHz     */
    const float *bwe_exc_extended,           /* i  : bandwidth extended exciatation         */
    const float voice_factors[],             /* i  : voicing factors                        */
    const float pitch_buf[],                 /* i  : pitch for each subframe                */
    const float voicing[]                    /* o  : OL maximum normalized correlation      */
);

void wb_tbe_dec(
    Decoder_State *st,                        /* i/o: decoder state structure                 */
    const short coder_type,                 /* i  : coding type                             */
    const float *bwe_exc_extended,          /* i  : bandwidth extended exciatation          */
    const float voice_factors[],            /* i  : voicing factors                         */
    float *synth                      /* i/o: ACELP core synthesis/final synthesis    */
);

void tbe_write_bitstream(
    Encoder_State *st                         /* i/o: encoder state structure                 */
);

void tbe_read_bitstream(
    Decoder_State *st                         /* i/o: decoder state structure                 */
);

void GenTransition(
    const float *input,                     /* i  : gain shape overlap buffer              */
    const float *old_hb_synth,              /* i  : synthesized HB from previous frame     */
    short length,                     /* i  : targeted length of transition signal   */
    float *output,                    /* o  : synthesized transitions signal         */
    float Hilbert_Mem[],              /* i/o: memory                                 */
    float state_lsyn_filt_shb_local[],/* i/o: memory                                 */
    short *syn_dm_phase,
    int   output_Fs,
    float *up_mem,
    int   rf_flag
    , int bitrate
);

void GenTransition_WB(
    const float *input,                         /* i  : gain shape overlap buffer            */
    const float *old_hb_synth,                  /* i  : synthesized HB from previous frame   */
    short length,                         /* i  : targeted length of transition signal */
    float *output,                        /* o  : synthesized transitions signal       */
    float state_lsyn_filt_shb1[],
    float state_lsyn_filt_shb2[],
    int   output_Fs,
    float *up_mem
);

void TBEreset_enc(
    Encoder_State *st,                        /* i/o: encoder state structure                 */
    short  bandwidth                  /* i  : bandwidth mode                          */
);

void TBEreset_dec(
    Decoder_State *st,                        /* i/o: decoder state structure                 */
    short  bandwidth                  /* i  : bandwidth mode                          */
);

short get_tbe_bits(                         /* o  : TBE bit consumption per frame           */
    short bitrate,                   /* i  : overall bitrate                         */
    short bandwidth,                 /* i  : bandwidht mode                          */
    short rf_mode                    /* i  : channel aware mode                      */
);

void fb_tbe_enc(
    Encoder_State *st,                        /* i/o: encoder state structure                 */
    const float new_input[],                /* i  : input speech at 48 kHz sample rate      */
    const float fb_exc[]                    /* i  : FB excitation from the SWB part         */
);

void fb_tbe_dec(
    Decoder_State *st,                        /* i/o: decoder state structure                 */
    const float fb_exc[],                   /* i  : FB excitation from the SWB part         */
    float *hb_synth                   /* i/o: high-band synthesis                     */
);

void calc_tilt_bwe(
    const float *sp,                        /* i  : input signal                            */
    float *tilt,                      /* o  : signal tilt                             */
    const short N                           /* i  : signal length                           */
);

void wtda_BWE(
    const float *new_audio,                 /* i  : input audio                             */
    float *old_wtda,                  /* i/o: windowed audio from previous frame      */
    const short L                           /* i  : length                                  */
);

void swb_pre_proc(
    Encoder_State *st,                         /* i/o: encoder state structure                 */
    const float *input,                     /* i  : original input signal                   */
    float *new_swb_speech,            /* o  : original input signal at 32kHz          */
    float *shb_speech,                /* o  : SHB target signal (6-14kHz) at 16kHz    */
    float realBuffer[CLDFB_NO_COL_MAX][CLDFB_NO_CHANNELS_MAX], /* i : real buffer     */
    float imagBuffer[CLDFB_NO_COL_MAX][CLDFB_NO_CHANNELS_MAX]  /* i : imag buffer     */
);

void wb_pre_proc(
    Encoder_State *st,                         /* i/o: encoder state structure                 */
    const float *new_inp_resamp16k,          /* i  : original input signal                   */
    float *hb_speech                   /* o  : HB target signal (6-8kHz) at 16kHz      */
);

void wb_bwe_enc(
    Encoder_State *st,                         /* i/o: encoder state structure                 */
    const float *new_wb_speech,              /* i  : original input signal at 16kHz          */
    short coder_type                   /* i  : coding type                             */
);

void wb_bwe_dec(
    float *synth,                      /* i/o: ACELP core synthesis/final synthesis    */
    float *hb_synth,                   /* o  : SHB synthesis/final synthesis           */
    const short output_frame,                /* i  : frame length                            */
    Decoder_State *st,                         /* i/o: decoder state structure                 */
    short coder_type,                  /* i  : coding type                             */
    float *voice_factors,              /* i  : voicing factors                         */
    const float pitch_buf[]                  /* i  : pitch buffer                            */
);

void swb_bwe_enc(
    Encoder_State *st,                        /* i/o: encoder state structure                 */
    const float *old_input_12k8,            /* i  : input signal @12.8kHz for SWB BWE       */
    const float *old_input_16k,             /* i  : input signal @16kHz for SWB BWE         */
    const float *old_syn_12k8_16k,          /* i  : ACELP core synthesis at 12.8kHz or 16kHz*/
    const float *new_swb_speech,            /* i  : original input signal at 32kHz          */
    const float *shb_speech,                /* i  : SHB target signal (6-14kHz) at 16kHz    */
    const short coder_type                  /* i  : coding type                             */
);

void swb_bwe_enc_hr(
    Encoder_State *st,                        /* i/o: encoder state structure                 */
    const float *new_input,                 /* i  : input signal                            */
    const short input_frame,                /* i  : frame length                            */
    const short coder_type,                 /* i  : coding type                             */
    const short unbits                      /* i  : number of core unused bits              */
);

void swb_bwe_dec(
    Decoder_State *st,                         /* i/o: decoder state structure                 */
    const float *synth,                      /* i  : ACELP core synthesis/final synthesis    */
    float *hb_synth,                   /* o  : SHB synthesis/final synthesis           */
    const short output_frame                 /* i  : frame length                            */
    ,short coder_type                        /* i  : coding type                             */
);

void swb_bwe_dec_hr(
    Decoder_State *st,                         /* i/o: decoder state structure                 */
    const float *syn_12k8_16k,               /* i  : ACELP core synthesis @16kHz             */
    float *hb_synth,                   /* o  : SHB synthesis                           */
    const short output_frame,                /* i  : frame length                            */
    const short unbits,                      /* i  : number of core unused bits              */
    const float pitch_buf[]                  /* i  : pitch buffer                            */
);

void swb_hr_noise_fill(
    const short is_transient,                 /* i  : transient flag                          */
    const short spect_start,                  /* i  : spectrum start point                    */
    const short spect_end,                    /* i  : spectrum end point                      */
    const float tilt_wb,                      /* i  : tilt of wideband signal                 */
    const float pitch,                        /* i  : pitch value                             */
    const short nq[],                         /* i  : AVQ nq index                            */
    short Nsv,                          /* i  : number of subband                       */
    short *bwe_highrate_seed,           /* i/o: seed of random noise                    */
    float *t_audio                      /* i/o: mdct spectrum                           */
);

float td_postprocess(                       /* o  : gain                                    */
    float hb_synth[],                 /* i/o: high-band synthesis                     */
    const short input_frame,                /* i  : frame length                            */
    const short last_extl                   /* i  : last extension layer                    */
);

void calc_normal_length(
    const short core,                       /* i  : core                                    */
    const float *sp,                        /* i  : input signal                            */
    const short mode,                       /* i  : input mode                              */
    const short extl,                       /* i  : extension layer                         */
    short *L_swb_norm,                /* o  : normalize length                        */
    short *prev_L_swb_norm            /* i/o: last normalize length                   */
);

void calc_norm_envelop(
    const float SWB_signal[],                /* i  : SWB spectrum                            */
    float *envelope,                   /* o  : normalized envelope                     */
    const short L_swb_norm,                  /* i  : length of envelope                      */
    const short SWB_flength,                 /* i  : Length of input/output                  */
    const short st_offset                    /* i  : offset                                  */
);

void time_envelop_shaping(
    float werr[],                      /* i/o: SHB synthesis                           */
    float SWB_tenv[],                  /* i/o: frequency envelope                      */
    const short L                            /* i  : frame length                            */
);

void time_reduce_pre_echo(
    const float *synth,                      /* i  : ACELP core synthesis                    */
    float *error,                      /* o  : SHB BWE synthesis                       */
    float prev_td_energy,              /* o  : last td energy                          */
    const short L                            /* i  : subframe length                         */
);

short WB_BWE_gain_pred(
    float *WB_fenv,                    /* o  : WB frequency envelopes                  */
    const float *core_dec_freq,              /* i  : Frequency domain core decoded signal    */
    const short coder_type,                  /* i  : coding type                             */
    short prev_code_type,              /* i  : coding type of last frame               */
    float prev_WB_fenv,                /* i  : envelope for last frame                 */
    float *voice_factors,              /* i  : voicing factors                         */
    const float pitch_buf[],                 /* i  : pitch buffer                            */
    long  last_core_brate,             /* i  : previous frame core bitrate             */
    float last_wb_bwe_ener             /* i  : previous frame wb bwe signal energy     */
    ,short last_extl                    /* i  : extl. layer for last frame              */
    ,float tilt
);

void WB_BWE_decoding(
    const float *core_dec_freq,              /* i  : Frequency domain core decoded signal    */
    float *WB_fenv,                    /* i  : WB frequency envelopes                  */
    float *WB_signal,                  /* o  : WB signal in MDCT domain                */
    const short WB_flength,                  /* i  : Length of input/output                  */
    const short mode,                        /* i  : classification for WB signal            */
    const short last_extl,                   /* i  : extl. layer for last frame              */
    float *prev_Energy,                /* i/o: energy for last frame                   */
    float *prev_WB_fenv,               /* i/o: envelope for last frame                 */
    short *prev_L_wb_norm,             /* i/o: length for last frame wb norm           */
    const short extl,                        /* i  : extension layer                         */
    const short coder_type,                  /* i  : coding type                             */
    const long  total_brate,                 /* i  : core layer bitrate                      */
    short *Seed,                       /* i/o: random generator seed                   */
    short *prev_flag,                  /* i/o: attenu flag of last frame               */
    short prev_coder_type              /* i  : coding type of last frame               */
);

void SWB_BWE_decoding(
    const float *core_dec_freq,              /* i  : Frequency domain core decoded signal    */
    float *SWB_fenv,                   /* i/o: SWB frequency envelopes                 */
    float *SWB_signal,                 /* o  : SWB signal in MDCT domain               */
    const short SWB_flength,                 /* i  : Length of input/output                  */
    const short mode,                        /* i  : classification for SWB signal           */
    short *frica_flag,                 /* o  : fricative signal flag                   */
    float *prev_Energy,                /* i/o: energy for last frame                   */
    float *prev_SWB_fenv,              /* i/o: envelope for last frame                 */
    short *prev_L_swb_norm,            /* i/o: length for last frame wb norm           */
    const float tilt_nb,                     /* i  : tilt of synthesis wb signal             */
    short *Seed,                       /* i/o: random generator seed                   */
    const short st_offset,                   /* i  : offset value due to different core      */
    float *prev_weight,                /* i/o: excitation weight value of last frame   */
    const short extl                         /* i  : extension layer                         */
    ,const short last_extl                    /* i  : extension layer of last frame           */
);

void CNG_reset_enc(
    Encoder_State *st,                        /* i/o: encoder state structure                         */
    LPD_state *mem,                       /* i/o: encoder memories                                */
    float *pitch_buf,                 /* o  : floating pitch for each subframe                */
    float *voice_factors              /* o  : voicing factors                                 */
    ,short VBR_cng_reset_flag
);

void a2isp(
    const float *a,                         /* i  : LP filter coefficients                          */
    float *isp,                       /* o  : Immittance spectral pairs                       */
    const float *old_isp                    /* i  : ISP vector from past frame                      */
);

void a2isf(
    float *a,
    float *isf,
    const float *old_isf,
    short lpcOrder);

void a2rc (
    const float *a,                         /* i  : LPC coefficients                                */
    float *refl,                      /* o  : Reflection co-efficients                        */
    const short lpcorder                    /* i  : LPC order                                       */
);

short lp_filt_exc_enc(
    const short codec_mode,                 /* i  : codec mode                                      */
    const long  core_brate,                 /* i  : core bitrate                                    */
    const short Opt_AMR_WB,                 /* i  : flag indicating AMR-WB IO mode                  */
    const short coder_type,                 /* i  : coding type                                     */
    const short i_subfr,                    /* i  : subframe index                                  */
    float *exc,                       /* i/o: pointer to excitation signal frame              */
    const float *h1,                        /* i  : weighted filter input response                  */
    const float *xn,                        /* i  : target vector                                   */
    float *y1,                        /* o  : zero-memory filtered adaptive excitation        */
    float *xn2,                       /* o  : target vector for innovation search             */
    const short L_subfr,                    /* i  : length of vectors for gain quantization         */
    const short L_frame,                    /* i  : frame size                                      */
    float *g_corr,                    /* o  : ACELP correlation values                        */
    const short clip_gain,                  /* i  : adaptive gain clipping flag                     */
    float *gain_pit,                  /* o  : adaptive excitation gain                        */
    short *lp_flag                    /* i/o  : mode selection                                */
);

void updt_tar(
    const float *x,                         /* i  : old target (for pitch search)                   */
    float *x2,                        /* o  : new target (for codebook search)                */
    const float *y,                         /* i  : filtered adaptive codebook vector               */
    const float gain,                       /* i  : adaptive codebook gain                          */
    const short L                           /* i  : subframe size                                   */
);

void analy_sp(
    float *speech,                    /* i  : speech buffer                                    */
    float *Bin_E,                     /* o  : per bin log energy spectrum                      */
    float *Bin_E_old,                 /* o  : per bin log energy spectrum for mid-frame        */
    float *fr_bands,                  /* o  : per band energy spectrum (2 analyses)            */
    float lf_E[],                     /* o  : per bin E for first VOIC_BINS bins (without DC)  */
    float *Etot,                      /* o  : total input energy                               */
    const short min_band,                   /* i  : minimum critical band                            */
    const short max_band,                   /* i  : maximum critical band                            */
    float *band_ener,                 /* o: energy in critical frequency bands without minimum noise floor E_MIN */
    float *PS                         /* o  : Per bin energy spectrum                          */
    ,float *fft_buff                   /* o  : FFT coefficients                                 */
);

void CNG_enc(
    Encoder_State *st,                        /* i/o: State structure                                 */
    const short L_frame,                    /* i  : length of the frame                             */
    float Aq[],                       /* o  : LP coefficients                                 */
    const float *speech,                    /* i  : pointer to current frame input speech buffer    */
    float enr,                        /* i  : frame energy output from Levinson recursion     */
    float *lsp_new,                   /* i/o: current frame LSPs                              */
    float *lsf_new,                   /* i/o: current frame LSFs                              */
    short *allow_cn_step,             /* o  : allow CN step                                   */
    short burst_ho_cnt,               /* i  : hangover frames at end of speech burst          */
    float *q_env,
    short *sid_bw,
    float *exc_mem2
);

void swb_CNG_enc(
    Encoder_State *st,                        /* i/o: State structure                                 */
    const float *shb_speech,                /* i  : SHB target signal (6-14kHz) at 16kHz            */
    const float *syn_12k8_16k               /* i  : ACELP core synthesis at 12.8kHz or 16kHz        */
);

void lsf_enc(
    Encoder_State *st,                        /* i/o: state structure                             */
    const short L_frame,                    /* i  : length of the frame                         */
    const short coder_type,                 /* i  : coding type                                 */
    float *lsf_new,                   /* o  : quantized LSF vector                        */
    float *lsp_new,                   /* i/o: LSP vector to quantize/quantized            */
    float *lsp_mid,                   /* i  : mid-frame LSP vector                        */
    float *Aq,                        /* o  : quantized A(z) for 4 subframes              */
    float *stab_fac,                  /* o  : ISF stability factor                        */
    const short Nb_ACELP_frames
);

void isf_enc_amr_wb(
    Encoder_State *st,                        /* i/o: state structure                             */
    float *isf_new,                   /* o  : quantized ISF vector                        */
    float *isp_new,                   /* i/o: ISP vector to quantize/quantized            */
    float *Aq,                        /* o  : quantized A(z) for 4 subframes              */
    float *stab_fac                   /* o  : ISF stability factor                        */
);

void find_targets(
    const float *speech,                    /* i  : pointer to the speech frame                      */
    const float *mem_syn,                   /* i  : memory of the synthesis filter                   */
    const short i_subfr,                    /* i  : subframe index                                   */
    float *mem_w0,                    /* i/o: weighting filter denominator memory              */
    const float *p_Aq,                      /* i  : interpolated quantized A(z) filter               */
    const float *res,                       /* i  : residual signal                                  */
    const short L_subfr,                    /* i  : length of vectors for gain quantization          */
    const float *Ap,                        /* i  : unquantized A(z) filter with bandwidth expansion */
    const float tilt_fac,                   /* i  : tilt factor                                      */
    float *xn,                        /* o  : Close-loop Pitch search target vector            */
    float *cn,                        /* o  : target vector in residual domain                 */
    float *h1                         /* o  : impulse response of weighted synthesis filter    */
);

void inov_encode(
    Encoder_State *st,                        /* i/o: encoder state structure                         */
    const  long core_brate,                 /* i  : core bitrate                                    */
    const short Opt_AMR_WB,                 /* i  : flag indicating AMR-WB IO mode                  */
    const short L_frame,                    /* i  : length of the frame                             */
    const short last_L_frame,               /* i  : length of the last frame                        */
    const short coder_type,                 /* i  : coding type                                     */
    const short bwidth,                     /* i  : input signal bandwidth                          */
    const short sharpFlag,                  /* i  : formant sharpening flag                         */
    const short i_subfr,                    /* i  : subframe index                                  */
    short const tc_subfr,                   /* i  : TC subframe index                               */
    const float *p_Aq,                      /* i  : LP filter coefficients                          */
    const float gain_pit,                   /* i  : adaptive excitation gain                        */
    float *cn,                        /* i/o: target vector in residual domain                */
    const float *exc,                       /* i  : pointer to excitation signal frame              */
    float *h1,                        /* i/o: weighted filter input response                  */
    const float tilt_code,                  /* i  : tilt of of the excitation of previous subframe  */
    const float pt_pitch,                   /* i  : pointer to current subframe fractional pitch    */
    const float *xn2,                       /* i  : target vector for innovation search             */
    float *code,                      /* o  : algebraic excitation                            */
    float *y2,                        /* o  : zero-memory filtered algebraic excitation       */
    short *unbits                     /* o  : number of unused bits for  EVS_PI               */
);

void acelp_1t64(
    Encoder_State *st,                        /* i/o: encoder state structure                       */
    const float dn[],                       /* i  : corr. between target and h[].                 */
    const float h[],                        /* i  : impulse response of weighted synthesis filter */
    float code[],                     /* o  : algebraic (fixed) codebook excitation         */
    float y[]                         /* o  : filtered fixed codebook excitation            */
);

void acelp_2t32(
    Encoder_State *st,                        /* i/o: encoder state structure                       */
    const float dn[],                       /* i  : corr. between target and h[].                 */
    const float h[],                        /* i  : impulse response of weighted synthesis filter */
    float code[],                     /* o  : algebraic (fixed) codebook excitation         */
    float y[]                         /* o  : filtered fixed codebook excitation            */
);

short acelp_4t64(
    Encoder_State *st,                        /* i/o: encoder state structure                       */
    float dn[],                       /* i : corr. between target and h[].                  */
    const float cn[],                       /* i : residual after long term prediction            */
    const float H[],                        /* i : impulse response of weighted synthesis filter  */
    float R[],                        /* i  : autocorrelation values                        */
    const short acelpautoc,                 /* i  : autocorrealtion flag                          */
    float code[],                     /* o : algebraic (fixed) codebook excitation          */
    float y[],                        /* o : filtered fixed codebook excitation             */
    short nbbits,                     /* i : number of bits per codebook                    */
    const short cmpl_flag,                  /* i  : coomplexity reduction flag                    */
    const short Opt_AMR_WB                  /* i  : flag indicating AMR-WB IO mode                */
);

void corr_xh(
    const float *x,                         /* i  : target signal                                 */
    float *y,                         /* o  : correlation between x[] and h[]               */
    const float *h,                         /* i  : impulse response of weighted synthesis filter */
    const int   L_subfr                     /* i  : length of the subframe                        */
);

void find_tilt(
    const float fr_bands[],                 /* i  : energy in frequency bands                  */
    const float bckr[],                     /* i  : per band background noise energy estimate  */
    float ee[2],                      /* o  : lf/hf E ration for present frame           */
    const short pitch[3],                   /* i  : open loop pitch values for 3 half-frames   */
    const float voicing[3],                 /* i  : normalized correlation for 3 half-frames   */
    const float *lf_E,                      /* i  : per bin energy  for low frequencies        */
    const float corr_shift,                 /* i  : normalized correlation correction          */
    const short bwidth,                     /* i  : input signal bandwidth                     */
    const short max_band,                   /* i  : maximum critical band                      */
    float hp_E[],                     /* o  : energy in HF                               */
    const short codec_mode,                 /* i  : Mode 1 or 2                                */
    float *bckr_tilt_lt,              /* i/o: lf/hf E ratio of background noise          */
    short Opt_vbr_mode
);

void init_gp_clip(
    float mem[]                       /* o  : memory of gain of pitch clipping algorithm */
);

short gp_clip(
    const long  core_brate,                 /* i  : core bitrate                               */
    const float *voicing,                   /* i  : normalized correlations (from OL pitch)    */
    const short i_subfr,                    /* i  : subframe index                             */
    const short coder_type,                 /* i  : coding type                                */
    const float xn[],                       /* i  : target vector                              */
    float mem[]                       /* i/o: memory of gain of pitch clipping algorithm */
);

void gp_clip_test_lsf(
    const long  core_brate,                 /* i  : core bitrate                               */
    const float lsf[],                      /* i  : LSF vector                                 */
    float mem[],                      /* i/o: memory of gain of pitch clipping algorithm */
    const short Opt_AMR_WB                  /* i  : flag indicating AMR-WB IO mode             */
);

void gp_clip_test_gain_pit(
    const long  core_brate,                 /* i  : core bitrate                                 */
    const float gain_pit,                   /* i  :   gain of quantized pitch                    */
    float mem[]                       /* i/o: memory of gain of pitch clipping algorithm   */
);

void analy_lp(
    const float speech[],                   /* i  : pointer to the denoised speech frame                */
    const short L_frame,                    /* i  : length of the frame                                 */
    const short L_look,                     /* i  : look-ahead length                                   */
    float *ener,                      /* o  : residual signal energy                              */
    float A[],                        /* o  : A(z) filter coefficients                            */
    float epsP[],                     /* o  : LP analysis residual energies for each iteration    */
    float lsp_new[],                  /* o  : current frame ISPs                                  */
    float lsp_mid[],                  /* o  : current mid-frame ISPs                              */
    float lsp_old[],                  /* i/o: previous frame unquantized ISPs                     */
    const short Top[2],                     /* i  : open loop pitch lag                                 */
    const float Tnc[2],                     /* i  : open loop pitch gain                                */
    const float sr_core                     /* i  : internal sampling rate                              */
);

void analy_lp_AMR_WB(
    const float speech[],                   /* i  : pointer to the speech frame                      */
    float *ener,                      /* o  : residual energy from Levinson-Durbin             */
    float A[],                        /* o  : A(z) filter coefficients                         */
    float epsP[],                     /* o  : LP analysis residual energies for each iteration */
    float isp_new[],                  /* o  : current frame ISPs                               */
    float isp_old[],                  /* i/o: previous frame unquantized ISPs                  */
    float isf_new[],                  /* o  : current frame ISFs                               */
    const int   Top,                        /* i  : open loop pitch lag                              */
    const float Tnc                         /* i  : open loop pitch gain                             */
);

void noise_est_init(
    float *totalNoise,                /* o  : noise estimate over all critical bands     */
    short *first_noise_updt,          /* o  : noise update initialization flag           */
    float bckr[],                     /* o  : per band background noise energy estimate  */
    float enrO[],                     /* o  : per band old input energy                  */
    float ave_enr[],                  /* o  : per band long-term average energies        */
    short *pitO,                      /* o  : open-loop pitch values from preceed. frame */
    short *aEn,                       /* o  : noise adaptation hangover counter          */
    short *st_harm_cor_cnt,           /* o  : harm and correlation timer                 */
    short *bg_cnt,                    /* o  : pause length counter                       */
    float *lt_tn_track,
    float *lt_tn_dist,
    float *lt_Ellp_dist,
    float *lt_haco_ev,
    short *low_tn_track_cnt
    ,float *Etot_st_est,
    float *Etot_sq_st_est

);

void long_enr(
    Encoder_State *st,                        /* i/o: encoder state structure                  */
    const float Etot,                       /* i  : total channel energy                     */
    const short localVAD_HE_SAD,            /* i  : HE-SAD flag without hangover             */
    short high_lpn_flag
);

void noise_est_pre(
    const float Etot,                       /* i  : Energy of current frame                  */
    const short ini_frame,                  /* i  : Frame number (init)                      */
    float *Etot_l,                    /* i/o: Track energy from below                  */
    float *Etot_h,                    /* i/o: Track energy from above                  */
    float *Etot_l_lp,                 /* i/o: Smoothed low energy                      */
    float *Etot_last,                 /* i/o: Energy of last frame                     */
    float *Etot_v_h2,                 /* i/o: Energy variations                        */
    float *sign_dyn_lp,               /* i/o: Smoother signal dynamics                 */
    short harm_cor_cnt,
    float *Etot_lp
);

void noise_est_down(
    const float fr_bands[],                 /* i  : per band input energy (contains 2 vectors)  */
    float bckr[],                     /* i/o: per band background noise energy estimate   */
    float tmpN[],                     /* o  : temporary noise update                      */
    float enr[],                      /* o  : averaged energy over both subframes         */
    const short min_band,                   /* i  : minimum critical band                       */
    const short max_band,                   /* i  : maximum critical band                       */
    float *totalNoise,                /* o  : noise estimate over all critical bands      */
    const float Etot,                       /* i  : Energy of current frame                     */
    float *Etot_last,                 /* i/o: Energy of last frame                        */
    float *Etot_v_h2                  /* i/o: Energy variaions of noise frames            */
);

void noise_est(
    Encoder_State *st,                        /* i/o: encoder state structure                                */
    const float tmpN[],                     /* i  : temporary noise update                                 */
    const short *pitch,                     /* i  : open-loop pitch values for each half-frame             */
    const float *voicing,                   /* i  : normalized correlation for all half-frames             */
    const float *epsP,                      /* i  : LP prediction error energies                           */
    const float Etot,                       /* i  : total channel E                                        */
    const float relE,                       /* i  : relative frame energy                                  */
    const float corr_shift,                 /* i  : normalized correlation correction                      */
    const float enr[],                      /* i  : averaged energy over both subframes                    */
    float fr_bands[],                 /* i  : spectrum per critical bands of the current frame       */
    float *cor_map_sum,               /* o  : sum of correlation map from mult-harm analysis         */
    float *sp_div,                    /* o  : soectral diversity feature                             */
    float *non_staX,                  /* o  : non-stationarity for sp/mus classifier                 */
    short *loc_harm,                  /* o  : multi-harmonicity flag for UV classifier               */
    const float *lf_E,                      /* i  : per bin energy  for low frequencies                    */
    short *st_harm_cor_cnt,           /* i  : 1st harm correlation timer                             */
    const float Etot_l_lp,                  /* i  : Smoothed low energy                                    */
    float *sp_floor                   /* o  : noise floor estimate                                   */
);

void vad_param_updt(
    Encoder_State *st,                        /* i/o: encoder state structure                                */
    const short pitch[3],                   /* i  : open loop pitch lag for each half-frame                */
    const float voicing[3],                 /* i  : maximum normalized correlation for each half-frame     */
    const float corr_shift,                 /* i  : correlation shift                                      */
    const float A[]                         /* i  : A(z) unquantized for the 4 subframes                   */
);

short multi_harm(                           /* o  : frame multi-harmonicity (1-harmonic, 0-not)     */
    const float Bin_E[],                    /* i  : log energy spectrum of the current frame        */
    float old_S[],                    /* i/o: prev. log-energy spectrum w. subtracted floor   */
    float cor_map_LT[],               /* i/o: LT correlation map                              */
    float *multi_harm_limit,          /* i/o: multi harminic threshold                        */
    const long  total_brate,                /* i  : total bitrate                                   */
    const short bwidth,                     /* i  : input signal bandwidth                          */
    short *cor_strong_limit,          /* i/o: HF correlation indicator                        */
    float *st_mean_avr_dyn,           /* i/o: long term average dynamic                       */
    float *st_last_sw_dyn,            /* i/o: last dynamic                                    */
    float *cor_map_sum,
    float *sp_floor                   /* o  : noise floor estimate                            */
);

void lp_gain_updt(
    const short i_subfr,                    /* i  :  subframe number                                */
    const float gain_pit,                   /* i  : Decoded gain pitch                              */
    const float norm_gain_code,             /* i  : Normalised gain code                            */
    float *lp_gainp,                  /* i/o: LP-filtered pitch gain(FEC)                     */
    float *lp_gainc,                  /* i/o: LP-filtered code gain (FEC)                     */
    const short L_frame                     /* i  : length of the frame                             */
);

void enc_pit_exc(
    Encoder_State *st,                        /* i/o: state structure                                 */
    LPD_state *mem,                       /* i/o: encoder memories                                */
    const float *speech,                    /* i  : Input speech                                    */
    const float Aw[],                       /* i  : weighted A(z) unquantized for subframes         */
    const float *Aq,                        /* i  : 12k8 Lp coefficient                             */
    const float Es_pred,                    /* i  : predicted scaled innov. energy                  */
    const short *T_op,                      /* i  : open loop pitch                                 */
    const float *voicing,                   /* i  : voicing                                         */
    const float *res,                       /* i  : residual signal                                 */
    float *synth,                     /* i/o: core synthesis                                  */
    float *exc,                       /* i/o: current non-enhanced excitation                 */
    short *T0,                        /* i/o: close loop integer pitch                        */
    short *T0_frac,                   /* i/o: close-loop pitch period - fractional part       */
    float *pitch_buf,                 /* i/o: Fractionnal per subframe pitch                  */
    const short nb_subfr,                   /* i  : Number of subframe considered                   */
    float *gpit                       /* o  : pitch gain per subframe                         */
);

void encod_audio(
    Encoder_State *st,                        /* i/o: state structure                                  */
    LPD_state *mem,                       /* i/o: encoder memories                                 */
    const float speech[],                   /* i  : input speech                                     */
    const float Aw[],                       /* i  : weighted A(z) unquantized for subframes          */
    const float Aq[],                       /* i  : 12k8 Lp coefficient                              */
    const short T_op[],                     /* i  : open loop pitch                                  */
    const float voicing[],                  /* i  : voicing                                          */
    const float *res,                       /* i  : residual signal                                  */
    float *synth,                     /* i/o: core synthesis                                   */
    float *exc,                       /* i/o: current non-enhanced excitation                  */
    float *pitch_buf,                 /* i/o: floating pitch values for each subframe          */
    float *voice_factors,             /* o  : voicing factors                                  */
    float *bwe_exc,                   /* o  : excitation for SWB TBE                           */
    const short attack_flag,                /* i  : Flag that point to an attack coded with AC mode (GSC)    */
    const short coder_type,                 /* i  : core coding type                                 */
    float *lsf_new,                   /* i  : current frame ISF vector                         */
    float *tmp_noise                  /* o  : long-term noise energy                           */
);

short Pit_exc_contribution_len(             /* o  : index of the last band where pitch contribution is significant */
    Encoder_State *st,                        /* i/o: state structure                                 */
    const float *dct_res,                   /* i  : DCT of residual                                 */
    float *dct_pitex,                 /* i/o: DCT of pitch contribution                       */
    float *pitch_buf,                 /* i/o: Pitch per subframe                              */
    short *hangover,                  /* i  : Hangover for the time contribution switching    */
    const short coder_type                  /* i  : coding type                                     */
);

short stab_est(
    float etot,                       /* i  : Total energy of the current frame   */
    float *lt_diff_etot,              /* i/o: Long term total energy variation    */
    float *mem_etot,                  /* i/o: Total energy memory                 */
    short *nb_thr_3,                  /* i/o: Number of consecutives frames of level 3 */
    short *nb_thr_1,                  /* i/o: Number of consecutives frames of level 1 */
    float *thresh,                    /* i/o: Detection thresold                 */
    short *last_music_flag,           /* i/o: Previous music detection ouptut    */
    short vad_flag
);

float gsc_gainQ(
    Encoder_State *st,                         /* i/o: encoder state structure           */
    const float y_gain4[],                   /* i  : gain per band                     */
    float y_gainQ[],                   /* o  : quantized gain per band           */
    const long  core_brate,                  /* i  : Core rate                         */
    const short coder_type,                  /* i  : coding type                       */
    const short bwidth                       /* i  : input signal bandwidth            */
);

void Ener_per_band_comp(
    const float exc_diff[],                  /* i  : gain per band                     */
    float y_gain4[],                   /* o  : gain per band to quantize         */
    const short Mband,                       /* i  : Max band                          */
    const short Eflag                        /* i  : flag of highest band              */
);

void Comp_and_apply_gain(
    float exc_diffQ[],                 /* i/o: gain per band                     */
    float Ener_per_bd_iQ[],            /* o  : Quant Ener per band               */
    float Ener_per_bd_yQ[],            /* o  : Ener per band for quantize y      */
    short Mbands_gn,                   /* i  : number of bands                   */
    const short ReUseGain                    /* i  : Reuse the gain in Ener_per_bd_yQ  */
);

void bands_and_bit_alloc(
    const short cor_strong_limit,             /* i  : HF correlation                                        */
    const short noise_lev,                    /* i  : dwn scaling factor                                    */
    const long  core_brate,                   /* i  : core bit rate                                         */
    const short Diff_len,                     /* i  : Lenght of the difference signal (before pure spectral)*/
    const short bits_used,                    /* i  : Number of bit used before frequency Q                 */
    short *bit,                         /* i/o: Number of bit allowed for frequency quantization      */
    float *ener_vec,                    /* i/o: Quantized energy vector                               */
    short *max_ener_band,               /* o  : Sorted order                                          */
    short *bits_per_bands_s,            /* i/o: Number of bit allowed per allowed subband (Q3)        */
    short *nb_subbands,                 /* o  : Number of subband allowed                             */
    const float *exc_diff,                    /* i  : Difference signal to quantize (encoder side only)     */
    float *concat_in,                   /* o  : Concatened PVQ's input vector (encoder side only)     */
    short *pvq_len,                     /* o  : Number of bin covered with the PVQ                    */
    const short coder_type,                   /* i  : coding type                                           */
    const short bwidth,                       /* i  : input signal bandwidth                                */
    const short GSC_noisy_speech              /* i  : GSC noisy speech flag                                 */
);

float gsc_gaindec(                           /* o  : average frequency gain                  */
    Decoder_State *st,                          /* i/o: decoder state structure                 */
    float y_gainQ[],                    /* o  : quantized gain per band                 */
    const long  core_brate,                   /* i  : core used                               */
    float old_y_gain[],                 /* i/o: AR gain quantizer for low rate          */
    const short coder_type,                   /* i  : coding type                             */
    const short bwidth                        /* i  : input signal bandwidth                  */
);

void freq_dnw_scaling(
    const short cor_strong_limit,             /* i  : HF correlation                          */
    const short coder_type,                   /* i  : coder type                              */
    const short noise_lev,                    /* i  : Noise level                             */
    const long  core_brate,                   /* i  : Core bitrate                            */
    float fy_norm[]                     /* i/o: Frequency quantized parameter           */
);

void decod_audio(
    Decoder_State *st,                        /* i/o: decoder static memory                     */
    float dct_epit[],                 /* o  : GSC excitation in DCT domain              */
    const float *Aq,                        /* i  : LP filter coefficient                     */
    const short coder_type,                 /* i  : coding type                               */
    float *tmp_noise,                 /* o  : long term temporary noise energy          */
    float *pitch_buf,                 /* o  : floating pitch values for each subframe   */
    float *voice_factors,             /* o  : voicing factors                           */
    float *exc_dct_in,                /* i/o: adapt. excitation exc                     */
    float *exc2,                      /* i/o: adapt. excitation/total exc               */
    float *bwe_exc,                   /* o  : excitation for SWB TBE                    */
    float *lsf_new,                   /* i  : current frame ISF vector                  */
    float *gain_buf
);

void gsc_dec(
    Decoder_State *st,                        /* i/o: State structure                           */
    float exc_dct_in[],               /* i/o: dct of pitch-only/total excitation        */
    const short pit_band_idx,               /* i  : pitch band index                          */
    const short Diff_len,                   /* i  :  */
    const short bits_used,                  /* i  : total number of bits used                 */
    const short nb_subfr,                   /* i  : Number of subframe considered             */
    const short coder_type,                 /* i  : coding type                               */
    short *last_bin,                  /* i  : last bin of bit allocation                */
    float *lsf_new,                   /* i  : ISFs at the end of the frame              */
    float *exc_wo_nf,                 /* o  : excitation (in f domain) without noisefill*/
    float *tmp_noise                  /* o  : long-term noise energy                    */
);

void dec_pit_exc(
    Decoder_State *st,                        /* i/o: decoder static memory                     */
    const short L_frame,                    /* i  : length of the frame                       */
    const float *Aq,                        /* i  : LP filter coefficient                     */
    const float Es_pred,                    /* i  : predicted scaled innov. energy            */
    float *pitch_buf,                 /* o  : floating pitch values for each subframe   */
    float *code,                      /* o  : innovation                                */
    float *exc,                       /* i/o: adapt. excitation exc                     */
    const short nb_subfr,                   /* i  : Number of subframe considered             */
    float *gain_buf
);

void highband_exc_dct_in(
    const long core_brate,                 /* i  : core bitrate                            */
    const short *mfreq_bindiv_loc,         /* i  : bin per bands tables                    */
    short last_bin,                  /* i  : last bin of bit allocation              */
    short Diff_len,                  /* i  : number of bin before cut-off frequency  */
    short noise_lev,                 /* i  : pulses dynamic                          */
    short pit_band_idx,              /* i  : bin position of the cut-off frequency   */
    float *exc_diffQ,                /* i  : frequency coefficients of per band      */
    short *seed_tcx,                 /* i  : Seed for noise                          */
    float *Ener_per_bd_iQ,           /* i  : Quantized energy of targeted vector     */
    short nb_subfr,                  /* i  : Number of subframe considered           */
    float *exc_dct_in,               /* o  : dct of residual signal                  */
    short last_coder_type,           /* i  : coding type of last frame               */
    short *bitallocation_band,       /* i  : bit allocation flag of each band        */
    float *lsf_new,                  /* i  : ISFs at the end of the frame            */
    float *last_exc_dct_in,          /* i  : dct of residual signal of last frame    */
    float *last_ener,                /* i  : frequency energy  of last frame         */
    short *last_bitallocation_band,  /* i  : bit allocation flag of each band  of last frame   */
    short *bitallocation_exc,        /* i  : flag of decoded coefficients            */
    short bfi,                       /* i  : bad frame indicator                     */
    const short coder_type,                /* i  : coder type                              */
    short bwidth,
    float *exc_wo_nf,                /* o  : excitation (in f domain) without noisefill   */
    const short GSC_noisy_speech
    ,float *lt_ener_per_band_fx       /* i/o: Average per band energy */
);

void inact_switch_ematch(
    float exc2[],                     /* i/o: CELP/GSC excitation buffer                  */
    float dct_exc_tmp[],              /* i  : GSC excitation in DCT domain                */
    float lt_ener_per_band[],         /* i/o: long-term energy per band                   */
    const short coder_type,                 /* i  : coding mode                                 */
    const short L_frame,                    /* i  : frame lenght                                */
    const long  core_brate                  /* i : core bit rate                                */
    ,const short bfi,                        /* i  : frame lost indicator                        */
    const short last_core,                  /* i  : Last core used                              */
    const short last_codec_mode             /* i  : Last codec mode                             */
);

void LD_music_post_filter(
    const float dtc_in[],                   /* i  : input synthesis                             */
    float dtc_out[],                  /* o  : output synthesis                            */
    const  long core_brate,                 /* i  : core bitrate                                */
    short *last_music_flag,           /* i/o: Previous music detection ouptut             */
    float *thresh,                    /* i/o: Detection thresold                          */
    short *nb_thr_1,                  /* i/o: Number of consecutives frames of level 1    */
    short *nb_thr_3,                  /* i/o: Number of consecutives frames of level 3    */
    float *lt_diff_etot,              /* i/o: Long term total energy variation            */
    float *mem_etot,                  /* i/o: Total energy memory                         */
    const float min_ns_gain,                /* i  : minimum gain for inter-harm noise red.      */
    float bckr[],                     /* i/o: per band bckgnd. noise energy estimate      */
    float lf_EO[],                    /* i/o: old per bin E for previous half frame       */
    float lp_gbin[],                  /* i/o: smoothed suppression gain, per FFT bin      */
    float *filt_lfE,                  /* i  : post filter weighting coefficient           */
    short *last_nonfull_music,        /* i  : Coder type : -1 in case of IO               */
    const short coder_type                  /* i  : Coder type : -1 in case of IO               */
    ,const short Last_coder_type             /* i   : last Coder type                            */
);

void Post_music_postP(
    float dct_buffer_in[],            /* i/o: excitation buffer                           */
    float exc_buffer_out[],           /* o  : DCT output buffer                           */
    float *exc2,                      /* i/o: Current excitation to be overwriten         */
    const float *mem_tmp,                   /* i  : previous frame synthesis memory             */
    float *st_mem_syn2,               /* i/o: current frame synthesis memory              */
    const float *Aq,                        /* i  : LPC filter coefficients                     */
    float *syn                        /* i/o: 12k8 synthesis                              */
);

void Prep_music_postP(
    float exc_buffer_in[],            /* i/o: excitation buffer                           */
    float dct_buffer_out[],           /* o  : DCT output buffer                           */
    float filt_lfE[],                 /* i/o: long term spectrum energy                   */
    const short last_core,                  /* i  : last core                                   */
    const float *pitch_buf,                 /* i  : current frame pitch information             */
    float *LDm_enh_lp_gbin            /* o  : smoothed suppression gain, per bin FFT      */
);

void speech_music_classif(
    Encoder_State *st,                        /* i/o: encoder state structure                         */
    short *sp_aud_decision0,          /* o  : 1st stage speech/music decision                 */
    short *sp_aud_decision1,          /* o  : 1st stage speech/music decision                 */
    short *sp_aud_decision2,          /* o  : 2nd stage speech/music decision                 */
    const float *new_inp,                   /* i  : new input signal                                */
    const float *inp,                       /* i  : input signal to locate attach position          */
    const short vad_flag,
    const short localVAD,
    const short localVAD_HE_SAD,            /* i  : HE-SAD flag without hangover                    */
    const short pitch[3],                   /* i  : open-loop pitch estimate in three subframes     */
    const float voicing[3],                 /* i  : voicing estimate in three subframes             */
    const float lsp_new[M],                 /* i  : LSPs in current frame                           */
    const float cor_map_sum,                /* i  : correlation map sum (from multi-harmonic anal.) */
    const float epsP[M+1],                  /* i  : LP prediciton error                             */
    const float PS[],                       /* i  : energy spectrum                                 */
    const float Etot,                       /* i  : total frame energy                              */
    const float old_cor,                    /* i  : max correlation from previous frame             */
    short *coder_type,                /* i/o: coding type                                     */
    short *attack_flag,               /* o  : flag to indicate if attack is to be treated by TC or GSC */
    const float non_staX,                   /* i  : unbound non-stationarity for sp/mus classifier  */
    const float relE,                       /* i  : relative frame energy                           */
    short *high_lpn_flag,
    const short flag_spitch                 /* i  : flag to indicate very short stable pitch        */
);

void find_wsp(
    const short L_frame,                    /* i  : length of the frame                         */
    const short L_subfr,                    /* i  : length of subframe                          */
    const short nb_subfr,                   /* i  : number of subframes                         */
    const float *A,                         /* i  : A(z) filter coefficients                    */
    float *Aw,                        /* o  : weighted A(z) filter coefficients           */
    const float *speech,                    /* i  : pointer to the denoised speech frame        */
    const float tilt_fact,                  /* i  : tilt factor                                 */
    float *wsp,                       /* o  : poitnter to the weighted speech frame       */
    float *mem_wsp,                   /* i/o: W(Z) denominator memory                     */
    const float gamma,                      /* i  : weighting factor                            */
    const short L_look                      /* i  : look-ahead                                  */
);

void pitch_ol_init(
    float *old_thres,                 /* o  : threshold for reinforcement of past pitch influence */
    short *old_pitch,                 /* o  : pitch  of the 2nd half-frame of previous frame      */
    short *delta_pit,                 /* o  : pitch evolution extrapolation                       */
    float *old_corr                   /* o  : correlation                                         */
);

void pitch_ol(
    short pitch[3],                   /* o  : open loop pitch lag for each half-frame                        */
    float voicing[3],                 /* o  : maximum normalized correlation for each half-frame             */
    short *old_pitch,                 /* i/o: OL pitch of the 2nd half-frame of the last frame               */
    float *old_corr,                  /* i/o: correlation                                                    */
    float corr_shift,                 /* i  : normalized correlation correction                              */
    float *old_thres,                 /* i/o: maximum correlation weighting with respect to past frame pitch */
    short *delta_pit,                 /* i/o: old pitch extrapolation correction (added to old pitch)        */
    float *st_old_wsp2,               /* i/o: weighted speech memory                                         */
    const float *wsp,                       /* i  : weighted speech for current frame and look-ahead               */
    float mem_decim2[3],              /* i/o: wsp decimation filter memory                                   */
    const float relE,                       /* i  : relative frame energy                                          */
    const short L_look,                     /* i  : look-ahead                                                     */
    const short last_class,                 /* i  : frame classification of last frame                             */
    const short bwidth,                     /* i  : bandwidth                                                      */
    const short Opt_SC_VBR                  /* i : SC-VBR flag                                                     */
);

void pitch_ol2(
    const short pit_min,                    /* i  : pit_min value                                                   */
    const short pitch_ol,                   /* i  : pitch to be improved                                            */
    float *pitch_fr,                  /* o  : adjusted 1/4 fractional pitch                                   */
    float *voicing_fr,                /* o  : adjusted 1/4 fractional voicing                                 */
    const short pos,                        /* i  : position in frame where to calculate the improv.                */
    const float *wsp,                       /* i  : weighted speech for current frame and look-ahead                */
    const short delta                       /* i  : delta for pitch search                                          */
);

void StableHighPitchDetect(
    short *flag_spitch,               /* o  : flag to indicate very short stable pitch */
    short pitch[],                    /* i/o: OL pitch buffer                         */
    const float voicing[],                  /* i  : OL pitch gains                          */
    const float Bin_E[],                    /* i  : per bin log energy spectrum                                     */
    const float wsp[],                      /* i  : weighted speech                                                 */
    const short localVAD,
    float *voicing_sm,                /* i/o: smoothed open-loop pitch gains                                  */
    float *voicing0_sm,               /* i/o: smoothed high pitch gains                                       */
    float *LF_EnergyRatio_sm,         /* i/o: smoothed [0, 300Hz] relative peak energy                        */
    short *predecision_flag,          /* i/o: predecision flag                                                */
    float *diff_sm,                   /* i/o: smoothed pitch frequency difference                             */
    float *energy_sm                  /* i/o: smoothed energy around pitch frequency                          */
);

void pitchDoubling_det(
    float *wspeech,
    short *pitch_ol,
    float *pitch_fr,
    float *voicing_fr
);

void gain_enc_amr_wb(
    Encoder_State *st,                        /* i/o: encoder state structure                                         */
    const float *xn,                        /* i  : target vector                                                   */
    const float *y1,                        /* i  : zero-memory filtered adaptive excitation                        */
    const float *y2,                        /* i  : zero-memory filtered algebraic codebook excitation              */
    const float *code,                      /* i  : algebraic excitation                                            */
    const long  core_brate,                 /* i  : core bitrate                                                    */
    float *gain_pit,                  /* i/o: Pitch gain / Quantized pitch gain                               */
    float *gain_code,                 /* o  : Quantized codebook gain                                         */
    float *gain_inov,                 /* o  : innovation gain                                                 */
    float *norm_gain_code,            /* o  : norm. gain of the codebook excitation                           */
    float *coeff,                     /* i/o: correlations <y1,y1>, -2<xn,y1>,<y2,y2>, -2<xn,y2> and 2<y1,y2> */
    const short clip_gain,                  /* i  : gain pitch clipping flag (1 = clipping)                         */
    float *past_qua_en                /* i/o: gain quantization memory (4 words)                              */

);

void gain_enc_lbr(
    Encoder_State *st,                        /* i/o: encoder state structure                                         */
    const long  core_brate,                 /* i  : core bitrate                                                    */
    const short coder_type,                 /* i  : coding type                                                     */
    const short i_subfr,                    /* i  : subframe index                                                  */
    const float *xn,                        /* i  : target vector                                                   */
    const float *y1,                        /* i  : zero-memory filtered adaptive excitation                        */
    const float *y2,                        /* i  : zero-memory filtered algebraic codebook excitation              */
    const float *code,                      /* i  : algebraic excitation                                            */
    float *gain_pit,                  /* o  : quantized pitch gain                                            */
    float *gain_code,                 /* o  : quantized codebook gain                                         */
    float *gain_inov,                 /* o  : gain of the innovation (used for normalization)                 */
    float *norm_gain_code,            /* o  : norm. gain of the codebook excitation                           */
    float *g_corr,                    /* i/o: correlations <y1,y1>, -2<xn,y1>,<y2,y2>, -2<xn,y2> and 2<y1,y2> */
    float gains_mem[],                /* i/o: pitch gain and code gain from previous subframes                */
    const short clip_gain                   /* i  : gain pitch clipping flag (1 = clipping)                         */
);

void gain_enc_mless(
    Encoder_State *st,                        /* i/o: encoder state structure                                         */
    const long  core_brate,                 /* i  : core bitrate                                                    */
    const short L_frame,                    /* i  : length of the frame                                             */
    const short coder_type,                 /* i  : coding type                                                     */
    const short i_subfr,                    /* i  : subframe index                                                  */
    const short tc_subfr,                   /* i  : TC subframe index                                               */
    const float *xn,                        /* i  : target vector                                                   */
    const float *y1,                        /* i  : zero-memory filtered adaptive excitation                        */
    const float *y2,                        /* i  : zero-memory filtered algebraic codebook excitation              */
    const float *code,                      /* i  : algebraic excitation                                            */
    const float Es_pred,                    /* i  : predicted scaled innovation energy                              */
    float *gain_pit,                  /* o  : quantized pitch gain                                            */
    float *gain_code,                 /* o  : quantized codebook gain                                         */
    float *gain_inov,                 /* o  : innovation gain                                                 */
    float *norm_gain_code,            /* o  : norm. gain of the codebook excitation                           */
    float *coeff,                     /* i/o: correlations <y1,y1>, -2<xn,y1>,<y2,y2>, -2<xn,y2> and 2<y1,y2> */
    const short clip_gain                   /* i  : gain pitch clipping flag (1 = clipping)                         */
);

void gain_enc_SQ(
    Encoder_State *st,                       /* i/o: encoder state structure                                         */
    const long  core_brate,                /* i  : core bitrate                                                    */
    const short coder_type,                /* i  : coding type                                                     */
    const short i_subfr,                   /* i  : subframe index                                                  */
    const short tc_subfr,                  /* i  : TC subframe index                                               */
    const float *xn,                       /* i  : target vector                                                   */
    const float *yy1,                      /* i  : zero-memory filtered adaptive excitation                        */
    const float *y2,                       /* i  : zero-memory filtered algebraic codebook excitation              */
    const float *code,                     /* i  : algebraic excitation                                            */
    const float Es_pred,                   /* i  : predicted scaled innovation energy                              */
    float *gain_pit,                 /* o  : quantized pitch gain                                            */
    float *gain_code,                /* o  : quantized codebook gain                                         */
    float *gain_inov,                /* o  : gain of the innovation (used for normalization)                 */
    float *norm_gain_code,           /* o  : norm. gain of the codebook excitation                           */
    float *g_corr,                   /* i/o: correlations <y1,y1>, -2<xn,y1>,<y2,y2>, -2<xn,y2> and 2<y1,y2> */
    const short clip_gain                  /* i  : gain pitch clipping flag (1 = clipping)                         */
);

short gain_enc_gaus(                        /* o  : Return index of quantization                */
    float *gain,                      /* i/o: Code gain to quantize                       */
    const short bits,                       /* i  : number of bits to quantize                  */
    const float lowBound,                   /* i  : lower bound of quantizer (dB)               */
    const float topBound                    /* i  : upper bound of quantizer (dB)               */
);

void E_corr_xy2(
    const float xn[],                       /* i  : target vector                               */
    const float y1[],                       /* i  : filtered excitation components 1            */
    const float y2[],                       /* i  : filtered excitation components 2            */
    float g_corr[],                   /* o  : correlations between x, y1, y2, y3, y4      */
    const short L_subfr                     /* i : subframe size                                */
);

float pit_encode(                           /* o  : Floating pitch for each subframe            */
    Encoder_State *st,                        /* i/o: encoder state structure                     */
    const long  core_brate,                 /* i  : core bitrate                                */
    const short Opt_AMR_WB,                 /* i  : flag indicating AMR-WB IO mode              */
    const short L_frame,                    /* i  : length of the frame                         */
    const short coder_type,                 /* i  : coding type                                 */
    short *limit_flag,                /* i/o: restrained(0) or extended(1) Q limits       */
    const short i_subfr,                    /* i  : subframe index                              */
    float *exc,                       /* i/o: pointer to excitation signal frame          */
    const short L_subfr,                    /* i  : subframe length                             */
    const short *T_op,                      /* i  : open loop pitch estimates in current frame  */
    short *T0_min,                    /* i/o: lower limit for close-loop search           */
    short *T0_max,                    /* i/o: higher limit for close-loop search          */
    short *T0,                        /* i/o: close loop integer pitch                    */
    short *T0_frac,                   /* i/o: close loop fractional part of the pitch     */
    const float *h1,                        /* i  : weighted filter input response              */
    const float *xn                         /* i  : target vector                               */
);

short find_uv(                              /* o  : coding type                                   */
    Encoder_State *st,                        /* i/o: encoder state structure                       */
    const float *pitch_fr,                  /* i  : pointer to adjusted fractional pitch (4 val.) */
    const float *voicing_fr,                /* i  : refined correlation for each subframes        */
    const float *voicing,                   /* i  : correlation for 3 half-frames                 */
    const float *speech,                    /* i  : pointer to speech signal for E computation    */
    const short localVAD,                   /* i  : vad without hangover                          */
    const float *ee,                        /* i  : lf/hf Energy ratio for present frame          */
    const float corr_shift,                 /* i  : normalized correlation correction in noise    */
    const float relE,                       /* i  : relative frame energy                         */
    const float Etot,                       /* i  : total energy                                  */
    const float hp_E[],                     /* i  : energy in HF                                  */
    short *flag_spitch,               /* i/o: flag to indicate very short stable pitch and high correlation */
    float voicing_sm                  /* i/o: smoothed open-loop pitch gains                */
    , const short last_core_orig              /* i  : original last core                            */
);

short signal_clas(                          /* o  : classification for current frames             */
    Encoder_State *st,                        /* i/o: encoder state structure                           */
    short *coder_type,                /* i/o: coder type                                        */
    const float voicing[3],                 /* i  : normalized correlation for 3 half-frames          */
    const float *speech,                    /* i  : pointer to speech signal for E computation        */
    const short localVAD,                   /* i  : vad without hangover                              */
    const short pit[3],                     /* i  : open loop pitch values for 3 half-frames          */
    const float *ee,                        /* i  : lf/hf E ration for 2 half-frames                  */
    const float relE,                       /* i  : frame relative E to the long term average         */
    const short L_look,                     /* i  : look-ahead                                        */
    short *uc_clas                    /* o  : flag for unvoiced class used in sp/mus classifier */
);

void wb_vad_init(
    short *nb_active_frames,          /* o  : nb of consecutive active speech frames                  */
    short *hangover_cnt,
    float *lp_speech,                 /* o  : long-term active speech level                           */
    short *nb_active_frames_he,       /* o  : nb of consecutive active speech frames                  */
    short *hangover_cnt_he,
    float *bcg_flux,                  /* o  : background noise fluctuation                            */
    short *soft_hangover,             /* o  : soft hangover counter                                   */
    short *voiced_burst,              /* o  : consecutive voiced speech counter                       */
    short *bcg_flux_init,             /* o  : initialization period for noise fluctuation estimation  */
    short *nb_active_frames_he1,      /* o  : nb of consecutive active speech frames 1                */
    short *hangover_cnt_he1,
    long  *vad_flag_reg_H,
    long  *vad_flag_reg_L,            /* o  :                                                      lower 31 bits */
    long  *vad_prim_reg,
    short *vad_flag_cnt_50,
    short *vad_prim_cnt_16,
    short *hangover_cnt_dtx,
    short *flag_noisy_speech_snr
    ,short *hangover_cnt_music         /* o  : counter of VAD DTX Music hangover frames                */
);

short dtx_hangover_addition(
    Encoder_State *st,                        /* i/o: encoder state structure                    */
    const short localVAD,                   /* i  : Primary vad decision                       */
    const short vad_flag,
    const float snr,                        /* i  : input single SNR estimate                  */
    const short cldfb_subtraction,          /* i  :                                            */
    short *vad_hover_flag
);

short wb_vad(
    Encoder_State *st,                        /* i/o: encoder state structure                    */
    const float fr_bands[],                 /* i  : per band input energy (contains 2 vectors) */
    short *localVAD,
    short *noisy_speech_HO,           /* o  : SC-VBR noisy speech HO flag                */
    short *clean_speech_HO,           /* o  : SC-VBR clean speech HO flag                */
    short *NB_speech_HO,              /* o  : SC-VBR NB speech HO flag                   */
    float *snr_sum_he,                /* i  : voicing metric from SAD                    */
    short *localVAD_HE_SAD,           /* o  : HE_SAD decision without hangovers          */
    short *flag_noisy_speech_snr      /* o  :  */
);

void bw_detect(
    Encoder_State *st,                        /* i/o: Encoder State                               */
    const float signal_in[],                /* i  : input signal                                */
    const short localVAD,
    float *enerBuffer
);

float gaus_encode(
    Encoder_State *st,                        /* i/o: encoder state structure                     */
    const short i_subfr,                    /* i  : subframe index                              */
    const float *h1,                        /* i  : weighted filter input response              */
    const float *xn,                        /* i  : target vector                               */
    float *exc,                       /* o  : pointer to excitation signal frame          */
    float *mem_w0,                    /* o  : weighting filter denominator memory         */
    float *gp_clip_mem,               /* o  : memory of gain of pitch clipping algorithm  */
    float *tilt_code,                 /* o  : synthesis excitation spectrum tilt          */
    float *code,                      /* o  : algebraic excitation                        */
    float *gain_code,                 /* o  : Code gain.                                  */
    float *y2,                        /* o  : zero-memory filtered adaptive excitation    */
    float *gain_inov,                 /* o  : innovation gain                             */
    float *voice_fac,                 /* o  : voicing factor                              */
    float *gain_pit,                  /* o  : adaptive excitation gain                    */
    float *norm_gain_code,            /* o  : normalized innovative cb. gain              */
    const long  core_brate                  /* i  : core bitrate                                */
);

void dtx(
    Encoder_State *st,                        /* i/o: encoder state structure                     */
    const short vad,                        /* i  : vad flag                                    */
    const float speech[]                    /* i  : Pointer to the speech frame                 */
);

void dtx_hangover_control(
    Encoder_State *st,                      /* i/o: encoder state structure                     */
    const float lsp_new[M]                /* i  : current frame LSPs                          */
);

void updt_enc(
    Encoder_State *st,                        /* i/o: state structure                             */
    const short L_frame,                    /* i  : length of the frame                         */
    const short coder_type,                 /* i  : coding type                                 */
    const float *old_exc,                   /* i  : buffer of excitation                        */
    const float *pitch_buf,                 /* i  : Floating pitch   for each subframe          */
    const float Es_pred,                    /* i  : predicited scaled innovation energy         */
    const float *Aq,                        /* i  : A(z) quantized for all subframes            */
    const float *lsf_new,                   /* i  : current frame LSF vector                    */
    const float *lsp_new,                   /* i  : current frame LSP vector                    */
    const float *old_bwe_exc                /* o  : buffer of excitation for SWB TBE            */
);

void updt_IO_switch_enc(
    Encoder_State *st,                      /* i/o: state structure                             */
    const short input_frame                 /* i  : input frame length                          */
);

void transition_enc(
    Encoder_State *st,                        /* i/o: encoder state structure                     */
    const long  core_brate,                 /* i  : core bitrate                                */
    const short L_frame,                    /* i  : length of the frame                         */
    const short coder_type,                 /* i  : coding type                                 */
    const short i_subfr,                    /* i  : subframe index                              */
    short *tc_subfr,                  /* i/o: TC subframe index                           */
    short *Jopt_flag,                 /* i  : joint optimization flag                     */
    short *position,                  /* i/o: maximum of residual signal index            */
    const float voicing[],                  /* i  : normalized correlations (from OL pitch)     */
    const short T_op[],                     /* i  : open loop pitch estimates in current frame  */
    short *T0,                        /* i/o: close loop integer pitch                    */
    short *T0_frac,                   /* i/o: close loop fractional part of the pitch     */
    short *T0_min,                    /* i/o: lower limit for close-loop search           */
    short *T0_max,                    /* i/o: higher limit for close-loop search          */
    float *exc,                       /* i/o: pointer to excitation signal frame          */
    float *y1,                        /* o  : zero-memory filtered adaptive excitation    */
    const float *res,                       /* i  : pointer to the LP residual signal frame     */
    const float *h1,                        /* i  : weighted filter input response              */
    const float *xn,                        /* i  : target vector                               */
    float *xn2,                       /* o  : target vector for innovation search         */
    float *gp_cl,                     /* i/o: memory of gain of pitch clipping algorithm  */
    float *gain_pit,                  /* o  : adaptive excitation gain                    */
    float *g_corr,                    /* o  : ACELP correlation values                    */
    short *clip_gain,                 /* i/o: adaptive gain clipping flag                 */
    float **pt_pitch,                 /* o:   floating pitch values                       */
    float *bwe_exc                    /* i/o: excitation for SWB TBE                      */
);

void set_impulse(
    const float xn[],                       /* i  : target signal                                   */
    const float h_orig[],                   /* i  : impulse response of weighted synthesis filter   */
    float exc[],                      /* o  : adaptive codebook excitation                    */
    float y1[],                       /* o  : filtered adaptive codebook excitation           */
    short *imp_shape,                 /* o  : adaptive codebook index                         */
    short *imp_pos,                   /* o  : position of the glotal impulse center index     */
    float *gain_trans                 /* o  : transition gain                                 */
);

void gain_enc_tc(
    Encoder_State *st,                        /* i/o: encoder state structure                                 */
    const long  core_brate,                 /* i  : core bitrate                                            */
    const short L_frame,                    /* i  : length of the frame                                     */
    const short i_subfr,                    /* i  : subframe index                                          */
    const short tc_subfr,                   /* i  : TC subframe index                                       */
    const float xn[],                       /* i  : target vector                                           */
    const float y2[],                       /* i  : zero-memory filtered algebraic codebook excitation      */
    const float code[],                     /* i  : algebraic excitation                                    */
    const float Es_pred,                    /* i  : predicted scaled innovation energy                      */
    float *gain_pit,                  /* o  : pitch gain / Quantized pitch gain                       */
    float *gain_code,                 /* o  : quantized codebook gain                                 */
    float *gain_inov,                 /* o  : innovation gain                                         */
    float *norm_gain_code             /* o  : norm. gain of the codebook excitation                   */
);

float corr_xy1(                             /* o  : pitch gain  (0..GAIN_PIT_MAX)               */
    const float xn[],                       /* i  : target signal                               */
    const float y1[],                       /* i  : filtered adaptive codebook excitation       */
    float g_corr[],                   /* o  : correlations <y1,y1>  and -2<xn,y1>         */
    const short L_subfr,
    const short norm_flag                   /* i : flag for constraining pitch contribution     */
);

void norm_corr(
    const float exc[],                      /* i  : excitation buffer                          */
    const float xn[],                       /* i  : target signal                              */
    const float h[],                        /* i  : weighted synthesis filter impulse response */
    const short t_min,                      /* i  : minimum value of searched range            */
    const short t_max,                      /* i  : maximum value of searched range            */
    float corr_norm[],                /* o  : normalized correlation                     */
    const short L_subfr                     /* i  : subframe size                              */
);

short pitch_fr4(                            /* o  : chosen integer pitch lag                   */
    const float exc[],                      /* i  : excitation buffer                          */
    const float xn[],                       /* i  : target signal                              */
    const float h[],                        /* i  : weighted synthesis filter impulse response */
    const short t0_min,                     /* i  : minimum value in the searched range.       */
    const short t0_max,                     /* i  : maximum value in the searched range.       */
    short *pit_frac,                  /* o  : chosen fraction (0, 1, 2 or 3)             */
    const short i_subfr,                    /* i  : flag to first subframe                     */
    const short limit_flag,                 /* i  : flag for limits (0=restrained, 1=extended) */
    const short t0_fr2,                     /* i  : minimum value for resolution 1/2           */
    const short t0_fr1,                     /* i  : minimum value for resolution 1             */
    const short L_frame,                    /* i  : length of the frame                        */
    const short L_subfr                     /* i  : size of subframe                           */
);

void pit_Q_enc(
    Encoder_State *st,                        /* i/o: encoder state structure                 */
    const short Opt_AMR_WB,                 /* i  : flag indicating AMR-WB IO mode          */
    const short nBits,                      /* i  : # of Q bits                             */
    const short delta,                      /* i  : Half the CL searched interval           */
    const short pit_flag,                   /* i  : absolute(0) or delta(1) pitch Q         */
    const short limit_flag,                 /* i  : restrained(0) or extended(1) Q limits   */
    const short T0,                         /* i  : integer pitch lag                       */
    const short T0_frac,                    /* i  : pitch fraction                          */
    short *T0_min,                    /* i/o: delta search min                        */
    short *T0_max                     /* o  : delta search max                        */
);

void pit16k_Q_enc(
    Encoder_State *st,                        /* i/o: encoder state structure                 */
    const short nBits,                      /* i  : # of Q bits                             */
    const short limit_flag,                 /* i  : restrained(0) or extended(1) Q limits   */
    const short T0,                         /* i  : integer pitch lag                       */
    const short T0_frac,                    /* i  : pitch fraction                          */
    short *T0_min,                    /* i/o: delta search min                        */
    short *T0_max                     /* o  : delta search max                        */
);

short abs_pit_enc(                          /* o  : pitch index                             */
    const short fr_steps,                   /* i  : fractional resolution step              */
    const short limit_flag,                 /* i  : restrained(0) or extended(1) limits     */
    const short T0,                         /* i  : integer pitch lag                       */
    const short T0_frac                     /* i  : pitch fraction                          */
);

short delta_pit_enc(                        /* o  : pitch index                             */
    const short fr_steps,                   /* i  : fractional resolution steps (2 or 4)    */
    const short T0,                         /* i  : integer pitch lag                       */
    const short T0_frac,                    /* i  : pitch fraction                          */
    const short T0_min                      /* i  : delta search min                        */
);

float AVQ_cod(                              /* o:   comfort noise gain factor               */
    const float xri[],                      /* i:   vector to quantize                      */
    int   xriq[],                     /* o:   quantized normalized vector (assuming the bit budget is enough) */
    const short nb_bits,                    /* i:   number of allocated bits                */
    const short Nsv                         /* i:   number of subvectors (lg=Nsv*8)         */
);

void AVQ_encmux(
    Encoder_State *st,                        /* i/o: encoder state structure                         */
    const short extl,                       /* i  : extension layer                                 */
    int   xriq[],                     /* i/o: rounded subvectors [0..8*Nsv-1] followed        */
    /*      by rounded bit allocations [8*Nsv..8*Nsv+Nsv-1] */
    short *nb_bits,                   /* i/o: number of allocated bits                        */
    const short Nsv,                        /* i:   number of subvectors                            */
    short nq[]                        /* o  : AVQ nq index                                    */
);

void re8_cod(
    int x[],                          /* i  : point in RE8 (8-dimensional integer vector)                         */
    int *n,                           /* i  : codebook number (*n is an integer defined in {0,2,3,4,..,n_max})    */
    long *I,                          /* o  : index of c (pointer to unsigned 16-bit word)                        */
    int k[]                           /* o  : index of v (8-dimensional vector of binary indices) = Voronoi index */
);

void pre_exc(
    const short Opt_AMR_WB,                 /* i  : flag indicating AMR-WB IO mode                       */
    const short L_frame,                    /* i  : frame length                                         */
    const float *speech,                    /* i  : input speech                                         */
    const float *p_Aq,                      /* i  : 12k8 Lp coefficient                                  */
    const float *p_A,                       /* i  : unquantized A(q) filter with bandwidth expansion     */
    const short coder_type,                 /* i  : coding type                                          */
    const short i_subfr,                    /* i  : current sub frame indicator                          */
    float *Ap,                        /* o  : weighted LP filter coefficients                      */
    const float *res,                       /* i  : residual signal                                      */
    float *h1,                        /* o  : impulse response of weighted synthesis filter        */
    float *xn,                        /* o  : close-loop Pitch search target vector                */
    float *cn,                        /* o  : target vector in residual domain                     */
    float *mem_syn,                   /* i/o: memory of the synthesis filter                       */
    float *mem_w0,                    /* i/o: weighting filter denominator memory                  */
    const short L_subfr                     /* i  : subframe length                                      */
);

void encod_unvoiced(
    Encoder_State *st,                        /* i/o: state structure                                 */
    LPD_state *mem,                       /* i/o: encoder memories                                */
    const float *speech,                    /* i  : input speech                                    */
    const float Aw[],                       /* i  : weighted A(z) unquantized for subframes         */
    const float *Aq,                        /* i  : LP coefficients                                 */
    const short vad_flag,
    const float *res,                       /* i  : residual signal                                 */
    float *syn,                       /* o  : core synthesis                                  */
    float *tmp_noise,                 /* o  : long-term noise energy                          */
    float *exc,                       /* i/o: current non-enhanced excitation                 */
    float *pitch_buf,                 /* o  : floating pitch values for each subframe         */
    float *voice_factors,             /* o  : voicing factors                                 */
    float *bwe_exc                    /* i/o: excitation for SWB TBE                          */
);

void encod_gen_voic(
    Encoder_State *st,                        /* i/o: state structure                                 */
    LPD_state    *mem,                        /* i/o: encoder memories                                */
    const short L_frame,                    /* i  : length of the frame                             */
    const short sharpFlag,                  /* i  : formant sharpening flag                         */
    const float speech[],                   /* i  : input speech                                    */
    const float Aw[],                       /* i  : weighted A(z) unquantized for subframes         */
    const float Aq[],                       /* i  : LP coefficients                                 */
    const short coder_type,                 /* i  : coding type                                     */
    const float Es_pred,                    /* i  : predicted scaled innov. energy                  */
    const short T_op[],                     /* i  : open loop pitch                                 */
    const float voicing[],                  /* i  : voicing                                         */
    const float *res,                       /* i  : residual signal                                 */
    float *syn,                       /* o  : core synthesis                                  */
    float *exc,                       /* i/o: current non-enhanced excitation                 */
    float *exc2,                      /* i/o: current enhanced excitation                     */
    float *pitch_buf,                 /* o  : floating pitch values for each subframe         */
    float *voice_factors,             /* o  : voicing factors                                 */
    float *bwe_exc,                   /* i/o: excitation for SWB TBE                          */
    short *unbits                     /* i/o: number of unused bits                           */
);

short encod_tran(
    Encoder_State *st,                        /* i/o: state structure                                 */
    LPD_state *mem,                       /* i/o: encoder memories                                */
    const short L_frame,                    /* i  : length of the frame                             */
    const float speech[],                   /* i  : input speech                                    */
    const float Aw[],                       /* i  : weighted A(z) unquantized for subframes         */
    const float Aq[],                       /* i  : LP coefficients                                 */
    const short coder_type,                 /* i  : coding type                                     */
    const float Es_pred,                    /* i  : predicted scaled innov. energy                  */
    const short T_op[],                     /* i  : open loop pitch                                 */
    const float voicing[],                  /* i  : voicing                                         */
    const float *res,                       /* i  : residual signal                                 */
    float *syn,                       /* o  : synthesis                                       */
    float *exc,                       /* i/o: current non-enhanced excitation                 */
    float *exc2,                      /* i/o: current enhanced excitation                     */
    float *pitch_buf,                 /* o  : floating pitch values for each subframe         */
    float *voice_factors,             /* o  : voicing factors                                 */
    float *bwe_exc,                   /* i/o: excitation for SWB TBE                          */
    const short attack_flag,                /* i  : Flag to indicate when attack is deal with TM    */
    short *unbits,                    /* i/o: number of unused bits                           */
    const short sharpFlag                   /* i  : formant sharpening flag                         */
);

void encod_amr_wb(
    Encoder_State *st,                        /* i/o: state structure                                 */
    LPD_state *mem,                       /* i/o: encoder state structure                         */
    const float speech[],                   /* i  : input speech                                    */
    const float Aw[],                       /* i  : weighted A(z) unquantized for subframes         */
    const float Aq[],                       /* i  : 12k8 Lp coefficient                             */
    const short pitch[3],                   /* i  : open-loop pitch values for quantiz.             */
    const float voicing[],                  /* i  : voicing                                         */
    const float *res,                       /* i  : residual signal                                 */
    float *syn,                       /* i/o: core synthesis                                  */
    float *exc,                       /* i/o: current non-enhanced excitation                 */
    float *exc2,                      /* i/o: current enhanced excitation                     */
    float *pitch_buf,                 /* i/o: floating pitch values for each subframe         */
    short *hf_gain,                   /* o  : decoded HF gain                                 */
    const float *speech16k                  /* i  : input speech @16kHz                             */
);

void stat_noise_uv_enc(
    Encoder_State *st,                        /* i/o: state structure                                 */
    const short coder_type,                 /* i  : coding type                                     */
    const float *epsP,                      /* i  : LP prediction errors                            */
    float *lsp_new,                   /* i  : end-frame LSP vector                            */
    float *lsp_mid,                   /* i  : mid-frame LSP vector                            */
    float *Aq,                        /* i  : A(z) quantized for the 4 subframes              */
    float *exc2                       /* i/o: excitation buffer                               */
);

void re8_compute_base_index(
    const int *x,                         /* i  : Elemen of Q2, Q3 or Q4                          */
    const int ka,                         /* i  : Identifier of the absolute leader related to x  */
    long *I                          /* o  : index                                           */
);

void transf_cdbk_enc(
    Encoder_State *st,                        /* i/o: encoder state structure                     */
    const long  core_brate,                 /* i  : core bitrate                                */
    const short extl,                       /* i  : extension layer                             */
    const short coder_type,                 /* i  : coding type                                 */
    const short harm_flag_acelp,            /* i  : harmonic flag for higher rates ACELP        */
    const short i_subfr,                    /* i  : subframe index                              */
    const short tc_subfr,                   /* i  : TC subframe index                           */
    float cn[],                       /* i/o: target vector in residual domain            */
    float exc[],                      /* i/o: pointer to excitation signal frame          */
    const float *p_Aq,                      /* i  : 12k8 Lp coefficient                         */
    const float Ap[],                       /* i  : weighted LP filter coefficients             */
    const float h1[],                       /* i  : weighted filter input response              */
    float xn[],                       /* i/o: target vector                               */
    float xn2[],                      /* i/o: target vector for innovation search         */
    float y1[],                       /* i/o: zero-memory filtered adaptive excitation    */
    const float y2[],                       /* i  : zero-memory filtered innovative excitation  */
    const float Es_pred,                    /* i  : predicited scaled innovation energy         */
    float *gain_pit,                  /* i/o: adaptive excitation gain                    */
    const float gain_code,                  /* i  : innovative excitation gain                  */
    float g_corr[],                   /* o  : ACELP correlation values                    */
    const short clip_gain,                  /* i  : adaptive gain clipping flag                 */
    float *mem_deemp,                 /* i/o: prequantizer deemhasis memory               */
    float *mem_preemp,                /* i/o: prequantizer preemhasis memory              */
    float *gain_preQ,                 /* o  : prequantizer excitation gain                */
    float code_preQ[],                /* o  : prequantizer excitation                     */
    short *unbits                     /* i/o: number of AVQ unused bits                   */
);

short gain_quant(                           /* o:   quantization index                          */
    float *gain,                      /* i/o: quantized gain                              */
    const float min,                        /* i:   value of lower limit                        */
    const float max,                        /* i:   value of upper limit                        */
    const short bits                        /* i:   number of bits to quantize                  */
);

void deemph_lpc(
    float *p_Aq_cuerr,                /* i : LP coefficients current frame                */
    float *p_Aq_old,                  /* i : LP coefficients previous frame               */
    float *LPC_de_curr,               /* o : De-emphasized LP coefficients current frame  */
    float *LPC_de_old                 /* o : De-emphasized LP coefficients previous frame */
    ,short deemph_old
);

void Interpol_delay(
    float *out,                       /* o : pitch interpolation output                   */
    float *last,                      /* i : last frame pitch lag                         */
    float *current,                   /* i : current frame pitch lag                      */
    short SubNum,                     /* i : subframe number                              */
    const float *frac                       /* i : interpolation constant                       */
);

float dequantize_uvg(
    int   iG1,                        /* i: gain 1 index                                  */
    int   *iG2,                       /* i: gain 2 index                                  */
    float *G,                         /* o: quantized gain                                */
    short bwidth
);

void generate_nelp_excitation(
    short *seed,                      /* i/o: random number seed                          */
    float *Gains,                     /* i  : excitation gains                            */
    float *output,                    /* o  : excitation output                           */
    float gain_fac                    /* i  : gain factor                                 */
);

void nelp_encoder(
    Encoder_State *st,                        /* i/o: encoder state                               */
    float *in,                        /* i  : residual signal                             */
    float *exc                        /* o  : NELP quantized excitation signal            */
    ,short reduce_gains
);

void encod_nelp(
    Encoder_State *st,                        /* i/o: state structure                             */
    LPD_state *mem,                       /* i/o: encoder memories                            */
    const float *speech,                    /* i  : input speech                                */
    const float Aw[],                       /* i  : weighted A(z) unquantized for subframes     */
    const float *Aq,                        /* i  : 12k8 Lp coefficient                         */
    float *res,                       /* o  : residual signal                             */
    float *synth,                     /* o  : core synthesis                              */
    float *tmp_noise,                 /* o  : long-term noise energy                      */
    float *exc,                       /* i/o: current non-enhanced excitation             */
    float *exc2,                      /* i/o: current enhanced excitation                 */
    float *pitch_buf,                 /* o  : floating pitch values for each subframe     */
    float *voice_factors,             /* o  : voicing factors                             */
    float *bwe_exc                    /* o  : excitation for SWB TBE                      */
);

void realft(
    float *data,                      /* i/o: data array   .......... */
    short n,                          /* i  : length of data array    */
    short isign                       /* i  : sign +1 or -1           */
);

DTFS_STRUCTURE *DTFS_new(
    void
);

void DTFS_copy(
    DTFS_STRUCTURE *Xout,                      /* o: DTFS       */
    DTFS_STRUCTURE Xinp                        /* i: DTFS       */
);

DTFS_STRUCTURE DTFS_sub(
    DTFS_STRUCTURE X1,                         /* i: X1 DTFS */
    DTFS_STRUCTURE X2                          /* i: X2 DTFS */
    /* o: X1 - X2 */
);

void DTFS_to_fs(
    const float *x,                         /* i: Time domain signal            */
    int N,                          /* i: Length of input vector        */
    DTFS_STRUCTURE *X,                         /* o: DTFS structure with a, b, lag */
    short Fs,
    short FR_flag                     /* i: FR flag                       */
);

void DTFS_fs_inv(
    DTFS_STRUCTURE *X,                         /* i : DTFS             */
    float *x,                         /* o : time domain sig  */
    int   N,                          /* i : Output length    */
    float ph0                         /* i : Input phase      */
);

void DTFS_fast_fs_inv(
    DTFS_STRUCTURE *X1_DTFS,                   /* i : X1 DTFS                  */
    float *out,                       /* o : time domain              */
    int   N                           /* i : Number of output samples */
);

void DTFS_car2pol(
    DTFS_STRUCTURE *X                          /* i/o : DTFS structure a, b, lag      */
    /* input in Cartesion, output in Polar */
);

void DTFS_pol2car(
    DTFS_STRUCTURE *X                          /* i/o : DTFS structure a, b, lag        */
    /* input in Polar, output in Cartesian   */
);

float DTFS_setEngyHarm(
    float f1,                         /* i  : lower band freq of input to control energy   */
    float f2,                         /* i  : upper band freq of input to control energy   */
    float g1,                         /* i  : lower band freq of output to control energy  */
    float g2,                         /* i  : upper band freq of output to control energy  */
    float en2,                        /* i  : Target Energy to set the DTFS to             */
    DTFS_STRUCTURE *X                          /* i/o: DTFS to adjust the energy of                 */
    /* o  : Return Input RMS between f1/f2 b4 scaling    */
);

void DTFS_to_erb(
    DTFS_STRUCTURE X,                 /* i : DTFS input   */
    float *out                        /* o : ERB output   */
);

void DTFS_zeroPadd(
    int   N,                          /* i  : Target lag  */
    DTFS_STRUCTURE *X                          /* i/o: DTFS        */
);

float DTFS_getEngy(
    DTFS_STRUCTURE X                           /* i : DTFS to compute energy of  */
    /* o : Energy                     */
);

float DTFS_getEngy_band(
    DTFS_STRUCTURE X,                          /* i : DTFS to compute energy of     */
    float lband,                      /* i : low end of band of interest   */
    float hband                       /* i : high end of band of interest  */
);

float DTFS_getEngy_band_wb(
    DTFS_STRUCTURE X,                          /* i : DTFS to compute energy of    */
    float lband,                      /* i : low end of band of interest  */
    float hband                       /* i : high end of band of interest */
);

double DTFS_freq_corr(
    DTFS_STRUCTURE X1,                         /* i : X1 DTFS      */
    DTFS_STRUCTURE X2,                         /* i : X2 DTFS      */
    float lband,                      /* i : low cutoff   */
    float hband                       /* i : high cutoff  */
    /* o : Correlation  */
);

float DTFS_setEngy(
    DTFS_STRUCTURE *X_DTFS,                    /* i/o: DTFS structure to set engy */
    float en2                         /* i  : Energy to set to           */
);

void DTFS_adjustLag(
    DTFS_STRUCTURE *X_DTFS,                    /* i/o: DTFS to adjust lag for   */
    int   N                           /* i  : Target lag               */
);

void DTFS_poleFilter(
    DTFS_STRUCTURE *X,                         /* i/o: DTFS to poleFilter inplace  */
    const float *LPC,                       /* i  : LPCs                        */
    int   N                           /* i  : LPCORDER                    */
);

void DTFS_zeroFilter(
    DTFS_STRUCTURE *X,                         /* i/o: DTFS to zeroFilter inplace   */
    const float *LPC,                       /* i  : LPCs                         */
    int   N                           /* i  : LPCORDER                     */
);

float DTFS_alignment_full(
    DTFS_STRUCTURE X1_DTFS,                    /* i : reference DTFS               */
    DTFS_STRUCTURE X2_DTFS,                    /* i : DTFS to shift                */
    int   num_steps                   /* i : resolution                   */
);

float DTFS_alignment_extract(
    DTFS_STRUCTURE refX1_DTFS,                 /* i : X1 the reference DTFS to keep fixed           */
    DTFS_STRUCTURE X2_DTFS,                    /* i : X2 the test DTFS to shift to find best match  */
    float Eshift,                     /* i : Expected shift - coarse value                 */
    const float *LPC2                       /* i : LPC to filter to find correlation in speech   */
    /* o : shift value to shift X2 by                    */
);

float DTFS_alignment_weight(
    DTFS_STRUCTURE refX1_DTFS,                 /* i : X1 the reference DTFS to keep fixed             */
    DTFS_STRUCTURE X2_DTFS,                    /* i : X2 the test DTFS to shift to find best match    */
    float Eshift,                     /* i : Expected shift - coarse value                   */
    const float *LPC1,                      /* i : LPC to filter to find correlation in spch       */
    const float *LPC2                       /* i : LPC to filter to find correlation in spch       */
    /* o : shift value to shift X2 by                      */
);

float DTFS_alignment_fine_new(
    DTFS_STRUCTURE X1_DTFS,                    /* i : X1 the reference DTFS to keep fixed            */
    DTFS_STRUCTURE X2_DTFS,                    /* i : X2 the test DTFS to shift to find best match   */
    float Eshift                      /* i : Expected shift - coarse value                  */
    /* o : shift value to shift X2 by                     */
);

void DTFS_phaseShift(
    DTFS_STRUCTURE *X,                         /* i : DTFS to shift                    */
    float ph                          /* i : phase to shift                   */
);

void DTFS_erb_inv(
    float *in,                        /* i : ERB inpt                        */
    int   *slot,                      /* i : ERB slots filled based on lag   */
    float *mfreq,                     /* i : erb frequence edges             */
    DTFS_STRUCTURE * X,                        /* o : DTFS after erb-inv              */
    short num_erb                     /* i : Number of ERB bands             */
);

short DTFS_quant_cw(
    DTFS_STRUCTURE *X,                         /* i/o: DTFS unquant inp, quant out        */
    int   pl,                         /* i  : Previous lag                       */
    const float *curr_lpc,                  /* i  : LPC                                */
    int   *POWER_IDX,                 /* o  : Power index                        */
    int   *AMP_IDX,                   /* o  : Amplitude index                    */
    float *lastLgainE,                /* i/o: last frame lowband gain            */
    float *lastHgainE,                /* i/o: last frame highband gain           */
    float *lasterbE                   /* i/o: last frame ERB vector              */
    /* Flag - amp quant performance pass/fail  */
);

void DTFS_dequant_cw(
    int   pl,                         /* i  : Previous lag                */
    int   POWER_IDX,                  /* i  : POWER index                 */
    const int   *AMP_IDX,                   /* i  : Amp Shape index             */
    float *lastLgainD,                /* i/o: low band last gain          */
    float *lastHgainD,                /* i/o: high band last gain         */
    float *lasterbD,                  /* i/o: last frame ERB vector       */
    DTFS_STRUCTURE *X,                         /* o  : DTFS structure dequantized  */
    short num_erb                     /* i  : Number of ERB bands         */
);

void DTFS_transform(
    DTFS_STRUCTURE X,                          /* i : Starting DTFS to use in WI    */
    DTFS_STRUCTURE X2,                         /* i : Ending DTFS to use in WI      */
    const float *phase,                     /* i : Phase contour                 */
    float *out,                       /* o : Output time domain waveform   */
    int   N,                          /* i : Number of samples to generate */
    short FR_flag                     /* i : Flag to indicate called in FR context */
);

float DTFS_getSpEngyFromResAmp(
    DTFS_STRUCTURE X,                          /* i : DTFS                             */
    float lband,                      /* i : Low band end to get energy from  */
    float hband,                      /* i : High band end to get energy from */
    const float *curr_lsp                   /* i : LPCs                             */
    /* o : speech energy                    */
);

void DTFS_peaktoaverage(
    DTFS_STRUCTURE X,                          /* i : DTFS                  */
    float *pos,                       /* o : positive peak to ave  */
    float *neg                        /* o : negative peak to ave  */
);

short ppp_extract_pitch_period(
    const float *in,                        /* i : input residual       */
    float *out,                       /* o : output residual      */
    int   l,                          /* i : lag                  */
    short *out_of_bound               /* o : out of bound flag    */
);

short ppp_quarter_encoder(
    Encoder_State *st,                        /* i/o: encoder state structure                  */
    DTFS_STRUCTURE *CURRCW_Q,                  /* o  : Quantized (amp/phase) DTFS               */
    DTFS_STRUCTURE *TARGETCW,                  /* o  : DTFS with quant phase but unquant Amp    */
    int   prevCW_lag,                 /* i  : previous lag                             */
    DTFS_STRUCTURE vCURRCW_NQ,                 /* i  : Unquantized DTFS                         */
    const float *curr_lpc,                  /* i  : LPCS                                     */
    float *lastLgainE,                /* i/o: last low band gain                       */
    float *lastHgainE,                /* i/o: last high band gain                      */
    float *lasterbE,                  /* i/o: last ERB vector                          */
    DTFS_STRUCTURE PREV_CW_E                   /* i  : past DTFS                                */
);

void WIsyn(
    DTFS_STRUCTURE PREVCW,                     /* i  : Prev frame DTFS                          */
    DTFS_STRUCTURE *CURR_CW_DTFS,              /* i/o: Curr frame DTFS                          */
    const float *curr_lpc,                  /* i  : LPC                                      */
    float *ph_offset,                 /* i/o: Phase offset to line up at end of frame  */
    float *out,                       /* o  : Waveform Interpolated time domain signal */
    int   N,                          /* i  : Number of samples of output to generate  */
    int FR_flag                       /* i  : called for post-smoothing in FR          */
);

void set_ppp_mode(
    Encoder_State *st,                      /* i/o: encoder state structure                 */
    short *coder_type,              /* i/o: coder type                              */
    const short noisy_speech_HO,          /* i  : SC-VBR noisy speech HO flag             */
    const short clean_speech_HO,          /* i  : SC-VBR clean speech HO flag             */
    const short NB_speech_HO,             /* i  : SC-VBR NB speech HO flag                */
    const short localVAD,
    const short localVAD_he,              /* i  : HE-SAD flag without hangover            */
    short *vad_flag,
    short T_op[],                   /* i  :open loop pitch lag                      */
    short sp_aud_decision1          /* i : Speech Audio Decision                    */
);

void lsf_syn_mem_backup(
    Encoder_State *st,                        /* i: state structure                                       */
    LPD_state* LPDmem,                /* i: LPD state memory structure                            */
    float *btilt_code,                /* i: tilt code                                             */
    float *bgc_threshold,             /* i:                                                       */
    float *clip_var_bck,              /* o:                                                       */
    short *next_force_sf_bck,         /* o:                                                       */
    float *lsp_new,                   /* i: LSP vector to quantize                                */
    float *lsp_mid,                   /* i: mid-frame LSP vector                                  */
    float *clip_var,                  /* o: pitch clipping state var                              */
    float *mem_AR,                    /* o: quantizer memory for AR model                         */
    float *mem_MA,                    /* o: quantizer memory for AR model                         */
    float *lsp_new_bck,               /* o: LSP vector to quantize- backup                        */
    float *lsp_mid_bck,               /* o: mid-frame LSP vector - backup                         */
    short *mCb1,                      /* o: counter for stationary frame after a transition frame */
    float *Bin_E,                     /* o: FFT Bin energy 128 *2 sets                            */
    float *Bin_E_old,                 /* o: FFT Bin energy 128 sets                               */
    float *mem_syn_bck,               /* o: synthesis filter memory                               */
    float *mem_w0_bck,                /* o: memory of the weighting filter                        */
    float *streaklimit,
    short *pstreaklen
);

void lsf_syn_mem_restore(
    Encoder_State *st,                        /* o: state structure                                        */
    LPD_state* LPDmem,                /* o: LPD_state vewctor                                      */
    float btilt_code,                 /* i:                                                        */
    float gc_threshold,               /* i:                                                        */
    float *clip_var_bck,              /* i:                                                        */
    short next_force_sf_bck,          /* i:                                                        */
    float *lsp_new,                   /* o: LSP vector to quantize                                 */
    float *lsp_mid,                   /* o: mid-frame LSP vector                                   */
    float clip_var,                   /* i: pitch clipping state var                               */
    float *mem_AR,                    /* i: quantizer memory for AR model                          */
    float *mem_MA,                    /* i: quantizer memory for AR model                          */
    float *lsp_new_bck,               /* i: LSP vector to quantize- backup                         */
    float *lsp_mid_bck,               /* i: mid-frame LSP vector - backup                          */
    short mCb1,                       /* i: counter for stationary frame after a transition frame  */
    float *Bin_E,                     /* i: FFT Bin energy 128 *2 sets                             */
    float *Bin_E_old,                 /* i: FFT Bin energy 128 sets                                */
    float *mem_syn_bck,               /* i: synthesis filter memory                                */
    float mem_w0_bck,                 /* i: memory of the weighting filter                         */
    float streaklimit,
    short pstreaklen
);
void ppp_voiced_encoder(
    Encoder_State *st,                        /* i/o: state structure                             */
    float *in,                        /* i  : residual signal                             */
    float *out,                       /* o  : Quantized residual signal                   */
    short delay,                      /* i  : open loop pitch                             */
    float *lpc1,                      /* i  : prev frame de-emphasized LPC                */
    float *lpc2,                      /* i  : current frame de-emphasized LPC             */
    float *exc,                       /* i  : previous frame quantized excitation         */
    float *pitch,                     /* o  : floating pitch values for each subframe     */
    float vadsnr                      /* i  : current frame SNR                           */
);

void encod_ppp(
    Encoder_State *st,                        /* i/o: state structure                                 */
    LPD_state *mem,                       /* i/o: encoder memories                                */
    const float speech[],                   /* i  : input speech                                    */
    const float Aw[],                       /* i  : weighted A(z) unquantized for subframes         */
    const float Aq[],                       /* i  : 12k8 Lp coefficient                             */
    short *coder_type,                /* i/o: coding type                                     */
    const short sharpFlag,                  /* i  : formant sharpening flag                         */
    const short T_op[],                     /* i  : open loop pitch                                 */
    const float voicing[],                  /* i  : voicing                                         */
    float *res,                       /* i/o: residual signal                                 */
    float *synth,                     /* i/o: core synthesis                                  */
    float *exc,                       /* i/o: current non-enhanced excitation                 */
    float *exc2,                      /* i/o: current enhanced excitation                     */
    float *pitch_buf,                 /* i/o: floating pitch values for each subframe         */
    float *voice_factors,             /* o  : voicing factors                                 */
    float *bwe_exc                    /* o  : excitation for SWB TBE                          */
);

void reset_rf_indices(
    Encoder_State *st                       /* i: state structure - contains partial RF indices     */
);

void signalling_enc_rf(
    Encoder_State *st                         /* i/o: encoder state structure             */
);

void acelp_core_dec(
    Decoder_State *st,                        /* i/o: Decoder state structure             */
    float synth_out[],                /* o  : synthesis                           */
    float bwe_exc_extended[],         /* i/o: bandwidth extended excitation       */
    float *voice_factors,             /* o  : voicing factors                     */
    float old_syn_12k8_16k[],         /* o  : intermediate ACELP synthesis at 12.8kHz or 16kHz to be used by SWB BWE */
    short coder_type,                 /* i  : coder type                          */
    short sharpFlag,                  /* i  : formant sharpening flag             */
    float pitch_buf[NB_SUBFR16k],     /* o  : floating pitch for each subframe    */
    short *unbits,                    /* o  : number of unused bits               */
    short *sid_bw                     /* o  : 0-NB/WB, 1-SWB SID                  */
);

void bass_psfilter_init(
    float old_syn[],                  /* o  : Old synthesis buffer 1              */
    float *mem_deemph_err,            /* o  : Error deemphasis memory             */
    float *lp_ener                    /* o  : long_term error signal energy       */
);

void bass_psfilter(
    const short Opt_AMR_WB,                 /* i  : AMR-WB IO flag                      */
    const float synth_in[],                 /* i  : synthesis (at 16kHz)                */
    const short L_frame,                    /* i  : length of the last frame            */
    const float pitch_buf[],                /* i  : pitch for every subfr [0,1,2,3]     */
    float old_syn[],                  /* i/o: NBPSF_PIT_MAX                       */
    float *mem_deemph_err,            /* o  : Error deemphasis memory             */
    float *lp_ener,                   /* o  : long_term error signal energy       */
    const short bpf_off,                    /* i  : do not use BPF when set to 1        */
    float v_stab,                     /* i  : stability factor                    */
    float *v_stab_smooth,             /* i  : smoothed stability factor           */
    float *mem_mean_pit,              /* i/o: average pitch memory                */
    short *Track_on_hist,             /* i/o: History of half frame usage         */
    short *vibrato_hist,              /* i/o: History of frames declared as vibrato*/
    float *psf_att,                   /* i/o: Post filter attenuation factor      */
    const short coder_type,                 /* i  : coder_type                          */
    float bpf_noise_buf[]             /* o  : BPF error signal (at int_fs)        */
);

void CNG_reset_dec(
    Decoder_State *st,                        /* i/o: decoder state structure            */
    float *pitch_buf,                 /* o  : floating pitch for each subframe   */
    float *voice_factors              /* o  : voicing factors                    */
);

void updt_dec(
    Decoder_State *st,                      /* i/o: state structure                         */
    const short L_frame,                  /* i  : length of the frame                     */
    const short coder_type,               /* i  : coding type                             */
    const float *old_exc,                 /* i  : buffer of excitation                    */
    const float *pitch_buf,               /* i  : floating pitch values for each subframe */
    const float Es_pred,                  /* i  : predicited scaled innovation energy     */
    const float *Aq,                      /* i  : A(z) quantized for all subframes        */
    const float *lsf_new,                 /* i  : current frame LSF vector                */
    const float *lsp_new,                 /* i  : current frame LSP vector                */
    const float voice_factors[],          /* i  : voicing factors                         */
    const float *old_bwe_exc,             /* i  : buffer of excitation                    */
    const float *gain_buf
);

void updt_IO_switch_dec(
    const short output_frame,               /* i  : output frame length                     */
    Decoder_State *st                         /* i/o: state structure                         */
);

void updt_dec_common(
    Decoder_State *st,                        /* i/o: decoder state structure                 */
    const short hq_core_type,               /* i  : HQ core type                            */
    const float *synth                      /* i  : decoded synthesis                       */
);

void CNG_dec(
    Decoder_State *st,                        /* i/o: State structure                         */
    const short L_frame,                    /* i  : length of the frame                     */
    float Aq[],                       /* o  : LP coefficients                         */
    const long  core_brate,                 /* i  : core bitrate                            */
    float *lsp_new,                   /* i/o: current frame LSPs                      */
    float *lsf_new,                   /* i/o: current frame LSFs                      */
    short *allow_cn_step,             /* o  : allow cn step                           */
    short *sid_bw,                    /* o  : 0-NB/WB, 1-SWB SID                      */
    float *q_env
);

void swb_CNG_dec(
    Decoder_State *st,                        /* i/o: State structure                         */
    const float *synth,                     /* i  : ACELP core synthesis at 32kHz           */
    float *shb_synth,                 /* o  : high-band CNG synthesis                 */
    const short sid_bw                      /* i  : 0-NB/WB, 1-SWB SID                      */
);

void lsf_dec(
    Decoder_State *st,                        /* i/o: State structure                             */
    const short tc_subfr,                   /* i  : TC subframe index                           */
    const short L_frame,                    /* i  : length of the frame                         */
    const short coder_type,                 /* i  : coding type                                 */
    const short bwidth,                     /* i  : input signal bandwidth                      */
    float *Aq,                        /* o  : quantized A(z) for 4 subframes              */
    short *LSF_Q_prediction,          /* o  : LSF prediction mode                         */
    float *lsf_new,                   /* o  : de-quantized LSF vector                     */
    float *lsp_new,                   /* o  : de-quantized LSP vector                     */
    float *lsp_mid                    /* o  : de-quantized mid-frame LSP vector           */
);

void isf_dec_amr_wb(
    Decoder_State *st,                      /* i/o: State structure                             */
    float *Aq,                        /* o  : quantized A(z) for 4 subframes              */
    float *isf_new,                   /* o  : de-quantized ISF vector                     */
    float *isp_new                    /* o  : de-quantized ISP vector                     */
);

void Es_pred_dec(
    float *Es_pred,                   /* o  : predicited scaled innovation energy         */
    const int   enr_idx,                    /* i  : indice                                      */
    const short nb_bits,                    /* i  : number of bits                              */
    const short no_ltp                      /* i  : no LTP flag                                 */
);

void gaus_dec(
    Decoder_State *st,                        /* i/o: decoder state structure                     */
    const long  core_brate,                 /* i  : core bitrate                                */
    const short i_subfr,                    /* i  : subframe index                              */
    float *code,                      /* o  : gaussian excitation                         */
    float *norm_gain_code,            /* o  : gain of the normalized gaussian excitation  */
    float *lp_gainp,                  /* i/o: lp filtered pitch gain(FER)                 */
    float *lp_gainc,                  /* i/o: lp filtered code gain (FER)                 */
    float *gain_inov,                 /* o  : unscaled innovation gain                    */
    float *tilt_code,                 /* o  : synthesis excitation spectrum tilt          */
    float *voice_fac,                 /* o  : estimated voicing factor                    */
    float *gain_pit,                  /* o  : reset pitch gain                            */
    float *pt_pitch,                  /* o  : reset floating pitch buffer                 */
    float *exc,                       /* o  : excitation signal frame                     */
    float *gain_code,                 /* o  : gain of the gaussian excitation             */
    float *exc2                       /* o  : scaled excitation signal frame              */
);

void gain_dec_amr_wb(
    Decoder_State *st,                        /* i/o: decoder state structure               */
    const long  core_brate,                 /* i  : core bitrate                          */
    float *gain_pit,                  /* o  : Quantized pitch gain                  */
    float *gain_code,                 /* o  : Quantized codeebook gain              */
    float *past_qua_en,               /* i/o: gain quantization memory (4 words)    */
    float *gain_inov,                 /* o  : unscaled innovation gain              */
    const float *code,                      /* i  : algebraic code excitation             */
    float *norm_gain_code             /* o  : norm. gain of the codebook excitation */
);

void gain_dec_lbr(
    Decoder_State *st,                        /* i/o: decoder state structure                          */
    const long  core_brate,                 /* i  : core bitrate                                     */
    const short coder_type,                 /* i  : coding type                                      */
    const short i_subfr,                    /* i  : subframe index                                   */
    const float *code,                      /* i  : algebraic excitation                             */
    float *gain_pit,                  /* o  : quantized pitch gain                             */
    float *gain_code,                 /* o  : quantized codebook gain                          */
    float *gain_inov,                 /* o  : gain of the innovation (used for normalization)  */
    float *norm_gain_code,            /* o  : norm. gain of the codebook excitation            */
    float gains_mem[]                 /* i/o: pitch gain and code gain from previous subframes */
);

void gain_dec_mless(
    Decoder_State *st,                        /* i/o: decoder state structure               */
    const long  core_brate,                 /* i  : core bitrate                          */
    const short L_frame,                    /* i  : length of the frame                   */
    const short coder_type,                 /* i  : coding type                           */
    const short i_subfr,                    /* i  : subframe number                       */
    const short tc_subfr,                   /* i  : TC subframe index                     */
    const float *code,                      /* i  : algebraic code excitation             */
    const float Es_pred,                    /* i  : predicted scaled innov. energy        */
    float *gain_pit,                  /* o  : Quantized pitch gain                  */
    float *gain_code,                 /* o  : Quantized codeebook gain              */
    float *gain_inov,                 /* o  : unscaled innovation gain              */
    float *norm_gain_code             /* o  : norm. gain of the codebook excitation */
);

void gain_dec_SQ(
    Decoder_State *st,                        /* i/o: decoder state structure               */
    const long  core_brate,                 /* i  : core bitrate                          */
    const short coder_type,                 /* i  : coding type                           */
    const short i_subfr,                    /* i  : subframe number                       */
    const short tc_subfr,                   /* i  : TC subframe index                     */
    const float *code,                      /* i  : algebraic code excitation             */
    const float Es_pred,                    /* i  : predicted scaled innov. energy        */
    float *gain_pit,                  /* o  : Quantized pitch gain                  */
    float *gain_code,                 /* o  : Quantized codeebook gain              */
    float *gain_inov,                 /* o  : unscaled innovation gain              */
    float *norm_gain_code             /* o  : norm. gain of the codebook excitation */
);

float gain_dec_gaus(                        /* o  : quantized codebook gain                 */
    const short index,                      /* i  : quantization index                      */
    const short bits,                       /* i  : number of bits to quantize              */
    const float lowBound,                   /* i  : lower bound of quantizer (dB)           */
    const float topBound,                   /* i  : upper bound of quantizer (dB)           */
    const float gain_inov,                  /* i  : unscaled innovation gain                */
    float *norm_gain_code             /* o  : gain of normalized gaus. excit.         */
);

float pit_decode(                           /* o  : floating pitch value                    */
    Decoder_State *st,                        /* i/o: decoder state structure                 */
    const long  core_brate,                 /* i  : core bitrate                            */
    const short Opt_AMR_WB,                 /* i  : flag indicating AMR-WB IO mode          */
    const short L_frame,                    /* i  : length of the frame                     */
    short i_subfr,                    /* i  : subframe index                          */
    const short coder_type,                 /* i  : coding type                             */
    short *limit_flag,                /* i/o: restrained(0) or extended(1) Q limits   */
    short *T0,                        /* o  : close loop integer pitch                */
    short *T0_frac,                   /* o  : close loop fractional part of the pitch */
    short *T0_min,                    /* i/o: delta search min for sf 2 & 4           */
    short *T0_max,                    /* i/o: delta search max for sf 2 & 4           */
    const short L_subfr                     /* i  : subframe length                         */
);

void abs_pit_dec(
    const short fr_steps,                   /* i  : fractional resolution steps (0, 2, 4)   */
    short pitch_index,                      /* i  : pitch index                             */
    const short limit_flag,                 /* i  : restrained(0) or extended(1) Q limits   */
    short *T0,                              /* o  : integer pitch lag                       */
    short *T0_frac                          /* o  : pitch fraction                          */
);

void delta_pit_dec(
    const short fr_steps,                   /* i  : fractional resolution steps (0, 2, 4)   */
    const short pitch_index,                /* i  : pitch index                             */
    short *T0,                        /* o  : integer pitch lag                       */
    short *T0_frac,                   /* o  : pitch fraction                          */
    const short T0_min                      /* i  : delta search min                        */
);

void pit_Q_dec(
    const short Opt_AMR_WB,                 /* i  : flag indicating AMR-WB IO mode          */
    const short pitch_index,                /* i  : pitch index                             */
    const short nBits,                      /* i  : # of Q bits                             */
    const short delta,                      /* i  : Half the CL searched interval           */
    const short pit_flag,                   /* i  : absolute(0) or delta(1) pitch Q         */
    const short limit_flag,                 /* i  : restrained(0) or extended(1) Q limits   */
    short *T0,                        /* o  : integer pitch lag                       */
    short *T0_frac,                   /* o  : pitch fraction                          */
    short *T0_min,                    /* i/o: delta search min                        */
    short *T0_max                     /* i/o: delta search max                        */
    ,short *BER_detect                 /* o  : BER detect flag                         */
);

void pit16k_Q_dec(
    const short pitch_index,                /* i  : pitch index                             */
    const short nBits,                      /* i  : # of Q bits                             */
    const short limit_flag,                 /* i  : restrained(0) or extended(1) Q limits   */
    short *T0,                        /* o  : integer pitch lag                       */
    short *T0_frac,                   /* o  : pitch fraction                          */
    short *T0_min,                    /* i/o: delta search min                        */
    short *T0_max                     /* i/o: delta search max                        */
    ,short *BER_detect   /* o  : BER detect flag                         */
);

void lp_filt_exc_dec(
    Decoder_State *st,                        /* i/o: decoder state structure                         */
    const short codec_mode,                 /* i  : codec mode                                      */
    const long  core_brate,                 /* i  : core bitrate                                    */
    const short Opt_AMR_WB,                 /* i  : flag indicating AMR-WB IO mode                  */
    const short coder_type,                 /* i  : coding type                                     */
    const short i_subfr,                    /* i  : subframe index                                  */
    const short L_subfr,                    /* i  : subframe size                                   */
    const short L_Frame,                    /* i  : frame size                                      */
    short lp_flag,                    /* i  : operation mode signalling                       */
    float *exc                        /* i/o: pointer to the excitation signal frame          */
);

void inov_decode(
    Decoder_State *st,                        /* i/o: decoder state structure                         */
    const long  core_brate,                 /* i  : core bitrate                                    */
    const short Opt_AMR_WB,                 /* i  : flag indicating AMR-WB IO mode                  */
    const short L_frame,                    /* i  : length of the frame                             */
    const short coder_type,                 /* i  : coding type                                     */
    const short sharpFlag,                  /* i  : formant sharpening flag                         */
    const short i_subfr,                    /* i  : subframe index                                  */
    const short tc_subfr,                   /* i  : TC subframe index                               */
    const float *p_Aq,                      /* i  : LP filter coefficients                          */
    const float tilt_code,                  /* i  : tilt of of the excitation of previous subframe  */
    const float pt_pitch,                   /* i  : pointer to current subframe fractional pitch    */
    float *code                       /* o  : algebraic excitation                            */
);

void dec_acelp_1t64(
    Decoder_State *st,                        /* i/o: decoder state structure                         */
    float code[]                      /* o  : algebraic (fixed) codebook excitation           */
);

void dec_acelp_2t32(
    Decoder_State *st,                        /* i/o: decoder state structure                         */
    float code[]                      /* o  : algebraic (fixed) codebook excitation           */
);

void dec_acelp_4t64(
    Decoder_State *st,                        /* i/o: decoder state structure                         */
    short nbbits,                     /* i  : number of bits per codebook                     */
    float code[],                     /* o  : algebraic (fixed) codebook excitation           */
    const short Opt_AMR_WB             /* i  : flag indicating AMR-WB IO mode                  */
);

unsigned int syn_output(
    float *synth,                     /* i/o: float synthesis signal                          */
    const short output_frame,               /* i  : output frame length                             */
    short *synth_out                  /* o  : integer 16 bits synthesis signal                */
);

void FEC_exc_estim(
    Decoder_State *st,                        /* i/o: Decoder static memory                           */
    const short L_frame,                    /* i  : length of the frame                             */
    float *old_exc,                   /* i/o: excitation buffer                               */
    float *exc2,                      /* o  : excitation buffer (for synthesis)               */
    float *exc_dct_in,                /* o  : GSC excitation in DCT domain                    */
    float *pitch_buf,                 /* o  : Floating pitch   for each subframe              */
    float *tmp_tc,                    /* o  : FEC pitch                                       */
    float *voice_factors,             /* o  : voicing factors                                 */
    float *bwe_exc,                   /* i/o: excitation for SWB TBE                          */
    float *lsf_new,                   /* i  : ISFs at the end of the frame                    */
    float *tmp_noise                  /* o  : long-term noise energy                          */
);

void FEC_lsf2lsp_interp(
    Decoder_State *st,                        /* i/o: Decoder static memory                        */
    const short L_frame,                    /* i  : length of the frame                          */
    float *Aq,                        /* o  : calculated A(z) for 4 subframes              */
    float *lsf,                       /* o  : estimated LSF vector                         */
    float *lsp                        /* o  : estimated LSP vector                         */
);

void FEC_lsf_estim_enc(
    Encoder_State *st,                        /* i  : Encoder static memory                           */
    const short L_frame,                    /* i  : length of the frame                             */
    float *lsf                        /* o  : estimated LSF vector                            */
);

float frame_energy(
    const short L_frame,                    /* i  : length of the frame                             */
    const float *pitch,                     /* i  : pitch values for each subframe                  */
    const float *speech,                    /* i  : pointer to speech signal for E computation      */
    const float lp_speech,                  /* i  : long term active speech energy average          */
    float *frame_ener                 /* o  : pitch-synchronous energy at frame end           */
);

void FEC_SinOnset(
    float *exc,                       /* i/o : exc vector to modify                                           */
    short puls_pos,                   /* i   : Last pulse position desired                                    */
    short T0,                         /* i   : decoded first frame pitch                                      */
    float enr_q,                      /* i   : energy provide by the encoder                                  */
    float *Aq,                        /* i   : Lsp coefficient                                                */
    const short L_frame                     /* i   : Frame lenght                                                   */
);

short FEC_enhACB(
    const short L_frame,                    /* i   : Frame lenght                                                              */
    const short last_L_frame,               /* i   : frame length of last frame                                                */
    float *exc_io,                    /* i/o : Adaptive codebook memory                                                  */
    const short new_pit,                    /* i   : decoded first frame pitch                                                 */
    const short puls_pos,                   /* i   : decoder position of the last glottal pulses decoded in the previous frame */
    const float bfi_pitch                   /* i   : Pitch used for concealment                                                */
);

void FEC_scale_syn(
    const short L_frame,                    /* i  : length of the frame                     */
    short  clas,                      /* i/o: frame classification                    */
    const short last_good,                  /* i  : last good frame classification          */
    float *synth,                     /* i/o: synthesized speech at Fs = 12k8 Hz      */
    const float *pitch,                     /* i  : pitch values for each subframe          */
    float enr_old,                    /* i  : energy at the end of prvious frame      */
    float enr_q,                      /* i  : transmitted energy for current frame    */
    const short coder_type,                 /* i  : coding type                             */
    const short LSF_Q_prediction,           /* i  : LSF prediction mode                     */
    short *scaling_flag,              /* i/o: flag to indicate energy control of syn  */
    float *lp_ener_FEC_av,            /* i/o: averaged voiced signal energy           */
    float *lp_ener_FEC_max,           /* i/o: averaged voiced signal energy           */
    const short bfi,                        /* i:   current  frame BFI                      */
    const long total_brate,                 /* i:   total bitrate                           */
    const short prev_bfi,                   /* i  : previous frame BFI                      */
    const long  last_core_brate,            /* i  : previous frame core bitrate             */
    float *exc,                       /* i/o: excitation signal without enhancement   */
    float *exc2,                      /* i/o: excitation signal with enhancement      */
    const float Aq[],                       /* i  : LP filter coefs                         */
    float *old_enr_LP,                /* i/o: LP filter E of last good voiced frame   */
    const float *mem_tmp,                   /* i  : temp. initial synthesis filter states   */
    float *mem_syn                    /* o  : initial synthesis filter states         */
    , int   avoid_lpc_burst_on_recovery /* i  : if true the excitation energy is limited if LP has big gain */
    , short force_scaling               /* i: force scaling                             */
);

void FEC_pitch_estim(
    const short Opt_AMR_WB,                 /* i  : flag indicating AMR-WB IO mode          */
    const short last_core,                  /* i  : last core                               */
    const short L_frame,                    /* i  : length of the frame                     */
    const short clas,                       /* i  : current frame classification            */
    const short last_good,                  /* i  : last good clas information              */
    const float pitch_buf[],                /* i  : Floating pitch   for each subframe      */
    const float old_pitch_buf[],            /* i  : buffer of old subframe pitch values     */
    float *bfi_pitch,                 /* i/o: update of the estimated pitch for FEC   */
    short *bfi_pitch_frame,           /* o  : frame length when pitch was updated     */
    short *upd_cnt,                   /* i/o: update counter                          */
    const short coder_type
);

void FEC_encode(
    Encoder_State *st,                        /* i/o: encoder state structure                         */
    const float *synth,                     /* i  : pointer to synthesized speech for E computation */
    const short coder_type,                 /* i  : type of coder                                   */
    short clas,                       /* i  : signal clas for current frame                   */
    const float *fpit,                      /* i  : close loop fractional pitch buffer              */
    const float *res,                       /* i  : LP residual signal frame                        */
    short *last_pulse_pos,            /* i/o: Position of the last pulse                      */
    const short L_frame,                    /* i  : Frame length                                    */
    const long  total_brate,                /* i  : total codec bitrate                             */
    const long  core_brate                  /* i  : core codec bitrate                              */
);

short FEC_pos_dec(
    Decoder_State *st,                        /* i/o: decoder state structure                 */
    const short coder_type,                 /* i  : coder type                              */
    const short last_good,                  /* i  : Last good classfication                 */
    short *last_pulse_pos,            /* o  : Last glotal pulse position in the lost ACB */
    short *clas,                      /* o  : Decoded classification                  */
    float *enr_q,                     /* o  : Decoded energy                          */
    const long  core_brate                  /* i  : Decoded bitrate                         */
);

void improv_amr_wb_gs(
    const short clas,                     /* i  : bitrate allocated to the core           */
    const short coder_type,               /* i  : coder_type                              */
    const long  core_brate,               /* i  : bitrate allocated to the core           */
    short *seed_tcx,                /* i/o: Seed used for noise generation          */
    float *old_Aq,                  /* i/o: old LPC filter coefficient              */
    float *mem_syn2,                /* i/o: synthesis memory                        */
    const float lt_voice_fac,             /* i/o: long term voice factor                  */
    const short locattack,                /* i  : Flag for a detected attack              */
    float *Aq,                      /* i/o: Decoded LP filter coefficient           */
    float *exc2,                    /* i/o: Decoded complete excitation             */
    float *mem_tmp,                 /* i/o: synthesis temporary memory              */
    float *syn,                     /* i/o: Decoded synthesis to be updated         */
    const float *pitch_buf,               /* i  : Decoded pitch buffer                    */
    const float Last_ener                 /* i  : Last energy                             */
    ,const short rate_switching_reset      /* i  : rate switching reset flag               */
    ,const short last_coder_type                    /* i  : Last coder_type */
);

short tc_classif(
    Decoder_State *st,                        /* i/o: decoder state structure                 */
    const short L_frame                     /* i  : length of the frame                     */
);

void transition_dec(
    Decoder_State *st,                        /* i/o: decoder state structure                 */
    const long  core_brate,                 /* i  : core bitrate                            */
    const short Opt_AMR_WB,                 /* i  : flag indicating AMR-WB IO mode          */
    const short L_frame,                    /* i  : length of the frame                     */
    const short i_subfr,                    /* i  : subframe index                          */
    const short coder_type,                 /* i  : coding type                             */
    const short tc_subfr,                   /* i  : TC subframe index                       */
    short *Jopt_flag,                 /* i  : joint optimization flag                 */
    float *exc,                       /* i/o: current frame excitation signal         */
    short *T0,                        /* o  : close loop integer pitch                */
    short *T0_frac,                   /* o  : close loop fractional part of the pitch */
    short *T0_min,                    /* i/o: delta search min for sf 2 & 4           */
    short *T0_max,                    /* i/o: delta search max for sf 2 & 4           */
    float **pt_pitch,                 /* o  : floating pitch values                   */
    short *position,                  /* i/o: first glottal impulse position in frame */
    float *bwe_exc                    /* i/o: excitation for SWB TBE                  */
);

void gain_dec_tc(
    Decoder_State *st,                        /* i/o: decoder state structure              */
    const long  core_brate,                 /* i  : core bitrate                         */
    const short L_frame,                    /* i  : length of the frame                  */
    const short i_subfr,                    /* i  : subframe number                      */
    const short tc_subfr,                   /* i  : TC subframe index                    */
    const float Es_pred,                    /* i  : predicted scaled innov. energy       */
    const float *code,                      /* i  : algebraic code excitation            */
    float *gain_pit,                  /* o  : pitch gain                           */
    float *gain_code,                 /* o  : Quantized codeebook gain             */
    float *gain_inov,                 /* o  : unscaled innovation gain             */
    float *norm_gain_code             /* o  : norm. gain of the codebook excit.    */
);

void stat_noise_uv_dec(
    Decoder_State *st,                        /* i/o: decoder static memory                */
    const short coder_type,                 /* i  : coding type                          */
    float *lsp_new,                   /* i  : end-frame LSP vector                 */
    float *lsp_mid,                   /* i  : mid-frame LSP vector                 */
    float *Aq,                        /* o  : A(z) quantized for the 4 subframes   */
    float *exc2                       /* i/o: excitation buffer                    */
);

void decod_nelp(
    Decoder_State *st,                        /* i/o: decoder static memory                   */
    const short coder_type,                 /* i  : coding type                             */
    float *tmp_noise,                 /* o  : long term temporary noise energy        */
    float *pitch_buf,                 /* o  : floating pitch values for each subframe */
    float *exc,                       /* o  : adapt. excitation exc                   */
    float *exc2,                      /* o  : adapt. excitation/total exc             */
    float *voice_factors,             /* o  : voicing factors                         */
    float *bwe_exc,                   /* o  : excitation for SWB TBE                  */
    const short bfi,                        /* i  : bad frame indicator                     */
    float *gain_buf
);

void nelp_decoder(
    Decoder_State *st,                        /* i/o: decoder static memory       */
    float *exc_nelp,                  /* o  : adapt. excitation/total exc */
    float *exc,                       /* o  : adapt. excitation exc       */
    short bfi,                        /* i  : frame error rate            */
    const short coder_type,                 /* i : coding type                  */
    float *gain_buf
);

void decod_ppp(
    Decoder_State *st,                        /* i/o: state structure                           */
    const float Aq[],                       /* i  : 12k8 Lp coefficient                       */
    float *pitch_buf,                 /* i/o: floating pitch values for each subframe   */
    float *exc,                       /* i/o: current non-enhanced excitation           */
    float *exc2,                      /* i/o: current enhanced excitation               */
    float *voice_factors,             /* o  : voicing factors                           */
    float *bwe_exc,                   /* o  : excitation for SWB TBE                    */
    float *gain_buf
    , short bfi
);

void ppp_quarter_decoder(
    Decoder_State *st,                        /* i/o: decoder state structure                 */
    DTFS_STRUCTURE *CURRCW_Q_DTFS,             /* i/o: Current CW DTFS                         */
    int   prevCW_lag,                 /* i  : Previous lag                            */
    float *lastLgainD,                /* i/o: Last gain lowband                       */
    float *lastHgainD,                /* i/o: Last gain highwband                     */
    float *lasterbD,                  /* i/o: Last ERB vector                         */
    short bfi,					    /* i  : FER flag */
    DTFS_STRUCTURE PREV_CW_D                   /* i  : Previous DTFS                           */
);

void ppp_voiced_decoder(
    Decoder_State *st,                        /* i/o: state structure                         */
    float *out,                       /* o : residual signal                          */
    const float *lpc2,                      /* i : current frame LPC                        */
    float *exc,                       /* i : previous frame excitation                */
    float *pitch                      /* o : floating pitch values for each subframe  */
    ,short bfi
);

void AVQ_demuxdec(
    Decoder_State *st,                        /* i/o: decoder state structure         */
    int   xriq[],                     /* o:   decoded subvectors [0..8*Nsv-1] */
    short *nb_bits,                   /* i/o: number of allocated bits        */
    const short Nsv,                        /* i:   number of subvectors            */
    short nq_out[]                    /* i/o: AVQ nq index                    */
);

void re8_dec(
    int nq,
    long I,
    int kv[],
    int y[]
);

void re8_decode_base_index(
    int n,                            /* i  : codebook number (*n is an integer defined in {0,2,3,4,..,n_max}) */
    long I,                           /* i  : index of c (pointer to unsigned 16-bit word)                     */
    int *x                            /* o  : point in RE8 (8-dimensional integer vector)                      */
);

void Init_post_filter(
    PFSTAT *pfstat                     /* i  : post-filter state memories              */
);

void nb_post_filt(
    const short L_frame,                    /* i  : frame length                            */
    const short L_subfr,                    /* i  : sub-frame length                        */
    PFSTAT *pfstat,                    /* i/o: Post filter related memories            */
    float *lp_noise,                  /* i/o: long term noise energy                  */
    const float tmp_noise,                  /* i  : noise energy                            */
    float *synth,                     /* i/o: synthesis                               */
    const float *Aq,                        /* i  : LP filter coefficient                   */
    const float *pitch_buf,                 /* i  : Floating pitch   for each subframe      */
    const short coder_type,                 /* i  : coder_type -> deactivated in AUDIO      */
    const short BER_detect,                 /* i  : BER detect flag                         */
    const short disable_hpf                 /* i  : flag to diabled HPF                     */
);

void decod_unvoiced(
    Decoder_State *st,                        /* i/o: decoder static memory                   */
    const float *Aq,                        /* i  : LP filter coefficient                   */
    const short coder_type,                 /* i  : coding type                             */
    float *tmp_noise,                 /* o  : long term temporary noise energy        */
    float *pitch_buf,                 /* o  : floating pitch values for each subframe */
    float *voice_factors,             /* o  : voicing factors                         */
    float *exc,                       /* o  : adapt. excitation exc                   */
    float *exc2,                      /* o  : adapt. excitation/total exc             */
    float *bwe_exc,                   /* i/o: excitation for SWB TBE                  */
    float *gain_buf
);

void decod_tran(
    Decoder_State *st,                        /* i/o: decoder static memory                   */
    const short L_frame,                    /* i  : length of the frame                     */
    const short tc_subfr,                   /* i  : TC subframe index                       */
    const float *Aq,                        /* i  : LP filter coefficient                   */
    const short coder_type,                 /* i  : coding type                             */
    const float Es_pred,                    /* i  : predicted scaled innov. energy          */
    float *pitch_buf,                 /* o  : floating pitch values for each subframe */
    float *voice_factors,             /* o  : voicing factors                         */
    float *exc,                       /* i/o: adapt. excitation exc                   */
    float *exc2,                      /* i/o: adapt. excitation/total exc             */
    float *bwe_exc,                   /* i/o: excitation for SWB TBE                  */
    short *unbits,                    /* i/o: number of unused bits                   */
    const short sharpFlag,                  /* i  : formant sharpening flag                  */
    float *gain_buf
);

void decod_gen_voic(
    Decoder_State *st,                        /* i/o: decoder static memory                   */
    const short L_frame,                    /* i  : length of the frame                     */
    const short sharpFlag,                  /* i  : formant sharpening flag                 */
    const float *Aq,                        /* i  : LP filter coefficient                   */
    const short coder_type,                 /* i  : coding type                             */
    const float Es_pred,                    /* i  : predicted scaled innov. energy          */
    const short do_WI,                      /* i  : FEC fast recovery flag                  */
    float *pitch_buf,                 /* o  : floating pitch   for each subframe      */
    float *voice_factors,             /* o  : voicing factors                         */
    float *exc,                       /* i/o: adapt. excitation exc                   */
    float *exc2,                      /* i/o: adapt. excitation/total exc             */
    float *bwe_exc,                   /* i/o: excitation for SWB TBE                  */
    short *unbits,                    /* i/o: number of unused bits                   */
    float *gain_buf
);

void decod_amr_wb(
    Decoder_State *st,                        /* i/o: decoder static memory                     */
    const float *Aq,                        /* i  : LP filter coefficients                    */
    float *pitch_buf,                 /* o  : floating pitch values for each subframe   */
    float *exc,                       /* i/o: adapt. excitation exc                     */
    float *exc2,                      /* i/o: adapt. excitation/total exc               */
    short *hf_gain,                   /* o  : decoded HF gain                           */
    float *voice_factors,              /* o  : voicing factors                           */
    float *gain_buf
);

void io_ini_dec(
    const int argc,                       /* i  : command line arguments number             */
    char *argv[],                    /* i  : command line arguments                    */
    FILE **f_stream,                 /* o  : input bitstream file                      */
    FILE **f_synth,                  /* o  : output synthesis file                     */
    short *quietMode,                 /* o  : limited printouts                         */
    short *noDelayCmp,                /* o  : turn off delay compensation               */
    Decoder_State *st,                        /* o  : Decoder static variables structure        */
#ifdef SUPPORT_JBM_TRACEFILE
    char **jbmTraceFileName,         /* o  : VOIP tracefilename                        */
#endif
    char **jbmFECoffsetFileName       /* o  : Output file for Optimum FEC offset        */
);

void init_decoder(
    Decoder_State *st                         /* o  : Decoder static variables structure      */
);

void destroy_decoder(
    Decoder_State *st                         /* o  : Decoder static variables structure      */
);

void evs_dec(
    Decoder_State *st,                       /* i/o: Decoder state structure                 */
    float *output,                   /* o  : output synthesis signal                 */
    frameMode frameMode                  /* i  : Decoder frame mode                      */
);

int decodeVoip(
    Decoder_State *st,
    FILE *f_stream,
    FILE *f_synth,
#ifdef SUPPORT_JBM_TRACEFILE
    const char *jbmTraceFileName,
#endif
    const char *jbmFECoffsetFileName, /* : Output file  for Optimum FEC offset        */
    const short quietMode
);

void get_next_frame_parameters(
    Decoder_State *st                       /* i/o: Decoder state structure */
);

void amr_wb_dec(
    Decoder_State *st,                      /* i/o: decoder state structure         */
    float *output_sp                /* o  : synthesis output                */
);

void transf_cdbk_dec(
    Decoder_State *st,                        /* i/o: decoder state structure   */
    const long  core_brate,                 /* i  : core bitrate                                    */
    const short coder_type,                 /* i  : coding type                                     */
    const short harm_flag_acelp,            /* i  : harmonic flag for higher rates ACELP            */
    const short i_subfr,                    /* i  : subframe index                                  */
    const short tc_subfr,                   /* i  : TC subframe index                               */
    const float Es_pred,                    /* i  : predicited scaled innovation energy             */
    const float gain_code,                  /* i  : innovative excitation gain                      */
    float *mem_preemp_preQ,           /* i/o: prequantizer preemhasis memory                  */
    float *gain_preQ,                 /* o  : prequantizer excitation gain                    */
    float *norm_gain_preQ,            /* o  : normalized prequantizer excitation gain         */
    float code_preQ[],                /* o  : prequantizer excitation                         */
    short *unbits                     /* o  : number of AVQ unused bits                       */
);

float gain_dequant(                         /* o:   decoded gain                  */
    short index,                      /* i:   quantization index            */
    const float min,                        /* i:   value of lower limit          */
    const float max,                        /* i:   value of upper limit          */
    const short bits                        /* i:   number of bits to dequantize  */
);

void hq_core_enc(
    Encoder_State *st,                        /* i/o: encoder state structure             */
    const float *audio,                     /* i  : input audio signal                  */
    const short input_frame,                /* i  : frame length                        */
    const short hq_core_type,               /* i  : HQ core type                        */
    const short Voicing_flag
);

short detect_transient(
    const float *in,                        /* i  : input signal                        */
    Encoder_State *st,                        /* i/o: encoder state structure             */
    const short L,                          /* i  : length                              */
    const short coder_type                  /* i  : coding type                         */
);

void wtda(
    const float *new_audio,                 /* i  : input audio                         */
    float *wtda_audio,                /* o  : windowed audio                      */
    float *old_wtda,                  /* i/o: windowed audio from previous frame  */
    const short left_mode,
    const short right_mode,                 /* window overlap of current frame (0: full, 2: none, or 3: half) */
    const short L                           /* i  : length                              */
);


void tcx_get_windows_mode1(
    const short left_mode,                  /* i: overlap mode of left window half          */
    const short right_mode,                 /* i: overlap mode of right window half         */
    float *left_win,                  /* o: left overlap window                       */
    float *right_win,                 /* o: right overlap window                      */
    float *left_win_int,              /* o: left overlap window                       */
    float *right_win_int,             /* o: right overlap window                      */
    short const L
);

void direct_transform(
    const float *in32,                      /* i  : input signal                        */
    float *out32,                     /* o  : output transformation               */
    const short is_transient,               /* i  : transient flag                      */
    const short L                           /* i  : length                              */
);

short noise_adjust(                         /* o  : index of noise attenuation          */
    const float *coeffs_norm,               /* i  : normalized coefficients             */
    const short *bitalloc,                  /* i  : bit allocation                      */
    const short *sfm_start,                 /* i  : Start of bands                      */
    const short *sfm_end,                   /* i  : End of bands                        */
    const short core_sfm                    /* i  : index of the end band for core      */
);

void interleave_spectrum(
    float *coefs,                     /* i/o: input and output coefficients       */
    const short length                      /* i  : length of spectrum                  */
);

void hq_hr_enc(
    Encoder_State *st,                        /* i/o: encoder state structure             */
    float *coefs,                     /* i/o: transform-domain coefficients       */
    const short length,                     /* i  : length of spectrum                  */
    short *num_bits,                  /* i  : number of available bits            */
    const short is_transient                /* i  : transient flag                      */
);

void logqnorm(
    const float *x,                         /* i  : coefficient vector                  */
    short *k,                         /* o  : index                               */
    const short L,                          /* i  : codebook length                     */
    const short N,                          /* i  : sub-vector size                     */
    const float *thren
);

void huff_dec(
    Decoder_State *st,                        /* i/o: decoder state structure                         */
    const short N,                          /* i  : Number of codewords to decode                   */
    const short buffer_len,                 /* i  : Number of bits to read                          */
    const short num_lengths,                /* i  : Number of different huffman codeword lengths    */
    const short *thres,                     /* i  : Threshold of first codeword of each length      */
    const short *offset,                    /* i  : Offset for first codeword                       */
    const short *huff_tab,                  /* i  : Huffman table order by codeword lengths         */
    short *index                      /* o  : Decoded index                                   */
);

void calc_norm(
    const float *x,                         /* i  : Input vector.                       */
    short *norm,                      /* o  : Quantization indices for norms      */
    short *normlg,                    /* o  : Quantized norms in log2             */
    const short start_band,                 /* i  : Indice of band to start coding      */
    const short num_bands,                  /* i  : Number of bands                     */
    const short *band_len,                  /* i  : Length of bands                     */
    const short *band_start                 /* i  : Start of bands                      */
);

void reordernorm(
    const short *ynrm,                      /* i  : quantization indices for norms      */
    const short *normqlg2,                  /* i  : quantized norms                     */
    short *idxbuf,                    /* o  : reordered quantization indices      */
    short *normbuf,                   /* o  : reordered quantized norms           */
    const short nb_sfm                      /* i  : number of bands                     */
);

void diffcod(
    const short N,                          /* i  : number of sub-vectors               */
    short *y,                         /* i/o: indices of quantized norms          */
    short *difidx                     /* o  : differential code                   */
);

void diffcod_lrmdct(
    const short N,
    const int be_ref,
    int *y,
    int *difidx,
    const short is_transient
);

void normalizecoefs(
    float *coefs,                     /* i/o: MDCT coefficients                   */
    const short *ynrm,                      /* i  : quantization indices for norms      */
    const short num_bands,                  /* i  : Number of bands                     */
    const short *band_start,                /* i  : Start of bands                      */
    const short *band_end                   /* i  : End of bands                        */
);

void bitallocsum(
    short *R,                         /* i  : bit-allocation vector               */
    const short nb_sfm,                     /* i  : number of sub-vectors               */
    short *sum,                       /* o  : total number of bits allocated      */
    short *Rsubband,                  /* o  : rate per subband (Q3)               */
    const short  v,                         /* i  : bit rate                            */
    const short  length,                    /* i  : length of spectrum                  */
    const short *sfmsize                    /* i  : Length of bands                     */
);

void hq_generic_hf_encoding(
    const float *coefs,                     /* i  : MDCT coefficients of weighted original      */
    float *hq_generic_fenv,           /* i/o: energy of SWB envelope                      */
    const short hq_generic_offset,          /* i  : frequency offset for extracting energy      */
    Encoder_State *st,                        /* i/o: encoder state structure                     */
    short *hq_generic_clas            /* o  : bwe excitation class                        */
);

short swb_bwe_gain_deq(                     /* o  : BWE class                                   */
    Decoder_State *st,                        /* i/o: decoder state structure                     */
    const short core,                       /* i  : core                                        */
    float *SWB_tenv,                  /* o  : time-domain BWE envelope                    */
    float *SWB_fenv,                  /* o  : frequency-domain BWE envelope               */
    const short hr_flag,                    /* i  : high rate flag                              */
    const short hqswb_clas                  /* i  : HQ BWE class                                */
);

void save_old_syn(
    const short L_frame,                    /* i  : frame length                        */
    const float syn[],                      /* i  : ACELP synthesis                     */
    float old_syn[],                  /* o  : old synthesis buffer                */
    float old_syn_12k8_16k[],         /* i/o: old synthesis buffer                */
    const float preemph_fac,                /* i  : preemphasis factor                  */
    float *mem_deemph                 /* i/o: deemphasis filter memory            */
);

void hq_generic_hf_decoding(
    const short HQ_mode,                    /* i  : HQ mode                                     */
    float *coeff_out1,                /* i/o : BWE input & temporary buffer               */
    const float *hq_generic_fenv,           /* i  : SWB frequency envelopes                     */
    float *coeff_out,                 /* o  : SWB signal in MDCT domain                   */
    const short hq_generic_offset,          /* i  : frequency offset for representing hq swb bwe*/
    short *prev_L_swb_norm,           /* i/o: last normalize length                       */
    const short hq_swb_bwe_exc_clas,         /* i  : bwe excitation class                        */
    const short *R
);

void hq_core_dec(
    Decoder_State *st,                        /* i/o: decoder state structure            */
    float out[],                      /* o  : output synthesis                   */
    const short output_frame,               /* i  : output frame length                */
    const short hq_core_type,               /* i  : HQ core type                       */
    const short core_switching_flag         /* i  : ACELP->HQ switching frame flag     */
);

void hq_hr_dec(
    Decoder_State *st,                        /* i/o: decoder state structure                 */
    float *coefsq,                    /* o  : transform-domain coefficients           */
    const short length,                     /* i  : frame length                            */
    short num_bits,                   /* i  : number of available bits                */
    short *ynrm,                      /* o  : norm quantization index vector          */
    short *is_transient,              /* o  : transient flag                          */
    short *hqswb_clas,                /* o  : HQ SWB class                            */
    float *SWB_fenv                   /* o  : SWB frequency envelopes                 */
);

void hdecnrm_context(
    Decoder_State *st,                        /* i/o: decoder state structure                 */
    const short N,                          /* i  : number of norms                         */
    short *index,                     /* o  : indices of quantized norms              */
    short *n_length                   /* o  : decoded stream length                   */
);

void hdecnrm_tran(
    Decoder_State *st,                        /* i/o: decoder state structure                 */
    const short N,
    short *index
);

void hdecnrm_resize(
    Decoder_State *st,                        /* i/o: decoder state structure                 */
    const short N,                          /* i  :number of SFMs                           */
    short *index                      /* o  : norm quantization index vector          */
);

void hdecnrm(
    Decoder_State *st,                        /* i/o: decoder state structure                 */
    const short N,                          /* i  : number of norms                         */
    short *index                      /* o  : indices of quantized norms              */
);

short find_last_band(                       /* o  : index of last band                      */
    const short *bitalloc,                  /* i  : bit allocation                          */
    const short nb_sfm                      /* i  : number of possibly coded bands          */
);

void fill_spectrum(
    float *coeff,                     /* i/o: normalized MLT spectrum / nf spectrum                 */
    short *R,                         /* i  : number of pulses per band                             */
    const short is_transient,               /* i  : transient flag                                        */
    short norm[],                     /* i  : quantization indices for norms                        */
    const float *hq_generic_fenv,           /* i  : HQ GENERIC envelope                                   */
    const short hq_generic_offset,          /* i  : HQ GENERIC offset                                     */
    const short nf_idx,                     /* i  : noise fill index                                      */
    const short length,                     /* i  : Length of spectrum (32 or 48 kHz)                     */
    const float env_stab,                   /* i  : Envelope stability measure [0..1]                     */
    short *no_att_hangover,           /* i/o: Frame counter for attenuation hangover                */
    float *energy_lt,                 /* i/o: Long-term energy measure for transient detection      */
    short *bwe_seed,                  /* i/o: random seed for generating BWE input                  */
    const short hq_generic_exc_clas,        /* i  : HF excitation class                                  */
    const short core_sfm,                   /* i  : index of the end band for core                        */
    short HQ_mode,                    /* i  : HQ mode                                               */
    float noise_level[],              /* i  : noise level for harmonic modes                        */
    long  core_brate,                 /* i  : target bit-rate                                       */
    float prev_noise_level[],         /* i/o: noise factor in previous frame                        */
    short *prev_R,                    /* i/o: bit allocation info. in previous frame                */
    float *prev_coeff_out,            /* i/o: decoded spectrum in previous frame                    */
    const short *peak_idx,                  /* i  : peak positions                                        */
    const short Npeaks,                     /* i  : number of peaks                                       */
    const short *npulses,                   /* i  : Number of assigned pulses per band                    */
    short prev_is_transient,
    float *prev_normq,
    float *prev_env,
    short prev_bfi,
    const short *sfmsize,                   /* i  : Length of bands                                       */
    const short *sfm_start,                 /* i  : Start of bands                                        */
    const short *sfm_end,                   /* i  : End of bands                                          */
    short *prev_L_swb_norm,
    short prev_hq_mode,
    const short num_sfm,
    const short num_env_bands
);

void env_stab_transient_detect(
    const short is_transient,               /* i:   Transient flag                                        */
    const short length,                     /* i  : Length of spectrum (32 or 48 kHz)                     */
    const short norm[],                     /* i  : quantization indices for norms                        */
    short *no_att_hangover,           /* i/o: Frame counter for attenuation hangover                */
    float *energy_lt,                 /* i/o: Long-term energy measure for transient detection      */
    const short HQ_mode,                    /* i  : HQ coding mode                                        */
    const short bin_th,                     /* i  : HVQ cross-over frequency bin                          */
    const float *coeff                      /* i  : Coded spectral coefficients                           */
);

void de_interleave_spectrum(
    float *coefs,                     /* i/o: input and output coefficients   */
    short length                      /* i  : length of spectrum              */
);

void inverse_transform(
    const float *InMDCT,                    /* i  : input MDCT vector               */
    float *Out,                       /* o  : output vector                   */
    const short IsTransient,                /* i  : transient flag                  */
    const short L,                          /* i  : output frame length             */
    const short L_inner                     /* i  : length of the transform         */
);

void window_ola(
    const float *ImdctOut,                  /* i  : input                           */
    float *auOut,                     /* o  : output audio                    */
    float *OldauOut,                  /* i/o: audio from previous frame       */
    const short L,                          /* i  : length                          */
    const short left_mode,
    const short right_mode, /* window overlap of current frame (0: full, 2: none, or 3: half) */
    const short use_bfi_win,                /* i  : use BFI windowing               */
    const short oldHqVoicing,               /* i  : previous HqVoicing              */
    float *oldgapsynth                /* i  : previous gapsynth               */
);

void map_quant_weight(
    const short normqlg2[],                 /* i  : quantized norms                 */
    short wnorm[],                    /* o  : weighted norm                   */
    const short is_transient                /* i  : transient flag                  */
);

void recovernorm(
    short *idxbuf,                    /* i  : reordered quantization indices  */
    short *ynrm,                      /* o  : recovered quantization indices  */
    short *normqlg2,                  /* o  : recovered quantized norms       */
    short nb_sfm                      /* i  : number of SFMs                  */
);

void reordvct(
    short *y,                         /* i/o: vector to rearrange             */
    short N,                          /* i  : dimensions                      */
    short *idx                        /* o  : reordered vector index          */
);

void bitalloc(
    short *y,                         /* i  : reordered norm of sub-vectors           */
    short *idx,                       /* i  : reordered sub-vector indices            */
    short sum,                        /* i  : number of available bits                */
    short N,                          /* i  : number of norms                         */
    short K,                          /* i  : maximum number of bits per dimension    */
    short *r,                         /* o  : bit-allacation vector                   */
    const short *sfmsize,                   /* i  : Length of bands                         */
    const short hqswb_clas                  /* i  : signal classification flag              */
);

short BitAllocF(
    short *y,                         /* i  : norm of sub-vectors                     */
    long  bit_rate,                   /* i  : bitrate                                 */
    short B,                          /* i  : number of available bits                */
    short N,                          /* i  : number of sub-vectors                   */
    short *R,                         /* o  : bit-allocation indicator                */
    short *Rsubband,                  /* o  : sub-band bit-allocation vector (Q3)     */
    const short hqswb_clas,           /* i  : hq swb class                            */
    const short num_env_bands         /* i  : Number sub bands to be encoded for HQ_SWB_BWE  */
);

short BitAllocWB(
    short *y,
    short B,
    short N,
    short *R,
    short *Rsubband
);

short hvq_pvq_bitalloc(
    short num_bits,               /* i/o: Number of available bits (including gain bits) */
    const short brate,                  /* i  : bitrate                     */
    const short bwidth,                 /* i  : Encoded bandwidth           */
    const short *ynrm,                  /* i  : Envelope coefficients       */
    const int   manE_peak,              /* i  : Peak energy mantissa        */
    const short expE_peak,              /* i  : Peak energy exponent        */
    short *Rk,                    /* o  : bit allocation for concatenated vector */
    short *R,                     /* i/o: Global bit allocation       */
    short *sel_bands,             /* o  : Selected bands for encoding */
    short *n_sel_bands            /* o  : No. of selected bands for encoding */
);

void floating_point_add(
    int   *mx,                        /* io: mantissa of the addend Q31 */
    short *ex,                        /* io: exponent of the addend Q0  */
    const int    my,                        /* i:  mantissa of the adder Q31  */
    const short  ey                         /* i:  exponent of the adder Q0   */
);

short rc_get_bits2(                         /* o: Number of bits needed         */
    const short N,                          /* i: Number of bits currently used */
    const unsigned int range                /* i: Range of range coder          */
);



void rc_enc_init(
    Encoder_State *st,                      /* i/o: Encoder state       */
    short tot_bits                          /* i  : Total bit budget    */
);

void rc_encode(
    Encoder_State *st,                      /* i/o: Encoder state                       */
    unsigned int cum_freq,                  /* i  : Cumulative frequency up to symbol   */
    unsigned int sym_freq,                  /* i  : Symbol probability                  */
    unsigned int tot                        /* i  : Total cumulative frequency          */
);

void rc_enc_finish(
    Encoder_State *st                       /* i/o: Encoder state       */
);

void rc_enc_bits(
    Encoder_State *st,                      /* i/o: Encoder state       */
    unsigned int value,                     /* i  : Value to encode     */
    short bits                              /* i  : Number of bits used */
);

void rc_enc_uniform(
    Encoder_State *st,                      /* i/o: Encoder state       */
    unsigned int value,                     /* i  : Value to encode     */
    unsigned int tot                        /* i  : Maximum value       */
);

void rc_dec_init(
    Decoder_State *st,                      /* i/o: Decoder State       */
    short tot_bits                          /* i  : Total bit budget    */
);

unsigned int rc_decode(                     /* o  : Decoded cumulative frequency    */
    Decoder_State *st,                      /* i/o: Decoder State                   */
    unsigned int tot                        /* i  : Total cumulative frequency      */
);

void rc_dec_update(
    Decoder_State *st,                      /* i/o: Decoder State           */
    unsigned int cum_freq,                  /* i  : Cumulative frequency    */
    unsigned int sym_freq                   /* i  : Symbol frequency        */
);

unsigned int rc_dec_bits(                   /* i  : Decoded value   */
    Decoder_State *st,                      /* i/o: Decoder State   */
    short bits                              /* i  : Number of bits  */
);

unsigned int rc_dec_uniform(                /* i  : Decoded value   */
    Decoder_State *st,                      /* i/o: Decoder State   */
    unsigned int tot                        /* i  : Maximum value   */
);

void rc_dec_finish(
    Decoder_State *st                         /* i/o: deocder state structure */
);

short pvq_core_enc (
    Encoder_State *st,
    float coefs_norm[],
    float coefs_quant[],
    short bit_budget,                 /* number of bits */
    short nb_sfm,
    const short *sfm_start,
    const short *sfm_end,
    const short *sfmsize,
    short *R,
    short *Rs,
    short *npulses,
    short *maxpulse,
    const short core
);

short pvq_core_dec (
    Decoder_State *st,
    const short *band_start,
    const short *band_end,
    const short *band_width,
    float coefs_quant[],              /* o  : output   MDCT     */
    short bits_tot,
    short nb_sfm,
    short *R,
    short *Rs,
    short *npulses,
    short *maxpulse,
    const short core
);

void pvq_encode(
    Encoder_State *st,
    const float *x,                         /* i:   vector to quantize               */
    short *y,                               /* o:   quantized vector (non-scaled int)*/
    float *xq,                        /* o:   quantized vector (scaled float)  */
    const short pulses,                     /* i:   number of allocated pulses       */
    const short N,                          /* i:   Length of vector                 */
    const float gain                        /* i:   Gain                             */
);

void pvq_decode(
    Decoder_State *st,
    float *xq,                        /* o:   decoded vector (scaled float)    */
    short *y,                         /* o:   decoded vector (non-scaled short)*/
    const short K,                          /* i:   number of allocated pulses       */
    const short N,                          /* i:   Length of vector                 */
    const float gain                        /* i:   Gain                             */
);


void rangeCoderFinalizationFBits(
    short Brc,
    unsigned int INTrc,
    short *FBits
);

void bandBitsAdjustment(
    short Brc,
    unsigned int  INTrc,
    short  Bavail,
    short Nbands,
    short D,
    short L,
    short Bband,
    short Breserv,
    short *Bband_adj,  /* o */
    short *Brem,       /* o */
    short *Breservplus /* o */
) ;
void conservativeL1Norm(
    short L,
    short Qvec,
    short Fcons,
    short Qavail,
    short Qreserv,
    short Dspec,
    short *Dvec,        /*o*/
    short *Qspare,      /*o*/
    short *Qreservplus, /*o*/
    short *Dspecplus    /*o*/
);

void NearOppSplitAdjustment(
    short qband,
    short qzero,
    const short Qac,
    const unsigned int  INTac,
    short qglobal,
    short FlagCons,
    short Np,
    short Nhead,
    short Ntail,
    short Nnear,
    short Nopp,
    short oppRQ3,
    short *qnear,
    short *qopp,
    short *qglobalupd
) ;
void   dsDirac2Dirac(
    short dsDiracIndex,
    short *diracs
);
int intLimCDivPos(
    int  NUM,
    short DEN
) ;

int   intLimCDivSigned(
    int   NUM,
    short DEN
) ;


short  shrtCDivSignedApprox(
    const short  num,
    const short  den
);

void dsDiracPerQuanta(short  td,
                      short  t_quanta,
                      short  dsm,
                      const unsigned char* const *frQuanta,
                      short *DsIdx
                     );
void    QuantaPerDsDirac(
    short td,
    short  dsDiracIndex,
    const unsigned char* const*  dimFrQuanta,
    short *Quanta
) ;

void   Ratio_rQ3(
    short opp,
    short near,
    short *result
);

void  obtainEnergyQuantizerDensity(
    short L,
    short R,
    short *Density
) ;

void densityAngle2RmsProjDec(
    short D,         /* i: density */
    short indexphi,  /* i: decoded  index from AR dec */
    short *oppQ15,   /* o: proj*/
    short *nearQ15,   /* o: */
    short *oppRatioQ3 /* o: */
) ;

void densityAngle2RmsProjEnc (
    short D,
    short phiQ14uq,
    short *indexphi, /* o: index to send to AR */
    short *oppQ15,   /* o: proj*/
    short *nearQ15,
    short *oppRatioQ3
) ;


void subband_gain_bits(
    const short *Rk,                        /* i  : bit allocation per band (Q3)*/
    const short N,                          /* i  : number of bands         */
    short *bits,                      /* o  : gain bits per band      */
    const short *sfmsize                    /* i  : Size of bands           */
);

void env_adj(
    const short *pulses,                    /* i  : number of pulses per band           */
    const short length,                     /* i  : length of spectrum                      */
    const short last_sfm,                   /* i  : Index of last band                      */
    float *adj,                       /* o  : Adjustment factors for the envelope     */
    const float env_stab,
    const short *sfmsize                    /* i  : Length of bands                         */
);

float env_stability(
    const short *ynrm,                      /*i  : Norm vector for current frame          */
    const short nb_sfm,                     /*i  : Number of sub-bands                    */
    short *mem_norm,                        /*i/o: Norm vector memory from past frame     */
    float *mem_env_delta                    /*i/o: Envelope stability memory for smoothing*/
);

float env_stab_smo(
    float env_stab,                         /*i  : env_stab value                         */
    float *env_stab_state_p,                /*i/o: env_stab state probabilities           */
    short *ho_cnt                           /*i/o: hangover counter for speech state      */
);

void core_switching_pre_enc(
    Encoder_State *st,                      /* i/o: encoder state structure                 */
    LPD_state *mem,                     /* i/o: encoder state structure                 */
    const float *old_inp_12k8,            /* i  : old input signal @12.8kHz               */
    const float *old_inp_16k              /* i  : old input signal @16kHz                 */
);

void core_switching_post_enc(
    Encoder_State *st,                      /* i/o: encoder state structure                 */
    const float *old_inp_12k8,            /* i  : old input signal @12.8kHz               */
    const float *old_inp_16k,             /* i  : old input signal @16kHz                 */
    const short pitch[3],                 /* i  : open-loop pitch values for quantiz.     */
    const float voicing[3],               /* i  : Open-loop pitch gains                   */
    const float A[]                       /* i  : unquant LP filter coefs.                */
);

void core_switching_post_dec(
    Decoder_State *st,                      /* i/o: decoder state structure                 */
    float *synth,                   /* i/o: output synthesis                        */
    const short output_frame,             /* i  : frame length                            */
    const short core_switching_flag,      /* i  : ACELP->HQ switching frame flag          */
    const short coder_type
);

void core_switching_pre_dec(
    Decoder_State *st,                       /* i/o: decoder state structure                 */
    const short output_frame               /* i  : frame length                            */
);

void bandwidth_switching_detect(
    Decoder_State *st                         /* i/o: encoder state structure                 */
);

void bw_switching_pre_proc(
    Decoder_State *st,                       /* i/o: decoder state structure                 */
    const float *old_syn_12k8_16k          /* i  : ACELP core synthesis @ 12.8kHz or 16kHz */
);

void updt_bw_switching(
    Decoder_State *st,                       /* i/o: decoder state structure                 */
    const float *synth,                    /* i  : float synthesis signal                  */
    const short *inner_frame_tbl           /* i  : HQ inner_frame signallisation table     */
);

void swb_tbe_reset(
    float mem_csfilt[],
    float mem_genSHBexc_filt_down_shb[],
    float state_lpc_syn[],
    float syn_overlap[],
    float state_syn_shbexc[],
    float *tbe_demph,
    float *tbe_premph,
    float mem_stp_swb[],
    float *gain_prec_swb
);

void swb_tbe_reset_synth(
    float genSHBsynth_Hilbert_Mem[],
    float genSHBsynth_state_lsyn_filt_shb_local[]
);


void fb_tbe_reset_enc(
    float elliptic_bpf_2_48k_mem[][4],
    float *prev_fb_energy
);

void fb_tbe_reset_synth(
    float fbbwe_hpf_mem[][4],
    float *prev_fbbwe_ratio
);

void wb_tbe_extras_reset(
    float mem_genSHBexc_filt_down_wb2[],
    float mem_genSHBexc_filt_down_wb3[]
);

void wb_tbe_extras_reset_synth(
    float state_lsyn_filt_shb[],
    float state_lsyn_filt_dwn_shb[],
    float mem_resamp_HB[]
);

void prep_tbe_exc(
    const short L_frame,                   /* i  : length of the frame                     */
    const short i_subfr,                   /* i  : subframe index                          */
    const float gain_pit,                  /* i  : Pitch gain                              */
    const float gain_code,                 /* i  : algebraic codebook gain                 */
    const float code[],                    /* i  : algebraic excitation                    */
    const float voice_fac,                 /* i  : voicing factor                          */
    float *voice_factors,            /* o  : TBE voicing factor                      */
    float bwe_exc[],                 /* i/o: excitation for TBE                      */
    const float gain_preQ,                 /* i  : prequantizer excitation gain            */
    const float code_preQ[],               /* i  : prequantizer excitation                 */
    const short T0,                        /* i  : integer pitch variables                 */
    const short coder_type,                /* i  : coding type                             */
    const long  core_brate                 /* i  : core bitrate                            */
);

void synthesise_fb_high_band(
    const float excitation_in[],            /* i  : full band excitation                    */
    float output[],                   /* o  : high band speech - 14.0 to 20 kHz       */
    const float fb_exc_energy,              /* i  : full band excitation energy             */
    const float ratio,                      /* i  : energy ratio                            */
    const short L_frame,                    /* i  : ACELP frame length                      */
    const short bfi,                        /* i  : fec flag                                */
    float *prev_fbbwe_ratio,          /* o  : previous frame energy for FEC           */
    float bpf_memory[][4]             /* i/o: memory for elliptic bpf 48k             */
);

void elliptic_bpf_48k_generic(
    const float input[],                   /* i  : input signal                            */
    float output[],                  /* o  : output signal                           */
    float memory[][4],               /* i/o: 4 arrays for memory                     */
    const float full_band_bpf[][5]         /* i  : filter coefficients b0,b1,b2,a0,a1,a2   */
);

void HQ_FEC_processing(
    Decoder_State *st,                       /* i/o: decoder state structure                          */
    float *t_audio_q,                        /* o  : MDCT coeffs. (for synthesis)                     */
    short is_transient,                      /* i  : Old flag for transient                           */
    float ynrm_values[][MAX_PGF],            /* i  : Old average Norm values for each group of bands  */
    float r_p_values[][MAX_ROW],             /* i  : Computed y-intercept and slope by Regression     */
    short num_Sb,                            /* i  : Number of sub-band group                         */
    short nb_sfm,                            /* i  : Number of sub-band                               */
    short *Num_bands_p,                      /* i  : Number of coeffs. for each sub-band              */
    short output_frame,                      /* i  : Frame size                                       */
    const short *sfm_start,                  /* i  : Start of bands                                   */
    const short *sfm_end                     /* i  : End of bands                                     */
);

void HQ_FEC_Mem_update(
    Decoder_State *st,                          /* i/o: decoder state structure            */
    float *t_audio_q,
    float *normq,
    short *ynrm,
    short *Num_bands_p,
    short is_transient,
    short hqswb_clas,
    short c_switching_flag,
    short nb_sfm,
    short num_Sb,
    float *mean_en_high,
    short hq_core_type,                  /* i : normal or low-rate MDCT(HQ) core */
    short output_frame
);

void time_domain_FEC_HQ(
    Decoder_State *st,                           /* i  : Decoder State                           */
    float *wtda_audio,                  /* i  : input                                   */
    float *out,                         /* o  : output audio                            */
    float mean_en_high,                 /* i  : transient flag                          */
    const short output_frame
);

void Next_good_after_burst_erasures(
    const float *ImdctOut,                  /* i  : input                                   */
    float *auOut,                     /* o  : output audio                            */
    float *OldauOut,                  /* i/o: audio from previous frame               */
    const short ol_size                     /* i  : overlap size                            */
);

void update_average_rate(
    Encoder_State *st                         /* i/o: encoder state structure                */
);

void reset_preecho_dec(
    Decoder_State *st                         /* i/o: decoder state structure   */
);

void preecho_sb(
    const long  brate,                       /* i  : core bit-rate                                           */
    const float wtda_audio[],                /* i  : imdct signal                                            */
    float *rec_sig,                    /* i  : reconstructed signal, output of the imdct transform     */
    const short output_frame,                /* i  : output frame length                                     */
    float *memfilt_lb,                 /* i/o: memory                                                  */
    float *mean_prev_hb,               /* i/o: memory                                                  */
    float *smoothmem,                  /* i/o: memory                                                  */
    float *mean_prev,                  /* i/o: memory                                                  */
    float *mean_prev_nc,               /* i/o: memory                                                  */
    float *wmold_hb,                   /* i/o: memory                                                  */
    short *prevflag,                   /* i/o: flag                                                    */
    short *pastpre,                    /* i/o: flag                                                    */
    const short bwidth
);

void hq2_core_configure(
    const short frame_length,                /* i  : frame length                            */
    const short num_bits,                    /* i  : number of bits                          */
    const short is_transient,                /* i  : transient flag                          */
    short *bands,
    short *length,
    short band_width[],
    short band_start[],
    short band_end[],
    Word32 *L_qint,                    /* o  : Q29 */
    Word16 *eref_fx,                   /* o  : Q10 */
    Word16 *bit_alloc_weight_fx,       /* o  : Q13 */
    short *gqlevs,
    short *Ngq,
    short *p2a_bands,
    float *p2a_th,
    float *pd_thresh,
    float *ld_slope,
    float *ni_coef,
    float *ni_pd_th,
    long  bwe_br
);

void hq_lr_enc(
    Encoder_State *st,                          /* i/o: encoder state structure                 */
    float t_audio[],                   /* i/o: transform-domain coefs.                 */
    const short inner_frame,                 /* i  : inner frame length                      */
    short *num_bits,                   /* i/o: number of available bits                */
    const short is_transient                 /* i  : transient flag                          */
);

void hq_lr_dec(
    Decoder_State *st,                         /* i/o: decoder state structure                 */
    float yout[],                      /* o  : transform-domain output coefs.          */
    const short inner_frame,                 /* i  : inner frame length                      */
    short num_bits,                    /* i  : number of available bits                */
    short *is_transient                /* o  : transient flag                          */
);

void hq2_bit_alloc (
    const float  band_energy[],              /* i  : band energy of each subband                 */
    const short  bands,                      /* i  : total number of subbands in a frame         */
    Word32  L_Rk[],                    /* i/o: Bit allocation/Adjusted bit alloc.          */
    short  *bit_budget,                /* i/o: bit bugdet                                  */
    short  *p2a_flags,                 /* i  : HF tonal indicator                          */
    const Word16 weight_fx,                  /* i  : weight (Q13)                                */
    const short  band_width[],               /* i  : Sub band bandwidth                          */
    const short  num_bits,                   /* i  : available bits                              */
    const short  hqswb_clas,                 /* i  : HQ2 class information                       */
    const short  bwidth,                     /* i  : input bandwidth                             */
    const short  is_transient                /* i  : indicator HQ_TRANSIENT or not               */
);

void hq2_noise_inject(
    float y2hat[],
    const short band_start[],
    const short band_end[],
    const short band_width[],
    float Ep[],
    float Rk[],
    const int   npulses[],
    short ni_seed,
    const short bands,
    const short ni_start_band,
    const short bw_low,
    const short bw_high,
    const float enerL,
    const float enerH,
    float last_ni_gain[],
    float last_env[],
    short *last_max_pos_pulse,
    short *p2a_flags,
    short p2a_bands,
    const short hqswb_clas,
    const short bwidth,
    const long  bwe_br
);

void mdct_spectrum_denorm(
    const int   inp_vector[],
    float y2[],
    const short band_start[],
    const short band_end[],
    const short band_width[],
    const float band_energy[],
    const int   npulses[],
    const short bands,
    const float ld_slope,
    const float pd_thresh
);

void reverse_transient_frame_energies(
    float band_energy[],               /* o  : band energies                           */
    const short bands                        /* i  : number of bands                         */
);

short peak_vq_enc(
    Encoder_State *st,                         /* i/o: encoder state structure     */
    const float *coefs,                      /* i  : Input coefficient vector    */
    float *coefs_out,                  /* o  : Quantized output vector     */
    const short brate,                       /* i  : Core bitrate                */
    const short num_bits,                    /* i  : Number of bits for HVQ      */
    const short vq_peaks,                    /* i  : Number of identified peaks  */
    const short *ynrm,                       /* i  : Envelope coefficients       */
    short *R,                          /* i/o: Bit allocation/updated bit allocation */
    short *vq_peak_idx,                /* i  : Peak index vector           */
    float *nf_gains                    /* i  : Estimated noise floor gains */
);

void hvq_dec(
    Decoder_State *st,                         /* i/o: decoder state structure   */
    const short num_bits,                    /* i : Number of available bits       */
    const long  core_brate,                  /* i : Core bit-rate                  */
    const short *ynrm,                       /* i : Envelope coefficients        */
    short *R,                          /* i/o: Bit allocation/updated bit allocation */
    float *noise_level,                /* o : Noise level                    */
    short *peak_idx,                   /* o : Peak position vector           */
    short *Npeaks,                     /* o : Total number of peaks          */
    float *coefsq_norm,                /* o : Output vector                  */
    const short core
);

void hq_configure_bfi(
    short *nb_sfm,                     /* o  : Number of sub bands               */
    short *num_Sb,                     /* o  : Number of FEC sub bands ?         */
    short *num_bands_p,                /* o  : FEC sub bands                     */
    short const **sfmsize,             /* o  : Subband bandwidths                */
    short const **sfm_start,           /* o  : Subband start coefficients        */
    short const **sfm_end              /* o  : Subband end coefficients          */
);

void peak_vq_dec(
    Decoder_State *st,                         /* i/o: decoder state structure     */
    float *coefs_out,                  /* o  : Output coefficient vetor    */
    const short brate,                       /* i  : Core bitrate                */
    const short num_bits,                    /* i  : Number of bits for HVQ      */
    const short *ynrm,                       /* i  : Envelope coefficients       */
    short *R,                          /* i/o: Bit allocation/updated bit allocation */
    short *vq_peak_idx,                /* o  : Peak position vector        */
    short *Npeaks,                     /* o  : Number of peaks             */
    const short core                         /* i  : Core type                   */
);

void swb_bwe_enc_lr(
    Encoder_State *st,                         /* i/o: encoder state structure     */
    const float m_core[],                    /* i  : core synthesis (MDCT)       */
    const float m_orig[],                    /* i/o: scaled orig signal (MDCT)   */
    float m[],                         /* o  : output, SWB part (MDCT)     */
    const long  total_brate,                 /* i  : total bitrate for selecting subband pattern */
    short BANDS,
    short *band_start,
    short *band_end,
    float *band_energy,
    short *p2a_flags,
    const short hqswb_clas,
    short lowlength,
    short highlength,
    short *prev_frm_index,
    const short har_bands,
    short *prev_frm_hfe2,
    short *prev_stab_hfe2,
    short band_width[],
    const float y2_ni[],
    short *ni_seed
);

void swb_bwe_dec_lr(
    Decoder_State *st,                         /* i/o: decoder state structure                     */
    const float m_core[],                    /* i  : lowband synthesis                           */
    float m[],                         /* o  : highband synthesis with lowband zeroed      */
    const long  total_brate,                 /* i  : total bitrate for selecting subband pattern */
    short BANDS,
    short *band_start,
    short *band_end,
    float *band_energy,
    short *p2a_flags,
    const short hqswb_clas,
    short lowlength,
    short highlength,
    const short har_bands,
    short *prev_frm_hfe2,
    short *prev_stab_hfe2,
    short band_width[],
    const float  y2_ni[],
    short *ni_seed
);

int get_usebit_npswb(
    short hqswb_clas
);

void GetPredictedSignal(
    const float *predBuf,                    /* i  : prediction buffer        */
    float *outBuf,                     /* o  : output buffer            */
    const short lag,                         /* i  : prediction buffer offset */
    const short fLen,                        /* i  : length of loop (output)  */
    const float gain                         /* i  : gain to be applied       */
);

void convert_lagIndices_pls2smp(
    short lagIndices_in[],
    short nBands_search,
    short lagIndices_out[],
    const float sspectra[],
    const short sbWidth[],
    const short fLenLow
);

void FindNBiggest2_simple(
    const float *inBuf,                      /* i  : input buffer (searched)                     */
    GainItem *g,                       /* o  : N biggest components found                  */
    const short nIdx,                        /* i  : search length                               */
    short *n,                          /* i  : number of components searched (N biggest)   */
    short N_NBIGGESTSEARCH
);

void updat_prev_frm(
    float y2[],
    float t_audio[],
    long bwe_br,
    short length,
    const short inner_frame,
    short bands,
    short bwidth,
    const short is_transient,
    short hqswb_clas,
    short *prev_hqswb_clas,
    short *prev_SWB_peak_pos,
    short prev_SWB_peak_pos_tmp[],
    short *prev_frm_hfe2,
    short *prev_stab_hfe2,
    short bws_cnt
);

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
);

float spectrumsmooth_noiseton(
    float spectra[],
    const float spectra_ni[],
    float sspectra[],
    float sspectra_diff[],
    float sspectra_ni[],
    const short fLenLow,
    short *ni_seed
);

void noiseinj_hf(
    float xSynth_har[],
    float th_g[],
    float band_energy[],
    float *prev_En_sb,
    const short p2a_flags[],
    short BANDS,
    short band_start[],
    short band_end[],
    const short fLenLow
);

void noise_extr_corcod(
    float spectra[],
    const float spectra_ni[],
    float sspectra[],
    float sspectra_diff[],
    float sspectra_ni[],
    const short fLenLow,
    short prev_hqswb_clas,
    float *prev_ni_ratio
);

void genhf_noise(
    float noise_flr[],
    float xSynth_har[],
    float *predBuf,
    short bands,                     /* i  : total number of subbands in a frame         */
    short harmonic_band,             /* i  : Number of LF harmonic frames                */
    short har_freq_est2,
    short pos_max_hfe2,
    short *pul_res,
    GainItem pk_sf[],
    const short fLenLow,
    const short fLenHigh,
    const short sbWidth[],
    const short lagIndices[],
    const short subband_offsets[],
    const short subband_search_offset[]
);

void ton_ene_est(
    float xSynth_har[],
    float be_tonal[],
    float band_energy[],
    short band_start[],
    short band_end[],
    short band_width[],
    const short fLenLow,
    const short fLenHigh,
    short bands,
    short har_bands,
    float ni_lvl,
    GainItem pk_sf[],
    short *pul_res
);

void Gettonl_scalfact(
    float *outBuf,                   /* o  : synthesized spectrum                        */
    const float *codbuf,                   /* i  : core coder                                  */
    const short fLenLow,                   /* i  : lowband length                              */
    const short fLenHigh,                  /* i  : highband length                             */
    short harmonic_band,             /* i  : Number of LF harmonic frames                */
    short bands,                     /* i  : total number of subbands in a frame         */
    float *band_energy,              /* i  : band energy of each subband                 */
    short *band_start,               /* i  : subband start indices                       */
    short *band_end,                 /* i  : subband end indices                         */
    const short p2aflags[],
    float be_tonal[],
    GainItem *pk_sf,
    short *pul_res
);

void SpectrumSmoothing_nss(
    float *inBuf,
    float *outBuf,
    int fLen
);

void SpectrumSmoothing(
    float *inBuf,
    float *outBuf,
    const short fLen,
    const float th_cut
);

void hq2_bit_alloc_har (
    float  *y,                          /* i  : band energy of sub-vectors                    */
    int    B,                           /* i  : number of available bits                      */
    short  N,                           /* i  : number of sub-vectors                         */
    Word32 *L_Rsubband,
    short  p2a_bands,
    long   core_brate,                  /* i  : core bit rate                                 */
    short  p2a_flags[],
    short  band_width[]
);

void GetSynthesizedSpecThinOut(
    const float *predBuf,
    float       *outBuf,
    const short nBands,
    const short *sbWidth,
    const short *lagIndices,
    const float *lagGains,
    const short predBufLen
);

void return_bits_normal2(
    short *bit_budget,
    const short p2a_flags[],
    const short bands,
    const short bits_lagIndices[]
);

void GetlagGains(
    const float *predBuf,
    const float *band_energy,
    const short nBands,
    const short *sbWidth,
    const short *lagIndices,
    const short predBufLen,
    float *lagGains
);

void preset_hq2_swb(
    const short hqswb_clas,
    const short band_end[],
    short *har_bands,
    short  p2a_bands,
    const short length,
    const short bands,
    short *lowlength,
    short *highlength,
    float m[]
);

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
);

void har_denorm_pulcnt(
    float spectra[],                   /* i/o: MDCT domain spectrum                        */
    short band_start[],                /* i  : Number subbands/Frame                       */
    short band_end[],                  /* i  : Band Start of each SB                       */
    float band_energy[],               /* i  : Band end of each SB                         */
    short band_width[],
    int   npulses[],
    const short har_bands                    /* i: No. of harmonic bands                         */
);

short har_est(
    float spectra[],
    short N,
    short *har_freq_est1,
    short *har_freq_est2,
    short *flag_dis,
    short *prev_frm_hfe2,
    const short subband_search_offset[],
    const short sbWidth[],
    short *prev_stab_hfe2
);

void spt_shorten_domain_pre(
    const short band_start[],
    const short band_end[],
    const short prev_SWB_peak_pos[],
    const short BANDS,
    const long  bwe_br,
    short new_band_start[],
    short new_band_end[],
    short new_band_width[]
);

void spt_shorten_domain_band_save(
    const short bands,
    const short band_start[],
    const short band_end[],
    const short band_width[],
    short org_band_start[],
    short org_band_end[],
    short org_band_width[]
);

void spt_shorten_domain_band_restore(
    const short bands,
    short band_start[],
    short band_end[],
    short band_width[],
    const short org_band_start[],
    const short org_band_end[],
    const short org_band_width[]
);

void spt_swb_peakpos_tmp_save(
    const float y2[],
    const short bands,
    const short band_start[],
    const short band_end[],
    short prev_SWB_peak_pos_tmp[]
);

void hq_ecu(
    const float *prevsynth,                 /* i  : buffer of previously synthesized signal   */
    float *ecu_rec,                   /* o  : reconstructed frame in tda domain         */
    short *time_offs,                 /* i/o: Sample offset for consecutive frame losses*/
    float *X_sav,                     /* i/o: Stored spectrum of prototype frame        */
    short *num_p,                     /* i/o: Number of identified peaks                */
    short *plocs,                     /* i/o: Peak locations                            */
    float *plocsi,                    /* i/o: Interpolated peak locations               */
    const float env_stab,                   /* i  : Envelope stability parameter              */
    short *last_fec,                  /* i/o: Flag for usage of pitch dependent ECU     */
    const short ph_ecu_HqVoicing,           /* i  : HQ Voicing flag                           */
    short *ph_ecu_active,             /* i  : Phase ECU active flag                     */
    float *gapsynth,                  /* o  : Gap synthesis                             */
    const short prev_bfi,                   /* i  : indicating burst frame error              */
    const short old_is_transient[2],        /* i  : flags indicating previous transient frames*/
    float *mag_chg_1st,               /* i/o: per band magnitude modifier for transients*/
    float Xavg[LGW_MAX],              /* i/o: Frequency group average gain to fade to   */
    float *beta_mute,                 /* o   : Factor for long-term mute                */
    const short output_frame,               /* i  : frame length                              */
    Decoder_State *st                 /* i/o: decoder state structure                   */
);

void hq_timedomain_conc(
    float *ecu_rec,                       /* o  : reconstructed frame in tda domain         */
    float *gapsynth,                      /* o  : Gap synthesis                             */
    const short output_frame,                   /* i  : frame length                              */
    const float *prevsynth,                     /* i  : buffer of previously synthesized signal   */
    Decoder_State *st                     /* i/o: decoder state structure                   */
);

void fft3(
    const float X[],                         /* i : input frame                                */
    float Y[],                         /* o : DFT of input frame                         */
    const short n                            /* i : block length (must be radix 3)             */
);

void ifft3(
    const float X[],                         /* i : input frame                                */
    float Y[],                         /* o : iDFT of input frame                        */
    const short n                            /* i : block length (must be radix 3)             */
);

void minimumStatistics(                     /* return: updated estimate of background noise */
    float*      noiseLevelMemory,           /* internal state */
    int*        noiseLevelIndex,            /* internal state */
    int*        currLevelIndex,             /* internal state (circular buffer) */
    float*      noiseEstimate,              /* previous estimate of background noise */
    float*      lastFrameLevel,             /* level of the last frame */
    float       currentFrameLevel,          /* level of the current frame */
    float const minLev,                     /* minimum level */
    int   const buffSize                    /* buffer size */
);

void E_LPC_int_lpc_tcx(
    const float lsf_old[],                  /* input : LSFs from past frame              */
    const float lsf_new[],                  /* input : LSFs from present frame           */
    float a[]                               /* output: interpolated LP coefficients      */
);

Word32 E_GAIN_closed_loop_search(
    Float32 exc[],
    Float32 xn[],
    Float32 h[],
    Word32 t0_min,
    Word32 t0_min_frac,
    Word32 t0_max,
    Word32 t0_max_frac,
    Word32 t0_min_max_res,
    Word32 *pit_frac,
    Word32 *pit_res,
    Word32 pit_res_max,
    Word32 i_subfr,
    Word32 pit_min,
    Word32 pit_fr2,
    Word32 pit_fr1,
    Word32 L_subfr
);

void E_ACELP_toeplitz_mul(
    float R[],
    float c[],
    float d[]
);

void acelp_pulsesign(
    const float cn[],
    float dn[],
    float dn2[],
    float sign[],
    float vec[],
    float alp
);

void E_ACELP_4tsearch(
    Float32 dn[],
    const Float32 cn[],
    const Float32 H[],
    float code[],
    PulseConfig *config,
    Word16 ind[],
    Float32 y[]
);

void E_ACELP_4tsearchx(
    Float32 dn[],
    const Float32 cn[],
    Float32 Rw[],
    float code[],
    PulseConfig *config,
    Word16 ind[]
);

short E_ACELP_indexing(
    Float32 code[],
    PulseConfig config,
    int num_tracks,
    int prm[]
);

void acelp_findcandidates(
    float dn2[],
    short dn2_pos[],
    short pos_max[],
    int L_subfr,
    int tracks
);

void E_ACELP_innovative_codebook(
    float *exc,        /* i  : pointer to the excitation frame                  */
    int T0,          /* i  : integer pitch lag                                */
    int T0_frac,     /* i  : fraction of lag                                  */
    int T0_res,      /* i  : pitch resolution                                 */
    float pitch_gain,  /* i  : adaptive codebook gain                           */
    float tilt_code,   /* i  : tilt factor                                      */
    int mode,        /* i  : innovative codebook mode                         */
    int pre_emphasis,/* i  : use pre_emphasis                                 */
    int pitch_sharpening, /* i  : use pitch sharpening                        */
    int phase_scrambling, /* i  : use phase scrambling                        */
    int formant_enh, /* i  : use formant enhancement                          */
    int formant_tilt,/* i  : use tilt of formant enhancement                  */
    float formant_enh_num, /* i  : formant enhancement numerator weighting factor*/
    float formant_enh_den, /* i  : formant enhancement denominator weighting factor*/
    const short i_subfr,     /* i  : subframe index                                   */
    const float *Aq,         /* i  : quantized LPC coefficients                       */
    float *h1,         /* i  : impulse response of weighted synthesis filter    */
    float *xn,         /* i  : Close-loop Pitch search target vector            */
    float *cn,         /* i  : Innovative codebook search target vector         */
    float *y1,         /* i  : zero-memory filtered adaptive excitation         */
    float *y2,         /* o  : zero-memory filtered algebraic excitation        */
    int acelpautoc,  /* i  : autocorrelation mode enabled                     */
    int **pt_indice, /* i/o: quantization indices pointer                     */
    float *code        /* o  : innovative codebook                              */
    ,const short L_frame,     /* i  : length of the frame                              */
    const short last_L_frame,/* i  : length of the last frame                         */
    const long total_brate   /* i  : total bit-rate                                   */
);

short E_ACELP_code43bit(
    const float code[],
    long unsigned *ps,
    int *p,
    unsigned short idxs[]
);

void fcb_pulse_track_joint(
    unsigned short *idxs,
    int wordcnt,
    unsigned long *index_n,
    int *pulse_num,
    int track_num
);

void D_ACELP_indexing(
    Float32 code[],
    PulseConfig config,
    int num_tracks,
    int prm[],
    short *BER_detect
);

void D_ACELP_decode_43bit(
    unsigned short idxs[],
    float code[],
    int *pulsestrack
);

void fcb_pulse_track_joint_decode(
    unsigned short *idxs,
    int wordcnt,
    long unsigned *index_n,
    int *pulse_num,
    int track_num
);

void lag_wind(
    float r[],            /* i/o: autocorrelations                                       */
    const short m,        /* i  : order of LP filter                                     */
    const int sr,         /* i  : sampling rate                                          */
    const short strength  /* i  : LAGW_WEAK, LAGW_MEDIUM, or LAGW_STRONG                 */
);

void adapt_lag_wind(
    float r[],            /* i/o: autocorrelations                                       */
    int m,                /* i  : order of LP filter                                     */
    const int Top,        /* i  : open loop pitch lags from curr. frame (or NULL if n/a) */
    const float Tnc,      /* i  : open loop pitch gains from curr. frame (NULL if n/a)   */
    int sr                /* i  : sampling rate                                          */
);

void hp20(
    Float32 signal[],
    Word32 lg,
    Float32 mem[],
    Word32 fs
);

void ham_cos_window(
    float *fh,
    int n1,
    int n2
);

float correlation_shift(                /* o  : noise dependent voicing correction      */
    const float totalNoise              /* i  : noise estimate over all critical bands  */
);

void init_coder_ace_plus(
    Encoder_State *st                   /* i/o: encoder state structure             */
);

void core_coder_reconfig(
    Encoder_State *st                   /* i/o: encoder state structure             */
);

void core_coder_mode_switch(
    Encoder_State *st,                    /* i/o: encoder state structure             */
    int bandwidth_in,           /* i  : bandwidth                           */
    int bitrate                 /* i  : bitrate                             */
);

void enc_acelp_tcx_main(
    const float new_samples[],          /* i  : new samples                         */
    Encoder_State *st,                    /* i/o: encoder state structure             */
    const short coder_type,             /* i  : coding type                         */
    const short pitch[3],               /* i  : open-loop pitch values for quantiz. */
    const float voicing[3],             /* i  : open-loop pitch gains               */
    float Aw[NB_SUBFR16k*(M+1)],  /* i  : weighted A(z) unquant. for subframes*/
    const float lsp_new[M],             /* i  : LSPs at the end of the frame        */
    const float lsp_mid[M],             /* i  : LSPs at the middle of the frame     */
    HANDLE_FD_CNG_ENC hFdCngEnc,        /* i/o: CNG handle                      */
    float bwe_exc_extended[],     /* i/o: bandwidth extended excitation       */
    float *voice_factors,         /* o  : voicing factors                     */
    float pitch_buf[],            /* o  : floating pitch for each subframe    */
    short vad_hover_flag
);

void pitch_pred_linear_fit(
    const short nbLostCmpt,             /* i:   bfi counter                                     */
    const short last_good,              /* i:   last classification type                        */
    float *old_pitch_buf,          /* i:   pitch lag buffer                                */
    float *old_fpitch,             /* i/o: pitch used for initial ACB generation           */
    float *T0_out,                 /* o:   estimated close loop pitch                      */
    int pit_min,                /* i:   Minimum pitch lag                               */
    int pit_max,                /* i:   Maximum pitch lag                               */
    float *mem_pitch_gain,         /* i:   lag pitch gain [0] is the most recent subfr lag */
    int limitation,
    short plc_use_future_lag,     /* i: number of subframes to predict                    */
    short *extrapolationFailed   /* o: flag if extrap decides not to change the pitch    */
    ,int nb_subfr                 /* i:   number of ACELP subframes                       */
);

void get_subframe_pitch(
    int nSubframes,             /* i:   number of subframes                             */
    float  pitchStart,             /* i:   starting pitch lag (in subframe -1)             */
    float  pitchEnd,               /* i:   ending pitch lag (in subframe nSubframes-1)     */
    float *pitchBuf                /* o:   interpolated pitch lag per subframe             */
);


void core_encode_openloop(
    Encoder_State *st,                    /* i/o: encoder state structure             */
    const short coder_type,             /* i  : coding type                         */
    const short pitch[3],               /* i  : open-loop pitch values for quantiz. */
    const float voicing[3],             /* i  : open-loop pitch gains               */
    const float Aw[NB_SUBFR16k*(M+1)],  /* i  : weighted A(z) unquant. for subframes*/
    const float lsp_new[M],             /* i  : LSPs at the end of the frame        */
    const float lsp_mid[M],             /* i  : LSPs at the middle of the frame     */
    float *pitch_buf,             /* i/o: floating pitch values for each subfr*/
    float *voice_factors,         /* o  : voicing factors                     */
    float *bwe_exc,               /* o  : excitation for SWB TBE              */
    const short vad_hover_flag
);

void core_acelp_tcx20_switching(
    Encoder_State *st,                    /* i/o: encoder state structure             */
    const short vad_flag,
    short sp_aud_decision0,
    float non_staX,
    short *pitch,                 /* i  : open-loop pitch values for quantiz. */
    float *pitch_fr,              /* i/o: fraction pitch values               */
    float *voicing_fr,            /* i/o: fractional voicing values           */
    const float currTempFlatness,       /* i  : flatness                            */
    const float lsp_mid[M],             /* i  : LSPs at the middle of the frame     */
    const float stab_fac                /* i  : LP filter stability                 */
);

void core_encode_twodiv(
    const float new_samples[],          /* i  : new samples                         */
    Encoder_State *st,                    /* i/o: encoder state structure             */
    const short coder_type,             /* i  : coding type                         */
    const short pitch[3],               /* i  : open-loop pitch values for quantiz. */
    const float voicing[3],             /* i  : open-loop pitch gains               */
    float Aw[NB_SUBFR16k*(M+1)]   /* i  : weighted A(z) unquant. for subframes*/
);

void core_encode_update(
    Encoder_State *st                     /* i/o: encoder state structure             */
);

void core_encode_update_cng(
    Encoder_State *st,                    /* i/o: encoder state structure             */
    float *timeDomainBuffer,
    float *A,
    const float Aw[]                    /* i  : weighted A(z) unquant. for subframes*/
);

void core_signal_analysis_high_bitrate(
    const float *new_samples,
    const short T_op[3],                /* i  : open-loop pitch values for quantiz. */
    const float voicing[3],             /* i  : open-loop pitch gains               */
    const short pitch_12k8[2],          /* i  : open-loop pitch @12.8kHz for adapt. lag windowing */
    float lsp[],
    float lsp_mid[],
    Encoder_State *st,
    int pTnsSize[],
    int pTnsBits[],
    int param_core[],
    int *ltpBits,
    int L_frame
    ,int L_frameTCX
);


float get_gain(     /* output: codebook gain (adaptive or fixed)    */
    float x[],      /* input : target signal                        */
    float y[],      /* input : filtered codebook excitation         */
    int n,          /* input : segment length                       */
    float *en_y     /* output: energy of y (sum of y[]^2, optional) */
);

void encode_acelp_gains(
    float *code,
    int gains_mode,
    float mean_ener_code,
    short clip_gain,
    ACELP_CbkCorr *g_corr,
    float *gain_pit,
    float *gain_code,
    int **pt_indice,
    float *past_gcode,
    float *gain_inov,
    int L_subfr,
    float *code2,
    float *gain_code2,
    short noisy_speech_flag
);

int Mode2_gain_enc_mless(
    const float *code,      /* i  : algebraic excitation                                            */
    int lcode,              /* i  : Subframe size                                                   */
    float *gain_pit,        /* o  : quantized pitch gain                                            */
    float *gain_code,       /* o  : quantized codebook gain                                         */
    ACELP_CbkCorr *coeff,   /* i/o: correlations <y1,y1>, -2<xn,y1>,<y2,y2>, -2<xn,y2> and 2<y1,y2> */
    float mean_ener,        /* i  : mean_ener defined in open-loop (3 bits)                         */
    const short clip_gain,  /* i  : gain pitch clipping flag (1 = clipping)                         */
    float *past_gcode,      /* i/o: past gain of code                                               */
    float *gain_inov,       /* o  : unscaled innovation gain                                        */
    const short coder_type  /* i  : type of coder                                                   */
);

void decode_acelp_gains(
    float *code,
    int gains_mode,
    float mean_ener_code,
    float *gain_pit,
    float *gain_code,
    int **pt_indice,
    float *past_gpit,
    float *past_gcode,
    float *gain_inov,
    int L_subfr,
    float *code2,
    float *gain_code2
);

void Es_pred_dec(
    float *Es_pred,       /* o  : predicited scaled innovation energy     */
    const int   enr_idx,        /* i  : indice                                  */
    const short nb_bits,        /* i  : number of bits                          */
    const short no_ltp          /* i  : no LTP flag                             */
);

void Mode2_pit_encode(
    short coder_type,   /* i  : coding model                               */
    short i_subfr,      /* i  : subframe index                             */
    int **pt_indice,  /* i/o: quantization indices pointer               */
    float *exc,         /* i/o: pointer to excitation signal frame         */
    const short *T_op,        /* i  : open loop pitch estimates in current frame */
    int *T0_min,      /* i/o: lower limit for close-loop search          */
    int *T0_min_frac, /* i/o: lower limit for close-loop search          */
    int *T0_max,      /* i/o: higher limit for close-loop search         */
    int *T0_max_frac, /* i/o: higher limit for close-loop search         */
    int *T0,          /* i/o: close loop integer pitch                   */
    int *T0_frac,     /* i/o: close loop fractional part of the pitch    */
    int *T0_res,      /* i/o: close loop pitch resolution                */
    float *h1,          /* i  : weighted filter impulse response           */
    float *xn,          /* i  : target vector                              */
    int  pit_min,
    int  pit_fr1,
    int  pit_fr1b,
    int  pit_fr2,
    int  pit_max,
    int  pit_res_max
);

void limit_T0_voiced(
    int nbits,
    int res,
    int T0,           /* i  : rough pitch estimate around which the search is done */
    int T0_frac,      /* i  : pitch estimate fractional part                       */
    int T0_res,       /* i  : pitch resolution                                     */
    int *T0_min,      /* o  : lower pitch limit                                    */
    int *T0_min_frac, /* o  : lower pitch limit                                    */
    int *T0_max,      /* o  : higher pitch limit                                   */
    int *T0_max_frac, /* o  : higher pitch limit                                   */
    int pit_min,      /* i  : Minimum pitch lag                                    */
    int pit_max       /* i  : Maximum pitch lag                                    */
);

void Mode2_abs_pit_enc(
    short T0,           /* i  : integer pitch lag              */
    int T0_frac,      /* i  : pitch fraction                 */
    int **pt_indice,  /* i/o: pointer to Vector of Q indexes */
    short pit_min,
    short pit_fr1,
    short pit_fr2,
    short pit_res_max
);

void Mode2_delta_pit_enc(
    short T0,               /* i  : integer pitch lag              */
    int T0_frac,          /* i  : pitch fraction                 */
    int T0_res,           /* i  : pitch resolution               */
    short T0_min,           /* i/o: delta search min               */
    short T0_min_frac,      /* i/o: delta search min               */
    int **pt_indice         /* i/o: pointer to Vector of Q indexes */
);

float Mode2_pit_decode(        /* o:   floating pitch value                    */
    const short coder_type, /* i:   coding model                            */
    short i_subfr,          /* i:   subframe index                          */
    int L_subfr,
    int **pt_indice,        /* i/o: quantization indices pointer            */
    int *T0,                /* o:   close loop integer pitch                */
    int *T0_frac,           /* o:   close loop fractional part of the pitch */
    int *T0_res,            /* i/o: pitch resolution                        */
    int *T0_min,            /* i/o: lower limit for close-loop search       */
    int *T0_min_frac,       /* i/o: lower limit for close-loop search       */
    int *T0_max,            /* i/o: higher limit for close-loop search      */
    int *T0_max_frac,       /* i/o: higher limit for close-loop search      */
    int pit_min,
    int pit_fr1,
    int pit_fr1b,
    int pit_fr2,
    int pit_max,
    int pit_res_max
);

void Mode2_abs_pit_dec(
    int *T0,                  /* o:   integer pitch lag              */
    int *T0_frac,             /* o:   pitch fraction                 */
    int *T0_res,              /* o:   pitch resolution               */
    int **pt_indice,          /* i/o: pointer to Vector of Q indexes */
    int pit_min,
    int pit_fr1,
    int pit_fr2,
    int pit_res_max
);

void Mode2_delta_pit_dec(
    int       *T0,          /* o:   integer pitch lag              */
    int       *T0_frac,     /* o:   pitch fraction                 */
    int       T0_res,       /* i:   pitch resolution               */
    int       *T0_min,      /* i: delta search min                 */
    int       *T0_min_frac, /* i: delta search min                 */
    int       **pt_indice   /* i/o: pointer to Vector of Q indexes */
);

void formant_post_filt(
    PFSTAT *pfstat,        /* i/o: Post filter related memories    */
    float *synth_in,      /* i  : 12k8 synthesis                  */
    const float *Aq,            /* i  : LP filter coefficient           */
    float *synth_out,     /* i/o: input signal                    */
    const short L_frame,        /* i  : frame length                    */
    const short L_subfr,        /* i  : sub-frame length                */
    const float lp_noise,       /* i  : background noise energy         */
    const long  rate,           /* i  : bit-rate                        */
    const short off_flag        /* i  : Off flag                        */
);

void qlpc_avq(
    const float *lsp,            /* (i) Input LSF vectors              */
    const float *lspmid,
    float *lsf_q,                /* (o) Quantized LFS vectors          */
    float *lsfmid_q,
    int *index,                  /* (o) Quantization indices           */
    int *nb_indices,             /* (o) Number of quantization indices */
    int *nbbits,                 /* (o) Number of quantization bits    */
    int core,
    float sr_core
);

int encode_lpc_avq(
    Encoder_State *st,
    int numlpc,
    int *param_lpc,
    int mode
);

int dlpc_avq(
    int *index,          /* (i)   Quantization indices                       */
    float *LSF_Q,        /* (o)   Quantized LSF vectors                      */
    int numlpc,          /* (i) Number of sets of lpc */
    float sr_core
);

int decode_lpc_avq( Decoder_State *st, int numlpc, int *param_lpc );

int vlpc_1st_cod(       /* output: codebook index                  */
    const float *lsf,     /* input:  vector to quantize              */
    float *lsfq,          /* i/o:    i:prediction   o:quantized lsf  */
    float sr_core
    ,float *w              /* o: lsf weights */
);

int vlpc_2st_cod(       /* output: number of allocated bits        */
    const float *lsf,     /* input:  normalized vector to quantize   */
    float *lsfq,          /* i/o:    i:1st stage   o:1st+2nd stage   */
    int *indx,            /* output: index[] (4 bits per words)      */
    int mode,             /* input:  0=abs, >0=rel                   */
    float sr_core
);

void vlpc_2st_dec(
    float *lsfq,          /* i/o:    i:1st stage   o:1st+2nd stage   */
    int *indx,            /* input:  index[] (4 bits per words)      */
    int mode,             /* input:  0=abs, >0=rel                   */
    float sr_core
);

void lsf_weight_2st(
    const float *lsfq,
    float *w,
    int mode,
    float sr_core
);

void mdct_window_sine(
    float *window,
    int n
);

void mdct_window_aldo(
    float *window1,
    float *window2,
    int n
);

void AVQ_cod_lpc(
    const float nvec[],     /* i:   vector to quantize              */
    int   nvecq[],    /* o:   quantized normalized vector (assuming the bit budget is enough) */
    int   *indx,      /* o:   index[] (4 bits per words)      */
    const short Nsv         /* i:   number of subvectors (lg=Nsv*8) */
);

void AVQ_dec_lpc(
    const int   indx[],     /* i  : index[] (4 bits per words)      */
    int   nvecq[],    /* o  : vector quantized                */
    const short Nsv         /* i  : number of subvectors (lg=Nsv*8) */
);

void vlpc_1st_dec(
    int index,                /* input:  codebook index                  */
    float *lsfq,              /* i/o:    i:prediction   o:quantized lsf  */
    float sr_core
);

void WindowSignal(
    TCX_config const *tcx_cfg,                /* input: configuration of TCX              */
    int offset,                               /* input: left folding point offset relative to the input signal pointer */
    const short left_overlap_mode,            /* input: overlap mode of left window half  */
    const short right_overlap_mode,           /* input: overlap mode of right window half */
    int * left_overlap_length,                /* output: TCX window left overlap length   */
    int * right_overlap_length,               /* output: TCX window right overlap length  */
    float const in[],                         /* input: input signal                      */
    int * L_frame,                            /* input/output: frame length               */
    float out[],                              /* output: output windowed signal           */
    int fullband                            /* input: fullband flag                     */
);

void HBAutocorrelation(
    TCX_config *tcx_cfg, /* input: configuration of TCX               */
    int left_mode,       /* input: overlap mode of left window half   */
    int right_mode,      /* input: overlap mode of right window half  */
    float speech[],      /* input: speech                             */
    int L_frame_glob,    /* input: frame length                       */
    float *r,            /* output: autocorrelations vector           */
    int m                /* input : order of LP filter                */
);

void TNSAnalysis(
    TCX_config *tcx_cfg, /* input: configuration of TCX */
    int L_frame,         /* input: frame length */
    int L_spec,
    const short tcxMode, /* input: TCX mode for the frame/subframe - TCX20 | TCX10 | TCX 5 (meaning 2 x TCX 5) */
    int isAfterACELP,    /* input: Flag indicating if the last frame was ACELP. For the second TCX subframe it should be 0  */
    float spectrum[],    /* input: MDCT spectrum of the subframe */
    STnsData * pTnsData, /* output: Tns data */
    int * pfUseTns,      /* output: Flag indicating if TNS is used */
    float* predictionGain
);

void ShapeSpectrum(
    TCX_config *tcx_cfg,/*input: configuration of TCX*/
    float A[],          /* input: quantized coefficients NxAz_q[M+1] */
    float gainlpc[],    /* output: MDCT gains for the previous frame */
    int L_frame_glob,   /* input: frame length             */
    int L_spec,
    float spectrum[],   /* i/o: MDCT spectrum */
    int fUseTns,        /* output: Flag indicating if TNS is used */
    Encoder_State *st
);

void QuantizeSpectrum(
    TCX_config *tcx_cfg,/*input: configuration of TCX*/
    float A[],          /* input: quantized coefficients NxAz_q[M+1] */
    Word16 Aqind[],     /* input: frame-independent quantized coefficients (M+1) */
    float gainlpc[],    /* input: MDCT gains of the previous frame */
    float synth[],
    int L_frame_glob,   /* input: frame length             */
    int L_frameTCX_glob,
    int L_spec,
    int nb_bits,        /*input: bit budget*/
    int tcxonly,        /*input: only TCX flag*/
    float spectrum[],   /* i/o: MDCT spectrum, input is shaped MDCT spectrum */
    STnsData * pTnsData,/* input: Tns data */
    int fUseTns,        /* input: Flag indicating if TNS is used */
    int tnsSize,        /* input: number of tns parameters put into prm */
    LPD_state *LPDmem,  /*i/o: memories*/
    int prm[],          /* output: tcx parameters          */
    int frame_cnt,       /* input: frame counter in the super_frame */
    Encoder_State *st,
    CONTEXT_HM_CONFIG *hm_cfg
);

/* Returns: index of next coefficient */
int get_next_coeff_mapped(
    int ii[2],                /* i/o: coefficient indexes       */
    int *pp,                  /* o  : peak(1)/hole(0) indicator */
    int *idx,                 /* o  : index in unmapped domain  */
    CONTEXT_HM_CONFIG *hm_cfg /* i  : HM configuration          */
);

/* Returns: index of next coefficient */
int get_next_coeff_unmapped(
    int *ii,                  /* i/o: coefficient index         */
    int *idx                  /* o  : index in unmapped domain  */
);

int update_mixed_context(
    int ctx,
    int a
);

void ACcontextMapping_encode2_no_mem_s17_LC(
    Encoder_State *st,
    int *x,
    int nt,
    int lastnz,
    int nbbits,
    int resQMaxBits,
    CONTEXT_HM_CONFIG *hm_cfg
);

long ACcontextMapping_decode2_no_mem_s17_LC(
    Decoder_State *st,          /* i/o: decoder state                           */
    int *x,                     /* o: decoded spectrum                          */
    long nt,                    /* i: size of spectrum                          */
    int nbbits,                 /* i: bit budget                                */
    int resQMaxBits,            /* i: residual coding maximum bits              */
    CONTEXT_HM_CONFIG *hm_cfg   /* i: context-based harmonic model configuration*/
);

int ACcontextMapping_encode2_estimate_no_mem_s17_LC(
    const int *x,
    int nt,
    int *lastnz,
    int *nEncoded,
    int target,
    int *stop,
    CONTEXT_HM_CONFIG *hm_cfg
);

void tcx_windowing_analysis(
    float const *signal,        /* i: signal vector                              */
    int L_frame,                /* i: frame length                               */
    int left_overlap,           /* i: left overlap length                        */
    float const *left_win,      /* i: left overlap window                        */
    int right_overlap,          /* i: right overlap length                       */
    float const *right_win,     /* i: right overlap window                       */
    float *output               /* o: windowed signal vector                     */
);

void tcx_windowing_synthesis_current_frame(
    float *signal,              /* i/o: signal vector                            */
    float *window,              /* i: TCX window vector                          */
    float *window_half,         /* i: TCX window vector for half-overlap window  */
    float *window_min,          /* i: TCX minimum overlap window                 */
    int window_length,          /* i: TCX window length                          */
    int window_half_length,     /* i: TCX half window length                     */
    int window_min_length,      /* i: TCX minimum overlap length                 */
    int left_rect,              /* i: left part is rectangular                   */
    int left_mode,              /* i: overlap mode of left window half           */
    float *acelp_zir,           /* i: acelp ZIR                                  */
    float *old_syn,             /* i: old synthesis                              */
    float *syn_overl,           /* i: overlap synthesis                          */
    float *A_zir,
    float *window_trans,        /* i: window for transition from ACELP           */
    int acelp_zir_len,
    int acelp_mem_len,
    int last_core_bfi,          /* i : last mode                                  */
    int last_is_cng,
    int fullbandScale
);

void tcx_windowing_synthesis_past_frame(
    float *signal,              /* i/o: signal vector                            */
    float *window,              /* i: TCX window vector                          */
    float *window_half,         /* i: TCX window vector for half-overlap window  */
    float *window_min,          /* i: TCX minimum overlap window                 */
    int window_length,          /* i: TCX window length                          */
    int window_half_length,     /* i: TCX half window length                     */
    int window_min_length,      /* i: TCX minimum overlap length                 */
    int right_mode              /* i: overlap mode (left_mode of current frame)  */
);

/* tcx_utils_enc.c */
void ProcessIGF(
    IGF_ENC_INSTANCE_HANDLE         const hInstance,          /**< in: instance handle of IGF Encoder */
    Encoder_State                        *st,                 /**< in: Encoder state */
    float                                *pMDCTSpectrum,      /**< in: MDCT spectrum */
    float                                *pPowerSpectrum,     /**< in: MDCT^2 + MDST^2 spectrum, or estimate */
    int                                   isTCX20,            /**< in: flag indicating if the input is TCX20 or TCX10/2xTCX5 */
    int                                   isTNSActive,        /**< in: flag indicating if the TNS is active */
    int                                   isTransition,       /**< in: flag indicating if the input is the transition from from ACELP to TCX20/TCX10 */
    int                                   frameno             /**< in: flag indicating index of current subframe */
);

/* tcx_utils_enc.c */
void AnalyzePowerSpectrum(
    Encoder_State *st,              /* i/o: encoder states                                  */
    int L_frame,                    /* input: frame length                                  */
    int L_frameTCX,                 /* input: full band frame length                        */
    int left_overlap,               /* input: left overlap length                           */
    int right_overlap,              /* input: right overlap length                          */
    float const mdctSpectrum[],     /* input: MDCT spectrum                                 */
    float const signal[],           /* input: windowed signal corresponding to mdctSpectrum */
    float powerSpec[]               /* output: Power spectrum                               */
);

void lpc2mdct(
    float *lpcCoeffs,
    int lpcOrder,
    float *mdct_gains
);

void mdct_preShaping(
    float x[],
    int lg,
    const float gains[]
);

void mdct_noiseShaping(
    float x[],
    int lg,
    const float gains[]
);

void mdct_noiseShaping_interp(
    float x[],
    int lg,
    const float gains[]
);

void AdaptLowFreqEmph(
    float x[],
    int xq[],
    float invGain,
    short tcx_lpc_shaped_ari,
    const float lpcGains[],
    const int lg
);

void PsychAdaptLowFreqEmph(
    float x[],
    const float lpcGains[]
);

void PsychAdaptLowFreqDeemph(
    float x[],
    const float lpcGains[],
    float lf_deemph_factors[]
);

void AdaptLowFreqDeemph(
    float x[],
    short tcx_lpc_shaped_ari,
    const float lpcGains[],
    const int lg,
    float lf_deemph_factors[]
);

float SQ_gain(          /* output: SQ gain                   */
    float x[],          /* input:  vector to quantize        */
    int nbitsSQ,      /* input:  number of bits targeted   */
    int lg            /* input:  vector size (2048 max)    */
);

void tcx_scalar_quantization(
    float *x,                 /* i: input coefficients            */
    int *xq,                  /* o: quantized coefficients        */
    int L_frame,              /* i: frame length                  */
    float gain,               /* i: quantization gain             */
    float offset,             /* i: rounding offset (deadzone)    */
    int *memQuantZeros,       /* o: coefficients set to 0         */
    const int alfe_flag
);

int tcx_scalar_quantization_rateloop(
    float *x,                 /* i  : input coefficients            */
    int *xq,                  /* o  : quantized coefficients        */
    int L_frame,              /* i  : frame length                  */
    float *gain,              /* i/o: quantization gain             */
    float offset,             /* i  : rounding offset (deadzone)    */
    int *memQuantZeros,       /* o  : coefficients set to 0         */
    int *lastnz_out,          /* i/o: last nonzero coeff index      */
    int target,               /* i  : target number of bits         */
    int *nEncoded,            /* o  : number of encoded coeff       */
    int *stop,                /* i/o: stop param                    */
    int sqBits_in_noStop,     /* i  : number of sqBits as determined in prev. quant. stage, w/o using stop mechanism (ie might exceed target bits) */
    int sqBits_in,            /* i  : number of sqBits as determined in prev. quant. stage, using stop mechanism (ie always <= target bits) */
    int tcxRateLoopOpt,       /* i  : turns on/off rateloop optimization */
    const int tcxonly,
    CONTEXT_HM_CONFIG *hm_cfg
);

void QuantizeGain(
    int n,
    float * pGain,
    int * pQuantizedGain
);

void tcx_noise_factor(
    float *x_orig,          /* i: unquantized mdct coefficients             */
    float *sqQ,             /* i: quantized mdct coefficients               */
    int iFirstLine,         /* i: first coefficient to be considered        */
    int lowpassLine,        /* i: last nonzero coefficients after low-pass  */
    int nMinHoleSize,       /* i: minimum size of hole to be checked        */
    int L_frame,            /* i: frame length                              */
    float gain_tcx,         /* i: tcx gain                                  */
    float tiltCompFactor,   /* i: LPC tilt compensation factor              */
    float *fac_ns,          /* o: noise factor                              */
    int *quantized_fac_ns   /* o: quantized noise factor                    */
);

void tcx_noise_filling(
    float *Q,
    const int noiseFillSeed,
    const int iFirstLine,
    const int lowpassLine,
    const int nTransWidth,
    const int L_frame,
    float tiltCompFactor,
    float fac_ns,
    unsigned char *infoTCXNoise
);

void tcx_encoder_memory_update(
    const float *wsig,      /* i : target weighted signal                     */
    float *xn_buf,          /* i/o: mdct output buffer/TD weigthed synthesis  */
    int L_frame_glob,       /* i: global frame length                         */
    const float *Ai,        /* i: Unquantized (interpolated) LPC coefficients */
    float *A,               /* i: Quantized LPC coefficients                  */
    float preemph,          /* i: preemphasis factor*/
    LPD_state *LPDmem,      /* i/o: coder memory state                        */
    Encoder_State *st,
    int m,
    float *synth
);

void tcx_decoder_memory_update(
    float *xn_buf,          /* i: mdct output buffer                          */
    float *synth,           /* i/o: synth                                     */
    int L_frame_glob,       /* i: global frame length                         */
    float *A,               /* i: Quantized LPC coefficients                  */
    Decoder_State *st,      /* i/o : decoder memory state                     */
    float *syn              /* o: st->syn                                     */
);



int tcx_ari_res_Q_spec( /* Returns: number of bits used (including "bits") */
    const float x_orig[], /* i: original spectrum                  */
    const int signs[],    /* i: signs (x_orig[.]<0)                */
    float x_Q[],          /* i/o: quantized spectrum               */
    int L_frame,          /* i: number of lines                    */
    float gain,           /* i: TCX gain                           */
    int prm[],            /* o: bit-stream                         */
    int target_bits,      /* i: number of bits available           */
    int bits,             /* i: number of bits used so far         */
    float deadzone,       /* i: quantizer deadzone                 */
    const float x_fac[]   /* i: spectrum post-quantization factors */
);

int tcx_ari_res_invQ_spec(/* Returns: number of bits used (including "bits") */
    float x_Q[],          /* i/o: quantized spectrum               */
    int L_frame,          /* i: number of lines                    */
    const int prm[],      /* i: bit-stream                         */
    int target_bits,      /* i: number of bits available           */
    int bits,             /* i: number of bits used so far         */
    float deadzone,       /* i: quantizer deadzone                 */
    const float x_fac[]   /* i: spectrum post-quantization factors */
);

int tcx_res_Q_gain(
    float sqGain,
    float *gain_tcx,
    int *prm,
    int sqTargetBits
);

int tcx_res_Q_spec(
    float *x_orig,
    float *x_Q,
    int L_frame,
    float sqGain,
    int *prm,
    int sqTargetBits,
    int bits,
    float sq_round,
    const float lf_deemph_factors[]
);

int tcx_res_invQ_gain(
    float *gain_tcx,
    const int *prm,
    int resQBits
);

int tcx_res_invQ_spec(
    float *x,
    int L_frame,
    const int *prm,
    int resQBits,
    int bits,
    float sq_round,
    const float lf_deemph_factors[]
);

void InitTnsConfigs(
    int nSampleRate,
    int L_frame,
    STnsConfig tnsConfig[2][2],
    int igfStopFreq,
    int bitrate
);

void SetTnsConfig(
    TCX_config * tcx_cfg,
    int isTCX20,
    int isAfterACELP
);

void ari_copy_states(
    Tastat *source,
    Tastat *dest
);

long mul_sbc_14bits(
    long r,
    long c
);

void ari_start_encoding_14bits(
    Tastat *s
);

long ari_encode_14bits(
    short *ptr,
    long bp,
    Tastat *s,
    long symbol,
    unsigned short *cum_freq
);

long ari_encode_14bits_ext(
    int *ptr,
    long bp,
    Tastat *s,
    long symbol,
    const unsigned short *cum_freq
);

long ari_done_encoding_14bits(
    int *ptr,
    long bp,
    Tastat *s
);

short ari_encode_check_budget_14bits(
    long bp,
    Tastat *s,
    long nb_bits
);

void ari_start_decoding_14bits(
    Decoder_State *st,
    Tastat *s
);

long ari_start_decoding_14bits_prm(
    const int *ptr,
    long bp,
    Tastat *s
);

void ari_decode_14bits_s17_ext(
    Decoder_State *st,
    int *res,
    Tastat *s,
    const unsigned short *cum_freq
);

void ari_decode_14bits_s27_ext(
    Decoder_State *st,
    int *res,
    Tastat *s,
    const unsigned short *cum_freq
);

void ari_decode_14bits_bit_ext(
    Decoder_State *st,
    int *res,
    Tastat *s
);

Word16 expfp(         /* o: Q15 */
    Word16 x,           /* i: mantissa  Q15-e */
    Word16 x_e          /* i: exponent  Q0 */
);

void powfp_odd2(
    Word16 base,     /* Q15 */
    Word16 exp,      /* Q0  */
    Word16 *pout1,   /* Q15 */
    Word16 *pout2    /* Q15 */
);

void tcx_arith_scale_envelope(
    Word16 L_spec_core,         /* i: number of lines to scale    Q0 */
    Word16 L_frame,             /* i: number of lines             Q0 */
    Word32 env[],               /* i: unscaled envelope           Q16 */
    Word16 target_bits,         /* i: number of available bits    Q0 */
    Word16 low_complexity,      /* i: low-complexity flag         Q0 */
    Word16 s_env[],             /* o: scaled envelope             Q15-e */
    Word16 *s_env_e             /* o: scaled envelope exponent    Q0 */
);

void tcx_arith_render_envelope(
    const Word16 A_ind[],         /* i: LPC coefficients of signal envelope        */
    Word16 L_frame,               /* i: number of spectral lines                   */
    Word16 L_spec,
    Word16 preemph_fac,           /* i: pre-emphasis factor                        */
    Word16 gamma_w,               /* i: A_ind -> weighted envelope factor          */
    Word16 gamma_uw,              /* i: A_ind -> non-weighted envelope factor      */
    Word32 env[]                  /* o: shaped signal envelope                     */
);

long ari_encode_14bits_range(
    int *ptr,
    long bp,
    long bits,
    Tastat *s,
    unsigned short cum_freq_low,
    unsigned short cum_freq_high
);

long ari_encode_14bits_sign(
    int *ptr,
    long bp,
    long bits,
    Tastat *s,
    long sign
);

long ari_done_cbr_encoding_14bits(
    int *ptr,
    long bp,
    long bits,
    Tastat *s
);

long ari_decode_14bits_pow(
    const int *ptr,
    long bp,
    long bits,
    int *res,
    Tastat *s,
    unsigned base
);

long ari_decode_14bits_sign(
    const int *ptr,
    long bp,
    long bits,
    int *res,
    Tastat *s
);

void tcx_arith_encode_envelope(
    float spectrum[],                       /* i/o: MDCT coefficients           */
    int signs[],                            /* o: signs (spectrum[.]<0)         */
    int L_frame,                            /* i: frame or MDCT length          */
    int L_frame_orig,                       /* i: length w/o BW limitation      */
    Encoder_State *st,                      /* i/o: coder state                 */
    const Word16 A_ind[],                   /* i: quantised LPC coefficients    */
    int target_bits,                        /* i: number of available bits      */
    int prm[],                              /* o: bitstream parameters          */
    int use_hm,                             /* i: use HM in current frame?      */
    int prm_hm[],                           /* o: HM parameter area             */
    short tcxltp_pitch,                     /* i: TCX LTP pitch in FD, -1 if n/a*/
    int *arith_bits,                        /* o: bits used for ari. coding     */
    int *signaling_bits                     /* o: bits used for signaling       */
    ,int low_complexity                     /* i: low-complexity flag           */
);

void tcx_arith_decode_envelope(
    float q_spectrum[],                     /* o: quantised MDCT coefficients   */
    int L_frame,                            /* i: frame or MDCT length          */
    int L_frame_orig,                       /* i: length w/o BW limitation      */
    Decoder_State *st,                      /* i/o: coder state                 */
    const short coder_type,                 /* i  : coder type                  */
    const Word16 A_ind[],                   /* i: quantised LPC coefficients    */
    float tcxltp_gain,                      /* i: TCX LTP gain                  */
    int target_bits,                        /* i: number of available bits      */
    const int prm[],                        /* i: bitstream parameters          */
    int use_hm,                             /* i: use HM in current frame?      */
    const int prm_hm[],                     /* i: HM parameter area             */
    short tcxltp_pitch,                     /* i: TCX LTP pitch in FD, -1 if n/a*/
    int *arith_bits,                        /* o: bits used for ari. coding     */
    int *signaling_bits                     /* o: bits used for signaling       */
    ,int low_complexity                     /* i: low-complexity flag           */
);

void UnmapIndex(
    int PeriodicityIndex,
    int Bandwidth,
    short LtpPitchLag,
    int SmallerLags,
    int *FractionalResolution,
    int *Lag
);

int SearchPeriodicityIndex(            /* Returns: PeriodicityIndex */
    const float Mdct[],                /* (I) Coefficients, Mdct[0..NumCoeffs-1]                      */
    const float UnfilteredMdct[],      /* (I) Unfiltered coefficients, UnfilteredMdct[0..NumCoeffs-1] */
    int NumCoeffs,                     /* (I) Number of coefficients                                  */
    int TargetBits,                    /* (I) Target bit budget (excl. Done flag)                     */
    short LtpPitchLag,
    float LtpGain,                     /* (I) LTP gain                                                */
    float *RelativeScore               /* (O) Energy concentration factor                             */
);

void ConfigureContextHm(
    int NumCoeffs,                      /* (I) Number of coefficients                         */
    int TargetBits,                     /* (I) Target bit budget (excl. Done flag)            */
    int PeriodicityIndex,               /* (I) Pitch related index                            */
    short LtpPitchLag,                  /* (I) TCX-LTP pitch in F.D.                          */
    CONTEXT_HM_CONFIG *hm_cfg           /* (O) Context-based harmonic model configuration     */
);

int EncodeIndex(
    int Bandwidth,                     /* 0: NB, 1: (S)WB */
    int PeriodicityIndex,
    Encoder_State *st
);

int CountIndexBits(
    int Bandwidth,                     /* 0: NB, 1: (S)WB */
    int PeriodicityIndex
);

int DecodeIndex(
    Decoder_State *st,
    int Bandwidth,                     /* 0: NB, 1: (S)WB */
    int *PeriodicityIndex
);


#define GET_ADJ(T,L)    GET_ADJ2(T,L,*FractionalResolution)
#define GET_ADJ2(T,L,F) (((L) << (F)) - (T))

int tcx_hm_render(
    int lag,              /* i: pitch lag                             */
    int fract_res,        /* i: fractional resolution of the lag      */
    float LtpGain,        /* i: LTP gain                              */
    Word16 p[]            /* o: harmonic model (Q13)                  */
);

void tcx_hm_modify_envelope(
    Word16 gain,          /* i: HM gain (Q11)                         */
    int lag,
    int fract_res,
    Word16 p[],           /* i: harmonic model (Q13)                  */
    Word32 env[],         /* i/o: envelope (Q16)                      */
    int L_frame           /* i: number of spectral lines              */
);

void tcx_hm_analyse(
    const float abs_spectrum[], /* i: absolute spectrum               */
    int L_frame,          /* i: number of spectral lines              */
    Word32 env[],         /* i/o: envelope shape (Q16)                */
    int targetBits,       /* i: target bit budget                     */
    int coder_type,       /* i: GC/VC mode                            */
    int prm_hm[],         /* o: HM parameters                         */
    short LtpPitchLag,    /* i: LTP pitch lag or -1 if none           */
    float LtpGain,        /* i: LTP gain                              */
    int *hm_bits          /* o: bit consumption                       */
);

void tcx_hm_decode(
    int L_frame,          /* i: number of spectral lines              */
    Word32 env[],         /* i/o: envelope shape (Q16)                */
    int targetBits,       /* i: target bit budget                     */
    int coder_type,       /* i: GC/VC mode                            */
    const int prm_hm[],   /* i: HM parameters                         */
    short LtpPitchLag,    /* i: LTP pitch lag or -1 if none           */
    float LtpGain ,       /* i: LTP gain                              */
    int *hm_bits          /* o: bit consumption                       */
);

void coder_tcx(
    int n,
    TCX_config *tcx_cfg,  /* input: configuration of TCX               */
    float A[],            /* input: quantized coefficients NxAz_q[M+1] */
    Word16 Aqind[],       /* input: frame-independent quantized coefficients (M+1) */
    float synth[],
    int L_frame_glob,     /* input: frame length                       */
    int L_frameTCX_glob,
    int L_spec,
    int nb_bits,          /* input: bit budget                         */
    int tcxonly,          /* input: only TCX flag                      */
    float spectrum[],     /* i/o: MDCT spectrum                        */
    LPD_state *LPDmem,    /* i/o: memories                             */
    int prm[],            /* output: tcx parameters                    */
    Encoder_State *st,
    CONTEXT_HM_CONFIG *hm_cfg
);

void coder_tcx_post(
    Encoder_State *st,
    LPD_state *LPDmem,
    TCX_config *tcx_cfg,
    float *synth,
    float *A,
    const float *Ai,
    float *wsig
);

void decoder_tcx(
    TCX_config *tcx_cfg,  /* input: configuration of TCX              */
    int prm[],            /* input:  parameters                       */
    float A[],            /* input:  coefficients NxAz[M+1]           */
    Word16 Aind[],        /* input: frame-independent coefficients Az[M+1] */
    int L_frame,          /* input:  frame length                     */
    int L_frameTCX,
    int L_spec,
    float synth[],        /* in/out: synth[-M..lg]                    */
    float synthFB[],
    Decoder_State *st,    /* i/o : coder memory state                 */
    const short coder_type,
    int bfi,              /* input:  Bad frame indicator              */
    int frame_cnt,        /* input: frame counter in the super_frame  */
    float stab_fac        /* input: stability of isf                  */
);

void decoder_tcx_post(
    Decoder_State *st,
    float *synth,
    float *synthFB,
    float *A,
    int bfi
);

void coder_acelp(
    ACELP_config *acelp_cfg,      /* i/o: configuration of the ACELP  */
    const short coder_type,     /* i  : coding type                 */
    const float A[],            /* i  : coefficients 4xAz[M+1]      */
    const float Aq[],           /* i  : coefficients 4xAz_q[M+1]    */
    const float speech[],       /* i  : speech[-M..lg]              */
    float synth[],        /* o  : synthesis                   */
    LPD_state *LPDmem,         /* i/o: ACELP memories              */
    const float voicing[],      /* i  : open-loop LTP gain          */
    const short T_op[],         /* i  : open-loop LTP lag           */
    int *prm,           /* o  : acelp parameters            */
    const float stab_fac,
    Encoder_State *st,            /* i/o : coder memory state         */
    HANDLE_PLC_ENC_EVS hPlc_Ext,
    const short target_bits,
    float *gain_pitch_buf,/* o  : gain pitch values           */
    float *gain_code_buf, /* o  : gain code values            */
    float *pitch_buf,     /* o  : pitch values for each subfr.*/
    float *voice_factors, /* o  : voicing factors             */
    float *bwe_exc        /* o  : excitation for SWB TBE      */
);

void coder_acelp_rf(
    const short target_bits,    /* i:   target bits                 */
    const float speech[],       /* i  : speech[-M..lg]              */
    const short coder_type,     /* i  : coding type                 */
    const short rf_frame_type,  /* i  : rf_frame_type               */
    const float A[],            /* i  : coefficients 4xAz[M+1]      */
    const float Aq[],           /* i  : coefficients 4xAz_q[M+1]    */
    const float voicing[],      /* i  : open-loop LTP gain          */
    const short T_op[],         /* i  : open-loop LTP lag           */
    const float stab_fac,       /* i  : LP stability factor         */
    Encoder_State *st,            /* i/o: coder memory state          */
    ACELP_config *acelp_cfg,     /* i/o: configuration of the ACELP  */
    float *exc_rf,        /* i/o: pointer to RF excitation    */
    float *syn_rf         /* i/o: pointer to RF synthesis     */
);

void decoder_acelp(
    Decoder_State *st,            /* i/o:  coder memory state         */
    const short coder_type,     /* i  : coder type                  */
    int prm[],            /* i  : parameters                  */
    const float A[],            /* i  : coefficients NxAz[M+1]      */
    ACELP_config acelp_cfg,      /* i  : ACELP config                */
    float synth[],        /* i/o:   synthesis                 */
    int   *pT,            /* o  :   pitch for all subframe    */
    float *pgainT,        /* o  :   pitch gain for all subfr  */
    const float stab_fac        /* i  : stability of isf            */
    ,float *pitch_buffer   /* o  : pitch values for each subfr.*/
    ,float *voice_factors  /* o  : voicing factors             */
    ,const short LSF_Q_prediction/* i  : LSF prediction mode         */
    ,float *bwe_exc        /* o  : excitation for SWB TBE      */
);

void enc_prm(
    const short coder_type,       /* i  : coding type                 */
    int param[],            /* i  : parameters                  */
    int param_lpc[],        /* i  : LPC parameters              */
    Encoder_State *st,            /* i/o: quantization Analysis values*/
    const short L_Frame,          /* i  : frame length                */
    CONTEXT_HM_CONFIG hm_cfg[],
    short * bits_param_lpc,
    short no_param_lpc
);

void enc_prm_rf(
    Encoder_State *st,
    const short rf_frame_type,
    const short fec_offset
);

void dec_prm(
    short *core,                  /* (0) : current frame mode            */
    short *last_core,             /* (0) : last frame mode               */
    short *coder_type,            /* o  : coder type                     */
    int param[],                  /* (o) : decoded parameters            */
    int param_lpc[],              /* (i) : LPC parameters                */
    Word16 *total_nbbits,         /* i/o : number of bits / decoded bits */
    Decoder_State *st,            /* i/o:  decoder memory state          */
    const int L_frame,
    short * bitsRead
);

void dec_prm_core(
    Decoder_State *st
);

void gauss_L2(
    const float h[],         /* i  : weighted LP filter impulse response     */
    float code[],            /* o  : gaussian excitation                     */
    float y2[],              /* i  : zero-memory filtered code. excitation   */
    float y11[],             /* o  : zero-memory filtered gauss. excitation  */
    float *gain,             /* o  : excitation gain                         */
    float g_corr[],          /* i/o : correlation structure for gain coding  */
    float gain_pit,          /* i : unquantized gain of code                 */
    float tilt_code,         /* i : tilt of code                             */
    const float *Aq,         /* i : quantized LPCs                           */
    float formant_enh_num,   /* i: formant enhancement factor                */
    short *seed_acelp        /*i/o: random seed                              */
);

void gaus_L2_dec(
    float *code,                /* o:   decoded gaussian codevector */
    float tilt_code,
    const float *A,
    float formant_enh_num,
    short *seed_acelp           /*i/o: random seed */
);

float interpolation(      /* o:   interpolated value   */
    const float *x,       /* i:   input vector         */
    const float *win,     /* i:   interpolation window */
    const short frac,     /* i:   fraction             */
    const short up_samp,  /* i:   upsampling factor    */
    const short nb_coef   /* i:   nb of filter coef    */
);

void predict_signal(
    const float excI[],  /* i  : input excitation buffer  */
    float excO[],  /* o  : output excitation buffer */
    const short T0,      /* i  : integer pitch lag        */
    short frac,    /* i  : fraction of lag          */
    const short frac_max,/* i  : max fraction             */
    const short L_subfr  /* i  : subframe size            */
);

void tcx_ltp_encode(
    int tcxltp_on,
    int tcxOnly,
    int tcxMode,
    int L_frame,
    int L_subfr,
    float *speech,
    float *speech_ltp,
    float *wsp,
    int Top,
    int *ltp_param,
    int *ltp_bits,
    int *pitch_int,
    int *pitch_fr,
    float *gain,
    int *pitch_int_past,
    int *pitch_fr_past,
    float *gain_past,
    float *norm_corr_past,
    int last_core,
    int pitmin,
    int pitfr1,
    int pitfr2,
    int pitmax,
    int pitres,
    struct TransientDetection const * pTransientDetection,
    int SideInfoOnly,
    float *A,
    int lpcorder
);

void tcx_ltp_post( int tcxltp_on,
                   short core,
                   int L_frame,
                   int L_frame_core,
                   int delay,
                   float *sig,
                   float *tcx_buf,
                   short tcx_buf_len,
                   int bfi,
                   int pitch_int,
                   int pitch_fr,
                   float gain,
                   int *pitch_int_past,
                   int *pitch_fr_past,
                   float *gain_past,
                   int *filtIdx_past,
                   int pitres,
                   int *pitres_past,
                   float damping,
                   int SideInfoOnly,
                   float *mem_in,
                   float *mem_out,
                   int bitrate
                 );

int tcx_ltp_decode_params(
    int *ltp_param,
    int *pitch_int,
    int *pitch_fr,
    float *gain,
    int pitmin,
    int pitfr1,
    int pitfr2,
    int pitmax,
    int pitres
);


int enc_lsf_tcxlpc(
    int **indices,             /* i  : Ptr to VQ indices */
    Encoder_State *st          /* i/o: Encoder state     */
);

void msvq_enc
(
    const float *const *cb,   /* i  : Codebook (indexed cb[*stages][levels][p])            */
    const int dims[],         /* i  : Dimension of each codebook stage (NULL: full dim.)   */
    const int offs[],         /* i  : Starting dimension of each codebook stage (NULL: 0)  */
    float u[],                /* i  : Vector to be encoded (prediction and mean removed)   */
    const int *levels,        /* i  : Number of levels in each stage                       */
    int maxC,                 /* i  : Tree search size (number of candidates kept from     */
    /*      one stage to the next == M-best)                     */
    int stages,               /* i  : Number of stages                                     */
    float w[],                /* i  : Weights                                              */
    int N,                    /* i  : Vector dimension                                     */
    int maxN,                 /* i  : Codebook dimension                                   */
    int Idx[]                 /* o  : Indices                                              */
);

void msvq_dec
(
    const float *const *cb,   /* i  : Codebook (indexed cb[*stages][levels][p])            */
    const int dims[],         /* i  : Dimension of each codebook stage (NULL: full dim.)   */
    const int offs[],         /* i  : Starting dimension of each codebook stage (NULL: 0)  */
    int stages,               /* i  : Number of stages                                     */
    int N,                    /* i  : Vector dimension                                     */
    int maxN,                 /* i  : Codebook dimension                                   */
    const int Idx[],          /* i  : Indices                                              */
    float *uq,                /* o  : quantized vector                                     */
    Word16 *uq_ind            /* o  : quantized vector (fixed point)                       */
);


void PulseResynchronization(
    float const * src_exc,
    float * dst_exc,
    int nFrameLength,
    int nSubframes,
    float pitchStart,
    float pitchEnd
);


void con_acelp(
    float A[],                    /* i  : coefficients NxAz[M+1]      */
    int coder_type,               /* i  : ACELP coder type            */
    float synth[],                /* i/o: synthesis                   */
    int *pT,                      /* o  :   pitch for all subframe    */
    float *pgainT,                /* o  :   pitch gain for all subfr  */
    float stab_fac,               /* i  : stability of isf            */
    Decoder_State *st,            /* i/o:  coder memory state         */
    float pitch_buffer[],         /* i/o: floating pitch values for each subframe          */
    float *voice_factors,         /* o  : voicing factors                                  */
    float *bwe_exc                /* o  : excitation for SWB TBE                           */
);

void con_tcx(
    Decoder_State* st,            /* i/o: coder memory state          */
    float          synth[]        /* i/o: synth[]                     */
);

int tcxlpc_get_cdk(             /* o  : codebook index              */
    int coder_type                /* (I) GC/VC indicator              */
);

int lsf_msvq_ma_encprm(
    Encoder_State *st,
    int *param_lpc,
    int core,
    int acelp_mode,
    int acelp_midLpc,
    short * bits_param_lpc,
    short no_indices
);

int lsf_msvq_ma_decprm(
    Decoder_State *st,
    int *param_lpc,
    int core,
    int acelp_mode,
    int acelp_midLpc,
    int narrowBand,
    int sr_core
);

int dec_lsf_tcxlpc(
    Decoder_State *st,         /* (I/O) Decoder state   */
    int **indices,             /* (O) Ptr to VQ indices */
    int narrowband,            /* (I) narrowband flag   */
    int cdk                    /* (I) codebook selector */
);

int D_lsf_tcxlpc(
    const int indices[],       /* (I) VQ indices        */
    float lsf_q[],             /* (O) quantized lsf     */
    Word16 lsp_q_ind[],        /* (O) quantized lsp (w/o MA prediction) */
    int narrowband,            /* (I) narrowband flag   */
    int cdk,                   /* (I) codebook selector */
    float mem_MA[]             /* (I) MA memory         */
);

void lsf_update_memory(
    int narrowband,           /* i  : narrowband flag                             */
    const float qlsf[],       /* i  : quantized lsf coefficients                  */
    float old_mem_MA[],       /* i  : MA memory                                   */
    float mem_MA[]            /* o  : updated MA memory                           */
);

int Q_lsf_tcxlpc(
    /* const */ float lsf[],   /* (I) original lsf      */
    float lsf_q[],             /* (O) quantized lsf     */
    Word16 lsp_q_ind[],        /* (O) quantized lsp (w/o MA prediction) */
    int indices[],             /* (O) VQ indices        */
    int narrowband,            /* (I) narrowband flag   */
    int cdk,                   /* (I) codebook selector */
    float mem_MA[],            /* (I) MA memory         */
    int coder_type,            /* (I) acelp extended mode*/
    float *Bin_Ener            /* (I) Spectrum energy*/
);

int E_LPC_lsp_unweight(
    /* const */ float lsp_w[], /* (I): weighted lsp             */
    float lsp_uw[],            /* (O): unweighted lsp           */
    float lsf_uw[],            /* (O): unweighted lsf           */
    float inv_gamma            /* (I): inverse weighting factor */
);

int lsf_ind_is_active(
    const Word16 lsf_q_ind[],
    const float means[],
    int bandwidth,
    int cdk
);

void midlsf_enc(
    float qlsf0[],
    float qlsf1[],
    const float lsf[],
    short *idx,
    int N,
    float *Bin_Ener,
    int narrowBand,
    int sr_core,
    int coder_type
);

void midlsf_dec(
    float qlsf0[],
    float qlsf1[],
    short idx,
    float qlsf[],
    int N,
    int coder_type,
    short *mid_lsf_int,
    short prev_bfi,
    short safety_net
);

void lsf_end_enc(
    Encoder_State *st,
    const float *lsf,
    float *qlsf,
    float *mem_AR,
    float *mem_MA,
    const short nBits,
    const short coder_type_org,
    const short bwidth,
    float *Bin_Ener,
    const float int_fs,
    long core_brate,
    float *streaklimit,
    short *pstreaklen,
    short force_sf,
    short rf_flag,
    short mode2_flag,
    int * lpc_param,
    short * no_stages,
    short * bits_param_lpc,
    short coder_type_raw
);

void lsf_end_dec(
    Decoder_State *st,
    const short coder_type_org,
    const short bwidth,
    const short nBits,
    float *qlsf,
    float *mem_AR,
    float *mem_MA,
    const float int_fs,
    long core_brate,
    unsigned int *p_offset_scale1,
    unsigned int *p_offset_scale2,
    unsigned int *p_offset_scale1_p,
    unsigned int *p_offset_scale2_p,
    short *p_no_scales,
    short *p_no_scales_p,
    short *safety_net,
    int *lpc_param,
    short *LSF_Q_prediction,  /* o  : LSF prediction mode                     */
    int * nb_indices
);

short find_pred_mode(
    const short coder_type,
    const short bwidth,
    const float int_fs,
    short * p_mode_lvq,
    short * p_mode_lvq_p,
    short core_brate
);

void lpc_quantization(
    Encoder_State * st,
    int core,
    int lpcQuantization,
    float lsf_old[],
    const float lsp[],
    const float lspmid[],
    float lsp_q[],
    float lsf_q[],
    float lspmid_q[],
    float mem_MA[],
    float mem_AR[],
    int narrowBand,
    short coder_type,
    int acelp_midLpc,
    int param_lpc[],
    int nbits_lpc[],
    short *seed_acelp,
    int sr_core,
    float *Bin_Ener,
    float *Bin_Ener_old,
    short * bits_param_lpc,
    short * no_param_lpc
);

void lpc_unquantize(
    Decoder_State * st,
    float *lsfold,
    float *lspold,
    float *lsf,
    float *lsp,
    int lpcQuantization,
    int *param_lpc,
    int numlpc,
    int core,
    float *mem_MA,
    float *lspmid,
    float *lsfmid,
    short coder_type,
    int acelp_midLpc,
    int narrow_band,
    short *seed_acelp,
    int sr_core,
    short *mid_lsf_int,
    short prev_bfi,
    short *LSF_Q_prediction,  /* o  : LSF prediction mode                     */
    short *safety_net
);

void dlpc_bfi(
    int L_frame,
    float *lsf_q,            /* o  : quantized lsfs                         */
    const float *lsfold,           /* i  : past quantized lsf                     */
    const short last_good,         /* i  : last good received frame               */
    const short nbLostCmpt,        /* i  : counter of consecutive bad frames      */
    float mem_MA[],          /* i/o: quantizer memory for MA model          */
    float mem_AR[],          /* i/o: quantizer memory for MA model          */
    float *stab_fac,         /* i  : lsf stability factor                   */
    float *lsf_adaptive_mean,/* i  : lsf adaptive mean, updated when BFI==0 */
    int   numlpc,            /* i  : Number of division per superframe      */
    float lsf_cng[],
    int   plcBackgroundNoiseUpdated,
    float *lsf_q_cng,        /* o  : quantized lsfs of background noise      */
    float *old_lsf_q_cng,    /* o  : old quantized lsfs for background noise */
    const float lsfBase[]          /* i  : base for differential lsf coding        */
);

int decode_lpc_stoch(short *ptr,
                     int numlpc,
                     int *param_lpc,
                     int mode,
                     int lpc_quant_type,
                     short *isf_model,
                     int *isf_end_nbits,
                     int *isf_mid_nbits
                    );

void lsfi_enc(
    float qisf0[],
    float qisf1[],
    float isf[],
    short *idx
);

void lsfi_dec(
    float qisf0[],
    float qisf1[],
    short idx,
    float qisf[]
);

void lsf_dec_bfi(
    short codec_mode,                 /* i: : codec_mode: MODE1 | MODE2              */
    float *lsf,                       /* o  : estimated LSF vector                   */
    const float *lsfold,              /* i  : past quantized lsf                     */
    float *lsf_adaptive_mean,         /* i  : lsf adaptive mean, updated when BFI==0 */
    const float lsfBase[],            /* i  : base for differential lsf coding       */
    float *mem_MA,                    /* i/o: quantizer memory for MA model          */
    float *mem_AR,                    /* o  : quantizer memory for AR model          */
    float stab_fac,                   /* i  : lsf stability factor                   */
    short last_coder_type,            /* i  : last coder type                        */
    short L_frame,                    /* i  : frame length                           */
    const short last_good,            /* i  : last good received frame               */
    const short nbLostCmpt,           /* i  : counter of consecutive bad frames      */
    int   plcBackgroundNoiseUpdated,  /* i  : background noise already updated?      */
    float *lsf_q_cng,                 /* o  : quantized lsfs of background noise     */
    float *lsf_cng,                   /* i  : long term target for fading to bg noise*/
    float *old_lsf_q_cng,             /* o  : old quantized lsfs for background noise*/
    short Last_GSC_pit_band_idx,      /* i  : AC mode (GSC) - Last pitch band index  */
    short Opt_AMR_WB,                 /* i  : IO flag                                */
    const short MODE1_bwidth          /* i  : coded bandwidth                        */
);

float const * PlcGetlsfBase(
    int const lpcQuantization,
    int const narrowBand,
    int const sr_core
);

void Unified_weighting(
    float Bin_Ener_128[], /* i  : FFT Bin energy 128 bins in two sets     */
    const float lsf[],          /* i  : LSF vector                              */
    float w[],            /* o  : LP weighting filter (numerator)         */
    const short narrowBand,     /* i  : flag for Narrowband                     */
    const short unvoiced,       /* i  : flag for Unvoiced frame                 */
    const short sr_core,        /* i  : sampling rate of core-coder             */
    const int   order           /* i  : LP order                                */
);

short vad_init(
    T_CldfbVadState *st
);

short vad_proc(
    float realValues[16][60],   /* CLDFB real values */
    float imagValues[16][60],   /* CLDFB imag values */
    float *sb_power,            /* Energy of CLDFB data */
    int numBands,               /* number of input bands */
    T_CldfbVadState *vad_st,      /* VAD state */
    short *cldfb_addition,
    short vada_flag
);

void subband_FFT(
    float Sr[16][60],            /*(i) real part */
    float Si[16][60],            /*(i) imag part */
    float *spec_amp              /*(o) spectral amplitude*/
);

int update_decision(
    T_CldfbVadState *st,
    float snr,                         /*(i) frequency domain SNR  */
    float tsnr,                        /*(i) time domain SNR       */
    float frame_energy,                /*(i) current frame energy  */
    float high_eng,                    /*(i) current frame high frequency energy*/
    int vad_flag,
    int music_backgound_f              /*(i) background music flag */
);

void frame_spec_dif_cor_rate(
    float spec_amp[],                  /*(i) spectral amplitude    */
    float pre_spec_low_dif[],          /*(io) low spectrum different*/
    float f_tonality_rate[]            /*(o) tonality rate*/
);

void ltd_stable(
    float frames_power[],               /*(io) energy of several frames*/
    float ltd_stable_rate[],            /*(o) time-domain stable rate*/
    float frame_energy,                 /*(i) current frame energy*/
    int frameloop                       /*(i) number of  frames*/
);

void SNR_calc(
    float frame_sb_energy[],            /*(i) energy of sub-band divided non-uniformly*/
    float sb_bg_energy[],               /*(i) sub-band background energy*/
    float t_bg_energy,                  /*(i) time background energy of several frames*/
    float *snr,                         /*(o) frequency domain SNR */
    float *tsnr,                        /*(o) time domain SNR */
    float frame_energy,                 /*(i) current frame energy */
    int bandwidth                       /*(i) band width*/
);




void background_update(
    T_CldfbVadState *st,
    float frame_energy,            /*(i) current frame energy 2*/
    int   update_flag,             /*(i) current frame update flag*/
    int   music_backgound_f        /*(i) background music flag*/
    , float snr
);

void bg_music_decision(
    T_CldfbVadState *st,
    int   *music_backgound_f,      /*(i) background music flag*/
    float frame_energy             /*(i) current frame energy 1*/
);

void est_energy(
    float sb_power[],                   /*(o) energy of sub-band divided uniformly*/
    float frame_sb_energy[],            /*(o) energy of sub-band divided non-uniformly*/
    float *p_frame_energy,              /*(o) frame energy 1*/
    float *p_frame_energy2,             /*(o) frame energy 2*/
    float *p_high_energy,               /*(o) high frequency energy*/
    int bw                              /*(i) band width*/
);

void spec_center(
    float spec_power[],                 /*(i) energy of sub-band divided uniformly*/
    float sp_center[],                  /*(o) spectral center*/
    int bandwidth                       /*(i) band width*/
);

void spec_flatness(
    float spec_amp[],                   /*(i) spectral amplitude*/
    float smooth_spec_amp[],            /*(i) smoothed spectral amplitude*/
    float sSFM[]                        /*(o) spectral flatness rate*/
);

int vad_decision(
    T_CldfbVadState *st,
    float snr ,                        /*(i) frequency domain SNR */
    float tsnr,                        /*(i) time domain SNR */
    float snr_flux,                    /*(i) average tsnr of several frames*/
    float lt_snr,                      /*(i)long time SNR calculated by fg_energy and bg_energy*/
    float lt_snr_org,                  /*(i)original long time SNR*/
    float lf_snr,                      /*(i) long time frequency domain SNR calculated by l_speech_snr and l_silence_snr*/
    float frame_energy,                /*(i) current frame energy */
    int  music_backgound_f             /*(i) background music flag*/
);

int comvad_decision(
    T_CldfbVadState *st,
    float snr,					            /*(i) frequency domain SNR */
    float tsnr,					            /*(i) time domain SNR */
    float snr_flux,                         /*(i) average tsnr of several frames*/
    float lt_snr,                           /*(i) long time SNR calculated by fg_energy and bg_energy*/
    float lt_snr_org,	                    /*(i) original long time SNR*/
    float lf_snr,                           /*(i) long time frequency domain
												   SNR calculated by l_speech_snr and l_silence_snr*/
    float frame_energy,                     /*(i) current frame energy */
    int   music_backgound_f,                /*(i) background music flag*/
    short *cldfb_addition,
    short vada_flag
);

void calc_snr_flux(
    float tsnr,                         /*(i) time-domain SNR*/
    float pre_snr[],                    /*(io)  time-domain SNR storage*/
    float *snr_flux                     /*(o) average tsnr*/
);

void calc_lt_snr(
    float *lt_snr_org,                 /*(o)original long time SNR*/
    float *lt_snr,                     /*(o)long time SNR calculated by fg_energy and bg_energy*/
    float fg_energy,                   /*(i)foreground energy sum  */
    int   fg_energy_count,             /*(i) number of the foreground energy frame */
    float bg_energy,                   /*(i)background energy sum  */
    int   bg_energy_count,             /*(i) number of the background energy frame */
    int   bw_index,                    /*(i) band width index*/
    float lt_noise_sp_center0          /*(i)long time noise spectral center by 0*/
);

void calc_lf_snr(
    float *lf_snr_smooth,              /*(o) smoothed lf_snr*/
    float *lf_snr,                     /*(o) long time frequency domain SNR calculated by l_speech_snr and l_silence_snr*/
    float l_speech_snr,                /*(i) sum of active frames snr */
    int   l_speech_snr_count,          /*(i) number of the active frame  */
    float l_silence_snr,               /*(i) sum of the nonactive frames snr*/
    int   l_silence_snr_count,         /*(i) number of the nonactive frame */
    int   fg_energy_count,             /*(i) number of the foreground energy frame */
    int   bg_energy_count,             /*(i) number of the background energy frame */
    int   bw_index                     /*(i) band width index*/
);

float construct_snr_thresh(
    float sp_center[],                 /*(i) spectral center*/
    float snr_flux,                    /*(i) snr flux*/
    float lt_snr,                      /*(i) long time time domain snr*/
    float lf_snr,                      /*(i) long time frequency domain snr*/
    int continuous_speech_num,         /*(i) continuous speech number*/
    int continuous_noise_num,          /*(i) continuous noise number*/
    int fg_energy_est_start,           /*(i) whether if estimated energy*/
    int bw_index                       /*(i) band width index*/
);

void createFdCngCom(
    HANDLE_FD_CNG_COM* hFdCngCom
);

void deleteFdCngCom(
    HANDLE_FD_CNG_COM * hFdCngCom
);

void initFdCngCom(
    HANDLE_FD_CNG_COM hs,              /* i/o: Contains the variables related to the FD-based CNG process */
    float scale
);

void initPartitions(
    const int * part_in,
    int npart_in,
    int startBand,
    int stopBand,
    int * part_out,
    int * npart_out,
    int * midband,
    float * psize,
    float * psize_inv,
    int stopBandFR
);

void minimum_statistics(
    int len,                        /* i  : Vector length */
    int lenFFT,                     /* i  : Length of the FFT part of the vectors */
    float * psize,
    float * msPeriodog,             /* i  : Periodograms */
    float * msNoiseFloor,
    float * msNoiseEst,             /* o  : Noise estimates */
    float * msAlpha,
    float * msPsd,
    float * msPsdFirstMoment,
    float * msPsdSecondMoment,
    float * msMinBuf,
    float * msBminWin,
    float * msBminSubWin,
    float * msCurrentMin,
    float * msCurrentMinOut,
    float * msCurrentMinSubWindow,
    int   * msLocalMinFlag,
    int   * msNewMinFlag,
    float * msPeriodogBuf,
    int   * msPeriodogBufPtr,
    HANDLE_FD_CNG_COM st           /* i/o: FD_CNG structure containing all buffers and variables */
);

void generate_comfort_noise_enc(
    Encoder_State *stcod
);

void generate_comfort_noise_dec(
    float ** bufferReal,            /* o  : Real part of input bands */
    float ** bufferImag,            /* o  : Imaginary part of input bands */
    Decoder_State *stdec
);

void generate_comfort_noise_dec_hf(
    float ** bufferReal,                /* o  : Real part of input bands */
    float ** bufferImag,                /* o  : Imaginary part of input bands */
    Decoder_State *stdec
);

void generate_masking_noise(
    float * timeDomainBuffer,           /* i/o: time-domain signal */
    HANDLE_FD_CNG_COM st,               /* i/o: FD_CNG structure containing all buffers and variables */
    short length,
    short core
);

void generate_masking_noise_update_seed(
    HANDLE_FD_CNG_COM st          /* i/o: FD_CNG structure containing all buffers and variables */
);

void generate_masking_noise_mdct(
    float * mdctBuffer,                 /* i/o: time-domain signal */
    HANDLE_FD_CNG_COM st                /* i/o: FD_CNG structure containing all buffers and variables */
);

void apply_scale(
    float *scale,
    int bandwidth,
    int bitrate
);

void compress_range(
    float* in,
    float* out,
    int len
);

/* Apply some dynamic range expansion to undo the compression */
void expand_range(
    float* in,
    float* out,
    int len
);

void bandcombinepow(
    float* bandpow,                 /* i  : Power for each band */
    int    nband,                   /* i  : Number of bands */
    int*   part,                    /* i  : Partition upper boundaries (band indices starting from 0) */
    int    npart,                   /* i  : Number of partitions */
    float* psize_inv,               /* i  : Inverse partition sizes */
    float* partpow                  /* o  : Power for each partition */
);

void scalebands(
    float* partpow,                 /* i  : Power for each partition */
    int*   part,                    /* i  : Partition upper boundaries (band indices starting from 0) */
    int    npart,                   /* i  : Number of partitions */
    int*   midband,                 /* i  : Central band of each partition */
    int    nFFTpart,                /* i  : Number of FFT partitions */
    int    nband,                   /* i  : Number of bands */
    float* bandpow,                 /* o  : Power for each band */
    short  flag_fft_en
);

void getmidbands(
    int*   part,                    /* i  : Partition upper boundaries (band indices starting from 0) */
    int    npart,                   /* i  : Number of partitions */
    int*   midband,                 /* o  : Central band of each partition */
    float* psize,                   /* o  : Partition sizes */
    float* psize_inv                /* o  : Inverse of partition sizes */
);

void AnalysisSTFT(
    const float *  timeDomainInput,
    float *  fftBuffer,             /* o  : FFT bins */
    HANDLE_FD_CNG_COM st           /* i/o: FD_CNG structure containing all buffers and variables */
);

void SynthesisSTFT(
    float * fftBuffer,
    float * timeDomainOutput,
    float * olapBuffer,
    const float * olapWin,
    int tcx_transition, HANDLE_FD_CNG_COM st
);

void rand_gauss(
    float *x,
    short *seed
);

void lpc_from_spectrum(
    float* powspec,
    int start,
    int stop,
    int fftlen,
    const float *fftSineTab,
    float *A,
    float preemph_fac
);

void createFdCngDec(
    HANDLE_FD_CNG_DEC* hFdCngDec
);

void deleteFdCngDec(
    HANDLE_FD_CNG_DEC * hFdCngDec
);

void initFdCngDec(
    HANDLE_FD_CNG_DEC hs,            /* i/o: Contains the variables related to the FD-based CNG process */
    float scale                    /* i: scaling factor */
);
void configureFdCngDec(
    HANDLE_FD_CNG_DEC hs,         /* i/o: Contains the variables related to the FD-based CNG process */
    short bandwidth,
    int bitrate,
    short L_frame
);

void configure_noise_estimation_dec(
    HANDLE_FD_CNG_DEC st   /* i/o: FD_CNG structure containing all buffers and variables */
);

void ApplyFdCng(
    float * timeDomainInput,
    float ** realBuffer,            /* i/o: Real part of the buffer */
    float ** imagBuffer,            /* i/o: Imaginary part of the buffer */
    HANDLE_FD_CNG_DEC st,           /* i/o: FD_CNG structure containing all buffers and variables */
    unsigned char m_frame_type,     /* i  : Type of frame at the decoder side */
    Decoder_State *stdec,
    const int concealWholeFrame,    /* i  : binary flag indicating frame loss                     */
    short is_music
);

void perform_noise_estimation_dec(
    const float * timeDomainInput,
    HANDLE_FD_CNG_DEC st           /* i/o: FD_CNG structure containing all buffers and variables */
);

void FdCng_decodeSID(
    Decoder_State *st               /* i/o: decoder state structure */
);

void FdCng_exc(
    HANDLE_FD_CNG_COM hs,
    short *CNG_mode,
    short L_frame,
    float *lsp_old,
    short first_CNG,
    float *lsp_CNG,
    float *Aq,                    /* o:   LPC coeffs */
    float *lsp_new,               /* o:   lsp  */
    float *lsf_new,               /* o:   lsf  */
    float *exc,                   /* o:   LP excitation   */
    float *exc2,                  /* o:   LP excitation   */
    float *bwe_exc                /* o:   LP excitation for BWE */
);

void noisy_speech_detection(
    const short vad,
    const float * ftimeInPtr,           /* i  : input time-domain frame                  */
    const int     frameSize,            /* i  : frame size                               */
    const float * msNoiseEst,           /* i  : noise estimate over all critical bands   */
    const float * psize,                /* i  : partition sizes                          */
    const int     nFFTpart,             /* i  : Number of partitions taken into account  */
    float *lp_noise,              /* i/o: long term total Noise energy average     */
    float *lp_speech,             /* i/o: long term active speech energy average   */
    short *flag_noisy_speech
);


void createFdCngEnc(
    HANDLE_FD_CNG_ENC* hFdCngEnc
);

void deleteFdCngEnc(
    HANDLE_FD_CNG_ENC * hFdCngEnc
);

void configureFdCngEnc(
    HANDLE_FD_CNG_ENC hs,         /* i/o: Contains the variables related to the FD-based CNG process */
    short bandwidth,
    int bitrate
);

void initFdCngEnc(
    HANDLE_FD_CNG_ENC hs,            /* i/o: Contains the variables related to the FD-based CNG process */
    int input_Fs,                     /* i: input signal sampling frequency in Hz */
    float scale                    /* i: scaling factor */
);

void resetFdCngEnc(
    Encoder_State * st
);

void perform_noise_estimation_enc(
    float *band_energies,               /* i: energy in critical bands without minimum noise floor E_MIN */
    float *enerBuffer,
    HANDLE_FD_CNG_ENC st
);

void AdjustFirstSID(
    int npart,
    float * msPeriodog,
    float * energy_ho,
    float * msNoiseEst,
    float * msNoiseEst_old,
    short * active_frame_counter,
    Encoder_State *stcod
);

void FdCng_encodeSID(
    HANDLE_FD_CNG_ENC st,                      /* i/o: FD_CNG structure containing all buffers and variables */
    Encoder_State *corest,
    float  preemph_fac
);

void GetParameters(
    ParamsBitMap const * paramsBitMap,
    int nParams,
    void const * pParameter,
    int ** pStream,
    int * pnSize,
    int * pnBits
);

void SetParameters(
    ParamsBitMap const * paramsBitMap,
    int nParams,
    void * pParameter,
    int const ** pStream,
    int * pnSize
);

void WriteToBitstream(
    ParamsBitMap const * paramsBitMap,
    int nParams,
    int const ** pStream,
    int * pnSize, Encoder_State *st,
    int * pnBits
);

void ReadFromBitstream(
    ParamsBitMap const * paramsBitMap,
    int nArrayLength,
    Decoder_State *st,
    int ** pStream,
    int * pnSize
);

void const * GetTnsFilterOrder(void const * p, int index, int * pValue);
void * SetTnsFilterOrder(void * p, int index, int value);
void const * GetNumOfTnsFilters(void const * p, int index, int * pValue);
void * SetNumOfTnsFilters(void * p, int index, int value);
void const * GetTnsEnabled(void const * p, int index, int * pValue);
void * SetTnsEnabled(void * p, int index, int value);
void const * GetTnsEnabledSingleFilter(void const * p, int index, int * pValue);
void * SetTnsEnabledSingleFilter(void * p, int index, int value);
void const * GetTnsFilterCoeff(void const * p, int index, int * pValue);
void * SetTnsFilterCoeff(void * p, int index, int value);

int GetSWBTCX10TnsFilterCoeffBits(int value, int index);
int EncodeSWBTCX10TnsFilterCoeff(int value, int index);
int DecodeSWBTCX10TnsFilterCoeff(Decoder_State *st, int index, int * pValue);
int GetSWBTCX20TnsFilterCoeffBits(int value, int index);
int EncodeSWBTCX20TnsFilterCoeff(int value, int index);
int DecodeSWBTCX20TnsFilterCoeff(Decoder_State *st, int index, int * pValue);

int GetWBTCX20TnsFilterCoeffBits(int value, int index);
int EncodeWBTCX20TnsFilterCoeff(int value, int index);
int DecodeWBTCX20TnsFilterCoeff(Decoder_State *st, int index, int * pValue);

int GetTnsFilterOrderBitsSWBTCX10(int value, int index);
int EncodeTnsFilterOrderSWBTCX10(int value, int index);
int DecodeTnsFilterOrderSWBTCX10(Decoder_State *st, int index, int * pValue);
int GetTnsFilterOrderBitsSWBTCX20(int value, int index);
int EncodeTnsFilterOrderSWBTCX20(int value, int index);
int DecodeTnsFilterOrderSWBTCX20(Decoder_State *st, int index, int * pValue);
int GetTnsFilterOrderBits(int value, int index);
int EncodeTnsFilterOrder(int value, int index);
int DecodeTnsFilterOrder(Decoder_State *st, int index, int * pValue);

void ResetTnsData(
    STnsData * pTnsData
);

void ClearTnsFilterCoefficients(
    STnsFilter * pTnsFilter
);

TNS_ERROR InitTnsConfiguration(
    int nSampleRate,
    int frameLength,
    STnsConfig * pTnsConfig,
    int igfStopFreq,
    int bitrate
);

int DetectTnsFilt(
    STnsConfig const * pTnsConfig,
    float const pSpectrum[],
    STnsData * pTnsData,
    float* predictionGain
);

TNS_ERROR ApplyTnsFilter(
    STnsConfig const * pTnsConfig,
    STnsData const * pTnsData,
    float spectrum[],
    int fIsAnalysis
);

int ITF_Detect(
    float const pSpectrum[],
    short int startLine,
    short int stopLine,
    int maxOrder,
    float* A,
    float* predictionGain,
    int* curr_order
);

TNS_ERROR ITF_Apply(
    float spectrum[],
    short int startLine,
    short int stopLine,
    const float* A,
    int curr_order
);

TNS_ERROR EncodeTnsData(
    STnsConfig const * pTnsConfig,
    STnsData const * pTnsData,
    int * stream,
    int * pnSize,
    int * pnBits
);

int DecodeTnsData(
    STnsConfig const * pTnsConfig,
    int const * stream,
    int * pnSize,
    STnsData * pTnsData
);

TNS_ERROR WriteTnsData(
    STnsConfig const * pTnsConfig,
    int const * stream,
    int * pnSize,
    Encoder_State *st,
    int * pnBits
);

TNS_ERROR ReadTnsData(
    STnsConfig const * pTnsConfig,
    Decoder_State * st,
    int * pnBits,
    int * stream,
    int * pnSize
);

void cldfbAnalysis (
    const float                 *timeIn,              /* i  : time buffer */
    float                       **realBuffer,         /* o  : real value buffer */
    float                       **imagBuffer,         /* o  : imag value buffer */
    int                         samplesToProcess,     /* i  : number of input samples */
    HANDLE_CLDFB_FILTER_BANK    h_cldfb               /* i  : filterbank state */
);

void cldfbSynthesis (
    float                      **realBuffer,         /* i  : real values */
    float                      **imagBuffer,         /* i  : imag values */
    float                       *timeOut,            /* o  : synthesized output */
    int                          samplesToProcess,   /* i  : number of samples */
    HANDLE_CLDFB_FILTER_BANK     h_cldfb             /* i  : filter bank state */
);

void analysisCldfbEncoder (
    Encoder_State *st,                    /* i/o: encoder state structure                    */
    const float *timeIn,
    int samplesToProcess,
    float realBuffer[CLDFB_NO_COL_MAX][CLDFB_NO_CHANNELS_MAX],
    float imagBuffer[CLDFB_NO_COL_MAX][CLDFB_NO_CHANNELS_MAX],
    float *ppBuf_Ener
);

int openCldfb (
    HANDLE_CLDFB_FILTER_BANK *h_cldfb,   /* i/o : filter bank handle */
    CLDFB_TYPE type,                     /* i   : analysis or synthesis */
    int samplerate                       /* i   : max samplerate to oeprate */
);

void resampleCldfb (
    HANDLE_CLDFB_FILTER_BANK hs,  /* i/o : filter bank handle */
    int newSamplerate             /* i   : new samplerate to operate */
);

int cldfb_save_memory(
    HANDLE_CLDFB_FILTER_BANK hs  /* i/o : filter bank handle */
);

int cldfb_restore_memory(
    HANDLE_CLDFB_FILTER_BANK hs   /* i/o : filter bank handle */
);

int cldfb_reset_memory (
    HANDLE_CLDFB_FILTER_BANK hs   /* i/o : filter bank handle */
);

void deleteCldfb(
    HANDLE_CLDFB_FILTER_BANK *h_cldfb  /* i/o : filter bank handle */
);

void fft_cldfb (
    float *data,                 /* i/o: input/output vector */
    int size                     /* size of fft operation */
);

void BITS_ALLOC_init_config_acelp(
    int bit_rate,
    int narrowBand,
    int nb_subfr,
    ACELP_config *acelp_cfg        /*o:  configuration structure of ACELP*/
);

int BITS_ALLOC_config_acelp(
    const int   bits_frame,         /* i  : remaining bit budget for the frame  */
    const short coder_type,         /* i  : acelp extended mode index           */
    ACELP_config *acelp_cfg,         /* i/o: configuration structure of ACELP    */
    const short narrowband,         /* i  : narrowband flag                     */
    const short nb_subfr            /* i  : number of subframes                 */
);

void FEC_clas_estim(
    const float *syn,
    const float *pitch,               /* i  : pitch values for each subframe                 */
    const short L_frame,              /* i  : length of the frame                            */
    const short coder_type,           /* i  : coder type                                     */
    const short codec_mode,           /* i  : codec mode                                     */
    float *mem_syn_clas_estim,  /* i/o: memory of the synthesis signal for frame class estimation */
    short *clas,                /* i/o: frame classification                           */
    float *lp_speech,           /* i/o: long term active speech energy average         */
    long  const bitrate,              /* i  : Decoded bitrate                                */
    short const Opt_AMR_WB,           /* i  : flag indicating AMR-WB IO mode                 */
    short *decision_hyst,       /* i/o: hysteresis of the music/speech decision        */
    short *locattack,           /* i/o: detection of attack (mainly to localized speech burst) */
    short *UV_cnt,              /* i/o: number of consecutives frames classified as UV */
    float *LT_UV_cnt,           /* i/o: long term consecutives frames classified as UV */
    float *Last_ener,           /* i/o: last_energy frame                              */
    short *amr_io_class,        /* i/o: classification for AMR-WB IO mode              */
    float *lt_diff_etot,        /* i/o: long-term total energy variation               */
    float *class_para,          /* o  : classification para. fmerit1                   */
    const float LTP_Gain,             /* i  :                                                */
    const int   narrowBand,           /* i  :                                                */
    const SIGNAL_CLASSIFIER_MODE mode,/* i  :                                                */
    const int   bfi,                  /* i  :                                                */
    const float preemph_fac,          /* i  :                                                */
    const int   tcxonly,              /* i  :                                                */
    const int   last_core_brate       /* i  : last core bitrate                              */
);

void InitTransientDetection(
    int nFrameLength,
    int nTCXDelay,
    struct TransientDetection * pTransientDetection
);

void RunTransientDetection(
    float const * input,
    int nSamplesAvailable,
    struct TransientDetection * pTransientDetection
);

float GetTCXAvgTemporalFlatnessMeasure(
    struct TransientDetection const * pTransientDetection,
    int nCurrentSubblocks,
    int nPrevSubblocks
);

float GetTCXMaxenergyChange(
    struct TransientDetection const * pTransientDetection,
    const int isTCX10,
    const int nCurrentSubblocks,
    const int nPrevSubblocks
);

void SetTCXModeInfo(
    Encoder_State *st,                  /* i/o: encoder state structure                  */
    struct TransientDetection const * pTransientDetection,
    short * tcxModeOverlap
);

void TCX_MDCT(
    float const *x,
    float *y,
    int l,
    int m,
    int r
);

void TCX_MDST(
    float const *x,
    float *y,
    int l,
    int m,
    int r
);

void TCX_MDCT_Inverse(
    float *x,
    float *y,
    int l,
    int m,
    int r
);

void post_decoder(
    Decoder_State *st,
    const short coder_type,           /* i  : coder type                  */
    float synth_buf[],
    const float pit_gain[],
    const int   pitch[],
    float signal_out[],
    float bpf_noise_buf[]
);

float bass_pf_enc(
    const float *orig,      /* (i) : 12.8kHz original signal                     */
    const float *syn,       /* (i) : 12.8kHz synthesis to postfilter             */
    const float pitch_buf[],/* (i) : Pitch gain for all subframes (gainT_sf[16]) */
    const float gainT_sf[], /* (i) : Pitch gain for all subframes (gainT_sf[16]) */
    const short l_frame,    /* (i) : frame length (should be multiple of l_subfr)*/
    const short l_subfr_in, /* (i) : sub-frame length (60/64)                    */
    float mem_bpf[],        /* i/o : memory state [2*L_FILT]                     */
    float mem_error_bpf[],  /* i/o : memory state [2*L_FILT]                     */
    int *gain_factor_param, /* (o) : quantized gain factor                       */
    const short mode,       /* (i) : coding mode of adapt bpf                    */
    float *mem_deemph_err,  /* i/o : Error deemphasis memory                     */
    float *lp_ener          /* i/o : long_term error signal energy               */
);

void cldfb_synth_set_bandsToZero(
    Decoder_State *st,
    float **rAnalysis,
    float **iAnalysis,
    const short nTimeSlots
);

void longadd(
    unsigned short a[],
    unsigned short b[],
    int lena,
    int lenb
);

void longshiftright(
    unsigned short a[],
    int b,
    unsigned short d[],
    int lena,
    int lend
);

void longshiftleft(
    unsigned short a[],
    int b,
    unsigned short d[],
    int len
);

void open_decoder_LPD(
    Decoder_State *st,
    const int bit_rate,
    const int bandwidth
);

void update_decoder_LPD_cng(
    Decoder_State *st,
    const short coder_type,           /* i  : coder type                  */
    float *timeDomainBuffer,
    float *A,
    float *bpf_noise_buf
);

void reconfig_decoder_LPD(
    Decoder_State *st,
    int bits_frame,
    int bandwidth,
    int bitrate,
    int L_frame_old
);

void mode_switch_decoder_LPD(
    Decoder_State *st,
    int bandwidth_in,
    int bitrate,
    int frame_size_index
);

void dec_acelp_tcx_frame(
    Decoder_State *st,                /* i/o: encoder state structure             */
    short *coder_type,        /* o  : coder type                          */
    short *concealWholeFrame, /* i/o: concealment flag                    */
    float *output,            /* o  : synthesis                           */
    float *bpf_noise_buf,     /* i/o: BPF noise buffer                    */
    float * pcmbufFB,
    float bwe_exc_extended[], /* i/o: bandwidth extended excitation       */
    float *voice_factors,     /* o  : voicing factors                     */
    float pitch_buf[]         /* o  : floating pitch for each subframe    */
);

void decoder_LPD(
    float signal_out[],               /* output: signal with LPD delay (7 subfrs) */
    float signal_outFB[],
    short *total_nbbits,        /* i/o:    number of bits / decoded bits    */
    Decoder_State *st,                /* i/o:    decoder memory state pointer     */
    float *bpf_noise_buf,       /* i/o: BPF noise buffer                    */
    short bfi,                  /* i  : BFI flag                            */
    short *bitsRead,            /* o  : number of read bits                 */
    short *coder_type,          /* o  : coder type                          */
    int   param[],              /* o  : buffer of parameters                */
    float *pitch_buf,           /* i/o: floating pitch values for each subfr*/
    float *voice_factors,       /* o  : voicing factors                     */
    float *ptr_bwe_exc          /* o  : excitation for SWB TBE              */
);


int tcxGetNoiseFillingTilt(
    float A[],
    int L_frame,
    int mode,
    float *noiseTiltFactor
);

void tcxFormantEnhancement(
    float xn_buf[],
    float gainlpc[],
    float spectrum[],
    int L_frame
);

void tcxInvertWindowGrouping(
    TCX_config *tcx_cfg,
    float xn_buf[],
    float spectrum[],
    int L_frame,
    int fUseTns,
    int last_core,
    int index,
    int frame_cnt,
    int bfi
);

void lerp(
    float *f,
    float *f_out,
    int bufferNewSize,
    int bufferOldSize
);

void coderLookAheadInnovation(
    const float A[],                   /* input: coefficients NxAz[M+1]   */
    int *pT,                           /* out:   pitch   */
    HANDLE_PLC_ENC_EVS hPlc_Ext,
    float *speechLookAhead,
    float *old_exc,
    int L_subfr,
    int L_frame
);

void encoderSideLossSimulation(
    Encoder_State *st,
    HANDLE_PLC_ENC_EVS hPlc_Ext,
    float *isf_q,
    float stab_fac,
    int calcOnlyISF,
    int L_frame
);

void enc_prm_side_Info(
    HANDLE_PLC_ENC_EVS hPlc_Ext,
    Encoder_State *st
);

void GplcTcxEncSetup(
    Encoder_State *st,
    HANDLE_PLC_ENC_EVS hPlc_Ext
);

short encSideSpecPowDiffuseDetector(
    float *isf_ref,
    float *isf_con,
    int sr_core,
    float *prev_isf4_mean,
    short sw
    ,short coder_type
);

void updateSpecPowDiffuseIdx(
    Encoder_State *st,
    const float gain_pitch_buf[],   /* i  : gain pitch values   */
    const float gain_code_buf[]     /* i  : gain pitch values   */
);

void getLookAheadResSig(
    float *speechLookAhead,
    const float *A,
    float *res,
    int L_frame,
    int L_subfr,
    int m,
    int numSubFrame
);

void updatelsfForConcealment(
    HANDLE_PLC_ENC_EVS decState,
    float *lsf
);

void getConcealedLP(
    HANDLE_PLC_ENC_EVS memDecState,
    float *AqCon,
    const float xsfBase[],
    const float sr_core,
    int last_good,
    int L_frame
);


void getConcealedlsf(
    HANDLE_PLC_ENC_EVS memDecState,
    const float lsfBase[],
    int L_frame, int last_good
);

void RecLpcSpecPowDiffuseLc(
    float *ispq,
    float *isp_old,
    float *isfq,
    Decoder_State *st
    , int reset_q
);

void modify_lsf(
    float *lsf,
    const short n,
    const int sr_core
    , int reset_q
);

void open_PLC_ENC_EVS(
    HANDLE_PLC_ENC_EVS hPlcExt,
    int sampleRate
);

void gPLC_encInfo(
    HANDLE_PLC_ENC_EVS self,
    const int brate,
    const int bwidth,
    const short last_clas,
    const int coder_type
);

void resetTecDec(
    HANDLE_TEC_DEC hTecDec
);

void updateTecBuf(
    float** pCldfbRealSrc,
    float** pCldfbImagSrc,
    float** pCldfbRealTgt,
    float** pCldfbImagTgt,
    int noCols,
    int lowSubband
);

void calcGainTemp_TBE(
    float** pCldfbRealSrc,
    float** pCldfbImagSrc,
    float* loBuffer,
    int startPos,           /*!<  Start position of the current envelope. */
    int stopPos,            /*!<  Stop position of the current envelope. */
    int lowSubband,         /* lowSubband */
    float* pGainTemp,
    short code
);

void procTecTfa_TBE(
    float *hb_synth,
    float *gain,
    short flat_flag,
    short last_core,
    int L_subfr,
    short code
);

void resetTecEnc(
    HANDLE_TEC_ENC hTecEnc,
    int flag
);

void calcHiEnvLoBuff(
    int noCols,
    const int* pFreqBandTable,  /*!<  freqbandTable. */
    int nSfb,                   /*!<  Number of scalefactors. */
    float** pYBuf,
    float* loBuf,
    float* hiTempEnv
);

void calcLoEnvCheckCorrHiLo(
    int noCols,
    const int* pFreqBandTable,  /*!<  freqbandTable. */
    float* loBuf,
    float* loTempEnv,
    float* loTempEnv_ns,
    float* hiTempEnv,
    int* corr_flag              /* 0 for original,  1 for TEC */
);

void tfaCalcEnv(
    const float* shb_speech,
    float* enr
);

short tfaEnc_TBE(
    float* enr,
    short last_core,
    float* voicing,
    float* pitch_buf
);

void tecEnc_TBE(
    int* corrFlag,
    const float* voicing,
    short coder_type
);

void set_TEC_TFA_code(
    const short corrFlag,
    short* tec_flag,
    short* tfa_flag
);

float Damping_fact(
    const short coder_type,/* i  : ACELP core coding mode                           */
    int nbLostCmpt,        /* i  : compt for number of consecutive lost frame       */
    short last_good,       /* i  : class of last good received frame                */
    float stab_fac,        /* i  : LSF stability factor                             */
    float *lp_gainp,       /* i/o: low passed pitch gain used for concealment       */
    int core               /* i  : current core: ACELP = 0, TCX20 = 1, TCX10 = 2    */
);

void fer_energy(
    const int L_frame,    /* i  : frame length                           */
    const short clas,     /* i  : frame classification                   */
    const float synth[],  /* i  : synthesized speech at Fs = 12k8 Hz     */
    const float pitch,    /* i  : pitch period                           */
    float *enr,     /* o  : pitch-synchronous or half_frame energy */
    const short useOffset /* i  : speech pointer offset (0 or L_FRAME)   */
);

float getLevelSynDeemph(
    float const h1Init[],     /* i: input value or vector to be processed */
    float const A[],          /* i: LPC coefficients                      */
    int   const lenLpcExc,    /* i: length of the LPC excitation buffer   */
    float const preemph_fac,  /* i: preemphasis factor                    */
    int   const numLoops      /* i: number of loops                       */
);

void genPlcFiltBWAdap(
    int   const sr_core,     /* i: core sampling rate                                         */
    float*      lpFiltAdapt, /* o: filter coefficients for filtering codebooks in case of flc */
    int   const type,        /* i: type of filter, either 0 : lowpass or 1 : highpass         */
    float const alpha        /* i: fade out factor [0 1) used decrease filter tilt            */
);

void highPassFiltering(
    const short last_good,      /* i:   last classification type                           */
    const int   L_buffer,       /* i:   buffer length                                      */
    float       exc2[],         /* i/o: unvoiced excitation before the high pass filtering */
    const float hp_filt[],      /* i:   high pass filter coefficients                      */
    const int   l_fir_fer       /* i:   high pass filter length                            */
);

int GetPLCModeDecision(
    Decoder_State *st                   /* i/o:    decoder memory state pointer */
);

void addBassPostFilter (
    const float *harm_timeIn,
    int samplesToProcess,
    float **rAnalysis,
    float **iAnalysis,
    HANDLE_CLDFB_FILTER_BANK cldfb
);

TONALMDCTCONCEAL_ERROR TonalMDCTConceal_Init(
    TonalMDCTConcealPtr self,
    unsigned int samplesPerBlock,
    unsigned int nSamplesCore,
    unsigned int nScaleFactors,
    TCX_config * tcx_cfg
);

TONALMDCTCONCEAL_ERROR TonalMDCTConceal_SaveFreqSignal(
    TonalMDCTConcealPtr self,
    float const * mdctSpectrum,
    unsigned int numSamples,
    unsigned int nNewSamplesCore,
    float const *scaleFactors
);

TONALMDCTCONCEAL_ERROR TonalMDCTConceal_UpdateState(
    TonalMDCTConcealPtr self,
    int numSamples,
    float pitchLag,
    int badBlock,
    int tonalConcealmentActive
);

TONALMDCTCONCEAL_ERROR TonalMDCTConceal_SaveTimeSignal(
    TonalMDCTConcealPtr self,
    float* timeSignal,
    unsigned int numSamples
);

TONALMDCTCONCEAL_ERROR TonalMDCTConceal_Detect(
    TonalMDCTConcealPtr const self,     /*IN */
    float const pitchLag,               /*IN */
    int * const umIndices               /*OUT*/
);

TONALMDCTCONCEAL_ERROR TonalMDCTConceal_Apply(
    TonalMDCTConcealPtr self,     /*IN */
    float* mdctSpectrum           /*OUT*/
);

TONALMDCTCONCEAL_ERROR TonalMDCTConceal_InsertNoise(
    TonalMDCTConcealPtr self,     /*IN */
    float* mdctSpectrum,          /*OUT*/
    int tonalConcealmentActive,
    short * pSeed,                /*IN/OUT*/
    float tiltCompFactor,
    float crossfadeGain,
    int crossOverFreq
);

void DetectTonalComponents(
    unsigned short int indexOfTonalPeak[],
    unsigned short int lowerIndex[],
    unsigned short int upperIndex[],
    unsigned int * pNumIndexes,
    float lastPitchLag,
    float currentPitchLag,
    float const lastMDCTSpectrum[],
    float const scaleFactors[],
    float const secondLastPowerSpectrum[],
    unsigned int nSamples,
    unsigned int nSamplesCore,
    float floorPowerSpectrum
);

void RefineTonalComponents(
    unsigned short int indexOfTonalPeak[],
    unsigned short int lowerIndex[],
    unsigned short int upperIndex[],
    float phaseDiff[],
    float phases[],
    unsigned int * pNumIndexes,
    float lastPitchLag,
    float currentPitchLag,
    float const lastMDCTSpectrum[],
    float const scaleFactors[],
    float const secondLastPowerSpectrum[],
    unsigned int nSamples,
    unsigned int nSamplesCore,
    float floorPowerSpectrum
);

void concealment_init(
    int N,
    void *_plcInfo
);

void concealment_decode(
    int pre_transient,
    float *invkoef,
    void *_plcInfo
);

void concealment_update(
    int bfi,
    int core,
    int harmonic,
    float *invkoef,
    void *_plcInfo
);

void concealment_update2(
    float *outx_new,
    void *_plcInfo,
    int FrameSize
);

void concealment_signal_tuning(
    int bfi,
    int curr_mode,
    float *outx_new,
    void *_plcInfo,
    int nbLostCmpt,
    int pre_bfi,
    float *OverlapBuf,
    int past_core_mode,
    float *outdata2,
    Decoder_State *st
);

void waveform_adj2(
    float *overlapbuf,
    float *outx_new,
    float *data_noise,
    float *outx_new_n1,
    float *nsapp_gain,
    float *nsapp_gain_n,
    float *recovery_gain,
    float step_concealgain,
    int    pitch,
    int    Framesize,
    int    delay,
    int    bfi_cnt,
    int    bfi
);


float SFM_Cal(
    float const fcoef[],
    int n
);

void set_state(
    int *state,
    int num,
    int N
);

void Shellsort_float(
    float *in,
    int n
);

int RFFTN(
    float *afftData,
    const float* trigPtr,
    int len,
    int isign
);


void DoFFT(
    float * re2,
    float * im2,
    short length
);


short getTcxonly(
    const int bitrate
);

short getTnsAllowed(
    const int bitrate,
    const short igf
);

short getCtxHm(
    const int bitrate,
    const short rf_flag
);

short getResq(
    const int bitrate
);

short getRestrictedMode(
    const int bitrate,
    const short Opt_AMR_WB
);

short getMdctWindowLength(
    const float fscale
);

short sr2fscale(
    const int sr
);

int getCoreSamplerateMode2(
    const int bitrate,
    const int bandwidth,
    const short rf_mode
);

float getTcxBandwidth(
    const int bandwidth
);

short getIgfPresent(
    const int bitrate,
    const int bandwidth,
    const short rf_mode
);

short getCnaPresent(
    const int bitrate,
    const int bandwidth
);

short getTcxLtp(
    const int sr_core
);

short initPitchLagParameters(
    const int sr_core,
    int *pit_min,
    int *pit_fr1,
    int *pit_fr1b,
    int *pit_fr2,
    int *pit_max
);

void attenuateNbSpectrum(
    int L_frame,
    float *spectrum
);

void SetModeIndex(
    Encoder_State *st,
    const long total_brate,
    const short bwidth
);

short getNumTcxCodedLines(
    const short bwidth
);

short getTcxLpcShapedAri(
    const int total_brate,
    const short bwidth,
    const short rf_mode
);

void IGFEncApplyMono(
    const IGF_ENC_INSTANCE_HANDLE                   hInstance,          /**< in:     | instance handle of IGF Encoder                         */
    const int                                       igfGridIdx,         /**< in:     | IGF grid index                                         */
    Encoder_State                                  *st,                 /**< in:     | Encoder state                                          */
    float                                          *pMDCTSpectrum,      /**< in:     | MDCT spectrum                                          */
    float                                          *pPowerSpectrum,     /**< in:     | MDCT^2 + MDST^2 spectrum, or estimate                  */
    int                                             isTCX20,            /**< in:     | flag indicating if the input is TCX20 or TCX10/2xTCX5  */
    int                                             isTNSActive,        /**< in:     | flag indicating if the TNS is active                   */
    int                                             last_core_acelp     /**< in:     | indictaor if last frame was acelp coded                */
);

void IGFEncConcatenateBitstream(
    const IGF_ENC_INSTANCE_HANDLE                   hInstance,          /**< in:     | instance handle of IGF Encoder                 */
    short                                           bsBits,             /**< in:     | number of IGF bits written to list of indices  */
    short                                          *next_ind,           /**< in/out: | pointer to actual bit indice                   */
    short                                          *nb_bits,            /**< in/out: | total number of bits already written           */
    Indice                                         *ind_list            /**< in:     | pointer to list of indices                     */
);

void IGFEncResetTCX10BitCounter(
    const IGF_ENC_INSTANCE_HANDLE                   hInstance           /**< in:     | instance handle of IGF Encoder */
);

void IGFEncSetMode(
    const IGF_ENC_INSTANCE_HANDLE                   hInstance,          /**< in:     | instance handle of IGF Encoder */
    const int                                       bitRate,            /**< in:     | encoder bitrate                */
    const int                                       mode                /**< in:     | encoder bandwidth mode         */
    , const int                                       f_mode              /**< in:     | flag to signal the RF mode */
);

int IGFEncWriteBitstream(                                               /**< out:    | number of bits written per frame                                             */
    const IGF_ENC_INSTANCE_HANDLE                   hInstance,          /**< in:     | instance handle of IGF Encoder                                               */
    void                                           *st,                 /**< in:     | encoder state                                                                */
    int                                            *pBitOffset,         /**< in:     | ptr to bitOffset counter                                                     */
    const int                                       igfGridIdx,         /**< in:     | igf grid index see declaration of IGF_GRID_IDX for details                   */
    const int                                       isIndepFlag         /**< in:     | if 1 frame is independent, 0 = frame is coded with data from previous frame  */
);

int IGFEncWriteConcatenatedBitstream(                                   /**< out:    | total number of bits written   */
    const IGF_ENC_INSTANCE_HANDLE                   hInstance,          /**< in:     | instance handle of IGF Encoder */
    void                                            *st                 /**< in:     | encoder state                  */
);

void IGFDecApplyMono(
    const IGF_DEC_INSTANCE_HANDLE                   hInstance,          /**< in:     | instance handle of IGF Decoder                       */
    float                                          *spectrum,           /**< in/out: | MDCT spectrum                                        */
    const int                                       igfGridIdx,         /**< in:     | in case of CELP->TCX switching, use 1.25 framelength */
    int                                             bfi                 /**< in:     | frame loss == 1, frame good == 0                     */
);

void IGFDecCopyLPCFlatSpectrum(
    const IGF_DEC_INSTANCE_HANDLE                   hInstance,          /**< in:     | instance handle of IGF Decoder     */
    const float                                    *pSpectrumFlat,      /**< in:     | LPC flattend spectrum from TCX dec */
    const int                                       igfGridIdx          /**< in:     | IGF grid index                     */
);

void IGFDecReadData(
    const IGF_DEC_INSTANCE_HANDLE                    hInstance,          /**< in:     | instance handle of IGF Deccoder                      */
    Decoder_State                                   *st,                 /**< in:     | decoder state                                        */
    const int                                        igfGridIdx,         /**< in:     | in case of CELP->TCX switching, use 1.25 framelength */
    const int                                        isIndepFrame        /**< in:     | if 1: arith dec force reset, if 0: no reset          */
);

int IGFDecReadLevel(                                                     /**< out:    | return igfAllZero flag indicating if no envelope is transmitted  */
    const IGF_DEC_INSTANCE_HANDLE                    hInstance,          /**< in:     | instance handle of IGF Deccoder                                  */
    Decoder_State                                   *st,                 /**< in:     | decoder state                                                    */
    const int                                        igfGridIdx,         /**< in:     | in case of CELP->TCX switching, use 1.25 framelength             */
    const int                                        isIndepFrame        /**< in:     | if 1: arith dec force reset, if 0: no reset                      */
);

void IGFDecRestoreTCX10SubFrameData(
    const IGF_DEC_INSTANCE_HANDLE                   hInstance,          /**< in:     | instance handle of IGF Decoder */
    const int                                       subFrameIdx         /**< in:     | index of subframe              */
);

void IGFDecSetMode(
    const IGF_DEC_INSTANCE_HANDLE                   hInstance,          /**< in:     | instance handle of IGF Decoder */
    const int                                       bitRate,            /**< in:     | bitrate                        */
    const int                                       mode,               /**< in:     | bandwidth mode                 */
    const int                                       defaultStartLine,   /**< in:     | default start subband index    */
    const int                                       defaultStopLine,    /**< in:     | default stop subband index     */
    const int                                       rf_mode             /**< in:     | flag to signal the RF mode */
);

void IGFDecStoreTCX10SubFrameData(
    const IGF_DEC_INSTANCE_HANDLE                   hInstance,          /**< in:     | instance handle of IGF Decoder */
    const int                                       subFrameIdx         /**< in:     | index of subframe              */
);

void IGFDecUpdateInfo(
    const IGF_DEC_INSTANCE_HANDLE                   hInstance,          /**< in:     | instance handle of IGF Decoder */
    const int                                       igfGridIdx          /**< in:     | IGF grid index                 */
);

int IGFCommonFuncsIGFConfiguration(                                                      /**< out:    | error value: 0 -> error, 1 -> ok   */
    int                                             bitRate,            /**< in:     | bitrate in bs e.g. 9600 for 9.6kbs */
    int                                             mode,               /**< in:     | bandwidth mode                     */
    H_IGF_INFO                                      hIGFInfo,           /**< out:    | IGF info handle                    */
    int                                             rf_mode             /**< in:     | flag to signal the RF mode */
);

int IGFCommonFuncsIGFGetCFTables(                                       /**< out:    | error value: 0 -> error, 1 -> ok     */
    int                                             bitRate,            /**< in:     | bitrate in bs e.g. 9600 for 9.6kbs   */
    int                                             mode,               /**< in:     | bandwidth mode                       */
    int                                             rf_mode,            /**< in:     | flag to signal the RF mode */
    const unsigned short                          **cf_se00,            /**< out:    | CF table for t == 0 and f == 0       */
    const unsigned short                          **cf_se01,            /**< out:    | CF table for t == 0 and f == 1       */
    short                                          *cf_off_se01,        /**< out:    | offset for CF table above            */
    const unsigned short                          **cf_se02,            /**< out:    | CF tables for t == 0 and f >= 2      */
    const short                                   **cf_off_se02,        /**< out:    | offsets for CF tables above          */
    const unsigned short                          **cf_se10,            /**< out:    | CF table for t == 1 and f == 0       */
    short                                          *cf_off_se10,        /**< out:    | offset for CF table above            */
    const unsigned short                          **cf_se11,            /**< out:    | CF tables for t == 1 and f >= 1      */
    const short                                   **cf_off_se11         /**< out:    | offsets for CF tables above          */
);

void IGFCommonFuncsWriteSerialBit(
    void                                           *st,                 /**< in:     | encoder/decoder state structure  */
    int                                            *pBitOffset,         /**< out:    | bit offset                       */
    int                                             bit                 /**< in:     | value of bit                     */
);

void IGFCommonFuncsWriteSerialBit(
    void                                           *st,                 /**< in:     | encoder/decoder state structure  */
    int                                            *pBitOffset,         /**< out:    | bit offset                       */
    int                                             bit                 /**< in:     | value of bit                     */
);

void IGFSCFEncoderOpen(
    IGFSCFENC_INSTANCE_HANDLE                       hPublicData,        /**< inout: handle to public data */
    int                                             scfCountLongBlock,
    int                                             bitRate,
    int                                             mode
    , int                                             rf_mode
);

int IGFSCFEncoderReset(
    IGFSCFENC_INSTANCE_HANDLE                       hPublicData         /**< inout: handle to public data or NULL in case there was no instance created */
);

int IGFSCFEncoderEncode(
    IGFSCFENC_INSTANCE_HANDLE                       hPublicData,        /**< inout: handle to public data or NULL in case there was no instance created */
    Encoder_State                                  *st,
    int                                             bitCount,           /**< in: offset to the first bit in bitbuffer which should be readed by iisArithDecoderDecode function */
    int                                            *sfe,                /**< in: ptr to an array which contain quantized scalefactor energies */
    int                                             indepFlag,          /**< in: if  1 on input the encoder will be forced to reset,
                                                                                 if  0 on input the encodder will be forced to encode without a reset */
    int                                             doRealEncoding      /**< in: whether the real encoding is needed, otherwise only the number of bits is used */
);

void IGFSCFEncoderSaveContextState(
    IGFSCFENC_INSTANCE_HANDLE                       hPublicData         /**< inout: handle to public data or NULL in case there was no instance created */
);

void IGFSCFEncoderRestoreContextState(
    IGFSCFENC_INSTANCE_HANDLE                       hPublicData         /**< inout: handle to public data or NULL in case there was no instance created */
);

void IGFSCFDecoderOpen(
    IGFSCFDEC_INSTANCE_HANDLE                       hPublicData,        /**< inout: handle to public data */
    int                                             scfCountLongBlock,
    int                                             bitRate,
    int                                             mode,
    int                                             rf_mode
);

void IGFSCFDecoderReset(
    IGFSCFDEC_INSTANCE_HANDLE                       hPublicData         /**< inout: handle to public data or NULL in case there was no instance created */
);

void IGFSCFDecoderDecode(
    IGFSCFDEC_INSTANCE_HANDLE                       hPublicData,        /**< inout: handle to public data or NULL in case there was no instance created */
    Decoder_State                                  *st,                 /**< inout: pointer to decoder state */
    int                                            *sfe,                /**< out: ptr to an array which will contain the decoded quantized coefficients */
    int                                             indepFlag           /**< in: if  1 on input the encoder will be forced to reset,
                                                                                 if  0 on input the encodder will be forced to encode without a reset */
);

short tbe_celp_exc_offset(
    const short T0,               /* i  : Integer pitch */
    const short T0_frac           /* i  : Fractional part of the pitch */
);

void blend_subfr2(
    float *sigIn1,                /* i  : input signal for fade-out */
    float *sigIn2,                /* i  : input signal for fade-in  */
    float *sigOut                 /* o  : output signal             */
);



#endif
