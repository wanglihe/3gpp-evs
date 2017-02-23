/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

/** \file jbm_jb4_circularbuffer.h circular buffer (FIFO) with fixed capacity */

#ifndef JBM_JB4_CIRCULARBUFFER_H
#define JBM_JB4_CIRCULARBUFFER_H JBM_JB4_CIRCULARBUFFER_H

/** handle for circular buffer (FIFO) with fixed capacity */
typedef struct JB4_CIRCULARBUFFER *JB4_CIRCULARBUFFER_HANDLE;
/** type of circular buffer elements */
typedef int JB4_CIRCULARBUFFER_ELEMENT;

/** Creates a circular buffer (FIFO)
 * @param[out] ph pointer to created handle
 * @return 0 if succeeded */
int JB4_CIRCULARBUFFER_Create( JB4_CIRCULARBUFFER_HANDLE *ph );
/** Destroys the circular buffer (FIFO) */
void JB4_CIRCULARBUFFER_Destroy( JB4_CIRCULARBUFFER_HANDLE *ph );
/** Initializes a circular buffer (FIFO) with a fixed maximum allowed number of elements
 * @param[in] capacity maximum allowed number of elements
 * @return 0 if succeeded */
int JB4_CIRCULARBUFFER_Init( JB4_CIRCULARBUFFER_HANDLE h, unsigned int capacity );

int JB4_CIRCULARBUFFER_Enque( JB4_CIRCULARBUFFER_HANDLE h, JB4_CIRCULARBUFFER_ELEMENT element );
int JB4_CIRCULARBUFFER_Deque( JB4_CIRCULARBUFFER_HANDLE h, JB4_CIRCULARBUFFER_ELEMENT *pElement );

/** Returns the first element. */
JB4_CIRCULARBUFFER_ELEMENT JB4_CIRCULARBUFFER_Front( const JB4_CIRCULARBUFFER_HANDLE h );
/** Returns the last element. */
JB4_CIRCULARBUFFER_ELEMENT JB4_CIRCULARBUFFER_Back( const JB4_CIRCULARBUFFER_HANDLE h );

int JB4_CIRCULARBUFFER_IsEmpty( const JB4_CIRCULARBUFFER_HANDLE h );
int JB4_CIRCULARBUFFER_IsFull( const JB4_CIRCULARBUFFER_HANDLE h );
unsigned int JB4_CIRCULARBUFFER_Size( const JB4_CIRCULARBUFFER_HANDLE h );

/** Calculates statistics over all elements: min element
 * @param[out] pMin minimum element */
void JB4_CIRCULARBUFFER_Min( const JB4_CIRCULARBUFFER_HANDLE h, JB4_CIRCULARBUFFER_ELEMENT *pMin );
/** Calculates statistics over all elements: max element
 * @param[out] pMax maximum element */
void JB4_CIRCULARBUFFER_Max( const JB4_CIRCULARBUFFER_HANDLE h, JB4_CIRCULARBUFFER_ELEMENT *pMax );
/** Calculates statistics over a considered fraction of all elements: min element and percentile
 * @param[in]  nElementsToIgnore number of highest elements to ignore for percentile calculation
 * @param[out] pMin minimum element
 * @param[out] pPercentile consideredFraction percentile of the highest elements */
void JB4_CIRCULARBUFFER_MinAndPercentile( const JB4_CIRCULARBUFFER_HANDLE h, unsigned int nElementsToIgnore,
        JB4_CIRCULARBUFFER_ELEMENT *pMin, JB4_CIRCULARBUFFER_ELEMENT *pPercentile );

#endif /* JBM_JB4_CIRCULARBUFFER_H */
