/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

/** \file jbm_jb4_inputbuffer.h RTP input buffer with fixed capacity. */

#ifndef JBM_JB4_INPUTBUFFER_H
#define JBM_JB4_INPUTBUFFER_H JBM_JB4_INPUTBUFFER_H

#include "jbm_types.h"

/** Handle for RTP input buffer with fixed capacity. */
/** Implemented as priority queue using an array based sorted circular buffer. */
typedef struct JB4_INPUTBUFFER *JB4_INPUTBUFFER_HANDLE;
/** type of sorted circular buffer elements */
typedef void* JB4_INPUTBUFFER_ELEMENT;

/** Creates a input buffer
 * @param[out] ph pointer to created handle
 * @return 0 if succeeded */
int JB4_INPUTBUFFER_Create( JB4_INPUTBUFFER_HANDLE *ph );
/** Destroys the input buffer */
void JB4_INPUTBUFFER_Destroy( JB4_INPUTBUFFER_HANDLE *ph );
/** Initializes a input buffer with a fixed maximum allowed number of elements
 * @param[in] capacity maximum allowed number of elements
 * @param[in] function to compare two elements: newElement==arrayElement ? 0 : (newElement>arrayElement ? +1 : -1)
 * @return 0 if succeeded */
int JB4_INPUTBUFFER_Init( JB4_INPUTBUFFER_HANDLE h, unsigned int capacity,
                          int (*compareFunction)( const JB4_INPUTBUFFER_ELEMENT newElement, const JB4_INPUTBUFFER_ELEMENT arrayElement,
                                  bool_t *replaceWithNewElementIfEqual ) );

/** Add an element to the buffer.
 * @return 0 if succeeded, -1 if buffer full, +1 if element with same index already stored. */
int JB4_INPUTBUFFER_Enque( JB4_INPUTBUFFER_HANDLE h, JB4_INPUTBUFFER_ELEMENT element, JB4_INPUTBUFFER_ELEMENT *replacedElement );
int JB4_INPUTBUFFER_Deque( JB4_INPUTBUFFER_HANDLE h, JB4_INPUTBUFFER_ELEMENT *pElement );

/** Returns the first element. */
JB4_INPUTBUFFER_ELEMENT JB4_INPUTBUFFER_Front( const JB4_INPUTBUFFER_HANDLE h );
/** Returns the last element. */
JB4_INPUTBUFFER_ELEMENT JB4_INPUTBUFFER_Back( const JB4_INPUTBUFFER_HANDLE h );
/** Returns the element with the given index (0 means front element). */
JB4_INPUTBUFFER_ELEMENT JB4_INPUTBUFFER_Element( const JB4_INPUTBUFFER_HANDLE h, unsigned int index );

int JB4_INPUTBUFFER_IsEmpty( const JB4_INPUTBUFFER_HANDLE h );
int JB4_INPUTBUFFER_IsFull( const JB4_INPUTBUFFER_HANDLE h );
unsigned int JB4_INPUTBUFFER_Size( const JB4_INPUTBUFFER_HANDLE h );

#endif /* JBM_JB4_INPUTBUFFER_H */
