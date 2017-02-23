/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/

/** \file jbm_jb4_circularbuffer.c circular buffer (FIFO) with fixed capacity */

/* system includes */
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "options.h"
/* local includes */
#include "jbm_jb4_circularbuffer.h"
/* instrumentation */


/** Calculates percentile by selecting greatest elements.
 * This function partial sorts all given elements in the given buffer.
 * @param[in,out] elements ascending sorted buffer of selected greatest elements
 * @param[in,out] size size of elements buffer
 * @param[in]     capacity maximum number of elements to buffer
 * @param[in]     newElement element to insert in buffer if great enough */
static void JB4_CIRCULARBUFFER_calcPercentile( JB4_CIRCULARBUFFER_ELEMENT *elements,
        unsigned int *size, unsigned int capacity, JB4_CIRCULARBUFFER_ELEMENT newElement );

/** circular buffer (FIFO) with fixed capacity */
struct JB4_CIRCULARBUFFER
{
    /** elements of circular buffer */
    JB4_CIRCULARBUFFER_ELEMENT *data;
    /** maximum allowed number of elements plus one free element (to decide between full/empty buffer) */
    unsigned int capacity;
    /** position of next enque operation */
    unsigned int writePos;
    /** position of next deque operation */
    unsigned int readPos;
};


/* Creates a circular buffer (FIFO) */
int JB4_CIRCULARBUFFER_Create( JB4_CIRCULARBUFFER_HANDLE *ph )
{
    JB4_CIRCULARBUFFER_HANDLE h = malloc( sizeof( struct JB4_CIRCULARBUFFER ) );


    h->data     = NULL;
    h->capacity = 0;
    h->writePos = 0;
    h->readPos  = 0;


    *ph = h;


    return 0;
}

/* Destroys the circular buffer (FIFO) */
void JB4_CIRCULARBUFFER_Destroy( JB4_CIRCULARBUFFER_HANDLE *ph )
{
    JB4_CIRCULARBUFFER_HANDLE h;


    if( !ph )
        return;
    h = *ph;
    if( !h )
        return;

    if( h->data )
        free( h->data );
    free( h );
    *ph = NULL;

}

/* Initializes a circular buffer (FIFO) with a fixed maximum allowed number of elements */
int JB4_CIRCULARBUFFER_Init( JB4_CIRCULARBUFFER_HANDLE h, unsigned int capacity )
{


    /* keep one element free to be able to decide between full/empty buffer */
    ++capacity;
    h->data     = malloc( capacity * sizeof( JB4_CIRCULARBUFFER_ELEMENT ) );

    h->capacity = capacity;
    h->writePos = 0;
    h->readPos  = 0;



    return 0;
}

int JB4_CIRCULARBUFFER_Enque( JB4_CIRCULARBUFFER_HANDLE h, JB4_CIRCULARBUFFER_ELEMENT element )
{

    if( JB4_CIRCULARBUFFER_IsFull( h ) )
    {
        return -1;
    }

    h->data[h->writePos] = element;
    ++h->writePos;
    if( h->writePos == h->capacity )
    {
        h->writePos = 0;
    }



    return 0;
}

int JB4_CIRCULARBUFFER_Deque( JB4_CIRCULARBUFFER_HANDLE h, JB4_CIRCULARBUFFER_ELEMENT *pElement )
{

    if( JB4_CIRCULARBUFFER_IsEmpty( h ) )
    {
        return -1;
    }

    *pElement = h->data[h->readPos];
    ++h->readPos;
    if( h->readPos == h->capacity )
    {
        h->readPos = 0;
    }



    return 0;
}

/* Returns the first element. */
JB4_CIRCULARBUFFER_ELEMENT JB4_CIRCULARBUFFER_Front( const JB4_CIRCULARBUFFER_HANDLE h )
{
    JB4_CIRCULARBUFFER_ELEMENT ret;


    ret =  h->data[h->readPos];

    return ret;
}

/* Returns the last element. */
JB4_CIRCULARBUFFER_ELEMENT JB4_CIRCULARBUFFER_Back( const JB4_CIRCULARBUFFER_HANDLE h )
{
    JB4_CIRCULARBUFFER_ELEMENT ret;


    if(h->writePos != 0U)
    {
        ret = h->data[h->writePos - 1];
    }
    else
    {
        ret = h->data[h->capacity - 1];
    }


    return ret;
}

int JB4_CIRCULARBUFFER_IsEmpty( const JB4_CIRCULARBUFFER_HANDLE h )
{
    int ret;


    if(h->readPos == h->writePos)
    {
        ret = 1;
    }
    else
    {
        ret = 0;
    }



    return ret;
}

int JB4_CIRCULARBUFFER_IsFull( const JB4_CIRCULARBUFFER_HANDLE h )
{
    int ret;



    if(((h->writePos + 1) % h->capacity) == h->readPos)
    {
        ret = 1;
    }
    else
    {
        ret = 0;
    }

    return ret;
}

unsigned int JB4_CIRCULARBUFFER_Size( const JB4_CIRCULARBUFFER_HANDLE h )
{
    unsigned int ret;


    if( h->readPos <= h->writePos )
    {
        ret = h->writePos - h->readPos;
    }
    else
    {
        /* wrap around */
        ret = h->writePos + h->capacity - h->readPos;
    }


    return ret;
}

/* Calculates statistics over all elements: min element */
void JB4_CIRCULARBUFFER_Min( const JB4_CIRCULARBUFFER_HANDLE h, JB4_CIRCULARBUFFER_ELEMENT *pMin )
{
    unsigned int i;
    JB4_CIRCULARBUFFER_ELEMENT minEle;


    /* init output variable */
    minEle = h->data[h->readPos];

    if( h->readPos <= h->writePos )
    {
        /* no wrap around */
        /* calc statistics for [readPos;writePos[ */
        for( i = h->readPos; i != h->writePos; ++i )
        {
            if( h->data[i] < minEle )
            {
                minEle = h->data[i];
            }
        }
    }
    else
    {
        /* wrap around */
        /* calc statistics for [readPos;capacity[ */
        for( i = h->readPos; i != h->capacity; ++i )
        {
            if( h->data[i] < minEle )
            {
                minEle = h->data[i];
            }
        }
        /* calc statistics for [0;writePos[ */
        for( i = 0; i != h->writePos; ++i )
        {
            if( h->data[i] < minEle )
            {
                minEle = h->data[i];
            }
        }
    }

    *pMin = minEle;

}

/* Calculates statistics over all elements: max element */
void JB4_CIRCULARBUFFER_Max( const JB4_CIRCULARBUFFER_HANDLE h, JB4_CIRCULARBUFFER_ELEMENT *pMax )
{
    unsigned int i;
    JB4_CIRCULARBUFFER_ELEMENT maxEle;


    /* init output variable */
    maxEle = h->data[h->readPos];
    if( h->readPos <= h->writePos )
    {
        /* no wrap around */
        /* calc statistics for [readPos;writePos[ */
        for( i = h->readPos; i != h->writePos; ++i )
        {
            if( h->data[i] > maxEle )
            {
                maxEle = h->data[i];
            }
        }
    }
    else
    {
        /* wrap around */
        /* calc statistics for [readPos;capacity[ */
        for( i = h->readPos; i != h->capacity; ++i )
        {
            if( h->data[i] > maxEle )
            {
                maxEle = h->data[i];
            }
        }
        /* calc statistics for [0;writePos[ */
        for( i = 0; i != h->writePos; ++i )
        {
            if( h->data[i] > maxEle )
            {
                maxEle = h->data[i];
            }
        }
    }

    *pMax = maxEle;

}

/* Calculates statistics over a considered fraction of all elements: min element and percentile */
void JB4_CIRCULARBUFFER_MinAndPercentile( const JB4_CIRCULARBUFFER_HANDLE h, unsigned int nElementsToIgnore,
        JB4_CIRCULARBUFFER_ELEMENT *pMin, JB4_CIRCULARBUFFER_ELEMENT *pPercentile )
{
    unsigned int i;
    JB4_CIRCULARBUFFER_ELEMENT maxElements[100];
    unsigned int maxElementsSize;
    unsigned int maxElementsCapacity;
    JB4_CIRCULARBUFFER_ELEMENT minEle;


    /* init output variables */
    minEle = h->data[h->readPos];

    /* To calculate the percentile, a number of elements with the highest values are collected in maxElements in
     * ascending sorted order. This array has a size of nElementsToIgnore plus one. This additional element is the
     * lowest of all maxElements, and is called the percentile of all elements. */

    maxElementsSize     = 0;
    maxElementsCapacity = nElementsToIgnore + 1;
    assert( maxElementsCapacity <= sizeof(maxElements) / sizeof(maxElements[0]) );
    if( h->readPos <= h->writePos )
    {
        /* no wrap around */
        /* calc statistics for [readPos;writePos[ */
        for( i = h->readPos; i != h->writePos; ++i )
        {
            if( h->data[i] < minEle )
            {
                minEle = h->data[i];
            }
            JB4_CIRCULARBUFFER_calcPercentile( maxElements, &maxElementsSize, maxElementsCapacity, h->data[i] );
        }
    }
    else
    {
        /* wrap around */
        /* calc statistics for [readPos;capacity[ */
        for( i = h->readPos; i != h->capacity; ++i )
        {
            if( h->data[i] < minEle )
            {
                minEle = h->data[i];
            }
            JB4_CIRCULARBUFFER_calcPercentile( maxElements, &maxElementsSize, maxElementsCapacity, h->data[i] );
        }
        /* calc statistics for [0;writePos[ */
        for( i = 0; i != h->writePos; ++i )
        {
            if( h->data[i] < minEle )
            {
                minEle = h->data[i];
            }
            JB4_CIRCULARBUFFER_calcPercentile( maxElements, &maxElementsSize, maxElementsCapacity, h->data[i] );
        }
    }

    *pPercentile = maxElements[0];
    *pMin = minEle;

}

/* Calculates percentile by selecting greatest elements. */
static void JB4_CIRCULARBUFFER_calcPercentile( JB4_CIRCULARBUFFER_ELEMENT *elements,
        unsigned int *size, unsigned int capacity, JB4_CIRCULARBUFFER_ELEMENT newElement )
{
    unsigned int i;


    /* insert newElement if elements buffer is not yet full */
    if( *size < capacity )
    {
        for( i = 0; i != *size; ++i )
        {
            if( newElement <= elements[i] )
            {
                /* insert newElement at index i */
                memmove( elements + i + 1, elements + i, ( *size - i ) * sizeof( JB4_CIRCULARBUFFER_ELEMENT ) );
                elements[i] = newElement;
                ++*size;
                return;
            }
        }
        /* newElement is maximum, just append it */
        elements[*size] = newElement;
        ++*size;
        return;
    }

    /* check if newElement is too small to be inserted in elements buffer */
    if( newElement <= elements[0] )
    {
        return;
    }

    /* select position to insert newElement to elements */
    for( i = *size - 1; i != 0; --i )
    {
        if( newElement >= elements[i] )
        {
            /* insert newElement at index i */
            memmove( elements, elements + 1, i * sizeof( JB4_CIRCULARBUFFER_ELEMENT ) );
            elements[i] = newElement;
            return;
        }
    }
    /* newElement is just greater than first on in elements buffer */
    elements[0] = newElement;

}

