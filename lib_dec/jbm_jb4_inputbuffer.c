/*====================================================================================
    EVS Codec 3GPP TS26.443 Aug 30, 2016. Version 12.7.0 / 13.3.0
  ====================================================================================*/


/** \file jbm_jb4_inputbuffer.c RTP input buffer with fixed capacity. */

/* system includes */
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "options.h"
/* instrumentation */
/* local includes */
#include "jbm_jb4_inputbuffer.h"


/** input buffer with fixed capacity */
struct JB4_INPUTBUFFER
{
    /** elements of input buffer */
    JB4_INPUTBUFFER_ELEMENT *data;
    /** maximum allowed number of elements plus one free element (to decide between full/empty buffer) */
    unsigned int capacity;
    /** position of next enque operation */
    unsigned int writePos;
    /** position of next deque operation */
    unsigned int readPos;
    /** function to compare two elements */
    int (*compareFunction)( const JB4_INPUTBUFFER_ELEMENT first, const JB4_INPUTBUFFER_ELEMENT second,
                            bool_t *replaceWithNewElementIfEqual );
};


/* Creates a input buffer */
int JB4_INPUTBUFFER_Create( JB4_INPUTBUFFER_HANDLE *ph )
{
    JB4_INPUTBUFFER_HANDLE h = malloc( sizeof( struct JB4_INPUTBUFFER ) );


    h->data            = NULL;
    h->capacity        = 0;
    h->writePos        = 0;
    h->readPos         = 0;
    h->compareFunction = NULL;

    *ph = h;

    return 0;
}

/* Destroys the input buffer */
void JB4_INPUTBUFFER_Destroy( JB4_INPUTBUFFER_HANDLE *ph )
{
    JB4_INPUTBUFFER_HANDLE h;

    if( !ph )
    {
        return;
    }
    h = *ph;
    if( !h )
    {
        return;
    }
    if( h->data )
        free( h->data );
    free( h );
    *ph = NULL;

}

/* Initializes a input buffer with a fixed maximum allowed number of elements */
int JB4_INPUTBUFFER_Init( JB4_INPUTBUFFER_HANDLE h, unsigned int capacity,
                          int (*compareFunction)( const JB4_INPUTBUFFER_ELEMENT first, const JB4_INPUTBUFFER_ELEMENT second,
                                  bool_t *replaceWithNewElementIfEqual ) )
{

    /* keep one element free to be able to decide between full/empty buffer */
    ++capacity;
    h->data            = malloc( capacity * sizeof( JB4_INPUTBUFFER_ELEMENT ) );
    h->capacity        = capacity;
    h->writePos        = 0;
    h->readPos         = 0;
    h->compareFunction = compareFunction;

    return 0;
}

int JB4_INPUTBUFFER_Enque( JB4_INPUTBUFFER_HANDLE h, JB4_INPUTBUFFER_ELEMENT element,
                           JB4_INPUTBUFFER_ELEMENT *replacedElement )
{
    unsigned int size;
    int low, high, middle, diff;
    unsigned int insertPos;
    unsigned int canMoveRight;
    unsigned int canMoveLeft;
    bool_t replace;
    *replacedElement = NULL;

    size = JB4_INPUTBUFFER_Size( h );
    if(size >= h->capacity - 1)
    {
        return -1;
    }

    /* appending the first element is straight forward */
    if( size == 0U )
    {
        h->data[h->writePos] = element;
        ++h->writePos;
        if( h->writePos == h->capacity )
        {
            h->writePos = 0;
        }
        return 0;
    }

    /* there's a high probability that the new element can be appended at the back */
    if( h->compareFunction( element, JB4_INPUTBUFFER_Back( h ), &replace ) > 0 )
    {
        h->data[h->writePos] = element;
        ++h->writePos;
        if( h->writePos == h->capacity )
        {
            h->writePos = 0;
        }
        return 0;
    }

    /* out of order: use binary search to get the position to insert */
    low  = 0;
    high = size - 1;
    while( low <= high )
    {
        middle = low + ( high - low ) / 2;
        diff = h->compareFunction( element, JB4_INPUTBUFFER_Element( h, middle ), &replace );
        if( diff < 0 )
        {
            high = middle - 1;
        }
        else if( diff > 0 )
        {
            low = middle + 1;
        }
        else   /* an element with same index is already stored */
        {
            if(replace != 0)
            {
                *replacedElement = h->data[( h->readPos + middle ) % h->capacity];
                h->data[( h->readPos + middle ) % h->capacity] = element;
                return 0;
            }
            return 1;
        }
    }

    assert( h->compareFunction( element, JB4_INPUTBUFFER_Element( h, low ), &replace ) != 0 );
    if( low > 0 )
        assert( h->compareFunction( element, JB4_INPUTBUFFER_Element( h, low - 1 ), &replace ) > 0 );
    assert( h->compareFunction( element, JB4_INPUTBUFFER_Element( h, low ), &replace ) < 0 );
    if( (unsigned int)(low + 1) < size )
        assert( h->compareFunction( element, JB4_INPUTBUFFER_Element( h, low + 1 ), &replace ) < 0 );

    insertPos = ( h->readPos + low ) % h->capacity;
    if( h->readPos < h->writePos )
    {
        canMoveRight = 1;
        canMoveLeft  = h->readPos > 0;
    }
    else
    {
        canMoveRight = insertPos < h->writePos;
        canMoveLeft  = insertPos > h->writePos;
    }

    assert( canMoveRight != 0 || canMoveLeft != 0 );
    (void)canMoveLeft;
    if( canMoveRight )
    {
        /* move higher elements to the right and insert at insertPos */
        memmove( h->data + insertPos + 1, h->data + insertPos,
                 ( h->writePos - insertPos ) * sizeof( JB4_INPUTBUFFER_ELEMENT ) );
        h->data[insertPos] = element;
        ++h->writePos;
        if( h->writePos == h->capacity )
        {
            h->writePos = 0;
        }
    }
    else
    {
        /* move lower elements to the left and insert before insertPos */
        memmove( h->data + h->readPos - 1, h->data + h->readPos,
                 low * sizeof( JB4_INPUTBUFFER_ELEMENT ) );
        h->data[insertPos-1] = element;
        --h->readPos;
        assert( h->readPos < 99999 );
    }

    return 0;
}

int JB4_INPUTBUFFER_Deque( JB4_INPUTBUFFER_HANDLE h, JB4_INPUTBUFFER_ELEMENT *pElement )
{
    if( JB4_INPUTBUFFER_IsEmpty( h ) )
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
JB4_INPUTBUFFER_ELEMENT JB4_INPUTBUFFER_Front( const JB4_INPUTBUFFER_HANDLE h )
{
    JB4_INPUTBUFFER_ELEMENT ret;


    ret = h->data[h->readPos];

    return ret;
}

/* Returns the last element. */
JB4_INPUTBUFFER_ELEMENT JB4_INPUTBUFFER_Back( const JB4_INPUTBUFFER_HANDLE h )
{
    JB4_INPUTBUFFER_ELEMENT ret;

    if( h->writePos != 0U )
    {
        ret = h->data[h->writePos - 1];
    }
    else
    {
        ret = h->data[h->capacity - 1];
    }

    return ret;
}

/* Returns the element with the given index (0 means front element). */
JB4_INPUTBUFFER_ELEMENT JB4_INPUTBUFFER_Element( const JB4_INPUTBUFFER_HANDLE h, unsigned int index )
{
    JB4_INPUTBUFFER_ELEMENT ret;


    /* return h->data[(h->readPos + index) % h->capacity] without error handling */
    if( h->readPos + index < h->capacity )
    {
        ret = h->data[h->readPos + index];
    }
    else
    {
        /* wrap around */
        ret = h->data[h->readPos + index - h->capacity];
    }

    return ret;
}

int JB4_INPUTBUFFER_IsEmpty( const JB4_INPUTBUFFER_HANDLE h )
{
    int ret;


    ret = h->readPos == h->writePos;

    return ret;
}

int JB4_INPUTBUFFER_IsFull( const JB4_INPUTBUFFER_HANDLE h )
{
    int ret;
    ret = 0;
    if( JB4_INPUTBUFFER_Size( h ) == h->capacity - 1 )
    {
        ret = 1;
    }
    return ret;
}

unsigned int JB4_INPUTBUFFER_Size( const JB4_INPUTBUFFER_HANDLE h )
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

