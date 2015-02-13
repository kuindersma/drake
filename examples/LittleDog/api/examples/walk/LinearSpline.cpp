#include <math.h>
#include <stddef.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include <bduLog.h>

#include "LinearSpline.h"

#define LINEAR_SPLINE_IS_EQUAL(a,b,tol)  ((((a)+(tol))>=(b)) && (((a)-(tol))<=(b)))

//--------------------------------------------------------------------
LinearSplineSet::LinearSplineSet( const char* name,
	int nPoints, LinearSplinePoint** pointArray ) :
	m_numPoints( nPoints ),
	m_t( NULL ),
	m_value( NULL ),
	m_name( NULL )
{
	m_name = strdup( name );

	// Assign space to copy point data
	m_t      = new float[ m_numPoints ];
	m_value  = new float[ m_numPoints ];

	// Process points and bail out if we get unsorted pairs 
	// Copy data to the x and y arrays
	int i;
	for ( i = 0; i < m_numPoints; i++ ) 
	{
		if ( (i < m_numPoints - 1) && (pointArray[i]->m_t > pointArray[i+1]->m_t) ) 
		{
			bdu_log_printf( BDU_LOG_LEVEL_WARN, "LinearSplinePoint constructor: "
				"Unsorted x values %f, %f from passed points %d, %d", 
				pointArray[i]->m_t, pointArray[i+1]->m_t, i, i+1); 
			return;
		}

		m_t[i] = pointArray[i]->m_t;
		m_value[i] = pointArray[i]->m_value;
	}
}

//--------------------------------------------------------------------
LinearSplineSet::~LinearSplineSet()
{
	delete [] m_t;
	m_t = NULL;

	delete [] m_value;
	m_value = NULL;

	free(m_name);
	m_name = NULL;
}

//--------------------------------------------------------------------
bool
LinearSplineSet::getIsInBounds(float t)
{
	return ((!m_numPoints) || (t < m_t[0]) || (t > m_t[m_numPoints-1])) ? 
		false : true; 
}

//--------------------------------------------------------------------
void
LinearSplineSet::setPoint( int idx, const LinearSplinePoint &p )
{
	setT( idx, p.m_t );
	setValue( idx, p.m_value );
}

//--------------------------------------------------------------------
void
LinearSplineSet::setPoint( int idx, float x, float y )
{
	setT( idx, x );
	setValue( idx, y );
}

//--------------------------------------------------------------------
int
LinearSplineSet::solve(float t, float* value, float* dvalue)
{
	if ( getIsInBounds(t) == false )
	{
		bdu_log_printf( BDU_LOG_LEVEL_WARN, 
			"LinearSplineSet::solve, given t (%f) is for spline %s is "
			"out of range!\n",  t,  m_name);
		return LINEAR_SPLINE_INPUT_OUT_OF_RANGE;
	}

	// find the segment that corresponds to t
	int seg;
	for ( seg = 1; seg < m_numPoints; seg++ )
		if ( t <= m_t[seg] )
			break;

	seg--;

	// solve the value at time t for the found segment
	float t0 = m_t[seg];
	float t1 = m_t[seg + 1];
	float v0 = m_value[seg];
	float v1 = m_value[seg + 1];

	if ( LINEAR_SPLINE_IS_EQUAL( t1 - t0, 0, 1.0e-8 ) )
	{
		bdu_log_printf(BDU_LOG_LEVEL_WARN, 
			"LinearSplineSet::solve, spline points for spline %s, seg %d"
			" are poorly defined (t1 == t0) (%f)!\n",m_name, seg,  t0);
		return LINEAR_SPLINE_POORLY_DEFINED;
	}

	if ( LINEAR_SPLINE_IS_EQUAL( v1 - v0, 0, 1.0e-8 ) )
		*value = v0;
	else
		*value = v0 + (((t - t0)/(t1 - t0)) * (v1 - v0));

	// compute the velocity
	*dvalue = (v1 - v0) / (t1 - t0);

	return LINEAR_SPLINE_OK;
}


//--------------------------------------------------------------------
