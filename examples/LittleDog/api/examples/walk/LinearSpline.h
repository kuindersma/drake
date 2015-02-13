#ifndef __LINEAR_SPLINE_H
#define __LINEAR_SPLINE_H

#include <float.h>

//--------------------------------------------------------------------
//
// structure LinearSplinePoint
//
// A (t,value) pair for populating a bdiLinearSplineSet.
//
//--------------------------------------------------------------------
struct LinearSplinePoint 
{
	LinearSplinePoint( float t = 0, float v = 0) { m_t = t; m_value = v; }

	float m_t;
	float m_value; 
};


//--------------------------------------------------------------------
//
// class LinearSplineSet
//
// A set of (t,value) values.
//
//--------------------------------------------------------------------
class LinearSplineSet 
{
	enum SetupErrors 
	{
		LINEAR_SPLINE_OK = 0,
		LINEAR_SPLINE_HAS_LESS_THAN_TWO_POINTS,
		LINEAR_SPLINE_INPUT_OUT_OF_RANGE,
		LINEAR_SPLINE_POORLY_DEFINED
	};

 public:

	// Constructor Args: num of points & the array of points itself
	
	LinearSplineSet(const char* name, int nPoints, LinearSplinePoint** pointArray );
	~LinearSplineSet();
	
	// Accessors for t and value

	void        setPoint( int idx, float t, float value );
	void        setPoint( int idx, const LinearSplinePoint & );
	void        setT( int idx, float value) { m_t[idx] = value; }
	void        setValue( int idx, float value) { m_value[idx] = value; }

	const char* getName( void ) const { return m_name; }
	float       getT( int idx ) { return m_t[idx]; }
	float       getValue( int idx ) { return m_value[idx]; }
	int         getNumPoints( void ) const {return m_numPoints;}
	bool        getIsInBounds(float t);
	int         solve(float t, float* value, float* dvalue); 	// returns SetupErrors

 private:

	int      m_numPoints;      // Number of points in the set
	float*   m_t;              // The array of t values
	float*   m_value;          // The array of values
	char*    m_name;           // optional name for user convenience

	LinearSplineSet( const LinearSplineSet & );
	LinearSplineSet & operator = ( const LinearSplineSet & );
};

#endif  /* __LINEAR_SPLINE_H */
