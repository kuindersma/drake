
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include <bduLog.h>

#include "SwingStanceGenerator.h"

/**************************/
/* SwingStanceSpline code */
/**************************/

/******************************************************************************/
SwingStanceSpline::SwingStanceSpline(int num_points, const char *name)
{
	m_spline_points = new LinearSplinePoint*[num_points];
	m_name = strdup(name );	

	int point_index;
	for (point_index = 0; point_index < num_points; point_index++)
		m_spline_points[point_index] = new LinearSplinePoint(0,0);
	
	m_spline = new LinearSplineSet(name, num_points,
		(LinearSplinePoint**)m_spline_points );
}


/******************************************************************************/
SwingStanceSpline::SwingStanceSpline( const SwingStanceSpline & rhs )
{
	int num_points = rhs.m_spline->getNumPoints();
	m_spline_points = new LinearSplinePoint*[num_points];
	m_name = strdup( rhs.m_spline->getName() );

	int point_index;
	for (point_index = 0; point_index < num_points; point_index++)
		m_spline_points[point_index] = new LinearSplinePoint(
			rhs.m_spline->getT(point_index), rhs.m_spline->getValue(point_index));
		
	m_spline = new LinearSplineSet( (char *)rhs.m_spline->getName(), num_points, 
		(LinearSplinePoint**)m_spline_points );
}


/******************************************************************************/
SwingStanceSpline:: ~SwingStanceSpline() 
{
	int point_index;
	
	for (point_index = 0; point_index < m_spline->getNumPoints(); point_index++ )
		delete m_spline_points[point_index];
	
	delete [] m_spline_points;		
	m_spline_points = NULL;

	delete m_spline;
	m_spline = NULL;

	free((char*)m_name);
	m_name = NULL;
}


/*****************************/
/* SwingStanceGenerator code */
/*****************************/


/******************************************************************************/
SwingStanceGenerator::SwingStanceGenerator(int num_swing_points, int num_stance_points)
{
	// (n1, n2)
	m_in_stance         = false;
	m_y_mirrored        = false;
	m_t_offset          = 0;
	m_xyz_offset        = bduVec3f(0,0,0);
	m_swing_spline_x    = new SwingStanceSpline(num_swing_points,  "swing_x spline");
	m_swing_spline_y    = new SwingStanceSpline(num_swing_points,  "swing_y spline");
	m_swing_spline_z    = new SwingStanceSpline(num_swing_points,  "swing_z spline");
	m_stance_spline_x   = new SwingStanceSpline(num_stance_points, "stance_x spline");
	m_stance_spline_y   = new SwingStanceSpline(num_stance_points, "stance_y spline");
	m_stance_spline_z   = new SwingStanceSpline(num_stance_points, "stance_z spline");
}


/******************************************************************************/
SwingStanceGenerator::SwingStanceGenerator(const SwingStanceGenerator & rhs )
{
	m_in_stance         = rhs.m_in_stance;
	m_y_mirrored        = rhs.m_y_mirrored;
	m_t_offset          = rhs.m_t_offset;
	m_xyz_offset        = rhs.m_xyz_offset;
	m_swing_spline_x    = new SwingStanceSpline( *(rhs.m_swing_spline_x) );
	m_swing_spline_y    = new SwingStanceSpline( *(rhs.m_swing_spline_y) );
	m_swing_spline_z    = new SwingStanceSpline( *(rhs.m_swing_spline_z) );
	m_stance_spline_x   = new SwingStanceSpline( *(rhs.m_stance_spline_x) );
	m_stance_spline_y   = new SwingStanceSpline( *(rhs.m_stance_spline_y) );
	m_stance_spline_z   = new SwingStanceSpline( *(rhs.m_stance_spline_z) );
}


/******************************************************************************/
SwingStanceGenerator::~SwingStanceGenerator() 
{
	delete m_swing_spline_x;
	delete m_swing_spline_y;
	delete m_swing_spline_z;
	delete m_stance_spline_x;
	delete m_stance_spline_y;
	delete m_stance_spline_z;
	
	m_y_mirrored       = false;
	m_t_offset         = 0;
	m_xyz_offset       = bduVec3f(0,0,0);
	m_swing_spline_x   = m_swing_spline_y  = m_swing_spline_z  = NULL;
	m_stance_spline_x  = m_stance_spline_y = m_stance_spline_z = NULL;
}


/******************************************************************************/
int 
SwingStanceGenerator::get_num_swing_points()  
{ 
	return m_swing_spline_x->get_spline()->getNumPoints(); 
}


/******************************************************************************/
int 
SwingStanceGenerator::get_num_stance_points() 
{ 
	return m_stance_spline_x->get_spline()->getNumPoints(); 
}


/******************************************************************************/
void
SwingStanceGenerator::set_spline_point_v(int point_index, float value, 
	SwingStanceSpline* spline)
{		
	spline->get_spline()->setValue(point_index,value);
}


/******************************************************************************/
void
SwingStanceGenerator::set_spline_point_t(int point_index, float t, 
	SwingStanceSpline* spline)
{
	spline->get_spline()->setT(point_index,t);
}


/******************************************************************************/
void
SwingStanceGenerator::get_spline_point_v(int point_index, float* value, 
	SwingStanceSpline* spline)
{
	*value = spline->get_spline()->getValue(point_index);
}


/******************************************************************************/
void
SwingStanceGenerator::get_spline_point_t(int point_index, float* t, 
	SwingStanceSpline* spline)
{
	*t = spline->get_spline()->getT(point_index);
}


/******************************************************************************/
int 
SwingStanceGenerator::set_swing_spline_point(int point_index, const bduVec3f & pos)
{
	if ( (point_index < 0) || 
		(point_index >= m_swing_spline_x->get_spline()->getNumPoints()) )
	{
		bdu_log_printf( BDU_LOG_LEVEL_WARN, 
			"SwingStanceGenerator::set_swing_spline_point(): index %d out of range.\n", 
			point_index);
		return -1;
	}

	bdu_log_printf( BDU_LOG_LEVEL_INFO, 
		"adding swing point[%d] to generator: (%f,%f,%f)\n", 
		point_index, pos.n[0], pos.n[1], pos.n[2]);
	
	set_spline_point_v(point_index, pos.n[0], m_swing_spline_x);
	set_spline_point_v(point_index, pos.n[1], m_swing_spline_y);
	set_spline_point_v(point_index, pos.n[2], m_swing_spline_z);
	
	return 0;
}


/******************************************************************************/
int 
SwingStanceGenerator::get_swing_spline_point(int point_index, bduVec3f* pos)
{
	if ( (point_index < 0) || 
		(point_index >= m_swing_spline_x->get_spline()->getNumPoints()) )
	{
		bdu_log_printf( BDU_LOG_LEVEL_WARN, 	
			"SwingStanceGenerator::get_swing_spline_point(): index %d out of range.\n", 
			point_index);
		return -1;
	}
	
	get_spline_point_v(point_index, &((*pos).n[0]), m_swing_spline_x);
	get_spline_point_v(point_index, &((*pos).n[1]), m_swing_spline_y);
	get_spline_point_v(point_index, &((*pos).n[2]), m_swing_spline_z);

	return 0;
}


/******************************************************************************/
int 
SwingStanceGenerator::set_stance_spline_point(int point_index, const bduVec3f & pos)
{
	if ( (point_index < 0) || 
		(point_index >= m_stance_spline_x->get_spline()->getNumPoints()) )
	{
		bdu_log_printf( BDU_LOG_LEVEL_WARN, 
			"SwingStanceGenerator::set_stance_spline_point(): index %d out of range.\n", 
			point_index);
		return -1;
	}
	
	bdu_log_printf( BDU_LOG_LEVEL_INFO, 
		"adding stance point[%d] to generator: (%f,%f,%f)\n", 
		point_index, pos.n[0], pos.n[1], pos.n[2]);

	set_spline_point_v(point_index, pos.n[0], m_stance_spline_x);
	set_spline_point_v(point_index, pos.n[1], m_stance_spline_y);
	set_spline_point_v(point_index, pos.n[2], m_stance_spline_z);
	
	return 0;
}

    
/******************************************************************************/
int 
SwingStanceGenerator::get_stance_spline_point(int point_index, bduVec3f* pos)
{
	if ( (point_index < 0) || 
		(point_index >= m_stance_spline_x->get_spline()->getNumPoints()) )
	{
		bdu_log_printf( BDU_LOG_LEVEL_WARN, 
			"SwingStanceGenerator::get_stance_spline_point(): index %d out of range.\n", 
			point_index);
		return -1;
	}

	get_spline_point_v(point_index, &((*pos).n[0]), m_stance_spline_x);
	get_spline_point_v(point_index, &((*pos).n[1]), m_stance_spline_y);
	get_spline_point_v(point_index, &((*pos).n[2]), m_stance_spline_z);

	return 0;
}


/******************************************************************************/
int 
SwingStanceGenerator::set_swing_spline_time(int point_index, float t)
{
	if ( (point_index < 0) || 
		(point_index >= m_swing_spline_x->get_spline()->getNumPoints()) )
	{
		bdu_log_printf( BDU_LOG_LEVEL_WARN, 
			"SwingStanceGenerator::set_swing_spline_time(): index %d out of range.\n", 
			point_index);
		return -1;
	}	

	bdu_log_printf( BDU_LOG_LEVEL_WARN, 
		"setting swing time for point[%d] to generator: (%f)\n", point_index, t );

	set_spline_point_t(point_index, t, m_swing_spline_x);
	set_spline_point_t(point_index, t, m_swing_spline_y);
	set_spline_point_t(point_index, t, m_swing_spline_z);

	return 0;
}


/******************************************************************************/
int 
SwingStanceGenerator::get_swing_spline_time(int point_index, float* t)
{
	if ( (point_index < 0) || 
		(point_index >= m_swing_spline_x->get_spline()->getNumPoints()) )
	{
		bdu_log_printf( BDU_LOG_LEVEL_WARN, 
			"SwingStanceGenerator::get_swing_spline_time(): index %d out of range.\n", 
			point_index);
		return -1;
	}	

	get_spline_point_t(point_index, t, m_swing_spline_x);

	return 0;
}


/******************************************************************************/
int 
SwingStanceGenerator::set_stance_spline_time(int point_index, float t)
{
	if ( (point_index < 0) || 
		(point_index >= m_stance_spline_x->get_spline()->getNumPoints()) )
	{
		bdu_log_printf( BDU_LOG_LEVEL_WARN, 
			"SwingStanceGenerator::set_stance_spline_time(): index %d out of range.\n", 
			point_index);
		return -1;
	}	

	bdu_log_printf( BDU_LOG_LEVEL_INFO, 
		"setting stance time for point[%d] to generator: (%f)\n", point_index, t );

	set_spline_point_t(point_index, t, m_stance_spline_x);
	set_spline_point_t(point_index, t, m_stance_spline_y);
	set_spline_point_t(point_index, t, m_stance_spline_z);

	return 0;
}


/******************************************************************************/
int 
SwingStanceGenerator::get_stance_spline_time(int point_index, float* t)
{
	if ( (point_index < 0) || 
		(point_index >= m_stance_spline_x->get_spline()->getNumPoints()) )
	{
		bdu_log_printf( BDU_LOG_LEVEL_WARN, 
			"SwingStanceGenerator::get_stance_spline_time(): index %d out of range.\n", 
			point_index);
		return -1;
	}	
	
	get_spline_point_t(point_index, t, m_stance_spline_x);

	return 0;
}


/******************************************************************************/
void   
SwingStanceGenerator::set_xyz_offset(const bduVec3f & xyz_offset)     
{
	m_xyz_offset = xyz_offset;
}


/******************************************************************************/
int
SwingStanceGenerator::solve(float t, bduVec3f* xyz, bduVec3f* dxdydz)
{
	float x, y, z, dx, dy, dz;
	int rval = solve(t, &x, &y, &z, &dx, &dy, &dz);

	(*xyz).n[ 0 ] = x;
	(*xyz).n[ 1 ] = y;
	(*xyz).n[ 2 ] = z;
	(*dxdydz).n[ 0 ] = dx;
	(*dxdydz).n[ 1 ] = dy;
	(*dxdydz).n[ 2 ] = dz;

	return rval;	
}


/******************************************************************************/
int 
SwingStanceGenerator::solve(float t, float* x, float* y, float* z, 
	float *dx, float* dy, float* dz)
{
	// Input t must fall in the range of s1 t(0) and s2 t(n2-1)
	// Note also that the point s2(n2-1) must equal s1(0).
	// Returns 0 on success, -1 on failure

	int stance_points = m_stance_spline_x->get_spline()->getNumPoints();

	if ( (m_swing_spline_x->get_spline()->getT(0) > t) || 
		(t > m_stance_spline_x->get_spline()->getT(stance_points-1) ) )
	{
		bdu_log_printf( BDU_LOG_LEVEL_WARN,
			"SwingStanceGenerator::get_xyz(): t (%f) out of range (%f,%f)\n", t, 
			m_swing_spline_x->get_spline()->getT(0),  
			m_stance_spline_x->get_spline()->getT(stance_points-1));

		return -1;
	}

	// perform t-offseting (wrap t and t_offset)
	float cycle_duration = m_stance_spline_x->get_spline()->getT(stance_points-1) -
		m_swing_spline_x->get_spline()->getT(0);
	float used_t = t + m_t_offset;

	if ( used_t > m_stance_spline_x->get_spline()->getT(stance_points-1) )
		used_t -= cycle_duration;
	else if ( used_t < m_swing_spline_x->get_spline()->getT(0) )
		used_t += cycle_duration;

	// solve for x,y,z given used_t
	int error_code_x, error_code_y, error_code_z;
	bduVec3f result, dresult;

	// don't log here since LinearSplineSet.cpp seems to do so for us.
	if ( used_t < m_stance_spline_x->get_spline()->getT(0) )
	{
		m_in_stance = false;
		error_code_x = m_swing_spline_x->get_spline()->solve(used_t, &result.n[0], &dresult.n[0]);
		error_code_y = m_swing_spline_y->get_spline()->solve(used_t, &result.n[1], &dresult.n[1]);
	  error_code_z = m_swing_spline_z->get_spline()->solve(used_t, &result.n[2], &dresult.n[2]);
	}
	else
	{
		m_in_stance = true;
		error_code_x = m_stance_spline_x->get_spline()->solve(used_t, &result.n[0], &dresult.n[0]);
		error_code_y = m_stance_spline_y->get_spline()->solve(used_t, &result.n[1], &dresult.n[1]);
		error_code_z = m_stance_spline_z->get_spline()->solve(used_t, &result.n[2], &dresult.n[2]);
	}

	// don't log here since LinearSplineSet.cpp already does so for us
	if ( (error_code_x != 0) || (error_code_y != 0) || (error_code_z != 0) )
		return -1;

	// y-mirror if needed
	if ( m_y_mirrored )
		result.n[1] *= -1;

	// offset xyz as needed
	*x = result.n[0] + m_xyz_offset.n[0];
	*y = result.n[1] + m_xyz_offset.n[1];
	*z = result.n[2] + m_xyz_offset.n[2];

	// there are no velocity offsets, so copy directly
	*dx = dresult.n[0];
	*dy = dresult.n[1];
	*dz = dresult.n[2];

	return 0;
}


/******************************************************************************/
float
SwingStanceGenerator::get_cycle_duration()
{
	// returns duration of entire swing-stance

	int num_stance_points = m_stance_spline_x->get_spline()->getNumPoints();

	return m_stance_spline_x->get_spline()->getT(num_stance_points-1) -
		m_swing_spline_x->get_spline()->getT(0);
}

