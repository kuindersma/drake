
#ifndef __SWING_STANCE_GENERATOR_H
#define __SWING_STANCE_GENERATOR_H

#include <bduVec3f.h>
#include "LinearSpline.h"

//
// This class contains persistent spline point memory.  m_spline_points can be
// set and then the m_spline should be used to retrieve the solution.
//

class SwingStanceSpline
{
 public:
	~SwingStanceSpline();
  
	SwingStanceSpline(int num_points, const char *name);
	SwingStanceSpline(const SwingStanceSpline & );
	const char* get_name() const {return m_name;}

 private:
	friend class SwingStanceGenerator;
	class LinearSplineSet* get_spline() { return m_spline; }
  
 private:
  
	// WARNING - this class has a clone constructor.  Addition or removal
	// of class members will require updates to the clone constructor.
  
	struct LinearSplinePoint**  m_spline_points;
	class LinearSplineSet*      m_spline;
	const char*                 m_name;
};


/******************************************************************************/

// From the user's point of view, this class consist of two end-point
// matched splines (s1 and s2).  The first spline (s1 == the swing spline) 
// has n1 points and the second spline (s2 == the stance spline) has n2 
// points.  Because endpoints match, the point at s1(0) == s2(n2-1) and 
// the point at s1(n1-1) == s2(0).  time runs from s1 t(0) through 
// s1 t(n1-1) on s1.  time runs from s1 t(n1-1) (== s2 t(0)) through 
// s2 t(n2-1) on s2.  Users should call the various set_swing_spline_*() and 
// set_stance_spline_*() functions.  When done, call solve() as often as needed.
//
// y_mirror of 1 implies that all m_*_spline_y values will negated.
//
// t_offset (whose magnitude must be within +/- a cycle duration) is 
// internally added to t in get_xyz() to produce the returned xyz.
//
// the xyz_offset is also internally added in get_xyz() prior to returning.
//
// is_stance indicates whether or not the last call to solve(t, ...) used a t
// (internally with the t_offset added to it) that fell in the stance time range.
//

class SwingStanceGenerator
{
 public:

	SwingStanceGenerator(int num_Swing_points, int num_stance_points); // (n1, n2)
	SwingStanceGenerator( const SwingStanceGenerator & );
	~SwingStanceGenerator();

	// GENERAL
	void   set_y_mirrored(bool enable_y_mirror)    {m_y_mirrored = enable_y_mirror;}
	void   set_t_offset(float t_offset)           {m_t_offset = t_offset;}
	void   set_xyz_offset(const bduVec3f & xyz_offset);
	
	int    solve(float t, float* x, float* y, float* z, float* dx, float* dy, float* dz);
	int    solve(float t, bduVec3f* xyz, bduVec3f* dxdydz);
	
	int    get_num_swing_points( void );
	int    get_num_stance_points( void );
	float  get_cycle_duration( void );
	float  get_t_offset( void ) { return this->m_t_offset; }
	bool   is_y_mirrored( void )                         {return m_y_mirrored;}
	bool   in_stance( void )                             {return m_in_stance;}

	const bduVec3f & get_xyz_offset( void )              {return m_xyz_offset;}

	// SWING
	int    set_swing_spline_point(int point_index, const bduVec3f & pos);
	int    set_swing_spline_time(int point_index, float t);

	int    get_swing_spline_point(int point_index, bduVec3f* pos);
	int    get_swing_spline_time(int point_index, float* t);

	// STANCE
	int    set_stance_spline_point(int point_index, const bduVec3f & pos);
	int    set_stance_spline_time(int point_index, float t);

	int    get_stance_spline_point(int point_index, bduVec3f* pos);
	int    get_stance_spline_time(int point_index, float* t);

 private:

	void   set_spline_point_v(int point_index, float value, SwingStanceSpline* spline);
	void   set_spline_point_t(int point_index, float t, SwingStanceSpline* spline);

	void   get_spline_point_v(int point_index, float* value, SwingStanceSpline* spline);
	void   get_spline_point_t(int point_index, float* t, SwingStanceSpline* spline);

	void   update_offset_matrix( void );

private:

	// WARNING - this class has a clone constructor.  Addition or removal
	// of class members will require updates to the clone constructor.

	bool                m_y_mirrored;
	bool                m_in_stance;
	float               m_t_offset;        // time offset used internally
	bduVec3f            m_xyz_offset;      // position offset used internally
	SwingStanceSpline*  m_swing_spline_x;
	SwingStanceSpline*  m_swing_spline_y;
	SwingStanceSpline*  m_swing_spline_z;
	SwingStanceSpline*  m_stance_spline_x;
	SwingStanceSpline*  m_stance_spline_y;
	SwingStanceSpline*  m_stance_spline_z;
};


#endif /* __SWING_STANCE_GENERATOR_H */

