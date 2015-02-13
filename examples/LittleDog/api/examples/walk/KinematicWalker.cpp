
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <memory.h>
#include <math.h>

#include <bduLog.h>
#include <bduVec3f.h>
#include <LittleDogDataSetFileStreamWriter.h>

#include "SwingStanceGenerator.h"
#include "KinematicWalker.h"
#include "MyLittleDog.h"

//--------------------------------------------------------------------

struct SplineDefaults
{
  float swing_p0[3];
  float swing_p0_t;
  float swing_p1[3];
  float swing_p1_t;
  float swing_p2[3];
  float swing_p2_t;
  float swing_p3[3];
  float swing_p3_t;
  float swing_p4[3];
  float swing_p4_t;
  float stance_p0_t;
  float stance_p1_t;
  float t_offset[LittleDog::NUM_LEGS];
  float xyz_offset[3];
  float t_scale;
};

//--------------------------------------------------------------------

#define MY_MAX_LOGGED     1000
#define MY_VARS_PER_LEG   8
#define MY_NON_LEG_VARS   3
#define MY_VAR_COUNT      MY_NON_LEG_VARS + (LittleDog::NUM_LEGS * MY_VARS_PER_LEG)
#define MAX_X_SCALE_WARP  0.25f
#define MAX_Y_SWAY_WARP   0.03   // meters

//--------------------------------------------------------------------

struct KinematicWalkerPrivate
{
	KinematicWalkerPrivate():
		m_spline_t( 0 ),
		m_prev_t( 0 ),
		m_turn_amount( 0 ),
		m_dataset( NULL ),
		m_log_count( 0 )
	{
		for ( int l = 0; l < LittleDog::NUM_LEGS; l++ )
			m_splines[ l ] = NULL;

		for ( int v = 0; v < MY_VAR_COUNT; v++ )
			m_vars[ v ] = 0;
	}

	float                               m_spline_t;
	float                               m_prev_t;
	float                               m_turn_amount;
	SwingStanceGenerator *              m_splines[ LittleDog::NUM_LEGS ];
	struct SplineDefaults               m_spline_defaults;
	LittleDogDataSetFileStreamWriter*   m_dataset;
	bduDataSetVarID                     m_vars[MY_VAR_COUNT];
	int                                 m_log_count;
	bduVec3f                            m_pos_foot_center_rt_body[ LittleDog::NUM_LEGS ];
};

//--------------------------------------------------------------------

KinematicWalker::~KinematicWalker()
{
	for ( int l = 0; l < LittleDog::NUM_LEGS; ++l )
		delete p->m_splines[ l ];

	if ( p->m_dataset )
	{
		bdu_log_printf(BDU_LOG_LEVEL_INFO, 
			"My data logged to my_dataset.data\n");
		p->m_dataset->save();
		delete p->m_dataset;
	}		

	delete p;
	p = NULL;
}


//--------------------------------------------------------------------

KinematicWalker::KinematicWalker( MyLittleDog & l ):
	p( new KinematicWalkerPrivate ),
	dog( l )
{
	assert( p );
}


//--------------------------------------------------------------------
void
KinematicWalker::initialize( void )
{
	bdu_log_printf( BDU_LOG_LEVEL_INFO, "KinematicWalker::initialize\n" );
	p->m_log_count = 0;	
	
	// learn about the robot
	initialize_foot_centers();

	// spline stuff
	initialize_params();

	int leg, joint;
	for ( leg = 0; leg < LittleDog::NUM_LEGS; ++leg )
		initialize_spline_for_leg( (LittleDog::LegIndex)leg );

	dog.getDataTimestamp(&p->m_prev_t);
	p->m_spline_t = p->m_spline_defaults.swing_p0_t;
	
	// robot stuff.
	for ( leg = 0; leg < LittleDog::NUM_LEGS; ++leg )
	{
		for ( joint = 0; joint < LittleDog::NUM_JOINTS; ++joint )
		{
			float k[] = { 7.5, 7.5, 6.7 };
			float b[] = { 0.32, 0.32, 0.32 };
			
			dog.setJointServoType( (LittleDog::LegIndex)leg, 
				(LittleDog::JointIndex)joint, LittleDog::SERVO_PD );
			dog.setJointPDServoGains( (LittleDog::LegIndex)leg, 
				(LittleDog::JointIndex)joint, k[ joint ], b[ joint ] );
		}

		dog.setLegCommandMode((LittleDog::LegIndex)leg, LittleDog::LEG_MODE_CARTESIAN );
	}

	// dataset stuff
	initialize_dataset();

	// change our internal state
	m_state = &KinematicWalker::stateInit;	

	bdu_log_printf( BDU_LOG_LEVEL_INFO, "kinematicWalker::initialize done.\n");
}

//--------------------------------------------------------------------
void
KinematicWalker::uninitialize( void )
{
	m_state = (void ( KinematicWalker:: * )() )NULL;

	for ( int leg = 0; leg < LittleDog::NUM_LEGS; ++leg )
		for ( int joint = 0; joint < LittleDog::NUM_JOINTS; ++joint )
			dog.setJointServoType( (LittleDog::LegIndex)leg, 
				(LittleDog::JointIndex)joint, LittleDog::SERVO_LIMP );

	if ( p->m_dataset )
	{
		p->m_dataset->save();
		delete p->m_dataset;
		p->m_dataset = NULL;
	}

	bdu_log_printf( BDU_LOG_LEVEL_INFO, "KinematicWalker::deactivate done.\n");
}

//--------------------------------------------------------------------
// one-time setup stuff
//--------------------------------------------------------------------
void 
KinematicWalker::initialize_foot_centers( void )
{
	// Computes the zero configuration location of each foot.
	// This is the xyz position that the foot would be relative to the
	// robot's body when all leg joints are set to zero.

	for ( int leg = 0; leg < LittleDog::NUM_LEGS; leg++ )
	{
		bduVec3f pos_hip_rx_rt_body;
		bduVec3f pos_hip_ry_rt_hip_rx;
		bduVec3f pos_knee_ry_rt_hip_ry;
		bduVec3f pos_foot_center_rt_knee_ry;
		
		dog.getJointOffset( (LittleDog::LegIndex)leg, 
			LittleDog::HIP_RX_OFFSET, &pos_hip_rx_rt_body );
		
		dog.getJointOffset( (LittleDog::LegIndex)leg, 
			LittleDog::HIP_RY_OFFSET, &pos_hip_ry_rt_hip_rx );
		
		dog.getJointOffset( (LittleDog::LegIndex)leg, 
			LittleDog::KNEE_RY_OFFSET, &pos_knee_ry_rt_hip_ry );
		
		dog.getJointOffset( (LittleDog::LegIndex)leg, 
			LittleDog::FOOT_OFFSET, &pos_foot_center_rt_knee_ry );
		
		for ( int i = 0; i < 3; i++ )
			p->m_pos_foot_center_rt_body[ leg ].n[i] = 
				pos_hip_rx_rt_body.n[i] + pos_hip_ry_rt_hip_rx.n[i] +
				pos_knee_ry_rt_hip_ry.n[i] + pos_foot_center_rt_knee_ry.n[i];
		
		bdu_log_printf( BDU_LOG_LEVEL_INFO,
			"The center of foot[%d] is (%f,%f,%f)\n", leg, 
			p->m_pos_foot_center_rt_body[leg].n[0],
			p->m_pos_foot_center_rt_body[leg].n[1],
			p->m_pos_foot_center_rt_body[leg].n[2] );
	}
}

//--------------------------------------------------------------------
void 
KinematicWalker::initialize_params( void )
{
	bdu_log_printf( BDU_LOG_LEVEL_INFO, 
		"KinematicWalker>> Starting parameters initialization\n");

	// these parameters are with respect to the center of the foot expressed
	// in the body coordinate frame (x forward, z up) when the foot is 
	// assumed to be in the zero-configuration location.  For convenience
	// all Y offsets are for the left side of the robot.  When these spline
	// defaults are used for the right side of the robot, the y's are negated.

	float horizontal_front_radius           =  0.040;
	float horizontal_hind_radius            =  0.040;
	float vertical_height_in_swing          =  0.075;
	float vertical_height_in_stance         =  0.05;
	float swing_duration                    =  0.500;
	float stance_duration                   =  swing_duration * 4.0;

	// Note that ground speed = ((horizontal_front_radius + horizontal_hind_radius) / (stance_duration)

	p->m_spline_defaults.swing_p0[0]        = -horizontal_hind_radius;
	p->m_spline_defaults.swing_p0[1]        = 0;
	p->m_spline_defaults.swing_p0[2]        = vertical_height_in_stance;
	p->m_spline_defaults.swing_p0_t         = 0;
	p->m_spline_defaults.swing_p1[0]        = 0;
	p->m_spline_defaults.swing_p1[1]        = 0;
	p->m_spline_defaults.swing_p1[2]        = vertical_height_in_swing;
	p->m_spline_defaults.swing_p1_t         = 0.2 * swing_duration;
	p->m_spline_defaults.swing_p2[0]        = 1.3 * horizontal_front_radius;
	p->m_spline_defaults.swing_p2[1]        = 0;
	p->m_spline_defaults.swing_p2[2]        = vertical_height_in_swing;
	p->m_spline_defaults.swing_p2_t         = 0.4 * swing_duration;
	p->m_spline_defaults.swing_p3[0]        = 1.4 * horizontal_front_radius;
	p->m_spline_defaults.swing_p3[1]        = 0;
	p->m_spline_defaults.swing_p3[2]        = 0.5 * ( vertical_height_in_swing + vertical_height_in_stance);
	p->m_spline_defaults.swing_p3_t         = 0.7 * swing_duration;
	p->m_spline_defaults.swing_p4[0]        = horizontal_front_radius;
	p->m_spline_defaults.swing_p4[1]        = 0;
	p->m_spline_defaults.swing_p4[2]        = vertical_height_in_stance;
	p->m_spline_defaults.swing_p4_t         = swing_duration;

	p->m_spline_defaults.stance_p0_t        = swing_duration;
	p->m_spline_defaults.stance_p1_t        = swing_duration + stance_duration;

	p->m_spline_defaults.t_offset[0]        = (-3.0/4.0) * (swing_duration + stance_duration);
	p->m_spline_defaults.t_offset[1]        = (-1.0/4.0) * (swing_duration + stance_duration);
	p->m_spline_defaults.t_offset[2]        = (-2.0/4.0) * (swing_duration + stance_duration);
	p->m_spline_defaults.t_offset[3]        = 0;

	p->m_spline_defaults.xyz_offset[0]      = 0;
	p->m_spline_defaults.xyz_offset[1]      = 0;
	p->m_spline_defaults.xyz_offset[2]      = 0;

	p->m_spline_defaults.t_scale            = 1.0;

	bdu_log_printf( BDU_LOG_LEVEL_INFO, 
		"KinematicWalker: done with parameters initialization\n");
}


//--------------------------------------------------------------------
void
KinematicWalker::initialize_spline_for_leg( LittleDog::LegIndex leg_num )
{
	SwingStanceGenerator * & s = p->m_splines[ leg_num ];
	
	if ( s != NULL )
	{
		bdu_log_printf( BDU_LOG_LEVEL_WARN,  
			"initialize_spline_for_leg() called redundantly.\n" );
		return;
	}

	s =	new SwingStanceGenerator( 5, 2 );
	initialize_spline_from_params( s, leg_num );
	s->set_t_offset( p->m_spline_defaults.t_offset[ leg_num ] );
}

#define ADD_OFFSETS(p1,p2,p3) bduVec3f(p1[0]+p2.n[0]+p3.n[0], p1[1]+p2.n[1]+p3.n[1], p1[2]+p2.n[2]+p3.n[2])

//--------------------------------------------------------------------
void
KinematicWalker::initialize_spline_from_params( SwingStanceGenerator *s, 
	LittleDog::LegIndex leg_num )
{
	// reads from p->m_spline_defaults and copies applicable values as is.

	// foot center position
	bduVec3f & foot_center = p->m_pos_foot_center_rt_body[ leg_num ];
	bduVec3f foot_offset;

	// note that we don't spread here using s->set_xyz_offset() because a
	// sinusoidal sway mechanism uses s->set_xyz_offset() elsewhere ...
	switch (leg_num)
	{
	case LittleDog::FL:
		s->set_xyz_offset( bduVec3f(0.005,0,0) );  // a constant spline offset
		foot_offset = bduVec3f( 0,0.025,0 );        // spread
		break;
		
	case LittleDog::FR:
		s->set_xyz_offset( bduVec3f(0.005,0,0) );  // a constant spline offset
		foot_offset = bduVec3f( 0,-0.025,0 );       // spread
		break;
		
	case LittleDog::HR:
		s->set_xyz_offset( bduVec3f(-0.05,0,0) ); // a constant spline offset
		foot_offset = bduVec3f( 0,-0.025,0 );       // spread
		break;

	case LittleDog::HL:
		s->set_xyz_offset( bduVec3f(-0.05,0,0) ); // a constant spline offset
		foot_offset = bduVec3f( 0,0.025,0 );        // spread
		break;

	default:
		break;
	}

	// swing point 0
	float *point = &(p->m_spline_defaults.swing_p0[0]);

	s->set_swing_spline_point(0, ADD_OFFSETS(point,foot_center,foot_offset) );
	s->set_swing_spline_time(0, p->m_spline_defaults.swing_p0_t);

	// swing point 1
	point = &(p->m_spline_defaults.swing_p1[0]);

	s->set_swing_spline_point(1, ADD_OFFSETS(point,foot_center,foot_offset) );
	s->set_swing_spline_time(1, p->m_spline_defaults.swing_p1_t);

	// swing point 2
	point = &(p->m_spline_defaults.swing_p2[0]);

	s->set_swing_spline_point(2, ADD_OFFSETS(point,foot_center,foot_offset) );
	s->set_swing_spline_time(2, p->m_spline_defaults.swing_p2_t);

	// swing point 3
	point = &(p->m_spline_defaults.swing_p3[0]);

	s->set_swing_spline_point(3, ADD_OFFSETS(point,foot_center,foot_offset) );
	s->set_swing_spline_time(3, p->m_spline_defaults.swing_p3_t);

	// swing point 4
	point = &(p->m_spline_defaults.swing_p4[0]);

	s->set_swing_spline_point(4, ADD_OFFSETS(point,foot_center,foot_offset) );
	s->set_swing_spline_time(4, p->m_spline_defaults.swing_p4_t);

	// stance point 0 (same as swing point 4)
	point = &(p->m_spline_defaults.swing_p4[0]);

	s->set_stance_spline_point(0, ADD_OFFSETS(point,foot_center,foot_offset) );
	s->set_stance_spline_time(0, p->m_spline_defaults.stance_p0_t);

	// stance point 1 (same as swing point 0)
	point = &(p->m_spline_defaults.swing_p0[0]);

	s->set_stance_spline_point(1, ADD_OFFSETS(point,foot_center,foot_offset) );
	s->set_stance_spline_time(1, p->m_spline_defaults.stance_p1_t);
}


//--------------------------------------------------------------------
void
KinematicWalker::initialize_dataset( void )
{
	p->m_dataset = new LittleDogDataSetFileStreamWriter(true, "my_dataset.data");
	int var_num  = 0;

	p->m_vars[var_num++] = p->m_dataset->addVariable( "robot_time" );
	p->m_vars[var_num++] = p->m_dataset->addVariable( "spline_time" );
	p->m_vars[var_num++] = p->m_dataset->addVariable( "sway" );

	// x_d, y_d, and z_d are what xyz we want to ask for.  	the fkik 
	// version of these tells us what xyz we "effectively" asked for
	// see logging area.
	p->m_vars[var_num++] = p->m_dataset->addVariable( "FL:foot_x_d" );
	p->m_vars[var_num++] = p->m_dataset->addVariable( "FL:foot_y_d" );
	p->m_vars[var_num++] = p->m_dataset->addVariable( "FL:foot_z_d" );
	p->m_vars[var_num++] = p->m_dataset->addVariable( "FL:foot_fkik_x_d" );
	p->m_vars[var_num++] = p->m_dataset->addVariable( "FL:foot_fkik_y_d" );
	p->m_vars[var_num++] = p->m_dataset->addVariable( "FL:foot_fkik_z_d" );
	p->m_vars[var_num++] = p->m_dataset->addVariable( "FL:servo_foot_pos_res");
	p->m_vars[var_num++] = p->m_dataset->addVariable( "FL:sway_warp" );

	p->m_vars[var_num++] = p->m_dataset->addVariable( "FR:foot_x_d" );
	p->m_vars[var_num++] = p->m_dataset->addVariable( "FR:foot_y_d" );
	p->m_vars[var_num++] = p->m_dataset->addVariable( "FR:foot_z_d" );
	p->m_vars[var_num++] = p->m_dataset->addVariable( "FR:foot_fkik_x_d" );
	p->m_vars[var_num++] = p->m_dataset->addVariable( "FR:foot_fkik_y_d" );
	p->m_vars[var_num++] = p->m_dataset->addVariable( "FR:foot_fkik_z_d" );
	p->m_vars[var_num++] = p->m_dataset->addVariable( "FR:servo_foot_pos_res");
	p->m_vars[var_num++] = p->m_dataset->addVariable( "FR:sway_warp" );

	p->m_vars[var_num++] = p->m_dataset->addVariable( "HL:foot_x_d" );
	p->m_vars[var_num++] = p->m_dataset->addVariable( "HL:foot_y_d" );
	p->m_vars[var_num++] = p->m_dataset->addVariable( "HL:foot_z_d" );
	p->m_vars[var_num++] = p->m_dataset->addVariable( "HL:foot_fkik_x_d" );
	p->m_vars[var_num++] = p->m_dataset->addVariable( "HL:foot_fkik_y_d" );
	p->m_vars[var_num++] = p->m_dataset->addVariable( "HL:foot_fkik_z_d" );
	p->m_vars[var_num++] = p->m_dataset->addVariable( "HL:servo_foot_pos_res");
	p->m_vars[var_num++] = p->m_dataset->addVariable( "HL:sway_warp" );

	p->m_vars[var_num++] = p->m_dataset->addVariable( "HR:foot_x_d" );
	p->m_vars[var_num++] = p->m_dataset->addVariable( "HR:foot_y_d" );
	p->m_vars[var_num++] = p->m_dataset->addVariable( "HR:foot_z_d" );
	p->m_vars[var_num++] = p->m_dataset->addVariable( "HR:foot_fkik_x_d" );
	p->m_vars[var_num++] = p->m_dataset->addVariable( "HR:foot_fkik_y_d" );
	p->m_vars[var_num++] = p->m_dataset->addVariable( "HR:foot_fkik_z_d" );
	p->m_vars[var_num++] = p->m_dataset->addVariable( "HR:servo_foot_pos_res");
	p->m_vars[var_num++] = p->m_dataset->addVariable( "HR:sway_warp" );
}

//--------------------------------------------------------------------
// CONTROL CODE - below source implements actual control           
//--------------------------------------------------------------------
void
KinematicWalker::update( void )
{
	( this->*m_state )();
}

//--------------------------------------------------------------------
void
KinematicWalker::update_time( void )
{
	// keep m_spline_t in the range of [0,cycle_duration]

	float current_t;
	dog.getDataTimestamp(&current_t);
	float dt        = p->m_spline_defaults.t_scale * (current_t - p->m_prev_t);
	p->m_prev_t     = current_t;

	float cycle_duration = p->m_splines[0]->get_cycle_duration();

	p->m_spline_t += dt;
	if ( p->m_spline_t > cycle_duration ) // time is forward
		p->m_spline_t -= cycle_duration;
	else if ( p->m_spline_t  < 0 )        // time is backwards
		p->m_spline_t += cycle_duration;

	// sanity check
	if ( (p->m_spline_t > p->m_spline_defaults.stance_p1_t) || 
		(p->m_spline_t < p->m_spline_defaults.swing_p0_t) )
	{
		bdu_log_printf( BDU_LOG_LEVEL_INFO,  
			"KinematicWalker::update_time() at t = %2.3f is resetting t\n: %2.3f > %2.3f > %2.3f\n", 
			current_t,
			p->m_spline_defaults.swing_p0_t,
			p->m_spline_t,
			p->m_spline_defaults.stance_p1_t);
	  
		// Reset p->m_spline_t to p->m_spline_defaults.swing_p0_t
		p->m_spline_t = p->m_spline_defaults.swing_p0_t;
		p->m_prev_t = current_t;
	}
}

//--------------------------------------------------------------------
//
// Determine the knee bend on a particular leg
//
bool
getKneeBendForward( LittleDog::LegIndex leg )
{
	// returns true if the knee should bend forward, false otherwise

	return ( (leg == LittleDog::FL) || (leg == LittleDog::FR) ) ? false : true;
}

//--------------------------------------------------------------------
void
KinematicWalker::update_splines( void )
{
	float sway_amplitude = 0.035;
	float sway = sway_amplitude * sin (2.0f * 3.14159 * 
		(p->m_spline_t/p->m_spline_defaults.stance_p1_t) + 3.14159 * 215.0 /180.0);

	p->m_dataset->setData( p->m_vars[ 2 ], sway );

	for ( int leg_num = 0; leg_num < LittleDog::NUM_LEGS; ++leg_num )
	{
		bduVec3f offset = p->m_splines[leg_num]->get_xyz_offset();
		offset.n[1] = sway;
		p->m_splines[leg_num]->set_xyz_offset( offset );
	}
}

//--------------------------------------------------------------------
//
// Input "index" is the leg number.
//
void
KinematicWalker::move( float t, LittleDog::LegIndex index )
{
	SwingStanceGenerator * spline = p->m_splines[ index ];
	float x, y, z, dx, dy, dz;
	
	// spline->solve() logs internally
	if ( spline->solve( t, &x, &y, &z, &dx, &dy, &dz ) != 0 )
		return;

	const bduVec3f & offset = spline->get_xyz_offset();
	float x_offset_delta = x - offset.n[0];
	float x_scale = 0;
	float t_scale = 1.0;
	
	// turn amount can be ramped up when walking fast since it makes walking more stable
	if ( p->m_spline_defaults.t_scale > 1 )
		t_scale = p->m_spline_defaults.t_scale;
	
	if ( (index == LittleDog::FL) || (index == LittleDog::HL) )
		x_scale = t_scale * MAX_X_SCALE_WARP * p->m_turn_amount;
	else
		x_scale = - t_scale * MAX_X_SCALE_WARP * p->m_turn_amount;


	// without turning we would only do a moveTo( index, x, y, z, dx, dy, dz );
	float sway_warp = compute_sway_turn_warp( index, t );

	int var_num = MY_NON_LEG_VARS + (index * MY_VARS_PER_LEG) + MY_VARS_PER_LEG - 1;
	p->m_dataset->setData( p->m_vars[var_num], sway_warp );

	moveTo( index, x + (x_scale * x_offset_delta), y + sway_warp, z, dx + (x_scale * dx), dy, dz );
}

//--------------------------------------------------------------------
//
// Input "index" is the leg number.
// Input "t" is global spline time
//
float
KinematicWalker::compute_sway_turn_warp( LittleDog::LegIndex index, 
	float spline_t )
{
	// normalize time from [0 to cycle_duration] where swing is at t = 0
	const float cycle_duration  = p->m_spline_defaults.stance_p1_t;
	const float swing_duration  = p->m_spline_defaults.stance_p0_t;
	const float stance_duration = cycle_duration - swing_duration;
	float t = spline_t + p->m_spline_defaults.t_offset[index];
	float warp_amount;
	if ( t < 0 )
		t += cycle_duration;

	if ( t <= swing_duration )
	{
		// the leg is in swing from t = 0 to swing_duration,
		// so gradually ramp down the Y warp magnitude

		warp_amount = ((swing_duration - t) / swing_duration) * 
			p->m_turn_amount * MAX_Y_SWAY_WARP;
	}
	else
	{
		// the leg is in stance from t = swing_duration to cycle_duration,
		// so gradually ramp up the Y warp magnitude
		
		warp_amount = ((t - swing_duration) / stance_duration) * 
			p->m_turn_amount * MAX_Y_SWAY_WARP;
	}

	if ( (index == LittleDog::FL) || (index == LittleDog::FR) )
		warp_amount = -warp_amount;

	return -warp_amount;
}

//--------------------------------------------------------------------
//
// Input "index" is the leg number.
//
void
KinematicWalker::moveTo( LittleDog::LegIndex index, float x, float y, float z, 
	float dx, float dy, float dz )
{
	bool bend_forward = getKneeBendForward(index);

	// servo the foot
	LD_ERROR_CODE res = dog.servoFootPosition( index, bend_forward,
		bduVec3f(x, y, z), bduVec3f(dx, dy, dz) );
	
	// for logging
	if ( p->m_log_count < MY_MAX_LOGGED )
	{
		int var_num = MY_NON_LEG_VARS + (index * MY_VARS_PER_LEG);

		p->m_dataset->setData( p->m_vars[ var_num++ ], x ); // foot_xd
		p->m_dataset->setData( p->m_vars[ var_num++ ], y ); // foot_yd
		p->m_dataset->setData( p->m_vars[ var_num++ ], z ); // foot_zd

		// internally, servoFootPosition() above may return errors, logging the results of IK and FK so as a test
		bduVec3f ikfk, forward, rearward;
		bduVec4f orient;
		dog.computeFootIK( index, bduVec3f(x,y,z), &forward, &rearward );
		dog.computeFootFK( index, (bend_forward == true) ? forward : rearward, &ikfk, &orient );

		p->m_dataset->setData( p->m_vars[ var_num++ ], ikfk.n[0] ); // foot_fkik_xd
		p->m_dataset->setData( p->m_vars[ var_num++ ], ikfk.n[1] ); // foot_fkik_yd
		p->m_dataset->setData( p->m_vars[ var_num++ ], ikfk.n[2] ); // foot_fkik_zd

		p->m_dataset->setData( p->m_vars[ var_num++ ], res ); // servo_foot_pos_res

		p->m_dataset->saveSample();
	}	

	// check results
	if ( res != LD_OKAY )
		bdu_log_printf( BDU_LOG_LEVEL_WARN, 
			"servoFootPosition for leg %d returned %d\n", index, res );
}

//--------------------------------------------------------------------
void 
KinematicWalker::setTimeScale( float scale )
{
	if ( scale > 1.5 )
		scale = 1.5;
	else if ( scale < 0.25 )
		scale = 0.25;

	p->m_spline_defaults.t_scale = scale;
}

//--------------------------------------------------------------------
float
KinematicWalker::getTimeScale( void )
{
	return p->m_spline_defaults.t_scale;
}

//--------------------------------------------------------------------
void
KinematicWalker::setTurnAmount(float turn)
{
  // the inputs should fall within range -1 to +1.

	if ( turn > 1 )
		turn = 1;
	if ( turn < -1)
		turn = -1;
	
	p->m_turn_amount = turn;
}

//--------------------------------------------------------------------
float
KinematicWalker::getTurnAmount( void )
{
	return p->m_turn_amount;
}

//--------------------------------------------------------------------
float 
KinematicWalker::getWalkVelocity( void )
{
	float dist = p->m_spline_defaults.swing_p4[0] - p->m_spline_defaults.swing_p0[0];
	float duration = (p->m_spline_defaults.stance_p1_t - p->m_spline_defaults.stance_p0_t) / 
		p->m_spline_defaults.t_scale;

	return dist / duration;  // meters per second
}

//--------------------------------------------------------------------
bool 
KinematicWalker::inStance( void )
{
	bool in_stance = true;
	
	for ( int leg = 0; leg < LittleDog::NUM_LEGS; leg++ )
		in_stance &= p->m_splines[ leg ]->in_stance();
	
	return in_stance;
}

//--------------------------------------------------------------------
void 
KinematicWalker::getStartingFootPos(bduVec3f pos[ LittleDog::NUM_LEGS ])
{
	// learn the desired location of each foot at the beginning of walk
	// and store them in the m_targeted_foot_pos.  This function should not
	// be called until after the walker has been initialized().  This function
	// should not be called while the trial is actively doing something.

	float current_t;
	dog.getDataTimestamp(&current_t);
	
	p->m_spline_t = p->m_spline_defaults.swing_p0_t;
	p->m_prev_t = current_t;
	
	update_splines();

	for ( int leg = 0; leg < LittleDog::NUM_LEGS; leg++ )
	{
		SwingStanceGenerator * spline = p->m_splines[ leg ];
		float x, y, z, dx, dy, dz;
		if ( spline->solve( p->m_spline_t, &x, &y, &z, &dx, &dy, &dz ) != 0 )
			pos[ leg ] = p->m_pos_foot_center_rt_body[ leg ];
		else
			pos[ leg ] = bduVec3f(x,y,z);
	}
}

//--------------------------------------------------------------------
void
KinematicWalker::stateInit( void )
{
	m_state = &KinematicWalker::stateWalk;
	bdu_log_printf( BDU_LOG_LEVEL_INFO, 
		"KinematicWalker::Transitioning to WALK state\n" );
}

//--------------------------------------------------------------------
//
// In this state, time progresses according to the t_scale * m_dt.
// The splines are solved the for the t, and the leg's desired
// positions are updated each dt.
//
void
KinematicWalker::stateWalk( void )
{
	update_splines();
	update_time();
  
	// for logging
	float robot_t;
	dog.getDataTimestamp(&robot_t);
	
	p->m_dataset->setData( p->m_vars[ 0 ], robot_t );
	p->m_dataset->setData( p->m_vars[ 1 ], p->m_spline_t );

	// move the legs
	int leg;
	for ( leg = 0; leg < LittleDog::NUM_LEGS; ++leg )
		move( p->m_spline_t, (LittleDog::LegIndex)leg );

	if ( p->m_log_count < MY_MAX_LOGGED )
		p->m_log_count++;
}

//====================================================================
