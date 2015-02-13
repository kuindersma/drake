#ifndef __KINEMATIC_WALKER_H__
#define __KINEMATIC_WALKER_H__

#include <littledog.h>

class KinematicWalker
{
 public:
	KinematicWalker( class MyLittleDog & );
	virtual ~KinematicWalker();
  
  virtual void initialize( void );
  virtual void update( void );
  virtual void uninitialize( void );
	
	// returns the starting [x,y,z] location of each foot, in body coordinates.
	void  getStartingFootPos(bduVec3f m_pos_foot_center_rt_body[ LittleDog::NUM_LEGS ]);

	// sets the time scale of the walk.  The ranges is [0.25 to 1.5] with 1.0 being nominal.
	// A time scale of 0.25 walks one fourth as fast as a time scale of 1.0.
	void  setTimeScale( float scale );

	// returns the current time scale.  See setTimeScale()
	float getTimeScale( void );

	// returns the theoretical velocity of the walk for the current time scale.
	float getWalkVelocity( void );

	// sets the turn amount of the walk.  The range is [-1.0 to +1.0] with 0.0 being
	// nominal (straight), -1 as full strength left turn, and +1 as full strength right.
	void setTurnAmount( float amount );

	// returns the current turn amount.  See setTurnAmount()
	float getTurnAmount( void );

	// returns true if each of the robot's legs is currently in stance.
	bool inStance( void );

 private:
	void stateInit(void);
	void stateWalk(void);

	void move( float t,  LittleDog::LegIndex index );
	void update_time( void );
	void moveTo(  LittleDog::LegIndex index, float x, float y, float z, 
		float dx, float dy, float dz );
	
	void initialize_spline_from_params( class SwingStanceGenerator *s, 
		LittleDog::LegIndex leg_num );
	void initialize_spline_for_leg( LittleDog::LegIndex index );
	void initialize_params();
	void initialize_dataset();
	void initialize_foot_centers();
	void update_splines();

	float compute_sway_turn_warp( LittleDog::LegIndex index, float spline_t );

	void (KinematicWalker:: *  m_state )(void);
	struct KinematicWalkerPrivate * p;

	class MyLittleDog & dog;
};


#endif  /*  __KINEMATIC_WALKER_H__  */
