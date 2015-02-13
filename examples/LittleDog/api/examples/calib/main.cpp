
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <errno.h>
#include <math.h>
#include <pthread.h>
#include <signal.h>

#include <littledog.h>

//====================================================================

class MyLittleDog : public LittleDog 
{
public:
	MyLittleDog() 
		: shutdown( false ) { }

	bool shutdown;
protected:


	virtual bool initControl( void );
	virtual void updateControl( void );
};

static MyLittleDog * dog = NULL;

#define CHK_ERROR(exp)              \
   { int r = LD_OKAY;             \
      if( (r = exp) != LD_OKAY ) {  \
         printf("CalCheck: error %d on line %d\n", r, __LINE__ ); \
         abortConnection(); }; }

static void signal_handler( int )
{
	printf( "Stopping control process\n" );
	assert( dog );
	dog->shutdown = true;
	dog->abortConnection();
}

int 
main( int, char * ) 
{
	signal( SIGINT, signal_handler );

	MyLittleDog l;
	dog = &l;

	printf( "Starting control \n" );
	l.initializeRobot();
	l.calibrate( true, false ); // args: force calibrate, end up standing

	if ( l.shutdown )
		return 0;

	printf("\n\n-----------------------------------------------------------------------\n");
	printf( "We're now going to validate your robot's calibration. Please place the\n");
	printf( "LittleDog on level ground, and press any key to continue.\n");
	printf("-----------------------------------------------------------------------\n\n");
	getchar();

	l.runTrial();

	// For now just hit CTRL+C to terminate
	printf( "Exiting main.\n" );
	fflush( stdout );
	return 0;
}


#define RAD2DEG(exp) ((180.0f/3.14159)*(exp))

bool
MyLittleDog::initControl( void )
{
	donePlanning();
	return true;
}


/// The calibration example assumes that you've just performed a calibration
/// and have set the robot on a level table before you entered updateControl.
/// Update control proceeds through two sets of checks:
///  1. it waits for the robot to be settled dynamically and for all the
///     joints to be in the 0 position.
///  2. it checks that the robot is level based on the gyro angles and that
///     the force is reasonably divided between the feet.
void
MyLittleDog::updateControl( void )
{
	// maximum allowable angular rate of the body
	static const float MAX_ROTATIONAL_MAG  = 0.25;   // rps

	// maximum acceleration deviation from gravity
	static const float MAX_ACCEL_MAG       = 3.0;    // m/sec^2

	// maximum combined roll,pitch deviation from level
	static const float MAX_ANG_NOT_LEVEL   = 0.05;  // rad = ~3.0deg

	// maximum difference in leg forces (since these are poorly calibrated
	// this is very loose).
	static const float MAX_FORCE_IMBALANCE = 0.1;  // plus/minus of required weight

	// how close to the zero position the robot's joints must servo.
	static const float NEAR_ZERO_TOLERANCE = 0.001;  // rad

	// gravity on earth
	static const float GRAVITY = 9.812; //m/sec^2

	// check for the shutdown flag, if true stop the robot (ctrl-c handler)
	if( shutdown )
		abortConnection();

	// 1.  setup joints to servo to the 0 position, wait for them to get
	//     suitably close.
	float k[] = { 7.5, 7.5, 6.7 };
	float b[] = { 0.32, 0.32, 0.32 };
	bool servoing_done = true;

	for(int l = 0; l < NUM_LEGS; l++)
	{
		LegIndex leg_index = static_cast< LegIndex >( l );
		CHK_ERROR( setLegCommandMode( leg_index, LEG_MODE_INDEPENDENT ) );
		
		for( int j = 0; j < NUM_JOINTS; j++ )
		{
			JointIndex joint_index = static_cast< JointIndex >( j );
			CHK_ERROR( setJointServoType( leg_index, joint_index, SERVO_PD ) );
			CHK_ERROR( setJointPDServoGains( leg_index, joint_index, k[ j ], b[ j ] ) );
			CHK_ERROR( setJointPDServoSetPoints( leg_index, joint_index, 0.0 ) );

			float q = getJointState( leg_index, joint_index, &q );

			servoing_done &= ( fabs( q ) < NEAR_ZERO_TOLERANCE );
		}
	}

	if( !servoing_done )
		return;

	// 2. grab the foot forces
	bduVec3f force;
	float forces[NUM_LEGS];
	float sum_forces = 0.0f;
	for(int i = 0 ; i < NUM_LEGS; i++)
	{
		CHK_ERROR( getFootForce( (LegIndex)i, &force ) );
		forces[i] = force.n[2];
		sum_forces += forces[i];
	}
	printf("Forces: %f %f %f %f = %f\n", forces[0], forces[1], forces[2], forces[3], sum_forces );

	// 3. grab the imu data
	bduVec4f orientation;
	bduVec3f angrates;
	bduVec3f accel;
	CHK_ERROR( getIMUReadings( &orientation, &angrates, &accel ) );


	// 4a. Check that we are still... is the body of the robot undergoing rotation?
	if( sqrt( angrates.n[ 0 ] * angrates.n[ 0 ] 
			+ angrates.n[ 1 ] * angrates.n[ 1 ]
			+ angrates.n[ 2 ] * angrates.n[ 2 ]) > MAX_ROTATIONAL_MAG )
	{
		printf("CalCheck: Undergoing rotation. %f %f %f. Wait... \n", angrates.n[0], angrates.n[1], angrates.n[2] );
		return;
	}


	// 4b. Check that we are still... are we accelerating linearly?
	if( fabs( sqrt( accel.n[ 0 ] * accel.n[ 0 ] 
			+ accel.n[ 1 ] * accel.n[ 1 ]
			+ accel.n[ 2 ] * accel.n[ 2 ]) - GRAVITY ) > MAX_ACCEL_MAG )
	{
		printf("CalCheck: Undergoing acceleration. %f %f %f. Wait... \n", accel.n[0], accel.n[1], accel.n[2] );
		return;
	}
	
	bool cal_failure   = false;


	// 5. If execution reaches here, the robot is dynamically still. Is it level?
	{
		const float angle = 2 * acos(orientation.n[3]);
		bduVec3f axis(orientation.n[0],	orientation.n[1],	orientation.n[2]);
		const float mag = sqrt(axis.n[0]*axis.n[0] + axis.n[1]*axis.n[1] + axis.n[2]*axis.n[2]);
		if( mag > 0.01 )
		{
			for(int i=0; i<3; i++) axis.n[i] *= angle/mag;
			const float xy_mag = sqrt(axis.n[0]*axis.n[0] + axis.n[1]*axis.n[1]);
			
			if( xy_mag > MAX_ANG_NOT_LEVEL )
			{
				// out of level is a calibration failure.
				printf( "CalCheck FAILURE: the robot is out of level in roll/pitch by %2.1f deg.\n",
					RAD2DEG(xy_mag));
				cal_failure = true;
			}
		}
		else
		{  // mag will only be zero near level, so count this as success.
			
		}
	}

	// 6. ...and are the foot forces balanced sufficiently?
	for(int i = 0; i < NUM_LEGS; i++)
	{
		if( ((forces[i] / sum_forces ) - 0.25) > MAX_FORCE_IMBALANCE )
		{
			// too much weight on a foot is a calibration failure.
			printf( "CalCheck FAILURE: the robot is carrying too much weight on foot %d\n", i );
			cal_failure = true;
		}
		else if ((( forces[i]/sum_forces ) + 0.25) < -MAX_FORCE_IMBALANCE )
		{
			// too little weight on a foot is a calibration failure.
			printf( "CalCheck FAILURE: the robot is carrying too little weight on foot %d\n", i );
			cal_failure = true;
		}
	}

	// 7. if execution reaches here, then we have performed all tests and can determine if
	// we've passed or failed.
	if( cal_failure )
	{		
		printf("\n\n-----------------------------------------------------------------------\n");	
		printf( "Calibration failed!\n" );
		printf("-----------------------------------------------------------------------\n");	
		stopRobot();
	}
	else
	{
		printf("\n\n-----------------------------------------------------------------------\n");	
		printf( "Calibration was SUCCESSFUL!\n" );
		printf("-----------------------------------------------------------------------\n");	
		stopRobot();
	}
	
}
//====================================================================
