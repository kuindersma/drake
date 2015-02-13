
//
//  System includes.
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <errno.h>
#include <math.h>
#include <signal.h>
#include <unistd.h>

//
//  LittleDog API includes.
//
#include <littledog.h>
#include <bduLog.h>

//
//  Forward declarations.
//
static void handle_result( const char* func, LD_ERROR_CODE result, 
	bool exit_on_error);
static void signal_handler( int );


//===================================================================

class MyLittleDog : public LittleDog
{
public:
  
	bool first_time;

	virtual bool initControl( void ) 
	{
		printf( "initControl(): called\n");

		first_time = true;

		return true; // tells the caller that we are ok
	}

	virtual void updateControl( void )
	{
		LD_ERROR_CODE result;

		// one-time only
		if ( first_time )
		{
			first_time = false;
			
			// Initialize the PD servos the first_time time through the loop
			
			printf( "updateControl(): Initializing all joints to use PD servos\n" );
			
			float k[] = { 7.5, 7.5, 6.7 };
			float b[] = { 0.32, 0.32, 0.32 };
			
			for ( int l=0; l<NUM_LEGS; ++l )
				for ( int j=0; j<NUM_JOINTS; ++j )
				{
					result = setJointServoType( (LegIndex)l, (JointIndex)j, SERVO_PD );
					handle_result("setJointServoType", result, true);

					result = setJointPDServoGains( (LegIndex)l, (JointIndex)j, k[j], b[j] );
					handle_result("setJointPDServoGains", result, true);
				}			

			donePlanning();
		}

		// Do control for this particular timestamp

		// For now just hit CTRL+C to terminate

		float t;
		getDataTimestamp( &t );
		
		// hip_rx_limits are changed from (-0.6, 0.6) so that the test
		// can be run while the dog is on its test stand (legs are spread)
		float hip_rx_limits[2]  = {0, 0.6};    
		float hip_ry_limits[2]  = {-3.5, 2.4};
		float knee_ry_limits[2] = {-3.1, 1.0};

		// averages
		float hip_rx_avg  = (hip_rx_limits[1]  + hip_rx_limits[0])/2.0f;
		float hip_ry_avg  = (hip_ry_limits[1]  + hip_ry_limits[0])/2.0f;
		float knee_ry_avg = (knee_ry_limits[1] + knee_ry_limits[0])/2.0f;

		// amplitudes
		float hip_rx_amp  = (hip_rx_limits[1]  - hip_rx_limits[0])/2.0f;
		float hip_ry_amp  = (hip_ry_limits[1]  - hip_ry_limits[0])/2.0f;
		float knee_ry_amp = (knee_ry_limits[1] - knee_ry_limits[0])/2.0f;

		// desired values
		float q_d1  = hip_rx_amp  * sin( 3.14 * 0.5 * t );
		float q_d2  = hip_ry_amp  * cos( 4.14 * 0.5 * t );
		float q_d3  = knee_ry_amp * sin( 5.14 * 0.5 * t );

		for ( int l=0; l<NUM_LEGS; ++l )
		{
			// exceptions
			float avg_k1 = 1.0, avg_k2 = 1.0, avg_k3 = 1.0;

			if ( (l == 1) || (l == 3) )
			{
				q_d1 *= -1;
				avg_k1 = -1;
			}

			if ( (l == 2) || (l == 3) )
			{
				avg_k2 = -1;
				avg_k3 = -1;
			}

			// set desired values
			result = setJointPDServoSetPoints( (LegIndex)l, HIP_RX, 
				(0.75f * q_d1) + (avg_k1 * hip_rx_avg) ); 
			handle_result("setJointPDServoPoints", result, false);

			result = setJointPDServoSetPoints( (LegIndex)l, HIP_RY, 
				(0.50f * q_d2) + (avg_k2 * hip_ry_avg) );
			handle_result("setJointPDServoPoints", result, false);

			result = setJointPDServoSetPoints( (LegIndex)l, KNEE_RY, 
				(0.75f * q_d3) + (avg_k3 * knee_ry_avg) );
			handle_result("setJointPDServoPoints", result, false);
		}
	}

	virtual void uninitControl( void ) 
	{
		printf( "uninitControl(): Cleaning up ...\n" );
	}

};

//===================================================================

static MyLittleDog * dog = NULL;

static void signal_handler( int )
{
	printf( "Stopping control process\n" );
	fflush( stdout );
  
	// We cannot delete 'dog' right here because it may be in use
	// in, say, the run or calibrate process.   By aborting the 
	// connection, it assures that these processes end quickly.

	LD_ERROR_CODE result = dog->abortConnection();

	handle_result("signal_handler::abortConnection", result, false); // must be false
}

//===================================================================

static void handle_result( const char* func, LD_ERROR_CODE result, bool exit_on_error)
{
	// if an error occurs, print it out and (optionally) exit

	if ( result == LD_OKAY )
		return;
	
	printf( "An error occurred calling '%s()'.\n\tError was: '%s'\n", 
		func, dog->getErrorCodeString(result));
	
	if ( exit_on_error )
	{
		delete dog;
		dog = NULL;
		exit(EXIT_FAILURE);
	}
}


//===================================================================
int 
main( int argc, char * argv[] ) 
{
	signal( SIGINT, signal_handler );

	dog = new MyLittleDog();

	bool forced = false;

	// pass in "-f" if you want to force calibration always
	if ( (argc == 2) && (!strcmp(argv[1],"-f")) )
		forced = true;

	LD_ERROR_CODE result;

	printf( "Connecting to Mocap ...\n" );
	result = dog->initializeMocap();
	handle_result("initializeMocap", result, false);

	printf( "Connecting to robot ...\n" );
	result = dog->initializeRobot();
	handle_result("initializeRobot", result, true);

	printf( "%s calibration of robot ...\n", forced ? "Forced " : "Non-forced" );
	result = dog->calibrate(forced);
	handle_result("calibrate", result, true);

	printf( "Starting user control ...\n" );
	result = dog->runTrial();
	handle_result("runTrial", result, true);

	printf( "Exiting main.\n" );
	fflush( stdout );

	delete dog;
	dog = NULL;
	return 0;
}

//====================================================================
