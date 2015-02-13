
#ifndef __LITTLEDOG_API_H__
#define __LITTLEDOG_API_H__

#include <bduVec3f.h>
#include <bduVec4f.h>
#include <bduMat3f.h>

//!
//!  \addtogroup ld_api  LittleDog Primary API
//!

//
//  LittleDog software version numbers.
//
#define LITTLEDOG_SOFTWARE_VERSION_STRING "1.1.2"
#define LITTLEDOG_SOFTWARE_VERSION_MAJOR 1
#define LITTLEDOG_SOFTWARE_VERSION_MINOR 1
#define LITTLEDOG_SOFTWARE_VERSION_POINT 2

//!
//!  \enum    LD_ERROR_CODE
//!
//!  \brief   Error codes used by most API functions.
//!
//!  \ingroup ld_api  LittleDog Primary API
//!  @{
//!
enum LD_ERROR_CODE
{
	LD_OKAY = 0,                           //!<  No error
	LD_TIMED_OUT,                          //!<  Operation timed out
	LD_INVALID_FUNCTION_CALL,              //!<  Function can only be called from updateControl()
	LD_INVALID_REQUEST,                    //!<  Requested function is out of context or order
	LD_ROBOT_COMMS_UNINITIALIZED,          //!<  No active host/robot communications
	LD_ROBOT_COMMS_INITIALIZED,            //!<  Host/robot comms currently active
	LD_ROBOT_COMMS_FAILED,                 //!<  Unable to setup host/robot communications
	LD_MOCAP_COMMS_FAILED,                 //!<  Mocap system communication not established
	LD_BAD_LEG_INDEX,                      //!<  Invalid leg index
	LD_BAD_JOINT_INDEX,                    //!<  Invalid joint index
	LD_BAD_JOINT_OFFSET,                   //!<  Invalid joint offset
	LD_BAD_BODY_INDEX,                     //!<  Invalid body index
	LD_BAD_MARKER_INDEX,                   //!<  Invalid mocap marker index
	LD_BAD_TERRAIN_INDEX,                  //!<  Invalid terrain index
	LD_BAD_LEG_SERVO_TYPE,                 //!<  Invalid leg servo type
	LD_BAD_JOINT_SERVO_TYPE,               //!<  Invalid joint servo type
	LD_BAD_PROPORTIONAL_GAIN,              //!<  Invalid gain (< 0.0) (in function setJointPDServoGains())
	LD_BAD_DERIVATIVE_GAIN,                //!<  Invalid gain (< 0.0) (in function setJointPDServoGains())
	LD_ANGLE_OUT_OF_RANGE,                 //!<  Passed angle is out of accepted range. Value clipped.
	LD_SINGULAR_MATRIX,                    //!<  A mathematical operation produced a singularity
	LD_INSUFFICIENT_ARGUMENTS,             //!<  Not enough arguments
	LD_CANNOT_REACH_IK_KNEE_BENT_FORWARD,  //!<  No knee forward IK solution for the requested position
	LD_CANNOT_REACH_IK_KNEE_BENT_REARWARD, //!<  No knee backward IK solution for the requested position
	LD_CANNOT_REACH_IK_KNEE_BENT_EITHER,   //!<  No IK solution for the requested position
	LD_USER_INIT_CONTROL_FAILED,           //!<  subclass initControl() returned false
	LD_TERMINATED,                         //!<  Call was terminated
	LD_BAD_QUATERNION_VECTOR,              //!<  Invalid unit quaternion (not unit length)
	LD_ROBOT_STOP_IS_BROKEN,               //!<  One or more robot leg stops is broken!  Cease using the Robot

	LD_NUM_ERROR_CODES              //!<  Number of error codes in LD_ERROR_CODE
};
//!  @}


/****************************************************************************/
//!  \class   LittleDog   littledog.h
//!
//!    The LittleDog class is the primary programming interface to the
//!    LittleDog robot and the Vicon motion capture system.
//!
//!    This class is an abstract base class, so a subclass \e must
//!    be derived from it (typically MyLittleDog in provided examples)
//!    to instantiate an object of this type.
//!
//!    Once an object of a LittleDog subclass is created, there is
//!    a specific progression of states the system must proceed
//!    through to run trials.
//!
//!    A typical session would proceed as follows:
//!
//!      - user creates MyLittleDog object
//!      - user calls initializeRobot() to set up host/robot communications
//!      - user calls calibrate()
//!          - this will cause the robot to find the zero position for all joints
//!          - if the robot has already been calibrated it will return immediately, unless it is a forced calibration)
//!          - an option argument can be provided to force recalibration of a previously calibrated robot
//!      - user (optionally) calls initializeMocap() or stopMocap() to change the host/vicon communication status.
//!          - The mocap status will remain unchanged from trial to trial.
//!          - If enabled, the latests available data from the mocap will be recorded each control dt
//!      - user calls runTrial()
//!          - the LittleDog state changes to indicate the planning phase of a trial has begun.
//!          - LittleDog calls virtual initControl() once.  Users will now have access to all API calls.
//!          - LittleDog makes multiple calls to virtual updateControl()
//!              - user calls donePlanning() from within updateControl() or from some other thread, 
//!                ending the planning phase and starting the running phase
//!              - user calls stopRobot() from within updateControl() or from some other thread
//!          - LittleDog calls virtual uninitControl() once
//!      - user calls stopRobot() from within the trial.
//!      - user repeats paired calls to runTrial() and stopRobot() as often as desired
//!      - user destroys MyLittleDog object
//!
//!    initializeMocap() is optional and only needs to be called if 
//!    motion capture system will be used during the trial(s).  If a trial has
//!    been started without initializing the motion capture system, 
//!    the trial must be terminated prior to calling initializeMocap().
//!    Users can call stopMocap() or initializeMocap() prior to calling runTrial()
//!    to change the mocap connection status.
//!
//!    \e Note: If the host computer detects a loss of communication with the robot
//!             it will automatically cease communications with the robot and transition
//!             to STATE_ABORTING_CONNECTION.  Once in this state the best recourse is 
//!             to try to initializeRobot() again, although you could also terminiate your
//!             application and start over again.  You should consider adjusting 
//!             your wireless communications settings to improve network performance if you
//!             observe this happening frequently.
//!
//!  \ingroup ld_api  LittleDog Primary API
//!  @{
//!

class LittleDog
{

 /////////////////////////////////////////////////////////////////////////
 //
 //  Enumerations and Data Types
 //
 /////////////////////////////////////////////////////////////////////////

 public:

	//!
	//!  \enum
	//!
	//!  \brief   LittleDog API execution states.
	//!
	enum State
	{
		STATE_UNINITIALIZED = 0,  //!<  Initial state after object construction
		STATE_ROBOT_INITIALIZED,  //!<  Robot communications are initialized
		STATE_ROBOT_CALIBRATING,  //!<  Robot is in the process of calibrating
		STATE_ROBOT_CALIBRATED,   //!<  Robot has been calibrated
		STATE_TRIAL_PLANNING,     //!<  Trial is in the planning phase
		STATE_TRIAL_RUNNING,      //!<  Trial is in the running phase
		STATE_STOPPING_ROBOT,     //!<  Stopping calibration or trial run
		STATE_ABORTING_CONNECTION //!<  Aborting trials and connections
	};
  
	//!
	//!  \enum
	//!
	//!  \brief   Leg indices.
	//!
	enum LegIndex
	{
		FL = 0,      //!<  Front left
		FR,          //!<  Front right
		HL,          //!<  Hind left
		HR,          //!<  Hind right
		NUM_LEGS     //!<  Number of legs
	};
  
	//!
	//!  \enum
	//!
	//!  \brief   Identifies offset types for parts of the robot
	//!
	//!  \_Description
	//!
	//!    There is no 'BODY_OFFSET' only because it would result in (0,0,0).
	//!    See the kinematic documentation in the LittleDog User Manual to 
	//!    relate these offsets to the robot.  All offsets
	//!    are expressed relative to their parent body.  The parent of HIP_RX_OFFSET is
	//!    the body (or trunk).  The parent of HIP_RY_OFFSET is HIP_RX_OFFSET.
	//!    The parent of KNEE_RY_OFFSET is HIP_RY_OFFSET.  And finally, the
	//!    parent of FOOT_OFFSET is KNEE_RY_OFFSET.  These offsets fully specify
	//!    the robot kinematics when taken together with the frame conventions 
	//!    outlined in the User Manual.
	//!
	enum JointOffset
	{
		HIP_RX_OFFSET = 0, //!< Offset to hip_rx from body reference frame
		HIP_RY_OFFSET,     //!< Offset to hip_ry from hip_rx joint
		KNEE_RY_OFFSET,    //!< Offset to knee from hip_ry joint
		FOOT_OFFSET,       //!< Offset to foot from knee joint
		NUM_OFFSETS
	};

	//!
	//!  \enum
	//!
	//!  \brief   Joint indices.
	//!
	enum JointIndex
	{
		HIP_RX = 0,    //!<  Hip joint, roll about X axis
		HIP_RY,        //!<  Hip joint, pitch about Y axis
		KNEE_RY,       //!<  Knee joint, pitch about Y axis
		NUM_JOINTS     //!<  Number of joints
	};
  
	//!
	//!  \enum
	//!
	//!  \brief   Body indices.
	//!
	//!  \_Description
	//!
	//!   These can be used to index into the array of MoCap BodyInfo structures.
	//!
	enum BodyIndex
	{
		B_GROUND =  -1,    //!<  Ground (not actually part of robot)
		B_TRUNK  =   0,    //!<  Trunk
		B_FL_ULEG,         //!<  Front left upper leg
		B_FL_LLEG,         //!<  Front left lower leg
		B_FR_ULEG,         //!<  Front right upper leg
		B_FR_LLEG,         //!<  Front right lower leg
		B_HL_ULEG,         //!<  Hind left upper leg
		B_HL_LLEG,         //!<  Hind left lower leg
		B_HR_ULEG,         //!<  Hind right upper leg
		B_HR_LLEG,         //!<  Hind right lower leg
		NUM_BODIES         //!<  Number of bodies
	};

	//!
	//!  \enum
	//!
	//!  \brief   Joint servo types.
	//!
	enum JointServoType
	{
		SERVO_LIMP = 0,    //!<  Limp servo; joint motors off
		SERVO_PD,          //!<  PD joint level servo + optional foot force control
		SERVO_TORQUE,      //!<  Open-Loop torque commands
		NUM_SERVOS         //!<  Number of servo types
	};

	//!
	//!  \struct  JointInfo
	//!
	//!  \brief   Information about a single joint
	//!
	struct JointInfo
	{
		// Actual 
		float q;            //!<  Current joint angle in radians
		float qd;           //!<  Current joint velocity in radians/sec
		
		// Desired
		float q_d;          //!<  Desired (PD) joint angle in radians
		float qd_d;         //!<  Desired (PD) joint velocity in radians/sec
		float tau_d;        //!<  Currently commanded torque in Newton meters
		float k;            //!<  PD proportional gain
		float b;            //!<  PD differential gain
		float ff_d;         //!<  Desired (PD) feed-forward torque in Newton meters
		float user_min_tau; //!<  Desired min joint torque limit, in Newton meters
		float user_max_tau; //!<  Desired max joint torque limit, in Newton meters

		JointInfo( void );  //!<  Default constructor
	};

	//! 
	//!  \struct  ForceControlInfo
	//!
	//!  \brief   Information about a single leg (gains in hip_rx, hip_ry, knee_ry order)
	//!
	struct ForceControlInfo
	{
		// Desired
		bduMat3f gains;            //!< Control gain on force errors, in Newton meters/Newton
		bduVec3f f_d;              //!< Desired force at the foot [x,y,z], in Newtons

		ForceControlInfo( void );  //!<  Default constructor
	};

	//!
	//!  \struct  DefaultUserTauLimits
	//!
	//!  \brief   Tau limits used when Calibrating, or Limp, and by default
	//!
	enum DefaultUserTauLimits
 	{
		USER_DEFAULT_MIN_TAU	= -1000,
		USER_DEFAULT_MAX_TAU	=  1000
	};

	//!
	//!  \struct  LegInfo
	//!
	//!  \brief   Information about a single leg
	//!
	//!  \_Description
	//!
	//!    Each leg has information about the joints that comprise the leg as well as 
	//!    information about the ground contact force sensor included in the lower leg.
	//!
	struct LegInfo 
	{ 
		JointInfo         joints[ NUM_JOINTS ];  //!<  Information for each leg joint
		ForceControlInfo  force_control;         //!<  Desired Information for the whole leg
		float             force[3];              //!<  Measured force at the foot [x,y,z]

		LegInfo( void );                         //!<  Default constructor
	};

	//!
	//!  \struct  IMUInfo
	//!
	//!  \brief   Information from the inertia measurement unit (IMU)
	//!
	//!  \_Description
	//!
	//!    The inertia measurement unit (IMU) provides estimated information
	//!    about the orientation of the robot body as well as a measurment of the 
	//!    angular rates of the body and the acceleration experienced by the body.
	//!    Note that the acceleration vector includes the gravitational acceleration.
	//!
	struct IMUInfo
	{
		bduVec3f orientation;   //!<  Body orientation as estimated by the IMU in yaw-roll-pitch Euler angles
		bduVec3f rates;         //!<  Angular rate measurment from the IMU
		bduVec3f accel;         //!<  Body acceleration as measured by the IMU (including gravity) 

		IMUInfo( void );        //!<  Default constructor
	};

	//!
	//!  \struct  Status
	//!
	//!  \brief   System status information (robot and host)
	//!
	//!  \_Description
	//!
	//!    Status provides flags for each joint in each leg, and
	//!    includes 'torque saturation' (both kinds) and 'at limit' flags as
	//!    well as count of control cycles skipped and packets
	//!    queued.  These latter statistics can be used to evaluate
	//!    the performance of the host-robot communication system and 
	//!    to detect when user control software is running excessively slow.
	//!
	struct Status
	{
		// Robot side information
		enum {
			CALIBRATED   = 0x00000001,  //!< The robot has been calibrated
			FROZEN       = 0x00000002,  //!< The robot is frozen due to a watchdog timeout
			UNDERVOLTAGE = 0x00000010,  //!< The battery voltage is low
			RESERVED0    = 0x00000020,  //!< Reserved for future use
			KILLSWITCH   = 0x00000040,  //!< The kill switch has been tripped
			WATCHDOG     = 0x00000080   //!< The watchdog has timed out
		};

		int  hardware_flags;                                 //!< Hardware status flags
		bool joint_saturated[ NUM_LEGS ][ NUM_JOINTS ];      //!< Torque saturation detected
		bool joint_user_saturated[ NUM_LEGS ][ NUM_JOINTS ]; //!< User requested torque limit reached
		bool joint_at_limit[ NUM_LEGS ][ NUM_JOINTS ];       //!< Joint at motion limits

		// Host side information
		int  host_skipped_control;                       //!< Count of missed control cycles by the host computer
		int  host_skipped_recv;                          //!< Count of late packets received by the host from the robot
		int  host_queued_packets;                        //!< Count of unprocessed data packets from the robot

		Status(void);
	};

	//!
	//!  \struct  RobotReadInfo
	//!
	//!  \brief   Complete robot state as seen by the host computer
	//!
	//!  \_Description
	//!
	//!    RobotReadInfo contains all the information describing the state of the 
	//!    robot as seen by the host computer.  The contents include information 
	//!    about battery state, each joint, each leg and all of the sensors, 
	//!    current status of the communication system as well as robot health.
	//!
	struct RobotReadInfo
	{
		float     timestamp;        //!< timestamp for this data sample from the robot in seconds   
		Status    status;           //!< Status structure
		
		float     battery_voltage;  //!< battery voltage in volts
		float     battery_current;  //!< battery current in amps
		IMUInfo   imu;              //!< inertial measurement unit structure
		float     temperature;      //!< Robot CPU temperature in deg Celcius

		LegInfo   legs[ NUM_LEGS ]; //!< array of Leg info structures
		
		float     prox_sensor;      //!< proximity sensor measurment in meters 

		RobotReadInfo( void );      //!<  Default constructor
	};

	//!
	//!  \enum
	//!
	//!  \brief   Leg servo mode can be either joint independent or Cartesian.
	//!
	//!  \_Description
	//!
	//!    Leg servo type is either independent or cartesian and are
	//!    specified on a per leg basis.  A leg can only be set to Cartesian mode
	//!    if all the joints of that leg are already set to PD servo mode.
	//!
	enum LegServoCommandMode
	{
		LEG_MODE_INDEPENDENT,  //!< The joints are controlled independently
		LEG_MODE_CARTESIAN     //!< Leg joints are coordinated to servo the foot's cartesian position
	};

	//------------------------------------------------------------------
	//
	//  Mocap Enumerations and Data Types
	//
	//------------------------------------------------------------------

	//!
	//!  \enum
	//!
	//!  \brief   Mocap constants.
	//!
	enum MocapConstants
	{
		MAX_MARKERS  = 80,        //!<  Maximum number of markers supported by the mocap system
		MAX_TERRAINS =  8         //!<  Maximum number of terrain boards in a trial
	};

	//!
	//!  \struct  MocapBodyInfo
	//!
	//!  \brief   Motion capture information reported on a per-body basis.
	//!
	//!  \_Description
	//!
	//!    Each rigid section of the robot (e.g., trunk,
	//!    upper leg, lower leg, etc.) is referred to
	//!    as a 'body' by the motion capture system.
	//!
	//!    Body information includes a position, orientation, joint angles,
	//!    and an associated age variable that indicates
	//!    how many frames 'old' the data is.
	//!
	//!    The position and orientation are expressed in 
	//!    global coordinates of the motion capture system.
	//!
	//!    The joint angles correspond to the joint angles of the robot 
	//!    kinematics. The
	//!    angles are associated with the proximal joints of the body.
	//!    e.g. the joints of the upper leg body are the hip joints
	//!    of the front left leg.  Note that upperlegs are associated with 2-dof
	//!    joints while lower legs are associated with the 1-dof knees 
	//!    and only have a 'primary angle'.
	//!
	//!    The age indicates how many frames 'old' the data is.
	//!    A value of zero (0.0) indicates the data is up-to-date.
	//!    A positive number indicates time in the past.
	//!    A negative number (-1.0) indicates the body has
	//!    never been seen.
	//!
	struct MocapBodyInfo
	{
		bduVec3f   m_position;             //!<  3d vector of global position, in meters
		bduVec3f   m_orientation;          //!<  3d vector of orientation, in yaw-roll-pitch Euler angles
		float      m_primary_angle;        //!<  primary angle to attached body, in radians
		float      m_secondary_angle;      //!<  secondary angle (used only for hips), in radians
		int        m_age;                  //!<  number of frames this data is out of date
	};

	//!
	//!  \struct  MocapMarkerInfo
	//!
	//!  \brief   Motion capture on a per marker basis.
	//!
	//!  \_Description
	//!
	//!    The LittleDog robot has a number of reflective
	//!    markers on it that the motion capture system
	//!    uses to track the position of the robot and
	//!    calculate joint angles of the legs.
	//!
	//!    Each marker has a position and an 'age'.
	//!
	//!    The position is expressed in global coordinates of the
	//!    motion capture system.
	//!
	//!    The age indicates how many frames 'old' the data is.
	//!    A value of zero (0.0) indicates the data is up-to-date.
	//!    A positive number indicates time in the past.
	//!    A negative number (-1.0) indicates the marker has
	//!    never been seen.
	//!
	struct MocapMarkerInfo
	{
		bduVec3f   m_position;             //!<  3d vector of global position, in meters
		int        m_age;                  //!<  number of frames this data is out of date
	};

	//!
	//!  \struct  MocapReadInfo
	//!
	//!  \brief   Complete robot motion capture information.
	//!
	//!  \_Description
	//!
	//!    Includes all marker and body information as well
	//!    as marker count and current frame number for the
	//!    robot only.
	//!
	struct MocapReadInfo
	{
		MocapBodyInfo     m_bodies[ NUM_BODIES ];    //!<  mocap body info
		MocapMarkerInfo   m_markers[ MAX_MARKERS ];  //!<  mocap marker info
		int               m_frame_number;            //!<  latest frame number from the mocap system
		int               m_num_markers;             //!<  number of mocap markers being tracked
		float             m_timestamp;               //!<  timestamp on host when data arrived

		MocapReadInfo( void );
	};
  
	//!
	//!  \enum
	//!
	//!  \brief   Terrain board ID constants.
	//!
	enum TerrainID
	{
		TERRAIN_NONE = 0,             //!<  No terrain board visible
		TERRAIN_A,                    //!<  Terrain board A
		TERRAIN_B,                    //!<  Terrain board B
		TERRAIN_C,                    //!<  Terrain board C
		TERRAIN_D,                    //!<  Terrain board D
		TERRAIN_E,                    //!<  Terrain board E
		TERRAIN_F,                    //!<  Terrain board F
		TERRAIN_G,                    //!<  Terrain board G
		TERRAIN_H,                    //!<  Terrain board H
		TERRAIN_I,                    //!<  Terrain board I
		TERRAIN_J,                    //!<  Terrain board J
		TERRAIN_K,                    //!<  Terrain board K
		TERRAIN_L,                    //!<  Terrain board L
		TERRAIN_M,                    //!<  Terrain board M
		TERRAIN_N,                    //!<  Terrain board N
		TERRAIN_O,                    //!<  Terrain board O
		TERRAIN_P,                    //!<  Terrain board P
		TERRAIN_Q,                    //!<  Terrain board Q
		TERRAIN_R,                    //!<  Terrain board R
		TERRAIN_S,                    //!<  Terrain board S
		TERRAIN_T,                    //!<  Terrain board T
		TERRAIN_U,                    //!<  Terrain board U
		TERRAIN_V,                    //!<  Terrain board V
		TERRAIN_W,                    //!<  Terrain board W
		TERRAIN_X,                    //!<  Terrain board X
		TERRAIN_Y,                    //!<  Terrain board Y
		TERRAIN_Z,                    //!<  Terrain board Z
		TERRAIN_AA,                   //!<  Terrain board AA
		TERRAIN_AB,                   //!<  Terrain board AB
		TERRAIN_AC,                   //!<  Terrain board AC
		TERRAIN_AD,                   //!<  Terrain board AD
		TERRAIN_AE,                   //!<  Terrain board AE
		TERRAIN_AF,                   //!<  Terrain board AF
		TERRAIN_AG,                   //!<  Terrain board AG
		TERRAIN_AH,                   //!<  Terrain board AH
		TERRAIN_AI,                   //!<  Terrain board AI
		TERRAIN_AJ,                   //!<  Terrain board AJ
		TERRAIN_AK,                   //!<  Terrain board AK
		TERRAIN_AL,                   //!<  Terrain board AL
		TERRAIN_AM,                   //!<  Terrain board AM
		TERRAIN_AN,                   //!<  Terrain board AN
		TERRAIN_AO,                   //!<  Terrain board AO
		TERRAIN_AP,                   //!<  Terrain board AP
		TERRAIN_AQ,                   //!<  Terrain board AQ
		TERRAIN_AR,                   //!<  Terrain board AR
		TERRAIN_AS,                   //!<  Terrain board AS
		TERRAIN_AT,                   //!<  Terrain board AT
		TERRAIN_AU,                   //!<  Terrain board AU
		TERRAIN_AV,                   //!<  Terrain board AV
		TERRAIN_AW,                   //!<  Terrain board AW
		TERRAIN_AX,                   //!<  Terrain board AX
		TERRAIN_AY,                   //!<  Terrain board AY
		TERRAIN_AZ,                   //!<  Terrain board AZ

                // ... Note: NUM_TERRAINS covers up to TERRAIN_ZZ

		NUM_TERRAINS = 703            //!<  Total number of possible terrain boards
                                              //!<  1 (for TERRAIN_NONE) + 26*27 = 703
	};

	//!
	//!  \struct  TerrainInfo
	//!
	//!  \brief   Motion capture information on a per terrain board basis.
	//!
	//!  \_Description
	//!
	//!    Each terrain board that is tracked has an ID tag,  position, orientation and
	//!    an associated age variable that indicates
	//!    how many frames 'old' the data is.
	//!
	struct TerrainInfo
	{
		enum TerrainID m_id;           //!<  ID of this terrain board
		bduVec3f       m_position;     //!<  3d position in global coordinates, in meters
		bduVec3f       m_orientation;  //!<  3d position in global coordinates, in yaw-roll-pitch Euler angles
		int            m_age;          //!<  number of frames this data is out of date
	};

	//!
	//!  \struct  TerrainReadInfo
	//!
	//!  \brief   Complete terrain motion capture information.
	//!
	//!    Includes all body information for all terrain boards
	//!    present in the current trial.
	//!
	struct TerrainReadInfo
	{
		TerrainInfo m_terrains[ MAX_TERRAINS ];
	};

  
 /////////////////////////////////////////////////////////////////////////
 //!
 //! @name    LittleDog Public Interface
 //!
 //!    Following is the public interface of the LittleDog class.
 //!    By calling these functions the user can begin a session,
 //!    move it through the various initialization steps, run
 //!    trials, and shut down the session.
 //!
 //! @{
 /////////////////////////////////////////////////////////////////////////

 public:

	//!
	//!  \brief    Primary constructor.
	//!
	//!  \_Description
	//!
	//!    This is the public constructor of the LittleDog class.
	//!
	//!  \_State_Information
	//!
	//!    The initial system state will be STATE_UNINITIALIZED.
	//!
	//!    Call initializeRobot() to advance the system state.
	//!
	LittleDog();

	//!
	//!  \brief    Virtual destructor.
	//!
	//!  \_Description
	//!
	//!    This is the virtual destructor of the LittleDog class.
	//!
	virtual ~LittleDog();
  
	//!
	//!  \brief    Query system for current state.
	//!
	//!  \_Description
	//!
	//!    This function returns the current state of the system.
	//!    Different function calls can be made from any given
	//!    state.
	//!
	//!  \_State_Information
	//!
	//!    This function can be called when the system is in any
	//!    state.  It does not change the state of the system.
	//!
	//!  \_Returns State enumeration result
	//!
	//!    Returned result will be one of:
	//!
	//!      - STATE_UNINITIALIZED       -  initial state after object construction
	//!      - STATE_ROBOT_INITIALIZED   -  robot is initialized
	//!      - STATE_ROBOT_CALIBRATING   -  robot is calibrating
	//!      - STATE_ROBOT_CALIBRATED    -  robot is calibrated and ready to run a trial
	//!      - STATE_TRIAL_PLANNING      -  trial is in the planning phase
	//!      - STATE_TRIAL_RUNNING       -  trial is in the running phase
	//!      - STATE_STOPPING_ROBOT      -  system is attempting to stop calibration or trial run
	//!      - STATE_ABORTING_CONNECTION -  system has abortied the connection to the robot (this is a terminal state)
	//!      - STATE_STOPPING            -  system is returning to STATE_ROBOT_INITIALIZED or
	//!                                     STATE_ROBOT_CALIBRATED, depending on last getState()
	//!
	State getState( void ) const;

	//!
	//!  \brief    Query system for the performer name, as derived from bdi_rt.cfg
	//!
	//!  \_Description
	//!
	//!    This function returns the processed user.performer.name found
	//!    in the bdi_rt.cfg file.  If the name or the config file was
	//!    missing, the LittleDog constructor will exit(EXIT_FAILURE)
	//!    the program after printing out instructions.  The performer
	//!    name should have no spaces in it and be less than 64 characters
	//!    in length.  It will be truncated to 64 characters and spaces
	//!    will be converted to underscores ('_').  This name will be
	//!    used when creating LittleDog dataset files.
	//!
	//!  \_Returns     NULL terminated string
	//!
	const char* getPerformerName( void ) const;

	//!
	//!  \brief    Query system for the performer name, as derived from bdi_rt.cfg
	//!
	//!  \_Description
	//!
	//!    This function returns the processed user.vicon.modelname found
	//!    in the bdi_rt.cfg file.  If the name or the config file was
	//!    missing, the LittleDog constructor will exit(EXIT_FAILURE)
	//!    the program after printing out instructions.  This name will
	//!    be used when connecting to the mocap system.
	//!
	//!  \_Returns     NULL terminated string
	//!
	const char* getMocapModelName( void ) const;

	//!
	//!  \brief    Convert an LD_ERROR_CODE to a string
	//!
	//!  \_Description
	//!
	//!    This function returns a NULL terminated string that
	//!    matches an LD_ERROR_CODE.  For example, passing
	//!    LD_TIMED_OUT will return the string "Operation timed out".
	//!    This function is guaranteed  never to return NULL.
	//!    If an invalid error code is given, a string to that 
	//!    effect is returned instead of NULL.
	//!
	//!  \_Parameters
	//!
	//!    \_in   error_code    - error code for which to get string
	//!
	//!  \_Returns     NULL terminated string
	//!
	const char* getErrorCodeString( LD_ERROR_CODE error_code ) const ;

	//!
	//!  \brief    Abort all communication with the robot.
	//!
	//!  \_Description
	//!
	//!    This function shuts down the communications.  This means that:
	//!      - any running trial will be aborted
	//!      - all communications with the robot will be disconnected
	//!      - all communications with the vicon will be disconnected
	//!      - the state will transition through STATE_ABORTING_CONNECTION
	//!      - the state will end up in STATE_UNINITIALIZED
	//!
	//!    This should \e not be the standard method of ending
	//!    trials or shutting down the system.  This function is
	//!    intended for use only in exceptional circumstances.
	//!    For a more graceful termination of a trial use stopRobot().
	//!
	//!  \_State_Information
	//!
	//!    This function can be called when the system is in any
	//!    state except STATE_UNINITIALIZED or
	//!    STATE_ABORTING_CONNECTION.  It is legal to call this
	//!    function from within updateControl().
	//!
	//!    This function does \e not block, but returns immediately.
	//!
	//!  \_Returns LD_ERROR_CODE result
	//!
	//!    Returned result will be one of:
	//!
	//!      - LD_OKAY              - abort succeeded
	//!
	LD_ERROR_CODE abortConnection( void );

	//!
	//!  \brief    Initiate communication with the mocap system.
	//!
	//!  \_Description
	//!
	//!    This function initializes communication with the
	//!    motion capture system.  Use of the motion capture system is
	//!    optional, but once initialized the motion capture system
	//!    will be used until stopMocap() or abortConnection() is called.
	//!
	//!  \_State_Information
	//!
	//!    This function can be called only when the system is in one
	//!    of the following states:
	//!
	//!       STATE_UNINITIALIZED, STATE_ROBOT_INITIALIZED, STATE_ROBOT_CALIBRATED
	//!
	//!    If this function succeeds, the call returns LD_OKAY and system state
	//!    remains in its current state.  The return value for getMocapInitialized()
	//!    will be true on success.  See stopMocap() for closing the mocap connection.
	//!
	//!    If this function fails, the system state will remain as it was before
	//!    this function call.
	//!
	//!    This function will block until connection is established or
	//!    until it is determined that the connection can not be established.
	//!
	//!  \_Returns LD_ERROR_CODE result
	//!
	//!    Returned result will be one of:
	//!
	//!      - LD_OKAY                     - connection succeeded
	//!      - LD_MOCAP_COMMS_FAILED       - connection failed
	//!      - LD_INVALID_REQUEST          - request is out of context 
	//!
	LD_ERROR_CODE initializeMocap( void );


	//!
	//!  \brief    Closes communications with the mocap system.
	//!
	//!  \_Description
	//!
	//!    This function shuts down communication with the
	//!    motion capture system.  Use of the motion capture system is
	//!    optional, but once initialized the motion capture system
	//!    will be used until stopMocap() or abortConnection() is called.
	//!
	//!  \_State_Information
	//!
	//!    This function can be called only when the system is in one
	//!    of the following states:
	//!
	//!       STATE_UNINITIALIZED, STATE_ROBOT_INITIALIZED, STATE_ROBOT_CALIBRATED
	//!
	//!    If this function succeeds, the call returns LD_OKAY and system state
	//!    remains in its current state.  The return value for getMocapInitialized()
	//!    will be false on success.  See initializeMocap() for starting the mocap connection.
	//!
	//!    If this function fails, the system state will remain as it was before
	//!    this function call.
	//!
	//!    This function will block until the connection has been shut down.
	//!
	//!  \_Returns LD_ERROR_CODE result
	//!
	//!    Returned result will be one of:
	//!
	//!      - LD_OKAY                     - connection shutdown succeeded
	//!      - LD_INVALID_REQUEST          - request is out of context 
	//!
	LD_ERROR_CODE stopMocap( void );

	//!
	//!  \brief    Initiate the communication with the robot.
	//!
	//!  \_Description
	//!
	//!    This function initializes communication with the
	//!    robot.
	//!
	//!  \_State_Information
	//!
	//!    This function can be called when the system is in
	//!    STATE_UNINITIALIZED.
	//!
	//!    If this function succeeds, the system state will become
	//!    STATE_ROBOT_INITIALIZED.
	//!
	//!    If this function fails, the system state will remain
	//!    in its current state (STATE_UNINITIALIZED).
	//!
	//!    Call calibrate() to advance the system state.
	//!
	//!  \_Returns LD_ERROR_CODE result
	//!
	//!    Returned result will be one of:
	//!
	//!      - LD_OKAY                     - connection succeeded
	//!      - LD_INVALID_REQUEST          - connection is aborting or operation is stopping
	//!      - LD_TIMEOUT                  - timed out trying to connect
	//!      - LD_ROBOTS_COMMS_FAILED      - comm errors trying to connect
	//!      - LD_INVALID_REQUEST          - request is out of context (connection aborting or robot stopping)
	//!
	LD_ERROR_CODE initializeRobot( void );

	//!
	//!  \brief    Initiate the robot calibration procedure.
	//!
	//!  \_Description
	//!
	//!    This function puts the system into a calibration phase.
	//!
	//!    The robot can be calibrated multiple times during a
	//!    session, between trials.
	//!
	//!    This function will block until calibration succeeds
	//!    or fails.  The calibration can be interrupted by
	//!    an asynchronous call to stopRobot() or abortConnection().
	//!    It is recommended that calibration happens in a
	//!    thread separate from any UI thread so that the
	//!    calibration can be stopped or aborted if necessary.
	//!
	//!    If input argument "force" is false and the robot has already
	//!    been calibrated, then this function returns immediately.  
	//!    If force is true, calibration will always occur.
	//!
	//!    If a calibration occurs, a final posture will be taken by
	//!    the robot post-calibration.  That posture is either a standing
	//!    posture (the zero position) or a crawl position.  Input
	//!    argument "crawl" dictates which posture to take post-calibration.
	//!
	//!  \_State_Information
	//!
	//!    This function can be called when the system is in
	//!    STATE_ROBOT_INITIALIZED or STATE_CALIBRATED.
	//!
	//!    While this function is running the system state will be
	//!    STATE_ROBOT_CALIBRATING.
	//!
	//!    If this function succeeds, the system state will become
	//!    STATE_ROBOT_CALIBRATED.
	//!
	//!    If this function fails, the system state will become
	//!    STATE_ROBOT_INITIALIZED.
	//!
	//!    Call runTrial() or calibrate() again to advance the system
	//!    state.
	//!
	//!  \_Parameters
	//!
	//!    \_in  force - forces recalibration if the robot has already been calibrated
	//!    \_in  crawl - post calibration, the robot assumes a crawl (if true) or stand (if false) posture.
	//!
	//!  \_Returns LD_ERROR_CODE result
	//!
	//!    Returned result will be one of:
	//!
	//!      - LD_OKAY                     - calibration succeeded
	//!      - LD_INVALID_REQUEST          - request is out of context (connection aborting or robot stopping)
	//!      - LD_TIMED_OUT                - calibration timed out
	//!      - LD_TERMINATED               - calibration terminated before completion 
	//!      - LD_ROBOT_STOP_IS_BROKEN     - the robot needs repair.  Cease using it.  The robot
	//!                                      should be packaged up and shipped to Boston Dynamics.
	LD_ERROR_CODE calibrate( bool force=false, bool crawl=true );

	//!
	//!  \brief    Begin a robot trial planning phase.
	//!
	//!  \_Description
	//!
	//!    This function begins a robot trial run during which user
	//!    control functions will be called and data will be collected.
	//!    The planning phase of the trial is entered automatically, and is to be
	//!    followed by the trial running phase, triggered by a user call to donePlanning().
	//!
	//!    This function will block until the trial is stopped by
	//!    calling stopRobot(), or aborted by calling abortConnection().
	//!
	//!    It is recommended that trial runs happen in a thread
	//!    separate from any UI thread, so that the trial can
	//!    be stopped or aborted if necessary.  Trials can be
	//!    stopped asychronously by calling stopRobot() or
	//!    abortConnection().
	//!
	//!    Multiple trials can be run during a session.  i.e. runTrial()
	//!    can be called multiple times before the LittleDog object
	//!    is destroyed.
	//!
	//!    Calling the calibrate() function between successive trials
	//!    is optional.
	//!
	//!    When the trial is stopped by a call to stopRobot(), data
	//!    will be saved before this call returns.  Data for both the planning 
	//!    and running phase of the trial will be saved in the produced dataset.
	//!
	//!  \_State_Information
	//!
	//!    This function can be called only when the system is in
	//!    STATE_CALIBRATED.
	//!
	//!    While the trial is running the system state will be either
	//!    STATE_TRIAL_PLANNING or STATE_TRIAL_RUNNING, depending upon whether
	//!    the user has called donePlanning() yet or not.  In both of these
	//!    states the user has access to essentially the same LittleDog API calls.
	//!
	//!    The robot state will initially enter into STATE_TRIAL_PLANNING and remain there 
	//!    until the user calls donePlanning() or the trial stops.
	//!    If donePlanning() is called successfully, the state will become STATE_TRIAL_RUNNING.
	//!    Once donePlanning() is called, it can not be called successfully again until the 
	//!    next trial begins.
	//!
	//!    When this function returns after a successful trial the system state will be
	//!    STATE_ROBOT_INITIALIZED.
	//!
	//!  \_Returns LD_ERROR_CODE result
	//!
	//!    Returned result will be one of:
	//!
	//!      - LD_OKAY                       - success
	//!      - LD_INVALID_REQUEST            - request is out of context (connection aborting or robot stopping)
	//!      - LD_ROBOT_COMMS_UNINITIALIZED  - communications not yet initialized
	//!      - LD_USER_INIT_CONTROL_FAILED   - user initialization failed
	//!      - LD_TERMINATED                 - trial terminated 
	//!
	LD_ERROR_CODE runTrial( void );

	//!
	//!  \brief    Begin robot trial running phase.
	//!
	//!  \_Description
	//!
	//!    This function ends the planning phase of a trial and forwards the LittleDog
	//!    state machine into the running phase of a trial.  During the running phase 
	//!    (as in the planning phase) the user control functions will be called and data
	//!    will continue to be collected.
	//!
	//!    It is legal to call this function from within updateControl().
	//!
	//!  \_State_Information
	//!
	//!    This function can be called only when the system is in STATE_TRIAL_PLANNING.
	//!
	//!    While the trial is in the running phase, the system state will be STATE_TRIAL_RUNNING.
	//!
	//!    This function does \e not block, but returns immediately.
	//!    If LD_OKAY is returned, then the state should become STATE_TRIAL_RUNNING.
	//!
	//!  \_Returns LD_ERROR_CODE result
	//!
	//!    Returned result will be one of:
	//!
	//!      - LD_OKAY                       - success, STATE_TRIAL_RUNNING begins
	//!      - LD_INVALID_REQUEST            - can not stop planning from the current state
	//!
	LD_ERROR_CODE donePlanning( void );

	//!
	//!  \brief    Stop current operation being performed on robot.
	//!
	//!  \_Description
	//!
	//!    This function attempts to stop the robot from moving.
	//!    There are two times during which the LittleDog software
	//!    will attempt to actuate robot joints: during the calibration
	//!    procedure, and during a trial.
	//!
	//!    This function is the preferred way to prematurely stop a
	//!    calibration or a trial.
	//!
	//!    \sa  abortConnection().
	//!
	//!  \_State_Information
	//!
	//!    This function can be called asychronously during calls
	//!    to calibrate() or runTrial() when the state is
	//!    STATE_ROBOT_CALIBRATING, STATE_TRIAL_PLANNING, or STATE_TRIAL_RUNNING.  It
	//!    is legal to call this function from within updateControl().
	//!
	//!    This function does \e not block, but returns immediately.
	//!    The stop is complete when the system state changes from
	//!    STATE_ROBOT_CALIBRATING to STATE_ROBOT_INITIALIZED, or from
	//!    STATE_TRIAL_PLANNING to STATE_ROBOT_CALIBRATED, or from
	//!    STATE_TRIAL_RUNNING to STATE_ROBOT_CALIBRATED.
	//! 
	//!  \_Returns LD_ERROR_CODE result
	//!
	//!    Returned result will be one of:
	//!
	//!      - LD_OKAY             - request to stop successfully posted
	//!      - LD_INVALID_REQUEST  - can not stop from the current state
	//!
	LD_ERROR_CODE stopRobot( void );


	//------------------------------------------------------------------
	//! @}
	//! @name    Status Get Functions
	//!
	//!    Accessor functions for retrieving robot state and status.
	//!
	//!    Unless otherwise specified, these functions should be called
	//!    only after communication has been established with the robot.
	//!    (i.e., the system state has passed through
	//!    STATE_ROBOT_INITIALIZED.)
	//!
	//! @{
	//------------------------------------------------------------------

	//!
	//!  \brief    Get the timestamp from the robot.
	//!
	//!  \_Description
	//!
	//!    This function returns the timestamp of the last
	//!    communication with the robot via the passed pointer.
	//!
	//!  \_Parameters
	//!
	//!    \_out  timestamp - time measured in seconds
	//!
	//!  \_Returns         LD_ERROR_CODE result
	//!
	//!     LD_INVALID_REQUEST
	//!     LD_OKAY
	//!
	LD_ERROR_CODE getDataTimestamp( float * timestamp) const;

	//!
	//!  \brief    Get status information from the robot.
	//!
	//!  \_Description
	//!
	//!    This funtion fills in an empty status structure provided by
	//!    the caller.  The Status struct includes (on a per leg/joint
	//!    basis) joint saturated and joint at limit status.  The
	//!    Status struct also includes counts of the number of host
	//!    skipped control packets and host queued packets.
	//!
	//!  \_Parameters
	//!
	//!    \_out s - a structure that provides several status indicators
	//!
	//!  \_Returns         LD_ERROR_CODE result
	//!
	//!     LD_INVALID_REQUEST
	//!     LD_OKAY
	//!
	LD_ERROR_CODE getStatus( Status * s ) const;

	//!
	//!  \brief    Get the robot battery voltage.
	//!
	//!  \_Parameters
	//!
	//!    \_out  voltage - a float for the battery voltage (in volts)
	//!
	//!  \_Returns         LD_ERROR_CODE result
	//!
	//!     LD_INVALID_REQUEST
	//!     LD_OKAY
	//!
	LD_ERROR_CODE getBatteryVoltage( float * voltage ) const;

	//!
	//!  \brief    Get the robot battery current.
	//!
	//!  \_Parameters
	//!
	//!    \_out  current - a float for the battery current (in Amps)
	//!
	//!  \_Returns         LD_ERROR_CODE result
	//!
	//!     LD_INVALID_REQUEST
	//!     LD_OKAY
	//!
	LD_ERROR_CODE getBatteryCurrent( float * current ) const;

	//!
	//!  \brief    Get the robot CPU temperature.
	//!
	//!  \_Parameters
	//!
	//!    \_out  voltage - a float pointer for the CPU temperature 
	//!                     (in deg C)
	//!
	//!  \_Returns         LD_ERROR_CODE result
	//!
	//!     LD_INVALID_REQUEST
	//!     LD_OKAY
	//!
	LD_ERROR_CODE getCPUTemperature( float * voltage ) const;

	//!
	//!  \brief    Check to see if mocap is initialized.  
	//!
	//!            This method can be called in any state.
	//!
	//!  \_Returns TRUE if mocap initialized, otherwise FALSE
	//!
	bool getMocapInitialized( void ) const;

	//!
	//!  \brief    Get the percentage of lost packets from robot.
	//!
	//!  \_Description
	//!
	//!    This function returns the computed average packet loss of robot-source
	//!    data over a running time window of 1 second.  It returns this
	//!    pre-computed data immediately.
	//!
	//!  \_Parameters
	//!
	//!    \_out  percentage - a float pointer for the percentage package loss.
	//!
	//!  \_Returns         LD_ERROR_CODE result
	//!
	//!     LD_INVALID_REQUEST
	//!     LD_OKAY
	//!
	LD_ERROR_CODE getPacketLoss( float * percentage ) const;

	//!
	//!  \brief    Indicates whether the motion capture system is functioning.
	//!
	//!  \_Description
	//!
	//!    This function returns a boolean indicating whether at least one
	//!    packet of motion capture data has arrived within the last second.
	//!
	//!  \_Returns         
	//!
	//!      - true if data is being received from the mocap system
	//!      - false if no data is being received from the mocap system
	//!
	bool getMocapDataReceived( void ) const;

	//------------------------------------------------------------------
	//! @}
	//! @name    Kinematic Transformations
	//!
	//!    Coordinate systems are attached to each body on LittleDog.
	//!    Coordinate systems are available fixed to the trunk, to the
	//!    upper leg, and to the lower leg.  A global coordinate system
	//!    is also defined by the motion capture system. The LittleDog API provides
	//!    transformations to convert back and forth between these
	//!    various coordinate systems.  
	//!
	//!    When the robot's legs are perpendicular to the ground, all the
	//!    link coordinate systems are aligned and the robot is said to
	//!    be in the "reference" configuration. 
	//!
	//!    All coordinate systems have the x-axis pointing forward, the 
	//!    y axis pointing to the left, and the z axis completing the right 
	//!    hand set.  The origin of the robot is fixed to the center of the 
	//!    robot's hips.  For all other bodies on the robot, (ie. upper/lower 
	//!    legs, etc), the origin of the frame of reference is at each body's 
	//!    inboard joint. 
	//!
	//!    Many functions in this section make use of unit quaternion 4
	//!    vectors to specify orientations.  These vectors are stored in
	//!    vector-scalar format (qx, qy, qz, w). Many of the underlying
	//!    structures contain orientation in yaw, roll, pitch Euler
	//!    coordinates.
	//!
	//! @{
	//------------------------------------------------------------------

	//!
	//!  \brief    Compute forward kinematics of a foot.
	//!
	//!  \_Description
	//!
	//!    This function computes the forward kinematics of
	//!    the foot associated with the selected leg.  Checks for 
	//!    reachability are not performed.  Angle range validation is
	//!    not performed.
	//!
	//!  \_Parameters
	//!
	//!    \_in   leg    - leg index
	//!    \_in   q      - 3 vector of joint angles
	//!    \_out  pos    - 3 vector describing foot position in the robot body coordinate system
	//!    \_out  orient - unit quaternion representing the orientation of the foot (vector, scaler) w.r.t. the robot body
	//!
	//!    leg angles (parameter q) are in order hip_x, hip_y and knee
	//!
	//!  \_Returns        LD_ERROR_CODE result
	//!
	//!    - LD_OKAY
	//!    - LD_BAD_LEG_INDEX
	//!    - LD_BAD_QUATERNION_VECTOR
	//!    - LD_INSUFFICIENT_ARGUMENTS
	//!
	LD_ERROR_CODE computeFootFK( LegIndex leg,
		const bduVec3f & q, 
		bduVec3f * pos,
		bduVec4f * orient );

	//!
	//!  \brief    Compute inverse kinematics of the foot associated with the selected leg.
	//!
	//!  \_Description
	//!
	//!    This function computes inverse kinematics of the foot
	//!    associated with the selected leg.
	//!
	//!    The leg index, as well as the target position, pos, for 
	//!    the foot expressed relative to the robot body coordinate 
	//!    system.  On success, the vectors q_bent_forward and
	//!    q_bent_rearward return the joint angles required to make
	//!    the center of the foot reach the specified location.
	//!
	//!  \_Parameters
	//!
	//!    \_in   leg             - leg index
	//!    \_in   pos             - 3 vector for desired foot position relative to the body
	//!    \_out  q_bent_forward  - joint angle 3 vector for 'forward knee' solution
	//!    \_out  q_bent_rearward - joint angle 3 vector for 'backward knee' solution
	//!
	//!    leg angles (q vector) are in hip_x, hip_y and knee order.
	//!
	//!  \_Returns        LD_ERROR_CODE result
	//!
	//!    Returned result will be one of:
	//!
	//!      - LD_OKAY
	//!      - LD_BAD_LEG_INDEX
	//!      - LD_INSUFFICIENT_ARGUMENTS
	//!      - LD_CANNOT_REACH_IK_KNEE_BENT_FORWARD
	//!      - LD_CANNOT_REACH_IK_KNEE_BENT_REARWARD
	//!      - LD_CANNOT_REACH_IK_KNEE_BENT_EITHER
	//!
	LD_ERROR_CODE computeFootIK( LegIndex leg,
		const bduVec3f & pos, 
		bduVec3f * q_bent_forward,
		bduVec3f * q_bent_rearward );

	//!
	//!  \brief    Compute the (point) differential kinematics of a foot.
	//!
	//!  \_Description
	//!
	//!     This function computes the (point) differential kinematics
	//!     of the foot associated with the selected leg.
	//!
	//!  \_Parameters
	//!
	//!    \_in   leg             - leg index
	//!    \_in   q               - 3 vector of joint angles
	//!    \_out  jac             - on success, contains the Jacobian expressed in robot body coordinates
	//!
	//!  \_Returns        LD_ERROR_CODE result
	//!
	//!    Returned result will be one of:
	//!
	//!      - LD_OKAY
	//!      - LD_BAD_LEG_INDEX
	//!      - LD_INSUFFICIENT_ARGUMENTS
	//!
	LD_ERROR_CODE computeFootDK( LegIndex leg,
		const bduVec3f & q,
		bduMat3f * jac );

	//!
	//!  \brief    Transform a position and/or orientation.
	//!
	//!  \_Description
	//!
	//!    This function transforms optional inputs pos0 and orient0 
	//!    according to the transform represented by required inputs
	//!    t_pos and t_orient.  The results are stored in option outputs
	//!    pos1 and orient1.
	//!
	//!  \_Parameters
	//!
	//!    \_in   t_pos      - translational component of the transform to be applied
	//!    \_in   t_orient   - rotational component of the transform to be applied
	//!    \_in   pos0       - position of initial frame        (optional)
	//!    \_in   orient0    - orientation of initial frame     (optional)
	//!    \_out  pos1       - position of transformed frame    (optional)
	//!    \_out  orient1    - orientation of transformed frame (optional)
	//!
	//!    Any of pos0, orient0, pos1, or orient1, can be NULL (unlike
	//!    most functions in the API pointers are used here for some inputs).
	//!    A NULL parameter will not be used or filled out, making it possible
	//!    to just transform points or orientations.    
	//!
	//!  \_Returns        LD_ERROR_CODE result
	//!
	//!    Returned result will be one of:
	//!
	//!      - LD_OKAY
	//!      - LD_BAD_QUATERNION_VECTOR
	//!
	LD_ERROR_CODE transform(
		const bduVec3f & t_pos,
		const bduVec4f & t_orient, 
		const bduVec3f * pos0,
		const bduVec4f * orient0, 
		bduVec3f * pos1,
		bduVec4f * orient1 ); 

	//!
	//!  \brief    Compute an inverse transform.
	//!
	//!  \_Description
	//!
	//!    The inverse of the transform represented by pos0 and 
	//!    orient0 is computed and the resulting transform is stored 
	//!    in pos1 and orient1.
	//!
	//!  \_Parameters
	//!
	//!    \_in   pos0         - translational component of the transform
	//!    \_in   orient0      - rotational component of the transform
	//!    \_out  pos1         - translational component of the inverse transform
	//!    \_out  orient1      - rotational component of the inverse transform
	//!              
	//!    All pointer parameters must be non NULL.
	//!              
	//!  \_Returns        LD_ERROR_CODE result
	//!
	//!    Returned result will be one of:
	//!
	//!      - LD_OKAY
	//!      - LD_INSUFFICIENT_ARGUMENTS
	//!      - LD_BAD_QUATERNION_VECTOR
	//!
	LD_ERROR_CODE invTransform( const bduVec3f & pos0,
		const bduVec4f & orient0, 
		bduVec3f * pos1,
		bduVec4f * orient1 ); 

	//!
	//!  \brief    Get the skeletal joint offsets
	//!
	//!  \_Description
	//!
	//!    Allows the caller to get the joint offset programatically instead
	//!    of consulting documentation or hard-coding the values.  These
	//!    skeletal offsets never change and are measured with respect to
	//!    the joint's parent.  The body is the root of the kinematic tree, 
	//!    and is followed by hip rx, hip ry, knee ry, and foot.  One
	//!    way to visualize this is to consider the robot standing in the
	//!    zero configuration (i.e. all legs perpendicular to the ground).
	//!    In this configuration these offsets describe the location of the
	//!    frames' origins with respect to its parent
	//!
	//!  \_Parameters
	//!
	//!    \_in   leg          - leg index
	//!    \_in   joint_offset - joint offset type
	//!    \_out  offset       - a vector (x,y,z), in meters, describing
	//!                          the offset for this joint type relative to the
	//!                          joint's parent, expressed in the body coordinate 
	//!                          system, when the robot is in the zero configuration.
	//!
	//!  \_Returns         LD_ERROR_CODE result
	//!
	LD_ERROR_CODE getJointOffset( LegIndex leg, 
		JointOffset joint_offset, 
		bduVec3f * offset ) const;

	//------------------------------------------------------------------
	//! @}
	//! @name    Get Functions
	//!
	//!    Unless otherwise specified, these functions should be called
	//!    only when the system is at least in STATE_ROBOT_INITIALIZED.
	//!
	//! @{
	//------------------------------------------------------------------

	//!
	//!  \brief    Read the joint state for the specified leg and joint.
	//!
	//!  \_Parameters
	//!
	//!    \_in   leg    - leg index
	//!    \_in   joint  - joint index
	//!    \_out  angle  - current joint angle
	//!    \_out  vel    - current joint velocity (optional)
	//!    \_out  tau_d  - current joint torque command (optional)
	//!
	//!    If optional output variables are NULL pointers the information 
	//!    will not be returned.
	//!
	//!  \_Returns         LD_ERROR_CODE result
	//!
	LD_ERROR_CODE getJointState( LegIndex leg,
		JointIndex joint, 
		float * angle,
		float * vel=0,
		float * tau_d=0 ) const;

	//!
	//!  \brief    Get the position of the robot feet in robot coordinates.
	//!
	//!  \_Description
	//!
	//!    Read the leg angles and call computeFootFK() to determine the position and orientation.
	//!
	//!  \_Parameters
	//!
	//!    \_in   leg    - leg index
	//!    \_out  pos    - 3d position of foot in robot coordinates
	//!    \_out  orient - unit quaternion representing the orientation of the foot in robot coordinates
	//!
	//!  \_Returns         LD_ERROR_CODE result
	//!
	LD_ERROR_CODE getFootPosition( LegIndex leg, bduVec3f* pos, bduVec4f* orient );

	//!
	//!  \brief    Get the force on a single leg
	//!
	//!  \_Description
	//!
	//!    Read the foot contact sensor for the specified leg.
	//!
	//!  \_Parameters
	//!
	//!    \_in   leg    - leg index
	//!    \_out  force  - returned force in Newtons in [x,y,z] order
	//!
	//!  \_Returns         LD_ERROR_CODE result
	//!
	LD_ERROR_CODE getFootForce( LegIndex leg, bduVec3f* force ) const;
	
	//!
	//!
	//!  \brief   Get readings from the onboard inertia measurement unit (IMU)
	//!
	//!  \_Description
	//!
	//!    This function returns the readings from the onboard inertial measurement unit
	//!    (IMU) for orientation, angular rate and acceleration.
	//!    The information is returned via the passed vector pointers.
	//!    
	//!    The orientation of the robot is reported as a unit quaternion in a 4-vector
	//!    (qx, qy, qz, w).
	//!
	//!    Angle rates are in rad/sec and are instantaneous rates about the body X, Y,
	//!    and Z axes.  The vector order is [X rate, Y rate, Z rate] in radians/sec.
	//!
	//!    Accelerations are in the body frame and reported as [X, Y, Z] and in 
	//!    meters/sec^2. The measurements include gravity, so when the robot rests
	//!    on a table, the vector reported will be [0, 0, 9.8] m/sec^2.
	//!
	//!  \_Parameters
	//!
	//!    \_out  orient    - a unit quaternion orientation (vector, scalar)
	//!    \_out  angrates  - a three dimensional angular rate (in rad/sec)   
	//!    \_out  acc       - a three dimensional acceleration (in m/sec^2)
	//!
	//!  \_Returns         LD_ERROR_CODE result
	//!
	LD_ERROR_CODE getIMUReadings( bduVec4f * orient,
		bduVec3f * angrates,
		bduVec3f * acc ) const;

	//!
	//!  \brief    Get the proximity sensor value from the robot
	//!
	//!  \_Description
	//!
	//!    This function returns the value from the proximity sensor via
	//!    the passed pointer.
	//!
	//!  \_Parameters
	//!
	//!    \_out  proximity - a distance measure, in meters
	//!
	//!  \_Returns         LD_ERROR_CODE result
	//!
	LD_ERROR_CODE getProximityReading( float * proximity) const;

	//!
	//!  \brief    Get the latest sensor readings from the robot.
	//!
	//!  \_Parameters
	//!
	//!    \_out  read_info  - latest unprocessed robot sensor readings
	//!
	//!  \_Returns         LD_ERROR_CODE result
	//!
	LD_ERROR_CODE getBulkRobotInputs(RobotReadInfo* read_info) const;

	//------------------------------------------------------------------
	//! @}
	//! @name    Motion Capture System Functions
	//!
	//!    All motion capture position and orientation readings are returned
	//!    in world coordinates.  To convert between reference
	//!    frames use the kinematic transform functions.
	//!
	//!    Unless otherwise specified, these functions should be called
	//!    only when the system is at least in STATE_INTIALIZED.
	//!
	//! @{
	//------------------------------------------------------------------

	//!
	//!  \brief    Get the mocap system joint angle for a specific joint.
	//!
	//!  \_Parameters
	//!
	//!    \_in   leg    - leg index
	//!    \_in   joint  - joint index
	//!    \_out  angle  - angle (in radians) of joint
	//!
	//!  \_Returns        LD_ERROR_CODE result
	//!
	LD_ERROR_CODE getMocapJointReading( LegIndex leg,
		JointIndex joint,
		float * angle) const;

	//!
	//!  \brief    Get the mocap system body for a specific body.
	//!
	//!  \_Parameters
	//!
	//!    \_in   body    - body index 
	//!    \_out  pos     - 3 dimensional vector describing position
	//!    \_out  orient  - unit quaternion (vector, scalar) describing orientation
	//!    \_out  age     - frames since this data was updated (in frames)
	//!
	//!  \_Returns        LD_ERROR_CODE result
	//!
	LD_ERROR_CODE getMocapBodyReading( BodyIndex body,
		bduVec3f * pos,
		bduVec4f * orient,   
		int * age ) const;

	//!
	//!  \brief    Get the mocap system's marker count.
	//!
	//!  \_Parameters
	//!
	//!    \_out  count   - total number of markers reported by the motion capture system
	//!
	//!  \_Returns        LD_ERROR_CODE result
	//!
	LD_ERROR_CODE getMocapMarkerCount( int * count ) const;
 
	//!
	//!  \brief    Get the mocap system marker position for a specific marker.
	//!
	//!  \_Parameters
	//!
	//!    \_in   inx     - marker index
	//!    \_out  v       - 3 dimensional vector describing marker position
	//!    \_out  age     - frames since this data was updated (in frames)
	//!
	//!  \_Returns        LD_ERROR_CODE result
	//!
	LD_ERROR_CODE getMocapMarkerPosition( int inx,
		bduVec3f * v,
		int * age ) const;

	//!
	//!  \brief    Get complete (single frame) of mocap system robot information.
	//!
	//!  \_Parameters
	//!
	//!    \_out  info    - pointer to empty MocapReadInfo structure, to be filled in
	//!
	//!  \_Returns        LD_ERROR_CODE result
	//!
	LD_ERROR_CODE getBulkMocapReadings( MocapReadInfo *info ) const;

	//!
	//!  \brief    Get the latest terrain reading from one of the terrain boards in use 
	//!
	//!  \_Parameters
	//!
	//!    \_in   index   - index of terrain to read ( 0 - MAX_TERRAINS-1 )
	//!    \_out  id      - terrain identifier ( TERRAIN_NONE - TERRAIN_ZZ )
	//!    \_out  pos     - position of terrain board in mocap system coords
	//!    \_out  orient  - unit quaternion orientation of terrain board in mocap system coords
	//!    \_out  age     - frames since this data was updated 
	//!
	//!  \_Returns        LD_ERROR_CODE result
	//!
	LD_ERROR_CODE getTerrainReadings( int index,
		TerrainID * id,
		bduVec3f * pos, 
		bduVec4f * orient, 
		int * age ) const;

	//!
	//!  \brief    Get complete (single frame) of mocap system terrain information.
	//!
	//!  \_Parameters
	//!
	//!    \_out  info    - pointer to empty TerrainReadInfo structure, to be filled in
	//!
	//!  \_Returns        LD_ERROR_CODE result
	//!
	LD_ERROR_CODE getBulkTerrainReadings( TerrainReadInfo * info ) const;


 /////////////////////////////////////////////////////////////////////////
 //! @}
 //
 //  LittleDog Protected Interface
 //
 //    The following functions are callable from a LittleDog object
 //    or an object derived from LittleDog.
 //
 //
 /////////////////////////////////////////////////////////////////////////

 protected:

	//------------------------------------------------------------------
	//!
	//! @name    User-Definable Virtual Functions
	//!
	//!    Virtual function 'callbacks' that provide hooks into
	//!    the control system.  These functions will be called by
	//!    LittleDog at the appropriate times when it is in
	//!    STATE_TRIAL_PLANNING and STATE_TRIAL_RUNNING.
	//!
	//! @{
	//------------------------------------------------------------------

	//!
	//!  \brief    Virtual function 'callback' for subclass control initialization.
	//!
	//!  \_Description
	//!
	//!    This virtual function is called at the beginning of a
	//!    trial when the system is moving from STATE_ROBOT_CALIBRATED
	//!    to STATE_TRIAL_PLANNING.
	//!
	//!    This function should not be called directly by the user.
	//!
	//!    The default implementation does nothing; any necessary
	//!    initialization of the subclass control system should be
	//!    done during the subclass override of this call.
	//!
	//!    If initialization of the subclass control system fails
	//!    for some reason (failure to find a file containing control
	//!    parameters, for example), the function should return false.
	//!    If false is returned, the system will abort the trial and
	//!    revert to STATE_ROBOT_INITIALIZED.
	//!   
	//!    \e Note: Users should use this function to initialize their
	//!    internal state.  They \e MUST not attempt to actively control
	//!    the robot or adjust robot operational paramaters from this function.
	//!
	//!  \_Returns         true if initialization succeeded, false if not
	//!
	virtual bool initControl( void ) {return true;}

	//!
	//!  \brief    Virtual function 'callback' for subclass control update.
	//!
	//!  \_Description
	//!
	//!    This is a pure virtual function, and \e must be overriden
	//!    by every subclass.
	//!
	//!    This virtual function is called periodically during the
	//!    trial while the system is in both STATE_TRIAL_PLANNING and STATE_TRIAL_RUNNING.
	//!
	//!    This function should not be called directly by the user.
	//!
	virtual void updateControl( void ) = 0;

	//!
	//!  \brief    Virtual function 'callback' for subclass control shutdown.
	//!
	//!  \_Description
	//!
	//!    This virtual function is called at the end of every trial.  
	//!    Usually this happens when the system is moving from
	//!    STATE_TRIAL_RUNNING to STATE_ROBOT_CALIBRATED.  But it is called
	//!    also when a trial is aborted or a connection is lost during a trial.
	//!
	//!    This function should not be called directly by the user.
	//!
	//!    The default implementation does nothing; any necessary
	//!    shutdown of the subclass control system should be
	//!    done during the subclass override of this call.
	//!
	virtual void uninitControl( void ) {;}

	//------------------------------------------------------------------
	//! @}
	//! @name    Set Functions
	//!
	//!    Unless otherwise specified, these functions should be called
	//!    only when the system is in STATE_TRIAL_PLANNING or STATE_TRIAL_RUNNING.
	//!
	//! @{
	//------------------------------------------------------------------

	//!  \brief    Set the servo type on a per joint basis
	//!
	//!            This function is only valid to use when the legs
	//!            are configured in LEG_MODE_INDEPENDENT
	//!
	//!  \_Parameters
	//!
	//!    \_in   leg    - leg index
	//!    \_in   joint  - joint index
	//!    \_in   type   - joint servo type
	//!
	//!    \e type is one of:
	//!      - SERVO_LIMP
	//!      - SERVO_PD
	//!      - SERVO_TORQUE
	//!
	//!  \_Returns         LD_ERROR_CODE result
	//!
	//!      - LD_OKAY
	//!      - LD_BAD_LEG_INDEX
	//!      - LD_BAD_JOINT_INDEX
	//!      - LD_BAD_LEG_SERVO_TYPE
	//!
	LD_ERROR_CODE setJointServoType( LegIndex leg,
		JointIndex joint,
		JointServoType type );

	//!
	//!  \brief    Specify the joint level PD servo gains for a joint in PD Servo mode
	//!
	//!  \_Description
	//!
	//!    \e Note:  This function only allows servo gain values to be set if the 
	//!              joint is already in SERVO_PD servo mode.
	//!
	//!  \_Parameters
	//!
	//!    \_in   leg       - leg index
	//!    \_in   joint     - joint index
	//!    \_in   k         - proportional gain
	//!    \_in   b         - derivative gain
	//!
	//!  \_Returns         LD_ERROR_CODE result
	//!
	//!      - LD_OKAY
	//!      - LD_BAD_LEG_INDEX
	//!      - LD_BAD_JOINT_INDEX
	//!      - LD_BAD_LEG_SERVO_TYPE
	//!      - LD_BAD_PROPORTIONAL_GAIN
	//!      - LD_BAD_DERIVATIVE_GAIN
	//!      - LD_BAD_JOINT_SERVO_TYPE
	//!
	LD_ERROR_CODE setJointPDServoGains( LegIndex leg,
		JointIndex joint,
		float k,
		float b );

	//!
	//!  \brief    Specify the leg level proportional force control gains for each joint
	//!
	//!  \_Description
	//!
	//!    This force control is an optional extension of SERVO_PD control.  The measured
	//!    foot force errors are multiplied by the corresponding gains on a per leg and joint
	//!    basis, resulting in joint-specific proportional force control.  To opt-out of 
	//!    using the force control (opting-out is the default), leave the force control 
	//!    gains set to zero.  
	//!
	//!    gains, a 3x3 matrix, has one row per joint (HIP_RX, HIP_RY, KNEE_RY) and one column
	//!    per force component (x, y, z).  Calling this method with no second argument
	//!    passed in sets all of the proportional force gains to zero.
	//!
	//!    \e Note:  This function only allows servo gain values to be set if the entire leg
	//!              has all of it's joints is already in SERVO_PD servo mode.
	//!
	//!  \_Parameters
	//!
	//!    \_in   leg       - leg index
	//!    \_in   gains     - proportional gains
	//!
	//!  \_Returns         LD_ERROR_CODE result
	//!
	//!      - LD_OKAY
	//!      - LD_BAD_LEG_INDEX
	//!      - LD_INVALID_FUNCTION_CALL
	//!      - LD_BAD_JOINT_SERVO_TYPE
	//!
	LD_ERROR_CODE setLegForceControlGains( LegIndex leg,
		const bduMat3f & gains = bduMat3f(bduMat3f::INITIALIZE_TO_ZERO) );

	//!
	//!  \brief    Specify the joint level proportional force control gains for a joint.
	//!
	//!  \_Description
	//!
	//!    \e Note:  This function is a joint-specific version of setLegForceControlGains().
	//!
	//!  \_Parameters
	//!
	//!    \_in   joint     - joint index
	//!    \_in   gains     - proportional gains
	//!
	//!  \_Returns         LD_ERROR_CODE result
	//!
	//!      - LD_OKAY
	//!      - LD_BAD_JOINT_INDEX
	//!      - LD_BAD_LEG_INDEX
	//!      - LD_INVALID_FUNCTION_CALL
	//!      - LD_BAD_JOINT_SERVO_TYPE
	//!
	LD_ERROR_CODE setJointForceControlGains( LegIndex leg, 
		JointIndex joint,
		const bduVec3f & gains = bduVec3f(0) );

	//!
	//!  \brief    Specify the joint level PD servo set points for a joint in PD Servo mode.
	//!
	//!  \_Description
	//!
	//!    \e Note:  This function only allows set point values to be set if the 
	//!              joint is already in SERVO_PD servo mode.
	//!
	//!  \_Parameters
	//!
	//!    \_in   leg       - leg index
	//!    \_in   joint     - joint index
	//!    \_in   q_d       - desired position
	//!    \_in   qd_d      - desired velocity
	//!    \_in   ff_d      - desired feedforward torque
	//!
	//!  \_Returns         LD_ERROR_CODE result
	//!
	LD_ERROR_CODE setJointPDServoSetPoints( LegIndex leg,
		JointIndex joint,
		float q_d,
		float qd_d = 0.0f,
		float ff_d = 0.0f );

	//!
	//!  \brief    Specify the desired joint level output torque for a joint in TAU Servo mode.
	//!
	//!  \_Description
	//!
	//!    \e Note:  This function only allows a torque to be set if the 
	//!              joint is already in SERVO_TORQUE servo mode
	//!
	//!  \_Parameters
	//!
	//!    \_in   leg       - leg index
	//!    \_in   joint     - joint index
	//!    \_in   tau_d     - desired torque in Newton meters
	//!
	//!  \_Returns         LD_ERROR_CODE result
	//!
	LD_ERROR_CODE setJointTorqueServoSetPoints( LegIndex leg,
		JointIndex joint,
		float tau_d );

	//!   
	//!  \brief    Specify the desired foot force for a leg in SERVO_PD mode.
	//!
	//!  \_Description
	//!
	//!    The f_d vector contains desired force components in Newtons for a particular leg.
	//!    The vector has (x, y, z) order.  See setLegForceControlGains().
	//!
	//!  \_Parameters
	//!
	//!    \_in   leg       - leg index
	//!    \_in   f_d       - desired force vector in Newtons (x,y,z)
	//!
	//!  \_Returns         LD_ERROR_CODE result
	//!
	//!      - LD_OKAY
	//!      - LD_BAD_LEG_INDEX
	//!      - LD_INVALID_FUNCTION_CALL
	//!
	LD_ERROR_CODE setLegForceControlSetPoints( LegIndex leg,
		const bduVec3f & f_d = bduVec3f(0) );

	//!
	//!  \brief    Specify foot position to servo to via a 3d location.
	//!
	//!  \_Description
	//!
	//!    For a leg in LEG_MODE_CARTESIAN servo command mode, this function
	//!    specifies the foot position in body coordinates that should
	//!    be reached.  If the leg is in another command mode, this function
	//!    will fail.
	//!
	//!    Note that input "bend_forward" is shorthand for KNEE_BENT_FORWARD.  Typically
	//!    for LittleDog the front-leg knees are bent rearward while the hind-leg
	//!    knees are bent forward.  To clarify, assume that the leg in question
	//!    has its foot directly beneath the hip.  There can be two kinematic solutions 
	//!    to this when the hip_ry is not zero-angled that depend upon the sign of the
	//!    hip_ry.  In this case, a positive hip_ry corresponds to knee rearward, 
	//!    and a negative hip_ry corresponds to knee forward.
	//!
	//!  \_Parameters
	//!
	//!    \_in   leg                    - leg index
	//!    \_in   bend_forward           - true or false.  false signifies bend_rearward
	//!    \_in   pos                    - 3d position of foot in robot coordinates
	//!    \_in   vel                    - 3d velocity of foot in robot coordinates
	//!    \_in   ignore_ik_unreachables - true of false
	//!    \_in   det_singularity        - used to test if matrix's determinant
	//!
	//!    The det_singularity argument is a number used to test if a
	//!    matrix's determinant is too small, thus singular.  If this
	//!    turns out to be the case LD_SINGULAR_MATRIX will be returned.
	//!
	//!    When the ignore_ik_unreachables parameter is false and the
	//!    requested position is unreachable, this function does
	//!    not servo the leg and returns an error.
	//!
	//!    When this parameter is true and the requested position is 
	//!    unreachable, the servo will be done to the best  possible
	//!    reachable position (which turns out to be the value that
	//!    computeFootIK() returns whether computeFootIK() succeeds or
	//!    fails).  The velocity input is ignored when IK is unreachable
	//!    and bduVec3f(0,0,0) is used instead.
	//!
	//!    In either case, when the requested position is not reachable,
	//!    the error code LD_CANNOT_REACH_IK_KNEE_BENT_FORWARD or
	//!    LD_CANNOT_REACH_IK_KNEE_BENT_REARWARD is returned.
	//!
	//!  \_Returns         LD_ERROR_CODE result
	//!
	//!      - LD_OKAY
	//!      - LD_BAD_LEG_INDEX
	//!      - LD_CANNOT_REACH_IK_KNEE_BENT_FORWARD
	//!      - LD_CANNOT_REACH_IK_KNEE_BENT_REARWARD
	//!      - LD_SINGULAR_MATRIX
	//!      - LD_BAD_LEG_SERVO_TYPE
	//!      - Errors produced via computeFootIK()
	//!
	LD_ERROR_CODE servoFootPosition ( LegIndex leg,
		bool bend_forward, 
		const bduVec3f & pos,
		const bduVec3f & vel, 
		bool ignore_ik_unreachables = true,
		float det_singularity = 1.0e-8 );

	//!  \brief    Set the leg command mode.
	//!
	//!  \_Description
	//!
	//!    This function sets the current servo command mode for each leg.
	//!
	//!    This function only allows a leg to be set to CARTESIAN if all
	//!    the joints of the leg are currently in PD mode.
	//!
	//!  \_Parameters
	//!
	//!    \_in   leg          - leg index
	//!    \_in   command_mode - LEG_MODE_INDEPENDENT or LEG_MODE_CARTESIAN
	//!
	//!  \_Returns          LD_ERROR_CODE result
	//!
	LD_ERROR_CODE setLegCommandMode( LegIndex leg,
		LegServoCommandMode command_mode);

	//!
	//!  \brief    Specify the joint output torque limits.
	//!
	//!  \_Description
	//!
	//!    \e Note:  This function sets torque limits that are in effect when
	//!              the servos are enabled.  The default limits are USER_DEFAULT_MIN_TAU
	//!              to USER_DEFAULT_MAX_TAU.  If the given min parameter is greater than
	//!              the given max parameter, the parameters are switched.  If no user limits 
	//!              are required, users can call this function with no limits, accepting
	//!              the defaults (see prototype).
	//!
	//!              These limits are reset to the default values when:
	//!
	//!                1) The LittleDog class is instanciating
	//!                2) The JointServoType is switched to SERVO_LIMP mode
	//!                3) Calibration is performed
	//!
	//!  \_Parameters
	//!
	//!    \_in   leg       - leg index
	//!    \_in   joint     - joint index
	//!    \_in   min_tau   - desired minimum torque limit in Newton meters
	//!    \_in   max_tau   - desired maximum torque limit in Newton meters
	//!
	//!  \_Returns         LD_ERROR_CODE result
	//!
	LD_ERROR_CODE setJointTorqueUserLimits( LegIndex leg,
		JointIndex joint,
		float min_tau = USER_DEFAULT_MIN_TAU, 
		float max_tau = USER_DEFAULT_MAX_TAU );


	//------------------------------------------------------------------
	//! @}


 /////////////////////////////////////////////////////////////////////////
 //
 //  Private Interface
 //
 //    The following functions are not callable outside the base
 //    LittleDog object class.
 //
 /////////////////////////////////////////////////////////////////////////

 private:

	//
	//  The copy consructor and operator= are declared private
	//  so that compilers will not generate default versions.
	//
	LittleDog( const LittleDog & );
	LittleDog & operator = ( const LittleDog & );
	
	//
	//  Pointer to internal data.
	//
	class LittleDogPrivate * p;

	//
	//  Internal utility functions.
	//
	LD_ERROR_CODE checkStandardErrors(void) const;
	LD_ERROR_CODE checkConnectedErrors(void) const;

	//
	//  Friend types.
	//
	friend class LittleDogLegKin;
	friend class LittleDogCalibrator;
	friend class LittleDogPrivate;
};

//!  @}

#endif  // __LITTLEDOG_API_H__
