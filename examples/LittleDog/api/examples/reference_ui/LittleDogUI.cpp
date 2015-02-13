
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <qapplication.h>
#include <qbrush.h>
#include <qcheckbox.h>
#include <qcursor.h>
#include <qgroupbox.h>
#include <qlabel.h>
#include <qlayout.h>
#include <qlineedit.h>
#include <qlcdnumber.h>
#include <qmessagebox.h>
#include <qpainter.h>
#include <qpen.h>
#include <qpushbutton.h>
#include <qthread.h>
#include <qtextedit.h>
#include <qtimer.h>

#include "LittleDogUI.h"
#include "DialogCalibratePose_uic.h"

static const int TIMER_UPDATE_LOG_PERIOD = 100;      // milliseconds
static const int TIMER_UPDATE_READINGS_PERIOD = 500; // milliseconds
static const int TIMER_UPDATE_CANVAS_PERIOD = 100;   // milliseconds
static const int MAX_THREAD_JOIN_WAIT_TIME = 10000;  // milliseconds

static const float MIN_ACCEPTABLE_VOLTAGE = 13.0f;
static const float MAX_ACCEPTABLE_TEMP = 85.0f;
static const float MAX_ACCEPTABLE_PACKET_LOSS = 5.0f;


/****************************************************************************/
/*
 *  Since the LittleDog::calibrate() function blocks we need to
 *   run it in a separate thread so that the GUI will remain
 *   responsive.
 */
class LittleDogCalibrateThread : public QThread
{
public:

	LittleDogCalibrateThread(LittleDog& dog, LittleDogUI& ui) : m_dog(dog), m_ui(ui) {}
	LittleDog&           m_dog;
	LittleDogUI&         m_ui;
	static LD_ERROR_CODE m_return_value;

	/*
	 *  Execution of thread ends when this function returns.
	 */
	virtual void run() 
	{
		// passing TRUE forces calibration
		m_return_value = m_dog.calibrate(TRUE);

		if ( m_return_value != LD_OKAY )
			bdu_log_printf(BDU_LOG_LEVEL_ERROR, "CALIBRATE ERROR:  %s\n",
				m_dog.getErrorCodeString(m_return_value));
	}
};

LD_ERROR_CODE LittleDogCalibrateThread::m_return_value = LD_OKAY;


/****************************************************************************/
/*
 *  Since the LittleDog::runTrial() function blocks we need to
 *   run it in a separate thread so that the GUI will remain
 *   responsive.
 */
class LittleDogTrialThread : public QThread
{
public:

	LittleDogTrialThread(LittleDog& dog, LittleDogUI& ui) : m_dog(dog), m_ui(ui) {}
	LittleDog&           m_dog;
	LittleDogUI&         m_ui;
	static LD_ERROR_CODE m_return_value;

	/*
	 *  Execution of thread ends when this function returns.
	 */
	virtual void run()
	{
		// Optionally start a timer.  It is only accurately canceled at 
		// frequency (1 / TIMER_UPDATE_READINGS_PERIOD * 0.001)
		m_ui.start_planning_timer();

		m_return_value = m_dog.runTrial();

		if ( !m_ui.get_goal_reached() )
			m_ui.generate_trial_summary();
	}
};

LD_ERROR_CODE LittleDogTrialThread::m_return_value = LD_OKAY;

/****************************************************************************/
/*
 *  Helper class that shows the hourglass when it is created, and
 *   restores the default cursor when it is destroyed.  Helps prevent
 *   missed restores by using local scope variables.
 */
class Hourglass
{
public:

	Hourglass()  {QApplication::setOverrideCursor(QCursor(Qt::WaitCursor));}
	~Hourglass() {QApplication::restoreOverrideCursor();}
};


/****************************************************************************/
/*
 *  The bduLog class can take a callback that has this exact signature.
 *   Forward the print to the LittleDogUI class's queue_to_log_widget()
 *   function.
 */
int
bdu_log_callback_function(int notify_level, const char* string, void* user_data)
{
	LittleDogUI* ui = (LittleDogUI*) user_data;

	if (ui)
		ui->queue_to_log_widget(notify_level, string);

	return 0;
}


/****************************************************************************/
/*
 *  LittleDogUI constructor.
 */
LittleDogUI::LittleDogUI(LittleDog& dog, int max_planning_duration_sec) :
	m_dog(dog),
	m_trial_thread(NULL),
	m_calibrate_thread(NULL),
	m_timer_planning_exceeded(NULL),
	m_last_elapsed_trial_run_duration(0),
	m_last_elapsed_trial_plan_duration(0),
	m_last_elapsed_idle_duration(0),
	m_max_planning_duration_sec(max_planning_duration_sec),
	m_flash_on(FALSE),
	m_temp_went_too_high(FALSE)
{
	/*
	 *  Set up a log that prints to an on-screen widget.
	 */
	m_log_callback = new bduLog(BDU_LOG_LEVEL_INFO, // WARN,
		bdu_log_callback_function,
		this);

	/*
	 *  Set up a timer that periodically updates the log widget.
	 */
	m_timer_update_log = new QTimer(this);
	connect(m_timer_update_log, SIGNAL(timeout()), this, SLOT(slot_update_log()));
	m_timer_update_log->start(TIMER_UPDATE_LOG_PERIOD);

	/*
	 *  Set up a timer that periodically updates readings, buttons,
	 *   and state labels.
	 */
	m_timer_update_readings = new QTimer(this);
	connect(m_timer_update_readings, SIGNAL(timeout()), this, SLOT(slot_update_readings()));
	connect(m_timer_update_readings, SIGNAL(timeout()), this, SLOT(slot_update_buttons_and_state_labels()));
	m_timer_update_readings->start(TIMER_UPDATE_READINGS_PERIOD);

	if ( m_max_planning_duration_sec != -1 )
	{
		m_timer_planning_exceeded = new QTimer(this);
		connect(m_timer_planning_exceeded, SIGNAL(timeout()), this, SLOT(slot_planning_exceeded()));
	}

	/*
	 *  Create an overhead view in m_frame_overhead_view.
	 */
	m_overhead_view = new OverheadViewWidget(m_frame_overhead_view, "overhead_view");

	QVBoxLayout* overhead_view_layout = new QVBoxLayout(m_frame_overhead_view, 4, 4, "overhead_view_layout");
	overhead_view_layout->addWidget(m_overhead_view);
	overhead_view_layout->activate();

	m_timer_update_overhead_view = new QTimer(this);
	connect(m_timer_update_overhead_view, SIGNAL(timeout()), this, SLOT(slot_update_overhead_view()));
	m_timer_update_overhead_view->start(TIMER_UPDATE_CANVAS_PERIOD);
	m_last_state = m_dog.getState();

	/*
	 *  Initial update of buttons and labels.
	 */
	slot_update_buttons_and_state_labels();
}


/****************************************************************************/
/*
 *  LittleDogUI destructor.
 */
LittleDogUI::~LittleDogUI()
{
	if (m_calibrate_thread)
	{
		wait_for_calibrate_thread(TRUE);
		delete m_calibrate_thread;
		m_calibrate_thread = NULL;
	}
	if (m_trial_thread)
	{
		wait_for_trial_thread(TRUE);
		delete m_trial_thread;
		m_trial_thread = NULL;
	}

	m_timer_update_log->stop();
	m_timer_update_readings->stop();
	m_timer_update_overhead_view->stop();

	stop_planning_timer();

	// Programming note: it is always ok to call delete on NULL pointers
	delete m_timer_update_readings;
	delete m_timer_planning_exceeded;
	delete m_log_callback;
}


/****************************************************************************/
bool
LittleDogUI::get_goal_reached() const
{
	if ( m_check_goal_reached->isChecked() )
		return true;

	if ( m_overhead_view->get_goal() == false )
		return false;

	float distance = 100.0 * m_overhead_view->get_distance_to_goal(m_dog);  // m -> cm

	return (distance <= 5.0) ? true : false;
}


/****************************************************************************/
/*
 *  This function optionally starts a single-shot timer that is used to 
 *   estimate how long the runTrial() call has been in STATE_TRIAL_PLANNING.
 */	
void
LittleDogUI::start_planning_timer()
{
	if ( m_timer_planning_exceeded )
		m_timer_planning_exceeded->start( 
			m_max_planning_duration_sec * 1000, TRUE ); // sec->msec
}


/****************************************************************************/
/*
 *  This function stops the planning timer, if it has been instantiated.
 *   It is ok to stop a timer that is already stopped.
 */
void
LittleDogUI::stop_planning_timer()
{
	if ( m_timer_planning_exceeded )
		m_timer_planning_exceeded->stop();
}


/****************************************************************************/
void
LittleDogUI::generate_trial_summary()
{
	m_last_elapsed_idle_duration = m_time_elapsed_idle.elapsed();
	m_time_elapsed_idle.start();

	/*
	 *  To facilitate testing, print out some information.
	 */
	
	float end_x, end_y, start_x, start_y, goal_x, goal_y;
	
	if ( (m_overhead_view->get_position(m_dog, &end_x, &end_y) == true) &&
		(m_overhead_view->get_start_position(&start_x, &start_y) == true) &
		(m_overhead_view->get_goal(&goal_x, &goal_y) == true ) )
	{
		float duration               = 0.001f * m_time_elapsed_trial_run.elapsed();
		float end_distance_to_goal   = sqrt( ((end_x - goal_x) * (end_x - goal_x)) + 
			((end_y - goal_y) * (end_y - goal_y)) );
		float start_distance_to_goal = sqrt( ((start_x - goal_x) * (start_x - goal_x)) + 
			((start_y - goal_y) * (start_y - goal_y)) );
		bool goal_was_reached        =  get_goal_reached();
		
		bdu_log_printf( BDU_LOG_LEVEL_WARN, "*** \n");
		bdu_log_printf( BDU_LOG_LEVEL_WARN, "*** Goal location     (x,y) (m): (%f,%f)\n",
			goal_x, goal_y );
		bdu_log_printf( BDU_LOG_LEVEL_WARN, "*** Started at        (x,y) (m): (%f,%f)\n", 
			start_x, start_y );
		bdu_log_printf( BDU_LOG_LEVEL_WARN, "*** Ended at          (x,y) (m): (%f,%f)\n", 
			end_x, end_y );
		bdu_log_printf( BDU_LOG_LEVEL_WARN, "*** Start dist. to goal     (m): %f\n", 
			start_distance_to_goal );
		bdu_log_printf( BDU_LOG_LEVEL_WARN, "*** End dist. to goal       (m): %f\n",
			end_distance_to_goal );
		bdu_log_printf( BDU_LOG_LEVEL_WARN, "*** Trial duration       (sec.): %f\n",
			duration );
		bdu_log_printf( BDU_LOG_LEVEL_WARN, "*** Goal reached         (bool): %s\n",
			goal_was_reached ? "Y" : "N" );
		
		bdu_log_printf( BDU_LOG_LEVEL_WARN, "*** \n");
	}
}


/****************************************************************************/
/*
 *  This function doesn't immediately print to m_edit_log
 *   as it may be called from a number of threads.  Instead
 *   it appends the passed string to a string that will
 *   be appended to m_edit_log's contents when
 *   m_timer_update_log times out.
 */
void
LittleDogUI::queue_to_log_widget(int /*level*/,
	const char* string)
{
	/*
	 *  - Temporarily suspend log update timer.
	 *  - Get lock for m_log_append_string.
	 *  - Append the string to m_log_append_string.
	 *  - Release lock for m_log_append_string.
	 *  - Resume log update timer.
	 */
	m_timer_update_log->stop();
	m_log_append_string_mutex.lock();

	m_log_append_string += string;

	m_log_append_string_mutex.unlock();
	m_timer_update_log->start(TIMER_UPDATE_LOG_PERIOD);
}


/****************************************************************************/
void
LittleDogUI::slot_clicked_connect_to_mocap_system()
{
	Hourglass hg;

	LD_ERROR_CODE err = m_dog.initializeMocap();
	if (err != LD_OKAY)
	{
		bdu_log_printf(BDU_LOG_LEVEL_ERROR, "ERROR:  %s\n",
			m_dog.getErrorCodeString(err));
		return;
	}

	bdu_log_print(BDU_LOG_LEVEL_WARN, "Successfully connected to mocap system.\n");
	m_overhead_view->get_goal();

	slot_update_buttons_and_state_labels();
}


/****************************************************************************/
void
LittleDogUI::slot_clicked_disconnect_from_mocap_system()
{
	Hourglass hg;

	LD_ERROR_CODE err = m_dog.stopMocap();
	if (err != LD_OKAY)
	{
		bdu_log_printf(BDU_LOG_LEVEL_ERROR, "ERROR:  %s\n",
			m_dog.getErrorCodeString(err));
		return;
	}

	bdu_log_print(BDU_LOG_LEVEL_WARN, "Successfully disconnected from mocap system.\n");

	slot_update_buttons_and_state_labels();
}


/****************************************************************************/
void
LittleDogUI::slot_clicked_connect_to_robot()
{
	Hourglass hg;

	LD_ERROR_CODE err = m_dog.initializeRobot();
	if (err != LD_OKAY)
	{
		bdu_log_printf(BDU_LOG_LEVEL_ERROR, "ERROR:  %s\n",
			m_dog.getErrorCodeString(err));
		return;
	}

	bdu_log_print(BDU_LOG_LEVEL_WARN, "Successfully connected to robot.\n");

	m_last_elapsed_idle_duration = m_time_elapsed_idle.elapsed();
	m_time_elapsed_idle.start();
	slot_update_buttons_and_state_labels();
}


/****************************************************************************/
void
LittleDogUI::slot_clicked_calibrate_robot()
{
	/*
	 *  Show a calibration dialog.
	 */
	DialogCalibratePose_uic dlg;

	QString pixmap_filename = "";
	const char* littledog_env_var = getenv("LITTLEDOG");
	if (littledog_env_var)
	{
		pixmap_filename = littledog_env_var;
		pixmap_filename += "/images/littledog_calibration_pose.jpg";

		QPixmap pixmap(pixmap_filename);
		if (!pixmap.isNull())
			dlg.m_label_pixmap->setPixmap(pixmap);
	}

	if (dlg.exec() != QDialog::Accepted)
		return;

	/*
	 *  May take awhile.  Show an hourglass.
	 */
	Hourglass hg;

	start_calibrate_thread();
}


/****************************************************************************/
void
LittleDogUI::slot_clicked_run_trial()
{
	bdu_log_printf( BDU_LOG_LEVEL_WARN, "*** Begin Trial Planning requested!\n" );

	if ( m_max_planning_duration_sec != -1 )
		bdu_log_printf( BDU_LOG_LEVEL_WARN, "*** Max Planning time is %d sec.\n",
			m_max_planning_duration_sec );
	else
		bdu_log_printf( BDU_LOG_LEVEL_WARN, "*** No Planning time limit requested.\n");

	/*
	 *  To facitate reporting, set the start_postion if we can.
	 */

	float x, y;

	if ( m_overhead_view->get_position(m_dog, &x, &y) == true )
		m_overhead_view->set_start_position( true, x, y );
	else
		m_overhead_view->set_start_position( false, 0, 0 );
	
	/*
	 *  Start the trial.
	 */

	start_trial_thread();
}


/****************************************************************************/
void
LittleDogUI::slot_clicked_done_planning()
{
	bdu_log_printf( BDU_LOG_LEVEL_WARN, "*** Begin Trial Running ... (planning done)\n" );

	m_dog.donePlanning();
}


/****************************************************************************/
void
LittleDogUI::slot_clicked_stop_robot()
{
	bdu_log_printf( BDU_LOG_LEVEL_WARN, "*** Stop Trial requested!\n" );

	/*
	 *  Stop threads.
	 */
	if (m_calibrate_thread)
		stop_calibrate_thread();
	if (m_trial_thread)
		stop_trial_thread();
}


/****************************************************************************/
void
LittleDogUI::slot_clicked_quit()
{
	/*
	 *  Tell the system to abort connections.  This should
	 *   stop everything and put the system into
	 *   STATE_ABORTING_CONNECTION.
	 */
	m_dog.abortConnection();

	/*
	 *  Immediately update buttons, etc.
	 */
	slot_update_buttons_and_state_labels();
	slot_update_readings();

	/*
	 *  Wait for threads to join.
	 */
	if (m_calibrate_thread)
		wait_for_calibrate_thread(TRUE);
	if (m_trial_thread)
		wait_for_trial_thread(TRUE);

	/*
	 *  Let closeEvent() handle shutting down.
	 */
	close();
}


/****************************************************************************/
void
LittleDogUI::slot_clicked_abort_connection()
{
	bdu_log_print(BDU_LOG_LEVEL_WARN, "*** Abort requested!\n");

	/*
	 *  Tell the system to abort connections.  This should
	 *   stop everything and put the system into
	 *   STATE_ABORTING_CONNECTION.
	 */
	m_dog.abortConnection();

	/*
	 *  Immediately update buttons, etc.
	 */
	slot_update_buttons_and_state_labels();
	slot_update_readings();

	/*
	 *  Wait for threads to join.
	 */
	if (m_calibrate_thread)
		wait_for_calibrate_thread(TRUE);
	if (m_trial_thread)
		wait_for_trial_thread(TRUE);
}


/****************************************************************************/
void
LittleDogUI::slot_pressed_button()
{
	/*
	 *  Qt "loses" the release event of button presses if
	 *   the timer goes off and the readings labels are
	 *   updated, so turn off the timer temporarily until
	 *   the button is released.
	 */
	if (m_timer_update_readings->isActive())
		m_timer_update_readings->stop();

	if (m_timer_update_overhead_view->isActive())
		m_timer_update_overhead_view->stop();
}


/****************************************************************************/
void
LittleDogUI::slot_released_button()
{
	/*
	 *  Restart update timer.
	 */
	if (!m_timer_update_readings->isActive())
		m_timer_update_readings->start(TIMER_UPDATE_READINGS_PERIOD);

	if (!m_timer_update_overhead_view->isActive())
		m_timer_update_overhead_view->start(TIMER_UPDATE_CANVAS_PERIOD);
}


/****************************************************************************/
/*
 *  This function is called when close() has been called on this widget.
 *   This can be due to a direct call (as in the slot_clicked_quit()
 *   function) or indirectly (user pressed a Close button on a
 *   titlebar).  Handle either case the same way.
 */
void
LittleDogUI::closeEvent(QCloseEvent* e)
{
	bdu_log_print(BDU_LOG_LEVEL_WARN, "Closing LittleDogUI...\n");

	if (m_timer_update_readings->isActive())
		m_timer_update_readings->stop();

	/*
	 *  Check here to see if a trial is running, and if so, ask
	 *   users if they really want to quit.
	 */
	LittleDog::State dog_state = m_dog.getState();
	if ( (dog_state == LittleDog::STATE_TRIAL_PLANNING) ||
		(dog_state == LittleDog::STATE_TRIAL_RUNNING) )
	{
		switch(QMessageBox::warning(this,
			"LittleDog",
			"Trial is in progress.  Really quit?",
			"&Yes", "&No", QString::null,
			0,      // Enter == Yes (button 0)
			1))     // Escape == No (button 1)
		{
			case 0:
				break;

			default:
				m_timer_update_readings->start(TIMER_UPDATE_READINGS_PERIOD);
				bdu_log_print(BDU_LOG_LEVEL_WARN, "... Quit aborted.\n");
				return;
		}
	}

	/*
	 *  If a trial is running or the robot is calibrating, call
	 *   abortConnection() first.
	 *
	 *  Otherwise we should just shut down.
	 */
	if ( (dog_state == LittleDog::STATE_TRIAL_PLANNING) ||
		(dog_state == LittleDog::STATE_TRIAL_RUNNING) ||
		(dog_state == LittleDog::STATE_ROBOT_CALIBRATING) )
	{
		/*
		 *  Call slot_clicked_abort_connection() instead of
		 *   directly calling m_dog.abortConnection(), as the
		 *   slot call will correctly wait for threads to end.
		 */
		slot_clicked_abort_connection();
	}

	bdu_log_print(BDU_LOG_LEVEL_WARN, "... done.\n");

	LittleDogUI_uic::closeEvent(e);
}


/****************************************************************************/
/*
 *  Helper function for starting m_calibrate_thread.
 */
void
LittleDogUI::start_calibrate_thread()
{
	if (m_calibrate_thread)
	{
		if (m_calibrate_thread->running() ||
			!m_calibrate_thread->finished())
		{
			bdu_log_print(BDU_LOG_LEVEL_WARN,
				"WARNING:  Attempted to start calibrate thread, but one is already running.\n");
			return;
		}
	}

	bdu_log_print(BDU_LOG_LEVEL_WARN, "Starting calibration...\n");

	/*
	 *  Create the calibrate thread.
	 */
	if (!m_calibrate_thread)
		m_calibrate_thread = new LittleDogCalibrateThread(m_dog,*this);

	/*
	 *  Start it.  The start() function will fork the thread
	 *   and immediately return after that.  The thread's run()
	 *   function will be executed, and the thread will terminate
	 *   when the run() function 'returns'.
	 */
	m_calibrate_thread->start();

	m_last_elapsed_idle_duration = m_time_elapsed_idle.elapsed();
	m_last_elapsed_trial_run_duration = 0;
	m_last_elapsed_trial_plan_duration = 0;

	bdu_log_print(BDU_LOG_LEVEL_WARN, "... done.\n");

	slot_update_buttons_and_state_labels();
}


/****************************************************************************/
/*
 *  Helper function for stopping m_calibrate_thread.
 */
void
LittleDogUI::stop_calibrate_thread()
{
	if (!m_calibrate_thread)
		return;

	if (!m_calibrate_thread->running() ||
		m_calibrate_thread->finished())
	{
		return;
	}

	LittleDog::State dog_state = m_dog.getState();
	if (dog_state != LittleDog::STATE_ROBOT_CALIBRATING)
		return;

	bdu_log_print(BDU_LOG_LEVEL_WARN, "Stopping calibration...\n");

	/*
	 *  Call stopRobot() to signal the calibration to end.  This
	 *   call immediately returns.  At some point the calibrate()
	 *   function called in LittleDogCalibrateThread::run()
	 *   should return, ending the thread.
	 */
	m_dog.stopRobot();
	m_last_elapsed_idle_duration = 0;
	m_time_elapsed_idle.start();

	/*
	 *  Wait for the calibrate thread to join.
	 */
	wait_for_calibrate_thread(FALSE);

	bdu_log_print(BDU_LOG_LEVEL_WARN, "... done.\n");

	slot_update_buttons_and_state_labels();
}


/****************************************************************************/
/*
 *  Helper function for joining m_calibrate_thread.
 */
void
LittleDogUI::wait_for_calibrate_thread(bool force_termination)
{
	if (!m_calibrate_thread)
		return;

	if (!m_calibrate_thread->running() ||
		m_calibrate_thread->finished())
	{
		return;
	}

	Hourglass hg;

	/*
	 *  Wait for thread to end.  wait() returns false
	 *   if the wait timed-out.  The wait should
	 *   end when the calibrate() function returns.
	 */
	unsigned long wait_time = force_termination ? MAX_THREAD_JOIN_WAIT_TIME : ULONG_MAX;

	if (!m_calibrate_thread->wait(wait_time))  // wait returns FALSE if it timed out
	{
		bdu_log_print(BDU_LOG_LEVEL_ERROR, "ERROR:  Calibrate thread not ending!\n");

		if (force_termination)
		{
			bdu_log_print(BDU_LOG_LEVEL_ERROR, "Forcing termination of calibrate thread!\n");

			m_calibrate_thread->terminate();
			m_calibrate_thread->wait(MAX_THREAD_JOIN_WAIT_TIME);
		}
	}
	else
	{
		if (m_calibrate_thread->m_return_value != LD_OKAY)
		{
			bdu_log_printf(BDU_LOG_LEVEL_ERROR,
				"ERROR:  Calibrate ended with error '%s'.\n",
				m_dog.getErrorCodeString(m_calibrate_thread->m_return_value));
		}
	}
}


/****************************************************************************/
/*
 *  Helper function for starting m_trial_thread.
 */
void
LittleDogUI::start_trial_thread()
{
	if (m_trial_thread)
	{
		if (m_trial_thread->running() ||
			!m_trial_thread->finished())
		{
			bdu_log_print(BDU_LOG_LEVEL_WARN,
				"WARNING:  Attempted to start trial thread, but one is already running.\n");
			return;
		}
	}

	bdu_log_print(BDU_LOG_LEVEL_WARN, "Starting trial...\n");

	/*
	 *  Create the trial thread.
	 */
	if (!m_trial_thread)
		m_trial_thread = new LittleDogTrialThread(m_dog, *this);

	/*
	 *  Start it.  The start() function will fork the thread
	 *   and immediately return after that.  The thread's run()
	 *   function will be executed, and the thread will terminate
	 *   when the run() function 'returns'.
	 */
	m_trial_thread->start();

	m_last_elapsed_idle_duration = m_time_elapsed_idle.elapsed();
	m_last_elapsed_trial_run_duration = 0;
	m_last_elapsed_trial_plan_duration = 0;

	m_time_elapsed_trial_plan.start();

	m_check_goal_reached->setChecked(FALSE);

	bdu_log_print(BDU_LOG_LEVEL_WARN, "... done.\n");

	slot_update_buttons_and_state_labels();
}


/****************************************************************************/
/*
 *  Helper function for stopping m_trial_thread.
 */
void
LittleDogUI::stop_trial_thread()
{
	if (!m_trial_thread)
		return;

	if (!m_trial_thread->running() ||
		m_trial_thread->finished())
	{
		return;
	}

	LittleDog::State dog_state = m_dog.getState();
	if ( (dog_state != LittleDog::STATE_TRIAL_PLANNING) &&
		(dog_state != LittleDog::STATE_TRIAL_RUNNING) )
		return;

	bdu_log_print(BDU_LOG_LEVEL_WARN, "Stopping trial...\n");

	/*
	 *  Call stopRobot() to signal the trial to end.  This
	 *   call immediately returns.  At some point the runTrial()
	 *   function called in LittleDogTrialThread::run()
	 *   should return, ending the thread.
	 */
	m_dog.stopRobot();

	if ( dog_state == LittleDog::STATE_TRIAL_PLANNING )
		m_last_elapsed_trial_plan_duration = m_time_elapsed_trial_plan.elapsed();
	else
		m_last_elapsed_trial_run_duration = m_time_elapsed_trial_run.elapsed();

	m_last_elapsed_idle_duration = 0;
	m_time_elapsed_idle.start();

	/*
	 *  Wait for the trial thread to join.
	 */
	wait_for_trial_thread(FALSE);

	bdu_log_print(BDU_LOG_LEVEL_WARN, "... done.\n");

	slot_update_buttons_and_state_labels();
}


/****************************************************************************/
/*
 *  Helper function for joining m_trial_thread.
 */
void
LittleDogUI::wait_for_trial_thread(bool force_termination)
{
	if (!m_trial_thread)
		return;

	if (!m_trial_thread->running() ||
		m_trial_thread->finished())
	{
		return;
	}

	Hourglass hg;

	/*
	 *  Wait for thread to end.  wait() returns false
	 *   if the wait timed-out.  The wait should
	 *   end when the runTrial() function returns.
	 */
	unsigned long wait_time = force_termination ? MAX_THREAD_JOIN_WAIT_TIME : ULONG_MAX;

	if (!m_trial_thread->wait(wait_time))  // wait returns FALSE if it timed out
	{
		bdu_log_print(BDU_LOG_LEVEL_ERROR, "ERROR:  Trial thread not ending.\n");

		if (force_termination)
		{
			bdu_log_print(BDU_LOG_LEVEL_ERROR, "Forcing termination of trial thread!\n");

			m_trial_thread->terminate();
			m_trial_thread->wait(MAX_THREAD_JOIN_WAIT_TIME);
		}
	}
	else
	{
		if (m_trial_thread->m_return_value != LD_OKAY)
		{
			bdu_log_printf(BDU_LOG_LEVEL_ERROR,
				"ERROR:  Trial ended with error '%s'.\n",
				m_dog.getErrorCodeString(m_trial_thread->m_return_value));
		}
	}
}


/****************************************************************************/
/*
 *  This function updates buttons and state labels based on the
 *   current state of the m_dog object.  It is called periodically
 *   when m_timer_update_readings times out.
 */
void
LittleDogUI::slot_update_buttons_and_state_labels()
{
	m_check_calibrate_thread_running->setChecked(FALSE);
	if (m_calibrate_thread && m_calibrate_thread->running())
		m_check_calibrate_thread_running->setChecked(TRUE);

	m_check_trial_thread_running->setChecked(FALSE);
	if (m_trial_thread && m_trial_thread->running())
		m_check_trial_thread_running->setChecked(TRUE);

	m_check_mocap_thread_running->setChecked(m_dog.getMocapInitialized());

	m_check_mocap_data->setChecked(FALSE);
	if ( m_dog.getMocapDataReceived() )
		m_check_mocap_data->setChecked(TRUE);

	/*
	 *  Start with a known state.
	 */
	m_button_connect_to_mocap_system->setEnabled(FALSE);
	m_button_disconnect_from_mocap_system->setEnabled(FALSE);
	m_button_connect_to_robot->setEnabled(FALSE);
	m_button_calibrate_robot->setEnabled(FALSE);
	m_button_run_trial->setEnabled(FALSE);
	m_button_done_planning->setEnabled(FALSE);
	m_button_stop_robot->setEnabled(FALSE);
	m_button_abort_connection->setEnabled(FALSE);
	m_button_quit->setEnabled(TRUE);

	bool mocap_initialized = m_dog.getMocapInitialized();

	/*
	 *  Update buttons based on m_dog state.
	 */
	LittleDog::State state = m_dog.getState();

	switch (state)
	{
		case LittleDog::STATE_UNINITIALIZED:
			m_button_connect_to_mocap_system->setEnabled(!mocap_initialized);
			m_button_disconnect_from_mocap_system->setEnabled(mocap_initialized);
			m_button_connect_to_robot->setEnabled(TRUE);
			break;

		case LittleDog::STATE_ROBOT_INITIALIZED:
			m_button_connect_to_mocap_system->setEnabled(!mocap_initialized);
			m_button_disconnect_from_mocap_system->setEnabled(mocap_initialized);
			m_button_calibrate_robot->setEnabled(TRUE);
			m_button_abort_connection->setEnabled(TRUE);
			break;

		case LittleDog::STATE_ROBOT_CALIBRATING:
			m_button_stop_robot->setEnabled(TRUE);
			m_button_abort_connection->setEnabled(TRUE);
			break;

		case LittleDog::STATE_ROBOT_CALIBRATED:
			m_button_connect_to_mocap_system->setEnabled(!mocap_initialized);
			m_button_disconnect_from_mocap_system->setEnabled(mocap_initialized);
			m_button_calibrate_robot->setEnabled(TRUE);
			m_button_run_trial->setEnabled(TRUE);
			m_button_abort_connection->setEnabled(TRUE);
			break;

		case LittleDog::STATE_TRIAL_PLANNING:
			m_button_done_planning->setEnabled(TRUE);
			m_button_abort_connection->setEnabled(TRUE);
			m_button_stop_robot->setEnabled(TRUE);
			break;

		case LittleDog::STATE_TRIAL_RUNNING:
			m_button_abort_connection->setEnabled(TRUE);
			m_button_stop_robot->setEnabled(TRUE);
			break;

		case LittleDog::STATE_STOPPING_ROBOT:
			break;

		case LittleDog::STATE_ABORTING_CONNECTION:
			break;
	}

	/*
	 *  Update labels based on m_dog state.
	 */
	QColor color_fg_disabled(255, 255, 255);
	QColor color_bg_disabled(96, 96, 96);

	QColor color_fg_enabled(255, 255, 255);
	QColor color_bg_enabled(0, 128, 0);

	QColor color_fg_abort(255, 255, 255);
	QColor color_bg_abort(128, 0, 0);

	QColor color_fg_current_on(0, 0, 0);
	QColor color_bg_current_on(255, 255, 255);

	QColor color_fg_current_off(255, 255, 255);
	QColor color_bg_current_off(0, 0, 128);

	QColor color_fg_current;
	QColor color_bg_current;

	m_label_STATE_UNINITIALIZED->setPaletteForegroundColor(color_fg_disabled);
	m_label_STATE_UNINITIALIZED->setPaletteBackgroundColor(color_bg_disabled);

	m_label_STATE_ROBOT_INITIALIZED->setPaletteForegroundColor(color_fg_disabled);
	m_label_STATE_ROBOT_INITIALIZED->setPaletteBackgroundColor(color_bg_disabled);

	m_label_STATE_ROBOT_CALIBRATING->setPaletteForegroundColor(color_fg_disabled);
	m_label_STATE_ROBOT_CALIBRATING->setPaletteBackgroundColor(color_bg_disabled);

	m_label_STATE_ROBOT_CALIBRATED->setPaletteForegroundColor(color_fg_disabled);
	m_label_STATE_ROBOT_CALIBRATED->setPaletteBackgroundColor(color_bg_disabled);

	m_label_STATE_TRIAL_PLANNING->setPaletteForegroundColor(color_fg_disabled);
	m_label_STATE_TRIAL_PLANNING->setPaletteBackgroundColor(color_bg_disabled);

	m_label_STATE_TRIAL_RUNNING->setPaletteForegroundColor(color_fg_disabled);
	m_label_STATE_TRIAL_RUNNING->setPaletteBackgroundColor(color_bg_disabled);

	m_label_STATE_STOPPING_ROBOT->setPaletteForegroundColor(color_fg_disabled);
	m_label_STATE_STOPPING_ROBOT->setPaletteBackgroundColor(color_bg_disabled);

	m_label_STATE_ABORTING_CONNECTION->setPaletteForegroundColor(color_fg_disabled);
	m_label_STATE_ABORTING_CONNECTION->setPaletteBackgroundColor(color_bg_disabled);

	if (m_flash_on)
	{
		color_fg_current = color_fg_current_on;
		color_bg_current = color_bg_current_on;
		m_flash_on = false;
	}
	else
	{
		color_fg_current = color_fg_current_off;
		color_bg_current = color_bg_current_off;
		m_flash_on = true;
	}

	switch (state)
	{
		case LittleDog::STATE_UNINITIALIZED:
			m_label_STATE_UNINITIALIZED->setPaletteForegroundColor(color_fg_current);
			m_label_STATE_UNINITIALIZED->setPaletteBackgroundColor(color_bg_current);
			break;

		case LittleDog::STATE_ROBOT_INITIALIZED:
			m_label_STATE_UNINITIALIZED->setPaletteForegroundColor(color_fg_enabled);
			m_label_STATE_UNINITIALIZED->setPaletteBackgroundColor(color_bg_enabled);

			m_label_STATE_ROBOT_INITIALIZED->setPaletteForegroundColor(color_fg_current);
			m_label_STATE_ROBOT_INITIALIZED->setPaletteBackgroundColor(color_bg_current);
			break;

		case LittleDog::STATE_ROBOT_CALIBRATING:
			m_label_STATE_UNINITIALIZED->setPaletteForegroundColor(color_fg_enabled);
			m_label_STATE_UNINITIALIZED->setPaletteBackgroundColor(color_bg_enabled);

			m_label_STATE_ROBOT_INITIALIZED->setPaletteForegroundColor(color_fg_enabled);
			m_label_STATE_ROBOT_INITIALIZED->setPaletteBackgroundColor(color_bg_enabled);

			m_label_STATE_ROBOT_CALIBRATING->setPaletteForegroundColor(color_fg_current);
			m_label_STATE_ROBOT_CALIBRATING->setPaletteBackgroundColor(color_bg_current);

			m_time_elapsed_idle.start();  // while calibrating, idle should not progress
			break;

		case LittleDog::STATE_ROBOT_CALIBRATED:
			m_label_STATE_UNINITIALIZED->setPaletteForegroundColor(color_fg_enabled);
			m_label_STATE_UNINITIALIZED->setPaletteBackgroundColor(color_bg_enabled);

			m_label_STATE_ROBOT_INITIALIZED->setPaletteForegroundColor(color_fg_enabled);
			m_label_STATE_ROBOT_INITIALIZED->setPaletteBackgroundColor(color_bg_enabled);

			m_label_STATE_ROBOT_CALIBRATING->setPaletteForegroundColor(color_fg_enabled);
			m_label_STATE_ROBOT_CALIBRATING->setPaletteBackgroundColor(color_bg_enabled);

			m_label_STATE_ROBOT_CALIBRATED->setPaletteForegroundColor(color_fg_current);
			m_label_STATE_ROBOT_CALIBRATED->setPaletteBackgroundColor(color_bg_current);
			break;

		case LittleDog::STATE_TRIAL_PLANNING:
			m_label_STATE_UNINITIALIZED->setPaletteForegroundColor(color_fg_enabled);
			m_label_STATE_UNINITIALIZED->setPaletteBackgroundColor(color_bg_enabled);

			m_label_STATE_ROBOT_INITIALIZED->setPaletteForegroundColor(color_fg_enabled);
			m_label_STATE_ROBOT_INITIALIZED->setPaletteBackgroundColor(color_bg_enabled);

			m_label_STATE_ROBOT_CALIBRATING->setPaletteForegroundColor(color_fg_enabled);
			m_label_STATE_ROBOT_CALIBRATING->setPaletteBackgroundColor(color_bg_enabled);

			m_label_STATE_ROBOT_CALIBRATED->setPaletteForegroundColor(color_fg_enabled);
			m_label_STATE_ROBOT_CALIBRATED->setPaletteBackgroundColor(color_bg_enabled);

			m_label_STATE_TRIAL_PLANNING->setPaletteForegroundColor(color_fg_current);
			m_label_STATE_TRIAL_PLANNING->setPaletteBackgroundColor(color_bg_current);
			break;

		case LittleDog::STATE_TRIAL_RUNNING:
			m_label_STATE_UNINITIALIZED->setPaletteForegroundColor(color_fg_enabled);
			m_label_STATE_UNINITIALIZED->setPaletteBackgroundColor(color_bg_enabled);

			m_label_STATE_ROBOT_INITIALIZED->setPaletteForegroundColor(color_fg_enabled);
			m_label_STATE_ROBOT_INITIALIZED->setPaletteBackgroundColor(color_bg_enabled);

			m_label_STATE_ROBOT_CALIBRATING->setPaletteForegroundColor(color_fg_enabled);
			m_label_STATE_ROBOT_CALIBRATING->setPaletteBackgroundColor(color_bg_enabled);

			m_label_STATE_ROBOT_CALIBRATED->setPaletteForegroundColor(color_fg_enabled);
			m_label_STATE_ROBOT_CALIBRATED->setPaletteBackgroundColor(color_bg_enabled);

			m_label_STATE_TRIAL_PLANNING->setPaletteForegroundColor(color_fg_enabled);
			m_label_STATE_TRIAL_PLANNING->setPaletteBackgroundColor(color_bg_enabled);

			m_label_STATE_TRIAL_RUNNING->setPaletteForegroundColor(color_fg_current);
			m_label_STATE_TRIAL_RUNNING->setPaletteBackgroundColor(color_bg_current);
			break;

		case LittleDog::STATE_STOPPING_ROBOT:
			m_label_STATE_UNINITIALIZED->setPaletteForegroundColor(color_fg_enabled);
			m_label_STATE_UNINITIALIZED->setPaletteBackgroundColor(color_bg_enabled);

			m_label_STATE_ROBOT_INITIALIZED->setPaletteForegroundColor(color_fg_enabled);
			m_label_STATE_ROBOT_INITIALIZED->setPaletteBackgroundColor(color_bg_enabled);

			m_label_STATE_STOPPING_ROBOT->setPaletteForegroundColor(color_fg_current);
			m_label_STATE_STOPPING_ROBOT->setPaletteBackgroundColor(color_bg_current);
			break;

		case LittleDog::STATE_ABORTING_CONNECTION:
			m_label_STATE_UNINITIALIZED->setPaletteForegroundColor(color_fg_abort);
			m_label_STATE_UNINITIALIZED->setPaletteBackgroundColor(color_bg_abort);

			m_label_STATE_ROBOT_INITIALIZED->setPaletteForegroundColor(color_fg_abort);
			m_label_STATE_ROBOT_INITIALIZED->setPaletteBackgroundColor(color_bg_abort);

			m_label_STATE_ROBOT_CALIBRATING->setPaletteForegroundColor(color_fg_abort);
			m_label_STATE_ROBOT_CALIBRATING->setPaletteBackgroundColor(color_bg_abort);

			m_label_STATE_ROBOT_CALIBRATED->setPaletteForegroundColor(color_fg_abort);
			m_label_STATE_ROBOT_CALIBRATED->setPaletteBackgroundColor(color_bg_abort);

			m_label_STATE_TRIAL_PLANNING->setPaletteForegroundColor(color_fg_abort);
			m_label_STATE_TRIAL_PLANNING->setPaletteBackgroundColor(color_bg_abort);

			m_label_STATE_TRIAL_RUNNING->setPaletteForegroundColor(color_fg_abort);
			m_label_STATE_TRIAL_RUNNING->setPaletteBackgroundColor(color_bg_abort);

			m_label_STATE_STOPPING_ROBOT->setPaletteForegroundColor(color_fg_abort);
			m_label_STATE_STOPPING_ROBOT->setPaletteBackgroundColor(color_bg_abort);

			m_label_STATE_ABORTING_CONNECTION->setPaletteForegroundColor(color_fg_current);
			m_label_STATE_ABORTING_CONNECTION->setPaletteBackgroundColor(color_bg_current);
			break;
	}

	bool enable_goal = m_dog.getMocapInitialized() &&
		m_overhead_view->get_goal() &&
		( (state == LittleDog::STATE_TRIAL_PLANNING) || 
			(state == LittleDog::STATE_TRIAL_RUNNING) || 
			(state == LittleDog::STATE_ROBOT_INITIALIZED) ||
			(state == LittleDog::STATE_ROBOT_CALIBRATED) );

	m_distance_label->setEnabled(enable_goal);
	m_lcd_goal_distance->setEnabled(enable_goal);
}


/****************************************************************************/
/*
 *  This function updates readings when communications have
 *   been established with the robot.  It is called periodically
 *   when m_timer_update_readings times out.
 */
void
LittleDogUI::slot_update_readings()
{
	LittleDog::State state = m_dog.getState();

	if ( (state != m_last_state) && (state == LittleDog::STATE_TRIAL_RUNNING) )
	{
		stop_planning_timer();
		m_last_elapsed_trial_plan_duration = m_time_elapsed_trial_plan.elapsed();

		m_last_elapsed_trial_run_duration = 0;
		m_time_elapsed_trial_run.start();
	}

	if ( (state == LittleDog::STATE_UNINITIALIZED) ||
		(state == LittleDog::STATE_ABORTING_CONNECTION)	)
	{
		m_last_elapsed_trial_run_duration = 0;
		m_last_elapsed_trial_plan_duration = 0;
		m_last_elapsed_idle_duration = 0;
		return;
	}

	m_lcd_battery_voltage->display("0");
	m_lcd_battery_current->display("0");
	m_lcd_cpu_temp->display("0");
	m_lcd_data_timestamp->display("0");
	m_lcd_packet_loss->display("0");
	m_lcd_trial_run_time->display("0");
	m_lcd_trial_plan_time->display("0");
	m_lcd_idle_time->display("0");

	float battery_voltage; m_dog.getBatteryVoltage(&battery_voltage);
	float battery_current; m_dog.getBatteryCurrent(&battery_current);
	float cpu_temp;        m_dog.getCPUTemperature(&cpu_temp);
	float data_timestamp;  m_dog.getDataTimestamp(&data_timestamp);
	float packet_loss;     m_dog.getPacketLoss(&packet_loss);

	float trial_run_time   = (state == LittleDog::STATE_TRIAL_RUNNING) ?
		0.001f * m_time_elapsed_trial_run.elapsed() : 
		0.001f * m_last_elapsed_trial_run_duration;
	float trial_plan_time   = (state == LittleDog::STATE_TRIAL_PLANNING) ?
		0.001f * m_time_elapsed_trial_plan.elapsed() : 
		0.001f * m_last_elapsed_trial_plan_duration;
	float idle_time    = ((state == LittleDog::STATE_ROBOT_INITIALIZED) 
		|| (state == LittleDog::STATE_ROBOT_CALIBRATED)) ?
		0.001f * m_time_elapsed_idle.elapsed() : 
		0.001f * m_last_elapsed_idle_duration;

	m_lcd_battery_voltage->display(QString::number(battery_voltage, 'f', 2));
	m_lcd_battery_current->display(QString::number(battery_current, 'f', 2));
	m_lcd_cpu_temp->display(QString::number(cpu_temp, 'f', 2));
	m_lcd_data_timestamp->display(QString::number(data_timestamp, 'f', 1));
	m_lcd_packet_loss->display(QString::number(packet_loss, 'f', 1));
	m_lcd_trial_run_time->display(QString::number(trial_run_time, 'f', 0));
	m_lcd_trial_plan_time->display(QString::number(trial_plan_time, 'f', 0));
	m_lcd_idle_time->display(QString::number(idle_time, 'f', 0));

	/*
	 *  Change color of voltage and packet loss widget if critical
	 */
	if (battery_voltage >= MIN_ACCEPTABLE_VOLTAGE)
		m_lcd_battery_voltage->setPaletteForegroundColor(QColor(85, 255, 0));
	else
		m_lcd_battery_voltage->setPaletteForegroundColor(QColor(255, 0, 0));

	if (packet_loss < MAX_ACCEPTABLE_PACKET_LOSS)
		m_lcd_packet_loss->setPaletteForegroundColor(QColor(85, 255, 0));
	else
		m_lcd_packet_loss->setPaletteForegroundColor(QColor(255, 0, 0));
	
	/*
	 *  Check to see if CPU temperature has gone too high.  Once
	 *   it has, the CPU temp widget will stay red, even if it
	 *   drops back down.
	 */
	if (!m_temp_went_too_high &&
		(cpu_temp >= MAX_ACCEPTABLE_TEMP))
	{
		m_temp_went_too_high = TRUE;
		m_lcd_cpu_temp->setPaletteForegroundColor(QColor(255, 0, 0));
	}

	m_last_state = state;
}


/****************************************************************************/
/*
 *  This function is called when the m_timer_planning_exceeded timer expires.
 *   It moves the LittleDog state machine into STATE_TRIAL_RUNNING.
 */
void
LittleDogUI::slot_planning_exceeded()
{
	bdu_log_print(BDU_LOG_LEVEL_WARN, 
		"*** Planning duration exceeded! ... calling donePlanning()\n" );

	LD_ERROR_CODE err = m_dog.donePlanning();

	if (err != LD_OKAY)
		bdu_log_printf(BDU_LOG_LEVEL_ERROR, "ERROR:  donePlanning() %s\n",
			m_dog.getErrorCodeString(err));
}


/****************************************************************************/
/*
 *  This function checks to see if there's something that's supposed to
 *   be append to the log.  If so, get the m_log_append_string mutex
 *   and do it.
 */
void
LittleDogUI::slot_update_log()
{
	m_timer_update_log->stop();
	m_log_append_string_mutex.lock();

	if (m_log_append_string != "")
	{
		m_edit_log->append(m_log_append_string);
		m_log_append_string = "";
	}

	m_log_append_string_mutex.unlock();
	m_timer_update_log->start(TIMER_UPDATE_LOG_PERIOD);
}


/****************************************************************************/
/*
 *  This function should be overridden in a subclass if the overhead
 *   view is to be useful.
 *
 *  This default version just updates the distance, checkbox, and draws the
 *  goal circle in the correct color.
 */
void
LittleDogUI::slot_update_overhead_view()
{
	m_lcd_goal_distance->display("0");

	/*
	 * if we are connected to mocap and the robot, then update the goal distance
	 * and check whether the goal has been reached.
	 */

	LittleDog::State state = m_dog.getState();

	if ( m_dog.getMocapInitialized() && 
		( (state == LittleDog::STATE_TRIAL_PLANNING) || 
			(state == LittleDog::STATE_TRIAL_RUNNING) || 
			(state == LittleDog::STATE_ROBOT_INITIALIZED) ||
			(state == LittleDog::STATE_ROBOT_CALIBRATED)) )
	{
		float distance = 100.0 * m_overhead_view->get_distance_to_goal(m_dog);  // m -> cm
		
		m_lcd_goal_distance->display(QString::number(distance, 'f', 1));
		
		if ( distance <= 5.0 )
		{
			bool was_checked = m_check_goal_reached->isChecked();
			m_check_goal_reached->setChecked(TRUE);
			
			if ( !was_checked )
				generate_trial_summary();

			goal_reached();
		}

		/*
		 *  Draw the goal red if we've never reached it and green otherwise.
		 */

		float goal_x, goal_y;
		m_overhead_view->get_goal(&goal_x, &goal_y);

		if ( m_check_goal_reached->isChecked() )
			m_overhead_view->draw_circle( goal_x, goal_y, 0.05, 0.05, 1, 0, 255, 0 );
		else
			m_overhead_view->draw_circle( goal_x, goal_y, 0.05, 0.05, 1, 255, 0, 0 );
	}
}


/****************************************************************************/
/*
 *  OverHeadView widget class
 */
/****************************************************************************/


/****************************************************************************/
OverheadViewWidget::OverheadViewWidget(QWidget* parent, const char* name) :
	QWidget(parent, name),
	m_qpainter(NULL),
	m_goal_file_state(GOAL_DATA_UNINITIALIZED)
{
	m_center_x = 0.0f;
	m_center_y = 0.0f;
	m_goal_x = 0.0f;
	m_goal_y = 0.0f;
	m_ppm = 100;
	m_start_x = 0.0f;
	m_start_y = 0.0f;
	m_start_pos_valid = false;
}


/****************************************************************************/
OverheadViewWidget::~OverheadViewWidget()
{
	if (m_qpainter)
		delete m_qpainter;
}


/****************************************************************************/
bool
OverheadViewWidget::get_goal(float* goal_x, float* goal_y)
{
	switch ( m_goal_file_state )
	{
	case GOAL_DATA_UNINITIALIZED:
		/*
		 *   Attempt to load the goal file and learn the goal.
		 */
		{
			/*
			 * compute the filename
			 */

			const char* LITTLEDOG_env_var = getenv("LITTLEDOG");
			
			if ( LITTLEDOG_env_var == NULL )
			{
				bdu_log_printf(BDU_LOG_LEVEL_ERROR, 
					"ERROR: get_goal_location() can't getenv(LITTLEDOG).\n" );
				m_goal_file_state = GOAL_DATA_MISSING;
				return false;
			}
			
			char data_file_name[512];
			sprintf( data_file_name, "%s/test/goal.txt", LITTLEDOG_env_var );
			
			/*
			 * open the file for text reading
			 */

			FILE* goal_file = fopen( data_file_name, "rt" );

			if ( goal_file == NULL )
			{
				bdu_log_printf(BDU_LOG_LEVEL_ERROR, 
					"ERROR: get_goal() can't open file \"%s\".\n", data_file_name );
				m_goal_file_state = GOAL_DATA_MISSING;
				return false;
			}

			/*
			 * scan file content into m_goal_x and m_goal_y (meters)
			 */

			if ( fscanf(goal_file, "%f %f", &m_goal_x, &m_goal_y) != 2 )
			{
				bdu_log_printf(BDU_LOG_LEVEL_ERROR, 
					"ERROR: get_goal() file \"%s\" has invalid contents.\n", data_file_name );
				m_goal_file_state = GOAL_DATA_MISSING;
				fclose(goal_file);
				return false;
			}

			m_goal_file_state = GOAL_DATA_ACQUIRED;
			fclose(goal_file);

			// fall through returning the data
		}

	case GOAL_DATA_ACQUIRED:
		if ( goal_x ) 
			*goal_x = m_goal_x;
		if ( goal_y )
			*goal_y = m_goal_y;
		return true;

	default:
	case GOAL_DATA_MISSING:
		if ( goal_x )
			*goal_x = 0.0;
		if ( goal_y ) 
			*goal_y = 0.0;
		return false;
	}
}

/****************************************************************************/
float 
OverheadViewWidget::get_distance_to_goal(const LittleDog& dog)
{
	/*
	 *  This function returns the distance of the robot's body center to
	 *  the goal location, in meters.
	 */

	if ( m_goal_file_state == GOAL_DATA_ACQUIRED )
	{
		float x, y;
		if ( get_position(dog, &x,&y) == true )
		{
			double dx = fabs(x - m_goal_x);
			double dy = fabs(y - m_goal_y);
			return sqrt( (dx*dx) + (dy*dy));
		}
	}

	// The tragedy of life lies in having no goals to reach.
	return 0.0f;
}


/****************************************************************************/
bool
OverheadViewWidget::get_position(const LittleDog& dog, float *x, float *y) const
{
	// gets the mocap position of the robot body.  returns true on success

	LittleDog::State state = dog.getState();

	if ( !dog.getMocapInitialized() ||
		( (state != LittleDog::STATE_TRIAL_PLANNING) &&
			(state != LittleDog::STATE_TRIAL_RUNNING) &&
			(state != LittleDog::STATE_ROBOT_INITIALIZED) &&
			(state != LittleDog::STATE_ROBOT_CALIBRATED) && 
			(state != LittleDog::STATE_STOPPING_ROBOT) ) )
		return false;

	int age;
	bduVec3f pos;
	bduVec4f orient;
	LD_ERROR_CODE err = dog.getMocapBodyReading(LittleDog::B_TRUNK, &pos, &orient, &age);
	
	if (err != LD_OKAY)
	{
		bdu_log_printf(BDU_LOG_LEVEL_ERROR, "ERROR:  get_position() %s\n",
			dog.getErrorCodeString(err));
		return false;
	}

	*x = pos.n[0];
	*y = pos.n[1];

	return true;
}


/****************************************************************************/
void
OverheadViewWidget::set_center(float center_x, float center_y)
{
	m_center_x = center_x;
	m_center_y = center_y;

	update_qpainter();
}


/****************************************************************************/
void
OverheadViewWidget::set_pixels_per_meter(int ppm)
{
	m_ppm = ppm;

	update_qpainter();
}


/****************************************************************************/
void
OverheadViewWidget::set_start_position(bool valid, float x, float y)
{
	m_start_pos_valid = valid;
	m_start_x = x;
	m_start_y = y;
}


/****************************************************************************/
bool
OverheadViewWidget::get_start_position(float* x, float* y) const
{
	if ( m_start_pos_valid == false )
		return false;

	*x = m_start_x;
	*y = m_start_y;
	return true;
}


/****************************************************************************/
void
OverheadViewWidget::clear(bool draw_grid, float grid_spacing)
{
	QPainter p(this);
	int w = width();
	int h = height();

	/*
	 *  Clear to bg color.
	 */
	QBrush brush_bg(QColor(0, 0, 0));
	p.fillRect(0, 0, w, h, brush_bg);

	/*
	 *  Update m_qpainter.
	 */
	update_qpainter();

	if (draw_grid)
	{
		/*
		 *  Draw horizontal and vertical grid lines.
		 */
		QColor grid_color(0, 96, 0);

		m_qpainter->setPen(QPen(grid_color, 1));

		static const int GRID_RADIUS = 10000;

		int grid_spacing_pix = (int) (grid_spacing * 1000.0f);

		int i;
		for (i = grid_spacing_pix; i <= GRID_RADIUS; i += grid_spacing_pix)
		{
			m_qpainter->drawLine( i, -GRID_RADIUS,  i, GRID_RADIUS);
			m_qpainter->drawLine(-i, -GRID_RADIUS, -i, GRID_RADIUS);

			m_qpainter->drawLine(-GRID_RADIUS,  i, GRID_RADIUS,  i);
			m_qpainter->drawLine(-GRID_RADIUS, -i, GRID_RADIUS, -i);
		}

		m_qpainter->setPen(QPen(grid_color, 2));

		m_qpainter->drawLine(0, -GRID_RADIUS, 0, GRID_RADIUS);
		m_qpainter->drawLine(-GRID_RADIUS, 0, GRID_RADIUS, 0);
	}
}


/****************************************************************************/
void
OverheadViewWidget::draw_visual(const LittleDogVisual& visual)
{
	if (!visual.visible)
		return;

	if (!m_qpainter)
		update_qpainter();

	m_qpainter->setPen(QPen(QColor(visual.color_r,
		visual.color_g,
		visual.color_b),
		visual.line_width));

	m_qpainter->save();

	int trunk_x_mm = (int) (1000.0f * visual.trunk_x);
	int trunk_y_mm = (int) (1000.0f * visual.trunk_y);
	float yaw_deg = visual.trunk_yaw * RAD2DEG_FLOAT;

	m_qpainter->translate(trunk_x_mm, trunk_y_mm);
	m_qpainter->rotate(yaw_deg);

	switch (visual.shape)
	{
		case VISUAL_SHAPE_ELLIPSE:
			draw_circle(0.0f, 0.0f,
				visual.trunk_length, visual.trunk_width,
				visual.line_width,
				visual.color_r,
				visual.color_g,
				visual.color_b);
			break;

		case VISUAL_SHAPE_RECT:
			draw_rect(0.0f, 0.0f,
				visual.trunk_length, visual.trunk_width,
				visual.line_width,
				visual.color_r,
				visual.color_g,
				visual.color_b);
			break;

		case VISUAL_SHAPE_X:
			draw_x(0.0f, 0.0f,
				visual.trunk_length, visual.trunk_width,
				visual.line_width,
				visual.color_r,
				visual.color_g,
				visual.color_b);
			break;
	}

	draw_line(0.0f, 0.0f,
		visual.trunk_length/2.0f, 0.0f,
		visual.line_width,
		visual.color_r,
		visual.color_g,
		visual.color_b);

	if (visual.feet_in_world_coords)
		m_qpainter->restore();

	int foot;
	for (foot = 0; foot < (int) LittleDog::NUM_LEGS; foot++)
	{
		if (!visual.feet[foot].visible)
			continue; // with for

		m_qpainter->save();
		{
			m_qpainter->setPen(QPen(QColor(visual.feet[foot].color_r,
				visual.feet[foot].color_g,
				visual.feet[foot].color_b),
				visual.feet[foot].line_width));

			int foot_x_mm = (int) (1000.0f * visual.feet[foot].foot_x);
			int foot_y_mm = (int) (1000.0f * visual.feet[foot].foot_y);

			m_qpainter->translate(foot_x_mm, foot_y_mm);

			if (visual.feet_in_world_coords)
				m_qpainter->rotate(yaw_deg);

			switch (visual.feet[foot].shape)
			{
				case VISUAL_SHAPE_ELLIPSE:
					draw_circle(0.0f, 0.0f,
						visual.foot_length, visual.foot_width,
						visual.feet[foot].line_width,
						visual.feet[foot].color_r,
						visual.feet[foot].color_g,
						visual.feet[foot].color_b);
					break;

				case VISUAL_SHAPE_RECT:
					draw_rect(0.0f, 0.0f,
						visual.foot_length, visual.foot_width,
						visual.feet[foot].line_width,
						visual.feet[foot].color_r,
						visual.feet[foot].color_g,
						visual.feet[foot].color_b);
					break;

				case VISUAL_SHAPE_X:
					draw_x(0.0f, 0.0f,
						visual.foot_length, visual.foot_width,
						visual.feet[foot].line_width,
						visual.feet[foot].color_r,
						visual.feet[foot].color_g,
						visual.feet[foot].color_b);
					break;
			}
		}
		m_qpainter->restore();
	}

	if (!visual.feet_in_world_coords)
		m_qpainter->restore();
}


/****************************************************************************/
void
OverheadViewWidget::draw_x(float x,
	float y,
	float height,
	float width,
	int line_width,
	int color_r,
	int color_g,
	int color_b)
{
	if (!m_qpainter)
		update_qpainter();

	m_qpainter->setPen(QPen(QColor(color_r, color_g, color_b), line_width));

	m_qpainter->save();
	{
		int x_mm = (int) (1000.0f * x);
		int y_mm = (int) (1000.0f * y);
		int width_mm = (int) (1000.0f * width);
		int height_mm = (int) (1000.0f * height);

		m_qpainter->drawLine(x_mm - height_mm/2, y_mm - width_mm/2,
			x_mm + height_mm/2, y_mm + width_mm/2);
		m_qpainter->drawLine(x_mm + height_mm/2, y_mm - width_mm/2,
			x_mm - height_mm/2, y_mm + width_mm/2);
	}
	m_qpainter->restore();
}


/****************************************************************************/
void
OverheadViewWidget::draw_circle(float x,
	float y,
	float height,
	float width,
	int line_width,
	int color_r,
	int color_g,
	int color_b)
{
	if (!m_qpainter)
		update_qpainter();

	m_qpainter->setPen(QPen(QColor(color_r, color_g, color_b), line_width));

	m_qpainter->save();
	{
		int x_mm = (int) (1000.0f * x);
		int y_mm = (int) (1000.0f * y);
		int width_mm = (int) (1000.0f * width);
		int height_mm = (int) (1000.0f * height);

		m_qpainter->drawEllipse(x_mm - height_mm/2, y_mm - width_mm/2,
			height_mm, width_mm);
	}
	m_qpainter->restore();
}


/****************************************************************************/
void
OverheadViewWidget::draw_rect(float x,
	float y,
	float height,
	float width,
	int line_width,
	int color_r,
	int color_g,
	int color_b)
{
	if (!m_qpainter)
		update_qpainter();

	m_qpainter->setPen(QPen(QColor(color_r, color_g, color_b), line_width));

	m_qpainter->save();
	{
		int x_mm = (int) (1000.0f * x);
		int y_mm = (int) (1000.0f * y);
		int width_mm = (int) (1000.0f * width);
		int height_mm = (int) (1000.0f * height);

		m_qpainter->drawRect(x_mm - height_mm/2, y_mm - width_mm/2,
			height_mm, width_mm);
	}
	m_qpainter->restore();
}


/****************************************************************************/
void
OverheadViewWidget::draw_line(float x0,
	float y0,
	float x1,
	float y1,
	int line_width,
	int color_r,
	int color_g,
	int color_b)
{
	if (!m_qpainter)
		update_qpainter();

	m_qpainter->setPen(QPen(QColor(color_r, color_g, color_b), line_width));

	m_qpainter->save();
	{
		int x0_mm = (int) (1000.0f * x0);
		int y0_mm = (int) (1000.0f * y0);
		int x1_mm = (int) (1000.0f * x1);
		int y1_mm = (int) (1000.0f * y1);

		m_qpainter->drawLine(x0_mm, y0_mm, x1_mm, y1_mm);
	}
	m_qpainter->restore();
}


/****************************************************************************/
void
OverheadViewWidget::set_font(const char* family, int point_size)
{
	if (!m_qpainter)
		update_qpainter();

	m_qpainter->setFont(QFont(family, point_size));
}


/****************************************************************************/
void
OverheadViewWidget::draw_text_widget_coords(int x,
	int y,
	const char* string,
	int color_r,
	int color_g,
	int color_b)
{
	if (!m_qpainter)
		update_qpainter();

	QPainter p(this);

	p.setPen(QPen(QColor(color_r, color_g, color_b), 1));
	p.setFont(m_qpainter->font());

	p.drawText(x, y, string);
}


/****************************************************************************/
void
OverheadViewWidget::draw_text_view_coords(float x,
	float y,
	const char* string,
	int color_r,
	int color_g,
	int color_b)
{
	if (!m_qpainter)
		update_qpainter();

	QPainter p(this);

	p.setPen(QPen(QColor(color_r, color_g, color_b), 1));
	p.setFont(m_qpainter->font());

	int widget_x;
	int widget_y;
	convert_view_coord_to_widget_coord(x, y, &widget_x, &widget_y);

	int w = width();
	int h = height();

	int center_x_pix = (int) (m_center_x * m_ppm);
	int center_y_pix = (int) (m_center_y * m_ppm);

	int ox = -w/2 - center_y_pix;
	int oy = -h/2 - center_x_pix;

	p.drawText(widget_x - ox, widget_y - oy, string);
}


/****************************************************************************/
void
OverheadViewWidget::convert_view_coord_to_widget_coord(float view_x, float view_y,
	int* widget_x, int* widget_y)
{
	/*
	 *  This transformation assumes that X is up, Y is to the left.
	 */
	*widget_x = (int) (-view_y * m_ppm);
	*widget_y = (int) (-view_x * m_ppm);
}


/****************************************************************************/
void
OverheadViewWidget::get_widget_size(int* w, int* h)
{
	*w = width();
	*h = height();
}


/****************************************************************************/
void
OverheadViewWidget::update_qpainter()
{
	/*
	 *  Create a new QPainter.
	 */
	if (m_qpainter)
		delete m_qpainter;
	m_qpainter = new QPainter(this);

	/*
	 *  Set the center of the view.
	 */
	int w = width();
	int h = height();

	int center_x_pix = (int) m_center_x * m_ppm;
	int center_y_pix = (int) m_center_y * m_ppm;

	int ox = -w/2 - center_y_pix;
	int oy = -h/2 - center_x_pix;
	m_qpainter->setWindow(ox, oy, w, h);

	/*
	 *  Set transformation matrix of painter.
	 *
	 *  Use a scale to get positive Y coords going up instead of down,
	 *   and to change units to mm.
	 *
	 *  Use a rotate to make positive X coords go up.
	 */
	m_qpainter->scale(m_ppm / 1000.0, -m_ppm / 1000.0);
	m_qpainter->rotate(90.0);
}


/****************************************************************************/
QPainter*
OverheadViewWidget::get_qpainter()
{
	if (!m_qpainter)
		update_qpainter();

	return m_qpainter;
}

