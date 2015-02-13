
#ifndef  __LittleDogUI_H
#define  __LittleDogUI_H

//
//  Qt includes.
//
#include <qmutex.h>
#include <qstring.h>
#include <qwidget.h>
#include <qdatetime.h>

//
//  LittleDog API includes.
//
#include <littledog.h>
#include <bduLog.h>

//
//  Local includes.
//
#include "LittleDogUI_uic.h"

//
//  Forward declarations.
//
class QTimer;

class LittleDogTrialThread;
class LittleDogCalibrateThread;
class OverheadViewWidget;

static const float PI_FLOAT = 3.1415927f;
static const float DEG2RAD_FLOAT = PI_FLOAT / 180.0f;
static const float RAD2DEG_FLOAT = 180.0f / PI_FLOAT;


/****************************************************************************/
class LittleDogUI : public LittleDogUI_uic
{
	Q_OBJECT
    
public:

	/*
	 *  The constructor for LittleDogUI accepts a LittleDog and a maximum 
	 *  allowed duration for planning.  After this duration (in seconds) of 
	 *  planning is reached the UI will automatically call donePlanning().
	 *  This feature is disabled when the max planning duration is set to -1.
	 */

	LittleDogUI(LittleDog& dog, int max_planning_duration_sec = -1);
	virtual ~LittleDogUI();

	/*
	 *  goal_reached() is automaticalaly called repeatedly by the GUI while
	 *  the robot's body crosses into the goal cylinder during a trial with 
	 *  mocap enabled and having valid goal data.  This function should be 
	 *  overridden in a subclass that uses the goal.
	 */
	virtual void goal_reached() {}

	/*
	 *  This helper function uses the results of get_goal() and the mocap
	 *  location of the LittleDog body to compute a distance, in meters,
	 *  between the two.  If get_goal() returns false, false is returned.
	 */
	bool get_goal_reached() const;

	void queue_to_log_widget(int l, const char* string);
	void wait_for_calibrate_thread(bool force_termination);
	void wait_for_trial_thread(bool force_termination);
	void generate_trial_summary();

	/*
	 *  These helper functions control the planning timer, if it has been
	 *  enabled - otherwise they do nothing.
	 */
	void start_planning_timer();
	void stop_planning_timer();

protected slots:

	virtual void slot_clicked_connect_to_mocap_system();
	virtual void slot_clicked_disconnect_from_mocap_system();
	virtual void slot_clicked_connect_to_robot();
	virtual void slot_clicked_calibrate_robot();
	virtual void slot_clicked_run_trial();
	virtual void slot_clicked_done_planning();
	virtual void slot_clicked_stop_robot();
	virtual void slot_clicked_quit();
	virtual void slot_clicked_abort_connection();

	virtual void slot_pressed_button();
	virtual void slot_released_button();

	virtual void closeEvent(QCloseEvent* e);

	virtual void slot_update_log();
	virtual void slot_update_buttons_and_state_labels();
	virtual void slot_update_readings();
	virtual void slot_planning_exceeded();

	/*
	 *  This function should be overridden in a subclass if the overhead
	 *   view is to be useful.
	 */
	virtual void slot_update_overhead_view();

protected:

	void start_calibrate_thread();
	void stop_calibrate_thread();

	void start_trial_thread();
	void stop_trial_thread();

	OverheadViewWidget*   m_overhead_view;
	QTimer*               m_timer_update_overhead_view;

	QTimer*               m_timer_update_log;
	bduLog*               m_log_callback;
	QMutex                m_log_append_string_mutex;
	QString               m_log_append_string;

	LittleDog&            m_dog;
	LittleDog::State      m_last_state;

	LittleDogTrialThread*     m_trial_thread;
	LittleDogCalibrateThread* m_calibrate_thread;

	QTimer*               m_timer_update_readings;
	QTimer*               m_timer_planning_exceeded;

	QTime                 m_time_elapsed_trial_run;
	QTime                 m_time_elapsed_trial_plan;
	QTime                 m_time_elapsed_idle;

	int                   m_last_elapsed_trial_run_duration;
	int                   m_last_elapsed_trial_plan_duration;
	int                   m_last_elapsed_idle_duration;

	int                   m_max_planning_duration_sec;

	bool                  m_flash_on;
	bool                  m_temp_went_too_high;
};



/****************************************************************************/
/*
 *  LittleDogVisual and LittleDogFootVisual are helper classes
 *   that make it easy to draw a set of shapes that represent
 *   a LittleDog.  It is possible to set the color, trunk position,
 *   and feet positions of these objects.
 */
/****************************************************************************/


/****************************************************************************/
/*
 *  Options for shape of various parts of LittleDogVisual.
 */
enum {
	VISUAL_SHAPE_ELLIPSE = 0,
	VISUAL_SHAPE_RECT,
	VISUAL_SHAPE_X
};


/****************************************************************************/
class LittleDogFootVisual
{
public:

	LittleDogFootVisual() :
		visible(true),
		shape(VISUAL_SHAPE_RECT),
		line_width(1),
		color_r(255),
		color_g(255),
		color_b(255),
		foot_x(0.0f),
		foot_y(0.0f)
		{};

	bool  visible;       // whether foot is visible

	int   shape;
	int   line_width;    // width of lines used to draw foot
	int   color_r;       // color of foot
	int   color_g;       // color of foot
	int   color_b;       // color of foot

	float foot_x;  // in meters
	float foot_y;  // in meters
};


/****************************************************************************/
/*
 *  NOTE: Bacause of QPainter considerations, units of length are in
 *   millimeters rather than meters.
 */
class LittleDogVisual
{
public:

	LittleDogVisual(bool vis, int r, int g, int b) :
		visible(vis),
		shape(VISUAL_SHAPE_RECT),
		line_width(1),
		color_r(r),
		color_g(g),
		color_b(b),
		trunk_length(0.33f),  // rough approximation
		trunk_width (0.08f),  // rough approximation
		trunk_x(0.0f),
		trunk_y(0.0f),
		trunk_yaw(0.0f),
		foot_length(0.08f),   // rough approximation
		foot_width(0.04f),    // rough approximation
		feet_in_world_coords(false)
	{
		feet[LittleDog::FL].foot_x = + 0.12f;
		feet[LittleDog::FL].foot_y = - 0.07f;
		feet[LittleDog::FL].color_r = color_r;
		feet[LittleDog::FL].color_g = color_g;
		feet[LittleDog::FL].color_b = color_b;

		feet[LittleDog::FR].foot_x = + 0.12f;
		feet[LittleDog::FR].foot_y = + 0.07f;
		feet[LittleDog::FR].color_r = color_r;
		feet[LittleDog::FR].color_g = color_g;
		feet[LittleDog::FR].color_b = color_b;

		feet[LittleDog::HL].foot_x = - 0.12f;
		feet[LittleDog::HL].foot_y = - 0.07f;
		feet[LittleDog::HL].color_r = color_r;
		feet[LittleDog::HL].color_g = color_g;
		feet[LittleDog::HL].color_b = color_b;

		feet[LittleDog::HR].foot_x = - 0.12f;
		feet[LittleDog::HR].foot_y = + 0.07f;
		feet[LittleDog::HR].color_r = color_r;
		feet[LittleDog::HR].color_g = color_g;
		feet[LittleDog::HR].color_b = color_b;
	}

	bool  visible;       // whether entire visual (including feet) is visible

	int   shape;
	int   line_width;    // width of lines used to draw trunk
	int   color_r;       // color of trunk
	int   color_g;       // color of trunk
	int   color_b;       // color of trunk

	float trunk_length;  // in meters
	float trunk_width;   // in meters

	float trunk_x;       // world coords, in meters
	float trunk_y;       // world coords, in meters
	float trunk_yaw;     // world coords, in radians

	float foot_length;   // in meters
	float foot_width;    // in meters

	LittleDogFootVisual feet[LittleDog::NUM_LEGS];  // feet objects
	bool                feet_in_world_coords;       // whether foot positions are in world or body coord
};


/****************************************************************************/
class OverheadViewWidget : public QWidget
{
	Q_OBJECT

public:

	OverheadViewWidget(QWidget* parent = 0, const char* name = 0);
	~OverheadViewWidget();

	/*
	 *  This function will consult the $LITTLEDOG/test/goal.txt file to learn
	 *  the goal location, and will return that location, converted to units of
	 *  meters from the mocap origin.  "true" is returned if the file was able
	 *  to be located and read successfully.
	 */
	bool get_goal(float* goal_x = NULL, float* goal_y = NULL);

	/*
	 *  This helper function reads the mocap position of the Littledog body.
	 *  True is returned if the request is valid and the (x,y) are populated.
	 */
	bool get_position(const LittleDog& dog, float *x, float *y) const;

	/*
	 *  This helper function uses the results of get_goal() and the mocap
	 *  location of the LittleDog body to compute a distance, in meters,
	 *  between the two.  If get_goal() returns false, 0 is returned.
	 */
	float get_distance_to_goal(const LittleDog& dog);

	/*
	 *  This helper function returns the position (in meters) of the dog when the 
	 *  "Start Trial" button was pushed, provided the robot was connected to
	 *  the mocap system.  True is returned if the request is valid and the 
	 *  (x,y) are populated.
	 */
	bool get_start_position( float *x, float *y) const;

	/*
	 *  Stores the start x and y values (the mocap position of the robot's body)
	 *  in meters.  Set valid to false if mocap is not available, and pass in
	 *  zero for x and y in that case.  Otherwise, pass in the values from
	 *  get_position() when the "Run Trial"  button is pressed.
	 */
	void set_start_position( bool valid, float x, float y );

	/*
	 *  This function sets where the center of the view
	 *   will be.  Coordinates are in meters.
	 *
	 *  For example, if -1.0, 2.0 are passed as arguments,
	 *   the point (-1.0, 2.0) in world coordinates will
	 *   be in the center of the view.
	 *
	 *  The default values are 0, 0.
	 */
	void set_center(float center_x, float center_y);

	/*
	 *  This function sets how many pixels per meter there
	 *   are in the view.
	 *
	 *  The default value is 100 pixels per meter.
	 */
	void set_pixels_per_meter(int ppm);

	/*
	 *  This function clears the view, and optionally draws
	 *   a grid.
	 */
	void clear(bool draw_grid = TRUE, float grid_spacing = 1.0f);

	/*
	 *  This function draws a LittleDogVisual.
	 */
	void draw_visual(const LittleDogVisual& visual);

	/*
	 *  This function draws an 'X'.
	 */
	void draw_x(float x,  // center in x
		float y,          // center in y
		float height,     // diameter in x
		float width,      // diameter in y
		int line_width,
		int color_r,
		int color_g,
		int color_b);

	/*
	 *  This function draws a circle.
	 */
	void draw_circle(float x, // center in x
		float y,              // center in y
		float height,         // diameter in x
		float width,          // diameter in y
		int line_width,
		int color_r,
		int color_g,
		int color_b);

	/*
	 *  This function draws a rectangle.
	 */
	void draw_rect(float x,  // center in x
		float y,             // center in y 
		float height,        // diameter in x
		float width,         // diameter in y
		int line_width,
		int color_r,
		int color_g,
		int color_b);

	/*
	 *  This function draws a line.
	 */
	void draw_line(float x0,
		float y0,
		float x1,
		float y1,
		int line_width,
		int color_r,
		int color_g,
		int color_b);

	/*
	 *  This function sets the font used in draw_text_*() functions.
	 *   Useful fonts are: "Times", "Courier".
	 */
	void set_font(const char* family, int point_size = 10);

	/*
	 *  This function renders the passed string in widget pixel
	 *   coordinates: (0, 0) in upper left corner.
	 */
	void draw_text_widget_coords(int x,
		int y,
		const char* string,
		int color_r,
		int color_g,
		int color_b);

	/*
	 *  This function renders the passed string in view
	 *   coordinates: x and y are in meters with the same
	 *   origin as draw_circle(), etc.
	 */
	void draw_text_view_coords(float x,
		float y,
		const char* string,
		int color_r,
		int color_g,
		int color_b);

	/*
	 *  This function converts passed widget coordinates
	 *   (in meters with X up, Y to left) to view coordinates
	 *   (in pixels with origin in upper left corner).
	 */
	void convert_view_coord_to_widget_coord(float view_x, float view_y,
		int* widget_x, int* widget_y);

	/*
	 *  This function returns the widget's size in pixels.
	 */
	void get_widget_size(int* width, int* height);

	/*
	 *  Accessor to QPainter object used for drawing,
	 *   if the above primitives are insufficient.
	 *
	 *   NOTE:  Units of the QPainter are in millimeters,
	 *   not meters.
	 */
	QPainter* get_qpainter();

protected:

	void update_qpainter();

	enum GoalFileState {
		GOAL_DATA_UNINITIALIZED,
		GOAL_DATA_ACQUIRED,
		GOAL_DATA_MISSING,
	};

	bool  m_start_pos_valid;
	float m_start_x;
	float m_start_y;
	float m_goal_x;
	float m_goal_y;
	float m_center_x;
	float m_center_y;
	int   m_ppm;

	QPainter* m_qpainter;
	GoalFileState m_goal_file_state;
};


#endif   // __LittleDogUI_H
