
#include <math.h>

#include <littledog.h>

#include "MyLittleDogUI.h"

static const int DRAW_FAKE_DESIRED_DATA = 0;


/****************************************************************************/
MyLittleDogUI::MyLittleDogUI(MyLittleDog& dog) :
	LittleDogUI(dog,-1),
	m_vis_current(TRUE, 0, 255, 0),
	m_vis_desired(TRUE, 255, 255, 0)
{
	m_overhead_view->set_pixels_per_meter(200);
}


/****************************************************************************/
void
MyLittleDogUI::slot_update_overhead_view()
{
	/*
	 *  Clear view and draw grid.
	 */
	m_overhead_view->clear(TRUE, 0.5f);

	/*
	 *  Call the base class since it updates goal distance information.
	 */
	LittleDogUI::slot_update_overhead_view();

	/*
	 *  Draw what we want to draw
	 */
	m_overhead_view->set_font("Times", 12);
	m_overhead_view->draw_text_widget_coords(20, 20,
		"LittleDog Example walk",
		255, 255, 255);

	if (!m_dog.getMocapInitialized())
	{
		m_overhead_view->draw_text_widget_coords(20, 40,
			"mocap not initialized...",
			255, 255, 255);

		return;
	}

	if (DRAW_FAKE_DESIRED_DATA)
	{
		/*
		 *  (For now test desired visuals by manually moving them,
		 *   instead of asking the robot for information.)
		 */
		m_vis_desired.line_width = 2;
		m_vis_desired.trunk_yaw += (3.0 * DEG2RAD_FLOAT);

		m_vis_desired.trunk_x =  sin(m_vis_desired.trunk_yaw);
		m_vis_desired.trunk_y = -cos(m_vis_desired.trunk_yaw);

		m_vis_desired.feet[LittleDog::FL].foot_x =  0.1f + 0.05f * sin(m_vis_desired.trunk_yaw / (4.0 * DEG2RAD_FLOAT));
		m_vis_desired.feet[LittleDog::FR].foot_x =  0.1f - 0.05f * sin(m_vis_desired.trunk_yaw / (4.0 * DEG2RAD_FLOAT));
		m_vis_desired.feet[LittleDog::HL].foot_x = -0.1f - 0.05f * sin(m_vis_desired.trunk_yaw / (4.0 * DEG2RAD_FLOAT));
		m_vis_desired.feet[LittleDog::HR].foot_x = -0.1f + 0.05f * sin(m_vis_desired.trunk_yaw / (4.0 * DEG2RAD_FLOAT));

		/*
		 *  Draw some random Xs and Os.
		 */
		m_overhead_view->draw_x(0.2f, 0.2f, 0.2f, 0.1f, 1,
			0, 85, 255);
		m_overhead_view->draw_circle(0.2f, -0.4f, 0.2f, 0.3f, 2,
			0, 85, 255);
		m_overhead_view->draw_line(0.2f, 0.2f, 0.2f, -0.4f, 1,
			0, 255, 255);

		/*
		 *  Draw the artificial trajectory m_vis_desired is following.
		 */
		m_overhead_view->draw_circle(0.0f, 0.0f, 2.0f, 2.0f, 2,
			128, 0, 0);
	}

	/*
	 *  Draw body only?
	 */
	bool body_only = false;
	if (strcmp(m_dog.getMocapModelName(),"littledog_body")==0)
		body_only = true;

	/*
	 *  Read mocap data.
	 */
	LittleDog::MocapReadInfo read_info;

	LD_ERROR_CODE result = m_dog.getBulkMocapReadings(&read_info);
	if (result != LD_OKAY)
	{
		bdu_log_printf(BDU_LOG_LEVEL_DEBUG,
			"An error occurred calling function '%s()'.\n"
			"    Error was: '%s'.\n", 
			"getBulkMocapReadings",
			m_dog.getErrorCodeString(result));
	}
	else
	{
		LittleDog::MocapBodyInfo* body;

		/*
		 *  Set m_vis_current trunk position and orientation.
		 */
		body = &read_info.m_bodies[LittleDog::B_TRUNK];
		m_vis_current.trunk_x = body->m_position.n[0];
		m_vis_current.trunk_y = body->m_position.n[1];
		m_vis_current.trunk_yaw = body->m_orientation.n[0];

		/*
		 *  Set m_vis_current feet positions.
		 */
		body = &read_info.m_bodies[LittleDog::B_FL_LLEG];
		m_vis_current.feet[LittleDog::FL].foot_x = body->m_position.n[0];
		m_vis_current.feet[LittleDog::FL].foot_y = body->m_position.n[1];
		m_vis_current.feet[LittleDog::FL].visible = !body_only;

		body = &read_info.m_bodies[LittleDog::B_FR_LLEG];
		m_vis_current.feet[LittleDog::FR].foot_x = body->m_position.n[0];
		m_vis_current.feet[LittleDog::FR].foot_y = body->m_position.n[1];
		m_vis_current.feet[LittleDog::FR].visible = !body_only;

		body = &read_info.m_bodies[LittleDog::B_HL_LLEG];
		m_vis_current.feet[LittleDog::HL].foot_x = body->m_position.n[0];
		m_vis_current.feet[LittleDog::HL].foot_y = body->m_position.n[1];
		m_vis_current.feet[LittleDog::HL].visible = !body_only;

		body = &read_info.m_bodies[LittleDog::B_HR_LLEG];
		m_vis_current.feet[LittleDog::HR].foot_x = body->m_position.n[0];
		m_vis_current.feet[LittleDog::HR].foot_y = body->m_position.n[1];
		m_vis_current.feet[LittleDog::HR].visible = !body_only;

		m_vis_current.feet_in_world_coords = 1;

		/*
		 *  Draw 'x' at the positions of all bodies.
		 */
		int i, max_i = LittleDog::NUM_BODIES;

		if (body_only)
			max_i = 1;

		for (i=0; i<max_i; i++)
		{
			body = &read_info.m_bodies[i];
			m_overhead_view->draw_x(body->m_position.n[0], body->m_position.n[1],
				0.05f, 0.05f, 1,
				255, 255, 255);
		}

		///*
		// *  Draw positions of all markers.
		// */
		//int i;
		//for (i=0; i<read_info.m_num_markers; i++)
		//{
		//	LittleDog::MocapMarkerInfo* marker = &read_info.m_markers[i];
		//	m_overhead_view->draw_x(marker->m_position.n[0], marker->m_position.n[1],
		//		0.05f, 0.05f, 1,
		//		0, 85, 255);
		//}
	}

	/*
	 *  Draw dog visuals.
	 */
	if (DRAW_FAKE_DESIRED_DATA)
		m_overhead_view->draw_visual(m_vis_desired);
	m_overhead_view->draw_visual(m_vis_current);
}
