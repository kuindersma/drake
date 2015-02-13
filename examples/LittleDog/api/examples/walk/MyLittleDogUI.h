
#ifndef  __MyLittleDogUI_H
#define  __MyLittleDogUI_H

#include <LittleDogUI.h>

#include "MyLittleDog.h"


/****************************************************************************/
class MyLittleDogUI : public LittleDogUI
{
	Q_OBJECT
    
public:

	MyLittleDogUI(MyLittleDog& dog);

protected slots:

	virtual void slot_update_overhead_view();

private:

	/*
	 *  A LittleDogVisual is a helper class that makes
	 *   it easy to draw a set of shapes that represent
	 *   a LittleDog.  It is possible to set the color,
	 *   trunk position, and feet positions of these
	 *   objects.
	 */
	LittleDogVisual m_vis_current;
	LittleDogVisual m_vis_desired;
};

#endif   // __MyLittleDogUI_H
