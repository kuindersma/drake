
#ifndef  __MyLittleDog_H
#define  __MyLittleDog_H

#include <littledog.h>

#include <LittleDogUI.h>


/****************************************************************************/
class MyLittleDog : public LittleDog
{
public:

	MyLittleDog();
  
	/*
	 *  LittleDog virtuals.
	 */
	virtual bool initControl( void );
	virtual void updateControl( void );
	virtual void uninitControl( void );
};

#endif // __MyLittleDog_H

