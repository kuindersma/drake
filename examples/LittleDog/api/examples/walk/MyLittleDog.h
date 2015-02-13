
#ifndef __MYLITTLEDOG_H__
#define __MYLITTLEDOG_H__

#include <littledog.h>

class MyLittleDog : public LittleDog
{
	friend class KinematicWalker;
  
 public:

	MyLittleDog();
	~MyLittleDog();
	
	virtual bool initControl( void );
	virtual void updateControl( void );
	virtual void uninitControl( void );

 private:
	
	class KinematicWalker * m_walker;
	bool m_walker_initialized;
};

#endif
