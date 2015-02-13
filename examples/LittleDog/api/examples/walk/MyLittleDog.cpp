
#include <memory.h>
#include <string.h>
#include <stddef.h>
#include <stdio.h>

#include <bduLog.h>
#include "MyLittleDog.h"
#include "KinematicWalker.h"

MyLittleDog::MyLittleDog() :
	m_walker(NULL),
	m_walker_initialized(false)
{
}

MyLittleDog::~MyLittleDog()
{
}


bool
MyLittleDog::initControl( void )
{
	bdu_log_print(BDU_LOG_LEVEL_WARN, "Initializing control.\n");

	m_walker = new KinematicWalker( *this );
	m_walker_initialized = false;
	donePlanning();
	return true;
}


void
MyLittleDog::updateControl( void )
{
	if ( m_walker_initialized == false )
	{
		m_walker->initialize();
		m_walker_initialized = true;
	}

	m_walker->update();
}


void
MyLittleDog::uninitControl( void )
{
	bdu_log_print(BDU_LOG_LEVEL_WARN, "Shutting down control.\n");

	m_walker->uninitialize();
	delete m_walker;
	m_walker = NULL;
}


