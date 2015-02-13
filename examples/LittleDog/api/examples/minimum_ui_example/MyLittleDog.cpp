
//
//  System includes.
//
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

//
//  LittleDog API includes.
//
#include <littledog.h>
#include <bduLog.h>

//
//  Local includes.
//
#include "MyLittleDog.h"


/****************************************************************************/
MyLittleDog::MyLittleDog()
{
}


/****************************************************************************/
//
//  This utility function checks returned error codes.
//
void
handle_result(LittleDog* dog,
	const char* func,
	LD_ERROR_CODE result,
	bool exit_on_error)
{
	// if an error occurs, print it out and (optionally) exit

	if (result == LD_OKAY)
		return;
	
	bdu_log_printf(BDU_LOG_LEVEL_WARN,
		"An error occurred calling function '%s()'.\n"
		"    Error was: '%s'.\n", 
		func,
		dog->getErrorCodeString(result));
	
	if (exit_on_error)
		exit(EXIT_FAILURE);
}


/****************************************************************************/
bool
MyLittleDog::initControl(void)
{
	bdu_log_printf(BDU_LOG_LEVEL_WARN, "initControl() called.\n");

	//
	// if we do all of our planning in initControl() we can call this here.
	//
	donePlanning(); // this switches us immediately to STATE_TRIAL_RUNNING

	//
	//  Return true to say control init was successful.
	//
	return true;
}


/****************************************************************************/
void
MyLittleDog::updateControl(void)
{
	//
	//  This minimum example does nothing.
	//
}


/****************************************************************************/
void
MyLittleDog::uninitControl(void)
{
	bdu_log_printf(BDU_LOG_LEVEL_WARN, "uninitControl() called.\n");
}
