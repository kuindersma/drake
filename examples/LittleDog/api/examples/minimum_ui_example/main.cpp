
//
//  System includes.
//
#include <stdlib.h>
#include <signal.h>

//
//  Qt includes.
//
#include <qapplication.h>

//
//  LittleDog API includes.
//
#include <littledog.h>
#include <bduLog.h>

//
//  Local includes.
//
#include "MyLittleDog.h"
#include "LittleDogUI.h"

//
//  Forward declarations.
//
static void signal_handler( int );

//===================================================================

static MyLittleDog * dog = NULL;
static QApplication * qt_app = NULL;

//===================================================================
static void signal_handler( int )
{
	static bool busy = false;
	if ( busy )
		return;
	busy = true;

	printf( "Stopping control process\n" );
	fflush( stdout );
  
	// We cannot delete 'dog' right here because it may be in use
	// in, say, the run or calibrate process.   By aborting the 
	// connection, it assures that these processes end quickly.

	LD_ERROR_CODE result = dog->abortConnection();

	if ( result != LD_OKAY )
		printf( "An error occurred calling 'abortConnection()'."
			"\n\tError was: '%s'\n", dog->getErrorCodeString(result));

	fprintf( stdout, "Waiting for state transition..." );

	while ( dog->getState() != LittleDog::STATE_UNINITIALIZED) 
	{
		sleep(1);
		fprintf(stdout, ".");
		fflush(stdout);
	}
	
	fprintf( stdout, "done!\n");
	qt_app->quit();
}

/****************************************************************************/
int 
main(int argc, char* argv[]) 
{
	signal( SIGINT, signal_handler );

	//
	// Check for max plan duration arg.  atoi returns 0 if it fails, not -1
	//
	int max_plan_time = -1;

	if ( argc > 1 )
		max_plan_time = atoi(argv[1]);

	//
	//  Set up a log that prints to a file.
	//
	bduLog* log_file = new bduLog(BDU_LOG_LEVEL_WARN,
		"littledog_log.txt",
		true);

	//
	//  Create the MyLittleDog object.
	//
	dog = new MyLittleDog;

	//
	//  Run the GUI.
	//
	QApplication app(argc, argv);
	qt_app = &app;

	LittleDogUI gui(*dog,max_plan_time);
	gui.show();

	app.connect(&app, SIGNAL(lastWindowClosed()), &app, SLOT(quit()));
	int exec_retval = app.exec();

	//
	//  Destroy the MyLittleDog object.
	//
	delete dog;

	//
	//  Close the log file.
	//
	delete log_file;

	return exec_retval;
}
