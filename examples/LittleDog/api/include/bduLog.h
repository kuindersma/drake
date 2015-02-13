
/*
 *  Copyright (C) 2006-2008 Boston Dynamics
 *  ALL RIGHTS RESERVED.
 *
 *  These coded instructions, statements, and computer programs
 *  contain unpublished proprietary information of Boston Dynamics
 *  and are protected by Copyright Laws of the United States.
 *  They may not be used, duplicated, or disclosed in any form, in
 *  whole or in part, without the prior written consent from Boston
 *  Dynamics.
 *
 *  RESTRICTED RIGHTS LEGEND
 *  Use, duplication, or disclosure by the government is subject
 *  to restrictions as set forth in FAR 52.227.19(c)(2) or
 *  subparagraph (c)(1)(ii) of the Rights in Technical Data and
 *  Computer Software clause at DFARS 252.227-7013 and/or in
 *  similar or successor clauses in the FAR, or the DOD or NASA
 *  FAR Supplement, or to subparagraphs (c)(1) and (c)(2) of the
 *  Commercial Computer Software--Restricted Rights at 48 CFR
 *  52.227-19, as applicable.  Unpublished-rights reserved under
 *  the Copyright Laws of the United States.
 *  Contractor/Manufacturer is:
 *  Boston Dynamics/78 Fourth Avenue/Waltham MA 02451.
 */

#ifndef __bduLog_H__
#define __bduLog_H__

//!
//!  \addtogroup bdu   Utility Classes and Functions
//!

class bdiLog;


typedef int bduLogCallbackFunction(int notify_level, const char* string, void* user_data);


//!
//!  \brief   Enumeration of log types supported by bduLog.
//!
enum
{

	BDU_LOG_TYPE_STDOUT = 0,     //!<  Send output to stdout
	BDU_LOG_TYPE_STDERR,         //!<  Send output to stderr
	BDU_LOG_TYPE_FILE,           //!<  Send output to file
	BDU_LOG_TYPE_CALLBACK        //!<  Send output to provided callback function

};


//!
//!  \brief   Enumeration of log levels supported by bduLog.
//!
enum
{
	BDU_LOG_LEVEL_FATAL  = 1,    //!<  Message indicates a fatal condition
	BDU_LOG_LEVEL_ERROR,         //!<  Message indicates something is seriously wrong
	BDU_LOG_LEVEL_WARN,          //!<  Message indicates of potential concern happened
	BDU_LOG_LEVEL_INFO,          //!<  Message is for information purposed
	BDU_LOG_LEVEL_DEBUG          //!<  Message is for debugging info
};


//!
//!  \brief    Log print function.
//!
//!  \_Description
//!
//!    This function prints the passed string to *all*
//!    open logs.  All logs that print at least the
//!    given log_level will print the string.
//!
//!  \_Parameters
//!
//!    \_in   log_level  - level of message
//!    \_in   string     - message to print
//!
void bdu_log_print(int log_level, const char* string);

//!
//!  \brief    Log printf function.
//!
//!  \_Description
//!
//!    This function prints the passed string to *all*
//!    open logs.  All logs that print at least the
//!    given log_level will print the string.
//!
//!    This function can be used essentially like fprintf().
//!
//!  \_Parameters
//!
//!    \_in   log_level  - level of message
//!    \_in   format     - printf-style format string
//!    \_in   ...        - varargs for format
//!
void bdu_log_printf(int log_level, const char* format, ...);


/****************************************************************************/
//!
//!  \class   bduLog bduLog.h
//!
//!  \brief   The bduLog class manages the Boston Dynamics log system.
//!
//!  \ingroup bdu
//!  @{
//!
class bduLog
{

public:

	//!
	//!  \brief    General constructor.
	//!
	//!  \_Description
	//!
	//!    This function opens a new output log.  The type of
	//!    log and verbosity of output depend on the passed
	//!    parameters.
	//!
	//!    More than one log can be active at any given time.
	//!    e.g., an application can have one log printing messages
	//!    to stdout, and another printing messages to a file.
	//!
	//!  \_Parameters
	//!
	//!    \_in   log_type    - type of log to open
	//!    \_in   log_level   - maximem level of messages to print
	//!
	//!    If BDU_LOG_TYPE_FILE is passed as the log_type, the
	//!    created file will be called "bdu_log.txt".
	//!
	bduLog(int log_type = BDU_LOG_TYPE_STDERR,
		int log_level = BDU_LOG_LEVEL_WARN);

	//!
	//!  \brief    File log constructor.
	//!
	//!  \_Description
	//!
	//!    This function opens a new file output log.
	//!
	//!    More than one log can be active at any given time.
	//!    e.g., an application can have one log printing messages
	//!    to stdout, and another printing messages to a file.
	//!
	//!  \_Parameters
	//!
	//!    \_in   log_level  - maximem level of messages to print
	//!    \_in   filename   - name of file to print messages into
	//!    \_in   clear_file - whether to clear file (1), or append to it (0)
	//!
	bduLog(int  log_level,
		const char*  filename,
		bool         clear_file);

	//!
	//!  \brief    Callback log constructor.
	//!
	//!  \_Description
	//!
	//!    This function opens a log that prints messages to
	//!    to the provided callback.
	//!
	//!    The passes user_data pointer will not be directly used
	//!    by the log.  It will simply be passed back in the callback
	//!    function.  One common use of the user_data pointer is to
	//!    have it hold a pointer to UI widget that the callback
	//!    function can then print to.
	//!
	//!  \_Parameters
	//!
	//!    \_in   log_level          - maximem level of messages to print
	//!    \_in   callback           - function to be called for log output
	//!    \_in   callback_user_data - generic pointer that will be returned in callback
	//!
	bduLog(int              log_level,
		bduLogCallbackFunction*  callback,
		void*                    callback_user_data);

	//!
	//!  \brief    Destructor.
	//!
	//!  \_Description
	//!
	//!    This function closes the log.
	//!
	~bduLog();

private:

	//
	//  Private, internal objects.
	//

	int   m_log_type;
	int   m_log_level;
	int   m_log_handle;

};

//! @}


//
//  Debugging macros.
//
#ifndef NDEBUG

#define bdu_log_assert(assertion) if (!(assertion)) {bdu_log_printf(BDU_LOG_LEVEL_ERROR, "Assertion '%s' failed in file %s line %d.\n", #assertion, __FILE__, __LINE__);}
#define bdu_log_debug(string) bdu_log_print(BDU_LOG_LEVEL_DEBUG, string)

#else

#define bdu_log_assert(assertion) /* */
#define bdu_log_debug(string) /* */

#endif


#endif //__bduLog_H__

/*
 *  Copyright (C) 2006-2008 Boston Dynamics
 *  ALL RIGHTS RESERVED.
 *
 *  These coded instructions, statements, and computer programs
 *  contain unpublished proprietary information of Boston Dynamics
 *  and are protected by Copyright Laws of the United States.
 *  They may not be used, duplicated, or disclosed in any form, in
 *  whole or in part, without the prior written consent from Boston
 *  Dynamics.
 *
 *  RESTRICTED RIGHTS LEGEND
 *  Use, duplication, or disclosure by the government is subject
 *  to restrictions as set forth in FAR 52.227.19(c)(2) or
 *  subparagraph (c)(1)(ii) of the Rights in Technical Data and
 *  Computer Software clause at DFARS 252.227-7013 and/or in
 *  similar or successor clauses in the FAR, or the DOD or NASA
 *  FAR Supplement, or to subparagraphs (c)(1) and (c)(2) of the
 *  Commercial Computer Software--Restricted Rights at 48 CFR
 *  52.227-19, as applicable.  Unpublished-rights reserved under
 *  the Copyright Laws of the United States.
 *  Contractor/Manufacturer is:
 *  Boston Dynamics/78 Fourth Avenue/Waltham MA 02451.
 */
