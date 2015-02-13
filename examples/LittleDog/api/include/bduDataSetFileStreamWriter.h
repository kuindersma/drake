
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

#ifndef __bduDataSetFileStreamWriter__
#define __bduDataSetFileStreamWriter__

//!
//!  \addtogroup bdu   Utility Classes and Functions
//!

typedef void* bduDataSetVarID;


/****************************************************************************/
//!
//!  \class bduDataSetFileStreamWriter bduDataSetFileStreamWriter.h
//!
//!  \brief This class is a low-memory, light-weight dataset file writer.
//!
//!  This class allocates little memory in its lifetime and writes each 
//!  sample out "to file" as the data is sampled.  The file is closed when
//!  save() is called.  All memory is freed when the destructor is called.
//!
//!  \ingroup bdu
//!  @{
//!

class bduDataSetFileStreamWriter
{
 public:

	//!
	//! \_Description
	//!
	//!    Destructor.   Internal memory is freed.   save() should be
	//!    called before deleting the object or the file will be
	//!    deleted by the destructor.
	//!
	~bduDataSetFileStreamWriter();

	//!
	//! \_Description
	//!
	//!    Finishes saving a dataset to the file.  The file is "polished"
	//!    for reading and then closed.  If save_as_binary was true in
	//!    the constructor data was stored in a more compact binary format.
	//!
	//!  \_Returns    true if save succeeded, otherwise false
	//!
	bool save( void );

	//!
	//! \_Description
	//!
	//!    Writes the sample values (set with calls to setData()) to the
	//!    file specified in the constructor.  The sample's memory is
	//!    then zeroed in preparation for the next set of calls to
	//!    setData().  Once all samples have been saved, save() can be
	//!    called to close and "polish" the output file for reading.
	//!
	//!  \_Returns    true on success, otherwise false
	//!
	bool saveSample( void );


	//!
	//! \_Description
	//!
	//!    Adds a new variable to the dataset.  All variables should be added
	//!    before setting the variable data.  The ID returned should be
	//!    stored for future use (or see getVariableID()).  The ID is used
	//!    as input to the setData() function.  Each variable has a unique ID.
	//!    Each variable should also have a unique name.  It is suggested
	//!    that all variable names be in lower case and have no spaces or
	//!    unusual characters in their names.
	//!  
	//!  \_Parameters
	//!
	//!    \_in   var_name - the name of the added variable
	//!
	//!  \_Returns    ID of added variable
	//!
	bduDataSetVarID  addVariable(const char* var_name);

	//!
	//! \_Description
	//!
	//!    Returns the bduDataSetVarID associated with a variable name.
	//!    If the variable is not found, this function returns 0.
	//!
	//!  \_Parameters
	//!
	//!    \_in   var_name - name used in addVariable() call
	//!
	//!  \_Returns    ID of variable with given name
	//!
	bduDataSetVarID  getVariableID(const char* var_name) const;

	//!
	//! \_Description
	//!
	//!    Returns the variable name assoicated with a bduDataSetVarID.
	//!    If the variable is not found, this function returns NULL.
	//!
	//!  \_Parameters
	//!
	//!    \_in   var_id - ID of added variable
	//!
	//!  \_Returns    name of variable with given ID
	//!
	const char* getVariableName(bduDataSetVarID var_id) const;
	
	//!
	//! \_Description
	//!
	//!    Sets the value for a particular variable for the current sample.
	//!    This function should not be called until all variables are added.
	//!    The var_id identifies the variable to change.  After all variables
	//!    have been set with setData() for a sample, then saveSample()
	//!    should be called.  This process repeats for each variable and sample.
	//!
	//!    Sample data is zeroed after each sample is written with saveSample().
	//!
	//!  \_Parameters
	//!
	//!    \_in   var_id       - ID of variable as returned by addVariable()
	//!    \_in   sample_value - the value that the variable should take
	//!
	//!    var_id must be a valid ID, as returned by addVariable().
	//!
	//!  \_Returns    true on success, otherwise false
	//!
	bool setData(bduDataSetVarID var_id, float sample_value);


 private:

	//
	//  Private internal objects.
	//
	friend class BDU_PRIVILEGED_DATA_SET_FILE_STREAM_WRITER_CLASS;

	bduDataSetFileStreamWriter(float sample_duration_in_seconds, bool save_as_binary, 
		const char * filename);
	bduDataSetFileStreamWriter(const bduDataSetFileStreamWriter &);
	bduDataSetFileStreamWriter& operator = (const bduDataSetFileStreamWriter &);

	struct DataSetPrivate* p;
};

//! @}

#endif // __bduDataSetFileStreamWriter__

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
