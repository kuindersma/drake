
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

#ifndef __bduDataSetFileWriter__
#define __bduDataSetFileWriter__

//!
//!  \addtogroup bdu   Utility Classes and Functions
//!

typedef void* bduDataSetVarID;


/****************************************************************************/
//!
//!  \class bduDataSetFileWriter bduDataSetFileWriter.h
//!
//!  \brief This class is an in-memory, light-weight dataset file writer.  
//!
//!  Memory for all of the samples is allocated up-front before
//!  the first call to setData(), which populates the data.  The file 
//!  is written when save() is is called.  All memory is freed when the 
//!  destructor is called.
//!
//!  \ingroup bdu
//!  @{
//!

class bduDataSetFileWriter
{
 public:

	//!
	//! \_Description
	//!
	//!    Destructor.   Internal memory is freed.  save() should be
	//!    called before deleting the object or data will not be saved.
	//!
	~bduDataSetFileWriter();

	//!
	//! \_Description
	//!
	//!    Saves a dataset to a file.  The data is written to a file
	//!    with the given name.
	//!
	//!  \_Parameters
	//!
	//!    \_in   file_name      - file to save to; can be relative or absolute
	//!    \_in   sample_count   - number of samples to save
	//!    \_in   save_as_binary - whether to save in binary format
	//!
	//!    Only sample_count samples are written.
	//!
	//!    sample_count must be less than  the amount allocated in the
	//!    constructor.
	//!
	//!    Sample data is originally cleared to zero, so use sample_count to
	//!    indicate how much of the data (see setData()) is saved.
	//!
	//!    If save_as_binary is true, the data is stored in a more
	//!    compact binary format.
	//!
	//!  \_Returns    true if save succeeded, otherwise false
	//!
	bool save(const char* file_name,
		int sample_count,
		bool save_as_binary = true) const;
	
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
	//!    Sets the value for a particular variable in a particular sample.
	//!    The var_id identifies the variable to change.  The sample_num
	//!    identifies which sample the sample_value should be stored in.
	//!
	//!    Sample data is allocated in the constructor and is initialized
	//!    to zero.
	//!
	//!  \_Parameters
	//!
	//!    \_in   var_id       - ID of variable as returned by addVariable()
	//!    \_in   sample_num   - number of samples to log
	//!    \_in   sample_value - the value that the variable should take
	//!
	//!    var_id must be a valid ID, as returned by addVariable().
	//!
	//!    sample_num should fall between 0 and (sample_count -1), where 
	//!    sample_count is the allocated_sample_count passed into the
	//!    constructor.
	//!
	//!  \_Returns    true on success, otherwise false
	//!
	bool setData(bduDataSetVarID var_id, int sample_num, float sample_value);

 private:

	//
	//  Private internal objects.
	//
	friend class BDU_PRIVILEGED_DATA_SET_FILE_WRITER_CLASS;

	bduDataSetFileWriter(float sample_duration_in_seconds, int allocated_sample_count);
	bduDataSetFileWriter(const bduDataSetFileWriter &);
	bduDataSetFileWriter& operator = (const bduDataSetFileWriter &);

	struct DataSetPrivate* p;
};

//! @}

#endif // __bduDataSetFileWriter__

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
