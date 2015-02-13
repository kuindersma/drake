
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

#ifndef __bduDataSetFileReader__
#define __bduDataSetFileReader__

//!
//!  \addtogroup bdu   Utility Classes and Functions
//!

/****************************************************************************/
//!
//!  \class bduDataSetFileReader bduDataSetFileReader.h
//!
//!  \brief The bduDataSetFileReader class is a light-weight dataset file reader
//!
//!  \ingroup bdu
//!  @{
//!
class bduDataSetFileReader
{
 public:

	//!
	//! \_Description
	//!
	//!    Destructor.  Internal memory is freed.
	//!
	~bduDataSetFileReader();

	//!
	//!  \_Returns    true if file open, otherwise false
	//!
	bool isFileOpen() const;

	//!
	//! \_Description
	//!
	//!    Returns the number of variables found in the file opened.
	//!    Returns 0 if the file is not opened successfully, or empty.
	//!
	//!  \_Returns    number of variables in file
	//!
	int getNumberOfVariables() const;

	//!
	//! \_Description
	//!
	//!    Returns the number of samples found in the file opened.
	//!    Returns 0 if the file is not opened successfully, or empty.
	//!    Each variable will have this same number of samples.
	//!
	//!  \_Returns    number of data samples in file
	//!
	int getNumberOfSamples() const;

	//!
	//! \_Description
	//!
	//!    Returns the duration of samples, in seconds.
	//!    The samples are presumed to have a fixed sample rate.
	//!    Returns 0 if the file is not opened successfully, or empty.
	//!
	//!  \_Returns    duration of samples, in seconds
	//!
	float getSampleDurationInSeconds() const;

	//!
	//! \_Description
	//!
	//!    Returns the samples for a particular variable.
	//!
	//!    The variable is identified with an index from 0 to the number
	//!    of variables - 1. The returned value is a pointer to a float
	//!    array of sample count length.  Returns NULL if the file was
	//!    not opened successfully, or empty, or if the var_index is out
	//!    of range.
	//!
	//! \_Parameters
	//!
	//!    \_in   var_index - index of desired variable
	//!
	//!  \_Returns    pointer to data of sample at given index
	//!	
	const float* getSamplesForVariable(int var_index) const;

	//!
	//! \_Description
	//!
	//!    This function returns the name of the variable at the given
	//!    index.
	//!
	//! \_Parameters
	//!
	//!    \_in   var_index - index of desired variable
	//!
	//!  \_Returns    name of variable at given index
	//!	
	const char* getVariableName(int var_index) const;

 private:

	//
	//  Private internal objects.
	//

	friend class BDU_PRIVILEGED_DATA_SET_FILE_READER_CLASS;

	bduDataSetFileReader(const char* file_name);
	bduDataSetFileReader(const bduDataSetFileReader &);
	bduDataSetFileReader & operator = (const bduDataSetFileReader &);

	struct        DataSetPrivate* p;
};

//! @}

#endif // __bduDataSetFileReader__

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
