
#ifndef __LittleDogDataSetFileWriter_H__
#define __LittleDogDataSetFileWriter_H__

#define BDU_PRIVILEGED_DATA_SET_FILE_WRITER_CLASS LittleDogDataSetFileWriter

#include <bduDataSetFileWriter.h>

//!
//!  \addtogroup ld_api  LittleDog Primary API
//!

/****************************************************************************/
//!
//!  \class LittleDogDataSetFileWriter LittleDogDataSetFileWriter.h
//!
//!  \brief a dataset writer customized for LittleDog
//!
//!  \ingroup bdu
//!  @{
//!
class LittleDogDataSetFileWriter : public bduDataSetFileWriter
{
 public:

	//!
	//! \_Description
	//!
	//!    Constructor.  A dataset is created in memory, pre-allocating the
	//!    the requested number of samples.  A constant sample duration of 
	//!    0.01 seconds is assumed and writen to the output file.  All sample
	//!    values will be stored internally as floats.
	//!
	//!    See bduDataSetFileWriter for complete documentation.
	//!
	//! \_Parameters
	//!
	//!    \_in   allocated_sample_count - the number of samples (arrays of
	//!           floats) allocated for each variable added (see addVariable()).
	//!           If less than this number of samples is logged (see setData()),
	//!           then some of the allocated memory is just never filled in and
	//!           is not written out.
	//!
	LittleDogDataSetFileWriter( int allocated_sample_count );
};

//! @}

#endif // __LittleDogDataSetFileWriter_H__
