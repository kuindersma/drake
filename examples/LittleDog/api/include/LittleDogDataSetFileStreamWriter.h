
#ifndef __LittleDogDataSetFileStreamWriter_H__
#define __LittleDogDataSetFileStreamWriter_H__

#define BDU_PRIVILEGED_DATA_SET_FILE_STREAM_WRITER_CLASS LittleDogDataSetFileStreamWriter
#include <bduDataSetFileStreamWriter.h>

//!
//!  \addtogroup ld_api  LittleDog Primary API
//!

/****************************************************************************/
//!
//!  \class LittleDogDataSetFileStreamWriter LittleDogDataSetFileStreamWriter.h
//!
//!  \brief a dataset writer customized for LittleDog
//!
//!  \ingroup bdu
//!  @{
//!
class LittleDogDataSetFileStreamWriter : public bduDataSetFileStreamWriter
{
 public:

	//!
	//! \_Description
	//!
	//!    Constructor.  A dataset is created in memory, but only memory
	//!    for one sample is allocated.  A constant sample duration of 
	//!    0.01 seconds is assumed and writen to the output file.  All sample
	//!    values will be stored internally as floats.  As each sample is
	//!    completed, the sample is streamed (flushed) to the file.  Finally,
	//!    when save() is called, the file is polished for reading and closed.
	//!
	//!    See bduDataSetFileStreamWriter for complete documentation.
	//!
	//! \_Parameters
	//!
	//!    \_in   save_as_binary - whether to save in binary format
	//!    \_in   file_name      - file to save to; can be relative or absolute
	//!
	//!    If save_as_binary is true, the data is stored in a more
	//!    compact binary format.
	//!
	LittleDogDataSetFileStreamWriter( bool save_as_binary, const char * file_name);
};

//! @}

#endif // __LittleDogDataSetFileStreamWriter_H__
