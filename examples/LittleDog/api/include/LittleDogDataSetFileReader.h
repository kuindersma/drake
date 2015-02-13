
#ifndef __LittleDogDataSetFileReader_H__
#define __LittleDogDataSetFileReader_H__

#define BDU_PRIVILEGED_DATA_SET_FILE_READER_CLASS LittleDogDataSetFileReader
#include <bduDataSetFileReader.h>

//!
//!  \addtogroup ld_api  LittleDog Primary API
//!

/****************************************************************************/
//!
//!  \class LittleDogDataSetFileReader LittleDogDataSetFileReader.h
//!
//!  \brief a dataset reader customized for LittleDog
//!
//!  \ingroup bdu
//!  @{
//!
class LittleDogDataSetFileReader : public bduDataSetFileReader
{
 public:

	//!
	//! \_Description
	//!
	//!    See bduDataSetFileReader for documentation.
	//!
	//!    This class is a LittleDog-specific dataset reader.
	//!
	LittleDogDataSetFileReader(const char* file_name);
};

//! @}

#endif // __LittleDogDataSetFileReader_H__
