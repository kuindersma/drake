
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

#ifndef __bduMat4f_H__
#define __bduMat4f_H__

//!
//!  \addtogroup bdu   Utility Classes and Functions
//!

/****************************************************************************/
//!
//!   \class bduMat4f bduMat4f.h
//!
//!   \brief    a light-weight 3x3 matrix class
//!
//!  \ingroup bdu
//!  @{
//!
class bduMat4f
{

public:

	enum InitialValues {
		INITIALIZE_TO_ZERO = 0,  //!<  Initialize all values to 0
		INITIALIZE_TO_IDENTITY   //!<  Initialize to the identity matrix
	};

	//!
	//! \_Description
	//!
	//!    Default constructor.  Matrix values are uninitialized.
	//!
	bduMat4f() {}

	//!
	//! \_Description
	//!
	//!    Constructor.  Matrix values are all initialized based
	//!    on initial_values parameter, which should be ZERO or
	//!    IDENTITY.
	//!
	//! \_Parameters
	//!
	//!    \_in   initial_values - how to initialize matrix values
	//!
	bduMat4f(InitialValues initial_values);

	//!
	//! \_Description
	//!
	//!    Current matrix values.  Can be accessed directly.
	//!
	float m[4][4];

};
//! @}

#endif //__bduMat4f_H__

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
