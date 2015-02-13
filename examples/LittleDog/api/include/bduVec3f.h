
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

#ifndef __bduVec3f_H__
#define __bduVec3f_H__

//!
//!  \addtogroup bdu   Utility Classes and Functions
//!

/****************************************************************************/
//!
//!  \class bduVec3f bduVec3f.h
//!
//!  \brief The bduVec3f class is a light-weight 3 vector.
//!
//!  \ingroup bdu
//!  @{
//!
class bduVec3f
{

public:

	//!
	//! \_Description
	//!
	//!    Default constructor.  Vector values are uninitialized.
	//!
	bduVec3f() {}

	//!
	//! \_Description
	//!
	//!    Constructor.  Vector values are all initialized to the
	//!    passed value.
	//!
	//! \_Parameters
	//!
	//!    \_in   initial_value - initial value of all vector values
	//!
	bduVec3f(float initial_value);

	//!
	//! \_Description
	//!
	//!    Constructor.  Vector values are initialized set to passed
	//!    values.
	//!
	//! \_Parameters
	//!
	//!    \_in   x  - 1st vector value
	//!    \_in   y  - 2nd vector value
	//!    \_in   z  - 3rd vector value
	//!
	bduVec3f(float x, float y, float z);

	//!
	//! \_Description
	//!
	//!    Current vector values.  Can be accessed directly.
	//!
	float n[3];

};
//! @}


#endif //__bduVec3f_H__

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
