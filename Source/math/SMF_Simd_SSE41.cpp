/*
  SMF -  Super Math Fabric  
  Copyright (C) 2014 Rasputtim <Rasputtim@hotmail.com>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the freeMem Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the freeMem Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

  */

//#pragma hdrstop

#include "math/SMF_SimdGeneric.h"
#include "math/SMF_SimdMMX.h"
#include "math/SMF_SimdSSE.h"
#include "math/SMF_SimdSSE2.h"
#include "math/SMF_SimdSSE3.h"
#include "math/SMF_SimdSSE41.h"
#include "math/SMF_SimdSSE42.h"
#include "math/SMF_Matriz.h"
#include "math/SMF_Vector.h"
#include "geometry/SMF_DrawVert.h"
#include "math/SMF_JointTransform.h"

namespace SMF{
namespace MATH{
//===============================================================
//
//	SSE3 implementation of CSIMDProcessor
//
//===============================================================

#if defined(MACOS_X) && defined(__i386__)

/*
============
CSIMD_SSE3::getName
============
*/
const char * CSIMD_SSE41::getName() const {
	return "MMX & SSE & SSE2 & SSE3 & SSE4.1";
}

#elif defined(WIN32)


/*
============
CSIMD_SSE3::getName
============
*/
const char * CSIMD_SSE41::getName() const {
	return "MMX & SSE & SSE2 & SSE3 & SSE4.1";
}


#endif /* _WIN32 */

} //end MATH
}//end SMF
