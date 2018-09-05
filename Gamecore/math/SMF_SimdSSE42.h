/*
  SMF - Super Math Fabric  (https://sourceforge.net/projects/sgfabric/)
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

#ifndef _SMF__MATH_SIMD_SSE42_H__
#define _SMF__MATH_SIMD_SSE42_H__

#include "../SMF_Config.h"
#include "SMF_SimdSSE41.h"
namespace SMF {
namespace MATH{
/*
===============================================================================

	SSE42 implementation of CSIMDProcessor
	
===============================================================================
*/
/**
 * \class CSIMD_SSE42
 *
 * \ingroup SMF_Math
 *
 * \brief Implementação SSE4.2 do processador SIMD
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 
 * \note 
 **/
class SMF_API CSIMD_SSE42 : public CSIMD_SSE41 {
public:
	CSIMD_SSE42():CSIMD_SSE41(CPUID_SSE41){};
	CSIMD_SSE42(cpuid_t id):CSIMD_SSE41(id){};
	~CSIMD_SSE42(){};
#if defined(MACOS_X) && defined(__i386__)
	virtual const char * VPCALL getName() const;

#elif defined(_WIN32)
	virtual const char * VPCALL getName() const;

#endif
};

} //end MATH
} //end SMF
#endif /* !__MATH_SIMD_SSE3_H__ */
