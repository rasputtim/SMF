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

#ifndef _SMF__MATH_SIMD_MMX_H__
#define _SMF__MATH_SIMD_MMX_H__
#include "../SMF_Config.h"
#include "SMF_SimdGeneric.h"
namespace SMF {
	namespace MATH{
/*
===============================================================================

	MMX implementation of CSIMDProcessor

===============================================================================
*/
/**
 * \class CSIMD_MMX
 *
 * \ingroup SMF_Math
 *
 * \brief Implementação do processador SIMD MMX 
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 
 * \note http://pt.wikipedia.org/wiki/MMX
 **/
class SMF_API CSIMD_MMX : public CSIMD_Generic {
public:
	CSIMD_MMX():CSIMD_Generic(CPUID_MMX){};
	CSIMD_MMX(cpuid_t id):CSIMD_Generic(id){};
	~CSIMD_MMX(){};

#if defined(MACOS_X) && defined(__i686__)
	virtual const char * VPCALL getName() const;

#else
	virtual const char * VPCALL getName() const;

	virtual void VPCALL memCopy( void *dst,			const void *src,		const int count );
	virtual void VPCALL memSet( void *dst,			const int val,			const int count );
///\todo implement this methods
	virtual void VPCALL quat_to_mat3x4(sf_m128 q, sf_m128 t, sf_m128 *m){};
	virtual sf_m128 VPCALL quat_transform_vec4(sf_m128 quat, sf_m128 vec){sf_m128 teste; return teste;};
	virtual sf_m128 VPCALL quat_mul_quat(sf_m128 q1, sf_m128 q2){sf_m128 teste; return teste;};
	virtual sf_m128 VPCALL quat_div_quat(sf_m128 q1, sf_m128 q2){sf_m128 teste; return teste;};

#endif
};

} //end MATH
} //end SMF
#endif /* !__MATH_SIMD_MMX_H__ */
