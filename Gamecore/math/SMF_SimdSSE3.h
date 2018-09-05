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

#ifndef _SMF__MATH_SIMD_SSE3_H__
#define _SMF__MATH_SIMD_SSE3_H__

#include "../SMF_Config.h"
#include "SMF_SimdSSE2.h"
namespace SMF {
namespace MATH{
/*
===============================================================================

	SSE3 implementation of CSIMDProcessor
	http://pt.wikipedia.org/wiki/SSE3
===============================================================================
*/
/**
 * \class CSIMD_SSE3
 *
 * \ingroup SMF_Math
 *
 * \brief Implementação SSE3 do processador SIMD
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 
 * \note http://pt.wikipedia.org/wiki/SSE3
 **/
class SMF_API CSIMD_SSE3 : public CSIMD_SSE2 {
public:
	CSIMD_SSE3():CSIMD_SSE2(CPUID_SSE2){};
	CSIMD_SSE3(cpuid_t id):CSIMD_SSE2(id){};
	~CSIMD_SSE3(){};
#if defined(MACOS_X) && defined(__i386__)
	virtual const char * VPCALL getName() const;

#elif defined(_WIN32)
	virtual const char * VPCALL getName() const;
	virtual void   VPCALL transformVerts( CVertex *verts, const int numVerts, const CMatJoint3x4 *joints, const CVec4D *weights, const int *index, const int numWeights );
	virtual float  VPCALL CSIMD_SSE3::vector3D_Dot(const CVec3D* pSrc1, const CVec3D* pSrc2);
	virtual float  VPCALL CSIMD_SSE3::vector4D_Dot(const CVec4D* pSrc1, const CVec4D* pSrc2);
	virtual float  VPCALL CSIMD_SSE3::vector3D_LengthSq(const CVec3D* pSrc1);
	virtual float  VPCALL CSIMD_SSE3::vector3D_Length(const CVec3D* pSrc1);
	virtual void   VPCALL CSIMD_SSE3::vector3D_Normalize(CVec3D* pVec);
	virtual void   VPCALL CSIMD_SSE3::vector3D_NormalizeOf(CVec3D* pOut, const CVec3D* pVec);
	virtual float  VPCALL CSIMD_SSE3::vector3D_Distance(const CVec3D* pVec1, const CVec3D* pVec2);

	virtual sf_m128  VPCALL mat4x4_mul_sse(const sf_m128 *matrix, sf_m128 vector);

#endif


};

} //end MATH
} //end SMF
#endif /* !__MATH_SIMD_SSE3_H__ */
