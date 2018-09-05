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

#ifndef _SMF__MATH_CPUID_SSE2_H__
#define _SMF__MATH_CPUID_SSE2_H__

#include "../SMF_Config.h"
#include "SMF_SimdSSE.h"
namespace SMF {
namespace GEO{
class CTriangleMesh;
}
	namespace MATH{
/*
===============================================================================

	SSE2 implementation of CSIMDProcessor

===============================================================================
*/
/**
 * \class CSIMD_SSE2
 *
 * \ingroup SMF_Math
 *
 * \brief Implementação SSE2 do processador SIMD
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 
 * \note http://pt.wikipedia.org/wiki/SSE2
 **/
class SMF_API CSIMD_SSE2 : public CSIMD_SSE {
public:

	CSIMD_SSE2():CSIMD_SSE(CPUID_SSE2){};
	CSIMD_SSE2(cpuid_t id):CSIMD_SSE(id){};
	~CSIMD_SSE2(){};

#if defined(MACOS_X) && defined(__i386__)
	virtual const char * VPCALL getName() const;
	virtual void VPCALL cmpLT( sf_u8 *dst,			const sf_u8 bitNum,		const float *src0,		const float constant,	const int count );

#else
	virtual const char * VPCALL getName() const;
#if 0
	virtual void VPCALL matX_LowerTriangularsolve( const CMatXD &L, float *x, const float *b, const int n, int skip = 0 );
	virtual void VPCALL matX_LowerTriangularsolveTranspose( const CMatXD &L, float *x, const float *b, const int n );
#endif
	virtual void VPCALL  mixedSoundToSamples( short *samples, const float *mixBuffer, const int numSamples );

#endif
	virtual sf_m128  VPCALL colmajor_mat4x4_mul_sse(const sf_m128 *matrix, sf_m128 vector);
	virtual  void   VPCALL mat4x4_mul_sse(sf_m128 *out, const sf_m128 *m1, const sf_m128 *m2);
static float intersectRay(const CTriangleMesh &trig, const CRay &ray);
static float intersectRay_TriangleIndex(const CTriangleMesh &trig, const CRay &ray, int &outTriangleIndex);
static float intersectRay_TriangleIndex_UV(const CTriangleMesh &trig, const CRay &ray, int &outTriangleIndex, float &outU, float &outV);

};

} //end MATH
} //end SMF
#endif /* !__MATH_CPUID_SSE2_H__ */
