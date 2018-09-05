/*
  SMF -  Salvathor Math Fabric  (http://smfabric.sourceforge.net)
  Copyright (C) 2014 Salvatore Giannotta Filho <a_materasu@hotmail.com>

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


#include "../../Gamecore/SMF_Config.h"
#include "../../Gamecore/exceptions/all.h"

/** @file TriangleMesh_IntersectRay_SSE.inl
	@author Jukka Jylänki
	@brief SSE implementation of ray-mesh intersection routines. */
//namespace SMF{
//namespace MATH{

#if defined(MATH_GEN_SSE2) && !defined(MATH_GEN_TRIANGLEINDEX)
float CSIMD_SSE2::intersectRay(const CTriangleMesh &trig, const CRay &ray)
#elif defined(MATH_GEN_SSE2) &&defined(MATH_GEN_TRIANGLEINDEX) && !defined(MATH_GEN_UV)
float CSIMD_SSE2::intersectRay_TriangleIndex(const CTriangleMesh &trig, const CRay &ray, int &outTriangleIndex)
#elif defined(MATH_GEN_SSE2) && defined(MATH_GEN_TRIANGLEINDEX) && defined(MATH_GEN_UV)
float CSIMD_SSE2::intersectRay_TriangleIndex_UV(const CTriangleMesh &trig, const CRay &ray, int &outTriangleIndex, float &outU, float &outV)
#elif defined (MATH_GEN_SSE41) && !defined(MATH_GEN_TRIANGLEINDEX)
float CSIMD_SSE41::intersectRay(const CTriangleMesh &trig, const CRay &ray)
#elif defined (MATH_GEN_SSE41) && defined(MATH_GEN_TRIANGLEINDEX) && !defined(MATH_GEN_UV)
float CSIMD_SSE41::intersectRay_TriangleIndex(const CTriangleMesh &trig, const CRay &ray, int &outTriangleIndex)
#elif defined (MATH_GEN_SSE41) && defined(MATH_GEN_TRIANGLEINDEX) && defined(MATH_GEN_UV)
float CSIMD_SSE41::intersectRay_TriangleIndex_UV(const CTriangleMesh &trig, const CRay &ray, int &outTriangleIndex, float &outU, float &outV)
#endif
{
//	std::cout << numTris << " tris: ";
//	TRACESTART(RayTriMeshIntersectSSE);

	SMF_ASSERT(sizeof(CVec3D) == 3*sizeof(float));
	SMF_ASSERT(sizeof(CTriangle) == 3*sizeof(CVec3D));
#ifdef _DEBUG
	SMF_ASSERT(trig.getVertexLayout() == 1); // Must be SoA4 structured!
#endif

	const float inf = CMath::INFINITY_FLOAT;
	sf_m128 nearestD = _mm_set1_ps(inf);
#ifdef MATH_GEN_UV
	sf_m128 nearestU = _mm_set1_ps(inf);
	sf_m128 nearestV = _mm_set1_ps(inf);
#endif
#ifdef MATH_GEN_TRIANGLEINDEX
	__m128i nearestIndex = _mm_set1_epi32(-1);
#endif

	const sf_m128 lX = _mm_load1_ps(&ray.pos.x);
	const sf_m128 lY = _mm_load1_ps(&ray.pos.y);
	const sf_m128 lZ = _mm_load1_ps(&ray.pos.z);

	const sf_m128 dX = _mm_load1_ps(&ray.dir.x);
	const sf_m128 dY = _mm_load1_ps(&ray.dir.y);
	const sf_m128 dZ = _mm_load1_ps(&ray.dir.z);

	const sf_m128 epsilon = _mm_set1_ps(1e-4f);
	const sf_m128 zero = _mm_setzero_ps();
	const sf_m128 one = _mm_set1_ps(1.f);

    const sf_m128 sign_mask = _mm_set1_ps(-0.f); // -0.f = 1 << 31

	SMF_ASSERT(((uintptr_t)trig.toFloatPtr() & 0xF) == 0);

	const float *tris = reinterpret_cast<const float*>(trig.toFloatPtr());

	for(int i = 0; i+4 <= trig.getCount(); i += 4)
	{
		sf_m128 v0x = _mm_load_ps(tris);
		sf_m128 v0y = _mm_load_ps(tris+4);
		sf_m128 v0z = _mm_load_ps(tris+8);

#ifdef SOA_HAS_EDGES
		sf_m128 e1x = _mm_load_ps(tris+12);
		sf_m128 e1y = _mm_load_ps(tris+16);
		sf_m128 e1z = _mm_load_ps(tris+20);

		sf_m128 e2x = _mm_load_ps(tris+24);
		sf_m128 e2y = _mm_load_ps(tris+28);
		sf_m128 e2z = _mm_load_ps(tris+32);
#else
		sf_m128 v1x = _mm_load_ps(tris+12);
		sf_m128 v1y = _mm_load_ps(tris+16);
		sf_m128 v1z = _mm_load_ps(tris+20);

		sf_m128 v2x = _mm_load_ps(tris+24);
		sf_m128 v2y = _mm_load_ps(tris+28);
		sf_m128 v2z = _mm_load_ps(tris+32);

		// Edge vectors
		sf_m128 e1x = _mm_sub_ps(v1x, v0x);
		sf_m128 e1y = _mm_sub_ps(v1y, v0y);
		sf_m128 e1z = _mm_sub_ps(v1z, v0z);

		sf_m128 e2x = _mm_sub_ps(v2x, v0x);
		sf_m128 e2y = _mm_sub_ps(v2y, v0y);
		sf_m128 e2z = _mm_sub_ps(v2z, v0z);
#endif
//		_mm_prefetch((const char *)(tris+36), _MM_HINT_T0);

		// begin calculating determinant - also used to calculate U parameter
		sf_m128 px = _mm_sub_ps(_mm_mul_ps(dY, e2z), _mm_mul_ps(dZ, e2y));
		sf_m128 py = _mm_sub_ps(_mm_mul_ps(dZ, e2x), _mm_mul_ps(dX, e2z));
		sf_m128 pz = _mm_sub_ps(_mm_mul_ps(dX, e2y), _mm_mul_ps(dY, e2x));

		// If det < 0, intersecting backfacing tri, > 0, intersecting frontfacing tri, 0, parallel to plane.
		sf_m128 det = _mm_add_ps(_mm_add_ps(_mm_mul_ps(e1x, px), _mm_mul_ps(e1y, py)), _mm_mul_ps(e1z, pz));

		// If determinant is near zero, ray lies in plane of triangle.

//		if (fabs(det) <= epsilon)
//			return FLOAT_INF;
		sf_m128 recipDet = _mm_rcp_ps(det);

		sf_m128 absdet = _mm_andnot_ps(sign_mask, det);
		sf_m128 out = _mm_cmple_ps(absdet, epsilon);

		// Calculate distance from v0 to ray origin
		sf_m128 tx = _mm_sub_ps(lX, v0x);
		sf_m128 ty = _mm_sub_ps(lY, v0y);
		sf_m128 tz = _mm_sub_ps(lZ, v0z);

		// Output barycentric u
		sf_m128 u = _mm_mul_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(tx, px), _mm_mul_ps(ty, py)), _mm_mul_ps(tz, pz)), recipDet);

//		if (u < 0.f || u > 1.f)
//			return FLOAT_INF; // Barycentric U is outside the triangle - early out.
		sf_m128 out2 = _mm_cmplt_ps(u, zero);
		out = _mm_or_ps(out, out2);
		out2 = _mm_cmpgt_ps(u, one);
		out = _mm_or_ps(out, out2);

		// Prepare to test V parameter
		sf_m128 qx = _mm_sub_ps(_mm_mul_ps(ty, e1z), _mm_mul_ps(tz, e1y));
		sf_m128 qy = _mm_sub_ps(_mm_mul_ps(tz, e1x), _mm_mul_ps(tx, e1z));
		sf_m128 qz = _mm_sub_ps(_mm_mul_ps(tx, e1y), _mm_mul_ps(ty, e1x));

		// Output barycentric v
		sf_m128 v = _mm_mul_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(dX, qx), _mm_mul_ps(dY, qy)), _mm_mul_ps(dZ, qz)), recipDet);

//		if (v < 0.f || u + v > 1.f) // Barycentric V or the combination of U and V are outside the triangle - no intersection.
//			return FLOAT_INF;
		out2 = _mm_cmplt_ps(v, zero);
		out = _mm_or_ps(out, out2);
		sf_m128 uv = _mm_add_ps(u, v);
		out2 = _mm_cmpgt_ps(uv, one);
		out = _mm_or_ps(out, out2);

		// Output signed distance from ray to triangle.
		sf_m128 t = _mm_mul_ps(_mm_add_ps(_mm_add_ps(_mm_mul_ps(e2x, qx), _mm_mul_ps(e2y, qy)), _mm_mul_ps(e2z, qz)), recipDet);

		// t < 0?
		out2 = _mm_cmplt_ps(t, zero);
		out = _mm_or_ps(out, out2);

		// Worse than previous result?
		out2 = _mm_cmpge_ps(t, nearestD);
		out = _mm_or_ps(out, out2);

		// The mask 'out' now contains 0xFF in all indices which are worse than previous, and
		// 0x00 in indices which are better.

#ifdef MATH_GEN_SSE41
		nearestD = _mm_blendv_ps(t, nearestD, out);
#else
		// If SSE 4.1 is not available:
		nearestD = _mm_and_ps(out, nearestD);
		t = _mm_andnot_ps(out, t);
		nearestD = _mm_or_ps(t, nearestD);
#endif

#ifdef MATH_GEN_UV

#ifdef MATH_GEN_SSE41
		nearestU = _mm_blendv_ps(u, nearestU, out); // 'blend' requires SSE4.1!
		nearestV = _mm_blendv_ps(v, nearestV, out); // 'blend' requires SSE4.1!
#else
		// If SSE 4.1 is not available:
		nearestU = _mm_and_ps(out, nearestU);
		nearestV = _mm_and_ps(out, nearestV);
		u = _mm_andnot_ps(out, u);
		v = _mm_andnot_ps(out, v);
		nearestU = _mm_or_ps(u, nearestU);
		nearestV = _mm_or_ps(v, nearestV);
#endif

#endif

#ifdef MATH_GEN_TRIANGLEINDEX
		__m128i hitIndex = _mm_set1_epi32(i);
#ifdef MATH_GEN_SSE41
		nearestIndex = _mm_castps_si128(_mm_blendv_ps(_mm_castsi128_ps(hitIndex), _mm_castsi128_ps(nearestIndex), out)); // 'blend' requires SSE4.1!
#else
		// If SSE 4.1 is not available:
		// Store the index of the triangle that was hit.
		nearestIndex = _mm_and_si128(_mm_castps_si128(out), nearestIndex);
		hitIndex = _mm_andnot_si128(_mm_castps_si128(out), hitIndex);
		nearestIndex = _mm_or_si128(hitIndex, nearestIndex);
#endif

#endif

		tris += 36;
	}

	float ds[16];
	float *alignedDS = (float*)(((uintptr_t)ds + 0xF) & ~0xF);

#ifdef MATH_GEN_UV
	float su[16];
	float *alignedU = (float*)(((uintptr_t)su + 0xF) & ~0xF);

	float sv[16];
	float *alignedV = (float*)(((uintptr_t)sv + 0xF) & ~0xF);

	_mm_store_ps(alignedU, nearestU);
	_mm_store_ps(alignedV, nearestV);
#endif

#ifdef MATH_GEN_TRIANGLEINDEX
	sf_u32 ds2[16];
	sf_u32 *alignedDS2 = (sf_u32*)(((uintptr_t)ds2 + 0xF) & ~0xF);

	_mm_store_si128((__m128i*)alignedDS2, nearestIndex);
#endif

	_mm_store_ps(alignedDS, nearestD);

	float smallestT = CMath::INFINITY_FLOAT;
//	float u = FLOAT_NAN, v = FLOAT_NAN;
	for(int i = 0; i < 4; ++i)
		if (alignedDS[i] < smallestT)
		{
			smallestT = alignedDS[i];
#ifdef MATH_GEN_TRIANGLEINDEX
			outTriangleIndex = alignedDS2[i]+i;
#endif

#ifdef MATH_GEN_UV
			outU = alignedU[i];
			outV = alignedV[i];
#endif
		}

//	TRACEEND(RayTriMeshIntersectSSE);

//	static double avgtimes = 0.f;
//	static double nAvgTimes = 0;
//	static double processedBytes;

//	processedBytes += numTris * 3 * 4;

//	avgtimes += Clock::TicksToMillisecondsD(time_RayTriMeshIntersectSSE);
//	++nAvgTimes;
//	std::cout << "Total avg (SSE): " << (avgtimes / nAvgTimes) << std::endl;
//	std::cout << "Hit distance (SSE): " << smallestT << ", index: " << hitTriangleIndex << ", UV: (" << u << ", " << v << ")" << std::endl;
//	std::cout << "(SSE) " << processedBytes / avgtimes * 1000.0 / 1024.0 / 1024.0 / 1024.0 << "GB/sec." << std::endl;

	return smallestT;
}


//float TriangleMesh::IntersectRay_AVX(const Ray &ray) const
//float TriangleMesh::IntersectRay_TriangleIndex_AVX(const Ray &ray, int &outIndex) const
//float TriangleMesh::IntersectRay_TriangleIndex_UV_AVX(const Ray &ray, int &outIndex, float &outU, float &outV) const




#ifdef MATH_GEN_SSE2
#undef MATH_GEN_SSE2
#endif
#ifdef MATH_GEN_SSE41
#undef MATH_GEN_SSE41
#endif
#ifdef MATH_GEN_TRIANGLEINDEX
#undef MATH_GEN_TRIANGLEINDEX
#endif
#ifdef MATH_GEN_UV
#undef MATH_GEN_UV
#endif

//} //end MATH
//} //end SMF
