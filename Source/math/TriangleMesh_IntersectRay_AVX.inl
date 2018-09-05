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


/** \file TriangleMesh_IntersectRay_AVX.inl
	\author Jukka Jyl�nki
	\brief AVX implementation of ray-mesh intersection routines.
*/


#if !defined(MATH_GEN_TRIANGLEINDEX)
float CSIMD_AVX::intersectRay(const CTriangleMesh &trig, const CRay &ray)
#elif defined(MATH_GEN_TRIANGLEINDEX) && !defined(MATH_GEN_UV)
float CSIMD_AVX::intersectRay_TriangleIndex(const CTriangleMesh &trig, const CRay &ray, int &outTriangleIndex)
#elif defined(MATH_GEN_TRIANGLEINDEX) && defined(MATH_GEN_UV)
float CSIMD_AVX::intersectRay_TriangleIndex_UV(const CTriangleMesh &trig, const CRay &ray, int &outTriangleIndex, float &outU, float &outV)
#endif
{
//	std::cout << numTris << " tris: ";
//	TRACESTART(RayTriMeshIntersectAVX);

	SMF_ASSERT(sizeof(CVec3D) == 3*sizeof(float));
	SMF_ASSERT(sizeof(CTriangle) == 3*sizeof(CVec3D));
#ifdef _DEBUG
	SMF_ASSERT(trig.getVertexLayout() == 2); // Must be SoA8 structured!
#endif

//	hitTriangleIndex = -1;
//	CVec3D pt;
	const float inf = CMath::INFINITY_FLOAT;
	__m256 nearestD = _mm256_set1_ps(inf);
#ifdef MATH_GEN_UV
	__m256 nearestU = _mm256_set1_ps(inf);
	__m256 nearestV = _mm256_set1_ps(inf);
#endif
#ifdef MATH_GEN_TRIANGLEINDEX
	__m256i nearestIndex = _mm256_set1_epi32(-1);
#endif

	const __m256 lX = _mm256_broadcast_ss(&ray.pos.x);
	const __m256 lY = _mm256_broadcast_ss(&ray.pos.y);
	const __m256 lZ = _mm256_broadcast_ss(&ray.pos.z);

	const __m256 dX = _mm256_broadcast_ss(&ray.dir.x);
	const __m256 dY = _mm256_broadcast_ss(&ray.dir.y);
	const __m256 dZ = _mm256_broadcast_ss(&ray.dir.z);

	const __m256 epsilon = _mm256_set1_ps(1e-4f);
	const __m256 zero = _mm256_setzero_ps();
	const __m256 one = _mm256_set1_ps(1.f);

	SMF_ASSERT(((uintptr_t)trig.toFloatPtr() & 0x1F) == 0);

	const float *tris = reinterpret_cast<const float*>(trig.toFloatPtr());

	for(int i = 0; i+8 <= trig.getCount(); i += 8)
	{
		__m256 v0x = _mm256_load_ps(tris);
		__m256 v0y = _mm256_load_ps(tris+8);
		__m256 v0z = _mm256_load_ps(tris+16);

#ifdef SOA_HAS_EDGES
		// Edge vectors
		__m256 e1x = _mm256_load_ps(tris+24);
		__m256 e1y = _mm256_load_ps(tris+32);
		__m256 e1z = _mm256_load_ps(tris+40);

		__m256 e2x = _mm256_load_ps(tris+48);
		__m256 e2y = _mm256_load_ps(tris+56);
		__m256 e2z = _mm256_load_ps(tris+64);
#else
		__m256 v1x = _mm256_load_ps(tris+24);
		__m256 v1y = _mm256_load_ps(tris+32);
		__m256 v1z = _mm256_load_ps(tris+40);

		__m256 v2x = _mm256_load_ps(tris+48);
		__m256 v2y = _mm256_load_ps(tris+56);
		__m256 v2z = _mm256_load_ps(tris+64);

		// Edge vectors
		__m256 e1x = _mm256_sub_ps(v1x, v0x);
		__m256 e1y = _mm256_sub_ps(v1y, v0y);
		__m256 e1z = _mm256_sub_ps(v1z, v0z);

		__m256 e2x = _mm256_sub_ps(v2x, v0x);
		__m256 e2y = _mm256_sub_ps(v2y, v0y);
		__m256 e2z = _mm256_sub_ps(v2z, v0z);
#endif
		// begin calculating determinant - also used to calculate U parameter
		__m256 px = _mm256_sub_ps(_mm256_mul_ps(dY, e2z), _mm256_mul_ps(dZ, e2y));
		__m256 py = _mm256_sub_ps(_mm256_mul_ps(dZ, e2x), _mm256_mul_ps(dX, e2z));
		__m256 pz = _mm256_sub_ps(_mm256_mul_ps(dX, e2y), _mm256_mul_ps(dY, e2x));

		// If det < 0, intersecting backfacing tri, > 0, intersecting frontfacing tri, 0, parallel to plane.
		__m256 det = _mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(e1x, px), _mm256_mul_ps(e1y, py)), _mm256_mul_ps(e1z, pz));

		// If determinant is near zero, ray lies in plane of triangle.

//		if (fabs(det) <= epsilon)
//			return CMath::INFINITY_FLOAT;
		__m256 recipDet = _mm256_rcp_ps(det);

		__m256 absdet = abs_ps256(det);
		__m256 out = _mm256_cmp_ps(absdet, epsilon, _CMP_LT_OQ);

		// Calculate distance from v0 to ray origin
		__m256 tx = _mm256_sub_ps(lX, v0x);
		__m256 ty = _mm256_sub_ps(lY, v0y);
		__m256 tz = _mm256_sub_ps(lZ, v0z);

		// Output barycentric u
		__m256 u = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(tx, px), _mm256_mul_ps(ty, py)), _mm256_mul_ps(tz, pz)), recipDet);

//		if (u < 0.f || u > 1.f)
//			return CMath::INFINITY_FLOAT; // Barycentric U is outside the triangle - early out.
		__m256 out2 = _mm256_cmp_ps(u, zero, 1);
		out = _mm256_or_ps(out, out2);
		out2 = _mm256_cmp_ps(u, one, _CMP_GT_OQ);
		out = _mm256_or_ps(out, out2);

		// Prepare to test V parameter
		__m256 qx = _mm256_sub_ps(_mm256_mul_ps(ty, e1z), _mm256_mul_ps(tz, e1y));
		__m256 qy = _mm256_sub_ps(_mm256_mul_ps(tz, e1x), _mm256_mul_ps(tx, e1z));
		__m256 qz = _mm256_sub_ps(_mm256_mul_ps(tx, e1y), _mm256_mul_ps(ty, e1x));

		// Output barycentric v
		__m256 v = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(dX, qx), _mm256_mul_ps(dY, qy)), _mm256_mul_ps(dZ, qz)), recipDet);

//		if (v < 0.f || u + v > 1.f) // Barycentric V or the combination of U and V are outside the triangle - no intersection.
//			return CMath::INFINITY_FLOAT;
		out2 = _mm256_cmp_ps(v, zero, _CMP_LT_OQ);
		out = _mm256_or_ps(out, out2);
		__m256 uv = _mm256_add_ps(u, v);
		out2 = _mm256_cmp_ps(uv, one, _CMP_GT_OQ);
		out = _mm256_or_ps(out, out2);

		// Output signed distance from ray to triangle.
		__m256 t = _mm256_mul_ps(_mm256_add_ps(_mm256_add_ps(_mm256_mul_ps(e2x, qx), _mm256_mul_ps(e2y, qy)), _mm256_mul_ps(e2z, qz)), recipDet);

		// t < 0?
		out2 = _mm256_cmp_ps(t, zero, _CMP_LT_OQ);
		out = _mm256_or_ps(out, out2);

		// Worse than previous result?
		out2 = _mm256_cmp_ps(t, nearestD, _CMP_GE_OQ);
		out = _mm256_or_ps(out, out2);

		// Store the index of the triangle that was hit.
#ifdef MATH_GEN_TRIANGLEINDEX
		__m256i hitIndex = _mm256_set1_epi32(i);
		nearestIndex = _mm256_castps_si256(_mm256_blendv_ps(_mm256_castsi256_ps(hitIndex), _mm256_castsi256_ps(nearestIndex), out)); // 'blend' requires SSE4.1!
#endif

#ifdef MATH_GEN_UV
		nearestU = _mm256_blendv_ps(u, nearestU, out);
		nearestV = _mm256_blendv_ps(v, nearestV, out);
#endif

		// The mask out now contains 0xFF in all indices which are worse than previous, and
		// 0x00 in indices which are better.
		nearestD = _mm256_blendv_ps(t, nearestD, out);

		tris += 72;
	}

	float ds[32];
	float *alignedDS = (float*)(((uintptr_t)ds + 0x1F) & ~0x1F);

#ifdef MATH_GEN_UV
	float su[32];
	float *alignedU = (float*)(((uintptr_t)su + 0x1F) & ~0x1F);

	float sv[32];
	float *alignedV = (float*)(((uintptr_t)sv + 0x1F) & ~0x1F);

	_mm256_store_ps(alignedU, nearestU);
	_mm256_store_ps(alignedV, nearestV);
#endif

#ifdef MATH_GEN_TRIANGLEINDEX
	sf_u32 ds2[32];
	sf_u32 *alignedDS2 = (sf_u32*)(((uintptr_t)ds2 + 0x1F) & ~0x1F);

	_mm256_store_si256((__m256i*)alignedDS2, nearestIndex);
#endif

	_mm256_store_ps(alignedDS, nearestD);

	float smallestT = CMath::INFINITY_FLOAT;
//	float u = FLOAT_NAN, v = FLOAT_NAN;
	for(int i = 0; i < 8; ++i)
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

//	TRACEEND(RayTriMeshIntersectAVX);

//	static double avgtimes = 0.f;
//	static double nAvgTimes = 0;
//	static double processedBytes;

//	processedBytes += numTris * 3 * 4;

//	avgtimes += Clock::TicksToMillisecondsD(time_RayTriMeshIntersectAVX);
//	++nAvgTimes;
//	std::cout << "Total avg (AVX): " << (avgtimes / nAvgTimes) << std::endl;
//	std::cout << "Hit distance (AVX): " << smallestT << ", index: " << hitTriangleIndex << ", UV: (" << u << ", " << v << ")" << std::endl;
//	std::cout << "(AVX) " << processedBytes / avgtimes * 1000.0 / 1024.0 / 1024.0 / 1024.0 << "GB/sec." << std::endl;

	return smallestT;
}


#ifdef MATH_GEN_TRIANGLEINDEX
#undef MATH_GEN_TRIANGLEINDEX
#endif
#ifdef MATH_GEN_UV
#undef MATH_GEN_UV
#endif

