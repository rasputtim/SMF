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

#include "SMF_Config.h"
#include "geometry/SMF_TriangleMesh.h"
#include "geometry/SMF_Circle.h"
#include "math/all.h"
#include "math/SMF_SimdAVX.h"
#include "geometry/all.h"

namespace SMF {
using namespace MATH;
using namespace Util;
namespace GEO{


// If defined, we preprocess our CTriangleMesh data structure to contain (v0, v1-v0, v2-v0)
// instead of (v0, v1, v2) triplets for faster ray-triangle mesh intersection.
#define SOA_HAS_EDGES

const int simdCapability = System::getProcessorId();

CTriangleMesh::CTriangleMesh()
:data(0), numTriangles(0)
{

}

CTriangleMesh::~CTriangleMesh()
{
	mem_alignedFree(data);
}

void CTriangleMesh::set(const CPolyhedron &polyhedron)
{
	std::vector<CTriangle> tris = polyhedron.triangulate();
	if (!tris.empty())
	{
		int alignment = (SMF_CONTAINS(simdCapability,CPUID_AVX)) ? 8 : ((SMF_CONTAINS(simdCapability,CPUID_SSE41) || SMF_CONTAINS(simdCapability,CPUID_SSE2)) ? 4 : 1);
		CVec3D degen(-CMath::INFINITY_FLOAT, -CMath::INFINITY_FLOAT, -CMath::INFINITY_FLOAT);
		CTriangle degent(degen, degen, degen);
		while(tris.size() % alignment != 0)
			tris.push_back(degent);
		set(&tris[0], tris.size());
	}
}

void CTriangleMesh::set(const float *triangleMesh, int numTriangles)
{
	if (SMF_CONTAINS(simdCapability, CPUID_AVX))
		SetSoA8(triangleMesh, numTriangles);
	else if (SMF_CONTAINS(simdCapability,CPUID_SSE41) || SMF_CONTAINS(simdCapability,CPUID_SSE2))
		SetSoA4(triangleMesh, numTriangles);
	else
		SetAoS(triangleMesh, numTriangles);
}

float CTriangleMesh::IntersectRay(const CRay &ray) const
{
if(CMath::MATH_AUTOMATIC_SSE){

	if (SMF_CONTAINS(simdCapability,CPUID_AVX))
		return CSIMD_AVX::intersectRay(*this,ray);
	if (SMF_CONTAINS(simdCapability,CPUID_SSE41))
		return CSIMD_SSE41::intersectRay(*this,ray);
	if (SMF_CONTAINS(simdCapability,CPUID_SSE2))
		return CSIMD_SSE2::intersectRay(*this,ray);

}else{

	int triangleIndex;
	float u, v;
	return IntersectRay_TriangleIndex_UV_CPP(ray, triangleIndex, u, v);
}
}

float CTriangleMesh::IntersectRay_TriangleIndex(const CRay &ray, int &outTriangleIndex) const
{
if(CMath::MATH_AUTOMATIC_SSE){
	if (SMF_CONTAINS(simdCapability ,CPUID_AVX))
		return CSIMD_AVX::intersectRay_TriangleIndex(*this,ray, outTriangleIndex);
	if (SMF_CONTAINS(simdCapability , CPUID_SSE41))
		return CSIMD_SSE41::intersectRay_TriangleIndex(*this,ray, outTriangleIndex);
	if (SMF_CONTAINS(simdCapability , CPUID_SSE2))
		return CSIMD_SSE2::intersectRay_TriangleIndex(*this,ray, outTriangleIndex);
}else{	
	float u, v;
	return IntersectRay_TriangleIndex_UV_CPP(ray, outTriangleIndex, u, v);
}
}

float CTriangleMesh::IntersectRay_TriangleIndex_UV(const CRay &ray, int &outTriangleIndex, float &outU, float &outV) const
{
if(CMath::MATH_AUTOMATIC_SSE){
	if (SMF_CONTAINS(simdCapability, CPUID_AVX))
		return CSIMD_AVX::intersectRay_TriangleIndex_UV(*this,ray, outTriangleIndex, outU, outV);
	if (SMF_CONTAINS(simdCapability,CPUID_SSE41))
		return CSIMD_SSE41::intersectRay_TriangleIndex_UV(*this,ray, outTriangleIndex, outU, outV);
	if (SMF_CONTAINS(simdCapability,CPUID_SSE2))
		return CSIMD_SSE2::intersectRay_TriangleIndex_UV(*this,ray, outTriangleIndex, outU, outV);
}else{

	return IntersectRay_TriangleIndex_UV_CPP(ray, outTriangleIndex, outU, outV);
}
}
void CTriangleMesh::ReallocVertexBuffer(int numTris)
{
	mem_alignedFree(data);
	data = (float*)mem_alignedMalloc(numTris*3*3*sizeof(float), 32);
	numTriangles = numTris;
}

void CTriangleMesh::SetAoS(const float *vertexData, int numTriangles)
{
	ReallocVertexBuffer(numTriangles);
#ifdef _DEBUG
	vertexDataLayout = 0; // AoS
#endif

	memcpy(data, vertexData, numTriangles*3*3*4);
}

void CTriangleMesh::SetSoA4(const float *vertexData, int numTriangles)
{
	ReallocVertexBuffer(numTriangles);
#ifdef _DEBUG
	vertexDataLayout = 1; // SoA4
#endif

	SMF_ASSERT(numTriangles % 4 == 0); // We must have an evenly divisible amount of triangles, so that the SoA swizzling succeeds.

	// From (xyz xyz xyz) (xyz xyz xyz) (xyz xyz xyz) (xyz xyz xyz)
	// To xxxx yyyy zzzz xxxx yyyy zzzz xxxx yyyy zzzz

	float *o = data;
	for(int i = 0; i + 4 <= numTriangles; i += 4)
	{
		for(int j = 0; j < 9; ++j)
		{
			*o++ = vertexData[0];
			*o++ = vertexData[9];
			*o++ = vertexData[18];
			*o++ = vertexData[27];
			++vertexData;
		}
		vertexData += 9 * 3;
	}

#ifdef SOA_HAS_EDGES
	o = data;
	for(int i = 0; i + 4 <= numTriangles; i += 4)
	{
		for(int j = 12; j < 24; ++j)
			o[j] -= o[j-12];
		for(int j = 24; j < 36; ++j)
			o[j] -= o[j-24];
		o += 36;
	}
#endif
}

void CTriangleMesh::SetSoA8(const float *vertexData, int numTriangles)
{
	ReallocVertexBuffer(numTriangles);
#ifdef _DEBUG
	vertexDataLayout = 2; // SoA8
#endif

	SMF_ASSERT(numTriangles % 8 == 0); // We must have an evenly divisible amount of triangles, so that the SoA swizzling succeeds.

	// From (xyz xyz xyz) (xyz xyz xyz) (xyz xyz xyz) (xyz xyz xyz)
	// To xxxxxxxx yyyyyyyy zzzzzzzz xxxxxxxx yyyyyyyy zzzzzzzz xxxxxxxx yyyyyyyy zzzzzzzz

	float *o = data;
	for(int i = 0; i + 8 <= numTriangles; i += 8)
	{
		for(int j = 0; j < 9; ++j)
		{
			*o++ = vertexData[0];
			*o++ = vertexData[9];
			*o++ = vertexData[18];
			*o++ = vertexData[27];
			*o++ = vertexData[36];
			*o++ = vertexData[45];
			*o++ = vertexData[54];
			*o++ = vertexData[63];
			++vertexData;
		}
		vertexData += 9 * 7;
	}

#ifdef SOA_HAS_EDGES
	o = data;
	for(int i = 0; i + 8 <= numTriangles; i += 8)
	{
		for(int j = 24; j < 48; ++j)
			o[j] -= o[j-24];
		for(int j = 48; j < 72; ++j)
			o[j] -= o[j-48];
		o += 72;
	}
#endif
}

float CTriangleMesh::IntersectRay_TriangleIndex_UV_CPP(const CRay &ray, int &outTriangleIndex, float &outU, float &outV) const
{
	SMF_ASSERT(sizeof(CVec3D) == 3*sizeof(float));
	SMF_ASSERT(sizeof(CTriangle) == 3*sizeof(CVec3D));
#ifdef _DEBUG
	SMF_ASSERT(vertexDataLayout == 0); // Must be AoS structured!
#endif

	float nearestD = CMath::INFINITY_FLOAT;
	float u, v, d;

	const CTriangle *tris = reinterpret_cast<const CTriangle*>(data);
	for(int i = 0; i < numTriangles; ++i)
	{
		d = CTriangle::intersectLineTri(ray.pos, ray.dir, tris->a, tris->b, tris->c, u, v);
		if (d >= 0.f && d < nearestD)
		{
			nearestD = d;
			outU = u;
			outV = v;
			outTriangleIndex = i;
		}
		++tris;
	}

	return nearestD;
}





} //end GEO
}  //end SMF
