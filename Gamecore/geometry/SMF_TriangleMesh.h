#ifndef __SMF_TRIANGLEMESH_
#define __SMF_TRIANGLEMESH_

#include "../SMF_Config.h"
#include "SMF_TriangleMesh.h"
#include "SMF_Triangle.h"
#include "SMF_Polyhedron.h"
#include "SMF_Polygon.h"
#include "../math/all.h"

namespace SMF {
using namespace MATH;
using namespace Util;
namespace MATH{
class CSIMD_SSE2;
class CSIMD_SSE41;
class CSIMD_AVX;

}
namespace GEO{

class CTriangle;
class CPolyhedron;
class CRay;


/**
 * \class CTriangleMesh
 *
 * \ingroup SMF_Geometric
 *
 * \if pt_br
 * \brief   Armazena uma malha de triângulos em uma matriz plana, otimizada para intersecção de raios
 * \elseif us_en
 * \brief 	This class stores a triangle mesh as flat array, optimized for ray intersections.
   \note    Represents an unindiced triangle mesh.
 * \endif
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
class SMF_API CTriangleMesh
{
public:
	CTriangleMesh();
	~CTriangleMesh();

	/// Specifies the vertex data of this triangle mesh. Replaces any old
	/// specified geometry.
	void set(const float *triangleMesh, int numTriangles);
	void set(const CVec3D *triangleMesh, int numTriangles) { set(reinterpret_cast<const float *>(triangleMesh), numTriangles); }
	void set(const CTriangle *triangleMesh, int numTriangles) { set(reinterpret_cast<const float *>(triangleMesh), numTriangles); }

	void set(const CPolyhedron &polyhedron);

	float IntersectRay(const CRay &ray) const;
	float IntersectRay_TriangleIndex(const CRay &ray, int &outTriangleIndex) const;
	float IntersectRay_TriangleIndex_UV(const CRay &ray, int &outTriangleIndex, float &outU, float &outV) const;

	void SetAoS(const float *vertexData, int numTriangles);
	void SetSoA4(const float *vertexData, int numTriangles);
	void SetSoA8(const float *vertexData, int numTriangles);

	float IntersectRay_TriangleIndex_UV_CPP(const CRay &ray, int &outTriangleIndex, float &outU, float &outV) const;
#ifdef _DEBUG
    SMF_INLINE int getVertexLayout()const {return vertexDataLayout;}
#endif
    SMF_INLINE float *toFloatPtr()  const {return data;}
    SMF_INLINE float *toFloatPtr() {return data;}
    SMF_INLINE int getCount()const {return numTriangles;}
private:
	float *data;
#ifdef _DEBUG
	int vertexDataLayout; // 0 - AoS, 1 - SoA4, 2 - SoA8
#endif
	int numTriangles;
	void ReallocVertexBuffer(int numTriangles);
};


} //end GEO
}  //end SMF

#endif // __SMF_TRIANGLEMESH_
