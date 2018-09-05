/*
  SMF -  Super Math Fabric  
  Copyright (C) 2014 Rasputtim <Rasputtim@hotmail.com>

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

  */

#include "math/SMF_JointTransform.h"
#include "math/SMF_Transforms.h"
//#pragma hdrstop

namespace SMF{

namespace MATH{

void CMatJoint3x4::set3x3Part(const CMat3D &r)
{
	//s assume(r.rsFinite());
	mat[0*4+0] = r[0][0]; mat[0*4+1] = r[0][1]; mat[0*4+2] = r[0][2];
	mat[1*4+0] = r[1][0]; mat[1*4+1] = r[1][1]; mat[1*4+2] = r[1][2];
	mat[2*4+0] = r[2][0]; mat[2*4+1] = r[2][1]; mat[2*4+2] = r[2][2];
}

void CMatJoint3x4::setMatrixRotatePart(const CQuaternion &q)
{
	// See e.g. http://www.geometrictools.com/Documentation/LinearAlgebraicQuaternions.pdf .

	//assume(q.isNormalized());
	const float x = q.x; const float y = q.y; const float z = q.z; const float w = q.w;
	mat[0*4+0] = 1 - 2*(y*y + z*z); mat[0*4+1] =     2*(x*y - z*w); mat[0*4+2] =     2*(x*z + y*w);
	mat[1*4+0] =     2*(x*y + z*w); mat[1*4+1] = 1 - 2*(x*x + z*z); mat[1*4+2] =     2*(y*z - x*w);
	mat[2*4+0] =     2*(x*z - y*w); mat[2*4+1] =     2*(y*z + x*w); mat[2*4+2] = 1 - 2*(x*x + y*y);
}


void CMatJoint3x4::set3x3PartRotateX(float angle)
{
	/*
	 ³   1   0   0 ³
	 ³   0  cz -sz ³
	 ³   0  sz  cz ³
	*/

	const float cosz = CMath::cos(angle);
	const float sinz = CMath::sin(angle);

	mat[0*4+0] = 1.f; mat[0*4+1] =  0.f; mat[0*4+2] =   0.f;
	mat[1*4+0] = 0.f; mat[1*4+1] = cosz; mat[1*4+2] = -sinz;
	mat[2*4+0] = 0.f; mat[2*4+1] = sinz; mat[2*4+2] =  cosz;
}

void CMatJoint3x4::set3x3PartRotateY(float angle)
{
	/*
	 ³  cz   0  sz ³
	 ³   0   1   0 ³
	 ³ -sz   0  cz ³
	*/

	const float cosz = CMath::cos(angle);
	const float sinz = CMath::sin(angle);

	mat[0*4+0] =  cosz; mat[0*4+1] = 0.f; mat[0*4+2] = sinz;
	mat[1*4+0] =   0.f; mat[1*4+1] = 1.f; mat[1*4+2] =  0.f;
	mat[2*4+0] = -sinz; mat[2*4+1] = 0.f; mat[2*4+2] = cosz;
}

void CMatJoint3x4::set3x3PartRotateZ(float angle)
{
	/*
	 ³  cz -sz   0 ³
	 ³  sz  cz   0 ³
	 ³   0   0   1 ³
	*/

	const float cosz = CMath::cos(angle);
	const float sinz = CMath::sin(angle);

	mat[0*4+0] = cosz; mat[0*4+1] = -sinz; mat[0*4+2] = 0.f;
	mat[1*4+0] = sinz; mat[1*4+1] =  cosz; mat[1*4+2] = 0.f;
	mat[2*4+0] =  0.f; mat[2*4+1] =   0.f; mat[2*4+2] = 1.f;
}

void CMatJoint3x4::setRotatePartX(float angle)
{
	set3x3PartRotateX(angle);
}

void CMatJoint3x4::setRotatePartY(float angle)
{
	set3x3PartRotateY(angle);
}

void CMatJoint3x4::setRotatePartZ(float angle)
{
	set3x3PartRotateZ(angle);
}

void CMatJoint3x4::setRotatePart(const CVec3D &axisDirection, float angle)
{
	setRotatePart(CQuaternion(axisDirection, angle));
}

void CMatJoint3x4::setRotatePart(const CQuaternion &q)
{
	setMatrixRotatePart(q);
}


CMatJoint3x4::CMatJoint3x4(const CQuaternion &orientation)
{
	setRotatePart(orientation);
	setTranslatePart(0, 0, 0);
}

CMatJoint3x4::CMatJoint3x4(const CQuaternion &orientation, const CVec3D &translation)
{
	setRotatePart(orientation);
	setTranslatePart(translation);
}
CJointQuaternion CMatJoint3x4::ToJointQuat( void ) const {
	CJointQuaternion	jq;
	float		trace;
	float		s;
	float		t;
	int     	i;
	int			j;
	int			k;

	static int 	next[3] = { 1, 2, 0 };

	trace = mat[0 * 4 + 0] + mat[1 * 4 + 1] + mat[2 * 4 + 2];

	if ( trace > 0.0f ) {

		t = trace + 1.0f;
		s = MATH::CMath::invSqrt( t ) * 0.5f;

		jq.q[3] = s * t;
		jq.q[0] = ( mat[1 * 4 + 2] - mat[2 * 4 + 1] ) * s;
		jq.q[1] = ( mat[2 * 4 + 0] - mat[0 * 4 + 2] ) * s;
		jq.q[2] = ( mat[0 * 4 + 1] - mat[1 * 4 + 0] ) * s;

	} else {

		i = 0;
		if ( mat[1 * 4 + 1] > mat[0 * 4 + 0] ) {
			i = 1;
		}
		if ( mat[2 * 4 + 2] > mat[i * 4 + i] ) {
			i = 2;
		}
		j = next[i];
		k = next[j];

		t = ( mat[i * 4 + i] - ( mat[j * 4 + j] + mat[k * 4 + k] ) ) + 1.0f;
		s = MATH::CMath::invSqrt( t ) * 0.5f;

		jq.q[i] = s * t;
		jq.q[3] = ( mat[j * 4 + k] - mat[k * 4 + j] ) * s;
		jq.q[j] = ( mat[i * 4 + j] + mat[j * 4 + i] ) * s;
		jq.q[k] = ( mat[i * 4 + k] + mat[k * 4 + i] ) * s;
	}

	jq.t[0] = mat[0 * 4 + 3];
	jq.t[1] = mat[1 * 4 + 3];
	jq.t[2] = mat[2 * 4 + 3];

	return jq;
}

CVec3D CMatJoint3x4::transformPos(const CVec3D &pointVector) const
{
if(CMath::MATH_AUTOMATIC_SSE){
	return CSIMD::getProcessor()->mat3x4_mul_vec(toM128Ptr(), MATH::SET_PS(1.f, pointVector.z, pointVector.y, pointVector.x));
}else{
	return transformPos(pointVector.x, pointVector.y, pointVector.z);
}
}

CVec3D CMatJoint3x4::transformPos(float x, float y, float z) const
{
if(CMath::MATH_AUTOMATIC_SSE){
	return CSIMD::getProcessor()->mat3x4_mul_vec(toM128Ptr(), MATH::SET_PS(1, z, y, x));
}else{
	CVec3D lin0(mat[0*4+1],mat[0*4+2],mat[0*4+3]);
	CVec3D lin1(mat[1*4+1],mat[1*4+2],mat[1*4+3]);
	CVec3D lin2(mat[2*4+1],mat[2*4+2],mat[2*4+3]);

	return CVec3D(DOT3_xyz(lin0, x,y,z) + mat[0*4+3],
				  DOT3_xyz(lin1, x,y,z) + mat[1*4+3],
				  DOT3_xyz(lin2, x,y,z) + mat[2*4+3]);
}
}

CVec3D CMatJoint3x4::transformDir(const CVec3D &directionVector) const
{
if(CMath::MATH_AUTOMATIC_SSE){

	return CSIMD::getProcessor()->mat3x4_mul_vec(&q[0], MATH::SET_PS(0.f, directionVector.z, directionVector.y, directionVector.x));
}else{
	return transformDir(directionVector.x, directionVector.y, directionVector.z);
}
}

CVec3D CMatJoint3x4::transformDir(float x, float y, float z) const
{
if(CMath::MATH_AUTOMATIC_SSE){
	return CSIMD::getProcessor()->mat3x4_mul_vec(toM128Ptr(), MATH::SET_PS(0, z, y, x));
}else{
	CVec3D lin0(mat[0*4+1],mat[0*4+2],mat[0*4+3]);
	CVec3D lin1(mat[1*4+1],mat[1*4+2],mat[1*4+3]);
	CVec3D lin2(mat[2*4+1],mat[2*4+2],mat[2*4+3]);

	return CVec3D(DOT3_xyz(lin0, x,y,z),
				  DOT3_xyz(lin1, x,y,z),
				  DOT3_xyz(lin2, x,y,z));
}
}

CVec4D CMatJoint3x4::transform(const CVec4D &vector) const
{
	if(CMath::MATH_AUTOMATIC_SSE) {
		return CVec4D(CSIMD::getProcessor()->mat3x4_mul_sse(Row(0).toM128Ptr(), vector.v));
	}else{
		return CVec4D(CSIMD::getGenProcessor()->mat3x4_mul_sse(Row(0).toM128Ptr(), vector.v));

	}
}

void CMatJoint3x4::setRow(int row, float m_r0, float m_r1, float m_r2, float m_r3)
{
	SMF_ASSERT(row >= 0);
	SMF_ASSERT(row < 3);
	SMF_ASSERT(MATH::isFinite(m_r0));
	SMF_ASSERT(MATH::isFinite(m_r1));
	SMF_ASSERT(MATH::isFinite(m_r2));
	SMF_ASSERT(MATH::isFinite(m_r3));
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (row < 0 || row >= 3)
		return; // Benign failure
#endif

if(CMath::MATH_AUTOMATIC_SSE){
///todo: put this code into correct CIMD Class	

	q[row] = MATH::SET_PS(m_r3, m_r2, m_r1, m_r0);
}else{
	mat[row*4+0] = m_r0;
	mat[row*4+1] = m_r1;
	mat[row*4+2] = m_r2;
	mat[row*4+3] = m_r3;
}
}

void CMatJoint3x4::setRow(int row, const CVec3D &rowVector, float w)
{
	setRow(row, rowVector.x, rowVector.y, rowVector.z, w);
}

void CMatJoint3x4::setRow(int row, const CVec4D &rowVector)
{
if(CMath::MATH_AUTOMATIC_SSE){
///todo: put this code into correct CIMD Class	

#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (row < 0 || row >= 3)
		return; // Benign failure
#endif
	q[row] = rowVector.v;
}else{
	setRow(row, rowVector.x, rowVector.y, rowVector.z, rowVector.w);
}
}

void CMatJoint3x4::setRow(int row, const float *data)
{
	SMF_ASSERT(data);
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (!data)
		return;
#endif

if(CMath::MATH_AUTOMATIC_SSE){
///todo: put this code into correct CIMD Class	

#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (row < 0 || row >= 3)
		return; // Benign failure
#endif
	this->q[row] = _mm_loadu_ps(data); // Assume unaligned load, since we don't know if data is 16-byte-aligned.
}else{
	setRow(row, data[0], data[1], data[2], data[3]);
}
}

void CMatJoint3x4::setCol(int column, float m_0c, float m_1c, float m_2c)
{
	SMF_ASSERT(column >= 0);
	SMF_ASSERT(column < 4);
	SMF_ASSERT(MATH::isFinite(m_0c));
	SMF_ASSERT(MATH::isFinite(m_1c));
	SMF_ASSERT(MATH::isFinite(m_2c));
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (column < 0 || column >= 4)
		return; // Benign failure
#endif
	mat[0*4+column] = m_0c;
	mat[1*4+column] = m_1c;
	mat[2*4+column] = m_2c;
}

void CMatJoint3x4::setCol(int column, const CVec3D &columnVector)
{
	setCol(column, columnVector.x, columnVector.y, columnVector.z);
}

void CMatJoint3x4::setCol(int column, const float *data)
{
	SMF_ASSERT(data);
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (!data)
		return;
#endif
	setCol(column, data[0], data[1], data[2]);
}
void CMatJoint3x4::set(float _00, float _01, float _02, float _03,
				   float _10, float _11, float _12, float _13,
				   float _20, float _21, float _22, float _23)
{
if(CMath::MATH_AUTOMATIC_SSE){
///todo: put this code into correct CIMD Class	
	q[0] = MATH::SET_PS(_03, _02, _01, _00);
	q[1] = MATH::SET_PS(_13, _12, _11, _10);
	q[2] = MATH::SET_PS(_23, _22, _21, _20);
}else{
	mat[0*4+0] = _00; mat[0*4+1] = _01; mat[0*4+2] = _02; mat[0*4+3] = _03;
	mat[1*4+0] = _10; mat[1*4+1] = _11; mat[1*4+2] = _12; mat[1*4+3] = _13;
	mat[2*4+0] = _20; mat[2*4+1] = _21; mat[2*4+2] = _22; mat[2*4+3] = _23;
}
}
CMatJoint3x4 &CMatJoint3x4::operator *=(float scalar)
{
if(CMath::MATH_AUTOMATIC_SSE){
//todo: put this code into CIMD Class
	sf_m128 s = _mm_set1_ps(scalar);
	q[0] = _mm_mul_ps(Row(0).v, s);
	q[1] = _mm_mul_ps(Row(1).v, s);
	q[2] = _mm_mul_ps(Row(2).v, s);
}else{
	for(int y = 0; y < 3; ++y)
		for(int x = 0; x < 4; ++x)
			mat[y * 4 + x] *= scalar;
}

	return *this;
}

CMatJoint3x4 &CMatJoint3x4::operator /=(float scalar)
{
	SMF_ASSERT(!CMath::equalsAbs(scalar, 0));

if(CMath::MATH_AUTOMATIC_SSE){
///todo: put this code into correct CIMD Class

	sf_m128 s = _mm_set1_ps(scalar);
	sf_m128 one = _mm_set1_ps(1.f);
	s = _mm_div_ps(one, s);
	q[0] = _mm_mul_ps(Row(0).v, s);
	q[1] = _mm_mul_ps(Row(1).v, s);
	q[2] = _mm_mul_ps(Row(2).v, s);
}else{
	float invScalar = 1.f / scalar;
	for(int y = 0; y < 3; ++y)
		for(int x = 0; x < 4; ++x)
			mat[y * 4 + x] *= invScalar;
}

	return *this;
}
CMatJoint3x4 CMatJoint3x4::operator *(float scalar) const
{
	CMatJoint3x4 r;
if(CMath::MATH_AUTOMATIC_SSE){
///todo: put this code into correct CIMD Class
	
	sf_m128 s = _mm_set1_ps(scalar);
	r.q[0] = _mm_mul_ps(q[0], s);
	r.q[1] = _mm_mul_ps(q[1], s);
	r.q[2] = _mm_mul_ps(q[2], s);
}else{
	r = *this;
	r *= scalar;
}

	return r;
}
bool CMatJoint3x4::isFinite() const
{
	for(int y = 0; y < 3; ++y)
		for(int x = 0; x < 4; ++x)
			if (!MATH::isFinite(mat[y*4+x]))
				return false;
	return true;
}

bool CMatJoint3x4::isIdentity(float epsilon) const
{
	for(int y = 0; y < 3; ++y)
		for(int x = 0; x < 4; ++x)
			if (!CMath::equalsAbs(mat[y*4+x], (x == y) ? 1.f : 0.f, epsilon))
				return false;

	return true;
}

bool CMatJoint3x4::isLowerTriangular(float epsilon) const
{
	return CMath::equalsAbs(mat[0*4+1], 0.f, epsilon)
		&& CMath::equalsAbs(mat[0*4+2], 0.f, epsilon)
		&& CMath::equalsAbs(mat[0*4+3], 0.f, epsilon)
		&& CMath::equalsAbs(mat[1*4+2], 0.f, epsilon)
		&& CMath::equalsAbs(mat[1*4+3], 0.f, epsilon)
		&& CMath::equalsAbs(mat[2*4+3], 0.f, epsilon);
}

bool CMatJoint3x4::isUpperTriangular(float epsilon) const
{
	return CMath::equalsAbs(mat[1*4+0], 0.f, epsilon)
		&& CMath::equalsAbs(mat[2*4+0], 0.f, epsilon)
		&& CMath::equalsAbs(mat[2*4+1], 0.f, epsilon);
}

bool CMatJoint3x4::isInvertible(float epsilon) const
{
	float d = determinant();
	bool isSingular = CMath::equalsAbs(d, 0.f, epsilon);
#ifdef TEST_FOR_CORRECTNESS
	CMat3D temp = CMat3DPart();
	SMF_ASSERT(temp.inverseSelf() != isSingular); // isInvertible() and Inverse() must match!
#endif
	return !isSingular;
}

bool CMatJoint3x4::isSymmetric(float epsilon) const
{
	return CMath::equalsAbs(mat[0*4+1], mat[1*4+0], epsilon) &&
		CMath::equalsAbs(mat[0*4+2], mat[2*4+0], epsilon) &&
		CMath::equalsAbs(mat[1*4+2], mat[2*4+1], epsilon);
}

bool CMatJoint3x4::isSkewSymmetric(float epsilon) const
{
	return CMath::equalsAbs(mat[0*4+0], 0.f, epsilon) &&
		CMath::equalsAbs(mat[1*4+1], 0.f, epsilon) &&
		CMath::equalsAbs(mat[2*4+2], 0.f, epsilon) &&
		CMath::equalsAbs(mat[0*4+1], -mat[1*4+0], epsilon) &&
		CMath::equalsAbs(mat[0*4+2], -mat[2*4+0], epsilon) &&
		CMath::equalsAbs(mat[1*4+2], -mat[2*4+1], epsilon);
}

bool CMatJoint3x4::hasUnitaryScale(float epsilon) const
{
	CVec3D scale = extractScale();
	return scale.compare(CVec3D(1.f, 1.f, 1.f), epsilon);
}

bool CMatJoint3x4::hasNegativeScale() const
{
	return determinant() < 0.f;
}

bool CMatJoint3x4::hasUniformScale(float epsilon) const
{
	CVec3D scale = extractScale();
	return CMath::equalsAbs(scale.x, scale.y, epsilon) && CMath::equalsAbs(scale.x, scale.z, epsilon);
}

bool CMatJoint3x4::isRowOrthogonal(float epsilon) const
{
	return Row(0).isPerpendicular3(Row(1), epsilon)
		&& Row(0).isPerpendicular3(Row(2), epsilon)
		&& Row(1).isPerpendicular3(Row(2), epsilon);
}

bool CMatJoint3x4::isColOrthogonal(float epsilon) const
{
	return Col(0).isPerpendicular(Col(1), epsilon)
		&& Col(0).isPerpendicular(Col(2), epsilon)
		&& Col(1).isPerpendicular(Col(2), epsilon);
}

bool CMatJoint3x4::isOrthonormal(float epsilon) const
{
	///\todo Epsilon magnitudes don't match.
	return isColOrthogonal(epsilon) && Row3(0).isNormalized(epsilon) && Row3(1).isNormalized(epsilon) && Row3(2).isNormalized(epsilon);
}

float CMatJoint3x4::determinant() const
{
	SMF_ASSERT(CMat3DPart().isFinite());
	if (CMath::MATH_AUTOMATIC_SSE) {
		return CSIMD::getProcessor()->mat3x4_determinant(Row(0).toM128Ptr());
	}else{
#if 0
		const float a = mat[0*4+0];
		const float b = mat[0*4+1];
		const float c = mat[0*4+2];
		const float d = mat[1*4+0];
		const float e = mat[1*4+1];
		const float f = mat[1*4+2];
		const float g = mat[2*4+0];
		const float h = mat[2*4+1];
		const float i = mat[2*4+2];

		return a*e*i + b*f*g + c*d*h - a*f*h - b*d*i - c*e*g;
#endif
	return CSIMD::getGenProcessor()->mat3x4_determinant(Row(0).toM128Ptr());
	}
}

bool CMatJoint3x4::inverseColOrthogonal()
{
	///\todo SSE.
#ifdef TEST_FOR_CORRECTNESS
	CMatJoint3x4 orig = *this;
#endif
	SMF_ASSERT(isColOrthogonal());
	float s1 = CVec3D(mat[0*4+0], mat[1*4+0], mat[2*4+0]).getLengthSqr();
	float s2 = CVec3D(mat[0*4+1], mat[1*4+1], mat[2*4+1]).getLengthSqr();
	float s3 = CVec3D(mat[0*4+2], mat[1*4+2], mat[2*4+2]).getLengthSqr();
	if (s1 < 1e-8f || s2 < 1e-8f || s3 < 1e-8f)
		return false;
	s1 = 1.f / s1;
	s2 = 1.f / s2;
	s3 = 1.f / s3;
	Swap(mat[0*4+1], mat[1*4+0]);
	Swap(mat[0*4+2], mat[2*4+0]);
	Swap(mat[1*4+2], mat[2*4+1]);

	mat[0*4+0] *= s1; mat[0*4+1] *= s1; mat[0*4+2] *= s1;
	mat[1*4+0] *= s2; mat[1*4+1] *= s2; mat[1*4+2] *= s2;
	mat[2*4+0] *= s3; mat[2*4+1] *= s3; mat[2*4+2] *= s3;

	setTranslatePart(transformDir(-mat[0*4+3], -mat[1*4+3], -mat[2*4+3]));
#ifdef TEST_FOR_CORRECTNESS
	SMF_ASSERT(!orig.isInvertible()|| (orig * *this).isIdentity());
	SMF_ASSERT(isRowOrthogonal());
#endif
	return true;
}

bool CMatJoint3x4::inverseOrthogonalUniformScale()
{
	///\todo SSE.
	SMF_ASSERT(isColOrthogonal(1e-3f));
	SMF_ASSERT(hasUniformScale());
	Swap(mat[0*4+1], mat[1*4+0]);
	Swap(mat[0*4+2], mat[2*4+0]);
	Swap(mat[1*4+2], mat[2*4+1]);
	float scale = CVec3D(mat[0*4+0], mat[1*4+0], mat[2*4+0]).getLengthSqr();
	if (scale == 0.f)
		return false;
	scale = 1.f / scale;

	mat[0*4+0] *= scale; mat[0*4+1] *= scale; mat[0*4+2] *= scale;
	mat[1*4+0] *= scale; mat[1*4+1] *= scale; mat[1*4+2] *= scale;
	mat[2*4+0] *= scale; mat[2*4+1] *= scale; mat[2*4+2] *= scale;

	setTranslatePart(transformDir(-mat[0*4+3], -mat[1*4+3], -mat[2*4+3]));

	return true;
}

void CMatJoint3x4::inverseOrthonormal()
{
	SMF_ASSERT(isOrthonormal());

	if (CMath::MATH_AUTOMATIC_SSE) {
		CSIMD::getProcessor()->mat3x4_inverse_orthonormal(Row(0).toM128Ptr(),Row(0).toM128Ptr());
	}else{
	CSIMD::getGenProcessor()->mat3x4_inverse_orthonormal(Row(0).toM128Ptr(),Row(0).toM128Ptr());

	}
}

void CMatJoint3x4::batchTransformPos(CVec3D *pointArray, int numPoints) const
{
	SMF_ASSERT(pointArray);
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (!pointArray)
		return;
#endif
	for(int i = 0; i < numPoints; ++i)
		pointArray[i] = MulPos(pointArray[i]);
}

void CMatJoint3x4::batchTransformPos(CVec3D *pointArray, int numPoints, int stride) const
{
	SMF_ASSERT(pointArray);
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (!pointArray)
		return;
#endif
	SMF_ASSERT(stride >= (int)sizeof(CVec3D));
	sf_u8 *data = reinterpret_cast<sf_u8*>(pointArray);
	for(int i = 0; i < numPoints; ++i)
	{
		CVec3D *v = reinterpret_cast<CVec3D*>(data + stride*i);
		*v = MulPos(*v);
	}
}

void CMatJoint3x4::batchTransformDir(CVec3D *dirArray, int numVectors) const
{
	SMF_ASSERT(dirArray);
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (!dirArray)
		return;
#endif
	for(int i = 0; i < numVectors; ++i)
		dirArray[i] = MulPos(dirArray[i]);
}

void CMatJoint3x4::batchTransformDir(CVec3D *dirArray, int numVectors, int stride) const
{
	SMF_ASSERT(dirArray);
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (!dirArray)
		return;
#endif
	SMF_ASSERT(stride >= (int)sizeof(CVec3D));
	sf_u8 *data = reinterpret_cast<sf_u8*>(dirArray);
	for(int i = 0; i < numVectors; ++i)
	{
		CVec3D *v = reinterpret_cast<CVec3D*>(data + stride*i);
		*v = MulDir(*v);
	}
}

void CMatJoint3x4::batchTransform(CVec4D *vectorArray, int numVectors) const
{
	SMF_ASSERT(vectorArray);
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (!vectorArray)
		return;
#endif
	for(int i = 0; i < numVectors; ++i)
		vectorArray[i] = (*this).mul(vectorArray[i]);
}

void CMatJoint3x4::batchTransform(CVec4D *vectorArray, int numVectors, int stride) const
{
	SMF_ASSERT(vectorArray);
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (!vectorArray)
		return;
#endif
	SMF_ASSERT(stride >= (int)sizeof(CVec4D));
	sf_u8 *data = reinterpret_cast<sf_u8*>(vectorArray);
	for(int i = 0; i < numVectors; ++i)
	{
		CVec4D *v = reinterpret_cast<CVec4D*>(data + stride*i);
		*v = (*this).mul(*v);
	}
}
CVec4D CMatJoint3x4::mul(const CVec4D &rhs) const
{
	return transform(rhs);
}

CTranslateOp CMatJoint3x4::translate(float tx, float ty, float tz)
{
	return CTranslateOp(tx, ty, tz);
}

CTranslateOp CMatJoint3x4::translate(const CVec3D &offset)
{
	return CTranslateOp(offset);
}

CScaleOp CMatJoint3x4::scale(float sx, float sy, float sz)
{
	return CScaleOp(sx, sy, sz);
}

CScaleOp CMatJoint3x4::scale(const CVec3D &scale)
{
	return CScaleOp(scale);
}
CMatJoint3x4 CMatJoint3x4::scale(const CVec3D &scale, const CVec3D &scaleCenter)
{
	return CMatJoint3x4::translate(scaleCenter) * CMatJoint3x4::scale(scale) * CMatJoint3x4::translate(-scaleCenter);
}

CONST_WIN32 CVec3D CMatJoint3x4::translatePart() const
{
	return Col(3);
}

CONST_WIN32 CMat3D CMatJoint3x4::rotatePart() const
{
	return CMat3DPart();
}


void CMatJoint3x4::transpose3()
{
	///\todo SSE.
	Swap(mat[0*4+1], mat[1*4+0]);
	Swap(mat[0*4+2], mat[2*4+0]);
	Swap(mat[1*4+2], mat[2*4+1]);
}

CMatJoint3x4 CMatJoint3x4::transposed3() const
{
	CMatJoint3x4 copy = *this;
	copy.transpose3();
	return copy;
}

bool CMatJoint3x4::inverse(float epsilon)
{
	///\todo SSE.
#ifdef TEST_FOR_CORRECTNESS
	CMatJoint3x4 orig = *this;
#endif
	CMat4D temp(*this); ///\todo It is possible optimize to avoid copying here by writing the inverse function specifically for float3x4.
	bool success = temp.inverseSelf();
	*this = temp.mat3x4Part();
#ifdef TEST_FOR_CORRECTNESS
	SMF_ASSERT(!success || (orig * *this).isIdentity(1e-1f));
#endif
	return success;
}


bool CMatJoint3x4::inverseTranspose()
{
	bool success = inverse();
	transpose3();
	// float3x4 cannot represent the translation element as the fourth row after transposing.
	// Since inverse transposes are used mainly to transform direction vectors, we can discard the translation component.
	setTranslatePart(0,0,0);
	return success;
}

CMatJoint3x4 CMatJoint3x4::inverseTransposed() const
{
	CMatJoint3x4 copy = *this;
	copy.transpose3();
	copy.inverse();
	// float3x4 cannot represent the translation element as the fourth row after transposing.
	// Since inverse transposes are used mainly to transform direction vectors, we can discard the translation component.
	copy.setTranslatePart(0,0,0);
	return copy;
}

float CMatJoint3x4::trace() const
{
	SMF_ASSERT(isFinite());
	return mat[0*4+0] + mat[1*4+1] + mat[2*4+2];
}

CMat4D::CMat4D(const CMatJoint3x4 &m)
{

	set(m.Row(0).x, m.Row(0).y, m.Row(0).z, m.Row(0).w,
		m.Row(1).x, m.Row(1).y, m.Row(1).z, m.Row(1).w,
		m.Row(2).x, m.Row(2).y, m.Row(2).z, m.Row(2).w,
		      0.f,       0.f,       0.f,     1.f);

}




} //end MATH
}//end SMF
