/*
  SMF -  Rasputtim Game Fabric  (https://sourceforge.net/projects/sgfabric/)
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

#ifndef _SMF__MATH_MATRIX_H__
#define _SMF__MATH_MATRIX_H__
#include "../SMF_Config.h"
#include "SMF_MathDefs.h"
#include "SMF_Mat3D.h"
#include "../exceptions/all.h"
//#include "../configuration/SMF_Configuration.h"
#include "SMF_Vector.h"
#include "SMF_Simd.h"
#include "SMF_Quaternion.h"
//#include "math/SMF_JointTransform.h"

namespace SMF{
namespace MATH {
/*
===============================================================================

  Matrix classes, all matrices are row-major except CMat3D

===============================================================================
*/
template<typename Matrix>
void setMatrixRotatePart(Matrix &m, const CQuaternion &q)
{
	// See e.g. http://www.geometrictools.com/Documentation/LinearAlgebraicQuaternions.pdf .

	assume(q.isNormalized());
	const float x = q.x; const float y = q.y; const float z = q.z; const float w = q.w;
	m[0][0] = 1 - 2*(y*y + z*z); m[0][1] =     2*(x*y - z*w); m[0][2] =     2*(x*z + y*w);
	m[1][0] =     2*(x*y + z*w); m[1][1] = 1 - 2*(x*x + z*z); m[1][2] =     2*(y*z - x*w);
	m[2][0] =     2*(x*z - y*w); m[2][1] =     2*(y*z + x*w); m[2][2] = 1 - 2*(x*x + y*y);
}

/** Sets the top-left 3x3 area of the matrix to the rotation matrix about the X-axis. Elements
	outside the top-left 3x3 area are ignored. This matrix rotates counterclockwise if multiplied
	in the order M*v, and clockwise if rotated in the order v*M.
	\param m The matrix to store the result.
	\param angle The rotation angle in radians.
*/
template<typename Matrix>
void set3x3PartRotateX(Matrix &m, float angle)
{
	/*
	 ³   1   0   0 ³
	 ³   0  cz -sz ³
	 ³   0  sz  cz ³
	*/

	const float cosz = cos(angle);
	const float sinz = sin(angle);

	m[0][0] = 1.f; m[0][1] =  0.f; m[0][2] =   0.f;
	m[1][0] = 0.f; m[1][1] = cosz; m[1][2] = -sinz;
	m[2][0] = 0.f; m[2][1] = sinz; m[2][2] =  cosz;
}

/** Sets the top-left 3x3 area of the matrix to the rotation matrix about the Y-axis. Elements
	outside the top-left 3x3 area are ignored. This matrix rotates counterclockwise if multiplied
	in the order M*v, and clockwise if rotated in the order v*M.
	\param m The matrix to store the result.
	\param angle The rotation angle in radians.
*/
template<typename Matrix>
void set3x3PartRotateY(Matrix &m, float angle)
{
	/*
	 ³  cz   0  sz ³
	 ³   0   1   0 ³
	 ³ -sz   0  cz ³
	*/

	const float cosz = cos(angle);
	const float sinz = sin(angle);

	m[0][0] =  cosz; m[0][1] = 0.f; m[0][2] = sinz;
	m[1][0] =   0.f; m[1][1] = 1.f; m[1][2] =  0.f;
	m[2][0] = -sinz; m[2][1] = 0.f; m[2][2] = cosz;
}

/** Sets the top-left 3x3 area of the matrix to the rotation matrix about the Z-axis. Elements
	outside the top-left 3x3 area are ignored. This matrix rotates counterclockwise if multiplied
	in the order M*v, and clockwise if rotated in the order v*M.
	@param m The matrix to store the result.
	@param angle The rotation angle in radians. */
template<typename Matrix>
void set3x3PartRotateZ(Matrix &m, float angle)
{
	/*
	 ³  cz -sz   0 ³
	 ³  sz  cz   0 ³
	 ³   0   0   1 ³
	*/

	const float cosz = cos(angle);
	const float sinz = sin(angle);

	m[0][0] = cosz; m[0][1] = -sinz; m[0][2] = 0.f;
	m[1][0] = sinz; m[1][1] =  cosz; m[1][2] = 0.f;
	m[2][0] =  0.f; m[2][1] =   0.f; m[2][2] = 1.f;
}


#define MATRIX_INVERSE_EPSILON		1e-14
#define MATRIX_EPSILON				1e-6

class CEulerAngles;
class CQuaternion;
class CCompQuaternion;
class CRotation;
class CMat4D;
class CMatJoint3x4;

/** para conseguir o melhor resultado de performance a biblioteca deve ser compilada com a opção /Ox (msvc) full otimization

**/
void PII_inverse_4x4(float* mat);
void PIII_inverse_4x4(float* src);

void PII_inverse_6x6(float* mat);
void PIII_InverseG_6x6(float *src);
void PIII_InverseS_6x6(float *src);



/**
 * \class CMat2D
 *
 * \ingroup SMF_Math
 *
 * \if pt_br
 * \brief Matriz 2x2
 * \elseif us_en
 * \brief 2x2 matrix
 * \endif
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
class SMF_API CMat2D {
public:
					CMat2D();
					explicit CMat2D( const CVec2D &x, const CVec2D &y );
					explicit CMat2D( const float xx, const float xy, const float yx, const float yy );
					explicit CMat2D( const float src[ 2 ][ 2 ] );

	const CVec2D &	operator[]( int index ) const;
	CVec2D &		operator[]( int index );
	CMat2D			operator-() const;
	CMat2D			operator*( const float a ) const;
	CVec2D			operator*( const CVec2D &vec ) const;
	CMat2D			operator*( const CMat2D &a ) const;
	CMat2D			operator+( const CMat2D &a ) const;
	CMat2D			operator-( const CMat2D &a ) const;
	CMat2D &		operator*=( const float a );
	CMat2D &		operator*=( const CMat2D &a );
	CMat2D &		operator+=( const CMat2D &a );
	CMat2D &		operator-=( const CMat2D &a );

	friend CMat2D	operator*( const float a, const CMat2D &mat );
	friend CVec2D	operator*( const CVec2D &vec, const CMat2D &mat );
	friend CVec2D &	operator*=( CVec2D &vec, const CMat2D &mat );
	CVec2D mul(const CVec2D &rhs) const { return *this * rhs; };
	CVec2D MulPos(const CVec2D &rhs) const { return mul(rhs); }
	CVec2D MulDir(const CVec2D &rhs) const { return mul(rhs); }
		/// Returns the given row. [noscript]
	/** \param row The zero-based index [0, 1] of the row to get. */
	CVec2D &Row(int row);
	const CVec2D &Row(int row) const;

	CVec2D &Row2(int row) { return Row(row); }
	const CVec2D &Row2(int row) const { return Row(row); }

	/// Returns the given column.
	/** \param col The zero-based index [0, 1] of the column to get. */
	CONST_WIN32 CVec2D Col(int col) const;
	CONST_WIN32 CVec2D Col2(int col) const { return Col(col); }

	/// Returns the main diagonal.
	/** The main diagonal consists of the elements at m[0][0], m[1][1]. */
	CONST_WIN32 CVec2D diagonal() const;

	bool			compare( const CMat2D &a ) const;						// exact compare, no epsilon
	bool			compare( const CMat2D &a, const float epsilon ) const;	// compare with epsilon
	bool			operator==( const CMat2D &a ) const;					// exact compare, no epsilon
	bool			operator!=( const CMat2D &a ) const;					// exact compare, no epsilon

	void			toZero();
	void			toIdentity();
	bool			isIdentity( const float epsilon = MATRIX_EPSILON ) const;
	bool			isSymmetric( const float epsilon = MATRIX_EPSILON ) const;
	bool			isDiagonal( const float epsilon = MATRIX_EPSILON ) const;

	float			trace() const;
	float			determinant() const;
    /// returns transpose
	CMat2D			transpose() const;
	CMat2D &		transposeSelf();
	/**
	\brief interte a matriz  ( m * m.inverse() = identity)
	\return retorna a matriz inversa ou throw exception, caso não haja inversa
	**/
	CMat2D		inverse() const throw(Exception::CMathException) ;		// returns the inverse ( m * m.inverse() = identity )
	/**
	\brief calcula a matriz inversa ( m * m.Inverse() = identity )
	\return retorna falso caso não haja inversa e verdadeiro caso haja
	**/
	bool			inverseSelf();		// returns false if determinant is zero
	/**
	\brief interte a matriz  ( m * m.inverse() = identity)
	\return retorna a matriz inversa ou throw exception, caso não haja inversa
	**/
	CMat2D			inverseFast() const throw(Exception::CMathException);	// returns the inverse ( m * m.inverse() = identity )
	/**
	\brief calcula a matriz inversa ( m * m.Inverse() = identity )
	\return retorna falso caso não haja inversa e verdadeiro caso haja
	**/
	bool			inverseFastSelf();	// returns false if determinant is zero

	int				getDimension() const;

	const float *	toFloatPtr() const;
	float *			toFloatPtr();
	const char *	toString( int precision = 2 ) const;

private:
	CVec2D			mat[ 2 ];
};

extern CMat2D mat2_zero;
extern CMat2D mat2_identity;
#define mat2_default	mat2_identity

SMF_INLINE_FORCED CMat2D::CMat2D() {
}

SMF_INLINE_FORCED CMat2D::CMat2D( const CVec2D &x, const CVec2D &y ) {
	mat[ 0 ].x = x.x; mat[ 0 ].y = x.y;
	mat[ 1 ].x = y.x; mat[ 1 ].y = y.y;
}

SMF_INLINE_FORCED CMat2D::CMat2D( const float xx, const float xy, const float yx, const float yy ) {
	mat[ 0 ].x = xx; mat[ 0 ].y = xy;
	mat[ 1 ].x = yx; mat[ 1 ].y = yy;
}

SMF_INLINE_FORCED CMat2D::CMat2D( const float src[ 2 ][ 2 ] ) {
	memcpy( mat, src, 2 * 2 * sizeof( float ) );
}

SMF_INLINE_FORCED const CVec2D &CMat2D::operator[]( int index ) const {
	SMF_ASSERT( ( index >= 0 ) && ( index < 2 ) );
	if (( index < 0 ) || ( index >= 2 )) {
		string error("invalid index");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}
	return mat[ index ];
}

SMF_INLINE_FORCED CVec2D &CMat2D::operator[]( int index ) {
	SMF_ASSERT( ( index >= 0 ) && ( index < 2 ) );
	if (( index < 0 ) || ( index >= 2 )) {
		string error("invalid index");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}
	return mat[ index ];
}

SMF_INLINE_FORCED CMat2D CMat2D::operator-() const {
	return CMat2D(	-mat[0][0], -mat[0][1],
					-mat[1][0], -mat[1][1] );
}

SMF_INLINE_FORCED CVec2D CMat2D::operator*( const CVec2D &vec ) const {
	return CVec2D(
		mat[ 0 ].x * vec.x + mat[ 0 ].y * vec.y,
		mat[ 1 ].x * vec.x + mat[ 1 ].y * vec.y );
}

SMF_INLINE_FORCED CMat2D CMat2D::operator*( const CMat2D &a ) const {
	return CMat2D(
		mat[0].x * a[0].x + mat[0].y * a[1].x,
		mat[0].x * a[0].y + mat[0].y * a[1].y,
		mat[1].x * a[0].x + mat[1].y * a[1].x,
		mat[1].x * a[0].y + mat[1].y * a[1].y );
}

SMF_INLINE_FORCED CMat2D CMat2D::operator*( const float a ) const {
	return CMat2D(
		mat[0].x * a, mat[0].y * a,
		mat[1].x * a, mat[1].y * a );
}

SMF_INLINE_FORCED CMat2D CMat2D::operator+( const CMat2D &a ) const {
	return CMat2D(
		mat[0].x + a[0].x, mat[0].y + a[0].y,
		mat[1].x + a[1].x, mat[1].y + a[1].y );
}

SMF_INLINE_FORCED CMat2D CMat2D::operator-( const CMat2D &a ) const {
	return CMat2D(
		mat[0].x - a[0].x, mat[0].y - a[0].y,
		mat[1].x - a[1].x, mat[1].y - a[1].y );
}

SMF_INLINE_FORCED CMat2D &CMat2D::operator*=( const float a ) {
	mat[0].x *= a; mat[0].y *= a;
	mat[1].x *= a; mat[1].y *= a;

    return *this;
}

SMF_INLINE_FORCED CMat2D &CMat2D::operator*=( const CMat2D &a ) {
	float x, y;
	x = mat[0].x; y = mat[0].y;
	mat[0].x = x * a[0].x + y * a[1].x;
	mat[0].y = x * a[0].y + y * a[1].y;
	x = mat[1].x; y = mat[1].y;
	mat[1].x = x * a[0].x + y * a[1].x;
	mat[1].y = x * a[0].y + y * a[1].y;
	return *this;
}

SMF_INLINE_FORCED CMat2D &CMat2D::operator+=( const CMat2D &a ) {
	mat[0].x += a[0].x; mat[0].y += a[0].y;
	mat[1].x += a[1].x; mat[1].y += a[1].y;

    return *this;
}

SMF_INLINE_FORCED CMat2D &CMat2D::operator-=( const CMat2D &a ) {
	mat[0].x -= a[0].x; mat[0].y -= a[0].y;
	mat[1].x -= a[1].x; mat[1].y -= a[1].y;

    return *this;
}

SMF_INLINE_FORCED CVec2D operator*( const CVec2D &vec, const CMat2D &mat ) {
	return mat * vec;
}

SMF_INLINE_FORCED CMat2D operator*( const float a, CMat2D const &mat ) {
	return mat * a;
}

SMF_INLINE_FORCED CVec2D &operator*=( CVec2D &vec, const CMat2D &mat ) {
	vec = mat * vec;
	return vec;
}
SMF_INLINE_FORCED CVec2D &CMat2D::Row(int row)
{
	SMF_ASSERT(row >= 0);
	SMF_ASSERT(row < 2);
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (row < 0 || row >= 3)
		row = 0; // Benign failure, just give the first row.
#endif
	return reinterpret_cast<CVec2D &>(mat[row]);
}

SMF_INLINE_FORCED const CVec2D &CMat2D::Row(int row) const
{
	SMF_ASSERT(row >= 0);
	SMF_ASSERT(row < 2);
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (row < 0 || row >= 2)
		row = 0; // Benign failure, just give the first row.
#endif
	return CVec2D(mat[0][row], mat[1][row]);
}

SMF_INLINE_FORCED CONST_WIN32 CVec2D CMat2D::Col(int col) const
{
	SMF_ASSERT(col >= 0);
	SMF_ASSERT(col < 2);
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (col < 0 || col >= 2)
		return CVec2D::nan;
#endif
	return reinterpret_cast<const CVec2D &>(mat[col]);

}

SMF_INLINE_FORCED CONST_WIN32 CVec2D CMat2D::diagonal() const
{
	return CVec2D(mat[0][0], mat[1][1]);
}

SMF_INLINE_FORCED bool CMat2D::compare( const CMat2D &a ) const {
	if ( mat[0].compare( a[0] ) &&
		mat[1].compare( a[1] ) ) {
		return true;
	}
	return false;
}

SMF_INLINE_FORCED bool CMat2D::compare( const CMat2D &a, const float epsilon ) const {
	if ( mat[0].compare( a[0], epsilon ) &&
		mat[1].compare( a[1], epsilon ) ) {
		return true;
	}
	return false;
}

SMF_INLINE_FORCED bool CMat2D::operator==( const CMat2D &a ) const {
	return compare( a );
}

SMF_INLINE_FORCED bool CMat2D::operator!=( const CMat2D &a ) const {
	return !compare( a );
}

SMF_INLINE_FORCED void CMat2D::toZero() {
	mat[0].toZero();
	mat[1].toZero();
}

SMF_INLINE_FORCED void CMat2D::toIdentity() {
	*this = mat2_identity;
}

SMF_INLINE_FORCED bool CMat2D::isIdentity( const float epsilon ) const {
	return compare( mat2_identity, epsilon );
}

SMF_INLINE_FORCED bool CMat2D::isSymmetric( const float epsilon ) const {
	return ( CMath::fabs( mat[0][1] - mat[1][0] ) < epsilon );
}

SMF_INLINE_FORCED bool CMat2D::isDiagonal( const float epsilon ) const {
	if ( CMath::fabs( mat[0][1] ) > epsilon ||
		CMath::fabs( mat[1][0] ) > epsilon ) {
		return false;
	}
	return true;
}

SMF_INLINE_FORCED float CMat2D::trace() const {
	return ( mat[0][0] + mat[1][1] );
}

SMF_INLINE_FORCED float CMat2D::determinant() const {
	return mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];
}

SMF_INLINE_FORCED CMat2D CMat2D::transpose() const {
	return CMat2D(	mat[0][0], mat[1][0],
					mat[0][1], mat[1][1] );
}

SMF_INLINE_FORCED CMat2D &CMat2D::transposeSelf() {
	float tmp;

	tmp = mat[0][1];
	mat[0][1] = mat[1][0];
	mat[1][0] = tmp;

	return *this;
}

SMF_INLINE_FORCED CMat2D CMat2D::inverse() const throw(Exception::CMathException){
	CMat2D invMat;

	invMat = *this;
	int r = invMat.inverseSelf();
	if (r ==  0) {
		string error("matriz have no inverse");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}
	return invMat;
}

SMF_INLINE_FORCED CMat2D CMat2D::inverseFast() const throw(Exception::CMathException){
	CMat2D invMat;

	invMat = *this;
	int r = invMat.inverseFastSelf();
	if (r ==  0) {
		string error("matriz have no inverse");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}
	return invMat;
}

SMF_INLINE_FORCED int CMat2D::getDimension() const {
	return 4;
}

SMF_INLINE_FORCED const float *CMat2D::toFloatPtr() const {
	return mat[0].toFloatPtr();
}

SMF_INLINE_FORCED float *CMat2D::toFloatPtr() {
	return mat[0].toFloatPtr();
}



/**
 * \class CMat4D
 *
 * \ingroup SMF_Math
 *
 * \if pt_br
 * \brief Matriz 4x4
 *
 * \elseif us_en
 * \brief 4x4 matrix
 *
 * \endif
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
class SMF_API CMat4D {
public:
	friend class CSIMD_Generic;
	CMat4D();
	explicit CMat4D( const CVec4D &x, const CVec4D &y, const CVec4D &z, const CVec4D &w );
	explicit CMat4D(const float xx, const float xy, const float xz, const float xw,
									const float yx, const float yy, const float yz, const float yw,
									const float zx, const float zy, const float zz, const float zw,
									const float wx, const float wy, const float wz, const float ww );
	explicit CMat4D( const CMat3D &rotation, const CVec3D &translation );
	explicit CMat4D( const float src[ 4 ][ 4 ] );

	/// Constructs this Mat4D from the given quaternion.
	explicit CMat4D(const CQuaternion  &orientation);

	/**
	\brief Constructs this float4x4 to represent the same transformation as the given float3x4.
	 The last row is set to [0 0 0 1].
	 */
	CMat4D(const CMatJoint3x4 &other);


	/**
	\brief Constructs this CMat4D from the given quaternion and translation.
	 Logically, the translation occurs after the rotation has been performed.
	 **/
	CMat4D(const CQuaternion  &orientation, const CVec3D &translation);

	const CVec4D &	operator[]( int index ) const throw(Exception::CMathException);
	CVec4D &		operator[]( int index ) throw(Exception::CMathException);
	CMat4D			operator*( const float a ) const;
	/// Performs standard matrix multiplication.
	CVec4D			operator*( const CVec4D &vec ) const;
	CVec3D			operator*( const CVec3D &vec ) const;
	CMat4D			operator*( const CMat4D &a ) const;
	/// Treats the CMat3D as a 4-by-4 matrix with the last row and column as identity, and multiplies the two matrices.
	CMat4D operator *(const CMat3D &rhs) const;

	/// Treats the CMatJoint3x4 as a 4-by-4 matrix with the last row as identity, and multiplies the two matrices.
	CMat4D operator *(const CMatJoint3x4 &rhs) const;


	/// Converts the quaternion to a CMat4D and multiplies the two matrices together.
	CMat4D operator *(const CQuaternion &rhs) const;

	CMat4D			operator+( const CMat4D &a ) const;
	CMat4D			operator-( const CMat4D &a ) const;
	CMat4D &		operator*=( const float a );
	CMat4D &		operator*=( const CMat4D &a );
	CMat4D &		operator+=( const CMat4D &a );
	CMat4D &		operator-=( const CMat4D &a );

	friend CMat4D	operator*( const float a, const CMat4D &mat );
	friend CVec4D	operator*( const CVec4D &vec, const CMat4D &mat );
	friend CVec3D	operator*( const CVec3D &vec, const CMat4D &mat );
	friend CVec4D &	operator*=( CVec4D &vec, const CMat4D &mat );
	friend CVec3D &	operator*=( CVec3D &vec, const CMat4D &mat );

	//CONST_WIN32 float at(int row, int col) const;

	/**
	\brief Returns the given row. [noscript]
	\param row The zero-based index [0, 3] of the row to get.
	*/
	CVec4D &Row(int row);
	const CVec4D &Row(int row) const;
	/// Returns the three first entries of the given row. [similarOverload: Row] [hideIndex]
	CVec3D &Row3(int row); ///< [noscript]
	const CVec3D &Row3(int row) const;

	/**
	\brief Returns the given column.
	 \param col The zero-based index [0, 3] of the column to get.
	 */
	CONST_WIN32 CVec4D Col(int col) const;
	/// Returns the three first entries of the given column. [similarOverload: Column] [hideIndex]
	CONST_WIN32 CVec3D Col3(int col) const;

	/**
	\brief Returns the main diagonal.
	 The main diagonal consists of the elements at m[0][0], m[1][1], m[2][2] and m[3][3]. */
	CONST_WIN32 CVec4D diagonal() const;
	/// Returns the three first entries of the main diagonal. [similarOverload: MainDiagonal] [hideIndex]
	CONST_WIN32 CVec3D diagonal3() const;

	/// Sets all values of this matrix.
	void set(float _00, float _01, float _02, float _03,
			 float _10, float _11, float _12, float _13,
			 float _20, float _21, float _22, float _23,
			 float _30, float _31, float _32, float _33);



	bool			compare( const CMat4D &a ) const;						// exact compare, no epsilon
	bool			compare( const CMat4D &a, const float epsilon ) const;	// compare with epsilon
	bool			operator==( const CMat4D &a ) const;					// exact compare, no epsilon
	bool			operator!=( const CMat4D &a ) const;					// exact compare, no epsilon

	void			toZero();
	void			toIdentity();
	/**
	\brief Tests if this is the identity matrix.
	\return Returns true if this matrix is the identity matrix, up to the given epsilon.
	**/
	bool			isIdentity( const float epsilon = MATRIX_EPSILON ) const;
	/**
	\brief Tests if this matrix is symmetric (M == M^T).
	 The test compares the elements for equality, up to the given epsilon. A matrix is symmetric if it is its own transpose.
	 **/
	bool			isSymmetric( const float epsilon = MATRIX_EPSILON ) const;
	bool			isDiagonal( const float epsilon = MATRIX_EPSILON ) const;
	bool			isRotated() const;

	/**
	\brief Tests if this matrix does not contain any NaNs or infs.
	\return Returns true if the entries of thisCMat4D are all finite, and do not contain NaN or infs.
	**/
	bool isFinite() const;


	/**
	\brief Tests if this matrix is in lower triangular form.
	\return Returns true if this matrix is in lower triangular form, up to the given epsilon.
	**/
	bool isLowerTriangular(float epsilon =CMath::EPSILON_SuperLow) const;

	/**
	\brief Tests if this matrix is in upper triangular form.
	\return Returns true if this matrix is in upper triangular form, up to the given epsilon.
	**/
	bool isUpperTriangular(float epsilon =CMath::EPSILON_SuperLow) const;

	/**
	\brief Tests if this matrix has an inverse.
	\return Returns true if this matrix can be inverted, up to the given epsilon.
	**/
	bool isInvertible(float epsilon =CMath::EPSILON_SuperLow) const;


	/**
	\brief Tests if this matrix is skew-symmetric (M == -M^T).
	The test compares the floating point elements of this matrix up to the given epsilon. A matrix M is skew-symmetric
	the identity M=-M^T holds.
	**/
	bool isSkewSymmetric(float epsilon =CMath::EPSILON_SuperLow) const;

	/**
	\brief Tests if this matrix is an idempotent matrix.
	An idempotent matrix is one for which the equality M*M=M holds. Projection matrices are commonly idempotent.
	**/
	bool isIdempotent(float epsilon =CMath::EPSILON_SuperLow) const;

	/**
	\breif Returns true if this matrix does not perform any scaling.
	A matrix does not do any scaling if the column vectors of this
	matrix are normalized in length, compared to the given epsilon. Note that this matrix may still perform
	reflection, i.e. it has a -1 scale along some axis.
	\note This function only examines the upper 3-by-3 part of this matrix.
	\note This function assumes that this matrix does not contain projection (the fourth row of this matrix is [0 0 0 1]).
	**/
	bool hasUnitaryScale(float epsilon =CMath::EPSILON_SuperLow) const;

	/**
	\brief Returns true if this matrix performs a reflection along some plane.
	 In 3D space, an even number of reflections corresponds to a rotation about some axis, so a matrix consisting of
	an odd number of consecutive mirror operations can only reflect about one axis. A matrix that contains reflection reverses
	the handedness of the coordinate system. This function tests if this matrix
	does perform mirroring. This occurs iff this matrix has a negative determinant. */
	bool hasNegativeScale() const;

	/**
	\brief Returns true if this matrix contains only uniform scaling, compared to the given epsilon.
	\note If the matrix does not really do any scaling, this function returns true (scaling uniformly by a factor of 1).
	\note This function only examines the upper 3-by-3 part of this matrix.
	\note This function assumes that this matrix does not contain projection (the fourth row of this matrix is [0 0 0 1]).
	**/
	bool hasUniformScale(float epsilon =CMath::EPSILON_SuperLow) const;

	/// Returns true if the row vectors of 3x3 top-left submatrix are all perpendicular to each other.
	bool isRowOrthogonal3(float epsilon =CMath::EPSILON_SuperLow) const;

	/// Returns true if the column vectors of 3x3 top-left submatrix are all perpendicular to each other.
	bool isColOrthogonal3(float epsilon =CMath::EPSILON_SuperLow) const;

	/**
	\brief Returns true if the column and row vectors of the 3x3 top-left submatrix form an orthonormal set.
	\note In math terms, there does not exist such a thing as 'orthonormal matrix'. In math terms,
	a matrix is orthogonal iff its column and row vectors are orthogonal *unit* vectors.
	In the terms of this library however, a matrix is orthogonal iff its column and row vectors are orthogonal (no need to be unitary),
	and a matrix is orthonormal if the column and row vectors are orthonormal.
	**/
	bool IsOrthonormal3(float epsilon =CMath::EPSILON_SuperLow) const;



	void			projectVector( const CVec4D &src, CVec4D &dst ) const;
	void			unprojectVector( const CVec4D &src, CVec4D &dst ) const;

	float			trace() const;
	/**
	\brief Computes the determinant of this matrix.
	If the determinant is nonzero, this matrix is invertible.
	If the determinant is negative, this matrix performs reflection about some axis.
	From http://msdn.microsoft.com/en-us/library/bb204853(VS.85).aspx :
	"If the determinant is positive, the basis is said to be "positively" oriented (or right-handed).
	If the determinant is negative, the basis is said to be "negatively" oriented (or left-handed)."
	**/
	float			determinant() const;
	/// Computes the determinant of the upper-left 3x3 submatrix of this matrix.
	float Determinant3() const;

	/// Returns the upper-left 3-by-3 part.
	CONST_WIN32 CMat3D CMat3DPart() const;

	/**
	\brief Returns the upper-left 3-by-4 part. [noscript]
	\note The CMatJoint3x4 and CMat4D are bit-compatible, so this function simply casts.
	**/

//	CMatJoint3x4 &CMatJoint3x4Part();
//	const CMatJoint3x4 &CMatJoint3x4Part() const;

	/// returns transpose
	CMat4D		transpose() const;
	CMat4D &	transposeSelf();
	/**
	\brief interte a matriz  ( m * m.inverse() = identity)
	\return retorna a matriz inversa ou throw exception, caso não haja inversa
	**/
	CMat4D		inverse() const throw(Exception::CMathException) ;		// returns the inverse ( m * m.inverse() = identity )


	/**
	\brief calcula a matriz inversa ( m * m.inverse() = identity )
	\return retorna falso caso não haja inversa e verdadeiro caso haja
	\note http://download.intel.com/design/PentiumIII/sml/24504301.pdf
	**/
	bool			inverseSelf();		// returns false if determinant is zero
	/**
	\brief interte a matriz  ( m * m.inverse() = identity)
	\return retorna a matriz inversa ou throw exception, caso não haja inversa
	**/
	CMat4D			inverseFast() const throw(Exception::CMathException);	// returns the inverse ( m * m.inverse() = identity )
	/**
	\brief calcula a matriz inversa ( m * m.inverse() = identity )
	\return retorna falso caso não haja inversa e verdadeiro caso haja
	**/
	bool		inverseFastSelf();	// returns false if determinant is zero
	CMat4D		transposeMultiply( const CMat4D &b ) const;

	int				getDimension() const;

	const float *	toFloatPtr() const;
	float *			toFloatPtr();
	const char *	toString( int precision = 2 ) const;

	/// Transforms the given point vector by this matrix M , i.e. returns M * (x, y, z, 1).
	CVec3D transformPos(const CVec3D &pointVector) const;
	CVec3D transformPos(float x, float y, float z) const;

	/// Transforms the given direction vector by this matrix M , i.e. returns M * (x, y, z, 0).
	CVec3D transformDir(const CVec3D &directionVector) const;
	CVec3D transformDir(float x, float y, float z) const;

	/**
	\brief Transforms the given 4-vector by this matrix M, i.e. returns M * (x, y, z, w).
	Does not perform a perspective divide afterwards, so remember to divide by w afterwards
	at some point, if this matrix contained a projection.
	**/
	CVec4D transform(const CVec4D &vector) const;

	CMat4D mul(const CMat3D &rhs) const { return *this * rhs; };
	CMat4D mul(const CMatJoint3x4 &rhs) const { return *this * rhs; };
	CMat4D mul(const CMat4D &rhs) const { return *this * rhs; };
	CMat4D mul(const CQuaternion &rhs) const { return *this * rhs; };
	CVec3D MulPos(const CVec3D &pointVector) const { return this->transformPos(pointVector); };
	CVec3D MulDir(const CVec3D &directionVector) const { return this->transformDir(directionVector); };
	CVec4D mul(const CVec4D &vector) const { return *this * vector; };
	/**
	\brief Returns the scale components of this matrix.
	This function decomposes this matrix M into a form M = M' * S, where M' has unitary column vectors and S is a diagonal matrix.
	\return extractScale returns the diagonal entries of S, i.e. the scale of the columns of this matrix . If this matrix
	represents a local->world space transformation for an object, then this scale represents a 'local scale', i.e.
	scaling that is performed before translating and rotating the object from its local coordinate system to its world
	position.
	\note This function assumes that this matrix does not contain projection (the fourth row of this matrix is [0 0 0 1]).
	\note This function does not detect and return reflection (-1 scale along some axis).
	*/
	CVec3D extractScale() const;

	/**
	\brief Returns the translation part.
	The translation part is stored in the fourth column of this matrix.
	This is equivalent to decomposing this matrix in the form M = T * M', i.e. this translation is applied last,
	after applying rotation and scale. If this matrix represents a local->world space transformation for an object,
	then this gives the world space position of the object.
	\note This function assumes that this matrix does not contain projection (the fourth row of this matrix is [0 0 0 1]).
	**/
	CONST_WIN32 CVec3D translatePart() const;

	/// Returns the top-left 3x3 part of this matrix. This stores the rotation part of this matrix (if this matrix represents a rotation).
	CONST_WIN32 CMat3D rotatePart() const;

	/// Returns the upper-left 3-by-3 part.
	CONST_WIN32 CMat3D mat3Part() const;

	/** Returns the upper-left 3-by-4 part. [noscript]
	\note The float3x4 and float4x4 are bit-compatible, so this function simply casts.
	**/
	CMatJoint3x4 & mat3x4Part();
	const CMatJoint3x4 & mat3x4Part() const;

	/**
	\brief Sets the translation part of this matrix.
	 This function sets the translation part of this matrix. These are the three first elements of the fourth column.
	All other entries are left untouched.
	**/
	void setTranslatePart(float tx, float ty, float tz);
	void setTranslatePart(const CVec3D &offset);

	/**
	\brief Sets the 3-by-3 part of this matrix to perform rotation about the positive X axis which passes through
	the origin. Leaves all other entries of this matrix untouched. [similarOverload: setRotatePart] [hideIndex]
	**/
	void setRotatePartX(float angleRadians);
	/**
	\brief Sets the 3-by-3 part of this matrix to perform rotation about the positive Y axis. Leaves all other
	entries untouched. [similarOverload: setRotatePart] [hideIndex]
	**/
	void setRotatePartY(float angleRadians);
	/**
	\brief Sets the 3-by-3 part of this matrix to perform rotation about the positive Z axis. Leaves all other
	entries untouched. [similarOverload: setRotatePart] [hideIndex]
	**/
	void setRotatePartZ(float angleRadians);

	/**
	\brief Sets the 3-by-3 part of this matrix to perform rotation about the given axis and angle (in radians). Leaves all other
	entries of this matrix untouched. [indexTitle: setRotatePart/X/Y/Z]
	**/
	void setRotatePart(const CVec3D &axisDirection, float angleRadians);
	/**
	\brief Sets the 3-by-3 part of this matrix to perform the rotation expressed by the given quaternion.
	Leaves all other entries of this matrix untouched.
	**/
	void setRotatePart(const CQuaternion &orientation);

	/**
	\brief  Sets the 3-by-3 part of this matrix.
	\note This is a convenience function which calls set3x3Part.
	\note This function erases the previous top-left 3x3 part of this matrix (any previous rotation, scaling and shearing, etc.). Translation is unaffected.
	**/
	void setRotatePart(const CMat3D &rotation) { set3x3Part(rotation); }

	void set3x3Part(const CMat3D &rotation);
	void set3x4Part(const CMatJoint3x4 &rotateTranslate);

	/**
	\brief Sets the three first elements of the given column. The fourth element is left unchanged.
	\param column The index of the column to set, in the range [0-3].
	\param data A pointer to an array of 3 floats that contain the new x, y and z values for the column.
	**/
	void setCol3(int column, const float *data);
	void setCol3(int column, const CVec3D &columnVector);
	void setCol3(int column, float m_0c, float m_1c, float m_2c);
	/**
	\brief Sets the values of the given row.
	\param row The index of the row to set, in the range [0-3].
	\param data A pointer to an array of 4 floats that contain the new x, y, z and w values for the row.
	**/
	void setRow(int row, const float *data);
	void setRow(int row, const CVec3D &rowVector, float m_r3);
	void setRow(int row, const CVec4D &rowVector);
	void setRow(int row, float m_r0, float m_r1, float m_r2, float m_r3);
	/**
	\brief Sets the values of the given column.
	\param column The index of the column to set, in the range [0-3].
	\param data A pointer to an array of 4 floats that contain the new x, y, z and w values for the column.
	**/
	void setCol(int column, const float *data);
	void setCol(int column, const CVec3D &columnVector, float m_3c);
	void setCol(int column, const CVec4D &columnVector);
	void setCol(int column, float m_0c, float m_1c, float m_2c, float m_3c);
	/**
	\brief Returns true if this matrix is seen to contain a "projective" part, i.e. whether the last row of this
	 matrix differs from [0 0 0 1].
	 **/
	bool ContainsProjection(float epsilon = 1e-3f) const;

	/// Returns a string representation of form "(m00, m01, m02, m03; m10, m11, m12, m13; ... )".
	std::string ToString() const;

	std::string ToString2() const;



public:
	CVec4D			mat[ 4 ];

};

extern CMat4D mat4_zero;
extern CMat4D mat4_identity;
#define mat4_default	mat4_identity




SMF_INLINE_FORCED CMat4D::CMat4D() {
}

SMF_INLINE_FORCED CMat4D::CMat4D( const CVec4D &x, const CVec4D &y, const CVec4D &z, const CVec4D &w ) {
	mat[ 0 ] = x;
	mat[ 1 ] = y;
	mat[ 2 ] = z;
	mat[ 3 ] = w;
}

SMF_INLINE_FORCED CMat4D::CMat4D( const float xx, const float xy, const float xz, const float xw,
							const float yx, const float yy, const float yz, const float yw,
							const float zx, const float zy, const float zz, const float zw,
							const float wx, const float wy, const float wz, const float ww ) {
	mat[0][0] = xx; mat[0][1] = xy; mat[0][2] = xz; mat[0][3] = xw;
	mat[1][0] = yx; mat[1][1] = yy; mat[1][2] = yz; mat[1][3] = yw;
	mat[2][0] = zx; mat[2][1] = zy; mat[2][2] = zz; mat[2][3] = zw;
	mat[3][0] = wx; mat[3][1] = wy; mat[3][2] = wz; mat[3][3] = ww;
}

SMF_INLINE_FORCED CMat4D::CMat4D( const CMat3D &rotation, const CVec3D &translation ) {
	// NOTE: CMat3D is transposed because it is column-major
	mat[ 0 ][ 0 ] = rotation[0][0];
	mat[ 0 ][ 1 ] = rotation[1][0];
	mat[ 0 ][ 2 ] = rotation[2][0];
	mat[ 0 ][ 3 ] = translation[0];
	mat[ 1 ][ 0 ] = rotation[0][1];
	mat[ 1 ][ 1 ] = rotation[1][1];
	mat[ 1 ][ 2 ] = rotation[2][1];
	mat[ 1 ][ 3 ] = translation[1];
	mat[ 2 ][ 0 ] = rotation[0][2];
	mat[ 2 ][ 1 ] = rotation[1][2];
	mat[ 2 ][ 2 ] = rotation[2][2];
	mat[ 2 ][ 3 ] = translation[2];
	mat[ 3 ][ 0 ] = 0.0f;
	mat[ 3 ][ 1 ] = 0.0f;
	mat[ 3 ][ 2 ] = 0.0f;
	mat[ 3 ][ 3 ] = 1.0f;
}

SMF_INLINE_FORCED CMat4D::CMat4D( const float src[ 4 ][ 4 ] ) {
	memcpy( mat, src, 4 * 4 * sizeof( float ) );
}

SMF_INLINE_FORCED const CVec4D &CMat4D::operator[]( int index ) const throw(Exception::CMathException){
	SMF_ASSERT( ( index >= 0 ) && ( index < 4 ) );
	if (( index < 0 ) || ( index >= 4 )) {
		string error("invalid index");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}
	return mat[ index ];
}

SMF_INLINE_FORCED CVec4D &CMat4D::operator[]( int index ) throw(Exception::CMathException){
	SMF_ASSERT( ( index >= 0 ) && ( index < 4 ) );
	if (( index < 0 ) || ( index >= 4 )) {
		string error("invalid index");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}
	return mat[ index ];
}

SMF_INLINE_FORCED CMat4D CMat4D::operator*( const float a ) const {
	return CMat4D(
		mat[0].x * a, mat[0].y * a, mat[0].z * a, mat[0].w * a,
		mat[1].x * a, mat[1].y * a, mat[1].z * a, mat[1].w * a,
		mat[2].x * a, mat[2].y * a, mat[2].z * a, mat[2].w * a,
		mat[3].x * a, mat[3].y * a, mat[3].z * a, mat[3].w * a );
}

SMF_INLINE_FORCED CVec4D CMat4D::operator*( const CVec4D &vec ) const {
	return CVec4D(
		mat[ 0 ].x * vec.x + mat[ 0 ].y * vec.y + mat[ 0 ].z * vec.z + mat[ 0 ].w * vec.w,
		mat[ 1 ].x * vec.x + mat[ 1 ].y * vec.y + mat[ 1 ].z * vec.z + mat[ 1 ].w * vec.w,
		mat[ 2 ].x * vec.x + mat[ 2 ].y * vec.y + mat[ 2 ].z * vec.z + mat[ 2 ].w * vec.w,
		mat[ 3 ].x * vec.x + mat[ 3 ].y * vec.y + mat[ 3 ].z * vec.z + mat[ 3 ].w * vec.w );
}

SMF_INLINE_FORCED CVec3D CMat4D::operator*( const CVec3D &vec ) const {
	float s = mat[ 3 ].x * vec.x + mat[ 3 ].y * vec.y + mat[ 3 ].z * vec.z + mat[ 3 ].w;
	if ( s == 0.0f ) {
		return CVec3D( 0.0f, 0.0f, 0.0f );
	}
	if ( s == 1.0f ) {
		return CVec3D(
			mat[ 0 ].x * vec.x + mat[ 0 ].y * vec.y + mat[ 0 ].z * vec.z + mat[ 0 ].w,
			mat[ 1 ].x * vec.x + mat[ 1 ].y * vec.y + mat[ 1 ].z * vec.z + mat[ 1 ].w,
			mat[ 2 ].x * vec.x + mat[ 2 ].y * vec.y + mat[ 2 ].z * vec.z + mat[ 2 ].w );
	}
	else {
		float invS = 1.0f / s;
		return CVec3D(
			(mat[ 0 ].x * vec.x + mat[ 0 ].y * vec.y + mat[ 0 ].z * vec.z + mat[ 0 ].w) * invS,
			(mat[ 1 ].x * vec.x + mat[ 1 ].y * vec.y + mat[ 1 ].z * vec.z + mat[ 1 ].w) * invS,
			(mat[ 2 ].x * vec.x + mat[ 2 ].y * vec.y + mat[ 2 ].z * vec.z + mat[ 2 ].w) * invS );
	}
}

SMF_INLINE_FORCED CMat4D CMat4D::operator*( const CMat4D &a ) const {
	CMat4D dst;
	if(CMath::MATH_AUTOMATIC_SSE){
		CSIMD::getProcessor()->mat4x4_mul_sse(dst.Row(0).toM128Ptr(),this->Row(0).toM128Ptr(),a.Row(0).toM128Ptr());
	}else{
		CSIMD::getGenProcessor()->mat4x4_mul_sse(dst.Row(0).toM128Ptr(),this->Row(0).toM128Ptr(),a.Row(0).toM128Ptr());
	}
	return dst;
}

SMF_INLINE_FORCED CMat4D CMat4D::operator *(const CMat3D &rhs) const
{
	///\todo SSE.
	CMat4D r;
	CVec3D col0=rhs.Col(0);
	CVec3D col1=rhs.Col(1);
	CVec3D col2=rhs.Col(2);

	r[0][0] = mat[0].toVec3()*col0;
	r[0][1] = mat[0].toVec3()*col1;
	r[0][2] = mat[0].toVec3()*col2;
	r[0][3] = mat[0][3];

	r[1][0] = mat[1].toVec3()*col0;
	r[1][1] = mat[1].toVec3()*col1;
	r[1][2] = mat[1].toVec3()*col2;
	r[1][3] = mat[1][3];

	r[2][0] = mat[2].toVec3()*col0;
	r[2][1] = mat[2].toVec3()*col1;
	r[2][2] = mat[2].toVec3()*col2;
	r[2][3] = mat[2][3];

	r[3][0] = mat[3].toVec3()*col0;
	r[3][1] = mat[3].toVec3()*col1;
	r[3][2] = mat[3].toVec3()*col2;
	r[3][3] = mat[3][3];

	return r;
}

SMF_INLINE_FORCED bool CMat4D::ContainsProjection(float epsilon) const
{
	return Row(3).compare(0.f, 0.f, 0.f, 1.f, epsilon) == false;
}

SMF_INLINE_FORCED CONST_WIN32 CMat3D CMat4D::CMat3DPart() const
{
	return CMat3D(mat[0][0], mat[0][1], mat[0][2],
					mat[1][0], mat[1][1], mat[1][2],
					mat[2][0], mat[2][1], mat[2][2]);
}

#if 0
CMatJoint3x4 &CMat4D::CMatJoint3x4Part()
{
	return CMatJoint3x4(mat[0][0], mat[0][1], mat[0][2],mat[0][3],
					mat[1][0], mat[1][1], mat[1][2],mat[1][3],
					mat[2][0], mat[2][1], mat[2][2],mat[2][3]);
}

const CMatJoint3x4 &CMat4D::CMatJoint3x4Part() const
{
	return CMatJoint3x4(mat[0][0], mat[0][1], mat[0][2],mat[0][3],
					mat[1][0], mat[1][1], mat[1][2],mat[1][3],
					mat[2][0], mat[2][1], mat[2][2],mat[2][3]);
}
#endif


SMF_INLINE_FORCED CMat4D CMat4D::operator+( const CMat4D &a ) const {
	return CMat4D(
		mat[0].x + a[0].x, mat[0].y + a[0].y, mat[0].z + a[0].z, mat[0].w + a[0].w,
		mat[1].x + a[1].x, mat[1].y + a[1].y, mat[1].z + a[1].z, mat[1].w + a[1].w,
		mat[2].x + a[2].x, mat[2].y + a[2].y, mat[2].z + a[2].z, mat[2].w + a[2].w,
		mat[3].x + a[3].x, mat[3].y + a[3].y, mat[3].z + a[3].z, mat[3].w + a[3].w );
}

SMF_INLINE_FORCED CMat4D CMat4D::operator-( const CMat4D &a ) const {
	return CMat4D(
		mat[0].x - a[0].x, mat[0].y - a[0].y, mat[0].z - a[0].z, mat[0].w - a[0].w,
		mat[1].x - a[1].x, mat[1].y - a[1].y, mat[1].z - a[1].z, mat[1].w - a[1].w,
		mat[2].x - a[2].x, mat[2].y - a[2].y, mat[2].z - a[2].z, mat[2].w - a[2].w,
		mat[3].x - a[3].x, mat[3].y - a[3].y, mat[3].z - a[3].z, mat[3].w - a[3].w );
}

SMF_INLINE_FORCED CMat4D &CMat4D::operator*=( const float a ) {
	mat[0].x *= a; mat[0].y *= a; mat[0].z *= a; mat[0].w *= a;
	mat[1].x *= a; mat[1].y *= a; mat[1].z *= a; mat[1].w *= a;
	mat[2].x *= a; mat[2].y *= a; mat[2].z *= a; mat[2].w *= a;
	mat[3].x *= a; mat[3].y *= a; mat[3].z *= a; mat[3].w *= a;
    return *this;
}

SMF_INLINE_FORCED CMat4D &CMat4D::operator*=( const CMat4D &a ) {
	*this = (*this) * a;
	return *this;
}

SMF_INLINE_FORCED CMat4D &CMat4D::operator+=( const CMat4D &a ) {
	mat[0].x += a[0].x; mat[0].y += a[0].y; mat[0].z += a[0].z; mat[0].w += a[0].w;
	mat[1].x += a[1].x; mat[1].y += a[1].y; mat[1].z += a[1].z; mat[1].w += a[1].w;
	mat[2].x += a[2].x; mat[2].y += a[2].y; mat[2].z += a[2].z; mat[2].w += a[2].w;
	mat[3].x += a[3].x; mat[3].y += a[3].y; mat[3].z += a[3].z; mat[3].w += a[3].w;
    return *this;
}

SMF_INLINE_FORCED CMat4D &CMat4D::operator-=( const CMat4D &a ) {
	mat[0].x -= a[0].x; mat[0].y -= a[0].y; mat[0].z -= a[0].z; mat[0].w -= a[0].w;
	mat[1].x -= a[1].x; mat[1].y -= a[1].y; mat[1].z -= a[1].z; mat[1].w -= a[1].w;
	mat[2].x -= a[2].x; mat[2].y -= a[2].y; mat[2].z -= a[2].z; mat[2].w -= a[2].w;
	mat[3].x -= a[3].x; mat[3].y -= a[3].y; mat[3].z -= a[3].z; mat[3].w -= a[3].w;
    return *this;
}

SMF_INLINE_FORCED CMat4D operator*( const float a, const CMat4D &mat ) {
	return mat * a;
}

SMF_INLINE_FORCED CVec4D operator*( const CVec4D &vec, const CMat4D &mat ) {
	return mat * vec;
}

SMF_INLINE_FORCED CVec3D operator*( const CVec3D &vec, const CMat4D &mat ) {
	return mat * vec;
}

SMF_INLINE_FORCED CVec4D &operator*=( CVec4D &vec, const CMat4D &mat ) {
	vec = mat * vec;
	return vec;
}

SMF_INLINE_FORCED CVec3D &operator*=( CVec3D &vec, const CMat4D &mat ) {
	vec = mat * vec;
	return vec;
}

SMF_INLINE_FORCED CVec4D &CMat4D::Row(int row)
{
	SMF_ASSERT(row >= 0);
	SMF_ASSERT(row < 4);

#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (row < 0 || row >= 4)
		row = 0; // Benign failure, just give the first row.
#endif

	return reinterpret_cast<CVec4D &>(mat[row]);
}

SMF_INLINE_FORCED const CVec4D &CMat4D::Row(int row) const
{
	SMF_ASSERT(row >= 0);
	SMF_ASSERT(row < 4);

#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (row < 0 || row >= 4)
		row = 0; // Benign failure, just give the first row.
#endif
	return reinterpret_cast<const CVec4D &>(mat[row]);
}

SMF_INLINE_FORCED CVec3D &CMat4D::Row3(int row)
{
	SMF_ASSERT(row >= 0);
	SMF_ASSERT(row < 4);

#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (row < 0 || row >= 4)
		row = 0; // Benign failure, just give the first row.
#endif
	return reinterpret_cast<CVec3D &>(mat[row]);
}

SMF_INLINE_FORCED const CVec3D &CMat4D::Row3(int row) const
{
	SMF_ASSERT(row >= 0);
	SMF_ASSERT(row < 4);

#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (row < 0 || row >= 4)
		row = 0; // Benign failure, just give the first row.
#endif
	return reinterpret_cast<const CVec3D &>(mat[row]);
}

SMF_INLINE_FORCED CONST_WIN32 CVec4D CMat4D::Col(int col) const
{
	SMF_ASSERT(col >= 0);
	SMF_ASSERT(col < 4);
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (col < 0 || col >= 4)
		return CVec4D::nan;
#endif
	return CVec4D(mat[0][col], mat[1][col], mat[2][col], mat[3][col]);
}

SMF_INLINE_FORCED CONST_WIN32 CVec3D CMat4D::Col3(int col) const
{
	SMF_ASSERT(col >= 0);
	SMF_ASSERT(col < 4);
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (col < 0 || col >= 4)
		return CVec3D::nan;
#endif
	return CVec3D(mat[0][col], mat[1][col], mat[2][col]);
}

SMF_INLINE_FORCED CONST_WIN32 CVec4D CMat4D::diagonal() const
{
	return CVec4D(mat[0][0], mat[1][1], mat[2][2], mat[3][3]);
}

SMF_INLINE_FORCED CONST_WIN32 CVec3D CMat4D::diagonal3() const
{
	return CVec3D(mat[0][0], mat[1][1], mat[2][2]);
}

SMF_INLINE_FORCED bool CMat4D::compare( const CMat4D &a ) const {
	sf_u32 i;
	const float *ptr1, *ptr2;

	ptr1 = reinterpret_cast<const float *>(mat);
	ptr2 = reinterpret_cast<const float *>(a.mat);
	for ( i = 0; i < 4*4; i++ ) {
		if ( ptr1[i] != ptr2[i] ) {
			return false;
		}
	}
	return true;
}

SMF_INLINE_FORCED bool CMat4D::compare( const CMat4D &a, const float epsilon ) const {
	sf_u32 i;
	const float *ptr1, *ptr2;

	ptr1 = reinterpret_cast<const float *>(mat);
	ptr2 = reinterpret_cast<const float *>(a.mat);
	for ( i = 0; i < 4*4; i++ ) {
		if ( CMath::fabs( ptr1[i] - ptr2[i] ) > epsilon ) {
			return false;
		}
	}
	return true;
}

SMF_INLINE_FORCED bool CMat4D::operator==( const CMat4D &a ) const {
	return compare( a );
}

SMF_INLINE_FORCED bool CMat4D::operator!=( const CMat4D &a ) const {
	return !compare( a );
}

SMF_INLINE_FORCED void CMat4D::toZero() {
	memset( mat, 0, sizeof( CMat4D ) );
}

SMF_INLINE_FORCED void CMat4D::toIdentity() {
	*this = mat4_identity;
}

SMF_INLINE_FORCED bool CMat4D::isIdentity( const float epsilon ) const {
	return compare( mat4_identity, epsilon );
#if 0
		for(int y = 0; y < 4; ++y)
		for(int x = 0; x < 4; ++x)
			if (!CMath::equalsAbs(mat[y][x], (x == y) ? 1.f : 0.f, epsilon))
				return false;

	return true;


#endif

}

SMF_INLINE_FORCED bool CMat4D::isSymmetric( const float epsilon ) const {
	for ( int i = 1; i < 4; i++ ) {
		for ( int j = 0; j < i; j++ ) {
			if ( CMath::fabs( mat[i][j] - mat[j][i] ) > epsilon ) {
				return false;
			}
		}
	}
	return true;

}

SMF_INLINE_FORCED bool CMat4D::isDiagonal( const float epsilon ) const {
	for ( int i = 0; i < 4; i++ ) {
		for ( int j = 0; j < 4; j++ ) {
			if ( i != j && CMath::fabs( mat[i][j] ) > epsilon ) {
				return false;
			}
		}
	}
	return true;
}

SMF_INLINE_FORCED bool CMat4D::isRotated() const {
	if ( !mat[ 0 ][ 1 ] && !mat[ 0 ][ 2 ] &&
		!mat[ 1 ][ 0 ] && !mat[ 1 ][ 2 ] &&
		!mat[ 2 ][ 0 ] && !mat[ 2 ][ 1 ] ) {
		return false;
	}
	return true;
}

SMF_INLINE_FORCED bool CMat4D::isFinite() const
{
	for(int y = 0; y < 4; ++y)
		for(int x = 0; x < 4; ++x)
			if (!MATH::isFinite(mat[y][x]))
				return false;
	return true;
}



SMF_INLINE_FORCED bool CMat4D::isLowerTriangular(float epsilon) const
{
	return CMath::equalsAbs(mat[0][1], 0.f, epsilon)
		&& CMath::equalsAbs(mat[0][2], 0.f, epsilon)
		&& CMath::equalsAbs(mat[0][3], 0.f, epsilon)
		&& CMath::equalsAbs(mat[1][2], 0.f, epsilon)
		&& CMath::equalsAbs(mat[1][3], 0.f, epsilon)
		&& CMath::equalsAbs(mat[2][3], 0.f, epsilon);
}

SMF_INLINE_FORCED bool CMat4D::isUpperTriangular(float epsilon) const
{
	return CMath::equalsAbs(mat[1][0], 0.f, epsilon)
		&& CMath::equalsAbs(mat[2][0], 0.f, epsilon)
		&& CMath::equalsAbs(mat[3][0], 0.f, epsilon)
		&& CMath::equalsAbs(mat[2][1], 0.f, epsilon)
		&& CMath::equalsAbs(mat[3][1], 0.f, epsilon)
		&& CMath::equalsAbs(mat[3][2], 0.f, epsilon);
}

SMF_INLINE_FORCED bool CMat4D::isInvertible(float epsilon) const
{
	///\todo Optimize.
	CMat4D copy = *this;
	return copy.inverseFastSelf();
}

SMF_INLINE_FORCED bool CMat4D::isSkewSymmetric(float epsilon) const
{
	for(int y = 0; y < 4; ++y)
		for(int x = y; x < 4; ++x)
			if (!CMath::equalsAbs(mat[y][x], -mat[x][y], epsilon))
				return false;
	return true;
}

SMF_INLINE_FORCED bool CMat4D::isIdempotent(float epsilon) const
{
	CMat4D m2 = *this * *this;
	return this->compare(m2, epsilon);
}


SMF_INLINE_FORCED void CMat4D::projectVector( const CVec4D &src, CVec4D &dst ) const {
	dst.x = src * mat[ 0 ];
	dst.y = src * mat[ 1 ];
	dst.z = src * mat[ 2 ];
	dst.w = src * mat[ 3 ];
}

SMF_INLINE_FORCED void CMat4D::unprojectVector( const CVec4D &src, CVec4D &dst ) const {
	dst = mat[ 0 ] * src.x + mat[ 1 ] * src.y + mat[ 2 ] * src.z + mat[ 3 ] * src.w;
}

SMF_INLINE_FORCED float CMat4D::trace() const {
	return ( mat[0][0] + mat[1][1] + mat[2][2] + mat[3][3] );
}



SMF_INLINE_FORCED CMat4D CMat4D::inverseFast() const throw(Exception::CMathException){



	CMat4D invMat;

	invMat = *this;
	int r = invMat.inverseFastSelf();
	if (r ==  0) {
		string error("matriz have no inverse");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}

	return invMat;
}

SMF_INLINE_FORCED CMat4D CMat3D::toMat4() const {
	// NOTE: CMat3D is transposed because it is column-major
	return CMat4D(	mat[0][0],	mat[1][0],	mat[2][0],	0.0f,
					mat[0][1],	mat[1][1],	mat[2][1],	0.0f,
					mat[0][2],	mat[1][2],	mat[2][2],	0.0f,
					0.0f,		0.0f,		0.0f,		1.0f );
}

SMF_INLINE_FORCED int CMat4D::getDimension() const {
	return 16;
}

SMF_INLINE_FORCED const float *CMat4D::toFloatPtr() const {
	return mat[0].toFloatPtr();
}

SMF_INLINE_FORCED float *CMat4D::toFloatPtr() {
	return mat[0].toFloatPtr();
}


static CVec3D mat3x4_mul_vec(const CVec4D *matrix, sf_m128 vector){

	CVec3D teste;
	return teste;
}

SMF_INLINE_FORCED CVec3D CMat4D::transformPos(const CVec3D &pointVector) const
{
	 SMF_ASSERT(!this->ContainsProjection()); // This function does not divide by w or output it, so cannot have projection.
if (CMath::MATH_AUTOMATIC_SSE) {
	return CSIMD::getProcessor()->mat3x4_mul_vec(mat[0].toM128Ptr(), MATH::SET_PS(1.f, pointVector.z, pointVector.y, pointVector.x));
}else{
	//return transformPos(pointVector.x, pointVector.y, pointVector.z);
	return CSIMD::getGenProcessor()->mat3x4_mul_vec(mat[0].toM128Ptr(), MATH::SET_PS(1.f, pointVector.z, pointVector.y, pointVector.x));
}
}

SMF_INLINE_FORCED CVec3D CMat4D::transformPos(float x, float y, float z) const
{
	SMF_ASSERT(!this->ContainsProjection()); // This function does not divide by w or output it, so cannot have projection.
	if (CMath::MATH_AUTOMATIC_SSE) {
		return CSIMD::getProcessor()->mat3x4_mul_vec(mat[0].toM128Ptr(), MATH::SET_PS(1.f, z, y, x));
	}else{
		return CSIMD::getProcessor()->mat3x4_mul_vec(mat[0].toM128Ptr(), MATH::SET_PS(1.f, z, y, x));

		//return CVec3D(Row(0).toVec3()* CVec3D(x,y,z),
		//			  Row(1).toVec3()* CVec3D(x,y,z),
		//			  Row(2).toVec3()* CVec3D(x,y,z));
	}
}

SMF_INLINE_FORCED CVec3D CMat4D::transformDir(const CVec3D &directionVector) const
{
	SMF_ASSERT(!this->ContainsProjection()); // This function does not divide by w or output it, so cannot have projection.
	if (CMath::MATH_AUTOMATIC_SSE) {
		return CSIMD::getProcessor()->mat3x4_mul_vec(mat[0].toM128Ptr(), MATH::SET_PS(0.f, directionVector.z, directionVector.y, directionVector.x));
	}else{
		return transformDir(directionVector.x, directionVector.y, directionVector.z);
	}
}

SMF_INLINE_FORCED CVec3D CMat4D::transformDir(float x, float y, float z) const
{
	SMF_ASSERT(!this->ContainsProjection()); // This function does not divide by w or output it, so cannot have projection.
	if (CMath::MATH_AUTOMATIC_SSE) {
		return CSIMD::getProcessor()->mat3x4_mul_vec(mat[0].toM128Ptr(), MATH::SET_PS(0.f, z, y, x));
	}else{
		return CVec3D(Row(0).toVec3()*CVec3D( x,y,z),
					  Row(0).toVec3()* CVec3D(x,y,z),
					  Row(0).toVec3()* CVec3D(x,y,z));
	}
}

SMF_INLINE_FORCED CVec4D CMat4D::transform(const CVec4D &vector) const
{
#if !defined(ANDROID) ///\bug Android GCC 4.6.6 gives internal compiler error!
if(CMath::MATH_AUTOMATIC_SSE){
	return CVec4D(CSIMD::getProcessor()->mat4x4_mul_sse(mat[0].toM128Ptr(), vector.v));
}else{
	return CVec4D(Row(0)* vector,
	              Row(1)* vector,
	              Row(2)* vector,
	              Row(3)* vector);
}
#endif
}

SMF_INLINE_FORCED CONST_WIN32 CVec3D CMat4D::translatePart() const
{
	return Col3(3);
}

SMF_INLINE_FORCED CONST_WIN32 CMat3D CMat4D::rotatePart() const
{
	return CMat3DPart();
}
SMF_INLINE_FORCED CONST_WIN32 CMat3D CMat4D::mat3Part() const
{
	return CMat3D(mat[0][0], mat[0][1], mat[0][2],
					mat[1][0], mat[1][1], mat[1][2],
					mat[2][0], mat[2][1], mat[2][2]);
}

SMF_INLINE_FORCED CMatJoint3x4 &CMat4D::mat3x4Part()
{
	return reinterpret_cast<CMatJoint3x4 &>(*this);
}

SMF_INLINE_FORCED const CMatJoint3x4 &CMat4D::mat3x4Part() const
{
	return reinterpret_cast<const CMatJoint3x4 &>(*this);
}
/**
 * \class CMat5D
 *
 * \ingroup SMF_Math
 *
 * \if pt_br
 * \brief Matriz 5x5
 *
 * \elseif us_en
 * \brief 5x5 matrix
 *
 * \endif
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
class SMF_API CMat5D {
public:
					CMat5D();
					explicit CMat5D( const CVec5D &v0, const CVec5D &v1, const CVec5D &v2, const CVec5D &v3, const CVec5D &v4 );
					explicit CMat5D( const float src[ 5 ][ 5 ] );

	const CVec5D &	operator[]( int index ) const throw(Exception::CMathException);
	CVec5D &		operator[]( int index ) throw(Exception::CMathException);
	CMat5D			operator*( const float a ) const;
	CVec5D			operator*( const CVec5D &vec ) const;
	CMat5D			operator*( const CMat5D &a ) const;
	CMat5D			operator+( const CMat5D &a ) const;
	CMat5D			operator-( const CMat5D &a ) const;
	CMat5D &		operator*=( const float a );
	CMat5D &		operator*=( const CMat5D &a );
	CMat5D &		operator+=( const CMat5D &a );
	CMat5D &		operator-=( const CMat5D &a );

	friend CMat5D	operator*( const float a, const CMat5D &mat );
	friend CVec5D	operator*( const CVec5D &vec, const CMat5D &mat );
	friend CVec5D &	operator*=( CVec5D &vec, const CMat5D &mat );
	/// exact compare, no epsilon
	bool			compare( const CMat5D &a ) const;
	/// compare with epsilon
	bool			compare( const CMat5D &a, const float epsilon ) const;
	/// exact compare, no epsilon
	bool			operator==( const CMat5D &a ) const;
	/// exact compare, no epsilon
	bool			operator!=( const CMat5D &a ) const;

	void			toZero();
	void			toIdentity();
	bool			isIdentity( const float epsilon = MATRIX_EPSILON ) const;
	bool			isSymmetric( const float epsilon = MATRIX_EPSILON ) const;
	bool			isDiagonal( const float epsilon = MATRIX_EPSILON ) const;

	float			trace() const;
	float			determinant() const;
	/// returns transpose
	CMat5D			transpose() const;
	CMat5D &		transposeSelf();
	/**
	\brief interte a matriz  ( m * m.inverse() = identity)
	\return retorna a matriz inversa ou throw exception, caso não haja inversa
	**/
	CMat5D		inverse() const throw(Exception::CMathException) ;
	/**
	\brief calcula a matriz inversa ( m * m.inverse() = identity )
	\return retorna falso caso não haja inversa e verdadeiro caso haja
	**/
	bool			inverseSelf();		// returns false if determinant is zero
	/**
	\brief interte a matriz  ( m * m.inverse() = identity)
	\return retorna a matriz inversa ou throw exception, caso não haja inversa
	**/
	CMat5D			inverseFast() const throw(Exception::CMathException);	// returns the inverse ( m * m.inverse() = identity )
	/**
	\brief calcula a matriz inversa ( m * m.inverse() = identity )
	\return retorna falso caso não haja inversa e verdadeiro caso haja
	**/
	bool			inverseFastSelf();	// returns false if determinant is zero


	int				getDimension() const;

	const float *	toFloatPtr() const;
	float *			toFloatPtr();
	const char *	toString( int precision = 2 ) const;

private:
	CVec5D			mat[ 5 ];
};

extern CMat5D mat5_zero;
extern CMat5D mat5_identity;
#define mat5_default	mat5_identity

SMF_INLINE_FORCED CMat5D::CMat5D() {
}

SMF_INLINE_FORCED CMat5D::CMat5D( const float src[ 5 ][ 5 ] ) {
	memcpy( mat, src, 5 * 5 * sizeof( float ) );
}

SMF_INLINE_FORCED CMat5D::CMat5D( const CVec5D &v0, const CVec5D &v1, const CVec5D &v2, const CVec5D &v3, const CVec5D &v4 ) {
	mat[0] = v0;
	mat[1] = v1;
	mat[2] = v2;
	mat[3] = v3;
	mat[4] = v4;
}

SMF_INLINE_FORCED const CVec5D &CMat5D::operator[]( int index ) const throw(Exception::CMathException){
	SMF_ASSERT( ( index >= 0 ) && ( index < 5 ) );
	if (( index < 0 ) || ( index >= 5 )) {
		string error("invalid index");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}
	return mat[ index ];
}

SMF_INLINE_FORCED CVec5D &CMat5D::operator[]( int index ) throw(Exception::CMathException){
	SMF_ASSERT( ( index >= 0 ) && ( index < 5 ) );
	if (( index < 0 ) || ( index >= 5 )) {
		string error("invalid index");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}
	return mat[ index ];
}

SMF_INLINE_FORCED CMat5D CMat5D::operator*( const CMat5D &a ) const {
	int i, j;
	const float *m1Ptr, *m2Ptr;
	float *dstPtr;
	CMat5D dst;

	m1Ptr = reinterpret_cast<const float *>(this);
	m2Ptr = reinterpret_cast<const float *>(&a);
	dstPtr = reinterpret_cast<float *>(&dst);

	for ( i = 0; i < 5; i++ ) {
		for ( j = 0; j < 5; j++ ) {
			*dstPtr = m1Ptr[0] * m2Ptr[ 0 * 5 + j ]
					+ m1Ptr[1] * m2Ptr[ 1 * 5 + j ]
					+ m1Ptr[2] * m2Ptr[ 2 * 5 + j ]
					+ m1Ptr[3] * m2Ptr[ 3 * 5 + j ]
					+ m1Ptr[4] * m2Ptr[ 4 * 5 + j ];
			dstPtr++;
		}
		m1Ptr += 5;
	}
	return dst;
}

SMF_INLINE_FORCED CMat5D CMat5D::operator*( const float a ) const {
	return CMat5D(
		CVec5D( mat[0][0] * a, mat[0][1] * a, mat[0][2] * a, mat[0][3] * a, mat[0][4] * a ),
		CVec5D( mat[1][0] * a, mat[1][1] * a, mat[1][2] * a, mat[1][3] * a, mat[1][4] * a ),
		CVec5D( mat[2][0] * a, mat[2][1] * a, mat[2][2] * a, mat[2][3] * a, mat[2][4] * a ),
		CVec5D( mat[3][0] * a, mat[3][1] * a, mat[3][2] * a, mat[3][3] * a, mat[3][4] * a ),
		CVec5D( mat[4][0] * a, mat[4][1] * a, mat[4][2] * a, mat[4][3] * a, mat[4][4] * a ) );
}

SMF_INLINE_FORCED CVec5D CMat5D::operator*( const CVec5D &vec ) const {
	return CVec5D(
		mat[0][0] * vec[0] + mat[0][1] * vec[1] + mat[0][2] * vec[2] + mat[0][3] * vec[3] + mat[0][4] * vec[4],
		mat[1][0] * vec[0] + mat[1][1] * vec[1] + mat[1][2] * vec[2] + mat[1][3] * vec[3] + mat[1][4] * vec[4],
		mat[2][0] * vec[0] + mat[2][1] * vec[1] + mat[2][2] * vec[2] + mat[2][3] * vec[3] + mat[2][4] * vec[4],
		mat[3][0] * vec[0] + mat[3][1] * vec[1] + mat[3][2] * vec[2] + mat[3][3] * vec[3] + mat[3][4] * vec[4],
		mat[4][0] * vec[0] + mat[4][1] * vec[1] + mat[4][2] * vec[2] + mat[4][3] * vec[3] + mat[4][4] * vec[4] );
}

SMF_INLINE_FORCED CMat5D CMat5D::operator+( const CMat5D &a ) const {
	return CMat5D(
		CVec5D( mat[0][0] + a[0][0], mat[0][1] + a[0][1], mat[0][2] + a[0][2], mat[0][3] + a[0][3], mat[0][4] + a[0][4] ),
		CVec5D( mat[1][0] + a[1][0], mat[1][1] + a[1][1], mat[1][2] + a[1][2], mat[1][3] + a[1][3], mat[1][4] + a[1][4] ),
		CVec5D( mat[2][0] + a[2][0], mat[2][1] + a[2][1], mat[2][2] + a[2][2], mat[2][3] + a[2][3], mat[2][4] + a[2][4] ),
		CVec5D( mat[3][0] + a[3][0], mat[3][1] + a[3][1], mat[3][2] + a[3][2], mat[3][3] + a[3][3], mat[3][4] + a[3][4] ),
		CVec5D( mat[4][0] + a[4][0], mat[4][1] + a[4][1], mat[4][2] + a[4][2], mat[4][3] + a[4][3], mat[4][4] + a[4][4] ) );
}

SMF_INLINE_FORCED CMat5D CMat5D::operator-( const CMat5D &a ) const {
	return CMat5D(
		CVec5D( mat[0][0] - a[0][0], mat[0][1] - a[0][1], mat[0][2] - a[0][2], mat[0][3] - a[0][3], mat[0][4] - a[0][4] ),
		CVec5D( mat[1][0] - a[1][0], mat[1][1] - a[1][1], mat[1][2] - a[1][2], mat[1][3] - a[1][3], mat[1][4] - a[1][4] ),
		CVec5D( mat[2][0] - a[2][0], mat[2][1] - a[2][1], mat[2][2] - a[2][2], mat[2][3] - a[2][3], mat[2][4] - a[2][4] ),
		CVec5D( mat[3][0] - a[3][0], mat[3][1] - a[3][1], mat[3][2] - a[3][2], mat[3][3] - a[3][3], mat[3][4] - a[3][4] ),
		CVec5D( mat[4][0] - a[4][0], mat[4][1] - a[4][1], mat[4][2] - a[4][2], mat[4][3] - a[4][3], mat[4][4] - a[4][4] ) );
}

SMF_INLINE_FORCED CMat5D &CMat5D::operator*=( const float a ) {
	mat[0][0] *= a; mat[0][1] *= a; mat[0][2] *= a; mat[0][3] *= a; mat[0][4] *= a;
	mat[1][0] *= a; mat[1][1] *= a; mat[1][2] *= a; mat[1][3] *= a; mat[1][4] *= a;
	mat[2][0] *= a; mat[2][1] *= a; mat[2][2] *= a; mat[2][3] *= a; mat[2][4] *= a;
	mat[3][0] *= a; mat[3][1] *= a; mat[3][2] *= a; mat[3][3] *= a; mat[3][4] *= a;
	mat[4][0] *= a; mat[4][1] *= a; mat[4][2] *= a; mat[4][3] *= a; mat[4][4] *= a;
	return *this;
}

SMF_INLINE_FORCED CMat5D &CMat5D::operator*=( const CMat5D &a ) {
	*this = *this * a;
	return *this;
}

SMF_INLINE_FORCED CMat5D &CMat5D::operator+=( const CMat5D &a ) {
	mat[0][0] += a[0][0]; mat[0][1] += a[0][1]; mat[0][2] += a[0][2]; mat[0][3] += a[0][3]; mat[0][4] += a[0][4];
	mat[1][0] += a[1][0]; mat[1][1] += a[1][1]; mat[1][2] += a[1][2]; mat[1][3] += a[1][3]; mat[1][4] += a[1][4];
	mat[2][0] += a[2][0]; mat[2][1] += a[2][1]; mat[2][2] += a[2][2]; mat[2][3] += a[2][3]; mat[2][4] += a[2][4];
	mat[3][0] += a[3][0]; mat[3][1] += a[3][1]; mat[3][2] += a[3][2]; mat[3][3] += a[3][3]; mat[3][4] += a[3][4];
	mat[4][0] += a[4][0]; mat[4][1] += a[4][1]; mat[4][2] += a[4][2]; mat[4][3] += a[4][3]; mat[4][4] += a[4][4];
	return *this;
}

SMF_INLINE_FORCED CMat5D &CMat5D::operator-=( const CMat5D &a ) {
	mat[0][0] -= a[0][0]; mat[0][1] -= a[0][1]; mat[0][2] -= a[0][2]; mat[0][3] -= a[0][3]; mat[0][4] -= a[0][4];
	mat[1][0] -= a[1][0]; mat[1][1] -= a[1][1]; mat[1][2] -= a[1][2]; mat[1][3] -= a[1][3]; mat[1][4] -= a[1][4];
	mat[2][0] -= a[2][0]; mat[2][1] -= a[2][1]; mat[2][2] -= a[2][2]; mat[2][3] -= a[2][3]; mat[2][4] -= a[2][4];
	mat[3][0] -= a[3][0]; mat[3][1] -= a[3][1]; mat[3][2] -= a[3][2]; mat[3][3] -= a[3][3]; mat[3][4] -= a[3][4];
	mat[4][0] -= a[4][0]; mat[4][1] -= a[4][1]; mat[4][2] -= a[4][2]; mat[4][3] -= a[4][3]; mat[4][4] -= a[4][4];
	return *this;
}

SMF_INLINE_FORCED CVec5D operator*( const CVec5D &vec, const CMat5D &mat ) {
	return mat * vec;
}

SMF_INLINE_FORCED CMat5D operator*( const float a, CMat5D const &mat ) {
	return mat * a;
}

SMF_INLINE_FORCED CVec5D &operator*=( CVec5D &vec, const CMat5D &mat ) {
	vec = mat * vec;
	return vec;
}

SMF_INLINE_FORCED bool CMat5D::compare( const CMat5D &a ) const {
	sf_u32 i;
	const float *ptr1, *ptr2;

	ptr1 = reinterpret_cast<const float *>(mat);
	ptr2 = reinterpret_cast<const float *>(a.mat);
	for ( i = 0; i < 5*5; i++ ) {
		if ( ptr1[i] != ptr2[i] ) {
			return false;
		}
	}
	return true;
}

SMF_INLINE_FORCED bool CMat5D::compare( const CMat5D &a, const float epsilon ) const {
	sf_u32 i;
	const float *ptr1, *ptr2;

	ptr1 = reinterpret_cast<const float *>(mat);
	ptr2 = reinterpret_cast<const float *>(a.mat);
	for ( i = 0; i < 5*5; i++ ) {
		if ( CMath::fabs( ptr1[i] - ptr2[i] ) > epsilon ) {
			return false;
		}
	}
	return true;
}

SMF_INLINE_FORCED bool CMat5D::operator==( const CMat5D &a ) const {
	return compare( a );
}

SMF_INLINE_FORCED bool CMat5D::operator!=( const CMat5D &a ) const {
	return !compare( a );
}

SMF_INLINE_FORCED void CMat5D::toZero() {
	memset( mat, 0, sizeof( CMat5D ) );
}

SMF_INLINE_FORCED void CMat5D::toIdentity() {
	*this = mat5_identity;
}

SMF_INLINE_FORCED bool CMat5D::isIdentity( const float epsilon ) const {
	return compare( mat5_identity, epsilon );
}

SMF_INLINE_FORCED bool CMat5D::isSymmetric( const float epsilon ) const {
	for ( int i = 1; i < 5; i++ ) {
		for ( int j = 0; j < i; j++ ) {
			if ( CMath::fabs( mat[i][j] - mat[j][i] ) > epsilon ) {
				return false;
			}
		}
	}
	return true;
}

SMF_INLINE_FORCED bool CMat5D::isDiagonal( const float epsilon ) const {
	for ( int i = 0; i < 5; i++ ) {
		for ( int j = 0; j < 5; j++ ) {
			if ( i != j && CMath::fabs( mat[i][j] ) > epsilon ) {
				return false;
			}
		}
	}
	return true;
}

SMF_INLINE_FORCED float CMat5D::trace() const {
	return ( mat[0][0] + mat[1][1] + mat[2][2] + mat[3][3] + mat[4][4] );
}

SMF_INLINE_FORCED CMat5D CMat5D::inverse() const throw(Exception::CMathException){
	CMat5D invMat;

	invMat = *this;
	int r = invMat.inverseSelf();
	if (r ==  0) {
		string error("matriz have no inverse");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}

	return invMat;
}

SMF_INLINE_FORCED CMat5D CMat5D::inverseFast() const throw(Exception::CMathException){
	CMat5D invMat;

	invMat = *this;
	int r = invMat.inverseFastSelf();
	if (r ==  0) {
		string error("matriz have no inverse");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}

	return invMat;
}

SMF_INLINE_FORCED int CMat5D::getDimension() const {
	return 25;
}

SMF_INLINE_FORCED const float *CMat5D::toFloatPtr() const {
	return mat[0].toFloatPtr();
}

SMF_INLINE_FORCED float *CMat5D::toFloatPtr() {
	return mat[0].toFloatPtr();
}


/**
 * \class CMat6D
 *
 * \ingroup SMF_Math
 *
 * \if pt_br
 * \brief Matriz 6x6
 *
 * \elseif us_en
 * \brief 6x6 matrix
 *
 * \endif
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
class SMF_API CMat6D {
public:
					CMat6D();
					explicit CMat6D( const CVec6D &v0, const CVec6D &v1, const CVec6D &v2, const CVec6D &v3, const CVec6D &v4, const CVec6D &v5 );
					explicit CMat6D( const CMat3D &m0, const CMat3D &m1, const CMat3D &m2, const CMat3D &m3 );
					explicit CMat6D( const float src[ 6 ][ 6 ] );

	const CVec6D &	operator[]( int index ) const throw(Exception::CMathException);
	CVec6D &		operator[]( int index ) throw(Exception::CMathException);
	CMat6D			operator*( const float a ) const;
	CVec6D			operator*( const CVec6D &vec ) const;
	CMat6D			operator*( const CMat6D &a ) const;
	CMat6D			operator+( const CMat6D &a ) const;
	CMat6D			operator-( const CMat6D &a ) const;
	CMat6D &		operator*=( const float a );
	CMat6D &		operator*=( const CMat6D &a );
	CMat6D &		operator+=( const CMat6D &a );
	CMat6D &		operator-=( const CMat6D &a );

	friend CMat6D	operator*( const float a, const CMat6D &mat );
	friend CVec6D	operator*( const CVec6D &vec, const CMat6D &mat );
	friend CVec6D &	operator*=( CVec6D &vec, const CMat6D &mat );
	/// exact compare, no epsilon
	bool			compare( const CMat6D &a ) const;
	/// compare with epsilon
	bool			compare( const CMat6D &a, const float epsilon ) const;
	/// exact compare, no epsilon
	bool			operator==( const CMat6D &a ) const;
	/// exact compare, no epsilon
	bool			operator!=( const CMat6D &a ) const;

	void			toZero();
	void			toIdentity();
	bool			isIdentity( const float epsilon = MATRIX_EPSILON ) const;
	bool			isSymmetric( const float epsilon = MATRIX_EPSILON ) const;
	bool			isDiagonal( const float epsilon = MATRIX_EPSILON ) const;

	CMat3D			subMat3( int n ) const;
	float			trace() const;
	float			determinant() const;
	CMat6D			transpose() const;	// returns transpose
	CMat6D &		transposeSelf();
	/**
	\brief interte a matriz  ( m * m.inverse() = identity)
	\return retorna a matriz inversa ou throw exception, caso não haja inversa
	**/
	CMat6D		inverse() const throw(Exception::CMathException) ;		// returns the inverse ( m * m.inverse() = identity )
	/**
	\brief calcula a matriz inversa ( m * m.inverse() = identity )
	\return retorna falso caso não haja inversa e verdadeiro caso haja
	**/
	bool			inverseSelf();		// returns false if determinant is zero
	/**
	\brief interte a matriz  ( m * m.inverse() = identity)
	\return retorna a matriz inversa ou throw exception, caso não haja inversa
	**/
	CMat6D			inverseFast() const throw(Exception::CMathException);	// returns the inverse ( m * m.inverse() = identity )
	/**
	\brief calcula a matriz inversa ( m * m.inverse() = identity )
	\return retorna falso caso não haja inversa e verdadeiro caso haja
	**/
	bool			inverseFastSelf();	// returns false if determinant is zero

	int				getDimension() const;

	const float *	toFloatPtr() const;
	float *			toFloatPtr();
	const char *	toString( int precision = 2 ) const;

private:
	CVec6D			mat[ 6 ];
};

extern CMat6D mat6_zero;
extern CMat6D mat6_identity;
#define mat6_default	mat6_identity

SMF_INLINE_FORCED CMat6D::CMat6D() {
}

SMF_INLINE_FORCED CMat6D::CMat6D( const CMat3D &m0, const CMat3D &m1, const CMat3D &m2, const CMat3D &m3 ) {
	mat[0] = CVec6D( m0[0][0], m0[0][1], m0[0][2], m1[0][0], m1[0][1], m1[0][2] );
	mat[1] = CVec6D( m0[1][0], m0[1][1], m0[1][2], m1[1][0], m1[1][1], m1[1][2] );
	mat[2] = CVec6D( m0[2][0], m0[2][1], m0[2][2], m1[2][0], m1[2][1], m1[2][2] );
	mat[3] = CVec6D( m2[0][0], m2[0][1], m2[0][2], m3[0][0], m3[0][1], m3[0][2] );
	mat[4] = CVec6D( m2[1][0], m2[1][1], m2[1][2], m3[1][0], m3[1][1], m3[1][2] );
	mat[5] = CVec6D( m2[2][0], m2[2][1], m2[2][2], m3[2][0], m3[2][1], m3[2][2] );
}

SMF_INLINE_FORCED CMat6D::CMat6D( const CVec6D &v0, const CVec6D &v1, const CVec6D &v2, const CVec6D &v3, const CVec6D &v4, const CVec6D &v5 ) {
	mat[0] = v0;
	mat[1] = v1;
	mat[2] = v2;
	mat[3] = v3;
	mat[4] = v4;
	mat[5] = v5;
}

SMF_INLINE_FORCED CMat6D::CMat6D( const float src[ 6 ][ 6 ] ) {
	memcpy( mat, src, 6 * 6 * sizeof( float ) );
}

SMF_INLINE_FORCED const CVec6D &CMat6D::operator[]( int index ) const throw(Exception::CMathException) {
	SMF_ASSERT( ( index >= 0 ) && ( index < 6 ) );
	if (( index < 0 ) || ( index >= 6 )) {
		string error("invalid index");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}
	return mat[ index ];
}

SMF_INLINE_FORCED CVec6D &CMat6D::operator[]( int index ) throw(Exception::CMathException) {
	SMF_ASSERT( ( index >= 0 ) && ( index < 6 ) );
	if (( index < 0 ) || ( index >= 6 )) {
		string error("invalid index");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}
	return mat[ index ];
}

SMF_INLINE_FORCED CMat6D CMat6D::operator*( const CMat6D &a ) const {
	int i, j;
	const float *m1Ptr, *m2Ptr;
	float *dstPtr;
	CMat6D dst;

	m1Ptr = reinterpret_cast<const float *>(this);
	m2Ptr = reinterpret_cast<const float *>(&a);
	dstPtr = reinterpret_cast<float *>(&dst);

	for ( i = 0; i < 6; i++ ) {
		for ( j = 0; j < 6; j++ ) {
			*dstPtr = m1Ptr[0] * m2Ptr[ 0 * 6 + j ]
					+ m1Ptr[1] * m2Ptr[ 1 * 6 + j ]
					+ m1Ptr[2] * m2Ptr[ 2 * 6 + j ]
					+ m1Ptr[3] * m2Ptr[ 3 * 6 + j ]
					+ m1Ptr[4] * m2Ptr[ 4 * 6 + j ]
					+ m1Ptr[5] * m2Ptr[ 5 * 6 + j ];
			dstPtr++;
		}
		m1Ptr += 6;
	}
	return dst;
}

SMF_INLINE_FORCED CMat6D CMat6D::operator*( const float a ) const {
	return CMat6D(
		CVec6D( mat[0][0] * a, mat[0][1] * a, mat[0][2] * a, mat[0][3] * a, mat[0][4] * a, mat[0][5] * a ),
		CVec6D( mat[1][0] * a, mat[1][1] * a, mat[1][2] * a, mat[1][3] * a, mat[1][4] * a, mat[1][5] * a ),
		CVec6D( mat[2][0] * a, mat[2][1] * a, mat[2][2] * a, mat[2][3] * a, mat[2][4] * a, mat[2][5] * a ),
		CVec6D( mat[3][0] * a, mat[3][1] * a, mat[3][2] * a, mat[3][3] * a, mat[3][4] * a, mat[3][5] * a ),
		CVec6D( mat[4][0] * a, mat[4][1] * a, mat[4][2] * a, mat[4][3] * a, mat[4][4] * a, mat[4][5] * a ),
		CVec6D( mat[5][0] * a, mat[5][1] * a, mat[5][2] * a, mat[5][3] * a, mat[5][4] * a, mat[5][5] * a ) );
}

SMF_INLINE_FORCED CVec6D CMat6D::operator*( const CVec6D &vec ) const {
	return CVec6D(
		mat[0][0] * vec[0] + mat[0][1] * vec[1] + mat[0][2] * vec[2] + mat[0][3] * vec[3] + mat[0][4] * vec[4] + mat[0][5] * vec[5],
		mat[1][0] * vec[0] + mat[1][1] * vec[1] + mat[1][2] * vec[2] + mat[1][3] * vec[3] + mat[1][4] * vec[4] + mat[1][5] * vec[5],
		mat[2][0] * vec[0] + mat[2][1] * vec[1] + mat[2][2] * vec[2] + mat[2][3] * vec[3] + mat[2][4] * vec[4] + mat[2][5] * vec[5],
		mat[3][0] * vec[0] + mat[3][1] * vec[1] + mat[3][2] * vec[2] + mat[3][3] * vec[3] + mat[3][4] * vec[4] + mat[3][5] * vec[5],
		mat[4][0] * vec[0] + mat[4][1] * vec[1] + mat[4][2] * vec[2] + mat[4][3] * vec[3] + mat[4][4] * vec[4] + mat[4][5] * vec[5],
		mat[5][0] * vec[0] + mat[5][1] * vec[1] + mat[5][2] * vec[2] + mat[5][3] * vec[3] + mat[5][4] * vec[4] + mat[5][5] * vec[5] );
}

SMF_INLINE_FORCED CMat6D CMat6D::operator+( const CMat6D &a ) const {
	return CMat6D(
		CVec6D( mat[0][0] + a[0][0], mat[0][1] + a[0][1], mat[0][2] + a[0][2], mat[0][3] + a[0][3], mat[0][4] + a[0][4], mat[0][5] + a[0][5] ),
		CVec6D( mat[1][0] + a[1][0], mat[1][1] + a[1][1], mat[1][2] + a[1][2], mat[1][3] + a[1][3], mat[1][4] + a[1][4], mat[1][5] + a[1][5] ),
		CVec6D( mat[2][0] + a[2][0], mat[2][1] + a[2][1], mat[2][2] + a[2][2], mat[2][3] + a[2][3], mat[2][4] + a[2][4], mat[2][5] + a[2][5] ),
		CVec6D( mat[3][0] + a[3][0], mat[3][1] + a[3][1], mat[3][2] + a[3][2], mat[3][3] + a[3][3], mat[3][4] + a[3][4], mat[3][5] + a[3][5] ),
		CVec6D( mat[4][0] + a[4][0], mat[4][1] + a[4][1], mat[4][2] + a[4][2], mat[4][3] + a[4][3], mat[4][4] + a[4][4], mat[4][5] + a[4][5] ),
		CVec6D( mat[5][0] + a[5][0], mat[5][1] + a[5][1], mat[5][2] + a[5][2], mat[5][3] + a[5][3], mat[5][4] + a[5][4], mat[5][5] + a[5][5] ) );
}

SMF_INLINE_FORCED CMat6D CMat6D::operator-( const CMat6D &a ) const {
	return CMat6D(
		CVec6D( mat[0][0] - a[0][0], mat[0][1] - a[0][1], mat[0][2] - a[0][2], mat[0][3] - a[0][3], mat[0][4] - a[0][4], mat[0][5] - a[0][5] ),
		CVec6D( mat[1][0] - a[1][0], mat[1][1] - a[1][1], mat[1][2] - a[1][2], mat[1][3] - a[1][3], mat[1][4] - a[1][4], mat[1][5] - a[1][5] ),
		CVec6D( mat[2][0] - a[2][0], mat[2][1] - a[2][1], mat[2][2] - a[2][2], mat[2][3] - a[2][3], mat[2][4] - a[2][4], mat[2][5] - a[2][5] ),
		CVec6D( mat[3][0] - a[3][0], mat[3][1] - a[3][1], mat[3][2] - a[3][2], mat[3][3] - a[3][3], mat[3][4] - a[3][4], mat[3][5] - a[3][5] ),
		CVec6D( mat[4][0] - a[4][0], mat[4][1] - a[4][1], mat[4][2] - a[4][2], mat[4][3] - a[4][3], mat[4][4] - a[4][4], mat[4][5] - a[4][5] ),
		CVec6D( mat[5][0] - a[5][0], mat[5][1] - a[5][1], mat[5][2] - a[5][2], mat[5][3] - a[5][3], mat[5][4] - a[5][4], mat[5][5] - a[5][5] ) );
}

SMF_INLINE_FORCED CMat6D &CMat6D::operator*=( const float a ) {
	mat[0][0] *= a; mat[0][1] *= a; mat[0][2] *= a; mat[0][3] *= a; mat[0][4] *= a; mat[0][5] *= a;
	mat[1][0] *= a; mat[1][1] *= a; mat[1][2] *= a; mat[1][3] *= a; mat[1][4] *= a; mat[1][5] *= a;
	mat[2][0] *= a; mat[2][1] *= a; mat[2][2] *= a; mat[2][3] *= a; mat[2][4] *= a; mat[2][5] *= a;
	mat[3][0] *= a; mat[3][1] *= a; mat[3][2] *= a; mat[3][3] *= a; mat[3][4] *= a; mat[3][5] *= a;
	mat[4][0] *= a; mat[4][1] *= a; mat[4][2] *= a; mat[4][3] *= a; mat[4][4] *= a; mat[4][5] *= a;
	mat[5][0] *= a; mat[5][1] *= a; mat[5][2] *= a; mat[5][3] *= a; mat[5][4] *= a; mat[5][5] *= a;
	return *this;
}

SMF_INLINE_FORCED CMat6D &CMat6D::operator*=( const CMat6D &a ) {
	*this = *this * a;
	return *this;
}

SMF_INLINE_FORCED CMat6D &CMat6D::operator+=( const CMat6D &a ) {
	mat[0][0] += a[0][0]; mat[0][1] += a[0][1]; mat[0][2] += a[0][2]; mat[0][3] += a[0][3]; mat[0][4] += a[0][4]; mat[0][5] += a[0][5];
	mat[1][0] += a[1][0]; mat[1][1] += a[1][1]; mat[1][2] += a[1][2]; mat[1][3] += a[1][3]; mat[1][4] += a[1][4]; mat[1][5] += a[1][5];
	mat[2][0] += a[2][0]; mat[2][1] += a[2][1]; mat[2][2] += a[2][2]; mat[2][3] += a[2][3]; mat[2][4] += a[2][4]; mat[2][5] += a[2][5];
	mat[3][0] += a[3][0]; mat[3][1] += a[3][1]; mat[3][2] += a[3][2]; mat[3][3] += a[3][3]; mat[3][4] += a[3][4]; mat[3][5] += a[3][5];
	mat[4][0] += a[4][0]; mat[4][1] += a[4][1]; mat[4][2] += a[4][2]; mat[4][3] += a[4][3]; mat[4][4] += a[4][4]; mat[4][5] += a[4][5];
	mat[5][0] += a[5][0]; mat[5][1] += a[5][1]; mat[5][2] += a[5][2]; mat[5][3] += a[5][3]; mat[5][4] += a[5][4]; mat[5][5] += a[5][5];
	return *this;
}

SMF_INLINE_FORCED CMat6D &CMat6D::operator-=( const CMat6D &a ) {
	mat[0][0] -= a[0][0]; mat[0][1] -= a[0][1]; mat[0][2] -= a[0][2]; mat[0][3] -= a[0][3]; mat[0][4] -= a[0][4]; mat[0][5] -= a[0][5];
	mat[1][0] -= a[1][0]; mat[1][1] -= a[1][1]; mat[1][2] -= a[1][2]; mat[1][3] -= a[1][3]; mat[1][4] -= a[1][4]; mat[1][5] -= a[1][5];
	mat[2][0] -= a[2][0]; mat[2][1] -= a[2][1]; mat[2][2] -= a[2][2]; mat[2][3] -= a[2][3]; mat[2][4] -= a[2][4]; mat[2][5] -= a[2][5];
	mat[3][0] -= a[3][0]; mat[3][1] -= a[3][1]; mat[3][2] -= a[3][2]; mat[3][3] -= a[3][3]; mat[3][4] -= a[3][4]; mat[3][5] -= a[3][5];
	mat[4][0] -= a[4][0]; mat[4][1] -= a[4][1]; mat[4][2] -= a[4][2]; mat[4][3] -= a[4][3]; mat[4][4] -= a[4][4]; mat[4][5] -= a[4][5];
	mat[5][0] -= a[5][0]; mat[5][1] -= a[5][1]; mat[5][2] -= a[5][2]; mat[5][3] -= a[5][3]; mat[5][4] -= a[5][4]; mat[5][5] -= a[5][5];
	return *this;
}

SMF_INLINE_FORCED CVec6D operator*( const CVec6D &vec, const CMat6D &mat ) {
	return mat * vec;
}

SMF_INLINE_FORCED CMat6D operator*( const float a, CMat6D const &mat ) {
	return mat * a;
}

SMF_INLINE_FORCED CVec6D &operator*=( CVec6D &vec, const CMat6D &mat ) {
	vec = mat * vec;
	return vec;
}

SMF_INLINE_FORCED bool CMat6D::compare( const CMat6D &a ) const {
	sf_u32 i;
	const float *ptr1, *ptr2;

	ptr1 = reinterpret_cast<const float *>(mat);
	ptr2 = reinterpret_cast<const float *>(a.mat);
	for ( i = 0; i < 6*6; i++ ) {
		if ( ptr1[i] != ptr2[i] ) {
			return false;
		}
	}
	return true;
}

SMF_INLINE_FORCED bool CMat6D::compare( const CMat6D &a, const float epsilon ) const {
	sf_u32 i;
	const float *ptr1, *ptr2;

	ptr1 = reinterpret_cast<const float *>(mat);
	ptr2 = reinterpret_cast<const float *>(a.mat);
	for ( i = 0; i < 6*6; i++ ) {
		if ( CMath::fabs( ptr1[i] - ptr2[i] ) > epsilon ) {
			return false;
		}
	}
	return true;
}

SMF_INLINE_FORCED bool CMat6D::operator==( const CMat6D &a ) const {
	return compare( a );
}

SMF_INLINE_FORCED bool CMat6D::operator!=( const CMat6D &a ) const {
	return !compare( a );
}

SMF_INLINE_FORCED void CMat6D::toZero() {
	memset( mat, 0, sizeof( CMat6D ) );
}

SMF_INLINE_FORCED void CMat6D::toIdentity() {
	*this = mat6_identity;
}

SMF_INLINE_FORCED bool CMat6D::isIdentity( const float epsilon ) const {
	return compare( mat6_identity, epsilon );
}

SMF_INLINE_FORCED bool CMat6D::isSymmetric( const float epsilon ) const {
	for ( int i = 1; i < 6; i++ ) {
		for ( int j = 0; j < i; j++ ) {
			if ( CMath::fabs( mat[i][j] - mat[j][i] ) > epsilon ) {
				return false;
			}
		}
	}
	return true;
}

SMF_INLINE_FORCED bool CMat6D::isDiagonal( const float epsilon ) const {
	for ( int i = 0; i < 6; i++ ) {
		for ( int j = 0; j < 6; j++ ) {
			if ( i != j && CMath::fabs( mat[i][j] ) > epsilon ) {
				return false;
			}
		}
	}
	return true;
}

SMF_INLINE_FORCED CMat3D CMat6D::subMat3( int n ) const {
	SMF_ASSERT( n >= 0 && n < 4 );
	if (( n < 0 ) || ( n >= 4 )) {
		string error("invalid index");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}

	int b0 = ((n & 2) >> 1) * 3;
	int b1 = (n & 1) * 3;
	return CMat3D(
		mat[b0 + 0][b1 + 0], mat[b0 + 0][b1 + 1], mat[b0 + 0][b1 + 2],
		mat[b0 + 1][b1 + 0], mat[b0 + 1][b1 + 1], mat[b0 + 1][b1 + 2],
		mat[b0 + 2][b1 + 0], mat[b0 + 2][b1 + 1], mat[b0 + 2][b1 + 2] );
}

SMF_INLINE_FORCED float CMat6D::trace() const {
	return ( mat[0][0] + mat[1][1] + mat[2][2] + mat[3][3] + mat[4][4] + mat[5][5] );
}

SMF_INLINE_FORCED CMat6D CMat6D::inverse() const throw(Exception::CMathException){
	CMat6D invMat;

	invMat = *this;
	int r = invMat.inverseSelf();
	if (r ==  0) {
		string error("matriz have no inverse");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}

	return invMat;
}

SMF_INLINE_FORCED CMat6D CMat6D::inverseFast() const throw(Exception::CMathException){
	CMat6D invMat;

	invMat = *this;
	int r = invMat.inverseFastSelf();
	if (r ==  0) {
		string error("matriz have no inverse");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}

	return invMat;
}

SMF_INLINE_FORCED int CMat6D::getDimension() const {
	return 36;
}

SMF_INLINE_FORCED const float *CMat6D::toFloatPtr() const {
	return mat[0].toFloatPtr();
}

SMF_INLINE_FORCED float *CMat6D::toFloatPtr() {
	return mat[0].toFloatPtr();
}



/**
 * \class CMatXD
 *
 * \ingroup SMF_Math
 *
 * \if pt_br
 * \brief   Matriz densa de tamanho arbitrário
 * \note	A matriz é armazenada em memória alinhada a 16 bytes e com preenchimento (pad)
 * \warning XMatXD não pode ser utilizada por multiplos threads devido ao pool de memória temporário
 * \elseif us_en
 * \brief  arbitrary sized dense real matrix
 * \note  The matrix lives on 16 byte aligned and 16 byte padded memory.
 * \warning due to the temporary memory pool CMatXD cannot be used by multiple threads.
 *
 * \endif
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
class SMF_API CMatXD {
public:
					CMatXD();
					explicit CMatXD( int rows, int columns );
					explicit CMatXD( int rows, int columns, float *src );
					~CMatXD();

	void			set( int rows, int columns, const float *src );
	void			set( const CMat3D &m1, const CMat3D &m2 );
	void			set( const CMat3D &m1, const CMat3D &m2, const CMat3D &m3, const CMat3D &m4 );

	const float *	operator[]( int index ) const throw(Exception::CMathException);
	float *			operator[]( int index ) throw(Exception::CMathException);
	CMatXD &		operator=( const CMatXD &a );
	CMatXD			operator*( const float a ) const;
	CVecXD			operator*( const CVecXD &vec ) const throw(Exception::CMathException);
	CMatXD			operator*( const CMatXD &a ) const throw(Exception::CMathException);
	CMatXD			operator+( const CMatXD &a ) const throw(Exception::CMathException);
	CMatXD			operator-( const CMatXD &a ) const throw(Exception::CMathException);
	CMatXD &		operator*=( const float a );
	CMatXD &		operator*=( const CMatXD &a );
	CMatXD &		operator+=( const CMatXD &a ) throw(Exception::CMathException);
	CMatXD &		operator-=( const CMatXD &a ) throw(Exception::CMathException);

	friend CMatXD	operator*( const float a, const CMatXD &m );
	friend CVecXD	operator*( const CVecXD &vec, const CMatXD &m );
	friend CVecXD &	operator*=( CVecXD &vec, const CMatXD &m );
	/// exact compare, no epsilon
	bool			compare( const CMatXD &a ) const throw(Exception::CMathException);
	/// compare with epsilon
	bool			compare( const CMatXD &a, const float epsilon ) const throw(Exception::CMathException);
	/// exact compare, no epsilon
	bool			operator==( const CMatXD &a ) const;
	/// exact compare, no epsilon
	bool			operator!=( const CMatXD &a ) const;
	/// set the number of rows/columns
	void			setSize( int rows, int columns ) throw(Exception::CMathException);
	/// change the size keeping data intact where possible
	void			changeSize( int rows, int columns, bool makeZero = false );
	/// get the number of rows
	int				getNumRows() const { return numRows; }
	/// get the number of columns
	int				getNumColumns() const { return numColumns; }
	/// set float array pointer
	void			setData( int rows, int columns, float *data ) throw(Exception::CMathException);
	/// clear matrix
	void			toZero();
	/// set size and clear matrix
	void			zero( int rows, int columns );
	/// clear to identity matrix
	void			toIdentity() throw(Exception::CMathException);
	/// set size and clear to identity matrix
	void			identity( int rows, int columns ) throw(Exception::CMathException);
	/// create diagonal matrix from vector
	void			diag( const CVecXD &v );
	/// fill matrix with random values
	void			random( int seed, float l = 0.0f, float u = 1.0f );
	void			random( int rows, int columns, int seed, float l = 0.0f, float u = 1.0f );
	/// (*this) = - (*this)
	void			negate();
	/// clamp all values
	void			clamp( float min, float max );
	/// swap rows
	CMatXD &		swapRows( int r1, int r2 );
	/// swap columns
	CMatXD &		swapColumns( int r1, int r2 );
	/// swap rows and columns
	CMatXD &		swapRowsColumns( int r1, int r2 );
	/// remove a row
	CMatXD &		removeRow( int r );
	/// remove a column
	CMatXD &		removeColumn( int r );
	/// remove a row and column
	CMatXD &		removeRowColumn( int r );
	/// clear the upper triangle
	void			clearUpperTriangle() throw(Exception::CMathException);
	/// clear the lower triangle
	void			clearLowerTriangle() throw(Exception::CMathException);
	/// get square sub-matrix from 0,0 to size,size
	void			squareSubMatrix( const CMatXD &m, int size ) throw(Exception::CMathException);
	/// return maximum element difference between this and m
	float			maxDifference( const CMatXD &m ) const throw(Exception::CMathException);

	bool			isSquare() const { return ( numRows == numColumns ); }
	bool			isZero( const float epsilon = MATRIX_EPSILON ) const;
	bool			isIdentity( const float epsilon = MATRIX_EPSILON ) const throw(Exception::CMathException);
	bool			isDiagonal( const float epsilon = MATRIX_EPSILON ) const throw(Exception::CMathException);
	bool			isTriDiagonal( const float epsilon = MATRIX_EPSILON ) const;
	bool			isSymmetric( const float epsilon = MATRIX_EPSILON ) const;
	bool			isOrthogonal( const float epsilon = MATRIX_EPSILON ) const;
	bool			isOrthonormal( const float epsilon = MATRIX_EPSILON ) const;
	bool			isPMatrix( const float epsilon = MATRIX_EPSILON ) const;
	bool			isZMatrix( const float epsilon = MATRIX_EPSILON ) const;
	bool			isPositiveDefinite( const float epsilon = MATRIX_EPSILON ) const;
	bool			isSymmetricPositiveDefinite( const float epsilon = MATRIX_EPSILON ) const;
	bool			isPositiveSemiDefinite( const float epsilon = MATRIX_EPSILON ) const;
	bool			isSymmetricPositiveSemiDefinite( const float epsilon = MATRIX_EPSILON ) const;
	/// returns product of diagonal elements
	float			trace() const throw(Exception::CMathException);
	/// returns determinant of matrix
	float			determinant() const throw(Exception::CMathException);
	/// returns transpose
	CMatXD			transpose() const;
	/// transposes the matrix itself
	CMatXD &		transposeSelf();
	/**
	\brief interte a matriz  ( m * m.inverse() = identity)
	\return retorna a matriz inversa ou throw exception, caso não haja inversa
	**/
	CMatXD		inverse() const throw(Exception::CMathException) ;		// returns the inverse ( m * m.inverse() = identity )
	/**
	\brief calcula a matriz inversa ( m * m.inverse() = identity )
	\return retorna falso caso não haja inversa e verdadeiro caso haja
	**/
	bool			inverseSelf() throw(Exception::CMathException);		// returns false if determinant is zero
	/**
	\brief interte a matriz  ( m * m.inverse() = identity)
	\return retorna a matriz inversa ou throw exception, caso não haja inversa
	**/
	CMatXD			inverseFast() const throw(Exception::CMathException);	// returns the inverse ( m * m.inverse() = identity )
	/**
	\brief calcula a matriz inversa ( m * m.inverse() = identity )
	\return retorna falso caso não haja inversa e verdadeiro caso haja
	**/
	bool			inverseFastSelf() throw(Exception::CMathException);	// returns false if determinant is zero
	/// in-place inversion, returns false if determinant is zero
	bool			lowerTriangularInverse();
	/// in-place inversion, returns false if determinant is zero
	bool			upperTriangularInverse();
	/// (*this) * vec
	CVecXD			multiply( const CVecXD &vec ) const throw(Exception::CMathException);
	/// this->transpose() * vec
	CVecXD			transposeMultiply( const CVecXD &vec ) const throw(Exception::CMathException);
	/// (*this) * a
	CMatXD			multiply( const CMatXD &a ) const throw(Exception::CMathException);
	/// this->transpose() * a
	CMatXD			transposeMultiply( const CMatXD &a ) const throw(Exception::CMathException);
	/// dst = (*this) * vec
	void			multiply( CVecXD &dst, const CVecXD &vec ) const throw(Exception::CMathException);
	/// dst += (*this) * vec
	void			multiplyAdd( CVecXD &dst, const CVecXD &vec ) const;
	/// dst -= (*this) * vec
	void			multiplySub( CVecXD &dst, const CVecXD &vec ) const;
	/// dst = this->transpose() * vec
	void			transposeMultiply( CVecXD &dst, const CVecXD &vec ) const;
	/// dst += this->transpose() * vec
	void			TransposeMultiplyAdd( CVecXD &dst, const CVecXD &vec ) const;
	/// dst -= this->transpose() * vec
	void			TransposeMultiplySub( CVecXD &dst, const CVecXD &vec ) const;
	/// dst = (*this) * a
	void			multiply( CMatXD &dst, const CMatXD &a ) const throw(Exception::CMathException);
	/// dst = this->transpose() * a
	void			transposeMultiply( CMatXD &dst, const CMatXD &a ) const throw(Exception::CMathException);
	/// returns total number of values in matrix
	int				getDimension() const;
	/// interpret beginning of row as a const CVec6D
	const CVec6D &	subVec6( int row ) const throw(Exception::CMathException);
	/// interpret beginning of row as an CVec6D
	CVec6D &		subVec6( int row ) throw(Exception::CMathException);
	/// interpret complete row as a const CVecXD
	const CVecXD	subVecX( int row ) const throw(Exception::CMathException);
	/// interpret complete row as an CVecXD
	CVecXD			subVecX( int row ) throw(Exception::CMathException);
	/// pointer to const matrix float array
	const float *	toFloatPtr() const;
	/// pointer to matrix float array
	float *			toFloatPtr();
	const char *	toString( int precision = 2 ) const;

	void			update_RankOne( const CVecXD &v, const CVecXD &w, float alpha );
	void			update_RankOneSymmetric( const CVecXD &v, float alpha );
	void			update_RowColumn( const CVecXD &v, const CVecXD &w, int r );
	void			update_RowColumnSymmetric( const CVecXD &v, int r );
	void			update_Increment( const CVecXD &v, const CVecXD &w );
	void			update_IncrementSymmetric( const CVecXD &v );
	void			update_Decrement( int r );
	/// invert in-place with Gauss-Jordan elimination
	bool			inverse_GaussJordan();
	bool			inverse_UpdateRankOne( const CVecXD &v, const CVecXD &w, float alpha );
	bool			inverse_UpdateRowColumn( const CVecXD &v, const CVecXD &w, int r );
	bool			inverse_UpdateIncrement( const CVecXD &v, const CVecXD &w );
	bool			inverse_UpdateDecrement( const CVecXD &v, const CVecXD &w, int r );
	void			inverse_solve( CVecXD &x, const CVecXD &b ) const;
	/// factor in-place: L * U
	bool			lu_Factor( int *index, float *det = NULL );
	bool			lu_UpdateRankOne( const CVecXD &v, const CVecXD &w, float alpha, int *index );
	bool			lu_UpdateRowColumn( const CVecXD &v, const CVecXD &w, int r, int *index );
	bool			lu_UpdateIncrement( const CVecXD &v, const CVecXD &w, int *index );
	bool			lu_UpdateDecrement( const CVecXD &v, const CVecXD &w, const CVecXD &u, int r, int *index );
	void			lu_solve( CVecXD &x, const CVecXD &b, const int *index ) const;
	void			lu_Inverse( CMatXD &inv, const int *index ) const;
	void			lu_UnpackFactors( CMatXD &L, CMatXD &U ) const;
	void			lu_MultiplyFactors( CMatXD &m, const int *index ) const;
	/// factor in-place: Q * R
	bool			qr_Factor( CVecXD &c, CVecXD &d );
	bool			qr_UpdateRankOne( CMatXD &R, const CVecXD &v, const CVecXD &w, float alpha );
	bool			qr_UpdateRowColumn( CMatXD &R, const CVecXD &v, const CVecXD &w, int r );
	bool			qr_UpdateIncrement( CMatXD &R, const CVecXD &v, const CVecXD &w );
	bool			qr_UpdateDecrement( CMatXD &R, const CVecXD &v, const CVecXD &w, int r );
	void			qr_solve( CVecXD &x, const CVecXD &b, const CVecXD &c, const CVecXD &d ) const;
	void			qr_solve( CVecXD &x, const CVecXD &b, const CMatXD &R ) const;
	void			qr_Inverse( CMatXD &inv, const CVecXD &c, const CVecXD &d ) const;
	void			qr_UnpackFactors( CMatXD &Q, CMatXD &R, const CVecXD &c, const CVecXD &d ) const;
	void			qr_MultiplyFactors( CMatXD &m, const CVecXD &c, const CVecXD &d ) const;
	/// factor in-place: U * diag(w) * V.transpose()
	bool			svd_Factor( CVecXD &w, CMatXD &V );
	void			svd_solve( CVecXD &x, const CVecXD &b, const CVecXD &w, const CMatXD &V ) const;
	void			svd_Inverse( CMatXD &inv, const CVecXD &w, const CMatXD &V ) const;
	void			svd_MultiplyFactors( CMatXD &m, const CVecXD &w, const CMatXD &V ) const;
	/// factor in-place: L * L.transpose()
	bool			cholesky_Factor();
	bool			cholesky_UpdateRankOne( const CVecXD &v, float alpha, int offset = 0 );
	bool			cholesky_UpdateRowColumn( const CVecXD &v, int r );
	bool			cholesky_UpdateIncrement( const CVecXD &v );
	bool			cholesky_UpdateDecrement( const CVecXD &v, int r );
	/**
	\brief calcula o fator
	Matrix inversion techniques based on Cholesky decomposition and the related LDL decomposition
	are efficient techniques widely used for inversion of positive-definite/symmetric matrices across
	multiple fields.
	Existing matrix inversion algorithms based on Cholesky decomposition
	use either equation solving [3] or triangular matrix operations [4] with most efficient implementation
	requiring operations.
	\note  http://en.wikipedia.org/wiki/cholesky_decomposition
	\note  http://arxiv.org/ftp/arxiv/papers/1111/1111.4144.pdf
	**/
	void			cholesky_solve( CVecXD &x, const CVecXD &b ) const;
	void			cholesky_Inverse( CMatXD &inv ) const;
	void			cholesky_MultiplyFactors( CMatXD &m ) const;
	/// factor in-place: L * D * L.transpose()
	bool			ldlt_Factor();
	bool			ldlt_UpdateRankOne( const CVecXD &v, float alpha, int offset = 0 );
	bool			ldlt_UpdateRowColumn( const CVecXD &v, int r );
	bool			ldlt_UpdateIncrement( const CVecXD &v );
	bool			ldlt_UpdateDecrement( const CVecXD &v, int r );
	void			ldlt_solve( CVecXD &x, const CVecXD &b ) const;
	void			ldlt_Inverse( CMatXD &inv ) const;
	void			ldlt_UnpackFactors( CMatXD &L, CMatXD &D ) const;
	void			ldlt_MultiplyFactors( CMatXD &m ) const;

	void			triDiagonal_ClearTriangles();
	bool			triDiagonal_solve( CVecXD &x, const CVecXD &b ) const;
	void			triDiagonal_Inverse( CMatXD &inv ) const;

	bool			eigen_solveSymmetricTriDiagonal( CVecXD &eigenValues );
	bool			eigen_solveSymmetric( CVecXD &eigenValues );
	bool			eigen_solve( CVecXD &realEigenValues, CVecXD &imaginaryEigenValues );
	void			eigen_SortIncreasing( CVecXD &eigenValues );
	void			eigen_SortDecreasing( CVecXD &eigenValues );

	static void		Test();

private:
	/// number of rows
	int				numRows;
	/// number of columns
	int				numColumns;
	/// floats allocated, if -1 then mat points to data set with setData
	int				alloced;
	/// memory the matrix is stored
	float *			mat;
	/// used to store intermediate results
	static float	temp[MATX_MAX_TEMP+4];
	/// pointer to 16 byte aligned temporary memory
	static float *	tempPtr;
	/// index into memory pool, wraps around
	static int		tempIndex;
	bool			useSIMD;
private:
	void			setTempSize( int rows, int columns ) throw(Exception::CMathException);
	float			DeterminantGeneric() const;
	bool			InverseSelfGeneric();
	void			qr_Rotate( CMatXD &R, int i, float a, float b );
	float			Pythag( float a, float b ) const;
	void			svd_BiDiag( CVecXD &w, CVecXD &rv1, float &anorm );
	void			svd_InitialWV( CVecXD &w, CMatXD &V, CVecXD &rv1 );
	void			HouseholderReduction( CVecXD &diag, CVecXD &subd );
	bool			QL( CVecXD &diag, CVecXD &subd );
	void			HessenbergReduction( CMatXD &H );
	void			ComplexDivision( float xr, float xi, float yr, float yi, float &cdivr, float &cdivi );
	bool			HessenbergToRealSchur( CMatXD &H, CVecXD &realEigenValues, CVecXD &imaginaryEigenValues );
};

SMF_INLINE_FORCED CMatXD::CMatXD():useSIMD(true) {
	//useSIMD=CGlobalConfiguration::useSIMDprocessor.getBool();
	useSIMD=CSIMD::globalUseSIMD();
	if(useSIMD)  CSIMD::initHeap();
	numRows = numColumns = alloced = 0;
	mat = NULL;
}

SMF_INLINE_FORCED CMatXD::~CMatXD() {
	// if not temp memory
	if ( mat != NULL && ( mat < CMatXD::tempPtr || mat > CMatXD::tempPtr + MATX_MAX_TEMP ) && alloced != -1 ) {
		mem_Free16( mat );
	}
}

SMF_INLINE_FORCED CMatXD::CMatXD( int rows, int columns ):useSIMD(true) {
	useSIMD=CSIMD::globalUseSIMD();
	if(useSIMD)  CSIMD::initHeap();
	numRows = numColumns = alloced = 0;
	mat = NULL;
	setSize( rows, columns );
}

SMF_INLINE_FORCED CMatXD::CMatXD( int rows, int columns, float *src ):useSIMD(true) {
	useSIMD=CSIMD::globalUseSIMD();
	if(useSIMD)  CSIMD::initHeap();
	numRows = numColumns = alloced = 0;
	mat = NULL;
	setData( rows, columns, src );
}

SMF_INLINE_FORCED void CMatXD::set( int rows, int columns, const float *src ) {
	setSize( rows, columns );
	memcpy( this->mat, src, rows * columns * sizeof( float ) );
}

SMF_INLINE_FORCED void CMatXD::set( const CMat3D &m1, const CMat3D &m2 ) {
	int i, j;

	setSize( 3, 6 );
	for ( i = 0; i < 3; i++ ) {
		for ( j = 0; j < 3; j++ ) {
			mat[(i+0) * numColumns + (j+0)] = m1[i][j];
			mat[(i+0) * numColumns + (j+3)] = m2[i][j];
		}
	}
}

SMF_INLINE_FORCED void CMatXD::set( const CMat3D &m1, const CMat3D &m2, const CMat3D &m3, const CMat3D &m4 ) {
	int i, j;

	setSize( 6, 6 );
	for ( i = 0; i < 3; i++ ) {
		for ( j = 0; j < 3; j++ ) {
			mat[(i+0) * numColumns + (j+0)] = m1[i][j];
			mat[(i+0) * numColumns + (j+3)] = m2[i][j];
			mat[(i+3) * numColumns + (j+0)] = m3[i][j];
			mat[(i+3) * numColumns + (j+3)] = m4[i][j];
		}
	}
}

SMF_INLINE_FORCED const float *CMatXD::operator[]( int index ) const throw(Exception::CMathException){
	SMF_ASSERT( ( index >= 0 ) && ( index < numRows ) );
	if (( index < 0 ) || ( index >= numRows )) {
		string error("invalid index");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}

	return mat + index * numColumns;
}

SMF_INLINE_FORCED float *CMatXD::operator[]( int index ) throw(Exception::CMathException){
	SMF_ASSERT( ( index >= 0 ) && ( index < numRows ) );
	if (( index < 0 ) || ( index >= numRows )) {
		string error("invalid index");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}
	return mat + index * numColumns;
}

SMF_INLINE_FORCED CMatXD &CMatXD::operator=( const CMatXD &a ) {
	setSize( a.numRows, a.numColumns );
if(useSIMD){
	SIMDProcessor->copy16( mat, a.mat, a.numRows * a.numColumns );
}else{
	memcpy( mat, a.mat, a.numRows * a.numColumns * sizeof( float ) );
}
	CMatXD::tempIndex = 0;
	return *this;
}

SMF_INLINE_FORCED CMatXD CMatXD::operator*( const float a ) const {
	CMatXD m;

	m.setTempSize( numRows, numColumns );
if(useSIMD){
	SIMDProcessor->mul16( m.mat, mat, a, numRows * numColumns );
}else{
	int i, s;
	s = numRows * numColumns;
	for ( i = 0; i < s; i++ ) {
		m.mat[i] = mat[i] * a;
	}
}
	return m;
}

SMF_INLINE_FORCED CVecXD CMatXD::operator*( const CVecXD &vec ) const throw(Exception::CMathException){
	CVecXD dst;

	SMF_ASSERT( numColumns == vec.getSize() );
	if (numColumns != vec.getSize()) {
		string error("invalid index");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}

	dst.setTempSize( numRows );
if(useSIMD){
	SIMDProcessor->matX_MultiplyVecX( dst, *this, vec );
}else{
	multiply( dst, vec );
}
	return dst;
}

SMF_INLINE_FORCED CMatXD CMatXD::operator*( const CMatXD &a ) const throw(Exception::CMathException){
	CMatXD dst;

	SMF_ASSERT( numColumns == a.numRows );
	if (numColumns != a.numRows) {
		string error("invalid index");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}

	dst.setTempSize( numRows, a.numColumns );
if(useSIMD){
	SIMDProcessor->matX_MultiplyMatX( dst, *this, a );
}else{
	multiply( dst, a );
}
	return dst;
}

SMF_INLINE_FORCED CMatXD CMatXD::operator+( const CMatXD &a ) const throw(Exception::CMathException){
	CMatXD m;

	SMF_ASSERT( numRows == a.numRows && numColumns == a.numColumns );
	if (numColumns != a.numRows || numColumns != a.numColumns) {
		string error("invalid index");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}

	m.setTempSize( numRows, numColumns );
if(useSIMD){
	SIMDProcessor->add16( m.mat, mat, a.mat, numRows * numColumns );
}else{
	int i, s;
	s = numRows * numColumns;
	for ( i = 0; i < s; i++ ) {
		m.mat[i] = mat[i] + a.mat[i];
	}
}
	return m;
}

SMF_INLINE_FORCED CMatXD CMatXD::operator-( const CMatXD &a ) const throw(Exception::CMathException){
	CMatXD m;

	SMF_ASSERT( numRows == a.numRows && numColumns == a.numColumns );
	if (numColumns != a.numRows || numColumns != a.numColumns) {
		string error("invalid index");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}

	m.setTempSize( numRows, numColumns );
if(useSIMD){
	SIMDProcessor->sub16( m.mat, mat, a.mat, numRows * numColumns );
}else{
	int i, s;
	s = numRows * numColumns;
	for ( i = 0; i < s; i++ ) {
		m.mat[i] = mat[i] - a.mat[i];
	}
}
	return m;
}

SMF_INLINE_FORCED CMatXD &CMatXD::operator*=( const float a ) {
if(useSIMD){
	SIMDProcessor->mulAssign16( mat, a, numRows * numColumns );
}else{
	int i, s;
	s = numRows * numColumns;
	for ( i = 0; i < s; i++ ) {
		mat[i] *= a;
	}
}
	CMatXD::tempIndex = 0;
	return *this;
}

SMF_INLINE_FORCED CMatXD &CMatXD::operator*=( const CMatXD &a ) {
	*this = *this * a;
	CMatXD::tempIndex = 0;
	return *this;
}

SMF_INLINE_FORCED CMatXD &CMatXD::operator+=( const CMatXD &a ) throw(Exception::CMathException) {
	SMF_ASSERT( numRows == a.numRows && numColumns == a.numColumns );
	if (numColumns != a.numRows || numColumns != a.numColumns) {
		string error("invalid index");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}

if(useSIMD){
	SIMDProcessor->addAssign16( mat, a.mat, numRows * numColumns );
}else{
	int i, s;
	s = numRows * numColumns;
	for ( i = 0; i < s; i++ ) {
		mat[i] += a.mat[i];
	}
}
	CMatXD::tempIndex = 0;
	return *this;
}

SMF_INLINE_FORCED CMatXD &CMatXD::operator-=( const CMatXD &a ) throw(Exception::CMathException){
	SMF_ASSERT( numRows == a.numRows && numColumns == a.numColumns );
	if (numColumns != a.numRows || numColumns != a.numColumns) {
		string error("invalid index");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}

if(useSIMD){
	SIMDProcessor->subAssign16( mat, a.mat, numRows * numColumns );
}else{
	int i, s;
	s = numRows * numColumns;
	for ( i = 0; i < s; i++ ) {
		mat[i] -= a.mat[i];
	}
}
	CMatXD::tempIndex = 0;
	return *this;
}

SMF_INLINE_FORCED CMatXD operator*( const float a, CMatXD const &m ) {
	return m * a;
}

SMF_INLINE_FORCED CVecXD operator*( const CVecXD &vec, const CMatXD &m ) {
	return m * vec;
}

SMF_INLINE_FORCED CVecXD &operator*=( CVecXD &vec, const CMatXD &m ) {
	vec = m * vec;
	return vec;
}

SMF_INLINE_FORCED bool CMatXD::compare( const CMatXD &a ) const throw(Exception::CMathException){
	int i, s;

	SMF_ASSERT( numRows == a.numRows && numColumns == a.numColumns );
	if (!( numRows == a.numRows && numColumns == a.numColumns )) {
		string error("invalid index");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}

	s = numRows * numColumns;
	for ( i = 0; i < s; i++ ) {
		if ( mat[i] != a.mat[i] ) {
			return false;
		}
	}
	return true;
}

SMF_INLINE_FORCED bool CMatXD::compare( const CMatXD &a, const float epsilon ) const throw(Exception::CMathException){
	int i, s;

	SMF_ASSERT( numRows == a.numRows && numColumns == a.numColumns );
	if (!( numRows == a.numRows && numColumns == a.numColumns )) {
		string error("invalid index");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}

	s = numRows * numColumns;
	for ( i = 0; i < s; i++ ) {
		if ( CMath::fabs( mat[i] - a.mat[i] ) > epsilon ) {
			return false;
		}
	}
	return true;
}

SMF_INLINE_FORCED bool CMatXD::operator==( const CMatXD &a ) const {
	return compare( a );
}

SMF_INLINE_FORCED bool CMatXD::operator!=( const CMatXD &a ) const {
	return !compare( a );
}

SMF_INLINE_FORCED void CMatXD::setSize( int rows, int columns ) throw(Exception::CMathException){
	SMF_ASSERT( mat < CMatXD::tempPtr || mat > CMatXD::tempPtr + MATX_MAX_TEMP );
	if ( !(mat < CMatXD::tempPtr || mat > CMatXD::tempPtr + MATX_MAX_TEMP)) {
		string error("invalid size");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}

	int alloc = ( rows * columns + 3 ) & ~3;
	if ( alloc > alloced && alloced != -1 ) {
		if ( mat != NULL ) {
			mem_Free16( mat );
		}
		mat = (float *) mem_Alloc16( alloc * sizeof( float ) );
		alloced = alloc;
	}
	numRows = rows;
	numColumns = columns;
	MATX_CLEAREND();
}

SMF_INLINE_FORCED void CMatXD::setTempSize( int rows, int columns ) throw(Exception::CMathException){
	int newSize;

	newSize = ( rows * columns + 3 ) & ~3;
	SMF_ASSERT( newSize < MATX_MAX_TEMP );
	if ( !(newSize < MATX_MAX_TEMP)) {
		string error("invalid index");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}

	if ( CMatXD::tempIndex + newSize > MATX_MAX_TEMP ) {
		CMatXD::tempIndex = 0;
	}
	mat = CMatXD::tempPtr + CMatXD::tempIndex;
	CMatXD::tempIndex += newSize;
	alloced = newSize;
	numRows = rows;
	numColumns = columns;
	MATX_CLEAREND();
}

SMF_INLINE_FORCED void CMatXD::setData( int rows, int columns, float *data ) throw(Exception::CMathException){
	SMF_ASSERT( mat < CMatXD::tempPtr || mat > CMatXD::tempPtr + MATX_MAX_TEMP );
	if ( !( mat < CMatXD::tempPtr || mat > CMatXD::tempPtr + MATX_MAX_TEMP )) {
		string error("invalid index");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}


	if ( mat != NULL && alloced != -1 ) {
		mem_Free16( mat );
	}
#ifdef _WIN32
	SMF_ASSERT( ( ( (int) data ) & 15 ) == 0 ); // data must be 16 byte aligned
	if ( !( ( ( (int) data ) & 15 ) == 0 )) {
		string error("invalid alignment. data must be 16 byte aligned");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}

#else
	SMF_ASSERT( ( ( (uintptr_t) data ) & 15 ) == 0 ); // data must be 16 byte aligned
	if ( !( ( ( (uintptr_t) data ) & 15 ) == 0 )) {
		string error("invalid alignment. data must be 16 byte aligned");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}
#endif
	mat = data;
	alloced = -1;
	numRows = rows;
	numColumns = columns;
	MATX_CLEAREND();
}

SMF_INLINE_FORCED void CMatXD::toZero() {
if(useSIMD){
	SIMDProcessor->zero16( mat, numRows * numColumns );
}else{
	memset( mat, 0, numRows * numColumns * sizeof( float ) );
}
}

SMF_INLINE_FORCED void CMatXD::zero( int rows, int columns ) {
	setSize( rows, columns );
if(useSIMD){
	SIMDProcessor->zero16( mat, numRows * numColumns );
}else{
	memset( mat, 0, rows * columns * sizeof( float ) );
}
}

SMF_INLINE_FORCED void CMatXD::toIdentity() throw(Exception::CMathException){
	SMF_ASSERT( numRows == numColumns );
	if ( numRows != numColumns) {
		string error("numColumns != numRows");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}
if(useSIMD){
	SIMDProcessor->zero16( mat, numRows * numColumns );
}else{
	memset( mat, 0, numRows * numColumns * sizeof( float ) );
}
	for ( int i = 0; i < numRows; i++ ) {
		mat[i * numColumns + i] = 1.0f;
	}
}

SMF_INLINE_FORCED void CMatXD::identity( int rows, int columns ) throw(Exception::CMathException){
	SMF_ASSERT( rows == columns );
	if ( numRows != numColumns) {
		string error("numRows != numColums");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}
	setSize( rows, columns );
	CMatXD::toIdentity();
}

SMF_INLINE_FORCED void CMatXD::diag( const CVecXD &v ) {
	zero( v.getSize(), v.getSize() );
	for ( int i = 0; i < v.getSize(); i++ ) {
		mat[i * numColumns + i] = v[i];
	}
}

SMF_INLINE_FORCED void CMatXD::random( int seed, float l, float u ) {
	int i, s;
	float c;
	CRandom rnd(seed);

	c = u - l;
	s = numRows * numColumns;
	for ( i = 0; i < s; i++ ) {
		mat[i] = l + rnd.randomFloat() * c;
	}
}

SMF_INLINE_FORCED void CMatXD::random( int rows, int columns, int seed, float l, float u ) {
	int i, s;
	float c;
	CRandom rnd(seed);

	setSize( rows, columns );
	c = u - l;
	s = numRows * numColumns;
	for ( i = 0; i < s; i++ ) {
		mat[i] = l + rnd.randomFloat() * c;
	}
}

SMF_INLINE_FORCED void CMatXD::negate() {
if(useSIMD){
	SIMDProcessor->negate16( mat, numRows * numColumns );
}else{
	int i, s;
	s = numRows * numColumns;
	for ( i = 0; i < s; i++ ) {
		mat[i] = -mat[i];
	}
}
}

SMF_INLINE_FORCED void CMatXD::clamp( float min, float max ) {
	int i, s;
	s = numRows * numColumns;
	for ( i = 0; i < s; i++ ) {
		if ( mat[i] < min ) {
			mat[i] = min;
		} else if ( mat[i] > max ) {
			mat[i] = max;
		}
	}
}

SMF_INLINE_FORCED CMatXD &CMatXD::swapRows( int r1, int r2 ) {
	float *ptr;

	//
	//
#ifdef _WIN32
	ptr = (float *) _alloca16( numColumns * sizeof( float ) );
#else
	ptr = ((float *)((((uintptr_t)alloca( (numColumns * sizeof( float ))+15 )) + 15) & ~15));

#endif
	memcpy( ptr, mat + r1 * numColumns, numColumns * sizeof( float ) );
	memcpy( mat + r1 * numColumns, mat + r2 * numColumns, numColumns * sizeof( float ) );
	memcpy( mat + r2 * numColumns, ptr, numColumns * sizeof( float ) );

	return *this;
}

SMF_INLINE_FORCED CMatXD &CMatXD::swapColumns( int r1, int r2 ) {
	int i;
	float tmp, *ptr;

	for ( i = 0; i < numRows; i++ ) {
		ptr = mat + i * numColumns;
		tmp = ptr[r1];
		ptr[r1] = ptr[r2];
		ptr[r2] = tmp;
	}

	return *this;
}

SMF_INLINE_FORCED CMatXD &CMatXD::swapRowsColumns( int r1, int r2 ) {

	swapRows( r1, r2 );
	swapColumns( r1, r2 );
	return *this;
}

SMF_INLINE_FORCED void CMatXD::clearUpperTriangle() throw(Exception::CMathException){
	SMF_ASSERT( numRows == numColumns );
	if ( numRows != numColumns) {
		string error("numRows != numColums");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}
	for ( int i = numRows-2; i >= 0; i-- ) {
		memset( mat + i * numColumns + i + 1, 0, (numColumns - 1 - i) * sizeof(float) );
	}
}

SMF_INLINE_FORCED void CMatXD::clearLowerTriangle() throw(Exception::CMathException){
	SMF_ASSERT( numRows == numColumns );
	if ( numRows != numColumns) {
		string error("numcolumns != numRows");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}
	for ( int i = 1; i < numRows; i++ ) {
		memset( mat + i * numColumns, 0, i * sizeof(float) );
	}
}

SMF_INLINE_FORCED void CMatXD::squareSubMatrix( const CMatXD &m, int size ) throw(Exception::CMathException){
	int i;
	SMF_ASSERT( size <= m.numRows && size <= m.numColumns );
	if ( !( size <= m.numRows && size <= m.numColumns )) {
		string error("invalid size");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}
	setSize( size, size );
	for ( i = 0; i < size; i++ ) {
		memcpy( mat + i * numColumns, m.mat + i * m.numColumns, size * sizeof( float ) );
	}
}

SMF_INLINE_FORCED float CMatXD::maxDifference( const CMatXD &m ) const throw(Exception::CMathException){
	int i, j;
	float diff, maxDiff;

	SMF_ASSERT( numRows == m.numRows && numColumns == m.numColumns );
	if ( !( numRows == m.numRows && numColumns == m.numColumns )) {
		string error("invalid size");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}

	maxDiff = -1.0f;
	for ( i = 0; i < numRows; i++ ) {
		for ( j = 0; j < numColumns; j++ ) {
			diff = CMath::fabs( mat[ i * numColumns + j ] - m[i][j] );
			if ( maxDiff < 0.0f || diff > maxDiff ) {
				maxDiff = diff;
			}
		}
	}
	return maxDiff;
}

SMF_INLINE_FORCED bool CMatXD::isZero( const float epsilon ) const {
	// returns true if (*this) == zero
	for ( int i = 0; i < numRows; i++ ) {
		for ( int j = 0; j < numColumns; j++ ) {
			if ( CMath::fabs( mat[i * numColumns + j] ) > epsilon ) {
				return false;
			}
		}
	}
	return true;
}

SMF_INLINE_FORCED bool CMatXD::isIdentity( const float epsilon ) const throw(Exception::CMathException){
	// returns true if (*this) == identity
	SMF_ASSERT( numRows == numColumns );
	if ( !( numRows == numRows )) {
		string error("numColumns != NumRows");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}

	for ( int i = 0; i < numRows; i++ ) {
		for ( int j = 0; j < numColumns; j++ ) {
			if ( CMath::fabs( mat[i * numColumns + j] - (float)( i == j ) ) > epsilon ) {
				return false;
			}
		}
	}
	return true;
}

SMF_INLINE_FORCED bool CMatXD::isDiagonal( const float epsilon ) const throw(Exception::CMathException){
	// returns true if all elements are zero except for the elements on the diagonal
	SMF_ASSERT( numRows == numColumns );
	if ( !( numRows == numColumns )) {
		string error("invalid size");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}

	for ( int i = 0; i < numRows; i++ ) {
		for ( int j = 0; j < numColumns; j++ ) {
			if ( i != j && CMath::fabs( mat[i * numColumns + j] ) > epsilon ) {
				return false;
			}
		}
	}
	return true;
}

SMF_INLINE_FORCED bool CMatXD::isTriDiagonal( const float epsilon ) const {
	// returns true if all elements are zero except for the elements on the diagonal plus or minus one column

	if ( numRows != numColumns ) {
		return false;
	}
	for ( int i = 0; i < numRows-2; i++ ) {
		for ( int j = i+2; j < numColumns; j++ ) {
			if ( CMath::fabs( (*this)[i][j] ) > epsilon ) {
				return false;
			}
			if ( CMath::fabs( (*this)[j][i] ) > epsilon ) {
				return false;
			}
		}
	}
	return true;
}

SMF_INLINE_FORCED bool CMatXD::isSymmetric( const float epsilon ) const {
	// (*this)[i][j] == (*this)[j][i]
	if ( numRows != numColumns ) {
		return false;
	}
	for ( int i = 0; i < numRows; i++ ) {
		for ( int j = 0; j < numColumns; j++ ) {
			if ( CMath::fabs( mat[ i * numColumns + j ] - mat[ j * numColumns + i ] ) > epsilon ) {
				return false;
			}
		}
	}
	return true;
}

SMF_INLINE_FORCED float CMatXD::trace() const throw(Exception::CMathException){
	float trace = 0.0f;

	SMF_ASSERT( numRows == numColumns );
	if ( !( numRows == numColumns )) {
		string error("invalid size");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}

	// sum of elements on the diagonal
	for ( int i = 0; i < numRows; i++ ) {
		trace += mat[i * numRows + i];
	}
	return trace;
}

SMF_INLINE_FORCED float CMatXD::determinant() const throw(Exception::CMathException){

	SMF_ASSERT( numRows == numColumns );
	if ( !( numRows == numColumns )) {
		string error("invalid size");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}

	switch( numRows ) {
		case 1:
			return mat[0];
		case 2:
			return reinterpret_cast<const CMat2D *>(mat)->determinant();
		case 3:
			return reinterpret_cast<const CMat3D *>(mat)->determinant();
		case 4:
			return reinterpret_cast<const CMat4D *>(mat)->determinant();
		case 5:
			return reinterpret_cast<const CMat5D *>(mat)->determinant();
		case 6:
			return reinterpret_cast<const CMat6D *>(mat)->determinant();
		default:
			return DeterminantGeneric();
	}
	return 0.0f;
}

SMF_INLINE_FORCED CMatXD CMatXD::transpose() const {
	CMatXD transpose;
	int i, j;

	transpose.setTempSize( numColumns, numRows );

	for ( i = 0; i < numRows; i++ ) {
		for ( j = 0; j < numColumns; j++ ) {
			transpose.mat[j * transpose.numColumns + i] = mat[i * numColumns + j];
		}
	}

	return transpose;
}

SMF_INLINE_FORCED CMatXD &CMatXD::transposeSelf() {
	*this = transpose();
	return *this;
}

SMF_INLINE_FORCED CMatXD CMatXD::inverse() const throw(Exception::CMathException){
	CMatXD invMat;

	invMat.setTempSize( numRows, numColumns );
	memcpy( invMat.mat, mat, numRows * numColumns * sizeof( float ) );
	int r = invMat.inverseSelf();
	if (r ==  0) {
		string error("matriz have no inverse");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}
	return invMat;
}

SMF_INLINE_FORCED bool CMatXD::inverseSelf() throw(Exception::CMathException){

	SMF_ASSERT( numRows == numColumns );
	if (!(numRows == numColumns)) {
		string error("matriz have no inverse");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}

	switch( numRows ) {
		case 1:
			if ( CMath::fabs( mat[0] ) < MATRIX_INVERSE_EPSILON ) {
				return false;
			}
			mat[0] = 1.0f / mat[0];
			return true;
		case 2:
			return reinterpret_cast<CMat2D *>(mat)->inverseSelf();
		case 3:
			return reinterpret_cast<CMat3D *>(mat)->inverseSelf();
		case 4:
			return reinterpret_cast<CMat4D *>(mat)->inverseSelf();
		case 5:
			return reinterpret_cast<CMat5D *>(mat)->inverseSelf();
		case 6:
			return reinterpret_cast<CMat6D *>(mat)->inverseSelf();
		default:
			return InverseSelfGeneric();
	}
}

SMF_INLINE_FORCED CMatXD CMatXD::inverseFast() const throw(Exception::CMathException){
	CMatXD invMat;

	invMat.setTempSize( numRows, numColumns );
	memcpy( invMat.mat, mat, numRows * numColumns * sizeof( float ) );
	int r = invMat.inverseFastSelf();
	if (r ==  0) {
		string error("matriz have no inverse");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}

	return invMat;
}

SMF_INLINE_FORCED bool CMatXD::inverseFastSelf() throw(Exception::CMathException){

	SMF_ASSERT( numRows == numColumns );
	if (numRows != numColumns) {
		string error("matriz have no inverse");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}
	switch( numRows ) {
		case 1:
			if ( CMath::fabs( mat[0] ) < MATRIX_INVERSE_EPSILON ) {
				return false;
			}
			mat[0] = 1.0f / mat[0];
			return true;
		case 2:
			return reinterpret_cast<CMat2D *>(mat)->inverseFastSelf();
		case 3:
			return reinterpret_cast<CMat3D *>(mat)->inverseFastSelf();
		case 4:
			return reinterpret_cast<CMat4D *>(mat)->inverseFastSelf();
		case 5:
			return reinterpret_cast<CMat5D *>(mat)->inverseFastSelf();
		case 6:
			return reinterpret_cast<CMat6D *>(mat)->inverseFastSelf();
		default:
			return InverseSelfGeneric();
	}
	return false;
}

SMF_INLINE_FORCED CVecXD CMatXD::multiply( const CVecXD &vec ) const throw(Exception::CMathException){
	CVecXD dst;

	SMF_ASSERT( numColumns == vec.getSize() );
	if ( !( numColumns == vec.getSize() )) {
		string error("invalid size");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}

	dst.setTempSize( numRows );
if(useSIMD){
	SIMDProcessor->matX_MultiplyVecX( dst, *this, vec );
}else{
	multiply( dst, vec );
}
	return dst;
}

SMF_INLINE_FORCED CMatXD CMatXD::multiply( const CMatXD &a ) const throw(Exception::CMathException){
	CMatXD dst;

	SMF_ASSERT( numColumns == a.numRows );
	if ( !( numColumns == a.numRows )) {
		string error("invalid size");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}
	dst.setTempSize( numRows, a.numColumns );
if(useSIMD){
	SIMDProcessor->matX_MultiplyMatX( dst, *this, a );
}else{
	multiply( dst, a );
}
	return dst;
}

SMF_INLINE_FORCED CVecXD CMatXD::transposeMultiply( const CVecXD &vec ) const throw(Exception::CMathException){
	CVecXD dst;

	SMF_ASSERT( numRows == vec.getSize() );
	if ( !( numRows == vec.getSize() )) {
		string error("invalid numRows");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}

	dst.setTempSize( numColumns );
if(useSIMD){
	SIMDProcessor->matX_TransposeMultiplyVecX( dst, *this, vec );
}else{
	transposeMultiply( dst, vec );
}
	return dst;
}

SMF_INLINE_FORCED CMatXD CMatXD::transposeMultiply( const CMatXD &a ) const throw(Exception::CMathException){
	CMatXD dst;

	SMF_ASSERT( numRows == a.numRows );
	if ( !( numRows == a.numRows )) {
		string error("invalid size");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}

	dst.setTempSize( numColumns, a.numColumns );
if(useSIMD){
	SIMDProcessor->matX_TransposeMultiplyMatX( dst, *this, a );
}else{
	transposeMultiply( dst, a );
}
	return dst;
}

SMF_INLINE_FORCED void CMatXD::multiply( CVecXD &dst, const CVecXD &vec ) const throw(Exception::CMathException){
if(useSIMD){
	SIMDProcessor->matX_MultiplyVecX( dst, *this, vec );
}else{
	int i, j;
	const float *mPtr, *vPtr;
	float *dstPtr;

	mPtr = mat;
	vPtr = vec.toFloatPtr();
	dstPtr = dst.toFloatPtr();
	for ( i = 0; i < numRows; i++ ) {
		float sum = mPtr[0] * vPtr[0];
		for ( j = 1; j < numColumns; j++ ) {
			sum += mPtr[j] * vPtr[j];
		}
		dstPtr[i] = sum;
		mPtr += numColumns;
	}
}
}

SMF_INLINE_FORCED void CMatXD::multiplyAdd( CVecXD &dst, const CVecXD &vec ) const {
if(useSIMD){
	SIMDProcessor->matX_MultiplyAddVecX( dst, *this, vec );
}else{
	int i, j;
	const float *mPtr, *vPtr;
	float *dstPtr;

	mPtr = mat;
	vPtr = vec.toFloatPtr();
	dstPtr = dst.toFloatPtr();
	for ( i = 0; i < numRows; i++ ) {
		float sum = mPtr[0] * vPtr[0];
		for ( j = 1; j < numColumns; j++ ) {
			sum += mPtr[j] * vPtr[j];
		}
		dstPtr[i] += sum;
		mPtr += numColumns;
	}
}
}

SMF_INLINE_FORCED void CMatXD::multiplySub( CVecXD &dst, const CVecXD &vec ) const {
if(useSIMD){
	SIMDProcessor->matX_MultiplySubVecX( dst, *this, vec );
}else{
	int i, j;
	const float *mPtr, *vPtr;
	float *dstPtr;

	mPtr = mat;
	vPtr = vec.toFloatPtr();
	dstPtr = dst.toFloatPtr();
	for ( i = 0; i < numRows; i++ ) {
		float sum = mPtr[0] * vPtr[0];
		for ( j = 1; j < numColumns; j++ ) {
			sum += mPtr[j] * vPtr[j];
		}
		dstPtr[i] -= sum;
		mPtr += numColumns;
	}
}
}

SMF_INLINE_FORCED void CMatXD::transposeMultiply( CVecXD &dst, const CVecXD &vec ) const {
if(useSIMD){
	SIMDProcessor->matX_TransposeMultiplyVecX( dst, *this, vec );
}else{
	int i, j;
	const float *mPtr, *vPtr;
	float *dstPtr;

	vPtr = vec.toFloatPtr();
	dstPtr = dst.toFloatPtr();
	for ( i = 0; i < numColumns; i++ ) {
		mPtr = mat + i;
		float sum = mPtr[0] * vPtr[0];
		for ( j = 1; j < numRows; j++ ) {
			mPtr += numColumns;
			sum += mPtr[0] * vPtr[j];
		}
		dstPtr[i] = sum;
	}
}
}

SMF_INLINE_FORCED void CMatXD::TransposeMultiplyAdd( CVecXD &dst, const CVecXD &vec ) const {
if(useSIMD){
	SIMDProcessor->matX_TransposeMultiplyAddVecX( dst, *this, vec );
}else{
	int i, j;
	const float *mPtr, *vPtr;
	float *dstPtr;

	vPtr = vec.toFloatPtr();
	dstPtr = dst.toFloatPtr();
	for ( i = 0; i < numColumns; i++ ) {
		mPtr = mat + i;
		float sum = mPtr[0] * vPtr[0];
		for ( j = 1; j < numRows; j++ ) {
			mPtr += numColumns;
			sum += mPtr[0] * vPtr[j];
		}
		dstPtr[i] += sum;
	}
}
}

SMF_INLINE_FORCED void CMatXD::TransposeMultiplySub( CVecXD &dst, const CVecXD &vec ) const {
if(useSIMD){
	SIMDProcessor->matX_TransposeMultiplySubVecX( dst, *this, vec );
}else{
	int i, j;
	const float *mPtr, *vPtr;
	float *dstPtr;

	vPtr = vec.toFloatPtr();
	dstPtr = dst.toFloatPtr();
	for ( i = 0; i < numColumns; i++ ) {
		mPtr = mat + i;
		float sum = mPtr[0] * vPtr[0];
		for ( j = 1; j < numRows; j++ ) {
			mPtr += numColumns;
			sum += mPtr[0] * vPtr[j];
		}
		dstPtr[i] -= sum;
	}
}
}

SMF_INLINE_FORCED void CMatXD::multiply( CMatXD &dst, const CMatXD &a ) const throw(Exception::CMathException){
if(useSIMD){
	SIMDProcessor->matX_MultiplyMatX( dst, *this, a );
}else{
	int i, j, k, l, n;
	float *dstPtr;
	const float *m1Ptr, *m2Ptr;
	double sum;

	SMF_ASSERT( numColumns == a.numRows );
	if ( !( numColumns == a.numRows )) {
		string error("invalid size");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}

	dstPtr = dst.toFloatPtr();
	m1Ptr = toFloatPtr();
	m2Ptr = a.toFloatPtr();
	k = numRows;
	l = a.getNumColumns();

	for ( i = 0; i < k; i++ ) {
		for ( j = 0; j < l; j++ ) {
			m2Ptr = a.toFloatPtr() + j;
			sum = m1Ptr[0] * m2Ptr[0];
			for ( n = 1; n < numColumns; n++ ) {
				m2Ptr += l;
				sum += m1Ptr[n] * m2Ptr[0];
			}
			*dstPtr++ = sum;
		}
		m1Ptr += numColumns;
	}
}
}

SMF_INLINE_FORCED void CMatXD::transposeMultiply( CMatXD &dst, const CMatXD &a ) const throw(Exception::CMathException){
if(useSIMD){
	SIMDProcessor->matX_TransposeMultiplyMatX( dst, *this, a );
}else{
	int i, j, k, l, n;
	float *dstPtr;
	const float *m1Ptr, *m2Ptr;
	double sum;

	SMF_ASSERT( numRows == a.numRows );
	if ( !( numRows == a.numRows )) {
		string error("invalid size");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}

	dstPtr = dst.toFloatPtr();
	m1Ptr = toFloatPtr();
	k = numColumns;
	l = a.numColumns;

	for ( i = 0; i < k; i++ ) {
		for ( j = 0; j < l; j++ ) {
			m1Ptr = toFloatPtr() + i;
			m2Ptr = a.toFloatPtr() + j;
			sum = m1Ptr[0] * m2Ptr[0];
			for ( n = 1; n < numRows; n++ ) {
				m1Ptr += numColumns;
				m2Ptr += a.numColumns;
				sum += m1Ptr[0] * m2Ptr[0];
			}
			*dstPtr++ = sum;
		}
	}
}
}

SMF_INLINE_FORCED int CMatXD::getDimension() const {
	return numRows * numColumns;
}

SMF_INLINE_FORCED const CVec6D &CMatXD::subVec6( int row ) const throw(Exception::CMathException){
	SMF_ASSERT( numColumns >= 6 && row >= 0 && row < numRows );
	if ( !( numColumns >= 6 && row >= 0 && row < numRows )) {
		string error("invalid size");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}

	return *reinterpret_cast<const CVec6D *>(mat + row * numColumns);
}

SMF_INLINE_FORCED CVec6D &CMatXD::subVec6( int row ) throw(Exception::CMathException){
	SMF_ASSERT( numColumns >= 6 && row >= 0 && row < numRows );
	if ( !( numColumns >= 6 && row >= 0 && row < numRows )) {
		string error("invalid size");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}

	return *reinterpret_cast<CVec6D *>(mat + row * numColumns);
}

SMF_INLINE_FORCED const CVecXD CMatXD::subVecX( int row ) const throw(Exception::CMathException){
	CVecXD v;
	SMF_ASSERT( row >= 0 && row < numRows );
	if ( !( row >= 0 && row < numRows )) {
		string error("invalid size");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}
	v.setData( numColumns, mat + row * numColumns );
	return v;
}

SMF_INLINE_FORCED CVecXD CMatXD::subVecX( int row ) throw(Exception::CMathException){
	CVecXD v;
	SMF_ASSERT( row >= 0 && row < numRows );
	if ( !( row >= 0 && row < numRows )) {
		string error("invalid size");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}
	v.setData( numColumns, mat + row * numColumns );
	return v;
}

SMF_INLINE_FORCED const float *CMatXD::toFloatPtr() const {
	return mat;
}

SMF_INLINE_FORCED float *CMatXD::toFloatPtr() {
	return mat;
}
} //end MATH
} //end SMF
#endif /* !__MATH_MATRIX_H__ */
