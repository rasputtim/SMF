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

#ifndef _SMF__MATH_MATRIX_3D_H___
#define _SMF__MATH_MATRIX_3D_H___
#include "../SMF_Config.h"
#include "SMF_MathDefs.h"
#include "../exceptions/all.h"
#include "SMF_Vector.h"
#include "SMF_Simd.h"


namespace SMF{
namespace MATH {
/*
===============================================================================

  Matrix classes, all matrices are row-major except CMat3D

===============================================================================
*/

#define MATRIX_INVERSE_EPSILON		1e-14
#define MATRIX_EPSILON				1e-6

class CEulerAngles;
class CQuaternion;
class CCompQuaternion;
class CRotation;
class CMat4D;
class CMatJoint3x4;




/**
 * \class CMat3D
 *
 * \ingroup SMF_Math
 *
 * \if pt_br
 * \brief Matriz 3x3
 * \note	a matriz é column-major para manter compatibilidade com OpenGL e outros softwares de renderização
 * \elseif us_en
 * \brief 3x3 matrix
 * \note	matrix is column-major
 * \endif
 *
 * \see http://en.wikipedia.org/wiki/Row-major_order
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: Rasputtim@hotmail.com
 *
 */
class SMF_API CMat3D {
public:
					CMat3D();
					/**
					\if pt_br
					\brief Constrói uma Matriz através dos vetores das colunas
					\param x coluna 0
					\param y coluna 1
					\param z coluna 2
					\warning  \b IMPORTANTE: \b Esta matriz é column-major então os vetores representam as colunas da matriz
					\elseif us_en
					\brief
					\note
					\enif

					**/
					explicit CMat3D( const CVec3D &x, const CVec3D &y, const CVec3D &z );
					explicit CMat3D( const float xx, const float xy, const float xz, const float yx, const float yy, const float yz, const float zx, const float zy, const float zz );
					/**
					\if pt_br
					\brief Constrói uma Matriz através dos pontos das colunas
					\param xx,xy,xz = coluna 0
					\param yx,yy,yz = coluna 1
					\param zx,zy,zz - coluna 2
					\warning  \b IMPORTANTE: \b Esta matriz é column-major então os pontos devem estar ordenados corretamente.
					\elseif us_en
					\brief
					\note
					\enif

					**/
					explicit CMat3D( const float src[ 3 ][ 3 ] );

	/** 
	\brief Returns the given row. [noscript]
	 \param row The zero-based index [0, 2] of the row to get. 
	 */
	CVec3D &Row(int row);
	const CVec3D &Row(int row) const;

	CVec3D &Row3(int row) { return Row(row); }
	const CVec3D &Row3(int row) const { return Row(row); }

	/** 
	\brief Returns the given column.
	 \param col The zero-based index [0, 2] of the column to get. 
	 **/
	CONST_WIN32 CVec3D Col(int col) const;
	CONST_WIN32 CVec3D Col3(int col) const { return Col(col); }

	/** 
	\brief Returns the main diagonal.
	 The main diagonal consists of the elements at m[0][0], m[1][1], m[2][2]. 
	 **/
	CONST_WIN32 CVec3D diagonal() const;
	/** 
	\brief Sets the values of the given row.
	\param row The index of the row to set, in the range [0-2].
	\param data A pointer to an array of 3 floats that contain the new x, y and z values for the row. 
	*/
	void setRow(int row, const float *data);
	void setRow(int row, float x, float y, float z);
	void setRow(int row, const CVec3D &rowVector);

	/** 
	\brief Sets the values of the given column.
	\param column The index of the column to set, in the range [0-2].
	\param data A pointer to an array of 3 floats that contain the new x, y and z values for the column. 
	**/
	void setCol(int column, const float *data);
	void setCol(int column, float x, float y, float z);
	void setCol(int column, const CVec3D &columnVector);

	/// Scales the given row by a scalar.
	void scaleRow(int row, float scalar);

	/// Scales the given column by a scalar.
	void scaleCol(int col, float scalar);

	CMat3D			operator-() const;
	CMat3D			operator*( const float a ) const;
	CVec3D			operator*( const CVec3D &vec ) const;

	CMat3D			operator*( const CMat3D &a ) const;
		/// Transforms the given vector by this matrix (in the order M * v).
	/// This function ignores the w component of the given input vector. This component is assumed to be either 0 or 1.
	CVec4D operator *(const CVec4D &rhs) const;

	/// Converts the quaternion to a CMat3D and multiplies the two matrices together.
	CMat3D			operator *(const CQuaternion &quat) const;

	CMat3D mul(const CMat3D &rhs) const { return *this * rhs; };
	CMat3D mul(const CQuaternion &rhs) const { return *this * rhs; };
	CVec3D mul(const CVec3D &rhs) const { return *this * rhs; };
	CVec3D MulPos(const CVec3D &rhs) const { return mul(rhs); }
	CVec3D MulDir(const CVec3D &rhs) const { return mul(rhs); }

	CMat3D			operator+( const CMat3D &a ) const;
	CMat3D			operator-( const CMat3D &a ) const;
	CMat3D &		operator*=( const float a );
	CMat3D &		operator*=( const CMat3D &a );
	CMat3D &		operator+=( const CMat3D &a );
	CMat3D &		operator-=( const CMat3D &a );

	friend CMat3D	operator*( const float a, const CMat3D &mat );
	friend CVec3D	operator*( const CVec3D &vec, const CMat3D &mat );
	friend CVec3D &	operator*=( CVec3D &vec, const CMat3D &mat );

	bool			compare( const CMat3D &a ) const;						// exact compare, no epsilon
	bool			compare( const CMat3D &a, const float epsilon ) const;	// compare with epsilon
	bool			operator==( const CMat3D &a ) const;					// exact compare, no epsilon
	bool			operator!=( const CMat3D &a ) const;					// exact compare, no epsilon

	void			toZero();
	/** 
	\brief Tests if this is the identity matrix.
	\return Returns true if this matrix is the identity matrix, up to the given epsilon. 
	**/
	void			toIdentity();
	bool			isIdentity( const float epsilon = MATRIX_EPSILON ) const;
	/** 
	\brief Tests if this matrix is symmetric (M == M^T).
	 The test compares the elements for equality, up to the given epsilon. A matrix is symmetric if it is its own transpose. 
	**/
	bool			isSymmetric( const float epsilon = MATRIX_EPSILON ) const;
	bool			isDiagonal( const float epsilon = MATRIX_EPSILON ) const;
	bool			isRotated() const;
	/** 
	\brief Tests if this matrix is in lower triangular form.
	 \return Returns true if this matrix is in lower triangular form, up to the given epsilon. 
	 */
	bool isLowerTriangular(float epsilon =CMath::EPSILON_SuperLow) const;

	/**
	\brief Tests if this matrix is in upper triangular form.
	\return Returns true if this matrix is in upper triangular form, up to the given epsilon. 
	*/
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
	\brief Returns true if this matrix does not perform any scaling.
	A matrix does not do any scaling if the column vectors of this
	matrix are normalized in length, compared to the given epsilon. Note that this matrix may still perform
	reflection, i.e. it has a -1 scale along some axis.
	\note This function only examines the upper 3-by-3 part of this matrix.
	\note This function assumes that this matrix does not contain projection (the fourth row of this matrix is [0 0 0 1]). 
	*/
	bool hasUnitaryScale(float epsilonSq = 1e-6f) const;

	/** 
	\brief Returns true if this matrix performs a reflection along some plane.
	In 3D space, an even number of reflections corresponds to a rotation about some axis, so a matrix consisting of
	an odd number of consecutive mirror operations can only reflect about one axis. A matrix that contains reflection reverses
	the handedness of the coordinate system. This function tests if this matrix
	does perform mirroring. This occurs iff this matrix has a negative determinant. 
	*/
	bool hasNegativeScale() const;

	/**
	\brief Returns true if this matrix contains only uniform scaling, compared to the given epsilon.
	\note If the matrix does not really do any scaling, this function returns true (scaling uniformly by a factor of 1).
	**/
	bool hasUniformScale(float epsilon =CMath::EPSILON_SuperLow) const;

	/// Returns true if the row vectors of this matrix are all perpendicular to each other.
	bool isRowOrthogonal(float epsilon =CMath::EPSILON_SuperLow) const;

	/// Returns true if the column vectors of this matrix are all perpendicular to each other.
	bool isColOrthogonal(float epsilon =CMath::EPSILON_SuperLow) const;
	bool isColOrthogonal3(float epsilon =CMath::EPSILON_SuperLow) const { return isColOrthogonal(epsilon); }

	/**
	\brief Returns true if the column and row vectors of this matrix form an orthonormal set.
	\note In math terms, there does not exist such a thing as 'orthonormal matrix'. In math terms, a matrix
	is orthogonal iff its column and row vectors are orthogonal *unit* vectors.
	In the terms of this library however, a matrix is orthogonal iff its column and row vectors are orthogonal (no need to be unitary),
	and a matrix is orthonormal if the column and row vectors are orthonormal.
	**/
	bool isOrthonormal(float epsilon =CMath::EPSILON_SuperLow) const;


	/** 
	\brief Tests if this vector contains valid finite elements.
	 \see isNormalized(), isZero(), isPerpendicular(). 
	 */
	bool			isFinite() const;

	void			projectVector( const CVec3D &src, CVec3D &dst ) const;
	void			unprojectVector( const CVec3D &src, CVec3D &dst ) const;
	/// fix degenerate axial cases
	bool			FixDegeneracies();
	/// change tiny numbers to zero
	bool			FixDenormals();		

	float			trace() const;
	float			determinant() const;
	CMat3D			orthoNormalize() const;
	CMat3D &		orthonormalizeSelf();
	/// returns transpose
	CMat3D			transpose() const;	
	CMat3D &		transposeSelf();
	/**
	\brief interte a matriz  ( m * m.inverse() = identity)
	\return retorna a matriz inversa ou throw exception, caso não haja inversa
	**/
	CMat3D		inverse() const throw(Exception::CMathException) ;		// returns the inverse ( m * m.inverse() = identity )
	/**
	\brief calcula a matriz inversa ( m * m.inverse() = identity )
	\return retorna falso caso não haja inversa e verdadeiro caso haja
	\note calculo da inversa: http://www.infoescola.com/matematica/matriz-inversa-inversao-por-matriz-adjunta/
	\note calculo do cofator: http://www.mundoeducacao.com/matematica/calculando-cofator-uma-matriz.htm
	**/
	bool			inverseSelf();		// returns false if determinant is zero
	/**
	\brief interte a matriz  ( m * m.inverse() = identity)
	\return retorna a matriz inversa ou throw exception, caso não haja inversa
	**/
	CMat3D			inverseFast() const throw(Exception::CMathException);	// returns the inverse ( m * m.inverse() = identity )
	/**
	\brief calcula a matriz inversa ( m * m.inverse() = identity )
	\return retorna falso caso não haja inversa e verdadeiro caso haja
	**/
	bool			inverseFastSelf();	// returns false if determinant is zero
	CMat3D			transposeMultiply( const CMat3D &b ) const;

	CMat3D			InertiaTranslate( const float mass, const CVec3D &centerOfMass, const CVec3D &translation ) const;
	CMat3D &		InertiaTranslateSelf( const float mass, const CVec3D &centerOfMass, const CVec3D &translation );
	CMat3D			InertiaRotate( const CMat3D &rotation ) const;
	CMat3D &		InertiaRotateSelf( const CMat3D &rotation );



	int				getDimension() const;

	CEulerAngles		toAngles() const;
	CQuaternion			toQuat() const;
	CCompQuaternion		toCQuat() const;
	CRotation		toRotation() const;
	CMat4D			toMat4() const;
	CVec3D			toAngularVelocity() const;

	const float *	toFloatPtr() const;
	float *			toFloatPtr();
	const char *	toString( int precision = 2 ) const;

	friend void		transposeMultiply( const CMat3D &inv, const CMat3D &b, CMat3D &dst );
	friend CMat3D	skewSymmetric( CVec3D const &src );
	/// Transforms the given 3-vector by this matrix M, i.e. returns M * (x, y, z).
	CVec3D transform(const CVec3D &vector) const;
	CVec3D transform(float x, float y, float z) const;

	/**
	\brief Transforms the given 3-vector by this matrix M so that the vector occurs on the left-hand side, i.e.
	computes the product lhs * M. This is against the common convention used by this class when transforming
	geometrical objects, but this operation is still occasionally useful for other purposes.
	(Remember that M * v != v * M in general). 
	**/
	CVec3D transformLeft(const CVec3D &lhs) const;

	/**
	\brief Transforms the given 4-vector by this matrix M, i.e. returns M * (x, y, z, w).
	This function ignores the w component of the given input vector. This component is assumed to be either 0 or 1.
	**/
	CVec4D transform(const CVec4D &vector) const;

	/** 
	\brief Returns the scale components of this matrix.
	This function decomposes this matrix M into a form M = M' * S, where M' has unitary column vectors and S is a diagonal matrix.
	\return extractScale returns the diagonal entries of S, i.e. the scale of the columns of this matrix . If this matrix
	represents a local->world space transformation for an object, then this scale represents a 'local scale', i.e.
	scaling that is performed before translating and rotating the object from its local coordinate system to its world
	position.
	\note This function assumes that this matrix does not contain projection (the fourth row of this matrix is [0 0 0 1]).
	\note This function does not detect and return reflection (-1 scale along some axis).
	**/
	CVec3D extractScale() const;
	/**
	\brief solves the linear equation Ax=b.
	The matrix A in the equations is this matrix. 
	**/
	bool solveAxb(CVec3D b, CVec3D &x) const;

	/// Performs a batch transform of the given array.
	void batchTransform(CVec3D *pointArray, int numPoints) const;

	/// Performs a batch transform of the given array.
	void batchTransform(CVec3D *pointArray, int numPoints, int stride) const;

	/// Performs a batch transform of the given array.
	/// This function ignores the w component of the input vectors. These components are assumed to be either 0 or 1.
	void batchTransform(CVec4D *vectorArray, int numVectors) const;

	/// Performs a batch transform of the given array.
	/// This function ignores the w component of the input vectors. These components are assumed to be either 0 or 1.
	void batchTransform(CVec4D *vectorArray, int numVectors, int stride) const;
/**  to test the use of operator []
	friend class CMatJoint3x4;
	friend class CMat4D;
	friend class CMat5D;
	friend class CMat6D;
	friend class CMatXD;
	friend class CEulerAngles;
	friend class CRotation;
	friend class CVec3D;
	friend class CAABBox;
	friend class CQuaternion;
**/
	const CVec3D &	operator[]( int index ) const throw(Exception::CMathException); //use with care because this is column major ordered
	CVec3D &		operator[]( int index ) throw(Exception::CMathException);   //use with care because this is column major ordered

private:

	CVec3D			mat[ 3 ];
};

extern CMat3D mat3_zero;
extern CMat3D mat3_identity;
#define mat3_default	mat3_identity

SMF_INLINE_FORCED CMat3D::CMat3D() {
}

SMF_INLINE_FORCED CMat3D::CMat3D( const CVec3D &x, const CVec3D &y, const CVec3D &z ) {
	mat[ 0 ].x = x.x; mat[ 0 ].y = x.y; mat[ 0 ].z = x.z;
	mat[ 1 ].x = y.x; mat[ 1 ].y = y.y; mat[ 1 ].z = y.z;
	mat[ 2 ].x = z.x; mat[ 2 ].y = z.y; mat[ 2 ].z = z.z;
}

SMF_INLINE_FORCED CMat3D::CMat3D( const float xx, const float xy, const float xz, const float yx, const float yy, const float yz, const float zx, const float zy, const float zz ) {
	mat[ 0 ].x = xx; mat[ 0 ].y = xy; mat[ 0 ].z = xz;
	mat[ 1 ].x = yx; mat[ 1 ].y = yy; mat[ 1 ].z = yz;
	mat[ 2 ].x = zx; mat[ 2 ].y = zy; mat[ 2 ].z = zz;
}

SMF_INLINE_FORCED CMat3D::CMat3D( const float src[ 3 ][ 3 ] ) {
	memcpy( mat, src, 3 * 3 * sizeof( float ) );
}

SMF_INLINE_FORCED const CVec3D &CMat3D::operator[]( int index ) const throw(Exception::CMathException){
	SMF_ASSERT( ( index >= 0 ) && ( index < 3 ) );
	if (( index < 0 ) || ( index >= 3 )) {
		string error("invalid index");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}
	return mat[ index ];
}


SMF_INLINE_FORCED CVec3D &CMat3D::operator[]( int index )throw(Exception::CMathException) {
	SMF_ASSERT( ( index >= 0 ) && ( index < 3 ) );
	if (( index < 0 ) || ( index >= 3 )) {
		string error("invalid index");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}
	return mat[ index ];
}


SMF_INLINE_FORCED CVec3D &CMat3D::Row(int row)
{
	SMF_ASSERT(row >= 0);
	SMF_ASSERT(row < 3);
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (row < 0 || row >= 3)
		row = 0; // Benign failure, just give the first row.
#endif
	return reinterpret_cast<CVec3D &>(mat[row]);
}

SMF_INLINE_FORCED const CVec3D &CMat3D::Row(int row) const
{
	SMF_ASSERT(row >= 0);
	SMF_ASSERT(row < 3);
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (row < 0 || row >= 3)
		row = 0; // Benign failure, just give the first row.
#endif
	return CVec3D(mat[0][row], mat[1][row], mat[2][row]);
}

SMF_INLINE_FORCED CONST_WIN32 CVec3D CMat3D::Col(int col) const
{
	SMF_ASSERT(col >= 0);
	SMF_ASSERT(col < 3);
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (col < 0 || col >= 3)
		return CVec3D::nan;
#endif
	return reinterpret_cast<const CVec3D &>(mat[col]);

}

SMF_INLINE_FORCED CONST_WIN32 CVec3D CMat3D::diagonal() const
{
	return CVec3D(mat[0][0], mat[1][1], mat[2][2]);
}

SMF_INLINE_FORCED void CMat3D::setRow(int row, float x, float y, float z)
{
	SMF_ASSERT(row >= 0);
	SMF_ASSERT(row < 3);
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (row < 0 || row >= 3)
		return;
#endif
	SMF_ASSERT(MATH::isFinite(x));
	SMF_ASSERT(MATH::isFinite(y));
	SMF_ASSERT(MATH::isFinite(z));
	mat[0][row] = x;
	mat[1][row] = y;
	mat[2][row] = z;
}

SMF_INLINE_FORCED void CMat3D::setRow(int row, const CVec3D &rowVector)
{
	setRow(row, rowVector.x, rowVector.y, rowVector.z);
}

SMF_INLINE_FORCED void CMat3D::setRow(int row, const float *data)
{
	SMF_ASSERT(data);
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (!data)
		return;
#endif
	setRow(row, data[0], data[1], data[2]);
}

SMF_INLINE_FORCED void CMat3D::setCol(int column, float x, float y, float z)
{
	SMF_ASSERT(column >= 0);
	SMF_ASSERT(column < 3);
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (column < 0 || column >= 3)
		return;
#endif
	SMF_ASSERT(MATH::isFinite(x));
	SMF_ASSERT(MATH::isFinite(y));
	SMF_ASSERT(MATH::isFinite(z));
	mat[column][0] = x;
	mat[column][1] = y;
	mat[column][2] = z;
}

SMF_INLINE_FORCED void CMat3D::setCol(int column, const CVec3D &columnVector)
{
	setCol(column, columnVector.x, columnVector.y, columnVector.z);
}

SMF_INLINE_FORCED void CMat3D::setCol(int column, const float *data)
{
	SMF_ASSERT(data);
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (!data)
		return;
#endif
	setCol(column, data[0], data[1], data[2]);
}

SMF_INLINE_FORCED CMat3D CMat3D::operator-() const {
	return CMat3D(	-mat[0][0], -mat[0][1], -mat[0][2],
					-mat[1][0], -mat[1][1], -mat[1][2],
					-mat[2][0], -mat[2][1], -mat[2][2] );
}

SMF_INLINE_FORCED CVec3D CMat3D::operator*( const CVec3D &vec ) const {
	return CVec3D(
		mat[ 0 ].x * vec.x + mat[ 1 ].x * vec.y + mat[ 2 ].x * vec.z,
		mat[ 0 ].y * vec.x + mat[ 1 ].y * vec.y + mat[ 2 ].y * vec.z,
		mat[ 0 ].z * vec.x + mat[ 1 ].z * vec.y + mat[ 2 ].z * vec.z );
}
SMF_INLINE_FORCED CVec4D CMat3D::operator *(const CVec4D &rhs) const
{
	return CVec4D((Row(0)* rhs.toVec3()),
	              (Row(1)* rhs.toVec3()),
	              (Row(2)* rhs.toVec3()),
	              rhs.w);
}

SMF_INLINE_FORCED CMat3D CMat3D::operator*( const CMat3D &a ) const {
	SMF_ASSERT(isFinite());
	SMF_ASSERT(a.isFinite());
    CMat3D Mat1 = *this;
	CMat3D Mat2 = a;

	int i, j;
	float *m1Ptr, *m2Ptr;
	float *dstPtr;
	CMat3D dst;

	m1Ptr = Mat1.toFloatPtr();
	m2Ptr = Mat2.toFloatPtr();

	float temp = 0;
	float one,two;
   int b, c;

   for(i = 0; i < 3; i++)
   {
       for(b = 0; b < 3; b++)
       {
           for(c = 0; c < 3; c++)
           {
               one = Mat1[c][b];
			   two = Mat2[i][c];
			   temp += one * two;
           }
           dst[i][b] = temp;
           temp = 0.0f;
       }
   }

	return dst;
#if 0
	   int i,j;

   for (i = 0;i < 3;i++)
      for (j = 0;j < 3;j++)
	 Result[i*3+j] = Mat1[i*3+0]*Mat2[0*3+j] +
	    Mat1[i*3+1]*Mat2[1*3+j] +
	    Mat1[i*3+2]*Mat2[2*3+j];
#endif
}


SMF_INLINE_FORCED CMat3D CMat3D::operator*( const float a ) const {
	return CMat3D(
		mat[0].x * a, mat[0].y * a, mat[0].z * a,
		mat[1].x * a, mat[1].y * a, mat[1].z * a,
		mat[2].x * a, mat[2].y * a, mat[2].z * a );
}

SMF_INLINE_FORCED CMat3D CMat3D::operator+( const CMat3D &a ) const {
	return CMat3D(
		mat[0].x + a[0].x, mat[0].y + a[0].y, mat[0].z + a[0].z,
		mat[1].x + a[1].x, mat[1].y + a[1].y, mat[1].z + a[1].z,
		mat[2].x + a[2].x, mat[2].y + a[2].y, mat[2].z + a[2].z );
}

SMF_INLINE_FORCED CMat3D CMat3D::operator-( const CMat3D &a ) const {
	return CMat3D(
		mat[0].x - a[0].x, mat[0].y - a[0].y, mat[0].z - a[0].z,
		mat[1].x - a[1].x, mat[1].y - a[1].y, mat[1].z - a[1].z,
		mat[2].x - a[2].x, mat[2].y - a[2].y, mat[2].z - a[2].z );
}

SMF_INLINE_FORCED CMat3D &CMat3D::operator*=( const float a ) {
	mat[0].x *= a; mat[0].y *= a; mat[0].z *= a;
	mat[1].x *= a; mat[1].y *= a; mat[1].z *= a;
	mat[2].x *= a; mat[2].y *= a; mat[2].z *= a;

    return *this;
}

SMF_INLINE_FORCED CMat3D &CMat3D::operator*=( const CMat3D &a ) {
	int i, j;
	const float *m2Ptr;
	float *m1Ptr, dst[3];

	m1Ptr = reinterpret_cast<float *>(this);
	m2Ptr = reinterpret_cast<const float *>(&a);

	for ( i = 0; i < 3; i++ ) {
		for ( j = 0; j < 3; j++ ) {
			dst[j]  = m1Ptr[0] * m2Ptr[ 0 * 3 + j ]
					+ m1Ptr[1] * m2Ptr[ 1 * 3 + j ]
					+ m1Ptr[2] * m2Ptr[ 2 * 3 + j ];
		}
		m1Ptr[0] = dst[0]; m1Ptr[1] = dst[1]; m1Ptr[2] = dst[2];
		m1Ptr += 3;
	}
	return *this;
}

SMF_INLINE_FORCED CMat3D &CMat3D::operator+=( const CMat3D &a ) {
	mat[0].x += a[0].x; mat[0].y += a[0].y; mat[0].z += a[0].z;
	mat[1].x += a[1].x; mat[1].y += a[1].y; mat[1].z += a[1].z;
	mat[2].x += a[2].x; mat[2].y += a[2].y; mat[2].z += a[2].z;

    return *this;
}

SMF_INLINE_FORCED CMat3D &CMat3D::operator-=( const CMat3D &a ) {
	mat[0].x -= a[0].x; mat[0].y -= a[0].y; mat[0].z -= a[0].z;
	mat[1].x -= a[1].x; mat[1].y -= a[1].y; mat[1].z -= a[1].z;
	mat[2].x -= a[2].x; mat[2].y -= a[2].y; mat[2].z -= a[2].z;

    return *this;
}

SMF_INLINE_FORCED CVec3D operator*( const CVec3D &vec, const CMat3D &mat ) {
	return mat * vec;
}

SMF_INLINE_FORCED CMat3D operator*( const float a, const CMat3D &mat ) {
	return mat * a;
}

SMF_INLINE_FORCED CVec3D &operator*=( CVec3D &vec, const CMat3D &mat ) {
	float x = mat[ 0 ].x * vec.x + mat[ 1 ].x * vec.y + mat[ 2 ].x * vec.z;
	float y = mat[ 0 ].y * vec.x + mat[ 1 ].y * vec.y + mat[ 2 ].y * vec.z;
	vec.z = mat[ 0 ].z * vec.x + mat[ 1 ].z * vec.y + mat[ 2 ].z * vec.z;
	vec.x = x;
	vec.y = y;
	return vec;
}

SMF_INLINE_FORCED bool CMat3D::compare( const CMat3D &a ) const {
	if ( mat[0].compare( a[0] ) &&
		mat[1].compare( a[1] ) &&
		mat[2].compare( a[2] ) ) {
		return true;
	}
	return false;
}
SMF_INLINE_FORCED CVec3D CMat3D::extractScale() const
{
	return CVec3D(Col(0).getLenght(), Col(1).getLenght(), Col(2).getLenght());
}
SMF_INLINE_FORCED bool CMat3D::compare( const CMat3D &a, const float epsilon ) const {
	if ( mat[0].compare( a[0], epsilon ) &&
		mat[1].compare( a[1], epsilon ) &&
		mat[2].compare( a[2], epsilon ) ) {
		return true;
	}
	return false;
}

SMF_INLINE_FORCED bool CMat3D::operator==( const CMat3D &a ) const {
	return compare( a );
}

SMF_INLINE_FORCED bool CMat3D::operator!=( const CMat3D &a ) const {
	return !compare( a );
}

SMF_INLINE_FORCED void CMat3D::toZero() {
	memset( mat, 0, sizeof( CMat3D ) );
}

SMF_INLINE_FORCED void CMat3D::toIdentity() {
	*this = mat3_identity;
}

SMF_INLINE_FORCED bool CMat3D::isIdentity( const float epsilon ) const {
	return compare( mat3_identity, epsilon );
}
SMF_INLINE_FORCED bool CMat3D::isFinite() const
{
	for(int y = 0; y < 3; ++y)
		for(int x = 0; x < 3; ++x)
			if (!MATH::isFinite(mat[y][x]))
				return false;
	return true;
}
SMF_INLINE_FORCED bool CMat3D::isSymmetric( const float epsilon ) const {
	if ( CMath::fabs( mat[1][0] - mat[0][1] ) > epsilon ) {
		return false;
	}
	if ( CMath::fabs( mat[2][0] - mat[0][2] ) > epsilon ) {
		return false;
	}
	if ( CMath::fabs( mat[2][1] - mat[1][2] ) > epsilon ) {
		return false;
	}
	return true;
}

SMF_INLINE_FORCED bool CMat3D::isDiagonal( const float epsilon ) const {
	if ( CMath::fabs( mat[1][0] ) > epsilon ||
		CMath::fabs( mat[2][0] ) > epsilon ||
		CMath::fabs( mat[0][1] ) > epsilon ||
		CMath::fabs( mat[2][1] ) > epsilon ||
		CMath::fabs( mat[0][2] ) > epsilon ||
		CMath::fabs( mat[1][2] ) > epsilon ) {
		return false;
	}
	return true;
}
SMF_INLINE_FORCED bool CMat3D::isLowerTriangular(float epsilon) const
{
	return CMath::equalsAbs(mat[1][0], 0.f, epsilon)
	    && CMath::equalsAbs(mat[2][0], 0.f, epsilon)
	    && CMath::equalsAbs(mat[2][1], 0.f, epsilon);
}

SMF_INLINE_FORCED bool CMat3D::isUpperTriangular(float epsilon) const
{
	return CMath::equalsAbs(mat[0][1], 0.f, epsilon)
	    && CMath::equalsAbs(mat[0][2], 0.f, epsilon)
	    && CMath::equalsAbs(mat[1][2], 0.f, epsilon);
}

SMF_INLINE_FORCED bool CMat3D::isInvertible(float epsilon) const
{
	float d = determinant();
	bool isSingular = CMath::equalsAbs(d, 0.f, epsilon);
	//SMF_ASSERT(CMat3D(*this).inverse(epsilon) != isSingular); // isInvertible() and Inverse() must match!
	return !isSingular;
}


SMF_INLINE_FORCED bool CMat3D::isSkewSymmetric(float epsilon) const
{
	return CMath::equalsAbs(mat[0][0], 0.f, epsilon) &&
		CMath::equalsAbs(mat[1][1], 0.f, epsilon) &&
		CMath::equalsAbs(mat[2][2], 0.f, epsilon) &&
		CMath::equalsAbs(mat[1][0], -mat[1][0], epsilon) &&
		CMath::equalsAbs(mat[2][0], -mat[2][0], epsilon) &&
		CMath::equalsAbs(mat[2][1], -mat[2][1], epsilon);
}

SMF_INLINE_FORCED bool CMat3D::hasUnitaryScale(float epsilon) const
{
	CVec3D scale = extractScale();
	return scale.compare(CVec3D(1.f, 1.f, 1.f), epsilon);
}

SMF_INLINE_FORCED bool CMat3D::hasNegativeScale() const
{
	return determinant() < 0.f;
}

SMF_INLINE_FORCED bool CMat3D::hasUniformScale(float epsilon) const
{
	CVec3D scale = extractScale();
	return CMath::equalsAbs(scale.x, scale.y, epsilon) && CMath::equalsAbs(scale.x, scale.z, epsilon);
}

SMF_INLINE_FORCED bool CMat3D::isRowOrthogonal(float epsilon) const
{
	return Row(0).isPerpendicular(Row(1), epsilon)
	    && Row(0).isPerpendicular(Row(2), epsilon)
	    && Row(1).isPerpendicular(Row(2), epsilon);
}

SMF_INLINE_FORCED bool CMat3D::isColOrthogonal(float epsilon) const
{
	return Col(0).isPerpendicular(Col(1), epsilon)
	    && Col(0).isPerpendicular(Col(2), epsilon)
	    && Col(1).isPerpendicular(Col(2), epsilon);
}

SMF_INLINE_FORCED bool CMat3D::isOrthonormal(float epsilon) const
{
	///\todo Epsilon magnitudes don't match.
	return isColOrthogonal(epsilon) && Row(0).isNormalized(epsilon) && Row(1).isNormalized(epsilon) && Row(2).isNormalized(epsilon);
}

SMF_INLINE_FORCED bool CMat3D::isRotated() const {
	return !compare( mat3_identity );
}

SMF_INLINE_FORCED void CMat3D::projectVector( const CVec3D &src, CVec3D &dst ) const {
	dst.x = src * mat[ 0 ];
	dst.y = src * mat[ 1 ];
	dst.z = src * mat[ 2 ];
}

SMF_INLINE_FORCED void CMat3D::unprojectVector( const CVec3D &src, CVec3D &dst ) const {
	dst = mat[ 0 ] * src.x + mat[ 1 ] * src.y + mat[ 2 ] * src.z;
}

SMF_INLINE_FORCED bool CMat3D::FixDegeneracies() {
	bool r = mat[0].fixDegenerateNormal();
	r |= mat[1].fixDegenerateNormal();
	r |= mat[2].fixDegenerateNormal();
	return r;
}

SMF_INLINE_FORCED bool CMat3D::FixDenormals() {
	bool r = mat[0].fixDenormals();
	r |= mat[1].fixDenormals();
	r |= mat[2].fixDenormals();
	return r;
}

SMF_INLINE_FORCED float CMat3D::trace() const {
	return ( mat[0][0] + mat[1][1] + mat[2][2] );
}

SMF_INLINE_FORCED CMat3D CMat3D::orthoNormalize() const {
	CMat3D ortho;

	ortho = *this;
	ortho[ 0 ].toNormal();
	ortho[ 2 ].cross( mat[ 0 ], mat[ 1 ] );
	ortho[ 2 ].toNormal();
	ortho[ 1 ].cross( mat[ 2 ], mat[ 0 ] );
	ortho[ 1 ].toNormal();
	return ortho;
}

SMF_INLINE_FORCED CMat3D &CMat3D::orthonormalizeSelf() {
	mat[ 0 ].toNormal();
	mat[ 2 ].cross( mat[ 0 ], mat[ 1 ] );
	mat[ 2 ].toNormal();
	mat[ 1 ].cross( mat[ 2 ], mat[ 0 ] );
	mat[ 1 ].toNormal();
	return *this;
}

SMF_INLINE_FORCED CMat3D CMat3D::transpose() const {
	return CMat3D(	mat[0][0], mat[1][0], mat[2][0],
					mat[0][1], mat[1][1], mat[2][1],
					mat[0][2], mat[1][2], mat[2][2] );
}

SMF_INLINE_FORCED CMat3D &CMat3D::transposeSelf() {
	float tmp0, tmp1, tmp2;

	tmp0 = mat[0][1];
	mat[0][1] = mat[1][0];
	mat[1][0] = tmp0;
	tmp1 = mat[0][2];
	mat[0][2] = mat[2][0];
	mat[2][0] = tmp1;
	tmp2 = mat[1][2];
	mat[1][2] = mat[2][1];
	mat[2][1] = tmp2;

	return *this;
}

SMF_INLINE_FORCED CMat3D CMat3D::inverse() const throw(Exception::CMathException){
	CMat3D invMat;

	invMat = *this;
	int r = invMat.inverseSelf();
	if (r ==  0) {
		string error("matriz have no inverse");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}

	return invMat;
}

SMF_INLINE_FORCED CMat3D CMat3D::inverseFast() const throw(Exception::CMathException){
	CMat3D invMat;

	invMat = *this;
	int r = invMat.inverseFastSelf();
	if (r ==  0) {
		string error("matriz have no inverse");
		throw Exception::CMathException(__FUNCTION__,__LINE__,error);
	}

	return invMat;
}

SMF_INLINE_FORCED CMat3D CMat3D::transposeMultiply( const CMat3D &b ) const {
	return CMat3D(	mat[0].x * b[0].x + mat[1].x * b[1].x + mat[2].x * b[2].x,
					mat[0].x * b[0].y + mat[1].x * b[1].y + mat[2].x * b[2].y,
					mat[0].x * b[0].z + mat[1].x * b[1].z + mat[2].x * b[2].z,
					mat[0].y * b[0].x + mat[1].y * b[1].x + mat[2].y * b[2].x,
					mat[0].y * b[0].y + mat[1].y * b[1].y + mat[2].y * b[2].y,
					mat[0].y * b[0].z + mat[1].y * b[1].z + mat[2].y * b[2].z,
					mat[0].z * b[0].x + mat[1].z * b[1].x + mat[2].z * b[2].x,
					mat[0].z * b[0].y + mat[1].z * b[1].y + mat[2].z * b[2].y,
					mat[0].z * b[0].z + mat[1].z * b[1].z + mat[2].z * b[2].z );
}

SMF_INLINE_FORCED void transposeMultiply( const CMat3D &transpose, const CMat3D &b, CMat3D &dst ) {
	dst[0].x = transpose[0].x * b[0].x + transpose[1].x * b[1].x + transpose[2].x * b[2].x;
	dst[0].y = transpose[0].x * b[0].y + transpose[1].x * b[1].y + transpose[2].x * b[2].y;
	dst[0].z = transpose[0].x * b[0].z + transpose[1].x * b[1].z + transpose[2].x * b[2].z;
	dst[1].x = transpose[0].y * b[0].x + transpose[1].y * b[1].x + transpose[2].y * b[2].x;
	dst[1].y = transpose[0].y * b[0].y + transpose[1].y * b[1].y + transpose[2].y * b[2].y;
	dst[1].z = transpose[0].y * b[0].z + transpose[1].y * b[1].z + transpose[2].y * b[2].z;
	dst[2].x = transpose[0].z * b[0].x + transpose[1].z * b[1].x + transpose[2].z * b[2].x;
	dst[2].y = transpose[0].z * b[0].y + transpose[1].z * b[1].y + transpose[2].z * b[2].y;
	dst[2].z = transpose[0].z * b[0].z + transpose[1].z * b[1].z + transpose[2].z * b[2].z;
}

SMF_INLINE_FORCED CMat3D skewSymmetric( CVec3D const &src ) {
	return CMat3D( 0.0f, -src.z,  src.y, src.z,   0.0f, -src.x, -src.y,  src.x,   0.0f );
}

SMF_INLINE_FORCED int CMat3D::getDimension() const {
	return 9;
}

SMF_INLINE_FORCED const float *CMat3D::toFloatPtr() const {
	return mat[0].toFloatPtr();
}

SMF_INLINE_FORCED float *CMat3D::toFloatPtr() {
	return mat[0].toFloatPtr();
}

} //end MATH
} //end SMF
#endif /* !__MATH_MATRIX_3D_H__ */
