/*
  SMF -  Super Math Fabric  
  Copyright (C) 2014 Rasputtim <carygrant@ig.com.br>

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

#ifndef _SMF__JOINT_TRANSFORM_H__
#define _SMF__JOINT_TRANSFORM_H__
#include "../math/SMF_Quaternion.h"
#include "../math/SMF_Vector.h"
#include "../math/SMF_Matriz.h"


namespace SMF {

namespace MATH {

class CScaleOP;
class CTranslateOp;
class CScaleOp;
class CTranslateOp;
	/**
 * \class CJointQuaternion
 *
 * \ingroup SMF_Math
 *
 * \if pt_br
 * \brief Classe utilizada na trnasformação de juntas com Quaternions (Joint transforms)
 * \elseif us_en
 *
 * \brief Used by Joint transforms with Quaternions
 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: carygrant@ig.com.br
 *
 */
class  SMF_API CJointQuaternion {
public:

	MATH::CQuaternion		q;
	MATH::CVec3D			t;
};


/**
 * \class CMatJoint3x4
 *
 * \ingroup SMF_Math
 *
 * \if pt_br
 * \brief Matriz 3x4 utilizada nas transformações de juntas (Joint transforms)
  \p MATH::CMat3D m;
  \p MATH::CVec3D t;

  \p m[0][0], m[1][0], m[2][0], t[0]
  \p m[0][1], m[1][1], m[2][1], t[1]
  \p m[0][2], m[1][2], m[2][2], t[2]
  \elseif us_en
 *
 * \brief  3x4 Matrix used for Joint transforms
  \p MATH::CMat3D m;
  \p MATH::CVec3D t;

  \p m[0][0], m[1][0], m[2][0], t[0]
  \p m[0][1], m[1][1], m[2][1], t[1]
  \p m[0][2], m[1][2], m[2][2], t[2]


 *
 * \author (last to touch it) $Autor: Rasputtim $
 *
 * \version 1.0 $Revision: 1.0 $
 *
 * Contact: carygrant@ig.com.br
 *
 */
class  SMF_API CMatJoint3x4 {
public:
	friend class CSIMD_Generic;
	/// Creates a new CMatJoint3x4 with uninitialized member values.
	CMatJoint3x4(){};
	/// Constructs a new CMatJoint3x4 by explicitly specifying all the matrix elements.
	/// The elements are specified in row-major format, i.e. the first row first followed by the second and third row.
	/// E.g. The element _10 denotes the scalar at second (index 1) row, first (index 0) column.
	explicit CMatJoint3x4(float _00, float _01, float _02, float _03,
			 float _10, float _11, float _12, float _13,
			 float _20, float _21, float _22, float _23);

	/// Constructs this CMatJoint3x4 from the given quaternion.
	explicit CMatJoint3x4(const CQuaternion &orientation);

	/// Constructs this CMatJoint3x4 from the given quaternion and translation.
	/// Logically, the translation occurs after the rotation has been performed.
	CMatJoint3x4(const CQuaternion &orientation, const CVec3D &translation);



	/** Returns the given row. [noscript]
	\param row The zero-based index [0, 2] of the row to get.
	**/
	CVec4D &Row(int row);
	const CVec4D &Row(int row) const;

	/// Returns the three first elements of the given row. [noscript]
	/** \param row The zero-based index [0, 2] of the row to get. */
	CVec3D &Row3(int row);
	const CVec3D &Row3(int row) const;

	/// Returns the given column.
	/** \param col The zero-based index [0, 3] of the column to get. */
	CONST_WIN32 CVec3D Col(int col) const;
	CONST_WIN32 CVec3D Col3(int col) const { return Col(col); }

	/// treat the last line as identity
	CONST_WIN32 CVec4D Col4(int col) const;
	/// Returns the main diagonal.
	/** The main diagonal consists of the elements at m[0][0], m[1][1], m[2][2]. */
	CONST_WIN32 CVec3D diagonal() const;

	/// Sets the values of the given row.
	/** \param row The index of the row to set, in the range [0-2].
		\param data A pointer to an array of 4 floats that contain the new x, y, z and w values for the row. */
	void setRow(int row, const float *data);
	void setRow(int row, const CVec3D &rowVector, float m_r3);
	void setRow(int row, const CVec4D &rowVector);
	void setRow(int row, float m_r0, float m_r1, float m_r2, float m_r3);

	/// Sets the values of the given column.
	/** \param column The index of the column to set, in the range [0-3].
		\param data A pointer to an array of 3 floats that contain the new x, y and z values for the column. */
	void setCol(int column, const float *data);
	void setCol(int column, const CVec3D &columnVector);
	void setCol(int column, float m_0c, float m_1c, float m_2c);

	/// Sets all values of this matrix.
	void set(float _00, float _01, float _02, float _03,
			 float _10, float _11, float _12, float _13,
			 float _20, float _21, float _22, float _23);

	void set3x3Part(const CMat3D &rotation);

	/**
	\brief Sets the top-left 3x3 area of the matrix to the rotation matrix about the X-axis. Elements
	outside the top-left 3x3 area are ignored. This matrix rotates counterclockwise if multiplied
	in the order M*v, and clockwise if rotated in the order v*M.
	\param m The matrix to store the result.
	\param angle The rotation angle in radians.
	*/
	void set3x3PartRotateX(float angle);
	/**
	\brief Sets the top-left 3x3 area of the matrix to the rotation matrix about the Y-axis. Elements
	outside the top-left 3x3 area are ignored. This matrix rotates counterclockwise if multiplied
	in the order M*v, and clockwise if rotated in the order v*M.
	\param m The matrix to store the result.
	\param angle The rotation angle in radians.
	*/
	void set3x3PartRotateY(float angle);
	/**
	\brief Sets the top-left 3x3 area of the matrix to the rotation matrix about the Z-axis. Elements
	outside the top-left 3x3 area are ignored. This matrix rotates counterclockwise if multiplied
	in the order M*v, and clockwise if rotated in the order v*M.
	\param m The matrix to store the result.
	\param angle The rotation angle in radians.
    */
	void set3x3PartRotateZ(float angle);

	void setMatrixRotatePart(const CQuaternion &q);
	/// Sets the 3-by-3 part of this matrix to perform rotation about the positive X axis which passes through
	/// the origin. Leaves all other entries of this matrix untouched. [similarOverload: SetRotatePart] [hideIndex]
	void setRotatePartX(float angleRadians);
	/// Sets the 3-by-3 part of this matrix to perform rotation about the positive Y axis. Leaves all other
	/// entries untouched. [similarOverload: SetRotatePart] [hideIndex]
	void setRotatePartY(float angleRadians);
	/// Sets the 3-by-3 part of this matrix to perform rotation about the positive Z axis. Leaves all other
	/// entries untouched. [similarOverload: SetRotatePart] [hideIndex]
	void setRotatePartZ(float angleRadians);

	/// Sets the 3-by-3 part of this matrix to perform rotation about the given axis and angle. Leaves all other
	/// entries of this matrix untouched. [indexTitle: SetRotatePart/X/Y/Z]
	void setRotatePart(const CVec3D &axisDirection, float angleRadians);
	/// Sets the 3-by-3 part of this matrix to perform the rotation expressed by the given quaternion.
	/// Leaves all other entries of this matrix untouched.
	void setRotatePart(const CQuaternion &orientation);
	/// Sets the 3-by-3 part of this matrix.
	/// @note This is a convenience function which calls Set3x3Part.
	/// @note This function erases the previous top-left 3x3 part of this matrix (any previous rotation, scaling and shearing, etc.). Translation is unaffected.
	void setRotatePart(const CMat3D &rotation) { set3x3Part(rotation); }

	void			setRotation( const MATH::CMat3D &m );





	void			setTranslation( const MATH::CVec3D &t );
	const float &   operator[]( int index ) const;

	CVec3D	operator*( const CVec3D &v ) const;							// only rotate

	CVec3D	operator*( const CVec4D &v ) const;							// rotate and translate

	/// Transforms the given vector by this matrix (in the order M * v).
	/// The fourth element of the vector is preserved.
	CVec4D mul(const CVec4D &rhs) const;

	/// Treats the CMat3D as a 4-by-4 matrix with the last row and column as identity, and multiplies the two matrices.
	CMatJoint3x4 operator *(const CMat3D &rhs) const;

	/// Treats the CMatJoint3x4 as a 4-by-4 matrix with the last row as identity, and multiplies the two matrices.
	CMatJoint3x4 operator *(const CMatJoint3x4 &rhs) const;

	/// Converts the quaternion to a CMatJoint3x4 and multiplies the two matrices together.
	CMatJoint3x4 operator *(const CQuaternion &rhs) const;

	/// Transforms the given vector by this matrix (in the order M * v).
	/// The fourth element of the vector is preserved.
	//CVec4D operator *(const CVec4D &rhs) const;

	CMatJoint3x4 operator *(float scalar) const;

	CMatJoint3x4 &	operator*=( const CMatJoint3x4 &a );							// transform
	CMatJoint3x4 &	operator/=( const CMatJoint3x4 &a );							// untransform

	CMatJoint3x4 &operator *=(float scalar);
	CMatJoint3x4 &operator /=(float scalar);
	/// exact compare, no epsilon
	bool			compare( const CMatJoint3x4 &a ) const;
	/// compare with epsilon
	bool			compare( const CMatJoint3x4 &a, const float epsilon ) const;
	/// exact compare, no epsilon
	bool			operator==(	const CMatJoint3x4 &a ) const;
	/// exact compare, no epsilon
	bool			operator!=(	const CMatJoint3x4 &a ) const;

	/// Transforms the given point vector by this matrix M , i.e. returns M * (x, y, z, 1).
	CVec3D transformPos(const CVec3D &pointVector) const;
	CVec3D transformPos(float x, float y, float z) const;

	/// Transforms the given direction vector by this matrix M , i.e. returns M * (x, y, z, 0).
	CVec3D transformDir(const CVec3D &directionVector) const;
	CVec3D transformDir(float x, float y, float z) const;

	/// Transforms the given 4-vector by this matrix M, i.e. returns M * (x, y, z, w).
	CVec4D transform(const CVec4D &vector) const;
		/**
		\brief Computes the determinant of this matrix.
		If the determinant is nonzero, this matrix is invertible.
		If the determinant is negative, this matrix performs reflection about some axis.
		"If the determinant is positive, the basis is said to be "positively" oriented (or right-handed).
		If the determinant is negative, the basis is said to be "negatively" oriented (or left-handed)."
		\see http://4msdn.microsoft.com/en-us/library/bb204853(VS.85).aspx :
		*/
	float determinant() const;

	/// Returns the upper-left 3-by-3 part.
	CONST_WIN32 CMat3D CMat3DPart() const;

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

	MATH::CMat3D			toMat3() const;
	MATH::CVec3D			toVec3() const;
	CJointQuaternion		ToJointQuat() const;
	const float *	toFloatPtr() const;
	float *			toFloatPtr();
	sf_m128 *		toM128Ptr();
	const sf_m128 *		toM128Ptr()const;

	CMatJoint3x4 mul(const CMat3D &rhs) const { return *this * rhs; };
	CMatJoint3x4 mul(const CMatJoint3x4 &rhs) const { return *this * rhs; };
//	CMat4D mul(const CMat4D &rhs) const { return *this * rhs; };
	CMatJoint3x4 mul(const CQuaternion &rhs) const { return *this * rhs; };
	CVec3D MulPos(const CVec3D &pointVector) const { return this->transformPos(pointVector); };
	CVec3D MulDir(const CVec3D &directionVector) const { return this->transformDir(directionVector); };
//	CVec4D mul(const CVec4D &vector) const { return *this * vector; };

	/**
	\brief Tests if this matrix does not contain any NaNs or infs.
	\return Returns true if the entries of this CMatJoint3x4 are all finite, and do not contain NaN or infs.
	**/
	bool isFinite() const;

	/**
	\brief Tests if this is the identity matrix.
	\return Returns true if this matrix is the identity matrix, up to the given epsilon.
	**/
	bool isIdentity(float epsilon =CMath::EPSILON_SuperLow) const;

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
	 This function treats this matrix as a square 4x4 matrix with the last row having the form [0 0 0 1].
	\return Returns true if this matrix can be inverted, up to the given epsilon.
	**/
	bool isInvertible(float epsilon =CMath::EPSILON_SuperLow) const;

	/**
	\brief Tests if this matrix is symmetric (M == M^T).
	 This function treats this matrix as a square 4x4 matrix with the last row having the form [0 0 0 1].
	The test compares the elements for equality, up to the given epsilon. A matrix is symmetric if it is its own transpose.
	**/
	bool isSymmetric(float epsilon =CMath::EPSILON_SuperLow) const;

	/**
	\brief Tests if this matrix is skew-symmetric (M == -M^T).
	 This function treats this matrix as a square 4x4 matrix with the last row having the form [0 0 0 1].
		The test compares the floating point elements of this matrix up to the given epsilon. A matrix M is skew-symmetric
		the identity M=-M^T holds.
	*/
	bool isSkewSymmetric(float epsilon =CMath::EPSILON_SuperLow) const;

	/**
	\brief Returns true if this matrix does not perform any scaling.
	 A matrix does not do any scaling if the column vectors of this
		matrix are normalized in length, compared to the given epsilon. Note that this matrix may still perform
		reflection, i.e. it has a -1 scale along some axis.
		\note This function only examines the upper 3-by-3 part of this matrix.
		\note This function assumes that this matrix does not contain projection (the fourth row of this matrix is [0 0 0 1]). */
	bool hasUnitaryScale(float epsilonSq = 1e-6f) const;

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

	/**
	\briefReturns true if the row vectors of 3x3 top-left submatrix are all perpendicular to each other.
	**/
	bool isRowOrthogonal(float epsilon =CMath::EPSILON_SuperLow) const;

	/**
	/brief Returns true if the column vectors of 3x3 top-left submatrix are all perpendicular to each other.
	*/
	bool isColOrthogonal(float epsilon =CMath::EPSILON_SuperLow) const;
	bool isColOrthogonal3(float epsilon =CMath::EPSILON_SuperLow) const { return isColOrthogonal(epsilon); }

	/**
	/brief Returns true if the column and row vectors of the 3x3 top-left submatrix form an orthonormal set.
	 \note In math terms, there does not exist such a thing as 'orthonormal matrix'. In math terms, a matrix
	 is orthogonal iff its column and row vectors are orthogonal *unit* vectors.
	 In the terms of this library however, a matrix is orthogonal iff its column and row vectors are orthogonal (no need to be unitary),
	 and a matrix is orthonormal if the column and row vectors are orthonormal.
	 **/
	bool isOrthonormal(float epsilon =CMath::EPSILON_SuperLow) const;


	/**
	\brief Sets the translation part of this matrix.
	This function sets the translation part of this matrix. These are the three first elements of the fourth column.
	All other entries are left untouched.
	*/
	void setTranslatePart(float tx, float ty, float tz) { setCol(3, tx, ty, tz); }
	void setTranslatePart(const CVec3D &offset) { setCol(3, offset); }

	/**
	\brief Inverts a column-orthogonal matrix.
	If a matrix is of form M=T*R*S, where T is an affine translation matrix,
	R is a rotation matrix and S is a diagonal matrix with non-zero but potentially non-uniform scaling
	factors (possibly mirroring), then the matrix M is column-orthogonal and this function can be used to compute the inverse.
	Calling this function is faster than the calling the generic matrix Inverse() function.
	Returns true on success. On failure, the matrix is not modified. This function fails if any of the
	elements of this vector are not finite, or if the matrix contains a zero scaling factor on X, Y or Z.
	\note The returned matrix will be row-orthogonal, but not column-orthogonal in general.
	The returned matrix will be column-orthogonal iff the original matrix M was row-orthogonal as well.
	(in which case S had uniform scale, inverseOrthogonalUniformScale() could have been used instead)
	**/
	bool inverseColOrthogonal();

	/**
	\brief Inverts a matrix that is a concatenation of only translate, rotate and uniform scale operations.
	If a matrix is of form M=T*R*S, where T is an affine translation matrix,
	R is a rotation matrix and S is a diagonal matrix with non-zero and uniform scaling factors (possibly mirroring),
	then the matrix M is both column- and row-orthogonal and this function can be used to compute the inverse.
	This function is faster than calling inverseColOrthogonal() or the generic Inverse().
	Returns true on success. On failure, the matrix is not modified. This function fails if any of the
	elements of this vector are not finite, or if the matrix contains a zero scaling factor on X, Y or Z.
	This function may not be called if this matrix contains any shearing or nonuniform scaling.
	**/
	bool inverseOrthogonalUniformScale();

	/**
	\brief Inverts a matrix that is a concatenation of only translate and rotate operations.
	If a matrix is of form M=T*R*S, where T is an affine translation matrix, R is a rotation
	matrix and S is either identity or a mirroring matrix, then the matrix M is orthonormal and this function can be used to compute the inverse.
	This function is faster than calling inverseOrthogonalUniformScale(), inverseColOrthogonal() or the
	generic Inverse().
	This function may not be called if this matrix contains any scaling or shearing, but it may contain mirroring.
	**/
	void inverseOrthonormal();

	/// Performs a batch transform of the given array.
	void batchTransformPos(CVec3D *pointArray, int numPoints) const;

	/// Performs a batch transform of the given array.
	void batchTransformPos(CVec3D *pointArray, int numPoints, int stride) const;

	/// Performs a batch transform of the given array.
	void batchTransformDir(CVec3D *dirArray, int numVectors) const;

	/// Performs a batch transform of the given array.
	void batchTransformDir(CVec3D *dirArray, int numVectors, int stride) const;

	/// Performs a batch transform of the given array.
	void batchTransform(CVec4D *vectorArray, int numVectors) const;

	/// Performs a batch transform of the given array.
	void batchTransform(CVec4D *vectorArray, int numVectors, int stride) const;
		/// Creates a new transformation matrix that scales by the given factors.

	/**
	\brief Returns the translation part.
	The translation part is stored in the fourth column of this matrix.
	This is equivalent to decomposing this matrix in the form M = T * M', i.e. this translation is applied last,
	after applying rotation and scale. If this matrix represents a local->world space transformation for an object,
	then this gives the world space position of the object.
	\note This function assumes that this matrix does not contain projection (the fourth row of this matrix is [0 0 0 1]).
	*/
	CONST_WIN32 CVec3D translatePart() const;

	/// Returns the upper-left 3x3 part of this matrix. This part stores the rotation of this transform.
	CONST_WIN32 CMat3D rotatePart() const;


	/**
	\brief Transposes the top-left 3x3 part of this matrix in-place. The fourth column (translation part) will
	remain intact.
	This operation swaps all elements with respect to the diagonal.
	**/
	void transpose3();

	/// Returns a copy of this matrix which has the top-left 3x3 part transposed.
	CMatJoint3x4 transposed3() const;

	/**
	\brief Computes the inverse transpose of this matrix in-place.
	 Use the inverse transpose to transform covariant vectors (normal vectors).
	\note This function resets the translation part of this matrix to zero.
	*/
	bool inverseTranspose();
	/**
	\brief Inverts this matrix using the generic Gauss's method.
	\return Returns true on success, false otherwise.
	**/
	bool inverse(float epsilon = 1e-3f);

	/**
	\brief Returns the inverse transpose of this matrix.
	 Use that matrix to transform covariant vectors (normal vectors).
	\note This function resets the translation part of this matrix to zero.
	*/
	CMatJoint3x4 inverseTransposed() const;

	/// Returns the sum of the diagonal elements of this matrix.
	float trace() const;


	/// Creates a new transformation matrix that scales by the given factors.
	/// This matrix scales with respect to origin.
	static CScaleOp scale(float sx, float sy, float sz);
	static CScaleOp scale(const CVec3D &scale);
	/// Creates a new transformation matrix that translates by the given offset.
	static CTranslateOp translate(float tx, float ty, float tz);
	static CTranslateOp translate(const CVec3D &offset);



	/** Creates a new CMatJoint3x4 that scales with respect to the given center point.
	\param scale The amount of scale to apply to the x, y and z directions.
	\param scaleCenter The coordinate system center point for the scaling. If omitted, the origin (0,0,0) will
	be used as the origin for the scale operation.
	**/
	static CMatJoint3x4 scale(const CVec3D &scale, const CVec3D &scaleCenter);





private:
	/**
	\ note
	[lin][col]
	[0][0] = mat[0 * 4 + 0]
	[0][1] = mat[0 * 4 + 1]
	[0][2] = mat[0 * 4 + 2]
	[0][3] = mat[0 * 4 + 3]
	[1][0] = mat[1 * 4 + 0]
	[1][1] = mat[1 * 4 + 1]
	[1][2] = mat[1 * 4 + 2]
	[1][3] = mat[1 * 4 + 3]
	[2][0] = mat[2 * 4 + 0]
	[2][1] = mat[2 * 4 + 1]
	[2][2] = mat[2 * 4 + 2]
	[2][3] = mat[2 * 4 + 3]

	**/
	union {
	float			mat[3*4];
	sf_m128	        q[3];
	};
};

SMF_INLINE_FORCED CMatJoint3x4::CMatJoint3x4(float _00, float _01, float _02, float _03,
		 float _10, float _11, float _12, float _13,
		 float _20, float _21, float _22, float _23)
{
	set(_00, _01, _02, _03,
		_10, _11, _12, _13,
		_20, _21, _22, _23);
}

SMF_INLINE_FORCED const float &CMatJoint3x4::operator[]( int index ) const {
	SMF_ASSERT( ( index < 0 ) && ( index >= 4 ) );
	return mat[ index ];
}

SMF_INLINE_FORCED void CMatJoint3x4::setRotation( const MATH::CMat3D &m ) {
	// NOTE: MATH::CMat3D is transposed because it is column-major
	mat[0 * 4 + 0] = m[0][0];
	mat[0 * 4 + 1] = m[1][0];
	mat[0 * 4 + 2] = m[2][0];
	mat[1 * 4 + 0] = m[0][1];
	mat[1 * 4 + 1] = m[1][1];
	mat[1 * 4 + 2] = m[2][1];
	mat[2 * 4 + 0] = m[0][2];
	mat[2 * 4 + 1] = m[1][2];
	mat[2 * 4 + 2] = m[2][2];
}

SMF_INLINE_FORCED void CMatJoint3x4::setTranslation( const CVec3D &t ) {
	mat[0 * 4 + 3] = t[0];
	mat[1 * 4 + 3] = t[1];
	mat[2 * 4 + 3] = t[2];
}

SMF_INLINE_FORCED MATH::CVec3D CMatJoint3x4::operator*( const CVec3D &v ) const {
	return MATH::CVec3D(	mat[0 * 4 + 0] * v[0] + mat[0 * 4 + 1] * v[1] + mat[0 * 4 + 2] * v[2],
					mat[1 * 4 + 0] * v[0] + mat[1 * 4 + 1] * v[1] + mat[1 * 4 + 2] * v[2],
					mat[2 * 4 + 0] * v[0] + mat[2 * 4 + 1] * v[1] + mat[2 * 4 + 2] * v[2] );
}

SMF_INLINE_FORCED CVec3D CMatJoint3x4::operator*( const CVec4D &v ) const {
	return MATH::CVec3D(	mat[0 * 4 + 0] * v[0] + mat[0 * 4 + 1] * v[1] + mat[0 * 4 + 2] * v[2] + mat[0 * 4 + 3] * v[3],
					mat[1 * 4 + 0] * v[0] + mat[1 * 4 + 1] * v[1] + mat[1 * 4 + 2] * v[2] + mat[1 * 4 + 3] * v[3],
					mat[2 * 4 + 0] * v[0] + mat[2 * 4 + 1] * v[1] + mat[2 * 4 + 2] * v[2] + mat[2 * 4 + 3] * v[3] );
}
SMF_INLINE_FORCED CMatJoint3x4 CMatJoint3x4::operator *(const CMat3D &rhs) const
{

	CMatJoint3x4 r;
	CVec3D lin0_3(mat[0 * 4 + 0],mat[0 * 4 + 1],mat[0 * 4 + 2]);
	CVec3D lin1_3(mat[1 * 4 + 0],mat[1 * 4 + 1],mat[1 * 4 + 2]);
	CVec3D lin2_3(mat[2 * 4 + 0],mat[2 * 4 + 1],mat[2 * 4 + 2]);

	CVec3D rhsCol0_3=rhs.Col(0);
	CVec3D rhsCol1_3=rhs.Col(1);
	CVec3D rhsCol2_3=rhs.Col(2);

	//linha 0
	r.mat[0 * 4 + 0] = lin0_3*rhsCol0_3;
	r.mat[0 * 4 + 1] = lin0_3*rhsCol1_3;
	r.mat[0 * 4 + 2] = lin0_3*rhsCol2_3;
	r.mat[0 * 4 + 3] = mat[0 * 4 + 3];
	//linha 1
	r.mat[1 * 4 + 0] = lin1_3*rhsCol0_3;
	r.mat[1 * 4 + 1] = lin1_3*rhsCol1_3;
	r.mat[1 * 4 + 2] = lin1_3*rhsCol2_3;
	r.mat[1 * 4 + 3] = mat[1 * 4 + 3];
	//linha 2
	r.mat[2 * 4 + 0] = lin2_3*rhsCol0_3;
	r.mat[2 * 4 + 1] = lin2_3*rhsCol1_3;
	r.mat[2 * 4 + 2] = lin2_3*rhsCol2_3;
	r.mat[2 * 4 + 3] = mat[2 * 4 + 3];

	return r;
}




SMF_INLINE_FORCED CVec4D &CMatJoint3x4::Row(int row)
{
	SMF_ASSERT(row >= 0);
	SMF_ASSERT(row < 3);

#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (row < 0 || row >= 3)
		row = 0; // Benign failure, just give the first row.
#endif
	return reinterpret_cast<CVec4D &>(mat[row*4]);
}

SMF_INLINE_FORCED const CVec4D &CMatJoint3x4::Row(int row) const
{
	SMF_ASSERT(row >= 0);
	SMF_ASSERT(row < 3);

#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (row < 0 || row >= 3)
		row = 0; // Benign failure, just give the first row.
#endif
	return reinterpret_cast<const CVec4D &>(mat[row*4]);
}

SMF_INLINE_FORCED CVec3D &CMatJoint3x4::Row3(int row)
{
	SMF_ASSERT(row >= 0);
	SMF_ASSERT(row < 3);

#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (row < 0 || row >= 3)
		row = 0; // Benign failure, just give the first row.
#endif
	return reinterpret_cast<CVec3D &>(mat[row*4]);
}

SMF_INLINE_FORCED const CVec3D &CMatJoint3x4::Row3(int row) const
{
	SMF_ASSERT(row >= 0);
	SMF_ASSERT(row < 3);

#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (row < 0 || row >= 3)
		row = 0; // Benign failure, just give the first row.
#endif
	return reinterpret_cast<const CVec3D &>(mat[row*4]);
}

SMF_INLINE_FORCED CONST_WIN32 CVec3D CMatJoint3x4::Col(int col) const
{
	SMF_ASSERT(col >= 0);
	SMF_ASSERT(col < 4);
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (col < 0 || col >= 4)
		return CVec3D::nan;
#endif
	return CVec3D(mat[0 * 4 + col], mat[1 * 4 + col], mat[2 * 4 + col]);
}

SMF_INLINE_FORCED CONST_WIN32 CVec4D CMatJoint3x4::Col4(int col) const
{
	SMF_ASSERT(col >= 0);
	SMF_ASSERT(col < 4);
#ifndef MATH_ENABLE_INSECURE_OPTIMIZATIONS
	if (col < 0 || col >= 4)
		return CVec4D::nan;
#endif
	if (col == 2){
	return CVec4D(mat[0 * 4 + col], mat[1 * 4 + col], mat[2 * 4 + col],1.0f);
	}else
	return CVec4D(mat[0 * 4 + col], mat[1 * 4 + col], mat[2 * 4 + col],0);
}

SMF_INLINE_FORCED CONST_WIN32 CVec3D CMatJoint3x4::diagonal() const
{
	return CVec3D(mat[0 * 4 + 0], mat[1 * 4 + 1], mat[2 * 4 + 2]);
}


SMF_INLINE_FORCED CMatJoint3x4 CMatJoint3x4::operator *(const CMatJoint3x4 &rhs) const
{
	CMatJoint3x4 r;
if(CMath::MATH_AUTOMATIC_SSE){
	CSIMD::getProcessor()->mat3x4_mul_sse(r.Row(0).toM128Ptr(), Row(0).toM128Ptr(), rhs.Row(0).toM128Ptr());
}else{
	CSIMD::getProcessor()->mat3x4_mul_sse(r.Row(0).toM128Ptr(), Row(0).toM128Ptr(), rhs.Row(0).toM128Ptr());
}

	return r;
}

SMF_INLINE_FORCED CMatJoint3x4 CMatJoint3x4::operator *(const CQuaternion &rhs) const
{
	CMat3D rot=rhs.toMat3();
	return *this * rot;
}
/*
SMF_INLINE_FORCED CVec4D CMatJoint3x4::operator *(const CVec4D &rhs) const
{
	return transform(rhs);
}
*/

SMF_INLINE_FORCED CMatJoint3x4 &CMatJoint3x4::operator*=( const CMatJoint3x4 &a ) {
	float dst[3];

	dst[0] = mat[0 * 4 + 0] * a.mat[0 * 4 + 0] + mat[1 * 4 + 0] * a.mat[0 * 4 + 1] + mat[2 * 4 + 0] * a.mat[0 * 4 + 2];
	dst[1] = mat[0 * 4 + 0] * a.mat[1 * 4 + 0] + mat[1 * 4 + 0] * a.mat[1 * 4 + 1] + mat[2 * 4 + 0] * a.mat[1 * 4 + 2];
	dst[2] = mat[0 * 4 + 0] * a.mat[2 * 4 + 0] + mat[1 * 4 + 0] * a.mat[2 * 4 + 1] + mat[2 * 4 + 0] * a.mat[2 * 4 + 2];
	mat[0 * 4 + 0] = dst[0];
	mat[1 * 4 + 0] = dst[1];
	mat[2 * 4 + 0] = dst[2];

	dst[0] = mat[0 * 4 + 1] * a.mat[0 * 4 + 0] + mat[1 * 4 + 1] * a.mat[0 * 4 + 1] + mat[2 * 4 + 1] * a.mat[0 * 4 + 2];
	dst[1] = mat[0 * 4 + 1] * a.mat[1 * 4 + 0] + mat[1 * 4 + 1] * a.mat[1 * 4 + 1] + mat[2 * 4 + 1] * a.mat[1 * 4 + 2];
	dst[2] = mat[0 * 4 + 1] * a.mat[2 * 4 + 0] + mat[1 * 4 + 1] * a.mat[2 * 4 + 1] + mat[2 * 4 + 1] * a.mat[2 * 4 + 2];
	mat[0 * 4 + 1] = dst[0];
	mat[1 * 4 + 1] = dst[1];
	mat[2 * 4 + 1] = dst[2];

	dst[0] = mat[0 * 4 + 2] * a.mat[0 * 4 + 0] + mat[1 * 4 + 2] * a.mat[0 * 4 + 1] + mat[2 * 4 + 2] * a.mat[0 * 4 + 2];
	dst[1] = mat[0 * 4 + 2] * a.mat[1 * 4 + 0] + mat[1 * 4 + 2] * a.mat[1 * 4 + 1] + mat[2 * 4 + 2] * a.mat[1 * 4 + 2];
	dst[2] = mat[0 * 4 + 2] * a.mat[2 * 4 + 0] + mat[1 * 4 + 2] * a.mat[2 * 4 + 1] + mat[2 * 4 + 2] * a.mat[2 * 4 + 2];
	mat[0 * 4 + 2] = dst[0];
	mat[1 * 4 + 2] = dst[1];
	mat[2 * 4 + 2] = dst[2];

	dst[0] = mat[0 * 4 + 3] * a.mat[0 * 4 + 0] + mat[1 * 4 + 3] * a.mat[0 * 4 + 1] + mat[2 * 4 + 3] * a.mat[0 * 4 + 2];
	dst[1] = mat[0 * 4 + 3] * a.mat[1 * 4 + 0] + mat[1 * 4 + 3] * a.mat[1 * 4 + 1] + mat[2 * 4 + 3] * a.mat[1 * 4 + 2];
	dst[2] = mat[0 * 4 + 3] * a.mat[2 * 4 + 0] + mat[1 * 4 + 3] * a.mat[2 * 4 + 1] + mat[2 * 4 + 3] * a.mat[2 * 4 + 2];
	mat[0 * 4 + 3] = dst[0];
	mat[1 * 4 + 3] = dst[1];
	mat[2 * 4 + 3] = dst[2];

	mat[0 * 4 + 3] += a.mat[0 * 4 + 3];
	mat[1 * 4 + 3] += a.mat[1 * 4 + 3];
	mat[2 * 4 + 3] += a.mat[2 * 4 + 3];

	return *this;
}

SMF_INLINE_FORCED CMatJoint3x4 &CMatJoint3x4::operator/=( const CMatJoint3x4 &a ) {
	float dst[3];

	mat[0 * 4 + 3] -= a.mat[0 * 4 + 3];
	mat[1 * 4 + 3] -= a.mat[1 * 4 + 3];
	mat[2 * 4 + 3] -= a.mat[2 * 4 + 3];

	dst[0] = mat[0 * 4 + 0] * a.mat[0 * 4 + 0] + mat[1 * 4 + 0] * a.mat[1 * 4 + 0] + mat[2 * 4 + 0] * a.mat[2 * 4 + 0];
	dst[1] = mat[0 * 4 + 0] * a.mat[0 * 4 + 1] + mat[1 * 4 + 0] * a.mat[1 * 4 + 1] + mat[2 * 4 + 0] * a.mat[2 * 4 + 1];
	dst[2] = mat[0 * 4 + 0] * a.mat[0 * 4 + 2] + mat[1 * 4 + 0] * a.mat[1 * 4 + 2] + mat[2 * 4 + 0] * a.mat[2 * 4 + 2];
	mat[0 * 4 + 0] = dst[0];
	mat[1 * 4 + 0] = dst[1];
	mat[2 * 4 + 0] = dst[2];

	dst[0] = mat[0 * 4 + 1] * a.mat[0 * 4 + 0] + mat[1 * 4 + 1] * a.mat[1 * 4 + 0] + mat[2 * 4 + 1] * a.mat[2 * 4 + 0];
	dst[1] = mat[0 * 4 + 1] * a.mat[0 * 4 + 1] + mat[1 * 4 + 1] * a.mat[1 * 4 + 1] + mat[2 * 4 + 1] * a.mat[2 * 4 + 1];
	dst[2] = mat[0 * 4 + 1] * a.mat[0 * 4 + 2] + mat[1 * 4 + 1] * a.mat[1 * 4 + 2] + mat[2 * 4 + 1] * a.mat[2 * 4 + 2];
	mat[0 * 4 + 1] = dst[0];
	mat[1 * 4 + 1] = dst[1];
	mat[2 * 4 + 1] = dst[2];

	dst[0] = mat[0 * 4 + 2] * a.mat[0 * 4 + 0] + mat[1 * 4 + 2] * a.mat[1 * 4 + 0] + mat[2 * 4 + 2] * a.mat[2 * 4 + 0];
	dst[1] = mat[0 * 4 + 2] * a.mat[0 * 4 + 1] + mat[1 * 4 + 2] * a.mat[1 * 4 + 1] + mat[2 * 4 + 2] * a.mat[2 * 4 + 1];
	dst[2] = mat[0 * 4 + 2] * a.mat[0 * 4 + 2] + mat[1 * 4 + 2] * a.mat[1 * 4 + 2] + mat[2 * 4 + 2] * a.mat[2 * 4 + 2];
	mat[0 * 4 + 2] = dst[0];
	mat[1 * 4 + 2] = dst[1];
	mat[2 * 4 + 2] = dst[2];

	dst[0] = mat[0 * 4 + 3] * a.mat[0 * 4 + 0] + mat[1 * 4 + 3] * a.mat[1 * 4 + 0] + mat[2 * 4 + 3] * a.mat[2 * 4 + 0];
	dst[1] = mat[0 * 4 + 3] * a.mat[0 * 4 + 1] + mat[1 * 4 + 3] * a.mat[1 * 4 + 1] + mat[2 * 4 + 3] * a.mat[2 * 4 + 1];
	dst[2] = mat[0 * 4 + 3] * a.mat[0 * 4 + 2] + mat[1 * 4 + 3] * a.mat[1 * 4 + 2] + mat[2 * 4 + 3] * a.mat[2 * 4 + 2];
	mat[0 * 4 + 3] = dst[0];
	mat[1 * 4 + 3] = dst[1];
	mat[2 * 4 + 3] = dst[2];

	return *this;
}

SMF_INLINE_FORCED bool CMatJoint3x4::compare( const CMatJoint3x4 &a ) const {
	int i;

	for ( i = 0; i < 12; i++ ) {
		if ( mat[i] != a.mat[i] ) {
			return false;
		}
	}
	return true;
}

SMF_INLINE_FORCED bool CMatJoint3x4::compare( const CMatJoint3x4 &a, const float epsilon ) const {
	int i;

	for ( i = 0; i < 12; i++ ) {
		if ( MATH::CMath::fabs( mat[i] - a.mat[i] ) > epsilon ) {
			return false;
		}
	}
	return true;
}

SMF_INLINE_FORCED bool CMatJoint3x4::operator==( const CMatJoint3x4 &a ) const {
	return compare( a );
}

SMF_INLINE_FORCED bool CMatJoint3x4::operator!=( const CMatJoint3x4 &a ) const {
	return !compare( a );
}

SMF_INLINE_FORCED MATH::CMat3D CMatJoint3x4::toMat3() const {
	return MATH::CMat3D(	mat[0 * 4 + 0], mat[1 * 4 + 0], mat[2 * 4 + 0],
					mat[0 * 4 + 1], mat[1 * 4 + 1], mat[2 * 4 + 1],
					mat[0 * 4 + 2], mat[1 * 4 + 2], mat[2 * 4 + 2] );
}

SMF_INLINE_FORCED MATH::CVec3D CMatJoint3x4::toVec3() const {
	return MATH::CVec3D( mat[0 * 4 + 3], mat[1 * 4 + 3], mat[2 * 4 + 3] );
}

SMF_INLINE_FORCED const float *CMatJoint3x4::toFloatPtr() const {
	return mat;
}

SMF_INLINE_FORCED float *CMatJoint3x4::toFloatPtr() {
	return mat;
}
SMF_INLINE_FORCED const sf_m128 *CMatJoint3x4::toM128Ptr() const {
	return reinterpret_cast<const sf_m128*>(&mat);
}

SMF_INLINE_FORCED sf_m128 *CMatJoint3x4::toM128Ptr() {
	return reinterpret_cast<sf_m128*>(&mat);
}


SMF_INLINE_FORCED CONST_WIN32 CMat3D CMatJoint3x4::CMat3DPart() const
{
	return CMat3D(mat[0*4+0], mat[0*4+1], mat[0*4+2],
					mat[1*4+0], mat[1*4+1], mat[1*4+2],
					mat[2*4+0], mat[2*4+1], mat[2*4+2]);
}

SMF_INLINE_FORCED CVec3D CMatJoint3x4::extractScale() const
{
	return CVec3D(Col(0).getLenght(), Col(1).getLenght(), Col(2).getLenght());
}

} //end MATH
} //end SMF
#endif /* !__JOINTTRANSFORM_H__ */
