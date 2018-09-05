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
#ifndef __SMF_TRANLATE_OPER__
#define __SMF_TRANLATE_OPER__

#include "../SMF_Config.h"


namespace SMF{
namespace MATH{
class CVec3D;
class CMat3D;
class CVec4D;
class CMat4D;
class CMatJoint3x4;

// A structure that represents the translate operation for 3D objects.
// This structure is used to optimize special cases of 3D transformation concatenations. The use of this
//	class occurs transparently to the user. You do not need to instantiate new CTranslateOp objects in your code. */
class SMF_API CTranslateOp
{
public:
	/// The x offset of translation.
	float x;
	/// The y offset of translation.
	float y;
	/// The z offset of translation.
	float z;

	/// Constructs an uninitialized CTranslateOp.
	CTranslateOp() {}

	/// Constructs a CTranslateOp that translates the given amount.
	explicit CTranslateOp(const CVec3D &offset);
	CTranslateOp(float x, float y, float z);

	/// Returns the translation offset (x, y, z).
	CVec3D Offset() const;

	/// Converts this CTranslateOp object to a matrix.
	CMatJoint3x4 toMat3x4() const;
	/// Converts this CTranslateOp object to a matrix.
	CMat4D toMat4D() const;

	/// Converts this CTranslateOp object to a matrix.
	operator CMatJoint3x4() const;
	/// Converts this CTranslateOp object to a matrix.
	operator CMat4D() const;
};

CMatJoint3x4 operator *(const CTranslateOp &lhs, const CMatJoint3x4 &rhs);
CMatJoint3x4 operator *(const CMatJoint3x4 &lhs, const CTranslateOp &rhs);
CMat4D operator *(const CTranslateOp &lhs, const CMat4D &rhs);
CMat4D operator *(const CMat4D &lhs, const CTranslateOp &rhs);

// A structure that represents the scale operation for 3D objects.
// This structure is used to optimize special cases of 3D transformation concatenations. The use of this
//	class occurs transparently to the user. You do not need to instantiate new CScaleOp objects in your code. */
class SMF_API CScaleOp
{
public:
	/// The scale factor along the x axis.
	float x;
	/// The scale factor along the y axis.
	float y;
	/// The scale factor along the z axis.
	float z;

	/// Constructs an uninitialized CScaleOp.
	CScaleOp() {}

	/// Constructs a CScaleOp with the given scale factors.
	explicit CScaleOp(const CVec3D &scale);
	CScaleOp(float sx, float sy, float sz);

	/// Returns the scale factors (x, y, z).
	CVec3D Offset() const;

	/// Converts this CScaleOp to a matrix.
	operator CMat3D() const;
	/// Converts this CScaleOp to a matrix.
	operator CMatJoint3x4() const;
	/// Converts this CScaleOp to a matrix.
	operator CMat4D() const;

	/// Converts this CScaleOp to a matrix.
	CMat3D ToFloat3x3() const;
	/// Converts this CScaleOp to a matrix.
	CMatJoint3x4 toMat3x4() const;
	/// Converts this CScaleOp to a matrix.
	CMat4D toMat4D() const;
};

CMat3D operator *(const CScaleOp &lhs, const CMat3D &rhs);
CMat3D operator *(const CMat3D &lhs, const CScaleOp &rhs);
CMatJoint3x4 operator *(const CScaleOp &lhs, const CMatJoint3x4 &rhs);
CMatJoint3x4 operator *(const CMatJoint3x4 &lhs, const CScaleOp &rhs);
CMat4D operator *(const CScaleOp &lhs, const CMat4D &rhs);
CMat4D operator *(const CMat4D &lhs, const CScaleOp &rhs);

CMatJoint3x4 operator *(const CScaleOp &lhs, const CTranslateOp &rhs);
CMatJoint3x4 operator *(const CTranslateOp &lhs, const CScaleOp &rhs);

#ifdef MATH_QT_INTEROP
Q_DECLARE_METATYPE(CTranslateOp)
Q_DECLARE_METATYPE(CTranslateOp*)
Q_DECLARE_METATYPE(CScaleOp)
Q_DECLARE_METATYPE(CScaleOp*)
#endif

}//end MATH
} //end SMF
#endif // __SMF_TRANF_OPER__