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
#include "math/SMF_Transforms.h"
#include "math/all.h"
#include "math/SMF_Vector.h"


namespace SMF{
namespace MATH{

CTranslateOp::CTranslateOp(float tx, float ty, float tz)
:x(tx), y(ty), z(tz)
{
}

CTranslateOp::CTranslateOp(const CVec3D &trans)
:x(trans.x), y(trans.y), z(trans.z)
{
}

CTranslateOp::operator CMatJoint3x4() const
{
	return toMat3x4();
}

CTranslateOp::operator CMat4D() const
{
	return toMat4D();
}


CMatJoint3x4 CTranslateOp::toMat3x4() const
{
	CMatJoint3x4 m;
	m.setRow(0, 1, 0, 0, x);
	m.setRow(1, 0, 1, 0, y);
	m.setRow(2, 0, 0, 1, z);
	return m;
}

CMat4D CTranslateOp::toMat4D() const
{
	CMat4D m;
	m.setRow(0, 1, 0, 0, x);
	m.setRow(1, 0, 1, 0, y);
	m.setRow(2, 0, 0, 1, z);
	m.setRow(3, 0, 0, 0, 1.f);
	return m;
}

CVec3D CTranslateOp::Offset() const
{
	return CVec3D(x, y, z);
}

CMatJoint3x4 operator *(const CTranslateOp &lhs, const CMatJoint3x4 &rhs)
{
	CMatJoint3x4 r = rhs;
	r.setTranslatePart(r.translatePart() + lhs.Offset());

	// Our optimized form of multiplication must be the same as this.
	SMF_ASSERT(r.compare((CMatJoint3x4)lhs * rhs));
	return r;
}

CMatJoint3x4 operator *(const CMatJoint3x4 &lhs, const CTranslateOp &rhs)
{
	CMatJoint3x4 r = lhs;
	r.setTranslatePart(lhs.transformPos(rhs.Offset()));

	// Our optimized form of multiplication must be the same as this.
	SMF_ASSERT(r.compare(lhs * (CMatJoint3x4)rhs));
	return r;
}

CMat4D operator *(const CTranslateOp &lhs, const CMat4D &rhs)
{
	CMat4D r = rhs;
	r.setTranslatePart(r.translatePart() + lhs.Offset());

	// Our optimized form of multiplication must be the same as this.
	SMF_ASSERT(r.compare(lhs.toMat4D() * rhs));
	return r;
}

CMat4D operator *(const CMat4D &lhs, const CTranslateOp &rhs)
{
	CMat4D r = lhs;
	r.setTranslatePart(lhs.transformPos(rhs.Offset()));

	// Our optimized form of multiplication must be the same as this.
	SMF_ASSERT(r.compare(lhs * rhs.toMat4D()));
	return r;
}

CScaleOp::CScaleOp(float sx, float sy, float sz)
:x(sx), y(sy), z(sz)
{
}

CScaleOp::CScaleOp(const CVec3D &scale)
:x(scale.x), y(scale.y), z(scale.z)
{
}

CScaleOp::operator CMat3D() const
{
	return ToFloat3x3();
}

CScaleOp::operator CMatJoint3x4() const
{
	return toMat3x4();
}

CScaleOp::operator CMat4D() const
{
	return toMat4D();
}

CMat3D CScaleOp::ToFloat3x3() const
{
	CMat3D m;
	m.setRow(0, x, 0, 0);
	m.setRow(1, 0, y, 0);
	m.setRow(2, 0, 0, z);
	return m;
}

CMatJoint3x4 CScaleOp::toMat3x4() const
{
	CMatJoint3x4 m;
	m.setRow(0, x, 0, 0, 0);
	m.setRow(1, 0, y, 0, 0);
	m.setRow(2, 0, 0, z, 0);
	return m;
}

CMat4D CScaleOp::toMat4D() const
{
	CMat4D m;
	m.setRow(0, x, 0, 0, 0);
	m.setRow(1, 0, y, 0, 0);
	m.setRow(2, 0, 0, z, 0);
	m.setRow(3, 0, 0, 0, 1.f);
	return m;
}

CMat3D operator *(const CScaleOp &lhs, const CMat3D &rhs)
{
	CMat3D ret = rhs;
	ret.scaleRow(0, lhs.x);
	ret.scaleRow(1, lhs.y);
	ret.scaleRow(2, lhs.z);

	// Our optimized form of multiplication must be the same as this.
	SMF_ASSERT(ret.compare((CMat3D)lhs * rhs));
	return ret;
}

CMat3D operator *(const CMat3D &lhs, const CScaleOp &rhs)
{
	CMat3D ret = lhs;
	ret.scaleCol(0, rhs.x);
	ret.scaleCol(1, rhs.y);
	ret.scaleCol(2, rhs.z);

	// Our optimized form of multiplication must be the same as this.
	SMF_ASSERT(ret.compare(lhs * (CMat3D)rhs));
	return ret;
}

CMatJoint3x4 operator *(const CScaleOp &lhs, const CMatJoint3x4 &rhs)
{
	CMatJoint3x4 ret;
	ret.setRow(0,rhs[0*4+0] * lhs.x,rhs[0*4+1] * lhs.x, rhs[0*4+2] * lhs.x, rhs[0*4+3] * lhs.x);
	ret.setRow(1, rhs[1*4+0] * lhs.y, rhs[1*4+1] * lhs.y, rhs[1*4+2] * lhs.y, rhs[1*4+3] * lhs.y);
	ret.setRow(2, rhs[2*4+0] * lhs.z, rhs[2*4+1] * lhs.z, rhs[2*4+2] * lhs.z, rhs[2*4+3] * lhs.z);

	SMF_ASSERT(ret.compare(lhs.toMat3x4() * rhs));
	return ret;
}

CMatJoint3x4 operator *(const CMatJoint3x4 &lhs, const CScaleOp &rhs)
{
	CMatJoint3x4 ret;
	ret.setRow(0, lhs[0*4+0] * rhs.x, lhs[0*4+1] * rhs.y, lhs[0*4+2] * rhs.z, lhs[0*4+3]);
	ret.setRow(1, lhs[1*4+0] * rhs.x, lhs[1*4+1] * rhs.y, lhs[1*4+2] * rhs.z, lhs[1*4+3]);
	ret.setRow(2, lhs[2*4+0] * rhs.x, lhs[2*4+1] * rhs.y, lhs[2*4+2] * rhs.z, lhs[2*4+3]);

	SMF_ASSERT(ret.compare(lhs * rhs.toMat3x4()));
	return ret;
}

CMat4D operator *(const CScaleOp &lhs, const CMat4D &rhs)
{
	CMat4D ret;
	ret[0][0] = rhs[0][0] * lhs.x; ret[0][1] = rhs[0][1] * lhs.x; ret[0][2] = rhs[0][2] * lhs.x; ret[0][3] = rhs[0][3] * lhs.x;
	ret[1][0] = rhs[1][0] * lhs.y; ret[1][1] = rhs[1][1] * lhs.y; ret[1][2] = rhs[1][2] * lhs.y; ret[1][3] = rhs[1][3] * lhs.y;
	ret[2][0] = rhs[2][0] * lhs.z; ret[2][1] = rhs[2][1] * lhs.z; ret[2][2] = rhs[2][2] * lhs.z; ret[2][3] = rhs[2][3] * lhs.z;
	ret[3][0] = rhs[3][0];		 ret[3][1] = rhs[3][1];		 ret[3][2] = rhs[3][2];		 ret[3][3] = rhs[3][3];

	SMF_ASSERT(ret.compare(lhs.toMat4D() * rhs));
	return ret;
}

CMat4D operator *(const CMat4D &lhs, const CScaleOp &rhs)
{
	CMat4D ret;
	ret[0][0] = lhs[0][0] * rhs.x; ret[0][1] = lhs[0][1] * rhs.y; ret[0][2] = lhs[0][2] * rhs.z; ret[0][3] = lhs[0][3];
	ret[1][0] = lhs[1][0] * rhs.x; ret[1][1] = lhs[1][1] * rhs.y; ret[1][2] = lhs[1][2] * rhs.z; ret[1][3] = lhs[1][3];
	ret[2][0] = lhs[2][0] * rhs.x; ret[2][1] = lhs[2][1] * rhs.y; ret[2][2] = lhs[2][2] * rhs.z; ret[2][3] = lhs[2][3];
	ret[3][0] = lhs[3][0] * rhs.x; ret[3][1] = lhs[3][1] * rhs.y; ret[3][2] = lhs[3][2] * rhs.z; ret[3][3] = lhs[3][3];

	SMF_ASSERT(ret.compare(lhs * rhs.toMat4D()));
	return ret;
}

CMatJoint3x4 operator *(const CScaleOp &lhs, const CTranslateOp &rhs)
{
	CMatJoint3x4 ret;
	ret.setRow(0, lhs.x,	 0,	 0, lhs.x * rhs.x);
	ret.setRow(1,	 0, lhs.y,	 0, lhs.y * rhs.y);
	ret.setRow(2,	 0,	 0, lhs.z, lhs.z * rhs.z);

	SMF_ASSERT(ret.compare(lhs.toMat3x4() * rhs));
	return ret;
}

CMatJoint3x4 operator *(const CTranslateOp &lhs, const CScaleOp &rhs)
{
	CMatJoint3x4 ret;
	ret.setRow(0, rhs.x,	 0,	 0, lhs.x);
	ret.setRow(1,	 0, rhs.y,	 0, lhs.y);
	ret.setRow(2,	 0,	 0, rhs.z, lhs.z);

	SMF_ASSERT(ret.compare(lhs.toMat3x4() * rhs));
	return ret;
}

CVec3D CScaleOp::Offset() const
{
	return CVec3D(x, y, z);
}

} //end MATH
} //end SMF